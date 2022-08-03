
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "mwArray.h"
#include "Properties.h"
#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

/** Processes the results from variationalMinVarThreshold to compute
    rdm of subchains and compare to the infinite temperature ones.
*/

int findTmpFile(int D1,int Dmax,int incrD,const string& mpsdir,const string& baseName,string& mpsfile);
const string mpsFileName(const string& mpsdir,const string& baseName,int D);

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  //  srandom(time(NULL));
  const char* infile=argv[++cntr];
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  int L = props.getIntProperty("L");  
  double J = props.getDoubleProperty("J");  
  double g = props.getDoubleProperty("g");  
  double h = props.getDoubleProperty("h");  
  int D1 = props.getIntProperty("D1");  
  int D2 = props.getIntProperty("D2");
  int incrD = props.getIntProperty("incrD");  
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir"); // where are the MPS
  string mpsfile=props.getProperty("mpsfile"); // basename of MPS files
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  
  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% N\t J\t g\t h\t D\t <H2>\t <H>\t L\t dist"<<endl;
    out->close();
    delete out;
  }
  cout<<setprecision(10);

  
  int d=2;
  IsingHamiltonian hamI(L,d,J,g,h); 
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamI.getHMPO(0.);

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int D=D1;
  
  MPS gs(L,D,d);
  complex_t valH,valH2;
  bool done=0;int failed=0;
  while(!done){
    string mpsDfile=mpsFileName(mpsdir,mpsfile,D);
    if(!file_exists(mpsDfile)){    
      cout<<"File with MPS for D="<<D<<" not found!: "<<mpsDfile<<endl;
      failed++;
    }
    else{
      gs.importMPS(mpsDfile.data());
      // Read the state
      // probably not needed
      gs.gaugeCond('R',1);
      //gs.gaugeCond('L',1);

      cout<<"Imported file for D="<<D<<endl;
      // Just for cross checkign, compute again <H> and <H2>
      valH=contractor.contract(gs,hamil,gs);
      valH2=contractor.contract2(hamil,gs);

      // Instead of doing this, I could use the strategy of
      // Contractor::getEntropy(mps,pos1,pos2), where the only thing we
      // need is the (DxD)^2 dimensional non-orthogonal matrix of
      // overlaps to find the eigenvalues of rho_Lc. Then substract
      // 1/2^Lc to get the trace distance needed. Nevertheless, if D is
      // large, it might be worse. Actually, it will only be
      // advantageous if D^2<d^Lc

    
      MPS auxRho(L,1,d*d);
      // The density operator, to trace out sites
      DMfromMPS(gs,auxRho);
      gs.clear();
      auxRho.gaugeCond('R',false);
      cout<<"Ceated the MPs for rho"<<endl;
      // Now look at the middle L sites, for L=2,4,6
      int Lc=8;int curL=L;
      int posL=(L-Lc)/2;int posR=posL+Lc-1; // outside tracing them
      while(Lc>=2){
	vector<int> toTrace;
	for(int k=0;k<curL;k++){
	  if(k<posL||k>posR)
	    toTrace.push_back(k);
	}
	auxRho.traceOutSites(toTrace);
	//int exactD=pow(d,Lc); // exact bond dimension for this=> can truncate
	auxRho.gaugeCond('R',0); // not normalizing
	cout<<"Traced out sites, now Lc="<<Lc<<endl;
	// convert to operator, expand substract the identity and compute SVD
	// Now, more economic: diag instead of svd
	MPO rdm(Lc);
	MPOfromMPS(auxRho,rdm);
	double singvals=0.;
	double idFc=1./pow(2,Lc);
	{
	  mwArray oper;
	  expandOper(rdm,oper);
	  rdm.clear();
	  cout<<"Built the full operator "<<endl;
	  mwArray U;vector<complex_t> Dv;
	  wrapper::eig(oper,Dv,U); // only eigenvalues
	  // now sum the absolute values
	  for(int k=0;k<Dv.size();k++){
	    double aux=real(Dv[k]);
	    //if(aux>=minEV)
	    singvals+=abs(aux-idFc);
	    //else
	    //singvals+=idFc;
	  }
	  if(Dv.size()<pow(2,Lc)){
	    cout<<"WARNING! Seems the dimension of the operator for D="<<D
		<<", Lc="<<Lc<<"is not full: dim(oper)="<<Dv.size()
		<<", 2^Lc="<<pow(2,Lc)<<endl;
	    singvals+=idFc*(pow(2,Lc)-Dv.size());
	  }
	}
	// mwArray oper;
	// expandOper(rdm,oper);
	// oper=oper-(1./pow(d,Lc))*identityMatrix(pow(d,Lc));
	// mwArray W,S,Vdagger;
	// wrapper::svd(oper,W,S,Vdagger);
	// oper.clear();
	// double singvals=real(S.trace());
	// open the file and write to it
	out=new ofstream(outfname.data(),ios::app);
	*out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"<<real(valH2)<<"\t"
	    <<real(valH)<<"\t"<<Lc<<"\t"<<singvals<<endl;
	out->close();delete out;    
	curL=Lc;Lc-=2;posL=1;posR=curL-2;
      }
    }
    D+=incrD;
    done=(D>D2)||(failed>10);
  }
    
}

const string mpsFileName(const string& mpsdir,const string& baseName,int cnt){
  stringstream s;
  s<<mpsdir<<"/"<<baseName;
  if(cnt>0) s<<"_"<<cnt;
  return s.str();
}

int findTmpFile(int D1,int Dmax,int incrD,const string& mpsdir,const string& baseName,string& mpsfile){
  int D=Dmax;
  bool found=0;
  while(D>=D1&&!found){
    mpsfile=mpsFileName(mpsdir,baseName,D);
    if(file_exists(mpsfile)){
      cout<<"findTmpFile found mpsfile="<<mpsfile<<endl;
      found=true;
      return D;
    }
    else{
      mpsfile="";
      D-=incrD;
    }
  }
  if(!found) return 0;
  cout<<"ERROR: Should not be here"<<endl;
  exit(1);
}
