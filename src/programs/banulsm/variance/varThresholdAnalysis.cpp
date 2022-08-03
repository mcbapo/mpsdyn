
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "mwArray.h"
#include "Properties.h"
#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"
#include "JoinedOperator.h"

using namespace shrt;

/** Processes the results from variationalMinVarThreshold to compute
    rdm of subchains and compare to the infinite temperature ones.
*/

void buildMixedMPO(const MPS& mps,int posL,int posR,MPO& mpo);
int findTmpFile(int D1,int Dmax,int incrD,const string& mpsdir,const string& baseName,string& mpsfile);
const string mpsFileName(const string& mpsdir,const string& baseName,int D);
void expandRDM(const MPS& mps,int posL,int posR,mwArray& A);
double computeRDMdistanceToId(const mwArray& rdm,int Lc);

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
  bool avrg=1; // whether to compute average distances for each Lc
  avrg=!(props.getIntProperty("average")==0); // replace only if 0 explicitly given (ow -1)
  
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
      cout<<"File with MPS for D="<<D<<" found and read: "<<mpsDfile<<endl;
      // Read the state
      gs.gaugeCond('R',1);
      // Just for cross checking, compute again <H> and <H2>
      valH=contractor.contract(gs,hamil,gs);
      valH2=contractor.contract2(hamil,gs);
      int minPosR=L/2; // assuming L even and smallest Lc is 2
      for(int k=L-1;k>=minPosR+1;k--){
        gs.gaugeCond(k,'L',true);
      }
      //cout<<"Imported file for D="<<D<<" and imposed proper gauge "<<endl;

      // Instead of doing this, I could use the strategy of
      // Contractor::getEntropy(mps,pos1,pos2), where the only thing we
      // need is the (DxD)^2 dimensional non-orthogonal matrix of
      // overlaps to find the eigenvalues of rho_Lc. Then substract
      // 1/2^Lc to get the trace distance needed. Nevertheless, if D is
      // large, it might be worse. Actually, it will only be
      // advantageous if D^2<d^Lc

      int Lcmax=8;
      int Lb=6; // boundary discarded
      if(!avrg){
	int Lc=Lcmax;
	int posL=(L-Lc)/2;int posR=posL+Lc-1; // outside tracing them
	MPO rdm(Lc);
	while(Lc>=2){
	  //mwArray oper;
	  //expandRDM(gs,posL,posR,oper);
	  mwArray oper=contractor.getRDM(gs,posL,posR);
	  double singvals=0.;
	  double idFc=1./pow(2,Lc);
	  {	  
	    mwArray U;vector<complex_t> Dv;
	    wrapper::eig(oper,Dv,U); // only eigenvalues
	    // now sum the absolute values
	    for(int k=0;k<Dv.size();k++){
	      double aux=real(Dv[k]);
	      singvals+=abs(aux-idFc);
	    }
	    if(Dv.size()<pow(2,Lc)){
	      cout<<"WARNING! Seems the dimension of the operator for D="<<D
		  <<", Lc="<<Lc<<"is not full: dim(oper)="<<Dv.size()
		  <<", 2^Lc="<<pow(2,Lc)<<endl;
	      singvals+=idFc*(pow(2,Lc)-Dv.size());
	    }
	  }
	  out=new ofstream(outfname.data(),ios::app);
	  *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"<<real(valH2)<<"\t"
	      <<real(valH)<<"\t"<<Lc<<"\t"<<singvals<<endl;
	  out->close();delete out;    
	  posL++;posR--;Lc=(posR-posL+1);
	}
      }
      else{ // average for various subsystems
	vector<int> Lcs;
	for(int lc=1;lc<=Lcmax;lc++){Lcs.push_back(lc);} // position lc-1
	vector<double> avrgDist(Lcs.size(),0.);
	vector<int> avrgCnt(Lcs.size(),0.);
	int posL=Lb;
	while(posL<L-Lb-1){ // at least there is a 2-site RDM to compute
	  int lc=Lcmax;
	  int posR=posL+lc-1;
	  while(posR>=L){ // cannot compute this one, need to reduce the size
	    lc--;
	    posR=posL+lc-1;
	  }
	  int dimL=pow(d,lc);
	  int dimR=dimL;
	  mwArray oper=contractor.getRDM(gs,posL,posR);
	  while(lc>=1){ // Compute the distance and trace another site
	    if(posR<L-Lb){
	      cout<<"Computing distance for rmd("<<lc<<") sites, pos:"<<posL<<"-"<<posR<<endl;
	      double dist_lc=computeRDMdistanceToId(oper,lc);
	      avrgDist[lc-1]+=dist_lc;
	      avrgCnt[lc-1]++; // to take the average later
	    }
	    if(lc>=1){
	      // trace out two sites
	      oper.reshape(Indices(dimL/d,d,dimR/d,d));
	      oper.permute(Indices(1,3,2,4));
	      dimL=dimL/d;dimR=dimR/d;
	      oper.reshape(Indices(dimL*dimR,d*d));
	      oper.multiplyRight(reshape(identityMatrix(d),Indices(d*d,1)));
	      oper.reshape(Indices(dimL,dimR));
	      posR--;
	    }
	    lc--;
	  }	    
	  posL++;
	}
	// Now write averages to the file
	out=new ofstream(outfname.data(),ios::app);
	for(int lc=1;lc<=Lcmax;lc++){
	  *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"<<real(valH2)<<"\t"
	      <<real(valH)<<"\t"<<lc<<"\t"<<avrgDist[lc-1]/avrgCnt[lc-1]<<endl;
	}
	out->close();delete out;    

      }    
    }
    D+=incrD;
    done=(D>D2)||(failed>10);
  }
    
}

double computeRDMdistanceToId(const mwArray& rdm,int Lc){
  double singvals=0.;
  double idFc=1./pow(2,Lc);
  mwArray U;vector<complex_t> Dv;
  wrapper::eig(rdm,Dv,U); // only eigenvalues
  // now sum the absolute values
  for(int k=0;k<Dv.size();k++){
    double aux=real(Dv[k]);
    singvals+=abs(aux-idFc);
  }
  if(Dv.size()<pow(2,Lc)){
    cout<<"WARNING! Seems the dimension of the operator for Lc="<<Lc
	<<"is not full: dim(oper)="<<Dv.size()
	<<", 2^Lc="<<pow(2,Lc)<<endl;
    singvals+=idFc*(pow(2,Lc)-Dv.size());
  }
  return singvals;
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



void buildMixedMPO(const MPS& mps,int posL,int posR,MPO& mpo){
  // Need to have gauge cond in mps: from left up to posL-1 and to right up to posR+1
  // I assume it is so
  int length=mps.getLength();
  if(posL<0||posR<0||posL>=length||posR>=length||posL>posR){
    cout<<"Error: trying to compute entropy between sites "<<posL
	<<" and "<<posR<<" when length is "<<length<<endl;
    exit(212);
  }
  int Lc=posR-posL+1;
  mpo.initLength(Lc);
  // Edge tensors
  {//left
    mwArray A=mps.getA(posL).getA(); // dxDlxDr
    int d=A.getDimension(0);int Dl=A.getDimension(1);int Dr=A.getDimension(2);
    A.reshape(Indices(d,1,Dl,Dr));
    Operator opKet(A);
    Operator opBra(A,Indices(3,2,1,4),true);
    const Operator* ops[2]={&opKet,&opBra};
    mpo.setOp(0,new JoinedOperator(2,ops),true);
  }
  {//right
    mwArray A=mps.getA(posR).getA(); // dxDlxDr
    int d=A.getDimension(0);int Dl=A.getDimension(1);int Dr=A.getDimension(2);
    A.reshape(Indices(d,Dl,Dr,1));
    Operator opKet(A);
    Operator opBra(A,Indices(3,2,1,4),true);
    const Operator* ops[2]={&opKet,&opBra};
    mpo.setOp(Lc-1,new JoinedOperator(2,ops),true);
  }
  // in between
  for(int k=1;k<Lc-1;k++){
    mwArray A=mps.getA(posL+k).getA(); // dxDlxDr
    int d=A.getDimension(0);int Dl=A.getDimension(1);int Dr=A.getDimension(2);
    A.reshape(Indices(d,Dl,1,Dr));
    Operator opKet(A);
    Operator opBra(A,Indices(3,2,1,4),true);
    const Operator* ops[2]={&opKet,&opBra};
    mpo.setOp(k,new JoinedOperator(2,ops),true);
  }
}


void expandRDM(const MPS& mps,int posL,int posR,mwArray& A){
  A=mps.getA(posL).getA();
  int d=A.getDimension(0);int Dl0=A.getDimension(1);int Dr=A.getDimension(2);
  A.permute(Indices(2,1,3));A.reshape(Indices(Dl0*d,Dr));
  int dphys=d; // temporary dimensions of the vector
  int Dv=Dr;
  for(int k=posL+1;k<=posR;k++){
    mwArray aux=mps.getA(k).getA();
    //cout<<"Multiplying tensor for site "<<k<<", A="<<aux.getDimensions()<<" with tmp result "<<A.getDimensions()<<endl;
    Indices dims=aux.getDimensions();
    int d=dims[0];int Dl=dims[1];int Dr=dims[2];
    if(Dl!=Dv){
      cout<<"Error: Dimensions do not agree in expandVec!"<<endl;
      exit(1);
    }
    aux.permute(Indices(2,1,3));
    aux.reshape(Indices(Dl,d*Dr));
    A.multiplyRight(aux);
    dphys=dphys*d;
    Dv=Dr;
    A.reshape(Indices(Dl0*dphys,Dv));
  }
  A.reshape(Indices(Dl0,dphys,Dv));A.permute(Indices(2,1,3));A.reshape(Indices(dphys,Dl0*Dv));
  // if the mps has a normalization factor, we need to include it!!
  A=mps.getNormFact()*A;
  A.multiplyRight(Hconjugate(A));
}  
