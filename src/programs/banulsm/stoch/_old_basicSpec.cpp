
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#include "StochasticHamiltonian.h"
#include "Properties.h"
#include "misc.h"

#define FREQSAVE 18000 // frequency of tmp saving: default 5 hrs 


/** 
Generates a file name for the MPS corresponding to a certain excitation.
*/
const string mpsfilename(int L,int D,double s,int level,const string mpsdir);

/** Compute individual occupations of sites */
void computeOccupations(const MPS& gs,int L,int d,vector<double>& occ);


/** Check whether there is a file with smaller bond dimension that can be used as initial guess. */
int findInitFile(string& str,int L,int D,double s,int level,const string initmpsdir);

/** basicSpec tries to find the ground state and excitations of the stochastic Hamiltonian
    \ref <StochasticHamiltonian>

    Receives as argument a Properties file, such as
    config/stochspec.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <s> (double) parameter \f$s\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)

*/

using namespace shrt;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  int L=props.getIntProperty("L");
  double s=props.getDoubleProperty("s");
  int D=props.getIntProperty("D");
  int nlev=props.getIntProperty("nlevel");
  string outfname=props.getProperty("outputfile");
  string mpsdir=props.getProperty("mpsdir");
  string initmpsdir=props.getProperty("initmpsdir");
  bool app=(props.getIntProperty("append")!=0);
  bool sectorOne=(props.getIntProperty("firstN")!=0);
  cout<<"Occupation of fist site fixed to "<<(sectorOne?1:0)<<endl;
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  int freqSave=props.getIntProperty("tmpSavingFreq");
  if(freqSave<0) freqSave=FREQSAVE;
  double noise=props.getDoubleProperty("noise");
  if(noise<0) noise=0.01;  

  bool rewrite=false;
  ofstream* out;
  if(!app||!file_exists(outfname)){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% Spectrum of stochastic model L="<<L<<", D="<<D<<endl;
    *out<<"% s="<<s<<endl;
    *out<<"% nlev\t <H>\t<H^2>\t <n(0)>\t <n(1)>\t ..."<<endl;
    out->close();delete out;
    rewrite=true;
  }
  //out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
  //out->close();delete out;
    
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);
  int d=2;

  double offset=-L/2;
  StochasticHamiltonian hSto(L-1,s,offset,sectorOne);
  const MPO& hamil=hSto.getHMPO();
  //  hamil.exportForMatlab("testSto.m");//exit(1);
  
  vector<MPS*> levels;

  // First, check how many MPS I already have in mpsdir
  bool oldMPS=true;
  int lval=0;
  MPS* nextLev;
  //  double lastE=0.;
  //double offset=0.;
  while(oldMPS&&lval<=nlev-1){
    const string filename=mpsfilename(L,D,s,lval,mpsdir);
    // Try to open file for MPS level lval
    if(file_exists(filename)){
      nextLev=new MPS(L-1,D,d);
      nextLev->importMPS(filename.data());
      nextLev->gaugeCond('R',1);
      nextLev->gaugeCond('L',1);
      levels.push_back(nextLev);
      cout<<"Recovered MPS file for level "<<lval<<" from file "<<filename<<endl;
      if(rewrite){
	complex_t Ek=contractor.contract(*nextLev,hamil,*nextLev)-offset*ONE_c;	
	complex_t Ek2=contractor.contract2(hamil,*nextLev)-2*offset*Ek-offset*offset*ONE_c;
	vector<double> occ;
	computeOccupations(*nextLev,L-1,d,occ);
	out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
	*out<<lval<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t";
	*out<<(sectorOne?1:0)<<"\t";
	for(int pos=0;pos<L-1;pos++)
	  *out<<occ[pos]<<"\t";
	*out<<endl;
	out->close();delete out;
	//	lastE=real(Ek);
      }
      // else{
      // 	lastE=real(contractor.contract(*nextLev,hamil,*nextLev));
      // }
      lval++;
    }
    else{
      oldMPS=false; // break the loop 
    }
  }
  //  offset=-lastE+2/L;
  cout<<"Recovered "<<lval<<" already computed MPS levels from disk. "
      <<"Computing next level now."<<endl;
  //cout<<"Last energy was "<<lastE<<" so now I will apply an offset "<<offset<<endl;

  while(lval<=nlev-1){
    MPS* exc=new MPS(L-1,D,d);
    const string filename=mpsfilename(L,D,s,lval,mpsdir); // where I will save it
    // First check if a temporary file for this state exists in the initmpsdir directory
    const string tmpfilename=filename+".tmp";
    bool found=0;
    if(file_exists(tmpfilename)){
      cout<<"There is a tmp file for same D: "<<tmpfilename<<endl;
      exc->importMPS(tmpfilename.data());
      cout<<"Imported tmp file for same D!"<<endl;
      found =true;
    }
    else{  // no tmp file
      cout<<"There is no tmp file for same D: "<<tmpfilename<<endl;    
      // First check if a valid initial guess exists in the initmpsdir directory
      if(initmpsdir.size()>0){
	string filenameInit;
	int Dinit=findInitFile(filenameInit,L,D,s,lval,initmpsdir);
	if(Dinit>0){
	  cout<<"Using initial guess for D="<<D<<" from result for "<<Dinit<<" in file "
	      <<filenameInit<<endl;
	  exc->importMPS(filenameInit.data());
	  if(Dinit<D)
	    exc->increaseBondDimensionWithNoise(D,noise);
	  found=true;
	}
      }
    }
    if(!found){
      cout<<"Initializing the search for level "<<lval<<" with a random MPS"<<endl;
      exc->setRandomState(); // the initial state, random -> not to the GS!!
    }
    exc->gaugeCond('R',1);
    exc->gaugeCond('L',1);
    double lambda=0.;

    if(lval==0){
      contractor.findGroundState(hamil,D,&lambda,*exc,0.,1,tmpfilename,freqSave);
      cout<<"Ground state found with eigenvalue "<<lambda<<endl;
    }
    else{
      contractor.findNextExcitedState(hamil,D,levels,&lambda,*exc,0.,1,tmpfilename,freqSave);
      cout<<"Next("<<lval<<") excited state found with E="<<lambda<<endl;
    }
    exc->gaugeCond('L',1);
    // Save to file
    exc->exportMPS(filename.data());

    // And delete tmp file, if existing
    if(file_exists(tmpfilename)){
      remove(tmpfilename.data());
    }
    
    levels.push_back(exc);
    complex_t Ek=contractor.contract(*exc,hamil,*exc)-offset*ONE_c;
    //    lastE=real(Ek);offset=lastE+2.;
    complex_t Ek2=contractor.contract2(hamil,*exc)-2*offset*Ek-offset*offset*ONE_c;
    vector<double> occ;
    computeOccupations(*exc,L-1,d,occ);
    out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
    *out<<lval<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t";
    *out<<(sectorOne?1:0)<<"\t";
    for(int pos=0;pos<L-1;pos++)
      *out<<occ[pos]<<"\t";
    *out<<endl;
    out->close();delete out;
    cout<<"Computed and saved occupations"<<endl;
    lval++;
  }

}


int findInitFile(string& str,int L,int D,double s,int level,const string initmpsdir){
  cout<<"Searching for init file for level "<<level<<" in dir "<<initmpsdir<<endl;
  for(int myD=D-1;myD>0;myD--){
    // TRy with one smaller D
    str=mpsfilename(L,myD,s,level,initmpsdir);
    if(file_exists(str)){
      cout<<"Found file "<<str<<endl;
      return myD;
    }
  }
  // if nothing found
  str="";
  return 0; 
}

const string mpsfilename(int L,int D,double s,int level,const string mpsdir){
  stringstream str;
  str<<mpsdir<<"/MPS_L"<<L<<"_s"<<s<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return str.str();
}

void computeOccupations(const MPS& gs,int L,int d,vector<double>& occ){
  // Construct a MPO for the particle number on a given site
  static MPO mpoId(L); // Basic MPO: all Identity
  static Operator idOp(reshape(identityMatrix(d),Indices(d,1,d,1)));
  static Operator* nOp; // individual operator number
  static bool init(false);
  if(!init){
    mwArray nSite(Indices(d,d));
    for(int k=0;k<d;k++){
      nSite.setElement(k*ONE_c,Indices(k,k));
    }
    nSite.reshape(Indices(d,1,d,1));
    nOp=new Operator(nSite);
    for(int k=0;k<L;k++)
      mpoId.setOp(k,&idOp,0);
    init=true;
  }
  Contractor& contractor=Contractor::theContractor();
  occ.clear();
  for(int pos=0;pos<L;pos++){
    mpoId.setOp(pos,nOp,false);
    complex_t nk=contractor.contract(gs,mpoId,gs);
    mpoId.setOp(pos,&idOp,false);
    occ.push_back(real(nk));
  }
}


