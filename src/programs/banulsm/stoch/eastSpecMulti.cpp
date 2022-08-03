
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#include "StochasticHamiltonian.h"
#include "Properties.h"
#include "misc.h"

#include "quicksort.h"

#define FREQSAVE 18000 // frequency of tmp saving: default 5 hrs 


/** 
Generates a file name for the MPS corresponding to a certain excitation.
*/
const string mpsfilename(int L,int D,double s,double c,int level,const string mpsdir);
const string jobfilename(int L,int D,double s,double c,const string jobsdir);

/** Compute individual occupations of sites */
void computeOccupations(const MPS& gs,int L,int d,vector<double>& occ);

/** Same for polarizations */
void computePolarizations(const MPS& gs,int L,int d,vector<double>& occ);


/** Check whether there is a file with smaller bond dimension that can be used as initial guess. */
int findInitFile(string& str,int L,int D,double s,double c,int level,const string initmpsdir);

/** Auxiliary object, level-energy, to sort */
class Level{
public:
  int lev;
  double energy;
  Level(int l, double E):lev(l),energy(E){};
  Level(const Level& l2):lev(l2.lev),energy(l2.energy){};
  bool operator<=(const Level& l2){return energy<=l2.energy;}
  bool operator<(const Level& l2){return energy<l2.energy;}
  bool operator==(const Level& l2){return energy==l2.energy;}
  friend bool operator<=(const Level& l1,const Level& l2){return l1.energy<=l2.energy;}
  friend bool operator<(const Level& l1,const Level& l2){return l1.energy<l2.energy;}
  friend bool operator==(const Level& l1,const Level& l2){return l1.energy==l2.energy;}
};



/** eastSpecMulti tries to find the ground state and excitations of
    the stochastic Hamiltonian \ref <StochasticHamiltonian> but,
    instead of a single state, as in basic, to find the GS it tries to
    enforce orthogonality to other states, too. Namely, we search at
    the same time (iterating) for 10 mutually orthogonal states.

    Receives as argument a Properties file, such as
    config/stochspec.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.

    Receives, among others, arguments:
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
  int maxD=props.getIntProperty("maxD");
  int incrD=props.getIntProperty("incrD");
  if(incrD<0) incrD=20;
  double c=props.getDoubleProperty("c");
  bool limitS=0;
  if(c<0||c>1){
    cout<<"Error: value of c not allowed (c in[0,1])"<<endl;
    exit(1);
  }
  int nlev=props.getIntProperty("nlevel");
  int knr=props.getIntProperty("knr");
  if(knr<=0) knr=1;
  string outfname=props.getProperty("outputfile");
  string polfname=props.getProperty("polfile");
  bool computePols=(polfname.length()>0);
  string mpsdir=props.getProperty("mpsdir");
  string initmpsdir=props.getProperty("initmpsdir");
  if(initmpsdir.size()<=0){
    cout<<"ERROR: Need the initmpsdir argument"<<endl;
    exit(1);
  }
  string jobsdir=props.getProperty("jobsdir");
  bool sectorOne=(props.getIntProperty("firstN")!=0);
  cout<<"Occupation of fist site fixed to "<<(sectorOne?1:0)<<endl;
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  double eigtol=props.getDoubleProperty("eigtol");
  if(eigtol<0) eigtol=0.;
  int freqSave=props.getIntProperty("tmpSavingFreq");
  if(freqSave<0) freqSave=FREQSAVE;
  double noise=props.getDoubleProperty("noise");
  if(noise<0) noise=0.01;  
  bool refine=(props.getIntProperty("refine")==1); // even if the level exists, run again
  
  ofstream* out;
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
  if(computePols){
    out=new ofstream(polfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<polfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% Polarization of stochastic model L="<<L<<", D="<<D<<endl;
    *out<<"% s="<<s<<endl;
    *out<<"% nlev\t <H>\t<H^2>\t <sigx(0)>\t <sigx(1)>\t ..."<<endl;
    out->close();delete out;
  }

  //out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
  //out->close();delete out;
  cout<<setprecision(15);
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;  
  contractor.setConvTol(tol);
  // if(D>20){
  //     contractor.setDebugLevel(2);
  //   contractor.setEigenSolver(dmrg);
  // }
  // else{
  if(eigtol>0)contractor.setEigTol(eigtol);
  contractor.setEigenSolver(primme);
    //  }
  int d=2;

  double offset=-L/2;
  StochasticHamiltonian hSto(L-1,s,offset,sectorOne,c);
  const MPO& hamil=hSto.getHMPO();
  //hamil.exportForMatlab("testSto.m");//exit(1);
  
  vector<MPS*> levels;
  vector<Level> unsrtLevs;
  vector<int> Dinit(nlev,D);
  
  // First, try to recover initial guesses for each of the desired
  // states and decide where to start the refining iteration. If there
  // is no file for a level, I stop, so that before running this, I
  // need to run the basic one.
  MPS* nextLev;
  for(int lval=0;lval<nlev;lval++){
    nextLev=new MPS(L-1,D,d);
    const string filename=mpsfilename(L,D,s,c,lval,mpsdir);
    // Try to open file for MPS level lval
    if(file_exists(filename)){
      nextLev->importMPS(filename.data());
    }
    else{ // try tmp file for the level
      const string tmpfilename=filename+".tmp";
      if(file_exists(tmpfilename)){
	nextLev->importMPS(tmpfilename.data());
	cout<<"Imported tmp file for level "<<lval<<" and D="<<D<<endl;
      }
      else{  // no tmp file
	// First check if a valid initial guess exists in the initmpsdir directory
	string filenameInit;
	int Dinit_=findInitFile(filenameInit,L,D,s,c,lval,initmpsdir);
	if(Dinit_>0){
	  cout<<"Using initial guess for (lval,D)=("<<lval<<","<<D<<") from result for "<<Dinit_<<endl;
	  nextLev->importMPS(filenameInit.data());
	  Dinit[lval]=Dinit_;
	}
	else{
	  cout<<"ERROR: Cannot find initial state for level "<<lval<<endl;
	  exit(1);
	}
      }
    }
    nextLev->gaugeCond('R',1);
    nextLev->gaugeCond('L',1);
    levels.push_back(nextLev);
  } // loop for lval
    
  // I will now iterate over the nlev levels, each time making the state orth to all
  if(refine){
    bool done=false;
    while(!done){
      // pass once over each level
      for(int lval=0;lval<nlev;lval++){
	MPS* exc=levels[lval];
	const string filename=mpsfilename(L,D,s,c,lval,mpsdir);
	// For orthogonality, take all the other states
	vector<MPS*> others;
	for(int k=0;k<nlev;k++)
	  if(k!=lval) others.push_back(levels[k]);
	if(Dinit[lval]<D)
	  exc->increaseBondDimensionWithNoise(D,noise);
	// search
	const string tmpfilename=filename+".tmp";
	double lambda=0.;
	contractor.findNextExcitedStateLM(hamil,D,others,&lambda,*exc,0.,knr,tmpfilename,freqSave);
	Dinit[lval]=D;
	exc->gaugeCond('L',1);
	// Save to file
	exc->exportMPS(filename.data());
	// And delete tmp file, if existing
	if(file_exists(tmpfilename)){
	  remove(tmpfilename.data());
	}
	complex_t Ek=contractor.contract(*exc,hamil,*exc)-offset*ONE_c;	
	unsrtLevs.push_back(Level(lval,real(Ek)));
      }   
      // decide if we are done?
      done=true;
    }
  }
  // Now reorder, and write observables  

  // I reorder the states found, rename the mpsfiles and rewrite the expectation values!
  quicksort(unsrtLevs,0,unsrtLevs.size()-1);
  bool chgOrder=0;
  for(int k=0;k<unsrtLevs.size();k++){
    Level& aux=unsrtLevs[k];
    int oldL=aux.lev;
    if(k!=oldL) chgOrder=true;
    cout<<"Level nr "<<k<<" was labeled "<<oldL<<", with energy E="<<aux.energy<<endl;
    const string oldfilename=mpsfilename(L,D,s,c,oldL,mpsdir);
    const string newfilename=mpsfilename(L,D,s,c,k,mpsdir)+".sorting";    
    cout<<"Renaming "<<oldfilename<<"->"<<newfilename<<endl;
    rename(oldfilename.data(),newfilename.data());
  }
  for(int k=0;k<unsrtLevs.size();k++){ // rename again to remove the extension
    const string oldfilename=mpsfilename(L,D,s,c,k,mpsdir)+".sorting";    
    const string newfilename=mpsfilename(L,D,s,c,k,mpsdir);    
    cout<<"Renaming "<<oldfilename<<"->"<<newfilename<<endl;
    rename(oldfilename.data(),newfilename.data());      
    // and write the data to file
    Level& aux=unsrtLevs[k];
    const MPS* exc=levels[aux.lev];
    complex_t Ek=contractor.contract(*exc,hamil,*exc)-offset*ONE_c;	
    complex_t Ek2=contractor.contract2(hamil,*exc)-2*offset*Ek-offset*offset*ONE_c;
    vector<double> occ;
    computeOccupations(*exc,L-1,d,occ);
    out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
    *out<<k<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t";
    *out<<(sectorOne?1:0)<<"\t";
    for(int pos=0;pos<L-1;pos++)
      *out<<occ[pos]<<"\t";
    *out<<endl;
    out->close();delete out;
    if(computePols){
      computePolarizations(*exc,L-1,d,occ);
      out=new ofstream(polfname.data(),ios::app);*out<<setprecision(15);
      *out<<k<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t";
      *out<<0<<"\t";
      for(int pos=0;pos<L-1;pos++)
	*out<<occ[pos]<<"\t";
      *out<<endl;
      out->close();delete out;
    }
  }
}



int findInitFile(string& str,int L,int D,double s,double c,int level,const string initmpsdir){
  cout<<"Searching for init file for level "<<level<<" in dir "<<initmpsdir<<endl;
  for(int myD=D-1;myD>0;myD--){
    // TRy with one smaller D
    str=mpsfilename(L,myD,s,c,level,initmpsdir);
    if(file_exists(str)){
      cout<<"Found file "<<str<<endl;
      return myD;
    }
  }
  // if nothing found
  str="";
  return 0; 
}

const string mpsfilename(int L,int D,double s,double c,int level,const string mpsdir){
  stringstream str;
  str<<mpsdir<<"/MPS_L"<<L<<"_s"<<s<<"_c"<<c<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return str.str();
}

const string jobfilename(int L,int D,double s,double c,const string jobsdir){
  stringstream str;
  str<<jobsdir<<"/job_L"<<L<<"_s"<<s<<"_c"<<c<<"_D"<<D;
  str<<".dat";
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

void computePolarizations(const MPS& gs,int L,int d,vector<double>& polX){
  // Construct a MPO for the sigma_x on a given site
  static MPO mpoId(L); // Basic MPO: all Identity
  static Operator idOp(reshape(identityMatrix(d),Indices(d,1,d,1)));
  static Operator* xOp; // individual operator sigmax
  static bool init(false);
  if(!init){
    complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    mwArray sigX(Indices(d,1,d,1),dataX);
    xOp=new Operator(sigX);
    for(int k=0;k<L;k++)
      mpoId.setOp(k,&idOp,0);
    init=true;
  }
  Contractor& contractor=Contractor::theContractor();
  polX.clear();
  for(int pos=0;pos<L;pos++){
    mpoId.setOp(pos,xOp,false);
    complex_t nk=contractor.contract(gs,mpoId,gs);
    mpoId.setOp(pos,&idOp,false);
    polX.push_back(real(nk));
  }
}


