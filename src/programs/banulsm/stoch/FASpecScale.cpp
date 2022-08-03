
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#include "FAHamiltonian.h"
#include "Properties.h"
#include "misc.h"

#include "quicksort.h"

#define FREQSAVE 18000 // frequency of tmp saving: default 5 hrs 


/** 
Generates a file name for the MPS corresponding to a certain excitation.
*/
const string mpsfilename(int L,int D,double sL,double c,int level,const string mpsdir);
const string jobfilename(int L,int D,double sL,double c,const string jobsdir);

/** Compute individual occupations of sites */
void computeOccupations(const MPS& gs,int L,int d,vector<double>& occ);

/** Same for polarizations */
void computePolarizations(const MPS& gs,int L,int d,vector<double>& occ);

void mpoNN(MPO& mpo,int L);

/** Check whether there is a file with smaller bond dimension that can be used as initial guess. */
int findInitFile(string& str,int L,int D,double sL,double c,int level,const string initmpsdir);

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



/** basicSpec tries to find the ground state and excitations of the stochastic Hamiltonian
    \ref <StochasticHamiltonian>

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
  double sL=props.getDoubleProperty("sL");
  double s=sL/(1.*L);
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
  string jobsdir=props.getProperty("jobsdir");
  bool app=(props.getIntProperty("append")!=0);
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  double eigtol=props.getDoubleProperty("eigtol");
  if(eigtol<0) eigtol=0.;
  int freqSave=props.getIntProperty("tmpSavingFreq");
  if(freqSave<0) freqSave=FREQSAVE;
  double noise=props.getDoubleProperty("noise");
  if(noise<0) noise=0.01;  
  bool reorder=(props.getIntProperty("reorder")==1); // if set, levels for this D are read and ordered in increasing E
  bool origApp=app;
  bool refine=(props.getIntProperty("refine")==1); // even if the level exists, run again
  double penZero=props.getDoubleProperty("penaltyZero");
  if(penZero<0) penZero=0.;
  
  bool rewrite=false;
  ofstream* out;
  if(!app||!file_exists(outfname)){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% Spectrum of stochastic FA model L="<<L<<", D="<<D<<endl;
    *out<<"% sL="<<sL<<endl;
    *out<<"% nlev\t <H>\t<H^2>\t <n(0)>\t <n(1)>\t ..."<<endl;
    out->close();delete out;
    rewrite=true;
    if(computePols){
      out=new ofstream(polfname.data());
      if(!out->is_open()){
	cout<<"Error: impossible to open file "<<polfname<<
	  " for output"<<endl;
	exit(1);
      }
      *out<<"% Polarization of FA stochastic model L="<<L<<", D="<<D<<endl;
      *out<<"% sL="<<sL<<endl;
      *out<<"% nlev\t <H>\t<H^2>\t <sigx(0)>\t <sigx(1)>\t ..."<<endl;
      out->close();delete out;
    }

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

  double offset=0.;
  FAHamiltonian hSto(L,s,offset,c,penZero);
  const MPO& hamil=hSto.getHMPO();
  //hamil.exportForMatlab("testFA.m");//exit(1);
  cout<<"Init Hamil"<<endl;
  
  // Want to compute nn
  MPO opNN(L);
  mpoNN(opNN,L);
  //  opNN.exportForMatlab("testNNop.m");exit(1);

  
  vector<MPS*> levels;
  vector<Level> unsrtLevs;
  
  // First, check how many MPS I already have in mpsdir
  bool oldMPS=true;
  int lval=0;
  MPS* nextLev;

  //  double lastE=0.;
  //double offset=0.;
  if(!refine){ // otherwise, rerun)
    while(oldMPS&&lval<=nlev-1){
      const string filename=mpsfilename(L,D,sL,c,lval,mpsdir);
      // Try to open file for MPS level lval
      if(file_exists(filename)){
	nextLev=new MPS(L,D,d);
	nextLev->importMPS(filename.data());
	nextLev->gaugeCond('R',1);
	nextLev->gaugeCond('L',1);
	levels.push_back(nextLev);
	complex_t Ek=contractor.contract(*nextLev,hamil,*nextLev)-offset*ONE_c;	
	complex_t nnk=contractor.contract(*nextLev,opNN,*nextLev);	
	complex_t Ek2=contractor.contract2(hamil,*nextLev)-2*offset*Ek-offset*offset*ONE_c;
	unsrtLevs.push_back(Level(lval,real(Ek)));
	cout<<"Recovered MPS file for level "<<lval<<" from file "<<filename<<endl;
	if(rewrite){
	  vector<double> occ;
	  computeOccupations(*nextLev,L,d,occ);
	  out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
	  *out<<lval<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t"<<real(nnk)<<"\t";
	  for(int pos=0;pos<L;pos++)
	    *out<<occ[pos]<<"\t";
	  *out<<endl;
	  out->close();delete out;
	  if(computePols){
	    computePolarizations(*nextLev,L,d,occ);
	    out=new ofstream(polfname.data(),ios::app);*out<<setprecision(15);
	    *out<<lval<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t";
	    *out<<0<<"\t";
	    for(int pos=0;pos<L;pos++)
	      *out<<occ[pos]<<"\t";
	    *out<<endl;
	    out->close();delete out;
	  }
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
  }
  
  //  offset=-lastE+2/L;
  cout<<"Recovered "<<lval<<" already computed MPS levels from disk. "
      <<"Computing next level now."<<endl;
  //cout<<"Last energy was "<<lastE<<" so now I will apply an offset "<<offset<<endl;

  while(lval<=nlev-1){
    MPS* exc=new MPS(L,D,d);
    double absEk=0.;
    const string filename=mpsfilename(L,D,sL,c,lval,mpsdir); // where I will save it
    // First check if a temporary file for this state exists in the initmpsdir directory
    const string tmpfilename=filename+".tmp";
    bool found=0;
    if(refine){ // first check for the same name
      if(file_exists(filename)){
	exc->importMPS(filename.data());
	cout<<"Imported initial file for same D!"<<endl;
	found=true;
	absEk=real(contractor.contract(*exc,hamil,*exc))-offset;
      }
    }
    if(!found){
      if(file_exists(tmpfilename)){
	cout<<"There is a tmp file for same D: "<<tmpfilename<<endl;
	exc->importMPS(tmpfilename.data());
	cout<<"Imported tmp file for same D!"<<endl;
	found =true;
	absEk=real(contractor.contract(*exc,hamil,*exc))-offset;
      }
      else{  // no tmp file
	cout<<"There is no tmp file for same D: "<<tmpfilename<<endl;    
	// First check if a valid initial guess exists in the initmpsdir directory
	if(initmpsdir.size()>0){
	  string filenameInit;
	  int Dinit=findInitFile(filenameInit,L,D,sL,c,lval,initmpsdir);
	  if(Dinit>0){
	    cout<<"Using initial guess for D="<<D<<" from result for "<<Dinit<<" in file "
		<<filenameInit<<endl;
	    exc->importMPS(filenameInit.data());
	    absEk=real(contractor.contract(*exc,hamil,*exc))-offset;
	    if(Dinit<D)
	      exc->increaseBondDimensionWithNoise(D,noise);
	    found=true;
	  }
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

    // Maybe try the rescaled H
    MPO hamilSc(L);
    double scale=1.;
    const MPO* hamilPtr=&hamil;
    if(found&&absEk!=0){
      if(lval>0||(absEk>1E-8)){
	if(absEk*tol<1E-13*L){
	  scale=1./max(1E-8,absEk);
	  double factorSc=pow(scale,1./L);
	  for(int k=0;k<L;k++){
	    hamilSc.setOp(k,new Operator(factorSc*hamil.getOp(k).getFullData()),true);
	  }
	  hamilPtr=&hamilSc;
	}
      }
    }
    if(lval==0){
      if(s<1){
	contractor.findGroundStateLM(*hamilPtr,D,&lambda,*exc,-L*scale,knr,tmpfilename,freqSave);
	//*	//	contractor.findGroundStateTwoSite(hamil,D,&lambda,*exc,-L,knr,tmpfilename,freqSave,true);
	//*	//	contractor.findGroundStateMultiSite(2,hamil,D,&lambda,*exc,-L,knr,tmpfilename,freqSave,true);
      }
      else
	contractor.findGroundState(*hamilPtr,D,&lambda,*exc,0.,knr,tmpfilename,freqSave);
      cout<<"Ground state found with eigenvalue "<<lambda<<endl;
    }
    else{
      contractor.findNextExcitedState(*hamilPtr,D,levels,&lambda,*exc,0.,knr,tmpfilename,freqSave);
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
    unsrtLevs.push_back(Level(lval,real(Ek)));
    complex_t nnk=contractor.contract(*exc,opNN,*exc);	
    //    lastE=real(Ek);offset=lastE+2.;
    complex_t Ek2=contractor.contract2(hamil,*exc)-2*offset*Ek-offset*offset*ONE_c;
    vector<double> occ;
    computeOccupations(*exc,L,d,occ);
    out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
    *out<<lval<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t"<<real(nnk)<<"\t";
    for(int pos=0;pos<L;pos++)
      *out<<occ[pos]<<"\t";
    *out<<endl;
    out->close();delete out;
    cout<<"For level "<<lval<<" E="<<Ek<<", DeltaH="<<real(Ek2)/real(Ek*Ek)-1.<<endl;
    if(computePols){
      computePolarizations(*exc,L,d,occ);
      out=new ofstream(polfname.data(),ios::app);*out<<setprecision(15);
      *out<<lval<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t";
      *out<<0<<"\t";
      for(int pos=0;pos<L;pos++)
	*out<<occ[pos]<<"\t";
      *out<<endl;
      out->close();delete out;
    }
    lval++;
  }

  if(reorder){
    // I reorder the states found, rename the mpsfiles and rewrite the expectation values!
    quicksort(unsrtLevs,0,unsrtLevs.size()-1);
    bool chgOrder=1;  // Now force always rewriting 0;
    for(int k=0;k<unsrtLevs.size();k++){
      Level& aux=unsrtLevs[k];
      int oldL=aux.lev;
      if(k!=oldL) chgOrder=true;
      cout<<"Level nr "<<k<<" was labeled "<<oldL<<", with energy E="<<aux.energy<<endl;
      const string oldfilename=mpsfilename(L,D,sL,c,oldL,mpsdir);
      const string newfilename=mpsfilename(L,D,sL,c,k,mpsdir)+".sorting";    
      cout<<"Renaming "<<oldfilename<<"->"<<newfilename<<endl;
      rename(oldfilename.data(),newfilename.data());
    }
    if(chgOrder){ // recreate the file for expectation values
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
	*out<<"% Polarization of FA stochastic model L="<<L<<", D="<<D<<endl;
	*out<<"% sL="<<sL<<endl;
	*out<<"% nlev\t <H>\t<H^2>\t <sigx(0)>\t <sigx(1)>\t ..."<<endl;
	out->close();delete out;
      }
    }
    for(int k=0;k<unsrtLevs.size();k++){ // rename again to remove the extension
      const string oldfilename=mpsfilename(L,D,sL,c,k,mpsdir)+".sorting";    
      const string newfilename=mpsfilename(L,D,sL,c,k,mpsdir);    
      cout<<"Renaming "<<oldfilename<<"->"<<newfilename<<endl;
      rename(oldfilename.data(),newfilename.data());      
      if(chgOrder){ // and rewrite the data if needed
	Level& aux=unsrtLevs[k];
	const MPS* exc=levels[aux.lev];
	complex_t Ek=contractor.contract(*exc,hamil,*exc)-offset*ONE_c;	
	complex_t nnk=contractor.contract(*exc,opNN,*exc);	
	complex_t Ek2=contractor.contract2(hamil,*exc)-2*offset*Ek-offset*offset*ONE_c;
	vector<double> occ;
	computeOccupations(*exc,L,d,occ);
	out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
	*out<<k<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t"<<real(nnk)<<"\t";
	for(int pos=0;pos<L;pos++)
	  *out<<occ[pos]<<"\t";
	*out<<endl;
	out->close();delete out;
	if(computePols){
	  computePolarizations(*exc,L,d,occ);
	  out=new ofstream(polfname.data(),ios::app);*out<<setprecision(15);
	  *out<<k<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t";
	  *out<<0<<"\t";
	  for(int pos=0;pos<L;pos++)
	    *out<<occ[pos]<<"\t";
	  *out<<endl;
	  out->close();delete out;
	}
      }
    }
  }



  // Prepare the script with the new job
  if(jobsdir.length()>0&&D<maxD){
    int newD=D+incrD;
    string newOutput(outfname);
    stringstream newSuff;newSuff<<"_D"<<newD; 
    newOutput.replace(newOutput.find("_D"),string::npos,newSuff.str());
    stringstream commandstr;
    commandstr<<argv[0]; // exec
    commandstr<<" "<<argv[1]; // config file
    commandstr<<" -L="<<L<<" -sL="<<sL;
    commandstr<<" -D="<<newD;
    commandstr<<" -maxD="<<maxD;
    commandstr<<" -incrD="<<incrD;
    commandstr<<" -c="<<c<<" -nlevel="<<nlev;
    commandstr<<" -knr="<<knr;
    commandstr<<" -outputfile="<<newOutput;
    if(computePols){
      string newPolOutput(polfname);
      newPolOutput.replace(newPolOutput.find("_D"),string::npos,newSuff.str());
      commandstr<<" -polfile="<<newPolOutput;
    }
    commandstr<<" -mpsdir="<<mpsdir<<" -initmpsdir="<<mpsdir;
    commandstr<<" -jobsdir="<<jobsdir;
    commandstr<<" -append="<<origApp;
    commandstr<<" -tol="<<tol;
    if(eigtol>0) commandstr<<" -eigtol="<<eigtol;
    commandstr<<" -tmpSavingFreq="<<freqSave;
    commandstr<<" -noise="<<noise;
    commandstr<<" -reorder="<<reorder;
    commandstr<<" -refine="<<refine;
    commandstr<<" -penaltyZero="<<penZero;
    string jobfile=jobfilename(L,newD,sL,c,jobsdir);

    ofstream* ojob=new ofstream(jobfile.data());
    if(!ojob->is_open()){
      cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
      cout<<commandstr.str()<<endl;
      cout<<endl;
    }
    else{
      int nPar=6;
      if(newD>=80) nPar=10;
      if(newD>=100) nPar=12;
      if(newD>=200) nPar=20;

      *ojob<<"msub_slurm -N fa.c"<<c<<".N"<<L<<".D"<<newD<<".s"<<s;
      *ojob<<" -P "<<nPar;
      *ojob<<" -- ";
      *ojob<<commandstr.str()<<endl;
      ojob->close();
    }
    delete ojob;

  }

  
}


int findInitFile(string& str,int L,int D,double sL,double c,int level,const string initmpsdir){
  cout<<"Searching for init file for level "<<level<<" in dir "<<initmpsdir<<endl;
  for(int myD=D-1;myD>0;myD--){
    // TRy with one smaller D
    str=mpsfilename(L,myD,sL,c,level,initmpsdir);
    if(file_exists(str)){
      cout<<"Found file "<<str<<endl;
      return myD;
    }
  }
  // if nothing found
  str="";
  return 0; 
}

const string mpsfilename(int L,int D,double sL,double c,int level,const string mpsdir){
  stringstream str;
  str<<mpsdir<<"/MPS_L"<<L<<"_sL"<<sL<<"_eps"<<c<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return str.str();
}

const string jobfilename(int L,int D,double sL,double c,const string jobsdir){
  stringstream str;
  str<<jobsdir<<"/job_fa_L"<<L<<"_sL"<<sL<<"_c"<<c<<"_D"<<D;
  str<<".dat";
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return str.str();
}

void mpoNN(MPO& mpo,int L){
  mpo.initLength(L);
  int d=2;
  // basic spin operators appearing
  mwArray sig0=identityMatrix(d);
  complex_t dataN[]={ZERO_c,ZERO_c,ZERO_c,ONE_c};
  mwArray operN(Indices(d,d),dataN);
  mwArray Z(Indices(2,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(operN.getElement(Indices(i1,i2)),Indices(1,i1,i2));
    }
  Z.reshape(Indices(2,d*d));

  int D=3;
  for(int k=0;k<L;k++){
    int Dl=(k==0)?1:D;
    int Dr=(k==L-1)?1:D;
    mwArray C(Indices(Dl,Dr,2));
    if(k<L-1){
      C.setElement(ONE_c,Indices(0,0,0)); // nothing yet
      C.setElement(ONE_c,Indices(0,1,1)); // 1st of pair
    }
    if(k!=0){
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // after everything
      C.setElement(ONE_c,Indices(1,Dr-1,1)); // 2nd of pair
    }
    // Reshape, multiply operators and set in MPO
    C.reshape(Indices(Dl*Dr,2));
    C.multiplyRight(Z);
    C.reshape(Indices(Dl,Dr,d,d));
    C.permute(Indices(3,1,4,2));
    mpo.setOp(k,new Operator(C),true);
  }
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
    //    cout<<"Computed occ ("<<pos<<") "<<nk<<endl;
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

void constructProductThetaPlus(MPS& mps,int L,int d,double c){
  mps=MPS(L,1,d);

  mwArray polC(Indices(d,1));
  polC.setElement(sqrt(1-c)*ONE_c,Indices(0,0));
  polC.setElement(sqrt(c)*ONE_c,Indices(1,0));
  polC.reshape(Indices(2,1,1));
  for(int k=0;k<L;k++)
      mps.replaceSite(k,polC,false);
  //  cout<<"Constructed exact product state for s=0"<<endl;

}

