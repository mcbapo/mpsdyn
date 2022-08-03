
#include <math.h>
#include <iomanip>

#include "MPS.h"
//#include "FoldedOperator.h"
#include "DoubleOperator.h"
#include "Contractor.h"
#include "Properties.h"
#include "misc.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

#define THRESHOLD 1E-6 // default threshold for probabilities

//#define DEBUGINFO
#undef DEBUGINFO
#define DEBUGPDF
//#undef DEBUGPDF

/** 
    MC sampling for the canonical expectation values. It is analogous
    to microcanonicalSampling, only for fixed beta instead of energy.
    Different to that program, though, this one requires the lowest
    and highest energies in the spectrum of the Hamiltonian. They
    could be computed as MPS by an independent program, and saved, so
    tht I can now read them. If not available, they will be
    approximated before starting.

*/


/** Auxiliary class to store, update the spin configuration */
class SamplingConfig{
  MPS state;
  int oldSpinId,spinId; // to ease recovering the state index
  int moved; // just to move one spin(to be changed for more general)
  mwArray sigX; // the rotation
private:
  SamplingConfig();
  SamplingConfig(const SamplingConfig& s);
  //void prepareMPS();
public:
  ~SamplingConfig(){};
  static SamplingConfig& config(); // I make it a singleton
  void resetWithMPS(const MPS& mps); // initialize with a certain config
  // propose a random change in configuration
  void proposeMove();
  void confirmLast();
  void discardLast();
  void getMPS(MPS& mps) const {mps=state;}// a copy
  int getIdx(){return spinId;}
  void proposeGlobalMove();
  void discardLastGlobal();
};

// Return a (random) product state from the computational basis (I need it to be  real for the optimized calculation)
void generateProductState(int N,int d,int D,MPS& state);


// Avoid optimization in each step and,instead, proceed as in TEBD, by local truncation
double computeProbCosFast(const MPS& state,int D,double J,double g,double h,double beta,
			  double deltaVar,double dt,double x,double r,double Emin,double Emax);

// Same but optimizing every step
double computeProbCos(const MPS& state,int D,double J,double g,double h,double beta,
		      double deltaVar,double dt,double x,double r,double Emin,double Emax);

// For debugging purposes!
double computeProbCosExact(const MPS& state,int D,double J,double g,double h,double E,
			   double deltaVar,double dt,double x,double r);

// Decide a spin configuration (0 or 1) such that the energy is close to E
void findProductStateWithEnergy(int N,MPS& mps,double E,double J,double g,double h);
//void findProductStateWithEnergy(int N,vector<bool>& spins,double E,const MPO& hamil);

void mapToComputationalBasis(MPS& spins);

// the middle Z, ZZ
void getObservablesZ(const MPS& mps,vector<double>& values);

void getObservablesZall(const MPS& mps,vector<double>& values);

// To get the extreme energies and save (recover) the files
void computeGS(MPS& state,double* ener,const MPO& hamil,int Dgs,const string& mpsdir,
	       const string modelstr,bool maxE=false);
const string mpsfilename(int N,int D,const string& mpsdir,const string& modelstr,
			 bool maxE=false);

double estimateEnergy(double beta,int N,double J,double g,double h,double Emin,double Emax);

int readTmpFile(MPS& prodConf,const string& tmpfile);

int d=2;

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

  int N = props.getIntProperty("L");  
  double J = props.getDoubleProperty("J");  
  double g = props.getDoubleProperty("g");  
  double h = props.getDoubleProperty("h");  
  int D = props.getIntProperty("D");  
  int Dgs = props.getIntProperty("Dgs");
  if(Dgs<=0) Dgs=D;
  int Ns = (int)props.getDoubleProperty("Ns");// nr of samples to take  
  double beta = props.getDoubleProperty("beta");  
  double thresh = props.getDoubleProperty("threshold");
  if(thresh<0) thresh=THRESHOLD;
  double deltaVar = props.getDoubleProperty("deltaVar");  
  double dt = props.getDoubleProperty("deltaTrotter");  
  string outfname=props.getProperty("outputFile");
  bool app=(props.getIntProperty("append")!=0);
  double xval=props.getDoubleProperty("x");  
  double rval=props.getDoubleProperty("r");  
  string outfnameGS=props.getProperty("outputFileGS");
  string mpsDir=props.getProperty("mpsdir"); // for GS and max exctd states
  string tmpFile=props.getProperty("tmpfile"); // for the data that allows to continue the calculation
  
  // String to identify states of this model
  stringstream str;
  str<<"J"<<J<<"_g"<<g<<"_h"<<g;

  //  srandom(17); // for tests
  srandom((int)time(0));

  // I need an estimate of the lowest and highest energies in the spectrum
  double Emin=0.;
  double Emax=0.;
  {
    MPS gs(N,Dgs,d);
    MPS gsNeg(N,Dgs,d); // max. excited state
        
    IsingHamiltonian isingH(N,d,J,g,h);const MPO& mpoH=isingH.getHMPO();
    computeGS(gs,&Emin,mpoH,Dgs,mpsDir,str.str(),0);
    cout<<"Ground state found with eigenvalue "<<Emin
	<<" (Energy per particle="<<Emin/N<<")"<<endl;
    computeGS(gsNeg,&Emax,mpoH,Dgs,mpsDir,str.str(),1);
    cout<<"Highest energy state found with eigenvalue "<<Emax
	<<" (Energy per particle="<<Emax/N<<")"<<endl;

    if(outfnameGS.length()!=0){
      ofstream* outGS;
      if(!file_exists(outfnameGS.data())||!app){
	outGS=new ofstream(outfnameGS);
	*outGS<<"% N\t J\t g\t h\t D\t level(0/1)\t Energy"<<endl;
	*outGS<<setprecision(12);
      }
      else
	outGS=new ofstream(outfnameGS,ios::app);
      if(!outGS->is_open()){
	cout<<"Error: impossible to open file "<<outfnameGS<<
	  " for output"<<endl;
	exit(1);
      }
      *outGS<<setprecision(12);
      *outGS<<N<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<Dgs<<"\t"<<0<<"\t"
	    <<Emin<<endl;
      *outGS<<N<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<Dgs<<"\t"<<1<<"\t"
	    <<Emax<<endl;
      outGS->close();
      delete outGS;
    }
  }
    
  // Estimate which energy to look for for the first state (estimate
  // average energy from Hartmann's slope at beta=0, or the extremes): -beta*sigma0^2
  double E0=estimateEnergy(beta,N,J,g,h,Emin,Emax);

  cout<<"Estimated energy "<<E0<<" for beta="<<beta<<endl;
  
  MPS prodConf(N,1,d);
  findProductStateWithEnergy(N,prodConf,E0,J,g,h);

  mapToComputationalBasis(prodConf);
  // Just for debugging
  IsingHamiltonian isingH(N,d,J,g,h);const MPO& hamil=isingH.getHMPO();
  Contractor& contractor=Contractor::theContractor();
  cout<<"Mapped to computational basis, now energy "<<contractor.contract(prodConf,hamil,prodConf)<<endl;
    
  SamplingConfig& config=SamplingConfig::config();
  config.resetWithMPS(prodConf);

  int cnt=0;
  double Pdelta=computeProbCosFast(prodConf,D,J,g,h,beta,deltaVar,dt,xval,rval,Emin,Emax);
  //double Pdelta=computeProbCos(prodConf,D,J,g,h,beta,deltaVar,dt,xval,rval,Emin,Emax);
  double PdeltaOld=Pdelta;
  
  vector<double> values;
  getObservablesZ(prodConf,values);

  vector<double> results(values);cnt++; // let's take this one too
  //for(int p=0;p<values.size();p++) results[p]*=Pdelta; // for the uniform test!
  //double totNorm=Pdelta; // for the uniform test!!
  
  // Just for debugging small systems:reconstruct the pdf
  vector<int>pdf(pow(2,N),0);
  
  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  ofstream* out;
  //if(!app){
  out=new ofstream(outfname.data());
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<"%cnt\t <ZI> \t <IZ>\t <ZZ(middle)> "<<endl;
  //}
  *out<<cnt<<"\t"<<real(values[0])<<"\t"<<real(values[1])<<"\t"<<real(values[2])<<"\t";
#ifdef DEBUGPDF
  //*out<<totNorm<<"\t";
  for(int k=0;k<pow(2,N);k++) *out<<pdf[k]<<"\t";
#endif
  *out<<endl;
  out->close();delete out;
    
  cout<<"Initial Pdelta="<<Pdelta<<" values "<<values<<" results="<<results<<endl;
  clock_t start=clock();
  clock_t finish;

  bool accept=0;
  while(cnt<Ns){
    cnt++;
    //cout<<"Starting cnt="<<cnt<<endl;
    // Try a move
    //config.proposeGlobalMove();
    config.proposeMove();
    config.getMPS(prodConf);
    Pdelta=computeProbCosFast(prodConf,D,J,g,h,beta,deltaVar,dt,xval,rval,Emin,Emax);
    //Pdelta=computeProbCos(prodConf,D,J,g,h,beta,deltaVar,dt,xval,rval,Emin,Emax);
    //if(0){ // test accept all
    if(Pdelta<thresh){
      Pdelta=0;
      if(PdeltaOld<thresh){
	cout<<"Why am I here??"<<endl;exit(1);
	accept=((double)random()/RAND_MAX)<.5; // from 0 to 1
      }
      else
	accept=0;
    }
    else{
      double rate=Pdelta/PdeltaOld;
      if(rate>1) accept=true;
      else
	accept=((double)random()/RAND_MAX)<rate; // from 0 to 1
    }
    //}
    //  accept=1;totNorm+=Pdelta;// for the uniform test!
    //    cout<<"New Pdelta="<<Pdelta<<" accept="<<accept<<" results="<<results<<", values="<<values<<endl;
    //cout<<"\tPdelta="<<Pdelta<<" accept="<<accept<<endl;
    if(accept){
      config.confirmLast();
      // do operators
      getObservablesZ(prodConf,values);
      PdeltaOld=Pdelta;
    }
    else{
      //config.discardLastGlobal();
      config.discardLast();
      // do operators: nothing changes
    }
    // Update the pdf
#ifdef DEBUGPDF
    pdf[config.getIdx()]++;
#endif
    for(int p=0;p<values.size();p++){
      results[p]+=values[p]; // for importance sampling!
      //results[p]+=values[p]*Pdelta; // for the uniform test!
    }
    //cout<<" cnt="<<cnt<<" values="<<values<<", results="<<results<<endl;
    if(cnt%100==0){ // write results to file
      //finish=clock();
      //cout<<"Time for cnt="<<cnt<<" = "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;  
      out=new ofstream(outfname.data(),ios::app);
      *out<<cnt<<"\t"<<real(results[0])<<"\t"<<real(results[1])<<"\t"<<real(results[2])<<"\t";
      //*out<<totNorm<<"\t";
#ifdef DEBUGPDF
      for(int k=0;k<pow(2,N);k++) *out<<pdf[k]<<"\t";
#endif
      *out<<endl;
      out->close();delete out;
      //start=clock();
    }
  }  

  cout<<"At the end, cnt="<<cnt<<", values="<<results<<endl;
  for(int k=0;k<results.size();k++){
    cout<<results[k]/cnt<<"\t";
    //cout<<results[k]/totNorm<<"\t";
  }
  cout<<endl;
  //cout<<"normalization: "<<totNorm<<endl;
  
}


void findProductStateWithEnergy(int N,MPS& spins,double E,double J,double g,double h){
  // I could easily find a general product state, but I want a real one, and additionally
  // one in the computational basis
  // I take the MPS solution for prod states and the select the largest overlap site by site
  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  IsingHamiltonian isingH(N,d,J,g,h,-E); // H-E
  const MPO& hamil=isingH.getHMPO();
  const MPO* ptrs[]={&hamil,&hamil};
  MPO ham2(N);
  MPO::join(2,ptrs,ham2);// MPO for (H-E)^2
  spins=MPS(N,1,d);double value=3;
  contractor.findGroundState(ham2,1,&value,spins);
  // this should indeed be real =>could actually use to define a local
  // basis that substitutes the computational one, and on each site,
  // there will be a different direction z and x, such taht this is
  // the 000.. state, and flipping means rotating to the 1, also real. TODO!!!!
  cout<<"Found a product state with energy "<<real(contractor.contract(spins,hamil,spins))+E<<endl;
}

void mapToComputationalBasis(MPS& spins){
  mwArray v0(Indices(d,1,1));
  mwArray v1(v0);
  v0.setElement(ONE_c,Indices(0,0,0));
  v1.setElement(ONE_c,Indices(1,0,0));
  for(int k=0;k<spins.getLength();k++){
    mwArray aux=spins.getA(k).getA();
    complex_t c0=aux.getElement(Indices(0,0,0));
    complex_t c1=aux.getElement(Indices(1,0,0));
    if(abs(c0)>abs(c1))
      spins.replaceSite(k,v0,0);
    else
      spins.replaceSite(k,v1,0);
  }
}

void getObservablesZ(const MPS& mps,vector<double>& values){
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  static mwArray sigZ(Indices(d,d),dataz);    
  int N=mps.getLength();
  values.clear();
  vector<int> Ps;Ps.push_back(N/2-1);Ps.push_back(N/2);
  for(int pos=0;pos<Ps.size();pos++){
    mwArray Aaux=mps.getA(Ps[pos]).getA(); // they are products
    Aaux.reshape(Indices(d,1));
    values.push_back(real((Hconjugate(Aaux)*sigZ*Aaux).getElement(Indices(0,0))));
  }    
  values.push_back(values[0]*values[1]);
  //cout<<"values="<<values<<endl;
}

void getObservablesZall(const MPS& mps,vector<double>& values){
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  static mwArray sigZ(Indices(d,d),dataz);    
  int N=mps.getLength();
  values.clear();
  for(int k=0;k<N;k++){
    mwArray Aaux=mps.getA(k).getA(); // they are products
    Aaux.reshape(Indices(d,1));
    values.push_back(real((Hconjugate(Aaux)*sigZ*Aaux).getElement(Indices(0,0))));
  }
}

void generateProductState(int N,int d,int D,MPS& state){
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  static mwArray sigX(Indices(d,d),datax);    
  state=MPS(N,D,d); // this is the 0000
  cout<<"State: ";
  for(int k=0;k<N;k++){
    //int r=(rand()%N ); // fromm 0 to N
    double r=(double)random()/RAND_MAX; // from 0 to 1
    if(r<.5){ // flip the site
      state.applyLocalOperator(k,sigX); // flipped to 1
      cout<<1;
    }
    else{
      cout<<0;
    }
  }
  cout<<endl;  
}

double coefficientC0log(int M){
  double logc0=0.;
  for(int k=1;k<=M;k++){
    logc0+=log2(k);
    if(k<=M/2) logc0-=2*log2(k);
  }
  logc0-=M;
  return logc0;
}

void prepareCoefficientsCm(vector<double>& cmvec,int M,int R){
  double nextC=0;  // log of component 0
  for(int m=1;m<=R;m++){
    //cmvec.push_back(pow(2,nextC));
    cmvec.push_back(nextC);
    nextC+=log2(M-2*(m-1))-log2(M+2*m);
  }
  //cmvec.push_back(pow(2,nextC)); //the R-th one
  cmvec.push_back(nextC); //the R-th one
  //cout<<"Coeffs cm (M="<<M<<"):"<<cmvec<<endl;
}

void prepareCoefficientsBm(vector<complex_t>& bmvec,double beta,double Emin,double Emax,int N,double r,int R){
  bmvec.clear();
  //cout<<"prepareCoefficientsBm, beta="<<beta<<", Emin="<<Emin<<", Emax="<<Emax<<", N="<<N<<", r="<<r<<", R="<<R<<endl;
  for(int m=0;m<=R;m++){
    complex_t expo=beta*ONE_c+(2*m/(r*sqrt(N)))*I_c;
    bmvec.push_back((exp(-Emin*expo)-exp(-Emax*expo))/expo);
    //cout<<"For m="<<m<<", beta+2mi/r sqrt(N)="<<expo<<",  exp(-Emin*this)="<<exp(-Emin*expo)
    //	<<",  exp(-Emax*this)="<<exp(-Emax*expo)<<", 1/this="<<ONE_c/expo<<endl;
  }
  cout<<"Coeffs bm (R="<<R<<"):"<<bmvec<<endl;//exit(1);
}

double computeProbCosExact(const MPS& state,int D,double J,double g,double h,double E,
			   double deltaVar,double dt,double x,double r){
  static bool init=0;
  int dim=pow(2,state.getLength());
  static mwArray pDelta(Indices(dim,dim));
  if(!init){
    int N=state.getLength();
    IsingHamiltonian isingH(N,d,J,g,h,-E); // I put the scale in the time
    const MPO& hamil=isingH.getHMPO();
    mwArray fullH;
    int Mr=ceil(r*r*N/(deltaVar*deltaVar)); //  this is M
  //  should be even
    Mr=Mr+(Mr%2); 
    expandOper(hamil,fullH);
    mwArray U;vector<complex_t> D;
    wrapper::eig(fullH,D,U,true);
    //cout<<"Mr="<<Mr<<", Eigs: "<<D<<endl;
    for(int k=0;k<D.size();k++){
      pDelta.setElement(pow(cos(real(D[k])/(r*sqrt(N))),Mr)*ONE_c,Indices(k,k));
    }
    pDelta.multiplyLeft(U);U.Hconjugate();pDelta.multiplyRight(U);
    init=true;
  }
  mwArray fullV;
  expandVec(state,fullV);
  return real((Hconjugate(fullV)*pDelta*fullV).getElement(0));
}


double computeProbCos(const MPS& state,int D,double J,double g,double h,double beta,
			  double deltaVar,double dt,double x,double r,double Emin,double Emax){
  Contractor& contractor=Contractor::theContractor();
  int N=state.getLength();
  //Determine the numberof steps, etc
  int M=ceil(r*r*N/(deltaVar*deltaVar)); //  this is M
  //  should be even
  M=M+(M%2); 
  int R=min(ceil(x*r*sqrt(N)/(deltaVar)),floor(M/2)); //max step nr =x*sqrt(M)
  static vector<double> cmvec;
  static vector<complex_t> bmvec;
  static bool init=false;
  if(!init){
    cout<<"First call to computeProbCosFast, Mr="<<M<<", R="<<R<<endl;
    prepareCoefficientsCm(cmvec,M,R);
    prepareCoefficientsBm(bmvec,beta,Emin,Emax,N,r,R);
    init=true;
    cout<<"Prepared cm coefficients "<<cmvec<<endl;
    cout<<"Prepared bm coefficients "<<bmvec<<endl;
  }
  
  // the rescaled Hamiltonian
  double scale=r*sqrt(N);
  IsingHamiltonian isingH(N,d,J,g,h); // I put the scale in the time
  MPO Uevol(N);
  double stepT=1./scale;
  int nrTr=ceil(stepT/dt);
  // cout<<"Original Trotter step "<<dt<<" interval="<<stepT<<" (r*sqrt(N)="<<scale
  //     <<"), step-nr="<<nrTr;
  dt=stepT/nrTr; //correct
  //  cout<<"(corrected dt="<<dt<<"), M="<<M<<" R="<<R<<endl;
  isingH.getUMPO(Uevol,dt*I_c,0.); // evolution operator  
  vector<complex_t> aMs; // store <p|exp(2*m*i*H/scale)|p> for m=0 until R
  MPS aux(state); // the one I go evolving
  aMs.push_back(contractor.contract(aux,state)); // should be 1
  // I need an auxiliary one to apply the MPO
  MPS placeholder(aux);
  MPS* next=&aux;MPS* old=&placeholder;
  for(int m=1;m<=R;m++){
    next=&aux;old=&placeholder;
    for(int k=0;k<nrTr;k++){
      if(k>0){MPS* tmp=old;old=next;next=tmp;} // exchange places
      contractor.optimize(Uevol,*old,*next,D); // last one is in next
      //next->approximate(Uevol,*old,D);
    }
    if(next!=&aux){aux=*next;}
    else{placeholder=*next;} // both have the latest one, now
    placeholder.conjugateMPS();
    aMs.push_back(contractor.contract(aux,placeholder)); // <exp(2miH/scale)> 
    placeholder.conjugateMPS(); // Return to previous
  }

  // Now should have all the vectors => combine for Pdelta (except c0 coeff)
  double Pdelta=real(bmvec[0]);
  for(int m=1;m<=R;m++){
    Pdelta+=pow(2,cmvec[m]+1)*real(bmvec[m]*aMs[m]);
  }
  double c0=coefficientC0log(M);
  return pow(2,c0)*Pdelta;
}


double computeProbCosFast(const MPS& state,int D,double J,double g,double h,double beta,
			  double deltaVar,double dt,double x,double r,double Emin,double Emax){
  Contractor& contractor=Contractor::theContractor();
  int N=state.getLength();
  //Determine the numberof steps, etc
  int M=ceil(r*r*N/(deltaVar*deltaVar)); //  this is M
  //  should be even
  M=M+(M%2); 
  int R=min(ceil(x*r*sqrt(N)/(deltaVar)),floor(M/2)); //max step nr =x*sqrt(M)
  static vector<double> cmvec;
  static vector<complex_t> bmvec;
  static bool init=false;
  if(!init){
    cout<<"First call to computeProbCosFast, Mr="<<M<<", R="<<R<<endl;
    prepareCoefficientsCm(cmvec,M,R);
    prepareCoefficientsBm(bmvec,beta,Emin,Emax,N,r,R);
    init=true;
    cout<<"Prepared cm coefficients "<<cmvec<<endl;
    cout<<"Prepared bm coefficients "<<bmvec<<endl;
  }
  
  // the rescaled Hamiltonian
  double scale=r*sqrt(N);
  IsingHamiltonian isingH(N,d,J,g,h); // I put the scale in the time
  MPO Uevol(N);
  double stepT=1./scale;
  int nrTr=ceil(stepT/dt);
  // cout<<"Original Trotter step "<<dt<<" interval="<<stepT<<" (r*sqrt(N)="<<scale
  //     <<"), step-nr="<<nrTr;
  dt=stepT/nrTr; //correct
  //  cout<<"(corrected dt="<<dt<<"), M="<<M<<" R="<<R<<endl;
  isingH.getUMPO(Uevol,dt*I_c,0.); // evolution operator  
  vector<complex_t> aMs; // store <p|exp(2*m*i*H/scale)|p> for m=0 until R
  MPS aux(state); // the one I go evolving
  aMs.push_back(contractor.contract(aux,state)); // should be 1
  // I need an auxiliary one to apply the MPO
  MPS placeholder(aux);
  MPS* next=&aux;MPS* old=&placeholder;
  for(int m=1;m<=R;m++){
    next=&aux;old=&placeholder;
    for(int k=0;k<nrTr;k++){
      if(k>0){MPS* tmp=old;old=next;next=tmp;} // exchange places
      //contractor.optimize(Uevol,*old,*next,D); // last one is in next
      next->approximate(Uevol,*old,D);
    }
    if(next!=&aux){aux=*next;}
    else{placeholder=*next;} // both have the latest one, now
    placeholder.conjugateMPS();
    aMs.push_back(contractor.contract(aux,placeholder)); // <exp(2miH/scale)> 
    placeholder.conjugateMPS(); // Return to previous
  }

  // Now should have all the vectors => combine for Pdelta (except c0 coeff)
  double Pdelta=real(bmvec[0]);
  for(int m=1;m<=R;m++){
    Pdelta+=pow(2,cmvec[m]+1)*real(bmvec[m]*aMs[m]);
  }
  double c0=coefficientC0log(M);
  return pow(2,c0)*Pdelta;
}


//private:
SamplingConfig::SamplingConfig():state(1,1,d),moved(-1),spinId(0){
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  sigX=mwArray(Indices(d,d),dataX);
  cout<<"Created default config"<<endl;
  oldSpinId=0;
}

//public:

SamplingConfig& SamplingConfig::config(){
  static SamplingConfig theConfig;
  return theConfig;
} 


void printZconfig(const MPS& state){
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  static mwArray sigZ=mwArray(Indices(d,d),dataZ);
  MPS aux(state);Contractor& contractor=Contractor::theContractor();
  for(int k=0;k<state.getLength();k++){
    aux.applyLocalOperator(k,sigZ);
    complex_t val=contractor.contract(aux,state);
    cout<<(real(val)<0);
    aux.applyLocalOperator(k,sigZ);
  }
  //cout<<endl;
}

int getZconfig(const MPS& state){
  int spinId=0;
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  static mwArray sigZ=mwArray(Indices(d,d),dataZ);
  MPS aux(state);Contractor& contractor=Contractor::theContractor();
  int N=state.getLength();
  for(int k=0;k<N;k++){
    aux.applyLocalOperator(k,sigZ);
    complex_t val=contractor.contract(aux,state);
    if(real(val)<0){
      spinId=(spinId|(1<<(N-k-1)));
      //cout<<"\t Site  "<<k<<" at 1 , now spinId="<<spinId<<endl;
    }
    aux.applyLocalOperator(k,sigZ);
  }
  return spinId;
}

void SamplingConfig::resetWithMPS(const MPS& mps){
  state=mps; // copy it
  moved=-1;
  spinId=getZconfig(mps);
  cout<<"Reset config with MPS("<<spinId<<"):";printZconfig(state);cout<<endl;
  oldSpinId=spinId;
}

// void SamplingConfig::prepareMPS(){
//   int N=spins.size();
//   if(state!=NULL) delete state;
//   state=new MPS(N,1,d);
//   for(int k=0;k<N;k++){
//     if(spins[k]){ // flip the site to 1
//       state.applyLocalOperator(k,sigX); // flipped to 1
//     }
//   }
// }

// propose a random change in configuration
void SamplingConfig::proposeMove(){
  oldSpinId=spinId;
  int N=state.getLength();
  moved=(random()% N); // from 0 to N
  state.applyLocalOperator(moved,sigX); // flipped
  spinId^=(1<<(N-moved-1));
#ifdef DEBUGINFO
  cout<<"\tProposing move to("<<spinId<<"):";printZconfig(state);cout<<endl;
#endif
}

// // propose a random change in configuration
void SamplingConfig::proposeGlobalMove(){
  //  To test: pick a state at random
  oldSpinId=spinId;
  int N=state.getLength();
  while(spinId==oldSpinId) spinId=(random()%(long)pow(2,N));
  moved=oldSpinId^spinId; // the changed ones  
  for(int k=0;k<N;k++){
    if(moved&(1<<(N-k-1))) state.applyLocalOperator(k,sigX); // flipped
  }
#ifdef DEBUGINFO
  cout<<"\tProposing move to("<<spinId<<"):";printZconfig(state);cout<<endl;
#endif
}


void SamplingConfig::confirmLast(){
  moved=-1;
#ifdef DEBUGINFO
  cout<<"Accepted move to("<<spinId<<"):";printZconfig(state);
#endif
}

void SamplingConfig::discardLast(){
  int N=state.getLength();
  state.applyLocalOperator(moved,sigX); // flip it back
  spinId^=(1<<(N-moved-1));
  moved=-1;
#ifdef DEBUGINFO
  cout<<"Rejected, back to("<<spinId<<"):";printZconfig(state);
#endif
}

void SamplingConfig::discardLastGlobal(){
  int N=state.getLength();
  for(int k=0;k<N;k++){
    if(moved&(1<<(N-k-1))) state.applyLocalOperator(k,sigX); // flip back
  }  
  spinId=oldSpinId;
  moved=-1;
#ifdef DEBUGINFO
  cout<<"Rejected, back to("<<spinId<<"):";printZconfig(state);
#endif
}

double estimateEnergy(double beta,int N,double J,double g,double h,double Emin,double Emax){
  double  sigma02=J*J*(N-1)+g*g*N+h*h*N;
  double E0=-beta*sigma02;
  //  double xmin=(Emin+beta*sigma02)/(sqrt(2*sigma02));
  //double xmax=(Emax+beta*sigma02)/(sqrt(2*sigma02));
  // double integr=int_xmin^xmax exp(-x^2) dx; // how do I evaluate here??
  // E0-=sqrt(sigma02/2)*(exp(-xmax*xmax)-exp(-xmin*xmin))/(integr);
  if(E0<Emin)  return Emin;
  if(E0>Emax) return Emax;
  return E0;
}

const string mpsfilename(int N,int D,const string& mpsdir,const string& modelstr,bool maxE){
  stringstream str;
  int level=maxE?1:0;
  str<<mpsdir<<"/MPS_L"<<N<<"_"<<modelstr<<"_l"<<level<<"_D"<<D;
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return str.str();
}

void computeGS(MPS& state,double* ener,const MPO& hamil,int Dgs,const string& mpsdir,const string modelstr,bool maxE){
  int N=hamil.getLength();
  Contractor& contractor=Contractor::theContractor();
  MPO hamilNeg(N);
  if(maxE){
    // need to change sign of H
    hamilNeg.setOp(0,new Operator(-1.*hamil.getOp(0).getFullData()),true);
    for(int k=1;k<N;k++)
      hamilNeg.setOp(k,&hamil.getOp(k),0);
  }
  bool found=0;
  // Now try to check for the MPS file
  const string mpsfile=mpsfilename(N,Dgs,mpsdir,modelstr,maxE);
  if(file_exists(mpsfile)){
    state.importMPS(mpsfile.data());
    cout<<"Recovered existing file for D="<<Dgs<<"=> no optimization"<<endl;
    state.gaugeCond('R',1);
    found=1;
    // compute the energy and return  
    if(maxE){
      *ener=real(contractor.contract(state,hamil,state));
    }
    else{
      *ener=real(contractor.contract(state,hamil,state));
    }
    return;
  }
  else{
    // try any smaller D
    int D_=Dgs-1;
    while(!found&&D_>0){ 
      const string mpsfileD=mpsfilename(N,D_,mpsdir,modelstr,maxE);
      if(file_exists(mpsfileD)){
	state.importMPS(mpsfileD.data());
	cout<<"Imported file for D="<<D_<<endl;
	found=1;
      }
      else D_--;
    }
    
    if(!found){
      cout<<"No file found for "<<(!maxE)<<"; Start with random state."<<endl;
      state.setRandomState(); // the intial state, random
    }
    state.gaugeCond('R',1);
    state.gaugeCond('L',1);
    // Now try to find the ground state
    if(maxE){
      contractor.findGroundState(hamilNeg,Dgs,ener,state);
      *ener=-*ener;
    }
    else{
      contractor.findGroundState(hamil,Dgs,ener,state);      
    }
    // Now save it
    state.exportMPS(mpsfile.data());
  }
    
}


void writeTmpFile(int cnt,const MPS& prodConf,const string& tmpfile){
  // Trick: save the cnt in the norm factor
  MPS
  int cnt=0;
  // Read
}


int readTmpFile(MPS& prodConf,const string& tmpfile){
  if(tmpfile.length()==0) return 0;
  int cnt=0;
  // Read
}
