
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#ifndef NEWIMPLEM
// To run in Mac Os X
#ifdef __APPLE_CC__
#include <CoreFoundation/CoreFoundation.h>
#endif
#include "libdiminf.h"
#include "diminf.h"

typedef OperatorRow MPO;
#define CHECKDIMS 1

#else
#include "SchwingerHamiltonian.h"
#endif

#include <vector>
using namespace std;
using namespace shrt;

#ifdef MACOSX
#define LOCALDIR "/Users/banulsm/Desktop/SVN/projects/mpsdyn/src/bin/"
#else
#define LOCALDIR "/home/banulsm/ALaPorra/bin/"
#endif

/** SchwingerPD runs the findGroundState routine with the MPO for the
    Schwinger Hamiltonian \ref <SchwingerHamiltonian>, and saves
    values of order parameters, to construct a phase diagram. 

    To be used as test in the new implementation, compile make test9

    Receives arguments:
    \param <L> (int) length of the chain
    \param <mg> (double) parameter \f$m/g\f$ of the Hamiltonian
    \param <x> (double) parameter \f$x\f$ of the Hamiltonian
    \param <L0> (double) parameter \f$l_0\f$ of the Hamiltonian
    \param <alpha> (double) parameter \f$\alpha\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)

    For it to work in Mac Os X, a thread has to be created which
    receives the main function as argument.
*/

int runSchwingerPD(int argc,const char* argv[]);

int main(int argc,const char* argv[]){
#if (!defined NEWIMPLEM)&&(defined __APPLE_CC__)
  //#ifdef __APPLE_CC__
  cout<<"Calling runMain inside the loop for Mac OS X"<<endl;
  mclmcrInitialize();
  return mclRunMain((mclMainFcnType)runSchwinger,argc,argv);
#else
  runSchwingerPD(argc,argv);
#endif
}

#include "time.h"

/** An auxiliary function which computes the coefficients of the given MPS
    in the computational basis. It takes as argument the state, and the 
    coefficient to be computed, as an integer from 0 to 2^L-1, and returns 
    the complex coefficient.
*/
complex_t computeCoefficient(const MPS& state,long int elem);
complex_t computeCoefficient(const MPS& state,const vector<int>& binZ);
vector<int> intToBin(long int z,int L);
long int binToLong(const vector<int>& refStr);
void getSzMPO(int L,MPO& Smpo);

int runSchwingerPD(int argc,const char* argv[]){

  char filename[120];
  int error=0;
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double mg=atof(argv[++cntr]);
  double x=atof(argv[++cntr]);
  double L0=atof(argv[++cntr]);
  double alpha=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  double mu=2*mg*sqrt(x);
  // WARNING!!! If x=0, this won't work, so then take mu=mg, but this is not right either!
  if(x==0.) {
    cout<<"WARNING!!!! The coefficients in case x=0 can be a bt funny"<<endl;
    mu=mg;
  }

  cout<<"Initialized arguments: L="<<L
      <<", mg="<<mg
      <<", mu="<<mu
      <<", x="<<x
      <<", l0="<<L0
      <<", alpha="<<alpha
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<setprecision(10)<<"# N\t x\t mg\t alpha\t l0\t D\t Energy \t Energy/(2Nx)\t <Psi-bar Psi>/g\t GammaA \t Gamma5"<<endl;
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  
#ifndef NEWIMPLEM
  initializelib(&error);
  try{
    initGlobals(mwArray(CHECKDIMS),mwArray(SEEDVAL));
  }
  catch(mwException& e){
    cout<<"Exception caught: "<<e.what()<<"in timedep"<<endl;
    exit(212);
  }
  cout<<"Initialized global constants"<<endl;
  try{
#endif
    Contractor& contractor=Contractor::theContractor();
    //    contractor.setConvTol(1E-8);
    cout<<"Initialized Contractor?"<<endl;
    int d=2;
    MPS init(L,D,d);
    MPS gs(L,D,d);
    //    cout<<"Initialized empty MPS: "<<init<<endl;
    gs.setRandomState(); // the intial state, random
    gs.gaugeCond('R',1);
    gs.gaugeCond('L',1);
    //cout<<"Set random MPS: "<<init<<endl;
    
    // First: put H and find the GS
    int Dop=5; // bond dimension of Schwinger H
#ifdef NEWIMPLEM
    double offset=0,nu=0,y=1.33;
    SchwingerHamiltonian hSch(L,mu,x,L0,alpha,offset,nu,y);
    const MPO& hamil=hSch.getHMPO();
    //    hamil.exportForMatlab("schwHamil.m");
    // exit(1);
#else
    MPO hamil(L,Dop,d); // it is empty    
    mwArray mpo;
    constructSchwingerHMPO(1,mpo,mwArray(L*1.0),mwArray(mu*1.0),
			   mwArray(x*1.0),mwArray(L0*1.0),mwArray(alpha*1.0));
    hamil.setOperatorArray(mpo);
#endif    
    
    //hamil.saveOperatorArray("/home/banulsm/SVN/peps/projects/folding/src/tests/hmpo.mat");
    double lambda=0.;
    clock_t start=clock();
    contractor.findGroundState(hamil,D,&lambda,gs);
    cout<<"Ground state found with eigenvalue "<<lambda<<endl;
    //cout<<"Checking: <H(0)>="<<contractor.contract(init,hamil,init)<<endl;
    //    strcpy(filename,LOCALDIR);
    // gs.exportForMatlab(strcat(filename,"mpsGS_N6.m"));
    // hamil.exportForMatlab(strcat(filename,"hamilN6.m"));

    clock_t final=clock();
    cout<<"Total time in findGroundState="<<(final-start)*1./CLOCKS_PER_SEC<<endl;
    cout<<"And expectation value of H="<<contractor.contract(gs,hamil,gs)<<endl;
    cout<<"Norm="<<contractor.contract(gs,gs)<<endl;

    bool printing=true;
    if(printing){
      // Use the MPOs for Gamma^alpha and Gamma^5
#ifdef NEWIMPLEM
      MPO GammaA(L);
      MPO Gamma5(L);
      MPO fCond(L);
      hSch.constructGammaAlphaMPO(GammaA);
      hSch.constructGamma5MPO(Gamma5);
      hSch.constructCondensateMPO(fCond);
      complex_t gammaAlpha=contractor.contract(gs,GammaA,gs);
      complex_t gamma5=contractor.contract(gs,Gamma5,gs);
      complex_t condensate=contractor.contract(gs,fCond,gs);
      cout<<"Condensate expectation value found to be "<<condensate<<endl;
      complex_t gA=gammaAlpha+(complex_t){L0+alpha,0.};
      complex_t g5=(sqrt(x)/L)*gamma5;
#else      
      MPO GammaA(L,3,d);
      MPO Gamma5(L,3,d);
      constructSchwingerGammaA(1,mpo,mwArray(L*1.0));
      GammaA.setOperatorArray(mpo);
      constructSchwingerGamma5(1,mpo,mwArray(L*1.0));
      Gamma5.setOperatorArray(mpo);
      // Now compute and save all the expectation values
      mwArray gammaAlpha=contractor.contract(gs,GammaA,gs);
      mwArray gamma5=contractor.contract(gs,Gamma5,gs);
      mwArray g5;
      mwArray coeff(sqrt(x)/L);
      m_times(1,g5,gamma5,coeff);
      mwArray gA;
      m_add(1,gA,gammaAlpha,mwArray(L0+alpha));
#endif                 
      *out<<setprecision(10)<<L<<"\t"<<x<<"\t"<<mg<<"\t"<<alpha<<"\t"<<L0<<"\t"
	  <<D<<"\t"<<lambda<<"\t"<<lambda/(2*L*x)<<"\t";
      *out<<condensate<<"\t";
      *out<<gA<<"\t";
      *out<<g5<<"\t";
      *out<<endl;


#ifdef COEFFCOMP
    // Now also check the expectation value of Sz and how close it is to an eigenvector
    MPO Sz(L);
    getSzMPO(L,Sz);
    //hamil.exportForMatlab("hamil.m");
    //Sz.exportForMatlab("Sz.m");
    complex_t sz=contractor.contract(gs,Sz,gs);
    MPO Sz2(L);
    const MPO* ptrs[2]={&Sz,&Sz};
    MPO::join(2,ptrs,Sz2);
    complex_t sz2=contractor.contract(gs,Sz2,gs);
    cout<<"<Sz>="<<sz<<", <Sz^2>="<<sz2<<", <Sz^2>-<Sz>^2="<<sz2-sz*sz<<endl;
    *out<<"% <Sz>\t <Sz^2>-<Sz>^2"<<endl;
    *out<<sz<<"\t"<<sz2-sz*sz<<endl<<endl;

      // In this new version, the coefficients of the GS in the computational basis are computed one by one,
      // if the chain is up to 20 sites long
      *out<<"% Basis element \t Coefficient"<<endl;
      // I want to use as reference the |0> state from the SC expansion, which in my case is 01010101....
      complex_t reference=ONE_c; // to check signs, I rephase them all by dividing by the first non-zero one
      vector<int> scGS(L,0); for(int k=1;k<L;k=k+2) scGS[k]=1;
      reference=computeCoefficient(gs,scGS);
      reference=(1/abs(reference))*conjugate(reference); // the phase I want to remove
      if(L<=20){
	for(long int z=0;z<pow(2,L);z++){
	  vector<int> binZ=intToBin(z,L);
	  for(int k=0;k<L;k++) *out<<binZ[k];
	  complex_t coeff=computeCoefficient(gs,z);
	  if(abs(coeff)>1E-12){
	    coeff=coeff*reference;
	  }
	  else coeff=ZERO_c;
	  int sumZ=0; for(int k=0;k<L;k++) sumZ+=pow(-1,binZ[k]);
	  if(sumZ==0) *out<<"*";
	  *out<<"\t"<<coeff<<endl;
	}
      }

#endif  // COEFFCOMP


    } // end of if(printing)
#ifndef NEWIMPLEM
  }
  catch(mwException& e){
    cout<<"Exception caught: "<<e.what()<<"in timedep"<<endl;
    exit(212);
  }
#endif
  out->close();
  delete out;
#ifndef NEWIMPLEM
  closelib();
#endif
return error;

}

complex_t computeCoefficient(const MPS& state,long int elem){
  int L=state.getLength();
  vector<int> binZ=intToBin(elem,L);
  return computeCoefficient(state,binZ);
}
complex_t computeCoefficient(const MPS& state,const vector<int>& binZ){
  int L=state.getLength();
  complex_t data0[]={ONE_c,ZERO_c};
  complex_t data1[]={ZERO_c,ONE_c};
  static mwArray A0(Indices(2,1,1),data0);
  static mwArray A1(Indices(2,1,1),data1);
  MPS basic(L,1,2);
  for(int k=0;k<L;k++){
    if(binZ[k]==0)
      basic.setA(k,A0);
    else if(binZ[k]==1)
      basic.setA(k,A1);
    else{
      cout<<"Error: binary sequence "<<binZ<<" not valid at pos "<<k+1<<endl;
      exit(1);
    }
  }
  Contractor& contractor=Contractor::theContractor();
  return contractor.contract(state,basic);
}



vector<int> intToBin(long int z,int L){
  long int tmp=z;
  vector<int> result(L,0);
  for(int k=L-1;k>=0;k--){
    result[k]=(int)(tmp%2);
    tmp=(tmp-result[k])/2;
  }
  return result;
}


long int binToLong(const vector<int>& refStr){
  int L=refStr.size();
  long int z=0;
  for(int k=0;k<L;k++)
    z+=refStr[k]*pow(2,L-k-1);
  return z;
}


void getSzMPO(int L,MPO& Smpo){
  Smpo.initLength(L);
  int d=2;
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  mwArray idOp=identityMatrix(d);
  mwArray spinOp(Indices(d,d),dataZ);

  mwArray Z=mwArray(Indices(2,d,d));
  Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(idOp.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(spinOp.getElement(Indices(i1,i2)),Indices(1,i1,i2));
    }
  Z.reshape(Indices(2,d*d));

  int D=2;
  mwArray C(Indices(D,D,2));
  mwArray Cl(Indices(1,D,2)),Cr(Indices(D,1,2));
  Cl.setElement(ONE_c,Indices(0,0,0));
  Cr.setElement(ONE_c,Indices(D-1,0,0));
  C.setElement(ONE_c,Indices(0,0,0));
  C.setElement(ONE_c,Indices(D-1,D-1,0));
  C.setElement(ONE_c,Indices(0,D-1,1));
  Cl.setElement(ONE_c,Indices(0,D-1,1));
  Cr.setElement(ONE_c,Indices(0,0,1));
  C.reshape(Indices(D*D,2));
  Cl.reshape(Indices(D,2));  Cr.reshape(Indices(D,2));
  mwArray term=C*Z;
  term.reshape(Indices(D,D,d,d));term.permute(Indices(3,1,4,2));
  mwArray termL=Cl*Z;
  termL.reshape(Indices(1,D,d,d));termL.permute(Indices(3,1,4,2));
  mwArray termR=Cr*Z;
  termR.reshape(Indices(D,1,d,d));termR.permute(Indices(3,1,4,2));

  Smpo.setOp(0,new Operator(termL),true);
  Smpo.setOp(L-1,new Operator(termR),true);
  if(L>2){
    Smpo.setOp(1,new Operator(term),true);
    for(int k=2;k<L-1;k++){
      Smpo.setOp(k,&Smpo.getOp(1),false);
    }
  }
}
