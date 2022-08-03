#include <math.h>
#include <iomanip>
#include <sstream>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "Properties.h"
#include "time.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

//#define SYMM 1 // If SYMM=1, only left part is computed explicitly, and right is just the reflection

#include "IsingHamiltonian.h"

using namespace shrt;

/** To try to understand what happens when slower particles are in the
    chain, I implement a different contraction strategy, in which time
    and space steps are alternated in a certain manner.  

    I work with pure states here (although it would be generalizable
    to mixed ones).

  */


void constructInitialProductState(int d,int D,int initSt,mwArray& W);
int getInitialSize(int M,int Ls,int Nx);
int getNrTSteps(int m);
int getInitialSizeQuadratic(int M,int Ls,int Nx);
int getNrTStepsQuadratic(int m);
void saveObservables(const MPS& mps,int cnt,double delta,const string& outfname);
// Prepare a product of two single site ops to be contracted between rho^1/2
void prepareOperatorProduct(mwArray& result,const mwArray& opA,const mwArray& opB);



int d=2;
double tolSVD=1E-12;      


int main(int argc,const char* argv[]){
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

  double J_=props.getDoubleProperty("J");
  double h_=props.getDoubleProperty("h");
  double g_=props.getDoubleProperty("g");
  int initSt=props.getIntProperty("initSt");
  double delta=props.getDoubleProperty("delta");
  int M=props.getIntProperty("M"); // max nr of time steps to be simulated
  int Nx=props.getIntProperty("Nx"); // nr of sites absorbed in env by X step
  if(Nx<1) Nx=1; // default is one
  int Ls=props.getIntProperty("Ls"); // sites to keep in the system at the end
  if(Ls<0) Ls=2; // absolute min
  int rate=props.getIntProperty("rate");
  if(rate<=0) rate=1; // default frequency
  int Dcut=props.getIntProperty("D"); // the truncation  (which is the dimension of the effective boundary site)
  //  int Dcontr=props.getIntProperty("Drho"); // the dimension used for the effective rho+boundary as it evolves
  string outfname=props.getProperty("outputfile");
  int savingfreq=props.getIntProperty("savingfreq"); // frequency with
						     // which
						     // intermediate
						     // results are
						     // saved to disk
  bool savingTmp=savingfreq>0;
  // string tmpdir=props.getProperty("tmpdir"); // directory to save
  // 					     // temporary tensors
  // 					     // (big!)
  // if(savingTmp&&tmpdir.empty()){
  //   cout<<"ERROR! No tmpdir specified to save temporary data"<<endl;
  //   exit(1);
  // }

  // 1) Prepare initial size. I assume it is a product, so the boundary D at the beginning is 1
  int L0=getInitialSize(M,Ls,Nx);

  MPS mps(L0,1,d);
  {mwArray A;
    constructInitialProductState(d,1,initSt,A);
    for(int k=0;k<L0;k++) mps.replaceSite(k,A,false);
  }
  mps.gaugeCond('R',1);mps.gaugeCond('L');
  
  // My Hamiltonian
  IsingHamiltonian hamH(3,d,J_,g_,h_);
  MPO Uevol3(3); // to hold the different operators
  {
    MPO Uevol3_single(3);
    hamH.getUMPO(Uevol3_single,-delta*I_c,0); // the single layer normal evolution
    // But I need to double it
    doubleMPO(Uevol3_single,Uevol3,true);
    cout<<"Created the double Uevol"<<endl;
  }
  
  Contractor& contractor=Contractor::theContractor();
  
  int cnt=0; // total nr of steps (real time)
  int cntX=0; // nr of times I've folded in boundary sites (evol in x)
  int Lm=L0; // current size of the MPS (incl. boundary)
  int Db=1; // dimension of the virtual leg at the boundary, representing the environment
  
  while(cnt<M){
    // How many normal time steps do I need to apply for this stage?
    int m=getNrTSteps(cntX);
    //    m=1; // trick to do linear
    // //    m=M; // trick to keep hard boundaries!!
    // Prepare the evolution MPO (which is the normal one, but with
    // the identity on the extra legs at the boundaries)
    MPO evolU(Lm);
    mwArray idDb=identityMatrix(Db);idDb.reshape(Indices(Db,1,Db,1));
    DoubleOperator UL(idDb,Uevol3.getOpData(0));
    DoubleOperator UR(Uevol3.getOpData(2),idDb);
    evolU.setOp(0,&UL,false);
    evolU.setOp(Lm-1,&UR,false);
    for(int l=1;l<Lm-1;l++) evolU.setOp(l,&Uevol3.getOp(1),false);
    cout<<"Prepared evol U with length "<<Lm<<endl;
    for(int k=0;k<m;k++){
      // Apply the step of usual evolution
      MPS aux(mps);
      contractor.optimize(evolU,aux,mps,Dcut);
      cnt++;
      if(cnt%rate==0){ // if needed, compute and save observables
	mps.gaugeCond('R',1);mps.gaugeCond('L',1);
	saveObservables(mps,cnt,delta,outfname); 
      }
    }
    cout<<"After evolving for "<<m<<" steps, about to trace out "<<Nx<<" sites per boundary "<<endl;
    // Now contract the boundaries and prepare the new effective system
    // the contraction of the boundary is equivalent to just tracing out those sites,
    // but imposing the right gauge conditions, that yields just the identity
    MPS aux(mps);
    for(int l=0;l<Nx;l++) aux.gaugeCond(l,'R',true);
    for(int l=Lm-1;l>Lm-1-Nx;l--) aux.gaugeCond(l,'L',true);
    mps=MPS(Lm-2*Nx,Dcut,d*d);
    //cout<<"Prepared place for new MPS of length "<<mps.getLength()<<endl;
    for(int l=1;l<Lm-2*Nx-1;l++) // place the old tensors
      mps.replaceSite(l,aux.getA(l+Nx),false);
    // the edges are a bit different
    mwArray AL=aux.getA(Nx).getA();
    mwArray AR=aux.getA(Lm-1-Nx).getA();
    // New boundary dimension is the left bond dim of site Nx (right of Lc-1-Nx)
    Indices dimsL=AL.getDimensions();
    Indices dimsR=AR.getDimensions();
    Db=dimsL[1]; // check with existing? with right?
    //cout<<"For boundaries, AL:"<<AL.getDimensions()<<", AR:"<<AR.getDimensions()<<"; new Db="<<Db<<endl;
    AL.permute(Indices(2,1,3));AL.reshape(Indices(Db*d*d,1,dimsL[2]));
    AR.permute(Indices(1,3,2));AR.reshape(Indices(d*d*Db,dimsR[1],1));
    mps.replaceSite(0,AL,false);
    mps.replaceSite(Lm-2*Nx-1,AR,false);
    mps.gaugeCond('R',1);mps.gaugeCond('L',1);
    cntX++;  Lm=Lm-2*Nx;
  }
  
}



int getInitialSize(int M,int Ls,int Nx){
  return getInitialSizeQuadratic(M,Ls,Nx);
}

int getNrTSteps(int m){
  return getNrTStepsQuadratic(m);
}

void saveObservables(const MPS& mps,int cnt,double delta,const string& outfname){
  static bool first=true;
  // Operators needed (only local, but doubled)
  static mwArray sigX;
  static mwArray sigY;
  static mwArray sigZ;

  ofstream* out;

  int L=mps.getLength();
  Contractor& contractor=Contractor::theContractor();  
  if(first){
    // Only the first time, if needed, write the header of the file
    if(cnt==0||!file_exists(outfname.data())){
      out=new ofstream(outfname.data());
      *out<<"% cnt\t t\t <sigX(L/2-1)>\t <sigY(L/2-1)>\t <sigZ(L/2-1)>"
	  <<"\t <sigX(L/2)>\t <sigY(L/2)>\t <sigZ(L/2)>"
	  <<"\t <sigZ(L/2-1) sigZ(L/2)>"<<endl;
      out->close();delete out;
    }
    // And initialize the actual (doubled) operators
    complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    complex_t dataY[]={ZERO_c,I_c,-1*I_c,ZERO_c};
    complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c}; //Z
    mwArray _sig0=identityMatrix(d);
    mwArray _sigX=mwArray(Indices(d,d),dataX);
    mwArray _sigY(Indices(d,d),dataY);
    mwArray _sigZ(Indices(d,d),dataZ);
    prepareOperatorProduct(sigX,_sigX,_sig0);
    prepareOperatorProduct(sigY,_sigY,_sig0);
    prepareOperatorProduct(sigZ,_sigZ,_sig0);
    first =false;
  }
  // Compute and write expectation values
  // Open the file to write them
  out=new ofstream(outfname.data(),ios::app);
  *out<<setprecision(15);
  *out<<cnt<<"\t"<<cnt*delta<<"\t";
  // the single site on site L/2-1
  MPS aux(mps);
  aux.applyLocalOperator(L/2-1,sigX,false);
  complex_t valX=contractor.contract(aux,mps);
  cout<<"cnt: "<<cnt<<" <s_x>="<<valX<<endl;
  aux.applyLocalOperator(L/2-1,sigY*sigX,false);
  complex_t valY=contractor.contract(aux,mps);
  aux.applyLocalOperator(L/2-1,sigZ*sigY,false);
  complex_t valZ=contractor.contract(aux,mps);
  aux.applyLocalOperator(L/2,sigZ,false);
  *out<<real(valX)<<"\t"<<real(valY)<<"\t"<<real(valZ)<<"\t";
  complex_t valZZ=contractor.contract(aux,mps);
  aux.applyLocalOperator(L/2-1,sigZ,false);
  valZ=contractor.contract(aux,mps);
  aux.applyLocalOperator(L/2,sigX*sigZ,false);
  valX=contractor.contract(aux,mps);
  aux.applyLocalOperator(L/2,sigY*sigX,false);
  valY=contractor.contract(aux,mps);
  *out<<real(valX)<<"\t"<<real(valY)<<"\t"<<real(valZ)<<"\t";
  *out<<real(valZZ)<<"\t";
  *out<<endl;  
  out->close();delete out;
}

int getInitialSizeQuadratic(int M,int Ls,int Nx){
  return Ls+2*Nx*ceil(sqrt(M));
}

int getNrTStepsQuadratic(int m){
  return 2*m+1;
}

void constructInitialProductState(int d,int D,int initSt,mwArray& W){
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sig0=identityMatrix(d);
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray sigX(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  
  switch(initSt){
  case 1: // X+
    {
      W=.5*(sig0+sigX);
      break;
    }
  case 2: // Y+
    {
      W=.5*(sig0+sigY);
      break;
    }
  case 3: // Z+
    {
      W=.5*(sig0+sigZ);
      break;
    }
  default:
    cout<<"Error: unknown type of intSt="<<initSt<<endl;
    exit(1);
  }
  W.reshape(Indices(d*d,1,1));
  W.resize(Indices(d*d,D,D));
}

void prepareOperatorProduct(mwArray& result,const mwArray& opA,const mwArray& opB){
  result=opA; 
  result.reshape(Indices(d*d,1));
  result.multiplyRight(reshape(opB,Indices(1,d*d)));
  result.reshape(Indices(d,d,d,d)); // order here: ud ud
  result.permute(Indices(1,3,2,4));
  result.reshape(Indices(d*d,d*d));
}



