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


int getInitialSize(int M,int Ls,double alpha);
int getNrTSteps(int m,double alpha);
int getInitialSizeQuadratic(int M,int Ls,double alpha);
int getNrTStepsQuadratic(int m,double alpha);
void saveObservables(const MPS& mps,int cnt,double delta,const string& outfname);

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
  double alpha=props.getDoubleProperty("alpha");
  if(alpha<0) alpha=1.;
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
  int Nx=1;
  int DcutS=props.getIntProperty("D"); // the truncation INSIDE the system
  if(DcutS<=0) DcutS=Dcut; // default like the other
  bool hardB=(props.getIntProperty("hardB")>0); // if set, no reduction of the LC is done, as if hard boundary
  // string tmpdir=props.getProperty("tmpdir"); // directory to save
  // 					     // temporary tensors
  // 					     // (big!)
  // if(savingTmp&&tmpdir.empty()){
  //   cout<<"ERROR! No tmpdir specified to save temporary data"<<endl;
  //   exit(1);
  // }

  // 1) Prepare initial size. I assume it is a product, so the boundary D at the beginning is 1
  int L0=getInitialSize(M,Ls,alpha);

  MPS mps(L0,1,d);
  switch(initSt){
  case 1:
    mps.setProductState(p_xplus);break;
  case 2:
    mps.setProductState(p_yplus);break;
  case 3:
    mps.setProductState(p_zero);break;
  default:
    cout<<"Error: unknown type of initSt="<<initSt<<endl;
    exit(1);
  }

  // My Hamiltonian
  IsingHamiltonian hamH(3,d,J_,g_,h_);
  MPO Uevol3(3); // to hold the different operators
  hamH.getUMPO(Uevol3,-delta*I_c,0); // the single layer normal evolution

  Contractor& contractor=Contractor::theContractor();
  
  int cnt=0; // total nr of steps (real time)
  int cntX=0; // nr of times I've folded in boundary sites (evol in x)
  int Lm=L0; // current size of the MPS (incl. boundary)
  int Db=1; // dimension of the virtual leg at the boundary, representing the environment

  saveObservables(mps,cnt,delta,outfname); 
  
  while(cnt<M){
    // How many normal time steps do I need to apply for this stage?
    int m=getNrTSteps(cntX,alpha);
    //    m=1; // trick to do linear
    if(hardB){
      cout<<"WARNING: using hard boundaries: no cone reduction"<<endl;
      m=M; // trick to keep hard boundaries!! All steps at once
    }
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
    if(!hardB){
      cout<<"After evolving for "<<m<<" steps, about to trace out "<<Nx<<" site(s) per boundary "<<endl;
      // Now contract the boundaries and prepare the new effective system
      // the contraction of the boundary is equivalent to just tracing out those sites,
      // but imposing the right gauge conditions, that yields just the identity
      MPS aux(mps);
      for(int l=0;l<Nx;l++) aux.gaugeCond(l,'R',true);
      for(int l=Lm-1;l>Lm-1-Nx;l--) aux.gaugeCond(l,'L',true);
      mps=MPS(Lm-2*Nx,Dcut,d);
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
      AL.permute(Indices(2,1,3));AL.reshape(Indices(Db*d,1,dimsL[2]));
      AR.permute(Indices(1,3,2));AR.reshape(Indices(d*Db,dimsR[1],1));
      mps.replaceSite(0,AL,false);
      mps.replaceSite(Lm-2*Nx-1,AR,false);
      mps.gaugeCond('R',1);mps.gaugeCond('L',1);
      cntX++;  Lm=Lm-2*Nx;
    }
  }
  
}



int getInitialSize(int M,int Ls,double alpha){
  return getInitialSizeQuadratic(M,Ls,alpha);
}

int getNrTSteps(int m,double alpha){
  return getNrTStepsQuadratic(m,alpha);
}

void saveObservables(const MPS& mps,int cnt,double delta,const string& outfname){
  static bool first=true;
  // Operators needed (only local)
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1*I_c,ZERO_c};
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c}; //Z
  mwArray sig0=identityMatrix(d);
  mwArray sigX=mwArray(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  mwArray sigZ(Indices(d,d),dataZ);

  ofstream* out;

  int L=mps.getLength();
  Contractor& contractor=Contractor::theContractor();  
  if(first){
    // Only the first time, if needed, write the header of the file
    if(cnt==0||!file_exists(outfname.data())){
      out=new ofstream(outfname.data());
      if(!out->is_open()){
	cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
	exit(1);
      }
      *out<<"% cnt\t t\t <sigX(L/2-1)>\t <sigY(L/2-1)>\t <sigZ(L/2-1)>"
	  <<"\t <sigX(L/2)>\t <sigY(L/2)>\t <sigZ(L/2)>"
	  <<"\t <sigZ(L/2-1) sigZ(L/2)>"<<endl;
      out->close();delete out;
    }    
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

int getInitialSizeQuadratic(int M,int Ls,double alpha){
  return Ls+2*ceil(sqrt(M)/alpha);
}

int getNrTStepsQuadratic(int m,double alpha){
  return alpha*alpha*(2*m+1);
}
