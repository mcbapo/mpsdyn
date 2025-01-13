#include "mwArray.h"
#include "misc.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "Properties.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>

using namespace std;
using namespace shrt;

/**
Compute the zero norm if the input MPS is a (classical) binary
string. Otherwise, the corresponding expectation value of number of
non-zero "components", the latter understood as having at least a
non-vanishing element (i.e. a particle) in a group of k neighboring sites.
 */

/* Prepare a MPO that implements the projector onto any total number
   of particles, for a chain of length L, with an extra tensor (whose
   index indicates the total N) to be contracted in position posIdx
   (between posIdx-1 and posIdx) */

void prepareCounterMPS(MPS& mps,int M);
void prepareNmodProjector(MPO& mpo,int N);
void prepareCompositeMPO(MPO& mpo,int L,int k);
void setRandomComputational(MPS& mps,stringstream& s,int k);

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  string directory="";
  // if(argc>2){
  //   const char* indir=argv[++cntr];
  //   directory=string(indir)+"/";
  // }
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  //const string outfile = props.getProperty("output");
  int L=props.getIntProperty("L"); // number of "components"
  int k=props.getIntProperty("k"); // number of bits per component=> size will be k*L
  int nrTests=props.getIntProperty("Ntest"); // number of random test strings
  int d=2;

  MPO mpo(k*L);
  prepareCompositeMPO(mpo,L,k);

  // Now I need to run tests
  Contractor& contractor=Contractor::theContractor();

  for(int g=0;g<nrTests;g++){
    // choose a random prod state
    MPS test(k*L,1,d);
    stringstream s;
    setRandomComputational(test,s,k);    

    complex_t value=contractor.contract(test,mpo,test);
    cout<<"State s("<<s.str()<<") value:"<<real(value)<<endl;

  }
  
  
}

//void prepareCounterMPS(MPS& mps,int M){
void prepareNmodProjector(MPO& mpo,int N){
  mpo.initLength(N);
  // left chain
  int D=2;int d=2;
  for(int k=0;k<N;k++){
    int Dl=(k==0)?1:D;
    mwArray opN(Indices(d,Dl,d,D));
    opN.fillWithZero();
    for(int p=0;p<Dl;p++){
      opN.setElement(ONE_c,Indices(0,p,0,p)); // no particle: stays the same
      opN.setElement(ONE_c,Indices(1,p,1,1)); // one particle: always returns 1
    }
    mpo.setOp(k,new Operator(opN),true);
  }
  // At the end, it is not really an MPO, as the last site has a hanging leg
  cout<<"Prepared Nmod projector for "<<N<<" sites "<<mpo<<endl;
}

void prepareCounterMPS(MPS& mps,int M){
  int d=2;int D=2;
  mps=MPS(M,D,d);
  mwArray basisChg(Indices(d,2));
  for(int p=0;p<d;p++){
    basisChg.setElement(ONE_c,Indices(p,0));
    int signP=p==0?1:-1;
    basisChg.setElement((1-signP)*.5*ONE_c,Indices(p,1));
  }
  for(int k=0;k<M;k++){
    int Dl=(k==0)?1:D;
    int Dr=(k==M-1)?1:D;
    mwArray C(Indices(d,Dl,Dr));
    if(k<M-1)
      C.setElement(ONE_c,Indices(0,0,0));
    if(k>0)
      C.setElement(ONE_c,Indices(0,Dl-1,Dr-1));
    C.setElement(ONE_c,Indices(1,0,Dr-1));
    // chg of basis
    C.reshape(Indices(2,Dl*Dr));
    C.multiplyLeft(basisChg);
    C.reshape(Indices(2,Dl,Dr));
    mps.replaceSite(k,C,false);
  }
  cout<<"Prepared counter MPS for "<<M<<" blocks: "<<mps<<endl;
}

void prepareCompositeMPO(MPO& mpo,int L,int k){
  MPO projN(k);
  MPS mpsC(L,1,1);
  int d=2;
  
  prepareNmodProjector(projN,k);
  prepareCounterMPS(mpsC,L);

  mpo.initLength(L*k);
  mwArray idOpMid_=identityMatrix(d*2); idOpMid_.reshape(Indices(d,2,d,2));
  Operator idOpMid(idOpMid_);
  mwArray idOp_=identityMatrix(d);idOp_.reshape(Indices(d,1,d,1));
  Operator idOp(idOp_);

  MPO aux1(k*L); // prod of "projectors", with last sites increasing phys dimension
  for(int l=0;l<L;l++){
    for(int s=0;s<k-1;s++){
      aux1.setOp(l*k+s,&projN.getOp(s),false);
    }
    mwArray data=projN.getOp(k-1).getFullData(); // dx2xdx2
    data.permute(Indices(1,4,2,3));
    data.reshape(Indices(d*2,2,d,1));
    aux1.setOp(l*k+k-1,new Operator(data),true);
  }
  cout<<"Set aux1 "<<aux1<<endl;
  // A second layer has just identities and the MPS components contracting the additional indices
  MPO aux2(k*L);
  for(int l=0;l<L;l++){
    for(int s=0;s<k-1;s++){
      if(l==0)
	aux2.setOp(l*k+s,&idOp,false);      
      else
	aux2.setOp(l*k+s,&idOpMid,false);      
    }
    mwArray auxA=mpsC.getA(l).getA(); // 2xDlxDr
    Indices dimsA=auxA.getDimensions();
    auxA.permute(Indices(2,1,3)); // dims1,dims0,dims2
    auxA.reshape(Indices(1,dimsA[1],dimsA[0],dimsA[2]));
    aux2.setOp(l*k+k-1,new DoubleOperator(idOp_,auxA),true);
  }
  cout<<"Set aux2 "<<aux2<<endl;
  const MPO* ptrs[]={&aux2,&aux1};
  MPO::join(2,ptrs,mpo);
  cout<<"Set mpo "<<mpo<<endl;
}

void setRandomComputational(MPS& mps,stringstream& s,int M){
  int L=mps.getLength();
  int d=2;
  mps=MPS(L,1,d);
  mwArray state0(Indices(d,1,1));state0.fillWithZero();state0.setElement(ONE_c,Indices(0,0,0));
  mwArray state1(Indices(d,1,1));state1.fillWithZero();state1.setElement(ONE_c,Indices(1,0,0));
  for(int k=0;k<L;k++){
    // random choice
    int state=((double)random()/RAND_MAX)<.25; //0 or 1
    if(state<.5){
      mps.replaceSite(k,state0,false);s<<"0 ";
    }
    else{
      mps.replaceSite(k,state1,false);s<<"1 ";
    }
    if((k+1)%M==0&&k<L-1) s<<":";
  }

}
