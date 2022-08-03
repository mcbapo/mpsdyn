#ifndef NEWIMPLEM
#include "libdiminf.h"
#include "diminf.h"
#else
#include "mwArray.h"
#include "MPO.h"
#endif
#include "Operator.h"
#include "Contractor.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>

#include "time.h"
clock_t start,finish;

#ifndef NEWIMPLEM
typedef OperatorRow MPO;
#endif


#ifdef MACOSX
#ifdef MCLAPTOP
#define LOCALDIR "/Users/banuls/SVN/peps/projects/mpsdyn/src/bin/"
#else
#define LOCALDIR "/Users/banulsm/Desktop/SVN/projects/mpsdyn/src/bin/"
#endif
#else
#define LOCALDIR "/home/banulsm/ALaPorra/bin/"
#endif


/** 
    Test program Nr. 7: Check the correct behaviour of Contractor class
    To compile: make test7
*/

using namespace std;
#ifdef NEWIMPLEM
using namespace shrt;
#endif

int main(){

  char filename[120];

#ifndef NEWIMPLEM
  int error;
  initializelib(&error);
  try{
    initGlobals(mwArray(CHECKDIMS),mwArray(SEEDVAL));
#endif

    int d=2, D=40, N=40, Dmax=(int)pow(2.,N/2.);Dmax=6;
    vector<int> dims(N);
    for(int k=0;k<N;k++){
      //dims[k]=k==0?1:d;
      dims[k]=d;
    }
  MPS psi1(N,D,dims);
  MPS psi2(N,D,dims);
  Contractor& theContractor=Contractor::theContractor();
#ifndef NEWIMPLEM
  mwArray result;
#else
  complex_t result;
#endif
  ProductState st1=p_yplus;
  psi1.setProductState(p_yplus);
  for(int k=3;k<4;k++){
    psi2.setProductState(ProductState(k));
    start=clock();
    result=theContractor.contract(psi1,psi2);
    cout<<"Contracting MPS@"<<&psi1<<", norm "<<
      theContractor.contract(psi1,psi1)<<endl;
    cout<<"Contracting MPS@"<<&psi2<<", norm "<<
      theContractor.contract(psi2,psi2)<<endl;
    cout<<"Result of contraction "<<st1<<"X"<<k<<"="<<result<<endl;
    cout<<"Result of contractionR "<<st1<<"X"<<k<<"="
	<<theContractor.contract(psi1,psi2,'R');
    finish=clock();
    cout<<", time "<<finish-start<<endl;
  }

  // Try now with random MPS/MPO
#ifdef NEWIMPLEM
  strcpy(filename,LOCALDIR);
  psi1.setRandomState();
  //psi1.exportMPStext(strcat(filename,"Psi1.mat"));
  psi1.importMPStext(strcat(filename,"Psi1.mat"));
  psi2.setRandomState();
  strcpy(filename,LOCALDIR);
  //psi2.exportMPStext(strcat(filename,"Psi2.mat"));
  psi2.importMPStext(strcat(filename,"Psi2.mat"));

#else
  strcpy(filename,LOCALDIR);
  psi1.importMPStext(strcat(filename,"Psi1.mat"));
  strcpy(filename,LOCALDIR);
  psi2.importMPStext(strcat(filename,"Psi2.mat"));
#endif

  cout<<"Result of contraction(L) ="
      <<theContractor.contract(psi1,psi2,'L')<<endl;
  cout<<"Result of contraction(R) ="
      <<theContractor.contract(psi1,psi2,'R')<<endl;

  //Now the MPO
  int Dop=3;
#ifdef NEWIMPLEM
  MPO ops(N);
  // Assign random Ops
  for(int k=0;k<N;k++){
    int Dlop=Dop;int Drop=Dop;
    if(k==0) Dlop=1;
    if(k==N-1) Drop=1;
    int du=dims[k];int dd=dims[k];
    mwArray aux(Indices(du,Dlop,dd,Drop));aux.fillRandom();
    Operator* op=new Operator(aux);
    ops.setOp(k,op);
  }
  strcpy(filename,LOCALDIR);
  //ops.exportMPOtext(strcat(filename,"MPO.mat"));
  ops.importMPOtext(strcat(filename,"MPO.mat"));
#else
  MPO ops(N,Dop,d);
  strcpy(filename,LOCALDIR);
  ops.importMPOtext(strcat(filename,"MPO.mat"));
#endif
  
  start=clock();
  cout<<"Result from contraction <Psi1|Ops|Psi2>="
      <<theContractor.contract(psi2,ops,psi1)<<endl;
  cout<<"Result from contractionR <Psi2|Ops|Psi1>="
      <<theContractor.contractR(psi2,ops,psi1)<<endl;
  finish=clock();
  cout<<" ===== Time for both "<<finish-start<<endl;

  start=clock();
  cout<<"Result from contraction <Psi2|Ops+ Ops|Psi2>="
      <<theContractor.contract2(ops,psi2)<<endl;
  cout<<"Result from contraction <Psi1|Ops Ops+|Psi1>="
      <<theContractor.contract2(psi1,ops)<<endl;
  cout<<"Result from contraction <Psi1|Ops+ Ops|Psi2>="
      <<theContractor.contract2(psi1,ops,psi2)<<endl;
  finish=clock();
  cout<<" ===== Time for three last "<<finish-start<<endl;

  psi1.gaugeCond('R',1);
  double entropyPsi1;
  int posi=0;//posi=N-1;
  //  for(int posi=1;posi<=N/2;posi++){
#ifdef NEWIMPLEM
  entropyPsi1=theContractor.getEntropy(psi1,posi);
#else
  entropyPsi1=psi1.getEntropy(posi);
#endif
  cout<<"Entropy of normalized Psi1(@"<<posi<<")="<<entropyPsi1<<endl;
  // }
  cout<<"And after imposing gauge to Psi1, <Psi1|Psi1>="
      <<theContractor.contract(psi1,psi1)<<endl;
  cout<<" and <Psi1|Psi2>="
      <<theContractor.contract(psi1,psi2,'L')<<endl;

#ifdef NEWIMPLEM
  // Test the Schmidt values in the middle
  std::vector<complex_t> lambdas;
  theContractor.getSchmidtValues(psi1,lambdas);
  cout<<"Schmidt values in the middle="<<lambdas<<endl;
  // Test the entropy of the middle subchain
  int posL=15;int posR=25;
  psi1.exportMPStext(strcat(filename,"Psi1.mat"));
  double entropyMPsi1=theContractor.getEntropy(psi1,posL,posR);
  cout<<"Entropy of the middle 10 sites ("<<posL<<"-"<<posR<<") "
      <<entropyMPsi1<<endl;
  //int n=50;
  //entropyMPsi1=theContractor.approximateEntropy(psi1,posL,posR,n);
  //cout<<"Approximated entropy of the middle "<<posR-posL+1<<" sites for n="<<n<<" : "
  //  <<entropyMPsi1<<endl;
  //  exit(1);
#endif

  //  exit(0);

  // Now try the optimization routines
  MPS final(psi2);
  start=clock();
  theContractor.optimize(ops,psi2,final);
  finish=clock();
  cout<<"After optimization Ops Psi2 ->> <result|result>="
      <<theContractor.contract(final,final)<<endl;
  cout<<"After optimization Ops Psi2 ->> <result|Ops|result>="
      <<theContractor.contract(final,ops,final)<<endl;
  cout<<"After optimization Ops Psi2 ->> <result|Ops|Psi2>="
      <<theContractor.contract(psi2,ops,final)<<endl;
  cout<<"And the state is "<<final<<endl;
  cout<<" ===== Optimize time="<<finish-start<<endl;

  start=clock();
  theContractor.optimizeL(ops,psi2,final);
  finish=clock();
  cout<<"After optimization Ops to L Psi2 ->> <result|result>="
      <<theContractor.contract(final,final)<<endl;
  cout<<"After optimization to L Ops Psi2 ->> <result|Ops|result>="
      <<theContractor.contract(final,ops,final)<<endl;
  cout<<"After optimization to L Ops Psi2 ->> <result|Ops|Psi2>="
      <<theContractor.contract(psi2,ops,final)<<endl;
  cout<<" ===== Optimize time="<<finish-start<<endl;

  // test dominant eigenvectors
  double lambda=0.;
  start=clock();
  theContractor.findRightEigenvector(ops,Dmax,&lambda,final);
  finish=clock();
  cout<<"After findRightEigenvector, <result|ops|result>="
     <<theContractor.contract(final,ops,final)<<endl;
  cout<<" Norm of result="<<theContractor.contract(final,final)<<endl;
  cout<<" Lambda="<<lambda<<endl;
  cout<<"State final="<<final<<endl;
  cout<<" ===== findRightEigenvector time="<<finish-start<<endl;


  MPS finalL(final);
  start=clock();
  theContractor.findLeftEigenvector(ops,Dmax,&lambda,finalL);
  finish=clock();
  cout<<"After findLeftEigenvector, <result|ops|result>="
     <<theContractor.contract(final,ops,finalL)<<endl;
  cout<<" Norm of result="<<theContractor.contract(finalL,finalL)<<endl;
  cout<<" LambdaL="<<lambda<<endl;
  cout<<" Overlap <L|R>="<<theContractor.contract(final,finalL)<<endl;
  cout<<"State finalL="<<finalL<<endl;
  cout<<" ===== findLeftEigenvector time="<<finish-start<<endl;

  // Now gauge condition between both!
  MPS::gaugeCond(final,finalL,'R');
  cout<<" After gaugeCond <L|R>="<<theContractor.contract(final,finalL)<<endl;
  cout<<" Norm of L="<<theContractor.contract(finalL,finalL)<<endl;
  cout<<" Norm of R="<<theContractor.contract(final,final)<<endl;

  // Test findGroundState (construct Hermitian Ops)
#ifdef NEWIMPLEM
  MPO opsH(N);
  // Assign random Ops
  for(int k=0;k<N;k++){
    int Dlop=Dop;int Drop=Dop;
    if(k==0) Dlop=1;
    if(k==N-1) Drop=1;
    int du=dims[k];
    mwArray aux(Indices(du,Dlop,du,Drop));aux.fillRandom();
    aux=aux+permute(conjugate(aux),Indices(3,2,1,4));
    Operator* op=new Operator(aux);
    opsH.setOp(k,op);
  }
  opsH.exportMPOtext(strcat(filename,"MPO_h.mat"));
#else
  MPO opsH(N,Dop,d);
  opsH.importMPOtext(strcat(filename,"MPO_h.mat"));
#endif
  cout<<"Read or created opsH "<<opsH<<endl;
  cout<<"<L|opsH|R>="<<theContractor.contract(final,opsH,finalL)<<endl;

  MPS finalGS(final);
  start=clock();
  theContractor.findGroundState(opsH,Dmax,&lambda,finalGS);
  finish=clock();
  cout<<"After findGroundState, <result|opsH|result>="
     <<theContractor.contract(finalGS,opsH,finalGS)<<endl;
  cout<<" Norm of result="<<theContractor.contract(finalGS,finalGS)<<endl;
  cout<<" Lambda="<<lambda<<endl;
  cout<<" ===== findǴroundState time="<<finish-start<<endl;

  // Test optimizeInv
  MPS finalI(final);MPS psi3(psi2);
  theContractor.optimize(opsH,psi2,psi3,Dmax);
  cout<<"opsH Psi2=Psi3, norm="<<theContractor.contract(psi3,psi3)
      <<" and overlap with psi2="<<theContractor.contract(psi2,psi3)<<endl;
  start=clock();
  theContractor.optimizeInv(opsH,psi3,finalI,Dmax);
  finish=clock();
  cout<<"After optimizeInv, <psi3|opsH|result>="
     <<theContractor.contract(finalI,opsH,psi3)<<endl;
  cout<<" Norm of result="<<theContractor.contract(finalI,finalI)<<endl;
  cout<<" Overlap with true result Psi2="
      <<theContractor.contract(psi2,finalI)<<endl;
  cout<<" ===== optimizeInv time="<<finish-start<<endl;

  // Test optimizeResolvent
  MPS psi4(psi3);double eta=0.1;
  start=clock();
  theContractor.optimizeResolvent(opsH,psi3,psi4,eta,Dmax);
  finish=clock();
   cout<<"After optimizeInv, <psi3|(opsH+i eta)|psi4>="
#ifdef NEWIMPLEM
       <<theContractor.contract(psi4,opsH,psi3)+
     eta*I_c*theContractor.contract(psi4,psi3)<<endl;
#else
   <<theContractor.contract(psi4,opsH,psi3)<<"+eta*i*"
   <<theContractor.contract(psi4,psi3)<<endl;
#endif
   cout<<" Norm of result="<<theContractor.contract(psi4,psi4)<<endl;
  cout<<" ===== optimizeResolvent time="<<finish-start<<endl;

   // Test approximateExponential
   MPS Phi0(psi1);
   MPS PhiDelta(psi1);double delta=0.02;
   theContractor.approximateExponential(opsH,Phi0,PhiDelta,delta,Dmax);
   cout<<"After approximateExponential, <Phi'|Phi'>="
       <<theContractor.contract(PhiDelta,PhiDelta)<<endl;
   cout<<"After approximateExponential, <Phi'|opsH|Phi'>="
       <<theContractor.contract(PhiDelta,opsH,PhiDelta)<<endl;
   cout<<"After approximateExponential, <Phi'|Phi0>="
       <<theContractor.contract(Phi0,PhiDelta)<<endl;

   // Test summing MPSs
   vector<const MPS*> list;vector<complex_t> listBeta;
   list.push_back(&psi1);listBeta.push_back((complex_t){2.,0.});
   list.push_back(&psi1);listBeta.push_back((complex_t){0.,2.});
   MPS nPsi(N,D,d);nPsi.setRandomState();
   theContractor.optimizeSum(list,listBeta,nPsi,D);
   cout<<"After optimizing the sum, <result|result>="
       <<theContractor.contract(nPsi,nPsi)<<endl;
   cout<<"After optimizing the sum, <psi1|result>="
       <<theContractor.contract(nPsi,psi1)<<endl;

#ifndef NEWIMPLEM
  }
  catch(mwException& e){
    cout<<"Exception caught: "<<e.what()<<"in testContractor"<<endl;
    exit(212);
  }
  closelib();
  return error;
#endif

}
