/**
   \file LMGHamiltonian.cpp
   Hamiltonian for the LMG model 
   
   \author Mari-Carmen Banuls
   \date 21/12/2016

*/

#include "LMGHamiltonian.h"
#include "Indices.h"
//#include "complex.h"

using namespace std;
using namespace shrt;

LMGHamiltonian::LMGHamiltonian(int N_,double h0_):
  hamil(N_),N(N_),d(2),h0(h0_),Z(){
  if(d!=2){
    cout<<"Error: LMGHamiltonian only defined for d=2"<<endl;
    exit(212);
  }
  initZ();
  setHMPO(hamil,J0,h0);
  J=h=0;
}

LMGHamiltonian::~LMGHamiltonian(){};


void LMGHamiltonian::initZ(){
  // basic spin operators appearing
  mwArray sig0=identityMatrix(d);
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sig3(Indices(d,d),dataz);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig1(Indices(d,d),datax);
  Z=mwArray(Indices(3,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Z.setElement(sig3.getElement(Indices(i1,i2)),Indices(2,i1,i2));
    }
  Z.reshape(Indices(3,d*d));
}

void LMGHamiltonian::setHMPO(MPO& mpo,double h){
  int D=3;
  for(int k=0;k<N;k++){
    // first and last are different
    if(k>1&&k<N-1){
      mpo.setOp(k,&mpo.getOp(k-1));
    }
    else{
      //cout<<"Creating new Operator for site "<<k<<" of LMGH MPO"<<endl;
      int Dl=D; int Dr=D;
      if(k==0) Dl=1;
      if(k==N-1) Dr=1;
      mwArray C(Indices(Dl,Dr,3));C.fillWithZero();
      complex_t effH=-.5*h*ONE_c; //h/2
      double factor=1./sqrt(2*N)*ONE_c;
      //cout<<"Coefficients: effJ="<<effJ<<", effG="<<effG<<", effH="<<effH
      //  <<endl;
      // Identity
      if(k!=N-1)
	C.setElement(ONE_c,Indices(0,0,0)); // just pass 1 (before anything)
      else
	C.setElement(-ONE_c,Indices(0,0,0)); // global constant -1
      if(k!=0)
	C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // just pass 1 (after everything)
      // sigmaz localterms (in all of them)
      C.setElement(effH,Indices(0,Dr-1,2));
      // left sigx term (all but last)
      if(k!=N-1){
	C.setElement(factor,Indices(0,1,1));
	if(k!=0) // and pass it on, except for first site
	  C.setElement(factor,Indices(1,1,0));
      }
      // right sigz term (all but first)
      if(k!=0)
	C.setElement(-factor,Indices(1,Dr-1,1));
      
      // Now contract and give the operator
      // Now contract MPS and operators and set proper index ordering
      C.reshape(Indices(Dl*Dr,3));
      //cout<<"C for site "<<k<<"="<<C<<endl;
      mwArray res=C*Z;
      //cout<<"Computed res "<<res<<endl;
      res.reshape(Indices(Dl,Dr,d,d));
      res.permute(Indices(3,1,4,2));
      //cout<<"  with size: "<<res.getDimensions()<<endl;
      // Create and store an Operator
      mpo.setOp(k,new Operator(res),true);
    }
  }
}
