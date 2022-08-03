#include "BosonMPO.h"

using namespace std;
using namespace shrt;

void BosonMPO::getNumberMPO(int L,int Nmax,MPO& result){
  result.initLength(L);
  int d=Nmax+1;
  mwArray Z=mwArray(Indices(2,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++){
    Z.setElement(ONE_c,Indices(0,i1,i1)); // identity
    Z.setElement(i1*ONE_c,Indices(1,i1,i1)); // N
  }
  Z.reshape(Indices(2,d*d));

  // the coefficients
  int Dop=2;
  // Center site
  mwArray C(Indices(Dop,Dop,2));C.fillWithZero();
  C.setElement(ONE_c,Indices(0,0,0));
  C.setElement(ONE_c,Indices(0,1,1));
  C.setElement(ONE_c,Indices(1,1,0));
  C.reshape(Indices(Dop*Dop,2));

  // Edges
  mwArray CL(Indices(1,Dop,2));CL.fillWithZero();
  CL.setElement(ONE_c,Indices(0,0,0));
  CL.setElement(ONE_c,Indices(0,1,1));
  CL.reshape(Indices(Dop,2));
  mwArray CR(Indices(Dop,1,2));CR.fillWithZero();
  CR.setElement(ONE_c,Indices(0,0,1));
  CR.setElement(ONE_c,Indices(1,0,0));
  CR.reshape(Indices(Dop,2));

  mwArray Op=C*Z;
  Op.reshape(Indices(Dop,Dop,d,d));
  Op.permute(Indices(3,1,4,2));
  mwArray OpL=CL*Z;
  OpL.reshape(Indices(1,Dop,d,d));
  OpL.permute(Indices(3,1,4,2));
  mwArray OpR=CR*Z;
  OpR.reshape(Indices(Dop,1,d,d));
  OpR.permute(Indices(3,1,4,2));

  // set operators in place
  result.setOp(0,new Operator(OpL),true);
  result.setOp(L-1,new Operator(OpR),true);
  if(L>1){
    result.setOp(1,new Operator(Op),true);
    const Operator* basicOp=&result.getOp(1);
    for(int k=2;k<L-1;k++)
      result.setOp(k,basicOp);
  }
}

void BosonMPO::getParityMPO(int L,int Nmax,MPO& result,int i,int j,
			    double avrgN){
  result.initLength(L);
  int d=Nmax+1;
  mwArray id=identityMatrix(d);
  id.reshape(Indices(d,1,d,1));
  mwArray expN(Indices(d,d));expN.fillWithZero();
  for(int k=0;k<d;k++){
    expN.setElement(exp((k-avrgN)*M_PIl*I_c),Indices(k,k)); // exp(i pi N)
  }
  expN.reshape(Indices(d,1,d,1));
  // Now set the operators
  bool storedId=false;int ptrId=-1;
  bool storedExp=false;int ptrExp=-1;
  for(int k=0;k<L;k++){
    if(k<i||k>j){
      if(!storedId){
	result.setOp(k,new Operator(id),true);
	ptrId=k;storedId=true;
      }
      else{
	result.setOp(k,&result.getOp(ptrId));
      }
    }
    else{
      if(!storedExp){
	result.setOp(k,new Operator(expN),true);
	ptrExp=k;storedExp=true;
      }
      else{
	result.setOp(k,&result.getOp(ptrExp));
      }
    }
  }
}

void BosonMPO::getSingleSiteBosonNumber(int Nmax,mwArray& result){
  int d=Nmax+1;
  result=mwArray(Indices(d,d));result.fillWithZero();
  for(int i1=0;i1<d;i1++){
    result.setElement(i1*ONE_c,Indices(i1,i1)); // N
  }
}
