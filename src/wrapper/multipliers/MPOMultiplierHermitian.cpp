
#include "MPOMultiplierHermitian.h"
#include <cmath>

using namespace std;
using namespace shrt;

MPOMultiplierHermitian::MPOMultiplierHermitian(int length):MPOMultiplier(length),mpoH(length){
  init=false;
  totalDim=0;
}

MPOMultiplierHermitian::MPOMultiplierHermitian(int length,const Operator** oplist):MPOMultiplier(length,oplist),mpoH(length){
  // init the HErmitian conjugate
  for(int k=0;k<length;k++){
    mpoH.setRotatedOp(k,oplist[k]->getFullData(),Indices(3,2,1,4),true);
  }
  initDims();
}

MPOMultiplierHermitian::MPOMultiplierHermitian(const MPO& mpo):MPOMultiplier(mpo),mpoH(mpo.getLength()){
  // init the HErmitian conjugate
  for(int k=0;k<length;k++){
    mpoH.setRotatedOp(k,mpo.getOpData(k),Indices(3,2,1,4),true);
  }  
  initDims();
}


void MPOMultiplierHermitian::setOp(int pos,const Operator* op,bool myop){
  MPOMultiplier::setOp(pos,op,myop);
  mpoH.setRotatedOp(pos,op->getFullData(),Indices(3,2,1,4),true);
  // if(init){
  //   int replDim=dimsIn[pos];
  //   totalDim=totalDim*op->getdorig()/dimsIn[pos];
  //   dimsIn[pos]=op->getdorig();
  // }
  // else{
  //   opset[pos]=true;
  //   if(allSet()) initDims();
  // }
}

void MPOMultiplierHermitian::getFullTensor(mwArray& result){
  MPOMultiplier::getFullTensor(result);
  result=.5*(result+Hconjugate(result));
}

void MPOMultiplierHermitian::setRotatedOp(int pos,const mwArray& op,
				 const shrt::Indices& neworder,
				 bool conjugate){
  MPOMultiplier::setRotatedOp(pos,op,neworder,conjugate);
  mpoH.setRotatedOp(pos,this->getOpData(pos),Indices(3,2,1,4),true);
  // if(init){
  //   int replDim=dimsIn[pos];
  //   totalDim=totalDim*getOp(pos).getdorig()/dimsIn[pos];
  //   dimsIn[pos]=getOp(pos).getdorig();
  // }
  // else{
  //   opset[pos]=true;
  //   if(allSet()) initDims();
  // }
}

// int MPOMultiplierHermitian::getSize() const{
//   if(init) return totalDim;
//   cout<<"Error: MPOMultiplierHermitian::getsize() can only be called when "
//       <<"all operators have been set"<<endl;
//   exit(1);
// }


void MPOMultiplierHermitian::product(const mwArray& input,mwArray& result){
  MPOMultiplier::product(input,result);
  mwArray resultH;
  mpoH.product(input,resultH);
  result=.5*(result+resultH);
}

#ifdef TESTINGTIMES
#include "sys/time.h"
#endif

 MPOMultiplierHermitian::~MPOMultiplierHermitian(){
#ifdef TESTINGTIMES
  cout<<"Destroying MPOMultiplierHermitian, called "<<counter<<" times, with average time "
      <<timer*1./(counter)<<endl;
#endif
}

