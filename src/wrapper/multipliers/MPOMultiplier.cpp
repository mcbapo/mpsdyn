
#include "MPOMultiplier.h"
#include "Indices.h"
#include <cmath>

int MAXDIM=pow(2,14);

using namespace std;
using namespace shrt;

MPOMultiplier::MPOMultiplier(int length):MPO(length),opset(length,false){
  init=false;
  totalDim=0;
}

void MPOMultiplier::initDims(){
  dimsIn=getOriginalDimensions();
  totalDim=1;
  for(int k=0;k<getLength();k++){
    totalDim*=dimsIn[k];
  }
  init=true;
}

bool MPOMultiplier::allSet() const{
  for(int k=0;k<getLength();k++)
    if(!opset[k]) return false;
  return true;
}

MPOMultiplier::MPOMultiplier(int length,const Operator** oplist):
  MPO(length,oplist),opset(length,true){
  initDims();
}

MPOMultiplier::MPOMultiplier(const MPO& mpo):MPO(1,mpo.Ops),opset(length,true){  
  initDims();
}


void MPOMultiplier::setOp(int pos,const Operator* op,bool myop){
  MPO::setOp(pos,op,myop);
  if(init){
    int replDim=dimsIn[pos];
    totalDim=totalDim*op->getdorig()/dimsIn[pos];
    dimsIn[pos]=op->getdorig();
  }
  else{
    opset[pos]=true;
    if(allSet()) initDims();
  }
}

void MPOMultiplier::getFullTensor(mwArray& result){
  if(totalDim>MAXDIM){
    cout<<"ERROR: MPOMultiplier would return a too big matrix (dim"<<totalDim<<"). Aborting."<<endl;
    cout<<"To allow for larger matrices, change the value of MAXDIM in MPOMultiplier.cpp"<<endl;
    exit(1);
  }
  expandOper(*this,result);
}

void MPOMultiplier::setRotatedOp(int pos,const mwArray& op,
				 const shrt::Indices& neworder,
				 bool conjugate){
  MPO::setRotatedOp(pos,op,neworder,conjugate);
  if(init){
    int replDim=dimsIn[pos];
    totalDim=totalDim*getOp(pos).getdorig()/dimsIn[pos];
    dimsIn[pos]=getOp(pos).getdorig();
  }
  else{
    opset[pos]=true;
    if(allSet()) initDims();
  }
}

int MPOMultiplier::getSize() const{
  if(init) return totalDim;
  cout<<"Error: MPOMultiplier::getsize() can only be called when "
      <<"all operators have been set"<<endl;
  exit(1);
}

#ifdef TESTINGTIMES
#include "sys/time.h"
#endif

void MPOMultiplier::product(const mwArray& input,mwArray& result){
  if(!init){
    cout<<"Error: MPOMultiplier("<<this<<")::product() can only be called when "
	<<"all operators have been set"<<endl;
    exit(1);
  }
  // Perform the contraction---TODO:CHECK!
#ifdef TESTINGTIMES
  struct timeval start,final;
  gettimeofday(&start,NULL);
#endif
  //cout<<"Call to MPOMultiplier::product("<<&input<<","<<input.getDimensions()<<"), result:"<<&result<<","<<result.getDimensions()<<endl;
  //result=input;
  mwArray aux0(input);
  int dimL=1; //"open" right bond left from previous multiplication
  int finalDim=1;
  for(int k=0;k<getLength();k++){
    // one by one, multiply one term, and reshape adequately
    Indices dims=getOp(k).getDimensions();
    int dd=dims[0];int du=dims[2];
    int Dl=dims[1];int Dr=dims[3];
    aux0.reshape(Indices(dimL*dd,-1));
    mwArray aux(getOpData(k));
    //cout<<"MPOMultiplier("<<this<<") prod (before k="<<k<<"): "<<result.getDimensions()<<", op:"<<aux.getDimensions()<<endl;
    aux.permute(Indices(1,4,2,3));
    aux.reshape(Indices(du*Dr,Dl*dd));
    aux0.multiplyLeft(aux);
    aux0.reshape(Indices(du,-1));
    aux0.transpose();
    dimL=Dr;
    finalDim*=du;
  }
  result=aux0;
  result.reshape(Indices(finalDim,1));
#ifdef TESTINGTIMES
  counter++;
  gettimeofday(&final,NULL);
  timer+=(final.tv_sec-start.tv_sec)*1E3+(final.tv_usec-start.tv_usec)*1E-3;
#endif
}

MPOMultiplier::~MPOMultiplier(){
#ifdef TESTINGTIMES
  cout<<"Destroying MPOMultiplier, called "<<counter<<" times, with average time "
      <<timer*1./(counter)<<endl;
#endif
}

