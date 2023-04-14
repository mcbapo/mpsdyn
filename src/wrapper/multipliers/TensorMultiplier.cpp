
#include "TensorMultiplier.h"
#include "Indices.h"

using namespace std;
using namespace shrt;

TensorMultiplier::TensorMultiplier(const mwArray& AL_,const mwArray& Aop_,
				   const mwArray& AR_):
  AL(),Aop(),AR(),betab(AL_.getDimension(0)),alpha(AL_.getDimension(1)),
  betak(AL_.getDimension(2)),du(Aop_.getDimension(0)),al1(Aop_.getDimension(1)),
  dd(Aop_.getDimension(2)),al2(Aop_.getDimension(3)),
  betabr(AR_.getDimension(0)),alphar(AR_.getDimension(1)),
  betakr(AR_.getDimension(2)){
#ifdef TESTFUN
  cout<<"Creating a basic TensorMultiplier with AL="<<AL_.getDimensions()<<", Aop="<<Aop_.getDimensions()
      <<", AR="<<AR_.getDimensions()<<endl;
#endif
  AL=reshape(AL_,Indices(betab*alpha,betak));
  Aop=Aop_;Aop.permute(Indices(1,4,2,3));Aop.reshape(Indices(du*al2,al1*dd));
  AR=reshape(AR_,Indices(betabr,alphar*betakr));AR.permute(Indices(2,1));
#ifdef CHECKDIMS
  if((al1!=alpha)||(al2!=alphar)){
    cout<<"Error: wrong dimensions constructing TensorMultiplier"
	<<" AL="<<AL_.getDimensions()<<", oper="<<Aop_.getDimensions()
	<<", AR="<<AR_.getDimensions()<<endl;
    exit(212);
  }
  // if((betab!=betak)||(betabr!=betakr)){
  //   cout<<"WARNING: different input and output dimensions constructing TensorMultiplier"
  // 	<<" AL="<<AL_.getDimensions()<<", oper="<<Aop_.getDimensions()
  // 	<<", AR="<<AR_.getDimensions()<<endl;
  //   //exit(212);
  // }
#endif
#ifdef TESTINGTIMES
  counter=0;timer=0;
#endif
}


// return the tensors with the original dimensions to be used by MPO
mwArray TensorMultiplier::getOriginalLeftTerm() const{
  return reshape(AL,Indices(betab,alpha,betak));	//betab*alpha,beta
}
mwArray TensorMultiplier::getOriginalMiddleTerm() const{
  mwArray oldAop(Aop);
  oldAop.reshape(Indices(du,al2,al1,dd));
  oldAop.permute(Indices(1,3,4,2));
  return oldAop;
}
mwArray TensorMultiplier::getOriginalRightTerm() const{
  return reshape(permute(AR,Indices(2,1)),Indices(betabr,alphar,betakr));
}


#ifdef TESTINGTIMES
#include "sys/time.h"
#endif

void TensorMultiplier::product(const mwArray& input,mwArray& result){
  // Perform the contraction---TODO:CHECK!
#ifdef TESTINGTIMES
  struct timeval start,final;
  gettimeofday(&start,NULL);
#endif
  mwArray aux=reshape(input,Indices(dd,betak,betakr));
  aux.permute(Indices(2,1,3));
  aux.reshape(Indices(betak,dd*betakr));
  // First contract tmpL with A
  //  aux=AL*aux;
  aux.multiplyLeft(AL);
  aux.reshape(Indices(betab,alpha*dd,betakr));
  aux.permute(Indices(2,1,3));
  aux.reshape(Indices(alpha*dd,betab*betakr));
  //mwArray tmpres=AL*reshape(permute(A,Indices(2,1,3)),Indices(betak,dd*betakr));
  // Second, contract oper
  //  aux=Aop*aux;
  aux.multiplyLeft(Aop);
  aux.reshape(Indices(du,al2,betab,betakr));
  aux.permute(Indices(1,3,2,4));
  //Finally, contract termR
  aux.reshape(Indices(du*betab,al2*betakr));
  //  result=aux*AR;
  aux.multiplyRight(AR);
  result=aux; //TODO: Make more efficient avoiding this copy
  result.reshape(Indices(-1,1));
#ifdef TESTINGTIMES
  counter++;
  gettimeofday(&final,NULL);
  timer+=(final.tv_sec-start.tv_sec)*1E3+(final.tv_usec-start.tv_usec)*1E-3;
#endif
}

TensorMultiplier::~TensorMultiplier(){
#ifdef TESTINGTIMES
  cout<<"Destroying TensorMultiplier, called "<<counter<<" times, with average time "
      <<timer*1./(counter)<<endl;
#endif
}

void TensorMultiplier::getFullTensor(mwArray& result){
  //  std::cout<<"In TensorMultiplier::getFullTensor"<<std::endl;
  result=AL;
  result.reshape(Indices(betab,alpha,betak));
  result.permute(Indices(1,3,2));result.reshape(Indices(betab*betak,alpha));
  mwArray tmp(Aop);
  tmp.reshape(Indices(du,al2,al1,dd));tmp.permute(Indices(3,1,4,2));
  tmp.reshape(Indices(al1,du*dd*al2));
  result.multiplyRight(tmp);
  result.reshape(Indices(betab*betak*du*dd,al2));
  tmp=reshape(AR,Indices(alphar,betakr*betabr));
  result.multiplyRight(tmp);
  result.reshape(Indices(betab,betak,du,dd,betakr,betabr));
  result.permute(Indices(3,1,6,4,2,5));
  result.reshape(Indices(du*betab*betabr,dd*betak*betakr));
}
