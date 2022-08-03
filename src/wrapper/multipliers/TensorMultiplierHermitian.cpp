
#include "TensorMultiplierHermitian.h"
#include "Indices.h"

using namespace std;
using namespace shrt;

TensorMultiplierHermitian::TensorMultiplierHermitian(const mwArray& AL_,const mwArray& Aop_,
						     const mwArray& AR_):
  TensorMultiplier(AL_,Aop_,AR_),ALh(),Aoph(),ARh(){
#ifdef TESTFUN
  cout<<"Creating a basic TensorMultiplierHermitian with AL="<<AL_.getDimensions()<<", Aop="<<Aop_.getDimensions()
      <<", AR="<<AR_.getDimensions()<<endl;
#endif
#ifdef CHECKDIMS
  if((betak!=betab)||(betakr!=betabr)||(dd!=du)){
    cout<<"Wrong dimensions (not symmetric) for Hermitian TensorMultiplier, as AL->"<<AL_.getDimensions()
	<<", Aop->"<<Aop_.getDimensions()<<", AR->"<<AR_.getDimensions()<<endl;
    exit(212);
  }
#endif
  ALh=AL_;ALh.permute(Indices(3,2,1),true); // conjugate
  ALh.reshape(Indices(betab*alpha,betak));
  Aoph=Aop_;Aoph.permute(Indices(3,2,1,4),true); // conjugate
  Aoph.permute(Indices(1,4,2,3));Aoph.reshape(Indices(du*al2,al1*dd));
  ARh=AR_;ARh.permute(Indices(3,2,1),true); // conjugate
  ARh.reshape(Indices(betabr,alphar*betakr));ARh.permute(Indices(2,1));
#ifdef TESTINGTIMES
  counter=0;timer=0;
#endif
}

#ifdef TESTINGTIMES
#include "sys/time.h"
#endif

void TensorMultiplierHermitian::product(const mwArray& input,mwArray& result){
  // Perform the contraction---TODO:CHECK!
#ifdef TESTINGTIMES
  struct timeval start,final;
  gettimeofday(&start,NULL);
#endif
  TensorMultiplier::product(input,result);
  mwArray aux=reshape(input,Indices(dd,betak,betakr));
  aux.permute(Indices(2,1,3));
  aux.reshape(Indices(betak,dd*betakr));
  // First contract tmpL with A
  //  aux=AL*aux;
  aux.multiplyLeft(ALh);
  aux.reshape(Indices(betab,alpha*dd,betakr));
  aux.permute(Indices(2,1,3));
  aux.reshape(Indices(alpha*dd,betab*betakr));
  //mwArray tmpres=AL*reshape(permute(A,Indices(2,1,3)),Indices(betak,dd*betakr));
  // Second, contract oper
  //  aux=Aop*aux;
  aux.multiplyLeft(Aoph);
  aux.reshape(Indices(du,al2,betab,betakr));
  aux.permute(Indices(1,3,2,4));
  //Finally, contract termR
  aux.reshape(Indices(du*betab,al2*betakr));
  //  result=aux*AR;
  aux.multiplyRight(ARh);
  aux.reshape(Indices(-1,1));
  result=.5*(result+aux);
#ifdef TESTINGTIMES
  counter++;
  gettimeofday(&final,NULL);
  timer+=(final.tv_sec-start.tv_sec)*1E3+(final.tv_usec-start.tv_usec)*1E-3;
#endif
}

TensorMultiplierHermitian::~TensorMultiplierHermitian(){
#ifdef TESTINGTIMES
  cout<<"Destroying TensorMultiplierHermitian, called "<<counter<<" times, with average time "
      <<timer*1./(counter)<<endl;
#endif
}
