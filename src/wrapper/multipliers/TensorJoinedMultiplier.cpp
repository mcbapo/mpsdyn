
#include "TensorJoinedMultiplier.h"
#include "Indices.h"

using namespace std;
using namespace shrt;

TensorJoinedMultiplier::TensorJoinedMultiplier(const mwArray& AL_,const vector<const mwArray*>& Aops_,
				   const mwArray& AR_):
  AL(AL_),Aops(),AR(AR_),dimsL(AL_.getDimensions()),dimsR(AR_.getDimensions()),
  dims(),nrOp(Aops_.size()){
#ifdef TESTFUN
  cout<<"Creating a basic TensorJoinedMultiplier with AL="<<AL_.getDimensions()<<", Aops="<<Aops_.size()
      <<", AR="<<AR_.getDimensions()<<endl;
#endif
  betab=dimsL[0];betabr=dimsR[0];
  betak=AL.isScalar()?1:dimsL[nrOp+1];
  betakr=AR.isScalar()?1:dimsR[nrOp+1];
  for(int k=0;k<nrOp;k++){
    Aops.push_back(*Aops_[k]);
    Indices dimO=Aops[k].getDimensions(); // du,al1,dd,al2
    Aops[k].permute(Indices(2,3,1,4));Aops[k].reshape(Indices(dimO[1]*dimO[2],dimO[0]*dimO[3]));
    dims.push_back(dimO);
    if(k==nrOp-1) dd=dimO[2];
    if(k==0) du=dimO[0];
    //cout<<"Op "<<k<<":"<<dimO<<endl;
  }
  AL.reshape(Indices(-1,betak));
  AR.reshape(Indices(betabr,-1));AR.permute(Indices(2,1));
#ifdef CHECKDIMS
  bool wrong=false;  
  if((nrOp+2!=dimsL.size())||(nrOp+2!=dimsR.size())) wrong=true;
  if(!wrong){
    for(int k=0;k<nrOp;k++){
      Indices dimOk=Aops[k].getDimensions(); // du,al1,dd,al2
      if((dimOk[1]!=dimsL[nrOp-k])||(dimOk[3]!=dimsR[nrOp-k])) wrong=true;
    }
  }
  if(wrong||(betab!=betak)||(betabr!=betakr)){
    cout<<"Error: wrong dimensions constructing TensorJoinedMultiplier"
	<<" termL="<<AL<<", opers=";
    for(int k=0;k<nrOp;k++)
      cout<<Aops[k].getDimensions()<<", ";
    cout<<", AL="<<AL<<endl;
    exit(212);
  }
#endif
#ifdef TESTINGTIMES
  counter=0;timer=0;
#endif
}

#ifdef TESTINGTIMES
#include "sys/time.h"
#endif

TensorJoinedMultiplier::~TensorJoinedMultiplier(){
  dims.clear();Aops.clear();
}

void TensorJoinedMultiplier::permuteL(const Indices& newOrder){
  // AL was all x betak, without reordering. I reset original dims, permute, and reshape again
  AL.reshape(dimsL);
  AL.permute(newOrder);
  // need to permute the dimensions too
  dimsL=AL.getDimensions();
  betab=dimsL[0];
  betak=AL.isScalar()?1:dimsL[nrOp+1];
  AL.reshape(Indices(-1,betak));
}

void TensorJoinedMultiplier::permuteR(const Indices& newOrder){
  // AR was betabr x all transposed. I reset original dims, permute, and reshape again
  AR.permute(Indices(2,1));
  AR.reshape(dimsR);
  AR.permute(newOrder);
  // need to permute the dimensions too
  dimsR=AR.getDimensions();
  betabr=dimsR[0];
  betakr=AR.isScalar()?1:dimsR[nrOp+1];
  AR.reshape(Indices(betabr,-1));AR.permute(Indices(2,1));

}

void TensorJoinedMultiplier::product(const mwArray& input,mwArray& result_){
  // Perform the contraction---TODO:CHECK!
  //std::cout<<"In TensorJoinedMultiplier::product"<<std::endl;
#ifdef TESTINGTIMES
  struct timeval start,final;
  gettimeofday(&start,NULL);
#endif
  mwArray result=reshape(input,Indices(dd,betak,betakr));
  result.permute(Indices(2,1,3));
  result.reshape(Indices(betak,dd*betakr));
  // First contract tmpL with A
  //  aux=AL*aux;
  result.multiplyLeft(AL); // betab*dimsL x dd*betakr
  int nextR=betakr;
  for(int k=0;k<nrOp;k++){
    Indices dimO=dims[nrOp-k-1];
    result.reshape(Indices(-1,dimO[1]*dimO[2],nextR));
    result.permute(Indices(3,1,2));
    result.reshape(Indices(-1,dimO[1]*dimO[2]));
    result.multiplyRight(Aops[nrOp-1-k]);
    nextR=dimO[3];
  }
  // at the end,
  result.reshape(Indices(-1,betab*dd,nextR));
  result.permute(Indices(2,3,1));
  result.reshape(Indices(betab*dd,-1));
  result.multiplyRight(AR);
  result.reshape(Indices(betab,dd,betabr));
  result.permute(Indices(2,1,3));
  result.reshape(Indices(-1,1));
  result_=result;
#ifdef TESTINGTIMES
  counter++;
  gettimeofday(&final,NULL);
  timer+=(final.tv_sec-start.tv_sec)*1E3+(final.tv_usec-start.tv_usec)*1E-3;
#endif
}



void TensorJoinedMultiplier::getFullTensor(mwArray& result){
  std::cout<<"In TensorJoinedMultiplier::getFullTensor NOT YET SUPPORTED"<<std::endl;
  result=AL;
  result.reshape(Indices(-1,betak)); // betab*alphas*betak  
  result.permute(Indices(2,1)); // betak*betab x alphas

  int nextR=betak;
  for(int k=0;k<nrOp;k++){
    Indices dimO=dims[nrOp-k-1];
    if(k==0){
      result.reshape(Indices(-1,dimO[1],nextR));
      result.permute(Indices(3,1,2));
      result.reshape(Indices(-1,dimO[1]));
      mwArray aux=Aops[nrOp-1];
      aux.reshape(Indices(dimO[1],dimO[2]*dimO[0]*dimO[3]));
      result.multiplyRight(aux);
      result.reshape(Indices(-1,dimO[2],dimO[0]*dimO[3]));
      result.permute(Indices(2,1,3));
    }
    else{
      result.reshape(Indices(-1,dimO[1]*dimO[2],nextR));
      result.permute(Indices(3,1,2));
      result.reshape(Indices(-1,dimO[1]*dimO[2]));
      result.multiplyRight(Aops[nrOp-1-k]);
    }
    nextR=dimO[3];
  }
  // At the end, i have: all alphas x dd x betak x betab x du*last alpha
  result.reshape(Indices(-1,dd*betak*betab*dd,nextR));
  result.permute(Indices(2,3,1));
  result.reshape(Indices(dd*betak*betab*dd,-1));
  mwArray aux(AR); aux.reshape(Indices(-1,betakr*betabr));
  result.multiplyRight(aux);
  result.reshape(Indices(dd,betak,betab,dd,betakr,betabr));
  result.permute(Indices(4,3,6,1,2,5));
  result.reshape(Indices(dd*betab*betabr,dd*betak*betakr));
}
