#include "TensorMultiplierProj.h"

#include "Indices.h"

using namespace std;
using namespace shrt;

TensorMultiplierProj::TensorMultiplierProj(const mwArray& AL_,
					   const mwArray& Aop_,
					   const mwArray& AR_,
					   double weight_,
					   const mwArray& vect_,
					   bool hermitian):
  vect(vect_),weight(weight_),basic(NULL){
  basic=hermitian?new TensorMultiplierHermitian(AL_,Aop_,AR_):new TensorMultiplier(AL_,Aop_,AR_);
}

TensorMultiplierProj::TensorMultiplierProj(const TensorMultiplierProj& ):
  vect(),weight(0),basic(NULL){};

TensorMultiplierProj& TensorMultiplierProj::operator=(const TensorMultiplierProj& ){
  return *this;
}


TensorMultiplierProj::~TensorMultiplierProj(){
  delete basic;
}


void TensorMultiplierProj::product(const mwArray& input,mwArray& result){
  // The first contribution comes from the inner TensorMultiplier action
  basic->product(input,result);
  // and now we substract something
  result=result-weight*(Hconjugate(vect)*input)*vect;
}

void TensorMultiplierProj::getFullTensor(mwArray& result){
  basic->getFullTensor(result);
  result=result-weight*vect*Hconjugate(vect);
}  
