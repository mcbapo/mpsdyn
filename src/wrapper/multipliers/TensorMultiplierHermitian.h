#ifndef TENSORMULTIPLIERHERMITIAN_H
#define TENSORMULTIPLIERHERMITIAN_H

#include "TensorMultiplier.h"

/** 
    \class <TensorMultiplierHermitian> implements the Multiplier
    interface with an optimized operation: product of three connected
    tensors, plus the assumption that the Tensor is Hermitian, i.e.,
    instead of the product it will return .5*(A+A^H) v */

class TensorMultiplierHermitian:
public TensorMultiplier{
  mwArray ALh;
  mwArray Aoph;
  mwArray ARh;
  /* int betab,alpha,betak; //termL size */
  /* int du,al1,dd,al2; // oper size */
  /* int betabr,alphar,betakr; //termR size */
#ifdef TESTINGTIMES
  long long timer;int counter;
#endif
 public:
  TensorMultiplierHermitian(const mwArray& AL_,const mwArray& Aop_,const mwArray& AR_);
  ~TensorMultiplierHermitian();
  void product(const mwArray& input,mwArray& result);
  // Returns the effective dimension of the vector on which this Multiplier 
  // can act.
  int getSize() const{return dd*betak*betakr;} //OR {return du*betab*betabr;}
  friend class Contractor;
  // protected:
 public:
  mwArray getLeftTerm(bool adjoint=false) const{return adjoint?ALh:AL;}
  mwArray getRightTerm(bool adjoint=false) const {return adjoint?ARh:AR;}
  mwArray getMiddleTerm(bool adjoint=false) const {return adjoint?Aoph:Aop;}
  virtual void getFullTensor(mwArray& result){
    //std::cout<<"In TensorMultiplierHermitian::getFullTensor"<<std::endl;
    TensorMultiplier::getFullTensor(result);
    result=.5*(result+Hconjugate(result));
  }
};

#endif // TENSORMULTIPLIER_H
