#ifndef TENSORMULTIPLIERHERMITIANOFF_H
#define TENSORMULTIPLIERHERMITIANOFF_H

#include "TensorMultiplierHermitian.h"

/** 
    \class <TensorMultiplierHermitianOffset> as
    TensorMultiplierHermitian applies an optimized operation: product
    of three connected tensors, and symmetrization, but then adds an
    offset (should be probably of the order of the largest
    eigenvalue) */

class TensorMultiplierHermitianOffset:
public TensorMultiplierHermitian{
  double offset;
 public:
  TensorMultiplierHermitianOffset(const mwArray& AL_,const mwArray& Aop_,const mwArray& AR_,double offset_=0.);
  ~TensorMultiplierHermitianOffset();
  void product(const mwArray& input,mwArray& result);
  void getFullTensor(mwArray& result){
//	  std::cout<<"In TensorMultiplierHermitianOffset::getFullTensor"<<std::endl;
	  TensorMultiplierHermitian::getFullTensor(result);
	  result=result+offset*identityMatrix(getSize());
  }
};

#endif // TENSORMULTIPLIERHERMITIANOFF_H
