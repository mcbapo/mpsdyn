#ifndef TENSORMULTIPLIEROFF_H
#define TENSORMULTIPLIEROFF_H

#include "TensorMultiplier.h"

/** 
    \class <TensorMultiplierOffset> as TensorMultiplier applies an
    optimized operation: product of three connected tensors, but then
    adds an offset (should be probably of the order of the largest
    eigenvalue) */

class TensorMultiplierOffset:
public TensorMultiplier{
  double offset;
 public:
  TensorMultiplierOffset(const mwArray& AL_,const mwArray& Aop_,const mwArray& AR_,double offset_=0.);
  ~TensorMultiplierOffset();
  double getOffset() const {return offset;}
  void product(const mwArray& input,mwArray& result);
  virtual void getFullTensor(mwArray& result);
};

#endif // TENSORMULTIPLIEROFF_H
