#ifndef HERMTENSORMULTIPLIER_H
#define HERMTENSORMULTIPLIER_H

#include "TensorMultiplier.h"

/** 
    \class <HermitianTensorMultiplier> is a TensorMultiplier (product 
    of three connected tensors) which forces the hermiticity of
    the tensor by computing always 1/2(T+T')*v. It is therefore 
    twice as expensive as the normal TensorMultiplier. */

class HermitianTensorMultiplier:
public TensorMultiplier{
  TensorMultiplier adjT; // the adjoint tensor
 public:
  HermitianTensorMultiplier(const mwArray& AL_,const mwArray& Aop_,const mwArray& AR_);
  /** Create an Hermitian version from an existing TensorMultiplier */
  HermitianTensorMultiplier(const TensorMultiplier& orig);
  ~HermitianTensorMultiplier();
  void product(const mwArray& input,mwArray& result);
  // Returns the effective dimension of the vector on which this Multiplier 
  // can act. This is done through the underlying TensorMultiplier
  //int getSize() const{return dd*betak*betakr;} //OR {return du*betab*betabr;}
};

#endif // HERMTENSORMULTIPLIER_H
