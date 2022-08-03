#ifndef TENSORMULTIPLIERPROJ_H
#define TENSORMULTIPLIERPROJ_H

#include "TensorMultiplierHermitian.h"

/** 
    \class <TensorMultiplierProj> implements the Multiplier interface
    as \class<TensorMultiplier>, but before returning the product,
    substracts the projection of the argument onto a given vector,
    with specified weight. This is used by Contractor to find the best
    MPS approximation to an excited state, by projecting out the
    ground state. */

class TensorMultiplierProj:
public Multiplier{
  TensorMultiplier* basic;
  mwArray vect;
  double weight;
 public:
  TensorMultiplierProj(const mwArray& AL_,const mwArray& Aop_,
		       const mwArray& AR_,
		       double weight,const mwArray& vect,bool hermitian=false);
  ~TensorMultiplierProj();
  void product(const mwArray& input,mwArray& result);
  int getSize() const{return basic->getSize();} 
  void getFullTensor(mwArray& result); 
private:
  TensorMultiplierProj(const TensorMultiplierProj& );
  TensorMultiplierProj& operator=(const TensorMultiplierProj& );
};

#endif // TENSORMULTIPLIERPROJ_H
