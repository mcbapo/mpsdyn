#ifndef BASICMULTIPLIER_H
#define BASICMULTIPLIER_H

#include "Multiplier.h"

/** 
    \class <BasicMultiplier> implements the Multiplier interface with
    the most basic operation: straightforward product */

class BasicMultiplier:
public Multiplier{
  const mwArray& A;
 public:
  BasicMultiplier(const mwArray& A_);
  void product(const mwArray& vector,mwArray& result);
  int getSize() const{return A.getDimension(0);}
  void getFullTensor(mwArray& result){result=A;}
};

#endif // BASICMULTIPLIER_H
