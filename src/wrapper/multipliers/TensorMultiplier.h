#ifndef TENSORMULTIPLIER_H
#define TENSORMULTIPLIER_H

#include "Multiplier.h"

/** 
    \class <TensorMultiplier> implements the Multiplier interface with
    an optimized operation: product of three connected tensors */

class TensorMultiplier:
public Multiplier{
 protected:
  mwArray AL;
  mwArray Aop;
  mwArray AR;
  int betab,alpha,betak; //termL size
  int du,al1,dd,al2; // oper size
  int betabr,alphar,betakr; //termR size
#ifdef TESTINGTIMES
  long long timer;int counter;
#endif
 public:
  TensorMultiplier(const mwArray& AL_,const mwArray& Aop_,const mwArray& AR_);
  virtual ~TensorMultiplier();
  void product(const mwArray& input,mwArray& result);
  // Returns the effective dimension of the vector on which this Multiplier 
  // can act.
  int getSize() const{return dd*betak*betakr;} //OR {return du*betab*betabr;}
  friend class Contractor;
  // protected:
 public:
  mwArray getLeftTerm() const{return AL;}
  mwArray getRightTerm() const {return AR;}
  mwArray getMiddleTerm() const {return Aop;}

  virtual void getFullTensor(mwArray& result);

  // return the tensors with the original dimensions to be used by MPO
  mwArray getOriginalLeftTerm() const;
  mwArray getOriginalRightTerm() const;
  mwArray getOriginalMiddleTerm() const;


};

#endif // TENSORMULTIPLIER_H
