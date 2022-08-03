#ifndef TENSORJMULTIPLIER_H
#define TENSORJMULTIPLIER_H

#include "Multiplier.h"
#include <vector>

/** 
    \class <TensorJoinedMultiplier> implements the Multiplier interface with
    an optimized operation: product of three connected tensors */

class TensorJoinedMultiplier:
public Multiplier{
 protected:
  mwArray AL;
  std::vector<mwArray> Aops;
  mwArray AR;
  int betab,betak; //termL size
  int du,dd; // physical (output/input) dimension
  std::vector<int> alphas;
  //int du,al1,dd,al2;
  shrt::Indices dimsL,dimsR;
  std::vector<shrt::Indices> dims;// oper sizes
  int betabr,betakr; //termR size
  int nrOp;
#ifdef TESTINGTIMES
  long long timer;int counter;
#endif
 public:
  TensorJoinedMultiplier(const mwArray& AL_,const std::vector<const mwArray*>& Aops_,
			 const mwArray& AR_);
  virtual ~TensorJoinedMultiplier();
  void product(const mwArray& input,mwArray& result);
  // Returns the effective dimension of the vector on which this Multiplier 
  // can act.
  int getSize() const{return du*betab*betabr;} //OR {return du*betab*betabr;}
  friend class Contractor;
  // Manipulate the internal tensors to use a different order
  void permuteL(const shrt::Indices& newOrder);
  void permuteR(const shrt::Indices& newOrder);
  
  // protected:
 public:
  mwArray getLeftTerm() const{return AL;}
  mwArray getRightTerm() const {return AR;}
  //mwArray getMiddleTerm() const {return Aop;}

  virtual void getFullTensor(mwArray& result);
};

#endif // TENSORJMULTIPLIER_H
