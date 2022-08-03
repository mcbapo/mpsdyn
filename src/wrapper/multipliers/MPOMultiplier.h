#ifndef TENSORMULTIPLIERMPO_H
#define TENSORMULTIPLIERMPO_H

#include "Multiplier.h"
#include <vector>
#include "MPO.h"

/** 
    \class <MPOMultiplier> implements the Multiplier interface similar
    to \class<TensorMultiplier>, but from an MPO, which, instead of
    being applied using TN techniques, is multiplied exactly tensor by
    tensor. 
 */

class MPOMultiplier:
public Multiplier,public MPO{
 protected:
  // For convenience I initialize a vector with the dimensions
  std::vector<int> dimsIn;
  int totalDim; // and their product
  bool init; // I need to know if it was correctly initialized or not
  std::vector<bool> opset;
#ifdef TESTINGTIMES
  long long timer;int counter;
#endif
 public:
  // Construct empty, with just the length
  MPOMultiplier(int length);
  MPOMultiplier(int length,const Operator** oplist);
  virtual ~MPOMultiplier();
  void product(const mwArray& input,mwArray& result);
  // Returns the effective dimension of the vector on which this Multiplier 
  // can act.
  int getSize() const;
  /** This will fail if the size is too larg, as it would construct a
      large matrix */
  void getFullTensor(mwArray& result){expandOper(*this,result);}
  void setOp(int pos,const Operator* op,bool myop=0);
  void setRotatedOp(int pos,const mwArray& op,const shrt::Indices& neworder,
		    bool conjugate=0);
 private:
  void initDims();
  bool allSet() const;
};

#endif // TENSORMULTIPLIERMPO_H
