#ifndef TENSORMULTIPLIERMPOHERM_H
#define TENSORMULTIPLIERMPOHERM_H

#include "MPOMultiplier.h"

/** 
    \class <MPOMultiplierHermitian> implements the Multiplier
    interface as a \class<MPOMultiplier>, but forcing Hermiticity,
    i.e., the product with a vector v will always return
    .5(M+M^{\dagger})|v>
 */

class MPOMultiplierHermitian:public MPOMultiplier{
 protected:
  MPOMultiplier mpoH;
 public:
  // Construct empty, with just the length
  MPOMultiplierHermitian(int length);
  MPOMultiplierHermitian(int length,const Operator** oplist);
  MPOMultiplierHermitian(const MPO& mpo);
  virtual ~MPOMultiplierHermitian();
  void product(const mwArray& input,mwArray& result);
  // // Returns the effective dimension of the vector on which this Multiplier 
  // // can act.
  // int getSize() const;
  /** This will fail if the size is too large, as it would construct a
      large matrix. It returns the Hermitian part */
  void getFullTensor(mwArray& result);
  void setOp(int pos,const Operator* op,bool myop=0);
  void setRotatedOp(int pos,const mwArray& op,const shrt::Indices& neworder,
		    bool conjugate=0);
};

#endif // TENSORMULTIPLIERMPO_H
