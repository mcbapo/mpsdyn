#ifndef TENSORMULTIPLIERPROJM_H
#define TENSORMULTIPLIERPROJM_H

#include <vector>
#include "TensorMultiplierHermitian.h"

/** 
    \class <TensorMultiplierProj> implements the Multiplier interface
    as \class<TensorMultiplier>, but before returning the product,
    substracts the projection of the argument onto a given vector,
    with specified weight. This is used by Contractor to find the best
    MPS approximation to an excited state, by projecting out the
    ground state and other computed MPS for eigenvalues. */

class TensorMultiplierProjMulti:
public Multiplier{
  Multiplier* basic;
  std::vector<const mwArray*> vectors;
  std::vector<double> weights;
  complex_t offset;
  bool mine;
 public:
  /** Important: the effective vectors passed in vector vect_ are not
      copied, for the sake of economy, and only a pointer is stored,
      assuming that Contractor will not do anything to them until it is
      ready with TensorMultiplierProjMulti.  Moreover, the vectors
      come already in the bra version, as resulting from
      contractMbraId in Contractor. */
  TensorMultiplierProjMulti(const mwArray& AL_,const mwArray& Aop_,
			    const mwArray& AR_,
			    std::vector<double> weights,
			    const std::vector<mwArray>& vect,
			    bool hermitian=false,			    
			    complex_t offset=ZERO_c);
  /** Instead of creating my own Multiplier, I can receive one that
      comes prepared (f.i. a linear combination of multipliers) */
  TensorMultiplierProjMulti(Multiplier& multi,
			    std::vector<double> weights,
			    const std::vector<mwArray>& vect,
			    complex_t offset=ZERO_c);

  ~TensorMultiplierProjMulti();
  void product(const mwArray& input,mwArray& result);
  int getSize() const{return basic->getSize();} 
  void getFullTensor(mwArray& result); 
private:
  TensorMultiplierProjMulti(const TensorMultiplierProjMulti&);
  TensorMultiplierProjMulti& operator=(const TensorMultiplierProjMulti&);
};

#endif // TENSORMULTIPLIERPROJM_H
