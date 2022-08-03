#ifndef LCMULTIPLIER_H
#define LCMULTIPLIER_H

#include "Multiplier.h"
#include <vector>

/** 
    \class <LinearCombinationMultiplier> implements the Multiplier
    resulting from a list of Multipliers and complex coefficients.
*/

class LinearCombinationMultiplier:
public Multiplier{
  std::vector<complex_t> coeffs;
  std::vector<Multiplier*> multis;
  complex_t offset;
  int N; // explicit size of ALL 
 public:
  /** Construct a Multiplier as a linear combination of other
      Multipliers (each with its specific product method) plus
      possibly an independent term, so that the implemented operator
      is
      \f[
      A=\sum_k coeffs[k]*multis[k] + offset*Id
      \f]
  */
  LinearCombinationMultiplier(const std::vector<complex_t>& coeffs_,
			      const std::vector<Multiplier*>& multis_,
			      complex_t offset=ZERO_c);
  ~LinearCombinationMultiplier();
  int getSize() const{return N;}
  void product(const mwArray& input,mwArray& result);

};

#endif // LCMULTIPLIER_H
