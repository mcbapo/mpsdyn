#ifndef DOPERATOR_H
#define DOPERATOR_H

#include "FoldedOperator.h"

/** 
    As \class<FoldedOperator>, only without the initial flip of left
    and right indices for one of the components.  The implementation
    just uses the FoldedOperator constructor which does not flip left
    and right dimensions of the left mwArray, and leaves the rest
    unchanged. Actually, it would be more logical that this was the
    mother class and FoldedOperator was derived, but it is kept this
    way for historical reasons.
  */

class DoubleOperator:
public FoldedOperator{
 public:
  /** Receive two arrays, corresponding to two operators that are
      gruped together now. They will be stored like that, and then
      used efficiently in the contractions.  The mwArrays must come
      with indices sorted as usual, up, left, down, right. The
      effective dimensions will be the product of the ones of the
      first argument times the ones of the second, although the big
      tensor is never trully constructed.
  */
  DoubleOperator(const mwArray& dataL,const mwArray& dataR);

  /** Destructor */
  ~DoubleOperator();

  /** Copy constructor */
  DoubleOperator(const DoubleOperator& op);

  /** Assignment operator */
  DoubleOperator& operator=(const DoubleOperator& op);

  void put(std::ostream& os) const;

};




#endif // DOPERATOR_H
