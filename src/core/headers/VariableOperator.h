#ifndef VAROPERATOR_H
#define VAROPERATOR_H

#include "Operator.h"

/** A basic Operator with the added functionality to let the
    components of the mwArray being changed from outside. Instead of
    creating a new copy, the idea is to be able to set individual
    elements of the internal mwArray. */

class VariableOperator: public Operator{
 protected:
  //mwArray* data;
 public:
  VariableOperator(int d,int Dl,int Dr,int d2=0);
  /** Construct an operator, for a given array whose dimensions should
     be arranged as: du*dl*dd*dr. Saves an internal copy of the
     data.  */
  VariableOperator(const mwArray& oper);
  ~VariableOperator();

  /** I cannot just store the pointer, I have to copy it */
  void setData(const mwArray* oper);

  /** Change just one element in the stored tensor*/
  void setElement(complex_t value,const shrt::Indices& idx);

  /** Function that actually prints the Operator, and gets overloaded
      by daughter classes. */
  void put(std::ostream& os) const;

};

#endif // VAROPERATOR_H
