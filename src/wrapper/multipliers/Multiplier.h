#ifndef MULTIPLIER_H
#define MULTIPLIER_H

#include "mwArray.h"

#ifdef USING_PRIMME
extern "C" {
#include "primme.h"
}
#endif // USING_PRIMME

/** 
    \class <Multilplier> defines the interface of a class designed to
    contain a member function, which implements a particular
    matrix-vector contraction.
 */

class Multiplier{
 public:
  virtual ~Multiplier(){};
  virtual void product(const mwArray& vector,mwArray& result)=0;
  /** Return the dimension of the vector on which the Multiplier may
      act. */
  virtual int getSize() const =0;

  virtual void getFullTensor(mwArray& result){
    std::cout<<"ERROR: Calling getFullTensor on abstract Multiplier!"<<std::endl;
    exit(1);
  }

#ifdef USING_PRIMME 
  /** To adapt to the interface of PRIMME eigensolver, the UI of the
      product function has to be different. This method acts as wrapper
      to be able to use PRIMME instead of eigs. */
  virtual void product_primme(void* x,void* y){
    mwArray vector;
    mwArray result;
    int n=getSize();
    vector.setPointer(shrt::Indices(n,1),(double*)x);
    result.setPointer(shrt::Indices(n,1),(double*)y);
    //std::cout<<"About to call product on "<<vector<<" and "
    //     <<result<<std::endl;
    product(vector,result);
    //    std::cout<<"After product: result is "<<result<<std::endl;
  }

#endif
};
  
#endif // MULTIPLIER_H
