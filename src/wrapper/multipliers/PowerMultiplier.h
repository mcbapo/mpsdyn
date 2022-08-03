#ifndef PMULTIPLIER_H
#define PMULTIPLIER_H

#include "mwArray.h"
#include "Multiplier.h"

#ifdef USING_PRIMME
extern "C" {
#include "primme.h"
}
#endif // USING_PRIMME

/** 
    \class <PowerMultilplier> is the n-th power of a certain Multiplier.
 */

class PowerMultiplier: public Multiplier{
  int n;
  Multiplier& multi;
 public:
 PowerMultiplier(Multiplier& multi_,int k):n(k),multi(multi_){};
  void product(const mwArray& vector,mwArray& result){
    result=vector;
    for(int k=0;k<n;k++){
      mwArray aux(result);
      multi.product(aux,result);
    }
  }
  /** Return the dimension of the vector on which the Multiplier may
      act. */
  int getSize() const {return multi.getSize();}

  void getFullTensor(mwArray& result){multi.getFullTensor(result);}

/* #ifdef USING_PRIMME  */
  /** To adapt to the interface of PRIMME eigensolver, the UI of the
      product function has to be different. This method acts as wrapper
      to be able to use PRIMME instead of eigs. */
  /*   void product_primme(void* x,void* y){ */
/*     mwArray vector; */
/*     mwArray result; */
/*     int n=getSize(); */
/*     vector.setPointer(shrt::Indices(n,1),(double*)x); */
/*     result.setPointer(shrt::Indices(n,1),(double*)y); */
/*     //std::cout<<"About to call product on "<<vector<<" and " */
/*     //	     <<result<<std::endl; */
/*     product(vector,result); */
/*   } */

/* #endif */
};
  
#endif // MULTIPLIER_H
