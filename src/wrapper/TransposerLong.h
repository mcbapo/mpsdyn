#ifndef TRANSPOSERL_H
#define  TRANSPOSERL_H

#include "Transposer.h"

class TransposerLong:public Transposer{
  Indices perm;
 public:
  //Initializes the vector of new positions, as 1->
  TransposerLong(const shrt::Indices& dimensions,const shrt::Indices& permutation);
  ~TransposerLong();
  friend std::ostream& operator<<(std::ostream& os,const TransposerLong& transp);

 protected:
  Indices getKey() const{return Indices(dims,perm);}
  friend class TransposerMap;

};

#endif // TRANSPOSERL_H



