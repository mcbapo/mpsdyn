#ifndef TRANSPOSER_H
#define  TRANSPOSER_H

#include "Indices.h"
#include "mwArray.h"
#include <vector>
#include <fstream>

class Transposer{
// I could generalize to arbitrary permutation, but right now I will
// only do it for matrices.
  int count;
 protected:
  shrt::Indices dims;
  int d1,d2;
  //Indices perm;
  // Podria ser en lugar de esto, algo optimizado para la cache?
  std::vector<std::vector<int>*> cycles;
  int nrcycles;
 public:
  //Initializes the vector of new positions, as 1->
  Transposer(const shrt::Indices& dimensions,const shrt::Indices& permutation);
  ~Transposer();
  void transpose(mwArray& mwA,bool conj=false);
  friend std::ostream& operator<<(std::ostream& os,const Transposer& transp);

protected:
  // return just a pointe to the inner dimensions, to be used as
  // key by the map
  shrt::Indices* getDims(){return &dims;}
  friend class TransposerMap;

};

#endif // TRANSPOSER_H



