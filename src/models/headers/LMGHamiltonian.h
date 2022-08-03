/**
   \file LGMHamiltonian.h
   Hamiltonian for the Lipkin-Meshkov-Glick model 
   (which is as Ising with transverse magnetic field, with long range interactions).
   The class is used to implement the MPO for the Hamiltonian
   but no evolution operator is provided.

   \author Mari-Carmen Banuls
   \date 21/12/2016

*/


#ifndef LMGHAMIL_H
#define LMGHAMIL_H

#include "Hamiltonian.h"
#include "MPO.h"

typedef double (*tdepfun)(double);

/** 
    Implement the LMG Hamiltonian for fixed system site, N, and magnetic field, h
    \f[
    H=\frac{1}{4 N}\left(\sum_i \sigma_x^{[i]}\right)^2 + \frac{h}{2} \sum_i \sigma_z^{[i]} 
    \f]
*/
class LMGHamiltonian: 
public Hamiltonian {
 protected:
  int N; // length 
  int d=2; // physical dimension
  double h0; // initial values of the constants	
  mwArray Z; 
  MPO hamil; // MPO for the Hamiltonian

 public:
  LMGHamiltonian(int N,double h0);
  virtual ~LMGHamiltonian();
  bool hasHMPO() const {return true;}
  const MPO& getHMPO() const {return hamil;}  
  bool hasUMPO() const{return false;}

  //  /** Construct an MPo for the commutator with an operator, and return in in the argument.  */
  //void getCommutatorMPO(MPO& mpo);

 private:
  void initZ();
  void setHMPO(MPO& mpo,double h);


};


#endif // LMGHAMIL_H
