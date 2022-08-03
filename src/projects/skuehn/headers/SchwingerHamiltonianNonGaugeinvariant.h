/**
   \file SchwingerHamiltonianNonGaugeinvariant.h
   Schwinger Hamiltonian with simplified version of gauge fields
   
   \author Stefan Kühn
   \date 16/01/2013
*/

#ifndef SCHWINGER_HAMILTONIAN_NON_GAUGEINVARIANT_H
#define SCHWINGER_HAMILTONIAN_NON_GAUGEINVARIANT_H

#include "MPS.h"
#include "Hamiltonian.h"

/** 
    Implementation of the Schwinger Hamiltonian in the spin
    representation with simplified gauge fields that violate the gauge invariance of the Hamiltonian
    \todo Complete implementation of the MPO-representation

*/

class SchwingerHamiltonianNonGaugeinvariant: 
public Hamiltonian{
protected:
  int N_fermions;
  int N_total;
  int D;
  double mu,x,mass,d,e;
  MPO hamil;
  MPO charge;
  mwArray Z; 
  
  /** Construct the MPO representation of the Schwinger Hamiltonian*/
  void initMPO(void);
  
   /** Construct the MPO representation of the Charge Operator*/
  void initChargeoperator(void);

 public:
  /** Creates an instance of Schwinger Hamiltonian for a chain of N fermions, bosonic dimension d and parameters mu, x, e. The bond dimension is chosen is set to D. The gauge field on the links is simplified and represented by a \f$ d\times d \f$ Matrix and therefore the gauge invariance of the Hamiltonian is broken up. The full Hamiltonian is given by
   \f[ W =\frac{2}{ag^2} H = W_0+xV \f]
   with 
   \f{eqnarray*}{ 
     W_0 &=& \sum_n L^2(n) + \frac{\mu}{2}\sum_n(-1)^n\sigma_3(n)+N\frac{\mu}{2}\\
     V &=& \sum_n\left[\sigma^+(n)(1+i\theta(n))\sigma1^-(n+1) + \sigma^+(n+1)(1-i\theta(n))\sigma1^-(n)\right] \\
    \f}
    For the Pauli-Matrices the following convention is chosen
    \f{eqnarray*}{
      \sigma^+(n) &=& \frac{1}{2}(\sigma_x - i\sigma_2) \\ 
      \sigma^-(n) &=& \frac{1}{2}(\sigma_x + i\sigma_2)
     \f}
    The operators \f$ \theta(n) \f$ and \f$ L(n)\f$ are represented with bosonic annihilation and creation operators \f$ a(n), a^\dagger(n) \f$
    \f{eqnarray*}
      \theta(n) &=& e\left(a(n)+a^\dagger(n)\right)\\
      L(n) &=& -\frac{i}{2e}\left(a(n) - a^\dagger(n)\right)
    \f}
  */
  SchwingerHamiltonianNonGaugeinvariant(int N,int d,int D,double mu,double x,double e);

  ~SchwingerHamiltonianNonGaugeinvariant();
  
  /** Return the MPO representation of the Hamiltonian*/
  const MPO& getHMPO(){return hamil;}
  /** Return whether the MPO representation exists*/
  bool hasHMPO() const {return true;}
  
  /** Return the MPO representation of the Charge Operator*/
  const MPO& getChargeMPO(){return charge;}
  
  /** Generate an MPO for L(n) since it might be needed for a special site it is not stored in the object itself, a reference to an MPO has to be given which will be overwritten. The parameter n specifies the bosonic site for which the MPO will be constructed, where n corresponds to a DMRG site, not a spin site.*/
  void constructLopMPO(MPO &mpo, int n);
  
  /** Generate an MPO for \f$ \frac{1}{2}[\sigma_z(n) + (-1)^n] \f$  since it might be needed for a special site it is not stored in the object itself, a reference to an MPO has to be given which will be overwritten. The parameter n specifies the fermionic site for which the MPO will be constructed, where n corresponds to a DMRG site, not a spin site.*/
  void constructGaussMPO(MPO &mpo, int n);
  
  /** Generate an MPO for \f$ \sigma_z \f$  since it might be needed for a special site it is not stored in the object itself, a reference to an MPO has to be given which will be overwritten. The parameter n specifies the fermionic site for which the MPO will be constructed, where n corresponds to a DMRG site, not a spin site.*/
  void constructsigmazMPO(MPO &mpo, int n);
  
  /** Generate an MPO for  the condensat fraction \f$ \frac{\sqrt{x}}{N}\sum_{n=1}^N\left(-1\right)^n(1+\sigma_z(n)) \f$  since it might be needed for a special site it is not stored in the object itself, a reference to an MPO has to be given which will be overwritten.*/
  void constructCondensateMPO(MPO &mpo);
  
};

 

#endif


