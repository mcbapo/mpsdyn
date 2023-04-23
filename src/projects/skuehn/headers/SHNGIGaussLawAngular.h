/**
   \file SHNGIGaussLawAngular.h
   Schwinger Hamiltonian with simplified version of the link field
   
   \author Stefan Kühn
   \date 20/03/2013
*/


#ifndef SHNGI_GAUSSLAW_ANGULAR_H
#define SHNGI_GAUSSLAW_ANGULAR_H

#include "MPS.h"
#include "Hamiltonian.h"

/** 
    Implementation of the Schwinger Hamiltonian in the spin representation with simplified link fields that violate the gauge invariance of the Hamiltonian. This version employs \f$ L(n) = a_n^\dagger a_n \f$ and chooses \f$ \theta(n) \f$ accordingly so that in the continuum the commutator fulfills  \f$ [\theta(n),L(m)]=i\delta_{n,m} \f$

*/

class SHNGIGaussLawAngular: 
public Hamiltonian{
protected:
  int N_fermions;
  int N_total;
  int D;
  unsigned int exp_order;
  double mu,x,d_bose,d_fermi,e,lambda;
  MPO hamil;
  MPO charge;
  
   /** Construct the MPO representation of the Charge Operator*/
  void initChargeoperator(void);
  
   /** Generate some matrices used for single and few body MPOs to be saved in the object, therefore in the MPO the pointes to the matrices can be stored and do not memory*/
  void initOperators(void);

 public:
  /** Creates an instance of Schwinger Hamiltonian for a chain of N fermions, bosonic dimension d and parameters \f$ \mu, x, e, \lambda \f$ and the order exp_order up to which the exponential is approximated. The bond dimension is set to D. The gauge field on the links is simplified and represented by a \f$ d\times d \f$ Matrix and therefore the gauge invariance of the Hamiltonian is broken up. The full Hamiltonian is given by
   \f[ W =\frac{2}{ag^2} H = W_0+xV + \lambda P \f]
   with 
   \f{eqnarray*}{ 
     W_0 &=& \sum_n L^2(n) + \frac{\mu}{2}\sum_n(-1)^n\sigma_3(n)+N\frac{\mu}{2}\\
     V &=& \sum_n\left[\sigma^+(n)(1+i\theta(n))\sigma1^-(n+1) + \sigma^+(n+1)(1-i\theta(n))\sigma1^-(n)\right] \\
    \f}
    For the Pauli-Matrices the following convention is chosen
    \f{eqnarray*}{
      \sigma^+(n) &=& \frac{1}{2}(\sigma_x + i\sigma_y) \\ 
      \sigma^-(n) &=& \frac{1}{2}(\sigma_x - i\sigma_y)
     \f}
    The operators \f$ \theta(n) \f$ and \f$ L(n)\f$ are represented with bosonic annihilation and creation operators \f$ a(n), a^\dagger(n) \f$
    \f{eqnarray*}
      \theta(n) &=& e\left(a(n)+a^\dagger(n)\right)\\
      L(n) &=& -\frac{i}{2e}\left(a(n) - a^\dagger(n)\right)
    \f}
    The term \f$ P \f$ acts as penalty term and ensures the Gauss-law, it is given by
    \f[ P=\sum_{n=1}^N\left(L(n)^2 - L(n-1)^2 -\frac{1}{2}\left[\sigma_3(n)+(-1)^n\right]\right)^2, \f]
    where it is assumed that \f$ L(N) = 0 = L(0) \f$. 
  */
  SHNGIGaussLawAngular(int N,int d,int D,double mu,double x,double e,double lambda);
  
  ~SHNGIGaussLawAngular();
  
  /** Return the MPO representation of the Hamiltonian*/
  const MPO& getHMPO() const {return hamil;}
  /** Return whether the MPO representation exists*/
  bool hasHMPO() const {return true;}
  
  /** Return the MPO representation of the Charge Operator \f$ \sum_{n=1}^N \sigma_z(n). \f$*/
  const MPO& getChargeMPO(){return charge;}
  
  /** Generate an MPO for L(n) since it might be needed for a special site it is not stored in the object itself, a reference to an MPO has to be given which will be overwritten. The parameter n specifies the bosonic site for which the MPO will be constructed, where n corresponds to a DMRG site, not a spin site.*/
  void constructLopMPO(MPO &mpo, int n);
  
  /** Generate an MPO for \f$ \frac{1}{2}[\sigma_z(n) + (-1)^n] \f$  since it might be needed for a special site it is not stored in the object itself, a reference to an MPO has to be given which will be overwritten. The parameter n specifies the fermionic site for which the MPO will be constructed, where n corresponds to a DMRG site, not a spin site.*/
  void constructGaussMPO(MPO &mpo, int n);
  
  /** Generate an MPO for \f$ \sigma_z \f$  since it might be needed for a special site it is not stored in the object itself, a reference to an MPO has to be given which will be overwritten. The parameter n specifies the fermionic site for which the MPO will be constructed, where n corresponds to a DMRG site, not a spin site.*/
  void constructsigmazMPO(MPO &mpo, int n);
  
  /** Generate an MPO for  the condensate fraction \f$ \frac{\sqrt{x}}{N}\sum_{n=1}^N\left(-1\right)^n(1+\sigma_z(n)) \f$. A reference to an MPO has to be given which will be overwritten.*/
  void constructCondensateMPO(MPO &mpo);
  
  /** Generate an MPO for the penalty term in the Hamiltonian which ensures the Gauss Law. A reference to an MPO has to be given which will be overwritten*/
  void constructPenaltyMPO(MPO &mpo);
  
};

 

#endif
