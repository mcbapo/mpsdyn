/**
   \file SchwingerHamiltonianGIGaussLaw.h
   Schwinger Hamiltonian with simplified version of gauge fields
   
   \author Stefan Kühn
   \date 13/03/2013
*/


#ifndef SCHWINGER_HAMILTONIAN_GI_GAUSSLAW_H
#define SCHWINGER_HAMILTONIAN_GI_GAUSSLAW_H

#include "MPS.h"
#include "Hamiltonian.h"

/** 
    Implementation of the Schwinger Hamiltonian in the spin
    representation with simplified gauge fields that violate the gauge invariance of the Hamiltonian
    By contrast to the other version the exponential for the link fields is calculted exactly which renders this version "more gauge invariant"

*/

/** List of inital states that can be automatically generated, |u> means spin up, |d> spin down |i> with integer i gives the flux between two spins*/
enum InitState{
  is_updown, //State of the form |u>|0>|d>|u>|0>|d>...
  is_downup, //State of the form |d>|0>|u>|d>|0>|u>...
  is_upup, //State of the form |u>|0>|u>|u>|0>|u>...
  is_downdown, //State of the form |d>|0>|d>|d>|0>|d>...
  is_string, //State with a string inside it
  is_pzoller_string //Stat with |u>|-1>|d>|-1>|u>|-1>|d>...
};

class SchwingerHamiltonianGIGaussLaw: 
public Hamiltonian{
private:
  //sigma_plus Op1 sigma_minus
  mwArray hopping1;
  //sigma_minus Op2 sigma_plus
  mwArray hopping2;
  //mwArray id L^2  id
  mwArray lterm;
  // 0.5*idpsz id id
  mwArray idpszleft;
  // id id 0.5*idpsz
  mwArray idpszright;  
  // idpsz id id
  mwArray firstidpszleft;  
  // id id idpsz
  mwArray lastidpszright;
  // bosonic identity as mpo to fill up
  mwArray id_bose_mpo;
  // fermionic identity as mpo to fill up
  mwArray id_fermi_mpo;
  //Flag if the terms were already calculated once and then stored
  bool terms_saved;
protected:
  int N_fermions;
  int N_total;
  int D;
  double mu,x,d_bose,d_fermi,e,lambda;
  MPO hamil;
  MPO charge;
  
   /** Construct the MPO representation of the Charge Operator*/
  void initChargeoperator(void);
  
   /** Generate some matrices used for single and few body MPOs to be saved in the object, therefore in the MPO the pointes to the matrices can be stored and do not memory*/
  void initOperators(void);

 public:
  /** Creates an instance of Schwinger Hamiltonian for a chain of N fermions, bosonic dimension d and parameters \f$ \mu, x, e \f$. The bond dimension is set to D. The gauge field on the links is simplified and represented by a \f$ d\times d \f$ Matrix and therefore the gauge invariance of the Hamiltonian is broken up. The full Hamiltonian is given by
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
  SchwingerHamiltonianGIGaussLaw(int N,int d,int D,double mu,double x,double e,double lambda);
  
  ~SchwingerHamiltonianGIGaussLaw();
  
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
  
  /** Generate the two MPOs for time evolution (odd-even-splitting) with desired timestep dt. The Hamiltonian for time evolution does not include a penalty term since the generator of the Gauss Law commutes with the Hamiltonian and therefore the state stays in the Gauss Law fulfilling sector of the Hilbertspace as long as the initial state has been in this sector. The two MPOs look like this:*/
  void constructEvolutionMPO(MPO &Uodd, MPO &Ueven, double dt);
  
  void constructEvolutionMPOefficient(MPO &Uodd, MPO &Ueven, double dt);
  
  void constructUoddMPO(MPO &Uodd, double dt);
  
  void constructUoddMPO(MPO &Uodd, double dt, double x1, double x2);
  
  void constructUevenMPO(MPO &Ueven, double dt);
  
  /** Generate an MPS (which fulfills GaussLaw \f$ L(n)-L(n-1) = \frac{1}{2}\left[\sigma^z(n)+(-1)^n\right]\f$) and therefore can be used as initial state to do time evolution and stay in the GaussLaw fulfilling subspace. The possibilites up to now are
  \f{eqnarray*}
  \mathrm{is\_updown} &= |\uparrow\rangle\otimes |0\rangle\otimes |\downarrow\rangle \\
  \mathrm{is\_downup} &= |\downarrow\rangle\otimes |0\rangle\otimes |\uparrow\rangle \\
  \mathrm{is\_upup} &= |\uparrow\rangle\otimes |0\rangle\otimes |\uparrow\rangle \\
  \mathrm{is\_downdown} &= |\downarrow\rangle\otimes |0\rangle\otimes |\downarrow\rangle
  \f}
   */
  MPS* constructInitalMPS(InitState state=is_updown,int extent=-1);
  
  /** Getter for x, pretty much self-explainatory*/
  double get_x() {return this->x;};
  /** Setter for x, pretty much self-explainatory*/
  void set_x(double xvalue){this->x = xvalue;};
  
  /** Getter for mu, pretty much self-explainatory*/
  double get_mu() {return this->mu;};
  /** Setter for mu, pretty much self-explainatory*/
  void set_mu(double muvalue){this->mu = muvalue;};
  
  /** In case a setter for parameters in the Hamiltonian is called, this function generates a new MPO for the Hamiltonian. As the Hamiltonian might not be needed the new calculation of the Hamiltonian MPO has to be explicitly started by the user*/
  void updateHMPO(void){initOperators();};
  
  
  
};

 

#endif
