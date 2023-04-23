/**
   \file SNGIGaussLawAngular.h
   Schwinger Hamiltonian with simplified version of the link field, this time GAUGE INVARIANT
   
   \author Stefan Kühn
   \date 10/04/2013
*/


#ifndef SHGI_GAUSSLAW_ANGULAR_H
#define SHGI_GAUSSLAW_ANGULAR_H

#include "MPS.h"
#include "Hamiltonian.h"

/** 
    Implementation of the Schwinger Hamiltonian in the spin representation with simplified link fields that violate the gauge invariance of the Hamiltonian. This version employs \f$ L(n) = a_n^\dagger a_n \f$ and chooses \f$ \theta(n) \f$ accordingly so that in the continuum the commutator fulfills  \f$ [\theta(n),L(m)]=i\delta_{n,m} \f$

*/

//I use the inital states at various headers, such that I would get a clash, that's why I secure it with an extra define
#ifndef InitState_E
#define InitState_E
/** List of inital states that can be automatically generated, |u> means spin up, |d> spin down |i> with integer i gives the flux between two spins*/
enum InitState{
  is_updown, //State of the form |u>|0>|d>|u>|0>|d>...
  is_downup, //State of the form |d>|0>|u>|d>|0>|u>...
  is_upup, //State of the form |u>|0>|u>|u>|0>|u>...
  is_downdown, //State of the form |d>|0>|d>|d>|0>|d>...
  is_string, //State with a string inside it
  is_pzoller_string //Stat with |u>|-1>|d>|-1>|u>|-1>|d>...
};
#endif

class SHGIGaussLawAngular: 
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
     V &=& \sum_n\left[\sigma^+(n)\exp(i\theta(n))\sigma1^-(n+1) + \sigma^+(n+1)\exp(-i\theta(n))\sigma1^-(n)\right] \\
    \f}
    For the Pauli-Matrices the following convention is chosen
    \f{eqnarray*}{
      \sigma^+(n) &=& \frac{1}{2}(\sigma_x + i\sigma_y) \\ 
      \sigma^-(n) &=& \frac{1}{2}(\sigma_x - i\sigma_y)
     \f}
    The operators \f$ \theta(n) \f$ and \f$ L(n)\f$ are given by
    \f{eqnarray*}
      \theta(n) &=& \frac{1}{N}\sum_{k,r,s}\frac{2\pi}{N}k\cdot\exp\left[-i\frac{2\pi}{N}k(s-r)\right]\varphi_s\langle\varphi_r,\cdot\rangle \\
      L(n) &=& diag(-j,-j+1,\dots,j-1,j),
    \f}
    where \f$ |\varphi_i\rangle \f$ are the basis states in which \f$ L(n) \f$ is diagonal.
    The term \f$ P \f$ acts as penalty term and ensures the Gauss-law, it is given by
    \f[ P=\sum_{n=1}^N\left(L(n)^2 - L(n-1)^2 -\frac{1}{2}\left[\sigma_3(n)+(-1)^n\right]\right)^2, \f]
    where it is assumed that \f$ L(N) = 0 = L(0) \f$. 
  */
  SHGIGaussLawAngular(int N,int d,int D,double mu,double x,double e,double lambda);
  
  ~SHGIGaussLawAngular();
  
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
  
  /** Generate an MPO for the penalty term in the Hamiltonian which ensures the Gauss Law. A reference to an MPO has to be given which will be overwritten. If a strength greater than zero is given explicitly, then it is taken as constant in front of the operator, otherwise the constant for the penalty out of the Hamiltonian is taken.*/
  void constructPenaltyMPO(MPO &mpo, double strength=-1.0);
  
  /** Generate the two MPOs for time evolution (odd-even-splitting) with desired timestep dt. The Hamiltonian for time evolution does not include a penalty term since the generator of the Gauss Law commutes with the Hamiltonian and therefore the state stays in the Gauss Law fulfilling sector of the Hilbertspace as long as the initial state has been in this sector.   */
  void constructEvolutionMPO(MPO &Uodd, MPO &Ueven, double dt);
  
  /** Basically the same as constructEvolutionMPO but parts of the tensor products are stored in the object and reused if the another MPO is constructed. Therefore if many MPOs for time evoultion are needed, this version is more efficient. For only one MPO this function produces more overhead*/
  void constructEvolutionMPOefficient(MPO &Uodd, MPO &Ueven, double dt);
  
   /** Basically the same as constructEvolutionMPOefficient but this time only the MPO for the odd sites is constructed. This might be useful in case one wants to do more elaborte stuff while ramping certain parameters of the Hamiltonian und so on*/
  void constructUoddMPO(MPO &Uodd, double dt);
  
  /** Overloaded version of constructUoddMPO(MPO &Uodd, double dt), instead of the x-value stored in the Hamiltonian it takes two x-values and constructs the time evolution operator with (x1+x2)/2. This is useful if one wants to do second order trotter evolution and ramp x adiabiatically. This routine allows one to take two odd timesteps together and therefore save computational effort.*/
  void constructUoddMPO(MPO &Uodd, double dt, double x1, double x2);
  
  /** Same as constructUoddMPO(MPO &Uodd, double dt) but this time for the even case*/
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
  
  /** MPO Ignacio suggested to monitor in which subspace we are staying, needs some improvement*/
  void constructUnitarySubspaceMPO(MPO &Uodd, MPO &Ueven, double dt, double strength=1.0);
  
  /** Projection operator to project the state into the Gauss Law fulfilling subspace*/
  void GaussProjector(MPO &mpo, int n);
  
  /** Projection operator to project the state into the Gauss Law fulfilling subspace*/
  void GaussProjector2(MPO &mpo, int n);
  
  /** Operator to enforce the unitary penalty*/
  void UnitaryPenaltyMPO(MPO &mpo, int n);
  
  /** Projector, which allows to project back into the subspace of total charge equal to zero. It should be used with care, because its bond dimension is half the number of spin sites and therefore can be quite large*/
  void ChargeProjectorMPO(MPO &mpo);
  
  /** The time evolution operator can by approximated by its taylor series up to first order \f$ U(t)=\exp(-i\Delta t W)\approx 1 - i\Delta t H \f$. This function generates a single MPO for that approximation of the time evolution operator. It is not unitary contrary to the Padé approximant of Contractor::approximateExponential(const MPO&,const MPS&,MPS&,double,int,bool,double ,double), but it is numerically cheaper. The bond dimension of the resulting MPO is four.*/
  void EvolutionMPO(MPO &mpo,double dt);
  
};

 

#endif
