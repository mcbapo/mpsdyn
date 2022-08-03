/**
   \file CQEDBackgroundField.h
   Realisation of the fundamental Hamiltonian for 1+1d cQED as proposed in Ignacios review
   
   \author Stefan Kühn
   \date 31/05/2013
*/


#ifndef CQEDBackgroundField_H
#define CQEDBackgroundField_H

#include "MPS.h"
#include "Hamiltonian.h"

/** 
    Implementation of the Schwinger Hamiltonian in the spin representation with simplified link fields that violate the gauge invariance of the Hamiltonian. This version employs \f$ L(n) = a_n^\dagger a_n \f$ and chooses \f$ \theta(n) \f$ accordingly so that in the continuum the commutator fulfills  \f$ [\theta(n),L(m)]=i\delta_{n,m} \f$

*/

/** List of initial states that can be automatically generated, |u> means spin up, |d> spin down |i> with integer i gives the flux between two spins*/
enum InitState{
  is_updown, //State of the form |u>|0>|d>|u>|0>|d>...
  is_downup, //State of the form |d>|0>|u>|d>|0>|u>...
  is_upup, //State of the form |u>|0>|u>|u>|0>|u>...
  is_downdown, //State of the form |d>|0>|d>|d>|0>|d>...
  is_string, //State with a string inside it
  is_pzoller_string //Stat with |u>|-1>|d>|-1>|u>|-1>|d>...
};

class CQEDBackgroundField: 
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
  
   /** Generate some matrices used for single and few body MPOs to be saved in the object, therefore in the MPO the pointers to the matrices can be stored and do not memory*/
  void initOperators(void);

 public:
  /** Creates an instance of cQED from arXiv: 1303.5050 in 1+1d. The Hamiltonian is given by 
   \f[ H = M\sum_n (-1)^n\psi_n^\dagger\psi_n + \frac{\varepsilon}{\sqrt{l(l+1)}}\sum_n \left(\psi_n^\dagger L_{+,n}\psi_{n+1}+h.c.\right)+\frac{g^2}{8}\sum_n L_{z,n}^2\f]
   where 
   \f[L_+ = a^\dagger b,\quad L_- = b^\dagger a,\quad L_z = \frac{1}{2}(N_a - N_b),\quad l=\frac{1}{2}(N_a+N_b) = \frac{N_0}{2}.  \f]    
   and \f$ a\f$ and \f$ b \f$ are the annihilation operators for two different bosonic species. The total number of bosons per link, \f$ N_{0} \f$, is fixed to an even value. After a Jordan Wigner transformation and rendering the result dimensionless we obtain the spin Hamiltonian which is implemented in this class. The formula is the following:
   \f[H = \frac{8\varepsilon}{g^2\sqrt{l(l+1)}}\sum_n\left[\sigma^+_n iL_{+,n}\sigma_{n+1}^- + \sigma^-_n (-iL_{-,n})\sigma_{n+1}^+\right] + \frac{4M}{g^2}\sum_n(-1)^n\left[ \sigma^z_{n+1}+1\right] + \sum_n L_{z,n}^2  \f]
   Since the groundstate search could bring us out of the sector of the states which fulfill Gauss' law, we add an additional penalty:
   \f[ \lambda \sum_n\left( L_{z,n} - L_{z,n-1} -\frac{1}{2}\left(\sigma^z_n +(-1)^n\right)^2\right)\f]
  */
  CQEDBackgroundField(int N,int d,int D,double mu,double x,double e,double lambda);
  
  ~CQEDBackgroundField();
  
  /** Return the MPO representation of the Hamiltonian*/
  const MPO& getHMPO(){return hamil;}
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
  
  /** Generate an MPO for the penalty term \f$ \sum_{n=1}^N\left(L(n)-L(n-1)-\frac{1}{2}\left[\sigma^z(n)+(-1)^n\right]\right)^2 \f$ in the Hamiltonian which ensures the Gauss Law. A reference to an MPO has to be given which will be overwritten*/
  void constructPenaltyMPO(MPO &mpo);
  
  /** Generate the two MPOs for time evolution (odd-even-splitting) with desired time step dt. The Hamiltonian for time evolution does not include a penalty term since the generator of the Gauss Law commutes with the Hamiltonian and therefore the state stays in the Gauss Law fulfilling sector of the Hilbertspace as long as the initial state has been in this sector. The two MPOs look like this:
  \f{eqnarray*}{
  U_\mathrm{odd} &= x\left[\sigma^+_{n}\frac{iL_{+,n}}{\sqrt{l(l+1)}}\sigma^-_{n+1} +  \sigma^-_n \frac{-iL_{-,n}}{\sqrt{l(l+1)}}\sigma^+_{n+1} \right] + \frac{\mu}{2}(-1)\left[\sigma^z_n+1\right] + \frac{\mu}{2}(+1)\left[\sigma^z_{n+1}+1\right] + L_{z,n}^2\quad\quad\quad\mbox{for n odd} \\
  U_\mathrm{even} &= x\left[\sigma^+_{n}\frac{iL_{+,n}}{\sqrt{l(l+1)}}\sigma^-_{n+1} +  \sigma^-_n \frac{-iL_{-,n}}{\sqrt{l(l+1)}}\sigma^+_{n+1} \right] + \frac{\mu}{2}(+1)\left[\sigma^z_n+1\right] + \frac{\mu}{2}(-1)\left[\sigma^z_{n+1}+1\right] + L_{z,n}^2\quad\quad\quad\mbox{for n odd} \\
  \f}
  */
  void constructEvolutionMPO(MPO &Uodd, MPO &Ueven, double dt);
  
  /** Basically the same as constructEvolutionMPO but parts of the tensor products are stored in the object and reused if the another MPO is constructed. Therefore if many MPOs for time evolution are needed, this version is more efficient. For only one MPO this function produces more overhead*/
  void constructEvolutionMPOefficient(MPO &Uodd, MPO &Ueven, double dt);
  
  /** Basically the same as constructEvolutionMPOefficient but this time only the MPO for the odd sites is constructed. This might be useful in case one wants to do more elaborte stuff while ramping certain parameters of the Hamiltonian and so on*/
  void constructUoddMPO(MPO &Uodd, double dt);
  
  /** Overloaded version of constructUoddMPO(MPO &Uodd, double dt), instead of the x-value stored in the Hamiltonian it takes two x-values and constructs the time evolution operator with (x1+x2)/2. This is useful if one wants to do second order trotter evolution and ramp x adiabiatically. This routine allows one to take two odd time steps together and therefore save computational effort.*/
  void constructUoddMPO(MPO &Uodd, double dt, double x1, double x2);
  
  /** Same as constructUoddMPO(MPO &Uodd, double dt) but this time for the even case*/
  void constructUevenMPO(MPO &Ueven, double dt);
  
  /** Construct the MPO for \f$ \sigma^+_{n}\frac{iL_{+,n}}{\sqrt{l(l+1)}}\sigma^-_{n+1} \f$ which flips a spin-down spin-up pair*/
  void constructFlipduMPO(MPO &MPO, int n);
  
  /** Construct the MPO for \f$ \sigma^-_{n}\frac{iL_{-,n}}{\sqrt{l(l+1)}}\sigma^+_{n+1} \f$ which flips a spin-up spin-down pair*/
  void constructFlipudMPO(MPO &MPO, int n);
  
  /** Function to generate a string on top of the state mps at specified starting position pos_start and with specified length. This routine does not check, if the given parameters lead to a valid string, the user has to pay attention that they are reasonable. Furthermore if the value of pos_start is negative the routine automatically places the string approximately in the middle of the MPS.*/
  void makeString(MPS &mps, int pos_start, int length);
  
  /** Generate an MPS (which fulfills GaussLaw \f$ L(n)-L(n-1) = \frac{1}{2}\left[\sigma^z(n)+(-1)^n\right]\f$) and therefore can be used as initial state to do time evolution and stay in the GaussLaw fulfilling subspace. The possibilities up to now are
  \f{eqnarray*}
  \mathrm{is\_updown} &= |\uparrow\rangle\otimes |0\rangle\otimes |\downarrow\rangle \\
  \mathrm{is\_downup} &= |\downarrow\rangle\otimes |0\rangle\otimes |\uparrow\rangle \\
  \mathrm{is\_upup} &= |\uparrow\rangle\otimes |0\rangle\otimes |\uparrow\rangle \\
  \mathrm{is\_downdown} &= |\downarrow\rangle\otimes |0\rangle\otimes |\downarrow\rangle
  \f}
   */
  MPS* constructInitalMPS(InitState state=is_updown,int extent=-1);
  
  /** Getter for x, pretty much self-explanatory*/
  double get_x() {return this->x;};
  /** Setter for x, pretty much self-explanatory*/
  void set_x(double xvalue){this->x = xvalue;};
  
  /** Getter for mu, pretty much self-explanatory*/
  double get_mu() {return this->mu;};
  /** Setter for mu, pretty much self-explanatory*/
  void set_mu(double muvalue){this->mu = muvalue;};
  
  /** In case a setter for parameters in the Hamiltonian is called, this function generates a new MPO for the Hamiltonian. As the Hamiltonian might not be needed the new calculation of the Hamiltonian MPO has to be explicitly started by the user*/
  void updateHMPO(void){initOperators();};
};

 

#endif
