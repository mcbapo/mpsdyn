/**    
   \file KSHamiltonian.h
   Class for truncated Kogut-Susskind with the rising/lowering operators of the electric field replaced by the ones presented in Zohar et. al Phys. Rev. A, 2013, 88, 023617
   
   \author Stefan Kühn
   \date 05/22/2016
*/

#ifndef KS2_Hamiltonian_H
#define KS2_Hamiltonian_H

#include "MPS.h"
#include "Hamiltonian.h"

//I use the initial states at various headers, such that I would get a clash, that's why I secure it with an extra define
#ifndef InitState_E
#define InitState_E
/** List of initial states that can be automatically generated, |u> means spin up, |d> spin down |i> with integer i gives the flux between two spins*/
enum InitState{
  is_updown, //State of the form |u>|0>|d>|u>|0>|d>...
  is_downup, //State of the form |d>|0>|u>|d>|0>|u>...
  is_upup, //State of the form |u>|0>|u>|u>|0>|u>...
  is_downdown, //State of the form |d>|0>|d>|d>|0>|d>...
  is_string, //State with a string inside it
  is_pzoller_string, //Start with |u>|-1>|d>|-1>|u>|-1>|d>...
  is_bgfield, // State of the form |d>|-1>|d>|-1>|u>|-1>|d>|-1>|u>...|-1>|u>
  is_bgfield_superpos // State of the form is_bgfield + is_updown
};
#endif


/** This class allows simulate Kogut-Susskind type Hamiltonians, where the operators for the links are replaced by something else. */


class KSHamiltonian: 
public Hamiltonian{
private:  
  //Array to store all fermionic operators
  mwArray FermiOperators;
  //Array to store all bose operators
  mwArray LinkOperators;
  
  // bosonic identity as mpo to fill up
  mwArray id_link_mpo;
  // fermionic identity as mpo to fill up
  mwArray id_fermi_mpo;
  
   
  
protected:
  //Number of sites
  int N_sites;  
  //Dimensionless mass
  double mu;
  //Coupling constant
  double x;
  //Background field
  double alpha;
  //Dimension for links and fermionic Hilbertspace
  int d_link,d_fermi;
  //Penalty strength
  double lambda;  
  
  //Hamiltonian MPO
  MPO hamiltonian;
  
   /** Generate some matrices used for single and few body MPOs to be saved in the object, therefore in the MPO the pointers to the matrices can be stored and do not consume memory*/
  void initOperators(void);
  
   /** After providing the general operators according to the chosen model, I can generate the actual Hamiltonian (with Mari Carmen's sign convention)*/
  void initHamiltonian(void);
  
  

 public:
  /** Creates an instance of the Kogut-Susskind-Hamiltonian. The form of the Hamiltonian is given by 
   \f[ H = M\sum_n (-1)^n\psi_n^\dagger\psi_n - \frac{i}{2a}\sum_n \left(\psi_n^\dagger L_{+,n}\psi_{n+1}+h.c.\right)+\frac{g^2}{2}\sum_n L_{z,n}^2 \f]
    After a Jordan Wigner transformation and rendering the result dimensionless we obtain the spin Hamiltonian which is implemented in this class. The formula is the following:
   \f[W = x\sum_n\left[\sigma^+_n L_{+,n}\sigma_{n+1}^- + \sigma^-_n L_{-,n}\sigma_{n+1}^+\right] + \frac{\mu}{2}\sum_n(-1)^n\left[ \sigma^z_{n+1}+1\right] + \sum_n L_{z,n}^2  \f]
   Since the groundstate search could bring us out of the sector of the states which fulfill Gauss' law, we add an additional penalty:
   \f[P(\lambda)= \lambda \sum_n\left( L_{z,n} - L_{z,n-1} -\frac{1}{2}\left(\sigma^z_n +(-1)^n\right)\right)^2 \f]
  */
  KSHamiltonian(int N,int d,double mu,double x,double alpha,double lambda);  
  
  ~KSHamiltonian();
  
  /** Return the MPO representation of the Hamiltonian*/
  const MPO& getHMPO() const {return hamiltonian;}
  /** Return whether the MPO representation exists*/
  bool hasHMPO() const {return true;}
  
  /** Getter for x, pretty much self-explanatory.*/
  double get_x() const {return this->x;};
  /** Setter for x, pretty much self-explanatory.*/
  void set_x(double xvalue){this->x = xvalue;};
  
  /** Getter for mu, pretty much self-explanatory*/
  double get_mu() const {return this->mu;};
  /** Setter for mu, pretty much self-explanatory.*/
  void set_mu(double muvalue){this->mu = muvalue;};
  
  /** Getter for lambda, pretty much self-explanatory*/
  double get_lambda() const {return this->lambda;};
  /** Setter for lambda, pretty much self-explanatory*/
  void set_lambda(double lambdavalue){this->lambda = lambdavalue;};
  
  
  /** In case a setter for parameters in the Hamiltonian is called, this function allows to generate a new MPO for the Hamiltonian. 
   *\warning As the Hamiltonian MPO might not be needed after a parameter has been changed,  the new calculation of the Hamiltonian MPO has to be explicitly started by the user.*/
  void updateHMPO(void);  
  
  /** Prepare some initial state, by default the strong coupling state is prepared*/
  void constructInitialMPS(MPS &mps, InitState state=is_updown) const;
  
  /** Generate an MPO for the penalty term \f$ \sum_{n=1}^{N}\left(L(n)-L(n-1)-\frac{1}{2}\left[\sigma^z(n)+(-1)^n\right]\right)^2 \f$ in the Hamiltonian which ensures the Gauss Law. A reference to an MPO has to be given which will be overwritten. If a strength greater than zero is given explicitly, then it is taken as constant in front of the operator, otherwise the constant \f$ \lambda \f$ for the penalty out of the Hamiltonian is taken.*/
  void constructPenaltyMPO(MPO &mpo, double strength=-1.0) const;
  
  /** Get the MPO representation for 
  \f[
  \Gamma^5 = \frac{\sqrt{x}}{N}\sum_{n=0}^{N-2}(-1)^n\left(\sigma_n^+L^+_n\sigma^-_{n+1} + \sigma_n^-L^-_n\sigma^+_{n+1} \right).
  \f]
  If factor_off is set to true, the prefactor \f$ \sqrt{x}/N \f$ is omitted.
  */  
  void constructGamma5MPO(MPO &mpo, bool prefactor_off=false) const;
  
  /** Get the MPO representation for 
  \f[
  \Gamma^\alpha = \frac{1}{N}\sum_{n=0}^{N-2}\left(L_{z,n}+\alpha \right).
  \f]
  If factor_off is set to true, the prefactor \f$ 1/N \f$ is omitted.
  */  
  void constructGammaAlphaMPO(MPO &mpo, bool prefactor_off=false) const;
  
  /** Get the MPO representation for 
  \f[
  \ \frac{\sqrt{x}}{N}\sum_{n=0}^{N-1}(-1)^n\left(\frac{1+\sigma_n^z}{2}\right).
  \f]
  If factor_off is set to true, the prefactor \f$ \sqrt{x}/N \f$ is omitted.
  */  
  void constructCondensateMPO(MPO &mpo, bool prefactor_off=false) const;
  
   /** Construct the discretized version of the momentum operator
   * \f[
   * \hat{P}=\int dx \psi^\dagger(x)i\partial\psi(x)
  \f]
  which is given by 
  \f[
    \hat{O}_p=-ix\sum_n\left[\sigma_n^-\sigma_{n+1}^z\sigma_{n+2}^+ - \mathrm{h.c.}\right]
  \f]
  after making it dimensionless like the Hamiltonian. The constant \f$ x \f$ is not taken into account and has to be added by hand afterwards
  */
  void constructMomentumMPO(MPO& Pmpo) const;
  
  
  /** Construct the MPO representation of the total charge operator
   \f[
   Q_\mathrm{tot}= \sum_{n=1}^{N} \frac{1}{2}\left(\sigma^z_n +(-1)^n\right) 
   \f]
  */
  void constructChargeMPO(MPO &mpo) const;
  
  /** Time evolution operator for the odd part of the odd even time evolution*/
  void constructUoddMPO(MPO &Uodd, double dt, bool imag_time=false);
  
  /** Time evolution operator for the even part of the odd even time evolution*/
  void constructUevenMPO(MPO &Uodd, double dt, bool imag_time=false);
  
};

#endif
