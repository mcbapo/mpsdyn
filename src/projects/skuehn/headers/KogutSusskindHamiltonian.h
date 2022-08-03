/**    
   \file KogutSusskindHamiltonian.h
   Class for Kogut-Susskind type Hamiltonians
   
   \author Stefan Kühn
   \date 17/02/2014
*/

#ifndef KS_Hamiltonian_H
#define KS_Hamiltonian_H

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
  is_pzoller_string //Start with |u>|-1>|d>|-1>|u>|-1>|d>...
};
#endif

//List of all the currently supported models (these define how the link operators look later on)
enum Model{
	//cQED model from proposal of Erez, Benni and Ignacio
	cQED,
	//cQED model from proposal of Erez, Benni and Ignacio with some simple noise terms
	cQEDnoise,
	//Effective model wich in the low energy limit corresponds to the cQED model of Erez and Benni
	cQEDeffective,
	//Zn model, with the usual Gauss Law as penalty for variational ground state search
	Zn,
	//Zn model, with the a Gauss Law, which respects the cyclic nature of the operators
	Zncgl,
	//Zn model, with the a Gauss Law, which respects the cyclic nature of the operators and some noise inside
	Zncglnoise,
	//Real Zn, meaning that the kinetic term is the one of a real Zn theory (not what is used in Zn, Zncgl and Zncglnoise)
	realZn
};

/** This class allows simulate Kogut-Susskind type Hamiltonians, where the operators for the links are replaced by something else. */


class KogutSusskindHamiltonian: 
public Hamiltonian{
private:  
  //Pauli Matrices
  mwArray sigma_x,sigma_y,sigma_z,sigma_plus,sigma_minus;
  //The L-operators
  mwArray Lp,Lm,Lz,Lzsquare;
  //Identities
  mwArray id_fermi,id_bose;
  //id_fermi+sigma_z (often needed, therefore separately provided)
  mwArray idpsz;
  
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
  //Terms for the noise
  mwArray Lpnoise;
  mwArray Lmnoise;
  //Flag if the terms were already calculated once and then stored
  bool terms_saved;
   
  
protected:
  //Number of fermions
  int N_fermions;
  //Total number of sites (2*N_fermions-1)
  int N_total;
  //Bond dimension
  int D;
  //Mass
  double mu;
  //Coupling constant
  double x;
  //Dimension for bosonic and fermionic Hilbertspace
  int d_bose,d_fermi;
  //Penalty strength
  double lambda;
  //Which model?
  Model model;
  //Strength of the noise terms in case cQEDnoise model is chosen
  double noise_strength;
  //Strength of the interaction generating factor
  double lxprefactor;
  
  //Hamiltonian MPO
  MPO hamil;
  
   /** Generate some matrices used for single and few body MPOs to be saved in the object, therefore in the MPO the pointers to the matrices can be stored and do not consume memory*/
  void initOperators(void);
  
   /** After providing the general operators according to the chosen model, I can generate the actual Hamiltonian*/
  void initHamiltonian(void);
  
  /** For the Zn model with the Gauss Law which respects the cyclicity of the model, I need a Hamiltonian which is a bit different*/
  void initZncglHamiltonian(void);
  
   /** For the effective cQED model I need a Hamiltonian which is a bit different*/
  void initeffectiveHamiltonian(void);

 public:
  /** Creates an instance of the Kogut-Susskind-Hamiltonian. The form of the Hamiltonian is given by 
   \f[ H = M\sum_n (-1)^n\psi_n^\dagger\psi_n + x\sum_n \left(\psi_n^\dagger L_{+,n}\psi_{n+1}+h.c.\right)+\frac{g^2}{8}\sum_n L_{z,n}^2\f]
    After a Jordan Wigner transformation and rendering the result dimensionless we obtain the spin Hamiltonian which is implemented in this class. The formula is the following:
   \f[H = x\sum_n\left[\sigma^+_n L_{+,n}\sigma_{n+1}^- + \sigma^-_n L_{-,n}\sigma_{n+1}^+\right] + \frac{4M}{g^2}\sum_n(-1)^n\left[ \sigma^z_{n+1}+1\right] + \sum_n L_{z,n}^2  \f]
   Since the groundstate search could bring us out of the sector of the states which fulfill Gauss' law, we add an additional penalty:
   \f[ \lambda \sum_n\left( L_{z,n} - L_{z,n-1} -\frac{1}{2}\left(\sigma^z_n +(-1)^n\right)\right)^2\f]
   The exact form of the operators \f$ L_{+,n},L_{-,n} \f$ is determined by the model one wants to look at.
   
   If the effective model is chosen, the Hamiltonian (in spin formulation) is the following:
   \f{eqnarray*}{
     H &=& H_b + H_f + H_E + H_M + H_\lambda \\
     H_f &=& x\sum_n\left[\sigma_n^+\sigma_{n+1}^- + \sigma_n^-\sigma_{n+1}^+\right]\\
     H_M &=& \frac{\mu}{2}\sum_n(-1)^n\left(\sigma_n^z+1\right)\\
     H_E &=& \sum_n L_{z,n}^2 \\ 
     H_b &=& c\sum_n\frac{L_{x,n}}{\sqrt{l(l+1)}}\\
     H_\lambda &=& \lambda\sum_n\left[L_{z,n}-L_{z,n-1}-\frac{1}{2}\left(\sigma^z_n+(-1)^n\right)\right]^2.
     \f}
   This Hamiltonian has the a low energy limit, which coincides with the above form of \f$ H \f$.
   
   If "realZn" is chosen, the kinetic term for the gauge field is replaced by the one of a real Zn lattice gauge theory which gives as Hamiltonian
   \f[ H =  x\sum_n\left[\sigma^+_n L_{+,n}\sigma_{n+1}^- + \sigma^-_n L_{-,n}\sigma_{n+1}^+\right] + \frac{\mu}{2}\sum_n(-1)^n\left[ \sigma^z_{n+1}+1\right] + \sum_n \left[1-\cos\left(\frac{2\pi}{d}L_{z,n}\right)\right]\f]
  */
  KogutSusskindHamiltonian(int N,int d,int D,double mu,double x,double lambda, Model model,double noise_strength=1.0);
  
  ~KogutSusskindHamiltonian();
  
  /** Return the MPO representation of the Hamiltonian*/
  const MPO& getHMPO() const {return hamil;}
  /** Return whether the MPO representation exists*/
  bool hasHMPO() const {return true;}
  
  /** Return the MPO representation of operator for the total \f$ S_z \f$, \f$ \sum_{n=1}^N \sigma_z(n). \f$ */
  void constructTotalSzMPO(MPO &mpo) const;
  
  /** Return the MPO representation of the Charge Operator \f$ \frac{1}{2}\sum_{n=1}^N (\sigma_z(n)+(-1)^n) \f$ */
  void constructChargeMPO(MPO &mpo) const;
  
  /** Generate an MPO for L(n) since it might be needed for a special site it is not stored in the object itself, a reference to an MPO has to be given which will be overwritten. The parameter n specifies the bosonic site for which the MPO will be constructed, where n corresponds to a DMRG site, not a spin site.*/
  void constructLopMPO(MPO &mpo, int n) const;
  
  /** Generate an MPO for L²(n) since it might be needed for a special site it is not stored in the object itself, a reference to an MPO has to be given which will be overwritten. The parameter n specifies the bosonic site for which the MPO will be constructed, where n corresponds to a DMRG site, not a spin site.*/
  void constructLopsquareMPO(MPO &mpo, int n) const;
  
  /** Generate an MPO for \f$ \sum_nL(n) \f$. A reference to an MPO has to be given which will be overwritten.*/
  void constructLtotMPO(MPO &mpo) const;
  
  /** Generate an MPO for \f$ \frac{1}{2}[\sigma_z(n) + (-1)^n] \f$  since it might be needed for a special site it is not stored in the object itself, a reference to an MPO has to be given which will be overwritten. The parameter n specifies the fermionic site for which the MPO will be constructed, where n corresponds to a DMRG site, not a spin site.*/
  void constructGaussMPO(MPO &mpo, int n) const;
  
  /** Generate an MPO for \f$ \sigma_z(n) \f$  since it might be needed for a special site it is not stored in the object itself, a reference to an MPO has to be given which will be overwritten. The parameter n specifies the fermionic site for which the MPO will be constructed, where n corresponds to a DMRG site, not a spin site.*/
  void constructsigmazMPO(MPO &mpo, int n) const;
  
  /** Generate an MPO for  the condensate fraction \f$ \frac{\sqrt{x}}{N}\sum_{n=1}^N\left(-1\right)^n(1+\sigma_z(n)) \f$. A reference to an MPO has to be given which will be overwritten.*/
  void constructCondensateMPO(MPO &mpo) const;
  
  /** Generate an MPO for the penalty term \f$ \sum_{n=1}^N\left(L(n)-L(n-1)-\frac{1}{2}\left[\sigma^z(n)+(-1)^n\right]\right)^2 \f$ in the Hamiltonian which ensures the Gauss Law. A reference to an MPO has to be given which will be overwritten. If a strength greater than zero is given explicitly, then it is taken as constant in front of the operator, otherwise the constant for the penalty out of the Hamiltonian is taken.*/
  void constructPenaltyMPO(MPO &mpo, double strength=-1.0) const;
  
  /** Generate an MPO for the penalty term \f$ \sum_{n=1}^N\left(\exp\left(i\frac{2\pi}{N}(L(n)-L(n-1)-\frac{1}{2}\left[\sigma^z(n)+(-1)^n\right]\right)-1\right)\cdot \left(\exp\left(i\frac{2\pi}{N}(L(n)-L(n-1)-\frac{1}{2}\left[\sigma^z(n)+(-1)^n\right]\right)-1\right)^\dagger\f$ in the Hamiltonian which ensures the Gauss Law modulo the dimension of the bosonic Hilbertspaces. A reference to an MPO has to be given which will be overwritten. If a strength greater than zero is given explicitly, then it is taken as constant in front of the operator, otherwise the constant for the penalty out of the Hamiltonian is taken.*/
  void constructCglPenaltyMPO(MPO &mpo, double strength=-1.0) const;
  
  /** \f{eqnarray*}{
  U_\mathrm{odd} &= x\left[\sigma^+_{n}\frac{iL_{+,n}}{\sqrt{l(l+1)}}\sigma^-_{n+1} +  \sigma^-_n \frac{-iL_{-,n}}{\sqrt{l(l+1)}}\sigma^+_{n+1} \right] + \frac{\mu}{2}(-1)\left[\sigma^z_n+1\right] + \frac{\mu}{2}(+1)\left[\sigma^z_{n+1}+1\right] + L_{z,n}^2\quad\quad\quad\mbox{for n odd} \\
  U_\mathrm{even} &= x\left[\sigma^+_{n}\frac{iL_{+,n}}{\sqrt{l(l+1)}}\sigma^-_{n+1} +  \sigma^-_n \frac{-iL_{-,n}}{\sqrt{l(l+1)}}\sigma^+_{n+1} \right] + \frac{\mu}{2}(+1)\left[\sigma^z_n+1\right] + \frac{\mu}{2}(-1)\left[\sigma^z_{n+1}+1\right] + L_{z,n}^2\quad\quad\quad\mbox{for n odd} \\
  \f}
  Generate the two MPOs for time evolution (odd-even-splitting) with desired time step dt. The Hamiltonian for time evolution does not include a penalty term since the generator of the Gauss Law commutes with the Hamiltonian and therefore the state stays in the Gauss Law fulfilling sector of the Hilbertspace as long as the initial state has been in this sector. The two MPOs look like in the formula above.
  *\warning This routine is deprecated and only acts as a wrapper which calls \ref constructUoddMPO and \ref constructUevenMPO to ensure compatibility with older versions
  */
  void constructEvolutionMPO(MPO &Uodd, MPO &Ueven, double dt);
  
  /** Basically the same as constructEvolutionMPOefficient but this time only the MPO for the odd sites is constructed. This might be useful in case one wants to do more elaborate stuff while ramping certain parameters of the Hamiltonian and so on.*/
  void constructUoddMPO(MPO &Uodd, double dt);
  
  /** Overloaded version of constructUoddMPO(MPO &Uodd, double dt), instead of the x-value stored in the Hamiltonian it takes two x-values and constructs the time evolution operator with (x1+x2)/2. This is useful if one wants to do second order trotter evolution and ramp x adiabiatically. This routine allows one to take two odd time steps together and therefore save computational effort.*/
  void constructUoddMPO(MPO &Uodd, double dt, double x1, double x2);
  
  /** Same as constructUoddMPO(MPO &Uodd, double dt) but this time for the even case.*/
  void constructUevenMPO(MPO &Ueven, double dt);
  
  /** Construct the MPO for \f$ \sigma^+_{n}\frac{iL_{+,n}}{\sqrt{l(l+1)}}\sigma^-_{n+1} \f$ which flips a spin-down spin-up pair*/
  void constructFlipduMPO(MPO &MPO, int n);
  
  /** Construct the MPO for \f$ \sigma^-_{n}\frac{iL_{-,n}}{\sqrt{l(l+1)}}\sigma^+_{n+1} \f$ which flips a spin-up spin-down pair.*/
  void constructFlipudMPO(MPO &MPO, int n);
  
  /** Function to generate a string on top of the state mps at specified starting position pos_start and with specified length. This routine does not check, if the given parameters lead to a valid string, the user has to pay attention that they are reasonable. Furthermore if the value of pos_start is negative the routine automatically places the string approximately in the middle of the MPS.*/
  void makeString(MPS &mps, int pos_start, int length);
  
  /** Generate an MPS (which fulfills GaussLaw \f$ L(n)-L(n-1) = \frac{1}{2}\left[\sigma^z(n)+(-1)^n\right]\f$) and therefore can be used as initial state to do time evolution and stay in the GaussLaw fulfilling subspace. The possibilities up to now are
  \f{eqnarray*}
  \mathrm{is\_updown} &= |\uparrow\rangle\otimes |0\rangle\otimes |\downarrow\rangle\dots \\
  \mathrm{is\_downup} &= |\downarrow\rangle\otimes |0\rangle\otimes |\uparrow\rangle\dots \\
  \mathrm{is\_upup} &= |\uparrow\rangle\otimes |0\rangle\otimes |\uparrow\rangle \dots\\
  \mathrm{is\_downdown} &= |\downarrow\rangle\otimes |0\rangle\otimes |\downarrow\rangle\dots
  \f}
   */
  MPS* constructInitalMPS(InitState state=is_updown,int extent=-1);
  
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
  
  /** Getter for noise constant, pretty much self-explanatory*/
  double get_noise() const {return this->noise_strength;};
  /** Setter for noise constant, pretty much self-explanatory*/
  void set_noise(double noisevalue){this->noise_strength = noisevalue;};
  
  /** Getter for noise prefactor of the Lx-term, only relevant for the effective cQED model*/
  double get_lxprefactor() const {return this->lxprefactor;};
  /** Setter for the prefactor of the Lx-term, only relevant for the effective cQED model*/
  void set_lxprefactor(double value){this->lxprefactor = value;};
  
  /** In case a setter for parameters in the Hamiltonian is called, this function allows to generate a new MPO for the Hamiltonian. 
   *\warning As the Hamiltonian MPO might not be needed after a parameter has been changed,  the new calculation of the Hamiltonian MPO has to be explicitly started by the user.*/
  void updateHMPO(void);
  
  /** The time evolution operator can by approximated by its taylor series up to first order \f$ U(t)=\exp(-i\Delta t W)\approx 1 - i\Delta t H \f$. This function generates a single MPO for that approximation of the time evolution operator. It is not unitary contrary to the Padé approximant of Contractor::approximateExponential(const MPO&,const MPS&,MPS&,double,int,bool,double ,double), but it is numerically cheaper. The bond dimension of the resulting MPO is four.
   \warning Deprecated, works only for the cQED model!*/
  void EvolutionMPO(MPO &mpo,double dt);
  
   /** The time evolution operator can by approximated by its taylor series up to first order \f$ U(t)=\exp(-i\Delta t W)\approx 1 - i\Delta t H \f$. This function generates a single MPO for that approximation of the time evolution operator. It is not unitary contrary to the Padé approximant of Contractor::approximateExponential(const MPO&,const MPS&,MPS&,double,int,bool,double ,double), but it is numerically cheaper. The bond dimension of the resulting MPO is four.
    \warning Right now, it only works for the cQED, Zn and Zncyclic model, noisy models or the effective cQED models are not supported*/
  void EvolutionMPO2(MPO &mpo,double dt);
  
  
  /** Projection operator to project the state into the Gauss Law fulfilling subspace. This version is local and a spin site has to be given, which will be projected back into this subspace.
   */
  void GaussProjector(MPO &mpo, int n);
  
  /** Projection operator to project the state into the Gauss Law fulfilling subspace. This version is global and projects at every spin site at once. Basically it is the same as \f$ \sum_n \f$ GaussProjector2(MPO &mpo, int n) , however it is more efficient than doing this sum since it is only a single MPO.
   */
  void GaussProjector(MPO &mpo);
  
  /** Function to generate a string operator \f$ S_{nl} = -i\sigma^+_nL^+_n\dots L^+_{n+l-1}\sigma^-_{n+l} \f$ of length \f$ l \f$ starting at site \f$ n=1,\dots,N_\mathrm{fermions}\f$. If n is odd, the adjoint operator is constructed (such that the operators generate a string on the strong coupling vacuum and does not annihilate it). If the flag for the adjoint operator is set to true, then the adjoint operator is constructed. In case the value for n<=0 then the string is automatically placed approximately in the middle of the chain.*/
  void getStringOperator(MPO &mpo, int l, int n=0, bool adjoint=false) const;
  
  /** Function to generate a superposition of mesons (strings of length 1) peaked around site \f$ n=1,\dots,N_\mathrm{fermions}-1\f$ with momentum \f$ k \f$. This is done by forming a superposition \f$ \sum_n \exp(-I\cdot k\cdot n)\exp\left(-\frac{1}{2w^2}(n-n_0)^2\right) \f$. Here Mesons are only taken to start at odd sites.*/
  void getMesonOperator(MPO &mpo, int n, double k=1.0, double width=-1.0) const;
  
  /** Function to generate a superposition of mesons and antimesions (strings of length 1) peaked around site \f$ n=1,\dots,N_\mathrm{fermions}-1\f$ with momentum \f$ k \f$. This is done by forming a superposition \f$ \sum_n \exp(-I\cdot k\cdot n)\exp\left(-\frac{1}{2w^2}(n-n_0)^2\right) \f$. Here the antimesons (mesons starting at even sites) are taken into account).*/
  void getMesonAntiMesonOperator(MPO &mpo, int n, double k=1.0, double width=-1.0) const;
  
  /** Construct the discretized version of the momentum operator
   * \f[
   * \hat{P}=\int dx \psi^\dagger(x)i\partial\psi(x)
  \f]
  which is given by 
  \f[
   * \hat{O}_p=-ix\sum_n\left[\sigma_n^-\sigma_{n+1}^z\sigma_{n+2}^+ - \mathrm{h.c.}\right]
  \f]
  after making it dimensionless like the Hamiltonian. The constant \f$ x \f$ is not taken into account and has to be added by hand afterwards
  */
  void constructMomentumMPO(MPO& Pmpo) const;
  
};

 

#endif
