#ifndef SU2_H
#define SU2_H

#include "MPS.h"
#include "Hamiltonian.h"
#include "Contractor.h"
#include "LGTstring.h"

/**
   \file SU2Hamiltonian.h Definition of the class that implements a truncated SU2 Hamiltonian with staggered fermions from Erez, Ignacio and Benni
   
   \author Stefan Kühn
   \date 12/08/2014
*/


/** List of possible single body operators*/
enum SingleBodyOperator{
  sxOp,	//sigma_x
  syOp,	//sigma_y
  szOp,	//sigma_z
  spOp,	//sigma_plus
  smOp,	//sigma_minus
  LxOp,	//Lx 
  LyOp,	//Ly
  LzOp,	//Lz
  RxOp,	//Rx
  RyOp,	//Ry
  RzOp,	//Rz
  JsqOp	//Jsquared
};

/** Enum for operator components*/
enum Component{
  x_comp, y_comp, z_comp
};

/** List of states that can be prepared as initial states */
enum InitState{
  is_updown, //Strong coupling ground state
  is_downup, //Strong coupling ground state with inverted spin configuration
  is_up, //All Spins up, now flux
  is_down, //All Spins down, now flux
  is_updown_onsite, //Up down configuration on every site, no flux
  is_downup_onsite, //Down up configuration on every site, no flux
  str_teststate_1, //Up-down on every site +  finite flux in some region
  str_teststate_2, //Fully polarized + finite flux in some region
};

class SU2Hamiltonian: 
public Hamiltonian{
private:
  //Pauli Matrices
  mwArray sigma_x,sigma_y,sigma_z,sigma_plus,sigma_minus;
  //Fermionic idenitity
  mwArray id_fermi;
  //Id+s3
  mwArray idpsz;
  //Gauge Operators and adjoints
  mwArray U00,U01,U10,U11,U00ad,U01ad,U10ad,U11ad;
  //Generators for SU(2) algebras
  mwArray Lx,Ly,Lz,Rx,Ry,Rz;
  //Casimir operator of the SU(2) algebra
  mwArray Jsquared;
  //Identity for a link
  mwArray id_link;
  //Properly reshaped identites which can be used for single body MPOs
  mwArray id_link_mpo,id_fermi_mpo;    
  //Matrices for the correlation function MPO
  mwArray charge_left,charge_right,unitary_left,unitary_right;
  //Flag if the correlation MPO matrices are saved
  bool correlation_matrices_saved;
  
  /**Get the spin on site \f$ n\f$, where n is now an index ranging from \f$ 1\f$ to \f$ N_\mathrm{tot}\f$. The index given has to be a valid fermionic site, this is not checked in the routine!*/
  void getSingleBodyMPO(MPO &Spin, int index,SingleBodyOperator Op) const;
  
  /** Get the operator for \f$ G_i=\sum_nG^\dagger_i(n)G_i(n),\,\, i\in\{x,y\}\f$ as these two components are quite similar*/
  void getGaussLawMPOxy(MPO &GaussLaw, Component comp) const;
  
  /** Get the operator for \f$ G_z=\sum_nG^\dagger_z(n)G_z(n)\f$ as this one it is quite different from the other two components*/
  void getGaussLawMPOz(MPO &GaussLaw) const;
  
  /** Construct two MPOs which allow to generate a string on top of the strong coupling ground state by applying first Podd (which acts on the even sites inside the specified region) and then Peven (which acts on the odd sites inside the specified region)*/
  void getStringMPO(MPO &Podd, MPO &Peven, int startpos, int endpos) const;
  
  /** Returns an MPS which is the superposition of all strings where the left charges are centered around pos_left and the right charges are centered around pos_right and have opposite momenta such that the string should start to stretch and eventually break. If no positions are specified, then the charges will be centered approximately in the middle of the left half and the right half of the chain. */
  void getMovingString_particles_antiparticles(MPS &mps,int pos_left, int pos_right, double width_left, double width_right, double k_vec, bool factor_on) const;
  
  /** Same as before, but here only strings starting at odd sites are considered (meaning, that in one direction only the particles are moving and in the other the antiparticles */
  void getMovingString_particles(MPS &mps,int pos_left, int pos_right, double width_left, double width_right, double k_vec, bool factor_on) const;
  
  /** Sometimes I have to split a k-local operators \f$ Op=O(1)\otimes O(2) \times \dots \otimes O(k) \f$ via SVD to get it in MPO form. This routine takes an mwArray Op which is the result of a kronecker product of single site operators and splits it into an matrices suitable to construct an MPO. The vector dims has to contain the physical dimensions of the single sites in the order they appear. The reference to the vector matrices will be overwritten and contains the MPO matrices afterwards.
   * \warning The kronecker product has to be taken with the index ordering used in the code (which is different from the one from MATLAB or MAPLE), otherwise the results are wrong */
  void splitLocalOp(mwArray Op,vector<int> dims,vector<mwArray> &matrices) const;
  
  protected:
  //Number of fermions
  int N_fermions;
  //Total number of sites (3*N_fermions-1)
  int N_total;  
  //Mass
  double mu;
  vector <double> mu_site;
  //Coupling constant
  double g;
  //Hopping constant 
  double epsilon;  
  //Hamiltonian MPO
  MPO hamil;
  //Hilbert space dimension of the links
  int d_link;
  //Dimension of the fermionic Hilbert spaces
  int d_fermi;
  
   /** Generate some matrices used for single and few body MPOs to be saved in the object, therefore in the MPO the pointers to the matrices can be stored and do not consume memory*/
  void initOperators(void);
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously)*/
  void initHamiltonian(void);
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously) where the mass is now a site depended quantity*/
  void initHamiltonian_sitedependendmass(void);
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously) where the mass is now a site depended quantity and some charges are forced to sit at certain sites*/
  void initHamiltonian_sitedependendmass_charge(vector<int> charge_pos, double lambda);
  
  
  public:
  /** Creates an instance of the truncated SU(2)-Hamiltonian. The form of the Hamiltonian is given by 
   \f[ H = M\sum_n (-1)^n\psi_n^\dagger\psi_n +  \varepsilon\sum_n \left(\psi_n^\dagger U_{n}\psi_{n+1}+h.c.\right)+\frac{g^2}{2}\sum_n L_{z,n}^2,\f]
   where \f$ \psi_n \f$ is a two component color spinor containing two single component fermionic fields \f$ \phi_{1,2}\f$ where the index refers to the species of fermions.
   After a generalized Jordan Wigner transformation and we obtain the spin Hamiltonian which is implemented in this class. The formula is the following:
   \f{eqnarray*}
   H = &\varepsilon &\sum_n\left[\sigma^+_1(n)\sigma^z_2(n) U_{00,n}\sigma_{1}^-(n+1)\right.\\
      &+& i\sigma^+_1(n)\sigma^z_2(n) U_{01,n}\sigma_{1}^z(n+1)\sigma_{2}^-(n+1)\\
      &-&i \sigma^+_2(n) U_{10,n}\sigma_{1}^-(n+1)\\
      &+& \sigma^+_2(n) U_{11,n}\sigma_{1}^z(n+1)\sigma_{2}^-(n+1)\\
      &+& \left.h.c.\right]\\
      &+&\frac{M}{2}\sum_n(-1)^n\left[ (\sigma^z_{1}(n)+1)+(\sigma^z_{2}(n)+1)\right] \\
      &+&\frac{g^2}{2} \sum_n L_{z,n}^2.
   \f}
 */
  SU2Hamiltonian(int N,double epsilon,double mu,double g);
  
  /** Same as other constructor, however here it is possible to have a site depended mass. The vector with with masses has to have the length N*/
  SU2Hamiltonian(int N,double epsilon,vector<double> mu,double g);
  
  /** Same as other constructor, however here it is possible to have a site depended mass plus charges specified at the position charge_pos. The charge at position n is ensured via the penalty \f$ \lambda\left(\sum_\alpha Q_\alpha(n) Q_\alpha(n) - \frac{3}{4}\right) = \lambda\left(1+\sigma^z_{1}(n)\sigma^z_{2}(n)\right) \f$. As \f$ \sum_\alpha Q_\alpha(n) Q_\alpha(n) \f$ is a Casimir Operator, this is gauge invariant.*/
  SU2Hamiltonian(int N,double epsilon,vector<double> mu,double g,vector<int> charge_pos,double lambda);
  
  ~SU2Hamiltonian();
  
  /** Return the MPO representation of the Hamiltonian*/
  const MPO& getHMPO() const {return hamil;}
  /** Return whether the MPO representation exists*/
  bool hasHMPO() const {return true;}
  
  /** Getter for the value of the interaction \f$ \varepsilon \f$*/
  double get_epsilon() const {return this->epsilon;};
  
  /** Setter for the value of the interaction \f$ \varepsilon \f$*/
  void set_epsilon(double epsilon_value) {this->epsilon=epsilon_value;}; 
  
  /** Getter for the value of the mass \f$ \mu \f$*/
  double get_mu() const {return this->mu;};
  
  /** Setter for the value of the interaction \f$ \mu \f$*/
  void set_mu(double mu_value) {this->mu=mu_value;}; 
  
  /** Setter for site dependend values of the interaction \f$ \mu \f$*/
  void set_mu(vector<double> mu_values) {this->mu_site=mu_values;}; 
  
  /** As the Hamiltonian after changing \f$ \varepsilon \f$ might not be immediately needed, the user has to update the Hamiltonian MPO manually with this function
   \warning Only working for site independent mass */
  void updateHMPO(){this->initHamiltonian();};
  
  /** As the Hamiltonian after changing \f$ \varepsilon \f$ might not be immediately needed, the user has to update the Hamiltonian MPO manually with this function
   \warning Only working for real site depended mass, not the special case where only one mass is given */
  void updateHMPOSiteDepMass(){this->initHamiltonian_sitedependendmass();};
  
  /** This function sets the matrices responsible for the gauge interaction to identities. Therefore one can compare the full model to cases where one neglects the gauge interaction. The procedure of removing the gauge matrices cannot be undone. Therefore this function has to be used with care.*/
  void removeGaugeInteraction();  
  
  /** Generate the the approximated time evolution MPO   \f[U(t+\Delta t,t)=1-i\Delta t H\f]   with time step dt. As this approximation is not unitary, the state has to be normalized after every step. If the flag imag_time is set to true, than the MPO for imaginary time evolution instead of real time will be constructed.
  */
  void constructEvolutionMPO(MPO &Umpo, double dt, bool imag_time=false) const;
  
  /** Generate the the approximated time evolution MPO   \f[U(t+\Delta t,t)=1-i\Delta t H\f]   with time step dt. As this approximation is not unitary, the state has to be normalized after every step. This version takes into account a possible site depended mass. If the flag imag_time is set to true, than the MPO for imaginary time evolution instead of real time will be constructed.
  */
  void constructEvolutionMPOSiteDepMass(MPO &Umpo, double dt, bool imag_time=false) const;
  
  /** Generate the the approximated time evolution MPO   \f[U(t+\Delta t,t)=1-i\Delta t H\f]  with time step dt. As this approximation is not unitary, the state has to be normalized after every step. This version takes into account a possible site depended mass plus static charges. If the flag imag_time is set to true, than the MPO for imaginary time evolution instead of real time will be constructed. */  
  void constructEvolutionMPOSiteDepMassCharge(MPO &Umpo, double dt, vector<int> charge_pos, double lambda, bool imag_time=false) const;
  
  /**Get the spin on site \f$ n\f$ for species 1, where n has to be an index ranging from \f$ 1 \f$ to \f$ N_\mathrm{fermions}\f$, by default the z-component is selected.
   */
  void getSpin1MPO(MPO &Spin, int index, Component comp=z_comp) const;
  
  /**Get the spin on site \f$ n\f$ for species 2, where n has to be an index ranging from \f$ 1 \f$ to \f$ N_\mathrm{fermions}\f$, by default the z-component is selected.
   */
  void getSpin2MPO(MPO &Spin, int index, Component comp=z_comp) const;
  
  /**Get the \f$ L_i(n)\f$ , where n has to be an index ranging from \f$ 1 \f$ to \f$ N_\mathrm{fermions}\f$ and \f$ i \f$ labels the component, by default the z-component is selected.
   */
  void getLMPO(MPO &Spin, int index, Component comp=z_comp) const;
  
  /**Get the \f$ R_i(n)\f$ , where n has to be an index ranging from \f$ 1 \f$ to \f$ N_\mathrm{fermions}\f$ and \f$ i \f$ labels the component, by default the z-component is selected.
   */
  void getRMPO(MPO &Spin, int index, Component comp=z_comp) const;
  
  /**Get the \f$ J^2(n)\f$ , where n has to be an index ranging from \f$ 1 \f$ to \f$ N_\mathrm{fermions}\f$
   */
  void getJsquaredMPO(MPO &Spin, int index) const;
  
  /** Get the operator for \f$ G_i=\sum_nG^\dagger_i(n)G_i(n),\,\, i\in\{x,y,z\}\f$ which allows to monitor if Gauss' Law is fulfilled*/
  void getGaussLawMPO(MPO &GaussLaw, Component comp) const;
  
  /** Create an inital state, right now the options are
  \f{eqnarray*}
  \mathrm{is\_updown} &= |\uparrow\rangle |\uparrow\rangle \otimes |0\rangle\otimes |\downarrow\rangle|\downarrow\rangle\dots\\ 
  \mathrm{is\_downup} &= |\downarrow\rangle|\downarrow\rangle\otimes |0\rangle\otimes|\uparrow\rangle |\uparrow\rangle \dots\\
  \mathrm{is\_up} &= |\uparrow\rangle|\uparrow\rangle\otimes |0\rangle\otimes|\uparrow\rangle |\uparrow\rangle \dots\\
  \mathrm{is\_down} &= |\downarrow\rangle|\downarrow\rangle\otimes |0\rangle\otimes|\downarrow\rangle |\downarrow\rangle \dots\\
  \f} */
  void constructInitialMPS(MPS &mps, InitState state=is_updown) const;
  
  /** Get projection operators for all the even and all the odd sides which allow to project back into the Gauss Law fulfilling sector (with eigenvalue 0) 
   \warning This projection operators have varying bond dimension which goes up to 100 at certain sites, an application is costly!*/
  void constructProjectors(MPO& Podd, MPO& Peven) const;
  
  /** Get projection operators for all the even and all the odd sides which allow to project in the subspace \f$ G_\alpha(n)=0\,\forall n\neq pos \,\wedge\, n\neq pos+leng \f$ and \f$ G_\alpha(pos)=1/2 \f$, \f$ G_\alpha(pos+leng)=-1/2 \f$. The component \f$ \alpha \f$ is selected by comp. If a position smaller than zero is specified, then the automatic placing of the routines creating strings is used. */
  void constructProjectors(MPO& Podd, MPO& Peven, int pos, int leng, Component comp) const;
  
  /** Create a string on top of the strong coupling ground state given state (in principle other states are also possible but the resulting object will not be a string). If a position between 1 and the number of fermions is specified, a string with specified length is created starting from the given position. If a position smaller than zero specified, the string will be placed approximately in the middle of the chain*/
  void createString(MPS& mps, Contractor &contractor,int D, unsigned int length, int pos=-1);
  
  /** Returns an MPO which allows to create a string on top of an arbitrary state. If a position in between \f$ 1 \f$ and \f$ N_\mathrm{fermions}\f$ is specified, the resulting string will start at the desired position, if a zero is given, the string will be placed approximately in the middle of the chain. The operator is chosen such that on the strong coupling vacuum a valid string would be generated which means \f$ \Psi_nU_n\dots U_{n+l-1}\Psi^\dagger_n \f$ for strings starting at odd n and \f$ \Psi_n^\dagger U_n^\dagger\dots U_{n+l-1}^\dagger\Psi_n \f$ starting at even n. If the conjugate is set true, this is flipped around. */
  void getStringMPO2(MPO &mpo, unsigned int length, unsigned int pos=0,bool conjugate=false) const;
  
  /** Get the MPO representation of the charge operator \f$ Q_{n,\alpha}= \frac{1}{2} \Psi_n^\dagger\sigma_\alpha\Psi_n \f$. Here pos is an a site index ranging between \f$ 1 \f$ and \f$ N_\mathrm{fermions}\f$ and comp labels the color component. */
  void getChargeMPO(MPO &mpo, unsigned int pos, Component comp) const;
  
  /** Get the MPO representation of \f$\sum_n Q_{n,\alpha}^\dagger Q_{n,\alpha} \f$. Note that \f$ Q_{n,\alpha}^\dagger Q_{n,\alpha} \f$ does not seem to gauge invariant at first as the group indices are not contracted. However \f$ Q_{n,x}^\dagger Q_{n,x}=Q_{n,y}^\dagger Q_{n,y}=Q_{n,z}^\dagger Q_{n,z} \f$, therefore \f$ \sum_\alpha Q_{n,\alpha}^\dagger Q_{n,\alpha} = 3Q_{n,\alpha}^\dagger Q_{n,\alpha} \f$ and therefore \f$ Q_{n,\alpha}^\dagger Q_{n,\alpha} \f$ is nevertheless gauge invariant. */
  void getChargeSquareMPO(MPO &mpo) const;
  
  /** Get the MPO representation of \f$ Q_{n,\alpha}^\dagger Q_{n,\alpha} \f$. Note that \f$ Q_{n,\alpha}^\dagger Q_{n,\alpha} \f$ does not seem to gauge invariant at first as the group indices are not contracted. However \f$ Q_{n,x}^\dagger Q_{n,x}=Q_{n,y}^\dagger Q_{n,y}=Q_{n,z}^\dagger Q_{n,z} \f$, therefore \f$ \sum_\alpha Q_{n,\alpha}^\dagger Q_{n,\alpha} = 3Q_{n,\alpha}^\dagger Q_{n,\alpha} \f$ and therefore \f$ Q_{n,\alpha}^\dagger Q_{n,\alpha} \f$ is nevertheless gauge invariant.*/
  void getChargeSquareMPO(MPO &mpo, unsigned int pos) const;
    
  /** Get the MPO representation of \f$ \sum_n \Psi_{n}^\dagger \Psi_{n} \f$. Note that this is only the number of fermions for the non-interacting case. This quantity commutes with the Hamiltonian and therefore should be conserved. */
  void getFermionNumberMPO(MPO &mpo) const;
  
  /** Get the MPO representation of \f$ \Psi_{n}^\dagger \Psi_{n} \f$. Note that this is only the number of fermions for the non-interacting case.*/
  void getFermionNumberMPO(MPO &mpo, unsigned int pos) const;
  
  /** Get the Number of sites */
  int getSites(){return this->N_fermions;};
  
  /** Build the MPO the reduced MPO \f$ O=I_A\otimes Tr_A(S|0\rangle\langle 0|S^\dagger) \f$ such that \f$ \langle\Psi|O|\Psi\rangle \f$ gives the probability to have string in region B.*/
  void buildReducedStringMPO(MPO &mpo,const MPS &mps,int index_left_site,int index_right_site) const;
  
  /** Build the string projector \f$ P=\sum_{i_k\neq 2}|2\rangle |i_1\rangle \dots \ |i_{L-1}\rangle|2\rangle \langle 2| \langle i_1| \dots \ \langle i_{L-1}|\langle 2| \f$ which detects a string of length \f$ L \f$ starting at position pos. A string is here a stat with flux different from zero on the links in between two links with flux zero. The resulting MPO is a projector which projects all non strings state on zero. If fermions is set to a value different from zero, it also checks if the configuration inside the string is valid (meaning no single fermion on a site inside the string). If the value 1 is given, it is only checked if there is no site with a single fermion on the site, if the value 2 is given, occupation is right (particles on odd sites, antiparticles on even sites). */
  void getStringProjector(MPO &mpo,int length,int pos, unsigned int fermions=0) const;
  
  /** Function to generate an initial distribution which is a large superposition of strong coupling string states. The superposition is created in such a way, that that the left charges are distributed according to a gaussian peaked around pos_left with width width_left and the right charges likewise around pos_right with width width_right. Additionally the charges on the left end of the string are given momenta to the left and the charges at the right end of the string are given momenta to the right. The result is an initial state with a charge distribution having two peaks and the charges separating during the evolution.
   * 
  If negative values for position and/or width are specified, default values are used. If the flag antiparticles is set to false, then only strings with a particle on the left end and an antiparticle on the right hand are created. Otherwise also strings starting with an antiparticle are created. The value of k_vec controls the magnitude of the k-vector used, it is 2*Pi/(N+1)*k_vec. */
  void getMovingString(MPS &mps, bool antiparticles=false, int pos_left=-1, int pos_right=-1, double width_left=-1.0, double width_right=-1.0,double k_vec=1.0,bool factor_on=true) const;
  
  /** Function to generate an MPO which can be applied to any state and generates a superposition of string states. The superposition is created in such a way, that that the left charges are distributed according to a gaussian peaked around pos_left with width width_left and the right charges likewise around pos_right with width width_right. Additionally the charges on the left end of the string are given momenta to the left and the charges at the right end of the string are given momenta to the right. The result is an initial state with a charge distribution having two peaks and the charges separating during the evolution.
   * 
  If negative values for position and/or width are specified, default values are used. The value of k_vec controls the magnitude of the k-vector used, it is 2*Pi/N*k_vec. */
  void getMovingString(MPO &mpo, int pos_left=-1, int pos_right=-1, double width_left=-1.0, double width_right=-1.0,double k_vec=1.0,bool factor_on=true) const;
  
  /** Function to generate an MPO which can be applied to any state and generates a moving meson. The superposition is created in such a way, that the mesons are formed as a gaussian peaked around pos with width. 
  If negative values for position and/or width are specified, default values are used. The value of k_vec controls the magnitude of the k-vector used, it is 2*Pi/N*k_vec, in particular it allows to create left and right moving mesons. 
  \warning DEPRECATED!!!*/
  void getMovingMeson(MPO &mpo, int pos=-1,  double width=-1.0, double k_vec=1.0) const;
  
  /** Function to generate an MPO which creates a distribution of mesons (mesons starting at an odd site) and antimesions (mesons starting at an even site). This is achieved by building the superposition \f$ \sum_{n=1}^{N_\mathrm{fermions-1}} \exp(-I\cdot k\cdot n)\exp\left(-\frac{1}{2\sigma^2}(n-n_0)^2\right) S_n\f$, where \f$ \sigma \f$ is the given by the width and \f$ k \f$ is \f$\mathrm{k}_\mathrm{vec}\cdot \pi \f$. \f$ S_n \f$ is a string operator of length one (whose exact form depends on whether the site is even or odd). If antimesons is set to false, only mesons are taken into account, if mesons is set to false, only antimesons are taken into account. 
   */
  void getMovingAntimesonMeson(MPO &mpo,int pos=-1, double width=-1.0, double k_vec=1.0,  bool antimesons=true, bool mesons=true) const;
  
  /** Construct an MPO which makes all quarks moving to the right and all anti-quarks moving to the left with a given kvector*/
  void getMovingQuarks(MPO &mpo, double kvector) const;
  
  /** Construct an MPO which makes all spins of one species move to the left and to the right, depending if they are flipped on a quark or an antiquark site. The resulting object has bond dimension 1.
   \warning This allows in principle to create non gauge invariant distributions.*/
  void getMovingQuarks2(MPO &mpo, double kvector, int species) const;
  
  /** Create a (non gauge invariant) MPS having a superposition of all single particle states given momentum to the left and all antiparticle states given momentum to the right while the gauge links carry no flux.*/
  void getMovingParticles(MPS &mps,int pos_left=-1, int pos_right=-1, double width_left=-1.0, double width_right=-1.0, double k_vec=1.0);
  
  /** Simplified interface for void getStringMPS(MPS &mps, vector<string_t> str_config) const, where the string of the desired length is automatically placed approximately in the middle of the stated */
  void getHeavyString(MPS &mps, unsigned int length) const;
  
  /** Get a strong coupling MPS which has a string(s) on in the link sites specified by the vector str_confi(the fermionic sites are simply the usual strong coupling state). This can be used for imaginary time evolution or ground state search to bias the inital wave function in the right sector.*/
  void constructStringMPS(MPS &mps, vector<string_t> str_config) const;
  
  /** Get a a string operator without the fermionic part (so basically a sum over all combinations of \f$ U_{a_1a_2}U_{a_2a_3}U_{a_3a_4}\dots U_{a_{l-2}a_{l-1}}, \f$ where \f$ a_i\in\{0,1\} \f$). This can be used real time evolution to apply it to the interacting vaccum and evolve it in time to mimic realtime dynamics of string breaking in the presence of external charges.*/
  void getHeavyString(MPO &mpo, unsigned int length, unsigned int pos=0, bool conjugate=false) const;
  
  /** Construct an MPO to create a single charge at site n in the chain by flipping a single spin. The result is NOT GAUGE INVARIANT therefore not fulfilling Gauss Law anymore. The variable spin determines if the spin 1 (spin=1) or spin 2 (spin=2) is flipped.*/
  void createSingleCharge(MPO &mpo, int n, int spin=1) const;
  
  /** MPO to detect a string of length l starting at pos by looking at the particles (or charge square respectively). A string is said to be an object with a single particle at site pos and pos+leng with either no particles or two particles in the sites in between. */
  void getCountStringsMPO(MPO &mpo, int length, int pos) const;
  
  /** MPO to compute the charge squared correlation function \f$ \langle Q^\dagger_\alpha(\mathrm{pos})Q_\alpha(\mathrm{pos}) Q^\dagger_\alpha(\mathrm{pos+length})Q_\alpha(\mathrm{pos+length})\rangle \f$. If the flag projector_on is set to true, then configurations which are not strong coupling inside the string are projected out.
  \warning Compared to getChargeSquareMPO(MPO &mpo, unsigned int pos) and getChargeSquareMPO(MPO &mpo) the prefactor 1/8 is neglected!*/
  void getChargeCorrelatorMPO(MPO &mpo, int length, int pos, bool projector_on=false) const;
  
  /** MPO to compute the charge squared correlation function \f$ \langle Q^\dagger_\alpha(\mathrm{pos})Q_\alpha(\mathrm{pos})U_\mathrm{pos}\dots U_\mathrm{pos+length-1} Q^\dagger_\alpha(\mathrm{pos+length})Q_\alpha(\mathrm{pos+length})\rangle \f$ where \f$ U_\mathrm{n}=\exp\left(i\,\mathrm{factor}\,\pi Q^\dagger_\alpha(n)Q_\alpha(n)\right) \f$ are unitaries acting on the spins on site n to prevent the decay of correlations. The constant \f$ \mathrm{factor} \f$ can be chosen arbitraryly, setting it to zero corresponds to having identities everywhere in betwen.
   \warning Once factor is chosen it cannot be changed during the entire run! The prefactors are again neglected as in getChargeCorrelatorMPO*/
  void getChargeCorrelatorMPO2(MPO &mpo, int length, int pos, double factor=0.0);
  
  /** Get the MPO for the chiral condensate \f$ \frac{\sqrt{\varepsilon}}{N}\sum_{n=1}^N (-1)^n\phi^\dagger_n\phi_n =  \f$*/
  void getCondensateMPO(MPO &mpo, bool factor_off=false) const;
  
};

#endif //SU2_H
