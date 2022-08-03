#ifndef SU2_NO_LINKS_H
#define SU2_NO_LINKS_H

#include "MPS.h"
#include "Hamiltonian.h"
#include "Contractor.h"
#include "LGTstring.h"

/**
   \file SU2HamiltonianNoLinks.h Definition of the class that implements a truncated SU2 Hamiltonian following the reduced gauge invariant basis presented in the paper "Hamer, C. SU(2) Yang-Mills theory in (1 + 1) dimensions: A finite-lattice approach Nuclear Physics B , 1982, 195, 503 - 521". Contrary to the version in ReducedSU2Hamiltonian.h, this one does not explicitly contain the gauge links, but instead the fermionic sites have one degree of freedom added to encode the gauge field.
   
   \author Stefan Kühn
   \date 09/10/2015
*/



class SU2HamiltonianNoLinks: 
public Hamiltonian{
protected:
  //Fermionic operators and adjoints
  mwArray O1,O2,O3,O4,O1dagger,O2dagger,O3dagger,O4dagger;
  //Fermionic occupation number operator
  mwArray Noc;
  //Fermionic identity
  mwArray id_fermi;
  //Gauge Operators and adjoints
  //mwArray Lpjdep,Lmjdep,Lpjdepdagger,Lmjdepdagger,Lm,Lp,Lmdagger,Lpdagger,Lvalop;
  //Electric energy
  mwArray Jsquared;
  //Properly reshaped identities which can be used for single body MPOs
  mwArray id_fermi_mpo;      
  
  /**Get the spin on site \f$ n\f$, where n is now an index ranging from \f$ 1\f$ to \f$ N_\mathrm{tot}\f$. The index given has to be a valid fermionic site, this is not checked in the routine!*/
  //void getSingleBodyMPO(MPO &Spin, int index,SingleBodyOperator Op) const;
  
  /** Get the time evolution operator for the hopping part. If odd is set to true, then Uodd is constructed, if odd is set to false, then Ueven is constructed. If the flag imag_time is set to true, imaginary time is used.*/
  //void getUMPO(MPO &mpo, double dt, bool odd, bool imag_time) const;
  
  //Number of fermions
  int N_fermions;
  //Mass
  double mu;
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
  //Penalty strength
  double lambda;
  //Penalty strength for the penalty term to target the subspace of total charge equal to zero
  double sz_penalty;
  
   /** Generate some matrices used for single and few body MPOs to be saved in the object, therefore in the MPO the pointers to the matrices can be stored and do not consume memory*/
  void initOperators(void);
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously)*/
  virtual void initHamiltonian(void);
  
  /** Get the MPO representation for the J or J_squared (both are structure-wise equal, only the values filled in the last matrix differ*/
  void get_J_Jsquared_MPO(MPO &JsquaredMPO, int index, bool squared) const;
  
  /**Get the MPO representation of a single body operator acting on site \f$ n\f$, where n is now an index ranging from \f$ 1\f$ to \f$ N_\mathrm{fermions}\f$. The index given has to be a valid fermionic site, this is not checked in the routine!*/
  int getSingleBodyMPO(MPO &mpo, int index, const mwArray &Op) const;
  
  /** Get the coefficients for the transitions. The flag encodes which type of coefficient I get, a value of 0 means \f$ (-1)^{l2-l1-1/2}\sqrt{\frac{2l2+1}{2l1+1}} \f$, a value of 1 means +1 and a value of 2 means -1. Furthermore it is checked whether the transition is possible with the allowed range for the flux values. If the transition is not allowed, a zero is returned.*/
  double get_coefficient(int l1, int l2, int flag) const;
  
  /** Empty default constructor for daughter class.*/
  SU2HamiltonianNoLinks();
  
  
public:
  /** Creates an instance of the truncated SU(2)-Hamiltonian in the reduced basis. In this formulation the gauge field is not explicitly implemented. Using the fact that the field value on site \f$ n \f$ can be reconstructed from its predecessors value and the fermionic content of the site
   \f{eqnarray*}
   l_n &= \pm 1/2 \quad\quad  &\mbox{ if a single fermion is on site } n \\
   l_n &= l_{n-1}  \quad\quad &\mbox{ else, } 
   \f}
   we encode the four possible cases by \f$ |0\rangle, |1\rangle, |2\rangle, |3\rangle \f$, where \f$ |0\rangle \f$ (\f$|3\rangle \f$) means that the site is occupied by no fermion (two fermions), and \f$ |1\rangle \f$ (\f$|2\rangle \f$) mean that the site is occupied by a single fermion and the field to the right is going down (up) by one quanta. Thus, by knowing \f$ l_0 \f$ one can reconstruct the field values form the fermionic configuration. For all the implementation here \f$ l_0=0 \f$ is used.
   
   The bond dimension of the Hamiltonian MPO is hence growing with the dimension of the links \f$ d_\mathrm{link}\f$ and is given by \f$ D=17+d_\mathrm{link} \f$.
   */
  SU2HamiltonianNoLinks(int N, double epsilon, double mu, double g, int dlink, double lambda=0.0, double sz_penalty=0.0);
  
  ~SU2HamiltonianNoLinks();
  
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
  
  /** As the Hamiltonian after changing a parameter might not be immediately needed, the user has to update the Hamiltonian MPO manually with this function */
  void updateHMPO(){this->initHamiltonian();};  
  
  /** Construct an arbitrary (gauge invariant) product state as MPS, the vector sites specifies the fermionic  states \f$ |0\rangle, |1\rangle, |2\rangle, |3\rangle \f$, via the numbers 0, 1, 2 and 3. The resulting state will be normalized and right gauged. If an empty vector is given, the strong coupling ground state is generated. If within the given \f$ d_\mathrm{link}\f$ the state is not possible, a warning is given.*/
  void constructProductStateMPS(MPS &mps,vector<int> sites) const;
  
  /** Get the MPO representation for the electric energy operator*/
  void getJsquaredMPO(MPO &JsquaredMPO, int index) const;
  
  /** Get the MPO representation for the operator returning the flux value on the link*/
  void getJMPO(MPO &JMPO, int index) const;
  
  /** Get the MPO representation for the fermionic occupation number operator for site n*/
  void getNocMPO(MPO &mpo, int index) const;
  
  /** Get the MPO representation for the total fermionic occupation number operator \f$ \sum_n N_n \f$. */
  void getNocMPO(MPO &mpo) const;
  
  /** Get the charge square operator, which is in the occupation number basis essentially given by \f[ Q^2_n = \frac{1}{4} \left(|1\rangle\langle 1| + |2\rangle\langle 2|\right).\f]*/
  void getChargeSquareMPO(MPO &mpo, int pos) const;
  
  /** Get the total charge square operator, which is in the occupation number basis given by \f[ \sum_n Q^2_n .\f]*/
  void getChargeSquareMPO(MPO &mpo) const;
  
  /** Given a state which would imply negative flux or exceed the maximum electric flux allowed by the given link dimension this MPO has a finite expectation value. If no strength different from zero is given, than the penalty strength of the Hamiltonian is taken. This version does not check whether \f$ l_n=0 \f$ which implies total U(1) charge equal to zero.*/
  void getPenaltyMPO(MPO &mpo, double strength=0.0) const;
  
  /** Given a state which would imply \f$ l_i<0 \f$ or \f$ l_i\geq d_\mathrm{link}\f$ for any \f$ i=0,...,N \f$ this MPO has a finite expectation value. If no strength different from zero is given, than the penalty strength of the Hamiltonian is taken. If the flag check_ln is set to true, also states with \f$ l_n!=0 \f$ are penalized.*/
  void getPenaltyMPO2(MPO &mpo, bool check_ln=false,  double strength=0.0) const;
  
  /** Something like the U(1)-charge meaning \f$ \sum_n N_{f,n} - \left(1-(-1)^n\right) \f$, where \f$ N_{f,n} \f$ is the operator for the particle number on site n. This essentially measures of there is the same number of particles and antiparticles.*/
  void getParticleBalanceMPO(MPO &mpo) const;
  
  /** Get particle antiparticle penalty MPO which is given by \f$ P_Q= \bigl(\sum_n N_{f,n} - \left(1-(-1)^n\right) \bigr)^2 \f$ */
  void getParticleBalancePenaltyMPO(MPO &mpo) const;
  
  /** Print nonzero components of tensor*/
  void print_tensor(const mwArray &tensor,bool fix_physical=true) const;
  
  /** An attempt to construct the lattice version of the non gauge invariant operator \f$ P=-i\int dx \bar{\psi}(x)\partial_x\psi(x) \f$. This version constructs the operator \f$\sum_n \left(\phi^\dagger_n\phi_{n+2} + \mbox{h.c.}\right)\f$ in the regular basis with links. This is achieved by combining it with an MPO which first unfolds the basis with no links in the basis with links and then traces out the gauge field. The resulting MPO thus gives the expectation value of the previous operator in a state where the gauge field has been traced out.*/
  void getMomentumMPO(MPO &mpo) const;
  
  /** An attempt to construct the lattice version of the non gauge invariant operator \f$ P=-i\int dx \bar{\psi}(x)\partial_x\psi(x) \f$. This version also constructs the operator \f$\sum_n \left(\phi^\dagger_n\phi_{n+2} + \mbox{h.c.}\right)\f$ in the regular basis with links. This is achieved by combining it with an MPO which first maps the fermionic content to this basis and simply ignores the gauge field.*/
  void getMomentumMPO2(MPO &mpo) const;
  
  /** Attempt to construct \f$ P=\sum_n\left(\phi^\dagger_n\phi_{n+2} + \mbox{h.c.}\right) \f$ directly in the gauge invariant basis. The resulting MPO has bond dimension 22, and is thus costly.*/
  void getMomentumMPO3(MPO &mpo) const;
  
  /** Do a cyclic shift to the right by one site. The bond dimension of this MPO is \f$ d^2 \f$ and therefore an application is costly!*/
  void getCyclicShiftMPO(MPO &mpo) const;
  
  /** Do a cyclic shift to the right by two sites, this is more efficient than joining the two cyclic shifts by one site to a single operator
   \warning Not yet working*/
  void getCyclicShiftSquareMPO(MPO &mpo) const;
  
  /** Do a (non-cyclic) shift to the right, the first site will then be filled with a zero. If flip_particles is set to true, then particles and antiparticles are flipped on each site except the first one*/
  void getShiftMPO(MPO &mpo, bool flip_particles=false) const;
  
  /** Construct an MPO which should represent the lattice version of charge conjugation. This MPO first does a cyclic shift (as in getCyclicShiftMPO) and subsequently transforms particles and antiparticles (as in getParticleTrafoMPO). If the flag project is set to true, then states ending with \f$ |1\rangle \f$ are projected out before shifting and transforming the particles. Essentially the same could be done be constructing all three MPOs separately and and then joining them, however this version is more comfortable.*/
  void getChargeConjugationMPO(MPO &mpo, bool project=false) const;
  
  /* Attempt to firste get rid of some sites at the edges and then use the bulk to compute the charge conjugation expecation value. Offset speciefies how much sites are traced away at each end, thus in total 2*offset sites are taken away from the original MPS.*/
  complex_t ChargeConjugationValue(const MPS &mps, Contractor &contractor, int offset=5) const;
  
  /** Transform "particles" into "antiparticles" which is done by the operator \f$ \otimes_{n=1}^N \left( |0\rangle \langle 3| + |3\rangle \langle 0| + |1\rangle \langle 1| + |2\rangle \langle 2|\right) \f$  */
  void getParticleTrafoMPO(MPO &mpo) const;
  
  /** Get MPO representation of projector which projects out all components of a state ending with \f$ |1\rangle \f$ at the end, as these have no definite behavior under charge conjugation*/
  void getProjectorMPO(MPO &mpo) const;
  
  //For debugging only
  void getNocMPO_other_basis(MPO &mpo, int n) const;
  
private:
  void getChargeConjugationMPO(MPO &mpo, bool project, int N) const;
};

#endif //SU2_NO_LINKS_H
