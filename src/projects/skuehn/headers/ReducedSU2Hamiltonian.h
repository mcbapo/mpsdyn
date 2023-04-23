#ifndef SU2_REDUCEDBASIS_H
#define SU2_REDUCEDBASIS_H

#include "MPS.h"
#include "Hamiltonian.h"
#include "Contractor.h"
#include "LGTstring.h"

/**
   \file ReducedSU2Hamiltonian.h Definition of the class that implements a truncated SU2 Hamiltonian following the reduced gauge invariant basis presented in the paper "Hamer, C. SU(2) Yang-Mills theory in (1 + 1) dimensions: A finite-lattice approach Nuclear Physics B , 1982, 195, 503 - 521"
   
   \author Stefan Kühn
   \date 15/08/2015
*/


/** List of possible single body operators*/
enum SingleBodyOperator{
  NfOp,		//Fermionic occupation number
  JsqOp,	//Jsquared = diag(l(l+1))
  LOp		//Loperator = diag(l)
};



class ReducedSU2Hamiltonian: 
public Hamiltonian{
private:
  //Fermionic operators and adjoints
  mwArray O1,O2,O1dagger,O2dagger;
  //Fermionic occupation number operator
  mwArray Noc;
  //Fermionic identity
  mwArray id_fermi;
  //Gauge Operators and adjoints
  mwArray Lpjdep,Lmjdep,Lpjdepdagger,Lmjdepdagger,Lm,Lp,Lmdagger,Lpdagger,Lvalop;
  //Electric energy
  mwArray Jsquared;
  //Identity for a link
  mwArray id_link;
  //Properly reshaped identities which can be used for single body MPOs
  mwArray id_link_mpo,id_fermi_mpo;      
  
  /**Get the spin on site \f$ n\f$, where n is now an index ranging from \f$ 1\f$ to \f$ N_\mathrm{tot}\f$. The index given has to be a valid fermionic site, this is not checked in the routine!*/
  void getSingleBodyMPO(MPO &Spin, int index,SingleBodyOperator Op) const;
  
  /** Get the time evolution operator for the hopping part. If odd is set to true, then Uodd is constructed, if odd is set to false, then Ueven is constructed. If the flag imag_time is set to true, imaginary time is used.*/
  void getUMPO(MPO &mpo, double dt, bool odd, bool imag_time) const;
  
  protected:
  //Number of fermions
  int N_fermions;
  //Total number of sites (2*N_fermions-1)
  int N_total; 
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
  
   /** Generate some matrices used for single and few body MPOs to be saved in the object, therefore in the MPO the pointers to the matrices can be stored and do not consume memory*/
  void initOperators(void);
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously)*/
  void initHamiltonian(void);
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously) with a penalty term for non gauge invariant states*/
  void initHamiltonianPenalty(void);
  
  
  public:
  /** Creates an instance of the truncated SU(2)-Hamiltonian in the reduced basis.*/
  ReducedSU2Hamiltonian(int N,double epsilon,double mu,double g, int Nmax);
  
  /** Creates an instance of the truncated SU(2)-Hamiltonian in the reduced basis. This version allows to include a penalty with strength \f$ \lambda \f$. 
   \warning The bond dimension of the version with penalty is significantly higher, it should therefore be only used if a finite penalty is desired and not with \f$ \lambda = 0 \f$*/
  ReducedSU2Hamiltonian(int N,double epsilon,double mu,double g, int Nmax, double lambda);
  
  ~ReducedSU2Hamiltonian();
  
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
  
  /** This function sets the matrices responsible for the gauge interaction to identities. Therefore one can compare the full model to cases where one neglects the gauge interaction. The procedure of removing the gauge matrices cannot be undone. Therefore this function has to be used with care.*/
  void removeGaugeInteraction();  
  
  /** Get the MPO representation for the electric energy operator*/
  void getJsquaredMPO(MPO &Spin, int index) const;
  
  /** Get the flux value on a link*/
  void getLMPO(MPO &mpo, int index) const;
  
  /** Get the MPO representation for the fermionic occupation number operator*/
  void getNocMPO(MPO &Spin, int index) const;
  
  /** Get the part starting on odd sites of an odd-even-splitting for the time evolution operator. If the flag imag_time is set to true, imaginary time is used.*/
  void getUoddMPO(MPO &mpo, double dt, bool imag_time=false) const;
  
  /** Get the part starting on even sites of an odd-even-splitting for the time evolution operator. If the flag imag_time is set to true, imaginary time is used.*/  
  void getUevenMPO(MPO &mpo, double dt, bool imag_time=false) const;
  
  /** Get the first order Taylor approximation for the time evolution operator as an MPO. Compared to the odd-even splitting this should preserve the gauge symmetry. If the flag imag_time is set to true, imaginary time is used.*/  
  void getUMPO(MPO &mpo, double dt, bool imag_time=false) const;
  
  /** Construct the strong coupling ground state, the result will be right gauged and normalized.*/
  void constructInitialMPS(MPS &mps) const;
  
  /** Construct an arbitrary product state as MPS, the vector sites specifies the fermionic occupation numbers as well as the states of the links on the sites in the form \f$ |n_1\rangle|2l_1\rangle|n_2\rangle|2l_2\rangle ...\f$, where \f$ l_i\in\{0,\dots,(d_\mathrm{link}-1)/2\}\f$ . The resulting state will be normalized and right gauged.*/
  void constructProductStateMPS(MPS &mps,vector<int> sites) const;
  
  /** Get a strong coupling MPS which has a string on in the link sites (the fermionic sites are simply the usual strong coupling state). This can be used for imaginary time evolution or ground state search to bias the initial wave function in the right sector.*/
  void getHeavyString(MPS &mps, unsigned int length, unsigned int flux_quanta=1, unsigned int pos=0) const;
  
  /** Get the charge square operator, which is in the occupation number basis essentially given by
   *\f[ Q^2_\alpha = \frac{1}{4} |1\rangle\langle 1|.\f]*/
  void getChargeSquareMPO(MPO &mpo, int pos) const;
  
  /**Chiral condensate MPO \f$ \frac{\sqrt{\varepsilon}}{N}\sum_{n=1}^N (-1)^n\phi^\dagger_n\phi_n \f$ */
  void getCondensateMPO(MPO &mpo, bool factor_off=false) const;
  
  /**Check (locally on position n) if the Gauss Law is fulfilled via computing an MPO representation for the operator \f$left(L_n-L_{n-1}\right)^2-Q_n^2 \f$. The operator \f$ L_n \f$ is the one giving the value \f$ l \f$ for the flux on link \f$ n \f$ (and not \f$ l(l+1) \f$ contrary to the operator for the electric energy).*/
  void getGaussLawMPO(MPO &mpo, int pos) const;
  
  /**Check globally if the Gauss Law is fulfilled. This is done with the positive semidefinite operator \f$ \lambda\left(\left(L_n^2-L_{n-1}^2\right)-Q_n^2\right)^2 \f$. If the strength is set to 0, then the value specified in the Hamiltonian for \f$ \lambda\f$ is taken. If a finite value is given, this one is used for \f$ \lambda \f$.*/
  void getGaussLawPenaltyMPO(MPO &mpo, double strength=0.0) const;
  
  /**Get a projector on the Gauss Law fulfilling subspace which is given by states \f$ |l_1\rangle \otimes|n\rangle\otimes |l_2\rangle \f$ where \f$ l_2=l_1\pm \frac{1}{2} \f$ for \f$ n=1 \f$ and \f$ l_2=l_1 \f$ for \f$ n=0,2 \f$. 
   \warning The bond dimension of this MPO is locally up to \f$ d_\mathrm{link}^2 \f$, thus an application of the projectors is costly!*/
  void getProjectorMPO(MPO &Podd, MPO &Peven) const;
  
  /**Get a projector on the Gauss Law fulfilling subspace which is given by states \f$ |l_1\rangle \otimes|n\rangle\otimes |l_2\rangle \f$ where \f$ l_2=l_1\pm \frac{1}{2} \f$ for \f$ n=1 \f$ and \f$ l_2=l_1 \f$ for \f$ n=0,2 \f$. 
   \warning This is an experiment if I can build a projector with lower bond dimension \f$ d_\mathrm{link} \f$, thus an application should be way less costly then the previous version.*/
  void getProjectorMPO2(MPO &Podd, MPO &Peven) const;
  
  /**Get the sum of local projectors on the Gauss Law fulfilling subspace which is given by states \f$ |l_1\rangle \otimes|n\rangle\otimes |l_2\rangle \f$ where \f$ l_2=l_1\pm \frac{1}{2} \f$ for \f$ n=1 \f$ and \f$ l_2=l_1 \f$ for \f$ n=0,2 \f$. Thus one can check, if the Gauss Law is fulfilled during evolution. The bond dimension of the operator is given by \f$ d_\mathrm{link}+2 \f$ */
  void getProjectorSumMPO(MPO &Pgauss) const;
  
  /**Get a projector on the Gauss Law fulfilling subspace of strings between heavy external charges. This space is given by by states \f$ |l_1\rangle \otimes|n\rangle\otimes |l_2\rangle \f$ where \f$ l_2=l_1\pm \frac{1}{2} \f$ for \f$ n=1 \f$ and \f$ l_2=l_1 \f$ for \f$ n=0,2 \f$ except for the sites pos and pos+length, where there is an offset by quanta allowed. 
   \todo Could be more efficient with bond dimension \f$ d_\mathrm{link} \f$!*/
  void getProjectorMPO(MPO &Podd, MPO &Peven, int length, int flux_quanta=1, int pos=0) const;
  
  /** Simplified interface for void getProjectorMPO2(MPO &Podd, MPO &Peven, string_t str, int flux_quanta=1) const, automatically places the string approximately in the middle of the chain. */
  void getProjectorMPO2(MPO &Podd, MPO &Peven, int length, int flux_quanta=1) const;  
  
  /**Get a projector on the Gauss Law fulfilling subspace of strings between heavy external charges. This space is given by by states \f$ |l_1\rangle \otimes|n\rangle\otimes |l_2\rangle \f$ where \f$ l_2=l_1\pm \frac{1}{2} \f$ for \f$ n=1 \f$ and \f$ l_2=l_1 \f$ for \f$ n=0,2 \f$ except for the sites pos and pos+length, where there is an offset by quanta allowed. 
   \todo Could be more efficient with bond dimension \f$ d_\mathrm{link} \f$!*/
  void getProjectorMPO2(MPO &Podd, MPO &Peven,vector<string_t> str_config, int flux_quanta=1) const;
  
  /** Generate a state containing all strong coupling strings specified in the vector str_config*/
  void constructStringMPS(MPS &mps,vector<string_t> str_config,int quanta=1) const;
};

#endif //SU2_REDUCEDBASIS_H
