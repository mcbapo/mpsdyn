#ifndef MULTI_FLAVOR_SCHWINGER_H
#define MULTI_FLAVOR_SCHWINGER_H

#include "MPS.h"
#include "Hamiltonian.h"

/**
   \file MultiFlavorSchwingerHamiltonian.h Definition of the class that implements the multi flavor Schwinger Model including a chemical potential term
   
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
  chargeOp, //staggered charge
  condensateOp, //chiral condensate with prefactor sqrt(x)/N
  condensate_nofactor_Op //chiral condensate without prefactor sqrt(x)/N
};

/** Enum for operator components*/
enum Component{
  x_comp, y_comp, z_comp
};

/** List of states that can be prepared as initial states */
enum InitState{
  is_strong_coupling,	//Strong coupling ground state
  is_anti_strong_coupling,	//Strong coupling ground state with flipped spins
  is_downup,		//All spins alternating down-up-down-up-...
  is_updown,		//All spins alternating up-down-up-down-...
  is_down,		//All spins down
  is_up			//All spins up
};


class MultiFlavorSchwingerHamiltonian: 
public Hamiltonian{
private:
  //Array to store all kinds of Pauli matrices
  mwArray Z; 
  //Array for an identity MPO
  mwArray id_fermi_mpo;
  
  /** Get the time evolution operator for the hopping part. If odd is set to true, then Uodd is constructed, if odd is set to false, then Ueven is constructed. If the flag imag_time is set to true, imaginary time is used.
  \todo Make it memorywise more efficient*/
  void getUMPO(MPO &mpo, double dt, bool odd, bool imag_time) const;
  
protected:
  //Physical dimension of a fermionic site
  int d_fermi;
  //Number of sites
  int N_sites;
  //Total number of sites
  int N_total;
  //Number of flavors
  int N_flavor;
  //Hopping, background field, penalty for total charge, penalty for particle number of certain flavor
  double x,l0,lambda,eta;
  //Flavor where a particle number should be imposed and the particle number
  int penalized_flavor,Npenalized_flavor;
  vector<double> mu,nu;
  //Hamiltonian MPO
  MPO hamiltonian;
  
   /** Generate some matrices used for single and few body MPOs to be saved in the object, therefore in the MPO the pointers to the matrices can be stored and do not consume memory*/
  void initOperators(void);
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously)*/
  void initHamiltonian(void);  
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously). This version puts a penalty on the flavor chosen previously to fix it to the selected particle number.*/
  void initHamiltonianPenalty(void);
  
  /**Apply a spin operator on flavor \f$ f \f$ on site \f$ n \f$, where n is now an index ranging from \f$ 0\f$ to \f$ N_\mathrm{sites}-1\f$ and f \f$ 0\f$ to \f$ N_\mathrm{flavors}-1\f$. The index given has to be a valid fermionic site, this is not checked in the routine!*/
  void getSingleBodyMPO(MPO &mpo,int site,int flavor, SingleBodyOperator Op) const;

 public:
  /** Creates an instance of the multi flavor Schwinger Hamiltonian for a chain of length N_sites, and N_flavor per site. The parameter x is the hopping parameter, mu and nu have to be vectors of length N_flavor specifying the mass and the chemical potential for each flavor respectively. If \f$ \lambda\neq 0 \f$ is given, a penalty is used to penalize states which have total charge not equal to zero. The resulting Hamiltonian MPO has bond dimension 3+2*N_flavor.
   \warning Right now the penalty pushing the model in the charge zero sub-sector is only working for <B>EVEN</B> number of sites*/
  MultiFlavorSchwingerHamiltonian(int N_sites,int N_flavor,double x,vector<double> mu,vector<double> nu,double lambda=0.0,double l0=0.0);
  
  /** Creates an instance of the multi flavor Schwinger Hamiltonian for a chain of length N_sites, and N_flavor per site. The parameter x is the hopping parameter, mu and nu have to be vectors of length N_flavor specifying the mass and the chemical potential for each flavor respectively. If \f$ \lambda\neq 0 \f$ is given, a penalty is used to penalize states which have total charge not equal to zero. The parameter eta specifies the strength of the penalty fixing the flavor penalized_flavor to the particle number Npenalized_flavor. This version should only be used if the penalty is really needed, as it yields a Hamiltonian MPO with bond dimension 4+2*N_flavor.
   \warning Right now the penalty pushing the model in the charge zero sub-sector is only working for <B>EVEN</B> number of sites*/
  MultiFlavorSchwingerHamiltonian(int N_sites,int N_flavor,double x,vector<double> mu,vector<double> nu,double lambda0,double l0, int penalized_flavor, int Npenalized_flavor, double eta);

  ~MultiFlavorSchwingerHamiltonian();

  /** Get the number of sites */
  int getN() const{return N_sites;}
  /** Get the masses */
  vector<double> getMu() const{return mu;}
  /** Get hopping parameter */
  double getX() const{return x;}
  /** Get chemical potentials */
  vector<double> getNu() const{return nu;}
  /** Return a reference to the MPO for the Hamiltonian. */
  const MPO& getHMPO() const{return hamiltonian;}

  bool hasHMPO() const {return true;}
  
  /** Apply a Pauli Matrix to a single site and a single flavor*/
  void getSpinMPO(MPO &mpo, int site, int flavor, Component comp=z_comp) const;  
  /** Get the condensate MPO \f$ \frac{\sqrt{x}}{N}\sum_{n=0}^{N-1}\sum_{f=0}^{N_f-1}\bar{\psi}_{n,f}\psi_{n,f}=\frac{\sqrt{x}}{N}\sum_{n=0}^{N-1}\sum_{f=0}^{N_f-1}(-1)^n\left(\frac{1+\sigma_z^{nN_f+f}}{2}\right) \f$ (simply a wrapper for the more general routine kept for compatibility). If factor_off is set to true the prefactor \f$ \sqrt{x}/N\f$ is set to 1. */
  void getCondensateMPO(MPO &mpo,bool factor_off=false) const;
  /** Get the (flavor resolved) condensate MPO. If a flavor between 0 and \f$ N_\mathrm{flavor} -1 \f$ is specified, the condensate MPO is given for this flavor, resulting in \f$ \frac{\sqrt{x}}{N}\sum_{n=0}^{N-1}\bar{\psi}_{n,f}\psi_{n,f}= \frac{\sqrt{x}}{N}\sum_{n=0}^{N-1}(-1)^n\left(\frac{1+\sigma_z^{nN_f+f}}{2}\right)\f$, if a negative number is given, all flavors are summed resulting in  \f$ \frac{\sqrt{x}}{N}\sum_{n=0}^{N-1}\sum_{f=0}^{N_f-1}\bar{\psi}_{n,f}\psi_{n,f} = \frac{\sqrt{x}}{N}\sum_{n=0}^{N-1}\sum_{f=0}^{N_f-1}(-1)^n\left(\frac{1+\sigma_z^{nN_f+f}}{2}\right) \f$. If factor_off is set to true the prefactor \f$ \sqrt{x}/N\f$ is set to 1. */
  void getCondensateMPO(MPO &mpo, int flavor,bool factor_off=false) const;
  /** Get the condensate MPO \f$ \frac{\sqrt{x}}{N}(-1)^n\left(\frac{1+\sigma_z^{nN_f+f}}{2}\right) \f$. If factor_off is set to true the prefactor \f$ \sqrt{x}/N\f$ is set to 1.*/
  void getCondensateMPO(MPO &mpo,int site, int flavor, bool factor_off=false) const;
  /** Get the total charge MPO given by
   * \f[ \sum_{n=0}^{N-1} \sum_{f=0}^{N_f-1} Q_{nf} = \sum_{n=0}^{N-1} \sum_{f=0}^{N_f-1} \frac{1+\sigma_z^{nN_f+f}}{2}-\frac{1}{2}\left(1-(-1)^n\right) = \sum_{n=0}^{N-1} \sum_{f=0}^{N_f-1} \frac{1}{2}\left(-1\right)^n+\frac{1}{2}\sigma_z^{nN_f+f} \f] */
  void getChargeMPO(MPO &mpo) const;
  /** Get the (flavor resolved) charge MPO. If a flavor between 0 and \f$ N_\mathrm{flavor} -1 \f$ is specified, the charge MPO is given for this flavor, resulting in
   * \f[ \sum_{n=0}^{N-1}  Q_{nf} = \sum_{n=0}^{N-1} \frac{1+\sigma_z^{nN_f+f}}{2}-\frac{1}{2}\left(1-(-1)^n\right) = \sum_{n=0}^{N-1} \sum_{f=0}^{N_f-1} \frac{1}{2}\left(-1\right)^n+\frac{1}{2}\sigma_z^{nN_f+f}. \f]
   If  a negative number is given, all flavors are summed resulting in \f$ \sum_{n=0}^{N-1} \sum_{f=0}^{N_f-1} Q_{nf} \f$.*/
  void getChargeMPO(MPO &mpo, int flavor) const;
  /** Get the flavor resolved charge MPO \f$ Q_{nf} = \frac{1}{2}\left(-1\right)^n+\frac{1}{2}\sigma_z^{nN_f+f} \f$ at a single site n.*/
  void getChargeMPO(MPO &mpo, int site, int flavor) const;
  /** Get the penalty MPO for the charge
   * \f[ \lambda\left(\sum_{n=0}^{N-1} \sum_{f=0}^{N_f-1} Q_{nf}\right)^2 \f]
   * If the flag const_off is set to true, then \f$\lambda\f$ is set to 1 independent of the value which was set originally in the constructor.*/
  void getChargePenaltyMPO(MPO &mpo, bool const_off=false) const;
  /** Get the particle number MPO*/
  void getParticleNumberMPO(MPO &mpo, int flavor) const;
  
  /** Get the time evolution operator for the hopping part of the Hamiltonian starting at the odd sites. If the flag imag_time is set to true, imaginary time is used.*/
  void getUoddMPO(MPO &mpo, double dt, bool imag_time=false) const;
  /** Get the time evolution operator for the hopping part of the Hamiltonian starting at the odd sites. If the flag imag_time is set to true, imaginary time is used.*/
  void getUevenMPO(MPO &mpo, double dt, bool imag_time=false) const;
  /** Get the time evolution operator for the mass, chemical potential and electric field part. If the flag imag_time is set to true, imaginary time is used. If a value of lmax>0 is given, all field values greater than lmax are truncated, otherwise the bond dimension will be proportional to \f$ N_\mathrm{sites}\cdot N_\mathrm{flavor} \f$.
  \todo Implement truncation*/
  void getUzMPO(MPO &mpo, double dt, bool imag_time=false,int lmax=0) const;
  /** Prepare an initial state, right now the options are
   \f{eqnarray*}
  \mathrm{is\_strong\_coupling} =& \underbrace{|\downarrow\rangle \otimes \dots \otimes|\downarrow\rangle}_{n=0,\,N_f \mathrm{ times}} \otimes \underbrace{|\uparrow\rangle \otimes \dots |\uparrow\rangle}_{n=1,\,N_f \mathrm{ times}}\otimes \underbrace{|\downarrow\rangle \otimes \dots \otimes|\downarrow\rangle}_{n=2,\,N_f \mathrm{ times}} \dots \\ 
  \mathrm{is\_anti\_strong\_coupling} =& \underbrace{|\uparrow\rangle \otimes \dots \otimes|\uparrow\rangle}_{n=0,\,N_f \mathrm{ times}} \otimes \underbrace{|\downarrow\rangle \otimes \dots |\downarrow\rangle}_{n=1,\,N_f \mathrm{ times}}\otimes \underbrace{|\uparrow\rangle \otimes \dots \otimes|\uparrow\rangle}_{n=2,\,N_f \mathrm{ times}} \dots \\
  \mathrm{is\_updown} =& |\uparrow\rangle \otimes |\downarrow\rangle \otimes |\uparrow\rangle \otimes |\downarrow\rangle \dots \\ 
  \mathrm{is\_downup} =& |\downarrow\rangle \otimes |\uparrow\rangle \otimes |\downarrow\rangle \otimes |\uparrow\rangle \dots \\ 
  \mathrm{is\_up} =& |\uparrow\rangle|\otimes |\uparrow\rangle \otimes |\uparrow\rangle \dots \\
  \mathrm{is\_down} =& |\downarrow\rangle|\otimes |\downarrow\rangle \otimes|\downarrow\rangle \dots 
  \f} */
  void constructInitialMPS(MPS &mps, InitState state=is_strong_coupling) const;
  
  /** Get a projector which projects on a state with N_f fermions in flavor f.
   \warning The bond dimension will be increasing from 1 to \f$ \lfloor\frac{N}{2}\rfloor \f$ starting from each site, thus an application is costly*/
  void getParticleNumberProjector(MPO &mpo, int flavor, int Nf) const;
  
  /** Get the MPO representation for \f$ \eta \Bigl(\sum_{n=0}^{N_\mathrm{sites}}\frac{1}{2}\left(\sigma_z^{nN_\mathrm{flavor}+\mathrm{flavor}}+1\right) - N_f\Bigr)^2 \f$. If \c factor_off is set to true then \f$ \eta=1.0 \f$ is used, otherwise the value of \f$ \eta \f$ form the Hamiltonian is taken.*/  
  void getParticlePenaltyMPO(MPO &mpo, int flavor, int Nf, bool const_off=false) const;
  
  /** Get local part of the time evolution operator (everything except the long range electric field part which is not local) for time step \c dt. If \c imag_time is set to true, the imaginary time version is constructed. If a \c site between \f$ 0 \f$ and \f$ N_\mathrm{sites} - 1\f$ is specified then the operator starts at this site, otherwise it is placed approximately in the middle of the system. The purpose of this operator is to test an extrapolation method by Philippe Corboz where he extrapolates in the truncation error rather than in the bond dimension (<a href="http://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.045116"> see reference here</a>) */
  void getLocalUMPO(MPO &mpo, double dt, bool imag_time, int site=-1) const;
  
  /** Get a two body correlation function MPO \f$ \langle O_{n_1,f_1}O_{n_2,f_2} \rangle, where  \f$ n_i\in\{0,\dots,N-1\}\f$, specifies the site and (\f$ f_i\in\{0,\dots,N_f-1\} \f$) the corresponding flavor on site \f$ n_ii \f$. If the flag condensate is set to true the condensate operator \f$ O_{n_i,f_i}=(-1)^{n_i\left(\frac{1+\sigma_z^{n_iN_f+f_i}}{2}\right) \f$ (without \f$ \sqrt{x}/N \f$ is used, otherwise simply the z-component of the spin.*/ 
  void getCorrelationMPO(MPO &mpo, int n1, int f1, int n2, int f2, bool condensate=false) const;
  
  
};


#endif // MULTI_FLAVOR_SCHWINGER_H
