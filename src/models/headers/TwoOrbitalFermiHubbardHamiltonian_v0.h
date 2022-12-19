#ifndef TWOORBFH_H
#define TWOORBFH_H

#include "Hamiltonian.h"
#include <vector>

/** 
    Implementation of two-orbital Fermi-Hubbard Hamiltonian on a finite chain with L 
    sites containing fermions with spin, with the form
    \f{eqnarray*}{ 
    H=&-t\sum_{j,\alpha=g e,\sigma=\uparrow \downarrow}\left (c_{j\alpha\sigma}^{\dagger}c_{j+1\alpha\sigma}
+\mathrm{h.c.} \right )
+U\sum_j n_{j\alpha\uparrow }n_{j\alpha\downarrow }\\
&+\sum_{i \sigma \sigma'}\left(V-\delta_{\sigma \sigma'} V_{\mathrm{ex}}\right) n_{i e \sigma} n_{i g \sigma'} \\
&+V_{\mathrm{ex}}\sum_{i,\sigma\neq \sigma'} c_{i g \sigma}^{\dagger}
c_{i e \sigma'}^{\dagger}
c_{i g \sigma'} c_{i e \sigma}
-\sum_{i \alpha} \mu_{\alpha} n_{i\alpha}+\lambda \left(\sum_{i \alpha \sigma} n_{i\alpha\sigma}-N_{\tot}\right)^2
    \f}
    After Jordan Wigner transformation, each site contains four spin
    1/2 systems each, which in this implementation are kept separate,
    so that the spin chain has 4L chains.
    The penalty term at the end is optional. If \f$\lambda \neq 0\f$,
    it is used to enforce the desired total number of fermions in the
    ground state search (not used in the time evolution).

Notice that this version conserves Ne and Ng independently, and does not seem to agree with the papers!!!!
*/

class TwoOrbitalFermiHubbardHamiltonian:
public Hamiltonian{
protected:
  int d; // dimension of individual spin sites (2)
  int L; // chain length is L (actual, 4L)
  double tg0,tg1,te0,te1,Ug,Ue,V,Vex,mu_g0,mu_g1,mu_e0,mu_e1; // parameters in H
  double tge0,tge1; // to allow changes Ne->Ng
  MPO hamil;
  int nrOps;
  mwArray Z;
  // four different penalty terms can be chosen, after construction
  int Ntot;double penNtot;
  int Ne;double penNe;
  int Ng;double penNg;
  int Sz;double penSz;

 public:
  /** Create a two-orbital FermiHubbard Hamiltonian with L fermionic sites */
  TwoOrbitalFermiHubbardHamiltonian(int L,const std::vector<double>& t,
				    const std::vector<double>& U,double V,double Vex,
				    const std::vector<double>& mu_g,const std::vector<double>& mu_e,
				    int Ntot=0,double penalty=0.);

  ~TwoOrbitalFermiHubbardHamiltonian();

  /** The penalty terms can be set after construction (which only
      allows the total number of fermions). If this method is called,
      the penalties are ALL reset, and it supersedes the constructor
      parameters. Notice also that if both Ne and Ng are set, Ntot is
      superfluous. If one of the penalties is zero, that term is ignored.  */
  void setPenalties(int Ng,double penaltyNg,int Ne,double penaltyNe,int Sz,double penaltySz,int Ntot,double penaltyNtot);
  
  const MPO& getHMPO() const {return hamil;}

  void getHMPO(MPO& mpo) const;
  /* /\** Return the Hamiltonian MPO but including a */
  /*     penalty term (coefficient penalty) to enforce the total Sz to be */
  /*     the desired Starget. *\/ */
  /* void getHMPOwithPenalty(MPO& mpo,double penalty,double Starget) const; */
  bool hasHMPO() const {return true;}

  // Prepare and return in the argument an MPO for the total number of fermions
  void getFermionNrMPO(MPO& mpo) const; 
  // Return an MPO for the number of fermions in one orbital (default g)
  void getFermionNrOrbitalMPO(MPO& mpo,bool orbital0=true) const; 
  // Return an MPO for the spin operator Sz 
  void getTotalSzMPO(MPO& mpo) const; 
  // Idem for Sx (resp Sy) totla spin operators
  void getTotalSxMPO(MPO& mpo) const; 
  void getTotalSyMPO(MPO& mpo) const; 
  
  /** Prepare an MPO with the double occupancy of g (if orbital0==1) or e orbital */
  void getDoubleOccupancyOrbitalMPO(MPO& mpo,bool orbital0) const;

  /** Prepare a MPO for \f$\sum_n S_n^{\alpha}S_{n+1}^{\alpha}\f$, where $\alpha=x,y,z$, according to parameter type (1=x, 2=y, 3=z)*/
  void getSpinCorrelatorMPO(MPO& mpo,int type) const;

  /** Prepare a MPO for \f$\sum_n T_{n}^{z}T_{n+1}^{z}\f$ */
  void getOrbitalSpinCorrelatorMPO(MPO& mpo) const;

  /** Construct the MPO corresponding to the evolution operator with
      the splitting even-odd in the Hamiltonian. Take as argument the
      complex factor in the exponential, delta.  For real time
      evolution, delta is -i*tau, for imaginary time, it should be
      -tau. */
  virtual void getExponentialMPOeven(MPO& expHe,complex_t delta) const {std::cout<<"Not the right method here!"<<endl;exit(1);}
  virtual void getExponentialMPOodd(MPO& expHo,complex_t delta) const {std::cout<<"Not the right method here!"<<endl;exit(1);}

 protected:
  
  /* /\** Prepare the even-odd decomposition of the evolution operator: */
  /*     returns the two terms to be inserted in the MPO for the  */
  /*     exp(delta H12) (delta a complex number). */
  /*     TI is not assumed and the exponential of the term */
  /*     starting in pos is returned.   */
  /* *\/ */
  /* void getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta, */
  /* 				 int pos=0) const; */

  /* /\** Compute the two-body term on positions (\param<pos>,\param<pos>+1): */
  /*  *\/ */
  /* void computeTwoBodyTerm(mwArray& result,int pos) const; */

 private:

  void initZ();
  void initHMPO();

  void getNumberMPO(MPO& mpo,int type) const;

  /** Particular correlator in x or y direction */
  void getSpinCorrelatorMPOXY(MPO& mpo,int type) const;

  /** Common structure to create the MPO for total spin x,y operators*/
  void getTotalSxyMPO(MPO& mpo,bool isX) const; 
  
  //  /** Common computation for recovering the even and odd parts of the exp */
  //  void getExponentialMPO(MPO& expHo,complex_t delta,bool even) const;




};


#endif // TWOORBFH_H
