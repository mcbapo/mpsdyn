#ifndef KONDO_H
#define KONDO_H

#include "Hamiltonian.h"
#include <vector>

/** 
    Implementation of Kondo Hamiltonian with site 0 containing the
    spin 1/2 impurity, and L+1 sites to the right containing fermions
    with spin (after Jordan Wigner transformation, those sites contain
    two spin 1/2 systems each).
    \f[
    H=TODOTODO

    \f]
*/

class KondoHamiltonian:
public Hamiltonian{
protected:
  int d; // dimension of individual spin sites (2)
  int L; // chain length is L+2, with L+1 fermionic sites
  double Jpar,Jperp,t;
  double omega; // a local term for the impurity
  MPO hamil;
  int nrOps;
  mwArray Zsp,Zferm;

 public:
  /** Create a Kondo Hamiltonian with L+1 fermionic sites */
  KondoHamiltonian(int L,double J,double t);
  
  /** Create a Kondo Hamiltonian, but optionally include a local Sz
      term for the impurity (designed to be used to locate the Fermi
      sea of the bath avoiding the degeneracy wrt impurity spin when
      the coupling is zero).  */
  KondoHamiltonian(int L,double Jpar,double Jperp,double t,double omega=0.);
  
  ~KondoHamiltonian();

  const MPO& getHMPO() const {return hamil;}
  
  /** Return the Hamiltonian MPO but including a
      penalty term (coefficient penalty) to enforce the total Sz to be
      the desired Starget. */
  void getHMPOwithPenalty(MPO& mpo,double penalty,double Starget) const;
  bool hasHMPO() const {return true;}

  /** Return, in this spin language, the operators corresponding to 
      \f$\sum_{\sigma,\sigma'}f_{n\sigma}^{\dagger} \sigma_{\sigma \sigma'}^{\alpha}f_{n\sigma}\f$
      for \f$\alpha=x,\,y,\z\f$. They are site independent. */
  void getFermionicOperators(mwArray& Fx,mwArray& Fy,mwArray& Fz) const;

  void getFermionNrMPO(MPO& mpo) const;
  void getTotalSzMPO(MPO& mpo) const;

  /** Construct the MPO for the integrated even correlations: 
      \f[I_{\mathrm{even}}=\sum_{x=0,2,4,\ldots} \langle \vec{S}\cdot \vec{F}_x\rangle\f]
      where \f$\vec{F}\f$ is the vector with components the fermionic operators above. 
  */
    void getIntegratedEvenCorrelationsMPO(MPO& mpo) const;

  /** Construct the MPO corresponding to the evolution operator with
      the splitting even-odd in the Hamiltonian. Take as argument the
      complex factor in the exponential, delta.  For real time
      evolution, delta is -i*tau, for imaginary time, it should be
      -tau. */
  virtual void getExponentialMPOeven(MPO& expHe,complex_t delta) const ;
  virtual void getExponentialMPOodd(MPO& expHo,complex_t delta) const ;

 protected:
  
  /** Prepare the even-odd decomposition of the evolution operator:
      returns the two terms to be inserted in the MPO for the 
      exp(delta H12) (delta a complex number).
      TI is not assumed and the exponential of the term
      starting in pos is returned.  
      TODO: why is this method (and the following one) public??
  */
  void getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,
				 int pos=0) const;
  /** Compute the two-body term on positions (\param<pos>,\param<pos>+1):
   */
  void computeTwoBodyTerm(mwArray& result,int pos) const;

 private:
  void constructOperatorProduct(mwArray& result,const mwArray& opA,
				const mwArray& opB) const;
  void initZ();
  void initHMPO();
  //void initHMPO(double* tis=0,double* Uis=0,double* muis=0,
  //	 double* Vis=0);
  /** Common computation for recovering the even and odd parts of the exp */
  void getExponentialMPO(MPO& expHo,complex_t delta,bool even) const;

  void computeTwoBodyTermSF(mwArray& result) const;
  void computeTwoBodyTermFF(mwArray& result,bool first) const;



};


#endif // KONDO_H
