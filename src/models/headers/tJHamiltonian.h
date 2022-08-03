#ifndef TJ_H
#define TJ_H

#include "Hamiltonian.h"
#include <vector>

/** 
    Implementation of the Hamiltonian for the t-J model on a finite
    chain with L sites containing fermions with spin (after Jordan
    Wigner transformation, those sites contain two spin 1/2 systems
    each).
    \f[
    H=-t\sum_{j,\sigma=\uparrow \downarrow}\left (c_{j,\sigma}^{\dagger}c_{j+1,\sigma}
+\mathrm{h.c.} \right )
+J\sum_{j} \left ({\bf S}_{j}\cdot {\bf S}_{j+1} -\frac{1}{4}n_{j}n_{j+1} \right )
    \f]
*/

class tJHamiltonian:
public Hamiltonian{
protected:
  int d; // dimension of individual spin sites (2)
  int L; // chain length is L
  double t,J;
  MPO hamil;
  int nrOps;
  mwArray Z;

 public:
  /** Create a tJ Hamiltonian with L fermionic sites */
  tJHamiltonian(int L,double t,double J);
    
  ~tJHamiltonian();

  const MPO& getHMPO() const {return hamil;}

  void getHMPO(MPO& mpo) const;
  /* /\** Return the Hamiltonian MPO but including a */
  /*     penalty term (coefficient penalty) to enforce the total Sz to be */
  /*     the desired Starget. *\/ */
  /* void getHMPOwithPenalty(MPO& mpo,double penalty,double Starget) const; */
  bool hasHMPO() const {return true;}

  /** Return, in this spin language, the operators corresponding to 
      \f$\sum_{\sigma,\sigma'}f_{n\sigma}^{\dagger} \sigma_{\sigma \sigma'}^{\alpha}f_{n\sigma}\f$
      for \f$\alpha=x,\,y,\z\f$. They are site independent. */
  void getFermionicOperators(mwArray& Fx,mwArray& Fy,mwArray& Fz) const;

  void getFermionNrMPO(MPO& mpo) const;
  void getTotalSzMPO(MPO& mpo) const;

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
  */
  void getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,
				 int pos=0) const;

  /** Compute the two-body term on positions (\param<pos>,\param<pos>+1):
   */
  void computeTwoBodyTerm(mwArray& result,int pos) const;

 private:

  void initZ();
  void initHMPO();


  /** Common computation for recovering the even and odd parts of the exp */
  void getExponentialMPO(MPO& expHo,complex_t delta,bool even) const;




};


#endif // FERMIHUBBARD_H
