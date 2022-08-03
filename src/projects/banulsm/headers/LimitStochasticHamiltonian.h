#ifndef LIMSTOCHASTIC_H
#define LIMSTOCHASTIC_H

#include "Hamiltonian.h"
#include <vector>

/** Limit for \f$s\to\infty\f$ of the <StochasticHamiltonian>, which
    reduces to
    \f[
    H(s)=-\sum_{i=1}^{L-1} n^{[i]} \sigma_x^{[i+1]}
    \f]

    The Hamiltonian commutes with the number operator on the leftmost
    site.  and we can study the sectors of occupation 0 and 1
    independently. For the occupation 0, the result is the same
    Hamiltonian acting on N-1, and for occupation 1 (default) need to
    add the \f$\sigma_x^{[1]}\f$ operator 
*/

class LimitStochasticHamiltonian:
public Hamiltonian{
protected:
  int L; // chain length
  int d; // dimension (only 2 supported)
  MPO hamil;
  mwArray Z;
  double offset;
  bool sectorOne;
  
 public:
  LimitStochasticHamiltonian(int L,double offset=0.,bool sectorOne=true);

  ~LimitStochasticHamiltonian();

  const MPO& getHMPO() const {return hamil;}

  bool hasHMPO() const {return true;}

  /** Construct the MPO corresponding to the evolution operator with
      the splitting even-odd in the Hamiltonian. Take as argument the
      complex factor in the exponential, delta.  For real time
      evolution, delta is -i*tau, for imaginary time, it should be
      -tau. */
  virtual void getExponentialMPOeven(MPO& expHe,complex_t delta) const ;
  virtual void getExponentialMPOodd(MPO& expHo,complex_t delta) const ;

  /** Construct an MPO for the commutator with an operator, and return in in the argument.  */
  void getCommutatorMPO(MPO& mpo);

 protected:

  /** Prepare the even-odd decomposition of the evolution operator:
      returns the two terms to be inserted in the MPO for the 
      exp(delta H12) (delta a complex number).
      TI is not assumed and the exponential of the term
      starting in pos is returned.  
  */
  void getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,
				 int pos=0) const;
  void getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,
				 const mwArray& H12,complex_t delta) const;
  /** Compute the two-body term on positions (\param<pos>,\param<pos>+1):
      \f[
      h_{pos,pos+1}=n^{pos}(e^{-s} \sigma_x^{pos+1}-1)
      \f]
  */
  void computeTwoBodyTerm(mwArray& result,int pos) const;

 private:
  void initZ();
  void initHMPO();
  //void initHMPO(double* tis=0,double* Uis=0,double* muis=0,
  //	 double* Vis=0);
  /** Common computation for recovering the even and odd parts of the exp */
  void getExponentialMPO(MPO& expHo,complex_t delta,bool even) const;
  //void getDoubleExponentialMPO(MPO& expHo,complex_t delta,bool even) const;
  //void getExtendedExponentialMPO(MPO& expHo,complex_t delta,bool even) const;

  // This is copied from IsingHamiltonian, but it would be better to
  // have it somewhere else (common single place)
  /** 
      Auxiliary function for the construction of the commutator
      MPO. If this should be more general (not only for this
      IsingHamiltonian case), these functions might need to be
      somewhere else. Specially the double
   */
  void initZdoub(mwArray& Z);
  void constructOperatorProduct(mwArray& result,const mwArray& opA,
				const mwArray& opB) const;



};


#endif // STOCHASTIC_H
