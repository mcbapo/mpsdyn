#ifndef FA_H
#define FA_H

#include "Hamiltonian.h"
#include <vector>

/** 
    Implementation of a simple $H_s$ Hamiltonian used to compute statistics
    of steady state in a stochastic model.

    \f[ H(s)=-\sum_{i=1}^{L-1} \left [ e^{-s} \sqrt{c(1-c)} \left ( n^{[i-1]}  
    \sigma_x^{[i]} + \sigma_x^{[i-1]} n^{[i]}\right ) -2(1-2c)
    n^{[i-1]}  n^{[i]}  \right ] -2 c \sum_{i=1}^L n^{[i]}+c n^{[1]}    
    +c n^{[L]}  \f]

    TODO: Suppress the totally empty state, that gives zero energy,
    and we probably don't want

*/

class FAHamiltonian:
public Hamiltonian{
protected:
  int L; // chain length
  int d; // dimension (only 2 supported)
  MPO hamil;
  mwArray Z;
  double s;
  double c;
  double offset;
  double noZeroPenalty;
  bool pbc;
  
 public:
  /** Create the FA Hamiltonian with the specified parameters. To create the limit s->-infty directly, use c=0.5 and s<=-1111 */
  FAHamiltonian(int L,double s,double offset=0.,double c=.5,double noZeroPen=0.,bool pbc=false);

  ~FAHamiltonian();

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
  void initHMPOpbc();
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


#endif // FA_H
