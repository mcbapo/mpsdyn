#ifndef THIRRING_H
#define THIRRING_H

#include "Hamiltonian.h"
#include <vector>

/** 
    Implementation of Thirring model, which is a Heisenberg
    Hamiltonian with Z magnetic field and site-dependent coefficients,
    plus a penalty term to target a sector if total magnetization.
    \f{eqnarray*}{ 
    H=&\sum_{i=0}^{L-2} \left ( -S_x^{[i]}S_x^{[i+1]}-S_y^{[i]}S_y^{[i+1]} + \left(g^2+\tilde{g}^2[1+(-1)^i]\right) S_z^{[i]}S_z^{[i+1]} \right)\\
    &+ \sum_{i=0}^{L-1} \left( (-1)^i ma +\tilde{g}^2 +g^2(1-\frac{1}{2}[\delta_{i\ 0}+\delta_{i\ N-1}])+\mu \right) S_z^{[i]}\\
    &+\lambda \left ( \sum_{i=0}^{L-1} S_z^{[i]}-S_{\mathrm{tot}}\right)^2
    \f}
    Notice that everything is written in terms of spin (not Pauli)
    operators.
*/

class ThirringHamiltonian:
public Hamiltonian{
protected:
  int L; // chain length
  int d; // dimension (only 2 supported)
  MPO hamil;
  mwArray Z;
  double ma,g2,lambda,gt2,mu;
  int Starget;

 public:
  /** Create a Thirring Hamiltonian with given parameters, and
      dimension of the individual sites d (only d=2, i.e. S=1/2,
      supported) */
  ThirringHamiltonian(int L,double ma,double g2,double lambda,int Starget,int d=2,double gt2=0.,double mu=0.);

  ~ThirringHamiltonian();

  const MPO& getHMPO() const {return hamil;}

  bool hasHMPO() const {return true;}

    /** Construct the MPO corresponding to the evolution operator with
      the splitting even-odd in the Hamiltonian. Take as argument the
      complex factor in the exponential, delta.  For real time
      evolution, delta is -i*tau, for imaginary time, it should be
      -tau. */
  virtual void getExponentialMPOeven(MPO& expHe,complex_t delta) const ;
  virtual void getExponentialMPOodd(MPO& expHo,complex_t delta) const ;

  /** Idem for the evolution of mixed states */
  virtual void getDoubleExponentialMPOeven(MPO& expHe,complex_t delta) const ;
  virtual void getDoubleExponentialMPOodd(MPO& expHo,complex_t delta) const ;

  virtual void getExtendedExponentialMPOeven(MPO& expHe,complex_t delta) const ;
  virtual void getExtendedExponentialMPOodd(MPO& expHo,complex_t delta) const ;

 protected:
    /** Prepare the even-odd decomposition of the evolution operator:
      returns the two terms to be inserted in the MPO for the 
      exp(delta H(pos,pos+1)) (delta a complex number).
      The exponential of the term starting in pos is returned.  
  */
  void getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,
				 int pos=0) const;
  /** Compute the two-body term on positions (\param<pos>,\param<pos>+1):
      \f[
      h_{pos,pos+1}=- S_x^{pos} S_x^{pos+1} +
      - S_y^{pos} S_y^{pos+1}+
      \left(g^2+\tilde{g}^2[1+(-1)^i]\right) S_z^{pos} S_z^{pos+1}+
      \frac{h({pos})}{2} S_z^{pos} +
      \frac{h({pos+1})}{2} S_z^{pos+1}
      \f]
      (except for edges, where the \f$h(0),\ h(L-1)\f$ terms are not divided by two).
*/
  void computeTwoBodyTerm(mwArray& result,int pos) const;

 private:
  void initZ();
  void initHMPO();

    /** Common computation for recovering the even and odd parts of the exp */
  void getExponentialMPO(MPO& expHo,complex_t delta,bool even) const;
  void getDoubleExponentialMPO(MPO& expHo,complex_t delta,bool even) const;
  void getExtendedExponentialMPO(MPO& expHo,complex_t delta,bool even) const;

};


#endif // HEISENBERG_H
