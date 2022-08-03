#ifndef HEISENBERGANC_H
#define HEISENBERGANC_H

#include "HeisenbergHamiltonian.h"
#include <vector>

/** 
    Implementation of a chain with two types of spins, \f$\sigma_n,\,
    \tau_n\f$, corresponding to two Heisenberg chains coupled by a
    \f$\sigma_z^{[n]} \tau_z^{[n]}\f$ term. Different to
    HeisenbergHamiltonian, here all the coefficients are
    translationally invariant. The second type of spins, when traced
    out, realizes a random magnetic field for the first chain.
    \f[
    H=\sum_{i=1}^{L-1} \frac{1}{4}\left ( \vec{\sigma}^{[i]}\vec{\sigma}^{[i+1]} +
    J \vec{\tau}^{[i]}\vec{\tau}^{[i+1]} \right ) +    
    \sum_{i=1}^{L} B \frac{1}{4} \sigma_z^{[i]}\tau_z^{[i]} 
    \f]
    The \f$\frac{1}{4}\f$ are due to the fact that the implementation 
    is defined in terms of spin matrices, not sigmas, as in the parent 
    HeisenbergHamiltonian.
    The implementation is done in terms of a single chain, in which 
    \f$\sigma\f$ and \f$\tau\f$ sites alternate.
*/

class HeisenbergAncillaryNoiseHamiltonian:
public HeisenbergHamiltonian{
	double B; // I will use them instead of Jxi in parent
	double J1x,J1y,J1z; // for the first chain (default 1)
	double J2x,J2y,J2z; // for the second one in case they are different
    double offset;
	bool TI;
	std::vector<double> Bs;
 public:
  /** Create a Hamiltonian of this type with length L and coefficients
      1 for the first chain, J for the second and B for the Ising term, TI. */
  HeisenbergAncillaryNoiseHamiltonian(int L,double J,double B);
  HeisenbergAncillaryNoiseHamiltonian(int L,double J1x,double J1y,double J1z,
				      double J2x,double J2y,double J2z,double B,double offset);
  /** Create a Hamiltonian of this type with length L and coefficients
      J, and non-TI Bs */
  HeisenbergAncillaryNoiseHamiltonian(int L,double J,const std::vector<double>& B);

  ~HeisenbergAncillaryNoiseHamiltonian();

  /** Not yet supported */
  const MPO& getHMPO() const {return hamil;}

  bool hasHMPO() const {return true;} 

  /** Construct the MPOs corresponding to the evolution operator.  In
      particular, we use the even-odd splitting for the sigma and tau
      Heisenberg terms independently, and require an additional
      exponential for the mixed term. Each of the methods prepares 
      a particular term, with the given complex factor in the 
      exponential, delta. For real time evolution, delta is
      -i*tau, for imaginary time, it should be -tau. 
      The even/odd terms require a third (boolean) argument, sigma, 
      to determine whether the term for \f$\sigma\f$ (if true, default) 
      or for \f$tau\f$ (if false) is computed.
*/
  void getExponentialMPOeven(MPO& expHe,complex_t delta,bool sigma=true) const;
  void getExponentialMPOodd(MPO& expHo,complex_t delta,bool sigma=true) const;
  void getExponentialMPOmixed(MPO& expHb,complex_t delta) const;
/** To make things simpler, I also provide here the double versions
 * of these exponential operators, that are used on the density
 * matrix MPS. Then, I can keep just the needed new Operators, and
 * proper links to the rest. 
*/
  void getDoubleExponentialMPOeven(MPO& expHe,complex_t delta,bool sigma=true) const;
  void getDoubleExponentialMPOodd(MPO& expHo,complex_t delta,bool sigma=true) const;
  void getDoubleExponentialMPOmixed(MPO& expHb,complex_t delta) const;


 private:
  void initHMPO();
  void initHMPOnonTI();
  //void initHMPO(double* tis=0,double* Uis=0,double* muis=0,
  //	 double* Vis=0);
  /** Common computation for recovering the even and odd parts of the exp 
// TODO!!!! Private function that prepares the operators, because
// placing them is the same for Double or normal!
*/
  void getExponentialMPO(MPO& expHeo,complex_t delta,bool even,bool sigma=true) const;
  void getDoubleExponentialMPO(MPO& expHeo,complex_t delta,bool even,bool sigma=true) const;

  void setExponentialOperators(MPO& expH,Operator* Opl,Operator* Opr,
			       Operator* idPhys,Operator* idBoth,
			       bool even,bool sigma) const;
  // whether the exponential of a certain term contains any
  // non-identity term
  bool containsOperators(bool even,bool sigma) const;
};

#endif // HEISENBERGANC_H
