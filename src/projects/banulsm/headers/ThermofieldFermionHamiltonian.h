#ifndef THERMOFIELDFERMIONHAMIL_H
#define THERMOFIELDFERMIONHAMIL_H

#include <vector>
#include "Hamiltonian.h"
#include "ThermofieldHamiltonian.h"

/** 
    A Hamiltonian containing only nearest-neighbor terms, for the
    particular case of the thermofield simulations with fermionic
    bath. This means we have two fermionic modes in the center of the
    chain, and L fermionic modes coupled to it on each side. For the
    MPS implementation, a Jordan Wigner transformation is done to 
    write the Hamiltonian in spin language.

    We assume couplings of the form \f[b_{2,0}\left ( L^{\dagger}
    b_0+\mathrm{h.c.}\right)\f] for the coupling between the spin and
    the first bosonic mode on the right, with
    \f[L=d_{\uparrow}+d_{\downarrow}\f] (being \f$d_{\sigma}\f$ the
    fermionic modes of the system), and \f[b_{1,0}\left (
    L^{\dagger}c_0^{\dagger} + \mathrm{h.c.}  \right)\f] for the first
    mode on the right, and \f[b_{1 i}
    \left(b_{i+1}^{\dagger}b_i+b_{i}^{\dagger}b_{i+1}\right)\f] \f[-
    b_{2 i}
    \left(c_{i+1}^{\dagger}c_i+c_{i}^{\dagger}c_{i+1}\right)\f] for
    the right bosonic chain.  The system will also have a free
    Hamiltonian
    \f$H_{\mathrm{S}}=\frac{U}{2}n_{\uparrow}n_{\downarrow}+V\sum_{\sigma}
    n_{\sigma}\f$.  
*/

class ThermofieldFermionHamiltonian:
public ThermofieldHamiltonian {

  int d;

 public:
  /** Create with all the required fields as arguments. Different to
      ThermofieldHamiltonina, the free Hamiltonian of the system is
      constructed from paramevers U and V, and also the operator in
      the coupling term has fixed form. */
  ThermofieldFermionHamiltonian(int L,
				const std::vector<double>& a1n,
				const std::vector<double>& a2n,
				const std::vector<double>& b1n,
				const std::vector<double>& b2n,
				const double U,const double V);

  ~ThermofieldFermionHamiltonian(){};
  /** Construct the MPO corresponding to the evolution operator with
      the splitting even-odd in the Hamiltonian. Take as argument the
      complex factor in the exponential, delta.  For real time
      evolution, delta is -i*tau, for imaginary time, it should be
      -tau.
      The even term includes boson pairs 01, 23,... AND the free
      Hamiltonian of the central spin. The odd term includes the
      bosonic pairs 12,34... AND both coupling terms between spin and
      zero-modes, which we assume commute.*/
  //  void getExponentialMPOeven(MPO& expHe,complex_t delta);
  // void getExponentialMPOodd(MPO& expHo,complex_t delta);

 private:
  void constructOperatorProduct(mwArray& result,const mwArray& opA,
				const mwArray& opB) const;
  void constructFreeHamiltonian(const double U,const double V);

 protected:
  virtual void initHMPO();
  virtual void initZ();

  void getTwoBodyTermBathExponential(mwArray& Ol,mwArray& Or,complex_t delta,int k,
				 bool left) const;
  /** Compute the two body term corresponding to the h-th and (k+1)-th
      mode in the left or roght baths. */
  void computeTwoBodyTermBath(mwArray& result,int k,bool left) const;
  void getCouplingExponential(mwArray& Opl,mwArray& Ospin,mwArray& Opr,
			      complex_t delta) const;
  void setSingleModeTerm(MPO& expH,int pos,int dim,
			 complex_t delta,bool left) const;
  void split3term(mwArray& Ol,mwArray& Oc,mwArray& Or,mwArray& expH,
		  int dimL,int dimC,int dimR,int& Dl,int& Dr) const;
};

#endif // THERMOFIELDFERMIONHAMIL_H
