#ifndef SOFTCOUL_H
#define SOFTCOUL_H

#include "Hamiltonian.h"
#include <vector>

/** 
    Implementation of SoftCoulomb Hamiltonian for fixed number of
    sites, using an MPO approximation as a sum of a certain number of
    exponentials. We include the local potential corresponding to the
    nucleus (total charge Z and located in the center of the system,
    in position \f$l_Z\f$).
     \f[
    H=\sum_{i=1}^{L-1} \frac{1}{2 \Delta^2} \left ( c_{i,\sigma}^{\dagger} c_{i+1,\sigma} +
    c_{i+1,\sigma}^{\dagger} c_{i,\sigma} \right ) -
    \sum_{i=1}^{L} \frac{Z}{\sqrt{((l_Z-i) \Delta )^2}+(a \Delta)^2} n_i+
    \sum_{i=1}^{L-1}\sum_{j>i} \frac{n_i n_j}{\sqrt{(i \Delta -j \Delta )^2}+(a \Delta)^2}
    +\sum_{i=1}^{L} \mu n_i
    +\eta \left( \sum_{i=1}^{L} n_i -N_e\right)^2
    +\eta_S \vec{S}^2
    \f]
    Notice it is defined in terms of spin matrices, not sigmas.
*/

class SoftCoulombHamiltonian:
public Hamiltonian{
  int L; // chain length
  int d; // dimension
  MPO hamil;
  mwArray Z;
  int M,Ze; // number of exponentials
  std::vector<complex_t> lambdas; // exponential factors
  mwArray X; // coefficients of the exponentials
  double Delta,a; // input parameters
  bool chemicalPotential; // whether there is a chemical potential
  bool penaltyNe; // whether there is a penalty term (N-Ne)^2
  bool penaltyS; // whether there is a penalty term for total spin!=0
  double eta; // mu or eta
  double etaS; // penalty for spin
  int Ne; // total N wanted
 public:
  /** Create a SoftCoulomb Hamiltonian for a chain of length L, using
      lattice spacing Delta, nuclear charge Z (in the center),
      softening parameter a, and using a total of M exponentials to
      approximate the power law decaying interactions in the MPO.  If
      eta!=0 is given, and no Ne value is provided, a chemical
      potential term with mu=eta is addded to the Hamiltonian. If Ne
      is also given, then we add a penalty term to ensure total number
      of electrons equal to Ne, and eta is then the penalty
      coefficient. */
  SoftCoulombHamiltonian(int L,double Delta,int Z,double a,int M,double eta=0,int Ne=0,double etaS=0);

  ~SoftCoulombHamiltonian();

  const MPO& getHMPO() const {return hamil;}

  bool hasHMPO() const {return true;}

  /** Returns (in mpo) the MPO corresponding to the total number of electrons 
      (bond dimension 2). */
  void getNumberMPO(MPO& mpo) const;

/** Return the mpo corresponding to \f$\vec{S}^2\f$, i.e. total spin
 * (bond dimension 5). */
  void getSpinMPO(MPO& mpo) const;

/** Return the single body operator that computes the numer per site. */
  void getLocalNumberOp(mwArray& oper) const;

 private:
  // compute the coefficients for the exponential terms
  void computeCoefficients();
  void initZ();
  void initHMPO();
  // auxiliary function (copied from LiouvillianXYedges) to compose a double site operator
  void constructOperatorProduct(mwArray& result,const mwArray& opA,
				const mwArray& opB) const;
};

#endif // SOFTCOUL_H
