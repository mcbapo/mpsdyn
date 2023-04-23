#ifndef SU2_NO_LINKS_DIMLESS_H
#define SU2_NO_LINKS_DIMLESS_H

#include "MPS.h"
#include "Hamiltonian.h"
#include "Contractor.h"
#include "LGTstring.h"
#include "SU2HamiltonianNoLinks.h"

/**
   \file SU2HamiltonianNoLinksDimless.h Definition of the class that implements a truncated SU2 Hamiltonian following the reduced gauge invariant basis presented in the paper "Hamer, C. SU(2) Yang-Mills theory in (1 + 1) dimensions: A finite-lattice approach Nuclear Physics B , 1982, 195, 503 - 521". Contrary to the version in ReducedSU2Hamiltonian.h, this one does not explicitly contain the gauge links, but instead the fermionic sites have one degree of freedom added to encode the gauge field. This version is furthermore dimensionless such that it is suitable for calculations which should be extrapolated to the continuum
   
   \author Stefan Kühn
   \date 23/10/2015
*/

class SU2HamiltonianNoLinksDimless: public SU2HamiltonianNoLinks
{
  
protected:
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously)*/
  void initHamiltonian(void);
  
public:
  /** Creates an instance of the truncated SU(2)-Hamiltonian in the reduced basis. The dimensionless Hamiltonian is given by
   \f[ W = \frac{2}{ag^2}H = xV + \sum_n 2\frac{m}{g}\sqrt{x} \f]
   The parameter mg is specifying the value for \f$ \frac{m}{g}\f$.*/
  SU2HamiltonianNoLinksDimless(int N, double x, double mg, int dlink, double lambda, double sz_penalty);
  
  ~SU2HamiltonianNoLinksDimless();
  
  
};

#endif //SU2_NO_LINKS_DIMLESS_H
