#ifndef SU2_NO_LINKS_Q_H
#define SU2_NO_LINKS_Q_H

#include "MPS.h"
#include "Hamiltonian.h"
#include "Contractor.h"
#include "LGTstring.h"
#include "SU2HamiltonianNoLinks.h"

/**
   \file SU2HamiltonianNoLinksQ.h Definition of the class that implements a truncated SU2 Hamiltonian following the reduced gauge invariant basis presented in the paper "Hamer, C. SU(2) Yang-Mills theory in (1 + 1) dimensions: A finite-lattice approach Nuclear Physics B , 1982, 195, 503 - 521". Contrary to the version in ReducedSU2Hamiltonian.h, this one does not explicitly contain the gauge links, but instead the fermionic sites have one degree of freedom added to encode the gauge field. This version additionally penalizes states that do not have a total abelian charge of zero.
   
   \author Stefan Kühn
   \date 20/10/2015
*/

class SU2HamiltonianNoLinksQ: 
public SU2HamiltonianNoLinks
{
private:
  double charge_penalty;
  
protected:
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously)*/
  void initHamiltonian(void);
  
public:
  /** Creates an instance of the truncated SU(2)-Hamiltonian in the reduced basis.*/
  SU2HamiltonianNoLinksQ(int N, double epsilon, double mu, double g, int dlink, double lambda, double sz_penalty, double charge_penalty);
  
  ~SU2HamiltonianNoLinksQ();
  
  
};

#endif //SU2_NO_LINKS_Q_H
