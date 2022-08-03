#ifndef NMPO_H
#define NMPO_H

#include "MPO.h"

/** Construct some  special MPOs, for a (finite) bosonic
    lattice in the occupation number basis.
    Basis on each site: 
    \f$\{|0\rangle,|1\rangle,\ldots|N_{max}\rangle\} \f$. */

class BosonMPO{
 public: 
  /** Create, specifying the length of the chain and the 
      maximum occupation number per site, the operator 
      total number of bosons: 
      \f[ 
      \hat N=\sum_i \hat n_i= \sum_i b_i^{\dagger} b_i
      \f]
  */
  static void getNumberMPO(int L,int Nmax,MPO& numberMPO);

  /** Creates the MPO for the "string" operator \f[
      exp\left(i\pi\sum_{i\leq k \leq j} (\hat{n}_k-\bar{n})\right).  \f] 
      If \param<avrgN> is different from 0, it is assumed to contain 
      \f$\bar{n}\f$. If it is zero, however, the exponents are computed 
      without \f$\bar{n}\f$, which
      is a constant (depending on the state) and can be included after
      computing the expectation value. 
      The MPO will hold a copy of the identity operator and another one of 
      the single body exponential. If several MPOs are to be used with 
      different values of i,j, the memory usage could be slightly 
      optimized by reusing these pointers for all of them. 
  */
  static void getParityMPO(int L,int Nmax,MPO& parityMPO,int i,int j,
			   double avrgN=0.);


/**
   Return the matrix to be used in the operator to measure the
   occupation number on a bosonic site of given maximum occupation.
 */
  static void getSingleSiteBosonNumber(int Nmax,mwArray& result);

 private:
  BosonMPO(){};
  ~BosonMPO(){};


};

#endif // NMPO_H
