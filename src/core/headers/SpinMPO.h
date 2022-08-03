#ifndef SMPO_H
#define SMPO_H

#include "MPO.h"

/** Construct some  special MPOs, for a (finite) spin
    chain, given the length and the local spin.
    Currently, only for spin 1/2.
*/

class SpinMPO{
 public: 
  /** 
      Total z spin component of the chain, given the length and the
      local spin dimension 2*s+1 (s=1/2).
      If staggered==true, the staggered version \f$\sum_i (-1)^{i}
      S_{i}^{\alpha}\f$ is computed instead (notice that the first
      site is numbered 0).
  */
  static void getSxMPO(int L,int spindim,MPO& Sx,bool staggered=0);
  static void getSyMPO(int L,int spindim,MPO& Sy,bool staggered=0);
  static void getSzMPO(int L,int spindim,MPO& Sz,bool staggered=0);

  /** Operator Sz^2, in a more economic way than squaring the previous
      MPO for Sz. */
  static void getSz2MPO(int L,int spindim,MPO& Sz2);
  /** Idem for Sx^2 */
  static void getSx2MPO(int L,int spindim,MPO& Sx2);
  /** Idem for Sy^2 */
  static void getSy2MPO(int L,int spindim,MPO& Sy2);

  /** Sum of nearest-neighbor correlations \f$\sum_i S_i^{\alpha} S_{i+1}^{\alpha}\f$*/
  
  static void getSxSxMPO(int L,int spindim,MPO& SxSx);
  static void getSySyMPO(int L,int spindim,MPO& SySy);
  static void getSzSzMPO(int L,int spindim,MPO& SzSz);

  
  /** Idem for \f$S^2=Sx^2+Sy^2+Sz^2\f$. If a fourth parameter alpha
   * is given, \f$S^2-alpha\f$ is constructed instead */
  static void getS2MPO(int L,int spindim,MPO& S2,double alpha=0.);


/** Projector onto total Sz=0 */
  static void getProjectorTotalSz(int L,int spindim,MPO& PSz0);


  /** Projector onto total spin 1 of each pair of neighbors.  If
      even==true, starts on 0.*/
  static void getProjectorTotalSPair(int L,int spindim,MPO& PS1,
				     bool even=true);

  /** 
      Construct an MPO (bond dimension \f$d^2\f$) that represents a
      cyclic translation in the chain (i.e., everyone is moved one to
      the right, and the last site is translated to the first one).
      The result is returned in T.
   */
  static void getCyclicTranslationMPO(int L,int dim,MPO& T);

/** Construct the MPO that measures the ocupation number of a given
 * momentum mode k (\f$k=-L/2,\ldots, 0,\ldots L/2-1\f$), i.e. corresponding to
 * quasimomentum \f$\frac{2\pi}{L} k\f$. If the MPO was alredy
 * constructed, the flag prepared can be set to true, and instead of
 * constructing the MPO from scratch, only the appropriate components
 * will be changed. */
  static void getMomentumDistributionMPO(int L,int dim,int k,MPO& Nkmpo,
					 bool prepared=false);

  /** Construct the quasimomentum mode
      \f[S_k^-=\frac{1}{\sqrt{N}}\sum_m e^{-i k m} \sigma_m^-\f], with
      \f$k=\frac{2\pi z}{N}\f$, for \f$z=0,\ldots N-1\f$.
   */
  static void getSkMPO(int N,int dim,int z,MPO& Skmpo,bool prepared=false);


  /** Construct the quasimomentum mode (for open chain)
      \f[S_k^-=\sqrt{\frac{2}{N+1}}\sum_m sin({ k m}) \sigma_m^-\f], with
      \f$k=\frac{\pi z}{N+1}\f$, for \f$z=1,\ldots N\f$.
   */
  static void getOBCSkMPO(int N,int dim,int z,MPO& Skmpo,bool prepared=false);

  //private:
  /** Construct an MPO for the sum, with equl weights, of operator
      spinOp acting on each site of a chain with length L, with idOp
      acting on every other site, e.g. \f$\sum_i \sigma_x^{(i)}\f$.

      If staggered==true, the staggered version \f$\sum_i (-1)^{i}
      S_{i}^{\alpha}\f$ is computed instead (notice that the first
      site is numbered 0). The same could be obtained with getModulatedSingleBodyMPO.
 */
  static void getSingleBodyMPO(int L,const mwArray& idOp,
			       const mwArray& spinOp,MPO& res,bool staggered=0);

  /** Construct a superposition of the given operator on each site, 
      modulated with phase \f$e^{-i q \pi m}\f$ on each site.
      If obc is true, instead of the exponential, the coefficient is \f$sin({ q \pi m})\f$.
   */
  static void getModulatedSingleBodyMPO(int L,const mwArray& idOp,
					const mwArray& spinOp,double q,MPO& res,bool obc=false);


  /** Construct the MPO for a nearest-neighbor correlator (averaged),
      such as \f$\sum_i S_i^y S_{i+1}^z\f$. To keep it flexible, two
      different operators \f$A_i\f$ and \f$B_{i+1}\f$ are
      specified.  */
  static void getNearestNeighbourMPO(int L,const mwArray& idOp,
				    const mwArray& opA,const mwArray& opB,
				    MPO& res);

};

#endif // SMPO_H
