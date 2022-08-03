#ifndef NNHAMIL_H
#define NNHAMIL_H

#include "Hamiltonian.h"

/** 
    A Hamiltonian containing only nearest-neighbor terms, which can be
    different, for a finite chain of given length. The physical
    dimensions can also be different for each site.
*/

class NNHamiltonian:
  public Hamiltonian{
  int L; // Nr of sites
  std::vector<int> d; // Physical dimensions
  std::vector<mwArray> terms; // Each of the terms, from 0 to L-1
  MPO mpoH; // Hamiltonian MPO
  bool uptodate; // whether the MPO is ready

public:
  /** Create empty Hamiltonian (make private?) */
  NNHamiltonian();
/** Create with given length and dimensions */
  NNHamiltonian(int L,const std::vector<int>& dims);
  /** Create reading from text file that has written from Matlab */
  NNHamiltonian(const char* filename);
  ~NNHamiltonian();

  /** Return the length */
  int getLength() const {return L;}

  /** Return the vector of dimensions in the given container. */
  void getDimensions(std::vector<int>& dims) const;

  bool hasHMPO() const{return true;}
  const MPO& getHMPO() const;
  bool hasUMPO() const {return true;}
/** Specify one of the NN terms in the Hamiltonian, acting on site
 * pos and pos+1. The mwArray is assumed to be given with dimensions
 * d1 d2 x d1 d2 (for physical dimensions d1 and d2 of pos and pos+1). Afterwards, the MPO for H is not ready, but needs to be recomputed explicitly via recomputeMPO. */
  void setTerm(int pos,const mwArray& term);

/** Return a copy to a certain term (acting on sites pos and
 * pos+1). The mwArray will have dimensions d1 d2 x d1 d2 (being d1
 * the physical dimension of site pos) */
  mwArray getTerm(int pos) const;

  /** Construct the MPO corresponding to the evolution operator with
      the splitting even-odd in the Hamiltonian. Take as argument the
      complex factor in the exponential, delta.  For real time
      evolution, delta is -i*tau, for imaginary time, it should be
      -tau. */
  void getExponentialMPOeven(MPO& expHe,complex_t delta);
  void getExponentialMPOodd(MPO& expHo,complex_t delta);

/** Return a MPO for a single term in the Hamiltonian (i.e.,q
 * everywhere the identity except on two sites, pos and pos+1). */
  void getMPOforTerm(MPO& localTerm,int pos) const;

/** Recompute the MPO for the Hamiltonian after changing (or maybe assigning for the first time) the individual two-body terms.
 */
      void recomputeMPO();
      
  private:
  /** Prepare the even-odd decomposition of the evolution operator:
      returns the two terms to be inserted in the MPO for the 
      exp(delta H12) (delta a complex number).
      The exponential of the term starting in pos is returned.  
  */
  void getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,
				 int pos=0);
  void initHMPO();
  //void initHMPO(double* tis=0,double* Uis=0,double* muis=0,
  //	 double* Vis=0);
  /** Common computation for recovering the even and odd parts of the exp */
  void getExponentialMPO(MPO& expHo,complex_t delta,bool even);

/** Split a two body term, via a SVD, in sum of product of local
 * operators. d1 and d2 are the (left and right) physical dimensions
 * according to which the svd is performed. */
  void splitTwoBodyTerm(const mwArray& term,int d1,int d2,
			mwArray& Ol,mwArray& Or,int& nr) const;
};

#endif
