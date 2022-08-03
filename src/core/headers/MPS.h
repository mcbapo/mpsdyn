/**
   \file MPS.h
   Basic class containing an MPS 
   
   \author Mari-Carmen Banuls
   \date 16/03/2011
*/

#ifndef MPS_H
#define MPS_H

#include "Indices.h"
#include "Site.h"
#include "Operator.h"
#ifndef TESTINGMPS
#include "MPO.h"
#endif
#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <vector>

#ifdef DISKSTRG
#include "FileSite.h"
typedef FileSite Site_t ;
#else
typedef Site Site_t ;
#endif


/** List of product states that I can construct */

enum ProductState{
  p_zero=zero, // |0> 
  p_one,  // |1>
  p_xplus,  // |x+>=|0> + |1>
  p_xminus, // |x->=|0> - |1>
  p_yplus,  // |y+>=|0> + i |1>
  p_yminus,  // |y->=|0> - i |1>
  p_special,  // if d=d0^2, max entangled of subsystems 
  p_maxent  // if d=d0^2, max entangled |00>+|11>+...+|d0-1 d0-1> 
};

/** List of gauge conditions an MPS can satisfy. This is to keep track
    of a previously appplied gauge condition, in case it saves us some
    time. It is immediately lost by replacing a site by hand, or
    applying an individual gauge condition to a certain site. */
enum Gauge{
  gNone,
  gLeft,
  gRight,
  gBoth   /* in special cases (product states), both left and right at
	     the same time*/
};

/** 
    The class that contains a (general) MPS state. 
*/

class MPS{

  /** The length of the chain (\todo a negative value could indicate 
     infinite length) */
  int length;
  /** Maximum dimension of the bond (\todo negative value could be used to 
     indicate no limit). */
  int bond;
  /** Array of sites (one per physical position in the chain) */
  Site_t** A;

  /** The normalization factor */
  double normfact;

  /** Which kind of gauge condition (if any) this MPS satisfies.
     (no gauge condition, to the left, or to the right).*/
  Gauge gauge;

 public:
  /** Default constructor */
  MPS();

  /** Creates a default MPS, with all sites of dimension d, initialized 
     to product state (x)|0> */
  MPS(int length,int bond,int d);
  /** Creates a MPS where the physical dimension of sites is specified 
     in vector d (of length length) and the bond dimensions in bond 
     (length-1) */
  MPS(int length,int bond[],int d[]);
  /** Creates a default MPS, with dimension d[k] for each site, initialized 
     to product state (x)|0> */
  MPS(int length,int bond,int* d);
  // Idem, with an STL vector for the physical dimensions
  MPS(int length,int bond,const std::vector<int>& d);
  // Idem, with an STL vector for both the physical and bond  dimensions
  MPS(int length,const std::vector<int>& bond,const std::vector<int>& d);
  /** Creates a copy of the MPS given as argument. */
  MPS(const MPS& aMPS);
  /** Creates a copy of the MPS given as argument, using bond as maximum . */
  MPS(const MPS& aMPS,int bond);
#ifndef TESTINGMPS
  /** Creates an MPS corresponding to the action of ops on the MPS given 
     as argument, using bond as maximum. If dir='U', the operator is 
     assumed to act on the bra. */
  MPS(const MPO& ops,const MPS& aMPS,int bond,char dir='D',char gauge='L');
#endif
  ~MPS();

  /** Length of the chain */
  int getLength() const {return length;}
  /** Maximum bond dimension */
  int getBond() const {return bond;}
  /** Retrieve a particular Site_t (numbering of sites from 0 to N-1!) */
  const Site_t& getA(int k) const;
  /** Return true iff the gauge condition to the left is guaranteed. */
  bool isGaugeL() const{return (gauge==gLeft)||(gauge==gBoth);}
  /** Return true iff the gauge condition to the right is guaranteed. */
  bool isGaugeR() const{return (gauge==gRight)||(gauge==gBoth);}

  /** Return true if the MPS is empty (result of default construction or 
      cleared)*/
  bool isEmpty() const{return length==0;}

  /** Set the A tensor of the particular site pos (numbering from 0 to
     N-1) to the value in data (which is copied). Requires that the
     total dimensions of the argument agree with the already stored
     ones (i.e. same product d x Dl x Dr). The mwArray is then
     reshaped to keep the same dimensions on the Site.
  */
  void setA(int pos,const mwArray& data);

  /** Set the A tensor of the particular site pos (numbering from 0 to N-1) 
     to the value in data (which is copied) after permuting the indices 
     of the latter appropriately, as indicated by newdims. The values of
     d, Dl and Dr are obtained from this rotation. 
     If conjugate is set to true, the data are also conjugated.
     If checkdims==1, Dl and Dr are checked for
     compatibility with neighbouring sites.
  */
  void setRotatedA(int pos,const mwArray& data,const shrt::Indices& newdims,
		   bool conjugate=0,bool checkdims=1);

  /** Unless the former method, this one completely replaces an existing
     Site_t with the argument. If checkdim is set to 0, no check on the new 
     Dl, Dr or physical d is made, so that arbitrary new dimensions can be
     set by this method. If checkdim==1, it checks only that Dl and Dr 
     dimensions agree with those of the neighbouring sites, but allowing a 
     change of physical dimension (useful for constructing MPS from shorter 
     ones) */
  void replaceSite(int pos,const Site_t& data,bool checkdim=1);
  /** This second version allows a similar change, taking as argument 
     directly the data to be copied inside the new Site_t. */
  void replaceSite(int pos,const mwArray& data,bool checkdim=1);

  /** Complex conjugate all sites */
  void conjugateMPS();
  
  /** Apply a local operator on the specified site. This does not change
   * the bond dimension, but it changes the state (and potentially could
   * also change the physical dimension in this site). If normalize is
   * set to true, the resulting state is normalized again (by gaugeR). */
  void applyLocalOperator(int pos,const mwArray& oper,bool normalize=false);

  /** Functions to handle normalization factors and already calculated 
     contraction of the tensors. The capability of changing this values 
     from the outside is not safe. It should be changed. */
  void setNormFact(double value){normfact=value;}
  double getNormFact() const {return normfact;}

  /** Set the MPS to a particular product state */
  void setProductState(ProductState state);

  /** Set the MPS to a random one by assigning random A tensors to
     every site (with the already given dimensions). Then apply gauge
     condition to the right to the random MPS, so that at the end, the
     resulting vector is normalized to one. */
  void setRandomState();

  /** Debugging utilities */
  friend ostream& operator<<(ostream& os, const MPS& mps);

  /** Apply the gauge condition to the right (dir='R') or to the left
     (dir='L'), and if normalize=true, leave the state normalized
     (otherwise, the normalization factor at the end is the norm of
     the state, so that this operation does not change the norm).  If
     cutD different from zero is given, this value is used to cut off
     the SVD along the whole procedure, keeping a maximum of cutD
     values at each bond (so that bond<=cutD).
  */
  void gaugeCond(char dir,bool normalize=false,int cutD=0);
  /** Apply the condition to a single site. 
     next indicates whether the next site in the proper direction 
     has to be updated, i.e. multiplied by the matrix that results 
     from applying the local gauge condition to its neighbour 
     (it is always done if the dimension changes). */
  void gaugeCond(int pos,char dir,bool next=false);

  /** Apply a different gauge condition, that imposes normalization
     of the product of two states (as <bra|ket>=1). */
  static void gaugeCond(MPS& ket,MPS& bra,char dir);

  /** Special truncation with the transpose. */
  void gaugeCondT(char dir,int cutD);

  
  /** Overload assignment */
  MPS& operator=(const MPS& mps);

  /** Export this MPS to the given file, as binary data.
  */
  void exportMPS(const char* filename) const;
  void importMPS(const char* filename);
  /** Auxiliary function, to export as text so that I can import it
      from the other library. Also more portable than binary writing. */
  void exportMPStext(const char* filename) const;
  void importMPStext(const char* filename);

  /** To elp debugging, this function wirtes to a text file the Matlab
      instructions to construct a cell array with the MPS tensors
      (that latr can be expanded to a vector in Matlab with expandVec,
      for instance). The optional argument precision can be used to
      specify the precision with which values are stored */
  void exportForMatlab(const char* filename,int precision=0) const;


  /** Construct this MPS using the tensors of an existing (shorter) MPS such 
     that the edges are the same and the center is filled in with copies of 
     the tensor at reppos. */
  void stretchMPS(int len,const MPS& mps,int reppos=0);

  /** Increase maximal bond dimension */
  void increaseBondDimension(int D);

  /** Increase maximal bond dimension, but instead of padding with
      zeros, add random noise to the extra elements. */
  void increaseBondDimensionWithNoise(int D,double noise);

  /** Increase the physical dimensions of the sites by embedding the MPS 
     tensors into larger ones (padded with zeroes). If the provided 
     dimensions are smaller than the existing ones, no transformation 
     is done. */
  void adjustPhysDimensions(const vector<int>& d);

  /** Clear all contained data.  */
  void clear();

  /** Trace out a list of sites and produce a shorter MPS. It will
      fail if the dimension of some of the site to be traced out is
      not a perfect square. */
  void traceOutSites(const std::vector<int>& pos,bool normalize=false);

  /** Insert a site in the MPS in position pos, increasing the length
      in one.  If checkdim is true (default), the bond dimensions are
      checked for consistency, i.e. new Dl and Dr should agree with
      the bond dimension beetween previous sites at pos-1 and pos. If
      normalize=true (default is false) the gauge (L) is imposed after
      the insertion and the MPS is normalized. */
  void insertSite(int pos,const mwArray& data,bool checkdim=1,bool normalize=false);

#ifndef TESTINGMPS
  /** Assignment of an approximation to the OperatorRow ops acting on
     the MPS aMPS, with max bond X (as in constructon, but can be
     applied at any time, substituting any prior content of MPS. Dir
     indicates whether the operator acts downwards (dir='D', default)
     so that the result is \f$|\Psi\rangle\approx Ops|\Phi_0\rangle\f$
     or upwards (dir='U'), so that \f$|\Psi\rangle\approx
     Ops^{\dagger}|\Phi_0\rangle\f $*/
  void approximate(const MPO& ops,const MPS& aMPS,int X,char dir='D',char gaugeDir='L');
  /* /\** Same, without conjugating. Notice that in this case, the direction indicates how to seew (L meaning from left to right) *\/ */
  /* void approximateT(const MPO& ops,const MPS& aMPS,int X,char dir='L'); */
#endif
  /** As the previous one, but approximates only a MPS */
  void approximate(const MPS& aMPS,int X);

  /** With the given probability, change the values of some of the
     tensors in the MPS, so that the state is randomly changed. */
  void perturb(double probability);

  /** 
      From this MPS, construct another one which is equivalent
      to the folded version of the original. If the original length
      was odd, the central tensor will be now in position 0 and
      length will be (L+1)/2, being L the original length. If L was
      even, new length will be L/2 with every tensor a folded one.
      Returns the new MPS (with copies of the tensors) in \param<mps>.
      In the first implementation, these folded tensors are just the 
      block form of the product. TODO: Support structure! 
  */
  void fold(MPS& mps) const;

  void applyExactly(const MPO& ops,bool dagger=false,
		char dir='L');
 protected:
  /** Return the current gauge, just for the sake of copying state. */
  Gauge getGauge() const{return gauge;}

 private:
#ifndef TESTINGMPS
  /** Auxiliary function to implement the approximation from a
      truncation of MPO on MPS. This is obsolete, as it is betterto
      first apply gauge to exact contraction */
  void truncate(const MPO& ops,const MPS& aMPS,int X,bool dagger=false,
		char dir='L');
  void truncateT(int X,char dir='L');
  //  void approxopL2();
#endif
};

/** Utility function which, given a MPS (the product of physical
 * dimensions should be smaller than 2^20, or it will crash) computes
 * the full vector in a d^L space, with indices ordered d0,d1,d2...
 */
void expandVec(const MPS& mps,mwArray& vec);

/** Given a vector vec, corresponding to a multidimensional array,
 * with dimensions dims (if the current shape does not agree with
 * that, a copy will be done and reshaped) transform it into a MPS
 * with L sites of corresponding physical dimensions
 * dims[0],dims[1],... , and maximum bond dimension D (just by local
 * truncation). */
void vecToMPS(const mwArray& vec,const shrt::Indices& dims,MPS& mps,int D);

/** Utility functions that allow the manipulations of MPO into MPS,
    and viceversa, and the composition of MPO as acting on vectorized
    mixed states, where the physical indices are doubled. */

#ifndef TESTINGMPS
/** Vectorize an MPO to make an MPS out of it. If up==true (default)
    the first component of the double physical index will be the one
    which was going up.*/
void MPSfromMPO(const MPO& mpo,MPS& mps,bool up=true);
/** Reverse operation. If conj=true, the tensors are complex conjugated.*/
void MPOfromMPS(const MPS& mps,MPO& mpo,bool up=true,bool conj=false);

/** doubleMPO taxes a normal MPO, acting on (system) physical
    dimensions, and constructs a double version, as in the evolution
    of a thermal state, for instance, where a second (transposed) copy
    of the MPO acts on the (ancillary) double indices. 
    Namely, this is the operation of folding:  
    MPO (.) MPO -> MPO x MPO^T
    If conjugate==true, the operation folds 
    MPO (.) MPO^{dagger} -> MPO x conj(MPO)   
*/
void doubleMPO(const MPO& simpleMPO,MPO& doubleMPO,bool conjugate=false);
/** extendMPO takes a normal MPO, and transforms it into an MPO acting
    on a mixed state only on one side (the system indices). On the
    rest of the double index (whose dimension must be given in dimA,
    but is typically the same as for the physical) only the identity
    acts.*/
void extendMPO(const MPO& simpleMPO,MPO& doubleMPO,int dimA);
/** Analogous to the one above, but puts the transposed MPO onto the
    ancillary indices*/
void extendTransposeMPO(const MPO& simpleMPO,MPO& doubleMPO,int dimA);

/** Froma a MPS, write a MPO corresponding to an operator with the MPS
    components on the diagonal. This is done by unfolding the physical
    index with a 3-legged delta. */
void diagonalMPOfromMPS(const MPS& mps,MPO& mpo,bool conj);
#endif


#endif // MPS_H
