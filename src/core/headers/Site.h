/**
   \file Site.h
   Definition of the elementary class containing the tensor for a
   single site in a MPS 
   
   \author Mari-Carmen Banuls
   \date 16/03/2011
*/

#ifndef SITE_H
#define SITE_H

//#include "mpsA.h"
//#include "matlab/libdiminf.h"
//#include "libdiminf.h"

#include "Indices.h"
#include "wrapper.h"
#include "mwArray.h"
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>

//typedef int dimen_t;

using namespace std;

/** A list of single site states that I can construct */
enum SingleState {
  zero=0, // |0> 
  one,  // |1>
  xplus,  // |x+>=|0> + |1>
  xminus, // |x->=|0> - |1>
  yplus,  // |y+>=|0> + i |1>
  yminus,  // |y->=|0> - i |1>
  special,
  maxent,  // |00>+|11> for local dimension 4
  randmps   // Fully random state, with the same dimensions
};


/** 
    The class that contains the A tensor for a particular site 
    in the chain, with dimensions (d,Dl,Dr).
    \image html Site.png
*/

class Site {

 protected:
  /** The position of the site */
  int pos;
  /** The dimensions, as d (physical), Dl (left bond), Dr (right bond) */ 
  //  vector<int> dims;
  /** Dimensions of left and right bonds. 
     TODO: These data are redundant (already in mwArray)=> Eliminate! */

 protected:

  /** The A tensor   */
  mwArray A; 

  friend class MPS;
  friend class TIMPS;
  friend class uMPS;
  friend class Contractor;

 public:
  /** Create a Site (default is state |0>) */
  Site(int pos,int d,int Dl,int Dr);
  /** Create a Site passing an array with dimensions. 
      \warning The length of the array MUST BE 3 */
  Site(int pos,int* dimensions);
  /** Copy constructor. Allow a different position to be set */
  Site(const Site& s,int pos=0);
  /** Create a Site for the given position, with a copy of this array as
     data. The dimensions of the array are taken as dxDlxDr */
  Site(int pos,const mwArray& A);
  virtual ~Site();

  // Return values
  virtual int getd() const{return A.getDimension(0);}
  virtual int getDl() const{return A.getDimension(1);}
  virtual int getDr() const{return A.getDimension(2);}

  int getPos() const {return pos;}
  virtual shrt::Indices getDimensions() const {return A.getDimensions();}

#ifndef DISKSTRG  
  const mwArray& getA() const {return A;}
#else
  // Then i need to override it by fileSite
  virtual const mwArray getA() const {return A;}
#endif

  // Copy this data.
  void setSite(const mwArray& newA);
  // Copy this data but permuting its dimensions as newdims indicates
  void setRotatedSite(const mwArray& newA,const shrt::Indices& newdims,
		      bool conjugate=0);
  // Set the tensor to all zeroes, keeping the dimensions 
  void setToZero();

  // TODO?: Ability just to store a given array (reduce number of copies!)

  /** Debugging utilities */
  friend ostream& operator<<(ostream& os, const Site& site);

  /** Set the array so that it represents the state of one particle given 
     as argument. */
  void setState(SingleState state);
  /** Overload assignment */
  Site& operator=(const Site& s);

  /** Read from a binary file */
  virtual void load(ifstream& data);
  /** Save to a binary file */
  virtual void save(ofstream& data) const;
  /** Read from a text file */
  virtual void loadtext(ifstream& data);
  /** Save to a text file */
  virtual void savetext(ofstream& data) const;

#ifndef TESTINGSITE
 protected:
#endif
  // Empty constructor for derived class FileSite
  // Almost empty constructor for 2D-Site
Site(int pos_):pos(pos_),A(){};
/** Auxiliary method for output*/
  virtual void put(ostream& os) const;

  /** Apply the gauge condition to the right. This should only be done
      by MPS, as it has no sense on an isolated tensor.  Argument
      rightterm is a reference that upon return contains the tensor to
      be multiplied on the site to right, and leftterm is the same
      term received from the site to the left.  If cutD is received,
      it is used to cut off the SVD, so that the resulting Site will
      have Dr<=D.  If cutD is not specified, instead of the full
      SVD,only a QR decomposition is done, and the new tensor will be
      Q, while the rightterm will be R.
  */
  virtual void gaugeR(mwArray& rightterm,const mwArray& leftterm,int cutD=0);

  /** The same functionality for gauge condition to the left.  If cutD
      is received, it is used to cut off the SVD, so that the
      resulting Site will have Dl<=D. If no cut is required, only a LQ
      decomposition is done, and the full SVD is not computed.
  */
  virtual void gaugeL(mwArray& leftterm,const mwArray& rightterm,int cutD=0);

  // Special truncation with the transpose (Hastings 2014)
  virtual void gaugeRT(mwArray& rightterm,const mwArray& leftterm,mwArray& lambda,int cutD);
  virtual void gaugeLT(mwArray& leftterm,const mwArray& rightterm,mwArray& lambda,int cutD);

  
  /**
     Given tensors A and B, impose a gauge condition for <B|A> as
     moving to the right (B+ A =I). This is applied by multiplying
     these and the next tensors by an appropriate X (Y) and its
     inverse, respectively: tensorA (B) returns the factor to be
     multiplied on the right, whereas lefttermA (B) should contain the
     corresponding one from the gaugeABr operation on the site to the
     left.  For the first site, leftterm is not received, while for
     the last one (Dr=1) tensor only contains a global normalization
     factor.  \warning IMPORTANT: For the function to work as
     expected, both A and B have to satisfy gaugeCond to the right,
     prior to this call (they are assumed both to be full rank)
 */
  static void gaugeR(Site* ket,Site* bra,mwArray& righttermK,
		     mwArray& righttermB,const mwArray& lefttermK,
		     const mwArray& lefttermB);
  /** Same but moving from right to left (ket and bra need to satisfy gauge cond to left) */
  static void gaugeL(Site* ket,Site* bra,mwArray& lefttermK,
		     mwArray& lefttermB,const mwArray& righttermK,
		     const mwArray& righttermB);

  
  // Apply the corresponding dimension reduction if the gauge condition 
  // requires it.
  virtual void reduceDimR(const mwArray& leftterm);
  virtual void reduceDimL(const mwArray& rightterm);

  // To keep track of a pi phase when converting from an extern set of 
  // tensors that is then subject to gauge condition. This is only useful 
  // for test programs.
  void changeSign();

  // Complex conjugate the tensor. Has to be added to the file version.
  virtual void conjugateA();
  
  /** Allow increasing the bond dimension, but only from MPS (so that it 
     is consistent along the whole state) */
  virtual void increaseBondDim(int Dl,int Dr);
  /** We could also cut off part of the virtual indices */
  virtual void decreaseBondDim(int Dl,int Dr);
  /** Idem for the physical dimension */
  virtual void increasePhysDim(int d);
  /** Cut the higher indices of d, in order to adjust the Site to a given 
     physical dimension. This should be used ONLY before any optimization 
     routine, in order to produce an initial MPS of the proper dimensions. */
  virtual void decreasePhysDim(int d);

  /** Increase the bond dimension, adding random noise to the extra elements */
  virtual void increaseBondDimWithNoise(int Dl,int Dr,double noise);


 protected:
  void setd(int x){exit(1);}
  void setDl(int x){exit(1);}
  void setDr(int x){exit(1);}
  void setPos(int x){pos=x;}

  void changeDimensions(int d,int Dl,int Dr,double noise=0.);

  /** Apply gauge condition to the left, at the same time an operator
      (or its adjoint, if dagger==true) is applied and the resulting
      bond dimension is truncated according to cutD.  If cutD is
      received, it is used to cut off the SVD, so that the resulting
      Site will have Dl<=D.
  */
  virtual void gaugeLOp(mwArray& leftterm,const mwArray& op,
		const mwArray& rightterm,bool dagger=false,
		int cutD=0);
  // Associated contraction of a remaining term
  virtual void reduceDimLOp(const mwArray& op,const mwArray& rightterm,
		    bool dagger=false);

  /** Same to the right
  */
  virtual void gaugeROp(mwArray& rightterm,const mwArray& op,
		const mwArray& leftterm,bool dagger=false,
		int cutD=0);
  // Associated contraction of a remaining term
  virtual void reduceDimROp(const mwArray& op,const mwArray& leftterm,
		    bool dagger=false);

};


#endif //SITE_H
