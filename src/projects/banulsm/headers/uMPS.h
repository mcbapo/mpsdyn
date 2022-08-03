/**
   \file uMPS.h
   Basic implementation of a uniform MPS.
   Should replace completely TIMPS.
   
   \author Mari-Carmen Banuls
   \date 27/12/2019
*/

#ifndef UMPS_H
#define UMPS_H

//#include "Indices.h"
//#include "Site.h"
#include "MPS.h"
#include "TensorMultiplier.h"

#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <vector>

#define TOLSVD 1E-8


/** Uniform MPS, made up of repetitions of a unit cell with any number
    of tensors, as in the schematic picture.

    \image html uMPS.png "Unit cell of a uMPS" width=500px

    The uMPS allows a canonical form where the left fixed point (for
    any cut between two tensors) is the identity, and the right fixed
    point is a diagonal matrix. The class imposes this automatically
    before computing an expectation value.

    It supports real or imaginary time evolution by applying an MPO
    (in the form of a vector of tensors that correspond to the
    operators) and truncating again into the canonical form.
*/


class uMPS{
  /** Number of different tensors (i.e. periodicity) */
  int nA;
  
  /** Array of sites (one or more, depending on nA) */
  Site** W;

  /** Array of corresponding diagonal matrices for the fixed point of
      the transfer operators, with left canonical form (Lambda[k] is
      the right fixed point for the transfer operator of tensors
      W[k]-W[k-1]) */
  mwArray** Lambda;

  bool isCanonical_; // to keep internal consistency
  
 public:
  /** Construct specifying number of sites (default 1) */
  uMPS(int d,int D);
  /** For several sites, with uniform physical and bond dimensions  */
  uMPS(int nA,int d,int D);
  /** General: each site can have different phyiscal and bond
      dimension (last right is the same as first left, so the vector
      of bond dimensions is supposed to give the left ones.*/
  uMPS(int nA,std::vector<int>& ds,std::vector<int>& Ds);
  /** Initialize with a list of tensors. Their dimensions must be in
the order (d,D,D), and be consistent with each other. */
  uMPS(std::vector<mwArray>& tensors);
  /** Copy constructor */
  uMPS(const uMPS& state);
  ~uMPS();

  /** Return number of tensors in the unit cell */  
  int getN() const {return nA;}
  /** Return the physical dimension of k-th tensor in the unit cell
      (numbered from 0).*/
  int getPhysicalDimension(int k) const {return W[k]->getd();}
  /** Return the left virtual dimension of k-th tensor in the unit cell
      (numbered from 0).*/
  int getLeftVirtualDimension(int k) const {return W[k]->getDl();}
  /** Return the right virtual dimension of k-th tensor in the unit cell
      (numbered from 0).*/
  int getRightVirtualDimension(int k) const {return W[k]->getDr();}

  /** Return (a copy of) the k-th tensor in the unit cell
      (numbered from 0).*/
  mwArray getA(int pos) const {return W[pos]->getA();}

  /** Return the square root of the right fixed point in the canonical
      form.  Notice that if the state is not in canonical form, this
      fails */
  mwArray getLambda(int pos) const;

  /** Change one of the tensors by hand. */
  void replaceSite(int pos,const mwArray& A);
  
  /** Print data */
  void put(ostream& os) const;
  
  /** Initialize with a given product state */
  void setProductState(ProductState state);

  /** Initialize with random tensors */
  void setRandomState();
  
  /** \todo Initialize with specific tensors?? */
  
  /** Read from a file */
  void importuMPS(const string& filename);
  
  /** Write to a file */
  void exportuMPS(const string& filename) const;

  /** Impose canonical form (if Dcut is 0, w.o. truncation) */
  void canonical(int Dcut=0,double tolSVD=TOLSVD); 
 
  /** Apply a MPO and impose canonical form again. The MPO is
      specified as a vector of tensors, each of them 4-dimensional,
      with legs ordered as in Operator. 

      \image html uMPS_applyMPO.png "Action of MPO on uMPS" width=500px
*/
  void applyMPO(std::vector<const mwArray*>& midU,double tolSVD=TOLSVD,int Dcut=0);

  /** Same as in applyMPO, but without any truncation: only exact
      application of the operators, increasing the virtual dimensions. */
  void applyMPOnoCanonical(std::vector<const mwArray*>& midU);

  /** Compute the expectation value of a product of local operators at
      particular sites (0 to nA-1). The first value in the vector of
      positions indicates where to start the block: the others have to
      be to the right (e.g. if nA=3, with tensors ABC, and pos={2,1},
      this means that the block taken is CAB, with ops at C and B)
  
      \todo Add the possibility of estimating long range
      correlations which go beyond the neighbouring unit cell. */
  complex_t computeExpectationValue(int nOp,const std::vector<int>& pos,
				    const std::vector<mwArray*>& op);
  

  /** For single site expectation values it is more convenient to
      recover the RDM for each single site.*/
  void getSingleSiteRDM(int pos,mwArray& rdm);

  /** Return the entropy for the cut to the left of pos (i.e. using
      the corresponding Lambda) */
  double getEntropy(int pos=0);

  /** In case of a purification, return the largest eigenvalue of the
      transfer operator for rho^2*/
  void getPurity(complex_t& eta);


  /** Get tr(rho_L^2) for the RDM of L sites */
  void getPurityBlock(complex_t& etaL,int L);

 protected:

  bool isCanonical() const{return isCanonical_;}

  // Get a copy of the lambda, irrespective of canonicity
  mwArray copyLambda(int pos) const{return *Lambda[pos];}

 private:
  void clear();
  /** Auxiliary functions used inside the application of the canonical form */
  void applyCyclicCanonical(int k0,std::vector<Multiplier*>& E,std::vector<Multiplier*>& El,int Dcut,double tolSVD);
  /** Compute the largest right eigenvalue of a product of
      multipliers, but assuming that the first one is anywhere (k0) in
      the vector, then wrapping around. If reverse==true, the wrapping
      is done backwards, as corresponds to the vector for El, but the
      multiplication is always applied from the rightmost one to the left.  */
  void getLargestEigenvalueCyclic(std::vector<Multiplier*>& M,int k0,mwArray& V,complex_t& lambda,bool reverse=false) const;
  complex_t getFirstNonZeroDiag(const mwArray& X) const;

  /** Special case: if only one tensor, I allow all kind of
      correlations to be computed. */
  complex_t computeExpectationValueSingle(int nOp,const std::vector<int>& pos,
				    const std::vector<mwArray*>& op) const;

};

ostream& operator<<(ostream& os, const uMPS& mps);

#endif // UMPS_H
