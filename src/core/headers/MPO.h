/**
   \file MPO.h
   Basic class containing a Matrix Product Operator 
   
   \author Mari-Carmen Banuls
   \date 16/03/2011
*/

#ifndef MPO_H
#define MPO_H

#include <vector>
#include "Indices.h"
#include "Operator.h"

// Help for the constant pointer declaration
typedef const Operator* optr_t;


/** 
    This class contains a row of Operators, that is to be applied to an MPS. 
    It could contain a different operator acting on each site. 
    \class <Contractor> allows MPO to be contracted between two MPS, to 
    find its ground state (for Hermitian MPOs), dominant eigenvectors,
    etc.

    \todo  I will be more interested in cases where the same 
    operator is applied to several sites. For this, I should extend the class, 
    instead of storing multiple copies of the same operator.
*/


class MPO{

 protected:
  /** The length of the chain */
  int length;

  /** Array of operators (one per physical position in the chain) */
  //  Operator** Ops;
  optr_t* Ops;

  // Whether there are Operators that I own, and have to be deleted at the end
  bool myown;
  bool* myops; // which ones

 public:
  /** Default constructor */
  MPO();
  /** Creates a default MPO, with all sites empty
   */
  MPO(int length);
  /** Creates an MPO, with operator o[k] on site k+1. Only the pointers are stored, not a copy of the operators. */
  MPO(int length,const Operator** oplist);
  /** 
      Construct a new MPO from joining together a number of
      MPOs, in the order they are passed as arguments. The
      resulting MPO has as elements basic Operator objects,
      whose dimensions (left and right) are the product of all the
      original ones. The new Operator objects are stored in the MPO.
      \warning For some optimization (like remembering the
      structure of multiple operators and using it in the
      contractions) the method join should be used instead.
  */
  MPO(int nr,const MPO* oprs[]);

  virtual ~MPO();

  /** Discard all contents */
  void clear();
  /** Clear the MPO and restart it, as a fresh, empty MPO with 
      the given length */
  void initLength(int len);
  /** Length of the chain */
  int getLength() const {return length;}
  /** Retrieve a particular Operator (numbering of sites from 0 to
      N-1!) */
  const Operator& getOp(int pos) const;
  /** Return true iff the corresponding operator is empty */
  bool isEmpty(int pos) const {return Ops[pos]==0;}
  /** Delete a stored Operator (USE WITH CARE!)*/
  void deleteOp(int pos);

  /** Retrieve data for a particular Operator (sites from 0 to
      N-1!). Returns a copy to the stored data.  */
  const mwArray getOpData(int pos) const;

  /** Set the operator in position pos (running from 0 to length-1).
     If flag myop is set to one, the MPO assumes ownership of the
     Operator pointer and will delete the object upon its own
     destruction. */
  void setOp(int pos,const Operator* op,bool myop=0);

  /** Set the Operator in position pos to a copy of op,
      after complex conjugation. */
  void setConjugatedOp(int pos,const Operator* op);

  /** Creates an Operator in position pos (0 to length-1) so that it
     saves a COPY of the data in op (which can be given as Operator or
     mwArray), but with dimensions permuted as indicated by newdorder.
     If conjugate==true, the data are also complex conjugated.  Since
     the operator (data) in op must be copied to be rotated, this
     method creates a new Operator.
  */
  void setRotatedOp(int pos,const Operator* op,const shrt::Indices& neworder,
		    bool conjugate=0);
  void setRotatedOp(int pos,const mwArray& op,const shrt::Indices& neworder,
		    bool conjugate=0);

  /** Print contents to a stream */
  friend std::ostream& operator<<(std::ostream& os,const MPO& ops);

  // To be overloaded by inherited classes
  virtual MPO* getBasicRow() const;

  /** Save to a binary file */
  void saveMPO(const char* filename) const;
  /** Read an MPO from a binary file */
  void loadMPO(const char* filename);

  /** Export to a text file, to make it easier to check with older
      version (also more portable). */
  void exportMPOtext(const char* filename) const;
  /** Import from a text file, to make it easier to check with older
      version */
  void importMPOtext(const char* filename);

  /** Export in text mode, to be able to copy into Matlab directly. The optional argument precision can be used to specify the precision with which values are stored  */
  void exportForMatlab(const char* filename, int precision=0) const;

  /** Set the MPO from a vector of mwArray */
  void setOperatorArray(const std::vector<const mwArray*>& opers);

  /** Returns a vector with the (final) physical dimensions of 
     every site in the chain. */
  const std::vector<int> getDimensions() const;

  /** Fill in the argument with the original) physical dimensions of 
      every site in the chain. */
  const std::vector<int> getOriginalDimensions() const;

  //#ifndef TESTINGMPO
    //protected:
  //#endif
  MPO& operator=(const MPO& opr);


 private:
  /** Creates a copy. Not to be used */
 MPO(const MPO& opr):length(0),Ops(NULL),myown(1),myops(NULL){};

  /** Just the repeated code that marks a given operator as own */
  void markAsOwn(int pos);

 public:
  /** 
      From this MPO, construct another one which is equivalent
      to the folded version of the original. If the original length
      was odd, the central operator will be now in position 0 and
      length will be (L+1)/2, being L the original length. If L was
      even, new length will be L/2 with every operator a folded one.
      Returns the new MPO (with copies of the tensors) in \param<mpo>.
  */
  void fold(MPO& mpo) const;

  /** 
      Construct a new MPO from joining together a number of MPOs, in
      the order they are passed as arguments. The elements in the
      resulting row are not simple Operator but JoinedOperator
      objects.  Notice that the individual MPOs (actually, the
      Operators) are copied.  The order of the MPOs should be such
      taht the last in the list is the one acting first on the ket.
  */
  static void join(int nr,const MPO* oprs[],MPO& joined);
};


/** Utility function which, given a MPO (the product of physical
 * dimensions should be smaller than 2^10, or it will crash) computes
 * the full matrix in a d^L space, with indices ordered
 * (d0,d1,d2...)x(d0'd1'd2'...)
 */
void expandOper(const MPO& mpo,mwArray& oper);

// For compatibility with former name
typedef MPO OperatorRow;

#endif  // MPO
