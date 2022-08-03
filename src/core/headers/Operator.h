/**
   \file Operator.h 

   Definition of the basic class that provides an Operator to be part
   of a MPO, and act on MPS. 
   
   \author Mari-Carmen Banuls
   \date 16/03/2011
*/

#ifndef OPERATOR_H
#define OPERATOR_H

#include "Indices.h"
#include "wrapper.h"
#include "mwArray.h"
#include "Site.h"

/** A general operator (from an operator row) which can act on a 
    given site of a MPS chain. Particular cases may have specific
    ways of being contracted. Therefore this could be mainly the general
    interface, that should be extended appropriately. */

class Operator{
  
 protected:
  shrt::Indices dims;
  //int d_; // resulting physical dimension
  //int d2_; // original physical dimension (that of the vector we act on)
  //int Dl_;
  //int Dr_;
  //mwArray data__;
  const mwArray* data;
  bool mydata; // whether I must destroy it at the end

 public:
  /** Construct an empty operator of dimensions d x Dl x d2 x Dr (if d2 not specified, it is assumed to be d) */
  Operator(int d,int Dl,int Dr,int d2=0);
  /** Construct an empty operator of dimensions dims */
  Operator(const shrt::Indices& dims);
  /** Construct an operator, for a given array whose dimensions should
     be arranged as: du*dl*dd*dr. Saves an internal copy of the
     data.  */
  Operator(const mwArray& oper);
  /** Special case: Construct an operator, from a rotated array whose 
     dimensions should be rearranged (permuted) as neworder indicates,  
     so that at the end they are du*dl*dd*dr 
     If conjugate==1, the date is complex conjugated before saving it. */
  Operator(const mwArray& oper,const shrt::Indices& neworder,
	   bool conjugate=0);

  /** 
      Construct a new Operator from joining together a number of
      Operators, in the order they are passed as arguments. The
      last operator in the array should be the first one acting on the
      ket. The resulting dimensions (left and right) are the product 
      of all the original ones. For some optimization (like remembering the
      structure of multiple operators and using it in the
      contractions) a JoinedOperator should be used instead.
  */
  Operator(int nr,const Operator* oprs[]);

  virtual ~Operator();

  /** Return the identity operator for a particular physical dimension */
  static const Operator* identity(int d);

  /** Copy a certain array into the local data, replacing the existing
      one, if any. */
  // virtual void copyData(const mwArray& oper);

  /** If I do not want to copy the array, just store the pointer */
  void setData(const mwArray* oper);

  /** Return (reference to) the full mwArray representing this
      Operator. It has to be implemented by composite operators. */
  virtual const mwArray getFullData() const {return *data;}

  /** Return a copy of the dimension vector */
  const shrt::Indices getDimensions() const{return dims;}

  /** Conjugate all stored data. If the data was not mine, it is
      copied now, because I am not supposed to change some external
      array. */
  virtual void conjugateOp();

  /** Contract the right index of data with the left one of rightOp
      and reshape the result, so that the internal data now is the 
      product of both. 
      \warning This might change the physical dimensions!!
      Inherited classes with structure should implement this to 
      keep their properties.
  */
  virtual void contractRight(const Operator* rightOp);

  /** Permute indices in a given manner. It takes as argument
      \param<perm> a length four \class<Indices>, and should be
      overloaded by any specialized Operator composed of several
      pieces (as JoinedOperator or FoldedOperator)
      If the data was not mine, it is
      copied now, because I am not supposed to change some external
      array.
  */
  virtual void permuteOp(shrt::Indices perm);


  /** Constructs a new Operator saving a COPY of the data in oper, but 
     with dimensions permuted as indicated by neworder. If conjugate==true, 
     the data are also complex conjugated. The values of d, d2, Dl and Dr 
     are set according to the newly stored data. */
  void setRotatedData(const mwArray& oper,const shrt::Indices& neworder,
		      bool conjugate=0);

  /** Contract the operator with a left term, between bra and ket, 
      and return a new "left" term. */
  virtual void contractL(mwArray& result,const mwArray& termL,const Site& ket,const Site& bra,bool dagger=0) const ;
  /** Contract the operator with a right term, between bra and ket, 
      and return a new "right" term. */
  virtual void contractR(mwArray& result,const mwArray& termR,const Site& ket,const Site& bra,bool dagger=0) const ;

  /** Contract the operator+ operator with a left term, between bra and ket, 
      and return a new "left" term. If dagger==1, instead of oper+ oper, 
      oper oper+ is calculated. */
  virtual void contractL2(mwArray& result,const mwArray& termL,const Site& ket,const Site& bra,bool dagger=0) const ;

  /** Contract the operator+ operator with a right term, between bra and ket, 
      and return a new "right" term. */
  virtual void contractR2(mwArray& result,const mwArray& termR,const Site& ket,const Site& bra,bool dagger=0) const ;

  /** Contract operator with a ket, between left and right terms */
  virtual void contractMket(mwArray& result,const mwArray& termL,const Site& ket,const mwArray& termR,bool dagger=0) const ;

  /** Contract operator with a bra, between left and right terms */
  virtual void contractMbra(mwArray& result,const mwArray& termL,const Site& bra,const mwArray& termR,bool dagger=0) const ;

  /** Contract operator between left and right terms to produce a N matrix */
  virtual void contractN(mwArray& result,const mwArray& termL,const mwArray& termR,bool dagger=0) const ;

  /** Contract the operator+ operator with left and right terms, 
      and return a N matrix. If dagger==1, instead of oper+ oper, 
      oper oper+ is calculated. */
  virtual void contractN2(mwArray& result,const mwArray& termL,const mwArray& termR,bool dagger=0) const ;

  friend ostream& operator<<(ostream& os,const Operator& sigma);
  friend class MPO;

  int getd() const {return dims[0];}
  int getDl() const {return dims[1];}
  int getDr() const {return dims[3];}
  int getdorig() const {return dims[2];}

  /** For exchange with other programs: save to a text file */
  virtual void savetext(std::ofstream& outfile) const;
  virtual void loadtext(std::ifstream& outfile);

  /** Function that actually prints the Operator, and gets overloaded
      by daughter classes. */
  virtual void put(std::ostream& os) const;

#ifndef TESTINGOBJ
 protected:
#endif
  // Added for some special operators which do not have a single data array
  Operator();
  // Do I need this? Just in case...
  Operator(const Operator& op);

  void clear();
  virtual Operator& operator=(const Operator& op);
  /** Construct an operator, for a given array but ignore the dimensions of 
     it and set d, Dl and Dr (d2) as specified. This is for derived classes 
     that do not store the operator as a block, like UOperator. */
  Operator(const mwArray& oper,int d,int Dl,int Dr,int d2=0);
  /** If I do not want to copy the array, just store the pointer. This 
     version is to be used by derived classes, as the former method. */
  void setData(const mwArray* oper,int d,int Dl,int Dr);

  /** Let friends obtain a COPY of  my internal data */
  //virtual const mwArray getData() const {return *data;} // BY VALUE?!?!
  // Return the corresponding dimension of the actual stored mwArray
  virtual int getInnerDimension(int k) const;


#ifndef TESTINGOBJ
 private:
#endif
  /** Auxiliary functions */
  /** Contract one term from left, assuming operator oper
   * ket is a Site (from the ket), dimension [d,Dl,Dr]
   * bra is the corresponding A from bra (but NOT CONJUGATED yet!)
   * oper is the operator matrix (du x alphal x dd x alphar)
   * termL is the temporary result, from terms to the left 
   * (dimensions [beta,alpha,beta]) */
  void contractleftop(mwArray& result,const mwArray& termL,const Site& ket,
		      const mwArray& oper,const Site& bra) const;
  /** Contract one term from right, assuming operator oper
   * ket is a Site (from the ket), dimension [d,Dl,Dr]
   * bra is the corresponding A from bra (but NOT CONJUGATED yet!)
   * oper is the operator (d1 x alphaL x d2 x alphaR)
   * termR is the temporary result, from terms to the right 
   * (dimensions [beta,alphaR,beta] */
  void contractrightop(mwArray& result,const mwArray& termR,const Site& ket,
		      const mwArray& oper,const Site& bra) const;
  /** Return the effective operator oper+ oper (d2 x alphaL^2 x d2 x
      alphaR^2). If dagger==1, returns oper oper+ instead
  */
  const mwArray doubleop(bool dagger) const;

  void contractNop(mwArray& result,const mwArray& termL,const mwArray& oper,
		   const mwArray& termR) const;

};




#endif // OPERATOR_H
