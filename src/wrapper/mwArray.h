#ifndef MWARRAY_H
#define MWARRAY_H

#include "complex.h"
#include "wrapper.h"
#include "Indices.h"
#include <fstream>
#include <cstring> // for memcpy
#include <cstdlib> // for random, should be sth better
#include <cmath>
#include <vector>
#include "stdarg.h" // for va_arg


/** 

    This is an intermediate class, providing mwArray interface (as Matlab) to
    tensors that contain my numerical data.
    These tensors are always complex (if data are real, zeros are stored).
    
*/

class mwArray {
  int rank; // nr of dimensions
  shrt::Indices dims; // actual vector of dimensions
  int nrcomponents; // total number of components 
  shrt::Indices cumdims; // auxiliary for getIndices
  bool computed;
 protected:
  double* components;
  bool mine;
 public:
  static const mwArray empty;
  /** Create an empty array */
  mwArray();
  /** Particular case: create a single component array */
  mwArray(complex_t c);
  /** Creates a new tensor with given number of dimensions and
      corresponding number of (complex) components. Initialized to zero. */
  mwArray(const shrt::Indices& dimensions);
  /** The same constructor can receive the dimensions as a simple
      array of integers */
  mwArray(int rank,const int* dimensions);
  
  //  /** Same, just passing a list of integers */
  //mwArray(int rank,...);

  /** Create a new tensor, filling it in with the given buffer.
      \param <dimensions> specifies the dimensions of the tensor
      \param <real> (default false) should be true if the data buffer
      contains only the real parts of the elements (complex parts are zero)
      \param <data> contains the sequence of element values (if
      real==false real and imaginary double per element)
  */
  mwArray(const shrt::Indices& dimensions,const double* data,bool real=false);
  /** Create a new tensor, filling it in with the given buffer
      of complex numbers.
  */
  mwArray(const shrt::Indices& dimensions,const complex_t* data);
  /** Create from a file */
  mwArray(std::ifstream& data);
  /** Copy constructor: performs a deep copy of the argument. */
  mwArray(const mwArray& orig);
  /** Special case: transform a vector into a mwArray, with the
   * dimensions of a vector (default is column order). This version
   * takes real components. */
  mwArray(const std::vector<double>& revals,bool column=true);
  /** Same as before, but taking complex components.
   */
  mwArray(const std::vector<complex_t>& cvals,bool column=true);
  /** Standard destructor */
  ~mwArray();
  /** Reduce it to empty mwArray */
  void clear();

  /** Load the mwArray from a binary file. */
  void load(std::ifstream& file);
  /** Save the mwArray to a binary file. */
  void save(std::ofstream& file) const;
  /** Save in text mode for exchange with other programs.  
      \param real should be true if the data are real, so that only one component
      is saved (but then it must be read with the flag real set to
      true!) 
  */
  void savetext(std::ofstream& data,bool real=false) const;
  /** Load in text mode for exchange with other programs */
  void loadtext(std::ifstream& data,bool real=false);
  /** Returns true iff the number of dimensions is two. */
  bool isMatrix() const {return (!isVector())&&(!isScalar())&&rank==2;}
  /** Returns true iff the number of dimensions is two and they are equal. */
  bool isSquare() const {return (isMatrix())&&(dims[0]==dims[1]);}
  /** Returns true iff the array is empty. */
  bool isEmpty() const {return rank==0;}
  /** Returns true iff the array is 2-dimensional and diagonal. */
  bool isDiagonal() const;
  /** Returns true iff the object has only one component. */
  bool isScalar() const {return nrcomponents==1;}
  /** Returns true iff the number of dimensions is one, or if it is
      two and one of them has size 1. */
  bool isVector() const {return (nrcomponents>1)&&
      (rank==1||(rank==2&&(dims[0]==1||dims[1]==1)));}
  /** Returns true iff the tensor corresponds to an Hermitian matrix
      (or real scalar), i.e. iff \f$\frac{A-A^{\dagger}}{2}=0\f$,
      according to the syntax of isNull. */
  bool isHermitian(double tol=0) const;
  /** Returns true if all elements are zero. If optional argument
      tol>0 is given, elements smaller than tol in absolute value are
      considered to be 0. */
  bool isNull(double tol=0) const;

  /** Return the largest element (in absolute value) and its position*/
  void getMaxElement(complex_t& value,shrt::Indices& indices) const;

  /** Used to manually set the value of a single component.
      (revalue,imvalue) are the real and imaginary parts of the
      component with indices given by the vector.
  */
  void setElement(double revalue,double imvalue,const shrt::Indices& indices);

  /** Used to manually set the value of a single complex component.
      with indices given by the vector.
  */
  void setElement(complex_t value,const shrt::Indices& indices);
  /** Idem from the linear index */
  void setElement(complex_t value,int linidx);

  /** Set the real part of all elements in the array from a list of
   doubles.  \param <nr> must be exactly the number of components in
   the array, and this is the number of elements read afterwards from
   the list of arguments.  The list of values are saved COLUMN MAJOR
   */
  void setRealData(int nr,...);
  /** Set the imaginary part of all elements in the array from a list
   of doubles.  \param <nr> must be exactly the number of components
   in the array, and this is the number of elements read afterwards
   from the list of arguments.  The list of values are saved COLUMN
   MAJOR */
  void setImagData(int nr,...);

  /** Set the data from an array of nr complex values saved in COLUMN
      MAJOR order */
  void setData(int nr,const complex_t data[]);

  /** 
      Returns, in revalue and imvalue, the real and imaginary part of 
      the element corresponding to the given indices.
   */
  inline void getElement(double& revalue,double& imvalue,
			 const shrt::Indices& indices) const;
  /** 
      Returns, in a complex_t, the real and imaginary part of 
      the element corresponding to the given indices.
   */
  inline complex_t getElement(const shrt::Indices& indices) const;
  /** Idem from a position */
  inline complex_t getElement(int linidx) const;

  /** Fill in the vector with the array dimensions */
  void getDimensions(shrt::Indices& dimensions) const{dimensions=dims;}

  /** Return a copy of the dimensions */
  shrt::Indices getDimensions() const{return dims;}

  /** Return the number of dimensions */
  int getRank() const {return rank;}

  /** Return the value of the k-th dimension of this array 
      (for k in [0...rank]) */
  inline int getDimension(int k) const;

  /** Return the number of components */
  int getNrComponents() const {return nrcomponents;}

  /** Redefine the dimensions, keeping the total number of elements
      constant. It only affects the way of accessing elements, but not
      the actual data.  One negative index can be included in the list
      of new dimensions, to provide the syntax of [] in Matlab,
      i.e. that index will have the required dimension to be
      compatible with the size of the array.
  */
  void reshape(const shrt::Indices& newdims);

  /** Reorder the dimensions of the matrix, as indicated by neworder.
      This requires copying the data in the new order. 
      \param <neworder> contains the new order of indices, from 1 to rank. 
      If \param<conj>==true, the result is conjugated at the same time.
      \todo Check cache usage. Check Ding's Vacancy tracking algorithm
      to optimize memory usage.  */
  void permute(const shrt::Indices& neworder,bool conj=false);
  /** Particular optimized case, for a transposition on a matrix */
  void transpose(bool conj=false);
  //void permute2(const std::vector<int>& neworder);

  /** Substitute the data by its complex conjugate */
  void conjugate();

  /** Substitute the data by its Hermitian conjugate (only for rank two arrays) */
  void Hconjugate();

  /** Change the sign of the current data, avoiding a copy */
  void changeSign();

  /** Fill in the mwArray with random values,discarding the previous
      buffer, if present. */
  void fillRandom();

  /** Fill with zeroes */
  void fillWithZero();

  /** Fill with (real) ones */
  void fillWithOne();

  /** Create a subarray by fixing some of the indices in the current
      one. The argument, ind, must have the same length as the number
      of dimensions of the array. The elements which are <0 mean that
      dimension is not touched, while any other value fixes the index
      corresponding to that dimension. A new mwArray is returned (as
      copy: could be done more efficiently) with the dimensions of the
      fixed indices collapsed.  */
  const mwArray subArray(const shrt::Indices& ind) const;

  /** 
      Change the dimensions to the ones in the argument. The vector
must have the same length as the current rank, and individual
dimensions are either expanded (padded with zeroes, or with random
numbers between 0 and noise) or reduced (first element kept, rest
discarded).  */
  void resize(const shrt::Indices& newdims,double noise=0.);

  /** Assignment operator: performs a deep copy of the rhs */
  mwArray& operator=(const mwArray& right);

  /** Replace the data by the result of multiplying with the argument
      on the right. If the argument has a single component, the
      current dimensions are kept, irrespective of right having a
      particular specification of (singleton) dimensions. */
  void multiplyRight(const mwArray& right);

  /** Replace the data by the result of multiplying with the argument
      on the left. If the argument has a single component, the
      current dimensions are kept, irrespective of left having a
      particular specification of (singleton) dimensions. */

  void multiplyLeft(const mwArray& left);


  /** Compute the trace of the (complex) rank 2 array, if square. */
  void trace(double& re,double& im) const;

  /** Return the trace in a complex value */
  complex_t trace() const;

  /** Returns the diagonal, if the tensor is a matrix, in the argument diag */
  void getDiagonal(mwArray& diag) const;
  /** Same, returns as a vector */
  void getDiagonal(std::vector<complex_t>& diag) const;

  /* Returns the real part of the diagonal as a (real) vector */
  void getRealDiagonal(std::vector<double>& diag) const;


  /** Contract leg nr k1 of this Tensor with leg nr k2 of the *
  argument B (k1 and k2 starting on 1). If rank of A was m and that of
  B, n, after returning, the tensor contains the result of the
  contraction with dimensions in the following order:
  [d(1)...d(k1-1),d(k1+1)...d(m),p(0)...p(k2-1),p(k2+1)...p(n)], where
  [d1...d(m)] and [p1...p(n)] where the original dimensions of A
  and B.  So far, no optimization is done regarding the order of the
  factors in the contraction, to make best use of matrix
  multiplication.  */
  void contractLeg(int k1,mwArray& B,int k2);
// As the previous one, but taking a constant mwArray as argument and
// making a temporary copy (less efficient) 
  void contractLeg(int k1,const mwArray& B,int k2);
  /* Same as the previous two, but returns the result in a new
   * tensor, and does not alter the */
  friend const mwArray contractLeg(const mwArray& A,int k1,
				   const mwArray& B,int k2);

/** Contract a list of legs */
  void contractLegs(const shrt::Indices& k1,mwArray& B,const shrt::Indices& k2);
  void contractLegs(const shrt::Indices& k1,const mwArray& B,
		    const shrt::Indices& k2);
  friend const mwArray contractLegs(const mwArray& A,const shrt::Indices& k1,
				    const mwArray& B,const shrt::Indices& k2);


  /** Comparison operator. Returns true iff arguments are identical,
      elementwise. */
  friend bool operator==(const mwArray& left,const mwArray& right);

  /** Product operator, contracting last dimension of argument left
      with first dimension of argument right.  */
  friend const mwArray operator*(const mwArray& left,const mwArray& right);

  /** Addition operator */
  friend const mwArray operator+(const mwArray& left,const mwArray& right);

  /** Substraction operator */
  friend const mwArray operator-(const mwArray& left,const mwArray& right);

  /** Multiplication by a (complex) scalar. I assume that it comes as
      an array with two components, real and imaginary parts. */
  friend const mwArray operator*(const double* alpha,const mwArray& mwA);

  /** Multiplication by a (complex) scalar. Includes the real and
      purely imaginary cases.*/
  friend const mwArray operator*(const complex_t& alpha,const mwArray& mwA);

  /** Multiplication by a real scalar.*/
  friend const mwArray operator*(const double alpha,const mwArray& mwA);

  /** Elementwise square root of complex elements.  \warning
      {Multivalued result not taken into account. Trivial choice of
      phase} */
  friend const mwArray sqrt(const mwArray& mwA);

  /** Invert a square diagonal matrix by inverting the diagonal
      elements one by one. If a tolerance is specified, the elements
      such that \f$z/trace(A)<tol\f$ are treated as zero, and the
      pseudoinverse is computed. */
  friend const mwArray invertDiag(const mwArray& mwA,int& nr);
  friend const mwArray invertDiag(const mwArray& mwA,int& nr,
				  const double tol);

  /** Construct a real diagonal matrix of dimensions nrxnr
      from an array of scalars */
  friend const mwArray realDiag(int nr,const double* values);
  friend const mwArray realDiag(int nr,const std::vector<double>& values);

  /** Construct a diagonal matrix of size nrxnr 
      from an array of complex scalars, given as pairs re,im */
  friend const mwArray diag(int nr,const double* values);

  /** Construct a diagonal matrix of size nrxnr 
      from an array of complex scalars, given as complex_t */
  friend const mwArray diag(int nr,const complex_t* values);
  friend const mwArray diag(const std::vector<complex_t>& values);
  /** Accept also a vector in mwArray form (only one non-trivial dimension) */
  friend const mwArray diag(const mwArray& values);

  /** Instead of reshaping in place, we can return a reshaped copy of
      the array */
  friend const mwArray reshape(const mwArray& A,const shrt::Indices& newdims);

  /** Instead of resizing in place, we can return a resized copy of
      the array */
  friend const mwArray resize(const mwArray& A,const shrt::Indices& newdims);

  /** Instead of reshaping in place, we can return a permuted copy of
      the array */
  friend const mwArray permute(const mwArray& A,const shrt::Indices& neworder);
  /** Instead of reshaping in place, we can return a conjugated copy of
      the array */
  friend const mwArray conjugate(const mwArray& A);
  /** Return the Hermitian conjugate of the matrix */
  friend const mwArray Hconjugate(const mwArray& A);

  /** Supported only for vectors and scalars: returns the absolute
   value or the Euclidean norm, respectively. TODO: implement for
   matrices */
  friend double norm(const mwArray& A);

  /** (Pseudo-)Invert a square matrix by diagonalizing and inverting the
eigenvalues above the specified tolerance.  
\todo Improve, using zgtref, zgetri */
  friend const mwArray inverse(const mwArray& mwA,int& nr);

  friend const mwArray inverse(const mwArray& mwA,int& nr,
			       const double tol);
  /** Utility to get the tensor product of two matrices, of any
      dimensions. Notice that, if arguments habe dimensions dA1xdA2
      and dB1xdB2, the result will be a rank-2 tensor with dimensions
      (dA1*dB1,dA2*dB2), and simply reshaping (without permutation)
      will produce the leg ordering (dA1,dB1,dA2,dB2)
 */ 
  friend void constructOperatorProduct(mwArray& result,const mwArray& opA,
				       const mwArray& opB);

  /** Operator to print some information about a mwArray. Mainly
      intended for debugging purposes, may or may not write the whole
      content of the buffer. 
  */
  friend std::ostream& operator<<(std::ostream& os,const mwArray& mwA);

  /** Idem, for Matlab friendly output, as we would input data in Matlab GUI, using name as the name of the array */
  friend std::ostream& putForMatlab(std::ostream& os,const mwArray& mwA,const char* name);
  friend std::ostream& putForMatlab(std::ostream& os,const mwArray& mwA,const char* name,int precision);


  /** 
      Return a pointer to the actual data buffer. It won't work on a
      const mwArray&, therefore there is a method to make a copy of
      the buffer without actually changing it.
   */
  const double* getComponents() const{return components;}

  /** Given two vectors, take the scalar product and return a scalar
      (in this case, as a complex). This method ignores
      the order in which the dimensions of each of the vectors are
      specified, and always contracts the long dimension. To compute
      something like \f$v \cdot v^T \f$ giving a matrix, use outerproduct. */
  friend const complex_t scalarproduct(const mwArray& A,const mwArray& B);
  /** Computes the outer product of two vectors, to produce the matrix
      \f$A \cdot B^T \f$, of dimensions dimA x dimB. */
  friend const mwArray outerproduct(const mwArray& A,const mwArray& B);

  /** Computes the tensor product of two matrices, of dimensions
      (da1,da2) and (db1,db2), to give a result (da1*db1,da2*db2). */
  friend const mwArray kron(const mwArray& A,const mwArray& B);
  /** A more efficient version: allows for reshape of the matrices */
  friend const mwArray kron(mwArray& A,mwArray& B);

  /** Hadamard product, or element-wise multiplication, of two arrays,
      which must have identical dimensions */
  friend const mwArray hadamard(const mwArray& A,const mwArray& B);


#ifndef TESTINGMWA
 protected:
#endif

  /** Gives back a deep copy of the data buffer, without changing the
      current object. Useful for copy construction and similar. */
  double* copyComponents() const;
  /** Get the address of a given element of the array. Should be
      called with as many integer arguments as the number of
      dimensions.  \todo Enforce this behavior.  */
  double* element(int d1,...);
  /** To be used by things like eigs, which need to set the inner
      pointer to a given address. */
  void setPointer(const shrt::Indices& dims,double* ptr);
  bool isOwned() const{return mine;}
  /** Check if a list of indices corresponds to a valid element in the
      array */
  bool checkIndices(const shrt::Indices& elem) const;

 private:
  /** Construct as result of a certain operation */
  mwArray(const mwArray& left,const mwArray& right,const char oper);
  /** Construct from product with scalar */
  mwArray(const double alpha[2],const mwArray& mwA);
  /** Construct from product with scalar */
  mwArray(const complex_t& alpha,const mwArray& mwA);
  /** Construct from product with real scalar */
  mwArray(const double alpha,const mwArray& mwA);

  /** Linear index that corresponds to a vector of indices. 
      This returns the position of the element in an array, but since 
      the elements are complex, the actual positions will be 
      2*idx and 2*idx+1 for the real and imaginary parts.
  */
  int linearIdx(const shrt::Indices& indices) const;

  /** Return the sequence of indices in the multidimensional array 
      that correspond to the linear index idx. */
  void getIndices(int idx,shrt::Indices& result) const;

  /** In general, given a vector of dimensions, compute the linear index. 
      \todo Should be out of here. Maybe a global function. */
  int linearIdx(const shrt::Indices& indices,
		const shrt::Indices& dimensions) const;

  friend void wrapper::svd(const mwArray& A,const double& tol,int& nr,
			   mwArray& U,mwArray& S,mwArray& Vdagger,bool full);

  friend void wrapper::product(const mwArray& A,const mwArray& B,mwArray& C);

  friend void wrapper::eig(const mwArray& A,std::vector<complex_t>& D,
			   mwArray& U,bool eigv);
  friend void wrapper::eigNH(const mwArray& A,std::vector<complex_t>& D,
			     mwArray& U,bool eigv,bool left);
  friend int wrapper::eigs(Multiplier& multi,const int k,
		   const char* which,std::vector<complex_t>& Dval,mwArray& U,
			    bool eigv,double tol);
  friend void wrapper::lsd(const mwArray& A,const mwArray& B,mwArray& X,
	   const double& tol);
  friend void wrapper::lss(const mwArray& A,const mwArray& B,mwArray& X,
	   const double& tol);
  friend void wrapper::lu(const mwArray& A,mwArray& L,mwArray& U,mwArray& P);
  friend void wrapper::lslu(const mwArray& A,const mwArray& B,mwArray& X);
  friend void wrapper::qr(const mwArray& A,mwArray& Q,mwArray& R);
  friend void wrapper::lq(const mwArray& A,mwArray& L,mwArray& Q);

 private:

  /** Return the sequence of indices in the multidimensional array
      that correspond to the linear index idx. For internal use: drops
      checking of dimensions and assumes vector is the right size.  */
  void _getindices_(int idx,shrt::Indices& result) const;
  void cumDims();
  int _linearidx_(const shrt::Indices& indices) const;
  int _linearidx_(const shrt::Indices& indices,
		  const shrt::Indices& dims_) const;

  friend class Transposer;
  friend class Multiplier;
};

// Inline functions

int mwArray::getDimension(int k) const{
  if(k<0){
    std::cout<<"Invalid dimension number("<<k<<") for mwArray of rank "
	     <<rank<<std::endl;
    exit(212);
  }
  if(k>=rank) return 1;
  return dims[k];
}

void mwArray::getElement(double& revalue,double& imvalue,
			 const shrt::Indices& indices) const{
  int linidx=linearIdx(indices);
  //cout<<"mwArray:getElement("<<indices<<")->"<<linidx<<endl;
  revalue=components[2*linidx];
  imvalue=components[2*linidx+1];
}

complex_t mwArray::getElement(const shrt::Indices& indices) const{
  complex_t c;
  getElement(c.re,c.im,indices);
  return c;
}

complex_t mwArray::getElement(int linidx) const{
  return (complex_t){components[2*linidx],components[2*linidx+1]};
}

/** Useful global functions */

const mwArray realDiag(int nr,const double* values);
const mwArray diag(int nr,const double* values);
const mwArray diag(const std::vector<complex_t>& values);
const mwArray diag(int nr,const complex_t* values);

/** Construct the identity matrix NxN */
const mwArray identityMatrix(int N);

//const mwArray& mwArray::empty=mwArray();

#endif // MWARRAY_H
