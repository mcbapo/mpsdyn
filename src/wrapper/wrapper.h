#ifndef WRAPPER_H
#define WRAPPER_H

#include <vector>
#include <iostream>
#include <fstream>
#include "complex.h"
#include "Indices.h"


#ifdef USING_PRIMME
#ifndef PRIMMEv2
extern "C" {
#endif
#include "primme.h"
#ifndef PRIMMEv2
}
#endif
#endif

/** 
\file Some functions that encapsulate the calls to LAPACK or BLAS
routines, with a user friendly interface.

Also utility functions to manipulate dimension vectors, among others. 
  */

class mwArray;
class Multiplier;


/** Namespace wrapper hides the direct calls to Lapack, blas and
    Arpack routines under a more friendly interface.  */

namespace wrapper{

  /** Singular Value decomposition of matrix A, returns U, S, V+ such
      that \f$A=U\cdot S\cdot V^{\dagger}\f$. If flag full is set to
      true, full unitaries are returned.
 */
  void svd(const mwArray& A,mwArray& U,mwArray& S,mwArray& Vdagger,bool full=false);

  /** Computes the singular Value decomposition of matrix A, returns
      U, S, V+ such that \f$A=U\cdot S\cdot V^{\dagger}\f$. Then
      analyzes how many singular values are bigger than the provided
      tol, and returns this number in nr. If the original matrix was
      MxN, the sizes of the returned arrays will be U(Mxnr), S(nrxnr)
      and V(nrxnr). 
      IMPORTANT: If tol is to be used to cut the small values, set nr to 0.
      Otherwise, nr values are kept at most.
      \todo Call the lapack routine once before, to assess optimal size!
  */
  void svd(const mwArray& A,const double& tol,int& nr,
	   mwArray& U,mwArray& S,mwArray& Vdagger,bool full=false);

  /** Product of two matrices, returning the result in C  */
  void product(const mwArray& A,const mwArray& B,mwArray& C);

  /** Full eigenvalue decomposition of matrix A. On exit, D contains
      the eigenvalues, and, if eigv==true, U contains the
      eigenvectors, as columns. */
  void eig(const mwArray& A,std::vector<complex_t>& D,mwArray& U,bool eigv=0);

  /** Eigenvalue calculation for general (non-Hermitian) matrices. */
  void eigNH(const mwArray& A,std::vector<complex_t>& Dval,mwArray& U,
	     bool eigv=0,bool left=0);

  /** Use iterative method to compute k eigenvalues and eigenvectors of
      the matrix A. On exit, D contains the eigenvalues, and U contains the
      corresponding eigenvectors, as columns, if \param<eigv> was true.
      \param <k> decides how many eigenvalues are computed
      \param <which> should be a string specifying the criterion:
      The possible options are
			   "LM": largest magnitude
			   "SM": smallest magnitude
			   "LA": largest real component
			   "SA": smallest real compoent
			   "LI": largest imaginary component
			   "SI": smallest imaginary component
      If U is not empty, the first column is taken as the initial vector.
 */
  int eigs(const mwArray& A,const int k,const char* which,
	    std::vector<complex_t>& D,mwArray& U,bool eigv=0,
	    double tol=0.);

  /** Use iterative method to compute k eigenvalues and eigenvectors of
      the matrix A. On exit, D contains the eigenvalues, and U contains the
      corresponding eigenvectors, as columns. 
      \see Multiplier
  */
  int eigs(Multiplier& multi,const int k,const char* which,
	    std::vector<complex_t>& D,mwArray& U,bool eigv=0,double tol=0.);


#ifdef USING_PRIMME
  #ifdef PRIMMEv2
  #define DYNAMIC_ver PRIMME_DYNAMIC
  #define Arnoldi_ver PRIMME_Arnoldi
  #define JDQR_ver PRIMME_JDQMR_ETol
  #else
  #define DYNAMIC_ver DYNAMIC
  #define Arnoldi_ver Arnoldi
  #define JDQR_ver JDQR
  #endif
  /** Compute k eigenpairs of the matrix using PRIMME, an 
      alternative eigensolver to ARPACK. 
      All parameters are as in eigs, except the following:
      \param <which> indicates which eigenvalues to compute. Should be
      one of the following:
               primme_smallest (default): Smallest algebraic
	       primme_largest: Largest algebraic
	       NOT YET SUPPORTED (require set of shifts)
	       primme_closest_geq: Closest to (>=) set of shifts
	       primme_closest_leq: Closest to (<=) set of shifts
	       primme_closest_abs: Closest (in abs value) to set of shifts  
  */
  int eigs_primme(const mwArray& A,const int k,primme_target which,
		  std::vector<complex_t>& D,mwArray& U,double tol=0.,
		  primme_preset_method method=DYNAMIC_ver);
  int eigs_primme(Multiplier& multi,const int nev,
		  primme_target which,std::vector<complex_t>& Dval,mwArray& U,
		  double tol=0.,
		  primme_preset_method method=DYNAMIC_ver,
		  int numShifts=0,double* targetShifts=NULL);
  int eigsgen_primme(Multiplier& multiA,Multiplier& multiB,const int nev,
		     primme_target which,std::vector<complex_t>& Dval,
		     mwArray& U,
		     double tol=0.,
		     primme_preset_method method=DYNAMIC_ver,
		     int numShifts=0,double* targetShifts=NULL);
  int call_eigs_primme(primme_params* primme,std::vector<complex_t>& Dval,
		       mwArray& U);
#endif

  /** Solve a system of linear equations by a linear least squares
      routine.  A*x=B is solved here using SVD decomposition.  If
      argument tol>0 is specified, it is used as a (relative value!)
      cutoff to decide on zero singular values of A. */
  void lsd(const mwArray& A,const mwArray& B,mwArray& X,
	   const double& tol=0.);

  /** Solve a system of linear equations by a linear least squares
      routine.  A*x=B is solved here using SVD decomposition.  If
      argument tol>0 is specified, it is used as a (relative value!)
      cutoff to decide on zero singular values of A. 
      Introduced to void crash of zgelsd for large matrices.
  */
  void lss(const mwArray& A,const mwArray& B,mwArray& X,
	   const double& tol=0.);

  /** Compute the LU decomposition of a given matrix, A.
      The result is returned in the arguments L (lower triangular
      with 1 in the diagonal), U (upper triangular) and P (permutation matrix), 
      such that A=P*L*U */
  void lu(const mwArray& A,mwArray& L,mwArray& U,mwArray& P);

  /** Solve a system of linear equations A*x=B, for a general square matrix A,
      using the LU decomposition. */
  void lslu(const mwArray& A,const mwArray& B,mwArray& X);
  
  /** Solve a system of linear equations A*x=B, for a general square matrix A,
      using the LU decomposition. This version employs the solver from the 
      Gnu Scientific Library. */
  void lsgsl(const mwArray& A,const mwArray& B,mwArray& X);


  /** QR decomposition of a general matrix A=Q*R If A has dimensions
      d1 x d2, the function returns a unitary (or isometry) Q of
      dimensions d1 x min(d1,d2) and an upper triangular (trapezoidal)
      R min(d1,d2) x d2 */
  void qr(const mwArray& A,mwArray& Q,mwArray& R);

  void lq(const mwArray& A,mwArray& L,mwArray& Q);

  /** Compute the matrix exponential by diagonalizing a full matrix
      and exponentiating the diagonal of eigenvalues:
      \f$ exp(A) = U exp(D) U^{\dagger}\f$. 
      The matrix should be diagonalizable. The coefficient coeff
      can contain a factor to put in the exponent.
*/
  void expm(const mwArray& A, mwArray& result,complex_t coeff=ONE_c);

/** Exponential of a square matrix, not necessarily diagonalizable, by
 * computing the Taylor series up to a certain order. */
  void expmTaylor(const mwArray& A, mwArray& result,int order,
		  complex_t coeff=ONE_c);
}


/** Auxiliary function to print a vector of indices */
std::ostream& operator<<(std::ostream& os,const std::vector<int>& dims);

/** Auxiliary function to print a vector of doubles */
std::ostream& operator<<(std::ostream& os,const std::vector<double>& dvals);

/** Auxiliary function to print a vector of complex values */
std::ostream& operator<<(std::ostream& os,const std::vector<complex_t>& cvals);

/** Auxiliary function to compare two vectors of indices and return
    true onl if every index in the first vector is smaller than the
    corresponding index in the second vector. */
//bool operator<(const shrt::Indices& dimsA,const shrt::Indices& dimsB);


#endif // WRAPPER_H
