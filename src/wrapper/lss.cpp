/**
   \file lsd.cpp
   
   \author Mari-Carmen Banuls
   \date 23/10/2013
*/

#include "wrapper.h"
#include "mwArray.h"
#include "Indices.h"
#include <cstdio>
#include <vector>

using namespace std;
using namespace shrt;

// Default precision for LSS cutoff
//#define PRECLSS 1E-15
// Machine precision
#define PRECLSS -1.

extern "C" { 
 
  /** Solution of LLS (linear least squares) problem A*x =b
   using SVD decomposition. */
  void zgelss_(int* m,int* n,int* nrhs,double* A,
	       int* lda,double* B,int* ldb,
	       double* S, double* rcond,
	       int* rank, 
	       double* work,int* ldwork,double* work2,
	       int* info);
}

void wrapper::lss(const mwArray& A,const mwArray& B,mwArray& X,
		  const double& tol){
  //cout<<"wrapper::lsd(A"<<A.getDimensions()<<", B"<<B.getDimensions()<<", tol"<<tol<<")"<<endl;
  // Prepare arguments to call the Lapack routine
  int M=A.getDimension(0);
  int N=A.getDimension(1);
  mwArray auxA(A); // fortran routine destroys this
  int nrhs=1;int lda=M;int ldb=max(M,N);
  double* result=B.copyComponents(); // copy the values in B
  //X=B;
  int rank=0;
  double rcond=tol==0.?PRECLSS:tol;
  double* S=new double[min(M,N)];
  int ldwork=6*max(M,N);
  double* work=new double[2*ldwork];
  double* work2=new double[2*5*min(M,N)];
  int info;
  zgelss_(&M,&N,&nrhs,auxA.components,&lda,result,&ldb,
	  S,&rcond,&rank,work,&ldwork,
	  work2,&info);
  if(info!=0){
    if(info<0)
      cout<<"Error: illegal input argument for zgelss at position "
	  <<-info<<endl;
    else
      cout<<"Error: failed to converge in zgelss, info="<<info<<endl;
    exit(212);
  }
  // Rearrange things
  X=mwArray(Indices(N,1),result);
  //contraction1=Hconjugate(result)*A*result;
  //contraction2=Hconjugate(result)*B;
  delete[] S;delete[] work;delete[] work2;
  delete[] result;
}
