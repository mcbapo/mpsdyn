/**
   \file lsd.cpp
   
   \author Mari-Carmen Banuls
   \date 16/03/2011
*/

#include "wrapper.h"
#include "mwArray.h"
#include "Indices.h"
#include <cstdio>
#include <vector>

//Library which I need in case to use the GSL-solver
#ifdef USING_GSL
#include<gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#endif

using namespace std;
using namespace shrt;

// Default precision for LSD cutoff
//#define PRECLSD 1E-15
// Machine precision
#define PRECLSD -1.

extern "C" { 
 
  /** Solution of LLS (linear least squares) problem A*x =b
   using SVD decomposition. */
  void zgelsd_(int* m,int* n,int* nrhs,double* A,
	       int* lda,double* B,int* ldb,
	       double* S, double* rcond,
	       int* rank, 
	       double* work,int* lwork,double* rwork,int* iwork,
	       int* info);
}

void wrapper::lsd(const mwArray& A,const mwArray& B,mwArray& X,
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
  double rcond=tol==0.?PRECLSD:tol;
  double* S=new double[min(M,N)];
  int lwork=6*max(M,N);
  double* work=new double[2*lwork];
  int smlsiz=30;
  int nlvl=(int)(1+log2(min(M,N)/(smlsiz+1)));
  int lrwork=10*N+2*N*smlsiz+8*N*nlvl+3*smlsiz+max((smlsiz*smlsiz),2*N+4);
  double* rwork=new double[lrwork];
  int* iwork=new int[min(M,N)*2*(3*nlvl+11)];
  int info;
  zgelsd_(&M,&N,&nrhs,auxA.components,&lda,result,&ldb,
	  S,&rcond,&rank,work,&lwork,
	  rwork,iwork,&info);
  if(info!=0){
    if(info<0)
      cout<<"Error: illegal input argument for zgelsd at position "
	  <<-info<<endl;
    else
      cout<<"Error: failed to converge in zgelsd, info="<<info<<endl;
    exit(212);
  }
  //Option to monitor condition number (might be helpful to detect stability issues)
  /*double minimum;
  double maximum;
  for(int k=0;k<min(M,N); k++)
  {
    if(k==0)
    {
      minimum=S[0];
      maximum=S[0];
    }
    else if(S[k]<minimum)
      minimum=S[k];
    else if(S[k]>maximum)
      maximum=S[k];
      
  }
  //Possibly ill conditioned?
  if(maximum/minimum>2.0)
    cout << "Condition number="<< maximum/minimum << endl;*/
  
  
  // Rearrange things
  X=mwArray(Indices(N,1),result);
  //contraction1=Hconjugate(result)*A*result;
  //contraction2=Hconjugate(result)*B;
  delete[] S;delete[] work;delete[] rwork;delete[] iwork;
  delete[] result;
}
