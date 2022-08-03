/**
   \file lslugsl.cpp
   
   \author Stefan Kühn
   \date 25/10/2013
*/

#include "wrapper.h"
#include "mwArray.h"
#include "Indices.h"
#include <cstdio>
#include <vector>
//Libraries which are needed to use the GSL-solver
#include<gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>


using namespace std;
using namespace shrt;


void wrapper::lsgsl(const mwArray& A,const mwArray& B,mwArray& X){
   //All the variables needed to solve the system
   gsl_matrix_complex *Amat;
   gsl_vector_complex *bvec,*xvec;
   gsl_permutation *p;
   int signum;
   
   int M=A.getDimension(0);
   int N=A.getDimension(1);  
   complex_t temp;
   Amat = gsl_matrix_complex_alloc(M,N);
   xvec = gsl_vector_complex_alloc(N);
   bvec = gsl_vector_complex_alloc(M);
   p = gsl_permutation_alloc(N);
   
   for(int i=0; i<M; i++)
   {
     for(int j=0;j<N;j++)
     {
       temp=A.getElement(Indices(i,j));
       gsl_matrix_complex_set(Amat,i,j,gsl_complex_rect(temp.re,temp.im));
     }
   }
   if(M!=N)
   {
     cout << "Error: A is not a rectangular matrix, instead it has dimensions "<< A.getDimensions()<<", lsgsl is not able to solve the system" << endl;
     exit(212);
   }   
   
   M=B.getDimension(0);
   N=B.getDimension(1);
   
   if(N!=1)
   {
     cout << "Error: B is not a vector, instead it has dimensions "<< ", lsgsl is not able to solve the system" << B.getDimensions() << endl;
     exit(212);
   }
   
   for(int i=0; i<M; i++)
   {
     temp=B.getElement(Indices(i,0));
     gsl_vector_complex_set(bvec,i,gsl_complex_rect(temp.re,temp.im));
   }
   
   //Solve linear system
   gsl_linalg_complex_LU_decomp(Amat,p,&signum);
   gsl_linalg_complex_LU_solve (Amat,p,bvec,xvec);
   
   //Transfer result to X
   X=mwArray(Indices(M,1));
   gsl_complex temp2;
   for(int i=0; i<M; i++)
   {
     temp2 = gsl_vector_complex_get(xvec,i);
     X.setElement(GSL_REAL(temp2),GSL_IMAG(temp2),Indices(i,0));
   }
   //Free memory
   gsl_matrix_complex_free(Amat);
   gsl_vector_complex_free(bvec);
   gsl_vector_complex_free(xvec);
   gsl_permutation_free(p);
}
