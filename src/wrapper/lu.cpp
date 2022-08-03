#include "wrapper.h"
#include "mwArray.h"
#include <cstdio>
#include <vector>

using namespace std;
using namespace shrt;

extern "C" {
  // LU decomposition for general matrices  
  void zgetrf_(int* m,int* n,
	       double* A,int* lda,
	       int* ipiv,int* info);
 // LS A*x=B for general MxM matrices  
  void zgetrs_(char* trans,int* m,int* nrhs,
	       double* A,int* lda,
	       int* ipiv,double* B,int* ldb,
	       int* info);
}

void wrapper::lu(const mwArray& A,mwArray& L,mwArray& U,mwArray& P){
  // Particular case: a scalar, save the call
  if(A.isScalar()){
    complex_t elem=A.getElement(vector<int>(A.getRank(),0)); // single element
    L=mwArray(ONE_c);
    U=A;
    P=mwArray(ONE_c);
    return;
  }
  if(A.getRank()!=2){
    cout<<"Error: Cannot do LU on an array if not rank 2, and here dims="
	<<A.getDimensions()<<endl;
    exit(212);
  }
  // Arguments for the fortran routine
  int info;
  int dimL=A.getDimension(0);int dimR=A.getDimension(1);
  U=A; // copied: replaced when coming back
  // Space for temporary computations (pivoting)
  int lwork=min(dimL,dimR);
  int* ipiv=new int[lwork];

  zgetrf_(&dimL,&dimR,U.components,&dimL,ipiv,&info);
  // Rearrange results
  L=U; 
  if(dimL<dimR) L.resize(Indices(dimL,dimL));
  if(dimL>dimR) U.resize(Indices(dimR,dimR));
  // Now the diagonal of L is 1 and the upper triang is zero
  for(int k=0;k<min(dimL,dimR);k++){
    L.setElement(ONE_c,Indices(k,k));
    for(int j=k+1;j<max(dimL,dimR);j++){
      if(j<L.getDimension(1))
	L.setElement(ZERO_c,Indices(k,j));
      if(j<U.getDimension(0))
	U.setElement(ZERO_c,Indices(j,k));
    }
  }

  // what about the permutation? According to Lapack:
  //   A = P * L * U
  //IPIV:  The pivot indices; for 1 <= i <= min(M,N), row i of the matrix was
  //	  interchanged with row	IPIV(i).
  /// But it seems to come as a row of transpostitions=>I transform it to a permutation of indices 1:dimL
  vector<int> tmp;for(int k=0;k<dimL;k++) tmp.push_back(k);
  cout<<"IPIV contains [";
  for(int k=0;k<lwork;k++){
    int trgt=ipiv[k]-1; // from 0
    cout<<trgt+1<<",";
    if(trgt!=k){
      int aux=tmp[k];
      tmp[k]=tmp[trgt];
      tmp[trgt]=aux;
    }
  }  
  cout<<"];\n After processing IPIV, the permutation vector looks like "<<tmp<<endl;
  P=mwArray(Indices(dimL,dimL));
  for(int k=0;k<dimL;k++)
    P.setElement(ONE_c,Indices(tmp[k],k));

  delete[] ipiv;
}


void wrapper::lslu(const mwArray& A,const mwArray& B,mwArray& X){
  // First, check that A is not a scalar, which is a particular case
  if(A.isScalar()){
    complex_t elem=A.getElement(vector<int>(A.getRank(),0)); // single element
    if(elem==ZERO_c){
      cout<<"Error, trying to invert zero in lslu!"<<endl;
      exit(212);
    }
    X=(ONE_c/elem)*B;
    return;
  }
  if(!A.isSquare()){
    cout<<"Error: Cannot do LSLU on a matrix that is not square. "
	<<A.getDimensions()<<endl;
    cout<<"Try lsd instead. "<<endl;
    exit(212);
  }
  if(B.getDimension(0)!=A.getDimension(0)){
    cout<<"Error, trying to solve linear system A x=B for A"<<A.getDimensions()
	<<" and B"<<B.getDimensions()<<endl;
    exit(212);
  }

  // I need to call the zgetrs routine, which requires the zgetrf factorization

  // Arguments for the fortran routine zgetrs
  int info;
  int dimL=A.getDimension(0);int dimR=A.getDimension(1);
  mwArray U=A; // copied: replaced when coming back
  // Space for temporary computations (pivoting)
  int lwork=min(dimL,dimR);
  int* ipiv=new int[lwork];
  zgetrf_(&dimL,&dimR,U.components,&dimL,ipiv,&info);
  if(info<0){
    cout<<"Problem solving lslu, after zgetrf info="<<info<<endl;
    exit(212);
  }
  // Arguments of zgetrs
  char trans='N'; // no transpose
  int nrhs=1; // nr of columns of B
  int lda=dimL;int ldb=dimL;
  double* result=B.copyComponents(); // copy the values in B
  zgetrs_(&trans,&dimL,&nrhs,U.components,&lda,ipiv,result,&ldb,&info);

  // Now return result in X
  if(info!=0){
    if(info<0)
      cout<<"Error: illegal input argument for zgelsd at position "
	  <<-info<<endl;
    else
      cout<<"Error: failed to converge in zgelsd, info="<<info<<endl;
    exit(212);
  }
  // Rearrange things
  X=mwArray(Indices(dimL,1),result);
  delete [] ipiv; delete[] result;
}

