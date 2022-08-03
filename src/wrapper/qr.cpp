#include "wrapper.h"
#include "mwArray.h"
#include <cstdio>
#include <vector>

using namespace std;
using namespace shrt;

extern "C" {
  // QR decomposition for general matrices  
  void zgeqrf_(int* m,int* n,
	       double* A,int* lda,
	       double* tau,double* work,
	       int* lwork,int* info);
  // Get the unitary Q matrix from the output of the above
  void zungqr_(int* m,int* n,
	       int* k,double* A,int* lda,
	       double* tau,double* work,
	       int* lwork,int* info);

  // LQ decomposition for general matrices  
  void zgelqf_(int* m,int* n,
	       double* A,int* lda,
	       double* tau,double* work,
	       int* lwork,int* info);
  // Get the unitary (isometric) Q matrix from the output of the above
  void zunglq_(int* m,int* n,
	       int* k,double* A,int* lda,
	       double* tau,double* work,
	       int* lwork,int* info);
}

void wrapper::qr(const mwArray& A,mwArray& Q,mwArray& R){
  // Particular case: a scalar, save the call
  if(A.isScalar()){
    Q=mwArray(ONE_c);Q.reshape(Indices(1,1));
    R=A;R.reshape(Indices(1,1));
    //cout<<"QR decomposition for a scalar "<<A<<
    //", Q="<<Q<<", R="<<R<<endl;
    return;
  }
  if(A.getRank()!=2){
    cout<<"Error: Cannot do QR on an array if not rank 2, and here dims="
	<<A.getDimensions()<<endl;
    exit(212);
  }
  // Arguments for the fortran routine
  int info;
  int M=A.getDimension(0);int N=A.getDimension(1);
  R=A; // copied: replaced when coming back
  int lda=M;
  // Space for temporary computations and Lapack return values
  double* tau=new double[2*min(M,N)];
  int lwork=-1; // will ask about size
  double workspacesize[2]; // place for the optimal size
  // cout<<"Calling zgeqrf_ with lwork="<<lwork<<", workspacesize="<<*workspacesize<<endl;
  zgeqrf_(&M,&N,R.components,&lda,tau,workspacesize,&lwork,&info); // query about size
  if(info<0){
    cout<<"Error in query for size to zgeqrf_, illegal value of parameter "<<-info<<endl;
    exit(212);
  }
  lwork=(int)(workspacesize[0]);
   // cout<<"After first call to zgeqrf, desired working space "<<lwork
   //     <<" and the other params: M="<<M<<", N="<<N<<", lda="<<lda
   //     <<endl;
  double* work=new double[2*lwork];
  zgeqrf_(&M,&N,R.components,&lda,tau,work,&lwork,&info);
  if(info<0){
    cout<<"Error in call to zgeqrf_, illegal value of parameter "<<-info<<endl;
    exit(212);
  }
  // cout<<"After Lapack zgeqrf, R="<<R
  //   <<" and the other params: M="<<M<<", N="<<N<<", lda="<<lda
  //   <<endl;
  // Now call zungqr to reconstruct the unitary
  Q=R; // will be destroyed
  int K=min(M,N);
  // cout<<"Will call now zungqr with K="<<K<<endl;
  // the unitary is partially reconstruted to M x min(M,N)
  lwork=-1; // again query for size
  zungqr_(&M,&K,&K,Q.components,&lda,tau,workspacesize,&lwork,&info);
  if(info<0){
    cout<<"Error in query for size to zungqr_, illegal value of parameter "<<-info<<endl;
    exit(212);
  }
  lwork=(int)(workspacesize[0]);
  // cout<<"After first call to zungqr, desired working space "<<lwork
  //     <<" and the other params: M="<<M<<", N="<<N<<", lda="<<lda<<", K="<<K
  //    <<endl;
  K=min(M,N); /** Contrary to the documentation in Lapack, K is changed after the first call! */
  // cout<<"Will call now zungqr with K="<<K<<endl;
  delete [] work;
  work=new double[2*lwork];
  zungqr_(&M,&K,&K,Q.components,&lda,tau,work,&lwork,&info);
  if(info<0){
    cout<<"Error in call to zungqr_, illegal value of parameter "<<-info<<endl;
    exit(212);
  }
  // cout<<"After Lapack zungqr, K="<<K<<" R="<<R<<", Q="<<Q<<endl;
  // now resize matrices and cleanup
  Q.resize(Indices(M,K));
  R.resize(Indices(K,N));
  // cout<<"After resizing to K="<<K<<", R="<<R<<", Q="<<Q<<endl;
  // set zeros below the diagonal of R
  for(int k=1;k<K;k++)
    for(int j=0;j<k;j++)
      R.setElement(ZERO_c,Indices(k,j));

  delete [] work;
  delete [] tau;
}

void wrapper::lq(const mwArray& A,mwArray& L,mwArray& Q){
  // Particular case: a scalar, save the call
  if(A.isScalar()){
    Q=mwArray(ONE_c);
    L=A;
    return;
  }
  if(A.getRank()!=2){
    cout<<"Error: Cannot do LQ on an array if not rank 2, and here dims="
	<<A.getDimensions()<<endl;
    exit(212);
  }
  // Arguments for the fortran routine
  int info;
  int M=A.getDimension(0);int N=A.getDimension(1);
  L=A; // copied: replaced when coming back
  int lda=M;
  // Space for temporary computations and Lapack return values
  double* tau=new double[2*min(M,N)];
  int lwork=-1; // will ask about size
  double workspacesize[2]; // place for the optimal size
  zgelqf_(&M,&N,L.components,&lda,tau,workspacesize,&lwork,&info); // query about size
  if(info<0){
    cout<<"Error in query for size to zgelqf_, illegal value of parameter "<<-info<<endl;
    exit(212);
  }
  lwork=(int)workspacesize[0];
  // cout<<"After first call to zgelqf, desired working space "<<lwork
  //     <<" and the other params: M="<<M<<", N="<<N<<", lda="<<lda
  //     <<endl;
  double* work=new double[2*lwork];
  zgelqf_(&M,&N,L.components,&lda,tau,work,&lwork,&info);
  if(info<0){
    cout<<"Error in call to zgelqf_, illegal value of parameter "<<-info<<endl;
    exit(212);
  }
  // cout<<"After Lapack zgelqf, L="<<L
  //     <<" and the other params: M="<<M<<", N="<<N<<", lda="<<lda
  //     <<endl;
  // Now call zungqr to reconstruct the unitary
  Q=L; // will be destroyed
  int K=min(M,N);
  //cout<<"Will call now zungqr with K="<<K<<endl;
  // the unitary is partially reconstruted to M x min(M,N)
  lwork=-1; // again query for size
  zunglq_(&K,&N,&K,Q.components,&lda,tau,workspacesize,&lwork,&info);
  if(info<0){
    cout<<"Error in query for size to zunglq_, illegal value of parameter "<<-info<<endl;
    exit(212);
  }
  lwork=(int)workspacesize[0];
  // cout<<"After first call to zunglq, desired working space "<<lwork
  //     <<" and the other params: M="<<M<<", N="<<N<<", lda="<<lda<<", K="<<K
  //    <<endl;
  K=min(M,N); /** Contrary to the documentation in Lapack, K is changed after the first call! */
  // cout<<"Will call now zunglq with K="<<K<<endl;
  delete [] work;
  work=new double[2*lwork];
  zunglq_(&K,&N,&K,Q.components,&lda,tau,work,&lwork,&info);
  if(info<0){
    cout<<"Error in call to zunglq_, illegal value of parameter "<<-info<<endl;
    exit(212);
  }
  // cout<<"After Lapack zunglq, K="<<K<<" L="<<L<<", Q="<<Q<<endl;
  // now resize matrices and cleanup
  Q.resize(Indices(K,N));
  L.resize(Indices(M,K));
  // cout<<"After resizing to K="<<K<<", L="<<L<<", Q="<<Q<<endl;
  // set zeros below the diagonal of R
  for(int k=1;k<K;k++)
    for(int j=0;j<k;j++)
      L.setElement(ZERO_c,Indices(j,k));

  delete [] work;
  delete [] tau;
}
