#include "wrapper.h"
#include "mwArray.h"
#include "Indices.h"
#include <cstdio>
#include <vector>

// What I consider to be zero when checking hermiticity (absolute value of matrix elements)
#define PREC_HER 1E-14

using namespace std;
using namespace shrt;

extern "C" { 

  /** Lapack routine. Computes all eigenvalues and, optionally,
      eigenvectors of a symmetric tridiagonal matrix using the divide
      and conquer method. For the eigenvalues of Hermitian matrix,
      zhetrd has to be used before. */
  void zstedc_(char* compz,int* N,double* D,double* E,double* Z,
	       int* ldz,double* WORK,int* lwork,double* RWORK,
	       int* lrwork,int* iwork,int* liwork,int* info);

  /** Lapack routine. Reduces a complex Hermitian matrix A to real
      symmetric tridiagonal form T by a unitary similarity
      transformation: \f[Q^{\dagger} * A * Q = T\f].
  */
  void zhetrd_(char* UPLO,int* N,double* A,int* lda,
	       double* D,double* E,
	       double* TAU,double* WORK,
	       int* lwork,int* info);

  /** Lapack routine. Computes the unitary generating the tridiagonal
      transformation of zhetrd, so that it can be passed to zstedc to
      compute full eigensystem. */
  void zungtr_(char* UPLO,int* N,double* A,int* lda,double* TAU,double* WORK,
	       int* lwork,int* info);

  /** Lapack routine: computes eigenvalues of a non-Hermitian
      matrix. Can also compute right or left eigenvectors. */
  void zgeev_(char* JOBVL,char* JOBVR,int* N,double* A,int* LDA,double* W,
	      double* VL,int* LDVL,double* VR,int* LDVR,double* WORK,
	      int* LWORK,double* RWORK,int* info);
}
  

void wrapper::eig(const mwArray& A,std::vector<complex_t>& Dval,mwArray& U,
		  bool eigv){
  if(A.isScalar()){
    vector<int> indx(A.getRank(),0);
    complex_t val=A.getElement(indx);
    Dval.clear();Dval.push_back(val);U=mwArray(ONE_c);
    return;
  }
  if(!A.isMatrix()){
    cout<<"Error: cannot use eig on array "<<A<<endl;
    exit(212);
  }
  //cout<<"Calling eig(A)"<<endl;
  Dval.clear();U.clear();
  if(A.isHermitian(PREC_HER)){
    mwArray AH=.5*(A+Hconjugate(A));
    //cout<<"Before calling LAPACK, AH:"<<AH<<endl;
    // Prepare to call Lapack zhetrd to reduce to tridiagonal symmetric
    char UPLO='U'; // really does not matter
    int N=AH.getDimension(0);
    double* D=new double[N]; // place for diagonal elements
    double* E=new double[N-1]; // place for off-diagonal elements
    double* tau=new double[2*(N-1)];
    int lwork=-1; // will query the routine for required value
    int info=0;
    double workT=0.;
    //cout<<"Calling zhetrd with N="<<N<<")"<<endl;
    //Query for optimal size of work space
    zhetrd_(&UPLO,&N,AH.components,&N,D,E,tau,&workT,&lwork,&info);
    //cout<<"After 1st zhetrd_"<<endl;
    if(info!=0){
      cout<<"Error: exit value info="<<info<<" when querying zhetrd on "
	  <<AH.getDimensions()<<endl;
      exit(212);
    }
    // Prepare space and call the real operation
    lwork=(int)workT;
    double* work=new double[2*lwork];
    zhetrd_(&UPLO,&N,AH.components,&N,D,E,tau,work,&lwork,&info);
    //cout<<"After 2nd zhetrd_"<<endl;
    if(info!=0){
      cout<<"Error: exit value info="<<info<<" when trying zhetrd on "
	  <<AH.getDimensions()<<endl;
      exit(212);
    }
    // some cleaning up
    delete[] work;
    //cout<<"After deleting work"<<endl;
    if(eigv){ // I want also the eigenvectors=> need the unitary used for tridiagonal reduction
      // query the size of workspace
      lwork=-1;
      zungtr_(&UPLO,&N,AH.components,&N,tau,&workT,&lwork,&info);
      // prepare workspace and proceed
      lwork=(int)workT;
      double* work=new double[2*lwork];
      zungtr_(&UPLO,&N,AH.components,&N,tau,work,&lwork,&info);
      if(info!=0){
	cout<<"Error: exit value info="<<info<<" when trying zungtr on "<<A.getDimensions()<<endl;
	exit(212);
      }
      delete[] work;
    }
    delete[] tau;
    // Now prepare for diagonalization
    char COMPZ=eigv?'V':'N';
    int ldz=AH.getDimension(0); //simply N
    lwork=-1;int lrwork=-1;int liwork=-1; // query for optimal sizes
    double rworkT=0.;int iworkT=0.;
    // Query for optimal work space
    zstedc_(&COMPZ,&N,D,E,AH.components,&ldz,&workT,&lwork,&rworkT,&lrwork,
	    &iworkT,&liwork,&info);
    // Prepare proper work space and diagonalize
    lwork=workT;lrwork=rworkT;liwork=iworkT;
    work=new double[2*lwork];
    double* rwork=new double[2*lrwork];
    int* iwork=new int[2*liwork];
    for(int ikk=0;ikk<liwork;ikk++){iwork[2*ikk]=0;iwork[2*ikk+1]=0;}
    for(int ikk=0;ikk<lrwork;ikk++){rwork[2*ikk]=0.;rwork[2*ikk+1]=0.;}
    zstedc_(&COMPZ,&N,D,E,AH.components,&ldz,work,&lwork,rwork,&lrwork,
	    iwork,&liwork,&info);  
    // Process output and clean up
    if(info<0){
      cout<<"Error: exit value info="<<info<<" when trying zstedc on "<<A.getDimensions()<<endl;
      exit(212);
    }
    //cout<<"AH="<<AH<<endl;
    if(eigv) 
      U=AH;
    for(int k=0;k<N;k++) Dval.push_back((complex_t){D[k],0.});
    delete[] work;delete[] rwork;delete[] iwork;
    delete[] D;delete[] E;
  }
  else{
    cout<<"Error: right now I cannot diagonalize a non-Hermitian matrix"<<endl;
    // cout<<"Error in eig, called for A:"<<endl;
    // putForMatlab(cout,A,"A");
    exit(212);
  }

}

void wrapper::eigNH(const mwArray& A,std::vector<complex_t>& Dval,mwArray& U,
		    bool eigv,bool left){
  if(A.isHermitian(PREC_HER)){
    // call the Hermitian one
    return wrapper::eig(A,Dval,U,eigv);
  }
  // If non-Hermitian, call the zeev routine from Lapack
  // I need to copy the matrix, because the argument was const
  mwArray Anh(A);
  char JOBVL=left&eigv?'V':'N'; //left EV computed
  char JOBVR=(!left)&eigv?'V':'N'; //right EV computed
  int N=Anh.getDimension(0);
  int LDA=Anh.getDimension(0);
  double* W=new double[2*N]; // place for eigenvalues
  mwArray VL,VR;
  int LDVL(N),LDVR(N);
  if(left&eigv) VL=mwArray(Indices(LDVL,N));
  if((!left)&eigv) VR=mwArray(Indices(LDVR,N));
  int LWORK=3*N+2;
  double* work=new double[2*LWORK];
  double* rwork=new double[2*2*N];
  int info=0;
  zgeev_(&JOBVL,&JOBVR,&N,Anh.components,&LDA,W,VL.components,&LDVL,VR.components,&LDVR,
	 work,&LWORK,rwork,&info);
  // Process output and clean up
  if(info!=0){
    cout<<"Error: exit value info="<<info<<" when trying zgeev on "<<A<<endl;
    exit(212);
  }
  // Copy eigenvalues in place
  Dval.clear();
  for(int k=0;k<N;k++){
     Dval.push_back((complex_t){W[2*k],W[2*k+1]});
  }
  // If eigenvalues were computed, return them now
  if(eigv){
    if(left) U=VL;
    else U=VR;
  }
  // clean up
  delete [] work;delete [] rwork; delete [] W;
}
