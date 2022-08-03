#include "wrapper.h"
#include "mwArray.h"

#ifdef MKL_VER
//#include "mkl.h"
#endif

extern "C" {
  // Product matrix vector
  void zgemv_(char* trans,const int* m,const int* n,const double* alpha,
	      const double* A,const int* lda,const double* X,
	      const int* incx,const double* beta,double* Y,const int*incy);

  // Product of matrices
  void zgemm_(char* transA,char* transB,int* dimL,int* dimR,int* dim,
	     double* alpha,const double* A,int* lda,const double* B,int* ldb,
	     double* beta,double* C,int* ldc);
}

using namespace std;
using namespace shrt;

/*
bool operator<(const Indices& dimA,const Indices& dimB){
  if(dimA.size()!=dimB.size()) return 0;
  for(int k=0;k<dimA.size();k++){
    if(dimA[k]>=dimB[k]) return 0;
  }
  return 1;
}
*/

ostream& operator<<(ostream& os,const vector<int>& dims){
  os<<"(";
  for(int li=0;li<dims.size();li++){
    os<<dims[li];
    if(li<dims.size()-1)
      os<<",";
  }
  return os<<")";
}

ostream& operator<<(ostream& os,const vector<double>& dvals){
  os<<"(";
  for(int li=0;li<dvals.size();li++){
    os<<dvals[li];
    if(li<dvals.size()-1)
      os<<",";
  }
  return os<<")";
}

ostream& operator<<(ostream& os,const vector<complex_t>& cvals){
  os<<"(";
  for(int li=0;li<cvals.size();li++){
    os<<cvals[li];
    if(li<cvals.size()-1)
      os<<",";
  }
  return os<<")";
}

void wrapper::product(const mwArray& A,const mwArray& B,mwArray& C){
  if(A.getNrComponents()==0||B.getNrComponents()==0){
    cout<<"Error: calling wrapper::product with an argument of zero components!: A:"
	<<A.getDimensions()<<" ("<<A.getNrComponents()<<"), B:"
	<<B.getDimensions()<<" ("<<B.getNrComponents()<<")"<<endl;
    //cout<<"A="<<A<<endl;
    //cout<<"B="<<B<<endl;
    exit(1);
  }
  if(A.isScalar()||B.isScalar()){
    cout<<"Error: shouldn't be calling wrapper::product for multiplication "
	<<"with a scalar"<<endl;
    exit(212);
  }
  if(A.getDimension(A.getRank()-1)==B.getDimension(0)){
    // Special case: matrix vector, vector matrix or vector vector products
    if(A.isVector()||B.isVector()){
      if(A.isMatrix()){ // matrix-vector
	//cout<<"Product matrix vector"<<endl;
	// Prepare place for the result
	int dimL=A.getDimension(0);int dimR=A.getDimension(1);
	if(C.nrcomponents!=dimL)
	  if(C.mine){
	    C.clear();
	    C.nrcomponents=dimL;
	    C.components=new double[2*C.nrcomponents];
	  }
	  else{
	    cout<<"Cannot put the result of product in shallow mwArray "
		<<C.dims<<endl;
	    exit(212);
	  }
	C.dims=Indices(dimL,1);
	C.rank=2;
	char trans='N';
	double alpha[]={1.,0.}; // coeff to call ZGEMV
	double beta[]={0.,0.};
	int lda=dimL;
	int incx=1;int incy=1;
	zgemv_(&trans,&dimL,&dimR,alpha,A.getComponents(),
	       &lda,B.getComponents(),&incx,beta,C.components,&incy);
      }
      else if(B.isMatrix()){ // vector matrix 
	//cout<<"Product vector matrix"<<endl;
	int dimL=B.getDimension(0);int dimR=B.getDimension(1);
	if(C.nrcomponents!=dimR)
	  if(C.mine){
	    C.clear();
	    C.nrcomponents=dimR;
	    C.components=new double[2*C.nrcomponents];
	  }
	  else{
	    cout<<"Cannot put the result of product in shallow mwArray "
		<<C.dims<<endl;
	    exit(212);
	  }
	C.dims=Indices(1,dimR);
	C.rank=2;
	char trans='T';
	double alpha[]={1.,0.}; // coeff to call ZGEMV
	double beta[]={0.,0.};
	int lda=dimL;
	int incx=1;int incy=1;
	zgemv_(&trans,&dimL,&dimR,alpha,B.getComponents(),
	       &lda,A.getComponents(),&incx,beta,C.components,&incy);
      }
      else{ // vector-vector -> either scalar or tensor product
	if(B.getDimension(0)==1)// outer product
	  C=outerproduct(A,B);
	else
	  C=scalarproduct(A,B);
      }
    }
    else{ // true matrix-matrix multiplication
      if(!A.isMatrix()||!B.isMatrix()){
	cout<<"ERROR: product should only be called on 2-dimensional tensors, but here A:"<<A.getDimensions()<<", B:"<<B.getDimensions()<<endl;
	exit(212);
      }
      int dimL=A.getDimension(0);int dim=A.getDimension(1);
      int dimR=B.getDimension(1);
      if(dim!=B.getDimension(0)){
	cout<<"Error: operands to A*B do not have the same dimension "
	    <<"(A_lastD="<<dim<<", B_firstD="<<B.getDimension(0)<<")"<<endl;
	exit(212);
      }
      // Prepare place for the result
      if(C.nrcomponents!=dimL*dimR)
	if(C.mine){
	  C.clear();
	  C.nrcomponents=dimL*dimR;
	  C.components=new double[2*C.nrcomponents];
	}
	else{
	  cout<<"Cannot put the result of product in shallow mwArray "
	      <<C.dims<<endl;
	  exit(212);
	}
      C.dims=Indices(dimL,dimR);
      C.rank=2;
      // TEST CASE - see precision for small matrices
      if(0&&dimL*dim*dimR<500){
	// do it brute force, instead of calling lapack
	for(int i_=0;i_<dimL;i_++)
	  for(int j_=0;j_<dimR;j_++){
	    complex_t val=ZERO_c;
	    for(int k_=0;k_<dim;k_++)
	      val=val+A.getElement(Indices(i_,k_))*B.getElement(Indices(k_,j_));
	    C.setElement(val,Indices(i_,j_));
	  }
	return;
      }
      // TEST CASE
      char trans='N';
      double alpha[]={1.,0.}; // coeff to call ZGEMM
      double beta[]={0.,0.};
      int lda=dimL;
      int ldb=dim;
      int ldc=dimL;
      // Call the fortran routine
      zgemm_(&trans,&trans,&dimL,&dimR,&dim,alpha,A.getComponents(),
	     &lda,B.getComponents(),&ldb,beta,C.components,&ldc);
    }
  }
  else{
    cout<<"Error: invalid dimensions "<<A.getDimensions()
	<<" and "<<B.getDimensions()<<" for mwArray product"<<endl;
    exit(212);
  }
}


void wrapper::expm(const mwArray& A, mwArray& result,complex_t coeff){
  //cout<<"In expm for matrix "<<A.getDimensions()<<", coeff="<<coeff<<endl; 
  mwArray U;
  vector<complex_t> D;
  // 1) Diagonalize the (full) matrix
  eig(A,D,U,true);
  //cout<<"After diagonalization, D "<<D.size()<<endl; 
  // 2) exponentiate the diagonal form
  int N=D.size();
  mwArray Dmat(Indices(N,N));
  for(int k=0;k<N;k++){
    Dmat.setElement(exp(coeff*D[k]),Indices(k,k));
  }
  // 3) go back to original basis: exp(A) = U exp(D) U'
  result=U;
  // TODO: use a trick: Multiply U*exp(D) column by column
  result.multiplyRight(Dmat);
  U.Hconjugate();
  result.multiplyRight(U);
}

void wrapper::expmTaylor(const mwArray& A, mwArray& result,int order,
			 complex_t coeff){
  //cout<<"In expm for matrix "<<A.getDimensions()<<", coeff="<<coeff<<endl; 
  int k=0;
  result=identityMatrix(A.getDimension(0)); // order 0
  complex_t coeffK=coeff;
  mwArray AK=A;
  bool converged=false;
  while(k<order){ //&&!converged){
    //if(abs(coeffK)<1E-100) converged=true;
    //else{
    //    result=result+coeffK*AK;
    result=result+coeff*AK;
    k++;
    //cout<<"expmTaylor order "<<k<<" -> coeff="<<coeffK<<endl;
    //coeffK=coeffK*(coeff*ONE_c/(k+1));
    AK=(coeff*(1./(k+1))*A)*AK;
    //}
  }
}
