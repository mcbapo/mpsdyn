#include "wrapper.h"
#include "mwArray.h"
#include <cstdio>
#include <vector>

using namespace std;

extern "C" {
  // SVD  
  void zgesvd_(char* JOBU,char* JOBVT,int* m,int* n,
	       double* A,int* lda,double* S,double* U,
	       int* ldu,double* VT,int* ldvt,
	       double* work,int* lwork, 
	       double* rwork,int* info);
}

void wrapper::svd(const mwArray& A,mwArray& U,mwArray& S,mwArray& Vdagger,bool full){
  int nr=0;
  double tol=0;
  svd(A,tol,nr,U,S,Vdagger,full);
}

void wrapper::svd(const mwArray& A,const double& tol,int& Xi,
		  mwArray& U,mwArray& S,mwArray& Vdagger,bool full){
  // Particular case: a scalar
  if(A.isScalar()){
    complex_t elem=A.getElement(vector<int>(A.getRank(),0)); // single element
    S=mwArray((complex_t){abs(elem),0.});
    Vdagger=mwArray(ONE_c);
    U=mwArray(elem*(1./abs(elem)));
    return;
  }
  if(A.getRank()!=2){
    cout<<"Error: Cannot do SVD on an array if not rank 2, and here dims="
	<<A.getDimensions()<<endl;
    exit(212);
  }
  // cout<<"SVD of A="<<A<<" with tol="<<tol<<", Xi="<<Xi
  //  <<", U="<<U<<", S="<<S<<", Vdagger="<<Vdagger<<endl;
  // Arguments for the fortran routine
  int info;
  char jobu='o';char jobt='s';
  int dimL=A.getDimension(0);int dimR=A.getDimension(1);
  double* D=new double[min(dimL,dimR)];
  vector<int> auxV(2);auxV[0]=min(dimL,dimR);auxV[1]=dimR;
  if(full) auxV[0]=dimR;
  Vdagger.clear();Vdagger.resize(auxV);
  U=A; // copied: replaced when coming back in the economic version
  // Space for temporary computations
  int lwork=2*(2*min(dimL,dimR)+max(dimL,dimR));
  double* work=new double[2*lwork];
  double* rwork=new double[5*min(dimL,dimR)];
  double auxDU=0.;
  if(full){
    jobu='a';jobt='a';
    mwArray Uf(A);U.resize(shrt::Indices(dimL,dimL));
    zgesvd_(&jobu,&jobt,&dimL,&dimR,Uf.components,&dimL,D,U.components,&dimL,
	    Vdagger.components,&auxV[0],
	    work,&lwork,rwork,&info);
    S=realDiag(min(dimL,dimR),D);
    S.resize(shrt::Indices(dimL,dimR)); // to match full U and Vd
    if(tol==0)
      Xi=min(dimL,dimR); // all of them
    else {
      // first decide if tol has to cut off some values
      Xi=0;bool highval=1;
      while(Xi<min(dimL,dimR)&&highval){
	if(D[Xi]>tol) Xi++;
	else highval=0;
      }
    }
    return;    
  }
  // cout<<"Calling ZGESVD with arguments JOBU="<<jobu<<", JOBVT="
  //  <<jobt<<", m="<<dimL<<", n="<<dimR<<" lda="<<dimL
  //  <<", ldu="<<dimL<<", ldtv="<<auxV[0]<<endl;
  zgesvd_(&jobu,&jobt,&dimL,&dimR,U.components,&dimL,D,&auxDU,&dimL,
	  Vdagger.components,&auxV[0],
	  work,&lwork,rwork,&info);
  // fix dimensions and contents
  if(Xi==0){   // If passed as argument, I use that as cutoff
    if(tol==0)
      Xi=min(dimL,dimR); // all of them
    else {
      // first decide if tol has to cut off some values
      Xi=0;bool highval=1;
      while(Xi<min(dimL,dimR)&&highval){
	if(D[Xi]>tol) Xi++;
	else highval=0;
      }
    }
    if(Xi==0){
      cout<<"Error: after svd, no singular value over tol="<<tol<<endl;
      exit(212);
    }
  }
  else{ // if Xi given is larger, reduce it to the proper limit (all values)
    if(Xi>min(dimL,dimR)) Xi=min(dimL,dimR);
    int Xitmp=0; // but maybe I still need to cut
    bool highval=1;
    while(Xitmp<Xi&&highval){
      if(D[Xitmp]>tol) Xitmp++;
      else highval=0;
    }
    Xi=Xitmp;
  }
  vector<int> auxU(2);auxU[1]=Xi;auxU[0]=dimL;
  U.resize(auxU);
  S=realDiag(Xi,D);
  if(Xi!=min(dimL,dimR)){
    auxV[0]=Xi; Vdagger.resize(auxV);
  }
  // TODO: Sth if the svd fails??

  delete[] work;
  delete[] rwork;
  delete[] D;
}
