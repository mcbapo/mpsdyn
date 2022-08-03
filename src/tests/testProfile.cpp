#ifndef NEWIMPLEM
#include "libdiminf.h"
#include "diminf.h"
#else
#include "mwArray.h"
#include "Indices.h"
#endif
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision

using namespace std;
#ifdef NEWIMPLEM
using namespace shrt;
#endif

/** 
    Test program Nr. 10: Check the performance of permute
    and product on mwArray.
    Compile with:    make test10
*/

#include "time.h"
int main(){

#ifndef NEWIMPLEM
  int error;
  initializelib(&error);
  try{
    initGlobals(mwArray(1),mwArray(SEEDVAL));
#endif

  cout<<"Testing permutations"<<endl;
  for(int D=10;D<50;D=D+10){
    int D1(D),D2(D),D3(D),D4(D);
#ifdef NEWIMPLEM
    mwArray A(Indices(D1,D2,D3,D4));
    Indices permutation(2,4,3,1);A.fillRandom();
    A.fillRandom();
#else
    mwArray A;
    mwArray args(1,4,mxCELL_CLASS);
    args(1)=mwArray(D1);args(2)=mwArray(D2);
    args(3)=mwArray(D3);args(4)=mwArray(D4);
    m_rand(1,A,mwArray(1),args);
    mwArray permutation(1,4,mxDOUBLE_CLASS);
    permutation(1)=mwArray(2);permutation(2)=mwArray(4);
    permutation(3)=mwArray(3);permutation(4)=mwArray(1);
#endif

    clock_t start=clock();
    int kmax=20;
    for(int k=0;k<kmax;k++){
#ifdef NEWIMPLEM
      mwArray B=permute(A,permutation);
#else
      mwArray B;
      m_permute(1,B,A,permutation);
#endif
    }
    clock_t finish=clock();
    cout<<"Time for "<<kmax<<" permutations when D="<<D<<" (total size "<<
      D1*D2*D3*D4<<") = "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
  }

  cout<<"Testing asymmetric transpositions"<<endl;
  for(int D=10;D<60;D=D+10){
    int D1(D),D2(D),D3(D),D4(D);
#ifdef NEWIMPLEM
    mwArray A(Indices(D1,D2*D3*D4));
    Indices permutation(2,1);A.fillRandom();
#else
    mwArray A;
    mwArray args(1,2,mxCELL_CLASS);
    args(1)=mwArray(D1);args(2)=mwArray(D2*D3*D4);
    m_rand(1,A,mwArray(1),args);
    mwArray permutation(1,2,mxDOUBLE_CLASS);
    permutation(1)=mwArray(2);permutation(2)=mwArray(1);
#endif

    clock_t start=clock();
    int kmax=20;
    for(int k=0;k<kmax;k++){
#ifdef NEWIMPLEM
      mwArray B=permute(A,permutation);
#else
      mwArray B;
      m_permute(1,B,A,permutation);
#endif
    }
    clock_t finish=clock();
    cout<<"Time for "<<kmax<<" transpositions when D="<<D<<" (total size "<<
      D1*D2*D3*D4<<") = "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
  }

  cout<<"Testing reshape"<<endl;
  for(int D=20;D<20;D=D+10){
    int D1(D),D2(D),D3(D),D4(D);
#ifdef NEWIMPLEM
    mwArray A(Indices(D1,D2,D3,D4));
    Indices newdims(D1*D2,D3*D4);A.fillRandom();
    A.fillRandom();
#else
    mwArray A;
    mwArray args(1,4,mxCELL_CLASS);
    args(1)=mwArray(D1);args(2)=mwArray(D2);
    args(3)=mwArray(D3);args(4)=mwArray(D4);
    m_rand(1,A,mwArray(1),args);
    mwArray newdims(1,2,mxCELL_CLASS);
    newdims(1)=mwArray(D1*D2);newdims(2)=mwArray(D3*D4);
#endif

    clock_t start=clock();
    for(int k=0;k<100;k++){
#ifdef NEWIMPLEM
      mwArray B=reshape(A,newdims);
#else
      mwArray B;
      m_reshape(1,B,A,newdims);
#endif
    }
    clock_t finish=clock();
    cout<<"Time for 100 reshapes when D="<<D<<" (total size "<<
      D1*D2*D3*D4<<") = "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
  }

  cout<<"Testing copying"<<endl;
  for(int D=20;D<20;D=D+10){
    int D1(D),D2(D),D3(D),D4(D);
#ifdef NEWIMPLEM
    mwArray A(Indices(D1,D2,D3,D4));
    Indices newdims(D1*D2,D3*D4);A.fillRandom();
    A.fillRandom();
#else
    mwArray A;
    mwArray args(1,4,mxCELL_CLASS);
    args(1)=mwArray(D1);args(2)=mwArray(D2);
    args(3)=mwArray(D3);args(4)=mwArray(D4);
    m_rand(1,A,mwArray(1),args);
#endif

    clock_t start=clock();
    for(int k=0;k<100;k++){
#ifdef NEWIMPLEM
      mwArray B=A;
#else
      mwArray B;
      m_copy(1,B,A);
#endif
    }
    clock_t finish=clock();
    cout<<"Time for 100 copies when D="<<D<<" (total size "<<
      D1*D2*D3*D4<<") = "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
  }


  cout<<"Testing MM product"<<endl;
  for(int D=24;D<24;D=D+2){
    int D1(D),D2(D);
#ifdef NEWIMPLEM
    mwArray M(Indices(D1,D2));M.fillRandom();
#else
    mwArray M;
    mwArray args(1,2,mxCELL_CLASS);
    args(1)=mwArray(D1);args(2)=mwArray(D2);
    m_rand(1,M,mwArray(1),args);
#endif
    clock_t start=clock();
    for(int k=0;k<100;k++){
#ifdef NEWIMPLEM
      mwArray B=M*M;
#else
      mwArray B;
      m_times(1,B,M,M);
#endif
    }
    clock_t finish=clock();
    cout<<"Time for 100 products M*M when D="<<D<<" (total size "<<
      D1*D2<<") = "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
  }

  cout<<"Testing MV product"<<endl;
  for(int D=24;D<24;D=D+2){
    int D1(D),D2(D);
#ifdef NEWIMPLEM
    mwArray M(Indices(D1,D2));M.fillRandom();
    mwArray V(Indices(D2,1));V.fillRandom();
#else
    mwArray M,V;
    mwArray args(1,2,mxCELL_CLASS);
    args(1)=mwArray(D1);args(2)=mwArray(D2);
    m_rand(1,M,mwArray(1),args);
    args(2)=1;
    m_rand(1,V,mwArray(1),args);
#endif
    clock_t start=clock();
    for(int k=0;k<100;k++){
#ifdef NEWIMPLEM
      mwArray B=M*V;
#else
      mwArray B;
      m_times(1,B,M,V);
#endif
    }
    clock_t finish=clock();
    cout<<"Time for 100 products M*V when D="<<D<<" (total size "<<
      D1*D2<<") = "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
  }

  cout<<"Testing reshape-permute-reshape"<<endl;
  for(int D=60;D<60;D=D+10){
    int D1(D),D2(D),D3(D),D4(D);
#ifdef NEWIMPLEM
    mwArray A(Indices(D1*D2,D3*D4));
    Indices permutation(2,4,3,1);
    Indices newdims(D2,D4,D3,D1);
    Indices olddims(D1,D2,D3,D4);
    A.fillRandom();
#else
    mwArray A;
    mwArray args(1,4,mxCELL_CLASS);
    args(1)=mwArray(D1);args(2)=mwArray(D2);
    args(3)=mwArray(D3);args(4)=mwArray(D4);
    mwArray newargs(1,4,mxCELL_CLASS);
    newargs(1)=mwArray(D1);newargs(2)=mwArray(D2);
    newargs(3)=mwArray(D3);newargs(4)=mwArray(D4);
    mwArray argsDim(1,2,mxCELL_CLASS);
    argsDim(1)=mwArray(D1*D2);argsDim(2)=mwArray(D3*D4);
    m_rand(1,A,mwArray(1),argsDim);
    mwArray permutation(1,4,mxDOUBLE_CLASS);
    permutation(1)=mwArray(2);permutation(2)=mwArray(4);
    permutation(3)=mwArray(3);permutation(4)=mwArray(1);
#endif

    clock_t start=clock();
    int kmax=100;
    for(int k=0;k<kmax;k++){
#ifdef NEWIMPLEM
      //mwArray B=reshape(permute(reshape(A,olddims),permutation),newdims);
      mwArray B=reshape(A,olddims);
      B.permute(permutation);
      B.reshape(newdims);
#else
      mwArray B;
      testfun(1,B,A,args,newargs);
#endif
    }
    clock_t finish=clock();
    cout<<"Time for "<<kmax<<" reshape-permute-reshape when D="<<D<<" (total size "<<
      D1*D2*D3*D4<<") = "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
  }

  
#ifndef NEWIMPLEM
  }
  catch(mwException& e){
    cout<<"Exception caught: "<<e.what()<<" in testProfile"<<endl;
    exit(212);
  }
  closelib();
  return error;
#endif

}
