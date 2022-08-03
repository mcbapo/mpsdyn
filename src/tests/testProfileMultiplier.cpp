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
    Test program Nr. 10: Check the performance of 
    products on mwArray.
    Compile with:    make test11
*/

#include "time.h"
int main(){

#ifndef NEWIMPLEM
  int error;
  initializelib(&error);
  try{
    initGlobals(mwArray(1),mwArray(SEEDVAL));
#endif

  cout<<"Testing square matrix multiplication "<<endl;
  for(int D=10;D<50;D=D+10){
    int D1(D),D2(D),D3(D),D4(D);
#ifdef NEWIMPLEM
    mwArray A(Indices(D1*D2,D3*D4));
    A.fillRandom();
#else
    mwArray A;
    mwArray args(1,2,mxCELL_CLASS);
    args(1)=mwArray(D1*D2);args(2)=mwArray(D3*D4);
    m_rand(1,A,mwArray(1),args);
#endif

    clock_t start=clock();
    int kmax=20;
    for(int k=0;k<kmax;k++){
#ifdef NEWIMPLEM
      mwArray B=A*A;
#else
      mwArray B;
      m_times(1,B,A,A);
#endif
    }
    clock_t finish=clock();
    cout<<"Time for "<<kmax<<" products when D="<<D<<" (total size "<<
      D1*D2<<") = "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
  }

  cout<<"Testing non-square multiplication "<<endl;
  for(int D=10;D<60;D=D+10){
    int D1(D),D2(D),D3(D),D4(D);
#ifdef NEWIMPLEM
    mwArray A(Indices(D1,D2*D3*D4));
    A.fillRandom();
    mwArray B(Indices(D1*D2*D3,D4));
    B.fillRandom();
#else
    mwArray A,B;
    mwArray args(1,2,mxCELL_CLASS);
    args(1)=mwArray(D1);args(2)=mwArray(D2*D3*D4);
    m_rand(1,A,mwArray(1),args);
    args(1)=mwArray(D2*D3*D4);args(2)=mwArray(D1);
    m_rand(1,B,mwArray(1),args);
#endif

    clock_t start=clock();
    int kmax=20;
    for(int k=0;k<kmax;k++){
#ifdef NEWIMPLEM
      mwArray C=A*B;
#else
      mwArray C;
      m_times(1,C,A,B);
#endif
    }
    clock_t finish=clock();
    cout<<"Time for "<<kmax<<" products when A is ="<<D1<<"x"<<
      D2*D3*D4<<") = "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
  }



  cout<<"Testing MM product"<<endl;
  for(int D=2;D<24;D=D+2){
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
  for(int D=2;D<24;D=D+2){
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
