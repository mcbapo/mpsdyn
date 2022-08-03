#ifndef NEWIMPLEM
#include "libdiminf.h"
#include "diminf.h"
#else
#include "mwArray.h"
#endif
#include "MPS.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>

#ifdef MACOSX
#define LOCALDIR "/Users/banulsm/Desktop/SVN/projects/mpsdyn/src/bin/"
//#define LOCALDIR "/Users/banuls/SVN/peps/projects/mpsdyn/src/bin/"
#else
#define LOCALDIR "/home/banulsm/ALaPorra/bin/"
#endif

using namespace std;

#include "time.h"
clock_t start,finish;

/** 
    Test program Nr. 4: Check the correct behaviour of MPS class
    Compile with make test4
*/

int main(){

 char filename[120];

#ifndef NEWIMPLEM
  int error;
  initializelib(&error);
  try{
    initGlobals(mwArray(CHECKDIMS),mwArray(SEEDVAL));
#endif

  MPS A;
  cout<<"Constructed default MPS A="<<A<<endl;

  int len=30;int bond=10;int d=4;
  MPS B(len,bond,d);
  cout<<"Constructed an MPS B="<<B<<endl;

  B.setProductState(p_xplus);
  cout<<"Set product state X+ B="<<B<<endl;
  
  B.gaugeCond('R',0);
  cout<<"Applied gauge condition to the right "<<B<<endl;

  B.gaugeCond('L',0);
  cout<<"Applied gauge condition to the left "<<B<<endl;

  MPS C(len,bond,d);MPS D;MPS F;
#ifndef NEWIMPLEM
  //C.importMPStext("/home/banulsm/ALaPorra/bin/miRandomB.mat");
  strcpy(filename,LOCALDIR);
  C.importMPStex(strcat(filename,"miRandomB.mat"));
  cout<<"Read from text file C="<<C<<endl;
#else
  C.setRandomState();
  cout<<"Filled with random values C="<<C<<endl;
  //C.exportMPStext("/home/banulsm/ALaPorra/bin/miRandomB.mat");
  strcpy(filename,LOCALDIR);
  C.exportMPStext(strcat(filename,"miRandomB.mat"));
  cout<<"C exported to text file!"<<endl;
#endif
  D=C;
  cout<<"Copied D="<<D<<endl;

  // Check gauge condition
  C.gaugeCond('L',0);
  cout<<"Applied gauge condition to the left "<<C<<endl;
  C.gaugeCond('R',0);
  cout<<"Applied gauge condition to the right "<<C<<endl;

  MPS E(D,2);
  cout<<"Created from cut E="<<E<<endl;

#ifdef NEWIMPLEM
  strcpy(filename,LOCALDIR);
  F.importMPStext(strcat(filename,"miRandomB.mat"));
  cout<<"F imported from text file! "<<endl;
  // Checking with origianl (D copied)
  for(int k=0;k<F.getLength();k++){
    cout<<"Diff Site "<<k<<" ";
    cout<<F.getA(k).getA()-D.getA(k).getA()<<endl;
  }
  cout<<"MPS F is "<<F<<endl;
  if(pow(d,len)<pow(2,10)){
    mwArray vec;
    expandVec(F,vec);
    cout<<"Now expanding to vector "<<vec<<endl;
  }
  C.exportMPS("miRandomB.mat");
  cout<<"C exported to binary file"<<endl;
  F.importMPS("miRandomB.mat");
  //F.importMPS("coherentOut/MPS_g1_gamma1.5_L30_D10_B.mat");
  cout<<"F imported from binary file! "<<endl;
  F.gaugeCond('R',1);
  cout<<"Applied gauge condition to F"<<endl;
  F.gaugeCond('L',1);

  // Now iterate and measure gauge condition
  MPS mps(len,bond,d);
  
  int nrtests=100;
  start=clock();
  for(int k=0;k<nrtests;k++){
    mps.setRandomState();
  }
  clock_t initOffset=clock()-start;
  start=clock();
  for(int k=0;k<nrtests;k++){
    mps.setRandomState();
    mps.gaugeCond('L',0);
  }
  finish=clock();
  clock_t timeQR=finish-start-initOffset;
  start=clock();
  for(int k=0;k<nrtests;k++){
    mps.setRandomState();
    mps.gaugeCond('L',0,bond);
  }
  finish=clock();
  clock_t timeSVD=finish-start-initOffset;

  cout<<"Average time for gauge with SVD: "<<timeSVD/nrtests
      <<" and with QR: "<<timeQR/nrtests<<endl;

  // Test for increaseBondDimension
  if(0){
    vector<int> dims;
    for(int k=0;k<125;k++){
      dims.push_back(2);
      dims.push_back(3);
    }
    dims.push_back(2);
    MPS testInc(251,10,dims);testInc.setRandomState();
    testInc.gaugeCond('R',1);
    testInc.gaugeCond('L',1);
    
    cout<<"Before increasing Bond dimension, testInc="<<testInc<<endl;
    testInc.increaseBondDimension(20);
    cout<<"After increasing Bond dimension, testInc="<<testInc<<endl;
  }

  // Test the partial trace feature.
  //  MPS rho(8,6,4);rho.setRandomState();
  MPS rho(4,6,4);rho.setRandomState();
  rho.exportForMatlab("initRho.m");
  std::vector<int> sites;
  //  sites.push_back(0);sites.push_back(1);sites.push_back(6);sites.push_back(7);
  sites.push_back(0);sites.push_back(1);sites.push_back(2);sites.push_back(3);
  rho.traceOutSites(sites,0);
  rho.exportForMatlab("tracedRho.m");
  cout<<"Norm of tracedRho="<<rho.getNormFact()<<endl;


#endif

#ifndef NEWIMPLEM
  }
  catch(mwException& e){
    cout<<"Exception caught: "<<e.what()<<"in testMPS"<<endl;
    exit(212);
  }
  closelib();
  return error;
#endif

}
