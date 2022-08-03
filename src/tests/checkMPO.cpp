#ifndef NEWIMPLEM
#include "libdiminf.h"
#include "diminf.h"
#else
#include "mwArray.h"
#include "MPO.h"
#endif
#include "Operator.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>


#ifndef NEWIMPLEM
typedef OperatorRow MPO;
#endif

/** 
    An auxiliary program that takes an MPO (of length L, bond
    dimension D and physical dimension d, specified as input
    arguments) from a file, and saves it (old implementation)
    in the matlab format, to be used for checking the correctness
    of an MPO.
    Receives arguments:
    \param <L> (int) length of the chain
    \param <D> (int) bond dimension of the MPO
    \param <d> (int) physical dimension
    \param <fname> (char*) name of the (text) file containing the MPO
    \param <outfname> (char*) name of the (Matlab) file where the MPO 
    is to be written 

    To compile:
    make checkMPO
*/

using namespace std;
#ifdef NEWIMPLEM
using namespace shrt;
#endif

int main(int argc,const char* argv[]){
  int cntr=0;
  int L=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  int d=atoi(argv[++cntr]);
  const char* fname=argv[++cntr];
  const char* outfname=argv[++cntr];

  cout<<"Starting chekMPO program"<<endl;


#ifndef NEWIMPLEM
  int error;
  initializelib(&error);
  try{
    initGlobals(mwArray(CHECKDIMS),mwArray(SEEDVAL));
    MPO mpo(L,D,d);
    cout<<"About to open file "<<fname<<endl;
    mpo.importMPOtext(fname);
    cout<<"MPO mpo loaded from text file "<<fname<<endl;
    for(int k=0;k<mpo.getLength();k++){
      cout<<"In the read mpo, size of op "<<k<<" is "<<mpo.getOpData(k).GetDimensions()<<endl;
    }
    mpo.saveOperatorArray(outfname);
    cout<<"MPO mpo saved to text file "<<outfname<<endl;

  }
  catch(mwException& e){
    cout<<"Exception caught: "<<e.what()<<"in testMPO"<<endl;
    exit(212);
  }
  closelib();
  return error;
#else
  cout<<"Error: checkMPO only to be used with old code!"<<endl;
  exit(212);
#endif
}
