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


#ifdef MACOSX
//#define LOCALDIR "/Users/banulsm/Desktop/SVN/projects/mpsdyn/src/bin/"
#define LOCALDIR "/Users/banuls/SVN/peps/projects/mpsdyn/src/bin/"
#else
#define LOCALDIR "/home/banulsm/ALaPorra/bin/"
#endif

/** 
    Test program Nr. 5: Check the correct behaviour of MPO class
    To compile: make test5
*/

#ifdef NEWIMPLEM
using namespace std;
using namespace shrt;
#endif

int main(){

  char filename[120];

#ifndef NEWIMPLEM
  int error;
  initializelib(&error);
  try{
    initGlobals(mwArray(CHECKDIMS),mwArray(SEEDVAL));
#endif

    int len=6;
    int d=2;int Dl=2;int Dr=2;
#ifdef NEWIMPLEM
    MPO mpo1(len);
    cout<<"Constructed MPO length "<<len<<", "<<mpo1<<endl;

    mwArray data1(Indices(d,1,d,Dr));data1.fillRandom();
    mwArray data2(Indices(d,Dl,d,Dr));data2.fillRandom();
    mwArray data3(Indices(d,Dl,d,1));data3.fillRandom();
    Operator oper1(data1);
    Operator oper2(data2);
    Operator oper3(data3);
    const Operator *tmpOp[6]={&oper1,&oper2,&oper2,&oper2,&oper2,&oper3};
    MPO mpo2(len,tmpOp);
    cout<<"Constructed MPO from array of opers "<<mpo2<<endl;
    
    vector<const mwArray*> tensors;
    tensors.push_back(&data1);tensors.push_back(&data2);    
    tensors.push_back(&data3);
    MPO mpo3(3);
    mpo3.setOperatorArray(tensors);
    cout<<"Constructed MPO from array of mwArray "<<mpo3<<endl;

    cout<<"Dimensions of that MPO "<<mpo3.getDimensions()<<endl;
    cout<<"Original dimensions of that MPO "<<mpo3.getOriginalDimensions()
	<<endl;

    // Test saving and recovering
    strcpy(filename,LOCALDIR);
    mpo2.saveMPO(strcat(filename,"mpo3.dat"));
    cout<<"MPO mpo3 saved to binary file "<<endl;
    MPO mpo4(len);
    strcpy(filename,LOCALDIR);
    mpo4.loadMPO(strcat(filename,"mpo3.dat"));
    cout<<"MPO mpo4 load from binary file "<<mpo4<<endl;

    // Same from text file
    strcpy(filename,LOCALDIR);
    mpo2.exportMPOtext(strcat(filename,"mpo3.txt"));
    cout<<"MPO mpo3 saved to text file "<<endl;
    MPO mpo5(len);
#else
    MPO mpo5(len,Dl,d);
#endif
    strcpy(filename,LOCALDIR);
    mpo5.importMPOtext(strcat(filename,"mpo3.txt"));
    cout<<"MPO mpo5 load from text file "<<mpo5<<endl;

    // Now do some test
#ifdef NEWIMPLEM
    MPO folded(len/2);
    mpo5.fold(folded);
    // For len odd case
    //int dims[]={d,d*d,d*d,d*d,d*d}; 
    // for even len
    int dims[]={1,d*d,d*d,d*d,d*d}; 
    MPS test1(folded.getLength(),2,dims);
    test1.setRandomState();
    test1.gaugeCond('R',1);
    strcpy(filename,LOCALDIR);
    test1.exportMPStext(strcat(filename,"miRandomB.mat"));
#else
    const MPO* foldedptr=mpo5.fold();
    const MPO& folded=*foldedptr;
    MPS test1(folded.getLength(),2,d*d);
    strcpy(filename,LOCALDIR);
    test1.importMPStext(strcat(filename,"miRandomB.mat"));
#endif
    Contractor& contractor=Contractor::theContractor();
    cout<<"Contraction of the folded MPO with a random MPS="
	<<contractor.contract(test1,folded,test1)<<endl;


#ifndef NEWIMPLEM
  }
  catch(mwException& e){
    cout<<"Exception caught: "<<e.what()<<"in testMPO"<<endl;
    exit(212);
  }
  closelib();
  return error;
#endif

}
