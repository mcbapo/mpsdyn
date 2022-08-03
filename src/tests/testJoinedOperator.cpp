#ifndef NEWIMPLEM
#include "libdiminf.h"
#include "diminf.h"
#else
#include "mwArray.h"
#endif
#include "JoinedOperator.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>

using namespace std;

/** 
    Test program Nr. 5: Check the correct behaviour of JoinedOperator class
*/

using namespace std;
using namespace shrt;

int main(){

#ifndef NEWIMPLEM
  int error;
  initializelib(&error);
  try{
    initGlobals(mwArray(CHECKDIMS),mwArray(SEEDVAL));
#endif
    int d=2;int Dl=2;int Dr=2;int d2=2;
    mwArray data(Indices(d,Dl,d,Dr));data.fillRandom();
    Operator oper(data);
    cout<<"Constructed from random data "<<oper<<endl;

    mwArray A(Indices(d,Dl,d,Dr));A.fillRandom();
    Operator op1(A), op2(conjugate(permute(A,Indices(3,2,1,4))));
    const Operator* ops[2]={&op2,&op1};

    cout<<"Basic operator op1"<<op1<<endl;
    cout<<"Basic operator op2(op1+)"<<op2<<endl;


    Operator C(2,ops);
    cout<<"Constructed block operator "<<C<<endl;
    JoinedOperator jC(2,ops);
    cout<<"Constructed joined operator "<<jC.getFullData()<<endl;

    // test CONTRACTIONS!
    int betabl=3;int betabr=3; //left right dims of bra
    int betakl=3;int betakr=3; //left right dims of ket
    Site ket(1,d2,betakl,betakr);ket.setState(randmps);
    Site bra(1,d,betabl,betabr);bra.setState(randmps);
    cout<<"************ Contractions **********"<<endl;
    cout<<"ket="<<ket<<endl;
    cout<<"bra="<<bra<<endl;

    mwArray result,resultJ;
    mwArray tmpL2(Indices(betabl,Dl*Dl,betakl));tmpL2.fillRandom();
    mwArray tmpR2(Indices(betabr,Dr*Dr,betakr));tmpR2.fillRandom();
    cout<<"tmpL2="<<tmpL2<<endl;
    cout<<"tmpR2="<<tmpR2<<endl;
    C.contractL(result,tmpL2,ket,bra);
    cout<<"\nOper.contractL->"<<result;
    jC.contractL(resultJ,tmpL2,ket,bra);
    cout<<"\nJoinedOper.contractL->"<<resultJ;
    C.contractR(result,tmpR2,ket,bra);
    cout<<"\nOper.contractR2->"<<result;
    jC.contractR(resultJ,tmpR2,ket,bra);
    cout<<"\nJoinedOper.contractR2->"<<resultJ;
    //cout<<"The double operator is "<<oper.doubleop(0)<<endl;

    C.contractMket(result,tmpL2,ket,tmpR2);
    cout<<"\nOper.contractMket->"<<result<<endl;
    jC.contractMket(resultJ,tmpL2,ket,tmpR2);
    cout<<"\nJoinedOper.contractMket->"<<resultJ<<endl;

    C.contractMbra(result,tmpL2,bra,tmpR2);
    cout<<"\nOper.contractMbra->"<<result<<endl;
    jC.contractMbra(resultJ,tmpL2,bra,tmpR2);
    cout<<"\nJoinedOper.contractMbra->"<<resultJ<<endl;

    C.contractN(result,tmpL2,tmpR2);
    cout<<"\nOper.contractN->"<<result<<endl;
    jC.contractN(resultJ,tmpL2,tmpR2);
    cout<<"\nJoinedOper.contractN->"<<resultJ<<endl;

    cout<<"Dimensions of C="<<C.getDimensions()<<endl;
    cout<<"Dimensions of Joined C="<<jC.getDimensions()<<endl;
    cout<<"And Joined C is "<<jC<<endl;

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
