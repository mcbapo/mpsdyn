#ifndef NEWIMPLEM
#include "libdiminf.h"
#include "diminf.h"
#else
#include "Indices.h"
#include "mwArray.h"
#endif
#include "Operator.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>


/** 
    Test program Nr. 3: Check the correct behaviour of Operator class
    Compilation: make test3
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
    int d=2;int Dl=1;int Dr=4;int d2=4;
    Operator op(d,Dl,Dr,d2);
    cout<<"Constructed Operator op="<<op<<endl;

    mwArray data(Indices(d,Dl,d,Dr));data.fillRandom();
    Operator oper(data);
    cout<<"Constructed from random data "<<oper<<endl;

    mwArray A(Indices(2,2,2,2));A.fillRandom();
    mwArray B(Indices(2,2,2,2));B.fillRandom();

    cout<<"From A="<<A<<endl;
    cout<<"And B="<<B<<endl;

    Operator op1(A), op2(B);
    const Operator* ops[2]={&op1,&op2};

    Operator C(2,ops);
    cout<<"Constructed joined operator "<<C<<endl;

    // Assignment (rotated)
    Operator D(A,Indices(1,3,4,2),true);
    cout<<"New operator D "<<D<<endl;
    D.setRotatedData(B,Indices(1,4,3,2),true);
    cout<<"Set data D="<<D<<endl;
    // test CONTRACTIONS!
    Site ket(1,d,3,3);ket.setState(randmps);
    Site bra(1,d,3,3);bra.setState(randmps);
    mwArray tmpL(Indices(3,Dl,3));tmpL.fillRandom();
    mwArray tmpR(Indices(3,Dr,3));tmpR.fillRandom();
    cout<<"************ Contractions **********"<<endl;
    cout<<"ket="<<ket<<endl;
    cout<<"bra="<<bra<<endl;
    cout<<"tmpL="<<tmpL<<endl;
    cout<<"tmpR="<<tmpR<<endl;
    cout<<"oper="<<oper<<endl;

    mwArray result;
    oper.contractL(result,tmpL,ket,bra);
    cout<<"contractL->"<<result;
    oper.contractR(result,tmpR,ket,bra);
    cout<<"contractR->"<<result;
    mwArray tmpL2(Indices(3,Dl*Dl,3));tmpL2.fillRandom();
    mwArray tmpR2(Indices(3,Dr*Dr,3));tmpR2.fillRandom();
    cout<<"tmpL2="<<tmpL2<<endl;
    cout<<"tmpR2="<<tmpR2<<endl;
    oper.contractL2(result,tmpL2,ket,bra);
    cout<<"contractL2->"<<result;
    oper.contractR2(result,tmpR2,ket,bra);
    cout<<"contractR2->"<<result;
    //cout<<"The double operator is "<<oper.doubleop(0)<<endl;

    oper.contractMket(result,tmpL,ket,tmpR);
    cout<<"contractMket->"<<result<<endl;

    oper.contractMbra(result,tmpL,bra,tmpR);
    cout<<"contractMbra->"<<result<<endl;

    oper.contractN(result,tmpL,tmpR);
    cout<<"contractN->"<<result<<endl;

    oper.contractN2(result,tmpL2,tmpR2);
    cout<<"contractN2->"<<result<<endl;


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
