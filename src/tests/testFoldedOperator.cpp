#ifndef NEWIMPLEM
#include "libdiminf.h"
#include "diminf.h"
#else
#include "mwArray.h"
#endif
#include "FoldedOperator.h"
#include "DoubleOperator.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>

using namespace std;

/** 
    Test program Nr. 8: Check the correct behaviour of FoldedOperator class
    Compile by
    make [WITHMKL=1] test8
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
    int dl=2;int Dll=2;int Drl=2;int d2l=2;
    int dr=2;int Dlr=2;int Drr=2;int d2r=2;
    //int dr=1;int Dlr=1;int Drr=1;int d2r=1;
    //mwArray data(Indices(d,Dl,d,Dr));data.fillRandom();
    //Operator oper(data);
    //cout<<"Constructed operator from random data: "<<oper<<endl;

    mwArray AR(Indices(dr,Dlr,dr,Drr));AR.fillRandom();
    mwArray AL(Indices(dl,Dll,dl,Drl));AL.fillRandom();
    mwArray AL_aux=permute(AL,Indices(1,4,3,2));
    FoldedOperator opF(AL,AR);
    DoubleOperator opD(AL_aux,AR);
    Operator fullOp(opF.getFullData());
    Operator fullOpD(opD.getFullData());

    cout<<"Basic piece opL "<<AL<<endl;
    cout<<"Basic piece opR "<<AR<<endl;

    //cout<<"Constructed block operator "<<fullOp<<endl;
    mwArray both=reshape(AL_aux,Indices(-1,1))*reshape(AR,Indices(1,-1));
    both.reshape(Indices(dl,Drl,d2l,Dll,dr,Dlr,d2r,Drr));
    both.permute(Indices(1,5,2,6,3,7,4,8));
    both.reshape(Indices(dl*dr,Drl*Dlr,d2l*d2r,Dll*Drr));
    Operator C(both);
    //cout<<"Constructed full operator by hand "<<C<<endl;

    // test CONTRACTIONS!
    int betabl=3;int betabr=3; //left right dims of bra
    int betakl=3;int betakr=3; //left right dims of ket
    Site ket(1,d2l*d2r,betakl,betakr);ket.setState(randmps);
    Site bra(1,dl*dr,betabl,betabr);bra.setState(randmps);
    cout<<"************ Contractions **********"<<endl;
    //cout<<"ket="<<ket<<endl;
    //cout<<"bra="<<bra<<endl;

    mwArray result,resultC,resultF,resultD;
    mwArray tmpL2(Indices(betabl,Drl*Dlr,betakl));tmpL2.fillRandom();
    mwArray tmpR2(Indices(betabr,Dll*Drr,betakr));tmpR2.fillRandom();
    //cout<<"tmpL2="<<tmpL2<<endl;
    //cout<<"tmpR2="<<tmpR2<<endl;
    C.contractL(resultC,tmpL2,ket,bra);
    //cout<<"\nOper.contractL->"<<resultC;
    opF.contractL(resultF,tmpL2,ket,bra);
    opD.contractL(resultD,tmpL2,ket,bra);
    //cout<<"\nFoldedOper.contractL->"<<resultF;
    fullOp.contractL(result,tmpL2,ket,bra);
    //cout<<"\nFullOper.contractL->"<<result;
    cout<<"**** contractL differences ********"<<endl;
    cout<<"Oper-Folded="<<resultC-resultF<<endl;
    cout<<"Full-Folded="<<result-resultF<<endl;
    cout<<"Double-Folded="<<resultD-resultF<<endl;
    cout<<"Oper-Full="<<resultC-result<<endl;

    C.contractR(resultC,tmpR2,ket,bra);
    //    cout<<"\nOper.contractR->"<<resultC;
    opF.contractR(resultF,tmpR2,ket,bra);
    opD.contractR(resultD,tmpR2,ket,bra);
    //cout<<"\nFoldedOper.contractR->"<<resultF;
    fullOp.contractR(result,tmpR2,ket,bra);
    //cout<<"\nFullOper.contractL->"<<result;
    cout<<"**** contractR differences ********"<<endl;
    cout<<"Oper-Folded="<<resultC-resultF<<endl;
    cout<<"Full-Folded="<<result-resultF<<endl;
    cout<<"Double-Folded="<<resultD-resultF<<endl;
    cout<<"Oper-Full="<<resultC-result<<endl;

    C.contractMket(resultC,tmpL2,ket,tmpR2);
    //cout<<"\nOper.contractMket->"<<resultC<<endl;
    opF.contractMket(resultF,tmpL2,ket,tmpR2);
    opD.contractMket(resultD,tmpL2,ket,tmpR2);
    //cout<<"\nFoldedOper.contractMket->"<<resultF<<endl;
    fullOp.contractMket(result,tmpL2,ket,tmpR2);
    //cout<<"\nFullOper.contractMket->"<<result<<endl;
    cout<<"**** contractMket differences ********"<<endl;
    cout<<"Oper-Folded="<<resultC-resultF<<endl;
    cout<<"Full-Folded="<<result-resultF<<endl;
    cout<<"Double-Folded="<<resultD-resultF<<endl;
    cout<<"Oper-Full="<<resultC-result<<endl;

    C.contractMbra(resultC,tmpL2,bra,tmpR2);
    //cout<<"\nOper.contractMbra->"<<resultC<<endl;
    opF.contractMbra(resultF,tmpL2,bra,tmpR2);
    opD.contractMbra(resultD,tmpL2,bra,tmpR2);
    //cout<<"\nFoldedOper.contractMbra->"<<resultF<<endl;
    fullOp.contractMbra(result,tmpL2,bra,tmpR2);
    //cout<<"\nFullOper.contractMbra->"<<result<<endl;
    cout<<"**** contractMbra differences ********"<<endl;
    cout<<"Oper-Folded="<<resultC-resultF<<endl;
    cout<<"Full-Folded="<<result-resultF<<endl;
    cout<<"Double-Folded="<<resultD-resultF<<endl;
    cout<<"Oper-Full="<<resultC-result<<endl;

    C.contractN(resultC,tmpL2,tmpR2);
    //cout<<"\nOper.contractN->"<<resultC<<endl;
    opF.contractN(resultF,tmpL2,tmpR2);
    opD.contractN(resultD,tmpL2,tmpR2);
    //cout<<"\nFoldedOper.contractN->"<<resultF<<endl;
    fullOp.contractN(result,tmpL2,tmpR2);
    //cout<<"\nFullOper.contractN->"<<result<<endl;
    cout<<"**** contractN differences ********"<<endl;
    cout<<"Oper-Folded="<<resultC-resultF<<endl;
    cout<<"Full-Folded="<<result-resultF<<endl;
    cout<<"Double-Folded="<<resultD-resultF<<endl;
    //cout<<"Oper-Full="<<resultC-result<<endl;

    cout<<"Dimensions of C="<<C.getDimensions()<<endl;
    cout<<"Dimensions of Folded C="<<opF.getDimensions()<<endl;
    cout<<"And Full C is "<<fullOp.getDimensions()<<endl;

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
