#ifndef NEWIMPLEM
#include "libdiminf.h"
#include "diminf.h"
#else
#include "mwArray.h"
#include "MPO.h"
#include "MPOMultiplier.h"
#endif
#include "Operator.h"
#include "Contractor.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>
#include <sys/time.h>

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

using namespace std;
using namespace shrt;

#define NUMREPE 100

int main(){

  char filename[120];

  int len=10;
  int d=2;int Dl=2;int Dr=2;
  MPOMultiplier mpo1(len);
  cout<<"Constructed MPO length "<<len<<", "<<mpo1<<endl;
  
  mwArray data1(Indices(d,1,d,Dr));data1.fillRandom();
  mwArray data2(Indices(d,Dl,d,Dr));data2.fillRandom();
  mwArray data3(Indices(d,Dl,d,1));data3.fillRandom();
  Operator oper1(data1);
  Operator oper2(data2);
  Operator oper3(data3);
  //  const Operator *tmpOp[6]={&oper1,&oper2,&oper2,&oper2,&oper2,&oper3};
  //MPO mpo2(len,tmpOp);
  mpo1.setRotatedOp(0,data1,Indices(1,2,3,4));
  for(int k=1;k<len-1;k++)
    mpo1.setRotatedOp(k,data2,Indices(1,2,3,4));
  mpo1.setRotatedOp(len-1,data3,Indices(1,2,3,4));
  cout<<"Constructed MPOMultiplier from several opers "<<mpo1<<endl;
  
    MPS test1(len,10,d);
  test1.setRandomState();
  test1.gaugeCond('R',1);
  cout<<"Constructed random MPS"<<endl;

  Contractor& contractor=Contractor::theContractor();
  cout<<"1) Contraction of the MPO with the random MPS="<<endl;
  MPS aux(test1);
  struct timeval start,final;
  long long timer=0;
  gettimeofday(&start,NULL);
  complex_t value;
  for(int l=0;l<NUMREPE;l++)
    contractor.optimize(mpo1,test1,aux);
  gettimeofday(&final,NULL);
  value=contractor.contract(aux,aux);
  timer=(final.tv_sec-start.tv_sec)*1E3+(final.tv_usec-start.tv_usec)*1E-3;
  cout<<"\t"<<value<<endl;
  cout<<"\t TIME:"<<timer<<endl;

  cout<<"2) Contraction of the MPO with the expanded MPS vector"<<endl;
  mwArray vect,result;
  expandVec(test1,vect);
  gettimeofday(&start,NULL);
  for(int l=0;l<NUMREPE;l++)
    mpo1.product(vect,result);
  gettimeofday(&final,NULL);
  timer=(final.tv_sec-start.tv_sec)*1E3+(final.tv_usec-start.tv_usec)*1E-3;
  cout<<"\t"<<conjugate(permute(result,Indices(2,1)))*result<<endl;
  cout<<"\t TIME:"<<timer<<endl;

  cout<<"3) Contraction of the expanded operator with the expanded vector "<<endl;
  mwArray bigOp;
  mpo1.getFullTensor(bigOp);
  gettimeofday(&start,NULL);
  for(int l=0;l<NUMREPE;l++)
    result=bigOp*vect;
  gettimeofday(&final,NULL);
  timer=(final.tv_sec-start.tv_sec)*1E3+(final.tv_usec-start.tv_usec)*1E-3;
  mwArray vecT(result);vecT.transpose(true);
  cout<<"\t"<<vecT*result<<endl;
  cout<<"\t TIME:"<<timer<<endl;

}
