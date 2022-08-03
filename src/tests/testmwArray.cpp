#include "Indices.h"
#include "mwArray.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include "Indices.h"

#ifdef USING_PRIMME
extern "C" {
#include "primme.h"
}
#endif

using namespace std;
using namespace shrt;

/** 
    Test program Nr. 1: Check the correct behaviour of mwArray class
    and wrapper routines.
    Compile with:     make test1
*/

int main(){
  // Create a default mwArray
  mwArray A1;
  cout<<setprecision(4)<<"Created default array A1="<<A1<<endl;
  // Create with given size
  int rank=3;
  Indices dims(2,4,3);
  cout<<"My Indices object is "<<dims<<endl;
  mwArray A2(dims);
  cout<<"And create a matrix of size "<<A2.getDimensions()<<endl;
  //exit(1);
  cout<<"Created sized array A2="<<A2<<endl;
  // Create from a file
  char filename[]="data.dat";
  ifstream infile(filename,ios::binary);
  if(!infile.is_open()){
    cout<<"Error: unable to open file "<<filename<<endl;
  }
  else{
    mwArray A3(infile);
    cout<<"Created array from file A3="<<A3<<endl;
    infile.close();
  }
  // Fill with random values between 0 and 1
  A2.fillRandom();
  cout<<"Filled array randomly A2="<<A2<<endl;
  // Save in a file
  ofstream outfile(filename,ios::binary);
  if(!outfile.is_open()){
    cout<<"Error: unable to open file to write"<<filename<<endl;
  }
  else{
    A2.save(outfile);
    outfile.close();
    cout<<"Saved A2="<<A2<<endl;
  }
  // SubArray extraction
  Indices subIdx(-1,2,-1);
  mwArray subA2=A2.subArray(subIdx);
  cout<<"Extracted subArray "<<subIdx<<": "<<subA2<<endl;

  // Copy constructor
  mwArray A4(A2);
  cout<<"Copied mwArray A2->A4="<<A4<<endl;

  // Permute indices
  Indices permut(2,3,1);
  A4.permute(permut);
  cout<<"Permuted mwArray A4="<<A4<<endl;

  // Construct from literal
  double data[]={1., 1., 2., 2., 3., 3., 4., 4.};
  complex_t dataC[]={{1., 1.}, {2., 2.}, {3., 3.}, {4., 4.}};
  //complex_t dataC[]={{1., 1.}, {2., 2.}, {3., 3.}, ONE_c};
  //double data[]={1., 0., 2., 0., 3., 0., 4., 0.};
  Indices size(2,2);
  mwArray A5(size,data);
  cout<<"Constructed from literal A5="<<A5<<endl;
  mwArray A5C(size,dataC);
  cout<<"Constructed from literal complex A5C="<<A5C<<endl;
  //vector<int> indices(2,1);indices[0]=2;
  // A5.permute(indices);
  // cout<<"Permuted mwArray A5="<<A5<<endl;
  
  // Try products
  mwArray A6=A5*A5;
  cout<<fixed<<"Computed product A6=A5*A5="<<A6<<endl;

  // Try sums
  mwArray A7=A5+A6;
  cout<<fixed<<"Computed sum A7=A5+A6="<<A7<<endl;

  // Try sums
  //mwArray A8=A5-A5;
  //cout<<fixed<<"Computed sum A8=A5-A5="<<A8<<endl;

  // Try product with scalar
  double alpha[2]={0.,1.}; //i
  mwArray A9=alpha*A5;
  cout<<"Computed A9=alpha*A5="<<A9<<endl;

  mwArray A9C=(complex_t){0.,1.}*A5;
  cout<<"Computed A9C=alphaC*A5="<<A9C<<endl;

  // Set isolated component
  Indices pos(0,0);
  pos[0]=1;
  A9.setElement(-3.,-5.,pos);
  cout<<"Changed component A9="<<A9<<endl;
 
  A9.setElement((complex_t){-6.,-10.},pos);
  cout<<"Changed component with complex A9="<<A9<<endl;
 
  complex_t eX;
  eX=A9.getElement(pos);
  //  cout<<"eX="<<eX<<endl;
  // cout<<"Saliendo voluntariamente"<<endl;  exit(1);
  cout<<"Read complex element @"<<pos<<"="<<A9.getElement(pos)<<endl;

  A9.setRealData(4,1.,2.,3.,4.);
  cout<<"setRealData A9="<<A9<<endl;

  // Try resize to bigger
  size[0]=4;size[1]=3;
  A9.resize(size);
  cout<<"Resized A9="<<A9<<endl;

  A9.setImagData(12,4.,3.,2.,1.,4.,3.,2.,1.,4.,3.,2.,1.);
  cout<<"setImagData A9="<<A9<<endl;
  
  complex_t tmpArray[12]={{4.,1.},
		 {3.,1.},
		   {2.,1.},
		     {1.,-1.},
		       {4.,-1.},
			 {3.,-2.},
			   {2.,.5},
			     {1.,0.},{4.,3.},{3.,.666},{0.,2.},{1.}};

  A9.setData(12,tmpArray);
  cout<<"setData A9="<<A9<<endl;

  // Try resize to smaller
  size[0]=3;size[1]=2;
  A9.resize(size);
  cout<<"Resized A9="<<A9<<endl;

  // Construct a diagonal matrix
  double data1[]={1., 2., 3.};
  mwArray A10=realDiag(3,data1);
  cout<<"Constructed diagonal A10="<<A10<<endl;

  double data2[]={1.,1., 2.,2., 3.,3.};
  mwArray A11=diag(3,data2);
  cout<<"Constructed diagonal A11="<<A11<<endl;

  vector<complex_t> dataCv(3,ONE_c);
  mwArray A11b=diag(dataCv);
  cout<<"Constructed diagonal A11b from vector="<<A11b<<endl;

  // Try the square root
  mwArray A12=sqrt(A9);
  cout<<"A12=sqrt(A9)="<<A12<<endl;
  
  // Try the Hermitian conjugate
  mwArray A13=Hconjugate(A12);
  cout<<"A13=sqrt(A9)+="<<A13<<endl;

  /**
  int d1=5;int d2=7;
  for(int k2=0;k2<d2;k2++)
    for(int k1=0;k1<d1;k1++){
      cout<<"El ("<<k1<<","<<k2<<") was at "<<k1+d1*k2
	  <<" -> "<<k2+d2*k1<<endl;
    }

    exit(1); 
*/
  // And the trace
  double reT,imT;
  A10.trace(reT,imT);
  cout<<"Trace of matrix A10 is "<<reT<<"+i"<<imT<<endl;
  cout<<"Trace of A11 is "<<A11.trace()<<endl;

  cout<<"Are they diagonal? A9->"<<A9.isDiagonal()<<endl;
  cout<<"   A10->"<<A10.isDiagonal()<<endl;
  cout<<"   A12->"<<A12.isDiagonal()<<endl;
  cout<<"   A6->"<<A6.isDiagonal()<<endl;

  // Test invertDiagonal
  int aux;
  mwArray A14=invertDiag(A10,aux);
  cout<<"Constructed inverse of A10: A14="<<A14<<endl;
  mwArray A15=invertDiag(A11,aux);
  cout<<"Constructed inverse of A11: A15="<<A15<<endl;

  // Test SVD

  mwArray A16(A9.getDimensions());A16.fillRandom();
  mwArray S,Vdagger,U;
  wrapper::svd(A16,U,S,Vdagger);
  cout<<"SVD of A16="<<A16<<endl;
  cout<<" gives U="<<U<<endl;
  cout<<"\t S="<<S<<endl;
  cout<<"\t Vdagger="<<Vdagger<<endl;
  cout<<"CHECKING: U S Vdagger-A16="<<U*S*Vdagger-A16<<endl;

  mwArray ident=identityMatrix(4);
  cout<<"Created identity matrix "<<ident<<endl;

  mwArray vectorR(Indices(1,4));vectorR.fillRandom();
  mwArray vectorC(Indices(4,1));vectorC.fillRandom();
  mwArray vector0(Indices(1,1));vector0.fillRandom();

  cout<<"Vector test: ident? "<<ident.isVector()
      <<", vectorR? "<<vectorR.isVector()
      <<", vectorC? "<<vectorC.isVector()<<", vector0? "
      <<vector0.isVector()<<endl;

  mwArray M0(Indices(4,4));M0.fillRandom();
  cout<<"VectorR="<<vectorR<<endl;
  cout<<"VectorC="<<vectorC<<endl;
  cout<<"Matrix="<<M0<<endl;
  cout<<"Product Matrix-vector:"<<M0*vectorC<<endl;
  cout<<"Product vector(T)-Matrix:"<<vectorR*M0<<endl;
  cout<<"Product vectorR vectorC:"<<vectorR*vectorC<<endl;

  // Test outerproduct too
  mwArray VVT=outerproduct(vectorR,vectorC);
  cout<<"Outer product VR.VCt="<<VVT<<endl;

  cout<<"Hconjugate(VVT)="<<Hconjugate(VVT)<<endl;
  mwArray herm=.5*(VVT+Hconjugate(VVT));
  // Test eig routine!
  mwArray Uvec;vector<complex_t> values;
  wrapper::eig(herm,values,Uvec,1);
  cout<<"After eig for VVT,, values are ";
  for(int k=0;k<4;k++)cout<<values[k]<<", ";
  cout<<"and vectors "<<Uvec<<endl;

  //Now try eigs
  mwArray herm2(Indices(4,4));herm2.fillRandom();
  herm2=herm2+Hconjugate(herm2);
  cout<<"Random hermitian "<<herm2<<endl;
  wrapper::eig(herm2,values,Uvec,1);
  cout<<"After eig for random hermitian, values are ";
  for(int k=0;k<4;k++)cout<<values[k]<<", ";
  cout<<"and vectors "<<Uvec<<endl;
  Uvec.clear();
  wrapper::eigs(herm2,2,"SM",values,Uvec,1);
  cout<<"After eigs, values are ";
  for(int k=0;k<2;k++)cout<<values[k]<<", ";
  cout<<"and vectors "<<Uvec<<endl;

#ifdef USING_PRIMME
  int err=wrapper::eigs_primme(herm2,2,primme_smallest,values,Uvec,1);
  if(err!=0){
    cout<<"ERROR in the diagonalization!"<<endl;
    exit(err);
  }
  cout<<"After eigs_primme, values are ";
  for(int k=0;k<2;k++)cout<<values[k]<<", ";
  cout<<"and vectors "<<Uvec<<endl;

#endif

  // Try solution of linear equations
  mwArray X(Indices(4,1));X.fillRandom(); // actual sol
  cout<<"A="<<herm2<<endl;
  cout<<"X="<<X<<endl;
  mwArray B=herm2*X; // rhs
  mwArray sol;
  wrapper::lsd(herm2,B,sol);
  cout<<"After lsd, solution: "<<sol<<endl;
  wrapper::lss(herm2,B,sol);
  cout<<"After lss, solution: "<<sol<<endl;
  wrapper::lslu(herm2,B,sol);
  cout<<"After lslu, solution: "<<sol<<endl;
  //exit(1);

  // Try the LU decomposition
  mwArray LUtest(Indices(4,3));LUtest.fillRandom();
  mwArray L,Uu,P;
  wrapper::lu(LUtest,L,Uu,P);
  cout<<"After LU of A="<<LUtest<<endl;
  cout<<"L="<<L<<endl;
  cout<<"U="<<Uu<<endl;
  cout<<"P="<<P<<endl;
  cout<<"P*L*U-A="<<P*L*Uu-LUtest<<endl;
  //exit(1);

  // Try the QR decomposition
  mwArray Q1,R1;
  wrapper::qr(LUtest,Q1,R1);
  cout<<"After QR of A="<<LUtest<<endl;
  cout<<"Q="<<Q1<<endl;
  cout<<"R="<<R1<<endl;
  cout<<"Q*R-A="<<Q1*R1-LUtest<<endl;
  //  exit(1);

  // And the LQ decomposition
  mwArray L1;
  wrapper::lq(LUtest,L1,Q1);
  cout<<"After LQ of A="<<LUtest<<endl;
  cout<<"Q="<<Q1<<endl;
  cout<<"L="<<L1<<endl;
  cout<<"L*Q-A="<<L1*Q1-LUtest<<endl;


  // Try now with defective rank
  cout<<"A="<<VVT<<endl;
  cout<<"X="<<X<<endl;
  B=VVT*X; // rhs
  wrapper::lsd(VVT,B,sol);
  cout<<"After lsd, solution: "<<sol<<endl;

  // Test permutation identity
  cout<<"Permute A4 "<<A4<<" as (2,1,3)"<<endl;
  cout<<permute(A4,Indices(2,1,3))<<endl;

  // Test shallow mwArray
  A4.permute(Indices(2,3,1));
  cout<<"Now A4 "<<A4<<" permuted as (2,3,1)"<<endl;
  mwArray Ashallow;
  Ashallow.setPointer(A4.getDimensions(),(double*)A4.getComponents());
  
  mwArray A80=permute(A4,Indices(3,1,2));
  cout<<"The copy is Ashallow="<<Ashallow<<endl;
  Ashallow.permute(Indices(3,1,2));
  cout<<"Now permute the shallow copy as (3,1,2) and get A4:"
      <<reshape(A4,A80.getDimensions())<<endl;
  cout<<"Diff wrt direct permuting: "<<A80-reshape(A4,A80.getDimensions())<<endl;

  cout<<"Test new trasposition!"<<endl;
  mwArray A81(Indices(2,5));A81.fillRandom();
  cout<<"Original A81="<<A81<<endl;
  mwArray A82=Hconjugate(A81);
  cout<<"Transposed A82="<<A82<<endl;
  mwArray A83=reshape(A81,Indices(2,1,5));
  A83.permute(Indices(3,2,1),true);A83.reshape(Indices(5,2));
  cout<<"Old permute gives"<<A83<<endl;
  cout<<"Difference "<<A83-A82<<endl;

  cout<<"Test expm!"<<endl;
  mwArray expH;
  wrapper::expm(herm2,expH);
  //cout<<setprecision(10);
  cout<<"Computed exponential of "<<herm2<<endl;
  cout<<"Result="<<expH<<endl;

  cout<<"Test inverse"<<endl;
  mwArray A90(Indices(5,5));A90.fillRandom();
  int nrVals;
  mwArray A91=inverse(A90,nrVals);
  //cout<<setprecision(10)<<endl;
  cout<<"Computed inverse of "<<A90<<endl;
  cout<<"with "<<nrVals<<" non vanishing values"<<endl;
  cout<<"Result: "<<A91<<endl;
  cout<<"Check: A*inv(A)-Id "<<A90*A91-identityMatrix(A90.getDimension(0))
      <<endl;
  cout<<"Check: inv(A)*A-Id "<<A90*A91-identityMatrix(A90.getDimension(0))
      <<endl;

  cout<<"Test single leg contraction"<<endl;
  int dimC=5;
  int d1(1),d2(4),d3(5),d4(5),d5(1),d6(4);
  Indices dims1(d1,d2,d3), order1(1,4,2,3); int dim1L=d1*d2*d3;
  Indices dims2(d4,d5,d6), order2(2,3,1,4); int dim2R=d4*d5*d6;
  mwArray A92(Indices(dims1,dimC));A92.fillRandom();
  mwArray A93(Indices(dimC,dims2));A93.fillRandom();
  mwArray A94=reshape(A92,Indices(dim1L,dimC))*reshape(A93,Indices(dimC,dim2R));
  A94.reshape(Indices(dims1,dims2));
  // now reshuffle and contract
  A92.permute(Indices(order1));
  A93.permute(Indices(order2));
  mwArray A95=contractLeg(A92,2,A93,3);
  cout<<"Difference between both forms of contracting:"<<endl;
  cout<<A94-A95<<endl;

  cout<<"Test multiple leg contraction"<<endl;
  Indices dimsC(2,3,2); int cTot=dimsC[0]*dimsC[1]*dimsC[2];
  mwArray A96(Indices(dim1L,cTot));A96.fillRandom();
  mwArray A97(Indices(cTot,dim2R));A97.fillRandom();
  mwArray A98=A96*A97;A98.reshape(Indices(dims1,dims2));
  A96.reshape(Indices(dims1,dimsC)); // ordered d1 d2 d3 dc1 dc2 dc3
  A97.reshape(Indices(dimsC,dims2)); // ordered dc1 dc2 dc3 d4 d5 d6
  order1=Indices(1,4,5,2,6,3);A96.permute(order1);
  order2=Indices(4,1,5,2,6,3);A97.permute(order2);
  mwArray A99=contractLegs(A96,Indices(2,3,5),A97,Indices(2,4,6));
  cout<<"Difference between both forms of contracting:"<<endl;
  cout<<A98-A99<<endl;
}
