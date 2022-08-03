#include <math.h>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "JoinedOperator.h"
#include "Properties.h"
#include "mwArray.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0


using namespace shrt;
using namespace std;

void commutatorMagnetization(MPO& mpo,int L);

// Copied from Liouvillian, should be somewhere else
void constructOperatorProduct(mwArray& result,const mwArray& opA,
			      const mwArray& opB);


const string jobfilename(int M,int D,double hp,const string mpsdir);
const string mpsfilename(int M,int D,const string mpsdir);

#define MAXITER 100
#define MAXPERT 10

int d=2;
double tau=0.8;
bool traceless=true;
bool norand=false;

int main(int argc,const char* argv[]){

  int cntr=0;
  const char* infile=argv[++cntr];

  Properties props(infile);

  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  stringstream commandstr;
  for(int c=0;c<argc;c++) commandstr<<argv[c]<<" ";


  int M=props.getIntProperty("M");
  int D=props.getIntProperty("D");
  string mpsdir=props.getProperty("mpsdir");
  string outfile=props.getProperty("outfile");
  string inputMPS=mpsfilename(M,D,mpsdir);
  const string tmpinfile = inputMPS+".tmp";  

  Contractor& contractor=Contractor::theContractor(); //inizialize contractor

  //Test outfile
  ofstream* out;
  out=new ofstream(outfile.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output"<<endl;
    exit(1);
  }
  out->close();delete out;

  
  MPO commutatorTotalSz(M);
  commutatorMagnetization(commutatorTotalSz,M);  

  //Import or define initial state
  MPS Omps(M,D,d*d);

  if(file_exists(tmpinfile)){
    Omps.importMPS(tmpinfile.data());
    cout<<"Imported tmp file for same D!"<<tmpinfile<<endl;
  }else{  // no tmp file
    cout<<"There is no tmp file for same D: "<<tmpinfile<<endl;
    if(inputMPS.empty()){ // no init file indicated in arguments
      cout<<"There is no input file for initD "<<endl;
      exit(1);
    }
    else{ // initMPS file read
      if(!file_exists(inputMPS)){
	cout<<"There is no input file for initD "<<endl;
	exit(1);
      }else{
      Omps.importMPS(inputMPS.data());
      cout<<"Imported initMPS file for initial state"<<inputMPS<<endl;
      }
    }
  }
  
  //Set gauge condition for the initial state
  Omps.gaugeCond('R',1);Omps.gaugeCond('L',1); 
  
  if(Omps.getBond()<D){
    forceHermitian(Omps);    
    Omps.increaseBondDimensionWithNoise(D,noise);
    substractTrace(Omps);
  }

  double norm=real(contractor.contract(Omps,Omps));
  double particleNumber=real(contractor.contract(Omps,commutatorTotalSz,Omps));

  cout<<"Magnetization commutator = "<<particleNumber/norm<<endl;


  out=new ofstream(outfile.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<"for output"<<endl;
    exit(1);
  }
  *out<<setprecision(15);
  *out<<D<<"\t"<<M<<"\t"<<lambda<<endl;
  out->close();
  delete out;
}




const string mpsfilename(int M,int D,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS"<<"_M"<<M<<"_D"<<D<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

const string jobfilename(int M,int D, double hP, const string jobsdir){
  stringstream s;
  s<<jobsdir<<"/job_H_M"<<M<<"_D"<<D<<"_h"<<hP;  
  s<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}


void constructOperatorProduct(mwArray& result,const mwArray& opA,
			      const mwArray& opB){
#ifdef CHECKDIMS
  if(opA.getDimensions()!=opB.getDimensions()){
    cout<<"Error in LiouvillianXYEdges::constructOperatorProduct for"
	<<" A="<<opA<<" and B="<<opB<<endl;
    exit(1);
  }
#endif
  result=opA; 
  result.reshape(Indices(d*d,1));
  result.multiplyRight(reshape(opB,Indices(1,d*d)));
  result.reshape(Indices(d,d,d,d)); // order here: iji'j'
  result.permute(Indices(1,3,2,4));  // order in the MPO: ii',jj'  
  result.reshape(Indices(d*d,d*d));
}


void commutatorMagnetization(MPO& mpo,int L){

  if(mpo.getLength()!=L)
    mpo.initLength(L);
  // Basic operators
  int nrOps=3;
  int dd=d*d;
  mwArray Z=mwArray(Indices(nrOps,dd,dd));
   // Basic pieces
  mwArray sig0=identityMatrix(d);//identity
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz

  // *** First of all, we need the identity on both
  mwArray term0=identityMatrix(dd);
  // (3) sigZ (x) Id
  mwArray term1;
  constructOperatorProduct(term1,sigZ,sig0);
  // (4)  Id (x) sigZ^T -> just a permutation of the previous one
  mwArray term2;
  constructOperatorProduct(term2,sig0,permute(sigZ,Indices(2,1)));
  // Fill in the operators in the Z array, to be contracted with C in
  // the construction of the MPO. The order is as listed above.
  for(int d1=0;d1<dd;d1++)
    for(int d2=0;d2<dd;d2++){
      Z.setElement(term0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
      Z.setElement(term1.getElement(Indices(d1,d2)),Indices(1,d1,d2));
      Z.setElement(term2.getElement(Indices(d1,d2)),Indices(2,d1,d2));
    }
  Z.reshape(Indices(nrOps,dd*dd));
  //cout<<"Finished Z:"<<Z.getDimensions()<<endl;
  // now the coefficients
  int D=2; // bond dimension of the MPO
  mwArray C1(Indices(1,D,nrOps)); // the first site
  mwArray C(Indices(D,D,nrOps)); // the middle site
  mwArray CN(Indices(D,1,nrOps)); // the last site

  // Identity terms when nothing has happened yet
  C.setElement(ONE_c,Indices(0,0,0));
  C1.setElement(ONE_c,Indices(0,0,0));
  // Identity terms after the end
  C.setElement(ONE_c,Indices(1,1,0));
  CN.setElement(ONE_c,Indices(1,0,0));
  // Single body terms (one for each operator, with proper weights)
  C.setElement(ONE_c,Indices(0,1,1)); // h sigZ x I
  C.setElement(-1.*ONE_c,Indices(0,1,2)); // -h I x sigZ^T
  C1.setElement(ONE_c,Indices(0,1,1)); // h sigZ x I
  C1.setElement(-1.*ONE_c,Indices(0,1,2)); // -h I x sigZ^T
  CN.setElement(ONE_c,Indices(0,0,1)); // h sigZ x I
  CN.setElement(-1.*ONE_c,Indices(0,0,2)); // -h I x sigZ^T

  // Now set the operators in place
  // Now I reshape, contract with Z and set the indices in proper order
  C.reshape(Indices(D*D,nrOps));
  C1.reshape(Indices(1*D,nrOps));
  CN.reshape(Indices(D*1,nrOps));
  C.multiplyRight(Z);C.reshape(Indices(D,D,dd,dd));C.permute(Indices(3,1,4,2));

  C1.multiplyRight(Z);C1.reshape(Indices(1*D*dd,dd));
  // Special situation: first and last sites are contracted with Id as traced out
  mwArray id2=1/sqrt(2.)*identityMatrix(d);id2.reshape(Indices(dd,1));
  C1.multiplyRight(id2);
  C1.reshape(Indices(1,D,dd,1));C1.permute(Indices(3,1,4,2));
  CN.multiplyRight(Z);CN.reshape(Indices(D*1*dd,dd));
  CN.multiplyRight(id2);
  CN.reshape(Indices(D,1,dd,1));CN.permute(Indices(3,1,4,2)); // dd x D x 1 x 1

  // Now I build the MPO with required ops.
  //1) on the edges, I contract C1 and CN with their adjoint
  C1.reshape(Indices(dd,D));
  CN.reshape(Indices(dd,D)); // dd D
  mwArray tmp(C1);tmp.Hconjugate(); 
  C1.multiplyLeft(tmp); // Dl,Dr
  tmp=CN;tmp.Hconjugate(); 
  CN.multiplyLeft(tmp); //Dlu,Dld
  // 2) These, I contract with the middle ones 
  tmp=C;tmp.permute(Indices(2,1,3,4));tmp.reshape(Indices(D,dd*dd*D));
  mwArray edgeL=C1*tmp;edgeL.reshape(Indices(D*dd,1,dd,D)); // Ddd 1 dd D
  tmp=C;tmp.reshape(Indices(dd*D*dd,D));
  CN.permute(Indices(2,1)); // put first the index for 2nd MPO
  mwArray edgeR=tmp*CN;edgeR.reshape(Indices(dd,D,dd,D));
  edgeR.permute(Indices(1,4,2,3));edgeR.reshape(Indices(dd*D,D,dd,1)); // ddD D dd 1
  // Now for the adjoint I also need modifications
  mwArray Cd(C);Cd.permute(Indices(3,2,1,4),true); // adjoint of normal C
  mwArray edgeLd(Cd);edgeLd.reshape(Indices(dd,1,D*dd,D));
  mwArray edgeRd(Cd);edgeRd.reshape(Indices(dd,D,dd*D,1));
  Operator CdOp(Cd),COp(C);
  Operator edgeLOp(edgeL),edgeROp(edgeR);
  Operator edgeLdOp(edgeLd),edgeRdOp(edgeRd);
  const Operator* arrL[2]={&edgeLdOp,&edgeLOp};
  const Operator* arrR[2]={&edgeRdOp,&edgeROp};
  const Operator* arrC[2]={&CdOp,&COp};
  //cout<<"About to place operators in MPO"<<endl;
  // And now I set the operators one by one
  mpo.setOp(0,new JoinedOperator(2,arrL),true);
  mpo.setOp(L-1,new JoinedOperator(2,arrR),true);
  if(L>2){
    mpo.setOp(1,new JoinedOperator(2,arrC),true);
    for(int k=2;k<L-1;k++)
      mpo.setOp(k,&mpo.getOp(1),false);
  }
}

void checkParticleNumber(MPO& mpo,int L){

  if(mpo.getLength()!=L)
    mpo.initLength(L);
  // Basic operators
  int nrOps=2;
  int D=2; // bond dimension of the MPO
  int dd=d*d;
  mwArray Z=mwArray(Indices(nrOps,dd,dd));
   // Basic pieces
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);sigZ.reshape(Indices(4,1));
  mwArray sigZT(Indices(d,d),dataz);sigZT.reshape(Indices(1,4));
  
  mwArray sig0=identityMatrix(d);//identity

  mwArray id4 = id2*id2T;
  mwArray Sz = sigZ*sigZT;
  

  Operator id4Op(id4), sigmaZ(Sz);

  for(int d1=0;d1<dd;d1++)
    for(int d2=0;d2<dd;d2++){
      Z.setElement(id4.getElement(Indices(d1,d2)),Indices(0,d1,d2));
      Z.setElement(Sz.getElement(Indices(d1,d2)),Indices(1,d1,d2));
    }
  Z.reshape(Indices(nrOps,dd*dd));
  //cout<<"Finished Z:"<<Z.getDimensions()<<endl;
  // now the coefficients

  mwArray C1(Indices(1,D,nrOps)); // the first site
  mwArray C(Indices(D,D,nrOps)); // the middle site
  mwArray CN(Indices(D,1,nrOps)); // the last site


  C.setElement(ONE_c,Indices(0,0,0));
  C1.setElement(ONE_c,Indices(0,0,0));

  C.setElement(ONE_c,Indices(1,1,0));
  CN.setElement(ONE_c,Indices(1,0,0));

  C.setElement(ONE_c,Indices(0,1,1)); 
  C1.setElement(ONE_c,Indices(0,1,1));
  CN.setElement(ONE_c,Indices(0,0,1));


  // Now set the operators in place
  // Now I reshape, contract with Z and set the indices in proper order
  C.reshape(Indices(D*D,nrOps));
  C1.reshape(Indices(1*D,nrOps));
  CN.reshape(Indices(D*1,nrOps));
  C.multiplyRight(Z);C.reshape(Indices(D,D,dd,dd));C.permute(Indices(3,1,4,2));
  C1.multiplyRight(Z);C1.reshape(Indices(1,D,dd,dd));C1.permute(Indices(3,1,4,2));
  CN.multiplyRight(Z);CN.reshape(Indices(D,1,dd,dd));CN.permute(Indices(3,1,4,2));

  Operator center(C);
  Operator first(C1);
  Operator last(CN);
  //cout<<"About to place operators in MPO"<<endl;
  // And now I set the operators one by one
  mpo.setOp(0,first,true);
  mpo.setOp(L-1,last,true);
  if(L>2){
    mpo.setOp(1,center,true);
    for(int k=2;k<L-1;k++)
      mpo.setOp(k,&mpo.getOp(1),false);
  }
}
