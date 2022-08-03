#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "JoinedOperator.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0


using namespace shrt;
using namespace std;

/** Compute numerically the smallest lambda for an MPO of given bond dimensions and range M */

//Construct Basic MPO, namely HxId-IdxH^T
//void constructCommutatorIsing(MPO& mpo,int L,double g,double h);
// Construct the real MPo that I need to use in the GS search, namely the square of the commutator, with the edge sites traced out
void constructAdjointMPO(MPO& HmpoD,const MPO& Hmpo);

// Construct the MPO for the double commutator, directly.
void constructDoubleCommutatorIsing(MPO& mpo,int L,double g,double h);


// Copied from Liouvillian, should be somewhere else
void constructOperatorProduct(mwArray& result,const mwArray& opA,
			      const mwArray& opB);
// Split it in two terms
void split2term(mwArray& Ol,mwArray& Or,mwArray& U,int& nr);

void formEvenDoubleMPO(MPO& mpo,int L,const mwArray& Olv,const mwArray& Orv,bool left=false);
void formOddDoubleMPO(MPO& mpo,int L,const mwArray& Olw,const mwArray& Orw,bool left=false);


void initializeState(MPS& state,bool dagger=false);
void initializeTracelessHermitianState(MPS& state);
void extendMPS(MPS& state);
void compressMPS(MPS& state);
void extendMPSwithId(MPS& state);

// substracts the trace by increasing the bond dimension in 1
// If I am not allowed to increase D, I will need to truncate back to the old one
void substractTrace(MPS& Omps);
complex_t computeTrace(const MPS& Omps);
void forceHermitian(MPS& Omps);

const string jobfilename(int M,int D,double g,double h,const string mpsdir);
const string mpsfilename(int M,int D,const string mpsdir);

void adapt(MPO& Hmpo,MPO& HmpoD);

#define MAXITER 100
#define MAXPERT 10

//#define ISINGCASE 1

int d=2;
//double gP=0.9045;
//double hP=0.8090;
double tau=0.8;

bool traceless=true;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  //  string directory="";
  // if(argc>2){
  //   const char* indir=argv[++cntr];
  //   directory=string(indir)+"/";
  // }
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  // Store the command line, for this is what should be relaunched if not all levels are done
  stringstream commandstr;
  for(int c=0;c<argc;c++) commandstr<<argv[c]<<" ";


  int M=props.getIntProperty("M");
  int D=props.getIntProperty("D");
  //  const char* outfnameG=argv[++cntr];
  // const char* outfnameL=argv[++cntr];
  string mpsdir=props.getProperty("mpsdir");
  string outfile=props.getProperty("outfile");
  int initD=props.getIntProperty("initD");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-5;
  double gP=props.getDoubleProperty("gP");
  double hP=props.getDoubleProperty("hP");
  double pen=props.getDoubleProperty("penalty");
  if(pen<0) pen=10.;
  string jobsdir=props.getProperty("jobsdir"); 
  double tolD=props.getDoubleProperty("tolD");
  if(tolD<0) tol=1E-4;
  //  int lastD=props.getIntProperty("lastD");
  double noise=props.getDoubleProperty("noise");
  if(noise<0) noise=0.;

  string outputMPS=mpsfilename(M,D,mpsdir);
  string inputMPS=mpsfilename(M,initD,mpsdir);
  if(initD<0) inputMPS="";

  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  cout<<"Initialized Contractor"<<endl;
  //  srandom(rSeed);
  
  ofstream* out;
  out=new ofstream(outfile.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output"<<endl;
    exit(1);
  }
  out->close();delete out;
  
  // *out<<"% D="<<D<<", M="<<M<<", rSeed="<<rSeed<<endl;

  // 1) Construct the MPO

  //  MPO Hmpo(M+2),HmpoD(M+2); // Special MPOs with traced out last site (on the input side)
  //  constructCommutatorIsing(Hmpo,M+2,gP,hP);
  //constructAdjointMPO(HmpoD,Hmpo);
  //adapt(Hmpo,HmpoD);


  MPO H2(M);
  //  const MPO* mpos[2]={&HmpoD,&Hmpo};
  //  MPO::join(2,mpos,H2);
  constructDoubleCommutatorIsing(H2,M,gP,hP);


  cout<<"constructed all MPOs:"<<endl;
  // cout<<"Hmpo "<<Hmpo<<endl;
  // cout<<"HmpoD "<<HmpoD<<endl;
  // Hmpo.exportForMatlab("Hmpo.m");
  // HmpoD.exportForMatlab("HmpoD.m");
  H2.exportForMatlab("H2.m");

  //Initialize the MPO O with random entries or existing MPS file
  //vector<int> dims(M+2,d*d);dims[0]=1;dims[M+1]=1;
  //  MPS Omps(M+2,D,dims);
  MPS Omps(M,D,d*d);
  cout<<"Initialized vector Omps="<<Omps<<endl;
  const string tmpfilename=outputMPS+".tmp";
  bool usingTmp=false;
  if(file_exists(tmpfilename)){
    cout<<"There is a tmp file for same D: "<<tmpfilename<<endl;
    Omps.importMPS(tmpfilename.data());
    usingTmp=true;
    cout<<"Imported tmp file for same D!"<<endl;
  }
  else{  // no tmp file
    cout<<"There is no tmp file for same D: "<<tmpfilename<<endl;
    if(inputMPS.empty()){ // no init file indicated in arguments
      cout<<"There is no input file for initD "<<endl;
      substractTrace(Omps);
      cout<<"substracted trace from random O"<<endl;
      //      Omps.setRandomState();
    }
    else{ // initMPS file read
      if(!file_exists(inputMPS)){
	cout<<"ERROR: Initial MPS not found in "<<inputMPS<<endl;
	exit(1);
      }
      Omps.importMPS(inputMPS.data());
      cout<<"Imported initMPS file for smaller D"<<endl;
    }
  }
  Omps.gaugeCond('R',1);Omps.gaugeCond('L',1); 
  cout<<"Imposed the gauge condition on Omps"<<endl;
  cout<<"Now going to contract with H2="<<H2<<endl;
  double oldLambda=real(contractor.contract(Omps,H2,Omps));
  cout<<"constructed (or recovered) initial MPS with lambda="<<oldLambda<<", trace "<<computeTrace(Omps)<<" and <O|O>="<<contractor.contract(Omps,Omps)<<endl;
  
  if(Omps.getBond()<D){
    forceHermitian(Omps);    
    Omps.increaseBondDimensionWithNoise(D,noise);
    substractTrace(Omps);
  }
  // Now apply the minimization with a penalty if trace is not zero.
  MPS idPen(M,1,d*d);
  for(int k=0;k<M;k++)
    idPen.setA(k,reshape(identityMatrix(d),Indices(d*d,1,1)));
  idPen.gaugeCond('R',0);
  // idPen.exportForMatlab("vId.m");
  vector<MPS*> orthog;orthog.push_back(&idPen);
  vector<double> penalty(1,-1.*pen);
  double lambda=real(contractor.contract(Omps,H2,Omps));
  contractor.setEigenSolver(primme);
  cout<<setprecision(15);
  contractor.findGroundStateWithProjectorPenalty(H2,D,penalty,orthog,&lambda,Omps,0.,1,tmpfilename);

  //  contractor.findGroundState(Htot,D,&lambda,Omps);
  cout<<"Found  eigenval="<<lambda<<endl;
  cout<<"Trace is "<<computeTrace(Omps)
      <<" and <O|O>="<<contractor.contract(Omps,Omps)<<endl;
  
  Omps.exportMPS(outputMPS.data());

  // Now, if it was saved, I can delete the temporary file
  if(file_exists(tmpfilename)){
    remove(tmpfilename.data());
  }

  out=new ofstream(outfile.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(15);
  *out<<D<<"\t"<<lambda<<"\t"<<M<<"\t"<<computeTrace(Omps)<<endl;
  out->close();
  delete out;

  if(usingTmp||abs((lambda-oldLambda)/(lambda))>=tolD&&D<260){ // not converged=> lauch again with increased D
    int newD=D+20; // in general
    // other cases
    if(D<20){
      if(D==1) newD=2;
      if(D==2) newD=4;
      if(D==4) newD=10;
      if(D==10) newD=D+10;
    }
    stringstream commandstr;
    commandstr<<argv[0];
    commandstr<<" "<<argv[1];
    commandstr<<" -gP="<<gP<<" -hP="<<hP;
    commandstr<<" -M="<<M<<" -outfile="<<outfile;
    commandstr<<" -jobsdir="<<jobsdir;
    commandstr<<" -noise="<<noise<<" -initD="<<D;
    commandstr<<" -penalty="<<pen<<" -tolD="<<tolD;
    commandstr<<" -tol="<<1E-6<<" -D="<<newD;
    commandstr<<" -mpsdir="<<mpsdir;

    const string jobfile=jobfilename(M,newD,gP,hP,jobsdir);
    ofstream* ojob=new ofstream(jobfile.data());
    if(!ojob->is_open()){
      cout<<"WARNING!!!! Could not open the JOB file "<<jobfile<<". Launch by hand with the following command:"<<endl;
      cout<<commandstr.str()<<endl;
      cout<<endl;
      // Do not exit, because I still need to save the MPS
    }
    else{
      *ojob<<"#!/bin/bash"<<endl;
      *ojob<<"msub_modified_tqo097 ";
      //if(D+20>60) *ojob<<" -P 4";
      *ojob<<" -N slowH."<<hP<<"."<<M<<"."<<newD;
      *ojob<<" -raw ";
      *ojob<<commandstr.str()<<endl;
      ojob->close();
    }
    delete ojob;
  }
  else{
    cout<<"Result converged within "<<tolD<<" at D="<<D<<" or D too large already"<<endl;
  }

}

const string mpsfilename(int M,int D,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS"<<"_M"<<M<<"_D"<<D<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

const string jobfilename(int M,int D,double gP,double hP,const string jobsdir){
  stringstream s;
  s<<jobsdir<<"/job_H_M"<<M<<"_D"<<D<<"_g"<<gP<<"_h"<<hP;  
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

void constructDoubleCommutatorIsing(MPO& mpo,int L,double g,double h){
  cout<<"constructDoubleCommutatorIsing(g="<<g<<",h="<<h<<")"<<endl;
  if(mpo.getLength()!=L)
    mpo.initLength(L);
  // Basic operators
  int nrOps=5;
  int dd=d*d;
  mwArray Z=mwArray(Indices(nrOps,dd,dd));
   // Basic pieces
  mwArray sig0=identityMatrix(d);//identity
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),datax);//sigmax
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz

  // *** First of all, we need the identity on both
  mwArray term0=identityMatrix(dd);
  // *** Now the ones appearing in single-body terms
  // (1) sigX (x) Id
  mwArray term1;
  constructOperatorProduct(term1,sigX,sig0);
  // (2)  Id (x) sigX -> just a permutation of the previous one
  mwArray term2;
  constructOperatorProduct(term2,sig0,permute(sigX,Indices(2,1)));
  // (3) sigZ (x) Id
  mwArray term3;
  constructOperatorProduct(term3,sigZ,sig0);
  // (4)  Id (x) sigZ^T -> just a permutation of the previous one
  mwArray term4;
  constructOperatorProduct(term4,sig0,permute(sigZ,Indices(2,1)));
  // Fill in the operators in the Z array, to be contracted with C in
  // the construction of the MPO. The order is as listed above.
  for(int d1=0;d1<dd;d1++)
    for(int d2=0;d2<dd;d2++){
      Z.setElement(term0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
      Z.setElement(term1.getElement(Indices(d1,d2)),Indices(1,d1,d2));
      Z.setElement(term2.getElement(Indices(d1,d2)),Indices(2,d1,d2));
      Z.setElement(term3.getElement(Indices(d1,d2)),Indices(3,d1,d2));
      Z.setElement(term4.getElement(Indices(d1,d2)),Indices(4,d1,d2));
    }
  Z.reshape(Indices(nrOps,dd*dd));
  //cout<<"Finished Z:"<<Z.getDimensions()<<endl;
  // now the coefficients
  int D=4; // bond dimension of the MPO
  mwArray C1(Indices(1,D,nrOps)); // the first site
  mwArray C(Indices(D,D,nrOps)); // the middle site
  mwArray CN(Indices(D,1,nrOps)); // the last site

  // Identity terms when nothing has happened yet
  C.setElement(ONE_c,Indices(0,0,0));
  C1.setElement(ONE_c,Indices(0,0,0));
  // Identity terms after the end
  C.setElement(ONE_c,Indices(D-1,D-1,0));
  CN.setElement(ONE_c,Indices(D-1,0,0));
  // Single body terms (one for each operator, with proper weights)
  C.setElement(h*ONE_c,Indices(0,D-1,3)); // h sigZ x I
  C.setElement(-h*ONE_c,Indices(0,D-1,4)); // -h I x sigZ^T
  C.setElement(g*ONE_c,Indices(0,D-1,1)); // g sigX x I
  C.setElement(-g*ONE_c,Indices(0,D-1,2)); // -g I x sigX^T
  C1.setElement(h*ONE_c,Indices(0,D-1,3)); // h sigZ x I
  C1.setElement(-h*ONE_c,Indices(0,D-1,4)); // -h I x sigZ^T
  C1.setElement(g*ONE_c,Indices(0,D-1,1)); // g sigX x I
  C1.setElement(-g*ONE_c,Indices(0,D-1,2)); // -g I x sigX^T
  CN.setElement(h*ONE_c,Indices(0,0,3)); // h sigZ x I
  CN.setElement(-h*ONE_c,Indices(0,0,4)); // -h I x sigZ^T
  CN.setElement(g*ONE_c,Indices(0,0,1)); // g sigX x I
  CN.setElement(-g*ONE_c,Indices(0,0,2)); // -g I x sigX^T
  //  Two body terms from H and H^T on the other side
  C.setElement(ONE_c,Indices(0,1,3)); //  sigZ x I
  C.setElement(ONE_c,Indices(1,D-1,3)); 
  C.setElement(-1.*ONE_c,Indices(0,2,4));  //  I x sigZ 
  C.setElement(ONE_c,Indices(2,D-1,4));
  C1.setElement(ONE_c,Indices(0,1,3)); //  sigZ x I
  C1.setElement(-1.*ONE_c,Indices(0,2,4));  //  I x sigZ 
  CN.setElement(ONE_c,Indices(1,0,3)); 
  CN.setElement(ONE_c,Indices(2,0,4));

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

void constructAdjointMPO(MPO& HmpoD,const MPO& Hmpo){
  // I just know that the only new Ops are first, second and last
  int L=Hmpo.getLength();
  HmpoD.initLength(L);
  if(L<2){
    cout<<"ERROR: Cannot work with MPO of length smaller than 2"<<endl;
    exit(1);
  }
  mwArray C1=Hmpo.getOpData(0);
  mwArray C=Hmpo.getOpData(1);
  mwArray CN=Hmpo.getOpData(L-1);
  // Transpose and conjugate them and create new Ops
  Operator* Op1=new Operator(C1,Indices(3,2,1,4),true);
  Operator* OpN=new Operator(CN,Indices(3,2,1,4),true);
  HmpoD.setOp(0,Op1,true);
  HmpoD.setOp(L-1,OpN,true);
  if(L>2){
    Operator* Op=new Operator(C,Indices(3,2,1,4),true);
    HmpoD.setOp(1,Op,true);
    for(int k=2;k<L-1;k++)
      HmpoD.setOp(k,Op,false);
  }
}

void adapt(MPO& Hmpo,MPO& HmpoD){
  int L=Hmpo.getLength();
  if(L<4){
    cout<<"ERROR: I cannot do this trick with only one site left!"<<endl;
    exit(1);
  }
  mwArray C1=Hmpo.getOpData(0);Indices dimsC=C1.getDimensions();
  mwArray C1d=HmpoD.getOpData(0);Indices dimsCd=C1d.getDimensions();
  C1d.permute(Indices(1,2,4,3));C1d.reshape(Indices(1*1*dimsCd[3],dimsCd[2])); 
  C1.reshape(Indices(dimsC[0],1*1*dimsC[3]));
  C1.multiplyLeft(C1d);
  C1.reshape(Indices(dimsCd[3],dimsC[3]));
  mwArray edgeL=Hmpo.getOpData(1);Indices dims=edgeL.getDimensions();
  edgeL.permute(Indices(2,1,3,4));edgeL.reshape(Indices(dims[1],dims[0]*dims[2]*dims[3]));


// first, I transfer the valid pointers
  

}


void initializeState(MPS& state,bool dagger){
  // Initialize the operator O, which is the projector on |0>
  mwArray Op(Indices(d,d));Op.fillWithZero();
  Op.setElement(ONE_c,Indices(0,0));
  if(dagger){
    Op.conjugate();
  }
  //  Op=identityMatrix(d);
  Op.reshape(Indices(d*d,1,1));
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  state=MPS(2,1,d*d);
  state.setA(0,Op);
  state.setA(1,(1./sqrt(2))*id2);
  //state.setA(1,id2);
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized state with norm "<<contractor.contract(state,state)<<endl;
}

void initializeTracelessState(MPS& state,bool dagger){
  // Initialize the operator O, which is the projector on |0>
  mwArray Op(Indices(d,d));Op.fillWithZero();
  Op.setElement(ONE_c,Indices(0,0));
  Op.setElement(-1.*ONE_c,Indices(1,1));
  if(dagger){
    Op.conjugate();
  }
  //  Op=identityMatrix(d);
  Op.reshape(Indices(d*d,1,1));
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  state=MPS(2,1,d*d);
  state.setA(0,Op);
  state.setA(1,(1./sqrt(2))*id2);
  //state.setA(1,id2);
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized state with norm "<<contractor.contract(state,state)<<endl;
}

void extendMPS(MPS& state){
  int L=state.getLength();
  MPS aux(state);
  state=MPS(L+2,state.getBond(),d*d);
  mwArray id2(ONE_c);id2.reshape(Indices(1,1,1));
  state.replaceSite(0,id2,false);
  state.replaceSite(L+1,id2,false);
  for(int k=0;k<L;k++)
    state.replaceSite(k+1,aux.getA(k),false);
  state.setNormFact(aux.getNormFact());
  //  cout<<"Extended MPS to "<<state<<endl;
}

void compressMPS(MPS& state){
  int L=state.getLength();
  if(state.getA(0).getDr()>1||state.getA(L-1).getDl()>1){
    cout<<"Cannot remove edges from MPS "<<state<<endl;
    exit(1);
  }
  MPS aux(state);
  state=MPS(L-2,state.getBond(),d*d);
  for(int k=0;k<L-2;k++)
    state.replaceSite(k,aux.getA(k+1),false);
  double oldNorm=aux.getNormFact();
  mwArray auxA=aux.getA(0).getA();
  if(auxA.getDimension(0)*auxA.getDimension(1)*auxA.getDimension(2)!=1){
    cout<<"Error: wrong dimensions of first site in compressMPS: A="<<auxA<<endl;
    exit(1);
  }
  oldNorm*=abs(auxA.getElement(0));
  auxA=aux.getA(L-1).getA();
  if(auxA.getDimension(0)*auxA.getDimension(1)*auxA.getDimension(2)!=1){
    cout<<"Error: wrong dimensions of last site in compressMPS: A="<<auxA<<endl;
    exit(1);
  }  
  oldNorm*=abs(auxA.getElement(0));
  state.setNormFact(oldNorm);
  // TODO: Check that sites 0 and L-1 are 1x1 and multiply also their values!
  //  cout<<"Extended MPS to "<<state<<endl;
}


void extendMPSwithId(MPS& state){
  int L=state.getLength();
  MPS aux(state);
  state=MPS(L+2,state.getBond(),d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  state.replaceSite(0,(1./sqrt(2))*id2,false);
  state.replaceSite(L+1,(1./sqrt(2))*id2,false);
  for(int k=0;k<L;k++)
    state.replaceSite(k+1,aux.getA(k),false);
  state.setNormFact(aux.getNormFact());
  //  cout<<"Extended MPS to "<<state<<endl;
}

void formEvenDoubleMPO(MPO& mpo,int L,const mwArray& Olv,const mwArray& Orv,bool left){
  // On the first and last sites, I have to contract the op with identity
  if(L%2!=0){
    cout<<"ERROR: Odd M not yet supported!"<<endl;
    exit(1);
  }
  Operator *edgeL,*edgeR;
  Operator *midL,*midR;
  //cout<<"Initializing double (even) op with Ol="<<Olv<<" and Or="<<Orv<<endl;
  // Prepare the ops the first time
  // Double operators for the center
  midL=new DoubleOperator(Olv,conjugate(Olv));
  midR=new DoubleOperator(Orv,conjugate(Orv));
  mwArray id2=(1./sqrt(2))*identityMatrix(d);id2.reshape(Indices(d*d,1));
  mwArray aux;
  Indices dims; 
  if(!left){
    aux=midL->getFullData();
    dims=aux.getDimensions(); // dd,1,dd,Dr2
    aux.permute(Indices(1,2,4,3));aux.reshape(Indices(dims[0]*dims[1]*dims[3],dims[2])); // dd*1*Dr2,dd
    aux.multiplyRight(id2); 
    aux.reshape(Indices(dims[0],dims[1],1,dims[3])); // dd,1,1,Dr2
    edgeL=new Operator(aux);
    aux=midR->getFullData();dims=aux.getDimensions(); // dd,Dl2,dd,1
    //cout<<"midR has dims "<<dims<<endl;
    aux.reshape(Indices(dims[0]*dims[1],dims[2])); // dd*Dl2,dd(*1)
    aux.multiplyRight(id2); 
    aux.reshape(Indices(dims[0],dims[1],1,dims[3])); // dd,Dl2,1,1
    edgeR=new Operator(aux);
  }
  else{
    aux=midL->getFullData();dims=aux.getDimensions(); // dd,1,dd,Dr2
    aux.reshape(Indices(dims[0],dims[1]*dims[2]*dims[3])); // dd,1*dd*Dr2
    id2.reshape(Indices(1,d*d));
    aux.multiplyLeft(id2); // 1,1*dd*Dr2
    aux.reshape(Indices(1,dims[1],dims[2],dims[3])); // 1,1,dd,Dr2
    edgeL=new Operator(aux);
    aux=midR->getFullData();dims=aux.getDimensions(); // dd,Dl2,dd,1
    //cout<<"midR has dims "<<dims<<endl;
    aux.reshape(Indices(dims[0],dims[1]*dims[2]*dims[3])); // dd,Dl2*dd*1
    aux.multiplyLeft(id2); // 1,Dl2*dd*1 
    aux.reshape(Indices(1,dims[1],dims[2],dims[3])); // 1,Dl2,dd,1
    edgeR=new Operator(aux);
  }
  // Set the operators
  mpo.initLength(L);
  if(L==2){
    mpo.setOp(0,midL,true);
    mpo.setOp(L-1,midR,true);
  }
  else{
    mpo.setOp(0,edgeL,true);
    mpo.setOp(L-1,edgeR,true);
  }
  int k=1;
  bool saved=false;
  while(k<L-1){
    if(!saved){
      mpo.setOp(k++,midR,true);
      mpo.setOp(k++,midL,true);    
      saved=true;
    }
    else{
      mpo.setOp(k++,midR,false);
      mpo.setOp(k++,midL,false);    
    }
  }
}

void formOddDoubleMPO(MPO& mpo,int L,const mwArray& Olv,const mwArray& Orv,bool left){
  // On the first and last sites, I have to contract the op with identity
  Operator *edgeL,*edgeR;
  Operator *midL,*midR;
  // Double operators for the center
  midL=new DoubleOperator(Olv,conjugate(Olv));
  midR=new DoubleOperator(Orv,conjugate(Orv));
  mwArray id2=(1./sqrt(2))*identityMatrix(d);id2.reshape(Indices(d*d,1));

  mwArray aux=midL->getFullData();
  Indices dims=aux.getDimensions(); // dd,1,dd,Dr2
  aux.permute(Indices(1,2,4,3));aux.reshape(Indices(dims[0]*dims[1]*dims[3],dims[2])); // dd*1*Dr2,dd
  aux.multiplyRight(id2); 
  aux.reshape(Indices(dims[0],dims[1]*1*dims[3])); // dd,1*1*Dr2
  id2.reshape(Indices(1,d*d));
  aux.multiplyLeft(id2); // 1,Dr2
  mwArray aux2=midR->getFullData();  // dd,Dl2,dd,1
  Indices dims2=aux2.getDimensions(); // dd,Dl2,dd,1
  if(!left){ // extra id2 for traicing out the extra site
    aux2.reshape(Indices(dims2[0],dims2[1]*dims2[2])); // dd,Dl2*dd(*1)
    aux2.multiplyLeft(id2); 
    aux2.reshape(Indices(dims2[1],dims2[2]*dims2[3])); // Dl2,dd*1
  }
  else{ // the trace is applied on the down indices
    aux2.reshape(Indices(dims2[0]*dims2[1],dims2[2])); // dd,Dl2*dd(*1)
    id2.reshape(Indices(d*d,1));
    aux2.multiplyRight(id2);
    aux2.reshape(Indices(dims2[0],dims2[1])); // dd,Dl2
    aux2.permute(Indices(2,1));
  }
  aux2.reshape(Indices(dims2[1],-1)); // Dl2,dd*1
  aux2.multiplyLeft(aux); // 1,dd*1
  if(!left)
    aux2.reshape(Indices(1,1,dims2[2],1)); // 1,1,dd,1
  else
    aux2.reshape(Indices(dims2[0],1,1,1)); // dd,1,1,1
  edgeL=new Operator(aux2);

  // Now the same for right side
  aux=midR->getFullData();dims=aux.getDimensions(); // dd,Dl2,dd,1
  aux.reshape(Indices(dims[0]*dims[1],dims[2])); // dd*Dl2,dd(*1)
  id2.reshape(Indices(d*d,1));
  aux.multiplyRight(id2); 
  aux.reshape(Indices(dims[0],dims[1]*1*dims[3])); // dd,Dl2*1*1
  id2.reshape(Indices(1,d*d));aux.multiplyLeft(id2);  //1,Dl2
  aux.reshape(Indices(dims[1],1));
  aux2=midL->getFullData(); dims2=aux2.getDimensions(); // dd,1,dd,Dr2
  if(!left){
    aux2.reshape(Indices(dims2[0],dims2[1]*dims2[2]*dims2[3])); // dd,1*dd*Dr2
    aux2.multiplyLeft(id2);  // 1,1*dd*Dr2
    aux2.reshape(Indices(dims2[2],dims2[3])); // 1*dd,Dr2
  }
  else{
    aux2.permute(Indices(3,1,2,4));
    aux2.reshape(Indices(dims2[2],dims2[0]*dims2[1]*dims2[3])); // dd,dd*1*Dr2
    aux2.multiplyLeft(id2);  // 1,1*dd*Dr2
    aux2.reshape(Indices(dims2[0],dims2[3])); // dd*1,Dr2
  }
  aux2.multiplyRight(aux); // dd,
  if(!left)
    aux2.reshape(Indices(1,1,dims[2],1));
  else
    aux2.reshape(Indices(dims[0],1,1,1));
  edgeR=new Operator(aux2);

  // Set the operators
  mpo.initLength(L);
  mpo.setOp(0,edgeL,true);
  mpo.setOp(L-1,edgeR,true);
  int k=1;bool stored=false;
  while(k<L-1){
    if(!stored){
      mpo.setOp(k++,midL,true);
      mpo.setOp(k++,midR,true);
      stored=true;    
    }
    else{
      mpo.setOp(k++,midL,false);
      mpo.setOp(k++,midR,false);    
    }
  }
  //cout<<"Initialized formOddDouble "<<mpo<<endl;
}

// void substractTrace(MPS& Omps){
//   int L=Omps.getLength();
//   // First compute trace
//   complex_t tr=computeTrace(Omps);
//   if(abs(imag(tr))>1E-10&&abs(imag(tr)/real(tr))>1E-3){
//     cout<<"Apparently there is a  problem, as the trace is not real! tr(O)="<<tr<<endl;
//     //exit(1);
//   }
//   double signTr=real(tr)>=0?1:-1;
//   double val=(1./d)*pow(abs(real(tr)),1./(L-2));
//   for(int k=1;k<L-1;k++){
//     mwArray A=Omps.getA(k).getA();
//     Indices dims=A.getDimensions(); // dd,Dl,Dr
//     int newDl=k==1?dims[1]:dims[1]+1;
//     int newDr=k==L-2?dims[2]:dims[2]+1;
//     A.resize(Indices(d*d,newDl,newDr));
//     A.reshape(Indices(d,d,newDl,newDr));
//     double sign=k==0?-1.*signTr:1.;
//     for(int j=0;j<d;j++)
//       A.setElement(sign*val*ONE_c,Indices(j,j,newDl-1,newDr-1)); 
//     A.reshape(Indices(d*d,newDl,newDr));
//     Omps.replaceSite(k,A,false);
//   }
// }

void substractTrace(MPS& Omps){
  int L=Omps.getLength();
  // First compute trace
  complex_t tr=computeTrace(Omps);
  cout<<"In substractTrace, trace is "<<tr<<endl;
  if(abs(imag(tr))>1E-10&&abs(imag(tr)/real(tr))>1E-3){
    cout<<"Apparently there is a  problem, as the trace is not real! tr(O)="<<tr<<endl;
    //exit(1);
  }
  double signTr=real(tr)>=0?1:-1;
  double val=(1./d)*pow(abs(real(tr)),1./L);
  for(int k=0;k<L;k++){
    mwArray A=Omps.getA(k).getA();
    Indices dims=A.getDimensions(); // dd,Dl,Dr
    int newDl=k==0?dims[1]:dims[1]+1;
    int newDr=k==L-1?dims[2]:dims[2]+1;
    A.resize(Indices(d*d,newDl,newDr));
    A.reshape(Indices(d,d,newDl,newDr));
    //    cout<<"Site "<<k<<" going from dims="<<dims<<" to "<<A.getDimensions()<<endl;
    double sign=k==0?-1.*signTr:1.;
    for(int j=0;j<d;j++)
      A.setElement(sign*val*ONE_c,Indices(j,j,newDl-1,newDr-1)); 
    A.reshape(Indices(d*d,newDl,newDr));
    Omps.replaceSite(k,A,false);
  }
}

// complex_t computeTrace(const MPS& Omps){
//   mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
//   int len=Omps.getLength();
//   MPS aux(len,1,1);
//   for(int k=1;k<len-1;k++){
//     aux.replaceSite(k,id2,false);
//   }
//   Contractor& contractor=Contractor::theContractor();
//   return contractor.contract(Omps,aux);
// }
 
complex_t computeTrace(const MPS& Omps){
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  int len=Omps.getLength();
  MPS aux(len,1,d*d);
  for(int k=0;k<len;k++){
    aux.setA(k,id2);
  }
  Contractor& contractor=Contractor::theContractor();
  return contractor.contract(Omps,aux);
}
 
void forceHermitian(MPS& Omps){
  int len=Omps.getLength();
  for(int k=0;k<len;k++){
    mwArray A=Omps.getA(k).getA();
    Indices dims=A.getDimensions(); // dd,Dl,Dr
    A.reshape(Indices(d,d,dims[1],dims[2]));
    mwArray B(A);
    B.permute(Indices(2,1,3,4),true);
    A=.5*(A+B);    
    A.reshape(Indices(d*d,dims[1],dims[2]));
    Omps.setA(k,A);
  }
}
