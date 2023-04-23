#include <math.h>
#include <iomanip>
#include <fstream>
#include <deque>

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

/** Compute numerically the smallest lambda for an MPO of given bond dimensions and range M */

//Construct Basic MPO, namely HxId-IdxH^T



// Construct the MPO for the double commutator, directly.
void constructDoubleCommutatorHeis(MPO& mpo,int L,double h, mwArray instances, int rowInstances, vector<double> &randFields);
void constructDoubleCommutatorHeisembergAncilla(MPO& mpo,int L,double J1x,double J1y,double J2x, double J2y, double J2z,double Bond);

// Copied from Liouvillian, should be somewhere else
void constructOperatorProduct(mwArray& result,const mwArray& opA,
			      const mwArray& opB);

// substracts the trace by increasing the bond dimension in 1
// If I am not allowed to increase D, I will need to truncate back to the old one
void substractTrace(MPS& Omps);
complex_t computeTrace(const MPS& Omps);
void forceHermitian(MPS& Omps);

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
  int maxD=props.getIntProperty("maxD");
  int rowInstances=props.getIntProperty("rowInstances");
  double Bond=props.getDoubleProperty("hP");
  string mpsdir=props.getProperty("mpsdir");
  //string mpsdirin=props.getProperty("mpsdirin");
  string inf=props.getProperty("instances");
  string outfile=props.getProperty("outfile");
  int initD=props.getIntProperty("initD");
  double threshold=props.getDoubleProperty("threshold");
  double numrounds=props.getDoubleProperty("numrounds");
  double lengthinterval=props.getDoubleProperty("lengthinterval");
  double MaxSteps=props.getDoubleProperty("MaxSteps");
  double tol=props.getDoubleProperty("tol");
  double tolEis=props.getDoubleProperty("tolEis");
  double tolZero=props.getDoubleProperty("tolZero");
  double tolLong=props.getDoubleProperty("tolLong");
  if(tol<0) tol=1E-4;
  double tolD=props.getDoubleProperty("tolD");
  if(tolD<0) tolD=1E-4;
  double pen=props.getDoubleProperty("penalty");
  if(pen<0) pen=10.;

  double noise=props.getDoubleProperty("noise");
  if(noise<0) noise=0.;

  string outputMPS=mpsfilename(M,D,mpsdir);   
  //string inputMPS=mpsfilename(M,initD,mpsdirin);
  string inputMPS=mpsfilename(M,initD,mpsdir);
  
  const string tmpoutfile = outputMPS+".tmp";
  const string tmpinfile = inputMPS+".tmp";
  if(initD<0) inputMPS="";
 
  //Import random instances of the hamiltonian
  const char* infname = inf.c_str();
  mwArray instances;
  ifstream in(infname);
  if(!in.is_open()){
    cout<<"Error: couldn't open file "<<infname<<" for output"<<endl;
    exit(1);
  }
  instances.load(in);
  in.close();


  //Initailize contractor

  Contractor& contractor=Contractor::theContractor(); //inizialize contractor
  contractor.setConvTol(tol);
  
  contractor.setEigTol(tolEis);
  contractor.setZeroTol(tolZero);
  contractor.setLongTol(tolLong);

  cout<<"Parameters: M="<<M<<", h="<<hP<<", D="<<D<<endl;
  cout<<"eigtol: "<<contractor.getEigTol()<<endl;
  cout<<"zerotol: "<<contractor.getZeroTol()<<endl;
  cout<<"long tolerance: "<<contractor.getConvLong()<<endl;
  cout<<"threshold: "<<threshold<<endl;
  cout<<"tolerance: "<<tol<<endl;
  cout<<"numrounds: "<<numrounds<<endl;
  cout<<"lengthinterval: "<<lengthinterval<<endl;
  


  //Test outfile
  
  ofstream* out;
  out=new ofstream(outfile.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output"<<endl;
    exit(1);
  }
  out->close();delete out;


  MPO H2(2*M);
  //  const MPO* mpos[2]={&HmpoD,&Hmpo};
  //  MPO::join(2,mpos,H2);
  double J1x=1,J1y=1,J2x=1,J2y=1,J2z=1;
  constructDoubleCommutatorHeisembergAncilla(H2,M,J1x,J1y,J2x,J2y,J2z,Bond);
  
  //Define penalities
  MPS idPen(M,1,d*d);
  for(int k=0;k<M;k++)
    idPen.setA(k,reshape(identityMatrix(d),Indices(d*d,1,1)));
  idPen.gaugeCond('R',0);
  vector<MPS*> orthog;orthog.push_back(&idPen);
  vector<double> penalty(1,-1.*pen);
  

  //Set the parameter for convergence:
  //first element: threshold value to reach(the power of 10: 1E-${value here})
  //second element: least number of rounds to make 
  //third element: length of the interval to check
  vector<double> convergenceTools;
  convergenceTools.push_back(threshold);
  convergenceTools.push_back(numrounds);
  convergenceTools.push_back(lengthinterval);
  convergenceTools.push_back(MaxSteps);
  

  //Import or define initial state
  MPS Omps(M,D,d*d);

  if(file_exists(tmpinfile)){
    Omps.importMPS(tmpinfile.data());
    cout<<"Imported tmp file for same D!"<<tmpinfile<<endl;
  }else{  // no tmp file
    cout<<"There is no tmp file for same D: "<<tmpinfile<<endl;
    if(inputMPS.empty()){ // no init file indicated in arguments
      cout<<"There is no input file for initD "<<endl;
      if(D<20){
	Omps.setRandomState();
	cout<<"Random State imported!"<<endl;
      }else{
	cout<<"No initial state imported: probably wrong parameters."<<endl;
	exit (1);
      }
      substractTrace(Omps);
      cout<<"substracted trace from random O"<<endl;

    }
    else{ // initMPS file read
      if(!file_exists(inputMPS)){
	if(D<20){
	  Omps.setRandomState();
	  cout<<"Random State imported!"<<endl;
	}else{
	  cout<<"No initial state imported: probably wrong parameters."<<endl;
	  exit (1);
	}
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

  double lambda=real(contractor.contract(Omps,H2,Omps));
  double oldLambda = lambda; 
  contractor.setEigenSolver(primme);
 
  contractor.findGroundStateWithProjectorPenalty(H2,D,penalty, convergenceTools,orthog,&lambda,Omps,0.,1,tmpoutfile);
 
  Omps.exportMPS(outputMPS.data());
  
  if(file_exists(tmpoutfile)){
    remove(tmpoutfile.data());
  }
  
  out=new ofstream(outfile.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<"for output"<<endl;
    exit(1);
  }
  *out<<setprecision(15);
  *out<<D<<"\t"<<M<<"\t"<<lambda<<"\t"<<computeTrace(Omps)<<endl;
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


// void constructDoubleCommutatorHeis(MPO& mpo,int L,double h, mwArray instances, int rowInstances, vector<double> randFields){

//   if(mpo.getLength()!=L)
//     mpo.initLength(L);
//   // Basic operators
//   double hi=h;
//   int nrOps=7;
//   int dd=d*d;
//   mwArray Z=mwArray(Indices(nrOps,dd,dd));
//   // Basic pieces
//   mwArray sig0=identityMatrix(d);//identity
//   complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
//   mwArray sigX(Indices(d,d),datax);//sigmax
//   complex_t datay[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
//   mwArray sigY(Indices(d,d),datay);//sigmay
//   complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
//   mwArray sigZ(Indices(d,d),dataz);//sigmaz

//   // *** First of all, we need the identity on both
//   mwArray term0=identityMatrix(dd); //id*id
//   // *** Now the ones appearing in single-body terms
//   // (1) sigX (x) Id
//   mwArray term1;
//   constructOperatorProduct(term1,sigX,sig0); //sigX*id
//   // (2)  Id (x) sigX -> just a permutation of the previous one
//   mwArray term2;
//   constructOperatorProduct(term2,sig0,permute(sigX,Indices(2,1))); //id*sigX^T
//   // (3) sigY (x) Id
//   mwArray term3;
//   constructOperatorProduct(term3,sigY,sig0); //sigY*id
//   // (4)  Id (x) sigY^T -> just a permutation of the previous one
//   mwArray term4;
//   constructOperatorProduct(term4,sig0,permute(sigY,Indices(2,1))); //id*sigY^T
//   // (5) sigZ (x) Id
//   mwArray term5;
//   constructOperatorProduct(term5,sigZ,sig0); //sigZ*id
//   // (6)  Id (x) sigZ^T -> just a permutation of the previous one
//   mwArray term6;
//   constructOperatorProduct(term6,sig0,permute(sigZ,Indices(2,1))); //id*sigZ^T
  
//   // Fill in the operators in the Z array, to be contracted with C in
//   // the construction of the MPO. The order is as listed above.
//   for(int d1=0;d1<dd;d1++)
//     for(int d2=0;d2<dd;d2++){
//       Z.setElement(term0.getElement(Indices(d1,d2)),Indices(0,d1,d2)); //id*id
//       Z.setElement(term1.getElement(Indices(d1,d2)),Indices(1,d1,d2)); //sigX*id
//       Z.setElement(term2.getElement(Indices(d1,d2)),Indices(2,d1,d2)); //id*sigX^T
//       Z.setElement(term3.getElement(Indices(d1,d2)),Indices(3,d1,d2)); //sigY*id
//       Z.setElement(term4.getElement(Indices(d1,d2)),Indices(4,d1,d2)); //id*sigY^T
//       Z.setElement(term5.getElement(Indices(d1,d2)),Indices(5,d1,d2)); //sigZ*id
//       Z.setElement(term6.getElement(Indices(d1,d2)),Indices(6,d1,d2)); //id*sigZ^T
//     }
//   Z.reshape(Indices(nrOps,dd*dd));
//   //cout<<"Finished Z:"<<Z.getDimensions()<<endl;
//   // now the coefficients
//   int D=8; // bond dimension of the MPO
//   mwArray C(Indices(D,D,nrOps)); // the middle site
//   mwArray C1(Indices(1,D,nrOps)); // the first site
//   mwArray CN(Indices(D,1,nrOps)); // the last site
  
//   // Single body terms (one for each operator, with proper weights)
//   C1.setElement(ONE_c,Indices(0,0,0));      //  I x I
//   C1.setElement(ONE_c,Indices(0,1,1));      //  sigX x I
//   C1.setElement(-1.*ONE_c,Indices(0,2,2));  //  I x sigX 
//   C1.setElement(ONE_c,Indices(0,3,3));      //  sigY x I
//   C1.setElement(-1.*ONE_c,Indices(0,4,4));  //  I x sigY
//   C1.setElement(ONE_c,Indices(0,5,5));      //  sigZ x I
//   C1.setElement(-1.*ONE_c,Indices(0,6,6));  //  I x sigZ 
  
//   double h1 = h*real(instances.getElement(Indices(rowInstances,L+1)));
//   //  h1 = 1.;
//   C1.setElement(h1*ONE_c,Indices(0,7,5));  // h sigZ x I
//   C1.setElement(-h1*ONE_c,Indices(0,7,6)); // -h I x sigZ^T
//   //cout<<h1<<",";

//   CN.setElement(ONE_c,Indices(D-1,0,0));
//   CN.setElement(ONE_c,Indices(1,0,1)); 
//   CN.setElement(ONE_c,Indices(2,0,2));
//   CN.setElement(ONE_c,Indices(3,0,3)); 
//   CN.setElement(ONE_c,Indices(4,0,4));
//   CN.setElement(ONE_c,Indices(5,0,5)); 
//   CN.setElement(ONE_c,Indices(6,0,6));
  
  
//   double hN = h*real(instances.getElement(Indices(rowInstances,L)));
//   //  hN = 1.;
//   CN.setElement(hN*ONE_c,Indices(0,0,5)); // h sigZ x I
//   CN.setElement(-hN*ONE_c,Indices(0,0,6)); // -h I x sigZ^T
  

//   deque<mwArray> totOps;
//   int ind=0;
//   for(int k=0; k<L; k++){
    
//     hi = h*real(instances.getElement(Indices(rowInstances,k)));
//     //    if (k!=22) hi = 1.;

//     if(k%2 == 0){
//       totOps.push_back(mwArray(Indices(D,D,nrOps)));
//       ind = totOps.size()-1; //cout<<hi<<"\t"<<k<<"\t"<<ind<<endl;  
//       }else{ 
//       totOps.push_front(mwArray(Indices(D,D,nrOps)));
//       ind = 0; //cout<<hi<<"\t"<<k<<"\t"<<ind<<endl;  
//     }

//     totOps[ind].setElement(ONE_c,Indices(0,0,0));      //I x I
//     totOps[ind].setElement(ONE_c,Indices(D-1,D-1,0));  //I x I
    
//     for(int s=1; s<D-1; s++){
//       totOps[ind].setElement(pow(-1.,s+1)*ONE_c,Indices(0,s,s));
//       totOps[ind].setElement(ONE_c,Indices(s,D-1,s)); 
//     }
    
//     totOps[ind].setElement(hi*ONE_c,Indices(0,D-1,5)); // hi sigZ x I
//     totOps[ind].setElement(-hi*ONE_c,Indices(0,D-1,6)); // -hi I x sigZ^T
    
//     totOps[ind].reshape(Indices(D*D,nrOps));
//     totOps[ind].multiplyRight(Z);
//     totOps[ind].reshape(Indices(D,D,dd,dd));
//     totOps[ind].permute(Indices(3,1,4,2));
    
//   }
  
//   double rf;
//   // cout<<endl<<"random magnetic fields well ordered:"<<endl;
//   // cout<<real(C1.getElement(Indices(0,7,5)))<<",";
//   for(int r=1;r<L-1;r++){
//     rf = real(totOps[r].getElement(Indices(2,0,2,D-1)));
//     randFields[r] = (rf/2)*(cos(M_PI*(r)/(2*L)) + cos(M_PI*(r-1)/(2*L)));
//   }
//   rf = real(totOps[0].getElement(Indices(2,0,2,D-1)));
//   randFields[0] = (rf/2)*cos(M_PI/(2*L));
//   rf = real(totOps[L-1].getElement(Indices(2,0,2,D-1)));
//   randFields[L-1] = (rf/2)*cos(M_PI*(L+1)/(2*L));
//   // cout<<real(CN.getElement(Indices(0,0,5)))<<endl<<endl;

//   C1.reshape(Indices(1*D,nrOps));
//   CN.reshape(Indices(D*1,nrOps));
  
//   C1.multiplyRight(Z);C1.reshape(Indices(1*D*dd,dd));
//   // Special situation: first and last sites are contracted with Id as traced out
//   mwArray id2=1/sqrt(2.)*identityMatrix(d);id2.reshape(Indices(dd,1));
//   C1.multiplyRight(id2);
//   C1.reshape(Indices(1,D,dd,1));C1.permute(Indices(3,1,4,2)); //dd x 1 x 1 x D
  
//   CN.multiplyRight(Z);CN.reshape(Indices(D*1*dd,dd));
//   CN.multiplyRight(id2);
//   CN.reshape(Indices(D,1,dd,1));CN.permute(Indices(3,1,4,2)); // dd x D x 1 x 1

//   // Now I build the MPO with required ops.
//   //1) on the edges, I contract C1 and CN with their adjoint
//   C1.reshape(Indices(dd,D));
//   CN.reshape(Indices(dd,D)); // dd D
//   mwArray tmp(C1);tmp.Hconjugate(); 
//   C1.multiplyLeft(tmp); // Dl,Dr
//   tmp=CN;tmp.Hconjugate(); 
//   CN.multiplyLeft(tmp); //Dlu,Dld
  
//   // 2) These, I contract with the middle ones 
//   tmp=totOps[0];tmp.permute(Indices(2,1,3,4));tmp.reshape(Indices(D,dd*dd*D));
//   mwArray edgeL=C1*tmp;edgeL.reshape(Indices(D*dd,1,dd,D)); // Ddd 1 dd D
//   tmp=totOps[L-1];tmp.reshape(Indices(dd*D*dd,D));
//   CN.permute(Indices(2,1)); // put first the index for 2nd MPO
//   mwArray edgeR=tmp*CN;edgeR.reshape(Indices(dd,D,dd,D));
//   edgeR.permute(Indices(1,4,2,3));edgeR.reshape(Indices(dd*D,D,dd,1)); // ddD D dd 1
  
//   // Now for the adjoint I also need modifications
//   mwArray Cdl(totOps[0]);Cdl.permute(Indices(3,2,1,4),true); // adjoint of normal C
//   mwArray Cdr(totOps[L-1]);Cdr.permute(Indices(3,2,1,4),true); // adjoint of normal C
//   mwArray edgeLd(Cdl);edgeLd.reshape(Indices(dd,1,D*dd,D));
//   mwArray edgeRd(Cdr);edgeRd.reshape(Indices(dd,D,dd*D,1));
  
//   Operator edgeLOp(edgeL),edgeROp(edgeR);
//   Operator edgeLdOp(edgeLd),edgeRdOp(edgeRd);
//   const Operator* arrL[2]={&edgeLdOp,&edgeLOp};
//   const Operator* arrR[2]={&edgeRdOp,&edgeROp};

//   vector<Operator*> CdOp,COp;
  
//   mpo.setOp(0,new JoinedOperator(2,arrL),true);
//   mpo.setOp(L-1,new JoinedOperator(2,arrR),true);
//   deque<mwArray> Cd(totOps);
//   //cdp.push_back(new Operator(*(Cd[3])));

  
//   if(L>2){
    
//     for(int k=1;k<L-1;k++){
      
//       Cd[k].permute(Indices(3,2,1,4),true);
//       CdOp.push_back(new Operator(Cd[k]));
//       COp.push_back(new Operator(totOps[k]));
//       const Operator* arrC[2] = {CdOp[k-1],COp[k-1]};
//       mpo.setOp(k,new JoinedOperator(2,arrC),true);  
//     }
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

void constructDoubleCommutatorHeis(MPO& mpo,int L,double h, mwArray instances, int rowInstances, vector<double> &randFields){

  if(mpo.getLength()!=L)
    mpo.initLength(L);
  // Basic operators

  int nrOps=7;
  int dd=d*d;
  mwArray Z=mwArray(Indices(nrOps,dd,dd));
  // Basic pieces
  mwArray sig0=identityMatrix(d);//identity
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),datax);//sigmax
  complex_t datay[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),datay);//sigmay
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz

  // *** First of all, we need the identity on both
  mwArray term0=identityMatrix(dd); //id*id
  // *** Now the ones appearing in single-body terms
  // (1) sigX (x) Id
  mwArray term1;
  constructOperatorProduct(term1,sigX,sig0); //sigX*id
  // (2)  Id (x) sigX -> just a permutation of the previous one
  mwArray term2;
  constructOperatorProduct(term2,sig0,permute(sigX,Indices(2,1))); //id*sigX^T
  // (3) sigY (x) Id
  mwArray term3;
  constructOperatorProduct(term3,sigY,sig0); //sigY*id
  // (4)  Id (x) sigY^T -> just a permutation of the previous one
  mwArray term4;
  constructOperatorProduct(term4,sig0,permute(sigY,Indices(2,1))); //id*sigY^T
  // (5) sigZ (x) Id
  mwArray term5;
  constructOperatorProduct(term5,sigZ,sig0); //sigZ*id
  // (6)  Id (x) sigZ^T -> just a permutation of the previous one
  mwArray term6;
  constructOperatorProduct(term6,sig0,permute(sigZ,Indices(2,1))); //id*sigZ^T
  
  // Fill in the operators in the Z array, to be contracted with C in
  // the construction of the MPO. The order is as listed above.
  for(int d1=0;d1<dd;d1++)
    for(int d2=0;d2<dd;d2++){
      Z.setElement(term0.getElement(Indices(d1,d2)),Indices(0,d1,d2)); //id*id
      Z.setElement(term1.getElement(Indices(d1,d2)),Indices(1,d1,d2)); //sigX*id
      Z.setElement(term2.getElement(Indices(d1,d2)),Indices(2,d1,d2)); //id*sigX^T
      Z.setElement(term3.getElement(Indices(d1,d2)),Indices(3,d1,d2)); //sigY*id
      Z.setElement(term4.getElement(Indices(d1,d2)),Indices(4,d1,d2)); //id*sigY^T
      Z.setElement(term5.getElement(Indices(d1,d2)),Indices(5,d1,d2)); //sigZ*id
      Z.setElement(term6.getElement(Indices(d1,d2)),Indices(6,d1,d2)); //id*sigZ^T
    }
  Z.reshape(Indices(nrOps,dd*dd));
  //cout<<"Finished Z:"<<Z.getDimensions()<<endl;
  // now the coefficients
  int D=8; // bond dimension of the MPO
  mwArray C(Indices(D,D,nrOps)); // the middle site
  mwArray C1(Indices(1,D,nrOps)); // the first site
  mwArray CN(Indices(D,1,nrOps)); // the last site

  vector<mwArray> totOps;
  deque<double> hi;
  for(int k=0; k<=L+1; k++){
    
    double tmphi = h*real(instances.getElement(Indices(rowInstances,k)));
    if(norand){
      hi.push_back(.5);
    }
    if(k%2 == 0 && !norand){
      hi.push_back(tmphi);
      }else{ 
      hi.push_front(tmphi);
    }
    
  }
  
  // Single body terms (one for each operator, with proper weights)
  C1.setElement(ONE_c,Indices(0,0,0));      //  I x I
  C1.setElement(ONE_c,Indices(0,1,1));      //  sigX x I
  C1.setElement(-1.*ONE_c,Indices(0,2,2));  //  I x sigX 
  C1.setElement(ONE_c,Indices(0,3,3));      //  sigY x I
  C1.setElement(-1.*ONE_c,Indices(0,4,4));  //  I x sigY
  C1.setElement(ONE_c,Indices(0,5,5));      //  sigZ x I
  C1.setElement(-1.*ONE_c,Indices(0,6,6));  //  I x sigZ 
  C1.setElement(hi[0]*ONE_c,Indices(0,7,5));   // h1 sigZ x I
  C1.setElement(-hi[0]*ONE_c,Indices(0,7,6));  // -h1 I x sigZ^T

  CN.setElement(ONE_c,Indices(D-1,0,0));    //  I x I
  CN.setElement(ONE_c,Indices(1,0,1));      //  sigX x I
  CN.setElement(ONE_c,Indices(2,0,2));      //  I x sigX
  CN.setElement(ONE_c,Indices(3,0,3));      //  sigY x I 
  CN.setElement(ONE_c,Indices(4,0,4));      //  I x sigY
  CN.setElement(ONE_c,Indices(5,0,5));      //  sigZ x I 
  CN.setElement(ONE_c,Indices(6,0,6));      //  I x sigZ
  CN.setElement(hi[L+1]*ONE_c,Indices(0,0,5));   // hN sigZ x I
  CN.setElement(-hi[L+1]*ONE_c,Indices(0,0,6));  // -hN I x sigZ^T

  
  for(int ind = 0; ind<L; ind++){
    totOps.push_back(mwArray(Indices(D,D,nrOps)));
    totOps[ind].setElement(ONE_c,Indices(0,0,0));      //I x I
    totOps[ind].setElement(ONE_c,Indices(D-1,D-1,0));  //I x I
    
    for(int s=1; s<D-1; s++){
      totOps[ind].setElement(pow(-1.,s+1)*ONE_c,Indices(0,s,s));
      totOps[ind].setElement(ONE_c,Indices(s,D-1,s)); 
    }
    
    totOps[ind].setElement(hi[ind+1]*ONE_c,Indices(0,D-1,5)); // hi sigZ x I
    totOps[ind].setElement(-hi[ind+1]*ONE_c,Indices(0,D-1,6)); // -hi I x sigZ^T
    
    totOps[ind].reshape(Indices(D*D,nrOps));
    totOps[ind].multiplyRight(Z);
    totOps[ind].reshape(Indices(D,D,dd,dd));
    totOps[ind].permute(Indices(3,1,4,2));
    
  }

  for(int r=0;r<=L+2;r++)
    randFields[r] = hi[r];//transform deque into vector
  
  C1.reshape(Indices(1*D,nrOps));
  CN.reshape(Indices(D*1,nrOps));
  
  C1.multiplyRight(Z);C1.reshape(Indices(1*D*dd,dd));
  // Special situation: first and last sites are contracted with Id as traced out
  mwArray id2=1/sqrt(2.)*identityMatrix(d);id2.reshape(Indices(dd,1));
  C1.multiplyRight(id2);
  C1.reshape(Indices(1,D,dd,1));C1.permute(Indices(3,1,4,2)); //dd x 1 x 1 x D
  
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
  tmp=totOps[0];tmp.permute(Indices(2,1,3,4));tmp.reshape(Indices(D,dd*dd*D));
  mwArray edgeL=C1*tmp;edgeL.reshape(Indices(D*dd,1,dd,D)); // Ddd 1 dd D
  tmp=totOps[L-1];tmp.reshape(Indices(dd*D*dd,D));
  CN.permute(Indices(2,1)); // put first the index for 2nd MPO
  mwArray edgeR=tmp*CN;edgeR.reshape(Indices(dd,D,dd,D));
  edgeR.permute(Indices(1,4,2,3));edgeR.reshape(Indices(dd*D,D,dd,1)); // ddD D dd 1
  
  // Now for the adjoint I also need modifications
  mwArray Cdl(totOps[0]);Cdl.permute(Indices(3,2,1,4),true); // adjoint of normal C
  mwArray Cdr(totOps[L-1]);Cdr.permute(Indices(3,2,1,4),true); // adjoint of normal C
  mwArray edgeLd(Cdl);edgeLd.reshape(Indices(dd,1,D*dd,D));
  mwArray edgeRd(Cdr);edgeRd.reshape(Indices(dd,D,dd*D,1));
  
  Operator edgeLOp(edgeL),edgeROp(edgeR);
  Operator edgeLdOp(edgeLd),edgeRdOp(edgeRd);
  const Operator* arrL[2]={&edgeLdOp,&edgeLOp};
  const Operator* arrR[2]={&edgeRdOp,&edgeROp};

  vector<Operator*> CdOp,COp;
  
  mpo.setOp(0,new JoinedOperator(2,arrL),true);
  mpo.setOp(L-1,new JoinedOperator(2,arrR),true);
  vector<mwArray> Cd(totOps);
  //cdp.push_back(new Operator(*(Cd[3])));

  
  if(L>2){
    
    for(int k=1;k<L-1;k++){
      
      Cd[k].permute(Indices(3,2,1,4),true);
      CdOp.push_back(new Operator(Cd[k]));
      COp.push_back(new Operator(totOps[k]));
      const Operator* arrC[2] = {CdOp[k-1],COp[k-1]};
      mpo.setOp(k,new JoinedOperator(2,arrC),true);  
    }
  }
}

 
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


void constructDoubleCommutatorHeisembergAncilla(MPO& mpo,int L,double J1x,double J1y,double J2x, double J2y, double J2z, double Bond){
  //cout<<"constructDoubleCommutatorIsing(g="<<g<<",h="<<h<<")"<<endl;
  // Basic operators
  int nrOps=7;
  int dd=d*d;
  mwArray Z=mwArray(Indices(nrOps,dd,dd));
   // Basic pieces
  mwArray sig0=identityMatrix(d);//identity
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),datax);//sigmax
  complex_t datay[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),datay);//sigmay
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz

  complex_t valJ1x=sqrt(J1x)*ONE_c;double sgnJ1x=J1x<0?-1.:1.;
  complex_t valJ1y=sqrt(J1y)*ONE_c;double sgnJ1y=J1y<0?-1.:1.;
  complex_t valJ2x=sqrt(J2x)*ONE_c;double sgnJ2x=J2x<0?-1.:1.;
  complex_t valJ2y=sqrt(J2y)*ONE_c;double sgnJ2y=J2y<0?-1.:1.;
  complex_t valJ2z=sqrt(J2z)*ONE_c;double sgnJ2z=J2z<0?-1.:1.;
  complex_t valBond=sqrt(Bond)*ONE_c;double sgnBond=Bond<0?-1.:1.;


  // *** First of all, we need the identity on both
  mwArray term0=identityMatrix(dd); //id*id
  // *** Now the ones appearing in single-body terms
  // (1) sigX (x) Id
  mwArray term1;
  constructOperatorProduct(term1,sigX,sig0); //sigX*id
  // (2)  Id (x) sigX -> just a permutation of the previous one
  mwArray term2;
  constructOperatorProduct(term2,sig0,permute(sigX,Indices(2,1))); //id*sigX^T
  // (3) sigY (x) Id
  mwArray term3;
  constructOperatorProduct(term3,sigY,sig0); //sigY*id
  // (4)  Id (x) sigY^T -> just a permutation of the previous one
  mwArray term4;
  constructOperatorProduct(term4,sig0,permute(sigY,Indices(2,1))); //id*sigY^T
  // (5) sigZ (x) Id
  mwArray term5;
  constructOperatorProduct(term5,sigZ,sig0); //sigZ*id
  // (6)  Id (x) sigZ^T -> just a permutation of the previous one
  mwArray term6;
  constructOperatorProduct(term6,sig0,permute(sigZ,Indices(2,1))); //id*sigZ^T

  // Fill in the operators in the Z array, to be contracted with C in
  // the construction of the MPO. The order is as listed above.
  for(int d1=0;d1<dd;d1++)
    for(int d2=0;d2<dd;d2++){
      Z.setElement(term0.getElement(Indices(d1,d2)),Indices(0,d1,d2)); //id*id
      Z.setElement(term1.getElement(Indices(d1,d2)),Indices(1,d1,d2)); //sigX*id
      Z.setElement(term2.getElement(Indices(d1,d2)),Indices(2,d1,d2)); //id*sigX^T
      Z.setElement(term3.getElement(Indices(d1,d2)),Indices(3,d1,d2)); //sigY*id
      Z.setElement(term4.getElement(Indices(d1,d2)),Indices(4,d1,d2)); //id*sigY^T
      Z.setElement(term5.getElement(Indices(d1,d2)),Indices(5,d1,d2)); //sigZ*id
      Z.setElement(term6.getElement(Indices(d1,d2)),Indices(6,d1,d2)); //id*sigZ^T
    }
  Z.reshape(Indices(nrOps,dd*dd));
  //cout<<"Finished Z:"<<Z.getDimensions()<<endl;
  // now the coefficients
  int D=14; // bond dimension of the MPO
  
  mwArray Ce(Indices(D,D,nrOps)); // the middle site even
  mwArray Co(Indices(D,D,nrOps)); // the middle site odd
  mwArray CeL(Indices(D,1,nrOps)); // the last even site
  mwArray CoL(Indices(D,1,nrOps)); // the last site
  mwArray Ce0(Indices(1,D,nrOps)); // the first site

  /*  CREATE EVEN PART */
  Ce.setElement(ONE_c,Indices(0,0,0)); 
  Ce.setElement(sgnJ1x*valJ1x,Indices(0,1,1));
  Ce.setElement(valJ1x,Indices(1,D-1,1)); 
  Ce.setElement(sgnJ1x*valJ1x,Indices(0,2,2));
  Ce.setElement(-ONE_c*valJ1x,Indices(2,D-1,2));
  
  Ce.setElement(sgnJ1y*valJ1y,Indices(0,3,3)); 
  Ce.setElement(valJ1y,Indices(3,D-1,3)); 
  Ce.setElement(sgnJ1y*valJ1y,Indices(0,4,4)); 
  Ce.setElement(-ONE_c*valJ1y,Indices(4,D-1,4)); 

  Ce.setElement(ONE_c,Indices(0,5,5));
  Ce.setElement(ONE_c,Indices(5,D-1,5)); 
  Ce.setElement(ONE_c,Indices(0,6,6));
  Ce.setElement(-ONE_c,Indices(6,D-1,6));

  Ce.setElement(ONE_c,Indices(7,7,0));
  Ce.setElement(ONE_c,Indices(8,8,0));
  Ce.setElement(ONE_c,Indices(9,9,0));
  Ce.setElement(ONE_c,Indices(10,10,0));
  Ce.setElement(ONE_c,Indices(11,11,0));
  Ce.setElement(ONE_c,Indices(12,12,0));
  Ce.setElement(ONE_c,Indices(D-1,D-1,0));

  /*  CREATE SECOND LAST PART */
  CeL=Ce;
  CeL.setElement(ZERO_c,Indices(0,0,0));
  CeL.setElement(ZERO_c,Indices(0,1,1));
  CeL.setElement(ZERO_c,Indices(0,2,2));
  CeL.setElement(ZERO_c,Indices(0,3,3));
  CeL.setElement(ZERO_c,Indices(0,4,4));

  /*  CREATE ODD PART */ 
  Co.setElement(Bond*ONE_c,Indices(5,D-1,5));  
  Co.setElement(-ONE_c*Bond,Indices(6,D-1,6)); 

  Co.setElement(sgnJ2x*valJ2x,Indices(0,7,1)); 
  Co.setElement(valJ2x,Indices(7,D-1,1)); 
  Co.setElement(sgnJ2x*valJ2x,Indices(0,8,2)); 
  Co.setElement(-ONE_c*valJ2x,Indices(8,D-1,2)); 
  
  Co.setElement(sgnJ2y*valJ2y,Indices(0,9,3)); 
  Co.setElement(valJ2y,Indices(9,D-1,3)); 
  Co.setElement(sgnJ2y*valJ2y,Indices(0,10,4)); 
  Co.setElement(-ONE_c*valJ2y,Indices(10,D-1,4));

  Co.setElement(sgnJ2z*valJ2z,Indices(0,11,5));
  Co.setElement(valJ2z,Indices(11,D-1,5)); 
  Co.setElement(sgnJ2z*valJ2z,Indices(0,12,6));
  Co.setElement(-ONE_c*valJ2z,Indices(12,D-1,6));

  Co.setElement(ONE_c,Indices(0,0,0)); 
  Co.setElement(ONE_c,Indices(1,1,0)); 
  Co.setElement(ONE_c,Indices(2,2,0)); 
  Co.setElement(ONE_c,Indices(3,3,0)); 
  Co.setElement(ONE_c,Indices(4,4,0)); 
  Co.setElement(ONE_c,Indices(5,5,0)); 
  Co.setElement(ONE_c,Indices(6,6,0)); 
  Co.setElement(ONE_c,Indices(D-1,D-1,0));

  /*  CREATE FIRST PART */
  Ce0.setElement(ONE_c,Indices(0,0,0)); 
  Ce0.setElement(sgnJ1x*valJ1x,Indices(0,1,1));
  Ce0.setElement(sgnJ1x*valJ1x,Indices(0,2,2));
  Ce0.setElement(sgnJ1y*valJ1y,Indices(0,3,3)); 
  Ce0.setElement(sgnJ1y*valJ1y,Indices(0,4,4)); 
  Ce0.setElement(ONE_c,Indices(0,5,5));
  Ce0.setElement(ONE_c,Indices(0,6,6));

  /*  CREATE LAST PART */
  CoL.setElement(valJ1x,Indices(7,0,1)); 
  CoL.setElement(-ONE_c*valJ1x,Indices(8,0,2)); 
  CoL.setElement(valJ1y,Indices(9,0,3)); 
  CoL.setElement(-ONE_c*valJ1y,Indices(10,0,4));
  CoL.setElement(ONE_c,Indices(11,0,5)); 
  CoL.setElement(-ONE_c,Indices(12,0,6));
  CoL.setElement(ONE_c,Indices(D-1,0,0));

  // Now set the operators in place
  // Now I reshape, contract with Z and set the indices in proper order

  Ce.reshape(Indices(D*D,nrOps));Ce.multiplyRight(Z);
  Ce.reshape(Indices(D,D,dd,dd));Ce.permute(Indices(3,1,4,2));

  Co.reshape(Indices(D*D,nrOps));Co.multiplyRight(Z);
  Co.reshape(Indices(D,D,dd,dd));Co.permute(Indices(3,1,4,2));
  
  CeL.reshape(Indices(D*D,nrOps));CeL.multiplyRight(Z);
  CeL.reshape(Indices(D,D,dd,dd));CeL.permute(Indices(3,1,4,2));
  
  Ce0.reshape(Indices(D,nrOps));Ce0.multiplyRight(Z);
  mwArray id2=1/sqrt(2.)*identityMatrix(d);id2.reshape(Indices(dd,1));
  Ce0.multiplyRight(id2);
  Ce0.reshape(Indices(1,D,dd,1));Ce0.permute(Indices(3,1,4,2));
  
  CoL.reshape(Indices(D*1,nrOps));CoL.multiplyRight(Z);
  CoL.reshape(Indices(D*1*dd,dd));
  CoL.multiplyRight(id2);
  CoL.reshape(Indices(D,1,dd,1));CoL.permute(Indices(3,1,4,2)); 

  // Now I build the MPO with required ops.
  //1) on the edges, I contract C1 and CN with their adjoint
  Ce0.reshape(Indices(dd,D));
  CoL.reshape(Indices(dd,D)); // dd D
  mwArray tmp(Ce0);tmp.Hconjugate(); 
  Ce0.multiplyLeft(tmp); // Dl,Dr
  tmp=CoL;tmp.Hconjugate(); 
  CoL.multiplyLeft(tmp); //Dlu,Dld
  // 2) These, I contract with the middle ones 
  tmp=Co;tmp.permute(Indices(2,1,3,4));tmp.reshape(Indices(D,dd*dd*D));
  mwArray edgeL=Ce0*tmp;edgeL.reshape(Indices(D*dd,1,dd,D)); // Ddd 1 dd D
  tmp=CeL;tmp.reshape(Indices(dd*D*dd,D));
  CoL.permute(Indices(2,1)); // put first the index for 2nd MPO
  mwArray edgeR=tmp*CoL;edgeR.reshape(Indices(dd,D,dd,D));
  edgeR.permute(Indices(1,4,2,3));edgeR.reshape(Indices(dd*D,D,dd,1)); // ddD D dd 1
  // Now for the adjoint I also need modifications
  mwArray Cde(Ce);Cde.permute(Indices(3,2,1,4),true); // adjoint of normal C
  mwArray Cdo(Co);Cdo.permute(Indices(3,2,1,4),true); // adjoint of normal C
  mwArray CdeL(CeL);CdeL.permute(Indices(3,2,1,4),true); // adjoint of normal C
 

  mwArray edgeLd(Cdo);edgeLd.reshape(Indices(dd,1,D*dd,D));
  mwArray edgeRd(CdeL);edgeRd.reshape(Indices(dd,D,dd*D,1));
  
  Operator CdOpE(Cde),COpE(Ce);
  Operator CdOpO(Cdo),COpO(Co);
  
  Operator edgeLOp(edgeL),edgeROp(edgeR);
  Operator edgeLdOp(edgeLd),edgeRdOp(edgeRd);
  const Operator* arrL[2]={&edgeLdOp,&edgeLOp};
  const Operator* arrR[2]={&edgeRdOp,&edgeROp};
  const Operator* arrCe[2]={&CdOpE,&COpE};
  const Operator* arrCo[2]={&CdOpO,&COpO};
  //cout<<"About to place operators in MPO"<<endl;
  // And now I set the operators one by one
  mpo.initLength(2*L);
  mpo.setOp(0,new JoinedOperator(2,arrL),true);
  mpo.setOp(2*L-1,new JoinedOperator(2,arrR),true);
  if(L>2){
    mpo.setOp(1,new JoinedOperator(2,arrCo),true);
    mpo.setOp(2,new JoinedOperator(2,arrCe),true);
    for(int k=3 ; k<2*L-1 ; k++){
      mpo.setOp(k++,&mpo.getOp(1),false);
      mpo.setOp(k,&mpo.getOp(2),false);
      
    }
  }
}
