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
void constructDoubleCommutatorHeisembergAncillaTrace(MPO& mpo,int L,double J1x,double J1y, double J1z,double J2x, double J2y, double J2z, double Bond);

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

  //Hamiltonian Parameters
  double hP=props.getDoubleProperty("hP");
  double J1x=props.getDoubleProperty("J1x");
  double J1y=props.getDoubleProperty("J1y");
  double J1z=props.getDoubleProperty("J1z");
  double J2x=props.getDoubleProperty("J2x");
  double J2y=props.getDoubleProperty("J2y");
  double J2z=props.getDoubleProperty("J2z");

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
  //double J1x=1.,J1y=1.,J1z=1.,J2x=0.1,J2y=0.1,J2z=0.1;
  //double Bond=.5;
  
  constructDoubleCommutatorHeisembergAncillaTrace(H2,M,J1x,J1y,J1z,J2x,J2y,J2z,hP);

  // H2.exportForMatlab("H2.m");
  // mwArray Mham;
  // expandOper(H2, Mham);
  // //if(Mham.isMatrix()) cout<<"iuppiii"<<endl;
  // vector<complex_t> eigvector; 
  // mwArray U; 
  // bool b;
  // double tool;
  // //wrapper::eig(Mham, eigvector, U, b = false);
  // cout<<"I am computing the eigenvalues..."<<endl;
  // int eigenval = wrapper::eigs(Mham, 30, "LM", eigvector, U, b = false, tool=0);
  // cout<<endl<<eigvector<<endl<<endl;

  int phys[2*M];
  for(int i=0;i<2*M;i++){
    if (i%2 == 0) {phys[i]=d*d;}
    else{phys[i]=1;}
  }
 
  //Define penalities
  MPS idPen(2*M,1,phys);
  for(int k=0;k<2*M;k++){
    if(k%2==0){
      idPen.setA(k,reshape(identityMatrix(d),Indices(d*d,1,1)));
    }else{
      idPen.setA(k,ONE_c);
    }
    
  }
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
  MPS Omps(2*M,D,phys);

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
    int dd = dims[0];
    int ds = sqrt(dd);
    int newDl=k==0?dims[1]:dims[1]+1;
    int newDr=k==L-1?dims[2]:dims[2]+1;
    A.resize(Indices(dd,newDl,newDr));
    A.reshape(Indices(ds,ds,newDl,newDr));
    //    cout<<"Site "<<k<<" going from dims="<<dims<<" to "<<A.getDimensions()<<endl;
    double sign=k==0?-1.*signTr:1.;
    for(int j=0;j<ds;j++)
      A.setElement(sign*val*ONE_c,Indices(j,j,newDl-1,newDr-1)); 
    A.reshape(Indices(dd,newDl,newDr));
    Omps.replaceSite(k,A,false);
  }
}

 
complex_t computeTrace(const MPS& Omps){
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  int len=Omps.getLength();
  int dim[len];
  for(int k=0; k<len ; k++) {
        mwArray A=Omps.getA(k).getA();
	Indices dims=A.getDimensions();
	dim[k] = dims[0];
  }
  MPS aux(len,1,dim);
  for(int k=0;k<len;k++){
    aux.setA(k,id2);
    k++;
  }
  Contractor& contractor=Contractor::theContractor();
  return contractor.contract(Omps,aux);
}
 
void forceHermitian(MPS& Omps){
  int len=Omps.getLength();
  for(int k=0;k<len;k++){
    mwArray A=Omps.getA(k).getA();
    Indices dims=A.getDimensions(); // dd,Dl,Dr
    int ds=sqrt(dims[0]);
    A.reshape(Indices(ds,ds,dims[1],dims[2]));
    mwArray B(A);
    B.permute(Indices(2,1,3,4),true);
    A=.5*(A+B);    
    A.reshape(Indices(dims[0],dims[1],dims[2]));
    Omps.setA(k,A);
  }
}


void constructDoubleCommutatorHeisembergAncillaTrace(MPO& mpo,int L,double J1x,double J1y, double J1z,double J2x, double J2y, double J2z, double Bond){
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

  //Create the 5 building blocks of the MPO
  //Si and Ti stay for Pauli matrix i

  /*  CREATE EVEN PART */

  // |1 Sx1 1Sx Sy1 1Sy Sz1 1Sz 0 0 0 0 0 0 0   |      
  // |                                     Sx1  |         
  // |                                     -1Sx |       
  // |                                     Sy1  |      
  // |                                     -1Sy |       
  // |                                     Sz1  |       
  // |                                     -1Sz |       
  // |                          1           0   |      
  // |                            1         0   |      
  // |                              1       0   |      
  // |                                1     0   |      
  // |                                  1   0   |      
  // |                                    1 0   |      
  // |                                      1   |      


  Ce.setElement(ONE_c,Indices(0,1,1)); 
  Ce.setElement(ONE_c*J1x,Indices(1,D-1,1)); 
  Ce.setElement(ONE_c,Indices(0,2,2));
  Ce.setElement(-ONE_c*J1x,Indices(2,D-1,2));
  
  Ce.setElement(ONE_c,Indices(0,3,3)); 
  Ce.setElement(ONE_c*J1y,Indices(3,D-1,3)); 
  Ce.setElement(ONE_c,Indices(0,4,4)); 
  Ce.setElement(-ONE_c*J1y,Indices(4,D-1,4)); 

  Ce.setElement(ONE_c,Indices(0,5,5));
  Ce.setElement(ONE_c*J1z,Indices(5,D-1,5)); 
  Ce.setElement(ONE_c,Indices(0,6,6));
  Ce.setElement(-ONE_c*J1z,Indices(6,D-1,6));

  Ce.setElement(ONE_c,Indices(0,0,0));  
  Ce.setElement(ONE_c,Indices(7,7,0));
  Ce.setElement(ONE_c,Indices(8,8,0));
  Ce.setElement(ONE_c,Indices(9,9,0));
  Ce.setElement(ONE_c,Indices(10,10,0));
  Ce.setElement(ONE_c,Indices(11,11,0));
  Ce.setElement(ONE_c,Indices(12,12,0));
  Ce.setElement(ONE_c,Indices(D-1,D-1,0));

  /*  CREATE SECOND LAST PART */

  // | 0   0   0   0   0  Sz1 1Sz 0 0 0 0 0 0 0   |      
  // |                                       Sx1  |         
  // |                                       -1Sx |       
  // |                                       Sy1  |      
  // |                                       -1Sy |       
  // |                                       Sz1  |       
  // |                                       -1Sz |       
  // |                            1           0   |      
  // |                              1         0   |      
  // |                                1       0   |      
  // |                                  1     0   |      
  // |                                    1   0   |      
  // |                                      1 0   |      
  // |                                        1   |      


  CeL=Ce;
  CeL.setElement(ZERO_c,Indices(0,0,0));
  CeL.setElement(ZERO_c,Indices(0,1,1));
  CeL.setElement(ZERO_c,Indices(0,2,2));
  CeL.setElement(ZERO_c,Indices(0,3,3));
  CeL.setElement(ZERO_c,Indices(0,4,4));

  /*  CREATE ODD PART */ 
 
  // | 1             Tx1 1Tx Ty1 1Ty Tz1 1Tz  0   |           
  // |   1                                    0   |      
  // |     1                                  0   |      
  // |       1                                0   |      
  // |         1                              0   |           
  // |           1                         Tz1*B  |      
  // |             1                      -Tz1*B  |      
  // |                                       Tx1  |         
  // |                                      -1Tx  |       
  // |                                       Ty1  |      
  // |                                      -1Ty  |       
  // |                                       Tz1  |       
  // |                                      -1Tz  |       
  // |                                        1   |      

  Co.setElement(Bond*ONE_c,Indices(5,D-1,5));  
  Co.setElement(-ONE_c*Bond,Indices(6,D-1,6)); 

  Co.setElement(ONE_c,Indices(0,7,1)); 
  Co.setElement(ONE_c*J2x,Indices(7,D-1,1)); 
  Co.setElement(ONE_c,Indices(0,8,2)); 
  Co.setElement(-ONE_c*J2x,Indices(8,D-1,2)); 
  
  Co.setElement(ONE_c,Indices(0,9,3)); 
  Co.setElement(ONE_c*J2y,Indices(9,D-1,3)); 
  Co.setElement(ONE_c,Indices(0,10,4)); 
  Co.setElement(-ONE_c*J2y,Indices(10,D-1,4));

  Co.setElement(ONE_c,Indices(0,11,5));
  Co.setElement(ONE_c*J2z,Indices(11,D-1,5)); 
  Co.setElement(ONE_c,Indices(0,12,6));
  Co.setElement(-ONE_c*J2z,Indices(12,D-1,6));

  Co.setElement(ONE_c,Indices(0,0,0)); 
  Co.setElement(ONE_c,Indices(1,1,0)); 
  Co.setElement(ONE_c,Indices(2,2,0)); 
  Co.setElement(ONE_c,Indices(3,3,0)); 
  Co.setElement(ONE_c,Indices(4,4,0)); 
  Co.setElement(ONE_c,Indices(5,5,0)); 
  Co.setElement(ONE_c,Indices(6,6,0)); 
  Co.setElement(ONE_c,Indices(D-1,D-1,0));

  /*  CREATE FIRST PART */

    // |1 Sx1 1Sx Sy1 1Sy Sz1 1Sz 0 0 0 0 0 0 0   |      

  Ce0.setElement(ONE_c,Indices(0,0,0)); 
  Ce0.setElement(ONE_c,Indices(0,1,1));
  Ce0.setElement(ONE_c,Indices(0,2,2));
  Ce0.setElement(ONE_c,Indices(0,3,3)); 
  Ce0.setElement(ONE_c,Indices(0,4,4)); 
  Ce0.setElement(ONE_c,Indices(0,5,5));
  Ce0.setElement(ONE_c,Indices(0,6,6));

  /*  CREATE LAST PART */

  // |   0    |           
  // |   0    |      
  // |   0    |      
  // |   0    |      
  // |   0    |           
  // |  Tz1*B |      
  // | -Tz1*B |      
  // |  Tx1   |         
  // | -1Tx   |       
  // |  Ty1   |      
  // | -1Ty   |       
  // |  Tz1   |       
  // | -1Tz   |       
  // |   1    |      


  CoL.setElement(Bond*ONE_c,Indices(5,0,5));  
  CoL.setElement(-ONE_c*Bond,Indices(6,0,6)); 
  CoL.setElement(ONE_c,Indices(D-1,0,0));
  CoL.setElement(ONE_c*J2x,Indices(7,0,1)); 
  CoL.setElement(-ONE_c*J2x,Indices(8,0,2)); 
  CoL.setElement(ONE_c*J2y,Indices(9,0,3)); 
  CoL.setElement(-ONE_c*J2y,Indices(10,0,4));
  CoL.setElement(ONE_c*J2z,Indices(11,0,5)); 
  CoL.setElement(-ONE_c*J2z,Indices(12,0,6));


  // Now set the operators in place
  // Now I reshape, contract with Z and set the indices in proper order
  mwArray id2=1/sqrt(2.)*identityMatrix(d);id2.reshape(Indices(dd,1));

  Ce.reshape(Indices(D*D,nrOps));Ce.multiplyRight(Z);
  Ce.reshape(Indices(D,D,dd,dd));Ce.permute(Indices(3,1,4,2));
  
  //I create the right operators @ the borders
  //and then I contract them with the identities 


  CeL.reshape(Indices(D*D,nrOps));CeL.multiplyRight(Z);
  CeL.reshape(Indices(D*D*dd,dd)); 
  CeL.multiplyRight(id2);
  CeL.reshape(Indices(D,D,dd,1));CeL.permute(Indices(3,1,4,2));//dd D 1 D

  Ce0.reshape(Indices(D,nrOps));Ce0.multiplyRight(Z);
  Ce0.reshape(Indices(1*D*dd,dd));
  Ce0.multiplyRight(id2);
  Ce0.reshape(Indices(1,D,dd,1));Ce0.permute(Indices(3,1,4,2));//dd 1 1 D

  CoL.reshape(Indices(D*1,nrOps));CoL.multiplyRight(Z);
  CoL.reshape(Indices(D*1*dd,dd));
  CoL.multiplyRight(id2);
  CoL.reshape(Indices(D,1,dd*1));CoL.permute(Indices(1,3,2));  

  Co.reshape(Indices(D*D,nrOps));Co.multiplyRight(Z);
  Co.reshape(Indices(D*D*dd,dd));
  Co.multiplyRight(id2);
  Co.reshape(Indices(D,D,dd,1));Co.permute(Indices(3,1,4,2));//dd D 1 D

  int d5 = dd*dd;

  mwArray tmp(Co);
  tmp.permute(Indices(2,1,3,4));
  tmp.reshape(Indices(D,dd*1*D));

  
  Ce0.reshape(Indices(dd,D));
  Ce0.multiplyRight(tmp);
  Ce0.reshape(Indices(d5,1,1,D)); 

  tmp=CoL; tmp.reshape(Indices(D,dd*1*1));
  CeL.reshape(Indices(dd*D*1,D));
  CeL.multiplyRight(tmp);
  CeL.reshape(Indices(dd,D,1,dd*1,1)); CeL.permute(Indices(1,4,2,3,5));
  CeL.reshape(Indices(d5,D,1,1));

  // NOW EXACLY THE SAME PROCEDURE AS ISING
  //C1 = Ce0 , CN = CeL and C = Ce 
  Ce0.reshape(Indices(d5,D));
  CeL.reshape(Indices(d5,D)); // dd D
  tmp=Ce0;tmp.Hconjugate(); 
  Ce0.multiplyLeft(tmp); // D,D
  tmp=CeL;tmp.Hconjugate(); 
  CeL.multiplyLeft(tmp); //D,Dl

  // 2) These, I contract with the middle ones 
  tmp=Ce;tmp.permute(Indices(2,1,3,4));tmp.reshape(Indices(D,dd*dd*D));
  mwArray edgeL=Ce0*tmp;edgeL.reshape(Indices(D*dd,1,dd,D)); // Ddd 1 dd D
  tmp=Co;tmp.reshape(Indices(dd*D*1,D));
  CeL.permute(Indices(2,1)); // put first the index for 2nd MPO
  mwArray edgeR=tmp*CeL;edgeR.reshape(Indices(dd,D,1,D));
  edgeR.permute(Indices(1,4,2,3));edgeR.reshape(Indices(dd*D,D,1,1)); // ddD D dd 1

  // Now for the adjoint I also need modifications
 
  // mwArray CdeL(CeL);CdeL.permute(Indices(3,2,1,4),true); // adjoint of normal 
  // mwArray edgeRdeL(CdeL);edgeRdeL.reshape(Indices(dd,D,dd*D,1));
 
  mwArray Cde(Ce);Cde.permute(Indices(3,2,1,4),true); // adjoint of normal C
  mwArray edgeLde(Cde);edgeLde.reshape(Indices(dd,1,D*dd,D));

  mwArray Cdo(Co);Cdo.permute(Indices(3,2,1,4),true); // adjoint of normal 
  mwArray edgeRdo(Cdo);edgeRdo.reshape(Indices(1,D,dd*D,1));
  
  Operator CdOpe(Cde),COpe(Ce);
  Operator CdOpo(Cdo),COpo(Co);
  Operator edgeLOp(edgeL),edgeROp(edgeR);
  Operator edgeLdOp(edgeLde),edgeRdOp(edgeRdo);

  const Operator* arrCe[2]={&CdOpe,&COpe};
  const Operator* arrCo[2]={&CdOpo,&COpo};
 
  const Operator* arrL[2]={&edgeLdOp,&edgeLOp};
  const Operator* arrR[2]={&edgeRdOp,&edgeROp};

  //cout<<"About to place operators in MPO"<<endl;
  // And now I set the operators one by one
  mpo.initLength(2*L);
  mpo.setOp(0,new JoinedOperator(2,arrL),true);
  mpo.setOp(2*L-1,new JoinedOperator(2,arrR),true);
  if(L>2){
    mpo.setOp(1,new JoinedOperator(2,arrCo),true);
    mpo.setOp(2,new JoinedOperator(2,arrCe),true);
    int k = 3;
    while(k < 2*L-1 ){
      mpo.setOp(k,&mpo.getOp(1),false); k += 1;
      mpo.setOp(k,&mpo.getOp(2),false); k += 1;
      
    }
  }
}

