#include <math.h>
#include <iomanip>
#include <fstream>
#include <deque>
#include <vector>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "JoinedOperator.h"
#include "Properties.h"
#include "mwArray.h"
#include "HeisenbergHamiltonian.h"
#include "IsingHamiltonian.h"
#include "Indices.h"
#include "Operator.h"

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
void constructDoubleCommutatorHeis(MPO& mpo,int L,double h, mwArray instances, int rowInstances, vector<double> &randFields);
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

const string jobfilename(int M,int D,double hp,const string mpsdir);
const string mpsfilename(int M,int D,const string mpsdir);

void adapt(MPO& Hmpo,MPO& HmpoD);

#define MAXITER 100
#define MAXPERT 10

//#define ISINGCASE 1

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
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  // Store the command line, for this is what should be relaunched if not all levels are done
  stringstream commandstr;
  for(int c=0;c<argc;c++) commandstr<<argv[c]<<" ";


  int M=props.getIntProperty("M");

  int rowInstances=props.getIntProperty("rowInstances");
  double tol = props.getDoubleProperty("tol");
  double hP=props.getDoubleProperty("hP");
  string outfile=props.getProperty("outfile");
  string inf=props.getProperty("instances");
  const char* infname = inf.c_str();
  
  mwArray instances;
  ifstream in(infname);
  if(!in.is_open()){
    cout<<"Error: couldn't open file "<<infname<<" for output"<<endl;
    exit(1);
  }
  instances.load(in);
  in.close();
  
  Contractor& contractor=Contractor::theContractor(); //inizialize contractor
  
  //  cout<<instances.getDimensions()<<endl;

  ofstream* out;
  out=new ofstream(outfile.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output"<<endl;
    exit(1);
  }
  out->close();delete out;
  double cosN,cosD, lambda2;

  //  while(M<=80){
  
    vector<double> Jx(M-1,1.);
    vector<double> Jy(M-1,1.);
    vector<double> Jz(M-1,1.);
    vector<double> randFields(M+2,0.);
    vector<double> fullRandFields(M,0.);
    vector<double> randsforMPS(M,0.);
    double N=M-1;
    MPO H2(M);
    
    constructDoubleCommutatorHeis(H2,M,10.,instances,rowInstances,randFields);

    // if(M<10){
    //   cout<<"instances=[";
    //   for(int j=0; j<M;j++) cout<<80.*real(instances.getElement(Indices(rowInstances,j)))<<",";
    //   cout<<"]"<<endl;
    //   cout<<"h=[";
    //   for(int j=0; j<M;j++) cout<<randFields[j]<<",";
    //   cout<<"]"<<endl;
    //   cout<<"Jx=[";
    //   for(int j=0; j<M-1;j++) cout<<Jx[j]<<",";
    //   cout<<"]"<<endl;
    // }


    //randsforMPS=randFields;
    for(int k=0;k<M-1;k++){
      //Jx[k]=cos(M_PI*(k-N/2)/N)/2;
      //Jx[k]=1.;
      //randFields[k] = randFields[k]*2;
      //randsforMPS[k] = randFields[k]*cos(M_PI*(k-N/2)/N);
      fullRandFields[k] = randFields[k+1];
      randsforMPS[k] = cos(M_PI*(k-N/2)/N);
    }
    //randFields[M-1] = randFields[M-1]*2;
    //randsforMPS[M-1]=randFields[M-1]*cos(M_PI*(M-1-N/2)/N);
    randsforMPS[M-1] = cos(M_PI*(M-1-N/2)/N);
    HeisenbergHamiltonian hamH(M,0.,0.,0.,randsforMPS,2);
    //HeisenbergHamiltonian hamH(M,0.,0.,0.,fullRandFields,2);
 
   
    MPS myMPS(M,5,4);
    MPSfromMPO(hamH.getHMPO(), myMPS, true);
    
    double lambda=real(contractor.contract(myMPS,H2,myMPS));
    double norm=real(contractor.contract(myMPS,myMPS));

    //if(M<10){
    //   cout<<"randsforMPS=[";
    //   for(int j=0; j<M;j++) cout<<randsforMPS[j]<<",";
    //   cout<<"]"<<endl;
    //   cout<<"Jx=[";
    //   for(int j=0; j<M-1;j++) cout<<Jx[j]<<",";
    //   cout<<"]"<<endl;
    // }
    cout<<"h=[";
    for(int j=0; j<M+1;j++) cout<<randFields[j]<<",";
    cout<<randFields[M+1]<<"]"<<endl;

    // cout<<"lambda: "<<lambda<<" norm: "<<norm<<endl;
    // //hamH.getHMPO().exportForMatlab("hamH.m");
    myMPS.exportForMatlab("myMPS.m");
    H2.exportForMatlab("H2.m");
    mwArray Mham;
    expandOper(H2, Mham);
    if(Mham.isMatrix()) cout<<"iuppiii"<<endl;
    vector<complex_t> eigvector; mwArray U;double tool;
    //wrapper::eig(Mham, eigvector, U, b = false);
    cout<<"I am computing the eigenvalues...the trace(H*H) = "<<(Mham*Mham).trace()<<endl;
    int eigenval = wrapper::eigs(Mham, 100, "LM", eigvector, U, false, tool=0);
    cout<<endl<<eigvector<<endl<<endl;


    

    // cosN=0;cosD=0;
    // for(int j=0; j<M-1; j++){
    //   cosN += cos(M_PI*(j-N/2)/N)*cos(M_PI*(j-N/2)/N)*(12+2*(randFields[j]*randFields[j]+randFields[j+1]*randFields[j+1]));
    //   cosN += -cos(M_PI*(j-N/2)/N)*cos(M_PI*(j+1-N/2)/N)*(6+2*randFields[j+1]*randFields[j+1]);
    //   cosN += -cos(M_PI*(j-N/2)/N)*cos(M_PI*(j-1-N/2)/N)*(6+2*randFields[j]*randFields[j]);
    //   cosD += cos(M_PI*(j-M/2)/M)*cos(M_PI*(j-M/2)/M)*(3+(randFields[j]*randFields[j]+randFields[j+1]*randFields[j+1])/4);
    // }
    
    // lambda2 =cosN/(cosD);
    // cout<<endl<<"lambda/norm: "<<lambda/norm<<" cosN/cosD: "<<lambda2<<endl;
    
    cosN=0;cosD=0;
    for(int j=0; j<M-1; j++){
      cosN += 2*(cos(M_PI*(j-N/2)/N)-cos(M_PI*(j-1-N/2)/N))*(cos(M_PI*(j-N/2)/N)-cos(M_PI*(j-1-N/2)/N));
      cosD += cos(M_PI*(j-N/2)/N)*cos(M_PI*(j-N/2)/N);
    }
    
    lambda2=4*cosN/cosD;
    cout<<endl<<"lambda/norm: "<<lambda/norm<<" cosN/cosD: "<<lambda2<<endl;
    

    out=new ofstream(outfile.data(),ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfile<<"for output"<<endl;
      exit(1);
    }
    *out<<setprecision(15);
    *out<<M<<"\t"<<lambda/norm<<"\t"<<lambda2<<endl;
    out->close();
    delete out;
    
    M=M+1;
    Jx.clear();Jy.clear();Jz.clear();randFields.clear();
    myMPS.clear();

    //}
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


  //MPS Omps(M,D,d*d);
  //Omps.setRandomState();
  // while(D<120){

  //   string outputMPS=mpsfilename(M,D,mpsdir); 
  //   const string tmpfilename=outputMPS+".tmp";
    
  //   substractTrace(Omps);
  //   cout<<"substracted trace from random Omps"<<endl;
    
  //   Omps.gaugeCond('R',1);Omps.gaugeCond('L',1); 
    
  //   //Omps.exportForMatlab("OmpsInit.m");
  //   cout<<"Imposed the gauge condition on Omps"<<endl;
  //   cout<<"Now going to contract with H2="<<H2<<endl;
  //   double oldLambda=real(contractor.contract(Omps,H2,Omps));
  //   cout<<"constructed (or recovered) initial MPS with lambda="<<oldLambda<<", trace "<<computeTrace(Omps)<<" and <O|O>="<<contractor.contract(Omps,Omps)<<endl;
    
  //   if(Omps.getBond()<D){
  //     forceHermitian(Omps);    
  //     Omps.increaseBondDimensionWithNoise(D,noise);
  //     substractTrace(Omps);
  //   }
  //   // Now apply the minimization with a penalty if trace is not zero.
  //   MPS idPen(M,1,d*d);
  //   for(int k=0;k<M;k++)
  //     idPen.setA(k,reshape(identityMatrix(d),Indices(d*d,1,1)));
  //   idPen.gaugeCond('R',0);
  //   vector<MPS*> orthog;orthog.push_back(&idPen);
  //   vector<double> penalty(1,-1.*pen);
  //   double lambda=real(contractor.contract(Omps,H2,Omps));
  //   contractor.setEigenSolver(primme);
  //   cout<<setprecision(15);
  //   contractor.findGroundStateWithProjectorPenalty(H2,D,penalty,orthog,&lambda,Omps,0.,1,tmpfilename);
  //   cout<<"Found  eigenval="<<lambda<<endl;
  //   cout<<"Trace is "<<computeTrace(Omps)
  // 	<<" and <O|O>="<<contractor.contract(Omps,Omps)<<endl;
    
  //   Omps.exportMPS(outputMPS.data());

  //   if(file_exists(tmpfilename)){
  //     remove(tmpfilename.data());
  //   }

  //   out=new ofstream(outfile.data(),ios::app);
  //   if(!out->is_open()){
  //     cout<<"Error: impossible to open file "<<outfile<<
  // 	" for output"<<endl;
  //     exit(1);
  //   }
  //   *out<<setprecision(15);
  //   *out<<D<<"\t"<<lambda<<"\t"<<M<<"\t"<<computeTrace(Omps)<<endl;
  //   out->close();
  //   delete out;

  //   D += 20;
  //   if(D<20){
  //     MPS Omps(M,D,d*d);
  //     string inputMPS=mpsfilename(M,D-2,mpsdir);
  //     Omps.importMPS(inputMPS.data());
  //   }
  // }
