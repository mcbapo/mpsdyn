#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "JoinedOperator.h"
#include "Properties.h"
#include "BasicMultiplier.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0


using namespace shrt;
using namespace std;

complex_t computeTrace(const MPS& Omps);

const string mpsfilename(int M,int D,const string mpsdir);
const string jobfilename(int M,int D,double gP,double hP,const string jobsdir);
void preparejobfile(Properties& props);
void constructMPOrange(MPO& mpoL,int range,int len);
int d=2;


//Construct Basic MPO, namely HxId-IdxH^T
//void constructCommutatorIsing(MPO& mpo,int L,double g,double h);
// Construct the real MPo that I need to use in the GS search, namely the square of the commutator, with the edge sites traced out
void constructAdjointMPO(MPO& HmpoD,const MPO& Hmpo);

// Construct the MPO for the double commutator, directly.
void constructDoubleCommutatorIsing(MPO& mpo,int L,double g,double h);



/**     
	From the results (MPS) saved by slowMPOham, construct smaller
	range MPOs by projecting the edges on the identity and compute
	the Frobenius norm of the corresponding commutator with H. In
	this case I use all possible projections onto M-2.
  */

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
  string outfile=props.getProperty("outputfile");
  double gP=props.getDoubleProperty("gP");
  double hP=props.getDoubleProperty("hP");

  string fileMPS=mpsfilename(M,D,mpsdir);

  Contractor& contractor=Contractor::theContractor();
  //  contractor.setConvTol(tol);
  cout<<"Initialized Contractor"<<endl;
  //  srandom(rSeed);

  ofstream* out;
  out=new ofstream(outfile.data());
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output"<<endl;
    exit(1);
  }
  // 1) Read the already computed state
  MPS Omps(M,D,d*d);
  if(file_exists(fileMPS)){
    Omps.importMPS(fileMPS.data());
    //cout<<"Imported MPS from file "<<fileMPS<<endl;
    //cout<<"Omps="<<Omps<<endl;
  }
  else{  // no tmp file
    cout<<"There is no MPS file "<<fileMPS<<" for specified parameters! NOTHING DONE !"<<endl;
    // Pass the required info
    // props.addProperty("command",argv[0]);
    // props.addProperty("propFile",argv[1]);
    // preparejobfile(jobsdir,props);
    exit(15);
  }

  // Starting on the MPS I just read, and going down in M, until 2
  // compute the commutator
  MPO H2(M);
  constructDoubleCommutatorIsing(H2,M,gP,hP);
  double lambda=real(contractor.contract(Omps,H2,Omps));
  double norm=real(contractor.contract(Omps,Omps)); // has to be one, but just in case
  // write out the old result
  *out<<M<<"\t"<<lambda<<endl;
  cout<<"Result for M "<<M<<" -> "<<lambda<<endl;
  while(M>2){
    // trace out things
    MPS auxL(Omps),auxR(Omps),auxC(Omps);
    auxL.traceOutSites(Indices(0,1),true);
    auxC.traceOutSites(Indices(0,M-1),true);
    auxR.traceOutSites(Indices(M-2,M-1),true);
    // save in an array of pointers just for convenience
    vector<const MPS*> auxPtr;
    auxPtr.push_back(&auxL);
    auxPtr.push_back(&auxC);
    auxPtr.push_back(&auxR);
    M=M-2;
    // compute the commutator
    constructDoubleCommutatorIsing(H2,M,gP,hP);
    // and now the matrix elements in the basis of these three vectors
    mwArray effH2(Indices(3,3));
    mwArray effN(Indices(3,3));
    for(int i=0;i<3;i++){
      for(int j=i;j<3;j++){
	complex_t Hij=contractor.contract(*auxPtr[j],H2,*auxPtr[i]);
	complex_t Nij=contractor.contract(*auxPtr[j],*auxPtr[i]);
	if(j==i) effH2.setElement(real(Hij),0,Indices(i,j));
	else{
	  effH2.setElement(Hij,Indices(i,j));
	  effH2.setElement(conjugate(Hij),Indices(j,i));
	}
	if(j==i) effN.setElement(real(Nij),0,Indices(i,j));
	else{
	  effN.setElement(Nij,Indices(i,j));
	  effN.setElement(conjugate(Nij),Indices(j,i));
	}
      }
    }
    // Now the solution should be the min gen. eigenvalue, but since
    // it is so small, I will do it by hand:
    // H v = lambda N v
    // N=Un S Un^dagger
    // [S^(-1/2) Un^dagger H U S^(-1/2)] [S^(1/2) U^dagger v] = lambda [S^(1/2) U^dagger v]
    //    cout<<"Solving the GEVP A="<<effH2<<", B="<<effN<<endl;
    // 1) diagonalize the norm:
    vector<complex_t> eigN,eigV; mwArray Un,U;
    wrapper::eig(effN,eigN,Un,true);
    //cout<<"Diagonalization of effN gives S="<<eigN<<", Un="<<Un<<endl;
    // 2) construct the effective matrix to be diagonalized next
    mwArray sqrtS(Indices(eigN.size(),eigN.size())); // sqrt(diag(eigN)); //
    for(int i=0;i<eigN.size();i++){
      sqrtS.setElement(sqrt(real(eigN[i])),0,Indices(i,i));
    }
    int nr=0;double tol=1E-8;
    mwArray sqrtSinv=invertDiag(sqrtS,nr,tol);
    //cout<<"sqrt(S)^-1 gives "<<sqrtSinv<<endl;
    effH2.multiplyLeft(sqrtSinv*Hconjugate(Un));
    effH2.multiplyRight(Un*sqrtSinv);
    effH2=.5*(effH2+Hconjugate(effH2));
    //cout<<"And transformed effH2="<<effH2<<endl;
    // 3) diagonalize [S^(-1/2) Un^dagger H U S^(-1/2)] 
    wrapper::eig(effH2,eigV,U,true);
    // 4) pick up the minimum eigenvalue and this is [S^(1/2) U^dagger v]
    //cout<<"Diagonalized [S^(-1/2) Un^dagger H U S^(-1/2)] with eigenvalues "<<eigV<<endl;
    double lambda=real(eigV[0]);
    // 5) recover the corresponding v
    U.multiplyLeft(Un*sqrtSinv); // first column contains the desired coefficients
    // write out the result
    *out<<M<<"\t"<<lambda<<endl;
    cout<<"Result for M "<<M<<" -> "<<lambda<<endl;
    if(M>2){
      // Find the best MPO approximation to the same sum!
      Omps=auxC;
      vector<complex_t> alpha;
      for(int i=0;i<3;i++)
	alpha.push_back(U.getElement(Indices(i,0)));
      contractor.optimizeSum(auxPtr,alpha,Omps,D);
      Omps.gaugeCond('r',true); // normalize
    }
  }


  out->close();delete out;
  // if(D<260){
  //   // Check if result was converged
  //   int newD=D+20; // in general
  //   // other cases
  //   if(D<20){
  //     if(D==1) newD=2;
  //     if(D==2) newD=4;
  //     if(D==4) newD=10;
  //     if(D==10) newD=D+10;
  //   }
  //   // Once everything is done, launch the next job, for next D
  //   // First, manipulate the properties, to include the new M 
  //   props.addProperty("command",argv[0]);
  //   props.addProperty("propFile",argv[1]);
  //   props.addProperty("D",newD);
  //   preparejobfile(jobsdir,props);
  // }


}


void preparejobfile(Properties& props){
  int M=props.getIntProperty("M");
  int D=props.getIntProperty("D");
  //  const char* outfnameG=argv[++cntr];
  // const char* outfnameL=argv[++cntr];
  string mpsdir=props.getProperty("mpsdir");
  string jobsdir=props.getProperty("jobsdir"); 
  string outfile=props.getProperty("weightsfile");
  string command=props.getProperty("command");
  string propFile=props.getProperty("propFile");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-5;
  double gP=props.getDoubleProperty("gP");
  double hP=props.getDoubleProperty("hP");
  stringstream commandstr;
  commandstr<<command;
  commandstr<<" "<<propFile;
    commandstr<<" -gP="<<gP<<" -hP="<<hP;
    commandstr<<" -M="<<M<<" -weightsfile="<<outfile;
    commandstr<<" -jobsdir="<<jobsdir;
    commandstr<<" -tol="<<tol<<" -D="<<D;
    commandstr<<" -mpsdir="<<mpsdir;
    const string jobfile=jobfilename(M,D,gP,hP,jobsdir);
    ofstream* ojob=new ofstream(jobfile.data());
    if(!ojob->is_open()){
      cout<<"WARNING!!!! Could not open the JOB file "<<jobfile<<". Launch by hand with the following command:"<<endl;
      cout<<commandstr.str()<<endl;
      cout<<endl;
      exit(15);
    }
    else{
      *ojob<<"#!/bin/bash"<<endl;
      *ojob<<"msub_modified_tqo097 ";
      //if(D+20>60) *ojob<<" -P 4";
      *ojob<<" -N slowHloc."<<hP<<"."<<M<<"."<<D;
      *ojob<<" -raw ";
      *ojob<<commandstr.str()<<endl;
      ojob->close();
    }
    delete ojob;
}


const string mpsfilename(int M,int D,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS"<<"_M"<<M<<"_D"<<D<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

const string jobfilename(int M,int D,double gP,double hP,const string jobsdir){
  stringstream s;
  s<<jobsdir<<"/job_Hloc_M"<<M<<"_D"<<D<<"_g"<<gP<<"_h"<<hP;  
  s<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}


void constructMPOrange(MPO& mpoL,int range,int len){
  //cout<<"In MPO(L="<<len<<"), range="<<range<<endl;
  if(len!=mpoL.getLength())
    mpoL.initLength(len);
  if(range>len){
    cout<<"ERROR! range cannot be larger than MPO!"<<endl;
    exit(1);
  }
  // Basic operators, construct them once, then reuse
  static bool init=false;
  static mwArray Z(Indices(3,d*d,d*d));
  if(!init){
    mwArray O0=.5*identityMatrix(d*d);
    O0.reshape(Indices(d,d,d,d));
    O0.permute(Indices(2,4,1,3));
    O0.reshape(Indices(d*d,d*d));
    mwArray O1=identityMatrix(d*d);
    // construct the Z mwArray
    for(int i1=0;i1<d*d;i1++)
      for(int i2=0;i2<d*d;i2++){
	Z.setElement(O0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
	Z.setElement(O1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
	Z.setElement(O1.getElement(Indices(i1,i2))-O0.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      }
    Z.reshape(Indices(3,d*d*d*d));
    init=true;
    //cout<<"Initialized the Z operator for the MPO"<<endl;
  }

  // DEBUGGING ONLY: I build the MPS alone to check coeffs
  MPS aux(len,range+1,3);

  // MPO terms. Special case, k=1=L
  int D=range+1;int Dl(D),Dr(D);
  for(int p=0;p<len;p++){
    if(p==0) Dl=1; else Dl=D;
    if(p==len-1) Dr=1; else Dr=D;
    mwArray C(Indices(Dl,Dr,3));
    if(len==1){ // special case: just one term
      C.setElement(ONE_c,Indices(Dl,Dr,2));
    }
    else{
      if(p<len-range)
	C.setElement(ONE_c,Indices(0,0,0));
      if(p>range-1)
	C.setElement(ONE_c,Indices(Dl-1,Dr-1,0));
      if(p<=len-range)
	if(p<len-1)
	  C.setElement(ONE_c,Indices(0,1,2)); // left operator
	else
	  C.setElement(ONE_c,Indices(0,0,2)); // left operator is also last
      for(int ik=1;ik<range-1;ik++){
	if(p>=ik&&(len-1-p)>=range-ik-1) // place for the rest of the term on both sides
	  C.setElement(ONE_c,Indices(ik,ik+1,1)); // in between
      }
      if(p>=range-1)
	C.setElement(ONE_c,Indices(range-1,Dr-1,2)); // in between
    }
    aux.replaceSite(p,permute(C,Indices(3,1,2)));
    //cout<<"Coefficient of site "<<p<<", C="<<C<<endl;
    C.reshape(Indices(Dl*Dr,3));
    C.multiplyRight(Z);
    C.reshape(Indices(Dl,Dr,d*d,d*d));
    C.permute(Indices(3,1,4,2));
    mpoL.setOp(p,new Operator(C),true); // no reusing operators!
  }
  // stringstream nameMPS;
  // nameMPS<<"mpsForRange"<<range<<"_L"<<len<<".m";
  // aux.exportForMatlab(nameMPS.str().data());
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


// Copied from Liouvillian, should be somewhere else
void constructOperatorProduct(mwArray& result,const mwArray& opA,
			      const mwArray& opB){
#ifdef CHECKDIMS
  if(opA.getDimensions()!=opB.getDimensions()){
    cout<<"Error in constructOperatorProduct for"
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
  // cout<<"constructDoubleCommutatorIsing(g="<<g<<",h="<<h<<")"<<endl;
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
