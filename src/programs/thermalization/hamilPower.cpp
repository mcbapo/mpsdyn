#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "JoinedOperator.h"
#include "Properties.h"
#include "IsingHamiltonian.h"

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

// Construct the MPO for the energy modulation with the specified coefficients
void constructModulatedEnergy(MPO& mpo,int M,const vector<double>& Js,const vector<double>& gs,const vector<double>& hs);


void computeAllTermsHn(int n,const MPO& H12,const MPS& Omps,int M,ofstream* out);
bool incrementRecursive(vector<int>& indices,int M,int& ptr);
void prepareMPOHn(const vector<int>& indices,const MPO& H12,MPO& Hn,int M);

void solveGEV(const mwArray& effH2,const mwArray& effN,double* lambda);

/**     
	Compute the commmutator of concrete Hamiltonian functions
	(e.g. modulated energy density and its powers)
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
  //  int D=props.getIntProperty("D");
  //  const char* outfnameG=argv[++cntr];
  // const char* outfnameL=argv[++cntr];
  string mpsdir=props.getProperty("mpsdir");
  string outfile=props.getProperty("outputfile");
  //double J=props.getDoubleProperty("J");
  double gP=props.getDoubleProperty("gP");
  double hP=props.getDoubleProperty("hP");
  int maxN=props.getIntProperty("maxN");
  int minN=props.getIntProperty("minN");
  if(minN<0) minN=1; // default
  //  string fileMPS=mpsfilename(M,D,mpsdir);

  //  contractor.setConvTol(tol);
  // cout<<"Initialized Contractor"<<endl;
  //  srandom(rSeed);
  int D=3;
  MPS Omps(M,D,d*d);
  MPO modH(M);
  vector<double> Js,gs,hs;
  double J=1;
  // Modulation cos(pi*(i-M/2)/M);
  for(int k=0;k<M;k++){
    double valk=cos(M_PIl*(k-M/2+1)/(M-1));
    double valkl=cos(M_PIl*(k-1-M/2+1)/(M-1));
    int last=k==M-1?1:0;
    int first=k==0?1:0;
    double Jmod=J*valk;
    double gmod=(1-last)*gP*.5*valk+(1-first)*gP*.5*valkl;
    double hmod=(1-last)*hP*.5*valk+(1-first)*hP*.5*valkl;
    if(!last) Js.push_back(Jmod);
    gs.push_back(gmod);
    hs.push_back(hmod);
  }
  cout<<"Modulated coefficients: Js "<<Js<<", gs "<<gs<<", hs "<<hs<<endl;
  constructModulatedEnergy(modH,M,Js,gs,hs);
  // Now transform in MPS
  MPSfromMPO(modH,Omps,true);
  Omps.gaugeCond('R',true);
  // And also prepare the MPO x Id for larger powers
  MPO modHdoub(M);
  extendMPO(modH,modHdoub,d);

  // The basic commutator
  MPO Hcom(M);
  constructDoubleCommutatorIsing(Hcom,M,gP,hP);

  // Abasis for the ops
  vector<MPS*> opBasis;
  MPS* aux=new MPS(Omps);
  opBasis.push_back(aux); // n=1
  int nTerms=opBasis.size();
  mwArray effN(Indices(nTerms,nTerms));effN.fillWithZero();
  mwArray effH2(Indices(nTerms,nTerms));effH2.fillWithZero();
  
  Contractor& contractor=Contractor::theContractor();
  MPS mpsId(M,1,d*d);
  {  
    MPO idMPO(M);
    mwArray id2=identityMatrix(2);
    for(int k=0;k<M;k++){
      idMPO.setOp(k,new Operator(reshape(id2,Indices(2,1,2,1))),true);
    }
    MPSfromMPO(idMPO,mpsId,d);
    mpsId.gaugeCond('R',true); // normalize properly
  }
  for(int n=minN;n<=maxN;n++){ // Repeat for n-th power of the Hamiltonian
    // Compute the commutator of Hmod^n
    complex_t commN=contractor.contract(Omps,Hcom,Omps);
    complex_t norm=contractor.contract(Omps,Omps); //should be 1
    complex_t overlapId=contractor.contract(Omps,mpsId); //if not 0
    if(n>minN){
      nTerms++;
      effN.resize(Indices(nTerms,nTerms));
      effH2.resize(Indices(nTerms,nTerms));
    }
    effH2.setElement(commN,Indices(nTerms-1,nTerms-1)); //diagonal terms
    effN.setElement(norm-conjugate(overlapId)*overlapId,Indices(nTerms-1,nTerms-1));
    for(int i=0;i<nTerms-1;i++){
      MPS* vi=opBasis[i];
      complex_t Nin=contractor.contract(Omps,*vi); //<vi|vn>
      complex_t Hin=contractor.contract(Omps,Hcom,*vi); //<vi|H2|vn>
      complex_t viId=contractor.contract(mpsId,*vi); //<vi|Id>
      effN.setElement(Nin-conjugate(viId)*overlapId,Indices(i,nTerms-1));//<vi|vn>-<vi|Id><Id|vn>
      effN.setElement(conjugate(Nin-conjugate(viId)*overlapId),Indices(nTerms-1,i));
      effH2.setElement(Hin,Indices(i,nTerms-1));
      effH2.setElement(conjugate(Hin),Indices(nTerms-1,i));
    }
    //cout<<"==n="<<n<<"====="<<endl;
    //putForMatlab(cout,effH2,"effH2");
    //putForMatlab(cout,effN,"effN");
    double lambda=10;
    solveGEV(effH2,effN,&lambda);
    ofstream* out;
    if(n==minN){
      out=new ofstream(outfile.data());
      *out<<"%n\t M\t lambda(Hn)\t min of H1...Hn\n";
    }
    else
      out=new ofstream(outfile.data(),ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfile<<
	" for output"<<endl;
      exit(1);
    }
    *out<<setprecision(10);
    *out<<n<<"\t"<<M<<"\t"<<real(commN)/real(norm)<<"\t"<<lambda<<endl;
    //cout<<"For n="<<n<<" and M="<<M<<", lambda(Hn)="<<real(commN)/real(norm)<<" and best LCfrom H1-Hn->"<<lambda<<endl;
    out->close();
    delete out;
    // now, if not at the end, calculate exactly the next power
    if(n<maxN){
      MPS* aux=new MPS(Omps);
      contractor.optimize(modHdoub,*aux,Omps,D*3);
      Omps.gaugeCond('R',true);
      *aux=Omps;      
      opBasis.push_back(aux);
    }
  }
}

void solveGEV(const mwArray& effH2_,const mwArray& effN_,double* lambda){
  // Copy the matrices formanipulation
  mwArray effH2(effH2_);mwArray effN(effN_);
    // 1) diagonalize the norm:
    vector<complex_t> eigN,eigV; mwArray Un,U;
    mwArray effH2b(effH2); // the original
    effN=.5*(effN+Hconjugate(effN));
    wrapper::eig(effN,eigN,Un,true);
    //cout<<"Diagonalization of effN gives S="<<eigN<<endl; //", Un="<<Un<<endl;
    // 2) construct the effective matrix to be diagonalized next I need
    // sqrt(S) and its inverse, but I will truncate to non-vanishing
    // eigenvalues only
    // how many?
    int nTerms=effH2.getDimension(0);
    int nr=nTerms;
    vector<complex_t>::iterator it=eigN.begin();
    while((real(*it)<1E-6)&&(++it!=eigN.end())) nr--;
    // cout<<"Nr of non-vanishingeigenvalues then="<<nr<<" (out of "<<nTerms<<")"<<endl;
    //mwArray sqrtS(Indices(nr,nr));sqrtS.fillWithZero();
    mwArray projNoNull(Indices(nTerms,nr)); projNoNull.fillWithZero();
    mwArray sqrtSinv(Indices(nr,nr));sqrtSinv.fillWithZero();
    for(int i=0;i<nr;i++){ //eigN.size()-nr;i<eigN.size();i++){
      double valI=real(eigN[nTerms-nr+i]);
      // sqrtS.setElement(sqrt(valI)*ONE_c,Indices(i,i));
      sqrtSinv.setElement((1./sqrt(valI))*ONE_c,Indices(i,i));
      projNoNull.setElement(ONE_c,Indices(nTerms-nr+i,i));
      //    cout<<"Setting element "<<Indices(i,i)<<" of sqrtS to "<<sqrt(real(eigN[i]))<<endl;
    }
    //    cout<<"sqrt(S)^-1 gives "<<sqrtSinv<<endl;
    effH2.multiplyLeft(sqrtSinv*Hconjugate(projNoNull)*Hconjugate(Un));
    effH2.multiplyRight(Un*projNoNull*sqrtSinv);
    effH2=.5*(effH2+Hconjugate(effH2));
    //cout<<"And transformed effH2="<<effH2<<endl;//.getDimensions()<<endl;
    // 3) diagonalize [S^(-1/2) Un^dagger H U S^(-1/2)] 
    wrapper::eig(effH2,eigV,U,false);
    // 4) pick up the minimum eigenvalue and this is [S^(1/2) U^dagger v]
    //cout<<"Diagonalized [S^(-1/2) Un^dagger H U S^(-1/2)] with eigenvalues "<<eigV<<endl;
    *lambda=real(eigV[0]);
    // 5) recover the corresponding v
    //  U.multiplyLeft(Un*sqrtSinv); // first column contains the desired coefficients
    // write out the result
    //*out<<M<<"\t"<<lambda<<endl;
}

void constructModulatedEnergy(MPO& mpo,int M,const vector<double>& Js,const vector<double>& gs,const vector<double>& hs){
  // Out of laziness, I just copy from IsingHamiltonian
  // basic spin operators appearing
  mwArray sig0=identityMatrix(d);
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sig3(Indices(d,d),dataz);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig1(Indices(d,d),datax);
  mwArray Z(Indices(3,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      //Z.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Z.setElement(sig3.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      //Z.setElement(sig3.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      Z.setElement(sig1.getElement(Indices(i1,i2)),Indices(2,i1,i2));
    }
  Z.reshape(Indices(3,d*d));
  int D=3;
  // But now all terms are different, because coefficients depend on
  // site
  mpo.initLength(M);
  for(int k=0;k<M;k++){
    //cout<<"Creating new Operator for site "<<k<<" of IsingH MPO"<<endl;
    int Dl=D; int Dr=D;
    if(k==0) Dl=1;
    if(k==M-1) Dr=1;
    int redefK=k+1-M/2; // from -M/2+1 until M/2
    mwArray C(Indices(Dl,Dr,3));C.fillWithZero();
    complex_t effJl=(Js[k]<0)?sqrt(abs(Js[k]))*I_c:sqrt(Js[k])*ONE_c;
    complex_t effJr=k==0?ZERO_c:((Js[k-1]<0)?sqrt(abs(Js[k-1]))*I_c:sqrt(Js[k-1])*ONE_c);
    complex_t effG=gs[k]*ONE_c;
    complex_t effH=hs[k]*ONE_c;
    // Identity
    if(k!=M-1)
      C.setElement(ONE_c,Indices(0,0,0)); // just pass 1 (before anything)
    if(k!=0)
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // just pass 1 (after)
    // sigmax terms (in all of them)
    C.setElement(effG,Indices(0,Dr-1,2));
    // sigmaz terms (in all of them)
    C.setElement(effH,Indices(0,Dr-1,1));
    // left sigz term (all but last)
    if(k!=M-1)
      C.setElement(effJl,Indices(0,1,1));
    // right sigz term (all but first)
    if(k!=0)
      C.setElement(effJr,Indices(1,Dr-1,1));
      
    // Now contract and give the operator
    // Now contract MPS and operators and set proper index ordering
    C.reshape(Indices(Dl*Dr,3));
    //cout<<"C for site "<<k<<"="<<C<<endl;
    mwArray res=C*Z;
    //cout<<"Computed res "<<res<<endl;
    res.reshape(Indices(Dl,Dr,d,d));
    res.permute(Indices(3,1,4,2));
    //cout<<"  with size: "<<res.getDimensions()<<endl;
    // Create and store an Operator
    mpo.setOp(k,new Operator(res),true);
  }
}

void computeAllTermsHn(int n,const MPO& H12,const MPS& Omps,int M,ofstream* out){
    MPO Hn(M);
    vector<int> indices(n,0); // all zeros
    bool done=false;
    Contractor& contractor=Contractor::theContractor();
    //double sumLambda=0;
    while(!done){
      prepareMPOHn(indices,H12,Hn,M);
      //cout<<"computeAllTermsHn indices"<<indices<<", Hn="<<Hn<<endl;
      // transform in MPS to compute overlap
      MPS aux(M,1,1);
      MPSfromMPO(Hn,aux,true);
      aux.gaugeCond('R',true);
      double lambda=abs(contractor.contract(Omps,aux));
      // cout<<"Hn normalized as "<<contractor.contract(aux,aux)<<endl;
      // cout<<"Operator normalized as "<<contractor.contract(Omps,Omps)<<endl;
      // cout<<"Overlap "<<contractor.contract(Omps,aux)<<", lambda="<<lambda<<endl;
      //      sumLambda+=lambda*lambda;
      //      cout<<" **** So far, cumulative sum of lambda^2="<<sumLambda<<endl;
      //double normAux=real(contractor.contract(aux,aux));
      for(int k=0;k<n;k++) *out<<indices[k]<<"\t";
      *out<<lambda<<endl;
      int ptr=n-1;
      done=!incrementRecursive(indices,M-2,ptr);
      //      cout<<"incrementRecursive returns "<<indices<<" (done="<<done<<")"<<endl;
    }

}

// Return false in case it is at the end (all max - M-1)
// I increment with the constraint that i1<=i2<=i3<=...iN
bool incrementRecursive(vector<int>& indices,int M,int& ptr){
  // cout<<"incrementRecursive("<<indices<<") pos "<<ptr<<endl;
  if(indices[ptr]<M){
    indices[ptr]++;
    for(int k=ptr+1;k<indices.size();k++){
      indices[k]=indices[ptr];
    }
    return true;
  }
  else{ // this one cannot be incremented => go left
    if(ptr>0)
      return incrementRecursive(indices,M,--ptr);
    else
      return false;
  }
}

void prepareMPOHn(const vector<int>& indices,const MPO& H12,MPO& Hn,int M){
  //  cout<<"Preparing Hn for indices "<<indices<<endl;
  Hn.clear();Hn.initLength(M);
  // initialize with all Ids
  mwArray id2=identityMatrix(2);
  for(int k=0;k<M;k++){
    Hn.setOp(k,new Operator(reshape(id2,Indices(2,1,2,1))),true);
  }
  vector<bool> occupied(M,false);
  int n=indices.size(); // how many terms
  for(int in=0;in<n;in++){
    for(int l=0;l<=1;l++){
      int pos=indices[in]+l;
      mwArray newOp=H12.getOp(l).getFullData();
      if(occupied[pos]){
	const Operator* toJoin[]={&Hn.getOp(pos),&H12.getOp(l)};
	JoinedOperator auxOp(2,toJoin);
	newOp=auxOp.getFullData();
      }
      Hn.setOp(pos,new Operator(newOp),true);
      occupied[pos]=true;
    }
  }
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

