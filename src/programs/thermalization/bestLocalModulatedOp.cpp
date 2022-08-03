#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "JoinedOperator.h"
#include "Properties.h"
#include "IsingHamiltonian.h"
#include <strstream>

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0


using namespace shrt;
using namespace std;

complex_t computeTrace(const MPS& Omps);

void constructMPOrange(MPO& mpoL,int range,int len);
int d=2;


//Construct Basic MPO, namely HxId-IdxH^T
//void constructCommutatorIsing(MPO& mpo,int L,double g,double h);
// Construct the real MPo that I need to use in the GS search, namely the square of the commutator, with the edge sites traced out
void constructAdjointMPO(MPO& HmpoD,const MPO& Hmpo);

// Construct the MPO for the double commutator, directly.
void constructDoubleCommutatorIsing(MPO& mpo,int L,double g,double h);

// Return false in case it is at the end (all max - M-1)
bool incrementRecursive(vector<int>& indices,int M,int& ptr);
// I increment WITHOUT the constraint that i1<=i2<=i3<=...iN
bool incrementRecursiveAll(vector<int>& indices,int M,int& ptr);

// Construct h12, with flexibility to set the single body terms balanced or not
void constructIsingHamiltonianTwoSites(MPO& h12,double J,double g1,double h1,double g2,double h2);

void addModulatedMPOtoBasis(vector<int>& indices,int M,vector<MPS*>& opBasis);
double modulation(int M,int R,int pos,const vector<int>& indices);
double modulatedFunction(int M,int R,int pos);


int modulMode;

/**     
	Construct explicitly the best linear combinations of local operators,
	of the form 
	\f[
	\sum_n f(n) C_{\alpha_1 \ldots \alpha_k} \sigma_{\alpha_1} \ldots \sigma_{\alpha_k}
	\f]
	for a given range \f$k\f$, where the function \f$f(n)\f$ could
	in principle be a cosine modulation.
	The program should optimize the coefficients \f$C_{\alpha_1
	\ldots \alpha_k}\f$
	for fixed \f$k\f$.
	Here we construct an MPO for each modulated sum 
	\f[
	\sum_n f(n) \sigma_{\alpha_1} \ldots \sigma_{\alpha_k}
	\f]
	and then proceed as before.

	Othermodulations can be tried by providing argument 
	-modulation=mode
	mode: 0=cosine (default); 1=Gaussian; 2=const; 3=Gaussian
	with different width
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

  int M=props.getIntProperty("M");
  //  int D=props.getIntProperty("D");
  //  const char* outfnameG=argv[++cntr];
  // const char* outfnameL=argv[++cntr];
  //  string mpsdir=props.getProperty("mpsdir");
  string outfile=props.getProperty("outputfile");
  double gP=props.getDoubleProperty("gP");
  double hP=props.getDoubleProperty("hP");
  int maxK=props.getIntProperty("maxK");
  modulMode=props.getIntProperty("modulation");
  if(modulMode<0) modulMode=0;

  //  string outfileIm=outfile+"_im";

  //  contractor.setConvTol(tol);
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  //  srandom(rSeed);

  MPO H2(M);
  constructDoubleCommutatorIsing(H2,M,gP,hP);
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
  // Now construct all the terms: for each set (a1,...ak), the MPO
  // with the proper modulation
  vector<MPS*> opBasis;

  // loop over 4^maxK terms
  bool done=false;
  vector<int> indices(maxK,0); // all of them Identity shouldn't
			       // contribute, so I startwith a 1 in
			       // the last position
  indices[maxK-1]=1;
  while(!done){
    //cout<<"Constructing MPS for set "<<indices<<endl;
    //construct the modulated MPO and add it to the basis
    addModulatedMPOtoBasis(indices,M,opBasis);
    // increment the indices
    int ptr=maxK-1; //try to increment the last one
    done=!incrementRecursiveAll(indices,3,ptr);
  }
  // Now compute the effective H and N matrices and solve
  int nTerms=opBasis.size();
  mwArray effN(Indices(nTerms,nTerms));effN.fillWithZero();
  mwArray effH2(Indices(nTerms,nTerms));effH2.fillWithZero();
  for(int i=0;i<nTerms;i++){
    MPS* vi=opBasis[i];
    if(i==0) vi->gaugeCond('R',true); 
    complex_t viId=contractor.contract(mpsId,*vi); //<vi|Id>
    effN.setElement(ONE_c-conjugate(viId)*viId,Indices(i,i));
    effH2.setElement(contractor.contract(*vi,H2,*vi),Indices(i,i));
    for(int j=i+1;j<nTerms;j++){
      MPS* vj=opBasis[j];
      if(i==0) vj->gaugeCond('R',true); 
      complex_t vjId=contractor.contract(mpsId,*vj); //<vj|Id>	
      complex_t Nij=contractor.contract(*vj,*vi); //<vi|vj>
      complex_t Hij=contractor.contract(*vj,H2,*vi); //<vi|H2|vj>
      //if(abs(Hij)>.1)
      //cout<<"Computed <v"<<i<<"|H2|v"<<j<<">="<<Hij<<endl;
      effN.setElement(Nij-conjugate(viId)*vjId,Indices(i,j));
      effN.setElement(conjugate(Nij)-conjugate(vjId)*viId,Indices(j,i));
      effH2.setElement(Hij,Indices(i,j));
      effH2.setElement(conjugate(Hij),Indices(j,i));
    }
  }
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
  wrapper::eig(effH2,eigV,U,true);
  // 4) pick up the minimum eigenvalue and this is [S^(1/2) U^dagger v]
  //cout<<"Diagonalized [S^(-1/2) Un^dagger H U S^(-1/2)] with eigenvalues "<<eigV<<endl;
  double lambda=real(eigV[0]);
  // 5) recover the corresponding v
  U.multiplyLeft(Un*projNoNull*sqrtSinv); // first column contains the desired coefficients
  // write out the result
  //*out<<M<<"\t"<<lambda<<endl;
  cout<<"Result for M="<<M<<", maxK="<<maxK<<" -> "<<lambda<<endl; //" all eigV"<<eigV<<endl;
  
  ofstream* out;
  out=new ofstream(outfile.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);
  *out<<M<<"\t"<<maxK<<"\t"<<lambda<<"\t";
  // I also want the coefficients!!
  // WARNING: BY HAND I only write up to 4^4
  int cutCoeff=pow(4,maxK);
  for(int k=0;k<cutCoeff;k++)
    if(k<U.getDimension(0)) *out<<real(U.getElement(Indices(k,0)))<<"\t"<<imag(U.getElement(Indices(k,0)))<<"\t";    
    else *out<<0<<"\t"<<0<<"\t";    
  *out<<endl;
  out->close();
  // ofstream* outI;
  // outI=new ofstream(outfileIm.data(),ios::app);
  // if(!outI->is_open()){
  //   cout<<"Error: impossible to open file "<<outfileIm<<
  //     " for output"<<endl;
  //   exit(1);
  // }
  // *outI<<setprecision(10);
  // *outI<<M<<"\t"<<maxK<<"\t"<<lambda<<"\t";
  // // I also want the coefficients!!
  // // WARNING: BY HAND I only write up to 4^n
  // for(int k=0;k<cutCoeff;k++)
  //   if(k<U.getDimension(0)) *outI<<imag(U.getElement(Indices(k,0)))<<"\t";    
  //   else *outI<<0<<"\t";    
  // *outI<<endl;
  // outI->close();
  // delete outI;
}

double modulatedFunctionCos(int M,int R,int pos){
    return cos(M_PIl*(pos+1.-(M-R+1)*.5)/(M-R));
}

double modulatedFunctionGauss(int M,int R,int pos){
  double n=pos+1-(M-R+1)*.5;
  double sigma=M*.25;
  return exp(-n*n/(2*sigma*sigma));
}

double modulatedFunctionGaussMultiSigma(int M,int R,int pos){
  double n=pos+1-(M-R+1)*.5;
  double sigma=(M-R+1)*.25;
  return exp(-n*n/(2*sigma*sigma));
}

double modulatedFunctionConst(int M,int R,int pos){
  return 1.;
}

double modulatedFunctionExp(int M,int R,int pos){
  double val=cos(M_PIl*(pos+1.-(M-R+1)*.5)/(M-R));
  double dist=(pos+1.-(M-R+1)*.5);
  return val*val*exp(-dist*dist/(M-R+1));
}

double modulatedFunctionDoub(int M,int R,int pos){
  return .5+.5*cos(2*M_PIl*(pos+1.-(M-R+1)*.5)/(M-R));
}

double (*funMod[])(int M,int R,int pos)={modulatedFunctionCos,modulatedFunctionGauss,
					 modulatedFunctionConst,
					 modulatedFunctionGaussMultiSigma};

double modulatedFunction(int M,int R,int pos){
  return (*funMod[modulMode])(M,R,pos);
}

double modulation(int M,int R,int pos,const vector<int>& indices){
  //cout<<"modulation M="<<M<<" R="<<R<<" pos="<<pos;
  int posE=pos;int Meff=M;int Reff=R;
  if(indices[0]*indices[R]==0){// only if some of the edges is Id, the
			       // effective range or pos may change
    int idL=0;int idR=R;
    while(idL<R&&indices[idL]==0){
      posE++;Reff--;idL++;// first non-null indices element
    }
    while(idR>idL&&indices[idR]==0){
      Reff--;idR--;
    }
    // cout<<"modulation of M="<<M<<", R="<<R<<", pos="<<pos<<", indices "
    // 	<<indices<<" Reff="<<Reff<<" posE="<<posE<<endl;
  }
  //  return modulatedFunctionDoub(Meff,Reff,posE);
  return modulatedFunction(Meff,Reff,posE);
}

double modulation(int M,int R,int pos,int ind0,int indR){
  int posE=ind0==0?pos+1:pos;
  int Re=R;
  Re=ind0==0?Re-1:Re;
  Re=indR==0?Re-1:Re;
  return modulatedFunction(M,Re,posE);
}

void addModulatedMPOtoBasis(vector<int>& indices,int M,vector<MPS*>& opBasis){
  // Have the common operators ready (don't know if it saves much)
  static bool init=false;
  static mwArray Z(Indices(4,d,d));
  if(!init){
    Z.fillWithZero();
    // The operators I need
    mwArray sig0=identityMatrix(d);
    complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    mwArray sigX(Indices(d,d),dataX);
    complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
    mwArray sigY(Indices(d,d),dataY);
    complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
    mwArray sigZ(Indices(d,d),dataZ);
    for(int i1=0;i1<d;i1++)
      for(int i2=0;i2<d;i2++){
	Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
	Z.setElement(sigX.getElement(Indices(i1,i2)),Indices(1,i1,i2));
	Z.setElement(sigY.getElement(Indices(i1,i2)),Indices(2,i1,i2));
	Z.setElement(sigZ.getElement(Indices(i1,i2)),Indices(3,i1,i2));
      }
    Z.reshape(Indices(4,d*d));
    init=true;
  }
  int range=indices.size();
  int newD=range+1;
  mwArray Abas(Indices(d*d,newD,newD));Abas.fillWithZero();
  // the one with the modulation, so I can multiply and sum
  mwArray Avar(Indices(d*d,newD,newD));Avar.fillWithZero();
  // set the proper elements according to vector indices
  for(int i1=0;i1<d*d;i1++){
    Avar.setElement(Z.getElement(Indices(indices[0],i1)),Indices(i1,0,1)); // the first one, which is variable
    Abas.setElement(Z.getElement(Indices(0,i1)),Indices(i1,0,0));
    Abas.setElement(Z.getElement(Indices(0,i1)),Indices(i1,newD-1,newD-1));
  }
  for(int p=1;p<indices.size();p++){
    for(int i1=0;i1<d*d;i1++){
      Abas.setElement(Z.getElement(Indices(indices[p],i1)),Indices(i1,p,p+1));
    }
  }

  MPS* mps=new MPS(M,newD,d*d);
  for(int k=0;k<M;k++){
    double cosFac=modulatedFunction(M,range-1,k);
    //double cosFac=modulation(M,range-1,k,indices);
    //cout<<"For site "<<k<<" coeff "<<cosFac<<endl;
    mwArray newA=Abas+cosFac*Avar;    
    if(k==0)
      mps->replaceSite(k,reshape(newA.subArray(Indices(-1,0,-1)),Indices(d*d,1,newD)),false);
    else
      if(k==M-1)
	mps->replaceSite(k,newA.subArray(Indices(-1,-1,newD-1)),false);
      else
	mps->replaceSite(k,newA,false);
  }
  opBasis.push_back(mps);
  //cout<<"Added mps for "<<indices<<", range="<<range<<": "<<*mps<<endl;
}

void addModulatedMPOtoBasis_ineff(vector<int>& indices,int M,vector<MPS*>& opBasis){
  // Construct the MPO
  MPO mpoSigmas(M);
  // It has four operators
  mwArray Z(Indices(4,d,d));Z.fillWithZero();
  // The operators I need
  mwArray sig0=identityMatrix(d);
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(sigX.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Z.setElement(sigY.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      Z.setElement(sigZ.getElement(Indices(i1,i2)),Indices(3,i1,i2));
    }
  Z.reshape(Indices(4,d*d));
  int maxK=indices.size();
  int D=maxK+1;

  for(int k=0;k<M;k++){
    //    cout<<"Setting coefficients for site "<<k<<endl;
    int Dl=D,Dr=D;
    if(k==0) Dl=1;
    if(k==M-1) Dr=1;
    mwArray C(Indices(Dl,Dr,4));C.fillWithZero();
    // Set elements
    if(k<M-1){ // the last one not, to avoid all Ids
      C.setElement(ONE_c,Indices(0,0,0)); // nothing yet
      //cout<<"Set 000"<<C<<endl;
    }
    if(k>0){
      C.setElement(ONE_c,Indices(maxK,Dr-1,0)); // after everything
      //cout<<"Set XX0"<<C<<endl;
    }
    if(k<M-1){ // the last one cannot start the set
      //      C.setElement(modulation(M,maxK-1,k,indices)*ONE_c,Indices(0,1,indices[0]));
      C.setElement(modulation(M,maxK-1,k,1,1)*ONE_c,Indices(0,1,indices[0])); // trick
									      // to
									      // ignore Ids
      //C.setElement(modulation(M,maxK-1,k,indices[0],indices[maxK-1])*ONE_c,Indices(0,1,indices[0]));
      //cout<<"Setting element 0,1,ind0="<<indices[0]<<" to "<<modulation(M,maxK-1,k)<<" ("<<C.getElement(Indices(0,1,indices[0]))<<")"<<endl;
    }
    for(int ik=1;ik<maxK;ik++){    
      //cout<<"Preparing left index "<<ik<<endl;
      if(k<M-1){
	//cout<<"site k("<<k<<")<M-1("<<M-1<<")"<<endl;
	if(k>0){
	  //cout<<"site k>0"<<endl;
	  C.setElement(ONE_c,Indices(ik,ik+1,indices[ik]));
	  //cout<<"Set "<<ik<<","<<ik+1<<"-ind("<<ik<<")="<<indices[ik]<<" "<<C<<endl;
	}
      }
      else{ // the last one only accepts the last one
	//cout<<"site k("<<k<<")!<M-1("<<M-1<<")"<<endl;
	if(ik==maxK-1){ // could remove this if
	  C.setElement(ONE_c,Indices(maxK-1,Dr-1,indices[maxK-1]));
	  //cout<<"Set "<<maxK-1<<"X-ind("<<maxK-1<<")="<<indices[maxK-1]<<" "<<C<<endl;
	}
      }
    }
    //cout<<"C="<<C<<endl;
    C.reshape(Indices(Dl*Dr,4));
    C.multiplyRight(Z);
    C.reshape(Indices(Dl,Dr,d,d));
    C.permute(Indices(3,1,4,2));
    mpoSigmas.setOp(k,new Operator(C),true);
  }
  // // To check
  // ostringstream nameDeb;
  // nameDeb<<"mpo_";for(int l=0;l<maxK;l++) nameDeb<<indices[l];nameDeb<<".m";
  // mpoSigmas.exportForMatlab(nameDeb.str().data());
  // Now mat to MPs and store
  MPS* aux=new MPS(M,1,1);
  MPSfromMPO(mpoSigmas,*aux,true);
  opBasis.push_back(aux);
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

// Return false in case it is at the end (all max - M-1)
// I increment WITHOUT the constraint that i1<=i2<=i3<=...iN
bool incrementRecursiveAll(vector<int>& indices,int M,int& ptr){
  // cout<<"incrementRecursive("<<indices<<") pos "<<ptr<<endl;
  if(indices[ptr]<M){
    indices[ptr]++;
    for(int k=ptr+1;k<indices.size();k++){
      indices[k]=0;
    }
    return true;
  }
  else{ // this one cannot be incremented => go left
    if(ptr>0)
      return incrementRecursiveAll(indices,M,--ptr);
    else
      return false;
  }
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

void constructIsingHamiltonianTwoSites(MPO& h12,double J,double g1,double h1,double g2,double h2){
  // basic spin operators appearing
  int d=2;
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
  int D; // bond dimension, depends on case
  mwArray C1,C2;
  if(h2==0&&g2==0){
    // then there is only ZZ+ (gX+hZ)x Id (bond dim 2)
    D=2;
    C1=mwArray(Indices(1,D,3));C1.fillWithZero();
    C2=mwArray(Indices(D,1,3));C2.fillWithZero();
    // single body sigmax terms 
    C1.setElement(g1*ONE_c,Indices(0,D-1,2));
    // sigmaz terms 
    C1.setElement(h1*ONE_c,Indices(0,D-1,1));
    // first half of ZZ term
    C1.setElement(J*ONE_c,Indices(0,0,1));
    // second half of ZZ
    C2.setElement(ONE_c,Indices(0,0,1));
    // or just identity if the operators were on site 1
    C2.setElement(ONE_c,Indices(D-1,0,0));
    //    cout<<"Case only g1,h1, C1="<<C1<<", C2="<<C2<<endl;
  }
  else{
    if(h1==0&&g1==0){ // all in site 2
      D=2;
      C1=mwArray(Indices(1,D,3));C1.fillWithZero();
      C2=mwArray(Indices(D,1,3));C2.fillWithZero();
      // single body sigmax terms 
      C2.setElement(g2*ONE_c,Indices(D-1,0,2));
      // sigmaz terms 
      C2.setElement(h2*ONE_c,Indices(D-1,0,1));
      // or just identity if the operators were on site 2
      C1.setElement(ONE_c,Indices(0,D-1,0));
      // first half of ZZ term
      C1.setElement(J*ONE_c,Indices(0,0,1));
      // second half of ZZ
      C2.setElement(ONE_c,Indices(0,0,1));
      //cout<<"Case only g2,h2, C1="<<C1<<", C2="<<C2<<endl;
    }
    else{ // all of them different from zero
      D=3;
      C1=mwArray(Indices(1,D,3));C1.fillWithZero();
      C2=mwArray(Indices(D,1,3));C2.fillWithZero();
      // single body sigmax terms 
      C1.setElement(g1*ONE_c,Indices(0,D-1,2));
      C2.setElement(g2*ONE_c,Indices(0,0,2));
      // sigmaz terms 
      C1.setElement(h1*ONE_c,Indices(0,D-1,1));
      C2.setElement(h2*ONE_c,Indices(0,0,1));
      // or identity on the first (before) or second (after single body)
      C1.setElement(ONE_c,Indices(0,0,0));
      C2.setElement(ONE_c,Indices(D-1,0,0));
      // first half of ZZ term
      C1.setElement(J*ONE_c,Indices(0,1,1));
      // second half of ZZ
      C2.setElement(ONE_c,Indices(1,0,1));
      //cout<<"Case all, C1="<<C1<<", C2="<<C2<<endl;
    }
  }
  // Finally, contract MPO terms and set them
  C1.reshape(Indices(1*D,3));
  C2.reshape(Indices(D*1,3));
  //cout<<"C for site "<<k<<"="<<C<<endl;
  mwArray res1=C1*Z;
  mwArray res2=C2*Z;
  //cout<<"Computed res "<<res<<endl;
  res1.reshape(Indices(1,D,d,d));
  res2.reshape(Indices(D,1,d,d));
  res1.permute(Indices(3,1,4,2));
  res2.permute(Indices(3,1,4,2));
  //cout<<"  with size: "<<res.getDimensions()<<endl;
  // Create and store an Operator
  h12.setOp(0,new Operator(res1),true);
  h12.setOp(1,new Operator(res2),true);
}
