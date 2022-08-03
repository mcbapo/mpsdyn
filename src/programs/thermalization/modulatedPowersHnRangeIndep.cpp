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

void computeAllTermsHn(int n,const MPO& H12,const MPS& Omps,int M,ofstream* out);
bool incrementRecursive(vector<int>& indices,int M,int& ptr);
void prepareMPOHn(const vector<int>& indices,const MPO& H12,MPO& Hn,int M);

// Construct h12, with flexibility to set the single body terms balanced or not
void constructIsingHamiltonianTwoSites(MPO& h12,double J,double g1,double h1,double g2,double h2);
// Construct all the MPS for valid operators from H^n
void prepareOperatorsHn(int n,const MPO& H12,int M,vector<MPS*>& terms);
void prepareOperatorsHnrange(int n,int r,const MPO& H12,int M,vector<MPS*>& terms,vector<double>& coeffs);
void prepareCoefficients(int n, int M, vector<double>& coeffs);
void prepareCoefficientsrange(int n, int r,int M, vector<double>& coeffs);
bool indicesWithinLimits(int exactOrder,int exactR,const vector<int>& indices);


double modulation(int M,int R,int pos);

void translateMPS(const MPS& mpsR,int M,vector<MPS*>& opBasis,vector<double>& coeffs);


/**     
	Compute best commutator from a LC of Hamiltonian terms, where
	products of few h12 at short range are modulated as cos.
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
  string outfile=props.getProperty("outputfile");
  double gP=props.getDoubleProperty("gP");
  double hP=props.getDoubleProperty("hP");
  int maxN=props.getIntProperty("maxN");
  //int minN=props.getIntProperty("minN");
  //if(minN<0) minN=1; // default
  //  int maxR=props.getIntProperty("maxR");
  //if(maxR<0) maxR=M; // default
  // I will just use terms which are shifted by one

  //  contractor.setConvTol(tol);
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  //  srandom(rSeed);

  // The basic piece is the two-body Hamiltonian term h12, as a MPO.
  IsingHamiltonian hamil2(2,2,1.,.5*gP,.5*hP);
  const MPO& H12=hamil2.getHMPO(0.);

  MPO H2(M);
  constructDoubleCommutatorIsing(H2,M,gP,hP);
  //  complex_t lambdaOr=contractor.contract(Omps,H2,Omps)/contractor.contract(Omps,Omps);

  // How many terms do I consider?
  // For a given n and r, I have 

  vector<vector<MPS*> > opBasis_; //(maxN*(maxN+1)/2,vector<MPS*>());
  vector<vector<double> > coeffs;//(maxN*(maxN+1)/2,vector<double>());
				 //// the modulation for each of them
  //I want to cunt where each n ends
  vector<int> nLim;
  for(int n=1;n<=maxN;n++){ // Repeat for n-th power of the
			       // Hamiltonian
    for(int r=0;r<n;r++){
      // First, I need to determine the independent range r terms with
      // n h12 factors (I use a fake length N)
      if(r+2<=M){
	vector<MPS*> _opBasis;
	vector<double> _coeffN; //ignored
	prepareOperatorsHnrange(n,r+1,H12,r+2,_opBasis,_coeffN);
	// For each of the, I include an independent set in opBasis,with
	// coefficients modulated by cosine (or whatever)
	for(int k=0;k<_opBasis.size();k++){
	  vector<MPS*> opBasis;
	  vector<double> coeffsK; 
	  translateMPS(*_opBasis[k],M,opBasis,coeffsK);
	  opBasis_.push_back(opBasis);
	  coeffs.push_back(coeffsK);
	  // cout<<"Prepared coefficients of a term n="<<n<<", r="<<r<<" (global nr "
	  //     <<opBasis_.size()<<"): "<<coeffsK<<endl;
	}
      }
    }
    // end of n
    nLim.push_back(opBasis_.size());
  }
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
  // Now compute the matrix of overlaps and the effective H for each
  // pair n,m
  int szNr=opBasis_.size();
  mwArray Heff_nm(Indices(szNr,szNr));
  mwArray Neff_nm(Indices(szNr,szNr));
  for(int posnr=0;posnr<szNr;posnr++){
    vector<MPS*>& opBasisN=opBasis_[posnr];
    vector<double>& coeffN=coeffs[posnr];
    mwArray vecN(coeffN,false);
    for(int posmr=posnr;posmr<szNr;posmr++){
      vector<MPS*>& opBasisM=opBasis_[posmr];
      vector<double>& coeffM=coeffs[posmr];
      mwArray vecM(coeffM,true);
      int nTerms=opBasisN.size();
      int mTerms=opBasisM.size();
      mwArray effN(Indices(nTerms,mTerms));effN.fillWithZero();
      mwArray effH2(Indices(nTerms,mTerms));effH2.fillWithZero();
      for(int i=0;i<nTerms;i++){
	MPS* vi=opBasisN[i];vi->gaugeCond('R',true); 
	complex_t viId=contractor.contract(mpsId,*vi); //<vi|Id>
	for(int j=0;j<mTerms;j++){
	  MPS* vj=opBasisM[j];
	  if(i==0) vj->gaugeCond('R',true); 
	  complex_t Nij=contractor.contract(*vj,*vi); //<vi|vj>
	  complex_t Hij=contractor.contract(*vj,H2,*vi); //<vi|H2|vj>
	  complex_t vjId=contractor.contract(mpsId,*vj); //<vj|Id>
	  effN.setElement(Nij-conjugate(viId)*vjId,Indices(i,j));//<vi|vj>-<vi|Id><Id|vj>
	  effH2.setElement(Hij,Indices(i,j));
	}
      }
      mwArray valH=vecN*effH2*vecM;
      mwArray valN=vecN*effN*vecM;
      Heff_nm.setElement(valH.getElement(0),Indices(posnr,posmr));
      Neff_nm.setElement(valN.getElement(0),Indices(posnr,posmr));
      if(posmr!=posnr){
	Heff_nm.setElement(conjugate(valH.getElement(0)),Indices(posmr,posnr));
	Neff_nm.setElement(conjugate(valN.getElement(0)),Indices(posmr,posnr));
      }
    }
  }

  //cout<<"Prepared effH="<<Heff_nm<<endl;
  //cout<<"And effN="<<Neff_nm<<endl;
  // And now run the optimization of the coefficients

  int totSz=Heff_nm.getDimension(0);
  for(int in=0;in<maxN;in++){ // for each value
    // 1) diagonalize the norm:
    vector<complex_t> eigN,eigV; mwArray Un,U;
    mwArray effN=.5*(Neff_nm+Hconjugate(Neff_nm));
    mwArray effH2=.5*(Heff_nm+Hconjugate(Heff_nm));
    int szIn=nLim[in];
    effN.resize(Indices(szIn,szIn));
    effH2.resize(Indices(szIn,szIn));
    wrapper::eig(effN,eigN,Un,true);
    //    cout<<"Diagonalization of effN gives S="<<eigN<<endl; //", Un="<<Un<<endl;
    // 2) construct the effective matrix to be diagonalized next I need
    // sqrt(S) and its inverse, but I will truncate to non-vanishing
    // eigenvalues only
    // how many?
    int nTerms=szIn;
    int nr=nTerms;
    vector<complex_t>::iterator it=eigN.begin();
    while((real(*it)<1E-6)&&(++it!=eigN.end())) nr--;
    //cout<<"Nr of non-vanishingeigenvalues then="<<nr<<" (out of "<<nTerms<<")"<<endl;
    //mwArray sqrtS(Indices(nr,nr));sqrtS.fillWithZero();
    mwArray projNoNull(Indices(nTerms,nr)); projNoNull.fillWithZero();
    mwArray sqrtSinv(Indices(nr,nr));sqrtSinv.fillWithZero();
    for(int i=0;i<nr;i++){ //eigN.size()-nr;i<eigN.size();i++){
      double valI=real(eigN[nTerms-nr+i]);
      // sqrtS.setElement(sqrt(valI)*ONE_c,Indices(i,i));
      sqrtSinv.setElement((1./sqrt(valI))*ONE_c,Indices(i,i));
      projNoNull.setElement(ONE_c,Indices(nTerms-nr+i,i));
      //cout<<"Setting element "<<Indices(i,i)<<" of sqrtS to "<<sqrt(real(eigN[nTerms-nr+i]))<<endl;
    }
    // cout<<"sqrt(S)^-1 gives "<<sqrtSinv<<endl;
    // cout<<"and Un is "<<Un<<endl;
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
    cout<<"Result for M="<<M<<", maxN="<<in+1<<" -> "<<lambda<<endl; //" all eigV"<<eigV<<endl;
    
    ofstream* out;
    if(in==0)
      out=new ofstream(outfile.data());
    else
      out=new ofstream(outfile.data(),ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfile<<
	" for output"<<endl;
      exit(1);
    }
    if(in==0) *out<<"% n \t M \t lambda\t coeffs"<<endl;
    *out<<setprecision(10);
    *out<<in+1<<"\t"<<M<<"\t"<<lambda<<"\t";
    for(int k=0;k<totSz;k++)
      if(k<effH2.getDimension(0)) *out<<abs(U.getElement(Indices(k,0)))<<"\t";    
      else *out<<0<<"\t";    
    *out<<endl;
    out->close();
    delete out;
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

// same if H12 is different in the last position
void prepareMPOHn(const vector<int>& indices,const MPO& H12,const MPO& H12last,MPO& Hn,int M){
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
    const MPO* miH12=&H12;
    if(indices[in]==M-2) miH12=&H12last;
    for(int l=0;l<=1;l++){
      int pos=indices[in]+l;
      mwArray newOp=miH12->getOp(l).getFullData();
      if(occupied[pos]){
	const Operator* toJoin[]={&Hn.getOp(pos),&miH12->getOp(l)};
	JoinedOperator auxOp(2,toJoin);
	newOp=auxOp.getFullData();
      }
      Hn.setOp(pos,new Operator(newOp),true);
      occupied[pos]=true;
    }
  }
}


bool prepareHermitianMPOHn(const vector<int>& indices,const MPO& H12,const MPO& H12last,MPO& Hn,int M){
  //  cout<<"Preparing hermitian Hn for indices "<<indices<<endl;
  Hn.clear();Hn.initLength(M);
  // initialize with all Ids
  mwArray id2=identityMatrix(2);
  for(int k=0;k<M;k++){
    Hn.setOp(k,new Operator(reshape(id2,Indices(2,1,2,1))),true);
  }
  vector<bool> occupied(M,false);
  bool overlaps=false; // whether some terms overlap, then I need to construct the 
  int n=indices.size(); // how many terms
  for(int in=0;in<n;in++){
    const MPO* miH12=&H12;
    if(indices[in]==M-2) miH12=&H12last;
    for(int l=0;l<=1;l++){
      int pos=indices[in]+l;
      mwArray newOp=miH12->getOp(l).getFullData();
      if(occupied[pos]){
	overlaps=true;
	const Operator* toJoin[]={&Hn.getOp(pos),&miH12->getOp(l)};
	// cout<<setprecision(10);
	// cout<<"Joining for pos "<<pos<<" Ops: "<<endl;
	// putForMatlab(cout,Hn.getOp(pos).getFullData(),"op1");
	// putForMatlab(cout,miH12->getOp(l).getFullData(),"op2");
	// cout<<endl;
	JoinedOperator auxOp(2,toJoin);
	newOp=auxOp.getFullData();
	// cout<<"Result; ";
	// putForMatlab(cout,newOp,"joined");
	// cout<<endl;
      }
      Hn.setOp(pos,new Operator(newOp),true);
      occupied[pos]=true;
    }
  }
  // If there were overlaps, I just compose with the hermitian conjugate doubling the bond dimension of the MPO. 
  if(overlaps){
    // ostringstream name;name<<"Hn_";for(int k=0;k<n;k++)name<<indices[k]<<"_";name<<"NonHerm.m";
    // Hn.exportForMatlab(name.str().data());
    // cout<<"Changing the MPO for Hn of "<<indices<<" to its hermitian conjugate. Saving first to "<<name.str()<<endl;
    for(int k=0;k<M;k++){
      mwArray aux=Hn.getOp(k).getFullData();
      Indices dims=aux.getDimensions();
      int newDl=k==0?1:2*dims[1];
      int newDr=k==M-1?1:2*dims[3];
      mwArray aux2(aux);
      aux2.resize(Indices(dims[0],newDl,dims[2],newDr)); // for both blocks
      for(int d1=0;d1<dims[0];d1++)
	for(int d2=0;d2<dims[0];d2++){
	  int offsetl=k==0?0:dims[1];
	  int offsetr=k==M-1?0:dims[3];
	  //cout<<"For pos "<<k<<" in the MPO, offsets("<<offsetl<<","<<offsetr<<")"<<endl;
	  for(int Dl=0;Dl<dims[1];Dl++)
	    for(int Dr=0;Dr<dims[3];Dr++)
	      aux2.setElement(conjugate(aux.getElement(Indices(d2,Dl,d1,Dr))),Indices(d1,Dl+offsetl,d2,Dr+offsetr));
	}
      Hn.setOp(k,new Operator(aux2),true);
    }
  }
  return overlaps;
}




void prepareOperatorsHn(int n,const MPO& H12,int M,vector<MPS*>& terms){ 
  //  int M=Omps.getLength();
  MPO Hn(M);
  mwArray id2=identityMatrix(2);
  MPS mpsId(M,1,4);mpsId.setProductState(p_maxent);
  mpsId.gaugeCond('R',true);
  vector<int> indices(n,0); // all zeros: this term has range 0, as all terms act on the same (first) position
  bool done=false;
  bool lastOne=false; // when it cannot increment indices any more
  Contractor& contractor=Contractor::theContractor();
  //double sumLambda=0;
  for(int k1=0;k1<M-n;k1++){ // only consecutive operators, e.g. 123
    indices[0]=k1;
    for(int ki=1;ki<n;ki++)indices[ki]=k1+ki;
    prepareHermitianMPOHn(indices,H12,H12,Hn,M);
    // transform in MPS to compute overlap
    MPS* aux=new MPS(M,1,1);
    MPSfromMPO(Hn,*aux,true);
    aux->gaugeCond('R',true);
    terms.push_back(aux);
  }
}
 
 void prepareCoefficients(int n, int M, vector<double>& coeffs){
   coeffs.clear();
   for(int k=0;k<M-n;k++)
     coeffs.push_back(cos(M_PIl*(k+1-(M-n+1)*.5)/(M-n)));
 }

void prepareCoefficientsrange(int n, int r,int M, vector<double>& coeffs){
   coeffs.clear();
   for(int k=0;k<M-r-1;k++)
     coeffs.push_back(cos(M_PIl*(k+1-(M-r)*.5)/(M-r-1)));
 }

// decide whether the given indices are within the 
bool indicesWithinLimits(int exactOrder,int exactR,const vector<int>& indices){
  //cout<<"Checking indices "<<indices<<" for n="<<exactOrder<<" r="<<exactR<<endl;
  if(indices.size()!=exactOrder) return false;
  vector<int>::const_iterator it=indices.begin();
  int leftMost=*it;int rightMost=*it+1;
  //cout<<"leftM="<<leftMost<<", rightM="<<rightMost<<endl;
  while(++it!=indices.end()){
    if(*it<leftMost) leftMost=*it;
    if(*it>=rightMost) rightMost=*it+1;
    //cout<<"After term "<<*it<<"leftM="<<leftMost<<", rightM="<<rightMost<<endl;
  }
  //  if(indices.back()-indices.front()+1>maxR) return false;
  if(rightMost-leftMost!=exactR) return false;
  return true; 
}


void prepareOperatorsHnrange(int n,int exactR,const MPO& H12,int M,vector<MPS*>& terms,
			     vector<double>& coeffs){ 
  //  int M=Omps.getLength();
  cout<<"Preparing ops for n="<<n<<", r="<<exactR<<endl;
  MPO Hn(M);
  mwArray id2=identityMatrix(2);
  MPS mpsId(M,1,4);mpsId.setProductState(p_maxent);
  mpsId.gaugeCond('R',true);
  vector<int> indices(n,0); // all zeros: this term has range 0, as all terms act on the same (first) position
  bool done=false;bool overlaps=false;
  bool lastOne=false; // when it cannot increment indices any more
  Contractor& contractor=Contractor::theContractor();
  //double sumLambda=0;
  while(!done){
    bool valid=indicesWithinLimits(n,exactR,indices);
    // need to loop over all possible terms with UP TO this range and order
    //      prepareMPOHn(indices,H12,Hn,M);
    if(valid){    
      overlaps=prepareHermitianMPOHn(indices,H12,H12,Hn,M);
      // cout<<"computed Hermitian Hn indices"<<indices<<endl;
      // ostringstream name;name<<"Hn_";for(int k=0;k<n;k++)name<<indices[k]<<"_";name<<".m";
      // Hn.exportForMatlab(name.str().data());
      // transform in MPS to compute overlap
      MPS* aux=new MPS(M,1,1);
      MPSfromMPO(Hn,*aux,true);
      // And I would like to remove the identity operator! It would only be needed if there are overlaps
      // if(overlaps){
      // 	complex_t comp=contractor.contract(*aux,mpsId);
      // 	if(abs(comp)>1E-8){
      // 	  //cout<<"Removing identity component ("<<abs(comp)<<") from term "<<indices<<endl;
      // 	  MPS aux2(*aux);
      // 	  vector<const MPS*> kets;kets.push_back(&aux2);kets.push_back(&mpsId);
      // 	  vector<complex_t> coeffs;coeffs.push_back(ONE_c);coeffs.push_back(-1*comp);
      // 	  contractor.optimizeSum(kets,coeffs,*aux);
      // 	}
      // }
      aux->gaugeCond('R',true);
      terms.push_back(aux);
      //what's the leftmost site?
      int leftmost=indices[0];
      for(int p=0;p<indices.size();p++) if(indices[p]<leftmost) leftmost=indices[p];
      double expon=(leftmost+1-(M-exactR+1)*.5)*(leftmost+1-(M-exactR+1)*.5);
      //coeffs.push_back(cos(M_PIl*(indices[0]+1-(M-exactR+1)*.5)/(M-exactR))*cos(M_PIl*(indices[0]+1-(M-exactR+1)*.5)/(M-exactR))*exp(-(expon)/((M-exactR)*.25)));
      coeffs.push_back(cos(M_PIl*(leftmost+1-(M-exactR+1)*.5)/(M-exactR)));
      //cout<<"Added term "<<indices<<") pos "<<terms.size()<<"= for n="<<n<<" r="<<exactR<<endl;    
    }
    valid=false;
    while(!lastOne&&!valid){
      int ptr=n-1;
      //lastOne=!incrementRecursive(indices,M-2,ptr);
      lastOne=!incrementRecursiveAll(indices,M-2,ptr);
      if(!lastOne) valid=indicesWithinLimits(n,exactR,indices);
      //if(!lastOne&&!valid) cout<<"---indices "<<indices<<" not valid"<<endl;
    }
    done=lastOne;
    //cout<<"incrementRecursiveAll+condition (maxR<="<<maxR<<") returns "<<indices<<" (done="<<done<<")"<<endl;
  }
}

double modulation(int M,int R,int pos){
  //cout<<"modulation M="<<M<<" R="<<R<<" pos="<<pos;
  return cos(M_PIl*(pos+1.-(M-R+1)*.5)/(M-R));
}

void translateMPS(const MPS& mpsR,int M,vector<MPS*>& opBasis,vector<double>& coeffs){
  int r=mpsR.getLength();
  //cout<<"Translations of MPS of range "<<r<<" over M="<<M<<"sites"<<endl;
  int d=2;
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  MPS allId(M,1,d*d);
  for(int k=0;k<M;k++){
    allId.setA(k,id2);
  }
  for(int k=0;k<M-r+1;k++){ // where I can start placing my MPS
    //cout<<" Starting on site "<<k<<", replacing: ";
    MPS* term=new MPS(allId); // first copy, the replace
    for(int l=0;l<r;l++){
      term->replaceSite(k+l,mpsR.getA(l),0);
      //cout<<k+l<<", ";
    }
    double val=modulation(M,r-1,k);
    //cout<<"- coeff "<<val<<endl;
    opBasis.push_back(term);
    coeffs.push_back(val);
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
