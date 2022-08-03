
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#include "SpinMPO.h"
#include "FAHamiltonian.h"
#include "Properties.h"
#include "misc.h"

//#include "quicksort.h"

#define FREQSAVE 18000 // frequency of tmp saving: default 5 hrs 


/** 
Generates a file name for the MPS corresponding to a certain excitation.
*/
const string mpsfilename(int L,int D,double sL,double c,int level,const string mpsdir);


void getProjectorTotalOccup(int L,MPO& PNz,mwArray& match);
void getProjectorTotalOccupX(int L,MPO& PNx,mwArray& matchX);

/** faDistributions takes the states computed by FASpec and
    computes observables on them.  In this case, the probability
    distributions of occupations in Z and X, and the entropy of half chain.

    Receives as argument a Properties file, such as
    config/stochspec.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.

    Receives, among others, arguments:
    \param <L> (int) length of the chain
    \param <c> (double) parameter \f$c\f$ of the Hamiltonian
    \param <s> (double) parameter \f$s\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)

*/

using namespace shrt;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  int L=props.getIntProperty("L");
  double sL=props.getDoubleProperty("sL");
  double s=sL/(1.*L);
  int D0=props.getIntProperty("D");
  int maxD=props.getIntProperty("maxD");
  int incrD=props.getIntProperty("incrD");
  if(incrD<0) incrD=20;
  double c=props.getDoubleProperty("c");
  if(c<0||c>1){
    cout<<"Error: value of c not allowed (c in [0,1])"<<endl;
    exit(1);
  }
  int nlev=props.getIntProperty("nlevel");
  string outfname=props.getProperty("outputfile");
  string outfnameS=props.getProperty("outputfileS"); // for the entropy
  string mpsdir=props.getProperty("mpsdir");
  int append=props.getIntProperty("append");
  bool app=(append!=0);
  cout<<"Appending to output?"<<app<<endl;  
  double penZero=props.getDoubleProperty("penaltyZero");
  if(penZero<0) penZero=0.;
  bool pbc=(props.getIntProperty("pbc")==1);

  bool origApp=app;
  ofstream* out;
  if(!app||!file_exists(outfname)){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% Distribution of occupations of FA model L="<<L<<", c="<<c<<endl;
    *out<<"% sL="<<sL<<endl;
    *out<<"% nlev\t D\t (1=X/3=Z)\t<H>\t<H^2>\t P(m=0)\t P(m=1)..."<<endl;
    out->close();delete out;
  }

  ofstream* outS;
  if(!app||!file_exists(outfnameS)){
    outS=new ofstream(outfnameS.data());
    if(!outS->is_open()){
      cout<<"Error: impossible to open file "<<outfnameS<<
	" for output"<<endl;
      exit(1);
    }
    *outS<<"% Entropy of FA model L="<<L<<endl;
    *outS<<"% c="<<c<<endl;
    *outS<<"% sL="<<sL<<endl;
    *outS<<"% nlev\t D\t <H>\t<H^2>\t pos(L/2)\t S(L/2)"<<endl;
    outS->close();delete outS;
  }

  //out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
  //out->close();delete out;
    
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;

  int d=2;

  double offset=0.;
  FAHamiltonian hSto(L,s,offset,c,penZero,pbc);
  const MPO& hamil=hSto.getHMPO();
  //hamil.exportForMatlab("testSto.m");//exit(1);

  // Operators I want to compute:
  MPO PNz(L);
  mwArray match;
  getProjectorTotalOccup(L,PNz,match);
  MPO PNx(L);
  mwArray matchX; // this is the same one!
  getProjectorTotalOccupX(L,PNx,matchX);

  // Now loop over Ds
  int D=D0;
  while(D<=maxD){
  
    // REad the existing computed levels
    bool oldMPS=true;
    int lval=0;

    while(oldMPS&&lval<=nlev-1){
      const string filename=mpsfilename(L,D,sL,c,lval,mpsdir);
      // Try to open file for MPS level lval
      if(file_exists(filename)){
	MPS nextLev(L,D,d);
	nextLev.importMPS(filename.data());
	cout<<"Recovered MPS file for level "<<lval<<" from file "<<filename<<endl;
	nextLev.gaugeCond('R',1);
	nextLev.gaugeCond('L',1);
	complex_t Ek=contractor.contract(nextLev,hamil,nextLev);
	complex_t Ek2=contractor.contract2(hamil,nextLev);
	double entropy=contractor.getEntropy(nextLev,L/2);
	
	outS=new ofstream(outfnameS.data(),ios::app);*outS<<setprecision(15);
	*outS<<lval<<"\t"<<D<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t";
	*outS<<L/2-1<<"\t"<<entropy<<"\t";      
	*outS<<endl;
	outS->close();delete outS;

	
	for(int ind=0;ind<2;ind++){
	  MPO& PN=(ind==0)?PNz:PNx;
	  // place for the results
	  mwArray contrL,contrR;
	  int posC=(L%2==0)?L/2:(L+1)/2;
	  contrL=contractor.contract(nextLev,PN,nextLev,posC,'L'); // Du x Xi x Dd
	  contrR=contractor.contract(nextLev,PN,nextLev,posC-1,'R'); // Du x Xi x Dd
	  Indices dimsL=contrL.getDimensions();
	  Indices dimsR=contrR.getDimensions();
	  // Check??
	  contrL.permute(Indices(2,1,3)); // XiL Du Dd
	  contrL.reshape(Indices(dimsL[1],dimsL[0]*dimsL[2]));
	  contrR.permute(Indices(1,3,2)); // Du Dd XiR
	  contrR.reshape(Indices(dimsR[0]*dimsR[2],dimsR[1]));
	  contrL.multiplyRight(contrR); // XiL x XiR
	  contrL.reshape(Indices(1,dimsL[1]*dimsR[1]));
	  contrL.multiplyRight(match);
	  
	  cout<<"Computed distribution for ind="<<ind<<" contrL="<<contrL.getDimensions()<<endl;
	  
	  // Write out results      
	  out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
	  *out<<lval<<"\t"<<D<<"\t";
	  if(ind==0) *out<<3<<"\t";
	  if(ind==1) *out<<1<<"\t";
	  *out<<real(Ek)<<"\t"<<real(Ek2)<<"\t";
	  for(int m=0;m<L+1;m++)
	    *out<<real(contrL.getElement(Indices(0,m)))<<"\t";
	  *out<<endl;
	  out->close();delete out;
	  
	}	
      
	lval++;
      }
      else{
	oldMPS=false; // break the loop 
      }
    }
    D+=incrD;
  }  
}


const string mpsfilename(int L,int D,double sL,double c,int level,const string mpsdir){
  stringstream str;
  str<<mpsdir<<"/MPS_L"<<L<<"_sL"<<sL<<"_eps"<<c<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return str.str();
}

void computeOccupations(const MPS& gs,int L,int d,vector<double>& occ){
  // Construct a MPO for the particle number on a given site
  static MPO mpoId(L); // Basic MPO: all Identity
  static Operator idOp(reshape(identityMatrix(d),Indices(d,1,d,1)));
  static Operator* nOp; // individual operator number
  static bool init(false);
  if(!init){
    mwArray nSite(Indices(d,d));
    for(int k=0;k<d;k++){
      nSite.setElement(k*ONE_c,Indices(k,k));
    }
    nSite.reshape(Indices(d,1,d,1));
    nOp=new Operator(nSite);
    for(int k=0;k<L;k++)
      mpoId.setOp(k,&idOp,0);
    init=true;
  }
  Contractor& contractor=Contractor::theContractor();
  occ.clear();
  for(int pos=0;pos<L;pos++){
    mpoId.setOp(pos,nOp,false);
    complex_t nk=contractor.contract(gs,mpoId,gs);
    mpoId.setOp(pos,&idOp,false);
    occ.push_back(real(nk));
  }
}

void computePolarizations(const MPS& gs,int L,int d,vector<double>& polX){
  // Construct a MPO for the sigma_x on a given site
  static MPO mpoId(L); // Basic MPO: all Identity
  static Operator idOp(reshape(identityMatrix(d),Indices(d,1,d,1)));
  static Operator* xOp; // individual operator sigmax
  static bool init(false);
  if(!init){
    complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    mwArray sigX(Indices(d,1,d,1),dataX);
    xOp=new Operator(sigX);
    for(int k=0;k<L;k++)
      mpoId.setOp(k,&idOp,0);
    init=true;
  }
  Contractor& contractor=Contractor::theContractor();
  polX.clear();
  for(int pos=0;pos<L;pos++){
    mpoId.setOp(pos,xOp,false);
    complex_t nk=contractor.contract(gs,mpoId,gs);
    mpoId.setOp(pos,&idOp,false);
    polX.push_back(real(nk));
  }
}


void getProjectorTotalOccup(int L,MPO& PNz,mwArray& match){
  int d=2;
  PNz.initLength(L);
  // In this case, I construct the tensors by hand
  for(int k=0;k<L/2;k++){
    mwArray termK(Indices(d,k+1,d,k+2));
    for(int j=0;j<k+1;j++){
      termK.setElement(ONE_c,Indices(0,j,0,j));
      termK.setElement(ONE_c,Indices(1,j,1,j+1));
    }
    // I have to take proper care of the even and odd lengths
    PNz.setOp(k,new Operator(termK),true);
    if(L%2==0||k<(L-1)/2){
      termK.permute(Indices(1,4,3,2));
      PNz.setOp(L-k-1,new Operator(termK),true);
    }
  }
  // dimensions of the left and right parts that have to be joined by the special matrix
  int dimL=(L%2==0)?L/2+1:(L+1)/2;
  int dimR=(L%2==0)?L/2+1:(L-1)/2;
  match=mwArray(Indices(dimL,dimR,L+1));
  // When I have to connect both sides, an auxiliary intermediate
  // matrix ensures that indices match (actually, it sums them and
  // gives an extra index with the sum, so tht I get all values from
  // one contraction)
  for(int k=0;k<dimL;k++)
    for(int l=0;l<dimR;l++)
      match.setElement(ONE_c,Indices(k,l,k+l));
  // To be able to use it with the tmp contractions:
  match.reshape(Indices(dimL*dimR,L+1));
}


void getProjectorTotalOccupX(int L,MPO& PNx,mwArray& matchX){
  int d=2;
  getProjectorTotalOccup(L,PNx,matchX);
  // I only need to change the basis:
  // Since exp(i*pi/4*sigy)*sigx*exp(-i*pi/4*sigy)=sigz, I can act with
  // exp(-i*pi/4*sigy)^L PNz exp(i*pi/4*sigy)^L

  complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  mwArray sig0=identityMatrix(d);

  mwArray rotY=(1./sqrt(2))*(sig0+I_c*sigY);
  mwArray rotYdag=Hconjugate(rotY);
  for(int k=0;k<L;k++){
    mwArray aux=PNx.getOpData(k); // d Dl d Dr
    Indices dimOp=aux.getDimensions();
    aux.reshape(Indices(d,-1));
    aux.multiplyLeft(rotYdag);
    aux.reshape(dimOp);
    aux.permute(Indices(1,2,4,3));
    aux.reshape(Indices(-1,d));
    aux.multiplyRight(rotY);
    aux.reshape(Indices(dimOp[0],dimOp[1],dimOp[3],dimOp[2]));
    aux.permute(Indices(1,2,4,3));
    PNx.setOp(k,new Operator(aux),true);
  }
}
