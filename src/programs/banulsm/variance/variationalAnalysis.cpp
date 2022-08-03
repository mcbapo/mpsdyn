
#include <math.h>
#include <iomanip>

#include "mwArray.h"
#include "MPS.h"
#include "misc.h"
#include "Contractor.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"
#include "HeisenbergHamiltonian.h"

using namespace shrt;

/** Runs the minimizeVariance routine with the MPO for the Ising
    Hamiltonian \ref <IsingHamiltonian> at a given energy, and saves
    the results (MPS) in a file.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <E0> (double) energy value
    \param <penH> (double) penalty fort the energy term
    \param <outfname> (char*) name of the output file for the MPS
*/

const string mpsFileName(const string& mpsdir,const int N,const double E0,const int D);
int findInitFile(string& str,const string& initmpsdir,const int N,const double E0,const int D);
const string jobfilename(int L,int D,double E0,const string& strMod,const string jobsdir,int parNr);
double getSchmidtVals(vector<double>& schmidt,const MPS& mps,int pos,char gaugeDir);
void expandRDM(const MPS& mps,int posL,int posR,mwArray& A);
double computeRDMdistanceToId(const mwArray& rdm,int Lc);

int computeBlockEntropiesExact(MPS& mps,vector<double>& trDistRDM,
			   vector<double>& entropyRDM,vector<int>& sizeRDM);

void computeTwoBodyTerm(mwArray& result,double Jx,double Jy,double Jz,
			double hx1,double hx2,double hy1,double hy2,double hz1,double hz2);

void splitTwoBodyTerms(const mwArray& h12,mwArray& Ol,mwArray& Or,mwArray& Omid12,mwArray& Omid21);


void computeEnergyEnergy(const MPS& mps,mwArray& h12,vector<int>& pos,vector<double>& h0hm,vector<double>& hm,int offset=0);


int d=2;
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
  double J_=props.getDoubleProperty("J");
  double h_=props.getDoubleProperty("h");
  double g_=props.getDoubleProperty("g");
  bool ising=(props.getIntProperty("ising")!=0); // only if -ising=0 is given it will be Heisenberg (or XY)  
  // These are only used if Heisenberg model
  double Jx_ = props.getDoubleProperty("Jx");  
  double Jy_ = props.getDoubleProperty("Jy");  
  double Jz_ = props.getDoubleProperty("Jz");
  // It can also include x and y local fields Heisenberg model
  double hx_ = props.getDoubleProperty("hx");  
  double hy_ = props.getDoubleProperty("hy");  

  double E0=props.getDoubleProperty("E");
  double maxDevE=props.getDoubleProperty("toleranceE");
  if(maxDevE<0) maxDevE=max(1E-2,.1*E0);
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  int D_=props.getIntProperty("D");
  int incrD=props.getIntProperty("incrD");
  int maxD=props.getIntProperty("maxD");
  string Sredfname=props.getProperty("outputSblock"); // for the entropy of the middle blocks
  string corrfname=props.getProperty("outputCorr"); // for the energy correlations
  string mpsdir=props.getProperty("mpsdir");
  // int Lb=props.getIntProperty("Lb");
  // if(Lb<0) Lb=6;
  // bool avrg=1;
  // avrg=!(props.getIntProperty("average")==0); // whether to average for each Lc subsystem
  bool app=(props.getIntProperty("append")!=0);

  bool isXY=(!ising&&Jz_==0);
  int parNr=props.getIntProperty("parNr"); // Just for job identifier
  cout<<"Initialized arguments: L="<<L;
  if(ising){
    cout<<"; ISING: J="<<J_
	<<", g="<<g_
	<<", h="<<h_;
  }
  else{
    if(!isXY){
      cout<<"; HEISENBERG: Jx="<<Jx_
	  <<", Jy="<<Jy_
	  <<", Jz="<<Jz_
	  <<", h="<<h_;
    }
    else{ // XY case
      cout<<"; XY: Jx="<<Jx_
	  <<", Jy="<<Jy_
	  <<", hx="<<hx_
	  <<", hy="<<hy_
	  <<", hz="<<h_;
    }    
  }
  cout<<", E="<<E0
      <<", app="<<app
      <<", D="<<D_<<endl;

  if((!file_exists(Sredfname.data()))||(!file_exists(corrfname.data()))) app=0; // also if it does not exist (new)
  ofstream* outSred;
  ofstream* outEcorr;
  if(!app){
    outSred=new ofstream(Sredfname.data());
    if(!outSred->is_open()){
      cout<<"Error: impossible to open file "<<Sredfname<<
	" for output"<<endl;
      exit(1);
    }
    if(ising){
      *outSred<<"% N\t J\t g\t h\t D\t Lc\t Sred(Lc)\t trDist"<<endl;
    }
    else{
      if(!isXY)
	*outSred<<"% N\t Jx\t Jy\t Jz\t h\t D\t Lc\t Sred\t trDist"<<endl;
      else
	*outSred<<"% N\t Jx\t Jy\t hx\t hy\t hz\t D\t Lc\t Sred(Lc)\t trDist"<<endl;
    }
    outSred->close();
    delete outSred;

    outEcorr=new ofstream(corrfname.data());
    if(!outEcorr->is_open()){
      cout<<"Error: impossible to open file "<<corrfname<<
	" for output of energy correls"<<endl;
      exit(1);
    }
    if(ising){
      *outEcorr<<"% N\t J\t g\t h\t D\t E\t <H^2>\t x\t <h0 h(x)>\t <h(x)>"<<endl;
    }
    else{
      if(!isXY)
	*outEcorr<<"% N\t Jx\t Jy\t Jz\t h\t D\t E\t <H^2>\t x\t <h0 h(x)>\t <h(x)>"<<endl;
      else
	*outEcorr<<"% N\t Jx\t Jy\t hx\t hy\t hz\t D\t E\t <H^2>\t x\t <h0 h(x)>\t <h(x)>"<<endl;
    }
    outEcorr->close();
    delete outEcorr;        
  }

  cout<<setprecision(10);

  Hamiltonian* ham0;
  if(ising){
    ham0=new IsingHamiltonian(L,d,J_,g_,h_);
  }
  else{
    ham0=new HeisenbergHamiltonian(L,Jx_,Jy_,Jz_,hx_,hy_,h_);
  }
  //IsingHamiltonian hamH0(L,d,J_,g_,h_);
  const MPO& hamil=ham0->getHMPO();
  cout<<"Constructed the hamil MPO"<<endl;
  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  cout<<"Initialized Contractor"<<endl;
  mwArray h12,hlast;
  if(ising){
    ham0=new IsingHamiltonian(L,d,J_,g_,h_);
    //    computeTwoBodyTerm(h12,0.,0.,J_,.5*g_,.5*g_,0.,0.,.5*h_,.5*h_);
    computeTwoBodyTerm(h12,0.,0.,J_,g_,0.,0.,0.,h_,0.); // all single site on left
  }
  else{
    ham0=new HeisenbergHamiltonian(L,Jx_,Jy_,Jz_,hx_,hy_,h_);
    //    computeTwoBodyTerm(h12,.25*Jx_,.25*Jy_,.25*Jz_,.25*hx_,.25*hx_,.25*hy_,.25*hy_,.25*h_,.25*h_);
    computeTwoBodyTerm(h12,.25*Jx_,.25*Jy_,.25*Jz_,.5*hx_,0.,.5*hy_,0.,.5*h_,0.);
  }

  
  // Iterate for all Ds
  for(int D=D_;D<=maxD;D+=incrD){  
    string mpsfile=mpsFileName(mpsdir,L,E0,D);
    MPS mps(L,D,d);
    // FIRST: Check if the same file exisits (then refine)
    if(!file_exists(mpsfile)){
      cout<<"No MPS file saved for parameters L="<<L<<" D="<<D<<endl;
    }
    else{
      mps.importMPS(mpsfile.data());
      mps.gaugeCond('R',1);mps.gaugeCond('L',1); // just in case
      cout<<"Imported initial state from file "<<mpsfile<<endl;
  
      // place for other params
      complex_t E=contractor.contract(mps,hamil,mps);      
      complex_t E2=contractor.contract2(hamil,mps);
      //double varE=real(E2)-real(E)*real(E);
      cout<<"Variance="<<E2-E*E<<" <H^2>="<<E2<<", <H>="<<E<<endl;
      // Compute and save same things as with Cheby

      double Shalf=contractor.getEntropy(mps);

      // entropies and distances to thermal of central blocks 
  
      vector<double> trDistRDM;
      vector<double> entropyRDM;
      vector<int> sizeRDM;
  
      // Also energy-energy correlations
      vector<double> h0hm(1,0.);vector<double> hm(1,0.);
      vector<int> posM(1,0);


      int offset=0;
      if(!ising){//&&abs(initSt)==10){
	if(L/2%4!=0)offset=2;
      }
      computeEnergyEnergy(mps,h12,posM,h0hm,hm,offset);
      // Write in file
      outEcorr=new ofstream(corrfname.data(),ios::app);
      if(!outEcorr->is_open()){
	cout<<"Error: impossible to open file "<<corrfname<<
	  " for output"<<endl;
	exit(1);
      }
      *outEcorr<<setprecision(15);
      for(int k=0;k<posM.size();k++){
	// OJO!!!!! This is a mistake, doing it wrong for Heisenberg!!!!
	if(ising)
	  *outEcorr<<L<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t";
	else
	  if(!isXY)
	    *outEcorr<<L<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t";
	  else
	    *outEcorr<<L<<"\t"<<Jx_<<"\t"<<Jy_<<"\t"<<hx_<<"\t"<<hy_<<"\t"<<h_<<"\t"<<D<<"\t";
	*outEcorr<<real(E)<<"\t"<<real(E2)<<"\t";
	*outEcorr<<posM[k]<<"\t"<<h0hm[k]<<"\t"<<hm[k]<<"\t";
	*outEcorr<<endl;
      }
      outEcorr->close(); delete outEcorr;

      //  int Dcut=computeBlockEntropiesSlow(mps,trDistRDM,entropyRDM,sizeRDM);
      int Dcut=D;computeBlockEntropiesExact(mps,trDistRDM,entropyRDM,sizeRDM);
      // Write in file
      outSred=new ofstream(Sredfname.data(),ios::app);
      if(!outSred->is_open()){
	cout<<"Error: impossible to open file "<<Sredfname<<
	  " for output"<<endl;
	exit(1);
      }
      *outSred<<setprecision(15);
      for(int k=0;k<sizeRDM.size();k++){
	// OJO!!!!! This is a mistake, doing it wrong for Heisenberg!!!!
	if(ising)
	  *outSred<<L<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t";
	else
	  if(!isXY)
	    *outSred<<L<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t"; // wrong!!!!
	  else
	    *outSred<<L<<"\t"<<Jx_<<"\t"<<Jy_<<"\t"<<hx_<<"\t"<<hy_<<"\t"<<h_<<"\t"<<D<<"\t";
	*outSred<<sizeRDM[k]<<"\t"<<entropyRDM[k]<<"\t"<<trDistRDM[k];
	*outSred<<endl;
      }
      outSred->close(); delete outSred;
    }
  }
}





int findInitFile(string& str,const string& initmpsdir,const int N,const double E0,const int D){
  cout<<"Searching for init file for D "<<D<<" in dir "<<initmpsdir<<endl;
  for(int myD=D-1;myD>0;myD--){
    // TRy with one smaller D
    str=mpsFileName(initmpsdir,N,E0,myD);
    if(file_exists(str)){
      cout<<"Found file "<<str<<endl;
      return myD;
    }
  }
  // if nothing found
  str="";
  return 0; 
}

const string mpsFileName(const string& mpsdir,const int N,const double E0,const int D){
  stringstream s;
  s<<mpsdir<<"/MPS_N"<<N<<"_E0"<<E0<<"_D"<<D;
  return s.str();
}

double getSchmidtVals(vector<double>& schmidt,const MPS& mps,int pos,char gaugeDir){
  Contractor& contractor=Contractor::theContractor();
  int len=mps.getLength();
  mwArray rdm=contractor.contract(mps,mps,pos,gaugeDir);
  // cout<<"getSchmidtVals of mps:"<<mps<<" at cut pos="<<pos<<"\nComputed the rdm: "<<rdm.getDimensions()<<endl;
  // Now diagonalize
  rdm.permute(Indices(2,1)); // here it does not matter
  vector<complex_t> Dval;mwArray U; //placeholder:ignored
  int D=rdm.getDimension(0); //mps.getA(pos).getDl();
  wrapper::eig(rdm,Dval,U);
  schmidt.clear();
  for(int k=0;k<D;k++) schmidt.push_back(real(Dval[k]));
  double entropy=0.;
  for(int k=0;k<D;k++){
    double tmp=real(Dval[k]);
    if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
  }
  //cout<<"getSchmidtVals computed vals="<<Dval<<"\n and Entropy="<<entropy<<endl;
  return entropy;
  
}

 double computeRDMdistanceToId(const mwArray& rdm,int Lc){
   double singvals=0.;
   double idFc=1./pow(2,Lc);
   mwArray U;vector<complex_t> Dv;
   wrapper::eig(rdm,Dv,U); // only eigenvalues
   // now sum the absolute values
   for(int k=0;k<Dv.size();k++){
     double aux=real(Dv[k]);
     singvals+=abs(aux-idFc);
   }
   if(Dv.size()<pow(2,Lc)){
     cout<<"WARNING! Seems the dimension of the operator for Lc="<<Lc
	 <<"is not full: dim(oper)="<<Dv.size()
	 <<", 2^Lc="<<pow(2,Lc)<<endl;
     singvals+=idFc*(pow(2,Lc)-Dv.size());
   }
   return singvals;
 }


void computeTwoBodyTerm(mwArray& result,double Jx,double Jy,double Jz,
			double hx1,double hx2,double hy1,double hy2,double hz1,double hz2){
  mwArray sig0=identityMatrix(d);//identity
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),datax);//sigmax
  complex_t datay[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigY(Indices(d,d),datay);//sigmay
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz
  
  mwArray XX,YY,ZZ,XI,IX,YI,IY,ZI,IZ;
  constructOperatorProduct(XX,sigX,sigX);
  constructOperatorProduct(YY,sigY,sigY);
  constructOperatorProduct(ZZ,sigZ,sigZ);
  constructOperatorProduct(XI,sigX,sig0);
  constructOperatorProduct(YI,sigY,sig0);
  constructOperatorProduct(ZI,sigZ,sig0);
  constructOperatorProduct(IX,sig0,sigX);
  constructOperatorProduct(IY,sig0,sigY);
  constructOperatorProduct(IZ,sig0,sigZ);

  result=Jx*XX+Jy*YY+Jz*ZZ+hx1*XI+hx2*IX+hy1*YI+hy2*IY+hz1*ZI+hz2*IZ;
  result.reshape(Indices(d,d,d,d));
  result.permute(Indices(1,3,2,4));
  result.reshape(Indices(d*d,d*d));
}

void splitTwoBodyTerms(const mwArray& h12,mwArray& Ol,mwArray& Or,mwArray& Omid12,mwArray& Omid21){
  mwArray S; // place for the singular values
  int nr=0;
  double tol=0.;
  //cout<<"Before decomposition, expH="<<expH<<endl;
  wrapper::svd(h12,tol,nr,Ol,S,Or);
  //  cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<", S="<<S<<endl;
  // redistribute the singular values
  mwArray sqrS=sqrt(S);
  // TODO: check no NaN???
  int xi=S.getDimension(0);
  Ol.multiplyRight(sqrS); // d*d*Xi
  Or.multiplyLeft(sqrS); // Xi*d*d
  // reshape for operators (i'x 1 x i x alpha ))
  Ol.reshape(Indices(d,1,d,-1));
  Or.reshape(Indices(-1,d,d,1));
  Or.permute(Indices(2,1,3,4));  

  Omid12=mwArray(Or); // d' xi d 1
  Omid12.reshape(Indices(d*xi,d));
  Omid12.multiplyRight(reshape(Ol,Indices(d,d*xi))); // d' *xi x d * xi
  Omid12.reshape(Indices(d,xi,d,xi));

  Omid21=mwArray(Ol); // d'*1*d*xi
  Omid21.permute(Indices(1,4,3,2)); //d' xi d 1
  Omid21.reshape(Indices(d*xi,d));
  Omid21.multiplyRight(reshape(Or,Indices(d,xi*d))); // d'*xi x xi(l)*d
  Omid21.reshape(Indices(d,xi,xi,d));
  Omid21.permute(Indices(1,3,4,2));
  
}

void computeEnergyEnergy(const MPS& mps,mwArray& h12,vector<int>& pos,vector<double>& h0hm,vector<double>& hm,int offset){
  cout<<"computeEnergyEnergy with h12="<<h12<<endl;
  h0hm.clear();hm.clear();pos.clear();
  // split the two body op as MPO
  mwArray Ol,Or,Omid12,Omid21;
  splitTwoBodyTerms(h12,Ol,Or,Omid12,Omid21);

   // // cout<<"After splitTwoBodyTerms,"<<endl;
   //  putForMatlab(cout,Ol,"Ol");
   //  putForMatlab(cout,Or,"Or");
   //  putForMatlab(cout,Omid12,"Omid12");
   //  //exit(1);

  Operator* Olop=new Operator(Ol);
  Operator* Orop=new Operator(Or);
  Operator* Omid=new Operator(Omid12);
  
  int L=mps.getLength();
  MPO corrHH(L);
  MPO corrHH1(L); // just for the first term
  MPO singleH(L);
  Operator* idOp=new Operator(reshape(identityMatrix(d),Indices(d,1,d,1)));
  for(int k=0;k<L;k++){
    corrHH.setOp(k,idOp,0);
    corrHH1.setOp(k,idOp,0);
    singleH.setOp(k,idOp,0);
  }
  // Now all MPOs are just identity
  corrHH1.setOp(L/2-1-offset,Olop,0);
  corrHH1.setOp(L/2-offset,Omid,0);
  corrHH1.setOp(L/2+1-offset,Orop,0);

  Contractor& contractor=Contractor::theContractor();

  // first the e.v. of h12 on the middle pair (L/2-1 and L/2)
  corrHH.setOp(L/2-1-offset,Olop,0);
  corrHH.setOp(L/2-offset,Orop,0);
  double hm_=real(contractor.contract(mps,corrHH,mps));
  double h0hm_=real(contractor.contract2(corrHH,mps));

  // cout<<"In this state, <h12>="<<hm_<<" and <h12^2>="<<h0hm_<<endl;
  // cout<<"State is "<<mps<<endl;
  
  pos.push_back(L/2-1-offset);h0hm.push_back(h0hm_);hm.push_back(hm_);
  // the site L/2 is special, too
  singleH.setOp(L/2-offset,Olop,0);
  singleH.setOp(L/2+1-offset,Orop,0);
  hm_=real(contractor.contract(mps,singleH,mps));
  singleH.setOp(L/2-offset,idOp,0);
  singleH.setOp(L/2+1-offset,idOp,0);

  h0hm_=real(contractor.contract(mps,corrHH1,mps));
  pos.push_back(L/2);
  hm.push_back(hm_);
  h0hm.push_back(h0hm_);

  for(int k=L/2+1-offset;k<L-1;k++){
    singleH.setOp(k,Olop,0);
    singleH.setOp(k+1,Orop,0);
    hm_=real(contractor.contract(mps,singleH,mps));    
    singleH.setOp(k,idOp,0);
    singleH.setOp(k+1,idOp,0);
    corrHH.setOp(k,Olop,0);
    corrHH.setOp(k+1,Orop,0);
    h0hm_=real(contractor.contract(mps,corrHH,mps));
    corrHH.setOp(k,idOp,0);
    corrHH.setOp(k+1,idOp,0);
    pos.push_back(k);hm.push_back(hm_);h0hm.push_back(h0hm_);
 
  }
  cout<<"At the end of computeEnergyEnergy, pos="<<pos<<endl;
  cout<<"At the end of computeEnergyEnergy, hm="<<hm<<endl;
  cout<<"At the end of computeEnergyEnergy, hmho="<<h0hm<<endl;
  delete idOp,Olop,Orop,Omid;
}

// for at most Lc=12 (max I can do exactly, I think), compute the
// exact RDM of the central Lc sites, its entropy and distance to the maximally mixed state.
// It could be possible to use a numerical approximation to the entropy for something larger

// Assuming proper "central" gauge
void expandRDM(const MPS& mps,int posL,int posR,mwArray& A){
  A=mps.getA(posL).getA();
  int d=A.getDimension(0);int Dl0=A.getDimension(1);int Dr=A.getDimension(2);
  A.permute(Indices(2,1,3));A.reshape(Indices(Dl0*d,Dr));
  int dphys=d; // temporary dimensions of the vector
  int Dv=Dr;
  for(int k=posL+1;k<=posR;k++){
    mwArray aux=mps.getA(k).getA();
    //cout<<"Multiplying tensor for site "<<k<<", A="<<aux.getDimensions()<<" with tmp result "<<A.getDimensions()<<endl;
    Indices dims=aux.getDimensions();
    int d=dims[0];int Dl=dims[1];int Dr=dims[2];
    if(Dl!=Dv){
      cout<<"Error: Dimensions do not agree in expandVec!"<<endl;
      exit(1);
    }
    aux.permute(Indices(2,1,3));
    aux.reshape(Indices(Dl,d*Dr));
    A.multiplyRight(aux);
    dphys=dphys*d;
    Dv=Dr;
    A.reshape(Indices(Dl0*dphys,Dv));
  }
  A.reshape(Indices(Dl0,dphys,Dv));A.permute(Indices(2,1,3));A.reshape(Indices(dphys,Dl0*Dv));
  // if the mps has a normalization factor, we need to include it!!
  A=mps.getNormFact()*A;
  A.multiplyRight(Hconjugate(A));
}  


int computeBlockEntropiesExact(MPS& mps_,vector<double>& trDistRDM,
			   vector<double>& entropyRDM,vector<int>& sizeRDM){
  int Lcmax=10;
  cout<<"computeBlockEntropiesExact "<<endl;
  int L=mps_.getLength();
  int minPosR=L/2; // assuming L even and smallest Lc is 2
  mps_.gaugeCond('R',1);
  for(int k=L-1;k>=minPosR+1;k--){
    mps_.gaugeCond(k,'L',true);
  } // now both sides are identity
  int Lc=Lcmax;
  int posL=(L-Lc)/2;int posR=posL+Lc-1; // outside tracing them
  while(Lc>=2){
    mwArray oper;
    expandRDM(mps_,posL,posR,oper);
    cout<<"Created RDM for Lc="<<Lc<<endl;
    double singvals=0.;
    double idFc=1./pow(2,Lc);
    mwArray U;vector<complex_t> Dv;
    wrapper::eig(oper,Dv,U); // only eigenvalues
    double entropy=0.;
    for(int ik=0;ik<Dv.size();ik++){
      double tmp=real(Dv[ik]);
      if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
    }
    cout<<"Diagonalized the rho(Lc)"<<endl;
    // now sum the absolute values
    for(int k=0;k<Dv.size();k++){
      double aux=real(Dv[k]);
      singvals+=abs(aux-idFc);
    }
    if(Dv.size()<pow(2,Lc)){
      cout<<"WARNING! Seems the dimension of the RDM operator for D="<<mps_.getA(L/2).getDl()
	  <<", Lc="<<Lc<<"is not full: dim(oper)="<<Dv.size()
	  <<", 2^Lc="<<pow(2,Lc)<<endl;
      singvals+=idFc*(pow(2,Lc)-Dv.size());
    }
    sizeRDM.push_back(Lc);
    trDistRDM.push_back(singvals);
    entropyRDM.push_back(entropy);
    posL++;posR--;Lc=(posR-posL+1);
  }
  
}

