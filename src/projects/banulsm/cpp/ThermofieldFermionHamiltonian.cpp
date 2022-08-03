#include "ThermofieldFermionHamiltonian.h"
#include "Indices.h"

using namespace std;
using namespace shrt;

//#define TRUNCESP 1

ThermofieldFermionHamiltonian::ThermofieldFermionHamiltonian(int L_,
							     const vector<double>& A1n_,
							     const vector<double>& A2n_,
							     const vector<double>& B1n_,
							     const vector<double>& B2n_,
							     const double U,const double V):
  ThermofieldHamiltonian(L_,vector<int>(2*L_+1,2),A1n_,A2n_,B1n_,B2n_),d(2){  
  cout<<"Creating TFFermionHamiltonian with arguments:"<<endl;
  // cout<<"dims="<<dims_<<endl;
   cout<<"A1n="<<A1n<<endl;
   cout<<"A2n="<<A2n<<endl;
   cout<<"B1n="<<B1n<<endl;
   cout<<"B2n="<<B2n<<endl;
  // dimensions in the parent ThermofieldHamiltonian are fixed, as these are fermionic modes
  // Probably not needed at all!
  dims[L_]=4;
  // for(int k=0;k<L;k++)
  //   dims.push_back(2);
  // dims.push_back(4);
  // for(int k=0;k<L;k++)
  //   dims.push_back(2);
  // Introduce the signs
  for(int k=0;k<L;k++){
    if(k>0)
      B2n[k]=-B2n[k];
    A2n[k]=-A2n[k];
  }
  // Now I create Hfree by hand, because form is fixed
  constructFreeHamiltonian(U,V);

  initZ();
  initHMPO();
}

void ThermofieldFermionHamiltonian::constructFreeHamiltonian(const double U,const double V){
  // The basic (single site) operators I need
  //  cout<<"Constructing free Hamiltonian"<<endl;
  mwArray sig0=identityMatrix(d);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  constructOperatorProduct(Hfree,sig0+sigZ,sig0+sigZ);
  cout<<"Hfree="<<Hfree<<endl;
  Hfree=.125*U*Hfree; //*U/8
  mwArray term1; // auxiliary for single body terms
  constructOperatorProduct(term1,.5*sigZ,sig0);
  Hfree=Hfree+V*term1;
  constructOperatorProduct(term1,sig0,.5*sigZ);
  Hfree=Hfree+V*term1+V*identityMatrix(d*d);
}

void ThermofieldFermionHamiltonian::constructOperatorProduct(mwArray& result,const mwArray& opA,
							     const mwArray& opB) const{
  Indices dimsA=opA.getDimensions();
  Indices dimsB=opB.getDimensions();
// #ifdef CHECKDIMS
//   if(opA.getDimensions()!=opB.getDimensions()){
//     cout<<"Error in ThermofieldFermionHamiltonian::constructOperatorProduct for"
// 	<<" A="<<opA<<" and B="<<opB<<endl;
//     exit(1);
//   }
// #endif
  result=opA; 
  result.reshape(Indices(dimsA[0]*dimsA[1],1));
  result.multiplyRight(reshape(opB,Indices(1,dimsB[0]*dimsB[1])));
  result.reshape(Indices(dimsA[0],dimsA[1],dimsB[0],dimsB[1])); // order here: iji'j'
  result.permute(Indices(1,3,2,4));  // order in the MPO: ii',jj'  
  result.reshape(Indices(dimsA[0]*dimsB[0],dimsA[1]*dimsB[1]));
}

void ThermofieldFermionHamiltonian::computeTwoBodyTermBath(mwArray& result,int k,
						      bool left) const{
  if(k>=L-1){
    cout<<"Error:cannot compute two body term for modes "<<k<<" and "<<k+1
	<<" of the "<<(left?"left":"right")<<" bath, when L="<<L<<endl;
    exit(1);
  }
  double alpha=left?B1n[k+1]:B2n[k+1];
  double gamma1=left?A1n[k]*.5:A2n[k]*.5;
  double gamma2=left?A1n[k+1]*.5:A2n[k+1]*.5;
  if(k+1==L-1) gamma2=2*gamma2;
  if(k==0) gamma1=2*gamma1; // the 0-th term goes here, unless only
			    // one mode!
  // Indices posi=getTwoBodyTermPositions(k,left);
  // int pos1=posi[0];
  // int pos2=posi[1];
  // int dimL=dims[min(pos1,pos2)];
  // int dimR=dims[max(pos1,pos2)];
  double gammaL=left?gamma2:gamma1;
  double gammaR=left?gamma1:gamma2;
  // cout<<"computeTwoBodyTermBath(k="<<k<<(left?",L)":",R)")<<" pos1="<<posi[0]
  //     <<" pos2="<<posi[1]<<", coeff NL="<<gammaL<<", coeff NR="<<gammaR
  //     <<", coeff inter="<<alpha<<endl;
  int d=2;
  mwArray sig0=identityMatrix(d);
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);

  // The two-body term is alpha/2 (sigma_x[k] sigma_x[k+1]+sigma_y[k]
  // sigma_y[k+1])+gamma1/2 (1+sigma_z[k])/2+gamma2/2 (1+sigma_z[k])/2
  // unless one of them is the edge of the chain, when gamma1 or 2 is
  // not divided by two
  constructOperatorProduct(result,sigX,sigX);
  mwArray tmp;
  constructOperatorProduct(tmp,sigY,sigY);
  result=-.5*alpha*(result+tmp);
  constructOperatorProduct(tmp,.5*(sig0+sigZ),sig0);
  result=result+gammaL*tmp;
  constructOperatorProduct(tmp,sig0,.5*(sig0+sigZ));
  result=result+gammaR*tmp;
}

void ThermofieldFermionHamiltonian::getTwoBodyTermBathExponential(mwArray& Ol,
								  mwArray& Or,
								  complex_t delta,
								  int k,
								  bool left) const{
  mwArray H12;
  Indices posi=getTwoBodyTermPositions(k,left);
  int pos1=posi[0];
  int pos2=posi[1];
  computeTwoBodyTermBath(H12,k,left);
   // cout<<"ThermofieldFermionHamiltonian::getTwoBodyTermExponential for "
   //     <<(left?"left":"right")<<" mode "<<k<<" corresponds to positions "
   //     <<pos1<<" and "<<pos2<<endl;
   // cout<<"The 2-body Hamiltonian is "<<H12<<endl;
  //H12.multiplyLeft(mwArray(delta));
  // Now take the matrix exponential
  mwArray expH;
  wrapper::expm(H12,expH,delta);
  // cout<<"Exponential of the h12="<<expH<<endl;
  // CAMBIO AQUI
  int nr=0;
  int d=2;
  split2term(Ol,Or,expH,d,d,nr);
}

void ThermofieldFermionHamiltonian::getCouplingExponential(mwArray& Opl,
							   mwArray& Ospin,mwArray& Opr,
							   complex_t delta) const{
  cout<<"ThermofieldFermionHamiltonian::getCouplingExponential"<<endl;
  // Construct the non-trivial operators for each of them 
  int d2=d*d;

  // Basic pieces
  mwArray sig0=identityMatrix(d);
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP); // sigma^+
  mwArray sigM=Hconjugate(sigP); //sigma^-
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataZ); //sigma_z

  // The term on the left (Bs=> with A1, B1): L^+ c_0^+ + H.c.
  mwArray tmp,termL,termR;
  constructOperatorProduct(tmp,sigZ,sigP);termL=tmp;
  constructOperatorProduct(tmp,sigP,sig0);termL=termL+tmp;
  constructOperatorProduct(tmp,sigM,termL);
  //termL=termL+tmp;
  termL=tmp+Hconjugate(tmp);
  cout<<"termL has size "<<termL.getDimensions()<<" and value "<<termL<<endl;
  constructOperatorProduct(tmp,sigM,sigZ);termR=tmp;
  constructOperatorProduct(tmp,sig0,sigM);termR=termR+tmp;
  constructOperatorProduct(tmp,termR,sigM);
  //termR=termR+tmp;
  termR=tmp+Hconjugate(tmp);
  cout<<"termR has size "<<termR.getDimensions()<<endl;
  mwArray H13;
  constructOperatorProduct(H13,-1.*B1n[0]*termL,sig0);
  cout<<"H13 has size "<<H13.getDimensions()<<endl;
  constructOperatorProduct(tmp,sig0,B2n[0]*termR);
  H13=H13+tmp;
  // This is the operator on three sites: now exp(delta*H13)
  wrapper::expm(H13,tmp,delta);
  // Now the SVD on two cuts!
  int Dl,Dr;
  split3term(Opl,Ospin,Opr,tmp,d,d2,d,Dl,Dr);
}

void ThermofieldFermionHamiltonian::split3term(mwArray& Ol,mwArray& Oc,mwArray& Or,mwArray& expH,
					       int dimL,int dimC,int dimR,int& Dl,int& Dr) const {
  cout<<"split3term of expH "<<expH.getDimensions()<<", dimL="<<dimL<<", dimC="<<dimC
      <<", dimR="<<dimR<<endl;
  // First the left cut with split2term
  mwArray tmp;
  split2term(Ol,tmp,expH,dimL,dimC*dimR,Dl);
  // Now tm is d2*d,Dl,d2*d
  tmp.reshape(Indices(dimC,dimR,Dl,dimC,dimR));
  tmp.permute(Indices(1,3,4,2,5));
  tmp.reshape(Indices(dimC*Dl*dimC,dimR*dimR));
  // And compute the SVD
  mwArray S; // place for the singular values
  Dr=0;
  double tol=0.;
  wrapper::svd(tmp,tol,Dr,Oc,S,Or);
  // cout<<"After SVD"<<endl;
  // redistribute the singular values
  S=sqrt(S);
  // TODO: check no NaN???
  Oc.multiplyRight(S);
  Or.multiplyLeft(S);
  // Now reshape adequately
  Oc.reshape(Indices(dimC,Dl,dimC,Dr));
  Or.reshape(Indices(Dr,dimR,dimR,1));
  Or.permute(Indices(2,1,3,4));  
}

void ThermofieldFermionHamiltonian::setSingleModeTerm(MPO& expH,int pos,int d,
						      complex_t delta,bool left) const{
  complex_t data[]={ONE_c,ZERO_c,ZERO_c,ZERO_c}; // (1+sigma_z)/2
  mwArray op(Indices(d,d),data);
  double coeff=left?A1n[0]:A2n[0];
  mwArray op2;
  wrapper::expm(op,op2,coeff*delta);
  op2.reshape(Indices(d,1,d,1));
  expH.setOp(pos,new Operator(op2),true);
}

void  ThermofieldFermionHamiltonian::initHMPO(){
  cout<<"ThermofieldFermionHamiltonian::initHMPO() Not yet implemented"<<endl;
}

void  ThermofieldFermionHamiltonian::initZ(){
  cout<<"ThermofieldFermionHamiltonian::initZ() Not yet implemented"<<endl;
}
