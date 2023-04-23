#include "ThermofieldHamiltonian.h"
#include "Indices.h"

using namespace std;
using namespace shrt;

//#define TRUNCESP 1

ThermofieldHamiltonian::ThermofieldHamiltonian(int L_,vector<int> dims_,
					       const vector<double>& A1n_,
					       const vector<double>& A2n_,
					       const vector<double>& B1n_,
					       const vector<double>& B2n_,
					       mwArray Hfree_,mwArray OL_):
  L(L_),hamil(2*L_+1),dims(dims_),B2n(B2n_),B1n(B1n_),
  A2n(A2n_),A1n(A1n_),Hfree(Hfree_),OL(OL_){
   // cout<<"Creating TFHamiltonian with arguments:"<<endl;
   // cout<<"dims="<<dims_<<endl;
   // cout<<"A1n="<<A1n<<endl;
   // cout<<"A2n="<<A2n<<endl;
   // cout<<"B1n="<<B1n<<endl;
   // cout<<"B2n="<<B2n<<endl;
   // cout<<"Hfree="<<Hfree<<endl;
   // cout<<"OL="<<OL<<endl;
  // Introduce the signs
  for(int k=0;k<L;k++){
    if(k>0)
      B2n[k]=-B2n[k];
    A2n[k]=-A2n[k];
  }
  initZ();
  initHMPO();
}

ThermofieldHamiltonian::ThermofieldHamiltonian(int L_,vector<int> dims_,
					       const vector<double>& A1n_,
					       const vector<double>& A2n_,
					       const vector<double>& B1n_,
					       const vector<double>& B2n_):L(L_),hamil(2*L_+1),dims(dims_),B2n(B2n_),B1n(B1n_),
									   A2n(A2n_),A1n(A1n_){};

ThermofieldHamiltonian::~ThermofieldHamiltonian(){};

Indices ThermofieldHamiltonian::getTwoBodyTermPositions(int k,bool left) const{
  int pos1=left?L-k-1:L+k+1; // real pos in the chain
  int pos2=left?pos1-1:pos1+1;
  return Indices(pos1,pos2);
}

void ThermofieldHamiltonian::computeTwoBodyTermBath(mwArray& result,int k,
						    bool left) const{
  if(k>=L-1){
    cout<<"Error:cannot compute two body bosonic term on the "
	<<(left?"left":"right")<<" for modes "<<k<<" and "<<k+1
	<<", when L="<<L<<endl;
    exit(1);
  }
  double alpha=left?B1n[k+1]:B2n[k+1];
  double gamma1=left?A1n[k]*.5:A2n[k]*.5;
  double gamma2=left?A1n[k+1]*.5:A2n[k+1]*.5;
  if(k+1==L-1) gamma2=2*gamma2;
  if(k==0) gamma1=2*gamma1; // the 0-th term goes here, unless only
			    // one mode!
  Indices posi=getTwoBodyTermPositions(k,left);
  int pos1=posi[0];
  int pos2=posi[1];
  //int dim1=dims[pos1];
  //int dim2=dims[pos2];
  int dimL=dims[min(pos1,pos2)];
  int dimR=dims[max(pos1,pos2)];
  double gammaL=left?gamma2:gamma1;
  double gammaR=left?gamma1:gamma2;
  // cout<<"computeTwoBodyTermBath(k="<<k<<(left?",L)":",R)")<<" pos1="<<pos1
  //      <<" pos2="<<pos2<<", coeff NL="<<gammaL<<", coeff NR="<<gammaR
  //      <<", coeff inter="<<alpha<<endl;

  mwArray aCreatL(Indices(dimL,dimL));
  for(int d=0;d<dimL-1;d++){
    aCreatL.setElement(sqrt(d+1),0.,Indices(d+1,d));
  }
  mwArray aDestrL(aCreatL);aDestrL.Hconjugate();
  mwArray aCreatR(Indices(dimR,dimR));
  for(int d=0;d<dimR-1;d++){
    aCreatR.setElement(sqrt(d+1),0.,Indices(d+1,d));
  }
  mwArray aDestrR(aCreatR);aDestrR.Hconjugate();
  mwArray NopL=aCreatL*aDestrL;
  mwArray NopR=aCreatR*aDestrR;
  mwArray IdL=identityMatrix(dimL);
  mwArray IdR=identityMatrix(dimR);
#ifdef TRUNCESP
  aCreatL.setElement(sqrt(dimL),0.,Indices(dimL-1,dimL-1));
  aDestrL.setElement(sqrt(dimL),0.,Indices(dimL-1,dimL-1));
  aCreatR.setElement(sqrt(dimR),0.,Indices(dimR-1,dimR-1));
  aDestrR.setElement(sqrt(dimR),0.,Indices(dimR-1,dimR-1));
#endif
  // The two-body term is alpha(ck+ c(k+1)+h.c.)+gamma1/2 N1 +gamma2/2 N2
  // unless one of them is the edge of the chain, when gamma1 or 2 is
  // not divided by two
  aCreatL.reshape(Indices(dimL*dimL,1));
  aDestrR.reshape(Indices(1,dimR*dimR));
  result=aCreatL;
  result.multiplyRight(aDestrR);
  result.reshape(Indices(dimL,dimL,dimR,dimR));
  result.permute(Indices(1,3,2,4));
  result.reshape(Indices(dimL*dimR,dimL*dimR));
  NopL.reshape(Indices(dimL*dimL,1));
  IdL.reshape(Indices(dimL*dimL,1));
  NopR.reshape(Indices(1,dimR*dimR));
  IdR.reshape(Indices(1,dimR*dimR));
  NopL.multiplyRight(IdR);
  NopR.multiplyLeft(IdL);
  mwArray diag=gammaL*NopL+gammaR*NopR;
  diag.reshape(Indices(dimL,dimL,dimR,dimR));
  diag.permute(Indices(1,3,2,4));
  diag.reshape(Indices(dimL*dimR,dimL*dimR));
  result=result+Hconjugate(result);
  result=alpha*result+diag;
  // cout<<"Boson term nr "<<k<<" (pos "<<min(pos1,pos2)<<")"<<endl;
  // // debugging
  // ofstream* debOut=new ofstream("debugQF.m",ios::app);
  // if(left) putForMatlab(*debOut,result,"term"); else putForMatlab(*debOut,result,"termR");
  // debOut->close();delete debOut;
}

void ThermofieldHamiltonian::getTwoBodyTermBathExponential(mwArray& Ol,
							     mwArray& Or,
							     complex_t delta,
							     int k,
							     bool left) const{
  mwArray H12;
  Indices posi=getTwoBodyTermPositions(k,left);
  int pos1=posi[0];
  int pos2=posi[1];
  //  int dimL=left?dim2:dim1;
  // int dimR=left?dim1:dim2;
  int dimL=dims[min(pos1,pos2)];
  int dimR=dims[max(pos1,pos2)];
  computeTwoBodyTermBath(H12,k,left);
   // cout<<"ThermofieldHamiltonian::getTwoBodyTermBathExponential for "
   //     <<(left?"left":"right")<<" mode "<<k<<" corresponds to positions "
   //     <<pos1<<" and "<<pos2<<", with dimensions "<<dimL<<" and "<<dimR<<endl;
   // cout<<"The 2-body Hamiltonian is "<<H12<<endl;
  //H12.multiplyLeft(mwArray(delta));
  // Now take the matrix exponential
  mwArray expH;
  wrapper::expm(H12,expH,delta);
  // cout<<"Exponential of the h12="<<expH<<endl;
  // CAMBIO AQUI
  int nr=0;
  split2term(Ol,Or,expH,dimL,dimR,nr);
}

void ThermofieldHamiltonian::split2term(mwArray& Ol,mwArray& Or,mwArray& expH,int dimL,int dimR,int& nr) const {
  //  cout<<"ThermofieldHamiltonian::split2term expH="<<expH<<", dimL="<<dimL<<", dimR="<<dimR<<endl;
  expH.reshape(Indices(dimL,dimR,dimL,dimR));
  expH.permute(Indices(1,3,2,4));
  expH.reshape(Indices(dimL*dimL,dimR*dimR));
  // And compute the SVD
  mwArray S; // place for the singular values
  nr=0;
  double tol=0.;
  //  cout<<"Before SVD expH="<<expH<<endl;
  wrapper::svd(expH,tol,nr,Ol,S,Or);
  // cout<<"After SVD"<<endl;
  // redistribute the singular values
  S=sqrt(S);
  // TODO: check no NaN???
  Ol.multiplyRight(S);
  Or.multiplyLeft(S);
  // reshape for operators (i'x 1 x i x alpha ))
  Ol.reshape(Indices(dimL,1,dimL,-1));
  Or.reshape(Indices(-1,dimR,dimR,1));
  Or.permute(Indices(2,1,3,4));  
}

void ThermofieldHamiltonian::getCouplingExponential(mwArray& Opl,
						    mwArray& Ospin,mwArray& Opr,
						    complex_t delta) const{
  //  cout<<"ThermofieldHamiltonian::getCouplingExponential"<<endl;
  // All Hamiltonian term involving the system and its neighbouring bosons, and the exponential, split in 3-site MPO
  int dimL=dims[L-1];
  int dimR=dims[L+1];
  int d=dims[L];
  // exp(delta B)-1

  mwArray bosonL(Indices(dimL,dimL));
  for(int id=0;id<dimL-1;id++){
   bosonL.setElement(sqrt(id+1),0.,Indices(id+1,id));
  }
  bosonL.Hconjugate(); // the destruct. is the one with OL
  
  mwArray bosonR(Indices(dimR,dimR));
  for(int id=0;id<dimR-1;id++){
    bosonR.setElement(sqrt(id+1),0.,Indices(id+1,id));
  }
  bosonR.Hconjugate(); // the destruct. is the one with OL
  
  // term is prop to OL* boson + h.c. on each side
  mwArray termL,termR;
  constructOperatorProduct(termL,bosonL,OL);termL=termL+Hconjugate(termL);
  constructOperatorProduct(termR,OL,bosonR);termR=termR+Hconjugate(termR);
  mwArray H13;
  constructOperatorProduct(H13,identityMatrix(dimL),B2n[0]*termR);
  mwArray tmp;
  constructOperatorProduct(tmp,B1n[0]*termL,identityMatrix(dimR));
  H13=H13+tmp;

  constructOperatorProduct(termL,identityMatrix(dimL),Hfree);
  constructOperatorProduct(tmp,termL,identityMatrix(dimR));
  H13=H13+tmp;
  // ofstream* debOut=new ofstream("debugQF.m",ios::app);
  // *debOut<<"dtExp="<<delta<<endl;
  // putForMatlab(*debOut,H13,"H13");
  // debOut->close();delete debOut;
  wrapper::expm(H13,tmp,delta);
    // Now the SVD on two cuts!
  int nl,nr;
  split3term(Opl,Ospin,Opr,tmp,dimL,d,dimR,nl,nr);
  // debOut=new ofstream("debugQF.m",ios::app);
  // putForMatlab(*debOut,Opl,"Opl");
  // putForMatlab(*debOut,Ospin,"Ospin");
  // putForMatlab(*debOut,Opr,"Opr");
  // debOut->close();delete debOut;
  }


void ThermofieldHamiltonian::split3term(mwArray& Ol,mwArray& Oc,mwArray& Or,mwArray& expH,
					int dimL,int dimC,int dimR,int& nl,int& nr) const {
  // cout<<"split3term of expH "<<expH.getDimensions()<<", dimL="<<dimL<<", dimC="<<dimC
  //     <<", dimR="<<dimR<<endl;
  // First the left cut with split2term
  mwArray tmp;
  split2term(Ol,tmp,expH,dimL,dimC*dimR,nl);
  // Now tm is d2*d,Dl,d2*d
  tmp.reshape(Indices(dimC,dimR,nl,dimC,dimR));
  tmp.permute(Indices(1,3,4,2,5));
  tmp.reshape(Indices(dimC*nl*dimC,dimR*dimR));
  // And compute the SVD
  mwArray S; // place for the singular values
  nr=0;
  double tol=0.;
  wrapper::svd(tmp,tol,nr,Oc,S,Or);
  // cout<<"After SVD"<<endl;
  // redistribute the singular values
  S=sqrt(S);
  // TODO: check no NaN???
  Oc.multiplyRight(S);
  Or.multiplyLeft(S);
  // Now reshape adequately
  Oc.reshape(Indices(dimC,nl,dimC,nr));
  Or.reshape(Indices(nr,dimR,dimR,1));
  Or.permute(Indices(2,1,3,4));  
}



void ThermofieldHamiltonian::getFreeSpinExponential(mwArray& Op,
						    complex_t delta) const{
  //  int d=dims[L];
  wrapper::expm(Hfree,Op,delta);
}

void ThermofieldHamiltonian::setSingleModeTerm(MPO& expHe,int pos,int dim,complex_t delta,bool left) const{
  mwArray op1(Indices(dim,dim));
  for(int d=0;d<dim-1;d++){
    op1.setElement(sqrt(d+1),0.,Indices(d+1,d));
  }
  mwArray op2(op1);op2.Hconjugate();
  op1.multiplyRight(op2);
  double coeff=left?A1n[0]:A2n[0];
  wrapper::expm(op1,op2,coeff*delta);
  op2.reshape(Indices(dim,1,dim,1));
  expHe.setOp(pos,new Operator(op2),true);
}

void ThermofieldHamiltonian::getExponentialMPOeven(MPO& expHe,complex_t delta) const{
  //  cout<<"ThermofieldHamiltonian::getExponentialMPOeven"<<endl;
  expHe.initLength(2*L+1);
  int d=dims[L];
  mwArray tmpOp=identityMatrix(d); // in this implementation, the free system Hamiltonian is including with the coupling
  //  getFreeSpinExponential(tmpOp,delta); 
  //int d=tmpOp.getDimension(0);
  tmpOp.reshape(Indices(d,1,d,1));
  expHe.setOp(L,new Operator(tmpOp),true);
  //  cout<<"Set operator at "<<L<<" to "<<tmpOp.getDimensions()<<endl;
  // Now set bosonic operators for 01,23...
  for(int k=0;k<L-1;k+=2){
    mwArray Ol,Or;
    // Left
    getTwoBodyTermBathExponential(Ol,Or,delta,k,true);
    Indices posi=getTwoBodyTermPositions(k,true);
    int posL=min(posi[0],posi[1]);
    int posR=max(posi[0],posi[1]);
    expHe.setOp(posL,new Operator(Ol),true);
    expHe.setOp(posR,new Operator(Or),true);
    //    cout<<"Set operator at "<<posL<<" to "<<Ol.getDimensions()<<endl;
    //    cout<<"Set operator at "<<posR<<" to "<<Or.getDimensions()<<endl;
    // Right
    getTwoBodyTermBathExponential(Ol,Or,delta,k,false);
    posi=getTwoBodyTermPositions(k,false);
    posL=min(posi[0],posi[1]);
    posR=max(posi[0],posi[1]);
    expHe.setOp(posL,new Operator(Ol),true);
    expHe.setOp(posR,new Operator(Or),true);
    //    cout<<"Set operator at "<<posL<<" to "<<Ol.getDimensions()<<endl;
    //    cout<<"Set operator at "<<posR<<" to "<<Or.getDimensions()<<endl;
  }
  // Now it can be that there is an identity at each edge!
  if(expHe.isEmpty(0)){
    int dim0=dims[0];
    // If there was only one mode, I also need the diagonal term here!
    if(L==1){
      setSingleModeTerm(expHe,0,dim0,delta,true); 
      // mwArray op1(Indices(dim0,dim0));
      // for(int d=0;d<dim0-1;d++){
      // 	op1.setElement(sqrt(d+1),0.,Indices(d+1,d));
      // }
      // mwArray op2(op1);op2.Hconjugate();
      // op1.multiplyRight(op2);
      // wrapper::expm(op1,op2,A1n[0]*delta);
      // op2.reshape(Indices(dim0,1,dim0,1));
      // expHe.setOp(0,new Operator(op2),true);
      // //cout<<"Set operator at "<<0<<" to diagonal term for mode 0 "<<op2.getDimensions()<<endl;
    }
    else{
      mwArray id0=identityMatrix(dim0);
      id0.reshape(Indices(dim0,1,dim0,1));
      expHe.setOp(0,new Operator(id0),true);
      cout<<"Set operator at "<<0<<" to Id "<<id0.getDimensions()<<endl;
    }
  }
  if(expHe.isEmpty(2*L)){
    int dimL=dims[2*L];
    if(L==1){
      setSingleModeTerm(expHe,2*L,dimL,delta,false); 
      // mwArray op1(Indices(dimL,dimL));
      // for(int d=0;d<dimL-1;d++){
      // 	op1.setElement(sqrt(d+1),0.,Indices(d+1,d));
      // }
      // mwArray op2(op1);op2.Hconjugate();
      // op1.multiplyRight(op2);
      // wrapper::expm(op1,op2,A2n[0]*delta);
      // op2.reshape(Indices(dimL,1,dimL,1));
      // expHe.setOp(2*L,new Operator(op2),true);
      // //      cout<<"Set operator at "<<2*L<<" to diagonal term for mode 0 "<<op2.getDimensions()<<endl;
    }
    else{    
      mwArray idL=identityMatrix(dimL);
      idL.reshape(Indices(dimL,1,dimL,1));
      expHe.setOp(2*L,new Operator(idL),true);
      cout<<"Set operator at "<<2*L<<" to Id "<<idL.getDimensions()<<endl;
    }
  }
}

void ThermofieldHamiltonian::getExponentialMPOodd(MPO& expHo,complex_t delta) const{
  expHo.initLength(2*L+1);
  // The term in the center, from the couplings of the spin and
  // zero-modes
  //  cout<<"ThermofieldHamiltonian::getExponentialMPOodd"<<endl;
  mwArray Ol,Ospin,Or;
  getCouplingExponential(Ol,Ospin,Or,delta);
  expHo.setOp(L-1,new Operator(Ol),true);
  expHo.setOp(L,new Operator(Ospin),true);
  expHo.setOp(L+1,new Operator(Or),true);
  // cout<<"Set operator at "<<L-1<<" to "<<Ol.getDimensions()<<endl;
  // cout<<"Set operator at "<<L<<" to "<<Ospin.getDimensions()<<endl;
  // cout<<"Set operator at "<<L+1<<" to "<<Or.getDimensions()<<endl;
  // Now set bosonic operators for 12,34...
  for(int k=1;k<L-1;k+=2){
    // Left
    getTwoBodyTermBathExponential(Ol,Or,delta,k,true);
    Indices posi=getTwoBodyTermPositions(k,true);
    int posL=min(posi[0],posi[1]);
    int posR=max(posi[0],posi[1]);
    expHo.setOp(posL,new Operator(Ol),true);
    expHo.setOp(posR,new Operator(Or),true);
    //cout<<"Set operator at "<<posL<<" to "<<Ol.getDimensions()<<endl;
    //cout<<"Set operator at "<<posR<<" to "<<Or.getDimensions()<<endl;
    // Right
    getTwoBodyTermBathExponential(Ol,Or,delta,k,false);
    posi=getTwoBodyTermPositions(k,false);
    posL=min(posi[0],posi[1]);
    posR=max(posi[0],posi[1]);
    expHo.setOp(posL,new Operator(Ol),true);
    expHo.setOp(posR,new Operator(Or),true);
    //cout<<"Set operator at "<<posL<<" to "<<Ol.getDimensions()<<endl;
    //cout<<"Set operator at "<<posR<<" to "<<Or.getDimensions()<<endl;
  }
  // Now it can be that there is an identity at each edge!
  if(expHo.isEmpty(0)){
    int dim0=dims[0];
    mwArray id0=identityMatrix(dim0);
    id0.reshape(Indices(dim0,1,dim0,1));
    expHo.setOp(0,new Operator(id0),true);
    cout<<"Set operator at "<<0<<" to Id "<<id0.getDimensions()<<endl;
  }
  if(expHo.isEmpty(2*L)){
    int dimL=dims[2*L];
    mwArray idL=identityMatrix(dimL);
    idL.reshape(Indices(dimL,1,dimL,1));
    expHo.setOp(2*L,new Operator(idL),true);
    cout<<"Set operator at "<<2*L<<" to Id "<<idL.getDimensions()<<endl;
  }
}

void  ThermofieldHamiltonian::initHMPO(){
  cout<<"ThermofieldHamiltonian::initHMPO() Not yet implemented"<<endl;
}

void  ThermofieldHamiltonian::initZ(){
  cout<<"ThermofieldHamiltonian::initZ() Not yet implemented"<<endl;
}
