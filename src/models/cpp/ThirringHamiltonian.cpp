#include "ThirringHamiltonian.h"
#include "Indices.h"
#include "DoubleOperator.h"

using namespace std;
using namespace shrt;

ThirringHamiltonian::ThirringHamiltonian(int L_,double ma_,double g2_,
					 double lambda_,int Starget_,int d_,double gt2_,double mu_):
  L(L_),ma(ma_),g2(g2_),lambda(lambda_),Starget(Starget_),d(d_),gt2(gt2_),mu(mu_),hamil(L_){
  // Checks
  if(d!=2){
    cout<<"Error: ThirringHamiltonian only supports d=2"<<endl;
    exit(212);
  }
  if(lambda<0){
    cout<<"WARNING! penalty term lambda="<<lambda<<"<0!!"<<endl;
  }
  initZ();
  initHMPO();
}

ThirringHamiltonian::~ThirringHamiltonian(){};


void ThirringHamiltonian::initZ(){
  Z=mwArray(Indices(4,d,d));Z.fillWithZero();
  // The operators I need. NOTICE:SPIN OPERATORS (factor 1/2 included)
  mwArray sig0=identityMatrix(d);
  complex_t dataX[]={ZERO_c,.5*ONE_c,.5*ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  complex_t dataY[]={ZERO_c,.5*I_c,-.5*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(sigX.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Z.setElement(sigY.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      Z.setElement(sigZ.getElement(Indices(i1,i2)),Indices(3,i1,i2));
    }
  Z.reshape(Indices(4,d*d));
}

void ThirringHamiltonian::initHMPO(){
  int D=(lambda!=0.)?6:5;
  //  int signL=lambda>=0.?1:-1;
  for(int k=0;k<L;k++){
    int signK=k%2==0?1:-1; // (-1)^k
    int Dl=D,Dr=D;
    if(k==0) Dl=1;
    if(k==L-1) Dr=1;
    mwArray C(Indices(Dl,Dr,4));
    // Set elements
    if(k!=L-1){
      C.setElement(ONE_c,Indices(0,0,0)); // nothing yet
      // First of the pair
      C.setElement(-ONE_c,Indices(0,1,1));
      C.setElement(-ONE_c,Indices(0,2,2));
      C.setElement((g2+gt2*(1+signK))*ONE_c,Indices(0,3,3));
      // if penalty
      if(lambda!=0.){
	C.setElement(2*lambda*ONE_c,Indices(0,4,3));
	if(k!=0) C.setElement(ONE_c,Indices(4,4,0));
      }
    }
    if(k!=0){ 
      // After everything
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0));
      // Second of the pair
      C.setElement(ONE_c,Indices(1,Dr-1,1));
      C.setElement(ONE_c,Indices(2,Dr-1,2));
      C.setElement(ONE_c,Indices(3,Dr-1,3));      
      if(lambda!=0.){
	C.setElement(ONE_c,Indices(4,Dr-1,3));
      }
    }
    // The individual term
    double single=signK*ma+g2+mu+gt2-2.*lambda*Starget;
    if(k==0||k==L-1)
      single-=g2*.5;
    C.setElement(single*ONE_c,Indices(0,Dr-1,3));
    C.setElement(lambda*(Starget*Starget*(1./L)+.25)*ONE_c,Indices(0,Dr-1,0)); // just
								   // a constant!
    // Reshape, multiply operators and set in MPO
    C.reshape(Indices(Dl*Dr,4));
    C.multiplyRight(Z);
    C.reshape(Indices(Dl,Dr,d,d));
    C.permute(Indices(3,1,4,2));
    hamil.setOp(k,new Operator(C),true);
  }
}

void ThirringHamiltonian::getExponentialMPOeven(MPO& expHe,complex_t delta) const {
 getExponentialMPO(expHe,delta,true);
}

void ThirringHamiltonian::getExponentialMPOodd(MPO& expHo,complex_t delta) const {
 getExponentialMPO(expHo,delta,false);
}

void ThirringHamiltonian::getDoubleExponentialMPOeven(MPO& expHe,complex_t delta) const {
 getDoubleExponentialMPO(expHe,delta,true);
}
void ThirringHamiltonian::getDoubleExponentialMPOodd(MPO& expHo,complex_t delta) const {
 getDoubleExponentialMPO(expHo,delta,false);
}

void ThirringHamiltonian::getExtendedExponentialMPOeven(MPO& expHe,complex_t delta) const {
 getExtendedExponentialMPO(expHe,delta,true);
}
void ThirringHamiltonian::getExtendedExponentialMPOodd(MPO& expHo,complex_t delta) const {
 getExtendedExponentialMPO(expHo,delta,false);
}

void ThirringHamiltonian::getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,
						      int pos) const {
  mwArray H12;
  computeTwoBodyTerm(H12,pos);
  //  cout<<"****** POS "<<pos<<")"<<endl;
  // Now take the matrix exponential
  mwArray expH;
  wrapper::expm(H12,expH,delta);
  //putForMatlab(cout,expH,"expH12");
  //cout<<"Exponential of the h12="<<expH<<endl;

  // Obtained with indices (i'j')(ij) => permute to (i'i)(j'j) for svd
  expH.reshape(Indices(d,d,d,d));
  expH.permute(Indices(1,3,2,4));
  expH.reshape(Indices(d*d,d*d));
  // And compute the SVD
  mwArray S; // place for the singular values
  int nr=0;
  double tol=0.;
  //cout<<"Before decomposition, expH="<<expH<<endl;
  wrapper::svd(expH,tol,nr,Ol,S,Or);
  //cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<", S="<<S<<endl;
  // redistribute the singular values
  S=sqrt(S);
  // TODO: check no NaN???
  Ol.multiplyRight(S);
  Or.multiplyLeft(S);
  // reshape for operators (i'x 1 x i x alpha ))
  Ol.reshape(Indices(d,1,d,-1));
  Or.reshape(Indices(-1,d,d,1));
  Or.permute(Indices(2,1,3,4));  
  //cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<endl;
  // putForMatlab(cout,Ol,"Ol1");
  // putForMatlab(cout,Or,"Or1");
}


void ThirringHamiltonian::computeTwoBodyTerm(mwArray& result,int pos) const{
  if(pos>=L-1){
    cout<<"ERROR: ThirringHamiltonian::computeTwoBodyTerm for site "<<pos<<" when L="<<L<<endl;
    exit(212);
  }
  double fact1=pos==0?1.:.5;
  double fact2=pos==L-2?1.:.5;
  static mwArray XX,YY,ZZ,XI,IX,YI,IY,ZI,IZ;
  static bool firstCall(true);
  if(firstCall){
    mwArray Id=identityMatrix(d);
    complex_t dataX[]={ZERO_c,.5*ONE_c,.5*ONE_c,ZERO_c};
    mwArray SX(Indices(d,d),dataX);
    complex_t dataY[]={ZERO_c,.5*I_c,-.5*I_c,ZERO_c};
    mwArray SY(Indices(d,d),dataY);
    complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
    mwArray SZ(Indices(d,d),dataZ);
    constructOperatorProduct(XX,SX,SX);
    constructOperatorProduct(YY,SY,SY);
    constructOperatorProduct(ZZ,SZ,SZ);
    constructOperatorProduct(XI,SX,Id);
    constructOperatorProduct(IX,Id,SX);
    constructOperatorProduct(YI,SY,Id);
    constructOperatorProduct(IY,Id,SY);
    constructOperatorProduct(ZI,SZ,Id);
    constructOperatorProduct(IZ,Id,SZ);
    firstCall=false;
  }
  int signi=pos%2==0?1:-1;
  double Jz=g2+gt2*(1+signi);
  double hz1=signi*ma+gt2+.5*g2+mu;
  if(pos==0)hz1-=.5*g2;
  double hz2=-signi*ma+gt2+.5*g2+mu;
  if(pos==(L-2))hz2-=.5*g2;
  result=-1.*XX-1.*YY+Jz*ZZ+fact1*hz1*ZI+fact2*hz2*IZ;
}


void ThirringHamiltonian::getExponentialMPO(MPO& expH,complex_t delta,bool even) const {
  mwArray Ol,Or;
  expH.initLength(L);
  int kfirst=even?0:1;
  int klast=(L-1)-(L-kfirst)%2; // last occupied by the loop
  cout<<"getExponentialMPO even="<<even<<" kfirst="<<kfirst<<" klast="<<klast<<endl;
  // if(L%2==0){
  //   if(even) klast=L-2;
  //   else klast=L-3;
  // }
  // else{
  //   if(even) klast=L-3;
  //   else klast=L-2;
  // }
  for(int k=kfirst;k<klast;k+=2){
    getTwoBodyTermExponential(Ol,Or,delta,k);
    expH.setOp(k,new Operator(Ol),true);
    expH.setOp(k+1,new Operator(Or),true);
    //cout<<"ExponentialMPO set operators on "<<k<<" and "<<k+1<<endl;
  }
  // Fill in the edges with identity operators
  mwArray ident=identityMatrix(d);
  ident.reshape(Indices(d,1,d,1)); // just dx1xdx1
  int ptrId=-1; // where there is already an Id stored in the MPO
  if(kfirst==1){ // the first site is the identity
    expH.setOp(0,new Operator(ident),true);
    ptrId=0;
    //cout<<"Set pos 0 to Id"<<endl;
  }
  if(klast<L-1){ // sth to be filled at the end
    for(int k=klast+1;k<L;k++){
      if(ptrId>=0)
	expH.setOp(k,&expH.getOp(ptrId));
      else{
	expH.setOp(k,new Operator(ident),true);
	ptrId=k; // in case there were more than one, which cannot be!
      }
      //cout<<"Set pos "<<k<<" to Id"<<endl;
    }
  }
}

void ThirringHamiltonian::getDoubleExponentialMPO(MPO& expH,complex_t delta,bool even) const {
  mwArray Ol,Or,OlA,OrA;
  expH.initLength(L);
  int kfirst=even?0:1;
  int klast=(L-1)-(L-kfirst)%2; // last occupied by the loop
  // if(L%2==0){
  //   if(even) klast=L-2;
  //   else klast=L-3;
  // }
  // else{
  //   if(even) klast=L-3;
  //   else klast=L-2;
  // }
  for(int k=kfirst;k<klast;k+=2){
    getTwoBodyTermExponential(Ol,Or,delta,k);
    getTwoBodyTermExponential(OlA,OrA,conjugate(delta),k);
    // mwArray H12;
    // computeTwoBodyTerm(H12,k);
    // getTwoBodyTermExponential(Ol,Or,H12,delta);
    // getTwoBodyTermExponential(OlA,OrA,H12,conjugate(delta));
    DoubleOperator* operOl=new DoubleOperator(Ol,permute(OlA,Indices(3,2,1,4)));
    DoubleOperator* operOr=new DoubleOperator(Or,permute(OrA,Indices(3,2,1,4)));
    expH.setOp(k,operOl,true);
    expH.setOp(k+1,operOr,true);
    //cout<<"ExponentialMPO set operators on "<<k<<" and "<<k+1<<endl;
  }
  mwArray idPhys=identityMatrix(d);
  idPhys.reshape(Indices(d,1,d,1)); // for the idOp
  DoubleOperator* idOpPhys=new DoubleOperator(idPhys,idPhys);
  int ptrId=-1; // where there is already an Id stored in the MPO
  if(kfirst==1){ // the first site is the identity
    expH.setOp(0,idOpPhys,true);
    ptrId=0;
    //cout<<"Set pos 0 to Id"<<endl;
  }
  if(klast<L-1){ // sth to be filled at the end
    for(int k=klast+1;k<L;k++){
      if(ptrId>=0)
	expH.setOp(k,idOpPhys);
      else{
	expH.setOp(k,idOpPhys,true);
	ptrId=k; // in case there were more than one, which cannot be!
      }
      //cout<<"Set pos "<<k<<" to Id"<<endl;
    }
  }

}

void ThirringHamiltonian::getExtendedExponentialMPO(MPO& expH,complex_t delta,bool even) const {
  mwArray Ol,Or;
  expH.initLength(L);
  int kfirst=even?0:1;
  int klast=(L-1)-(L-kfirst)%2; // last occupied by the loop
  // if(L%2==0){
  //   if(even) klast=L-2;
  //   else klast=L-3;
  // }
  // else{
  //   if(even) klast=L-3;
  //   else klast=L-2;
  // }
  mwArray idPhys=identityMatrix(d);
  idPhys.reshape(Indices(d,1,d,1)); // for the idOp
  for(int k=kfirst;k<klast;k+=2){
    getTwoBodyTermExponential(Ol,Or,delta,k);
    // mwArray H12;
    // computeTwoBodyTerm(H12,k);
    // getTwoBodyTermExponential(Ol,Or,H12,delta);
    DoubleOperator* operOl=new DoubleOperator(Ol,idPhys);
    DoubleOperator* operOr=new DoubleOperator(Or,idPhys);
    expH.setOp(k,operOl,true);
    expH.setOp(k+1,operOr,true);
    //cout<<"ExponentialMPO set operators on "<<k<<" and "<<k+1<<endl;
  }
  DoubleOperator* idOpPhys=new DoubleOperator(idPhys,idPhys);
  int ptrId=-1; // where there is already an Id stored in the MPO
  if(kfirst==1){ // the first site is the identity
    expH.setOp(0,idOpPhys,true);
    ptrId=0;
    //cout<<"Set pos 0 to Id"<<endl;
  }
  if(klast<L-1){ // sth to be filled at the end
    for(int k=klast+1;k<L;k++){
      if(ptrId>=0)
	expH.setOp(k,idOpPhys);
      else{
	expH.setOp(k,idOpPhys,true);
	ptrId=k; // in case there were more than one, which cannot be!
      }
      //cout<<"Set pos "<<k<<" to Id"<<endl;
    }
  }

}

