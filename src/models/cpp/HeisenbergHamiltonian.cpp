#include "HeisenbergHamiltonian.h"
#include "Indices.h"
#include "DoubleOperator.h"

using namespace std;
using namespace shrt;

HeisenbergHamiltonian::HeisenbergHamiltonian(int L_,int d_):L(L_),d(d_),hamil(L_),offset(0.),labda(0.),Starget(0){initZ();}

HeisenbergHamiltonian::HeisenbergHamiltonian(int L_,vector<double> Jx,vector<double> Jy,
					     vector<double> Jz,vector<double> h,int d_,double offset_):
  L(L_),Jxi(Jx),Jyi(Jy),Jzi(Jz),hi(h),hxi(L_,0.),hyi(L_,0.),d(d_),hamil(L_),offset(offset_),lambda(0.),Starget(0){
  // Checks
  if(d!=2){
    cout<<"Error: HeisenbergHamiltonian only supports d=2"<<endl;
    exit(212);
  }
  if(Jx.size()!=L-1||Jy.size()!=L-1||Jz.size()!=L-1||h.size()!=L){
    cout<<"Error: length of coefficient vectors is not adequate for "
	<<"HeisenbergHamiltonian constructor"<<endl;
    exit(212);
  }
  initZ();
  initHMPO(offset);
}

HeisenbergHamiltonian::HeisenbergHamiltonian(int L_,double Jx,double Jy,double Jz,
					     vector<double> h,int d_,double offset_):
  L(L_),Jxi(L_-1,Jx),Jyi(L_-1,Jy),Jzi(L_-1,Jz),hi(h),hxi(L_,0.),hyi(L_,0.),d(d_),hamil(L_),lambda(0.),Staget(0),offset(offset_){
  // Checks
  if(d!=2){
    cout<<"Error: HeisenbergHamiltonian only supports d=2"<<endl;
    exit(212);
  }
  if(h.size()!=L){
    cout<<"Error: length of coefficient vector h is not adequate for "
	<<"HeisenbergHamiltonian constructor"<<endl;
    exit(212);
  }
  initZ();
  initHMPO(offset);
}

HeisenbergHamiltonian::HeisenbergHamiltonian(int L_,double Jx,double Jy,double Jz,
					     double h,int d_,double offset_):
  L(L_),Jxi(L_-1,Jx),Jyi(L_-1,Jy),Jzi(L_-1,Jz),hi(L_,h),hxi(L_,0.),hyi(L_,0.),d(d_),hamil(L_),lambda(0.),Staget(0),offset(offset_){
  // Checks
  if(d!=2){
    cout<<"Error: HeisenbergHamiltonian only supports d=2"<<endl;
    exit(212);
  }
  initZ();
  initHMPO(offset);
}

HeisenbergHamiltonian::HeisenbergHamiltonian(int L_,vector<double> Jx,vector<double> Jy,
					     vector<double> Jz,vector<double> hx,
					     vector<double> hy,vector<double> hz,int d_,double offset_):
  L(L_),Jxi(Jx),Jyi(Jy),Jzi(Jz),hxi(hx),hyi(hy),hi(hz),d(d_),hamil(L_),lambda(0.),Staget(0),offset(offset_){
  // Checks
  if(d!=2){
    cout<<"Error: HeisenbergHamiltonian only supports d=2"<<endl;
    exit(212);
  }
  if(Jx.size()!=L-1||Jy.size()!=L-1||Jz.size()!=L-1||hz.size()!=L||hx.size()!=L||hy.size()!=L){
    cout<<"Error: length of coefficient vectors is not adequate for "
	<<"HeisenbergHamiltonian constructor"<<endl;
    exit(212);
  }
  initZ();
  initHMPO(offset);
}

HeisenbergHamiltonian::HeisenbergHamiltonian(int L_,double Jx,double Jy,double Jz,
					     vector<double> hx,vector<double> hy,vector<double> hz,int d_,double offset_):
  L(L_),Jxi(L_-1,Jx),Jyi(L_-1,Jy),Jzi(L_-1,Jz),hxi(hx),hyi(hy),hi(hz),d(d_),hamil(L_),lambda(0.),Staget(0),offset(offset_){
  // Checks
  if(d!=2){
    cout<<"Error: HeisenbergHamiltonian only supports d=2"<<endl;
    exit(212);
  }
  if(hz.size()!=L||
     (hx.size()!=L&&hx.size()!=0)||
     (hy.size()!=L&&hy.size()!=0)){
    cout<<"Error: length of coefficient vector(s) h is not adequate for "
	<<"HeisenbergHamiltonian constructor"<<endl;
    exit(212);
  }
  initZ();
  initHMPO(offset);
}

HeisenbergHamiltonian::HeisenbergHamiltonian(int L_,double Jx,double Jy,double Jz,
					     double hx,double hy,double hz,int d_,double offset_):
  L(L_),Jxi(L_-1,Jx),Jyi(L_-1,Jy),Jzi(L_-1,Jz),hxi(L_,hx),hyi(L_,hy),hi(L_,hz),d(d_),hamil(L_),lambda(0.),Staget(0),offset(offset_){
  // Checks
  if(d!=2){
    cout<<"Error: HeisenbergHamiltonian only supports d=2"<<endl;
    exit(212);
  }
  initZ();
  initHMPO(offset);
}

void HeisenbergHamiltonian::setTargetSz(int Starget_,double penalty){
  Starget=Starget_;
  lambda=penalty;
  initHMPO(offset);
}


HeisenbergHamiltonian::~HeisenbergHamiltonian(){
  Jxi.clear();Jyi.clear();Jzi.clear();hi.clear();hxi.clear();hyi.clear();
}

void HeisenbergHamiltonian::getExponentialMPOeven(MPO& expHe,complex_t delta) const {
 getExponentialMPO(expHe,delta,true);
}
void HeisenbergHamiltonian::getExponentialMPOodd(MPO& expHo,complex_t delta) const {
 getExponentialMPO(expHo,delta,false);
}

void HeisenbergHamiltonian::getDoubleExponentialMPOeven(MPO& expHe,complex_t delta) const {
 getDoubleExponentialMPO(expHe,delta,true);
}
void HeisenbergHamiltonian::getDoubleExponentialMPOodd(MPO& expHo,complex_t delta) const {
 getDoubleExponentialMPO(expHo,delta,false);
}

void HeisenbergHamiltonian::getExtendedExponentialMPOeven(MPO& expHe,complex_t delta) const {
 getExtendedExponentialMPO(expHe,delta,true);
}
void HeisenbergHamiltonian::getExtendedExponentialMPOodd(MPO& expHo,complex_t delta) const {
 getExtendedExponentialMPO(expHo,delta,false);
}

void HeisenbergHamiltonian::getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,
						      int pos) const {
  mwArray H12;
  computeTwoBodyTerm(H12,pos);
   cout<<"****** POS "<<pos<<")"<<endl;
  getTwoBodyTermExponential(Ol,Or,H12,delta);
  // putForMatlab(cout,H12,"H12");
  // // Now take the matrix exponential
  // mwArray expH;
  // wrapper::expm(H12,expH,delta);
  // //putForMatlab(cout,expH,"expH12");
  // //cout<<"Exponential of the h12="<<expH<<endl;

  // // Obtained with indices (i'j')(ij) => permute to (i'i)(j'j) for svd
  // expH.reshape(Indices(d,d,d,d));
  // expH.permute(Indices(1,3,2,4));
  // expH.reshape(Indices(d*d,d*d));
  // // And compute the SVD
  // mwArray S; // place for the singular values
  // int nr=0;
  // double tol=0.;
  // //cout<<"Before decomposition, expH="<<expH<<endl;
  // wrapper::svd(expH,tol,nr,Ol,S,Or);
  // //cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<", S="<<S<<endl;
  // // redistribute the singular values
  // S=sqrt(S);
  // // TODO: check no NaN???
  // Ol.multiplyRight(S);
  // Or.multiplyLeft(S);
  // // reshape for operators (i'x 1 x i x alpha ))
  // Ol.reshape(Indices(d,1,d,-1));
  // Or.reshape(Indices(-1,d,d,1));
  // Or.permute(Indices(2,1,3,4));  
  // //cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<endl;
  // //putForMatlab(cout,Ol,"Ol");
  // //putForMatlab(cout,Or,"Or");
}

void HeisenbergHamiltonian::getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,
						      const mwArray& H12,complex_t delta) const {
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

void HeisenbergHamiltonian::computeTwoBodyTerm(mwArray& result,int pos) const{
  if(pos>=L-1){
    cout<<"ERROR: HeisenbergHamiltonian::computeTwoBodyTerm for site "<<pos<<" when L="<<L<<endl;
    exit(212);
  }
  double h1=pos==0?1.:.5;
  double h2=pos==L-2?1.:.5;
  computeTwoBodyTerm(result,Jxi[pos],Jyi[pos],Jzi[pos],h1*hxi[pos],h2*hxi[pos+1],h1*hyi[pos],h2*hyi[pos+1],h1*hi[pos],h2*hi[pos+1]);
}

void HeisenbergHamiltonian::computeTwoBodyTerm(mwArray& result,double Jx,double Jy,double Jz,
					       double h1,double h2) const{
  static mwArray XX,YY,ZZ,ZI,IZ;
  static bool firstCall(true);
  if(firstCall){
    complex_t dataXX[]={ZERO_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ZERO_c};
    XX=.25*mwArray(Indices(d*d,d*d),dataXX);
    complex_t dataYY[]={ZERO_c,ZERO_c,ZERO_c,-ONE_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,-ONE_c,ZERO_c,ZERO_c,ZERO_c};
    YY=.25*mwArray(Indices(d*d,d*d),dataYY);
    complex_t dataZZ[]={ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,-ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,-ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,ONE_c};
    ZZ=.25*mwArray(Indices(d*d,d*d),dataZZ);
    complex_t dataIZ[]={ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,-ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,-ONE_c};
    IZ=.5*mwArray(Indices(d*d,d*d),dataIZ);
    complex_t dataZI[]={ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,-ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,-ONE_c};
    ZI=.5*mwArray(Indices(d*d,d*d),dataZI);
    firstCall=false;
  }
  result=Jx*XX+Jy*YY+Jz*ZZ+h1*ZI+h2*IZ;
}

void HeisenbergHamiltonian::computeTwoBodyTerm(mwArray& result,double Jx,double Jy,double Jz,
					       double hx1,double hx2,double hy1,double hy2,double hz1,double hz2) const{
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
  result=Jx*XX+Jy*YY+Jz*ZZ+hx1*XI+hx2*IX+hy1*YI+hy2*IY+hz1*ZI+hz2*IZ;
}



void HeisenbergHamiltonian::initZ(){
  Z=mwArray(Indices(4,d,d));Z.fillWithZero();
  // The operators I need
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

void HeisenbergHamiltonian::initHMPO(double offset){
  int D=(lambda!=0)?6:5;
  for(int k=0;k<L;k++){
    int Dl=D,Dr=D;
    if(k==0) Dl=1;
    if(k==L-1) Dr=1;
    mwArray C(Indices(Dl,Dr,4));
    // Set elements
    if(k!=L-1){
      C.setElement(ONE_c,Indices(0,0,0)); // nothing yet
      // First of the pair
      C.setElement(Jxi[k]*ONE_c,Indices(0,1,1));
      C.setElement(Jyi[k]*ONE_c,Indices(0,2,2));
      C.setElement(Jzi[k]*ONE_c,Indices(0,3,3));
      // if penalty
      if(lambda!=0.){
	C.setElement(2*lambda*ONE_c,Indices(0,4,3));//first sigz of penalty
	if(k!=0) C.setElement(ONE_c,Indices(4,4,0)); //in between
      }
    }
    if(k!=0){ 
      // After everything
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0));
      // Second of the pair
      C.setElement(ONE_c,Indices(1,Dr-1,1));
      C.setElement(ONE_c,Indices(2,Dr-1,2));
      C.setElement(ONE_c,Indices(3,Dr-1,3));
      if(lambda!=0.){ // second of the penalty
	C.setElement(ONE_c,Indices(4,Dr-1,3));
      }
    }
    // The individual term(s)
    C.setElement((hi[k]-2*lambda*Starget)*ONE_c,Indices(0,Dr-1,3));
    if(hxi[k]!=0) C.setElement(hxi[k]*ONE_c,Indices(0,Dr-1,1));
    if(hyi[k]!=0) C.setElement(hyi[k]*ONE_c,Indices(0,Dr-1,2));
    // The offset term, if present
    if(offset!=0)
      C.setElement((offset/L)*ONE_c,Indices(0,Dr-1,0));
    // Reshape, multiply operators and set in MPO
    C.reshape(Indices(Dl*Dr,4));
    C.multiplyRight(Z);
    C.reshape(Indices(Dl,Dr,d,d));
    C.permute(Indices(3,1,4,2));
    hamil.setOp(k,new Operator(C),true);
  }
}

void HeisenbergHamiltonian::getExponentialMPO(MPO& expH,complex_t delta,bool even) const {
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

void HeisenbergHamiltonian::getDoubleExponentialMPO(MPO& expH,complex_t delta,bool even) const {
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
    mwArray H12;
    computeTwoBodyTerm(H12,k);
    getTwoBodyTermExponential(Ol,Or,H12,delta);
    getTwoBodyTermExponential(OlA,OrA,H12,conjugate(delta));
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

void HeisenbergHamiltonian::getExtendedExponentialMPO(MPO& expH,complex_t delta,bool even) const {
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
    mwArray H12;
    computeTwoBodyTerm(H12,k);
    getTwoBodyTermExponential(Ol,Or,H12,delta);
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

void HeisenbergHamiltonian::getCommutatorMPO(MPO& mpo){
  int D=8; // bond dimension of the MPO
  mwArray Z; // The operators
  initZdoub(Z);
  int nrOps=Z.getDimension(0); // Nr of operators (7)

  // Since the coefficients may depend on the site, I need a non-TI MPO
  for(int k=0;k<L;k++){
    int Dl=D,Dr=D;
    if(k==0) Dl=1;
    if(k==L-1) Dr=1;
    mwArray C(Indices(Dl,Dr,4));
    // Set elements
    if(k!=L-1){
      C.setElement(ONE_c,Indices(0,0,0)); // nothing yet
      // First of the pair (system-system) or (ancilla-ancilla)
      C.setElement(Jxi[k]*ONE_c,Indices(0,1,1));
      C.setElement(-Jxi[k]*ONE_c,Indices(0,2,2));
      C.setElement(Jyi[k]*ONE_c,Indices(0,3,3));
      C.setElement(-Jyi[k]*ONE_c,Indices(0,4,4));
      C.setElement(Jzi[k]*ONE_c,Indices(0,5,5));
      C.setElement(-Jzi[k]*ONE_c,Indices(0,6,6));
    }
    if(k!=0){ 
      // After everything
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0));
      // Second of the pair
      C.setElement(ONE_c,Indices(1,Dr-1,1));
      C.setElement(ONE_c,Indices(2,Dr-1,2));
      C.setElement(ONE_c,Indices(3,Dr-1,3));      
      C.setElement(ONE_c,Indices(4,Dr-1,4));      
      C.setElement(ONE_c,Indices(5,Dr-1,5));      
      C.setElement(ONE_c,Indices(6,Dr-1,6));      
    }
    // The individual term(s)
    C.setElement(hi[k]*ONE_c,Indices(0,Dr-1,5));
    C.setElement(-hi[k]*ONE_c,Indices(0,Dr-1,6));

    // Reshape, multiply operators and set in MPO
    C.reshape(Indices(Dl*Dr,nrOps));
    C.multiplyRight(Z);
    C.reshape(Indices(Dl,Dr,d,d));
    C.permute(Indices(3,1,4,2));
    mpo.setOp(k,new Operator(C),true);
  }
}


void HeisenbergHamiltonian::initZdoub(mwArray& Z){
  // Basic pieces of the commutator MPO
  mwArray sig0=identityMatrix(d);//identity
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),datax);//sigmax
  complex_t datay[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigY(Indices(d,d),datax);//sigmax
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz

  // Because my sites are double, the actual individual operators will
  // be tensor products of these, or of these with the identity So, I
  // construct the basic products one by one.  (Notice that not all of
  // them are independent, and I could also do a more efficient
  // contraction, but this cost is negligible compared to anything
  // else, as I only need to do this onec when initializing the
  // problem, so I opt for the clearerst construction here).

  // *** First of all, we need the identity on both
  mwArray term0=identityMatrix(d*d);
  // *** Now the ones acting on system and ancilla
  // (1) sigX (x) Id
  mwArray term1;
  constructOperatorProduct(term1,sigX,sig0);
  // (2)  Id (x) sigX -> just a permutation of the previous one
  mwArray term2;
  constructOperatorProduct(term2,sig0,sigX);
  // (3) sigY (x) Id
  mwArray term3;
  constructOperatorProduct(term3,sigY,sig0);
  // (4)  Id (x) sigY -> just a permutation of the previous one
  mwArray term4;
  constructOperatorProduct(term4,sig0,sigY);
  // (5) sigZ (x) Id
  mwArray term5; 
  constructOperatorProduct(term5,sigZ,sig0);
  // (6)  Id (x) sigZ -> just a permutation of the previous one
  mwArray term6;
  constructOperatorProduct(term6,sig0,sigZ);

  // Fill in the operators in the Z array, to be contracted with C in
  // the construction of the MPO. The order is as listed above.
  int nrOps=5;
  Z=mwArray(Indices(nrOps,d*d,d*d));
  for(int d1=0;d1<d*d;d1++)
    for(int d2=0;d2<d*d;d2++){
      Z.setElement(term0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
      Z.setElement(term1.getElement(Indices(d1,d2)),Indices(1,d1,d2));
      Z.setElement(term2.getElement(Indices(d1,d2)),Indices(2,d1,d2));
      Z.setElement(term3.getElement(Indices(d1,d2)),Indices(3,d1,d2));
      Z.setElement(term4.getElement(Indices(d1,d2)),Indices(4,d1,d2));
      Z.setElement(term5.getElement(Indices(d1,d2)),Indices(5,d1,d2));
      Z.setElement(term6.getElement(Indices(d1,d2)),Indices(6,d1,d2));
    }
}

