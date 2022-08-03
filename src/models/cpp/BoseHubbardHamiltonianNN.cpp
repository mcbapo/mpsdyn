#include "BoseHubbardHamiltonian.h"
#include "Indices.h"

using namespace std;
using namespace shrt;

BoseHubbardHamiltonian::BoseHubbardHamiltonian(int L_,int N_,double t_,
					       double U_,double mu_,double V_,
					       double Vint_):
  hamil(L_),tis(1,t_),Uis(1,.5*U_),muis(1,mu_),Vis(1,V_),Vintis(1,Vint_){
  L=L_;N=N_;
  //t=t_;U=U_;mu=mu_;V=V_;
  // Construct the MPO
  initZ();
  initHMPO();
}

BoseHubbardHamiltonian::BoseHubbardHamiltonian(int L_,int N_,double* tis_,
					       double* Uis_,double* muis_,
					       double* Vis_,double* Vintis_):
  hamil(L_){
  L=L_;N=N_;
  //t=U=mu=V=0; // no constant values
  for(int k=0;k<L;k++){
    if(k<L-1)
      tis.push_back(tis_[k]);
    Uis.push_back(.5*Uis_[k]);
    if(Vis_!=0)
      Vis.push_back(Vis_[k]);
    if(muis_!=0)
      muis.push_back(muis_[k]);
    if(Vintis_!=0)
      Vintis.push_back(Vintis_[k]);
  }
  // Construct the MPO
  initZ();
  initHMPO();
}

void auxFillVector(vector<double>& vec,int L,const string name){
  //cout<<"auxFillVector "<<name<<" "<<vec<<" up to "<<L<<endl;
  if(vec.size()<L){
    int len=vec.size();
    double t_=0.;
    if(len!=0) t_=*vec.end();
    cout<<"BoseHubbardHamiltonian constructor, filling in vector "<<name<<", length "
	<<len<<" with last value "<<t_<<endl;
    for(int k=0;k<L-len;k++)vec.push_back(t_);
  }
}

BoseHubbardHamiltonian::BoseHubbardHamiltonian(int L_,int N_,const std::vector<double>& tis_,
					       const std::vector<double>& Uis_,
					       const std::vector<double>& muis_,
					       const std::vector<double>& Vis_,
					       const std::vector<double>& Vintis_):
  hamil(L_),L(L_),N(N_),tis(tis_),Uis(Uis_),Vis(Vis_),muis(muis_),Vintis(Vintis_){
  // check if some vector has to be stretched
  //  cout<<"BoseHubbardHamiltonian constructor with vectors!"<<endl;
  if(tis.size()<L&&tis.size()>1){ // length 1 is special: treated as homogeneous
    auxFillVector(tis,L-1,"t");
  }
  auxFillVector(Uis,L,"U");
  for(int k=0;k<Uis.size();k++) Uis[k]=.5*Uis[k];
  if(muis.size()>0) auxFillVector(muis,L,"mu");
  if(Vis.size()>0) auxFillVector(Vis,L,"V");
  if(Vintis.size()>0) auxFillVector(Vintis,L,"Vint");
  // Construct the MPO
  initZ();
  initHMPO();
}

void BoseHubbardHamiltonian::initZ(){
  int d=N+1; // physical dimension
  // basic spin operators appearing
  // op0 -> identity on d
  // op1 -> n(n-1)
  // op2 -> b+
  // op3 -> b
  // op4 -> n (for mu and penalty term)
  // if penalty, then n^2 too
  int nrOps=5;
  Z=mwArray(Indices(nrOps,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++){
    // fill in
    Z.setElement(ONE_c,Indices(0,i1,i1)); // identity
    Z.setElement(ONE_c*(i1*(i1-1)),Indices(1,i1,i1)); // diagonal N(N-1)
    if(i1>0)
      Z.setElement(ONE_c*sqrt(i1),Indices(2,i1,i1-1)); // b+
    if(i1<N)
      Z.setElement(ONE_c*sqrt(i1+1),Indices(3,i1,i1+1)); // b
    Z.setElement(ONE_c*(i1),Indices(4,i1,i1)); // diagonal N
  }
  Z.reshape(Indices(nrOps,d*d));
  //cout<<"Initialized Z:"<<Z<<endl;
}

//void BoseHubbardHamiltonian::initHMPO(double* tis,double* Uis,
//				      double* muis,double* Vis){
void BoseHubbardHamiltonian::initHMPO(){
  cout<<"initHMPO()"<<endl;
  int D=4; // bond dimension of the MPO
  if(Vintis.size()!=0&&Vintis[0]!=0) D++; // one more if NN interactions
  int d=N+1; // physical dimension

  double t_(0.),U_,mu_,V_,Vint_;
  double t_prev(0.); // since I split the t term left to right,
  //I need to check both
  
  for(int k=0;k<L;k++){
    // decide coefficients, whether they depend on position or not
    if(tis.size()==1){ // they are all constant
      t_=tis[0];t_prev=tis[0]; // don't care that the last one does not enter
      V_=Vis[0];mu_=muis[0];U_=Uis[0];Vint_=Vintis[0];
    }
    else{
      if(k<L-1) t_=tis[k];
      if(k>0) t_prev=abs(tis[k-1]); // sign always on left
      U_=Uis[k];
      if(muis.empty()) mu_=0.;
      else{
	if(muis.size()!=1) mu_=muis[k];
	else mu_=muis[0];
      }
      if(Vis.size()!=1) V_=Vis[k];
      else V_=Vis[0];
      if(Vintis.empty()) Vint_=0;
      else{
	if(Vintis.size()!=1) Vint_=Vintis[k];
	else Vint_=Vintis[0];
      }
    }
    int signt=1;
    int signVint=1;
    if(t_<0){
      t_=-t_;signt=-1;
    }
    if(Vint_<0){
      Vint_=-Vint_;signVint=-1;
    }
    int Dl=D,Dr=D;
    if(k==0) Dl=1;
    if(k==L-1) Dr=1;
    cout<<"For site k="<<k<<" size should be "<<Indices(Dl,Dr,5)<<
    "; using t="<<t_<<", U="<<U_<<endl;
    mwArray C(Indices(Dl,Dr,5));C.fillWithZero();
    if(k==0){
      C.setElement(ONE_c,Indices(0,0,0));
      //C.setElement(ONE_c*V_,Indices(0,D-1,4));
      C.setElement(ONE_c*U_,Indices(0,D-1,1));
      C.setElement(ONE_c*(-1.)*mu_+ONE_c*V_,Indices(0,D-1,4));
      C.setElement(I_c*signt*sqrt(t_),Indices(0,1,2));
      C.setElement(I_c*signt*sqrt(t_),Indices(0,2,3));
      if(Vint_!=0)
	C.setElement(signVint*sqrt(Vint_)*ONE_c,Indices(0,3,4));
    }
    else if(k==L-1){
      C.setElement(ONE_c,Indices(D-1,Dr-1,0));
      //C.setElement(ONE_c*V_,Indices(0,Dr-1,4));
      C.setElement(ONE_c*U_,Indices(0,Dr-1,1));
      C.setElement(ONE_c*(-1.)*mu_+ONE_c*V_,Indices(0,Dr-1,4));
      C.setElement(I_c*sqrt(t_prev),Indices(1,Dr-1,3));
      C.setElement(I_c*sqrt(t_prev),Indices(2,Dr-1,2));
      if(Vint_!=0)
	C.setElement(sqrt(Vint_)*ONE_c,Indices(3,Dr-1,4));
    }
    else{ // intermediate ones
      C.setElement(ONE_c,Indices(0,0,0));
      C.setElement(ONE_c,Indices(D-1,D-1,0));
      //C.setElement(ONE_c*V_,Indices(0,D-1,4));
      C.setElement(ONE_c*U_,Indices(0,D-1,1));
      C.setElement(ONE_c*((-1)*mu_+V_),Indices(0,D-1,4));
      C.setElement(I_c*signt*sqrt(t_),Indices(0,1,2));
      C.setElement(I_c*signt*sqrt(t_),Indices(0,2,3));
      C.setElement(I_c*sqrt(t_prev),Indices(1,D-1,3));
      C.setElement(I_c*sqrt(t_prev),Indices(2,D-1,2));
      if(Vint_!=0){
	C.setElement(signVint*sqrt(Vint_)*ONE_c,Indices(0,3,4));
	C.setElement(sqrt(Vint_)*ONE_c,Indices(3,D-1,4));
      }
    }
    // Now contract MPS and operators and set proper index ordering
    C.reshape(Indices(Dl*Dr,5));
    //cout<<"C for site "<<k<<"="<<C<<endl;
    mwArray res=C*Z;
    //cout<<"Computed res "<<res<<endl;
    res.reshape(Indices(Dl,Dr,d,d));
    // Create and store an Operator
    hamil.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);

  }
}

void BoseHubbardHamiltonian::getTwoBodyTermExponential(mwArray& Ol,
						       mwArray& Or,
						       complex_t delta,
						       int pos){
  mwArray H12;
  int d=N+1;
  computeTwoBodyTerm(H12,pos);
  //H12.multiplyLeft(mwArray(delta));
  // Now take the matrix exponential
  mwArray expH;
  wrapper::expm(H12,expH,delta);
  //cout<<"Exponential of the h12="<<expH<<endl;

  // Obtained with indices (i'j')(ij) => permute to (i'i)(j'j) for svd
  expH.reshape(Indices(d,d,d,d));
  expH.permute(Indices(1,3,2,4));
  expH.reshape(Indices(d*d,d*d));
  // And compute the SVD
  mwArray S; // place for the singular values
  int nr=0;
  double tol=0.;
  wrapper::svd(expH,tol,nr,Ol,S,Or);
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

}

void BoseHubbardHamiltonian::computeTwoBodyTerm(mwArray& result,int pos) 
  const{
  double t_,U1_,U2_,mu1_,mu2_,V1_,V2_,Vint_;  
  if(tis.size()>1){
    t_=tis[pos];
    U1_=Uis[pos];U2_=Uis[pos+1];
    if(Vis.size()>1){
      V1_=Vis[pos];V2_=Vis[pos+1];
    }
    else{
      V1_=V2_=Vis[0];
    }
    if(muis.size()>1){
      mu1_=muis[pos];mu2_=muis[pos+1];
    }
    else{
      mu1_=mu2_=muis[0];
    }
    if(Vintis.size()>1){
      Vint_=Vintis[pos];
    }
    else{
      Vint_=Vintis[0];
    }
  }
  else{ // TI case
    t_=tis[0];U1_=U2_=Uis[0];V1_=V2_=Vis[0];mu1_=mu2_=muis[0];
    Vint_=Vintis[0];
  }
  int d=N+1;
  static bool firstCall(true);
  static mwArray twobody,onebodyU1,onebodyU2;
  static mwArray onebodymu1,onebodymu2;
  static mwArray twobodyInt;
  if(firstCall){
    //cout<<"First time computeTwoBodyTerm is called"<<endl;
    mwArray aCreat(Indices(d,d));
    for(int k=0;k<N;k++){
      aCreat.setElement(sqrt(k+1),0.,Indices(k+1,k));
    }
    mwArray aDestr(Hconjugate(aCreat));
    //cout<<"Created aCreat"<<aCreat.getDimensions()<<", and aDestr"
    //<<aDestr.getDimensions()<<endl;
    //cout<<"Created a+="<<aCreat<<", and a="<<aDestr<<endl;
    mwArray Nop=aCreat*aDestr;
    mwArray Id=identityMatrix(N+1);
    mwArray NopNopminus1=Nop*(Nop-Id);
    //    cout<<"Created Nop"<<Nop.getDimensions()<<" and Id"<<Id.getDimensions()<<endl;
    //cout<<"Created Nop="<<Nop<<" and N(N-1)="<<NopNopminus1<<endl;
    // left one-body term with U
    onebodyU1=kron(NopNopminus1,Id);
    // right one-body term with U
    onebodyU2=kron(Id,NopNopminus1);
    //cout<<"Created onebodyU1"<<onebodyU1.getDimensions()
    //	<<" and onebodyU2"<<onebodyU2.getDimensions()<<endl;
    // left and right one-body term with mu and V
    onebodymu1=kron(Nop,Id);
    onebodymu2=kron(Id,Nop);
    //cout<<"Created onebodymu1"<<onebodymu1.getDimensions()
    //	<<" and onebodymu2"<<onebodymu2.getDimensions()<<endl;
    // two-body term
    mwArray aDaggera=kron(aCreat,aDestr);
    twobodyInt=kron(Nop,Nop);
    //cout<<"Created a+ a"<<aDaggera.getDimensions()<<endl;
    twobody=aDaggera;
    aDaggera.Hconjugate();
    twobody=twobody+aDaggera;
    //cout<<"Created a+ a+a a+"<<twobody.getDimensions()<<endl;
  }
  firstCall=false;

  // For the first and last sites, the single body terms
  // do not have the factor .5, as there is no other term including them
  double factorL=.5; 
  if(pos==0) factorL=1.;
  double factorR=.5;
  if(pos+1==L-1) factorR=1.;

  result=-t_*twobody+factorL*U1_*onebodyU1+factorR*U2_*onebodyU2
    +factorL*(-mu1_+V1_)*onebodymu1+factorR*(-mu2_+V2_)*onebodymu2
    +Vint_*twobodyInt;

}

void BoseHubbardHamiltonian::getExponentialMPO(MPO& expH,complex_t delta,bool even){
  mwArray Ol,Or;
  int kfirst=even?0:1;
  int klast=(L-1)-(L-kfirst)%2; // last occupied by the loop
  for(int k=kfirst;k<klast;k+=2){
    bool change=(Ol.isEmpty())|| // the absolute first one
      (k==2)|| // the second one is different from the first one!!!!
      (k==L-2)||  // the last of the chain is always special (.5 factor on right)
      (tis.size()>1&& // parameters are not constant
       (tis[k]!=tis[k-2]||Uis[k]!=Uis[k-2]|| // and any of them is different
	Vis[k]!=Vis[k-2]||muis[k]!=muis[k-2]||
	Uis[k+1]!=Uis[k-1]||Vis[k+1]!=Vis[k-1]||muis[k+1]!=muis[k-1]));
    if(change){ // I need to recompute the Operators
      cout<<"Recomputing exponential operator for pos "<<k<<endl;
      getTwoBodyTermExponential(Ol,Or,delta,k);
      expH.setOp(k,new Operator(Ol),true);
      expH.setOp(k+1,new Operator(Or),true);
    }
    else{ // I can just copy the pointers to the ones before
      expH.setOp(k,&expH.getOp(k-2));
      expH.setOp(k+1,&expH.getOp(k-1));
    }
    //cout<<"Set pos "<<k<<" to Ol and pos "<<k+1<<" to Or"<<endl;
  }
  // Fill in the edges with identity operators
  mwArray ident=identityMatrix(N+1);
  ident.reshape(Indices(N+1,1,N+1,1)); // just dx1xdx1
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

void BoseHubbardHamiltonian::getExponentialMPOeven(MPO& expHe,complex_t delta){
  getExponentialMPO(expHe,delta,true);
}
void BoseHubbardHamiltonian::getExponentialMPOodd(MPO& expHo,complex_t delta){
  getExponentialMPO(expHo,delta,false);
}
