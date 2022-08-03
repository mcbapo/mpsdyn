/**
   \file FAHamiltonian.cpp
   Hamiltonian for a stochastic model
   
   \author Mari-Carmen Banuls
   \date 12/9/2018

*/

#include "FAHamiltonian.h"

using namespace std;
using namespace shrt;

FAHamiltonian::FAHamiltonian(int L_,double s_,double offset_,double c_,double noZero_,bool pbc_):d(2),L(L_),hamil(L),s(s_),c(c_),offset(offset_),noZeroPenalty(noZero_),pbc(pbc_){
  initZ();
  if(pbc)
    initHMPOpbc();
  else
    initHMPO();
}

FAHamiltonian::~FAHamiltonian(){};

void FAHamiltonian::initZ(){
  // basic spin operators appearing
  mwArray sig0=identityMatrix(d);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig1(Indices(d,d),datax);
  complex_t dataN[]={ZERO_c,ZERO_c,ZERO_c,ONE_c};
  mwArray operN(Indices(d,d),dataN);
  mwArray oper0=sig0-operN; // projector onto 0, in case it is suppressed
  //  mwArray operB=exp(-s)*sqrt(c*(1-c))*sig1-c*sig0-(1-2*c)*operN;
  Z=mwArray(Indices(4,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(operN.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Z.setElement(sig1.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      Z.setElement(oper0.getElement(Indices(i1,i2)),Indices(3,i1,i2));
    }
  Z.reshape(Indices(4,d*d));
}

void FAHamiltonian::initHMPO(){
  int D=noZeroPenalty!=0.?5:4;
  double factX=-exp(-s)*sqrt(c*(1-c)); // the factor with sigX
  bool limSneg=0;
  // A special case is the limit of very large s!
  if(c==0.5&&s<=-1111){
    cout<<"Creating FA Hamiltonian in the limit s->-infty! (no s or c dependence)"<<endl;
    limSneg=1;factX=-1;
  }

  for(int k=0;k<L;k++){
    // only first, second and last are distinct
    if(k>1&&k<L-1){
      hamil.setOp(k,&hamil.getOp(k-1));
    }
    else{
      //cout<<"Creating new Operator for site "<<k<<" of FAH MPO"<<endl;
      int Dl=D; int Dr=D;
      if(k==0) Dl=1;
      if(k==L-1) Dr=1;
      mwArray C(Indices(Dl,Dr,4));C.fillWithZero();
      // Set elements
      if(k!=L-1){ // all but last: first in pair
	C.setElement(ONE_c,Indices(0,0,0)); // nothing yet
	// First of the pair
	C.setElement(ONE_c,Indices(0,1,1)); // term with n first
	C.setElement(factX*ONE_c,Indices(0,2,2)); // term with sigx first
	if(k!=0){
	  if(!limSneg)
	    C.setElement(2*c*ONE_c,Indices(0,Dr-1,1)); // single n term
	}
      }
      else{
	if(!limSneg)
	  C.setElement(c*ONE_c,Indices(0,Dr-1,1)); // single n term
      }
      if(k!=0){ 
	// After everything
	C.setElement(ONE_c,Indices(Dl-1,Dr-1,0));
	// Second of the pair
	C.setElement(factX*ONE_c,Indices(1,Dr-1,2)); // nX
	C.setElement(ONE_c,Indices(2,Dr-1,1)); // Xn
	if(!limSneg)
	  C.setElement(2*(1-2*c)*ONE_c,Indices(1,Dr-1,1)); // nn
      }
      else{
	if(!limSneg)
	  C.setElement(c*ONE_c,Indices(0,Dr-1,1)); // single n term
      }
      if(offset!=0)
	C.setElement((offset/L)*ONE_c,Indices(0,Dr-1,0)); // the offset term
      if(noZeroPenalty!=0){
	if(k==0)
	  C.setElement(noZeroPenalty*ONE_c,Indices(0,3,3));
	else
	  if(k==L-1)
	    C.setElement(ONE_c,Indices(3,0,3));
	  else
	    C.setElement(ONE_c,Indices(3,3,3));
      }
      // Reshape, multiply operators and set in MPO
      C.reshape(Indices(Dl*Dr,4));
      C.multiplyRight(Z);
      C.reshape(Indices(Dl,Dr,d,d));
      C.permute(Indices(3,1,4,2));
      hamil.setOp(k,new Operator(C),true);
    }
  }
}

// To keep the other untouched, I do this extra method, although would not need it (could modify the other one with more cases)
void FAHamiltonian::initHMPOpbc(){
  int D=noZeroPenalty!=0.?7:6;
  double factX=-exp(-s)*sqrt(c*(1-c)); // the factor with sigX
  bool limSneg=0;
  // A special case is the limit of very large s!
  if(c==0.5&&s<=-1111){
    cout<<"Creating FA Hamiltonian in the limit s->-infty! (no s or c dependence)"<<endl;
    limSneg=1;
  }

  // basic spin operators appearing
  mwArray sig0=identityMatrix(d);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig1(Indices(d,d),datax);
  complex_t dataN[]={ZERO_c,ZERO_c,ZERO_c,ONE_c};
  mwArray operN(Indices(d,d),dataN);
  mwArray oper0=sig0-operN; // projector onto 0, in case it is suppressed
  mwArray oper1=factX*sig1+(1-2*c)*operN+c*sig0; // includes the minus sign
  if(limSneg) oper1=-1.*sig1; // includes the minus sign
  //  mwArray operB=exp(-s)*sqrt(c*(1-c))*sig1-c*sig0-(1-2*c)*operN;
  Z=mwArray(Indices(4,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(operN.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Z.setElement(oper1.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      Z.setElement(oper0.getElement(Indices(i1,i2)),Indices(3,i1,i2));
    }
  Z.reshape(Indices(4,d*d));

  
  for(int k=0;k<L;k++){
    // only first, second and last are distinct
    if(k>1&&k<L-1){
      hamil.setOp(k,&hamil.getOp(k-1));
    }
    else{
      //cout<<"Creating new Operator for site "<<k<<" of FAH MPO"<<endl;
      int Dl=D; int Dr=D;
      if(k==0) Dl=1;
      if(k==L-1) Dr=1;
      mwArray C(Indices(Dl,Dr,4));C.fillWithZero();
      // Set elements
      if(k!=L-1){ // all but last: first in pair
	C.setElement(ONE_c,Indices(0,0,0)); // nothing yet
	// First of the pair
	C.setElement(ONE_c,Indices(0,1,1)); // term with n first
	C.setElement(ONE_c,Indices(0,2,2)); // term with sigx first
	if(k==0){ // also the boundary terms
	  C.setElement(ONE_c,Indices(0,3,1)); // term with n second
	  C.setElement(ONE_c,Indices(0,4,2)); // term with sigx second
	}
	else{ // pass those on
	  C.setElement(ONE_c,Indices(3,3,0)); // term with n second
	  C.setElement(ONE_c,Indices(4,4,0)); // term with sigx second
	}
      }
      else{ // on the last, put the first part of boundary
	  C.setElement(ONE_c,Indices(3,Dr-1,2)); // term with n second
	  C.setElement(ONE_c,Indices(4,Dr-1,1)); // term with sigx second
      }
      if(k!=0){ 
	// After everything
	C.setElement(ONE_c,Indices(Dl-1,Dr-1,0));
	// Second of the pair
	C.setElement(ONE_c,Indices(1,Dr-1,2)); // nX
	C.setElement(ONE_c,Indices(2,Dr-1,1)); // Xn
      }
      if(offset!=0)
	C.setElement((offset/L)*ONE_c,Indices(0,Dr-1,0)); // the offset term
      if(noZeroPenalty!=0){
	if(k==0)
	  C.setElement(noZeroPenalty*ONE_c,Indices(0,5,3));
	else
	  if(k==L-1)
	    C.setElement(ONE_c,Indices(5,0,3));
	  else
	    C.setElement(ONE_c,Indices(5,5,3));
      }
      // Reshape, multiply operators and set in MPO
      C.reshape(Indices(Dl*Dr,4));
      C.multiplyRight(Z);
      C.reshape(Indices(Dl,Dr,d,d));
      C.permute(Indices(3,1,4,2));
      hamil.setOp(k,new Operator(C),true);
    }
  }
}

void FAHamiltonian::getExponentialMPOeven(MPO& expHe,complex_t delta) const {
 getExponentialMPO(expHe,delta,true);
}
void FAHamiltonian::getExponentialMPOodd(MPO& expHo,complex_t delta) const {
 getExponentialMPO(expHo,delta,false);
}

void FAHamiltonian::getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,
						      int pos) const {
  mwArray H12;
  computeTwoBodyTerm(H12,pos);
  getTwoBodyTermExponential(Ol,Or,H12,delta);
}

void FAHamiltonian::getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,
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
  //putForMatlab(cout,Ol,"Ol");
  //putForMatlab(cout,Or,"Or");
}

void FAHamiltonian::getExponentialMPO(MPO& expH,complex_t delta,bool even) const {
  mwArray Ol,Or;
  expH.initLength(L);
  int kfirst=even?0:1;
  int klast=(L-1)-(L-kfirst)%2; // last occupied by the loop
  for(int k=kfirst;k<klast;k+=2){
    getTwoBodyTermExponential(Ol,Or,delta,k);
    expH.setOp(k,new Operator(Ol),true);
    expH.setOp(k+1,new Operator(Or),true);
    // TODO: could save some space using pointers
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



void FAHamiltonian::getCommutatorMPO(MPO& mpo){
  cout<<"Not ready yet!"<<endl;
  exit(1);
  // int D=4; // bond dimension of the MPO
  // mwArray Z; // The operators
  // initZdoub(Z);
  // int nrOps=Z.getDimension(0); // Nr of operators (5)
  // // coefficients, in the specified time
  // double J_(J0),g_(g0),h_(h0);
  // if(timedep){
  //   J_=J(t);
  //   g_=g(t);
  //   h_=h(t);
  // }
  // complex_t effJ=(J_<0)?sqrt(abs(J_))*I_c:sqrt(J_)*ONE_c;

  // // Now construct the coefficients for every term. I need different
  // // ones for the first and last sites, but the rest are all the same
  // mwArray C1(Indices(1,D,nrOps)); // the first site
  // mwArray C(Indices(D,D,nrOps)); // the middle site
  // mwArray CN(Indices(D,1,nrOps)); // the last site

  // // Identity terms when nothing has happened yet
  // C.setElement(ONE_c,Indices(0,0,0));
  // C1.setElement(ONE_c,Indices(0,0,0));
  // // Identity terms after the end
  // C.setElement(ONE_c,Indices(D-1,D-1,0));
  // CN.setElement(ONE_c,Indices(D-1,0,0));
  // // Single body terms (one for each operator, with proper weights)
  // C.setElement(h_*ONE_c,Indices(0,D-1,1));
  // C.setElement(-h_*ONE_c,Indices(0,D-1,2));
  // C.setElement(g_*ONE_c,Indices(0,D-1,3));
  // C.setElement(-g_*ONE_c,Indices(0,D-1,4));
  // // The same for first site 
  // C1.setElement(h_*ONE_c,Indices(0,D-1,1));
  // C1.setElement(-h_*ONE_c,Indices(0,D-1,2));
  // C1.setElement(g_*ONE_c,Indices(0,D-1,3));
  // C1.setElement(-g_*ONE_c,Indices(0,D-1,4));
  // // And also for the last site
  // CN.setElement(h_*ONE_c,Indices(0,0,1));
  // CN.setElement(-h_*ONE_c,Indices(0,0,2));
  // CN.setElement(g_*ONE_c,Indices(0,0,3));
  // CN.setElement(-g_*ONE_c,Indices(0,0,4));
  // // And now the two-body terms
  // // J ( sigma_z x Id ) (x) ( sigma_z x Id ) 
  // C.setElement(effJ,Indices(0,2,1));
  // C.setElement(effJ,Indices(2,D-1,1));
  // C1.setElement(effJ,Indices(0,2,1));
  // CN.setElement(effJ,Indices(2,0,1));
  // // -J ( Id x sigma_z ) (x) ( Id x sigma_z ) 
  // C.setElement(effJ*I_c,Indices(0,1,2));
  // C.setElement(effJ*I_c,Indices(1,D-1,2));
  // C1.setElement(effJ*I_c,Indices(0,1,2));
  // CN.setElement(effJ*I_c,Indices(1,0,2));

  // // Now I reshape, contract with Z and set the indices in proper order
  // C.reshape(Indices(D*D,nrOps));
  // C1.reshape(Indices(1*D,nrOps));
  // CN.reshape(Indices(D*1,nrOps));
  // Z.reshape(Indices(nrOps,d*d*d*d));
  // C.multiplyRight(Z);C.reshape(Indices(D,D,d*d,d*d));C.permute(Indices(3,1,4,2));
  // C1.multiplyRight(Z);C1.reshape(Indices(1,D,d*d,d*d));C1.permute(Indices(3,1,4,2));
  // CN.multiplyRight(Z);CN.reshape(Indices(D,1,d*d,d*d));CN.permute(Indices(3,1,4,2));

  // // Construct the three operators
  // Operator* Op1=new Operator(C1);
  // Operator* Op=new Operator(C);
  // Operator* OpN=new Operator(CN);

  // // And set them to the adequate positions in the MPO!
  // mpo.setOp(0,Op1,true); // mpo is now resposible of this pointer
  // mpo.setOp(N-1,OpN,true);
  // if(N>2){
  //   mpo.setOp(1,Op,true);
  //   for(int k=2;k<L-1;k++){
  //     // For all the remaining middle sites, I do not need to create
  //     // new Operators, but just point to the same one!
  //     mpo.setOp(k,Op,false);      
  //   }
  // }
  // else{
  //   // delete the extra Op, to avoid memory leaks, in this very silly case of only two sites
  //   delete Op;
  // }
}




void FAHamiltonian::initZdoub(mwArray& Z){
  cout<<"Not ready yet!"<<endl;
  exit(1);
  // // Basic pieces
  // mwArray sig0=identityMatrix(d);//identity
  // complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  // mwArray sigX(Indices(d,d),datax);//sigmax
  // complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  // mwArray sigZ(Indices(d,d),dataz);//sigmaz

  // // Because my sites are double, the actual individual operators will
  // // be tensor products of these, or of these with the identity So, I
  // // construct the basic products one by one.  (Notice that not all of
  // // them are independent, and I could also do a more efficient
  // // contraction, but this cost is negligible compared to anything
  // // else, as I only need to do this onec when initializing the
  // // problem, so I opt for the clearerst construction here).

  // // *** First of all, we need the identity on both
  // mwArray term0=identityMatrix(d*d);
  // // *** Now the ones appearing in single-body terms
  // // (1) sigZ (x) Id
  // mwArray term1; 
  // constructOperatorProduct(term1,sigZ,sig0);
  // // (2)  Id (x) sigZ -> just a permutation of the previous one
  // mwArray term2;
  // constructOperatorProduct(term2,sig0,sigZ);
  // // (3) sigX (x) Id
  // mwArray term3;
  // constructOperatorProduct(term3,sigX,sig0);
  // // (4)  Id (x) sigX -> just a permutation of the previous one
  // mwArray term4;
  // constructOperatorProduct(term4,sig0,sigX);

  // // Fill in the operators in the Z array, to be contracted with C in
  // // the construction of the MPO. The order is as listed above.
  // int nrOps=5;
  // Z=mwArray(Indices(nrOps,d*d,d*d));
  // for(int d1=0;d1<d*d;d1++)
  //   for(int d2=0;d2<d*d;d2++){
  //     Z.setElement(term0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
  //     Z.setElement(term1.getElement(Indices(d1,d2)),Indices(1,d1,d2));
  //     Z.setElement(term2.getElement(Indices(d1,d2)),Indices(2,d1,d2));
  //     Z.setElement(term3.getElement(Indices(d1,d2)),Indices(3,d1,d2));
  //     Z.setElement(term4.getElement(Indices(d1,d2)),Indices(4,d1,d2));
  //   }
}

void FAHamiltonian::constructOperatorProduct(mwArray& result,const mwArray& opA,
						const mwArray& opB) const{
#ifdef CHECKDIMS
  if(opA.getDimensions()!=opB.getDimensions()){
    cout<<"Error in FAHamiltonian::constructOperatorProduct for"
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

void FAHamiltonian::computeTwoBodyTerm(mwArray& result,int pos) const{
  if(pos>=L-1){
    cout<<"ERROR: FAHamiltonianHamiltonian::computeTwoBodyTerm for site "<<pos<<" when L="<<L<<endl;
    exit(212);
  }
  static mwArray NX,XN,NI,IN,NN;
  static bool firstCall(true);
  if(firstCall){
    // basic spin operators appearing
    mwArray sig0=identityMatrix(d);
    complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    mwArray sig1(Indices(d,d),datax);
    complex_t dataN[]={ZERO_c,ZERO_c,ZERO_c,ONE_c};
    mwArray operN(Indices(d,d),dataN);
    constructOperatorProduct(NX,operN,sig1);
    constructOperatorProduct(XN,sig1,operN);
    constructOperatorProduct(NI,operN,sig0);
    constructOperatorProduct(IN,sig0,operN);
    constructOperatorProduct(NN,operN,operN);
  }
  result=-exp(-s)*sqrt(c*(1-c))*(NX+XN)+2*c*NI+2*(1-2*c)*NN;
  if(pos==0) result=result-c*NI;  
  if(pos==L-1) result=result-c*IN;
}
