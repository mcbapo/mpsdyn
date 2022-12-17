/**
   \file IsingHamiltonian.cpp
   Hamiltonian for the Ising model (with parallel and transverse magnetic field)
   
   \author Mari-Carmen Banuls
   \date 20/10/2011

*/

#include "IsingHamiltonian.h"
#include "Indices.h"
//#include "complex.h"

using namespace std;
using namespace shrt;

IsingHamiltonian::IsingHamiltonian(int N_,int d_,double J0_,double g0_,
				   double h0_,double K_):
  hamil(N_),N(N_),d(d_),J0(J0_),g0(g0_),h0(h0_),t(0),timedep(0),
  J(NULL),g(NULL),h(NULL),Z(),K(K_){
  if(d!=2){
    cout<<"Error: IsingHamiltonian only defined for d=2"<<endl;
    exit(212);
  }
  initZ();
  t=0.;timedep=false;
  setHMPO(hamil,J0,g0,h0);
  J=g=h=0;
}

IsingHamiltonian::~IsingHamiltonian(){};

void IsingHamiltonian::registerTimeDependence(tdepfun J_,tdepfun g_,
					      tdepfun h_){
  if(J_!=0){
    cout<<"Registering J(t) dependency "<<endl;
    J=J_;timedep=1;
    J0=(*J)(0);
  }
  if(g_!=0){
    g=g_;timedep=1;
    g0=(*g)(0);
  }
  if(h_!=0){
    cout<<"Registering h(t) dependency "<<endl;
    h=h_;timedep=1;
    h0=(*h)(0);
  }
}

void IsingHamiltonian::initZ(){
  // basic spin operators appearing
  mwArray sig0=identityMatrix(d);
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sig3(Indices(d,d),dataz);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig1(Indices(d,d),datax);
  Z=mwArray(Indices(3,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      //Z.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Z.setElement(sig3.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      //Z.setElement(sig3.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      Z.setElement(sig1.getElement(Indices(i1,i2)),Indices(2,i1,i2));
    }
  Z.reshape(Indices(3,d*d));
}

void IsingHamiltonian::setHMPO(MPO& mpo,double J,double g,double h){
  int D=3;
  for(int k=0;k<N;k++){
    // only first, second and last are different
    if(k>1&&k<N-1){
      mpo.setOp(k,&mpo.getOp(k-1));
    }
    else{
      //cout<<"Creating new Operator for site "<<k<<" of IsingH MPO"<<endl;
      int Dl=D; int Dr=D;
      if(k==0) Dl=1;
      if(k==N-1) Dr=1;
      int redefK=k+1-N/2;
      mwArray C(Indices(Dl,Dr,3));C.fillWithZero();
      complex_t effJ=(J<0)?sqrt(abs(J))*I_c:sqrt(J)*ONE_c;
      complex_t effG=g*ONE_c;
      complex_t effH=h*ONE_c;
      complex_t effK=(K/N)*ONE_c;
      //cout<<"Coefficients: effJ="<<effJ<<", effG="<<effG<<", effH="<<effH
      //  <<endl;
      // Identity
      if(k!=N-1)
	C.setElement(ONE_c,Indices(0,0,0)); // just pass 1 (before anything)
      if(k!=0)
	C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // just pass 1 (after)
      // sigmax terms (in all of them)
      C.setElement(effG,Indices(0,Dr-1,2));
      // sigmaz terms (in all of them)
      C.setElement(effH,Indices(0,Dr-1,1));
      // independent offset term, in all, too
      C.setElement(effK,Indices(0,Dr-1,0));
      // left sigz term (all but last)
      if(k!=N-1)
	C.setElement(effJ,Indices(0,1,1));
      // right sigz term (all but first)
      if(k!=0)
	C.setElement(effJ,Indices(1,Dr-1,1));
      
      // Now contract and give the operator
      // Now contract MPS and operators and set proper index ordering
      C.reshape(Indices(Dl*Dr,3));
      //cout<<"C for site "<<k<<"="<<C<<endl;
      mwArray res=C*Z;
      //cout<<"Computed res "<<res<<endl;
      res.reshape(Indices(Dl,Dr,d,d));
      res.permute(Indices(3,1,4,2));
      //cout<<"  with size: "<<res.getDimensions()<<endl;
      // Create and store an Operator
      mpo.setOp(k,new Operator(res),true);
    }
  }
}

const MPO& IsingHamiltonian::getHMPO(double t_){
  if(timedep&&t!=t_){
    t=t_;
    setHMPO(hamil,J(t),g(t),h(t));
  }
  return hamil;
}

void IsingHamiltonian::getUMPO(MPO& mpo,double delta,double t,bool imag,
			       int orderT,int nrOp){
  mpo.initLength(N);
  if(orderT!=2){
    cout<<"Error: order Trotter "<<orderT<<" nor supported yet (only 2, for now)"<<endl;
    exit(212);
  }
  // Set the edge operators, which are different
  complex_t eps=imag?-delta*ONE_c:-delta*I_c;
  mpo.setOp(0,new Operator(constructU2T(t,eps,0)),true);
  mpo.setOp(N-1,new Operator(constructU2T(t,eps,N-1)),true);
  for(int k=1;k<N-1;k++){
    if(k==1)
      mpo.setOp(k,new Operator(constructU2T(t,eps,k)),true);
    else
      mpo.setOp(k,&mpo.getOp(k-1));
  }
}  


void IsingHamiltonian::getUMPO(MPO& mpo,complex_t deltaC,double t,
			       int orderT,int nrOp){
  mpo.initLength(N);
  if(orderT!=2){
    cout<<"Error: order Trotter "<<orderT<<" nor supported yet (only 2, for now)"<<endl;
    exit(212);
  }
  // Set the edge operators, which are different
  mpo.setOp(0,new Operator(constructU2T(t,deltaC,0)),true);
  mpo.setOp(N-1,new Operator(constructU2T(t,deltaC,N-1)),true);
  for(int k=1;k<N-1;k++){
    if(k==1)
      mpo.setOp(k,new Operator(constructU2T(t,deltaC,k)),true);
    else
      mpo.setOp(k,&mpo.getOp(k-1));
  }
}  

mwArray IsingHamiltonian::isingU2MPO(complex_t eps,int pos) const{
  int D=2;int Dl=D;int Dr=D;
  if(eps==ZERO_c){
    mwArray res=identityMatrix(d);
    res.reshape(Indices(d,1,d,1));
    return res;
  }
  mwArray C(Indices(D,D,3));
  if(pos!=0&&pos!=N-1){// middle term
    C.setElement(cosh(eps),Indices(0,0,0));
    C.setElement(sinh(eps),Indices(1,1,0));
    C.setElement(sqrt(cosh(eps)*sinh(eps)),Indices(0,1,1));
    C.setElement(sqrt(cosh(eps)*sinh(eps)),Indices(1,0,1));
  }
  else{
    C.resize(Indices(1,D,3));
    C.setElement(sqrt(cosh(eps)),Indices(0,0,0));
    C.setElement(sqrt(sinh(eps)),Indices(0,1,1));
    if(pos==0)
      Dl=1;
    else
      Dr=1;
  }
  //cout<<"Coefficients in the exponent for site "<<pos<<" eps="<<eps<<endl;
  C.reshape(Indices(Dl*Dr,3));
  //cout<<"C tensor "<<reshape(resize(C,Indices(Dl*Dr,2)),Indices(Dl,Dr,2))<<endl;
  mwArray res=C*Z;
  res.reshape(Indices(Dl,Dr,d,d));
  res.permute(Indices(3,1,4,2));
  //cout<<" IsingHamiltonian::isingU2MPO("<<pos<<")->"<<res.getDimensions()
  //  <<res<<endl;
  return res;
}

mwArray IsingHamiltonian::isingU1MPO(complex_t epsilon,double g,
				     double h) const{
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigX(Indices(d,d),datax);
  mwArray sigZ(Indices(d,d),dataz);
  mwArray id=identityMatrix(d);
  double norm=sqrt(g*g+h*h);
  mwArray res=cosh(norm*epsilon)*id;
  if(norm!=0) // else only identity
    res=res+(h/norm)*sinh(norm*epsilon)*sigZ+
      (g/norm)*sinh(norm*epsilon)*sigX;
  res.reshape(Indices(d,1,d,1));
  //cout<<" IsingHamiltonian::isingU1MPO->"<<res.getDimensions()<<endl;
  return res;
}
 
mwArray IsingHamiltonian::constructU2T(double t,complex_t delta,int pos){
  //cout<<"IsingHamiltonian::constructU2T("<<delta<<","<<pos<<")"<<endl;
  double Jt=J==0?J0:J(t);
  double gt=g==0?g0:g(t);
  double ht=h==0?h0:h(t);
  mwArray U1=isingU1MPO(.5*delta,gt,ht);
  mwArray U2=isingU2MPO(Jt*delta,pos);
  int Dl=U2.getDimension(1);
  int Dr=U2.getDimension(3);
  U1.reshape(Indices(d,d));
  //cout<<"U1 is "<< U1<<endl;
  U2.reshape(Indices(d,Dl*d*Dr));
  //cout<<"U2 is "<< U2<<endl;
  mwArray res=U1*U2;
  res.reshape(Indices(d*Dl,d,Dr));
  res.permute(Indices(1,3,2));
  res.reshape(Indices(d*Dl*Dr,d));
  res=res*U1;
  res.reshape(Indices(d*Dl,Dr,d));
  res.permute(Indices(1,3,2));
  res.reshape(Indices(d,Dl,d,Dr));
  //cout<<"IsingHamiltonian::constructU2T("<<delta<<","<<pos<<")->"
  //  <<res.getDimensions()<<endl;
  return res;
}

mwArray IsingHamiltonian::getUOperator(double delta,double t,bool imag,
				       int orderT,int nrOp){
  if(orderT!=2){
    cout<<"Error: order Trotter "<<orderT
	<<" not supported yet (only 2, for now)"<<endl;
    exit(212);
  }
  else if(nrOp>0){
    cout<<"Incorrect nrOp("<<nrOp<<") for order Trotter 2, ignoring"<<endl;
  }
  // Set the edge operators, which are different
  complex_t eps=imag?-delta*ONE_c:-delta*I_c;
  return constructU2T(t,eps);
}

void IsingHamiltonian::getCommutatorMPO(MPO& mpo,double t){
  int D=4; // bond dimension of the MPO
  mwArray Z; // The operators
  initZdoub(Z);
  int nrOps=Z.getDimension(0); // Nr of operators (5)
  // coefficients, in the specified time
  double J_(J0),g_(g0),h_(h0);
  if(timedep){
    J_=J(t);
    g_=g(t);
    h_=h(t);
  }
  complex_t effJ=(J_<0)?sqrt(abs(J_))*I_c:sqrt(J_)*ONE_c;

  // Now construct the coefficients for every term. I need different
  // ones for the first and last sites, but the rest are all the same
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
  C.setElement(h_*ONE_c,Indices(0,D-1,1));
  C.setElement(-h_*ONE_c,Indices(0,D-1,2));
  C.setElement(g_*ONE_c,Indices(0,D-1,3));
  C.setElement(-g_*ONE_c,Indices(0,D-1,4));
  // The same for first site 
  C1.setElement(h_*ONE_c,Indices(0,D-1,1));
  C1.setElement(-h_*ONE_c,Indices(0,D-1,2));
  C1.setElement(g_*ONE_c,Indices(0,D-1,3));
  C1.setElement(-g_*ONE_c,Indices(0,D-1,4));
  // And also for the last site
  CN.setElement(h_*ONE_c,Indices(0,0,1));
  CN.setElement(-h_*ONE_c,Indices(0,0,2));
  CN.setElement(g_*ONE_c,Indices(0,0,3));
  CN.setElement(-g_*ONE_c,Indices(0,0,4));
  // And now the two-body terms
  // J ( sigma_z x Id ) (x) ( sigma_z x Id ) 
  C.setElement(effJ,Indices(0,2,1));
  C.setElement(effJ,Indices(2,D-1,1));
  C1.setElement(effJ,Indices(0,2,1));
  CN.setElement(effJ,Indices(2,0,1));
  // -J ( Id x sigma_z ) (x) ( Id x sigma_z ) 
  C.setElement(effJ*I_c,Indices(0,1,2));
  C.setElement(effJ*I_c,Indices(1,D-1,2));
  C1.setElement(effJ*I_c,Indices(0,1,2));
  CN.setElement(effJ*I_c,Indices(1,0,2));

  // Now I reshape, contract with Z and set the indices in proper order
  C.reshape(Indices(D*D,nrOps));
  C1.reshape(Indices(1*D,nrOps));
  CN.reshape(Indices(D*1,nrOps));
  Z.reshape(Indices(nrOps,d*d*d*d));
  C.multiplyRight(Z);C.reshape(Indices(D,D,d*d,d*d));C.permute(Indices(3,1,4,2));
  C1.multiplyRight(Z);C1.reshape(Indices(1,D,d*d,d*d));C1.permute(Indices(3,1,4,2));
  CN.multiplyRight(Z);CN.reshape(Indices(D,1,d*d,d*d));CN.permute(Indices(3,1,4,2));

  // Construct the three operators
  Operator* Op1=new Operator(C1);
  Operator* Op=new Operator(C);
  Operator* OpN=new Operator(CN);

  // And set them to the adequate positions in the MPO!
  mpo.setOp(0,Op1,true); // mpo is now resposible of this pointer
  mpo.setOp(N-1,OpN,true);
  if(N>2){
    mpo.setOp(1,Op,true);
    for(int k=2;k<N-1;k++){
      // For all the remaining middle sites, I do not need to create
      // new Operators, but just point to the same one!
      mpo.setOp(k,Op,false);      
    }
  }
  else{
    // delete the extra Op, to avoid memory leaks, in this very silly case of only two sites
    delete Op;
  }
}




void IsingHamiltonian::initZdoub(mwArray& Z){
  // Basic pieces
  mwArray sig0=identityMatrix(d);//identity
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),datax);//sigmax
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
  // *** Now the ones appearing in single-body terms
  // (1) sigZ (x) Id
  mwArray term1; 
  constructOperatorProduct(term1,sigZ,sig0);
  // (2)  Id (x) sigZ -> just a permutation of the previous one
  mwArray term2;
  constructOperatorProduct(term2,sig0,sigZ);
  // (3) sigX (x) Id
  mwArray term3;
  constructOperatorProduct(term3,sigX,sig0);
  // (4)  Id (x) sigX -> just a permutation of the previous one
  mwArray term4;
  constructOperatorProduct(term4,sig0,sigX);

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
    }
}

void IsingHamiltonian::constructOperatorProduct(mwArray& result,const mwArray& opA,
						const mwArray& opB){
#ifdef CHECKDIMS
  if(opA.getDimensions()!=opB.getDimensions()){
    cout<<"Error in CoherentDissipationLiouvillian::constructOperatorProduct for"
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
