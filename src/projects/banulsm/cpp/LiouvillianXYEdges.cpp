
#include "LiouvillianXYEdges.h"
#include "DoubleOperator.h"

using namespace std;
using namespace shrt;


LiouvillianXYEdges::LiouvillianXYEdges(int N_,double gamma_,double h_,
				       double GammaL1_,
				       double GammaL2_,double GammaR3_,
				       double GammaR4_,double scale_):
  d(2),N(N_),h(h_),gamma(gamma_),GammaL1(GammaL1_),GammaL2(GammaL2_),
  GammaR3(GammaR3_),GammaR4(GammaR4_),
  Lmpo(N_),LAdjmpo(N_),trLmpo(N_),isInit(false),scale(scale_){
  initL();
  initLAdjoint();
  initTrL();
}

// LiouvillianXYEdges::LiouvillianXYEdges(int N_,double g_,double gamma_,bool priv):
//   d(2),N(N_),g(g_),gamma(gamma_),Lmpo(N_),LAdjmpo(N_),trLmpo(N_),isInit(false){};

LiouvillianXYEdges::~LiouvillianXYEdges(){
  clear();
}

void LiouvillianXYEdges::clear(){
  Lmpo.clear();
  isInit=false;
  LAdjmpo.clear();
}

void LiouvillianXYEdges::initLAdjoint(){
  cout<<"LiouvillianXYEdges::initLAdjoint()"<<endl;
  // I can only do this if the Lmpo is already prepared
  initL();

  // I will create copies of the edge and middle operators with the
  // adequate rotation of indices (up<->down) and complex conjugation
  // and for the rest of the chain, I just point to the operator in
  // the second site, as L itself does

  LAdjmpo.setRotatedOp(0,&Lmpo.getOp(0),Indices(3,2,1,4),true);
  LAdjmpo.setRotatedOp(N-1,&Lmpo.getOp(N-1),Indices(3,2,1,4),true);
  if(N>2){
    LAdjmpo.setRotatedOp(1,&Lmpo.getOp(1),Indices(3,2,1,4),true);
    for(int k=2;k<N-1;k++)
      LAdjmpo.setOp(k,&LAdjmpo.getOp(1),false);            
  }
}


void LiouvillianXYEdges::initL(){
  if(isInit) return;
  cout<<"LiouvillianXYEdges::initL()"<<endl;
  int D=6; // bond dimension of the MPO
  mwArray Z; // The operators
  initZ(Z);
  int nrOps=Z.getDimension(0); // Nr of operators (13)
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
  // magnetic field in H
  C.setElement(-h*I_c,Indices(0,D-1,5));
  C.setElement(h*I_c,Indices(0,D-1,6));
  C1.setElement(-h*I_c,Indices(0,D-1,5));
  C1.setElement(h*I_c,Indices(0,D-1,6));
  CN.setElement(-h*I_c,Indices(0,0,5));
  CN.setElement(h*I_c,Indices(0,0,6));
  // Dissipation terms (only in C1 and CN!)
  C1.setElement(2.*GammaL1*ONE_c,Indices(0,D-1,11));
  C1.setElement(-GammaL1*ONE_c,Indices(0,D-1,7));
  C1.setElement(-GammaL1*ONE_c,Indices(0,D-1,8));
  CN.setElement(2.*GammaR3*ONE_c,Indices(0,0,11));
  CN.setElement(-GammaR3*ONE_c,Indices(0,0,7));
  CN.setElement(-GammaR3*ONE_c,Indices(0,0,8));

  C1.setElement(2.*GammaL2*ONE_c,Indices(0,D-1,12));
  C1.setElement(-GammaL2*ONE_c,Indices(0,D-1,9));
  C1.setElement(-GammaL2*ONE_c,Indices(0,D-1,10));
  CN.setElement(2.*GammaR4*ONE_c,Indices(0,0,12));
  CN.setElement(-GammaR4*ONE_c,Indices(0,0,9));
  CN.setElement(-GammaR4*ONE_c,Indices(0,0,10));

  // Two body terms only from H but on both sides
  // sigmax sigmax
  C.setElement(-.5*(1+gamma)*I_c,Indices(0,1,1)); 
  C1.setElement(-.5*(1+gamma)*I_c,Indices(0,1,1)); 
  C.setElement(ONE_c,Indices(1,D-1,1)); 
  CN.setElement(ONE_c,Indices(1,0,1)); 
  // taux taux
  C.setElement(.5*(1+gamma)*I_c,Indices(0,2,2)); 
  C1.setElement(.5*(1+gamma)*I_c,Indices(0,2,2)); 
  C.setElement(ONE_c,Indices(2,D-1,2)); 
  CN.setElement(ONE_c,Indices(2,0,2)); 

  // sigmay sigmay
  C.setElement(-.5*(1-gamma)*I_c,Indices(0,3,3)); 
  C1.setElement(-.5*(1-gamma)*I_c,Indices(0,3,3)); 
  C.setElement(ONE_c,Indices(3,D-1,3)); 
  CN.setElement(ONE_c,Indices(3,0,3)); 
  // tauy tauy
  C.setElement(.5*(1-gamma)*I_c,Indices(0,4,4)); 
  C1.setElement(.5*(1-gamma)*I_c,Indices(0,4,4)); 
  C.setElement(ONE_c,Indices(4,D-1,4)); 
  CN.setElement(ONE_c,Indices(4,0,4)); 

  // Now I reshape, contract with Z and set the indices in proper order
  C.reshape(Indices(D*D,nrOps));
  C1.reshape(Indices(1*D,nrOps));
  CN.reshape(Indices(D*1,nrOps));
  Z.reshape(Indices(nrOps,d*d*d*d));
  C.multiplyRight(Z);C.reshape(Indices(D,D,d*d,d*d));C.permute(Indices(3,1,4,2));
  C1.multiplyRight(Z);C1.reshape(Indices(1,D,d*d,d*d));C1.permute(Indices(3,1,4,2));
  CN.multiplyRight(Z);CN.reshape(Indices(D,1,d*d,d*d));CN.permute(Indices(3,1,4,2));

  double factor=scale!=1.?pow(scale,1/((double)N)):1.; // scale factor
  // Construct the three operators
  Operator* Op1=new Operator(factor*C1);
  Operator* Op=new Operator(factor*C);
  Operator* OpN=new Operator(factor*CN);

  // And set them to the adequate positions in the MPO!
  Lmpo.setOp(0,Op1,true); // Lmpo is now resposible of this pointer
  Lmpo.setOp(N-1,OpN,true);
  if(N>2){
    Lmpo.setOp(1,Op,true);
    for(int k=2;k<N-1;k++){
      // For all the remaining middle sites, I do not need to create
      // new Operators, but just point to the same one!
      Lmpo.setOp(k,Op,false);      
    }
  }
  else{
    // delete the extra Op, to avoid memory leaks, in this very silly case of only two sites
    delete Op;
  }

  // to avoid repeating all this, set a flag that indicates that L is already prepared
  isInit=true;
}


void LiouvillianXYEdges::initZ(mwArray& Z){
  // Basic pieces
  mwArray sig0=identityMatrix(d);//identity
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),datax);//sigmax
  complex_t datay[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigY(Indices(d,d),datay);//sigmay
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz

  complex_t dataplus[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigPlus(Indices(d,d),dataplus);//sigma plus
  complex_t dataminus[]={ZERO_c,ONE_c,ZERO_c,ZERO_c};
  mwArray sigMinus(Indices(d,d),dataminus);//sigma minus

  // Because my sites are double, the actual individual operators will
  // be tensor products of these, or of these with the identity So, I
  // construct the basic products one by one.  (Notice that not all of
  // them are independent, and I could also do a more efficient
  // contraction, but this cost is negligible compared to anything
  // else, as I only need to do this once when initializing the
  // problem, so I opt for the clearerst construction here).

  // *** First of all, we need the identity on both
  mwArray term0=identityMatrix(d*d);
  // *** Now the ones appearing in single-body terms
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
  // (7) sigma_plus sigma_minus (x) Id
  mwArray term7;
  constructOperatorProduct(term7,sigPlus*sigMinus,sig0);
  // (8) Id (x) sigma_plus sigma_minus
  mwArray term8;
  constructOperatorProduct(term8,sig0,sigPlus*sigMinus);
  // (9) sigma_minus sigma_plus (x) Id
  mwArray term9;
  constructOperatorProduct(term9,sigMinus*sigPlus,sig0);
  // (10) Id (x) sigma_minus sigma_plus
  mwArray term10;
  constructOperatorProduct(term10,sig0,sigMinus*sigPlus);
  // (11) sigma_minus  (x) sigma_minus
  mwArray term11;
  constructOperatorProduct(term11,sigMinus,sigMinus);
  // (12) sigma_plus  (x) sigma_plus
  mwArray term12;
  constructOperatorProduct(term12,sigPlus,sigPlus);

  // Fill in the operators in the Z array, to be contracted with C in
  // the construction of the MPO. The order is as listed above.
  int nrOps=13;
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
      Z.setElement(term7.getElement(Indices(d1,d2)),Indices(7,d1,d2));
      Z.setElement(term8.getElement(Indices(d1,d2)),Indices(8,d1,d2));
      Z.setElement(term9.getElement(Indices(d1,d2)),Indices(9,d1,d2));
      Z.setElement(term10.getElement(Indices(d1,d2)),Indices(10,d1,d2));
      Z.setElement(term11.getElement(Indices(d1,d2)),Indices(11,d1,d2));
      Z.setElement(term12.getElement(Indices(d1,d2)),Indices(12,d1,d2));
    }
}

void LiouvillianXYEdges::constructOperatorProduct(mwArray& result,const mwArray& opA,
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

void LiouvillianXYEdges::initTrL(){
  // I can only do this if the Lmpo is already prepared
  initL();

  // I will create copies of the edge and middle operators with the
  // adequate contraction of indices (upi,upj) and complex conjugation
  // and for the rest of the chain, I just point to the operator in
  // the second site, as L itself does

  mwArray auxI=identityMatrix(d);
  auxI.reshape(Indices(1,d*d));

  // First of the chain
  mwArray op=Lmpo.getOp(0).getFullData();
  Indices dimOp=op.getDimensions();
  op.reshape(Indices(dimOp[0],-1));
  op.multiplyLeft(auxI);
  op.reshape(Indices(dimOp[1],d,d,dimOp[3]));
  op.permute(Indices(3,1,2,4));
  trLmpo.setOp(0,new Operator(op),true);

  // Last of the chain
  op=Lmpo.getOp(N-1).getFullData();
  dimOp=op.getDimensions();
  op.reshape(Indices(dimOp[0],-1));
  op.multiplyLeft(auxI);
  op.reshape(Indices(dimOp[1],d,d,dimOp[3]));
  op.permute(Indices(3,1,2,4));
  trLmpo.setOp(N-1,new Operator(op),true);

  if(N>2){
    op=Lmpo.getOp(1).getFullData();
    dimOp=op.getDimensions();
    op.reshape(Indices(dimOp[0],-1));
    op.multiplyLeft(auxI);
    op.reshape(Indices(dimOp[1],d,d,dimOp[3]));
    op.permute(Indices(3,1,2,4));
    trLmpo.setOp(1,new Operator(op),true);
    for(int k=2;k<N-1;k++)
      trLmpo.setOp(k,&trLmpo.getOp(1),false);            
  }
}

void LiouvillianXYEdges::getTrLMPO(int dpur,MPO& result){
  result.initLength(N);
  mwArray auxId=identityMatrix(dpur);
  auxId.reshape(Indices(dpur,1,dpur,1));

  // Set operator 0
  mwArray op=trLmpo.getOp(0).getFullData();
  result.setOp(0,new DoubleOperator(op,auxId),true);
  // Set operator N-1
  op=trLmpo.getOp(N-1).getFullData();
  result.setOp(N-1,new DoubleOperator(op,auxId),true);
  if(N>2){
    op=trLmpo.getOp(1).getFullData();
    result.setOp(1,new DoubleOperator(op,auxId),true);
    for(int k=2;k<N-1;k++)
      result.setOp(k,&result.getOp(1),false);            
  }
}
