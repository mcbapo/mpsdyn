
#include "CoherentDissipationLiouvillian.h"
#include "DoubleOperator.h"

using namespace std;
using namespace shrt;


CoherentDissipationLiouvillian::CoherentDissipationLiouvillian(int N_,double g_,double gamma_,complex_t epsilon_):
  d(2),N(N_),g(g_),gamma(gamma_),Lmpo(N_),LAdjmpo(N_),trLmpo(N_),darkSubspace(N_),isInit(false),epsilon(epsilon_){
  initL();
  initLAdjoint();
  initTrL();
  initProjectorDS();
}

CoherentDissipationLiouvillian::CoherentDissipationLiouvillian(int N_,double g_,double gamma_,bool priv):
  d(2),N(N_),g(g_),gamma(gamma_),Lmpo(N_),LAdjmpo(N_),trLmpo(N_),darkSubspace(N_),isInit(false){};

CoherentDissipationLiouvillian::~CoherentDissipationLiouvillian(){
  clear();
}

void CoherentDissipationLiouvillian::clear(){
  Lmpo.clear();
  isInit=false;
  LAdjmpo.clear();
  darkSubspace.clear();
}

void CoherentDissipationLiouvillian::initLAdjoint(){
  cout<<"CoherentDissipationLiouvillian::initLAdjoint()"<<endl;
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



void CoherentDissipationLiouvillian::initL(){
  if(isInit) return;
  cout<<"CoherentDissipationLiouvillian::initL()"<<endl;
  int D=6; // bond dimension of the MPO
  mwArray Z; // The operators
  initZ(Z);
  int nrOps=Z.getDimension(0); // Nr of operators (10)
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
  C.setElement(-g*I_c,Indices(0,D-1,1));
  C.setElement(g*I_c,Indices(0,D-1,2));
  C.setElement(-gamma*ONE_c,Indices(0,D-1,3));
  C.setElement(-gamma*ONE_c,Indices(0,D-1,4));
  C.setElement(2*gamma*ONE_c,Indices(0,D-1,5));
  // The same for first site (notice some coefficients are halved)
  C1.setElement(-g*I_c,Indices(0,D-1,1));
  C1.setElement(g*I_c,Indices(0,D-1,2));
  C1.setElement(-.5*gamma*ONE_c,Indices(0,D-1,3));
  C1.setElement(-.5*gamma*ONE_c,Indices(0,D-1,4));
  C1.setElement(gamma*ONE_c,Indices(0,D-1,5));
  // And also for the last site
  CN.setElement(-g*I_c,Indices(0,0,1));
  CN.setElement(g*I_c,Indices(0,0,2));
  CN.setElement(-.5*gamma*ONE_c,Indices(0,0,3));
  CN.setElement(-.5*gamma*ONE_c,Indices(0,0,4));
  CN.setElement(gamma*ONE_c,Indices(0,0,5));
  // And now the two-body terms
  // gamma ( sigma_minus x Id ) (x) ( Id x sigma_minus ) 
  C.setElement(ONE_c,Indices(0,1,6));
  C.setElement(gamma*epsilon*ONE_c,Indices(1,D-1,8));
  C1.setElement(ONE_c,Indices(0,1,6));
  CN.setElement(gamma*epsilon*ONE_c,Indices(1,0,8));
  // gamma ( Id x sigma_minus ) (x) ( sigma_minus x Id ) 
  C.setElement(ONE_c,Indices(0,3,8));
  C.setElement(gamma*epsilon*ONE_c,Indices(3,D-1,6));
  C1.setElement(ONE_c,Indices(0,3,8));
  CN.setElement(gamma*epsilon*ONE_c,Indices(3,0,6));
  // -0.5*gamma ( sigma_plus x Id ) (x) ( sigma_minus x Id ) 
  C.setElement(ONE_c,Indices(0,2,7));
  C.setElement(-.5*epsilon*gamma*ONE_c,Indices(2,D-1,6));
  C1.setElement(ONE_c,Indices(0,2,7));
  CN.setElement(-.5*epsilon*gamma*ONE_c,Indices(2,0,6));
  // -.5*gamma ( sigma_minus x Id ) (x) ( sigma_plus x Id ) 
  C.setElement(-.5*epsilon*gamma*ONE_c,Indices(1,D-1,7));
  CN.setElement(-.5*epsilon*gamma*ONE_c,Indices(1,0,7));
  // -0.5*gamma ( Id x sigma_plus ) (x) ( Id x sigma_minus ) 
  C.setElement(ONE_c,Indices(0,4,9));
  C.setElement(-.5*epsilon*gamma*ONE_c,Indices(4,D-1,8));
  C1.setElement(ONE_c,Indices(0,4,9));
  CN.setElement(-.5*epsilon*gamma*ONE_c,Indices(4,0,8));
  // -0.5*gamma ( Id x sigma_minus ) (x) ( Id x sigma_plus ) 
  C.setElement(-.5*epsilon*gamma*ONE_c,Indices(3,D-1,9));
  CN.setElement(-.5*epsilon*gamma*ONE_c,Indices(3,0,9));

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


void CoherentDissipationLiouvillian::initZ(mwArray& Z){
  // Basic pieces
  mwArray sig0=identityMatrix(d);//identity
  //  complex_t dataz[]={ZERO_c,ZERO_c,ZERO_c,ONE_c}; // sigPlus sigMinus (1-sigmaz)/2
  //mwArray sigZ(Indices(d,d),dataz);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),datax);//sigmax
  complex_t dataplus[]={ZERO_c,ONE_c,ZERO_c,ZERO_c};
  mwArray sigPlus(Indices(d,d),dataplus);//sigma plus
  complex_t dataminus[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigMinus(Indices(d,d),dataminus);//sigma minus
  mwArray sigZ=sigPlus*sigMinus; // Not really sigma_z, but the product of these two

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
  // (1) sigX (x) Id
  mwArray term1;
  constructOperatorProduct(term1,sigX,sig0);
  // (2)  Id (x) sigX -> just a permutation of the previous one
  mwArray term2;
  constructOperatorProduct(term2,sig0,sigX);
  // (3) (1+sigZ)/2 (x) Id
  mwArray term3; 
  constructOperatorProduct(term3,sigZ,sig0);
  // (4)  Id (x) (1+sigZ)/2 -> just a permutation of the previous one
  mwArray term4;
  constructOperatorProduct(term4,sig0,sigZ);
  // (5) sigma_minus (x) sigma_minus
  mwArray term5;
  constructOperatorProduct(term5,sigMinus,sigMinus);
  // *** Now the ones appearing in two-body terms
  // (6) sigma_minus (x) Id
  mwArray term6;
  constructOperatorProduct(term6,sigMinus,sig0);
  // (7) sigma_plus (x) Id
  mwArray term7;
  constructOperatorProduct(term7,sigPlus,sig0);
  // (8) Id (x) sigma_minus
  mwArray term8;
  constructOperatorProduct(term8,sig0,sigMinus);
  // (9) Id (x) sigma_plus
  mwArray term9;
  constructOperatorProduct(term9,sig0,sigPlus);

  // Fill in the operators in the Z array, to be contracted with C in
  // the construction of the MPO. The order is as listed above.
  int nrOps=10;
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
    }
}

void CoherentDissipationLiouvillian::initProjectorDS(){
  // Basic pieces
  complex_t dataP0[]={ONE_c,ZERO_c,ZERO_c,ZERO_c};
  complex_t dataP1[]={ZERO_c,ZERO_c,ZERO_c,ONE_c};
  complex_t dataP01[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  complex_t dataP10[]={ZERO_c,ONE_c,ZERO_c,ZERO_c};
  mwArray sigP0(Indices(d,d),dataP0);
  mwArray sigP1(Indices(d,d),dataP1);
  mwArray sigP01(Indices(d,d),dataP01);
  mwArray sigP10(Indices(d,d),dataP10);
  // Now construct the Z array
  int nrOps=4;
  mwArray Zd(Indices(nrOps,d,d));
  for(int d1=0;d1<d;d1++)
    for(int d2=0;d2<d;d2++){
      Zd.setElement(sigP0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
      Zd.setElement(sigP01.getElement(Indices(d1,d2)),Indices(1,d1,d2));
      Zd.setElement(sigP10.getElement(Indices(d1,d2)),Indices(2,d1,d2));
      Zd.setElement(sigP1.getElement(Indices(d1,d2)),Indices(3,d1,d2));      
    }
  Zd.reshape(Indices(nrOps,d*d));
  // And now the C coefficients
  int D=4; // bond dimension
  for(int k=0;k<N;k++){
    int Dl=k==0?1:D;
    int Dr=k==N-1?1:D;
    complex_t signk=k%2==0?ONE_c:-1.*ONE_c;
    if(epsilon!=ONE_c)
      signk=signk*power(epsilon,1.*k);
    mwArray C(Indices(Dl,Dr,nrOps));
    // Just |0><0|
    C.setElement(ONE_c,Indices(0,0,0)); // before anything
    if(k>0&&k<N-1){
      C.setElement(ONE_c,Indices(1,1,0)); // in the middle
      C.setElement(ONE_c,Indices(2,2,0)); // in the middle
    }
    if(k>0)
      C.setElement(ONE_c,Indices(3,Dr-1,0)); // after everything
    //  |0><1|
    if(k<N-1) // first in the pair
      C.setElement(signk*(1./sqrt(N))*ONE_c,Indices(0,1,1)); // (-1)^k
    if(k>0) // second in the pair
      C.setElement(signk*(1./sqrt(N))*ONE_c,Indices(2,Dr-1,1)); // (-1)^k
    //  |1><0|
    if(k<N-1) // first in the pair
      C.setElement(signk*(1./sqrt(N))*ONE_c,Indices(0,2,2)); // (-1)^k
    if(k>0) // second in the pair
      C.setElement(signk*(1./sqrt(N))*ONE_c,Indices(1,Dr-1,2)); // (-1)^k
    //  |1><1| // both together
    C.setElement(signk*signk*(1./N)*ONE_c,Indices(0,Dr-1,3)); 
    C.reshape(Indices(Dl*Dr,nrOps));
    C.multiplyRight(Zd);C.reshape(Indices(Dl,Dr,d,d));C.permute(Indices(3,1,4,2));

    Operator* Op=new Operator(C);
    darkSubspace.setOp(k,Op,true);
  }
}

void CoherentDissipationLiouvillian::constructOperatorProduct(mwArray& result,const mwArray& opA,
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

void CoherentDissipationLiouvillian::initTrL(){
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

void CoherentDissipationLiouvillian::getTrLMPO(int dpur,MPO& result){
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

