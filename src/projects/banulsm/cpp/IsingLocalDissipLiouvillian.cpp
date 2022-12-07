#include "IsingLocalDissipLiouvillian.h"
#include "DoubleOperator.h"

using namespace std;
using namespace shrt;

IsingLocalDissipLiouvillian::IsingLocalDissipLiouvillian(int N_,double V_,double Omega_,double Delta_,
							 double gamma_,double scale_):
  d(2),N(N_),V(V_),Omega(Omega_),Delta(Delta_),gamma(gamma_),scale(scale_){
  initZ();
  initL();
  initLAdjoint();
}

IsingLocalDissipLiouvillian::~IsingLocalDissipLiouvillian(){
  clear();
}

void IsingLocalDissipLiouvillian::clear(){
  Lmpo.clear();
  LAdjmpo.clear();
}

void IsingLocalDissipLiouvillian::initL(){
  cout<<"IsingLocalDissipLiouvillian::initL()"<<endl;
  Lmpo.initLength(N);
  int D=4; // bond dimension of the MPO
  int nrOps=Z.getDimension(0); // Nr of operators (8)
  // Now construct the coefficients for every term. I need different
  // ones for the first and last sites, but the rest are all the same
  for (int k=0;k<N;k++){
    int Dl=k==0?1:D;
    int Dr=k==N-1?1:D;
    if(k>1&&k<N-1){
      // simply point to the already constructed operator
      Lmpo.setOp(k,&Lmpo.getOp(1),false);
    }
    else{
      //cout<<"Creating Operator for site "<<k<<" with size "<<Indices(Dl,Dr,nrOps)<<endl;
      mwArray C(Indices(Dl,Dr,nrOps));
      // set coefficients
      // Identity terms when nothing has happened yet
      if(k<N-1) C.setElement(ONE_c,Indices(0,0,0));
      // Identity terms after the end
      if(k>0) C.setElement(ONE_c,Indices(Dl-1,Dr-1,0));
      // Single body terms (one for each operator, with proper weights)
      C.setElement(-.5*Omega*I_c,Indices(0,Dr-1,1)); //sigx
      C.setElement(.5*Omega*I_c,Indices(0,Dr-1,2)); // taux

      double coeffZ=-.5*(V-Delta);
      if(k==0||k==N-1) coeffZ+=.25*V;
      C.setElement(-coeffZ*I_c,Indices(0,Dr-1,3)); //sigz
      C.setElement(coeffZ*I_c,Indices(0,Dr-1,4)); //tauz

      C.setElement(gamma*ONE_c,Indices(0,Dr-1,5)); // L(x)L^*
      C.setElement(-.5*gamma*ONE_c,Indices(0,Dr-1,6)); // L^dag L
      C.setElement(-.5*gamma*ONE_c,Indices(0,Dr-1,7)); // L^T L^*

      // two-body terms
      if(k<N-1){ // first in the pair
	C.setElement(-.25*V*I_c,Indices(0,1,3));
	C.setElement(.25*V*I_c,Indices(0,2,4));
      }
      if(k>0){ // second in the pair
	C.setElement(ONE_c,Indices(1,Dr-1,3));
	C.setElement(ONE_c,Indices(2,Dr-1,4));
      }
      // Now reshape, combine with Z and set in place

      C.reshape(Indices(Dl*Dr,nrOps));
      Z.reshape(Indices(nrOps,d*d*d*d));
      C.multiplyRight(Z);C.reshape(Indices(Dl,Dr,d*d,d*d));C.permute(Indices(3,1,4,2));

      double factor=scale!=1.?pow(scale,1/((double)N)):1.; // scale factor
      Lmpo.setOp(k,new Operator(factor*C),true); // Lmpo is now resposible of this pointer
    }

  }

}

void IsingLocalDissipLiouvillian::initLAdjoint(){
  cout<<"IsingLocalDissipLiouvillian::initLAdjoint()"<<endl;
  // I can only do this if the Lmpo is already prepared
  initL();
  LAdjmpo.initLength(N);

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

void IsingLocalDissipLiouvillian::initZ(){
  int nrOps=8;
  Z=mwArray(Indices(nrOps,d*d,d*d));

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

  // the double ops used here (called sigma on stm, tau on acilla)

  // *** First of all, we need the identity on both
  mwArray term0=identityMatrix(d*d);
  // *** Now the ones appearing in single-body terms
  // (1) sigX (x) Id
  mwArray term1;
  constructOperatorProduct(term1,sigX,sig0);
  // (2)  Id (x) tauX -> just a permutation of the previous one
  mwArray term2;
  constructOperatorProduct(term2,sig0,sigX);
  // (3) sigZ (x) Id
  mwArray term3;
  constructOperatorProduct(term3,sigZ,sig0);
  // (4)  Id (x) tauZ -> just a permutation of the previous one
  mwArray term4;
  constructOperatorProduct(term4,sig0,sigZ);
  // (5) sigma_plus  (x) tau_plus
  mwArray term5;
  constructOperatorProduct(term5,sigPlus,sigPlus);
  // (6) sigma_minus sigma_plus (x) Id
  mwArray term6;
  constructOperatorProduct(term6,sigMinus*sigPlus,sig0);
  // (7) Id (x) tau_minus tau_plus
  mwArray term7;
  constructOperatorProduct(term7,sig0,sigMinus*sigPlus);

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
    }  
}

