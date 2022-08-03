# include "SoftCoulombHamiltonian.h"

#include "Indices.h"

using namespace std;
using namespace shrt;

SoftCoulombHamiltonian::SoftCoulombHamiltonian(int L_,double Delta_,int Ze_,double a_,int n_,
					       double eta_,int Ne_,double penS_):
  L(L_),Delta(Delta_),Ze(Ze_),a(a_),M(n_),eta(eta_),Ne(Ne_),etaS(penS_),d(2),hamil(L_){
  chemicalPotential=eta!=0&&Ne==0;
  penaltyNe=eta!=0&&Ne!=0;
  penaltyS=etaS!=0;
  computeCoefficients(); // first, determine the coefficients and weights of exponentials
  // then initialize the MPO
  initZ();
  //  cout<<"Initialized Z operators"<<endl;
  initHMPO();
}

SoftCoulombHamiltonian::~SoftCoulombHamiltonian(){
  Z.clear();
  lambdas.clear();
  X.clear();
}

void SoftCoulombHamiltonian::computeCoefficients(){
  // 1) Construct the indiidual coefficients of the terms in the true Hamiltonian
  //  mwArray fVals(Indices(L-1,1));
  mwArray fVals(Indices(L,1));
  //  for(int k=1;k<=L-1;k++){
  for(int k=1;k<=L;k++){
    fVals.setElement((1./Delta)/sqrt(k*k+a*a)*ONE_c,Indices(k-1,0));
  }
  //cout<<"fVals="<<fVals<<endl;
  // Now the Hankel matrix F
  mwArray Fmat(Indices(L-M+1,M));
  for(int k=0;k<M;k++) // construct column by column
    for(int j=0;j<L-M+1;j++)
	Fmat.setElement(fVals.getElement(Indices(j+k,0)),Indices(j,k));
  //cout<<"F="<<Fmat<<endl;
  // QR decomposition of F
  mwArray Q,R;
  wrapper::qr(Fmat,Q,R);
  // cout<<"Q="<<Q<<", R="<<R<<endl;
  // Define U1,U2
  mwArray U1(Q);
  U1.resize(Indices(L-M,Q.getDimension(1)));
  //cout<<"U1="<<U1<<endl;
  mwArray U2(U1); // same size
  for(int k=0;k<Q.getDimension(1);k++)
    for(int j=0;j<L-M;j++)
      U2.setElement(Q.getElement(Indices(j+1,k)),Indices(j,k));
  // cout<<"U2="<<U2<<endl;
  // Diagonalize pseudoinv(U1)*U2
  double tol=1E-8;int nr=0;
  // cannot obtain the inverse (non-Hermitian)
  mwArray Wu1,Su1,Vdu1;
  wrapper::svd(U1,tol,nr,Wu1,Su1,Vdu1);
  // int nr=U1.getDimension(0);
  //  mwArray invU1=inverse(U1,nr,tol);
  mwArray invU1=invertDiag(Su1,nr,tol);
  Vdu1.Hconjugate();Wu1.Hconjugate();
  invU1.multiplyLeft(Vdu1);
  invU1.multiplyRight(Wu1);
  invU1.multiplyRight(U2);
  // But this guy is also non-Hermitian and I need to diagonalize it!! 
  mwArray U;
  wrapper::eigNH(invU1,lambdas,U,0,0);
  //cout<<"Obtained eigenvalues "<<lambdas<<endl;
  // Now compute the coefficients x by solving min sum_k|f(k)-sum_i(x_i lambda_i^k)|^2
  // Matrix of coefficients
  mwArray _M(Indices(L,M));
  for(int k=0;k<L;k++)
    for(int j=0;j<M;j++)
      _M.setElement(pow(real(lambdas[j]),k+1)*ONE_c,Indices(k,j));

  // And solve the system
  wrapper::lsd(_M,fVals,X);

  cout<<"Lambda values "<<lambdas<<endl;
  cout<<"Solved the system: found coefficients "<<X<<endl;
}

void SoftCoulombHamiltonian::initZ(){
  int nrOps=12; // should be one less if no penaltyNe and more if
		// penaltyS
  if(penaltyS) nrOps=16;
  Z=mwArray(Indices(nrOps,d*d,d*d));Z.fillWithZero();
  // The basic (single site) operators I need
  mwArray sig0=identityMatrix(d);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP);
  complex_t dataM[]={ZERO_c,ONE_c,ZERO_c,ZERO_c};
  mwArray sigM(Indices(d,d),dataM);
  // Now the "true" operatros, corresponding to pairs of spins
  mwArray S0=identityMatrix(d*d);
  // trick: put them all in a temporary vector
  vector<mwArray> opers(nrOps,S0);
  constructOperatorProduct(opers[1],sigP,sigZ); // S(+z)
  constructOperatorProduct(opers[2],sigM,sigZ); // S(-z)
  constructOperatorProduct(opers[3],sig0,sigP); // S(1+)
  constructOperatorProduct(opers[4],sig0,sigM); // S(1-)
  constructOperatorProduct(opers[5],sigM,sig0); // S(-1)
  constructOperatorProduct(opers[6],sigP,sig0); // S(+1)
  constructOperatorProduct(opers[7],sigZ,sigM); // S(z-)
  constructOperatorProduct(opers[8],sigZ,sigP); // S(z+)
  // For N_i and the square, I need to operate
  mwArray Sz0;constructOperatorProduct(Sz0,sigZ,sig0); // S(zz)
  mwArray S0z;constructOperatorProduct(S0z,sig0,sigZ); // S(zz)
  opers[9]=S0+.5*Sz0+.5*S0z; //N_i
  opers[10]=opers[9]*opers[9]; // N_i^2
  mwArray Szz;constructOperatorProduct(Szz,sigZ,sigZ); // S(zz)
  opers[11]=.25*(S0+S0z+Sz0+Szz);
  if(penaltyS){ // I need also the spin ones
    mwArray Spm;constructOperatorProduct(Spm,sigP,sigM); // S(+-)
    mwArray Smp;constructOperatorProduct(Smp,sigM,sigP); // S(+-)
    opers[12]=-.5*I_c*(Spm-Smp); //Sigma_x
    opers[13]=(-.5)*(Spm+Smp);  //Sigma_y
    opers[14]=.25*(Sz0-S0z); //Sigma_z
    opers[15]=.125*(opers[0]-Szz);
  }
  // Now set them all in place
  for(int i1=0;i1<d*d;i1++)
    for(int i2=0;i2<d*d;i2++){
      for(int k=0;k<nrOps;k++)
	Z.setElement(opers[k].getElement(Indices(i1,i2)),Indices(k,i1,i2));
    }
  Z.reshape(Indices(nrOps,d*d*d*d));
}

void SoftCoulombHamiltonian::constructOperatorProduct(mwArray& result,const mwArray& opA,
						      const mwArray& opB) const{
#ifdef CHECKDIMS
  if(opA.getDimensions()!=opB.getDimensions()){
    cout<<"Error in SoftCoulombHamiltonian::constructOperatorProduct for"
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

void SoftCoulombHamiltonian::initHMPO(){
  int D=penaltyNe?M+7:M+6;
  int DnoS=D; // without the S penalty
  if(penaltyS) D+=3;
  int nrOps=Z.getDimension(0); // Nr of operators (11 or 15)
  //cout<<"initHMPO(); Z has dims "<<Z.getDimensions()<<endl;
  // Now construct the coefficients for every term. I need different
  // ones for the first and last sites, and in the middle, I need to
  // change BY HAND one of the terms before assigning, but most of
  // them are TI.
  mwArray C1(Indices(1,D,nrOps)); // the first site
  mwArray C(Indices(D,D,nrOps)); // the middle site
  mwArray CN(Indices(D,1,nrOps)); // the last site

  // Identity before anything was set
  C.setElement(ONE_c,Indices(0,0,0));
  C1.setElement(ONE_c,Indices(0,0,0));
  // If penaltyNe, then there is an identity term (which I could remove)
  if(penaltyNe)
    CN.setElement(eta*Ne*Ne*ONE_c,Indices(0,0,0));
  // Identity terms after the end
  C.setElement(ONE_c,Indices(D-1,D-1,0));
  CN.setElement(ONE_c,Indices(D-1,0,0));
  // TI Single body term, only if penaltyNe
  if(penaltyNe){
    C.setElement(eta*ONE_c,Indices(0,D-1,10));
    C1.setElement(eta*ONE_c,Indices(0,D-1,10));
    CN.setElement(eta*ONE_c,Indices(0,0,10));
  }
  // And the site dependent for first and last
  complex_t addi=penaltyNe?-2*eta*Ne*ONE_c:eta*ONE_c;
  int mz=L%2!=0?(L-1)/2:L/2;
  complex_t coeff0=(-Ze/(Delta*sqrt(mz*mz+a*a))+1./(Delta*Delta))*ONE_c+addi;
  C1.setElement(coeff0,Indices(0,D-1,9));
  //cout<<"Set local term of first site to "<<-Ze/(Delta*sqrt((mz*mz+a*a)))<<", coeff0="<<coeff0<<endl;
  complex_t coeffN=(-Ze/(Delta*sqrt((L-1-mz)*(L-1-mz)+a*a))+1./(Delta*Delta))*ONE_c+addi;
  CN.setElement(coeffN,Indices(0,0,9));
  //cout<<"Set local term of last site to "<<-(double)Ze/(Delta*sqrt((double)((L-1-mz)*(L-1-mz)+a*a)))<<", coeffN="<<coeffN<<endl;

  // The TI local term that comes form interaction
  C1.setElement(1./(Delta*a)*ONE_c,Indices(0,D-1,11));
  C.setElement(1./(Delta*a)*ONE_c,Indices(0,D-1,11));
  CN.setElement(1./(Delta*a)*ONE_c,Indices(0,0,11));
  
  // Two-body terms with the exponentials
  for(int k=0;k<M;k++){
    complex_t xk=X.getElement(Indices(k,0));
    C.setElement(xk,Indices(0,k+1,9));
    C1.setElement(xk,Indices(0,k+1,9));
    // C.setElement((1./Delta)*xk,Indices(0,k+1,9));
    // C1.setElement((1./Delta)*xk,Indices(0,k+1,9));
    // intermediate terms
    C.setElement(lambdas[k],Indices(k+1,k+1,0));
    // second half
    C.setElement(lambdas[k],Indices(k+1,D-1,9));
    CN.setElement(lambdas[k],Indices(k+1,0,9));
  }
  // kinetic terms
  complex_t coeffT=1./(sqrt(2.)*Delta)*I_c;
  for(int p=1;p<=4;p++){
    C.setElement(coeffT,Indices(0,M+p,p));
    C.setElement(coeffT,Indices(M+p,D-1,p+4));
    C1.setElement(coeffT,Indices(0,M+p,p));
    CN.setElement(coeffT,Indices(M+p,0,p+4));
  }
  // Special two-body term, only of penaltyNe
  if(penaltyNe){
    complex_t sq2eta=eta>=0?sqrt(2.*eta)*ONE_c:sqrt(2.*abs(eta))*I_c;
    C.setElement(sq2eta,Indices(0,M+5,9));
    C.setElement(sq2eta,Indices(M+5,D-1,9));
    C.setElement(ONE_c,Indices(M+5,M+5,0));
    C1.setElement(sq2eta,Indices(0,M+5,9));
    CN.setElement(sq2eta,Indices(M+5,0,9));
  }
  // Special components if there is also a spin penalty
  if(penaltyS){
    // TI local term
    C.setElement(3*etaS*ONE_c,Indices(0,D-1,15));
    C1.setElement(3*etaS*ONE_c,Indices(0,D-1,15));
    CN.setElement(3*etaS*ONE_c,Indices(0,0,15));
    // two body terms
    int Xi=DnoS-1; // the first bond value these terms can use
    complex_t coeffS=etaS>=0?sqrt(2*abs(etaS))*ONE_c:sqrt(2*abs(etaS))*I_c;
    C.setElement(coeffS,Indices(0,Xi,12)); // XX
    C.setElement(ONE_c,Indices(Xi,Xi,0));
    C.setElement(coeffS,Indices(Xi,D-1,12));
    C1.setElement(coeffS,Indices(0,Xi,12));
    CN.setElement(coeffS,Indices(Xi,0,12));
    C.setElement(coeffS,Indices(0,Xi+1,13)); // YY
    C.setElement(ONE_c,Indices(Xi+1,Xi+1,0));
    C.setElement(coeffS,Indices(Xi+1,D-1,13));
    C1.setElement(coeffS,Indices(0,Xi+1,13));
    CN.setElement(coeffS,Indices(Xi+1,0,13));
    C.setElement(coeffS,Indices(0,Xi+2,14)); // ZZ
    C.setElement(ONE_c,Indices(Xi+2,Xi+2,0));
    C.setElement(coeffS,Indices(Xi+2,D-1,14));
    C1.setElement(coeffS,Indices(0,Xi+2,14));
    CN.setElement(coeffS,Indices(Xi+2,0,14));
  }
  // Set the operators on the edges
  C1.reshape(Indices(1*D,nrOps));
  CN.reshape(Indices(D*1,nrOps));
  C1.multiplyRight(Z);C1.reshape(Indices(1,D,d*d,d*d));C1.permute(Indices(3,1,4,2));
  CN.multiplyRight(Z);CN.reshape(Indices(D,1,d*d,d*d));CN.permute(Indices(3,1,4,2));
  Operator* Op1=new Operator(C1);
  Operator* OpN=new Operator(CN);
  hamil.setOp(0,Op1,true); // mpo is now resposible of this pointer
  //cout<<"Set Op(0) to "<<*Op1<<endl;
  hamil.setOp(L-1,OpN,true);
  //cout<<"Set Op("<<L-1<<") to "<<*OpN<<endl;
  // One by one, because of the non-TI part, set the operators in the middle
  for(int k=1;k<L-1;k++){
    // set the site dependent term
    complex_t coeffk=(-Ze/(Delta*sqrt((k-mz)*(k-mz)+a*a))+1./(Delta*Delta))*ONE_c+addi;
    C.setElement(coeffk,Indices(0,D-1,9));
    // reshape and set the operator
    mwArray Caux=reshape(C,Indices(D*D,nrOps));
    Caux.multiplyRight(Z);Caux.reshape(Indices(D,D,d*d,d*d));Caux.permute(Indices(3,1,4,2));
    Operator* Op=new Operator(Caux);
    hamil.setOp(k,Op,true);
    //cout<<"Set Op("<<k<<") to "<<*Op<<endl;
  }
}

void SoftCoulombHamiltonian::getNumberMPO(MPO& mpo) const{
  mpo.initLength(L);
  mwArray _Z(Indices(2,d*d,d*d));_Z.fillWithZero();
  // The basic (single site) operators I need
  mwArray sig0=identityMatrix(d);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray Sz0;constructOperatorProduct(Sz0,sigZ,sig0); // S(zz)
  mwArray S0z;constructOperatorProduct(S0z,sig0,sigZ); // S(zz)
  mwArray S0;constructOperatorProduct(S0,sig0,sig0);
  mwArray S1=S0+.5*Sz0+.5*S0z; //N_i
  for(int i1=0;i1<d*d;i1++)
    for(int i2=0;i2<d*d;i2++){
	_Z.setElement(S0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
	_Z.setElement(S1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
    }
  _Z.reshape(Indices(2,d*d*d*d));
  // the coefficients
  int Dop=2;
  // Center site
  mwArray C(Indices(Dop,Dop,2));C.fillWithZero();
  C.setElement(ONE_c,Indices(0,0,0));
  C.setElement(ONE_c,Indices(0,1,1));
  C.setElement(ONE_c,Indices(1,1,0));
  C.reshape(Indices(Dop*Dop,2));

  // Edges
  mwArray CL(Indices(1,Dop,2));CL.fillWithZero();
  CL.setElement(ONE_c,Indices(0,0,0));
  CL.setElement(ONE_c,Indices(0,1,1));
  CL.reshape(Indices(Dop,2));
  mwArray CR(Indices(Dop,1,2));CR.fillWithZero();
  CR.setElement(ONE_c,Indices(0,0,1));
  CR.setElement(ONE_c,Indices(1,0,0));
  CR.reshape(Indices(Dop,2));

  mwArray Op=C*_Z;
  Op.reshape(Indices(Dop,Dop,d*d,d*d));
  Op.permute(Indices(3,1,4,2));
  mwArray OpL=CL*_Z;
  OpL.reshape(Indices(1,Dop,d*d,d*d));
  OpL.permute(Indices(3,1,4,2));
  mwArray OpR=CR*_Z;
  OpR.reshape(Indices(Dop,1,d*d,d*d));
  OpR.permute(Indices(3,1,4,2));

  // set operators in place
  mpo.setOp(0,new Operator(OpL),true);
  mpo.setOp(L-1,new Operator(OpR),true);
  if(L>2){
    mpo.setOp(1,new Operator(Op),true);
    const Operator* basicOp=&mpo.getOp(1);
    for(int k=2;k<L-1;k++)
      mpo.setOp(k,basicOp);
  }

}

void SoftCoulombHamiltonian::getSpinMPO(MPO& mpo) const{
  mpo.initLength(L);
  int nrOps=5;
  mwArray _Z(Indices(nrOps,d*d,d*d));_Z.fillWithZero();
  // The basic (single site) operators I need
  mwArray sig0=identityMatrix(d);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP);
  complex_t dataM[]={ZERO_c,ONE_c,ZERO_c,ZERO_c};
  mwArray sigM(Indices(d,d),dataM);
  mwArray Sz0;constructOperatorProduct(Sz0,sigZ,sig0); // S(zz)
  mwArray S0z;constructOperatorProduct(S0z,sig0,sigZ); // S(zz)
  mwArray S0;constructOperatorProduct(S0,sig0,sig0);

  vector<mwArray> opers(nrOps,S0);
  mwArray Spm;constructOperatorProduct(Spm,sigP,sigM); // S(+-)
  mwArray Smp;constructOperatorProduct(Smp,sigM,sigP); // S(+-)
  opers[1]=-.5*I_c*(Spm-Smp); //Sigma_x
  opers[2]=-.5*(Spm+Smp);  //Sigma_y
  opers[3]=.25*(Sz0-S0z); //Sigma_z
  mwArray Szz;constructOperatorProduct(Szz,sigZ,sigZ); // S(zz)
  opers[4]=.125*(opers[0]-Szz); // Sigma_(x/y/z)^2
  // Now set them all in place
  for(int i1=0;i1<d*d;i1++)
    for(int i2=0;i2<d*d;i2++){
      for(int k=0;k<nrOps;k++)
	_Z.setElement(opers[k].getElement(Indices(i1,i2)),Indices(k,i1,i2));
    }
  _Z.reshape(Indices(nrOps,d*d*d*d));
  // the coefficients
  int Dop=5;
  // Center site
  mwArray C(Indices(Dop,Dop,nrOps));C.fillWithZero();
  // Edges
  mwArray CL(Indices(1,Dop,nrOps));CL.fillWithZero();
  mwArray CR(Indices(Dop,1,nrOps));CR.fillWithZero();
  // Identity terms before any op
  C.setElement(ONE_c,Indices(0,0,0));
  CL.setElement(ONE_c,Indices(0,0,0));
  // Identity terms after everything
  C.setElement(ONE_c,Indices(Dop-1,Dop-1,0));
  CR.setElement(ONE_c,Indices(Dop-1,0,0));
  // local term
  C.setElement(3*ONE_c,Indices(0,Dop-1,4));
  CL.setElement(3*ONE_c,Indices(0,Dop-1,4));
  CR.setElement(3*ONE_c,Indices(0,0,4));
  // two body terms
  int Xi=1; // the first bond value these terms can use
  C.setElement(sqrt(2)*ONE_c,Indices(0,Xi,1)); // XX
  C.setElement(ONE_c,Indices(Xi,Xi,0));
  C.setElement(sqrt(2)*ONE_c,Indices(Xi,Dop-1,1));
  CL.setElement(sqrt(2)*ONE_c,Indices(0,Xi,1));
  CR.setElement(sqrt(2)*ONE_c,Indices(Xi,0,1));
  C.setElement(sqrt(2)*ONE_c,Indices(0,Xi+1,2)); // YY
  C.setElement(ONE_c,Indices(Xi+1,Xi+1,0));
  C.setElement(sqrt(2)*ONE_c,Indices(Xi+1,Dop-1,2));
  CL.setElement(sqrt(2)*ONE_c,Indices(0,Xi+1,2));
  CR.setElement(sqrt(2)*ONE_c,Indices(Xi+1,0,2));
  C.setElement(sqrt(2)*ONE_c,Indices(0,Xi+2,3)); // ZZ
  C.setElement(ONE_c,Indices(Xi+2,Xi+2,0));
  C.setElement(sqrt(2)*ONE_c,Indices(Xi+2,Dop-1,3));
  CL.setElement(sqrt(2)*ONE_c,Indices(0,Xi+2,3));
  CR.setElement(sqrt(2)*ONE_c,Indices(Xi+2,0,3));

  C.reshape(Indices(Dop*Dop,nrOps));
  CL.reshape(Indices(Dop,nrOps));
  CR.reshape(Indices(Dop,nrOps));

  mwArray Op=C*_Z;
  Op.reshape(Indices(Dop,Dop,d*d,d*d));
  Op.permute(Indices(3,1,4,2));
  mwArray OpL=CL*_Z;
  OpL.reshape(Indices(1,Dop,d*d,d*d));
  OpL.permute(Indices(3,1,4,2));
  mwArray OpR=CR*_Z;
  OpR.reshape(Indices(Dop,1,d*d,d*d));
  OpR.permute(Indices(3,1,4,2));

  // set operators in place
  mpo.setOp(0,new Operator(OpL),true);
  mpo.setOp(L-1,new Operator(OpR),true);
  if(L>2){
    mpo.setOp(1,new Operator(Op),true);
    const Operator* basicOp=&mpo.getOp(1);
    for(int k=2;k<L-1;k++)
      mpo.setOp(k,basicOp);
  }

}

void SoftCoulombHamiltonian::getLocalNumberOp(mwArray& oper) const{
  // The basic (single site) operators I need
  mwArray sig0=identityMatrix(d);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray Sz0;constructOperatorProduct(Sz0,sigZ,sig0); // S(zz)
  mwArray S0z;constructOperatorProduct(S0z,sig0,sigZ); // S(zz)
  mwArray S0;constructOperatorProduct(S0,sig0,sig0);
  oper=S0+.5*Sz0+.5*S0z; //N_i
}
