#include "SchwingerHamiltonian.h"
#include "Indices.h"
#include "Contractor.h"

using namespace std;
using namespace shrt;

/** Empty constructor only for the daughter */
SchwingerHamiltonian::SchwingerHamiltonian():
  hamil(1),N(0),mu(0),x(0),l0(0),alpha(0),nu(0),gweight(0),
  offset(0),Z(){};

SchwingerHamiltonian::SchwingerHamiltonian(int N_,double mu_,double x_,
					   double l0_,double alpha_,
					   double offset_,double nu_,double y_): 
  hamil(N_),N(N_),mu(mu_),x(x_),l0(l0_),alpha(alpha_),nu(nu_),gweight(y_),
  offset(offset_),Z(){
  // Construct the MPO
  initZ();
  initHMPO();
}

void SchwingerHamiltonian::initZ(){
  int d=2; // physical dimension
  // basic spin operators appearing
  mwArray sig0=identityMatrix(d);
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sig3(Indices(d,d),dataz);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig1(Indices(d,d),datax);
  complex_t datay[]={ZERO_c,-1*I_c,I_c,ZERO_c};
  mwArray sig2(Indices(d,d),datay);
  Z=mwArray(Indices(4,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Z.setElement(sig2.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      Z.setElement(sig3.getElement(Indices(i1,i2)),Indices(3,i1,i2));
    }
}

void SchwingerHamiltonian::initHMPO(){
  int D=5; // bond dimension of the MPO
  int d=2; // physical dimension
  double alpha_= alpha + l0;
  // Auxiliary: the values that depend on the site

  // even-odd alpha
  mwArray alphas(Indices(1,N));
  for(int k=0;k<N;k++){
    double signo=(k+1)%2==0?+1.:-1.; // since this k starts in 0
    alphas.setElement((complex_t){alpha_-.25*(1.-signo)},Indices(0,k));
  }

  // sums
  mwArray sumA(Indices(1,N));sumA.fillWithZero();
  complex_t sumTot=ZERO_c;
  for(int k=N-2;k>=0;k--){
    sumTot=sumTot+alphas.getElement(Indices(0,k));
    sumA.setElement(sumTot,Indices(0,k));
  }
  // \TODO To optimize the performance I should only change the C elements
  // when calling the function several times
  for(int k=0;k<N;k++){
    // the MPS containing all the coefficients. But they will depend on the
    // position
    int Dl=D,Dr=D;
    if(k==0) Dl=1;
    if(k==N-1) Dr=1;
    mwArray C(Indices(Dl,Dr,4));C.fillWithZero();
    // identities
    if(k!=N-1) // all but the last one can pass on a D=0
      C.setElement(ONE_c,Indices(0,0,0));

    if(k!=0) // in the first one I canot put 1 and finish
      C.setElement(ONE_c,Indices(5-1,Dr-1,0));

    if((k>0)&&(k<N-2)) // all but the last and first to last can be intermediate
      C.setElement(ONE_c,Indices(3,3,0));

    complex_t alphak=alphas.getElement(Indices(0,k));
    double signok=(k+1)%2==0?+1.:-1.;
    if(k!=N-1){
      C.setElement(gweight*alphak*alphak+(complex_t){gweight*(k+1)/4.+mu/2.+nu/2.+(double)offset/(double)N,0.},
		   Indices(0,Dr-1,0)); // fk' 1
      C.setElement((complex_t){signok*mu/2.+nu/2.,0.}+
		   gweight*sumA.getElement(Indices(0,k)),Indices(0,Dr-1,3));
    }
    else{
      C.setElement((complex_t){mu/2.+nu/2.+(double)offset/(double)N,0.},Indices(0,Dr-1,0)); // fk' 1
      C.setElement((complex_t){signok*mu/2.+nu/2.,0.},Indices(0,Dr-1,3));
    }
    //C(1,Dr,4)=(mu/2)*(-1)^k + sum(alphas(k:N)); // fk Z

    if(k!=(N-1)){
      C.setElement((complex_t){sqrt(x/2),0.},Indices(0,1,1)); // XX
      C.setElement((complex_t){sqrt(x/2),0.},Indices(0,2,2)); // YY
      if(k<N-2)
	C.setElement(gweight*ONE_c,Indices(0,3,3)); // k ZZ
      if(k!=0){
	C.setElement((complex_t){sqrt(x/2),0.},Indices(1,4,1)); // XX
	C.setElement((complex_t){sqrt(x/2),0.},Indices(2,4,2)); // YY
	C.setElement((complex_t){(N-k-1)/2.,0.},Indices(3,4,3)); // N-k ZZ
	//if(k<(N-2)) C.setElement(ONE_c,Indices(3,3,0)); // 
      }
    }
    else{ // in the last site I can only have the second half
	C.setElement((complex_t){sqrt(x/2),0.},Indices(1,0,1)); // XX
	C.setElement((complex_t){sqrt(x/2),0.},Indices(2,0,2)); // YY
	//C.setElement((complex_t){(N-k)/2.,0.},Indices(3,0,3)); // N-k+1 ZZ
    }

    // Now contract MPS and operators and set proper index ordering
    mwArray res=reshape(reshape(C,Indices(Dl*Dr,4))*
			reshape(Z,Indices(4,d*d)),Indices(Dl,Dr,d,d));
    // Create and store an Operator
   hamil.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
  }

}

SchwingerHamiltonian::~SchwingerHamiltonian(){
    hamil.clear();
}

void SchwingerHamiltonian::constructGamma5MPO(MPO& gamma5){
  gamma5.initLength(N);
  int D=4,d=2;
  double val=sqrt(sqrt(x)/N);
  mwArray sig0=identityMatrix(2);
  complex_t data2[]={ZERO_c,ZERO_c,ONE_c,ZERO_c}; // sigma-
  mwArray sig2(Indices(d,d),data2); 
  mwArray sig1=permute(sig2,Indices(2,1)); // sigma+
  mwArray Z(Indices(3,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Z.setElement(sig2.getElement(Indices(i1,i2)),Indices(2,i1,i2));
    } 
  Z.reshape(Indices(3,d*d));

  int Dl=D,Dr=D;
  mwArray C0(Indices(1,Dr,3)); // first site
  mwArray CL(Indices(Dl,1,3)); // last site
  mwArray Ce(Indices(Dl,Dr,3)); // even (middle) sites
  mwArray Co(Indices(Dl,Dr,3)); // odd sites

  // Only identity
  C0.setElement(ONE_c,Indices(0,0,0));
  Ce.setElement(ONE_c,Indices(0,0,0));
  Co.setElement(ONE_c,Indices(0,0,0));
  Ce.setElement(ONE_c,Indices(Dl-1,Dr-1,0));
  Co.setElement(ONE_c,Indices(Dl-1,Dr-1,0));
  CL.setElement(ONE_c,Indices(Dl-1,0,0));
  // Left part
  C0.setElement(val*ONE_c,Indices(0,1,1));
  C0.setElement(val*ONE_c,Indices(0,2,2));
  Ce.setElement(val*ONE_c,Indices(0,1,1));
  Ce.setElement(val*ONE_c,Indices(0,2,2));
  Co.setElement(-val*ONE_c,Indices(0,1,1));
  Co.setElement(-val*ONE_c,Indices(0,2,2));
  // Right part
  Ce.setElement(val*ONE_c,Indices(1,Dr-1,2));
  Ce.setElement(val*ONE_c,Indices(2,Dr-1,1));
  Co.setElement(val*ONE_c,Indices(1,Dr-1,2));
  Co.setElement(val*ONE_c,Indices(2,Dr-1,1));
  CL.setElement(val*ONE_c,Indices(1,0,2));
  CL.setElement(val*ONE_c,Indices(2,0,1));

  C0.reshape(Indices(1*Dr,3));
  Ce.reshape(Indices(Dl*Dr,3));
  Co.reshape(Indices(Dl*Dr,3));
  CL.reshape(Indices(Dl*1,3));

  mwArray res0=C0*Z;res0.reshape(Indices(1,Dr,d,d));res0.permute(Indices(3,1,4,2));
  mwArray resL=CL*Z;resL.reshape(Indices(Dl,1,d,d));resL.permute(Indices(3,1,4,2));
  mwArray rese=Ce*Z;rese.reshape(Indices(Dl,Dr,d,d));rese.permute(Indices(3,1,4,2));
  mwArray reso=Co*Z;reso.reshape(Indices(Dl,Dr,d,d));reso.permute(Indices(3,1,4,2));
  Operator* OpEven=new Operator(rese);
  Operator* OpOdd=new Operator(reso);
  Operator* Op0=new Operator(res0);
  Operator* OpL=new Operator(resL);

  gamma5.setOp(0,Op0,true);
  gamma5.setOp(N-1,OpL,true);
  if(N>2){
    gamma5.setOp(1,OpOdd,true);
    if(N>3)
      gamma5.setOp(2,OpEven,true);
    else delete OpEven;
    for(int k=3;k<N-1;k++)
      if(k%2==0)
	gamma5.setOp(k,OpEven,false);
      else
	gamma5.setOp(k,OpOdd,false);
  }    
  else delete OpOdd; // Not necessary, as N=2, only 0 and last
  }

// void SchwingerHamiltonian::constructVecMPO(MPO& Vmpo) const{
//   //  cout<<" SchwingerHamiltonian::constructVMPO"<<endl;
//   Vmpo.initLength(N);
//   int D=4,d=2;
//   mwArray sig0=identityMatrix(d);
//   complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
//   mwArray sig1(Indices(d,d),datax);
//   complex_t datay[]={ZERO_c,-1*I_c,I_c,ZERO_c};
//   mwArray sig2(Indices(d,d),datay);
//   mwArray Zv=mwArray(Indices(3,d,d));Zv.fillWithZero();
//   for(int i1=0;i1<d;i1++)
//     for(int i2=0;i2<d;i2++){
//       Zv.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
//       Zv.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
//       Zv.setElement(sig2.getElement(Indices(i1,i2)),Indices(2,i1,i2));
//     }
//   Zv.reshape(Indices(3,d*d));
//   // I could just do:
//   //   mwArray Zv=Z.subArray(Indices(3,-1,-1));

//   // Middle site:
//   mwArray C(Indices(D,D,3)),Cl(Indices(1,D,3)),Cr(Indices(D,1,3));
//   C.setElement(ONE_c,Indices(0,0,0));
//   Cl.setElement(ONE_c,Indices(0,0,0));
//   Cr.setElement(ONE_c,Indices(0,0,0));
//   C.setElement(ONE_c,Indices(3,3,0));
//   Cr.setElement(ONE_c,Indices(3,1,0));
//   C.setElement(-1.*sqrt(.5)*ONE_c,Indices(0,1,1));
//   C.setElement(sqrt(.5)*I_c,Indices(0,2,2));
//   Cl.setElement(-1*sqrt(.5)*ONE_c,Indices(0,1,1));
//   Cl.setElement(sqrt(.5)*I_c,Indices(0,2,2));
//   C.setElement(sqrt(.5)*I_c,Indices(1,3,2));
//   C.setElement(sqrt(.5)*ONE_c,Indices(2,3,1));
//   Cr.setElement(sqrt(.5)*I_c,Indices(1,0,2));
//   Cr.setElement(sqrt(.5)*ONE_c,Indices(2,0,1));
//   C.reshape(Indices(D*D,3));
//   Cl.reshape(Indices(D,3));
//   Cr.reshape(Indices(D,3));

//   mwArray Vterm=C*Zv;
//   mwArray VtermL=Cl*Zv;
//   mwArray VtermR=Cr*Zv;
//   Vterm.reshape(Indices(D,D,d,d));Vterm.permute(Indices(3,1,4,2));
//   //cout<<"Coonstructed Vterm, dims "<<Vterm.getDimensions()<<endl;
//   VtermL.reshape(Indices(1,D,d,d));VtermL.permute(Indices(3,1,4,2));
//   //cout<<"Coonstructed VtermL, dims "<<VtermL.getDimensions()<<endl;
//   VtermR.reshape(Indices(D,1,d,d));VtermR.permute(Indices(3,1,4,2));
//   //  cout<<"Coonstructed VtermR, dims "<<VtermR.getDimensions()<<endl;

//   // fill in the MPO
//   // edges: these are the MPO's to destroy
//   Vmpo.setOp(0,new Operator(VtermL),true);
//   Vmpo.setOp(N-1,new Operator(VtermR),true);
//   // the first intermediate one is also its own, the rest are just
//   // copies
//   Vmpo.setOp(1,new Operator(Vterm),true);
//   for(int k=2;k<N-1;k++){
//     Vmpo.setOp(k,&Vmpo.getOp(1),false);
//   }
// }


void SchwingerHamiltonian::constructVMPO(MPO& Vmpo) const{
  //  cout<<" SchwingerHamiltonian::constructVMPO"<<endl;
  Vmpo.initLength(N);
  int D=4,d=2;
  mwArray sig0=identityMatrix(d);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig1(Indices(d,d),datax);
  complex_t datay[]={ZERO_c,-1*I_c,I_c,ZERO_c};
  mwArray sig2(Indices(d,d),datay);
  mwArray Zv=mwArray(Indices(3,d,d));Zv.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Zv.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Zv.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Zv.setElement(sig2.getElement(Indices(i1,i2)),Indices(2,i1,i2));
    }
  Zv.reshape(Indices(3,d*d));
  // I could just do:
  //   mwArray Zv=Z.subArray(Indices(3,-1,-1));

  // Middle site:
  mwArray C(Indices(D,D,3)),Cl(Indices(1,D,3)),Cr(Indices(D,1,3));
  C.setElement(ONE_c,Indices(0,0,0));
  Cl.setElement(ONE_c,Indices(0,0,0));
  //Cr.setElement(.01*ONE_c,Indices(0,0,0)); // REMOVE!!!! This is 1
  C.setElement(ONE_c,Indices(3,3,0));
  //Cr.setElement(ONE_c,Indices(3,1,0)); //MAL! (was (3,0,1)!! )
  Cr.setElement(ONE_c,Indices(3,0,0)); //BIEN!
  C.setElement(sqrt(.5)*ONE_c,Indices(0,1,1));
  C.setElement(sqrt(.5)*ONE_c,Indices(0,2,2));
  Cl.setElement(sqrt(.5)*ONE_c,Indices(0,1,1));
  Cl.setElement(sqrt(.5)*ONE_c,Indices(0,2,2));
  C.setElement(sqrt(.5)*ONE_c,Indices(1,3,1));
  C.setElement(sqrt(.5)*ONE_c,Indices(2,3,2));
  Cr.setElement(sqrt(.5)*ONE_c,Indices(1,0,1));
  Cr.setElement(sqrt(.5)*ONE_c,Indices(2,0,2));
  C.reshape(Indices(D*D,3));
  Cl.reshape(Indices(D,3));
  Cr.reshape(Indices(D,3));

  mwArray Vterm=C*Zv;
  mwArray VtermL=Cl*Zv;
  mwArray VtermR=Cr*Zv;
  Vterm.reshape(Indices(D,D,d,d));Vterm.permute(Indices(3,1,4,2));
  //cout<<"Coonstructed Vterm, dims "<<Vterm.getDimensions()<<endl;
  VtermL.reshape(Indices(1,D,d,d));VtermL.permute(Indices(3,1,4,2));
  //cout<<"Coonstructed VtermL, dims "<<VtermL.getDimensions()<<endl;
  VtermR.reshape(Indices(D,1,d,d));VtermR.permute(Indices(3,1,4,2));
  //  cout<<"Coonstructed VtermR, dims "<<VtermR.getDimensions()<<endl;

  // fill in the MPO
  // edges: these are the MPO's to destroy
  Vmpo.setOp(0,new Operator(VtermL),true);
  Vmpo.setOp(N-1,new Operator(VtermR),true);
  // the first intermediate one is also its own, the rest are just
  // copies
  Vmpo.setOp(1,new Operator(Vterm),true);
  for(int k=2;k<N-1;k++){
    Vmpo.setOp(k,&Vmpo.getOp(1),false);
  }

}

void SchwingerHamiltonian::constructGammaAlphaMPO(MPO& gammaA){
  gammaA.initLength(N);
  int D=2,d=2;
  mwArray sig0=identityMatrix(2);
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigz(Indices(d,d),dataz);

  double constant=N%2==0?l0+alpha-.25:l0+alpha-(N-1.)/(4.*N); // the independent term

  mwArray Z(Indices(2,d,d));Z.fillWithZero();

  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(sigz.getElement(Indices(i1,i2)),Indices(1,i1,i2));
    }

  for(int k=0;k<N;k++){
    int Dl=D,Dr=D;
    if(k==0) Dl=1;
    if(k==N-1) Dr=1;
    mwArray C(Indices(Dl,Dr,2));

    // identities
    if(k!=N-1) // all but the last one can pass on a D=0
      C.setElement(ONE_c,Indices(0,0,0));
    else // the last one carries the constant
      C.setElement(constant*ONE_c,Indices(0,0,0));
    if(k<N-1)
      C.setElement((complex_t){(N-k-1)/(2.*N),0.},Indices(0,Dr-1,1));
    if(k!=0) // everything could have been included before
      C.setElement(ONE_c,Indices(1,Dr-1,0));
    Z.reshape(Indices(2,d*d));
    C.reshape(Indices(Dl*Dr,2));
    // Now contract MPS and operators and set proper index ordering
    mwArray res=C*Z;
    res.reshape(Indices(Dl,Dr,d,d));
    gammaA.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
  }
}

/** Another MPO: this one for staggered polarization */

void SchwingerHamiltonian::constructCondensateMPO(MPO& gammaC){
  gammaC.initLength(N);
  int D=2,d=2;
  double val=sqrt(x)/N;
  mwArray sig0=identityMatrix(2);
  //  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,ZERO_c};
  mwArray sig3(Indices(d,d),dataz);
  mwArray Z=mwArray(Indices(2,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(sig3.getElement(Indices(i1,i2)),Indices(1,i1,i2));
    }
  mwArray C0(Indices(1,D,2));C0.fillWithZero(); // first one
  mwArray Ce(Indices(D,D,2));Ce.fillWithZero(); // even site
  mwArray Co(Indices(D,D,2));Co.fillWithZero(); // odd site
  mwArray CN(Indices(D,1,2));CN.fillWithZero(); // last one
  
  // Nothing, just identity
  C0.setElement(ONE_c,Indices(0,0,0));
  Ce.setElement(ONE_c,Indices(0,0,0));
  Co.setElement(ONE_c,Indices(0,0,0));
  // Operator, with sign
  C0.setElement(val*ONE_c,Indices(0,1,1));
  Ce.setElement(val*ONE_c,Indices(0,1,1));
  Co.setElement((-1*val)*ONE_c,Indices(0,1,1));
  if(N%2==0) // then N-1 odd
    CN.setElement((-1*val)*ONE_c,Indices(0,0,1));
  else
    CN.setElement(val*ONE_c,Indices(0,0,1));    
  // After the operator, just identity
  Ce.setElement(ONE_c,Indices(1,1,0));
  Co.setElement(ONE_c,Indices(1,1,0));
  CN.setElement(ONE_c,Indices(1,0,0));

  // Four operators needed: edges, even and odd
  mwArray res0=reshape(reshape(C0,Indices(D,2))*
			reshape(Z,Indices(2,d*d)),Indices(1,D,d,d));
  mwArray resN=reshape(reshape(CN,Indices(D,2))*
			reshape(Z,Indices(2,d*d)),Indices(D,1,d,d));
  mwArray rese=reshape(reshape(Ce,Indices(D*D,2))*
			reshape(Z,Indices(2,d*d)),Indices(D,D,d,d));
  mwArray reso=reshape(reshape(Co,Indices(D*D,2))*
		       reshape(Z,Indices(2,d*d)),Indices(D,D,d,d));
  Operator* op0=new Operator(permute(res0,Indices(3,1,4,2)));
  Operator* opN=new Operator(permute(resN,Indices(3,1,4,2)));
  Operator* ope=new Operator(permute(rese,Indices(3,1,4,2)));
  Operator* opo=new Operator(permute(reso,Indices(3,1,4,2)));
  // Fill in the MPO
  gammaC.setOp(0,op0,true);
  gammaC.setOp(N-1,opN,true);
  bool savede=0;bool savedo=0;
  for(int k=1;k<N-1;k++){
    if(k%2==0)
      if(savede)// already saved
	gammaC.setOp(k,ope,false);
      else{
	gammaC.setOp(k,ope,true);
	savede=true;
      }
    else //odd site
      if(savedo)// already saved
	gammaC.setOp(k,opo,false);
      else{
	gammaC.setOp(k,opo,true);
	savedo=true;
      }
  }
}


void SchwingerHamiltonian::constructMomentumMPO(MPO& Pmpo) const{
  if(N<=2){
    cout<<"Error: cannot construct the momentum MPO (three body) for a chain of length "<<N<<endl;
    exit(1);
  }
  Pmpo.initLength(N);
  int d=2;
  int D=6; 
  mwArray sig0=identityMatrix(2);
  //  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigz(Indices(d,d),dataz);
  mwArray sigP(Indices(d,d));sigP.setElement(ONE_c,Indices(0,1));
  mwArray sigM(sigP);sigM.Hconjugate();
  mwArray Z(Indices(4,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(sigM.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Z.setElement(sigz.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      Z.setElement(sigP.getElement(Indices(i1,i2)),Indices(3,i1,i2));
    }
  Z.reshape(Indices(4,d*d));
  mwArray C0(Indices(1,D,4));C0.fillWithZero(); // first one
  mwArray C(Indices(D,D,4));C.fillWithZero(); // middle site
  mwArray CN(Indices(D,1,4));CN.fillWithZero(); // last one
  // Set the elements
  C0.setElement(ONE_c,Indices(0,0,0)); // nothing yet  
  C.setElement(ONE_c,Indices(0,0,0)); 
  C.setElement(ONE_c,Indices(D-1,D-1,0));  // all finished
  CN.setElement(ONE_c,Indices(D-1,0,0)); 
  C0.setElement(ONE_c,Indices(0,1,1)); // first term, sigma- sigmaz sigma+ 
  C.setElement(ONE_c,Indices(0,1,1)); 
  C0.setElement(ONE_c,Indices(0,3,3)); // first term, sigma+ sigmaz sigma- 
  C.setElement(ONE_c,Indices(0,3,3)); 
  C.setElement(ONE_c,Indices(2,D-1,3)); // last term, sigma- sigmaz sigma+ 
  CN.setElement(ONE_c,Indices(2,0,3)); 
  C.setElement(ONE_c,Indices(4,D-1,1)); // last term, sigma+ sigmaz sigma-
  CN.setElement(ONE_c,Indices(4,0,1)); 
  C.setElement(-1*I_c,Indices(1,2,2)); // middle term, sigma- sigmaz sigma+ 
  C.setElement(I_c,Indices(3,4,2)); // middle term, sigma+ sigmaz sigma- 

  C0.reshape(Indices(1*D,4));
  C.reshape(Indices(D*D,4));
  CN.reshape(Indices(D*1,4));

  mwArray res0=C0*Z;res0.reshape(Indices(1,D,d,d));res0.permute(Indices(3,1,4,2));
  mwArray res=C*Z;res.reshape(Indices(D,D,d,d));res.permute(Indices(3,1,4,2));
  mwArray resN=CN*Z;resN.reshape(Indices(D,1,d,d));resN.permute(Indices(3,1,4,2));

  Operator* op0=new Operator(res0);
  Operator* op=new Operator(res);
  Operator* opN=new Operator(resN);

  Pmpo.setOp(0,op0,true);
  Pmpo.setOp(N-1,opN,true);
  if(N>2)
    Pmpo.setOp(1,op,true);
  for(int k=2;k<N-1;k++)
    Pmpo.setOp(k,op,false);
}


void SchwingerHamiltonian::constructMomentumSquaredMPO(MPO& Pmpo) const{
  if(N<=2){
    cout<<"Error: cannot construct the momentum squared MPO (three body) for a chain of length "<<N<<endl;
    exit(1);
  }
  Pmpo.initLength(N);
  int d=2;
  int D=6; 
  mwArray sig0=identityMatrix(2);
  //  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigz(Indices(d,d),dataz);
  mwArray sigP(Indices(d,d));sigP.setElement(ONE_c,Indices(0,1));
  mwArray sigM(sigP);sigM.Hconjugate();
  mwArray locS=sig0+sigz;
  mwArray Z(Indices(5,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(sigM.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Z.setElement(sigz.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      Z.setElement(sigP.getElement(Indices(i1,i2)),Indices(3,i1,i2));
      Z.setElement(locS.getElement(Indices(i1,i2)),Indices(4,i1,i2));
    }
  Z.reshape(Indices(5,d*d));
  mwArray C0(Indices(1,D,5));C0.fillWithZero(); // first one
  mwArray C(Indices(D,D,5));C.fillWithZero(); // middle site
  mwArray CN(Indices(D,1,5));CN.fillWithZero(); // last one
  // Set the elements
  C0.setElement(ONE_c,Indices(0,0,0)); // nothing yet  
  C.setElement(ONE_c,Indices(0,0,0)); 
  C.setElement(ONE_c,Indices(D-1,D-1,0));  // all finished
  CN.setElement(ONE_c,Indices(D-1,0,0)); 
  C0.setElement(ONE_c,Indices(0,1,1)); // first term, sigma- sigmaz sigma+ 
  C.setElement(ONE_c,Indices(0,1,1)); 
  C0.setElement(ONE_c,Indices(0,3,3)); // first term, sigma+ sigmaz sigma- 
  C.setElement(ONE_c,Indices(0,3,3)); 
  C.setElement(ONE_c,Indices(2,D-1,3)); // last term, sigma- sigmaz sigma+ 
  CN.setElement(ONE_c,Indices(2,0,3)); 
  C.setElement(ONE_c,Indices(4,D-1,1)); // last term, sigma+ sigmaz sigma-
  CN.setElement(ONE_c,Indices(4,0,1)); 
  C.setElement(x*ONE_c,Indices(1,2,2)); // middle term, sigma- sigmaz sigma+ 
  C.setElement(x*ONE_c,Indices(3,4,2)); // middle term, sigma+ sigmaz sigma- 
  C0.setElement(x*ONE_c,Indices(0,D-1,4)); // single body term 1+sigz
  C.setElement(x*ONE_c,Indices(0,D-1,4)); 
  CN.setElement(x*ONE_c,Indices(0,0,4)); 

  C0.reshape(Indices(1*D,5));
  C.reshape(Indices(D*D,5));
  CN.reshape(Indices(D*1,5));

  mwArray res0=C0*Z;res0.reshape(Indices(1,D,d,d));res0.permute(Indices(3,1,4,2));
  mwArray res=C*Z;res.reshape(Indices(D,D,d,d));res.permute(Indices(3,1,4,2));
  mwArray resN=CN*Z;resN.reshape(Indices(D,1,d,d));resN.permute(Indices(3,1,4,2));

  Operator* op0=new Operator(res0);
  Operator* op=new Operator(res);
  Operator* opN=new Operator(resN);

  Pmpo.setOp(0,op0,true);
  Pmpo.setOp(N-1,opN,true);
  if(N>2)
    Pmpo.setOp(1,op,true);
  for(int k=2;k<N-1;k++)
    Pmpo.setOp(k,op,false);
}



/** Construct a translation operator which is cyclic */

void SchwingerHamiltonian::constructCyclicTranslationMPO(MPO& T,bool right) const{
  T.initLength(N);
  int d=2;
  int D=d*d;  
  // The basic operator on any middle site
  mwArray basicOp=identityMatrix(d*d*d);
  basicOp.reshape(Indices(d,d,d,d,d,d));
  basicOp.permute(Indices(5,1,2,3,4,6));
  basicOp.reshape(Indices(d,d*d,d,d*d));
  // The operator on the left
  mwArray leftOp=identityMatrix(d*d);
  leftOp.reshape(Indices(d,1,d,d*d));
  mwArray rightOp=identityMatrix(d*d);
  rightOp.reshape(Indices(d,d,d,d));
  rightOp.permute(Indices(4,1,2,3));
  rightOp.reshape(Indices(d,d*d,d,1));
  if(!right){ // Take the adjoint, as this is the translation to left
    leftOp.permute(Indices(3,2,1,4),true);
    basicOp.permute(Indices(3,2,1,4),true);
    rightOp.permute(Indices(3,2,1,4),true);
  }
  // Now set the operators in the chain
  T.setOp(0,new Operator(leftOp),true);
  T.setOp(N-1,new Operator(rightOp),true);
  if(N>2)
    T.setOp(1,new Operator(basicOp),true);
  for(int k=2;k<N-1;k++)
    T.setOp(k,&T.getOp(1),false);
}

void SchwingerHamiltonian::constructCyclicTranslationRotateMPO(MPO& T,bool right) const{
  T.initLength(N);
  int d=2;
  int D=d*d;  
  // The basic operator on any middle site
  mwArray basicOp=identityMatrix(d*d*d);
  basicOp.reshape(Indices(d,d,d,d,d,d));
  basicOp.permute(Indices(5,1,2,3,4,6));
  basicOp.reshape(Indices(d,d*d,d,d*d));
  // The operator on the left
  mwArray leftOp=identityMatrix(d*d);
  leftOp.reshape(Indices(d,1,d,d*d));
  mwArray rightOp=identityMatrix(d*d);
  rightOp.reshape(Indices(d,d,d,d));
  rightOp.permute(Indices(4,1,2,3));
  rightOp.reshape(Indices(d,d*d,d,1));
  if(!right){ // Take the adjoint, as this is the translation to left
    leftOp.permute(Indices(3,2,1,4),true);
    basicOp.permute(Indices(3,2,1,4),true);
    rightOp.permute(Indices(3,2,1,4),true);
  }
  // for the odd sites, I apply a sigmax rotation now
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigx(Indices(d,d),datax);
  mwArray basicOpOdd(basicOp);
  basicOpOdd.reshape(Indices(d,D*d*D));
  basicOpOdd.multiplyLeft(sigx);
  basicOpOdd.reshape(Indices(d*D,d,D));
  basicOpOdd.permute(Indices(1,3,2));
  basicOpOdd.reshape(Indices(d*D*D,d));
  basicOpOdd.multiplyRight(sigx);
  basicOpOdd.reshape(Indices(d,D,D,d));
  basicOpOdd.permute(Indices(1,2,4,3));
  // Also to the last one, which must be odd
  rightOp.reshape(Indices(d,D*d));
  rightOp.multiplyLeft(sigx);
  rightOp.reshape(Indices(d*D,d));
  rightOp.multiplyRight(sigx);
  rightOp.reshape(Indices(d,D,d,1));

  // Now set the operators in the chain
  T.setOp(0,new Operator(leftOp),true);
  T.setOp(N-1,new Operator(rightOp),true);
  if(N>2){
    T.setOp(1,new Operator(basicOpOdd),true);
    T.setOp(2,new Operator(basicOp),true);
    for(int k=3;k<N-1;k=k+2){
      T.setOp(k,&T.getOp(1),false);
      T.setOp(k+1,&T.getOp(2),false);
    }
  }
}


void SchwingerHamiltonian::constructShiftRotateMPO(MPO& T,bool right) const{
  T.initLength(N);
  int d=2;
  int D=d;  
  // The basic operator on any middle site
  mwArray basicOp=identityMatrix(d*d);
  basicOp.reshape(Indices(d,d,d,d));
  basicOp.permute(Indices(3,1,2,4));
  // The basic operator for odd sites
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigx(Indices(d,d),datax);
  sigx.reshape(Indices(d*d,1));
  mwArray basicOpOdd(sigx);
  sigx.permute(Indices(2,1));
  basicOpOdd.multiplyRight(sigx);
  basicOpOdd.reshape(Indices(d,d,d,d));
  // The leftmost operator
  mwArray leftOp(Indices(d,1,d,D));
  leftOp.setElement(ONE_c,Indices(0,0,0,0));
  leftOp.setElement(ONE_c,Indices(0,0,1,1));
  // The rightmost operator
  mwArray rightOp(Indices(d,D,d));
  rightOp.setElement(ONE_c,Indices(0,0,0));
  //rightOp.setElement(ONE_c,Indices(1,0,1));
  rightOp.setElement(ONE_c,Indices(0,0,1)); // for ignoring the last one
  rightOp.setElement(ONE_c,Indices(1,1,0));
  rightOp.setElement(ONE_c,Indices(1,1,1));
  // if(!right){
  //   // If going left everything looks different
  //   // exchange left and right edges, before applying the sigmax factors!!
  //   mwArray auxOp(leftOp);
  //   leftOp=rightOp;leftOp.permute(Indices(1,3,2));leftOp.reshape(Indices(d,1,d,D));
  //   rightOp=auxOp;rightOp.permute(Indices(1,4,3,2));
  //   // and THIS rightOp has now the sigmax
  //   // and reorder the legs of the others, too
  //   basicOp.permute(Indices(3,2,1,4));
  //   basicOpOdd.permute(Indices(3,2,1,4));
  // } 
  // and because the last is odd, apply also sigmax
  sigx.reshape(Indices(d,d));
  rightOp.reshape(Indices(d,D*d));
  rightOp.multiplyLeft(sigx);
  rightOp.reshape(Indices(d*D,d));
  rightOp.multiplyRight(sigx);
  rightOp.reshape(Indices(d,D,d,1));

  if(!right){
    // If going left, I do the adjoint operator (could decide to do it differently)
    leftOp.permute(Indices(3,2,1,4),true);
    rightOp.permute(Indices(3,2,1,4),true);
    basicOp.permute(Indices(3,2,1,4),true);
    basicOpOdd.permute(Indices(3,2,1,4),true);
  }

  // And now set the operators in the MPO
  T.setOp(0,new Operator(leftOp),true);
  T.setOp(N-1,new Operator(rightOp),true);
  if(N>2){ // at least 4
    T.setOp(1,new Operator(basicOpOdd),true);
    T.setOp(2,new Operator(basicOp),true);
    for(int k=3;k<N-1;k=k+2){
      T.setOp(k,&T.getOp(1),false);
      T.setOp(k+1,&T.getOp(2),false);
    }
  }
}


/*** Implementation of the attempt for unitaries ***/
void SchwingerHamiltonian::getExponentialMPOx(MPO& expH,complex_t delta,
					      bool even) const{
  cout<<"In SchwingerHamiltonian::getExponentialMPOx with delta="<<delta
      <<" and even="<<even<<endl;
  int d=2;
  mwArray Ol,Or;
  int kfirst=even?0:1;
  int klast=(N-1)-(N-kfirst)%2; // last occupied by the loop
  for(int k=kfirst;k<klast;k+=2){
    if(Ol.isEmpty()){ // I need to compute the Operators the first time
      getTwoBodyTermExponential(Ol,Or,delta);
      expH.setOp(k,new Operator(Ol),true);
      expH.setOp(k+1,new Operator(Or),true);
    }
    else{ // I can just copy the pointers to the ones before
      expH.setOp(k,&expH.getOp(k-2));
      expH.setOp(k+1,&expH.getOp(k-1));
    }
    //    cout<<"Set pos "<<k<<" to Ol "<<Ol.getDimensions()<<" and pos "<<k+1<<" to Or "<<Or.getDimensions()<<endl;
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
  if(klast<N-1){ // sth to be filled at the end
    for(int k=klast+1;k<N;k++){
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

/** The version that produces a single MPO for exp(H_e) and exp(H_o) */
void SchwingerHamiltonian::getExponentialMPOx(MPO& expHx,complex_t delta) const{
  bool imag=(delta.im==0); // imaginary time => delta is a real number (-tau)
  bool transverse=0; // not doing transverse contraction, here
  // I need three operators, for first and last sites, and any one in the middle of the chain
  mwArray op1,opN,opk;
  double delta_=abs(delta);
  constructOperatorUxy(op1,delta_,0,imag,transverse);
  constructOperatorUxy(opN,delta_,N-1,imag,transverse);
  if(N>1)  // otherwise, this is not needed
    constructOperatorUxy(opk,delta_,1,imag,transverse);
  expHx.setOp(0,new Operator(op1),true);
  expHx.setOp(N-1,new Operator(opN),true);
  bool saved=false;
  for(int k=1;k<N-1;k++){
    if(!saved){
      expHx.setOp(k,new Operator(opk),true);
      saved=true; // now I already have a pointer
    }
    else{
      expHx.setOp(k,&expHx.getOp(k-1));
    }      
  }
}

void SchwingerHamiltonian::getExponentialMPOL(MPO& expHL,complex_t delta) const {
  cout<<"SchwingerHamiltonian::getExponentialMPOL, delta="<<delta<<endl;
  // Approximated by Id+delta H_L
  // This is copied from initHMPO, and should probably be elsewhere
  int D=3; // bond dimension of the MPO
  int d=2; // physical dimension
  double alpha_= alpha + l0;
  // Auxiliary: the values that depend on the site

  // even-odd alpha
  mwArray alphas(Indices(1,N));
  for(int k=0;k<N;k++){
    double signo=(k+1)%2==0?+1.:-1.; // since this k starts in 0
    alphas.setElement((complex_t){alpha_-.25*(1.-signo)},Indices(0,k));
  }

  // sums
  mwArray sumA(Indices(1,N));sumA.fillWithZero();
  complex_t sumTot=ZERO_c;
  for(int k=N-2;k>=0;k--){
    sumTot=sumTot+alphas.getElement(Indices(0,k));
    sumA.setElement(sumTot,Indices(0,k));
  }
  // \TODO To optimize the performance I should only change the C elements
  // when calling the function several times
  for(int k=0;k<N;k++){
    // the MPS containing all the coefficients. But they will depend on the
    // position
    int Dl=D,Dr=D;
    if(k==0) Dl=1;
    if(k==N-1) Dr=1;
    mwArray C(Indices(Dl,Dr,4));C.fillWithZero();
    // identities
    // all but the last one can pass on a D=0 (last one corresponding to Id)
    C.setElement(ONE_c,Indices(0,0,0));
    if(k!=0) // in the first one I canot put 1 and finish
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0));
    if((k>0)&&(k<N-2)) // all but the last and first to last can be intermediate
      C.setElement(ONE_c,Indices(1,1,0));
    complex_t alphak=alphas.getElement(Indices(0,k));
    double signok=(k+1)%2==0?+1.:-1.;
    if(k!=N-1){
      complex_t value=alphak*alphak+(complex_t){(k+1)/4.+mu/2.+nu/2.+(double)offset/(double)N,0.};
      C.setElement(delta*value,Indices(0,Dr-1,0)); // fk' 1
      value=(complex_t){signok*mu/2.+nu/2.,0.}+sumA.getElement(Indices(0,k));
      C.setElement(delta*value,Indices(0,Dr-1,3));
    }
    else{
      complex_t value=(complex_t){mu/2.+nu/2.+(double)offset/(double)N,0.};
      C.setElement(delta*value+C.getElement(Indices(0,Dr-1,0)),Indices(0,Dr-1,0)); // fk' 1
      value=(complex_t){signok*mu/2.+nu/2.,0.};
      C.setElement(delta*value,Indices(0,Dr-1,3));
    }
    if(k!=(N-1)){
      if(k<N-2)
	C.setElement(delta*ONE_c,Indices(0,1,3)); // k ZZ
      if(k!=0){
	C.setElement((complex_t){(N-k-1)/2.,0.},Indices(1,Dr-1,3)); // N-k ZZ
      }
    }
    // Now contract MPS and operators and set proper index ordering
    //cout<<"expHL about to setOp("<<k<<") from C "<<C.getDimensions()<<" and Z "<<Z.getDimensions()<<endl;
    mwArray res=reshape(reshape(C,Indices(Dl*Dr,4))*
			reshape(Z,Indices(4,d*d)),Indices(Dl,Dr,d,d));
    // Create and store an Operator
    expHL.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
  }
}

void SchwingerHamiltonian::getMPOHL(MPO& HL) const{
  cout<<"SchwingerHamiltonian::getMPOHL"<<endl;
  // This is copied from getExponentialMPOL, and should probably be elsewhere
  int D=3; // bond dimension of the MPO
  int d=2; // physical dimension
  double alpha_= alpha + l0;
  // Auxiliary: the values that depend on the site

  // even-odd alpha
  mwArray alphas(Indices(1,N));
  for(int k=0;k<N;k++){
    double signo=(k+1)%2==0?+1.:-1.; // since this k starts in 0
    alphas.setElement((complex_t){alpha_-.25*(1.-signo)},Indices(0,k));
  }

  // sums
  mwArray sumA(Indices(1,N));sumA.fillWithZero();
  complex_t sumTot=ZERO_c;
  for(int k=N-2;k>=0;k--){
    sumTot=sumTot+alphas.getElement(Indices(0,k));
    sumA.setElement(sumTot,Indices(0,k));
  }
  for(int k=0;k<N;k++){
    // the MPS containing all the coefficients. But they will depend on the
    // position
    int Dl=D,Dr=D;
    if(k==0) Dl=1;
    if(k==N-1) Dr=1;
    mwArray C(Indices(Dl,Dr,4));C.fillWithZero();
    // identities
    // all but the last one can pass on a D=0 
    if(k!=N-1)
      C.setElement(ONE_c,Indices(0,0,0));
    if(k!=0) // in the first one I canot put 1 and finish
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0));
    if((k>0)&&(k<N-2)) // all but the last and first to last can be intermediate
      C.setElement(ONE_c,Indices(1,1,0));
    complex_t alphak=alphas.getElement(Indices(0,k));
    double signok=(k+1)%2==0?+1.:-1.;
    if(k!=N-1){
      complex_t value=alphak*alphak+(complex_t){(k+1)/4.+mu/2.+nu/2.+(double)offset/(double)N,0.};
      C.setElement(value,Indices(0,Dr-1,0)); // fk' 1
      value=(complex_t){signok*mu/2.+nu/2.,0.}+sumA.getElement(Indices(0,k));
      C.setElement(value,Indices(0,Dr-1,3));
    }
    else{
      complex_t value=(complex_t){mu/2.+nu/2.+(double)offset/(double)N,0.};
      C.setElement(value+C.getElement(Indices(0,Dr-1,0)),Indices(0,Dr-1,0)); // fk' 1
      value=(complex_t){signok*mu/2.+nu/2.,0.};
      C.setElement(value,Indices(0,Dr-1,3));
    }
    if(k!=(N-1)){
      if(k<N-2)
	C.setElement(ONE_c,Indices(0,1,3)); // k ZZ
      if(k!=0){
	C.setElement((complex_t){(N-k-1)/2.,0.},Indices(1,Dr-1,3)); // N-k ZZ
      }
    }
    // Now contract MPS and operators and set proper index ordering
    //cout<<"expHL about to setOp("<<k<<") from C "<<C.getDimensions()<<" and Z "<<Z.getDimensions()<<endl;
    mwArray res=reshape(reshape(C,Indices(Dl*Dr,4))*
			reshape(Z,Indices(4,d*d)),Indices(Dl,Dr,d,d));
    // Create and store an Operator
    HL.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
  }
}

void SchwingerHamiltonian::getBasicUOperators(int pos,double initT,
					      double delta,int orderT,
					      vector<Operator&>& opSet,
					      bool imag,bool transverse){
  if(orderT!=2){
    cout<<"ERROR: SchwingerHamiltonian not supporting other than Trotter "
	<<"order 2 (with a trick!)"<<endl;
    exit(212);
  }
  //opSet.clear();
  //opSet.push_back(oper1);
  //opSet.push_back(oper2);
}

int SchwingerHamiltonian::getUOpersPerStep(int orderT) const{
  if(orderT!=2){
    cout<<"ERROR: SchwingerHamiltonian not supporting other than Trotter "
	<<"order 2 (with a trick!)"<<endl;
    exit(212);
  }
  /* If I use a single MPO to approximate the exponential of H_x as
     even x odd and another MPO for the 1+delta H_L */
  return 2;

  /* If the exponential of H_x is better approximated by second order
     Trotter, then I have three MPOS (even odd even or the other way
     round) for it, plus one for the therm, 1+delta H_L */
  //  return 4;

}

Indices SchwingerHamiltonian::getUOperDimensions(int orderT,int opernr,
						 bool transverse) const{
  if(orderT!=2){
    cout<<"ERROR: SchwingerHamiltonian not supporting other than Trotter "
	<<"order 2 (with a trick!)"<<endl;
    exit(212);
  }
  if(opernr<1||opernr>2){
    cout<<"ERROR: SchwingerHamiltonian has only two operators in the unitary "
	<<"expansion, asking for "<<opernr<<endl;
    exit(212);
  }
  int d=2;
  switch(opernr){
  case 1: // first one is sigx sigx ...
    if(!transverse)
      return Indices(d,4,d,4);
    else
      return Indices(4,d,4,d);
    break;
  case 2: // 1+i delta H_L
    if(!transverse)
      return Indices(d,3,d,3);
    else
      return Indices(3,d,3,d);
    break;
  }
}

void SchwingerHamiltonian::constructDforUxy(mwArray& D,double delta,int pos,bool imag) const{
  int Dl=pos==0?1:4;
  int Dr=pos==N-1?1:4;
  //  int Dr=pos==N-1?4:1;
  bool even=(pos%2==0);
  D.fillWithZero();
  complex_t eps=imag?I_c*x*delta:ONE_c*(-x)*delta;
  // auxiliary values
  complex_t cd=cos(eps);complex_t sd=sin(eps);
  complex_t phPI4=sqrt(2)*.5*(complex_t){1.,1.}; // exp(i pi/4)=sqrt(i)
  if(pos!=0&&pos!=N-1){
    complex_t a=.5*sqrt((ONE_c+cd)*sd)*phPI4;
    complex_t b=.5*sqrt((ONE_c-cd)*sd)*phPI4;
    // D(0)
    D.setElement(.5*(ONE_c+cd),Indices(0,0,0));
    D.setElement(.5*I_c*sd,Indices(1,1,0));
    D.setElement(.5*I_c*sd,Indices(2,2,0));
    D.setElement(.5*(ONE_c-cd),Indices(3,3,0));
    // D(1)
    D.setElement(a,Indices(0,1,1));
    D.setElement(a,Indices(1,0,1));
    if(!even){
      D.setElement(I_c*b,Indices(2,3,1));
      D.setElement(I_c*(-b),Indices(3,2,1));
    }
    else{
      D.setElement(I_c*(-b),Indices(2,3,1));
      D.setElement(I_c*b,Indices(3,2,1));
    }
    // D(2)
    D.setElement(a,Indices(0,2,2));
    D.setElement(a,Indices(2,0,2));
    if(!even){
      D.setElement(-I_c*b,Indices(1,3,2));
      D.setElement(I_c*b,Indices(3,1,2));
    }
    else{
      D.setElement(I_c*b,Indices(1,3,2));
      D.setElement(-I_c*b,Indices(3,1,2));
    }
    // D(3)
    D.setElement(.5*sqrt(ONE_c-cd*cd),Indices(0,3,3));
    D.setElement(.5*sqrt(ONE_c-cd*cd),Indices(3,0,3));
    if(!even){
      D.setElement(-.5*sd,Indices(1,2,3));
      D.setElement(.5*sd,Indices(2,1,3));
    }
    else{
      D.setElement(.5*sd,Indices(1,2,3));
      D.setElement(-.5*sd,Indices(2,1,3));
    }
  }
  else{ // one of the edges
    complex_t a=sqrt(.5*sd)*phPI4;
    // D(0)
    D.setElement(sqrt(.5*(ONE_c+cd)),Indices(0,0,0));
    if(pos==0){
      D.setElement(a,Indices(0,1,1));
      D.setElement(a,Indices(0,2,2));
      D.setElement(sqrt(.5*(ONE_c-cd)),Indices(0,3,3));
    }
    else{
      D.setElement(a,Indices(1,0,1));
      D.setElement(a,Indices(2,0,2));
      D.setElement(sqrt(.5*(ONE_c-cd)),Indices(3,0,3));
    }
  }
}

void SchwingerHamiltonian::constructOperatorUxy(mwArray& res,
						double delta,
						int pos,bool imag,
						bool transverse) const {
  int d=2;
  int Dl=pos==0?1:4;
  int Dr=pos==N-1?4:1;
  mwArray D(Indices(Dl,Dr,4));
  constructDforUxy(D,delta,pos,imag);
  res=reshape(reshape(D,Indices(Dl*Dr,4))*
		      reshape(Z,Indices(4,d*d)),Indices(Dl,Dr,d,d));
  if(transverse) // orientation for transverse MPO
    res.permute(Indices(2,3,1,4));
  else //normal orientation
    res.permute(Indices(3,1,4,2));
}

void SchwingerHamiltonian::getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,
						     complex_t delta) const{
  //cout<<"In  SchwingerHamiltonian::getTwoBodyTermExponential("<<
  //delta<<")"<<endl;
  mwArray H12;
  int d=2;
  computeTwoBodyTerm(H12);
  //cout<<"twobody="<<H12<<endl;
  //H12.multiplyLeft(mwArray(delta));
  // Now take the matrix exponential
  mwArray expH;
  wrapper::expm(H12,expH,delta);
  // Obtained with indices (i'j')(ij) => permute to (i'i)(j'j) for svd
  expH.reshape(Indices(d,d,d,d));
  expH.permute(Indices(1,3,2,4));
  expH.reshape(Indices(d*d,d*d));
  //cout<<"exp(delta H12)="<<expH<<endl;
  // And compute the SVD
  mwArray S; // place for the singular values
  int nr(0);
  double tol=0.;
  wrapper::svd(expH,tol,nr,Ol,S,Or);
  //cout<<"Ol="<<Ol<<"; Or="<<Or<<"; S="<<S<<endl;
  // redistribute the singular values
  S=sqrt(S);
  // TODO: check no NaN???
  Ol.multiplyRight(S);
  Or.multiplyLeft(S);
  // reshape for operators (i'x 1 x i x alpha ))
  Ol.reshape(Indices(d,1,d,-1));
  Or.reshape(Indices(-1,d,d,1));
  Or.permute(Indices(2,1,3,4));  
}


void SchwingerHamiltonian::computeTwoBodyTerm(mwArray& result) 
  const{
  int d=2;
  static bool firstCall(true);
  static mwArray twobody;
  if(firstCall){
    //cout<<"First time computeTwoBodyTerm is called"<<endl;
    complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    const mwArray sigX(Indices(d,d),datax);
    complex_t datay[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
    const mwArray sigY(Indices(d,d),datay);
    twobody=kron(sigX,sigX)+kron(sigY,sigY);
    //cout<<"Computed twobody term "<<twobody<<endl;
    firstCall=false;
  }
  //cout<<"SchwingerHamiltonian::computeTwoBodyTerm with twobody already"<<endl;
  result=(x/2)*twobody;
}


void SchwingerHamiltonian::constructVecMPO(MPO& Vmpo) const{
  cout<<" SchwingerHamiltonian::constructVecMPO"<<endl;
  Vmpo.initLength(N);
  int d=2;
  mwArray sig0=identityMatrix(d);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig1(Indices(d,d),datax);
  complex_t datay[]={ZERO_c,-1*I_c,I_c,ZERO_c};
  mwArray sig2(Indices(d,d),datay);
  mwArray Zv=mwArray(Indices(3,d,d));Zv.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Zv.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Zv.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Zv.setElement(sig2.getElement(Indices(i1,i2)),Indices(2,i1,i2));
    }
  Zv.reshape(Indices(3,d*d));
  // I could just do:
  //   mwArray Zv=Z.subArray(Indices(3,-1,-1));

  // Distinguish even (including 0-th) and odd sites
  int D=4;
  mwArray Cl(Indices(1,D,3)),Cr(Indices(D,1,3));
  mwArray C(Indices(D,D,3));
  C.setElement(ONE_c,Indices(0,0,0));
  Cl.setElement(ONE_c,Indices(0,0,0));
  C.setElement(ONE_c,Indices(D-1,D-1,0));
  Cr.setElement(ONE_c,Indices(D-1,0,0));

  C.setElement(sqrt(.5)*I_c,Indices(0,1,1));
  C.setElement(-1.*sqrt(.5)*I_c,Indices(0,2,2));
  Cl.setElement(sqrt(.5)*I_c,Indices(0,1,1));
  Cl.setElement(-1.*sqrt(.5)*I_c,Indices(0,2,2));

  C.setElement(sqrt(.5)*ONE_c,Indices(1,D-1,2));
  C.setElement(sqrt(.5)*ONE_c,Indices(2,D-1,1));
  Cr.setElement(sqrt(.5)*ONE_c,Indices(1,0,2));
  Cr.setElement(sqrt(.5)*ONE_c,Indices(2,0,1));

  C.reshape(Indices(D*D,3));
  Cl.reshape(Indices(D,3));
  Cr.reshape(Indices(D,3));

  mwArray Vterm=C*Zv;
  mwArray VtermL=Cl*Zv;
  mwArray VtermR=Cr*Zv;
  Vterm.reshape(Indices(D,D,d,d));Vterm.permute(Indices(3,1,4,2));
  //cout<<"Coonstructed Vterm, dims "<<Vterm.getDimensions()<<endl;
  VtermL.reshape(Indices(1,D,d,d));VtermL.permute(Indices(3,1,4,2));
  //cout<<"Coonstructed VtermL, dims "<<VtermL.getDimensions()<<endl;
  VtermR.reshape(Indices(D,1,d,d));VtermR.permute(Indices(3,1,4,2));
  //  cout<<"Coonstructed VtermR, dims "<<VtermR.getDimensions()<<endl;

  // fill in the MPO
  // edges: these are the MPO's to destroy
  Vmpo.setOp(0,new Operator(VtermL),true);
  Vmpo.setOp(N-1,new Operator(VtermR),true);
  // the first intermediate one is also its own, the rest are just
  // copies
  if(N>2){
    Vmpo.setOp(1,new Operator(Vterm),true);
    for(int k=2;k<N-1;k++){
      Vmpo.setOp(k,&Vmpo.getOp(1),false);
    }
  }
}

void SchwingerHamiltonian::constructVecMPOobc(MPO& Vmpo,int k0) const{
  cout<<" SchwingerHamiltonian::constructVecMPOobc"<<endl;
  Vmpo.initLength(N);
  int d=2;
  mwArray sig0=identityMatrix(d);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig1(Indices(d,d),datax);
  complex_t datay[]={ZERO_c,-1*I_c,I_c,ZERO_c};
  mwArray sig2(Indices(d,d),datay);
  mwArray Zv=mwArray(Indices(3,d,d));Zv.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Zv.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Zv.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Zv.setElement(sig2.getElement(Indices(i1,i2)),Indices(2,i1,i2));
    }
  Zv.reshape(Indices(3,d*d));
  // I could just do:
  //   mwArray Zv=Z.subArray(Indices(3,-1,-1));

  // Distinguish even (including 0-th) and odd sites
  int D=4;
  mwArray Cl(Indices(1,D,3)),Cr(Indices(D,1,3));
  Cl.setElement(ONE_c,Indices(0,0,0));
  Cr.setElement(ONE_c,Indices(D-1,0,0));
  Cl.setElement(sqrt(.5)*sin(k0*M_PIl/(N+1))*I_c,Indices(0,1,1));
  Cl.setElement(-1.*sqrt(.5)*sin(k0*M_PIl/(N+1))*I_c,Indices(0,2,2));
  Cr.setElement(sqrt(.5)*ONE_c,Indices(1,0,2));
  Cr.setElement(sqrt(.5)*ONE_c,Indices(2,0,1));
  Cl.reshape(Indices(D,3));
  Cr.reshape(Indices(D,3));
  mwArray VtermL=Cl*Zv;
  mwArray VtermR=Cr*Zv;
  VtermL.reshape(Indices(1,D,d,d));VtermL.permute(Indices(3,1,4,2));
  //cout<<"Coonstructed VtermL, dims "<<VtermL.getDimensions()<<endl;
  VtermR.reshape(Indices(D,1,d,d));VtermR.permute(Indices(3,1,4,2));
  //  cout<<"Coonstructed VtermR, dims "<<VtermR.getDimensions()<<endl;
  // fill in the MPO
  // edges: these are the MPO's to destroy
  Vmpo.setOp(0,new Operator(VtermL),true);
  Vmpo.setOp(N-1,new Operator(VtermR),true);

  mwArray C(Indices(D,D,3));
  C.setElement(ONE_c,Indices(0,0,0));
  C.setElement(ONE_c,Indices(D-1,D-1,0));
  C.setElement(sqrt(.5)*I_c,Indices(0,1,1));
  C.setElement(-1.*sqrt(.5)*I_c,Indices(0,2,2));
  C.setElement(sqrt(.5)*ONE_c,Indices(1,D-1,2));
  C.setElement(sqrt(.5)*ONE_c,Indices(2,D-1,1));

  if(N>2)
    for(int k=1;k<N-1;k++){
      C.setElement(sqrt(.5)*sin(k0*M_PIl*(k+1)/(N+1))*I_c,Indices(0,1,1));
      C.setElement(-1.*sqrt(.5)*sin(k0*M_PIl*(k+1)/(N+1))*I_c,Indices(0,2,2));	       
      mwArray Vterm=reshape(C,Indices(D*D,3))*Zv;
      Vterm.reshape(Indices(D,D,d,d));Vterm.permute(Indices(3,1,4,2));
      Vmpo.setOp(k,new Operator(Vterm),true);
    }
}

void SchwingerHamiltonian::constructVecMPOe(MPO& Vmpo) const{
  cout<<" SchwingerHamiltonian::constructVecMPO"<<endl;
  Vmpo.initLength(N);
  int d=2;
  mwArray sig0=identityMatrix(d);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig1(Indices(d,d),datax);
  complex_t datay[]={ZERO_c,-1*I_c,I_c,ZERO_c};
  mwArray sig2(Indices(d,d),datay);
  mwArray Zv=mwArray(Indices(3,d,d));Zv.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Zv.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Zv.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Zv.setElement(sig2.getElement(Indices(i1,i2)),Indices(2,i1,i2));
    }
  Zv.reshape(Indices(3,d*d));
  // I could just do:
  //   mwArray Zv=Z.subArray(Indices(3,-1,-1));

  // Distinguish even (including 0-th) and odd sites
  int Der=4,Dor=2;
  mwArray Cl(Indices(1,Der,3)),Cr(Indices(Der,1,3));
  mwArray Ce(Indices(Dor,Der,3)),Co(Indices(Der,Dor,3));
  Ce.setElement(ONE_c,Indices(0,0,0));
  Co.setElement(ONE_c,Indices(0,0,0));
  Cl.setElement(ONE_c,Indices(0,0,0));
  Ce.setElement(ONE_c,Indices(Dor-1,Der-1,0));
  Co.setElement(ONE_c,Indices(Der-1,Dor-1,0));
  Cr.setElement(ONE_c,Indices(Der-1,0,0));

  Ce.setElement(sqrt(.5)*I_c,Indices(0,1,1));
  Ce.setElement(-1.*sqrt(.5)*I_c,Indices(0,2,2));
  Cl.setElement(sqrt(.5)*I_c,Indices(0,1,1));
  Cl.setElement(-1.*sqrt(.5)*I_c,Indices(0,2,2));

  Co.setElement(sqrt(.5)*ONE_c,Indices(1,Dor-1,2));
  Co.setElement(sqrt(.5)*ONE_c,Indices(2,Dor-1,1));
  Cr.setElement(sqrt(.5)*ONE_c,Indices(1,0,2));
  Cr.setElement(sqrt(.5)*ONE_c,Indices(2,0,1));

  Ce.reshape(Indices(Dor*Der,3));
  Co.reshape(Indices(Der*Dor,3));
  Cl.reshape(Indices(Der,3));
  Cr.reshape(Indices(Der,3));

  mwArray Vterme=Ce*Zv;
  mwArray Vtermo=Co*Zv;
  mwArray VtermL=Cl*Zv;
  mwArray VtermR=Cr*Zv;
  Vterme.reshape(Indices(Dor,Der,d,d));Vterme.permute(Indices(3,1,4,2));
  Vtermo.reshape(Indices(Der,Dor,d,d));Vtermo.permute(Indices(3,1,4,2));
  //cout<<"Coonstructed Vterm, dims "<<Vterm.getDimensions()<<endl;
  VtermL.reshape(Indices(1,Der,d,d));VtermL.permute(Indices(3,1,4,2));
  //cout<<"Coonstructed VtermL, dims "<<VtermL.getDimensions()<<endl;
  VtermR.reshape(Indices(Der,1,d,d));VtermR.permute(Indices(3,1,4,2));
  //  cout<<"Coonstructed VtermR, dims "<<VtermR.getDimensions()<<endl;

  // fill in the MPO
  // edges: these are the MPO's to destroy
  Vmpo.setOp(0,new Operator(VtermL),true);
  Vmpo.setOp(N-1,new Operator(VtermR),true);
  // the first intermediate one is also its own, the rest are just
  // copies
  if(N>2){
    Vmpo.setOp(1,new Operator(Vtermo),true);
    Vmpo.setOp(2,new Operator(Vterme),true);
    for(int k=3;k<N-1;k=k+2){
      Vmpo.setOp(k,&Vmpo.getOp(1),false);
      Vmpo.setOp(k+1,&Vmpo.getOp(2),false);
    }
  }
}



void SchwingerHamiltonian::constructScalMPO(MPO& Vmpo) const{
    cout<<" SchwingerHamiltonian::constructScalMPO"<<endl;
  Vmpo.initLength(N);
  int d=2;
  mwArray sig0=identityMatrix(d);
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  mwArray sig3(Indices(d,d),dataz);
  mwArray Zv=mwArray(Indices(2,d,d));Zv.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Zv.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Zv.setElement(sig3.getElement(Indices(i1,i2)),Indices(1,i1,i2));
    }
  Zv.reshape(Indices(2,d*d));

  // Distinguish even (including 0-th) and odd sites by the sign
  int D=2;
  mwArray Cl(Indices(1,D,2)),Cr(Indices(D,1,2));
  mwArray Ce(Indices(D,D,2));mwArray Co(Indices(D,D,2));
  Ce.setElement(ONE_c,Indices(0,0,0));Co.setElement(ONE_c,Indices(0,0,0));
  Cl.setElement(ONE_c,Indices(0,0,0));
  Ce.setElement(ONE_c,Indices(D-1,D-1,0));
  Co.setElement(ONE_c,Indices(D-1,D-1,0));
  Cr.setElement(ONE_c,Indices(D-1,0,0));

  Ce.setElement(ONE_c,Indices(0,1,1));
  Cl.setElement(ONE_c,Indices(0,1,1));
  Co.setElement(-1.*ONE_c,Indices(0,1,1));
  Cr.setElement(-1.*ONE_c,Indices(0,0,1));

  Ce.reshape(Indices(D*D,2));
  Co.reshape(Indices(D*D,2));
  Cl.reshape(Indices(D,2));
  Cr.reshape(Indices(D,2));

  mwArray Vterme=Ce*Zv;
  mwArray Vtermo=Co*Zv;
  mwArray VtermL=Cl*Zv;
  mwArray VtermR=Cr*Zv;
  Vterme.reshape(Indices(D,D,d,d));Vterme.permute(Indices(3,1,4,2));
  Vtermo.reshape(Indices(D,D,d,d));Vtermo.permute(Indices(3,1,4,2));
  //cout<<"Coonstructed Vterm, dims "<<Vterm.getDimensions()<<endl;
  VtermL.reshape(Indices(1,D,d,d));VtermL.permute(Indices(3,1,4,2));
  //cout<<"Coonstructed VtermL, dims "<<VtermL.getDimensions()<<endl;
  VtermR.reshape(Indices(D,1,d,d));VtermR.permute(Indices(3,1,4,2));
  //  cout<<"Coonstructed VtermR, dims "<<VtermR.getDimensions()<<endl;

  // fill in the MPO
  // edges: these are the MPO's to destroy
  Vmpo.setOp(0,new Operator(VtermL),true);
  Vmpo.setOp(N-1,new Operator(VtermR),true);
  // the first intermediate one is also its own, the rest are just
  // copies
  if(N>2){
    Vmpo.setOp(1,new Operator(Vtermo),true);
    Vmpo.setOp(2,new Operator(Vterme),true);
    for(int k=3;k<N-1;k=k+2){
      Vmpo.setOp(k,&Vmpo.getOp(1),false);
      Vmpo.setOp(k+1,&Vmpo.getOp(2),false);
    }
  }
}

void SchwingerHamiltonian::constructScalMPOobc(MPO& Vmpo,int k0) const{
    cout<<" SchwingerHamiltonian::constructScalMPO"<<endl;
  Vmpo.initLength(N);
  int d=2;
  mwArray sig0=identityMatrix(d);
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  mwArray sig3(Indices(d,d),dataz);
  mwArray Zv=mwArray(Indices(2,d,d));Zv.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Zv.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Zv.setElement(sig3.getElement(Indices(i1,i2)),Indices(1,i1,i2));
    }
  Zv.reshape(Indices(2,d*d));

  // Distinguish even (including 0-th) and odd sites by the sign
  int D=2;
  mwArray Cl(Indices(1,D,2)),Cr(Indices(D,1,2));
  mwArray Ce(Indices(D,D,2));mwArray Co(Indices(D,D,2));
  Ce.setElement(ONE_c,Indices(0,0,0));Co.setElement(ONE_c,Indices(0,0,0));
  Cl.setElement(ONE_c,Indices(0,0,0));
  Ce.setElement(ONE_c,Indices(D-1,D-1,0));
  Co.setElement(ONE_c,Indices(D-1,D-1,0));
  Cr.setElement(ONE_c,Indices(D-1,0,0));
  Ce.setElement(ONE_c,Indices(0,1,1));
  Cl.setElement(ONE_c*sin(k0*M_PIl/(N+1)),Indices(0,1,1));
  Co.setElement(-1.*ONE_c,Indices(0,1,1));
  Cr.setElement(-1.*sin(k0*M_PIl*(N)/(N+1))*ONE_c,Indices(0,0,1));
  Cl.reshape(Indices(D,2));
  Cr.reshape(Indices(D,2));
  mwArray VtermL=Cl*Zv;
  mwArray VtermR=Cr*Zv;
  VtermL.reshape(Indices(1,D,d,d));VtermL.permute(Indices(3,1,4,2));
  //cout<<"Coonstructed VtermL, dims "<<VtermL.getDimensions()<<endl;
  VtermR.reshape(Indices(D,1,d,d));VtermR.permute(Indices(3,1,4,2));
  //  cout<<"Coonstructed VtermR, dims "<<VtermR.getDimensions()<<endl;

  // fill in the MPO
  // edges: these are the MPO's to destroy
  Vmpo.setOp(0,new Operator(VtermL),true);
  Vmpo.setOp(N-1,new Operator(VtermR),true);

  if(N>2)
    for(int k=1;k<N-1;k=k+2){  
      Co.setElement(-1.*sin(k0*M_PIl*(k+1)/(N+1))*ONE_c,Indices(0,1,1));
      Ce.setElement(ONE_c*sin(k0*M_PIl*(k+2)/(N+1)),Indices(0,1,1));
      Ce.reshape(Indices(D*D,2));
      Co.reshape(Indices(D*D,2));
      mwArray Vterme=Ce*Zv;
      mwArray Vtermo=Co*Zv;
      Vterme.reshape(Indices(D,D,d,d));Vterme.permute(Indices(3,1,4,2));
      Vtermo.reshape(Indices(D,D,d,d));Vtermo.permute(Indices(3,1,4,2));
      //cout<<"Coonstructed Vterm, dims "<<Vterm.getDimensions()<<endl;
      Vmpo.setOp(k,new Operator(Vtermo),true);
      Vmpo.setOp(k+1,new Operator(Vterme),true);
    }
}

void SchwingerHamiltonian::constructVecTwoSiteMPO(MPO& Vmpo) const{
  cout<<" SchwingerHamiltonian::constructVecTwoSiteMPO"<<endl;
  Vmpo.initLength(2);
  int d=2;
  mwArray sig0=identityMatrix(d);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig1(Indices(d,d),datax);
  complex_t datay[]={ZERO_c,-1*I_c,I_c,ZERO_c};
  mwArray sig2(Indices(d,d),datay);
  mwArray Zv=mwArray(Indices(3,d,d));Zv.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Zv.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Zv.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Zv.setElement(sig2.getElement(Indices(i1,i2)),Indices(2,i1,i2));
    }
  Zv.reshape(Indices(3,d*d));
  // I could just do:
  //   mwArray Zv=Z.subArray(Indices(3,-1,-1));

  // Distinguish even (including 0-th) and odd sites
  int Der=4;
  mwArray Cl(Indices(1,Der,3)),Cr(Indices(Der,1,3));
  Cl.setElement(ONE_c,Indices(0,0,0));
  Cr.setElement(ONE_c,Indices(Der-1,0,0));
  Cl.setElement(sqrt(.5)*I_c,Indices(0,1,1));
  Cl.setElement(-1.*sqrt(.5)*I_c,Indices(0,2,2));
  Cr.setElement(sqrt(.5)*ONE_c,Indices(1,0,2));
  Cr.setElement(sqrt(.5)*ONE_c,Indices(2,0,1));
  Cl.reshape(Indices(Der,3));
  Cr.reshape(Indices(Der,3));
  mwArray VtermL=Cl*Zv;
  mwArray VtermR=Cr*Zv;
  VtermL.reshape(Indices(1,Der,d,d));VtermL.permute(Indices(3,1,4,2));
  //cout<<"Coonstructed VtermL, dims "<<VtermL.getDimensions()<<endl;
  VtermR.reshape(Indices(Der,1,d,d));VtermR.permute(Indices(3,1,4,2));
  //  cout<<"Coonstructed VtermR, dims "<<VtermR.getDimensions()<<endl;

  // fill in the MPO
  // edges: these are the MPO's to destroy
  Vmpo.setOp(0,new Operator(VtermL),true);
  Vmpo.setOp(1,new Operator(VtermR),true);
}

void SchwingerHamiltonian::constructScalTwoSiteMPO(MPO& Vmpo) const{
  cout<<" SchwingerHamiltonian::constructScalTwoSiteMPO"<<endl;
  Vmpo.initLength(2);
  int d=2;
  mwArray sig0=identityMatrix(d);
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  mwArray sig3(Indices(d,d),dataz);
  mwArray Zv=mwArray(Indices(2,d,d));Zv.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Zv.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Zv.setElement(sig3.getElement(Indices(i1,i2)),Indices(1,i1,i2));
    }
  Zv.reshape(Indices(2,d*d));

  // Distinguish even (including 0-th) and odd sites by the sign
  int D=2;
  mwArray Cl(Indices(1,D,2)),Cr(Indices(D,1,2));
  Cl.setElement(ONE_c,Indices(0,0,0));
  Cr.setElement(ONE_c,Indices(D-1,0,0));
  Cl.setElement(ONE_c,Indices(0,1,1));
  Cr.setElement(-1.*ONE_c,Indices(0,0,1));
  Cl.reshape(Indices(D,2));
  Cr.reshape(Indices(D,2));
  mwArray VtermL=Cl*Zv;
  mwArray VtermR=Cr*Zv;
  VtermL.reshape(Indices(1,D,d,d));VtermL.permute(Indices(3,1,4,2));
  //cout<<"Coonstructed VtermL, dims "<<VtermL.getDimensions()<<endl;
  VtermR.reshape(Indices(D,1,d,d));VtermR.permute(Indices(3,1,4,2));
  //  cout<<"Coonstructed VtermR, dims "<<VtermR.getDimensions()<<endl;

  // fill in the MPO
  // edges: these are the MPO's to destroy
  Vmpo.setOp(0,new Operator(VtermL),true);
  Vmpo.setOp(1,new Operator(VtermR),true);
}

complex_t SchwingerHamiltonian::computeParity(const MPS& orig) const{
  cout<<"SchwingerHamiltonian::computeParity"<<endl;
  int d=orig.getA(0).getd();
  int L=orig.getLength();
  if(L%2!=0){
    cout<<"ERROR: Length should be even!"<<endl;
    exit(1);
  }
  // First flip, then apply the parity set of MPOs
  MPS aux(orig);
  MPS result(orig);
  for(int k=0;k<orig.getLength()/2;k++){
    // Each pair gets transformed: moved to sym position, but relative 
    // pos does not change, thus the virtual bonds have to be swapped,
    // to connect to the right neighbors.
    mwArray auxE=aux.getA(2*k).getA(); // goes to N-2k-2
    mwArray auxO=aux.getA(2*k+1).getA(); // goes to N-2k-1
    Indices dimsE=auxE.getDimensions();
    Indices dimsO=auxO.getDimensions();
    auxE.reshape(Indices(dimsE[0]*dimsE[1],dimsE[2]));
    auxO.permute(Indices(2,1,3));
    auxO.reshape(Indices(dimsO[1],dimsO[0]*dimsO[2]));
    auxE.multiplyRight(auxO);
    // swap the virtual indices
    auxE.reshape(Indices(dimsE[0],dimsE[1],dimsO[0],dimsO[2]));
    auxE.permute(Indices(1,4,3,2));
    auxE.reshape(Indices(dimsE[0]*dimsO[2],dimsO[0]*dimsE[1]));
    // SVD 
    mwArray U,S,Vdagger;
    wrapper::svd(auxE,U,S,Vdagger);
    //cout<<"Checking SVD: "<<auxE-U*S*Vdagger<<endl;
    auxE=U*sqrt(S);
    auxO=sqrt(S)*Vdagger;
    int Dl=auxE.getDimension(1);
    auxE.reshape(Indices(dimsE[0],dimsO[2],Dl));
    auxO.reshape(Indices(Dl,dimsO[0],dimsE[1]));
    auxO.permute(Indices(2,1,3));
    // now set the tensors in place
    result.replaceSite(N-2-2*k,auxE,0);
    result.replaceSite(N-2-2*k+1,auxO,0);
  }
  Contractor& contractor=Contractor::theContractor();
  cout<<"Flipped sites, now overlap with original is "
      <<contractor.contract(result,orig)<<endl;
  // Finally, apply the set of MPOs that take care of the (-1)^m(m-1)/2
  MPO Pmpo(L);
  constructParityMPO(Pmpo);
  cout<<"Constructed single body MPO "<<endl;
  aux=result;
  contractor.optimize(Pmpo,aux,result);    
  cout<<"Optimized after local terms"<<endl;
  constructSignParityMPO(Pmpo);
  return  contractor.contract(result,Pmpo,orig);

  // //TODO: instead of repeating the construction, I could just swap the operators
  // // (only for pos=0 and N-1 would then be different)
  // for(int k=0;k<L;k++){
  //   aux=result;
  //   constructParityMPO(Pmpo,k);
  //   contractor.optimize(Pmpo,aux,result);    
  // }


  // //result.gaugeCond('R',0);
}

void SchwingerHamiltonian::constructParityMPO(MPO& Pmpo) const{
  Pmpo.initLength(N);
  int d=2;
  double signo=pow(-1,N/2.); // i^L
  complex_t dataE[]={-1.*(signo)*ONE_c,ZERO_c,ZERO_c,ONE_c};
  complex_t dataO[]={signo*ONE_c,ZERO_c,ZERO_c,ONE_c};
  mwArray sigE(Indices(d,d),dataE);
  mwArray sigO(Indices(d,d),dataO);
  sigE.reshape(Indices(d,1,d,1));
  sigO.reshape(Indices(d,1,d,1));
  for(int k=0;k<N/2;k++){
    if(k==0){
      Pmpo.setOp(2*k,new Operator(sigE),true);
      Pmpo.setOp(2*k+1,new Operator(sigO),true);
    }
    else{
      Pmpo.setOp(2*k,&Pmpo.getOp(0),false);
      Pmpo.setOp(2*k+1,&Pmpo.getOp(1),false);
    }
  }
}

void SchwingerHamiltonian::constructSignParityMPO(MPO& Pmpo) const{
  //cout<<"SchwingerHamiltonian::constructSignParityMPO"<<endl;
  Pmpo.initLength(N);
  int d=2;int L=N;
  // edges
  mwArray Cl(Indices(d,1,d,2));
  Cl.setElement(ONE_c,Indices(0,0,0,1));
  Cl.setElement(ONE_c,Indices(1,0,1,0));
  mwArray Cr(Indices(d,2,d,1));
  Cr.setElement(ONE_c,Indices(0,1,0,0));
  Cr.setElement(ONE_c,Indices(1,0,1,0));
  Pmpo.setOp(0,new Operator(Cl),true);
  Pmpo.setOp(N-1,new Operator(Cr),true);
  for(int k=1;k<N/2;k++){
    mwArray C(Indices(d,k+1,d,k+2));
    for(int m=0;m<k+1;m++){
      C.setElement(ONE_c,Indices(0,m,0,m+1));
      C.setElement(ONE_c,Indices(1,m,1,m));
    }
    Pmpo.setOp(k,new Operator(C),true);
    //cout<<"Set term "<<k<<" in sign MPO to "<<C<<endl;
  }
  mwArray Cmid(Indices(d,N/2+1,d,N-N/2));
  for(int m=0;m<N/2+1;m++){
    for(int mp=0;mp<N/2;mp++){
      Cmid.setElement(pow(-1.,(m+mp)*(m+mp+1)/2)*ONE_c,Indices(0,m,0,mp));
      Cmid.setElement(pow(-1.,(m+mp)*(m+mp-1)/2)*ONE_c,Indices(1,m,1,mp));
    }
  }
  Pmpo.setOp(N/2,new Operator(Cmid),true);
  for(int k=N/2+1;k<N-1;k++){
    mwArray C(Indices(d,L-k+1,d,L-k));
    for(int m=0;m<L-k;m++){
      C.setElement(ONE_c,Indices(0,m+1,0,m));
      C.setElement(ONE_c,Indices(1,m,1,m));
    }
    Pmpo.setOp(k,new Operator(C),true);
    //cout<<"Set term "<<k<<" in sign MPO to "<<C<<endl;
  }
}


void SchwingerHamiltonian::constructParityMPO(MPO& Pmpo,int pos) const{
  Pmpo.initLength(N);
  int d=2;
  double signo=pow(-1,N/2.); // i^L
  complex_t dataP[]={I_c,ZERO_c,ZERO_c,ONE_c}; // local terms e^(i pi/2 O)
  complex_t dataPos0[]={signo*I_c,ZERO_c,ZERO_c,ZERO_c};
  complex_t dataPos1[]={ZERO_c,ZERO_c,ZERO_c,ONE_c};
  mwArray sigP(Indices(d,d),dataP);
  mwArray idO=identityMatrix(d);
  mwArray sigPos0(Indices(d,d),dataPos0);
  mwArray sigPos1(Indices(d,d),dataPos1);
  // sigP.reshape(Indices(d,1,d,1));
  // idO.reshape(Indices(d,1,d,1));
  // sigPos0.reshape(Indices(d,1,d,1));
  // sigPos1.reshape(Indices(d,1,d,1));
  mwArray Zv=mwArray(Indices(4,d,d));Zv.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Zv.setElement(idO.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Zv.setElement(sigP.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Zv.setElement(sigPos0.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      Zv.setElement(sigPos1.getElement(Indices(i1,i2)),Indices(3,i1,i2));
    }
  Zv.reshape(Indices(4,d*d));

  // For each position (not edges and not pos)
  mwArray Ck(Indices(2,2,4));Ck.fillWithZero();
  Ck.setElement(ONE_c,Indices(0,0,0));
  Ck.setElement(ONE_c,Indices(1,1,1));
  Ck.reshape(Indices(2*2,4));
  mwArray Cpos(Indices(2,2,4));Cpos.fillWithZero();
  Cpos.setElement(ONE_c,Indices(0,0,3));
  if(pos%2==0)
    Cpos.setElement(ONE_c,Indices(1,1,2));
  else
    Cpos.setElement(-1.*ONE_c,Indices(1,1,2));
  Cpos.reshape(Indices(2*2,4));
  mwArray Cl(Indices(1,2,4));Cl.fillWithZero();
  if(pos!=0){
    Cl.setElement(ONE_c,Indices(0,0,0));  
    Cl.setElement(ONE_c,Indices(0,1,1));  
  }
  else{
    Cl.setElement(ONE_c,Indices(0,0,3));  
    Cl.setElement(ONE_c,Indices(0,1,2));  
  }
  mwArray Cr(Indices(2,1,4));Cr.fillWithZero();
  if(pos!=N-1){
    Cr.setElement(ONE_c,Indices(0,0,0));  
    Cr.setElement(ONE_c,Indices(1,0,1));  
  }
  else{
    Cr.setElement(ONE_c,Indices(0,0,3));  
    Cr.setElement(ONE_c,Indices(1,0,2));  
  }
  Cl.reshape(Indices(2,4));Cr.reshape(Indices(2,4));

  mwArray VtermL=Cl*Zv;
  mwArray VtermR=Cr*Zv;
  VtermL.reshape(Indices(1,2,d,d));VtermL.permute(Indices(3,1,4,2));
  VtermR.reshape(Indices(2,1,d,d));VtermR.permute(Indices(3,1,4,2));
  mwArray Vterm=Ck*Zv;
  Vterm.reshape(Indices(2,2,d,d));Vterm.permute(Indices(3,1,4,2));
  mwArray VtermPos=Cpos*Zv;
  VtermPos.reshape(Indices(2,2,d,d));VtermPos.permute(Indices(3,1,4,2));
  Pmpo.setOp(0,new Operator(VtermL),true);
  Pmpo.setOp(N-1,new Operator(VtermR),true);
  bool used=false;int ref=1; // whether I have already employed Vterm
  for(int k=1;k<N-1;k++){
    if(k==pos)
      Pmpo.setOp(k,new Operator(VtermPos),true);
    else
      if(!used){
	Pmpo.setOp(k,new Operator(Vterm),true);
	used=true;ref=k;
      }
      else
    	Pmpo.setOp(k,&Pmpo.getOp(ref),false);
  }
}


void SchwingerHamiltonian::applyParity(const MPS& orig,MPS& result) const{
  int d=orig.getA(0).getd();
  int L=orig.getLength();
  if(L%2!=0){
    cout<<"ERROR: Length should be even!"<<endl;
    exit(1);
  }
  // First apply the product operator, then flip
  MPO signs(L);
  constructParityMPO(signs);
  MPS aux(orig);
  Contractor& contractor=Contractor::theContractor();
  contractor.optimize(signs,orig,aux);
  result=orig;
  for(int k=0;k<orig.getLength()/2;k++){
    // Each pair gets transformed: moved to sym position, but relative 
    // pos does not change, thus the virtual bonds have to be swapped,
    // to connect to the right neighbors.
    mwArray auxE=aux.getA(2*k).getA(); // goes to N-2k-2
    mwArray auxO=aux.getA(2*k+1).getA(); // goes to N-2k-1
    Indices dimsE=auxE.getDimensions();
    Indices dimsO=auxO.getDimensions();
    auxE.reshape(Indices(dimsE[0]*dimsE[1],dimsE[2]));
    auxO.permute(Indices(2,1,3));
    auxO.reshape(Indices(dimsO[1],dimsO[0]*dimsO[2]));
    auxE.multiplyRight(auxO);
    // swap the virtual indices
    auxE.reshape(Indices(dimsE[0],dimsE[1],dimsO[0],dimsO[2]));
    auxE.permute(Indices(1,4,3,2));
    auxE.reshape(Indices(dimsE[0]*dimsO[2],dimsO[0]*dimsE[1]));
    // SVD 
    mwArray U,S,Vdagger;
    wrapper::svd(auxE,U,S,Vdagger);
    //cout<<"Checking SVD: "<<auxE-U*S*Vdagger<<endl;
    auxE=U*sqrt(S);
    auxO=sqrt(S)*Vdagger;
    int Dl=auxE.getDimension(1);
    auxE.reshape(Indices(dimsE[0],dimsO[2],Dl));
    auxO.reshape(Indices(Dl,dimsO[0],dimsE[1]));
    auxO.permute(Indices(2,1,3));
    // now set the tensors in place
    // the sign, I put in the first one)
    if(k==0)
      result.replaceSite(N-2-2*k,-1.*auxE,0);
    else
      result.replaceSite(N-2-2*k,auxE,0);
    result.replaceSite(N-2-2*k+1,auxO,0);
  }
  //result.gaugeCond('R',0);
}

void SchwingerHamiltonian::applyReflectionRotation(const MPS& orig,MPS& result) const{
  int d=orig.getA(0).getd();
  int L=orig.getLength();
  if(L%2!=0){
    cout<<"ERROR: Length should be even!"<<endl;
    exit(1);
  }
  complex_t data1[]={ZERO_c,ONE_c,ONE_c,ZERO_c}; // sigma_x
  mwArray sigmax(Indices(d,d),data1);

  result=orig;
  for(int k=0;k<orig.getLength();k++){
    // Each Site gets transformed: moved to sym position, N-k-1
    // and multiplied by sigmax
    mwArray auxA=orig.getA(k).getA(); 
    Indices dimsA=auxA.getDimensions();
    auxA.reshape(Indices(dimsA[0],dimsA[1]*dimsA[2]));
    //if(k%2!=0)
    auxA.multiplyLeft(sigmax);
    auxA.reshape(Indices(dimsA[0],dimsA[1],dimsA[2]));
    auxA.permute(Indices(1,3,2)); // swap virtual indices
    result.replaceSite(N-k-1,auxA,0);
  }
  //result.gaugeCond('R',0);
}


// void SchwingerHamiltonian::constructZZMPO(MPO& Vmpo,bool even) const{
//   //  cout<<" SchwingerHamiltonian::constructVMPO"<<endl;
//   Vmpo.initLength(N);
//   int D=3,d=2;
//   mwArray sig0=identityMatrix(d);
//   complex_t data1[]={ONE_c,ZERO_c,ZERO_c,ZERO_c};
//   mwArray sig1(Indices(d,d),data1);
//   complex_t data2[]={ZERO_c,ZERO_c,ZERO_c,ONE_c};
//   mwArray sig2(Indices(d,d),data2);
//   mwArray Zv=mwArray(Indices(3,d,d));Zv.fillWithZero();
//   for(int i1=0;i1<d;i1++)
//     for(int i2=0;i2<d;i2++){
//       Zv.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
//       Zv.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
//       Zv.setElement(sig2.getElement(Indices(i1,i2)),Indices(2,i1,i2));
//     }
//   Zv.reshape(Indices(3,d*d));
//   // Middle sites:
//   mwArray Ce(Indices(D,D,3)),Co(Indices(D,D,3));
//   // edges
//   mwArray Cl(Indices(1,D,3)),Cr(Indices(D,1,3));
//   Ce.setElement(ONE_c,Indices(0,0,0));
//   Co.setElement(ONE_c,Indices(0,0,0));
//   Cl.setElement(ONE_c,Indices(0,0,0));
//   Ce.setElement(ONE_c,Indices(D-1,D-1,0));
//   Co.setElement(ONE_c,Indices(D-1,D-1,0));
//   Cr.setElement(ONE_c,Indices(D-1,0,0));
//   // first component
//   if(even){ // the first member of a pair is even
//     Ce.setElement(ONE_c,Indices(0,1,1));
//     Co.setElement(ONE_c,Indices(1,D-1,2));
//     if(N%2!=0)
//       Cr.setElement(ONE_c,Indices(1,0,2)); // Last is ODD??
//   }
//   else{
//     Co.setElement(ONE_c,Indices(0,1,2));
//     Ce.setElement(ONE_c,Indices(1,D-1,1));
//     Cl.setElement(ONE_c,Indices(0,1,2)); // First is ODD????
//     if(N%2==0)
//       Cr.setElement(ONE_c,Indices(1,0,1)); // Last is EVEN??
//   }
//   Ce.reshape(Indices(D*D,3));
//   Co.reshape(Indices(D*D,3));
//   Cl.reshape(Indices(D,3));
//   Cr.reshape(Indices(D,3));

//   mwArray Vterme=Ce*Zv;
//   mwArray Vtermo=Co*Zv;
//   mwArray VtermL=Cl*Zv;
//   mwArray VtermR=Cr*Zv;
//   Vterme.reshape(Indices(D,D,d,d));Vterme.permute(Indices(3,1,4,2));
//   Vtermo.reshape(Indices(D,D,d,d));Vtermo.permute(Indices(3,1,4,2));
//   VtermL.reshape(Indices(1,D,d,d));VtermL.permute(Indices(3,1,4,2));
//   VtermR.reshape(Indices(D,1,d,d));VtermR.permute(Indices(3,1,4,2));

//   Vmpo.setOp(0,new Operator(VtermL),true);
//   Vmpo.setOp(N-1,new Operator(VtermR),true);
//   // the first intermediate one of each parity is also its own, the rest are just
//   // copies
//   Operator* ope=new Operator(Vterme);
//   Operator* opo=new Operator(Vtermo);
//   bool savede=0;bool savedo=0;
//   for(int k=1;k<N-1;k++){
//     if((k+1)%2==0)
//       if(savede)// already saved
// 	Vmpo.setOp(k,ope,false);
//       else{
// 	Vmpo.setOp(k,ope,true);
// 	savede=true;
//       }
//     else //odd site
//       if(savedo)// already saved
// 	Vmpo.setOp(k,opo,false);
//       else{
// 	Vmpo.setOp(k,opo,true);
// 	savedo=true;
//       }
//   }

// }

// Correlated shift, two 1's to left (hc to right)
void SchwingerHamiltonian::constructShiftPairMPO(MPO& Vmpo) const{
  Vmpo.initLength(N);
  int D=4,d=2,nOp=4;
  mwArray sig0=identityMatrix(d); // identity
  complex_t data1[]={ZERO_c,ONE_c,ZERO_c,ZERO_c}; // sigPlus
  mwArray sig1(Indices(d,d),data1);
  complex_t data2[]={ZERO_c,ZERO_c,ONE_c,ZERO_c}; // sigMinus
  mwArray sig2(Indices(d,d),data2);
  complex_t data3[]={ONE_c,ZERO_c,ZERO_c,ZERO_c}; // (1+sigZ)/2
  mwArray sig3(Indices(d,d),data3);
  mwArray Zv=mwArray(Indices(nOp,d,d));Zv.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Zv.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Zv.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Zv.setElement(sig2.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      Zv.setElement(sig3.getElement(Indices(i1,i2)),Indices(3,i1,i2));
    }
  Zv.reshape(Indices(nOp,d*d));
  // Middle sites:
  mwArray C(Indices(D,D,nOp));
  // edges
  mwArray Cl(Indices(1,D,nOp)),Cr(Indices(D,1,nOp));
  C.setElement(ONE_c,Indices(0,0,0));
  Cl.setElement(ONE_c,Indices(0,0,0));
  C.setElement(ONE_c,Indices(0,1,1));
  Cl.setElement(ONE_c,Indices(0,1,1));
  C.setElement(ONE_c,Indices(1,2,3));
  C.setElement(ONE_c,Indices(2,3,2));
  Cr.setElement(ONE_c,Indices(2,0,2));
  C.setElement(ONE_c,Indices(3,3,0));
  Cr.setElement(ONE_c,Indices(3,0,0));

  C.reshape(Indices(D*D,nOp));
  Cl.reshape(Indices(D,nOp));
  Cr.reshape(Indices(D,nOp));
  mwArray Vterm=C*Zv;
  mwArray VtermL=Cl*Zv;
  mwArray VtermR=Cr*Zv;
  Vterm.reshape(Indices(D,D,d,d));Vterm.permute(Indices(3,1,4,2));
  VtermL.reshape(Indices(1,D,d,d));VtermL.permute(Indices(3,1,4,2));
  VtermR.reshape(Indices(D,1,d,d));VtermR.permute(Indices(3,1,4,2));

  Vmpo.setOp(0,new Operator(VtermL),true);
  Vmpo.setOp(N-1,new Operator(VtermR),true);
  if(N>2){
    Vmpo.setOp(1,new Operator(Vterm),true);
    for(int k=2;k<N-1;k++)
      Vmpo.setOp(k,&Vmpo.getOp(1),false);
  }
}


void SchwingerHamiltonian::constructZZMPO(MPO& Vmpo,bool even) const{
  //  cout<<" SchwingerHamiltonian::constructVMPO"<<endl;
  Vmpo.initLength(N);
  int D=3,d=2;
  mwArray sig0=identityMatrix(d);
  complex_t data1[]={ONE_c,ZERO_c,ZERO_c,ZERO_c};
  mwArray sig1(Indices(d,d),data1);
  complex_t data2[]={ZERO_c,ZERO_c,ZERO_c,ONE_c};
  mwArray sig2(Indices(d,d),data2);
  mwArray Zv=mwArray(Indices(2,d,d));Zv.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Zv.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Zv.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
    }
  Zv.reshape(Indices(2,d*d));
  // Middle sites:
  mwArray C1(Indices(D,D,2)),C2(Indices(D,D,2));
  // edges
  mwArray Cl(Indices(1,D,2)),Cr(Indices(D,1,2));
  C1.setElement(ONE_c,Indices(0,0,0));
  C2.setElement(ONE_c,Indices(0,0,0));
  Cl.setElement(ONE_c,Indices(0,0,0));
  C1.setElement(ONE_c,Indices(D-1,D-1,0));
  C2.setElement(ONE_c,Indices(D-1,D-1,0));
  Cr.setElement(ONE_c,Indices(D-1,0,0));
  // first component
  C1.setElement(ONE_c,Indices(0,1,1));
  if(!even) // start on odd (like 1)
    Cl.setElement(ONE_c,Indices(0,1,1));
  // second component
  C2.setElement(ONE_c,Indices(1,D-1,1));
  if(!even&&(N%2==0)) // only when the second half is in the last
    Cr.setElement(ONE_c,Indices(1,0,1));
  C1.reshape(Indices(D*D,2));
  C2.reshape(Indices(D*D,2));
  Cl.reshape(Indices(D,2));
  Cr.reshape(Indices(D,2));

  mwArray Vterm1=C1*Zv;
  mwArray Vterm2=C2*Zv;
  mwArray VtermL=Cl*Zv;
  mwArray VtermR=Cr*Zv;
  Vterm1.reshape(Indices(D,D,d,d));Vterm1.permute(Indices(3,1,4,2));
  Vterm2.reshape(Indices(D,D,d,d));Vterm2.permute(Indices(3,1,4,2));
  VtermL.reshape(Indices(1,D,d,d));VtermL.permute(Indices(3,1,4,2));
  VtermR.reshape(Indices(D,1,d,d));VtermR.permute(Indices(3,1,4,2));

  Vmpo.setOp(0,new Operator(VtermL),true);
  Vmpo.setOp(N-1,new Operator(VtermR),true);
  // the first intermediate one of each parity is also its own, the rest are just
  // copies
  Operator* ope=even?new Operator(Vterm1):new Operator(Vterm2);
  Operator* opo=even?new Operator(Vterm2):new Operator(Vterm1);
  bool savede=0;bool savedo=0;
  for(int k=1;k<N-1;k++){
    if((k+1)%2==0)
      if(savede)// already saved
	Vmpo.setOp(k,ope,false);
      else{
	Vmpo.setOp(k,ope,true);
	savede=true;
      }
    else //odd site
      if(savedo)// already saved
	Vmpo.setOp(k,opo,false);
      else{
	Vmpo.setOp(k,opo,true);
	savedo=true;
      }
  }

}

void SchwingerHamiltonian::getExactExponentialMPOL(MPO& expHL,complex_t delta,int cutoff,double tol,int lmax,double normFact) const {
  cout<<"SchwingerHamiltonian::getExactExponentialMPOL, delta="<<delta<<" for N="<<N<<endl;
  expHL.initLength(N);
  int d=2; // physical dimension
  double alpha_= alpha + l0;
  /** I will implement the cutoff simply by applying SVD to the resulting MPO tensors! */
  mwArray tmpL(ONE_c);
  if(cutoff!=0){
    if(lmax!=0) cutoff=2*lmax+1;
    cout<<"Applying a cutoff ("<<cutoff<<") to the dimension of the expHL MPO"<<endl;
  }
  // The dimension of the operator will grow site by site, as I move to the right of the chain.
  // So, I need to construct a new operator for every site
  int Dl=1;int Dr=2;
  for(int k=0;k<N-1;k++){
    if(k==N-2) Dr=1;
    mwArray Ck(Indices(d,Dl,d,Dr));
    // set the non-vanishing elements
    for(int xi=0;xi<Dl;xi++){
      int valL0=k%2==0?alpha_-(k+2)/2:alpha_-(k+1)/2; // l corresponding to xip=0 (right index)
      for(int s=0;s<d;s++){
	int xip=s==0?xi+1:xi;
	int lp=xip+valL0; // value of L_n
	if(k==N-2) xip=0; // then I don't care any more
	// mu/2 (1 + (-1)^(k+1) sigma_z)
	double signK=k%2==0?-1.:1.; //(-1)^(k+1)
	int signS=s==0?1:-1; //sig_z=(-1)^(s)
	double massExp=(1.+signK*signS)*.5*mu;
	// Chemical potential term is nu*(1+sigma_z)/2=> for s=0 it is nu, otherwise 0
	double nuExp=s==0?nu:0;
	//int massExp=k%2==0?(s==0?0:1):(s==0?1:0);
	//	complex_t elem=exp(delta*lp*lp)*exp(delta*massExp);
	// Added term of chemical potential. TODO: CHECK!!!!!!
	complex_t elem=exp(delta*lp*lp)*exp(delta*massExp)*exp(delta*nuExp);
	if(lmax!=0&&abs(lp)>lmax) elem=0*elem;
	Ck.setElement(elem,Indices(s,xi,s,xip));
      }
    }
    if(cutoff!=0){
      // cout<<"Cutting off C("<<k<<") which was "<<Ck.getDimensions()<<" with cutoff="
      // 	  <<cutoff<<" and tol="<<tol<<endl;
      // if cutting off, apply left matrix
      Ck.permute(Indices(2,1,3,4));
      Ck.reshape(Indices(Dl,d*d*Dr));
      int newDl=tmpL.getDimension(0);
      Ck.multiplyLeft(tmpL);
      Ck.reshape(Indices(newDl*d*d,Dr));
      mwArray aux(Ck),S;
      int newDr=cutoff;
      wrapper::svd(aux,tol,newDr,Ck,S,tmpL);
      Ck.multiplyRight(S);
      //int newDr=Ck.getDimension(1);
      Ck.reshape(Indices(newDl,d,d,newDr));
      Ck.permute(Indices(2,1,3,4));
    }
    // Another strategy is to cut just the largest values of l, thus not allowing the electric flux to grow beyond something
    //    cout<<"Setting operator "<<k<<" with dims "<<Ck.getDimensions()<<endl;
    expHL.setOp(k,new Operator(normFact*Ck),true);
    Dl=Dr;Dr++;
  }
  mwArray Clast=identityMatrix(d);
  if(mu!=0){
    double signN=N%2==0?1.:-1.; //(-1)^(N)
    Clast.setElement(exp(delta*(1.+signN)*.5*mu),Indices(0,0));
    Clast.setElement(exp(delta*(1.-signN)*.5*mu),Indices(1,1));
  }
  Clast.reshape(Indices(d,1,d,1));
  expHL.setOp(N-1,new Operator(Clast),true);
}

// void SchwingerHamiltonian::getPartExponentialMPOL(MPO& expHL,complex_t delta,
// 						  int second) const{
//   if(second<1||second>N-2){
//     cout<<"Error: getPartExponentialMPOL can only be used with last sigz on"
// 	<<" sites [1,...N-2]"<<endl;
//     exit(212);
//   }
//   complex_t deltaE=(N-second-1.)*.5*delta; // including the factor in front
//   int D=2;
//   int d=2;
//   // basic operators
//   mwArray Ak=identityMatrix(d);
//   mwArray Bk=identityMatrix(d);
//   Ak.setElement(exp(deltaE),Indices(0,0));
//   Ak.setElement(exp(-1.*deltaE),Indices(1,1));
//   Bk.setElement(exp(-1.*deltaE),Indices(0,0));
//   Bk.setElement(exp(deltaE),Indices(1,1));
//   mwArray sigP(Indices(d,d));sigP.setElement(ONE_c,Indices(0,1));
//   mwArray sigM(sigP);sigM.Hconjugate();
//   mwArray sig0=identityMatrix(d);
//   mwArray _Z(Indices(5,d,d));
//   for(int d1=0;d1<d;d1++)
//     for(int d2=0;d2<d;d2++){
//       Z.setElement(sig0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
//       Z.setElement(Ak.getElement(Indices(d1,d2)),Indices(1,d1,d2));
//       Z.setElement(Bk.getElement(Indices(d1,d2)),Indices(2,d1,d2));
//       Z.setElement(sigP.getElement(Indices(d1,d2)),Indices(3,d1,d2));
//       Z.setElement(sigM.getElement(Indices(d1,d2)),Indices(4,d1,d2));
//     }


// }
