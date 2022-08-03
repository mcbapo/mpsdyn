#include "KondoHamiltonian.h"
#include "Indices.h"

using namespace shrt;
using namespace std;

KondoHamiltonian::KondoHamiltonian(int L_,double J_,double t_):d(2),L(L_),Jpar(J_),Jperp(J_),t(t_),omega(0.),hamil(L+2){
  initZ();
  initHMPO();
}
KondoHamiltonian::KondoHamiltonian(int L_,double Jpar_,double Jperp_,double t_,double omega_):d(2),L(L_),Jpar(Jpar_),Jperp(Jperp_),t(t_),omega(omega_),hamil(L+2){
  initZ();
  initHMPO();
}

KondoHamiltonian::~KondoHamiltonian(){};

// TODO: This should be soemwhere general, and work with any
// dimensions of the opers
void KondoHamiltonian::constructOperatorProduct(mwArray& result,
						const mwArray& opA,
						const mwArray& opB) const{
#ifdef CHECKDIMS
  if(opA.getDimensions()!=opB.getDimensions()){
    cout<<"Error in KondoHamiltonian::constructOperatorProduct for"
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

void KondoHamiltonian::initZ(){
    // Basic pieces
  mwArray sig0=identityMatrix(d);//identity
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP);//sigmaPlus
  mwArray sigM=Hconjugate(sigP); //sigmaMinus
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz

  // The spin site uses "normal" Z for a spin Hamiltonian
  Zsp=mwArray(Indices(4,d,d));Zsp.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Zsp.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Zsp.setElement(sigP.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Zsp.setElement(sigM.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      Zsp.setElement(sigZ.getElement(Indices(i1,i2)),Indices(3,i1,i2));
    }
  Zsp.reshape(Indices(4,d*d));

  // Initialize the operators for fermionic sites: these are like two
  // spins: all double operators
  nrOps=13;
  vector<mwArray> terms(nrOps,mwArray::empty);
  // *** First of all, we need the identity on both
  terms[0]=identityMatrix(d*d);
  // *** Now the ones appearing in fermion hopping terms
  // (1) sigP (x) tauZ
  constructOperatorProduct(terms[1],sigP,sigZ);  
  // (2) sigM (x) tauZ
  constructOperatorProduct(terms[2],sigM,sigZ);  
  // (3) Id (x) tauP
  constructOperatorProduct(terms[3],sig0,sigP);  
  // (4) Id (x) tauM
  constructOperatorProduct(terms[4],sig0,sigM);  
  // (5) sigM (x) Id
  constructOperatorProduct(terms[5],sigM,sig0);
  // (6) sigP (x) Id
  constructOperatorProduct(terms[6],sigP,sig0);
  // (7) sigZ (x) tauM
  constructOperatorProduct(terms[7],sigZ,sigM);
  // (8) sigZ (x) tauP
  constructOperatorProduct(terms[8],sigZ,sigP);

  // *** And the local ones forthe first site
  // (9) F_0^+= sigP (x) tauM
  constructOperatorProduct(terms[9],sigP,sigM);
  // (10) F_0^-= sigM (x) tauP
  constructOperatorProduct(terms[10],sigM,sigP);
  // (11) F_0^z= .5*(1+sigZ) (x) Id - Id (x) (1+tauZ)/2 // actually, defd here without the .5
  mwArray aux;
  constructOperatorProduct(terms[11],sigZ,sig0);
  constructOperatorProduct(aux,sig0,sigZ);
  terms[11]=terms[11]-aux; // simply sigZ-tauZ
  // // (12) sigZ (x) tauZ (for penalty only)
  // constructOperatorProduct(terms[12],sigZ,sigZ);
  terms[12]=terms[11]*terms[11];
  
  // Fill in the operators in the Z array, to be contracted with C in
  // the construction of the MPO. The order is as listed above.
  Zferm=mwArray(Indices(nrOps,d*d,d*d));
  for(int d1=0;d1<d*d;d1++){
    for(int d2=0;d2<d*d;d2++){
      for(int k=0;k<nrOps;k++){
	Zferm.setElement(terms[k].getElement(Indices(d1,d2)),Indices(k,d1,d2));
      }
    }
  }
  Zferm.reshape(Indices(nrOps,d*d*d*d));
}

void KondoHamiltonian::getHMPOwithPenalty(MPO& mpo,double penalty,double Starget) const{
  // Basically repeat the construction of H and include the penalty
  mpo.initLength(hamil.getLength());
  int D=7; // bond dimension of the MPO
  int Dr_s=D-1; // D-2
  {
    // I create the spin site first, for that's unique
    mwArray C(Indices(1,Dr_s,4)); // two less than the others
    C.setElement(-.5*Jpar*ONE_c,Indices(0,1,1)); // J/2 S^+
    C.setElement(-.5*Jpar*ONE_c,Indices(0,2,2)); // J/2 S^-
    C.setElement(.25*.5*Jperp*ONE_c,Indices(0,3,3)); // J/2 S^z *2 (spin matrix)
    C.setElement(ONE_c,Indices(0,0,0)); // just pass on Id
    C.setElement((.25+Starget*Starget)*penalty*ONE_c,Indices(0,Dr_s-1,0)); // the global constant from the penalty
    C.setElement(omega*.5*ONE_c-penalty*Starget*ONE_c,Indices(0,Dr_s-1,3)); // single site term for the impurity
    C.setElement(penalty*ONE_c,Indices(0,4,3)); // long range term
    C.reshape(Indices(Dr_s,4));
    C.multiplyRight(Zsp);
    C.reshape(Indices(1,Dr_s,d,d));
    C.permute(Indices(3,1,4,2));
    mpo.setOp(0,new Operator(C),true);
  }
  // Now construct the tensors for the L+1 fermionic sites
  for(int k=0;k<=L;k++){
    int Dl=(k==0)?Dr_s:D;
    int Dr=(k==L)?1:D;
    mwArray C(Indices(Dl,Dr,nrOps)); 
    if(k<L){
      C.setElement(ONE_c,Indices(0,0,0)); // nothing yet, only Id
      C.setElement(ONE_c,Indices(Dl-2,Dr-2,0)); // in the middle of long range
    }
    C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // after everything
    C.setElement(-.5*Starget*penalty*ONE_c,Indices(0,Dr-1,11)); // local part of the penalty
    //    C.setElement(-.5*.25*penalty*ONE_c,Indices(0,Dr-1,12)); // local part of the penalty
    C.setElement(.25*.25*penalty*ONE_c,Indices(0,Dr-1,12)); // local part of the penalty
    complex_t tcompl=k==0?sqrt(2)*t*ONE_c:t*ONE_c; // hopping parameter
    if(k==0){ // the first one is special
      C.setElement(ONE_c,Indices(1,Dr-1,10)); // second part of: S+ F-
      C.setElement(ONE_c,Indices(2,Dr-1,9));  // S- F+
      C.setElement(ONE_c,Indices(3,Dr-1,11)); // Sz Fz 
      C.setElement(.25*ONE_c,Indices(4,Dr-1,11)); // the long range term from penalty
    }
    else{
      // second part of the hopping term
      C.setElement(ONE_c,Indices(1,Dr-1,5)); 
      C.setElement(ONE_c,Indices(2,Dr-1,6)); 
      C.setElement(ONE_c,Indices(3,Dr-1,7)); 
      C.setElement(ONE_c,Indices(4,Dr-1,8)); 
      // second part of the penalty long range term
      C.setElement(.25*ONE_c,Indices(5,Dr-1,11)); 
    }
    if(k<L){
      // first part of the hopping term
      C.setElement(tcompl,Indices(0,1,1)); 
      C.setElement(tcompl,Indices(0,2,2)); 
      C.setElement(tcompl,Indices(0,3,3)); 
      C.setElement(tcompl,Indices(0,4,4));
      // first part of the penalty long range term
      C.setElement(.5*penalty*ONE_c,Indices(0,5,11)); 
    }
    C.reshape(Indices(Dl*Dr,nrOps));
    C.multiplyRight(Zferm);
    C.reshape(Indices(Dl,Dr,d*d,d*d));
    C.permute(Indices(3,1,4,2));
    mpo.setOp(k+1,new Operator(C),true);    
  }

}

void KondoHamiltonian::initHMPO(){
  int D=6; // bond dimension of the MPO
  int Dr_s=D-1; // D-2
  {
    // I create the spin site first, for that's unique
    mwArray C(Indices(1,Dr_s,4)); // two less than the others
    C.setElement(-.5*Jpar*ONE_c,Indices(0,1,1)); // J/2 S^+
    C.setElement(-.5*Jpar*ONE_c,Indices(0,2,2)); // J/2 S^-
    C.setElement(.25*.5*Jperp*ONE_c,Indices(0,3,3)); // J/2 S^z *2 (spin matrix)
    C.setElement(ONE_c,Indices(0,0,0)); // just pass on Id
    C.setElement(omega*.5*ONE_c,Indices(0,Dr_s-1,3)); // single site term for the impurity
    C.reshape(Indices(Dr_s,4));
    C.multiplyRight(Zsp);
    C.reshape(Indices(1,Dr_s,d,d));
    C.permute(Indices(3,1,4,2));
    hamil.setOp(0,new Operator(C),true);
  }
  // Now construct the tensors for the L+1 fermionic sites
  for(int k=0;k<=L;k++){
    int Dl=(k==0)?Dr_s:D;
    int Dr=(k==L)?1:D;
    mwArray C(Indices(Dl,Dr,nrOps)); 
    if(k<L)
      C.setElement(ONE_c,Indices(0,0,0)); // nothing yet, only Id
    //if(k>0)
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // after everything
    complex_t tcompl=k==0?sqrt(2)*t*ONE_c:t*ONE_c; // hopping parameter
    if(k==0){ // the first one is special
      C.setElement(ONE_c,Indices(1,Dr-1,10)); // second part of: S+ F-
      C.setElement(ONE_c,Indices(2,Dr-1,9));  // S- F+
      C.setElement(ONE_c,Indices(3,Dr-1,11)); // Sz Fz 
    }
    else{
      // second part of the hopping term
      C.setElement(ONE_c,Indices(1,Dr-1,5)); 
      C.setElement(ONE_c,Indices(2,Dr-1,6)); 
      C.setElement(ONE_c,Indices(3,Dr-1,7)); 
      C.setElement(ONE_c,Indices(4,Dr-1,8)); 
    }
    if(k<L){
      // first part of the hopping term
      C.setElement(tcompl,Indices(0,1,1)); 
      C.setElement(tcompl,Indices(0,2,2)); 
      C.setElement(tcompl,Indices(0,3,3)); 
      C.setElement(tcompl,Indices(0,4,4)); 
    }
    C.reshape(Indices(Dl*Dr,nrOps));
    C.multiplyRight(Zferm);
    C.reshape(Indices(Dl,Dr,d*d,d*d));
    C.permute(Indices(3,1,4,2));
    hamil.setOp(k+1,new Operator(C),true);    
  }
}

void KondoHamiltonian::getFermionicOperators(mwArray& Fx,mwArray& Fy,mwArray& Fz) const{
  mwArray sig0=identityMatrix(d);//identity
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP);//sigmaPlus
  mwArray sigM=Hconjugate(sigP); //sigmaMinus
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz

  mwArray auxPM;
  constructOperatorProduct(auxPM,sigP,sigM);

  Fx=(-1.)*auxPM-Hconjugate(auxPM);
  Fy=I_c*auxPM-I_c*Hconjugate(auxPM);

  constructOperatorProduct(Fz,sigZ,sig0);
  mwArray aux0Z;
  constructOperatorProduct(aux0Z,sig0,sigZ);
  Fz=.5*(Fz-aux0Z);
}

void KondoHamiltonian::getIntegratedEvenCorrelationsMPO(MPO& mpo) const{
  // Instead of creating again the ops, could use the ones in
  // Zsp,Zferm, but they are written in terms of Splus, Sminus
  mpo.initLength(L+2);  
  // spin operators
  complex_t dataX[]={ZERO_c,.5*ONE_c,.5*ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,.5*I_c,-.5*I_c,ZERO_c};
  complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
  mwArray Sx(Indices(d,d),dataX);//Sx 
  mwArray Sy(Indices(d,d),dataY);//Sy
  mwArray Sz(Indices(d,d),dataZ);//Sz
  mwArray Fx,Fy,Fz; // fermion spin operators
  getFermionicOperators(Fx,Fy,Fz);
  Fx=.5*Fx;Fy=.5*Fy;Fz=.5*Fz; // proper spin ops
  mwArray F0=identityMatrix(d*d);
  
  int D=4; // bond dimension of the MPO
  {
    // Construct the oper for impurity site
    mwArray Zi(Indices(3,d,d));
    for(int d1=0;d1<d;d1++){
      for(int d2=0;d2<d;d2++){
	Zi.setElement(Sx.getElement(Indices(d1,d2)),Indices(0,d1,d2));
	Zi.setElement(Sy.getElement(Indices(d1,d2)),Indices(1,d1,d2));
	Zi.setElement(Sz.getElement(Indices(d1,d2)),Indices(2,d1,d2));
      }
    }
    Zi.reshape(Indices(3,d*d));
    mwArray C(Indices(1,D,3));C.fillWithZero(); 
    C.setElement(ONE_c,Indices(0,0,0)); // Sx
    C.setElement(ONE_c,Indices(0,1,1)); // Sy
    C.setElement(ONE_c,Indices(0,2,2)); // Sz
    C.reshape(Indices(D,3));
    C.multiplyRight(Zi);
    C.reshape(Indices(1,D,d,d));
    C.permute(Indices(3,1,4,2));
    mpo.setOp(0,new Operator(C),true);
  }
  // Now the opers for even and odd sites (index of ferm)
  mwArray Zf(Indices(4,d*d,d*d));
  for(int d1=0;d1<d*d;d1++){
    for(int d2=0;d2<d*d;d2++){
      Zf.setElement(F0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
      Zf.setElement(Fx.getElement(Indices(d1,d2)),Indices(1,d1,d2));
      Zf.setElement(Fy.getElement(Indices(d1,d2)),Indices(2,d1,d2));
      Zf.setElement(Fz.getElement(Indices(d1,d2)),Indices(3,d1,d2));
    }
  }
  Zf.reshape(Indices(4,d*d*d*d));
  int evOp(-1),odOp(-1);
  for(int k=0;k<L+1;k++){
    if(k==L||(k%2==0&&evOp<0)||(k%2!=0&&odOp<0)){ // need to create the op
      int Dl=D;
      int Dr=(k==L)?1:D;
      mwArray C(Indices(Dl,Dr,4));C.fillWithZero();
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // after everyting is done
      if(k<L){
	C.setElement(ONE_c,Indices(0,0,0)); // pass everything
	C.setElement(ONE_c,Indices(1,1,0)); // pass everything
	C.setElement(ONE_c,Indices(2,2,0)); // pass everything
      }
      if(k%2==0){ // if even, put second term
	C.setElement(ONE_c,Indices(0,Dr-1,1)); 
	C.setElement(ONE_c,Indices(1,Dr-1,2)); 
	C.setElement(ONE_c,Indices(2,Dr-1,3)); 
      }
      C.reshape(Indices(Dl*Dr,4));
      C.multiplyRight(Zf);
      C.reshape(Indices(Dl,Dr,d*d,d*d));
      C.permute(Indices(3,1,4,2));
      mpo.setOp(k+1,new Operator(C),true);
      if(k%2==0) evOp=k+1;
      else odOp=k+1;
    }
    else{
      if(k%2==0)
	mpo.setOp(k+1,&mpo.getOp(evOp),false);
      else
	mpo.setOp(k+1,&mpo.getOp(odOp),false);
    }    
  }  
}


void KondoHamiltonian::getFermionNrMPO(MPO& mpo) const {
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz
  mwArray sig0=identityMatrix(d);

  mwArray nrUp;
  constructOperatorProduct(nrUp,.5*(sigZ+sig0),sig0);
  mwArray nrDown;
  constructOperatorProduct(nrDown,sig0,.5*(sig0+sigZ));
  mwArray nrLoc=nrUp+nrDown;
  mwArray sig00=identityMatrix(d*d);

  mwArray Zf(Indices(2,d*d,d*d));
  for(int d1=0;d1<d*d;d1++){
    for(int d2=0;d2<d*d;d2++){
      Zf.setElement(sig00.getElement(Indices(d1,d2)),Indices(0,d1,d2));
      Zf.setElement(nrLoc.getElement(Indices(d1,d2)),Indices(1,d1,d2));
    }
  }
  Zf.reshape(Indices(2,d*d*d*d));
  
  mpo.initLength(L+2);
  sig0.reshape(Indices(d,1,d,1));
  mpo.setOp(0,new Operator(sig0),true);
  for(int k=0;k<L+1;k++){
    if(k<=1||k==L){ // the different ones
      int Dl=(k==0)?1:2;
      int Dr=(k==L)?1:2;
      mwArray C(Indices(Dl,Dr,2));
      if(k<L)
	C.setElement(ONE_c,Indices(0,0,0)); // except last one
      if(k>0)
	C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // except first one
      C.setElement(ONE_c,Indices(0,Dr-1,1));
      C.reshape(Indices(Dl*Dr,2));C.multiplyRight(Zf);
      C.reshape(Indices(Dl,Dr,d*d,d*d));
      C.permute(Indices(3,1,4,2));
      mpo.setOp(k+1,new Operator(C),true); 
    }
    else{
      mpo.setOp(k+1,&mpo.getOp(2),false); // could use a link
    }
  }
  
}

void KondoHamiltonian::getTotalSzMPO(MPO& mpo) const {
  mpo.initLength(L+2);
  // The first site is the impurity: Id or .5*sigZ
  {int Dl(1),Dr(2);
    mwArray Csp(Indices(Dl,Dr,4));
    Csp.setElement(ONE_c,Indices(0,0,0));
    Csp.setElement(.5*ONE_c,Indices(0,1,3));
    Csp.reshape(Indices(Dl*Dr,4));
    Csp.multiplyRight(Zsp);
    Csp.reshape(Indices(Dl,Dr,d,d));
    Csp.permute(Indices(3,1,4,2));
    mpo.setOp(0,new Operator(Csp),true);     
  }
  for(int k=0;k<L+1;k++){
    if(k==0||k==L){ // the different ones
      int Dl=2;
      int Dr=(k==L)?1:2;
      mwArray C(Indices(Dl,Dr,nrOps));
      if(k<L)
	C.setElement(ONE_c,Indices(0,0,0)); // except last one
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); 
      C.setElement(.25*ONE_c,Indices(0,Dr-1,11));
      C.reshape(Indices(Dl*Dr,nrOps));C.multiplyRight(Zferm);
      C.reshape(Indices(Dl,Dr,d*d,d*d));
      C.permute(Indices(3,1,4,2));
      mpo.setOp(k+1,new Operator(C),true); 
    }
    else{
      mpo.setOp(k+1,&mpo.getOp(1),false); // could use a link
    }
  }
  
}


void KondoHamiltonian::getExponentialMPOeven(MPO& expHe,complex_t delta) const {
   getExponentialMPO(expHe,delta,true);
}

void KondoHamiltonian::getExponentialMPOodd(MPO& expHo,complex_t delta) const {
  getExponentialMPO(expHo,delta,false);
}


void KondoHamiltonian::getExponentialMPO(MPO& expH,complex_t delta,bool even) const {
  mwArray Ol,Or;
  expH.initLength(L+2);
  // if even, I need to include the special case of impurity-first
  // fermion; if odd, the first term is also different, because the
  // coefficient of the hopping term is different. All the other
  // fermion terms are the same.
  //cout<<"KondoHamiltonian::getExponentialMPO(delta="<<delta<<",even="<<even<<")"<<endl;
  int kfirst=even?0:1;
  int len=L+2;
  int klast=(len-1)-(len-kfirst)%2; // last occupied by the loop
  // The first exponential is always special
  bool repeat=false;int posOl(0),posOr(0);
  for(int k=kfirst;k<klast;k+=2){
    if(!repeat){
      getTwoBodyTermExponential(Ol,Or,delta,k);
      expH.setOp(k,new Operator(Ol),true);
      expH.setOp(k+1,new Operator(Or),true);
      if(k>kfirst){ // from now on, repeat these operators
	repeat=true;
	posOl=k;posOr=k+1;
	//cout<<"Will repeat ops Ol("<<k<<") Or("<<k+1<<")"<<endl;
      }
      //cout<<"Set original operators on "<<k<<" and "<<k+1<<endl;
    }
    else{
      expH.setOp(k,&expH.getOp(posOl),false);
      expH.setOp(k+1,&expH.getOp(posOr),false);
      //cout<<"Set repeated operators on "<<k<<" and "<<k+1<<endl;
    }
  }
  // Fill in the edges with identity operators
  if(kfirst==1){ // the first site is the identity
    mwArray ident=identityMatrix(d);
    ident.reshape(Indices(d,1,d,1)); // just dx1xdx1
    expH.setOp(0,new Operator(ident),true);
    //cout<<"Set pos "<<0<<" to Id"<<endl;
  }
  if(klast<len-1){ // sth to be filled at the end
    mwArray ident=identityMatrix(d*d);
    ident.reshape(Indices(d*d,1,d*d,1)); // just dx1xdx1
    expH.setOp(len-1,new Operator(ident),true);
    //cout<<"Set pos "<<len-1<<" to Id"<<endl;
  }
}

void KondoHamiltonian::getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,
						 int pos) const {
  mwArray H12;
  computeTwoBodyTerm(H12,pos);
  // Now take the matrix exponential
  mwArray expH;
  //  if(pos==2) cout<<"h12("<<pos<<")="<<H12<<endl;
  wrapper::expm(H12,expH,delta);
  //putForMatlab(cout,expH,"expH12");
  //if(pos==2) cout<<"Exponential of the h12("<<pos<<")="<<expH<<endl;
  int d2=d*d;
  int d1=(pos==0)?d:d2;
  // Obtained with indices (i'j')(ij) => permute to (i'i)(j'j) for svd
  expH.reshape(Indices(d1,d2,d1,d2));
  expH.permute(Indices(1,3,2,4));
  expH.reshape(Indices(d1*d1,d2*d2));
  // And compute the SVD
  mwArray S; // place for the singular values
  int nr=0;
  double tol=0.;
  //cout<<"Before decomposition, expH="<<expH<<endl;
  wrapper::svd(expH,tol,nr,Ol,S,Or);
  //if(pos==2) cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<", S="<<S<<endl;
  // redistribute the singular values
  S=sqrt(S);
  // TODO: check no NaN???
  Ol.multiplyRight(S);
  Or.multiplyLeft(S);
  // reshape for operators (i'x 1 x i x alpha ))
  Ol.reshape(Indices(d1,1,d1,-1));
  Or.reshape(Indices(-1,d2,d2,1));
  Or.permute(Indices(2,1,3,4));  
  //cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<endl;
  //putForMatlab(cout,Ol,"Ol");
  //putForMatlab(cout,Or,"Or");
}

void KondoHamiltonian::computeTwoBodyTerm(mwArray& result,int pos) const{
  if(pos>=L+1){
    cout<<"ERROR: HeisenbergHamiltonian::computeTwoBodyTerm for site "<<pos<<" when L="<<L<<endl;
    exit(212);
  }
  if(pos==0){ // case spin-fermion
    computeTwoBodyTermSF(result);
  }
  else{
    bool firstFermion=(pos==1);
    computeTwoBodyTermFF(result,firstFermion);
  }
}

void KondoHamiltonian::computeTwoBodyTermSF(mwArray& result) const {
  // J/2 (S^+ F- + S^- F+) + J/2 S_z F_z
  // Basic pieces for the spin term
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP);//sigmaPlus
  mwArray sigM=Hconjugate(sigP); //sigmaMinus
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz); //sigmaz (S_z is .5 this)
  mwArray sig0=identityMatrix(d);
  // Basic pieces for the fermion term
  // F^+= sigP (x) tauM
  mwArray Fplus,Fminus,Fz;
  constructOperatorProduct(Fplus,sigP,sigM);
  // F^-= sigM (x) tauP
  constructOperatorProduct(Fminus,sigM,sigP);
  // F^z= .5*(1+sigZ) (x) Id - Id (x) (1+tauZ)/2
  mwArray aux;
  constructOperatorProduct(Fz,sigZ,sig0);
  constructOperatorProduct(aux,sig0,sigZ);
  Fz=.5*(Fz-aux);
  sigP.reshape(Indices(d*d,1));
  Fminus.reshape(Indices(1,d*d*d*d));
  sigP.multiplyRight(Fminus);
  sigM.reshape(Indices(d*d,1));
  Fplus.reshape(Indices(1,d*d*d*d));
  sigM.multiplyRight(Fplus);
  sigZ.reshape(Indices(d*d,1));
  mwArray singleTerm(sigZ);
  singleTerm.multiplyRight(reshape(identityMatrix(d*d),Indices(1,d*d*d*d)));
  Fz.reshape(Indices(1,d*d*d*d));
  sigZ.multiplyRight(Fz);
  result=-.5*Jpar*(sigP+sigM)+.25*Jperp*sigZ+omega*.5*singleTerm;
  result.reshape(Indices(d,d,d*d,d*d));
  result.permute(Indices(1,3,2,4));
  result.reshape(Indices(d*d*d,d*d*d));
}

void KondoHamiltonian::computeTwoBodyTermFF(mwArray& result,bool first) const{
  // Hopping term: four components. Construct again the fermion terms,
  // although I could save them from initHMPO().
  // Basic pieces for the spin terms
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP);//sigmaPlus
  mwArray sigM=Hconjugate(sigP); //sigmaMinus
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz); //sigmaz (S_z is .5 this)
  mwArray sig0=identityMatrix(d);
  vector<mwArray> terms(9,mwArray::empty);
  // *** First of all, we keep the identity on both
  int d2=d*d;
  terms[0]=identityMatrix(d2); // keep for clarity of numbering. 
  // *** Now the ones appearing in fermion hopping terms
  // (1) sigP (x) tauZ
  constructOperatorProduct(terms[1],sigP,sigZ);
  terms[1].reshape(Indices(d2*d2,1));
  // (2) sigM (x) tauZ
  constructOperatorProduct(terms[2],sigM,sigZ);  
  terms[2].reshape(Indices(d2*d2,1));
  // (3) Id (x) tauP
  constructOperatorProduct(terms[3],sig0,sigP);  
  terms[3].reshape(Indices(d2*d2,1));
  // (4) Id (x) tauM
  constructOperatorProduct(terms[4],sig0,sigM);  
  terms[4].reshape(Indices(d2*d2,1));
  // (5) sigM (x) Id
  constructOperatorProduct(terms[5],sigM,sig0);
  terms[5].reshape(Indices(1,d2*d2));
  // (6) sigP (x) Id
  constructOperatorProduct(terms[6],sigP,sig0);
  terms[6].reshape(Indices(1,d2*d2));
  // (7) sigZ (x) tauM
  constructOperatorProduct(terms[7],sigZ,sigM);
  terms[7].reshape(Indices(1,d2*d2));
  // (8) sigZ (x) tauP
  constructOperatorProduct(terms[8],sigZ,sigP);
  terms[8].reshape(Indices(1,d2*d2));

  result=terms[1]*terms[5]+terms[2]*terms[6]+terms[3]*terms[7]+terms[4]*terms[8];

  result=t*result;
  if(first) result=sqrt(2)*result;
  result.reshape(Indices(d2,d2,d2,d2));
  result.permute(Indices(1,3,2,4));
  result.reshape(Indices(d2*d2,d2*d2));
}


