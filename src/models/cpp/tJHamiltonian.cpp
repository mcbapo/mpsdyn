#include "tJHamiltonian.h"
#include "Indices.h"

using namespace shrt;
using namespace std;

tJHamiltonian::tJHamiltonian(int L_,double t_,double J_):d(2),L(L_),t(t_),J(J_),hamil(L){
  initZ();
  initHMPO();
}

tJHamiltonian::~tJHamiltonian(){};


void tJHamiltonian::initZ(){
    // Basic pieces
  mwArray sig0=identityMatrix(d);//identity
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP);//sigmaPlus
  mwArray sigM=Hconjugate(sigP); //sigmaMinus
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz

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
  // (9) Sx=(sigP (x) tauM + sigM (x) tauP)/2
  mwArray auxPM;
  constructOperatorProduct(auxPM,sigP,sigM);
  terms[9]=.5*(auxPM+Hconjugate(auxPM));
  // (10) Sy=i(sigP (x) tauM - sigM (x) tauP)/2
  terms[10]=.5*I_c*(auxPM-Hconjugate(auxPM));
  // (11) Sz=(sigZ (x) Id - Id (x) tauZ)/4
  mwArray ZId,IdZ;
  constructOperatorProduct(ZId,sigZ,sig0);
  constructOperatorProduct(IdZ,sig0,sigZ);
  terms[11]=.25*(ZId-IdZ);
  // (12) n=Id+(sigZ (x) Id + Id (x) tauZ)/2
  terms[12]=identityMatrix(d*d)+.5*(ZId+IdZ);
  // // (9) (1+sigZ)/2 (x) (1+tauZ)/2
  // mwArray projZ=.5*(identityMatrix(d)+sigZ);
  // constructOperatorProduct(terms[9],projZ,projZ);
  
  // Fill in the operators in the Z array, to be contracted with C in
  // the construction of the MPO. The order is as listed above.
  Z=mwArray(Indices(nrOps,d*d,d*d));
  for(int d1=0;d1<d*d;d1++){
    for(int d2=0;d2<d*d;d2++){
      for(int k=0;k<nrOps;k++){
	Z.setElement(terms[k].getElement(Indices(d1,d2)),Indices(k,d1,d2));
      }
    }
  }
  Z.reshape(Indices(nrOps,d*d*d*d));
}


void tJHamiltonian::initHMPO(){
  int D=10; // bond dimension of the MPO
  // Tensors site by site
  for(int k=0;k<L;k++){
    int Dl=(k==0)?1:D;
    int Dr=(k==L-1)?1:D;
    mwArray C(Indices(Dl,Dr,nrOps)); 
    if(k<L-1)
      C.setElement(ONE_c,Indices(0,0,0)); // nothing yet, only Id
    if(k>0)
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // after everything
    complex_t tcompl=t*ONE_c; // hopping parameter (notice the minus sign is absorbed in the spin ops)
    complex_t Jcompl=J*ONE_c; // interaction parameter
    if(k>0){
      // second part of the hopping term
      C.setElement(ONE_c,Indices(1,Dr-1,5)); 
      C.setElement(ONE_c,Indices(2,Dr-1,6)); 
      C.setElement(ONE_c,Indices(3,Dr-1,7)); 
      C.setElement(ONE_c,Indices(4,Dr-1,8));
      // second part of the interaction terms
      C.setElement(ONE_c,Indices(5,Dr-1,9)); 
      C.setElement(ONE_c,Indices(6,Dr-1,10)); 
      C.setElement(ONE_c,Indices(7,Dr-1,11)); 
      C.setElement(ONE_c,Indices(8,Dr-1,12));
    }
    if(k<L-1){
      // first part of the hopping term
      C.setElement(tcompl,Indices(0,1,1)); 
      C.setElement(tcompl,Indices(0,2,2)); 
      C.setElement(tcompl,Indices(0,3,3)); 
      C.setElement(tcompl,Indices(0,4,4));
      // first part of the interaction terms
      C.setElement(Jcompl,Indices(0,5,9)); 
      C.setElement(Jcompl,Indices(0,6,10)); 
      C.setElement(Jcompl,Indices(0,7,11)); 
      C.setElement(-.25*Jcompl,Indices(0,8,12));      
    }
    // No single body term
    C.reshape(Indices(Dl*Dr,nrOps));
    C.multiplyRight(Z);
    C.reshape(Indices(Dl,Dr,d*d,d*d));
    C.permute(Indices(3,1,4,2));
    hamil.setOp(k,new Operator(C),true);    
  }
}

void tJHamiltonian::getFermionicOperators(mwArray& Fx,mwArray& Fy,mwArray& Fz) const{
  mwArray sig0=identityMatrix(d);//identity
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP);//sigmaPlus
  mwArray sigM=Hconjugate(sigP); //sigmaMinus
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz

  mwArray auxPM;
  constructOperatorProduct(auxPM,sigP,sigM);

  Fx=.5*((-1.)*auxPM-Hconjugate(auxPM));
  Fy=.5*(I_c*auxPM-I_c*Hconjugate(auxPM));

  constructOperatorProduct(Fz,sigZ,sig0);
  mwArray aux0Z;
  constructOperatorProduct(aux0Z,sig0,sigZ);
  Fz=.25*(Fz-aux0Z);
}


void tJHamiltonian::getFermionNrMPO(MPO& mpo) const {
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
  
  mpo.initLength(L);
  for(int k=0;k<L;k++){
    if(k<=1||k==L-1){ // the different ones
      int Dl=(k==0)?1:2;
      int Dr=(k==L-1)?1:2;
      mwArray C(Indices(Dl,Dr,2));
      if(k<L-1)
	C.setElement(ONE_c,Indices(0,0,0)); // except last one
      if(k>0)
	C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // except first one
      C.setElement(ONE_c,Indices(0,Dr-1,1));
      C.reshape(Indices(Dl*Dr,2));C.multiplyRight(Zf);
      C.reshape(Indices(Dl,Dr,d*d,d*d));
      C.permute(Indices(3,1,4,2));
      mpo.setOp(k,new Operator(C),true);
    }
    else
      mpo.setOp(k,&mpo.getOp(1),false); // for the rest I can use the same link from 1 to L-2
  }
  
}

void tJHamiltonian::getTotalSzMPO(MPO& mpo) const {
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz
  mwArray sig0=identityMatrix(d);

  mwArray nrUp;
  constructOperatorProduct(nrUp,.5*(sigZ+sig0),sig0);
  mwArray nrDown;
  constructOperatorProduct(nrDown,sig0,.5*(sig0+sigZ));
  mwArray aux=nrUp-nrDown; // TODO!!!! CHECK?
  mwArray sig00=identityMatrix(d*d);
  // // F_0^z= .5*(1+sigZ) (x) Id - Id (x) (1+tauZ)/2 

  mwArray Zf(Indices(2,d*d,d*d));
  for(int d1=0;d1<d*d;d1++){
    for(int d2=0;d2<d*d;d2++){
      Zf.setElement(sig00.getElement(Indices(d1,d2)),Indices(0,d1,d2));
      Zf.setElement(aux.getElement(Indices(d1,d2)),Indices(1,d1,d2));
    }
  }
  Zf.reshape(Indices(2,d*d*d*d));
  
  mpo.initLength(L);
  for(int k=0;k<L;k++){
    if(k<=1||k==L-1){ // the different ones
      int Dl=(k==0)?1:2;
      int Dr=(k==L-1)?1:2;
      mwArray C(Indices(Dl,Dr,2));
      if(k<L-1)
	C.setElement(ONE_c,Indices(0,0,0)); // except last one
      if(k>0)
	C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // not first one
      C.setElement(ONE_c,Indices(0,Dr-1,1));
      C.reshape(Indices(Dl*Dr,2));C.multiplyRight(Zf);
      C.reshape(Indices(Dl,Dr,d*d,d*d));
      C.permute(Indices(3,1,4,2));
      mpo.setOp(k,new Operator(C),true); 
    }
    else{
      mpo.setOp(k,&mpo.getOp(1),false); // could use a link
    }
  }  
}


void tJHamiltonian::getExponentialMPOeven(MPO& expHe,complex_t delta) const {
   getExponentialMPO(expHe,delta,true);
}

void tJHamiltonian::getExponentialMPOodd(MPO& expHo,complex_t delta) const {
  getExponentialMPO(expHo,delta,false);
}


void tJHamiltonian::getExponentialMPO(MPO& expH,complex_t delta,bool even) const {
  mwArray Ol,Or;
  expH.initLength(L);
  //cout<<"tJHamiltonian::getExponentialMPO(delta="<<delta<<",even="<<even<<")"<<endl;
  int kfirst=even?0:1;
  int len=L;
  int klast=(len-1)-(len-kfirst)%2; // last occupied by the loop
  // The first exponential is always special (the rest are copies, as I put the U term on the left site)
  bool repeat=false;int posOl(0),posOr(0);
  for(int k=kfirst;k<klast;k+=2){
    if(!repeat){
      getTwoBodyTermExponential(Ol,Or,delta,k);
      expH.setOp(k,new Operator(Ol),true);
      expH.setOp(k+1,new Operator(Or),true);
      // from now on, repeat these operators
      repeat=true;
      posOl=k;posOr=k+1;
      //cout<<"Will repeat ops Ol("<<k<<") Or("<<k+1<<")"<<endl;
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
    mwArray ident=identityMatrix(d*d);
    ident.reshape(Indices(d*d,1,d*d,1)); // just dx1xdx1
    expH.setOp(0,new Operator(ident),true);
    //cout<<"Set pos "<<0<<" to Id"<<endl;
  }
  if(klast<len-1){ // sth to be filled at the end:
    mwArray ident=identityMatrix(d*d);
    ident.reshape(Indices(d*d,1,d*d,1)); // just dx1xdx1
    expH.setOp(len-1,new Operator(ident),true);
    //cout<<"Set pos "<<len-1<<" to Id"<<endl;
  }
}

void tJHamiltonian::getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,
						 int pos) const {
  mwArray H12;
  computeTwoBodyTerm(H12,pos);
  // Now take the matrix exponential
  mwArray expH;
  //  if(pos==2) cout<<"h12("<<pos<<")="<<H12<<endl;
  wrapper::expm(H12,expH,delta);
  //putForMatlab(cout,expH,"expH12");
  if(pos==L-1){// special case: singlebody term for the last site
    Ol=expH;
    Ol.reshape(Indices(d*d,1,d*d,1));
    return;
  }

  //if(pos==2) cout<<"Exponential of the h12("<<pos<<")="<<expH<<endl;
  int d2=d*d;
  int d1=d2;
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

void tJHamiltonian::computeTwoBodyTerm(mwArray& result,int pos) const{
  if(pos>=L-1){
    cout<<"ERROR: tJHamiltonian::computeTwoBodyTerm for site "<<pos<<" when L="<<L<<endl;
    exit(212);
  }
  // pos is mostly ignored, otherwise, as everything is TI

  // Hopping term: four components. Construct again the fermion terms,
  // although I could save them from initHMPO().
  // Basic pieces for the spin terms
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP);//sigmaPlus
  mwArray sigM=Hconjugate(sigP); //sigmaMinus
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz); //sigmaz (S_z is .5 this)
  mwArray sig0=identityMatrix(d);
  vector<mwArray> terms(13,mwArray::empty);
  // *** First of all, we keep the identity on both
  int d2=d*d;
  terms[0]=identityMatrix(d2); 
  terms[0].reshape(Indices(1,d2*d2));
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
  // (9) Sx=(sigP (x) tauM + sigM (x) tauP)/2
  mwArray auxPM;
  constructOperatorProduct(auxPM,sigP,sigM);
  terms[9]=.5*(auxPM+Hconjugate(auxPM));
  terms[9].reshape(Indices(d2*d2,1));
  // (10) Sy=i(sigP (x) tauM - sigM (x) tauP)/2
  terms[10]=.5*I_c*(auxPM-Hconjugate(auxPM));
  terms[10].reshape(Indices(d2*d2,1));
  // (11) Sz=(sigZ (x) Id - Id (x) tauZ)/4
  mwArray ZId,IdZ;
  constructOperatorProduct(ZId,sigZ,sig0);
  constructOperatorProduct(IdZ,sig0,sigZ);
  terms[11]=.25*(ZId-IdZ);
  terms[11].reshape(Indices(d2*d2,1));
  // (12) n=Id+(sigZ (x) Id + Id (x) tauZ)/2
  terms[12]=identityMatrix(d*d)+.5*(ZId+IdZ);
  terms[12].reshape(Indices(d2*d2,1));
  result=terms[1]*terms[5]+terms[2]*terms[6]+terms[3]*terms[7]+terms[4]*terms[8];
  result=t*result+J*(terms[9]*permute(terms[9],Indices(2,1))+
		     terms[10]*permute(terms[10],Indices(2,1))+
		     terms[11]*permute(terms[11],Indices(2,1))-
		     .25*terms[12]*permute(terms[12],Indices(2,1)));
  result.reshape(Indices(d2,d2,d2,d2));
  result.permute(Indices(1,3,2,4));
  result.reshape(Indices(d2*d2,d2*d2));
}


