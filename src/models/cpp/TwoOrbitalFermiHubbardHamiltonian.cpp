#include "TwoOrbitalFermiHubbardHamiltonian.h"
#include "Indices.h"

using namespace shrt;
using namespace std;

// TwoOrbitalFermiHubbardHamiltonian::TwoOrbitalFermiHubbardHamiltonian(int L_,double t_,double U_,double V_,double Vex_,double mu_g_,double mu_e_):
//   d(2),L(L_),t(t_),U(U_),V(V_),Vex(Vex_),mu_g(mu_g_),mu_e(mu_e_),hamil(4*L){
//   initZ();
//   initHMPO();
// }

TwoOrbitalFermiHubbardHamiltonian::TwoOrbitalFermiHubbardHamiltonian(int L_,const vector<double>& t,const vector<double>& U,
								     double V_,double Vex_,
								     const vector<double>& mu_g,const vector<double>& mu_e,
								     int Ntot_,double penalty_):
  d(2),L(L_),hamil(4*L),V(V_),Vex(Vex_),Ntot(Ntot_),penNtot(penalty_),penSz(0.){
  if(t.size()>=4){
    tg0=t[0];tg1=t[1];te0=t[2];te1=t[3];
    if(t.size()>4){tge0=t[4];tge1=t[5];}
    else {tge0=0.;tge1=0.;}
  }
  else{
    if(t.size()==1){
      tg0=t[0];
      tg1=tg0;
      te0=tg0;te1=tg0;
      tge0=0.;tge1=0.;
    }else{ cout<<"Invalid input t for TwoOrbitalFermiHubbardHamiltonian"<<endl;exit(1);}
  }
  if(U.size()==2){
    Ug=U[0];Ue=U[1];
  }
  else{
    if(U.size()==1){
      Ug=U[0];Ue=Ug;
    }else{ cout<<"Invalid input U for TwoOrbitalFermiHubbardHamiltonian"<<endl;exit(1);}
  }
  if(mu_g.size()==2){
    mu_g0=mu_g[0];mu_g1=mu_g[1];
  }
  else{
    if(mu_g.size()==1){
      mu_g0=mu_g[0];mu_g1=mu_g0;
    }else{ cout<<"Invalid input mu_g for TwoOrbitalFermiHubbardHamiltonian"<<endl;exit(1);}
  }
  if(mu_e.size()==2){
    mu_e0=mu_e[0];mu_e1=mu_e[1];
  }
  else{
    if(mu_e.size()==1){
      mu_e0=mu_e[0];mu_e1=mu_e0;
    }else{ cout<<"Invalid input mu_e for TwoOrbitalFermiHubbardHamiltonian"<<endl;exit(1);}
  }
  initZ();
  initHMPO();
}

TwoOrbitalFermiHubbardHamiltonian::~TwoOrbitalFermiHubbardHamiltonian(){};

void TwoOrbitalFermiHubbardHamiltonian::setPenalties(int Sz_,double penaltySz,
						     int Ntot_,double penaltyNtot){
  Ntot=Ntot_;penNtot=penaltyNtot;
  Sz=Sz_;penSz=penaltySz;  
  initHMPO();
}


void TwoOrbitalFermiHubbardHamiltonian::initZ(){
    // Basic pieces
  mwArray sig0=identityMatrix(d);//identity
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP);//sigmaPlus
  mwArray sigM=Hconjugate(sigP); //sigmaMinus
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz
  mwArray projZ0=.5*(sig0+sigZ);

  nrOps=5;

  // Fill in the operators in the Z array, to be contracted with C in
  // the construction of the MPO. The order is Id,sigPlus,sigMinus,sigZ,proj0
  Z=mwArray(Indices(nrOps,d,d));
  for(int d1=0;d1<d;d1++){
    for(int d2=0;d2<d;d2++){
	Z.setElement(sig0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
	Z.setElement(sigP.getElement(Indices(d1,d2)),Indices(1,d1,d2));
	Z.setElement(sigM.getElement(Indices(d1,d2)),Indices(2,d1,d2));
	Z.setElement(sigZ.getElement(Indices(d1,d2)),Indices(3,d1,d2));
	Z.setElement(projZ0.getElement(Indices(d1,d2)),Indices(4,d1,d2));
    }
  }
  Z.reshape(Indices(nrOps,d*d));
}

// TODO: optimize the mPO before running
//#include "Contractor.h"
//#include "MPS.h"

void TwoOrbitalFermiHubbardHamiltonian::initHMPO(){
  int D=17; // bond dimension of the MPO without penalty terms
  int Dz,Dn;
  if(penNtot!=0){
    Dn=D-1;D++; // take the last one for this and increase
    cout<<"Dn="<<Dn<<endl;
  }
  if(penSz!=0){
    Dz=D-1;D++; // take the last one for this and increase
    cout<<"Dz="<<Dz<<endl;
  }
  cout<<"Considering penalty terms, bond dimension of the MPO for the Hamiltonian is "<<D<<endl;
  // Tensors site by site, but sites numbered as (n,alpha,sigma), where n is the site (0 to L-1), alpha the elctronic state (g,e) and sigma the spin (0,1 for up down)
  for(int n=0;n<L;n++){
    for(int alpha=0;alpha<2;alpha++){
      for(int s=0;s<2;s++){
	// linear index from (alpha,s)
	int p=alpha*2+s;
	// linear index
	int k=4*n+p;//cout<<"Site "<<k<<"("<<n<<","<<alpha<<","<<s<<")"<<endl;
	int signS=(s==0)?1:-1; //(-1)^s
	int Dl=(k==0)?1:D;
	int Dr=(k==4*L-1)?1:D;
	mwArray C(Indices(Dl,Dr,nrOps)); 
	if(k<4*L-1)
	  C.setElement(ONE_c,Indices(0,0,0)); // nothing yet, only Id
	if(k>0)
	  C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // after everything
	// First part of hopping (or beginning of -Vex term)
	if(n<L-1||(p==0||(tge1!=0.&&p==1))){
	  C.setElement(ONE_c,Indices(0,p+1,1)); // sigP 
	  C.setElement(ONE_c,Indices(0,p+5,2)); // sigM
	}
	// In the middle of the hoping term
	for(int r=0;r<4;r++){
	  if((r<p&&n<L-1)||(r>p&&n>0)){ // first half for first phys site
	    C.setElement(ONE_c,Indices(r+1,r+1,3)); // sigZ 
	    C.setElement(ONE_c,Indices(r+5,r+5,3)); // sigZ
	  }
	}
	// Second part of hopping carries the subindex (alpha,s)
	complex_t tcompl;
	switch(p){
	case 0: tcompl=tg0*ONE_c;break;
	case 1: tcompl=tg1*ONE_c;break;
	case 2: tcompl=te0*ONE_c;break;
	case 3: tcompl=te1*ONE_c;break;
	}
	if(n>0){
	  C.setElement(tcompl,Indices(p+1,Dr-1,2)); // sigM 
	  C.setElement(tcompl,Indices(p+5,Dr-1,1)); // sigP
	}
	// if hopping between g and e is allowed, I have to include some second terms
	if(tge0!=0.||tge1!=0.){
	  // middle terms also in sites 1 and 2 of the same site, even for n=0 and n=L-1
	  if((n==0||n==L-1)&&(p==1||p==2)){
	    C.setElement(ONE_c,Indices(p+1,p+1,3)); // sigZ 
	    C.setElement(ONE_c,Indices(p+5,p+5,3)); // sigZ
	  }	  
	  complex_t tge_compl;
	  switch(p){ // only the ones that can put the second term
	  case 2: tge_compl=tge0*ONE_c;break;
	  case 3: tge_compl=tge1*ONE_c;break;
	  }
	  if(p>=2){
	    C.setElement(tge_compl+C.getElement(Indices(p+1,Dr-1,2)),
			 Indices(p+1,Dr-1,2)); // sigM 
	    C.setElement(tge_compl+C.getElement(Indices(p+5,Dr-1,1)),
			 Indices(p+5,Dr-1,1)); // sigP
	  }
	}	

	// Vex terms
	if(p==1){
 	  C.setElement(ONE_c,Indices(1,9,2)); // sigM
 	  C.setElement(ONE_c,Indices(5,10,1)); // sigP
 	  C.setElement(ONE_c,Indices(1,14,1)); // sigP
 	  C.setElement(ONE_c,Indices(5,15,2)); // sigM
	}
	if(p==2){
 	  C.setElement(ONE_c,Indices(9,9,2)); // sigM
 	  C.setElement(ONE_c,Indices(10,10,1)); // sigP
 	  C.setElement(ONE_c,Indices(15,15,1)); // sigP
 	  C.setElement(ONE_c,Indices(14,14,2)); // sigM
	}
	if(p==3){
 	  C.setElement(-Vex*ONE_c,Indices(9,Dr-1,1)); // sigP
 	  C.setElement(-Vex*ONE_c,Indices(10,Dr-1,2)); // sigM
 	  C.setElement(Vex*ONE_c,Indices(14,Dr-1,2)); // sigM
 	  C.setElement(Vex*ONE_c,Indices(15,Dr-1,1)); // sigP
	}
	// U term
	if(p==0) C.setElement(ONE_c,Indices(0,11,4)); // proj0
	if(p==1) C.setElement(Ug*ONE_c,Indices(11,Dr-1,4)); // proj0
	if(p==2) C.setElement(ONE_c,Indices(0,12,4)); // proj0
	if(p==3) C.setElement(Ue*ONE_c,Indices(12,Dr-1,4)); // proj0
	// V term
	if(p==1){
 	  C.setElement(ONE_c,Indices(11,11,0)); // Id
	  C.setElement(ONE_c,Indices(0,13,4)); // proj0
	}
	if(p==2){
	  C.setElement(ONE_c,Indices(11,11,0)); // Id
	  C.setElement((V-Vex)*ONE_c,Indices(11,Dr-1,4)); // proj0
	  C.setElement(ONE_c,Indices(13,13,0)); // Id
	  C.setElement(V*ONE_c,Indices(13,Dr-1,4)); // proj0
	}
	if(p==3){
	  C.setElement(V*ONE_c,Indices(11,Dr-1,4)); // proj0
	  C.setElement((V-Vex)*ONE_c,Indices(13,Dr-1,4)); // proj0
	}
	// if penalty terms are present, I need extra terms in the tensors
	if(penNtot!=0){
	  if(k<4*L-1){ // first of the pair
	    C.setElement(ONE_c,Indices(0,Dn,4)); // proj0
	  }
	  if(k>0){ // second of the pair
	    C.setElement(2*penNtot*ONE_c,Indices(Dn,Dr-1,4)); // proj0
	  }
	  if(k>0&&k<4*L-1) C.setElement(ONE_c,Indices(Dn,Dn,0)); // pass it
	}
	if(penSz!=0){
	  if(k<4*L-1){ // first of the pair
	    C.setElement(signS*ONE_c,Indices(0,Dz,4)); // proj0
	  }
	  if(k>0){ // second of the pair
	    C.setElement(2*signS*penSz*ONE_c,Indices(Dz,Dr-1,4)); // proj0
	  }
	  if(k>0&&k<4*L-1) C.setElement(ONE_c,Indices(Dz,Dz,0)); // pass it
	}

	// Single site chemical potential terms
	complex_t mucompl;
	switch(p){
	case 0: mucompl=-mu_g0*ONE_c;break;
	case 1: mucompl=-mu_g1*ONE_c;break;
	case 2: mucompl=-mu_e0*ONE_c;break;
	case 3: mucompl=-mu_e1*ONE_c;break;
	}
	// Add also the local part of penalty terms, if any of them is present
	if(penNtot!=0) mucompl=mucompl+penNtot*(1-2*Ntot)*ONE_c;
	if(penSz!=0) mucompl=mucompl+penSz*(1-2*signS*Sz)*ONE_c;
	C.setElement(mucompl,Indices(0,Dr-1,4)); // proj0

	
	// Add also the constant part of penalty terms, if any of them is present
	double termLoc=0.;
	if(penNtot!=0) termLoc+=penNtot*Ntot*Ntot;
	if(penSz!=0) termLoc+=penSz*Sz*Sz;
	  
	C.setElement((termLoc/(4*L))*ONE_c,Indices(0,Dr-1,0)); // constant term (could drop)

	// Now set the op in the MPO
	C.reshape(Indices(Dl*Dr,nrOps));
	C.multiplyRight(Z);
	C.reshape(Indices(Dl,Dr,d,d));
	C.permute(Indices(3,1,4,2));
	//cout<<"Setting operator of Hmpo at site "<<k<<"("<<n<<","<<alpha<<","<<s<<")"<<endl;
	hamil.setOp(k,new Operator(C),true);    
      }
    }
  }
}

// void TwoOrbitalFermiHubbardHamiltonian::getFermionicOperators(mwArray& Fx,mwArray& Fy,mwArray& Fz) const{
//   mwArray sig0=identityMatrix(d);//identity
//   complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
//   mwArray sigP(Indices(d,d),dataP);//sigmaPlus
//   mwArray sigM=Hconjugate(sigP); //sigmaMinus
//   complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
//   mwArray sigZ(Indices(d,d),dataz);//sigmaz

//   mwArray auxPM;
//   constructOperatorProduct(auxPM,sigP,sigM);

//   Fx=(-1.)*auxPM-Hconjugate(auxPM);
//   Fy=I_c*auxPM-I_c*Hconjugate(auxPM);

//   constructOperatorProduct(Fz,sigZ,sig0);
//   mwArray aux0Z;
//   constructOperatorProduct(aux0Z,sig0,sigZ);
//   Fz=.5*(Fz-aux0Z);
// }


void TwoOrbitalFermiHubbardHamiltonian::getNumberMPO(MPO& mpo,int type) const {
   complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
   mwArray sigZ(Indices(d,d),dataz);//sigmaz
   mwArray sig0=identityMatrix(d);
   mwArray proj0=.5*(sigZ+sig0);
   mwArray Zf(Indices(2,d,d));
   for(int d1=0;d1<d;d1++){
     for(int d2=0;d2<d;d2++){
       Zf.setElement(sig0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
       Zf.setElement(proj0.getElement(Indices(d1,d2)),Indices(1,d1,d2));
     }
   }
   Zf.reshape(Indices(2,d*d));

   mpo.initLength(4*L);
   for(int k=0;k<4*L;k++){
     int Dl=(k==0)?1:2;
     int Dr=(k==4*L-1)?1:2;
     mwArray C(Indices(Dl,Dr,2));
     if(k<4*L-1)
       C.setElement(ONE_c,Indices(0,0,0)); // except last one
     if(k>0)
       C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // except first one
     // in which sites I include the local n operator depends on type
     if(type==0){ // all of them
       C.setElement(ONE_c,Indices(0,Dr-1,1));
     }
     else{
       // what is the orbital corresponding to this k=4*n+2*alpha+s
       int j=k%4; // alpha*2+s
       int alpha=(j-j%2)/2;
       if((type==1&&alpha==0)||(type==2&&alpha==1)){
	 C.setElement(ONE_c,Indices(0,Dr-1,1));
       }
     }
     C.reshape(Indices(Dl*Dr,2));C.multiplyRight(Zf);
     C.reshape(Indices(Dl,Dr,d,d));
     C.permute(Indices(3,1,4,2));
     mpo.setOp(k,new Operator(C),true);
   }
   
}

void TwoOrbitalFermiHubbardHamiltonian::getFermionNrMPO(MPO& mpo) const {
  getNumberMPO(mpo,0);
}

void TwoOrbitalFermiHubbardHamiltonian::getFermionNrOrbitalMPO(MPO& mpo,bool orbital0) const{
  int type=orbital0?1:2;
  getNumberMPO(mpo,type);
}

void TwoOrbitalFermiHubbardHamiltonian::getDoubleOccupancyOrbitalMPO(MPO& mpo,bool orbital0) const{
   complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
   mwArray sigZ(Indices(d,d),dataz);//sigmaz
   mwArray sig0=identityMatrix(d);
   mwArray proj0=.5*(sigZ+sig0);
   mwArray Zf(Indices(2,d,d));
   for(int d1=0;d1<d;d1++){
     for(int d2=0;d2<d;d2++){
       Zf.setElement(sig0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
       Zf.setElement(proj0.getElement(Indices(d1,d2)),Indices(1,d1,d2));
     }
   }
   Zf.reshape(Indices(2,d*d));

   mpo.initLength(4*L);
   int D=3; // bond dim of the MPOs
   for(int k=0;k<4*L;k++){
     int Dl=(k==0)?1:D;
     int Dr=(k==4*L-1)?1:D;
     mwArray C(Indices(Dl,Dr,2));
     if(k<4*L-1)
       C.setElement(ONE_c,Indices(0,0,0)); // except last one
     if(k>0)
       C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // except first one
     // in which pairs of sites I include the local n operator depends on type
     // g-> proj on sites 0 and 1; // e->proj on sites 2 and 3
     if((k%4==0&&orbital0)||(k%4==2&&!orbital0)) // first of the pair
       C.setElement(ONE_c,Indices(0,1,1));
     if((k%4==1&&orbital0)||(k%4==3&&!orbital0)) // second of the pair
       C.setElement(ONE_c,Indices(1,Dr-1,1));
     C.reshape(Indices(Dl*Dr,2));C.multiplyRight(Zf);
     C.reshape(Indices(Dl,Dr,d,d));
     C.permute(Indices(3,1,4,2));
     mpo.setOp(k,new Operator(C),true);
   }
   
}


void TwoOrbitalFermiHubbardHamiltonian::getTotalSzMPO(MPO& mpo) const {
   complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
   mwArray sigZ(Indices(d,d),dataz);//sigmaz
   mwArray sig0=identityMatrix(d);
   mwArray proj0=.5*(sigZ+sig0);
   mwArray Zf(Indices(2,d,d));
   for(int d1=0;d1<d;d1++){
     for(int d2=0;d2<d;d2++){
       Zf.setElement(sig0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
       Zf.setElement(proj0.getElement(Indices(d1,d2)),Indices(1,d1,d2));
     }
   }
   Zf.reshape(Indices(2,d*d));

   mpo.initLength(4*L);
   for(int k=0;k<4*L;k++){
     int Dl=(k==0)?1:2;
     int Dr=(k==4*L-1)?1:2;
     mwArray C(Indices(Dl,Dr,2));
     if(k<4*L-1)
       C.setElement(ONE_c,Indices(0,0,0)); // except last one
     if(k>0)
       C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // except first one
     // in which sites I include the local n operator depends on spin
     // what is the orbitalvalue of s corresponding to this k=4*n+2*alpha+s
     int j=k%4; // alpha*2+s
     //   int alpha=(j-j%2)/2;
     int s=j%2;int signS=(s==0)?1:-1;
     C.setElement(signS*ONE_c,Indices(0,Dr-1,1));
     C.reshape(Indices(Dl*Dr,2));C.multiplyRight(Zf);
     C.reshape(Indices(Dl,Dr,d,d));
     C.permute(Indices(3,1,4,2));
     mpo.setOp(k,new Operator(C),true);
   }
   

}

void TwoOrbitalFermiHubbardHamiltonian::getSpinCorrelatorMPO(MPO& mpo,int type) const {
  switch(type){
  case 1:
  case 2:
    getSpinCorrelatorMPOXY(mpo,type);
    break;
  case 3:
    complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
    mwArray sigZ(Indices(d,d),dataz);//sigmaz
    mwArray sig0=identityMatrix(d);
    mwArray Zf(Indices(2,d,d));
    for(int d1=0;d1<d;d1++){
      for(int d2=0;d2<d;d2++){
	Zf.setElement(sig0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
	Zf.setElement(sigZ.getElement(Indices(d1,d2)),Indices(1,d1,d2));
      }
    }
    Zf.reshape(Indices(2,d*d));
    mpo.initLength(4*L);
    int D=4; 
    for(int k=0;k<4*L;k++){
      int Dl=(k==0)?1:D;
      int Dr=(k==4*L-1)?1:D;
      mwArray C(Indices(Dl,Dr,2));
      if(k<4*L-1)
	C.setElement(ONE_c,Indices(0,0,0)); // except last one
      if(k>0)
	C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // except first one
      // in which sites I include the local n operator depends on spin
      // what is the orbitalvalue of s corresponding to this k=4*n+2*alpha+s
      int j=k%4; // alpha*2+s
      int n=(k-j)/4; //n
      int signS=(j%2==0?)1:-1;
      if(n%2==0){ // even n
	if(n>0) // first of pair
	  C.setElement(signS*.5*ONE_c,Indices(0,1,1));
	if(n<L-1) //second of pair
	  C.setElement(signS*.5*ONE_c,Indices(2,Dr-1,1));
	if(j!=0) // let it through
	  C.setElement(ONE_c,Indices(1,1,0));
      }
      else{
	if(n<L-1) // first of pair
	  C.setElement(signS*.5*ONE_c,Indices(0,2,1));
	if(n>0) //second of pair
	  C.setElement(signS*.5*ONE_c,Indices(1,Dr-1,1));
	if(j!=0) // let it through
	  C.setElement(ONE_c,Indices(2,2,0));
      }
      C.reshape(Indices(Dl*Dr,2));C.multiplyRight(Zf);
      C.reshape(Indices(Dl,Dr,d,d));
      C.permute(Indices(3,1,4,2));
      mpo.setOp(k,new Operator(C),true);
    }
    break; 
  }
}

void TwoOrbitalFermiHubbardHamiltonian::getSpinCorrelatorMPOXY(MPO& mpo,int type) const {
  mwArray sig0=identityMatrix(d);//identity
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP);//sigmaPlus
  mwArray sigM=Hconjugate(sigP); //sigmaMinus
  mwArray Zf(Indices(3,d,d));
  for(int d1=0;d1<d;d1++){
    for(int d2=0;d2<d;d2++){
      Zf.setElement(sig0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
      Zf.setElement(sigP.getElement(Indices(d1,d2)),Indices(1,d1,d2));
      Zf.setElement(sigM.getElement(Indices(d1,d2)),Indices(3,d1,d2));
    }
  }
  Zf.reshape(Indices(3,d*d));
  mpo.initLength(4*L);
  int D=8; 
  for(int k=0;k<4*L;k++){
    int Dl=(k==0)?1:D;
    int Dr=(k==4*L-1)?1:D;
    mwArray C(Indices(Dl,Dr,2));
    if(k<4*L-1)
      C.setElement(ONE_c,Indices(0,0,0)); // except last one
    if(k>0)
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // except first one
    int j=k%4; // alpha*2+s
    int n=(k-j)/4; //n
    //int signS=(j%2==0?)1:-1;
    int signTyp=(type==1)?-1:+1;
    if(n<L-1){ // first one
      if(j==0||j==2){
	C.setElement(ONE_c,Indices(0,1,1)); // sigP
	C.setElement(signTyp*ONE_c,Indices(0,2,2)); // sigM
      }
      else{ // j=1 and 3
	C.setElement(ONE_c,Indices(1,3+n%2,2)); // sigM
	C.setElement(ONE_c,Indices(2,3+n%2,1)); // sigP

      }
    }
    if(n>0){ // second one
      if(j==0||j==2){
	C.setElement(signTyp*ONE_c,Indices(3+n%2,5,1)); // sigP
	C.setElement(ONE_c,Indices(3+n%2,6,2)); // sigM
      }
      else{ // j=1 and 3
	C.setElement(ONE_c,Indices(5,Dr-1,2)); // sigM
	C.setElement(ONE_c,Indices(6,Dr-1,1)); // sigP
      }
    }
    if(k>0&&k<4*L-1){ // could improve
      if(j==0||j==1){
	C.setElement(ONE_c,Indices(4-n%2,4-n%2,0)); // even lets 4 through, odd lets 3
      }
      else{ // j=2,3
	C.setElement(ONE_c,Indices(3+n%2,3+n%2,0)); // even lets 3 through, odd lets 4
      }
    }
    // HERE!!!!

    C.reshape(Indices(Dl*Dr,3));C.multiplyRight(Zf);
    C.reshape(Indices(Dl,Dr,d,d));
    C.permute(Indices(3,1,4,2));
    mpo.setOp(k,new Operator(C),true);
  }
}

void TwoOrbitalFermiHubbardHamiltonian::getOrbitalSpinCorrelatorMPO(MPO& mpo) const {
    complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
    mwArray sigZ(Indices(d,d),dataz);//sigmaz
    mwArray sig0=identityMatrix(d);
    mwArray Zf(Indices(2,d,d));
    for(int d1=0;d1<d;d1++){
      for(int d2=0;d2<d;d2++){
	Zf.setElement(sig0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
	Zf.setElement(sigZ.getElement(Indices(d1,d2)),Indices(1,d1,d2));
      }
    }
    Zf.reshape(Indices(2,d*d));
    mpo.initLength(4*L);
    int D=4; 
    for(int k=0;k<4*L;k++){
      int Dl=(k==0)?1:D;
      int Dr=(k==4*L-1)?1:D;
      mwArray C(Indices(Dl,Dr,2));
      if(k<4*L-1)
	C.setElement(ONE_c,Indices(0,0,0)); // except last one
      if(k>0)
	C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // except first one
      // in which sites I include the local n operator depends on spin
      // what is the orbitalvalue of s corresponding to this k=4*n+2*alpha+s
      int j=k%4; // alpha*2+s
      int n=(k-j)/4; //n
      int signS=(j<2)?1:-1;
      if(n%2==0){ // even n
	if(n>0) // first of pair
	  C.setElement(signS*.5*ONE_c,Indices(0,1,1));
	if(n<L-1) //second of pair
	  C.setElement(signS*.5*ONE_c,Indices(2,Dr-1,1));
	if(j!=0) // let it through
	  C.setElement(ONE_c,Indices(1,1,0));
      }
      else{
	if(n<L-1) // first of pair
	  C.setElement(signS*.5*ONE_c,Indices(0,2,1));
	if(n>0) //second of pair
	  C.setElement(signS*.5*ONE_c,Indices(1,Dr-1,1));
	if(j!=0) // let it through
	  C.setElement(ONE_c,Indices(2,2,0));
      }
      C.reshape(Indices(Dl*Dr,2));C.multiplyRight(Zf);
      C.reshape(Indices(Dl,Dr,d,d));
      C.permute(Indices(3,1,4,2));
      mpo.setOp(k,new Operator(C),true);
    }
}


//   complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
//   mwArray sigZ(Indices(d,d),dataz);//sigmaz
//   mwArray sig0=identityMatrix(d);

//   mwArray nrUp;
//   constructOperatorProduct(nrUp,.5*(sigZ+sig0),sig0);
//   mwArray nrDown;
//   constructOperatorProduct(nrDown,sig0,.5*(sig0+sigZ));
//   mwArray aux=nrUp-nrDown; // TODO!!!! CHECK?
//   mwArray sig00=identityMatrix(d*d);
//   // // F_0^z= .5*(1+sigZ) (x) Id - Id (x) (1+tauZ)/2 

//   mwArray Zf(Indices(2,d*d,d*d));
//   for(int d1=0;d1<d*d;d1++){
//     for(int d2=0;d2<d*d;d2++){
//       Zf.setElement(sig00.getElement(Indices(d1,d2)),Indices(0,d1,d2));
//       Zf.setElement(aux.getElement(Indices(d1,d2)),Indices(1,d1,d2));
//     }
//   }
//   Zf.reshape(Indices(2,d*d*d*d));
  
//   mpo.initLength(L);
//   for(int k=0;k<L;k++){
//     if(k<=1||k==L-1){ // the different ones
//       int Dl=(k==0)?1:2;
//       int Dr=(k==L-1)?1:2;
//       mwArray C(Indices(Dl,Dr,2));
//       if(k<L-1)
// 	C.setElement(ONE_c,Indices(0,0,0)); // except last one
//       if(k>0)
// 	C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // not first one
//       C.setElement(ONE_c,Indices(0,Dr-1,1));
//       C.reshape(Indices(Dl*Dr,2));C.multiplyRight(Zf);
//       C.reshape(Indices(Dl,Dr,d*d,d*d));
//       C.permute(Indices(3,1,4,2));
//       mpo.setOp(k,new Operator(C),true); 
//     }
//     else{
//       mpo.setOp(k,&mpo.getOp(1),false); // could use a link
//     }
//   }  
// }


// void TwoOrbitalFermiHubbardHamiltonian::getExponentialMPOeven(MPO& expHe,complex_t delta) const {
//    getExponentialMPO(expHe,delta,true);
// }

// void TwoOrbitalFermiHubbardHamiltonian::getExponentialMPOodd(MPO& expHo,complex_t delta) const {
//   getExponentialMPO(expHo,delta,false);
// }


// void TwoOrbitalFermiHubbardHamiltonian::getExponentialMPO(MPO& expH,complex_t delta,bool even) const {
//   mwArray Ol,Or;
//   expH.initLength(L);
//   //cout<<"TwoOrbitalFermiHubbardHamiltonian::getExponentialMPO(delta="<<delta<<",even="<<even<<")"<<endl;
//   int kfirst=even?0:1;
//   int len=L;
//   int klast=(len-1)-(len-kfirst)%2; // last occupied by the loop
//   // The first exponential is always special (the rest are copies, as I put the U term on the left site)
//   bool repeat=false;int posOl(0),posOr(0);
//   for(int k=kfirst;k<klast;k+=2){
//     if(!repeat){
//       getTwoBodyTermExponential(Ol,Or,delta,k);
//       expH.setOp(k,new Operator(Ol),true);
//       expH.setOp(k+1,new Operator(Or),true);
//       // from now on, repeat these operators
//       repeat=true;
//       posOl=k;posOr=k+1;
//       //cout<<"Will repeat ops Ol("<<k<<") Or("<<k+1<<")"<<endl;
//       //cout<<"Set original operators on "<<k<<" and "<<k+1<<endl;
//     }
//     else{
//       expH.setOp(k,&expH.getOp(posOl),false);
//       expH.setOp(k+1,&expH.getOp(posOr),false);
//       //cout<<"Set repeated operators on "<<k<<" and "<<k+1<<endl;
//     }
//   }
//   // Fill in the edges with identity operators
//   if(kfirst==1){ // the first site is the identity
//     mwArray ident=identityMatrix(d*d);
//     ident.reshape(Indices(d*d,1,d*d,1)); // just dx1xdx1
//     expH.setOp(0,new Operator(ident),true);
//     //cout<<"Set pos "<<0<<" to Id"<<endl;
//   }
//   if(klast<len-1){ // sth to be filled at the end: the single body term!
//     getTwoBodyTermExponential(Ol,Or,delta,len-1);
//     expH.setOp(len-1,new Operator(Ol),true);
//     //    mwArray ident=identityMatrix(d*d);
//     //ident.reshape(Indices(d*d,1,d*d,1)); // just dx1xdx1
//     //expH.setOp(len-1,new Operator(ident),true);
//     //cout<<"Set pos "<<len-1<<" to Id"<<endl;
//   }
// }

// void TwoOrbitalFermiHubbardHamiltonian::getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,
// 						 int pos) const {
//   mwArray H12;
//   computeTwoBodyTerm(H12,pos);
//   // Now take the matrix exponential
//   mwArray expH;
//   //  if(pos==2) cout<<"h12("<<pos<<")="<<H12<<endl;
//   wrapper::expm(H12,expH,delta);
//   //putForMatlab(cout,expH,"expH12");
//   if(pos==L-1){// special case: singlebody term for the last site
//     Ol=expH;
//     Ol.reshape(Indices(d*d,1,d*d,1));
//     return;
//   }

//   //if(pos==2) cout<<"Exponential of the h12("<<pos<<")="<<expH<<endl;
//   int d2=d*d;
//   int d1=d2;
//   // Obtained with indices (i'j')(ij) => permute to (i'i)(j'j) for svd
//   expH.reshape(Indices(d1,d2,d1,d2));
//   expH.permute(Indices(1,3,2,4));
//   expH.reshape(Indices(d1*d1,d2*d2));
//   // And compute the SVD
//   mwArray S; // place for the singular values
//   int nr=0;
//   double tol=0.;
//   //cout<<"Before decomposition, expH="<<expH<<endl;
//   wrapper::svd(expH,tol,nr,Ol,S,Or);
//   //if(pos==2) cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<", S="<<S<<endl;
//   // redistribute the singular values
//   S=sqrt(S);
//   // TODO: check no NaN???
//   Ol.multiplyRight(S);
//   Or.multiplyLeft(S);
//   // reshape for operators (i'x 1 x i x alpha ))
//   Ol.reshape(Indices(d1,1,d1,-1));
//   Or.reshape(Indices(-1,d2,d2,1));
//   Or.permute(Indices(2,1,3,4));  
//   //cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<endl;
//   //putForMatlab(cout,Ol,"Ol");
//   //putForMatlab(cout,Or,"Or");
// }

// void TwoOrbitalFermiHubbardHamiltonian::computeTwoBodyTerm(mwArray& result,int pos) const{
// }


