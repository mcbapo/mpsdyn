#include "TwoOrbitalFermiHubbardHamiltonian_v0.h"
#include "Indices.h"

using namespace shrt;
using namespace std;

/** Notice this version conserves Ne and Ng independently, and does not seem to agree with the papers!!! */

// TwoOrbitalFermiHubbardHamiltonian::TwoOrbitalFermiHubbardHamiltonian(int L_,double t_,double U_,double V_,double Vex_,double mu_g_,double mu_e_):
//   d(2),L(L_),t(t_),U(U_),V(V_),Vex(Vex_),mu_g(mu_g_),mu_e(mu_e_),hamil(4*L){
//   initZ();
//   initHMPO();
// }

TwoOrbitalFermiHubbardHamiltonian::TwoOrbitalFermiHubbardHamiltonian(int L_,const vector<double>& t,const vector<double>& U,
								     double V_,double Vex_,
								     const vector<double>& mu_g,const vector<double>& mu_e,
								     int Ntot_,double penalty_):
  d(2),L(L_),hamil(4*L),V(V_),Vex(Vex_),Ntot(Ntot_),penNtot(penalty_),penNg(0.),penNe(0.),penSz(0.){
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

void TwoOrbitalFermiHubbardHamiltonian::setPenalties(int Ng_,double penaltyNg,int Ne_,double penaltyNe,int Sz_,double penaltySz,
						     int Ntot_,double penaltyNtot){
  Ntot=Ntot_;penNtot=penaltyNtot;
  Ng=Ng_;penNg=penaltyNg;
  Ne=Ne_;penNe=penaltyNe;
  if(penNe>0 &&penNg>0){
    if(penNtot>0) cout<<"WARNING: Ignoring penalty for Ntot=Ng+Ne!"<<endl;
    penNtot=0.;
  }
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


void TwoOrbitalFermiHubbardHamiltonian::initHMPO(){
  int D=15; // bond dimension of the MPO without penalty terms
  int Dg,De,Dz,Dn;
  if(penNtot!=0){
    Dn=D-1;D++; // take the last one for this and increase
    cout<<"Dn="<<Dn<<endl;
  }
  if(penNg!=0){
    Dg=D-1;D++; // take the last one for this and increase
    cout<<"Dg="<<Dg<<endl;
  }
  if(penNe!=0){
    De=D-1;D++; // take the last one for this and increase
    cout<<"De="<<De<<endl;
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

	// Vex term (Vex tilde??)
	if(p==1){
 	  C.setElement(ONE_c,Indices(1,9,2)); // sigM
 	  C.setElement(ONE_c,Indices(5,10,1)); // sigP
	}
	if(p==2){
 	  C.setElement(ONE_c,Indices(9,9,2)); // sigM
 	  C.setElement(ONE_c,Indices(10,10,1)); // sigP
	}
	if(p==3){
 	  C.setElement(-Vex*ONE_c,Indices(9,Dr-1,1)); // sigP
 	  C.setElement(-Vex*ONE_c,Indices(10,Dr-1,2)); // sigM
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
	if(penNg!=0){
	  if(k<4*L-1&&alpha==0){ // first of the pair
	    C.setElement(ONE_c,Indices(0,Dg,4)); // proj0
	  }
	  if(k>0&&alpha==0){ // second of the pair
	    C.setElement(2*penNg*ONE_c,Indices(Dg,Dr-1,4)); // proj0
	  }
	  if(k>0&&k<4*L-1) C.setElement(ONE_c,Indices(Dg,Dg,0)); // pass it
	}
	if(penNe!=0){
	  if(k<4*L-1&&alpha==1){ // first of the pair
	    C.setElement(ONE_c,Indices(0,De,4)); // proj0
	  }
	  if(k>0&&alpha==1){ // second of the pair
	    C.setElement(2*penNe*ONE_c,Indices(De,Dr-1,4)); // proj0
	  }
	  if(k>0&&k<4*L-1) C.setElement(ONE_c,Indices(De,De,0)); // pass it
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
	if(penNg!=0&&alpha==0) mucompl=mucompl+penNg*(1-2*Ng)*ONE_c;
	if(penNe!=0&&alpha==1) mucompl=mucompl+penNe*(1-2*Ne)*ONE_c;
	if(penSz!=0) mucompl=mucompl+penSz*(1-2*signS*Sz)*ONE_c;
	C.setElement(mucompl,Indices(0,Dr-1,4)); // proj0

	
	// Add also the constant part of penalty terms, if any of them is present
	double termLoc=0.;
	if(penNtot!=0) termLoc+=penNtot*Ntot*Ntot;
	if(penNg!=0) termLoc+=penNg*Ng*Ng;
	if(penNe!=0) termLoc+=penNe*Ne*Ne;
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
    for(int k=0;k<4*L;k++){ // k=(n,alpha,s)=4*n+2*alpha+s
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
      int signS=(j%2==0)?1:-1; // (-1)^s
      if(n%2==0){ // even n
	if(n<L-1) // first of pair
	  C.setElement(signS*.5*ONE_c,Indices(0,1,1));
	if(n>0){ //second of pair
	  C.setElement(signS*.5*ONE_c,Indices(2,Dr-1,1));
	  if(j!=3) // let it through from site before (except last j)
	    C.setElement(ONE_c,Indices(2,2,0));
	}		
	if(j!=0&&k<4*L-1) // let it through within the same site
	  C.setElement(ONE_c,Indices(1,1,0));
      }
      else{
	if(n<L-1) // first of pair
	  C.setElement(signS*.5*ONE_c,Indices(0,2,1));
	if(n>0){ //second of pair
	  C.setElement(signS*.5*ONE_c,Indices(1,Dr-1,1));
	  if(j!=3) // let it through from site before (except last j)
	    C.setElement(ONE_c,Indices(1,1,0));
	}
	if(j!=0&&k<4*L-1) // let it through
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
      Zf.setElement(sigM.getElement(Indices(d1,d2)),Indices(2,d1,d2));
    }
  }
  Zf.reshape(Indices(3,d*d));
  mpo.initLength(4*L);
  int D=8; 
  for(int k=0;k<4*L;k++){
    int Dl=(k==0)?1:D;
    int Dr=(k==4*L-1)?1:D;
    mwArray C(Indices(Dl,Dr,3));
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
    for(int k=0;k<4*L;k++){ // k(n,alpha,s)=4*n+2*alpha+s
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
      int signS=(j<2)?1:-1; // (-1)^alpha
      if(n%2==0){ // even n
	if(n<L-1) // first of pair
	  C.setElement(signS*.5*ONE_c,Indices(0,1,1));
	if(n>0){ //second of pair
	  C.setElement(signS*.5*ONE_c,Indices(2,Dr-1,1));
	  if(j!=3) // let it through from site before (except last j)
	    C.setElement(ONE_c,Indices(2,2,0));
	}
	if(j!=0) // let it through
	  C.setElement(ONE_c,Indices(1,1,0));
      }
      else{
	if(n<L-1) // first of pair
	  C.setElement(signS*.5*ONE_c,Indices(0,2,1));
	if(n>0){ //second of pair
	  C.setElement(signS*.5*ONE_c,Indices(1,Dr-1,1));
	  if(j!=3) // let it through from site before (except last j)
	    C.setElement(ONE_c,Indices(1,1,0));
	}
	if(j!=0&&k<4*L-1) // let it through
	  C.setElement(ONE_c,Indices(2,2,0));
      }
      C.reshape(Indices(Dl*Dr,2));C.multiplyRight(Zf);
      C.reshape(Indices(Dl,Dr,d,d));
      C.permute(Indices(3,1,4,2));
      mpo.setOp(k,new Operator(C),true);
    }
}




