#include "TwoOrbitalFermiHubbardHamiltonianBlock.h"
#include "Indices.h"

using namespace shrt;
using namespace std;

/** Notice this version conserves Ne and Ng independently, and does not seem to agree with the papers!!! */

// TwoOrbitalFermiHubbardHamiltonianBlock::TwoOrbitalFermiHubbardHamiltonianBlock(int L_,double t_,double U_,double V_,double Vex_,double mu_g_,double mu_e_):
//   d(2),L(L_),t(t_),U(U_),V(V_),Vex(Vex_),mu_g(mu_g_),mu_e(mu_e_),hamil(4*L){
//   initZ();
//   initHMPO();
// }

TwoOrbitalFermiHubbardHamiltonianBlock::TwoOrbitalFermiHubbardHamiltonianBlock(int L_,const vector<double>& t,const vector<double>& U,
								     double V_,double Vex_,
								     const vector<double>& mu_g,const vector<double>& mu_e,
								     int Ntot_,double penalty_):
  dB(16),L(L_),hamil(L),V(V_),Vex(Vex_),Ntot(Ntot_),penNtot(penalty_),penNg(0.),penNe(0.),penSz(0.){
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
    }else{ cout<<"Invalid input t for TwoOrbitalFermiHubbardHamiltonianBlock"<<endl;exit(1);}
  }
  if(U.size()==2){
    Ug=U[0];Ue=U[1];
  }
  else{
    if(U.size()==1){
      Ug=U[0];Ue=Ug;
    }else{ cout<<"Invalid input U for TwoOrbitalFermiHubbardHamiltonianBlock"<<endl;exit(1);}
  }
  if(mu_g.size()==2){
    mu_g0=mu_g[0];mu_g1=mu_g[1];
  }
  else{
    if(mu_g.size()==1){
      mu_g0=mu_g[0];mu_g1=mu_g0;
    }else{ cout<<"Invalid input mu_g for TwoOrbitalFermiHubbardHamiltonianBlock"<<endl;exit(1);}
  }
  if(mu_e.size()==2){
    mu_e0=mu_e[0];mu_e1=mu_e[1];
  }
  else{
    if(mu_e.size()==1){
      mu_e0=mu_e[0];mu_e1=mu_e0;
    }else{ cout<<"Invalid input mu_e for TwoOrbitalFermiHubbardHamiltonianBlock"<<endl;exit(1);}
  }
  initLocalOps();
  initZ();
  initHMPO();
}

TwoOrbitalFermiHubbardHamiltonianBlock::~TwoOrbitalFermiHubbardHamiltonianBlock(){};

void TwoOrbitalFermiHubbardHamiltonianBlock::setPenalties(int Ng_,double penaltyNg,int Ne_,double penaltyNe,int Sz_,double penaltySz,
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


// auxiliary function to transform a single site op into a local
// operator on the 4-site, at the specified site, with or without Z
// string (which can be left or right)

void TwoOrbitalFermiHubbardHamiltonianBlock::promoteLocalOp(mwArray& result,const mwArray& localOp,int pos,bool Zs,bool Zsright){
  //  cout<<"promoteLocalOp(pos="<<pos<<", with stringZ="<<Zs<<", to right="<<Zsright<<")"<<endl;
  int d=2;
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz
  if(pos<0||pos>=4){
    cout<<"ERROR: wrong position in local op!"<<endl;
    exit(1);
  }
  static vector<mwArray> Zstr;
  static bool init=false;
  if(!init){
    mwArray ZZ,ZZZ;
    Zstr.push_back(sigZ);
    constructOperatorProduct(ZZ,sigZ,sigZ);
    constructOperatorProduct(ZZZ,ZZ,sigZ);
    Zstr.push_back(ZZ);
    Zstr.push_back(ZZZ);
    init=true;
  }

  mwArray aux;
  if(pos<3){// right part
    if(Zs&&Zsright){ // on the right there is a string
      constructOperatorProduct(aux,localOp,Zstr[3-pos-1]);
      //cout<<"\t string of Zs to right: "<<Zstr[3-pos-1].getDimensions()<<endl;
    }
    else{
      constructOperatorProduct(aux,localOp,identityMatrix(pow(d,3-pos)));
      //cout<<"\t identity to right: "<<pow(d,3-pos)<<endl;
    }
  }
  else{ // op on the last site
    aux=localOp;
  }  
  if(pos>0){ // left part
    if(Zs&&!Zsright){ // Zs to the left
      constructOperatorProduct(result,Zstr[pos-1],aux);
      //cout<<"\t string of Zs to left: "<<Zstr[pos-1].getDimensions()<<endl;
    }
    else{ // only identities
      constructOperatorProduct(result,identityMatrix(pow(d,pos)),aux);
      //cout<<"\t identity to left: "<<pow(d,pos)<<endl;
    }
  }
  else
    result=aux;
  //  cout<<"At the end of promoteLocalOp, result:"<<result.getDimensions()<<endl;
}

void TwoOrbitalFermiHubbardHamiltonianBlock::initLocalOps(){
  //  cout<<"initLocalOps()"<<endl;
  int d=2;
  mwArray sig0=identityMatrix(d);//identity
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP);//sigmaPlus
  mwArray sigM=Hconjugate(sigP); //sigmaMinus

  mwArray sigP0,sigP1,sigP2,sigP3; // sigPlus on each site, with no Zs
  promoteLocalOp(sigP0,sigP,0,0);
  promoteLocalOp(sigP1,sigP,1,0);
  promoteLocalOp(sigP2,sigP,2,0);
  promoteLocalOp(sigP3,sigP,3,0);

  // sigMinus on each site, with no Zs
  mwArray sigM0=Hconjugate(sigP0);
  mwArray sigM1=Hconjugate(sigP1);
  mwArray sigM2=Hconjugate(sigP2);
  mwArray sigM3=Hconjugate(sigP3);

  opSx=-0.5*I_c*(sigP0*sigM1+sigP2*sigM3);
  opSx=opSx+Hconjugate(opSx);

  opSy=-.5*(sigP0*sigM1+sigP2*sigM3);
  opSy=opSy+Hconjugate(opSy);  
  
}

void TwoOrbitalFermiHubbardHamiltonianBlock::initZ(){
  cout<<"initZ()"<<endl;
  // Basic pieces: identity plus 16 operators for hopping (+,- for each mode in each
  // site, but also left and right, due to different Z strings)
  // plus the single site operators
  int d=2;
  nrOps=1; // identity
  mwArray sig0=identityMatrix(d);//identity
  complex_t dataP[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
  mwArray sigP(Indices(d,d),dataP);//sigmaPlus
  mwArray sigM=Hconjugate(sigP); //sigmaMinus
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz
  mwArray projZ0=.5*(sig0+sigZ);
  mwArray sig0b=identityMatrix(dB);//identity in block

  // Construct the basic operators by composing the others
  mwArray sigP0l,sigP1l,sigP2l,sigP3l; // sigPlus on each site, with Zs to the right
  //  promoteLocalOp(rmwArray& result,const mwArray& localOp,int pos,bool Zs=0,bool Zsright=true){
  promoteLocalOp(sigP0l,sigP,0,1,1);
  promoteLocalOp(sigP1l,sigP,1,1,1);
  promoteLocalOp(sigP2l,sigP,2,1,1);
  promoteLocalOp(sigP3l,sigP,3,1,1);

  mwArray sigP0r,sigP1r,sigP2r,sigP3r; // sigPlus on each site, with Zs to the left
  promoteLocalOp(sigP0r,sigP,0,1,0);
  promoteLocalOp(sigP1r,sigP,1,1,0);
  promoteLocalOp(sigP2r,sigP,2,1,0);
  promoteLocalOp(sigP3r,sigP,3,1,0);

  // the other part of hopping is adjoint of these
  mwArray sigM0l=Hconjugate(sigP0l);
  mwArray sigM1l=Hconjugate(sigP1l);
  mwArray sigM2l=Hconjugate(sigP2l);
  mwArray sigM3l=Hconjugate(sigP3l);
  mwArray sigM0r=Hconjugate(sigP0r);
  mwArray sigM1r=Hconjugate(sigP1r);
  mwArray sigM2r=Hconjugate(sigP2r);
  mwArray sigM3r=Hconjugate(sigP3r);

  
  nrOps+=4*4;
  
  // now the local operators: n==projZ0

  mwArray n0,n1,n2,n3;
  promoteLocalOp(n0,projZ0,0,0);
  promoteLocalOp(n1,projZ0,1,0);
  promoteLocalOp(n2,projZ0,2,0);
  promoteLocalOp(n3,projZ0,3,0);

  nrOps+=4; // chemical potential
  
  // For convenience, I keep the local operators:
  opNg=n0+n1;
  opNe=n2+n3;
  opSz=.5*(n0+n2-n1-n3);
  opTz=.5*(n0+n1-n2-n3);
  // mwArray ng(n0+n1);
  // mwArray ne(n2+n3);
  // mwArray Sz(n0+n2-n1-n3);
  // mwArray Tz(n0+n1-n2-n3);
  mwArray ntot=opNg+opNe;
  // and their squares
  mwArray ng2(opNg*opNg);
  mwArray ne2(opNe*opNe);
  mwArray Sz2(opSz*opSz);
  mwArray Tz2(opTz*opTz);
  mwArray ntot2(ntot*ntot);
  
  // Now in H I have also:
  // n_k_up n_k_down (n0 n1 + n2 n3), term with U
  mwArray termUg=n0*n1;
  mwArray termUe=n2*n3;
  // V term ne_s ng_s' (n0n3+n1n2)
  mwArray Vterm=n0*n3+n1*n2;
  // V-Vex term nes ng's (n0n2+n1n3)
  mwArray VVexterm=n0*n2+n1*n3;
  // Vex term cgs^+ ces'^+ cgs' ces (sP sM sM sP + h.c.)
  mwArray Vexterm;
  {mwArray aux,aux2;
    constructOperatorProduct(aux,sigP,sigM);
    constructOperatorProduct(aux2,aux,sigM);
    constructOperatorProduct(aux,aux2,sigP);
    Vexterm=aux+Hconjugate(aux);
  }

  nrOps+=5; // local terms in H

  // for the penalty terms, potentially, also ne, ng, ntot, Sz (not Tz) as individual ops, but also their squares
  nrOps+=2*4; 
  

  
  // Fill in the operators in the Z array, to be contracted with C in
  // the construction of the MPO. The order is Id,sigPlusL,sigMinusL,sigPlusR,sigMinusR,local ops
  Z=mwArray(Indices(nrOps,dB,dB));
  for(int d1=0;d1<dB;d1++){
    for(int d2=0;d2<dB;d2++){
	Z.setElement(sig0b.getElement(Indices(d1,d2)),Indices(0,d1,d2));
	Z.setElement(sigP0l.getElement(Indices(d1,d2)),Indices(1+0,d1,d2));
	Z.setElement(sigP1l.getElement(Indices(d1,d2)),Indices(1+1,d1,d2));
	Z.setElement(sigP2l.getElement(Indices(d1,d2)),Indices(1+2,d1,d2));
	Z.setElement(sigP3l.getElement(Indices(d1,d2)),Indices(1+3,d1,d2));

	Z.setElement(sigM0l.getElement(Indices(d1,d2)),Indices(5+0,d1,d2));
	Z.setElement(sigM1l.getElement(Indices(d1,d2)),Indices(5+1,d1,d2));
	Z.setElement(sigM2l.getElement(Indices(d1,d2)),Indices(5+2,d1,d2));
	Z.setElement(sigM3l.getElement(Indices(d1,d2)),Indices(5+3,d1,d2));
	
	Z.setElement(sigP0r.getElement(Indices(d1,d2)),Indices(9+0,d1,d2));
	Z.setElement(sigP1r.getElement(Indices(d1,d2)),Indices(9+1,d1,d2));
	Z.setElement(sigP2r.getElement(Indices(d1,d2)),Indices(9+2,d1,d2));
	Z.setElement(sigP3r.getElement(Indices(d1,d2)),Indices(9+3,d1,d2));

	Z.setElement(sigM0r.getElement(Indices(d1,d2)),Indices(13+0,d1,d2));
	Z.setElement(sigM1r.getElement(Indices(d1,d2)),Indices(13+1,d1,d2));
	Z.setElement(sigM2r.getElement(Indices(d1,d2)),Indices(13+2,d1,d2));
	Z.setElement(sigM3r.getElement(Indices(d1,d2)),Indices(13+3,d1,d2));


	// Now all local ops
	Z.setElement(termUg.getElement(Indices(d1,d2)),Indices(17,d1,d2)); // Ug
	Z.setElement(termUe.getElement(Indices(d1,d2)),Indices(18,d1,d2)); // Ue
	Z.setElement(Vterm.getElement(Indices(d1,d2)),Indices(19,d1,d2)); // V
	Z.setElement(VVexterm.getElement(Indices(d1,d2)),Indices(20,d1,d2)); // V-Vex delta(ss')
	Z.setElement(Vexterm.getElement(Indices(d1,d2)),Indices(21,d1,d2)); // Vex

	// Projectors
	Z.setElement(n0.getElement(Indices(d1,d2)),Indices(22+0,d1,d2)); // n0
	Z.setElement(n1.getElement(Indices(d1,d2)),Indices(22+1,d1,d2)); // n1
	Z.setElement(n2.getElement(Indices(d1,d2)),Indices(22+2,d1,d2)); // n2
	Z.setElement(n3.getElement(Indices(d1,d2)),Indices(22+3,d1,d2)); // n3

	// Terms for penalty
	Z.setElement(opNg.getElement(Indices(d1,d2)),Indices(26,d1,d2)); // ng
	Z.setElement(opNe.getElement(Indices(d1,d2)),Indices(27,d1,d2)); // ne
	Z.setElement(2*opSz.getElement(Indices(d1,d2)),Indices(28,d1,d2)); // Sz
	Z.setElement(ntot.getElement(Indices(d1,d2)),Indices(29,d1,d2)); // ntot

	Z.setElement(ng2.getElement(Indices(d1,d2)),Indices(30,d1,d2)); // ng^2
	Z.setElement(ne2.getElement(Indices(d1,d2)),Indices(31,d1,d2)); // ne^2
	Z.setElement(4*Sz2.getElement(Indices(d1,d2)),Indices(32,d1,d2)); // Sz^2
	Z.setElement(ntot2.getElement(Indices(d1,d2)),Indices(33,d1,d2)); // ntot^2

    }
  }
  Z.reshape(Indices(nrOps,dB*dB));
  cout<<"Prepared Z:"<<Z.getDimensions()<<endl;
}


void TwoOrbitalFermiHubbardHamiltonianBlock::initHMPO(){
  cout<<"initHMPO()"<<endl;
  int D=10; // bond dimension of the MPO without penalty terms
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

  // Second part of hopping carries the subindex (alpha,s)
  vector<complex_t> tcompl;
  tcompl.push_back(tg0*ONE_c);
  tcompl.push_back(tg1*ONE_c);
  tcompl.push_back(te0*ONE_c);
  tcompl.push_back(te1*ONE_c);

  vector<complex_t> mucompl;
  mucompl.push_back(-mu_g0*ONE_c);
  mucompl.push_back(-mu_g1*ONE_c);
  mucompl.push_back(-mu_e0*ONE_c);
  mucompl.push_back(-mu_e1*ONE_c);

  for(int k=0;k<L;k++){
    int Dl=(k==0)?1:D;
    int Dr=(k==L-1)?1:D;
    mwArray C(Indices(Dl,Dr,nrOps)); 
    if(k<L-1)
      C.setElement(ONE_c,Indices(0,0,0)); // nothing yet, only Id
    if(k>0)
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // after everything
    if(k<L-1){
      // First part of hopping: 8 possibilities (sigPl or sigMl)
      for(int p=0;p<4;p++)
	C.setElement(ONE_c,Indices(0,p+1,p+1)); // sigPl 
      for(int p=0;p<4;p++)
	C.setElement(ONE_c,Indices(0,p+5,p+5)); // sigMl 
    }
    if(k>0){ // second term of hopping
      for(int p=0;p<4;p++){
       C.setElement(tcompl[p],Indices(p+1,Dr-1,p+13)); // sigMr 
       C.setElement(tcompl[p],Indices(p+5,Dr-1,p+9)); // sigPr 
      }
      if(tge0!=0.){
	C.setElement(tge0*ONE_c,Indices(1,Dr-1,15)); // sigP0l sigM2r
	C.setElement(tge0*ONE_c,Indices(5,Dr-1,11)); // sigM0l sigP2r
      }
      if(tge1!=0.){
	C.setElement(tge1*ONE_c,Indices(2,Dr-1,16)); // sigP1l sigM3r
	C.setElement(tge1*ONE_c,Indices(6,Dr-1,12)); // sigM1l sigP3r
      }
    }
    // Local terms (in all sites)
    // U
    C.setElement(Ug*ONE_c,Indices(0,Dr-1,17)); // U* sum nxup nxdown
    C.setElement(Ue*ONE_c,Indices(0,Dr-1,18)); // U* sum nxup nxdown
    // V
    C.setElement(V*ONE_c,Indices(0,Dr-1,19)); // V nes ngs'
    C.setElement((V-Vex)*ONE_c,Indices(0,Dr-1,20)); // V-Vex for s=s'
    // Vex
    C.setElement((-Vex)*ONE_c,Indices(0,Dr-1,21)); // Vex 

    // chemical potentials
    for(int p=0;p<3;p++)
      C.setElement(mucompl[p],Indices(0,Dr-1,22+p)); // n_p

    // constant shift (may ignore)
    complex_t termLoc(ZERO_c);
    // Terms from penalties
    if(penNg!=0){
      C.setElement(penNg*ONE_c,Indices(0,Dr-1,30)); // ng^2
      if(k<L-1) C.setElement(2*penNg*ONE_c,Indices(0,Dg,26)); // ng (first)
      if(k>0&&k<L-1) C.setElement(ONE_c,Indices(Dg,Dg,0)); // intermediate
      if(k>0) C.setElement(ONE_c,Indices(Dg,Dr-1,26)); // ng (second)
      C.setElement(-2*penNg*Ng*ONE_c,Indices(0,Dr-1,26)); // linear term
      termLoc+=penNg*Ng*Ng*ONE_c;
    }
    if(penNe!=0){
      C.setElement(penNe*ONE_c,Indices(0,Dr-1,31)); // ne^2
      if(k<L-1) C.setElement(2*penNe*ONE_c,Indices(0,De,27)); // ne (first)
      if(k>0&&k<L-1) C.setElement(ONE_c,Indices(De,De,0)); // intermediate
      if(k>0) C.setElement(ONE_c,Indices(De,Dr-1,27)); // ne (second)
      C.setElement(-2*penNe*Ne*ONE_c,Indices(0,Dr-1,27)); // linear term
      termLoc+=penNe*Ne*Ne*ONE_c;
    }
    if(penSz!=0){
      C.setElement(penSz*ONE_c,Indices(0,Dr-1,32)); // Sz^2
      if(k<L-1) C.setElement(2*penSz*ONE_c,Indices(0,Dz,28)); // ne (first)
      if(k>0&&k<L-1) C.setElement(ONE_c,Indices(Dz,Dz,0)); // intermediate
      if(k>0) C.setElement(ONE_c,Indices(Dz,Dr-1,28)); // ne (second)
      C.setElement(-2*penSz*Sz*ONE_c,Indices(0,Dr-1,28)); // linear term
      termLoc+=penSz*Sz*Sz*ONE_c;
    }
    if(penNtot!=0){
      C.setElement(penNtot*ONE_c,Indices(0,Dr-1,33)); // n^2
      if(k<L-1) C.setElement(2*penNtot*ONE_c,Indices(0,Dn,29)); // n (first)
      if(k>0&&k<L-1) C.setElement(ONE_c,Indices(Dn,Dn,0)); // intermediate
      if(k>0) C.setElement(ONE_c,Indices(Dn,Dr-1,29)); // n (second)
      C.setElement(-2*penNtot*Ntot*ONE_c,Indices(0,Dr-1,29)); // linear term
      termLoc+=penNtot*Ntot*Ntot*ONE_c;
    }

    C.setElement((1./L)*termLoc,Indices(0,Dr-1,0)); // Indentity (could skip)
    
    // Now set the op in the MPO
    C.reshape(Indices(Dl*Dr,nrOps));
    C.multiplyRight(Z);
    C.reshape(Indices(Dl,Dr,dB,dB));
    C.permute(Indices(3,1,4,2));
    //cout<<"Setting operator of Hmpo at site "<<k<<"("<<n<<","<<alpha<<","<<s<<")"<<endl;
    hamil.setOp(k,new Operator(C),true);    

  } // end loop over sites
}


void TwoOrbitalFermiHubbardHamiltonianBlock::getSingleSiteMPO(MPO& mpo,const mwArray& opN) const{
   mwArray sig0=identityMatrix(dB);
   mwArray Zf(Indices(2,dB,dB));
   for(int d1=0;d1<dB;d1++){
     for(int d2=0;d2<dB;d2++){
       Zf.setElement(sig0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
       Zf.setElement(opN.getElement(Indices(d1,d2)),Indices(1,d1,d2));
     }
   }
   Zf.reshape(Indices(2,dB*dB));

   mpo.initLength(L);
   for(int k=0;k<L;k++){
     int Dl=(k==0)?1:2;
     int Dr=(k==L-1)?1:2;
     mwArray C(Indices(Dl,Dr,2));
     if(k<L-1)
       C.setElement(ONE_c,Indices(0,0,0)); // except last one
     if(k>0)
       C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // except first one
     C.setElement(ONE_c,Indices(0,Dr-1,1));
     C.reshape(Indices(Dl*Dr,2));C.multiplyRight(Zf);
     C.reshape(Indices(Dl,Dr,dB,dB));
     C.permute(Indices(3,1,4,2));
     if(k<=1||k==L-1)
       mpo.setOp(k,new Operator(C),true);
     else{
       mpo.setOp(k,&mpo.getOp(1),false);
     }
   }   
}

void TwoOrbitalFermiHubbardHamiltonianBlock::getNNCorrelationMPO(MPO& mpo,const mwArray& opA,const mwArray& opB) const{
   mwArray sig0=identityMatrix(dB);
   mwArray Zf(Indices(3,dB,dB));
   for(int d1=0;d1<dB;d1++){
     for(int d2=0;d2<dB;d2++){
       Zf.setElement(sig0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
       Zf.setElement(opA.getElement(Indices(d1,d2)),Indices(1,d1,d2));
       Zf.setElement(opB.getElement(Indices(d1,d2)),Indices(2,d1,d2));
     }
   }
   Zf.reshape(Indices(3,dB*dB));

   mpo.initLength(L);
   for(int k=0;k<L;k++){
     int Dl=(k==0)?1:3;
     int Dr=(k==L-1)?1:3;
     mwArray C(Indices(Dl,Dr,3));
     if(k<L-1)
       C.setElement(ONE_c,Indices(0,0,0)); // except last one
     if(k>0)
       C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // except first one
     if(k<L-1) // first of pair
       C.setElement(ONE_c,Indices(0,1,1));
     if(k>0) // second of pair
       C.setElement(ONE_c,Indices(1,Dr-1,2));
     C.reshape(Indices(Dl*Dr,3));C.multiplyRight(Zf);
     C.reshape(Indices(Dl,Dr,dB,dB));
     C.permute(Indices(3,1,4,2));
     if(k<=1||k==L-1)
       mpo.setOp(k,new Operator(C),true);
     else
       mpo.setOp(k,&mpo.getOp(1),false);
   }   
}


void TwoOrbitalFermiHubbardHamiltonianBlock::getNumberMPO(MPO& mpo,int type) const {
   mwArray sig0=identityMatrix(dB);
   mwArray opN;
   switch(type){
   case 0:
     opN=opNg+opNe;
     break;
   case 1:
     opN=opNg;
     break;
   case 2:
     opN=opNe;
     break;
   default:
     cout<<"Wrong type of nr operator in TwoOrbitalFermiHubbardHamiltonianBlock::getNumberMPO"<<endl;
     exit(1);
     break;
   }
   getSingleSiteMPO(mpo,opN);
}

void TwoOrbitalFermiHubbardHamiltonianBlock::getFermionNrMPO(MPO& mpo) const {
  getNumberMPO(mpo,0);
}

void TwoOrbitalFermiHubbardHamiltonianBlock::getFermionNrOrbitalMPO(MPO& mpo,bool orbital0) const{
  int type=orbital0?1:2;
  getNumberMPO(mpo,type);
}

void TwoOrbitalFermiHubbardHamiltonianBlock::getDoubleOccupancyOrbitalMPO(MPO& mpo,bool orbital0) const{
   mwArray sig0=identityMatrix(dB);
   mwArray projOrb=orbital0==1?opNg:opNe;
   getNNCorrelationMPO(mpo,projOrb,projOrb);
}


void TwoOrbitalFermiHubbardHamiltonianBlock::getTotalSzMPO(MPO& mpo) const {
  getSingleSiteMPO(mpo,opSz);
}


void TwoOrbitalFermiHubbardHamiltonianBlock::getSpinCorrelatorMPO(MPO& mpo,int type) const {
  switch(type){
  case 1:
    getNNCorrelationMPO(mpo,opSx,opSx);
    break;
  case 2:
    getNNCorrelationMPO(mpo,opSy,opSy);
    break;
  case 3:
    getNNCorrelationMPO(mpo,opSz,opSz);
    break;
  }
}

void TwoOrbitalFermiHubbardHamiltonianBlock::getOrbitalSpinCorrelatorMPO(MPO& mpo) const {
  getNNCorrelationMPO(mpo,opTz,opTz);
}


// methods for time evolution, to be done

void TwoOrbitalFermiHubbardHamiltonianBlock::getExponentialMPOeven(MPO& expHe,complex_t delta) const {};
void TwoOrbitalFermiHubbardHamiltonianBlock::getExponentialMPOodd(MPO& expHo,complex_t delta) const {};

void TwoOrbitalFermiHubbardHamiltonianBlock::getDoubleExponentialMPOeven(MPO& expHe,complex_t delta) const {};
void TwoOrbitalFermiHubbardHamiltonianBlock::getDoubleExponentialMPOodd(MPO& expHo,complex_t delta) const {};

void TwoOrbitalFermiHubbardHamiltonianBlock::getExtendedExponentialMPOeven(MPO& expHe,complex_t delta) const {};
void TwoOrbitalFermiHubbardHamiltonianBlock::getExtendedExponentialMPOodd(MPO& expHo,complex_t delta) const {};




