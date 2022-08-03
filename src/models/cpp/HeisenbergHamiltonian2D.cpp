#include "HeisenbergHamiltonian2D.h"
#include "Indices.h"
#include "DoubleOperator.h"

using namespace std;
using namespace shrt;


HeisenbergHamiltonian2D::HeisenbergHamiltonian2D(int Lx_,int Ly_,int nL_,double Jx_,double Jy_,double Jz_,
						 double hx_,double hy_,double hz_,double offset):hamil(Lx_*Ly_/nL_){
  d=2;
  Jx=Jx_;Jy=Jy_;Jz=Jz_;
  hx=hx_;hy=hy_;hz=hz_;
  Lx=Lx_;Ly=Ly_;
  nL=nL_;

  if(Ly%nL!=0){
    cout<<"ERROR: Not yet ready for sites with different nr of spins!"<<endl;
    exit(1);
  }

  legNr=Ly/nL;
  cout<<"Creating HH ladder for ("<<Lx<<","<<Ly<<") with "<<legNr<<" legs and "<<nL<<" spins per site"<<endl;


  hamil.initLength(legNr*Lx);
  
  initZ();
  initHMPO(offset);
}


HeisenbergHamiltonian2D::~HeisenbergHamiltonian2D(){hamil.clear();Z.clear();}

void HeisenbergHamiltonian2D::getExponentialMPOevenH(MPO& expHe,complex_t delta) const {
 getExponentialMPOx(expHe,delta,true);
}
void HeisenbergHamiltonian2D::getExponentialMPOoddH(MPO& expHo,complex_t delta) const {
 getExponentialMPOx(expHo,delta,false);
}
void HeisenbergHamiltonian2D::getExponentialMPOevenV(MPO& expHe,complex_t delta) const {
 getExponentialMPOy(expHe,delta,true);
}
void HeisenbergHamiltonian2D::getExponentialMPOoddV(MPO& expHo,complex_t delta) const {
 getExponentialMPOy(expHo,delta,false);
}

void HeisenbergHamiltonian2D::getDoubleExponentialMPOevenH(MPO& expHe,complex_t delta) const {
 getDoubleExponentialMPOx(expHe,delta,true);
}
void HeisenbergHamiltonian2D::getDoubleExponentialMPOoddH(MPO& expHo,complex_t delta) const {
 getDoubleExponentialMPOx(expHo,delta,false);
}
void HeisenbergHamiltonian2D::getDoubleExponentialMPOevenV(MPO& expHe,complex_t delta) const {
 getDoubleExponentialMPOy(expHe,delta,true);
}
void HeisenbergHamiltonian2D::getDoubleExponentialMPOoddV(MPO& expHo,complex_t delta) const {
 getDoubleExponentialMPOy(expHo,delta,false);
}


void HeisenbergHamiltonian2D::getExtendedExponentialMPOevenH(MPO& expHe,complex_t delta) const {
 getExtendedExponentialMPOx(expHe,delta,true);
}
void HeisenbergHamiltonian2D::getExtendedExponentialMPOoddH(MPO& expHo,complex_t delta) const {
 getExtendedExponentialMPOx(expHo,delta,false);
}
void HeisenbergHamiltonian2D::getExtendedExponentialMPOevenV(MPO& expHe,complex_t delta) const {
 getExtendedExponentialMPOy(expHe,delta,true);
}
void HeisenbergHamiltonian2D::getExtendedExponentialMPOoddV(MPO& expHo,complex_t delta) const {
 getExtendedExponentialMPOy(expHo,delta,false);
}


void HeisenbergHamiltonian2D::getTwoBodyTermExponentialx(mwArray& Ol,mwArray& Or,complex_t delta,
						      int posX) const {
  mwArray H12;
  computeTwoBodyTermHorizontal(H12,posX);
  getTwoBodyTermExponential(Ol,Or,H12,delta);
}

void HeisenbergHamiltonian2D::getTwoBodyTermExponentialy(mwArray& Ol,mwArray& Or,complex_t delta,
						      int posY) const {
  mwArray H12;
  computeTwoBodyTermVertical(H12,posY);
  getTwoBodyTermExponential(Ol,Or,H12,delta);
}

void HeisenbergHamiltonian2D::getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,
						      const mwArray& H12,complex_t delta) const {
  int dloc=pow(d,nL); 
  // Now take the matrix exponential
  mwArray expH;
  wrapper::expm(H12,expH,delta);
  //putForMatlab(cout,expH,"expH12");
  //cout<<"Exponential of the h12="<<expH<<endl;

  // Obtained with indices (i'j')(ij) => permute to (i'i)(j'j) for svd
  expH.reshape(Indices(dloc,dloc,dloc,dloc));
  expH.permute(Indices(1,3,2,4));
  expH.reshape(Indices(dloc*dloc,dloc*dloc));
  // And compute the SVD
  mwArray S; // place for the singular values
  int nr=0;
  double tol=0.;
  //cout<<"Before decomposition, expH="<<expH<<endl;
  wrapper::svd(expH,tol,nr,Ol,S,Or);
  //cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<", S="<<S<<endl;
  // redistribute the singular values
  S=sqrt(S);
  // TODO: check no NaN???
  Ol.multiplyRight(S);
  Or.multiplyLeft(S);
  // reshape for operators (i'x 1 x i x alpha ))
  Ol.reshape(Indices(dloc,1,dloc,-1));
  Or.reshape(Indices(-1,dloc,dloc,1));
  Or.permute(Indices(2,1,3,4));  
  //cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<endl;
  //putForMatlab(cout,Ol,"Ol");
  //putForMatlab(cout,Or,"Or");
}

void HeisenbergHamiltonian2D::computeTwoBodyTermHorizontal(mwArray& result,int pos) const{
  if(pos>=Lx-1){
    cout<<"ERROR: HeisenbergHamiltonian2D::computeTwoBodyTermHorizontal for rung "<<pos<<" when Lx="<<Lx<<endl;
    exit(212);
  }
  // the local terms are included as 1/4 in general, and in edges 1/4 more, except if there is only one leg, in which case the factor is 1/2. 
  double factorY=(legNr==1&&nL==1)?.5:.25;
  double factor1=pos==0?2.:1.;factor1=factor1*factorY;
  double factor2=pos==Lx-2?2.:1.;factor2=factor2*factorY;
  mwArray Hx;
  computeTwoBodyTerm(Hx,Jx,Jy,Jz,factor1*hx,factor2*hx,factor1*hy,factor2*hy,factor1*hz,factor2*hz);
  //  cout<<"HeisenbergHamiltonian2D::computeTwoBodyTermHorizontal rung="<<pos<<", factor1="<<factor1<<", factor2="<<factor2<<endl;
  // I need to tensor nL of those
  result=Hx; // d*d x d*d
  int dtot=d;
  for(int k=1;k<nL;k++){
    mwArray tmp(result);
    constructOperatorProduct(result,tmp,identityMatrix(d*d));
    constructOperatorProduct(tmp,identityMatrix(dtot*dtot),Hx);
    result=result+tmp;
    result.reshape(Indices(dtot,dtot,d,d,dtot,dtot,d,d));
    result.permute(Indices(1,3,2,4,5,7,6,8));    
    dtot*=d;
    result.reshape(Indices(dtot*dtot,dtot*dtot));
  }
  if(legNr==1&&nL>1){
    cout<<"Special case: single leg, but multispin sites=> horizontal terms need"
	<<" to include internal vertical pairs as local "<<endl;
    mwArray term1=0.*identityMatrix(d); // term2 is the same
    for(int l=0;l<nL-1;l++){
      double factorX=.25;
      double factor1loc=(l==0)?2.:1.;factor1loc*=factorX;
      double factor2loc=(l==nL-2)?2.:1.;factor2loc*=factorX;
      //cout<<"Term ("<<pos<<","<<l<<")->f1="<<factor1loc<<" f2="<<factor2loc<<endl;
      mwArray Hy;
      computeTwoBodyTerm(Hy,Jx,Jy,Jz,factor1loc*hx,factor2loc*hx,factor1loc*hy,factor2loc*hy,factor1loc*hz,factor2loc*hz);
      mwArray tmp;
      int dleft=pow(d,l);
      //putForMatlab(cout,Hy,"Hy");
      if(l==0) tmp=Hy;
      else constructOperatorProduct(tmp,identityMatrix(dleft),Hy);
      mwArray aux(term1);
      constructOperatorProduct(term1,aux,identityMatrix(d));
      term1=term1+tmp;
    }
    mwArray aux;constructOperatorProduct(aux,term1,identityMatrix(dtot));
    factorY=.5;
    factor1=pos==0?2.:1.;factor1=factor1*factorY;
    factor2=pos==Lx-2?2.:1.;factor2=factor2*factorY;
    //    cout<<" for the (local) vertical terms at "<<pos<<", factor1="<<factor1<<", (next)factor2="<<factor2<<endl;
    mwArray tmp=factor1*aux;
    result=result+factor1*aux;
    constructOperatorProduct(aux,identityMatrix(dtot),term1);
    result=result+factor2*aux;
    tmp=tmp+factor2*aux;
  }
}


void HeisenbergHamiltonian2D::computeTwoBodyTermVertical(mwArray& result,int pos) const{
  if(pos>=legNr-1){
    cout<<"ERROR: HeisenbergHamiltonian2D::computeTwoBodyTermVertical for pos(leg) "<<pos<<" when nr.legs="<<legNr<<endl;
    exit(212);
  }
  int dloc=pow(d,nL);
  // the local terms are included as 1/4 in general, and in edges 1/4 more, except if there is only one rung, in which case the factor is 1/2
  double factorX=Lx>1?.25:.5;  
  // But in this construction, also two body terms completely included in a site are summed twice andd need their own factor
  double factor1=(pos==0&&nL==1)?2.:1.;factor1*=factorX;
  double factor2=(pos==legNr-2&&nL==1)?2.:1.;factor2*=factorX;
  mwArray h12;
  computeTwoBodyTerm(h12,Jx,Jy,Jz,factor1*hx,factor2*hx,factor1*hy,factor2*hy,factor1*hz,factor2*hz);
  if(nL==1){
    result=h12;
  }
  else{
    // I need 2*nL-1 terms
    // now the "local" terms include 2-body within the effective site
    factor1=pos==0?2.:1.;factor1*=factorX;
    factor2=pos==legNr-2?2.:1.;factor2*=factorX;

    mwArray term1(Indices(dloc,dloc)),term2(Indices(dloc,dloc));

    // Just add up all the contributions to site 1 (2)
    for(int l=0;l<nL-1;l++){
      mwArray auxOp;
      auxOp=Z.subArray(Indices(l*6+4,-1)); // XX
      term1=term1+factor1*Jx*reshape(auxOp,Indices(dloc,dloc));
      term2=term2+factor2*Jx*reshape(auxOp,Indices(dloc,dloc));
      auxOp=Z.subArray(Indices(l*6+5,-1)); // YY
      term1=term1+factor1*Jy*reshape(auxOp,Indices(dloc,dloc));
      term2=term2+factor2*Jy*reshape(auxOp,Indices(dloc,dloc));
      auxOp=Z.subArray(Indices(l*6+6,-1)); // ZZ
      term1=term1+factor1*Jz*reshape(auxOp,Indices(dloc,dloc));
      term2=term2+factor2*Jz*reshape(auxOp,Indices(dloc,dloc));

      // include also local terms, but for each of them I need the truly local factors,
      // whih depend on the site and the spin position
      // for term1
      double factor1loc=(pos==0&&l==0)?2.:1.;factor1loc*=factorX;
      double factor2loc=factorX; // never called for a single leg=> cannot be the last
      // for term2
      double factor1loc2=factorX; // pos+1 never 0
      double factor2loc2=(pos==legNr-2&&l==nL-2)?2.:1.;factor2loc2*=factorX;
      auxOp=Z.subArray(Indices(l*6+1,-1)); // XI
      term1=term1+factor1*factor1loc*hx*reshape(auxOp,Indices(dloc,dloc));    
      term2=term2+factor2*factor1loc2*hx*reshape(auxOp,Indices(dloc,dloc));    
      auxOp=Z.subArray(Indices((l+1)*6+1,-1)); // IX
      term1=term1+factor1*factor2loc*hx*reshape(auxOp,Indices(dloc,dloc));
      term2=term2+factor2*factor2loc2*hx*reshape(auxOp,Indices(dloc,dloc));
      auxOp=Z.subArray(Indices(l*6+2,-1)); // YI
      term1=term1+factor1*factor1loc*hy*reshape(auxOp,Indices(dloc,dloc));    
      term2=term2+factor2*factor1loc2*hy*reshape(auxOp,Indices(dloc,dloc));    
      auxOp=Z.subArray(Indices((l+1)*6+2,-1)); // IY
      term1=term1+factor1*factor2loc*hy*reshape(auxOp,Indices(dloc,dloc));
      term2=term2+factor2*factor2loc2*hy*reshape(auxOp,Indices(dloc,dloc));
      auxOp=Z.subArray(Indices(l*6+3,-1)); // ZI
      term1=term1+factor1*factor1loc*hz*reshape(auxOp,Indices(dloc,dloc));    
      term2=term2+factor2*factor1loc2*hz*reshape(auxOp,Indices(dloc,dloc));    
      auxOp=Z.subArray(Indices((l+1)*6+3,-1)); // IZ
      term1=term1+factor1*factor2loc*hz*reshape(auxOp,Indices(dloc,dloc));
      term2=term2+factor2*factor2loc2*hz*reshape(auxOp,Indices(dloc,dloc));
    }
    constructOperatorProduct(result,term1,identityMatrix(dloc));
    mwArray aux;
    constructOperatorProduct(aux,identityMatrix(dloc),term2);
    result=result+aux;
    constructOperatorProduct(aux,identityMatrix(dloc/d),h12);
    constructOperatorProduct(h12,aux,identityMatrix(dloc/d));
    result=result+h12;
  }
}

void HeisenbergHamiltonian2D::computeTwoBodyTerm(mwArray& result,double Jx,double Jy,double Jz,
					       double hx1,double hx2,double hy1,double hy2,double hz1,double hz2) const{
  static mwArray XX,YY,ZZ,XI,IX,YI,IY,ZI,IZ;
  static bool firstCall(true);
  if(firstCall){
    mwArray Id=identityMatrix(d);
    complex_t dataX[]={ZERO_c,.5*ONE_c,.5*ONE_c,ZERO_c};
    mwArray SX(Indices(d,d),dataX);
    complex_t dataY[]={ZERO_c,.5*I_c,-.5*I_c,ZERO_c};
    mwArray SY(Indices(d,d),dataY);
    complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
    mwArray SZ(Indices(d,d),dataZ);
    constructOperatorProduct(XX,SX,SX);
    constructOperatorProduct(YY,SY,SY);
    constructOperatorProduct(ZZ,SZ,SZ);
    constructOperatorProduct(XI,SX,Id);
    constructOperatorProduct(IX,Id,SX);
    constructOperatorProduct(YI,SY,Id);
    constructOperatorProduct(IY,Id,SY);
    constructOperatorProduct(ZI,SZ,Id);
    constructOperatorProduct(IZ,Id,SZ);
    firstCall=false;
  }
  result=Jx*XX+Jy*YY+Jz*ZZ+hx1*XI+hx2*IX+hy1*YI+hy2*IY+hz1*ZI+hz2*IZ;
}


void HeisenbergHamiltonian2D::initZ(){
  int dloc=pow(d,nL);
  int nOps=nL*6+1; //nL*3+(nL-1)*3+1;
  cout<<"initZ with nL="<<nL<<" dloc="<<dloc<<" nOps="<<nOps<<endl;
  Z=mwArray(Indices(nOps,dloc,dloc));Z.fillWithZero();
  // The very basic pieces for the operators I need
  mwArray sig0=identityMatrix(d);
  complex_t dataX[]={ZERO_c,.5*ONE_c,.5*ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  complex_t dataY[]={ZERO_c,.5*I_c,-.5*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray idLoc=identityMatrix(dloc);
  for(int i1=0;i1<dloc;i1++)
    for(int i2=0;i2<dloc;i2++){
      Z.setElement(idLoc.getElement(Indices(i1,i2)),Indices(0,i1,i2));
    }
  Indices dims;for(int p=0;p<nL;p++) dims.push_back(d); // indiv indices of ops
  // the local operators
  mwArray XI,YI,ZI,XX,YY,ZZ;
  if(nL==1){
    XI=sigX;YI=sigY;ZI=sigZ;
    XX=0.*XI;YY=XX;ZZ=XX;
  }
  else{
    constructOperatorProduct(XI,sigX,identityMatrix(pow(d,nL-1)));
    constructOperatorProduct(YI,sigY,identityMatrix(pow(d,nL-1)));
    constructOperatorProduct(ZI,sigZ,identityMatrix(pow(d,nL-1)));
    {mwArray XX_,YY_,ZZ_;
      constructOperatorProduct(XX_,sigX,sigX);
      constructOperatorProduct(YY_,sigY,sigY);
      constructOperatorProduct(ZZ_,sigZ,sigZ);
      if(nL==2){XX=XX_;YY=YY_;ZZ=ZZ_;}
      else{
	constructOperatorProduct(XX,XX_,identityMatrix(pow(d,nL-2)));
	constructOperatorProduct(YY,YY_,identityMatrix(pow(d,nL-2)));
	constructOperatorProduct(ZZ,ZZ_,identityMatrix(pow(d,nL-2)));
      }
    }
    XI.reshape(Indices(dims,dims));YI.reshape(Indices(dims,dims));ZI.reshape(Indices(dims,dims));
    XX.reshape(Indices(dims,dims));YY.reshape(Indices(dims,dims));ZZ.reshape(Indices(dims,dims));
  }
  Indices order;for(int p=1;p<=2*nL;p++)order.push_back(p);
  for(int k=0;k<nL;k++){
    mwArray XI_(XI);mwArray YI_(YI);mwArray ZI_(ZI);
    if(k>0){ // first one is ok
      // reshape the composite ops: permute one of the identities with the op
      Indices newOrder(order);
      newOrder[0]=newOrder[k];
      newOrder[nL]=newOrder[k+nL];
      newOrder[k]=1;
      newOrder[k+nL]=nL+1;
      //cout<<"order: "<<order<<" neworder: "<<newOrder<<endl;
      XI_.permute(newOrder);
      YI_.permute(newOrder);
      ZI_.permute(newOrder);
    }
    XI_.reshape(Indices(dloc,dloc));YI_.reshape(Indices(dloc,dloc));ZI_.reshape(Indices(dloc,dloc));        
    for(int i1=0;i1<dloc;i1++)
      for(int i2=0;i2<dloc;i2++){
	Z.setElement(XI_.getElement(Indices(i1,i2)),Indices(k*6+1,i1,i2));
	Z.setElement(YI_.getElement(Indices(i1,i2)),Indices(k*6+2,i1,i2));
	Z.setElement(ZI_.getElement(Indices(i1,i2)),Indices(k*6+3,i1,i2));
      }
    if(k<nL-1){
      Indices newOrder2(order);
      int auxInd=newOrder2[k+1]; // exchange with second site (2)
      newOrder2[k+1]=2;newOrder2[k+nL+1]=2+nL;
      newOrder2[1]=auxInd;newOrder2[nL+1]=auxInd+nL;
      auxInd=newOrder2[k]; // exchange with 1st site (1)  
      newOrder2[k]=1;newOrder2[k+nL]=1+nL;
      newOrder2[0]=auxInd;newOrder2[nL]=auxInd+nL;
      //      cout<<"order: "<<order<<" neworder2: "<<newOrder2<<endl;
      mwArray XX_(XX);XX_.permute(newOrder2);XX_.reshape(Indices(dloc,dloc));
      mwArray YY_(YY);YY_.permute(newOrder2);YY_.reshape(Indices(dloc,dloc));
      mwArray ZZ_(ZZ);ZZ_.permute(newOrder2);ZZ_.reshape(Indices(dloc,dloc));
      for(int i1=0;i1<dloc;i1++)
	for(int i2=0;i2<dloc;i2++){
	  Z.setElement(XX_.getElement(Indices(i1,i2)),Indices(k*6+4,i1,i2));
	  Z.setElement(YY_.getElement(Indices(i1,i2)),Indices(k*6+5,i1,i2));
	  Z.setElement(ZZ_.getElement(Indices(i1,i2)),Indices(k*6+6,i1,i2));
	}
    }      
  }
  Z.reshape(Indices(nOps,dloc*dloc));
}

void HeisenbergHamiltonian2D::initHMPO(double offset){
  int D=5; // but then insert another index for leg 
  int L=Lx*legNr; // length of MPO
  int nOps=Z.getDimension(0);
  int dloc=pow(d,nL);
  int Daux=legNr; // auxiliary counter for legs
  int Ds=nL; // and counter for spin within each site
  for(int k=0;k<Lx;k++){
    // rung by rung
    for(int l=0;l<legNr;l++){ // leg by leg
      int pos=k*legNr+l; // position in the MPO
      int Dl=D,Dr=D,Dl2=Daux,Dr2=Daux,Dl3=Ds,Dr3=Ds;
      if(pos==0) Dl=1;
      if(pos==L-1) Dr=1;
      if(pos==0) Dl2=1;
      if(pos==L-1) Dr2=1;
      if(pos==0) Dl3=1;
      if(pos==L-1) Dr3=1;
      // if(k==0) Dl2=1;
      // if(k==Lx-1) Dr2=1;
      mwArray C(Indices(Dl,Dl2,Dl3,Dr,Dr2,Dr3,nOps));
      cout<<"Preparing C tensor for site ("<<k<<","<<l<<") pos="<<pos<<" (of "<<L<<") dims="<<C.getDimensions()<<endl;
      // Set elements
      if(pos!=L-1){
	C.setElement(ONE_c,Indices(0,0,0,0,0,0,0)); // nothing yet
	// First of a pair going horizontally (or vertically for all but the last leg)
	for(int p=0;p<nL;p++){
	  C.setElement(Jx*ONE_c,Indices(0,0,0,1,l,p,p*6+1));
	  C.setElement(Jy*ONE_c,Indices(0,0,0,2,l,p,p*6+2));
	  C.setElement(Jz*ONE_c,Indices(0,0,0,3,l,p,p*6+3));
	  //cout<<"First of pair (leg"<<l<<",site"<<p<<") at element(x) "<<Indices(0,0,0,1,l,p,p*6+1)<<endl;
	}
	// Except first and last, I have to let two-body terms through
	if(pos!=0){
	  for(int kl=0;kl<legNr;kl++){
	    if(kl!=l){
	      	for(int p=0;p<nL;p++){
		  C.setElement(ONE_c,Indices(1,kl,p,1,kl,p,0));
		  C.setElement(ONE_c,Indices(2,kl,p,2,kl,p,0));
		  C.setElement(ONE_c,Indices(3,kl,p,3,kl,p,0));
		  //cout<<"Passing pair (leg"<<kl<<",site"<<p<<") (x) "<<Indices(1,kl,p,1,kl,p,0)<<endl;
		}
	    }
	  }
	}
      }
      if(pos!=0){ 
	// After everything	
	C.setElement(ONE_c,Indices(Dl-1,0,0,Dr-1,0,0,0));
	// Second of the pair
	if(k>0){ // in the first rung I cannot have the second of sth from left
	  for(int p=0;p<nL;p++){
	    C.setElement(ONE_c,Indices(1,l,p,Dr-1,0,0,p*6+1));
	    C.setElement(ONE_c,Indices(2,l,p,Dr-1,0,0,p*6+2));
	    C.setElement(ONE_c,Indices(3,l,p,Dr-1,0,0,p*6+3));
	    //cout<<"Second of pair (leg"<<l<<",site"<<p<<") at element(x) "<<Indices(1,l,p,Dr-1,0,0,p*6+1)<<endl;
	  }
	}
	// In all rungs, and all legs after the first, I can have the second half of the pair from below
	if(l>0){
	  C.setElement(ONE_c,Indices(1,l-1,nL-1,Dr-1,0,0,1));
	  C.setElement(ONE_c,Indices(2,l-1,nL-1,Dr-1,0,0,2));
	  C.setElement(ONE_c,Indices(3,l-1,nL-1,Dr-1,0,0,3));
	  //cout<<"Second of pair (from leg"<<l-1<<") at element(x) "<<Indices(1,l-1,nL-1,Dr-1,0,0,1)<<endl;
	}
      }
      // The single-site term(s)
      for(int p=0;p<nL;p++){
	C.setElement(hx*ONE_c,Indices(0,0,0,Dr-1,0,0,p*6+1));
	C.setElement(hy*ONE_c,Indices(0,0,0,Dr-1,0,0,p*6+2));
	C.setElement(hz*ONE_c,Indices(0,0,0,Dr-1,0,0,p*6+3));
	// the two-site terms within an effective site
	if(p<nL-1){
	  C.setElement(Jx*ONE_c,Indices(0,0,0,Dr-1,0,0,p*6+4));
	  C.setElement(Jy*ONE_c,Indices(0,0,0,Dr-1,0,0,p*6+5));
	  C.setElement(Jz*ONE_c,Indices(0,0,0,Dr-1,0,0,p*6+6));
	  //cout<<"Local two body term(xx) at element "<<Indices(0,0,0,Dr-1,0,0,p*6+4)<<endl;
	}
      }
      // The offset term, if present
      if(offset!=0)
	C.setElement((offset/L)*ONE_c,Indices(0,0,0,Dr-1,0,0,0));
      // Reshape, multiply operators and set in MPO
      C.reshape(Indices(Dl*Dl2*Dl3*Dr*Dr2*Dr3,nOps));
      C.multiplyRight(Z);
      C.reshape(Indices(Dl*Dl2*Dl3,Dr*Dr2*Dr3,dloc,dloc));
      C.permute(Indices(3,1,4,2));
      hamil.setOp(pos,new Operator(C),true);
    }
  }
}

void HeisenbergHamiltonian2D::getExponentialMPOx(MPO& expH,complex_t delta,bool even) const {
  // Notice that this implementation has many copies of the same operators and could be optimized.
  mwArray Ol,Or;
  int L=Lx*legNr; // length of the MPOs
  int dloc=pow(d,nL);
  // I will need local operator for identities
  mwArray opIdloc=identityMatrix(dloc);opIdloc.reshape(Indices(dloc,1,dloc,1));
  // I will make legNr MPOs and then join them in expH
  MPO ** mpos=new MPO*[legNr];
  for(int k=0;k<legNr;k++){
    mpos[k]=new MPO(L);
    // fill them initially with identities, and then just replace the operators that are different
    for(int p=0;p<L;p++) mpos[k]->setOp(p,new Operator(opIdloc),true);
  }
  // rung by rung I set in each MPO the terms that connect each effective site to the neighbor
  int firstRung=even?0:1;
  int lastRung=(Lx-1)-(Lx-firstRung)%2; // last occupied by the loop
  cout<<"Constructing ExponentialMPOx even="<<even<<" firstRung="<<firstRung<<" lastRung="<<lastRung<<endl;

  // if even, the first rung occupied on leg 0 is 0, o.w. it is 1
  for(int r=0;r<Lx-1;r++){
    int firstLeg=(even==(r%2==0))?0:1;
    getTwoBodyTermExponentialx(Ol,Or,delta,r); //r*legNr);
    int Dop=Ol.getDimension(3);
    // also the intermediate identity, that passes the bond    
    mwArray opId(opIdloc);
    opId.reshape(Indices(dloc*dloc,1));opId.multiplyRight(reshape(identityMatrix(Dop),Indices(1,Dop*Dop)));
    opId.reshape(Indices(dloc,dloc,Dop,Dop));opId.permute(Indices(1,3,2,4));
    // leg by leg, set the ops
    for(int p=firstLeg;p<legNr;p=p+2){
      mpos[p]->setOp(p+r*legNr,new Operator(Ol),true);
      for(int k=1;k<legNr;k++) mpos[p]->setOp(p+r*legNr+k,new Operator(opId),true);
      mpos[p]->setOp(p+(r+1)*legNr,new Operator(Or),true);
      //cout<<"Set Opl("<<p+r*legNr<<") Or("<<p+(r+1)*legNr<<")"<<endl;      
    }
    // putForMatlab(cout,Ol,"Ol");
    // putForMatlab(cout,Or,"Or");
  }
  // And finally join the MPOS
  const MPO* mpos_[legNr];for(int k=0;k<legNr;k++) mpos_[k]=mpos[k];
  MPO::join(legNr,mpos_,expH);
  //cout<<"Created exponentialMPOx(even="<<even<<"): "<<expH<<endl;
  for(int k=0;k<legNr;k++){
    mpos[k]->clear();
    delete mpos[k];
  }
  delete []mpos;
}

// void HeisenbergHamiltonian2D::getExponentialMPOx(MPO& expH,complex_t delta,bool even) const {
//   // Notice that this implementation has many copies of the same operators and could be optimized.
//   mwArray Ol,Or;
//   int L=Lx*legNr; // length of the MPOs
//   int dloc=pow(d,nL);
//   // I will need local operator for identities
//   mwArray opIdloc=identityMatrix(dloc);opIdloc.reshape(Indices(dloc,1,dloc,1));
//   // I will make legNr MPOs and then join them in expH
//   MPO ** mpos=new MPO*[legNr];
//   for(int k=0;k<legNr;k++){
//     mpos[k]=new MPO(L);
//     // fill them initially with identities, and then just replace the operators that are different
//     for(int p=0;p<L;p++) mpos[k]->setOp(p,new Operator(opIdloc),true);
//   }
//   // rung by rung I set in each MPO the terms that connect each effective site to the neighbor
//   int firstRung=even?0:1;
//   int lastRung=(Lx-1)-(Lx-firstRung)%2; // last occupied by the loop
//   cout<<"Constructing ExponentialMPOx even="<<even<<" firstRung="<<firstRung<<" lastRung="<<lastRung<<endl;
//   for(int r=firstRung;r<lastRung;r+=2){
//     // Because everything TI, it is the same operator for all horizontal
//     // terms (first and last sites/rungs have a difference in local
//     // terms) I have legNr copies, displaced by one site, so only need to get the ops for the first leg
//     getTwoBodyTermExponentialx(Ol,Or,delta,r); //r*legNr);
//     int Dop=Ol.getDimension(3);
//     // also the intermediate identity, that passes the bond    
//     mwArray opId(opIdloc);
//     opId.reshape(Indices(dloc*dloc,1));opId.multiplyRight(reshape(identityMatrix(Dop),Indices(1,Dop*Dop)));
//     opId.reshape(Indices(dloc,dloc,Dop,Dop));opId.permute(Indices(1,3,2,4));
//     // leg by leg, set the ops
//     for(int p=0;p<legNr;p++){
//       mpos[p]->setOp(p+r*legNr,new Operator(Ol),true);
//       for(int k=1;k<legNr;k++) mpos[p]->setOp(p+r*legNr+k,new Operator(opId),true);
//       mpos[p]->setOp(p+(r+1)*legNr,new Operator(Or),true);
//       //cout<<"Set Opl("<<p+r*legNr<<") Or("<<p+(r+1)*legNr<<")"<<endl;      
//     }
//     // putForMatlab(cout,Ol,"Ol");
//     // putForMatlab(cout,Or,"Or");
//   }
//   // And finally join the MPOS
//   const MPO* mpos_[legNr];for(int k=0;k<legNr;k++) mpos_[k]=mpos[k];
//   MPO::join(legNr,mpos_,expH);
//   //cout<<"Created exponentialMPOx(even="<<even<<"): "<<expH<<endl;
//   for(int k=0;k<legNr;k++){
//     mpos[k]->clear();
//     delete mpos[k];
//   }
//   delete []mpos;
// }

void HeisenbergHamiltonian2D::getExponentialMPOy(MPO& expH,complex_t delta,bool even) const {
  // Notice that this implementation has many copies of the same operators and could be optimized.
  mwArray Ol,Or;
  int L=Lx*legNr; // length of the MPOs
  int dloc=pow(d,nL);
  // I will need local operator for identities
  mwArray opIdloc=identityMatrix(dloc);opIdloc.reshape(Indices(dloc,1,dloc,1));
  expH.initLength(L);
  // fill initially with identities, and then just replace the operators taht are different
  for(int p=0;p<L;p++) expH.setOp(p,new Operator(opIdloc),true);

  int kfirst=even?0:1;
  int klast=(legNr-1)-(legNr-kfirst)%2; // last leg occupied by the loop

  // leg by leg I set in each MPO the terms that connect each effective site to the next neighbour
  for(int k=kfirst;k<klast;k+=2){
    // for each rung, I need almost the same pair of operators for the same legs.
    //The first and last rungs might be different, but to avoid that, the local terms are included in the horizontal terms    
    getTwoBodyTermExponentialy(Ol,Or,delta,k);
    for(int r=0;r<Lx;r++){
      // set the ops
      expH.setOp(k+r*legNr,new Operator(Ol),true);
      expH.setOp(k+1+r*legNr,new Operator(Or),true);
    }    
  }
}


void HeisenbergHamiltonian2D::getDoubleExponentialMPOx(MPO& expH,complex_t delta,bool even) const {
  // Notice that this implementation has many copies of the same operators and could be optimized.
  mwArray Ol,Or,OlA,OrA;
  int L=Lx*legNr; // length of the MPOs
  int dloc=pow(d,nL);
  // I will need local operator for identities
  mwArray opIdloc=identityMatrix(dloc);opIdloc.reshape(Indices(dloc,1,dloc,1));
  // I will make legNr MPOs and then join them in expH
  MPO** mpos=new MPO*[legNr];
  for(int k=0;k<legNr;k++){
    mpos[k]=new MPO(L);
    // fill them initially with identities, and then just replace the operators taht are different
    for(int p=0;p<L;p++) mpos[k]->setOp(p,new DoubleOperator(opIdloc,opIdloc),true);
  }
  // rung by rung I set in each MPO the terms that connect each effective site to the neighbor
  int firstRung=even?0:1;
  int lastRung=(Lx-1)-(Lx-firstRung)%2; // last occupied by the loop
  cout<<"Constructing DoubleExponentialMPOx even="<<even<<" firstRung="<<firstRung<<" lastRung="<<lastRung<<endl;
 
  for(int r=firstRung;r<lastRung;r+=2){
    mwArray H12;
    computeTwoBodyTermHorizontal(H12,r);
    getTwoBodyTermExponential(Ol,Or,H12,delta);
    getTwoBodyTermExponential(OlA,OrA,H12,conjugate(delta));

    int Dop=Ol.getDimension(3);
    // also the intermediate identity, that passes the bond    
    mwArray opId(opIdloc);
    opId.reshape(Indices(dloc*dloc,1));opId.multiplyRight(reshape(identityMatrix(Dop),Indices(1,Dop*Dop)));
    opId.reshape(Indices(dloc,dloc,Dop,Dop));opId.permute(Indices(1,3,2,4));
    // leg by leg, set the ops
    for(int p=0;p<legNr;p++){
      mpos[p]->setOp(p+r*legNr,new DoubleOperator(Ol,permute(OlA,Indices(3,2,1,4))),true);
      for(int k=1;k<legNr;k++) mpos[p]->setOp(p+r*legNr+k,new DoubleOperator(opId,permute(opId,Indices(3,2,1,4))),true);
      mpos[p]->setOp(p+(r+1)*legNr,new DoubleOperator(Or,permute(OrA,Indices(3,2,1,4))),true);
    }
  }
    // And finally join the MPOS
  const MPO* mpos_[legNr];for(int k=0;k<legNr;k++) mpos_[k]=mpos[k];
  MPO::join(legNr,mpos_,expH);
  for(int k=0;k<legNr;k++){
    mpos[k]->clear();
    delete mpos[k];
  }
  delete []mpos;
}

void HeisenbergHamiltonian2D::getDoubleExponentialMPOy(MPO& expH,complex_t delta,bool even) const {
  // Notice that this implementation has many copies of the same operators and could be optimized.
  mwArray Ol,Or,OlA,OrA;
  int L=Lx*legNr; // length of the MPO
  int dloc=pow(d,nL);
  // I will need local operator for identities
  mwArray opIdloc=identityMatrix(dloc);opIdloc.reshape(Indices(dloc,1,dloc,1));
  expH.initLength(L);
  // fill initially with identities, and then just replace the operators taht are different
  for(int p=0;p<L;p++) expH.setOp(p,new DoubleOperator(opIdloc,opIdloc),true);

  int kfirst=even?0:1;
  int klast=(legNr-1)-(legNr-kfirst)%2; // last leg occupied by the loop

  // leg by leg I set in each MPO the terms that connect each effective site to the next neighbour
  for(int k=kfirst;k<klast;k+=2){
    // for each rung, I need almost the same pair of operators 
    //The first and last rungs might be different, but to avoid that, the local terms are included in the horizontal terms
    mwArray H12;
    computeTwoBodyTermVertical(H12,k);
    getTwoBodyTermExponential(Ol,Or,H12,delta);
    getTwoBodyTermExponential(OlA,OrA,H12,conjugate(delta));
    for(int r=0;r<Lx;r++){ // each rung
      // set the ops
      expH.setOp(k+r*legNr,new DoubleOperator(Ol,permute(OlA,Indices(3,2,1,4))),true);
      expH.setOp(k+1+r*legNr,new DoubleOperator(Or,permute(OrA,Indices(3,2,1,4))),true);
    }    
  }
}


void HeisenbergHamiltonian2D::getExtendedExponentialMPOx(MPO& expH,complex_t delta,bool even) const {
  // Notice that this implementation has many copies of the same operators and could be optimized.
  mwArray Ol,Or,OlA,OrA;
  int L=Lx*legNr; // length of the MPOs
  int dloc=pow(d,nL);
  // I will need local operator for identities
  mwArray opIdloc=identityMatrix(dloc);opIdloc.reshape(Indices(dloc,1,dloc,1));
  // I will make legNr MPOs and then join them in expH
  MPO** mpos=new MPO*[legNr];
  for(int k=0;k<legNr;k++){
    mpos[k]=new MPO(L);
    // fill them initially with identities, and then just replace the operators taht are different
    for(int p=0;p<L;p++) mpos[k]->setOp(p,new DoubleOperator(opIdloc,opIdloc),true);
  }
  // rung by rung I set in each MPO the terms that connect each effective site to the neighbor
  int firstRung=even?0:1;
  int lastRung=(Lx-1)-(Lx-firstRung)%2; // last occupied by the loop
  
  for(int r=firstRung;r<lastRung;r+=2){
    mwArray H12;
    computeTwoBodyTermHorizontal(H12,r);
    getTwoBodyTermExponential(Ol,Or,H12,delta);
    getTwoBodyTermExponential(OlA,OrA,H12,conjugate(delta));

    int Dop=Ol.getDimension(3);
    // also the intermediate identity, that passes the bond    
    mwArray opId(opIdloc);
    opId.reshape(Indices(dloc*dloc,1));opId.multiplyRight(reshape(identityMatrix(Dop),Indices(1,Dop*Dop)));
    opId.reshape(Indices(dloc,dloc,Dop,Dop));opId.permute(Indices(1,3,2,4));
    // leg by leg, set the ops
    for(int p=0;p<legNr;p++){
      mpos[p]->setOp(p+r*legNr,new DoubleOperator(Ol,permute(OlA,Indices(3,2,1,4))),true);
      for(int k=1;k<legNr;k++) mpos[p]->setOp(p+r*legNr+k,new DoubleOperator(opId,permute(opId,Indices(3,2,1,4))),true);
      mpos[p]->setOp(p+(r+1)*legNr,new DoubleOperator(Or,permute(OrA,Indices(3,2,1,4))),true);
    }
  }
    // And finally join the MPOS
  const MPO* mpos_[legNr];for(int k=0;k<legNr;k++) mpos_[k]=mpos[k];
  MPO::join(legNr,mpos_,expH);
  for(int k=0;k<legNr;k++){
    mpos[k]->clear();
  }
  delete []mpos;
}

void HeisenbergHamiltonian2D::getExtendedExponentialMPOy(MPO& expH,complex_t delta,bool even) const {
  // Notice that this implementation has many copies of the same operators and could be optimized.
  mwArray Ol,Or,OlA,OrA;
  int L=Lx*legNr; // length of the MPO
  int dloc=pow(d,nL);
  // I will need local operator for identities
  mwArray opIdloc=identityMatrix(dloc);opIdloc.reshape(Indices(dloc,1,dloc,1));
  expH.initLength(L);
  // fill initially with identities, and then just replace the operators taht are different
  for(int p=0;p<L;p++) expH.setOp(p,new DoubleOperator(opIdloc,opIdloc),true);

  int kfirst=even?0:1;
  int klast=(legNr-1)-(legNr-kfirst)%2; // last leg occupied by the loop

  // leg by leg I set in each MPO the terms that connect each effective site to the next neighbour
  for(int k=kfirst;k<klast;k+=2){
    // for each rung, I need almost the same pair of operators 
    //The first and last rungs might be different, but to avoid that, the local terms are included in the horizontal terms
    mwArray H12;
    computeTwoBodyTermVertical(H12,k);
    getTwoBodyTermExponential(Ol,Or,H12,delta);
    for(int r=0;r<Lx;r++){ // each rung
      // set the ops
      expH.setOp(k+r*legNr,new DoubleOperator(Ol,opIdloc),true);
      expH.setOp(k+1+r*legNr,new DoubleOperator(Or,opIdloc),true);
    }    
  }
}


