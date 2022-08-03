/**
   \file Site.cpp
   Implementation of the elementary class containing the tensor for a
   single site in a MPS 
   
   \author Mari-Carmen Banuls
   \date 16/03/2011
*/

#include "Site.h"
#include <cmath>

using namespace shrt;


Site::Site(int posit,int dphys,int Dleft,int Dright):
  pos(posit),A(Indices(dphys,Dleft,Dright)){
  //dims[0]=dphys;dims[1]=Dleft;dims[2]=Dright;
  //A.resize(dims);
  A.setElement(1.,0.,Indices(0,0,0));
}

Site::Site(int posit,int* dimensions):
  pos(posit),A(3,dimensions) {
  //  dims[0]=(int)dimensions[0];
  //dims[1]=(int)dimensions[1];
  //dims[2]=(int)dimensions[2];
  A.setElement(1.,0.,Indices(0,0,0));
  //cout<<"Created Site("<<pos<<" with dims "<<A.getDimensions()<<endl;
}


Site::Site(const Site& s,int _pos):
  pos(_pos==0?s.pos:_pos),A(s.A) {}
  //  pos(_pos==0?s.pos:_pos),dims(s.dims),A(s.A) {}

Site::Site(int pos_,const mwArray& A_):
  pos(pos_),A(A_){
  Indices dims=A.getDimensions();
  if(A.getRank()<3){
    dims.push_back(1);
    if(A.getRank()<2) dims.push_back(1);
    A.reshape(dims);
  }
  if(A.getRank()==0||A.getRank()>3){
    cout<<"Error: cannot create a Site for pos "<<pos_
	<<" from array of dims "<<A_.getDimensions()<<endl;
    exit(212);
  }
}

Site::~Site(){}


void Site::setSite(const mwArray& newA){
  // Right now, this makes a copy of this data, but maybe it could
  // make it more efficiently! Clone()?
  Indices dims=A.getDimensions();
  if((newA.getRank()<=0||newA.getRank()>3)||
     (newA.getNrComponents()!=A.getNrComponents())){
    cout<<"Error Site::setSite trying to set matrix pos("<<pos<<
      ")of different or invalid dimensions!"
      " old:"<<A.getDimensions()<<", new:"<<newA.getDimensions()<<endl;
    exit(2);
  }
  A=reshape(newA,dims);
}

void Site::setRotatedSite(const mwArray& newA,const Indices& newdims,
			  bool conjugate){
  A=newA;
  A.permute(newdims);
  if(conjugate)
    A.conjugate();
  Indices dims=A.getDimensions();
  if(A.getRank()<3){
    cout<<"WARNING: setRotatedsite with rank "<<A.getRank()<<endl;
    dims.push_back(1);
    if(A.getRank()<2) dims.push_back(1);
  }
}


ostream& operator<<(ostream& os,const Site& site){
  site.put(os);
  return os;
  //  return os<<"Site("<<site.pos<<"), tensor:"<<site.A<<endl;
}

void Site::put(ostream& os) const{
  os<<"Site("<<pos<<"), tensor:"<<A.getDimensions()<<endl;
}


void Site::setToZero(){
  A.fillWithZero();
}

void Site::setState(SingleState state){
  setToZero();
  Indices pos1(0,0,0);
  Indices pos2(1,0,0);
  switch(state){
  case zero:
    A.setElement(1.,0.,pos1);
    break;
  case one:
    A.setElement(1.,0.,pos2);
    break;
  case xplus:
    A.setElement(1./sqrt(2.),0.,pos1);
    A.setElement(1./sqrt(2.),0.,pos2);
    break;
  case xminus:
    A.setElement(1./sqrt(2.),0.,pos1);
    A.setElement(-1./sqrt(2.),0.,pos2);
    break;
  case yplus:
    A.setElement(1./sqrt(2.),0.,pos1);
    A.setElement(0.,1./sqrt(2.),pos2);
    break;
  case yminus:
    A.setElement(1./sqrt(2.),0.,pos1);
    A.setElement(0.,-1./sqrt(2.),pos2);
    break;
  case special:
    {int dim0=A.getDimension(0);
    if(dim0<4){
      cout<<"Cannot create special state if dim="<<dim0<<endl;
      exit(1);
    }
    pos1[0]=1;
    pos2[0]=dim0-2; // the last one (|01>+|10>)
    A.setElement(1./sqrt(2.),0.,pos1);
    A.setElement(1./sqrt(2.),0.,pos2);}
    break;
  case maxent:
    {int dim0=A.getDimension(0);
      // Assume it is a perfect square
      int dbas=(int)sqrt(dim0*1.);
      //    if(dim0!=4){
      if(dbas*dbas!=dim0){
	cout<<"Cannot create maximally entangled state if dim="<<dim0<<endl;
	exit(1);
      }
      double val=1/sqrt(dbas*1.);
      int Dleft=A.getDimension(1);
      int Dright=A.getDimension(2);
      A.reshape(Indices(dbas,dbas,Dleft,Dright));
      for(int cntr=0;cntr<dbas;cntr++){
	A.setElement(val,0.,Indices(cntr,cntr,0,0));
      }
      A.reshape(Indices(dim0,Dleft,Dright));}
    //pos1[0]=0;
    //pos2[0]=dim0-1; // the last one (|01>+|10>)
    //A.setElement(1./sqrt(2.),0.,pos1);
    //A.setElement(1./sqrt(2.),0.,pos2);}
    break;
  case randmps:
    A.fillRandom();
    break;
  default:
    cout<<"Error: state "<<state<<"not known: state set to zero"<<endl;
    A.setElement(1.,0.,pos1);
  }
}

Site& Site::operator=(const Site& s){
  if(this!=&s){
    pos=s.pos;
    A=s.A;
    //A.getDimensions(dims);
  }
  return *this;
}

/****** OJO!!!!!!!!! *****/

void Site::gaugeR(mwArray& rightterm,const mwArray& leftterm,int cutD){
  reduceDimR(leftterm);
  int d=A.getDimension(0);int Dl=A.getDimension(1);int Dr=A.getDimension(2);
  Indices packedDims(d*Dl,Dr);
  mwArray aux=reshape(A,packedDims);
  if(cutD==0){
    //cout<<"In gaugeR, before QR, aux="<<aux<<endl;
    wrapper::qr(aux,A,rightterm);
    //cout<<"In gaugeR, right after QR, A="<<A<<endl;
    //cout<<"\t rightterm="<<rightterm<<endl;
  }
  else{
    mwArray S,Vdagger;
    //cout<<"In gaugeR, before SVD, aux="<<aux<<endl;
    double tol=0.;
    wrapper::svd(aux,tol,cutD,A,S,Vdagger);
    rightterm=S*Vdagger;
    // cout<<"In gaugeR, right after SVD, A="<<A<<endl;
    // cout<<"\t S="<<S<<endl;
    // cout<<"\t Vdagger="<<Vdagger<<endl;
  }
  Dr=A.getDimension(1); // new Dr
  A.reshape(Indices(d,Dl,Dr));
  if(Dr==1&&rightterm.getDimension(1)==1){//last site => do not include phase in the tensor
    complex_t aux=rightterm.getElement(0); // single element
    complex_t phase=(1./abs(aux))*aux;
    A.multiplyRight(mwArray(phase));
    rightterm.setElement((complex_t){abs(aux),0.},0);
  }
}

void Site::gaugeL(mwArray& leftterm,const mwArray& rightterm,int cutD){
  reduceDimL(rightterm);
  //Since A is dxDlxDr, I have to permute indices
  Indices neworder(2,1,3);
  A.permute(neworder);
  Indices dims=A.getDimensions(); // now Dl, d, Dr
  // And reshape to Dlx(d Dr)
  Indices packedDims(dims[0],dims[1]*dims[2]);
  if(cutD==0){
    wrapper::lq(reshape(A,packedDims),leftterm,A);
  }
  else{
    // Plain SVD works, but gets slower if dims are very different; move to eig?
    mwArray U,S;
    double tol=0.;
    wrapper::svd(reshape(A,packedDims),tol,cutD,U,S,A);
    leftterm=U*S;
  }
  int Dl=A.getDimension(0); // new Dl'
  dims[0]=Dl;
  A.reshape(dims); // Dl' d Dr
  A.permute(neworder);
  if(Dl==1&&leftterm.getDimension(0)==1){//fist site => do not include phase in the tensor
    complex_t aux=leftterm.getElement(0); // single element
    complex_t phase=(1./abs(aux))*aux;
    A.multiplyRight(mwArray(phase));
    leftterm.setElement((complex_t){abs(aux),0.},0);
  }
}

void Site::reduceDimR(const mwArray& leftterm){
  //if(leftterm.getRank()==1&&leftterm.getDimension(0)==1){
  if(leftterm.isScalar()){
    // multiply by scalar only
    A=leftterm.getComponents()*A;
  //A=leftterm.getElement(0)*A;
  }
  else{
    int newD=leftterm.getDimension(0);
    Indices neworder(2,1,3);
    Indices dims=A.getDimensions(); // d Dl Dr
    Indices packedDims(dims[1],dims[0]*dims[2]);
    A.permute(neworder);
    A.reshape(packedDims);
    A=leftterm*A; //newD dxDr
    //A=leftterm*reshape(permute(A,neworder),packedDims); //newD dxDr
    if(newD!=dims[1]) dims[1]=newD;
    Indices newdims(dims[1],dims[0],dims[2]);
    A.reshape(newdims);
    A.permute(neworder);
  }
}

void Site::reduceDimL(const mwArray& rightterm){
  //cout<<"Site::reduceDimL("<<rightterm<<"), dims "<<rightterm.getDimensions()<<endl;
  //if(rightterm.getRank()==1&&rightterm.getDimension(0)==1){
  if(rightterm.isScalar()){
    // multiply by scalar only
    //cout<<"Site::reduceDimL multiplying by scalar"<<endl;
    A=rightterm.getComponents()*A;
  }
  else{
    //  cout<<"Site::reduceDimL multiplying by a term dims "<<rightterm.getDimensions()<<endl;
    int newD=rightterm.getDimension(1); // new right dim of this tensor
    Indices dims=A.getDimensions();
    Indices packedDims(dims[0]*dims[1],dims[2]);
    //packedDims[0]=dims[0]*dims[1];packedDims[1]=dims[2];
    A.reshape(packedDims);
    A=A*rightterm;
    //A=reshape(A,packedDims)*rightterm;
    if(newD!=dims[2]) dims[2]=newD;
    A.reshape(dims);
  }
}

void Site::gaugeLOp(mwArray& leftterm,const mwArray& oper,
		    const mwArray& rightterm,bool dagger,int cutD){
  // contract the three tensors together
  //cout<<"Site::gaugeLOp"<<endl;
  reduceDimLOp(oper,rightterm,dagger);
  // and apply standard gauge condition to the left
  gaugeL(leftterm,mwArray(ONE_c),cutD);
}

void Site::gaugeROp(mwArray& rightterm,const mwArray& oper,
		    const mwArray& leftterm,bool dagger,int cutD){
  // contract the three tensors together
  //cout<<"Site::gaugeLOp"<<endl;
  reduceDimROp(oper,leftterm,dagger); // MISSING!! TODO!!
  // and apply standard gauge condition to the right
  gaugeR(rightterm,mwArray(ONE_c),cutD);
}

void Site::reduceDimLOp(const mwArray& oper_,const mwArray& rightterm,
			bool dagger){
  //  if(rightterm.getRank()==1&&rightterm.getDimension(0)==1){
  //A=rightterm.getComponents()*A;
  //}
  //cout<<"basic Site::reduceDimLOp on A("<<pos<<"):"<<*this<<" oper "<<oper_<<endl;
  // putForMatlab(cout,A,"A0");
  // putForMatlab(cout,oper_,"O");
  // putForMatlab(cout,rightterm,"tmpR");
  mwArray oper=oper_;
  if(dagger){
    oper.permute(Indices(3,2,1,4),true); // permute and conjugate at once
    //oper=permute(conjugate(oper),Indices(3,2,1,4));
  }
  int du=oper.getDimension(0);int Dlo=oper.getDimension(1);
  int dd=oper.getDimension(2);int Dro=oper.getDimension(3);
  int d=getd();int Dl=getDl();int Dr=getDr();
  int beta=rightterm.getDimension(rightterm.getRank()-1); // right bond
  mwArray tmp(rightterm);
  #ifdef CHECKDIMS
  if((dd!=d)||(tmp.getNrComponents()!=Dr*Dro*beta)){
    cout<<"Error in Site::reduceDimLOp dimensions: Site="<<getDimensions()
	<<", oper="<<oper.getDimensions()<<", tmpR="<<tmp.getDimensions()
	<<" dd="<<dd<<" d="<<d<<" Dr="<<Dr<<" Dro="<<Dro<<" beta="<<beta
	<<endl;
    exit(212);
  }
  #endif
  if(tmp.getRank()<3){
    // TODO: check dimensions
    // if dimensions are very wrong, this should fail here
    tmp.reshape(Indices(Dro,Dr,beta));
  }
  // A* rightterm -> Dl*d, Dro*beta
  A.permute(Indices(2,1,3));
  A.reshape(Indices(Dl*d,Dr));
  tmp.permute(Indices(2,1,3));
  tmp.reshape(Indices(Dr,Dro*beta));
  A=A*tmp;
  // oper * tmp -> du*Dlo, Dl*beta
  tmp=reshape(oper,Indices(du*Dlo,dd*Dro));
  A.reshape(Indices(Dl,d*Dro,beta));
  A.permute(Indices(2,1,3));
  A.reshape(Indices(d*Dro,Dl*beta));
  A=tmp*A;
  A.reshape(Indices(du,Dlo*Dl,beta));
  //  cout<<"After multiplying, A("<<pos<<")="<<endl;putForMatlab(cout,A,"A1");
  //setSite(reshape(tmp,Indices(du,Dlo*Dl,beta)));
}

void Site::reduceDimROp(const mwArray& oper_,const mwArray& leftterm,
			bool dagger){
  //  if(rightterm.getRank()==1&&rightterm.getDimension(0)==1){
  //A=rightterm.getComponents()*A;
  //}
  //  cout<<"basic Site::reduceDimLOp on A("<<pos<<"):"<<*this<<endl;
  // putForMatlab(cout,A,"A0");
  // putForMatlab(cout,oper_,"O");
  // putForMatlab(cout,rightterm,"tmpR");
  mwArray oper=oper_;
  if(dagger){
    oper.permute(Indices(3,2,1,4),true); // permute and conjugate at once
    //oper=permute(conjugate(oper),Indices(3,2,1,4));
  }
  int du=oper.getDimension(0);int Dlo=oper.getDimension(1);
  int dd=oper.getDimension(2);int Dro=oper.getDimension(3);
  int d=getd();int Dl=getDl();int Dr=getDr();
  int beta=leftterm.getDimension(0); // left bond
  mwArray tmp(leftterm);
#ifdef CHECKDIMS
  if((dd!=d)||(tmp.getNrComponents()!=Dl*Dlo*beta)){
    cout<<"Error in Site::reduceDimROp dimensions: Site="<<getDimensions()
	<<", oper="<<oper.getDimensions()<<", tmpR="<<tmp.getDimensions()<<endl;
    exit(212);
  }
#endif
  if(tmp.getRank()<3){
    // TODO: check dimensions
    // if dimensions are very wrong, this should fail here
    tmp.reshape(Indices(beta,Dlo,Dl));
  }
  //****** HERE!!!!!
  //****** HERE!!!!!
  // leftterm*A -> beta*Dlo, d*Dr
  A.permute(Indices(2,1,3));// D d Dr
  A.reshape(Indices(Dl,d*Dr));
  tmp.reshape(Indices(beta*Dlo,Dl));
  A.multiplyLeft(tmp); //beta*Dlo, d*Dr
  A.reshape(Indices(beta,Dlo*d,Dr));
  A.permute(Indices(2,1,3));
  A.reshape(Indices(Dlo*d,beta*Dr));
  A.multiplyLeft(reshape(permute(reshape(oper,Indices(du,Dlo*dd,Dro)),
				 Indices(1,3,2)),
			 Indices(du*Dro,Dlo*dd))); // du*Dro,beta*Dr
  A.reshape(Indices(du,Dro,beta,Dr));
  A.permute(Indices(1,3,2,4));
  A.reshape(Indices(du,beta,Dro*Dr));
  //  cout<<"After multiplying, A("<<pos<<")="<<endl;putForMatlab(cout,A,"A1");
  //setSite(reshape(tmp,Indices(du,Dlo*Dl,beta)));
}

void Site::changeSign(){
  A.changeSign();
}

void Site::conjugateA(){
  A.conjugate();
}

void Site::gaugeR(Site* ket,Site* bra,mwArray& righttermK,
		  mwArray& righttermB,const mwArray& lefttermK,
		  const mwArray& lefttermB){
  Indices transpose(2,1); // vector (2,1) for transp.
  // Multiply terms from previous site
  ket->reduceDimR(lefttermK);
  bra->reduceDimR(lefttermB);
  // We want to impose Sum B+ A=1 => contract
  Indices dimA=ket->A.getDimensions();
  Indices dimB=bra->A.getDimensions();
  Indices dim2A((ket->getd())*(ket->getDl()),ket->getDr());
  //dim2A[0]=(ket->getd())*(ket->getDl());dim2A[1]=(ket->getDr());
  Indices dim2B((bra->getd())*(bra->getDl()),bra->getDr());
  //dim2B[0]=(bra->getd())*(bra->getDl());dim2B[1]=(bra->getDr());
  mwArray A=reshape(ket->A,dim2A);
  mwArray B=reshape(bra->A,dim2B);
  mwArray M;
  M=Hconjugate(B)*A;
  mwArray U,S,Vdagger;
  double tol=1E-10;int Xi=0;
  wrapper::svd(M,0.,Xi,U,S,Vdagger);
  mwArray auxS=sqrt(S);
  mwArray invS_=invertDiag(auxS,Xi,tol); // pseudoinverse
  double trSre,trSim;
  auxS.trace(trSre,trSim);
  auxS=(1./trSre)*auxS;
  
  //cout<<"Inside gaugeR->auxS="<<auxS<<"\n\tinvS_="<<invS_<<endl;

  if(A.getDimension(1)>Xi&&B.getDimension(1)<=Xi){
    cout<<"BIG ERROR!"<<endl;
    exit(212);
    //ket->A=A*Hconjugate(Vdagger);
    //bra->A=B*U*conjugate(invS_)*conjugate(invS_);
    //righttermK=Vdagger;
    //righttermB=conjugate(auxS)*conjugate(auxS)*Hconjugate(U);
  }
  else if(A.getDimension(1)<=Xi&&B.getDimension(1)>Xi){
    cout<<"BIG ERROR!"<<endl;
    exit(212);
    //ket->A=A*Hconjugate(Vdagger)*invS_*invS_;
    //bra->A=B*U;
    //righttermK=auxS*auxS*Vdagger;
    //righttermB=Hconjugate(U);
  }
  else{ // No dim cut by Xi (also if both were)
    // WARNING: If both are cut?
    ket->A=A*Hconjugate(Vdagger)*invS_;
    bra->A=B*U*conjugate(invS_);
    righttermK=auxS*Vdagger;
    righttermB=conjugate(auxS)*Hconjugate(U);
  }
  ket->A.reshape(dimA);bra->A.reshape(dimB);
}

void Site::gaugeL(Site* ket,Site* bra,mwArray& lefttermK,
		  mwArray& lefttermB,const mwArray& righttermK,
		  const mwArray& righttermB){
  Indices transpose(2,1); // vector (2,1) for transp.
  // Multiply terms from previous site
  ket->reduceDimL(righttermK);
  bra->reduceDimL(righttermB);
  // We want to impose Sum A B+=1 => contract
  Indices dimA=ket->A.getDimensions();
  Indices dimB=bra->A.getDimensions();
  Indices dim2A(ket->getDl(),(ket->getd())*ket->getDr());
  //dim2A[0]=(ket->getd())*(ket->getDl());dim2A[1]=(ket->getDr());
  Indices dim2B(bra->getDl(),(bra->getd())*bra->getDr());
  //dim2B[0]=(bra->getd())*(bra->getDl());dim2B[1]=(bra->getDr());
  mwArray A=permute(ket->A,Indices(2,1,3));A.reshape(dim2A);
  mwArray B=permute(bra->A,Indices(2,1,3));B.reshape(dim2B);
  mwArray M;
  M=A*Hconjugate(B); // Dla x Dlb
  mwArray U,S,Vdagger;
  double tol=1E-10;int Xi=0;
  wrapper::svd(M,0.,Xi,U,S,Vdagger);
  mwArray auxS=sqrt(S);
  mwArray invS_=invertDiag(auxS,Xi,tol); // pseudoinverse
  double trSre,trSim;
  auxS.trace(trSre,trSim);
  auxS=(1./trSre)*auxS;
  
  //cout<<"Inside gaugeL->auxS="<<auxS<<"\n\tinvS_="<<invS_<<endl;
  // Q: the next assertions do not seem very clear now
  if(A.getDimension(0)>Xi&&B.getDimension(0)<=Xi){
    cout<<"BIG ERROR!"<<endl;
    exit(212);
    //ket->A=A*Hconjugate(Vdagger);
    //bra->A=B*U*conjugate(invS_)*conjugate(invS_);
    //righttermK=Vdagger;
    //righttermB=conjugate(auxS)*conjugate(auxS)*Hconjugate(U);
  }
  else if(A.getDimension(0)<=Xi&&B.getDimension(0)>Xi){
    cout<<"BIG ERROR!"<<endl;
    exit(212);
    //ket->A=A*Hconjugate(Vdagger)*invS_*invS_;
    //bra->A=B*U;
    //righttermK=auxS*auxS*Vdagger;
    //righttermB=Hconjugate(U);
  }
  else{ // No dim cut by Xi (also if both were)
    // WARNING: If both are cut?
    A.multiplyLeft(invS_*Hconjugate(U));
    A.reshape(Indices(dimA[1],dimA[0],dimA[2]));A.permute(Indices(2,1,3));    
    B.multiplyLeft(conjugate(invS_)*Vdagger);
    B.reshape(Indices(dimB[1],dimB[0],dimB[2]));B.permute(Indices(2,1,3));
    lefttermK=U*auxS;
    lefttermB=Hconjugate(Vdagger)*conjugate(auxS);
  }
  ket->A=A;
  bra->A=B;
}

void Site::changeDimensions(int d_,int Dl_,int Dr_,double noise){
  Indices dims(d_,Dl_,Dr_);
  //cout<<"Site changing the dimensions of A:"<<A<<", to "<<dims<<endl;
  A.resize(dims,noise);
  //cout<<"Site changed the dimensions of A:"<<A<<endl;
}

void Site::increaseBondDim(int Dl_,int Dr_){
  changeDimensions(A.getDimension(0),Dl_,Dr_);
}

void Site::increaseBondDimWithNoise(int Dl_,int Dr_,double noise){
  changeDimensions(A.getDimension(0),Dl_,Dr_,noise);
}


void Site::decreaseBondDim(int Dl_,int Dr_){
  changeDimensions(A.getDimension(0),Dl_,Dr_);
}  

void Site::increasePhysDim(int d_){
  if(d_>getd()){
    // cout<<"WARNING: Site increasing physical dimension "<<getd()
    // 	<<"->"<<d_<<endl;
    changeDimensions(d_,A.getDimension(1),A.getDimension(2));
    //    changeDimensions(getd(),A.getDimension(1),A.getDimension(2));
  }
}
void Site::decreasePhysDim(int d_){
  if(d_<getd())
    // cout<<"WARNING: Site decreasing physical dimension "<<getd()
    // 	<<"->"<<d_<<endl;
    changeDimensions(d_,A.getDimension(1),A.getDimension(2));
  //    changeDimensions(getd(),A.getDimension(1),A.getDimension(2));
}

void Site::load(ifstream& data){
  data.read((char*)&pos,sizeof(int));
  if(!data){
    cout<<"Error trying to read Site from file: could not read position"<<endl;
    exit(1);
  }
  A.load(data);
}

void Site::save(ofstream& data) const{
  data.write((char*)&pos,sizeof(int));
  if(!data){
    cout<<"Error trying to write Site("<<pos
	<<") to file: could not write position"<<endl;
    exit(1);
  }
  A.save(data);
}

void Site::loadtext(ifstream& data){
  data>>pos;
  A.loadtext(data);
}

void Site::savetext(ofstream& data) const{
  data<<pos<<" ";
  A.savetext(data);
}


// Testing

void Site::gaugeRT(mwArray& rightterm,const mwArray& leftterm,mwArray& lambda,int cutD){
  cout<<"ERROR! gaugeRT not ready!! Needs to be written as gaugeLT below!"<<endl;exit(1);
  reduceDimR(leftterm);
  int d=A.getDimension(0);int Dl=A.getDimension(1);int Dr=A.getDimension(2);
  A.permute(Indices(2,1,3)); // Dl d Dr
  A.reshape(Indices(Dl,d*Dr));
  mwArray lambdaOld(lambda);
  lambda.multiplyRight(A);
  lambda.reshape(Indices(Dl*d,Dr));
  A.reshape(Indices(Dl*d,Dr));
  lambda.multiplyLeft(permute(A,Indices(2,1))); //Dr Dr
  mwArray U,S,Vd;
  double tol=0.;
  wrapper::svd(lambda,tol,cutD,U,S,Vd);
  lambda=S;
  rightterm=Vd;
  Vd.Hconjugate();
  A.multiplyRight(Vd);
  Dr=A.getDimension(1); // new Dr
  A.reshape(Indices(Dl,d,Dr));
  lambda=lambdaOld;
  lambda.multiplyRight(reshape(A,Indices(Dl,d*Dr)));
  lambda.reshape(Indices(Dl*d,Dr));
  lambda.multiplyLeft(permute(reshape(conjugate(A),Indices(Dl*d,Dr)),Indices(2,1))); //Dr Dr
  // And diagonalize
  vector<complex_t> Dval;
  wrapper::eig(lambda,Dval,U,true);
  lambda=diag(Dval);
  A.reshape(Indices(Dl*d,Dr));
  A.multiplyRight(U);A.reshape(Indices(Dl,d,Dr));
  U.Hconjugate();
  rightterm.multiplyLeft(U);

  A.permute(Indices(2,1,3)); // d Dl Dr
  //A.reshape(Indices(d,Dl,Dr));
  // if(Dr==1&&rightterm.getDimension(1)==1){//last site => do not include phase in the tensor
  //   complex_t aux=rightterm.getElement(0); // single element
  //   complex_t phase=(1./abs(aux))*aux;
  //   A.multiplyRight(mwArray(phase));
  //   rightterm.setElement((complex_t){abs(aux),0.},0);
  // }
}

void Site::gaugeLT(mwArray& leftterm,const mwArray& rightterm,mwArray& lambda,int cutD){
  //cout<<"Site::gaugeLT("<<pos<<") received rightterm:"<<rightterm.getDimensions()<<" lambda "<<lambda<<", A:"<<A.getDimensions()<<endl;
  reduceDimL(rightterm);
  //Since A is dxDlxDr, I have to permute indices
  Indices neworder(2,1,3);
  A.permute(neworder);
  Indices dims=A.getDimensions(); // now Dl, d, Dr
  //mwArray lambdaOld(lambda);
  // multiply the lambda
  lambda.multiplyLeft(reshape(A,Indices(dims[0]*dims[1],dims[2])));
  // And reshape to Dlx(d Dr)
  lambda.reshape(Indices(dims[0],dims[1]*dims[2]));
  mwArray lambdaOld(lambda);
  // multiply by transpose of A
  lambda.multiplyRight(permute(reshape(A,Indices(dims[0],dims[1]*dims[2])),
			       Indices(2,1))); // now Dl x Dl
  //complex_t errSym;Indices exIn;
  //(lambda-permute(lambda,Indices(2,1))).getMaxElement(errSym,exIn);
  //cout<<"Is the new Lambda symmetric? "<<errSym<<endl;
  mwArray U,S,Vd;
  double tol=0.; //1E-12;//0.;
  wrapper::svd(lambda,tol,cutD,U,S,Vd); // returns at most cutD sing vals
  // if(pos<=1){
  //   cout<<"gaugeLT("<<pos<<") S: "<<S.getDimensions()<<endl;
  //   cout<<"lambda:"<<lambda<<", U:"<<U<<endl;
  // }
  // U is the proper isometry to reduce the dim of the next one
  // But I should fix the phases! // 
  // conj(Vd)*U should be it, but if there are degeneracies, there
  // would be a block => instead, make U+ lambda U* real (it is already diagonal)
  vector<complex_t> phases;
  (permute(conjugate(U),Indices(2,1))*lambda*conjugate(U)).getDiagonal(phases);
  // at this point, phases containes the diagonal elements WITH phases squared
  for(int k=0;k<phases.size();k++)
    if(abs(phases[k])>1E-20)
      phases[k]=sqrt(phases[k]/abs(phases[k]));
    else
      phases[k]=ONE_c;
  U=U*diag(phases);
  //cout<<"Now phases are "<<phases<<" and modified lambda "<<permute(conjugate(U),Indices(2,1))*lambda*conjugate(U)<<endl;
  // It may be convenient to force Lambdab to be prop to identity, by setting U=U*sqrt(S) and A->S^{-1/2}*U (U containing the phases)
  
  leftterm=U; 
  A.reshape(Indices(dims[0],dims[1]*dims[2])); // Dl x d*Dr
  U.Hconjugate();A.multiplyLeft(U);
  dims[0]=A.getDimension(0); // new Dl
  lambda=S;
  A.reshape(dims);
  A.permute(neworder);
  // // if(Dl==1&&leftterm.getDimension(0)==1){//fist site => do not include phase in the tensor
  // //   complex_t aux=leftterm.getElement(0); // single element
  // //   complex_t phase=(1./abs(aux))*aux;
  // //   A.multiplyRight(mwArray(phase));
  // //   leftterm.setElement((complex_t){abs(aux),0.},0);
  // // }
}

// void Site::gaugeLT_v0(mwArray& leftterm,const mwArray& rightterm,mwArray& lambda,int cutD){
//   //cout<<"Site::gaugeLT("<<pos<<") received rightterm:"<<rightterm.getDimensions()<<" lambda "<<lambda<<", A:"<<A.getDimensions()<<endl;
//   reduceDimL(rightterm);
//   //Since A is dxDlxDr, I have to permute indices
//   Indices neworder(2,1,3);
//   A.permute(neworder);
//   Indices dims=A.getDimensions(); // now Dl, d, Dr
//   //mwArray lambdaOld(lambda);
//   // multiply the lambda
//   lambda.multiplyLeft(reshape(A,Indices(dims[0]*dims[1],dims[2])));
//   // And reshape to Dlx(d Dr)
//   lambda.reshape(Indices(dims[0],dims[1]*dims[2]));
//   mwArray lambdaOld(lambda);
//   // multiply by transpose of A
//   lambda.multiplyRight(permute(reshape(A,Indices(dims[0],dims[1]*dims[2])),
// 			       Indices(2,1))); // now Dl x Dl
//   //complex_t errSym;Indices exIn;
//   //(lambda-permute(lambda,Indices(2,1))).getMaxElement(errSym,exIn);
//   //cout<<"Is the new Lambda symmetric? "<<errSym<<endl;
//   mwArray U,S,Vd;
//   double tol=1E-12;//0.;
//   wrapper::svd(lambda,tol,cutD,U,S,Vd); // returns at most cutD sing vals
//   // U is the proper isometry to reduce the dim of the next one
//   // But I should fix the phases! // 
//   // conj(Vd)*U should be it, but if there are degeneracies, there
//   // would be a block => instead, make U+ lambda U* real (it is already diagonal)
//   vector<complex_t> phases;
//   (permute(conjugate(U),Indices(2,1))*lambda*conjugate(U)).getDiagonal(phases);
//   // at this point, phases containes the diagonal elements WITH phases squared
  
//   // mwArray expIPhi=conjugate(Vd)*U; // should be the phase diagonal exp(-i Phi)
//   // cout<<"lambda:"<<lambda<<"U:"<<U<<" S:"<<S<<" Vd:"<<Vd<<" Vt*U:"<<expIPhi<<endl;
//   // vector<complex_t> phases;expIPhi.getDiagonal(phases);
//   // {complex_t elem;Indices inds;
//   //   expIPhi.reshape(Indices(phases.size(),phases.size()));
//   //   (diag(phases)-expIPhi).getMaxElement(elem,inds);
//   //   if(abs(elem)>1E-12){
//   //     cout<<"WARNING: Vt*U not diagonal:"<<expIPhi<<endl;
//   //     cout<<"Check Udag*lambda*U*"<<permute(conjugate(U),Indices(2,1))*lambda*conjugate(U)<<endl;
//   //     exit(1);

//   //   }
//   // }
//   for(int k=0;k<phases.size();k++)
//     if(abs(phases[k])>1E-10)
//       phases[k]=sqrt(phases[k]/abs(phases[k]));
//     else
//       phases[k]=ONE_c;
//   U=U*diag(phases);
//   cout<<"Now phases are "<<phases<<" and modified lambda "<<permute(conjugate(U),Indices(2,1))*lambda*conjugate(U)<<endl;
//   leftterm=U; 
//   A.reshape(Indices(dims[0],dims[1]*dims[2])); // Dl x d*Dr
//   U.Hconjugate();A.multiplyLeft(U);
//   dims[0]=A.getDimension(0); // new Dl
//   //  A.reshape(dims);
//   lambda=lambdaOld;
//   //  lambda.multiplyLeft(reshape(A,Indices(dims[0]*dims[1],dims[2])));
//   // And reshape to Dlx(d Dr)
//   //lambda.reshape(Indices(dims[0],dims[1]*dims[2]));
//   // multiply by transpose of A
//   lambda.multiplyRight(permute(reshape(conjugate(A),Indices(dims[0],dims[1]*dims[2])),
//   			       Indices(2,1))); // now Dl x Dl
//   lambda.multiplyLeft(U);
//   // And diagonalize lambda
//   vector<complex_t> Dval;
//   wrapper::eig(lambda,Dval,U,true);
//   lambda=diag(Dval);
//   leftterm.multiplyRight(U);
//   U.Hconjugate();
//   A.reshape(Indices(dims[0],dims[1]*dims[2]));A.multiplyLeft(U);
//   A.reshape(dims);
//   A.permute(neworder);
//   // if(Dl==1&&leftterm.getDimension(0)==1){//fist site => do not include phase in the tensor
//   //   complex_t aux=leftterm.getElement(0); // single element
//   //   complex_t phase=(1./abs(aux))*aux;
//   //   A.multiplyRight(mwArray(phase));
//   //   leftterm.setElement((complex_t){abs(aux),0.},0);
//   // }
// }

