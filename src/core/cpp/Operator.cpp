/**
   \file Operator.cpp 

   Implementation of the basic class that provides an Operator to be part
   of a MPO, and act on MPS. 
   
   \author Mari-Carmen Banuls
   \date 16/03/2011
*/

#include "math.h"
#include "Operator.h"

using namespace shrt;
using namespace std;

Operator::Operator(int d,int Dl,int Dr,int d2):dims(),data(NULL),mydata(1){
  dims.push_back(d);
  dims.push_back(Dl);
  dims.push_back(d);
  dims.push_back(Dr);
  if(d2!=0) dims[2]=d2;
  // TODO: Default is Identity??
}

Operator::Operator(const Indices& dims_):dims(dims_),data(NULL),mydata(1){
  if(dims.size()>4){
    cout<<"Error: cannot create Operator from mwArray of dims "<<dims
	<<endl;
    exit(2);
  }    
  // If rank was smaller than 4, fill in with ones
  while(dims.size()<4) dims.push_back(1);
}


Operator::Operator(const mwArray& data_):
  dims(data_.getDimensions()),data(NULL),mydata(1){
  if(data_.getRank()>4){
    cout<<"Error: cannot create Operator from mwArray of dims "<<dims
	<<endl;
    exit(2);
  }    
  // If rank was smaller than 4, fill in with ones
  while(dims.size()<4) dims.push_back(1);
  data=new mwArray(reshape(data_,dims));
  mydata=1;
}

Operator::Operator(const mwArray& oper,const Indices& neworder,
		   bool conjugate_):dims(),data(NULL),mydata(1){
  if((oper.getRank()!=neworder.size())||oper.getRank()>4){
    cout<<"Error: incompatible dimensions of oper="<<oper<<" and "
	<<"permutation vector "<<neworder<<" at Operator::Operator"<<endl;
    exit(212);
  }
  (permute(oper,neworder)).getDimensions(dims);
  while(dims.size()<4) dims.push_back(1);
  if(conjugate_) 
    data=new mwArray(reshape(permute(conjugate(oper),neworder),dims));
  else
    data=new mwArray(reshape(permute(oper,neworder),dims));
}

int Operator::getInnerDimension(int k)const{
  if(data==0) return 0;
  if(k<0||k>=4){
    cout<<"Error: cannot retrieve dimension "<<k<<" for an Operator"<<endl;
  }
  return data->getDimension(k);
}


// TODO!! Replace by invocation to JoinedOperator FullData??
// Or remove entirely!!
Operator::Operator(int nr,const Operator* oprs[]):
  dims(),data(NULL),mydata(1){
  cout<<"WARNING WARNING WARNING!!!!!!!"<<endl;
  cout<<"Calling Operator constructor from array of Operators!!! "<<endl;
  cout<<"Should probably be using JoinedOperator instead!"<<endl;
  int dd=oprs[0]->getInnerDimension(0);
  int Dl=oprs[0]->getInnerDimension(1);
  int ddp=oprs[0]->getInnerDimension(2);
  int Dr=oprs[0]->getInnerDimension(3);
  Indices curdims(dd,Dl,Dr);int prod=dd*Dl*Dr;
  mwArray tmp=permute(oprs[0]->getFullData(),Indices(1,2,4,3));
  tmp=reshape(tmp,Indices(prod,ddp));
  int dimL=Dl;int dimR=Dr; //dims to l and r so far
  int ddn,Dln,ddpn,Drn;
  for(int k=1;k<nr;k++){
    ddn=oprs[k]->getInnerDimension(0);
    Dln=oprs[k]->getInnerDimension(1);
    ddpn=oprs[k]->getInnerDimension(2);
    Drn=oprs[k]->getInnerDimension(3);
#ifdef CHECKDIMS
    if(ddn!=ddp){
      cout<<"Error: incompatible dimensions in Operator::Operator(int,"
	  <<"Operator*[]), on term "<<k<<endl;
      exit(2);
    }
    else ddp=ddpn;
#endif
    mwArray aux=permute(oprs[k]->getFullData(),Indices(1,2,4,3));
    aux=reshape(aux,Indices(ddn,Dln*Drn*ddpn));
    tmp=tmp*aux;
    curdims=Indices(curdims,Indices(Dln,Drn));prod=prod*Dln*Drn;
    tmp=reshape(tmp,Indices(prod,ddpn));
    dimL=dimL*Dln;dimR=dimR*Drn;
  }
  tmp.reshape(Indices(curdims,ddpn)); // dd1,Dl1,Dr1,Dl2,Dr2,...ddpn
  Indices permutation(1);
  int last=curdims.size()+1;
  for(int l=2;l<=last;l=l+2) permutation.push_back(l);
  for(int l=3;l<=last-1;l=l+2) permutation.push_back(l);
  tmp.permute(permutation);
  data=new mwArray(reshape(tmp,Indices(dd,dimL,ddpn,dimR)));
  mydata=1;
  dims=data->getDimensions();
}

Operator::Operator():dims(),data(NULL),mydata(1){};

Operator::~Operator(){
  clear();
}

void Operator::clear(){
  dims.clear();
  if(mydata)
    if(data!=0) delete data;
  mydata=0;
}

/* Protected version of the constructor */
Operator::Operator(const mwArray& data_,int d,int Dl,int Dr,int d2):
  dims(),data(NULL),mydata(1){
  dims=Indices(d,Dl,d2==0?d:d2,Dr); 
  data=new mwArray(data_);
  mydata=1;
  cout<<"Derived class creating Operator from mwArray "<<data_
      <<" with fake dims "<<dims<<endl;
}

Operator::Operator(const Operator& op):
  dims(op.getDimensions()),data(NULL),mydata(1){
  data=new mwArray(op.getFullData());
  mydata=1;
}

Operator& Operator::operator=(const Operator& op){
  if(this!=&op){
    clear();
    dims=op.getDimensions();
    mydata=op.mydata;
    if(mydata){
      data=new mwArray(op.getFullData());
    }
    else{
      data=op.data; // just pointer
    }
  cout<<"WARNING: Trying to copy basic Operator!"<<endl;
  //exit(222);
  }
  return *this;
}


const Operator* Operator::identity(int d){
  mwArray ident=identityMatrix(d);
  return new Operator(reshape(ident,Indices(d,1,d,1)));
}

void Operator::contractRight(const Operator* rightOp){
  if(rightOp->getDl()!=getDr()){
    cout<<"ERROR: Cannot contract the right index of "<<*this
	<<" with the left one of "<<*rightOp<<endl;
  }
  mwArray dataR(rightOp->getFullData());
  Indices dimsR=dataR.getDimensions();
  dataR.permute(Indices(2,1,3,4));
  dataR.reshape(Indices(dimsR[1],-1));
  mwArray* dataux=new mwArray(*data);
  dataux->reshape(Indices(-1,dims[3]));
  dataux->multiplyRight(dataR);
  dataux->reshape(Indices(dims[0],dims[1],dims[2],dimsR[0],dimsR[2],dimsR[3]));
  dataux->permute(Indices(1,4,2,3,5,6));
  dataux->reshape(Indices(dims[0]*dimsR[0],dims[1],dims[2]*dimsR[2],dimsR[3]));
  dims=dataux->getDimensions();
  if(mydata) delete data;
  data=dataux;mydata=1;
}

void Operator::permuteOp(shrt::Indices perm){
  if(perm.size()!=4){
    cout<<"ERROR: permute can only handle list of four indices!!"<<endl;
    exit(212);
  }
  mwArray* aux=new mwArray(*data);
  aux->permute(perm);
  if(mydata){
    delete data;
  }
  mydata=1;
  data=aux;
  dims=data->getDimensions();
}

void Operator::setData(const mwArray* oper){
  if(oper->getRank()!=4){
    cout<<"Error: Trying to setData in operator from mwArray pointer "
	<<" to "<<*oper<<endl;
    exit(2);
  }   
  clear();
  dims=oper->getDimensions();
  data=oper;
  mydata=0;
  //    cout<<"Operator set to data ptr "<<data<<endl;
}

/*************************
void Operator::copyData(const mwArray& oper){
  if(oper.getRank()!=4){
    cout<<"Error: Trying to copyData in operator from mwArray of dimensions "
	<<oper.getDimensions()<<", i.e. rank="<<oper.getRank()<<endl;
    exit(2);
  }   
  clear();
  data=new mwArray(oper);
  dims=data->getDimensions();
  mydata=1;
  cout<<"WARNING: copyData invoked from basic Operator"<<endl;
}
*/
void Operator::conjugateOp(){
  mwArray* aux=new mwArray(*data);
  aux->conjugate();
  if(mydata){
    delete data;
  }
  mydata=1;
  data=aux;
}

void Operator::setRotatedData(const mwArray& oper,const Indices& neworder,
			      bool conjugate_){
  if((oper.getRank()!=neworder.size())||oper.getRank()>4){
    cout<<"Error: incompatible dimensions of oper="<<oper<<" and "
	<<"permutation vector "<<neworder<<" at Operator::setRotatedData()"
	<<endl;
    exit(212);
  }
  clear();
  (permute(oper,neworder)).getDimensions(dims);
  mydata=1;
  while(dims.size()<4) dims.push_back(1);
  if(conjugate_) 
    data=new mwArray(reshape(permute(conjugate(oper),neworder),dims));
  else
    data=new mwArray(reshape(permute(oper,neworder),dims));
}


/* Protected version, ignoring dimensions of the data */
void Operator::setData(const mwArray* oper,int d,int Dl,int Dr){
  clear();
  dims=Indices(d,Dl,d,Dr); // assume initial and final dim are equal
  data=oper;
  mydata=0;
  //    cout<<"Operator set to ptr to "<<data<<" and dims "<<dims<<endl;
}

void Operator::contractleftop(mwArray& result,const mwArray& termL,
			      const Site& ket,
			      const mwArray& oper,const Site& bra) const {
  int d=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int betab=termL.getDimension(0);int alpha=termL.getDimension(1);
  int betak=termL.getDimension(2);
  int du=oper.getDimension(0);int al1=oper.getDimension(1);
  int dd=oper.getDimension(2);int al2=oper.getDimension(3);
#ifdef CHECKDIMS
  if((d!=du)||(dp!=dd)||(Dlp!=betak)||(Dl!=betab)||(al1!=alpha)){
    cout<<"Wrong dimensions in Operator::contractleftop: bra "
	<<bra.getDimensions()
	<<", ket "<<ket.getDimensions()<<", termL "<<termL.getDimensions()
	<<", oper "<<oper.getDimensions()<<endl;
    exit(212);
  }
#endif
  // TODO!!!! Change order, as this is suboptimal!!
  // result=reshape(termL,Indices(betab*alpha,betak));
  // mwArray aux=permute(ket.getA(),Indices(2,1,3));aux.reshape(Indices(Dlp,dp*Drp));
  // result=result*aux;
  // aux=permute(bra.getA(),Indices(3,1,2));aux.conjugate();
  // aux.reshape(Indices(Dr*d,Dl));
  // result.reshape(Indices(betab,alpha*dp*Drp));
  // result=aux*result;
  // result.reshape(Indices(Dr,d*alpha*dp,Drp));
  // result.permute(Indices(1,3,2));
  // result.reshape(Indices(Dr*Drp,d*alpha*dp));
  // // Now contract oper
  // result=result*reshape(oper,Indices(du*al1*dd,al2));
  // result.reshape(Indices(Dr,Drp,al2));
  // result.permute(Indices(1,3,2));

  // better order of contractions
  result=permute(ket.getA(),Indices(2,1,3));result.reshape(Indices(Dlp,dp*Drp));
  result.multiplyLeft(reshape(termL,Indices(betab*alpha,betak))); // betab*alpha x dp*Drp
  // now oper
  result.reshape(Indices(betab,alpha*dp,Drp));
  result.permute(Indices(2,1,3));result.reshape(Indices(alpha*dp,betab*Drp));
  result.multiplyLeft(reshape(
  			      permute(reshape(oper,Indices(du,al1*dd,al2)),
  				      Indices(1,3,2)),
  			      Indices(du*al2,al1*dd))); // du*al2 betab*Drp
  result.reshape(Indices(du,al2,betab,Drp));
  result.permute(Indices(1,3,2,4));
  result.reshape(Indices(du*betab,al2*Drp));
  result.multiplyLeft(Hconjugate(reshape(bra.getA(),Indices(d*Dl,Dr)))); // Dr x al2*Drp
  result.reshape(Indices(Dr,al2,Drp));
}

void Operator::contractL(mwArray& result,const mwArray& termL,const Site& ket,
			 const Site& bra,bool dagger) const {
  if(!dagger)
    contractleftop(result,termL,ket,getFullData(),bra);
  //    cout<<"Operator "<<data_<<" contractL with ket"<<
  //ket<<" and bra "<<bra<endl;
  else{
    mwArray adjO=permute(getFullData(),Indices(3,2,1,4));adjO.conjugate();
    contractleftop(result,termL,ket,adjO,bra);
  }
}

void Operator::contractrightop(mwArray& result,const mwArray& termR,const Site& ket,
			      const mwArray& oper,const Site& bra) const {
  int d=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int betab=termR.getDimension(0);int alpha=termR.getDimension(1);
  int betak=termR.getDimension(2);
  int d1=oper.getDimension(0);int al1=oper.getDimension(1);
  int d2=oper.getDimension(2);int al2=oper.getDimension(3);
#ifdef CHECKDIMS
  if((dp!=d2)||(d!=d1)||(Drp!=betak)||(Dr!=betab)||(al2!=alpha)){
    cout<<"Wrong dimensions in Operator::contractrightop: bra "
	<<bra.getDimensions()<<", ket "
	<<ket.getDimensions()<<", termL "<<termR.getDimensions()
	<<", oper "<<oper.getDimensions()<<endl;
    exit(212);
  }
#endif
  // Notice, the original order of contraction was not optimal!!!! I keep it here, just in case some bug is discovered (3.11.2021)
  // // First step, contract ket with tmpR
  // mwArray aux=permute(ket.getA(),Indices(3,1,2));aux.reshape(Indices(Drp,dp*Dlp));
  // result=reshape(termR,Indices(betab*alpha,betak))*aux;
  // // Then Abra
  // aux=reshape(bra.getA(),Indices(d*Dl,Dr));aux.conjugate();
  // result.reshape(Indices(betab,alpha*dp*Dlp));
  // result=aux*result;
  // // Now oper
  // result.reshape(Indices(d,Dl,alpha*dp,Dlp));
  // result.permute(Indices(1,3,2,4));
  // result.reshape(Indices(d*alpha*dp,Dl*Dlp));
  // result=reshape(permute(oper,Indices(2,1,4,3)),Indices(al1,d1*al2*d2))*
  //   result;
  // result.reshape(Indices(al1,Dl,Dlp));
  // result.permute(Indices(2,1,3));
  
  // First step, contract ket with tmpR
  result=ket.getA();
  result.permute(Indices(3,1,2)); // Drp dp Dlp
  result.reshape(Indices(Drp,dp*Dlp));
  result.multiplyLeft(reshape(termR,Indices(betab*alpha,betak))); // betab*alpha x dp*Dlp
  // Now oper
  result.reshape(Indices(betab,alpha,dp,Dlp));
  result.permute(Indices(3,2,1,4)); // dp alpha betab Dlp
  result.reshape(Indices(dp*alpha,betab*Dlp));
  result.multiplyLeft(reshape(oper,Indices(d1*al1,d2*al2))); // d1*al1 x betab*Dlp
  // Then Abra
  result.reshape(Indices(d1,al1,betab,Dlp));
  result.permute(Indices(1,3,2,4)); // d1 betab al1 Dp
  result.reshape(Indices(d1*betab,al1*Dlp));
  result.multiplyLeft(reshape(permute(conjugate(bra.getA()),Indices(2,1,3)),
			      Indices(Dl,d*Dr))); // Dl (bra) x al1*Dlp
  result.reshape(Indices(Dl,al1,Dlp));
}


void Operator::contractR(mwArray& result,const mwArray& termR,const Site& ket,
			 const Site& bra,bool dagger) const {
  if(!dagger)
    contractrightop(result,termR,ket,getFullData(),bra);
  else{
    mwArray adjO=permute(getFullData(),Indices(3,2,1,4));adjO.conjugate();
    contractrightop(result,termR,ket,adjO,bra);
  }
}

const mwArray Operator::doubleop(bool dagger) const{
  // \todo: save a copy as static, so that I can reuse it?
  mwArray oper=getFullData();
  if(dagger) oper=permute(conjugate(oper),Indices(3,2,1,4));
  int d1=oper.getDimension(0);int al1=oper.getDimension(1);
  int d2=oper.getDimension(2);int al2=oper.getDimension(3);
  mwArray tmpres=permute(oper,Indices(3,2,4,1));tmpres.conjugate();
  tmpres.reshape(Indices(d2*al1*al2,d1));
  tmpres=tmpres*reshape(oper,Indices(d1,al1*d2*al2));
  tmpres.reshape(Indices(d2,al1,al2,al1,d2,al2));
  tmpres.permute(Indices(1,2,4,5,3,6));
  tmpres.reshape(Indices(d2,al1*al1,d2,al2*al2));
  return tmpres;
}

void Operator::contractL2(mwArray& result,const mwArray& termL,const Site& ket,
			  const Site& bra,bool dagger) const {
  mwArray aux=doubleop(dagger);
  contractleftop(result,termL,ket,aux,bra);
}
void Operator::contractR2(mwArray& result,const mwArray& termR,const Site& ket,
			  const Site& bra,bool dagger) const {
  mwArray aux=doubleop(dagger);
  contractrightop(result,termR,ket,aux,bra);
}

/****** AQUI!!!!!!!!!! ********/
void Operator::contractMket(mwArray& result,const mwArray& termL,const Site& ket,
			    const mwArray& termR,bool dagger) const {
  //  cout<<"basic Operator::contractMket of this:"<<*this<<endl;
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int al1=getInnerDimension(1);int al2=getInnerDimension(3);
  int du,dd;
  if(dagger){
    du=getInnerDimension(2);dd=getInnerDimension(0);
  }
  else{
    dd=getInnerDimension(2);du=getInnerDimension(0);
  }
  int betabr=termR.getDimension(0);int alphar=termR.getDimension(1);
  int betakr=termR.getDimension(2);
  int betabl=termL.getDimension(0);int alphal=termL.getDimension(1);
  int betakl=termL.getDimension(2);
#ifdef CHECKDIMS
  if((dp!=dd)||(betakl!=Dlp)||(betakr!=Drp)||(al1!=alphal)||(al2!=alphar)){
    cout<<"Wrong dimensions in contractMketop: ket "<<ket.getDimensions()
	<<", termL "<<termL.getDimensions()
	<<", termR "<<termR.getDimensions()<<", oper "<<data->getDimensions()
	<<", dagger="<<dagger<<endl;
    exit(212);
  }
#endif
  // First step: contract Aket with tmpR
  result=permute(ket.getA(),Indices(2,1,3));result.reshape(Indices(Dlp,dp*Drp));
  result=reshape(termL,Indices(betabl*alphal,betakl))*result;
  result.reshape(Indices(betabl,alphal*dp,Drp));
  result.permute(Indices(1,3,2));
  result.reshape(Indices(betabl*Drp,alphal*dp));
  // Now contract oper
  mwArray aux;
  if(dagger){
    aux=permute(getFullData(),Indices(2,1,3,4));aux.conjugate();
    //tmpres=tmpres*reshape(permute(permute(conjugate(*data),Indices(3,2,1,4)),
    //				  Indices(2,3,1,4)),Indices(al1*dd,du*al2));
  }  
  else{
    aux=permute(getFullData(),Indices(2,3,1,4));
    //tmpres=tmpres*reshape(permute(getFullData(),Indices(2,3,1,4)),Indices(al1*dd,du*al2));
  }
  aux.reshape(Indices(al1*dd,du*al2));
  result=result*aux;
  // Now add termR
  result.reshape(Indices(betabl,Drp,du,al2));
  result.permute(Indices(3,1,4,2));
  result.reshape(Indices(du*betabl,al2*Drp));
  aux=permute(termR,Indices(2,3,1));aux.reshape(Indices(alphar*betakr,betabr));
  result=result*aux;
  result.reshape(Indices(du*betabl*betabr,1));
}


void Operator::contractMbra(mwArray& result,const mwArray& termL,const Site& bra,
			    const mwArray& termR,bool dagger) const {
  //  cout<<"basic Operator::contractMbra on this="<<*this<<endl;
  int dp=bra.getd();int Dlp=bra.getDl();int Drp=bra.getDr();
  int al1=getInnerDimension(1);int al2=getInnerDimension(3);
  int du=getInnerDimension(0);int dd=getInnerDimension(2);
  int betabr=termR.getDimension(0);int alphar=termR.getDimension(1);
  int betakr=termR.getDimension(2);
  int betabl=termL.getDimension(0);int alphal=termL.getDimension(1);
  int betakl=termL.getDimension(2);
#ifdef CHECKDIMS
  if((dp!=du)||(betabl!=Dlp)||(betabr!=Drp)||(al1!=alphal)||(al2!=alphar)){
    cout<<"Wrong dimensions in contractMbraop: bra "<<bra.getDimensions()
	<<", termL "<<termL.getDimensions()
	<<", termR "<<termR.getDimensions()<<", oper "
	<<data->getDimensions()<<endl;
    exit(212);
  }
#endif
  // First step: contract Abra with tmpR
  result=conjugate(bra.getA());result.reshape(Indices(dp*Dlp,Drp));
  result=result*reshape(termR,Indices(betabr,alphar*betakr));
  // Now oper
  //  mwArray aux=permute(getFullData(),Indices(2,3,1,4));
  mwArray aux;
  if(dagger){
    aux=permute(getFullData(),Indices(2,1,3,4));aux.conjugate();
    //tmpres=tmpres*reshape(permute(permute(conjugate(*data),Indices(3,2,1,4)),
    //				  Indices(2,3,1,4)),Indices(al1*dd,du*al2));
  }  
  else{
    aux=permute(getFullData(),Indices(2,3,1,4));
    //tmpres=tmpres*reshape(permute(getFullData(),Indices(2,3,1,4)),Indices(al1*dd,du*al2));
  }
  aux.reshape(Indices(al1*dd,du*al2));
  result.reshape(Indices(dp,Dlp,alphar,betakr));
  result.permute(Indices(1,3,2,4));
  result.reshape(Indices(dp*alphar,Dlp*betakr));
  result=aux*result;
  // finally add termL
  aux=permute(termL,Indices(3,2,1));aux.reshape(Indices(betakl,alphal*betabl));
  result.reshape(Indices(al1,dd,Dlp,betakr));
  result.permute(Indices(1,3,2,4));
  result.reshape(Indices(al1*Dlp,dd*betakr));
  result=aux*result;
  result.reshape(Indices(betakl,dd,betakr));
  result.permute(Indices(2,1,3));
  result.reshape(Indices(dd*betakl*betakr,1));
}

void Operator::contractNop(mwArray& result,const mwArray& termL,
			   const mwArray& oper,
			   const mwArray& termR) const{
  int betab=termL.getDimension(0);int alpha=termL.getDimension(1);
  int betak=termL.getDimension(2);
  int betabr=termR.getDimension(0);int alphar=termR.getDimension(1);
  int betakr=termR.getDimension(2);
  int du=oper.getDimension(0);int al1=oper.getDimension(1);
  int dd=oper.getDimension(2);int al2=oper.getDimension(3);
#ifdef CHECKDIMS
  if((al1!=alpha)||(al2!=alphar)||(betab!=betak)||(betabr!=betakr)){
    cout<<"Wrong dimensions in contractNop: termL "<<termL.getDimensions()
	<<", oper "<<oper.getDimensions()
	<<", termR "<<termR.getDimensions()<<endl;
    exit(212);
  }
#endif
  // First step: contract oper with tmpL
  result=permute(termL,Indices(1,3,2));result.reshape(Indices(betab*betak,alpha));
  mwArray aux=permute(oper,Indices(2,1,3,4));aux.reshape(Indices(al1,du*dd*al2));
  result=result*aux;
  // now contract with tmpR
  result.reshape(Indices(betab*betak*du*dd,al2));
  aux=permute(termR,Indices(2,1,3));aux.reshape(Indices(alphar,betabr*betakr));
  result=result*aux;
  // finally restore indices to their place
  result.reshape(Indices(betab,betak,du,dd,betabr,betakr));
  result.permute(Indices(3,1,5,4,2,6));
  result.reshape(Indices(du*betab*betabr,dd*betak*betakr));
}


void Operator::contractN2(mwArray& result,const mwArray& termL,const mwArray& termR,
			  bool dagger) const {
    contractNop(result,termL,doubleop(dagger),termR);
}

void Operator::contractN(mwArray& result,const mwArray& termL,
			 const mwArray& termR,
			 bool dagger) const {
  if(!dagger)
    contractNop(result,termL,getFullData(),termR);
  else
    contractNop(result,termL,permute(conjugate(getFullData()),
				     Indices(3,2,1,4)),termR);
}

void Operator::put(ostream& os) const{
  os<<"basic Operator: ";
  if(data!=0)
    os<<getDimensions();
    //    os<<*(data);
  else
    os<<"Operator empty";
}

ostream& operator<<(ostream& os,const Operator& oper){
  oper.put(os);
  return os;
  //  if(oper.data!=0)
  // return os<<*(oper.data);
  //else
  //  return os<<"Operator empty";
}

void Operator::savetext(ofstream& outfile) const{
  getFullData().savetext(outfile);
}

void Operator::loadtext(ifstream& infile){
  clear();
  mwArray aux;
  aux.loadtext(infile);
  data=new mwArray(aux);
  dims=data->getDimensions();
  mydata=1;
}
