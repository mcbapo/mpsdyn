#include "FoldedOperator.h"

using namespace std;
using namespace shrt;

FoldedOperator::FoldedOperator(const mwArray& dataL_,
			       const mwArray& dataR,bool norev):
  dataL(NULL),mydataL(0){
  if(dataR.getRank()!=4||dataL_.getRank()!=4){
    cout<<"Error: invalid dimensions for the creation of FoldedOperator, "
	<<"dataL="<<dataL_<<", dataR="<<dataR<<endl;
  }
  data=new mwArray(dataR);
  if(norev)
    dataL=new mwArray(dataL_);
  else
    dataL=new mwArray(permute(dataL_,Indices(1,4,3,2)));
  updateDims();
  mydata=true;mydataL=true;
  //cout<<"Constructed foldedOp with dims "<<dims<<endl;
}

void FoldedOperator::updateDims(){
  dims.clear();
  for(int k=0;k<4;k++)
    dims.push_back(data->getDimension(k)*dataL->getDimension(k));
}

FoldedOperator::FoldedOperator(const FoldedOperator& op):
  Operator(op),dataL(NULL),mydataL(0){
  dataL=new mwArray(*op.dataL);
  updateDims();
  mydataL=true;
}

FoldedOperator& FoldedOperator::operator=(const FoldedOperator& op){
  if(this!=&op){
    clear();
    this->Operator::operator=(op);
    dataL=new mwArray(*op.dataL);
    mydataL=true;
    updateDims();
  }
  return *this;
}

void FoldedOperator::clear(){
  if(mydataL) if(dataL!=0) delete dataL;
  Operator::clear();
  dataL=0;
}

FoldedOperator::~FoldedOperator(){
  clear();
}

void FoldedOperator::conjugateOp(){
  Operator::conjugateOp(); // conjugate the right one
  if(dataL!=0){
    mwArray* auxL=new mwArray(*dataL);
    auxL->conjugate();
    if(mydataL) delete dataL;
    mydataL=true;dataL=auxL;
  }
}

void FoldedOperator::permuteOp(Indices perm){
  Operator::permuteOp(perm);
  if(dataL!=0){
    mwArray* auxL=new mwArray(*dataL);
    auxL->permute(perm);
    if(mydataL) delete dataL;
    mydataL=true;dataL=auxL;
  }
  updateDims();
}

void FoldedOperator::put(ostream& os) const{
  os<<"FoldedOperator with components L:"
    <<dataL->getDimensions()
    <<" and R:"<<data->getDimensions();
}

//ostream& operator<<(ostream& os,const FoldedOperator& oper){
//  return os<<"FoldedOperator with components L"
//	   <<oper.dataL->getDimensions()
//	   <<" and R"<<oper.data->getDimensions();
// /*if(oper.dataL!=0||oper.data!=0)
//  return os<<*(oper.dataL)<<*(oper.data);
//else
//return os<<"FoldedOperator empty";*/
//}



void FoldedOperator::contractRight(const FoldedOperator* rightOp){
  mwArray* data_[]={new mwArray(*dataL),new mwArray(*data)};
  mwArray dataR[]={rightOp->getDataLeft(),rightOp->getDataRight()};
  for(int k=0;k<2;k++){
    Indices dimsLoc=data_[k]->getDimensions();
    Indices dimsR=dataR[k].getDimensions();
    dataR[k].permute(Indices(2,1,3,4));
    dataR[k].reshape(Indices(dimsR[1],-1));
    data_[k]->reshape(Indices(-1,dimsLoc[3]));
    data_[k]->multiplyRight(dataR[k]);
    data_[k]->reshape(Indices(dimsLoc[0],dimsLoc[1],dimsLoc[2],dimsR[0],
			    dimsR[2],dimsR[3]));
    data_[k]->permute(Indices(1,4,2,3,5,6));
    data_[k]->reshape(Indices(dimsLoc[0]*dimsR[0],dimsLoc[1],
			    dimsLoc[2]*dimsR[2],dimsR[3]));
}
  if(mydata) delete data;
  if(mydataL) delete dataL;
  data=data_[1];dataL=data_[0];
  updateDims();
}

// #######
// When I just want to store the pointer to operR (L is copied)
void FoldedOperator::setData(const mwArray* operL,const mwArray* operR){
  Indices dimsL=operL->getDimensions();
  Indices dimsR=operR->getDimensions();
  if((dimsL.size()<2)||(dimsR.size()<2)||(dimsL[0]!=dimsR[0])||
     (dimsL[1]!=dimsR[3])||(dimsL[3]!=dimsR[3])||
     (dimsL[3]!=dimsR[1])){
    cout<<"Error: cannot set FoldedOperator to data with dimensions L:"
	<<dimsL<<", dimensions R:"<<dimsR<<endl;
    exit(2);
  }
  clear();
  data=operR; // not copied
  mydataL=1;
  dataL=new mwArray(permute(*operL,Indices(1,4,3,2))); mydataL=1; // copied
  cout<<"FoldedOperator holding alien data "<<endl;
  updateDims();
}

const mwArray FoldedOperator::getFullData() const{
  int r1=data->getDimension(0);int r2=data->getDimension(1);
  int r3=data->getDimension(2);int r4=data->getDimension(3);
  int l1=dataL->getDimension(0);int l2=dataL->getDimension(1);
  int l3=dataL->getDimension(2);int l4=dataL->getDimension(3);
  mwArray op=reshape(*dataL,Indices(l1*l2*l3*l4,1))*
    reshape(*data,Indices(1,r1*r2*r3*r4));
  op=permute(reshape(op,Indices(l1,l2,l3,l4,r1,r2,r3,r4)),
	     Indices(1,5,2,6,3,7,4,8));
  op=reshape(op,Indices(l1*r1,l2*r2,l3*r3,l4*r4));
  //  cout<<"Using specialized FoldedOperator::getFullData()! for this "<<*this;
  //cout<<" gives op "<<op.getDimensions()<<endl;
  return op;
}

void FoldedOperator::contractL(mwArray& result,const mwArray& termL,
			       const Site& ket,const Site& bra,
			       bool dagger) const {
  contractleftfop(result,termL,ket,bra,dagger);
}


void FoldedOperator::contractR(mwArray& result,const mwArray& termR,
			       const Site& ket,const Site& bra,
			       bool dagger) const {
  contractrightfop(result,termR,ket,bra,dagger);
}

void FoldedOperator::contractMket(mwArray& result,const mwArray& termL,
				  const Site& ket,const mwArray& termR,
				  bool dagger) const {
  //  cout<<"FoldedOperator::contractMket, this="<<*this<<endl;
  contractMketfop(result,termL,ket,termR,dagger);
}
 
void FoldedOperator::contractMbra(mwArray& result,const mwArray& termL,
				  const Site& bra,const mwArray& termR,
				  bool dagger) const {
  //cout<<"FoldedOperator::contractMbra, this="<<*this<<endl;
  contractMbrafop(result,termL,bra,termR,dagger);
}

void FoldedOperator::contractN(mwArray& result,const mwArray& termL,
			       const mwArray& termR,bool dagger) const {
  contractNfop(result,termL,termR,dagger);
}
 
void FoldedOperator::contractleftfop(mwArray& result,const mwArray& termL,
				     const Site& ket,
				     const Site& bra,bool dagger) const{
  int d2=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int dp2=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int betab=termL.getDimension(0);
  int alpha2=termL.getDimension(1);  
  int betak=termL.getDimension(2);
  int dul=dagger?dataL->getDimension(2):dataL->getDimension(0);
  int al1l=dataL->getDimension(1);
  int ddl=dagger?dataL->getDimension(0):dataL->getDimension(2);
  int al2l=dataL->getDimension(3);
  int dur=dagger?data->getDimension(2):data->getDimension(0);
  int al1r=data->getDimension(1);
  int ddr=dagger?data->getDimension(0):data->getDimension(2);
  int al2r=data->getDimension(3);
  const mwArray* operL=dataL;
  const mwArray* operR=data;
  if(dagger){
    operL=new mwArray(permute(conjugate(*dataL),Indices(3,2,1,4)));
    operR=new mwArray(permute(conjugate(*data),Indices(3,2,1,4)));
  }
#ifdef CHECKDIMS
  if((d2!=dul*dur)||(dp2!=ddl*ddr)||(Dlp!=betak)||(Dl!=betab)||
     (al1l*al1r!=alpha2)){
    //||(ddl!=ddr)||(dul!=dur)||(al1l!=al1r)||(al2l!=al2r)){
    cout<<"FoldedOperator::contractleftfop, wrong dimensions, for site "
	<<ket.getPos()
	<<", ket "<<ket.getA().getDimensions()
	<<", bra "<<bra.getA().getDimensions()
	<<", termL "<<termL.getDimensions()
	<<", operL "<<dataL->getDimensions()
	<<", operR "<<data->getDimensions()<<endl;
    exit(212);
  }
#endif
  if((dp2>Dlp||d2>Dl)){
    //First contract Aket with operL
    mwArray tmpres=reshape(ket.getA(),Indices(ddl,ddr,Dlp,Drp));
    tmpres.permute(Indices(2,3,4,1));
    tmpres.reshape(Indices(ddr*Dlp*Drp,ddl));
    mwArray aux=permute(*operL,Indices(3,1,2,4));aux.reshape(Indices(ddl,dul*al1l*al2l));
    tmpres=tmpres*aux;
    //Then with operR
    tmpres.reshape(Indices(ddr,Dlp,Drp,dul,al1l,al2l));
    tmpres.permute(Indices(2,3,4,5,6,1));
    tmpres.reshape(Indices(-1,ddr));
    aux=permute(*operR,Indices(3,1,2,4));aux.reshape(Indices(ddr,-1));
    tmpres=tmpres*aux;
    tmpres.reshape(Indices(Dlp,Drp,dul,al1l,al2l,dur,al1r,al2r));
    tmpres.permute(Indices(1,2,4,5,7,8,3,6));
    tmpres.reshape(Indices(-1,dul*dur));
    //and with Abra
    aux=conjugate(bra.getA());aux.reshape(Indices(dul,dur,Dl,Dr));aux.reshape(Indices(dul*dur,-1));
    result=tmpres*aux;
    //Finally contract everything with termL
    result.reshape(Indices(Dlp,Drp,al1l,al2l,al1r,al2r,Dl,Dr));
    result.permute(Indices(2,4,6,8,1,3,5,7));
    result.reshape(Indices(-1,Dlp*al1l*al1r*Dl));
    aux=permute(termL,Indices(3,2,1));aux.reshape(Indices(betab*al1l*al1r*betak,1));
    result=result*aux;
    result.reshape(Indices(Drp,al2l,al2r,Dr));
    result.permute(Indices(4,2,3,1));
    result.reshape(Indices(Dr,al2l*al2r,Drp));   
  }
  else{
    // First step: contract Aket with tmpL
    mwArray tmpres=permute(ket.getA(),Indices(2,1,3));tmpres.reshape(Indices(Dlp,dp2*Drp));
    tmpres=reshape(termL,Indices(betab*alpha2,betak))*tmpres;
    // Now contract Abra
    mwArray aux=conjugate(bra.getA());aux.permute(Indices(3,1,2));
    aux.reshape(Indices(Dr*d2,Dl));
    tmpres.reshape(Indices(betab,alpha2*dp2*Drp));
    tmpres=aux*tmpres;
    // Now contract operR
    tmpres.reshape(Indices(Dr,dul,dur,al1l,al1r,ddl,ddr,Drp));
    tmpres.permute(Indices(1,8,2,4,6,3,5,7));
    tmpres.reshape(Indices(Dr*Drp*dul*al1l*ddl,dur*al1r*ddr));
    tmpres=tmpres*reshape(*operR,Indices(dur*al1r*ddr,al2r));
    // Now operL
    tmpres.reshape(Indices(Dr*Drp,dul*al1l*ddl,al2r));
    tmpres.permute(Indices(1,3,2));
    tmpres.reshape(Indices(Dr*Drp*al2r,dul*al1l*ddl));
    result=tmpres*reshape(*operL,Indices(dul*al1l*ddl,al2l));
    result.reshape(Indices(Dr,Drp,al2r,al2l));
    result.permute(Indices(1,4,3,2));
    result.reshape(Indices(Dr,al2l*al2r,Drp));
  }
  if(dagger){
    delete operR,operL;
  }
}

void FoldedOperator::contractrightfop(mwArray& result,const mwArray& termR,
				     const Site& ket,
				     const Site& bra,bool dagger) const{
  int d2=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int dp2=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int betab=termR.getDimension(0);
  int alpha2=termR.getDimension(1);  
  int betak=termR.getDimension(2);
  int d1l=dagger?dataL->getDimension(2):dataL->getDimension(0);
  int al1l=dataL->getDimension(1);
  int d2l=dagger?dataL->getDimension(0):dataL->getDimension(2);
  int al2l=dataL->getDimension(3);
  int d1r=dagger?data->getDimension(2):data->getDimension(0);
  int al1r=data->getDimension(1);
  int d2r=dagger?data->getDimension(0):data->getDimension(2);
  int al2r=data->getDimension(3);
  const mwArray* operL=dataL;
  const mwArray* operR=data;
  if(dagger){
    operL=new mwArray(permute(conjugate(*dataL),Indices(3,2,1,4)));
    operR=new mwArray(permute(conjugate(*data),Indices(3,2,1,4)));
  }
#ifdef CHECKDIMS
  if((dp2!=d2l*d2r)||(d2!=d1l*d1r)||(Drp!=betak)||(Dr!=betab)||
     (al2l*al2r!=alpha2)){
       //||(d2l!=d2r)||(d1l!=d1r)||(al1l!=al1r)||(al2l!=al2r)){
    cout<<"FoldedOperator::contractrightfop, wrong dimensions, "
	<<"ket "<<ket.getA().getDimensions()
	<<"bra "<<bra.getA().getDimensions()
	<<"termR "<<termR.getDimensions()
	<<"operL "<<dataL->getDimensions()
	<<"operR "<<data->getDimensions()<<endl;
    exit(212);
  }
#endif

  if((dp2>Drp||d2>Dr)){
    //First step, contract Aket with operL
    mwArray tmpres=ket.getA();
    tmpres.reshape(Indices(d2l,d2r,Dlp,Drp));tmpres.permute(Indices(2,3,4,1));
    tmpres.reshape(Indices(d2r*Dlp*Drp,d2l));
    mwArray aux=*operL;aux.permute(Indices(3,1,2,4));aux.reshape(Indices(d2l,d1l*al1l*al2l));
    tmpres=tmpres*aux;			  
    //Then with operR
    tmpres.reshape(Indices(d2r,Dlp,Drp,d1l,al1l,al2l));
    tmpres.permute(Indices(2,3,4,5,6,1));
    tmpres.reshape(Indices(-1,d2r));
    tmpres=tmpres*reshape(permute(*operR,Indices(3,1,2,4)),Indices(d2r,-1));
    tmpres.reshape(Indices(Dlp,Drp,d1l,al1l,al2l,d1r,al1r,al2r));
    tmpres.permute(Indices(1,2,4,5,7,8,3,6));
    tmpres.reshape(Indices(-1,d1l*d1r));
    //Then Contract everything with Abra
    aux=bra.getA();aux.conjugate();aux.reshape(Indices(d1l,d1r,Dl,Dr));aux.reshape(Indices(d1l*d1r,-1));
    result=tmpres*aux;
    //Finally contract with termR
    result.reshape(Indices(Dlp,Drp,al1l,al2l,al1r,al2r,Dl,Dr));
    result.permute(Indices(1,3,5,7,2,4,6,8));
    result.reshape(Indices(-1,Drp*al2l*al2r*Dr));
    aux=termR;aux.permute(Indices(3,2,1));aux.reshape(Indices(betab*al2l*al2r*betak,1));
    result = result*aux;
    result.reshape(Indices(Dlp,al1l,al1r,Dl));
    result.permute(Indices(4,2,3,1));
    result.reshape(Indices(Dl,al1l*al1r,Dlp));
  }
  else{
    // First step, contract Aket with tmpR
    mwArray tmpres=termR;tmpres.reshape(Indices(betab*alpha2,betak));
    mwArray aux=ket.getA();aux.permute(Indices(3,1,2));aux.reshape(Indices(Drp,dp2*Dlp));
    tmpres=tmpres*aux;
    // Then Abra
    tmpres.reshape(Indices(betab,alpha2*dp2*Dlp));
    aux=bra.getA();aux.conjugate();aux.reshape(Indices(d2*Dl,Dr));
    tmpres=aux*tmpres;
    // Now operL
    tmpres.reshape(Indices(d1l,d1r,Dl,al2l,al2r,d2l,d2r,Dlp));
    tmpres.permute(Indices(1,4,6,2,5,7,3,8));
    tmpres.reshape(Indices(d1l*al2l*d2l,d1r*al2r*d2r*Dl*Dlp));
    aux=*operL;aux.permute(Indices(2,1,4,3));aux.reshape(Indices(al1l,d1l*al2l*d2l));
    tmpres=aux*tmpres;
    // And operR
    tmpres.reshape(Indices(al1l,d1r*al2r*d2r,Dl*Dlp));
    tmpres.permute(Indices(2,1,3));
    tmpres.reshape(Indices(d1r*al2r*d2r,al1l*Dl*Dlp));
    aux=*operR;aux.permute(Indices(2,1,4,3));aux.reshape(Indices(al1r,d1r*al2r*d2r));
    result=aux*tmpres;
    result.reshape(Indices(al1r,al1l,Dl,Dlp));
    result.permute(Indices(3,2,1,4));
    result.reshape(Indices(Dl,al1l*al1r,Dlp));
  }
  if(dagger) delete operR,operL;
}

void FoldedOperator::contractMketfop(mwArray& result,const mwArray& termL,
				     const Site& ket,
				     const mwArray& termR,bool dagger) const{

  int dp2=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int betabr=termR.getDimension(0);
  int alphar=termR.getDimension(1);  
  int betakr=termR.getDimension(2);
  int betabl=termL.getDimension(0);
  int alphal=termL.getDimension(1);  
  int betakl=termL.getDimension(2);
  int dul=dagger?dataL->getDimension(2):dataL->getDimension(0);
  int al1l=dataL->getDimension(1);
  int ddl=dagger?dataL->getDimension(0):dataL->getDimension(2);
  int al2l=dataL->getDimension(3);
  int dur=dagger?data->getDimension(2):data->getDimension(0);
  int al1r=data->getDimension(1);
  int ddr=dagger?data->getDimension(0):data->getDimension(2);
  int al2r=data->getDimension(3);
  const mwArray* operL=dataL;
  const mwArray* operR=data;
  if(dagger){
    operL=new mwArray(permute(conjugate(*dataL),Indices(3,2,1,4)));
    operR=new mwArray(permute(conjugate(*data),Indices(3,2,1,4)));
  }
#ifdef CHECKDIMS
  if((dp2!=ddl*ddr)||(betakl!=Dlp)||(betakr!=Drp)||(al1l*al1r!=alphal)||
            (al2l*al2r!=alphar)){
    //||(ddl!=ddr)||(dul!=dur)||(al1l!=al1r)||(al2l!=al2r)){
    cout<<"FoldedOperator::contractMketfop, wrong dimensions, "
	<<"ket "<<ket.getA().getDimensions()
	<<"termL "<<termL.getDimensions()
	<<"termR "<<termR.getDimensions()
	<<"operL "<<dataL->getDimensions()
	<<"operR "<<data->getDimensions()<<endl;
    exit(212);
}
#endif
  // First step: contract Aket with operR
  mwArray tmpres=*operR;tmpres.permute(Indices(1,2,4,3));tmpres.reshape(Indices(dur*al1r*al2r,ddr));
  mwArray aux=ket.getA();aux.reshape(Indices(ddl,ddr,Dlp,Drp));
  aux.permute(Indices(2,1,3,4));aux.reshape(Indices(ddr,ddl*Dlp*Drp));
  tmpres=tmpres*aux;
  // Now, termR
  tmpres.reshape(Indices(dur*al1r,al2r,ddl*Dlp,Drp));
  tmpres.permute(Indices(1,3,2,4));
  tmpres.reshape(Indices(dur*al1r*ddl*Dlp,al2r*Drp));
  aux=termR;aux.reshape(Indices(betabr*al2l,al2r*betakr));
  aux.permute(Indices(2,1)); 
  tmpres=tmpres*aux;
  // Contract operL
  tmpres.reshape(Indices(dur,al1r,ddl,Dlp,betabr,al2l));
  tmpres.permute(Indices(3,6,1,2,4,5));
  tmpres.reshape(Indices(ddl*al2l,dur*al1r*Dlp*betabr));
  tmpres=reshape(*operL,Indices(dul*al1l,ddl*al2l))*tmpres;
  // finally, termL
  tmpres=reshape(termL,Indices(betabl,al1l*al1r*betakl))*
    reshape(permute(reshape(tmpres,Indices(dul,al1l,dur,al1r*Dlp,betabr)),
		    Indices(2,4,1,3,5)),Indices(al1l*al1r*Dlp,dul*dur*betabr));
  // and reshape properly
  result=reshape(permute(reshape(tmpres,Indices(betabl,dul*dur,betabr)),
			 Indices(2,1,3)),Indices(dul*dur*betabl*betabr,1));
  if(dagger) delete operR,operL;
}

void FoldedOperator::contractMbrafop(mwArray& result,const mwArray& termL,
				     const Site& bra,
				     const mwArray& termR,bool dagger) const{
  int dp2=bra.getd();int Dlp=bra.getDl();int Drp=bra.getDr();
  int betabr=termR.getDimension(0);
  int alphar=termR.getDimension(1);  
  int betakr=termR.getDimension(2);
  int betabl=termL.getDimension(0);
  int alphal=termL.getDimension(1);  
  int betakl=termL.getDimension(2);
  int dul=dagger?dataL->getDimension(2):dataL->getDimension(0);
  int al1l=dataL->getDimension(1);
  int ddl=dagger?dataL->getDimension(0):dataL->getDimension(2);
  int al2l=dataL->getDimension(3);
  int dur=dagger?data->getDimension(2):data->getDimension(0);
  int al1r=data->getDimension(1);
  int ddr=dagger?data->getDimension(0):data->getDimension(2);
  int al2r=data->getDimension(3);
  const mwArray* operL=dataL;
  const mwArray* operR=data;
  if(dagger){
    operL=new mwArray(permute(conjugate(*dataL),Indices(3,2,1,4)));
    operR=new mwArray(permute(conjugate(*data),Indices(3,2,1,4)));
  }
#ifdef CHECKDIMS
  if((dp2!=dul*dur)||(betabl!=Dlp)||(betabr!=Drp)||(al1l*al1r!=alphal)||
     (al2l*al2r!=alphar)){
    //||(ddl!=ddr)||(dul!=dur)||(al1l!=al1r)||(al2l!=al2r)){
    cout<<"FoldedOperator::contractMbrafop, wrong dimensions, "
	<<"bra "<<bra.getA().getDimensions()
	<<"termL "<<termL.getDimensions()
	<<"termR "<<termR.getDimensions()
	<<"operL "<<dataL->getDimensions()
	<<"operR "<<data->getDimensions()<<endl;
    exit(212);
}
#endif
  // first: Abra with operL
  mwArray tmpres=*operL;
  tmpres.permute(Indices(2,3,4,1));
  tmpres.reshape(Indices(al1l*ddl*al2l,dul));
  mwArray aux=bra.getA();aux.conjugate();aux.reshape(Indices(dul,dur*Dlp*Drp));
  tmpres=tmpres*aux;
  // Second step: termL
  tmpres.reshape(Indices(al1l,ddl*al2l*dur,betabl,betabr));
  tmpres.permute(Indices(2,4,3,1));
  tmpres.reshape(Indices(ddl*al2l*dur*betabr,betabl*al1l));
  tmpres=tmpres*reshape(termL,Indices(betabl*al1l,al1r*betakl));
  // Now operR
  tmpres.reshape(Indices(ddl*al2l,dur,betabr,al1r,betakl));
  tmpres.permute(Indices(1,3,5,2,4));
  tmpres.reshape(Indices(ddl*al2l*betabr*betakl,dur*al1r));
  tmpres=tmpres*reshape(*operR,Indices(dur*al1r,ddr*al2r));
  // and termR
  tmpres.reshape(Indices(ddl,al2l,betabr,betakl,ddr,al2r));
  tmpres.permute(Indices(1,5,4,3,2,6));
  tmpres.reshape(Indices(ddl*ddr*betakl,betabr*al2l*al2r));
  result=tmpres*reshape(termR,Indices(betabr*al2l*al2r,betakr));
  // Finally, reshape 
  result.reshape(Indices(ddl*ddr*betakl*betakr,1));
  if(dagger) delete operR,operL;
}

void FoldedOperator::contractNfop(mwArray& result,const mwArray& termL,
				  const mwArray& termR,bool dagger) const{
  int betabr=termR.getDimension(0);
  int alphar=termR.getDimension(1);  
  int betakr=termR.getDimension(2);
  int betabl=termL.getDimension(0);
  int alphal=termL.getDimension(1);  
  int betakl=termL.getDimension(2);
  int dul=dagger?dataL->getDimension(2):dataL->getDimension(0);
  int al1l=dataL->getDimension(1);
  int ddl=dagger?dataL->getDimension(0):dataL->getDimension(2);
  int al2l=dataL->getDimension(3);
  int dur=dagger?data->getDimension(2):data->getDimension(0);
  int al1r=data->getDimension(1);
  int ddr=dagger?data->getDimension(0):data->getDimension(2);
  int al2r=data->getDimension(3);
  const mwArray* operL=dataL;
  const mwArray* operR=data;
  if(dagger){
    operL=new mwArray(permute(conjugate(*dataL),Indices(3,2,1,4)));
    operR=new mwArray(permute(conjugate(*data),Indices(3,2,1,4)));
  }
#ifdef CHECKDIMS
  if((al1l*al1r!=alphal)||(al2l*al2r!=alphar)){
    //||(ddl!=ddr)||(dul!=dur)||(al1l!=al1r)||(al2l!=al2r))  {
    cout<<"FoldedOperator::contractNfop, wrong dimensions, "
	<<"termL "<<termL.getDimensions()
	<<"termR "<<termR.getDimensions()
	<<"operL "<<dataL->getDimensions()
	<<"operR "<<data->getDimensions()<<endl;
    exit(212);
}
#endif
  // Contract termL and operL
  mwArray tmpres=termL;
  tmpres.reshape(Indices(betabl,al1l,al1r*betakl));
  tmpres.permute(Indices(1,3,2));
  tmpres.reshape(Indices(betabl*al1r*betakl,al1l));
  mwArray aux=*operL;
  aux.permute(Indices(2,1,3,4));
  aux.reshape(Indices(al1l,dul*ddl*al2l));
  tmpres=tmpres*aux;
  tmpres.reshape(Indices(betabl,al1r,betakl*dul*ddl,al2l));
  tmpres.permute(Indices(1,3,2,4));
  tmpres.reshape(Indices(betabl*betakl*dul*ddl,al1r*al2l));
  // contract termR and operR

  aux=termR;
  aux.reshape(Indices(betabr*al2l,al2r,betakr));
  aux.permute(Indices(2,1,3));
  aux.reshape(Indices(al2r,betabr*al2l*betakr));

  mwArray tmpres2=reshape(*operR,Indices(dur*al1r*ddr,al2r))*aux;
  
  tmpres2.reshape(Indices(dur,al1r,ddr*betabr,al2l,betakr));
  tmpres2.permute(Indices(2,4,1,3,5));
  tmpres2.reshape(Indices(al1r*al2l,dur*ddr*betabr*betakr));

  // Now contract both together
  result=tmpres*tmpres2;
  result.reshape(Indices(betabl,betakl,dul,ddl,dur,ddr,betabr,betakr));
  result.permute(Indices(3,5,1,7,4,6,2,8));
  result.reshape(Indices(dul*dur*betabl*betabr,ddl*ddr*betakl*betakr));

  if(dagger) delete operR,operL;
}
