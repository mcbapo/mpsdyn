/**
   \file MPO.cpp
   Basic class containing a Matrix Product Operator
   
   \author Mari-Carmen Banuls
   \date 16/03/2011
*/
#include "MPO.h"
#include "JoinedOperator.h"

using namespace shrt;

MPO::MPO(){
  Ops=0;
  myown=0;
  myops=0;
}

MPO::MPO(int len):length(len),Ops(NULL),myown(1),myops(NULL){
  Ops=new const Operator*[length];
  for(int k=0;k<length;k++)
    Ops[k]=0;
  myown=0;
  myops=0;
//  cout<<"Created Operator Row "<<this<<" of length="<<length<<" and bond="<<bond<<endl;
}

MPO::MPO(int len,const Operator** op):
  length(len),Ops(NULL),myown(1),myops(NULL){
  Ops=new const Operator*[length];
  for(int k=0;k<length;k++){
    Ops[k]=op[k];
  }
  myown=0;
  myops=0;// nothing is mine
}

MPO::MPO(int nr,const MPO* oprs[]):
  length(oprs[0]->getLength()),Ops(NULL),myown(1),myops(NULL){  
  //  bond=(int)pow((double)oprs[0]->getBond(),(double)nr); // temptative
  //Ops=new const JoinedOperator*[length];
  Ops=new const Operator*[length];
  myown=1; // everything is mine: creating the Ops
  myops=new bool[length];
  for(int pos=0;pos<length;pos++){
    // Construct the array of Operators to be fused in the new one
    const Operator** loc=new const Operator*[nr];
    for(int j=0;j<nr;j++){
      loc[j]=&(oprs[j]->getOp(pos));
    }
    Ops[pos]=new Operator(nr,loc);
    //Ops[pos]=new JoinedOperator(nr,loc);
    myops[pos]=1;
    delete []loc;
  }
}

void MPO::clear(){
  if(myown){
    for(int k=0;k<length;k++){
      if(myops[k]){
	if(Ops[k]!=0) delete Ops[k];
      }
    }
    delete []myops;
    myown=0;
  }
  if(Ops!=0)
    delete []Ops;
  Ops=0;
}

void MPO::initLength(int len){
  clear();
  length=len;
  Ops=new const Operator*[length];
  for(int k=0;k<length;k++)
    Ops[k]=0;
  myown=0;
  myops=0;  
}

MPO::~MPO(){
  clear();
}

MPO& MPO::operator=(const MPO& opr){
  if(this!=&opr){
    cout<<"Copy-assigning MPO!"<<endl;
    if(length!=opr.length){
      if(Ops!=0) delete []Ops;
      //      Ops=new const Operator*[length];
      cout<<"WARNING! Copy constructor MPO!"<<endl;
      Ops=new const Operator*[opr.length];
    }
    length=opr.length;
    for(int k=0;k<length;k++){
      Ops[k]=opr.Ops[k];
    }
    if(myown) delete []myops;
    myown=0; // nothing is mine
  }
  return *this;
}


const Operator& MPO::getOp(int k) const{
  //cout<<"MPO: retrieving "<<k<<"="<<*Ops[k]<<endl;
  if(Ops[k]==0){
    cout<<"Error: trying to retrieve empty Operator at "<<k<<endl;
    exit(2);
  }
  return *Ops[k];
}

void MPO::deleteOp(int k){
  if(myown){
    myops[k]=0;
    if(Ops[k]!=0)
      delete Ops[k];
  }
}

const mwArray MPO::getOpData(int pos) const{
  return Ops[pos]->getFullData();
}

void MPO::setOp(int pos,const Operator* op,bool myop_){
  // TODO: Check dimension compatibility wrt neighbours
  if(pos<0||pos>=length){
    cout<<"Error: MPO cannot set Operator at pos "<<pos<<endl;
    exit(212);
  }
  if(myown)
    if(myops[pos]){
      delete Ops[pos];
      myops[pos]=0;
    }
  Ops[pos]=op;
  if(myop_)
    markAsOwn(pos);
}

void MPO::setConjugatedOp(int pos,const Operator* op){
  if(pos<0||pos>=length){
    cout<<"Error: MPO cannot set Operator at pos "<<pos<<endl;
    exit(212);
  }
  if(myown)
    if(myops[pos]){
      delete Ops[pos];myops[pos]=0;
    }
  Operator* auxOp=new Operator(*op);
  auxOp->conjugateOp();
  Ops[pos]=(const Operator*)auxOp;
  markAsOwn(pos);
}

void MPO::markAsOwn(int pos){
  if(!myown){
    myown=1;
    myops=new bool[length];
    for(int i=0;i<length;i++)
      myops[i]=0;
  }
  myops[pos]=1;
}

void MPO::setRotatedOp(int pos,const Operator* op,
		       const Indices& neworder,bool conjugate){
  // cout<<"MPO::setRotatedOp("<<pos<<") from Operator: "
//       <<op->data_->GetDimensions()
//       <<" permuted to "<<newdims<<endl;
  setRotatedOp(pos,op->getFullData(),neworder,conjugate);
}

void MPO::setRotatedOp(int pos,const mwArray& data,
		       const Indices& neworder,bool conjugate){
  // Since I only keep the pointer to Operator, creation of new Op is 
  // done here
  // cout<<"MPO::setRotatedOp("<<pos<<") from mwArray "
  //        <<data.GetDimensions()
  //        <<" permuted to "<<newdims<<endl;
  // I need to remove the previous one, if it was there and mine!!!
  if(myown)
    if(myops[pos]){
      if(Ops[pos]!=0)
	delete Ops[pos];
    }
  Ops[pos]=new Operator(data,neworder,conjugate);
  //   cout<<"Warning: setRotatedOp creating a new Operator may cause memory "
  //       <<"leak. To avoid it, the operator should be deleted by deleteOp("
  //       <<pos<<")!"<<endl;
  markAsOwn(pos);
}


ostream& operator<<(ostream& os,const MPO& ops){
  os<<"MPO len="<<ops.getLength()<<" :";
  for(int k=0;k<ops.getLength();k++){
    os<<"pos["<<k<<"]: ";
    if(ops.isEmpty(k)) os<<"-empty-"<<endl;
    else os<<ops.getOp(k).getDimensions()<<endl;
    //else os<<ops.getOp(k)<<endl;
  }
  return os;
}

MPO* MPO::getBasicRow() const{
  cout<<"WARNING!! getBasicRow at MPO"<<endl;
  return (MPO*)this;
}

void MPO::saveMPO(const char* filename) const{
  /** TODO!!!!!! Proper save for different operators!! */
  ofstream outfile(filename,ios::binary);
  if(!outfile.is_open()){
    cout<<"Error: unable to open file "<<filename<<" to write MPO"<<endl;
    exit(2);
  }
  else{
    outfile.write((char*)&length,sizeof(int));
    if(!outfile){
      cout<<"Error trying to write MPO length to file "<<filename<<endl;
      exit(1);
    }
    for(int k=0;k<length;k++)
      Ops[k]->getFullData().save(outfile);
    outfile.close();
  }
}

#include <iomanip>
void MPO::exportMPOtext(const char* filename) const{
  ofstream outfile(filename);
  if(!outfile.is_open()){
    cout<<"Error: unable to open file "<<filename<<" to write"<<endl;
    exit(2);
  }
  else{
    outfile<<length<<" ";
    //outfile<<setprecision(10);
    for(int k=0;k<length;k++){
      //Ops[k]->savetext(outfile);
      Ops[k]->getFullData().savetext(outfile);
    }
    outfile.close();
 }
}
void MPO::importMPOtext(const char* filename){
  ifstream infile(filename);
  if(!infile.is_open()){
    cout<<"Error: unable to open file "<<filename<<" to read (text)"<<endl;
    exit(2);
  }
  else{
    clear(); // discard what I had before
    infile>>length;
    Ops=new const Operator*[length];
    myown=1;myops=new bool[length];
    for(int k=0;k<length;k++){
      mwArray aux;
      aux.loadtext(infile);
      Ops[k]=new Operator(aux);
      //Ops[k]=new Operator(1,1,1,1);
      //Ops[k]->loadtext(infile);
      myops[k]=1;
    }
    infile.close();
 }
}

void MPO::exportForMatlab(const char* filename,int precision) const{
  ofstream outfile(filename);
  if(!outfile.is_open()){
    cout<<"Error: unable to open file "<<filename<<" to write"<<endl;
    exit(2);
  }
  else{
    if(precision>0)
      outfile << setprecision(precision);
    outfile<<"opers={};clear aux;"<<endl;
    //outfile<<setprecision(10);
    for(int k=0;k<length;k++){
      //Ops[k]->savetext(outfile);
      putForMatlab(outfile,Ops[k]->getFullData(),"aux");
      outfile<<"opers{"<<k+1<<"}=aux;clear aux;"<<endl;
    }
    outfile.close();
 }
}

void MPO::loadMPO(const char* filename) {
  ifstream infile(filename,ios::binary);
  if(!infile.is_open()){
    cout<<"Error: unable to open file "<<filename<<" to read MPO from file "<<filename<<endl;
    exit(2);
  }
  else{
    clear();
    infile.read((char*)&length,sizeof(int));
    if(!infile){
      cout<<"Error trying to read MPO length from file!"<<endl;
      exit(1);
    }
    Ops=new const Operator*[length];
    myown=1;
    myops=new bool[length];
    for(int l=0;l<length;l++){
      Ops[l]=new Operator(mwArray(infile));
      myops[l]=1;
    }
  }
  infile.close();
}

const vector<int> MPO::getDimensions() const{
  //cout<<"MPO::getDimensions() for "<<*this<<endl;
  vector<int> dims(length,0);
  for(int k=0;k<length;k++){
    dims[k]=Ops[k]->getd();
    // cout<<"Ops["<<k<<"], dim="<<dims[k]<<endl;
  }
  return dims;
}

const vector<int> MPO::getOriginalDimensions() const{
  //cout<<"MPO::getOriginalDimensions() for "<<*this<<endl;
  vector<int> dims(length,0);
  for(int k=0;k<length;k++){
    dims[k]=Ops[k]->getdorig();
    //cout<<"Ops["<<k<<"],"<<Ops[k]->getInnerData().GetDimensions()
    //<<", d2="<<dims[k]<<endl;
  }
  return dims;
}

void MPO::setOperatorArray(const vector<const mwArray*>& opers){
  int len=opers.size();
//   cout<<"readOperArray read opers of size "<<size
//       <<", length="<<len<<endl;
  clear();
  length=len;
  myown=1;
  myops=new bool[length];
  Ops=new const Operator*[length];
  for(int l=0;l<length;l++){
    //cout<<"setOperArray created Operator for site "<<l<<" with dimensions "
    //<<d1<<", "<<d2<<", "<<d3<<", "<<d4<<endl;
    Ops[l]=new Operator(*opers[l]);
    myops[l]=1;
  }
}

#include "FoldedOperator.h"

void MPO::fold(MPO& newR) const {
  // Site 0 is now the former operator in the middle
  cout<<"Basic MPO::fold(), original length="<<length<<endl;
  bool isodd=length%2!=0;
  // length of the new chain
  //int newlen=length%2==0?length/2:(length+1)/2;
  int newlen=isodd?(length+1)/2:length/2;

  cout<<"new length="<<newlen<<endl;
  if(isodd) // if it was full, restart now with the proper length
    newR.initLength(newlen); 
  else
    newR.initLength(newlen+1); // auxiliary site 0 with only contraction

  // what is to be substracted from middlept, with k, for folding pair
  int middlept=length/2;
  int offset=isodd?0:1;
  int initpos=0; // Pos in the folded chain
  if(isodd){// the first operator is singular: only the middle one with rearranged indices  
    mwArray aux=getOpData(length/2);
    const Operator& op=getOp(length/2);
    aux.permute(Indices(1,3,2,4));
    aux.reshape(Indices(op.getd(),1,op.getdorig(),op.getDl()*op.getDr()));
    newR.setOp(0,new Operator(aux),1);
    cout<<"Odd case: Operator 0 set to middle op"<<endl;
    initpos++;
  }
  else{
    // In the even case, the "middle" operator is substituted by folding two operators in 
    // the center, but also contracting their common index
    // I will do this with an extra site which contracts both, and has all remaining 
    // dimensions equal to one.
    mwArray aux=identityMatrix(getOp(length/2).getDl());
    aux.reshape(Indices(1,1,1,-1));
    newR.setOp(0,new Operator(aux),1);
    cout<<"Even case: Operator 0 set to contraction"<<endl;
  }
  for(int k=initpos;k<newlen;k++){
    int pos1=middlept-k-offset; // left operator
    int pos2=middlept+k; // right operator
    int d=getOp(pos1).getd();
    int newpos=isodd?k:k+1;
    newR.setOp(newpos,new FoldedOperator(getOpData(pos1),getOpData(pos2)),1);
    //cout<<"Operator "<<k<<" set from "<<pos1<<" and "<<pos2<<endl;
  }
}


void MPO::join(int nr,const MPO* oprs[],MPO& joined){
  int len=oprs[0]->getLength();
  //cout<<"Basic MPO::join(), nr="<<nr<<", len="<<len<<endl;
  joined.initLength(len);
  for(int pos=0;pos<len;pos++){
    // Construct the array of Operators to be fused in the new one
    const Operator** loc=new const Operator*[nr];
    for(int j=0;j<nr;j++){
      loc[j]=&(oprs[j]->getOp(pos));
    }
    joined.setOp(pos,new JoinedOperator(nr,loc),true);
    //cout<<"Set JoinedOperator at pos "<<pos<<" dims "
    //	<<joined.getOp(pos).getDimensions()<<endl;
    delete []loc;
  }
}

void expandOper(const MPO& mpo,mwArray& oper){
  oper=mwArray(ONE_c);
  int dphys1=1; // temporary output dimensions 
  int dphys2=1; // temporary input (right) dimensions 
  int Do=1; // open bond dimension
  int N=mpo.getLength();
  oper.reshape(Indices(dphys1*dphys2,Do));
  for(int k=0;k<N;k++){
    mwArray aux=mpo.getOp(k).getFullData();
    //cout<<"Multiplying oper{"<<k<<"}="<<aux<<endl;
    Indices dims=aux.getDimensions();
    int d1=dims[0];int Dl=dims[1];int d2=dims[2];int Dr=dims[3];
    if(Dl!=Do){
      cout<<"Error: Dimensions do not agree in expandOper!"<<endl;
      exit(1);
    }
    aux.permute(Indices(2,1,3,4));
    aux.reshape(Indices(Dl,d1*d2*Dr));
    oper.multiplyRight(aux);
    oper.reshape(Indices(dphys1,dphys2,d1,d2,Dr));
    oper.permute(Indices(1,3,2,4,5));
    dphys1=dphys1*d1;
    dphys2=dphys2*d2;
    Do=Dr;
    oper.reshape(Indices(dphys1*dphys2,Do));
  }
  // Now, if Do at the end is not one, I do not get an operator, but
  // still hae one contractable index. This may be desired, so I allow
  // it
  if(Do==1)
    oper.reshape(Indices(dphys1,dphys2));
  else
    oper.reshape(Indices(dphys1,dphys2,Do));
}
