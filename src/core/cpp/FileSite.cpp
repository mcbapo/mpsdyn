/**
   \file FileSite.cpp
   Implementation of the class that stores site tensors in disk
   
   \author Mari-Carmen Banuls
   \date 31/01/2012
*/

#include "FileSite.h"
#include <cmath>
#include <stdlib.h>
#include <unistd.h>
#define MAXLEN 140

using namespace shrt;

//const char* defaultDir="/ptmp/mpq/banulsm";
bool FileSite::dirset(0);
char* FileSite::dir=0;//uninitialized

void FileSite::setDir(const char* dir_){
  if(!dirset){
    dir=new char[MAXLEN];  
    sprintf(dir,"%s",dir_);
    dirset=1;
    char command[MAXLEN];
    sprintf(command,"mkdir %s",dir);    
    int errcomm=system(command); // created dir
  }
  else{
    cout<<"Trying to set directory for FileSite, "<<
      "when it was already set?"<<endl;
  }
}

void FileSite::initFile(const char* dir_){
  if(!dirset)
    setDir(dir_);
  filename=new char[MAXLEN];  
  if(dir_==0){
    sprintf(filename,"%s/%pXXXXXX",dir,this);    
  }
  else{
    char command[MAXLEN];
    sprintf(command,"mkdir %s",dir_);    
    sprintf(filename,"%s/%pXXXXXX",dir_,this);        
    int errcomm=system(command); // created dir
  }
  int fd=mkstemp(filename); 
  close(fd); // peligro de colision?
  //tmpnam(filename);
  //sprintf(filename,"%s.mat",filename);
  //cout<<"Now my filename is "<<filename<<endl;
}

FileSite::FileSite(int posit,int dphys,int Dleft,int Dright,
		   const char* dir_):
  Site(posit,dphys,Dleft,Dright),inmem(1),dims(dphys,Dleft,Dright),
  filename(NULL){
   // cout<<"Constructed FileSite("<<posit<<","<<dphys<<","<<Dleft<<
   // ","<<Dright<<")- no A"<<endl;
  initFile(dir_);
  write();
}

FileSite::FileSite(int posit,int* dimensions,
		   const char* dir_):
  Site(posit,dimensions),inmem(1),filename(NULL),dims(){
  // cout<<"Constructed FileSite("<<posit<<","<<getDimensions()<<")- no A"<<endl;
  initFile(dir_);
  dims=Site::getDimensions();
  write();
}

FileSite::FileSite(const FileSite& s,int _pos):
  //Site(_pos,s.getA()),inmem(1){
  Site(_pos,s.getA()),inmem(1),dims(s.getDimensions()),filename(NULL){
  if(_pos==0) setPos(s.getPos());
   // cout<<"Constructed FileSite("<<getPos()<<" from copy (original dims:"<<
   //   s.getDimensions()<<", new("<<getDimensions()<<endl;
  initFile();
  write();
}

FileSite::FileSite(int pos_,const mwArray& A_):
  Site(pos_,A_),inmem(1),filename(NULL),dims(){
   // cout<<"Constructed FileSite("<<getPos()<<") at "<<this<<" from A (original:"<<
   //   A<<") new("<<getd()<<
   //   ","<<getDl()<<","<<getDr()<<")"<<endl;
  dims=Site::getDimensions();
  initFile();
  write();
}

FileSite::~FileSite(){
  remove(filename);
  delete []filename;
}


void FileSite::setSite(const mwArray& newA){
  int nrcomp=dims[0]*dims[1]*dims[2];
  if((newA.getRank()<=0||newA.getRank()>3)||
     (newA.getNrComponents()!=nrcomp)){
    cout<<"Error FileSite::setSite trying to set matrix pos("<<pos<<
      ")of different or invalid dimensions!"
      " old:"<<dims<<", new:"<<newA.getDimensions()<<endl;
    exit(2);
  }
  A=reshape(newA,dims);
  inmem=1;
   // cout<<"Set FileSite("<<getPos()<<" to newA:"
   //     <<A.getDimensions()<<endl;
  write();
}

void FileSite::setToZero(){
  A=mwArray(dims);inmem=1;
  A.fillWithZero();
  write();
}

const mwArray FileSite::getA() const{
    if(inmem) return A;
    ifstream infile(filename,ios::binary);
    if(!infile.is_open()){
      cout<<"Error: unable to open file "<<filename<<" for FileSite"<<endl;
      exit(1);
    }
    else{
      mwArray B(infile);
      infile.close();
      return B;
  }
//     cout<<"In getA() of pos "<<getPos()<<", at "<<this<<"("<<getd()<<
//       ","<<getDl()<<","<<getDr()<<"), read from file "<<
//       filename<<" returning "
// 	<<B.GetDimensions()<<endl;

}


void FileSite::setRotatedSite(const mwArray& newA,const Indices& newdims,
			  bool conjugate){
  Site::setRotatedSite(newA,newdims,conjugate);inmem=1;
  dims=Site::getDimensions();
  write();
}


void FileSite::load(ifstream& data){
  //read(); // ignore?
  Site::load(data);inmem=1;
  dims=Site::getDimensions();
  write();
}

void FileSite::save(ofstream& data) const{
   if(inmem)
     Site::save(data);
   else{ // read it
     data.write((char*)&pos,sizeof(int));
     mwArray B=getA();
     B.save(data);
   }
}

void FileSite::loadtext(ifstream& data){
  Site::loadtext(data);inmem=1;
  dims=Site::getDimensions();
  write();
}

void FileSite::savetext(ofstream& data) const{
  if(inmem)
    Site::savetext(data);
  else{ // read it
     data<<pos<<" ";
     mwArray B=getA();
     B.savetext(data);
  }
}

void FileSite::put(ostream& os) const{
  os<<"FileSite("<<getPos()<<"), A dims"<<dims<<", inmem "<<inmem
    <<", A:"<<A.getDimensions()<<endl;
    //<<" values "<<getA()<<endl;
}

// ostream& operator<<(ostream& os,const FileSite& site){
//   return os<<"FileSite("<<site.getPos()<<"), A:"
// 	   <<site.getA()<<endl;
// }

void FileSite::setState(SingleState state){
  //  read(); // need the dimensions!
  A=mwArray(dims);
  Site::setState(state);inmem=1;
  write();
}

FileSite& FileSite::operator=(const FileSite& s){
  if(this!=&s){
    dims=s.getDimensions();
    A=s.getA();
    pos=s.pos;
    //    this->Site::operator=(s);
    //dims=Site::getDimensions();
    write();
    // cout<<"Assigning Site "<<&s<<" to "<<this<<endl;
  }
  return *this;
}

void FileSite::gaugeR(mwArray& rightterm,const mwArray& leftterm,int cutD){
  bool wasinmem=inmem;
  if(!inmem) read();
  Site::gaugeR(rightterm,leftterm,cutD);
  dims=Site::getDimensions();
  if(!wasinmem) write();
}

void FileSite::gaugeL(mwArray& leftterm,const mwArray& rightterm,int cutD){
  //cout<<"FileSite::gaugeL("<<pos<<")"<<endl;
  bool wasinmem=inmem;
  if(!inmem) read();
  Site::gaugeL(leftterm,rightterm,cutD);
  dims=Site::getDimensions();
  if(!wasinmem) write();
}

void FileSite::reduceDimR(const mwArray& leftterm){
  bool wasinmem=inmem;
  if(!inmem) read();
  Site::reduceDimR(leftterm);
  dims=Site::getDimensions();
  if(!wasinmem) write();
}

void FileSite::reduceDimL(const mwArray& rightterm){
  bool wasinmem=inmem;
  if(!inmem) read();
  Site::reduceDimL(rightterm);
  dims=Site::getDimensions();
  if(!wasinmem) write();
}

void FileSite::gaugeLOp(mwArray& leftterm,const mwArray& op,
			const mwArray& rightterm,bool dagger,int cutD){
  bool wasinmem=inmem;
  if(!inmem) read();
  //cout<<"FileSite::gaugeLOp()"<<endl;
  Site::gaugeLOp(leftterm,op,rightterm,dagger,cutD);
  dims=Site::getDimensions();
  if(!wasinmem) write();
}
 
void FileSite::reduceDimLOp(const mwArray& op,const mwArray& rightterm,
			    bool dagger){
  bool wasinmem=inmem;
  if(!inmem) read();
  //cout<<"FileSite "<<*this<<" reduceDimLOp"<<endl;
  Site::reduceDimLOp(op,rightterm,dagger);
  dims=Site::getDimensions();
  if(!wasinmem) write();             
}


void FileSite::changeSign(){
  bool wasinmem=inmem;
  if(!inmem) read();
  Site::changeSign();
  if(!wasinmem) write();
}


void FileSite::gaugeR(FileSite* ket,FileSite* bra,mwArray& righttermK,
		  mwArray& righttermB,const mwArray& lefttermK,
		  const mwArray& lefttermB){
  bool inmemK=ket->isInMem();
  bool inmemB=ket->isInMem();
  if(!inmemK) ket->read();
  if(!inmemB) bra->read();
  Site::gaugeR(ket,bra,righttermK,righttermB,lefttermK,lefttermB);
  if(!inmemK) ket->write();
  if(!inmemB) bra->write();

    // // Maybe a problem using the same in and output argument=> copy)
    // mwArray auxK(ket->getA());mwArray K;
    // mwArray auxB(bra->getA());mwArray B;
    // gaugeABr(4,K,B,righttermK,righttermB,auxK,auxB,
    // 	     lefttermK,lefttermB);
    // ket->setSite(K);
    // bra->setSite(B);

}
  
void FileSite::increaseBondDim(int Dl_,int Dr_){
  // cout<<"FileSite increaseBondDim()"<<endl;
  bool wasinmem=inmem;
  if(!inmem) read();
  Site::changeDimensions(getd(),Dl_,Dr_);
  dims=Site::getDimensions();
  //cout<<"Before increasing Bond from "<<A.getDimensions()<<" to Dl="<<Dl_
  //  <<", Dr="<<Dr_<<", A="<<A<<endl;
  //Site::increaseBondDim(Dl_,Dr_);
  //cout<<"After changing Dl, Dr, A="<<A<<endl;
  if(!wasinmem) write();
}

void FileSite::decreaseBondDim(int Dl_,int Dr_){
  // cout<<"FileSite decreaseBondDim()"<<endl;
  bool wasinmem=inmem;
  if(!inmem) read();
  Site::changeDimensions(getd(),Dl_,Dr_);
  dims=Site::getDimensions();
  //Site::increaseBondDim(Dl_,Dr_);
  if(!wasinmem) write();
}

void FileSite::increasePhysDim(int d_){
  // cout<<"FileSite increasePhysDim()"<<endl;
  bool wasinmem=inmem;
  if(!inmem) read();
  if(d_<=A.getDimension(0))
    Site::changeDimensions(d_,getDl(),getDr());
  //  Site::increasePhysDim(d_);
  dims=Site::getDimensions();
  if(!wasinmem) write();
}

void FileSite::decreasePhysDim(int d_){
  // cout<<"FileSite decreasePhysDim()"<<endl;
  bool wasinmem=inmem;
  if(!inmem) read();
  Site::changeDimensions(d_,getDl(),getDr());
  //Site::decreasePhysDim(d_);
  dims=Site::getDimensions();
  if(!wasinmem) write();
}


void FileSite::read(){
  ifstream infile(filename,ios::binary);
  if(!infile.is_open()){
    cout<<"Error: unable to open file "<<filename<<" for FileSite"<<endl;
    exit(1);
  }
  else{
    // cout<<"Reading FileSite "<<getPos()<<" from "<<filename<<endl;
    A.load(infile);
    infile.close();
    dims=A.getDimensions();
    inmem=1;  
  }
}


void FileSite::write(){
  if(!inmem){
    cout<<"Error trying to save to file FileSite not in memory!"<<endl;
    exit(1);
  }
  ofstream outfile(filename,ios::binary);
  if(!outfile.is_open()){
    cout<<"Error: unable to open file to write FileSite "<<filename<<endl;
    exit(1);
  }
  else{
    // cout<<"Writing FileSite "<<getPos()<<" to "<<filename<<endl;
    A.save(outfile);
    outfile.close();
  }
  discard();
}

void FileSite::discard(){
  if(inmem){
    A.clear();
    inmem=0;
    // cout<<"FileSite("<<getPos()<<") discarded"<<endl;
  }
}
