/**
   \file MPS.cpp
   Basic class containing an MPS 
   
   \author Mari-Carmen Banuls
   \date 16/03/2011
*/
#include "MPS.h"
#include <iomanip>
using namespace shrt;
using namespace std;

MPS::MPS():length(0),bond(0),A(),normfact(0),gauge(gNone){};

MPS::MPS(int len,int bon,int d):length(len),bond(bon),A(0),
				normfact(1),gauge(gNone){
  //cout<<"MPS constructor with len "<<length<<", bond "<<bond<<<", d "
  //  <<d<<endl;
  A=new Site_t*[length];
  int dimensions[3]={d,1,bond};
  if(length==1) dimensions[2]=1;
  A[0]=new Site_t(1,dimensions);
  if(length>1) dimensions[1]=bond;
  for(int k=1;k<length-1;k++){
    dimensions[1]=bond;
    A[k]=new Site_t(k+1,dimensions);
  }
  dimensions[2]=1;
  A[length-1]=new Site_t(length,dimensions);
  normfact=1.; // default initialization!
  gauge=gNone; //gBoth; no gauge, because dims are not min, no full rank
  //cout<<"Created MPS "<<this<<endl;
}

MPS::MPS(int len,int bon[],int d[]):length(len),bond(1),
				    normfact(1.),A(0),gauge(gNone){
  A=new Site_t*[length];
  int Dl=1;bond=Dl; // the maximum
  for(int k=0;k<length;k++){    
    int Dr=k<length-1?bon[k]:1; 
    //cout<<"Creating for site "<<k+1<<" A("<<d[k]<<","<<Dl<<","<<Dr<<")"<<endl;
    A[k]=new Site_t(k+1,d[k],Dl,Dr);
    Dl=Dr;
    if(Dr>bond) bond=Dr;
  }
  normfact=1.; // default initialization!
  gauge=gNone;
}

MPS::MPS(int length_,int bond_,int* d):length(length_),bond(bond_),
				       normfact(1.),A(0),gauge(gNone){
  // As before, but all the same bond
  A=new Site_t*[length];
  int Dl=1;
  for(int k=0;k<length;k++){    
    int Dr=k<length-1?bond:1; 
    //cout<<"Creating for site "<<k+1<<" A("<<d[k]<<","<<Dl<<","<<Dr<<")"<<endl;
    A[k]=new Site_t(k+1,d[k],Dl,Dr);
    Dl=Dr;
  }
  normfact=1; // default initialization!
  gauge=gNone;
  //cout<<"Created MPS "<<this<<endl;
}

MPS::MPS(int length_,int bond_,const vector<int>& d):
  length(length_),bond(bond_),normfact(1.),A(0),gauge(gNone){
  if(d.size()!=length){
    cout<<"Error: length ("<<length<<") and number of components of dims "
	<<d<<" do not agree to create MPS"<<endl;
    exit(1);
  }
  // As before, but all the same bond
  A=new Site_t*[length];
  int Dl=1;
  for(int k=0;k<length;k++){    
    int Dr=k<length-1?bond:1; 
    //cout<<"Creating for site "<<k+1<<" A("<<d[k]<<","<<Dl<<","<<Dr<<")"<<endl;
    A[k]=new Site_t(k+1,d[k],Dl,Dr);
    Dl=Dr;
  }
  normfact=1; // default initialization!
  gauge=gNone; //gBoth;
  //cout<<"Created MPS "<<this<<endl;
}

MPS::MPS(int length_,const vector<int>& bond_,const vector<int>& d):
  length(length_),bond(1),normfact(1.),A(0),gauge(gNone){
  if((d.size()!=length)||(bond_.size()!=length-1)){
    cout<<"Error: length ("<<length<<") and number of components of dims "
	<<d<<" or bond dimes "<<bond_<<" do not agree to create MPS"<<endl;
    exit(1);
  } 
  A=new Site_t*[length];
  int Dl=1;bond=Dl; // the maximum
  for(int k=0;k<length;k++){    
    int Dr=k<length-1?bond_[k]:1; 
    //cout<<"Creating for site "<<k+1<<" A("<<d[k]<<","<<Dl<<","<<Dr<<")"<<endl;
    A[k]=new Site_t(k+1,d[k],Dl,Dr);
    Dl=Dr;
    if(Dr>bond) bond=Dr;
  }
  normfact=1.; // default initialization!
  gauge=gNone;
}

MPS::MPS(const MPS& aMPS):
  length(aMPS.getLength()),
  bond(aMPS.getBond()),A(0),normfact(1.),gauge(gNone){
  A=new Site_t*[length];
  for(int k=0;k<length;k++){
    A[k]=new Site_t(aMPS.getA(k));
  }
  normfact=aMPS.getNormFact();
  gauge=aMPS.gauge;
  //  cout<<"Created copy of MPS with D="<<bond<<endl;
//   cout<<" Old: "<<aMPS<<endl;
//   cout<<" New: "<<*this<<endl;
}

void MPS::clear(){
  //cout<<"Clearing MPS "<<this<<endl;
  if(A!=0){
    for(int k=0;k<length;k++)
      delete A[k];
    delete []A;
  }
  A=0;  
  length=0;bond=0;
  gauge=gNone;
}

MPS::~MPS(){
  //cout<<"Destruct MPS "<<this<<endl;
  clear();
}
   
const Site_t& MPS::getA(int k) const{
  if(k<0||k>=length){
    cout<<"Error getA trying to access position "<<k<<", while length="
	<<length<<endl;
    exit(1);
  }
  return *A[k];
}

void MPS::setA(int k,const mwArray& data){
  if(k<0||k>=length){
    cout<<"Error setA trying to access position "<<k<<", while length="
	<<length<<endl;
    exit(1);
  }
  A[k]->setSite(data);
  gauge=gNone;
}

void MPS::setRotatedA(int pos,const mwArray& data,const Indices& newdims,
		      bool conjugate,bool checkdims){
  //cout<<"MPS::setRotatedA("<<pos<<") to "<<newdims<<endl;
  if(pos>=length){
    cout<<"Error setRotatedA trying to access position "<<pos<<", while length="
	<<length<<endl;
    exit(1);
  }
  A[pos]->setRotatedSite(data,newdims,conjugate);
  // Check Dl, Dr!
  if(checkdims)
  if(((pos>0&&A[pos]->getDl()!=A[pos-1]->getDr()) // wrong Dl
      ||(pos==0&&A[pos]->getDl()!=1))|| // for the first site
     ((pos<length-1&&A[pos]->getDr()!=A[pos+1]->getDl()) // wrong Dr
      ||(pos==length-1&&A[pos]->getDr()!=1))){ // for the last site
    cout<<"Incompatible dimensions at MPS:setRotatedSite("<<pos<<"): d="
	<<A[pos]->getd()<<", Dl="<<A[pos]->getDl()<<", Dr="
	<<A[pos]->getDr()<<endl;
    exit(1);
  }
  gauge=gNone;
}


void MPS::replaceSite(int pos,const Site_t& data,bool checkdims){
  if(pos>=length){
    cout<<"Error replaceSite(pos,Site_t) trying to access position "<<pos<<", while length="
	<<length<<endl;
    exit(1);
  }
  if(checkdims){
    // check Dl and Dr dims
    if(A[pos]->getDl()!=data.getDl()){
      cout<<"MPS::replaceSite trying to set incompatible Dl dim="
	  <<data.getDl()<<" at pos "<<pos<<", which had oldDl="<<
	A[pos]->getDl()<<endl;
      exit(1);
    }
    if(A[pos]->getDr()!=data.getDr()){
      cout<<"MPS::replaceSite trying to set incompatible Dr dim="
	  <<data.getDr()<<" at pos "<<pos<<", which had oldDr="<<
	A[pos]->getDr()<<endl;
      exit(1);
    }
  }
  delete A[pos];
  A[pos]=new Site_t(data,pos+1);
  gauge=gNone;
}

void MPS::replaceSite(int pos,const mwArray& data,bool checkdims){
  if(pos<0||pos>=length){
    cout<<"Error replaceSite(pos,data) trying to access position "
	<<pos<<", while length="
	<<length<<endl;
    exit(1);
  }
  if(checkdims){
    // check Dl and Dr dims
    Indices dims=data.getDimensions();
    int rank=data.getRank();
    int Dl_=data.getDimension(1);int Dr_=rank<3?1:data.getDimension(2);
    if(A[pos]->getDl()!=Dl_){
      cout<<"MPS::replaceSite trying to set incompatible Dl dim="
	  <<data.getDimension(1)<<" at pos "<<pos<<", which had oldDl="<<
	A[pos]->getDl()<<endl;
      exit(1);
    }
    if(A[pos]->getDr()!=Dr_){
      cout<<"MPS::replaceSite trying to set incompatible Dr dim="
	  <<Dr_<<" at pos "<<pos<<", which had oldDr="<<
	A[pos]->getDr()<<endl;
      exit(1);
    }
  }
  // todo: instead of this, directly replace inside Site (with setSite)
  delete A[pos];
  A[pos]=new Site_t(pos+1,data);
  gauge=gNone;
}


void MPS::applyLocalOperator(int pos,const mwArray& oper,bool normalize){
  // Check dimensions
  if(pos<0||pos>=length){
    cout<<"Error in MPS::applyLocalOperator trying to apply operator on position "<<pos
	<<", when length is "<<length<<endl;
    exit(1);
  }
  Indices dims=oper.getDimensions();
  if(dims.size()!=2||dims[1]!=A[pos]->getd()){
    cout<<"Error in MPS::applyLocalOperator dimensions of operator "<<oper
	<<" incompatible with local site "<<A[pos]->getA()<<endl;
    exit(1);
  }
  mwArray locA=A[pos]->getA();
  locA.reshape(Indices(dims[1],-1));
  locA.multiplyLeft(oper);
  locA.reshape(Indices(dims[0],A[pos]->getDl(),A[pos]->getDr()));
  if(A[pos]->getd()!=dims[0]){
    //    cout<<"Local operator changing physical dimension of MPS site "<<pos<<" from "
    //	<<A[pos]->getd()<<" to "<<dims[0]<<endl;
    A[pos]->changeDimensions(dims[0],A[pos]->getDl(),A[pos]->getDr());
  }
  A[pos]->setSite(locA);
  if(normalize)
    gaugeCond('R',true);
}

void MPS::setProductState(ProductState state){
 for(int k=0;k<length;k++){
   //cout<<"setProductState("<<k<<")"<<endl;
   if(state==p_special&&(k==0||k==length-1)){
     A[k]->setState((SingleState)p_zero);     
   }
   else
     A[k]->setState((SingleState)state);
 }
 //gauge=gBoth; // check?
 gauge=gNone;
}

void MPS::setRandomState(){
  for(int k=0;k<length;k++){
    A[k]->setState(randmps);
    // Now apply gauge cond to the right
    int oldDr=A[k]->getDr();
    mwArray tmp;
    A[k]->gaugeR(tmp,mwArray(ONE_c));
    if((k<length-1)&&(A[k]->getDr()!=oldDr)){
      A[k+1]->reduceDimR(tmp);
    }
  }
  gauge=gRight;
}

ostream& operator<<(ostream& os, const MPS& mps){
  os<<"MPS{length="<<mps.length<<", bond="<<mps.bond<<"}, gauge="<<mps.gauge;
  os<<endl;
  for(int k=0;k<mps.length;k++){
    os<<"Site("<<k<<") "<<mps.A[k]->getDimensions()<<endl;
    //os<<"pos("<<k<<"): "<<mps.getA(k).getA()<<endl;
  }
  return os;
}

void MPS::gaugeCond(char dir,bool normal,int cutD){
  mwArray tmp(ONE_c); // First term
  int maxbond=1; // actual maximum bond
  switch(dir){
  case 'L':
  case'l':
    for(int pos=length-1;pos>=0;pos--){
      mwArray link(tmp); // copy: not efficient! TODO: Areference??
      // cout<<"About to apply Site("<<pos<<").gaugeL, before A="
      // 	  <<A[pos]->getA()<<" and temporary term "
      // 	  <<link<<endl;
      A[pos]->gaugeL(tmp,link,cutD);
      maxbond=max(maxbond,max(A[pos]->getDl(),A[pos]->getDr()));
    };
    gauge=gLeft;
    break;
  case 'R':
  case 'r':
    for(int pos=0;pos<length;pos++){
      mwArray right(tmp); // copy: not efficient! TODO: Areference??
      //cout<<"About to apply Site("<<pos<<").gaugeR, before A="
      //<<A[pos]->getA()<<" and temporary term "
      // <<right<<endl;
      A[pos]->gaugeR(tmp,right,cutD);
      //cout<<"After Site("<<pos<<").gaugeR, A="<<A[pos]->getA()
      //  <<" and temporary term "<<tmp<<endl;
      maxbond=max(maxbond,max(A[pos]->getDl(),A[pos]->getDr()));
      //cout<<"gaugeCondR@"<<pos<<" Dl="<<A[pos]->getDl()<<", Dr="
      //<<A[pos]->getDr()<<endl;
    };
    gauge=gRight;
    break;
  default:
    cout<<"Error unknown direction "<<dir<<" for gauge condition"<<endl;
    exit(2);
  };
  if(maxbond<bond){
    //cout<<"Reducing bond "<<bond<<"->"<<maxbond<<" in gaugeCond()"<<endl;
    bond=maxbond;
  }
  //cout<<"At the end of gaugeCond("<<dir<<"), the last factor is "<<tmp<<endl;
  // WARNING: If tmp<0, I am introducing a global pi phase !!
//   if((double)tmp<0){
//     cout<<"gaugeCond ends up with negative factor=> keep sign by "
// 	<<"multiplying the last modified tensor "<<lastpos<<endl;
//     A[lastpos]->changeSign();
//   }
  // At the end, the MPS is normalized, and the contraction is 1
  if(normal)
    normfact=1.;
  else{   // gaugeCond without normalization? (keeping the factor)
    if(tmp.getNrComponents()!=1){
      cout<<"Error! After gauge"<<dir<<", final tmp size "
	  <<tmp.getDimensions()<<endl;
      exit(212);
    }
    normfact=normfact*abs(tmp.getElement(0));
  }
}

void MPS::gaugeCond(int pos,char dir,bool next){
  mwArray tmp(ONE_c); // First term
  Indices oldD=A[pos]->getDimensions();
  switch(dir){
  case 'L':
  case'l':
    {
    // if(pos==0)
    // 	cout<<"gaugeL(0) from "<<A[pos]->getA()<<endl;
    A[pos]->gaugeL(tmp,mwArray(ONE_c));
    // If dimensions change, multiply the following term
    if(((A[pos]->getDl()!=oldD[1])||next)&&(pos>0)){   
      A[pos-1]->reduceDimL(tmp);
    }
    if(pos==0&&next){ // For the first site, apply normalization
      if(tmp.getNrComponents()!=1){
	cout<<"Error! After gaugeL on site 0, tmp size "
	    <<tmp.getDimensions()<<endl;
	exit(212);
      }
      normfact=normfact*abs(tmp.getElement(0));
    }
    //if(pos==0)
    // 	cout<<" to "<<A[pos]->getA()<<", with tmp="<<tmp<<endl;
    break;
    }
  case 'R':
  case 'r':
    {
    // if(pos==length-1)
    // 	cout<<"gaugeR(length-1) from "<<A[pos]->getA()<<endl;
    A[pos]->gaugeR(tmp,mwArray(ONE_c));
    if((pos<length-1)&&
       ((A[pos]->getDr()!=oldD[2])||next)){
      A[pos+1]->reduceDimR(tmp);
    }
    if(next&&pos==length-1){ // For the last site, apply normalization
      if(tmp.getNrComponents()!=1){
	cout<<"Error! After gaugeR on last site, tmp size "
	    <<tmp.getDimensions()<<endl;
	exit(212);
      }
      normfact=normfact*abs(tmp.getElement(0));
    }
    // if(pos==length-1)
    // 	cout<<" to "<<A[pos]->getA()<<", with tmp="<<tmp<<endl;
    break;
    }
  default:
    cout<<"Error unknown direction for gauge condition"<<endl;
    exit(2);
  };
  gauge=gNone;
}

MPS& MPS::operator=(const MPS& mps){
  if(this!=&mps){ // no self-assignment
    //cout<<"Assigning MPS "<<&mps<<" to "<<this<<endl;
    if(A!=0){
      for(int k=0;k<length;k++)
	delete A[k];
      delete []A;
    }
    length=mps.getLength();
    bond=mps.getBond();
    A=new Site_t*[length];
    for(int k=0;k<length;k++){
      A[k]=new Site_t(mps.getA(k));
    }
    normfact=mps.getNormFact();
    gauge=mps.getGauge();
  }
  return *this;
}


void MPS::exportMPS(const char* filename) const{
  ofstream outfile(filename,ios::binary);
  if(!outfile.is_open()){
    cout<<"Error: unable to open file "<<filename<<" to write"<<endl;
    exit(2);
  }
  else{
    outfile.write((char*)&length,sizeof(int));
    outfile.write((char*)&bond,sizeof(int));
    outfile.write((char*)&normfact,sizeof(double));
    outfile.write((char*)&gauge,sizeof(Gauge));
    if(!outfile){
      cout<<"Error trying to write MPS to file!"<<endl;
      exit(1);
    }
    for(int k=0;k<length;k++)
      A[k]->save(outfile);
    outfile.close();
  }
}

void MPS::exportForMatlab(const char* filename,int precision) const{
  //ofstream outfile(filename,ios::binary);
  ofstream outfile(filename);
  if(!outfile.is_open()){
    cout<<"Error: unable to open file "<<filename<<" to write"<<endl;
    exit(2);
  }
  else{
    if(precision>0)
      outfile<<setprecision(precision);
    outfile<<"tensors={};clear aux;"<<endl;
    //outfile<<setprecision(10);
    for(int k=0;k<length;k++){
      //Ops[k]->savetext(outfile);
      putForMatlab(outfile,A[k]->getA(),"aux",precision);
      outfile<<"tensors{"<<k+1<<"}=aux;clear aux;"<<endl;
    }
    outfile<<"normFactor="<<normfact<<";"<<endl;
    outfile.close();
  }
}

void MPS::importMPS(const char* filename){
  clear();
  ifstream infile(filename,ios::binary);
  if(!infile.is_open()){
    cout<<"Error: unable to open file "<<filename<<" to read"<<endl;
    exit(2);
  }
  else{
    clear();
    infile.read((char*)&length,sizeof(int));
    infile.read((char*)&bond,sizeof(int));
    infile.read((char*)&normfact,sizeof(double));
    infile.read((char*)&gauge,sizeof(Gauge));
    if(!infile){
      cout<<"Error trying to read MPS from file!"<<endl;
      exit(1);
    }
    A=new Site_t*[length];
    for(int k=0;k<length;k++){
      A[k]=new Site_t(k+1,mwArray(0));
      A[k]->load(infile);
    }
    infile.close();
  }
}


void MPS::stretchMPS(int len,const MPS& mps,int reppos){
  // clear what was here before
  for(int k=0;k<length;k++)
    delete A[k];
  delete []A;
  // length of the original MPS
  int len0=mps.getLength();
  if(reppos==0) reppos=len0/2;
  // assigned length
  length=len;
  int deltaL=length-len0; // difference
  bond=mps.getBond();
  A=new Site_t*[length];
  if(deltaL>=0){
    //cout<<"Stretching MPS("<<len0<<") to length="<<length<<endl;;     
    // copy the mps at the edges and fill in copies of the intermediate sites
    for(int k=0;k<len0;k++){
      int kp=k<len0/2?k:k+deltaL;
      //cout<<"Setting Site_t["<<kp<<"] from old ["<<k<<"]"<<endl;
      A[kp]=new Site_t(kp+1,mps.getA(k).getA());
    }
    for(int k=0;k<deltaL;k++){
      int kp=k+len0/2;
      A[kp]=new Site_t(kp+1,mps.getA(reppos).getA());
      //cout<<"Setting Site_t["<<kp<<"] from old ["<<reppos<<"]"<<endl;
    }
    // TODO! Comprobar que las dims del tensor intermedio cumplen Dl=Dr
    normfact=mps.getNormFact();
  }
  else{
    cout<<"WARNING: Cannot stretch a MPS ("<<len0<<") to a smaller length"
	<<length<<"!=> ignoring intermediate sites"
	<<endl;
    // copy the mps at the edges 
    for(int k=0;k<length;k++){
      int kp=k<length/2?k:k-deltaL;
	A[k]=new Site_t(k+1,mps.getA(kp).getA());
	//cout<<"Setting Site_t["<<kp<<"] from old ["<<k<<"]"<<endl;
    }
    normfact=mps.getNormFact();
  }
  gauge=gNone;
}

void MPS::gaugeCond(MPS& ket,MPS& bra,char dir){
  switch(dir){
  case 'R':
  case 'r':{
    //cout<<"Special gaugecond on ket="<<ket
    //	<<" and bra="<<bra<<endl;
    // start imposing proper gauge condition on both arguments
    //if(!ket.isGaugeR()) ket.gaugeCond('R',1);
    //if(!bra.isGaugeR()) bra.gaugeCond('R',1);
    mwArray tmpK(ONE_c); // First term ket
    mwArray tmpB(ONE_c); // First term bra
    int length=ket.length;
    if(length!=bra.length){
      cout<<"Error: gaugeCond(ket,bra) can act only of MPS of same length"
	  <<endl;
      exit(1);
    }
    for(int pos=0;pos<length;pos++){
      mwArray rightK(tmpK); // copy: not efficient! TODO: Areference??
      mwArray rightB(tmpB); // copy: not efficient! TODO: Areference??
      Site_t::gaugeR(ket.A[pos],bra.A[pos],tmpK,tmpB,rightK,rightB);
    };
    ket.normfact=1.;
    bra.normfact=1.;
    ket.gauge=gNone;
    bra.gauge=gNone;
  }
    break;
  case 'L':
  case 'l':{
    //cout<<"Special gaugecond on ket="<<ket
    //	<<" and bra="<<bra<<endl;
    // start imposing proper gauge condition on both arguments
    //    if(!ket.isGaugeL()) ket.gaugeCond('L',1);
    //if(!bra.isGaugeL()) bra.gaugeCond('L',1);
    mwArray tmpK(ONE_c); // First term ket
    mwArray tmpB(ONE_c); // First term bra
    int length=ket.length;
    if(length!=bra.length){
      cout<<"Error: gaugeCond(ket,bra) can act only of MPS of same length"
	  <<endl;
      exit(1);
    }
    for(int pos=length-1;pos>=0;pos--){
      mwArray rightK(tmpK); // copy: not efficient! TODO: Areference??
      mwArray rightB(tmpB); // copy: not efficient! TODO: Areference??
      Site_t::gaugeL(ket.A[pos],bra.A[pos],tmpK,tmpB,rightK,rightB);
    };
    ket.normfact=1.;
    bra.normfact=1.;
    ket.gauge=gNone;
    bra.gauge=gNone;
  }
    break;
  default:
    cout<<"Error unknown direction for gauge condition"<<endl;
    exit(2);
  };
}

// an auxiliary function for random numbers
class X{
public:
  X(){srand(time(NULL));}
};

double randdouble()
{
  static X x;
  static double scale=RAND_MAX+1.;
  double base=rand()/scale;
  double fine=rand()/scale;
  return base+fine/scale;
}

void MPS::perturb(double probability){
  if(probability>1){
    cout<<"ERROR: probability of perturbation must be <1"<<endl;
    exit(1);
  }
  int count=0;
  for(int k=0;k<length;k++){
    double p=randdouble();   
    if(p<probability){
      //cout<<"Kicking site "<<k+1<<" randomly!"<<endl;
      A[k]->setState(randmps);
      count++;
    }
  }
  cout<<"After MPS::perturb("<<probability<<") "<<count<<" out of "<<length
      <<" sites touched"<<endl;
  gauge=gNone;
}

void MPS::increaseBondDimension(int D){
  increaseBondDimensionWithNoise(D,0.);
}

void MPS::increaseBondDimensionWithNoise(int D,double noise){
  if(D<=bond){
    // Cannot be ignored if the physical dimensions have been altered!
    //cout<<"Ignoring increase bond("<<D<<")"<<endl;
    return;
  }
  //  double noisePS=pow(noise,1./(double)length); // noise persite
  double noisePS=noise/((double)(length*(D-bond))); // noise persite
  // cout<<"MPS::increaseBondDimension("<<bond<<"->"<<D<<") with noise "
  //     <<noisePS<<endl;
  bool cutsth=0; // whether some bond was reduced
  int maxbond=1; // actual maximum
  int cumd=1;  // take care of maximal phys dim!
  int cumdrLast=1; //max dimension on a right bond, before reaching sth larger than D
  int boundR=length-1; // pos where that happens
  bool done=false;
  while(!done&&boundR>0){
    int dim=A[boundR]->getd();    
    if(cumdrLast*dim<D){
      boundR--;cumdrLast*=dim;
    }
    else done=true;
  }
  int cumdr=boundR>0?D:cumdrLast*A[0]->getd(); // to keep comparisons simple 
  for(int k=0;k<length;k++){
    Indices oldD=A[k]->getDimensions();
    int newDr,newDl;
    //    if(oldDl<bond||k==0) newDl=oldDl;
    if(k==0)  //first one unchanged
      newDl=oldD[1];
    else {
      if(cumd<D||cumdr<D)
	newDl=min(cumd,cumdr);
      else newDl=D;
    }
    int posd=A[k]->getd();
    if(cumd<D){ // once reached D, I do not increase it any more
      cumd*=posd;
    }
    if(k==boundR) // exactly on the first site where the right dimension matters
      cumdr=cumdrLast; 
    else if(k>boundR)
      cumdr=cumdr/posd;
    if(k==length-1)  //last one
      newDr=oldD[2];
    else {
      if(cumd<D||cumdr<D)
	newDr=min(cumd,cumdr);
      else newDr=D;
    }
    if(newDl<oldD[1]||newDr<oldD[2]) cutsth=1;
    //  cout<<"@"<<k<<" d="<<posd<<" cumd="<<cumd<<", cumdr="<<cumdr
    //	<<", old ("<<oldD[1]
    //	<<", "<<oldD[2]<<") newDs("<<newDl<<", "<<newDr<<")"<<endl;

    A[k]->increaseBondDimWithNoise(newDl,newDr,noisePS);
    maxbond=max(max(newDl,newDr),maxbond);
  }
  //bond=D;
  bond=maxbond;
  //  cout<<"Now bond dimension="<<bond<<endl;
  // I am not sure if the gauge condition changes 
  // (in principle adding zeros shouldn't but cutting things may do it)
  //  if(cutsth) 
  gauge=gNone;
}


void MPS::adjustPhysDimensions(const vector<int>& d){
  if(d.size()!=length){
    cout<<"Error: canot adjust physical dimensions of this MPS: "<<*this
	<<" with dimension vector "<<d<<endl;
    exit(2);
  }
  bool sth=0; //debug
  for(int k=0;k<length;k++){    
    if(d[k]>A[k]->getd()){
      //cout<<"Increasing site "<<k+1<<" A from "<<A[k]->getd()
      // <<" to "<<d[k]<<endl;    
      //cout<<"A was "<<A[k]->getA()<<endl;
      A[k]->increasePhysDim(d[k]);
      sth=1;
//       cout<<"A is "<<A[k]->getA()<<endl;
    }
    else if(d[k]<A[k]->getd()){
      //cout<<"WARNING! Decreasing site "<<k+1<<" A from "<<A[k]->getd()
      //<<" to "<<d[k]<<endl;
      //       cout<<"A was "<<A[k]->getA()<<endl;
      A[k]->decreasePhysDim(d[k]);
      sth=1;
      //cout<<"A is "<<A[k]->getA()<<endl;

    }
    //else
    //cout<<"No change to A["<<k<<"]="<<A[k]->getA()<<endl;
  }  
  if(sth){
    // cout<<"MPS::adjustPhysDims changed sth!!=> reset bond! to "<<bond<<endl;
    int D=bond;
    bond=1;gauge=gNone;
    increaseBondDimension(D);
  }
}

// New creation modalities, by approximation (cut on bond) to given 
// MPS or Operator on MPS.

MPS::MPS(const MPS& aMPS,int D):
  length(aMPS.getLength()),bond(aMPS.getBond()),A(0),
  normfact(aMPS.getNormFact()),gauge(aMPS.gauge)
{
  // \todo! Repeated code from copy constructor: should integrate
  // will *this=aMPS work? 
  A=new Site_t*[length];
  for(int k=0;k<length;k++){
    A[k]=new Site_t(aMPS.getA(k));
  }
  gaugeCond('L',0,D); // keep norm 
}

// // Only for some checking!
// #ifdef TESTINGMPS
// #include "Contractor.h"
// #endif

#ifndef TESTINGMPS
MPS::MPS(const MPO& ops,const MPS& aMPS,int D,char dir,char gaugeDir):
  length(ops.getLength()),bond(D),A(0),gauge(gNone),normfact(1.){
  //  cout<<"MPS constructing from ops on MPS + cut bond; now this is "<<*this<<endl;
  *this=aMPS;
  if(D!=0)
    truncate(ops,aMPS,D,dir,gaugeDir);
  else
    applyExactly(ops,dir,gaugeDir);
}
#endif
// #ifdef TESTINGMPS
//   cout<<"At the end of constructing a MPs from Op on MPS, we have norm="
//       <<normfact<<endl;
//   Contractor& theContr=Contractor::theContractor();
//   complex_t overl=theContr.contract(aMPS,ops,*this);
//   complex_t thisNorm=theContr.contract(*this,*this);
//   cout<<"And overlap with the original is "<<
//     overl<<
//     " when norm was "<<theContr.contract(aMPS,aMPS)
//       <<" and new is "<<thisNorm<<
//     " ratio= "<< real(overl)/(real(thisNorm))<<endl;
//   //cout<<"and S="<<getEntropy()<<endl;
// #endif



// As the ones above, but for use on "live" MPS

void MPS::approximate(const MPS& aMPS,int D){
  if(this!=&aMPS){ // no self-assignment
    *this=aMPS; // Copy
    gaugeCond('L',false,D); // and cut (keep norm)
  }
  else{
    cout<<"approximate failed: trying to assign to same MPS!"<<endl;
  }
}

void MPS::exportMPStext(const char* filename) const{
  ofstream outfile(filename);
  if(!outfile.is_open()){
    cout<<"Error: unable to open file "<<filename<<" to write"<<endl;
    exit(2);
  }
  else{
    outfile<<setprecision(10)<<length<<" "<<bond<<" "<<normfact<<" ";
    for(int k=0;k<length;k++){
      A[k]->savetext(outfile);
    }
    outfile.close();
  }
}

void MPS::importMPStext(const char* filename){
  clear();
  ifstream infile(filename);
  if(!infile.is_open()){
    cout<<"Error: unable to open file "<<filename<<" to read text"<<endl;
    exit(2);
  }
  else{
    clear();
    infile>>length>>bond>>normfact;
    A=new Site_t*[length];
    for(int k=0;k<length;k++){
      A[k]=new Site_t(k+1,mwArray(Indices(1,1,1)));
      A[k]->loadtext(infile);
    }
    infile.close();
  }
}

#ifndef TESTINGMPS
void MPS::approximate(const MPO& ops,const MPS& aMPS,int D,char dir,char gaugeDir){
  bool dagger=dir=='U'?true:false;
  truncate(ops,aMPS,D,dagger,gaugeDir);
}

void MPS::applyExactly(const MPO& ops,bool dagger,
		       char gaugedir){
  mwArray tmp(ONE_c);
  int maxbond=1;
  switch(gaugedir){
  case 'L': // TO the left
  case'l':
    for(int pos=length-1;pos>=0;pos--){
      mwArray rightT(tmp);
      int cutD=A[pos]->getDl()*ops.getOp(pos).getDimensions()[1]; // max dim: no cut
      A[pos]->gaugeLOp(tmp,ops.getOpData(pos),rightT,dagger,cutD);
    }
    gauge=gLeft;
    break;
  case 'R': // TO the right
  case'r':
    for(int pos=0;pos<=length-1;pos++){
      mwArray leftT(tmp);
      int cutD=A[pos]->getDr()*ops.getOp(pos).getDimensions()[3]; // max dim: no cut
      A[pos]->gaugeROp(tmp,ops.getOpData(pos),leftT,dagger,cutD);
    }
    gauge=gRight;
    break;
  default:
    cout<<"Error unknown direction for applyExactly MPO!"<<endl;
    exit(2);
  }
}


void MPS::truncate(const MPO& ops,const MPS& aMPS,int cutD,
		   bool dagger,char dir){
  //cout<<"MPS::truncate (dir"<<dir<<") making a copy of an MPS "<<aMPS<<endl;
  //if(this!=&aMPS) // no self-assignment
  *this=aMPS; // copy
  //Now apply ops exactly and apply gauge before truncating
  mwArray tmp(ONE_c);
  int maxbond=1;
  switch(dir){
  case 'L':
  case'l':
    // Need gauge from right before truncating
    for(int pos=length-1;pos>=0;pos--){
      mwArray link(tmp); // copy: not efficient! TODO: Areference??
      //cout<<"MPS::truncate at pos "<<pos<<" before gaugeLOp, with link "
      //   <<"dimensions "<<link.getDimensions()<<endl;
      A[pos]->gaugeLOp(tmp,ops.getOpData(pos),link,dagger,cutD);
      //cout<<"MPS::truncate at pos "<<pos<<" after gaugeLOp, with tmp "
      //   <<"dimensions "<<tmp.getDimensions()<<endl;
      maxbond=max(maxbond,max(A[pos]->getDl(),A[pos]->getDr()));
    };
    gauge=gLeft;
    break;
  case 'R':
  case 'r':
    for(int pos=0;pos<length;pos++){
      mwArray link(tmp); // copy: not efficient! TODO: Areference??
      // cout<<"MPS::truncate at pos "<<pos<<" before gaugeLOp, with link "
      //   <<"dimensions "<<link.getDimensions()<<endl;
      A[pos]->gaugeROp(tmp,ops.getOpData(pos),link,dagger,cutD);
      // cout<<"MPS::truncate at pos "<<pos<<" after gaugeLOp, with tmp "
      //   <<"dimensions "<<tmp.getDimensions()<<endl;
      maxbond=max(maxbond,max(A[pos]->getDl(),A[pos]->getDr()));
    };
    gauge=gRight;
    break;
  default:
    cout<<"Error unknown direction for gauge condition with MPO!"<<endl;
    exit(2);
  }
  if(tmp.getNrComponents()!=1){
    cout<<"Error! After truncate"<<dir<<" with Ops and MPS, final tmp size "
	<<tmp.getDimensions()<<endl;
    exit(212);
  }  
  if(maxbond<bond) bond=maxbond;
  normfact=normfact*abs(tmp.getElement(0));
}
 
#endif


void expandVec(const MPS& mps,mwArray& vec){
  vec=mwArray(ONE_c);
  int dphys=1; // temporary dimensions of the vector
  int Dv=1;
  int N=mps.getLength();
  vec.reshape(Indices(dphys,Dv));
  for(int k=0;k<N;k++){
    mwArray aux=mps.getA(k).getA();
    //cout<<"Multiplying tensor for site "<<k<<", A="<<aux<<endl;
    Indices dims=aux.getDimensions();
    int d=dims[0];int Dl=dims[1];int Dr=dims[2];
    if(Dl!=Dv){
      cout<<"Error: Dimensions do not agree in expandVec!"<<endl;
      exit(1);
    }
    aux.permute(Indices(2,1,3));
    aux.reshape(Indices(Dl,d*Dr));
    vec.multiplyRight(aux);
    dphys=dphys*d;
    Dv=Dr;
    vec.reshape(Indices(dphys,Dv));
  }
  // if the mps has a normalization factor, we need to include it!!
  vec=mps.getNormFact()*vec;
}  

void vecToMPS(const mwArray& vec,const Indices& dims,MPS& mps,int D){
  Indices locdims=vec.getDimensions();
  mwArray remain(vec);
  if(locdims!=dims){
    remain.reshape(dims);
  }
  int dimL=1;
  Indices dimcur(dims);
  int length=dims.size();
  mps=MPS(length,D,dimcur);
  for(int k=0;k<length-1;k++){
    mwArray U,S,Vdagger;
    remain.reshape(Indices(dimL*dimcur[0],-1));
    wrapper::svd(remain,U,S,Vdagger);
    int dimR=U.getDimension(1); //exact Dright
    if(D<dimR){ // cut to D
      U.resize(Indices(U.getDimensions()[0],D));
      S.resize(Indices(D,S.getDimensions()[1]));
      dimR=D;
    }
    U.reshape(Indices(dimL,dimcur[0],dimR));
    U.permute(Indices(2,1,3));
    mps.replaceSite(k,U,false);
    dimcur.erase(dimcur.begin());
    dimL=dimR;
    remain=S*Vdagger;
    remain.reshape(Indices(dimL,-1));
  }
  remain.permute(Indices(2,1));
  cout<<"In vecToMPS setting pos "<<length-1<<" to dims "<<remain.getDimensions()<<endl;
  mps.replaceSite(length-1,remain,false);
}

void MPS::fold(MPS& mps) const {
  // Site 0 is now the former middle tensor
  //cout<<"Basic MPS::fold(), original length="<<length<<endl;
  bool isodd=length%2!=0;
  // length of the new chain
  //int newlen=length%2==0?length/2:(length+1)/2;
  int newlen=isodd?(length+1)/2:length/2;
  int d=A[0]->getd(); // ignored
  //  cout<<"A[0]="<<*A[0]<<endl;
  if(isodd)
    mps=MPS(newlen,bond,d*d);
  else
    mps=MPS(newlen+1,bond,d*d); // auxiliary site 0 with only contraction  
  //cout<<"new length="<<mps.getLength()<<endl;

  // what is to be substracted from middlept, with k, for folding pair
  int middlept=length/2;
  int offset=isodd?0:1;
  int initpos=0; // Pos in the folded chain
  if(isodd){// the first operator is singular: only the middle one with rearranged indices  
    mwArray aux=A[length/2]->getA();
    Indices dims=A[length/2]->getDimensions();
    aux.reshape(Indices(dims[0],1,dims[1]*dims[2]));
    mps.replaceSite(0,aux,0);
    //cout<<"Odd case: Tensor at 0 set to middle one"<<endl;
    initpos++;
  }
  else{
    // In the even case, the "middle" tensor is substituted by folding two tensors in 
    // the center, but also contracting their common index
    // I will do this with an extra site which contracts both, and has all remaining 
    // dimensions equal to one.
    mwArray aux=identityMatrix(A[length/2]->getDl());
    aux.reshape(Indices(1,1,-1));
    mps.replaceSite(0,aux,0);
    //cout<<"Even case: Tensor 0 set to contraction"<<endl;
  }
  for(int k=initpos;k<newlen;k++){
    int pos1=middlept-k-offset; // left operator
    int pos2=middlept+k; // right operator
    int d=A[pos1]->getd();
    int newpos=isodd?k:k+1;
    // Future?
    // FoldedSite newA(*A[pos1],*A[pos2]);
    //mps.setA(newpos,newA);
    // Meanwhile, simple construction
    mwArray leftA=A[pos1]->getA();
    mwArray rightA=A[pos2]->getA();
    Indices dims1=A[pos1]->getDimensions();
    Indices dims2=A[pos2]->getDimensions();
    leftA.reshape(Indices(-1,1)); // dl Dll Drl
    rightA.reshape(Indices(1,-1)); // dr Dlr Drr
    leftA.multiplyRight(rightA); 
    leftA.reshape(Indices(dims1[0],dims1[1],dims1[2],
			  dims2[0],dims2[1],dims2[2]));
    leftA.permute(Indices(1,4,3,5,2,6));
    leftA.reshape(Indices(dims1[0]*dims2[0],dims1[2]*dims2[1],
			  dims1[1]*dims2[2]));
    mps.replaceSite(newpos,leftA,0);
    //cout<<"Operator "<<k<<" set from "<<pos1<<" and "<<pos2<<endl;
  }
  cout<<"MPS.fold->"<<mps<<endl;
}

void MPS::insertSite(int pos,const mwArray& data,bool checkdims,bool normalize){
  if(pos>length||pos<0){
    cout<<"ERROR: trying to insert Site "<<pos<<" in MPS of length "<<length<<" (allowed 0-"
	<<length<<")"<<endl;
    exit(1);
  }
  // Check proper shape and (if applicalbe) dimensions of the array
  Indices dimsA=data.getDimensions();
  bool okdims=dimsA.size()==3;
  if(pos==0&&dimsA[1]!=1) okdims=0;
  if(pos==length&&dimsA[2]!=1) okdims=0;
  if(checkdims){ // more complete check of dims
    //cout<<"insertSite, checking dimensions of data "<<dimsA<<endl;
    if(pos>0&&dimsA[1]!=A[pos-1]->getDr()) okdims=0;
    if(pos<length&&dimsA[2]!=A[pos]->getDl()) okdims=0;
  }
  if(!okdims){
    cout<<"ERROR: Trying to insert site with array of wrong dimensions: "<<dimsA
	<<" in pos "<<pos<<" of MPS "<<*this<<endl;
    exit(1);
  }
  //cout<<"MPS::insertSite("<<pos<<"), dims "<<dimsA<<endl;
  int newL=length+1;
  Site_t** Along=new Site_t*[newL];
  int k2=0;
  for(int k=0;k<length;k++){
    if(k<pos) {Along[k]=A[k];Along[k]->pos=k;}
    else {Along[k+1]=A[k];Along[k+1]->pos=k+1;}
  }
  Along[pos]=new Site_t(pos,data);
  length=newL;
  delete []A;
  A=Along;
  gauge=gNone;
  if(normalize)
    gaugeCond('L',1);      
}
  
void MPS::traceOutSites(const vector<int>& pos,bool normalize){
  if(pos.size()==0){
    //cout<<"MPS::traceOutSites nothing to do"<<endl;
    return;
  }
  //cout<<"MPS::traceOutSites "<<pos<<endl;
  //  cout<<"Original: "<<*this<<endl;
  int newL=length;
  // First, I collapse the physical dimension of the traced sites:
  for(int k=0;k<pos.size();k++){
    if(pos[k]>=length) continue;
    //cout<<"Index k="<<k<<", pos[k]="<<pos[k]<<endl;
    // Check that the physical dimension is a perfect square
    int d2=A[pos[k]]->getd();
    int d0=sqrt(d2);
    if(d0*d0!=d2){
      cout<<"Error: cannot trace out site "<<pos[k]<<" as it has dimension "<<d2
	  <<" which does not seem to be  of the form d*d"<<endl;
      exit(1);
    }
    // Construct the identity to apply local trace
    mwArray id=identityMatrix(d0);id.reshape(Indices(1,d2));
    applyLocalOperator(pos[k],id);
    //cout<<"After applying the operator, dims of the site are "<<A[pos[k]]->getDimensions()<<endl;
    newL--;
    // If I am not at the end, I now absorb this term to the right
    if(pos[k]<length-1){
      //cout<<"Absorbing the resulting tensor for "<<k<<" in the next one"<<endl;
      mwArray leftT=A[pos[k]]->getA();Indices dims=A[pos[k]]->getDimensions();
      leftT.reshape(Indices(dims[1],dims[2]));
      //cout<<"site at "<<pos[k]<<": leftT="<<leftT<<endl;
      A[pos[k]+1]->reduceDimR(leftT);
      //cout<<"Now site at "<<pos[k]+1<<": "<<A[pos[k]+1]->getDimensions()<<endl;
      mwArray newA=identityMatrix(dims[1]);newA.reshape(Indices(1,dims[1],dims[1]));
      if(dims[2]!=dims[1]){
	A[pos[k]]->changeDimensions(1,dims[1],dims[1]);
	//cout<<"Changed dimensions of site "<<pos[k]<<": "<<A[pos[k]]->getDimensions()<<endl;    
      }
      A[pos[k]]->setSite(newA);
      //      cout<<"Now site at "<<pos[k]<<": "<<A[pos[k]]->getA()<<endl;    
    }
  }

  //  cout<<"After collapsing the physical dimensions: "<<*this<<endl;

  // Now I want to keep an MPS with the remaining number of sites. I
  // use the following trick: first I sweep over the chain from the
  // right, to absorb also possible terms to the right of the
  // rightmost physical site; then I can simply remove the traced out sites

  // First sweep R->L
  int k=length-1;bool traced=true;
  while(k>=0&&traced){
    traced=A[k]->getd()==1;
    if(traced){
      //cout<<"Absorbing the resulting tensor for "<<k<<" in the previous one"<<endl;
      mwArray rightT=A[k]->getA();
      //cout<<"Now dimensions of the site "<<pos[k]<<" are "<<A[k]->getDimensions()<<endl;
      Indices dims=A[k]->getDimensions();
      rightT.reshape(Indices(dims[1],dims[2]));
      //cout<<"site at "<<k<<": rightT="<<rightT<<endl;
      if(k>0){
	A[k-1]->reduceDimL(rightT);
	//cout<<"Now site at "<<k-1<<": "<<A[k-1]->getDimensions()<<endl;
      }
      else{ // special case: k=0, then everything was traced=> the resulting constant should end up in normFactor!
	if(dims[1]*dims[2]!=1){
	  cout<<"Apparently there is an error in traceOutsites, when I am tracing out "<<pos
	      <<" in a MPS of length "<<length<<" (so newL="<<newL<<"), when going to the right after tracing "
	      <<*this<<" and on site 0 cannot continue"<<endl;
	  exit(1);
	} 
	normfact=normfact*abs(A[k]->getA().getElement(Indices(0,0,0)));
      }
      mwArray newA=identityMatrix(dims[2]);newA.reshape(Indices(1,dims[2],dims[2]));
      // I might need to change the dimension by hand to be able to do this
      if(dims[2]!=dims[1]){
	A[k]->changeDimensions(1,dims[2],dims[2]);
	//cout<<"Changed dimensions of site "<<pos[k]<<": "<<A[k]->getDimensions()<<endl;    
      }
      A[k]->setSite(newA);
    }
    k--;
  }
  // And now all the traced sites are just identities, so I can simply remove them.
  Site_t** Ashort=new Site_t*[newL];
  int k2=0;
  for(int k=0;k<length;k++){
    traced=A[k]->getd()==1;
    if(traced) delete A[k];
    else Ashort[k2++]=A[k];
  }
  length=newL;
  delete []A;
  A=Ashort;
  for(int k=0;k<length;k++)A[k]->pos=k; // reset indices
  gauge=gNone;
  if(normalize)
    gaugeCond('R',1);
}


// void MPS::approximateT(const MPS& aMPS,int D,char dir){
//   if(this!=&aMPS){ // no self-assignment
//     *this=aMPS; // Copy
//     gaugeCondT(dir,false,D); // and cut (keep norm)
//   }
//   else{
//     cout<<"approximateT failed: trying to assign to same MPS!"<<endl;
//   }
// }

void MPS::gaugeCondT(char dir,int cutD){
  mwArray tmp(ONE_c); // First term
  int maxbond=1; // actual maximum bond
  switch(dir){
  case 'L': // towards L
  case'l':
    // first need gauge to R imposed
    //gaugeCond('R',false,0); // no cut
    {
    mwArray lambdaR(ONE_c); // A*A^+ product from right cut
    for(int pos=length-1;pos>=0;pos--){
      mwArray auxL(lambdaR);
      mwArray link(tmp); // copy: not efficient! TODO: Areference??
      // cout<<"About to apply Site("<<pos<<").gaugeL, before A="
      // 	  <<A[pos]->getA()<<" and temporary term "
      // 	  <<link<<endl;
      A[pos]->gaugeLT(tmp,link,lambdaR,cutD);
      //mwArray Apos(A[pos]->getA());Indices dims=Apos.getDimensions(); // d Dl Dr
      //auxL.multiplyLeft(reshape(Apos,Indices(dims[0]*dims[1],dims[2])));
      //auxL.reshape(dims);auxL.permute(Indices(2,1,3));auxL.reshape(Indices(dims[1],dims[0]*dims[2])); //Dl d Dr
      //auxL.multiplyRight(reshape(permute(Apos,Indices(1,3,2)),Indices(dims[0]*dims[2],dims[1])));
      //cout<<"After gaugeLT("<<pos //<<"), lambdaR is "<<lambdaR<<endl;
	  //  <<"), A L At="<<auxL<<endl;
      maxbond=max(maxbond,max(A[pos]->getDl(),A[pos]->getDr()));
    };}
    gauge=gNone;
    break;
  case 'R':
  case 'r':
    // This is not clear if it should work
    // first need gauge to L imposed
    //gaugeCond('L',false,0); // no cut
    {
    mwArray lambdaL(ONE_c); // A^+*A product
    for(int pos=0;pos<length;pos++){
      mwArray right(tmp); // copy: not efficient! TODO: Areference??
      //cout<<"About to apply Site("<<pos<<").gaugeR, before A="
      //<<A[pos]->getA()<<" and temporary term "
      // <<right<<endl;
      A[pos]->gaugeRT(tmp,right,lambdaL,cutD);
      //cout<<"After Site("<<pos<<").gaugeR, A="<<A[pos]->getA()
      //  <<" and temporary term "<<tmp<<endl;
      maxbond=max(maxbond,max(A[pos]->getDl(),A[pos]->getDr()));
      //cout<<"gaugeCondR@"<<pos<<" Dl="<<A[pos]->getDl()<<", Dr="
      //<<A[pos]->getDr()<<endl;
    };}
    gauge=gNone;
    break;
  default:
    cout<<"Error unknown direction "<<dir<<" for gauge condition"<<endl;
    exit(2);
  };
  if(maxbond<bond){
    //cout<<"Reducing bond "<<bond<<"->"<<maxbond<<" in gaugeCond()"<<endl;
    bond=maxbond;
  }
  //cout<<"At the end of gaugeCond("<<dir<<"), the last factor is "<<tmp<<endl;
  // WARNING: If tmp<0, I am introducing a global pi phase !!
  //   if((double)tmp<0){
  //     cout<<"gaugeCond ends up with negative factor=> keep sign by "
  // 	<<"multiplying the last modified tensor "<<lastpos<<endl;
  //     A[lastpos]->changeSign();
  //   }
  // At the end, the left term is kept in the norm factor
  if(tmp.getNrComponents()!=1){
    cout<<"Error! After gauge"<<dir<<", final tmp size "
	<<tmp.getDimensions()<<endl;
    exit(212);
  }
  normfact=normfact*abs(tmp.getElement(0));
}


#ifndef TESTINGMPS

#include "DoubleOperator.h"

void MPSfromMPO(const MPO& mpo,MPS& mps,bool up){
  int L=mpo.getLength();
  // vector<int> du=mpo.getDimensions();
  // vector<int> dd=mpo.getOriginalDimensions();
  // int D=mpo.getOp(0).getDr();

  // vector<int> newDims(du);
  // for(int k=0;k<L;k++) newDims[k]*=dd[k];

  mps=MPS(L,1,1);

  for(int k=0;k<L;k++){
    mwArray aux=mpo.getOp(k).getFullData();
    Indices dims=aux.getDimensions();
    if(up)
      aux.permute(Indices(1,3,2,4));
    else
      aux.permute(Indices(3,1,2,4));
    aux.reshape(Indices(dims[0]*dims[2],dims[1],dims[3]));
    mps.replaceSite(k,aux,0); // brute force
  }

}

void MPOfromMPS(const MPS& mps,MPO& mpo,bool up,bool conj){
  int L=mps.getLength();

  mpo.initLength(L);

  for(int k=0;k<L;k++){
    mwArray aux=mps.getA(k).getA();
    Indices dims=aux.getDimensions();
    int d0=sqrt(dims[0]);
    if(d0*d0!=dims[0]){
      cout<<"Error in MPOfromMPS: Dimension of site "<<k<<" ("<<dims[0]
	  <<") does not seem a square=> don't know how t divide it"<<endl;
      exit(1);
    }
    aux.reshape(Indices(d0,d0,dims[1],dims[2]));
    if(up)
      aux.permute(Indices(1,3,2,4),conj);
    else
      aux.permute(Indices(2,3,1,4),conj);
    // take into account the normalization of the MPS, by including it
    // in the first tensor (it could be distributed among all of them)
    if(k==0) aux=mps.getNormFact()*aux;
    mpo.setOp(k,new Operator(aux),true);
  }

}

void doubleMPO(const MPO& simpleMPO,MPO& doubleMPO,bool conj_){
  int L=simpleMPO.getLength();
  doubleMPO.initLength(L); // just in case
  for(int k=0;k<L;k++){
    mwArray aux=simpleMPO.getOp(k).getFullData();
    if(!conj_)
      doubleMPO.setOp(k,new DoubleOperator(aux,permute(aux,Indices(3,2,1,4))),true);
    else
      doubleMPO.setOp(k,new DoubleOperator(aux,conjugate(aux)),true);
  }
}

void extendMPO(const MPO& simpleMPO,MPO& doubleMPO,int dimA){
  int L=simpleMPO.getLength();
  doubleMPO.initLength(L); // just in case
  mwArray identPhys=identityMatrix(dimA);
  identPhys.reshape(Indices(dimA,1,dimA,1)); 
  for(int k=0;k<L;k++){  
    doubleMPO.setOp(k,new DoubleOperator(simpleMPO.getOp(k).getFullData(),identPhys),true);
  }
}

void extendMPO(const MPO& simpleMPO,MPO& doubleMPO){
  int L=simpleMPO.getLength();
  doubleMPO.initLength(L); // just in case
  const vector<int> dims=simpleMPO.getDimensions();
  for(int k=0;k<L;k++){
    int dimA=dims[k];
    mwArray identPhys=identityMatrix(dimA);
    identPhys.reshape(Indices(dimA,1,dimA,1)); 
    doubleMPO.setOp(k,new DoubleOperator(simpleMPO.getOp(k).getFullData(),identPhys),true);
  }
}

void extendTransposeMPO(const MPO& simpleMPO,MPO& doubleMPO,int dimA){
  int L=simpleMPO.getLength();
  doubleMPO.initLength(L); // just in case
  mwArray identPhys=identityMatrix(dimA);
  identPhys.reshape(Indices(dimA,1,dimA,1)); 
  for(int k=0;k<L;k++){  
    doubleMPO.setOp(k,new DoubleOperator(identPhys,permute(simpleMPO.getOp(k).getFullData(),Indices(3,2,1,4))),true);
  }
}

void extendTransposeMPO(const MPO& simpleMPO,MPO& doubleMPO){
  int L=simpleMPO.getLength();
  doubleMPO.initLength(L); // just in case
  const vector<int> dims=simpleMPO.getDimensions();
  for(int k=0;k<L;k++){
    int dimA=dims[k];
    mwArray identPhys=identityMatrix(dimA);
    identPhys.reshape(Indices(dimA,1,dimA,1)); 
    doubleMPO.setOp(k,new DoubleOperator(identPhys,permute(simpleMPO.getOp(k).getFullData(),Indices(3,2,1,4))),true);
  }
}

void diagonalMPOfromMPS(const MPS& mps,MPO& mpo,bool conj){
  int L=mps.getLength();
  mwArray auxId; // for the 3-legged delta
  int d0;
  mpo.initLength(L);

  for(int k=0;k<L;k++){
    mwArray aux=mps.getA(k).getA();
    Indices dims=aux.getDimensions(); // dims[0] is phys. dim
    if(d0!=dims[0]){
      d0=dims[0];
      auxId=mwArray(Indices(d0,d0,d0));
      for(int p=0;p<d0;p++) // the 3-legged delta
	auxId.setElement(ONE_c,Indices(p,p,p));
      auxId.reshape(Indices(d0*d0,d0));
    }
    aux.reshape(Indices(d0,dims[1]*dims[2]));
    aux.multiplyLeft(auxId);aux.reshape(Indices(d0,d0,dims[1],dims[2]));
    aux.permute(Indices(1,3,2,4));
    if(conj)
      aux.conjugate();
    // take into account the normalization of the MPS, by including it
    // in the first tensor (it could be distributed among all of them)
    if(k==0) aux=mps.getNormFact()*aux;
    mpo.setOp(k,new Operator(aux),true);
  }

}

void MPS::conjugateMPS(){
  for(int k=0;k<length;k++){
    A[k]->conjugateA();
  }
}

#endif
