#include "mwArray.h"
//#ifndef MKL_VER
#include "TransposerMap.h"
//#else
#ifdef MKL_VER
extern "C" {
  // in place transposition of a matrix
  void mkl_zimatcopy (char* ordering,char* trans,int* rows,int* cols,
		      complex_t* alpha,double* a,int* src_lda,int* dst_lda);
  }
#endif

using namespace std;
using namespace shrt;

// static member empty
const mwArray mwArray::empty=mwArray();

const mwArray identityMatrix(int N){
  mwArray A(Indices(N,N));
  for(int k=0;k<N;k++){
    A.setElement(1.,0.,Indices(k,k));
  }
  return A;
}

void mwArray::getDiagonal(mwArray& diag) const{
  //cout<<"mwArray::getDiagonal() of "<<*this<<endl;
  if(isMatrix()){
    //cout<<"mwArray::getDiagonal() is a matrix "<<endl;
    int minDim=min(dims[0],dims[1]);
    diag=mwArray(Indices(1,minDim));
    //cout<<"Constructed place for the diagonal "<<diag<<endl;
    for(int k=0;k<minDim;k++){
      //cout<<"Element "<<k<<" of the diagoonal should be: "<<getElement(Indices(k,k))<<endl;
      diag.setElement(getElement(Indices(k,k)),Indices(0,k));
    }
    //cout<<"At the end of getDiagonal, diag="<<diag<<endl;
  }
  else
    if(isScalar())
      diag=*this; // just a copy
    else
      if(isVector())
	diag=mwArray(getElement(0)); // just the first element
      else{
	cout<<"Error: cannot extract the diagonal of a mwArray "<<*this<<endl;
	exit(1);
      }
}

void mwArray::getRealDiagonal(std::vector<double>& diag) const{
  mwArray aux;getDiagonal(aux);
  diag.clear();  
  for(int k=0;k<aux.getDimension(1);k++)
    diag.push_back(real(aux.getElement(k)));
}

void mwArray::getDiagonal(std::vector<complex_t>& diag) const{
  diag.clear();
  if(isMatrix()){
    //cout<<"mwArray::getDiagonal() is a matrix "<<endl;
    int minDim=min(dims[0],dims[1]);

    for(int k=0;k<minDim;k++)
      diag.push_back(getElement(Indices(k,k)));
  }
  else{
    if(isScalar()||isVector())
      diag.push_back(getElement(0)); // just the first element
    else{
      cout<<"Error: cannot extract the diagonal of a mwArray "<<*this<<endl;
      exit(1);
    }
  }
}




const mwArray realDiag(int nr,const double* values){
  Indices size(nr,nr);
  mwArray A(size);
  for(int k=0;k<nr;k++){
    A.setElement(values[k],0.,Indices(k,k));
  }
  return A;
}

const mwArray realDiag(int nr,const vector<double>& values){
  Indices size(nr,nr);
  mwArray A(size);
  for(int k=0;k<nr;k++){
    A.setElement(values[k],0.,Indices(k,k));
  }
  return A;
}

const mwArray diag(int nr,const double* values){
  Indices size(nr,nr);
  mwArray A(size);
  for(int k=0;k<nr;k++){
    A.setElement(values[2*k],values[2*k+1],Indices(k,k));
  }
  return A;
}

const mwArray diag(int nr,const complex_t* values){
  Indices size(nr,nr);
  mwArray A(size);
  for(int k=0;k<nr;k++){
    A.setElement(values[k],Indices(k,k));
  }
  return A;
}

const mwArray diag(const vector<complex_t>& values){
  int nr=values.size();
  Indices size(nr,nr);
  mwArray A(size);
  for(int k=0;k<nr;k++){
    A.setElement(values[k],Indices(k,k));
  }
  return A;
}

const mwArray diag(const mwArray& values){
  if(values.isVector()){
    int nr=values.getNrComponents();
    Indices size(nr,nr);
    mwArray A(size);
    for(int k=0;k<nr;k++){
      A.setElement(values.getElement(k),Indices(k,k));
    }
    return A;
  }
  else{
    cout<<"ERROR: cannot construct a diagonal matrix from the non-vector tensor "<<values.getDimensions()<<endl;
    exit(1);
  }
}

const mwArray mwArray::subArray(const Indices& ind) const{
  if(dims.size()==1){ // special case: only one dimension (scalar or vector)
    // TODO: Check that any extra index in ind is -1. 
    cout<<"WARNING: subArray used in array of single dimension, returning element"<<endl;
    mwArray B(getElement(ind[0]));
    return B;
  }
  if(ind.size()!=dims.size()){
    cout<<"Error: cannot extract subArray with indices "<<ind<<
      "from mwArray"<<*this<<endl;
    exit(1);
  }
  Indices newdims;
  Indices oldInd(ind);
  // TODO: reverse: more efficient?
  int* ptr=new int[ind.size()]; // place to map new dims to old ones
  int cnt=0;
  for(int k=0;k<dims.size();k++){
    if(ind[k]<0){ 
      newdims.push_back(dims[k]);
      ptr[cnt++]=k;
    }
  }
  //cout<<"SubArray "<<ind<<", newdims="<<newdims<<endl;
  mwArray B(newdims);
  Indices newind;
  for(int k=0;k<B.getNrComponents();k++){
    B.getIndices(k,newind);
    for(int l=0;l<newdims.size();l++){
      oldInd[ptr[l]]=newind[l];
    }
    B.setElement(getElement(oldInd),k);
  }
  delete []ptr;
  return B;
}

mwArray::mwArray():rank(0),dims(),cumdims(),nrcomponents(0),computed(false),
		   mine(true),components(NULL){};

mwArray::mwArray(complex_t c):dims(Indices(1)),rank(1),nrcomponents(1),
			      cumdims(Indices(1)),mine(true),
			      computed(false),components(NULL){
  components=new double[2];
  *(complex_t*)components=c;
  computed=true;mine=true;
}

mwArray::mwArray(const Indices& dimensions):
  dims(dimensions),rank(dimensions.size()),
  nrcomponents(1),computed(false),mine(true),
  cumdims(),components(NULL){
  for(int k=0;k<rank;k++) nrcomponents=nrcomponents*dims[k];
  components=new double[nrcomponents*2]; // complex
  fillWithZero();
} 

mwArray::mwArray(const vector<double>& revals,bool column):
  nrcomponents(revals.size()),components(NULL),cumdims(),dims(),
  mine(true),computed(false),rank(2){
  dims=column?Indices(nrcomponents,1):Indices(1,nrcomponents);
  components=new double[nrcomponents*2]; // complex
  for(int k=0;k<nrcomponents;k++){
    components[2*k]=revals[k]; 
    components[2*k+1]=0;
  }
}

mwArray::mwArray(const vector<complex_t>& cvals,bool column):
  nrcomponents(cvals.size()),dims(),cumdims(),
  components(NULL),mine(true),rank(2),computed(false){
  dims=column?Indices(nrcomponents,1):Indices(1,nrcomponents);
  components=new double[nrcomponents*2]; // complex
  for(int k=0;k<nrcomponents;k++){
    *((complex_t*)&components[2*k])=cvals[k];
  }
}


void mwArray::cumDims(){
  cumdims=Indices(vector<int>(rank,1));
  int tot=1;
  for(int k=1;k<rank;k++){
    tot=tot*dims[k-1];
    cumdims[k]=tot;
  }
  computed=true;
}

mwArray::mwArray(int rank_,const int* dimensions):
  rank(rank_),dims(vector<int>(rank,0)),nrcomponents(1),
  components(NULL),mine(true),cumdims(),computed(false){
  for(int k=0;k<rank;k++){ 
    dims[k]=dimensions[k];
    nrcomponents=nrcomponents*dims[k];
  }
  components=new double[nrcomponents*2]; // complex
  fillWithZero();
  computed=false;mine=true;
}

/**
mwArray::mwArray(int rank_,...):rank(rank_),dims(rank_,1),nrcomponents(1){
  va_list dimensions;
  va_start(dimensions,rank_);
  for(int k=0;k<rank;k++){
    dims[k]=va_arg(dimensions,int);
    cout<<"Read dim["<<k<<"]="<<dims[k]<<endl;
    nrcomponents=nrcomponents*dims[k];
g  }
  cout<<"Nr components="<<nrcomponents<<", dims="<<dims<<", rank="<<rank<<endl;
  components=new double[nrcomponents*2]; // complex
  fillWithZero();
  //  cout<<"Constructed from variable list of input args "<<*this<<endl;
  }*/

mwArray::mwArray(const Indices& dimensions,const double* data,bool real):
  dims(dimensions),rank(dimensions.size()),nrcomponents(1),computed(false),
  mine(true),components(NULL),cumdims(){
  for(int k=0;k<rank;k++) nrcomponents=nrcomponents*dims[k];
  components=new double[nrcomponents*2]; // complex
  if(real)
    for(int k=0;k<nrcomponents;k++){      
      components[2*k]=data[k];
      components[2*k+1]=0.0;
    }
  else
    for(int k=0;k<2*nrcomponents;k++){
      components[k]=data[k];
    }
  mine=true;
}

mwArray::mwArray(const Indices& dimensions,const complex_t* data):
  dims(dimensions),rank(dimensions.size()),nrcomponents(1),computed(false),
  components(NULL),mine(true),cumdims(){
  for(int k=0;k<rank;k++) nrcomponents=nrcomponents*dims[k];
  components=new double[nrcomponents*2]; // complex
  for(int k=0;k<nrcomponents;k++){
    //components[2*k]=data[k].re; 
    //components[2*k+1]=data[k].im;
    // todo: can use built-in-like??
    *((complex_t*)&components[2*k])=data[k];
  }
  mine=true;
}


mwArray::mwArray(ifstream& data):
  rank(0),computed(false),mine(true),nrcomponents(0),
  components(NULL),dims(),cumdims(){
  this->load(data);
}

mwArray::mwArray(const mwArray& orig):
  dims(orig.dims),rank(orig.rank),nrcomponents(orig.nrcomponents),
  computed(false),components(NULL),mine(true),cumdims(){
  if(orig.nrcomponents!=0)
    components=orig.copyComponents();
  // CHECK!
}

void mwArray::setRealData(int nr,...){
  if(nr!=nrcomponents){
    cout<<"Error: trying to set real components for "<<nr
	<<" elements, when size is "<<nrcomponents<<endl;
    exit(212);
  }
  va_list values;
  va_start(values, nr);
  for(int k=0;k<nr;k++){
    components[2*k]=va_arg(values,double);
  }
  va_end(values);
}

void mwArray::setImagData(int nr,...){
  if(nr!=nrcomponents){
    cout<<"Error: trying to set imaginary components for "<<nr
	<<" elements, when size is "<<nrcomponents<<endl;
    exit(212);
  }
  va_list values;
  va_start(values, nr);
  for(int k=0;k<nr;k++){
    components[2*k+1]=va_arg(values,double);
  }
  va_end(values);
}

void mwArray::setData(int nr,const complex_t data[]){
  if(nr!=nrcomponents){
    cout<<"Error: trying to set complex components for "<<nr
	<<" elements, when size is "<<nrcomponents<<endl;
    exit(212);
  }
  for(int k=0;k<nrcomponents;k++){
    *(complex_t*)&components[2*k]=data[k];
  }
}

void mwArray::fillWithZero(){
  for(int k=0;k<nrcomponents;k++){
    components[2*k]=0.;
    components[2*k+1]=0.;
  }
}

void mwArray::fillWithOne(){
  for(int k=0;k<nrcomponents;k++){
    components[2*k]=1.;
    components[2*k+1]=0.;
  }
}

mwArray::~mwArray(){
  clear();
}

void mwArray::clear(){
  rank=0;
  dims.clear();
  if(mine)
    if(components!=0)
      delete[] components;
  components=0;
  cumdims.clear();computed=false;
  mine=true;
}

double* mwArray::copyComponents() const{
  double* _comp=new double[2*nrcomponents];
  // CHECK!
  //  memcpy((void*)_comp,(void*)components,sizeof(double)*2*nrcomponents); 
  memcpy((void*)_comp,(void*)components,sizeof(complex_t)*nrcomponents); 
  return _comp;
}


void mwArray::load(ifstream& data){
  // Discard previous contents
  dims.clear();
  cumdims.clear();computed=false;
  data.read((char*)&rank,sizeof(int));
  if(!data){ // sth happened when reading
    cout<<"Error trying to read mwArray rank from file "<<endl;
    exit(1);
  } 
  //cout<<"load read rank="<<rank<<endl;
  int* tmpDim=new int[rank];
  data.read((char*)tmpDim,sizeof(int)*rank);
  if(!data){ // sth happened when reading
    cout<<"Error trying to read mwArray dimensions from file ("
	<<sizeof(int)*rank<<" characters needed, only "
	<<data.gcount()<<" could be extracted)"<<endl;
    exit(1);
  } 
  //cout<<"load read dimensions (";
  //for(int k=0;k<rank;k++) cout<<tmpDim[k]<<",";
  //cout<<")"<<endl;
  //exit(1);
  int totalsize=1;
  for(int k=0;k<rank;k++){ 
    dims.push_back(tmpDim[k]);
    totalsize=totalsize*tmpDim[k];
  }
  delete[] tmpDim;
  if(!mine){
    if(totalsize!=nrcomponents){
      cout<<"Error: shallow mwArray can only load data of the same size"<<endl;
      exit(212);
    }
  }
  else{
    if(components!=0) delete[] components;
    components=new double[2*totalsize];
  }
  nrcomponents=totalsize;
  data.read((char*)components,sizeof(double)*2*totalsize); 
  if(!data){ // sth happened when reading
    cout<<"Error trying to read mwArray from file ("
	<<sizeof(double)*2*totalsize<<" elements needed, only "
	<<data.gcount()<<" could be extracted)"<<endl;
    exit(1);
  } 
}

void mwArray::save(ofstream& data) const {
  data.write((char*)&rank,sizeof(int));
  int* tmpDim=new int[rank];
  for(int k=0;k<rank;k++){
    tmpDim[k]=dims[k];
  }
  data.write((char*)tmpDim,sizeof(int)*rank);
  data.write((char*)components,sizeof(double)*2*nrcomponents);
  delete[] tmpDim;
  if(!data){
    cout<<"Error writing mwArray to file"<<endl;
    exit(1);
  }
}

bool mwArray::checkIndices(const Indices& indices) const{
  for(int k=rank-1;k>=0;k--)
    if(indices[k]<0||indices[k]>=dims[k])
      return false;
  return true;
}

int mwArray::linearIdx(const Indices& indices) const{
  // Fortran order: column major
  int idx=0;
  for(int k=rank-1;k>=0;k--){
    if(indices[k]<0||indices[k]>=dims[k]){
      cout<<"Error: in linearIdx invalid index["<<k<<"]="<<indices[k]
	  <<" for mwArray "<<*this<<" of corresponding dimension "
	  <<dims[k]<<endl;
      exit(212);
    }
    idx=dims[k]*idx+indices[k];
  }
  //cout<<"linearIdx"<<indices<<"->"<<idx<<endl;
  return idx; // complex numbers
}


int mwArray::_linearidx_(const Indices& indices) const{
  // Fortran order: column major
  int idx=0;
  if(!computed) {
    cout<<"ERROR: linearidx needs cumDim to be called first!"<<endl;
    exit(4);
  }
  for(int k=0;k<rank;k++)
    idx+=indices[k]*cumdims[k];
  //cout<<"NEW linearidx"<<indices<<"->"<<idx;
  //cout<<"CHECKING with old ->"<<linearIdx(indices)<<endl;
  return idx; // complex numbers
}

int mwArray::linearIdx(const Indices& indices,const Indices& dims_) const {
  // Fortran order: column major
  int idx=0;
  for(int k=rank-1;k>=0;k--){
    if(indices[k]<0||indices[k]>=dims_[k]){
      cout<<"Error: in linearIdx(vec,vec) invalid index["<<k<<"]"<<indices[k]
	  <<" for mwArray "<<*this<<" of corresponding dimension "
	  <<dims_[k]<<endl;
      exit(212);
    }
    idx=dims_[k]*idx+indices[k];
  }
  //cout<<"linearIdx"<<indices<<", for dims"<<dims<<"->"<<idx<<endl;
  return idx; // complex numbers
}


int mwArray::_linearidx_(const Indices& indices,const Indices& dims_) const {
  static Indices mydims=dims_;
  static Indices _cumdims(dims_.size(),1);
  static bool cumdimsOK=false;
  if(!cumdimsOK||mydims!=dims_){
    mydims=dims_;_cumdims=Indices(vector<int>(dims_.size(),1));int tot=1;
    for(int k=1;k<dims_.size();k++){
      tot=tot*dims_[k-1];_cumdims[k]=tot;
    }
    cumdimsOK=true;
  }
  // Fortran order: column major
  int idx=0;
  for(int k=0;k<rank;k++)
    idx+=indices[k]*_cumdims[k];
  //cout<<"NEW linearidx"<<indices<<",dims"<<dims_<<"->"<<idx;
  //cout<<"CHECKING with old ->"<<linearIdx(indices,dims_)<<endl;
  return idx; // complex numbers
}

double* mwArray::element(int d1,...){
  Indices pos(rank,d1);
  va_list indices;
  va_start(indices, d1);
  for(int k=1;k<rank;k++){
    pos[k]=va_arg(indices,int);
  }
  va_end(indices);
  //cout<<"Accessing element "<<pos<<endl;
  cumDims();
  int linIdx=_linearidx_(pos);
  return &components[2*linIdx];
}

void mwArray::setElement(double revalue,double imvalue,
			 const Indices& indices){
  //cout<<"Setting component "<<indices<<endl;
  if(!checkIndices(indices)){
    cout<<"Error: Element "<<indices<<" out of range for mwArray of size "<<dims<<endl;
    exit(1);
  }
  cumDims();
  int linidx=_linearidx_(indices);
  components[2*linidx]= revalue;
  components[2*linidx+1]= imvalue;
}

void mwArray::setElement(complex_t value,const Indices& indices){
  setElement(value.re,value.im,indices);
}

void mwArray::setElement(complex_t value,int linidx){
  components[2*linidx]= value.re;
  components[2*linidx+1]= value.im;
}

void mwArray::reshape(const Indices& newdims){
  if(newdims!=dims){ // otherwise ignore
    // First check that the total dimensions are right:
    int newtot=1;int newrank=newdims.size();
    int placeholder=-1; // whether there is a negative index
    for(int k=0;k<newrank;k++){ 
      if(newdims[k]<0){
	if(placeholder>=0){ // there can be only one
	  cout<<"Error: reshape receives wrong dimensions "<<newdims<<endl;
	  exit(212);
	}
	placeholder=k;
      }
      else newtot=newtot*newdims[k];
    }
    dims=newdims;
    // if there was an empty dimension, fill it in now
    if(placeholder>=0) 
      dims[placeholder]=nrcomponents/newtot;
    else if(newtot!=nrcomponents){
      cout<<"Error: reshape cannot change the total number of elements("
	  <<nrcomponents<<" to "<<newtot<<")"<<endl;
      exit(212);
    }
    rank=newrank;
    cumdims.clear();computed=false;
  }
}

bool incrementOne1(Indices& orig,const Indices& dims,int rank){
  bool done=false;
  int pos=0;
  //cout<<orig<<"+1 (DIMS"<<dims<<"!)->";
  while(!done&&pos<rank){
    if(orig[pos]<dims[pos]-1){
      orig[pos]++;done=true;
      for(int k=0;k<pos;k++) orig[k]=0;
    } 
    else{
      pos++;
    }
  }
  //cout<<orig<<endl;
  return done;
}


bool incrementOne(Indices& orig,const Indices& dims,int rank){
  bool done=false;
  int pos=0;
  //cout<<orig<<"+1 (DIMS"<<dims<<"!)->";
  while(!done&&pos<rank){
    if(orig[pos]==dims[pos]-1){
      orig[pos++]=0;
    }
    else{
      orig[pos]++;done=true;
    }
  }
  //cout<<orig<<endl;
  return done;
}

bool change(const Indices& neworder,const Indices& dims){
  int sz=dims.size();
  if(neworder.size()!=sz) return true;
  for(int k=0;k<sz;k++)
    if(neworder[k]!=k+1) return true;
  //cout<<"No change permuting "<<dims<<" as "<<neworder<<endl;
  return false;
}

bool isTransposition(const Indices& ind){
  return ind.size()==2&&ind[0]==2&&ind[1]==1;
}

bool isExchange(const Indices& neworder,int& i1,int& i2){
  int sz=neworder.size();
  int changed[3]={0,0,0};
  int nr=0;
  for(int k=0;k<sz;k++){
    if(neworder[k]!=k+1){
      changed[nr++]=k;
    }
    if(nr>2) return false;
  }
  if(nr==2){i1=changed[0];i2=changed[1]; return true;}
  return false;
}

void mwArray::permute(const Indices& neworder,bool conj){
  if(!change(neworder,dims)) return; // nothing to do
  if(isScalar()) 
    if(!conj) return;
    else{
      conjugate(); return;
    }
  // Particular case: transposition of 2D matrix
  //if(isTransposition(neworder)&&isSquare()) return transpose(conj);
  if(isTransposition(neworder)) return transpose(conj);
  // Particular case: extra indices, understood as 1
  if(neworder.size()>rank){
    Indices newdims(dims);
    while(newdims.size()<neworder.size()){
      newdims.push_back(1);
    }
    //cout<<"In permute, dims "<<dims<<" reshaped to "<<newdims<<endl;
    reshape(newdims);
    // todo: avoid
    permute(neworder);
  }
  if(neworder.size()!=dims.size()){
    cout<<"Error! Cannot permute "<<dims<<" as "<<neworder<<endl;
    exit(2);
  }
  double* newcomp=new double[nrcomponents*2];
  // auxidx holds indices from the second on, and increases one by one
  // (the first one is looped by hand) newdims holds the permuted
  // dimensions (which are the dims of the final vector) auxcum holds
  // the permutation of the cumulative product of dimensions for the
  // final vector, i.e., from the new dims, a permutation of newcumdims
  Indices auxidx(vector<int>(rank-1,0)),auxcum(vector<int>(rank,0));
  Indices newdims(dims); // new order
  Indices invorder(neworder); // reverse permutation
  Indices newcumdims(vector<int>(rank,1)); // new cumulative product of dimensions
  for(int p=0;p<rank;p++){
    int _ne=neworder[p]-1;
    newdims[p]=dims[_ne];
    if(p>0) newcumdims[p]=newcumdims[p-1]*newdims[p-1];
    invorder[_ne]=p+1;
    auxcum[_ne]=newcumdims[p];   //auxcum[p]=newcumdims[invorder[p]-1];
  }     
  //cout<<"Permuting old dims "<<dims<<" to "<<neworder<<", inverse perm="<<invorder<<endl;
  //cout<<"So, newdims="<<newdims<<", and newcumdims="<<newcumdims<<" and auxcum="<<auxcum<<endl;

  Indices auxdims(dims);auxdims.erase(auxdims.begin()); // remove the first one
  // loop one by one over the whole data and copy adequately
  bool done=0;int dims0=dims[0];int auxcum0=auxcum[0];
  // int k=0;
  //double re,im;
  double* ptr=&components[0];
  //  while(k<nrcomponents){
  while(!done){
    //cout<<"Looking at elements with indices "<<auxidx<<endl;
    int offset=0;
    for(int p=0;p<rank-1;p++) offset+=2*auxcum[p+1]*auxidx[p];
    for(int i1=0;i1<dims0;i1++){
      newcomp[offset]=*(ptr++); // real
      //cout<<"offset "<<offset<<endl;
      newcomp[++offset]=conj?-*(ptr++):*(ptr++); // imag  k++
      //cout<<"offset "<<offset<<endl;
      //cout<<"Copied element "<<k++<<" to position "<<(offset-1)/2
      //  <<", value "<<newcomp[offset-1]<<"+i "<<newcomp[offset]<<endl;
      offset+=2*auxcum0-1; // counting from the imaginary part!
    }
    done=!incrementOne(auxidx,auxdims,rank-1);
  }
  ptr=0;
  if(mine){
    delete[] components;
    components=newcomp;
  }
  else{ // copy in place and delete auxiliary
    memcpy((void*)components,(void*)newcomp,sizeof(double)*2*nrcomponents);
    delete[] newcomp;
  }
  dims.swap(newdims);
  computed=false;
}
#ifdef MKL_VER_INEFF
// use the MKL library to do in place transpositions?
// Less efficient: not used now
void mwArray::transpose(bool conj){
  if(isMatrix()){
    int d1=dims[0];int d2=dims[1];
    char ordering='C'; // column major
    char trans=conj?'C':'T'; // transpose-conjugate or transpose alone
    complex_t alpha={1.,0.}; // scaling factor
    mkl_zimatcopy(&ordering,&trans,&d1,&d2,&alpha,components,&d1,&d2);
    dims[0]=d2;dims[1]=d1;
    computed=false;
  }
  else{
    if(isScalar()){
      if(conj) conjugate();
    }
    else if(isVector()){
      if(conj) conjugate();
      return reshape(Indices(dims[1],dims[0]));
    }
    else{
      cout<<"Error: cannot transpose an array of dims "<<dims<<endl;
      exit(212);
    }
  }
}      
#else
void mwArray::transpose(bool conj){
  if(isMatrix()){
    int d1=dims[0];int d2=dims[1];
    if(d1==d2){
      int offset=0;
      double* ptr=&components[0];
      //int k=0;
      for(int i2=0;i2<d2;i2++){
	int offset=2*i2;
	//for(int i1=0;i1<d1;i1++){
	for(int i1=0;i1<i2;i1++){
	  //TODO: instead of an int offset, use the pointer directly?
	  complex_t aux=*(complex_t*)ptr;
	  //cout<<"Exchanging elements "<<aux<<" and " 
	  //  <<*(complex_t*)(&components[offset])<<endl;
	  if(conj){
	    *(complex_t*)ptr=::conjugate(*(complex_t*)(&components[offset]));
	    *(complex_t*)(&components[offset])=::conjugate(aux);
	  }
	  else{
	    *(complex_t*)ptr=*(complex_t*)(&components[offset]);
	    *(complex_t*)(&components[offset])=aux;
	  }
	  //	offset+=2*d2-1;
	  ptr+=2;
	  //cout<<"orig k="<<k<<"->"<<offset<<endl;k+=2;
	  offset+=2*d2;
	}
	if(conj){ // also conjugate the diagonal
	  *(complex_t*)ptr=::conjugate(*(complex_t*)ptr);
	}
	// Now I have to jump over the lower diagonal!
	ptr+=(d1-i2)*2;
      }
    }
    else{ // for non-square matrix I have to follow cycles.
      TransposerMap& theMap=TransposerMap::getTransposerMap();
      // check if the Transposer is known. If not, create and archive
      // it 
      Transposer* trans;
      //cout<<"Looking up map for Transposer of "<<dims<<endl;
      iter_t posIter=theMap.find(dims);
      if(posIter!=theMap.end())
	trans=posIter->second;
	//cout<<"Recovered from map transposer for "<<dims<<" as "<<*trans<<endl;}
      else{
	trans=new Transposer(dims,Indices(2,1));
	theMap.insert(make_pair(dims,trans));
	//cout<<"Inserted transposer for "<<dims<<endl;
      }
      // Use the transposer
      trans->transpose(*this,conj);
    }   
    dims[0]=d2;dims[1]=d1;
    computed=false;
  }
  else{
    if(isScalar()){
      if(conj) conjugate();
    }
    else if(isVector()){
      //cout<<"Transposing a vector (conj="<<conj<<")"<<endl;
      if(conj) conjugate();
      return reshape(Indices(dims[1],dims[0]));
    }
    else{
      cout<<"Error: cannot transpose an array of dims "<<dims<<endl;
      exit(212);
    }
  }
}
#endif   


/*
void mwArray::exchangeTwo(int i1,int i2){
  // Do a permutation of two indices corresponding to the same dimension

}
*/

void mwArray::conjugate(){
  for(int k=0;k<nrcomponents;k++)
    components[2*k+1]=-components[2*k+1];
}

void mwArray::Hconjugate(){
  if(isScalar()&&rank==1){
    conjugate(); 
    return;
  }
  if(rank==2){
    //transpose(true);
    //Indices neworder(2,1);
    permute(Indices(2,1),true);
  }
  else{
    cout<<"Error: cannot apply Hermitian conjugate to array of size"<<dims
	<<endl;
    exit(212);
  }
}


const mwArray reshape(const mwArray& A,const Indices& newdims){
  mwArray B(A);
  B.reshape(newdims);
  return B;
}

const mwArray resize(const mwArray& A,const Indices& newdims){
  mwArray B(A);
  B.resize(newdims);
  return B;
}

const mwArray permute(const mwArray& A,const Indices& neworder){
  // Could be faster if I avoid the copy, but repeat the code somehow!
  mwArray B(A);
  B.permute(neworder,false);
  return B;
}
const mwArray conjugate(const mwArray& A){
  mwArray B(A);
  B.conjugate();
  return B;
}
const mwArray Hconjugate(const mwArray& A){
  mwArray B(A);
  B.Hconjugate();
  return B;
}

double norm(const mwArray& A){
  if(A.isScalar()){
    complex_t singleEl=A.getElement(0);
    return abs(singleEl);
  }
  else if(A.isVector()){
    complex_t aux=scalarproduct(conjugate(A),A);
    return sqrt(abs(aux));
  }
  else{
    cout<<"Error: norm not supported for mwArray of dims "
	<<A.getDimensions()<<endl;
    exit(212);
  }
}

const mwArray sqrt(const mwArray& mwA){
  mwArray B(mwA);
  for(int k=0;k<B.nrcomponents;k++){
    double re=B.components[2*k];
    double im=B.components[2*k+1];
    double z=sqrt(re*re+im*im);
    double phi=im==0?0:(re==0?(im>0?M_PIl:(im<0?-M_PIl:0)):atan(im/re));
    B.components[2*k]=sqrt(z)*cos(phi/2);
    B.components[2*k+1]=sqrt(z)*sin(phi/2);
  }
  return B;
}

const mwArray invertDiag(const mwArray& mwA,int& nr){
  return invertDiag(mwA,nr,0);
}

const mwArray invertDiag(const mwArray& mwA,int& nr,const double tol){
  if(mwA.isScalar()){ // particular case: inverse of the complex element
    complex_t elem=mwA.getElement(0);
    if(abs(elem)<tol) return mwArray(ZERO_c);
    else return mwArray(ONE_c/elem);
  }
  if(!mwA.isMatrix()||mwA.dims[0]!=mwA.dims[1]){
    cout<<"Error: cannot apply invertDiag to matrix of size "<<mwA.dims<<endl;
    exit(212);
  }
  double reT,imT;mwA.trace(reT,imT);
  double absT=sqrt(reT*reT+imT*imT);
  mwArray B(mwA.dims); // All zeros
  nr=0;
  for(int k=0;k<B.dims[0];k++){
    Indices pos(k,k);
    complex_t elem=mwA.getElement(pos);
    if(tol==0||abs(elem)/absT>tol){
      B.setElement(ONE_c/elem,pos);nr++;
    }
    else{
      B.setElement(ZERO_c,pos);
    }
    // double re,im;
    // mwA.getElement(re,im,pos);
    // double z=sqrt(re*re+im*im);
    // if(tol==0||z/absT>tol){
    //   double phi=atan(im/re);
    //   B.setElement((1/z)*cos(phi),-(1/z)*sin(phi),pos);
    //   nr++;
    // }
    // else{
    //   B.setElement(0.,0.,pos);
    // }
  }
  return B;
}


const mwArray inverse(const mwArray& mwA,int& nr){
  return inverse(mwA,nr,0.);
}


const mwArray inverse(const mwArray& mwA,int& nr,
		      const double tol){
  if(!mwA.isMatrix()||!mwA.isSquare()){
    cout<<"Error: cannot apply inverse to tensor of size "<<mwA.dims<<endl;
    exit(212);
  }
  if(mwA.isScalar()||mwA.isDiagonal())
    return invertDiag(mwA,nr,tol);
  // If not diagonal, svd and invert the diagonal
  mwArray U,S,Vdagger; // place for the results
  wrapper::svd(mwA,tol,nr,U,S,Vdagger);
  cout<<"mwArray::inverse found "<<nr<<" singular values"<<endl;
  cout<<" and sizes U "<<U.getDimensions()<<", S "<<S.getDimensions()
      <<", V+ "<<Vdagger.getDimensions()<<endl;
  S=invertDiag(S,nr,tol);
  U.Hconjugate(); // U+
  Vdagger.Hconjugate(); // V
  Vdagger.multiplyRight(S); // V*S
  Vdagger.multiplyRight(U); // V*S*U+
  return Vdagger;
}

bool mwArray::isDiagonal()const{
  if(rank==2){
    mwArray B(*this);
    for(int k=0;k<min(dims[0],dims[1]);k++){
      B.setElement(ZERO_c,Indices(k,k));
    }
    //cout<<"In isDiagonal, set auxiliary matrix B "<<B<<endl;
    double sum=0.;
    for(int k=0;k<nrcomponents;k++) 
      sum+=(B.components[2*k]+B.components[2*k+1]);
    if(sum==0.) return 1;
  }
  return 0;
}

void mwArray::changeSign(){
  for(int k=0;k<nrcomponents;k++){
    components[2*k]=-components[2*k];
    components[2*k+1]=-components[2*k+1];
  }
}

void mwArray::getIndices(int idx,Indices& result) const{
  if(idx<0||idx>=nrcomponents){
    cout<<"Invalid index "<<idx<<" for mwArray "<<*this<<" of size "
	<<nrcomponents<<endl;
    exit(212);
  }
  // Very little efficient: should be improved. boost??
  // Maybe just hardcode for small nr of dimensions
  //result.clear();int idx_=idx;
  result=Indices(vector<int>(rank,0));int idx_=idx;
  for(int k=0;k<rank;k++){
    int val=(idx_)%dims[k];
    //cout<<"Mod(idx,dim["<<k<<"]="<<dims[k]<<")="<<val<<endl;
    //result.push_back(val);
    result[k]=val;
    idx_=(idx_-val)/dims[k];
  }
}

void mwArray::_getindices_(int idx,Indices& result) const {
  // Very little efficient: should be improved. boost??
  // Maybe just hardcode for small nr of dimensions
  //  Indices aux(result);getIndices(idx,aux); //REMOVE!
  if(!computed) {
    cout<<"ERROR: getindices needs cumDim to be called first!"<<endl;
    exit(4);
  }
  int rem=idx;
  for(int k=rank-1;k>=0;k--){
    rem=idx%cumdims[k];
    result[k]=(idx-rem)/cumdims[k];
    idx=rem;
  }
  //cout<<"NEW getIndices("<<idx<<")->"<<result<<"; OLD "<<aux<<endl;
  //cout<<"getIndices("<<idx<<")->"<<result<<endl;
}


mwArray& mwArray::operator=(const mwArray& orig){
  if(this!=&orig){ // no self-assignment allowed
    if(mine){
      if(nrcomponents!=orig.getNrComponents()||components==0){
	if(components!=0) delete[] components;
	components=orig.copyComponents();
	nrcomponents=orig.getNrComponents();
      }
      else{
	memcpy((void*)components,(void*)orig.components,sizeof(double)*2*nrcomponents);
      }
    }
    else{ // can only assign same nr of elements
      if(nrcomponents!=orig.getNrComponents()){
	cout<<"Error: shallow mwArray can only be assigned to "
	    <<"the same number of components"<<endl;
	exit(212);
      }
      memcpy((void*)components,(void*)orig.components,sizeof(double)*2*nrcomponents);
    }
    rank=orig.getRank();
    orig.getDimensions(dims);
    computed=false;
  }
  return *this;
}

mwArray::mwArray(const mwArray& left,const mwArray& right,
		 const char oper):rank(0),nrcomponents(0),
				  components(0),computed(false),
				  mine(true),cumdims(),dims(){
  switch(oper){
  case '*':
    // Particular case: one scalar!
    //if((left.getRank()==1)&&(left.getDimension(0)==1)){
    if(left.isScalar()){
      //cout<<"Multiplying by scalar "<<alpha[0]<<"+"<<alpha[1]<<endl;
      *this=mwArray(left.getElement(0),right);
    }
    //else if((right.getRank()==1)&&(right.getDimension(0)==1)){
    else if(right.isScalar()){
      //cout<<"Multiplying by scalar "<<alpha[0]<<"+"<<alpha[1]<<endl;
      *this=mwArray(right.getElement(0),left);
    }
    else{
      if(left.getRank()>2||right.getRank()>2){
	cout<<"Error: operands to A*B do not have valid rank(2) "
	    <<"(A.rank="<<left.getRank()<<", B.rank="
	    <<right.getRank()<<")"<<endl;
	exit(212);
      }
      // Prepare place for the result: not needed: happens in wrapper
      // Call the procedure
      wrapper::product(left,right,*this);
    }
    break;
  case '+':
    if(left.rank!=right.rank){
      if(!left.isEmpty()&&!right.isEmpty()){
	cout<<"Error: incompatible dimensions in mwArray sum "<<left.dims
	    <<" and "<<right.dims<<endl;
	exit(212);
      }
      else{
	const mwArray* copying;
	if(left.isEmpty())
	  copying=&right;
	else
	  copying=&left;
	rank=copying->rank;copying->getDimensions(dims);
	nrcomponents=copying->nrcomponents;
	components=copying->copyComponents();
	break;
      }
    }
    for(int k=0;k<left.rank;k++){
      if(left.dims[k]!=right.dims[k]){
	cout<<"Error: incompatible dimensions in mwArray sum "<<left.dims
	    <<" and "<<right.dims<<endl;
	exit(212);
      }
    }
    //cout<<"Sum of matrices with dims "<<dimsA<<endl;
    rank=left.rank;left.getDimensions(dims);
    nrcomponents=left.nrcomponents;
    components=right.copyComponents();
    for(int k=0;k<2*nrcomponents;k++){
      components[k]+=left.components[k];
    }
    break;
  case '-':
    if(left.rank!=right.rank){
      if(!left.isEmpty()&&!right.isEmpty()){
	cout<<"Error: incompatible dimensions in mwArray substraction "<<left.dims
	    <<" and "<<right.dims<<endl;
	exit(212);
      }
      else{ // some of them is empty, then return the other! but the sing might change
	const mwArray* copying;
	if(left.isEmpty())
	  copying=&right;
	else
	  copying=&left;
	rank=copying->rank;copying->getDimensions(dims);
	nrcomponents=copying->nrcomponents;
	components=copying->copyComponents();
	if(!right.isEmpty())
	  changeSign();
	break;
      }
    }
    for(int k=0;k<left.rank;k++){
      if(left.dims[k]!=right.dims[k]){
	cout<<"Error: incompatible dimensions in mwArray substraction "<<left.dims
	    <<" and "<<right.dims<<endl;
	exit(212);
      }
    }
    //cout<<"Sum of matrices with dims "<<dimsA<<endl;
    rank=left.rank;left.getDimensions(dims);
    nrcomponents=left.nrcomponents;
    components=left.copyComponents();
    for(int k=0;k<2*nrcomponents;k++){
      components[k]-=right.components[k];
    }
    break;
  default:
    cout<<"Unknown operator in construction of mwArray"<<endl;
    exit(212);
  }
}

void mwArray::multiplyRight(const mwArray& right){
  if(right.isScalar()){
    complex_t c=right.getElement(0);
    for(int k=0;k<nrcomponents;k++)
      *(complex_t*)&components[2*k]*=c;
  }
  else // right is a proper tensor
    if(isScalar()){ // this is the scalar, and dims of result are as right
      Indices olddims=right.getDimensions();
      *this=getElement(0)*right;
      reshape(olddims);
    }
    else // Proper multiplication
      wrapper::product(mwArray(*this),right,*this);
}

void mwArray::multiplyLeft(const mwArray& left){
  if(left.isScalar()){
    complex_t c=left.getElement(0);
    for(int k=0;k<nrcomponents;k++)
      *(complex_t*)&components[2*k]*=c;
  }
  else
    if(isScalar()){
      Indices olddims=left.getDimensions();
      *this=getElement(0)*left;
      reshape(olddims);
      return;
    }
    else
      wrapper::product(left,mwArray(*this),*this);
}

mwArray::mwArray(const double alpha[2],const mwArray& mwA):
  rank(mwA.rank),dims(mwA.dims),nrcomponents(mwA.nrcomponents),
  computed(false),mine(true),cumdims(),components(NULL){
  components=new double[2*nrcomponents];
  //mwA.copyComponents();
  complex_t c=(complex_t){alpha[0],alpha[1]};
  for(int k=0;k<nrcomponents;k++){
    *(complex_t*)&components[2*k]=c*(*(complex_t*)&mwA.components[2*k]);
    //double u=alpha[0]*components[2*k]-alpha[1]*components[2*k+1];
    //double v=alpha[1]*components[2*k]+alpha[0]*components[2*k+1];
    //components[2*k]=u;
    //components[2*k+1]=v;
  }
}

mwArray::mwArray(const double alpha,const mwArray& mwA):
  rank(mwA.rank),dims(mwA.dims),nrcomponents(mwA.nrcomponents),
  computed(false),mine(true),cumdims(),components(NULL){
  //  components=mwA.copyComponents();
  components=new double[2*nrcomponents];
  for(int k=0;k<nrcomponents;k++){
    *(complex_t*)&components[2*k]=alpha*(*(complex_t*)&mwA.components[2*k]);
    //    components[2*k]=alpha*components[2*k];
    //components[2*k+1]=alpha*components[2*k+1];
  }
}

mwArray::mwArray(const complex_t& alpha,const mwArray& mwA):
  rank(mwA.rank),dims(mwA.dims),nrcomponents(mwA.nrcomponents),
  computed(false),mine(true),cumdims(),components(NULL){
  //  components=mwA.copyComponents();
  components=new double[2*nrcomponents];
  for(int k=0;k<nrcomponents;k++){
    *(complex_t*)&components[2*k]=alpha*(*(complex_t*)&mwA.components[2*k]);
    //    complex_t old=*(complex_t*)&components[2*k];
    //*(complex_t*)&components[2*k]=alpha*old;
  }
}
// Multiply two matrices (only works for rank 2 arrays)
const mwArray operator*(const mwArray& left,const mwArray& right){
  return mwArray(left,right,'*');
}

// Add two arrays of the same dimensions
const mwArray operator+(const mwArray& left,const mwArray& right){
  return mwArray(left,right,'+');
}
// Substract two arrays of the same dimensions
const mwArray operator-(const mwArray& left,const mwArray& right){
  return mwArray(left,right,'-');
}

// Multiply by complex scalar
const mwArray operator*(const double* alpha,const mwArray& mwA){
  return mwArray(alpha,mwA);
}
//Multiply by real scalar
const mwArray operator*(const double alpha,const mwArray& mwA){
  return mwArray(alpha,mwA);
}

const mwArray operator*(const complex_t& alpha,const mwArray& mwA){
  return mwArray(alpha,mwA);
}

bool operator==(const mwArray& left,const mwArray& right){
  if((left.rank!=right.rank)||(left.nrcomponents!=right.nrcomponents)|| 
     (left.dims!=right.dims))
    return false;
  // if all metadata coincide, check elements one by one
  for(int k=0;k<2*left.nrcomponents;k++){
    if(left.components[k]!=right.components[k]) return false;
  }
  // if passed all checks, they are equal
  return true;
}

bool mwArray::isNull(double tol) const{
  for(int k=0;k<nrcomponents;k++){
    complex_t val=*((complex_t*)&components[2*k]);
    if(abs(val)>tol){
      //cout<<"isNull check failing, for tol "<<tol
      // <<", due to matrix element "<<val<<endl;
      return false;
    }
  } 
  return true;
}

void mwArray::getMaxElement(complex_t& value,shrt::Indices& indices) const{
  double maxV=0;int pos=0;
  for(int k=0;k<nrcomponents;k++){
    complex_t val=*((complex_t*)&components[2*k]);
    if(abs(val)>maxV){
      maxV=abs(val);value=val;pos=k;
    }
  } 
  getIndices(pos,indices);
}

bool mwArray::isHermitian(double tol) const{
  return (::Hconjugate(*this)-(*this)).isNull(2*tol);
}

const complex_t scalarproduct(const mwArray& A,const mwArray& B){
  if((!A.isVector()||!B.isVector())||
     (A.getNrComponents()!=B.getNrComponents())){
    cout<<"Error: cannot compute scalarproduct from dimensions "<<A.dims
	<<" and "<<B.dims<<endl;
    exit(212);
  }
  complex_t tmp=ZERO_c;
  for(int k=0;k<A.getNrComponents();k++){
    tmp=tmp+(*(complex_t*)&A.components[2*k])*(*(complex_t*)&B.components[2*k]);
  }
  return tmp;
}

#include "Indices.h"
 
const mwArray outerproduct(const mwArray& A,const mwArray& B){
  if(!A.isVector()||!B.isVector()){
    cout<<"Error: cannot compute outerproduct from dimensions "<<A.dims
	<<" and "<<B.dims<<endl;
    exit(212);
  }
  int dimA=A.getNrComponents();
  int dimB=B.getNrComponents();
  mwArray C(Indices(dimA,dimB));
  mwArray auxA=reshape(A,Indices(dimA));
  mwArray auxB=reshape(B,Indices(dimB));
  for(int kA=0;kA<dimA;kA++)
    for(int kB=0;kB<dimB;kB++){
      C.setElement(auxA.getElement(shrt::Indices(kA))*
		   auxB.getElement(shrt::Indices(kB)),
		   shrt::Indices(kA,kB));
    }
  return C;
}

const mwArray kron(const mwArray& A,const mwArray& B){
  if(A.getRank()!=2||B.getRank()!=2){
    cout<<"Error: cannot apply kron to A"<<A.getDimensions()
	<<" and B"<<B.getDimensions()<<", only matrix,matrix"<<endl;
    exit(212);
  }
  cout<<"kron A"<<A.getDimensions()<<"xB"<<B.getDimensions()<<endl;
  Indices dimsA=A.getDimensions();
  Indices dimsB=B.getDimensions();
  mwArray result=reshape(A,Indices(-1,1))*reshape(B,Indices(1,-1));
  result.reshape(Indices(dimsA[0],dimsA[1],dimsB[0],dimsB[1]));
  result.permute(Indices(1,3,2,4));
  result.reshape(Indices(dimsA[0]*dimsB[0],dimsA[1]*dimsB[1]));
  cout<<"result"<<result.getDimensions()<<endl;
  return result;
}

const mwArray kron(mwArray& A,mwArray& B){
  if(A.getRank()!=2||B.getRank()!=2){
    cout<<"Error: cannot apply kron to A"<<A.getDimensions()
	<<" and B"<<B.getDimensions()<<", only matrix,matrix"<<endl;
    exit(212);
  }
  if(&A==&B){
    cout<<"kron receiving the same argument twice, should copy one of them!"
	<<endl;
    return kron(A,mwArray(B));
  }
  //cout<<"Al llegar a kron: A->"<<A.getDimensions()<<", B->"
  //  <<B.getDimensions()<<endl;
  Indices dimsA=A.getDimensions();
  Indices dimsB=B.getDimensions();
  A.reshape(Indices(-1,1));
  B.reshape(Indices(1,-1));
  //  cout<<"kron ha hecho un reshape: A->"<<A.getDimensions()<<", B->"
  //  <<B.getDimensions()<<endl;
  mwArray result=A*B;
  result.reshape(Indices(dimsA[0],dimsA[1],dimsB[0],dimsB[1]));
  result.permute(Indices(1,3,2,4));
  result.reshape(Indices(dimsA[0]*dimsB[0],dimsA[1]*dimsB[1]));
  A.reshape(dimsA);
  B.reshape(dimsB);
  return result;
}

const mwArray hadamard(const mwArray& A,const mwArray& B){
  // Check dimensions
  if(A.getDimensions()!=B.getDimensions()){
    cout<<"Error: cannot do element-wise multiplication of two arrays "
	<<"with different dimensions A("<<A.getDimensions()<<"), B("
	<<B.getDimensions()<<")"<<endl;
    exit(1);
  }
  mwArray result=A;
  for(int k=0;k<result.getNrComponents();k++){
    *(complex_t*)&result.components[2*k]*=*(complex_t*)&B.components[2*k];
  }
  return result;
}


void mwArray::trace(double& re,double& im) const {
  complex_t c=trace();
  re=c.re;im=c.im;
  /*if(rank==2&&dims[0]==dims[1]){
    re=0.;im=0.;
    for(int k=0;k<dims[0];k++){
      double re_,im_;
      Indices pos(2,k);
      getElement(re_,im_,pos);
      re+=re_;im+=im_;
    }
  }
  else{
    cout<<"Error: cannot compute trace of array of size "<<dims<<endl;
    exit(212);
    }*/
}

complex_t mwArray::trace() const{
  if(isScalar()) return *((complex_t*)&components[0]);
  if(rank==2&&dims[0]==dims[1]){
    complex_t result(ZERO_c);
    for(int k=0;k<dims[0];k++){
      result=result+getElement(Indices(k,k));
    }
    return result;
  }
  else{
    cout<<"Error: cannot compute trace of array of size "<<dims<<endl;
    exit(212);
  }
}

#define LAST 111 // large enough that there are not such big rank mwArrays
// Auxiliary function to run over all the dimensions, increasing the
// rightmost index by one every time. Argument cur may be used to
// specify a different dimension to increase. If that's reached the
// maximum, the previous one is increased, and the new set of indices
// is contained in the vector, while the return value is 1. If the
// values in the vector correspond to the last component, the vector
// is not changed and -1 is returned.
int increaseIndices(Indices& indices,const Indices& dims,int cur=LAST){
  //cout<<"increaseIndices "<<indices<<"; dimensions"<<dims<<endl;
  if(indices.size()!=dims.size()){ 
    cout<<"Invalid indices size! "<<indices<<" for dimensions "
	<<dims<<endl;
    return -1;
  }
  if(cur==LAST) cur=indices.size()-1;
  if(cur<0||cur>dims.size()-1) return -1;
  if(indices[cur]==dims[cur]-1)
    return increaseIndices(indices,dims,cur-1);
  else{
    //cout<<"increaseIndices ("<<indices;
    indices[cur]=indices[cur]+1;
    for(int r=cur+1;r<dims.size();r++){
      indices[r]=0;
    }
    //cout<<")->("<<indices<<"), at index "<<cur+1<<endl;
    return 1;
  }

}

#include <iomanip>

ostream& putForMatlab(ostream& os,const mwArray& mwA,const char* name){
  return putForMatlab(os,mwA,name,0);
}

ostream& putForMatlab(ostream& os,const mwArray& mwA,const char* name,int precision){
  if(precision>0)
    os<<setprecision(precision);
  if(mwA.isScalar()){    
    return os<<name<<"="<<mwA.getElement(0)<<";";
  }
  if(mwA.nrcomponents==0) return os<<name<<"=[];";
  // Show components as matlab does
  Indices indices(vector<int>(mwA.rank,0));
  if(mwA.rank>=2){
    int curidx=2;
    int next=1;
    while(next>0){
      os<<name<<"(:,:,";
      for(int ll=2;ll<mwA.rank-1;ll++) os<<indices[ll]+1<<",";
      os<<indices[mwA.rank-1]+1<<")=[\n\t";
      // Print the 12 matrix
      // Reset the i1,i2 indices
      indices[0]=0;indices[1]=0;
      for(int p1=0;p1<mwA.dims[0];p1++){
	indices[0]=p1;
	for(int p2=0;p2<mwA.dims[1];p2++){
	  double re,im;indices[1]=p2;
	  mwA.getElement(re,im,indices);
	  //os<<"("<<indices[0]<<","<<indices[1]<<")";
	  os<<re;
	  if(im>=0)
	    os<<"+"<<im<<"i, ";
	  else
	    os<<"-"<<-im<<"i, ";
	}
	if(p1==mwA.dims[0]-1) 	os<<"];\n\t";
	else	os<<";\n\t";
      }
      os<<"\n";
      //cout<<"About to call increaseIndices from "<<indices<<endl;
      next=increaseIndices(indices,mwA.dims,LAST);
    }
  }
  else{ // raank is only one
    os<<name<<"=[";
    for(int p=0;p<mwA.dims[0];p++){
      indices[0]=p;double re,im;
      mwA.getElement(re,im,indices);
      os<<re;
      if(im>=0)
	os<<"+"<<im<<"i, ";
      else
	os<<"-"<<-im<<"i, ";
    }
    os<<"];\n";
  }
  return os;
}



ostream& operator<<(ostream& os,const mwArray& mwA){
  if(mwA.isScalar()){
    return os<<mwA.getElement(0);
  }
  os<<"mwArray dims"<<mwA.dims<<", N_elements="
    <<mwA.nrcomponents<<"\n";
  if(mwA.nrcomponents==0) return os;
  os<<" components=\n";
  // Show components as matlab does
  Indices indices(vector<int>(mwA.rank,0));
  if(mwA.rank>=2){
    int curidx=2;
    int next=1;
    while(next>0){
      os<<"(:,:,";
      for(int ll=2;ll<mwA.rank-1;ll++) os<<indices[ll]+1<<",";
      os<<indices[mwA.rank-1]+1<<")=[\n\t";
      // Print the 12 matrix
      // Reset the i1,i2 indices
      indices[0]=0;indices[1]=0;
      for(int p1=0;p1<mwA.dims[0];p1++){
	indices[0]=p1;
	for(int p2=0;p2<mwA.dims[1];p2++){
	  double re,im;indices[1]=p2;
	  mwA.getElement(re,im,indices);
	  //os<<"("<<indices[0]<<","<<indices[1]<<")";
	  os<<re;
	  if(im>=0)
	    os<<"+"<<im<<"i, ";
	  else
	    os<<"-"<<-im<<"i, ";
	}
	if(p1==mwA.dims[0]-1) 	os<<"];\n\t";
	else	os<<";\n\t";
      }
      os<<"\n";
      //cout<<"About to call increaseIndices from "<<indices<<endl;
      next=increaseIndices(indices,mwA.dims,LAST);
    }
  }
  else{
    for(int p=0;p<mwA.dims[0];p++){
      indices[0]=p;double re,im;
      mwA.getElement(re,im,indices);
      os<<re;
      if(im>=0)
	os<<"+"<<im<<"i, ";
      else
	os<<"-"<<-im<<"i, ";
    }
  }
  return os;
}

void mwArray::fillRandom(){
  if(rank>0){
    for(int k=0;k<nrcomponents;k++){
      components[2*k]=random()*1./RAND_MAX;
      components[2*k+1]=random()*1./RAND_MAX;
    }
  }
}

void mwArray::resize(const Indices& newdims,double noise){
  if(newdims==dims) return; // nothing to do
  if(!mine){
    cout<<"Error: cannot resize a shallow array "<<endl;
    exit(212);
  }
  if(rank==0){ // Create empty matrix and return
    rank=newdims.size();
    dims=newdims;
    nrcomponents=1;
    for(int k=0;k<rank;k++)nrcomponents=nrcomponents*dims[k];
    components=new double[2*nrcomponents];
    fillWithZero();
    cumdims.clear();computed=false;
    return;
  }
  if(newdims.size()!=rank){
    cout<<"Error: invalid vector of dimensions "<<newdims
	<<" for resize (current "<<dims<<")"<<endl;
    exit(212);
  }
  int Nnrcomponents=1;
  Indices newcumdims(vector<int>(rank,1));
  for(int k=0;k<rank;k++){
    Nnrcomponents=Nnrcomponents*newdims[k];
    if(k>0) newcumdims[k]=newcumdims[k-1]*newdims[k-1];
  }
  double* newcomp=new double[Nnrcomponents*2];
  for(int k=0;k<Nnrcomponents;k++){
    if(noise!=0){
      newcomp[2*k]=noise*random()*1./RAND_MAX;
      newcomp[2*k+1]=noise*random()*1./RAND_MAX;
    }
    else{
      newcomp[2*k]=0.;
      newcomp[2*k+1]=0.;
    }
  }  
  //cout<<"RESIZE: olddims "<<dims<<" -> "<<newdims<<endl;
  // read them all and write, one by one, to the new positions
  Indices oldidx(vector<int>(rank,0));
  //Indices _oldidx(rank,0);
  int k=0;
  while(k<nrcomponents){
    //for(int k=0;k<nrcomponents;k++){
    //_getindices_(k,_oldidx);
    if(oldidx<newdims){ // if within the new dimensions
      //int newpos=_linearidx_(oldidx,newdims);
      int newpos=0;
      for(int p=0;p<rank;p++)
	newpos+=oldidx[p]*newcumdims[p];
      //cout<<"RESIZE oldidx "<<_oldidx<<", goes to pos "<<_linearidx_(_oldidx,newdims)<<endl;
      //cout<<"NOW oldidx="<<oldidx<<", goes to "<<newpos<<endl;
      newcomp[2*newpos]=components[2*k];
      newcomp[2*newpos+1]=components[2*k+1]; 
    }
    //else cout<<"Element "<<oldidx<<" is cut out"<<endl;
    k++;incrementOne(oldidx,dims,rank);   
  }
  delete[] components;
  components=newcomp;
  dims=newdims;
  nrcomponents=Nnrcomponents;
}

void mwArray::savetext(ofstream& data,bool real) const {
  data<<rank<<" ";
  for(int k=0;k<rank;k++){
    data<<dims[k]<<" ";
  }
  if(!real)
    for(int k=0;k<2*nrcomponents;k++){
      data<<components[k]<<" ";
      //cout<<"written comp["<<k<<"]="<<components[k]<<endl;
    }
  else
    for(int k=0;k<nrcomponents;k++){
      data<<components[2*k]<<" ";
      //cout<<"written comp["<<k<<"]="<<components[k]<<endl;
    }
  //cout<<"Saved mwArray "<<dims<<endl;
}

void mwArray::loadtext(ifstream& data,bool real){
  // Discard previous contents
  dims.clear();int newnrcomponents=1;
  if(mine)
    if(components!=0) delete[] components;
  data>>rank;
  cout<<"load read rank="<<rank<<endl;
  for(int k=0;k<rank;k++){
    int tmpD;
    data>>tmpD;
    dims.push_back(tmpD);newnrcomponents=newnrcomponents*tmpD;
  }
  if(!mine){
    if(newnrcomponents!=nrcomponents){
      cout<<"Error: cannot load data from text into shallow mwArray "
	  <<"unless size is the same "<<endl;
      exit(212);
    }
  }
  else{
    nrcomponents=newnrcomponents;
    components=new double[2*nrcomponents];
  }
  if(!real)
    for(int k=0;k<2*nrcomponents;k++)
      data>>components[k];
  else
    for(int k=0;k<nrcomponents;k++){
      data>>components[2*k];
      components[2*k+1]=0.;
    }
}

void mwArray::setPointer(const shrt::Indices& dims_,double* ptr){
  clear();
  components=ptr;
  dims=dims_;
  rank=dims.size();
  nrcomponents=1;
  for(int k=0;k<rank;k++) nrcomponents*=dims[k];
  mine=false;
}

void mwArray::contractLeg(int k1,const mwArray& B,int k2){
  mwArray aux(B);
  contractLeg(k1,aux,k2);
}

void mwArray::contractLeg(int k1,mwArray& B,int k2){
  // read dimensions, to restore them
  Indices dimsB=B.getDimensions();
  Indices dimsA=dims;
  k1--;k2--; // they come from 0 to n-1
#ifdef CHECKDIMS
  if(dims[k1]!=B.getDimension(k2)){
    cout<<"Error in contractLeg: dimensions to be contracted do not agree"
	<<endl;
  }
#endif
  //cout<<"contractLeg with A->"<<dimsA<<", B->"<<dimsB<<", k1="
  //  <<k1<<", k2="<<k2<<endl;
  // First, reshape
  int dimsLA=1,dimsRA=1,dimsLB=1,dimsRB=1;
  for(int k=0;k<max(rank,B.getRank());k++){
    if(k<rank){
      if(k<k1) dimsLA*=dimsA[k];
      if(k>k1) dimsRA*=dimsA[k];
    }
    if(k<B.getRank()){
      if(k<k2) dimsLB*=dimsB[k];
      if(k>k2) dimsRB*=dimsB[k];
    }
  }
  //cout<<"Will reshape A to ("<<dimsLA<<","<<dimsA[k1]<<","<<dimsRA<<")"<<endl;
  reshape(Indices(dimsLA,dimsA[k1],dimsRA));
  //cout<<"Will reshape B to ("<<dimsLB<<","<<dimsB[k2]<<","<<dimsRB<<")"<<endl;
  B.reshape(Indices(dimsLB,dimsB[k2],dimsRB));
  // permute indices (there may be singleton dimensions, now, but that
  // should not be a problem
  // We should decide here how to do the product, but maybe we would
  // need extra permutations later to keep consistent ordering of the
  // indices.
  permute(Indices(1,3,2));
  B.permute(Indices(2,1,3));
  reshape(Indices(dims[0]*dims[1],dimsA[k1]));
  //cout<<"Now reshape A to "<<dims<<endl;
  B.reshape(Indices(dimsB[k2],B.getDimension(1)*B.getDimension(2)));
  //cout<<"And B to "<<B.getDimensions()<<endl;
  // Do the product
  multiplyRight(B);
  // Reshape again
  vector<int> newdimsA;
  for(int k=0;k<dimsA.size();k++)
    if(k!=k1) newdimsA.push_back(dimsA[k]);
  for(int k=0;k<dimsB.size();k++)
    if(k!=k2) newdimsA.push_back(dimsB[k]);
  reshape(Indices(newdimsA));
}

const mwArray contractLeg(const mwArray& A,int k1,
			  const mwArray& B,int k2){
  mwArray C(A);
  C.contractLeg(k1,B,k2);
  return C;
}

 const mwArray contractLegs(const mwArray& A,const Indices& k1,
			    const mwArray& B,const Indices& k2){
  mwArray C(A);
  C.contractLegs(k1,B,k2);
  return C;
 }

void mwArray::contractLegs(const Indices& k1,
		  const mwArray& B,const Indices& k2){
  mwArray aux(B);
  contractLegs(k1,aux,k2);
}

// Auxiliary function to reorder the list of indices 0 to rank-1 
// putting the selected ones (lista) at the end.
// It also returns the list of new dimensions according to that order,
// and the product of the segregated ones
// If front==true, the list dimensions are returned in the front 
//Needs algorithm library 
#include <algorithm>
void reSort(const Indices& oldDims,const Indices& lista,
	    Indices& newDims,Indices& newL,int& prodL,int& prodR,
	    bool front=false){
  //cout<<"reSort receives original "<<oldDims<<" of which "<<lista
  //  <<" are to be contracted"<<endl;
  int rank=oldDims.size();
  int sizeL=lista.size();
  newL.clear();newDims.clear();
  int ptr=0;
  int cnt=0;
  prodL=1;
  
  for(cnt=0;cnt<rank; cnt++)
  {
    //Check if current index is one which is going to be contracted, if not put it in the list of remaining indices
    if (std::find(lista.begin(),lista.end(),cnt)== lista.end())
    {
      newL.push_back(cnt+1);   
      newDims.push_back(oldDims[cnt]);
      prodL*=oldDims[cnt];
    }      
  }
  //Old code, can run into trouble if the entries in lista are not in ascending order
  /*while(cnt<rank){
    if(cnt==lista[ptr]){
      if(ptr<sizeL)
	ptr++;
    }
    else{
      newL.push_back(cnt+1);   
      newDims.push_back(oldDims[cnt]);
      prodL*=oldDims[cnt];
    }
    cnt++;
  }*/
  // now, add all the lista
  ptr=0;prodR=1;
  Indices aux;
  while(ptr<lista.size()){
    //newDims.push_back(oldDims[lista[ptr]]);
    prodR*=oldDims[lista[ptr]];
    if(!front)
      newL.push_back(lista[ptr++]+1);
    else
      aux.push_back(lista[ptr++]+1);
  }
  if(front){
    newL=Indices(aux,newL);
    int tmp=prodL;
    prodL=prodR;prodR=tmp;
  }
  //cout<<"reSort returns newdims "<<newDims<<" newList "<<newL
  //  <<"products(L="<<prodL<<",R="<<prodR<<")"<<endl;
}

void mwArray::contractLegs(const Indices& k1_,
			   mwArray& B,const Indices& k2_){
   // First, check lengths of indices to contract
   // TODO: would it be better to check only the total dimension?
  // Convert numbers to range 0-n-1
  Indices k1(k1_),k2(k2_);
  k1.increment(-1);
  k2.increment(-1);
#ifdef CHECKDIMS
   bool error=0;
   if(k1.size()!=k2.size()){
     cout<<"contractLegs ERROR: length of lists of Indices are not the same: "
	 <<k1<<", and "<<k2<<endl;
     exit(1);
   }
   for(int k=0;k<k1.size()&&error==0;k++){
     if(dims[k1[k]]!=B.getDimension(k2[k])) error=1;
   }
   if(error){
     cout<<"Incompatible lists of indices to contract "<<k1<<" and "<<k2<<endl;
     exit(1);
   }
#endif
   // Reorder indices, so that all the ones in k1 are together at the
   // end of A, and all the ones in k2 are together at the beginning
   // of B
   // TODO: better to do in several steps, so that it is always a
   // transposition? 
   Indices newOrA,newDimsA;
   Indices newOrB,newDimsB;
   int prodL1,prodL2;
   int prodR1,prodR2;
   reSort(dims,k1,newDimsA,newOrA,prodL1,prodR1);
   reSort(B.getDimensions(),k2,newDimsB,newOrB,prodL2,prodR2,true);
   // permute the tensors according to the list
   permute(newOrA);reshape(Indices(prodL1,prodR1));
   B.permute(newOrB);B.reshape(Indices(prodL2,prodR2));
   //   B.permute(newOrB);B.reshape(Indices(prodL2,prodR2));
   //B.permute(Indices(2,1));
   multiplyRight(B);
   reshape(Indices(newDimsA,newDimsB));
}


void constructOperatorProduct(mwArray& result,const mwArray& opA,
			      const mwArray& opB) {
  if(opA.getRank()!=2||opB.getRank()!=2){
    cout<<"Error: cannot only apply constructOperatorProduct to matrices. Here A"<<opA.getDimensions()
	<<" and B"<<opB.getDimensions()<<endl;
    exit(1);
  }
  Indices dimsA=opA.getDimensions();
  Indices dimsB=opB.getDimensions();
  result=opA; 
  result.reshape(Indices(dimsA[0]*dimsA[1],1));
  result.multiplyRight(reshape(opB,Indices(1,dimsB[0]*dimsB[1])));
  result.reshape(Indices(dimsA[0],dimsA[1],dimsB[0],dimsB[1])); // order here: iji'j'
  result.permute(Indices(1,3,2,4));  // order in the MPO: ii',jj'  
  result.reshape(Indices(dimsA[0]*dimsB[0],dimsA[1]*dimsB[1]));
}



