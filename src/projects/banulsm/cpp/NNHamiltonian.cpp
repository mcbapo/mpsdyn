#include "NNHamiltonian.h"
#include <iomanip>

// Tolerance for svd of Hamiltonian terms
#define TOLSVD 1E-10 

using namespace std;
using namespace shrt;

NNHamiltonian::NNHamiltonian():uptodate(0),mpoH(0){};

NNHamiltonian::NNHamiltonian(int L_,const vector<int>& dims):L(L_),uptodate(0),mpoH(L_){
  if(dims.size()!=L_){
    cout<<"Error: argument dims for creation of NNHamiltonian does not"
	<<" have the right length"<<endl;
    exit(1);
  }
  for(int k=0;k<L;k++){
    d.push_back(dims[k]);
    if(k<L-2) terms.push_back(mwArray());
  }
  mpoH.initLength(L);
};

NNHamiltonian::NNHamiltonian(const char* filename):mpoH(0){
  ifstream input(filename);
  if(!input.is_open()){
    cout<<"Error: cannot open file "<<filename<<" to read NNHamiltonian"<<endl;
    exit(1);
  }
  input>>L;
  input>>setprecision(20);
  d=vector<int>(L,0);
  mpoH.initLength(L);
  cout<<"Read number sites "<<L<<endl;
  for(int k=0;k<L;k++){
    input>>d[k];
    cout<<"Read dimension d["<<k<<"]="<<d[k]<<endl;
    if(k<L-1){ 
      terms.push_back(mwArray());
      //cout<<"pushed back term["<<k<<"]="<<terms[k]<<endl;
    }
  }
  for(int k=0;k<L-1;k++){
    //    cout<<"Going to read term "<<k<<" for dims "<<d[k]<<","<<d[k+1]<<endl;
    // // int ndim;
    // // input>>aux;
    // // ndim=(int)aux;
    // // Indices dims(1,ndim);
    // // for(int l=0;l<ndim;l++){
    // //   input>>aux;
    // //   dims[l]=(int)aux;
    // // }
    // // cout<<"Read dimensions of term "<<k<<", "<<dims<<endl;
    mwArray tmpT;
    tmpT.loadtext(input);
    // Correct for Hermiticity?
    // I assume term comes as d1xd2  d1xd2
    mwArray aux(tmpT);aux.permute(Indices(3,4,1,2),true);
    tmpT=.5*(aux+tmpT);
    terms[k]=tmpT;
    //cout<<"Read term "<<terms[k]<<endl;
    cout<<"Read term ("<<k<<") "<<terms[k].getDimensions()<<endl;
  }
  input.close();
  uptodate=false;
  mpoH.initLength(L);
  initHMPO();
}
  

NNHamiltonian::~NNHamiltonian(){
  d.clear();terms.clear();mpoH.clear();
}

void NNHamiltonian::getDimensions(std::vector<int>& dims) const{
  if(!d.empty()){
    dims=d;
  }
  else{
    cout<<"Warning: copying empyty dimensions"<<endl;
    dims.clear();
  }
}

const MPO& NNHamiltonian::getHMPO() const {
  if(d.empty()){
    cout<<"Error: cannot construct the MPO for empty NNHamiltonian"<<endl;
    exit(1);
  }
  if(!uptodate){
      cout<<"Error: MPO for NNHamiltonian is not computed. Function recomputeMPO should be called first!"<<endl;
      exit(1);
    //initHMPO();
  }
  return mpoH;
}

void NNHamiltonian::recomputeMPO(){
    if(d.empty()){
        cout<<"Error: cannot construct the MPO for empty NNHamiltonian"<<endl;
        exit(1);
    }
 //Should check also that no mwArray is empty!
    for (int l=0; l<L-1; l++) {
        if(terms[l].isEmpty()){
            cout<<"Error: cannot construct the MPO because term nr "<<l<<" is empty!"<<endl;
            exit(1);
        }
    }
    initHMPO();
}

mwArray NNHamiltonian::getTerm(int pos) const{
  // Check position
  if(pos>=L||d.size()<=pos+1){
    cout<<"Error at NNHamiltonian::getTerm->No term for site "<<pos<<endl;
    exit(1);
  }
  return terms[pos];
}

void NNHamiltonian::setTerm(int pos,const mwArray& term){
  // 1) check dimensions
  if(pos>=L||d.size()<=pos+1){
    cout<<"Error: NNHamiltonian::setTerm cannot set term for site "<<pos<<endl;
    exit(1);
  }
  int d1=d[pos];
  int d2=d[pos+1];
  Indices dimsT=term.getDimensions();
  if(dimsT.size()!=2||dimsT[0]!=d1*d2||dimsT[1]!=d1*d2){
    cout<<"Error: wrong dimensions for Hamiltonian term for "
	<<"position "<<pos<<"->"<<term<<endl;
    exit(1);
  }
  // 2) store the term in the list
  if(terms.size()<pos){ // not yet initialized?
    mwArray aux; // empty matrix
    while(terms.size()<=pos)
      terms.push_back(aux);
  }
  terms[pos]=term; // set the given term
  uptodate=0; // something changed, so the MPO is not ready
}

void NNHamiltonian::getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,int pos){
  //  cout<<"NNHamiltonian::getTwoBodyTermExponential(delta="<<delta<<", pos="<<pos<<")"<<endl;
  if(pos>=L||d.size()<=pos+1){
    cout<<"Error: Not enough dimensions specified for site "<<pos<<endl;
    exit(1);
  }
  mwArray term=terms[pos];
  int d1=d[pos];int d2=d[pos+1];
  // 1) prepare the operators for the mpoH
  // Reshape according to the svd to be applied
  term.reshape(Indices(d1*d2,d1*d2));
  if(!term.isHermitian(1E-15)){
    cout<<"WARNING: Term for pos "<<pos<<" is not exactly Hermitian!!!"<<endl;
    term=.5*(term+Hconjugate(term));
    //cout<<term-Hconjugate(term)<<endl;
  }
  //cout<<"getTwoBodyTermExponential for term ("<<pos<<") "<<term<<endl;exit(2);
  mwArray expH;
  wrapper::expm(term,expH,delta);
  cout<<"getTwoBodyTermExponential for expm ("<<pos<<") "<<endl;
  //  cout<<"After expm, got "<<expH.getDimensions()<<endl;
  // Obtained with indices (i'j')(ij) => permute to (i'i)(j'j) for svd
  expH.reshape(Indices(d1,d2,d1,d2));
  expH.permute(Indices(1,3,2,4));
  expH.reshape(Indices(d1*d1,d2*d2));
  int nr; // place for result of the SVD
  splitTwoBodyTerm(expH,d1,d2,Ol,Or,nr);
  // cout<<"getTwoBodyTermExponential(pos="<<pos<<", delta="<<delta<<") returns Ol "
  //  <<Ol.getDimensions()<<", Or "<<Or.getDimensions()<<endl;
}

// TODO: What if one term is not there!
void NNHamiltonian::initHMPO(){
  // Run over the terms splitting each one and setting the operators
  // I do not know what to do if a term is missing, so I check:
  for(int k=0;k<L-1;k++){
    if(terms[k].isEmpty()){
      cout<<"Error: empty term in the Hamiltonian not yet supported"<<endl;
      exit(1);
    }
  }
  mwArray Or_; // the right term from the one before
  int Dl=1;
  int nr; // nr of singular values in between two sites
  for(int k=0;k<L-1;k++){
    //cout<<"initHMPO(), term "<<k<<endl;
    mwArray term=terms[k];
    int d1=d[k];int d2=d[k+1];
    // 1) prepare the operators for the mpoH
    // Reshape according to the svd to be applied
    term.reshape(Indices(d1,d2,d1,d2));
    term.permute(Indices(1,3,2,4));
    term.reshape(Indices(d1*d1,d2*d2));
    mwArray Ol,Or;nr=0;
    splitTwoBodyTerm(term,d1,d2,Ol,Or,nr);
    // 2) Set the components of the MPO
    int Dr=nr+2; // right direction 
    mwArray oper(Indices(d1,Dl,d1,Dr));
    //cout<<"For site "<<k<<", place "<<oper.getDimensions()<<endl;
    for(int da=0;da<d1;da++){ // nothing or all
      if(k<L-1) oper.setElement(ONE_c,Indices(da,0,da,0));
      if(k>0) oper.setElement(ONE_c,Indices(da,Dl-1,da,Dr-1));
      for(int al=0;al<nr;al++)
	for(int dp=0;dp<d1;dp++)
	  oper.setElement(Ol.getElement(Indices(da,0,dp,al)),
			  Indices(da,0,dp,al+1));
      if(k!=0){ // then I have a Or term
	for(int al=0;al<Dl-2;al++)
	  for(int dp=0;dp<d1;dp++)
	  oper.setElement(Or_.getElement(Indices(da,al,dp,0)),
			  Indices(da,al+1,dp,Dr-1));
      }
    }
    // Set the operator in the MPO
    //cout<<"Setting term "<<k<<" to "<<oper.getDimensions()<<endl;
    mpoH.setOp(k,new Operator(oper),true);
    Or_=Or; // for the next one
    Dl=Dr;
  }
  // For the last one, only Or
  int d1=d[L-1];int Dr=1;
  mwArray oper(Indices(d1,Dl,d1,Dr));
  //cout<<"For last site, place "<<oper.getDimensions()<<endl;
  for(int da=0;da<d1;da++){
    oper.setElement(ONE_c,Indices(da,Dl-1,da,Dr-1));
    for(int al=0;al<nr;al++)
      for(int dp=0;dp<d1;dp++)
	oper.setElement(Or_.getElement(Indices(da,al,dp,0)),
			Indices(da,al+1,dp,Dr-1));
  }
  //cout<<"initHMPO(), setting last term "<<oper.getDimensions()<<endl;
  mpoH.setOp(L-1,new Operator(oper),true);
  uptodate=true;
}

void NNHamiltonian::splitTwoBodyTerm(const mwArray& term,int d1,int d2,
				     mwArray& Ol, mwArray& Or,int& nr) const {
  mwArray S; // place for the singular values
  double tol=TOLSVD;
  //cout<<"splitTwoBodyTerm received tensor "<<term.getDimensions()<<endl;
  wrapper::svd(term,tol,nr,Ol,S,Or);
  //cout<<"splitTwoBodyTerm found tensors Ol
  //"<<Ol.getDimensions()<<"\n Or "<<Or.getDimensions()<<"\n and S
  //"<<S.getDimensions()<<endl;
  //cout<<"splitTwoBodyTerm found singular values "<<S<<endl;
  // redistribute the singular values
  S=sqrt(S);
  // TODO: check no NaN???
  Ol.multiplyRight(S);
  Or.multiplyLeft(S);
  // reshape for operators (i'x 1 x i x alpha ))
  Ol.reshape(Indices(d1,1,d1,nr));
  Or.reshape(Indices(nr,d2,d2,1));
  Or.permute(Indices(2,1,3,4));  
  //cout<<"splitTwoBodyTerm returns Ol "<<Ol.getDimensions()<<", Or "<<Or.getDimensions()<<endl;
}

void NNHamiltonian::getExponentialMPO(MPO& expH,complex_t delta,bool even){
  cout<<"NNHamiltonian::getExponentialMPO(delta="<<delta<<", even="<<even<<")"<<endl;
  mwArray Ol,Or;
  int kfirst=even?0:1;
  int klast=(L-1)-(L-kfirst)%2; // last occupied by the loop
  for(int k=kfirst;k<klast;k+=2){
    getTwoBodyTermExponential(Ol,Or,delta,k);
    expH.setOp(k,new Operator(Ol),true);
    expH.setOp(k+1,new Operator(Or),true);
  }
  // Fill in the edges with identity operators
  if(kfirst==1){ // Identity in the first one
    int d1=d[0];
    mwArray ident=identityMatrix(d1);
    ident.reshape(Indices(d1,1,d1,1)); // just dx1xdx1
    expH.setOp(0,new Operator(ident),true);
  }
  if(klast<L-1){ // Identity at the end
    for(int k=klast+1;k<L;k++){
      int d1=d[k];
      mwArray ident=identityMatrix(d1);
      ident.reshape(Indices(d1,1,d1,1)); // just dx1xdx1
      expH.setOp(L-1,new Operator(ident),true);
    }
  }
}


void NNHamiltonian::getExponentialMPOeven(MPO& expHe,complex_t delta){
  getExponentialMPO(expHe,delta,true);
}

void NNHamiltonian::getExponentialMPOodd(MPO& expHo,complex_t delta){
  getExponentialMPO(expHo,delta,false);
}

void NNHamiltonian::getMPOforTerm(MPO& localTerm,int pos) const{
  if(pos>=L||d.size()<=pos+1){
    cout<<"Error: NNHamiltonian cannot getMPO for term "<<pos<<endl;
    exit(1);
  }
  localTerm.initLength(L);
  for(int k=0;k<L;k++){
    if(k<pos||k>pos+1){
      int dloc=d[pos];
      mwArray id=identityMatrix(dloc);
      id.reshape(Indices(dloc,1,dloc,1));
      localTerm.setOp(k,new Operator(id),true);
    }
    else if(k==pos){
      mwArray term=terms[pos];
      int d1=d[k];int d2=d[k+1];
      // Reshape according to the svd to be applied
      term.reshape(Indices(d1,d2,d1,d2));
      term.permute(Indices(1,3,2,4));
      term.reshape(Indices(d1*d1,d2*d2));
      mwArray Ol,Or;int nr=0;
      splitTwoBodyTerm(term,d1,d2,Ol,Or,nr);
      localTerm.setOp(pos,new Operator(Ol),true);
      localTerm.setOp(pos+1,new Operator(Or),true);      
    }
  }
}
