/**
   \file uMPS.h
   Basic implementation of a uniform MPS 
   
   \author Mari-Carmen Banuls
   \date 27/12/2019
*/
#include "uMPS.h"
#include "MPOMultiplier.h"
#include "DoubleOperator.h"

using namespace std;
using namespace shrt;

uMPS::uMPS(int d,int D):nA(1){
  W=new Site*[nA];
  Lambda=new mwArray*[nA];
  for(int k=0;k<nA;k++){
    W[k]=new Site(k,Indices(d,D,D));
    Lambda[k]=new mwArray(identityMatrix(D));
  }
  isCanonical_=false;
}

uMPS::uMPS(int nA_,int d,int D):nA(nA_){
  W=new Site*[nA];
  Lambda=new mwArray*[nA];
  for(int k=0;k<nA;k++){
    W[k]=new Site(k,Indices(d,D,D));
    Lambda[k]=new mwArray(identityMatrix(D));
  }
  isCanonical_=false;
}

uMPS::uMPS(int nA_,vector<int>& ds,vector<int>& Ds):nA(nA_){
  W=new Site*[nA];
  Lambda=new mwArray*[nA];
  for(int k=0;k<nA;k++){
    int Dl=Ds[k];
    int Dr=k<nA-1?Ds[k+1]:Ds[0];
    W[k]=new Site(k,Indices(ds[k],Dl,Dr));
    Lambda[k]=new mwArray(identityMatrix(Dl));
  }
  isCanonical_=false;
}

uMPS::uMPS(vector<mwArray>& tensors){
  nA=tensors.size();
  W=new Site*[nA];
  Lambda=new mwArray*[nA];
  for(int k=0;k<nA;k++){
    W[k]=new Site(k,tensors[k]);
    int Dl=W[k]->getDl();
    Lambda[k]=new mwArray(identityMatrix(Dl));
    if(k>=1){ // check left dimension with the previous right
      if(W[k]->getDl()!=W[k-1]->getDr()){
	cout<<"Error: incompatible dimensions in constructor of uMPS"<<endl;
	exit(1);
      }      
    }
    if(k==nA-1){ // check right dim with the initial left
      if(W[0]->getDl()!=W[k]->getDr()){
	cout<<"Error: incompatible dimensions in constructor of uMPS"<<endl;
	exit(1);
      }
    }
  }
  canonical();
}

uMPS::uMPS(const uMPS& state):nA(state.getN()){
  W=new Site*[nA];
  Lambda=new mwArray*[nA];
  for(int k=0;k<nA;k++){
    W[k]=new Site(k,state.getA(k));
    Lambda[k]=new mwArray(state.copyLambda(k));
  }
  isCanonical_=state.isCanonical();  
}

uMPS::~uMPS(){clear();}

void uMPS::clear(){
  nA=0;
  if(W!=0){
    for(int k=0;k<nA;k++){
      delete W[k];
      delete Lambda[k];
    }
    delete []W;
    delete []Lambda;
  }
  isCanonical_=false;
}

mwArray uMPS::getLambda(int pos) const{
  if(!isCanonical_){ //canonical(0);
    cout<<"ERROR: uMPS::getLambda invoked before imposing canonical form!"<<endl;
									    exit(1);
  }
  return *Lambda[pos];
}

double uMPS::getEntropy(int pos){
  mwArray lamb=getLambda(pos);
  double entr=0.;
  int Dl=lamb.getDimension(0);
  for(int k=0;k<Dl;k++){
    double valK=real(lamb.getElement(Indices(k,k)));
    if(abs(valK)>1E-12){
      entr+=-2*valK*valK*log2(valK);
    }
  }
  return entr;
}

/** Initialize with a given product state */
void uMPS::setProductState(ProductState state){
  for(int k=0;k<nA;k++)
    W[k]->setState((SingleState)state);
  canonical();
}

void uMPS::setRandomState(){
  for(int k=0;k<nA;k++)
    W[k]->setState(randmps);
  canonical();
}

void uMPS::replaceSite(int pos,const mwArray& A){
  if(pos<0||pos>=nA){
    cout<<"ERROR: Invalid index in uMPS::replaceSite"<<endl;
						       exit(1);
  }
  Indices dims=A.getDimensions();
  if(nA>1){//Check dimensions
    int posL=pos>0?pos-1:nA-1;
    int posR=pos<nA-1?pos+1:0;
    if(dims[0]!=W[posL]->getDr()||dims[2]!=W[posR]->getDr()){
      cout<<"ERROR: Invalid dimensions "<<dims<<" of the new tensor in uMPS::replaceSite"<<endl;
											   exit(1);
    }
  }
  W[pos]->setRotatedSite(A,Indices(1,2,3),false);
  Lambda[pos]->clear();
  cout<<"Tensor "<<pos<<" replaced by "<<A.getDimensions()<<endl;
  isCanonical_=false;
  canonical();
}

void uMPS::put(ostream& os) const{
  os<<"uMPS("<<nA<<"): ";
  for(int k=0;k<nA;k++){
    os<<"\tW["<<k<<"]="<<W[k]->getDimensions(); //<<*W[k];
    //os<<" Lambda["<<k<<"]="<<*Lambda[k];
    os<<endl;
  }
}

ostream& operator<<(ostream& os, const uMPS& mps){
  mps.put(os);
  return os;
}
  
  /** Read from a file */
void uMPS::importuMPS(const string& filename){
  clear();
  ifstream infile(filename,ios::binary);
  if(!infile.is_open()){
    cout<<"Error: unable to open file "<<filename<<" to read"<<endl;
    exit(2);
  }
  infile.read((char*)&nA,sizeof(int));
  if(!infile){ // sth happened when reading
    cout<<"Error trying to read uMPS from file "<<endl;
    exit(1);
  } 
  W=new Site*[nA];
  Lambda=new mwArray*[nA];
  for(int k=0;k<nA;k++){
    W[k]=new Site(k,mwArray(0));
    W[k]->load(infile);
    int Dl=W[k]->getDl();
    Lambda[k]=new mwArray(identityMatrix(Dl));
  }
  infile.close();
  isCanonical_=false;
  canonical();
}
  
  /** Write to a file */
void uMPS::exportuMPS(const string& filename) const{
  ofstream outfile(filename.data(),ios::binary);
  if(!outfile.is_open()){
    cout<<"Error: unable to open file "<<filename<<" to write uMPS"<<endl;
    exit(2);
  }
  else{
    outfile.write((char*)&nA,sizeof(int));
    for(int k=0;k<nA;k++){
      W[k]->save(outfile);
    }
    outfile.close();
  }
}

void uMPS::canonical(int Dcut,double tolSVD){
  //cout<<"uMPS::canonical(Dcut="<<Dcut<<")"<<endl;
  // For each cyclic permutation of the multipliers (tensors), need to
  // apply the gauge condition, with the multipliers in the right
  // order, and applying the changes to the left on the k-th tensor
  // and to the right on the (k-1)-th
  for(int k0=0;k0<nA;k0++){
    // Each tensor gives a TensorMultiplier
    vector<Multiplier*> E;
    vector<Multiplier*> El(nA,NULL);
    for(int k=0;k<nA;k++){
      mwArray W_=W[k]->getA();
      int dk=W[k]->getd();
      mwArray auxOp=identityMatrix(dk);auxOp.reshape(Indices(1,dk,1,dk)); // trick for TensorMultiplier
      E.push_back(new TensorMultiplier(permute(W_,Indices(2,1,3)),auxOp,permute(conjugate(W_),Indices(2,1,3))));
      El[nA-k-1]=new TensorMultiplier(permute(conjugate(W_),Indices(3,1,2)),permute(conjugate(auxOp),Indices(3,2,1,4)),permute(W_,Indices(3,1,2)));
    }
    applyCyclicCanonical(k0,E,El,Dcut,tolSVD);
    for(int k=0;k<nA;k++){
      delete E[k];delete El[k];
    }
    E.clear();El.clear();
  }
  isCanonical_=true;

  if(nA>1&&Dcut!=0){
    //cout<<"Repeating the gauge cond. after truncation, starting with uMPS "<<*this<<endl;
    isCanonical_=false;
    canonical(0,tolSVD);
  }
}

void uMPS::applyCyclicCanonical(int k0,vector<Multiplier*>& E,vector<Multiplier*>& El,int Dcut,double tolSVD){
  //  cout<<"uMPS::applyCyclicCanonical("<<k0<<")"<<endl;
  int D=W[k0]->getDl();
  if(Dcut==0) Dcut=D; // no truncation
  mwArray Xt,Yt,X,Y;
  if(!Lambda[k0]->isEmpty()){
    Xt=*Lambda[k0];
  }
  else{
    Xt=identityMatrix(D);
  }
  Yt=identityMatrix(D);
  complex_t etaL,etaR;
  getLargestEigenvalueCyclic(E,k0,Xt,etaR,false);
  //cout<<"Largest Right eigenvalue: "<<etaR<<", vec:"<<Xt.getDimensions()<<endl;
  getLargestEigenvalueCyclic(El,nA-1-k0,Yt,etaL,true); Yt.Hconjugate();  // 1 x xi*xi'*D*D'
  // cout<<"Found largest eigenvalues R:"<<etaR<<", L:"<<etaL<<", and <Yt|Xt>="<<Yt*Xt<<endl;
  //  cout<<"Yt="<<Yt<<endl;cout<<"Xt="<<Xt<<endl;
  Xt.reshape(Indices(D,D));
  // This guy should be Hermitian, but the global phase may be
  // arbitrary from the algorithm=> I find it looking at the diagonal,
  // and remove it.
  complex_t elX=getFirstNonZeroDiag(Xt);
  complex_t phase=conjugate(elX)/abs(elX);
  Xt=phase*Xt;
  if(!Xt.isHermitian(1E-10)){
    cout<<"% ERROR!??! Xt not Hermitian (should be Hermitian times a phase!): "<<endl;
    // exp(i 2* phi) X = X^dagger
    /// gt a value which is not zero 
    //putForMatlab(cout,Xt,"Xt");
  }
  Xt=.5*(Xt+Hconjugate(Xt));
  //cout<<" Xt="<<Xt<<endl;
  // cout<<"Empezando con Yt:"<<Yt.getDimensions()<<endl;
  Yt.reshape(Indices(D,D));
  Yt.permute(Indices(2,1));
  complex_t elY=getFirstNonZeroDiag(Yt);
  phase=conjugate(elY)/abs(elY);
  Yt.multiplyLeft(phase);
  if(!Yt.isHermitian(1E-10)){
    cout<<"% ERROR!??! Yt not Hermitian: "<<endl;
    //putForMatlab(cout,Yt,"Yt");
  }
  Yt=.5*(Yt+Hconjugate(Yt));
  complex_t trYX;
  {
    mwArray auxTr=Yt*Xt;
    trYX=auxTr.trace(); // orhtonormalize Yt*Xt
  }
  if(abs(imag(trYX))>1E-10){
    cout<<"ERROR!Why is the product of Y and X not a positive real number??"<<endl;
    exit(1);
  }
  Yt=1/sqrt(abs(trYX))*Yt;
  Xt=1/sqrt(abs(trYX))*Xt;
  //cout<<"Xt is "<<Xt<<endl;
  //cout<<"Yt is "<<Yt<<endl;
  vector<complex_t> Dx,Dy;
  mwArray Ux,Uy;
  wrapper::eig(Xt,Dx,Ux,true);
  //  cout<<"After eig(Xt) Dx="<<Dx<<endl;
  wrapper::eig(Yt,Dy,Uy,true);
  //cout<<"After eig(Yt) Dy="<<Dy<<endl;
  if(real(Dx[0])<-1E12){
    cout<<"X seems to be all negative=> changing sign"<<endl;
    X=Ux*sqrt(-ONE_c*diag(Dx))*Hconjugate(Ux);
  }
  else X=Ux*sqrt(diag(Dx))*Hconjugate(Ux);
  if(real(Dy[0])<-1E12){
    cout<<"Y seems to be all negative=> changing sign"<<endl;
    Y=Uy*sqrt(-ONE_c*diag(Dy))*Hconjugate(Uy);
  }
  else Y=Uy*sqrt(diag(Dy))*Hconjugate(Uy);
  //    cout<<"Created X:"<<X.getDimensions()<<" and Y:"<<Y.getDimensions()<<endl;
  mwArray U,S,V;
  wrapper::svd(Y*X,tolSVD,Dcut,U,S,V);
  *Lambda[k0]=(1./sqrt(real((S*S).trace())))*S;
  mwArray valsS;S.getDiagonal(valsS);
  // cout<<"After svd of YX, S="<<valsS<<endl;
  //cout<<"Lambda["<<k0<<"]="<<*Lambda[k0]<<endl;
  int nrS=0;
  mwArray auxS=invertDiag(S,nrS,1E-8);
  //cout<<"nrS="<<nrS<<endl;//cout<<"invS="<<auxS<<endl;
  mwArray W_l; // a copy!
  if(nA>1){
    if(k0==0) W_l=W[nA-1]->getA();
    else W_l=W[k0-1]->getA();
  }
  else{
    W_l=W[0]->getA();
  }
  W_l=(ONE_c/sqrt(etaR))*W_l;
  int d_l=W_l.getDimension(0);
  int Dl_l=W_l.getDimension(1);
  int Dr_l=W_l.getDimension(2);
  W_l.reshape(Indices(d_l*Dl_l,Dr_l));
  W_l.multiplyRight(X*Hconjugate(V)*auxS);//invertDiag(S,nrS));
  //W.multiplyRight(X*Hconjugate(V)*invertDiag(S,nrS));
  W_l.reshape(Indices(d_l,Dl_l,Dcut));
  // Now, reset in the right position. Use setRotated to avoid check of nr of components
  if(k0==0) W[nA-1]->setRotatedSite(W_l,Indices(1,2,3),false);
  else W[k0-1]->setRotatedSite(W_l,Indices(1,2,3),false);
  //cout<<"In applyCyclicCanonical("<<k0<<"), left W is k="<<(k0==0?nA-1:k0-1)<<", now "<<(k0==0?W[nA-1]->getDimensions():W[k0-1]->getDimensions())<<endl;
  // Transform the k-th tensor on the left
  mwArray W_r=W[k0]->getA();
  int d_r=W_r.getDimension(0);
  int Dl_r=W_r.getDimension(1);
  int Dr_r=W_r.getDimension(2);
  W_r.permute(Indices(2,1,3));
  W_r.reshape(Indices(Dl_r,d_r*Dr_r));
  W_r.multiplyLeft(Hconjugate(U)*Y);
  W_r.reshape(Indices(Dcut,d_r,Dr_r));
  W_r.permute(Indices(2,1,3));
  W[k0]->setRotatedSite(W_r,Indices(1,2,3),false);
  //cout<<"In applyCyclicCanonical("<<k0<<"), W["<<k0<<"]->"<<W[k0]->getDimensions()<<endl;
}
 
  /** Apply a MPO and impose canonical form again */
void uMPS::applyMPO(std::vector<const mwArray*>& midU,double tolSVD,int Dcut){
  //cout<<"uMPS::applyMPO()"<<endl;
  applyMPOnoCanonical(midU);
  canonical(Dcut,tolSVD);
}


void uMPS::applyMPOnoCanonical(std::vector<const mwArray*>& midU){
  //cout<<"uMPS::applyMPOnoCAnonical()"<<endl;
  for(int k=0;k<nA;k++){
    //cout<<"Applying operator on tensor "<<k<<endl;
    mwArray W_=W[k]->getA(); // copy
    int d=W_.getDimension(0);
    int Dl=W_.getDimension(1);
    int Dr=W_.getDimension(2);
    int xi_l=midU[k]->getDimension(1); // d xi_l d xi_r
    int xi_r=midU[k]->getDimension(3);
    W_.reshape(Indices(d,Dl*Dr));
    W_.multiplyLeft(reshape(permute(*midU[k],Indices(1,2,4,3)),Indices(d*xi_l*xi_r,d))); //d*xi*xiXD*D
    W_.reshape(Indices(d,xi_l,xi_r,Dl,Dr));
    W_.permute(Indices(1,2,4,3,5));
    W_.reshape(Indices(d,xi_l*Dl,xi_r*Dr));
    W[k]->setRotatedSite(W_,Indices(1,2,3),false);
    //    cout<<"After applying MPO, tensor at site "<<k<<" is "<<W[k]->getDimensions()<<endl;
    // Trick for the later application of canonical
    Lambda[k]->reshape(Indices(Dl*Dl,1));
    mwArray auxM=identityMatrix(xi_l);
    //mwArray auxM(Indices(xi_l,xi_l));auxM.fillWithOne();auxM=(.5)*auxM;
    Lambda[k]->multiplyRight(reshape(auxM,Indices(1,xi_l*xi_l)));
    Lambda[k]->reshape(Indices(Dl,Dl,xi_l,xi_l));
    Lambda[k]->permute(Indices(3,1,4,2));
    Lambda[k]->reshape(Indices(xi_l*Dl,xi_l*Dl));

  }
  isCanonical_=false;
}


void uMPS::getLargestEigenvalueCyclic(vector<Multiplier*>& M,int k0,mwArray& V,complex_t& lambda,bool reverse) const{
  //cout<<"uMPS::getLargestEigenvalueCyclic("<<k0<<") with init V:"<<V.getDimensions()<<endl;
  vector<int> ind;
  if(!reverse){
    for(int k=k0;k<nA;k++) ind.push_back(k);
    for(int k=0;k<k0;k++) ind.push_back(k);
  }
  else{
    for(int k=k0+1;k<nA;k++) ind.push_back(k);
    for(int k=0;k<=k0;k++) ind.push_back(k);
  }
  int D=M[ind[nA-1]]->getSize();
  //cout<<"getLargestEigenvalue M["<<ind[nA-1]<<"] size "<<D<<endl;
  bool conv=0;
  if(V.isEmpty()){
    V=mwArray(Indices(D,1));V.fillRandom();
    V=(1./norm(V))*V;
  }
  else{
    V.reshape(Indices(D,1));
  }
  int numIt=1;
  while(!conv&&numIt<10000){
    //cout<<"Iter: "<<numIt<<" V.nrComp "<<V.getNrComponents()<<endl;
    mwArray newV(V);
    for(int k=nA-1;k>=0;k--){
      mwArray aux(newV);
      M[ind[k]]->product(aux,newV);
    }
    mwArray aux=Hconjugate(newV)*V;
    complex_t new_lambda=aux.getElement(0);
    //if(abs(new_lambda-lambda)/abs(lambda)<1E-13) conv=true;    
    if(abs(new_lambda/lambda-ONE_c)<1E-13) conv=true;    
    V=(1./norm(newV))*newV;lambda=new_lambda;
    numIt++;
  }
  if(!conv) cout<<"**** getLargestEigenvalue not converged!"<<endl;
  //else cout<<"numIt="<<numIt<<endl;
  // cout<<"It nr "<<numIt<<". Found largest eigenvalue to be "<<lambda<<endl;
  // mwArray aux;M.product(V,aux);
  // cout<<"CHECK!! M*V-lambda*V"<<norm(aux-lambda*V)<<endl;
}

complex_t uMPS::getFirstNonZeroDiag(const mwArray& X) const {
  //cout<<"uMPS::getFirstNonZeroDiag()"<<endl;
  bool done=false;
  int D=X.getDimension(0);
  int cnt=0;
  //cout<<"getFirstNonZeroDiag of X:"<<X.getDimensions()<<endl;
  while(!done&&cnt<D){
    complex_t res=X.getElement(Indices(cnt,cnt));
    //cout<<"Looked al element ("<<cnt<<")="<<res<<endl;
    if(abs(res)>1E-3) return res;
    //cout<<"Looked al element ("<<cnt<<") was not big enough"<<endl;
    cnt++;
  }
  if(!done){
    cout<<"ERROR! All diagonal elements of X are too small!"<<endl;
    putForMatlab(cout,X,"X");
    exit(1);
  }
  exit(1); // never reached
}

complex_t uMPS::computeExpectationValue(int nOp,const vector<int>& pos,
					const vector<mwArray*>& op) {
  //  cout<<"uMPS::computeExpectationValue() pos="<<pos<<endl;
  if(!isCanonical_) canonical();
  if(nA==1){
    return computeExpectationValueSingle(nOp,pos,op);
  }
  if(nOp>nA){
    cout<<"Error: uMPS::computeExpectationValue can only handle up to "
	<<nA<<" operators"<<endl;
    exit(1);
  }
  // which positions have an operator, and which one
  vector<mwArray*> hasOp(nA,NULL);
  for(int k=0;k<nOp;k++){
    hasOp[pos[k]]=op[k];
  }
  // the block will be taken to start at pos[0]
  int k0=pos[0];
  vector<int> ind;
  for(int k=k0;k<nA;k++) ind.push_back(k);
  for(int k=0;k<k0;k++) ind.push_back(k);
  // on the right, the corresponding Lambda squared
  mwArray V(*Lambda[k0]);V=V*V;
  for(int k=nA-1;k>=0;k--){
    mwArray W_=W[ind[k]]->getA();
    int dk=W_.getDimension(0);
    mwArray auxOp=identityMatrix(dk);
    if(hasOp[ind[k]]!=0){
      auxOp=*hasOp[ind[k]]; // d x d      
    }
    auxOp.reshape(Indices(1,dk,1,dk)); // trick for TensorMultiplier
    TensorMultiplier E(permute(W_,Indices(2,1,3)),auxOp,permute(conjugate(W_),Indices(2,1,3)));
    mwArray aux(V);
    E.product(aux,V);
  }
  // At the very end, contract with identity
  mwArray idL=identityMatrix(W[k0]->getDl());
  idL.reshape(Indices(-1,1));idL.Hconjugate();
  return (idL*V).getElement(0);
}

void uMPS::getSingleSiteRDM(int pos,mwArray& rdm) {
  if(!isCanonical_) canonical();
  vector<int> ind;
  for(int k=pos;k<nA;k++) ind.push_back(k);
  for(int k=0;k<pos;k++) ind.push_back(k);
  //  cout<<"uMPS::getSingleSiteRDM("<<pos<<"); ind="<<ind<<endl;
  // on the right, the corresponding Lambda squared
  mwArray V(*Lambda[pos]); 
  V=V*V;
  for(int k=nA-1;k>=1;k--){
    mwArray W_=W[ind[k]]->getA();
    int dk=W_.getDimension(0);
    mwArray auxOp=identityMatrix(dk);
    auxOp.reshape(Indices(1,dk,1,dk)); // trick for TensorMultiplier
    TensorMultiplier E(permute(W_,Indices(2,1,3)),auxOp,permute(conjugate(W_),Indices(2,1,3)));
    mwArray aux(V);
    E.product(aux,V);
    //cout<<"After multiplier for "<<ind[k]<<", V="<<V<<endl;
  }
  mwArray Wpos=W[pos]->getA();
  rdm=Wpos;
  int d=Wpos.getDimension(0);
  int D=Wpos.getDimension(1);
  int Dr=Wpos.getDimension(2);
  V.reshape(Indices(Dr,Dr));
  rdm.reshape(Indices(d*D,Dr));
  rdm.multiplyRight(V);
  rdm.reshape(Indices(d,D*Dr));
  rdm.multiplyRight(Hconjugate(reshape(Wpos,Indices(d,D*Dr)))); // d x d
  rdm.transpose();
  //  rdm.reshape(Indices(d*d,1));
}


complex_t uMPS::computeExpectationValueSingle(int nOp,const vector<int>& pos,
					      const vector<mwArray*>& op) const{
  //  cout<<"uMPS::computeExpectationValueSingle("<<nOp<<") pos="<<pos<<endl;
  if(!isCanonical_){
    cout<<"ERROR: uMPS::computeExpectationValueSingle should only be called after imposing canonical form"<<endl;
    exit(1);
  }
  // in this case, the tensor is always the same one
  mwArray W_=W[0]->getA();
  int d=W_.getDimension(0);
  // which positions have an operator, and which one
  mwArray auxOp=identityMatrix(d);
  vector<mwArray*> hasOp(pos.size(),&auxOp);
  for(int k=0;k<nOp;k++){
    hasOp[pos[k]]=op[k];
  }
  //cout<<"hasOp:"<<hasOp<<endl;
  // on the right, the corresponding Lambda squared
  mwArray V(*Lambda[0]);V=V*V;
  for(int k=nOp-1;k>=0;k--){
    if(hasOp[k]!=0){
      auxOp=*hasOp[k]; // d x d
      //cout<<"Operator on position "<<k<<" is:"<<auxOp<<endl; 
    }
    auxOp.reshape(Indices(1,d,1,d)); // trick for TensorMultiplier
    TensorMultiplier E(permute(W_,Indices(2,1,3)),auxOp,permute(conjugate(W_),Indices(2,1,3)));
    mwArray aux(V);
    E.product(aux,V);
  }
  // At the very end, contract with identity
  mwArray idL=identityMatrix(W[0]->getDl());
  idL.reshape(Indices(-1,1));idL.Hconjugate();
  return (idL*V).getElement(0);
}

void uMPS::getPurity(complex_t& eta){
  // This ony works if this uMPS is a purification, i.e. d=d0*d0
  if(!isCanonical_) canonical(); // ensure canonical form
  // Construct the MPO(s) for the transfer operator in rho^2

  // Each tensor gives a MPOMultiplier
  vector<Multiplier*> E;
  for(int k=0;k<nA;k++){
    mwArray W_=W[k]->getA(); // dkxDlxDr
    int dk=W[k]->getd();
    int dk0=floor(sqrt(dk));
    if(dk0*dk0!=dk){
      cout<<"ERROR: The "<<k<<"-th tensor in this uMPO does not correspond to a purification"
	  <<", with d="<<dk<<endl;
      exit(1);
    }
    int Dl=W_.getDimension(1);
    int Dr=W_.getDimension(2);
    MPOMultiplier* Ek=new MPOMultiplier(4);
    Operator* A0=new Operator(reshape(W_,Indices(dk,Dl,Dr,1)),Indices(2,4,3,1),0);
    Ek->setOp(0,A0,true); // all mine
    Operator* A3=new Operator(reshape(W_,Indices(dk,Dl,Dr,1)),Indices(2,1,3,4),true);
    Ek->setOp(3,A3,true); // all mine
    // the middle ones are doubled with identities
    mwArray idPhys=identityMatrix(dk0);idPhys.reshape(Indices(1,dk0,1,dk0));
    W_.reshape(Indices(dk0,dk0,Dl,Dr));
    Ek->setOp(1,new DoubleOperator(permute(conjugate(W_),Indices(3,1,4,2)),idPhys),true);
    Ek->setOp(2,new DoubleOperator(permute(W_,Indices(3,2,4,1)),idPhys),true);
    E.push_back(Ek);
  }
  mwArray Rvec; // for now, not returning
  getLargestEigenvalueCyclic(E,0,Rvec,eta,false);

  for(int k=0;k<nA;k++){
    delete E[k];
  }
  E.clear();
}


void uMPS::getPurityBlock(complex_t& eta,int L){
  // This ony works if this uMPS is a purification, i.e. d=d0*d0
  if(!isCanonical_) canonical(); // ensure canonical form
  if(L%nA!=0){
    cout<<"ERROR: Cannot yet compute rhoL^2 for a block if the length is not multiple of the period"<<endl;
    exit(1);
  }
  int nRep=L/nA; // nr of times the multiplier will be applied
  // Construct the MPO(s) for the transfer operator in rho^2
  // Each tensor gives a MPOMultiplier
  vector<Multiplier*> E;
  for(int k=0;k<nA;k++){
    mwArray W_=W[k]->getA(); // dkxDlxDr
    int dk=W[k]->getd();
    int dk0=floor(sqrt(dk));
    if(dk0*dk0!=dk){
      cout<<"ERROR: The "<<k<<"-th tensor in this uMPO does not correspond to a purification"
	  <<", with d="<<dk<<endl;
      exit(1);
    }
    int Dl=W_.getDimension(1);
    int Dr=W_.getDimension(2);
    MPOMultiplier* Ek=new MPOMultiplier(4);
    Operator* A0=new Operator(reshape(W_,Indices(dk,Dl,Dr,1)),Indices(2,4,3,1),0);
    Ek->setOp(0,A0,true); // all mine
    Operator* A3=new Operator(reshape(W_,Indices(dk,Dl,Dr,1)),Indices(2,1,3,4),true);
    Ek->setOp(3,A3,true); // all mine
    // the middle ones are doubled with identities
    mwArray idPhys=identityMatrix(dk0);idPhys.reshape(Indices(1,dk0,1,dk0));
    W_.reshape(Indices(dk0,dk0,Dl,Dr));
    Ek->setOp(1,new DoubleOperator(permute(conjugate(W_),Indices(3,1,4,2)),idPhys),true);
    Ek->setOp(2,new DoubleOperator(permute(W_,Indices(3,2,4,1)),idPhys),true);
    E.push_back(Ek);
  }
  // Now I need the boundary vectors
  mwArray V(*Lambda[0]); // I'll take std block
  V=V*V; // Dr x Dr
  int Dr=V.getDimension(0);
  V.reshape(Indices(Dr*Dr,1));
  V.multiplyRight(permute(V,Indices(2,1)));
  V.reshape(Indices(Dr,Dr,Dr,Dr));
  V.permute(Indices(3,2,1,4));
  V.reshape(Indices(Dr*Dr*Dr*Dr,1));

  mwArray Vl=identityMatrix(Dr*Dr);
  Vl.reshape(Indices(Dr,Dr,Dr,Dr));
  Vl.permute(Indices(1,2,4,3));
  Vl.reshape(Indices(1,Dr*Dr*Dr*Dr));

  for(int kR=0;kR<nRep;kR++){
    for(int k=nA-1;k>=0;k--){
      mwArray aux(V);
      E[k]->product(aux,V);
    }
  }  
  mwArray aux=Vl*V;
  eta=aux.getElement(0);
}

