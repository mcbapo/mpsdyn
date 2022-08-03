#include "JoinedOperator.h"
#include "Indices.h"

using namespace std;
using namespace shrt;

/** */
JoinedOperator::JoinedOperator(int n_,const Operator* ops[]):
  n(n_),data(){
  int du=ops[0]->getd();int dd=ops[0]->getdorig();
  int Dl=1;int Dr=1;
  // cout<<"Creating JoinedOperator with "<<n<<" terms: ";
  // for(int k=0;k<n;k++)cout<<ops[k]->getDimensions()<<", ";
  // cout<<endl;
  for(int k=0;k<n;k++){   // Saves a copy of each data matrix
    if(k>0&&(ops[k]->getd()!=dd)){
      cout<<"Error: incompatible dimensions in JoinedOperator, component "
	  <<k<<": ops["<<k-1<<"]->dd="<<dd<<" and ops["<<k<<"]->du="
	  <<ops[k]->getd()<<endl;
      //cout<<ops[k-1]->getDimensions()<<endl;
      //cout<<ops[k]->getDimensions()<<endl;
      exit(212);
    }
    dd=ops[k]->getdorig();
    Dl*=ops[k]->getDl();Dr*=ops[k]->getDr();
    data.push_back(ops[k]->getFullData());
  }
  dims=Indices(du,Dl,dd,Dr);
}

JoinedOperator::~JoinedOperator(){
  dims.clear();
  data.clear();
}

JoinedOperator::JoinedOperator(const JoinedOperator& op):
  n(op.n),data(op.data){
  dims=op.getDimensions();
}

void JoinedOperator::conjugateOp(){
  for(int k=0;k<n;k++){
    data[k].conjugate();
  }
}

void JoinedOperator::permuteOp(Indices perm){
  if(perm[0]!=dims[0]||perm[2]!=dims[2]){
    cout<<"ERROR: JoinedOperator only allows permutation of side indices!"
	<<endl;
    exit(212);
  }
  for(int k=0;k<n;k++){
    data[k].permute(perm);
  }
  // permute the dimensions too
  Indices oldDims=dims;
  for(int k=0;k<4;k++){
    dims[k]=oldDims[perm[k]];
  }
}

void JoinedOperator::contractRight(const JoinedOperator* rightOp){
  if(rightOp->getNumber()!=n){
    cout<<"ERROR: Cannot use contractRight on JoinedOperators of different"
	<<" number of components"<<endl;
    exit(212);
  }
  for(int k=0;k<n;k++){
    Indices dimsLoc=data[k].getDimensions();
    mwArray dataR=rightOp->getDataN(k);
    Indices dimsR=dataR.getDimensions();
    dataR.permute(Indices(2,1,3,4));
    dataR.reshape(Indices(dimsR[1],-1));
    data[k].reshape(Indices(-1,dimsLoc[3]));
    data[k].multiplyRight(dataR);
    data[k].reshape(Indices(dimsLoc[0],dimsLoc[1],dimsLoc[2],dimsR[0],
			    dimsR[2],dimsR[3]));
    data[k].permute(Indices(1,4,2,3,5,6));
    data[k].reshape(Indices(dimsLoc[0]*dimsR[0],dimsLoc[1],
			    dimsLoc[2]*dimsR[2],dimsR[3]));
  }
  // set the global dimensions
  dims[0]=data[0].getDimension(0); //up
  dims[3]=rightOp->getDr(); // right (left won't change)
  dims[2]=data[n-1].getDimension(2); // down
}


JoinedOperator& JoinedOperator::operator=(const JoinedOperator& op){
  if(this!=&op){
    this->Operator::operator=(op);
    n=op.n;
    data=op.data;
    dims=op.getDimensions();
  }
  return *this;
}

void JoinedOperator::setData(int nr,const mwArray* opers[]){
  n=nr;
  data.clear();dims.clear();
  int du=opers[0]->getDimension(0);int dd=opers[0]->getDimension(2);
  int Dl=1;int Dr=1;
  for(int k=0;k<nr;k++){
    if(k>0&&(opers[k]->getDimension(0)!=dd)){
      cout<<"Error: incompatible dimensions in JoinedOperator, component "<<k<<endl;
      exit(212);
    }
    dd=opers[k]->getDimension(2);
    Dl*=opers[k]->getDimension(1);Dr*=opers[k]->getDimension(3);
    data.push_back(*opers[k]);
  }
  dims=Indices(du,Dl,dd,Dr);
}

void JoinedOperator::savetext(ofstream& outfile) const{
  outfile<<n<<" ";
  for(int k=0;k<n;k++){
    data[k].savetext(outfile);
  }
}

void JoinedOperator::loadtext(ifstream& infile){
  infile>>n;
  data.clear();data=vector<mwArray>(n);
  int Dl=1;int Dr=1;
  for(int k=0;k<n;k++){
    data[k].loadtext(infile);
    Dl*=data[k].getDimension(1);Dr*=data[k].getDimension(3);
 }
  dims=Indices(data[0].getDimension(0),Dl,data[n-1].getDimension(2),Dr);
}


/** ********************* CONTRACTIONS  ************************ */

void JoinedOperator::contractL(mwArray& result,const mwArray& termL,
			       const Site& ket,const Site& bra,
			       bool dagger) const{
  contractleftjop(result,termL,ket,bra,dagger);
}

void JoinedOperator::contractR(mwArray& result,const mwArray& termR,
			       const Site& ket,const Site& bra,
			       bool dagger) const{
  contractrightjop(result,termR,ket,bra,dagger);
}

void JoinedOperator::contractMket(mwArray& result,const mwArray& termL,
				  const Site& ket,const mwArray& termR,
				  bool dagger) const{
  contractMketjop(result,termL,ket,termR,dagger);
}

void JoinedOperator::contractMbra(mwArray& result,const mwArray& termL,
				  const Site& bra,const mwArray& termR) const{
  contractMbrajop(result,termL,bra,termR,false);
}

void JoinedOperator::contractN(mwArray& result,const mwArray& termL,
			       const mwArray& termR,bool dagger) const {
  contractNjop(result,termL,termR,dagger);
}


void JoinedOperator::put(std::ostream& os) const{
  os<<"JoinedOperator["<<n<<"]: ";
  for(int k=0;k<n;k++){
    os<<"op["<<k<<"]="<<data[k]<<endl;
  }
}


// ostream& operator<<(ostream& os,const JoinedOperator& oper){
//   os<<"JoinedOperator["<<oper.n<<"]: ";
//   for(int k=0;k<oper.n;k++){
//     os<<"op["<<k<<"]="<<oper.data[k]<<endl;
//   }
//   return os;
// }

const mwArray JoinedOperator::getFullData() const{
  int dd=data[0].getDimension(0);
  int Dl=data[0].getDimension(1);
  int ddp=data[0].getDimension(2);
  int Dr=data[0].getDimension(3);
  Indices curdims(dd,Dl,Dr);int prod=dd*Dl*Dr;
  mwArray tmp=permute(data[0],Indices(1,2,4,3));
  tmp.reshape(Indices(prod,ddp));
  int dimL=Dl;int dimR=Dr; //dims to l and r so far
  int ddn(dd),Dln(Dl),ddpn(ddp),Drn(Dr);
  for(int k=1;k<n;k++){
    ddn=data[k].getDimension(0);
    Dln=data[k].getDimension(1);
    ddpn=data[k].getDimension(2);
    Drn=data[k].getDimension(3);
#ifdef CHECKDIMS
    if(ddn!=ddp){
      cout<<"Error: incompatible dimensions in JoinedOperator::getFullData"
	  <<", on term "<<k<<endl;
      exit(2);
    }
    else ddp=ddpn;
#endif
    mwArray aux=permute(data[k],Indices(1,2,4,3));
    aux.reshape(Indices(ddn,Dln*Drn*ddpn));
    tmp.multiplyRight(aux);
    curdims=Indices(curdims,Indices(Dln,Drn));prod=prod*Dln*Drn;
    tmp.reshape(Indices(prod,ddpn));
    dimL=dimL*Dln;dimR=dimR*Drn;
  }
  tmp.reshape(Indices(curdims,ddpn)); // dd1,Dl1,Dr1,Dl2,Dr2,...ddpn
  vector<int> permutation(1,1);
  int last=curdims.size()+1;
  for(int l=2;l<=last;l=l+2) permutation.push_back(l);
  for(int l=3;l<=last-1;l=l+2) permutation.push_back(l);
  tmp.permute(permutation);
  tmp.reshape(Indices(dd,dimL,ddpn,dimR));
  return tmp; //reshape(tmp,Indices(dd,dimL,ddpn,dimR));
}

void JoinedOperator::contractleftjop(mwArray& result,
				     const mwArray& termL,const Site& ket,
				     const Site& bra,bool dagger) const{
  int d=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int betab=termL.getDimension(0);
  int alpha=termL.getDimension(1);  
  int betak=termL.getDimension(2);
  int du=dagger?data[0].getDimension(2):data[0].getDimension(0);
  int al1=data[0].getDimension(1);
  int dd=dagger?data[0].getDimension(0):data[0].getDimension(2);
  int al2=data[0].getDimension(3);
  int dun=dagger?data[n-1].getDimension(2):data[n-1].getDimension(0);
  int al1n=data[n-1].getDimension(1);
  int ddn=dagger?data[n-1].getDimension(0):data[n-1].getDimension(2);
  int al2n=data[n-1].getDimension(3);
#ifdef CHECKDIMS
  bool wrong=0;
  int ddlast=dd;
  int alphaL=al1;int alphaR=al2;
  for(int k=1;k<n;k++){  // k=2:n      
    int duk=dagger?data[k].getDimension(2):data[k].getDimension(0);
    int al1k=data[k].getDimension(1);
    int ddk=dagger?data[k].getDimension(0):data[k].getDimension(2);
    int al2k=data[k].getDimension(3);
    if(ddlast!=duk) wrong=1;
    alphaL=alphaL*al1k;alphaR=alphaR*al2k;      
    ddlast=ddk;
  }
  // After the last one, check left and right dimensions
  if((d!=du)||(dp!=ddn)||(Dlp!=betak)||(Dl!=betab)||
     (alphaL!=alpha))
    wrong=1;  
  if(wrong){
    cout<<"Wrong dimensions in JoinedOperator::contractleftjop:"<<endl;
    cout<<"Constraints one by one: d("<<d<<")=?du("<<du<<"), "
	<<"dp("<<dp<<")=?ddn("<<ddn<<"), "
	<<"Dlp("<<Dlp<<")=?betak("<<betak<<"), "
	<<"Dl("<<Dl<<")=?betab("<<betab<<"), "
	<<"alphaL("<<alphaL<<")=?alpha("<<alpha
	<<")"<<endl;
    cout<<" bra"<<bra.getDimensions()<<", ket"<<ket.getDimensions()
	<<", termL"<<termL.getDimensions()
	<<", n="<<n;
    for(int k=0;k<n;k++)
      cout<<", oper["<<k<<"]="<<data[k].getDimensions();
    cout<<endl;
    exit(212);
  }
#endif
  // First step: contract Aket with tmpL -> betab*alpha, dp*Drp
  mwArray aux=ket.getA();
  aux.permute(Indices(2,1,3));aux.reshape(Indices(Dlp,dp*Drp));
  mwArray tmpres(termL);tmpres.reshape(Indices(betab*alpha,betak));
  tmpres.multiplyRight(aux);
  tmpres.reshape(Indices(betab,alpha*dp*Drp));
  aux=bra.getA();aux.conjugate();
  aux.permute(Indices(3,1,2));
  aux.reshape(Indices(Dr*d,Dl));
  tmpres.multiplyLeft(aux);
  //  mwArray tmpres=reshape(termL,Indices(betab*alpha,betak))*
  //reshape(permute(ket.getA(),Indices(2,1,3)),Indices(Dlp,dp*Drp));
  // Now contract Abra -> Dr*d, alpha*dp*Drp
  //tmpres=reshape(permute(conjugate(bra.getA()),Indices(3,1,2)),
  //		 Indices(Dr*d,Dl))*
  // reshape(tmpres,Indices(betab,alpha*dp*Drp));
  // Now contract oper term by term
  // first reshape and sort indices
  tmpres.reshape(Indices(Dr*d*alpha*dp,Drp));
  tmpres.permute(Indices(2,1));  //
  //  tmpres=permute(reshape(tmpres,Indices(Dr*d*alpha*dp,Drp)),
  //		 Indices(2,1));  //
  int alphan=1; // new alpha
  //  mwArray aux;
  for(int k=n-1;k>0;k--){   //n:-1:2
    if(dagger)
      aux=permute(conjugate(data[k]),Indices(3,2,1,4));      
    else
      aux=data[k];
    dun=aux.getDimension(0);al1n=aux.getDimension(1);
    ddn=aux.getDimension(2);al2n=aux.getDimension(3);
    tmpres.reshape(Indices(-1,al1n*ddn));
    aux.permute(Indices(2,3,1,4));
    aux.reshape(Indices(al1n*ddn,dun*al2n));
    tmpres.multiplyRight(aux);
    tmpres.reshape(Indices(-1,al2n));
    tmpres.permute(Indices(2,1));
    // tmpres=reshape(tmpres,Indices(-1,al1n*ddn))*
    //  reshape(permute(aux,Indices(2,3,1,4)),Indices(al1n*ddn,dun*al2n));
      //tmpres=permute(reshape(tmpres,Indices(-1,al2n)),Indices(2,1));
    alphan=alphan*al2n; // new alpha
  }
  // The last contraction, with data[0]
  if(dagger)
    aux=permute(conjugate(data[0]),Indices(3,2,1,4));    
  else
    aux=data[0];
  //[du,al1,dd,al2]=size(oper{1});
  aux.reshape(Indices(du*al1*dd,al2));
  tmpres.reshape(Indices(-1,d*al1*dd));
  tmpres.multiplyRight(aux);
  //tmpres=reshape(tmpres,Indices(-1,d*al1*dd))*
  //reshape(aux,Indices(du*al1*dd,al2));
  alphan=alphan*al2; // new alpha
  tmpres.reshape(Indices(-1,Dr*al2));
  tmpres.permute(Indices(2,1));
  tmpres.reshape(Indices(Dr,alphan,Drp));
  result=tmpres;
  //  result=reshape(permute(reshape(tmpres,Indices(-1,Dr*al2)),Indices(2,1)),
  //		 Indices(Dr,alphan,Drp));
}

void JoinedOperator::contractrightjop(mwArray& result,const mwArray& termR,
				      const Site& ket,
				      const Site& bra,bool dagger) const{
  int d=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int betab=termR.getDimension(0);
  int alpha=termR.getDimension(1);  
  int betak=termR.getDimension(2);

  int d1=dagger?data[0].getDimension(2):data[0].getDimension(0);
  int al1=data[0].getDimension(1);
  int d2=dagger?data[0].getDimension(0):data[0].getDimension(2);
  int al2=data[0].getDimension(3);
  int d1n=dagger?data[n-1].getDimension(2):data[n-1].getDimension(0);
  int al1n=data[n-1].getDimension(1);
  int d2n=dagger?data[n-1].getDimension(0):data[n-1].getDimension(2);
  int al2n=data[n-1].getDimension(3);
#ifdef CHECKDIMS
  bool wrong=0;
  int alphaL=al1;int alphaR=al2;
  int dd=d2;
  for(int k=1;k<n;k++){  // k=2:n      
    int duk=dagger?data[k].getDimension(2):data[k].getDimension(0);
    int al1k=data[k].getDimension(1);
    int ddk=dagger?data[k].getDimension(0):data[k].getDimension(2);
    int al2k=data[k].getDimension(3);
    if(dd!=duk) wrong=1;
    alphaL=alphaL*al1k;alphaR=alphaR*al2k;      
    dd=ddk;
  }
  // After the last one, check left and right dimensions
  if((d!=d1)||(dp!=d2n)||(Drp!=betak)||(Dr!=betab)||(alphaR!=alpha))
    wrong=1;
  if(wrong){
    cout<<"Wrong dimensions in JoinedOperator::contractrightjop:"<<endl;
    cout<<"Constraints one by one: d("<<d<<")=?d1("<<d1<<"), "
	<<"dp("<<dp<<")=?d2n("<<d2n<<"), "
	<<"Drp("<<Drp<<")=?betak("<<betak<<"), "
	<<"Dr("<<Dr<<")=?betab("<<betab<<"), "
	<<"alphaR("<<alphaR<<")=?alpha("<<alpha
	<<")"<<endl;
    cout<<" bra="<<bra<<", ket="<<ket<<", termR="<<termR
	<<", n="<<n;
    for(int k=0;k<n;k++)
      cout<<", oper["<<k<<"]="<<data[k];
    cout<<endl;
    exit(212);
  }
#else 
  int dd;
#endif
  // First step, contract Aket with tmpR
  mwArray tmpres(termR);
  tmpres.reshape(Indices(betab*alpha,betak));
  mwArray aux(ket.getA());
  aux.permute(Indices(3,1,2));
  aux.reshape(Indices(Drp,dp*Dlp));
  tmpres.multiplyRight(aux);
  tmpres.reshape(Indices(betab,alpha*dp*Dlp));
  aux=bra.getA();aux.conjugate();
  aux.permute(Indices(2,1,3));
  aux.reshape(Indices(Dl*d,Dr));
  tmpres.multiplyLeft(aux);
  //  mwArray tmpres=reshape(termR,Indices(betab*alpha,betak))*
  //  reshape(permute(ket.getA(),Indices(3,1,2)),Indices(Drp,dp*Dlp));
  // Then Abra
  //  tmpres=reshape(permute(conjugate(bra.getA()),Indices(2,1,3)),
  //		 Indices(Dl*d,Dr))*
  //  reshape(tmpres,Indices(betab,alpha*dp*Dlp));
  // Now oper term by term
  // first reshape and sort indices  -> Dl to the bottom
  tmpres.reshape(Indices(Dl,d*alpha*dp*Dlp));
  tmpres.permute(Indices(2,1)); 
  //  tmpres=permute(reshape(tmpres,Indices(Dl,d*alpha*dp*Dlp)),Indices(2,1)); 
  int alphan=1; // new alpha
  //  mwArray aux;
  int du;
  for(int k=0;k<n-1;k++){ //k=1:n-1
    if(dagger)
      aux=permute(conjugate(data[k]),Indices(3,2,1,4));      
    else
      aux=data[k];
    du=aux.getDimension(0);int al1=aux.getDimension(1);
    dd=aux.getDimension(2);int al2=aux.getDimension(3);
    aux.permute(Indices(2,3,1,4));
    aux.reshape(Indices(al1*dd,du*al2));
    tmpres.reshape(Indices(du*al2,-1));
    tmpres.multiplyLeft(aux);
    //    tmpres=reshape(permute(aux,Indices(2,3,1,4)),Indices(al1*dd,du*al2))
    //*reshape(tmpres,Indices(du*al2,-1));
    tmpres.reshape(Indices(al1,-1));
    tmpres.permute(Indices(2,1));
    //    tmpres=permute(reshape(tmpres,Indices(al1,-1)),Indices(2,1));
    alphan=alphan*al1;
  }
  // The last contraction, with data[n-1]
  if(dagger)
    aux=permute(conjugate(data[n-1]),Indices(3,2,1,4));
  else
    aux=data[n-1];
  int dun=aux.getDimension(0);
  al1n=aux.getDimension(1);
  int ddn=aux.getDimension(2);
  al2n=aux.getDimension(3);

  tmpres.reshape(Indices(d*al2n*dd,-1));
  aux.permute(Indices(2,1,4,3));
  aux.reshape(Indices(al1n,dun*al2n*ddn));
  tmpres.multiplyLeft(aux);
  //tmpres=reshape(permute(aux,Indices(2,1,4,3)),Indices(al1n,dun*al2n*ddn))*
  // reshape(tmpres,Indices(d*al2n*dd,-1));
  alphan=alphan*al1n;
  tmpres.reshape(Indices(al1n*Dlp,-1));
  tmpres.permute(Indices(2,1));
  tmpres.reshape(Indices(Dl,alphan,Dlp));
  result=tmpres;
  //  result=reshape(permute(reshape(tmpres,Indices(al1n*Dlp,-1)),Indices(2,1)),
  //		 Indices(Dl,alphan,Dlp));
}

void JoinedOperator::contractMketjop(mwArray& result,const mwArray& termL,
				     const Site& ket,const mwArray& termR,
				     bool dagger) const{
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int betabr=termR.getDimension(0);int betabl=termL.getDimension(0);
  int alphar=termR.getDimension(1);int alphal=termL.getDimension(1);  
  int betakr=termR.getDimension(2);int betakl=termL.getDimension(2);
#ifdef CHECKDIMS
  bool wrong=0;
  int du=dagger?data[0].getDimension(2):data[0].getDimension(0);
  int al1=data[0].getDimension(1);
  int dd=dagger?data[0].getDimension(0):data[0].getDimension(2);
  int al2=data[0].getDimension(3);
  int alphaL=al1;int alphaR=al2;
  for(int k=1;k<n;k++){  // k=2:n      
    int duk=dagger?data[k].getDimension(2):data[k].getDimension(0);
    int al1k=data[k].getDimension(1);
    int ddk=dagger?data[k].getDimension(0):data[k].getDimension(2);
    int al2k=data[k].getDimension(3);
    if(dd!=duk) wrong=1;
    alphaL=alphaL*al1k;alphaR=alphaR*al2k;      
    du=duk;dd=ddk;
  }
  // After the last one, check left and right dimensions
  if((alphaL!=alphal)||(alphaR!=alphar)||(Dlp!=betakl)||
     (Drp!=betakr)||(dp!=dd))
    wrong=1;
  if(wrong){
    cout<<"Wrong dimensions in JoinedOperator::contractMketjop:"
	<<" termL="<<termL<<", ket="<<ket<<", termR="<<termR
	<<", n="<<n;
    for(int k=0;k<n;k++)
      cout<<", oper["<<k<<"]="<<data[k];
    cout<<endl;
    exit(212);
  }
#else
  int du,dd; // Just define them
#endif
  // First step: contract Aket with tmpL
  mwArray tmpres(termL);
  tmpres.reshape(Indices(betabl*alphal,betakl));
  mwArray aux(ket.getA()); // dp x Dlp x Drp
  aux.permute(Indices(2,1,3)); // Dlp x dp x Drp
  aux.reshape(Indices(Dlp,dp*Drp));
  tmpres.multiplyRight(aux); // betabl*alphal x dp*Drp
  //  mwArray tmpres=reshape(termL,Indices(betabl*alphal,betakl))*
  //  reshape(permute(ket.getA(),Indices(2,1,3)),Indices(Dlp,dp*Drp));
  // Now shift the first Xi to first place and contract alphas and d
  // adding operators one by one, starting in n
  int remain=alphal;int passed=1;
  tmpres.reshape(Indices(betabl*remain*dp,Drp));
  tmpres.permute(Indices(2,1)); // Drp x betabl*alphal*dp
  //  tmpres=permute(reshape(tmpres,Indices(betabl*remain*dp,Drp)),Indices(2,1));
  for(int k=n-1;k>=0;k--){  // k=n:-1:1
    mwArray op=data[k];
    if(dagger) op=permute(conjugate(data[k]),Indices(3,2,1,4));
    du=op.getDimension(0);int al1=op.getDimension(1);
    int d2=op.getDimension(2);int al2=op.getDimension(3);
    remain=remain/al1;
    tmpres.reshape(Indices(passed*Drp*betabl*remain,al1*d2));
    aux=op; // du x al1 x dd x alp1
    aux.permute(Indices(2,3,1,4));
    aux.reshape(Indices(al1*d2,du*al2));
    tmpres.multiplyRight(aux); // passed*Drp*betabl*remain x du*al2
    //tmpres=reshape(tmpres,Indices(passed*Drp*betabl*remain,al1*d2))*
    //reshape(permute(op,Indices(2,3,1,4)),Indices(al1*d2,du*al2));
    tmpres.reshape(Indices(passed*Drp*betabl*remain*du,al2));
    tmpres.permute(Indices(2,1));
    //    tmpres=permute(reshape(tmpres,Indices(passed*Drp*betabl*remain*du,al2)),
    //		   Indices(2,1));
    passed=passed*al2;    
  }
  if((passed!=alphar)||(remain!=1)){
    cout<<"Not all operators contracted in contractMketjop"<<endl;
    exit(212);
  }
  // Now contract with tmpR
  tmpres.reshape(Indices(passed*Drp,betabl*du));
  aux=termR;aux.reshape(Indices(betabr,alphar*betakr));
  tmpres.multiplyLeft(aux);
  //  tmpres=reshape(termR,Indices(betabr,alphar*betakr))*
  //  reshape(tmpres,Indices(passed*Drp,betabl*du));
  tmpres.reshape(Indices(betabr,betabl,du));
  tmpres.permute(Indices(3,2,1));
  tmpres.reshape(Indices(-1,1));
  result=tmpres;
  //  result=reshape(permute(reshape(tmpres,Indices(betabr,betabl,du)),
  //			 Indices(3,2,1)),Indices(-1,1));
}

void JoinedOperator::contractMbrajop(mwArray& result,const mwArray& termL,
				     const Site& bra,const mwArray& termR,
				     bool dagger) const{
  int d=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int betabr=termR.getDimension(0);int betabl=termL.getDimension(0);
  int alphar=termR.getDimension(1);int alphal=termL.getDimension(1);  
  int betakr=termR.getDimension(2);int betakl=termL.getDimension(2);
#ifdef CHECKDIMS
  bool wrong=0;
  int du=dagger?data[0].getDimension(2):data[0].getDimension(0);
  int al1=data[0].getDimension(1);
  int dd=dagger?data[0].getDimension(0):data[0].getDimension(2);
  int al2=data[0].getDimension(3);
  int alphaL=al1;int alphaR=al2;
  for(int k=1;k<n;k++){  // k=2:n      
    int duk=dagger?data[k].getDimension(2):data[k].getDimension(0);
    int al1k=data[k].getDimension(1);
    int ddk=dagger?data[k].getDimension(0):data[k].getDimension(2);
    int al2k=data[k].getDimension(3);
    if(dd!=duk) wrong=1;
    alphaL=alphaL*al1k;alphaR=alphaR*al2k;      
    du=duk;dd=ddk;
  }
  // After the last one, check left and right dimensions
  if((alphaL!=alphal)||(alphaR!=alphar)||(Dl!=betabl)||
     (Dr!=betabr)||(d!=du))
    wrong=1;
  if(wrong){
    cout<<"Wrong dimensions in JoinedOperator::contractMbrajop:"
	<<" termL="<<termL<<", bra="<<bra<<", termR="<<termR
	<<", n="<<n;
    for(int k=0;k<n;k++)
      cout<<", oper["<<k<<"]="<<data[k];
    cout<<endl;
    exit(212);
  }
#else
  // Define the appropriate dims
  int du,dd;
#endif
  // First step: contract Abra with tmpL
  mwArray tmpres(bra.getA());
  tmpres.conjugate();
  tmpres.permute(Indices(3,1,2));
  tmpres.reshape(Indices(Dr*d,Dl));
  mwArray aux(termL);
  aux.reshape(Indices(betabl,alphal*betakl));
  tmpres.multiplyRight(aux);
  //mwArray tmpres=reshape(permute(conjugate(bra.getA()),Indices(3,1,2)),
  //			 Indices(Dr*d,Dl))*
  //reshape(termL,Indices(betabl,alphal*betakl));
  // Now shift the first Xi to last place and contract alphas and d
  // adding operators one by one, starting in n
  int remain=alphal;int passed=1;int d2;
  tmpres.reshape(Indices(Dr,d*remain*betakl));
  tmpres.permute(Indices(2,1));
  //tmpres=permute(reshape(tmpres,Indices(Dr,d*remain*betakl)),Indices(2,1));
  for(int k=0;k<n;k++){
    mwArray op=data[k];
    if(dagger) op=permute(conjugate(data[k]),Indices(3,2,1,4));
    du=op.getDimension(0);int al1=op.getDimension(1);
    d2=op.getDimension(2);int al2=op.getDimension(3);
    remain=remain/al1;   
    tmpres.reshape(Indices(du*al1,remain*betakl*Dr*passed));
    aux=op;
    aux.permute(Indices(1,2,4,3));
    aux.reshape(Indices(du*al1,al2*d2));
    aux.permute(Indices(2,1));
    tmpres.multiplyLeft(aux);
    //tmpres=permute(reshape(permute(op,Indices(1,2,4,3)),
    //			   Indices(du*al1,al2*d2)),Indices(2,1))*
    //reshape(tmpres,Indices(du*al1,remain*betakl*Dr*passed));
    tmpres.reshape(Indices(al2,d2*remain*betakl*Dr*passed));
    tmpres.permute(Indices(2,1));
    //tmpres=permute(reshape(tmpres,Indices(al2,d2*remain*betakl*Dr*passed)),
    //		   Indices(2,1));
    passed=passed*al2; 
  }
  if((passed!=alphar)||(remain!=1)){
    cout<<"Not all operators contracted in contractMbrajop"<<endl;
    exit(212);
  }
  // Now contract with tmpR
  tmpres.reshape(Indices(d2*betakl,Dr*alphar));
  aux=termR;aux.reshape(Indices(betabr*alphar,betakr));
  tmpres.multiplyRight(aux);
  tmpres.reshape(Indices(d2*betakl*betakr,1));
  result=tmpres;
  //tmpres=reshape(tmpres,Indices(d2*betakl,Dr*alphar))*
  //reshape(termR,Indices(betabr*alphar,betakr));
  result=reshape(tmpres,Indices(d2*betakl*betakr,1));
}

void JoinedOperator::contractNjop(mwArray& result,const mwArray& termL,
				  const mwArray& termR,bool dagger) const{
  int betabr=termR.getDimension(0);int betab=termL.getDimension(0);
  int alphar=termR.getDimension(1);int alpha=termL.getDimension(1);  
  int betakr=termR.getDimension(2);int betak=termL.getDimension(2);
#ifdef CHECKDIMS
  bool wrong=0;
  int du=dagger?data[0].getDimension(2):data[0].getDimension(0);
  int al1=data[0].getDimension(1);
  int dd=dagger?data[0].getDimension(0):data[0].getDimension(2);
  int al2=data[0].getDimension(3);
  int alphaL=al1;int alphaR=al2;
  for(int k=1;k<n;k++){ // k=2:n      
    int duk=dagger?data[k].getDimension(2):data[k].getDimension(0);
    int al1k=data[k].getDimension(1);
    int ddk=dagger?data[k].getDimension(0):data[k].getDimension(2);
    int al2k=data[k].getDimension(3);
    if(dd!=duk) wrong=1;
    alphaL=alphaL*al1k;alphaR=alphaR*al2k;      
    du=duk;dd=ddk;
  }
  // After the last one, check left and right dimensions
  if((alphaL!=alpha)||(alphaR!=alphar))
      wrong=1;
  if(wrong){
    cout<<"Wrong dimensions in JoinedOperator::contractNjop:"
	<<" termL="<<termL<<", termR="<<termR
	<<", n="<<n;
    for(int k=0;k<n;k++)
      cout<<", oper["<<k<<"]="<<data[k];
    cout<<endl;
    exit(212);
  }
#else
  // Just define the dims
  int du,al1,dd,al2;
#endif
  // Now do the contractions:
  // Start with topmost operator
  mwArray tmpres=permute(termL,Indices(2,3,1));
  // The first one gets contracted only in alpha
  mwArray op=data[0];
  if(dagger) op=permute(conjugate(op),Indices(3,2,1,4));
  du=op.getDimension(0);al1=op.getDimension(1);
  dd=op.getDimension(2);al2=op.getDimension(3);
  tmpres=reshape(permute(op,Indices(1,3,4,2)),Indices(du*dd*al2,al1))
    *reshape(tmpres,Indices(al1,-1));
  int remain=alpha/al1;int passed=al2;int dphys=du;
  tmpres=reshape(permute(reshape(tmpres,Indices(du,dd,al2,remain*betak*betab)),
			 Indices(2,4,1,3)),Indices(dd,remain*betak*betab*dphys*passed));
  for(int k=1;k<n;k++){ // k=2:n
    op=data[k];
    if(dagger) op=permute(conjugate(op),Indices(3,2,1,4));
    du=op.getDimension(0);al1=op.getDimension(1);
    dd=op.getDimension(2);al2=op.getDimension(3);
    remain=remain/al1;
    tmpres=permute(reshape(op,Indices(du*al1,dd*al2)),Indices(2,1))
      *reshape(tmpres,Indices(du*al1,remain*betak*betab*dphys*passed));
    tmpres=reshape(permute(reshape(tmpres,
				   Indices(dd,al2,remain*betak*betab*dphys*passed)),
			   Indices(1,3,2)),Indices(dd,-1));
    passed=passed*al2;
  }
  if((remain!=1)||(passed!=alphar)){
    cout<<"Not all operators contracted in contractNjop"<<endl;
    exit(212);
  }
  // Now contract with tmpR
  tmpres=reshape(tmpres,Indices(dd*betak*betab*dphys,passed))*
    reshape(permute(termR,Indices(2,1,3)),Indices(alphar,betabr*betakr));
  result=reshape(permute(reshape(tmpres,Indices(dd,betak,betab,dphys,betabr,betakr)),
			 Indices(4,3,5,1,2,6)),Indices(dphys*betab*betabr,dd*betak*betakr));
}
