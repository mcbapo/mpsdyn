#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "HeisenbergHamiltonian.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0


using namespace shrt;
using namespace std;

/** Compute numerically the terms appearing in the coefficient of the slowl relaxing operator */


void constructDoubleCommutatorHeis(MPO& mpo,int L,double h, mwArray instances, int rowInstances, vector<double> &randFields);

//Construct random unitaries W on two sites
void constructRandomUnitary(mwArray& W);
void constructUnitariesIsing(mwArray& W,mwArray& V,double g,double h,double tau);
void constructUnitarySwap(mwArray& W);

void computeTwoBodyExHam(mwArray& Ol,mwArray& Or,double Jx,double Jy,double Jz,double h1,double h2, complex_t delta);
void formDoubleMPO(MPO& mpo, int L, vector<mwArray> &UniMPO, bool left=false);


// Split it in two terms
void split2term(mwArray& Ol,mwArray& Or,mwArray& U,int& nr);

void formEvenDoubleMPO(MPO& mpo,int L,const mwArray& Olv,const mwArray& Orv,bool left=false);
void formOddDoubleMPO(MPO& mpo,int L,const mwArray& Olw,const mwArray& Orw,bool left=false);

void initializeState(MPS& state,bool dagger=false);
void initializeTracelessState(MPS& state,bool dagger=false);
void initializeOptimTracelessState(MPS& state,bool dagger=false);
void initializeOptimTraceState(MPS& state,bool dagger=false);
void extendMPS(MPS& state);
void extendMPSwithId(MPS& state);

int d=2;
double gP=0.9045;
double hP=0.8090;


//bool traceless=true;
bool traceless=true;
bool norand=false;
bool noHam=false;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int M=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  int rowInstances=atoi(argv[++cntr]);
  double h=atof(argv[++cntr]);
  double tau=atof(argv[++cntr]);;
  const char* outfnameG=argv[++cntr];
  const char* outfnameL=argv[++cntr];
  const char* infname=argv[++cntr];
  int rSeed=atoi(argv[++cntr]);
  
  mwArray instances;
  ifstream in(infname);
  if(!in.is_open()){
    cout<<"Error: couldn't open file "<<infname<<" for output"<<endl;
    exit(1);
  }
  instances.load(in);
  in.close();

  // Indices dims;
  // instances.getDimensions(dims);
  // cout<<"Dimensions:"<<dims<<endl;
  // for(int i=0; i<60; i++)
  //   cout<<instances.getElement(Indices(0,i))<<endl;
  // ofstream fileInst("instFile.txt");
  // //fileInst= new ofstream("instFile.txt");
  // instances.savetext(fileInst, true);
  // fileInst.close(); 
   
  vector<double> randField;
  int c=0,newc=0;
  while(c <= 4*(M+1) ) {
    if(c%96==0 && c>0) {rowInstances+=1; newc+=96;}
    randField.insert(randField.begin(), h*real(instances.getElement(Indices(rowInstances,c-newc))));c+=1;
    randField.push_back(h*real(instances.getElement(Indices(rowInstances,c-newc))));c+=1;    
  }

  // cout<<"inst vector: ";
  // for(int h=0; h<randField.size(); h++)
  //   cout<<randField.at(h)<<", ";
  // cout<<endl;

  HeisenbergHamiltonian hamH(4*(M+1) + 2,1.,1.,1.,randField,d);
  MPO expHe(4*(M+1) + 2); MPO expHo(4*(M+1) + 2);
  hamH.getExponentialMPOeven(expHe,-tau*I_c);
  hamH.getExponentialMPOodd(expHo,-tau*I_c);
  int lenOp=expHe.getLength();

  vector<double> randFields(M+2);
  MPO H2(M);
  constructDoubleCommutatorHeis(H2,M,hP,instances,rowInstances, randFields);

  
  
  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(1E-8);
  cout<<"Initialized Contractor"<<endl;
  srandom(rSeed);
  ofstream* out;
  out=new ofstream(outfnameG);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfnameG<<
      " for output"<<endl;
    exit(1);
  }
  *out<<"# D="<<D<<", M="<<M<<", rSeed="<<rSeed<<endl;
  *out<<setprecision(10);

  vector<double> g; // results
  vector<mwArray> Even, Odd;
  vector<mwArray> EvenH, OddH;
  mwArray Ur,Ul;
  vector<double> hi;

   // 3) Initialize the "state" for the evolution: an MPS with double indices, started as OpxId/sqrt(2)
  int L=2;
  MPS state(L,1,d*d),stateL(L,1,d*d);
  int lenState=0;

  initializeOptimTracelessState(state);
  initializeOptimTracelessState(stateL,true);

  int s=0; 
  for(int k=1;k<=M+1;k++){
    // compute g(2*(k-1))
    MPS last(state); // length L, normalizes with sqrt for L-1
    if(k>1)
      extendMPSwithId(last); // to match left part
    
    complex_t gtmp=contractor.contract(last,stateL);
    // check for reality
    if(abs(imag(gtmp))>1E-10&&abs(real(gtmp)/imag(gtmp))<1E5){
      cout<<"Apparently there is an error, as g("<<2*(k-1)<<")="<<gtmp<<" is not real"<<endl;
      //exit(1);
    }
    g.push_back(real(gtmp));

    *out<<2*(k-1)<<"\t"<<real(gtmp)<<"\t"<<imag(gtmp)<<endl;
    //cout<<"g("<<2*(k-1)<<")\t"<<gtmp<<endl;
    if(k>1){
      // enlarge the MPS
      extendMPS(state);
      L=L+2;
    }
    //    cout<<"Starting step "<<k<<", extended state has <|>="<<contractor.contract(state,stateL)<<endl;
    // apply "even" part (V)
    MPO mpo(L), mpoH(L);
      
    lenState=state.getLength();
    for(int j=(lenOp-lenState)/2; j<(lenOp+lenState)/2; j++){
      EvenH.push_back(expHe.getOpData(j));
    }
    formDoubleMPO(mpo,L,EvenH);
   
    {
      MPS aux(state);
      contractor.optimize(mpo,aux,state,D);
      //   cout<<"After applying the even part (V) on O "<<endl;
    }
    EvenH.clear();
    mpo.clear();
    
    // enlarge the MPS
    extendMPS(state);
    L=L+2;
    
    lenState=state.getLength();
    for(int j=(lenOp-lenState)/2; j<(lenOp+lenState)/2; j++){
      OddH.push_back(expHo.getOpData(j));
    }
    cout<<"Length State: "<<lenState<<" length Op: "<<OddH.size()<<endl;
    
    formDoubleMPO(mpo,L,OddH);

    {
      MPS aux(state);
      contractor.optimize(mpo,aux,state,D);
      //cout<<"After applying the odd part (W) on O "<<endl;
    }
    mpo.clear();
    OddH.clear();
    // compute g(2*k-1) contracting with last
    last=stateL;extendMPSwithId(last);
    gtmp=contractor.contract(state,last);
    // check for reality

    if(abs(imag(gtmp))>1E-10&&abs(real(gtmp)/imag(gtmp))<1E5){
      cout<<"Apparently there is an error, as g("<<2*k-1<<")="<<gtmp<<" is not real"<<endl;
      //exit(1);
    }
    g.push_back(real(gtmp));
    *out<<2*(k-1)+1<<"\t"<<real(gtmp)<<"\t"<<imag(gtmp)<<endl;
    //cout<<"g("<<2*(k-1)+1<<")\t"<<gtmp<<endl;
    if(k<M+1){
    // And now I evolve the left side
    // enlarge the MPS
    extendMPS(stateL);

    lenState=stateL.getLength();
    for(int j=(lenOp-lenState)/2; j<(lenOp+lenState)/2; j++){
      OddH.push_back(expHo.getOpData(j));
    }
    formDoubleMPO(mpo,L,OddH,true);
   
    {
      MPS aux(stateL);
      contractor.optimizeL(mpo,aux,stateL,D);
      //      cout<<"After applying the odd part (W) on O+ "<<endl;
    }
    OddH.clear();
    mpo.clear();

    // enlarge the MPS
    extendMPS(stateL);
    // apply "odd" part
   

    lenState=stateL.getLength();
    for(int j=(lenOp-lenState)/2; j<(lenOp+lenState)/2; j++){
      EvenH.push_back(expHe.getOpData(j));
    }
    
    formDoubleMPO(mpo,stateL.getLength(),EvenH,true); 
    
    {
      MPS aux(stateL);
      contractor.optimizeL(mpo,aux,stateL,D);
      //      cout<<"After applying the even part (V) on O+ "<<endl;
    }
    EvenH.clear();
    mpo.clear();

    }
  }
                 
  out->close();
  delete out;
  
  cout<<g<<endl;

  // Now I can construct the matrices for the optimization problem, for m=1,..M
  out=new ofstream(outfnameL);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfnameL<<
      " for output"<<endl;
    exit(1);
  }
  *out<<"# D="<<D<<", M="<<M<<", rSeed="<<rSeed<<endl;
  *out<<setprecision(10);
  *out<<"# M\t min(lambda)"<<endl;

  for(int m=1;m<=M;m++){
    mwArray X(Indices(2*m+1,2*m+1)),Y(Indices(2*m+1,2*m+1));
    for(int i=0;i<2*m+1;i++){
      for(int j=0;j<2*m+1;j++){
	X.setElement(ONE_c*(g[abs(j-i)]-g[abs(j-i+1)]),Indices(i,j));
	if(traceless)
	  Y.setElement(ONE_c*(g[abs(j-i)]),Indices(i,j));
	else
	  Y.setElement(ONE_c*(g[abs(j-i)]-.5),Indices(i,j));
      }      
    }
    X=X+Hconjugate(X);
    
    
    mwArray Uy;vector<complex_t> Dy;
    if(m==5)
      cout<<"Y="<<Y<<endl;

    wrapper::eig(Y,Dy,Uy,true);
    
    mwArray eigv=Uy[0];

    if(m==5)
      cout<<"Y="<<Dy<<endl;

    int nr=0;double tol=1E-6;
    mwArray D2=invertDiag(sqrt(diag(Dy)),nr,tol);
    D2.multiplyLeft(Uy);
    X.multiplyRight(D2);
    X.multiplyLeft(Hconjugate(D2));
    mwArray Ux;vector<complex_t> Dx;
    X=X+Hconjugate(X);
    
    wrapper::eig(X,Dx,Ux,true);
    *out<<m<<"\t"<<real(Dx[0])<<endl;

  }
  out->close();
  delete out;
}



void constructRandomUnitary(mwArray& W){
  mwArray M(Indices(d*d,d*d));
  M.fillRandom();
  //  cout<<"Constructed random M "<<M<<endl;
  M=M+Hconjugate(M);
  vector<complex_t> D;
  wrapper::eig(M,D,W,true);
  // complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  // complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  // static bool first=true;
  // if(first){
  //   M=mwArray(Indices(d*d,1),dataX);first=false;
  // }
  // else
  //   M=mwArray(Indices(d*d,1),dataZ);
  // M.multiplyRight(permute(M,Indices(2,1)));
  // M.reshape(Indices(d,d,d,d));
  // M.permute(Indices(1,3,2,4));
  // M.reshape(Indices(d*d,d*d));
  // wrapper::expm(M,W,I_c);
  //cout<<"Constructed random unitary "<<W<<endl;
  //cout<<"W W+="<<W*Hconjugate(W)<<endl;
}

void constructUnitariesIsing(mwArray& W,mwArray& V,double g,double h,double tau){
  // The unitaries of the Ising model with transverse and parallel field
  // W is the product of single body X terms with z terms, and V is teh z term on contiguous sites
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray sigX(Indices(d,d),dataX);
  mwArray twoZ(sigZ+h*.5*identityMatrix(d));twoZ.reshape(Indices(d*d,1));
  twoZ.multiplyRight(permute(twoZ,Indices(2,1))); // (sigZ+h) x (sigZ+h)
  twoZ.reshape(Indices(d,d,d,d));twoZ.permute(Indices(1,3,2,4));
  twoZ.reshape(Indices(d*d,d*d));
  twoZ=twoZ-h*h*.25*identityMatrix(d*d);
  wrapper::expm(twoZ,V,-tau*I_c);
  wrapper::expm(sigX,W,-tau*g*I_c);
  W.reshape(Indices(d*d,1));
  W.multiplyRight(permute(W,Indices(2,1)));
  W.reshape(Indices(d,d,d,d));W.permute(Indices(1,3,2,4));
  W.reshape(Indices(d*d,d*d));
  W.multiplyRight(V);
}

void constructUnitarySwap(mwArray& W){
  W=identityMatrix(d*d);
  W.reshape(Indices(d,d,d,d));
  W.permute(Indices(1,2,4,3));
  W.reshape(Indices(d*d,d*d));
}

void split2term(mwArray& Ol,mwArray& Or,mwArray& U,int& nr){
  //cout<<"split2term U="<<U<<endl;
  int dimL(d),dimR(d);
  U.reshape(Indices(dimL,dimR,dimL,dimR));
  U.permute(Indices(1,3,2,4));
  U.reshape(Indices(dimL*dimL,dimR*dimR));
  // And compute the SVD
  mwArray S; // place for the singular values
  nr=0;
  double tol=0.;
  //  cout<<"Before SVD expH="<<expH<<endl;
  wrapper::svd(U,tol,nr,Ol,S,Or);
  // cout<<"After SVD"<<endl;
  // redistribute the singular values
  S=sqrt(S);
  // TODO: check no NaN???
  Ol.multiplyRight(S);
  Or.multiplyLeft(S);
  // reshape for operators (i'x 1 x i x alpha ))
  Ol.reshape(Indices(dimL,1,dimL,-1));
  Or.reshape(Indices(-1,dimR,dimR,1));
  Or.permute(Indices(2,1,3,4));  
  //cout<<"split2term Ol="<<Ol<<", Or="<<Or<<", nr="<<nr<<endl;
}

void initializeState(MPS& state,bool dagger){
  // Initialize the operator O, which is the projector on |0>
  mwArray Op(Indices(d,d));Op.fillWithZero();
  Op.setElement(ONE_c,Indices(0,0));
  // if(dagger){
  //   Op.conjugate();
  // }
  //  Op=identityMatrix(d);
  Op.reshape(Indices(d*d,1,1));
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  state=MPS(2,1,d*d);
  state.setA(0,Op);
  state.setA(1,(1./sqrt(2))*id2);
  //state.setA(1,id2);
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized state with norm "<<contractor.contract(state,state)<<endl;
}

void initializeTracelessState(MPS& state,bool dagger){
  // Initialize the operator O, which is the projector on |0>
  mwArray Op(Indices(d,d));Op.fillWithZero();
  Op.setElement(ONE_c,Indices(0,0));
  Op.setElement(-1.*ONE_c,Indices(1,1));
  //if(dagger){
  //Op.conjugate();
  //}
  //  Op=identityMatrix(d);
  Op.reshape(Indices(d*d,1,1));
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  state=MPS(2,1,d*d);
  state.setA(0,Op);
  state.setA(1,(1./sqrt(2))*id2);
  //state.setA(1,id2);
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized state with norm "<<contractor.contract(state,state)<<endl;
}

void initializeOptimTracelessState(MPS& state,bool dagger){
  // Initialize the operator O, which is a given combination of sigx sigy and sigz
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray sigX(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  //mwArray Op=.0301*sigX-.6618*sigY+.7491*sigZ;
  mwArray Op=sigZ;
  
  //mwArray Op=sigY;
  // No need to conjugate, as it is done automatically in the Contractor!!!!!
  //  if(dagger){
    //    Op.conjugate();
  //}
  //  Op=identityMatrix(d);
  Op.reshape(Indices(d*d,1,1));
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  state=MPS(2,1,d*d);
  state.setA(0,Op);
  state.setA(1,(1./sqrt(2))*id2);
  //state.setA(1,id2);
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized state with norm "<<contractor.contract(state,state)<<endl;
}

void initializeOptimTraceState(MPS& state,bool dagger){
  // Initialize the operator O, which is a given combination of sigx sigy and sigz
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray sigX(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  mwArray Op=.0301*sigX-.6618*sigY+.7491*sigZ;
  Op=.5*(Op+identityMatrix(d)); // convert to projector!  //mwArray Op=sigY;
  // No need to conjugate, as it is done automatically in the Contractor!!!!!
  //  if(dagger){
    //    Op.conjugate();
  //}
  //  Op=identityMatrix(d);
  Op.reshape(Indices(d*d,1,1));
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  state=MPS(2,1,d*d);
  state.setA(0,Op);
  state.setA(1,(1./sqrt(2))*id2);
  //state.setA(1,id2);
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized state with norm "<<contractor.contract(state,state)<<endl;
}

void extendMPS(MPS& state){
  int L=state.getLength();
  MPS aux(state);
  state=MPS(L+2,state.getBond(),d*d);
  mwArray id2(ONE_c);id2.reshape(Indices(1,1,1));
  state.replaceSite(0,id2,false);
  state.replaceSite(L+1,id2,false);
  for(int k=0;k<L;k++)
    state.replaceSite(k+1,aux.getA(k),false);
  state.setNormFact(aux.getNormFact());
  //  cout<<"Extended MPS to "<<state<<endl;
}

void extendMPSwithId(MPS& state){
  int L=state.getLength();
  MPS aux(state);
  state=MPS(L+2,state.getBond(),d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  state.replaceSite(0,(1./sqrt(2))*id2,false);
  state.replaceSite(L+1,(1./sqrt(2))*id2,false);
  for(int k=0;k<L;k++)
    state.replaceSite(k+1,aux.getA(k),false);
  state.setNormFact(aux.getNormFact());
  //  cout<<"Extended MPS to "<<state<<endl;
}

void formEvenDoubleMPO(MPO& mpo,int L,const mwArray& Olv,const mwArray& Orv,bool left){
  // On the first and last sites, I have to contract the op with identity
  static bool init(false);
  static Operator *edgeL,*edgeR,*edgeLleft,*edgeRleft;
  static Operator *midL,*midR;
  if(!init){
    //cout<<"Initializing double (even) op with Ol="<<Olv<<" and Or="<<Orv<<endl;
    // Prepare the ops the first time
    // Double operators for the center
    midL=new DoubleOperator(Olv,conjugate(Olv));
    midR=new DoubleOperator(Orv,conjugate(Orv));
    mwArray id2=(1./sqrt(2))*identityMatrix(d);id2.reshape(Indices(d*d,1));
    mwArray aux=midL->getFullData();
    Indices dims=aux.getDimensions(); // dd,1,dd,Dr2
    //cout<<"midL has dims "<<dims<<endl;
    aux.permute(Indices(1,2,4,3));aux.reshape(Indices(dims[0]*dims[1]*dims[3],dims[2])); // dd*1*Dr2,dd
    aux.multiplyRight(id2); 
    aux.reshape(Indices(dims[0],dims[1],1,dims[3])); // dd,1,1,Dr2
    edgeL=new Operator(aux);
    aux=midR->getFullData();dims=aux.getDimensions(); // dd,Dl2,dd,1
    //cout<<"midR has dims "<<dims<<endl;
    aux.reshape(Indices(dims[0]*dims[1],dims[2])); // dd*Dl2,dd(*1)
    aux.multiplyRight(id2); 
    aux.reshape(Indices(dims[0],dims[1],1,dims[3])); // dd,Dl2,1,1
    edgeR=new Operator(aux);
    //same but contracting identity on the other side, for left operator 
    id2.reshape(Indices(1,d*d));
    aux=midL->getFullData(); dims=aux.getDimensions(); // dd,1,dd,Dr2
    aux.reshape(Indices(dims[0],dims[1]*dims[2]*dims[3])); // dd,1*dd*Dr2
    aux.multiplyLeft(id2); 
    aux.reshape(Indices(1,dims[1],dims[2],dims[3]));
    edgeLleft=new Operator(aux);
    aux=midR->getFullData();dims=aux.getDimensions(); // dd,Dl2,dd,1
    //cout<<"midR has dims "<<dims<<endl;
    aux.reshape(Indices(dims[0],dims[1]*dims[2])); // dd,Dl2*dd(*1)
    aux.multiplyLeft(id2); 
    aux.reshape(Indices(1,dims[1],dims[2],dims[3])); // 1,Dl2,dd,1
    edgeRleft=new Operator(aux);
    init=true;
  }
  // Set the operators
  mpo.initLength(L);
  if(L==2){
    mpo.setOp(0,midL,false);
    mpo.setOp(L-1,midR,false);
  }
  else{
    if(!left){
      mpo.setOp(0,edgeL,false);
      mpo.setOp(L-1,edgeR,false);
    }
    else{
      mpo.setOp(0,edgeLleft,false);
      mpo.setOp(L-1,edgeRleft,false);
    }
  }
  int k=1;
  while(k<L-1){
    mpo.setOp(k++,midR,false);
    mpo.setOp(k++,midL,false);    
  }
}

void formOddDoubleMPO(MPO& mpo,int L,MPO& expHo,bool left){
  // On the first and last sites, I have to contract the op with identity
  static bool init(false);
  static Operator *edgeL,*edgeR,*edgeLleft,*edgeRleft;
  static Operator *midL,*midR, *mid;
  if(!init){
    //    cout<<"Initializing double (odd) op with Ol="<<Olv<<" and Or="<<Orv<<endl;
    // Prepare the ops the first time
    // Double operators for the center
    midL=new DoubleOperator(expHo.getOpData(0),conjugate(expHo.getOpData(0)));
    midR=new DoubleOperator(expHo.getOpData(L-1),conjugate(expHo.getOpData(L-1)));
    mwArray id2=(1./sqrt(2))*identityMatrix(d);id2.reshape(Indices(d*d,1,1));
    mwArray aux=midL->getFullData();
    Indices dims=aux.getDimensions(); // dd,1,dd,Dr2
    aux.permute(Indices(1,2,4,3));aux.reshape(Indices(dims[0]*dims[1]*dims[3],dims[2])); // dd*1*Dr2,dd
    aux.multiplyRight(id2); 
    aux.reshape(Indices(dims[0],dims[1],1,dims[3])); // dd,1,1,Dr2
    edgeL=new Operator(aux);
    aux=midR->getFullData();dims=aux.getDimensions(); // dd,Dl2,dd,1
    aux.reshape(Indices(dims[0]*dims[1],dims[2])); // dd*Dl2,dd(*1)
    aux.multiplyRight(id2); 
    aux.reshape(Indices(dims[0],dims[1],1,dims[3])); // dd,Dl2,1,1
    edgeR=new Operator(aux);
    

    //same but contracting identity on the other side, for left operator 
    id2.reshape(Indices(1,d*d));
    aux=midL->getFullData(); dims=aux.getDimensions(); // dd,1,dd,Dr2
    aux.reshape(Indices(dims[0],dims[1]*dims[2]*dims[3])); // dd,1*dd*Dr2
    aux.multiplyLeft(id2); 
    aux.reshape(Indices(1,dims[1],dims[2],dims[3]));
    edgeLleft=new Operator(aux);
    aux=midR->getFullData();dims=aux.getDimensions(); // dd,Dl2,dd,1
    //cout<<"midR has dims "<<dims<<endl;
    aux.reshape(Indices(dims[0],dims[1]*dims[2])); // dd,Dl2*dd(*1)
    aux.multiplyLeft(id2); 
    aux.reshape(Indices(1,dims[1],dims[2],dims[3])); // 1,Dl2,dd,1
    edgeRleft=new Operator(aux);
    init=true;
  }
  // Set the operators
  mpo.initLength(L);
  if(!left){
    mpo.setOp(0,edgeL,false);
    mpo.setOp(L-1,edgeR,false);
  }
  else{
    mpo.setOp(0,edgeLleft,false);
    mpo.setOp(L-1,edgeRleft,false);
  }
  int k=1;
  while(k<L-1){
    mid=new DoubleOperator(expHo.getOpData(k),conjugate(expHo.getOpData(k)));
    mpo.setOp(k++,mid,false);   
  }
  //cout<<"Initialized formOddDouble "<<mpo<<endl;
}


void computeTwoBodyExHam (mwArray& Or, mwArray& Ol,double Jx,double Jy,double Jz,double h1,double h2, complex_t delta) {
  
  static bool init(false);
  static mwArray XX,YY,ZZ,ZI,IZ;

  mwArray H12; mwArray expH;
  
  if(!init){
    
    complex_t dataXX[]={ZERO_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ZERO_c};
    XX=.25*mwArray(Indices(d*d,d*d),dataXX);
    complex_t dataYY[]={ZERO_c,ZERO_c,ZERO_c,-ONE_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,-ONE_c,ZERO_c,ZERO_c,ZERO_c};
    YY=.25*mwArray(Indices(d*d,d*d),dataYY);
    complex_t dataZZ[]={ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,-ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,-ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,ONE_c};
    ZZ=.25*mwArray(Indices(d*d,d*d),dataZZ);
    complex_t dataIZ[]={ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,-ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,-ONE_c};
    IZ=.5*mwArray(Indices(d*d,d*d),dataIZ);
    complex_t dataZI[]={ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,-ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,ONE_c,ZERO_c,ZERO_c,ZERO_c,ZERO_c,-ONE_c};
    ZI=.5*mwArray(Indices(d*d,d*d),dataZI);
    init=true;
  }
  H12=Jx*XX+Jy*YY+Jz*ZZ+h2*ZI+h1*IZ;
   
  wrapper::expm(H12,expH,delta);
  expH.reshape(Indices(d,d,d,d));
  expH.permute(Indices(1,3,2,4));
  expH.reshape(Indices(d*d,d*d));
  mwArray S; // place for the singular values
  int nr=0;
  double tol=0.;
  wrapper::svd(expH,tol,nr,Ol,S,Or);
  S=sqrt(S);
  Ol.multiplyRight(S);
  Or.multiplyLeft(S);
  Ol.reshape(Indices(d,1,d,-1));
  Or.reshape(Indices(-1,d,d,1));
  Or.permute(Indices(2,1,3,4));  
}

void formDoubleMPO(MPO& mpo, int L, vector<mwArray> &UniMPO, bool left){
  // On the first and last sites, I have to contract the op with identity

  static Operator *edgeL,*edgeR,*edgeLleft,*edgeRleft;
  static Operator *midL,*midR;
  
  mwArray Olv(UniMPO.front());
  mwArray Orv(UniMPO.back());

  //cout<<"Initializing double (even) op with Ol="<<Olv<<" and Or="<<Orv<<endl;
  // Prepare the ops the first time
  // Double operators for the center
  midL=new DoubleOperator(Olv,conjugate(Olv));
  midR=new DoubleOperator(Orv,conjugate(Orv));
  mwArray id2=(1./sqrt(2))*identityMatrix(d);id2.reshape(Indices(d*d,1));
  mwArray aux=midL->getFullData();
  Indices dims=aux.getDimensions(); // dd,1,dd,Dr2
  //cout<<"midL has dims "<<dims<<endl;

  aux.permute(Indices(1,2,4,3));aux.reshape(Indices(dims[0]*dims[1]*dims[3],dims[2])); // dd*1*Dr2,dd
  aux.multiplyRight(id2); 
  aux.reshape(Indices(dims[0],dims[1],1,dims[3])); // dd,1,1,Dr2
  edgeL=new Operator(aux);
  aux=midR->getFullData();dims=aux.getDimensions(); // dd,Dl2,dd,1
  //cout<<"midR has dims "<<dims<<endl;
  aux.reshape(Indices(dims[0]*dims[1],dims[2])); // dd*Dl2,dd(*1)
  aux.multiplyRight(id2); 
  aux.reshape(Indices(dims[0],dims[1],1,dims[3])); // dd,Dl2,1,1
  edgeR=new Operator(aux);
  //same but contracting identity on the other side, for left operator 
  id2.reshape(Indices(1,d*d));
  aux=midL->getFullData(); dims=aux.getDimensions(); // dd,1,dd,Dr2
  aux.reshape(Indices(dims[0],dims[1]*dims[2]*dims[3])); // dd,1*dd*Dr2
  aux.multiplyLeft(id2); 
  aux.reshape(Indices(1,dims[1],dims[2],dims[3]));
  edgeLleft=new Operator(aux);
  aux=midR->getFullData();dims=aux.getDimensions(); // dd,Dl2,dd,1
  //cout<<"midR has dims "<<dims<<endl;
  aux.reshape(Indices(dims[0],dims[1]*dims[2])); // dd,Dl2*dd(*1)
  aux.multiplyLeft(id2); 
  aux.reshape(Indices(1,dims[1],dims[2],dims[3])); // 1,Dl2,dd,1
  edgeRleft=new Operator(aux);

  // Set the operators
  mpo.initLength(L);
  if(L==2){
    mpo.setOp(0,midL,false);
    mpo.setOp(L-1,midR,false);
  }
  else{
    if(!left){
      mpo.setOp(0,edgeL,false);
      mpo.setOp(L-1,edgeR,false);
    }
    else{
      mpo.setOp(0,edgeLleft,false);
      mpo.setOp(L-1,edgeRleft,false);
    }
  }
  int k=1;
  vector<Operator*> mid;
  while(k<L-1){
    mid.push_back(new DoubleOperator(UniMPO[k],conjugate(UniMPO[k])));
    mpo.setOp(k++,mid[k-1],true);
  }
  
}


void constructDoubleCommutatorHeis(MPO& mpo,int L,double h, mwArray instances, int rowInstances, vector<double> &randFields){

  if(mpo.getLength()!=L)
    mpo.initLength(L);
  // Basic operators

  int nrOps=7;
  int dd=d*d;
  mwArray Z=mwArray(Indices(nrOps,dd,dd));
  // Basic pieces
  mwArray sig0=identityMatrix(d);//identity
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),datax);//sigmax
  complex_t datay[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),datay);//sigmay
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz

  // *** First of all, we need the identity on both
  mwArray term0=identityMatrix(dd); //id*id
  // *** Now the ones appearing in single-body terms
  // (1) sigX (x) Id
  mwArray term1;
  constructOperatorProduct(term1,sigX,sig0); //sigX*id
  // (2)  Id (x) sigX -> just a permutation of the previous one
  mwArray term2;
  constructOperatorProduct(term2,sig0,permute(sigX,Indices(2,1))); //id*sigX^T
  // (3) sigY (x) Id
  mwArray term3;
  constructOperatorProduct(term3,sigY,sig0); //sigY*id
  // (4)  Id (x) sigY^T -> just a permutation of the previous one
  mwArray term4;
  constructOperatorProduct(term4,sig0,permute(sigY,Indices(2,1))); //id*sigY^T
  // (5) sigZ (x) Id
  mwArray term5;
  constructOperatorProduct(term5,sigZ,sig0); //sigZ*id
  // (6)  Id (x) sigZ^T -> just a permutation of the previous one
  mwArray term6;
  constructOperatorProduct(term6,sig0,permute(sigZ,Indices(2,1))); //id*sigZ^T
  
  // Fill in the operators in the Z array, to be contracted with C in
  // the construction of the MPO. The order is as listed above.
  for(int d1=0;d1<dd;d1++)
    for(int d2=0;d2<dd;d2++){
      Z.setElement(term0.getElement(Indices(d1,d2)),Indices(0,d1,d2)); //id*id
      Z.setElement(term1.getElement(Indices(d1,d2)),Indices(1,d1,d2)); //sigX*id
      Z.setElement(term2.getElement(Indices(d1,d2)),Indices(2,d1,d2)); //id*sigX^T
      Z.setElement(term3.getElement(Indices(d1,d2)),Indices(3,d1,d2)); //sigY*id
      Z.setElement(term4.getElement(Indices(d1,d2)),Indices(4,d1,d2)); //id*sigY^T
      Z.setElement(term5.getElement(Indices(d1,d2)),Indices(5,d1,d2)); //sigZ*id
      Z.setElement(term6.getElement(Indices(d1,d2)),Indices(6,d1,d2)); //id*sigZ^T
    }
  Z.reshape(Indices(nrOps,dd*dd));
  //cout<<"Finished Z:"<<Z.getDimensions()<<endl;
  // now the coefficients
  int D=8; // bond dimension of the MPO
  mwArray C(Indices(D,D,nrOps)); // the middle site
  mwArray C1(Indices(1,D,nrOps)); // the first site
  mwArray CN(Indices(D,1,nrOps)); // the last site

  vector<mwArray> totOps;
  deque<double> hi;
  for(int k=0; k<=L+1; k++){
    
    double tmphi = h*real(instances.getElement(Indices(rowInstances,k)));
    if(norand){
      hi.push_back(.5);
    }
    if(k%2 == 0 && !norand){
      hi.push_back(tmphi);
      }else{ 
      hi.push_front(tmphi);
    }
    
  }
  
  // Single body terms (one for each operator, with proper weights)
  C1.setElement(ONE_c,Indices(0,0,0));      //  I x I
  C1.setElement(ONE_c,Indices(0,1,1));      //  sigX x I
  C1.setElement(-1.*ONE_c,Indices(0,2,2));  //  I x sigX 
  C1.setElement(ONE_c,Indices(0,3,3));      //  sigY x I
  C1.setElement(-1.*ONE_c,Indices(0,4,4));  //  I x sigY
  C1.setElement(ONE_c,Indices(0,5,5));      //  sigZ x I
  C1.setElement(-1.*ONE_c,Indices(0,6,6));  //  I x sigZ 
  C1.setElement(hi[0]*ONE_c,Indices(0,7,5));   // h1 sigZ x I
  C1.setElement(-hi[0]*ONE_c,Indices(0,7,6));  // -h1 I x sigZ^T

  CN.setElement(ONE_c,Indices(D-1,0,0));    //  I x I
  CN.setElement(ONE_c,Indices(1,0,1));      //  sigX x I
  CN.setElement(ONE_c,Indices(2,0,2));      //  I x sigX
  CN.setElement(ONE_c,Indices(3,0,3));      //  sigY x I 
  CN.setElement(ONE_c,Indices(4,0,4));      //  I x sigY
  CN.setElement(ONE_c,Indices(5,0,5));      //  sigZ x I 
  CN.setElement(ONE_c,Indices(6,0,6));      //  I x sigZ
  CN.setElement(hi[L+1]*ONE_c,Indices(0,0,5));   // hN sigZ x I
  CN.setElement(-hi[L+1]*ONE_c,Indices(0,0,6));  // -hN I x sigZ^T

  
  for(int ind = 0; ind<L; ind++){
    totOps.push_back(mwArray(Indices(D,D,nrOps)));
    totOps[ind].setElement(ONE_c,Indices(0,0,0));      //I x I
    totOps[ind].setElement(ONE_c,Indices(D-1,D-1,0));  //I x I
    
    for(int s=1; s<D-1; s++){
      totOps[ind].setElement(pow(-1.,s+1)*ONE_c,Indices(0,s,s));
      totOps[ind].setElement(ONE_c,Indices(s,D-1,s)); 
    }
    
    totOps[ind].setElement(hi[ind+1]*ONE_c,Indices(0,D-1,5)); // hi sigZ x I
    totOps[ind].setElement(-hi[ind+1]*ONE_c,Indices(0,D-1,6)); // -hi I x sigZ^T
    
    totOps[ind].reshape(Indices(D*D,nrOps));
    totOps[ind].multiplyRight(Z);
    totOps[ind].reshape(Indices(D,D,dd,dd));
    totOps[ind].permute(Indices(3,1,4,2));
    
  }

  for(int r=0;r<=L+2;r++)
    randFields[r] = hi[r];//transform deque into vector
  
  C1.reshape(Indices(1*D,nrOps));
  CN.reshape(Indices(D*1,nrOps));
  
  C1.multiplyRight(Z);C1.reshape(Indices(1*D*dd,dd));
  // Special situation: first and last sites are contracted with Id as traced out
  mwArray id2=1/sqrt(2.)*identityMatrix(d);id2.reshape(Indices(dd,1));
  C1.multiplyRight(id2);
  C1.reshape(Indices(1,D,dd,1));C1.permute(Indices(3,1,4,2)); //dd x 1 x 1 x D
  
  CN.multiplyRight(Z);CN.reshape(Indices(D*1*dd,dd));
  CN.multiplyRight(id2);
  CN.reshape(Indices(D,1,dd,1));CN.permute(Indices(3,1,4,2)); // dd x D x 1 x 1

  // Now I build the MPO with required ops.
  //1) on the edges, I contract C1 and CN with their adjoint
  C1.reshape(Indices(dd,D));
  CN.reshape(Indices(dd,D)); // dd D
  mwArray tmp(C1);tmp.Hconjugate(); 
  C1.multiplyLeft(tmp); // Dl,Dr
  tmp=CN;tmp.Hconjugate(); 
  CN.multiplyLeft(tmp); //Dlu,Dld
  
  // 2) These, I contract with the middle ones 
  tmp=totOps[0];tmp.permute(Indices(2,1,3,4));tmp.reshape(Indices(D,dd*dd*D));
  mwArray edgeL=C1*tmp;edgeL.reshape(Indices(D*dd,1,dd,D)); // Ddd 1 dd D
  tmp=totOps[L-1];tmp.reshape(Indices(dd*D*dd,D));
  CN.permute(Indices(2,1)); // put first the index for 2nd MPO
  mwArray edgeR=tmp*CN;edgeR.reshape(Indices(dd,D,dd,D));
  edgeR.permute(Indices(1,4,2,3));edgeR.reshape(Indices(dd*D,D,dd,1)); // ddD D dd 1
  
  // Now for the adjoint I also need modifications
  mwArray Cdl(totOps[0]);Cdl.permute(Indices(3,2,1,4),true); // adjoint of normal C
  mwArray Cdr(totOps[L-1]);Cdr.permute(Indices(3,2,1,4),true); // adjoint of normal C
  mwArray edgeLd(Cdl);edgeLd.reshape(Indices(dd,1,D*dd,D));
  mwArray edgeRd(Cdr);edgeRd.reshape(Indices(dd,D,dd*D,1));
  
  Operator edgeLOp(edgeL),edgeROp(edgeR);
  Operator edgeLdOp(edgeLd),edgeRdOp(edgeRd);
  const Operator* arrL[2]={&edgeLdOp,&edgeLOp};
  const Operator* arrR[2]={&edgeRdOp,&edgeROp};

  vector<Operator*> CdOp,COp;
  
  mpo.setOp(0,new JoinedOperator(2,arrL),true);
  mpo.setOp(L-1,new JoinedOperator(2,arrR),true);
  vector<mwArray> Cd(totOps);
  //cdp.push_back(new Operator(*(Cd[3])));

  
  if(L>2){
    
    for(int k=1;k<L-1;k++){
      
      Cd[k].permute(Indices(3,2,1,4),true);
      CdOp.push_back(new Operator(Cd[k]));
      COp.push_back(new Operator(totOps[k]));
      const Operator* arrC[2] = {CdOp[k-1],COp[k-1]};
      mpo.setOp(k,new JoinedOperator(2,arrC),true);  
    }
  }
}




  // 1) create the mbl Hamiltonian: even and odd MPOs included
  
  // vector<double> h;
  // deque<double> hi;
  // for(int k=0; k<=2*M+1; k++){
    
  //   double tmphi = h*real(instances.getElement(Indices(rowInstances,k)));
  //   if(norand){
  //     hi.push_back(1);
  //   }
  //   if(k%2 == 0 && !norand){
  //     hi.push_back(tmphi);
  //     }else{ 
  //     hi.push_front(tmphi);
  //   }
    
  // }
  // for(unsigned int k=0; k<hi.size(); k++) 
  //   h.push_back(hi.at(k));

  // HeisenbergHamiltonian hamH(2*M+1,1,1,1,h,d);
  // MPO expHe(2*M+1); MPO expHo(2*M+1);
  // hamH.getExponentialMPOeven(expHe,-tau*I_c);
  // hamH.getExponentialMPOodd(expHo,-tau*I_c);
     
