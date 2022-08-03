#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0


using namespace shrt;
using namespace std;

/** Compute numerically the terms appearing in the coefficient of the
    slowly relaxing operator, in the case of a TI superposition which
    is then filtered */


//Construct random unitaries W on two sites
void constructRandomUnitary(mwArray& W);
void constructUnitariesIsing(mwArray& W,mwArray& V,double g,double h,double tau);
void constructUnitarySwap(mwArray& W);

// Split it in two terms
void split2term(mwArray& Ol,mwArray& Or,mwArray& U,int& nr);

void formEvenDoubleMPO(MPO& mpo,int L,const mwArray& Olv,const mwArray& Orv,bool left=false);
void formOddDoubleMPO(MPO& mpo,int L,const mwArray& Olw,const mwArray& Orw,bool left=false);

void initializeState(MPS& state,bool even,bool dagger=false);
void initializeTracelessState(MPS& state,bool even,bool dagger=false);
void initializeOptimTracelessState(MPS& state,bool even,bool dagger=false);
void initializeOptimTraceState(MPS& state,bool even,bool dagger=false);

// Compute and summ all non-vanishing (valid) overlaps of a left and
// right vector, assuming the leftmost site of the right has the same
// parity as the rightmost of the left
complex_t sumOverlaps(const MPS& state,const MPS& stateL,ofstream* output=0);


void extendMPS(MPS& state);
void extendMPSwithId(MPS& state);
void extendMPSside(MPS& state,int n,bool left);

int d=2;
double gP=0.9045;
double hP=0.8090;
double tau=0.8;

//bool traceless=true;
bool traceless=false;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int M=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfnameG=argv[++cntr];
  int rSeed=atoi(argv[++cntr]);

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
  *out<<"% D="<<D<<", M="<<M<<", rSeed="<<rSeed<<endl;
  *out<<setprecision(15);

  vector<double> gEE,gEO,gOE,gOO; // results

  mwArray Olv,Orv,Olw,Orw;
  // 1) Construct random unitaries W and V
  mwArray W,V;
  int nv,nw;
  constructRandomUnitary(W);
  constructRandomUnitary(V);
  // //constructUnitarySwap(W);
  // //constructUnitarySwap(V);  
  //  constructUnitariesIsing(W,V,gP,hP,tau);
  cout<<"Chosen W="<<W<<endl;
  cout<<"Chosen V="<<V<<endl;
  // 2) Split them in left and right parts, so that I can put them in MPOs
  split2term(Olv,Orv,V,nv);
  split2term(Olw,Orw,W,nw);
  //  cout<<"Constructed and splitted random unitaries "<<endl;
  //  cout<<"Olv="<<Olv<<endl;
  // cout<<"Orv="<<Orv<<endl;

  // 3) Initialize the "state" for the evolution: an MPS with double indices, started as OpxId/sqrt(2)
  int L=2;
  MPS stateE(L,1,d*d),stateEL(L,1,d*d);
  MPS stateO(L,1,d*d),stateOL(L,1,d*d);
  if(!traceless){
    initializeOptimTraceState(stateE,true,false);
    initializeOptimTraceState(stateEL,true,true);
    initializeOptimTraceState(stateO,false,false);
    initializeOptimTraceState(stateOL,false,true);
  }
  else{
    initializeOptimTracelessState(stateE,true,false);
    initializeOptimTracelessState(stateEL,true,true);
    initializeOptimTracelessState(stateO,false,false);
    initializeOptimTracelessState(stateOL,false,true);
  }

  //    ofstream* debug=new ofstream("allGs.txt");
  ofstream* debug=0;

  // 4) Apply m steps of evolution
  for(int k=1;k<=M+1;k++){
    // compute g(2*(k-1))
    MPS lastE(stateE); // length L, normalizes with sqrt for L-1
    MPS lastO(stateO); // length L, normalizes with sqrt for L-1
    // if(k>1){
    //   extendMPSwithId(lastE); // to match left part
    //   extendMPSwithId(lastO); // to match left part
    // }
    complex_t gtmpEE,gtmpEO,gtmpOE,gtmpOO;
    //    *debug<<"% contribs "<<2*(k-1)<<endl;
    gtmpEE=sumOverlaps(lastE,stateEL,debug);
    gtmpEO=sumOverlaps(lastO,stateEL,debug);
    gtmpOE=sumOverlaps(lastE,stateOL,debug);
    gtmpOO=sumOverlaps(lastO,stateOL,debug);
    // check for reality
    if(abs(imag(gtmpEE))>1E-10&&abs(real(gtmpEE)/imag(gtmpEE))<1E5){
      cout<<"Apparently there is an error, as gEE("<<2*(k-1)<<")="<<gtmpEE<<" is not real"<<endl;
      //exit(1);
    }
    gEE.push_back(real(gtmpEE));
    gEO.push_back(real(gtmpEO));
    gOE.push_back(real(gtmpOE));
    gOO.push_back(real(gtmpOO));
    *out<<2*(k-1)<<"\t"<<real(gtmpEE)<<"\t"<<imag(gtmpEE)
	<<"\t"<<real(gtmpEO)<<"\t"<<imag(gtmpEO)
	<<"\t"<<real(gtmpOE)<<"\t"<<imag(gtmpOE)
	<<"\t"<<real(gtmpOO)<<"\t"<<imag(gtmpOO)
	<<endl;
    //cout<<"g("<<2*(k-1)<<")\t"<<gtmp<<endl;
    if(k>1){
      // enlarge the MPS
      extendMPS(stateE);
      extendMPS(stateO);
      L=L+2;
    }
    //    cout<<"Starting step "<<k<<", extended state has <|>="<<contractor.contract(state,stateL)<<endl;
    // apply "even" part (V)
    MPO mpo(L);
    formEvenDoubleMPO(mpo,L,Olv,Orv);
    {
      MPS aux(stateE);
      contractor.optimize(mpo,aux,stateE,D);
      aux=stateO;
      contractor.optimize(mpo,aux,stateO,D);
      //   cout<<"After applying the even part (V) on O "<<endl;
    }
    // enlarge the MPS
    extendMPS(stateE);extendMPS(stateO);
    L=L+2;
    //cout<<"And extended state has <|>="<<contractor.contract(state,state)<<endl;
    // apply "odd" part
    formOddDoubleMPO(mpo,L,Olw,Orw);
    {
      MPS aux(stateE);
      contractor.optimize(mpo,aux,stateE,D);
      aux=stateO;
      contractor.optimize(mpo,aux,stateO,D);
      //cout<<"After applying the odd part (W) on O "<<endl;
    }
    // compute g(2*k-1) contracting with last
    lastE=stateEL;//extendMPSwithId(lastE);
    lastO=stateOL;//extendMPSwithId(lastO);
    //    *debug<<"% contribs "<<2*(k-1)+1<<endl;
    gtmpEE=sumOverlaps(stateE,lastE,debug);
    gtmpEO=sumOverlaps(stateO,lastE,debug);
    gtmpOE=sumOverlaps(stateE,lastO,debug);
    gtmpOO=sumOverlaps(stateO,lastO,debug);
    // check for reality
    if(abs(imag(gtmpEE))>1E-10&&abs(real(gtmpEE)/imag(gtmpEE))<1E5){
      cout<<"Apparently there is an error, as gEE("<<2*k-1<<")="<<gtmpEE<<" is not real"<<endl;
      //exit(1);
    }
    gEE.push_back(real(gtmpEE));
    gEO.push_back(real(gtmpEO));
    gOE.push_back(real(gtmpOE));
    gOO.push_back(real(gtmpOO));
    *out<<2*(k-1)+1<<"\t"<<real(gtmpEE)<<"\t"<<imag(gtmpEE)
	<<"\t"<<real(gtmpEO)<<"\t"<<imag(gtmpEO)
	<<"\t"<<real(gtmpOE)<<"\t"<<imag(gtmpOE)
	<<"\t"<<real(gtmpOO)<<"\t"<<imag(gtmpOO)
	<<endl;
    //cout<<"g("<<2*(k-1)+1<<")\t"<<gtmp<<endl;
    if(k<M+1){
      // And now I evolve the left side
      // enlarge the MPS
      extendMPS(stateEL);extendMPS(stateOL);
      formOddDoubleMPO(mpo,stateEL.getLength(),Olw,Orw,true);
      {
	MPS aux(stateEL);
	contractor.optimizeL(mpo,aux,stateEL,D);
	aux=stateOL;
	contractor.optimizeL(mpo,aux,stateOL,D);
	//      cout<<"After applying the odd part (W) on O+ "<<endl;
      }
      // enlarge the MPS
      extendMPS(stateEL);extendMPS(stateOL);
      // apply "odd" part
      formEvenDoubleMPO(mpo,stateEL.getLength(),Olv,Orv,true);
      {
	MPS aux(stateEL);
	contractor.optimizeL(mpo,aux,stateEL,D);
	aux=stateOL;
	contractor.optimizeL(mpo,aux,stateOL,D);
	//      cout<<"After applying the even part (V) on O+ "<<endl;
      }
    }
  }
  // extendMPSwithId(state); // to match left part
  // complex_t gtmp=contractor.contract(state,stateL);
  // // check for reality
  // if(abs(imag(gtmp))>1E-10&&abs(real(gtmp)/imag(gtmp))<1E5){
  //   cout<<"Apparently there is an error, as g("<<2*(M)<<")="<<gtmp<<" is not real"<<endl;
  //   //exit(1);
  // }
  // g.push_back(real(gtmp));
  // *out<<2*(M)<<"\t"<<real(gtmp)<<"\t"<<imag(gtmp)<<endl;
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

void initializeState(MPS& state,bool even,bool dagger){
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
  if(even){
    state.setA(0,Op);
    state.setA(1,(1./sqrt(2))*id2);
  }
  else{
    state.setA(1,Op);
    state.setA(0,(1./sqrt(2))*id2);
  }
  //state.setA(1,id2);
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized state with norm "<<contractor.contract(state,state)<<endl;
}

void initializeTracelessState(MPS& state,bool even,bool dagger){
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
  if(even){
    state.setA(0,Op);
    state.setA(1,(1./sqrt(2))*id2);
  }
  else{
    state.setA(1,Op);
    state.setA(0,(1./sqrt(2))*id2);
  }
  //state.setA(1,id2);
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized state with norm "<<contractor.contract(state,state)<<endl;
}

void initializeOptimTracelessState(MPS& state,bool even,bool dagger){
  // Initialize the operator O, which is a given combination of sigx sigy and sigz
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray sigX(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  mwArray Op=.0301*sigX-.6618*sigY+.7491*sigZ;
  //mwArray Op=sigY;
  // No need to conjugate, as it is done automatically in the Contractor!!!!!
  //  if(dagger){
    //    Op.conjugate();
  //}
  //  Op=identityMatrix(d);
  Op.reshape(Indices(d*d,1,1));
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  state=MPS(2,1,d*d);
  if(even){
    state.setA(0,Op);
    state.setA(1,(1./sqrt(2))*id2);
  }
  else{
    state.setA(1,Op);
    state.setA(0,(1./sqrt(2))*id2);
  }
  //state.setA(1,id2);
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized state with norm "<<contractor.contract(state,state)<<endl;
}

void initializeOptimTraceState(MPS& state,bool even,bool dagger){
  // Initialize the operator O, which is a given combination of sigx sigy and sigz
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray sigX(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  complex_t dataN[]={0.0301*ONE_c,-0.6618*ONE_c,0.7491*ONE_c};
  mwArray vN(Indices(3,1),dataN);
  double normV=norm(vN);vN=1/sqrt(normV)*vN;
  mwArray Op=real(vN.getElement(0))*sigX+real(vN.getElement(1))*sigY+real(vN.getElement(2))*sigZ;
  Op=.5*(Op+identityMatrix(d)); // convert to projector!  //mwArray Op=sigY;
  // No need to conjugate, as it is done automatically in the Contractor!!!!!
  //  if(dagger){
    //    Op.conjugate();
  //}
  //  Op=identityMatrix(d);
  Op.reshape(Indices(d*d,1,1));
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  state=MPS(2,1,d*d);
  if(even){
    state.setA(0,Op);
    state.setA(1,(1./sqrt(2))*id2);
  }
  else{
    state.setA(1,Op);
    state.setA(0,(1./sqrt(2))*id2);
  }
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

void formOddDoubleMPO(MPO& mpo,int L,const mwArray& Olv,const mwArray& Orv,bool left){
  // On the first and last sites, I have to contract the op with identity
  static bool init(false);
  static Operator *edgeL,*edgeR,*edgeLleft,*edgeRleft;
  static Operator *midL,*midR;
  if(!init){
    //    cout<<"Initializing double (odd) op with Ol="<<Olv<<" and Or="<<Orv<<endl;
    // Prepare the ops the first time
    // Double operators for the center
    midL=new DoubleOperator(Olv,conjugate(Olv));
    midR=new DoubleOperator(Orv,conjugate(Orv));
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
    mpo.setOp(k++,midR,false);
    mpo.setOp(k++,midL,false);    
  }
  //cout<<"Initialized formOddDouble "<<mpo<<endl;
}

void extendMPSside(MPS& state,int n,bool left){
  int L=state.getLength();
  MPS aux(state);
  state=MPS(L+n,state.getBond(),d*d);
  mwArray id2=1/sqrt(2)*identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  int firstI=left?0:L;
  int lastI=left?n-1:L+n-1;
  for(int p=firstI;p<=lastI;p++)
    state.replaceSite(p,id2,false);
  int firstA=left?n:0;
  for(int k=0;k<L;k++)
    state.replaceSite(k+firstA,aux.getA(k),false);
  state.setNormFact(aux.getNormFact());
  //  cout<<"Extended MPS to "<<state<<endl;
}

void compressMPSside(MPS& state,int n,bool left){
  int L=state.getLength();
  MPS aux(state);
  state=MPS(L-n,state.getBond(),d*d);
  int shift=left?n:0;
  for(int k=0;k<L-n;k++)
    state.replaceSite(k,aux.getA(k+shift),false);
  double oldNorm=aux.getNormFact();
  int firstRem=left?0:L-n;
  int lastRem=left?n-1:L-1;
  for(int k=firstRem;k<=lastRem;k++){
    mwArray auxA=aux.getA(k).getA();
    if(auxA.getDimension(0)*auxA.getDimension(1)*auxA.getDimension(2)!=d*d){
      cout<<"Error: wrong dimensions of removed site ("<<k<<") in compressMPS: A="<<auxA<<endl;
      exit(1);
    }
    auxA.reshape(Indices(d,d));mwArray checkA=auxA-1/sqrt(2)*identityMatrix(d);
    if(!checkA.isNull(1E-8)){
      cout<<"Potential error: in this program, I should only remove sites which are a(normalized) identity!!"
	  <<" but here I have tensor "<<auxA<<" with trace "<<auxA.trace()
	  <<endl;
      exit(1);
      //oldNorm*=abs(factor);
    }
  }
  state.setNormFact(oldNorm);
  // TODO: Check that sites 0 and L-1 are 1x1 and multiply also their values!
  //  cout<<"Extended MPS to "<<state<<endl;
}

complex_t sumOverlaps(const MPS& state,const MPS& stateL,ofstream* out){
  int Lr=state.getLength();
  int Ll=stateL.getLength();
  Contractor& contractor=Contractor::theContractor();

  // special case: if both lengths are 2 (step 0)
  if(Ll==2&&Lr==2) return contractor.contract(state,stateL)-.5*ONE_c;

  int offsetR=Ll-1;
  int offsetL=Lr-1;
//   cout<<"sumOverlaps received state "<<Lr<<" and stateL "<<Ll<<endl;
//   cout<<"offsetR(<-state)="<<offsetR<<", offsetL(stateL->)="<<offsetL<<endl;
  MPS auxSt(state);
  MPS auxStL(stateL);
  extendMPSside(auxSt,offsetR,true); // state is extended on the left by offsetR sites
  extendMPSside(auxStL,offsetL,false);
//   cout<<" state extended to "<<auxSt.getLength()<<endl;
//   cout<<" stateL extended to "<<auxStL.getLength()<<endl;
  complex_t result=ZERO_c;
  bool stop=false;
  while(!stop){
    complex_t tmp=contractor.contract(auxSt,auxStL);
    if(!traceless) tmp=tmp-.5*ONE_c;
    if(out!=0)
      *out<<tmp<<"\t";
//     cout<<"  Contraction "<<tmp<<endl;
    result=result+tmp;
    if(auxSt.getLength()-Lr+auxStL.getLength()-Ll==4) stop=true; // no further cut
    if(!stop){ // the last time cannot cut two
      compressMPSside(auxSt,2,true);  
      compressMPSside(auxStL,2,false);  
//       cout<<" state reduced to "<<auxSt.getLength()<<endl;
//       cout<<" stateL reduced to "<<auxStL.getLength()<<endl;
    }
  }
  // Out of the loop I should be in the position to contract "normally"
  // so I take the original ones again and extend the shorter
  auxSt=state;auxStL=stateL;
  MPS& shorter=Lr<Ll?auxSt:auxStL; 
  extendMPSwithId(shorter);
  //  cout<<" Extended state"<<(Lr<Ll?" ":"L ")<<"to "<<shorter.getLength()<<endl;
  complex_t tmp=contractor.contract(auxSt,auxStL);
  if(!traceless) tmp=tmp-.5*ONE_c;
  //cout<<"  Contraction "<<tmp<<endl;
  if(out!=0)
    *out<<tmp<<"\t";
  result=result+tmp;
  // And now I have to extend the state on the right and the stateL on the left
  auxSt=state;auxStL=stateL;
  offsetR=Ll>Lr?3:1;
  offsetL=Ll>Lr?1:3;
  //cout<<"Now to the other direction: offsetR(state->)="<<offsetR<<", offsetL(<-stateL)="<<offsetL<<endl;
  extendMPSside(auxSt,offsetR,false); // state is extended on the left by offsetR sites
  extendMPSside(auxStL,offsetL,true);
  //cout<<" state extended to "<<auxSt<<endl;
  //cout<<" stateL extended to "<<auxStL<<endl;
  stop=false;
  while(!stop){
    complex_t tmp=contractor.contract(auxSt,auxStL);
    if(!traceless) tmp=tmp-.5*ONE_c;
    //cout<<"  Contraction "<<tmp<<endl;
    if(out!=0)
      *out<<tmp<<"\t";
    result=result+tmp;
    if(auxSt.getLength()-Lr==Ll-1) stop=true;
    if(!stop){
      extendMPSside(auxSt,2,false);
      extendMPSside(auxStL,2,true);
      //cout<<" state extended to "<<auxSt.getLength()<<endl;
      //cout<<" stateL extended to "<<auxStL.getLength()<<endl;
    }
  }
    if(out!=0)
      *out<<endl;
  return result;
}

