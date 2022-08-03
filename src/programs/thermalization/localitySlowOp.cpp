#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0


using namespace shrt;
using namespace std;

/** Compute numerically the terms appearing in the coefficient of the slowl relaxing operator */


//Construct random unitaries W on two sites
void constructRandomUnitary(mwArray& W);
void constructUnitariesIsing(mwArray& W,mwArray& V,double g,double h,double tau);
void constructUnitarySwap(mwArray& W);

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

void fitLength(MPS& variable,int len);


void constructMPOrange(MPO& mpoL,int range,int len);

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
  const char* outfnameR=argv[++cntr];
  int rSeed=atoi(argv[++cntr]);

  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(1E-8);
  cout<<"Initialized Contractor"<<endl;
  srandom(rSeed);
  ofstream* out;
  out=new ofstream(outfnameR);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfnameR<<
      " for output"<<endl;
    exit(1);
  }
  *out<<"% D="<<D<<", M="<<M<<", rSeed="<<rSeed<<endl;
  *out<<setprecision(15);

  vector<double> g; // results

  mwArray Olv,Orv,Olw,Orw;
  // 1) Construct random unitaries W and V
  mwArray W,V;
  int nv,nw;
  //constructRandomUnitary(W);
  //constructRandomUnitary(V);
  constructUnitarySwap(W);
  constructUnitarySwap(V);  
  //constructUnitariesIsing(W,V,gP,hP,tau);
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
  MPS state(L,1,d*d),stateL(L,1,d*d);
  if(!traceless){
    initializeOptimTraceState(state);
    initializeOptimTraceState(stateL,true);
    //    initializeState(state);
    //initializeState(stateL,true);
  }
  else{
    initializeOptimTracelessState(state);
    initializeOptimTracelessState(stateL,true);
  }
  //cout<<"Initialized state(s) "<<endl;
  //  state.exportForMatlab("/afs/ipp/u/banulsm/Op.m");

  vector<MPS> On(M+1,MPS(1,1,1)),OLn(M+1,MPS(1,1,1));
  // where to save the results
  mwArray R(Indices(4*(M+1),2*M+1,2*M+1));

  // 4) Apply m steps of evolution
  for(int k=0;k<=M;k++){
    On[k]=state;
    OLn[k]=stateL; // k=0 is in both, but easier to keep this way
    // Compute overlaps
    int len=state.getLength();
    int lenL=stateL.getLength();
    for(int r=len;r>=1;r--){
      MPO mpoR(len);
      constructMPOrange(mpoR,r,len); // the projector
      for(int n=0;n<=k;n++){
	// pos pos terms
	int len2=On[n].getLength();
	if(len2>=r){ // pos pos
	  //cout<<"Overlap r="<<r<<" of O["<<k<<"]("<<len<<") with O["<<n<<"]("<<len2<<")"<<endl;
	  MPS auxPos=On[n];
	  fitLength(auxPos,len);
	  complex_t value=contractor.contract(state,mpoR,auxPos);
	  R.setElement(value,Indices(r,n+M,k+M));
	  *out<<r<<"\t"<<n<<"\t"<<k<<"\t"<<real(value)<<endl;
	}
	// neg pos
	len2=OLn[n].getLength();
	if(len2>=r){
	  //cout<<"Overlap r="<<r<<" of O["<<k<<"]("<<len<<") with O["<<-n<<"]("<<len2<<")"<<endl;
	  if(len2<=len){
	    MPS auxPos=OLn[n]; // neg pos
	    fitLength(auxPos,len);
	    complex_t value=contractor.contract(state,mpoR,auxPos);
	    R.setElement(value,Indices(r,M-n,k+M));
	    *out<<r<<"\t"<<-n<<"\t"<<k<<"\t"<<real(value)<<endl;
	  }
	  else{
	    //cout<<" ******* O["<<-n<<"] has length "<<len2<<", larger than O["<<k<<"]("<<len<<")"<<endl;
	    // ENlarge the MPO and the state
	    MPS auxState(state);
	    fitLength(auxState,len2);
	    MPO mpoExtra(len2);
	    constructMPOrange(mpoExtra,r,len2);
	    //cout<<"Now contracting extended O["<<k<<"]: "<<auxState<<" with O["<<-n<<"]: "<<OLn[n]<<" and MPO "<<mpoExtra<<endl;
	    complex_t value=contractor.contract(auxState,mpoExtra,OLn[n]);
	    R.setElement(value,Indices(r,M-n,k+M));
	    *out<<r<<"\t"<<-n<<"\t"<<k<<"\t"<<real(value)<<endl;
	  }
	}
      }
    }
    if(k>0){
      for(int r=lenL;r>=1;r--){
	MPO mpoR(lenL);
	// For stateL I need MPO of length lenL, which is 2 longer!!
	constructMPOrange(mpoR,r,lenL); // the projector
	for(int n=0;n<=k;n++){ // pos neg
	  int len2=On[n].getLength();
	  if(len2>=r){
	    //	    cout<<"Overlap r="<<r<<" of O["<<-k<<"]("<<lenL<<") with O["<<n<<"]("<<len2<<")"<<endl;
	    MPS auxPos=On[n];
	    fitLength(auxPos,lenL);
	    complex_t value=contractor.contract(stateL,mpoR,auxPos);
	    R.setElement(value,Indices(r,n+M,M-k));	  
	    *out<<r<<"\t"<<n<<"\t"<<-k<<"\t"<<real(value)<<endl;
	  }
	  // Same for neg neg
	  if(n>0){
	    len2=OLn[n].getLength();
	    if(len2>=r){
	      //cout<<"Overlap r="<<r<<" of O["<<-k<<"]("<<lenL<<") with O["<<-n<<"]("<<len2<<")"<<endl;
	      MPS auxPos=OLn[n]; // neg pos
	      fitLength(auxPos,lenL);
	      complex_t value=contractor.contract(stateL,mpoR,auxPos);
	      R.setElement(value,Indices(r,M-n,M-k));
	      *out<<r<<"\t"<<-n<<"\t"<<-k<<"\t"<<real(value)<<endl;
	    }
	  }
	}
      }
    }
    // Evolve one more
    if(k<M){
    //    cout<<"Starting evolution step "<<k+1<<" state("<<state.getLength()<<"), stateL("<<stateL.getLength()<<")"<<endl;
    if(k>0){
      // enlarge the MPS
      extendMPS(state);
      L=L+2;
    }
    // apply "even" part (V)
    MPO mpo(L);
    formEvenDoubleMPO(mpo,L,Olv,Orv);
    {
      MPS aux(state);
      contractor.optimize(mpo,aux,state,D);
      //   cout<<"After applying the even part (V) on O "<<endl;
    }
    // enlarge the MPS
    extendMPS(state);
    L=L+2;
    //cout<<"And extended state has <|>="<<contractor.contract(state,state)<<endl;
    // apply "odd" part
    formOddDoubleMPO(mpo,L,Olw,Orw);
    {
      MPS aux(state);
      contractor.optimize(mpo,aux,state,D);
      //cout<<"After applying the odd part (W) on O "<<endl;
    }
    // And now I evolve the left side
    // enlarge the MPS
    extendMPS(stateL);
    formOddDoubleMPO(mpo,stateL.getLength(),Olw,Orw,true);
    {
      MPS aux(stateL);
      contractor.optimizeL(mpo,aux,stateL,D);
      //      cout<<"After applying the odd part (W) on O+ "<<endl;
    }
    // enlarge the MPS
    extendMPS(stateL);
    // apply "odd" part
    formEvenDoubleMPO(mpo,stateL.getLength(),Olv,Orv,true);
    {
      MPS aux(stateL);
      contractor.optimizeL(mpo,aux,stateL,D);
      //      cout<<"After applying the even part (V) on O+ "<<endl;
    }
    //    cout<<"After evolution step "<<k+1<<" state("<<state.getLength()<<"), stateL("<<stateL.getLength()<<")"<<endl;
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


void fitLength(MPS& state,int len){
  int L=state.getLength();
  if(L==len) return;
  int n=len-L;
  if(n%2!=0){
    cout<<"ERROR: Can only fitLength in a symmetric way!"<<endl;
    exit(1);
  }
  MPS aux(state);
  state=MPS(len,state.getBond(),d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  int firstI=n/2;
  for(int k=0;k<n/2;k++)
    state.replaceSite(k,(1./sqrt(2))*id2,false);
  for(int k=0;k<n/2;k++)
    state.replaceSite(n/2+L+k,(1./sqrt(2))*id2,false);
  for(int k=0;k<L;k++)
    state.replaceSite(k+n/2,aux.getA(k),false);
  state.setNormFact(aux.getNormFact());
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


void constructMPOrange(MPO& mpoL,int range,int len){
  //cout<<"In MPO(L="<<len<<"), range="<<range<<endl;
  if(len!=mpoL.getLength())
    mpoL.initLength(len);
  if(range>len){
    cout<<"ERROR! range cannot be larger than MPO!"<<endl;
    exit(1);
  }
  // Basic operators, construct them once, then reuse
  static bool init=false;
  static mwArray Z(Indices(3,d*d,d*d));
  if(!init){
    mwArray O0=.5*identityMatrix(d*d);
    O0.reshape(Indices(d,d,d,d));
    O0.permute(Indices(2,4,1,3));
    O0.reshape(Indices(d*d,d*d));
    mwArray O1=identityMatrix(d*d);
    // construct the Z mwArray
    for(int i1=0;i1<d*d;i1++)
      for(int i2=0;i2<d*d;i2++){
	Z.setElement(O0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
	Z.setElement(O1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
	Z.setElement(O1.getElement(Indices(i1,i2))-O0.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      }
    Z.reshape(Indices(3,d*d*d*d));
    init=true;
    cout<<"Initialized the Z operator for the MPO"<<endl;
  }

  // DEBUGGING ONLY: I build the MPS alone to check coeffs
  MPS aux(len,range+1,3);

  // MPO terms. Special case, k=1=L
  int D=range+1;int Dl(D),Dr(D);
  for(int p=0;p<len;p++){
    if(p==0) Dl=1; else Dl=D;
    if(p==len-1) Dr=1; else Dr=D;
    mwArray C(Indices(Dl,Dr,3));
    if(len==1){ // special case: just one term
      C.setElement(ONE_c,Indices(Dl,Dr,2));
    }
    else{
      if(p<len-range)
	C.setElement(ONE_c,Indices(0,0,0));
      if(p>range-1)
	C.setElement(ONE_c,Indices(Dl-1,Dr-1,0));
      if(p<=len-range)
	if(p<len-1)
	  C.setElement(ONE_c,Indices(0,1,2)); // left operator
	else
	  C.setElement(ONE_c,Indices(0,0,2)); // left operator is also last
      for(int ik=1;ik<range-1;ik++){
	if(p>=ik&&(len-1-p)>=range-ik-1) // place for the rest of the term on both sides
	  C.setElement(ONE_c,Indices(ik,ik+1,1)); // in between
      }
      if(p>=range-1)
	C.setElement(ONE_c,Indices(range-1,Dr-1,2)); // in between
    }
    aux.replaceSite(p,permute(C,Indices(3,1,2)));
    //cout<<"Coefficient of site "<<p<<", C="<<C<<endl;
    C.reshape(Indices(Dl*Dr,3));
    C.multiplyRight(Z);
    C.reshape(Indices(Dl,Dr,d*d,d*d));
    C.permute(Indices(3,1,4,2));
    mpoL.setOp(p,new Operator(C),true); // no reusing operators!
  }
  // stringstream nameMPS;
  // nameMPS<<"mpsForRange"<<range<<"_L"<<len<<".m";
  // aux.exportForMatlab(nameMPS.str().data());
}
