#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0


using namespace shrt;
using namespace std;

/** Compute numerically the largest lambda for an MPO of given bond dimensions and range M */

//Construct random unitaries W on two sites
void constructRandomUnitary(mwArray& W);
void constructUnitariesIsing(mwArray& W,mwArray& V,double g,double h,double tau);
void constructUnitarySwap(mwArray& W);

// Split it in two terms
void split2term(mwArray& Ol,mwArray& Or,mwArray& U,int& nr);

void formEvenDoubleMPO(MPO& mpo,int L,const mwArray& Olv,const mwArray& Orv,bool left=false);
void formOddDoubleMPO(MPO& mpo,int L,const mwArray& Olw,const mwArray& Orw,bool left=false);


void initializeState(MPS& state,bool dagger=false);
void initializeTracelessHermitianState(MPS& state);
void extendMPS(MPS& state);
void compressMPS(MPS& state);
void extendMPSwithId(MPS& state);

// substracts the trace by increasing the bond dimension in 1
// If I am not allowed to increase D, I will need to truncate back to the old one
void substractTrace(MPS& Omps);
complex_t computeTrace(const MPS& Omps);
void forceHermitian(MPS& Omps);

#define MAXITER 100
#define MAXPERT 10


int d=2;
double gP=0.9045;
double hP=0.8090;
double tau=0.8;

bool traceless=true;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int M=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  // const char* outfnameL=argv[++cntr];
  int rSeed=atoi(argv[++cntr]);

  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(1E-8);
  cout<<"Initialized Contractor"<<endl;
  srandom(rSeed);
   ofstream* out;
   out=new ofstream(outfname);
   if(!out->is_open()){
     cout<<"Error: impossible to open file "<<outfname<<
       " for output"<<endl;
     exit(1);
   }
   *out<<"% D="<<D<<", M="<<M<<", rSeed="<<rSeed<<endl;
   *out<<setprecision(10);
   *out<<"% iter\t eigenvalue\t lambda\t tr(O)\t <O|O>"<<endl;

  // 1) Construct random unitaries W and V
  mwArray W,V;
  int nv,nw;
  mwArray Olv,Orv,Olw,Orw;
  constructRandomUnitary(W);
  constructRandomUnitary(V);
  //  W=identityMatrix(d*d);
  //constructUnitarySwap(W);
  //constructUnitarySwap(V);
  //constructUnitariesIsing(W,V,gP,hP,tau);
  cout<<"Chosen W="<<W<<endl;
  cout<<"Chosen V="<<V<<endl;
  split2term(Olv,Orv,V,nv);
  split2term(Olw,Orw,W,nw);
  //  putForMatlab(cout,Olv,"Olv");  
  //putForMatlab(cout,Orv,"Orv");  
//  cout<<"Constructed and splitted random unitaries "<<endl;

  // I first construct the basic even and odd terms
  MPO UxUc(M+2),UdxUt(M+2); // the joined ones
  MPO mpoV(M+2),mpoW(M+2),mpoVd(M+2),mpoWd(M+2);
  formEvenDoubleMPO(mpoV,M+2,Olv,Orv); 
  formOddDoubleMPO(mpoW,M+2,Olw,Orw); // automatically 2 on each side are traced out
  const MPO* mpos[2]={&mpoW,&mpoV};
  MPO::join(2,mpos,UxUc);

  // And also the ones for the adjoint
  // Fast trick for the adjoints
  // V.Hconjugate();W.Hconjugate();
  // split2term(Olv,Orv,V,nv);
  // split2term(Olw,Orw,W,nw);
  // Now I conjugate and transpose by hand
  Olv.permute(Indices(3,2,1,4),true);
  Orv.permute(Indices(3,2,1,4),true);
  Olw.permute(Indices(3,2,1,4),true);
  Orw.permute(Indices(3,2,1,4),true);
  //putForMatlab(cout,Olv,"OlvD");  
  //putForMatlab(cout,Orv,"OrvD");  
  formEvenDoubleMPO(mpoVd,M+2,Olv,Orv,true);
  formOddDoubleMPO(mpoWd,M+2,Olw,Orw,true);
  const MPO* mposd[2]={&mpoVd,&mpoWd};
  MPO::join(2,mposd,UdxUt);
  // UdxUt=mpoVd;
  // UxUc=mpoV;

  cout<<"constructed all MPOs:"<<endl;
  // cout<<"mpoV looks like "<<mpoV<<endl;
  // cout<<"mpoW looks like "<<mpoW<<endl;
  // cout<<"U x Uc looks like "<<UxUc<<endl;
  //UxUc.exportForMatlab("UxUc.m");
  //UdxUt.exportForMatlab("UdxUt.m");
  //  cout<<"U x Uc looks like "<<UxUc<<endl;
  //Initialize the MPO O with random entries (make Hermitian and traceless)
  //  MPS Omps(M,D-1,d*d);
  MPS Omps(M,1,d*d);
  Omps.setRandomState();
  //cout<<"constructed initial MPS "<<Omps<<" with trace "<<computeTrace(Omps)<<" and <O|O>="<<contractor.contract(Omps,Omps)<<endl;
  forceHermitian(Omps);
  //cout<<"Initial MPS made Hermitian"<<Omps<<" with trace "<<computeTrace(Omps)<<" and <O|O>="<<contractor.contract(Omps,Omps)<<endl;
  substractTrace(Omps);
  //cout<<"Initial MPS made traceless ("<<computeTrace(Omps)<<")"<<endl;
  Omps.gaugeCond('R',1);Omps.gaugeCond('L',1); 
  //cout<<"After making it hermitian, substracting the trace and applying the gauge, trace is "<<computeTrace(Omps)
  //  <<" and <O|O>="<<contractor.contract(Omps,Omps)<<endl;
  //cout<<"Initialized state(s) "<<endl;
  //  state.exportForMatlab("/afs/ipp/u/banulsm/Op.m");

  // Now apply the power method
  bool conv=0;
  int count=0;int perturbcnt=0;
  double lambda=4;double convtol=1E-5;
  while(!conv){
    count++;
    // cout<<"Step "<<count<<" of the iteration, eigenval="<<lambda<<", lambda="<<2+2-lambda<<endl;
    extendMPS(Omps);
    //cout<<"Initial MPS extended to "<<Omps<<endl;
    Omps.gaugeCond('R',1);Omps.gaugeCond('L',1); 
    //cout<<"Initial MPS fulfilling gauge cond"<<endl;
    MPS aux1(Omps),aux2(Omps); 
    if(count<2){ // init guess
	aux1.approximate(UxUc,Omps,D);
	//cout<<"aux1 has been found that approximates the action of UxUc on Omps"<<endl;
	aux2.approximate(UdxUt,Omps,D);
	//cout<<"aux2 has been found that approximates the action of UdxUt on Omps"<<endl;
    }
    contractor.optimize(UxUc,Omps,aux1,D);
    //cout<<"aux1 has been found that optimizes the action of UxUc on Omps"<<endl;
    contractor.optimize(UdxUt,Omps,aux2,D);
    //cout<<"aux2 has been found that optimizes the action of UdxUt on Omps"<<endl;
    vector<const MPS*> mpss;
    mpss.push_back(&aux1);
    mpss.push_back(&aux2);
    vector<complex_t> coeffs(2,ONE_c);
    MPS aux3(Omps);
    mpss.push_back(&aux3);
    coeffs.push_back(2*ONE_c);
    //Omps.setRandomState();
    contractor.optimizeSum(mpss,coeffs,Omps,3*D);
    //  cout<<"Omps now optimized as the sum of aux1, aux2, 2 Omps"<<endl;

    double newLamb=Omps.getNormFact();
    // At the end, check convergence
    if(abs(newLamb-lambda)<abs((lambda)*convtol)){
      cout<<"Iteration apparently converged after "<<count
      	  <<" attempts with lambda="<<newLamb<<" and tol="<<convtol<<endl;
      conv=1;
    }
    else{
      //cout<<"last var="<<abs(newLamb-*lambda)/abs(*lambda)<<"; ";
      if(count>MAXITER){
	cout<<"No convergence after "<<count<<" trials"<<endl;
	Omps.perturb(.01); // Perturbation on 1 out of 100 sites
	count=0;
	perturbcnt++;
	if(perturbcnt>=MAXPERT){
	  cout<<"Could not find an eigenvector after "<<perturbcnt
	      <<" trials (each of them = perturbation-"<<MAXITER<<" applications"
	      <<" of the operator=> FAILED"<<endl;
	  exit(1); // TODO: Seguir de todas formas??
	}
      }    
    }
    lambda=newLamb;
    Omps.gaugeCond('R',1);Omps.gaugeCond('L',1); 
    compressMPS(Omps);
    //complex_t tr=computeTrace(Omps);
    //    cout<<"At the end of step "<<count<<", the Omps has trace "<<tr<<endl;
    forceHermitian(Omps);
    substractTrace(Omps);
    aux1=Omps;
    contractor.optimizeMPS(aux1,Omps,D);
    // cout<<"After step "<<count<< ", trace is "<<computeTrace(Omps)
    // 	<<" and <O|O>="<<contractor.contract(Omps,Omps)<<endl;
    *out<<count<<"\t"<<lambda<<"\t"<<4-lambda<<"\t"<<computeTrace(Omps)<<"\t"<<contractor.contract(Omps,Omps)<<endl;
  }

  cout<<"Final result lambda("<<M<<")="<<2+2-lambda<<endl;

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
  if(dagger){
    Op.conjugate();
  }
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
  if(dagger){
    Op.conjugate();
  }
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

void compressMPS(MPS& state){
  int L=state.getLength();
  if(state.getA(0).getDr()>1||state.getA(L-1).getDl()>1){
    cout<<"Cannot remove edges from MPS "<<state<<endl;
    exit(1);
  }
  MPS aux(state);
  state=MPS(L-2,state.getBond(),d*d);
  for(int k=0;k<L-2;k++)
    state.replaceSite(k,aux.getA(k+1),false);
  double oldNorm=aux.getNormFact();
  mwArray auxA=aux.getA(0).getA();
  if(auxA.getDimension(0)*auxA.getDimension(1)*auxA.getDimension(2)!=1){
    cout<<"Error: wrong dimensions of first site in compressMPS: A="<<auxA<<endl;
    exit(1);
  }
  oldNorm*=abs(auxA.getElement(0));
  auxA=aux.getA(L-1).getA();
  if(auxA.getDimension(0)*auxA.getDimension(1)*auxA.getDimension(2)!=1){
    cout<<"Error: wrong dimensions of last site in compressMPS: A="<<auxA<<endl;
    exit(1);
  }  
  oldNorm*=abs(auxA.getElement(0));
  state.setNormFact(oldNorm);
  // TODO: Check that sites 0 and L-1 are 1x1 and multiply also their values!
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
  if(L%2!=0){
    cout<<"ERROR: Odd M not yet supported!"<<endl;
    exit(1);
  }
  Operator *edgeL,*edgeR;
  Operator *midL,*midR;
  //cout<<"Initializing double (even) op with Ol="<<Olv<<" and Or="<<Orv<<endl;
  // Prepare the ops the first time
  // Double operators for the center
  midL=new DoubleOperator(Olv,conjugate(Olv));
  midR=new DoubleOperator(Orv,conjugate(Orv));
  mwArray id2=(1./sqrt(2))*identityMatrix(d);id2.reshape(Indices(d*d,1));
  mwArray aux;
  Indices dims; 
  if(!left){
    aux=midL->getFullData();
    dims=aux.getDimensions(); // dd,1,dd,Dr2
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
  }
  else{
    aux=midL->getFullData();dims=aux.getDimensions(); // dd,1,dd,Dr2
    aux.reshape(Indices(dims[0],dims[1]*dims[2]*dims[3])); // dd,1*dd*Dr2
    id2.reshape(Indices(1,d*d));
    aux.multiplyLeft(id2); // 1,1*dd*Dr2
    aux.reshape(Indices(1,dims[1],dims[2],dims[3])); // 1,1,dd,Dr2
    edgeL=new Operator(aux);
    aux=midR->getFullData();dims=aux.getDimensions(); // dd,Dl2,dd,1
    //cout<<"midR has dims "<<dims<<endl;
    aux.reshape(Indices(dims[0],dims[1]*dims[2]*dims[3])); // dd,Dl2*dd*1
    aux.multiplyLeft(id2); // 1,Dl2*dd*1 
    aux.reshape(Indices(1,dims[1],dims[2],dims[3])); // 1,Dl2,dd,1
    edgeR=new Operator(aux);
  }
  // Set the operators
  mpo.initLength(L);
  if(L==2){
    mpo.setOp(0,midL,true);
    mpo.setOp(L-1,midR,true);
  }
  else{
    mpo.setOp(0,edgeL,true);
    mpo.setOp(L-1,edgeR,true);
  }
  int k=1;
  bool saved=false;
  while(k<L-1){
    if(!saved){
      mpo.setOp(k++,midR,true);
      mpo.setOp(k++,midL,true);    
      saved=true;
    }
    else{
      mpo.setOp(k++,midR,false);
      mpo.setOp(k++,midL,false);    
    }
  }
}

void formOddDoubleMPO(MPO& mpo,int L,const mwArray& Olv,const mwArray& Orv,bool left){
  // On the first and last sites, I have to contract the op with identity
  Operator *edgeL,*edgeR;
  Operator *midL,*midR;
  // Double operators for the center
  midL=new DoubleOperator(Olv,conjugate(Olv));
  midR=new DoubleOperator(Orv,conjugate(Orv));
  mwArray id2=(1./sqrt(2))*identityMatrix(d);id2.reshape(Indices(d*d,1));

  mwArray aux=midL->getFullData();
  Indices dims=aux.getDimensions(); // dd,1,dd,Dr2
  aux.permute(Indices(1,2,4,3));aux.reshape(Indices(dims[0]*dims[1]*dims[3],dims[2])); // dd*1*Dr2,dd
  aux.multiplyRight(id2); 
  aux.reshape(Indices(dims[0],dims[1]*1*dims[3])); // dd,1*1*Dr2
  id2.reshape(Indices(1,d*d));
  aux.multiplyLeft(id2); // 1,Dr2
  mwArray aux2=midR->getFullData();  // dd,Dl2,dd,1
  Indices dims2=aux2.getDimensions(); // dd,Dl2,dd,1
  if(!left){ // extra id2 for traicing out the extra site
    aux2.reshape(Indices(dims2[0],dims2[1]*dims2[2])); // dd,Dl2*dd(*1)
    aux2.multiplyLeft(id2); 
    aux2.reshape(Indices(dims2[1],dims2[2]*dims2[3])); // Dl2,dd*1
  }
  else{ // the trace is applied on the down indices
    aux2.reshape(Indices(dims2[0]*dims2[1],dims2[2])); // dd,Dl2*dd(*1)
    id2.reshape(Indices(d*d,1));
    aux2.multiplyRight(id2);
    aux2.reshape(Indices(dims2[0],dims2[1])); // dd,Dl2
    aux2.permute(Indices(2,1));
  }
  aux2.reshape(Indices(dims2[1],-1)); // Dl2,dd*1
  aux2.multiplyLeft(aux); // 1,dd*1
  if(!left)
    aux2.reshape(Indices(1,1,dims2[2],1)); // 1,1,dd,1
  else
    aux2.reshape(Indices(dims2[0],1,1,1)); // dd,1,1,1
  edgeL=new Operator(aux2);

  // Now the same for right side
  aux=midR->getFullData();dims=aux.getDimensions(); // dd,Dl2,dd,1
  aux.reshape(Indices(dims[0]*dims[1],dims[2])); // dd*Dl2,dd(*1)
  id2.reshape(Indices(d*d,1));
  aux.multiplyRight(id2); 
  aux.reshape(Indices(dims[0],dims[1]*1*dims[3])); // dd,Dl2*1*1
  id2.reshape(Indices(1,d*d));aux.multiplyLeft(id2);  //1,Dl2
  aux.reshape(Indices(dims[1],1));
  aux2=midL->getFullData(); dims2=aux2.getDimensions(); // dd,1,dd,Dr2
  if(!left){
    aux2.reshape(Indices(dims2[0],dims2[1]*dims2[2]*dims2[3])); // dd,1*dd*Dr2
    aux2.multiplyLeft(id2);  // 1,1*dd*Dr2
    aux2.reshape(Indices(dims2[2],dims2[3])); // 1*dd,Dr2
  }
  else{
    aux2.permute(Indices(3,1,2,4));
    aux2.reshape(Indices(dims2[2],dims2[0]*dims2[1]*dims2[3])); // dd,dd*1*Dr2
    aux2.multiplyLeft(id2);  // 1,1*dd*Dr2
    aux2.reshape(Indices(dims2[0],dims2[3])); // dd*1,Dr2
  }
  aux2.multiplyRight(aux); // dd,
  if(!left)
    aux2.reshape(Indices(1,1,dims[2],1));
  else
    aux2.reshape(Indices(dims[0],1,1,1));
  edgeR=new Operator(aux2);

  // Set the operators
  mpo.initLength(L);
  mpo.setOp(0,edgeL,true);
  mpo.setOp(L-1,edgeR,true);
  int k=1;bool stored=false;
  while(k<L-1){
    if(!stored){
      mpo.setOp(k++,midL,true);
      mpo.setOp(k++,midR,true);
      stored=true;    
    }
    else{
      mpo.setOp(k++,midL,false);
      mpo.setOp(k++,midR,false);    
    }
  }
  //cout<<"Initialized formOddDouble "<<mpo<<endl;
}

void substractTrace(MPS& Omps){
  int L=Omps.getLength();
  // First compute trace
  complex_t tr=computeTrace(Omps);
  if(abs(imag(tr))>1E-10&&abs(imag(tr)/real(tr))>1E-3){
    cout<<"Apparently there is a  problem, as the trace is not real! tr(O)="<<tr<<endl;
    exit(1);
  }
  double signTr=real(tr)>=0?1:-1;
  double val=(1./d)*pow(abs(real(tr)),1./L);
  for(int k=0;k<L;k++){
    mwArray A=Omps.getA(k).getA();
    Indices dims=A.getDimensions(); // dd,Dl,Dr
    int newDl=k==0?dims[1]:dims[1]+1;
    int newDr=k==L-1?dims[2]:dims[2]+1;
    A.resize(Indices(d*d,newDl,newDr));
    A.reshape(Indices(d,d,newDl,newDr));
    double sign=k==0?-1.*signTr:1.;
    for(int j=0;j<d;j++)
      A.setElement(sign*val*ONE_c,Indices(j,j,newDl-1,newDr-1)); 
    A.reshape(Indices(d*d,newDl,newDr));
    Omps.replaceSite(k,A,false);
  }
}

complex_t computeTrace(const MPS& Omps){
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  int len=Omps.getLength();
  MPS aux(len,1,d*d);
  for(int k=0;k<len;k++){
    aux.setA(k,id2);
  }
  Contractor& contractor=Contractor::theContractor();
  return contractor.contract(Omps,aux);
}
 
void forceHermitian(MPS& Omps){
  int len=Omps.getLength();
  for(int k=0;k<len;k++){
    mwArray A=Omps.getA(k).getA();
    Indices dims=A.getDimensions(); // dd,Dl,Dr
    A.reshape(Indices(d,d,dims[1],dims[2]));
    mwArray B(A);
    B.permute(Indices(2,1,3,4),true);
    A=.5*(A+B);    
    A.reshape(Indices(d*d,dims[1],dims[2]));
    Omps.setA(k,A);
  }
}
