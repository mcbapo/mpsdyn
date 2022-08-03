#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "SpinMPO.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "HeisenbergAncillaryNoiseHamiltonian.h"

using namespace shrt;

/** mutualInfoTI simulates the evolution of a mixed state, which
    initially is maximally mixed for every site except for L0, in the
    first chain, initially in a product state (Z+ or X+).
    The program simulates the evolution and computes the 2-Renyi 
    entropies of different parts of the chain, to determine the mutual information.
    To keep track of how good the evolution is, I will also compute
    the energy along the way.

    Receives arguments:
    \param <L> (int) number of sites (in the original chain)
    \param <L0> (int) number of sites (in the original chain)
    \param <isXY> (bool) whether the model to be considered is only XY (no ZZ term) 
    \param <J> (double) parameter \f$J\f$ of the Heisenberg
                  Hamiltonian for the ancillary chain
    \param <B> (double) parameter \f$B\f$ of the interaction
                        between both chains
    \param <init> (char) initial state in the middle, can be 1 (X),
                        2(Y) or 3 (Z)
    \param <delta> (double) width of the time step used
    \param <M> (int) number of steps
    \param <rate> (int) frequency of recording data (nr of steps)
                        Results will be saved for 0,rate,2*rate...,M
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the results
                        (will be replaced if it exists, unless newInstance==0)
    \param <mpsfile> (char*) name of the file where to save the MPS 
    \param <newInstance> (bool) whether to start a new evolution 
			(if not, the MPS in mpsfile is used)
    \param <M0> (int)  only read if newInstance is false, it
                       is the number of steps already done in the
		       existing file.
*/


/** Product of all identities, but single site operator Op on site pos */
void constructLocal(int L,int pos,mwArray op,int d,MPS& result);
/** Vectorize an MPO to make an MPS out of it. If up==true (default)
    the first component of the double physical index will be the one
    which was going up.*/
void MPSfromMPO(const MPO& mpo,MPS& mps,bool up=true);
void MPOfromMPS(const MPS& mps,MPO& mpo,bool up=true);

/** doubleMPO taxes a normal MPO, acting on (system)
    physical dimensions, and constracts a double version, as in the
    evolution of a thermal state, for instance, where a second
    (transposed) copy of the MPO acts on the (ancillary) double indices. */
void doubleMPO(const MPO& simpleMPO,MPO& doubleMPO);
/** extendMPO takes a normal MPO, and transforms it into an MPO acting
    on a mixed state only on one side (the system indices). On the
    rest of the double index (whose dimension must be given in dimA,
    but is typically the same as for the physical) only the identity
    acts.*/
void extendMPO(const MPO& simpleMPO,MPO& doubleMPO,int dimA);
/** Analogous to the one above, but puts the transposed MPO onto the ancillary indices*/
void extendTransposeMPO(const MPO& simpleMPO,MPO& doubleMPO,int dimA);

/** Substitute a single site in the MPS by the one corresponding to
    the single site operator initSt indicates*/
void setSingleInitState(MPS& rho0,int L,int d,int initSt,int pos);

/** Substitute L0 sites in the MPS by the state that encodes initSt in L0 sites, 
    starting at pos0. */
void setExtendedInitState(MPS& rho0,int L,int d,int initSt,int L0,int pos0);

/** Compute and print to ofstrem the expectation values of X, Y, Z on
    pos */
void computeLocalOp(const MPS& rho,int pos,int d,ofstream* out);

/** Compute the 2-Renyi entropies of two pieces of the chain, namely
    A, defined to be from pos1A until pos2A (both included and
    numbered from 0 to L) and the complementary, B */
double compute2RenyiS(const MPS& redRho,complex_t traceR);
void compute2RenyiEntropiesAB(const MPS& rho,complex_t trR,int pos1A,int pos2A,double* SA,double*SB);
void computeEntropies(const MPS& rho,complex_t trR,ofstream* out,int L0,bool header=false);

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  int L0=atoi(argv[++cntr]);
  bool isXY=false;
  isXY=atoi(argv[++cntr])!=0;
  double J_=atof(argv[++cntr]);
  double B_=atof(argv[++cntr]);
  int initSt=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int rate=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  const char* mpsfile=argv[++cntr];
  int newInstance=atoi(argv[++cntr]);
  if(newInstance<0) newInstance=1;
  int M0=0;
  if(!newInstance&&argc>14) M0=atoi(argv[++cntr]);


  cout<<"Initialized arguments "
      <<", L="<<L
      <<", L0="<<L0
      <<", isXY?="<<isXY
      <<", J="<<J_
      <<", B="<<B_
      <<", delta="<<delta
      <<", M="<<M
      <<", D="<<D
      <<", outfile="<<outfname
      <<", mpsfile="<<mpsfile
      <<", newInstance="<<newInstance
      <<", M0="<<M0
      <<", rate="<<rate
      <<endl;

  ofstream* out;
  bool appendingFile=(!newInstance&&file_exists(outfname));
  if(appendingFile)
    out=new ofstream(outfname,ios::app);
  else // created again
    out=new ofstream(outfname);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  if(!appendingFile){
    *out<<"% L="<<L<<", J="<<J_<<", B="<<B_
	<<" D="<<D<<", init St="<<initSt<<endl;
    *out<<"% MPS in "<<mpsfile<<endl;
    *out<<"% M\t delta\t t\t D\t Re[tr(rho H)]\t Im[tr(rho H)]\t Re[tr(rho)] \tIm[tr(rho)] \t";
    computeEntropies(MPS(2*L,1,1),I_c,out,L0,true);
    *out<<endl;
  }
  *out<<setprecision(15);
  int Mtot=M0+M;

  Contractor& contractor=Contractor::theContractor();
  //contractor.setConvTol(1E-10);
  //cout<<"Initialized Contractor"<<endl;
  int d=2;

  // First: create the Hamiltonian
  double J1x(1.),J1y(1.),J1z(isXY?0.:1.);
  double J2x(J_),J2y(J_),J2z(isXY?0.:J_);
  HeisenbergAncillaryNoiseHamiltonian hamH(L,J1x,J1y,J1z,J2x,J2y,J2z,B_);
  //hamH.getHMPO().exportForMatlab("hamilHeisAnc.m");  
  //cout<<"Created the Hamiltonian"<<endl;

  // The evolution operator for a rho will result from double operators
  // Which MPOs do I need for the evolution
  MPO expHSe(2*L),expHSo(2*L),expHSo2(2*L),
    expHTe(2*L),expHTo(2*L),expHTo2(2*L),
    expHB(2*L),expHB2(2*L);
  hamH.getDoubleExponentialMPOeven(expHSe,-delta*I_c,true);
  hamH.getDoubleExponentialMPOodd(expHSo,-delta*I_c,true);
  hamH.getDoubleExponentialMPOodd(expHSo2,-delta*.5*I_c,true);
  hamH.getDoubleExponentialMPOeven(expHTe,-delta*I_c,false);
  hamH.getDoubleExponentialMPOodd(expHTo,-delta*I_c,false);
  hamH.getDoubleExponentialMPOodd(expHTo2,-delta*.5*I_c,false);
  hamH.getDoubleExponentialMPOmixed(expHB,-delta*I_c);
  hamH.getDoubleExponentialMPOmixed(expHB2,-delta*.5*I_c);

  // cout<<"Created all the exponential operators:"<<endl;
  // cout<<"expHeSigma:"<<expHSe<<endl;
  // cout<<"expHeTau:"<<expHTe<<endl;  
  // cout<<"expHoSigma:"<<expHSo<<endl;  
  // cout<<"expHoTau:"<<expHTo<<endl;  
  // cout<<"expHoSigma(1/2):"<<expHSo2<<endl;  
  // cout<<"expHoTau(1/2):"<<expHTo2<<endl;  
  // cout<<"expHB:"<<expHB<<endl; 

  // Now construct the initial rho MPO, in this case, the identity,
  // except for one site 
  MPS rho0(2*L,1,d*d);
  rho0.setProductState(p_maxent);
  int posCenter=L; // position in the chain where I set the different operator
  //  setSingleInitState(rho0,L,d,initSt,L);
  setExtendedInitState(rho0,L,d,initSt,L0,(L-L0)/2);

  if(!newInstance)
    rho0.importMPS(mpsfile);
  rho0.gaugeCond('R',1); // normalize
  
  // Now do the evolution step by step, and compute the energy,
  // and the renyi entropies
  MPO hamil(2*L);
  const MPO& _hamil=hamH.getHMPO();
  extendMPO(_hamil,hamil,d);

  // Also the identity, to compute the trace
  MPS Id(2*L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<2*L;k++)
    Id.setA(k,id2);

  int cnt=M0;
  double time=cnt*delta;
  complex_t Eev,trRho;
  Eev=contractor.contract(rho0,hamil,Id);
  trRho=contractor.contract(rho0,Id);
  double S2R=compute2RenyiS(rho0,trRho);
  //cout<<"Time evolution step nr "<<cnt<<", time="<<time
  //  <<", trace="<<trRho<<", <H>="<<Eev
  //  <<endl;
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"
      <<real(Eev)<<"\t"<<imag(Eev)<<"\t"<<real(trRho)<<"\t"<<imag(trRho)<<"\t";
  computeEntropies(rho0,trRho,out,L0);
  *out<<endl;
  // vector<complex_t> lambda2;
  // contractor.getSchmidtValues(rho0,lambda2); // lambda^2 in the middle
  // for(int id=0;id<D;id++)
  //   *out<<real(lambda2[id])<<"\t";
  // //  computeEntropies(rho0,out);
  // computeLocalSigZ(rho0,d,out);
  // //  computeErrors(rho0,out,Ds);
  *out<<endl;
  out->close(); // I probably don't want to keep it open for days!
  delete out;
  while(cnt<Mtot){
    MPS aux(rho0); // temporary copy
    contractor.optimize(expHSo2,rho0,aux,D);
    if(J_!=0)
      contractor.optimize(expHTo2,aux,rho0,D);
    else
      rho0=aux;
    // Now apply r-1 times the pair Ho He
    int cntLoop=0;
    while(cntLoop<rate&&cnt<Mtot){
      aux=rho0;  // AQUI!!!!*****
      contractor.optimize(expHB2,aux,rho0,D);
      contractor.optimize(expHSe,rho0,aux,D);
      if(J_!=0)
	contractor.optimize(expHTe,aux,rho0,D);
      else
	rho0=aux;
      contractor.optimize(expHB2,rho0,aux,D);
      if(cntLoop<rate-1){
	// contractor.optimize(expHSo,rho0,aux,D);
	// contractor.optimize(expHTo,aux,rho0,D);
	contractor.optimize(expHSo,aux,rho0,D);
	if(J_!=0)
	  contractor.optimize(expHTo,rho0,aux,D);
	else
	  aux=rho0;
      }
      else{
	// contractor.optimize(expHSo2,rho0,aux,D);
	// contractor.optimize(expHTo2,aux,rho0,D);
	contractor.optimize(expHSo2,aux,rho0,D);
	if(J_!=0)
	  contractor.optimize(expHTo2,rho0,aux,D);
	else
	  aux=rho0;
      }
      rho0=aux;
      //cout<<"Time evolution step nr "<<cnt+1<<", time="<<time<<", value "<<M1ev<<", trace "<<trR<<endl;
      cnt++;cntLoop++;time+=delta;
    }
    //    rho0.gaugeCond('R',1);
    Eev=contractor.contract(rho0,hamil,Id);
    trRho=contractor.contract(rho0,Id);
    //    S2R=compute2RenyiS(rho0,trRho);
    //cout<<"Time evolution step nr "<<cnt<<", time="<<time
    //	<<", trace="<<trRho<<", <H>="<<Eev
    //	<<endl;
    out=new ofstream(outfname,ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"
	<<real(Eev)<<"\t"<<imag(Eev)<<"\t"<<real(trRho)<<"\t"<<imag(trRho)<<"\t";
    computeEntropies(rho0,trRho,out,L0);
    *out<<endl;
    out->close();
    delete out;
  }
  rho0.exportMPS(mpsfile);
}


void constructLocal(int L,int pos,mwArray op,int d,MPS& result){
   mwArray sig0=identityMatrix(d);
   sig0.reshape(Indices(d*d,1,1));
   mwArray sigZ(op);sigZ.reshape(Indices(d*d,1,1));
   result=MPS(L,1,d*d);
   for(int k=0;k<L;k++){
     result.setA(k,sig0);
   }
   result.setA(pos,sigZ);
}


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

void setSingleInitState(MPS& rho0,int L,int d,int initSt,int pos){
  complex_t Xop[]={ZERO_c,ONE_c,ONE_c,ZERO_c}; 
  complex_t Yop[]={ZERO_c,I_c,-1.*I_c,ZERO_c}; 
  complex_t Zop[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c}; 
  mwArray X(Indices(d*d,1,1),Xop);
  mwArray Y(Indices(d*d,1,1),Yop);
  mwArray Z(Indices(d*d,1,1),Zop);
  switch(initSt){
  case 1:
    { cout<<"Setting init State to X"<<endl;
      rho0.setA(pos,X);
      break;}
  case 2:
    {  cout<<"Setting init State to Y"<<endl;
      rho0.setA(pos,Y);
      break;}
  case 3:
    {  cout<<"Setting init State to Z"<<endl;
      rho0.setA(pos,Z);
      break;}
  default:{
    cout<<"Error! Unknown initSt "<<initSt<<endl;
    exit(1);
  }
  }
}

void setExtendedInitState(MPS& rho0,int L,int d,int initSt,int L0,int pos0){
   mwArray sig0=identityMatrix(d);
   complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
   complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
   complex_t dataY[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
   mwArray sigZ(Indices(d,d),dataZ);
   mwArray sigX(Indices(d,d),dataX);
   mwArray sigY(Indices(d,d),dataY);

   mwArray state;
   switch(initSt){
   case 1: // X+
     {
     state=.5*(sig0+sigX);
     break;
     }
   case 2: // Y+
     {
     state=.5*(sig0+sigY);
     break;
     }
   case 3: // Z+
     {
     state=.5*(sig0+sigZ);
     break;
     }
   default:
     cout<<"Error: unknown type onf intSt="<<initSt<<endl;
     exit(1);
   }
   sig0.reshape(Indices(d*d,1,1));
   state.reshape(Indices(d*d,1,1));

   // Set the product of identities on the edges and all ancillas
   for(int k=0;k<pos0;k++){
     rho0.replaceSite(2*k,sig0,false);
     rho0.replaceSite(2*k+1,sig0,false);
   }
   for(int k=pos0+L0;k<L;k++){
     rho0.replaceSite(2*k,sig0,false);
     rho0.replaceSite(2*k+1,sig0,false);
   }
   for(int k=pos0;k<pos0+L0;k++){
     rho0.replaceSite(2*k,state,false);
     rho0.replaceSite(2*k+1,sig0,false);
   }
}

void doubleMPO(const MPO& simpleMPO,MPO& doubleMPO){
  int L=simpleMPO.getLength();
  doubleMPO.initLength(L); // just in case
  for(int k=0;k<L;k++){
    mwArray aux=simpleMPO.getOp(k).getFullData();
    doubleMPO.setOp(k,new DoubleOperator(aux,permute(aux,Indices(3,2,1,4))),true);
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

void extendTransposeMPO(const MPO& simpleMPO,MPO& doubleMPO,int dimA){
  int L=simpleMPO.getLength();
  doubleMPO.initLength(L); // just in case
  mwArray identPhys=identityMatrix(dimA);
  identPhys.reshape(Indices(dimA,1,dimA,1)); 
  for(int k=0;k<L;k++){  
    doubleMPO.setOp(k,new DoubleOperator(identPhys,permute(simpleMPO.getOp(k).getFullData(),Indices(3,2,1,4))),true);
  }
}

void computeLocalOp(const MPS& rho,int pos,int d,ofstream* out){
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d*d,1,1),dataX);
  complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigY(Indices(d*d,1,1),dataY);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d*d,1,1),dataZ);

  // cout<<"Computing local expectation values "<<endl;
  mwArray sig0=identityMatrix(d);
  sig0.reshape(Indices(d*d,1,1));
  int L=rho.getLength();
  MPS oper(L,1,d*d);
  for(int k=0;k<L;k++){
    oper.setA(k,sig0);
  }

  Contractor& contractor=Contractor::theContractor();
  complex_t resL;
  oper.setA(pos,sigX);
  resL=contractor.contract(rho,oper);
  *out<<real(resL)<<"\t"<<imag(resL)<<"\t";
  oper.setA(pos,sigY);
  resL=contractor.contract(rho,oper);
  *out<<real(resL)<<"\t"<<imag(resL)<<"\t";
  oper.setA(pos,sigZ);
  resL=contractor.contract(rho,oper);
  *out<<real(resL)<<"\t"<<imag(resL)<<"\t";
}






void constructLocal(int L,int pos,int initSt,int d,MPS& result){
   mwArray sig0=identityMatrix(d);
   sig0.reshape(Indices(d*d,1,1));
  complex_t Xop[]={ZERO_c,ONE_c,ONE_c,ZERO_c}; 
  complex_t Yop[]={ZERO_c,I_c,-1.*I_c,ZERO_c}; 
  complex_t Zop[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c}; 
  mwArray X(Indices(d*d,1,1),Xop);
  mwArray Y(Indices(d*d,1,1),Yop);
  mwArray Z(Indices(d*d,1,1),Zop);
  mwArray* op;
  switch(initSt){
  case 1:
    { op=&X;
      break;}
  case 2:
    { op=&Y;
      break;}
  case 3:
    { op=&Z;
      break;}
  default:{
    cout<<"Error! Unknown initSt "<<initSt<<endl;
    exit(1);
  }
  }
  op->reshape(Indices(d*d,1,1));
  result=MPS(L,1,d*d);
  for(int k=0;k<L;k++){
    result.setA(k,sig0);
  }
  result.setA(pos,*op);
}



void MPOfromMPS(const MPS& mps,MPO& mpo,bool up){
  int L=mps.getLength();

  mpo.initLength(L);

  for(int k=0;k<L;k++){
    mwArray aux=mps.getA(k).getA();
    Indices dims=aux.getDimensions();
    int d0=sqrt(dims[0]);
    if(d0*d0!=dims[0]){
      cout<<"Error: Dimension of site "<<k<<" ("<<dims[0]
	  <<") does not seem a square=> don't know how t divide it"<<endl;
      exit(1);
    }
    aux.reshape(Indices(d0,d0,dims[1],dims[2]));
    if(up)
      aux.permute(Indices(1,3,2,4));
    else
      aux.permute(Indices(2,3,1,4));
    mpo.setOp(k,new Operator(aux),true);
  }

}

double traceNorm(const MPS& rho){
  MPO rhoMPO(rho.getLength());
  MPOfromMPS(rho,rhoMPO);
  // if(rho.getLength()==4){
  //   rhoMPO.exportForMatlab("rhoTest4.m");
  //   exit(1);
  // }
  mwArray Op;
  //cout<<"Computing trace norm of mps, length "<<rho.getLength()<<endl;
  expandOper(rhoMPO,Op);
  // cout<<"Operator expanded to dims "<<Op.getDimensions()<<endl;
  // trace norm, is the sum of singular values
  mwArray U,S,Vdag;
  wrapper::svd(Op,U,S,Vdag);
  return real(S.trace());
}

int d=2;


/** Compute the 2-Renyi entropies of two pieces of the chain, namely
    A, defined to be from pos1A until pos2A (both included and
    numbered from 0 to L) and the complementary, B */
void compute2RenyiEntropiesAB(const MPS& rho,complex_t trR,int pos1A,int pos2A,double* SA,double*SB){
  int L=rho.getLength();
  // First trace out the complementary of A
  MPS rhoA(rho);
  vector<int> toTrace;
  for(int k=0;k<L;k++){
    if(k<pos1A||k>pos2A)
      toTrace.push_back(k);
  }
  rhoA.traceOutSites(toTrace);
  // Now compute its 2-Renyi entropy
  *SA=compute2RenyiS(rhoA,trR);
  // And now the same for the complementary
  toTrace.clear();
  for(int k=pos1A;k<=pos2A;k++)
    toTrace.push_back(k);
  rhoA=rho;
  rhoA.traceOutSites(toTrace);
  // Now compute its 2-Renyi entropy
  *SB=compute2RenyiS(rhoA,trR);  
}


double compute2RenyiS(const MPS& redRho,complex_t trR){
  // I need trace(rho^2)
  
  // If rho was Hermitian, as it ideally should,
  // trace(rho^2)=trace(rho+ rho) We will have problably some
  // numerical error, so that rho has an antihermitian part (hopefully
  // small). Then trace(rho+ rho) is second order in that error, while
  // trace(rho^2) is first order. With a bit more work, I could get
  // trace(rho^2) too, and then the trace of the hermitian part
  // squared is exactly .5*(Re[tr(rho^2)]+tr(rho+ rho))
  Contractor& contractor=Contractor::theContractor();
  complex_t trR2_c=contractor.contract(redRho,redRho);
  // check for real character?
  double trR2=real(trR2_c);
  if(abs(imag(trR2_c))>1E-10&&abs(imag(trR2_c))>1E-8*abs(trR2)){
    cout<<"WARNING (or error): the trace of rho^+ rho does not seem to be real!"<<endl;
    exit(1);
  }
  double normRho=real(conjugate(trR)*trR);
  return -log2(trR2/normRho); // log2(trace(rho^2))/(1-2)
  
  // For trace(rho^2) I construct an MPO that swaps up and doen
  // indices in the rho MPS, and compute the contraction of this in
  // between rho and the bra for its complex conjugate (to compensate
  // the conjugation that happens within contract)
  // static Operator swap;
  // static bool init=false;
  // if(!init){
  //   mwArray id2=identityMatrix(d*d);
  //   id2.reshape(Indices(d*d,d,d));
  //   id2.permute(Indices(1,3,2));
  //   id2.reshape(Indices(d*d,1,d*d,1));
  //   swap=Operator(id2);
  //   init=true;
  // }
  // int L=redRho.getL();
  // MPO aux(L);
  // for(int k=0;k<L;k++)
  //   aux.setOp(k,&swap,false);

}

void computeEntropies(const MPS& rho,complex_t trR,ofstream* out,int L0,bool header){
  *out<<setprecision(15);
  int L2=rho.getLength(); // 2L
  int L=L2/2; 
  double SA,SB;int pos1,pos2;
  // Now compute a series of predetermined partitions:
  vector<int> toTrace;
  MPS rho1(rho);double S2R;
  if(!header){
    // trace the ancillas first
    for(int k=0;k<L;k++)
      toTrace.push_back(2*k+1);
    rho1.traceOutSites(toTrace);
    S2R=compute2RenyiS(rho1,trR);
    *out<<S2R<<"\t";
  }
  else{    
    *out<<"S2_1(rho)\t";
  }
  // 1) central L0 (first chain) vs rest
  pos1=(L-L0)/2;pos2=pos1+L0-1;
  if(!header){
    compute2RenyiEntropiesAB(rho1,trR,pos1,pos2,&SA,&SB);
    *out<<SA<<"\t"<<SB<<"\t";
  }
  else *out<<"S2[A_1("<<pos1<<"-"<<pos2<<")]\t S2[(L-A)_1] \t";
  // 2) central 2*L0 vs rest
  pos1=(L-2*L0)/2;pos2=pos1+2*L0-1;
  if(!header){
    compute2RenyiEntropiesAB(rho1,trR,pos1,pos2,&SA,&SB);
    *out<<SA<<"\t"<<SB<<"\t";
  }
  else *out<<"S2[A_1("<<pos1<<"-"<<pos2<<")]\t S2[(L-A)_1] \t";
  // 3) left L/2 vs right (ancillas not included!)
  pos1=0;pos2=L/2-1;
  if(!header){
    compute2RenyiEntropiesAB(rho1,trR,pos1,pos2,&SA,&SB);
    *out<<SA<<"\t"<<SB<<"\t";
  }
  else *out<<"S2[A_1("<<pos1<<"-"<<pos2<<")]\t S2[(L-A)_1] \t";
  // 4) central L/2 vs rest (ancillas not included)
  pos1=(L-L/2)/2;pos2=pos1+L/2-1;
  if(!header){
    compute2RenyiEntropiesAB(rho1,trR,pos1,pos2,&SA,&SB);
    *out<<SA<<"\t"<<SB<<"\t";
  }
  else *out<<"S2[A_1("<<pos1<<"-"<<pos2<<")]\t S2[(L-A)_1] \t";
  // Now the same with ancillas also in the A and B subsystems
  if(!header){
    S2R=compute2RenyiS(rho,trR);
    *out<<S2R<<"\t";
  }
  else{    
    *out<<"S2(rho)\t";
  }  
  // 1) central L0 (and its ancillas) vs rest
  pos1=(L-L0)/2;pos2=pos1+L0-1;
  if(!header){
    compute2RenyiEntropiesAB(rho,trR,2*pos1,2*pos2+1,&SA,&SB);
    *out<<SA<<"\t"<<SB<<"\t";
  }
  else *out<<"S2[A("<<pos1<<"-"<<pos2<<")]\t S2[(L-A)] \t";
  // 2) central 2*L0 vs rest
  pos1=(L-2*L0)/2;pos2=pos1+2*L0-1;
  if(!header){
    compute2RenyiEntropiesAB(rho,trR,2*pos1,2*pos2+1,&SA,&SB);
    *out<<SA<<"\t"<<SB<<"\t";
  }
  else *out<<"S2[A("<<pos1<<"-"<<pos2<<")]\t S2[(L-A)] \t";
  // 3) left L/2 vs right
  pos1=0;pos2=L/2-1;
  if(!header){
    compute2RenyiEntropiesAB(rho,trR,2*pos1,2*pos2+1,&SA,&SB);
    *out<<SA<<"\t"<<SB<<"\t";
  }
  else *out<<"S2[A("<<pos1<<"-"<<pos2<<")]\t S2[(L-A)] \t";
  // 4) central L/2 vs rest (ancillas not included)
  pos1=(L-L/2)/2;pos2=pos1+L/2-1;
  if(!header){
    compute2RenyiEntropiesAB(rho,trR,2*pos1,2*pos2+1,&SA,&SB);
    *out<<SA<<"\t"<<SB<<"\t";
  }
  else *out<<"S2[A("<<pos1<<"-"<<pos2<<")]\t S2[(L-A)] \t";

}

