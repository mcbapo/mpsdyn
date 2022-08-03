#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "HeisenbergHamiltonian.h"

using namespace shrt;

// To activate Trotter order 4 
//#define ORDER4_C 1

/** mutualInfoSingleChain starts or continues the evolution using for the
    Hamiltonian the random parameters contained in the mwArray file
    provided as argument.
    It thus simulates the time evolution of an MPO under a
    disordered Heisenberg Hamiltonian \ref <HeisenbergHamiltonian>. 
    A number of sites in the center of the chain are initially in a 
    well defined product state (|0> or|0+1>)
    tensored with identities on both sides.
    The program simulates the evolution and computes the 2-Renyi 
    entropies of different parts of the chain, to determine the mutual information.
    To keep track of how good the evolution is, I will also compute
    the energy along the way.

    Receives arguments:
    \param <L> (int) number of sites
    \param <L0> (int) number of sites (in the original chain)
    \param <isXY> (bool) whether the model to be considered is only XY (no ZZ term) 
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian (Jx=Jy=Jz or Jx=Jy)
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian (local values 
                        randomly distributed between -h and h)			
    \param <isCont> (bool) whether the model to be considered has continuous randomness
    \param <paramFile> (char*) file with the coefficients for the magnetic field
    \param <nrInst> (int) Nr of row in the file to be used (if rows
                       are longer than L, the first L values of the row will be used).
    \param <init> (char) initial pure state in the middle L0 sites, can be 1 (X+),
                        2(Y+) or 3 (Z+)
    \param <delta> (double) width of the time step used
    \param <M> (int) number of steps
    \param <rate> (int) frequency of recording data (nr of steps)
                       Results will be saved for 0,rate,2*rate...,M
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the results
                   (will be appended if newInstance is 0 and it exists)
    \param <mpsfile> (char*) name of the file where to read/save the final MPS 
    \param <newInstance> (bool) whether to start a new evolution from
                         the inhomogeneous polarization distribution
			 (if not, the MPS in mpsfile is used)
    \param <M0> (int) OPTIONAL, only read if newInstance is false, it
                         is the number of steps already done in the
			 existing file.
*/

/** Construct the initial state and observable, but vectorized*/
//void constructM1(int L,int d,MPS& result);

/** Vectorize an MPO to make an MPS out of it. If up==true (default)
    the first component of the double physical index will be the one
    which was going up.*/
void MPSfromMPO(const MPO& mpo,MPS& mps,bool up=true);
void MPOfromMPS(const MPS& mps,MPO& mpo,bool up=true);

/** Set initial state to initSt (sigma_x,y or z), encoded in L0 spins (from pos0) */
void setExtendedInitState(MPS& rho0,int L,int d,int initSt,int L0,int pos0);


void computeLocalOp(const MPS& rho,int pos,int d,ofstream* out);
void computeTraceNorms(const MPS& rho0,const char* outfname,int L,int L0,bool isXY,double J_,double h_,int initSt,
		       double delta,int M,int D);
void constructLocal(int L,int pos,int initSt,int d,MPS& result);
double traceNorm(const MPS& rho);

/** Compute the 2-Renyi entropies of two pieces of the chain, namely
    A, defined to be from pos1A until pos2A (both included and
    numbered from 0 to L) and the complementary, B */
double compute2RenyiS(const MPS& redRho,complex_t traceR);
void compute2RenyiEntropiesAB(const MPS& rho,complex_t trR,int pos1A,int pos2A,double* SA,double*SB);
void computeEntropies(const MPS& rho,complex_t trR,ofstream* out,int L0,bool header=false);


int d=2;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  int L0=atoi(argv[++cntr]);
  bool isXY=false;
  isXY=atoi(argv[++cntr])!=0;
  double J_=atof(argv[++cntr]);
  double h_=atof(argv[++cntr]);
  bool isCont=false;
  isCont=atoi(argv[++cntr])!=0;
  const char* paramfname=argv[++cntr];
  int nrInst=atoi(argv[++cntr]);
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
  if(!newInstance&&argc>17) M0=atoi(argv[++cntr]);
    

  cout<<"Initialized arguments "
      <<", L="<<L
      <<", L0="<<L0
      <<", isXY?="<<isXY
      <<", isCont?="<<isCont
      <<", J="<<J_
      <<", h="<<h_
      <<", delta="<<delta
      <<", M="<<M
      <<", D="<<D
      <<", outfile="<<outfname
      <<", paramfile="<<paramfname
      <<", nrInst="<<nrInst
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
    *out<<"% L="<<L<<", L0="<<L0<<", J="<<J_<<", h="<<h_
	<<", isXY="<<isXY<<", initSt="<<initSt<<", D="<<D<<endl;
    *out<<"% hs in "<<paramfname<<", row "<<nrInst<<endl;
    *out<<"% MPS in "<<mpsfile<<endl;
  }

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  vector<double> J(L-1,J_);
  vector<double> Jz(L-1,isXY?0.:J_);
  vector<double> h(L,h_);
  
  if(!isCont)
    cout<<"WARNING!! Turning to discrete values of the random fields!"<<endl;

  // Now read the parameters
  {
    mwArray H;
    ifstream inData(paramfname);
    if(!inData.is_open()){
      cout<<"Error: Cannot open the mwArray file "<<paramfname<<" to read data "<<endl;
    }
    H.load(inData);
    double aux_(0.);
    for(int k=0;k<L;k++){
      double haux(0.);
      H.getElement(haux,aux_,Indices(nrInst,k));
      h[k]=h_*haux;
      if(!isCont)
	h[k]=haux>0?h_:-h_;
    }
    // TODO: Check for right dimensions
    inData.close();
  }

  cout<<"Read h="<<h<<endl;  
  if(!appendingFile){
    *out<<"% Random values h:"<<endl;
    *out<<setprecision(18);
    for(int k=0;k<L;k++) *out<<"% h("<<k+1<<")="<<h[k]<<endl;
    *out<<setprecision(15);
    *out<<"% M\t delta\t t\t D\t Re[tr(rho H)]\t Im[tr(rho H)]\t Re[tr(rho)] \tIm[tr(rho)] \t";
    *out<<"S2(rho)\t";
    computeEntropies(MPS(L,1,1),I_c,out,L0,true);
    *out<<endl;
  }
  int Mtot=M0+M;

  // First: create the Hamiltonian
  HeisenbergHamiltonian hamH(L,J,J,Jz,h,d);
  cout<<"Created the Hamiltonian"<<endl;
  //const MPO& hamilMPO=hamH.getHMPO();
  //hamilMPO.exportForMatlab("hamilMPO.m");
  //exit(1);

  // The evolution operator for a rho will result from double operators
  // Which MPOs do I need for the indices in both sides
#ifndef ORDER4_C
  MPO expHe2(L),expHe(L),expHo(L);
  // I construct directly the double operators
  hamH.getDoubleExponentialMPOeven(expHe2,-delta*.5*I_c);
  hamH.getDoubleExponentialMPOeven(expHe,-delta*I_c);
  hamH.getDoubleExponentialMPOodd(expHo,-delta*I_c);
#endif
  // Should I do Trotter 4???
#ifdef ORDER4_C
    cout<<"Running fourth Trotter order with complex parameters!"<<endl;
    complex_t Apar=(.5*ONE_c+(sqrt(3)/6.)*I_c);
    // For Trotter fourth order I need 5 different operators
    MPO expHe_a(L),expHe_aC(L),expHo_aC_2(L),expHo_a_2(L),expHo_2(L);
    complex_t deltaC=-delta*I_c;
    hamH.getDoubleExponentialMPOeven(expHe_a,Apar*deltaC);
    hamH.getDoubleExponentialMPOeven(expHe_aC,conjugate(Apar)*deltaC);
    hamH.getDoubleExponentialMPOodd(expHo_a_2,.5*Apar*deltaC);
    hamH.getDoubleExponentialMPOodd(expHo_2,.5*deltaC);
    hamH.getDoubleExponentialMPOodd(expHo_aC_2,.5*conjugate(Apar)*deltaC);
#endif

  // Now construct the initial rho MPO, in this case, the identity,
  // except for L0 sites 
  MPS rho0(L,1,d*d);
  // Now read the initial rho MPO from the file or construct it for new instance
  if(!newInstance)
    rho0.importMPS(mpsfile);
  else{
    //rho0.setProductState(p_maxent);
    setExtendedInitState(rho0,L,d,initSt,L0,(L-L0)/2);
  }

  // Also the identity, to compute the trace
  MPS Id(L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++)
    Id.setA(k,id2);
  
  MPS hmps(L,1,d*d);
  const MPO& hamil=hamH.getHMPO();
  MPSfromMPO(hamil,hmps);
  
  // Now do the evolution step by step, and compute the entropies
  int cnt=M0;
  double time=delta*cnt;
  complex_t trR=contractor.contract(rho0,Id);
  complex_t Eev=contractor.contract(rho0,hmps);
  double S2R=compute2RenyiS(rho0,trR);
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"
	<<real(Eev)<<"\t"<<imag(Eev)<<"\t"<<real(trR)<<"\t"<<imag(trR)<<"\t";
  *out<<S2R<<"\t";
  computeEntropies(rho0,trR,out,L0);
  *out<<endl;
  out->close();delete out;
  while(cnt<Mtot){
    //cout<<"Time evolution step nr "<<cnt<<", time="<<time<<endl;
    //", value "<<M1ev<<", trace "<<trR<<endl;
    MPS aux(rho0); // temporary copy
#ifdef ORDER4_C
    contractor.optimize(expHo_aC_2,aux,rho0,D);
#endif
#ifndef ORDER4_C
    contractor.optimize(expHe2,aux,rho0,D);
#endif
    // cout<<"After expHe2, v'*H*v="<<contractor.contract(imaGS,hamil,imaGS)
    // 	   <<", v*v="<<contractor.contract(imaGS,imaGS)<<endl;
    // Now apply r-1 times the pair Ho He
    int cntLoop=0;
    while(cntLoop<rate&&cnt<Mtot){
#ifndef ORDER4_C
      contractor.optimize(expHo,rho0,aux,D);
      // cout<<"After expHo, v'*H*v="<<contractor.contract(imaGS,hamil,imaGS)
      // 	   <<", v*v="<<contractor.contract(imaGS,imaGS)<<endl;
      if(cntLoop<rate-1&&cnt<Mtot-1){
	contractor.optimize(expHe,aux,rho0,D);
      }
      else{
	contractor.optimize(expHe2,aux,rho0,D);
      }
#endif
#ifdef ORDER4_C
      contractor.optimize(expHe_aC,rho0,aux,D);
      contractor.optimize(expHo_2,aux,rho0,D);
      contractor.optimize(expHe_a,rho0,aux,D);
      if(cntLoop<rate-1&&cnt<Mtot-1){
	contractor.optimize(expHo_2,aux,rho0,D);
      }
      else{
	contractor.optimize(expHo_a_2,aux,rho0,D);
      }
#endif
      cnt++;cntLoop++;time+=delta;
    }
    // rho0.gaugeCond('R',1);
    // cout<<"After expHe2, v'*H*v="<<contractor.contract(imaGS,hamil,imaGS)
    // 	   <<", v*v="<<contractor.contract(imaGS,imaGS)<<endl;
    //rho0.gaugeCond('R',true);
    trR=contractor.contract(rho0,Id);
    Eev=contractor.contract(rho0,hmps);
    out=new ofstream(outfname,ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"
	<<real(Eev)<<"\t"<<imag(Eev)<<"\t"<<real(trR)<<"\t"<<imag(trR)<<"\t";
    S2R=compute2RenyiS(rho0,trR);
    *out<<S2R<<"\t";
    computeEntropies(rho0,trR,out,L0);
    *out<<endl;
    out->close();delete out;
  }

  rho0.exportMPS(mpsfile);

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
  int L=rho.getLength();
  double SA,SB;int pos1,pos2;
  // Now compute a series of predetermined partitions:
  // 1) central L0 vs rest
  pos1=(L-L0)/2;pos2=pos1+L0-1;
  if(!header){
    compute2RenyiEntropiesAB(rho,trR,pos1,pos2,&SA,&SB);
    *out<<SA<<"\t"<<SB<<"\t";
  }
  else *out<<"S2[A("<<pos1<<"-"<<pos2<<")]\t S2[L-A] \t";
  // 2) central 2*L0 vs rest
  pos1=(L-2*L0)/2;pos2=pos1+2*L0-1;
  if(!header){
    compute2RenyiEntropiesAB(rho,trR,pos1,pos2,&SA,&SB);
    *out<<SA<<"\t"<<SB<<"\t";
  }
  else *out<<"S2[A("<<pos1<<"-"<<pos2<<")]\t S2[L-A] \t";
  // 3) left L/2 vs right
  pos1=0;pos2=L/2-1;
  if(!header){
    compute2RenyiEntropiesAB(rho,trR,pos1,pos2,&SA,&SB);
    *out<<SA<<"\t"<<SB<<"\t";
  }
  else *out<<"S2[A("<<pos1<<"-"<<pos2<<")]\t S2[L-A] \t";
  // 4) central L/2 vs rest
  pos1=(L-L/2)/2;pos2=pos1+L/2-1;
  if(!header){
    compute2RenyiEntropiesAB(rho,trR,pos1,pos2,&SA,&SB);
    *out<<SA<<"\t"<<SB<<"\t";
  }
  else *out<<"S2[A("<<pos1<<"-"<<pos2<<")]\t S2[L-A] \t";
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

// void constructM1(int L,int d,MPS& result){
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

   // Set the product of identities on the edges
   for(int k=0;k<pos0;k++)
     rho0.replaceSite(k,sig0,false);
   for(int k=pos0+L0;k<L;k++)
     rho0.replaceSite(k,sig0,false);
   for(int k=pos0;k<pos0+L0;k++)
     rho0.replaceSite(k,state,false);
}
