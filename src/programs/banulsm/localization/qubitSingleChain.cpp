#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "HeisenbergHamiltonian.h"

using namespace shrt;

/** qubitSingleChain starts or continues the evolution using for the
    Hamiltonian the random parameters contained in the mwArray file
    provided as argument.
    It thus simulates the time evolution of an MPO under a
    disordered Heisenberg Hamiltonian \ref <HeisenbergHamiltonian>.
    A number of sites in the center of the chain are used to encode a single qubit, so that 
    \f$| \tilde{0} \rangle \langle \tilde{0}|\f$ corresponds to 
    \f$\left ( | 0 \rangle \langle 0| \right )^{\otimes L_0}\f$
    tensored with identities on both sides.
    Starting from an MPO which encodes the \f$sigma_{x(y,z)}\f$, 
    the program simulates the evolution and computes the trace norm of the resulting state 
    when tracing out the edges, and keeping \f$L_0,L0+2,\ldots 10\f$
    spins. This will allow us to compute the recovery fidelity.
    computes the average and localized magnetization along the way.
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
    \param <init> (char) initial state in the middle, can be 1 (X),
                        2(Y) or 3 (Z)
    \param <delta> (double) width of the time step used
    \param <M> (int) number of steps
    \param <rate> (int) frequency of recording data (nr of steps)
                       Results will be saved for 0,rate,2*rate...,M
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the results
                   (will be appended if newInstance is 0 and it exists)
    \param <mpsfile> (char*) name of the file where to read/save the MPS 
    \param <jobsdir> (char*) name of the dir to write the new jobs
    \param <outfname_traces> (char*) name of the output file for the trace norms results
    \param <normFreq> (int) frequency (steps) to compute and write trace norms of reduced states.
                    If it is 0, non will be computed.
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
  const char* jobsdir=argv[++cntr];
  const char* outfname_Traces=argv[++cntr];
  int normFreq=atoi(argv[++cntr]);
  int newInstance=atoi(argv[++cntr]);
  if(newInstance<0) newInstance=1;
  int M0=0;
  if(!newInstance&&argc>12) M0=atoi(argv[++cntr]);
    

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
      <<", normFreq="<<normFreq
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
    *out<<"% M\t delta\t t\t D\t tr(rho H)\t tr(rho) \t <X>\t <Y>\t <Z>"<<endl;
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
  MPO expHe2(L),expHe(L),expHo(L);
  // I construct directly the double operators
  hamH.getDoubleExponentialMPOeven(expHe2,-delta*.5*I_c);
  hamH.getDoubleExponentialMPOeven(expHe,-delta*I_c);
  hamH.getDoubleExponentialMPOodd(expHo,-delta*I_c);

  // // I construct the simple ones first, and then from them the DoubleOperators
  //   {
  // Regular ones
  // MPO _expHe2(L),_expHe(L),_expHo(L);
  //   MPO _expHe2A(L),_expHeA(L),_expHoA(L); // for the second indices
  //     hamH.getExponentialMPOeven(_expHe2,-delta*.5*I_c);
  // hamH.getExponentialMPOeven(_expHe,-delta*I_c);
  // hamH.getExponentialMPOodd(_expHo,-delta*I_c);
  // // _expHe.exportForMatlab("expHe.m");
  // // _expHo.exportForMatlab("expHo.m");
  //   hamH.getExponentialMPOeven(_expHe2A,delta*.5*I_c);
  //   hamH.getExponentialMPOeven(_expHeA,delta*I_c);
  //   hamH.getExponentialMPOodd(_expHoA,delta*I_c);
  //   for(int k=0;k<L;k++){
  //     expHe2.setOp(k,new DoubleOperator(_expHe2.getOp(k).getFullData(),
  // 					permute(_expHe2A.getOp(k).getFullData(),Indices(3,2,1,4))),true);
  //     expHe.setOp(k,new DoubleOperator(_expHe.getOp(k).getFullData(),
  // 				       permute(_expHeA.getOp(k).getFullData(),Indices(3,2,1,4))),true);
  //     expHo.setOp(k,new DoubleOperator(_expHo.getOp(k).getFullData(),
  // 				       permute(_expHoA.getOp(k).getFullData(),Indices(3,2,1,4))),true);
  //   }
  //   // The simple ones will be destroyed now. I already have them copied
  // }
  
 
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
  //MPO auxRho(L);
  //MPOfromMPS(rho0,auxRho);
  //auxRho.exportForMatlab("rho0.m");
  //exit(1);

  // Also the identity, to compute the trace
  MPS Id(L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++)
    Id.setA(k,id2);
  
  MPS hmps(L,1,d*d);
  const MPO& hamil=hamH.getHMPO();
  MPSfromMPO(hamil,hmps);
  
  //  M1.exportForMatlab("M1.m");
  //Id.exportForMatlab("Id.m");
  
  // Now do the evolution step by step, and compute the expectation value of M1 every time
  int cnt=M0;
  int posCenter=L/2; // for local opers
  double time=delta*cnt;
  complex_t trR=contractor.contract(rho0,Id);
  complex_t Eev=contractor.contract(rho0,hmps);
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"
	<<real(Eev)<<"\t"<<imag(Eev)<<"\t"<<real(trR)<<"\t"<<imag(trR)<<"\t";
  computeLocalOp(rho0,posCenter,d,out);
  *out<<endl;
  out->close();delete out;
  computeTraceNorms(rho0,outfname_Traces,L,L0,isXY,J_,h_,initSt,delta,cnt,D);
  while(cnt<Mtot){
    //cout<<"Time evolution step nr "<<cnt<<", time="<<time<<endl;
    //", value "<<M1ev<<", trace "<<trR<<endl;
    MPS aux(rho0); // temporary copy
    contractor.optimize(expHe2,aux,rho0,D);
    // cout<<"After expHe2, v'*H*v="<<contractor.contract(imaGS,hamil,imaGS)
    // 	   <<", v*v="<<contractor.contract(imaGS,imaGS)<<endl;
    // Now apply r-1 times the pair Ho He
    int cntLoop=0;
    while(cntLoop<rate&&cnt<Mtot){
      contractor.optimize(expHo,rho0,aux,D);
      // cout<<"After expHo, v'*H*v="<<contractor.contract(imaGS,hamil,imaGS)
      // 	   <<", v*v="<<contractor.contract(imaGS,imaGS)<<endl;
      if(cntLoop<rate-1&&cnt<Mtot-1){
	contractor.optimize(expHe,aux,rho0,D);
      }
      else{
	contractor.optimize(expHe2,aux,rho0,D);
      }
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
    computeLocalOp(rho0,posCenter,d,out);
    *out<<endl;
    out->close();delete out;
    // export the MPS if cnt is multiple of savingFreq
    if(normFreq&&cnt%normFreq==0){
      computeTraceNorms(rho0,outfname_Traces,L,L0,isXY,J_,h_,initSt,delta,cnt,D);
    }
  }

  rho0.exportMPS(mpsfile);

}


// void setExtendedInitState(MPS& rho0,int L,int d,int initSt,int L0,int pos0){
//    mwArray sig0=identityMatrix(d);
//    sig0.reshape(Indices(d*d,1,1));
//    for(int k=0;k<pos0;k++)
//      rho0.replaceSite(k,sig0,false);
//    for(int k=pos0+L0;k<L;k++)
//      rho0.replaceSite(k,sig0,false);

//   int ind01=2; //complex_t op01[]={ZERO_c,ZERO_c,ONE_c,ZERO_c}; 
//   int ind10=1; //complex_t op10[]={ZERO_c,ONE_c,ZERO_c,ZERO_c}; 
//   int ind00=0; //complex_t op00[]={ONE_c,ZERO_c,ZERO_c,ZERO_c}; 
//   int ind11=3; //complex_t op11[]={ZERO_c,ZERO_c,ZERO_c,ONE_c}; 
//   int ind0,ind1;complex_t coeff0,coeff1;
//   switch(initSt){
//   case 1:
//     { cout<<"Setting init State to X, starting on pos "<<pos0<<", finishing "<<pos0+L0-1<<endl;
//       ind0=ind01;ind1=ind10;
//       coeff0=ONE_c;coeff1=ONE_c;
//       break;}
//   case 2:
//     {  cout<<"Setting init State to Y, starting on pos "<<pos0<<", finishing "<<pos0+L0-1<<endl;
//       ind0=ind01;ind1=ind10;
//       coeff0=-1.*I_c;coeff1=I_c;
//       break;}
//   case 3:
//     {  cout<<"Setting init State to Z, starting on pos "<<pos0<<", finishing "<<pos0+L0-1<<endl;
//       ind0=ind00;ind1=ind11;
//       coeff0=ONE_c;coeff1=-1.*ONE_c;
//       break;}
//   default:{
//     cout<<"Error! Unknown initSt "<<initSt<<endl;
//     exit(1);
//   }
//   }
//   // The original MPS is changed s.t. a shorter MPS is inserted (from pos0 to pos0+L0-1) with the right state.
//   //left edge 
//   int D=L0==1?1:2;
//   mwArray Al(Indices(d*d,max(1,rho0.getA(pos0).getDl()),D));Al.fillWithZero();
//   mwArray Am(Indices(d*d,D,D));Am.fillWithZero();
//   mwArray Ar(Indices(d*d,D,max(1,rho0.getA(pos0+L0-1).getDr())));Ar.fillWithZero();
//   Al.setElement(coeff0,Indices(ind0,0,0));Al.setElement(coeff1,Indices(ind1,0,D-1));
//   Am.setElement(ONE_c,Indices(ind0,0,0));Am.setElement(ONE_c,Indices(ind1,D-1,D-1));
//   Ar.setElement(ONE_c,Indices(ind0,0,0));Ar.setElement(ONE_c,Indices(ind1,D-1,0));
//   // And the special identity for the intermediate ancillas:
//   //rho0.setA(pos0,Al);
//   rho0.replaceSite(pos0,Al,false);
//   for(int k=1;k<L0-1;k++){
//     //rho0.setA((pos0+k),Am);
//     rho0.replaceSite((pos0+k),Am,false);
//   }
//   //  rho0.setA(pos0+L0-1,Ar);
//   rho0.replaceSite(pos0+L0-1,Ar,false);
// }



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


void computeTraceNorms(const MPS& rho0,const char* outfname,int L,int L0,bool isXY,double J_,double h_,int initSt,
		       double delta,int M,int D){
  ofstream* out;
  bool appendingFile=file_exists(outfname)&&(M!=0);
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
    *out<<"% L="<<L<<", J="<<J_<<", h="<<h_
	<<" D="<<D<<", init St="<<initSt<<endl;
    *out<<"% M\t delta\t t\t D\t tr(rho H)\t tr(rho(t)) \t tr(sigma"<<initSt<<")\t |rho(t)_2|_1"<<endl;
  }
  *out<<setprecision(15);

  Contractor& contractor=Contractor::theContractor();
  //cout<<"Initialized Contractor"<<endl;
  // Also the identity, to compute the trace
  MPS Id(L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++)
    Id.setA(k,id2);

  MPS rho(rho0); // copying
  rho.gaugeCond('R',1);

  complex_t trRho,Eev,sigev;
  {
    double Jx(J_),Jy(J_),Jz(isXY?0.:J_);
    HeisenbergHamiltonian hamH(L,Jx,Jy,Jz,h_,d);
    const MPO& hamil=hamH.getHMPO();
    MPS hmps(L,1,d*d);
    MPSfromMPO(hamil,hmps);
    Eev=contractor.contract(rho,hmps);
    trRho=contractor.contract(rho,Id);
    constructLocal(L,L/2,initSt,d,hmps);
    sigev=contractor.contract(rho,hmps);
  }

  *out<<M<<"\t"<<delta<<"\t"<<M*delta<<"\t"<<D<<"\t"<<real(Eev)<<"\t"<<imag(Eev)<<"\t"
      <<real(trRho)<<"\t"<<imag(trRho)<<"\t"
      <<real(sigev)<<"\t"<<imag(sigev)<<"\t"; // 10 columns for nothing (3 is time)

  // Now I need to start tracing out things, to construct reduced operators and to compute their trace norms.
  // I will keep 10, 8,... L0+2, L0
  // Column 11 is 10; 12 is 8; 13 is 6; 14 is 4; 15 is 2


  for(int cut=5;cut>=L0/2;cut--){
    vector<int> toTrace;
    // 1) trace all but the central block of 2*cut spins
    int lastTracedL=max(L/2-cut-1,-1);int firstTracedR=min(L/2+cut,L);
    for(int k=0;k<=lastTracedL;k++){
      toTrace.push_back(k);
    }
    for(int k=firstTracedR;k<L;k++){
      toTrace.push_back(k);
    }
    //cout<<"To get the central "<<cut*2<<" spins (L="<<L<<"), tracing out "<<toTrace<<endl;
    MPS redRho(rho);
    redRho.traceOutSites(toTrace);
    double tN10=traceNorm(redRho);
    // And finally write the results, preceeded by the number kept!
    *out<<2*cut<<"\t";
    *out<<tN10<<"\t";
  }
  //cout<<"Computation finished"<<endl;
  *out<<endl;
  out->close();
  delete out;
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
   sig0.reshape(Indices(d*d,1,1));
   complex_t data00[]={ONE_c,ZERO_c,ZERO_c,ZERO_c};
   mwArray op00(Indices(d*d,1,1),data00);
   complex_t data01[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
   mwArray op01(Indices(d*d,1,1),data01);
   complex_t data10[]={ZERO_c,ONE_c,ZERO_c,ZERO_c};
   mwArray op10(Indices(d*d,1,1),data10);
   complex_t data11[]={ZERO_c,ZERO_c,ZERO_c,ONE_c};
   mwArray op11(Indices(d*d,1,1),data11);
   mwArray Z(Indices(4,d*d));
   for(int j=0;j<d*d;j++){
     Z.setElement(op00.getElement(Indices(j,0,0)),Indices(0,j));
     Z.setElement(op01.getElement(Indices(j,0,0)),Indices(1,j));
     Z.setElement(op10.getElement(Indices(j,0,0)),Indices(2,j));
     Z.setElement(op11.getElement(Indices(j,0,0)),Indices(3,j));
   }
   // Set the product of identities on the edges
   for(int k=0;k<pos0;k++)
     rho0.replaceSite(k,sig0,false);
   for(int k=pos0+L0;k<L;k++)
     rho0.replaceSite(k,sig0,false);
   int ind1(0),ind2(3);
   complex_t coeff1(ONE_c),coeff2(-1*ONE_c);
   if(initSt<3){
     ind1=1;ind2=2;
     coeff1=ONE_c;coeff2=ONE_c;
     if(initSt==2){
       coeff1=-1*I_c;coeff2=I_c;
     }
   }
   int D=2; // Dimension of the MPS
   mwArray C(Indices(D,D,4)),C1(Indices(1,D,4)),CL(Indices(D,1,4));
   C.setElement(ONE_c,Indices(0,0,ind1));
   C.setElement(ONE_c,Indices(D-1,D-1,ind2));
   C1.setElement(coeff1,Indices(0,0,ind1));
   C1.setElement(coeff2,Indices(0,D-1,ind2));
   CL.setElement(ONE_c,Indices(0,0,ind1));
   CL.setElement(ONE_c,Indices(D-1,0,ind2));
   
   C1.reshape(Indices(1*D,4));C1.multiplyRight(Z);C1.reshape(Indices(1,D,d*d));C1.permute(Indices(3,1,2));
   rho0.replaceSite(pos0,C1,false);
   C.reshape(Indices(D*D,4));C.multiplyRight(Z);C.reshape(Indices(D,D,d*d));C.permute(Indices(3,1,2));
   for(int k=pos0+1;k<pos0+L0-1;k++)
     rho0.replaceSite(k,C,false);
   CL.reshape(Indices(D*1,4));CL.multiplyRight(Z);CL.reshape(Indices(D,1,d*d));CL.permute(Indices(3,1,2));
   rho0.replaceSite(pos0+L0-1,CL,false);
}
