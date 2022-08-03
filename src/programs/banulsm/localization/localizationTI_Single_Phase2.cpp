#include <math.h>
#include <iomanip>
#include <vector>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "SpinMPO.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "HeisenbergAncillaryNoiseHamiltonian.h"

using namespace shrt;

/** Starting from an MPs produced by localizationTI_Single (and taing basically the same parameters),
    this program computes trace norms of central regions around the intial polarized spin,
    using the evolved MPS produced by the other program.

    Receives arguments:
    \param <L> (int) number of sites (in the original chain)
    \param <isXY> (bool) whether the model to be considered is only XY (no ZZ term) 
    \param <J> (double) parameter \f$J\f$ of the Heisenberg
                        Hamiltonian for the ancillary chain
    \param <B> (double) parameter \f$B\f$ of the interaction
                        between both chains
    \param <initSt> (char) initial state in the middle, can be 1 (X),
                        2(Y) or 3 (Z)
    \param <delta> (double) width of the time step used
    \param <M> (int) number of steps
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the results
                        (will be appended)
    \param <mpsfile> (char*) name of the file from which to read the MPS 
*/


void constructLocal(int L,int pos,int initSt,int d,MPS& result);

/** Vectorize an MPO to make an MPS out of it. If up==true (default)
    the first component of the double physical index will be the one
    which was going up.*/
void MPSfromMPO(const MPO& mpo,MPS& mps,bool up=true);

/** Make an MPO from an MPS with double indices. If up==true (default)
    the first component of the double physical index will be the one
    which is now going up.*/
void MPOfromMPS(const MPS& mps,MPO& mpo,bool up=true);

double traceNorm(const MPS& rho);

int d=2;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  bool isXY=false;
  isXY=atoi(argv[++cntr])!=0;
  double J_=atof(argv[++cntr]);
  double B_=atof(argv[++cntr]);
  int initSt=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  const char* mpsfile=argv[++cntr];

  cout<<"Initialized arguments "
      <<", L="<<L
      <<", J="<<J_
      <<", B="<<B_
      <<", delta="<<delta
      <<", initSt="<<initSt
      <<", M="<<M
      <<", D="<<D
      <<", outfile="<<outfname
      <<", mpsfile="<<mpsfile
      <<endl;

  ofstream* out;
  bool appendingFile=file_exists(outfname);
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
    *out<<"% M\t delta\t t\t D\t tr(rho H)\t tr(rho(t)) \t tr(sigma"<<initSt<<")\t |rho(t)_12|_1"
	<<" \t |rho(t)_8|_1 \t |rho(t)_10s|_1 \t |rho(t)_6s|_1 \t |rho(t)_4s|_1"<<endl;
  }
  *out<<setprecision(15);

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;

  // Also the identity, to compute the trace
  MPS Id(2*L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<2*L;k++)
    Id.setA(k,id2);

  cout<<"Analyzing state in file "<<mpsfile<<endl;
  MPS rho(L,1,1);
  rho.importMPS(mpsfile);
  rho.gaugeCond('R',1);

  complex_t trRho,Eev,sigev;
  {
    // Create the Hamiltonian
    double J1x(1.),J1y(1.),J1z(isXY?0.:1.);
    double J2x(J_),J2y(J_),J2z(isXY?0.:J_);
    HeisenbergAncillaryNoiseHamiltonian hamH(L,J1x,J1y,J1z,J2x,J2y,J2z,B_);
  //    HeisenbergAncillaryNoiseHamiltonian hamH(L,J_,B_); 
    const MPO& hamil=hamH.getHMPO();
    MPS hmps(2*L,1,d*d);
    MPSfromMPO(hamil,hmps);
    Eev=contractor.contract(rho,hmps);
    trRho=contractor.contract(rho,Id);
    constructLocal(2*L,L,initSt,d,hmps);
    sigev=contractor.contract(rho,hmps);
  }
  // Now I need to start tracing out things, to construct reduced operators and to compute their trace norms.
  // If I started from the largest subsystem, I could reuse some results, but it is simpler to go one by one.
  // I should normalize everything with 1/tr(rho)

  vector<int> toTrace;

  // 1) trace all but the central block of 12 spins, including second
  // leg in the ladder
  int cut=3;
  int lastTracedL=max(L/2-cut-1,-1);int firstTracedR=min(L/2+cut,L);
  for(int k=0;k<=lastTracedL;k++){
    toTrace.push_back(2*k);
    toTrace.push_back(2*k+1);
  }
  for(int k=firstTracedR;k<L;k++){
    toTrace.push_back(2*k);
    toTrace.push_back(2*k+1);
  }
  cout<<"To get the central "<<cut*4<<" spins (L="<<L<<"), tracing out "<<toTrace<<endl;
  MPS redRho(rho);
  redRho.traceOutSites(toTrace);
  //   double tN12=traceNorm(redRho);
   double tN12=0;

  // // 2) keep 12 in total  // trace out 4 from each side
  // toTrace.clear();
  // int len=redRho.getLength();
  // lastTracedL=max(len/2-3-1,-1);firstTracedR=min(len/2+3,len);
  // for(int k=0;k<=lastTracedL;k++)
  //   toTrace.push_back(k);
  // for(int k=firstTracedR;k<len;k++)
  //   toTrace.push_back(k);
  // cout<<"To get the central 12 spins (L="<<len<<"), tracing out "<<toTrace<<endl;
  // redRho.traceOutSites(toTrace);
  // double tN12=traceNorm(redRho);

  // 3) keep 8 in total // trace out 2 from each side
  toTrace.clear();
  int len=redRho.getLength();
  lastTracedL=max(len/2-2-1,-1);firstTracedR=min(len/2+2,len);
  //lastTracedL=max(len/2-2,-1);firstTracedR=min(len/2+1,len);
  for(int k=0;k<=lastTracedL;k++){
    toTrace.push_back(2*k);
    toTrace.push_back(2*k+1);
  }
  for(int k=firstTracedR;k<len;k++){
    toTrace.push_back(2*k);
    toTrace.push_back(2*k+1);
  }
  cout<<"To get the central 8 spins (L="<<len<<"), tracing out "<<toTrace<<endl;
  redRho.traceOutSites(toTrace);
  double tN8=traceNorm(redRho);

  // 4) keep 10 spins from first leg
  toTrace.clear();
  lastTracedL=max(L/2-5-1,-1);firstTracedR=min(L/2+5,L);
  for(int k=0;k<=lastTracedL;k++){
    toTrace.push_back(2*k); // also the edges
    toTrace.push_back(2*k+1); // also the edges
  }
  for(int k=lastTracedL+1;k<firstTracedR;k++){
    toTrace.push_back(2*k+1); // all the ancillas
  }
  for(int k=firstTracedR;k<L;k++){
    toTrace.push_back(2*k);
    toTrace.push_back(2*k+1); // also the edges
  }
  cout<<"To get the central 10 spins of first chain (L="<<len
      <<"), tracing out "<<toTrace<<endl;  
  rho.traceOutSites(toTrace);
  double tN10s=traceNorm(rho);

  // 5) keep 6 spins // trace 2 from each side
  toTrace.clear();
  len=rho.getLength();
  lastTracedL=max(len/2-3-1,-1);firstTracedR=min(len/2+3,len);
  for(int k=0;k<=lastTracedL;k++)
    toTrace.push_back(k);
  for(int k=firstTracedR;k<L;k++)
    toTrace.push_back(k);
  cout<<"To get the central 6 spins of first chain (L="<<len
      <<"), tracing out "<<toTrace<<endl;  
  rho.traceOutSites(toTrace);
  double tN6s=traceNorm(rho);


  // 6) keep 4 (2) spins // trace one from each side
  toTrace.clear();
  len=rho.getLength();
  lastTracedL=max(len/2-2-1,-1);firstTracedR=min(len/2+2,len);
  //lastTracedL=max(len/2-2,-1);firstTracedR=min(len/2+1,len);
  for(int k=0;k<=lastTracedL;k++)
    toTrace.push_back(k);
  for(int k=firstTracedR;k<L;k++)
    toTrace.push_back(k);
  cout<<"To get the central 4 spins of first chain (L="<<len
      <<"), tracing out "<<toTrace<<endl;  
  rho.traceOutSites(toTrace);
  double tN4s=traceNorm(rho);

  // And finally write the results
  *out<<M<<"\t"<<delta<<"\t"<<M*delta<<"\t"<<D<<"\t"<<real(Eev)<<"\t"<<imag(Eev)<<"\t"
      <<real(trRho)<<"\t"<<imag(trRho)<<"\t"
      <<real(sigev)<<"\t"<<imag(sigev)<<"\t"
      <<tN12<<"\t"<<tN8<<"\t"
      <<tN10s<<"\t"<<tN6s<<"\t"<<tN4s<<"\t";

  cout<<"Computation finished"<<endl;
  *out<<endl;
  out->close();
  delete out;

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
  cout<<"Computing trace norm of mps, length "<<rho.getLength()<<endl;
  expandOper(rhoMPO,Op);
  cout<<"Operator expanded to dims "<<Op.getDimensions()<<endl;
  // trace norm, is the sum of singular values
  mwArray U,S,Vdag;
  wrapper::svd(Op,U,S,Vdag);
  return real(S.trace());
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
