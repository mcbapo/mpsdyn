
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "FoldedOperator.h"
#include "DoubleOperator.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

//#include "SchwingerHamiltonianSz.h"
#include "SchwingerHamiltonian.h"

using namespace shrt;

/** schwingerThermUnfold tries to find the thermal state of the
    Schwinger Hamiltonian \ref <SchwingerHamiltonianSz> using
    purification and imaginary time evolution. 
    Unless in schwingerThermUnfold, sites are double here, 
    without unfolding the system-ancilla dof.
    With this program, we may run Trotter or Taylor expansion.

    Receives arguments:
    \param <N> (int) length of the chain
    \param <mg> (double) parameter \f$m/g\f$ of the Hamiltonian
    \param <x> (double) parameter \f$x\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <beta> (double) inverse temperature
    \param <delta> (double) width of time step (actually, this is the 
                   coeff. in the exponential, and the time step turns 
		   out to be 2*delta, because the answer is rho^+ rho
    \param <stepBlk> (int) save results every stepBlk steps (means
            every delta*stepBlk*2 (4) in beta)
    \param <mode> (int) running Trotter(1) or Taylor(0) or exact exponential for HL (2)
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)
    \param <mpsfile> (char*) optional argument, with a file name where to store the 
                     final MPS.

*/

/** Basic functions to transform normal MPOs onto MPOs which act on
    the special MPS for mixed states, i.e. with squared physical
    dimension.  */
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
/** Analogous to the one ebove, but puts the transposed MPO onto the ancillary indices*/
void extendTransposeMPO(const MPO& simpleMPO,MPO& doubleMPO,int dimA);

/** Compute overlap with GS (<GS|rho|GS>) */
complex_t computeEV(const MPS& rho,const MPS& GS);


/** Depending on the method (Trotter/Taylor) a different series of operators is prepared 
   and returned as an array of MPO pointers. */
void prepareDoubleOperatorsTrotter(const SchwingerHamiltonian& hSchP,complex_t deltaC,vector<MPO*>& ops);
void prepareDoubleOperatorsTaylor(const SchwingerHamiltonian& hSchP,vector<MPO*>& ops);
void prepareDoubleOperatorsExpL(const SchwingerHamiltonian& hSchP,complex_t deltaC,vector<MPO*>& ops);

/** The MPOs constructed in prepareDoubleOperators are applied to simulate a certain number of steps of 
    imaginary time evolution */
void applyStepsTrotter(const vector<MPO*>& ops,MPS& state,int D,int nrSteps);
void applyStepsTrotterExpL(const vector<MPO*>& ops,MPS& state,int D,int nrSteps);
void applyStepsTaylor(const vector<MPO*>& ops,MPS& state,int D,double delta,int nrSteps,int order);

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int N=atoi(argv[++cntr]);
  double mg=atof(argv[++cntr]);
  double x=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  double beta=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int stepBlk=atoi(argv[++cntr]);
  int mode=atoi(argv[++cntr]);
  bool trotter=mode!=0;  // 1 iff the arg!=0
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);
  bool saveFinal=false;
  string mpsFile;
  if(argc>11){
    saveFinal=true;
    mpsFile=string(argv[++cntr]);
  }

  double mu=2*mg*sqrt(x);

  double alpha=0.;
  double L0P=0.;
  //double L0M=-1.;

  // Other parameters in the Hamiltonian
  double gweight=1.; // weight of gauge coupling (set to 0 for non-interacting case)
  double nu=0.; // chemical potential
  double zpen=0.; // penalty term for the Sz!=0 states
  double offset=0.; // global offset in energy
  double scale=1.;
  

  cout<<"Initialized arguments: L="<<N
      <<", mu="<<mu
      <<", x="<<x
      <<", Tmin="<<1/beta
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% N\t x\t mg\t D\t delta\t nrsteps\t beta\t T\t Energy"
	<<endl;
    *out<<setprecision(10);
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);

  int M=beta/(2*delta);

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d0=2;
  MPS thS(N,1,d0*d0);
  thS.setProductState(p_maxent);
  cout<<"Initialized identity state. Check norm (should be 2^N="<<pow(2,N)<<")="
      <<contractor.contract(thS,thS)<<endl;
  // I keep it for later (identity contraction)
  MPS idState(thS);
 
  //SchwingerHamiltonianSz hSchP(N,mu,x,L0P,alpha,zpen,offset,nu,gweight,scale);
  SchwingerHamiltonian hSchP(N,mu,x,L0P,alpha);
  const MPO& hamilP=hSchP.getHMPO();
  cout<<"Initialized basic Hamiltonian"<<endl;
  MPO condensate(N);
  hSchP.constructCondensateMPO(condensate);

  // just to check, values in the GS
  double E0=0.;complex_t Sigma0=ZERO_c; 
  // For tests, I can also keep the real GS and compute the overlap!
  MPS auxGS(N,D,d0);
  contractor.findGroundState(hamilP,D,&E0,auxGS);
  Sigma0=contractor.contract(auxGS,condensate,auxGS);
  cout<<"In the GS (computed with D="<<D<<") the Energy is "<<E0<<" and the condensate is "<<Sigma0<<endl;
  *out<<N<<"\t"<<x<<"\t"<<mg<<"\t"<<D<<"\t"
      <<delta<<"\t inf \t inf \t 0 \t"
      <<E0<<"\t"<<real(Sigma0)
      <<"\t 1"<<endl;
  

  // double MPO for the H and the condensate to be applied on rho
  MPO doubleH(N),doubleSigma(N);
  extendMPO(hamilP,doubleH,d0);
  extendMPO(condensate,doubleSigma,d0);

  vector<MPO*> ops;
  complex_t deltaC=(complex_t){-delta,0.}; // delta  as complex
  if(mode==1)
    prepareDoubleOperatorsTrotter(hSchP,deltaC,ops);
  else if(mode==0)
    prepareDoubleOperatorsTaylor(hSchP,ops);
  else if(mode==2)
    prepareDoubleOperatorsExpL(hSchP,deltaC,ops);


  //exit(1);

  cout<<"Before applying any evolution step, energy is E(Id)="
      <<contractor.contract(thS,doubleH,thS)/contractor.contract(thS,thS)
      <<" and condensate is Sigma(Id)="
      <<contractor.contract(thS,doubleSigma,thS)/contractor.contract(thS,thS)
      <<endl;

  // Now apply the evolution
  int leftSteps=M; // nr of steps left to apply
  double betai=0; // current value of beta
  int cnt=0;
    double norm2=real(contractor.contract(thS,thS)); //trace of rho^2
    complex_t energy2=contractor.contract(thS,doubleH,thS);
    complex_t cond2=contractor.contract(thS,doubleSigma,thS);
    complex_t overl1=computeEV(thS,auxGS);
    *out<<N<<"\t"<<x<<"\t"<<mg<<"\t"<<D<<"\t"
	<<delta<<"\t"<<2*cnt<<"\t"<<2*betai<<"\t"<<.5/betai<<"\t"
	<<real(energy2)/norm2<<"\t"<<real(cond2)/norm2
	<<"\t"<<overl1<<endl;
  while(leftSteps>0){
    int nrToApply=min(stepBlk,leftSteps);
    if(mode==1)
      applyStepsTrotter(ops,thS,D,nrToApply);
    else if (mode==0)
      applyStepsTaylor(ops,thS,D,delta,nrToApply,1); // Taylor 1st order
    else if (mode==2)
      applyStepsTrotterExpL(ops,thS,D,nrToApply);
    leftSteps-=nrToApply;
    betai+=2*delta*nrToApply;
    cnt+=nrToApply;
    // compute/record observables now
    double norm2=real(contractor.contract(thS,thS)); //trace of rho^2
    complex_t energy2=contractor.contract(thS,doubleH,thS);
    complex_t cond2=contractor.contract(thS,doubleSigma,thS);
    complex_t overl1=computeEV(thS,auxGS);
    *out<<N<<"\t"<<x<<"\t"<<mg<<"\t"<<D<<"\t"
	<<delta<<"\t"<<2*cnt<<"\t"<<2*betai<<"\t"<<.5/betai<<"\t"
	<<real(energy2)/norm2<<"\t"<<real(cond2)/norm2
	<<"\t"<<overl1<<endl;
    cout<<"At the end of step "<<cnt<<", 2*beta="<<2*betai<<", T="<<.5/betai
	<<", energy="<<real(energy2)/norm2
	<<", condensate="<<real(cond2)/norm2
	<<", overlap with GS="<<overl1
	<<endl;
  }
  // Finally, if needed, save the state
  if(saveFinal)
    thS.exportMPS(mpsFile.data());

}

void prepareDoubleOperatorsTrotter(const SchwingerHamiltonian& hSchP,complex_t deltaC,vector<MPO*>& ops){
  complex_t deltaC_2=.5*deltaC;
  int N=hSchP.getN();
  // Get the (simple) exponential MPOs
  MPO _expHxe_2(N),_expHxo(N),_expHL(N),_expHL_2(N);
  hSchP.getExponentialMPOx(_expHxe_2,deltaC_2,true);  
  _expHxe_2.exportForMatlab("Hxe2.m");
  cout<<"Constructed exponential MPO for Hxe_2"<<endl;
  hSchP.getExponentialMPOx(_expHxo,deltaC,false);
  _expHxo.exportForMatlab("Hxo.m");
  cout<<"Constructed exponential MPO for Hxo"<<endl;
  hSchP.getExponentialMPOL(_expHL,deltaC);  
  _expHL.exportForMatlab("HL.m");
  cout<<"Constructed exponential MPO for HL"<<endl;
  hSchP.getExponentialMPOL(_expHL_2,deltaC_2);
  cout<<"Constructed all exponential MPOs"<<endl;
  // Prepare the double ones, which will be returned in the vector
  MPO* expHxe_2=new MPO(N);
  MPO* expHxo=new MPO(N);
  MPO* expHL=new MPO(N);
  MPO* expHL_2=new MPO(N);
  int d=2;
  doubleMPO(_expHxe_2,*expHxe_2);
  doubleMPO(_expHxo,*expHxo);
  doubleMPO(_expHL_2,*expHL_2);
  doubleMPO(_expHL,*expHL);
  // Now put them in the proper vector to return
  ops.push_back(expHL_2);
  ops.push_back(expHxo);
  ops.push_back(expHxe_2);
  ops.push_back(expHL);
}

void applyStepsTrotter(const vector<MPO*>& ops,MPS& state,int D,int nrSteps){
  Contractor& contractor=Contractor::theContractor();
  MPS aux(state);
  contractor.optimize(*ops[0],aux,state,D); // the L*delta/2
  for(int k=0;k<nrSteps;k++){
    contractor.optimize(*ops[2],state,aux,D); // the Hxe*delta/2
    contractor.optimize(*ops[1],aux,state,D); // the Hxo*delta
    contractor.optimize(*ops[2],state,aux,D); // the Hxe*delta/2
    if(k<nrSteps-1)
      contractor.optimize(*ops[3],aux,state,D); // the HL*delta in between steps
    else
      contractor.optimize(*ops[0],aux,state,D); // the L*delta/2 at the end
    // Now I need to normalize
    state.gaugeCond('R',1);
  }

}

void prepareDoubleOperatorsExpL(const SchwingerHamiltonian& hSchP,complex_t deltaC,vector<MPO*>& ops){
  complex_t deltaC_2=.5*deltaC;
  int N=hSchP.getN();
  // Get the (simple) exponential MPOs
  MPO _expHxe_2(N),_expHxo(N),_expHL(N),_expHxo_2(N);
  hSchP.getExponentialMPOx(_expHxe_2,deltaC_2,true);  
  //_expHxe_2.exportForMatlab("Hxe2.m");
  cout<<"Constructed exponential MPO for Hxe_2"<<endl;
  hSchP.getExponentialMPOx(_expHxo,deltaC,false);
  cout<<"Constructed exponential MPO for Hxo"<<endl;
  hSchP.getExponentialMPOx(_expHxo_2,deltaC,false);
  //_expHxo_2.exportForMatlab("Hxo2.m");
  cout<<"Constructed exponential MPO for Hxo_2"<<endl;
  hSchP.getExactExponentialMPOL(_expHL,deltaC);  
  _expHL.exportForMatlab("expHL.m");
  cout<<"Constructed exponential MPO for HL"<<endl;
  cout<<_expHL<<endl;
  // Prepare the double ones, which will be returned in the vector
  MPO* expHxe_2=new MPO(N);
  MPO* expHxo=new MPO(N);
  MPO* expHxo_2=new MPO(N);
  MPO* expHL=new MPO(N);
  int d=2;
  doubleMPO(_expHxe_2,*expHxe_2);
  doubleMPO(_expHxo,*expHxo);
  doubleMPO(_expHxo_2,*expHxo_2);
  doubleMPO(_expHL,*expHL);
  // Now put them in the proper vector to return
  ops.clear();
  ops.push_back(expHxo_2);
  ops.push_back(expHxo);
  ops.push_back(expHxe_2);
  ops.push_back(expHL);
}

void applyStepsTrotterExpL(const vector<MPO*>& ops,MPS& state,int D,int nrSteps){
  Contractor& contractor=Contractor::theContractor();
  MPS aux(state);
  contractor.optimize(*ops[0],aux,state,D); // the Hxo*delta/2
  for(int k=0;k<nrSteps;k++){
    contractor.optimize(*ops[2],state,aux,D); // the Hxe*delta/2
    contractor.optimize(*ops[3],aux,state,D); // the HL*delta
    state.gaugeCond('R',1);
    contractor.optimize(*ops[2],state,aux,D); // the Hxe*delta/2
    if(k<nrSteps-1)
      contractor.optimize(*ops[1],aux,state,D); // the Hxo*delta in between steps
    else
      contractor.optimize(*ops[0],aux,state,D); // the Hxo*delta/2 at the end
    // Now I need to normalize
    state.gaugeCond('R',1);
  }

}


/** The evolution we want to simulate is e^(-delta H) rho e^(-delta H) => e^(-delta [HxId + IdxH^T])|rho>
on the vectorized operator, where H^T is the transposed H. We ocould construct an effective H_eff=HxId+IdxH^T
as a single MPO in two ways:
i) directly summing together the extended H and the corresponding one
for H^T on the second index, and creating a MPO with bond dimension
twice that of the H MPO (i.e. 10) 
ii) constructing the MPO explicitly, such that the bond dimension is
optimal (we would get 8, but would need to do it in the same way as
initH works inside SchwingerHamiltonian, so probably it should be
there).
iii) the fastest way (although it is not clear if it will be better
that the one before) is to just construct two MPOs for each of the
terms in the effective Hamiltonian and to apply them both for every step.
In this program, I just construct the last option. */


void prepareDoubleOperatorsTaylor(const SchwingerHamiltonian& hSchP,vector<MPO*>& ops){
  int d0=2;int N=hSchP.getN();
  const MPO& hamilP=hSchP.getHMPO();
  MPO* Hdoub=new MPO(N);
  MPO* HTdoub=new MPO(N);
  extendMPO(hamilP,*Hdoub,d0);
  extendTransposeMPO(hamilP,*HTdoub,d0);
  ops.clear();
  ops.push_back(Hdoub);
  ops.push_back(HTdoub);
}

/** There would also be other ways of applying the operators (e.g use
    both to optimize Id-delta HxId -delta Idx H^T, as if the full
    H_eff was constructed, but sicne Contractor can only deal with a
    single MPO fr approximateTaylorExponential, I will apply both
    exponentials sequentially. */
void applyStepsTaylor(const vector<MPO*>& ops,MPS& state,int D,double delta,int nrSteps,int order){
  Contractor& contractor=Contractor::theContractor();
  for(int k=0;k<nrSteps;k++){
    MPS aux(state);
    contractor.approximateTaylorExponential(*ops[0],state,aux,delta,order,D,true);
    contractor.approximateTaylorExponential(*ops[1],aux,state,delta,order,D,true);
    state.gaugeCond('R',1); // normalize
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


complex_t computeEV(const MPS& rho,const MPS& GS){
  Contractor& contractor=Contractor::theContractor();
  // transform the MPS from rho in a MPO
  int L=rho.getLength();
  MPO rhoOp(L);int d=2;
  for(int k=0;k<L;k++){
    mwArray aux=rho.getA(k).getA();
    Indices dims=aux.getDimensions();
    d=sqrt(dims[0]);
    int Dl=dims[1];int Dr=dims[2];
    aux.reshape(Indices(d,d,Dl,Dr));
    aux.permute(Indices(1,3,2,4));
    rhoOp.setOp(k,new Operator(aux),true); 
  }
  MPS id(L,1,d*d);
  mwArray idS=identityMatrix(d);
  idS.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++){
    id.setA(k,idS);
  }
  complex_t traz=contractor.contract(rho,id);
  return contractor.contract(GS,rhoOp,GS)/traz;
}
