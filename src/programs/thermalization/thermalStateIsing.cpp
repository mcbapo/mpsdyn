
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "FoldedOperator.h"
#include "DoubleOperator.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

/** Construct the MPO for the thermal state of a finite chain under
    the Ising Hamiltonian using the purification ansatz, i.e., an MPO
    for rho(beta/2) and rho(beta)=rho(beta/2)^dagger rho(beta/2).
    
    Arguments:
    \param <N> (int) length of the chain
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian (transverse field)
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian (parallel field)
    \param <D> (int) bond dimension for the half MPO
    \param <beta> (double) inverse temperature
    \param <delta> (double) width of time step: this is the increment in beta
                   as imaginary time evolution is simulated (meaning that each
		   exponential has coefficient delta/2)
                   The time step turns out to be 2*delta, because the answer 
		   is rho^+ rho
    \param <stepBlk> (int) save results every stepBlk steps (means
            every delta*stepBlk*2 in beta)
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)
    \param <mpsfile> (char*) optional argument, with a file name where to store the 
                     final MPS.


*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int N=atoi(argv[++cntr]);
  double g=atof(argv[++cntr]);
  double h=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  double beta=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int stepBlk=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);
  bool saveFinal=false;
  string mpsFile;
  if(argc>10){
    saveFinal=true;
    mpsFile=string(argv[++cntr]);
  }

  cout<<"Initialized arguments: L="<<N
      <<", g="<<g
      <<", h="<<h
      <<", Tmin="<<1/beta
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% N\t g\t h\t D\t delta\t nrsteps\t beta\t T\t Energy/N\t <h(L/2)>"
	<<endl;
    *out<<setprecision(12);
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);

  int M=beta/(2*delta); // nr of steps I need in each beta/2
  // Maybe I need an extra step, if beta/2 is not divisible by delta
  double remainder=beta/2-M*delta;
  // I will need a step of remainder on each half of the state


  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d0=2;
  MPS thS(N,1,d0*d0);
  thS.setProductState(p_maxent); // initial identity, with funny normalization
  cout<<"Initialized identity state. Check norm (unnormalized should be 2^N="<<pow(2,N)<<")="
      <<contractor.contract(thS,thS)<<endl;
  // I keep it for later (identity contraction)
  MPS idState(thS);
 
  IsingHamiltonian hamI(N,d0,-1.,-1.*g,-1.*h); // J=-1;
  const MPO& mpoH=hamI.getHMPO();
  cout<<"Initialized basic Hamiltonian"<<endl;

  // double MPO for the H and the condensate to be applied on rho
  MPO doubleH(N);
  extendMPO(mpoH,doubleH,d0);

  // To compute the local energy density in the center, I use a trick:
  // define a Hamiltonian on only two sites, with half the g and h
  IsingHamiltonian ham12(2,d0,-1.,-.5*g,-.5*h);
  MPO localH(N); Operator* idOp;
  {
    const MPO& h12=ham12.getHMPO();
    // fill in a long MPO with identities
    mwArray id=identityMatrix(d0);id.reshape(Indices(d0,1,d0,1));
    idOp=new Operator(id);
    MPO auxH12(N);
    for(int k=0;k<N/2-1;k++) 
      auxH12.setOp(k,idOp,false);
    auxH12.setOp(N/2-1,&h12.getOp(0),false);
    auxH12.setOp(N/2,&h12.getOp(1),false);
    for(int k=N/2+1;k<N;k++) 
      auxH12.setOp(k,idOp,false);
    extendMPO(auxH12,localH,d0);
  }


  MPO evolU(N); // MPO for imaginary time evolution
  {
    MPO auxU(N);
    hamI.getUMPO(auxU,delta/2,0.,true);
    doubleMPO(auxU,evolU); // the folded e^(-delta/2 H)xe^(-delta/2 H^T)
  }

  cout<<"Before applying any evolution step, energy is E(Id)="
      <<contractor.contract(thS,doubleH,thS)/contractor.contract(thS,thS)
      <<endl;

  // Now apply the evolution
  int leftSteps=M; // nr of steps left to apply
  double betai=0; // current value of beta
  int cnt=0;
  while(leftSteps>0){
    int nrToApply=min(stepBlk,leftSteps);
    for(int k=0;k<nrToApply;k++){
      MPS aux(thS);
      contractor.optimize(evolU,aux,thS,D); 
    }
    leftSteps-=nrToApply;
    betai+=2*delta*nrToApply;
    cnt+=nrToApply;
    // compute/record observables now
    double norm2=real(contractor.contract(thS,thS)); //trace of rho^2
    complex_t energy2=contractor.contract(thS,doubleH,thS);
    complex_t locEnergy=contractor.contract(thS,localH,thS);
    *out<<N<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"
	<<delta<<"\t"<<2*cnt<<"\t"<<betai<<"\t"<<1./betai<<"\t"
	<<real(energy2)/norm2<<"\t"
	<<real(locEnergy)/norm2
	<<endl;
    cout<<"At the end of step "<<cnt<<", beta="<<betai<<", T="<<1./betai
	<<", energy="<<real(energy2)/norm2<<", local energy="<<real(locEnergy)/norm2
	<<endl;
  }
  // If I am not exactly at the desired beta, I have to apply an extra step
  if(remainder!=0){
    MPO auxU(N);
    hamI.getUMPO(auxU,remainder/2,0.,true);
    doubleMPO(auxU,evolU); // the folded e^(-delta/2 H)xe^(-delta/2 H^T)
    MPS aux(thS);
    contractor.optimize(evolU,aux,thS,D); 
    betai+=remainder*2;cnt++;
    // compute/record observables now
    double norm2=real(contractor.contract(thS,thS)); //trace of rho^2
    complex_t energy2=contractor.contract(thS,doubleH,thS);
    complex_t locEnergy=contractor.contract(thS,localH,thS);
    *out<<N<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"
	<<delta<<"\t"<<2*cnt<<"\t"<<betai<<"\t"<<1./betai<<"\t"
	<<real(energy2)/norm2<<"\t"
	<<real(locEnergy)/norm2
	<<endl;
    cout<<"At the end of step "<<cnt<<", beta="<<betai<<", T="<<1./betai
	<<", energy="<<real(energy2)/norm2<<", local energy="<<real(locEnergy)/norm2
	<<endl;
  }
  // Finally, if needed, save the state
  if(saveFinal)
    thS.exportMPS(mpsFile.data());

}
