
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#include "BoseHubbardHamiltonian.h"
#include "BosonMPO.h"

// To stop right after GS computation using U
//#define ONLYGS 1
// To switch off the potential for the evolution
//#define OPENV 1

/** bosehubWell finds the ground state of the Bose-Hubbard
    Hamiltonian     \ref <BoseHubbardHamiltonian>
    in a certain potential well, with the expression 
    \f$V(x)=V_0 (1-exp(-(x-x_0)^2/a^2)) \f$, with interactions which
    may be different in the three sectors (L,C,R).
    The limit if the three regions is specified by arguments
    \param<xL> and \param<xR> (both between 0 and L-1),
    which are the last sites included in the well.
    Then, it makes it evolve with the Hamiltonian with interaction
    for a certain number of steps.

    Receives arguments:
    \param <L> (int) total length of the chain
    \param <N> (int) cutoff for the occupation number
    \param <xL> (int) site where the well starts
    \param <xR> (int) site where the well ends
    \param <t> (double) parameter \f$t\f$ of the Hamiltonian
    \param <UL> (double) parameter \f$U\f$ of the Hamiltonian
        in the left part of the well
    \param <UC> (double) parameter \f$U\f$ of the Hamiltonian
        in the center of the well
    \param <UR> (double) parameter \f$U\f$ of the Hamiltonian
        in the right part of the well
    \param <mu> (double) parameter \f$\mu\f$ of the Hamiltonian
    \param <V0> (double) parameter \f$V0\f$ of the Hamiltonian 
    \param <x0> (double) parameter \f$x0\f$ of the Hamiltonian 
    \param <a> (double) parameter \f$a\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <delta> (double) time step \f$\delta\f$ 
    \param <M> (int) maximum number of steps
    \param <stepBl> (int) periodicity for saving results
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)

*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  int N=atoi(argv[++cntr]);
  int xL=atoi(argv[++cntr]);
  int xR=atoi(argv[++cntr]);
  double t=atof(argv[++cntr]);
  double UL=atof(argv[++cntr]);
  double UC=atof(argv[++cntr]);
  double UR=atof(argv[++cntr]);
  double mu=atof(argv[++cntr]);
  double V0=atof(argv[++cntr]);
  int x0=atoi(argv[++cntr]);
  double a=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int stepBlk=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% L="<<L<<", N="<<N<<", xL="<<xL<<",xR="<<xR<<endl;
    *out<<"% D="<<D<<", EnergyGS= ";
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d=N+1;
  MPS gsP(L,D,d);
  gsP.setRandomState(); // the initial state, random -> not to the GS!!
  gsP.gaugeCond('R',1);
  gsP.gaugeCond('L',1);
  double lambda=0.;

  double* tis=new double[L];
  double* Uis=new double[L];
  double* Uis0=new double[L];
  double* muis=new double[L];
  double* muis0=new double[L];
  double* Vis=new double[L];
  double* VisE=new double[L];
  for(int k=0;k<L;k++){
    tis[k]=t;
    Uis0[k]=0;
    Uis[k]=k<xL?UL:(k<=xR?UC:UR);
    // another test!
    //    Uis[k]=k<xL?UL:(k<=xR?UC:UR*(1-exp(-(k-x0)*(k-x0)/(a*a))));
#ifdef ONLYGS 
    Uis0[k]=Uis[k]; // to check!
#endif
    muis0[k]=mu;
    muis[k]=0;
    Vis[k]=V0*(1-exp(-(k-x0)*(k-x0)/(a*a)));
#ifdef OPENV
    VisE[k]=0;
#else
    VisE[k]=Vis[k];
#endif
    cout<<"Site "<<k<<", U="<<Uis0[k]<<",V="<<Vis[k]<<endl;
  }

  BoseHubbardHamiltonian hBH0(L,N,tis,Uis0,muis0,Vis);
  const MPO& hamil=hBH0.getHMPO();

  cout<<"Before starting, expectation value is="
      <<contractor.contract(gsP,hamil,gsP)<<", and norm="
      <<contractor.contract(gsP,gsP)<<endl;

  contractor.findGroundState(hamil,D,&lambda,gsP);

  cout<<"Ground state found with eigenvalue "<<lambda<<endl;
  complex_t normGS=contractor.contract(gsP,gsP);
  complex_t energyGS=contractor.contract(gsP,hamil,gsP);
  cout<<"After optimizing, expectation value is="
      <<energyGS<<", and norm="
      <<normGS<<endl;
  *out<<energyGS/normGS<<endl;


  // Construct a MPO for the particle number on a given site
  mwArray identM=identityMatrix(d);
  identM.reshape(shrt::Indices(d,1,d,1));
  Operator ident(identM);
  double* values=new double[d];
  for(int k=0;k<d;k++){
    values[k]=k;
  }
  mwArray nrop=realDiag(d,values); // individual operator number
  nrop.reshape(shrt::Indices(d,1,d,1));
  Operator nrk(nrop);
  // Basic MPO: all Identity
  MPO mpoN(L);
  for(int k=0;k<L;k++)
    mpoN.setOp(k,&ident,0);
  // Now iterate for occupation number
  *out<<"%k\t t(k)\t U(k)\t mu(k)\t V(k)\t <n(k)>"<<endl;
  // and at the same time compute total number of particles
  double Ntot=0.;
  for(int k=0;k<L;k++){
    mpoN.setOp(k,&nrk,0);
    complex_t nrk=contractor.contract(gsP,mpoN,gsP);
    nrk=nrk/normGS;
    *out<<k<<"\t"<<tis[k]<<"\t"<<Uis0[k]<<"\t"<<muis0[k]<<"\t"<<Vis[k]
	<<"\t"<<nrk<<endl;
    // restore the identity
    mpoN.setOp(k,&ident,0);
    Ntot=Ntot+real(nrk);
  }

  cout<<"Total number of particles N="<<Ntot<<endl;

  // Now check with the number operator MPO
  MPO Nmpo(L);
  BosonMPO::getNumberMPO(L,N,Nmpo);
  complex_t totalN=contractor.contract(gsP,Nmpo,gsP);
  totalN=totalN/normGS;
  cout<<"Computed with the MPO, Ntot="<<totalN<<endl;
  complex_t totalN2=contractor.contract2(gsP,Nmpo,gsP);
  totalN2=totalN2/normGS;
  cout<<"Computed with the MPO, <Ntot^2>="<<totalN2<<endl;
  *out<<"% Ntot="<<totalN<<" N^2="<<totalN2<<", Delta N="<<totalN2-totalN*totalN<<endl;
#ifdef ONLYGS 
  exit(1);
#endif

  // MPO parMPO(L);
  // BosonMPO::getParityMPO(L,N,parMPO,xL,xR);
  // complex_t parit=contractor.contract(gsP,parMPO,gsP);
  // parit=parit/normGS;
  // cout<<"Found parity operator in the GS: "<<parit<<endl;

  // Now check how close we are to an eigenstate of the new Hamiltonian
  BoseHubbardHamiltonian hBH(L,N,tis,Uis,muis,VisE);
  const MPO& hamil1=hBH.getHMPO();
  complex_t variance1=contractor.contract2(gsP,hamil1,gsP);
  complex_t energy1=contractor.contract(gsP,hamil1,gsP);
  *out<<"% Delta H^2="<<variance1-energy1*energy1<<endl;

  // exit(1);
  // Now start the evolution
  // Get the exponential MPOs
  MPO expHe(L),expHo(L),expHe_2(L);
  hBH.getExponentialMPOeven(expHe,(complex_t){0,-delta});
  hBH.getExponentialMPOeven(expHe_2,(complex_t){0,-delta/2});
  hBH.getExponentialMPOodd(expHo,(complex_t){0,-delta});

  int cnt=0;
  *out<<"% Results of the time evolution "<<endl;
  *out<<"% cnt(delta="<<delta<<")\t ";
  for(int k=0;k<L;k++) *out<<"n("<<k<<")\t";
  *out<<"Ntot\t S("<<xR<<")"<<endl;

  MPS evolS(gsP);
  while(cnt<M){
    MPS aux(evolS);
    contractor.optimize(expHe_2,aux,evolS,D);
    int ks=0;
    while((ks<stepBlk-1)&&(cnt<M-1)){
      aux=evolS;
      contractor.optimize(expHo,aux,evolS,D);
      aux=evolS;
      contractor.optimize(expHe,aux,evolS,D);
      ks++;cnt++;
    }
    aux=evolS;
    contractor.optimize(expHo,aux,evolS,D);
    aux=evolS;
    contractor.optimize(expHe_2,aux,evolS,D);
    cnt++;
    double time=delta*cnt;
    // // MPS aux(evolS);
    // // contractor.optimize(expHe_2,aux,evolS,D);
    // // aux=evolS;
    // // contractor.optimize(expHo,aux,evolS,D);
    // // aux=evolS;
    // // contractor.optimize(expHe_2,aux,evolS,D);
    // After step of evolution, compute distrib of part.
    // Normalize, just in case
    evolS.gaugeCond('R',1);
    *out<<cnt<<"\t";
    double Ntot=0.;
    for(int k=0;k<L;k++){
      mpoN.setOp(k,&nrk,0);
      complex_t nrk=contractor.contract(evolS,mpoN,evolS);
      *out<<real(nrk)<<"\t";
      // restore the identity
      mpoN.setOp(k,&ident,0);
      Ntot=Ntot+real(nrk);
    }
    *out<<Ntot<<"\t";
    // Now I would also like to compute the entropy between the right 
    // part and the left one
    double entr=contractor.getEntropy(evolS,xR);
    *out<<entr<<endl;
  }

  out->close();

}
