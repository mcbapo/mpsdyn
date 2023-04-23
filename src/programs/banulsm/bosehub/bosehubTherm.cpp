
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "BoseHubbardHamiltonian.h"

/** bosehubTherm tries to find the thermal state of the Bose-Hubbard
    Hamiltonian \ref <BoseHubbardHamiltonian> using purification and
    imaginary time evolution.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <N> (int) cutoff for the occupation number
    \param <t> (double) parameter \f$t\f$ of the Hamiltonian
    \param <U> (double) parameter \f$U\f$ of the Hamiltonian
    \param <mu> (double) parameter \f$\mu\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <beta> (double) inverse temperature, \f$\beta\f$
    \param <delta> (double) width of time step
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)

*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  int N=atoi(argv[++cntr]);
  double t=atof(argv[++cntr]);
  double U=atof(argv[++cntr]);
  double mu=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  double beta=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  cout<<"Initialized arguments: L="<<L
      <<", t="<<t
      <<", U="<<U
      <<", mu="<<mu
      <<", Tmin="<<1/beta
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% L="<<L<<", t="<<t<<", U="<<U<<", mu="<<mu
	<<", Nmax="<<N<<", D="<<D<<endl;
    *out<<"% delta\t nrsteps\t beta\t T\t Energy"
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

  int M=beta/(delta);


  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d0=N+1;
  int d=d0*d0;
  MPS thS(L,D,d);
  thS.setProductState(p_maxent); // the initial state, maximally entangled
  cout<<"Initialized identity state"<<endl;

  //int Dop=5; // bond dimension of Bose-Hubbard H

  //  double offset=0.; // done by Contractor, if needed
  BoseHubbardHamiltonian hBH(L,N,t,U,mu);
  //SchwingerHamiltonian hSchM(L,mu,x,L0M,alpha);
  const MPO& hamil=hBH.getHMPO();
  cout<<"Initialized basic Hamiltonian"<<endl;

  // now an auxiliary MPO, just identities
  mwArray ident=identityMatrix(d0);
  ident.reshape(shrt::Indices(d0,1,d0,1)); // just dx1xdx1

  MPO foldedH(L);
  for(int k=0;k<L;k++){
    // Remember left one is reflected
    foldedH.setOp(k,new FoldedOperator(ident,hamil.getOpData(k)),true);
    cout<<"Set FoldedOperator at position "<<k<<endl;
  }

  // Check the ground state?
  //  double lambda=0.;
  // MPS gS(L,D,d);gS.setRandomState();
  //gS.gaugeCond('R',1);
  //contractor.findGroundState(foldedH,D,&lambda,gS);
  //cout<<"Ground state found with eigenvalue "<<lambda<<endl;

  // Now try an optimization
  int cnt=0;
  double betai=0;
  // For convergence of the method, it is often worth to use offset<~E_GS
  double offset=-16.610; 
  while(cnt<M){
    MPS aux(thS);
    contractor.approximateExponential(foldedH,aux,thS,delta,D,imag,offset);
    // Now I need to normalize
    thS.gaugeCond('R',1);
    betai+=delta;
    cnt++;
    complex_t energy=contractor.contract(thS,foldedH,thS);
    *out<<delta<<"\t"<<cnt<<"\t"<<betai<<"\t"<<1/betai<<"\t"
	<<energy<<endl;
    cout<<"At the end of step "<<cnt<<", beta="<<betai<<", T="<<1/betai
	<<", energy="<<energy<<endl;
  }
  

}
