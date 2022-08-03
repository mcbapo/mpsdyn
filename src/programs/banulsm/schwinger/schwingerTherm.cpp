
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "SchwingerHamiltonian.h"

/** schwingerTherm tries to find the thermal state of the Schwinger
    Hamiltonian \ref <SchwingerHamiltonian> using purification and
    imaginary time evolution.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <mg> (double) parameter \f$m/g\f$ of the Hamiltonian
    \param <x> (double) parameter \f$x\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <T> (double) temperature
    \param <delta> (double) width of time step (actually, this is the 
                   coeff. in the exponential)
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)

*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double mg=atof(argv[++cntr]);
  double x=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  double beta=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  double mu=2*mg*sqrt(x);

  double alpha=.5;
  double L0P=0.;
  //double L0M=-1.;


  cout<<"Initialized arguments: L="<<L
      <<", mu="<<mu
      <<", x="<<x
      <<", Tmin="<<1/beta
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"# N\t x\t mg\t D\t delta\t nrsteps\t beta\t T\t Energy"
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
  int d0=2;
  int d=d0*d0;
  MPS thS(L,D,d);
  thS.setProductState(p_maxent); // the intial state, maximally entangled
  cout<<"Initialized identity state"<<endl;

  int Dop=5; // bond dimension of Schwinger H

  double offset=0.;
  SchwingerHamiltonian hSchP(L,mu,x,L0P,alpha,offset);
  //SchwingerHamiltonian hSchM(L,mu,x,L0M,alpha);
  const MPO& hamilP=hSchP.getHMPO();
  cout<<"Initialized basic Hamiltonian"<<endl;

  // now an auxiliary MPO, just identities
  mwArray ident=identityMatrix(d0);
  ident.reshape(shrt::Indices(d0,1,d0,1)); // just dx1xdx1

  MPO foldedH(L);
  for(int k=0;k<L;k++){
    // Remember left one is reflected
    foldedH.setOp(k,new FoldedOperator(ident,hamilP.getOpData(k)),true);
    cout<<"Set FoldedOperator at position "<<k<<endl;
  }

  // Check the ground state?
  double lambda=0.;
  MPS gS(L,D,d);gS.setRandomState();
  contractor.findGroundState(foldedH,D,&lambda,gS);
  cout<<"Ground state (sector +1/2) found with eigenvalue "<<lambda<<endl;

  // Now try an optimization
  int cnt=0;
  double betai=0;
  while(cnt<M){
    MPS aux(thS);
    contractor.approximateExponential(foldedH,aux,thS,delta,D,imag,-580);//,-1175,2500);
    // Now I need to normalize
    thS.gaugeCond('R',1);
    betai+=delta;
    cnt++;
    complex_t energy=contractor.contract(thS,foldedH,thS);
    *out<<L<<"\t"<<x<<"\t"<<mg<<"\t"<<D<<"\t"
	<<delta<<"\t"<<cnt<<"\t"<<betai<<"\t"<<1/betai<<"\t"
	<<energy<<endl;
    cout<<"At the end of step "<<cnt<<", beta="<<betai<<", T="<<1/betai
	<<", energy="<<energy<<endl;
  }
  

}
