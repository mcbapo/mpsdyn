
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "BoseHubbardHamiltonian.h"

/** bosehubTherm2 tries to find the thermal state of the Bose-Hubbard
    Hamiltonian \ref <BoseHubbardHamiltonian> using purification and
    imaginary time evolution, but splitting the auxiliary and true sites.

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

    WARNING!!!! Different to bosehubTherm4(5), the beta computed is
    actually the argument given and not twice as much!  This is so
    because here I only apply the operaton on the system, and
    identities on the ancillas.

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
  MPS thS(2*L,D,d0);
  //cout<<"Initialized MPS "<<thS<<endl;
  // Tensor for the even sites (left part of the max. entgld)
  mwArray Aeven(shrt::Indices(d0,D,D));
  // Tensor for the odd sites (right part of the max. entgld)
  mwArray Aodd(shrt::Indices(d0,D,D));
  double val=1./sqrt(sqrt(d0));
  for(int k=0;k<d0;k++){
    Aeven.setElement(val,0.,shrt::Indices(k,0,k));
    Aodd.setElement(val,0.,shrt::Indices(k,k,0));
  }
  // Now set the sites in the MPS
  for(int k=0;k<L;k++){
    if(k==0){
      mwArray firstSite(Aeven);
      firstSite.resize(shrt::Indices(d0,1,D));
      thS.setA(2*k,firstSite);
    }
    else
      thS.setA(2*k,Aeven);
    if(k==L-1){
      mwArray lastSite(Aodd);
      lastSite.resize(shrt::Indices(d0,D,1));
      thS.setA(2*k+1,lastSite);
    }
    else
      thS.setA(2*k+1,Aodd);
  }
  cout<<"Initialized identity state. Check norm="
      <<contractor.contract(thS,thS)<<endl;

  int Dop=4; // bond dimension of Bose-Hubbard H

  //  double offset=0.; // done by Contractor, if needed
  BoseHubbardHamiltonian hBH(L,N,t,U,mu);
  //SchwingerHamiltonian hSchM(L,mu,x,L0M,alpha);
  const MPO& hamil=hBH.getHMPO();
  cout<<"Initialized basic Hamiltonian"<<endl;

  // Now combine this ops with identities, to write
  // the H on the doubled system
  MPO unfoldedH(2*L);

  // Operator on the even sites
  mwArray identPhys=identityMatrix(d0);
  mwArray identD=identityMatrix(Dop);
  identPhys.reshape(shrt::Indices(d0*d0,1)); 
  identD.reshape(shrt::Indices(1,Dop*Dop)); 
  mwArray _idOp=identPhys*identD;
  _idOp.reshape(shrt::Indices(d0,d0,Dop,Dop));
  _idOp.permute(shrt::Indices(1,3,2,4));
  Operator idOp(_idOp);
  // the last one is different, because Dop should be one now
  identPhys.reshape(shrt::Indices(d0,1,d0,1));
  Operator idOpL(identPhys);

  for(int k=0;k<L;k++){
    unfoldedH.setOp(2*k,&hamil.getOp(k));
    if(k<L-1)
      unfoldedH.setOp(2*k+1,&idOp);
    else
      unfoldedH.setOp(2*k+1,&idOpL);
  }  

  // Check the ground state?
  double lambda=0.;
  MPS gS(2*L,D,d0);gS.setRandomState();
  contractor.findGroundState(unfoldedH,D,&lambda,gS);
  cout<<"Ground state found with eigenvalue "<<lambda<<endl;

  // Now try an optimization
  int cnt=0;
  double betai=0;
  // For convergence of the method, it is often worth to use offset<~E_GS
  double offset=-16.610; 
  while(cnt<M){
    MPS aux(thS);
    contractor.approximateExponential(unfoldedH,aux,thS,delta,D,imag,offset);
    // Now I need to normalize
    thS.gaugeCond('R',1);
    betai+=delta;
    cnt++;
    complex_t energy=contractor.contract(thS,unfoldedH,thS);
    *out<<delta<<"\t"<<cnt<<"\t"<<betai<<"\t"<<1/betai<<"\t"
	<<energy<<endl;
    cout<<"At the end of step "<<cnt<<", beta="<<betai<<", T="<<1/betai
	<<", energy="<<energy<<endl;
  }
  

}
