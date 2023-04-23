
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include "BosonMPO.h"

#define SKIPPAR 1 // Not do parity calculations
#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "BoseHubbardHamiltonian.h"

/** bosehubTherm5 tries to find the thermal state of the Bose-Hubbard
    Hamiltonian \ref <BoseHubbardHamiltonian> using purification and
    imaginary time evolution, with the MPO for the exponential of H as
    second order Trotter, but splitting the auxiliary and true sites.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <N> (int) cutoff for the occupation number
    \param <t> (double) parameter \f$t\f$ of the Hamiltonian
    \param <U> (double) parameter \f$U\f$ of the Hamiltonian
    \param <mu> (double) parameter \f$\mu\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <beta> (double) inverse temperature, \f$\beta\f$
    \param <delta> (double) width of time step (actually, it is 
    the delta in the exponent, and the time step turns out to be 2*delta)
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

  int M=beta/(2*delta);


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


  //  double offset=0.; // done by Contractor, if needed
  BoseHubbardHamiltonian hBH(L,N,t,U,mu);
  int DopH=4; // bod dim of the BH Hamiltonian
  //SchwingerHamiltonian hSchM(L,mu,x,L0M,alpha);
  const MPO& hamil=hBH.getHMPO();
  cout<<"Initialized basic Hamiltonian"<<endl;

  complex_t deltaC=(complex_t){-delta,0.}; // delta  as complex
  complex_t deltaC_2=(complex_t){-delta*.5,0.}; // delta/2

  // Get the exponential MPOs
  MPO expHe(2*L),expHo(2*L),expHe_2(2*L);
  hBH.getExponentialMPOeven(expHe,deltaC);
  hBH.getExponentialMPOeven(expHe_2,deltaC_2);
  hBH.getExponentialMPOodd(expHo,deltaC);

  int Dop=d0*d0; // bond dimension of Bose-Hubbard expHe

  // Operator on the even sites
  mwArray identPhys=identityMatrix(d0);
  mwArray identD=identityMatrix(Dop);
  identPhys.reshape(shrt::Indices(d0*d0,1)); 
  identD.reshape(shrt::Indices(1,Dop*Dop)); 
  mwArray _idOp=identPhys*identD;
  _idOp.reshape(shrt::Indices(d0,d0,Dop,Dop));
  _idOp.permute(shrt::Indices(1,3,2,4));
  Operator idOp(_idOp);
  // Cut down the bond dim for the H operator!
  _idOp.resize(shrt::Indices(d0,DopH,d0,DopH));
  Operator idOpH(_idOp);
  int DopN=2;
  // Cut down the bond dim for the N operator!
  _idOp.resize(shrt::Indices(d0,DopN,d0,DopN));
  Operator idOpN(_idOp);
  // the last one is different, because Dop should be one now
  identPhys.reshape(shrt::Indices(d0,1,d0,1));
  Operator idOpL(identPhys);

  // The MPO for H:
  MPO unfoldedH(2*L);
  for(int k=0;k<L;k++){
    unfoldedH.setOp(2*k,&hamil.getOp(k));
    if(k<L-1)
      unfoldedH.setOp(2*k+1,&idOpH);
    else
      unfoldedH.setOp(2*k+1,&idOpL);
  }  

  // Now combine the operators, to write
  // the expH X conj(expH) on the doubled system
  MPO unexpHe(2*L),unexpHo(2*L),unexpHe_2(2*L);
  // first, without joining the MPOs, I construct the long versions of these ones
  // and the complex conjugates for the other side
  for(int k=0;k<L;k++){
    unexpHe.setOp(2*k,&expHe.getOp(k));
    unexpHo.setOp(2*k,&expHo.getOp(k));
    unexpHe_2.setOp(2*k,&expHe_2.getOp(k));
    if(k%2==0){
      if(k!=L-1) unexpHe.setOp(2*k+1,&idOp);
      else unexpHe.setOp(2*k+1,&idOpL);
      if(k!=L-1) unexpHe_2.setOp(2*k+1,&idOp);
      else unexpHe_2.setOp(2*k+1,&idOpL);
      unexpHo.setOp(2*k+1,&idOpL);
    }
    else{
      unexpHe.setOp(2*k+1,&idOpL);
      unexpHe_2.setOp(2*k+1,&idOpL);
      if(k!=L-1) unexpHo.setOp(2*k+1,&idOp);
      else unexpHo.setOp(2*k+1,&idOpL);
    }
  }
  /*  cout<<"Constructed exponentials on system, filled with ancillas "<<endl;
  cout<<"expHe->"<<unexpHe<<endl;
  cout<<"expHo->"<<unexpHo<<endl;
  cout<<"expHe2->"<<unexpHe_2<<endl;*/

  // And the conjugates acting on acillas (for these, I have to copy
  // all the operators to conjugate them!)
  MPO unexpHeA(2*L),unexpHoA(2*L),unexpHe_2A(2*L);
  for(int k=0;k<L;k++){
    unexpHeA.setConjugatedOp(2*k+1,&expHe.getOp(k));
    unexpHoA.setConjugatedOp(2*k+1,&expHo.getOp(k));
    unexpHe_2A.setConjugatedOp(2*k+1,&expHe_2.getOp(k));
    if(k%2!=0){
      unexpHeA.setOp(2*k,&idOp);
      unexpHe_2A.setOp(2*k,&idOp);
      unexpHoA.setOp(2*k,&idOpL);
    }
    else{
      unexpHeA.setOp(2*k,&idOpL);
      unexpHe_2A.setOp(2*k,&idOpL);
      if(k!=0) unexpHoA.setOp(2*k,&idOp);
      else unexpHoA.setOp(2*k,&idOpL);
    }
  }
  /*  cout<<"Constructed exponentials on ancillas, filled with id on stm "<<endl;
  cout<<"expHeA->"<<unexpHeA<<endl;
  cout<<"expHoA->"<<unexpHoA<<endl;
  cout<<"expHe2A->"<<unexpHe_2A<<endl; */

  // Right now, the joinedOperators contain copies!  this is very
  // inefficient here, because I will be making many unnecessary
  // copies of the arrays.
  // Place for the MPOs
  MPO unexpHeT(2*L),unexpHoT(2*L),unexpHe_2T(2*L);
  const MPO* aux[2];
  aux[0]=&unexpHeA;aux[1]=&unexpHe; 
  MPO::join(2,aux,unexpHeT);
  aux[0]=&unexpHe_2A;aux[1]=&unexpHe_2; 
  MPO::join(2,aux,unexpHe_2T);
  aux[0]=&unexpHoA;aux[1]=&unexpHo; 
  MPO::join(2,aux,unexpHoT);


  // Check the ground state?
  double lambda=0.;
  MPS gS(2*L,D,d0);gS.setRandomState();
  contractor.findGroundState(unfoldedH,D,&lambda,gS);
  cout<<"Ground state found with eigenvalue "<<lambda<<endl;
  *out<<"% EnergyGS="<<lambda<<"\t";

  exit(1);

  // MPO for number op
  MPO Nmpo_(L);
  MPO Nmpo(2*L);
  BosonMPO::getNumberMPO(L,N,Nmpo_);
  // MPO for different parity strings
  MPO parMPO_(L); // basic one
  MPO parMPO(2*L); // with ancillas, all holding identity
  for(int k=0;k<L;k++){
    parMPO.setOp(2*k+1,&idOpL);
    Nmpo.setOp(2*k,&Nmpo_.getOp(k));
    if(k==L-1) Nmpo.setOp(2*k+1,&idOpL);
    else Nmpo.setOp(2*k+1,&idOpN);
  }

  complex_t NtotGS=contractor.contract(gS,Nmpo,gS);
  *out<<NtotGS<<endl;

  int posCent=L/2; // central site

  // Now try an optimization
  int cnt=0;
  double betai=0;

  int stepBlk=5; //blocks of steps

  complex_t energy0=contractor.contract(thS,unfoldedH,thS);
  *out<<delta<<"\t"<<cnt<<"\t"<<betai<<"\t"<<1/betai<<"\t"
      <<energy0<<endl;
  complex_t Ntot0=contractor.contract(thS,Nmpo,thS);
  cout<<"Before any step, beta="<<betai<<", T="<<1/betai
      <<", energy="<<energy0<<"\t"<<Ntot0<<endl;
  while(cnt<M){
    MPS aux(thS);
    contractor.optimize(unexpHe_2T,aux,thS,D);
    cout<<"Optimized unexpHe_2 step "<<cnt<<endl;
    int k=0;
    while((k<stepBlk-1)&&(cnt<M-1)){
      aux=thS;
      contractor.optimize(unexpHoT,aux,thS,D);
      cout<<"Optimized unexpHo step "<<cnt<<endl;
      aux=thS;
      contractor.optimize(unexpHeT,aux,thS,D);
      cout<<"Optimized unexpHe step "<<cnt<<endl;
      betai+=2*delta;
      k++;cnt++;
    }
    aux=thS;
    contractor.optimize(unexpHoT,aux,thS,D);
    cout<<"Optimized unexpHo step "<<cnt<<endl;
    aux=thS;
    contractor.optimize(unexpHe_2T,aux,thS,D);
    cout<<"Optimized unexpHe_2 step "<<cnt<<endl;
    betai+=2*delta;
    cnt++;
    // Now I need to normalize
    thS.gaugeCond('R',1);
    complex_t energy=contractor.contract(thS,unfoldedH,thS);
    complex_t Ntot=contractor.contract(thS,Nmpo,thS);

    *out<<delta<<"\t"<<cnt<<"\t"<<betai<<"\t"<<1/betai<<"\t"
	<<energy<<"\t"<<Ntot<<"\t";
    cout<<"At the end of step "<<cnt<<", beta="<<betai<<", T="<<1/betai
	<<", energy="<<energy<<", N="<<Ntot<<endl;

#ifndef SKIPPAR
    // And all the parity values
    double nbar=real(Ntot/L); // average n
    for(int k=1;k<L/2;k++){
      BosonMPO::getParityMPO(L,N,parMPO_,posCent-k,posCent+k,nbar);
      // Now I need to fill in the identities on ancillas!
      for(int ka=0;ka<L;ka++)
	parMPO.setOp(2*ka,&parMPO_.getOp(ka));
      //cout<<parMPO<<endl;
      complex_t parit=contractor.contract(thS,parMPO,thS);
      *out<<posCent-k<<"\t"<<posCent+k<<"\t"<<parit<<"\t";
      cout<<"Parity("<<posCent-k<<","<<posCent+k<<")="<<parit<<endl;
    }
#endif
    *out<<endl;

  }
  *out<<delta<<"\t Inf\t Inf\t"<<0<<"\t"
      <<lambda<<"\t"<<NtotGS<<endl;
  cout<<"And in the limit beta=Inf, T=0,"
      <<" GS energy="<<lambda<<endl;

}
