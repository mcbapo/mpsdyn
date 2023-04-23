
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include "BosonMPO.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "BoseHubbardHamiltonian.h"

/** bosehubTrapTherm tries to find the thermal state of the Bose-Hubbard
    Hamiltonian \ref <BoseHubbardHamiltonian> inside a trap, using 
    purification and
    imaginary time evolution, with the MPO for the exponential of H as
    second order Trotter, but splitting the auxiliary and true sites.
    This version will create, in the output directory, a file every 5 steps,
    with the parity string and correlations, and the number of particles 
    per site.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <N> (int) cutoff for the occupation number
    \param <t> (double) parameter \f$t\f$ of the Hamiltonian
    \param <U> (double) parameter \f$U\f$ of the Hamiltonian
    \param <mu> (double) parameter \f$\mu\f$ of the Hamiltonian
    \param <V> (double) parameter \f$V_0\f$ of the Hamiltonian 
            (asume center in the middle, (L+1)/2 )
    \param <D> (int) maximum bond dimension
    \param <beta> (double) inverse temperature, \f$\beta\f$
    \param <delta> (double) width of time step (actually, it is 
    the delta in the exponent, and the time step turns out to be 2*delta)
    \param <stepBlk> (int) save results every stepBlk steps (means
            every delta*stepBlk*2(4) in beta)
    \param <outdname> (char*) name of the output directory

*/

void computeCorrelations(int Nmax,const MPS& state,const char* filename,
			 const MPS& idState,const MPO& unfoldedH,double beta);
void getMPON(int L,int N,int i,MPO& NmpoI);
void getMPOs(int L,int N,int i,int j,MPO& parityString,MPO& SiSj);
const MPO& getMPONtot(int L,int Nmax);


int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  int N=atoi(argv[++cntr]);
  double t=atof(argv[++cntr]);
  double U=atof(argv[++cntr]);
  double mu=atof(argv[++cntr]);
  double V=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  double beta=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int stepBlk=atoi(argv[++cntr]);
  const char* outdname=argv[++cntr];

  cout<<"Initialized arguments: L="<<L
      <<", t="<<t
      <<", U="<<U
      <<", mu="<<mu
      <<", V="<<V
      <<", Tmin="<<1/beta
      <<", outdir="<<outdname
      <<", D="<<D<<endl;

  double x0=(L+1)*.5;


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
  MPS idState(thS);

  /* With this implementation, BoseHubbardHamiltonian needs all
     parameters in all sites of the lattice to create a
     non-homogeneous system.*/
  double* tis=new double[L];
  double* Uis=new double[L];
  double* muis=new double[L];
  double* Vis=new double[L];

  // fill in values as needed
  for(int k=0;k<L;k++){
    tis[k]=t;Uis[k]=U;muis[k]=mu;
    Vis[k]=V*(k-x0)*(k-x0);
  }

  BoseHubbardHamiltonian hBH(L,N,tis,Uis,muis,Vis);

  // delete the auxiliary double vectors
  delete []tis;
  delete []Uis;
  delete []muis;
  delete []Vis;

  int DopH=hBH.getBondDimension(); // bond dim of the BH Hamiltonian
  const MPO& hamil=hBH.getHMPO();
  cout<<"Initialized basic Hamiltonian"<<endl;

  complex_t deltaC=(complex_t){-delta,0.}; // delta  as complex
  complex_t deltaC_2=(complex_t){-delta*.5,0.}; // delta/2

  // Get the exponential MPOs
  // WARNING!! Length should be L, half is empty, now!!
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
  /*  cout<<"Constructed exponentials on ancillas, filled with id on stm "
      <<endl;
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
  //  double lambda=0.;
  //MPS gS(2*L,D,d0);gS.setRandomState();
  //contractor.findGroundState(unfoldedH,D,&lambda,gS);
  //cout<<"Ground state found with eigenvalue "<<lambda<<endl;

  int posCent=L/2; // central site


  // Now try an optimization
  int cnt=0;
  double betai=0;

  char outfile[200];
  sprintf(outfile,"%s/beta%g.dat",outdname,betai*2);
  cout<<"Writing to file "<<outfile<<endl;
  computeCorrelations(N,thS,outfile,idState,unfoldedH,betai);

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

    // Create a name for the file where I want to write the results
    // for this time step
    char outfile[200];
    sprintf(outfile,"%s/beta%g.dat",outdname,betai);
    cout<<"Writing to file "<<outfile<<endl;
    computeCorrelations(N,thS,outfile,idState,unfoldedH,betai);
  }
}

const MPO& getMPONtot(int L,int Nmax){
  int d=Nmax+1;
  /** Keep the N perator for the next calls */
  // MPO for number op
  static MPO Nmpo(L),Nmpo_(L/2);
  static bool initialized=false;
  static Operator* idOpN,*idOpL;
  if(!initialized){
    //cout<<"Constructing Nmpo"<<endl;
    int Dop=2;
    mwArray identPhys=identityMatrix(d);
    mwArray identD=identityMatrix(Dop);
    identPhys.reshape(shrt::Indices(d*d,1)); 
    identD.reshape(shrt::Indices(1,Dop*Dop)); 
    mwArray _idOp=identPhys*identD;
    _idOp.reshape(shrt::Indices(d,d,Dop,Dop));
    _idOp.permute(shrt::Indices(1,3,2,4));
    idOpN=new Operator(_idOp);
    // the last one is different, because Dop should be one now
    identPhys.reshape(shrt::Indices(d,1,d,1));
    idOpL=new Operator(identPhys);
    BosonMPO::getNumberMPO(L/2,Nmax,Nmpo_);
    for(int k=0;k<L/2;k++){
      Nmpo.setOp(2*k,&Nmpo_.getOp(k));
      if(k==L/2-1) Nmpo.setOp(2*k+1,idOpL);
      else Nmpo.setOp(2*k+1,idOpN);
    }
    initialized=true;
  }
  //cout<<"Constructed Nmpo"<<Nmpo<<endl;
  return Nmpo;
}

using namespace shrt;

void getMPON(int L,int N,int i,MPO& NmpoI){
  int d=N+1;
  static bool initialized=false;
  static Operator* idOpPhys,* Nop;
  if(!initialized){
    mwArray identPhys=identityMatrix(d);
    identPhys.reshape(shrt::Indices(d,1,d,1));
    idOpPhys=new Operator(identPhys);
    mwArray locN(Indices(d,d));locN.fillWithZero();
    for(int k=0;k<d;k++){
      locN.setElement(k,0,Indices(k,k)); // exp(i pi N)
    }
    locN.reshape(Indices(d,1,d,1));
    Nop=new Operator(locN);
    initialized=true;
  }
  NmpoI.initLength(2*L);
  for(int k=0;k<L;k++){
    if(k!=i)
      NmpoI.setOp(2*k,idOpPhys);
    else
      NmpoI.setOp(2*k,Nop);
    NmpoI.setOp(2*k+1,idOpPhys);
  }
}

void getMPOs(int L,int N,int i,int j,MPO& parityString,MPO& SiSj){
  int d=N+1;
  static bool initialized=false;
  static Operator* idOpPhys,* Sop;
  if(!initialized){
    mwArray identPhys=identityMatrix(d);
    identPhys.reshape(shrt::Indices(d,1,d,1));
    idOpPhys=new Operator(identPhys);
    mwArray expN(Indices(d,d));expN.fillWithZero();
    for(int k=0;k<d;k++){
      expN.setElement(exp((k)*M_PIl*I_c),Indices(k,k)); // exp(i pi N)
    }
    expN.reshape(Indices(d,1,d,1));
    Sop=new Operator(expN);
    initialized=true;
  }
  parityString.initLength(2*L);
  SiSj.initLength(2*L);
  for(int k=0;k<L;k++){
    if(k<i||k>j){
      parityString.setOp(2*k,idOpPhys);
      SiSj.setOp(2*k,idOpPhys);
    }
    else{
      parityString.setOp(2*k,Sop);
      if((k==i)||(k==j))
	SiSj.setOp(2*k,Sop);
      else
	SiSj.setOp(2*k,idOpPhys);
    }
    parityString.setOp(2*k+1,idOpPhys);
    SiSj.setOp(2*k+1,idOpPhys);
  }

}


void computeCorrelations(int Nmax,const MPS& state,const char* filename,
			 const MPS& idState,const MPO& unfoldedH,
			 double beta){
  // Open file (append?)
  //  ofstream* out=new ofstream(filename,ios::app);
  ofstream* out=new ofstream(filename);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<filename<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);
  int L=state.getLength();
  int d=Nmax+1;

  /** Two ways of computing: Psi Op Psi will actually compute 
      trace(exp(-2 delta*(2k) H) Op), because Psi already contains 
      exp(-2*k*delta H). So, I can use Psi to compute the expectation
      value of beta= 2*k*delta and 4*k*delta in the same step */

  // Normalize the state
  Contractor& contractor=Contractor::theContractor();
  double norm2=real(contractor.contract(state,state)); // trace of rho2
  double norm1=real(contractor.contract(state,idState)); //trace of rho
  const MPO& Ntot=getMPONtot(L,Nmax);
  complex_t ntot2=contractor.contract(state,Ntot,state);
  complex_t ntot=contractor.contract(state,Ntot,idState);
  complex_t energy2=contractor.contract(state,unfoldedH,state);
  complex_t energy=contractor.contract(state,unfoldedH,idState);
  *out<<"% beta= "<<beta<<", Energy(beta)="<<real(energy)/norm1
      <<", Ntot(beta)="<<real(ntot)/norm1<<"\n"
      <<"% 2*beta="<<2*beta<<" Energy(2*beta)="<<real(energy2)/norm2
      <<", Ntot(2*beta)="<<real(ntot2)/norm2<<endl;

  mwArray Pij(Indices(L/2,L/2));Pij.fillWithZero();
  mwArray Sij(Indices(L/2,L/2));Sij.fillWithZero();
  mwArray Ni(Indices(L/2));Ni.fillWithZero();
  mwArray Pij2(Indices(L/2,L/2));Pij2.fillWithZero();
  mwArray Sij2(Indices(L/2,L/2));Sij2.fillWithZero();
  mwArray Ni2(Indices(L/2));Ni2.fillWithZero();

  // Compute and store all correlations
  MPO parityString(L),SiSj(L),NmpoI(L);
  for(int i=0;i<L/2;i++){ 
    // Compute expectation value of number
    getMPON(L/2,Nmax,i,NmpoI);
    complex_t nval2=contractor.contract(state,NmpoI,state);
    complex_t nval=contractor.contract(state,NmpoI,idState);
    Ni2.setElement(nval2/norm2,Indices(i));
    Ni.setElement(nval/norm1,Indices(i));
    for(int j=i;j<L/2;j++){
      getMPOs(L/2,Nmax,i,j,parityString,SiSj);      
      // compute and write
      complex_t corrPar2=contractor.contract(state,parityString,state);
      complex_t corrPar=contractor.contract(state,parityString,idState);
      complex_t corrSS2=contractor.contract(state,SiSj,state);
      complex_t corrSS=contractor.contract(state,SiSj,idState);
      Pij2.setElement(corrPar2/norm2,Indices(i,j));
      Sij2.setElement(corrSS2/norm2,Indices(i,j));
      Pij.setElement(corrPar/norm1,Indices(i,j));
      Sij.setElement(corrSS/norm1,Indices(i,j));
    }
  }
  // Now write to the file
  *out<<"% ******************* beta *********************** "<<endl;
  *out<<"\n\n % For beta="<<beta<<endl;
  // N operator for a site
  *out<<"% i\t <n_i>"<<endl;
  for(int i=0;i<L/2;i++){
    *out<<i<<"\t"<<real(Ni.getElement(Indices(i)))<<"\t"<<imag(Ni.getElement(Indices(i)))<<endl;
  }
  *out<<"\n\n"<<endl;
  *out<<"% beta="<<beta<<", Two body correlations, S_i=exp(i pi n_i)"<<endl;
  *out<<"% i\t j \t <S_i S_j>"<<endl;
  for(int i=0;i<L/2;i++){
    for(int j=i;j<L/2;j++){
      complex_t sel=Sij.getElement(Indices(i,j));
      *out<<i<<"\t"<<j<<"\t"<<real(sel)<<"\t"<<imag(sel)<<endl;
    }
  }
  // And the parity string operator
  *out<<"\n\n"<<endl;
  *out<<"% beta="<<beta<<",Parity string operator Pij:=<S_i S_[i+1]...S_j>" <<endl;
  *out<<"% i\t j \t <P_ij>"<<endl;
  for(int i=0;i<L/2;i++){
    for(int j=i;j<L/2;j++){
      complex_t pel=Pij.getElement(Indices(i,j));
      *out<<i<<"\t"<<j<<"\t"<<real(pel)<<"\t"<<imag(pel)<<endl;
    }
  }

  *out<<"% ******************* 2*beta *********************** "<<endl;
  *out<<"\n\n % For 2*beta="<<2*beta<<endl;
  // N operator for a site
  *out<<"% i\t <n_i>"<<endl;
  for(int i=0;i<L/2;i++){
    *out<<i<<"\t"<<real(Ni2.getElement(Indices(i)))<<"\t"<<imag(Ni2.getElement(Indices(i)))<<endl;
  }
  *out<<"\n\n"<<endl;
  *out<<"% 2*beta="<<2*beta<<", Two body correlations, S_i=exp(i pi n_i)"<<endl;
  *out<<"% i\t j \t <S_i S_j>"<<endl;
  for(int i=0;i<L/2;i++){
    for(int j=i;j<L/2;j++){
      complex_t sel=Sij2.getElement(Indices(i,j));
      *out<<i<<"\t"<<j<<"\t"<<real(sel)<<"\t"<<imag(sel)<<endl;
    }
  }
  // And the parity string operator
  *out<<"\n\n"<<endl;
  *out<<"% 2*beta="<<2*beta<<",Parity string operator Pij:=<S_i S_[i+1]...S_j>" <<endl;
  *out<<"% i\t j \t <P_ij>"<<endl;
  for(int i=0;i<L/2;i++){
    for(int j=i;j<L/2;j++){
      complex_t pel=Pij2.getElement(Indices(i,j));
      *out<<i<<"\t"<<j<<"\t"<<real(pel)<<"\t"<<imag(pel)<<endl;
    }
  }

  // Close file
  out->close();
}
