#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#include "BoseHubbardHamiltonian.h"
#include "BosonMPO.h"

using namespace shrt;

/** bosehubWithInt finds the ground state of the Bose-Hubbard
    Hamiltonian     \ref <BoseHubbardHamiltonian>
    with constant (homogeneous) parameters, including a NN interaction
    term. At the end, it outputs the energy of the ground state and
    the particle densities.
    // TODO: Then, add one particle in the center, evolve in real time
    // and plot entropy growth.

    Receives arguments:
    \param <L> (int) total length of the chain
    \param <N> (int) cutoff for the occupation number
    \param <t> (double) parameter \f$t\f$ of the Hamiltonian
    \param <U> (double) parameter \f$U\f$ of the Hamiltonian
    \param <mu> (double) parameter \f$mu\f$ of the Hamiltonian
    \param <V> (double) parameter \f$V\f$ of the Hamiltonian
    \param <Vint> (double) parameter \f$Vint\f$ of the Hamiltonian
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
  double t=atof(argv[++cntr]);
  double U=atof(argv[++cntr]);
  double mu=atof(argv[++cntr]);
  double V=atof(argv[++cntr]);
  double Vint=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int stepBlk=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% L="<<L<<", N="<<N<<endl;
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

  BoseHubbardHamiltonian hBH0(L,N,t,U,mu,V,Vint);
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
  *out<<"% t="<<t<<endl;
  *out<<"% U="<<U<<endl;
  *out<<"% mu="<<mu<<endl;
  *out<<"% V="<<V<<endl;
  *out<<"% Vint="<<Vint<<endl;
  *out<<"%k\t <n(k)> \t (Delta nk)^2"<<endl;
  // and at the same time compute total number of particles
  double Ntot=0.;
  for(int k=0;k<L;k++){
    mpoN.setOp(k,&nrk,0);
    complex_t nrk_=contractor.contract(gsP,mpoN,gsP);
    complex_t nrk2=contractor.contract2(gsP,mpoN,gsP);
    nrk_=nrk_/normGS;
    nrk2=nrk2/normGS;
    *out<<k  //<<"\t"<<t<<"\t"<<U<<"\t"<<mu<<"\t"<<V<<"\t"<<Vint
	<<"\t"<<nrk_<<"\t"<<nrk2-nrk_*nrk_<<endl;
    // restore the identity
    mpoN.setOp(k,&ident,0);
    Ntot=Ntot+real(nrk_);
  }

  cout<<"Total number of particles N="<<Ntot<<endl;
  *out<<"%\n Ntot="<<Ntot<<endl;

  // Now apply creation operator on the central site
  mwArray aDagger(Indices(d,d));
  for(int n=0;n<N;n++){
    aDagger.setElement(ONE_c*sqrt(n+1),Indices(n+1,n));
  }
  aDagger.reshape(Indices(d,1,d,1));
  int center=L%2==0?L/2:(L-1)/2;
  cout<<"Cut located in "<<center<<endl;
  *out<<"%Before perturbation, S("<<center<<")="<<contractor.getEntropy(gsP,center)<<endl;
  Operator* aDagop=new Operator(aDagger);
  mpoN.setOp(center,aDagop,false);
  MPS init(gsP);
  {
    MPS tmp(init);
    contractor.optimize(mpoN,tmp,init,D);
    mpoN.setOp(center,&ident,false);
  }
  delete aDagop;

  complex_t normInit=contractor.contract(init,init);
  complex_t energyInit=contractor.contract(init,hamil,init);

  *out<<"% Perturbation applied, now S("<<center<<")="<<contractor.getEntropy(init,center)<<endl;
  *out<<"%k\t <n(k)> \t (Delta nk)^2"<<endl;
  Ntot=0.;
  for(int k=0;k<L;k++){
    mpoN.setOp(k,&nrk,0);
    complex_t nrk_=contractor.contract(init,mpoN,init);
    complex_t nrk2=contractor.contract2(init,mpoN,init);
    nrk_=nrk_/normInit;
    nrk2=nrk2/normInit;
    *out<<k  //<<"\t"<<t<<"\t"<<U<<"\t"<<mu<<"\t"<<V<<"\t"<<Vint
	<<"\t"<<nrk_<<"\t"<<nrk2-nrk_*nrk_<<endl;
    // restore the identity
    mpoN.setOp(k,&ident,0);
    Ntot=Ntot+real(nrk_);
  }

  cout<<"Total number of particles after creating one, N="<<Ntot<<endl;
  cout<<"New energy "<<energyInit/normInit<<endl;


  // Now start the evolution
  // Get the exponential MPOs
  MPO expHe(L),expHo(L),expHe_2(L);
  hBH0.getExponentialMPOeven(expHe,(complex_t){0,-delta});
  hBH0.getExponentialMPOeven(expHe_2,(complex_t){0,-delta*.5});
  hBH0.getExponentialMPOodd(expHo,(complex_t){0,-delta});

  //  expHe.exportMPOtext("checkExpHe.txt");
  //expHe_2.exportMPOtext("checkExpHe_2.txt");
  //expHo.exportMPOtext("checkExpHo.txt");
  //  exit(1);

  int cnt=0;
  *out<<"% Results of the time evolution "<<endl;
  *out<<"% cnt(delta="<<delta<<")\t ";
  for(int k=0;k<L;k++) *out<<"n("<<k<<")\t";
  *out<<"Ntot \t E \t S(L/2)"<<endl;

  MPS evolS(init);
  while(cnt<M){
    cout<<"Time evolution step "<<cnt<<endl;
    MPS aux(evolS);
    double errE,errO;
    contractor.optimize(expHe_2,aux,evolS,D,&errE);
    cout<<"After step "<<cnt<<" rel errE="<<errE<<endl;
    int ks=0;
    while((ks<stepBlk-1)&&(cnt<M-1)){
      aux=evolS;
      contractor.optimize(expHo,aux,evolS,D,&errO);
      //contractor.optimize(expHo,aux,evolS,D);
      aux=evolS;
      contractor.optimize(expHe,aux,evolS,D,&errE);
      //contractor.optimize(expHe,aux,evolS,D);
      ks++;cnt++;
      cout<<"After step "<<cnt-1<<" rel errE="<<errE
        <<", rel errO="<<errO<<endl;
    }
    aux=evolS;
    contractor.optimize(expHo,aux,evolS,D,&errO);
    //contractor.optimize(expHo,aux,evolS,D);
    aux=evolS;
    contractor.optimize(expHe_2,aux,evolS,D,&errE);
    //contractor.optimize(expHe_2,aux,evolS,D);
    cout<<"After step "<<cnt<<" rel errO="<<errO
    	<<" rel errE="<<errE<<endl;
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
    *out<<Ntot<<"\t"
	<<real(contractor.contract(evolS,hamil,evolS))<<"\t";
    // Now I would also like to compute the entropy between the right 
    // part and the left one
    double entr=contractor.getEntropy(evolS,center);
    *out<<entr<<endl;
  }


  out->close();


}
