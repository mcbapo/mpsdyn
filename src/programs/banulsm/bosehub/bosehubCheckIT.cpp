
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#include "BoseHubbardHamiltonian.h"
#include "BosonMPO.h"


/** bosehubCheckIT finds the ground state of the Bose-Hubbard
    Hamiltonian \ref <BoseHubbardHamiltonian> in a certain potential
    well, with the expression \f[V_0(x)=-A sech^2(B(x-x_0))\f] and
    initial interactions equal to \f$U\f$.  To check the
    implementation, GS is found first with the direct search and then
    with imaginary time evolution.

    Receives arguments:
    \param <L> (int) total length of the chain
    \param <N> (int) cutoff for the occupation number
    \param <t> (double) parameter \f$t\f$ of the Hamiltonian
    \param <U> (double) parameter \f$U\f$ of the Hamiltonian
    \param <mu> (double) parameter \f$\mu\f$ of the Hamiltonian
    \param <A> (double) parameter \f$A\f$ of the initial potential
    \param <B> (double) parameter \f$B\f$ of the initial potential
    \param <x0> (double) parameter \f$x0\f$ of the initial potential 
                            (center)
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
  double A=atof(argv[++cntr]);
  double B=atof(argv[++cntr]);
  int x0=atoi(argv[++cntr]);
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

  double* tis=new double[L];
  double* Uis0=new double[L];
  double* muis0=new double[L];
  double* Vis0=new double[L];
  for(int k=0;k<L;k++){
    tis[k]=t;
    Uis0[k]=U;
    // another test!
    //    Uis[k]=k<xL?UL:(k<=xR?UC:UR*(1-exp(-(k-x0)*(k-x0)/(a*a))));
#ifdef ONLYGS 
    //Uis0[k]=Uis[k]; // to check!
#endif
    muis0[k]=mu;
    double aux=cosh(B*(k-x0));
    Vis0[k]=A/(aux*aux);
    cout<<"Site "<<k<<", U="<<Uis0[k]<<",V="<<Vis0[k]<<endl;
  }

  BoseHubbardHamiltonian hBH0(L,N,tis,Uis0,muis0,Vis0);
  const MPO& hamil=hBH0.getHMPO();

  cout<<"Before starting, expectation value is="
      <<contractor.contract(gsP,hamil,gsP)<<", and norm="
      <<contractor.contract(gsP,gsP)<<endl;

  contractor.findGroundState(hamil,D,&lambda,gsP);

  cout<<"Ground state found (direct minim) with eigenvalue "<<lambda<<endl;
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
    *out<<k<<"\t"<<tis[k]<<"\t"<<Uis0[k]<<"\t"<<muis0[k]<<"\t"<<Vis0[k]
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
  *out<<"% Delta H not relevant now"<<endl;
  cout<<"Now searching the GS by imaginary time evolution"<<endl;
  // Get the exponential MPOs
  MPO expHe(L),expHo(L),expHe_2(L);
  hBH0.getExponentialMPOeven(expHe,(complex_t){-delta,0});
  hBH0.getExponentialMPOeven(expHe_2,(complex_t){-delta/2,0});
  hBH0.getExponentialMPOodd(expHo,(complex_t){-delta,0});

  int cnt=0;
  *out<<"% Results of the time evolution "<<endl;
  *out<<"% cnt(delta="<<delta<<")\t ";
  for(int k=0;k<L;k++) *out<<"n("<<k<<")\t";
  *out<<"Ntot"<<"\t Ntot(MPO)\t <H>"<<endl;

  MPS evolS(gsP);
  evolS.setRandomState();
  while(cnt<M){
    cout<<"(Imaginary) Time evolution step "<<cnt<<endl;
    MPS aux(evolS);
    //double errE,errO;
    contractor.optimize(expHe_2,aux,evolS,D);
    //cout<<"After step "<<cnt<<" rel errE="<<errE<<endl;
    int ks=0;
    while((ks<stepBlk-1)&&(cnt<M-1)){
      aux=evolS;
      //contractor.optimize(expHo,aux,evolS,D,&errO);
      contractor.optimize(expHo,aux,evolS,D);
      aux=evolS;
      //      contractor.optimize(expHe,aux,evolS,D,&errE);
      contractor.optimize(expHe,aux,evolS,D);
      ks++;cnt++;
      //cout<<"After step "<<cnt-1<<" rel errE="<<errE
      //  <<", rel errO="<<errO<<endl;
    }
    aux=evolS;
    //contractor.optimize(expHo,aux,evolS,D,&errO);
    contractor.optimize(expHo,aux,evolS,D);
    aux=evolS;
    //contractor.optimize(expHo,aux,evolS,D,&errO);
    contractor.optimize(expHe_2,aux,evolS,D);
    //cout<<"After step "<<cnt<<" rel errO="<<errO
    //	<<" rel errE="<<errE<<endl;
    cnt++;
    double time=delta*cnt;
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
    *out<<Ntot<<"\t"<<real(contractor.contract(evolS,Nmpo,evolS))
	<<"\t"<<real(contractor.contract(evolS,hamil,evolS))
	<<"\t"<<contractor.contract(evolS,gsP)<<endl;
    // Now I would also like to compute the entropy between the right 
    // part and the left one
    //    double entr=contractor.getEntropy(evolS,xR);
    //*out<<entr<<endl;
  }

  // At the end I compute the overlap between gS as IT and variational
  out->close();

}
