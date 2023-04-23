
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#include "BoseHubbardHamiltonian.h"
#include "BosonMPO.h"

/** bosehub tries to find the ground state of the Bose-Hubbard Hamiltonian
    \ref <BoseHubbardHamiltonian>

    Receives arguments:
    \param <L> (int) length of the chain
    \param <N> (int) cutoff for the occupation number
    \param <t> (double) parameter \f$t\f$ of the Hamiltonian
    \param <U> (double) parameter \f$U\f$ of the Hamiltonian
    \param <mu> (double) parameter \f$\mu\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
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
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"# L="<<L<<", N="<<N<<", D="<<D<<", EnergyGS= ";
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
  BoseHubbardHamiltonian hBH(L,N,t,U,mu);
  const MPO& hamil=hBH.getHMPO();

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

  // Finally, compute total number of particles
  // A very inefficient way :

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
  // and at the same time compute total number of particles
  double Ntot=0.;
  *out<<"% k\t <n(k)>"<<endl;
  for(int k=0;k<L;k++){
    mpoN.setOp(k,&nrk,0);
    complex_t nrk=contractor.contract(gsP,mpoN,gsP);
    nrk=nrk/normGS;
    *out<<k<<"\t"<<nrk<<endl;
    cout<<k<<"\t"<<nrk<<endl;
    // restore the identity
    mpoN.setOp(k,&ident,0);
    Ntot=Ntot+real(nrk);
  }
  *out<<"% Ntot="<<Ntot<<endl;
  cout<<"Total energy="<<energyGS/normGS<<", and number of particles N="<<Ntot<<endl;

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

  double nbar=real(totalN/L);

  *out<<"% Parity calculations"<<endl;
  *out<<"i\t k\t P(i,k)"<<endl;
  MPO parMPO(L);
  int posCent=L/2; // central site
  //complex_t factorN=exp(-(totalN*(1./L))*M_PIl*I_c);
  for(int k=1;k<L/2;k++){
    BosonMPO::getParityMPO(L,N,parMPO,posCent-k,posCent+k,nbar);
    complex_t parit=contractor.contract(gsP,parMPO,gsP);
    //parit=factorN*parit/normGS;
    parit=parit/normGS;
    *out<<posCent-k<<"\t"<<posCent+k<<"\t"<<parit<<endl;
    cout<<"Parity("<<posCent-k<<","<<posCent+k<<")="<<parit<<endl;
  }

  /**
  *out<<"% Fr the test MPS"<<endl;
  // TEST with given MPS
  MPS Psi(L,1,d);
  mwArray Atest(shrt::Indices(d,1,1));Atest.fillWithZero();
  mwArray Atest0(shrt::Indices(d,1,1));Atest0.fillWithZero();
  Atest.setElement(1.,0.,shrt::Indices(3,0,0)); // 2 per site
  Atest0.setElement(1.,0.,shrt::Indices(0,0,0)); // 0 per site
  for(int k=0;k<L;k++)
    if(k%3==0)
      Psi.setA(k,Atest);
    else
      Psi.setA(k,Atest0);
  cout<<"Constructed Psi"<<endl;
  cout<<"with norm "<<contractor.contract(Psi,Psi)<<endl;
  cout<<"And nr of particles="<<contractor.contract(Psi,Nmpo,Psi)<<endl;
  for(int k=1;k<L/2;k++){
    BosonMPO::getParityMPO(L,N,parMPO,posCent-k,posCent+k,1);
    complex_t parit=contractor.contract(Psi,parMPO,Psi);
    //parit=factorN*parit/normGS;
    //parit=parit/normGS;
    *out<<posCent-k<<"\t"<<posCent+k<<"\t"<<parit<<endl;
    cout<<"Parity("<<posCent-k<<","<<posCent+k<<")="<<parit<<endl;
  }
  */


  out->close();

}
