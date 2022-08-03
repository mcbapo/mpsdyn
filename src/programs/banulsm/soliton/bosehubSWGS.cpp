
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#include "BoseHubbardHamiltonian.h"
#include "BosonMPO.h"

/** bosehubSW tries to find the ground state of the Bose-Hubbard
    Hamiltonian     \ref <BoseHubbardHamiltonian>
    in a square well potential, with interactions which
    may be different in the three sectors (L,C,R).
    The limit if the three regions is specified by arguments
    \param<xL> and \param<xR> (both between 0 and L-1),
    which are the last sites included in the well.

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
    \param <VL> (double) parameter \f$\V0\f$ of the Hamiltonian 
    (constant) on the left part of the well
    \param <VR> (double) parameter \f$\V0\f$ of the Hamiltonian
    (constant) on the right part of the well
    \param <D> (int) maximum bond dimension
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
  double VL=atof(argv[++cntr]);
  double VR=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
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
  double* muis=new double[L];
  double* Vis=new double[L];
  for(int k=0;k<L;k++){
    tis[k]=t;
    Uis[k]=k<xL?UL:(k<=xR?UC:UR);
    muis[k]=mu;
    Vis[k]=k<xL?VL:(k<=xR?0:VR);
    cout<<"Site "<<k<<", U="<<Uis[k]<<",V="<<Vis[k]<<endl;
  }

  BoseHubbardHamiltonian hBH(L,N,tis,Uis,muis,Vis);
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
    *out<<k<<"\t"<<tis[k]<<"\t"<<Uis[k]<<"\t"<<muis[k]<<"\t"<<Vis[k]
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

  MPO parMPO(L);
  BosonMPO::getParityMPO(L,N,parMPO,xL,xR);
  complex_t parit=contractor.contract(gsP,parMPO,gsP);
  parit=parit/normGS;
  cout<<"Found parity operator in the GS: "<<parit<<endl;

  out->close();

}
