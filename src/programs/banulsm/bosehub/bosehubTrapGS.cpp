
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include "BosonMPO.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "BoseHubbardHamiltonian.h"

using namespace shrt;

/** bosehub tries to find the ground state of the Bose-Hubbard Hamiltonian
    \ref <BoseHubbardHamiltonian> inside a trap.
    Receives arguments:
    \param <L> (int) length of the chain
    \param <N> (int) cutoff for the occupation number
    \param <t> (double) parameter \f$t\f$ of the Hamiltonian
    \param <U> (double) parameter \f$U\f$ of the Hamiltonian
    \param <mu> (double) parameter \f$\mu\f$ of the Hamiltonian
    \param <V> (double) parameter \f$V_0\f$ of the Hamiltonian 
            (asume center in the middle, (L+1)/2 )
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)

*/

void getMPOs(int L,int N,int i,int j,MPO& parityString,MPO& SiSj);

void computeCorrelations(int Nmax,const MPS& state,const char* filename,
			 const MPS& idState,const MPO& unfoldedH);

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
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  double x0=(L+1)*.5;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% L="<<L<<", J/U="<<t<<", mu="<<mu<<", V0="<<V<<"\n"
	<<"% Nmax="<<N<<", D="<<D<<", "<<endl;
    *out<<"% EnergyGS= ";
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

  const MPO& hamil=hBH.getHMPO();

  cout<<"Before starting, expectation value is="
      <<contractor.contract(gsP,hamil,gsP)<<", and norm="
      <<contractor.contract(gsP,gsP)<<endl;

  contractor.findGroundState(hamil,D,&lambda,gsP);

  cout<<"Ground state found with eigenvalue "<<lambda<<endl;

  complex_t normGS=contractor.contract(gsP,gsP);
  complex_t energyGS=contractor.contract(gsP,hamil,gsP);
  MPO Nmpo(L);
  BosonMPO::getNumberMPO(L,N,Nmpo);
  complex_t totalN=contractor.contract(gsP,Nmpo,gsP);
  cout<<"After optimizing, expectation value is="
      <<energyGS<<", and norm="
      <<normGS<<endl;
  *out<<energyGS/normGS<<", Ntot="<<totalN/normGS
      <<endl;

  mwArray identPhys=identityMatrix(d);
  mwArray locN(Indices(d,d));locN.fillWithZero();
  mwArray expN(Indices(d,d));expN.fillWithZero();
  for(int k=0;k<d;k++){
    locN.setElement(k,0,Indices(k,k)); // (N_k)
    expN.setElement(exp((k)*M_PIl*I_c),Indices(k,k)); // exp(i pi N)
  }
  identPhys.reshape(Indices(d,1,d,1));
  locN.reshape(Indices(d,1,d,1));
  expN.reshape(Indices(d,1,d,1));
  Operator* idOpPhys=new Operator(identPhys);
  Operator* Nop=new Operator(locN);
  Operator* Sop=new Operator(expN);

  // N operator for a site
  *out<<"% i\t <n_i>"<<endl;
  for(int i=0;i<L;i++){
    Nmpo.setOp(i,idOpPhys);
  } // Now it is the identity operator
  for(int k=0;k<L;k++){
    Nmpo.setOp(k,Nop);
    complex_t Nvalk=contractor.contract(gsP,Nmpo,gsP);
    Nvalk=Nvalk/normGS;
    *out<<k<<"\t"<<real(Nvalk)<<"\t"<<imag(Nvalk)<<endl;
    Nmpo.setOp(k,idOpPhys);
  }

  // Now the correlations Si Sj

  *out<<"\n\n"<<endl;
  *out<<"% Two body correlations, S_i=exp(i pi n_i)"<<endl;
  *out<<"% i\t j \t <S_i S_j>"<<endl;
  for(int i=0;i<L;i++){
    Nmpo.setOp(i,Sop);
    for(int j=i;j<L;j++){
      if(j!=i) Nmpo.setOp(j,Sop);
      complex_t Sij=contractor.contract(gsP,Nmpo,gsP);
      Sij=Sij/normGS;
      *out<<i<<"\t"<<j<<"\t"<<real(Sij)<<"\t"<<imag(Sij)<<endl;
      if(j!=i) Nmpo.setOp(j,idOpPhys);
    }
    Nmpo.setOp(i,idOpPhys);
  }

  // And the parity string operator
  *out<<"\n\n"<<endl;
  *out<<"% Parity string operator Pij:=<S_i S_[i+1]...S_j>" <<endl;
  *out<<"% i\t j \t <P_ij>"<<endl;

  for(int i=0;i<L;i++){
    for(int j=i+1;j<L;j++){
      BosonMPO::getParityMPO(L,N,Nmpo,i,j);
      complex_t Pij=contractor.contract(gsP,Nmpo,gsP);
      Pij=Pij/normGS;
      *out<<i<<"\t"<<j<<"\t"<<real(Pij)<<"\t"<<imag(Pij)<<endl;
    }
  }

}
