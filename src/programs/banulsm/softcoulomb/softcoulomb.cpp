
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#include "SoftCoulombHamiltonian.h"


/** softcoulomb tries to find the ground state of the soft Coulomb Hamiltonian
    \ref <SoftCoulombHamiltonian>
    // TODO: also try excited states

    Receives arguments:
    \param <L> (int) length of the chain
    \param <Delta> (double) lattice spacing
    \param <Z> (int) charge of the nucleus (located in the center 
                     of the system
    \param <a> (double) softening parameter \f$a\f$ of the Hamiltonian
    \param <M> (int) number of exponential terms used to approximate 
                     the Hamiltonian
    \param <Npen> (double) penalty term for number of particles
    \param <Ntot> (int) desired number of particles (if Npen==0, no effect)
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)
    \param <mpsfname> (char*) name of the file where to save the final MPS
    \param <initmpsfname> (char*) [optional] name of the file where
                         the initial MPS is saved

*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double Delta=atof(argv[++cntr]);
  int Z=atoi(argv[++cntr]);
  double a=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  double Npen=atof(argv[++cntr]);
  int Ntot=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);
  const char* mpsfname=argv[++cntr];
  const char* initmpsfname=0;
  if(argc>12) initmpsfname=argv[++cntr];

  // Create the Hamiltonian
  SoftCoulombHamiltonian H(L,Delta,Z,a,M,Npen,Ntot);
  // and get the MPO
  const MPO& hamil=H.getHMPO();
  cout<<"Initialized Hamiltonian"<<endl;
  //  hamil.exportForMatlab("HSC.m");

  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  cout<<"Initialized Contractor"<<endl;
  
  int d=4;
  MPS gsP(L,D,d);
  if(initmpsfname!=0)
    gsP.importMPS(initmpsfname);
  else
    gsP.setRandomState(); // the initial state, random -> not to the GS!!
  gsP.gaugeCond('R',1);
  gsP.gaugeCond('L',1);
  double lambda=0.;

  MPO Nmpo(L);
  H.getNumberMPO(Nmpo);
  cout<<"Before starting, expectation value is="
      <<contractor.contract(gsP,hamil,gsP)<<endl;
  cout<<", and norm="
      <<contractor.contract(gsP,gsP)<<endl;
  cout<<", and <N>="<<contractor.contract(gsP,Nmpo,gsP)<<endl;

  contractor.findGroundState(hamil,D,&lambda,gsP);

  cout<<"Ground state found with eigenvalue "<<lambda<<endl;
  complex_t normGS=contractor.contract(gsP,gsP);
  complex_t energyGS=contractor.contract(gsP,hamil,gsP);
  cout<<"After optimizing, expectation value is="
      <<energyGS<<", and norm="
      <<normGS<<endl;
  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% L\t D\t M\t Ze\t E(gs)\t N(gs)"<<endl;
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);
  *out<<L<<"\t"<<D<<"\t"<<M<<"\t"<<Z;
  *out<<"\t"<<energyGS/normGS;

  complex_t nGS=contractor.contract(gsP,Nmpo,gsP);
  cout<<"<N>="<<nGS<<endl;
  *out<<"\t"<<nGS/normGS<<endl;
  gsP.exportMPS(mpsfname);
  out->close();
  delete out;
}

