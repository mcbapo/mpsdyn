/**
   \file ZnPenalties.cpp
   At first, I used a simple penalty to monitor penalty energy and to push the variational calculations in the right subspace of the Hilbertspace. Later (25.03.2014) we realized, that a penalty which respects the cyclic nature of the Zn model might be more appropriate. To be able to reevaluate my old results with this new penalty, this little program is designed to read a state and compute the penalty energy in both ways.
   
  \param <N> (int) Number of Fermions
  \param <d> (int) Physical dimension of Bosons in between
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <x> (double) Coupling constant
  \param <e> (double) Constant of the Operator \f$ \Theta(n)\f$
  \param <lambda> (double) Constant for the penalty-term in the Hamiltonian
   
  \author Stefan Kühn
  \date 28/03/2014

*/


#include <cstdlib>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <string>
#include <algorithm>
//Lib which is needed for sleep command on *nix Systems
//#include "unistd.h"

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include "misc.h"

#include "KogutSusskindHamiltonian.h"

#define NUM_OF_PARAMS 8
#define FILENAME "KS.txt"
#define EPSILON 1.0E-9

using namespace std;

//bool file_exists(string filename);

int main(int argc,const char* argv[]){
  int cntr=0;
  if(argc < NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << " Model N d D mu x e la" << endl;
    cout << "N:     " << "Number of spins" << endl
	 << "d:     " << "Bosonic Hilbertspace dimension" << endl
	 << "D:     " << "Bond dimension of the MPS" << endl
	 << "µ:     " << "Dimensionless constant" << endl
	 << "x:     " << "Coupling constant" << endl
	 << "e:     " << "Constant for gauge field" << endl
	 << "la:    " << "Constant for penalty term" << endl
	 << "ig:    " << "Input state" << endl
	 << "ac:    " << "Accuracy for contractor (optional)" << endl;
    return -1;
  }
  int N=atoi(argv[++cntr]);
  int d=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  double mu=atof(argv[++cntr]);
  double x=atof(argv[++cntr]);
  double e=atof(argv[++cntr]);
  double lam=atof(argv[++cntr]);
  bool set_accuracy = false;
  double desired_accuracy;
  string ig_file = argv[++cntr];
  
  if(argc == NUM_OF_PARAMS + 2)
  {
    set_accuracy = true;
    desired_accuracy = atof(argv[++cntr]);
  }
  
    
  
  if(argc>=NUM_OF_PARAMS+1){
    cout << "Parameter: " << endl
	<< "------------------------------" << endl
	<< "N  = " << N << endl
	<< "d  = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "x  = " << x << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "ig = " << ig_file << endl;
  }
  if(argc>=NUM_OF_PARAMS+2)
     cout << "acc= " << desired_accuracy << endl;
  cout << "-------------------------------" << endl;
  
  //Now check, if the state which should be read exists, in case not, I can immediately kill the program and do not have to waste time constructing things
  if(!file_exists(ig_file))
  {
    cout << "Input state does not exist, program will be aborted..." << endl;
    exit(666);
  }
    
    
  int L=2*N-1; 
  
  //Create an array containing the physical dimensions (alternating between fermion and boson, ending with fermion)
  int *phys_dims = new int[L];
  int *bond_dim = new int[L-1];
  for(int i=0; i<L; ++i){
    phys_dims[i] = 2;
    if((i+1)<L){
      phys_dims[i+1] = d;
      i++;
    }
  }
  
  for (int i=0; i<L-1; i++){
    bond_dim[i] = D;
    //cout << "Eintrag Nr. " << i << " ist " << phys_dims[i] << endl;
  }
  
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  
  if(set_accuracy){
    cout << "Setting Convergence Tolerance to " << desired_accuracy<<endl;
    contractor.setConvTol(desired_accuracy);
  }
  
  //Prepare a MPS
  MPS mps(L,bond_dim,phys_dims);
  //Prepare the Hamiltonian, which is needed to get the MPOs for the penalty energies
  KogutSusskindHamiltonian Hsngi(N,d,D,mu,x,lam,Zncgl,0.0);
  
  
  mps.importMPS(ig_file.c_str());
  cout << "Imported state with D=" <<  mps.getBond() << " and L=" << mps.getLength() << endl ;

 
  complex_t penalty_energy,cgl_penalty_energy;
  MPO Penalty(1);
  MPO PenaltyCgl(1);
  Hsngi.constructPenaltyMPO(Penalty);
  Hsngi.constructCglPenaltyMPO(PenaltyCgl);
  penalty_energy = contractor.contract(mps,Penalty,mps); 
  cgl_penalty_energy = contractor.contract(mps,PenaltyCgl,mps);
  cout << endl << "Penalty energy:        "<< penalty_energy << endl;
  cout << "Cyclic Penalty energy: "<< cgl_penalty_energy << endl << endl;
  
  ofstream resfile;
  resfile.open("penalties.txt");
  resfile << penalty_energy << endl << cgl_penalty_energy << endl;
  resfile.close();
  
  //Free Memory
  delete[] phys_dims;
  delete[] bond_dim; 
  
  cout << "End of program" << endl;
  return 0;
}


