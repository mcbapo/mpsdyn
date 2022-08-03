/**
   \file third_test.cpp
   Full gauge invariant version of the Schwinger Model yet with truncated bosonic dimension (the link fileds are calulated exactly)
   
  \param <N> (int) Number of Fermions
  \param <d> (int) Physical dimension of Bosons in between
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <x> (double) Coupling constant
  \param <e> (double) Constant of the Operator \f$ \Theta(n)\f$
  \param <lambda> (double) Constant for the penalty-term in the Hamiltonian
   
  \author Stefan Kühn
  \date 13/03/2013

*/


#include <cstdlib>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <string>

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"

#include "SchwingerHamiltonianGIGaussLaw.h"

#define NUM_OF_PARAMS 7
#define FILENAME "SchwingerGI.txt"

using namespace std;

bool file_exists(string filename);

int main(int argc,const char* argv[]){
  int cntr=0;
  if(argc != NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << " N d D mu x e la" << endl;
    cout << "N: " << "Number of spins" << endl
	 << "d: " << "Physical dimension of bosonic lattice sites" << endl
	 << "D: " << "Bond dimension of the MPS" << endl
	 << "µ: " << "Dimensionless constant" << endl
	 << "x: " << "Coupling constant" << endl
	 << "e: " << "Constant for gauge field" << endl
	 << "la:" << "Constant for penalty term" << endl;
    return -1;
  }
  int N=atoi(argv[++cntr]);
  double d=atof(argv[++cntr]);
  //double d2=atof(argv[++cntr]);
  double D=atof(argv[++cntr]);
  double mu=atof(argv[++cntr]);
  double x=atof(argv[++cntr]);
  double e=atof(argv[++cntr]);
  double lam=atof(argv[++cntr]);
  
  cout << "Parameter: " << endl
       << "-----------" << endl
       << "N  = " << N << endl
       << "d  = " << d << endl
       << "D  = " << D << endl
       << "µ  = " << mu << endl
       << "x  = " << x << endl
       << "e  = " << e << endl
       << "la = " << lam << endl
       << "-----------" << endl << endl;
       
  //Measure time
  time_t start,end;
  
  time(&start);
  
  int L=2*N-1; 
  
  //Create an array containing the physical dimensions (alternating between fermion and boson, ending with fermion)
  int *phys_dims = new int[L];
  int *bond_dim = new int[L-1];
  for(int i=0; i<L; ++i){
    phys_dims[i] = 2;
    cout << "phys_dim[" << i << "]=" << phys_dims[i] << endl;
    if((i+1)<L){
      phys_dims[i+1] = d;
      cout << "phys_dim[" << i+1 << "]=" << phys_dims[i+1] << endl;
      i++;
    }
  }
  
  for (int i=0; i<L-1; i++){
    bond_dim[i] = D;
    //cout << "Eintrag Nr. " << i << " ist " << phys_dims[i] << endl;
  }
  
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  
  MPS mps(L,bond_dim,phys_dims);
  mps.setRandomState();
  mps.gaugeCond('R',1);
  mps.gaugeCond('L',1);
  SchwingerHamiltonianGIGaussLaw Hsngi(N,d,D,mu,x,e,lam);
  //cout << "SchwingerHamiltonianNonGaugeinvariant = " << Hsngi.getHMPO() << endl;
  //cout << "MPS davor " << mps << endl;
  
  //cout << "MPS is left-gauged " << mps.isGaugeL() << endl;
  //cout << "MPS is right-gauged " << mps.isGaugeR() << endl;
  
  //Try to find the groundstate
  const MPO& hamil=Hsngi.getHMPO();  
  double E0 = 0.0;
  complex_t wert;
  
  cout << "Starting to calculate Groundstate" << endl;
  contractor.findGroundState(hamil,D,&E0,mps);   
  cout << "Groundstate energy: " << E0 << endl;
  cout << "Groundstate energy via Contractor " << contractor.contract(mps,hamil,mps) << endl;
  //mps.exportForMatlab ("Gs.m");
  
  //Calculate the contribution of the penalty term
  complex_t penalty_energy;
  MPO Penalty(1);
  Hsngi.constructPenaltyMPO(Penalty);
  penalty_energy = contractor.contract(mps,Penalty,mps);
  cout << "Penalty energy: "<< penalty_energy << endl;
  
  cout << "Ratio Penalty in percent/Groundstate energy: " << penalty_energy.re/E0*100.0<<" %" << endl<<endl; 
  
  //Try to calculate the charge
  const MPO& chargeop=Hsngi.getChargeMPO();  
  complex_t charge_0=contractor.contract(mps,chargeop,mps);
  cout << "Charge Groundstate " << charge_0.re << endl;
  
  
  //Calculate the therms for the Gauss Law and see if it is fullfilled
  ofstream Lnfile;
  Lnfile.open("Ln.txt", ios_base::app);
  complex_t Lnvalue;
  MPO Ln(1);
  for(int i=1; i<L-1;i += 2){
    Hsngi.constructLopMPO(Ln,i);
    Lnvalue = contractor.contract(mps,Ln,mps);
    cout << "--> Expectation value L("<< i << ") = " << Lnvalue<< endl;
    Lnfile << Lnvalue.re << "\t";
  }
  Lnfile << endl;
  Lnfile.close();
  
  ofstream Gaussfile;
  Gaussfile.open("Gaussfile.txt", ios_base::app);
  complex_t Gaussvalue;
  MPO Gauss(1);
  for(int i=0; i<L;i += 2){
    Hsngi.constructGaussMPO(Gauss,i);
    Gaussvalue = contractor.contract(mps,Gauss,mps);
    cout << "--> Expectation value 0.5*[sigma_z+(-1)^"<< i/2+1 << "] = " << Gaussvalue<< endl;
    Gaussfile << Gaussvalue.re << "\t";
  }
  Gaussfile << endl;
  Gaussfile.close();
  
  //Calculate the z-Component of the spins
  ofstream spinfile;
  spinfile.open("Spin.txt", ios_base::app);
  complex_t Zcomponent;
  MPO Sigmaz(1);
  for(int i=0; i<L;i += 2){
    Hsngi.constructsigmazMPO(Sigmaz,i);
    Zcomponent = contractor.contract(mps,Sigmaz,mps);
    cout << "--> Expectation value sigma_z("<< i/2+1 << ") = " << Zcomponent << endl << endl;
    spinfile << Zcomponent.re << "\t";
  }
  spinfile << endl;
  spinfile.close(); 
  
  //Calculate the condensate farction
  complex_t cfraction;
  MPO Condensate(1);
  Hsngi.constructCondensateMPO(Condensate);
  cfraction = contractor.contract(mps,Condensate,mps);
  cout << "Condensate: "<< cfraction << endl << endl;  
  
  //Ouput the results up to now (sometimes the computation of the exicted state(s) fails and the Contractor routine quits the program, therefore all results which were not stored would be lost). To avoid this, save the first part of the results now
  ofstream resfile;  
  cout << "Saving results..." << endl;
  if(file_exists(FILENAME))
  {
    resfile.open(FILENAME,ios_base::app);  
    resfile << endl << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << x << "\t" << e << "\t" << lam << "\t" << E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" <<penalty_energy.re << "\t";
  }
  else
  {
    resfile.open(FILENAME,ios_base::app);  
    resfile << "#N" << "\t" << "d" << "\t" << "D" << "\t" << "mu" << "\t" << "x" << "\t" << "e" << "\t" << "lambda" << "\t"<< "E0" << "\t" << "C0" << "\t" << "Cfrac" << "\t" << "E1" << "\t" << "C1" << "\tE_pen"<< endl;
    resfile << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << x << "\t" << e << "\t" << lam << "\t"<< E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" << penalty_energy.re << "\t";
  }
  resfile.close();
  
  
  //Get the excited state(s) 
  /*
  vector<MPS*> excitedstates;
  double lambdak;
  MPS first_ex(L,bond_dim,phys_dims);
  
  excitedstates.push_back(&mps);
  
  cout << endl << "Computing excited states" << endl;
  contractor.findNextExcitedState(hamil,D,excitedstates,&lambdak,first_ex); 
  cout << "Excited state energy: " << lambdak << endl<< endl;
  complex_t charge_1=contractor.contract(first_ex,chargeop,first_ex);
  
  //Save results for the exicted state(s)
  cout << "Saving results..." << endl;
  resfile.open(FILENAME,ios_base::app);  
  resfile << lambdak << "\t" << charge_1.re ;
  resfile.close();*/  
  
  //Free Memory
  delete[] phys_dims;
  delete[] bond_dim; 
  
 
  
  //End Time
  time(&end);
  cout << "Run time: " << difftime(end,start) <<"s" <<  endl;  
  cout << "End of program" << endl;
  return 0;
}

bool file_exists(string filename)
{
  ifstream file;
  file.open(filename.c_str());
  if(file.good())
    return true;
  else
    return false;
}