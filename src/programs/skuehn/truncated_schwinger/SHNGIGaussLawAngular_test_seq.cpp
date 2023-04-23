/**
   \file SHNGIGaussLawAngular_test.cpp
   Full gauge invariant version of the Schwinger Model yet with truncated bosonic dimension (the link fileds are APPROXIMATED). This version allows a sequential run for different values of x and keeps the preliminary ground state in the memory to use it as a startstate. Therefore one is not dependent on the saved ground state form the disk (which is/was at that time problematic due to the errors with the filesystem)
   
  \param <N> (int) Number of Fermions
  \param <d> (int) Physical dimension of Bosons in between
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <xs> (double) Coupling constant to start with
  \param <xe> (double) Coupling constant to end with
  \param <dx> (double) Stepsize for the coupling constant
  \param <e> (double) Constant of the Operator \f$ \Theta(n)\f$
  \param <lambda> (double) Constant for the penalty-term in the Hamiltonian
   
  \author Stefan Kühn
  \date 03/07/2013

*/


#include <cstdlib>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <string>
//Lib which is needed for sleep command on *nix Systems
//#include "unistd.h"

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"

#include "SHNGIGaussLawAngular.h"

#define NUM_OF_PARAMS 9
#define FILENAME "SchwingerNGI.txt"

using namespace std;

bool file_exists(string filename);

int main(int argc,const char* argv[]){
  int cntr=0;
  if(argc < NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << " N d D mu x e la" << endl;
    cout << "N: " << "Number of spins" << endl
	 << "d: " << "Physical dimension of bosonic lattice sites" << endl
	 << "D: " << "Bond dimension of the MPS" << endl
	 << "µ: " << "Dimensionless constant" << endl
	 << "xs:" << "Coupling constant to start with" << endl
	 << "xe:" << "Coupling constant to end with" << endl
	 << "dx:" << "Stepsize for the coupling constant" << endl
	 << "e: " << "Constant for gauge field" << endl
	 << "la:" << "Constant factor for penalty term" << endl
	 << "ig:" << "True if inital guess for the ground state should be read form HDD (optional)" << endl
	 << "ac:" << "Accuracy for contractor (optional)" << endl;
    return -1;
  }
  int N=atoi(argv[++cntr]);
  double d=atof(argv[++cntr]);
  //double d2=atof(argv[++cntr]);
  double D=atof(argv[++cntr]);
  double mu=atof(argv[++cntr]);
  double xs=atof(argv[++cntr]);
  double xe=atof(argv[++cntr]);
  double deltax=atof(argv[++cntr]);
  double e=atof(argv[++cntr]);
  double lam_factor=atof(argv[++cntr]);
  bool ig = false;
  bool set_accuracy = false;
  double desired_accuracy;
  double lam;
  
  if(argc == NUM_OF_PARAMS + 2)
    ig = true;
  if(argc == NUM_OF_PARAMS + 3)
  {
    ig = true;
    set_accuracy = true;
    ++cntr;
    desired_accuracy = atof(argv[++cntr]);
  }
    
  
  if(argc==NUM_OF_PARAMS+1){
    cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "N  = " << N << endl
	<< "d  = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
	<< "dx = " << deltax << endl
	<< "e  = " << e << endl
	<< "la = " << lam_factor << endl
	<< "-----------" << endl << endl;
  }
  else if (argc==NUM_OF_PARAMS+2){
    cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "N  = " << N << endl
	<< "d  = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
	<< "dx = " << deltax << endl
	<< "e  = " << e << endl
	<< "la = " << lam_factor << endl
	<< "ig = " << ig << endl
	<< "-----------" << endl << endl;
  }
  else{
     cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "N  = " << N << endl
	<< "d  = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
	<< "dx = " << deltax << endl
	<< "e  = " << e << endl
	<< "la = " << lam_factor << endl
	<< "ig = " << ig << endl
	<< "acc= " << desired_accuracy << endl
	<< "-----------" << endl << endl;
  }
    
       
  //Measure time
  time_t start,end;
  time_t start_local,end_local;
  double readtime_gs=0.0,writetime_gs=0.0,writetime_gauss=0.0,writetime_spin=0.0,writetime_ln=0.0,writetime_total=0.0,writetime_data=0.0;
  
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
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  
    if(set_accuracy){
    cout << "Setting Convergence Tolareance to " << desired_accuracy<<endl<<endl;
    contractor.setConvTol(desired_accuracy);
  }
  
  
  MPS mps(L,bond_dim,phys_dims);
  mps.setRandomState();
  mps.gaugeCond('R',1);
  mps.gaugeCond('L',1);
  
  
  double E0 = 0.0;
  complex_t wert;
  
  ifstream inital_guess;
  inital_guess.open("GS_guess.dat");
  cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<        "                   %Starting to calculate Groundstate%                     " << endl
       <<        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
       
  //cout << "Starting to calculate Groundstate" << endl;
  if(inital_guess.is_open() && ig)
  {
    cout << "Found initial guess for the groundstate" << endl;
    time(&start_local);
    mps.importMPS("GS_guess.dat");
    time(&end_local);
    readtime_gs = difftime(end_local,start_local);
    cout << "Time for reading GS: " << readtime_gs << " s" << endl;
  }
  for(double x=xs; x<=xe; x += deltax)
  {
    lam = lam_factor*x;
    SHNGIGaussLawAngular Hsngi(N,d,D,mu,x,e,lam);    
    const MPO& hamil=Hsngi.getHMPO();  
    
    //Try to find the groundstate
    contractor.findGroundState(hamil,D,&E0,mps);   
    cout << "Groundstate energy: " << E0 << endl;
    cout << "Groundstate energy via Contractor " << contractor.contract(mps,hamil,mps) << endl;
    //mps.exportForMatlab ("Gs.m");
    
    //String for the name of the groundstate
    string GS_name="GS_";
    //String needed for the conversion
    ostringstream convert;  
    //Build the filename
    convert << N;
    GS_name = GS_name+"N"+convert.str();
    convert.str("");
    convert.clear();
    convert << x;
    GS_name = GS_name+"x"+convert.str(); 
    convert.str("");
    convert.clear();
    convert << mu;
    GS_name = GS_name+"mu"+convert.str()+".dat"; 
    time(&start_local);
    mps.exportMPS(GS_name.c_str());
    time(&end_local);
    writetime_gs = difftime(end_local,start_local);
    cout << "Time for writing GS: " << writetime_gs << " s" << endl;
    
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
    cout << endl;
    for(int i=1; i<L-1;i += 2){
      Hsngi.constructLopMPO(Ln,i);
      Lnvalue = contractor.contract(mps,Ln,mps);
      cout << "--> Expectation value L("<< (i+1)/2 << ") = " << Lnvalue<< endl;
      time(&start_local);
      Lnfile << Lnvalue.re << "\t";
      time(&end_local);
      writetime_ln += difftime(end_local,start_local);
    }
    time(&start_local);
    Lnfile << endl;
    Lnfile.close();
    time(&end_local);
    writetime_ln += difftime(end_local,start_local);
    
    cout << "Time for writing Ln: " << writetime_ln << " s" << endl;
    
    ofstream Gaussfile;
    Gaussfile.open("Gaussfile.txt", ios_base::app);
    complex_t Gaussvalue;
    MPO Gauss(1);
    cout << endl;
    for(int i=0; i<L;i += 2){
      Hsngi.constructGaussMPO(Gauss,i);
      Gaussvalue = contractor.contract(mps,Gauss,mps);
      cout << "--> Expectation value 0.5*[sigma_z+(-1)^"<< i/2+1 << "] = " << Gaussvalue<< endl;
      time(&start_local);
      Gaussfile << Gaussvalue.re << "\t";
      time(&end_local);
      writetime_gauss += difftime(end_local,start_local);
    }
    time(&start_local);
    Gaussfile << endl;
    Gaussfile.close();
    time(&end_local);
    writetime_gauss += difftime(end_local,start_local);
    
    cout << "Time for writing Gauss: " << writetime_gauss << " s" << endl;
    
    //Calculate the z-Component of the spins
    ofstream spinfile;
    spinfile.open("Spin.txt", ios_base::app);
    complex_t Zcomponent;
    MPO Sigmaz(1);
    cout << endl;
    for(int i=0; i<L;i += 2){
      Hsngi.constructsigmazMPO(Sigmaz,i);
      Zcomponent = contractor.contract(mps,Sigmaz,mps);
      cout << "--> Expectation value sigma_z("<< i/2+1 << ") = " << Zcomponent << endl;
      time(&start_local);
      spinfile << Zcomponent.re << "\t";
      time(&end_local);
      writetime_spin += difftime(end_local,start_local);
    }
    time(&start_local);
    spinfile << endl;
    spinfile.close(); 
    time(&end_local);
    writetime_spin += difftime(end_local,start_local);
    
    cout << "Time for writing Spin: " << writetime_spin << " s" << endl;
    
    //Calculate the condensate farction
    complex_t cfraction;
    MPO Condensate(1);
    cout << endl;
    Hsngi.constructCondensateMPO(Condensate);
    cfraction = contractor.contract(mps,Condensate,mps);
    cout << "Condensate: "<< cfraction << endl << endl;  
    
    //Ouput the results up to now (sometimes the computation of the exicted state(s) fails and the Contractor routine quits the program, therefore all results which were not stored would be lost). To avoid this, save the first part of the results now
    ofstream resfile;  
    cout << "Saving results..." << endl;
    if(file_exists(FILENAME))
    {
      time(&start_local);
      resfile.open(FILENAME,ios_base::app);  
      resfile << endl << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << x << "\t" << e << "\t" << lam << scientific << "\t" << E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" <<penalty_energy.re << "\t";
    }
    else
    {
      time(&start_local);
      resfile.open(FILENAME,ios_base::app);  
      resfile << "#N" << "\t" << "d" << "\t" << "D" << "\t" << "mu" << "\t" << "x" << "\t" << "e" << "\t" << "lambda" << "\t"<< "E0" << "\t" << "C0" << "\t" << "Cfrac" << "\t" << "E1" << "\t" << "C1" << "\tE_pen"<< endl;
      resfile << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << x << "\t" << e << "\t" << lam << scientific << "\t"<< E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" << penalty_energy.re << "\t";
    }
    resfile.close();
    time(&end_local);
    writetime_data = difftime(end_local,start_local);
    
    cout << "Time for writing Data: " << writetime_data << " s" << endl;
    writetime_total = writetime_data + writetime_gauss + writetime_gs + writetime_ln + writetime_spin;
    cout << "Time for writing in total: " << writetime_total << " s" << endl;
    cout << "Time on filesystem in total (reading + writing): " << writetime_total+readtime_gs << endl; 
  }
  
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