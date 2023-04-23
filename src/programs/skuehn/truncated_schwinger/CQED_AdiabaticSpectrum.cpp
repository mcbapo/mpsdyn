/**
  \file CQED_AdiabaticSpectrum.cpp
  Basically the same file as CQED_AdiabaticEvolution but  with some simplifications. The goal is to find out why at some point the simulation fails miserably for some parameters and for others not
   
  \param <N> (int) Number of Fermions
  \param <d> (int) Physical dimension of Bosons in between
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <x> (double) Coupling constant
  \param <e> (double) Constant of the Operator \f$ \Theta(n)\f$
  \param <lambda> (double) Constant for the penalty-term in the Hamiltonian
   
  \author Stefan Kühn
  \date 02/08/2013

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
#include "Indices.h"

#include "CQED.h"
#include "ramp_xvalue.h"

#define NUM_OF_PARAMS 9
#define FILENAME "cQED.txt"
#define RESOLUTION 100
#define EPS 1E-10

using namespace std;
using namespace shrt;

bool file_exists(string filename);
void print_MPS(MPS);

int main(int argc,const char* argv[]){
  int cntr=0;
  
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << " N d D mu x e la" << endl;
    cout << "N: " << "Number of spins" << endl
	 << "N0:" << "Total number of bosons" << endl
	 << "D: " << "Bond dimension of the MPS" << endl
	 << "µ: " << "Dimensionless constant" << endl
	 << "x: " << "Coupling constant" << endl
	 << "e: " << "Constant for gauge field" << endl
	 << "la:" << "Constant for penalty term" << endl
	 << "dt:" << "Timestep size" << endl
	 << "Nt:" << "Number of timesteps which should be calculated" << endl
	 << "ig:" << "True if inital guess for the ground state should be read form HDD (optional)" << endl
	 << "ac:" << "Accuracy for contractor (optional)" << endl;
    return -1;
  }
  
  //Number of spin sites
  int N=atoi(argv[++cntr]);
  //Total number fo bosons living on a link
  double d=atof(argv[++cntr]);
  //double d2=atof(argv[++cntr]);
  //Bond dimension
  double D=atof(argv[++cntr]);
  //Mass
  double mu=atof(argv[++cntr]);
  //Coupling constant
  double x=atof(argv[++cntr]);
  //Constant which can be used to manipulate the Hamiltonian
  double e=atof(argv[++cntr]);
  //Strength of the penalty
  double lam=atof(argv[++cntr]);
  //Size of timestep
  double dt=atof(argv[++cntr]);
  //Total number of timesteps
  int num_of_steps=atoi(argv[++cntr]);
  //Desired order of the evolution (by default 1)
  int evolution_order=1;
  //Inital state which should be used?
  bool ig = false;
  //Set the accuracy of contractor?
  bool set_accuracy = false;
  //Variable for the desired accuracy in case it should be set
  double desired_accuracy;
  //Name for file for initial state
  string ig_file;
  //Variable needed to verify if file for inital guess should have dummy name or a custom name
  size_t found_dat;
  
  //Only file for inital guess is given aditionally to all other parameters
  if(argc == NUM_OF_PARAMS + 2)
  {
    ig = true;
    ig_file = argv[++cntr];
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file: " << ig_file << endl;
  }
  //Custom accuracy for contractor and file for inital state
  if(argc == NUM_OF_PARAMS + 3)
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file " << ig_file << endl;
  }
    
  //Print parameter which were set
  if(argc==NUM_OF_PARAMS+1){
    cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "N  = " << N << endl
	<< "d  = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "x  = " << x << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "dt = " << dt << endl
	<< "Nt = " << num_of_steps << endl
	<< "-----------" << endl << endl;
  }
  else if (argc==NUM_OF_PARAMS+2){
    cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "N  = " << N << endl
	<< "d  = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "x  = " << x << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "dt = " << dt << endl
	<< "Nt = " << num_of_steps << endl
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
	<< "x  = " << x << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "dt = " << dt << endl
	<< "Nt = " << num_of_steps << endl
	<< "ig = " << ig << endl
	<< "acc= " << desired_accuracy << endl
	<< "-----------" << endl << endl;
  }
    
       
  //Measure time
  time_t start,end;
  time_t start_local,end_local;
  double readtime_gs=0.0,writetime_gs=0.0,writetime_gauss=0.0,writetime_spin=0.0,writetime_ln=0.0,writetime_total=0.0,writetime_data=0.0;
  
  time(&start);
  
  //Total number fo sites (fermionic and bosonic)
  int L=2*N-1; 
  
  //Create an array containing the physical dimensions (alternating between fermion and boson, ending with fermion)
  int *phys_dims = new int[L];
  int *bond_dim = new int[L-1];
  for(int i=0; i<L; ++i){
    phys_dims[i] = 2;
    cout << "phys_dim[" << i << "]=" << phys_dims[i] << endl;
    if((i+1)<L){
      phys_dims[i+1] = d+1;
      cout << "phys_dim[" << i+1 << "]=" << phys_dims[i+1] << endl;
      i++;
    }
  }
  for(int i=0; i<(L-1); i++)
  {
    bond_dim[i] = D;
  }
  
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  //Set accuracy of contractor if desired
  if(set_accuracy){
    cout << "Setting Convergence Tolareance to " << desired_accuracy<<endl<<endl;
    contractor.setConvTol(desired_accuracy);
  }
  
  
  
  //Pointer to MPS to hold the inital state which will be evolved in time
  MPS *mps;
  
  //Get the Hamiltonian object for the cQED model
  CQED Hsngi(N,d,D,mu,x,e,lam);  
  //Get Hamiltonian MPO
  const MPO& hamil = Hsngi.getHMPO();  
  //MPO for the odd and even part of the evolution operator
  MPO exp_even(1),exp_odd(1);
  //Energy values
  double E0 = 0.0; 
  double E0_direct = 0.0; 
  //Auxiliary MPS to store an old version before a timestep to be able to monitor changes
  MPS aux,mps_old;
  //Variable for the current energy value during time evolution
  complex_t curr_energy;
  //Variable for the current norm of the state during time evolution
  complex_t curr_norm;
  //Variable for the current overlap with the wavefunction of the previous step
  complex_t curr_overlap;
  //Variable for the value of the penalty energy
  complex_t penalty_energy;
  //MPO for the penalty term
  MPO Penalty(1);
  //Variable for the condensate
  complex_t cfraction;
  //MPO to determine the condensate
  MPO Condensate(1);
  //Variable for the current x value
  double xvalue=0.0;
  //File for energy values and other observables during the evolution
  ofstream energy_file;
  //Filepointer for the inital guess
  ifstream inital_guess;
  
  //Try to open file for inital state
  inital_guess.open(ig_file.c_str());
  //In case file for inital state exists read it and adjust bond dimension, otherwise start with predefined inital state
  if(inital_guess.is_open() && ig)
  {
    //Case in which an inital guess exists and will be read
    cout << "Found initial guess for the groundstate" << endl;
    time(&start_local);
    mps->importMPS(ig_file.c_str());
    cout << "Bond dimension after import: "<< mps->getBond() << endl;
    cout << "Length after import: " << mps->getLength() << endl;
    //In case an intial guess with bond dimension smaller than set is imported it has to be expanded again
    if(mps->getBond() < D)
    {
      //Increase bond dimension again, in case the inital guess has a smaller bond dimension
      //mps.increaseBondDimension(D);
      //Increase bond dimension again, in case the inital guess has a smaller bond dimension
      //This time with a little noise, maybe better to avoid stability issues?
      mps->increaseBondDimensionWithNoise(D,1E-2);
      cout << "Bond dimension after increase: "<< mps->getBond() << endl;
      
    }
    else if(mps->getBond() > D)
    {
      cout << "Error, the bond dimension of the inital guess is greater than the maximum bond dimension D=" << D << " set by the user" << endl;
      cout << "Bond dimension cannot be shrunk" << endl;
      exit(1);
    }
    time(&end_local);
    readtime_gs = difftime(end_local,start_local);
    cout << "Time for reading GS: " << readtime_gs << " s" << endl;
  }
  else
  {
    //Prepare some inital state
    mps = Hsngi.constructInitalMPS();
    //mps->setRandomState();
    //mps = Hsngi.constructInitalMPS(is_downup);
    //mps = Hsngi.constructInitalMPS(is_upup);
    //mps = Hsngi.constructInitalMPS(is_downdown);    
  } 
  
  //Determine the inital energy
  cout << "Energy before time evolution " << contractor.contract(*mps,hamil,*mps) << endl;  
  //Determine condensate before evolution
  Hsngi.constructCondensateMPO(Condensate);
  cfraction = contractor.contract(*mps,Condensate,*mps);  
  cout <<"Intial condensate: " << cfraction << " with x-value " << Hsngi.get_x() << endl;
  
  //Try to find the groundstate by time evolution of an intial state 
  
  
  
  //Second order time evolution (but this time a try with some more optimization)
  cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
      <<        "                  %Starting second order time evolution%                     " << endl
      <<        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  //Construct the two additional MPOs for the halfsteps
  MPO exp_even_half(1),exp_odd_half(1);
  
  energy_file.open("efile3.txt");
  energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\t1. Ex\tTrue Ex\tOverlap GS\tOverlap Ex (resolved every "<< RESOLUTION << " steps)" << endl;
  
  double xold;
  bool is_complete=false;
  
  //Stuff for excited states
  vector<MPS*> excitedstates;
  double lambdak=0.0,lambdak_direct=0.0;
  MPS first_ex(L,bond_dim,phys_dims);
  MPS gs_direct(L,bond_dim,phys_dims);
  MPS first_ex_direct(L,bond_dim,phys_dims);
  
  first_ex_direct.setProductState(p_zero);
  gs_direct.setProductState(p_zero);
  
  
  //Initial half-step
  xvalue =0.0;
  xold = xvalue;
  Hsngi.set_x(xvalue);
  Hsngi.constructUoddMPO(exp_odd_half,dt/2.0);
  contractor.optimize(exp_odd_half,*mps,aux,D);
  //Restore norm
  aux.setNormFact(1.0);
  
  //Actual time evolution
  //Note: Due to the fact that I group the odd timesteps in the loop and have do the correction at the end, this time the loop only runs to i<num_of_steps to have the desired number of timesteps at the end
  for(int i=1; i<num_of_steps; ++i)
  {
    mps_old = *mps;
    
    //Complete step for even part (make sure that inital guess is empty)
    Hsngi.constructUevenMPO(exp_even,dt);       
    mps->clear();
    //cout << "exp(-t U_even) " ;
    contractor.optimize(exp_even,aux,*mps,D);
    //Restore norm
    mps->setNormFact(1.0);
    
    //Get next xvalue and make a complete step with the mixture of the old xvalue und the new one
    xvalue = ramp_xvalue(0.0,x,num_of_steps-1,i);  
    //cout << "x = " << xvalue << ", xold = " << xold << endl;
    Hsngi.constructUoddMPO(exp_odd,dt,xold,xvalue);      
    aux.clear();
    //cout << "[exp(-t/2 U_odd) exp(-t/2 U_odd)] " ;
    contractor.optimize(exp_odd,*mps,aux,D);
    //Restore norm
    aux.setNormFact(1.0);
    
    //Update Hamiltonian
    Hsngi.set_x(xvalue);
    
    //Set x old for next iteration (needed to construct Ueven)
    xold = xvalue;
	  
    //make sure that the state is normalizes again, therefore I just apply the gauge condition and normalize it with this routine (gauging is not necessary but this is the simplest way to do it)
    //mps->gaugeCond('L',true);
    //aux.setNormFact(1.0);
    aux.gaugeCond('L',true);
    
    if(((i+1)%RESOLUTION)==0)
    {
      //Complete the step
      Hsngi.constructUevenMPO(exp_even,dt);
      mps->clear();
      contractor.optimize(exp_even,aux,*mps,D);
      //Restore norm
      mps->setNormFact(1.0);
      
      Hsngi.constructUoddMPO(exp_odd_half,dt/2.0);	
      aux.clear();
      contractor.optimize(exp_odd_half,*mps,aux,D);	
      //Restore norm
      aux.setNormFact(1.0);
      
      //cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
      
      //Save result in mps
      *mps = aux;
      
      //Get the values which should be monitored	
      Hsngi.updateHMPO();
      Hsngi.constructPenaltyMPO(Penalty);
      
      //Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
      curr_overlap = contractor.contract(mps_old,*mps);
      cout << "Overlap " << curr_overlap << "\t";
    
      //curr_energy = contractor.contract(*mps,hamil2,*mps);
      curr_energy = contractor.contract(*mps,hamil,*mps);
      curr_norm = contractor.contract(*mps,*mps);
      
      penalty_energy = contractor.contract(*mps,Penalty,*mps);
      
      //contractor.findGroundState(hamil2,D,&E0,aux);
      //contractor.findGroundState(hamil,D,&E0,aux);
      
      Hsngi.constructCondensateMPO(Condensate);
      cfraction = contractor.contract(*mps,Condensate,*mps);
      
      
      if (Hsngi.get_x() >= 0.5)
      {
	contractor.findGroundState(hamil,D,&E0,gs_direct);    
	
	cout << "Getting excited state" << endl;
	excitedstates.clear();
	excitedstates.push_back(mps);
	contractor.findNextExcitedState(hamil,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;
	
	cout << "Getting excited state from exact GS" << endl;
	excitedstates.clear();
	excitedstates.push_back(&gs_direct);
	contractor.findNextExcitedState(hamil,D,excitedstates,&lambdak_direct,first_ex_direct); 
	cout << "Excited state energy (from direct GS): " << lambdak_direct << endl<< endl;
      }
      
      cout<< "Norm in step " << (i+1) << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
      //energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
      energy_file << (i+1)*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<<"\t" << Hsngi.get_x() << "\t" << lambdak << "\t" << lambdak_direct <<" \t" << abs(contractor.contract(gs_direct,*mps)) << "\t" << abs(contractor.contract(first_ex_direct,first_ex)) << endl;
      
      
      //Do a halfstep again to be in the be able to use the grouped-together time evolution again
      //Raise the value of i, since I start an an additional timestep to get the values at the right time
      i++;
      //It can be, that I know already finished all timesteps, I had to calculate, so the look is not supposed to run any more and I don't need the following additional halfstep. For this case, exit the look and set a flag, that the correction after the loop to get a full timestep is not necessary any more
      if(i==num_of_steps)
      {
	is_complete = true;
	cout << "Exiting loop and setting flag to " << is_complete << endl;
	break;
      }
      
      xvalue = ramp_xvalue(0.0,x,num_of_steps-1,i-1);     
      xold = xvalue;
      Hsngi.set_x(xvalue);
      aux.clear();
      contractor.optimize(exp_odd_half,*mps,aux,D);
      //Restore norm
      aux.setNormFact(1.0);
    }
  }    
  //Complete the timestep if necessary
  if(!is_complete)
  {
    Hsngi.constructUevenMPO(exp_even,dt);
    mps->clear();
    cout << "exp(-t U_even) " ;
    contractor.optimize(exp_even,aux,*mps,D);
    //Restore norm
    mps->setNormFact(1.0);
    
    Hsngi.constructUoddMPO(exp_odd_half,dt/2.0);	
    aux.clear();
    cout << "exp(-t/2 U_odd) " << endl;
    contractor.optimize(exp_odd_half,*mps,aux,D);
    //Restore norm
    aux.setNormFact(1.0);
    
    *mps = aux;
  }
  
  energy_file.close();
  
  
  
  
  E0 = curr_energy.re;
//   //String for the name of the groundstate
//   string GS_name="GS_";
//   //String needed for the conversion
//   ostringstream convert;  
//   //Build the filename
//   convert << N;
//   GS_name = GS_name+"N"+convert.str();
//   convert.str("");
//   convert.clear();
//   convert << x;
//   GS_name = GS_name+"x"+convert.str(); 
//   convert.str("");
//   convert.clear();
//   convert << mu;
//   GS_name = GS_name+"mu"+convert.str()+".dat"; 
//   
//   time(&start_local);
//   mps.exportMPS(GS_name.c_str());
//   time(&end_local);
//   writetime_gs = difftime(end_local,start_local);
//   cout << "Time for writing GS: " << writetime_gs << " s" << endl;
   
  
  mps->exportForMatlab("MPS.m");
  
  
  
//   
//   cout << "Ratio Penalty in percent/Groundstate energy: " << penalty_energy.re/E0*100.0<<" %" << endl<<endl; 
  
  //Try to calculate the charge
  const MPO& chargeop=Hsngi.getChargeMPO();  
  complex_t charge_0=contractor.contract(*mps,chargeop,*mps);
  cout << "Charge Groundstate " << charge_0.re << endl;  
  
  //Calculate the therms for the Gauss Law and see if it is fullfilled
  ofstream Lnfile;
  Lnfile.open("Ln.txt", ios_base::app);
  complex_t Lnvalue;
  MPO Ln(1);
  cout << endl;
  for(int i=1; i<L-1;i += 2){
    Hsngi.constructLopMPO(Ln,i);
    Lnvalue = contractor.contract(*mps,Ln,*mps);
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
    Gaussvalue = contractor.contract(*mps,Gauss,*mps);
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
    Zcomponent = contractor.contract(*mps,Sigmaz,*mps);
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
  cout << endl;
  Hsngi.constructCondensateMPO(Condensate);
  cfraction = contractor.contract(*mps,Condensate,*mps);
  cout << "Condensate: "<< cfraction << endl << endl;  
  
  //Ouput the results up to now (sometimes the computation of the exicted state(s) fails and the Contractor routine quits the program, therefore all results which were not stored would be lost). To avoid this, save the first part of the results now
  ofstream resfile;  
  cout << "Saving results..." << endl;
  if(file_exists(FILENAME))
  {
    time(&start_local);
    resfile.open(FILENAME,ios_base::app);  
    resfile << endl << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << x << "\t" << e << "\t" << lam << scientific <<"\t" << E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" <<penalty_energy.re << "\t";
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
  mps->clear();
  delete mps;
  
 
  
  //End Time
  time(&end);
  cout << "Run time: " << difftime(end,start) <<"s" <<  endl;  
  writetime_total = writetime_data + writetime_gauss + writetime_gs + writetime_ln + writetime_spin;
  cout << "Time for writing in total: " << writetime_total << " s" << endl;
  cout << "Time on filesystem in total (reading + writing): " << writetime_total+readtime_gs << endl; 
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

void print_MPS(MPS mps)
{
  cout << "MPS of length " << mps. getLength() << endl;
  for(int i=0; i<mps.getLength(); i++)
  {
    cout << "[" << i << "]=" <<  (mps.getA(i)).getA() << endl;
  }
}