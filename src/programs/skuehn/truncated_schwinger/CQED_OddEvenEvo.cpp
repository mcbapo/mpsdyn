/**
   \file OddEvenEvo.cpp
   Realisation of a simulation of the time evolution under the fundamental Hamiltonian for 1+1d cQED as proposed in Ignacios review. This file is dedicated to the time evolution via Odd-Even splitting of the time evolution operator. Basically this files contains the versions from CQED_AdiabaticEvolution.cpp in a slightly optimized and cleaned up form.
   
  \param <N> (int) Number of Fermions
  \param <d> (int) Physical dimension of Bosons in between
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <x> (double) Coupling constant
  \param <e> (double) Constant of the Operator \f$ \Theta(n)\f$
  \param <lambda> (double) Constant for the penalty-term in the Hamiltonian
   
  \author Stefan Kühn
  \date 16/12/2013

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
#include "misc.h"

#include "CQED.h"
#include "SHGIGaussLawAngular.h"
#include "ramp_xvalue.h"


#define NUM_OF_PARAMS 10
#define FILENAME "cQED.txt"
//#define FILENAME "Zn.txt"
#define RESOLUTION 10
#define BU_RESOLUTION (100*RESOLUTION)
#define EPS 1E-10
#define DMIN 3
#define THRESHOLD 0.0




using namespace std;
using namespace shrt;

//bool file_exists(string filename);
void print_MPS(MPS);
double sum(vector<double>);
double sum_abs(vector<double>);
//Function to read a file containing some parameters
int read_params(string paramfile,int &, int &, int &, double &, double &, double &, double &, double &, int &, int &, double &, int &);
//Function to save a temporary results to be able to use it again
int save_temp_result(MPS &,int, int, int, double, double, double, double, double, int, int, double, int);

int main(int argc,const char* argv[]){
  int cntr=0;
  
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << " N N0 D mu xs xe e la dt Nt ig ac O" << endl;
    cout << "N: " << "Number of spins" << endl
	 << "N0:" << "Total number of bosons/bosonic Hilbertspace dimension" << endl
	 << "D: " << "Bond dimension of the MPS" << endl
	 << "µ: " << "Dimensionless constant" << endl
	 << "xs:" << "Coupling constant to start ramping with" << endl
	 << "xe:" << "Coupling constant to end ramping with" << endl
	 << "e: " << "Constant for gauge field" << endl
	 << "la:" << "Constant for penalty term" << endl
	 << "dt:" << "Timestep size" << endl
	 << "Nt:" << "Number of timesteps which should be calculated" << endl
	 << "ig:" << "True if inital guess for the ground state should be read form HDD (optional)" << endl
	 << "ac:" << "Accuracy for contractor (optional)" << endl
	 << "O: " << "Order of the time evolution (optional, by default first order)" << endl
	 << "dx:" << "Stepsize at which intermediate values should be taken (optional)" << endl
	 << "xt:" << "Threshold above which intermediate steps are taken (optional)" << endl
	 << "ov:" << "Compute overlap with given reference data (optional)" << endl;
    return -1;
  }
  
  //Number of spin sites
  int N=atoi(argv[++cntr]);
  //Total number fo bosons living on a link
  int d=atoi(argv[++cntr]);
  //double d2=atof(argv[++cntr]);
  //Bond dimension
  int D=atoi(argv[++cntr]);
  //Mass
  double mu=atof(argv[++cntr]);
  //Coupling constant
  double xs=atof(argv[++cntr]);
  double xe=atof(argv[++cntr]);
  //Coupling constant to start with
  double xjump=xs;
  //Stepsize for steps to take an intermediate value of x
  double dx;
  bool dx_flag = false;
  //Threshold value above which intermediate x-values are evaluated
  double xt = xs;
  //Time step to start with
  int i_start = 1;
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
  //Flag indicating, if the overlap between the state and a reference state should be computed and variable for the path to the reference
  bool get_overlap=false;  
  string path_to_reference;
  MPS reference_mps;
  //Name for file for initial state
  string ig_file;
  //Variable needed to verify if file for inital guess should have dummy name or a custom name
  size_t found_dat;
  //Flag I set in case I discover a jump
  bool no_jump_flag=true;
  //Variables for time measurement
  time_t start,end;
  time_t start_local,end_local;
  time_t start_estimation,end_estimation;
  double readtime_gs=0.0,writetime_gs=0.0,writetime_gauss=0.0,writetime_spin=0.0,writetime_ln=0.0,writetime_total=0.0,writetime_data=0.0;
  //String needed for conversion
  ostringstream convert;  
  
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
  //Custom accuracy for contractor, file for inital state and order fo time evolution is given
  if(argc == NUM_OF_PARAMS + 4)
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    evolution_order = atoi(argv[++cntr]);
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file " << ig_file << endl;
  }
  //Custom accuracy for contractor, file for inital state, order for time evolution and stepsize for intermediate snapshots is given
  if(argc == NUM_OF_PARAMS + 5)
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    evolution_order = atoi(argv[++cntr]);
    dx = atof(argv[++cntr]);
    dx_flag = true;
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file " << ig_file << endl;
  }
  
  //Custom accuracy for contractor, file for inital state, order for time evolution, xt  and stepsize for intermediate snapshots is given
  if(argc == NUM_OF_PARAMS + 6)
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    evolution_order = atoi(argv[++cntr]);
    dx = atof(argv[++cntr]);
    dx_flag = true;
    xt = atof(argv[++cntr]);
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file " << ig_file << endl;
  }
  
  //Custom accuracy for contractor, file for inital state, order for time evolution, xt, get_overlap and stepsize for intermediate snapshots is given
  if(argc == NUM_OF_PARAMS + 7)
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    evolution_order = atoi(argv[++cntr]);
    dx = atof(argv[++cntr]);
    dx_flag = true;
    xt = atof(argv[++cntr]);
    get_overlap=true;
    path_to_reference = argv[++cntr];
    
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
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
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
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "dt = " << dt << endl
	<< "Nt = " << num_of_steps << endl
	<< "ig = " << ig << endl
	<< "-----------" << endl << endl;
  }
  else if (argc==NUM_OF_PARAMS+3){
     cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "N  = " << N << endl
	<< "d  = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "dt = " << dt << endl
	<< "Nt = " << num_of_steps << endl
	<< "ig = " << ig << endl
	<< "acc= " << desired_accuracy << endl
	<< "-----------" << endl << endl;
  }
  else if (argc==NUM_OF_PARAMS+4){
     cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "N  = " << N << endl
	<< "d  = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "dt = " << dt << endl
	<< "Nt = " << num_of_steps << endl
	<< "ig = " << ig << endl
	<< "acc= " << desired_accuracy << endl
	<< "dx = " << dx << endl
	<< "-----------" << endl << endl;
  }
  else if(argc==NUM_OF_PARAMS+5){
     cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "N  = " << N << endl
	<< "d  = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "dt = " << dt << endl
	<< "Nt = " << num_of_steps << endl
	<< "ig = " << ig << endl
	<< "acc= " << desired_accuracy << endl
	<< "dx = " << dx << endl
	<< "xt= " << xt << endl
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
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "dt = " << dt << endl
	<< "Nt = " << num_of_steps << endl
	<< "ig = " << ig << endl
	<< "acc= " << desired_accuracy << endl
	<< "dx = " << dx << endl
	<< "xt = " << xt << endl
	<< "pwd= " << path_to_reference << endl
	<< "-----------" << endl << endl;
  }
    
       
  //Start time measurement
  time(&start);
  
  //Total number of sites (fermionic and bosonic)
  int L=2*N-1; 
  
  //Create an array containing the physical dimensions (alternating between fermion and boson, ending with fermion)
  int *phys_dims = new int[L];
  int *bond_dim = new int[L-1];
  for(int i=0; i<L; ++i){
    phys_dims[i] = 2;
    //cout << "phys_dim[" << i << "]=" << phys_dims[i] << endl;
    if((i+1)<L){
      phys_dims[i+1] = d+1;
      //cout << "phys_dim[" << i+1 << "]=" << phys_dims[i+1] << endl;
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
  
  //CQED-Hamiltonian
  CQED Hsngi(N,d,D,mu,xs,e,lam);
  
  //Zn Hamilton
  //SHGIGaussLawAngular Hsngi(N,d,D,mu,xs,e,lam);
  
  
  //Get Hamiltonian MPO
  const MPO& hamil = Hsngi.getHMPO();  
  //MPO for the odd and even part of the evolution operator
  MPO exp_even(1),exp_odd(1);
  //Energy value
  double E0 = 0.0; 
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
  complex_t cfraction, cfraction_old;
  //MPO to determine the condensate
  MPO Condensate(1);
  //Variable for the current x value
  double xvalue=0.0;
  //File for energy values and other observables during the evolution
  ofstream energy_file;
  //Filepointer for the inital guess
  ifstream inital_guess;
  //Charge MPO and variable for the charge
  const MPO& chargeop=Hsngi.getChargeMPO();  
  complex_t charge_0;
  //Ofstream for the file to hold the final results
  ofstream resfile;  
  //Status if backupfile is correctly read
  int status=-1;
  
  //Try to open file for inital state
  inital_guess.open(ig_file.c_str());
  cout << "ig_file" << ig_file << endl;
  //In case file for inital state exists read it and adjust bond dimension, otherwise start with predefined inital state
  if(inital_guess.is_open() && ig)
  {
    //Case in which an inital guess exists and will be read
    cout << "Found initial guess for the groundstate" << endl;
    time(&start_local);
    mps = new MPS();
    mps->importMPS(ig_file.c_str());
    cout << "Bond dimension after import: "<< mps->getBond() << endl;
    cout << "Length after import: " << mps->getLength() << endl;
    status = read_params("bu_params.txt",N,d,D,mu,xs,xe,lam,dt,num_of_steps,evolution_order,xjump,i_start);
    //Check if the inital guess is a backup from an older simulation and if I have to adjust the threshold value xt at which data is taken
    if(status==0)
    {
      //Get the old value of x
      double xold = ramp_xvalue_time(xs,xe,num_of_steps*dt,((double)i_start)*dt+dt/2.0);
      //Get last value at which data was taken befor a snapshot was saved and add dx, such that one has the next value at which data is required
      //It could still happen, that I get certain datapoints twice, since the backup intervall and the intervall at which data is taken are different, so after the last backup I might already have a number of new datapoints which are before the new backup. These are then computed twice
      if(xold > xt)
	xt = xold - fmod(xold-xt,dx)+dx;
      
      cout << "Adjusting threshold value to " << xt << endl;
    }
    i_start++;
   
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
    //mps = Hsngi.constructInitalMPS(is_string,5);        
    //mps = Hsngi.constructInitalMPS(is_pzoller_string);
  } 
  
  //Determine the inital energy
  cout << "Energy before time evolution " << contractor.contract(*mps,hamil,*mps) << endl;  
  //Determine condensate before evolution
  Hsngi.constructCondensateMPO(Condensate);
  cfraction = contractor.contract(*mps,Condensate,*mps);  
  cout <<"Intial condensate: " << cfraction << " with x-value " << Hsngi.get_x() << endl;
  
  
  //First order Time Evolution
  if(evolution_order==1)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<        "                   %Starting first order time evolution%                    " << endl
       <<        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    
    energy_file.open("efile.txt");        
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Actual time evolution
    for(int i=i_start; i<=num_of_steps; ++i)
    {      
      //Get the current value of x
      xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);  
      
      //Get the time evolution operators
      Hsngi.set_x(xvalue);
      //Hsngi.constructEvolutionMPO(exp_odd,exp_even,dt);
      Hsngi.constructEvolutionMPOefficient(exp_odd,exp_even,dt);
      
      
      mps_old = *mps;
      // Not very optimized: apply step by step
      contractor.optimize(exp_odd,*mps,aux,D);
      contractor.optimize(exp_even,aux,*mps,D);
      
      //make sure that the state is normalizes again, therefore I just apply the gauge condition and normalize it with this routine (gauging is not necessary but this is the simplest way to do it)
      mps->gaugeCond('L',true);
            
      //ONLY FOR TESTING!!!
      //contractor.optimize(exp_even,mps,aux,D);
      
      if((i%RESOLUTION)==0)
      {
	//Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
	Hsngi.updateHMPO();
	Hsngi.constructPenaltyMPO(Penalty);
	//I realized that getHMPO() returns a constant reference to the field for the Hamilton MPO of the Class CQED. As soon as I updathe the MPO field in the class the reference will always give me the updated version, so there is no need to generate a new one
	//const MPO& hamil2 = Hsngi.getHMPO();  
	
	curr_overlap = contractor.contract(mps_old,*mps);
	cout << "Overlap " << curr_overlap << "\t";
      
	//curr_energy = contractor.contract(*mps,hamil2,*mps);
	curr_energy = contractor.contract(*mps,hamil,*mps);
	curr_norm = contractor.contract(*mps,*mps);
	
	aux = *mps;
	//contractor.findGroundState(hamil2,D,&E0,aux);
	contractor.findGroundState(hamil,D,&E0,aux);
	
	Hsngi.constructCondensateMPO(Condensate);
	cfraction = contractor.contract(*mps,Condensate,*mps);
	
	
	penalty_energy = contractor.contract(*mps,Penalty,*mps);
		
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy <<"\t" << E0 << "\t" << cfraction << endl;
	if(abs(curr_overlap)< 0.99)
	{
	  //D += 10;
	  cout << "Increasing maximum bond dimension, current bond dimension is " << D << endl;
	}
	//Do a backup of the current result if desired 
	if((i%BU_RESOLUTION)==0)
	{
	  save_temp_result(*mps,N,d,D,mu,xs,xe,lam,dt,num_of_steps,evolution_order,xvalue,i);
	}
      }
    }
    energy_file.close();
  }
  
  
  //Second order time evolution
  if(evolution_order==2)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<        "                  %Starting second order time evolution%                     " << endl
       <<        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    
    energy_file.open("efile2.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Actual time evolution
    for(int i=i_start; i<=num_of_steps; ++i)
    {
      mps_old = *mps;
      
      //Get the current value of x
      xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);  
      //cout << "x =  " << xvalue << endl;
   
      //Get the time evolution operators and the current Hamiltonian
      Hsngi.set_x(xvalue);
      /*Hsngi.constructEvolutionMPO(exp_odd,exp_even,dt);
      Hsngi.constructEvolutionMPO(exp_odd_half,exp_even_half,dt/2.0);*/
      Hsngi.constructEvolutionMPOefficient(exp_odd,exp_even,dt);
      Hsngi.constructEvolutionMPOefficient(exp_odd_half,exp_even_half,dt/2.0);      
      
      //First half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      contractor.optimize(exp_odd_half,*mps,aux,D);
      //Complete step for even part (make sure that inital guess is empty)
      mps->clear();
      contractor.optimize(exp_even,aux,*mps,D);
      //Final half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      contractor.optimize(exp_odd_half,*mps,aux,D);
      
      //Save result in MPS 
      *mps = aux;
      
      //make sure that the state is normalizes again, therefore I just apply the gauge condition and normalize it with this routine (gauging is not necessary but this is the simplest way to do it)
      mps->gaugeCond('L',true);
      
      if((i%RESOLUTION)==0)
      {
	Hsngi.updateHMPO();
	Hsngi.constructPenaltyMPO(Penalty);
	//I realized that getHMPO() returns a constant reference to the field for the Hamilton MPO of the Class CQED. As soon as I updathe the MPO field in the class the reference will always give me the updated version, so there is no need to generate a new one
	//const MPO& hamil2 = Hsngi.getHMPO(); 
	
	cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
	
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
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<< "\t" << Hsngi.get_x() << endl;
	
	if(abs(curr_overlap)< 0.99)
	{
	  //D += 20;
	  cout << "Increasing maximum bond dimension, current bond dimension is " << D << endl;
	}
	//Do a backup of the current result if desired 
	if((i%BU_RESOLUTION)==0)
	{
	  save_temp_result(*mps,N,d,D,mu,xs,xe,lam,dt,num_of_steps,evolution_order,xvalue,i);
	}
      }		
    }
    energy_file.close();
  }
  
  //Still second order time evolution (but this time a try with some more optimization)
  if(evolution_order==3)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "                  %Starting second order time evolution%                   " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    
    //MPS for the true GS (in case it is needed)
    MPS GS_exact(L,bond_dim,phys_dims);
    //Since I deform the inital state slowly, the inital state should be a good guess for the ground state
    GS_exact.setRandomState();
    //Variable to monitor the overlap of the true GS with the state computed by time evolution
    complex_t gs_overlap=ZERO_c;
    
    energy_file.open("efile3.txt");
    energy_file <<"#Time\t"
                <<"Energy\t"
		<<"Norm\t"
		<<"Penalty\t"
		<<"True GS\t"
		<<"Condensate\t"
		<<"x\t"
		<<"Trunc Err\t"
		<<"C_0\t"
		<<"Dbond\t"
		<<"Overlap (resolved every "<< RESOLUTION << " steps)" << endl;
    
    double xold;
    bool is_complete=false;
    
    //Stuff for excited states
    vector<MPS*> excitedstates;
    double lambdak=0.0;
    MPS first_ex(L,bond_dim,phys_dims);
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Construct the penalty MPO (with fixed strength)
    Hsngi.constructPenaltyMPO(Penalty,1.0);
    
    //Initial half-step
    xvalue =xs;
    xold = xvalue;
    Hsngi.set_x(xvalue);
    //cout << "exp(-t/2 U_odd) " ;
    Hsngi.constructUoddMPO(exp_odd_half,dt/2.0);
    //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
    errN = -1.0;
    contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
    truncation_errors.push_back(errN);
    //Restore norm
    aux.setNormFact(1.0);
    
    //Actual time evolution
    //Note: Due to the fact that I group the odd timesteps in the loop and have do the correction at the end, this time the loop only runs to i<num_of_steps to have the desired number of timesteps at the end
    for(int i=i_start; i<num_of_steps; ++i)
    {
      mps_old = *mps;
      
      //Complete step for even part (make sure that inital guess is empty)
      Hsngi.constructUevenMPO(exp_even,dt);       
      mps->clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      //Restore norm
      mps->setNormFact(1.0);
      
      //Get next xvalue and make a complete step with the mixture of the old xvalue und the new one
      xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i);
      //cout << "x = " << xvalue << ", xold = " << xold << endl;
      Hsngi.constructUoddMPO(exp_odd,dt,xold,xvalue);      
      aux.clear();
      //cout << "[exp(-t/2 U_odd) exp(-t/2 U_odd)] " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd,*mps,aux,D,&errN);
      //Restore norm
      aux.setNormFact(1.0);
      
      //Update the x-value in the Hamiltonian
      Hsngi.set_x(xvalue);      
      //Set x old for next iteration (needed to construct Ueven)
      xold = xvalue;            
      
      //Take data every RESOLUTION steps
      if(((i+1)%RESOLUTION)==0)
      {
	//Complete the step
	Hsngi.constructUevenMPO(exp_even,dt);
	mps->clear();
	//I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
	errN = -1.0;
	contractor.optimize(exp_even,aux,*mps,D,&errN);
	truncation_errors.push_back(errN);
	//Restore norm
	mps->setNormFact(1.0);
	
	Hsngi.constructUoddMPO(exp_odd_half,dt/2.0);
	aux.clear();
	//I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
	errN = -1.0;
	contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
	truncation_errors.push_back(errN);
	//Restore norm
	aux.setNormFact(1.0);
	
	//Save result in mps
	*mps = aux;
	
	//Check truncation errors
	errtot = sum_abs(truncation_errors);
	//cout.unsetf(ios::fixed | ios::scientific);
	truncation_errors.clear();
	
	//Get the values which should be monitored, first of all get a Hamiltonian MPO with the current value of x
	Hsngi.updateHMPO();
	Hsngi.constructCondensateMPO(Condensate);
	
	//Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
	curr_overlap = contractor.contract(mps_old,*mps);      
	curr_energy = contractor.contract(*mps,hamil,*mps);
	curr_norm = contractor.contract(*mps,*mps);	
	penalty_energy = contractor.contract(*mps,Penalty,*mps);	
	cfraction = contractor.contract(*mps,Condensate,*mps);
	charge_0 = contractor.contract(*mps,chargeop,*mps);
	
	/*if (abs(curr_energy)>4.0*contractor.getConvTol())
	{
	  cout << "Getting excited state" << endl;
	  excitedstates.clear();
	  excitedstates.push_back(mps);

	  contractor.findNextExcitedState(hamil,D,excitedstates,&lambdak,first_ex); 
	  cout << "Excited state energy: " << lambdak << endl<< endl;
	}*/	
	
	/*if(xvalue >= 0.5)
	{
	  //GS_exact = *mps;
	  contractor.findGroundState(hamil,D,&E0,GS_exact);
	  gs_overlap=contractor.contract(GS_exact,*mps);
	  cout << "Overlap with true GS: " << abs(gs_overlap) << endl;
	}*/	
	
	energy_file << scientific << setprecision(11)
		    << (i+1)*dt << "\t" 
	            << curr_energy.re << "\t" 
		    << curr_norm.re << "\t" 
		    << penalty_energy << "\t" 
		    << E0 <<"\t"
		    << cfraction/sqrt(Hsngi.get_x())<<"\t" 
		    << Hsngi.get_x() << "\t" 
		    << errtot << "\t" 
		    << charge_0 << "\t" 
		    << D << "\t" 
		    << abs(curr_overlap) << endl;

      
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Calculate the terms for the Gauss Law and see if it is fullfilled
	/*ofstream Lnfile;
	Lnfile.open("Ln.txt", ios_base::app);
	complex_t Lnvalue;
	MPO Ln(1);
	cout << endl;
	cout << "Getting Ln" << endl;
	for(int ind=1; ind<L-1;ind += 2){
	  Hsngi.constructLopMPO(Ln,ind);
	  Lnvalue = contractor.contract(*mps,Ln,*mps);
	  //cout << "--> Expectation value L("<< (i+1)/2 << ") = " << Lnvalue<< endl;
	  Lnfile << Lnvalue.re << "\t";
	}
	Lnfile << endl;
	Lnfile.close();
	
	//Calculate the z-Component of the spins
	ofstream spinfile;
	spinfile.open("Spin.txt", ios_base::app);
	complex_t Zcomponent;
	MPO Sigmaz(1);
	cout << endl;
	cout << "Getting Spin" << endl;
	for(int ind=0; ind<L;ind += 2){
	  Hsngi.constructsigmazMPO(Sigmaz,ind);
	  Zcomponent = contractor.contract(*mps,Sigmaz,*mps);
	  //cout << "--> Expectation value sigma_z("<< i/2+1 << ") = " << Zcomponent << endl;
	  spinfile << Zcomponent.re << "\t";
	}
	spinfile << endl;
	spinfile.close(); */	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Output
	cout << scientific;
	cout << "Step:             " << i << endl
	     << "x:                " << xvalue << endl
	     << "D:                " << D << endl
	     << "Overlap:          " << curr_overlap << endl
	     << "Norm:             " << curr_norm << endl
	     << "Energy:           " << curr_energy << endl
	     << "Truncation Error: " << errtot << endl << endl;
	
	
	//Do a halfstep again to be in the be able to use the grouped-together time evolution again
	//Raise the value of i, since I start an an additional timestep to get the values at the right time
	i++;
	//Do a backup of the current result if desired 
	if((i%BU_RESOLUTION)==0)
	{
	  save_temp_result(*mps,N,d,D,mu,xs,xe,lam,dt,num_of_steps,evolution_order,xvalue,i);
	}
	//It can be, that I now already finished all timesteps, I had to calculate, so the loop is not supposed to run any more and I don't need the following additional halfstep. For this case, exit the loop and set a flag, that the correction after the loop to get a full timestep is not necessary any more
	if(i==num_of_steps)
	{
	  is_complete = true;
	  cout << "Exiting loop and setting flag to " << is_complete << endl;
	  break;
	}
	
	xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);
	xold = xvalue;
	Hsngi.set_x(xvalue);
	Hsngi.constructUoddMPO(exp_odd_half,dt/2.0);
	aux.clear();
	//I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
	errN = -1.0;
	contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
	truncation_errors.push_back(errN);
	//Restore norm
	aux.setNormFact(1.0);
      }
    }    
    //Complete the timestep if necessary
    if(!is_complete)
    {
      Hsngi.constructUevenMPO(exp_even,dt);
      mps->clear();
      //cout << "exp(-t U_even) " ;
      contractor.optimize(exp_even,aux,*mps,D);
      //Restore norm
      mps->setNormFact(1.0);
      
      Hsngi.constructUoddMPO(exp_odd_half,dt/2.0);
      aux.clear();
      //cout << "exp(-t/2 U_odd) " << endl;
      contractor.optimize(exp_odd_half,*mps,aux,D);
      //Restore norm
      aux.setNormFact(1.0);
      
      *mps = aux;
      
      //In case I had to complete the last timestep, my value for the energy is not up to date anymore and I have to recompute it
      Hsngi.updateHMPO();      
      curr_energy = contractor.contract(*mps,hamil,*mps);
    }
    
    energy_file.close();
  } 
  
  //Fifth order Trotter decomposition with the Forest-Ruth formula
  if(evolution_order==4)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "                  %Starting fifth order time evolution%                    " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
       
    //A bunch of additional MPOs because the Forest-Ruth splitting requires four different MPOs in total
    MPO exp_odd_2(1),exp_even_2(1);
    //Constant I need for the Forest-Ruth splitting (theta=1/(2-2^(1/3))
    double theta = 1.351207191959657634047687808971460826922;
 
    //MPS for the true GS (in case it is needed)
    MPS GS_exact(L,bond_dim,phys_dims);
    //Since I deform the inital state slowly, the inital state should be a good guess for the ground state
    GS_exact.setRandomState();
    //Variable to monitor the overlap of the true GS with the state computed by time evolution
    complex_t gs_overlap=ZERO_c;
    
    energy_file.open("efile4.txt");
    energy_file <<"#Time\t"
                <<"Energy\t"
		<<"Norm\t"
		<<"Penalty\t"
		<<"True GS\t"
		<<"Condensate\t"
		<<"x\t"
		<<"Trunc Err\t"
		<<"C_0\t"
		<<"Dbond\t"
		<<"Overlap (resolved every "<< RESOLUTION << " steps)" << endl;
    
    double xold;
    bool is_complete=false;
    
    //Stuff for excited states
    vector<MPS*> excitedstates;
    double lambdak=0.0;
    MPS first_ex(L,bond_dim,phys_dims);
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Construct the penalty MPO (with fixed strength)
    Hsngi.constructPenaltyMPO(Penalty,1.0);
        
    //Actual time evolution
    for(int i=i_start; i<=num_of_steps; ++i)
    {
      mps_old = *mps;
      
      //Get the current value of x
      xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);  
      
      //Get the time evolution operators
      Hsngi.set_x(xvalue);
      
      //Complete first step
      Hsngi.constructUoddMPO(exp_odd,dt/2.0*theta);       
      aux.clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd,*mps,aux,D,&errN);
      //Restore norm
      aux.setNormFact(1.0);
      
      //Complete second step
      Hsngi.constructUevenMPO(exp_even,dt*theta);       
      mps->clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      //Restore norm
      mps->setNormFact(1.0);
      
      //Complete third step
      Hsngi.constructUoddMPO(exp_odd_2,dt/2.0*(1.0-theta));
      aux.clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_2,*mps,aux,D,&errN);
      //Restore norm
      aux.setNormFact(1.0);
      
      //Complete fourth step
      Hsngi.constructUevenMPO(exp_even_2,dt*(1-2*theta));
      mps->clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even_2,aux,*mps,D,&errN);
      //Restore norm
      mps->setNormFact(1.0);
      
      //Complete fifth step (identical to the third step)
      aux.clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_2,*mps,aux,D,&errN);
      //Restore norm
      aux.setNormFact(1.0);
      
      //Complete sixth step (identical to the second step)
      mps->clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      //Restore norm
      mps->setNormFact(1.0);
      
      //Complete seventh step (identical to the first step)
      aux.clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd,*mps,aux,D,&errN);
      //Restore norm
      aux.setNormFact(1.0);
      
      //Save result in mps
      *mps = aux;
      
      //Take data every RESOLUTION steps
      if((i%RESOLUTION)==0)
      {
	
	//Check truncation errors
	errtot = sum_abs(truncation_errors);
	//cout.unsetf(ios::fixed | ios::scientific);
	truncation_errors.clear();
	
	//Get the values which should be monitored, first of all get a Hamiltonian MPO with the current value of x
	Hsngi.updateHMPO();
	Hsngi.constructCondensateMPO(Condensate);
	
	//Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
	curr_overlap = contractor.contract(mps_old,*mps);      
	curr_energy = contractor.contract(*mps,hamil,*mps);
	curr_norm = contractor.contract(*mps,*mps);	
	penalty_energy = contractor.contract(*mps,Penalty,*mps);	
	cfraction = contractor.contract(*mps,Condensate,*mps);
	charge_0 = contractor.contract(*mps,chargeop,*mps);
	
	/*if (abs(curr_energy)>4.0*contractor.getConvTol())
	{
	  cout << "Getting excited state" << endl;
	  excitedstates.clear();
	  excitedstates.push_back(mps);

	  contractor.findNextExcitedState(hamil,D,excitedstates,&lambdak,first_ex); 
	  cout << "Excited state energy: " << lambdak << endl<< endl;
	}*/	
	
	/*if(xvalue >= 0.5)
	{
	  //GS_exact = *mps;
	  contractor.findGroundState(hamil,D,&E0,GS_exact);
	  gs_overlap=contractor.contract(GS_exact,*mps);
	  cout << "Overlap with true GS: " << abs(gs_overlap) << endl;
	}*/	
	
	energy_file << scientific << setprecision(11)
		    << i*dt << "\t" 
	            << curr_energy.re << "\t" 
		    << curr_norm.re << "\t" 
		    << penalty_energy << "\t" 
		    << E0 <<"\t"
		    << cfraction/sqrt(Hsngi.get_x())<<"\t" 
		    << Hsngi.get_x() << "\t" 
		    << errtot << "\t" 
		    << charge_0 << "\t" 
		    << D << "\t" 
		    << abs(curr_overlap) << endl;

      
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Calculate the terms for the Gauss Law and see if it is fullfilled
	/*ofstream Lnfile;
	Lnfile.open("Ln.txt", ios_base::app);
	complex_t Lnvalue;
	MPO Ln(1);
	cout << endl;
	cout << "Getting Ln" << endl;
	for(int ind=1; ind<L-1;ind += 2){
	  Hsngi.constructLopMPO(Ln,ind);
	  Lnvalue = contractor.contract(*mps,Ln,*mps);
	  //cout << "--> Expectation value L("<< (i+1)/2 << ") = " << Lnvalue<< endl;
	  Lnfile << Lnvalue.re << "\t";
	}
	Lnfile << endl;
	Lnfile.close();
	
	//Calculate the z-Component of the spins
	ofstream spinfile;
	spinfile.open("Spin.txt", ios_base::app);
	complex_t Zcomponent;
	MPO Sigmaz(1);
	cout << endl;
	cout << "Getting Spin" << endl;
	for(int ind=0; ind<L;ind += 2){
	  Hsngi.constructsigmazMPO(Sigmaz,ind);
	  Zcomponent = contractor.contract(*mps,Sigmaz,*mps);
	  //cout << "--> Expectation value sigma_z("<< i/2+1 << ") = " << Zcomponent << endl;
	  spinfile << Zcomponent.re << "\t";
	}
	spinfile << endl;
	spinfile.close(); */	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Output
	cout << scientific;
	cout << "Step:             " << i << endl
	     << "x:                " << xvalue << endl
	     << "D:                " << D << endl
	     << "Overlap:          " << abs(curr_overlap) << endl
	     << "Norm:             " << curr_norm << endl
	     << "Energy:           " << curr_energy << endl
	     << "Truncation Error: " << errtot << endl << endl;
	     
	//Do a backup of the current result if desired 
	if((i%BU_RESOLUTION)==0)
	{
	  save_temp_result(*mps,N,d,D,mu,xs,xe,lam,dt,num_of_steps,evolution_order,xvalue,i);
	}
      }
    }        
    energy_file.close();
  } 
  
  //Fifth order Trotter decomposition with the Forest-Ruth formula (with a bit of optimization, still experimental)
  if(evolution_order==5)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "                  %Starting fifth order time evolution%                    " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
       
    //A bunch of additional MPOs because the Forest-Ruth splitting requires four different MPOs in total
    MPO exp_odd_2(1),exp_even_2(1);
    //Constant I need for the Forest-Ruth splitting (theta=1/(2-2^(1/3))
    double theta = 1.351207191959657634047687808971460826922;
 
    //MPS for the true GS (in case it is needed)
    MPS GS_exact(L,bond_dim,phys_dims);
    //Since I deform the inital state slowly, the inital state should be a good guess for the ground state
    GS_exact.setRandomState();
    //Variable to monitor the overlap of the true GS with the state computed by time evolution
    complex_t gs_overlap=ZERO_c;
    
    energy_file.open("efile5.txt");
    energy_file <<"#Time\t"
                <<"Energy\t"
		<<"Norm\t"
		<<"Penalty\t"
		<<"True GS\t"
		<<"Condensate\t"
		<<"x\t"
		<<"Trunc Err\t"
		<<"C_0\t"
		<<"Dbond\t"
		<<"Overlap (resolved every "<< RESOLUTION << " steps)" << endl;
    
    double xold;
    bool is_complete=false;
    
    //Stuff for excited states
    vector<MPS*> excitedstates;
    double lambdak=0.0;
    MPS first_ex(L,bond_dim,phys_dims);
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Construct the penalty MPO (with fixed strength)
    Hsngi.constructPenaltyMPO(Penalty,1.0);
        
    //Actual time evolution
    
    //Initial half-step
    xvalue =xjump;
    xold = xvalue;
    Hsngi.set_x(xvalue);
    
    //Complete first step
    Hsngi.constructUoddMPO(exp_odd,dt/2.0*theta);       
    aux.clear();
    //cout << "exp(-t U_even) " ;
    //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
    errN = -1.0;
    contractor.optimize(exp_odd,*mps,aux,D,&errN);
    //Restore norm
    aux.setNormFact(1.0);
    
    
    //Note: Due to the fact that I group the odd timesteps in the loop and have do the correction at the end, this time the loop only runs to i<num_of_steps to have the desired number of timesteps at the end
    for(int i=i_start; i<num_of_steps; ++i)
    {
      mps_old = *mps;      
      
      
      //Complete second step
      Hsngi.constructUevenMPO(exp_even,dt*theta);       
      mps->clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      //Restore norm
      mps->setNormFact(1.0);
      
      //Complete third step
      Hsngi.constructUoddMPO(exp_odd_2,dt/2.0*(1.0-theta));
      aux.clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_2,*mps,aux,D,&errN);
      //Restore norm
      aux.setNormFact(1.0);
      
      //Complete fourth step
      Hsngi.constructUevenMPO(exp_even_2,dt*(1.0-2.0*theta));
      mps->clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even_2,aux,*mps,D,&errN);
      //Restore norm
      mps->setNormFact(1.0);
      
      //Complete fifth step (identical to the third step)
      aux.clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_2,*mps,aux,D,&errN);
      //Restore norm
      aux.setNormFact(1.0);
      
      //Complete sixth step (identical to the second step)
      mps->clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      //Restore norm
      mps->setNormFact(1.0);
      
      //Get next xvalue and make a complete step with the mixture of the old xvalue und the new one
      xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i);    
      
      //Complete seventh step combined with the first step of the next round
      Hsngi.constructUoddMPO(exp_odd,dt*theta,xold,xvalue);
      aux.clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd,*mps,aux,D,&errN);
      //Restore norm
      aux.setNormFact(1.0);
      
      //Update the x-value in the Hamiltonian
      Hsngi.set_x(xvalue);      
      //Set x old for next iteration (needed to construct Ueven)
      xold = xvalue;            
      
      //Save result in mps
      *mps = aux;
      
      //Take data every RESOLUTION steps
      if(((i+1)%RESOLUTION)==0)
      {
	//Complete the step
	
	//Complete second step
	Hsngi.constructUevenMPO(exp_even,dt*theta);       
	mps->clear();
	//cout << "exp(-t U_even) " ;
	//I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
	errN = -1.0;
	contractor.optimize(exp_even,aux,*mps,D,&errN);
	//Restore norm
	mps->setNormFact(1.0);
	
	//Complete third step
	Hsngi.constructUoddMPO(exp_odd_2,dt/2.0*(1.0-theta));
	aux.clear();
	//cout << "exp(-t U_even) " ;
	//I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
	errN = -1.0;
	contractor.optimize(exp_odd_2,*mps,aux,D,&errN);
	//Restore norm
	aux.setNormFact(1.0);
	
	//Complete fourth step
	Hsngi.constructUevenMPO(exp_even_2,dt*(1-2*theta));
	mps->clear();
	//cout << "exp(-t U_even) " ;
	//I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
	errN = -1.0;
	contractor.optimize(exp_even_2,aux,*mps,D,&errN);
	//Restore norm
	mps->setNormFact(1.0);
	
	//Complete fifth step (identical to the third step)
	aux.clear();
	//cout << "exp(-t U_even) " ;
	//I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
	errN = -1.0;
	contractor.optimize(exp_odd_2,*mps,aux,D,&errN);
	//Restore norm
	aux.setNormFact(1.0);
	
	//Complete sixth step (identical to the second step)
	mps->clear();
	//cout << "exp(-t U_even) " ;
	//I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
	errN = -1.0;
	contractor.optimize(exp_even,aux,*mps,D,&errN);
	//Restore norm
	mps->setNormFact(1.0);
	
	//Complete the last half-step
	Hsngi.constructUoddMPO(exp_odd,dt/2.0*theta);
	aux.clear();
	//cout << "exp(-t U_even) " ;
	//I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
	errN = -1.0;
	contractor.optimize(exp_odd,*mps,aux,D,&errN);
	//Restore norm
	aux.setNormFact(1.0);
	
	*mps = aux;
	
	//Check truncation errors
	errtot = sum_abs(truncation_errors);
	//cout.unsetf(ios::fixed | ios::scientific);
	truncation_errors.clear();
	
	//Get the values which should be monitored, first of all get a Hamiltonian MPO with the current value of x
	Hsngi.updateHMPO();
	Hsngi.constructCondensateMPO(Condensate);
	
	//Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
	curr_overlap = contractor.contract(mps_old,*mps);      
	curr_energy = contractor.contract(*mps,hamil,*mps);
	curr_norm = contractor.contract(*mps,*mps);	
	penalty_energy = contractor.contract(*mps,Penalty,*mps);	
	cfraction = contractor.contract(*mps,Condensate,*mps);
	charge_0 = contractor.contract(*mps,chargeop,*mps);
	
	energy_file << scientific << setprecision(11)
		    << (i+1)*dt << "\t" 
	            << curr_energy.re << "\t" 
		    << curr_norm.re << "\t" 
		    << penalty_energy << "\t" 
		    << E0 <<"\t"
		    << cfraction/sqrt(Hsngi.get_x())<<"\t" 
		    << Hsngi.get_x() << "\t" 
		    << errtot << "\t" 
		    << charge_0 << "\t" 
		    << D << "\t" 
		    << abs(curr_overlap) << endl;

	//Output
	cout << scientific;
	cout << "Step:             " << i+1 << endl
	     << "x:                " << xvalue << endl
	     << "D:                " << D << endl
	     << "Overlap:          " << abs(curr_overlap) << endl
	     << "Norm:             " << curr_norm << endl
	     << "Energy:           " << curr_energy << endl
	     << "Truncation Error: " << errtot << endl << endl;
	     
	
	//Do a halfstep again to be in the be able to use the grouped-together time evolution again
	//Raise the value of i, since I start an an additional timestep to get the values at the right time
	i++;
	//Do a backup of the current result if desired 
	if((i%BU_RESOLUTION)==0)
	{
	  save_temp_result(*mps,N,d,D,mu,xs,xe,lam,dt,num_of_steps,evolution_order,xvalue,i);
	}
	//It can be, that I now already finished all timesteps, I had to calculate, so the loop is not supposed to run any more and I don't need the following additional halfstep. For this case, exit the loop and set a flag, that the correction after the loop to get a full timestep is not necessary any more
	if(i==num_of_steps)
	{
	  is_complete = true;
	  cout << "Exiting loop and setting flag to " << is_complete << endl;
	  break;
	}
	
	xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);
	xold = xvalue;
	Hsngi.set_x(xvalue);
	Hsngi.constructUoddMPO(exp_odd,dt/2.0*theta);
	aux.clear();
	//I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
	errN = -1.0;
	contractor.optimize(exp_odd,*mps,aux,D,&errN);
	truncation_errors.push_back(errN);
	//Restore norm
	aux.setNormFact(1.0);
      }
    }        
    energy_file.close();
  } 
  
  
  /*************************************************************************************************************************************************************
   * Basically same as above but here I take the time dependency of W into account while doing the splittng, which makes the splitting formluas a bit different
   * The choice of the order is the same as above but with an additional one in fornt of it
   * **********************************************************************************************************************************************************/
  
  //Second order time evolution
  if(evolution_order==12)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<        "                  %Starting second order time evolution%                     " << endl
       <<        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    
    energy_file.open("efile12.txt");
    energy_file <<"#Time\t"
                <<"Energy\t"
		<<"Norm\t"
		<<"Penalty\t"
		<<"True GS\t"
		<<"Condensate\t"
		<<"x\t"
		<<"Trunc Err\t"
		<<"C_0\t"
		<<"Dbond\t"
		<<"Overlap (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Construct the penalty MPO (with fixed strength)
    Hsngi.constructPenaltyMPO(Penalty,1.0);
    
    //Actual time evolution
    for(int i=i_start-1; i<num_of_steps; ++i)
    {
      mps_old = *mps;
      
      //Get the current value of x
      //xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);  
      xvalue = ramp_xvalue_time(xs,xe,num_of_steps*dt,((double)i)*dt+dt/2.0);
      //cout << "x =  " << xvalue << endl;
   
      //Get the time evolution operators and the current Hamiltonian
      Hsngi.set_x(xvalue);
      Hsngi.constructEvolutionMPOefficient(exp_odd,exp_even,dt);
      Hsngi.constructEvolutionMPOefficient(exp_odd_half,exp_even_half,dt/2.0);      
      
      //First half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      truncation_errors.push_back(errN);
      //Complete step for even part (make sure that inital guess is empty)
      mps->clear();
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      truncation_errors.push_back(errN);
      //Final half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      truncation_errors.push_back(errN);
      
      //Save result in MPS 
      *mps = aux;
      
      //make sure that the state is normalizes again
      mps->setNormFact(1.0);
      
      if(((i+1)%RESOLUTION)==0)
      {
	//Get the current MPO for the Hamilton
	Hsngi.updateHMPO();
	
	//Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
	curr_overlap = contractor.contract(mps_old,*mps);      
	curr_energy = contractor.contract(*mps,hamil,*mps);
	curr_norm = contractor.contract(*mps,*mps);	
	penalty_energy = contractor.contract(*mps,Penalty,*mps);
	charge_0 = contractor.contract(*mps,chargeop,*mps);
	
	Hsngi.constructCondensateMPO(Condensate);
	cfraction = contractor.contract(*mps,Condensate,*mps);
	
	//Check truncation errors
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	//Output
	cout << scientific;
	cout << "Step:             " << i+1 << endl
	     << "x:                " << xvalue << endl
	     << "D:                " << D << endl
	     << "Overlap:          " << abs(curr_overlap) << endl
	     << "Norm:             " << curr_norm << endl
	     << "Energy:           " << curr_energy << endl
	     << "Truncation Error: " << errtot << endl << endl;
	     
	energy_file << scientific << setprecision(11)
		    << (i+1)*dt << "\t" 
	            << curr_energy.re << "\t" 
		    << curr_norm.re << "\t" 
		    << penalty_energy << "\t" 
		    << E0 <<"\t"
		    << cfraction/sqrt(Hsngi.get_x())<<"\t" 
		    << Hsngi.get_x() << "\t" 
		    << errtot << "\t" 
		    << charge_0 << "\t" 
		    << D << "\t" 
		    << abs(curr_overlap) << endl;
	
	//Do a backup of the current result if desired 
	if(((i+1)%BU_RESOLUTION)==0)
	{
	  /*I save a snapshot for i+1 because:
	   * - I already completed the time step with i, the next to proceed is i+1
	   * - after reading the input file I add +1 to i
	   * - this loop starts with i_start-1, so effectively I start with value of i, I put in the function save_temp_result
	   * - since I don't want to do timesteps twice, I add +1 to i, to have it correct again*/
	  //Todo: a bit convoluted, should be cleared up a bit
	  save_temp_result(*mps,N,d,D,mu,xs,xe,lam,dt,num_of_steps,evolution_order,xvalue,i+1);
	}
      }		
    }
    energy_file.close();
  }
  
  //Second order time evolution
  //This version has the option to specify certain values of x which are then used to take data besides the final snapshot
  if(evolution_order==13)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<        "                  %Starting second order time evolution%                     " << endl
       <<        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    
    energy_file.open("efile13.txt");
    energy_file <<"#Time\t"
                <<"Energy\t"
		<<"Norm\t"
		<<"Penalty\t"
		<<"True GS\t"
		<<"Condensate\t"
		<<"x\t"
		<<"Trunc Err\t"
		<<"C_0\t"
		<<"Dbond\t"
		<<"Overlap (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Flag if I have to take some data
    bool take_data_flag=false;
    
    //Save the inital bond dimension (which will be the maximum bond dimension reached)
    int Dmax=D;
    
    //Construct the penalty MPO (with fixed strength)
    Hsngi.constructPenaltyMPO(Penalty,1.0);
    
    time(&start_estimation);
    //Actual time evolution
    for(int i=i_start-1; i<num_of_steps; ++i)
    {
      mps_old = *mps;
      
      //Get the current value of x
      //xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);  
      xvalue = ramp_xvalue_time(xs,xe,num_of_steps*dt,((double)i)*dt+dt/2.0,3.0);
      //xvalue = ramp_xvalue_time_hybrid(xs,xe,num_of_steps*dt,((double)i)*dt+dt/2.0,150.0);
      //cout << "x =  " << xvalue << endl;
      
      //Get the current value of bon dimension
      //D = ramp_bond(DMIN,Dmax,num_of_steps-1,i);
      //cout << "Current Bond: " << D << endl;
      
      //Case that the user speciefied some value of x at which he wants to have a snapshot
      if(dx_flag && xvalue>=xt)
      {
	xt += dx;
	take_data_flag = true;
      }
   
      //Get the time evolution operators and the current Hamiltonian
      Hsngi.set_x(xvalue);
      Hsngi.constructEvolutionMPOefficient(exp_odd,exp_even,dt);
      Hsngi.constructEvolutionMPOefficient(exp_odd_half,exp_even_half,dt/2.0);      
      
      //First half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      truncation_errors.push_back(errN);
      //Complete step for even part (make sure that inital guess is empty)
      mps->clear();
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      truncation_errors.push_back(errN);
      //Final half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      truncation_errors.push_back(errN);
      
      //Save result in MPS 
      *mps = aux;
      
      //make sure that the state is normalizes again
      mps->setNormFact(1.0);
      
      if((((i+1)%RESOLUTION)==0) || take_data_flag)
      {
	//Get the current MPO for the Hamilton
	Hsngi.updateHMPO();
	
	//Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
	curr_overlap = contractor.contract(mps_old,*mps);      
	curr_energy = contractor.contract(*mps,hamil,*mps);
	curr_norm = contractor.contract(*mps,*mps);	
	penalty_energy = contractor.contract(*mps,Penalty,*mps);
	charge_0 = contractor.contract(*mps,chargeop,*mps);
	
	Hsngi.constructCondensateMPO(Condensate);
	cfraction = contractor.contract(*mps,Condensate,*mps);
	
	//Check truncation errors
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	
	
	//Output
	cout << scientific;
	cout << "Step:             " << i+1 << endl
	     << "x:                " << xvalue << endl
	     << "D:                " << D << endl
	     << "Overlap:          " << abs(curr_overlap) << endl
	     << "Norm:             " << curr_norm << endl
	     << "Energy:           " << curr_energy << endl
	     << "Truncation Error: " << errtot << endl << endl;
	     
	//Estimate remainin runtime
	time(&end_estimation);
	cout << "Remaining time (estimated): " << difftime(end_estimation,start_estimation)/(i-i_start-1)*(num_of_steps-i-i_start-1)/3600.0 << "h " << endl << endl;
	     
	energy_file << scientific << setprecision(11)
		    << (i+1)*dt << "\t" 
	            << curr_energy.re << "\t" 
		    << curr_norm.re << "\t" 
		    << penalty_energy << "\t" 
		    << E0 <<"\t"
		    << cfraction/sqrt(Hsngi.get_x())<<"\t" 
		    << Hsngi.get_x() << "\t" 
		    << errtot << "\t" 
		    << charge_0 << "\t" 
		    << D << "\t" 
		    << abs(curr_overlap) << endl;
	
	//Do a backup of the current result if desired 
	if(((i+1)%BU_RESOLUTION)==0)
	{
	  /*I save a snapshot for i+1 because:
	   * - I already completed the time step with i, the next to proceed is i+1
	   * - after reading the input file I add +1 to i
	   * - this loop starts with i_start-1, so effectively I start with value of i, I put in the function save_temp_result
	   * - since I don't want to do timesteps twice, I add +1 to i, to have it correct again*/
	  //Todo: a bit convoluted, should be cleared up a bit	  
	  save_temp_result(*mps,N,d,D,mu,xs,xe,lam,dt,num_of_steps,evolution_order,xvalue,i+1);
	}
	//Case that I reached some value of x, where I want to have a snapshot of data
	if(take_data_flag)
	{
	  curr_overlap = -1.2*ONE_c;
	  //Case that an overlap with a reference state should be computed
	  if(get_overlap)
	  {
	    //String for the name of the groundstate
	    string reference_name = "GS_";
	  
	    //Build the filename
	    convert.str("");
	    convert.clear();
	    convert << N;
	    reference_name = reference_name+"N"+convert.str();
	    convert.str("");
	    convert.clear();
	    convert << d;
	    reference_name = reference_name+"N0"+convert.str(); 
	    convert.str("");
	    convert.clear();
	    convert <<xt-dx;
	    reference_name = reference_name+"x"+convert.str(); 	    
	    convert.str("");
	    convert.clear();
	    convert << mu;	    
	    reference_name =path_to_reference+reference_name+"mu"+convert.str()+".dat"; 
	    
	    cout << "Computing overlap with " << reference_name << endl;
	    
	    reference_mps.importMPS(reference_name.c_str());
	    
	    curr_overlap = contractor.contract(*mps,reference_mps);
	  }
	  
	  
	  //Take the data
	  if(file_exists(FILENAME))
	  {
	    resfile.open(FILENAME,ios_base::app);  
	    resfile << endl << scientific << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << Hsngi.get_x() << "\t" << e << "\t" << lam << scientific <<"\t" << curr_energy.re << "\t" << charge_0.re << "\t" << cfraction.re << "\t" <<penalty_energy.re << "\t" << curr_overlap << "\t";
	  }
	  else
	  {
	    resfile.open(FILENAME,ios_base::app);  
	    resfile << "#N" << "\t" << "d" << "\t" << "D" << "\t" << "mu" << "\t" << "x" << "\t" << "e" << "\t" << "lambda" << "\t"<< "E0" << "\t" << "C0" << "\t" << "Cfrac" << "\t" << "E1" << "\t" << "C1" << "\tE_pen"<< endl;
	    resfile << scientific << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << Hsngi.get_x() << "\t" << e << "\t" << lam << scientific << "\t"<< curr_energy.re << "\t" << charge_0.re << "\t" << cfraction.re << "\t" << penalty_energy.re << "\t" << curr_overlap << "\t";  
	    
	  }
	  resfile.close();
	  
	  //Reset the flag
	  take_data_flag=false;
	  
	  //Save snapshot
	  //String for the name of the groundstate
	  string GS_name="GS_";
	  
	  //Build the filename
	  convert.str("");
	  convert.clear();
	  convert << N;
	  GS_name = GS_name+"N"+convert.str();
	  convert.str("");
	  convert.clear();
	  convert <<xt-dx;
	  GS_name = GS_name+"x"+convert.str(); 
	  convert.str("");
	  convert.clear();
	  convert << d;
	  GS_name = GS_name+"N0"+convert.str(); 
	  convert.str("");
	  convert.clear();
	  convert << mu;
	  GS_name = GS_name+"mu"+convert.str()+".dat"; 
	  //Save the ground state
	  time(&start_local);
	  mps->exportMPS(GS_name.c_str());
	  time(&end_local);
	  writetime_gs = difftime(end_local,start_local);
	  cout << "Time for writing GS: " << writetime_gs << " s" << endl;
	  
	  //I kill the program now, in case my overlap is too small, to have at least the datapoint, at which I killed the program
	  if(get_overlap && (abs(curr_overlap)<THRESHOLD))
	  {
	    cout << "Overlap between current state and reference state is " << abs(curr_overlap) << " and therefore smaller than " << THRESHOLD << ", program will be aborted" << endl;
	    exit(666);
	  }
	  
	  cout << "Succesfully took data for x=" << Hsngi.get_x() << endl << endl;
	}
      }		
    }
    energy_file.close();
  }
  
  //In previous versions of the programm it might happen, that I have already completed the simulation and then the programm crashes (while the reference state for the overlap does not exist). In this case the loop won't run anymore, therefore I have to compute stuff by hand
  if (i_start > num_of_steps)
  {
    cout << endl << "WARNING, simulation is already terminated, computing last values directly" << endl << endl;
    xvalue = ramp_xvalue_time(xs,xe,num_of_steps*dt,((double)(num_of_steps-1))*dt+dt/2.0,3.0);
    
    //Set the x-value
    Hsngi.set_x(xvalue);    
    //Get the current MPO for the Hamilton
    Hsngi.updateHMPO();
      
    //Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
    curr_energy = contractor.contract(*mps,hamil,*mps);
    curr_norm = contractor.contract(*mps,*mps);	
    penalty_energy = contractor.contract(*mps,Penalty,*mps);
    charge_0 = contractor.contract(*mps,chargeop,*mps);
    
    Hsngi.constructCondensateMPO(Condensate);
    cfraction = contractor.contract(*mps,Condensate,*mps);
    
  }
  
  E0 = curr_energy.re;
  
  //String for the name of the groundstate
  string GS_name="GS_";
  //Build the filename
  convert.str("");
  convert.clear();
  convert << N;
  GS_name = GS_name+"N"+convert.str();
  convert.str("");
  convert.clear();
  convert << xe;
  GS_name = GS_name+"x"+convert.str(); 
  convert.str("");
  convert.clear();
  convert << d;
  GS_name = GS_name+"N0"+convert.str(); 
  convert.str("");
  convert.clear();
  convert << mu;
  GS_name = GS_name+"mu"+convert.str()+".dat"; 
  //Save the ground state
  time(&start_local);
  mps->exportMPS(GS_name.c_str());
  time(&end_local);
  writetime_gs = difftime(end_local,start_local);
  cout << "Time for writing GS: " << writetime_gs << " s" << endl;
   
  
  mps->exportForMatlab("MPS.m");
  cout << "Ratio Penalty in percent/Groundstate energy: " << penalty_energy.re/E0*100.0<<" %" << endl<<endl; 
  
  //Try to calculate the charge
  charge_0=contractor.contract(*mps,chargeop,*mps);
  cout << "Charge Groundstate " << charge_0.re << endl;  
  
  //Calculate the terms for the Gauss Law and see if it is fullfilled
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
  cout << "Saving results..." << endl;
  if(file_exists(FILENAME))
  {
    time(&start_local);
    resfile.open(FILENAME,ios_base::app);  
    resfile << endl << scientific << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << xe << "\t" << e << "\t" << lam << scientific <<"\t" << E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" <<penalty_energy.re << "\t";
  }
  else
  {
    time(&start_local);
    resfile.open(FILENAME,ios_base::app);  
    resfile << "#N" << "\t" << "d" << "\t" << "D" << "\t" << "mu" << "\t" << "x" << "\t" << "e" << "\t" << "lambda" << "\t"<< "E0" << "\t" << "C0" << "\t" << "Cfrac" << "\t" << "E1" << "\t" << "C1" << "\tE_pen"<< endl;
    resfile << scientific << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << xe << "\t" << e << "\t" << lam << scientific << "\t"<< E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" << penalty_energy.re << "\t" ;
  }
  resfile.close();
  
  //Get the final overlap in case it is desired
  curr_overlap = -1.2*ONE_c;
  //Case that an overlap with a reference state should be computed
  if(get_overlap)
  {
    //String for the name of the groundstate
    string reference_name = "GS_";
  
    //Build the filename
    convert.str("");
    convert.clear();
    convert << N;
    reference_name = reference_name+"N"+convert.str();
    convert.str("");
    convert.clear();
    convert << d;
    reference_name = reference_name+"N0"+convert.str(); 
    convert.str("");
    convert.clear();
    //convert << (int) Hsngi.get_x();
    //convert << setprecision(1) << Hsngi.get_x();
    convert << xe;
    reference_name = reference_name+"x"+convert.str(); 	    
    convert.str("");
    convert.clear();
    convert << mu;	    
    reference_name =path_to_reference+reference_name+"mu"+convert.str()+".dat"; 
    
    cout << "Computing overlap with " << reference_name << endl;
    
    reference_mps.importMPS(reference_name.c_str());
    
    curr_overlap = contractor.contract(*mps,reference_mps);
  }
  resfile.open(FILENAME,ios_base::app); 
  resfile << curr_overlap << "\t";
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
  double runtime_total = difftime(end,start);
  cout << "Run time: " << runtime_total <<"s" <<  endl;  
  writetime_total = writetime_data + writetime_gauss + writetime_gs + writetime_ln + writetime_spin;
  cout << "Time for writing in total: " << writetime_total << " s" << endl;
  cout << "Time on filesystem in total (reading + writing): " << writetime_total+readtime_gs << endl; 
  
  //Simply uncomment in case benchmarking file should be written
  ofstream benchmarking;
  if(file_exists("Benchmark.txt"))
  {
    benchmarking.open("Benchmark.txt",ios_base::app);
  }
  else
  {
    benchmarking.open("Benchmark.txt");
    benchmarking << "#N\tD\tN_0\tT_total(s)\tT_perstep(s)" << endl;
  }
  benchmarking << std::resetiosflags(std::ios::scientific) << N << "\t" << D << "\t" << d << scientific << setprecision(10) <<"\t" << runtime_total << "\t" << runtime_total/num_of_steps << endl;
  benchmarking.close();
  
  
  cout << "End of program" << endl;
  return 0;
}

/*bool file_exists(string filename)
{
  ifstream file;
  file.open(filename.c_str());
  if(file.good())
    return true;
  else
    return false;
}*/

void print_MPS(MPS mps)
{
  cout << "MPS of length " << mps. getLength() << endl;
  for(int i=0; i<mps.getLength(); i++)
  {
    cout << "[" << i << "]=" <<  (mps.getA(i)).getA() << endl;
  }
}

double sum(vector<double> numbers)
{
  double res=0.0;
  for (int i=0; i<numbers.size(); i++)
    res += numbers[i];
  return res;
}

double sum_abs(vector<double> numbers)
{
  double res=0.0;
  for (int i=0; i<numbers.size(); i++)
    res += fabs(numbers[i]);
  return res;
}

int read_params(string paramfile,int &N, int &d, int &D, double &mu, double &xs, double &xe, double &lam, double &dt, int &steps, int &evolution_order, double &xjump, int &ijump)
{
  ifstream pfile;
  string dummy;
  
  cout << "Trying to read file " << paramfile << endl;
  if(!file_exists(paramfile))
    return -1;
  else
  {
    pfile.open(paramfile.c_str());    
    pfile >> dummy; pfile >> dummy; pfile >> N;
    pfile >> dummy; pfile >> dummy; pfile >> d;
    pfile >> dummy; pfile >> dummy; pfile >> D;
    pfile >> dummy; pfile >> dummy; pfile >> mu;
    pfile >> dummy; pfile >> dummy; pfile >> xs;
    pfile >> dummy; pfile >> dummy; pfile >> xe;
    pfile >> dummy; pfile >> dummy; pfile >> lam;
    pfile >> dummy; pfile >> dummy; pfile >> dt;
    pfile >> dummy; pfile >> dummy; pfile >> steps;
    pfile >> dummy; pfile >> dummy; pfile >> evolution_order;
    pfile >> dummy; pfile >> dummy; pfile >> xjump;
    pfile >> dummy; pfile >> dummy; pfile >> ijump;
    
    cout << "Read N         = " << N << endl;
    cout << "Read N_0       = " << d<< endl;
    cout << "Read D         = " << D<< endl;
    cout << "Read mu        = " << mu<< endl;
    cout << "Read xs        = " << xs<< endl;
    cout << "Read xe        = " << xe<< endl;
    cout << "Read lamda     = " << lam<< endl;
    cout << "Read dt        = " << dt<< endl;
    cout << "Read steps     = " << steps<< endl;
    cout << "Read evo_order = " << evolution_order<< endl;
    cout << "Read xjump     = " << xjump<< endl;
    cout << "Read ijump     = " << ijump<< endl;
    
    pfile.close();
    return 0;
  }
}

int save_temp_result(MPS &mps,int N, int d, int D, double mu, double xs, double xe, double lam, double dt, int steps, int evolution_order, double xjump, int ijump)
{
  //Output file
  ofstream paramfile;
  //Variables for timestamp
  time_t timer;
  struct tm * timeinfo;
  
  time(&timer);
  timeinfo = localtime (&timer);
  
  cout << "Trying to write parameterfile bu_params.txt" << endl;
 
  paramfile.open("bu_params.txt");    
  paramfile << "N = " << N << endl;
  paramfile << "d = " << d << endl;
  paramfile << "D = " << D << endl;
  paramfile << "mu = " << mu << endl;
  paramfile << "xs = " << xs << endl;
  paramfile << "xe = " << xe << endl;
  paramfile << "lambda = " << lam << endl;
  paramfile << "dt = " << dt << endl;
  paramfile << "steps = " << steps << endl;
  paramfile << "evo_order = " << evolution_order << endl;
  paramfile << scientific << setprecision(11) ;
  paramfile << "x_save = " << xjump << endl;
  paramfile << "i_save = " << ijump << endl;
  paramfile << "Time: " << asctime(timeinfo) << endl;
  paramfile.close();
  
  mps.exportMPS("MPS.dat");
  
  cout << "Successfully saved a snapshot" << endl << endl;
  
  paramfile.close();
  
  return 0;
}

