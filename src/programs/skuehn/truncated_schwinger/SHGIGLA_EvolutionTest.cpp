/**
   \file SHGILA_EvolutionTest.cpp
   I want to test some stuff conerninc the Z_N model and not spoil my other codes
   
  \param <N> (int) Number of Fermions
  \param <d> (int) Physical dimension of Bosons in between
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <x> (double) Coupling constant
  \param <e> (double) Constant of the Operator \f$ \Theta(n)\f$
  \param <lambda> (double) Constant for the penalty-term in the Hamiltonian
   
  \author Stefan Kühn
  \date 09/09/2013

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

#include "SHGIGaussLawAngular.h"
#include "ramp_xvalue.h"

#define NUM_OF_PARAMS 10
#define FILENAME "Z2.txt"
#define RESOLUTION 10
#define EPS 1E-10
#define OFFSET 5192
//Tolerance for the error using the routine to apply an MPO and cut back the bond dimension
#define TOLERANCE 1.0E-6
//Maximum bond dimension which is allowed
#define DMAX 100
//Threshold for dectection of jumps
#define THRESHOLD 200
//Project back into Gauss Law fulfilling subspace?
#define PROJECTOR_FLAG true
//Project back into Gauss Law fulfilling subspace while sweeping right?
#define RIGHT_PROJECTOR_FLAG true
      
using namespace std;
using namespace shrt;

bool file_exists(string filename);
void print_MPS(MPS);
void print_MPO(const MPO&);
//double mean(vector<double>);
double sum(vector<double>);
double sum_abs(vector<double>);

int main(int argc,const char* argv[]){
  int cntr=0;
  
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << " N d D mu xs xe e la dt Nt ig ac O" << endl;
    cout << "N: " << "Number of spins" << endl
	 << "d: " << "Dimension of the bosonic Hilbert space" << endl
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
	 << "up:" << "Strength of the unitary penalty (optional, by default 1)" << endl;
    return -1;
  }
  
  //Number of spin sites
  int N=atoi(argv[++cntr]);
  //Total number fo bosons living on a link
  double d=atoi(argv[++cntr]);
  //double d2=atof(argv[++cntr]);
  //Bond dimension
  double D=atoi(argv[++cntr]);
  //Mass
  double mu=atof(argv[++cntr]);
  //Coupling constant
  double xs=atof(argv[++cntr]);
  double xe=atof(argv[++cntr]);
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
  //Strength of the unitary penalty
  double unitary_penalty_strength=1.0;
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
  
  //Custom accuracy for contractor, file for inital state, order for time evolution and unitary penalty strength is given
  if(argc == NUM_OF_PARAMS + 5)
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    evolution_order = atoi(argv[++cntr]);
    unitary_penalty_strength = atof(argv[++cntr]);
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
	<< "up = " << unitary_penalty_strength << endl
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
    //cout << "phys_dim[" << i << "]=" << phys_dims[i] << endl;
    if((i+1)<L){
      phys_dims[i+1] = d;
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
  
  cout << "Using Eigensolver " << contractor.getSolverName() << endl;
  
  
  //Pointer to MPS to hold the inital state which will be evolved in time
  MPS *mps;
  
  //Get the Hamiltonian object for the cQED model
  SHGIGaussLawAngular Hsngi(N,d,D,mu,xs,e,lam);  
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
  complex_t cfraction;
  //MPO to determine the condensate
  MPO Condensate(1);
  //Variable for the current x value
  double xvalue=0.0;
  //Variable for the current mass value
  double muvalue=0.0;
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
    mps = new MPS();
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
    //mps = Hsngi.constructInitalMPS(is_string,5);        
    //mps = Hsngi.constructInitalMPS(is_pzoller_string);
  } 
  
  //Determine the inital energy
  cout << "Energy before time evolution " << contractor.contract(*mps,hamil,*mps) << endl;  
  //Determine condensate before evolution
  Hsngi.constructCondensateMPO(Condensate);
  cfraction = contractor.contract(*mps,Condensate,*mps);  
  cout <<"Intial condensate: " << cfraction << " with x-value " << Hsngi.get_x() << endl;
  
  Hsngi.constructPenaltyMPO(Penalty);
  cout << "Initial Penalty energy " << contractor.contract(*mps,Penalty,*mps) << endl;
  
  //Try to find the groundstate by time evolution of an intial state   
  
  //First order Time Evolution
  if(evolution_order==1)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<        "                   %Starting first order time evolution%                    " << endl
       <<        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    
    energy_file.open("efile.txt");        
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Actual time evolution
    for(int i=1; i<=num_of_steps; ++i)
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
      }
    }
    energy_file.close();
  }
  
  
  //Second order time evolution
  if(evolution_order==2)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "      %Starting second order time evolution with unitary penalty%          " << endl
       <<        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    MPO pen_even(1),pen_odd(1);
    MPO pen_even_half(1),pen_odd_half(1);
    
    Hsngi.constructUnitarySubspaceMPO(pen_odd,pen_even,dt,unitary_penalty_strength);
    Hsngi.constructUnitarySubspaceMPO(pen_odd_half,pen_even_half,dt/2.0,unitary_penalty_strength);  
    
    energy_file.open("efile2.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\tTrunc Err\tC_0 (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Variables to monitor the intermediate charge
    //Try to calculate the charge
    const MPO& chargeoperator=Hsngi.getChargeMPO();
    complex_t curr_charge;
    
    
    //Actual time evolution
    for(int i=1; i<=num_of_steps; ++i)
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
      
      //cout << "Let's go" << endl;
      
      //First half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      truncation_errors.push_back(errN);
      //Renormalize artificially
      aux.setNormFact(1.0);
      //Complete step for even part (make sure that inital guess is empty)
      mps->clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      mps->setNormFact(1.0);
      //Final half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      aux.setNormFact(1.0);
      
      //Now the same game for the unitary penalty thing
      mps->clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(pen_odd_half,aux,*mps,D,&errN);
      truncation_errors.push_back(errN);
      //Renormalize artificially
      mps->setNormFact(1.0);
      
      aux.clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(pen_even,*mps,aux,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      aux.setNormFact(1.0);
      
      
      mps->clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(pen_odd_half,aux,*mps,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      mps->setNormFact(1.0);
      
      
      //Save result in MPS  (obsolete now)
      //*mps = aux;
      
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
	curr_charge = contractor.contract(*mps,chargeoperator,*mps);
	
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	cout << "Total truncation error " << errtot << endl;
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<< "\t" << Hsngi.get_x() << "\t" << errtot << "\t" << curr_charge << endl;
	
	if(abs(curr_overlap)< 0.99)
	{
	  //D += 20;
	  cout << "Increasing maximum bond dimension, current bond dimension is " << D << endl;
	}
      }		
    }
    energy_file.close();
  }
  
  //Still second order time evolution (but this time a try with some more optimization)
  if(evolution_order==3)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "%Starting second order time evolution with error monitoring and option to project back into GL fulfilling subspace%       " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    //MPO for projection into Gauss Law fulfilling subspace
    MPO projector(1);
    
    //MPS for the true GS (in case it is needed)
    MPS GS_exact(L,bond_dim,phys_dims);
    //Since I deform the inital state slowly, the inital state should be a good guess for the ground state
    GS_exact.setRandomState();
    //Variable to monitor the overlap of the true GS with the state computed by time evolution
    complex_t gs_overlap;
    gs_overlap.re = 0.0;
    gs_overlap.im = 0.0;
    
    energy_file.open("efile3.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\t1. Ex\tOverlap\tDbond\tTrunc Err\tC_0 (resolved every "<< RESOLUTION << " steps)" << endl;
    
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
    
    //Variables to monitor the intermediate charge
    //Try to calculate the charge
    const MPO& chargeoperator=Hsngi.getChargeMPO();
    complex_t curr_charge;
    
    //Initial half-step
    xvalue = xs;
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
    
    //Filenae for projector
    string pro_filename;
    
    //Actual time evolution
    //Note: Due to the fact that I group the odd timesteps in the loop and have do the correction at the end, this time the loop only runs to i<num_of_steps to have the desired number of timesteps at the end
    for(int i=1; i<num_of_steps; ++i)
    {
      mps_old = *mps;
      
      //Complete step for even part (make sure that inital guess is empty)
      Hsngi.constructUevenMPO(exp_even,dt);       
      mps->clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      
      truncation_errors.push_back(errN);
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
      
      truncation_errors.push_back(errN);
      //Restore norm
      aux.setNormFact(1.0);
      
      //Update Hamiltonian
      Hsngi.set_x(xvalue);
      
      //Set x old for next iteration (needed to construct Ueven)
      xold = xvalue;
      
      //Project back intor the Gauss Law fulfilling subspace if desired
      if(PROJECTOR_FLAG)
      {
	cout << "Projecting state into Gauss Law fulfilling subspace" << endl;
	for(int k=1; k<=N; k++)
	{
	  Hsngi.GaussProjector2(projector,k);
	  contractor.optimize(projector,aux,*mps,D);
	  
	  stringstream ss;//create a stringstream
	  ss << k;//add number to the stream
	  pro_filename = "gpro_";
	  pro_filename.append(ss.str());
	  pro_filename.append(".m");
	  cout << "pro_filename = " << pro_filename<<endl; 
	  projector.exportForMatlab(pro_filename.c_str());
	  
	  if(fabs(1.0-abs(contractor.contract(aux,*mps)))>1E-5)
		cout << "Overlap after projection " << contractor.contract(aux,*mps) << endl;
	  aux = *mps;		
	}
	  }
            
      //make sure that the state is normalizes again, therefore I just apply the gauge condition and normalize it with this routine (gauging is not necessary but this is the simplest way to do it)
      //mps->gaugeCond('L',true);
      //aux.setNormFact(1.0);
      //aux.gaugeCond('L',true);
      
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
	
	//cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
	
	//Save result in mps
	*mps = aux;
	
	//Ideally the truncation errors are all positive values but sometimes if the value is very small I encounter negative values. So I take the absolute values to avoid strange features
	errtot = sum_abs(truncation_errors);
	cout << scientific;
	cout << "Sum of errors " << errtot << endl;	
	cout.unsetf(ios::fixed | ios::scientific);
	truncation_errors.clear();
		
	/*cout << "TOLERANCE " << TOLERANCE << endl;
	cout << "RESOLUTION " << RESOLUTION << endl;
	cout << "Current bond " << D << endl;
	cout << "DMAX " << DMAX << endl;*/
	
	/*if((errtot >= (RESOLUTION*TOLERANCE)) && D<DMAX)
	{
	  D = (2*D > DMAX) ? DMAX : 2*D;
	  cout << "Increasing bond dimension, current bond dimension is " << D << endl;
	}*/
	
	
	//Get the values which should be monitored, first of all get a Hamiltonian MPO with the current value of x
	Hsngi.updateHMPO();
	//I think this is needless here, since I do not change the value of lambda any more
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
	curr_charge = contractor.contract(*mps,chargeoperator,*mps);
	
	/*if (abs(curr_energy)>4.0*contractor.getConvTol())
	{
	  cout << "Getting excited state" << endl;
	  excitedstates.clear();
	  excitedstates.push_back(mps);

	  contractor.findNextExcitedState(hamil,D,excitedstates,&lambdak,first_ex); 
	  cout << "Excited state energy: " << lambdak << endl<< endl;
	}*/
	
	//cout<< "Norm in step " << (i+1) << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	cout<< "x-value " << xvalue << "\tEnergy " << curr_energy << endl;
	
	/*if(xvalue >= 0.5)
	{
	  //GS_exact = *mps;
	  contractor.findGroundState(hamil,D,&E0,GS_exact);
	  gs_overlap=contractor.contract(GS_exact,*mps);
	  cout << "True GS Energy " << E0 << "\t Overlap with true GS: " << abs(gs_overlap) << endl;
	}*/
	
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << (i+1)*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<<"\t" << Hsngi.get_x() << "\t" << lambdak << "\t" << abs(gs_overlap) << "\t" << D << "\t" << errtot << "\t" << curr_charge <<endl;
	
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
	
	
	//Do a halfstep again to be in the be able to use the grouped-together time evolution again
	//Raise the value of i, since I start an an additional timestep to get the values at the right time
	i++;
	//It can be, that I now already finished all timesteps, I had to calculate, so the loop is not supposed to run any more and I don't need the following additional halfstep. For this case, exit the look and set a flag, that the correction after the loop to get a full timestep is not necessary any more
	if(i==num_of_steps)
	{
	  is_complete = true;
	  cout << "Exiting loop and setting flag to " << is_complete << endl;
	  break;
	}
	
	xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);
	xold = xvalue;
	Hsngi.set_x(xvalue);
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
  }
  
  //Again second order time evolution but this time a more elaborate in terms of timestep and stuff. The idea is to give an initial time step and a number of steps (which gives a desired total evolution time). During the evolution the time step / number of steps is adjusted accordingly to make the evolution as fast as possible but nevertheless reach the desired total time
  if(evolution_order==4)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "      %Starting second order time evolution with variable timestep%        " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    
    energy_file.open("efile4.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate (resolved every "<< RESOLUTION << " steps)" << endl;
    
    double xold,total_time,elapsed_time=0.0;
    int count = 1;
    
    
    //Calculate the desired total time
    total_time= num_of_steps*dt;
    
    //Initial half-step
    xvalue = xs;
    xold = xvalue;
    Hsngi.set_x(xvalue);
    Hsngi.constructUoddMPO(exp_odd_half,dt/2.0);
    contractor.optimize(exp_odd_half,*mps,aux,D);
    
    //Actual time evolution
    do{
      mps_old = *mps;
      
      //Complete step for even part (make sure that inital guess is empty)
      Hsngi.constructUevenMPO(exp_even,dt);       
      mps->clear();
      contractor.optimize(exp_even,aux,*mps,D);
      
      //Get next xvalue and make a complete step with the mixture of the old xvalue und the new one
      xvalue = ramp_xvalue_time(xs,xe,total_time,elapsed_time);      
      Hsngi.constructUoddMPO(exp_odd,dt,xold,xvalue);      
      aux.clear();
      contractor.optimize(exp_odd,*mps,aux,D);
      
      //Update Hamiltonian
      Hsngi.set_x(xvalue);
      
      //Set x old for next iteration (needed to construct Ueven)
      xold = xvalue;
            
      //make sure that the state is normalizes again, therefore I just apply the gauge condition and normalize it with this routine (gauging is not necessary but this is the simplest way to do it)
      //mps->gaugeCond('L',true);
      
      if((count%(RESOLUTION))==0)
      {
	//Complete the step
	Hsngi.constructUevenMPO(exp_even,dt);
	mps->clear();
	contractor.optimize(exp_even,aux,*mps,D);
	
	Hsngi.constructUoddMPO(exp_odd_half,dt/2.0);	
	aux.clear();
	contractor.optimize(exp_odd_half,*mps,aux,D);	
	
	cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
	
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
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm after time " << elapsed_time << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << elapsed_time << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x()) << endl;
	
	
	//Do a halfstep again to be in the be able to use the grouped-together time evolution again
	//Raise the value of i, since I start an an additional timestep to get the values at the right time
	elapsed_time += dt;
	xvalue = ramp_xvalue_time(xs,xe,total_time,elapsed_time);     
	xold = xvalue;
	Hsngi.set_x(xvalue);
	aux.clear();
	contractor.optimize(exp_odd_half,*mps,aux,D);
      }
      elapsed_time += dt;
      count++;
      cout << "Elapsed time " << elapsed_time << " = " << elapsed_time/dt << " of total time " << total_time <<endl;
      if(elapsed_time<total_time)
	cout << "Sollte nicht abbrechen, "<< total_time - elapsed_time <<" eps oder so "<< EPS <<endl;
      else
	cout << "Sollte abbrechen " << endl;
    
    }while(fabs(elapsed_time-total_time)>EPS);    
    //Complete the timestep
    Hsngi.constructUevenMPO(exp_even,dt);
    mps->clear();
    contractor.optimize(exp_even,aux,*mps,D);
    
    Hsngi.constructUoddMPO(exp_odd_half,dt/2.0);	
    aux.clear();
    contractor.optimize(exp_odd_half,*mps,aux,D);
    *mps = aux;
    
    energy_file.close();
  }
  
  //Still second order time evolution (but this time with exponential increase in D, currently D is doubled)
  if(evolution_order==5)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "%Starting second order time evolution with error monitoring and D-doubling%" << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    
    //MPS for the true GS (in case it is needed)
    MPS GS_exact(L,bond_dim,phys_dims);
    //Since I deform the inital state slowly, the inital state should be a good guess for the ground state
    GS_exact.setRandomState();
    //Variable to monitor the overlap of the true GS with the state computed by time evolution
    complex_t gs_overlap;
    gs_overlap.re = 0.0;
    gs_overlap.im = 0.0;
    
    energy_file.open("efile5.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\t1. Ex\tOverlap\tDbond\tTrunc Err (resolved every "<< RESOLUTION << " steps)" << endl;
    
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
    for(int i=1; i<num_of_steps; ++i)
    {
      mps_old = *mps;
      
      //Complete step for even part (make sure that inital guess is empty)
      Hsngi.constructUevenMPO(exp_even,dt);       
      mps->clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      truncation_errors.push_back(errN);
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
      truncation_errors.push_back(errN);
      //Restore norm
      aux.setNormFact(1.0);
      
      //Update Hamiltonian
      Hsngi.set_x(xvalue);
      
      //Set x old for next iteration (needed to construct Ueven)
      xold = xvalue;
            
      //make sure that the state is normalizes again, therefore I just apply the gauge condition and normalize it with this routine (gauging is not necessary but this is the simplest way to do it)
      //mps->gaugeCond('L',true);
      //aux.setNormFact(1.0);
      //aux.gaugeCond('L',true);
      
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
	
	//cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
	
	//Save result in mps
	*mps = aux;
	
	//Ideally the truncation errors are all positive values but sometimes if the value is very small I encounter negative values. So I take the absolute values to avoid strange features
	errtot = sum_abs(truncation_errors);
	cout << scientific;
	cout << "Sum of errors " << errtot << endl;	
	cout.unsetf(ios::fixed | ios::scientific);
	truncation_errors.clear();
		
	/*cout << "TOLERANCE " << TOLERANCE << endl;
	cout << "RESOLUTION " << RESOLUTION << endl;
	cout << "Current bond " << D << endl;
	cout << "DMAX " << DMAX << endl;*/
	
	if((errtot >= (RESOLUTION*TOLERANCE)) && D<DMAX)
	{
	  D = (2*D > DMAX) ? DMAX : 2*D;
	  cout << "Increasing bond dimension, current bond dimension is " << D << endl;
	}
	
	
	//Get the values which should be monitored, first of all get a Hamiltonian MPO with the current value of x
	Hsngi.updateHMPO();
	//I think this is needless here, since I do not change the value of lambda any more
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
	
	/*if (abs(curr_energy)>4.0*contractor.getConvTol())
	{
	  cout << "Getting excited state" << endl;
	  excitedstates.clear();
	  excitedstates.push_back(mps);

	  contractor.findNextExcitedState(hamil,D,excitedstates,&lambdak,first_ex); 
	  cout << "Excited state energy: " << lambdak << endl<< endl;
	}*/
	
	//cout<< "Norm in step " << (i+1) << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	cout<< "x-value " << xvalue << "\tEnergy " << curr_energy << endl;
	
	/*if(xvalue >= 0.5)
	{
	  //GS_exact = *mps;
	  contractor.findGroundState(hamil,D,&E0,GS_exact);
	  gs_overlap=contractor.contract(GS_exact,*mps);
	  cout << "True GS Energy " << E0 << "\t Overlap with true GS: " << abs(gs_overlap) << endl;
	}*/
	
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << (i+1)*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<<"\t" << Hsngi.get_x() << "\t" << lambdak << "\t" << abs(gs_overlap) << "\t" << D << "\t" << errtot <<endl;
	
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
	
	
	//Do a halfstep again to be in the be able to use the grouped-together time evolution again
	//Raise the value of i, since I start an an additional timestep to get the values at the right time
	i++;
	//It can be, that I now already finished all timesteps, I had to calculate, so the loop is not supposed to run any more and I don't need the following additional halfstep. For this case, exit the look and set a flag, that the correction after the loop to get a full timestep is not necessary any more
	if(i==num_of_steps)
	{
	  is_complete = true;
	  cout << "Exiting loop and setting flag to " << is_complete << endl;
	  break;
	}
	
	xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);
	xold = xvalue;
	Hsngi.set_x(xvalue);
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
  }
  
  //Still second order time evolution (but this time with linear increase in D)
  if(evolution_order==6)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          " %Starting second order time evolution with error monitoring and D-linear% " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    
    //MPS for the true GS (in case it is needed)
    MPS GS_exact(L,bond_dim,phys_dims);
    //Since I deform the inital state slowly, the inital state should be a good guess for the ground state
    GS_exact.setRandomState();
    //Variable to monitor the overlap of the true GS with the state computed by time evolution
    complex_t gs_overlap;
    gs_overlap.re = 0.0;
    gs_overlap.im = 0.0;
    
    energy_file.open("efile6.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\t1. Ex\tOverlap\tDbond\tTrunc Err (resolved every "<< RESOLUTION << " steps)" << endl;
    
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
    for(int i=1; i<num_of_steps; ++i)
    {
      mps_old = *mps;
      
      //Complete step for even part (make sure that inital guess is empty)
      Hsngi.constructUevenMPO(exp_even,dt);       
      mps->clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      truncation_errors.push_back(errN);
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
      truncation_errors.push_back(errN);
      //Restore norm
      aux.setNormFact(1.0);
      
      //Update Hamiltonian
      Hsngi.set_x(xvalue);
      
      //Set x old for next iteration (needed to construct Ueven)
      xold = xvalue;
            
      //make sure that the state is normalizes again, therefore I just apply the gauge condition and normalize it with this routine (gauging is not necessary but this is the simplest way to do it)
      //mps->gaugeCond('L',true);
      //aux.setNormFact(1.0);
      //aux.gaugeCond('L',true);
      
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
	
	//cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
	
	//Save result in mps
	*mps = aux;
	
	//Ideally the truncation errors are all positive values but sometimes if the value is very small I encounter negative values. So I take the absolute values to avoid strange features
	errtot = sum_abs(truncation_errors);
	cout << scientific;
	cout << "Sum of errors " << errtot << endl;	
	cout.unsetf(ios::fixed | ios::scientific);
	truncation_errors.clear();
		
	/*cout << "TOLERANCE " << TOLERANCE << endl;
	cout << "RESOLUTION " << RESOLUTION << endl;
	cout << "Current bond " << D << endl;
	cout << "DMAX " << DMAX << endl;*/
	
	if((errtot >= (RESOLUTION*TOLERANCE)) && D<DMAX)
	{
	  D = ((D+10) > DMAX) ? DMAX : D+10;
	  cout << "Increasing bond dimension, current bond dimension is " << D << endl;
	}
	
	
	//Get the values which should be monitored, first of all get a Hamiltonian MPO with the current value of x
	Hsngi.updateHMPO();
	//I think this is needless here, since I do not change the value of lambda any more
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
	
	/*if (abs(curr_energy)>4.0*contractor.getConvTol())
	{
	  cout << "Getting excited state" << endl;
	  excitedstates.clear();
	  excitedstates.push_back(mps);

	  contractor.findNextExcitedState(hamil,D,excitedstates,&lambdak,first_ex); 
	  cout << "Excited state energy: " << lambdak << endl<< endl;
	}*/
	
	//cout<< "Norm in step " << (i+1) << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	cout<< "x-value " << xvalue << "\tEnergy " << curr_energy << endl;
	
	/*if(xvalue >= 0.5)
	{
	  //GS_exact = *mps;
	  contractor.findGroundState(hamil,D,&E0,GS_exact);
	  gs_overlap=contractor.contract(GS_exact,*mps);
	  cout << "True GS Energy " << E0 << "\t Overlap with true GS: " << abs(gs_overlap) << endl;
	}*/
	
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << (i+1)*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<<"\t" << Hsngi.get_x() << "\t" << lambdak << "\t" << abs(gs_overlap) << "\t" << D << "\t" << errtot <<endl;
	
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
	
	
	//Do a halfstep again to be in the be able to use the grouped-together time evolution again
	//Raise the value of i, since I start an an additional timestep to get the values at the right time
	i++;
	//It can be, that I now already finished all timesteps, I had to calculate, so the loop is not supposed to run any more and I don't need the following additional halfstep. For this case, exit the look and set a flag, that the correction after the loop to get a full timestep is not necessary any more
	if(i==num_of_steps)
	{
	  is_complete = true;
	  cout << "Exiting loop and setting flag to " << is_complete << endl;
	  break;
	}
	
	xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);
	xold = xvalue;
	Hsngi.set_x(xvalue);
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
  }
  
  /*Second order time evolution
  While I ramp the value of x, I decrease an inital mass. The idea is to have a finite gap at the beginning which is due to the mass term. As the value for x grows I might not need the mass anymore to have a gap and therefore I can decrease it. I use second order time evolution without any optimization here because this way it is easier to ramp the values (I do not have to interpolate intermediate values). Currently I do not use the unitary penalty*/
  
  if(evolution_order==7)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "        %Starting second order time evolution with mass decrease%          " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    //MPO pen_even(1),pen_odd(1);
    MPO pen_even_half(1),pen_odd_half(1);
    
    //Hsngi.constructUnitarySubspaceMPO(pen_odd,pen_even,dt,unitary_penalty_strength);
    //Hsngi.constructUnitarySubspaceMPO(pen_odd_half,pen_even_half,dt/2.0,unitary_penalty_strength);  
    
    energy_file.open("efile7.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\tTrunc Err (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    double inital_mass = 0.5;
    ifstream mass_file;    
    if (file_exists("mass.txt"))
    {
      mass_file.open("mass.txt");
      mass_file >> inital_mass;
      cout << "Read inital value for µ: " << inital_mass << endl;
    }
    
    
    //Actual time evolution
    for(int i=1; i<=num_of_steps; ++i)
    {
      mps_old = *mps;
      
      //Get the current value of x
      xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1); 
      muvalue = ramp_muvalue(inital_mass,0.0,num_of_steps-1,i-1); 
      //cout << "x =  " << xvalue << endl;
      //cout << "mu =  " << muvalue << endl;
   
      //Get the time evolution operators and the current Hamiltonian
      Hsngi.set_x(xvalue);
      Hsngi.set_mu(muvalue);
      /*Hsngi.constructEvolutionMPO(exp_odd,exp_even,dt);
      Hsngi.constructEvolutionMPO(exp_odd_half,exp_even_half,dt/2.0);*/
      Hsngi.constructEvolutionMPOefficient(exp_odd,exp_even,dt);
      Hsngi.constructEvolutionMPOefficient(exp_odd_half,exp_even_half,dt/2.0);  
      
      
      //First half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      aux.setNormFact(1.0);
      //Complete step for even part (make sure that inital guess is empty)
      mps->clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      mps->setNormFact(1.0);
      //Final half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      aux.setNormFact(1.0);
      
      /*//Now the same game for the unitary penalty thing
      mps->clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(pen_odd_half,aux,*mps,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      mps->setNormFact(1.0);
      
      aux.clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(pen_even,*mps,aux,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      aux.setNormFact(1.0);
      
      
      mps->clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(pen_odd_half,aux,*mps,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      mps->setNormFact(1.0);*/
      
      //Save result in mps
      *mps = aux;
      
      if((i%RESOLUTION)==0)
      {
	Hsngi.updateHMPO();
	Hsngi.constructPenaltyMPO(Penalty);
	//I realized that getHMPO() returns a constant reference to the field for the Hamilton MPO of the Class CQED. As soon as I updathe the MPO field in the class the reference will always give me the updated version, so there is no need to generate a new one
	//const MPO& hamil2 = Hsngi.getHMPO(); 
	
	cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
	cout<< "µ-value for data calculation: " << Hsngi.get_mu() << endl;
	
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
	
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	cout << "Total truncation error " << errtot << endl;
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<< "\t" << Hsngi.get_x() << "\t" << errtot << endl;
      }		
    }
    energy_file.close();
  }
  
  //Still second order time evolution (this is a special version which detects jumps and has some other features which I need for debugging)
  //Currently it detects these jumps in the energy and kills the program once the jump occured
  if(evolution_order==10)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "        %Starting second order time evolution with debugging stuff%        " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    
    //MPS for the true GS (in case it is needed)
    MPS GS_exact(L,bond_dim,phys_dims);
    //Since I deform the inital state slowly, the inital state should be a good guess for the ground state
    GS_exact.setRandomState();
    //Variable to monitor the overlap of the true GS with the state computed by time evolution
    complex_t gs_overlap;
    gs_overlap.re = 0.0;
    gs_overlap.im = 0.0;
    
    complex_t Eold = contractor.contract(*mps,hamil,*mps);
    
    energy_file.open("efile10.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\t1. Ex\tOverlap\tDbond\tTrunc Err (resolved every "<< RESOLUTION << " steps)" << endl;
    
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
    for(int i=1; i<num_of_steps; ++i)
    {
      mps_old = *mps;
      
      //Complete step for even part (make sure that inital guess is empty)
      Hsngi.constructUevenMPO(exp_even,dt);       
      mps->clear();
      //cout << "exp(-t U_even) " ;
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      truncation_errors.push_back(errN);
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
      truncation_errors.push_back(errN);
      //Restore norm
      aux.setNormFact(1.0);
      
      //Update Hamiltonian
      Hsngi.set_x(xvalue);
      
      //Set x old for next iteration (needed to construct Ueven)
      xold = xvalue;
            
      //make sure that the state is normalizes again, therefore I just apply the gauge condition and normalize it with this routine (gauging is not necessary but this is the simplest way to do it)
      //mps->gaugeCond('L',true);
      //aux.setNormFact(1.0);
      //aux.gaugeCond('L',true);
      
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
	
	//cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
	
	//Save result in mps
	*mps = aux;
	
	//Ideally the truncation errors are all positive values but sometimes if the value is very small I encounter negative values. So I take the absolute values to avoid strange features
	errtot = sum_abs(truncation_errors);
	cout << scientific;
	cout << "Sum of errors " << errtot << endl;	
	cout.unsetf(ios::fixed | ios::scientific);
	truncation_errors.clear();
		
	/*cout << "TOLERANCE " << TOLERANCE << endl;
	cout << "RESOLUTION " << RESOLUTION << endl;
	cout << "Current bond " << D << endl;
	cout << "DMAX " << DMAX << endl;*/
	
	/*if((errtot >= (RESOLUTION*TOLERANCE)) && D<DMAX)
	{
	  D = (2*D > DMAX) ? DMAX : 2*D;
	  cout << "Increasing bond dimension, current bond dimension is " << D << endl;
	}*/
	
	
	//Get the values which should be monitored, first of all get a Hamiltonian MPO with the current value of x
	Hsngi.updateHMPO();
	//I think this is needless here, since I do not change the value of lambda any more
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
	
	/*if (abs(curr_energy)>4.0*contractor.getConvTol())
	{
	  cout << "Getting excited state" << endl;
	  excitedstates.clear();
	  excitedstates.push_back(mps);

	  contractor.findNextExcitedState(hamil,D,excitedstates,&lambdak,first_ex); 
	  cout << "Excited state energy: " << lambdak << endl<< endl;
	}*/
	
	//cout<< "Norm in step " << (i+1) << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	cout<< "x-value " << xvalue << "\tEnergy " << curr_energy << endl;
	
	if(fabs(curr_energy.re - Eold.re)> THRESHOLD)
	{
	  cout << "Jump detected, killing program " << endl;
	}
	
	/*if(xvalue >= 0.5)
	{
	  //GS_exact = *mps;
	  contractor.findGroundState(hamil,D,&E0,GS_exact);
	  gs_overlap=contractor.contract(GS_exact,*mps);
	  cout << "True GS Energy " << E0 << "\t Overlap with true GS: " << abs(gs_overlap) << endl;
	}*/
	
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << (i+1)*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<<"\t" << Hsngi.get_x() << "\t" << lambdak << "\t" << abs(gs_overlap) << "\t" << D << "\t" << errtot <<endl;
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Calculate the terms for the Gauss Law and see if it is fullfilled
	ofstream Lnfile;
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
	spinfile.close(); 
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	if(fabs(curr_energy.re - Eold.re)> THRESHOLD)
	{
	  cout << "Jump detected, killing program " << endl;
	  exit(666);
	}
	Eold = curr_energy;
	
	
	//Do a halfstep again to be in the be able to use the grouped-together time evolution again
	//Raise the value of i, since I start an an additional timestep to get the values at the right time
	i++;
	//It can be, that I now already finished all timesteps, I had to calculate, so the loop is not supposed to run any more and I don't need the following additional halfstep. For this case, exit the look and set a flag, that the correction after the loop to get a full timestep is not necessary any more
	if(i==num_of_steps)
	{
	  is_complete = true;
	  cout << "Exiting loop and setting flag to " << is_complete << endl;
	  break;
	}
	
	xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);
	xold = xvalue;
	Hsngi.set_x(xvalue);
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
  }
  
  //Second order time evolution, this time instead of a unitary penalty I use some projection back into the right subspace. Since the unitary penalty I used up to now is very costly, because its bond dimension is huge, this way would be cheaper and therefore favorable (in case it works). This version is a tryout of the projection thing
  if(evolution_order==12)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "      %Starting second order time evolution with unitary projector%        " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    MPO unitary_projector(1);
    
    energy_file.open("efile12.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\tTrunc Err (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Monitor number of iterations
    int iterations;
    ofstream iteration;
    iteration.open("Iterations.txt");
    
    int switch_unitary_projector_off=num_of_steps;
    
    if(file_exists("projector_off.txt"))
    {
      ifstream projectorfile;
      projectorfile.open("projector_off.txt");
      projectorfile >> switch_unitary_projector_off;
      projectorfile.close();
      cout << "# Will switch of projector after " << switch_unitary_projector_off << " steps of total " << num_of_steps << " steps" << endl;
    }
    
    //Filename for projector ofstream
    string pro_filename;
    
    //Actual time evolution
    for(int i=1; i<=num_of_steps; ++i)
    {
       mps_old = *mps;
      
      //Get the current value of x
      xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);  
      //cout << "x =  " << xvalue << endl;
   
      //Get the time evolution operators and the current Hamiltonian
      Hsngi.set_x(xvalue);
      Hsngi.constructEvolutionMPOefficient(exp_odd,exp_even,dt);
      Hsngi.constructEvolutionMPOefficient(exp_odd_half,exp_even_half,dt/2.0);  
      
      //cout << "Let's go" << endl;
      
      //First half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      aux.setNormFact(1.0);
      //Complete step for even part (make sure that inital guess is empty)
      mps->clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      mps->setNormFact(1.0);
      //Final half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      aux.setNormFact(1.0);     
      
      //Save result in MPS 
      *mps = aux;
      
      //Project back into unitary subspace
      if(i <= switch_unitary_projector_off)
      {
	for(int k=1; k<=N; k++)
	{
	  Hsngi.UnitaryPenaltyMPO(unitary_projector,k);
	  //In case I want the projectors to check them in matlab or so
	  /*stringstream ss;//create a stringstream
	  ss << k;//add number to the stream
	  pro_filename = "pro_";
	  pro_filename.append(ss.str());
	  pro_filename.append(".m");
	  cout << "pro_filename = " << pro_filename<<endl; 
	  unitary_projector.exportForMatlab(pro_filename.c_str());*/
	  aux.clear();
	  contractor.optimize(unitary_projector,*mps,aux,D);
	  *mps = aux;
	}
      }
	
      
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
	
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	cout << "Total truncation error " << errtot << endl;
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<< "\t" << Hsngi.get_x() << "\t" << errtot << endl;
	
	if(abs(curr_overlap)< 0.99)
	{
	  //D += 20;
	  cout << "Increasing maximum bond dimension, current bond dimension is " << D << endl;
	}
      }		
    }
    energy_file.close();
  }
  
/*Short explanation:
 * For debugging and stuff, I once prepared a function optimize() in contractor, which does not start sweeping from left to right, but the opposite direction. However, contractor is maintained by Mari Carmen and therefore the routine only exists in my local copy. This means that after every new version of Contractor I run into trouble, because optimize might have changed and I have to fix everything by hand. Furthermore I always have an incosistent local copy here on my machine and on the cluster. Therefore I disable this piece of code this way until I maybe have a better solution of doing this right sweeps.*/
#ifdef RIGHTSWEEP_ACTIVE
  //Still second order time evolution (since I still try to fight the jumps this version does the optimization sweeps from right to left to find out if somehow the sweep direction seems to bias the results)
  if(evolution_order==13)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          " %Starting second order time evolution with error monitoring, contractor is running from right to left this time% " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    //MPO for projection into Gauss Law fulfilling subspace
    MPO projector(1);
    
    //MPS for the true GS (in case it is needed)
    MPS GS_exact(L,bond_dim,phys_dims);
    //Since I deform the inital state slowly, the inital state should be a good guess for the ground state
    GS_exact.setRandomState();
    //Variable to monitor the overlap of the true GS with the state computed by time evolution
    complex_t gs_overlap;
    gs_overlap.re = 0.0;
    gs_overlap.im = 0.0;
    
    energy_file.open("efile13.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\t1. Ex\tOverlap\tDbond\tTrunc Err\tC_0 (resolved every "<< RESOLUTION << " steps)" << endl;
    
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
    
    //Variables to monitor the intermediate charge
    //Try to calculate the charge
    const MPO& chargeoperator=Hsngi.getChargeMPO();
    complex_t curr_charge;
    
    //Initial half-step
    xvalue = xs;
    xold = xvalue;
    Hsngi.set_x(xvalue);
    //cout << "exp(-t/2 U_odd) " ;
    Hsngi.constructUoddMPO(exp_odd_half,dt/2.0);
    //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
    errN = -1.0;
    contractor.optimize_starting_from_right(exp_odd_half,*mps,aux,D,&errN);
    
    truncation_errors.push_back(errN);
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
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize_starting_from_right(exp_even,aux,*mps,D,&errN);
      
      truncation_errors.push_back(errN);
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
      contractor.optimize_starting_from_right(exp_odd,*mps,aux,D,&errN);
      
      truncation_errors.push_back(errN);
      //Restore norm
      aux.setNormFact(1.0);
      
      //Update Hamiltonian
      Hsngi.set_x(xvalue);
      
      //Set x old for next iteration (needed to construct Ueven)
      xold = xvalue;
      
      //Project back intor the Gauss Law fulfilling subspace if desired
      if(RIGHT_PROJECTOR_FLAG)
      {
	cout << "Projecting state into Gauss Law fulfilling subspace" << endl;
	for(int k=1; k<=N; k++)
	{
	  Hsngi.GaussProjector2(projector,k);
	  contractor.optimize(projector,aux,*mps,D);
	  
	  /*stringstream ss;//create a stringstream
	  ss << k;//add number to the stream
	  pro_filename = "gpro_";
	  pro_filename.append(ss.str());
	  pro_filename.append(".m");
	  cout << "pro_filename = " << pro_filename<<endl; 
	  projector.exportForMatlab(pro_filename.c_str());*/
	  
	  if(fabs(1.0-abs(contractor.contract(aux,*mps)))>1E-5)
	    cout << "Overlap after projection " << contractor.contract(aux,*mps) << endl;
	  aux = *mps;
	}
    }
            
      //make sure that the state is normalizes again, therefore I just apply the gauge condition and normalize it with this routine (gauging is not necessary but this is the simplest way to do it)
      //mps->gaugeCond('L',true);
      //aux.setNormFact(1.0);
      //aux.gaugeCond('L',true);
      
      if(((i+1)%RESOLUTION)==0)
      {
	//Complete the step
	Hsngi.constructUevenMPO(exp_even,dt);	
	mps->clear();
	//I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
	errN = -1.0;
	contractor.optimize_starting_from_right(exp_even,aux,*mps,D,&errN);
	
	truncation_errors.push_back(errN);
	//Restore norm
	mps->setNormFact(1.0);
	
	Hsngi.constructUoddMPO(exp_odd_half,dt/2.0);	
	aux.clear();
	//I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
	errN = -1.0;
	contractor.optimize_starting_from_right(exp_odd_half,*mps,aux,D,&errN);
	
	truncation_errors.push_back(errN);
	//Restore norm
	aux.setNormFact(1.0);
	
	//cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
	
	//Save result in mps
	*mps = aux;
	
	//Ideally the truncation errors are all positive values but sometimes if the value is very small I encounter negative values. So I take the absolute values to avoid strange features
	errtot = sum_abs(truncation_errors);
	cout << scientific;
	cout << "Sum of errors " << errtot << endl;	
	cout.unsetf(ios::fixed | ios::scientific);
	truncation_errors.clear();
		
	/*cout << "TOLERANCE " << TOLERANCE << endl;
	cout << "RESOLUTION " << RESOLUTION << endl;
	cout << "Current bond " << D << endl;
	cout << "DMAX " << DMAX << endl;*/
	
	/*if((errtot >= (RESOLUTION*TOLERANCE)) && D<DMAX)
	{
	  D = (2*D > DMAX) ? DMAX : 2*D;
	  cout << "Increasing bond dimension, current bond dimension is " << D << endl;
	}*/
	
	
	//Get the values which should be monitored, first of all get a Hamiltonian MPO with the current value of x
	Hsngi.updateHMPO();
	//I think this is needless here, since I do not change the value of lambda any more
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
	curr_charge = contractor.contract(*mps,chargeoperator,*mps);
	
	/*if (abs(curr_energy)>4.0*contractor.getConvTol())
	{
	  cout << "Getting excited state" << endl;
	  excitedstates.clear();
	  excitedstates.push_back(mps);

	  contractor.findNextExcitedState(hamil,D,excitedstates,&lambdak,first_ex); 
	  cout << "Excited state energy: " << lambdak << endl<< endl;
	}*/
	
	//cout<< "Norm in step " << (i+1) << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	cout<< "x-value " << xvalue << "\tEnergy " << curr_energy << endl;
	
	/*if(xvalue >= 0.5)
	{
	  //GS_exact = *mps;
	  contractor.findGroundState(hamil,D,&E0,GS_exact);
	  gs_overlap=contractor.contract(GS_exact,*mps);
	  cout << "True GS Energy " << E0 << "\t Overlap with true GS: " << abs(gs_overlap) << endl;
	}*/
	
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << (i+1)*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<<"\t" << Hsngi.get_x() << "\t" << lambdak << "\t" << abs(gs_overlap) << "\t" << D << "\t" << errtot << "\t" << curr_charge <<endl;
	
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
	
	
	//Do a halfstep again to be in the be able to use the grouped-together time evolution again
	//Raise the value of i, since I start an an additional timestep to get the values at the right time
	i++;
	//It can be, that I now already finished all timesteps, I had to calculate, so the loop is not supposed to run any more and I don't need the following additional halfstep. For this case, exit the look and set a flag, that the correction after the loop to get a full timestep is not necessary any more
	if(i==num_of_steps)
	{
	  is_complete = true;
	  cout << "Exiting loop and setting flag to " << is_complete << endl;
	  break;
	}
	
	xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);
	xold = xvalue;
	Hsngi.set_x(xvalue);
	aux.clear();
	//I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
	errN = -1.0;
	contractor.optimize_starting_from_right(exp_odd_half,*mps,aux,D,&errN);
	
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
  }
#endif
  
  //Pure second order time evolution, which does nothing with projecting and so on and therefore is suitabel to compare to my matlab version of exact odd-even evolution
  if(evolution_order==14)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "      %Starting second order time evolution with unitary penalty%          " << endl
       <<        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    MPO pen_even(1),pen_odd(1);
    MPO pen_even_half(1),pen_odd_half(1);
    
    Hsngi.constructUnitarySubspaceMPO(pen_odd,pen_even,dt,unitary_penalty_strength);
    Hsngi.constructUnitarySubspaceMPO(pen_odd_half,pen_even_half,dt/2.0,unitary_penalty_strength);  
    
    energy_file.open("efile14.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\tTrunc Err\tC_0 (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Variables to monitor the intermediate charge
    //Try to calculate the charge
    const MPO& chargeoperator=Hsngi.getChargeMPO();
    complex_t curr_charge;
    
    //Monitor number of iterations
    int iterations;
    ofstream iteration;
    iteration.open("Iterations.txt");
    
    //Actual time evolution
    for(int i=1; i<=num_of_steps; ++i)
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
      
      //cout << "Let's go" << endl;
      
      //First half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      aux.setNormFact(1.0);
      //Complete step for even part (make sure that inital guess is empty)
      mps->clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      mps->setNormFact(1.0);
      //Final half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      aux.setNormFact(1.0);      
      
      //Save result in MPS  (obsolete now)
      *mps = aux;
      
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
	curr_charge = contractor.contract(*mps,chargeoperator,*mps);
	
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	cout << "Total truncation error " << errtot << endl;
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << i*dt << scientific << setprecision(10) << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<< "\t" << Hsngi.get_x() << "\t" << errtot << "\t" << curr_charge << endl;
      }		
    }
    energy_file.close();
  }
  
    
  //Second order time evolution, which does nothing with projecting and so on. However in this version I use a Joined Operator to have the time evolution more gauge invariant and see if it improves anything.
  if(evolution_order==15)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "                   %Starting second order time evolution %                 " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    
    energy_file.open("efile15.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\tTrunc Err\tC_0 (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Variables to monitor the intermediate charge
    //Try to calculate the charge
    const MPO& chargeoperator=Hsngi.getChargeMPO();
    complex_t curr_charge;
    
    //Variables I need to join the two operators for time evolution
    const MPO* ptrs[3];
    MPO TimeEvo(L);
    
    //Actual time evolution
    for(int i=1; i<=num_of_steps; ++i)
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
      
      ptrs[0]=&exp_odd_half;
      ptrs[1]=&exp_even;
      ptrs[2]=&exp_odd_half;
      
      MPO::join(3,ptrs,TimeEvo);
      
      //cout << "Let's go" << endl;
      
      //First half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(TimeEvo,*mps,aux,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      aux.setNormFact(1.0);
      
      //Save result in MPS  (obsolete now)
      *mps = aux;
      
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
	curr_charge = contractor.contract(*mps,chargeoperator,*mps);
	
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	cout << "Total truncation error " << errtot << endl;
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << i*dt << scientific << setprecision(10) << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<< "\t" << Hsngi.get_x() << "\t" << errtot << "\t" << curr_charge << endl;
      }		
    }
    energy_file.close();
  }
  
  //Time evolution directly via the Hamiltonian by using an approximation for the matrix exponential
  if(evolution_order==16)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "           %Starting time evolution with approximate exponential %         " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
 
    
    energy_file.open("efile16.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\tTrunc Err\tC_0 (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Variables to monitor the intermediate charge
    //Try to calculate the charge
    const MPO& chargeoperator=Hsngi.getChargeMPO();
    complex_t curr_charge;

    //Initialize aux
    aux = *mps;
    
    //aux.setRandomState();
    
    //mps->setRandomState();

    //Actual time evolution
    for(int i=1; i<=num_of_steps; i++)
    {
      mps_old = *mps;

      //Get the current value of x
      xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);  
      //xvalue = ramp_xvalue_mixed(xs,xe,num_of_steps-1,i-1,5000);  
      //cout << "x =  " << xvalue << endl;
   
      //Get the time evolution operators and the current Hamiltonian
      Hsngi.set_x(xvalue);
      Hsngi.updateHMPO();
      
      //I set the inital guess to something random because if I initialize it with *mps, I have something which is completly orthogonal to H|*mps> and therefore the approximation does not do anything
      
      //aux.setRandomState();
      aux = *mps;
      /*cout << "<aux|aux>  davor " << contractor.contract(aux,aux) << endl;
      cout << "<aux|mps>  davor " << contractor.contract(aux,*mps) << endl;*/
      contractor.approximateExponential(hamil,*mps,aux,dt,D,false,0.0,0.0);   
      //contractor.approximateExponential(hamil,*mps,aux,dt,0,true,0.0,0.0);  
      cout << endl;
      cout << "<aux|aux> " << contractor.contract(aux,aux) << endl;
      cout << "<mps|mps> " << contractor.contract(*mps,*mps) << endl;
      cout << "<aux|mps> " << abs(contractor.contract(aux,*mps)) << endl;
      aux.setNormFact(1.0);
      //Save result in MPS 
      mps->clear();
      *mps = aux;
      
      if((i%RESOLUTION)==0)
      {
	Hsngi.updateHMPO();
	//Since I might set the penalty strength in the Hamiltonian to zero but nevertheless want to monitor the penalty energy, I set the constant for the penalty MPO explicitly to one here
	Hsngi.constructPenaltyMPO(Penalty,1.0);
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
	curr_charge = contractor.contract(*mps,chargeoperator,*mps);
	
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	cout << "Total truncation error " << errtot << endl;
	cout << "Penalty Energy " << penalty_energy << endl;
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << scientific << setprecision(10) << i*dt  << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<< "\t" << Hsngi.get_x() << "\t" << errtot << "\t" << curr_charge << endl;
      }		
    }
    energy_file.close();
  }
  //First order version of the joined time evolution (although I already have a second order version). Maybe this thing is more gauge invariant since I have not done so much splitting to the operator
  if(evolution_order==17)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "          %Starting first order time evolution (joined version)%           " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    
    energy_file.open("efile17.txt");        
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Variables I need to join the two operators for time evolution
    const MPO* ptrs[2];
    MPO TimeEvo(L);
    
    //Actual time evolution
    for(int i=1; i<=num_of_steps; ++i)
    {      
      //Get the current value of x
      xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);  
      
      //Get the time evolution operators
      Hsngi.set_x(xvalue);
      //Hsngi.constructEvolutionMPO(exp_odd,exp_even,dt);
      Hsngi.constructEvolutionMPOefficient(exp_odd,exp_even,dt);     
      
      mps_old = *mps;
      
      ptrs[0]=&exp_odd;
      ptrs[1]=&exp_even;      
      MPO::join(2,ptrs,TimeEvo);
      
      // Not very optimized: apply step by step
      /*contractor.optimize(exp_odd,*mps,aux,D);
      contractor.optimize(exp_even,aux,*mps,D);*/
      
      contractor.optimize(TimeEvo,*mps,aux,D);
      
      //Save result again in *mps
      *mps=aux;
      
      //make sure that the state is normalizes again, therefore I just apply the gauge condition and normalize it with this routine (gauging is not necessary but this is the simplest way to do it)
      mps->setNormFact(1.0);
      
      
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
	//contractor.findGroundState(hamil,D,&E0,aux);
	
	Hsngi.constructCondensateMPO(Condensate);
	cfraction = contractor.contract(*mps,Condensate,*mps);
	
	
	penalty_energy = contractor.contract(*mps,Penalty,*mps);
		
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy <<"\t" << E0 << "\t" << cfraction << endl;
      }
    }
    energy_file.close();
  }
  
  //Pure second order time evolution, which does nothing with projecting and so on and therefore is suitabel to compare to my matlab version of exact odd-even evolution
  if(evolution_order==18)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "      %Starting second order time evolution with unitary penalty%          " << endl
       <<        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
    MPO pen_even(1),pen_odd(1);
    MPO pen_even_half(1),pen_odd_half(1);
    
    Hsngi.constructUnitarySubspaceMPO(pen_odd,pen_even,dt,unitary_penalty_strength);
    Hsngi.constructUnitarySubspaceMPO(pen_odd_half,pen_even_half,dt/2.0,unitary_penalty_strength);  
    
    energy_file.open("efile18.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\tTrunc Err\tC_0 (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Variables to monitor the intermediate charge
    //Try to calculate the charge
    const MPO& chargeoperator=Hsngi.getChargeMPO();
    complex_t curr_charge;
    
    MPO charge_zero_projector(1);
    Hsngi.ChargeProjectorMPO(charge_zero_projector);
    cout << "charge-zero-projector " << charge_zero_projector << endl;

    
    //Actual time evolution
    for(int i=1; i<=num_of_steps; ++i)
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
      
      //cout << "Let's go" << endl;
      
      //First half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      aux.setNormFact(1.0);
      //Complete step for even part (make sure that inital guess is empty)
      mps->clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_even,aux,*mps,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      mps->setNormFact(1.0);
      //Final half-step for odd part (make sure that inital guess is empty)
      aux.clear();
      //I have to set errN to a value different from zero, to guarantee to get the truncation error. Since the truncation error should be greater than zero, I set it to a negative value to be sure to obtain a valid truncation error
      errN = -1.0;
      contractor.optimize(exp_odd_half,*mps,aux,D,&errN);
      
      truncation_errors.push_back(errN);
      //Renormalize artificially
      aux.setNormFact(1.0);      
      
      //Project back into charge-zero-subspace
      contractor.optimize(charge_zero_projector,aux,*mps,D,&errN);
      
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
	curr_charge = contractor.contract(*mps,chargeoperator,*mps);
	
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	cout << "Total truncation error " << errtot << endl;
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << i*dt << scientific << setprecision(10) << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<< "\t" << Hsngi.get_x() << "\t" << errtot << "\t" << curr_charge << endl;
      }		
    }
    energy_file.close();
  }
  
  //Time evolution via expandind the time evolution operator as talyor series (up to first order) and then applying it to the MPS. This should be faster as approximateExponential (at least when I tweak it) and also have the same favorable property of gauge invariance. However the time evolution operator in this case is not unitary and will change my norm. Therefore it is restored by hand
  if(evolution_order==20)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "     		 %Time evolution with Taylor series%                   " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  
    energy_file.open("efile20.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\tTrunc Err\tC_0 (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Variables to monitor the intermediate charge
    //Try to calculate the charge
    const MPO& chargeoperator=Hsngi.getChargeMPO();
    complex_t curr_charge;
    
    vector <const MPS*> kets;
    vector <complex_t> coefficients;
    kets.resize(2);
    coefficients.resize(2);
    
    coefficients[0] = ONE_c;
    coefficients[1] = -1.0*dt*I_c;
    
    MPS res;
    
    res = *mps;
    aux = *mps;
    aux.setRandomState();
    
    //Actual time evolution
    for(int i=1; i<=num_of_steps; ++i)
    {
      mps_old = *mps;
      
      //Get the current value of x
      xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);  
      //cout << "x =  " << xvalue << endl;
   
      //Get the time evolution operators and the current Hamiltonian
      Hsngi.set_x(xvalue);
      Hsngi.updateHMPO();
      
      //The first order Taylor expansion of exp(-i dt H) = 1 - i*dt*H + O(dt²), therefore I prepare H|phi> and then use optimizeSum to optimize the result
      errN=-1;
      cout << "Optimizing |Phi>=H|Psi> for x="<< xvalue << endl;
      cout << "<aux|aux> = " << contractor.contract(aux,aux) << endl;
      aux.setRandomState();
      if(xvalue>0.0)
      {
	contractor.optimize(hamil,*mps,aux,D,&errN);
	//aux.setNormFact(1.0);
	aux.gaugeCond('L',true);
	truncation_errors.push_back(errN);
	kets[0] = mps;
	kets[1] = &aux;
	coefficients[0] = ONE_c;
	coefficients[1] = -1.0*dt*I_c;
      }
      else
      {
	kets[0] = mps;
	kets[1] = &aux;
	coefficients[0] = ONE_c;
	coefficients[1] = ZERO_c;
      }
	
      
      cout << "Optimizing |Psi>-i*dt|Phi>" << endl;
      contractor.optimizeSum(kets,coefficients,res,D);
      *mps = res;
      //Restore norm factor
      mps->setNormFact(1.0);
      if((i%RESOLUTION)==0)
      {
	Hsngi.updateHMPO();
	Hsngi.constructPenaltyMPO(Penalty,1.0);
	
	cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
	
	//Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
	curr_overlap = contractor.contract(mps_old,*mps);
	cout << "Overlap " << curr_overlap << "\t";
      
	//curr_energy = contractor.contract(*mps,hamil2,*mps);
	curr_energy = contractor.contract(*mps,hamil,*mps);
	curr_norm = contractor.contract(*mps,*mps);
	
	penalty_energy = contractor.contract(*mps,Penalty,*mps);
	
	Hsngi.constructCondensateMPO(Condensate);
	cfraction = contractor.contract(*mps,Condensate,*mps);
	curr_charge = contractor.contract(*mps,chargeoperator,*mps);
	
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	cout << "Total truncation error " << errtot << endl;
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << scientific << setprecision(10) << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<< "\t" << Hsngi.get_x() << "\t" << errtot << "\t" << curr_charge << endl;
      }		
    }
    energy_file.close();
  }
  //Again time evolution via taylor series, this time more optimized since I construct an operator for 1-idtH to only have to optimize once
  if(evolution_order==21)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "     		 %Time evolution with Taylor series%                   " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
 
    energy_file.open("efile21.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\tTrunc Err\tC_0 (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Variables to monitor the intermediate charge
    //Try to calculate the charge
    const MPO& chargeoperator=Hsngi.getChargeMPO();
    complex_t curr_charge;

   
    aux = *mps;
    
    MPO evolutionoperator(1);
    
    
    //Actual time evolution
    for(int i=1; i<=num_of_steps; ++i)
    {
      mps_old = *mps;
      
      //Get the current value of x
      xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);  
   
      //Get the time evolution operators and the current Hamiltonian
      Hsngi.set_x(xvalue);
      //Hsngi.updateHMPO();
      
      //The first order Taylor expansion of exp(-i dt H) = 1 - i*dt*H + O(dt²), therefore I prepare H|phi> and then use optimizeSum to optimize the result
      Hsngi.EvolutionMPO(evolutionoperator,dt);
      //cout << "EvolutionMPO " << evolutionoperator << endl;
      errN=-1.0;
      contractor.optimize(evolutionoperator,*mps,aux,D,&errN);
      aux.setNormFact(1.0);
      *mps = aux;
      //mps->setNormFact(1.0);
      
      if((i%RESOLUTION)==0)
      {
	Hsngi.updateHMPO();
	//Since I might set the penalty strength in the Hamiltonian to zero but nevertheless want to monitor the penalty energy, I set the constant for the penalty MPO explicitly to one here
	Hsngi.constructPenaltyMPO(Penalty,1.0);
	
	cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
	
	//Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
	curr_overlap = contractor.contract(mps_old,*mps);
	cout << "Overlap " << curr_overlap << "\t";
      
	//curr_energy = contractor.contract(*mps,hamil2,*mps);
	curr_energy = contractor.contract(*mps,hamil,*mps);
	curr_norm = contractor.contract(*mps,*mps);
	
	penalty_energy = contractor.contract(*mps,Penalty,*mps);
	
	Hsngi.constructCondensateMPO(Condensate);
	cfraction = contractor.contract(*mps,Condensate,*mps);
	curr_charge = contractor.contract(*mps,chargeoperator,*mps);
	
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	cout << "Total truncation error " << errtot << endl;
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << scientific << setprecision(10) << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<< "\t" << Hsngi.get_x() << "\t" << errtot << "\t" << curr_charge << endl;
      }		
    }
    energy_file.close();
  }
  
  E0 = curr_energy.re;
  //String for the name of the groundstate
  string GS_name="GS_";
  //String needed for the conversion
  ostringstream convert;  
  //Build the filename
  convert << N;
  GS_name = GS_name+"N"+convert.str();
  convert.str("");
  convert.clear();
  convert << xe;
  GS_name = GS_name+"x"+convert.str(); 
  convert.str("");
  convert.clear();
  convert << d;
  GS_name = GS_name+"d"+convert.str();
  convert.str("");
  convert.clear();
  convert << mu;
  GS_name = GS_name+"mu"+convert.str()+".dat"; 
  
  time(&start_local);
  mps->exportMPS(GS_name.c_str());
  time(&end_local);
  writetime_gs = difftime(end_local,start_local);
  cout << "Time for writing GS: " << writetime_gs << " s" << endl;
   
  
  mps->exportForMatlab("MPS.m");
  
  
  
//   
//   cout << "Ratio Penalty in percent/Groundstate energy: " << penalty_energy.re/E0*100.0<<" %" << endl<<endl; 
  
  //MPO up(1);
  //Hsngi.UnitaryPenaltyMPO(up,6);
  
  //print_MPO(up);
  
  //cout << "Unitary Penalty MPO " << up << endl;
 
  /*aux.clear();
  contractor.optimize(up,*mps,aux,D);
  *mps = aux;*/
  
  //cout << "Penalty expectation " << contractor.contract(*mps,up,*mps)<< endl;
  
  //Test the Gauss Law projector
  //cout << "Construction the Projector for the Gauss Law" << endl;
  /*MPO glp(1);
  Hsngi.GaussProjector2(glp,1);
  
  print_MPO(glp);
  
  cout << "Gauss Law Projector MPO " << glp << endl;
 
  aux.clear();
  contractor.optimize(glp,*mps,aux,D);
  *mps = aux;
  
  cout << " VOR SPINFLIP " << endl;
  for (int i =1; i<=N; i++)
  {
	Hsngi.GaussProjector2(glp,i);
	//cout << "MPO " << glp << endl;
	cout << "GLP expectation for spin site "<< i <<": " << contractor.contract(*mps,glp,*mps)<< endl;
  }
  
  complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
  mwArray sigma_plus(Indices(2,2),data_sp);
  
  complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
  mwArray sigma_minus(Indices(2,2),data_sm);
  
  //mps->applyLocalOperator(0,sigma_plus,true); 
  mps->applyLocalOperator(0,sigma_minus,true); 
  
  cout << " NACH SPINFLIP " << endl;
  for (int i =1; i<=N; i++)
  {
	Hsngi.GaussProjector2(glp,i);
	cout << "GLP expectation for spin site "<< i <<": " << contractor.contract(*mps,glp,*mps)<< endl;
  }
  
  //cout << "GLP expectation " << contractor.contract(*mps,glp,*mps)<< endl;*/
  
  //Try to calculate the charge
  const MPO& chargeop=Hsngi.getChargeMPO();  
  complex_t charge_0=contractor.contract(*mps,chargeop,*mps);
  cout << "Charge Groundstate " << charge_0.re << endl;  
  
  
  //cout << "Chareoperator " << chargeop <<  endl;
  //print_MPO(chargeop);
  
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
  ofstream resfile;  
  cout << "Saving results..." << endl;
  if(file_exists(FILENAME))
  {
    time(&start_local);
    resfile.open(FILENAME,ios_base::app);  
    resfile << endl << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << xe << "\t" << e << "\t" << lam << scientific <<"\t" << E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" <<penalty_energy.re << "\t";
  }
  else
  {
    time(&start_local);
    resfile.open(FILENAME,ios_base::app);  
    resfile << "#N" << "\t" << "d" << "\t" << "D" << "\t" << "mu" << "\t" << "x" << "\t" << "e" << "\t" << "lambda" << "\t"<< "E0" << "\t" << "C0" << "\t" << "Cfrac" << "\t" << "E1" << "\t" << "C1" << "\tE_pen"<< endl;
    resfile << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << xe << "\t" << e << "\t" << lam << scientific << "\t"<< E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" << penalty_energy.re << "\t";
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
    benchmarking << "#N\tD\td\tT_total(s)\tT_perstep(s)" << endl;
  }
  benchmarking << std::resetiosflags(std::ios::scientific) << N << "\t" << D << "\t" << d << scientific << setprecision(10) <<"\t" << runtime_total << "\t" << runtime_total/num_of_steps << endl;
  benchmarking.close();
  
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

void print_MPO(const MPO &mpo)
{
  cout << "MP0 of length " << mpo.getLength() << endl;
  for(int i=0; i<mpo.getLength(); i++)
  {
    cout << "[" << i << "]=" << (mpo.getOp(i)).getFullData() << endl;
  }
}

