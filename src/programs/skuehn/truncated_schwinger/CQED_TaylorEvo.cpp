/**
   \file CQED_TaylorEvo.cpp
   Realisation of a simulation of the time evolution under the fundamental Hamiltonian for 1+1d cQED as proposed in Ignacios review. This file specially contains the evolution via approximation of the evolution operator by a Taylor series up to first order. In comparison to CQED_AdiabaticEvolution.cpp this file is a bit cleaner and offers more options which might be useful especially on the cluster. 
   
  \param <N> (int) Number of spins
  \param <N0> (int) Total number of bosons an a link
  \param <D> (int) Bond dimension of the MPS
  \param <µ> (double) Constant for the mass term
  \param <xs> (double) Coupling constant to start ramping with" << endl
  \param <xe> (double) Coupling constant to end ramping with" << endl
  \param <e> (double) Constant for gauge field (deprecated, should be set to 1)
  \param <la> (double) Constant for penalty term in the Hamiltonian (does not affect the time evolution, only useful for energy monitoring)
  \param <dt> (double) Timestep size
  \param <Nt> (int) Number of timesteps which should be calculated
  \param <ig> (void) If an arbitrary argument ist given here, an inital guess for the ground state form HDD will be read. If the argmument ist a valid ground state filename, the corresponding file will be read (if existing), otherwise the standard file "GS_guess.dat" will be read. If there is no ground state file found, the program proceeds as there was no ground state given at all (optional)
  \param <ac> (double) Accuracy for contractor (optional)
  \param <O> (int) Select the type of evolution (optional)
  \param <ip> (string) Select the parameterfile corresponding to the inital guess (optional)
   
  \author Stefan Kühn
  \date 22/11/2013

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

#define NUM_OF_PARAMS 10
#define FILENAME "cQED.txt"
#define RESOLUTION 10
//Tolerance for the error using the routine to apply an MPO and cut back the bond dimension
#define TOLERANCE 1.0E-6
//On the cluster it might happen that I run out of walltime, and therefore would loose all of my results. If the SAFTEYMODE is turned on, I write my results from time to time and therefore I am in principle able to read them again and continue in a new run.
#define SAFTEYMODE false
//Resolution at which I want to save my GS in case I turned SAFTEYMODE on
#define INTERVAL 100

//Norm output for debugging
#define LOWER_TIME -0.01
#define UPPER_TIME 25.0


using namespace std;
using namespace shrt;

//Routine to check if a file with given filename exists
bool file_exists(string filename);
//Print the entire MPS (which means the content of every matrix!)
void print_MPS(MPS);
//Print the entire MPO (which means the content of every matrix!)
void print_MPO(MPO&);
//Sum up the entries of a vector of doubles
double sum(vector<double>);
//Sum up the absolute values of the entries of a vector of doubles
double sum_abs(vector<double>);
//Read parameter from file (in case I want to resume an old evolution
int read_params(string paramfile,int &, int &, int &, double &, double &, double &, double &, double &, int &, int &, double &, int &);

int main(int argc,const char* argv[]){
  int cntr=0;
  
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << " N N0 D mu xs xe e la dt Nt ig ac O" << endl;
    cout << "N: " << "Number of spins" << endl
	 << "N0:" << "Total number of bosons" << endl
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
	 << "O: " << "Type of evolution:" << endl
	 << "   1) Linear ramp, no projection" << endl
	 << "   2) Cubic ramp, no projection" << endl
	 << "   3) Linear ramp, projection" << endl
	 << "   4) Cubic ramp, projection " << endl
	 << "ip:" << "File for initial parameters (optional)" << endl;
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
  //Time step to start from (by default 1, but in case I want to read a snapshot I have to adjust this)
  int i_start = 1;
  //x-value for the input file
  double xsnap;
  //Inital state which should be used?
  bool ig = false;
  //Set the accuracy of contractor?
  bool set_accuracy = false;
  //Variable for the desired accuracy in case it should be set
  double desired_accuracy;
  //Should a parameterfile be read?
  bool ip = false;
  //Name for file for initial state
  string ig_file;
  //Variable needed to verify if file for inital guess should have dummy name or a custom name
  size_t found_dat;
  //Name for the parameter file
  string ip_file;
  //Flag I set in case I discover a jump
  bool no_jump_flag=true;
  //Vector to save the errors during truncation
  vector<double> truncation_errors;
  //Variable to hold the error for a single truncation
  double errN;
  //Variable to hold the summed error of all truncations since the last reset
  double errtot;
  //Variables for time measurement
  time_t start,end;
  time_t start_local,end_local;
  double readtime_gs=0.0,writetime_gs=0.0,writetime_gauss=0.0,writetime_spin=0.0,writetime_ln=0.0,writetime_total=0.0,writetime_data=0.0;
  //Energy value
  double E0 = 0.0; 
  //Variable for the current energy value during time evolution
  complex_t curr_energy;
  //Variable for the current norm of the state during time evolution
  complex_t curr_norm;
  //Variable for the current overlap with the wavefunction of the previous step
  complex_t curr_overlap;
  //Variable for the value of the penalty energy
  complex_t penalty_energy;
  //Variable for the condensate
  complex_t cfraction, cfraction_old;
  //Variable to monitor the intermediate charge
  complex_t curr_charge;
  //Variable for the current x value during the evolution
  double xvalue=0.0;
  //File for energy values and other observables during the evolution
  ofstream energy_file;
  //Filepointer for the inital guess
  ifstream inital_guess;
  
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
  //Custom accuracy for contractor, file for inital state and order of time evolution is given
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
  //Custom accuracy for contractor, file for inital state,order of time evolution and parameterfile is given
  if(argc == NUM_OF_PARAMS + 5)
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
    ip = true;
    ip_file = argv[++cntr];
    
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
      << "O  = " << evolution_order << endl
      << "-----------" << endl << endl;
  }
  else {
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
	<< "O  = " << evolution_order << endl
	<< "ip = " << ip_file << endl
	<< "-----------" << endl << endl;
  }
  
    
  //Start the time measurement
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
  //Auxiliary MPS to store an old version before a timestep to be able to monitor changes
  MPS aux,mps_old;
  //Get the Hamiltonian object for the cQED model
  CQED Hsngi(N,d,D,mu,xs,e,lam);  
  //Get Hamiltonian MPO
  const MPO& hamil = Hsngi.getHMPO();  
  //MPO for the the evolution operator
  MPO evolutionoperator(1);
  //MPO for the penalty term
  MPO Penalty(1);
  Hsngi.constructPenaltyMPO(Penalty,1.0);
  //MPO to determine the condensate
  MPO Condensate(1);
  //MPO for the chargeoperator
  const MPO& chargeoperator=Hsngi.getChargeMPO();
  //MPO for the projector
  MPO projector(1);
  //In case I want to project back, I have to initialize the projector
  if(evolution_order==3 || evolution_order==4)
  {
    Hsngi.GaussProjector2(projector);
  }
  
  //Bond dimension to start und to end with
  /*int D_start = 5;
  int D_end = D;
  D=D_start;*/
  
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
      mps->increaseBondDimension(D);
      //Increase bond dimension again, in case the inital guess has a smaller bond dimension
      //This time with a little noise, maybe better to avoid stability issues?
      //mps->increaseBondDimensionWithNoise(D,1E-2);     
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
    
    //See if an file with parameters is given
    if(ip && file_exists(ip_file))
    {
      read_params(ip_file,N,d,D,mu,xs,xe,lam,dt,num_of_steps,evolution_order,xsnap,i_start);
      Hsngi.set_x(xsnap);
      Hsngi.updateHMPO();
    }
      
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
  
  //Once the inital MPS is fixed, make a copy in the auxiliary MPS
  aux = *mps;
   
  //Get the right ramping-function such that I do not have to do if-queries over and over again during the evolution
  double (*ramping)(double,double,int,int); 
  if(evolution_order==1 || evolution_order==3)
    ramping = &ramp_xvalue_linear; 
  else
    ramping = &ramp_xvalue_cubic;   
    
  //Prepare the file to monitor energy and other observables during the evolution
  energy_file.open("efile_taylor.txt");
  energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\tTrunc Err\tC_0\tD (resolved every "<< RESOLUTION << " steps)" << endl;
  
  //File to monitor the norms during the "critical phase" (for debugging)
  ofstream normfile;
  normfile.open("Norm.txt");
  normfile << "#Time\tNorm after H|Psi>\tNorm after P|Psi>"<< endl;
  
  //Determine the inital energy
  cout << "Inital energy: " << contractor.contract(*mps,hamil,*mps) << endl;  
  //Determine the inital condensate
  Hsngi.constructCondensateMPO(Condensate);
  cfraction = contractor.contract(*mps,Condensate,*mps);  
  cout <<"Intial condensate: " << cfraction << " with x-value " << Hsngi.get_x() << endl;
  
  //Evolution via Taylor series with the dedicated operator for1-idtH to only have to optimize once per timestep
  cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
      <<          "     		 %Time evolution with Taylor series%                   " << endl
      <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl; 
  //Actual time evolution
  for(int i=i_start; i<=num_of_steps; ++i)
  {
    mps_old = *mps;
    
    //Get the current value of x according to the chosen ramping
    xvalue = (*ramping)(xs,xe,num_of_steps-1,i-1);
    
    //Get the current value of D according to the chosen ramping
    //D = ramp_bond(D_start,D_end,num_of_steps-1,i-1);
  
    //Get the time evolution operators and the current Hamiltonian
    Hsngi.set_x(xvalue);
    
    //Get the evolution MPo for the current value of x
    Hsngi.EvolutionMPO2(evolutionoperator,dt);
    //print_MPO(evolutionoperator);
    //evolutionoperator.exportForMatlab("EvoMPO.m");
    //Apply it to the MPS and truncate the bond dimension
    errN=-1.0;
    contractor.optimize(evolutionoperator,*mps,aux,D,&errN);
    truncation_errors.push_back(errN);
    /*if(i*dt>LOWER_TIME && i*dt<UPPER_TIME)
    {
      curr_norm=contractor.contract(aux,aux);
      cout << "Norm after Evolution operator: " << curr_norm.re  << endl;
      cout << "Truncation Error:              " << errN  << endl;
      normfile << scientific << setprecision(10) << i*dt << "\t" << curr_norm << "\t";
    }*/
    aux.setNormFact(1.0);
    *mps = aux;
    
    //Project back into the Gauss Law fulfilling subspace if desired
    if(evolution_order==3 || evolution_order==4)
    {
      //Reset errorflag and project back
      errN=-1.0;      
      contractor.optimize(projector,aux,*mps,D,&errN);
      //cout << "Norm aux after projection: " << contractor.contract(aux,aux) << endl;
      truncation_errors.push_back(errN);
      curr_norm = contractor.contract(*mps,*mps);
      /*if(curr_norm.re < 0.99 || (i*dt>LOWER_TIME && i*dt<UPPER_TIME))
      {
	 cout << "Warning, norm after projecting is: " << curr_norm << endl;
	 if(i*dt>LOWER_TIME && i*dt<UPPER_TIME)
	   normfile << curr_norm << endl;
      }*/
      //Renormalize
      mps->setNormFact(1.0); 
    }
    
    //if((i%RESOLUTION)==0 || (i*dt>LOWER_TIME && i*dt<UPPER_TIME))
    if((i%RESOLUTION)==0)
    {
      //Get Hamiltonian for current value of x
      Hsngi.updateHMPO();   
      
      //Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep, penalty energy, condensate fraction and charge
      curr_overlap = contractor.contract(mps_old,*mps);   
      curr_energy = contractor.contract(*mps,hamil,*mps);
      curr_norm = contractor.contract(*mps,*mps);      
      penalty_energy = contractor.contract(*mps,Penalty,*mps);
      //Since the condensate MPO is x-dependend, I have to reconstruct it every time
      Hsngi.constructCondensateMPO(Condensate);
      cfraction = contractor.contract(*mps,Condensate,*mps);
      curr_charge = contractor.contract(*mps,chargeoperator,*mps);
      
      //Evaluate truncation errors
      errtot = sum_abs(truncation_errors);
      truncation_errors.clear();
      
      //Some output for the user
      cout << "Step:             " << i << endl;
      cout << "x:                " << Hsngi.get_x() << endl;
      cout << "D:                " << D << endl;
      cout << "Overlap:          " << abs(curr_overlap) << endl;
      cout << "Norm:             " << curr_norm.re << endl;
      cout << "Energy:           " << curr_energy << endl;
      cout << "Truncation error: " << errtot << endl;
      cout << endl;
      
      //Save data 
      energy_file << scientific << setprecision(10) 
                  << i*dt << "\t" 
		  << curr_energy.re  << "\t" 
		  << curr_norm.re << "\t" 
		  << penalty_energy << "\t" 
		  << E0 <<"\t" 
		  << cfraction/sqrt(Hsngi.get_x())<< "\t" 
		  << Hsngi.get_x() << "\t" 
		  << errtot << "\t" 
		  << curr_charge << "\t"
		  << D << "\t"
		  << abs(curr_overlap) << endl;
      
      //In case I want to save intermediate results
      if(SAFTEYMODE && i<num_of_steps && ((i%(INTERVAL*RESOLUTION))==0))
      {
	//Save the ground state
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
	GS_name = GS_name+"N0"+convert.str(); 
	convert.str("");
	convert.clear();
	convert << mu;
	GS_name = GS_name+"mu"+convert.str(); 
	convert.str("");
	convert.clear();
	convert << i;
	GS_name = GS_name+"i"+convert.str()+".dat"; 
	
	mps->exportMPS(GS_name.c_str());
	
	//Furthermore save a file giving the exact values for this run, I would need to restart the evolution
	string param_name="parameter";
	convert.str("");
	convert.clear();
	convert << i;
	param_name = param_name+"i"+convert.str()+".dat"; 
	
	ofstream paramfile;
	paramfile.open(param_name.c_str());
	paramfile << "N = " << N << endl;
	paramfile << "d = " << d << endl;
	paramfile << "D = " << D << endl;
	paramfile << "mu = " << mu << endl;
	paramfile << "xs = " << xs << endl;
	paramfile << "xe = " << xe << endl;
	paramfile << "lambda = " << lam << endl;
	paramfile << "dt = " << dt << endl;
	paramfile << "steps = " << num_of_steps << endl;
	paramfile << "evo_order = " << evolution_order << endl;
	paramfile << scientific << setprecision(11) ;
	paramfile << "x_save = " << xvalue << endl;
	paramfile << "i_save = " << i << endl;	  
	paramfile.close();
	
	cout << "Successfully saved a snapshot" << endl << endl;
      }
    }		
  }
  //Close energy file
  energy_file.close();  
  
  //Close norm file
  normfile.close();
  
  //Get final energy value (in case num_of_steps is not a multiple of the chosen resolution)
  Hsngi.updateHMPO();
  E0 = (contractor.contract(*mps,hamil,*mps)).re;  
  
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
  GS_name = GS_name+"N0"+convert.str(); 
  convert.str("");
  convert.clear();
  convert << mu;
  GS_name = GS_name+"mu"+convert.str()+".dat"; 
  
  time(&start_local);
  mps->exportMPS(GS_name.c_str());
  time(&end_local);
  writetime_gs = difftime(end_local,start_local);
  cout << "Time for writing GS: " << writetime_gs << " s" << endl;
  
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
  
  //Try to calculate the charge
  complex_t charge_0=contractor.contract(*mps,chargeoperator,*mps);
  cout << "Charge: " << charge_0.re << endl << endl;  
  
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
    benchmarking << "#N\tD\tN_0\tT_total(s)\tT_perstep(s)" << endl;
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
  cout << "MPS of length " << mps.getLength() << endl;
  for(int i=0; i<mps.getLength(); i++)
  {
    cout << "[" << i << "]=" <<  (mps.getA(i)).getA() << endl;
  }
}
void print_MPO(MPO &mpo)
{
  cout << "MPO of length " << mpo.getLength() << endl;
  for(int i=0; i<mpo.getLength(); i++)
  {
    cout << "[" << i << "]=" <<  (mpo.getOp(i)).getFullData() << endl;
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
    
    cout << "Read N=" << N << endl;
    cout << "Read d=" << d<< endl;
    cout << "Read D=" << D<< endl;
    cout << "Read mu=" << mu<< endl;
    cout << "Read xs=" << xs<< endl;
    cout << "Read xe=" << xe<< endl;
    cout << "Read lamda=" << lam<< endl;
    cout << "Read dt=" << dt<< endl;
    cout << "Read steps=" << steps<< endl;
    cout << "Read evo_order=" << evolution_order<< endl;
    cout << "Read xjump=" << xjump<< endl;
    cout << "Read ijump=" << ijump<< endl << endl;
    
    pfile.close();
    return 0;
  }
}