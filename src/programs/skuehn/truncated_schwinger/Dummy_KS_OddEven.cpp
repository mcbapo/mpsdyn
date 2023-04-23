/**
   \file KS_OddEven.cpp
   Realisation of a simulation of the time evolution under the fundamental Hamiltonian for 1+1d cQED as proposed in Ignacios review. This file is dedicated to the time evolution via Odd-Even splitting of the time evolution operator. Basically this files contains the versions from CQED_AdiabaticEvolution.cpp in a slightly optimized and cleaned up form.
   
  \param <N> (int) Number of Fermions
  \param <d> (int) Physical dimension of Bosons in between
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <x> (double) Coupling constant
  \param <e> (double) Constant of the Operator \f$ \Theta(n)\f$
  \param <lambda> (double) Constant for the penalty-term in the Hamiltonian
   
  \author Stefan Kühn
  \date 08/04/2014

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
#include "Indices.h"
#include "misc.h"

#include "KogutSusskindHamiltonian.h"
#include "SHGIGaussLawAngular.h"
#include "ramp_xvalue.h"

#define NUM_OF_PARAMS 11
#define FILENAME "KS_dummy.txt"
#define KILLFILE "crash_dump.txt"
#define RESOLUTION 10
#define BU_RESOLUTION (100*RESOLUTION)
#define EPS 1E-10
#define DMIN 3
//#define THRESHOLD 0.99
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
    cout << "Model: " << "Select model Hamiltonian {cQED,Zn,Zncgl,cQEDnoise,Zncglnoise}" << endl
	 << "N:     " << "Number of spins" << endl
	 << "d:     " << "Bosonic Hilbertspace dimension" << endl
	 << "D:     " << "Bond dimension of the MPS" << endl
	 << "µ:     " << "Dimensionless constant" << endl
	 << "xs:    " << "Coupling constant to start ramping with" << endl
	 << "xe:    " << "Coupling constant to end ramping with" << endl
	 << "e:     " << "Constant for gauge field" << endl
	 << "la:    " << "Constant for penalty term" << endl
	 << "dt:    " << "Timestep size" << endl
	 << "Nt:    " << "Number of timesteps which should be calculated" << endl
	 << "ig:    " << "True if inital guess for the ground state should be read form HDD (optional)" << endl
	 << "ac:    " << "Accuracy for contractor (optional)" << endl
	 << "O:     " << "Order of the time evolution (optional, by default first order)" << endl
	 << "dx:    " << "Stepsize at which intermediate values should be taken (optional)" << endl
	 << "xt:    " << "Threshold above which intermediate steps are taken (optional)" << endl
	 << "noise: " << "Strength of the noise (optional)" << endl
	 << "ov:    " << "Compute overlap with given reference data (optional)" << endl;
    return -1;
  }
  
  //Model
  string model=argv[++cntr];
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
  //String needed for conversion
  ostringstream convert;  
  //Strength for the noise
  double noise_strength=1.0;
  
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
    noise_strength=atof(argv[++cntr]);
    
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file " << ig_file << endl;
  }
  
  //Custom accuracy for contractor, file for inital state, order for time evolution, xt, get_overlap and stepsize for intermediate snapshots is given
  if(argc == NUM_OF_PARAMS + 8)
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    evolution_order = atoi(argv[++cntr]);
    dx = atof(argv[++cntr]);
    dx_flag = true;
    xt = atof(argv[++cntr]);
    noise_strength=atof(argv[++cntr]);
    get_overlap=true;
    path_to_reference = argv[++cntr];
    
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file " << ig_file << endl;
  }
  
    
  //Print parameter which were set
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
       << "Nt = " << num_of_steps << endl;
  if (argc>=NUM_OF_PARAMS+2){
    cout << "ig = " << ig << endl;
  }
  if (argc>=NUM_OF_PARAMS+3){
     cout << "acc= " << desired_accuracy << endl;
  }
  if (argc>=NUM_OF_PARAMS+4){
     cout << "O  = " << evolution_order << endl;
  }
  if (argc>=NUM_OF_PARAMS+5){
     cout << "dx = " << dx << endl;
  }
  if(argc>=NUM_OF_PARAMS+6){
     cout << "xt = " << xt << endl;
  }
  if(argc>=NUM_OF_PARAMS+7){
     cout << "n =  " << noise_strength << endl;
  }
  if(argc>=NUM_OF_PARAMS+8){
     cout << "pwd= " << path_to_reference << endl;
  }
  cout << "-----------" << endl << endl;
    
       
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
  
  //Get the corresponding Hamiltonian, to avoid lower und upper case issues, convert everything to lower case
  transform(model.begin(),model.end(),model.begin(),::tolower);
  Model mod;
  if(model.compare("zn")==0)
    mod = Zn;
  else if(model.compare("cqed")==0)
    mod = cQED;
  else if(model.compare("cqednoise")==0)
    mod = cQEDnoise;
  else if(model.compare("zncgl")==0)
    mod = Zncgl;
  else if(model.compare("zncglnoise")==0)
    mod = Zncglnoise;
  else
  {
    cout << "Unknown model, program will be aborted..." << endl;
    exit(666);
  }
  KogutSusskindHamiltonian Hsngi(N,d,D,mu,xs,lam,mod,noise_strength);
  
  
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
  //Filepointer for the inital guess
  ifstream inital_guess;
  //Charge MPO and variable for the charge
  MPO chargeop(1);
  Hsngi.constructChargeMPO(chargeop);  
  complex_t charge_0;
  //Ofstream for the file to hold the final results
  ofstream resfile;  
  
  //Try to open file for inital state
  inital_guess.open(ig_file.c_str());
  cout << "ig_file" << ig_file << endl;
  //In case file for inital state exists read it and adjust bond dimension, otherwise start with predefined inital state
  if(inital_guess.is_open() && ig)
  {
    //Case in which an inital guess exists and will be read
    cout << "Found initial guess for the groundstate" << endl;
    mps = new MPS();
    mps->importMPS(ig_file.c_str());
    cout << "Bond dimension after import: "<< mps->getBond() << endl;
    cout << "Length after import: " << mps->getLength() << endl;
    int status = read_params("bu_params.txt",N,d,D,mu,xs,xe,lam,dt,num_of_steps,evolution_order,xjump,i_start);
    //Check if the inital guess is a backup from an older simulation and if I have to adjust the threshold value xt at which data is taken
    if(status==0)
    {
      //Get the old value of x
      double xold = ramp_xvalue_time(xs,xe,num_of_steps*dt,((double)i_start)*dt+dt/2.0);
      //Get last value at which data was taken before a snapshot was saved and add dx, such that one has the next value at which data is required
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
  
  //Second order time evolution
  //This version has the option to specify certain values of x which are then used to take data besides the final snapshot
  if(evolution_order==13)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<        "                  %Starting second order time evolution%                     " << endl
       <<        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    //Construct the two additional MPOs for the halfsteps
    MPO exp_even_half(1),exp_odd_half(1);
 
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Flag if I have to take some data
    bool take_data_flag=false;
    
    //Save the inital bond dimension (which will be the maximum bond dimension reached)
    int Dmax=D;
    
    //Construct the penalty MPO according to the chosen model(with fixed strength)
    if((mod==Zncgl)||(mod==Zncglnoise))
      Hsngi.constructCglPenaltyMPO(Penalty,1.0);
    else
      Hsngi.constructPenaltyMPO(Penalty,1.0);
    
    //Actual time evolution
    for(int i=i_start-1; i<num_of_steps; ++i)
    {
      mps_old = *mps;
      
      //Get the current value of x
      //xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);  
      xvalue = ramp_xvalue_time(xs,xe,num_of_steps*dt,((double)i)*dt+dt/2.0,3.0);
      //xvalue = ramp_xvalue_time_hybrid(xs,xe,num_of_steps*dt,((double)i)*dt+dt/2.0,150.0);
      //cout << "x =  " << xvalue << endl;
      
      //Case that the user speciefied some value of x at which he wants to have a snapshot
      if(dx_flag && xvalue>=xt)
      {
	xt += dx;
	take_data_flag = true;
      }
   
      //Get the time evolution operators and the current Hamiltonian
      Hsngi.set_x(xvalue);
      //Ramp the noise with the xvalue
      Hsngi.set_noise(noise_strength*xvalue);
      
      //Reconstruct the data points
      if(take_data_flag)
      {
	//Get the current MPO for the Hamilton
	Hsngi.updateHMPO();
	
	//State I have to read
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
	mps->importMPS(GS_name.c_str());
	cout << "Imported " << GS_name << endl;

	//Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep     
	curr_energy = contractor.contract(*mps,hamil,*mps);
	curr_norm = contractor.contract(*mps,*mps);	
	penalty_energy = contractor.contract(*mps,Penalty,*mps);
	charge_0 = contractor.contract(*mps,chargeop,*mps);
	
	Hsngi.constructCondensateMPO(Condensate);
	cfraction = contractor.contract(*mps,Condensate,*mps);
	
	//Check truncation errors
	
	//Output
	cout << scientific;
	cout << "Step:             " << i+1 << endl
	     << "x:                " << xvalue << endl
	     << "D:                " << D << endl
	     << "Norm:             " << curr_norm << endl
	     << "Energy:           " << curr_energy << endl
	     << "Noise:            " << Hsngi.get_noise()  << endl;
	 
	curr_overlap = -1.2*ONE_c;	
	
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
      }
    }
  }
   
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
  mps->importMPS(GS_name.c_str());
  cout << "Imported " << GS_name << endl;   
        
  //Get the current MPO for the Hamilton
  Hsngi.updateHMPO();
      
  //Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
  curr_energy = contractor.contract(*mps,hamil,*mps);
  curr_norm = contractor.contract(*mps,*mps);	
  penalty_energy = contractor.contract(*mps,Penalty,*mps);
  charge_0 = contractor.contract(*mps,chargeop,*mps);
  
  Hsngi.constructCondensateMPO(Condensate);
  cfraction = contractor.contract(*mps,Condensate,*mps);
  
  E0 = curr_energy.re;
   
  
  cout << "Ratio Penalty in percent/Groundstate energy: " << penalty_energy.re/E0*100.0<<" %" << endl<<endl; 
  
  //Try to calculate the charge
  charge_0=contractor.contract(*mps,chargeop,*mps);
  cout << "Charge Groundstate " << charge_0.re << endl;  
  
  //Calculate the terms for the Gauss Law and see if it is fullfilled
  ofstream Lnfile;
  Lnfile.open("Ln_dummy.txt", ios_base::app);
  complex_t Lnvalue;
  MPO Ln(1);
  cout << endl;
  for(int i=1; i<L-1;i += 2){
    Hsngi.constructLopMPO(Ln,i);
    Lnvalue = contractor.contract(*mps,Ln,*mps);
    cout << "--> Expectation value L("<< (i+1)/2 << ") = " << Lnvalue<< endl;
    Lnfile << Lnvalue.re << "\t";

  }
  Lnfile << endl;
  Lnfile.close();
  
  
  ofstream Gaussfile;
  Gaussfile.open("Gaussfile_dummy.txt", ios_base::app);
  complex_t Gaussvalue;
  MPO Gauss(1);
  cout << endl;
  for(int i=0; i<L;i += 2){
    Hsngi.constructGaussMPO(Gauss,i);
    Gaussvalue = contractor.contract(*mps,Gauss,*mps);
    cout << "--> Expectation value 0.5*[sigma_z+(-1)^"<< i/2+1 << "] = " << Gaussvalue<< endl;
    Gaussfile << Gaussvalue.re << "\t";
  }
  Gaussfile << endl;
  Gaussfile.close();
  
  
  //Calculate the z-Component of the spins
  ofstream spinfile;
  spinfile.open("Spin_dummy.txt", ios_base::app);
  complex_t Zcomponent;
  MPO Sigmaz(1);
  cout << endl;
  for(int i=0; i<L;i += 2){
    Hsngi.constructsigmazMPO(Sigmaz,i);
    Zcomponent = contractor.contract(*mps,Sigmaz,*mps);
    cout << "--> Expectation value sigma_z("<< i/2+1 << ") = " << Zcomponent << endl;
    spinfile << Zcomponent.re << "\t";
  }
  spinfile << endl;
  spinfile.close(); 
  
  
  //Calculate the condensate farction
  cout << endl;
  Hsngi.constructCondensateMPO(Condensate);
  cfraction = contractor.contract(*mps,Condensate,*mps);
  cout << "Condensate: "<< cfraction << endl << endl;  
  
  //Ouput the results up to now (sometimes the computation of the exicted state(s) fails and the Contractor routine quits the program, therefore all results which were not stored would be lost). To avoid this, save the first part of the results now
  cout << "Saving results..." << endl;
  if(file_exists(FILENAME))
  {
    resfile.open(FILENAME,ios_base::app);  
    resfile << endl << scientific << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << xe << "\t" << e << "\t" << lam << scientific <<"\t" << E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" <<penalty_energy.re << "\t";
  }
  else
  {
    resfile.open(FILENAME,ios_base::app);  
    resfile << "#N" << "\t" << "d" << "\t" << "D" << "\t" << "mu" << "\t" << "x" << "\t" << "e" << "\t" << "lambda" << "\t"<< "E0" << "\t" << "C0" << "\t" << "Cfrac" << "\t" << "E1" << "\t" << "C1" << "\tE_pen"<< endl;
    resfile << scientific << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << xe << "\t" << e << "\t" << lam << scientific << "\t"<< E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" << penalty_energy.re << "\t";
  }
  resfile.close();
  
  //Get the final overlap in case it is desired
  //The reason I do it in such a cumbersome way is the following: If something goes wrong while trying to open the last reference state to compute the overlap, my program crashes and I loose all the other data not just the last overlap. Right now importMPS does not allow me to catch the error if opening the reference state fails. Therefore I have no other choice than writing everything to a file before starting to compute the last overlap, to be sure to save my data.
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
    //New convetion, d is the true Hilbertspace dimension
    //convert << d;
    //Old convetion, d is the boson number
    convert << d-1;
    reference_name = reference_name+"N0"+convert.str(); 
    convert.str("");
    convert.clear();
    //convert << (int) Hsngi.get_x();
    convert << (int) round(Hsngi.get_x());
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
  

  
  
  //Free Memory
  delete[] phys_dims;
  delete[] bond_dim; 
  mps->clear();
  delete mps;
  
 
  
  //End Time
  time(&end);
  double runtime_total = difftime(end,start);
  cout << "Run time: " << runtime_total <<"s" <<  endl;    
  
  cout << "End of program" << endl;
  return 0;
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

