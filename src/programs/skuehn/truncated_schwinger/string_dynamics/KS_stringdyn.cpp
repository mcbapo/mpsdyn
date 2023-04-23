/**
   \file KS_stringdyn.cpp
   Program to look at certain string dynamics of the truncated Kogut Susskind Hamiltonian
   
  \param <N> (int) Number of Fermions
  \param <d> (int) Physical dimension of Bosons in between
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <x> (double) Coupling constant
  \param <lambda> (double) Constant for the penalty-term in the Hamiltonian
   
  \author Stefan Kühn
  \date 10/06/2015

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

#define NUM_OF_PARAMS 9
#define FILENAME "KS.txt"
#define KILLFILE "crash_dump.txt"
#define FILENAME_CONFIG "KS_config.txt"
#define FILENAME_COEFFIS "KScoefficients.txt"
#define RESOLUTION 5
#define BU_RESOLUTION (100*RESOLUTION)
#define EPS 1E-10
#define DMIN 3




using namespace std;
using namespace shrt;

//bool file_exists(string filename);
void print_MPS(MPS);
double sum(vector<double>);
double sum_abs(vector<double>);
//Function to read a file containing some parameters
int read_params(string paramfile,int &, int &, int &, double &, double &,  double &, double &, int &);
//Function to save a temporary results to be able to use it again
int save_temp_result(MPS &,int, int, int, double, double, double, double, int);

//Determine the current configuration
vector<complex_t> current_configuration(vector<MPO*> &Spin,vector<MPO*> &Gauge,MPS &mps,Contractor& contractor,int N);
/** Save the current configuration to a given file */
int save_configuration(string filename, vector<complex_t> config);

//Get the distribution of entanglement entropy
vector<complex_t> get_EntanglementEntropy(Contractor &contractor, MPS mps, int N);



int main(int argc,const char* argv[]){
  int cntr=0;
  
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << "Model N d D mu x la dt Nt ig ac noise" << endl;
    cout << "Model: " << "Select model Hamiltonian {cQED,Zn,Zncgl,cQEDnoise,Zncglnoise}" << endl
	 << "N:     " << "Number of spins" << endl
	 << "d:     " << "Bosonic Hilbertspace dimension" << endl
	 << "D:     " << "Bond dimension of the MPS" << endl
	 << "µ:     " << "Dimensionless constant" << endl
	 << "x:     " << "Interaction strength" << endl
	 << "la:    " << "Constant for penalty term" << endl
	 << "dt:    " << "Timestep size" << endl
	 << "Nt:    " << "Number of timesteps which should be calculated" << endl
	 << "ig:    " << "True if inital guess for the ground state should be read form HDD (optional)" << endl
	 << "ac:    " << "Accuracy for contractor (optional)" << endl
	 << "noise: " << "Strength of the noise (optional)" << endl;
    return -1;
  }
  
  //Model
  string model=argv[++cntr];
  //Number of spin sites
  int N=atoi(argv[++cntr]);
  //Total number of sites
  int L=2*N-1;
  //Total number fo bosons living on a link
  int d=atoi(argv[++cntr]);
  //Bond dimension
  int D=atoi(argv[++cntr]);
  //Mass
  double mu=atof(argv[++cntr]);
  //Coupling constant
  double x=atof(argv[++cntr]);
  //Time step to start with
  int i_start = 1;
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
  //Variables for time measurement
  time_t start,end;
  time_t start_local,end_local;
  double readtime_gs=0.0,writetime_gs=0.0,writetime_gauss=0.0,writetime_spin=0.0,writetime_ln=0.0,writetime_total=0.0,writetime_data=0.0;
  //String needed for conversion
  ostringstream convert;  
  //Strength for the noise
  double noise_strength=1.0;
  //Vector for the configuration
  vector<complex_t> configuration,coefficients;
  
  //Only file for inital guess is given aditionally to all other parameters
  if(argc >= NUM_OF_PARAMS + 2)
  {
    ig = true;
    ig_file = argv[++cntr];
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file: " << ig_file << endl;
  }
  //Custom accuracy for contractor and file for inital state
  if(argc >= NUM_OF_PARAMS + 3)
  {
    set_accuracy = true;    
    desired_accuracy = atof(argv[++cntr]);
  }  
  
  //Custom accuracy for contractor, file for inital state, order for time evolution, xt, get_overlap and stepsize for intermediate snapshots is given
  if(argc >= NUM_OF_PARAMS + 4)
  {
    noise_strength=atof(argv[++cntr]);
  }  
    
  //Print parameter which were set
  cout << "Parameter: " << endl
       << "-----------" << endl
       << "N  = " << N << endl
       << "d  = " << d << endl
       << "D  = " << D << endl
       << "µ  = " << mu << endl
       << "x  = " << x << endl
       << "la = " << lam << endl
       << "dt = " << dt << endl
       << "Nt = " << num_of_steps << endl;
  if (argc>=NUM_OF_PARAMS+2)
    cout << "ig = " << ig << endl;
  if (argc>=NUM_OF_PARAMS+3)
     cout << "acc= " << desired_accuracy << endl;
  if(argc>=NUM_OF_PARAMS+4)
     cout << "n =  " << noise_strength << endl;
  cout << "-----------" << endl << endl;
  
  //Initialize random seed
  srand (time(NULL));
    
       
  //Start time measurement
  time(&start);
  
  
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
  KogutSusskindHamiltonian Hsngi(N,d,D,mu,x,lam,mod,noise_strength);
  
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
  //File for results
  ofstream file;
  //Filepointer for the inital guess
  ifstream inital_guess;
  //File for reading kvector, string parameters and so on
  ifstream infile;
  //Charge MPO and variable for the charge
  MPO chargeop(1);
  Hsngi.constructChargeMPO(chargeop);  
  complex_t charge_0;
  //Ofstream for the file to hold the final results
  ofstream resfile;  
  //State
  vector<int> dims;
  dims.resize(L,2);
  for(int i=1;i<L; i+=2)
    dims[i]=d;
  MPS*  mps = new MPS(L,D,dims);
  mps->setRandomState();
  
  //Prepare the single body MPOs
  vector<MPO*> Spin_MPOs,Gauge_MPOs;  
  MPO *mpo;
  for(int k=1; k<=N; k++)
  {
    mpo=new MPO(1);
    Hsngi.constructsigmazMPO(*mpo,2*k-2);
    Spin_MPOs.push_back(mpo);    
    mpo=new MPO(1);       
    if(k<N)
    {
      mpo=new MPO(1);
      //Hsngi.constructLopMPO(*mpo,2*k-1);
      Hsngi.constructLopsquareMPO(*mpo,2*k-1);
      Gauge_MPOs.push_back(mpo);       
    }
  }

  //Construct the penalty MPO according to the chosen model(with fixed strength)
  if((mod==Zncgl)||(mod==Zncglnoise))
    Hsngi.constructCglPenaltyMPO(Penalty,1.0);
  else
    Hsngi.constructPenaltyMPO(Penalty,1.0);
  
  
  //Try to open file for inital state
  inital_guess.open(ig_file.c_str());
  cout << "ig_file" << ig_file << endl;
  //In case file for inital state exists read it and adjust bond dimension, otherwise start with predefined inital state
  if(inital_guess.is_open() && ig)
  {
    //Case in which an inital guess exists and will be read
    cout << "Found initial guess for the groundstate" << endl;
    time(&start_local);
    mps->importMPS(ig_file.c_str());
    cout << "Bond dimension after import: "<< mps->getBond() << endl;
    cout << "Length after import: " << mps->getLength() << endl;
    int status = read_params("bu_params.txt",N,d,D,mu,x,lam,dt,num_of_steps);
    //Check if the inital guess is a backup from an older simulation and if I have to adjust the threshold value xt at which data is taken
    i_start++;
   
    //In case an intial guess with bond dimension smaller than set is imported it has to be expanded again
    if(mps->getBond() < D)
    {
      //Increase bond dimension again, in case the inital guess has a smaller bond dimension
      mps->increaseBondDimension(D);
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
    //mps = Hsngi.constructInitalMPS();
    //mps->setRandomState();
    //mps = Hsngi.constructInitalMPS(is_downup);
    //mps = Hsngi.constructInitalMPS(is_upup);
    //mps = Hsngi.constructInitalMPS(is_downdown);
    //mps = Hsngi.constructInitalMPS(is_string,5);        
    //mps = Hsngi.constructInitalMPS(is_pzoller_string);
    
    MPS vacuum(L,D,dims);
    contractor.findGroundState(hamil,D,&E0,vacuum);
    
    //Build a string
    /*MPO StringOp(1);    
    Hsngi.getStringOperator(StringOp,5);
    cout << "Computed ground state" << endl;
    mps->approximate(StringOp,vacuum,D);
    contractor.optimize(StringOp,vacuum,*mps,D);  
    mps->gaugeCond('R',true);*/
    
    //Single Meson
    /*cout << "Starting to build meson" << endl;
    MPS *tmp_mps1;
    MPO MesonOp(1),Momentum(1);
    MPS aux;
    int pos_left=round(N/2.0);
    double kleft=-2.0;
    complex_t pen1, pen2, pen3;
    if(file_exists("kvecs.txt"))
    {
      infile.open("kvecs.txt");
      infile >> kleft ;
      cout << "Read k-vector:" << endl << "--> kleft=" << kleft << endl;
    }
    Hsngi.getMesonAntiMesonOperator(MesonOp,pos_left,kleft,2.0);
    tmp_mps1=new MPS(vacuum);	tmp_mps1->approximate(MesonOp,vacuum,D);
    contractor.optimize(MesonOp,vacuum,*tmp_mps1,D);
    tmp_mps1->gaugeCond('R',true);
    pen1=contractor.contract(*tmp_mps1,Penalty,*tmp_mps1);
    cout << "Penalty energy meson state 1: " << pen1 << endl;
    Hsngi.constructMomentumMPO(Momentum);
    cout << "Momentum meson state: " << contractor.contract(*tmp_mps1,Momentum,*tmp_mps1) << endl;
    cout << "Momentum vacuum: " << contractor.contract(vacuum,Momentum,vacuum) << endl;
    cout << "Momentum square meson state: " << contractor.contract2(Momentum,*tmp_mps1)/contractor.contract(*tmp_mps1,*tmp_mps1) << endl;
    cout << "Momentum square vacuum: " << contractor.contract2(Momentum,vacuum)/contractor.contract(vacuum,vacuum) << endl;
    *mps=MPS(*tmp_mps1);
    if(pen1.re>1E-3)
    {
      cout << "Error, approximation is in the wrong sector, will abort program..." << endl;
      exit(666);
    }    
    tmp_mps1->clear(); delete tmp_mps1;*/
    
    //Two Mesons with momentum
    cout << "Starting to build meson" << endl;
    MPS *tmp_mps1,*tmp_mps2;
    MPO MesonOp(1),Momentum(1);
    MPS aux;
    int pos_left=round(N/2.0)-10, pos_right=round(N/2.0)+10;
    double kleft=-2.0,kright=2.0;
    bool include_antimesons=false;
    complex_t pen1, pen2, pen3;
    if(file_exists("kvecs.txt"))
    {
      infile.open("kvecs.txt");
      infile >> kleft >> kright;
      if(!infile.eof())
	infile >> include_antimesons;
      cout << "Read k-vectors:" << endl << "--> kleft=" << kleft << endl << "--> kright=" << kright << endl;
      if(include_antimesons)	
	cout << "Will include anti mesons" << endl;
	
    }
    if(include_antimesons)
      Hsngi.getMesonAntiMesonOperator(MesonOp,pos_left,kleft,2.0);
    else
      Hsngi.getMesonOperator(MesonOp,pos_left,kleft,2.0);
    tmp_mps1=new MPS(vacuum);	tmp_mps1->approximate(MesonOp,vacuum,D);
    contractor.optimize(MesonOp,vacuum,*tmp_mps1,D);
    if(include_antimesons)
      Hsngi.getMesonAntiMesonOperator(MesonOp,pos_right,kright,2.0);
    else
      Hsngi.getMesonOperator(MesonOp,pos_right,kright,2.0);
    tmp_mps2=new MPS(vacuum);	tmp_mps2->approximate(MesonOp,vacuum,D);
    contractor.optimize(MesonOp,vacuum,*tmp_mps2,D);
    pen1=contractor.contract(*tmp_mps1,Penalty,*tmp_mps1);
    pen2=contractor.contract(*tmp_mps2,Penalty,*tmp_mps2);
    cout << "Penalty energy meson state 1: " << pen1 << endl;
    cout << "Penalty energy meson state 2: " << pen2 << endl; 
    Hsngi.constructMomentumMPO(Momentum);
    cout << "Momentum meson state 1: " << contractor.contract(*tmp_mps1,Momentum,*tmp_mps1) << endl;
    cout << "Momentum meson state 2: " << contractor.contract(*tmp_mps2,Momentum,*tmp_mps2) << endl;
    cout << "Momentum vacuum: " << contractor.contract(vacuum,Momentum,vacuum) << endl;
    cout << "Momentum square meson state 1: " << contractor.contract2(Momentum,*tmp_mps1)/contractor.contract(*tmp_mps1,*tmp_mps1) << endl;
    cout << "Momentum square meson state 2: " << contractor.contract2(Momentum,*tmp_mps2)/contractor.contract(*tmp_mps2,*tmp_mps2) << endl;
    cout << "Momentum square vacuum: " << contractor.contract2(Momentum,vacuum)/contractor.contract(vacuum,vacuum) << endl;
    vector<const MPS*> states; states.push_back(tmp_mps1); states.push_back(tmp_mps2);
    vector<complex_t> phases; 
    phases.push_back(ONE_c); phases.push_back(ONE_c);
    aux=MPS(vacuum); aux.setRandomState();
    contractor.optimizeSum(states,phases,aux,D);
    aux.gaugeCond('R',true);
    *mps=MPS(aux);
    pen3=contractor.contract(*mps,Penalty,*mps);
    cout << "Penalty energy two meson state: " << pen3 << endl;
    cout << "Momentum two meson state: " << contractor.contract(*mps,Momentum,*mps) << endl;
    cout << "Momentum square two meson state: " << contractor.contract2(Momentum,*mps) << endl;
    if(pen1.re>1E-3 || pen2.re>1E-3 || pen3.re>1E-3)
    {
      cout << "Error, approximation is in the wrong sector, will abort program..." << endl;
      exit(666);
    }    
    tmp_mps1->clear(); delete tmp_mps1;
    tmp_mps2->clear(); delete tmp_mps2;
    
    //Two strings of length one taken with relative phase
    /*cout << "Starting to build meson" << endl;
    MPS *tmp_mps1,*tmp_mps2;
    MPO StringOp(1),Momentum(1);
    MPS aux;
    int pos_left=round(N/4.0)+1, pos_right=round(3.0*N/4.0)-1;
    double kleft=-2.0,kright=2.0;
    complex_t pen1, pen2, pen3;
    if(file_exists("kvecs.txt"))
    {
      infile.open("kvecs.txt");
      infile >> kleft >> kright;
      cout << "Read k-vectors:" << endl << "--> kleft=" << kleft << endl << "--> kright=" << kright << endl;
    }
    Hsngi.getStringOperator(StringOp,1,pos_left);
    tmp_mps1=new MPS(vacuum);	tmp_mps1->approximate(StringOp,vacuum,D);
    contractor.optimize(StringOp,vacuum,*tmp_mps1,D);
    Hsngi.getStringOperator(StringOp,1,pos_right);
    tmp_mps2=new MPS(vacuum);	tmp_mps2->approximate(StringOp,vacuum,D);
    contractor.optimize(StringOp,vacuum,*tmp_mps2,D);
    pen1=contractor.contract(*tmp_mps1,Penalty,*tmp_mps1);
    pen2=contractor.contract(*tmp_mps2,Penalty,*tmp_mps2);
    cout << "Penalty energy string state 1: " << pen1 << endl;
    cout << "Penalty energy string state 2: " << pen2 << endl; 
    Hsngi.constructMomentumMPO(Momentum);
    cout << "Momentum string state 1: " << contractor.contract(*tmp_mps1,Momentum,*tmp_mps1) << endl;
    cout << "Momentum string state 2: " << contractor.contract(*tmp_mps2,Momentum,*tmp_mps2) << endl;
    cout << "Momentum vacuum: " << contractor.contract(vacuum,Momentum,vacuum) << endl;
    cout << "Momentum square string state 1: " << contractor.contract2(Momentum,*tmp_mps1)/contractor.contract(*tmp_mps1,*tmp_mps1) << endl;
    cout << "Momentum square string state 2: " << contractor.contract2(Momentum,*tmp_mps2)/contractor.contract(*tmp_mps2,*tmp_mps2) << endl;
    cout << "Momentum square vacuum: " << contractor.contract2(Momentum,vacuum)/contractor.contract(vacuum,vacuum) << endl;
    vector<const MPS*> states; states.push_back(tmp_mps1); states.push_back(tmp_mps2);
    vector<complex_t> phases; 
    phases.push_back(exp(-I_c*kleft*M_PIl)); phases.push_back(exp(I_c*kright*M_PIl));
    aux=MPS(vacuum); aux.setRandomState();
    contractor.optimizeSum(states,phases,aux,D);
    aux.gaugeCond('R',true);
    *mps=MPS(aux);
    pen3=contractor.contract(*mps,Penalty,*mps);
    cout << "Penalty energy two string state: " << pen3 << endl;
    cout << "Momentum two string state: " << contractor.contract(*mps,Momentum,*mps) << endl;
    cout << "Momentum square two string state: " << contractor.contract2(Momentum,*mps) << endl;
    if(pen1.re>1E-3 || pen2.re>1E-3 || pen3.re>1E-3)
    {
      cout << "Error, approximation is in the wrong sector, will abort program..." << endl;
      exit(666);
    }    
    tmp_mps1->clear(); delete tmp_mps1;
    tmp_mps2->clear(); delete tmp_mps2;*/
  } 
  
  
  //Determine the inital energy
  cout << "Energy before time evolution " << contractor.contract(*mps,hamil,*mps) << endl;  
  //Determine condensate before evolution
  Hsngi.constructCondensateMPO(Condensate);
  cfraction = contractor.contract(*mps,Condensate,*mps);  
  cout <<"Intial condensate: " << cfraction << " with x-value " << Hsngi.get_x() << endl;
  
  //Determine and save the configuration once before the evolution
  configuration=current_configuration(Spin_MPOs,Gauge_MPOs,*mps,contractor,N);  
  save_configuration(FILENAME_CONFIG,configuration);
  
  //Determine the entanglement entropy once before the evolution and save it
  coefficients = get_EntanglementEntropy(contractor,*mps,N);
  file.open(FILENAME_COEFFIS, ios_base::app);
  file << 0 << "\t" << 0 ;
  for(int ind=0; ind<coefficients.size(); ind++)
    file << "\t" << coefficients[ind];
  file << endl;
  file.close();
  
  
  
  //Second order time evolution
  //This version has the option to specify certain values of x which are then used to take data besides the final snapshot
  cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
      <<         "                  %Starting second order time evolution%                     " << endl
      <<         "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  //Construct the two additional MPOs for the halfsteps
  MPO exp_even_half(1),exp_odd_half(1);
  
  energy_file.open("KS_efile.txt");
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
  double errN,max_err,total_max_err=0.0;
  double errtot;
  
  
  //Save the inital bond dimension (which will be the maximum bond dimension reached)
  int Dmax=D;
  
  //Time evolution operators
  Hsngi.constructUevenMPO(exp_even,dt);
  Hsngi.constructUoddMPO(exp_odd_half,dt/2.0); 
  
  //Actual time evolution
  for(int i=i_start-1; i<num_of_steps; ++i)
  {
    mps_old = *mps;    
    
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
      //Hsngi.updateHMPO();
      
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
      max_err =*std::max_element(truncation_errors.begin(),truncation_errors.end());
      truncation_errors.clear();
      
      if(max_err>total_max_err)
	total_max_err=max_err;
      
      //Output
      cout << scientific;
      cout  << "Step:              " << i+1 << endl
	    << "x:                 " << x << endl
	    << "D:                 " << D << endl
	    << "Overlap:           " << abs(curr_overlap) << endl
	    << "Norm:              " << curr_norm << endl
	    << "Energy:            " << curr_energy << endl
	    << "Noise:             " << Hsngi.get_noise() << endl
	    << "Truncation Errors: " << errtot << endl
	    << "Maximum Error:     " << total_max_err << endl
	    << "Penalty Energy:    " << penalty_energy << endl << endl;
	    
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
      
      //Take the data
      if(file_exists(FILENAME))
      {
	resfile.open(FILENAME,ios_base::app);  
	resfile << endl << scientific << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << Hsngi.get_x() << "\t" << lam << scientific <<"\t" << curr_energy.re << "\t" << charge_0.re << "\t" << cfraction.re << "\t" <<penalty_energy.re << "\t" << curr_overlap << "\t";
      }
      else
      {
	resfile.open(FILENAME,ios_base::app);  
	resfile << "#N" << "\t" << "d" << "\t" << "D" << "\t" << "mu" << "\t" << "x" << "\t" << "e" << "\t" << "lambda" << "\t"<< "E0" << "\t" << "C0" << "\t" << "Cfrac" << "\t" << "E1" << "\t" << "C1" << "\tE_pen"<< endl;
	resfile << scientific << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << Hsngi.get_x() << "\t" << lam << scientific << "\t"<< curr_energy.re << "\t" << charge_0.re << "\t" << cfraction.re << "\t" << penalty_energy.re << "\t" << curr_overlap << "\t";  
	
      }
      resfile.close();
      
      configuration=current_configuration(Spin_MPOs,Gauge_MPOs,*mps,contractor,N);  
      save_configuration(FILENAME_CONFIG,configuration);
      
      coefficients = get_EntanglementEntropy(contractor,*mps,N);
      file.open(FILENAME_COEFFIS, ios_base::app);
      file << i << "\t" << i*dt ;
      for(int ind=0; ind<coefficients.size(); ind++)
	file << "\t" << coefficients[ind];
      file << endl;
      file.close();
      
      //Do a backup of the current result if desired 
      //Don't do a backup for the last step, the reason is if I do that, the loop wont run and therefore a lot of quantities will not be initialized properly which results in quite some mess in the output file
      //TODO: Get this cleaner and manage to also recompute quantities out of a finished simulation
      if(((i+1)%BU_RESOLUTION)==0)
      {
	/*I save a snapshot for i+1 because:
	  * - I already completed the time step with i, the next to proceed is i+1
	  * - after reading the input file I add +1 to i
	  * - this loop starts with i_start-1, so effectively I start with value of i, I put in the function save_temp_result
	  * - since I don't want to do timesteps twice, I add +1 to i, to have it correct again*/
	//TODO: a bit convoluted, should be cleared up a bit	  
	save_temp_result(*mps,N,d,D,mu,x,lam,dt,num_of_steps);
      }
	
	//Save snapshot
	//String for the name of the groundstate
	/*string GS_name="GS_";
	
	//Build the filename
	convert.str("");
	convert.clear();
	convert << N;
	GS_name = GS_name+"N"+convert.str();
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
	cout << "Time for writing GS: " << writetime_gs << " s" << endl;*/
	
	
	cout << "Succesfully took data for x=" << Hsngi.get_x() << endl << endl;
      }
  }
  energy_file.close();
  
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
    resfile << endl << scientific << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << x << "\t" << lam << scientific <<"\t" << E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" <<penalty_energy.re << "\t";
  }
  else
  {
    time(&start_local);
    resfile.open(FILENAME,ios_base::app);  
    resfile << "#N" << "\t" << "d" << "\t" << "D" << "\t" << "mu" << "\t" << "x" << "\t" << "e" << "\t" << "lambda" << "\t"<< "E0" << "\t" << "C0" << "\t" << "Cfrac" << "\t" << "E1" << "\t" << "C1" << "\tE_pen"<< endl;
    resfile << scientific << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << x << "\t" << lam << scientific << "\t"<< E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" << penalty_energy.re << "\t";
  }
  resfile.close();
  
  
  time(&end_local);
  writetime_data = difftime(end_local,start_local);
  
  cout << "Time for writing Data: " << writetime_data << " s" << endl;   
  
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
    benchmarking.open("Benchmark.txt",ios_base::app);
  else
  {
    benchmarking.open("Benchmark.txt");
    benchmarking << "#N\tD\tN_0\tT_total(s)\tT_perstep(s)" << endl;
  }
  benchmarking << std::resetiosflags(std::ios::scientific) << N << "\t" << D << "\t" << d << scientific << setprecision(10) <<"\t" << runtime_total << "\t" << runtime_total/num_of_steps << endl;
  benchmarking.close();
  
  //Free memory
  for(int k=0; k<Spin_MPOs.size(); k++)
  {
    Spin_MPOs[k]->clear();
    delete Spin_MPOs[k];
  }
  for(int k=0; k<Gauge_MPOs.size(); k++)
  {
    Gauge_MPOs[k]->clear();
    delete Gauge_MPOs[k];
  }
  
  mps->clear();
  delete mps;
  
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

int read_params(string paramfile,int &N, int &d, int &D, double &mu, double &x, double &lam, double &dt, int &steps)
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
    pfile >> dummy; pfile >> dummy; pfile >> x;
    pfile >> dummy; pfile >> dummy; pfile >> lam;
    pfile >> dummy; pfile >> dummy; pfile >> dt;
    pfile >> dummy; pfile >> dummy; pfile >> steps;
    
    cout << "Read N         = " << N << endl;
    cout << "Read N_0       = " << d<< endl;
    cout << "Read D         = " << D<< endl;
    cout << "Read mu        = " << mu<< endl;
    cout << "Read x        = " << x<< endl;
    cout << "Read lamda     = " << lam<< endl;
    cout << "Read dt        = " << dt<< endl;
    cout << "Read steps     = " << steps<< endl;
    
    pfile.close();
    return 0;
  }
}

int save_temp_result(MPS &mps,int N, int d, int D, double mu, double x, double lam, double dt, int steps)
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
  paramfile << "x = " << x << endl;
  paramfile << "lambda = " << lam << endl;
  paramfile << "dt = " << dt << endl;
  paramfile << "steps = " << steps << endl;
  paramfile << scientific << setprecision(11) ;
  paramfile << "Time: " << asctime(timeinfo) << endl;
  paramfile.close();
  
  mps.exportMPS("MPS.dat");
  
  cout << "Successfully saved a snapshot" << endl << endl;
  
  paramfile.close();
  
  return 0;
}

vector<complex_t> current_configuration(vector<MPO*> &Spin, vector<MPO*> &Gauge,MPS &mps,Contractor& contractor,int N)
{
  vector<complex_t> result;
  if(Spin.size()!=(Gauge.size()+1))
  {
    cout << "Error, cannot determine the configuration" << endl;
    exit(666);
  }
  else
  {
    result.resize(2*N-1);
    #pragma omp parallel for
    for(int k=1; k<=N; k++)
    {
      result[2*k-2]=contractor.contract(mps,*Spin[k-1],mps);
      if(k<N)
	result[2*k-1]=contractor.contract(mps,*Gauge[k-1],mps);
    }
  }
  return result;
}

int save_configuration(string filename, vector<complex_t> config)
{
  ofstream file;
  file.open(filename.c_str(),ios_base::app);
  for(int k=0; k<config.size()-1; k++)
    file << config[k] << "\t";
  file << config.back() << endl;
  file.close();
  return 0;
}



vector<complex_t> get_EntanglementEntropy(Contractor &contractor, MPS mps, int N)
{
  vector<complex_t> result;
  result.resize(2*N-2);

  #pragma omp parallel for 
  for (int i=1; i<2*N-1; i++)
  {
    result[i-1].re = contractor.getEntropy(mps,i);
    result[i-1].im = 0;
  }
  return result;
}

