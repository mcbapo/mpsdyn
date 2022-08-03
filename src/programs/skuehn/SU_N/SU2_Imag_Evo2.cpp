/**
   \file SU2_Imag_Evo.cpp
   Driver program for imaginary time evoltuion of the truncated SU(2) Hamiltonian with staggered fermions from Erez, Ignacio and Benni
   
  \param <N> (int) Number of Fermions
  \param <D> (int) Bond dimension
  \param <epsilon> (double) Hopping parameter
  \param <mu> (double) Mass
  \param <g> (double) Coupling constant
  \param <dt> (double) Time step size
  \param <name> (string) Name of a potential ground state file (optional)
  \param <leng> (int) length of the string, if 0 is given no string is placed
  
  \author Stefan Kühn
  \date 27/03/2015

*/


#include <cstdlib>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>
//Lib which is needed for sleep command on *nix Systems
//#include "unistd.h"

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include "Indices.h"
#include "misc.h"

#include "SU2Hamiltonian.h"


#define EPS 1E-10
#define NUM_OF_PARAMS 7
#define FILENAME "SU2.txt"
#define FILENAME_CONFIG "SU2_configuration.txt"
#define FILENAME_COEFFIS "SU2coefficients.txt"
#define CHARGE_FILE "SU2charge.txt"
#define EVOFILE "SU2efile.txt"
#define RESOLUTION 10
#define PROJECTION_INTERVAL 100
#define BACKUP_INTERVAL (100*RESOLUTION)
#define STRING_PARAM_FILE "string_params.dat"
#define TIMEFILE "time.txt"
#define PROJECTION_ON false

using namespace std;
using namespace shrt;

/**Struct representing a string */
struct str
{
  int pos;
  int leng;
};

/** Different methods to compute the string coefficients*/
enum StringCoefficients{
  OverlapSC,		//Overlap with a single strong coupling stringstate
  ReducedSC,		//Trace the region outside the string and compute local overlap
  SdaggerS,		//Expectation value of S^dagger S
  StringProjector,	//Project on states which have no flux on two sites and some finite flux in between
  ParticleProjector,	//Project on states which have a single fermion on two sites and two/zero fermions in between
  QsqCorrelation	//Determine the charge square correlations
};



/** Print MPS, in contrast to the regular streaming operator this shows explicitly the content of the matrices */
void print_MPS(MPS);

//Determine the current configuration
vector<complex_t> current_configuration(vector<MPO*> &Spin1,vector<MPO*> &Spin2,vector<MPO*> &Gauge,MPS &mps,Contractor& contractor,int N);

/** Save the current configuration to a given file */
int save_configuration(string filename, vector<complex_t> config);

/** Print the current configuration in a readable way */
void print_configuration(vector<complex_t> config);

/** Function to generate a backup of the current state during the time evolution*/
void generate_backup(const MPS &mps,int N,int D,double x,double mu, double g, double dt, int step, double curr_time);
int load_backup(MPS &mps,int &N,int &D,double &x,double &mu, double &g, double &dt, int &step, double &curr_time);

//Function to read the parameters for an initial superposition form a file
int read_string_parameter(string filename, bool &antiparticles, int &pos_left, int &pos_right, double &width_left, double &width_right, double &k_vec, bool &factor_on);

//Get full config of a state to check the Gauss Law manually afterwars
void get_full_config(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N);

//All kinds of functions to determine the string coefficients
void get_OverlapStrongCouplingStringstate(SU2Hamiltonian &H, vector<MPO*> &ops, int N, int &num_ops);
void get_SdaggerS(SU2Hamiltonian &H, vector<MPO*> &ops, int N, int &num_ops);
void get_ReducedOverlapStrongCouplingStringstates(SU2Hamiltonian &H,Contractor &contractor, vector<MPO*> &ops, int N, int &num_ops);
void get_LinkProjectorCoefficients(SU2Hamiltonian &H, vector<MPO*> &ops, int N, unsigned int method, int &num_ops);
void get_SpinProjectorCoefficients(SU2Hamiltonian &H, vector<MPO*> &ops, int N, int &num_ops);
void getChargeSquareCorrelations(SU2Hamiltonian &H, vector<MPO*> &ops, int N, int method, int &num_ops);
//Driver function to compute the string coefficients
vector<complex_t> get_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, vector<MPO*> ops, StringCoefficients method, MPS mps, int num_ops);
void set_up_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, vector<MPO*> &ops, int N, StringCoefficients method, int flag, int &num_ops);

//Function to read the parameters for the method to determine the string coefficients
int read_string_coefficients_parameter(string filename, StringCoefficients &method, int & flag);

//Conmpute the values for the charge square locally
vector <complex_t> get_ChargeSquare(Contractor &contractor, vector<MPO*> ops, MPS mps, int N);

//Simple factorial function to compute Taylor coefficients (as int might not be enough, I use double to prevent overflows)
double factorial(int n)
{
  double res=1.0;
  for(int i=n; i>0; i--)
    res *= (double) i;
  return res;
}

//Read all string configurations from a given file
void read_string_configs(string filename, vector<string_t> &configs);

int main(int argc,const char* argv[]){
  int cntr=0;
  
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage:   " << argv[0] << " N D µ g epsilon" << endl;
    cout << "N:       " << "Number of spins" << endl
	 << "D:       " << "Bond dimension of the MPS" << endl
	 << "epsilon: " << "Hopping parameter" << endl
	 << "µ:       " << "Mass" << endl
	 << "g:       " << "Coupling constant" << endl
	 << "n:       " << "Number of time steps to be computed" << endl
	 << "dt:      " << "Time step size" << endl
	 << "order:   " << "Order of the time evolution (optional)" << endl
	 << "leng:    " << "Length of the initial string (optional)" << endl
	 << "lambda:  " << "Strength of the penalty for the static charges (optional)" << endl
	 << "file:    " << "Name of the file for an existing state to be further evolved (optional)" << endl
	 << "offset:  " << "Offset to be used to evolve an existing state (optional)" << endl;
    return -1;
  }
  
  /***********************************************************************************
   * 					Variables				     *
   **********************************************************************************/
  
  //Number of spin sites
  int N=atoi(argv[++cntr]);
  //Bond dimension
  int D=atoi(argv[++cntr]);
  //Hopping parameter
  double epsilon=atof(argv[++cntr]);
  //Mass
  double mu=atof(argv[++cntr]);
  vector <double> masses(N,mu);
  //Coupling constant
  double g=atof(argv[++cntr]);
  //Number of time steps to be computed
  int steps=atoi(argv[++cntr]);
  //Time step size
  double dt=atof(argv[++cntr]);
  double curr_time=0.0;
  //Name of a possible state to import and go on involving
  string state_filename="";
  //Length of the initial string
  int leng=0;  
  //Strength of the penalty for the charges
  double lambda=0.0;
  //Variables for time measurement
  time_t start,end,start_local,end_local;  
  //File for results
  ofstream file;
  //Variables for Gauss Law Components
  complex_t Gx_val,Gy_val,Gz_val; 
  //Variable for the energy
  complex_t E0;
  //Vector to save a configuration
  vector<complex_t> configuration;
  //Vector to save the coefficients for the string states
  vector<complex_t> coefficients;
  //Variables for the error
  double curr_err,sum_of_errors=0;
  vector<double> errors;
  //Offset in case I want to continue evolving a state
  int offset=0;
  //Method to determine the coefficients
  StringCoefficients method=OverlapSC;
  int flag=0;
  vector<MPO*> string_observable_ops;
  int num_ops=0;
  //Variables for higher order Taylor expansion
  vector< const MPS* >  kets;
  vector< complex_t >  beta;
  vector< double >  normfactors;
  MPS *tmp_mps;
  int order=1;
  int backup_data_flag=-1;
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  contractor.setConvTol(1.0E-8);
  
  if(argc >= NUM_OF_PARAMS+2)
    order=atoi(argv[++cntr]);
  if(argc >= NUM_OF_PARAMS+3)
    leng=atoi(argv[++cntr]);
  if(argc >= NUM_OF_PARAMS+4)
    lambda=atof(argv[++cntr]);
  if(argc >= NUM_OF_PARAMS+5)
    state_filename=argv[++cntr];
  if(argc >= NUM_OF_PARAMS+6)
    offset=atoi(argv[++cntr]);
  
  //Look if there is a file describing which method with which flag should be used for determining the string coefficients
  if(file_exists("CoefficientMethod.txt"))
  {
    cout << "Reading parameters for string coefficient computation" << endl;
    read_string_coefficients_parameter("CoefficientMethod.txt",method,flag);
  }
  
  //Start time measurement
  time(&start);
  
  //Now manipulate the masses to have the starting and the end point more massive
  /*vector <int> chargepos;
  if(leng!=0)
  {
    int startpos = round(N/2)-round(leng/2);
    int endpos = startpos+leng;
    chargepos.push_back(startpos);
    chargepos.push_back(endpos);
    masses[startpos-1]=masses[startpos-1]*20;
    masses[endpos-1]=masses[endpos-1]*20;
    cout << "Masses: " << masses << endl;    
  }

  cout << "Lambda: " << lambda << endl;*/
  
  /***********************************************************************************
   * 				Set up MPOs and MPSs				     *
   **********************************************************************************/  

  //Get the Hamiltonian
  SU2Hamiltonian HNGI(N,epsilon,mu,g);
  //Get the Hamiltonian with site depended masses
  //SU2Hamiltonian HNGI(N,epsilon,masses,g);
  //Get the Hamiltonian with site depended masses and static charges imposed
  //SU2Hamiltonian HNGI(N,epsilon,masses,g,chargepos,lambda);
    
  //Turn gauge interaction off in case this is wanted
#ifdef GAUGEOFF
  HNGI.removeGaugeInteraction();
  HNGI.updateHMPO();
#endif
  
  const MPO& hamil=HNGI.getHMPO();
  
  //Prepare the single body MPOs
  vector<MPO*> Spin1_MPOs,Spin2_MPOs,Gauge_MPOs;  
  MPO *mpo;
  for(int k=1; k<=N; k++)
  {
    mpo=new MPO(1);
    HNGI.getSpin1MPO(*mpo,k);
    Spin1_MPOs.push_back(mpo);    
    mpo=new MPO(1);
    HNGI.getSpin2MPO(*mpo,k);
    Spin2_MPOs.push_back(mpo);       
    if(k<N)
    {
      mpo=new MPO(1);
      HNGI.getJsquaredMPO(*mpo,k);
      Gauge_MPOs.push_back(mpo);       
    }
  }  
  
  vector<MPO*> charge;
  for(int k=1; k<=N; k++)
  {
    mpo=new MPO(1);
    HNGI.getChargeSquareMPO(*mpo,k);
    charge.push_back(mpo);
  }
  
  //Set up the MPOs for the string coefficients
  set_up_StringCoefficients(HNGI,contractor,string_observable_ops,N,method,flag,num_ops);
  
  //Get the time evolution MPO
  MPO Umpo(1);
  HNGI.constructEvolutionMPO(Umpo,dt,true);
  //HNGI.constructEvolutionMPOSiteDepMass(Umpo,dt,true);
  //HNGI.constructEvolutionMPOSiteDepMassCharge(Umpo,dt,chargepos,lambda,true);
  
  //Get the projection operators for projecting into the gauss law fulfilling subspace
  MPO Podd(1),Peven(1);
  //HNGI.constructProjectors(Podd,Peven);
  //HNGI.constructProjectors(Podd,Peven,-1,leng,y_comp);
  
  //MPS for the state and auxiliary state and the interacting vacuum
  MPS mps,aux;
  
  //MPOs to monitor the Gauss Law
  MPO Gx(1), Gy(1), Gz(1);
  HNGI.getGaussLawMPO(Gx,x_comp);
  HNGI.getGaussLawMPO(Gy,y_comp);
  HNGI.getGaussLawMPO(Gz,z_comp);
  
  //MPO to monitor the charge
  MPO charge_square(1);  
  HNGI.getChargeSquareMPO(charge_square);
  
  //Take care of the initial MPS to start the evolution
  if(!state_filename.empty())
  {
    //Case that I want to import a state
    if(file_exists(state_filename))
      mps.importMPS(state_filename.c_str());
    else
    {
      cout << "Warning, could not import state" << endl;
      exit(666);
    }
  }
  else if (file_exists("BackupData.txt"))
  {
    //I have a backup, so read the parameters from the backup
    backup_data_flag=load_backup(mps,N,D,epsilon,mu,g,dt,offset,curr_time);
  }
  else if(leng!=0)
  {
    /*******************************************************
     * Generate a string (on top of the interacting vacuum)
     * ****************************************************/
     /*MPO StringOp(1);
     MPS aux;    
     HNGI.getStringMPO2(StringOp,leng,0,false);
     cout << "Creating a string with length " << leng << endl;
     aux.approximate(StringOp,vacuum,D);    
     contractor.optimize(StringOp,vacuum,aux,D);
     aux.gaugeCond('R',true);
     mps=MPS(aux);*/
     /*******************************************************
     * Generate a  strong coupling string
     * ****************************************************/
     /*MPO StringOp(1);
     MPS aux,tmp;    
     HNGI.getStringMPO2(StringOp,leng,0,false);
     cout << "Creating a string with length " << leng << endl;
     HNGI.constructInitialMPS(tmp);
     aux.approximate(StringOp,tmp,D);    
     contractor.optimize(StringOp,tmp,aux,D);
     aux.gaugeCond('R',true);
     mps=MPS(aux);*/
    /*******************************************************
     * Put some momentum in the (interacting) vacuum
     * ****************************************************/
    /*MPO Moving(1);    
    cout << "Generating moving quarks" << endl;*/
    //First option (only distinguish quarks and antiquarks)
    /*HNGI.getMovingQuarks(Moving,1.0);
    aux.approximate(Moving,vacuum,D);    
    contractor.optimize(Moving,vacuum,aux,D);
    aux.gaugeCond('R',true);
    mps=MPS(aux);*/
    
    /*******************************************************
     * Get a large superposition of STRONG COUPLING strings peaked around predifined values
     * ****************************************************/
    /*if(!file_exists(STRING_PARAM_FILE))
      HNGI.getMovingString(mps,false,N/2,N/2+1,2.0,2.0,1.0);
    else
    {
      int pos_left,pos_right;
      double width_left,width_right,k_vec;
      bool antiparticles,factor_on;
      read_string_parameter(STRING_PARAM_FILE,antiparticles,pos_left,pos_right,width_left,width_right,k_vec,factor_on);
      HNGI.getMovingString(mps,antiparticles,pos_left,pos_right,width_left,width_right,k_vec,factor_on);
    }*/
    
    /*******************************************************
     * Get a large superposition of strings peaked around predifined values
     * ****************************************************/
    /*MPO StringOp(1);
    MPS aux;
    if(!file_exists(STRING_PARAM_FILE))
      HNGI.getMovingString(StringOp,N/2,N/2+1,2.0,2.0,1.0);
    else
    {
      int pos_left,pos_right;
      double width_left,width_right,k_vec;
      bool antiparticles,factor_on;
      read_string_parameter(STRING_PARAM_FILE,antiparticles,pos_left,pos_right,width_left,width_right,k_vec,factor_on);
      HNGI.getMovingString(StringOp,pos_left,pos_right,width_left,width_right,k_vec,factor_on);
    }    
    aux.approximate(StringOp,vacuum,D);    
    contractor.optimize(StringOp,vacuum,aux,D);
    aux.gaugeCond('R',true);
    mps=MPS(aux);*/
    
    
    
    /*******************************************************
     * Non gauge invariant superposition for testing if the particles move
     * ****************************************************/    
    /*cout << "Creating a non gauge invariant superposition" << endl;
    cout << "WARNING: Projector should be switched off for these runs!" << endl;
    if(!file_exists(STRING_PARAM_FILE))
      HNGI.getMovingParticles(mps,(int) N/4,(int) N/4,2.0,2.0,1.0);
    else
    {
      int pos_left,pos_right;
      double width_left,width_right,k_vec;
      bool antiparticles,factor_on;
      
      read_string_parameter(STRING_PARAM_FILE,antiparticles,pos_left,pos_right,width_left,width_right,k_vec,factor_on);
      HNGI.getMovingParticles(mps,pos_left,pos_right,width_left,width_right,k_vec);
    }*/
    /*******************************************************
     * Drop a single charge on the interacting vacuum (non gauge invariant)
     * ****************************************************/    
    /*cout << "Creating a (non gauge invariant) single charge" << endl;
    cout << "WARNING: Projector should be switched off for these runs!" << endl;
    MPS aux;
    MPO SingleChargeMPO(1);
    HNGI.createSingleCharge(SingleChargeMPO,(int)(N/2),1);
    aux.approximate(SingleChargeMPO,vacuum,D);    
    contractor.optimize(SingleChargeMPO,vacuum,aux,D);
    aux.gaugeCond('R',true);
    mps=MPS(aux);*/
    /*******************************************************
     * Create a "heavy string" state as initial state for imaginary time evolution
     * ****************************************************/ 
    HNGI.getHeavyString(mps,leng);
    get_full_config(HNGI,contractor,mps,N);
    mps.increaseBondDimension(D);
    /*******************************************************
     * Create a strong coupling state containing multiple strings as initial state for imaginary time evolution
     * ****************************************************/ 
    /*vector<string_t> configs;  
    read_string_configs("strconfigs.dat",configs);  
    HNGI.constructStringMPS(mps,configs);
    get_full_config(HNGI,contractor,mps,N);
    mps.increaseBondDimension(D);*/
  }
  else
  {
    //I simply start with the strong coupling state
    HNGI.constructInitialMPS(mps);
  } 
  
  //Determine initial values of energy and Gauss Law
  E0=contractor.contract(mps,hamil,mps);   
  Gx_val=contractor.contract(mps,Gx,mps);
  Gy_val=contractor.contract(mps,Gy,mps);
  Gz_val=contractor.contract(mps,Gz,mps);
  
  cout << "Inital energy: " << contractor.contract(mps,hamil,mps) << endl;
  cout << "Gx: " << Gx_val << endl;
  cout << "Gy: " << Gy_val << endl;
  cout << "Gz: " << Gz_val << endl << endl;   
  
  if(backup_data_flag != 0)
  {
    //Save configuration once before evolution
    configuration=current_configuration(Spin1_MPOs,Spin2_MPOs,Gauge_MPOs,mps,contractor,N);
    save_configuration(FILENAME_CONFIG,configuration);
  
    //Print start configuration
    cout << "Start configuration:" << endl;
    print_configuration(configuration);  
    
    //Determine the charge
    coefficients=get_ChargeSquare(contractor,charge,mps,N);
    file.open(CHARGE_FILE, ios_base::app);
    file << 0 << "\t" << 0 ;
    for(int ind=0; ind<coefficients.size(); ind++)
      file << "\t" << coefficients[ind];
    file << endl;
    file.close(); 
    
    //Compute some string observables
    coefficients=get_StringCoefficients(HNGI,contractor,string_observable_ops,method,mps,num_ops);
    file.open(FILENAME_COEFFIS, ios_base::app);
    file << 0 << "\t" << 0;
    for(int ind=0; ind<coefficients.size(); ind++)
    {
      file << "\t" << coefficients[ind];
    }
    file << endl;
    file.close();
  }    
    
  //Save initial state of the time evolution
  MPS init_state(mps); 
    
  //Copy mps to have an MPS with right length and physical dimension, then increase the bond dimension and set it to a random state
  aux=MPS(mps);
  aux.increaseBondDimension(D);
  aux.setRandomState();  
  
  //Prepare the evolution file
  if(!file_exists(EVOFILE))
  {
    file.open(EVOFILE, ios_base::app);
    file << "#Step\tTime\tBond Dim\tepsilon\tmu\tg\tE0\tGx\tGy\tGz\tOverlap\tCharge square\tError" << endl;
    file.close();
  }
  //Save once before evolution in case I do not further evolve a backup
  if(backup_data_flag!=0)
  {
   file.open(EVOFILE, ios_base::app);
   file << 0 << "\t" << 0 << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << setprecision(10) << E0 << "\t" << Gx_val << "\t" << Gy_val << "\t" << Gz_val  << "\t" << abs(contractor.contract(init_state,mps)) << "\t" << contractor.contract(mps,charge_square,mps) << "\t" << sum_of_errors << endl;
   file.close();
  }
  
  /***********************************************************************************
   *				Actual evolution				     *
   **********************************************************************************/
  //Monitor the energy during evolution and decide whether it's time to stop or to reduce the time step or something else
  complex_t energy_old=E0;
  
  //Start time measurement
  time(&start_local);
  
  cout << endl;
  cout << "------------------------------------------------------------------------"<<endl;
  cout << "              Starting time evolution of order "<< order <<endl;
  cout << "------------------------------------------------------------------------"<<endl;
  for(int i=1+offset; i<=steps; i++)
  {
#ifdef RAMP
    if(((double) i)/2000.0<=1.0)
      //HNGI.set_epsilon(epsilon*(((double) i)/2000.0));
      HNGI.set_epsilon(epsilon*pow(((double) i)/2000.0,3.0));
    else
      HNGI.set_epsilon(epsilon);
    HNGI.constructEvolutionMPO(Umpo,dt);
#endif
    if(order==1)
    {
      ///////////////////////////////////////////////////////////
      //Simple first order Taylor expansion, not unitary anymore, treated separately, because I have an explicit MPO for the evolution operator which is more efficient
      ///////////////////////////////////////////////////////////
      curr_err=-1.0;
      contractor.optimize(Umpo,mps,aux,D,&curr_err);
      errors.push_back(curr_err);
      //As the Taylor evolution is not unitary, restore the norm
      //aux.setNormFact(1.0);
      aux.gaugeCond('R',true);
    }
    else
    {
      ///////////////////////////////////////////////////////////
      //Taylor expansion up to a given order
      ///////////////////////////////////////////////////////////
      tmp_mps=new MPS(mps);
      kets.push_back(tmp_mps);
      normfactors.push_back(1.0);
      beta.push_back(ONE_c);
      for(int k=1; k<=order; k++)
      {
	curr_err=-1.0;
	aux.approximate(hamil,mps,D);
	contractor.optimize(hamil,mps,aux,D,&curr_err);
	errors.push_back(curr_err);
	//contractor.optimize(hamil,mps,aux,D);
	//I save the norm of the state separately, this seems to make the approximation more stable
	normfactors.push_back(normfactors[k-1]*sqrt(contractor.contract(aux,aux).re));
	//Put the norm in the prefactors of the Taylor expansion
	beta.push_back(normfactors[k]*pow(-1.0*dt,k)/factorial(k)*ONE_c);
	aux.setNormFact(1.0);
	tmp_mps=new MPS(aux);
	kets.push_back(tmp_mps);
	mps=aux;
	aux.setRandomState();
      }    
      //Compute the resulting state and normalize it
      contractor.optimizeSum(kets,beta,aux,D);
      aux.gaugeCond('R',true);
      //Free memory
      for (int k=0; k<kets.size(); k++)
	delete kets[k];
      kets.clear();
      normfactors.clear();
      beta.clear();
    }   
    
    if(PROJECTION_ON && ((i%PROJECTION_INTERVAL)==0))
    {    
      cout << "Projecting" << endl;
      mps.approximate(Podd,aux,D);
      //WARNING: Error computation seems to be very expensive and memory consuming and might not be possible in larger simulations
      curr_err=-1.0;
      //contractor.optimize(Podd,aux,mps,D,&curr_err);
      contractor.optimize(Podd,aux,mps,D);
      sum_of_errors=curr_err;
      mps.gaugeCond('R',true); 
      aux.approximate(Podd,mps,D);
      //WARNING: Error computation seems to be very expensive and memory consuming and might not be possible in larger simulations
      curr_err=-1.0;
      //contractor.optimize(Peven,mps,aux,D,&curr_err);
      contractor.optimize(Peven,mps,aux,D);
      sum_of_errors+=curr_err;
      aux.gaugeCond('R',true);
      cout << "Errors projection: " << sum_of_errors << endl;
    }
    else if((i%PROJECTION_INTERVAL)==0)
      cout << "WARNING: Projector switched off!" << endl;
    
    //Time step complete, update variable for total time
    curr_time += dt;
    
    //Take some data
    if((i%RESOLUTION)==0)
    {
#ifdef RAMP
      cout << "Current value of epsilon: "<< HNGI.get_epsilon() << endl;
      HNGI.updateHMPO();
#endif
      E0=contractor.contract(aux,hamil,aux);
      Gx_val=contractor.contract(aux,Gx,aux);
      Gy_val=contractor.contract(aux,Gy,aux);
      Gz_val=contractor.contract(aux,Gz,aux);
      cout << "Step: " << i << ", Norm: " << contractor.contract(aux,aux) << ", Energy: " << E0 << ", Gx: " << Gx_val << ", Gy: " << Gy_val << ", Gz: " << Gz_val <<endl;
      configuration=current_configuration(Spin1_MPOs,Spin2_MPOs,Gauge_MPOs,aux,contractor,N);
      save_configuration(FILENAME_CONFIG,configuration);  
      
      sum_of_errors =std::accumulate(errors.begin(),errors.end(),0.0);
      cout << "Errors: " << sum_of_errors << endl;
      errors.clear();      
      
      file.open(EVOFILE, ios_base::app);
      file << i << "\t" << curr_time << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << setprecision(10) << E0 << "\t" << Gx_val << "\t" << Gy_val << "\t" << Gz_val  << "\t" << abs(contractor.contract(init_state,aux)) << "\t" << contractor.contract(aux,charge_square,aux) << "\t" << sum_of_errors << endl;
      file.close();      
      
      if(N<=6)
	print_configuration(configuration);
      if((i%BACKUP_INTERVAL)==0 && i!=0)
	generate_backup(aux,N,D,epsilon,mu,g,dt,i,curr_time);
      
      //Determine the charge
      coefficients=get_ChargeSquare(contractor,charge,aux,N);
      file.open(CHARGE_FILE, ios_base::app);
      file << i << "\t" << curr_time ;
      for(int ind=0; ind<coefficients.size(); ind++)
	file << "\t" << coefficients[ind];
      file << endl;
      file.close(); 
      
      //Compare to old energy value and act accordingly
//       if((abs((E0-energy_old)/energy_old)<1.0E-4) && (dt>=1.0E-3))
//       {
// 	dt /=4.0;
// 	HNGI.constructEvolutionMPO(Umpo,dt);
// 	cout << "Relative change: " << abs((E0-energy_old)/energy_old) <<", reducing time step, dt=" << dt << endl;	
// 	energy_old=E0;
//       }
      
      /*if(i==40)
	dt /=10.0;
      
      if(i==200)
	dt /=10.0;*/
      
      //Compute some string observables
      /*coefficients=get_StringCoefficients(HNGI,contractor,string_observable_ops,method,aux,num_ops);
      file.open(FILENAME_COEFFIS, ios_base::app);
      file << i << "\t" << curr_time ;
      for(int ind=0; ind<coefficients.size(); ind++)
      {
	file << "\t" << coefficients[ind];
      }
      file << endl;
      file.close();*/
    }
    mps=aux;
  }  
  
  //Get final results
  E0=contractor.contract(aux,hamil,aux);
  Gx_val=contractor.contract(aux,Gx,aux);
  Gy_val=contractor.contract(aux,Gy,aux);
  Gz_val=contractor.contract(aux,Gz,aux);
  
  time(&end_local);
  file.open(TIMEFILE,ios_base::app);
  file << "Time per step\t" << difftime(end_local,start_local)/(steps-offset) << "s" << endl;
  file.close();
  
  
  //Save the ground state
  string GS_name="GS_";
  //String needed for the conversion
  ostringstream convert;  
  //Build the filename
  convert << N;
  GS_name = GS_name+"N"+convert.str();
  convert.str("");
  convert.clear();
  convert << D;
  GS_name = GS_name+"D"+convert.str()+".dat";
  mps.exportMPS(GS_name.c_str());
  
  if(file_exists(FILENAME))
  {
    file.open(FILENAME, ios_base::app);
    file << setprecision(16) << N << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << E0 << "\t" << Gx_val << "\t" << Gy_val << "\t" << Gz_val << endl;
  }
  else
  {
    file.open(FILENAME, ios_base::app);
    file << setprecision(16) << "#N\tD\tepsilon\tmu\tg\tE0\tGx\tGy\tGz" << endl;
    file << setprecision(16) << N << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << E0 << "\t" << Gx_val << "\t" << Gy_val << "\t" << Gz_val << endl;
  }
  file.close();
  
  //Save the Gauss Law configuration to be able to compare to the inital one
  get_full_config(HNGI,contractor,mps,N);
  
  
  //Free memory
  for(int k=0; k<Spin1_MPOs.size(); k++)
  {
    Spin1_MPOs[k]->clear();
    delete Spin1_MPOs[k];
  }
  for(int k=0; k<Spin2_MPOs.size(); k++)
  {
    Spin2_MPOs[k]->clear();
    delete Spin2_MPOs[k];
  }
  for(int k=0; k<Gauge_MPOs.size(); k++)
  {
    Gauge_MPOs[k]->clear();
    delete Gauge_MPOs[k];
  }
  
  for(int k=0; k<charge.size(); k++)
  {
    charge[k]->clear();
    delete charge[k];
  } 
  
  for(int k=0; k<string_observable_ops.size(); k++)
  {
    string_observable_ops[k]->clear();
    delete string_observable_ops[k];
  }
  
  //End time measurement and output the runtime
  time(&end);
  double runtime_total = difftime(end,start);
  cout << "Runtime: " << runtime_total << "s=" << runtime_total/60.0 << "min="<< runtime_total/3600.0 << "h="<<runtime_total/(3600.0*24.0) << "d"<< endl;
  

  cout << "End of program" << endl;
  return 0; 
}


/***********************************************************************************
 *			Previously declared functions				   *
 **********************************************************************************/

//Explicitly output the MPS meaning all the matrices
void print_MPS(MPS mps)
{
  cout << "MPS of length " << mps. getLength() << endl;
  for(int i=0; i<mps.getLength(); i++)
  {
    cout << "[" << i << "]=" <<  (mps.getA(i)).getA() << endl;
  }
}

//Get the spin and flux configuration of a given MPS on every site
vector<complex_t> current_configuration(vector<MPO*> &Spin1, vector<MPO*> &Spin2, vector<MPO*> &Gauge,MPS &mps,Contractor& contractor,int N)
{
  vector<complex_t> result;
  if((Spin1.size()!=Spin2.size()) || (Spin1.size()!=(Gauge.size()+1)))
  {
    cout << "Error, cannot determine the configuration" << endl;
    exit(666);
  }
  else
  {
    result.resize(3*N-1);
    #pragma omp parallel for
    for(int k=1; k<=N; k++)
    {
      result[3*k-3]=contractor.contract(mps,*Spin1[k-1],mps);
      result[3*k-2]=contractor.contract(mps,*Spin2[k-1],mps);
      if(k<N)
	result[3*k-1]=contractor.contract(mps,*Gauge[k-1],mps);
    }
  }
  return result;
}


//Once a configuration vector has been obtained with the previous function, this function saves it
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
 
//Print a configuration vector on screen (resonably formatted)
void print_configuration(vector<complex_t> config)
{
  for(int k=0; k<config.size(); k+=3)
  {
    printf("%.2f ",config[k].re);
    printf("%.2f ",config[k+1].re);
    if(k+2<config.size())
      printf("%.2f | ",config[k+2].re);
  }
  printf("\n");
}

void get_OverlapStrongCouplingStringstate(SU2Hamiltonian &H, vector<MPO*> &ops, int N, int &num_ops)
{
  MPO* mpo;
  MPS strong_coupling;
  
  cout << "Using overlap with the strong coupling string state" << endl;
  
  //I want the strong coupling vacuum not the interacting vacuum I get as input
  H.constructInitialMPS(strong_coupling);
  
  num_ops=0;
  for(int leng=1; leng<=N; leng+=2)
  {
    for(int pos=1; pos<=N-leng; pos++)
    {
      mpo = new MPO(1);
      H.getStringMPO2(*mpo,leng,pos);
      //WARNING: contract is taking arguments in this order: Ket, MPO, Bra !!!
      ops.push_back(mpo);
      num_ops++;
    }
  }
}

void get_SdaggerS(SU2Hamiltonian &H, vector<MPO*> &ops, int N, int &num_ops)
{
  MPO* mpo;
  vector<complex_t> result;
  
  cout << "Using expectation value of S^dagger S" << endl;

  num_ops=0;
  //Loop over length
  for(int leng=1; leng<=N-1; leng+=2)
  {
    //Loop over position
    for(int pos=1; pos<=N-leng; pos++)
    {
      /*******************************************************************************
       * Compute <\Psi(t)|SS^\dagger|\Psi>
       * ****************************************************************************/      
      mpo = new MPO(1);
      //Get the adjoint string operator S^\dagger
      H.getStringMPO2(*mpo,leng,pos,true);      
      //result.push_back(contractor.contract2(mpo,mps));
      ops.push_back(mpo);
      num_ops++;
    }
  }
}

void get_ReducedOverlapStrongCouplingStringstates(SU2Hamiltonian &H, Contractor & contractor,vector<MPO*> &ops, int N, int &num_ops)
{
  MPO *mpo, StringOp(1);
  MPS strong_coupling_vacuum,op;
  
  H.constructInitialMPS(strong_coupling_vacuum);
  op=MPS(strong_coupling_vacuum);
  op.increaseBondDimension(2);
  op.setRandomState();
  
  cout << "Using reduced overlap with the strong coupling string state" << endl;
  
  num_ops=0;
  //Loop over length
  for(int leng=1; leng<=N-1; leng+=2)
  {
    //Loop over position
    for(int pos=1; pos<=N-leng; pos++)
    {
      /*******************************************************************************
       * Compute <\Psi(t)|Id_A \otimes Tr_A(S|0><0|S^\dagger) |\Psi(t)>
       * ****************************************************************************/ 
      mpo = new MPO(1);
      H.getStringMPO2(StringOp,leng,pos);
      op.approximate(StringOp,strong_coupling_vacuum,2);    
      contractor.optimize(StringOp,strong_coupling_vacuum,op,2);
      H.buildReducedStringMPO(*mpo,op,pos,pos+leng);
      ops.push_back(mpo);      
      num_ops++;
    }
  }
}

void get_LinkProjectorCoefficients(SU2Hamiltonian &H, vector<MPO*> &ops, int N, unsigned int method, int &num_ops)
{  
  MPO* Projector;
  
  if(method==0)
    cout << "Using string projector without taking the spins into account" << endl;
  else if(method==1)
    cout << "Using string projector and looking for parallel polarizations of the spins" << endl;
  if(method==2)
    cout << "Using string projector and looking for polarization in up (down) direction on odd (even) sites" << endl;  
 
  num_ops=0;
  for(int leng=1; leng<=N-1; leng++)
  {
    //Loop over position
    for(int pos=1; pos<=N-leng; pos++)
    {
      /*******************************************************************************
       * Compute <\Psi(t)|P|\Psi(t)>
       * ****************************************************************************/
      Projector=new MPO(1);
      H.getStringProjector(*Projector,leng,pos,method);
      ops.push_back(Projector);
      num_ops++;
    }
  }
}

void get_SpinProjectorCoefficients(SU2Hamiltonian &H, vector<MPO*> & ops, int N, int &num_ops)
{
  MPO* Projector;
  
  cout << "Using projector on spin sites" << endl;  
 
  num_ops=0;
  for(int leng=1; leng<=N-1; leng++)
  {
    //Loop over position
    for(int pos=1; pos<=N-leng; pos++)
    {
      /*******************************************************************************
       * Compute <\Psi(t)|P|\Psi(t)>
       * ****************************************************************************/
      Projector = new MPO(1);
      H.getCountStringsMPO(*Projector,leng,pos);
      ops.push_back(Projector);
      num_ops++;
    }
  }
}


void getChargeSquareCorrelations(SU2Hamiltonian &H, vector<MPO*> &ops, int N, int method, int &num_ops)
{
  MPO* correlation;
  
  if(method==0)
    cout << "Using old version with projector" << endl;
  else if(method==1)
    cout << "Using old version without projector" << endl;
  else
    cout << "Using new version with/without unitaries" << endl;
  
  num_ops=0;
  //Loop over length
  for(int leng=1; leng<=N-1; leng++)
  {
    //Loop over position
    for(int pos=1; pos<=N-leng; pos++)
    {
      correlation = new MPO(1);
      if(method==0)
	H.getChargeCorrelatorMPO(*correlation,leng,pos,true);
      else if(method==1)
	H.getChargeCorrelatorMPO(*correlation,leng,pos,false);
      else
	H.getChargeCorrelatorMPO2(*correlation,leng,pos,0.5);
      ops.push_back(correlation);
      num_ops++;
    }
  }
}

//Function to compute the coefficients (basically does nothing but choosing the right subroutine)
void set_up_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, vector<MPO*> &ops, int N, StringCoefficients method, int flag, int &num_ops)
{  
  if(method==OverlapSC)
   get_OverlapStrongCouplingStringstate(H,ops,N,num_ops);
  else if(method==ReducedSC)
    get_ReducedOverlapStrongCouplingStringstates(H,contractor,ops,N,num_ops);
  else if(method==SdaggerS)
    get_SdaggerS(H,ops,N,num_ops);
  else if(method==StringProjector)
    get_LinkProjectorCoefficients(H,ops,N,flag,num_ops);
  else if(method==ParticleProjector)
    get_SpinProjectorCoefficients(H,ops,N,num_ops);
  else if(method==QsqCorrelation)
    getChargeSquareCorrelations(H,ops,N,flag,num_ops);
}

vector<complex_t> get_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, vector<MPO*> ops, StringCoefficients method, MPS mps, int num_ops)
{
  vector<complex_t> result;
  result.resize(num_ops);
  
  if(method==OverlapSC)
  {
    MPS strong_coupling;
    H.constructInitialMPS(strong_coupling);
    #pragma omp parallel for 
    for (int i=0; i<num_ops; i++)
      result[i] = contractor.contract(strong_coupling,*ops[i],mps)/sqrt(abs(contractor.contract2(*ops[i],strong_coupling)));
  }
  else if(method==SdaggerS)
  {
    #pragma omp parallel for 
    for (int i=0; i<num_ops; i++)
      result[i] = contractor.contract2(*ops[i],mps);
  }
  else 
  {
    #pragma omp parallel for 
    for (int i=0; i<num_ops; i++)
      result[i] = contractor.contract(mps,*ops[i],mps);
  }
  return result;
}

//Determine all three components of L, R and Q on every site and check the GaussLaw
void get_full_config(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N)
{
  vector <complex_t> Lx,Ly,Lz,Rx,Ry,Rz,Qx,Qy,Qz,Gx,Gy,Gz;
  MPO mpo(1);
  
  ofstream file;
  
  for(int i=1; i<=N; i++)
  {
    if(i<N)
    {
      H.getLMPO(mpo,i,x_comp);
      Lx.push_back(contractor.contract(mps,mpo,mps));
      H.getLMPO(mpo,i,y_comp);
      Ly.push_back(contractor.contract(mps,mpo,mps));
      H.getLMPO(mpo,i,z_comp);
      Lz.push_back(contractor.contract(mps,mpo,mps));
      
      H.getRMPO(mpo,i,x_comp);
      Rx.push_back(contractor.contract(mps,mpo,mps));
      H.getRMPO(mpo,i,y_comp);
      Ry.push_back(contractor.contract(mps,mpo,mps));
      H.getRMPO(mpo,i,z_comp);
      Rz.push_back(contractor.contract(mps,mpo,mps));
    }
    H.getChargeMPO(mpo,i,x_comp);
    Qx.push_back(contractor.contract(mps,mpo,mps));
    H.getChargeMPO(mpo,i,y_comp);
    Qy.push_back(contractor.contract(mps,mpo,mps));
    H.getChargeMPO(mpo,i,z_comp);
    Qz.push_back(contractor.contract(mps,mpo,mps));
  }    
  for (int i=1;i<=N; i++)
  {
    if(i==1)
    {
      Gx.push_back(Lx[i-1]-Qx[i-1]);
      Gy.push_back(Ly[i-1]-Qy[i-1]);
      Gz.push_back(Lz[i-1]-Qz[i-1]);
    }
    else if(i==N)
    {
      Gx.push_back(-Rx[i-2]-Qx[i-1]);
      Gy.push_back(-Ry[i-2]-Qy[i-1]);
      Gz.push_back(-Rz[i-2]-Qz[i-1]);
    }
    else
    {
      Gx.push_back(Lx[i-1]-Rx[i-2]-Qx[i-1]);
      Gy.push_back(Ly[i-1]-Ry[i-2]-Qy[i-1]);
      Gz.push_back(Lz[i-1]-Rz[i-2]-Qz[i-1]);
    }
  }
  
  if(file_exists("GaussLaw.txt"))
    file.open("GaussLaw.txt", ios_base::app);
  else
  {
    file.open("GaussLaw.txt", ios_base::app);
    file << "#N\tLx\tLy\tLz\tRx\tRy\tRz\tQx\tQy\tQz\tGx\tGy\tGz" << endl;
  }
  for(int i=1; i<=N; i++)
  {
    if(i<N)
      file << setprecision(10)<< i << "\t" << Lx[i-1] << "\t" << Ly[i-1] << "\t" << Lz[i-1] << "\t" << Rx[i-1] << "\t" << Ry[i-1] << "\t" << Rz[i-1] << "\t" << Qx[i-1] << "\t" << Qy[i-1] << "\t" << Qz[i-1] << "\t" << Gx[i-1] << "\t" << Gy[i-1] << "\t" << Gz[i-1] << endl;
    else
      file << setprecision(10) << i << "\t" << "0" << "\t" << "0" << "\t" << "0" << "\t" << "0" << "\t" << "0" << "\t" << "0" << "\t" << Qx[i-1] << "\t" << Qy[i-1] << "\t" << Qz[i-1] << "\t" << Gx[i-1] << "\t" << Gy[i-1] << "\t" << Gz[i-1] << endl;
  }
  file.close();
}

vector <complex_t> get_ChargeSquare(Contractor &contractor, vector<MPO*> ops, MPS mps, int N)
{
  vector<complex_t> result;
  result.resize(N);
  
#pragma omp parallel for
  for(int i=0; i<N; i++)
    result[i]=contractor.contract(mps,*ops[i],mps);
  
  return result;
}


//Save an MPS and the data belonging to it
void generate_backup(const MPS &mps,int N,int D,double x,double mu, double g, double dt, int step, double curr_time)
{
  ofstream file;
  time_t timer;
  char filename[50];
  
  //Set the timer to see when backup has been generated
  time(&timer);
  
  file.open("BackupData.txt");
  file << asctime(localtime(&timer))
       << "N:    " << N << endl
       << "D:    " << D << endl
       << "x:    " << x << endl
       << "mu:   " << mu << endl
       << "g:    " << g << endl
       << "dt:   " << dt << endl
       << "step: " << step << endl
       << "time: " << curr_time << endl;  
  file.close();
  
  mps.exportMPS("MPS.dat");
}

int load_backup(MPS &mps,int &N,int &D,double &x,double &mu, double &g, double &dt, int &step, double &curr_time)
{
  ifstream file;
  string variable_name;
  
  //Check if BackupData is there, in case yes proceed with loading the previous snapshot
  if(!file_exists("BackupData.txt"))
    return -1;
  //In case there is some backupdata, see if the MPS file exists
  if(!file_exists("MPS.dat"))
  {
    cout << "Error while loading backup, file for MPS is missing" << endl;
    return -1;
  }
  //Everything seems to be there, now proceed with loading the file  
  file.open("BackupData.txt");  
  //File is not working
  if(!file.good())
    return -1;
  //Start reading the file
  //file >> variable_name;
  getline(file,variable_name);
  cout << "Loading backup from the " << variable_name << endl;
  if(!file.good())
    return -1;
  file >> variable_name >> N;
  cout << "Read: " << variable_name << "\t" << N << endl;  
  if(!file.good())
    return -1;
  file >> variable_name >> D;
  cout << "Read: " << variable_name << "\t" << D << endl;  
  if(!file.good())
    return -1;
  file >> variable_name >> x;
  cout << "Read: " << variable_name << "\t" << x << endl;  
  if(!file.good())
    return -1;
  file >> variable_name >> mu;
  cout << "Read: " << variable_name << "\t" << mu << endl;  
  if(!file.good())
    return -1;
  file >> variable_name >> g;
  cout << "Read: " << variable_name << "\t" << g << endl;  
  if(!file.good())
    return -1;
  file >> variable_name >> dt;
  cout << "Read: " << variable_name << "\t" << dt << endl;  
  if(!file.good())
    return -1;
  file >> variable_name >> step;
  cout << "Read: " << variable_name << "\t" << step << endl;   
  if(!file.good())
    return -1;
  file >> variable_name >> curr_time;
  cout << "Read: " << variable_name << "\t" << curr_time << endl;   
  
  mps.importMPS("MPS.dat");
  return 0;
}

//In case I want to create one of the large superpositions of string states, this allows to read the desired parameters form an external file
int read_string_parameter(string filename, bool &antiparticles, int &pos_left, int &pos_right, double &width_left, double &width_right, double &k_vec, bool &factor_on)
{
  ifstream file;
  string variable_name;
  file.open(filename.c_str());
  
  if(!file.good())
    return -1;
  
  file >> variable_name >> antiparticles;
  cout << "Read: " << variable_name << "\t" << antiparticles << endl;  
  file >> variable_name >> pos_left;
  cout << "Read: " << variable_name << "\t" << pos_left << endl;  
  file >> variable_name >> pos_right;
  cout << "Read: " << variable_name << "\t" << pos_right << endl;  
  file >> variable_name >> width_left;
  cout << "Read: " << variable_name << "\t" << width_left << endl;  
  file >> variable_name >> width_right;
  cout << "Read: " << variable_name << "\t" << width_left << endl;  
  file >> variable_name >> k_vec;  
  cout << "Read: " << variable_name << "\t" << k_vec << endl;  
  file >> variable_name >> factor_on;
  cout << "Read: " << variable_name << "\t" << factor_on << endl;
  
  file.close();  
  return 0;
}

//Function to read the parameters for the method to determine the string coefficients
int read_string_coefficients_parameter(string filename, StringCoefficients &method, int & flag)
{
  ifstream file;
  string variable_name;
  file.open(filename.c_str());
  int tmp;
  
  if(!file.good())
    return -1;
  
  file >> variable_name >> tmp;
  method = (StringCoefficients) tmp;
  cout << "Read: " << variable_name << "\t" << method << endl;  
  file >> variable_name >> flag;
  cout << "Read: " << variable_name << "\t" << flag << endl;  
  
  file.close();  
  return 0;
}

void read_string_configs(string filename, vector<string_t> &configs)
{
  configs.clear();
  ifstream file;
  string_t curr_confi;
  if(!file_exists(filename))
    cout << "Warning, file " << filename << " for the string configurations does not exist" << endl;
  else
  {
    file.open(filename.c_str());
    while(!file.eof())
    {
      file >> curr_confi.pos >> curr_confi.leng;
      configs.push_back(curr_confi);
    }
    
    file.close();
  }
}
