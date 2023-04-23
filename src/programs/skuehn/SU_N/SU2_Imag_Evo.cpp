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
  \date 16/03/2015

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
#define FILENAME_EX "SU2_ex.txt"
#define FILENAME_CONFIG "SU2_configuration.txt"
#define FILENAME_COEFFIS "SU2coefficients.txt"
#define CHARGE_FILE "SU2charge.txt"
#define EVOFILE "SU2efile.txt"
#define RESOLUTION 10
#define PROJECTION_INTERVAL 100
#define BACKUP_INTERVAL (500*RESOLUTION)
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
vector<complex_t> current_configuration(vector<MPO*> &Spin1,vector<MPO*> &Spin2,vector<MPO*> &Gauge,MPS &mps,Contractor& contractor);

/** Save the current configuration to a given file */
int save_configuration(string filename, vector<complex_t> config);

/** Print the current configuration in a readable way */
void print_configuration(vector<complex_t> config);

/** Determine the projection on all string states arising from the strong coupling vacuum specified in the vector configurations */
vector<complex_t> get_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, vector< vector<str> > configurations, int N);
/** Little help function, given a string of maximum length determine all configurations with strings on the state which are possible*/
void get_all_strings_helper(int start, int end, int leng, vector<str> tmp, vector< vector<str> > &res);
/** Given a start site and and end site and a maximum length, determine all possible string configurations on the strong coupling vacuum containing strings in the specified region with length at most the given maximum length.*/
vector< vector<str> > get_all_strings(int start, int end, int leng);

/** Function to generate a backup of the current state during the time evolution*/
void get_snapshot(const MPS &mps,int N,int D,double x,double mu, double g, double dt, int step);
//Function to read the parameters for an initial superposition form a file
int read_string_parameter(string filename, bool &antiparticles, int &pos_left, int &pos_right, double &width_left, double &width_right, double &k_vec, bool &factor_on);
//Get full config of a state to check the Gauss Law manually afterwars
void get_full_config(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N);

//All kinds of functions to determine the string coefficients
vector<complex_t> get_OverlapStrongCouplingStringstate(SU2Hamiltonian &H, Contractor &contractor, MPS& mps, int N);
vector<complex_t> get_SdaggerS(SU2Hamiltonian &H, Contractor &contractor, MPS &mps,int N);
vector<complex_t> get_ReducedOverlapStrongCouplingStringstates(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N);
vector<complex_t> get_LinkProjectorCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N, unsigned int method);
vector<complex_t> get_SpinProjectorCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N);
vector<complex_t> getChargeSquareCorrelations(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N, int method);
//Driver function to compute the string coefficients
vector<complex_t> get_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N, StringCoefficients method, int flag);

//Function to read the parameters for the method to determine the string coefficients
int read_string_coefficients_parameter(string filename, StringCoefficients &method, int & flag);

//Simple factorial function to compute Taylor coefficients (as int might not be enough, I use double to prevent overflows)
double factorial(int n)
{
  double res=1.0;
  for(int i=n; i>0; i--)
    res *= (double) i;
  return res;
}

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
	 << "GS:      " << "Name of the initial state which should be imported (optional)" << endl
	 << "leng:    " << "Length of the initial string (optional)" << endl
	 << "lambda:  " << "Strength of the penalty for the static charges (optional)" << endl 	 << "file:    " << "Name of the file for an exisiting state to be further evolved (optional)" << endl
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
  //Name of the ground state to import
  string gs_filename="";
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
  double curr_err,sum_of_errors;
  vector<double> errors;
  //Offset in case I want to continue evolving a state
  int offset=0;
  //Method to determine the coefficients
  StringCoefficients method=OverlapSC;
  int flag=0;
  
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
    gs_filename=argv[++cntr];
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
  
  const MPO& hamil=HNGI.getHMPO();
  
  //Turn gauge interaction off in case this is wanted
#ifdef GAUGEOFF
  HNGI.removeGaugeInteraction();
  HNGI.updateHMPO();
#endif
  
  //MPS for the state and auxiliary state and the interacting vacuum
  MPS mps,aux;
  
  //Check if it is hermitian
  /*HNGI.constructInitialMPS(mps);
  mps.increaseBondDimension(D);
  for (int i=0; i< 10000; i++)
  {
    mps.setRandomState();
    cout << i << ": <mps|H|mps>=" << contractor.contract(mps,hamil,mps) << ", <mps|mps>="<< contractor.contract(mps,mps) << endl;
  }*/
  
  //MPOs to monitor the Gauss Law
  MPO Gx(1), Gy(1), Gz(1);
  HNGI.getGaussLawMPO(Gx,x_comp);
  HNGI.getGaussLawMPO(Gy,y_comp);
  HNGI.getGaussLawMPO(Gz,z_comp);
  
  //In case I want a string on top of it
  if(leng!=0)
  {
    /*******************************************************
     * Generate a string (on top of the interacting vacuum)
     * ****************************************************/
//    HNGI.createString(mps,contractor,min(2,leng+1),leng);
     //MPO StringOp(1);
     //MPS aux;    
     //HNGI.getStringMPO2(StringOp,leng,0,false);
//     HNGI.getMovingString(StringOp);
//     cout << "WARNING will create a moving string!" << endl;
//     
//     cout << "Successfully retrieved String Operator" << endl;
//     /*cout << "Consistency check 1 String Operator: " << contractor.contract2(StringOp,vacuum) << endl;
//     cout << "Consistency check 2 String Operator: " << contractor.contract2(vacuum,StringOp) << endl;*/
     /*cout << "Creating a string with length " << leng << endl;
     aux.approximate(StringOp,vacuum,D);    
     contractor.optimize(StringOp,vacuum,aux,D);
     aux.gaugeCond('R',true);*/
//     /*cout << "WARNING: String switched off" << endl;
//     aux=MPS(vacuum);*/
     //mps=MPS(aux);
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
    
    //Second option (each color species differently)
    /*HNGI.getMovingQuarks2(Moving,3.0,1);
    aux.approximate(Moving,vacuum,D);    
    contractor.optimize(Moving,vacuum,aux,D);
    HNGI.getMovingQuarks2(Moving,0.0,2);
    mps.approximate(Moving,aux,D);
    contractor.optimize(Moving,aux,mps,D);
    mps.gaugeCond('R',true);*/
    
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
  }
  else
  {
    HNGI.constructInitialMPS(mps);
  }
  
  //Case that I want to import a state which replaces the MPS
  //TODO: Make it more efficient, avoid computing a possible string state in case it will be replaced anyway
  if(file_exists(state_filename))
    mps.importMPS(state_filename.c_str());
  else if(state_filename.length()>0)
    cout << "Warning, could not import state" << endl;
  
  cout << "------------------------------------------------------------------------"<<endl;
  cout << "------------------------Starting to evolve state------------------------"<<endl;
  cout << "------------------------------------------------------------------------"<<endl; 
  
  //Determine initial values of energy and Gauss Law
  E0=contractor.contract(mps,hamil,mps); 
  
  Gx_val=contractor.contract(mps,Gx,mps);
  Gy_val=contractor.contract(mps,Gy,mps);
  Gz_val=contractor.contract(mps,Gz,mps);
  
  cout << "Gx: " << Gx_val << endl;
  cout << "Gy: " << Gy_val << endl;
  cout << "Gz: " << Gz_val << endl;   
  
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
  
  //Get the time evolution MPO
  MPO Umpo(1);
  HNGI.constructEvolutionMPO(Umpo,dt,true);
  //HNGI.constructEvolutionMPOSiteDepMass(Umpo,dt,true);
  //HNGI.constructEvolutionMPOSiteDepMassCharge(Umpo,dt,chargepos,lambda,true);
  
  //Get the projection operators for projecting into the gauss law fulfilling subspace
  MPO Podd(1),Peven(1);
  HNGI.constructProjectors(Podd,Peven);
  
  //Save configuration once before evolution
  configuration=current_configuration(Spin1_MPOs,Spin2_MPOs,Gauge_MPOs,mps,contractor);
  save_configuration(FILENAME_CONFIG,configuration);
  
  //Print start configuration
  cout << "Start configuration" << endl;
  print_configuration(configuration); 
  
  cout << "Inital energy: " << contractor.contract(mps,hamil,mps) << endl;
    
  //Save initial state of the time evolution
  MPS init_state(mps);
  
  //MPO to monitor the charge
  MPO charge(1),charge_square(1);
  
  HNGI.getChargeSquareMPO(charge_square);
    
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
  
  //Possible string configurations
  vector< vector<str> > str_confis=get_all_strings(1,N,N-1);
    
  
  /***********************************************************************************
   *				Actual evolution				     *
   **********************************************************************************/
  //Monitor the energy during evolution and decide whether it's time to stop or to reduce the time step or something else
  complex_t energy_old=E0;
  double curr_time=0.0;
  //Variables for hihger order Taylor expansion
  vector< const MPS* >  kets;
  vector< complex_t >  beta;
  vector< double >  normfactors;
  int order;
  cout << "Order: " ;
  cin >> order;
  //Compute the taylor coefficients
  complex_t taylor_coefficient;
  MPS *tmp_mps;
  
  cout << "Taylor coefficients: " << beta << endl;
  time(&start_local);
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
    ///////////////////////////////////////////////////////////
    //Simple first order Taylor expansion, not unitary anymore
    ///////////////////////////////////////////////////////////
    if(order==1)
    {
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
      cout << "Building states" << endl;
      normfactors.push_back(1.0);
      beta.push_back(ONE_c);
      for(int k=1; k<=order; k++)
      {
	curr_err=-1.0;
	aux.approximate(hamil,mps,D);
	contractor.optimize(hamil,mps,aux,D,&curr_err);
	normfactors.push_back(normfactors[k-1]*sqrt(contractor.contract(aux,aux).re));
	beta.push_back(normfactors[k]*pow(-1.0*dt,k)/factorial(k)*ONE_c);
	aux.setNormFact(1.0);
	tmp_mps=new MPS(aux);
	kets.push_back(tmp_mps);
	cout << "<mps|aux>=" << contractor.contract(mps,aux) << endl;
	cout << "<mps|H|mps>=" << contractor.contract(mps,hamil,mps)<<", <mps|mps>=" << contractor.contract(mps,mps) << endl;
	mps=aux;
	aux.setRandomState();
      }    
      
      //contractor.optimizeSum(kets,beta,aux,D);
      contractor.optimizeSum(kets,beta,aux,D);
      aux.gaugeCond('R',true);
      //Free memory
      cout << "Freeing memory" << endl;
      for (int k=0; k<kets.size(); k++)
	delete kets[k];
      kets.clear();
      normfactors.clear();
      beta.clear();
    }
    
    //Padé approximant, compared to a Taylor expansion this is unitary
    //WARNING: Not working
    //contractor.approximateExponential(hamil,mps,aux,dt,D,0.0,0.0);
    
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
      configuration=current_configuration(Spin1_MPOs,Spin2_MPOs,Gauge_MPOs,mps,contractor);
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
	get_snapshot(aux,N,D,epsilon,mu,g,dt,i);
      
      //Determine the charge
      file.open(CHARGE_FILE, ios_base::app);
      file << i << "\t" << curr_time ;
      for(int pos=1; pos<=N; pos++)
      {	
	HNGI.getChargeSquareMPO(charge,pos);
	file << "\t" << contractor.contract(aux,charge,aux);
      }
      file << endl;
      file.close();     
      
      //Compare to old energy value and act accordingly
      //if((i%(10*RESOLUTION))==0)
      /*if(i==3000)
      {
	//if((abs((E0-energy_old)/energy_old)<1.0E-5) && (dt>=1.0E-4))
	//{
	  dt /=4.0;
	  HNGI.constructEvolutionMPO(Umpo,dt);
	  cout << "Relative change: " << abs((E0-energy_old)/energy_old) <<", reducing time step, dt=" << dt << endl;
	//}
	//energy_old=E0;
      } */
      
      if(i==40)
      {
	dt /=100.0;
	order=1;
	HNGI.constructEvolutionMPO(Umpo,dt);
      } 
      
      /*if(i==500)
      {
	dt /=10.0;
	HNGI.constructEvolutionMPO(Umpo,dt);	  
      }*/
      
      //Compute some string observables
      /*coefficients=get_StringCoefficients(HNGI,contractor,aux,N,(StringCoefficients) method,flag);
      //Determine the charge
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
vector<complex_t> current_configuration(vector<MPO*> &Spin1, vector<MPO*> &Spin2, vector<MPO*> &Gauge,MPS &mps,Contractor& contractor)
{
  complex_t tmp;
  vector<complex_t> res;
  if((Spin1.size()!=Spin2.size()) || (Spin1.size()!=(Gauge.size()+1)))
  {
    cout << "Error, cannot determine the configuration" << endl;
    exit(666);
  }
  else
  {
    for(int k=0; k<Spin1.size(); k++)
    {
      tmp = contractor.contract(mps,*Spin1[k],mps);
      res.push_back(tmp);
      tmp = contractor.contract(mps,*Spin2[k],mps);
      res.push_back(tmp);
      if(k<Gauge.size())
      {
	tmp = contractor.contract(mps,*Gauge[k],mps);
	res.push_back(tmp);
      }
    }
  }
  return res;
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

//Helper function to construct all possible string configurations
void get_all_strings_helper(int start, int end, int leng, vector<str> tmp, vector< vector<str> > &res)
{
  int curr_leng=0;
  //Local copy of current path
  vector<str> local_tmp;
  str curr_str;
  
  //1. case: new start site >= end site, no strings possible, therefore end
  if(start>=end)
  {
    //Did I put some strings, in case yes save result
    if(tmp.size()>0)
      res.push_back(tmp);
  }
  //2. case: Strings can still be placed
  else
  {
    //Loop over length
    curr_leng=0;
    do{
      //Treat every possibility separate
      local_tmp.clear();
      local_tmp=tmp;
      //Case that current string still fits on chain
      if(start+curr_leng<=end)
      {
	if(curr_leng>0)
	{
	  curr_str.pos=start;
	  curr_str.leng=curr_leng;
	  local_tmp.push_back(curr_str);
	}
	get_all_strings_helper(start+curr_leng+1,end,leng,local_tmp,res);	
      }
      //Current string does not fit anymore, as all others are even longer I can stop
      else
	break;
      (curr_leng==0) ? curr_leng++ : curr_leng+=2;
    }while(curr_leng<=leng);
  }
  
}

//Determine all possible string configurations in between start and end with a maximum length of leng
vector< vector<str> > get_all_strings(int start, int end, int leng)
{
  //Vector to keep the result
  vector< vector<str> > res;
  //Vector to keep track of the current path
  vector<str> tmp;
  //Now call the helper
  get_all_strings_helper(start,end,leng,tmp,res);
  return res;
}

//This routine computes all possible string configurations arising form the strong coupling ground state, which are EXPONENTIALLY MANY! Use with care for small systems only
vector<complex_t> get_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS& mps,vector< vector<str> > configuration, int N)
{
  MPO* mpo;
  MPS strong_coupling;
  vector<complex_t> res;
  vector<str> curr_confi;
  complex_t temp;
  
  const MPO* StrMPOs[N];
  
  //I want the strong coupling vacuum not the interacting vacuum I get as input
  H.constructInitialMPS(strong_coupling);
  
  for(int num_confis=0; num_confis<configuration.size(); num_confis++)
  {
    curr_confi=configuration[num_confis];    
    //StrMPOs = new MPO[curr_confi.size()];
    cout << "String:" ;
    for(int k=0; k<curr_confi.size(); k++)
    {
      cout << " (" << curr_confi[k].pos <<","<<curr_confi[k].leng<<")"; 
      mpo=new MPO(1);
      H.getStringMPO2(*mpo,curr_confi[k].leng,curr_confi[k].pos);
      StrMPOs[k]=mpo;
    }
    cout << endl;
    //For safety put NULL pointers in unused spaces;
    for(int k=curr_confi.size(); k<N; k++)
      StrMPOs[k]=NULL;
   
    //I join the operators which makes them behave like one in case I have more than one
    mpo=new MPO(1);
    MPO::join(curr_confi.size(),StrMPOs,*mpo);
    temp = contractor.contract(strong_coupling,*mpo,mps)/sqrt(abs(contractor.contract2(*mpo,strong_coupling)));    
    res.push_back(temp);
    //Free memory
    mpo->clear();
    delete mpo;  
    for(int k=0; k<curr_confi.size(); k++)
      delete StrMPOs[k];
    
    //cout << "coefficient: " << pow(abs(temp),2.0) << endl;
  }
  return res;
}

vector<complex_t> get_OverlapStrongCouplingStringstate(SU2Hamiltonian &H, Contractor &contractor, MPS& mps, int N)
{
  MPO mpo(1);
  MPS strong_coupling;
  vector<complex_t> res;
  complex_t temp;
  
  cout << "Using overlap with the strong coupling string state" << endl;
  
  //I want the strong coupling vacuum not the interacting vacuum I get as input
  H.constructInitialMPS(strong_coupling);
  
  for(int leng=1; leng<=N; leng+=2)
  {
    for(int pos=1; pos<=N-leng; pos++)
    {
      H.getStringMPO2(mpo,leng,pos);
      //WARNING: contract is taking arguments in this order: Ket, MPO, Bra !!!
      temp = contractor.contract(strong_coupling,mpo,mps)/sqrt(abs(contractor.contract2(mpo,strong_coupling)));
      res.push_back(temp);
    }
  }
  return res;
}

vector<complex_t> get_SdaggerS(SU2Hamiltonian &H, Contractor &contractor, MPS &mps,int N)
{
  MPO mpo(1);
  vector<complex_t> result;
  
  cout << "Using exectation value of S^dagger S" << endl;

  //Loop over length
  for(int leng=1; leng<=N-1; leng+=2)
  {
    //Loop over position
    for(int pos=1; pos<=N-leng; pos++)
    {
      /*******************************************************************************
       * Compute <\Psi(t)|SS^\dagger|\Psi>
       * ****************************************************************************/      
      //Get the adjoint string operator S^\dagger
      H.getStringMPO2(mpo,leng,pos,true);      
      result.push_back(contractor.contract2(mpo,mps));
    }
  }
  return result;
}

vector<complex_t> get_ReducedOverlapStrongCouplingStringstates(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N)
{
  vector<complex_t> result;
  
  MPO StringOp(1);
  MPS strong_coupling_vacuum,op;
  
  H.constructInitialMPS(strong_coupling_vacuum);
  op=MPS(strong_coupling_vacuum);
  op.increaseBondDimension(2);
  op.setRandomState();
  
  cout << "Using reduced overlap with the strong coupling string state" << endl;
  
  //Loop over length
  for(int leng=1; leng<=N-1; leng+=2)
  {
    //Loop over position
    for(int pos=1; pos<=N-leng; pos++)
    {
      /*******************************************************************************
       * Compute <\Psi(t)|Id_A \otimes Tr_A(S|0><0|S^\dagger) |\Psi(t)>
       * ****************************************************************************/ 
      H.getStringMPO2(StringOp,leng,pos);
      op.approximate(StringOp,strong_coupling_vacuum,2);    
      contractor.optimize(StringOp,strong_coupling_vacuum,op,2);
      H.buildReducedStringMPO(StringOp,op,pos,pos+leng);
      result.push_back(contractor.contract(mps,StringOp,mps));
    }
  }
  return result;
}

vector<complex_t> get_LinkProjectorCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N, unsigned int method)
{
  vector<complex_t> result;
  
  MPO Projector(1);
  
  if(method==0)
    cout << "Using string projector without taking the spins into account" << endl;
  else if(method==1)
    cout << "Using string projector and looking for parallel polarizations of the spins" << endl;
  if(method==2)
    cout << "Using string projector and looking for polarization in up (down) direction on odd (even) sites" << endl;  
 
  for(int leng=1; leng<=N-1; leng++)
  {
    //Loop over position
    for(int pos=1; pos<=N-leng; pos++)
    {
      /*******************************************************************************
       * Compute <\Psi(t)|P|\Psi(t)>
       * ****************************************************************************/
      H.getStringProjector(Projector,leng,pos,method);
      result.push_back(contractor.contract(mps,Projector,mps));
    }
  }
  return result;
}

vector<complex_t> get_SpinProjectorCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N)
{
  vector<complex_t> result;  
  MPO Projector(1);
  
  cout << "Using projector on spin sites" << endl;  
 
  for(int leng=1; leng<=N-1; leng++)
  {
    //Loop over position
    for(int pos=1; pos<=N-leng; pos++)
    {
      /*******************************************************************************
       * Compute <\Psi(t)|P|\Psi(t)>
       * ****************************************************************************/
      H.getCountStringsMPO(Projector,leng,pos);
      result.push_back(contractor.contract(mps,Projector,mps));
    }
  }
  return result;
}


vector<complex_t> getChargeSquareCorrelations(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N, int method)
{
  vector<complex_t> result;
  MPO correlation(1);
  
  if(method==0)
    cout << "Using old version with projector" << endl;
  else if(method==1)
    cout << "Using old version without projector" << endl;
  else
    cout << "Using new version with/without unitaries" << endl;
  
   //Loop over length
  for(int leng=1; leng<=N-1; leng++)
  {
    //Loop over position
    for(int pos=1; pos<=N-leng; pos++)
    {
      if(method==0)
	H.getChargeCorrelatorMPO(correlation,leng,pos,true);
      else if(method==1)
	H.getChargeCorrelatorMPO(correlation,leng,pos,false);
      else
	H.getChargeCorrelatorMPO2(correlation,leng,pos,0.5);
      result.push_back(contractor.contract(mps,correlation,mps));
    }
  }
  return result;
}

//Function to compute the coefficients (basically does nothing but choosing the right subroutine)
vector<complex_t> get_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N, StringCoefficients method, int flag)
{
  vector<complex_t> result;
  
  if(method==OverlapSC)
  {
    result=get_OverlapStrongCouplingStringstate(H,contractor,mps,N);
  }
  else if(method==ReducedSC)
  {
    result=get_ReducedOverlapStrongCouplingStringstates(H,contractor,mps,N);
  }
  else if(method==SdaggerS)
  {
    result=get_SdaggerS(H,contractor,mps,N);
  }
  else if(method==StringProjector)
  {
    result=get_LinkProjectorCoefficients(H,contractor,mps,N,flag);
  }
  else if(method==ParticleProjector)
  {
    result=get_SpinProjectorCoefficients(H,contractor,mps,N);
  }
  else if(method==QsqCorrelation)
  {
    result=getChargeSquareCorrelations(H,contractor,mps,N,flag);
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

//Save an MPS and the data belonging to it
void get_snapshot(const MPS &mps,int N,int D,double x,double mu, double g, double dt, int step)
{
  ofstream file;
  time_t timer;
  time(&timer);
  
  file.open("BackupData.txt");
  file << asctime(localtime(&timer))
       << "N:    " << N << endl
       << "D:    " << D << endl
       << "x:    " << x << endl
       << "mu:   " << mu << endl
       << "g:    " << g << endl
       << "dt:   " << dt << endl
       << "step: " << step << endl;  
  file.close();
  
  mps.exportMPS("MPS.dat");
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
