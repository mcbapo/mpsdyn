/**
   \file SU2_Evo.cpp
   Realisation of the truncated SU(2) Hamiltonian with staggered fermions from Erez, Ignacio and Benni
   
  \param <N> (int) Number of Fermions
  \param <D> (int) Bond dimension
  \param <epsilon> (double) Hopping parameter
  \param <mu> (double) Mass
  \param <g> (double) Coupling constant
  \param <dt> (double) Time step size
  \param <name> (string) Name of a potential ground state file (optional)
  \param <leng> (int) length of the string, if 0 is given no string is placed
  
  \author Stefan Kühn
  \date 30/03/2015

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
#define BACKUP_INTERVAL  (100*RESOLUTION)
#define STRING_PARAM_FILE "string_params.dat"
#define TIMEFILE "time.txt"
#define PROJECTION_ON true

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
  QsqCorrelation,	//Determine the charge square correlations
  Entanglement		//Entanglment entropy of each bipartition
};

/** Print MPS, in contrast to the regular streaming operator this shows explicitly the content of the matrices */
void print_MPS(MPS);

//Determine the current configuration
vector<complex_t> current_configuration(vector<MPO*> &Spin1,vector<MPO*> &Spin2,vector<MPO*> &Gauge,MPS &mps,Contractor& contractor,int N);
/** Save the current configuration to a given file */
int save_configuration(string filename, vector<complex_t> config);
int save_configuration_vacuum(string filename, vector<complex_t> config);

/** Function to generate a backup of the current state during the time evolution*/
void generate_backup(const MPS &mps,int N,int D,double x,double mu, double g, double dt, int step, double curr_time);
int load_backup(MPS &mps,int &N,int &D,double &x,double &mu, double &g, double &dt, int &step, double &curr_time);

/** Print the current configuration in a readable way */
void print_configuration(vector<complex_t> config);

//Function to read the parameters for an initial superposition form a file
int read_string_parameter(string filename, bool &antiparticles, int &pos_left, int &pos_right, double &width_left, double &width_right, double &k_vec, bool &factor_on);

//All kinds of functions to determine the string coefficients
void get_OverlapStrongCouplingStringstate(SU2Hamiltonian &H, vector<MPO*> &ops, int N, int &num_ops);
void get_SdaggerS(SU2Hamiltonian &H, vector<MPO*> &ops, int N, int &num_ops);
void get_ReducedOverlapStrongCouplingStringstates(SU2Hamiltonian &H,Contractor &contractor, vector<MPO*> &ops, int N, int &num_ops);
void get_LinkProjectorCoefficients(SU2Hamiltonian &H, vector<MPO*> &ops, int N, unsigned int method, int &num_ops);
void get_SpinProjectorCoefficients(SU2Hamiltonian &H, vector<MPO*> &ops, int N, int &num_ops);
void getChargeSquareCorrelations(SU2Hamiltonian &H, vector<MPO*> &ops, int N, int method, int &num_ops);
void getEntanglementEntropy(SU2Hamiltonian &H, vector<MPO*> &ops, int N, int method, int &num_ops);
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
	 << "file:    " << "Name of the file for an exisiting state to be further evolved (optional)" << endl
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
  //Name of the ground state to import
  string gs_filename="";
  //Name of a possible state to import and go on involving
  string state_filename="";
  //Length of the initial string
  int leng=0;  
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
  int num_ops;
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
    gs_filename=argv[++cntr];
  if(argc >= NUM_OF_PARAMS+3)
    leng=atoi(argv[++cntr]);
  if(argc >= NUM_OF_PARAMS+4)
    state_filename=argv[++cntr];
  if(argc >= NUM_OF_PARAMS+5)
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
  if(leng!=0)
  {
    /*int startpos = round(N/2)-round(leng/2);
    int endpos = startpos+leng;
    masses[startpos-1]=masses[startpos-1]*100;
    masses[endpos-1]=masses[endpos-1]*100;*/
    //Another option which Erez came up with, have a site dependend mass in such a was that it mimics an external electric field
    /*double lambda=0.3;
    if(file_exists("mass.txt"))
    {
      ifstream massfile;
      massfile.open("mass.txt");
      massfile >> lambda;
      massfile.close();
      cout << "Read lambda: " << lambda << endl;
    }
    
    for(int i=0; i<masses.size(); i++)
    {
      masses[i] = masses[i]-2.0*lambda*pow(-1.0,i+1)*i;
    }*/
    
    //Box around the mesons
    for(int i=0; i<masses.size(); i++)
      masses[i] = 20.0*masses[i];
    
    masses[round(N/4.0)-1] = mu;
    masses[round(N/4.0)] = mu;
    masses[round(N/4.0)+1] = mu;
    
    masses[round(3.0*N/4.0)-3] = mu;
    masses[round(3.0*N/4.0)-2] = mu;
    masses[round(3.0*N/4.0)-1] = mu;
    
    cout << "Masses: " << masses << endl;
    
  }
  
  /***********************************************************************************
   * 				Set up MPOs and MPSs				     *
   **********************************************************************************/  
  //Get the Hamiltonian
  SU2Hamiltonian HNGI(N,epsilon,mu,g);
  //Get the Hamiltonian with site depended masses
  //SU2Hamiltonian HNGI(N,epsilon,masses,g);
  
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
  HNGI.constructEvolutionMPO(Umpo,dt);
  //HNGI.constructEvolutionMPOSiteDepMass(Umpo,dt,false);
  //HNGI.constructEvolutionMPOSiteDepMassCharge(Umpo,dt,chargepos,lambda,true);
  
  //Get the projection operators for projecting into the gauss law fulfilling subspace
  MPO Podd(1),Peven(1);
  HNGI.constructProjectors(Podd,Peven);
  
  //MPS for the state and auxiliary state and the interacting vacuum
  MPS vacuum,mps,aux;
  
  //MPOs to monitor the Gauss Law
  MPO Gx(1), Gy(1), Gz(1);
  HNGI.getGaussLawMPO(Gx,x_comp);
  HNGI.getGaussLawMPO(Gy,y_comp);
  HNGI.getGaussLawMPO(Gz,z_comp);
  
  //MPO to monitor the charge
  MPO charge_square(1);  
  HNGI.getChargeSquareMPO(charge_square);
  
  //Construct an initial state
  double initial_energy;
  if((gs_filename.length()==0) || !file_exists(gs_filename)) 
  {
    HNGI.constructInitialMPS(vacuum);
    time(&start_local);
    contractor.findGroundState(hamil,D,&initial_energy,vacuum);    
    time(&end_local);
    file.open(TIMEFILE,ios_base::app);
    file << "Time for computing the ground state\t" << difftime(end_local,start_local) << "s " << endl;
    file.close();    
    vacuum.exportMPS("Groundstate.dat");
  }
  else
  {
    cout << "Importing ground state" << endl;
    vacuum.importMPS(gs_filename.c_str());
    initial_energy=contractor.contract(vacuum,hamil,vacuum).re;
  }
  
  cout << "Vacuum energy: " << initial_energy << endl;
  
  //In case I want a string on top of it
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
     * Get a meson superpostion around a predefined value with momentum
     * ****************************************************/
    MPS *tmp_mps1,*tmp_mps2;
    MPO MesonOp(1);
    MPS aux;
    //HNGI.constructInitialMPS(vacuum);
    int pos_left=round(N/4.0)-1, pos_right=round(3.0*N/4.0);
    //HNGI.getMovingMeson(MesonOp,pos_left,2.0,0.5);
    HNGI.getMovingAntimesonMeson(MesonOp,pos_left,2.0,0.5);
    tmp_mps1=new MPS(vacuum);	tmp_mps1->approximate(MesonOp,vacuum,D);
    contractor.optimize(MesonOp,vacuum,*tmp_mps1,D);
    //HNGI.getMovingMeson(MesonOp,pos_right,2.0,-0.5);
    HNGI.getMovingAntimesonMeson(MesonOp,pos_right,2.0,-0.5);
    tmp_mps2=new MPS(vacuum);	tmp_mps2->approximate(MesonOp,vacuum,D);
    contractor.optimize(MesonOp,vacuum,*tmp_mps2,D);
    vector<const MPS*> states; states.push_back(tmp_mps1); states.push_back(tmp_mps2);
    vector<complex_t> phases; 
    phases.push_back(ONE_c); phases.push_back(ONE_c);
    aux=MPS(vacuum); aux.setRandomState();
    contractor.optimizeSum(states,phases,aux,D);
    aux.gaugeCond('R',true);
    mps=MPS(aux);
    tmp_mps1->clear(); delete tmp_mps1;
    tmp_mps2->clear(); delete tmp_mps2;
    /*******************************************************
     * Single meson superpostion around a predefined value with momentum
     * ****************************************************/
    /*MPS *tmp_mps1;
    MPO MesonOp(1);
    //HNGI.constructInitialMPS(vacuum);
    int pos=round(N/2.0);
    HNGI.getMovingMeson(MesonOp,pos,0.6,0.8);
    tmp_mps1=new MPS(vacuum);	tmp_mps1->approximate(MesonOp,vacuum,D);
    contractor.optimize(MesonOp,vacuum,*tmp_mps1,D);
    tmp_mps1->gaugeCond('R',true);
    mps=MPS(*tmp_mps1);
    tmp_mps1->clear(); delete tmp_mps1;*/
    /*******************************************************
     * Simply two mesons
     * ****************************************************/
    /*MPS *tmp_mps1,*tmp_mps2;
    MPO StringOp(1);
    MPS aux;
    int pos_left=round(N/4.0), pos_right=round(3.0*N/4.0)-1;
    cout << "pos left: " << pos_left << endl;
    cout << "pos right: " << pos_right << endl;
    HNGI.getStringMPO2(StringOp,1,pos_left,false);
    tmp_mps1=new MPS(vacuum);	tmp_mps1->approximate(StringOp,vacuum,D);
    contractor.optimize(StringOp,vacuum,*tmp_mps1,D);
    HNGI.getStringMPO2(StringOp,1,pos_right,false);
    tmp_mps2=new MPS(vacuum);	tmp_mps2->approximate(StringOp,vacuum,D);
    contractor.optimize(StringOp,vacuum,*tmp_mps2,D);
    vector<const MPS*> states; states.push_back(tmp_mps1); states.push_back(tmp_mps2);
    vector<complex_t> phases; 
    phases.push_back(ONE_c); phases.push_back(ONE_c);
    aux=MPS(vacuum); aux.setRandomState();
    contractor.optimizeSum(states,phases,aux,D);
    aux.gaugeCond('R',true);
    mps=MPS(aux);
    tmp_mps1->clear(); delete tmp_mps1;
    tmp_mps2->clear(); delete tmp_mps2;*/
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
     * Create a "heavy string" state on top of the ground state
     * ****************************************************/ 
     /*MPO StringOp(1);  
     HNGI.getHeavyString(StringOp,leng,0,false);
     aux.approximate(StringOp,vacuum,D);    
     contractor.optimize(StringOp,vacuum,aux,D);
     aux.gaugeCond('R',true);
     mps=MPS(aux);*/
  }
  else
  {
    //I simply start with my (interacting) vacuum
    mps=MPS(vacuum);
  }
  
  //Only determine the vacuum configuration and the first configuration in case I'm not further evolving a backup
  if(backup_data_flag!=0)
  {  
    //Determine the vacuum configuration
    cout << "Configuration vacuum"<< endl;
    configuration=current_configuration(Spin1_MPOs,Spin2_MPOs,Gauge_MPOs,vacuum,contractor,N);  
    print_configuration(configuration); 
    save_configuration_vacuum(FILENAME_CONFIG,configuration);
    cout << "Charge squared vacuum: " << contractor.contract(vacuum,charge_square,vacuum) << endl;
    
    //Determine the charge square of the vacuum
    coefficients=get_ChargeSquare(contractor,charge,vacuum,N);
    file.open(CHARGE_FILE, ios_base::app);
    file << "#Vacuum:" << endl;
    file << 0 << "\t" << 0;
    for(int ind=0; ind<coefficients.size(); ind++)
      file << "\t" << coefficients[ind];
    file << endl;
    file << "#State:" << endl;
    file.close();
    
    //Compute the string observables of the vacuum
    coefficients=get_StringCoefficients(HNGI,contractor,string_observable_ops,method,vacuum,num_ops);
    file.open(FILENAME_COEFFIS, ios_base::app);
    file << "#Vacuum:" << endl;
    file << 0 << "\t" << 0 ;
    for(int ind=0; ind<coefficients.size(); ind++)
      file << "\t" << coefficients[ind];
    file << endl;
    file << "#State:" << endl;
    file.close();
  }
    
  //Determine initial values
  E0=contractor.contract(mps,hamil,mps); 
  
  Gx_val=contractor.contract(mps,Gx,mps);
  Gy_val=contractor.contract(mps,Gy,mps);
  Gz_val=contractor.contract(mps,Gz,mps);
  
  cout << "Inital energy: " << contractor.contract(mps,hamil,mps) << endl;
  cout << "Gx: " << Gx_val << endl;
  cout << "Gy: " << Gy_val << endl;
  cout << "Gz: " << Gz_val << endl;    
  
  if(backup_data_flag != 0)
  {
    //Save configuration once before evolution
    configuration=current_configuration(Spin1_MPOs,Spin2_MPOs,Gauge_MPOs,mps,contractor,N);
    save_configuration(FILENAME_CONFIG,configuration);
    
    //Print start configuration
    cout << "Start configuration" << endl;
    print_configuration(configuration); 
    
    
    //Determine the charge once before evolution
    coefficients=get_ChargeSquare(contractor,charge,mps,N);
    file.open(CHARGE_FILE, ios_base::app);
    file << 0 << "\t" << 0 ;
    for(int ind=0; ind<coefficients.size(); ind++)
      file << "\t" << coefficients[ind];
    file << endl;
    file.close();
    
    //Compute some string observables once before evolution
    coefficients=get_StringCoefficients(HNGI,contractor,string_observable_ops,method,mps,num_ops);
    file.open(FILENAME_COEFFIS, ios_base::app);
    file << 0 << "\t" << 0 ;
    for(int ind=0; ind<coefficients.size(); ind++)
      file << "\t" << coefficients[ind];
    file << endl;
    file.close();    
  }
    
  //Save initial state (which is the vacuum with probably a string on top)
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
  
  //Write to file once before evolution in case it is not coming from a backup
  if(backup_data_flag!=0)
  {
    file.open(EVOFILE, ios_base::app);
    file << 0 << "\t" << 0 << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << setprecision(10) << E0 << "\t" << Gx_val << "\t" << Gy_val << "\t" << Gz_val  << "\t" << abs(contractor.contract(init_state,mps)) << "\t" << contractor.contract(mps,charge_square,mps) << "\t" << sum_of_errors << endl;
    file.close();
  }
  
  bool mass_flag=1;
  int mass_counter=0;
  /***********************************************************************************
   *				Actual evolution				     *
   **********************************************************************************/
  time(&start_local);
  cout << "------------------------------------------------------------------------"<<endl;
  cout << "------------------------Starting to evolve state------------------------"<<endl;
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
      //Simple first order Taylor expansion, not unitary anymore
      curr_err=-1.0;
      contractor.optimize(Umpo,mps,aux,D,&curr_err);
      errors.push_back(curr_err);
      aux.setNormFact(1.0);
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
	beta.push_back(normfactors[k]*power(-I_c*dt,k)/factorial(k)*ONE_c);
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
    
    //Padé approximant, compared to a Taylor expansion this is unitary
    //WARNING: Not working
    //contractor.approximateExponential(hamil,mps,aux,dt,D,0.0,0.0);

    
    if(PROJECTION_ON && (i%PROJECTION_INTERVAL)==0)
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
      
      //Compute some string observables
      coefficients=get_StringCoefficients(HNGI,contractor,string_observable_ops,method,aux,num_ops);
      file.open(FILENAME_COEFFIS, ios_base::app);
      file << i << "\t" << curr_time ;
      for(int ind=0; ind<coefficients.size(); ind++)
	file << "\t" << coefficients[ind];
      file << endl;
      file.close();     
    }
    mps=aux;
    /*if(mass_flag && (curr_time>0.0) && (i%100)==0)
    {
      cout.flush(); 
      //mass_flag=0;
      cout << endl << "Switching mass at time " << curr_time << endl;
      mass_counter++;
      if(curr_time<=10.0E-1)
      {
	masses.assign(N,20.0*mu);
	masses[round(N/4.0)-1+mass_counter] = mu;
	masses[round(N/4.0)+mass_counter] = mu;
	masses[round(N/4.0)+1+mass_counter] = mu;
	
	masses[round(3.0*N/4.0)-3-mass_counter] = mu;
	masses[round(3.0*N/4.0)-2-mass_counter] = mu;
	masses[round(3.0*N/4.0)-1-mass_counter] = mu;
      }
      else
	masses.assign(N,mu);
      
      HNGI.set_mu(masses);
      HNGI.updateHMPOSiteDepMass();
      HNGI.constructEvolutionMPOSiteDepMass(Umpo,dt,false);      
    }*/
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

void print_MPS(MPS mps)
{
  cout << "MPS of length " << mps. getLength() << endl;
  for(int i=0; i<mps.getLength(); i++)
  {
    cout << "[" << i << "]=" <<  (mps.getA(i)).getA() << endl;
  }
}


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

int save_configuration_vacuum(string filename, vector<complex_t> config)
{
  ofstream file;
  file.open(filename.c_str(),ios_base::app);
  file << "#Vacuum:" << endl;
  for(int k=0; k<config.size()-1; k++)
    file << config[k] << "\t";
  file << config.back() << endl;
  file << "#State:" << endl;
  file.close();
  return 0;
}
  
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

/*vector<complex_t> get_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS& mps, MPS& vacuum, int N)
{
  MPO mpo(1);
  MPS strong_coupling;
  vector<complex_t> res;
  complex_t temp;
  
  //I want the strong coupling vacuum not the interacting vacuum I get as input
  H.constructInitialMPS(strong_coupling);
  
  for(int leng=1; leng<=N; leng+=2)
  {
    for(int pos=1; pos<=N-leng; pos++)
    {
      H.getStringMPO2(mpo,leng,pos);
      //WARNING: contract is taking arguments in this order: Ket, MPO, Bra !!!
      //temp = contractor.contract(vacuum,mpo,mps)/sqrt(abs(contractor.contract2(mpo,vacuum)));
      temp = contractor.contract(strong_coupling,mpo,mps)/sqrt(abs(contractor.contract2(mpo,strong_coupling)));
      res.push_back(temp);
      //cout << "Created string at position " << pos << " with length " << leng << " with overlap " <<  temp << endl;
    }
  }
  return res;
}*/


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
  else if(method==Entanglement)
    num_ops=3*H.getSites()-2;
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
  else if(method==Entanglement)
  {
    #pragma omp parallel for 
    for (int i=1; i<(3*H.getSites()-1); i++)
    {
      result[i-1].re = contractor.getEntropy(mps,i);
      result[i-1].im = 0;
    }
  }
  else 
  {
    #pragma omp parallel for 
    for (int i=0; i<num_ops; i++)
      result[i] = contractor.contract(mps,*ops[i],mps);
  }
  return result;
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