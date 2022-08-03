/**
   \file VacuumConfig.cpp
   I realized that I want to know some observables of the vacuum to be able to substract them afterwards. However, I did not compute them in several cases, I saved the ground state though. Therefore I have this dummy program which basically behaves like SU2_Evo2 but only computes the desired observables of the vacuum
   
  \param <N> (int) Number of Fermions
  \param <D> (int) Bond dimension
  \param <epsilon> (double) Hopping parameter
  \param <mu> (double) Mass
  \param <g> (double) Coupling constant
  \param <dt> (double) Time step size
  \param <name> (string) Name of a potential ground state file (optional)
  \param <leng> (int) length of the string, if 0 is given no string is placed
  
  \author Stefan Kühn
  \date 10/04/2015

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
#define BACKUP_INTERVAL  (100*RESOLUTION)
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
  double curr_err,sum_of_errors;
  vector<double> errors;
  //Offset in case I want to continue evolving a state
  int offset=0;
  //Method to determine the coefficients
  StringCoefficients method=OverlapSC;
  int flag=0;
  vector<MPO*> string_observable_ops;
  int num_ops;
  //Variables for higher order Taylor expansion
  
  
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
  if(argc >= NUM_OF_PARAMS+4)
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
    int startpos = round(N/2)-round(leng/2);
    int endpos = startpos+leng;
    masses[startpos-1]=masses[startpos-1]*100;
    masses[endpos-1]=masses[endpos-1]*100;
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
  //HNGI.constructEvolutionMPOSiteDepMass(Umpo,dt,true);
  //HNGI.constructEvolutionMPOSiteDepMassCharge(Umpo,dt,chargepos,lambda,true);
  
  //Get the projection operators for projecting into the gauss law fulfilling subspace
  MPO Podd(1),Peven(1);
  HNGI.constructProjectors(Podd,Peven);
  
  //MPS for the state and auxiliary state and the interacting vacuum
  MPS vacuum;
  
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
  
  //Save configuration of the vacuum
  cout << "Computing spin and flux configuration" << endl;
  configuration=current_configuration(Spin1_MPOs,Spin2_MPOs,Gauge_MPOs,vacuum,contractor,N);
  save_configuration(FILENAME_CONFIG,configuration);
  
  //Determine the charge
  cout << "Computing charge square configuration" << endl;
  coefficients=get_ChargeSquare(contractor,charge,vacuum,N);
  file.open(CHARGE_FILE, ios_base::app);
  file << 0 << "\t" << 0;
  for(int ind=0; ind<coefficients.size(); ind++)
    file << "\t" << coefficients[ind];
  file << endl;
  file.close();
  
  //Compute some string observables
  cout << "Computing string observables" << endl;
  coefficients=get_StringCoefficients(HNGI,contractor,string_observable_ops,method,vacuum,num_ops);
  file.open(FILENAME_COEFFIS, ios_base::app);
  file << 0 << "\t" << 0 ;
  for(int ind=0; ind<coefficients.size(); ind++)
    file << "\t" << coefficients[ind];
  file << endl;
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

vector <complex_t> get_ChargeSquare(Contractor &contractor, vector<MPO*> ops, MPS mps, int N)
{
  vector<complex_t> result;
  result.resize(N);
  
  #pragma omp prallel for
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