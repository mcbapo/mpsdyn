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
  \date 30/01/2014

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
#define NUM_OF_PARAMS 6
#define FILENAME_GAUGE "SU2_gauge.txt"
#define FILENAME_NON_GAUGE "SU2_non_gauge.txt"
#define FILENAME_CONFIG_GAUGE "SU2_gauge_config.txt"
#define FILENAME_CONFIG_NON_GAUGE "SU2_non_gauge_config.txt"
#define FILENAME_CONFFIS "SU2coefficients.txt"
#define CHARGE_GAUGE_FILE "SU2_gauge_charge.txt"
#define CHARGE_NON_GAUGE_FILE "SU2_non_gauge_charge.txt"
#define TIMEFILE "time.txt"

using namespace std;
using namespace shrt;

/**Struct representing a string */
struct str
{
  int pos;
  int leng;
};

//bool file_exists(string filename);

/** Print MPS, in contrast to the regular streaming operator this shows explicitly the content of the matrices */
void print_MPS(MPS);

//Determine the current configuration
vector<complex_t> current_configuration(vector<MPO*> &Spin1,vector<MPO*> &Spin2,vector<MPO*> &Gauge,MPS &mps,Contractor& contractor);

/** Save the current configuration to a given file */
int save_configuration(string filename, vector<complex_t> config);

/** Print the current configuration in a readable way */
void print_configuration(vector<complex_t> config);

//Determine the projection on all single string states
//vector<complex_t> get_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, MPS &vacuum, int N);
/** Determine the projection on all string states arising from the strong coupling vacuum specified in the vector configurations */
vector<complex_t> get_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, vector< vector<str> > configurations, int N);
/** Little help function, given a string of maximum length determine all configurations with strings on the state which are possible*/
void get_all_strings_helper(int start, int end, int leng, vector<str> tmp, vector< vector<str> > &res);
/** Given a start site and and end site and a maximum length, determine all possible string configurations on the strong coupling vacuum containing strings in the specified region with length at most the given maximum length.*/
vector< vector<str> > get_all_strings(int start, int end, int leng);
/** Determine if a a string on a certain position is present in the state (independently of the configuration of the other sites)*/
vector<complex_t> get_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS &mps,MPS &vacuum, int N);
/** Yet another variant. Determine if a a string on a certain position is present in the state (independently of the configuration of the other sites)**/
vector<complex_t> get_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N, unsigned int method);
/** Function to generate a backup of the current state during the time evolution*/
void get_snapshot(const MPS &mps,int N,int D,double x,double mu, double g, double dt, int step);
//Function to read the parameters for an initial superposition form a file
int read_string_parameter(string filename, bool &antiparticles, int &pos_left, int &pos_right, double &width_left, double &width_right, double &k_vec, bool &factor_on);
//Get all possible correlations
vector<complex_t> getChargeSquareCorrelations(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N);

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
	 << "leng:    " << "Length of the string to be placed" << endl
	 << "GS1:     " << "Name of the gauge invariant ground state which should be imported (optional)" << endl
	 << "GS2:     " << "Name of the  non gauge invariant ground state which should be imported (optional)" << endl;
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
  //Length of the initial string
  int leng=atoi(argv[++cntr]);
  int leng_short,leng_long;
  //Name of the ground state to import
  string gs_gauge_filename="";
  string gs_non_gauge_filename="";  
  //Variables for time measurement
  time_t start,end,start_local,end_local;  
  //File for results
  ofstream file_gauge,file_non_gauge,file;
  //Variables for Gauss Law Components
  complex_t Gx_val,Gy_val,Gz_val; 
  //Variables for the energies
  complex_t E_vacuum,E_nongauge_vaccuum,E_string,E_nongauge_string,E_string_short=ZERO_c,E_string_long=ZERO_c,E_nongauge_string_long=ZERO_c,E_nongauge_string_short=ZERO_c;
  double initial_energy;
  //Vector to save a configuration
  vector<complex_t> configuration;
  //Vector to save the coefficients for the string states
  vector<complex_t> coefficients;
  //Vector to save the correlation functions
  vector<complex_t> correlations,correlations_nongauge;
  //Variables for the error
  double curr_err,sum_of_errors;
  vector<double> errors;
  //Offset in case I want to continue evolving a state
  int offset=0;
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  contractor.setConvTol(1.0E-8);
  
  if(argc == NUM_OF_PARAMS+2)
    gs_gauge_filename=argv[++cntr];
  else if(argc == NUM_OF_PARAMS+3)
  {
    gs_gauge_filename=argv[++cntr];
    gs_non_gauge_filename=argv[++cntr];
  }

  
  //Start time measurement
  time(&start);
  
  //Provide site dependend mass
  //Another option which Erez came up with, have a site dependend mass in such a was that it mimics an external electric field
  double lambda=0.3;
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
  }
  cout << "Masses: " << masses << endl;
    
  
  
  /***********************************************************************************
   * 				Set up MPOs and MPSs				     *
   **********************************************************************************/  

  //Get the Hamiltonian
  SU2Hamiltonian HNGI(N,epsilon,mu,g);
  //Get the Hamiltonian with site depended masses
  //SU2Hamiltonian HNGI(N,epsilon,masses,g);
  
  const MPO& hamil=HNGI.getHMPO();    
  
  //MPS for the state and auxiliary state and the interacting vacuum
  MPS vacuum,mps,aux,nongauge_vacuum,nongauge_mps,state_short,state_long,state_nongauge_short,state_nongauge_long;
  
  //MPOs to monitor the Gauss Law
  MPO Gx(1), Gy(1), Gz(1);
  HNGI.getGaussLawMPO(Gx,x_comp);
  HNGI.getGaussLawMPO(Gy,y_comp);
  HNGI.getGaussLawMPO(Gz,z_comp);
  
  //MPO for the string operator
  MPO StringOp(1),StringOp2(1),JoinedString(1),StringOpdagger(1),StringOp2dagger(1),EnergyMPO(1);
  const MPO* StrMPOs[2];
  const MPO* AllMPOs[5];
  
  //MPO for total charge square and local charge
  MPO charge_square(1),charge(1);
  HNGI.getChargeSquareMPO(charge_square);
  
  cout << endl;
  
  /***********************************************************************************
   * 				Gauge invariant case				     *
   **********************************************************************************/
  
  //Construct an initial state
  HNGI.constructInitialMPS(vacuum);
  if((gs_gauge_filename.length()==0) || !file_exists(gs_gauge_filename)) 
  { 
    time(&start_local);
    contractor.findGroundState(hamil,D,&initial_energy,vacuum);    
    time(&end_local);
    file.open(TIMEFILE,ios_base::app);
    file << "Time for computing the ground state\t" << difftime(end_local,start_local) << "s " << endl;
    file.close();
    
    vacuum.exportMPS("GS_gauge.dat");
  }
  else
  {
    cout << "Importing gauge invariant ground state" << endl;
    vacuum.importMPS(gs_gauge_filename.c_str());
  }  
  E_vacuum=contractor.contract(vacuum,hamil,vacuum);
  cout << "Initial energy: " << E_vacuum << endl;

  
  if(leng!=0)
  {
    //Generate a string
     HNGI.getStringMPO2(StringOp,leng,0,false);
     cout << "Creating a string with length " << leng << endl;
     aux.approximate(StringOp,vacuum,D);    
     contractor.optimize(StringOp,vacuum,aux,2*D);
     aux.gaugeCond('R',true);
     mps=MPS(aux);
  }
  else
  {
    mps=MPS(vacuum);
  } 
  E_string=contractor.contract(mps,hamil,mps);
  cout << "Energy gauge invariant string state: " << E_string << endl;
  
  //Now look at the double configurations
  if(leng>1)
  {
    leng_short = ((((int) floor(leng/2.0))%2)==0) ? ((int) floor(leng/2.0))-1 : ((int) floor(leng/2.0));
    
    HNGI.getStringMPO2(StringOp,leng_short,round(N/2.0)-leng_short,false);
    HNGI.getStringMPO2(StringOp2,leng_short,round(N/2.0)+1,false);
    HNGI.getStringMPO2(StringOpdagger,leng_short,round(N/2.0)-leng_short,true);
    HNGI.getStringMPO2(StringOp2dagger,leng_short,round(N/2.0)+1,true);
    
    //Computing the state explicitly would be costly, but an exact contraction is easily possible when I join the operators
    StrMPOs[0]=&StringOp2;
    StrMPOs[1]=&StringOp;
    MPO::join(2,StrMPOs,JoinedString);
    
    AllMPOs[4]=&StringOp;
    AllMPOs[3]=&StringOp2;
    AllMPOs[2]=&hamil;
    AllMPOs[1]=&StringOp2dagger;
    AllMPOs[0]=&StringOpdagger;
    MPO::join(5,AllMPOs,EnergyMPO);
    
    //Keep in mind that after applying a string operator the state is not normalized anymore
    E_string_short = contractor.contract(vacuum,EnergyMPO,vacuum)/contractor.contract2(JoinedString,vacuum);
  }
  if(leng<(N-3))
  {
    leng_long = ((((int) ceil(leng/2.0))%2)==0) ? ((int) ceil(leng/2.0))+1 : ((int) ceil(leng/2.0));
    
    HNGI.getStringMPO2(StringOp,leng_long,round(N/2.0)-leng_long,false);
    HNGI.getStringMPO2(StringOp2,leng_long,round(N/2.0)+1,false);
    HNGI.getStringMPO2(StringOpdagger,leng_long,round(N/2.0)-leng_long,true);
    HNGI.getStringMPO2(StringOp2dagger,leng_long,round(N/2.0)+1,true);
    
    //Computing the state explicitly would be costly, but an exact contraction is easily possible when I join the operators
    StrMPOs[0]=&StringOp2;
    StrMPOs[1]=&StringOp;
    MPO::join(2,StrMPOs,JoinedString);
    
    AllMPOs[4]=&StringOp;
    AllMPOs[3]=&StringOp2;
    AllMPOs[2]=&hamil;
    AllMPOs[1]=&StringOp2dagger;
    AllMPOs[0]=&StringOpdagger;
    MPO::join(5,AllMPOs,EnergyMPO);
    
    //Keep in mind that after applying a string operator the state is not normalized anymore
    E_string_long = contractor.contract(vacuum,EnergyMPO,vacuum)/contractor.contract2(JoinedString,vacuum);
  }  
  cout << "E_string_short: " << E_string_short << endl;
  cout << "E_string_long: " << E_string_long << endl << endl;
  
  /***********************************************************************************
   * 			       Take gauge interaction away			     *
   **********************************************************************************/
  HNGI.removeGaugeInteraction();
  HNGI.updateHMPO();
  
  HNGI.constructInitialMPS(nongauge_vacuum);
  if((gs_non_gauge_filename.length()==0) || !file_exists(gs_non_gauge_filename)) 
  { 
    time(&start_local);
    contractor.findGroundState(hamil,D,&initial_energy,nongauge_vacuum);    
    time(&end_local);
    file.open(TIMEFILE,ios_base::app);
    file << "Time for computing the ground state\t" << difftime(end_local,start_local) << "s " << endl;
    file.close();
    
    nongauge_vacuum.exportMPS("GS_non_gauge.dat");
  }
  else
  {
    cout << "Importing non gauge invariant ground state" << endl;
    nongauge_vacuum.importMPS(gs_non_gauge_filename.c_str());
  }  
  E_nongauge_vaccuum=contractor.contract(nongauge_vacuum,hamil,nongauge_vacuum);
  cout << "Initial energy: " << E_nongauge_vaccuum << endl;
  
  if(leng!=0)
  {
    //Generate a string
    //Keep in mind that removing the U matrices also affects the string operators, so if this one is rebuilt, then it is a different one from the pervious case
    HNGI.getStringMPO2(StringOp,leng,0,false);
    cout << "Creating a string with length " << leng << endl;
    aux.approximate(StringOp,nongauge_vacuum,D);    
    contractor.optimize(StringOp,nongauge_vacuum,aux,2*D);
    aux.gaugeCond('R',true);
    nongauge_mps=MPS(aux);
  }
  else
  {
    nongauge_mps=MPS(nongauge_vacuum);
  }  
  E_nongauge_string=contractor.contract(nongauge_mps,hamil,nongauge_mps);
  cout << "Energy non gauge invariant string state: " << E_nongauge_string << endl;
  
  //Now look at the double configurations
  EnergyMPO.clear();
  if(leng>1)
  { 
    HNGI.getStringMPO2(StringOp,leng_short,round(N/2.0)-leng_short,false);
    HNGI.getStringMPO2(StringOp2,leng_short,round(N/2.0)+1,false);
    HNGI.getStringMPO2(StringOpdagger,leng_short,round(N/2.0)-leng_short,true);
    HNGI.getStringMPO2(StringOp2dagger,leng_short,round(N/2.0)+1,true);
    
    //Computing the state explicitly would be costly, but an exact contraction is easily possible when I join the operators
    StrMPOs[0]=&StringOp2;
    StrMPOs[1]=&StringOp;
    MPO::join(2,StrMPOs,JoinedString);
    
    AllMPOs[4]=&StringOp;
    AllMPOs[3]=&StringOp2;
    AllMPOs[2]=&hamil;
    AllMPOs[1]=&StringOp2dagger;
    AllMPOs[0]=&StringOpdagger;
    MPO::join(5,AllMPOs,EnergyMPO);
    
    //Keep in mind that after applying a string operator the state is not normalized anymore
    E_nongauge_string_short = contractor.contract(nongauge_vacuum,EnergyMPO,nongauge_vacuum)/contractor.contract2(JoinedString,nongauge_vacuum);
  }
  EnergyMPO.clear();
  if(leng<(N-3))
  {    
    HNGI.getStringMPO2(StringOp,leng_long,round(N/2.0)-leng_long,false);
    HNGI.getStringMPO2(StringOp2,leng_long,round(N/2.0)+1,false);
    HNGI.getStringMPO2(StringOpdagger,leng_long,round(N/2.0)-leng_long,true);
    HNGI.getStringMPO2(StringOp2dagger,leng_long,round(N/2.0)+1,true);
    
    //Computing the state explicitly would be costly, but an exact contraction is easily possible when I join the operators
    StrMPOs[0]=&StringOp2;
    StrMPOs[1]=&StringOp;
    MPO::join(2,StrMPOs,JoinedString);
    
    AllMPOs[4]=&StringOp;
    AllMPOs[3]=&StringOp2;
    AllMPOs[2]=&hamil;
    AllMPOs[1]=&StringOp2dagger;
    AllMPOs[0]=&StringOpdagger;
    MPO::join(5,AllMPOs,EnergyMPO);
    
    //Keep in mind that after applying a string operator the state is not normalized anymore
    E_nongauge_string_long = contractor.contract(nongauge_vacuum,EnergyMPO,nongauge_vacuum)/contractor.contract2(JoinedString,nongauge_vacuum);
  }
  
  cout << "E_nongauge_string_short: " << E_nongauge_string_short << endl;
  cout << "E_nongauge_string_long: " << E_nongauge_string_long << endl << endl;
  
  /***********************************************************************************
   * 				Look at some observables			     *
   **********************************************************************************/
  
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
  
  //Save configurations
  configuration=current_configuration(Spin1_MPOs,Spin2_MPOs,Gauge_MPOs,mps,contractor);
  save_configuration(FILENAME_CONFIG_GAUGE,configuration);
  
  configuration=current_configuration(Spin1_MPOs,Spin2_MPOs,Gauge_MPOs,nongauge_mps,contractor);
  save_configuration(FILENAME_CONFIG_NON_GAUGE,configuration);
  
  //Correlation functions
  correlations=getChargeSquareCorrelations(HNGI,contractor,mps,N);
  correlations_nongauge=getChargeSquareCorrelations(HNGI,contractor,nongauge_mps,N);
  file_gauge.open("SU2_correlations_gauge.txt", ios_base::app);
  file_non_gauge.open("SU2_correlations_non_gauge.txt", ios_base::app);
  file_gauge << leng ;
  file_non_gauge << leng ;
  for(int i=0; i<correlations.size(); i++)
  {	
    file_gauge << "\t" << correlations[i];
    file_non_gauge << "\t" << correlations_nongauge[i];
  }
  //cout << endl;
  file_gauge << endl;
  file_non_gauge << endl;
  file_gauge.close();
  file_non_gauge.close();
    
  //Charge squared
  file_gauge.open(CHARGE_GAUGE_FILE, ios_base::app);
  file_non_gauge.open(CHARGE_NON_GAUGE_FILE, ios_base::app);
  file_gauge << leng ;
  file_non_gauge << leng ;
  for(int pos=1; pos<=N; pos++)
  {	
    HNGI.getChargeSquareMPO(charge,pos);
    file_gauge << "\t" << contractor.contract(mps,charge,mps);
    file_non_gauge << "\t" << contractor.contract(nongauge_mps,charge,nongauge_mps);
  }
  //cout << endl;
  file_gauge << endl;
  file_non_gauge << endl;
  file_gauge.close();
  file_non_gauge.close();
  
  Gx_val=contractor.contract(mps,Gx,mps);
  Gy_val=contractor.contract(mps,Gy,mps);
  Gz_val=contractor.contract(mps,Gz,mps);
  
  if(!file_exists(FILENAME_GAUGE))
  {
    file.open(FILENAME_GAUGE);
    file << "#N\tx\tmu\tg\tLeng\tE_vac\tE_string\tQ²\tGx\tGy\tGz\tl_short\tE_short\tl_long\tE_long" << endl;
    file.close();
  }
  file_gauge.open(FILENAME_GAUGE,ios_base::app);
  file_gauge <<N << "\t" << epsilon  << "\t" << mu  << "\t" << g << "\t" << leng << "\t" << E_vacuum << "\t" <<E_string << "\t" <<contractor.contract(mps,charge_square,mps) << "\t" << Gx_val << "\t" << Gy_val << "\t" << Gz_val << "\t" << leng_short << "\t" <<E_string_short << "\t" << leng_long << "\t" << E_string_long << endl ;  
  file_gauge.close();
  
  Gx_val=contractor.contract(nongauge_mps,Gx,nongauge_mps);
  Gy_val=contractor.contract(nongauge_mps,Gy,nongauge_mps);
  Gz_val=contractor.contract(nongauge_mps,Gz,nongauge_mps);
  
  if(!file_exists(FILENAME_NON_GAUGE))
  {
    file.open(FILENAME_NON_GAUGE);
    file << "#N\tx\tmu\tg\tLeng\tE_vac\tE_string\tQ²\tGx\tGy\tGz\tl_short\tE_short\tl_long\tE_long" << endl;
    file.close();
  }
  file_non_gauge.open(FILENAME_NON_GAUGE,ios_base::app);
  file_non_gauge << N << "\t" << epsilon  << "\t" << mu  << "\t" << g << "\t" << leng << "\t" << E_nongauge_vaccuum << "\t" << E_nongauge_string << "\t" << contractor.contract(nongauge_mps,charge_square,nongauge_mps) << "\t" << Gx_val << "\t" << Gy_val << "\t" << Gz_val << "\t" << leng_short << "\t" << E_nongauge_string_short << "\t" << leng_long << "\t" << E_nongauge_string_long << endl;
  file_non_gauge.close();  
  
  
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
void print_MPS(MPS mps)
{
  cout << "MPS of length " << mps. getLength() << endl;
  for(int i=0; i<mps.getLength(); i++)
  {
    cout << "[" << i << "]=" <<  (mps.getA(i)).getA() << endl;
  }
}


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

//All kinds of variants for String coefficients
vector<complex_t> get_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, MPS &vacuum,int N)
{
  MPO mpo(1);
  vector<complex_t> result;
  complex_t tmp;
  
  MPS stringstate;
  
  MPS strong_coupling_test_state;
  H.constructInitialMPS(strong_coupling_test_state);
  
  /*cout << "My vacuum" << endl;
  print_MPS(vacuum);
  
  cout << "Strong coupling vacuum" << endl;
  print_MPS(strong_coupling_test_state);*/
  
  cout << "Overlap <My vacuum|Strong coupling vacuum>: " << abs(contractor.contract(vacuum,strong_coupling_test_state)) << endl;
  

  //Loop over length
  //cout << "Coeffis: (";
  for(int leng=1; leng<=N-1; leng+=2)
  {
    //Loop over position
    for(int pos=1; pos<=N-leng; pos++)
    {
      cout << " --> String at pos " << pos << " of length " << leng << endl;
      /*******************************************************************************
       * Compute <0|S^\dagger|\Psi(t)> which is something like a Wilson loop in real time
       * ****************************************************************************/      
      //Get the adjoint string operator S^\dagger
      //H.getStringMPO2(mpo,leng,pos,true);      
      //Attention: contractor takes argument in the order |ket> MPO <bra|
      //result.push_back(contractor.contract(mps,mpo,vacuum));
      
      /*******************************************************************************
       * Compute <\Psi(t)|SS^\dagger|\Psi>
       * ****************************************************************************/      
      //Get the adjoint string operator S^\dagger
      //H.getStringMPO2(mpo,leng,pos,true);      
      //result.push_back(contractor.contract2(mpo,mps));
      
      /*******************************************************************************
       * Compute <\Psi(t)|Id_A \otimes Tr_A(S|\Omega><\Omega|S^\dagger)) |\Psi(t)>
       * ****************************************************************************/      
      //Warning, if the bond dimension of the mps is large, this will be a huge object which eventually does not fit in the memory anymore
      //Get the adjoint string operator S^\dagger
      H.getStringMPO2(mpo,leng,pos);
      //cout << "Vacuum overlap: " << abs(contractor.contract(vacuum,strong_coupling_test_state)) << endl;
      //cout << "Vacuum norm: " << contractor.contract(vacuum,vacuum) << endl;
      //cout << "<mps|vacuum>: " << abs(contractor.contract(mps,vacuum)) << endl;
      //Get an initial guess for the stringstate
      stringstate.approximate(mpo,vacuum,vacuum.getBond());
      //cout << "Vacuum bond: " << vacuum.getBond() << endl;
      //Now compute the stringstate
      contractor.optimize(mpo,vacuum,stringstate);
      /*cout << "Warning replaced interacting vacuum state with something else" << endl;
      contractor.optimize(mpo,strong_coupling_test_state,stringstate);*/
      mpo.clear();
      H.buildReducedStringMPO(mpo,stringstate,pos,pos+leng);
      result.push_back(contractor.contract(mps,mpo,mps));
      cout << "<mps|stringstate>: " << abs(contractor.contract(mps,stringstate)) << ", Coefficient: " <<  contractor.contract(mps,mpo,mps) << endl;
    }
  }
  /*cout << ")" <<endl;
  cout << "Sum: " << sum << endl;*/
  return result;
}

vector<complex_t> get_StringCoefficients(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N, unsigned int method)
{
  vector<complex_t> result;
  
  MPS strong_coupling_vacuum, op;
  MPO StringOp(1);
  
  H.constructInitialMPS(strong_coupling_vacuum);
  op=MPS(strong_coupling_vacuum);
  op.increaseBondDimension(20);
  op.setRandomState();
  
  //Loop over length
  //int count=1;
  for(int leng=1; leng<=N-1; leng+=2)
  //for(int leng=1; leng<=N-1; leng++)
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
      //cout << "<mps|stringstate>: " << abs(contractor.contract(op,mps)) << ", Coefficient: " << contractor.contract(mps,StringOp,mps) << endl;
      /*******************************************************************************
       * Compute <\Psi(t)|P|\Psi(t)>
       * ****************************************************************************/
      //H.getStringProjector(StringOp,leng,pos,method);
      //result.push_back(contractor.contract(mps,StringOp,mps));
      //cout << "Number: " << count << ", Pos: " << pos << ", Length: " << leng << ", <Psi|P|Psi>=" << contractor.contract(mps,StringOp,mps)<< endl;
      //count++;
    }
  }
  return result;
}

vector<complex_t> getChargeSquareCorrelations(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N)
{
  vector<complex_t> result;
  MPO correlation(1);
  
   //Loop over length
  for(int leng=1; leng<=N-1; leng++)
  {
    //Loop over position
    for(int pos=1; pos<=N-leng; pos++)
    {
      H.getChargeCorrelatorMPO2(correlation,leng,pos,0.5);
      result.push_back(contractor.contract(mps,correlation,mps));
    }
  }
  return result;
}


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
