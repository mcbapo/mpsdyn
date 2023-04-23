/**
   \file ReducedSU2Stringstates.cpp
   
   
  \param <N> (int) Number of Fermions
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <epsilon> (double) Hopping parameter
  \param <mu> (double) Mass
  \param <g> (double) Coupling constant
  
  \author Stefan Kühn
  \date 01/09/2015
  
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
#include "MPO.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include "Indices.h"
#include "misc.h"
#include "LGTstring.h"

#include "ReducedSU2Hamiltonian.h"

#define EPS 1E-8
#define NUM_OF_PARAMS 10
#define RESOLUTION 10
#define CONFIGFILE "SU2redconfiguration.txt"
#define CHARGEFILE "SU2redcharge.txt"
#define EVOFILE "SU2redevo.txt"
#define GAUSSFILE "SU2redgausslaw.txt"
#define ENTROPYFILE "SU2redentropy.txt"

using namespace std;
using namespace shrt;

/** Little help function, given a string of maximum length determine all configurations with strings on the state which are possible*/
void get_all_strings_helper(int start, int end, int leng, vector<string_t> tmp, vector< vector<string_t> > &res);
/** Given a start site and and end site and a maximum length, determine all possible string configurations on the strong coupling vacuum containing strings in the specified region with length at most the given maximum length.*/
vector< vector<string_t> > get_all_strings(int start, int end, int leng);


void print_MPS(MPS);
//Get spin and flux configuration
vector<complex_t> get_full_config(vector<MPO*> NocOps, vector<MPO*> FluxOps, Contractor &contractor, MPS &mps, int N);
//Get charge configuration
vector<complex_t> get_config(vector<MPO*> Ops, Contractor &contractor, MPS &mps, int N);
//Save the current configuration to a given file
int save_configuration(string filename, vector<complex_t> config);
//Time Evolution
void compute_evolution(MPO &exp_even,MPO &exp_odd,MPO &exp_odd_half,Contractor &contractor, MPS &mps, int steps,bool normalize);
void compute_evolution(MPO &Umpo,Contractor &contractor, MPS &mps, int steps,bool normalize);
void print_config(vector<complex_t> confi);
//Compute the sum of the absolute values of the vector entries
double sum_abs(vector<complex_t> input);
//Given the projectors, project back into the Gauss Law fulfilling subspace
void project(MPO &Podd, MPO &Peven, Contractor &contractor, MPS &mps, bool normalize=false);
//Given the projector, project back into the Gauss Law fulfilling subspace
void project(MPO &Projector, Contractor &contractor, MPS &mps, bool normalize=false);

//Read all string configurations from a given file
void read_string_configs(string filename, vector<string_t> &configs);

int main(int argc,const char* argv[]){
  int cntr=0;
  
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1)
  {
    cout << "Program compiled on " << __DATE__ << ", " << __TIME__ << endl;
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage:   " << argv[0] << " N D µ g epsilon" << endl;
    cout << "N:       " << "Number of spins" << endl
	 << "D:       " << "Bond dimension of the MPS" << endl
	 << "epsilon: " << "Hopping parameter" << endl
	 << "µ:       " << "Mass" << endl
	 << "g:       " << "Coupling constant" << endl
	 << "n:       " << "Number of states on a link" << endl
	 << "steps:   " << "Number of steps for time evolution, a negative number means variational search" << endl
	 << "dt:      " << "Time step size" << endl
	 << "file:    " << "File containing string configurations" << endl
	 << "quanta:  " << "How much flux quanta the string should carry" << endl
	 << "itime:   " << "Evolution in imaginary time, default is yes (optional)" << endl
	 << "file:    " << "Name of the file for an existing state to be further evolved (optional)" << endl
	 << "offset:  " << "Offset to be used to evolve an existing state (optional)" << endl;
    return -1;
  }
  
  //Number of spin sites
  int N=atoi(argv[++cntr]);
  //Bond dimension
  int D=atoi(argv[++cntr]);
  //Hopping parameter
  double epsilon=atof(argv[++cntr]);
  //Mass
  double mu=atof(argv[++cntr]);
  //Coupling constant
  double g=atof(argv[++cntr]);
  //Number of excited states to be computed
  int dlink=atoi(argv[++cntr]);
  //Number of time steps to be computed
  int steps=atoi(argv[++cntr]);
  //Number of time steps to be computed
  double dt=atof(argv[++cntr]);
  //Length of the string in case time evolution is computed
  string config_file=argv[++cntr]; 
  //The number of flux quanta
  int quanta=atoi(argv[++cntr]); 
  //optional arguments
  bool imag_evo=true; 
  if(argc > NUM_OF_PARAMS + 1)
    imag_evo = atoi(argv[++cntr]);
  string init_state_file="";
  if(argc > NUM_OF_PARAMS + 2)
    init_state_file=argv[++cntr];
  int offset=0;
  if(argc > NUM_OF_PARAMS + 3)
    offset=atoi(argv[++cntr]);  
  
  //Entanglement entropy
  double entropy;
  //Revision
  string revision="$Rev: 1580 $";
  //Parse revision string
  revision=revision.substr(6,4);
  
  
  cout << "------------------------------------------------------------------------"<<endl;
  cout << "Generated with code revision " << revision << endl;
  cout << "Parameters: " << endl
       << " --> N        = " << N << endl
       << " --> D        = " << D << endl
       << " --> epsilon  = " << epsilon << endl
       << " --> mu       = " << mu << endl
       << " --> g        = " << g << endl
       << " --> n        = " << dlink << endl
       << " --> steps    = " << steps << endl
       << " --> dt       = " << dt << endl
       << " --> file     = " << config_file << endl
       << " --> quanta   = " << quanta << endl
       << " --> itime    = " << imag_evo << endl
       << " --> file     = " << init_state_file << endl
       << " --> offset   = " << offset << endl;
  cout << "------------------------------------------------------------------------"<<endl;
       
  
  //File for results
  ofstream file; 
  //Vector to save expectation values of operators
  vector<complex_t> result;
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  contractor.setConvTol(EPS);

  //Get the Hamiltonian
  ReducedSU2Hamiltonian HNGI(N,epsilon,mu,g,dlink);
  //ReducedSU2Hamiltonian HNGI(N,epsilon,mu,g,dlink,1000.0);
  const MPO& hamil=HNGI.getHMPO();    
  
  cout << "Successfully retrived ReducedSU2Hamiltonian" << endl;
  
  //Now prepare and MPS
  vector<int> dims(2*N-1,3);
  for(int i=1; i<(2*N-1); i+=2) 
    dims[i]=dlink;
  MPS mps(2*N-1,D,dims);
  double E0=0.0;
  
  
  //Build the MPOs for the observables
  vector<MPO*> FluxOps,NocOps,ChargeOps,GaussLawOps,LvalOps;
  MPO condensate(1),penalty(1);
  MPO* tmp_mpo;
  for (int i=1;i<=N; i++)
  {
    tmp_mpo = new MPO(1);
    HNGI.getNocMPO(*tmp_mpo,i);
    NocOps.push_back(tmp_mpo);
    
    tmp_mpo = new MPO(1);
    HNGI.getChargeSquareMPO(*tmp_mpo,i);
    ChargeOps.push_back(tmp_mpo);
    
    tmp_mpo = new MPO(1);
    HNGI.getGaussLawMPO(*tmp_mpo,i);
    GaussLawOps.push_back(tmp_mpo);
    if(i<N)
    {
      tmp_mpo = new MPO(1);
      HNGI.getJsquaredMPO(*tmp_mpo,i);
      FluxOps.push_back(tmp_mpo);
      
      tmp_mpo = new MPO(1);
      HNGI.getLMPO(*tmp_mpo,i);
      LvalOps.push_back(tmp_mpo);
    }
  } 
  HNGI.getCondensateMPO(condensate,true);
  HNGI.getGaussLawPenaltyMPO(penalty,1.0); 
  
  cout << "------------------------------------------------------------------------"<<endl;
  cout << "----------------Starting to compute time evolution----------------------"<<endl;
  cout << "------------------------------------------------------------------------"<<endl; 
  cout.flush();
  
  MPO Uodd(1),Uodd_half(1),Ueven(1),Umpo(1),Podd(1),Peven(1),Pgauss(1);
    
  
  HNGI.getUoddMPO(Uodd,dt,imag_evo);
  HNGI.getUoddMPO(Uodd_half,dt/2.0,imag_evo);
  HNGI.getUevenMPO(Ueven,dt,imag_evo);
  HNGI.getUMPO(Umpo,dt,imag_evo);   
  HNGI.getProjectorSumMPO(Pgauss);
  
  if(!file_exists(config_file))
  {
    if(file_exists(init_state_file.c_str()))
       mps.importMPS(init_state_file.c_str());
    else
    {
      cout << "No valid file for string configuration given, will start with strong coupling state" << endl;
      HNGI.constructInitialMPS(mps);      
    }
    HNGI.getProjectorMPO2(Podd,Peven);
  }
  else
  {
    //Old an deprecated
    /*HNGI.getHeavyString(mps,leng,quanta);
    HNGI.getProjectorMPO2(Podd,Peven,leng,quanta);*/
    
    //Read string configs from file and compute projectors as well as the string state accordingly
    vector<string_t> configs;  
    read_string_configs(config_file,configs);      
    HNGI.getProjectorMPO2(Podd,Peven,configs,quanta);
    
    if(file_exists(init_state_file.c_str()))
    {  
      cout << "Reading existing state " << init_state_file << endl;
      mps.importMPS(init_state_file.c_str());
    }
    else
    {
      cout << "Starting with strong coupling string " << endl;
      HNGI.constructStringMPS(mps,configs);
    }    
  }    
      
  if(!file_exists(EVOFILE))
  {
    file.open(EVOFILE, ios_base::app);
    file << "#Created with code revision " << revision << endl;
    file << "#Step\tD\tepsilon\tmu\tg\tE\tC\tP\tNorm" << endl;
    file.close();
  }
  else
  {
    file.open(EVOFILE, ios_base::app);
    file << "#Created with code revision " << revision << endl;
    file.close();
  }  
   
  //If I start a fresh evolution, save configuration once before I actually start doing anything
  if(offset==0)
  {
    result = get_full_config(NocOps,FluxOps,contractor,mps,N);
    save_configuration(CONFIGFILE,result); 
    result = get_config(ChargeOps,contractor,mps,N);
    save_configuration(CHARGEFILE,result); 
    result = get_config(NocOps,contractor,mps,N-1);
    cout << "Noc=" << result << endl; 
    result = get_config(LvalOps,contractor,mps,N-1);
    cout << "Lval=" << result << endl;      
    result = get_config(GaussLawOps,contractor,mps,N);
    save_configuration(GAUSSFILE,result);
    cout << "Gauss Law=" << result << endl;
    cout << "Gauss Law Penalty " << contractor.contract(mps,penalty,mps) << endl;
    cout << "Projector sum " << contractor.contract(mps,Pgauss,mps) << endl;
    result = get_full_config(NocOps,FluxOps,contractor,mps,N);
    save_configuration(CONFIGFILE,result); 
    
    file.open(EVOFILE, ios_base::app);
    file << 0 << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << setprecision(10) << contractor.contract(mps,hamil,mps) << "\t" << contractor.contract(mps,condensate,mps) << "\t" << contractor.contract(mps,penalty,mps) << "\t" << contractor.contract(mps,mps) <<endl;
    file.close();  
    cout << "Energy after 0 steps: " << contractor.contract(mps,hamil,mps) << endl;
    cout << "Condensate after 0 steps: " << contractor.contract(mps,condensate,mps) << endl;
    cout << "Norm after 0 steps: " << contractor.contract(mps,mps) << endl << endl;
  }
  
  //Adjust bond dimension 
  mps.increaseBondDimension(D);
  
  //Now compute the evolution
  for(int i=1+(offset/RESOLUTION); i<=((steps+offset)/RESOLUTION); i++)
  {
    compute_evolution(Ueven,Uodd,Uodd_half,contractor,mps,RESOLUTION,imag_evo);
    //compute_evolution(Umpo,contractor,mps,RESOLUTION,imag_evo);
    //cout << "Projector turned off" << endl;
    project(Podd,Peven,contractor,mps,false);
    
    result = get_full_config(NocOps,FluxOps,contractor,mps,N);
    save_configuration(CONFIGFILE,result);
    
    result = get_config(ChargeOps,contractor,mps,N);
    save_configuration(CHARGEFILE,result);       
    //cout << "Charge=" << result << endl;
    
    result = get_config(LvalOps,contractor,mps,N-1);
    //cout << "Lval=" << result << endl;
    
    result = get_config(GaussLawOps,contractor,mps,N);
    save_configuration(GAUSSFILE,result);
    
    entropy = contractor.getEntropy(mps,0);
    file.open(ENTROPYFILE, ios_base::app);
    file <<  i*RESOLUTION << "\t" << entropy << endl;
    file.close();
    
    //cout << "Gauss Law=" << result << endl;
    cout << "Gauss Law Penalty:          " << contractor.contract(mps,penalty,mps) << endl;
    cout << "Projector sum:              " << contractor.contract(mps,Pgauss,mps) << endl;
    cout << "Center entropy:             " << entropy << endl;
    
    file.open(EVOFILE, ios_base::app);
    file << i*RESOLUTION << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << setprecision(10) << contractor.contract(mps,hamil,mps) << "\t" << contractor.contract(mps,condensate,mps) << "\t" << contractor.contract(mps,penalty,mps) << "\t" << contractor.contract(mps,mps) <<endl;
    file.close();
    cout << "Energy after " << i*RESOLUTION << " steps:     " << contractor.contract(mps,hamil,mps) << endl;
    cout << "Condensate after " << i*RESOLUTION << " steps: " << contractor.contract(mps,condensate,mps) << endl;
    cout << "Norm after " << i*RESOLUTION << " steps:        " << contractor.contract(mps,mps) << endl << endl;
  }
  
  //Save the (ground) state
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
  
  
  //Free memory
  for(int k=0; k<N; k++)
  {
    NocOps[k]->clear();
    delete NocOps[k];
    
    ChargeOps[k]->clear();
    delete ChargeOps[k];
    
    GaussLawOps[k]->clear();
    delete GaussLawOps[k];
    if(k<(N-1))
    {
      FluxOps[k]->clear();
      delete FluxOps[k];
      
      LvalOps[k]->clear();
      delete LvalOps[k];
    }
  }
  
  cout << "End of program" << endl;
  return 0;
  
}

void get_all_strings_helper(int start, int end, int leng, vector<string_t> tmp, vector< vector<string_t> > &res)
{
  int curr_leng=0;
  //Local copy of current path
  vector<string_t> local_tmp;
  string_t curr_str;
  
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


vector< vector<string_t> > get_all_strings(int start, int end, int leng)
{
  //Vector to keep the result
  vector< vector<string_t> > res;
  //Vector to keep track of the current path
  vector<string_t> tmp;
  //Now call the helper
  get_all_strings_helper(start,end,leng,tmp,res);
  return res;
}

void print_MPS(MPS mps)
{
  cout << "MPS of length " << mps. getLength() << endl;
  for(int i=0; i<mps.getLength(); i++)
  {
    cout << "[" << i << "]=" <<  (mps.getA(i)).getA() << endl;
  }
}

vector<complex_t> get_full_config(vector<MPO*> NocOps, vector<MPO*> FluxOps, Contractor &contractor, MPS &mps, int N)
{
  vector <complex_t> res;
  res.resize(2*N-1);
  
  #pragma omp prallel for
  for(int i=1; i<=N; i++)
  {
    res[2*i-2]=contractor.contract(mps,*NocOps[i-1],mps);
    if(i<N)
      res[2*i-1]=contractor.contract(mps,*FluxOps[i-1],mps);    
  }    
  return res;
}

vector<complex_t> get_config(vector<MPO*> Ops, Contractor &contractor, MPS &mps, int N)
{
  vector <complex_t> res;
  res.resize(N);
  
  #pragma omp prallel for
  for(int i=1; i<=N; i++)
    res[i-1]=contractor.contract(mps,*Ops[i-1],mps);
    
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


//Little function to compute the trotter splitting, as it is quite complicated and I want to take two time steps together to improve performance
void compute_evolution(MPO &exp_even,MPO &exp_odd,MPO &exp_odd_half,Contractor &contractor, MPS &mps, int steps,bool normalize)
{
   MPS aux=MPS(mps);
   
  //Now with a bit more optimization to save time  
  //Inital half step
  contractor.optimize(exp_odd_half,mps,aux);
  for(int i=1; i<=(steps-1); i++)
  {
    //cout << "\rStep " << i <<"      " ;
    cout.flush();   
    mps.setRandomState();
    contractor.optimize(exp_even,aux,mps);
    aux.setRandomState();
    contractor.optimize(exp_odd,mps,aux);
    //cout << endl << "Norm=" << contractor.contract(mps,mps) << endl;
    if(normalize)
    {
      //cout << "Normalizing" << endl;
      mps.gaugeCond('R',true);      
    }
  }
  
  //Complete last step
  contractor.optimize(exp_even,aux,mps);
  contractor.optimize(exp_odd_half,mps,aux);  
  mps=aux;
  if(normalize)
  {
    //cout << "Normalizing" << endl;
    mps.gaugeCond('R',true);      
  }  
  //cout << endl;
}

void compute_evolution(MPO &Umpo,Contractor &contractor, MPS &mps, int steps,bool normalize)
{
  MPS aux=MPS(mps);
  for(int i=1; i<=steps; i++)
  {
    contractor.optimize(Umpo,mps,aux);
    mps=aux;
    if(normalize)
      mps.gaugeCond('R',true);
  }
}

void project(MPO &Podd, MPO &Peven, Contractor &contractor, MPS &mps,bool normalize)
{
  MPS aux=MPS(mps);
  contractor.optimize(Podd,mps,aux);
  contractor.optimize(Peven,aux,mps);
  if(normalize)
    mps.gaugeCond('R',true);
}

void project(MPO &Projector, Contractor &contractor, MPS &mps,bool normalize)
{
  MPS aux=MPS(mps);
  contractor.optimize(Projector,mps,aux);
  mps = aux;
  if(normalize)
    mps.gaugeCond('R',true);
}


void print_config(vector<complex_t> confi)
{
  cout << "Config=(";
  int count=0;
  for(int i=0; i<confi.size(); i++)
  {
    count ++;
    if(count !=2 )
      cout << confi[i] << ",";
    else
    {
      cout << confi[i] << "   ";
      count = 0;
    }
  }
  cout << ")" << endl;
}

double sum_abs(vector<complex_t> input)
{
  double result=0.0;
  for(int i=0; i<input.size(); i++)
    result += + abs(input[i]);
  return result;
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

