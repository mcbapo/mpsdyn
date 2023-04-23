/**
   \file ReducesSU2.cpp
   
   
  \param <N> (int) Number of Fermions
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <epsilon> (double) Hopping parameter
  \param <mu> (double) Mass
  \param <g> (double) Coupling constant
  
  \author Stefan Kühn
  \date 15/08/2015

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

#include "ReducedSU2Hamiltonian.h"


#define EPS 1E-8
#define NUM_OF_PARAMS 10
#define RESOLUTION 10
#define FILENAME "SU2red.txt"
#define CONFIGFILE "SU2redconfiguration.txt"
#define CHARGEFILE "SU2redcharge.txt"
#define EVOFILE "SU2redevo.txt"
#define GAUSSFILE "SU2redgausslaw.txt"

using namespace std;
using namespace shrt;

//bool file_exists(string filename);
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
	 << "n:       " << "Number of states on a link" << endl
	 << "steps:   " << "Number of steps for time evolution, a negative number means variational search" << endl
	 << "dt:      " << "Time step size" << endl
	 << "leng:    " << "Length of the inital strong coupling string (leng<=0 means no string)" << endl
	 << "quanta:  " << "How much flux quanta the string should carry" << endl;
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
  int leng=atoi(argv[++cntr]); 
  //The number of flux quanta
  int quanta=atoi(argv[++cntr]); 
  
  cout << "------------------------------------------------------------------------"<<endl;
  cout << "Parameters: " << endl
       << " --> N        = " << N << endl
       << " --> D        = " << D << endl
       << " --> epsilon  = " << epsilon << endl
       << " --> mu       = " << mu << endl
       << " --> g        = " << g << endl
       << " --> n        = " << dlink << endl
       << " --> steps    = " << steps << endl
       << " --> dt       = " << dt << endl
       << " --> leng     = " << leng << endl
       << " --> quanta   = " << quanta << endl;
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

  
  if(steps<0)
  {
    cout << "------------------------------------------------------------------------"<<endl;
    cout << "------------------Starting to compute ground state----------------------"<<endl;
    cout << "------------------------------------------------------------------------"<<endl; 
    cout.flush();
    
    HNGI.constructInitialMPS(mps);
    mps.increaseBondDimension(D);    
    contractor.findGroundState(hamil,D,&E0,mps);
    cout << "Condensate " << contractor.contract(mps,condensate,mps) << endl;
    cout << "Ground state energy " << contractor.contract(mps,hamil,mps) << endl;
    
    result=get_full_config(NocOps,FluxOps,contractor,mps,N);
    save_configuration(CONFIGFILE,result); 
    result.clear();
    result=get_config(ChargeOps,contractor,mps,N);
    save_configuration(CHARGEFILE,result);
    result = get_config(GaussLawOps,contractor,mps,N);
    cout << "Gauss Law=" << result << endl;
    cout << "Gauss Law Penalty " << contractor.contract(mps,penalty,mps) << endl;
    
    if(file_exists(FILENAME))
    {
      file.open(FILENAME, ios_base::app);
      //file << setprecision(16) << N << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << E0 << "\t" << Gx_val << "\t" << Gy_val << "\t" << Gz_val << endl;
      file << setprecision(16) << N << "\t" << D << "\t" << dlink << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << E0 << "\t" << contractor.contract(mps,condensate,mps) << "\t" << contractor.contract(mps,penalty,mps) << endl;
    }
    else
    {
      file.open(FILENAME, ios_base::app);
      file << setprecision(16) << "#N\tD\tdlink\tepsilon\tmu\tg\tE0\tC\tP" << endl;
      //file << setprecision(16) << N << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << E0 << "\t" << Gx_val << "\t" << Gy_val << "\t" << Gz_val << endl;
      file << setprecision(16) << N << "\t" << D << "\t" << dlink << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << E0 << "\t" << contractor.contract(mps,condensate,mps) << "\t" << contractor.contract(mps,penalty,mps) << endl;
    }
    file.close();
    }
  else
  {
    cout << "------------------------------------------------------------------------"<<endl;
    cout << "----------------Starting to compute time evolution----------------------"<<endl;
    cout << "------------------------------------------------------------------------"<<endl; 
    cout.flush();
    
    MPO Uodd(1),Uodd_half(1),Ueven(1),Umpo(1),Podd(1),Peven(1),Pgauss(1);
    bool imag_evo=true;   
    
    HNGI.getUoddMPO(Uodd,dt,imag_evo);
    HNGI.getUoddMPO(Uodd_half,dt/2.0,imag_evo);
    HNGI.getUevenMPO(Ueven,dt,imag_evo);
    HNGI.getUMPO(Umpo,dt,imag_evo);   
    
    if(leng<=0)
    {
      HNGI.constructInitialMPS(mps);
      //HNGI.getProjectorMPO2(Podd,Peven);
    }
    else
    {
      HNGI.getHeavyString(mps,leng,quanta);
      //HNGI.getProjectorMPO(Podd,Peven,leng,quanta);
    }
    /*HNGI.getProjectorMPO(Podd,Peven);
    HNGI.getProjectorMPO2(Podd,Peven);
    HNGI.getProjectorMPO(Podd,Peven,leng,quanta);*/
    HNGI.getProjectorMPO2(Podd,Peven,leng,quanta);
    HNGI.getProjectorSumMPO(Pgauss);
    
    result = get_full_config(NocOps,FluxOps,contractor,mps,N);
    save_configuration(CONFIGFILE,result); 
    
    mps.increaseBondDimension(D);    
    if(!file_exists(EVOFILE))
    {
      file.open(EVOFILE, ios_base::app);
      file << "#Step\tD\tepsilon\tmu\tg\tE\tC\tP\tNorm" << endl;
      file.close();
    }
    result = get_full_config(NocOps,FluxOps,contractor,mps,N);
    //print_config(result);
    save_configuration(CONFIGFILE,result); 
    result = get_config(ChargeOps,contractor,mps,N);
    save_configuration(CHARGEFILE,result); 
    //cout << "Charge=" << result << endl;
    result = get_config(LvalOps,contractor,mps,N-1);
    cout << "Lval=" << result << endl;      
    result = get_config(GaussLawOps,contractor,mps,N);
    save_configuration(GAUSSFILE,result);
    cout << "Gauss Law=" << result << endl;
    cout << "Gauss Law Penalty " << contractor.contract(mps,penalty,mps) << endl;
    cout << "Projector sum " << contractor.contract(mps,Pgauss,mps) << endl;
      
    file.open(EVOFILE, ios_base::app);
    file << 0 << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << setprecision(10) << contractor.contract(mps,hamil,mps) << "\t" << contractor.contract(mps,condensate,mps) << "\t" << contractor.contract(mps,penalty,mps) << "\t" << contractor.contract(mps,mps) <<endl;
    file.close();  
    cout << "Energy after 0 steps: " << contractor.contract(mps,hamil,mps) << endl;
    cout << "Condensate after 0 steps: " << contractor.contract(mps,condensate,mps) << endl;
    cout << "Norm after 0 steps: " << contractor.contract(mps,mps) << endl << endl;
    for(int i=1; i<=(steps/RESOLUTION); i++)
    {
      compute_evolution(Ueven,Uodd,Uodd_half,contractor,mps,RESOLUTION,imag_evo);
      //compute_evolution(Umpo,contractor,mps,RESOLUTION,imag_evo);
      project(Podd,Peven,contractor,mps,false);
      
      result = get_full_config(NocOps,FluxOps,contractor,mps,N);
      save_configuration(CONFIGFILE,result);
      
      result = get_config(ChargeOps,contractor,mps,N);
      save_configuration(CHARGEFILE,result);       
      //cout << "Charge=" << result << endl;
      
      result = get_config(LvalOps,contractor,mps,N-1);
      cout << "Lval=" << result << endl;
      
      result = get_config(GaussLawOps,contractor,mps,N);
      save_configuration(GAUSSFILE,result);
      //cout << "Gauss Law=" << result << endl;
      cout << "Gauss Law Penalty " << contractor.contract(mps,penalty,mps) << endl;
      cout << "Projector sum " << contractor.contract(mps,Pgauss,mps) << endl;
      
      file.open(EVOFILE, ios_base::app);
      file << i*RESOLUTION << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << setprecision(10) << contractor.contract(mps,hamil,mps) << "\t" << contractor.contract(mps,condensate,mps) << "\t" << contractor.contract(mps,penalty,mps) << "\t" << contractor.contract(mps,mps) <<endl;
      file.close();
      cout << "Energy after " << i*RESOLUTION << " steps: " << contractor.contract(mps,hamil,mps) << endl;
      cout << "Condensate after " << i*RESOLUTION << " steps: " << contractor.contract(mps,condensate,mps) << endl;
      cout << "Norm after " << i*RESOLUTION << " steps: " << contractor.contract(mps,mps) << endl << endl;
    }  
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
    
