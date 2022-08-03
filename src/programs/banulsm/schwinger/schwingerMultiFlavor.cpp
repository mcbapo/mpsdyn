/**
   \file schwingerMultiFlavor.cpp
   Multi Flavor Schwinger model
   
  \param <file> (string) Name of the driver file containing the parameters
  
  \author Stefan Kühn
  \date 30/06/2015

*/


#include <cstdlib>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <string>
#include <algorithm>
#include <vector>
//Lib which is needed for sleep command on *nix Systems
//#include "unistd.h"

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include "Indices.h"
#include "misc.h"

#include "SchwingerHamiltonian.h"
#include "MultiFlavorSchwingerHamiltonian.h"

#define NUM_OF_PARAMS 1
#define FILENAME "MultiFlavor.txt"
#define TMP_FILE "GS.tmp"


using namespace std;
using namespace shrt;

vector<complex_t> get_configuration(MultiFlavorSchwingerHamiltonian &H, Contractor &contractor, MPS &mps, int N, int F);
vector<complex_t> get_charge(vector<complex_t> spinconfiguration, int N, int F);
vector<complex_t> get_condensate(vector<complex_t> spinconfiguration, int N, int F);
void save_configuration(vector<complex_t> configuration, string filename);
void read_parameters(string filename,int &N, int &F, int &D, double &x, vector<double> &mu, vector<double> &nu, double &lambda, double &l0, int &bu_interval, string &init_file, double &tol, int &Npen, double &eta);
vector<double> process_data(string str, int F);

int main(int argc,const char* argv[])
{  
  //Code Revision (magic field, automatically filled in by subversion)
  string revision="$Rev: 1522 $";
  //Parse revision string
  revision=revision.substr(6,4);
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1)
  {
    cout << "Program compiled on " << __DATE__ << ", " << __TIME__ << " with code revision " << revision << endl;
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage:   " << argv[0] << " <driver_file>" << endl;
    return -1;
  }
  
  //Read input file
  string inputfile(argv[1]);  
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
    
  //Variables
  vector<double> mu,nu;
  vector<complex_t> particle_expectationvals;
  double E0;
  int N;
  int D;
  int F;
  double x,lambda,l0;
  ofstream file;
  string init_file;
  int bu_interval=50;
  double tol;
  int Npen;
  double eta=0.0;
  
  //Read input parameters from the input file
  read_parameters(inputfile,N,F,D,x,mu,nu,lambda,l0,bu_interval,init_file,tol,Npen,eta);
  
  //Set the convergence tolerance
  contractor.setConvTol(tol);
  
  //MPS for the ground state
  MPS mps;
  
  //Get the Hamiltonian
  MultiFlavorSchwingerHamiltonian *Hmulti;
  MPO Npenalty(1);
  if (eta==0.0)
    Hmulti = new MultiFlavorSchwingerHamiltonian(N,F,x,mu,nu,lambda,l0);
  else
  {
    Hmulti = new MultiFlavorSchwingerHamiltonian(N,F,x,mu,nu,lambda,l0,0,Npen,eta);
    Hmulti->getParticlePenaltyMPO(Npenalty,0,Npen);
  }
  const MPO& hamil=Hmulti->getHMPO();
  
  cout << "Successfully retrieved multi flavor Schwinger Hamiltonian" << endl;	cout.flush();   

  
  //Set up all kins of MPOs
  MPO charge(1),penalty(1);
  vector<MPO*> condensates,particles;
  MPO *condensate_flavor,*particles_flavor;
  Hmulti->getChargeMPO(charge);
  Hmulti->getChargePenaltyMPO(penalty,true);
  for(int i=0; i<F ; i++)
  {
    condensate_flavor=new MPO(1);
    particles_flavor=new MPO(1);
    Hmulti->getCondensateMPO(*condensate_flavor,i,true);
    Hmulti->getParticleNumberMPO(*particles_flavor,i);
    condensates.push_back(condensate_flavor);      
    particles.push_back(particles_flavor);
  }
  condensate_flavor=new MPO(1);
  particles_flavor=new MPO(1);
  Hmulti->getCondensateMPO(*condensate_flavor,-1,true);
  Hmulti->getParticleNumberMPO(*particles_flavor,-1);
  condensates.push_back(condensate_flavor);
  particles.push_back(particles_flavor);
  
  //Now prepare the initial state, see first if temporary results exist or if a file for an initial guess is specified, if nothing is there, start with a strong coupling ground state
  if(file_exists(TMP_FILE))
  {
    mps.importMPS(TMP_FILE);    
    cout << "Successfully loaded temporary result" << endl;
  }
  else if(file_exists(init_file))
  {
    mps.importMPS(init_file.c_str());
    cout << "Successfully loaded initial MPS" << endl;
  }
  else
    Hmulti->constructInitialMPS(mps);
  
  cout << "---------------------------------------------------------------------"<<endl;
  cout << "---------------Starting to compute ground state----------------------"<<endl;
  cout << "---------------------------------------------------------------------"<<endl; 
  cout.flush();
  contractor.findGroundState(hamil,D,&E0,mps,0.0,1,TMP_FILE,bu_interval);  
  mps.exportMPS("GS.dat");
  //Compute some observables 
  for (int i=0; i<=F; i++)
    particle_expectationvals.push_back(contractor.contract(mps,*particles[i],mps));
  cout << "Ground state energy: " << E0 << endl;
  cout << "Condensate:          " << contractor.contract(mps,*condensates[F],mps) << endl;
  cout << "Total charge:        " << contractor.contract(mps,charge,mps) << endl;
  cout << "Penalty energy:      " << contractor.contract(mps,penalty,mps) << endl;
  
  
  if(file_exists(FILENAME))
  {
    file.open(FILENAME,ios_base::app);
    file << N << "\t" << D << "\t" << F << "\t" << x << "\t";
    for(int i=0; i<mu.size(); i++)
      file << setprecision(15) << mu[i] << "\t";
    for(int i=0; i<nu.size(); i++)
      file << setprecision(15) << nu[i] << "\t";
    file << setprecision(15) << E0 << "\t";
    file << setprecision(15) << real(contractor.contract2(hamil,mps)) - E0*E0 << "\t";
    for(int i=0; i<=F; i++)
      file << setprecision(15) << contractor.contract(mps,*condensates[i],mps) << "\t";
    for(int i=0; i<=F; i++)
      file << setprecision(15) << particle_expectationvals[i] << "\t";
    for(int i=0; i<F; i++)
    file << setprecision(15) << real(contractor.contract2(*particles[i],mps)) - real(particle_expectationvals[i])*real(particle_expectationvals[i])  << "\t";
    file << setprecision(15) << contractor.contract(mps,charge,mps) << "\t" << contractor.contract(mps,penalty,mps);
    if(eta!=0.0)
      file << "\t" << contractor.contract(mps,Npenalty,mps);
    
    file << endl;
  }
  else
  {
    file.open(FILENAME,ios_base::app);
    file <<"#N\tD\tF\tx\tmu[0]..mu[F-1]\tnu[0]..nu[F-1]\tE0\tVar(E0)\tC[0]..C[F-1]\tC\tN[0]..N[F-1]\tNtot\tVar(N[0])..Var(N[F-1])\tC_0\tP\n";
    file << N << "\t" << D << "\t" << F << "\t" << x << "\t";
    for(int i=0; i<mu.size(); i++)
      file << setprecision(15) << mu[i] << "\t";
    for(int i=0; i<nu.size(); i++)
      file << setprecision(15) << nu[i] << "\t";
    file << setprecision(15) << E0 << "\t";
    file << setprecision(15) << real(contractor.contract2(hamil,mps)) - E0*E0 << "\t";
    for(int i=0; i<=F; i++)
      file << setprecision(15) << contractor.contract(mps,*condensates[i],mps) << "\t";
    for(int i=0; i<=F; i++)
      file << setprecision(15) << particle_expectationvals[i] << "\t";
    for(int i=0; i<F; i++)
      file << setprecision(15) << real(contractor.contract2(*particles[i],mps)) - real(particle_expectationvals[i])*real(particle_expectationvals[i])  << "\t";
    file << setprecision(15) << contractor.contract(mps,charge,mps) << "\t" << contractor.contract(mps,penalty,mps);
    if(eta!=0.0)
      file << "\t" << contractor.contract(mps,Npenalty,mps);
    file << endl;
  }
  file.close();
  
  //Compute the ground state configuration
  vector <complex_t> spin_vals,condensate_vals,charge_vals;
  spin_vals=get_configuration(*Hmulti,contractor,mps,N,F);
  save_configuration(spin_vals,"MultiFlavorConfiguration.txt");
  charge_vals=get_charge(spin_vals,N,F);
  save_configuration(charge_vals,"MultiFlavorCharge.txt");
  condensate_vals=get_condensate(spin_vals,N,F);  
  save_configuration(condensate_vals,"MultiFlavorCondensate.txt");
  
  //Free memory
  for(int i=0; i<condensates.size(); i++)
  {
    condensates[i]->clear();    
    particles[i]->clear();
    delete condensates[i];
    delete particles[i];
  }   
  delete Hmulti;

  cout << "End of program" << endl;
  return 0;
  
}

vector<complex_t> get_configuration(MultiFlavorSchwingerHamiltonian &H, Contractor &contractor, MPS &mps, int N, int F)
{
  MPO mpo(1);
  vector<complex_t> result;
  for(int n=0; n<N; n++)
    for(int f=0; f<F; f++)
    {
      H.getSpinMPO(mpo,n,f);
      result.push_back(contractor.contract(mps,mpo,mps));
    }
  return result;
}

vector<complex_t> get_charge(vector<complex_t> spinconfiguration, int N, int F)
{
  vector<complex_t> result;
  for(int n=0; n<N; n++)
    for(int f=0; f<F; f++)      
      result.push_back(0.5*(pow(-1.0,n)*ONE_c+spinconfiguration[n*F+f]));
  return result;
}

vector<complex_t> get_condensate(vector<complex_t> spinconfiguration, int N, int F)
{
  vector<complex_t> result;
  for(int n=0; n<N; n++)
    for(int f=0; f<F; f++)      
      result.push_back(0.5*pow(-1.0,n)*(ONE_c+spinconfiguration[n*F+f]));
  return result;
}

void save_configuration(vector<complex_t> configuration, string filename)
{
  ofstream file;
  file.open(filename.c_str(),ios_base::app);
  
  for(int i=0; i<configuration.size(); i++)
  {
    if(i<(configuration.size()-1))
      file << configuration[i] << "\t";
    else
      file << configuration[i] << endl;  
  }
  file.close();
}

void read_parameters(string filename,int &N, int &F, int &D, double &x, vector<double> &mu, vector<double> &nu,double &lambda, double &l0, int &bu_interval,string &init_file, double &tol, int &Npen, double &eta)
{
  if(!file_exists(filename))
  {
    cout << "Error, given input file " << filename << " does not exist" << endl;
    exit(1);
  }
  
  int status=0;
  ifstream file;
  file.open(filename.c_str());
  string dummy,mus,nus;
  
  file >> dummy >> N;
  cout << "Read N = " << N << endl;
  if (!file.good())
    status=1;
  file >> dummy >> F;
  cout << "Read F = " << F << endl;
  if (!file.good())
    status=1;
  file >> dummy >> D;
  cout << "Read D = " << D << endl;
  if (!file.good())
    status=1;
  file >> dummy >> x;
  cout << "Read x = " << x << endl;
  if (!file.good())
    status=1;
  file >> dummy >> mus;
  mu=process_data(mus,F);
  cout << "Read mu =" << mu << endl;
  if (!file.good())
    status=1;
  file >> dummy >> nus;
  nu=process_data(nus,F);
  cout << "Read nu =" << nu << endl;
  if (!file.good())
    status=1;
  file >> dummy >> lambda;
  cout << "Read la = " << lambda << endl;
  if (!file.good())
    status=1;
  file >> dummy >> l0;
  cout << "Read l_0 = " << l0 << endl;
  if (!file.good())
    status=1;
  init_file="";
  file >> dummy >> init_file;
  cout << "Read init_file = " << init_file << endl;
  if (!file.good())
    status=1;
  file >> dummy >> bu_interval;
  cout << "Read bu_interval = " << bu_interval << endl;
  if (!file.good())
    status=1;
  file >> dummy >> tol;
  cout << "Read tol = " << tol << endl;
  if (file.good())
  {
    file >> dummy >> Npen;
    cout << "Read Npen = " << Npen << endl;
    file >> dummy >> eta;
    cout << "Read eta = " << eta << endl;
  }
  
  
  
  if(status==1)
  {
    cout << "Error while reading the parameters, program will be aborted" << endl;
    exit(666);
  }
    
}
  
 
vector<double> process_data(string str, int F)
{
  vector<double> result;
  stringstream  data(str);
  string line;
  while(getline(data,line,','))
  {
    result.push_back(atof(line.c_str()));
  }
  //Check for possible incompatibilities 
  if(result.size()!=F)
  {
    //Number of entries does not match number of flavors
    if(result.size()==1)
    {
      //Only a single entry given, will assume the same for all flavors
      result.resize(F,result[0]);
    }
    else
    {
      cout << "Error number of masses/chemical potentials does not match number of flavors" << endl;
      exit(1);
    }
  }  
  return result;
}
