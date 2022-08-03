/**
   \file schwingerMultiFlavor.cpp
   Multi Flavor Schwinger model
   
  \param <file> (string) Name of the driver file containing the parameters
  
  \author Stefan Kühn
  \date 05/04/2015

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

#include "SchwingerHamiltonian.h"
#include "MultiFlavorSchwingerHamiltonian.h"

#define NUM_OF_PARAMS 1
#define FILENAME "MultiFlavorRefined.txt"
#define TMP_FILE "GS_refined.tmp"


using namespace std;
using namespace shrt;

vector<complex_t> get_configuration(MultiFlavorSchwingerHamiltonian &H, Contractor &contractor, MPS &mps, int N, int F);
vector<complex_t> get_charge(vector<complex_t> spinconfiguration, int N, int F);
vector<complex_t> get_condensate(vector<complex_t> spinconfiguration, int N, int F);
void save_configuration(vector<complex_t> configuration, string filename);
void read_parameters(string filename,int &N, int &F, int &D, double &x, vector<double> &mu, vector<double> &nu, double &lambda, double &l0, int &bu_interval, string &init_file, double &tol);
vector<double> process_data(string str, int F);

int main(int argc,const char* argv[])
{  
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage:   " << argv[0] << " <driver_file> <tolerance>" << endl;
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
  double E0;
  int N;
  int D;
  int F;
  double x,lambda,l0;
  ofstream file;
  string init_file;
  int bu_interval=50;
  double tol;
  double Ni_tol=1E-1;
  
  if(argc >(NUM_OF_PARAMS+1))
  {
    Ni_tol = atof(argv[2]);
    cout << "Using N_i tolerance " << Ni_tol << endl;
  }
  else
    cout << "Using default N_i tolerance " << Ni_tol << endl;
  
  //Read input parameters from the input file
  read_parameters(inputfile,N,F,D,x,mu,nu,lambda,l0,bu_interval,init_file,tol);

  contractor.setConvTol(tol);
  
  MPS mps;
  
  //Get the Hamiltonian
  MultiFlavorSchwingerHamiltonian Hmulti(N,F,x,mu,nu,lambda,l0);
  const MPO& hamil=Hmulti.getHMPO();    
  
  cout << "Successfully retrieved multi flavor Schwinger Hamiltonian" << endl;	cout.flush();   

  
  //Set up all kins of MPOs
  MPO charge(1),penalty(1);
  vector<MPO*> condensates,particles;
  MPO *condensate_flavor,*particles_flavor;
  Hmulti.getChargeMPO(charge);
  Hmulti.getChargePenaltyMPO(penalty,true);
  for(int i=0; i<F ; i++)
  {
    condensate_flavor=new MPO(1);
    particles_flavor=new MPO(1);
    cout << "Getting condensate for flavor " << i << endl;
    Hmulti.getCondensateMPO(*condensate_flavor,i,true);
    Hmulti.getParticleNumberMPO(*particles_flavor,i);
    condensates.push_back(condensate_flavor);      
    particles.push_back(particles_flavor);
  }
  condensate_flavor=new MPO(1);
  particles_flavor=new MPO(1);
  Hmulti.getCondensateMPO(*condensate_flavor,-1,true);
  Hmulti.getParticleNumberMPO(*particles_flavor,-1);
  condensates.push_back(condensate_flavor);
  particles.push_back(particles_flavor);
  
  //Read the state which should be refined
  if(file_exists(init_file))
  {
    mps.importMPS(init_file.c_str());
    cout << "Successfully loaded initial MPS" << endl;
  }
  else
  {
    cout << "Error initial state " << init_file << " does not exist, will abort program" << endl;
    exit(666);
  }
  
  //Now figure out the particle number and see if it needs refinement within the tolerance I allow  
  //\TODO: Make it more general, right now it only works for 2 flavors
  vector<double> Nf,Nfdiff;
  bool refine_flag=false;
  for(int i=0; i<F; i++)
    Nf.push_back(real(contractor.contract(mps,*particles[i],mps)));  
  Nfdiff.resize(Nf.size());  
  std::adjacent_difference(Nf.begin(), Nf.end(), Nfdiff.begin());
  
  /*Now check if refinement is necessary, thus check if:
   * - N2-N1 is close to an integer
   * - in case it is, if this integer is even
   * - N1 is close to an integer
   * - N2 is close to an integer */
  if(fabs(round(Nfdiff[1])-Nfdiff[1])>Ni_tol || (((int) round(Nfdiff[1]))%2) != 0 || fabs(round(Nf[0])-Nf[0])>Ni_tol || fabs(round(Nf[0])-Nf[0])>Ni_tol)
    refine_flag=true;  
  
  if(!refine_flag)
  {
    cout << "N_i values are within the given tolerance " << Ni_tol << ", no refining needed" << endl;
    goto stop;
  }
  else
    cout << "Going for refinement:" << endl;
  
  //Get the projectors for the relevant tuples of particle numbers and project the state
  {
  //Compute the possible values of N_1-N_2 according to the theory
  for(int i=0; i<Nf.size(); i++)
    cout << "- Nf[" << i << "]=" << Nf[i] << endl;
  cout << "- dN=" << Nfdiff[1] << endl;
    
  int dNlower = floor(fabs(Nfdiff[1])/2.0)*2;
  int dNupper = ceil(fabs(Nfdiff[1])/2.0)*2;
  
  
  //Now compute the resulting N1,N2 values for the projectors
  int N1lower,N1upper,N2lower,N2upper;
  
  N1lower = (int) ((N+dNlower)/2.0);
  N2lower = (int) ((N-dNlower)/2.0);
  
  N1upper = (int) ((N+dNupper)/2.0);
  N2upper = (int) ((N-dNupper)/2.0);
  
  cout << "- dNlower: " << dNlower << ", N1lower: " << N1lower << ", N2lower: " << N2lower << endl;
  cout << "- dNupper: " << dNupper << ", N1upper: " << N1upper << ", N2upper: " << N2upper <<endl << endl;
  
  MPO P1lower(1),P1upper(1),P2lower(1),P2upper(1),Plower(1),Pupper(1); 
  double val_lower=0.0,val_upper=0.0;
  complex_t tmp,tmp1,tmp2;
  
  Hmulti.getParticleNumberProjector(P1lower,0,N1lower);
  Hmulti.getParticleNumberProjector(P2lower,1,N2lower);
  const MPO* oprsLower[2];
  const MPO* oprsUpper[2];
  
  /*Plower.initLength(2*N);
  oprsLower[0]=&P1lower;
  oprsLower[1]=&P2lower;  
  MPO::join(2,oprsLower,Plower);
  tmp  = contractor.contract(mps,Plower,mps);
  cout << "Expectation value of Plower: " << tmp << endl;
  val_lower = real(tmp);*/
  
  //Different strategy, it should be cheaper to have both seperately
  tmp1 = contractor.contract(mps,P1lower,mps);
  tmp2 = contractor.contract(mps,P2lower,mps);
  val_lower = real(tmp1)*real(tmp2);
  cout << "- Expectation value of P1lower: " << tmp1 << endl;
  cout << "- Expectation value of P2lower: " << tmp2 << endl;
  cout << "- Val lower: " << val_lower << endl << endl;
  
  if(N1lower != N1upper && N2lower != N2upper)
  {
    Hmulti.getParticleNumberProjector(P1upper,0,N1upper);
    Hmulti.getParticleNumberProjector(P2upper,1,N2upper);
    /*Pupper.initLength(2*N);
    oprsUpper[0]=&P1upper;
    oprsUpper[1]=&P2upper;
    MPO::join(2,oprsUpper,Pupper); 
    tmp  = contractor.contract(mps,Pupper,mps);
    cout << "Expectation value of Pupper: " << tmp << endl;
    val_upper = real(tmp);*/
    
    //Different strategy, it should be cheaper to have both seperately
    tmp1 = contractor.contract(mps,P1upper,mps);
    tmp2 = contractor.contract(mps,P2upper,mps);
    val_upper = real(tmp1)*real(tmp2);
    cout << "- Expectation value of P2upper: " << tmp1 << endl;
    cout << "- Expectation value of P2upper: " << tmp2 << endl;
    cout << "- Val upper: " << val_upper << endl << endl;
  }
  
  
  //Now do actual projection, therefore get a copy of the original mps
  MPS init_mps(mps);
  if(val_upper>val_lower)
  {
    //contractor.optimize(Pupper,init_mps,mps,D);
    //Different strategy, as I anyway want it to be an eigenstate of P1 and P2 project one after another
    contractor.optimize(P1upper,mps,init_mps,D);
    init_mps.setNormFact(1.0);
    mps = init_mps;
    contractor.optimize(P2upper,init_mps,mps,D);
  }
  else
  {
    //contractor.optimize(Plower,init_mps,mps,D);
    //Different strategy, as I anyway want it to be an eigenstate of P1 and P2 project one after another
    contractor.optimize(P1lower,mps,init_mps,D);
    init_mps.setNormFact(1.0);
    mps = init_mps;
    contractor.optimize(P2lower,init_mps,mps,D);
  }
  mps.gaugeCond('R',true);  
  }
  
  
  cout << "---------------------------------------------------------------------"<<endl;
  cout << "--------------------Starting to refine state-------------------------"<<endl;
  cout << "---------------------------------------------------------------------"<<endl; 
  cout.flush();
  contractor.findGroundState(hamil,D,&E0,mps,0.0,1,TMP_FILE,bu_interval);  
  mps.exportMPS("GS_refined.dat");
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
    for(int i=0; i<=F; i++)
      file << setprecision(15) << contractor.contract(mps,*condensates[i],mps) << "\t";
    for(int i=0; i<=F; i++)
      file << setprecision(15) << contractor.contract(mps,*particles[i],mps) << "\t";
    file << setprecision(15) << contractor.contract(mps,charge,mps) << "\t" << contractor.contract(mps,penalty,mps) << endl;
  }
  else
  {
    file.open(FILENAME,ios_base::app);
    file <<"#N\tD\tF\tx\tmu[0]..mu[F-1]\tnu[0]..nu[F-1]\tE0\tC[0]..C[F-1]\tC\tN[0]..N[F-1]\tNtot\tC_0\tP\n";
    file << N << "\t" << D << "\t" << F << "\t" << x << "\t";
    for(int i=0; i<mu.size(); i++)
      file << setprecision(15) << mu[i] << "\t";
    for(int i=0; i<nu.size(); i++)
      file << setprecision(15) << nu[i] << "\t";
    file << setprecision(15) << E0 << "\t";
    for(int i=0; i<=F; i++)
      file << setprecision(15) << contractor.contract(mps,*condensates[i],mps) << "\t";
    for(int i=0; i<=F; i++)
      file << setprecision(15) << contractor.contract(mps,*particles[i],mps) << "\t";
    file << setprecision(15) << contractor.contract(mps,charge,mps) << "\t" << contractor.contract(mps,penalty,mps) << endl;
  }
  file.close();
  
  
  //Compute the ground state configuration
  {
  vector <complex_t> spin_vals,condensate_vals,charge_vals;
  spin_vals=get_configuration(Hmulti,contractor,mps,N,F);
  save_configuration(spin_vals,"MultiFlavorConfigurationRefined.txt");
  charge_vals=get_charge(spin_vals,N,F);
  save_configuration(charge_vals,"MultiFlavorChargeRefined.txt");
  condensate_vals=get_condensate(spin_vals,N,F);  
  save_configuration(condensate_vals,"MultiFlavorCondensateRefined.txt");
  }
  
  //Free memory
  stop:
  for(int i=0; i<condensates.size(); i++)
  {
    condensates[i]->clear();    
    particles[i]->clear();
    delete condensates[i];
    delete particles[i];
  }   

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

void read_parameters(string filename,int &N, int &F, int &D, double &x, vector<double> &mu, vector<double> &nu,double &lambda, double &l0, int &bu_interval,string &init_file, double &tol)
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
