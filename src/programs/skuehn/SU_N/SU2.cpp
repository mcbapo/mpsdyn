/**
   \file SU2.cpp
   Realisation of the truncated SU(2) Hamiltonian with staggered fermions from Erez, Ignacio and Benni
   
  \param <N> (int) Number of Fermions
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <epsilon> (double) Hopping parameter
  \param <mu> (double) Mass
  \param <g> (double) Coupling constant
  
  \author Stefan Kühn
  \date 12/08/2014

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
#define FILENAME "SU2.txt"
#define FILENAME_EX "SU2_ex.txt"
#define FILENAME_CONFIG "SU2_configuration.txt"
#define CHARGE_FILE "SU2charge.txt"

using namespace std;
using namespace shrt;

//bool file_exists(string filename);
void print_MPS(MPS);
//Get full config of a state to check the Gauss Law manually afterwars
void get_full_config(SU2Hamiltonian &H, Contractor &contractor, MPS &mps, int N);

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
	 << "n:       " << "Number of excited states (if a positive number is given)" << endl
	 << "ig_file: " << "Filname containing an initial guess (optional)" << endl
	 << "factor:  " << "Factor which is used to enlarge two masses to create a pair (optional)" << endl
	 << "leng:    " << "Seperation of the pairs in sites, should be odd (optional)" << endl
	 << "lambda:  " << "Penalty strength (optional)" << endl;
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
  int num_of_exstates=atoi(argv[++cntr]);
  //Possible input file
  string ig_file="";
  //Factor which is used to enlarge the masses for the pair
  double mass_factor;
  //Seperation of the charge-anticharge pair
  int separation;
  //Masses
  vector<double> masses;
  masses.push_back(mu);
  vector<int> chargepos;
  //Penalty for charge configuration
  double lambda=0.0;
  
  if(argc >= NUM_OF_PARAMS+2)
    ig_file=argv[++cntr];
  
  if(argc >= NUM_OF_PARAMS+4)
  {
    mass_factor=atof(argv[++cntr]);
    separation=atoi(argv[++cntr]);
    if((separation>N || (separation%2)!=1) && separation>0)
    {
      cout << "Seperation must be odd and smaller than given number of sites!" << endl;
      exit(666);
    }
    else if(separation>0)
    {
      masses.clear();
      masses.assign(N,mu);
      //Compute where the charge anticharge pair is sitting
      int start_index=round(N/2)-round(separation/2), end_index=start_index+separation;
      chargepos.push_back(start_index);
      chargepos.push_back(end_index);      
      masses[start_index-1]=1.0*mass_factor*mu;
      masses[end_index-1]=1.0*mass_factor*mu;
    }
    else
    {
      cout << "Taking homgeneous mass" << endl;
      masses.clear();
      masses.assign(N,mu);      
    }
  }  
  cout << "Masses: " << masses << endl;
  
  if(argc == NUM_OF_PARAMS+5)
    lambda = atof(argv[++cntr]);
  
  //File for results
  ofstream file;
  //Variables for Gauss Law Components
  complex_t Gx_val,Gy_val,Gz_val; 
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  contractor.setConvTol(1.0E-8);

  //Get the Hamiltonian
  SU2Hamiltonian HNGI(N,epsilon,mu,g);
  //SU2Hamiltonian HNGI(N,epsilon,masses,g);
  //chargepos.clear();
  //SU2Hamiltonian HNGI(N,epsilon,masses,g,chargepos,lambda);
  const MPO& hamil=HNGI.getHMPO();    
  
  cout << "Successfully retrived SU2Hamiltonian" << endl;	cout.flush();
  
  //Now compute the groundstate
  vector<int> dims(3*N-1,2);
  for(int i=2; i<(3*N-1); i+=3) 
    dims[i]=5;
  MPS mps(3*N-1,D,dims);
  double E0=0.0,E1=0.0;
  
  MPO Gx(1), Gy(1), Gz(1);
  HNGI.getGaussLawMPO(Gx,x_comp);
  HNGI.getGaussLawMPO(Gy,y_comp);
  HNGI.getGaussLawMPO(Gz,z_comp);
  
  MPO Condensate(1);
  HNGI.getCondensateMPO(Condensate,true);
  /*for(int i=1; i<=10; i++)
  {
    mps.setRandomState();
    cout << i << ": <Psi|H|Psi>=" << contractor.contract(mps,hamil,mps)<<"	<Psi|Psi>=" << contractor.contract(mps,mps)<< "	G_x=" << contractor.contract(mps,Gx,mps) <<"	G_y=" << contractor.contract(mps,Gy,mps) <<"	G_z=" << contractor.contract(mps,Gz,mps) <<endl;
  }*/
  
  //Choose inital state, random or strong coupling ground state, or given input file
  if(!ig_file.empty() && file_exists(ig_file))
  {
    cout << "Trying to import ground state from " << ig_file << endl;
    mps.importMPS(ig_file.c_str());
    cout << "Successfully read GS from file " << ig_file << endl;
  }
  else
  {
    //Create a string on the interacting vacuum
    /*HNGI.constructInitialMPS(mps);
    MPO StringOp(1);
    MPS aux = MPS(mps);
    HNGI.getStringMPO2(StringOp,separation,0,false);
    aux.approximate(StringOp,mps,D); 
    contractor.optimize(StringOp,mps,aux,D);
    mps = MPS(aux);
    cout << "MPS: " << mps << endl;*/
    
    //Random start state
    //mps.setRandomState();
    
    //Strong coupling ground state
    HNGI.constructInitialMPS(mps);
    
    //String betweeen heavy external charges
    //HNGI.getHeavyString(mps,separation);
    
    //Get the configuration once before the ground state search
    //get_full_config(HNGI,contractor,mps,N);
  }
  
  
  
  cout << "------------------------------------------------------------------------"<<endl;
  cout << "---------------Starting to compute ground state----------------------"<<endl;
  cout << "------------------------------------------------------------------------"<<endl; 
  cout.flush();
  contractor.findGroundState(hamil,D,&E0,mps);
  Gx_val=contractor.contract(mps,Gx,mps);
  Gy_val=contractor.contract(mps,Gy,mps);
  Gz_val=contractor.contract(mps,Gz,mps);
  
  cout << "Gx: " << Gx_val << endl;
  cout << "Gy: " << Gy_val << endl;
  cout << "Gz: " << Gz_val << endl;
  
  cout << "Chiral condensate: " << contractor.contract(mps,Condensate,mps) << endl;
  
  //Get the configuration once after the ground state search to make sure if I'm staying in the right Gauss Law sector (or if the search spoiled it)
  get_full_config(HNGI,contractor,mps,N);
  
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
    file << setprecision(16) << N << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << E0 << "\t" << Gx_val << "\t" << Gy_val << "\t" << Gz_val << "\t" << contractor.contract(mps,Condensate,mps) << endl;
  }
  else
  {
    file.open(FILENAME, ios_base::app);
    file << setprecision(16) << "#N\tD\tepsilon\tmu\tg\tE0\tGx\tGy\tGz\tC" << endl;
    file << setprecision(16) << N << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << E0 << "\t" << Gx_val << "\t" << Gy_val << "\t" << Gz_val << "\t" << contractor.contract(mps,Condensate,mps) << endl;
  }
  file.close();
  
  //Get the configuration and save it
  vector<complex_t> configuration;  
  MPO mpo(1),charge(1);
  file.open(FILENAME_CONFIG, ios_base::app);
  for(int k=1; k<=N; k++)
  {
    HNGI.getSpin1MPO(mpo,k);
    configuration.push_back(contractor.contract(mps,mpo,mps));
    
    HNGI.getSpin2MPO(mpo,k);
    configuration.push_back(contractor.contract(mps,mpo,mps));
    
    if(k<N)
    {
      /*HNGI.getLMPO(mpo,k);
      configuration.push_back(contractor.contract(mps,mpo,mps));
      
      HNGI.getRMPO(mpo,k);
      configuration.push_back(contractor.contract(mps,mpo,mps));*/
      
      HNGI.getJsquaredMPO(mpo,k);
      configuration.push_back(contractor.contract(mps,mpo,mps));
      
      //Write to file
      file <<configuration[3*k-3] << "\t" << configuration[3*k-2] << "\t" << configuration[3*k-1] << "\t";
    }
    else
    {
      file <<configuration[3*k-3] << "\t" << configuration[3*k-2] << endl;
    }
  }
  file.close();
  //cout << "Configuration: " << configuration << endl;
  
  //Determine the chargequare of the ground state
  file.open(CHARGE_FILE, ios_base::app);
  for(int pos=1; pos<=N; pos++)
  {	
    HNGI.getChargeSquareMPO(charge,pos);
    file << contractor.contract(mps,charge,mps) << "\t";
  }
  //cout << endl;
  file << endl;
  file.close();
  
  
  //Try to get an excited state
  MPS first_ex;
  first_ex = mps;
  
  vector<MPS*> computed_states;
  computed_states.push_back(&mps);
  
  cout << endl;
  vector<MPS> excited_states;
  vector<double> energies;
  vector<complex_t> Gx_vals_ex,Gy_vals_ex,Gz_vals_ex;
  
  //Save also the ground state
  energies.push_back(E0);
  Gx_vals_ex.push_back(Gx_val);	Gy_vals_ex.push_back(Gy_val);	Gz_vals_ex.push_back(Gz_val);
  
  for(int i=1; i<=num_of_exstates; i++)
  {
    first_ex.setRandomState();  
    contractor.findNextExcitedState(hamil,D,computed_states,&E1,first_ex);
    computed_states.push_back(new MPS(first_ex));
    cout << i <<". excited state with energy " << E1 << endl;
    Gx_val=contractor.contract(first_ex,Gx,first_ex);
    Gy_val=contractor.contract(first_ex,Gy,first_ex);
    Gz_val=contractor.contract(first_ex,Gz,first_ex);
    
    cout << "Gx: " << Gx_val << endl;
    cout << "Gy: " << Gy_val << endl;
    cout << "Gz: " << Gz_val << endl<< endl;
    
    energies.push_back(E1);
    Gx_vals_ex.push_back(Gx_val);	Gy_vals_ex.push_back(Gy_val);	Gz_vals_ex.push_back(Gz_val);
  }
  
  //Free the memory, do NOT delete the ground state!
  for(int i=1; i<computed_states.size(); i++)
  {
    computed_states[i]->clear();
    delete computed_states[i];
  }
  
  //Save observables for excited states
  file.open(FILENAME_EX, ios_base::app);
  file << setprecision(16) << "#N=" << N << "\t" << "D=" << D << "\t" << "epsilon=" << epsilon << "\t" << "mu=" << mu << "\t" << "g=" << g << "\t" <<"E0=" <<E0 << endl;
  for(int i=0; i<energies.size(); i++)
  {
    if(i==0)
      file << "Gs:\t" << energies[i] << "\t" << Gx_vals_ex[i] << "\t" << Gy_vals_ex[i] << "\t" << Gz_vals_ex[i] << endl;
    else
      file << i << ". Ex:\t" << energies[i] << "\t" << Gx_vals_ex[i] << "\t" << Gy_vals_ex[i] << "\t" << Gz_vals_ex[i] << endl;
  }
  file << endl;
  file.close();
  

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


    
    
