/**
   \file CQED_StringvsMeson.cpp
   Generate a string on top of the groundstate and the related Meson states and compare their energy. Rather quick and dirty!
   
  \param <N> (int) Number of Fermions
  \param <d> (int) Physical dimension of Bosons in between
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <xs> (double) Coupling constant to start with
  \param <xe> (double Coupling constant to end with
  \param <e> (double) Constant of the Operator \f$ \Theta(n)\f$
  \param <lambda> (double) Constant for the penalty-term in the Hamiltonian
  \param <dt> (int) Size of a single time step
  \param <Nt> (int) Number of timesteps to evolve 
  \param <n> (int) Number of the site where the string starts 
  \param <l> (int) Length of the string which will be generated
  \param <ig> (string) filename for the inital guess
  \param <ac> (double) accuracy for the contractor
  \param <O> (int) Picks the desired algorithm for time evolution
  \author Stefan Kühn
  \date 22/08/2013

*/


#include <cstdlib>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <string>
//Lib which is needed for sleep command on *nix Systems
//#include "unistd.h"

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include "Indices.h"

#include "CQED.h"
#include "ramp_xvalue.h"

#define NUM_OF_PARAMS 9
#define FILENAME "cQED.txt"
#define RESOLUTION 10
#define EPS 1E-10
#define OFFSET 5192

using namespace std;
using namespace shrt;

bool file_exists(string filename);
void print_MPS(MPS);
void makeStringandMeson(int,int,MPS&,MPS&,MPS&,int,int);


int main(int argc,const char* argv[]){
  int cntr=0;
  
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << " N d D mu xs xe e la dt Nt ig ac O" << endl;
    cout << "N: " << "Number of spins" << endl
	 << "N0:" << "Total number of bosons" << endl
	 << "D: " << "Bond dimension of the MPS" << endl
	 << "µ: " << "Dimensionless constant" << endl
	 << "x: " << "Coupling constant" << endl
	 << "e: " << "Constant for gauge field" << endl
	 << "la:" << "Constant for penalty term" << endl
	 << "n: " << "Site where the string starts (n counts all sites, not only spin sites)" << endl
	 << "l: " << "Length of the string (must be an odd number)" << endl
	 << "ig:" << "True if inital guess for the ground state should be read form HDD (optional)" << endl
	 << "ac:" << "Accuracy for contractor (optional)" << endl
	 << endl;
    return -1;
  }
  
  //Number of spin sites
  int N=atoi(argv[++cntr]);
  //Total number fo bosons living on a link
  double d=atoi(argv[++cntr]);
  //double d2=atof(argv[++cntr]);
  //Bond dimension
  double D=atoi(argv[++cntr]);
  //Mass
  double mu=atof(argv[++cntr]);
  //Coupling constant
  double x=atof(argv[++cntr]);
  //Constant which can be used to manipulate the Hamiltonian
  double e=atof(argv[++cntr]);
  //Strength of the penalty
  double lam=atof(argv[++cntr]);
  //Position where the string starts
  int string_pos=atoi(argv[++cntr]);
  //Length of the string
  int string_length=atoi(argv[++cntr]);
  //Desired order of the evolution (by default 1)
  int evolution_order=1;
  //Inital state which should be used?
  bool ig = false;
  //Set the accuracy of contractor?
  bool set_accuracy = false;
  //Variable for the desired accuracy in case it should be set
  double desired_accuracy;
  //Name for file for initial state
  string ig_file;
  //Variable needed to verify if file for inital guess should have dummy name or a custom name
  size_t found_dat;
  
  //Only file for inital guess is given aditionally to all other parameters
  if(argc == NUM_OF_PARAMS + 2)
  {
    ig = true;
    ig_file = argv[++cntr];
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file: " << ig_file << endl;
  }
  //Custom accuracy for contractor and file for inital state
  if(argc == NUM_OF_PARAMS + 3)
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file " << ig_file << endl;
  }
  //Custom accuracy for contractor, file for inital state and order fo time evolution is given
  if(argc == NUM_OF_PARAMS + 4)
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    evolution_order = atoi(argv[++cntr]);
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file " << ig_file << endl;
  }
    
  //Print parameter which were set
  if(argc==NUM_OF_PARAMS+1){
    cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "N  = " << N << endl
	<< "d  = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "x  = " << x << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "n  = " << string_pos << endl
	<< "l  = " << string_length << endl
	<< "-----------" << endl << endl;
  }
  else if (argc==NUM_OF_PARAMS+2){
    cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "N  = " << N << endl
	<< "d  = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "x = " << x << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "n  = " << string_pos << endl
	<< "l  = " << string_length << endl
	<< "ig = " << ig << endl
	<< "-----------" << endl << endl;
  }
  else{
     cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "N  = " << N << endl
	<< "d  = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "x  = " << x << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "n  = " << string_pos << endl
	<< "l  = " << string_length << endl
	<< "ig = " << ig << endl
	<< "acc= " << desired_accuracy << endl
	<< "-----------" << endl << endl;
  }
    
       
  //Measure time
  time_t start,end;
  time_t start_local,end_local;
  double readtime_gs=0.0,writetime_gs=0.0,writetime_gauss=0.0,writetime_spin=0.0,writetime_ln=0.0,writetime_total=0.0,writetime_data=0.0;
  
  time(&start);
  
  //Total number fo sites (fermionic and bosonic)
  int L=2*N-1; 
  
  //Create an array containing the physical dimensions (alternating between fermion and boson, ending with fermion)
  int *phys_dims = new int[L];
  int *bond_dim = new int[L-1];
  for(int i=0; i<L; ++i){
    phys_dims[i] = 2;
    //cout << "phys_dim[" << i << "]=" << phys_dims[i] << endl;
    if((i+1)<L){
      phys_dims[i+1] = d+1;
      //cout << "phys_dim[" << i+1 << "]=" << phys_dims[i+1] << endl;
      i++;
    }
  }
  for(int i=0; i<(L-1); i++)
  {
    bond_dim[i] = D;
  }
  
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  //Set accuracy of contractor if desired
  if(set_accuracy){
    cout << "Setting Convergence Tolareance to " << desired_accuracy<<endl<<endl;
    contractor.setConvTol(desired_accuracy);
  }
  
  
  
  //Pointer to MPS to hold the inital state which will be evolved in time
  MPS *mps;
  
  //Get the Hamiltonian object for the cQED model
  CQED Hsngi(N,d,D,mu,x,e,lam);  
  //Get Hamiltonian MPO
  const MPO& hamil = Hsngi.getHMPO();  
  //MPO for the odd and even part of the evolution operator
  MPO exp_even(1),exp_odd(1);
  //Energy value
  double E0 = 0.0; 
  //Auxiliary MPS to store an old version before a timestep to be able to monitor changes
  MPS aux,mps_old;
  //Variable for the current energy value during time evolution
  complex_t curr_energy;
  //Variable for the current norm of the state during time evolution
  complex_t curr_norm;
  //Variable for the current overlap with the wavefunction of the previous step
  complex_t curr_overlap;
  //Variable for the value of the penalty energy
  complex_t penalty_energy;
  //MPO for the penalty term
  MPO Penalty(1);
  //Variable for the condensate
  complex_t cfraction;
  //MPO to determine the condensate
  MPO Condensate(1);
  //Variable for the current x value
  double xvalue=0.0;
  //File for energy values and other observables during the evolution
  ofstream energy_file;
  //Filepointer for the inital guess
  ifstream inital_guess;
  
  //Try to open file for inital state
  inital_guess.open(ig_file.c_str());
  //In case file for inital state exists read it and adjust bond dimension, otherwise start with predefined inital state
  if(inital_guess.is_open() && ig)
  {
    //Case in which an inital guess exists and will be read
    cout << "Found initial guess for the groundstate" << endl;
    time(&start_local);
    mps = new MPS();
    mps->importMPS(ig_file.c_str());
    cout << "Bond dimension after import: "<< mps->getBond() << endl;
    cout << "Length after import: " << mps->getLength() << endl;
    //In case an intial guess with bond dimension smaller than set is imported it has to be expanded again
    if(mps->getBond() < D)
    {
      //Increase bond dimension again, in case the inital guess has a smaller bond dimension
      //mps.increaseBondDimension(D);
      //Increase bond dimension again, in case the inital guess has a smaller bond dimension
      //This time with a little noise, maybe better to avoid stability issues?
      mps->increaseBondDimensionWithNoise(D,1E-2);
      cout << "Bond dimension after increase: "<< mps->getBond() << endl;
      
    }
    else if(mps->getBond() > D)
    {
      cout << "Error, the bond dimension of the inital guess is greater than the maximum bond dimension D=" << D << " set by the user" << endl;
      cout << "Bond dimension cannot be shrunk" << endl;
      exit(1);
    }
    time(&end_local);
    readtime_gs = difftime(end_local,start_local);
    cout << "Time for reading GS: " << readtime_gs << " s" << endl;
  }
  else
  {
    //Prepare some inital state
    mps = Hsngi.constructInitalMPS();
    //mps->setRandomState();
    //mps = Hsngi.constructInitalMPS(is_downup);
    //mps = Hsngi.constructInitalMPS(is_upup);
    //mps = Hsngi.constructInitalMPS(is_downdown);
    //mps = Hsngi.constructInitalMPS(is_string,3);        
    //mps = Hsngi.constructInitalMPS(is_pzoller_string);
  } 
  
  //Determine the inital energy
  cout << "Energy before time evolution " << contractor.contract(*mps,hamil,*mps) << endl;  
  //Determine condensate before evolution
  Hsngi.constructCondensateMPO(Condensate);
  cfraction = contractor.contract(*mps,Condensate,*mps);  
  cout <<"Intial condensate: " << cfraction << " with x-value " << Hsngi.get_x() << endl;
  
  //Still second order time evolution (but this first evolve to the ground state and afterwards generate a string and evolve this guy, compared to the other variant I calculate the groundstate variationally and therefore save some time, which might be useful for exploring parameters)
 
    
   MPO Ltot(1);
   MPS meson1,meson2;
   
   Hsngi.constructLtotMPO(Ltot);
   cout << "Total Flux inital state: " << contractor.contract(*mps,Ltot,*mps) << endl;
   Hsngi.constructCondensateMPO(Condensate);
   cfraction = contractor.contract(*mps,Condensate,*mps);
   cout << "Condensate inital state " << cfraction << endl;
    
    
  contractor.findGroundState(hamil,D,&E0,*mps);      
  cout << "Ground state energy: " << contractor.contract(*mps,hamil,*mps) << endl;
  cout << "Total Flux ground state: " << contractor.contract(*mps,Ltot,*mps) << endl;
  cfraction = contractor.contract(*mps,Condensate,*mps);
  cout << "Condensate ground state " << cfraction << endl;
    
  //Now generate a string
  cout << "Generating string and some meson states" << endl;
  //Hsngi.makeString(*mps,string_pos,string_length);
  makeStringandMeson(N,d,*mps,meson1,meson2,string_pos,string_length);
  
  
  ofstream mesonfile;
  
  if(file_exists("Meson.txt"))
    mesonfile.open("Meson.txt",ios_base::app);
  else
  {
    mesonfile.open("Meson.txt",ios_base::app);
    mesonfile << "#x\tmu\tl\tEgs\tEs\tEm1\tEm2\n";
  }
  
  
  cout << endl;
  cout << "Energy ground state: " << E0 << endl;
  cout << "Energy string state: " << contractor.contract(*mps,hamil,*mps) << endl;
  cout << "Energy meson1 state: " << contractor.contract(meson1,hamil,meson1) << endl;
  cout << "Energy meson2 state: " << contractor.contract(meson2,hamil,meson2) << endl;
  cout << endl;
  
  mesonfile << x << "\t" << mu << "\t" << string_length << "\t" << E0 << "\t" << contractor.contract(*mps,hamil,*mps).re << "\t" << contractor.contract(meson1,hamil,meson1).re << "\t" << contractor.contract(meson2,hamil,meson2).re << endl;
  
  mesonfile.close();
  
  
  
  cout << "Norm after string generation is " << contractor.contract(*mps,*mps) << endl;
  cout << "Total Flux after string generation: " << contractor.contract(*mps,Ltot,*mps) << endl;
  cfraction = contractor.contract(*mps,Condensate,*mps);
  cout << "Condensate after string generation " << cfraction << endl;
    
  double lambdak=0.0;
  MPS first_ex(L,bond_dim,phys_dims);   
    
  
  
  E0 = contractor.contract(*mps,hamil,*mps).re;  
  
  //Try to calculate the charge
  const MPO& chargeop=Hsngi.getChargeMPO();  
  complex_t charge_0=contractor.contract(*mps,chargeop,*mps);
  cout << "Charge Groundstate " << charge_0.re << endl;  
  
  //Calculate the terms for the Gauss Law and see if it is fullfilled
  ofstream Lnfile;
  Lnfile.open("Ln.txt", ios_base::app);
  complex_t Lnvalue;
  MPO Ln(1);
  cout << endl;
  for(int i=1; i<L-1;i += 2){
    Hsngi.constructLopMPO(Ln,i);
    Lnvalue = contractor.contract(*mps,Ln,*mps);
    cout << "--> Expectation value L("<< (i+1)/2 << ") = " << Lnvalue<< endl;
    time(&start_local);
    Lnfile << Lnvalue.re << "\t";
    time(&end_local);
    writetime_ln += difftime(end_local,start_local);
  }
  time(&start_local);
  Lnfile << endl;
  Lnfile.close();
  time(&end_local);
  writetime_ln += difftime(end_local,start_local);
  
  cout << "Time for writing Ln: " << writetime_ln << " s" << endl;
  
  ofstream Gaussfile;
  Gaussfile.open("Gaussfile.txt", ios_base::app);
  complex_t Gaussvalue;
  MPO Gauss(1);
  cout << endl;
  for(int i=0; i<L;i += 2){
    Hsngi.constructGaussMPO(Gauss,i);
    Gaussvalue = contractor.contract(*mps,Gauss,*mps);
    cout << "--> Expectation value 0.5*[sigma_z+(-1)^"<< i/2+1 << "] = " << Gaussvalue<< endl;
    time(&start_local);
    Gaussfile << Gaussvalue.re << "\t";
    time(&end_local);
    writetime_gauss += difftime(end_local,start_local);
  }
  time(&start_local);
  Gaussfile << endl;
  Gaussfile.close();
  time(&end_local);
  writetime_gauss += difftime(end_local,start_local);
  
  cout << "Time for writing Gauss: " << writetime_gauss << " s" << endl;
  
  //Calculate the z-Component of the spins
  ofstream spinfile;
  spinfile.open("Spin.txt", ios_base::app);
  complex_t Zcomponent;
  MPO Sigmaz(1);
  cout << endl;
  for(int i=0; i<L;i += 2){
    Hsngi.constructsigmazMPO(Sigmaz,i);
    Zcomponent = contractor.contract(*mps,Sigmaz,*mps);
    cout << "--> Expectation value sigma_z("<< i/2+1 << ") = " << Zcomponent << endl;
    time(&start_local);
    spinfile << Zcomponent.re << "\t";
    time(&end_local);
    writetime_spin += difftime(end_local,start_local);
  }
  time(&start_local);
  spinfile << endl;
  spinfile.close(); 
  time(&end_local);
  writetime_spin += difftime(end_local,start_local);
   
  cout << "Time for writing Spin: " << writetime_spin << " s" << endl;
  
  cout << endl << "Meson 1" << endl;
  
  for(int i=1; i<L-1;i += 2){
    Hsngi.constructLopMPO(Ln,i);
    Lnvalue = contractor.contract(meson1,Ln,meson1);
    cout << "--> Expectation value L("<< (i+1)/2 << ") = " << Lnvalue<< endl;
    time(&start_local);
    Lnfile << Lnvalue.re << "\t";
    time(&end_local);
    writetime_ln += difftime(end_local,start_local);
  }
  cout << endl;
   for(int i=0; i<L;i += 2){
    Hsngi.constructsigmazMPO(Sigmaz,i);
    Zcomponent = contractor.contract(meson1,Sigmaz,meson1);
    cout << "--> Expectation value sigma_z("<< i/2+1 << ") = " << Zcomponent << endl;
    time(&start_local);
    spinfile << Zcomponent.re << "\t";
    time(&end_local);
    writetime_spin += difftime(end_local,start_local);
  }
  
  cout << endl << "Meson 2" << endl;
  for(int i=1; i<L-1;i += 2){
    Hsngi.constructLopMPO(Ln,i);
    Lnvalue = contractor.contract(meson2,Ln,meson2);
    cout << "--> Expectation value L("<< (i+1)/2 << ") = " << Lnvalue<< endl;
    time(&start_local);
    Lnfile << Lnvalue.re << "\t";
    time(&end_local);
    writetime_ln += difftime(end_local,start_local);
  }
  cout << endl;
   for(int i=0; i<L;i += 2){
    Hsngi.constructsigmazMPO(Sigmaz,i);
    Zcomponent = contractor.contract(meson2,Sigmaz,meson2);
    cout << "--> Expectation value sigma_z("<< i/2+1 << ") = " << Zcomponent << endl;
    time(&start_local);
    spinfile << Zcomponent.re << "\t";
    time(&end_local);
    writetime_spin += difftime(end_local,start_local);
  }
  
  
  //Calculate the condensate farction
  cout << endl;
  Hsngi.constructCondensateMPO(Condensate);
  cfraction = contractor.contract(*mps,Condensate,*mps);
  cout << "Condensate: "<< cfraction << endl << endl;  
  
  //Ouput the results up to now (sometimes the computation of the exicted state(s) fails and the Contractor routine quits the program, therefore all results which were not stored would be lost). To avoid this, save the first part of the results now
  ofstream resfile;  
  cout << "Saving results..." << endl;
  if(file_exists(FILENAME))
  {
    time(&start_local);
    resfile.open(FILENAME,ios_base::app);  
    resfile << endl << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << x << "\t" << e << "\t" << lam << scientific <<"\t" << E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" <<penalty_energy.re << "\t";
  }
  else
  {
    time(&start_local);
    resfile.open(FILENAME,ios_base::app);  
    resfile << "#N" << "\t" << "d" << "\t" << "D" << "\t" << "mu" << "\t" << "x" << "\t" << "e" << "\t" << "lambda" << "\t"<< "E0" << "\t" << "C0" << "\t" << "Cfrac" << "\t" << "E1" << "\t" << "C1" << "\tE_pen"<< endl;
    resfile << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << x << "\t" << e << "\t" << lam << scientific << "\t"<< E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" << penalty_energy.re << "\t";
  }
  resfile.close();
  time(&end_local);
  writetime_data = difftime(end_local,start_local);
  
  cout << "Time for writing Data: " << writetime_data << " s" << endl;  
  
  //Get the excited state(s) 
  /*
  vector<MPS*> excitedstates;
  double lambdak;
  MPS first_ex(L,bond_dim,phys_dims);
  
  excitedstates.push_back(&mps);
  
  cout << endl << "Computing excited states" << endl;
  contractor.findNextExcitedState(hamil,D,excitedstates,&lambdak,first_ex); 
  cout << "Excited state energy: " << lambdak << endl<< endl;
  complex_t charge_1=contractor.contract(first_ex,chargeop,first_ex);
  
  //Save results for the exicted state(s)
  cout << "Saving results..." << endl;
  resfile.open(FILENAME,ios_base::app);  
  resfile << lambdak << "\t" << charge_1.re ;
  resfile.close();*/  
  
  //Free Memory
  delete[] phys_dims;
  delete[] bond_dim; 
  mps->clear();
  delete mps;
  
 
  
  //End Time
  time(&end);
  cout << "Run time: " << difftime(end,start) <<"s" <<  endl;  
  writetime_total = writetime_data + writetime_gauss + writetime_gs + writetime_ln + writetime_spin;
  cout << "Time for writing in total: " << writetime_total << " s" << endl;
  cout << "Time on filesystem in total (reading + writing): " << writetime_total+readtime_gs << endl; 
  cout << "End of program" << endl;
  return 0;
}

bool file_exists(string filename)
{
  ifstream file;
  file.open(filename.c_str());
  if(file.good())
    return true;
  else
    return false;
}

void print_MPS(MPS mps)
{
  cout << "MPS of length " << mps. getLength() << endl;
  for(int i=0; i<mps.getLength(); i++)
  {
    cout << "[" << i << "]=" <<  (mps.getA(i)).getA() << endl;
  }
}

void makeStringandMeson(int N_fermions,int d_bose,MPS &stringmps, MPS &meson_1_mps, MPS &meson_2_mps, int pos_start, int length)
{
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  
  complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
  mwArray sigma_minus(Indices(2,2),data_sm);
  
  complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
  mwArray sigma_plus(Indices(2,2),data_sp);
 
  //Provide bosonic operators
  mwArray Op1(Indices(d_bose+1,d_bose+1));
  mwArray id_bosonic = identityMatrix(d_bose+1);
  mwArray Op2;
  
  Op1.fillWithZero();
  
  //Construct Op1 
  //Since I want to reuse my old implementation which has hopping terms like sp Op1 sm + sm op2 sp I absorb i in Op1 and -i in Op2
  double norm_factor = sqrt(d_bose/2.0*(d_bose/2.0+1));
  complex_t tmp;
  for(int i=0; i<=d_bose; ++i)
  {
    for(int j=0; j<=d_bose; ++j)
    {
      //Lp
      if(i==j+1)
      {
	Op1.setElement(0.0,sqrt(i*(d_bose-j))/norm_factor,Indices(i,j));
      }
    }
  }
  
  //Get Op2 as the hermitian conjugate of Op1
  Op2=Op1;
  Op2.transpose(true);
  
  //In case pos_start is a negative number just prepare a string approximately in the middle
  if(pos_start<0)
  {
    //Find the middle of the chain where the string should reside
    int middle = (int) (N_fermions/2.0);
    //Translate into absolute site number
    pos_start = 2*(middle-(length-1)/2)-2;    
    cout << "Construction string starting at site " << pos_start << endl;
  }
  
  //Check if string fits on chain in principle (not paying attention to Gauss Law and stuff which might be an additional limitation)
  if((stringmps.getLength()-pos_start)< 2*length)
  {
	  cerr << "MPS is not long enough to generate the desires string, program will be aborted" << endl;
	  exit(1);
  }
  
  //Check if startposition is a valid spin site and if it is a spin up or a spin down site
  if((pos_start%2)==1)
  {
	  cerr << "Cannot generate a string starting at a bosonic site, program will be aborted..." << endl;
	  exit(1);
  }
  int sign = pow(-1,pos_start/2+1);
  /*cout << "sign=" << sign << endl;
  cout << "length=" << length << endl;*/
  
  if(sign==1)
  {
	 stringmps.applyLocalOperator(pos_start,sigma_plus,true); 
	 for(int i=pos_start+1;i<pos_start+2*length;i+=2)
		stringmps.applyLocalOperator(i,Op1,true);
	 stringmps.applyLocalOperator(pos_start+2*length,sigma_minus,true);
  }
  else
  {
	 stringmps.applyLocalOperator(pos_start,sigma_minus,true); 
	 for(int i=pos_start+1;i<pos_start+2*length;i+=2)
		stringmps.applyLocalOperator(i,Op2,true);
	 stringmps.applyLocalOperator(pos_start+2*length,sigma_plus,true);	
  }
  
  meson_1_mps = stringmps;
  
  for(int i=pos_start+2;i<pos_start+2*length-2; i+=4)
  {
    if(sign==-1)
    {
      meson_1_mps.applyLocalOperator(i,sigma_plus,true); 
      meson_1_mps.applyLocalOperator(i+1,Op1,true);      
      meson_1_mps.applyLocalOperator(i+2,sigma_minus,true); 
    }
    else
    {
      meson_1_mps.applyLocalOperator(i,sigma_minus,true); 
      meson_1_mps.applyLocalOperator(i+1,Op2,true);      
      meson_1_mps.applyLocalOperator(i+2,sigma_plus,true); 
    }
  }
  meson_2_mps = meson_1_mps;
  for(int i=pos_start+4;i<pos_start+2*length-2; i+=4)
  {
    
    if(sign==-1)
    {
      meson_2_mps.applyLocalOperator(i,sigma_plus,true); 
      meson_2_mps.applyLocalOperator(i+1,Op1,true);      
      meson_2_mps.applyLocalOperator(i+2,sigma_minus,true); 
    }
    else
    {
      meson_2_mps.applyLocalOperator(i,sigma_minus,true); 
      meson_2_mps.applyLocalOperator(i+1,Op2,true);      
      meson_2_mps.applyLocalOperator(i+2,sigma_plus,true); 
    }    
  }
  
  
}