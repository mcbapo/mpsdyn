/**
   \file KS_Hamil_test.cpp
   At the moment file purely for testing the unified implementation for the Kogut Susskind type Hamiltonian (with arbitrary operators on the links)
   
  \param <N> (int) Number of Fermions
  \param <d> (int) Physical dimension of Bosons in between
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <x> (double) Coupling constant
  \param <e> (double) Constant of the Operator \f$ \Theta(n)\f$
  \param <lambda> (double) Constant for the penalty-term in the Hamiltonian
   
  \author Stefan Kühn
  \date 17/02/2014

*/


#include <cstdlib>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <string>
#include <algorithm>
//Lib which is needed for sleep command on *nix Systems
//#include "unistd.h"

//Use temporary storage for intermediate results
#define TMPSTORE 0

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include "misc.h"

#include "KogutSusskindHamiltonian.h"

#define NUM_OF_PARAMS 8
#define FILENAME "KS.txt"
#define EPSILON 1.0E-9

using namespace std;

//bool file_exists(string filename);

int main(int argc,const char* argv[]){
  int cntr=0;
  if(argc < NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << " Model N d D mu x e la" << endl;
    cout << "Model: " << "Select model Hamiltonian {cQED,cQEDnoise,Zn,Zncgl,Zncglnoise,cQEDeffective (experimental), realZn (experimental)}" << endl
	 << "N:     " << "Number of spins" << endl
	 << "d:     " << "Bosonic Hilbertspace dimension" << endl
	 << "D:     " << "Bond dimension of the MPS" << endl
	 << "µ:     " << "Dimensionless constant" << endl
	 << "x:     " << "Coupling constant" << endl
	 << "e:     " << "Constant for gauge field" << endl
	 << "la:    " << "Constant for penalty term" << endl
	 << "ig:    " << "True if inital guess for the ground state should be read form HDD (optional)" << endl
	 << "ac:    " << "Accuracy for contractor (optional)" << endl
	 << "n:     " << "Strength for the noise, if cQEDnoise/Zncglnoise model is chosen (optional)"<< endl
	 << "lx:    " << "Constant for the Lx-term, if cQEDeffective model is chosen (optional)"<< endl;
    return -1;
  }
  string model=argv[++cntr];
  int N=atoi(argv[++cntr]);
  int d=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  double mu=atof(argv[++cntr]);
  double x=atof(argv[++cntr]);
  double e=atof(argv[++cntr]);
  double lam=atof(argv[++cntr]);
  bool ig = false;
  bool set_accuracy = false;
  double desired_accuracy;
  double noise_strength=1.0;
  double lxprefactor=1.0;
  string ig_file;
  
  size_t found_dat;
  
  if(argc == NUM_OF_PARAMS + 2)
  {
    ig = true;
    ig_file = argv[++cntr];
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file: " << ig_file << endl;    
  }
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
  if(argc == NUM_OF_PARAMS + 4)
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    noise_strength = atof(argv[++cntr]);
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file " << ig_file << endl;
  }
  if(argc == NUM_OF_PARAMS + 5)
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    noise_strength = atof(argv[++cntr]);
    lxprefactor = atof(argv[++cntr]);
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file " << ig_file << endl;
  }
    
  
  if(argc>=NUM_OF_PARAMS+1){
    cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "N  = " << N << endl
	<< "d  = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "x  = " << x << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl;
  }
  if(argc>=NUM_OF_PARAMS+2){
    cout << "ig = " << ig << endl;
  }
  if(argc>=NUM_OF_PARAMS+3){
     cout << "acc= " << desired_accuracy << endl;
  }
  if(argc>=NUM_OF_PARAMS+4){
     cout << "n =  " << noise_strength << endl;
  }
  if(argc>NUM_OF_PARAMS+4){
     cout << "lx=  " << lxprefactor << endl;
  }
  cout << "-----------" << endl << endl;
    
       
  //Measure time
  time_t start,end;
  time_t start_local,end_local;
  double readtime_gs=0.0,writetime_gs=0.0,writetime_gauss=0.0,writetime_spin=0.0,writetime_ln=0.0,writetime_total=0.0,writetime_data=0.0;
  
  time(&start);
  
  int L=2*N-1; 
  
  //Create an array containing the physical dimensions (alternating between fermion and boson, ending with fermion)
  int *phys_dims = new int[L];
  int *bond_dim = new int[L-1];
  for(int i=0; i<L; ++i){
    phys_dims[i] = 2;
    //cout << "phys_dim[" << i << "]=" << phys_dims[i] << endl;
    if((i+1)<L){
      phys_dims[i+1] = d;
      //cout << "phys_dim[" << i+1 << "]=" << phys_dims[i+1] << endl;
      i++;
    }
  }
  
  for (int i=0; i<L-1; i++){
    bond_dim[i] = D;
    //cout << "Eintrag Nr. " << i << " ist " << phys_dims[i] << endl;
  }
  
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  
  if(set_accuracy){
    cout << "Setting Convergence Tolerance to " << desired_accuracy<<endl<<endl;
    contractor.setConvTol(desired_accuracy);
  }
  
  //Set the seed for the C++ rand, such that I get truly distinct states (in case I want that)
  srand (time(NULL));
  
  //Set the initial state to a random state
  MPS mps(L,bond_dim,phys_dims);
  mps.setRandomState();
  mps.gaugeCond('R',1);
  mps.gaugeCond('L',1);
  
  
  //Get the corresponding Hamiltonian, to avoid lower und upper case issues, convert everything to lower case
  transform(model.begin(),model.end(),model.begin(),::tolower);
  Model mod;
  if(model.compare("zn")==0)
    mod = Zn;
  else if(model.compare("cqed")==0)
    mod = cQED;
  else if(model.compare("cqednoise")==0)
    mod = cQEDnoise;
  else if(model.compare("cqedeffective")==0)
    mod = cQEDeffective;
  else if(model.compare("zncgl")==0)
    mod = Zncgl;
  else if(model.compare("zncglnoise")==0)
    mod = Zncglnoise;
  else if(model.compare("realzn")==0)
    mod = realZn;
  else
  {
    cout << "Unknown model, program will be aborted..." << endl;
    exit(666);
  }
  KogutSusskindHamiltonian Hsngi(N,d,D,mu,x,lam,mod,noise_strength);
  if(model.compare("cqedeffective")==0)
  {
    Hsngi.set_lxprefactor(lxprefactor);
    Hsngi.updateHMPO();
  }
  
  //Override the random state with a deterministic inital state
  MPS *start_state;
  start_state=Hsngi.constructInitalMPS();
  mps = *start_state;
  start_state->clear();
  delete start_state;
  
  //cout << "SchwingerHamiltonianNonGaugeinvariant = " << Hsngi.getHMPO() << endl;
  //cout << "MPS davor " << mps << endl;
  
  //cout << "MPS is left-gauged " << mps.isGaugeL() << endl;
  //cout << "MPS is right-gauged " << mps.isGaugeR() << endl;
  
  //Try to find the groundstate
  const MPO& hamil=Hsngi.getHMPO();  
  double E0 = 0.0;
  complex_t wert;
  
  ifstream inital_guess;
  inital_guess.open(ig_file.c_str());
  cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<        "                   %Starting to calculate Groundstate%                     " << endl
       <<        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
       
  //cout << "Starting to calculate Groundstate" << endl;
  if(inital_guess.is_open() && ig)
  {
    cout << "Found initial guess for the groundstate" << endl;
    time(&start_local);
    cout << "Bond dimension before import: "<< mps.getBond() << endl;
    cout << "Length before import: " << mps.getLength() << endl;
    mps.importMPS(ig_file.c_str());
    cout << "Bond dimension after import: "<< mps.getBond() << endl;
    cout << "Length after import: " << mps.getLength() << endl;
    cout << "Energy after import: " << contractor.contract(mps,hamil,mps) << endl;
    //Sometimes I have problems with corrupted files on the cluster filesystem, which obviously seem readable but give energy=infinity. Therefore I check the norm, after I read something from the filesystem just to make sure it is a valid state and I do not waste time computing nonsense.
    if(abs(contractor.contract(mps,mps))>(1.0+EPSILON))
    {
      //State is somehow corrupted
      cout << endl <<"INTIAL STATE IS CORRUPTED, NORM IS " << endl;
      cout << "Replacing it by a random inital guess" << endl << endl;
      mps.setRandomState();
    }
    
    //In case an intial guess with bond dimension smaller than set is imported it has to be expanded again
    if(mps.getBond() < D)
    {
      /*mps.gaugeCond('R',1);
      mps.gaugeCond('L',1);
      cout << "Energy after gauging: " << contractor.contract(mps,hamil,mps) << endl;*
      cout << "MPS: " << mps << endl;*/
      //Increase bond dimension again, in case the inital guess has a smaller bond dimension
      mps.increaseBondDimension(D);
      //Increase bond dimension again, in case the inital guess has a smaller bond dimension
      //This time with a little noise, maybe better to avoid stability issues?
      //The noise is a total noise for the entire MPS, however I want it to be per site such that I compute the required total noise on the fly
      //mps.increaseBondDimensionWithNoise(D,1.0E-2*(mps.getLength()*(D-mps.getBond())));
      //mps.increaseBondDimensionWithNoise(D,0.0);      
      
    }
    else if(mps.getBond() > D)
    {
      cout << "WARNING, the bond dimension of the inital guess is greater than the maximum bond dimension D=" << D << " set by the user" << endl;
      cout << "Bond dimension will be shrunk" << endl;
      MPS initial_mps;
      initial_mps = MPS(mps,D);
      contractor.optimizeMPS(initial_mps,mps,D);
      mps = initial_mps;
      initial_mps.clear();
      
    }
    cout << "Bond dimension after increase: "<< mps.getBond() << endl;
    time(&end_local);
    readtime_gs = difftime(end_local,start_local);
    cout << "Time for reading GS: " << readtime_gs << " s" << endl;
  }
  
  cout << "Bond dimension which is set: D=" << D << endl;
  contractor.findGroundState(hamil,D,&E0,mps);  
  cout << "Bond dimension of MPS: " << mps.getBond() << endl;
  cout << "Groundstate energy: " << E0 << endl;
  cout << "Groundstate energy via Contractor " << contractor.contract(mps,hamil,mps) << endl;
  //mps.exportForMatlab ("Gs.m");
  
  //String for the name of the groundstate
  string GS_name="GS_";
  //String needed for the conversion
  ostringstream convert;  
  //Build the filename
  convert << N;
  GS_name = GS_name+"N"+convert.str();
  convert.str("");
  convert.clear();
  convert << d;
  GS_name = GS_name+"d"+convert.str(); 
  convert.str("");
  convert.clear();
  convert << x;
  GS_name = GS_name+"x"+convert.str(); 
  convert.str("");
  convert.clear();
  convert << mu;
  GS_name = GS_name+"mu"+convert.str()+".dat"; 
  
  time(&start_local);
  mps.exportMPS(GS_name.c_str());
  time(&end_local);
  writetime_gs = difftime(end_local,start_local);
  cout << "Time for writing GS: " << writetime_gs << " s" << endl;
  
  //Calculate the contribution of the penalty term
  complex_t penalty_energy;
  MPO Penalty(1);

  if(model.compare("zncgl")==0 || model.compare("zncglnoise")==0 || model.compare("realzn")==0)
  {
    MPO PenaltyCgl(1);
    Hsngi.constructCglPenaltyMPO(PenaltyCgl);
    penalty_energy = contractor.contract(mps,PenaltyCgl,mps);
    cout << "Cyclic Penalty energy: "<< penalty_energy << endl;
  }
  else
  {
    Hsngi.constructPenaltyMPO(Penalty,1.0);
    penalty_energy = contractor.contract(mps,Penalty,mps);
    cout << "Penalty energy: "<< penalty_energy << endl;
  }
    
  
  cout << "Ratio Penalty in percent/Groundstate energy: " << penalty_energy.re/E0*100.0<<" %" << endl<<endl; 
  
  //Try to calculate the charge
  MPO chargeop(1);
  Hsngi.constructChargeMPO(chargeop);  
  complex_t charge_0=contractor.contract(mps,chargeop,mps);
  cout << "Charge Groundstate " << charge_0.re << endl; 
  
  //Calculate the therms for the Gauss Law and see if it is fullfilled
  ofstream Lnfile;
  Lnfile.open("Ln.txt", ios_base::app);
  complex_t Lnvalue;
  MPO Ln(1);
  cout << endl; 
  for(int i=1; i<L-1;i += 2){
    Hsngi.constructLopMPO(Ln,i);
    Lnvalue = contractor.contract(mps,Ln,mps);
    cout << "--> Expectation value L("<< (i+1)/2 << ") = " << Lnvalue<< endl;
    time(&start_local);
    Lnfile << setprecision(10) << scientific << Lnvalue.re << "\t";
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
    Gaussvalue = contractor.contract(mps,Gauss,mps);
    cout << "--> Expectation value 0.5*[sigma_z+(-1)^"<< i/2+1 << "] = " << Gaussvalue<< endl;
    time(&start_local);
    Gaussfile << setprecision(10) << scientific << Gaussvalue.re << "\t";
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
    Zcomponent = contractor.contract(mps,Sigmaz,mps);
    cout << "--> Expectation value sigma_z("<< i/2+1 << ") = " << Zcomponent << endl;
    time(&start_local);
    spinfile << setprecision(10) << scientific << Zcomponent.re << "\t";
    time(&end_local);
    writetime_spin += difftime(end_local,start_local);
  }
  time(&start_local);
  spinfile << endl;
  spinfile.close(); 
  time(&end_local);
  writetime_spin += difftime(end_local,start_local);
  
  cout << "Time for writing Spin: " << writetime_spin << " s" << endl;
  
  //Calculate the condensate farction
  complex_t cfraction;
  MPO Condensate(1);
  cout << endl;
  Hsngi.constructCondensateMPO(Condensate);
  cfraction = contractor.contract(mps,Condensate,mps);
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
  
 
  
  //End Time
  time(&end);
  cout << "Run time: " << difftime(end,start) <<"s" <<  endl;  
  writetime_total = writetime_data + writetime_gauss + writetime_gs + writetime_ln + writetime_spin;
  cout << "Time for writing in total: " << writetime_total << " s" << endl;
  cout << "Time on filesystem in total (reading + writing): " << writetime_total+readtime_gs << endl; 
  cout << "End of program" << endl;
  return 0;
}

// bool file_exists(string filename)
// {
//   ifstream file;
//   file.open(filename.c_str());
//   if(file.good())
//     return true;
//   else
//     return false;
// }
