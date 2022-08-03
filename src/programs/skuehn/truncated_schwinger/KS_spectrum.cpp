/**
   \file KS_spectrum.cpp
   Program to compute the spectrum of one of the Hamiltonians realized by KogutSusskindHamiltonian.cpp
   
  \param <N> (int) Number of Fermions
  \param <d> (int) Physical dimension of Bosons in between
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <x> (double) Coupling constant
  \param <e> (double) Constant of the Operator \f$ \Theta(n)\f$
  \param <lambda> (double) Constant for the penalty-term in the Hamiltonian
   
  \author Stefan Kühn
  \date 20/05/2014

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

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include "misc.h"

#include "KogutSusskindHamiltonian.h"

#define NUM_OF_PARAMS 10
#define FILENAME "KS.txt"

#define THRESHOLD 1E-3

using namespace std;

int main(int argc,const char* argv[]){
  int cntr=0;
  if((argc < NUM_OF_PARAMS + 1) || (argc > NUM_OF_PARAMS + 6)){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << " model N d D mu x e la" << endl;
    cout << "Model: " << "Select model Hamiltonian {cQED,Zn,Zncgl,cQEDnoise,Zncglnoise}" << endl
         << "N:     " << "Number of spins" << endl
	 << "N0:    " << "Total number of bosons" << endl
	 << "D:     " << "Bond dimension of the MPS" << endl
	 << "µ:     " << "Dimensionless constant" << endl
	 << "xs:    " << "Coupling constant to start with" << endl
	 << "xe:    " << "Coupling constant to end with" << endl
	 << "dx:    " << "Stepping in coupling constant" << endl
	 << "e:     " << "Constant for gauge field" << endl
	 << "la:    " << "Constant for penalty term" << endl
	 << "ig:    " << "True if inital guess for the ground state should be read form HDD (optional)" << endl
	 << "ac:    " << "Accuracy for contractor (optional)" << endl
	 << "lm:    " << "Mode for lambda, if an integer different from zero is given, lambda will be adjusted relatively to x" << endl
	 << "ne:    " << "Number of excited states that should be calculated (optional by default only one is calculated)" << endl
	 << "se:    " << "Desired sector, an integer different from zero restricts the computation to the Gauss Law fulfilling sector (in case a suitable penalty is set)" << endl;
    return -1;
  }
  
  //Model
  string model=argv[++cntr];
  //Number of fermionic sites
  int N=atoi(argv[++cntr]);
  //Number of bosons on a link
  int d=atoi(argv[++cntr]);
  //Bond dimension
  int D=atoi(argv[++cntr]);
  //Mass
  double mu=atof(argv[++cntr]);
  //Inital coupling strength
  double xs=atof(argv[++cntr]);
  //Coupling to end with
  double xe=atof(argv[++cntr]);
  //Steps in which the coupling constant should be increased
  double dx=atof(argv[++cntr]);
  //Constant in the Hamiltonian (I think it is not needed anymore and should be removed?)
  double e=atof(argv[++cntr]);
  //Constant for the penalty
  double lam_constant=atof(argv[++cntr]);
  //Actual penalty strength
  double lam;
  //Number of excited states which should be calculated
  int num_of_excited_states = 1;
  //Determine if lambda should be fixed or adjusted relatively to the value of x
  bool lambda_is_fixed = true;
  bool ig = false;
  bool set_accuracy = false;
  double desired_accuracy;
  //Restrict computation to Gauss Law fulfilling sector
  bool gl_sector=false;
  string ig_file;
  
  //Variable to find out of specific file is given or just the generic file "GS_guess.dat" should be read
  size_t found_dat;
  
  //Inital guess is given
  if(argc == NUM_OF_PARAMS + 2)
  {
    ig = true;
    ig_file = argv[++cntr];
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file: " << ig_file << endl;    
  }
  //Inital guess + accuracy is given
  else if(argc == NUM_OF_PARAMS + 3)
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
  //Inital guess + accuracy + mode for lambda is given
  else if(argc == NUM_OF_PARAMS + 4)
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    if(atoi(argv[++cntr]) != 0)
      lambda_is_fixed  = false;
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file " << ig_file << endl;
  }
  //Inital guess + accuracy + mode for lambda + number of excited states which have to be calculated is given
  else if(argc == NUM_OF_PARAMS + 5)
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    if(atoi(argv[++cntr]) != 0)
      lambda_is_fixed  = false;
    num_of_excited_states = atoi(argv[++cntr]);
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file " << ig_file << endl;
  }
  //Inital guess + accuracy + mode for lambda + number of excited states which have to be calculated + restriction to gauge invariant sector is given
  else
  {
    ig = true;
    set_accuracy = true;
    ig_file = argv[++cntr];
    desired_accuracy = atof(argv[++cntr]);
    if(atoi(argv[++cntr]) != 0)
      lambda_is_fixed  = false;
    num_of_excited_states = atoi(argv[++cntr]);
    found_dat = ig_file.find(".dat");
    if(found_dat == std::string::npos)
      ig_file = "GS_guess.dat";
    cout << "IG file " << ig_file << endl;
    if(atoi(argv[++cntr]) != 0)
      gl_sector = true;
  } 
    
  
  if(argc==NUM_OF_PARAMS+1){
    cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "Model: " << model << endl
	<< "N  = " << N << endl
	<< "N0 = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
	<< "e  = " << e << endl
	<< "la = " << lam_constant << endl
	<< "-----------" << endl << endl;
  }
  else if (argc==NUM_OF_PARAMS+2){
    cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "Model: " << model << endl
	<< "N  = " << N << endl
	<< "N0 = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
	<< "e  = " << e << endl
	<< "la = " << lam_constant << endl
	<< "ig = " << ig << endl
	<< "-----------" << endl << endl;
	
  }
  else if (argc==NUM_OF_PARAMS+3){
     cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "Model: " << model << endl
	<< "N  = " << N << endl
	<< "N0 = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
	<< "e  = " << e << endl
	<< "la = " << lam_constant << endl
	<< "ig = " << ig << endl
	<< "acc= " << desired_accuracy << endl
	<< "-----------" << endl << endl;
  }
  else if (argc==NUM_OF_PARAMS+4){
     cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "Model: " << model << endl
	<< "N  = " << N << endl
	<< "N0 = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
	<< "e  = " << e << endl
	<< "la = " << lam_constant << endl
	<< "ig = " << ig << endl
	<< "acc= " << desired_accuracy << endl
	<< "lm = " << lambda_is_fixed << endl
	<< "-----------" << endl << endl;
  }
  else if (argc==NUM_OF_PARAMS+5){
     cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "Model: " << model << endl
	<< "N  = " << N << endl
	<< "N0 = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
	<< "e  = " << e << endl
	<< "la = " << lam_constant << endl
	<< "ig = " << ig << endl
	<< "acc= " << desired_accuracy << endl
	<< "lm = " << lambda_is_fixed << endl
	<< "ne = " << num_of_excited_states << endl
	<< "-----------" << endl << endl;
  }
  else{
     cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "Model: " << model << endl
	<< "N  = " << N << endl
	<< "N0 = " << d << endl
	<< "D  = " << D << endl
	<< "µ  = " << mu << endl
	<< "xs = " << xs << endl
	<< "xe = " << xe << endl
	<< "e  = " << e << endl
	<< "la = " << lam_constant << endl
	<< "ig = " << ig << endl
	<< "acc= " << desired_accuracy << endl
	<< "lm = " << lambda_is_fixed << endl
	<< "ne = " << num_of_excited_states << endl
	<< "se = " << gl_sector << endl
	<< "-----------" << endl << endl;
  }    
       
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
    cout << "phys_dim[" << i << "]=" << phys_dims[i] << endl;
    if((i+1)<L){
      phys_dims[i+1] = d+1;
      cout << "phys_dim[" << i+1 << "]=" << phys_dims[i+1] << endl;
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
    cout << "Setting Convergence Tolareance to " << desired_accuracy<<endl<<endl;
    contractor.setConvTol(desired_accuracy);
  } 
  
  //MPSs for the gauge and non gauge invariant case
  MPS mps(L,bond_dim,phys_dims);
  mps.setRandomState();
  mps.gaugeCond('R',1);
  mps.gaugeCond('L',1);

  //Try to find the groundstate
  double E0 = 0.0,E0_nopenalty=0.0;
  complex_t wert;
  vector<MPS*> excitedstates;
  double lambdak;
  MPS *first_ex;
  
  ifstream inital_guess;
  inital_guess.open(ig_file.c_str());
  cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<        "                   %Starting to calculate Groundstate%                     " << endl
       <<        "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
       
  if(inital_guess.is_open() && ig)
  {
    cout << "Found initial guess for the groundstate" << endl;
    time(&start_local);
    cout << "Bond dimension before import: "<< mps.getBond() << endl;
    cout << "Length before import: " << mps.getLength() << endl;
    mps.importMPS(ig_file.c_str());
    cout << "Bond dimension after import: "<< mps.getBond() << endl;
    cout << "Length after import: " << mps.getLength() << endl;
    //In case an intial guess with bond dimension smaller than set is imported it has to be expanded again
    if(mps.getBond() < D)
    {
      //Increase bond dimension again, in case the inital guess has a smaller bond dimension
      //mps.increaseBondDimension(D);
      //Increase bond dimension again, in case the inital guess has a smaller bond dimension
      //This time with a little noise, maybe better to avoid stability issues?
      mps.increaseBondDimensionWithNoise(D,1E-2);
      
    }
    else if(mps.getBond() > D)
    {
      cout << "Error, the bond dimension of the inital guess is greater than the maximum bond dimension D=" << D << " set by the user" << endl;
      cout << "Bond dimension cannot be shrunk" << endl;
      exit(1);
    }
    cout << "Bond dimension after increase: "<< mps.getBond() << endl;
    time(&end_local);
    readtime_gs = difftime(end_local,start_local);
    cout << "Time for reading GS: " << readtime_gs << " s" << endl;
  }
  
  lam = lam_constant;
  double overlap;
  
  //Get the corresponding Hamiltonian, to avoid lower und upper case issues, convert everything to lower case
  transform(model.begin(),model.end(),model.begin(),::tolower);
  Model mod;
  if(model.compare("zn")==0)
    mod = Zn;
  else if(model.compare("cqed")==0)
    mod = cQED;
  else if(model.compare("cqednoise")==0)
    mod = cQEDnoise;
  else if(model.compare("zncgl")==0)
    mod = Zncgl;
  else if(model.compare("zncglnoise")==0)
    mod = Zncglnoise;
  else
  {
    cout << "Unknown model, program will be aborted..." << endl;
    exit(666);
  }
  cout << "WARNING: Until now no noise is supportet, will set it to zero" << endl;
  KogutSusskindHamiltonian Hsngi(N,d,D,mu,xs,lam,mod,0.0);
  const MPO& hamil=Hsngi.getHMPO();   
  MPO chargeop(1);
  Hsngi.constructChargeMPO(chargeop);
  complex_t charge_0; 
  MPO Penalty(1);
 //Construct the penalty MPO according to the chosen model(with fixed strength)
  if((mod==Zncgl)||(mod==Zncglnoise))
    Hsngi.constructCglPenaltyMPO(Penalty,1.0);
  else
    Hsngi.constructPenaltyMPO(Penalty,1.0);
  complex_t penalty_energy;
  
  //Set inital state in the right charge subsector
  MPS* init;
  init = Hsngi.constructInitalMPS();
  mps = *init;
  
  for (double x=xs; x<=xe; x+=dx)
  {
    //Update the paramater in the Hamiltonian
    if(!lambda_is_fixed)
    {
      lam = x*lam_constant;
      Hsngi.set_lambda(lam);
    }   
    Hsngi.set_x(x);
    //Update the Hamiltonian MPO according to the new parameters
    Hsngi.updateHMPO();
    
    cout << endl;
    contractor.findGroundState(hamil,D,&E0,mps);   
    penalty_energy = contractor.contract(mps,Penalty,mps);
    
    cout << "Groundstate energy:                          " << E0 << endl;
    cout << "Groundstate energy via Contractor:           " << contractor.contract(mps,hamil,mps) << endl;
    cout << "Penalty contribution:                        " << penalty_energy << endl;
    cout << "Ratio Penalty in percent/Groundstate energy: " << penalty_energy.re/E0*100.0<<" %" << endl; 
    //Calculate the charge
    charge_0=contractor.contract(mps,chargeop,mps);
    cout << "Charge Groundstate:                          " << charge_0 << endl;  
    //Calculate the condensate farction
    complex_t cfraction;
    MPO Condensate(1);
    Hsngi.constructCondensateMPO(Condensate);
    cfraction = contractor.contract(mps,Condensate,mps);
    cout << "Condensate:                                  "<< cfraction << endl << endl;  
    
    
    //String for the name of the groundstate
    string GS_name="GS_";
    //String needed for the conversion
    ostringstream convert;  
    //Build the filename
    convert << N;
    GS_name = GS_name+"N"+convert.str();
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
    //cout << "Time for writing GS: " << writetime_gs << " s" << endl;      
    
    
    
    //Calculate the therms for the Gauss Law and see if it is fullfilled
    ofstream Lnfile;
    Lnfile.open("Ln.txt", ios_base::app);
    complex_t Lnvalue;
    MPO Ln(1);
    for(int i=1; i<L-1;i += 2){
      Hsngi.constructLopMPO(Ln,i);
      Lnvalue = contractor.contract(mps,Ln,mps);
      //cout << "--> Expectation value L("<< (i+1)/2 << ") = " << Lnvalue<< endl;
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
    
    //cout << "Time for writing Ln: " << writetime_ln << " s" << endl;
    
    ofstream Gaussfile;
    Gaussfile.open("Gaussfile.txt", ios_base::app);
    complex_t Gaussvalue;
    MPO Gauss(1);
    for(int i=0; i<L;i += 2){
      Hsngi.constructGaussMPO(Gauss,i);
      Gaussvalue = contractor.contract(mps,Gauss,mps);
      //cout << "--> Expectation value 0.5*[sigma_z+(-1)^"<< i/2+1 << "] = " << Gaussvalue<< endl;
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
    
    //cout << "Time for writing Gauss: " << writetime_gauss << " s" << endl;
    
    //Calculate the z-Component of the spins
    ofstream spinfile;
    spinfile.open("Spin.txt", ios_base::app);
    complex_t Zcomponent;
    MPO Sigmaz(1);
    for(int i=0; i<L;i += 2){
      Hsngi.constructsigmazMPO(Sigmaz,i);
      Zcomponent = contractor.contract(mps,Sigmaz,mps);
      //cout << "--> Expectation value sigma_z("<< i/2+1 << ") = " << Zcomponent << endl;
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
    
    //cout << "Time for writing Spin: " << writetime_spin << " s" << endl;
    
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
    
    //cout << "Time for writing Data: " << writetime_data << " s" << endl;  
    
    
    ofstream spectrumfile;
    if(file_exists("SpectrumKS.txt"))
    {
      spectrumfile.open("SpectrumKS.txt",ios_base::app);
    }
    else
    { 
      spectrumfile.open("SpectrumKS.txt",ios_base::app);
      spectrumfile << "#x\tE0 \tE0_pen\tE1 \tE1_pen ..." << endl;
    }
    spectrumfile <<  scientific << setprecision(10) << x << "\t" << E0 << "\t" << penalty_energy << "\t" ;
    for(int k=1; k<=num_of_excited_states; k++)
    {
      first_ex = new MPS(L,bond_dim,phys_dims); 
      
      //Add groundstate in case it is the first round
      if(k==1)
	excitedstates.push_back(&mps);
      
      //Make output for debugging
      for(int iter=0; iter < excitedstates.size(); iter++)
	cout << "Exstate["<<iter<<"] with energy " << contractor.contract(*(excitedstates[iter]),hamil,*(excitedstates[iter])) << " and penalty "  <<  contractor.contract(*(excitedstates[iter]),Penalty,*(excitedstates[iter])) << " and charge " << contractor.contract(*(excitedstates[iter]),chargeop,*(excitedstates[iter])) << endl;
      
      //Set inital guess (do not take old state since we look for something orthogonal!!!)
      first_ex->setRandomState();
      cout << endl << "Computing excited state number " << k << endl;
      //In case I only want Gauss states in the Gauss law fulfilling sector
      if(gl_sector)
      {
	
	while(1){
	  contractor.findNextExcitedState(hamil,D,excitedstates,&lambdak,*first_ex);
	  penalty_energy = contractor.contract(*first_ex,Penalty,*first_ex);
	  if(penalty_energy.re > THRESHOLD)
	  {
	    //Obviously I could not find something gauge invariant, therefore I try to set a larger penalty to shift everything which is detrimental away
	    Hsngi.set_lambda(2.0*Hsngi.get_lambda());
	    Hsngi.updateHMPO();
	    cout << "Failed to compute Gauss Law fulfilling excited state, trying again with larger penalty " << Hsngi.get_lambda() << endl;
	  }
	  else
	    break;
	}
      }
      //Otherwise I take all states
      else
	contractor.findNextExcitedState(hamil,D,excitedstates,&lambdak,*first_ex);      
      
      excitedstates.push_back(first_ex);
      cout << "Excited state energy:  " << lambdak << endl;
      cout << "Excites state charge:  " << contractor.contract(*first_ex,chargeop,*first_ex) << endl;
      cout << "Excited state penalty: " << penalty_energy << endl<< endl;
      complex_t charge_1=contractor.contract(*first_ex,chargeop,*first_ex);
      
      //Save results for the exicted state(s)
      resfile.open(FILENAME,ios_base::app);  
      resfile << lambdak << "\t" << charge_1.re << "\t";
      resfile.close(); 
      
      //Write to dedicated spectrumfile
      spectrumfile << lambdak << "\t" << penalty_energy << "\t";
    }
    spectrumfile << endl;
    spectrumfile.close();
    
    //I want to identify if the excited states are really degenerate or if i simply find a bunch of them twice, therefore I calculate overlaps between them
    ofstream overlaps;
    if (!file_exists("Overlaps.txt"))
    {
      overlaps.open("Overlaps.txt");
      for(int k=1; k<=num_of_excited_states; k++)
	for(int l=k+1; l<=num_of_excited_states;l++)
	  overlaps<< "|<Psi_"<<k<<"|Psi_"<<l<<">|"<<"\t"; 
      overlaps << endl;
    }
    else
    {
      overlaps.open("Overlaps.txt",ios_base::app);
    }
      
    for(int k=1; k<=num_of_excited_states; k++)
    {
      for(int l=k+1; l<=num_of_excited_states;l++)
      {
	overlap = abs(contractor.contract(*(excitedstates[k]),*(excitedstates[l])));
	cout<< "|<Psi_"<<k<<"|Psi_"<<l<<">|="<< overlap << endl;
	overlaps << scientific << setprecision(6) <<overlap << "\t";
      }
    }
    overlaps << endl;
    overlaps.close();
    
    //Free the excited states (Entry at zero is the GS and is still needed)
    for(int k=1; k<excitedstates.size(); k++)
    {
      cout << "Clearing state " << k << endl;
      excitedstates[k]->clear();
      delete excitedstates[k];
    }
    excitedstates.clear();
    
     
  }
  
  //Free Memory
  delete[] phys_dims;
  delete[] bond_dim; 
  
  init->clear();
  delete init;
  
 
  
  //End Time
  time(&end);
  cout << "Run time: " << difftime(end,start) <<"s" <<  endl;  
  writetime_total = writetime_data + writetime_gauss + writetime_gs + writetime_ln + writetime_spin;
  cout << "Time for writing in total: " << writetime_total << " s" << endl;
  cout << "Time on filesystem in total (reading + writing): " << writetime_total+readtime_gs << endl; 
  cout << "End of program" << endl;
  return 0;
}