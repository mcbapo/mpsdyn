/**
   \file EstimateTruncationError.cpp
   Given an MPS of bond dimension D, this programm cuts back the bond dimension to a given value in certain steps specified by the user. After after every truncation the observables are recomputed to get an estimate for the truncation error
   
  \param <Dmin> (int) Minimum bond dimension
  \param <dD> (int) Step for bond dimension reduction
  \param <x> (double) coupling constant
  \param <µ> (double) Mass constant of the Hamiltonian
  \param <e> (double) Constant in thte Hamiltonian
  \param <lambda> (double) Constant for the penalty-term in the Hamiltonian
  \param <ig_file> (string) Filename of the MPS file
  \param <acc> (double) Contractor accuracy (optional)
  
   
  \author Stefan Kühn
  \date 16/01/2013

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

#define FILENAME "cQED.txt"
#define NUM_OF_PARAMS 7

using namespace std;
using namespace shrt;

bool file_exists(string filename);
void print_MPS(MPS);
double sum(vector<double>);
double sum_abs(vector<double>);

int main(int argc,const char* argv[]){
  
   //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << " D dD mu x e la ig ac" << endl;
    cout << "D: " << "Minimum bond dimension which should be reached by truncation" << endl
	 << "dD:" << "Step size for bond dimension truncation" << endl
	 << "µ: " << "Dimensionless constant" << endl
	 << "x: " << "Coupling constant" << endl
	 << "e: " << "Constant for gauge field" << endl
	 << "la:" << "Constant for penalty term" << endl
	 << "ig:" << "Filename of the MPS which should be truncated" << endl
	 << "ac:" << "Accuracy for contractor (optional)" << endl;
    return -1;
  }
  
  //Counter for the arguments
  int cntr=0;
  //Minimum bond dimension
  int Dmin=atoi(argv[++cntr]);
  //Stepsize for bond dimension reduction
  int dD=atoi(argv[++cntr]);
  //Mass
  double mu=atof(argv[++cntr]);
  //Coupling constant
  double x=atof(argv[++cntr]);
  //Constant which can be used to manipulate the Hamiltonian
  double e=atof(argv[++cntr]);
  //Strength of the penalty
  double lam=atof(argv[++cntr]);
  //Name for file for initial state
  const char* filename=argv[++cntr];
  string ig_file(filename);
  //Set the accuracy of contractor?
  bool set_accuracy = false;
  //Variable for the desired accuracy in case it should be set
  double desired_accuracy;
  //File for the results
  ofstream resfile;
  
  //Custom accuracy for contractor and file for inital state
  if(argc == NUM_OF_PARAMS + 2)
  {
    set_accuracy = true;
    desired_accuracy = atof(argv[++cntr]);    
  }  
    
  //Print parameter which were set
  if(argc==NUM_OF_PARAMS+1){
    cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "D  = " << Dmin << endl
	<< "dD = " << dD << endl
	<< "µ  = " << mu << endl
	<< "x = " << x << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "ig = " << ig_file << endl
	<< "-----------" << endl << endl;
  }
  else {
    cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "D  = " << Dmin << endl
	<< "dD = " << dD << endl
	<< "µ  = " << mu << endl
	<< "x = " << x << endl
	<< "e  = " << e << endl
	<< "la = " << lam << endl
	<< "ig = " << ig_file << endl
	<< "acc= " << desired_accuracy << endl
	<< "-----------" << endl << endl;
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
    cout << "Setting Convergence Tolerance to " << desired_accuracy<<endl<<endl;
    contractor.setConvTol(desired_accuracy);
  }  
  
  //Read the state
  MPS initial_mps,mps;
  initial_mps.importMPS(ig_file.c_str());
  
  //Get the parameters of the MPS
  int L = initial_mps.getLength();
  int N = (L+1)/2;
  int D = initial_mps.getBond();
  int d = (initial_mps.getA(1)).getd()-1;
  
  cout << "N = " << N << endl
       << "D = " << D << endl
       << "d = " << d << endl;  
  
  //Get the Hamiltonian object for the cQED model
  CQED Hsngi(N,d,D,mu,x,e,lam);  
  //Get Hamiltonian MPO
  const MPO& hamil = Hsngi.getHMPO();  
  complex_t E0;
  
  //Penalty MPO
  MPO Penalty(1);
  Hsngi.constructPenaltyMPO(Penalty);
  complex_t penalty_energy;
  //Electric field MPO and observable
  MPO Ln(1);
  complex_t Lnvalue;
  //Charge MPO and variable for the charge
  const MPO& chargeop=Hsngi.getChargeMPO();  
  complex_t charge_0; 
  ofstream Lnfile;
  //Spin MPO and observables
  MPO Sigmaz(1);
  complex_t Zcomponent;
  ofstream spinfile;
  //Gauss MPO
  MPO Gauss(1);
  complex_t Gaussvalue;
  ofstream Gaussfile;
  //Condensate MPO
  MPO Condensate(1);
  Hsngi.constructCondensateMPO(Condensate);
  complex_t cfraction;
  
  //Monitor the truncation error
  double error;
  
  
  for(int curr_bond=D-dD; curr_bond>=Dmin; curr_bond-= dD)
  {  
    //Build a new MPS out of the old one by truncating in bond dimension
    error = -1.0;
    mps = MPS(initial_mps,curr_bond);
    contractor.optimizeMPS(initial_mps,mps,curr_bond,&error);
    
    cout << endl << "Current bond dimension: " << mps.getBond() << endl;
    cout << "Current norm: " << contractor.contract(mps,mps) << endl;
    cout << "Truncation Error: " << error << endl;
    
//     //String for the name of the groundstate
//     string GS_name="GS_";
//     //String needed for the conversion
//     ostringstream convert;  
//     //Build the filename
//     convert << N;
//     GS_name = GS_name+"N"+convert.str();
//     convert.str("");
//     convert.clear();
//     convert << x;
//     GS_name = GS_name+"x"+convert.str(); 
//     convert.str("");
//     convert.clear();
//     convert << d;
//     GS_name = GS_name+"N0"+convert.str(); 
//     convert.str("");
//     convert.clear();
//     convert << mu;
//     GS_name = GS_name+"mu"+convert.str()+".dat"; 
//     //Save the ground state
//     mps.exportMPS(GS_name.c_str());
    
    //Calculate the energy
    E0 = contractor.contract(mps,hamil,mps);
    
    //Calculate the penalty energy
    penalty_energy = contractor.contract(mps,Penalty,mps);
    
    
    //Try to calculate the charge
    charge_0=contractor.contract(mps,chargeop,mps);
    cout << "Charge Groundstate " << charge_0.re << endl;  
    
    //Calculate the terms for the Gauss Law and see if it is fullfilled
    Lnfile.open("Ln.txt", ios_base::app);
    cout << endl;
    for(int i=1; i<L-1;i += 2){
      Hsngi.constructLopMPO(Ln,i);
      Lnvalue = contractor.contract(mps,Ln,mps);
      cout << "--> Expectation value L("<< (i+1)/2 << ") = " << Lnvalue<< endl;
      Lnfile << Lnvalue.re << "\t";
    }
    Lnfile << endl;
    Lnfile.close();
    
    Gaussfile.open("Gaussfile.txt", ios_base::app);
    cout << endl;
    for(int i=0; i<L;i += 2){
      Hsngi.constructGaussMPO(Gauss,i);
      Gaussvalue = contractor.contract(mps,Gauss,mps);
      cout << "--> Expectation value 0.5*[sigma_z+(-1)^"<< i/2+1 << "] = " << Gaussvalue<< endl;
      Gaussfile << Gaussvalue.re << "\t";
    }
    Gaussfile << endl;
    Gaussfile.close();
    

    
    //Calculate the z-Component of the spins
    spinfile.open("Spin.txt", ios_base::app);
    cout << endl;
    for(int i=0; i<L;i += 2){
      Hsngi.constructsigmazMPO(Sigmaz,i);
      Zcomponent = contractor.contract(mps,Sigmaz,mps);
      cout << "--> Expectation value sigma_z("<< i/2+1 << ") = " << Zcomponent << endl;
      spinfile << Zcomponent.re << "\t";
    }
    spinfile << endl;
    spinfile.close(); 
    
    
    //Calculate the condensate farction
    cout << endl;
    cfraction = contractor.contract(mps,Condensate,mps);
    cout << "Condensate: "<< cfraction << endl << endl;  
    
    //Ouput the results up to now (sometimes the computation of the exicted state(s) fails and the Contractor routine quits the program, therefore all results which were not stored would be lost). To avoid this, save the first part of the results now
    cout << "Saving results..." << endl;
    if(file_exists(FILENAME))
    {
      resfile.open(FILENAME,ios_base::app);  
      resfile << endl << setprecision(10) << N << "\t" << d << "\t" << curr_bond << "\t" << mu << "\t" << x << "\t" << e << "\t" << lam << scientific <<"\t" << E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" <<penalty_energy.re << "\t";
    }
    else
    {
      resfile.open(FILENAME,ios_base::app);  
      resfile << "#N" << "\t" << "d" << "\t" << "D" << "\t" << "mu" << "\t" << "x" << "\t" << "e" << "\t" << "lambda" << "\t"<< "E0" << "\t" << "C0" << "\t" << "Cfrac" << "\t" << "E1" << "\t" << "C1" << "\tE_pen"<< endl;
      resfile << setprecision(10) << N << "\t" << d << "\t" << curr_bond << "\t" << mu << "\t" << x << "\t" << e << "\t" << lam << scientific << "\t"<< E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" << penalty_energy.re << "\t";
    }
    resfile.close();
  
  }  
  
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

double sum(vector<double> numbers)
{
  double res=0.0;
  for (int i=0; i<numbers.size(); i++)
    res += numbers[i];
  return res;
}

double sum_abs(vector<double> numbers)
{
  double res=0.0;
  for (int i=0; i<numbers.size(); i++)
    res += fabs(numbers[i]);
  return res;
}


