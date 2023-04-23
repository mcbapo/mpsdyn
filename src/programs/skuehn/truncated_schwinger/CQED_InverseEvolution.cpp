/**
   \file CQED_InverseEvolution.cpp
   Realisation of a simulation of the time evolution under the fundamental Hamiltonian for 1+1d cQED as proposed in Ignacios review. This version takes an inital state and evolves it back in time. This might be useful for testing and debugging.
   
  \param <ig_file> (string) File from which the initial state should be read from
  \param <param_file> (string) File in which the parameters for the inital state are stored
   
  \author Stefan Kühn
  \date 14/11/2013

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

#define NUM_OF_PARAMS 2
#define FILENAME "INVcQED.txt"
#define RESOLUTION 10
#define EPS 1E-10
#define OFFSET 5192
//Tolerance for the error using the routine to apply an MPO and cut back the bond dimension
#define TOLERANCE 1.0E-6
//Maximum bond dimension which is allowed
#define DMAX 100
//Project back into Gauss Law fulfilling subspace?
#define PROJECTOR_FLAG false

using namespace std;
using namespace shrt;

bool file_exists(string filename);
void print_MPS(MPS);
double sum(vector<double>);
double sum_abs(vector<double>);
int read_params(string paramfile,int &, int &, int &, double &, double &, double &, double &, double &, int &, int &, double &, int &);

int main(int argc,const char* argv[]){
  int cntr=0;
  
  //Case that number of input arguments is to small, display some instruction
  if(argc != NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << " GS_fileparamfile" << endl;
    cout << "GS_file:   File containing the ground state which will be evolved backwards" << endl;
    cout << "paramfile: File containing the corresponding parameters" << endl;
    return -1;
  }
  
  //Number of spin sites
  int N;
  //Total number of sites
  int L;
  //Total number fo bosons living on a link
  int d;
  //Bond dimension
  int D;
  //Mass
  double mu;
  //Coupling constant
  double xs;
  double xe;
  //Constant which can be used to manipulate the Hamiltonian
  double e;
  //Strength of the penalty
  double lam;
  //Size of timestep
  double dt;
  //Total number of timesteps
  int num_of_steps;
  //Desired order of the evolution (by default 1)
  int evolution_order;
  //Name for file for initial state
  string ig_file=argv[1];
  //Name for the parameter file
  string paramfile=argv[2];
  //Flag I set in case I discover a jump
  bool no_jump_flag=true;
  //When did the jump occur
  int ijump;
  //x-value at the jump
  double xjump;    
       
  //Measure time
  time_t start,end;
  time_t start_local,end_local;
  double readtime_gs=0.0,writetime_gs=0.0,writetime_gauss=0.0,writetime_spin=0.0,writetime_ln=0.0,writetime_total=0.0,writetime_data=0.0;
  
  time(&start);
  
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(1.0E-10);
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  
  //Pointer to MPS to hold the inital state which will be evolved in time
  MPS *mps;
  
  //Read the inital state
  mps = new MPS();
  mps->importMPS(ig_file.c_str());
  //Since there is a MPS which I want to evolve backwards, I have to adjust all parameters according to these given in the parameter file
  int status = read_params(paramfile,N,d,D,mu,xs,xe,lam,dt,num_of_steps,evolution_order,xjump,ijump);
  if(status!=0)
  {
    cout << "Error while reading parameter file, program will be aborted" << endl;
    exit(666);
  }
  
  //Since I want to evolve backwards, replace dt by -dt
  dt = -1.0*dt;
  
  //Determine total number of sites
  L = 2*N-1;
  
  
 
  
  //Get the Hamiltonian object for the cQED model
  CQED Hsngi(N,d,D,mu,xjump,e,lam);  
  //Get Hamiltonian MPO
  const MPO& hamil = Hsngi.getHMPO();  
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
  complex_t cfraction, cfraction_old;
  //MPO to determine the condensate
  MPO Condensate(1);
  //Variable for the current x value
  double xvalue=0.0;
  //File for energy values and other observables during the evolution
  ofstream energy_file;
  //Filepointer for the inital guess
  ifstream inital_guess;

  
  
  
  //Determine the inital energy
  cout << "Energy before time evolution " << contractor.contract(*mps,hamil,*mps) << endl;  
  //Determine condensate before evolution
  Hsngi.constructCondensateMPO(Condensate);
  cfraction = contractor.contract(*mps,Condensate,*mps);  
  cout <<"Intial condensate: " << cfraction << " with x-value " << Hsngi.get_x() << endl;
  
  //Try evolve back in time some given intial state 
 
  //Again time evolution via taylor series, this time more optimized since I construct an operator for 1-idtH to only have to optimize once
  if(evolution_order==21)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "     		 %Time evolution with Taylor series%                   " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
 
    energy_file.open("invefile21.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\tTrunc Err\tC_0 (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Variables to monitor the intermediate charge
    //Try to calculate the charge
    const MPO& chargeoperator=Hsngi.getChargeMPO();
    complex_t curr_charge;
    
    //Monitor number of iterations
    int iterations;
    ofstream iteration;
    iteration.open("Iterations.txt");
   
    aux = *mps;
    
    MPO evolutionoperator(1);
    
    
    //Actual time evolution
    for(int i=ijump; i>0; i--)
    {
      mps_old = *mps;
      
      //Get the current value of x
      xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);  
   
      //Get the time evolution operators and the current Hamiltonian
      Hsngi.set_x(xvalue);
      //Hsngi.updateHMPO();
      
      //The first order Taylor expansion of exp(-i dt H) = 1 - i*dt*H + O(dt²), therefore I prepare H|phi> and then use optimizeSum to optimize the result
      Hsngi.EvolutionMPO(evolutionoperator,dt);
      //cout << "EvolutionMPO " << evolutionoperator << endl;
      errN=-1.0;
      contractor.optimize(evolutionoperator,*mps,aux,D,&errN);
      aux.setNormFact(1.0);
      *mps = aux;
      //mps->setNormFact(1.0);
      
      if((i%RESOLUTION)==0)
      {
	Hsngi.updateHMPO();
	//Since I might set the penalty strength in the Hamiltonian to zero but nevertheless want to monitor the penalty energy, I set the constant for the penalty MPO explicitly to one here
	Hsngi.constructPenaltyMPO(Penalty,1.0);
	
	cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
	
	//Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
	curr_overlap = contractor.contract(mps_old,*mps);
	cout << "Overlap " << curr_overlap << "\t";
      
	//curr_energy = contractor.contract(*mps,hamil2,*mps);
	curr_energy = contractor.contract(*mps,hamil,*mps);
	curr_norm = contractor.contract(*mps,*mps);
	
	penalty_energy = contractor.contract(*mps,Penalty,*mps);
	
	Hsngi.constructCondensateMPO(Condensate);
	cfraction = contractor.contract(*mps,Condensate,*mps);
	curr_charge = contractor.contract(*mps,chargeoperator,*mps);
	
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	cout << "Total truncation error " << errtot << endl;
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << scientific << setprecision(10) << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<< "\t" << Hsngi.get_x() << "\t" << errtot << "\t" << curr_charge << endl;
      }		
    }
    energy_file.close();
  }
  
  //Again time evolution via taylor series (using the explicit MPO for 1-itH). In this version I project back in the Gauss Law fulfilling subspace every timestep. 
  if(evolution_order==22)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "% Time evolution with Taylor series and projection into Gauss Law subspace %                   " << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
 
    energy_file.open("invefile22.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\tTrunc Err\tC_0 (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Variables to monitor the intermediate charge
    //Try to calculate the charge
    const MPO& chargeoperator=Hsngi.getChargeMPO();
    complex_t curr_charge;
    
    //Monitor number of iterations
    int iterations;
    ofstream iteration;
    iteration.open("Iterations.txt");
   
    aux = *mps;
    
    MPO evolutionoperator(1);
    MPO projector(1);
    cfraction_old = cfraction;
    
    
    //Actual time evolution
    for(int i=ijump; i>0; i--)
    {
      mps_old = *mps;
      
      //Get the current value of x
      xvalue = ramp_xvalue(xs,xe,num_of_steps-1,i-1);  
   
      //Get the time evolution operators and the current Hamiltonian
      Hsngi.set_x(xvalue);
      //Hsngi.updateHMPO();
      
      //The first order Taylor expansion of exp(-i dt H) = 1 - i*dt*H + O(dt²), therefore I prepare H|phi> and then use optimizeSum to optimize the result
      Hsngi.EvolutionMPO(evolutionoperator,dt);
      //cout << "EvolutionMPO " << evolutionoperator << endl;
      errN=-1.0;
      contractor.optimize(evolutionoperator,*mps,aux,D,&errN);
      aux.setNormFact(1.0);
      *mps = aux;
      //mps->setNormFact(1.0);
      
      //Project back intor the Gauss Law fulfilling subspace if desired
      if(PROJECTOR_FLAG)
      {
	//cout << "Projecting state into Gauss Law fulfilling subspace" << endl;
	Hsngi.GaussProjector2(projector);
	//cout << "Projector " << projector << endl;
	/*for(int k=1; k<=N; k++)
	{
	  Hsngi.GaussProjector2(projector,k);
	  iterations = contractor.optimize(projector,aux,*mps,D);
	  iteration << iterations << endl;
	  /*stringstream ss;//create a stringstream
	  ss << k;//add number to the stream
	  pro_filename = "gpro_";
	  pro_filename.append(ss.str());
	  pro_filename.append(".m");
	  cout << "pro_filename = " << pro_filename<<endl; 
	  projector.exportForMatlab(pro_filename.c_str());*/
	  
	  /*if(fabs(1.0-abs(contractor.contract(aux,*mps)))>1E-5)
		cout << "Overlap after projection " << contractor.contract(aux,*mps) << endl;
	  aux = *mps;		
	}*/
	contractor.optimize(projector,aux,*mps,D);
      }
      //*mps = aux;
      
      if((i%RESOLUTION)==0)
      {
	Hsngi.updateHMPO();
	//Since I might set the penalty strength in the Hamiltonian to zero but nevertheless want to monitor the penalty energy, I set the constant for the penalty MPO explicitly to one here
	Hsngi.constructPenaltyMPO(Penalty,1.0);
	
	cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
	
	//Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
	curr_overlap = contractor.contract(mps_old,*mps);
	cout << "Overlap " << curr_overlap << "\t";
      
	//curr_energy = contractor.contract(*mps,hamil2,*mps);
	curr_energy = contractor.contract(*mps,hamil,*mps);
	curr_norm = contractor.contract(*mps,*mps);
	
	penalty_energy = contractor.contract(*mps,Penalty,*mps);
	
	Hsngi.constructCondensateMPO(Condensate);
	cfraction = contractor.contract(*mps,Condensate,*mps);
	curr_charge = contractor.contract(*mps,chargeoperator,*mps);
	
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	cout << "Total truncation error " << errtot << endl;
	
	
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << scientific << setprecision(10) << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<< "\t" << Hsngi.get_x() << "\t" << errtot << "\t" << curr_charge << endl;
      }		
    }
    energy_file.close();
  }
  
  //Again time evolution via taylor series (using the explicit MPO for 1-itH). In this version I project back in the Gauss Law fulfilling subspace every timestep. Furthermore I do NOT start with x=0 as I do in the other version, but with a finite (though very small) value of x
  if(evolution_order==23)
  {
    cout << endl <<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
       <<          "% Time evolution with Taylor series and projection into Gauss Law subspace starting with finite x %" << endl
       <<          "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
 
    energy_file.open("invefile23.txt");
    energy_file <<"#Time\tEnergy\tNorm\tPenalty\tTrue GS\tCondensate\tx\tTrunc Err\tC_0 (resolved every "<< RESOLUTION << " steps)" << endl;
    
    //Vector to save the errors during truncation
    vector<double> truncation_errors;
    double errN;
    double errtot;
    
    //Variables to monitor the intermediate charge
    //Try to calculate the charge
    const MPO& chargeoperator=Hsngi.getChargeMPO();
    complex_t curr_charge;
    
    //Monitor number of iterations
    int iterations;
    ofstream iteration;
    iteration.open("Iterations.txt");
   
    aux = *mps;
    
    MPO evolutionoperator(1);
    MPO projector(1);
    
    
    //Actual time evolution
    for(int i=ijump; i>0; i--)
    {
      mps_old = *mps;
      
      //Get the current value of x
      xvalue = ramp_xvalue(xs,xe,num_of_steps,i);  
   
      //Get the time evolution operators and the current Hamiltonian
      Hsngi.set_x(xvalue);
      //Hsngi.updateHMPO();
      
      //The first order Taylor expansion of exp(-i dt H) = 1 - i*dt*H + O(dt²), therefore I prepare H|phi> and then use optimizeSum to optimize the result
      Hsngi.EvolutionMPO(evolutionoperator,dt);
      //cout << "EvolutionMPO " << evolutionoperator << endl;
      errN=-1.0;
      contractor.optimize(evolutionoperator,*mps,aux,D,&errN);
      aux.setNormFact(1.0);
      *mps = aux;
      //mps->setNormFact(1.0);
      
      //Project back intor the Gauss Law fulfilling subspace if desired
      if(PROJECTOR_FLAG)
      {
	//cout << "Projecting state into Gauss Law fulfilling subspace" << endl;
	Hsngi.GaussProjector2(projector);
	contractor.optimize(projector,aux,*mps,D);
      }
      
      if((i%RESOLUTION)==0)
      {
	Hsngi.updateHMPO();
	//Since I might set the penalty strength in the Hamiltonian to zero but nevertheless want to monitor the penalty energy, I set the constant for the penalty MPO explicitly to one here
	Hsngi.constructPenaltyMPO(Penalty,1.0);

	cout<< "x-value for data calculation: " << Hsngi.get_x() << endl;
	
	//Monitor norm and energy as well as overlap of the current state with the state at the beginning of the timestep
	curr_overlap = contractor.contract(mps_old,*mps);
	cout << "Overlap " << curr_overlap << "\t";
      
	//curr_energy = contractor.contract(*mps,hamil2,*mps);
	curr_energy = contractor.contract(*mps,hamil,*mps);
	curr_norm = contractor.contract(*mps,*mps);
	
	penalty_energy = contractor.contract(*mps,Penalty,*mps);
	
	Hsngi.constructCondensateMPO(Condensate);
	cfraction = contractor.contract(*mps,Condensate,*mps);
	curr_charge = contractor.contract(*mps,chargeoperator,*mps);
	
	errtot = sum_abs(truncation_errors);
	truncation_errors.clear();
	
	cout << "Total truncation error " << errtot << endl;
	
	/*vector<MPS*> excitedstates;
	double lambdak;
	MPS first_ex(L,bond_dim,phys_dims);
	excitedstates.push_back(mps);

	contractor.findNextExcitedState(hamil2,D,excitedstates,&lambdak,first_ex); 
	cout << "Excited state energy: " << lambdak << endl<< endl;*/
	
	cout<< "Norm in step " << i << " is " << contractor.contract(*mps,*mps) << "\tNorm aux is " << contractor.contract(aux,aux) << "\tEnergy " << curr_energy << endl;
	//energy_file << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << lambdak << endl;
	energy_file << scientific << setprecision(10) << i*dt << "\t" << curr_energy.re << "\t" << curr_norm.re << "\t" << penalty_energy << "\t" << E0 <<"\t" << cfraction/sqrt(Hsngi.get_x())<< "\t" << Hsngi.get_x() << "\t" << errtot << "\t" << curr_charge << endl;
      }		
    }
    energy_file.close();
  }
  
  
  E0 = curr_energy.re;
//   //String for the name of the groundstate
//   string GS_name="GS_";
//   //String needed for the conversion
//   ostringstream convert;  
//   //Build the filename
//   convert << N;
//   GS_name = GS_name+"N"+convert.str();
//   convert.str("");
//   convert.clear();
//   convert << x;
//   GS_name = GS_name+"x"+convert.str(); 
//   convert.str("");
//   convert.clear();
//   convert << mu;
//   GS_name = GS_name+"mu"+convert.str()+".dat"; 
//   
//   time(&start_local);
//   mps.exportMPS(GS_name.c_str());
//   time(&end_local);
//   writetime_gs = difftime(end_local,start_local);
//   cout << "Time for writing GS: " << writetime_gs << " s" << endl;
   
  
  mps->exportForMatlab("MPS.m");
  
  
  
//   
//   cout << "Ratio Penalty in percent/Groundstate energy: " << penalty_energy.re/E0*100.0<<" %" << endl<<endl; 
  
  //Try to calculate the charge
  const MPO& chargeop=Hsngi.getChargeMPO();  
  complex_t charge_0=contractor.contract(*mps,chargeop,*mps);
  cout << "Charge Groundstate " << charge_0.re << endl;  
  
  //Calculate the terms for the Gauss Law and see if it is fullfilled
  ofstream Lnfile;
  Lnfile.open("InvLn.txt", ios_base::app);
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
  Gaussfile.open("InvGaussfile.txt", ios_base::app);
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
  spinfile.open("InvSpin.txt", ios_base::app);
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
    resfile << endl << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << xe << "\t" << e << "\t" << lam << scientific <<"\t" << E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" <<penalty_energy.re << "\t";
  }
  else
  {
    time(&start_local);
    resfile.open(FILENAME,ios_base::app);  
    resfile << "#N" << "\t" << "d" << "\t" << "D" << "\t" << "mu" << "\t" << "x" << "\t" << "e" << "\t" << "lambda" << "\t"<< "E0" << "\t" << "C0" << "\t" << "Cfrac" << "\t" << "E1" << "\t" << "C1" << "\tE_pen"<< endl;
    resfile << setprecision(10) << N << "\t" << d << "\t" << D << "\t" << mu << "\t" << xe << "\t" << e << "\t" << lam << scientific << "\t"<< E0 << "\t" << charge_0.re << "\t" << cfraction.re << "\t" << penalty_energy.re << "\t";
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
  mps->clear();
  delete mps;
  
 
  
  //End Time
  time(&end);
  double runtime_total = difftime(end,start);
  cout << "Run time: " << runtime_total <<"s" <<  endl;  
  writetime_total = writetime_data + writetime_gauss + writetime_gs + writetime_ln + writetime_spin;
  cout << "Time for writing in total: " << writetime_total << " s" << endl;
  cout << "Time on filesystem in total (reading + writing): " << writetime_total+readtime_gs << endl; 
  
  //Simply uncomment in case benchmarking file should be written
  ofstream benchmarking;
  if(file_exists("Benchmark.txt"))
  {
    benchmarking.open("Benchmark.txt",ios_base::app);
  }
  else
  {
    benchmarking.open("Benchmark.txt");
    benchmarking << "#N\tD\tN_0\tT_total(s)\tT_perstep(s)" << endl;
  }
  benchmarking << std::resetiosflags(std::ios::scientific) << N << "\t" << D << "\t" << d << scientific << setprecision(10) <<"\t" << runtime_total << "\t" << runtime_total/num_of_steps << endl;
  benchmarking.close();
  
  
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
int read_params(string paramfile,int &N, int &d, int &D, double &mu, double &xs, double &xe, double &lam, double &dt, int &steps, int &evolution_order, double &xjump, int &ijump)
{
  ifstream pfile;
  string dummy;
  
  cout << "Trying to read file " << paramfile << endl;
  if(!file_exists(paramfile))
    return -1;
  else
  {
    pfile.open(paramfile.c_str());    
    pfile >> dummy; pfile >> dummy; pfile >> N;
    pfile >> dummy; pfile >> dummy; pfile >> d;
    pfile >> dummy; pfile >> dummy; pfile >> D;
    pfile >> dummy; pfile >> dummy; pfile >> mu;
    pfile >> dummy; pfile >> dummy; pfile >> xs;
    pfile >> dummy; pfile >> dummy; pfile >> xe;
    pfile >> dummy; pfile >> dummy; pfile >> lam;
    pfile >> dummy; pfile >> dummy; pfile >> dt;
    pfile >> dummy; pfile >> dummy; pfile >> steps;
    pfile >> dummy; pfile >> dummy; pfile >> evolution_order;
    pfile >> dummy; pfile >> dummy; pfile >> xjump;
    pfile >> dummy; pfile >> dummy; pfile >> ijump;
    
    cout << "Read N=" << N << endl;
    cout << "Read d=" << d<< endl;
    cout << "Read D=" << D<< endl;
    cout << "Read mu=" << mu<< endl;
    cout << "Read xs=" << xs<< endl;
    cout << "Read xe=" << xe<< endl;
    cout << "Read lamda=" << lam<< endl;
    cout << "Read dt=" << dt<< endl;
    cout << "Read steps=" << steps<< endl;
    cout << "Read evo_order=" << evolution_order<< endl;
    cout << "Read xjump=" << xjump<< endl;
    cout << "Read ijump=" << ijump<< endl;
    
    pfile.close();
    return 0;
  }
}