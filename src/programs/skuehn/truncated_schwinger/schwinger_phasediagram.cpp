/**
  \file schwinger_phasediagram.cpp
  Program to compute the phase diagram for the (truncated) Schwinger Model with background field
   
  \param <file> (string) Name of the driver file containing the parameters
  
  \author Stefan Kühn
  \date 05/22/2016

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

#include "Properties.h"
#include "SchwingerHamiltonian.h"
#include "SchwingerHamiltonianSz.h"
#include "KSHamiltonian.h"

#define NUM_OF_PARAMS 1
#define FILENAME "schwinger_pd.txt"
#define FILENAME_TRUNC "schwinger_truncated_pd.txt"
#define TMP_FILE "GS.tmp"



using namespace std;
using namespace shrt;

//Function to do some basic second order trotterized odd even time evolution
void compute_evolution(MPO &exp_even,MPO &exp_odd,MPO &exp_odd_half,Contractor &contractor, MPS &mps,int steps,bool normalize); 
void getSchwingerInitState(MPS &mps,int N,InitState state);
string getFilename(string model, int level, int dlink, int N, double x, double mg, int D);

int main(int argc,const char* argv[])
{  
  //Code Revision (magic field, automatically filled in by subversion)
  string revision="$Rev: 1542 $";
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
  Properties props(inputfile.c_str());
  
  //System size
  int N = props.getIntProperty("N");
  //Maximum bond dimension
  int D = props.getIntProperty("D");
  //Hamiltonian parameters
  double x=props.getDoubleProperty("x");
  double mg=props.getDoubleProperty("mg");  
  int dlink=props.getIntProperty("dlink");
  double alpha=props.getDoubleProperty("alpha");
  //From the ratio m/g and x compute the dimensionless mass
  double mu=2*mg*sqrt(x);
  //Penalty for unphysical/non Gauss Law fulfilling states
  double lambda=props.getDoubleProperty("lambda");
  //Decide which model to run
  string model=props.getProperty("model");
  //Tolerance for the contractor
  double tol=props.getDoubleProperty("tol");
  //Time step size and number of time steps
  double dt=props.getDoubleProperty("dt");
  int steps=props.getIntProperty("steps");
  //Type of start state
  string startstate=props.getProperty("startstate");
  //Noiselevel if there should be some noise involved
  double noiselevel = props.getDoubleProperty("noiselevel");
  string outputname;
  tol=tol==-1? 1.0E-8 : tol;
  noiselevel=noiselevel==-1? 1.0E-1 : noiselevel;
  //File for the results
  ofstream file;
  
  //Seed the random number generator
  srand (time(NULL));
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;  
  //In case PRIMME-support is possible, set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  //Set the convergence tolerance
  contractor.setConvTol(tol);   
  
  //Variables for the expectation value fo the observalbles
  complex_t penaltyval, gamma5val, gammaAval, condensateval, momentumval, chargeval;
  //Variable for the energies of the excite state and the ground state
  double E0,E1;
  
  //MPS for the ground state and the first excited state and vector to save computed states
  MPS *mps;
  vector <MPS*> computed_states;  
  
  //Gamma5 MPO 
  MPO gamma5(1);
  //GammaA MPO
  MPO gammaA(1);
  //Chiral Condensate MPO
  MPO condensate(1);
  //Momentum MPO
  MPO momentum(1);
  
  if(model.compare("truncated") == 0)
  {
    //Get the Hamiltonian
    KSHamiltonian KS(N,dlink,mu,x,alpha,lambda);  
    const MPO& hamil=KS.getHMPO();    
    
    KS.constructGamma5MPO(gamma5);    
    KS.constructGammaAlphaMPO(gammaA);    
    KS.constructCondensateMPO(condensate);    
    KS.constructMomentumMPO(momentum);
    //Penalty MPO to check if Gauss Law is fulfilled
    MPO penalty(1);
    KS.constructPenaltyMPO(penalty,1.0);
    //Total charge MPO
    MPO charge(1);
    KS.constructChargeMPO(charge); 
    
    //Now compute the ground state and the first excited state  
    cout << "---------------------------------------------------------------------"<<endl;
    cout << "---------------------------Truncated Schwinger-----------------------"<<endl;
    cout << "---------------------------------------------------------------------"<<endl;       
    if(dt<0 || steps<0)
    {
      cout << "--> Computing ground state for truncated Schwinger Hamiltonian" << endl;
      mps = new MPS(N,D,2);
      if(startstate.compare("sc") == 0)
	KS.constructInitialMPS(*mps,is_updown);
      else if(startstate.compare("bg") == 0)
	KS.constructInitialMPS(*mps,is_bgfield);      
      else if(startstate.compare("superpos") == 0)
	KS.constructInitialMPS(*mps,is_bgfield_superpos);
      else if(startstate.compare("superpos_noise") == 0)
      {
	KS.constructInitialMPS(*mps,is_bgfield_superpos);
	mps->increaseBondDimensionWithNoise(D,noiselevel);
      }
      else if(startstate.find(".dat") != string::npos && file_exists(startstate))
      {
	cout << "Reading state from file " << startstate << endl;
	mps->importMPS(startstate.c_str());
      }
      else      
      {
	cout << "Using random state" << endl; 
	KS.constructInitialMPS(*mps);
	mps->increaseBondDimension(D);
	mps->setRandomState();
      }
      contractor.findGroundState(hamil,D,&E0,*mps);
      outputname=getFilename("truncated",0,dlink,N,x,mg,D);
      //mps->exportMPS("GS_truncated.dat");
      mps->exportMPS(outputname.c_str());
      computed_states.push_back(mps);
      //Compute some observables
      penaltyval = contractor.contract(*mps,penalty,*mps);
      gamma5val = contractor.contract(*mps,gamma5,*mps);
      gammaAval = contractor.contract(*mps,gammaA,*mps);
      condensateval = contractor.contract(*mps,condensate,*mps);
      momentumval = contractor.contract(*mps,momentum,*mps);
      chargeval = contractor.contract(*mps,charge,*mps);
      cout << "<GS|H|GS> =           " << E0 << endl <<
	      "<GS|Penalty|GS> =     " << penaltyval << endl << 
	      "<GS|Gamma5|GS> =      " << gamma5val << endl << 
	      "<GS|GammaA|GS> =      " << gammaAval << endl <<
	      "<GS|Condensate|GS> =  " << condensateval << endl <<
	      "<GS|P|GS> =           " << momentumval << endl <<
	      "<GS|Charge|GS> =      " << chargeval << endl;
      //Save observables to file
      if(!file_exists(FILENAME_TRUNC))
      {
	file.open(FILENAME_TRUNC,ios_base::app);
	file <<  "#Level\tN\tD\tdlink\tx\tmg\talpha\tlambda\tE0\tCondensate\tGamma5\tGammaA\tMomentum\tCharge\tPenalty" << endl;
	file.close();
      }
      file.open(FILENAME_TRUNC,ios_base::app);
      file  << "0\t" 
	    << N << "\t" 
	    << D << "\t" 
	    << dlink << "\t" 
	    << x << "\t" 
	    << mg << "\t" 
	    << alpha << "\t" 
	    << lambda <<"\t" << setprecision(15)
	    << E0 << "\t" 
	    << condensateval << "\t"
	    << gamma5val << "\t"
	    << gammaAval << "\t"
	    << momentumval << "\t"
	    << chargeval << "\t"
	    << penaltyval << endl;
      file.close();
      
      //Proceed with the first excited state
      cout << endl <<"--> Computing first excited state for truncated Schwinger Hamiltonian" << endl;
      mps = new MPS(N,D,2);
      KS.constructInitialMPS(*mps);
      mps->setRandomState();
      contractor.findNextExcitedState(hamil,D,computed_states,&E1,*mps);
      outputname=getFilename("truncated",1,dlink,N,x,mg,D);
      //mps->exportMPS("1EX_truncated.dat");
      mps->exportMPS(outputname.c_str());      
      computed_states.push_back(mps);
      //Compute some observables
      penaltyval = contractor.contract(*mps,penalty,*mps);
      gamma5val = contractor.contract(*mps,gamma5,*mps);
      gammaAval = contractor.contract(*mps,gammaA,*mps);
      condensateval = contractor.contract(*mps,condensate,*mps);
      momentumval = contractor.contract(*mps,momentum,*mps);
      chargeval = contractor.contract(*mps,charge,*mps);
      cout << "<1EX|H|1EX> =           " << E1 << endl <<
	      "<1EX|Penalty|1EX> =     " << penaltyval << endl << 
	      "<1EX|Gamma5|1EX> =      " << gamma5val << endl << 
	      "<1EX|GammaA|1EX> =      " << gammaAval << endl <<
	      "<1EX|Condensate|1EX> =  " << condensateval << endl <<
	      "<1EX|P|1Ex> =           " << momentumval << endl <<
	      "<1EX|Charge|1EX> =      " << chargeval << endl;
      //Save observables to file
      file.open(FILENAME_TRUNC,ios_base::app);
      file  << "1\t" 
	    << N << "\t" 
	    << D << "\t" 
	    << dlink << "\t" 
	    << x << "\t" 
	    << mg << "\t" 
	    << alpha << "\t" 
	    << lambda <<"\t" << setprecision(15)
	    << E1 << "\t" 
	    << condensateval << "\t"
	    << gamma5val << "\t"
	    << gammaAval << "\t"
	    << momentumval << "\t"
	    << chargeval << "\t"
	    << penaltyval << endl;
      file.close();
    }
    else
    {
      cout << "--> Computing ground state for truncated Schwinger Hamiltonian via imaginary time evolution" << endl;
      mps = new MPS(N,D,2);
      KS.constructInitialMPS(*mps);
      mps->increaseBondDimension(D);
      
      MPO exp_even(1),exp_odd(1),exp_odd_half(1);
      KS.constructUoddMPO(exp_odd,dt,true);
      KS.constructUoddMPO(exp_odd_half,dt/2.0,true);
      KS.constructUevenMPO(exp_even,dt,true);
      
      compute_evolution(exp_even,exp_odd,exp_odd_half,contractor,*mps,steps,true);
      
      outputname=getFilename("truncated",0,dlink,N,x,mg,D);
      //mps->exportMPS("GS_truncated.dat");
      mps->exportMPS(outputname.c_str());
      computed_states.push_back(mps);
      //Compute some observables
      E0 = real(contractor.contract(*mps,hamil,*mps));
      penaltyval = contractor.contract(*mps,penalty,*mps);
      gamma5val = contractor.contract(*mps,gamma5,*mps);
      gammaAval = contractor.contract(*mps,gammaA,*mps);
      condensateval = contractor.contract(*mps,condensate,*mps);
      momentumval = contractor.contract(*mps,momentum,*mps);
      chargeval = contractor.contract(*mps,charge,*mps);
      cout << "<GS|H|GS> =           " << E0 << endl <<
	      "<GS|Penalty|GS> =     " << penaltyval << endl << 
	      "<GS|Gamma5|GS> =      " << gamma5val << endl << 
	      "<GS|GammaA|GS> =      " << gammaAval << endl <<
	      "<GS|Condensate|GS> =  " << condensateval << endl <<
	      "<GS|P|GS> =           " << momentumval << endl <<
	      "<GS|Charge|GS> =      " << chargeval << endl;
      //Save observables to file
      if(!file_exists(FILENAME_TRUNC))
      {
	file.open(FILENAME_TRUNC,ios_base::app);
	file <<  "#Level\tN\tD\tdlink\tx\tmg\talpha\tlambda\tE0\tCondensate\tGamma5\tGammaA\tMomentum\tCharge\tPenalty" << endl;
	file.close();
      }
      file.open(FILENAME_TRUNC,ios_base::app);
      file  << "0\t" 
	    << N << "\t" 
	    << D << "\t" 
	    << dlink << "\t" 
	    << x << "\t" 
	    << mg << "\t" 
	    << alpha << "\t" 
	    << lambda <<"\t" << setprecision(15)
	    << E0 << "\t" 
	    << condensateval << "\t"
	    << gamma5val << "\t"
	    << gammaAval << "\t"
	    << momentumval << "\t"
	    << chargeval << "\t"
	    << penaltyval << endl;
      file.close();
    }
  }
  else
  {  
    //Exact Schwinger Hamiltonian
    SchwingerHamiltonianSz  HS(N,mu,x,0.0,alpha,lambda);
    const MPO& hamil_schw=HS.getHMPO();

    HS.constructGamma5MPO(gamma5);
    HS.constructGammaAlphaMPO(gammaA);
    HS.constructCondensateMPO(condensate);
    HS.constructMomentumMPO(momentum);
      
    //Now compute the ground state
    cout << "---------------------------------------------------------------------"<<endl;
    cout << "-----------------------Full Schwinger Penalty------------------------"<<endl;
    cout << "---------------------------------------------------------------------"<<endl;
    cout << "--> Computing ground state for full Schwinger Hamiltonian" << endl;
    
    mps=new MPS(N,D,2);
    if(startstate.compare("sc") == 0)
      getSchwingerInitState(*mps,N,is_updown);
    else if(startstate.compare("bg") == 0)
      getSchwingerInitState(*mps,N,is_bgfield);
    else if(startstate.compare("superpos") == 0)
      getSchwingerInitState(*mps,N,is_bgfield_superpos);
    else if(startstate.compare("superpos_noise") == 0)
    {
      getSchwingerInitState(*mps,N,is_bgfield_superpos);
      mps->increaseBondDimensionWithNoise(D,noiselevel);
    }
    else if(startstate.find(".dat") != string::npos && file_exists(startstate))
    {
      cout << "Reading state from file " << startstate << endl;
      mps->importMPS(startstate.c_str());
    }
    else      
    {
      cout << "Using random state" << endl; 
      mps->setRandomState();
    }
    
    contractor.findGroundState(hamil_schw,D,&E0,*mps);
    outputname=getFilename("full",0,dlink,N,x,mg,D);
    //mps->exportMPS("GS_full.dat");
    mps->exportMPS(outputname.c_str());
    computed_states.push_back(mps);
    //Compute some observables  
    gamma5val = contractor.contract(*mps,gamma5,*mps);
    gammaAval = contractor.contract(*mps,gammaA,*mps);
    condensateval = contractor.contract(*mps,condensate,*mps);
    momentumval = contractor.contract(*mps,momentum,*mps);
    cout << "<GS|H|GS> =           " <<E0-0.5*N*mu << endl <<
	    "<GS|Gamma5|GS> =      " << gamma5val  << endl << 
	    "<GS|GammaA|GS> =      " << gammaAval << endl <<
	    "<GS|Condensate|GS> =  " << condensateval << endl <<
	    "<GS|P|GS> =           " << momentumval << endl;
    //Save observables to file
    if(!file_exists(FILENAME))
    {
      file.open(FILENAME,ios_base::app);
      file <<  "#Level\tN\tD\tx\tmg\talpha\tlambda\tE0\tCondensate\tGamma5\tGammaA\tMomentum" << endl;
      file.close();
    }
    file.open(FILENAME,ios_base::app);
    file  << "0\t"
	  << N << "\t" 
	  << D << "\t"
	  << x << "\t" 
	  << mg << "\t" 
	  << alpha << "\t" 
	  << lambda <<"\t" << setprecision(15)
	  << E0-0.5*N*mu  << "\t" 
	  << condensateval << "\t"
	  << gamma5val << "\t"
	  << gammaAval << "\t"
	  << momentumval << endl;
    file.close();
    
    //Proceed with the first excited state
    cout << endl << "--> Computing first excited state for full Schwinger Hamiltonian" << endl;
    mps=new MPS(N,D,2);
    mps->setRandomState();
    contractor.findNextExcitedState(hamil_schw,D,computed_states,&E1,*mps);
    outputname=getFilename("full",1,dlink,N,x,mg,D);
    //mps->exportMPS("1EX_full.dat");
    mps->exportMPS(outputname.c_str());
    computed_states.push_back(mps);
    //Compute some observables  
    gamma5val = contractor.contract(*mps,gamma5,*mps);
    gammaAval = contractor.contract(*mps,gammaA,*mps);
    condensateval = contractor.contract(*mps,condensate,*mps);
    momentumval = contractor.contract(*mps,momentum,*mps);
    cout << "<1EX|H|1EX> =           " <<E1-0.5*N*mu << endl <<
	    "<1EX|Gamma5|1EX> =      " << gamma5val  << endl << 
	    "<1EX|GammaA|1EX> =      " << gammaAval << endl <<
	    "<1EX|Condensate|1EX> =  " << condensateval << endl <<
	    "<1EX|P|1EX> =           " << momentumval << endl;
    //Save observables to file
    file.open(FILENAME,ios_base::app);
    file  << "1\t"
	  << N << "\t" 
	  << D << "\t" 
	  << x << "\t" 
	  << mg << "\t" 
	  << alpha << "\t" 
	  << lambda <<"\t" << setprecision(15)
	  << E1-0.5*N*mu  << "\t" 
	  << condensateval << "\t"
	  << gamma5val << "\t"
	  << gammaAval << "\t"
	  << momentumval << endl;
    file.close();
  }
  
  //Free memory
  for(int i=0; i<computed_states.size(); i++)
  {
    computed_states[i]->clear();
    delete computed_states[i];
  }
  
  cout << "End of program" << endl;
  return 0;
  
}

//Little function to compute the trotter splitting, as it is quite complicated and I want to take two time steps together to improve performance
void compute_evolution(MPO &exp_even,MPO &exp_odd,MPO &exp_odd_half,Contractor &contractor, MPS &mps, int steps,bool normalize)
{
   MPS aux=MPS(mps);

  //Initial half step
  contractor.optimize(exp_odd_half,mps,aux);
  for(int i=1; i<=(steps-1); i++)
  {
    cout.flush();
    mps.setRandomState();
    contractor.optimize(exp_even,aux,mps);
    aux.setRandomState();
    contractor.optimize(exp_odd,mps,aux);
    if(normalize)
      mps.gaugeCond('R',true);
  }

  //Complete last step
  contractor.optimize(exp_even,aux,mps);
  contractor.optimize(exp_odd_half,mps,aux);
  mps=aux;
  if(normalize)
    mps.gaugeCond('R',true);
} 

void getSchwingerInitState(MPS &mps,int N,InitState state)
{
  mwArray spinsiteup(Indices(2,1,1));		spinsiteup.fillWithZero();
  mwArray spinsitedown(Indices(2,1,1));		spinsitedown.fillWithZero();
  
  spinsiteup.setElement(ONE_c,Indices(0,0,0));
  spinsitedown.setElement(ONE_c,Indices(1,0,0));
  
  if(state == is_updown)
  {
    cout << "Constructing state |u>|0>|d>|0>..." << endl;
    
    mps.clear();
    mps=MPS(N,1,2);
    
    for(int i=0; i<N; i++)
    {
      if((i%2)!=0)
	mps.setA(i,spinsitedown);
      else
	mps.setA(i,spinsiteup);
    }
  }  
  else if(state == is_bgfield)
  {
    //Little error check
    if((N%2)!=0)
    {
      cout << "State |d>|-1>|d>|-1>|u>|-1>|d>|-1>|u>...|-1>|u> can only be constructed for an even number of sites" << endl;
      cout << "Program will be aborted..." << endl;
      exit(666);
    }
    cout << "Constructing state |d>|-1>|d>|-1>|u>|-1>|d>|-1>|u>...|-1>|u>" << endl;
    
    mps.clear();
    mps=MPS(N,1,2);
    
    mps.setA(0,spinsitedown);
    mps.setA(N-1,spinsiteup);
    
    for(int i=2; i<N; i++)
    {
      if((i%2)==0)
	mps.setA(i-1,spinsitedown);
      else
	mps.setA(i-1,spinsiteup);
    }
  }
  else if(state == is_bgfield_superpos)
  {
    //Little error check
    if((N%2)!=0)
    {
      cout << "State can only be constructed for an even number of sites" << endl;
      cout << "Program will be aborted..." << endl;
      exit(666);
    }
    cout << "Constructing superposition state" << endl;
    
    mps.clear();
    mps=MPS(N,2,2);
    
    mwArray First,Odd,Even,Last;    
    First=mwArray(Indices(2,1,2));	First.fillWithZero();
    Odd=mwArray(Indices(2,2,2));	Odd.fillWithZero();
    Even=mwArray(Indices(2,2,2));	Even.fillWithZero();
    Last=mwArray(Indices(2,2,1));	Last.fillWithZero();
    
    //Set elements
    First.setElement(ONE_c,Indices(0,0,0));
    First.setElement(ONE_c,Indices(1,0,1));
    
    Odd.setElement(ONE_c,Indices(0,0,0));
    Odd.setElement(ONE_c,Indices(0,1,1));
    
    Even.setElement(ONE_c,Indices(1,0,0));
    Even.setElement(ONE_c,Indices(1,1,1));
    
    Last.setElement(ONE_c,Indices(0,1,0));
    Last.setElement(ONE_c,Indices(1,0,0));
    
    mps.setA(0,First);    
    mps.setA(N-1,Last);       
    
    for(int i=2; i<N; i++)
    {
      if((i%2)==0)      
	mps.setA(i-1,Even);
      else	
	mps.setA(i-1,Odd);      
    }
  }
  else
    cout << "State not supported" << endl;
}

string getFilename(string model, int level, int dlink, int N, double x, double mg, int D)
{
  stringstream name; 
  name.str("");
  if(level==0)
    name << "GS_";
  else
    name << level << "EX_";
  
  if(model.compare("full") == 0)
    name << "full_";
  else
    name << "truncated_d" << dlink << "_";
  
  name << "D" << D <<"_x" << (int) (x*100.0) << "_N" << N << "_mg" << (int) (mg*10000.0) << ".dat";
        
  return name.str();
}
