/**
   \file Z2_phasediagram.cpp
   
   This program is designed to compute the phase structure of the Z2 Hamiltonian in spin formulation 
   
  \param <inputfile> (string) File containing the parameters for the run
  
  \author Stefan Kühn
  \date 14/12/2017

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
#include "unistd.h"

#include "MPS.h"
#include "MPO.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include "Indices.h"
#include "misc.h"

#include "Properties.h"
#include "Z2Hamiltonian.h"

#define NUM_OF_PARAMS 1
#define FILENAME "Z2pd.txt"


using namespace std;
using namespace shrt;


int main(int argc,const char* argv[])
{
  //Code Revision (magic field, automatically filled in by subversion)
  string revision="$Rev: 1579 $";
  //Parse revision string
  revision=revision.substr(6,4);
  
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1)
  {
    cout << "Program compiled on " << __DATE__ << ", " << __TIME__ << " with code revision " << revision << endl;
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage:   " << argv[0] << " <inputfile>" << endl;
    return -1;
  }
  
  /**************************************************************************************
				 * Set up parameters*
   * ***********************************************************************************/
  
  const char* infile=argv[1];
  Properties props(infile);
  
  //System size
  int L = props.getIntProperty("L");
  //Bond dimension
  int D=props.getIntProperty("D");
  //Hamiltonian parameters
  double x=props.getDoubleProperty("x");
  double mu=props.getDoubleProperty("mu");  
  double alpha=props.getDoubleProperty("alpha"); 
  //Penalty for unphysical states
  double penalty=props.getDoubleProperty("penalty"); 
  //Do I want to append to file
  int app=props.getIntProperty("app");
  
  cout << "------------------------------------------------------------------------"<<endl;
  cout << "Parameters: " << endl
       << " --> L           = " << L << endl
       << " --> D           = " << D << endl
       << " --> x           = " << x << endl
       << " --> mu          = " << mu << endl
       << " --> pen         = " << penalty << endl
       << " --> alpha       = " << alpha << endl
       << " --> app         = " << app << endl;
  cout << "------------------------------------------------------------------------"<<endl;
  
  /**************************************************************************************
			      * Set up MPOs & Contractor*
   * ***********************************************************************************/
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  contractor.setConvTol(1E-6);
  
  //Get the Hamiltonian
  Z2Hamiltonian HZ2(L,x,mu,alpha,penalty);
  const MPO& hamil=HZ2.getHMPO();        
  cout << "Successfully retrieved Z2-Hamiltonian" << endl;
  
  MPO PenaltyMPO(1);
  HZ2.getPenaltyMPO(PenaltyMPO);
  
  //Now prepare a MPS for the ground state
  MPS mps(2*L-1,D,2);  
  mps.setRandomState();
  
  //Variables for observables
  double E0,S;
  complex_t Epen;
  
  /**************************************************************************************
			      * Actual computation*
   * ***********************************************************************************/
  cout << "Computing ground state" << endl;

  contractor.findGroundState (hamil,D,&E0,mps); 
  
  //Save the ground state
  string GS_name="Z2GS_";
  //String needed for the conversion
  ostringstream convert;  
  //Build the filename
  convert << L;
  GS_name = GS_name+"L"+convert.str();
  convert.str("");
  convert.clear();
  convert << D;
  GS_name = GS_name+"D"+convert.str();
  convert.str("");
  convert.clear();
  convert << (int) x*1000.0;
  GS_name = GS_name+"x"+convert.str();
  convert.str("");
  convert.clear();
  convert << (int) mu*1000.0;
  GS_name = GS_name+"mu"+convert.str();
  convert.str("");
  convert.clear();
  convert << penalty;
  GS_name = GS_name+"la"+convert.str()+".dat";
  mps.exportMPS(GS_name.c_str());
  
  //Compute the entropy
  S=contractor.getEntropy(mps);
  
  //Compute the penalty contribution
  Epen=contractor.contract(mps,PenaltyMPO,mps);
  
  //Prepare the file for the results
  ofstream file; 
  if(app==1)
    file.open(FILENAME, ios_base::app);
  else
  {
    file.open(FILENAME);
    file << "#Created with code revision " << revision << endl;
    file << "#L\tx\tmu\tlambda\tD\tE0\tS\tPenalty" << endl;
  }
  file << setprecision(16) ;
  file << L << "\t" << x << "\t" << mu << "\t" << penalty << "\t"  << D << "\t" << E0 << "\t" << S << "\t" << Epen << endl;
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

