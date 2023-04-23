/**
   \file ContractMPS.cpp
   Given two MPS of the same length and physical dimension as input, this programm computes the overlap
   
  \param <mps1> (string) Name of the file containing the first MPS
  \param <mps2> (string) Name of the file containing the second MPS
  \param <acc> (double) Accuracy for the contractor (optional)
  
  \author Stefan Kühn
  \date 17/01/2013

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
#define NUM_OF_PARAMS 2

using namespace std;
using namespace shrt;

int main(int argc,const char* argv[]){
  int cntr=0;
   //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1){
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage: " << argv[0] << " mps1 mps2" << endl;
    cout << "mps1: " << "Name of the file containing the first MPS" << endl
	 << "mps2: " << "Name of the file containing the second MPS" << endl;
    return -1;
  }
  
  //Name for MPS1
  const char* filename_mps1=argv[++cntr];
  string mps1_file(filename_mps1);
  //Name for MPS2
  const char* filename_mps2=argv[++cntr];
  string mps2_file(filename_mps2);
  //Set the accuracy of contractor?
  bool set_accuracy = false;
  //Variable for the desired accuracy in case it should be set
  double desired_accuracy;
  //File for result
  ofstream result;
  
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
	<< "MPS 1: " << mps1_file << endl
	<< "MPS 2: " << mps2_file << endl
	<< "-----------" << endl << endl;
  }
  else {
    cout << "Parameter: " << endl
	<< "-----------" << endl
	<< "MPS 1: " << mps1_file << endl
	<< "MPS 2: " << mps2_file << endl
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
  MPS mps1,mps2;
  mps1.importMPS(mps1_file.c_str());
  mps2.importMPS(mps2_file.c_str());

  complex_t overlap = contractor.contract(mps1,mps2);
  
  cout << "<MPS1|MPS2>=" << scientific << setprecision(11) << overlap << endl;
  cout << "|<MPS1|MPS2>|=" << scientific << setprecision(11) << abs(overlap) << endl;
  
  result.open("res.txt");
  result << scientific << setprecision(10) << overlap << endl;
  result.close();
  
  cout << "End of program" << endl;
  return 0;
}

