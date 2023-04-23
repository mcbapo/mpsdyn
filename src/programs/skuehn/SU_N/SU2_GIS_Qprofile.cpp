/**
  \file SU2_GIS_Qprofile.cpp
  Driver program to compute the charge profile for states of the SU(2) Hamiltonian on the gauge invariant subspace as developed by Tao and Pablo
  
  \author Stefan Kühn
  \date 03/08/2018

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
#include "Contractor.h"
#include "Indices.h"
#include "misc.h"
#include "Properties.h"

#include "SU2HamiltonianGIS.h"


#define EPS 1E-10
#define NUM_OF_PARAMS 1
#define FILENAME "SU2QQvals.txt"

using namespace std;
using namespace shrt;

int main(int argc,const char* argv[])
{  

  //Read the properties file
  const char* infile=argv[1];

  //File for results
  ofstream file;
  
  //Mps
  MPS mps;
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
     
  //Now compute the groundstate starting from a random state or a given input state
  if(file_exists(infile))
  {
    cout << "Importing MPS from " << infile << endl;
    mps.importMPS(infile);
    
    int L = mps.getLength();
    mwArray Qsquare,Qx,Qy,Qz,sigma_z,sigma_plus,sigma_minus,id_fermi;
    
    id_fermi = identityMatrix(2);    
    complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
    sigma_z=mwArray(Indices(2,2),data_z);
    complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
    sigma_minus=mwArray(Indices(2,2),data_sm);  
    complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
    sigma_plus=mwArray(Indices(2,2),data_sp);    
    
    Qx = 0.5*I_c*(kron(sigma_minus,sigma_plus) - kron(sigma_plus,sigma_minus));
    Qy = -0.5*(kron(sigma_minus,sigma_plus) + kron(sigma_plus,sigma_minus));
    Qz = 0.25*(kron(sigma_z,id_fermi) - kron(id_fermi,sigma_z));
    Qsquare = Qz*Qz;
    
    MPS curr_mps;
    
    vector<complex_t> val_x,val_y,val_z;
    complex_t tmp;
    
    for(int i=0; i<L; i++)
    {
      for(int j=i; j<L; j++)
      {
	cout << "Computing <Q_" << i << " Q_" << j << ">" << endl;
	if(i==j)
	{
	  curr_mps = mps;
	  curr_mps.applyLocalOperator(i,Qsquare);
	  tmp=contractor.contract(mps,curr_mps);
	  val_x.push_back(tmp);
	  val_y.push_back(tmp);
	  val_z.push_back(tmp);
	}
	else
	{
	  curr_mps = mps;
	  curr_mps.applyLocalOperator(i,Qx);
	  curr_mps.applyLocalOperator(j,Qx);  
	  val_x.push_back(contractor.contract(mps,curr_mps));
	  
	  curr_mps = mps;
	  curr_mps.applyLocalOperator(i,Qy);
	  curr_mps.applyLocalOperator(j,Qy);
	  val_y.push_back(contractor.contract(mps,curr_mps));
	  
	  curr_mps = mps;
	  curr_mps.applyLocalOperator(i,Qz);
	  curr_mps.applyLocalOperator(j,Qz);
	  val_z.push_back(contractor.contract(mps,curr_mps));
	}
      }    
    }
    
    
    
    //Save the results
    file.open(FILENAME);
    for(int i=0; i<val_x.size(); i++)
    {
      if(i==0)
	file << val_x[i];
      else
	file << "\t" << val_x[i];
    }
    file << endl;
    for(int i=0; i<val_y.size(); i++)
    {
      if(i==0)
	file << val_y[i];
      else
	file << "\t" << val_y[i];
    }
    file << endl;
    for(int i=0; i<val_z.size(); i++)
    {
      if(i==0)
	file << val_z[i];
      else
	file << "\t" << val_z[i];
    }
    file << endl;
    file.close();
  }
  
 
  cout << "End of program" << endl;
  return 0; 
  
}


    
    
