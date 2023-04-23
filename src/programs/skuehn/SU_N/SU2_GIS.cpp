/**
  \file SU2_GIS.cpp
  Driver program to compute the ground state of the SU(2) Hamiltonian  on the gauge invariant subspace as developed by Tao and Pablo
   
  \param <N> (int) Number of Fermions
  \param <D> (int) Bond dimension
  \param <epsilon> (double) Hopping parameter 
  \param <mu> (double)  Mass for the Hamiltonian
  \param <g> (double) Coupling constant for the Hamiltonian
  
  \author Stefan Kühn
  \date 02/23/2018

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
#define FILENAME "SU2GIS.txt"
#define FILENAMECORR "SU2QQCorr.txt"

using namespace std;
using namespace shrt;

int main(int argc,const char* argv[])
{  
  //Code Revision (magic field, automatically filled in by subversion)
  string revision="$Rev: 1593 $";
  //Parse revision string
  revision=revision.substr(6,4);
  
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1)
  {
    cout << "Program compiled on " << __DATE__ << ", " << __TIME__ << " with code revision " << revision << endl;
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage:   " << argv[0] << " <param_file>" << endl;
    return -1;
  }
  
  //Read the properties file
  const char* infile=argv[1];
  Properties props(infile);
   
  //Number of spin sites
  int N=props.getIntProperty("N");
  //Bond dimension
  int D=props.getIntProperty("D");
  //Hopping parameter
  double epsilon=props.getDoubleProperty("t");
  //Mass
  double mu=props.getDoubleProperty("m");
  //Coupling constant
  double g=props.getDoubleProperty("g");
  //Penalty strength
  double lambda=props.getDoubleProperty("lambda");
  //Penalty for single occupancy
  double eta=props.getDoubleProperty("eta");
  //Positions for the external charges
  int k=props.getIntProperty("pc1");
  int l=props.getIntProperty("pc2");
  //File of the inital guess
  string ig_file=props.getProperty("ig_file");
  
  //File for results
  ofstream file;
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  contractor.setConvTol(EPS);

  //Get the Hamiltonian according to the case if external charges are involved or not
  SU2HamiltonianGIS *HSU2;
  MPS mps;
  if(k==-1 && l==-1)
  {
    cout << "Constructing SU(2) Hamiltonian without external charges" << endl;
    HSU2 = new SU2HamiltonianGIS(N,epsilon,mu,g,lambda);
    mps = MPS(N,D,4);
  }
  else
  {
    cout << "Constructing SU(2) Hamiltonian with external charges" << endl;
    HSU2 = new SU2HamiltonianGIS(N,epsilon,mu,g,k,l,lambda,eta);
    mps = MPS(N+2,D,4);
  }
  const MPO& hamil=HSU2->getHMPO();    
    
  //Now compute the groundstate starting from a random state or a given input state
  if(!ig_file.empty() && file_exists(ig_file))
  {
    cout << "Importing MPS from " << ig_file << endl;
    mps.importMPS(ig_file.c_str());
  }
  else
    mps.setRandomState();
  
  double E0=0.0;
    
  cout << "------------------------------------------------------------------------"<<endl;
  cout << "------------------Starting to compute ground state----------------------"<<endl;
  cout << "------------------------------------------------------------------------"<<endl; 
  cout.flush();
  
  contractor.findGroundState(hamil,D,&E0,mps);
   
  //Save the ground state
  
  string GS_name;
  if(k==-1 && l==-1)
    GS_name="SU2GS_";
  else
    GS_name="SU2GSExt_";
  //String needed for the conversion
  ostringstream convert;  
  //Build the filename
  convert << N;
  GS_name = GS_name+"N"+convert.str();
  convert.str("");
  convert.clear();
  convert << D;
  GS_name = GS_name+"_D"+convert.str();
  convert.str("");
  convert.clear();
  convert << epsilon;
  GS_name = GS_name+"_t"+convert.str();
  convert.str("");
  convert.clear();
  convert << mu;
  GS_name = GS_name+"_m"+convert.str();
  convert.str("");
  convert.clear();
  convert << g;
  GS_name = GS_name+"_g"+convert.str();
  convert.str("");
  convert.clear();
  convert << lambda;
  GS_name = GS_name+"_la"+convert.str();
  if(k!=-1 && l!=-1)
  {
    if(k<l)
    {
      convert.str("");
      convert.clear();
      convert << k;
      GS_name = GS_name+"_cleft"+convert.str();
      convert.str("");
      convert.clear();
      convert << l;
      GS_name = GS_name+"_cright"+convert.str();
    }
    else
    {
      convert.str("");
      convert.clear();
      convert << l;
      GS_name = GS_name+"_cleft"+convert.str();
      convert.str("");
      convert.clear();
      convert << k;
      GS_name = GS_name+"_cright"+convert.str();
    }
    convert.str("");
    convert.clear();
    convert << eta;
    GS_name = GS_name+"_eta"+convert.str();
  }
  GS_name = GS_name+".dat";  
  mps.exportMPS(GS_name.c_str());
  
  if(file_exists(FILENAME))
  {
    file.open(FILENAME, ios_base::app);
    if(k==-1 && l==-1)    
      file << setprecision(16) << N << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << E0 << endl;
    else
      file << setprecision(16) << N << "\t" << D << "\t"<< k << "\t" << l << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << E0 << endl;
  }
  else
  {
    file.open(FILENAME, ios_base::app);
    if(k==-1 && l==-1)    
    {
      file << "#N\tD\tt\tm\tg\tE0" << endl;
      file << setprecision(16) << N << "\t" << D << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << E0 << endl;
    }
    else
    {
      file << "#N\tD\tcleft\tcright\tt\tm\tg\tE0" << endl;
      file << setprecision(16) << N << "\t" << D << "\t"<< k << "\t" << l << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << E0 << endl;
    }
  }
  file.close();
  
  if(k!=-1 && l!=-1)
  {
    //Compute the correlation between the external charges if there are any
    MPO Qcorr(1);
    complex_t cQx,cQy,cQz;
    
    HSU2->getExtChargeCorrelationMPO(Qcorr,x_comp);
    cQx = contractor.contract(mps,Qcorr,mps);
    
    HSU2->getExtChargeCorrelationMPO(Qcorr,y_comp);
    cQy = contractor.contract(mps,Qcorr,mps);
    
    HSU2->getExtChargeCorrelationMPO(Qcorr,z_comp);
    cQz = contractor.contract(mps,Qcorr,mps);
    
    if(file_exists(FILENAMECORR))
    {
      file.open(FILENAMECORR, ios_base::app);
      file << setprecision(16) << N << "\t" << D << "\t"<< k << "\t" << l << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << cQx << "\t" << cQy << "\t" << cQz << endl;
    }
    else
    {
      file.open(FILENAMECORR, ios_base::app);
      file << "#N\tD\tcleft\tcright\tt\tm\tg\tQ^xQ^x\tQ^yQ^y\tQ^zQ^z" << endl;
      file << setprecision(16) << N << "\t" << D << "\t"<< k << "\t" << l << "\t" << epsilon << "\t" << mu << "\t" << g << "\t" << cQx  << "\t" << cQy << "\t" << cQz  << endl;
    }
    file.close();
    
  }
    
    
    
  delete HSU2;
  cout << "End of program" << endl;
  return 0; 
  
}


    
    
