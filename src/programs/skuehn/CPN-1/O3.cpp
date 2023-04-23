/**
   \file O3Hamiltonian.cpp
   Truncated O(3) rotor Hamiltonian
   
  \author Stefan Kühn
  \date 03/30/2019

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
#include "Properties.h"
#include "misc.h"

#include "O3Hamiltonian.h"

#define NUM_OF_PARAMS 1
#define FILENAME "O3.txt"
#define FILENAMECONFIG "O3config.txt"


using namespace std;
using namespace shrt;


int main(int argc,const char* argv[])
{  
  //Code Revision (magic field, automatically filled in by subversion)
  string revision="$Rev: 1522 $";
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
  
  //Read the properties file
  const char* infile=argv[1];
  Properties props(infile);
   
  //Number of spin sites
  int N=props.getIntProperty("N");
  //Bond dimension
  int D=props.getIntProperty("D");
  //Truncation parameters
  int lmax=props.getIntProperty("lmax");
  int mmax=props.getIntProperty("mmax");
  //Parameters for the Hamiltonian
  double g2=props.getDoubleProperty("g2");
  double au=props.getDoubleProperty("au");
  double offset=props.getDoubleProperty("offset");
  double lambda=props.getDoubleProperty("lambda");
  double eta=props.getDoubleProperty("eta");
  int q = props.getIntProperty("q");
  if(lambda < 0.)
    lambda = 0.;
  string ig_file=props.getProperty("ig_file");  
  //Determine the constants for the Hamiltonian
  double c1=g2/2.;
  double c2=au;
  double c3=-1./g2;
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
     
  //Get the Hamiltonian
  O3Hamiltonian *HO3;
  if(eta<0)
    HO3 = new O3Hamiltonian(N,lmax,mmax,c1,c2,c3,lambda);
  else
    HO3 = new O3Hamiltonian(N,lmax,mmax,c1,c2,c3,lambda,eta,q);    
  const MPO& hamil=HO3->getHMPO();
  cout << "Successfully retrieved O(3) Hamiltonian" << endl;
  cout.flush(); 
  
  //MPS for the ground state
  int d = hamil.getOp(0).getd();
  MPS mps(N,D,d),*first_ex;
  first_ex = new MPS(N,D,d);
  first_ex->setRandomState();
  
  //Now compute the groundstate starting from a random state or a given input state
  if(!ig_file.empty() && file_exists(ig_file))
  {
    cout << "Importing MPS from " << ig_file << endl;
    mps.importMPS(ig_file.c_str());
  }
  else
    mps.setRandomState();
  
  cout << "---------------------------------------------------------------------"<<endl;
  cout << "---------------Starting to compute ground state----------------------"<<endl;
  cout << "---------------------------------------------------------------------"<<endl; 
  cout.flush();
  double E0,E1;
  vector<MPS*> excitedstates;
  contractor.findGroundState(hamil,D,&E0,mps,offset);
  cout << "Ground state energy: " << E0 << endl; 
  excitedstates.push_back(&mps);
  contractor.findNextExcitedState(hamil,D,excitedstates,&E1,*first_ex,offset);
  cout << "First excited state energy: " << E1 << endl;
  string GS_name,EX_name,name;
  GS_name="O3GS_";
  EX_name="O31EX_";
  //String needed for the conversion
  ostringstream convert;  
  //Build the filename
  convert << N;
  name = name+"N"+convert.str();
  convert.str("");
  convert.clear();
  convert << D;
  name = name+"_D"+convert.str();
  convert.str("");
  convert.clear();
  convert << lmax;
  name = name+"_lmax"+convert.str();
  convert.str("");
  convert.clear();
  convert << mmax;
  name = name+"_mmax"+convert.str();
  convert.str("");
  convert.clear();
  convert << au;
  name = name+"_au"+convert.str();
  convert.str("");
  convert.clear();
  convert << g2;
  name = name+"_g2"+convert.str()+".dat";
  GS_name = GS_name + name;
  EX_name = EX_name + name;
  mps.exportMPS(GS_name.c_str());
  first_ex->exportMPS(EX_name.c_str());
  
  //Compute the entropy in the ground state and the first excited state
  double S_GS,S_EX;  
  S_GS = contractor.getEntropy(mps); 
  S_EX = contractor.getEntropy(*first_ex);
  
  //Compute the potential in the ground state and the first excited state
  complex_t pot_GS,pot_EX;
  MPO pot(1);
  HO3->getPotentialMPO(pot);
  pot_GS = contractor.contract(mps,pot,mps);
  pot_EX = contractor.contract(*first_ex,pot,*first_ex);
  
  ofstream file;
  
  //Compute site resolved expectation values of Jz and J^2
  complex_t Jsq_total_GS=ZERO_c, Jz_total_GS=ZERO_c, Jsq_total_EX=ZERO_c, Jz_total_EX=ZERO_c;
  vector<complex_t> Jsq_vals_GS,Jz_vals_GS,Jsq_vals_EX,Jz_vals_EX;
  Jsq_vals_GS.resize(N);
  Jz_vals_GS.resize(N);
  Jsq_vals_EX.resize(N);
  Jz_vals_EX.resize(N);
  MPO Op(1);
  for (int i=0;i<N; i++)
  {  
    HO3->getJsqMPO(Op,i);
    Jsq_vals_GS[i] = contractor.contract(mps,Op,mps);
    Jsq_vals_EX[i] = contractor.contract(*first_ex,Op,*first_ex);
    Jsq_total_GS += Jsq_vals_GS[i];
    Jsq_total_EX += Jsq_vals_EX[i];
    HO3->getJzMPO(Op,i);
    Jz_vals_GS[i] = contractor.contract(mps,Op,mps);
    Jz_vals_EX[i] = contractor.contract(*first_ex,Op,*first_ex);
    Jz_total_GS += Jz_vals_GS[i];
    Jz_total_EX += Jz_vals_EX[i];
  } 
  
  if(file_exists(FILENAMECONFIG))
  {
    file.open(FILENAMECONFIG,ios_base::app);
    for (int i=0; i<N; i++)
      file << setprecision(15) << Jsq_vals_GS[i] << "\t" << Jz_vals_GS[i] << "\t" << Jsq_vals_EX[i] << "\t" << Jz_vals_EX[i] << endl;
  }
  else
  {
    file.open(FILENAMECONFIG,ios_base::app);
    file <<"#<J^2>_GS\t<J_z>_GS\t<J^2>_EX\t<J_z>_EX\n";
    for (int i=0; i<N; i++)
      file << setprecision(15) << Jsq_vals_GS[i] << "\t" << Jz_vals_GS[i] << "\t" << Jsq_vals_EX[i] << "\t" << Jz_vals_EX[i] << endl;
  }
  file.close();
  
  //Compute the susceptibility
  MPO charge_square(1),charge(1);
  complex_t Qsq_GS,Qsq_EX;
  HO3->getChargeMPO(charge);  
  HO3->getChargeSquareMPO(charge_square); 
  
  Qsq_GS=contractor.contract(mps,charge_square,mps);
  Qsq_EX=contractor.contract(*first_ex,charge_square,*first_ex);
  
  if(file_exists(FILENAME))
  {
    file.open(FILENAME,ios_base::app);
    file << setprecision(15) << N << "\t" << D << "\t" << lmax << "\t" << mmax << "\t" << c1 << "\t" << c2 << "\t" << c3 << "\t" << E0 << "\t" << E1 << "\t" << pot_GS << "\t" << pot_EX << "\t" <<  S_GS << "\t" << S_EX << "\t" << Jsq_total_GS << "\t" << Jz_total_GS << "\t" << Jsq_total_EX << "\t" << Jz_total_EX << "\t" << Qsq_GS << "\t" << Qsq_EX << endl;
  }
  else
  {
    file.open(FILENAME,ios_base::app);
    file <<"#N\tD\tlmax\tmmax\tc1\tc2\tc3\tE0\tE1\t<GS|nn|GS>\t<EX|nn|EX>\tS_GS\tS_EX\t<GS|J^2|GS>\t<GS|J^z|GS>\t<EX|J^2|EX>\t<EX|J^z|EX>\t<GS|Q^2|GS>\t<EX|Q^2|EX>\n";
    file << setprecision(15) << N << "\t" << D << "\t" << lmax << "\t" << mmax << "\t" << c1 << "\t" << c2 << "\t" << c3 << "\t" << E0 << "\t" << E1 << "\t" << pot_GS << "\t" << pot_EX << "\t" << S_GS << "\t" << S_EX << "\t" << Jsq_total_GS << "\t" << Jz_total_GS << "\t" << Jsq_total_EX << "\t" << Jz_total_EX << "\t" << Qsq_GS << "\t" << Qsq_EX << endl;
  }
  file.close();
  
  //Free memory
  first_ex->clear();
  delete HO3;
  delete first_ex;
  
  cout << "End of program" << endl;
  return 0;
  
}
