/**
   \file SU2NoLinksVariance.cpp
   
   This program is designed to compute the variance in energy after the large scale production runs were completed
   
  \param <inputfile> (string) File containing the parameters for the run
  
  \author Stefan Kühn
  \date 02/05/2016ss

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
#include "SU2HamiltonianNoLinks.h"
#include "SU2HamiltonianNoLinksQ.h"
#include "SU2HamiltonianNoLinksDimless.h"
#include "SU2HamiltonianNoLinksDimlessQ.h"

#define NUM_OF_PARAMS 1
#define FILENAME "SU2variances.txt"

using namespace std;
using namespace shrt;


/**Given the parameters, construct the name of the state*/
string constructFileName(int N, int D,double x, double mg, int dlink, double tol, int level_nr);

int main(int argc,const char* argv[])
{
  //Code Revision (magic field, automatically filled in by subversion)
  string revision="$Rev: 1524 $";
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
  int N = props.getIntProperty("N");
  //Bond dimension
  int D=props.getIntProperty("D");
  //Hamiltonian parameters
  double x=props.getDoubleProperty("x");
  double mg=props.getDoubleProperty("mg");  
  int dlink=props.getIntProperty("dlink");
  //Penalty for unphysical states
  double penalty=props.getDoubleProperty("penalty");  
  //Penalty for states with non-vanishing flux at the site N
  double szpenalty=props.getDoubleProperty("szpenalty");  
  //Penalty for states with non-vanishing U(1) charge
  double qpenalty=props.getDoubleProperty("qpenalty");
  //Contractor tolerance
  double tol=props.getDoubleProperty("tol");
  
  cout << "------------------------------------------------------------------------"<<endl;
  cout << "Parameters: " << endl
       << " --> N           = " << N << endl
       << " --> D           = " << D << endl
       << " --> x           = " << x << endl
       << " --> mg          = " << mg << endl
       << " --> dlink       = " << dlink << endl
       << " --> pen         = " << penalty << endl
       << " --> szpen       = " << szpenalty << endl
       << " --> qpen        = " << qpenalty << endl
       << " --> tol         = " << tol << endl << endl;
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
  contractor.setConvTol(tol);

  //Get the Hamiltonian (dimensionless version with charge penalty)
  SU2HamiltonianNoLinksDimlessQ HNGI(N,x,mg,dlink,penalty,szpenalty,qpenalty);
  const MPO& hamil=HNGI.getHMPO();        
  cout << "Successfully retrieved SU2Hamiltonian" << endl;
 
  MPS mps;
  
  
  /**************************************************************************************
			      * Actual computation*
   * ***********************************************************************************/
  //Variable for the energy
  double Eval, Esquare_val=0.0;
  //Auxiliary variables for filenames and statuses
  string filename;
  //Level number
  int level_nr=0;
  //Vector for the energies and the variances
  vector<double> variances,energies;
  while(1)
  { 
    filename=constructFileName(N,D,x,mg,dlink,tol,level_nr);
    if(file_exists(filename.c_str()))
    {
      //Import state and compute energy as well as variance
      mps.importMPS(filename.c_str());
      cout << "Found level number " << level_nr << " in " << filename << endl;
      Eval=real(contractor.contract(mps,hamil,mps));
      Esquare_val=real(contractor.contract2(hamil,mps));
      energies.push_back(Eval);
      variances.push_back(Esquare_val-Eval*Eval);
    }
    else
    {
      break;
    }
    level_nr++;
  }
  
  ofstream file;
  file.open(FILENAME);
  file << "#Created with code revision " << revision << endl;
  file << setprecision(16) << "#N\tD\tdlink\tx\tmg\ttol\tpen\tszpen\tqpen" << endl;
  file << "#E0\tVar(E0)" << endl;
  file << "#E1\tVar(E1)" << endl;
  file << "#..." << endl;
  file << setprecision(15) << N << "\t" << D << "\t" << dlink << "\t" << x << "\t" << mg << "\t" << tol << "\t" << penalty << "\t" << szpenalty << "\t" << qpenalty << endl;
      
  for(int i=0; i<energies.size(); i++)
  {    
    file  << energies[i] << "\t"
	  << variances[i] << endl;
  } 
  file.close();  
  
  cout << "End of program" << endl;
  return 0;  
}

string constructFileName(int N, int D,double x, double mg, int dlink, double tol, int level_nr)
{
  stringstream sstm;    
  sstm.str("");
  if(level_nr==0)
    sstm << "GS_";
  else
    sstm<< level_nr << "EX_";
  sstm << "N" << N << "D" << D << "d" << dlink << "x" << (int) x << "mg" << (int) (mg*1000) << "tol" << tol << ".dat";
  return sstm.str();
}
