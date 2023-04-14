
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"

#include "IsingHamiltonian.h"
#include "MPOMultiplierHermitian.h"
#include "SpinMPO.h"
#include <cmath>
#include <unistd.h>

#include "misc.h"
//#include "quicksort.h"

using namespace std;
using namespace shrt;

#define MAXLEN 120

int d=2;
double tol=1E-8;

/** Read a GS MPS from a file (needed input)
    and compute the spectrum of the effective local Hamiltonian for a certain number of sites over the whole chain.
*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  string directory="";
  // if(argc>2){
  //   const char* indir=argv[++cntr];
  //   directory=string(indir)+"/";
  // }
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  const string mpsfile=props.getProperty("mpsfile"); // full path expected
  // Hamiltonian parameters
  double J=props.getDoubleProperty("J");
  double g=props.getDoubleProperty("g");
  double h=props.getDoubleProperty("h");  
  int block = props.getIntProperty("block");
  if(block<=0) block=2;
  const string outfile=props.getProperty("outfile");  
  bool app=(props.getIntProperty("append")!=0); // default: append
  int nrEV=props.getIntProperty("numberEV"); // max nr of eigenvalues to record
  

  // Read the MPS
  if(!file_exists(mpsfile)){
    cout<<"ERROR: File for MPS "<<mpsfile<<" not found!"<<endl;
    exit(1);
  }

  MPS state(1,1,1);
  state.importMPS(mpsfile.data());
  cout<<"Read GS file "<<endl;
  int L=state.getLength();

  state.gaugeCond('R',1);
  state.gaugeCond('L',1);

  // Prepare output file
  ofstream* out;
  if(!app||!file_exists(outfile))
    out=new ofstream(outfile.data());
  else{
    out=new ofstream(outfile.data(),ios::app);
    // Write the header
  }
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output (append="<<app<<")"<<endl;
    exit(1);
  }
  out->close();delete out;

  IsingHamiltonian hamIsing(L,d,J,g,h);
  const MPO& mpoH=hamIsing.getHMPO();
  
  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);

  cout<<"The energy of the read state is "<<contractor.contract(state,mpoH,state)<<endl;
  
  for(int k=0;k<L-block;k++){

    MPOMultiplierHermitian Heff(block+2); // maybe could use non-Hermitian, as it is in fact Hermitian (factor 2 cheaper)
    contractor.getEffectiveOperatorMPOMultiplier(state,mpoH,block,k,Heff);
    int totDim=Heff.getSize();
    cout<<"Hamiltonian for site "<<k<<" range "<<block<<" dimensions "<<totDim<<endl;
    // Now extract the exact spectrum (no eigenvectors)
    vector<complex_t> Dval;
    mwArray U;
    //wrapper::eig(Heff,Dval,U,false);
    //    wrapper::eigs(Heff,nrEV,"SR",Dval,U,false);
    wrapper::eigs_primme(Heff,min(nrEV,totDim),primme_smallest,Dval,U,false);
    out=new ofstream(outfile.data(),ios::app);
    *out<<setprecision(15);
    *out<<k<<"\t";
    for(int p=0;p<min(nrEV,(int)Dval.size());p++){
      *out<<real(Dval[p])<<"\t";
    }
    if(Dval.size()<nrEV) // fill in with 0
      for(int p=Dval.size();p<nrEV;p++) *out<<0<<"\t";
    *out<<endl;
    out->close();delete out;
  }

}
