/**
  \file schwinger_finite_rep_contribution.cpp
  Program to compute the phase diagram for the (truncated) Schwinger Model with background field
   
  \param <file> (string) Name of the driver file containing the parameters
  
  \author Stefan Kühn
  \date 04/19/2017

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
#include "FoldedOperator.h"
#include "Contractor.h"
#include "Indices.h"
#include "misc.h"

#include "Properties.h"
#include "SchwingerHamiltonian.h"
#include "SchwingerHamiltonianSz.h"
#include "KSHamiltonian.h"

#define NUM_OF_PARAMS 1

using namespace std;
using namespace shrt;

//Generate the name of the corresponding input file
string getFilename(string model, int level, int dlink, int N, double x, double mg, int D);
void getProjectorMPOFull(MPO &mpo, int N, int lcut);
void getProjectorMPOTruncated(MPO &mpo, int N, int lcut, int dlink);

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
  double mu=2.0*mg*sqrt(x);
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
  tol=tol==-1? 1.0E-8 : tol;
  noiselevel=noiselevel==-1? 1.0E-1 : noiselevel;
  //File for the results
  ofstream file;
  
  //For the full model, I might want to truncate lmax
  int lmax=N/2;
  if(argc >= NUM_OF_PARAMS + 2)
  {
    if(atoi(argv[2])<=0)
      lmax = N/2;
    else
      lmax = atoi(argv[2]);    
  }  
 
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;  
  //In case PRIMME-support is possible, set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  
 
  //Import the ground state
  MPS mps;
  string inputname;
  if(model.compare("truncated") == 0)
    inputname=getFilename("truncated",0,dlink,N,x,mg,D);
  else
    inputname=getFilename("full",0,dlink,N,x,mg,D);
  
  if(file_exists(inputname.c_str()))
    mps.importMPS(inputname.c_str());
  else
  {
    cout << "Input file " << inputname << " does not exist, program will be aborted..." << endl;
    exit(666);
  }
  
  //MPO for the projector
  MPO projector(1);
  vector<double> qvals;
  
  if(model.compare("truncated") == 0)
  { 
    lmax = (int)round((dlink-1.0)/2.0);
    file.open("contributions_sectors.txt");
    file << "#lcut\t<psi|P_{|l|<=lcut}|psi>\tNorm difference" << endl;
    double overlap,prev=0;
    for(int i=1; i<=lmax; i++)
    {
      getProjectorMPOTruncated(projector,N,i,dlink);
      overlap=real(contractor.contract(mps,projector,mps));
      file << setprecision(15) << i-1 << "\t" << overlap << "\t" <<sqrt(overlap)-prev << endl;
      prev=sqrt(overlap);
    }
    file << setprecision(15) << lmax << "\t" << 1.0 << "\t" <<1.0-prev << endl;
    file.close();
  }
  else
  { 
    if(lmax!=N/2)
      cout << "Checking up to lmax=" << lmax << endl;
    file.open("contributions_sectors.txt");
    file << "#lcut\t<psi|P_{|l|<=lcut}|psi>\tNorm difference" << endl;
    double overlap,prev=0;    
    for(int i=1; i<=lmax; i++)
    {
      getProjectorMPOFull(projector,N,i);
      overlap=real(contractor.contract(mps,projector,mps));
      file << setprecision(15) << i-1 << "\t" << overlap << "\t" <<sqrt(overlap)-prev << endl;
      prev=sqrt(overlap);
    }
    file << setprecision(15) << lmax << "\t" << 1.0 << "\t" <<1.0-prev << endl;
    file.close();
  }  
 
  cout << "End of program" << endl;
  return 0;
  
}

void getProjectorMPOFull(MPO &mpo, int N, int lcut)
{
  mpo.initLength(N);
  if(lcut<=0 || lcut>N/2)
    lcut=N/2;
  int D=2*lcut+1;
  
  mwArray Modd,Meven,Mfirst,Mlast;
  
  Modd = mwArray(Indices(2,D,2,D));	Modd.fillWithZero();
  Meven = mwArray(Indices(2,D,2,D));	Meven.fillWithZero();
  Mfirst = mwArray(Indices(2,1,2,D));	Mfirst.fillWithZero();
  Mlast = mwArray(Indices(2,D,2,1));	Mlast.fillWithZero();
  
  for(int k=0; k<D; k++)
  {
    //Spin up
    Modd.setElement(ONE_c,Indices(0,k,0,k));
    if(k+1<D)
      Meven.setElement(ONE_c,Indices(0,k,0,k+1));
    //Spin down
    if(k-1>=0)
      Modd.setElement(ONE_c,Indices(1,k,1,k-1));
    Meven.setElement(ONE_c,Indices(1,k,1,k));
  }
  //First site, always odd
  Mfirst.setElement(ONE_c,Indices(0,0,0,lcut));  
  if(lcut-1>=0)
    Mfirst.setElement(ONE_c,Indices(1,0,1,lcut-1));
  //Last site
  if((N%2)==0)
  {
    //Even number of sites, last one is even
    for(int k=0; k<D; k++)
    {
      if(k+1<D)
	Mlast.setElement(ONE_c,Indices(0,k,0,0));
      Mlast.setElement(ONE_c,Indices(1,k,1,0));
    }
  }
  else
  {
    //Odd number of sites, last one is odd
    for(int k=0; k<D; k++)
    {
      Mlast.setElement(ONE_c,Indices(0,k,0,0));
      if(k-1>=0)
	Mlast.setElement(ONE_c,Indices(1,k,1,0));
    }
  }
  mpo.setOp(0, new Operator(Mfirst),true);
  mpo.setOp(N-1,new Operator(Mlast),true);
  for(int i=1; i<N-1; i++)
  {
    if(i==1)
      mpo.setOp(i,new Operator(Meven),true);
    else if(i==2)
      mpo.setOp(i,new Operator(Modd),true);
    else if((i%2)!=0)
      mpo.setOp(i,&mpo.getOp(1),false);
    else
      mpo.setOp(i,&mpo.getOp(2),false);      
  }
}
void getProjectorMPOTruncated(MPO &mpo, int N, int lcut, int dlink)
{
  mpo.initLength(2*N-1);
  mwArray Ploc=identityMatrix(dlink);
  mwArray id=identityMatrix(2);
  int lmax = (int)round((dlink-1.0)/2.0);
  
  for(int i=0; i<=(lmax-lcut); i++)
  {
    Ploc.setElement(ZERO_c,Indices(i,i));
    Ploc.setElement(ZERO_c,Indices(2*lmax-i,2*lmax-i));
  }
  
  Ploc.reshape(Indices(dlink,1,dlink,1));
  id.reshape(Indices(2,1,2,1));
  
  mpo.setOp(0,new Operator(id),true);
  mpo.setOp(1,new Operator(Ploc),true);
  
  for(int i=2; i<2*N-1; i+=2)
  {
    mpo.setOp(i,&mpo.getOp(0),false);
    if((i+1)<2*N-2)
      mpo.setOp(i+1,&mpo.getOp(1),false);
  }
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