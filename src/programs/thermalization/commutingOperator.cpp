
//#include "mwArray.h"
//#include "MPO.h"
#include "IsingHamiltonian.h"
//#include "Operator.h"
#include "Contractor.h"
#include "Properties.h"
//#include "JoinedOperator.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <stdio.h>

#define MAXLEN 200  // Max nr of characters in a filename


/** 
    Program commutingOperator: Given a Hamiltonian (Ising is 
    the only one I have so far), find an operator A on L 
    sites that minimizes tr([H,A]^2).

    To compile: make commut

    Receives as argument a Properties file, such as
    config/commProp.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.
    To replace some of the values in the file by a command line
    option, the additional arguments have to be of the form
    -propName=propValue

    The properties that will be required:
    + L (int) length of the operator (support of A)
    + D (int) maximum bond dimension of A'S MPO
    + J (double) parameter $J$ of the Hamiltonian
    + g (double) parameter $g$ of the Hamiltonian
    + h (double) parameter $h$ of the Hamiltonian
    + tol (double) tolerance asked from Contractor
    + mpsfile (string) file where the resulting MPS will be stored.
    + oldmpsfile [OPTIONAL] (string) name of another (binary) file, containing 
                    an earlier MPS to be used as starting point

		    \todo: what do I compute afterwards. Infinite chain evolution?

*/

using namespace std;
using namespace shrt;

// /** Given a state, compute all the expectation values (on each
//     position) of a given single body operator. The results are
//     written to the ostream (could be a file, or cout). */
// void computeSingleBodyEV(const MPS& rho,const mwArray& op1,ostream& os);

// /** Same as before, but compute the expectation value of the product
//     of op1 acting on one site and op2 acting on the site plus l.  If
//     maxLen!=0 is provided, it writes up to maxLen elements in the
//     line, padding with 0 when needed.  
// */
// void computeTwoBodyEV(const MPS& rho,const mwArray& op1,
// 		      const mwArray& op2,int l,ostream& os,
// 		      int maxLen=0);

// /** Prepare header for a file */
// void writeHeader(ofstream& out,const CoherentDissipationLiouvillian& L,
// 		 double lambda0,
// 		 int D,double scal,double costT,const string oldMPSfile);

/** Global matrices */
int d=2;
complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
mwArray sigX(Indices(d,d),dataX);
complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
mwArray sigY(Indices(d,d),dataY);
complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
mwArray sigZ(Indices(d,d),dataZ);

//double noise=1E-60;

int main(int argc,const char* argv[]){
  int cntr=0;
  const char* infile=argv[++cntr];
  string directory="";
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  int L=props.getIntProperty("L");
  int D=props.getIntProperty("D");
  double J=props.getDoubleProperty("J");
  double g=props.getDoubleProperty("g");
  double h=props.getDoubleProperty("h");
  double tol=props.getDoubleProperty("tol");
  const string newMPSfile=props.getProperty("mpsfile");
  const string oldMPSfile=props.getProperty("oldmpsfile");
  bool usingOldMPS=(!oldMPSfile.empty()&&oldMPSfile.length()>1);
  // const string outdir=props.getProperty("outdir");
  // double noise=props.getDoubleProperty("noise");
  // if(noise<0.) noise=0.;

  cout<<"oldMPSfile is "<<oldMPSfile<<endl;
  cout<<"its size is "<<oldMPSfile.size()<<endl;
  cout<<"its length is "<<oldMPSfile.length()<<endl;


  // Construct the Hamiltonian (on L+ boundary)
  IsingHamiltonian H(L+2,d,J,g,h);
  // Commutator superoperator
  MPO commMPO(L+2);
  H.getCommutatorMPO(commMPO);
  
  cout<<"Constructed commutator"<<endl;

  // Prepare an MPS for the initial guess of the optimization
  // But now, the actual A has two extra sites, which are just constants, so I create an MPS with L+2 sites
  mwArray edge=identityMatrix(d);edge.reshape(Indices(d*d,1,1));
  vector<int> bonds(L+1,1);
  vector<int> dims(L+2,d*d);
  for(int k=0;k<L+1;k++) bonds[k]=D; // 0 and L are 1
  MPS A(L+2,bonds,dims);
  // And set the inital (L middle) tensors to the initial MPS
  if(usingOldMPS){ // if read from file, copy them
    MPS _A_(L,D,d*d);
    _A_.importMPS(oldMPSfile.data());
    cout<<"Read initial A from file "<<oldMPSfile<<endl;
    for(int k=0;k<L;k++)
      A.replaceSite(k+1,_A_.getA(k));
  }
  else{
    A.setRandomState();
    cout<<"Initial A set to random MPS"<<endl;
  }

  // And now try to optimize over A with Contractor::findGroundState
  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  if(D<=12)
    contractor.setEigenSolver(fulleig);
  contractor.setConvTol(tol);

  double lambda=0;
  contractor.findGroundState(commMPO,D,&lambda,A);

  cout<<"Finished GS search with lambda="<<lambda<<endl;

  if(!newMPSfile.empty()&&newMPSfile.length()>1)
    // Now I have the GS, I will write it to a file
    A.exportMPS(newMPSfile.data());


}

// void computeSingleBodyEV(const MPS& rho,const mwArray& op1,ostream& os){
//   int L=rho.getLength();
//   MPS aux(L,1,d*d);
//   mwArray basic=identityMatrix(d); // basic piece of the MPO
//   basic.reshape(Indices(d*d,1,1));
//   for(int k=0;k<L;k++){
//     aux.setA(k,basic);
//   }

//   Contractor& contractor=Contractor::theContractor();
//   double traceRho=real(contractor.contract(rho,aux));

//   // WARNING: No check of dimensions done!!
//   mwArray op1_=reshape(op1,Indices(d*d,1,1));
//   op1_.conjugate(); // It gets conjugated in the scalar product! 

//   // Now, site by site, replace the identity by the operator,
//   // contract, write the result to os, and replace the operator by the
//   // identity again
//   for(int k=0;k<L;k++){
//     aux.setA(k,op1_); // operator on site k
//     complex_t valK=contractor.contract(rho,aux);
//     os<<real(valK)/traceRho<<"\t";
//     aux.setA(k,basic); // operator removed
//   }
// }

// void computeTwoBodyEV(const MPS& rho,const mwArray& op1,
// 		      const mwArray& op2,int l,ostream& os,
// 		      int maxLen){
//   int L=rho.getLength();
//   MPS aux(L,1,d*d);
//   mwArray basic=identityMatrix(d); // basic piece of the MPO
//   basic.reshape(Indices(d*d,1,1));
//   for(int k=0;k<L;k++){
//     aux.setA(k,basic);
//   }

//   Contractor& contractor=Contractor::theContractor();
//   double traceRho=real(contractor.contract(rho,aux));

//   // WARNING: No check of dimensions done!!
//   mwArray op1_=reshape(op1,Indices(d*d,1,1));
//   op1_.conjugate(); // It gets conjugated in the scalar product! 
//   mwArray op2_=reshape(op2,Indices(d*d,1,1));
//   op2_.conjugate(); // It gets conjugated in the scalar product! 

//   // Now, site by site, replace the identity at k by the operator op1,
//   // the identity in k+l by op2,
//   // contract, write the result to os, and replace the operators by the
//   // identity again
//   int cnt=0; // nr of terms written
//   for(int k=0;k+l<L;k++){
//     aux.setA(k,op1_); // operator on site k
//     aux.setA(k+l,op2_); // operator on site k+l
//     complex_t valK=contractor.contract(rho,aux);
//     os<<real(valK)/traceRho<<"\t";
//     cnt++;
//     aux.setA(k,basic); // operator1 removed
//     aux.setA(k+l,basic); // operator2 removed
//   }
//   if(maxLen>0)
//     while(cnt<maxLen){
//       os<<0.<<"\t";
//       cnt++;
//     }

// }

// void writeHeader(ofstream& out,const CoherentDissipationLiouvillian& L,
// 		 double lambda0,
// 		 int D,double scal,double costT,const string oldMPSfile){
//   out<<"% L="<<L.getLength()<<endl;
//   out<<"% g="<<L.getG()<<endl;
//   out<<"% gamma="<<L.getGamma()<<endl;
//   out<<"% D="<<D<<endl;
//   out<<"% scale="<<scal<<endl;
//   Contractor& contr=Contractor::theContractor();
//   out<<"% Contractor tol="<<contr.getConvTol()<<endl;
//   out<<"% Contractor solver="<<contr.getSolverName()<<endl;
//   out<<"% Time (s) needed = "<<costT<<endl;
//   if(!oldMPSfile.empty())
//     out<<"% Using initial state from "<<oldMPSfile<<endl;
//   else
//     out<<"% Using random initial state "<<endl;
//   out<<"% Lowest eigenvalue found to be "<<lambda0<<endl;
// }

