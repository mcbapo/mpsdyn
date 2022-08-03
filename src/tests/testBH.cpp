#ifndef NEWIMPLEM
#include "libdiminf.h"
#include "diminf.h"
#else
#include "mwArray.h"
#include "MPO.h"
#include "BoseHubbardHamiltonian.h"
#endif
#include "Operator.h"
#include "Contractor.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>


#ifndef NEWIMPLEM
typedef OperatorRow MPO;
#endif


#ifdef MACOSX
#define LOCALDIR "/Users/banulsm/Desktop/SVN/projects/mpsdyn/src/bin/"
#else
#define LOCALDIR "/home/banulsm/transfer/"
#endif


/** 
    Test program Nr. 12: Check the correct behaviour of BH Hamiltonian
    To compile: make test12
*/

#ifdef NEWIMPLEM
using namespace std;
using namespace shrt;
#endif

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  int N=atoi(argv[++cntr]);
  int xL=atoi(argv[++cntr]);
  int xR=atoi(argv[++cntr]);
  double t=atof(argv[++cntr]);
  double UL=atof(argv[++cntr]);
  double UC=atof(argv[++cntr]);
  double UR=atof(argv[++cntr]);
  double mu=atof(argv[++cntr]);
  double VL=atof(argv[++cntr]);
  double VR=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);

  char filename[120];

  int d=N+1;

#ifndef NEWIMPLEM
  int error;
  initializelib(&error);
  try{
    initGlobals(mwArray(CHECKDIMS),mwArray(SEEDVAL));
#endif

#ifdef NEWIMPLEM
  double* tis=new double[L];
  double* Uis=new double[L];
  double* muis=new double[L];
  double* Vis=new double[L];
  for(int k=0;k<L;k++){
    tis[k]=t;
    Uis[k]=k<xL?UL:(k<=xR?UC:UR);
    muis[k]=mu;
    Vis[k]=k<xL?VL:(k<=xR?0:VR);
    cout<<"Site "<<k<<", U="<<Uis[k]<<",V="<<Vis[k]<<endl;
  }

  BoseHubbardHamiltonian hBH(L,N,tis,Uis,muis,Vis);
  const MPO& mpoBH=hBH.getHMPO();

  strcpy(filename,LOCALDIR);
  mpoBH.exportMPOtext(strcat(filename,"mpoBH.txt"));
#else
  MPO mpoBH(L,4,d);
  strcpy(filename,LOCALDIR);
  mpoBH.importMPOtext(strcat(filename,"mpoBH.txt"));

#endif

    MPS gsP(L,D,d);
    gsP.setRandomState(); // the initial state, random -> not to the GS!!
    gsP.setRandomState(); // the initial state, random -> not to the GS!!
    gsP.gaugeCond('R',1);
    gsP.gaugeCond('L',1);
    double lambda=0.;

#ifdef NEWIMPLEM    
    strcpy(filename,LOCALDIR);
    gsP.exportMPStext(strcat(filename,"mpsBH0.txt"));
#else
    strcpy(filename,LOCALDIR);
    gsP.importMPStext(strcat(filename,"mpsBH0.txt"));
#endif

    Contractor& contractor=Contractor::theContractor();
    contractor.findGroundState(mpoBH,D,&lambda,gsP);

    int Dop=d*d;
#ifdef NEWIMPLEM    
    // Test two body H12
    mwArray h12;
    int pos=0;complex_t delta={.27,0.};
    hBH.computeTwoBodyTerm(h12,pos);
    cout<<"Computed two-body term"<<h12<<endl;
    mwArray Op1,Op2;
    hBH.getTwoBodyTermExponential(Op1,Op2,delta,pos);
    cout<<"Computed exponential opers Op1="<<Op1
	<<", Op2="<<Op2<<endl;

    // Compute also the MPO
    MPO expHe(L),expHo(L);
    hBH.getExponentialMPOeven(expHe,delta);
    hBH.getExponentialMPOodd(expHo,delta);
    strcpy(filename,LOCALDIR);
    expHe.exportMPOtext(strcat(filename,"expMPOe.txt"));
    strcpy(filename,LOCALDIR);
    expHo.exportMPOtext(strcat(filename,"expMPOo.txt"));

#else
    // Could do this through checkMPO 
    MPO expHe(L,Dop,d),expHo(L,Dop,d);
    expHe.importMPOtext("/home/banulsm/transfer/expMPOe.txt");
    expHe.saveOperatorArray("/home/banulsm/transfer/matlab/expMPOe.mat");
    expHo.importMPOtext("/home/banulsm/transfer/expMPOo.txt");
    expHo.saveOperatorArray("/home/banulsm/transfer/matlab/expMPOo.mat");

#endif
    

    cout<<"Ground state found with eigenvalue "<<lambda<<endl;
#ifndef NEWIMPLEM
    mpoBH.saveOperatorArray("/home/banulsm/transfer/mpoBH.mat");
  }
  catch(mwException& e){
    cout<<"Exception caught: "<<e.what()<<"in testMPO"<<endl;
    exit(212);
  }
  closelib();
  return error;
#endif
}
