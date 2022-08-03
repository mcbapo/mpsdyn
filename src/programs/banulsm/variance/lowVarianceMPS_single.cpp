
#include <math.h>
#include <iomanip>

#include "mwArray.h"
#include "MPS.h"
#include "misc.h"
#include "Contractor.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

/** Runs the minimizeVariance routine with the MPO for the Ising
    Hamiltonian \ref <IsingHamiltonian> at a given energy, and saves
    the results (MPS) in a file.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <E0> (double) energy value
    \param <penH> (double) penalty fort the energy term
    \param <outfname> (char*) name of the output file for the MPS
*/


int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  int L=props.getIntProperty("L");
  double J=props.getDoubleProperty("J");
  double h=props.getDoubleProperty("h");
  double g=props.getDoubleProperty("g");
  double E0=props.getDoubleProperty("E");
  double penH=props.getDoubleProperty("penH");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  int D=props.getIntProperty("D");
  string mpsfile=props.getProperty("mpsfile");
  string initmpsfile=props.getProperty("initmpsfile");
  string outputfile=props.getProperty("outputfile");

  cout<<"Initialized arguments: L="<<L
      <<", J="<<J
      <<", g="<<g
      <<", h="<<h
      <<", initmpsfile="<<initmpsfile
      <<", D="<<D
      <<", E="<<E0<<endl;

  cout<<setprecision(10);

  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  contractor.setConvTol(tol);
  cout<<"Initialized Contractor"<<endl;
  string tmpfile(mpsfile);tmpfile.append(".tmp");
  int d=2;
  MPS gs(L,D,d);
  bool readInit=false; // whether I success in reading an initial guess
  // FIRST: Check for tmp file
  if(file_exists(tmpfile)){
    gs.importMPS(tmpfile.data());
    cout<<"Imported initial state from tmp file "<<tmpfile<<endl;
    readInit=true;
  }
  else{ // see if initmpsfile was specified
    if(initmpsfile.size()>0){
      if(file_exists(initmpsfile)){    
	gs.importMPS(initmpsfile.data());    readInit=true;
	cout<<"Imported initial state from file "<<initmpsfile<<endl;
      }
      else{// check for a tmp file of the initmpsfile
	string initTmp(initmpsfile);initTmp.append(".tmp");
	if(file_exists(initTmp)){
	  gs.importMPS(initTmp.data());    readInit=true;
	  cout<<"Imported initial state from file "<<initTmp<<endl;
	}
      }
    }
  }
  if(!readInit){
    gs.setRandomState(); // the intial state, random
    cout<<"Initialized random state"<<endl;
  }
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  
  // First: put H and find the GS
  IsingHamiltonian hamI(L,d,J,g,h,0.);//-E); 
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamI.getHMPO(0.);
  cout<<"Constructed the hamil MPO"<<endl;
  
  // place for other params
  complex_t E=contractor.contract(gs,hamil,gs);

  complex_t E2=contractor.contract2(hamil,gs);

  //hamil.exportForMatlab("isingHmpo.m");
  cout<<"Initial value of E, with initial state"<<endl;
  cout<<E<<endl;
  cout<<"and <H^2>="<<E2<<endl;

  // const MPO* ptrs[]={&hamil,&hamil};
  // MPO ham2(L);
  // MPO::join(2,ptrs,ham2);// MPO for H^2
  
  double varE=real(E2)-real(E)*real(E);
  cout<<"Initial variance: "<<varE<<endl;
  //  contractor.findGroundState(ham2,D,&expH2,gs);
  contractor.minimizeVariance(gs,D,hamil,E0,penH,&varE,tmpfile,500);

  E=contractor.contract(gs,hamil,gs);
  E2=contractor.contract2(hamil,gs);
  cout<<"After the optimization, cost function is "<<varE<<endl;
  cout<<"Variance="<<E2-E*E<<" <H^2>="<<E2<<", <H>="<<E<<endl;

  if(mpsfile.size()>0){
    gs.exportMPS(mpsfile.data());
    cout<<"Exported solution to file "<<mpsfile<<endl;
  }
  remove(tmpfile.data());
  
    //    *out<<"% N\t J\t g\t h\t D\t <H^2>\t <H_L^2>\t <H_R^2> \t <H>\t <HL>\t <HR>"<<endl;
  if(outputfile.size()>0){
    ofstream* out=new ofstream(outputfile.data(),ios::app);
    *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"
	<<varE<<"\t"<<real(E2)<<"\t"<<real(E)<<"\t"<<real(E2-E*E);
    *out<<endl;
    out->close();
    delete out;
  }
}
