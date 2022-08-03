
#include <math.h>
#include <iomanip>

#include "mwArray.h"
#include "MPS.h"
#include "misc.h"
#include "Contractor.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "HeisenbergHamiltonian.h"

using namespace shrt;

/** 
Generate the name of the file with the command in it
*/
const string jobfilename(int L,int D,double E,double Jx,double Jy,double Jz,double h,const string jobsdir);
/** 
Generates a file name for the MPS corresponding to a certain energy.
*/
const string mpsfilename(int L,int D,double E,double Jx,double Jy,double Jz,double h,const string mpsdir);

/** Prepare a job script with the instruction to launch another instance  */
void writeJob(const string jobsdir,const string exec,const string configfile,
	      int L,int D,double E,double Jx,double Jy,double Jz,double h,double tol,
	      double Emin,double Emax,int nrStates,const string mpsdir,double penH,
	      const string initmpsfile,int Dmin,int Dmax,int stepD,bool greedy=0);

/** Runs the minimizeVariance routine with the MPO for the Heisenberg
    Hamiltonian \ref <HeisenbergHamiltonian> at a given energy, and saves
    the results (MPS) in a file.
    It also launches the next one (creates a script that does) until
    the threshold E=0 is reached (so that we can start from both edges
    of the spectrum. 
    The first time it is called (without a value for the min and max
    energies) it locates the extrema of the spectrum and prepares two
    scripts.

*/

int d=2;

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
  double Jx=props.getDoubleProperty("Jx");
  double Jy=props.getDoubleProperty("Jy");
  double Jz=props.getDoubleProperty("Jz");
  double h=props.getDoubleProperty("h");
  double E0=props.getDoubleProperty("E");
  double Emin=props.getDoubleProperty("Emin");
  double Emax=props.getDoubleProperty("Emax");
  int nrStates=props.getDoubleProperty("nrStates"); // how many to compute in total (incl. extremes)
  // if nrStates is 0, the tower will be done only for the specified E, but not launched for other energies
  bool firstCall=(Emin==-1)&&(Emax==-1);
  double penH=props.getDoubleProperty("penH");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  int D=props.getIntProperty("D");
  int Dmin=props.getIntProperty("Dmin"); // the smallest D
  int Dmax=props.getIntProperty("Dmax"); // where to stop
  int stepD=props.getIntProperty("stepD");
  if(stepD<0) stepD=20;
  string mpsdir=props.getProperty("mpsdir");
  string initmpsfile=props.getProperty("initmpsfile");
  //  string outputfile=props.getProperty("outputfile");
  string jobsdir=props.getProperty("jobsdir");
  bool launchImmediately=(props.getIntProperty("greedy")==1); // whether to launch the next immediately!
  
  cout<<"Initialized arguments: L="<<L
      <<", Jx="<<Jx
      <<", Jy="<<Jy
      <<", Jz="<<Jz
      <<", h="<<h
      <<", initmpsfile="<<initmpsfile
      <<", D="<<D
      <<", E="<<E0<<endl;

  cout<<setprecision(10);

  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  contractor.setConvTol(tol);
  cout<<"Initialized Contractor"<<endl;
  

  // First: put H and find the GS
  HeisenbergHamiltonian hamH(L,Jx,Jy,Jz,h,d,0.);//-E); 
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamH.getHMPO();
  cout<<"Constructed the hamil MPO"<<endl;

  if(firstCall){
    // If the first call, locate edges of the spectrum
    Emin=0;
    MPS gs(L,D,d);
    gs.setRandomState();gs.gaugeCond('R');gs.gaugeCond('L',1);
    contractor.findGroundState(hamil,D,&Emin,gs);
    cout<<"Found GS of H with E="<<Emin<<endl;
    // Save as first state
    string gsfilename=mpsfilename(L,D,Emin,Jx,Jy,Jz,h,mpsdir);
    gs.exportMPS(gsfilename.data());
    
    // Same for largest eigenvalue:
    Emax=0;
    MPS exc(L,D,d);
    exc.setRandomState();exc.gaugeCond('R');exc.gaugeCond('L',1);
    {
      MPO hamilMinus(L);
      hamilMinus.setOp(0,new Operator(-1.*hamil.getOp(0).getFullData()),true);
      for(int k=1;k<L;k++){
	hamilMinus.setOp(k,&hamil.getOp(k),false);
      }
      contractor.findGroundState(hamilMinus,D,&Emax,exc);
      Emax=-Emax;
    }
    cout<<"Found max Exc of H with E="<<Emax<<endl;
    // Save and prepare jobs for both extremes
    string excfilename=mpsfilename(L,D,Emax,Jx,Jy,Jz,h,mpsdir);
    exc.exportMPS(excfilename.data());

    if(nrStates>1){
      // Now the step of energy
      double deltaE=(Emax-Emin)/(nrStates-1);
      // Launch four jobs: two inwards in E and two upwards in D
      // The first ones cannot use this as initial state (too close to eigenstate)
      writeJob(jobsdir,argv[0],argv[1],L,D,Emin+deltaE,Jx,Jy,Jz,h,tol,Emin,Emax,nrStates,
	       mpsdir,penH,"",Dmin,Dmax,stepD,launchImmediately);
      writeJob(jobsdir,argv[0],argv[1],L,D,Emax-deltaE,Jx,Jy,Jz,h,tol,Emin,Emax,nrStates,
	       mpsdir,penH,"",Dmin,Dmax,stepD,launchImmediately);
      // now the ones in D
      if(D+stepD<=Dmax){
	writeJob(jobsdir,argv[0],argv[1],L,D+stepD,Emin,Jx,Jy,Jz,h,tol,Emin,Emax,nrStates,
		 mpsdir,penH,gsfilename,Dmin,Dmax,stepD,launchImmediately);
	writeJob(jobsdir,argv[0],argv[1],L,D+stepD,Emax,Jx,Jy,Jz,h,tol,Emin,Emax,nrStates,
		 mpsdir,penH,excfilename,Dmin,Dmax,stepD,launchImmediately);
      }
    }
    else{ // launch just one, this time knowing the edges
      writeJob(jobsdir,argv[0],argv[1],L,D,E0,Jx,Jy,Jz,h,tol,Emin,Emax,nrStates,
	       mpsdir,penH,"",Dmin,Dmax,stepD,launchImmediately);
    }
  }
  else{ // Compute just one state
    MPS gs(L,D,d);
    string mpsfile=mpsfilename(L,D,E0,Jx,Jy,Jz,h,mpsdir);    
    string tmpfile(mpsfile);tmpfile.append(".tmp");

    if(launchImmediately){
      // Launch the next jobs immediately, without waiting for this one!
      // write the job script
      // increase D
      if(D+stepD<=Dmax){
	writeJob(jobsdir,argv[0],argv[1],L,D+stepD,E0,Jx,Jy,Jz,h,tol,Emin,Emax,nrStates,
		 mpsdir,penH,mpsfile,Dmin,Dmax,stepD,launchImmediately);
      }
      if(D==Dmin&&nrStates>1){// then move in energy, too }
	double deltaE=(Emax-Emin)/(nrStates-1);
	double nextE=E0<0?E0+deltaE:E0-deltaE;
	if(E0*nextE>=0||abs(nextE)<.5*deltaE){ // both on the same side. 
	  writeJob(jobsdir,argv[0],argv[1],L,D,nextE,Jx,Jy,Jz,h,tol,Emin,Emax,nrStates,
		   mpsdir,penH,mpsfile,Dmin,Dmax,stepD,launchImmediately);
	}    
      }  
    }

    bool readInit=false; // whether I succeeded in reading an initial guess
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
      cout<<"Initialized random state because I could not find the initfile "<<initmpsfile<<endl;
    }
    gs.gaugeCond('R',1);
    gs.gaugeCond('L',1);
  
  
    complex_t E=contractor.contract(gs,hamil,gs);
    complex_t E2=contractor.contract2(hamil,gs);

    //hamil.exportForMatlab("isingHmpo.m");
    cout<<"Initial value of E, with initial state"<<endl;
    cout<<E<<endl;
    cout<<"and <H^2>="<<E2<<endl;
    
    double varE=real(E2)-real(E)*real(E);
    cout<<"Initial variance: "<<varE<<endl;
    contractor.minimizeVariance(gs,D,hamil,E0,penH,&varE,tmpfile,200);

    E=contractor.contract(gs,hamil,gs);
    E2=contractor.contract2(hamil,gs);
    cout<<"After the optimization, cost function is "<<varE<<endl;
    cout<<"Variance="<<E2-E*E<<" <H^2>="<<E2<<", <H>="<<E<<endl;

    if(mpsfile.size()>0){
      gs.exportMPS(mpsfile.data());
      cout<<"Exported solution to file "<<mpsfile<<endl;
    }
    remove(tmpfile.data());
    if(!launchImmediately){
      // And write the job script
      // increase D
      if(D+stepD<=Dmax){
	writeJob(jobsdir,argv[0],argv[1],L,D+stepD,E0,Jx,Jy,Jz,h,tol,Emin,Emax,nrStates,
		 mpsdir,penH,mpsfile,Dmin,Dmax,stepD);
      }
      if(D==Dmin&&nrStates>1){// then move in energy, too }
	double deltaE=(Emax-Emin)/(nrStates-1);
	double nextE=E0<0?E0+deltaE:E0-deltaE;
	if(E0*nextE>=0){ // both on the same side. WARNING: 0 maybe left out
	  writeJob(jobsdir,argv[0],argv[1],L,D,nextE,Jx,Jy,Jz,h,tol,Emin,Emax,nrStates,
		   mpsdir,penH,mpsfile,Dmin,Dmax,stepD);
	}    
      }  
    }
  }
}



const string mpsfilename(int L,int D,double E,double Jx,double Jy,double Jz,double h,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/mps_N"<<L<<"_Jx"<<Jx<<"_Jy"<<Jy<<"_Jz"<<Jz<<"_h"<<h<<"_E"<<E<<"_D"<<D<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

const string jobfilename(int L,int D,double E,double Jx,double Jy,double Jz,double h,const string jobsdir){
  stringstream s;
  s<<jobsdir<<"/job_L"<<L<<"_Jx"<<Jx<<"_Jy"<<Jy<<"_Jz"<<Jz<<"_h"<<h<<"_E"<<E<<"_D"<<D<<".dat";
    //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

void writeJob(const string jobsdir,const string exec,const string configfile,
	      int L,int D,double E,double Jx,double Jy,double Jz,double h,double tol,
	      double Emin,double Emax,int nrStates,const string mpsdir,double penH,
	      const string initmpsfile,int Dmin,int Dmax,int stepD,bool greedy){
  string jobfile=jobfilename(L,D,E,Jx,Jy,Jz,h,jobsdir);
  stringstream commandstr;
  commandstr<<exec;
  commandstr<<" "<<configfile;
  commandstr<<" -L="<<L<<" -Jx="<<Jx<<" -Jy="<<Jy<<" -Jz="<<Jz<<" -h="<<h;
  commandstr<<" -D="<<D<<" -Dmin="<<Dmin<<" -Dmax="<<Dmax<<" -stepD="<<stepD;
  commandstr<<" -E="<<E<<" -penH="<<penH<<" -tol="<<tol;
  commandstr<<" -Emin="<<Emin<<" -Emax="<<Emax<<" -nrStates="<<nrStates;
  commandstr<<" -mpsdir="<<mpsdir<<" -jobsdir="<<jobsdir;
  if(initmpsfile.size()>0)
    commandstr<<" -initmpsfile="<<initmpsfile;
  if(greedy)
    commandstr<<" -greedy=1";

  ofstream* ojob=new ofstream(jobfile.data());
  if(!ojob->is_open()){
    cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
    cout<<commandstr.str()<<endl;
    cout<<endl;
  }
  else{
    *ojob<<"#!/bin/bash"<<endl;
    *ojob<<"msub_slurm";
    if(D<=60) *ojob<<" -P 2";
    else if(D>60&&D<100) *ojob<<" -P 4";
    else if(D>=100&&D<=140)  *ojob<<" -P 12";
    else *ojob<<" -P 14";
    *ojob<<" -N loVa.L"<<L<<".D"<<D<<".E"<<E<<".Jx"<<Jx<<".Jy"<<Jy<<".Jz"<<Jz<<".h"<<h<<" ";
    *ojob<<" -- ";
    *ojob<<commandstr.str()<<endl;
    ojob->close();
  }
  delete ojob;
}
