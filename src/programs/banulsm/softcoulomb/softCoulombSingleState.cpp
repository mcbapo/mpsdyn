
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"
//#include "HermitianTensorMultiplier.h"

#include "SoftCoulombHamiltonian.h"
//#include "SpinMPO.h"
#include <cmath>

using namespace std;
using namespace shrt;

#define MAXLEN 160

/** 
Generates a file name for the MPS corresponding to a certain excitation.
*/
const string mpsfilename(int L,int D,int M,int Ze,int level,const string mpsdir);
/** 
Generate the name of the file with the command in it
*/
const string jobfilename(int L,int D,int M,int Ze,const string mpsdir);


/** Computes a single eigenvector of the soft Coulomb Hamiltonian
    \ref <SoftCoulombHamiltonian>

    Receives as argument a Properties file, such as
    config/softCoulombProp.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.

    If property minnumber is set to a certain value, the program will
    run until it has those many levels computed, and then will invoke
    searchOrder to order them.

    It is possible to read existing MPS files from another directory,
    if they do not exist in the destination one, and use them as
    initial states. For this, specify the name of the old directory
    initmpsdir in Properties. If the files for a smaller D are to be
    used, specify the value in initD (in that case, initmpsdir could
    be the same as mpsdir). A noise parameter can be specified in the
    list of properties, and then the extra components of the tensors
    in the initial guess are filled with some noise.

    Moreover, the initial state for the search can be constructed from
    a shorter MPS, by stretching it (assuming homogeneous tensors in
    the center). To do this, set stretchmps=1 in Properties, and
    specify initL for the intial length of the MPS to be
    stretched. The directory (initmpsdir) and bond dimension (initD)
    as before.

    Property jobsdir will indicate where to write the file with the 
    command that has to be relaunched 

    Property onlygs makes the thing stop after the GS, and prepare 
    a job to run for the next bond dimension (+20)
    until a maximum one (taken to be 140)

    To replace some of the values in the file by a command line
    option, the additional arguments have to be of the form
    -propName=propValue
    (except for the vector properties)
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

  // Store the command line, for this is what should be relaunched if not all levels are done (and at the end for observables)
  //  stringstream commandstr;
  //for(int c=0;c<argc;c++) commandstr<<argv[c]<<" ";
  // Now, str contains the command line, that I need to relaunch, if this was not the last run

  int L = props.getIntProperty("L");
  double M=props.getDoubleProperty("M");
  double Z=props.getDoubleProperty("Ze");
  double a=props.getDoubleProperty("a");
  double Delta=props.getDoubleProperty("Delta");
  //  int D0 = props.getIntProperty("D0");
  //  int nLev = props.getIntProperty("nLev");
  directory=props.getProperty("outputdir");
  const string outfile=directory+"/"+props.getProperty("outfile");
  int app=props.getIntProperty("append");
  if(app<0) app=1; // default will be to append
  double tol=props.getDoubleProperty("tol");
  // int unif=1; // All same D 
  // int unif=props.getIntProperty("unif");
  double Npen=props.getDoubleProperty("penaltyNe");
  int Ntot=props.getIntProperty("Ne");
  double Spen=props.getDoubleProperty("penaltyS");
  if(Spen<0) Spen=0;
  int D=props.getIntProperty("D");
  double offset=props.getDoubleProperty("offset");
  string mpsdir=props.getProperty("mpsdir"); // default current
  if(mpsdir.empty()) 
    mpsdir=".";
  int nLev=props.getIntProperty("nLev");
  if(nLev<0) nLev=0; // where to stop
  int stopAtGS=(nLev==0);
  string initmpsdir=props.getProperty("initmpsdir"); // default none
  if(initmpsdir.empty()){
    initmpsdir=props.getProperty("oldmpsdir"); // default none
  }
  int initD=props.getIntProperty("initD");
  int stretchmps=props.getIntProperty("stretchmps");
  if(stretchmps==-1) stretchmps=0; 
  int initL=props.getIntProperty("initL");
  double noise=props.getDoubleProperty("noise");
  if(noise<0) noise=0.;
  if(!initmpsdir.empty()){
    if(initD<0){
      cout<<"WARNING: No initD specified, using the same one as for this run ("<<D<<")"<<endl;
      initD=D;
    }
    if(initL==-1){
      if(stretchmps)
	cout<<"WARNING: stretchmps set to true, but no initial length given. Assuming "
	    <<"the same one, so no stretching will take place "<<endl;
      initL=L;
    }
  }
  else{
    cout<<"WARNING: No initmpsdir specified => no initial guess for MPS will be used "<<endl;
  }


  if(offset==-1) offset=0;

  string jobsdir=props.getProperty("jobsdir"); 
  int maxTime=props.getIntProperty("maxTime");
  if(maxTime==-1) maxTime=0;
  int freqSv=props.getIntProperty("savingTmp");
  if(freqSv==-1) freqSv=5;

  // if(stopAtGS){ // If only computing GS, Iwill then use larger D
  //   commandstr.str(""); // empty again
  //   commandstr<<argv[0];
  //   commandstr<<" "<<argv[1];
  //   commandstr<<" -L="<<L<<" -M="<<M<<" -Delta="<<Delta;
  //   commandstr<<" -a="<<a<<" -D="<<D+20<<" -Ze="<<Z;
  //   commandstr<<" -append="<<app<<" -outfile=L"<<L<<"Z"<<Z<<"D"<<D+20;
  //   commandstr<<" -outputdir="<<directory<<" -tol="<<tol;
  //   commandstr<<" -penaltyNe="<<Npen<<" -Ne="<<Ntot;
  //   commandstr<<" -penaltyS="<<Spen;
  //   commandstr<<" -offset="<<offset;
  //   commandstr<<" -mpsdir="<<mpsdir;
  //   commandstr<<" -initmpsdir="<<mpsdir; // what I just computed is the new initmpsdir
  //   commandstr<<" -initD="<<D<<" ";
  //   commandstr<<" -noise="<<noise;
  //   commandstr<<" -jobsdir="<<jobsdir;
  //   commandstr<<" -maxTime="<<maxTime;
  //   commandstr<<" -nLev=0 ";
  //   int mem=ceil(L*(D+20)*(D+20)*(4*16.)*16E-6)+200;
  //   commandstr<<" -msubmem="<<mem<<"M";
  // }

  cout<<"Initialized arguments: L="<<L
      <<", Delta="<<Delta
      <<", Ze="<<Z
      <<", a="<<a
      <<", M="<<M
      <<", outfile="<<outfile
      <<", app="<<app
      <<", tol="<<tol
      <<", penaltyNe="<<Npen
      <<", Ne="<<Ntot
      <<", D0="<<D<<endl;

#ifdef DISKSTRG
  //  int errD=system("rm -rf /ptmp/mpq/banulsm/*");
  char tmpdir[150];
  sprintf(tmpdir,"/ptmp/mpq/banulsm/SC/Sitesx%dL%dD%dM%dXXXXXX",x,L,D,M);
  char* dirid=mkdtemp(tmpdir);
  if(dirid==0){
    cout<<"Error: couldn't apply mkdtemp "<<tmpdir<<endl;
    exit(1);
  }
  FileSite::setDir(tmpdir);
#endif

  ofstream* out;

  out=new ofstream(outfile.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output (append="<<app<<")"<<endl;
    exit(1);
  }
  *out<<setprecision(15);

  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);
  double lambda=0.;
  int d=2;
  SoftCoulombHamiltonian H(L,Delta,Z,a,M,Npen,Ntot,Spen);
  const MPO& hamil=H.getHMPO();
  // Observables I will compute
  // Compute the expectation value of N, as it commutes with H
  MPO Nmpo(L);
  H.getNumberMPO(Nmpo);
  MPO Smpo(L);
  H.getSpinMPO(Smpo);
  // if(L<=5){
  //   hamil.exportForMatlab("HSC.m");
  //   Smpo.exportForMatlab("Smpo.m");
  //   Nmpo.exportForMatlab("Nmpo.m");
  // }

  vector<MPS*> levels;

  // First, check how many MPS I already have in mpsdir
  bool oldMPS=true;
  int lval=0;
  MPS* nextLev;
  while(oldMPS&&(lval<=nLev)){
    const string filename=mpsfilename(L,D,M,Z,lval,mpsdir);
    cout<<"Checking for file ("<<lval<<") "<<filename<<endl;
    // Try to open file for MPS level lval
    if(file_exists(filename)){
      nextLev=new MPS(L,D,d);
      nextLev->importMPS(filename.data());
      levels.push_back(nextLev);
      cout<<"Recovered MPS file for level "<<lval<<" from file "<<filename<<endl;
      lval++;
    }
    else{
      oldMPS=false; // break the loop 
    }
  }

  if(stopAtGS&&lval>0){
    cout<<"Stopping at GS, which is already computed. Nothing to do"<<endl;
    out->close();
    delete out;
    exit(333);
  }
  if(lval>nLev){
    cout<<"Stopping at "<<nLev<<", which is already computed. Nothing to do"<<endl;
    out->close();
    delete out;
    exit(333);
  }


  // Now I already have read all existing MPS
  cout<<"Recovered "<<lval<<" already computed MPS levels from disk. "
      <<"Computing level "<<lval<<" now (GS is 0)."<<endl;


  // values for each level
  complex_t Hval,H2val,normN,Nval,Sval;
  int d2=d*d;int Dh=6+M; // dim of H mpo
  MPS* exc=new MPS(L,D,d2);
  bool tmpRes=false;
  // First check if there is any temporary result
  const string tmpfilename=mpsfilename(L,D,M,Z,lval,mpsdir)+".tmp";
  if(!file_exists(tmpfilename)){
    const string initfilename=mpsfilename(initL,initD,M,Z,lval,initmpsdir);
    if(initmpsdir.empty())
      exc->setRandomState();
    else{
      cout<<"Will search initial state for level "<<lval<<" in "<<initfilename<<endl;
      if(!file_exists(initfilename)){
	cout<<"ERROR: Initial MPS not found in "<<initfilename<<endl;
	exit(1);
      }
      if(initL<L){
	MPS auxMPS(initL,initD,d2);
	auxMPS.importMPS(initfilename.data());
	cout<<"Stretching initial state (l="<<lval<<") from "<<initfilename<<endl;
	exc->stretchMPS(L,auxMPS);
      }
      else{
	exc->importMPS(initfilename.data());
	cout<<"Recovered initial state for level "<<lval<<" from file "
	    <<initfilename<<endl;
      }
      if(initD<D&&noise!=0.)
	exc->increaseBondDimensionWithNoise(D,noise);
    }
  }
  else{
    cout<<"Importing temporary result from "<<tmpfilename<<endl;
    exc->importMPS(tmpfilename.data());
    tmpRes=true;
    cout<<"Recovered temporary state for level "<<lval<<" from file "
	<<tmpfilename<<endl;    
  }
  cout<<"After initializing my new MPS, exc, before imposing gauge"<<endl;
  exc->gaugeCond('R',1);
  cout<<"After imposing gauge R"<<endl;
  exc->gaugeCond('L',1);
  cout<<"After imposing gauge L"<<endl;
  if(lval==0){
    if(!tmpRes&&!initmpsdir.empty()&&initL<L){
      // In this case, I might try to optimize the added sites first, see if it helps
      contractor.sweepPart(initL/2,initL/2+L-initL,hamil,D,*exc);
    }
    //    contractor.findGroundState(hamil,D,&lambda,*exc);
    contractor.findGroundState(hamil,D,&lambda,*exc,offset,1,tmpfilename,freqSv);
    cout<<"Ground state found with eigenvalue "<<lambda<<endl;
    *out<<"% L="<<L<<endl;
    *out<<"% Delta="<<Delta<<endl;
    *out<<"% Ze="<<Z<<endl;
    *out<<"% a="<<a<<endl;
    *out<<"% M="<<M<<endl;
    *out<<"% D="<<D<<endl;
    *out<<"% Npen="<<Npen<<endl;
    *out<<"% Ntot="<<Ntot<<endl;
    *out<<"% offset="<<offset<<endl;
    *out<<"% E0="<<lambda-offset<<endl;
    cout<<"Now header was written to file "<<outfile<<endl;
  }
  else{
    if(!tmpRes&&!initmpsdir.empty()&&initL<L){
      contractor.sweepPartNextExcitedState(initL/2,initL/2+L-initL,hamil,D,levels,*exc);
    }
    cout<<"About to start findNextExcitedState"<<endl;
    contractor.findNextExcitedState(hamil,D,levels,&lambda,*exc,offset,1,tmpfilename,maxTime);
    // if(lambda>0){
    //   cout<<"Stopping iteration as I got an (excited) energy above zero!"<<endl;
    //   exit(212);
    // }
  }
  cout<<"computing expectation values"<<endl;
  Nval=contractor.contract(*exc,Nmpo,*exc);
  cout<<"First only Nval"<<endl;
  bool good=Npen==0||(abs(real(Nval)-Ntot)<=1E-4);
  // Save the state to disk. But watch out, bus errors??
  if(good){
    const string filename=mpsfilename(L,D,M,Z,lval,mpsdir);
    cout<<"About to write MPS to file "<<filename<<endl;
    exc->exportMPS(filename.data());
    cout<<"computing remaining expectation values"<<endl;
    Hval=contractor.contract(*exc,hamil,*exc);
    cout<<"computed Hval"<<endl;
    //    H2val=contractor.contract2(hamil,*exc);
    // cout<<"computed H2val"<<endl;
    normN=contractor.contract(*exc,*exc);
    cout<<"computed norm"<<endl;
    Sval=contractor.contract(*exc,Smpo,*exc);
    cout<<"computed Sval; now write to file"<<endl;
    // write out results to the file
    *out<<lval<<"\t"<<exc->getBond()<<"\t";
    *out<<real(normN)<<"\t";
    //*out<<real(H2val)<<"\t"; //<<imag(H2)<<"\t";
    *out<<0<<"\t"; //<<imag(H2)<<"\t";
    *out<<real(Hval)<<"\t"; //<<imag(Hval)<<"\t";
    *out<<real(Nval)<<"\t";
    *out<<real(Sval)<<"\t";
    *out<<endl;
    cout<<"Expectation values have been written to file"<<endl;
  }
  else{
    cout<<"Error: Level "<<lval<<" has been found with <N>="<<Nval
	<<", when desired was "<<Ntot<<": increasing penalty value "
	<<"(currently "<<Npen<<")"<<endl;
    if(file_exists(tmpfilename))
      remove(tmpfilename.data());
    Npen=10*Npen;lval--;
  }

  
  // Write the command to a file, if it needs to be relaunched
  const string jobfile=jobfilename(L,D,M,Z,jobsdir);
  cout<<"Preparing jobfile: "<<jobfile<<endl;
  stringstream commandstr;
  if(lval<nLev)
    commandstr<<argv[0];
  else //lval>=nLev => Launch the observable computation
    commandstr<<"./scsp2";
  commandstr<<" "<<argv[1];
  commandstr<<" -L="<<L<<" -M="<<M<<" -Delta="<<Delta;
  commandstr<<" -Ze="<<Z<<" -D="<<D<<" -a="<<a;
  commandstr<<" -append="<<app<<" -outfile="<<props.getProperty("outfile");
  commandstr<<" -outputdir="<<directory<<" -tol="<<tol;
  commandstr<<" -penaltyNe="<<Npen<<" -Ne="<<Ntot;
  commandstr<<" -penaltyS="<<Spen;
  commandstr<<" -offset="<<offset;
  commandstr<<" -mpsdir="<<mpsdir;
  if(!initmpsdir.empty())
    commandstr<<" -initmpsdir="<<initmpsdir; // what I just computed is the new initmpsdir
  if(initD>0){
    commandstr<<" -initD="<<initD<<" ";
    commandstr<<" -noise="<<noise;
  }
  commandstr<<" -jobsdir="<<jobsdir;
  commandstr<<" -maxTime="<<maxTime;
  commandstr<<" -nLev="<<nLev;
  //  int mem=ceil(L*D*D*(4*16.+4*4.*(lval+1))*16E-6)+200;
  int mem=ceil((L*D*D*d2*d2*(lval+1)+2*L*D*D*Dh+d2*d2*Dh*Dh*L+2*d2*d2*Dh*Dh*D*D)*16E-6)+200;
  commandstr<<" -msubmem="<<mem<<"M";
  
  cout<<"Line "<<commandstr<<" will be written to jobfile: "<<jobfile<<endl;

  ofstream* ojob=new ofstream(jobfile.data());
  if(!ojob->is_open()){
    cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
    cout<<commandstr.str()<<endl;
    cout<<endl;
    // Do not exit, because I still need to save the MPS
  }
  else{
    *ojob<<"#!/bin/bash"<<endl;
    *ojob<<"msub_modified_tqo097 -N sc."<<L<<"."<<Z<<"."<<M<<"."<<D<<" -l h_vmem="
	 <<mem<<"M -raw ";
    *ojob<<commandstr.str()<<endl;
  }

  // // Save the state to disk. But watch out, bus errors??
  // if(good){
  //   const string filename=mpsfilename(L,D,M,Z,lval,mpsdir);
  //   cout<<"About to write MPS to file "<<filename<<endl;
  //   exc->exportMPS(filename.data());
  // }

  // Now, if it was saved, I can delete the temporary file
  if(file_exists(tmpfilename)){
    remove(tmpfilename.data());
  }
  
  if(ojob->is_open()) ojob->close();
  delete ojob;
  out->close();
  delete out;

  // Clear the MPS vector
  while(levels.size()>1){
    delete levels.back();
    levels.pop_back();
  }  


   // for(int z=1;z<L;z++){
   //   delete Skmpos[z];
   // }
   // delete []Skmpos;
}

#include <sstream>    

const string mpsfilename(int L,int D,int M,int Ze,int level,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS_L"<<L<<"_M"<<M<<"_Z"<<Ze<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

const string jobfilename(int L,int D,int M,int Ze,const string jobsdir){
  stringstream s;
  s<<jobsdir<<"/job_L"<<L<<"_M"<<M<<"_Z"<<Ze<<"_D"<<D;
  //if(tol!=1E-7) s<<"_tol"<<-log10(tol);
  s<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

