
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"
//#include "HermitianTensorMultiplier.h"

#include "SchwingerHamiltonianSz.h"
#include "SpinMPO.h"
#include <cmath>

using namespace std;
using namespace shrt;

#define MAXLEN 120

/** 
Generates a file name for the MPS corresponding to a certain excitation.
*/
const string mpsfilename(int L,int D,double x,double mg,int level,const string mpsdir);
/** 
Generate the name of the file with the command in it
*/
const string jobfilename(int L,int D,double x,double mg,double tol,const string mpsdir);

/** Check if the filename exists */
bool file_exists(const string filename);


/** As in phase1_save, but compute just a single vector (After
    checking the already computed states). When the scalar is found
    and saved, a special file is created in the destination directory,
    so that future instances don't waste time.

    This program tries to find the vector and scalar mass gap by first
    finding the MPS approximation to the ground state, then finding
    another one for the first excited state, but in every step, the
    Hamiltonian is projected onto the Sz=0 subspace.  A certain number
    of (MPS approximated) eigenstates is computed, and for each of
    them the expectation value of the operator SR (translation of one
    to the right times sigmax rotation) is calculated, together with
    the momentum operator.  When the scalar candidate is located
    (value of SR, and- maybe- momentum), the program stops.  Here we
    save all the MPS computed to disk, to a directory whose name is
    uniquely related to the instance , so that at the beginning, the
    directory (mpsdir in properties file) is checked, and if MPS files
    exst, they are read and the program starts using them already.

    Instead of targetting a fixed number of eigenstates, it stops
    right after identifying the candidate scalar, via the sign
    of the real part of <SR> (>0) and its modulus (close to that of
    the GS) and the value (small) of the momentum operator P2.

    Receives as argument a Properties file, such as
    config/excitProp.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.

    If property onlyvec is set to 1, the program stops after finding
    the vector state (default is to compute until the scalar).

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

  // Store the command line, for this is what should be relaunched if not all levels are done
  stringstream commandstr;
  for(int c=0;c<argc;c++) commandstr<<argv[c]<<" ";
  // Now, str contains the command line, that I need to relaunch, if this was not the last run

  int L = props.getIntProperty("L");
  double mg=props.getDoubleProperty("mg");
  double alpha=props.getDoubleProperty("alpha");
  double x=props.getDoubleProperty("x");
  //  int D0 = props.getIntProperty("D0");
  //  int nLev = props.getIntProperty("nLev");
  directory=props.getProperty("outputdir");
  const string outfile=directory+"/"+props.getProperty("outfile");
  int app=props.getIntProperty("append");
  if(app<0) app=1; // default will be to append
  double tol=props.getDoubleProperty("tol");
  // int unif=1; // All same D 
  // int unif=props.getIntProperty("unif");
  int D=props.getIntProperty("D");
  double zpen=props.getDoubleProperty("penalty");
  double offset=props.getDoubleProperty("offset");
  string mpsdir=props.getProperty("mpsdir"); // default current
  if(mpsdir.empty()) 
    mpsdir=".";
  int stopAtVec=props.getIntProperty("onlyvec");
  if(stopAtVec<0) stopAtVec=0; // default searches also scalar
  int stopAtGS=props.getIntProperty("onlygs");
  if(stopAtGS<0) stopAtGS=0; // default searches also vector, at least
  else stopAtVec=1; // just in case, I stop after 1
  int minNumber=props.getIntProperty("minnumber");
  bool runUntil=(minNumber>0);  
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

  if(stopAtGS){ // If only computing GS, Iwill then use larger D
    commandstr.str(""); // empty again
    commandstr<<argv[0];
    commandstr<<" "<<argv[1];
    commandstr<<" -L="<<L<<" -mg="<<mg<<" -alpha="<<alpha;
    commandstr<<" -x="<<x<<" -D="<<D+20;
    commandstr<<" -append="<<app<<" -outfile=L"<<L<<"x"<<x<<"D"<<D+20;
    commandstr<<" -outputdir="<<directory<<" -tol="<<tol;
    commandstr<<" -penalty="<<zpen;
    commandstr<<" -offset="<<offset;
    commandstr<<" -mpsdir="<<mpsdir;
    commandstr<<" -initmpsdir="<<mpsdir; // what I just computed is the new initmpsdir
    commandstr<<" -initD="<<D<<" ";
    commandstr<<" -noise="<<noise;
    commandstr<<" -jobsdir="<<jobsdir;
    commandstr<<" -maxTime="<<maxTime;
    commandstr<<" -minnumber=0 -onlygs=1 ";
    int mem=ceil(L*(D+20)*(D+20)*(16.)*16E-6)+200;
    commandstr<<" -msubmem="<<mem<<"M";
  }


  // Parameters I do not use here: chemical potential and weight of the gauge term
  double nu=0;
  // I might want to switch off the gauge terms
  //double gweight=props.getDoubleProperty("gweight");
  //if(gweight==-1) gweight=1; // default: do not alter
  double gweight=1.; // Not touching the gauge terms (std)

  double mu=2*mg*sqrt(x);
  double L0=0.; // l0 parameter, which is irrelevant, as only appears
		// as alpha+l0

  //  double t0=.01; // todo: differently -> reference time

  cout<<"Initialized arguments: L="<<L
      <<", mu="<<mu
      <<", x="<<x
      <<", alpha="<<alpha
      <<", outfile="<<outfile
      <<", app="<<app
      <<", tol="<<tol
      <<", zpen="<<zpen
      <<", gweight="<<gweight
      <<", D0="<<D<<endl;

#ifdef DISKSTRG
  //  int errD=system("rm -rf /ptmp/mpq/banulsm/*");
  char tmpdir[150];
  sprintf(tmpdir,"/ptmp/mpq/banulsm/Schwinger/Sitesx%dL%dD%dmg%dXXXXXX",x,L,D,(int)(mg*1000));
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
  int d=2;
  // First: put H a
  int Dop=5; // bond dimension of Schwinger H
  double lambda=0.;
  SchwingerHamiltonianSz hSch(L,mu,x,L0,alpha,zpen,offset,nu,gweight);
  //  SchwingerHamiltonian hSch(L,mu,x,L0,alpha,offset,nu,gweight);
  const MPO& hamil=hSch.getHMPO();
  // Observables I will compute
  // Compute the expectation value of Sz, as it commutes with H
  MPO Szmpo(L);
  SpinMPO::getSzMPO(L,d,Szmpo);
  // Sz2 is the penalty term (we need to substract it) 
  MPO Sz2mpo(L);
  SpinMPO::getSz2MPO(L,d,Sz2mpo);
  // The momentum operator <P2> is what we want
  MPO Pmpo(L);
  hSch.constructMomentumMPO(Pmpo);
  const MPO* oprsP[2]={&Pmpo,&Pmpo};
  MPO P2mpo(L);
  MPO::join(2,oprsP,P2mpo);
  // // The momentum squared operator
  // MPO PSqmpo(L);
  // hSch.constructMomentumSquaredMPO(PSqmpo);

  // the shift by one and rotate sigmax
  MPO SRmpo(L);
  hSch.constructShiftRotateMPO(SRmpo);
  //  SRmpo.exportForMatlab("SRmpo.m");

  // The fermionic condensate and the Gamma_alpha, Gamma_5 parameters
  MPO condMPO(L),GammaA(L),Gamma5(L);
  hSch.constructCondensateMPO(condMPO);
  hSch.constructGammaAlphaMPO(GammaA);
  hSch.constructGamma5MPO(Gamma5);

  vector<MPS*> levels;

  // First, check how many MPS I already have in mpsdir
  bool oldMPS=true;
  int lval=0;
  MPS* nextLev;
  while(oldMPS&&(!stopAtVec||lval<=1)){
    const string filename=mpsfilename(L,D,x,mg,lval,mpsdir);
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

  // If I was instructed to stop after vector and it is already there,
  // something strange must have happened. I break now.
  if(stopAtVec&&lval>1){
    cout<<"Stopping at vector, which is already computed. Nothing to do"<<endl;
    out->close();
    delete out;
    exit(333);
  }

  if(stopAtGS&&lval>0){
    cout<<"Stopping at GS, which is already computed. Nothing to do"<<endl;
    out->close();
    delete out;
    exit(333);
  }


  // Now I already have read all existing MPS
  cout<<"Recovered "<<lval<<" already computed MPS levels from disk. Computing level "<<lval<<" now."<<endl;


  // values for each level
  complex_t Hval,H2val,normN,Szval,Sz2val,P2val,SRval,fC,gammaA,gamma5;
  bool found=0;

  MPS* exc=new MPS(L,D,d);
  bool tmpRes=false;
  // First check if there is any temporary result
  const string tmpfilename=mpsfilename(L,D,x,mg,lval,mpsdir)+".tmp";
  if(!file_exists(tmpfilename)){
    const string initfilename=mpsfilename(initL,initD,x,mg,lval,initmpsdir);
    if(initmpsdir.empty())
      exc->setRandomState();
    else{
      cout<<"Will search initial state for level "<<lval<<" in "<<initfilename<<endl;
      if(!file_exists(initfilename)){
	cout<<"ERROR: Initial MPS not found in "<<initfilename<<endl;
	exit(1);
      }
      if(initL<L){
	MPS auxMPS(initL,initD,d);
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
    contractor.findGroundState(hamil,D,&lambda,*exc);
    cout<<"Ground state found with eigenvalue "<<lambda<<endl;
    *out<<"% L="<<L<<endl;
    *out<<"% mg="<<mg<<endl;
    *out<<"% x="<<x<<endl;
    *out<<"% alpha="<<alpha<<endl;
    *out<<"% D="<<D<<endl;
    *out<<"% penalty="<<zpen<<endl;
    *out<<"% offset="<<offset<<endl;
    *out<<"% gweight="<<gweight<<endl;
    *out<<"% E0="<<lambda-offset<<endl;
  }
  else{
    if(!tmpRes&&!initmpsdir.empty()&&initL<L){
      contractor.sweepPartNextExcitedState(initL/2,initL/2+L-initL,hamil,D,levels,*exc);
    }
    cout<<"About to start findNextExcitedState"<<endl;
    contractor.findNextExcitedState(hamil,D,levels,&lambda,*exc,0.,1,tmpfilename,maxTime);
    if(lambda>0){
      cout<<"Stopping iteration as I got an (excited) energy above zero!"<<endl;
      exit(212);
    }
  }
  
  Hval=contractor.contract(*exc,hamil,*exc);
  H2val=contractor.contract2(hamil,*exc);
  normN=contractor.contract(*exc,*exc);
  Szval=contractor.contract(*exc,Szmpo,*exc);
  Sz2val=contractor.contract(*exc,Sz2mpo,*exc);
  P2val=contractor.contract(*exc,P2mpo,*exc);
  SRval=contractor.contract(*exc,SRmpo,*exc);
  fC=contractor.contract(*exc,condMPO,*exc); // fermion condensate
  gammaA=contractor.contract(*exc,GammaA,*exc); // Gamma_alpha
  gamma5=contractor.contract(*exc,Gamma5,*exc); // Gamma_5
  bool good=(abs(Szval)<=1E-5);
  if(good){
    // write out results to the file
    *out<<lval<<"\t"<<exc->getBond()<<"\t";
    *out<<real(normN)<<"\t";
    *out<<real(H2val)<<"\t"; //<<imag(H2)<<"\t";
    *out<<real(Hval-zpen*Sz2val)-offset<<"\t"; //<<imag(Hval)<<"\t";
    *out<<real(Szval)<<"\t";
    *out<<real(P2val)<<"\t";
    *out<<real(SRval)<<"\t"<<imag(SRval)<<"\t";
    *out<<real(fC)<<"\t";
    *out<<real(gammaA)<<"\t";
    *out<<real(gamma5)<<"\t";
    *out<<endl;
    // Decide if this was the scalar and I can thus stop.
    if(lval>0) // GS does not count
      if(real(SRval)>0){
	found=true;
	double Es=real(Hval-zpen*Sz2val)-offset;
	cout<<"** FIRST Candidate found, at l="<<lval<<" with properties: "<<endl; 
	cout<<"** E="<<Es<<", P2="<<P2val<<", SR="<<SRval<<endl;
	if(abs(phase(SRval))>M_PIl/4.){
	  cout<<"WARNING!! THis might not be it!!!"<<endl;
	}
      }
  }
  else{
    cout<<"Error: Level "<<lval<<" has been found with <Sz>="<<Szval<<": increasing penalty value "
	<<"(currently "<<zpen<<")"<<endl;
    if(file_exists(tmpfilename))
      remove(tmpfilename.data());
    zpen=10*zpen;lval--;
  }

  
  // Write the command to a file, if it needs to be relaunched
  if(runUntil){
    const string jobfile=jobfilename(L,D,x,mg,tol,jobsdir);
    stringstream commandstr;
    if(lval<minNumber){
      commandstr<<argv[0];
    }
    else{ // Now launch searchOrder
      commandstr<<"./searchO";
    }
    commandstr<<" "<<argv[1];
    commandstr<<" -L="<<L<<" -mg="<<mg<<" -alpha="<<alpha;
    commandstr<<" -x="<<x<<" -D="<<D;
    commandstr<<" -append="<<app<<" -outfile=L"<<L<<"x"<<x<<"D"<<D;
    commandstr<<" -outputdir="<<directory<<" -tol="<<tol;
    commandstr<<" -penalty="<<zpen;
    commandstr<<" -offset="<<offset;
    commandstr<<" -mpsdir="<<mpsdir;
    commandstr<<" -initmpsdir="<<initmpsdir; // what I just computed is the new initmpsdir
    commandstr<<" -initD="<<initD<<" ";
    commandstr<<" -noise="<<noise;
    commandstr<<" -jobsdir="<<jobsdir;
    commandstr<<" -maxTime="<<maxTime;
    commandstr<<" -minnumber="<<minNumber+1;
    int mem=ceil(L*D*D*(16.+4.*minNumber)*16E-6)+200;
    commandstr<<" -msubmem="<<mem<<"M";
    ofstream* ojob=new ofstream(jobfile.data());
    if(!ojob->is_open()){
      cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
      cout<<commandstr.str()<<endl;
      cout<<endl;
      // Do not exit, because I still need to save the MPS
    }
    else{
      *ojob<<"#!/bin/bash"<<endl;
      if(tol!=1E-7){
	*ojob<<"msub_modified_tqo097 -N tol"<<(int)(-log10(tol))<<".D"<<D<<".M"<<props.getProperty("mg")<<"."<<x<<"."<<L<<" -l h_vmem="
	     <<mem<<"M -raw ";
      }
      else
      *ojob<<"msub_modified_tqo097 -N D"<<D<<".M"<<props.getProperty("mg")<<"."<<x<<"."<<L<<" -l h_vmem="
	   <<mem<<"M -raw ";
      *ojob<<commandstr.str()<<endl;
      ojob->close();
    }
    delete ojob;
  }
  else{
    int nextD=stopAtGS?D+20:D;
    if((stopAtGS&&D<=140&&tol==1E-7)||(!stopAtGS&&!found&&(!stopAtVec||lval<1))){
      const string jobfile=jobfilename(L,nextD,x,mg,tol,jobsdir);
      ofstream* ojob=new ofstream(jobfile.data());
      if(!ojob->is_open()){
	cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
	cout<<commandstr.str()<<endl;
	cout<<endl;
	// Do not exit, because I still need to save the MPS
      }
      else{
	int nextLval=stopAtGS?0:lval;
	int mem=ceil(L*nextD*nextD*(16.+4.*(nextLval+1))*16E-6)+200;
	if(tol!=1E-7){
	  *ojob<<"msub_modified_tqo097 -N tol"<<(int)(-log10(tol))<<".D"<<nextD<<".M"<<props.getProperty("mg")<<"."<<x<<"."<<L<<" -l h_vmem="
	       <<mem<<"M -raw ";
	}
	else
	  *ojob<<"msub_modified_tqo097 -N D"<<nextD<<".M"<<props.getProperty("mg")<<"."<<x<<"."<<L<<" -l h_vmem="
	       <<mem<<"M -raw ";
	*ojob<<commandstr.str()<<endl;
	ojob->close();
      }
      delete ojob;
    }
  }
  // Save the state to disk. But watch out, bus errors??
  if(good){
    const string filename=mpsfilename(L,D,x,mg,lval,mpsdir);
    exc->exportMPS(filename.data());
  }

  // Now, if it was saved, I can delete the temporary file
  if(file_exists(tmpfilename)){
    remove(tmpfilename.data());
  }
  
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

const string mpsfilename(int L,int D,double x,double mg,int level,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS_L"<<L<<"_m"<<mg<<"_x"<<x<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

const string jobfilename(int L,int D,double x,double mg,double tol,const string jobsdir){
  stringstream s;
  s<<jobsdir<<"/job_L"<<L<<"_m"<<mg<<"_x"<<x<<"_D"<<D;
  if(tol!=1E-7) s<<"_tol"<<-log10(tol);
  s<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

bool file_exists(const string filename)
{
  ifstream ifile(filename.data());
  //cout<<"file_exists("<<filename<<")="<<ifile<<endl;
  return ifile;
}
