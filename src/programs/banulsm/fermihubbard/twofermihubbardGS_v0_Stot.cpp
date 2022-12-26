
#include <math.h>
#include <iomanip>

#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "misc.h"
#include <list>

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "Properties.h"
#include "TwoOrbitalFermiHubbardHamiltonian_v0.h"

using namespace shrt;

/** twofermihubbardGS_v0 runs the findGroundState routine with the MPO for the
    Hamiltonian \ref <TwoOrbitalFermiHubbardHamiltonian_v0>, and computes GS and
    maybe other things

    In this version (20.12.2022), the initial state is chosen to be
    another one with larger Sz, and a S^- operator is applied (or
    more) before normalizing and searching variationally. It gets an
    initial Sz, searches for one (or more less) and launches the next.
*/

int findTmpFile(int D,const string& mpsdir,const string& baseName,string& mpsfile);
const string mpsFileName(const string& mpsdir,const string& baseName,int D);
const string getMPSfilename(int L,const vector<double>& t,const vector<double>& U,double V,double Vex,const vector<double>& mu_g,const vector<double>& mu_e,bool GS=1);

const string getMPSfilenameConstraints(int L,double penNg,int Ng,double penNe,int Ne,double penSz,int Sz,double penNtot,int Ntot,const vector<double>& t,const vector<double>& U,double V,double Vex,const vector<double>& mu_g,const vector<double>& mu_e,bool GS=1);

bool checkQuantumNumber(const MPS& gs,const MPO& mpo,int targetV,double tolVal,double tolVar);


void prepareStateWithNfermions(MPS& mps,int Nf,int L,int D);

/** 
Generate the name of the file with the command in it
*/

const string jobFileName(const string& jobsdir,const string& baseName,int D);
void prepareNewJob(const string& jobsdir,const string& basempsfile,int nextSz,double nextpenNtot,double nextpenNg,double nextpenNe,double nextpenSz,
		   const string& exec,const string& propFile,const Properties& props);

void prepareNewJobNextD(const string& jobsdir,const string& basempsfile,int nextD,double nextpenNtot,double nextpenNg,double nextpenNe,double nextpenSz,
		   const string& exec,const string& propFile,const Properties& props);


int main(int argc,const char* argv[]){
  // Read input arguments
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
  double tg0=props.getDoubleProperty("tg0");
  double tg1=props.getDoubleProperty("tg1");
  double te0=props.getDoubleProperty("te0");
  double te1=props.getDoubleProperty("te1");
  double Ug=props.getDoubleProperty("Ug");
  double Ue=props.getDoubleProperty("Ue");
  double V=props.getDoubleProperty("V");
  double Vex=props.getDoubleProperty("Vex");
  double mu_g0=props.getDoubleProperty("mu_g0");
  double mu_g1=props.getDoubleProperty("mu_g1");
  double mu_e0=props.getDoubleProperty("mu_e0");
  double mu_e1=props.getDoubleProperty("mu_e1");
  double penNtot=props.getDoubleProperty("penaltyNtot");
  if(penNtot<0) penNtot=0;
  int Ntot=props.getIntProperty("Ntot");
  if(Ntot<0) Ntot=0;
  double penNg=props.getDoubleProperty("penaltyNg");
  if(penNg<0) penNg=0;
  int Ng=props.getIntProperty("Ng");
  if(Ng<0) Ng=0;
  double penNe=props.getDoubleProperty("penaltyNe");
  if(penNe<0) penNe=0;
  int Ne=props.getIntProperty("Ne");
  if(Ne<0) Ne=0;
  double penSz=props.getDoubleProperty("penaltySz");
  if(penSz<0){
    cout<<"ERROR: This version of twoFH GS needs a target Sz value, with a penalty"<<endl;
    exit(1);
  }
  int Sz=props.getIntProperty("Sz");

  if(penNe>0&&penNg>0) penNtot=0;
  
  int D=props.getIntProperty("D");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir"); // where are the MPS
  int Stot=props.getIntProperty("Stot");
  if(Stot<0) Stot=max(Ntot,Ne+Ng)/2;
  int Szprev=props.getIntProperty("Szprev");
  int incrSz=Sz-Szprev; // if negative, going down with Sminus, o.w. up with Splus (recall factor 2) 
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  int freqSv=props.getIntProperty("freqSv"); // save to file every so many sweeps
  if(freqSv<0) freqSv=10;
  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  // To prepare new jobs that rerun with increased D (or, if necessary, increased penalty)
  bool relaunch=props.getIntProperty("relaunch")==1;
  int incrD=props.getIntProperty("incrD");
  int maxD=props.getIntProperty("maxD");
  int minD=props.getIntProperty("minD"); // the one used to sweep the first time
  string jobsdir;
  if(relaunch){
    jobsdir=props.getProperty("jobsdir"); 
  }
  //  bool refine=(props.getIntProperty("refine")==1); // in that case, only runs if the same state (same D, no tmp) exists. the same is passed on
  // bool restart=(props.getIntProperty("restart")==1); // in that case, only runs if the state with D0 exists, but not passed on to new job
  // int Drestart=D;
  // if(restart){
  //   Drestart=props.getIntProperty("Drestart");
  //   if(Drestart<=0) restart=0;
  // }

  //refine=0;//restart=0;tryVs=0; // in this version, I don't do this
  // if(tryVs&&restart){
  //   cout<<"Option tryVs is incompatible with restart"<<endl;
  //   exit(1);
  // }  
  srandom(time(NULL));
  // cout<<"Initialized arguments: L="<<L
  //     <<", t="<<t
  //     <<", U="<<U
  //     <<", outfile="<<outfname
  //     <<", app="<<app
  //     <<", D="<<D<<endl;

  // To ease comparisons with matlab, for debugging
  cout<<"Parameters:"<<endl;
  cout<<"tg0="<<tg0<<";tg1="<<tg1<<";te0="<<te0<<";te1="<<te1<<";"<<endl;
  cout<<"Ug="<<Ug<<";Ue="<<Ue<<";"<<endl;
  cout<<"V="<<V<<";Vex="<<Vex<<";"<<endl;
  cout<<"mu_g0="<<mu_g0<<";mu_g1="<<mu_g1<<";mu_e0="<<mu_e0<<";mu_e1="<<mu_e1<<";"<<endl;
  if(penNtot>0)
    cout<<"targetNtot="<<Ntot<<"; lambda="<<penNtot<<";"<<endl;
  if(penNg>0) cout<<"targetNg="<<Ng<<"; lambda="<<penNg<<";"<<endl;
  if(penNe>0) cout<<"targetNe="<<Ne<<"; lambda="<<penNe<<";"<<endl;
  cout<<"targetSz="<<Sz<<"; lambda="<<penSz<<";"<<endl;
  cout<<"origSz="<<Szprev<<";"<<endl;
  cout<<"L="<<L<<";"<<endl;
  cout<<"[H,Nfg,Nfe,Sz]=FermiHubbardHamiltonianTwoOrbital(L,tg0,tg1,te0,te1,Ug,Ue,V,Vex,mu_g0,mu_g1,mu_e0,mu_e1);\nvalsEx=eig(H);\nmin(valsEx)"<<endl;
  
  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);
  //  cout<<"Initialized Contractor"<<endl;
  int d=2;
  //  int dF=d*d; // dimension of fermion sites (L+1)

  MPS gs(4*L,D,d);
  // First, check whether there is a tmp file

  // Prepare parameters
  vector<double> t;t.push_back(tg0);t.push_back(tg1);t.push_back(te0);t.push_back(te1);
  vector<double> U;U.push_back(Ug);U.push_back(Ue);
  vector<double> mu_g;mu_g.push_back(mu_g0);mu_g.push_back(mu_g1);
  vector<double> mu_e;mu_e.push_back(mu_e0);mu_e.push_back(mu_e1);

  TwoOrbitalFermiHubbardHamiltonian ham(L,t,U,V,Vex,mu_g,mu_e,Ntot,penNtot);
  ham.setPenalties(Ng,penNg,Ne,penNe,Sz,penSz,Ntot,penNtot);
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=ham.getHMPO();

  // Get also the operators
  MPO fermNr(4*L);
  ham.getFermionNrMPO(fermNr); // mpo with total nr of fermions
  MPO fermNr0(4*L);
  ham.getFermionNrOrbitalMPO(fermNr0,true); // mpo with total nr of fermions in g
  MPO fermNr1(4*L);
  ham.getFermionNrOrbitalMPO(fermNr1,false); // mpo with total nr of fermions in e
  MPO fermSz(4*L);
  ham.getTotalSzMPO(fermSz);
  MPO fermSx(4*L);
  ham.getTotalSxMPO(fermSx);
  MPO fermSy(4*L);
  ham.getTotalSyMPO(fermSy);
  cout<<"Constructed the hamil MPO and Nferm MPOs"<<endl;


  
  string basempsfile=getMPSfilenameConstraints(L,penNg,Ng,penNe,Ne,penSz,Sz,penNtot,Ntot,t,U,V,Vex,mu_g,mu_e,1); // basename of MPS files
  string basempsfileRef=getMPSfilenameConstraints(L,penNg,Ng,penNe,Ne,penSz,Szprev,penNtot,Ntot,t,U,V,Vex,mu_g,mu_e,1); // basename of MPS files
  string mpsfile=mpsFileName(mpsdir,basempsfile,D);
  string mpsfileRef=mpsFileName(mpsdir,basempsfileRef,D);
  string tmpMPSfile;

  // First check if it is already computed, or if there is a tmp file to continue
  if(file_exists(mpsfile)){
    cout<<"Found already state at "<<mpsfile<<endl;
    // Check if the S^2 is ok
    gs.importMPS(mpsfile.data());
    complex_t valSx2=contractor.contract2(fermSx,gs);
    complex_t valSy2=contractor.contract2(fermSy,gs);
    complex_t valSz2=contractor.contract2(fermSz,gs);
    double valS2=real(valSx2+valSy2+.25*valSz2);
    if(abs(valS2-Stot*(Stot+1))<.2*Stot*(Stot+1)){
      if(incrD<=0||maxD<=0||D==minD){
	cout<<"Launching for another Sz"<<endl;
	string basempsfileNext=getMPSfilenameConstraints(L,penNg,Ng,penNe,Ne,penSz,Sz+incrSz,penNtot,Ntot,t,U,V,Vex,mu_g,mu_e,1); // basename of MPS files
	prepareNewJob(jobsdir,basempsfileNext,Sz+incrSz,penNtot,penNg,penNe,penSz,argv[0],argv[1],props);	
      }
      if(incrD>0&&maxD>0&&D+incrD<=maxD){
	cout<<"Preparing (also) a larger D run"<<endl;
	prepareNewJobNextD(jobsdir,basempsfile,D+incrD,penNtot,penNg,penNe,penSz,argv[0],argv[1],props);	
      }
      exit(0);
    }
    cout<<"But the state does not have the desired S2 (found "<<valS2<<" wanted "<<Stot*(Stot+1)<<")"<<endl;
  }
  // check for other files
  bool foundFile=false;
  int tmpD=findTmpFile(D,mpsdir,basempsfile,tmpMPSfile);
  if(tmpD>0){
    cout<<"Recovered initial guess from file "<<tmpMPSfile<<endl;
    gs.importMPS(tmpMPSfile.data());
    if(gs.getLength()!=4*L){
      cout<<"Error! File "<<tmpMPSfile<<" does not contain valid MPS "<<endl;
      exit(1);
    }
    if(tmpD==D){// found some tmp with same D
      foundFile=true;
    }
    else{ // found some smaller D file
      complex_t valSx2=contractor.contract2(fermSx,gs);
      complex_t valSy2=contractor.contract2(fermSy,gs);
      complex_t valSz2=contractor.contract2(fermSz,gs);
      double valS2=real(valSx2+valSy2+.25*valSz2);
      if(abs(valS2-Stot*(Stot+1))<.2*Stot*(Stot+1)){
	foundFile=true;
	cout<<"Found smaller D with adequate S^2 "<<valS2<<endl;
      }
      else{
	cout<<"The lower D has S^2="<<abs(valS2)<<" (targetting "<<Stot*(Stot+1)<<")"<<endl;
      }
    }
  }
  if(!foundFile){ // no tmp file, neither a valid smaller D
    cout<<"Failed to recover a tmp file ("<<basempsfile<<")=> trying to use previous Sz"<<endl;  
  
    if(file_exists(mpsfileRef)){
      cout<<"Found state for reference Sz "<<mpsfileRef<<endl;
      gs.importMPS(mpsfileRef.data());
      // Now apply Sminus or plus a certain number of times
      MPO Sladder(4*L);
      
      if(incrSz>0)
	ham.getTotalSplusMPO(Sladder);
      else
	ham.getTotalSminusMPO(Sladder);

      int nrLad=abs((Szprev-Sz)/2);
      cout<<"Will need to apply a ladder operator "<<nrLad<<" times to go from "<<Szprev<<" to "<<Sz<<endl;
      for(int kl=0;kl<nrLad;kl++){
	MPS aux(gs);
	contractor.optimize(Sladder,aux,gs,D);
	// normalize
	gs.gaugeCond('R',1);
	gs.gaugeCond('L',1);
	// Check energy and expectation values
	cout<<"After "<<kl+1<<" application(s) of the ladder operator, the state has: "<<endl;
	complex_t valH=contractor.contract(gs,hamil,gs);
	cout<<" <H>="<<real(valH);
	complex_t valNf=contractor.contract(gs,fermNr,gs);
	cout<<"\t and nrs of fermions:"<<endl;
	cout<<"\tNtot "<<valNf<<endl;
	complex_t valNfg=contractor.contract(gs,fermNr0,gs);
	cout<<"\tNg "<<valNfg<<endl;
	complex_t valNfe=contractor.contract(gs,fermNr1,gs);
	cout<<"\tNe "<<valNfe<<endl;
	complex_t valSz=contractor.contract(gs,fermSz,gs);
	complex_t valSz2=contractor.contract2(fermSz,gs);
	cout<<"\tSz "<<valSz<<" variance:"<<real(valSz2-valSz*valSz)<<endl;
	complex_t valSx2=contractor.contract2(fermSx,gs);
	complex_t valSy2=contractor.contract2(fermSy,gs);
	cout<<"\tS^2 "<<valSx2+valSy2+.25*valSz2<<endl;
      }      
    }
    else{
      cout<<"ERROR: Failed to recover a file for previous Sz ("<<basempsfileRef<<")"<<endl;
      exit(1);
    }
  }
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  //  cout<<"Initial state has norm "<<contractor.contract(gs,gs)<<endl;
  
  // if(4*L<=12){
  //   hamil.exportForMatlab("fH2Hmpo.m");
  //   fermNr.exportForMatlab("Nfmpo.m");
  //   fermNr0.exportForMatlab("Nf0mpo.m");
  //   fermNr1.exportForMatlab("Nf1mpo.m");
  //   fermSz.exportForMatlab("Szmpo.m");
  // }

  // Check if the initial state (read or random) is too far from the targetted numbers, and, if so, perturb it
  if(penNtot+penNg+penNe+penSz>0){
    cout<<"Checking initial state for targeted quantum numbers"<<endl;
    bool goodsector=1;
    if(penNg>0){
      cout<<"Checking Ng "<<endl;
      goodsector=goodsector&&checkQuantumNumber(gs,fermNr0,Ng,.5,1E-3);
    }
    if(penNe>0){
      cout<<"Checking Ne "<<endl;
      goodsector=goodsector&&checkQuantumNumber(gs,fermNr1,Ne,.5,1E-3);
    }
    if(penSz>0){
      cout<<"Checking Sz "<<endl;
      goodsector=goodsector&&checkQuantumNumber(gs,fermSz,Sz,.5,1E-3);
    }
    if(penNtot>0){
      cout<<"Checking Ntot "<<endl;
      goodsector=goodsector&&checkQuantumNumber(gs,fermNr,Ntot,.5,1E-3);
    }
    if(!goodsector){
      cout<<"Some nr failed=> perturbing the initial state"<<endl;
      gs.perturb(.1);
    }
  }
  else{
    cout<<"No need to check Q Nrs? penalties:"<<penNg<<","<<penNe<<","<<penSz<<","<<penNtot<<endl;
  }
  
  MPO doubleG(4*L);
  ham.getDoubleOccupancyOrbitalMPO(doubleG,1);
  MPO doubleE(4*L);
  ham.getDoubleOccupancyOrbitalMPO(doubleE,0);

  
  double E0=0.;
  E0=real(contractor.contract(gs,hamil,gs));
  cout<<"Initial state, with initial energy ";
  cout<<E0<<endl;
  cout<<"\t and initial nr of fermions "<<contractor.contract(gs,fermNr,gs)<<endl;

   
  tmpMPSfile=mpsfile+"_tmp";//int freqSv=10; // save to file every so many sweeps
  double offset=0.;int nrEig=1;
  // // sice I am running, I will move a posssible final file to tmp status, in case this run is interrupted!
  // if(file_exists(mpsfile)){
  //   rename(mpsfile.data(),tmpMPSfile.data());
  //   cout<<"Renamed existing "<<mpsfile<<" to "<<tmpMPSfile<<endl;
  // }
  contractor.findGroundState(hamil,D,&E0,gs,offset,nrEig,tmpMPSfile,freqSv);
  cout<<"Ground state found with eigenvalue "<<E0<<endl;
  complex_t valH2=contractor.contract2(hamil,gs);
  complex_t valNf=contractor.contract(gs,fermNr,gs);
  complex_t valNfg=contractor.contract(gs,fermNr0,gs);
  complex_t valNfe=contractor.contract(gs,fermNr1,gs);
  complex_t valSz=contractor.contract(gs,fermSz,gs);
  complex_t valNf2=contractor.contract2(fermNr,gs);
  complex_t valNfg2=contractor.contract2(fermNr0,gs);
  complex_t valNfe2=contractor.contract2(fermNr1,gs);
  complex_t valSz2=contractor.contract2(fermSz,gs);
  complex_t valSx2=contractor.contract2(fermSx,gs);
  complex_t valSy2=contractor.contract2(fermSy,gs);
  double valS2=real(valSx2+valSy2+.25*valSz2);
  cout<<" <H^2>="<<real(valH2);
  cout<<"\t and nrs of fermions:"<<endl;
  cout<<"\tNtot "<<valNf<<" variance:"<<real(valNf2-valNf*valNf)<<endl;
  cout<<"\tNg "<<valNfg<<" variance:"<<real(valNfg2-valNfg*valNfg)<<endl;
  cout<<"\tNe "<<valNfe<<" variance:"<<real(valNfe2-valNfe*valNfe)<<endl;
  cout<<"\tSz "<<valSz<<" variance:"<<real(valSz2-valSz*valSz)<<endl;
  cout<<"\tS^2 "<<valS2<<endl;
  // save it
  gs.exportMPS(mpsfile.data());
  complex_t valDg=contractor.contract(gs,doubleG,gs);
  complex_t valDe=contractor.contract(gs,doubleE,gs);

  // Write results to txt file
  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output (append="<<app<<")"<<endl;
      exit(1);
    }
    *out<<"%L\tD\t<H>\t<H^2>\t<N>\t<N^2>\t<N(g)>\t<N(g)^2>\t<N(e)>\t<N(e)^2>\t<Sz>\t<Sz^2>\t<double g>\t<double e>\t";
    *out<<"params[ts,Us,V,Vex,mus]";
    *out<<endl;
    out->close();delete out;
  }
  out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output (append="<<app<<")"<<endl;
    exit(1);
  }
  *out<<setprecision(15);

  *out<<L<<"\t"<<D<<"\t"<<E0<<"\t"<<real(valH2)<<"\t";
  *out<<real(valNf)<<"\t"<<real(valNf2)<<"\t";
  *out<<real(valNfg)<<"\t"<<real(valNfg2)<<"\t";
  *out<<real(valNfe)<<"\t"<<real(valNfe2)<<"\t";
  *out<<real(valSz)<<"\t"<<real(valSz2)<<"\t";
  *out<<real(valDg)<<"\t"<<real(valDe)<<"\t";
  *out<<tg0<<"\t"<<tg1<<"\t"<<te0<<"\t"<<te1<<"\t";
  *out<<Ug<<"\t"<<Ue<<"\t"<<V<<"\t"<<Vex<<"\t";
  *out<<mu_g0<<"\t"<<mu_g1<<"\t"<<mu_e0<<"\t"<<mu_e1<<"\t";
  *out<<endl;

  out->close();delete out;
  if(file_exists(tmpMPSfile))
    remove(tmpMPSfile.data()); // remove potential tmp file

  if(relaunch){// prepare the next job
    int nextSz=Sz+incrSz;
    double nextpenNtot(penNtot),nextpenNg(penNg),nextpenNe(penNe),nextpenSz(penSz);
    bool increasePen=false;
    if(penNtot>0&&abs(real(valNf)-Ntot)>1E-2){ // got the wrong Nf
      cout<<"Result with wrong Ntot "<<valNf<<endl;
      increasePen=true;nextpenNtot+=abs(E0)/((real(valNf)-Ntot)*(real(valNf)-Ntot));}
    if(penNg>0&&abs(real(valNfg)-Ng)>1E-2){ // got the wrong Ng
      cout<<"Result with wrong Ng "<<valNfg<<endl;
      increasePen=true;nextpenNg*=2.;}//+=abs(E0)/((real(valNfg)-Ng)*(real(valNfg)-Ng));}
    if(penNe>0&&abs(real(valNfe)-Ne)>1E-2){ // got the wrong Ne
      cout<<"Result with wrong Ne "<<valNfe<<" (target "<<Ne<<", diff "<<abs(real(valNfe)-Ne)<<")"<<endl;
      increasePen=true;nextpenNe*=2.;} //nextpenNe+=abs(E0)/((real(valNfe)-Ne)*(real(valNfe)-Ne));}
    if(penNg>0&&abs(real(valSz)-Sz)>1E-2){ // got the wrong Sz
      cout<<"Result with wrong Sz "<<valSz<<endl;
      increasePen=true;nextpenSz*=2.;nextSz=Sz; // not changing Sz yet
    } 
    if(increasePen){ // should maybe erase the file and run the same with higher penalty!
      cout<<"WARNING: The computed state is probably not good enough!!!"<<endl;
    }
      
    //    string basempsfileNext=getMPSfilenameConstraints(L,penNg,Ng,penNe,Ne,penSz,nextSz,penNtot,Ntot,t,U,V,Vex,mu_g,mu_e,1); // basename of MPS files
    //prepareNewJob(jobsdir,basempsfileNext,nextSz,nextpenNtot,nextpenNg,nextpenNe,nextpenSz,argv[0],argv[1],props);

    if(abs(valS2-Stot*(Stot+1))<.1*Stot*(Stot+1)){
      if(incrD<=0||maxD<=0||D==minD){
	cout<<"Launching for another Sz"<<endl;
	string basempsfileNext=getMPSfilenameConstraints(L,penNg,Ng,penNe,Ne,penSz,Sz+incrSz,penNtot,Ntot,t,U,V,Vex,mu_g,mu_e,1); // basename of MPS files
	prepareNewJob(jobsdir,basempsfileNext,Sz+incrSz,penNtot,penNg,penNe,penSz,argv[0],argv[1],props);	
      }
      if(incrD>0&&maxD>0&&D+incrD<=maxD){
	cout<<"Preparing (also) a larger D run"<<endl;
	prepareNewJobNextD(jobsdir,basempsfile,D+incrD,penNtot,penNg,penNe,penSz,argv[0],argv[1],props);	
      }
    }
    else{
      cout<<"ERROR: Since the total found S2 ("<<valS2<<") differs from the desired one ("<<Stot*(Stot+1)<<"), I stop the jobs here!"<<endl;
      exit(2);
    }

  }
	  
}

int findTmpFile(int D,const string& mpsdir,const string& baseName,string& mpsfile){
  mpsfile=mpsFileName(mpsdir,baseName,D)+"_tmp";
  bool found=0;
  if(file_exists(mpsfile)){
    cout<<"findTmpFile found mpsfile="<<mpsfile<<endl;
    found=true;
    return D;
  }
  else{// search for smaller Ds
     cout<<"Not found tmp file? "<<mpsfile<<endl;
     int cnt=D;
     bool found=0;
     while(cnt>0&&!found){
       mpsfile=mpsFileName(mpsdir,baseName,cnt);
       if(file_exists(mpsfile)){
   	cout<<"findTmpFile found mpsfile="<<mpsfile<<endl;
   	found=true;
   	return cnt;
       }
       else{
   	mpsfile="";
   	cnt--;
       }
     }
  }
  if(!found) return 0;
  cout<<"ERROR: Should not be here"<<endl;
  exit(1);
}


const string jobFileName(const string& jobsdir,const string& baseName,int D){
  stringstream s;
  s<<jobsdir<<"/job_Stot_"<<baseName<<"_"<<D;
  return s.str();
}

const string mpsFileName(const string& mpsdir,const string& baseName,int D){
  stringstream s;
  s<<mpsdir<<"/"<<baseName;
  if(D>0) s<<"_"<<D;
  return s.str();
}

const string getMPSfilename(int L,const vector<double>& t,const vector<double>& U,double V,double Vex,const vector<double>& mu_g,const vector<double>& mu_e,bool GS){
  stringstream s;
  if(GS)
    s<<"gsFH";
  else
    s<<"excFH";
  s<<"_L"<<L<<"_t";
  for(int k=0;k<t.size();k++)s<<t[k]<<"_";
  s<<"_U";for(int k=0;k<U.size();k++)s<<U[k]<<"_";
  s<<"_V"<<V<<"_Vex"<<Vex<<"_mug";
  for(int k=0;k<mu_g.size();k++)s<<mu_g[k]<<"_";
  s<<"_mue";for(int k=0;k<mu_e.size();k++)s<<mu_e[k]<<"_";
  return s.str();
}


const string getMPSfilenameConstraints(int L,double penNg,int Ng,double penNe,int Ne,double penSz,int Sz,double penNtot,int Ntot,const vector<double>& t,const vector<double>& U,double V,double Vex,const vector<double>& mu_g,const vector<double>& mu_e,bool GS){
  stringstream s;
  if(GS)
    s<<"gsFH";
  else
    s<<"excFH";
  s<<"_L"<<L;
  s<<"_Ng";if(penNg>0) s<<Ng; else s<<"X";
  s<<"_Ne";if(penNe>0) s<<Ne; else s<<"X";
  s<<"_Sz"; if(penSz>0) s<<Sz; else s<<"X";
  s<<"_Ntot"; if(penNtot>0) s<<Ntot; else s<<"X";
  s<<"_t";
  for(int k=0;k<t.size();k++)s<<t[k]<<"_";
  s<<"_U";for(int k=0;k<U.size();k++)s<<U[k]<<"_";
  s<<"_V"<<V<<"_Vex"<<Vex<<"_mug";
  for(int k=0;k<mu_g.size();k++)s<<mu_g[k]<<"_";
  s<<"_mue";for(int k=0;k<mu_e.size();k++)s<<mu_e[k]<<"_";
  return s.str();
}


// Select a number of random Nsel indices from 0 to Nfrom-1
void selectRandomSet(int Nsel,int Nfrom,vector<int>& selection){
  selection.clear(); // just in case
  vector<int> inds;
  for(int p=0;p<Nsel;p++){
    inds.push_back(random()%(Nfrom-p));
  }
  //cout<<"Selected random positions "<<inds<<endl;
  // now make them correspond to indices
  // straightforward: take them off the list of initial indices
  list<int> allInds;
  for(int p=0;p<Nfrom;p++) allInds.push_back(p);
  for(int k=0;k<Nsel;k++){
    list<int>::iterator iter=allInds.begin();
    advance(iter,inds[k]); // not very efficient, as it traverses the list one by one!
    selection.push_back(*iter);
    allInds.erase(iter);
    //cout<<"Now selection:"<<selection<<" and list:";
    //iter=allInds.begin();while(iter!=allInds.end())cout<<*iter++<<",";
    //cout<<endl;
  }
  //cout<<"Corresponding to indices "<<selection<<endl;
  // Q: should I compute them, without creating the list? is it more efficient?
}

void prepareStateWithNfermions(MPS& mps,int Nf,int L,int D){
  complex_t datav0[]={ONE_c,ZERO_c};
  mwArray v0(Indices(2,1,1),datav0);//
  complex_t datav1[]={ZERO_c,ONE_c};
  mwArray v1(Indices(2,1,1),datav1);//
  mps=MPS(L,1,2);mps.setProductState(p_one); // all empty corresponds to 1s
  // now set Nf sites to 1
  vector<int> pos1;
  selectRandomSet(Nf,L,pos1);
  //cout<<"Selected the random set of positions "<<pos1<<endl;
  for(int k=0;k<Nf;k++){
    mps.setA(pos1[k],v0);
  }
  if(D>1) mps.increaseBondDimension(D);
}


void prepareNewJob(const string& jobsdir,const string& basempsfile,int nextSz,double nextpenNtot,double nextpenNg,double nextpenNe,double nextpenSz,const string& exec,const string& propFile,const Properties& props){
  int L=props.getIntProperty("L");
  double tg0=props.getDoubleProperty("tg0");
  double tg1=props.getDoubleProperty("tg1");
  double te0=props.getDoubleProperty("te0");
  double te1=props.getDoubleProperty("te1");
  double Ug=props.getDoubleProperty("Ug");
  double Ue=props.getDoubleProperty("Ue");
  double V=props.getDoubleProperty("V");
  double Vex=props.getDoubleProperty("Vex");
  double mu_g0=props.getDoubleProperty("mu_g0");
  double mu_g1=props.getDoubleProperty("mu_g1");
  double mu_e0=props.getDoubleProperty("mu_e0");
  double mu_e1=props.getDoubleProperty("mu_e1");
  int Ntot=props.getIntProperty("Ntot");
  if(Ntot<0) Ntot=0;
  int Ng=props.getIntProperty("Ng");
  if(Ng<0) Ng=0;
  int Ne=props.getIntProperty("Ne");
  if(Ne<0) Ng=0;
  int Sz=props.getIntProperty("Sz");
  int Stot=props.getIntProperty("Stot");
  int Szprev=props.getIntProperty("Szprev");
  if(abs(nextSz*.5)>Stot){
    cout<<"Extreme Sz reached, no job prepared"<<endl;
    return;
  }
  if(Sz==0){
    cout<<"Computed up to the middle, stopping here "<<endl;
    return;
  }
  int prevSz=nextSz==Sz?Szprev:Sz;
  int D=props.getIntProperty("D");
  int incrD=props.getIntProperty("incrD");
  int maxD=props.getIntProperty("maxD");
  int minD=props.getIntProperty("minD"); // the one used to sweep the first time
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  string outfname=props.getProperty("output"); // this has at the end the old Sz!
  int posSz=outfname.find_last_of("Sz");
  outfname=outfname.substr(0,posSz+1);
  outfname.append(to_string(nextSz));
  cout<<"Next output file will be "<<outfname<<endl;
  string mpsdir=props.getProperty("mpsdir");
  bool app=1;
  //  app=!(props.getIntProperty("append")==0); 
  bool relaunch=props.getIntProperty("relaunch")==1;
  // int incrD=props.getIntProperty("incrD");
  // int maxD=props.getIntProperty("maxD");

  int freqSv=props.getIntProperty("freqSv"); // save to file every so many sweeps  
  // bool refine=(props.getIntProperty("refine")==1||(nextD==props.getIntProperty("D"))||(L==24&&nextD==260)); // remove the last one! 
  // double incrVtry=props.getDoubleProperty("tryVsstep");

  const string jobfile=jobFileName(jobsdir,basempsfile,D);
  stringstream commandstr;
  commandstr<<exec<<" "<<propFile; // command and properties file
  commandstr<<" -L="<<L<<" -freqSv="<<freqSv
	    <<" -tg0="<<tg0<<" -tg1="<<tg1<<" -te0="<<te0<<" -te1="<<te1
	    <<" -Ug="<<Ug<<" -Ue="<<Ue<<" -V="<<V<<" -Vex="<<Vex
	    <<" -mu_g0="<<mu_g0<<" -mu_g1="<<mu_g1<<" -mu_e0="<<mu_e0<<" -mu_e1="<<mu_e1;
  commandstr<<" -tol="<<tol<<" -output="<<outfname<<" -mpsdir="<<mpsdir<<" -append="<<app;
  commandstr<<" -D="<<D;//<<" -Ntot="<<Ntot<<" -penalty="<<nextpen;
  if(nextpenNtot>0) commandstr<<" -Ntot="<<Ntot<<" -penaltyNtot="<<nextpenNtot;
  if(nextpenNg>0) commandstr<<" -Ng="<<Ng<<" -penaltyNg="<<nextpenNg;
  if(nextpenNe>0) commandstr<<" -Ne="<<Ne<<" -penaltyNe="<<nextpenNe;
  commandstr<<" -Sz="<<nextSz<<" -penaltySz="<<nextpenSz;
  commandstr<<" -Szprev="<<prevSz;
  commandstr<<" -Stot="<<Stot;
  commandstr<<" -relaunch="<<relaunch<<" -jobsdir="<<jobsdir;
  commandstr<<" -incrD="<<incrD<<" -maxD="<<maxD<<" -minD="<<minD;
    
  ofstream* ojob=new ofstream(jobfile.data());
  if(!ojob->is_open()){
    cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
    cout<<commandstr.str()<<endl;
    cout<<endl;
    // Do not exit, because I still need to save the MPS
  }
  else{
    int minNumber=1;
    *ojob<<"#!/bin/bash"<<endl;
    *ojob<<"msub_slurm -N D"<<D<<".L"<<L<<".Ng"<<Ng<<".Ne"<<Ne<<".V"<<V<<".Sz"<<nextSz<<".Stot -mem="
	 <<ceil(L*D*D*(12.+4.*minNumber)*16E-6)+200<<"M ";
    //if(L==12) *ojob<<"-P 2 ";
    //else{
    if(D<200&&L<12) *ojob<<"-P 8 ";
    else *ojob<<"-P 20 ";
    //}
    // time limit, to be removed!
    //    *ojob<<"-t 3:00:00 ";
    *ojob<<"-- ";
    *ojob<<commandstr.str()<<endl;
    ojob->close();
  }
  delete ojob;
}

bool checkQuantumNumber(const MPS& gs,const MPO& mpo,int targetV,double tolVal,double tolVar){
  Contractor& contractor=Contractor::theContractor();
  complex_t valOp=contractor.contract(gs,mpo,gs);
  complex_t valOp2=contractor.contract2(mpo,gs);
  if((abs(valOp-targetV*ONE_c)>tolVal)&&
     (abs(valOp2-valOp*valOp)<abs(tolVar*valOp*valOp))){
    cout<<"State does not have the desired quantum nr ("<<targetV<<"), but "<<real(valOp)
	<<", with variance "<<real(valOp2-valOp*valOp)<<endl;
    return false;
  }
  else{ // remove for real runs
    cout<<"State has acceptable quantum nr (target="<<targetV<<", value="<<real(valOp)
	<<"), with variance "<<real(valOp2-valOp*valOp)<<endl;
    return true;
  }
}

void prepareNewJobNextD(const string& jobsdir,const string& basempsfile,int nextD,double nextpenNtot,double nextpenNg,double nextpenNe,double nextpenSz,const string& exec,const string& propFile,const Properties& props){
  int L=props.getIntProperty("L");
  double tg0=props.getDoubleProperty("tg0");
  double tg1=props.getDoubleProperty("tg1");
  double te0=props.getDoubleProperty("te0");
  double te1=props.getDoubleProperty("te1");
  double Ug=props.getDoubleProperty("Ug");
  double Ue=props.getDoubleProperty("Ue");
  double V=props.getDoubleProperty("V");
  double Vex=props.getDoubleProperty("Vex");
  double mu_g0=props.getDoubleProperty("mu_g0");
  double mu_g1=props.getDoubleProperty("mu_g1");
  double mu_e0=props.getDoubleProperty("mu_e0");
  double mu_e1=props.getDoubleProperty("mu_e1");
  int Ntot=props.getIntProperty("Ntot");
  if(Ntot<0) Ntot=0;
  int Ng=props.getIntProperty("Ng");
  if(Ng<0) Ng=0;
  int Ne=props.getIntProperty("Ne");
  if(Ne<0) Ng=0;
  int Sz=props.getIntProperty("Sz");
  int Stot=props.getIntProperty("Stot");
  int Szprev=props.getIntProperty("Szprev");
  int D=nextD;
  int incrD=props.getIntProperty("incrD");
  int maxD=props.getIntProperty("maxD");
  int minD=props.getIntProperty("minD"); // the one used to sweep the first time
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  string outfname=props.getProperty("output"); // this has at the end the old Sz!
  // int posSz=outfname.find_last_of("Sz");
  // outfname=outfname.substr(0,posSz+1);
  // outfname.append(to_string(nextSz));
  //  cout<<"Next output file will be "<<outfname<<endl;
  string mpsdir=props.getProperty("mpsdir");
  bool app=1;
  //  app=!(props.getIntProperty("append")==0); 
  bool relaunch=props.getIntProperty("relaunch")==1;
  // int incrD=props.getIntProperty("incrD");
  // int maxD=props.getIntProperty("maxD");

  int freqSv=props.getIntProperty("freqSv"); // save to file every so many sweeps  
  // bool refine=(props.getIntProperty("refine")==1||(nextD==props.getIntProperty("D"))||(L==24&&nextD==260)); // remove the last one! 
  // double incrVtry=props.getDoubleProperty("tryVsstep");

  const string jobfile=jobFileName(jobsdir,basempsfile,D);
  stringstream commandstr;
  commandstr<<exec<<" "<<propFile; // command and properties file
  commandstr<<" -L="<<L<<" -freqSv="<<freqSv
	    <<" -tg0="<<tg0<<" -tg1="<<tg1<<" -te0="<<te0<<" -te1="<<te1
	    <<" -Ug="<<Ug<<" -Ue="<<Ue<<" -V="<<V<<" -Vex="<<Vex
	    <<" -mu_g0="<<mu_g0<<" -mu_g1="<<mu_g1<<" -mu_e0="<<mu_e0<<" -mu_e1="<<mu_e1;
  commandstr<<" -tol="<<tol<<" -output="<<outfname<<" -mpsdir="<<mpsdir<<" -append="<<app;
  commandstr<<" -D="<<D;//<<" -Ntot="<<Ntot<<" -penalty="<<nextpen;
  if(nextpenNtot>0) commandstr<<" -Ntot="<<Ntot<<" -penaltyNtot="<<nextpenNtot;
  if(nextpenNg>0) commandstr<<" -Ng="<<Ng<<" -penaltyNg="<<nextpenNg;
  if(nextpenNe>0) commandstr<<" -Ne="<<Ne<<" -penaltyNe="<<nextpenNe;
  commandstr<<" -Sz="<<Sz<<" -penaltySz="<<nextpenSz;
  commandstr<<" -Szprev="<<Szprev;
  commandstr<<" -Stot="<<Stot;
  commandstr<<" -relaunch="<<relaunch<<" -jobsdir="<<jobsdir;
  commandstr<<" -incrD="<<incrD<<" -maxD="<<maxD<<" -minD="<<minD;

    
  ofstream* ojob=new ofstream(jobfile.data());
  if(!ojob->is_open()){
    cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
    cout<<commandstr.str()<<endl;
    cout<<endl;
    // Do not exit, because I still need to save the MPS
  }
  else{
    int minNumber=1;
    *ojob<<"#!/bin/bash"<<endl;
    *ojob<<"msub_slurm -N D"<<D<<".L"<<L<<".Ng"<<Ng<<".Ne"<<Ne<<".V"<<V<<".Sz"<<Sz<<".Stot -mem="
	 <<ceil(L*D*D*(12.+4.*minNumber)*16E-6)+200<<"M ";
    //if(L==12) *ojob<<"-P 2 ";
    //else{
    if(D<200&&L<12) *ojob<<"-P 8 ";
    else *ojob<<"-P 20 ";
    //}
    // time limit, to be removed!
    //    *ojob<<"-t 3:00:00 ";
    *ojob<<"-- ";
    *ojob<<commandstr.str()<<endl;
    ojob->close();
  }
  delete ojob;


}
