
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

/** twofermihubbardGS runs the findGroundState routine with the MPO for the
    Hamiltonian \ref <TwoOrbitalFermiHubbardHamiltonian>, and computes GS and
    maybe other things... 

    Receives arguments:
    \param <L> (int) length of the chain
    \param <t> (double) parameter \f$t\f$ of the Hamiltonian
    \param <U> (double) parameter \f$U\f$ of the Hamiltonian
    \param <V> (double) parameter \f$V\f$ of the Hamiltonian
    \param <Vex> (double) parameter \f$V_{\mathrm ex}\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the energy

    In this version (28.7.2022) a first warmup sweep is done for D=2, and then, for larger dimensions (with a hardcoded increment between 
    D=2 and D0) the usual optimization is run starting from smaller D.
*/

int findTmpFile(int D,const string& mpsdir,const string& baseName,string& mpsfile);
const string mpsFileName(const string& mpsdir,const string& baseName,int D);
const string getMPSfilename(int L,const vector<double>& t,const vector<double>& U,double V,double Vex,const vector<double>& mu_g,const vector<double>& mu_e,bool GS=1);

const string getMPSfilenameConstrains(int L,double penNg,int Ng,double penNe,int Ne,double penSz,int Sz,double penNtot,int Ntot,const vector<double>& t,const vector<double>& U,double V,double Vex,const vector<double>& mu_g,const vector<double>& mu_e,bool GS=1);

bool checkQuantumNumber(const MPS& gs,const MPO& mpo,int targetV,double tolVal,double tolVar);


void prepareStateWithNfermions(MPS& mps,int Nf,int L,int D);

bool checkWithBestImbalance(double E0,const string& lowestEnergyFile);
const string getLowestEfile(const string& lowEdir,int L,int D,const vector<double>& t,const vector<double>& U,double V,double Vex,const vector<double>& mu_g,const vector<double>& mu_e);

/** 
Generate the name of the file with the command in it
*/

const string jobFileName(const string& jobsdir,const string& baseName,int D);
void prepareNewJob(const string& jobsdir,const string& basempsfile,int nextD,double nextpenNtot,double nextpenNg,double nextpenNe,double nextpenSz,
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
  if(penSz<0) penSz=0;
  int Sz=props.getIntProperty("Sz");

  if(penNe>0&&penNg>0) penNtot=0;
  
  int D=props.getIntProperty("D");
  int warmup=props.getIntProperty("warmup"); // how many rounds of warmup
  if(warmup<0) warmup=0;
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir"); // where are the MPS
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  int freqSv=props.getIntProperty("freqSv"); // save to file every so many sweeps
  if(freqSv<0) freqSv=10;
  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  // To prepare new jobs that rerun with increased D (or, if necessary, increased penalty)
  bool relaunch=props.getIntProperty("relaunch")==1;
  int incrD=0;int maxD=D;string jobsdir;
  if(relaunch){
    incrD=props.getIntProperty("incrD");
    maxD=props.getIntProperty("maxD");
    jobsdir=props.getProperty("jobsdir"); 
  }
  bool multipleImbalances=props.getIntProperty("multiImbalance")==1;
  string lowEdir;
  if(multipleImbalances){
    lowEdir=props.getProperty("lowEdir"); 
  }
  bool refine=(props.getIntProperty("refine")==1); // in that case, only runs if the same state (same D, no tmp) exists. the same is passed on
  bool restart=(props.getIntProperty("restart")==1); // in that case, only runs if the state with D0 exists, but not passed on to new job
  int Drestart=D;
  if(restart){
    Drestart=props.getIntProperty("Drestart");
    if(Drestart<=0) restart=0;
  }
  bool tryVs=(props.getIntProperty("tryVs")==1); // in that case, check the best solution with neighbouring Vs, and, if it exists and has better energy, start with that one
  double incrVtry=props.getDoubleProperty("tryVsstep");
  if(incrVtry<=0) incrVtry=1.;

  refine=0;restart=0;tryVs=0; // in this version, I don't do this
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
  if(penSz>0) cout<<"targetSz="<<Sz<<"; lambda="<<penSz<<";"<<endl;
  cout<<"L="<<L<<";"<<endl;
  cout<<"[H,Nfg,Nfe,Sz]=FermiHubbardHamiltonianTwoOrbital(L,tg0,tg1,te0,te1,Ug,Ue,V,Vex,mu_g0,mu_g1,mu_e0,mu_e1);\nvalsEx=eig(H);\nmin(valsEx)"<<endl;
  
  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);
  //  cout<<"Initialized Contractor"<<endl;
  int d=2;
  //  int dF=d*d; // dimension of fermion sites (L+1)

  MPS gs(4*L,D,d);gs.setRandomState();gs.gaugeCond('R',1);
  // First, check whether there is a tmp file

  // Prepare parameters
  vector<double> t;t.push_back(tg0);t.push_back(tg1);t.push_back(te0);t.push_back(te1);
  vector<double> U;U.push_back(Ug);U.push_back(Ue);
  vector<double> mu_g;mu_g.push_back(mu_g0);mu_g.push_back(mu_g1);
  vector<double> mu_e;mu_e.push_back(mu_e0);mu_e.push_back(mu_e1);

  string basempsfile=getMPSfilename(L,t,U,V,Vex,mu_g,mu_e,1); // basename of MPS files
  if(penNtot+penNg+penNe+penSz>0)
    basempsfile=getMPSfilenameConstrains(L,penNg,Ng,penNe,Ne,penSz,Sz,penNtot,Ntot,t,U,V,Vex,mu_g,mu_e,1); // basename of MPS files
  string mpsfile=mpsFileName(mpsdir,basempsfile,D);
  string tmpMPSfile;

  if(file_exists(mpsfile)){
    cout<<"Refining existing state "<<mpsfile<<endl;
    gs.importMPS(mpsfile.data());
  }
  else{
    // check for other files
    if(findTmpFile(D,mpsdir,basempsfile,tmpMPSfile)){// found some tmp or smaller D file
      cout<<"Recovered initial guess from file "<<tmpMPSfile<<endl;
      gs.importMPS(tmpMPSfile.data());
      if(gs.getLength()!=4*L){
	cout<<"Error! File "<<tmpMPSfile<<" does not contain valid MPS "<<endl;
	exit(1);
      }
      gs.increaseBondDimensionWithNoise(D,.001);
    }
    else{
      cout<<"Failed to recover a tmp file ("<<basempsfile<<")"<<endl;
      gs.setRandomState(); // the intial state, random
      //    cout<<"Initialized random state from scratch "<<endl;
      // if(penalty>0){
      //   prepareStateWithNfermions(gs,Ntot,4*L,1);
      //   gs.gaugeCond('R',1);
      //   gs.gaugeCond('L',1);
      //   cout<<"Constructed a state with (supposedly) Nferm="<<Ntot<<endl;
      //   cout<<"Checking:"<<contractor.contract(gs,fermNr,gs)<<endl;
      //   if(4*L<=12) gs.exportForMatlab("gsMPS.m");
      // }
    }
  }
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  //  cout<<"Initial state has norm "<<contractor.contract(gs,gs)<<endl;
  
  // First: put H and find the GS
  double E0=0.;
  TwoOrbitalFermiHubbardHamiltonian ham(L,t,U,V,Vex,mu_g,mu_e,Ntot,penNtot);
  if(penNg+penNe+penSz>0){
    ham.setPenalties(Ng,penNg,Ne,penNe,Sz,penSz,Ntot,penNtot);
  }
    cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=ham.getHMPO();
  E0=real(contractor.contract(gs,hamil,gs)); // current value of the energy


  if(warmup){
    if(1){
    cout<<"Applying "<<warmup<<" symmetric warm up rounds"<<endl;
    for(int pp=0;pp<warmup;pp++){
      for(int pos=0;pos<(4*L-1)/2+1;pos++){
    	contractor.sweepPart(pos,pos+1,hamil,D,gs,0.); // one right
    	contractor.sweepPart(4*L-1-pos,4*L-1-pos-1,hamil,D,gs,0.); // one left
      }
      double Ewarm=real(contractor.contract(gs,hamil,gs));
      cout<<"After symmetric sweep nr "<<pp+1<<" energy value "<<Ewarm<<endl;
    }
    }
    else{
    cout<<"Applying "<<warmup<<" centered warm up rounds"<<endl;
    for(int pp=0;pp<warmup;pp++){
      //      for(int pos=0;pos<L/2-1;pos++){
      for(int pos=0;pos<(4*L-1)/2+1;pos++){
	int leftB=4*L/2-1-pos;int rightB=4*L/2+pos;
    	//cout<<"Sweep from "<<leftB<<"-"<<rightB<<endl;
    	contractor.sweepPart(leftB,rightB,hamil,D,gs,0.); // 2*(pos+1) right
    	contractor.sweepPart(rightB,leftB,hamil,D,gs,0.); // same left
      }
      double Ewarm=real(contractor.contract(gs,hamil,gs));
      cout<<"After centered sweep nr "<<pp+1<<" energy value "<<Ewarm<<endl;
    }
    }
  }
    

  
  // // if I am trying to refine, I should also check for better guesses as smaller Ds
  // if(refine){ // maybe only for L=50?
  //   int Ds=D-1;
  //   while(Ds>0){
  //     string mpsfileDs=mpsFileName(mpsdir,basempsfile,Ds);
  //     if(file_exists(mpsfileDs)){
  // 	MPS mpsDs;
  // 	mpsDs.importMPS(mpsfileDs.data());
  // 	double EDs=real(contractor.contract(mpsDs,hamil,mpsDs));
  // 	if(EDs<E0){
  // 	  gs=mpsDs;E0=EDs;
  // 	  cout<<"Lower energy found for smaller D="<<Ds<<endl;
  // 	  gs.increaseBondDimensionWithNoise(D,.001);
  // 	}
  //     }
  //     Ds--;
  //   }
  // }

  
  MPO fermNr(4*L);
  ham.getFermionNrMPO(fermNr); // mpo with total nr of fermions
  MPO fermNr0(4*L);
  ham.getFermionNrOrbitalMPO(fermNr0,true); // mpo with total nr of fermions in g
  MPO fermNr1(4*L);
  ham.getFermionNrOrbitalMPO(fermNr1,false); // mpo with total nr of fermions in e
  MPO fermSz(4*L);
  ham.getTotalSzMPO(fermSz);
  cout<<"Constructed the hamil MPO and Nferm MPOs"<<endl;
  if(4*L<=12){
    hamil.exportForMatlab("fH2Hmpo.m");
    fermNr.exportForMatlab("Nfmpo.m");
    fermNr0.exportForMatlab("Nf0mpo.m");
    fermNr1.exportForMatlab("Nf1mpo.m");
    fermSz.exportForMatlab("Szmpo.m");
  }

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

  
  cout<<"Initial state, with initial energy ";
  cout<<contractor.contract(gs,hamil,gs)<<endl;
  cout<<"\t and initial nr of fermions "<<contractor.contract(gs,fermNr,gs)<<endl;

   
  tmpMPSfile=mpsfile+"_tmp";//int freqSv=10; // save to file every so many sweeps
  double offset=0.;int nrEig=1;
  // sice I am running, I will move a posssible final file to tmp status, in case this run is interrupted!
  if(file_exists(mpsfile)){
    rename(mpsfile.data(),tmpMPSfile.data());
    cout<<"Renamed wxisting "<<mpsfile<<" to "<<tmpMPSfile<<endl;
  }
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
  cout<<" <H^2>="<<real(valH2);
  cout<<"\t and nrs of fermions:"<<endl;
  cout<<"\tNtot "<<valNf<<" variance:"<<real(valNf2-valNf*valNf)<<endl;
  cout<<"\tNg "<<valNfg<<" variance:"<<real(valNfg2-valNfg*valNfg)<<endl;
  cout<<"\tNe "<<valNfe<<" variance:"<<real(valNfe2-valNfe*valNfe)<<endl;
  cout<<"\tSz "<<valSz<<" variance:"<<real(valSz2-valSz*valSz)<<endl;
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
    int nextD;double nextpenNtot(penNtot),nextpenNg(penNg),nextpenNe(penNe),nextpenSz(penSz);
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
      increasePen=true;nextpenSz*=2.;} //+=abs(E0)/((real(valSz)-Sz)*(real(valSz)-Sz));}
    if(!increasePen){
      nextD=D+incrD;
    }
    else{
      nextD=D;
    }
    if(nextD<=maxD){
      // if this energy was larger than the best one obtained, maybe stop this imbalance
      bool doLaunch=true;
      if(!increasePen&&multipleImbalances){
	string lowestEfile=getLowestEfile(lowEdir,L,D,t,U,V,Vex,mu_g,mu_e);
	doLaunch=checkWithBestImbalance(E0,lowestEfile);
      }
      if(doLaunch)
	prepareNewJob(jobsdir,basempsfile,nextD,nextpenNtot,nextpenNg,nextpenNe,nextpenSz,argv[0],argv[1],props);
    }
  }
	  
}

int findTmpFile(int D,const string& mpsdir,const string& baseName,string& mpsfile){
  mpsfile=mpsFileName(mpsdir,baseName,D)+"_tmp";
  bool found=0;
  if(file_exists(mpsfile)){
    cout<<"findTmpFile found mpsfile="<<mpsfile<<endl;
    found=true;
    return 1;
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
  s<<jobsdir<<"/job_"<<baseName<<"_"<<D;
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

const string getLowestEfile(const string& lowEdir,int L,int D,const vector<double>& t,const vector<double>& U,double V,double Vex,const vector<double>& mu_g,const vector<double>& mu_e){
  stringstream s;
  s<<lowEdir<<"/lowE";
  s<<"_L"<<L<<"_D"<<D<<"_t";
  for(int k=0;k<t.size();k++)s<<t[k]<<"_";
  s<<"_U";for(int k=0;k<U.size();k++)s<<U[k]<<"_";
  s<<"_V"<<V<<"_Vex"<<Vex<<"_mug";
  for(int k=0;k<mu_g.size();k++)s<<mu_g[k]<<"_";
  s<<"_mue";for(int k=0;k<mu_e.size();k++)s<<mu_e[k]<<"_";
  return s.str();
}

const string getMPSfilenameConstrains(int L,double penNg,int Ng,double penNe,int Ne,double penSz,int Sz,double penNtot,int Ntot,const vector<double>& t,const vector<double>& U,double V,double Vex,const vector<double>& mu_g,const vector<double>& mu_e,bool GS){
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


void prepareNewJob(const string& jobsdir,const string& basempsfile,int nextD,double nextpenNtot,double nextpenNg,double nextpenNe,double nextpenSz,const string& exec,const string& propFile,const Properties& props){
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
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir");
  bool app=1;
  //  app=!(props.getIntProperty("append")==0); 
  bool relaunch=props.getIntProperty("relaunch")==1;
  int incrD=props.getIntProperty("incrD");
  int maxD=props.getIntProperty("maxD");
  // by hand, manipulate the increments
  if(nextD==maxD){
    if(maxD==50){
      incrD=10;maxD=100;
    }
    if(maxD==100){
      incrD=20;maxD=160;
    }
  }

  int freqSv=props.getIntProperty("freqSv"); // save to file every so many sweeps  
  bool multipleImbalances=props.getIntProperty("multiImbalance")==1;
  string lowEdir;
  if(multipleImbalances){
    lowEdir=props.getProperty("lowEdir"); 
  }
  bool refine=(props.getIntProperty("refine")==1||(nextD==props.getIntProperty("D"))||(L==24&&nextD==260)); // remove the last one! 
  double incrVtry=props.getDoubleProperty("tryVsstep");

  const string jobfile=jobFileName(jobsdir,basempsfile,nextD);
  stringstream commandstr;
  commandstr<<exec<<" "<<propFile; // command and properties file
  commandstr<<" -L="<<L<<" -freqSv="<<freqSv
	    <<" -tg0="<<tg0<<" -tg1="<<tg1<<" -te0="<<te0<<" -te1="<<te1
	    <<" -Ug="<<Ug<<" -Ue="<<Ue<<" -V="<<V<<" -Vex="<<Vex
	    <<" -mu_g0="<<mu_g0<<" -mu_g1="<<mu_g1<<" -mu_e0="<<mu_e0<<" -mu_e1="<<mu_e1;
  commandstr<<" -tol="<<tol<<" -output="<<outfname<<" -mpsdir="<<mpsdir<<" -append="<<app;
  commandstr<<" -D="<<nextD;//<<" -Ntot="<<Ntot<<" -penalty="<<nextpen;
  if(nextpenNtot>0) commandstr<<" -Ntot="<<Ntot<<" -penaltyNtot="<<nextpenNtot;
  if(nextpenNg>0) commandstr<<" -Ng="<<Ng<<" -penaltyNg="<<nextpenNg;
  if(nextpenNe>0) commandstr<<" -Ne="<<Ne<<" -penaltyNe="<<nextpenNe;
  if(nextpenSz>0) commandstr<<" -Sz="<<Sz<<" -penaltySz="<<nextpenSz;
  commandstr<<" -relaunch="<<relaunch<<" -incrD="<<incrD<<" -maxD="<<maxD<<" -jobsdir="<<jobsdir;
  if(multipleImbalances) commandstr<<" -multiImbalance=1 -lowEdir="<<lowEdir;
  commandstr<<" -refine="<<refine;
  commandstr<<" -tryVs=1 -tryVsstep="<<incrVtry;
  
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
    *ojob<<"msub_slurm -N D"<<nextD<<".L"<<L<<".Ng"<<Ng<<".Ne"<<Ne<<".V"<<V<<" -mem="
	 <<ceil(L*nextD*nextD*(12.+4.*minNumber)*16E-6)+200<<"M ";
    //if(L==12) *ojob<<"-P 2 ";
    //else{
    if(nextD>80&&L>12) *ojob<<"-P 8 ";
    else *ojob<<"-P 4 ";
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

bool checkWithBestImbalance(double E0,const string& lowestEnergyFile){
  // A lock file
  
  string stringlockFileName=lowestEnergyFile+"_lock";
  bool checked=true; return checked; // ignore this and run always
  FILE* fp = nullptr;
  do {        
    fp = fopen(stringlockFileName.c_str(), "wx");
  } while(!fp);
  fclose(fp);
  cout<<"The lock file "<<stringlockFileName<<" is created"<<endl;
  // Now I open the other file to read the best energy so far
  ifstream f_bestE(lowestEnergyFile);
  double bestE;
  f_bestE>>bestE;
  f_bestE.close();
  if(bestE-E0>=.2*E0){
    ofstream newE(lowestEnergyFile);
    newE<<E0<<endl;
    newE.close();    
    checked=true;
  }
  else{
    cout<<"This energy is larger than the current best:"<<bestE<<endl;
    checked=false;
  }
  // remove the lock
  remove(stringlockFileName.data());
  return checked;
}
