
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

/** twofermihubbardObs reads states computed by the GS search adn computes observables.
twofermihubbardGS runs the findGroundState routine with the MPO for the
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
*/

//int findTmpFile(int D,const string& mpsdir,const string& baseName,string& mpsfile);
const string mpsFileName(const string& mpsdir,const string& baseName,int D);
const string getMPSfilename(int L,const vector<double>& t,const vector<double>& U,double V,double Vex,const vector<double>& mu_g,const vector<double>& mu_e,bool GS=1);

const string getMPSfilenameConstrains(int L,double penNg,int Ng,double penNe,int Ne,double penSz,int Sz,double penNtot,int Ntot,const vector<double>& t,const vector<double>& U,double V,double Vex,const vector<double>& mu_g,const vector<double>& mu_e,bool GS=1);

//bool checkQuantumNumber(const MPS& gs,const MPO& mpo,int targetV,double tolVal,double tolVar);

//void prepareStateWithNfermions(MPS& mps,int Nf,int L,int D);

//bool checkWithBestImbalance(double E0,const string& lowestEnergyFile);
//const string getLowestEfile(const string& lowEdir,int L,int D,const vector<double>& t,const vector<double>& U,double V,double Vex,const vector<double>& mu_g,const vector<double>& mu_e);

/** 
Generate the name of the file with the command in it
*/

//const string jobFileName(const string& jobsdir,const string& baseName,int D);
//void prepareNewJob(const string& jobsdir,const string& basempsfile,int nextD,double nextpenNtot,double nextpenNg,double nextpenNe,double nextpenSz,
//		   const string& exec,const string& propFile,const Properties& props);


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
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir"); // where are the MPS
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  // int freqSv=props.getIntProperty("freqSv"); // save to file every so many sweeps
  // if(freqSv<0) freqSv=10;
  // if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  // // To prepare new jobs that rerun with increased D (or, if necessary, increased penalty)
  // bool relaunch=props.getIntProperty("relaunch")==1;
  int incrD=0;int maxD=D;
  //string jobsdir;
  //if(relaunch){
  incrD=props.getIntProperty("incrD");
  maxD=props.getIntProperty("maxD");
    //jobsdir=props.getProperty("jobsdir"); 
    //}
  // bool multipleImbalances=props.getIntProperty("multiImbalance")==1;
  // string lowEdir;
  // if(multipleImbalances){
  //   lowEdir=props.getProperty("lowEdir"); 
  // }
  // bool refine=(props.getIntProperty("refine")==1); // in that case, only runs if the same state (same D, no tmp) exists. the same is passed on
  // bool restart=(props.getIntProperty("restart")==1); // in that case, only runs if the state with D0 exists, but not passed on to new job
  // int Drestart=D;
  // if(restart){
  //   Drestart=props.getIntProperty("Drestart");
  //   if(Drestart<=0) restart=0;
  // }
  // bool tryVs=(props.getIntProperty("tryVs")==1); // in that case, check the best solution with neighbouring Vs, and, if it exists and has better energy, start with that one
  // double incrVtry=props.getDoubleProperty("tryVsstep");
  // if(incrVtry<=0) incrVtry=1.;
  // if(tryVs&&restart){
  //   cout<<"Option tryVs is incompatible with restart"<<endl;
  //   exit(1);
  // }  
    //  srandom(time(NULL));
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

  // Prepare parameters
  vector<double> t;t.push_back(tg0);t.push_back(tg1);t.push_back(te0);t.push_back(te1);
  vector<double> U;U.push_back(Ug);U.push_back(Ue);
  vector<double> mu_g;mu_g.push_back(mu_g0);mu_g.push_back(mu_g1);
  vector<double> mu_e;mu_e.push_back(mu_e0);mu_e.push_back(mu_e1);

  TwoOrbitalFermiHubbardHamiltonian ham(L,t,U,V,Vex,mu_g,mu_e,Ntot,penNtot);
  if(penNg+penNe+penSz>0){
    ham.setPenalties(Ng,penNg,Ne,penNe,Sz,penSz,Ntot,penNtot);
  }
  //  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=ham.getHMPO();

  MPO fermNr(4*L);
  ham.getFermionNrMPO(fermNr); // mpo with total nr of fermions
  MPO fermNr0(4*L);
  ham.getFermionNrOrbitalMPO(fermNr0,true); // mpo with total nr of fermions in g
  MPO fermNr1(4*L);
  ham.getFermionNrOrbitalMPO(fermNr1,false); // mpo with total nr of fermions in e
  MPO fermSz(4*L);
  ham.getTotalSzMPO(fermSz);
  MPO doubleG(4*L);
  ham.getDoubleOccupancyOrbitalMPO(doubleG,1);
  MPO doubleE(4*L);
  ham.getDoubleOccupancyOrbitalMPO(doubleE,0);
  cout<<"Constructed the hamil MPO, fermion nrs and double occupancy MPOs"<<endl;


  MPO corrXX(4*L);
  ham.getSpinCorrelatorMPO(corrXX,1);
  cout<<"Constructed spin corr XX MPO"<<endl;
  MPO corrYY(4*L);
  ham.getSpinCorrelatorMPO(corrYY,2);
  cout<<"Constructed spin corr YY MPO"<<endl;
  MPO corrZZ(4*L);
  ham.getSpinCorrelatorMPO(corrZZ,3);
  cout<<"Constructed spin corr ZZ MPO"<<endl;
  MPO corrTz(4*L);
  ham.getOrbitalSpinCorrelatorMPO(corrTz);
  cout<<"Constructed orbital spin corr ZZ MPO"<<endl;

  // if(4*L<=12){
  //   hamil.exportForMatlab("fH2Hmpo.m");
  //   fermNr.exportForMatlab("Nfmpo.m");
  //   fermNr0.exportForMatlab("Nf0mpo.m");
  //   fermNr1.exportForMatlab("Nf1mpo.m");
  //   fermSz.exportForMatlab("Szmpo.m");
  // }
  
  while(D<=maxD){
  
    MPS gs(4*L,D,d); //gs.setRandomState();gs.gaugeCond('R',1);

    string basempsfile=getMPSfilename(L,t,U,V,Vex,mu_g,mu_e,1); // basename of MPS files
    if(penNtot+penNg+penNe+penSz>0)
      basempsfile=getMPSfilenameConstrains(L,penNg,Ng,penNe,Ne,penSz,Sz,penNtot,Ntot,t,U,V,Vex,mu_g,mu_e,1); // basename of MPS files
    string mpsfile=mpsFileName(mpsdir,basempsfile,D);
    if(file_exists(mpsfile)){
      cout<<"Read state "<<mpsfile<<endl;
      gs.importMPS(mpsfile.data());
      gs.gaugeCond('R',1);
      gs.gaugeCond('L',1);
  
      // Compute observables
      double E0=0.;
      E0=real(contractor.contract(gs,hamil,gs)); // current value of the energy
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
      
      complex_t valDg=contractor.contract(gs,doubleG,gs);
      complex_t valDe=contractor.contract(gs,doubleE,gs);

      complex_t valSSx=contractor.contract(gs,corrXX,gs);
      complex_t valSSy=contractor.contract(gs,corrYY,gs);
      complex_t valSSz=contractor.contract(gs,corrZZ,gs);

      complex_t valTz=contractor.contract(gs,corrTz,gs);
    
      // Write results to txt file
      ofstream* out;
      if(!app||!file_exists(outfname.data())){
	out=new ofstream(outfname.data());
	if(!out->is_open()){
	  cout<<"Error: impossible to open file "<<outfname<<
	    " for output (append="<<app<<")"<<endl;
	  exit(1);
	}
	*out<<"%L\tD\t<H>\t<H^2>\t<N>\t<N^2>\t<N(g)>\t<N(g)^2>\t<N(e)>\t<N(e)^2>\t<Sz>\t<Sz^2>\t<double g>\t<double e>\t";
	*out<<"<XX>\t<YY>\t<ZZ>\t<TzTz>\t";
	*out<<"params[ts,Us,V,Vex,mus]";
	*out<<endl;
	app=1; // already initialized
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
      *out<<real(valSSx)<<"\t";
      *out<<real(valSSy)<<"\t";
      *out<<real(valSSz)<<"\t";
      *out<<real(valTz)<<"\t";
      *out<<tg0<<"\t"<<tg1<<"\t"<<te0<<"\t"<<te1<<"\t";
      *out<<Ug<<"\t"<<Ue<<"\t"<<V<<"\t"<<Vex<<"\t";
      *out<<mu_g0<<"\t"<<mu_g1<<"\t"<<mu_e0<<"\t"<<mu_e1<<"\t";
      *out<<endl;

      out->close();delete out;
    }
    else
      cout<<"No file for D="<<D<<endl;
    D+=incrD;
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

