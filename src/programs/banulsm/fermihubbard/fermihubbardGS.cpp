
#include <math.h>
#include <iomanip>

#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "misc.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "Properties.h"
#include "FermiHubbardHamiltonian.h"

using namespace shrt;

/** fermihubbardGS runs the findGroundState routine with the MPO for the
    Ising Hamiltonian \ref <FermiHubbardHamiltonian>, and computes GS and
    maybe other things... 

    Receives arguments:
    \param <L> (int) length of the chain
    \param <t> (double) parameter \f$t\f$ of the Hamiltonian
    \param <U> (double) parameter \f$U\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the energy
*/

int findTmpFile(int D,const string& mpsdir,const string& baseName,string& mpsfile);
const string mpsFileName(const string& mpsdir,const string& baseName,int D);
const string getMPSfilename(int L,double t,double U,bool GS=1);

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
  double t=props.getDoubleProperty("t");
  double U=props.getDoubleProperty("U");
  int D=props.getIntProperty("D");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir"); // where are the MPS
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  
  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)

  // cout<<"Initialized arguments: L="<<L
  //     <<", t="<<t
  //     <<", U="<<U
  //     <<", outfile="<<outfname
  //     <<", app="<<app
  //     <<", D="<<D<<endl;

  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);
  //  cout<<"Initialized Contractor"<<endl;
  int d=2;
  //  int dF=d*d; // dimension of fermion sites (L+1)

  MPS gs(L,D,d*d);
  // First, check whether there is a tmp file

  string mpsfile=getMPSfilename(L,t,U,1); // basename of MPS files
  string tmpMPSfile;
  if(findTmpFile(D,mpsdir,mpsfile,tmpMPSfile)){// found some tmp or smaller D file
    cout<<"Recovered initial guess from file "<<tmpMPSfile<<endl;
    gs.importMPS(tmpMPSfile.data());
    if(gs.getLength()!=L){
      cout<<"Error! File "<<tmpMPSfile<<" does not contain valid MPS "<<endl;
      exit(1);
    }
  }
  else{
    gs.setRandomState(); // the intial state, random
    //    cout<<"Initialized random state from scratch "<<endl;
  }
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  //  cout<<"Initial state has norm "<<contractor.contract(gs,gs)<<endl;
  
  // First: put H and find the GS
  double E0=0.;
  FermiHubbardHamiltonian ham(L,t,U);
  //  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=ham.getHMPO();
  // cout<<"Constructed the hamil MPO"<<endl;
  //   cout<<"Initial value, with initial state"<<endl;
  // cout<<contractor.contract(gs,hamil,gs)<<endl;

  contractor.findGroundState(hamil,D,&E0,gs);
  //cout<<"Ground state found with eigenvalue "<<E0<<endl;
  // save it
  tmpMPSfile=mpsFileName(mpsdir,mpsfile,D);
  gs.exportMPS(tmpMPSfile.data());
  //  cout<<"Ground state saved to "<<tmpMPSfile<<endl;

  MPO nrFerm(L);
  ham.getFermionNrMPO(nrFerm);
  complex_t nF=contractor.contract(gs,nrFerm,gs);
  //cout<<"Computing fermion number in this state: "<<nF<<endl;
  //cout<<"And the corresponding variance is: "<<contractor.contract2(nrFerm,gs)-nF*nF<<endl;
  
  MPO totSz(L);
  ham.getTotalSzMPO(totSz);
  complex_t totZ=contractor.contract(gs,totSz,gs);
  //cout<<"Computing total Sz in this state: "<<totZ<<endl;
  //cout<<"And the corresponding variance is: "<<contractor.contract2(totSz,gs)-totZ*totZ<<endl;
  cout<<"*** GS: E="<<E0<<", Nf="<<nF<<", Sz="<<totZ<<endl;


  // Same for maximally excited state
  MPS exc(L,D,d*d);
  string mpsfile1=getMPSfilename(L,t,U,0); // basename of MPS files
  string tmpMPSfile1;
  if(findTmpFile(D,mpsdir,mpsfile1,tmpMPSfile1)){// found some tmp or smaller D file
    //    cout<<"Recovered initial guess for exc from file "<<tmpMPSfile1<<endl;
    exc.importMPS(tmpMPSfile1.data());
    if(exc.getLength()!=L){
      cout<<"Error! File "<<tmpMPSfile1<<" does not contain valid MPS "<<endl;
      exit(1);
    }
  }
  else{
    exc.setRandomState(); // the intial state, random
    //cout<<"Initialized random state for exc from scratch "<<endl;
  }
  exc.gaugeCond('R',1);
  exc.gaugeCond('L',1);

  double Emax=0;
  {
    MPO hamilMinus(L);
    hamilMinus.setOp(0,new Operator(-1.*hamil.getOp(0).getFullData()),true);
    for(int k=1;k<L;k++){
      hamilMinus.setOp(k,&hamil.getOp(k),false);
    }
    contractor.findGroundState(hamilMinus,D,&Emax,exc);
    Emax=-Emax;

    tmpMPSfile1=mpsFileName(mpsdir,mpsfile1,D);
    exc.exportMPS(tmpMPSfile1.data());
    //cout<<"Max excited state saved to "<<tmpMPSfile1<<endl;
  }
  //cout<<"Found max Exc of H with E="<<Emax<<endl;
  complex_t nFexc=contractor.contract(exc,nrFerm,exc);
  //cout<<"Computing fermion number in this state: "<<nFexc<<endl;
  //cout<<"And the corresponding variance is: "<<contractor.contract2(nrFerm,exc)-nFexc*nFexc<<endl;
  
  complex_t totZexc=contractor.contract(exc,totSz,exc);
  //cout<<"Computing total Sz in this state: "<<totZexc<<endl;
  //cout<<"And the corresponding variance is: "<<contractor.contract2(totSz,exc)-totZexc*totZexc<<endl;
  cout<<"*** maxExc: Emax="<<Emax<<", Nf="<<nFexc<<", Sz="<<totZexc<<endl;

  
  // // For debugging. Seems ok (vs FermiHubbardHamiltonian.m in matlab) 
  //   if(L<=6){
  //    hamil.exportForMatlab("fermihubHmpo.m");
  //    totSz.exportForMatlab("fermihubSzmpo.m");
  //    nrFerm.exportForMatlab("fermihubNfmpo.m");
  //  }

  // Now check energy of a given initial state
  MPS psi(L,1,d*d);
  {    
    mwArray A00(Indices(d*d,1,1));
    mwArray A11(Indices(d*d,1,1));
    A00.setElement(ONE_c,Indices(0,0,0));
    A11.setElement(ONE_c,Indices(d*d-1,0,0));
    for(int k=0;k<L;k++){
      psi.replaceSite(k,A00);
    }
    cout<<"Z+ (all full state) has: "
	<<"E="<<contractor.contract(psi,hamil,psi)
	<<", Nf="<<contractor.contract(psi,nrFerm,psi)
	<<", Sz="<<contractor.contract(psi,totSz,psi)<<endl;

    for(int fr=2;fr<=L;fr++){
      for(int k=0;k<L;k++){
	if(k%fr==0)
	  psi.replaceSite(k,A00);
	else
	  psi.replaceSite(k,A11);
      }
      cout<<"State with 1/"<<fr<<" singlets has: "
	<<"E="<<contractor.contract(psi,hamil,psi)
	<<", Nf="<<contractor.contract(psi,nrFerm,psi)
	<<", Sz="<<contractor.contract(psi,totSz,psi)<<endl;

    }

    for(int fr=3;fr<=L;fr++){
      for(int k=0;k<L;k++){
	if(k%fr!=0)
	  psi.replaceSite(k,A00);
	else
	  psi.replaceSite(k,A11);
      }
      cout<<"State with "<<fr-1<<"/"<<fr<<" singlets has: "
	<<"E="<<contractor.contract(psi,hamil,psi)
	<<", Nf="<<contractor.contract(psi,nrFerm,psi)
	<<", Sz="<<contractor.contract(psi,totSz,psi)<<endl;

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

const string getMPSfilename(int L,double t,double U,bool GS){
  stringstream s;
  if(GS)
    s<<"gsFH";
  else
    s<<"excFH";
  s<<"_L"<<L<<"_t"<<t<<"_U"<<U;
  return s.str();
}


