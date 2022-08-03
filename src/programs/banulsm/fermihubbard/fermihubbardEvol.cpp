
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

/** fermihubbardEvol runs the evolution with the
    Hamiltonian \ref <FermiHubbardHamiltonian>, starting from a state with a given fraction of 
    fully occupied sites, and computes the Loschmidt echo.

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
  double dt=props.getDoubleProperty("dt");
  double timeF=props.getDoubleProperty("time"); // final time
  vector<int> posSing=props.getVectorIntProperty("posSinglet");
  vector<int> posNoSing=props.getVectorIntProperty("posNoSinglet");
  int Nsingl;
  if(posSing.size()>0) Nsingl=posSing.size();
  else
    if(posNoSing.size()>0) Nsingl=L-posNoSing.size();
    else{
      cout<<"Error: no position given for the singlets"<<endl;
      exit(1);
    }
  // if(posSing.size()==0){
  //   cout<<"Positions for the singlets not explicit, using the complementary"<<endl;

  //   Nsingl=posSing.size();    

  // }
  
  int stepBlock=props.getIntProperty("stepBlock"); // saving frequency
  if(stepBlock<=0) stepBlock=1;
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  string outfname=props.getProperty("output");
  //string mpsdir=props.getProperty("mpsdir"); // where are the MPS
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  
  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)

  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    *out<<"%time\t real(<Psi0|Psi(t)>)\t imag(<Psi0|Psi(t)>)\t E\t Nf\t Sz"<<endl;;
    out->close();
    delete out;
  }
  
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
  FermiHubbardHamiltonian ham(L,t,U);
  //  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=ham.getHMPO();
  MPO nrFerm(L);
  ham.getFermionNrMPO(nrFerm);
  MPO totSz(L);
  ham.getTotalSzMPO(totSz);

  MPO evolUe(L),evolUo(L),evolUe_2(L);
  ham.getExponentialMPOeven(evolUe,-dt*I_c);
  ham.getExponentialMPOeven(evolUe_2,-dt*.5*I_c);
  ham.getExponentialMPOodd(evolUo,-dt*I_c);

  // how many actual singlets do I need to put?
  //  int Nsingl=floor(L*fr);
  cout<<"Init state with "<<Nsingl<<" singlets"<<endl;
  //int NsinglFr=floor(L/Nsingl);
  MPS psi(L,1,d*d);
  mwArray A00(Indices(d*d,1,1));
  mwArray A11(Indices(d*d,1,1));
  A00.setElement(ONE_c,Indices(0,0,0));
  A11.setElement(ONE_c,Indices(d*d-1,0,0));
  if(posSing.size()>0){
    for(int k=0;k<L;k++){
      psi.replaceSite(k,A11);
    }
    for(int k=0;k<Nsingl;k++){
      cout<<"Setting singlet in pos "<<posSing[k]<<endl;
      psi.replaceSite(posSing[k],A00);
    }
  }
  else{
    for(int k=0;k<L;k++){
      psi.replaceSite(k,A00);
    }
    for(int k=0;k<posNoSing.size();k++){
      cout<<"Removing singlet from pos "<<posNoSing[k]<<endl;
      psi.replaceSite(posNoSing[k],A11);
    }
  }
  psi.gaugeCond('R',1);
  complex_t E=contractor.contract(psi,hamil,psi);
  complex_t Nf=contractor.contract(psi,nrFerm,psi);
  complex_t Sz=contractor.contract(psi,totSz,psi);
  cout<<"Initial state has: "
      <<"E="<<E<<", Nf="<<Nf<<", Sz="<<Sz<<endl;
  MPS psi0(psi); // for reference
  int M=ceil(timeF/dt); // nr of Trotter steps
  double time=0;
  int remainingSteps=M;

  // write the first line in the file
  out=new ofstream(outfname.data(),ios::app);*out<<setprecision(10);
  *out<<0<<"\t"<<1<<"\t"<<0<<"\t"<<real(E)<<"\t"<<real(Nf)<<"\t"<<real(Sz)<<endl;
  out->close();
  delete out;
  int cnt=0;
  while(remainingSteps>0){
    int nrToApply=min(stepBlock,remainingSteps);
    MPS aux(psi);
    contractor.optimize(evolUe_2,aux,psi,D);
    for(int k=0;k<nrToApply;k++){
      contractor.optimize(evolUo,psi,aux,D);
      //cout<<"Step "<<++cnt<<endl;
      if(k<nrToApply-1)
	contractor.optimize(evolUe,aux,psi,D);
      else
	contractor.optimize(evolUe_2,aux,psi,D);
    }
    psi.gaugeCond('R',1);
    remainingSteps-=nrToApply;
    time+=dt*nrToApply;
    complex_t overl=contractor.contract(psi,psi0);
    E=contractor.contract(psi,hamil,psi);
    Nf=contractor.contract(psi,nrFerm,psi);
    Sz=contractor.contract(psi,totSz,psi);

    out=new ofstream(outfname.data(),ios::app);*out<<setprecision(10);
    *out<<time<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t"<<real(E)<<"\t"<<real(Nf)<<"\t"<<real(Sz)<<endl;
    out->close();
    delete out;
    //cout<<time<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t"<<real(E)<<"\t"<<real(Nf)<<"\t"<<real(Sz)<<endl;
  }



  // // For debugging. Seems ok (vs FermiHubbardHamiltonian.m in matlab) 
  //   if(L<=6){
  //    hamil.exportForMatlab("fermihubHmpo.m");
  //    totSz.exportForMatlab("fermihubSzmpo.m");
  //    nrFerm.exportForMatlab("fermihubNfmpo.m");
  //  }

  
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


