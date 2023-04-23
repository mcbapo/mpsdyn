#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "SpinMPO.h"
#include "Properties.h"
#include "misc.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

#include "IsingHamiltonian.h"

using namespace shrt;

/** Construct the initial state */
void constructInitialRho(int L,int d,int initSt,MPS& rho0);

/** Product of all identities, but sigz on site pos */
void constructLocal(int L,int pos,int d,MPS& result);

void computeLocalSigZ(const MPS& rho,int d,ofstream* out);
void computeEntropies(const MPS& rho,ofstream* out);

/** Prepare a name for the mps correspondng to step cnt */
const string mpsFileName(const string& mpsdir,const string& baseName,int cnt);

/** Save a tmp mps file, replacing the previous one */
void tmpSave(MPS& rho,const string& baseDir,const string& baseName,
	     int cnt,string& tmpFile);

/** Check for tmp files */
int findTmpFile(int M,const string& mpsdir,const string& baseName,string& mpsfile);

/**
   longTimeAverage (compiled as make longTav) evolves the time
   averaged density matrix, to find the long time limit. 
   The intermediate times are not supposed to be accurate, but the
   infinite time limit should converge to the true diagonal 
   ensemble.

   Parameters are provided as a Properties file, where, in particular,
   the following should be included
   \param <L> (int) number of sites
   \param <delta> (double) step width
   \param <J> (double) parameter \f$J\f$ of the Hamiltonian 
   \param <g> (double) parameter \f$g\f$ of the Hamiltonian 
   \param <h> (double) parameter \f$h\f$ of the Hamiltonian
   \param <Gamma> (double) weight of the double commutator term
   \param <initSt> (int) which initial rho to consider. Options are:
   1 - |X+>
   2 - |Y+>
   3 - |Z+>
   \param <delta> (double) width of the time step used
   \param <M> (int) max number of steps
   \param <rate> (int) frequency (steps) to write results
   \param <D> (int) maximum bond dimension
   \param <Dpsi> (int) [optional] bond dimension used for the pure 
                     state (default D)
   \param <outfname> (char*) name of the output file for the results
                     (will be replaced if it exists)
   \param <mpsdir> (char*) directory where to store the tmp (or final)
                     state, so that I can continue the evolution
   \param <mpsfile> (char*) base of the filename where to save or from 
                     which to read (only if newInstance is false) the
		     tmp MPS. The full name will be changed to include 
		     the time step it corresponds to. The basis of
		     thename should be such that we do not mix
		     evolutions with different parameters (e.g. cahnging step width).
   \param <savingFreq> (int) frequency (steps) to save the tmp MPS file
   \param <newInstance> (bool) whether to start a new evolution from
                         the inhomogeneous polarization distribution
			 (if not, the MPS in mpsfile is used)

 */

int d=2;

int main(int argc,const char* argv[]){
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
  double J_=props.getDoubleProperty("J");
  double h_=props.getDoubleProperty("h");
  double g_=props.getDoubleProperty("g");
  int initSt=props.getIntProperty("initSt");
  double delta=props.getDoubleProperty("delta");
  double Gamma=props.getDoubleProperty("Gamma");
  int M=props.getIntProperty("M");
  int rate=props.getIntProperty("rate");
  int D=props.getIntProperty("D");
  string outfname=props.getProperty("outputfile");
  string mpsdir=props.getProperty("mpsdir");
  string mpsfile=props.getProperty("mpsfile");
  int savingFreq=props.getIntProperty("savingFreq");
  int _aux_=props.getIntProperty("newInstance");
  bool newInstance=_aux_>0;
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-5;    
  bool useHam=true;_aux_=props.getIntProperty("noHam");
  if(_aux_>0) useHam=false;

  // For output
  ofstream* out;
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  IsingHamiltonian hamH(L,d,J_,g_,h_);
  //hamH.registerTimeDependence();
  cout<<"Created the Hamiltonian"<<endl;

  //Now I need to prepare the evolution operators. Since [H,.] and
  //[H,[H,.]] commute, I will use the normal ones for the first
  //(unitary) term, and Taylor for the second

  // As first order approach, I just take the MPO for [H,.] (D=4), 
  // and square it, then take Taylor expansion for the exponential
  // A more efficient strategy exists, by designing the MPO for H^2
  // and do terms separately

  // 1) MPO for the commutator
  MPO mpoComm(L),mpoComm2(L);
  hamH.getCommutatorMPO(mpoComm,0.);
  cout<<"Created the commutator MPO"<<endl;
  // And I need to square it
  const MPO* mpos[2]={&mpoComm,&mpoComm};
  MPO::join(2,mpos,mpoComm2);
  cout<<"Created the commutator squared MPO"<<endl;

  // 2) Basic terms for the exponential: just one
  MPO Uevol(L); // this is the double one
  if(useHam){
    MPO _Uevol_1(L); // held just temprarily
    hamH.getUMPO(_Uevol_1,-delta*I_c,0); // the single layer normal evolution
    // But I need to double it
    doubleMPO(_Uevol_1,Uevol,true);
    cout<<"Created the double Uevol"<<endl;
  }

  // Prepare the observables to compute
  // The identity, to compute the trace
  MPS Id(L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++)
    Id.setA(k,id2);
  // The Hamiltonian to get the energy  
  MPS hmps(L,1,d*d);
  const MPO& hamil=hamH.getHMPO();
  MPSfromMPO(hamil,hmps);

  // Sum of sigmaz/sigmax operators
  MPS SxSum(L,2,d*d),SzSum(L,2,d*d);
  {
    MPO Saux(L);
    SpinMPO::getSxMPO(L,d,Saux);
    MPSfromMPO(Saux,SxSum);
    SpinMPO::getSzMPO(L,d,Saux);
    MPSfromMPO(Saux,SzSum);
  }
  // And I will also need the sigmaz/sigmax operator on central site, which I could compute ony by one, but I assemble it here
  MPS SxCenter(Id);
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  sigX.reshape(Indices(d*d,1,1));
  SxCenter.setA(L/2,sigX);sigX.reshape(Indices(d,d));

  // Initial state rho0
  MPS rho0(L,1,d*d);
  // Now if we start from scratch, initialize the MPS. Otherwise, read
  // the tmp file
  string tmpMPSfile;
  int cnt=0;
  double time=0.;
  if(!newInstance){
    cnt=findTmpFile(M,mpsdir,mpsfile,tmpMPSfile);
    if(cnt>0)
      cout<<"Recovered state for cnt="<<cnt<<" at file "<<tmpMPSfile<<endl;
  }
  if(newInstance||cnt==0){ // init from scratch
    constructInitialRho(L,d,initSt,rho0);    
    cout<<"Created state from scratch"<<endl;
    out=new ofstream(outfname.data());
    // Write header
    *out<<"% cnt\t time\t D\t tr(rho)\t <E>\t <Sxav>\t <Szav>\t <Sx(L/2)>"<<endl;
    out->close();
    delete out;
  }
  else{ 
    rho0.importMPS(tmpMPSfile.data());
  }
  rho0.increaseBondDimension(D);
  rho0.gaugeCond('R',1);
  rho0.gaugeCond('L',1);
  time=cnt*delta;


  // And now iterate the evolution
  cout<<"=== STARTING TIME EVOLUTION FROM STEP "<<cnt<<" ==="<<endl;
  complex_t trR=contractor.contract(rho0,Id);
  complex_t Eev=contractor.contract(rho0,hmps);
  complex_t Sxev=contractor.contract(rho0,SxSum);
  complex_t Szev=contractor.contract(rho0,SzSum);
  complex_t SxCev=contractor.contract(rho0,SxCenter);
  out=new ofstream(outfname.data(),ios::app);
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D
      <<"\t"<<real(trR)<<"\t"<<imag(trR)
      <<"\t"<<real(Eev)<<"\t"<<imag(Eev)
      <<"\t"<<real(Sxev)<<"\t"<<imag(Sxev)
      <<"\t"<<real(Szev)<<"\t"<<imag(Szev)
      <<"\t"<<real(SxCev)<<"\t"<<imag(SxCev)
      <<endl;
  out->close();


  int stepsLeft=M-cnt;
  while(stepsLeft>0){
    int nrToApply=min(stepsLeft,rate);
    for(int k=0;k<nrToApply;k++){
      MPS aux(rho0);
      if(useHam)
	contractor.optimize(Uevol,rho0,aux,D);
      contractor.approximateTaylorExponential(mpoComm2,aux,rho0,Gamma,1,D,true);
      cnt++;stepsLeft--;time+=delta;
      rho0.gaugeCond('R',1);
    }
    // Now save results
    trR=contractor.contract(rho0,Id);
    Eev=contractor.contract(rho0,hmps)/trR;
    Sxev=contractor.contract(rho0,SxSum)/trR;
    Szev=contractor.contract(rho0,SzSum)/trR;
    SxCev=contractor.contract(rho0,SxCenter)/trR;
    out=new ofstream(outfname.data(),ios::app);
    *out<<setprecision(12);
    *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D
	<<"\t"<<real(trR)<<"\t"<<imag(trR)
	<<"\t"<<real(Eev)<<"\t"<<imag(Eev)
	<<"\t"<<real(Sxev)<<"\t"<<imag(Sxev)
	<<"\t"<<real(Szev)<<"\t"<<imag(Szev)
	<<"\t"<<real(SxCev)<<"\t"<<imag(SxCev)
	<<endl;
    out->close();
    // And if necessary, store tmpfile

    if(cnt>0&&(cnt%savingFreq==0||cnt==M)){
      cout<<"Saving backup at cnt="<<cnt<<endl;
      tmpSave(rho0,mpsdir,mpsfile,cnt,tmpMPSfile);
    }
  }

}



const string mpsFileName(const string& mpsdir,const string& baseName,int cnt){
  stringstream s;
  s<<mpsdir<<"/"<<baseName;
  if(cnt>0) s<<"_"<<cnt;
  return s.str();
}

void tmpSave(MPS& rho,const string& baseDir,const string& baseName,int cnt,string& tmpFile){
  // first create a new tmpFile name
  string newTmp=mpsFileName(baseDir,baseName,cnt);
  rho.exportMPS(newTmp.data());
  // now remove the old one
  if(tmpFile.length()>0&&file_exists(tmpFile)){
    remove(tmpFile.data());
  }
  // and replace the string
  tmpFile=newTmp;
}

int findTmpFile(int M,const string& mpsdir,const string& baseName,string& mpsfile){
  int cnt=M;
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
  if(!found) return 0;
}

void constructInitialRho(int L,int d,int initSt,MPS& rho0){
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sig0=identityMatrix(d);
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray sigX(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  
  mwArray state;
  switch(initSt){
  case 1: // X+
    {
      state=.5*(sig0+sigX);
      break;
    }
  case 2: // Y+
    {
      state=.5*(sig0+sigY);
      break;
    }
  case 3: // Z+
    {
      state=.5*(sig0+sigZ);
      break;
    }
  default:
    cout<<"Error: unknown type onf intSt="<<initSt<<endl;
    exit(1);
  }
  state.reshape(Indices(d*d,1,1));

  rho0=MPS(L,1,d*d);
  for(int k=0;k<L;k++)
    rho0.setA(k,state);
}
