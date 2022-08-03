
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

/** Check if the filename exists */
bool file_exists(const string filename);


/** This program tries to find the vector and scalar mass gap by
    first finding the MPS approximation to the ground state, then
    finding another one for the first excited state, but in every
    step, the Hamiltonian is projected onto the Sz=0 subspace.  A
    certain number of (MPS approximated) eigenstates is computed, and
    for each of them the expectation value of the operator SR
    (translation of one to the right times sigmax rotation) is
    calculated, together with the momentum operator.  When the scalar
    candidate is located (value of SR, and- maybe- momentum), the
    program stops.  Here we save all the MPS computed to disk, to a
    directory whose name is uniquely related to the instance , so that
    at the beginning, the directory (mpsdir in properties file) is
    checked, and if MPS files exst, they are read and the program
    starts using them already.

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
  string initmpsdir=props.getProperty("initmpsdir"); // default none
  int initD=props.getIntProperty("initD");
  int stretchmps=props.getIntProperty("stretchmps");
  if(stretchmps=-1) stretchmps=0; 
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

  if(!app){
    out=new ofstream(outfile.data());
  }
  else{
    out=new ofstream(outfile.data(),ios::app);
    *out<<"%%%%%%%%%%%%%%%%%%%"<<endl;
  }
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

  // Now I already have read all existing MPS
  cout<<"Recovered "<<lval<<" already computed MPS levels from disk."<<endl;
  if(lval==0){
    nextLev=new MPS(L,D,d);
    const string initfilename=mpsfilename(initL,initD,x,mg,lval,initmpsdir);
    cout<<"Will search initial GS in "<<initfilename<<endl;
    if(initmpsdir.empty()||!file_exists(initfilename)) 
      nextLev->setRandomState();
    else{
      if(initL<L){
	MPS auxMPS(initL,initD,d);
	auxMPS.importMPS(initfilename.data());
	nextLev->stretchMPS(L,auxMPS);
	cout<<"Stretched initial state (l="<<lval<<")from "<<initfilename<<endl;
      }
      else{
	nextLev->importMPS(initfilename.data());
	cout<<"Recovered initial state for level "<<lval<<" from file "
	    <<initfilename<<endl;
      }
      if(initD<D&&noise!=0.)
	nextLev->increaseBondDimensionWithNoise(D,noise);
    }
    nextLev->gaugeCond('R',1);
    nextLev->gaugeCond('L',1);
    contractor.findGroundState(hamil,D,&lambda,*nextLev);
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
    // try to compute the expectation value of the condensate, too
    complex_t fC=contractor.contract(*nextLev,condMPO,*nextLev);
    *out<<"% <condensate>="<<real(fC)<<endl;
    levels.push_back(nextLev); 
    const string filename=mpsfilename(L,D,x,mg,lval,mpsdir);
    nextLev->exportMPS(filename.data());
    lval++;
  }

  // Some values I want to keep
  complex_t P20,SR0; // ground state reference
  complex_t P2s,SRs; double Es; // best scalar so far
  int l=0; // nr of level
  bool found=0; // whether I found the state I looked for
  int lS=0; // index of the scalar state found
  // For each excited state, repeat the following procedure
  while(!found&&(!stopAtVec||l<=1)){
  //  for(int l=0;l<=nLev;l++){
    if(l>=lval){ // compute one new excitation)
      cout<<"Initializing computation of excitation nr "<<l<<" with D="<<D<<endl;
      // Now find the first excited state
      MPS* exc=new MPS(L,D,d);
      const string initfilename=mpsfilename(initL,initD,x,mg,l,initmpsdir);
      if(initmpsdir.empty()||!file_exists(initfilename)) 
	exc->setRandomState(); // the initial state, random
      else{
	if(initL<L){
	  MPS auxMPS(initL,initD,d);
	  auxMPS.importMPS(initfilename.data());
	  exc->stretchMPS(L,auxMPS);
	  cout<<"Stretched initial state (l="<<lval<<")from "<<initfilename<<endl;
	}
	else{	
	  exc->importMPS(initfilename.data());
	  cout<<"Recovered initial state for level "<<l<<" from file "
	      <<initfilename<<endl;
	}
	if(initD<D&&noise!=0.)
	  exc->increaseBondDimensionWithNoise(D,noise);
      }
      double lambda1=0.;
      contractor.findNextExcitedState(hamil,D,levels,&lambda1,*exc);
      if(lambda1>0){
	cout<<"Stopping iteration as I got an energy above zero!"<<endl;
	break;
      }
      levels.push_back(exc);
      const string filename=mpsfilename(L,D,x,mg,l,mpsdir);
      exc->exportMPS(filename.data());
    }
    // Compute the expectation values
    const MPS& exc=*levels[l];
    //MPS reflected(exc);
    //hSch.applyReflectionRotation(exc,reflected);
    complex_t Hval=contractor.contract(exc,hamil,exc);
    //cout<<"Computed <H>="<<Hval<<endl;
    complex_t H2val=contractor.contract2(hamil,exc);
    //cout<<"Computed <H2>="<<H2val<<endl;
    complex_t normN=contractor.contract(exc,exc);
    //cout<<"Computed norm="<<normN<<endl;
    complex_t Szval=contractor.contract(exc,Szmpo,exc);
    //cout<<"Computed <Sz>="<<Szval<<endl;
    complex_t Sz2val=contractor.contract(exc,Sz2mpo,exc);
    //cout<<"Computed <Sz^2>="<<Sz2val<<endl;
    //    complex_t Pval=contractor.contract(exc,Pmpo,exc);
    complex_t P2val=contractor.contract(exc,P2mpo,exc);
    //cout<<"Computed <P^2>="<<P2val<<endl;
    //complex_t PSqval=contractor.contract(exc,PSqmpo,exc);
    complex_t SRval=contractor.contract(exc,SRmpo,exc);
    complex_t fC=contractor.contract(exc,condMPO,exc); // fermion condensate
    complex_t gammaA=contractor.contract(exc,GammaA,exc); // Gamma_alpha
    complex_t gamma5=contractor.contract(exc,Gamma5,exc); // Gamma_5
 
   // Retain the GS value for reference
    if(l==0){ 
      P20=P2val;SR0=SRval;
    }
    //cout<<"Computed <SR>="<<SRval<<endl;
    //complex_t Parity=contractor.contract(exc,reflected);
    //    complex_t SCEcoeff=contractor.contract(exc,sceVac);
    // complex_t SCEcoeff2=contractor.contract(exc,sceVacC);

    *out<<l<<"\t"<<exc.getBond()<<"\t";
    *out<<real(normN)<<"\t";
    *out<<real(H2val)<<"\t"; //<<imag(H2)<<"\t";
    *out<<real(Hval-zpen*Sz2val)-offset<<"\t"; //<<imag(Hval)<<"\t";
    //*out<<real(H2/normN-Hval*Hval/(normN*normN))<<"\t";
    *out<<real(Szval)<<"\t";
    //*out<<real(Pval)<<"\t";
    *out<<real(P2val)<<"\t";
    *out<<real(SRval)<<"\t"<<imag(SRval)<<"\t";
    *out<<real(fC)<<"\t";
    *out<<real(gammaA)<<"\t";
    *out<<real(gamma5)<<"\t";
    // for(int s=0;s<=nLev;s++)
    //   if(s<l)
    // 	*out<<contractor.contract(exc,*levels[s])<<"\t";
    //   else if(s==l)
    // 	*out<<1<<"\t";
    //   else
    // 	*out<<0<<"\t";
    *out<<endl;

    // Decide whether I have found the state
    if(l>0) // GS does not count
      if(real(SRval)>0){
	found=true;lS=l;
	Es=real(Hval-zpen*Sz2val)-offset;P2s=P2val;SRs=SRval;
	cout<<"** FIRST Candidate found, at l="<<lS<<" with properties: "<<endl; 
	cout<<"** E="<<Es<<", P2="<<P2s<<", SR="<<SRs<<endl;
	if(abs(phase(SRval))>M_PIl/4.){
	  cout<<"WARNING!! THis might not be it!!!"<<endl;
	}
      }
    l++;
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

bool file_exists(const string filename)
{
  ifstream ifile(filename.data());
  //cout<<"file_exists("<<filename<<")="<<ifile<<endl;
  return ifile;
}
