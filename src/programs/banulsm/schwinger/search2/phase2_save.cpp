
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"
//#include "HermitianTensorMultiplier.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "SchwingerHamiltonianSz.h"
//#include "SchwingerHamiltonianSzSR.h"
#include "SpinMPO.h"
#include <cmath>

using namespace std;
using namespace shrt;

#ifdef MACOSX
#define LOCALDIR "/Users/banulsm/Desktop/SVN/projects/mpsdyn/src/bin/"
#else
#define LOCALDIR "/home/banulsm/ALaPorra/bin/"
#endif
#define MAXLEN 120


/** SchwingerExcSzSR3 tries to find the scalar mass gap in two steps:
    1) first, it finds the MPS approximation to the ground state, then
    another one for the first excited state (vector). The energy of the 
    scalar is then estimated as ~2 mV (or a bit less, let's see)
    2) After the energy of the scalar state is estimated in that way,
    its value, E, is used to construct a new Hamiltonian, (H-E)^2,
    whose ground state is found again as an MPS. Several excited
    states are then found, from this Hamiltonian, and the first one
    with the proper label will be identified with the scalar state.

    As results, for each the GS and the excited states, the expectation
    value of the operator SR (translation of one to the right times
    sigmax rotation) is calculated, together with the momentum
    operator. They will later serve to identify the vector and scalar
    energies. 

    Receives as argument a Properties file, such as
    config/excitProp.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.

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
  int nLev = props.getIntProperty("nLev");
  directory=props.getProperty("outputdir");
  const string outfile=directory+"/"+props.getProperty("outfile");
  int app=props.getIntProperty("append");
  double tol=props.getDoubleProperty("tol");
  int D=props.getIntProperty("Dcommon");
  double zpen=props.getDoubleProperty("penalty");
  //  double lamb=props.getDoubleProperty("penaltySR");
  double offset=props.getDoubleProperty("offset");
  if(offset==-1) offset=0;
  double scale=props.getDoubleProperty("scale");
  if(scale==-1) scale=1;
  //double offset2=props.getDoubleProperty("offset2");
  //if(offset2==-1) offset2=0;
  //  double estimE=props.getDoubleProperty("estimatedS");

  const string filestate=props.getProperty("mpsfile");
  if(filestate.empty()){
    cout<<"No file specified containing a previously calculated scalar "
	<<"candidate => ABORTING!"<<endl;
    cout<<"USE phase1 first"<<endl;
    exit(1);
  }


  // Parameters I do not use here: chemical potential and weight of the gauge term
  double nu=0;
  // I might want to switch off the gauge terms
  //double gweight=props.getDoubleProperty("gweight");
  //  if(gweight==-1) gweight=1; // default: do not alter
  double gweight=1.;

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
      <<", D="<<D<<endl;

#ifdef DISKSTRG
  //  int errD=system("rm -rf /ptmp/mpq/banulsm/*");
  char tmpdir[150];
  sprintf(tmpdir,"/ptmp/mpq/banulsm/Schwinger/Sitesx%dL%dD%dmg%dXXXXXX",x,L,D,(int)(mg*1000));
  char* dirid=mkdtemp(tmpdir);
  if(dirid==0){
    cout<<"Error: couldn't apply mkdtemp "<<tmpdir<<endl;
    exit(1);
  }
  cout<<"Using FileSites in "<<tmpdir<<endl;
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

  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);
  int d=2;
  MPS gs(L,D,d);
  gs.importMPS(filestate.data());
  cout<<"Imported MPS candidate from file "<<filestate<<endl;
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
   
  // First: put H and find the GS
  double lambda=0.;
  double estimE;
  SchwingerHamiltonianSz hSchOrig(L,mu,x,L0,alpha,zpen,offset,nu,gweight);
  const MPO& hamilOrig=hSchOrig.getHMPO();
  estimE=real(contractor.contract(gs,hamilOrig,gs))-offset;
  cout<<"The MPS read had energy "<<estimE<<endl;
  
  // Now create [alpha(H-E)][alpha(H-E)] and find its lowest states
  SchwingerHamiltonianSz hSchScal(L,mu,x,L0,alpha,zpen,offset-estimE,nu,gweight,scale);
  const MPO& hamilScal=hSchScal.getHMPO();
  // And my "Hamiltonian" is the square of that!
  const MPO* ptr[2]={&hamilScal,&hamilScal};
  MPO hamil(L);
  MPO::join(2,ptr,hamil);


  //hamilOrig.exportForMatlab("schwHamil.m");
  // To check!!
  contractor.findGroundState(hamil,D,&lambda,gs);
  cout<<"Ground state of (H-estimE) found with eigenvalue "<<sqrt(lambda)/scale<<endl;
  *out<<"% L="<<L<<endl;
  *out<<"% mg="<<mg<<endl;
  *out<<"% x="<<x<<endl;
  *out<<"% alpha="<<alpha<<endl;
  *out<<"% D="<<D<<endl;
  *out<<"% penalty="<<zpen<<endl;
  *out<<"% offset="<<offset<<endl;
  *out<<"% gweight="<<gweight<<endl;
  *out<<"% scale="<<scale<<endl;
  *out<<"% estimE="<<estimE<<endl;
  *out<<setprecision(10);
  *out<<"% lambda[alpha(H-E)]^2="<<lambda<<endl;
  *out<<"% E="<<sqrt(lambda)/scale-offset+estimE<<endl;

  // try to compute the expectation value of the condensate, too
  MPO cond(L);
  hSchScal.constructCondensateMPO(cond);
  complex_t fC=contractor.contract(gs,cond,gs);
  *out<<"% <condensate>="<<real(fC)<<endl;
  
  // Observables I will compute
  // Compute the expectation value of Sz, as it commutes with H
  MPO Szmpo(L);
  SpinMPO::getSzMPO(L,d,Szmpo);
  MPO Sz2mpo(L);
  SpinMPO::getSz2MPO(L,d,Sz2mpo);
  // The momentum operator
  MPO Pmpo(L);
  hSchScal.constructMomentumMPO(Pmpo);
  const MPO* oprsP[2]={&Pmpo,&Pmpo};
  MPO P2mpo(L);
  MPO::join(2,oprsP,P2mpo);
  //  P2mpo.exportForMatlab("P2mpo.m");
  // // The momentum squared operator
  // MPO PSqmpo(L);
  // hSch.constructMomentumSquaredMPO(PSqmpo);

  // the shift by one and rotate sigmax
  MPO SRmpo(L);
  //hSchOrig.constructShiftRotateMPO(SRmpo,false);
  //SRmpo.exportForMatlab("SRDmpo.m");
  hSchScal.constructShiftRotateMPO(SRmpo);
  //SRmpo.exportForMatlab("SRmpo.m");

  // The fermionic condensate and the Gamma_alpha, Gamma_5 parameters
  MPO condMPO(L),GammaA(L),Gamma5(L);
  hSchOrig.constructCondensateMPO(condMPO);
  hSchOrig.constructGammaAlphaMPO(GammaA);
  hSchOrig.constructGamma5MPO(Gamma5);

  vector<MPS*> levels;
  //  vector<MPS*> levelsExc; // The ones wrt which I want to orthogonalize
  levels.push_back(&gs); 

  for(int l=0;l<nLev;l++){
    if(l>0){ // compute one new excitation)
      int D1=D;
      cout<<"Initializing computation of excitation nr "<<l<<" with D="<<D1<<endl;
      // Now find the first excited state
      MPS* exc=new MPS(L,D1,d);
      exc->setRandomState(); // the initial state, random
      double lambda1=0.;
      contractor.findNextExcitedState(hamil,D1,levels,&lambda1,*exc);
      if(lambda1>0){
	cout<<"Stopping iteration as I got an energy above zero!"<<endl;
	break;
      }
      levels.push_back(exc);
    }
   // Compute the expectation values
    const MPS& exc=*levels[l];
 
   complex_t Hval=contractor.contract(exc,hamilOrig,exc);
    //cout<<"Computed <H>="<<Hval<<endl;
    complex_t H2val=contractor.contract2(hamilOrig,exc);
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

  }


  out->close();
  delete out;

   // for(int z=1;z<L;z++){
   //   delete Skmpos[z];
   // }
   // delete []Skmpos;
}


