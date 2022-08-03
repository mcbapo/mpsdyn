
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"
//#include "HermitianTensorMultiplier.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "SchwingerHamiltonianSz.h"
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


/** SchwingerExcSz tries to find the vector and scalar mass gap by
    first finding the MPS approximation to the ground state, then
    finding another one for the first excited state, but in every
    step, the Hamiltonian is projected onto the Sz=0 subspace.  A
    certain number of (MPS approximated) eigenstates is computed, and
    for each of them the expectation value of the operator SR
    (translation of one to the right times sigmax rotation) is
    calculated, together with the momentum operator.  When the scalar
    candidate is located (value of SR, and- maybe- momentum), a couple
    more states are computed (in some cases we could find a higher
    state earlier) and the program stops. 

    \todo The scalar candidate can saved to disk, and the results,
    including the name of the file where the MPS is saved, should be
    written to a special output file, which could be used by another
    program to start the search on an already approximated state.
    This could be run as a first search, with low D, to locate the
    approximate scalar candidate, and feed the corresponding energy to
    a second program, which used the findClosestEigenstate methods.

    Instead of targetting a fixed number of eigenstates, it stops
    three states after identifying the candidate scalar, via the sign
    of the real part of <SR> (>0) and its modulus (close to that of
    the GS) and the value (small) of the momentum operator P2.

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
  //  int D0 = props.getIntProperty("D0");
  //  int nLev = props.getIntProperty("nLev");
  directory=props.getProperty("outputdir");
  const string outfile=directory+"/"+props.getProperty("outfile");
  int app=props.getIntProperty("append");
  double tol=props.getDoubleProperty("tol");
  // int unif=1; // All same D 
  // int unif=props.getIntProperty("unif");
  int D=props.getIntProperty("Dcommon");
  double zpen=props.getDoubleProperty("penalty");
  double offset=props.getDoubleProperty("offset");

  const string filestate=props.getProperty("mpsfile");
  if(filestate.empty()){
    cout<<"No file specified to save the scalar candidate => not saving it!"<<endl;
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

  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);
  int d=2;
  MPS gs(L,D,d);
  gs.setRandomState(); // the intial state, random
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
   
  // First: put H and find the GS
  int Dop=5; // bond dimension of Schwinger H
  double lambda=0.;
  SchwingerHamiltonianSz hSch(L,mu,x,L0,alpha,zpen,offset,nu,gweight);
  //  SchwingerHamiltonian hSch(L,mu,x,L0,alpha,offset,nu,gweight);
  const MPO& hamil=hSch.getHMPO();
  // To check!!
  //  hamil.exportForMatlab("schwHamil.m");
  //  exit(1);
  contractor.findGroundState(hamil,D,&lambda,gs);
  cout<<"Ground state found with eigenvalue "<<lambda<<endl;
  *out<<"% L="<<L<<endl;
  *out<<"% mg="<<mg<<endl;
  *out<<"% x="<<x<<endl;
  *out<<"% alpha="<<alpha<<endl;
  *out<<"% D="<<D<<endl;
  *out<<"% penalty="<<zpen<<endl;
  *out<<"% offset="<<offset<<endl;
  *out<<"% gweight="<<gweight<<endl;
  *out<<setprecision(10);
  *out<<"% E0="<<lambda-offset<<endl;

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

 // try to compute the expectation value of the condensate, too
  complex_t fC=contractor.contract(gs,condMPO,gs);
  *out<<"% <condensate>="<<real(fC)<<endl;


  vector<MPS*> levels;
  levels.push_back(&gs); 
  // Some values I want to keep
  complex_t P20,SR0; // ground state reference
  complex_t P2s,SRs; double Es; // best scalar so far
  int l=0; // nr of level
  bool found=0; // whether I found the state I looked for
  int lS=0; // index of the scalar state found
  // For each excited state, repeat the following procedure
  while(!found||l<=lS+2){
  //  for(int l=0;l<=nLev;l++){
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
      if(real(SRval)>0)
	if(abs(phase(SRval))<M_PIl/4.){
	  // Compare to SR and P2 of GS
	  	  if(abs(real(SRval)-real(SR0))/real(SR0)<.2){// I am asking that it is very close to the GS value
		    //if(abs(real(SRval))>.4){ // I am asking  a minimum value
	    if(!found){
	      found=true;lS=l;
	      Es=real(Hval-zpen*Sz2val)-offset;P2s=P2val;SRs=SRval;
	      cout<<"** FIRST Candidate found, at l="<<lS<<" with properties: "<<endl; 
	      cout<<"** E="<<Es<<", P2="<<P2s<<", SR="<<SRs<<endl;
	      if(!filestate.empty()){
		cout<<"** Saving the state to file "<<filestate<<endl;
		exc.exportMPS(filestate.data());
	      }
	    }
	    else{
	      // Compare with previous candidate
	      if(real(SRval)>real(SRs)){// Accepting it
		lS=l;
		Es=real(Hval-zpen*Sz2val)-offset;P2s=P2val;SRs=SRval;
		cout<<"** NEW Candidate found, at l="<<lS<<" with properties: "<<endl; 
		cout<<"** E="<<Es<<", P2="<<P2s<<", SR="<<SRs<<endl;
		if(!filestate.empty()){
		  cout<<"** Saving the state to file "<<filestate<<endl;
		  exc.exportMPS(filestate.data());
		}
	      }
	    }
	  }
	  else{ // right phase, but seems very small, so it might not be the good one, or maybe the D is too small
	    if(!found){ // This could be the first candidate

	    }
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


