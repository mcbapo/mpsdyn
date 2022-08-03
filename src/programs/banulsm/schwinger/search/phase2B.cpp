
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
    we continue trying to compute excitations closest to the estimated value (?), E.

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
  int D0 = props.getIntProperty("D0");
  int nLev = props.getIntProperty("nLev");
  directory=props.getProperty("outputdir");
  const string outfile=directory+"/"+props.getProperty("outfile");
  int app=props.getIntProperty("append");
  double tol=props.getDoubleProperty("tol");
  int unif=props.getIntProperty("unif");
  int Dall=props.getIntProperty("Dcommon");
  double zpen=props.getDoubleProperty("penalty");
  //  double lamb=props.getDoubleProperty("penaltySR");
  double offset=props.getDoubleProperty("offset");
  if(offset==-1) offset=0;
  double offset2=props.getDoubleProperty("offset2");
  if(offset2==-1) offset2=0;
  double estimE=props.getDoubleProperty("estimatedS");

  int D=Dall;

  // Parameters I do not use here: chemical potential and weight of the gauge term
  double nu=0;
  // I might want to switch off the gauge terms
  double gweight=props.getDoubleProperty("gweight");
  if(gweight==-1) gweight=1; // default: do not alter

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
      <<", D="<<D<<endl;

#ifdef DISKSTRG
  //  int errD=system("rm -rf /ptmp/mpq/banulsm/*");
  char tmpdir[150];
  sprintf(tmpdir,"/ptmp/mpq/banulsm/Schwinger/Sitesx%dL%dD%dmg%dXXXXXX",x,L,D0,(int)(mg*1000));
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
  gs.setRandomState(); // the intial state, random
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
   
  // First: put H and find the GS
  double lambda=0.;
  SchwingerHamiltonianSz hSchOrig(L,mu,x,L0,alpha,zpen,offset,nu,gweight);
  const MPO& hamilOrig=hSchOrig.getHMPO();
  //hamilOrig.exportForMatlab("schwHamil.m");
  // To check!!
  contractor.findGroundState(hamilOrig,D,&lambda,gs);
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
  
   // try to compute the expectation value of the condensate, too
    MPO cond(L);
    hSchOrig.constructCondensateMPO(cond);
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
  hSchOrig.constructMomentumMPO(Pmpo);
  const MPO* oprsP[2]={&Pmpo,&Pmpo};
  MPO P2mpo(L);
  MPO::join(2,oprsP,P2mpo);
  P2mpo.exportForMatlab("P2mpo.m");
  // // The momentum squared operator
  // MPO PSqmpo(L);
  // hSch.constructMomentumSquaredMPO(PSqmpo);

  // the shift by one and rotate sigmax
  MPO SRmpo(L);
  //hSchOrig.constructShiftRotateMPO(SRmpo,false);
  //SRmpo.exportForMatlab("SRDmpo.m");
  hSchOrig.constructShiftRotateMPO(SRmpo);
  SRmpo.exportForMatlab("SRmpo.m");

  // // the shift by one and rotate sigmax, but to left
  // MPO SRmpoL(L);
  // hSch.constructShiftRotateMPO(SRmpoL,false);
  // SRmpo.exportForMatlab("SRmpoL.m");

  // // the shift by one and rotate sigmax, to the right but cyclic
  // MPO TRmpo(L);
  // hSch.constructCyclicTranslationRotateMPO(TRmpo,true);
  // MPO TRmpoL(L); // idem to left
  // hSch.constructCyclicTranslationRotateMPO(TRmpoL,false);
  // TRmpo.exportForMatlab("TRmpo.m");



  // // A special (product) state: the vacuum for x=0, to check the coefficients
  // MPS sceVac(L,1,d);sceVac.setProductState(p_zero);
  // MPS sceVacC(L,1,d);sceVacC.setProductState(p_one); // the complementary one
  // for(int k=0;k<L;k++){
  //   if(k%2!=0)
  //     sceVac.setA(k,sceVacC.getA(k).getA());
  //   else
  //     sceVacC.setA(k,sceVac.getA(k).getA());
  // }


  vector<MPS*> levels;
  vector<MPS*> levelsExc; // The ones wrt which I want to orthogonalize
  levels.push_back(&gs); 
  // Find also the first excited state, which will be the vector (Q: find a second one, too?)
  cout<<"Initializing computation of first excitation with D="<<D1<<endl;
  // Now find the first excited state
  MPS* vec=new MPS(L,D,d);
  vec->setRandomState(); // the initial state, random
  double lambda1=0.;
  contractor.findNextExcitedState(hamilOrig,D,levels,&lambda1,*vec);
  cout<<"First excited (vector) state found with eigenvalue "<<lambda1<<"=> mass gap="<<lambda1-lambda<<endl;
  levels.push_back(vec);
  cout<<"Substracting estimation of the scalar energy "<<estimE<<" (gap~"<<estimE-lambda<<")"<<endl;
  // NEw auxiliary Hamiltonian
  SchwingerHamiltonianSz hSchOrigAux(L,mu,x,L0,alpha,zpen,offset-estimE,nu,gweight);
  const MPO& hamilOrigAux=hSchOrigAux.getHMPO();
  cout<<"Value of H-E: GS="<<contractor.contract(gs,hamilOrigAux,gs)
      <<";\n\t V="<<contractor.contract(*vec,hamilOrigAux,*vec)<<endl;
  // Find the eigenstate closest to zero (offset)
  // The new "Hamiltonian"
  const MPO* mpoPtr[2];mpoPtr[0]=&hamilOrigAux;mpoPtr[1]=&hamilOrigAux;
  MPO hamil2(L);
  MPO::join(2,mpoPtr,hamil2);
  MPS* gs2=new MPS(L,D,d);
  gs2->setRandomState();
  contractor.setEigenSolver(primme_JDQR); // seems to work better for squared things
  contractor.findGroundState(hamil2,D,&lambda1,*gs2,offset2);
  levels.push_back(gs2);
  levelsExc.push_back(gs2);
  // For a certain number of excited states repeat the following procedure
  for(int l=0;l<=nLev+1;l++){
    if(l>2){ // compute one new excitation of hamil2
      cout<<"Initializing computation of excitation nr "<<l-2<<" of the H2 with D="<<D1<<endl;
      // Now find the first excited state
      MPS* exc=new MPS(L,D,d);
      exc->setRandomState(); // the initial state, random
      contractor.findNextExcitedState(hamil2,D,levelsExc,&lambda1,*exc,offset2);
      //contractor.findNextExcitedState(hamil2,D1,levels,&lambda1,*exc,offset2);
      //complex_t SRval=contractor.contract(*exc,SRmpo,*exc);
      //complex_t estimE=contractor.contract(*exc,hamilOrig,*exc);
      //contractor.findNextExcitedState(hamil2,D1,levels,&lambda1,*exc);
      levels.push_back(exc);
      levelsExc.push_back(exc);
    }
    // Compute the expectation values
    const MPS& exc=*levels[l];
    //complex_t Hval=contractor.contract(exc,hamil,exc);
    complex_t Hval=contractor.contract(exc,hamilOrig,exc);
    //cout<<"Computed <H>="<<Hval<<endl;
    complex_t H2val=contractor.contract2(hamilOrig,exc);
    //    complex_t H2val=contractor.contract2(hamil,exc);
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


    //cout<<"Computed <SR>="<<SRval<<endl;
    //complex_t Parity=contractor.contract(exc,reflected);
    //    complex_t SCEcoeff=contractor.contract(exc,sceVac);
    // complex_t SCEcoeff2=contractor.contract(exc,sceVacC);

    *out<<l<<"\t"<<exc.getBond()<<"\t";
    *out<<real(normN)<<"\t";
    *out<<real(H2val)<<"\t"; //<<imag(H2)<<"\t";
    //*out<<real(Hval-zpen*Sz2val)-lamb*(1-real(SRval))<<"\t"; //<<imag(Hval)<<"\t";
    *out<<real(Hval)<<"\t"; //<<imag(Hval)<<"\t";
    //*out<<real(H2/normN-Hval*Hval/(normN*normN))<<"\t";
    *out<<real(Szval)<<"\t";
    //*out<<real(Pval)<<"\t";
    *out<<real(P2val)<<"\t";
    *out<<real(SRval)<<"\t"<<imag(SRval)<<"\t";
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


