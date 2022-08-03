
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

bool constructInitialRho(int L,int d,int initSt,MPS& rho0);

/** Prepare a name for the mps correspondng to step cnt */
const string mpsFileName(const string& mpsdir,const string& baseName,int cnt);

/** Save a tmp mps file, replacing the previous one */
void tmpSave(MPS& rho,const string& baseDir,const string& baseName,
	     int cnt,string& tmpFile);

/** Check for tmp files */
int findTmpFile(int M,const string& mpsdir,const string& baseName,string& mpsfile);


/** Try to approximate the diagonal ensemble by projecting the initial
    state onto the null space of the commutator. For this, we
    implement an iterative approximation to the pseudoinverse, and use
    that the projector we want is actually
    /f[P=Id-A A^P\f]
    where \f$A^P\f$ is the Moore-Penrose pseudoinverse.

    Optimizations could be to keep positivity or use a more efficient
    approximation for the projector. 

    NOTE: Seems to be working in the exact case (precise comparison to
    MAtlab still missing) but this can only be achieved for L=2 or 3
    sites! (since the effective physical dimension for the projector,
    which is living in the same space as the commutator, is d^4=16, so
    for L=4 we need already D=256 for the exact result!).
    Q: Optimizations??? A weaker version acting only on the state??
*/


int d=2;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  srandom(time(NULL));
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
  double delta=props.getDoubleProperty("delta");
  if(delta<0) delta=0.02; 
  int D = props.getIntProperty("D");  
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir");
  string mpsfile=props.getProperty("mpsfile");
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  int M=props.getIntProperty("M"); // LARGEST STEP NUMBER
  double Emin=props.getDoubleProperty("Emin");
  double Emax=props.getDoubleProperty("Emax");
  string efname=props.getProperty("energiesfile");
  int initSt=props.getIntProperty("initSt"); // initial state to use:
  //                               +1=X+; -1=X-; (2,3 for Y, Z)
  //                               +4,+5,+6 are the staggered(+-) for X, Y and Z and -4,-5,-6 the -+
  //                               +7 is a (complex) random state and +8 a real random (in which cases, the random seed may be given)
  int randSeed=props.getIntProperty("seed");

  int rate=props.getIntProperty("rate");
  if(rate<0) rate=1;
  int saveFreq=props.getIntProperty("saveFreq");
  if(saveFreq<0) saveFreq=M;
  saveFreq=lcm(rate,saveFreq);
  int lastSavedCnt=-1;

  cout<<"Initialized arguments: L="<<L
      <<", J="<<J_
      <<", g="<<g_
      <<", h="<<h_
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  if(!file_exists(outfname.data())){
    cout<<"File *"<<outfname.data()<<"* not found: not appending"<<endl;
    app=0; // also if it does not exist (new)
  }
  ofstream* out;
  if(!app){
    cout<<"REcreating output file *"<<outfname.data()<<"*"<<endl;
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% cnt\t D\t real(tr)\t im(tr)\t re(E)\t im(E)\t re(Sx)\t re(Sy)\t re(Sz)\tcomm\tlog2(normP)";
    *out<<endl;
    out->close();
    delete out;
  }
  cout<<setprecision(10);

  
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  
  // The original Hamiltonian 
  IsingHamiltonian hamH(L,d,J_,g_,h_);
  const MPO& hamil=hamH.getHMPO();
    // The original Hamiltonian to get the energy
  MPS hmps(L,1,d*d);
  MPSfromMPO(hamil,hmps);


  double scale;
  if(Emin==-1&&Emax==-1){// not provided
    cout<<"First need to estimate the energy band to rescale the term Hc Hc^+"<<endl;
    // First rescale the Hamiltonian:
    //hamH.registerTimeDependence();
    const MPO& hamil=hamH.getHMPO();
    // Get ground and highest state energy, for rescaling
    Emin=0;
    MPS gs(L,D,d);
    gs.setRandomState();gs.gaugeCond('R');gs.gaugeCond('L',1);
      contractor.findGroundState(hamil,D,&Emin,gs);
      cout<<"Found GS of H with E="<<Emin<<endl;
      // Same for largest eigenvalue:
      Emax=0;
      MPS exc(L,D,d);
      exc.setRandomState();exc.gaugeCond('R');exc.gaugeCond('L',1);
      {
	MPO hamilMinus(L);
	hamilMinus.setOp(0,new Operator(-1.*hamil.getOp(0).getFullData()),true);
	for(int k=1;k<L;k++){
	  hamilMinus.setOp(k,&hamil.getOp(k),false);
	}
	contractor.findGroundState(hamilMinus,D,&Emax,exc);
	Emax=-Emax;
      }
      cout<<"Found max Exc of H with E="<<Emax<<endl;
      if(efname.size()!=0){ // save for other cases
	ofstream* outE=new ofstream(efname.data());
	*outE<<setprecision(15);
	*outE<<Emin<<endl;
	*outE<<Emax<<endl;
	outE->close();delete outE;
      }
    }
  else{
    cout<<"Using for extreme energies the provided arguments Emin="<<Emin
	<<" Emax="<<Emax<<endl;
  }
  scale=(1.-delta)/((Emax-Emin)*(Emax-Emin));
  cout<<"Scale factor thus: "<<scale<<endl;
  
  int Dh=3; // bond dimension of H, just in case

  //Commutator superoperator
  MPO commMPO(L);
  hamH.getCommutatorMPO(commMPO,0.); //commutator constructed
  cout<<"Commutator constructed"<<endl;

    // Prepare the observables to compute
  // The identity, to compute the trace
  MPS idMPS(L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++)
    idMPS.setA(k,id2);
  // And I will also need the sigmaz/sigmay/sigmax operator on central site, which I could compute ony by one, but I assemble it here
  MPS SxCenter(idMPS);
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  sigX.reshape(Indices(d*d,1,1));
  SxCenter.setA(L/2,sigX);sigX.reshape(Indices(d,d));

  MPS SyCenter(idMPS);
  complex_t dataY[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  sigY.reshape(Indices(d*d,1,1));
  SyCenter.setA(L/2,sigY);sigY.reshape(Indices(d,d));

  MPS SzCenter(idMPS);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  sigZ.reshape(Indices(d*d,1,1));
  SzCenter.setA(L/2,sigZ);sigZ.reshape(Indices(d,d));


    // Initial state rho0
  MPS rho0(L,D,d*d);
  bool isReal=constructInitialRho(L,d,initSt,rho0);
  rho0.gaugeCond('R',0);rho0.gaugeCond('L',0);
  double normMPS=log2(rho0.getNormFact()); // hopefully 0 (pure states)
  rho0.setNormFact(1.);
  complex_t E0=contractor.contract(rho0,hmps);
  complex_t tr0=contractor.contract(rho0,idMPS);
  complex_t valx0=contractor.contract(rho0,SxCenter);
  complex_t valy0=contractor.contract(rho0,SyCenter);
  complex_t valz0=contractor.contract(rho0,SzCenter);
  complex_t valC2=contractor.contract2(commMPO,rho0);
  cout<<"Initial state has trace="<<tr0<<", <H>="<<E0/tr0
      <<", |[H,rho0]|^2="<<pow(2,2*normMPS)*valC2<<endl;

  out=new ofstream(outfname.data(),ios::app);
  *out<<0<<"\t"<<D<<"\t"<<real(tr0)<<"\t"<<imag(tr0)<<"\t"
      <<real(E0)<<"\t"<<imag(E0)<<"\t"
      <<real(valx0)<<"\t"<<real(valy0)<<"\t"<<real(valz0)<<"\t"
      <<real(valC2)<<"\t"<<0.
      <<endl;
  out->close();delete out;

  // Now, to start the iteration, the first approx to the projector is Id-scale*Hc*Hc
  MPS idMPSdoub(L,1,d*d*d*d); // vectorized superoperator Identity
  mwArray id22=identityMatrix(d*d);id22.reshape(Indices(d*d*d*d,1,1));
  for(int k=0;k<L;k++)
    idMPSdoub.setA(k,id22);
  MPS projNull(idMPSdoub);
  double normP=0;
  {
    MPS aux(L,D,d*d);
    MPSfromMPO(commMPO,aux);
    MPO comm2(L);
    extendMPO(commMPO,comm2,d*d);
    //    cout<<"comm2="<<comm2<<endl;
    MPS aux2(aux);
    contractor.optimize(comm2,aux2,aux,D);
    cout<<"Approximated square of commutator"<<endl; //: "<<aux<<endl;
    //    cout<<" to be summed to idMPS: "<<idMPSdoub<<endl;
    
    vector<const MPS*> vecs;vector<complex_t> coeffs;
    vecs.push_back(&idMPSdoub);
    vecs.push_back(&aux);
    coeffs.push_back(ONE_c); 
    coeffs.push_back(-scale*ONE_c); // coeff of the Hc Hc
    contractor.optimizeSum(vecs,coeffs,projNull,D); // Do I want this to be exact? Dh*Dh+1
    cout<<"Approximated Id-alpha* square of commutator"<<endl;
    projNull.gaugeCond('R',0);projNull.gaugeCond('L',0);
    normP+=log2(projNull.getNormFact());
    projNull.setNormFact(1.);
  }

  // Now iterate:
  // make MPO out of proj
  // apply to initial state (compute exp values)
  // apply to proj MPS for next order

  int cnt=1;
  // TODO: read temporary results and restart from intermediate step
  while(cnt<=M){
    MPO projMPO(L);
    MPOfromMPS(projNull,projMPO);
    complex_t trM=contractor.contract(rho0,projMPO,idMPS);
    complex_t enerM=contractor.contract(rho0,projMPO,hmps);
    complex_t valxM=contractor.contract(rho0,projMPO,SxCenter);
    complex_t valyM=contractor.contract(rho0,projMPO,SyCenter);
    complex_t valzM=contractor.contract(rho0,projMPO,SzCenter);
    // would like the commutator squared, too, but it may be too expensive: remove?
    complex_t valC2;
    {
      MPO Pcomm2P(L);
      {
	const MPO* mpos[4]={&projMPO,&commMPO,&commMPO,&projMPO};
	MPO::join(4,mpos,Pcomm2P);
      }
      valC2=contractor.contract(rho0,Pcomm2P,rho0);
    }

    // write results
    out=new ofstream(outfname.data(),ios::app);
    *out<<cnt<<"\t"<<D<<"\t"<<real(trM)<<"\t"<<imag(trM)<<"\t"
	<<real(enerM)<<"\t"<<imag(enerM)<<"\t"
	<<real(valxM)<<"\t"<<real(valyM)<<"\t"<<real(valzM)<<"\t"
	<<real(valC2)<<"\t"<<normP<<endl;
    out->close();delete out;
    cout<<"step "<<cnt<<" <Sx>="<<real(valxM)/real(trM)
	<<" <Sy>="<<real(valyM)/real(trM)
	<<" <Sz>="<<real(valzM)/real(trM)
	<<" |[Hc,rho]|^2="<<real(valC2)*pow(2,2*normP)
	<<endl;
    // one step further
    if(cnt<M){
      MPO extProj(L);
      extendMPO(projMPO,extProj,d*d);
      MPS aux(projNull);
      contractor.optimize(extProj,aux,projNull,D);
      projNull.gaugeCond('R',0);projNull.gaugeCond('L',0);
      normP+=log2(projNull.getNormFact());
      projNull.setNormFact(1.);
    }
    cnt++;
  }
}



bool constructInitialRho(int L,int d,int initSt,MPS& rho0){
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sig0=identityMatrix(d);
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray sigX(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  bool isReal=false;  

  // Product states
  if(initSt<9){
    mwArray state;
    switch(initSt){
    case 1: // X+
      {
	state=.5*(sig0+sigX);
	isReal=true;
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
	isReal=true;
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
  if(initSt==9){
    // Special case: |Y+>+|Y->
    MPS auxP(L,1,d);auxP.setProductState(p_yplus);
    MPS auxM(L,1,d);auxM.setProductState(p_yminus);
    Contractor& contr=Contractor::theContractor();
    MPS tmp(L,2,d);tmp.setRandomState();
    vector<const MPS*>vecs;vector<complex_t> coeffs(2,ONE_c);
    vecs.push_back(&auxP);vecs.push_back(&auxM);
    contr.optimizeSum(vecs,coeffs,tmp,2);
    tmp.gaugeCond('L',1);
    tmp.gaugeCond('R',1);
    // And now transform to double MPS
    for(int k=0;k<L;k++){
      mwArray tens=tmp.getA(k).getA(); //d x Dl xDr
      int Dl=tens.getDimension(1);int Dr=tens.getDimension(2);
      tens.reshape(Indices(d*Dl*Dr,1));
      tens.multiplyRight(Hconjugate(tens));
      tens.reshape(Indices(d,Dl,Dr,d,Dl,Dr));
      tens.permute(Indices(1,4,2,5,3,6));
      tens.reshape(Indices(d*d,Dl*Dl,Dr*Dr));
      rho0.replaceSite(k,tens,0);
    }
    isReal=true;
  }
  
  return isReal;
}
