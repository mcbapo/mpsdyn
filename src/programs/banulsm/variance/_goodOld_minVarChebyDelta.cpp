
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"
#include "HeisenbergHamiltonian.h"

using namespace shrt;

void sumMPO(MPO& result,const MPO& mpo1,const MPO& mpo2);

void computeValues(const MPS& state,const MPO& hamil,vector<complex_t>& valH2,vector<complex_t>& valH);
void computeLocalEnergy(const MPS& mps,double J_,double g_,double h_,vector<double>& locE);

void expandRDM(const MPS& mps,int posL,int posR,mwArray& A); // obsolete!
double computeRDMdistanceToId(const mwArray& rdm,int Lc);
void initialState(MPS& mps,int intSt,int D);
void initialStepState(MPS& mps,double thetaL,double phiL,double thetaR,double phiR,int D);
double getJacksonKPMcoeff(int k,int M);

/** Construct a file name to save (tmp) results */
const string mpsFileName(const string& mpsdir,const string& baseName,const string& label,int cnt);

/** Save the tmp files */
void tmpSave(MPS& mps,double norm_MPS,MPS& chebyT_Nm1,double norm_Nm1,MPS& chebyT_N,
	     double norm_N,const string& baseDir,const string& baseName,int cnt,vector<string>& tmpFiles);

/** Find tmp files and restore the status */
int tmpRestore(MPS& mps,double& norm_MPS,MPS& chebyT_Nm1,double& norm_Nm1,MPS& chebyT_N,
	       double& norm_N,const string& baseDir,const string& baseName,int M,vector<string>& tmpFiles);

/** SHOULD BE IN Contractor */
double getSchmidtVals(vector<double>& schmidt,const MPS& mps,int pos,char gaugeDir='R');
  
/** minVarChebyDelta constructs a Chebyshev poly approx to the delta
    function and uses it to filter out large energies, and thus
    decrease the variance.

    In this version, the expansion approximates directly the product
    of the Chebysev polynomial times the initial state, since the recursion is
    the same.

    I will use first order polynomials now

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
  bool ising=(props.getIntProperty("ising")!=0); // only if -ising=0 is given it will be Heisenberg (or XY)  
  // These are only used if Heisenberg model
  double Jx_ = props.getDoubleProperty("Jx");  
  double Jy_ = props.getDoubleProperty("Jy");  
  double Jz_ = props.getDoubleProperty("Jz");
  // It can also include x and y local fields Heisenberg model
  double hx_ = props.getDoubleProperty("hx");  
  double hy_ = props.getDoubleProperty("hy");  
  double delta=props.getDoubleProperty("delta");
  if(delta<0) delta=0.02; 
  int D = props.getIntProperty("D");  
  string outfname=props.getProperty("output");
  string Sfname=props.getProperty("outputS");
  string mpsdir=props.getProperty("mpsdir");
  int saveFr = props.getIntProperty("savingFreq");
  if(saveFr>0&&mpsdir.size()==0){
    cout<<"ERROR: Cannot save tmp results because no dir is given as argument!"<<endl;
    exit(1);
  }
  string mpsfileroot=props.getProperty("mpsfileroot");
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  int M=props.getIntProperty("M"); // LARGEST ORDER
  double Emin=props.getDoubleProperty("Emin");
  double Emax=props.getDoubleProperty("Emax");
  string efname=props.getProperty("energiesfile");
  int initSt=props.getIntProperty("initSt"); // initial state to use:
  //                               +1=X+; -1=X-; (2,3 for Y, Z)
  //                               +4,+5,+6 are the staggered(+-) for X, Y and Z and -4,-5,-6 the -+
  //                               +7 is a (complex) random state and +8 a real random (in which cases, the random seed may be given)
  //                               +10 is the staggered with period 4 (++--) in the Z direction, and -10 is (--++)
  //                               +11 is a step configuration (valid right now only for Ising)
  //                                   then requires Thleft and Thright to fix the excess of energy left-right
  //                                   and also a file to save the distribution of energy density
  double thetaL=props.getDoubleProperty("thetaL");
  double thetaR=props.getDoubleProperty("thetaR");
  double phiL=props.getDoubleProperty("phiL");
  double phiR=props.getDoubleProperty("phiR");
  string Edistfname=props.getProperty("outputEdist");
  int randSeed=props.getIntProperty("seed");

  int rate=props.getIntProperty("rate");
  bool checkRDM=(props.getIntProperty("checkRDM")>0); // whether to check how close RDM are to thermal
  if(rate<0) rate=1;
  int saveFreq=props.getIntProperty("saveFreq");
  if(saveFreq<0) saveFreq=M;
  saveFreq=lcm(rate,saveFreq);
  int lastSavedCnt=-1;
  int Lcmax=props.getIntProperty("Lcmax");
  if(Lcmax<0||Lcmax>10) Lcmax=6;
  int Lb=props.getIntProperty("Lb");
  if(Lb<0) Lb=6;
  bool avrg=1;
  avrg=!(props.getIntProperty("average")==0); // whether to average for each Lc subsystem

  bool isXY=(!ising&&Jz_==0);
  cout<<"Initialized arguments: L="<<L;
  if(ising){
    cout<<"; ISING: J="<<J_
	<<", g="<<g_
	<<", h="<<h_;
  }
  else{
    if(!isXY){
      cout<<"; HEISENBERG: Jx="<<Jx_
	  <<", Jy="<<Jy_
	  <<", Jz="<<Jz_
	  <<", h="<<h_;
    }
    else{ // XY case
      cout<<"; XY: Jx="<<Jx_
	  <<", Jy="<<Jy_
	  <<", hx="<<hx_
      	  <<", hy="<<hy_
	  <<", hz="<<h_;
    }    
  }
  cout<<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    if(ising){
      *out<<"% N\t J\t g\t h\t D\t cnt\t <H2>\t <H>";
    }
    else{
      if(!isXY)
	*out<<"% N\t Jx\t Jy\t Jz\t h\t D\t cnt\t <H2>\t <H>";
      else
	*out<<"% N\t Jx\t Jy\t hx\t hy\t hz\t D\t cnt\t <H2>\t <H>";
    }
    if(checkRDM){
      for(int kl=Lcmax;kl>=2;kl-=2){
	*out<<"\t"<<kl<<"\t dist(rho_"<<kl<<")";
      }
    }
    *out<<"\t S(L/2)";
    *out<<endl;
    out->close();
    delete out;
  }

  ofstream* outS;
  if(!file_exists(Sfname.data())){
    outS=new ofstream(Sfname.data());
    if(!outS->is_open()){
      cout<<"Error: impossible to open file "<<Sfname<<
	" for output"<<endl;
      exit(1);
    }
    if(ising){
      *outS<<"% N\t J\t g\t h\t D\t cnt\t lambda[0]^2\t lamdba[1]^2..."<<endl;
    }
    else{
      if(!isXY)
	*outS<<"% N\t Jx\t Jy\t Jz\t h\t D\t cnt\t lambda[0]^2\t lamdba[1]^2..."<<endl;
      else
	*outS<<"% N\t Jx\t Jy\t hx\t hy\t hz\t D\t cnt\t lambda[0]^2\t lamdba[1]^2..."<<endl;
    }
    *outS<<endl;
    outS->close();
    delete outS;
  }

  ofstream* outEd;
  if(ising&&initSt==11&&Edistfname.size()>0){
    if(!file_exists(Edistfname.data())){
      outEd=new ofstream(Edistfname.data());
      if(!outEd->is_open()){
	cout<<"Error: impossible to open file "<<Edistfname<<
	  " for output"<<endl;
	exit(1);
      }
      *outEd<<"% N\t J\t g\t h\t D\t cnt\t <h01> \t <h12>..."<<endl;
      outEd->close();
      delete outEd;
    }
  }

  
  cout<<setprecision(10);

  
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;

  MPS mps(L,D,d);
  if(ising&&initSt==11){
    initialStepState(mps,thetaL,phiL,thetaR,phiR,D);
  }
  else{
      initialState(mps,initSt,D);
  }
  cout<<"Initialized state ("<<initSt<<") from scratch, norm "<<contractor.contract(mps,mps)<<endl;
  mps.gaugeCond('R',1);
  mps.gaugeCond('L',1);
  
  
  // The original Hamiltonian (I'll also need the rescaled one for the Chebyshev expansion)
  // First: H
  //  IsingHamiltonian hamI(L,d,J,g,h); 
  Hamiltonian* ham0;
  if(ising){
    ham0=new IsingHamiltonian(L,d,J_,g_,h_);
  }
  else{
    ham0=new HeisenbergHamiltonian(L,Jx_,Jy_,Jz_,hx_,hy_,h_);
  }
  //IsingHamiltonian hamH0(L,d,J_,g_,h_);
  const MPO& hamil0=ham0->getHMPO();

  double scale,offset;
  if(Emin==-1&&Emax==-1){// not provided
    cout<<"First need to estimate the energy band to rescale "
	<<"and shift the Hamiltonian"<<endl;
    // First rescale the Hamiltonian:
    //hamH.registerTimeDependence();
    const MPO& hamil=ham0->getHMPO();
    // Get ground and highest state energy, for rescaling, but use always Dgs=60
    Emin=0;
    int Dgs=60;
    MPS gs(L,Dgs,d);
    gs.setRandomState();gs.gaugeCond('R');gs.gaugeCond('L',1);
    contractor.setEigenSolver(primme);
    contractor.findGroundState(hamil,Dgs,&Emin,gs);
      cout<<"Found GS of H with E="<<Emin<<endl;
      // Same for largest eigenvalue:
      Emax=0;
      MPS exc(L,Dgs,d);
      exc.setRandomState();exc.gaugeCond('R');exc.gaugeCond('L',1);
      {
	MPO hamilMinus(L);
	hamilMinus.setOp(0,new Operator(-1.*hamil.getOp(0).getFullData()),true);
	for(int k=1;k<L;k++){
	  hamilMinus.setOp(k,&hamil.getOp(k),false);
	}
	contractor.findGroundState(hamilMinus,Dgs,&Emax,exc);
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
  // Now the scale has to be such that E fits in [-1+delta,1+delta]
  scale=(1.-delta)/max(abs(Emax),abs(Emin));
  offset=0.;
  cout<<"Rescaling with scale="<<scale<<endl;
  
  // Now the rescaled H+alpha, such that the spec is in [-1+delta,1+delta]
  Hamiltonian* hamH;
  if(ising){
    hamH=new IsingHamiltonian(L,d,J_*scale,g_*scale,h_*scale,offset);
  }
  else{
    hamH=new HeisenbergHamiltonian(L,Jx_*scale,Jy_*scale,Jz_*scale,hx_*scale,hy_*scale,h_*scale);
  }
  //IsingHamiltonian hamH0(L,d,J_,g_,h_);
  const MPO& hamil_=hamH->getHMPO();

  //  IsingHamiltonian hamH(L,d,J_*scale,g_*scale,h_*scale,offset);
  //hamH.registerTimeDependence();
  cout<<"Created the rescaled Hamiltonian: scale="<<scale<<", offset="<<offset<<endl;
  // const MPO& hamil_=hamH.getHMPO();
  //  int Dh=3; // bond dimension of H, just in case

  cout<<"Initial state has <H>="<<contractor.contract(mps,hamil0,mps)
      <<", <H^2>="<<contractor.contract2(hamil0,mps)<<endl;
    

  // Now start the Chebyshev expansion computing T-1(0) and T0(Id) on the state
  MPS* chebyT_Nm1=new MPS(mps); // poly(-1)=0
  // Set the first tensor to 0=> ignore norm factor in this case
  chebyT_Nm1->gaugeCond('R',0);chebyT_Nm1->gaugeCond('L',1);
  chebyT_Nm1->setA(0,ZERO_c*chebyT_Nm1->getA(0).getA());
  double norm_Nm1=0.;
  
  MPS* chebyT_N=new MPS(mps); // normalized already
  double norm_N=0.;

  double coeffK=getJacksonKPMcoeff(0,M);
  double normMPS=log2(coeffK)-log2(M_PIl); // +1 for type II
  
  complex_t valsE,valsE2;
  double Shalf;

  // If intermediate steps are saved, I need mps with norm_MPS, chebyT_N with norm_N and chebyT_Nm1 with norm_Nm1
  vector<string> tmpFiles;
  int cnt=0;
  if(mpsdir.size()>0)
    cnt=tmpRestore(mps,normMPS,*chebyT_Nm1,norm_Nm1,*chebyT_N,norm_N,mpsdir,mpsfileroot,M,tmpFiles);
  if(cnt>0){
    cnt++; // do not want to repeat the iteration for cnt
    cout<<"After reading tmp files for step "<<cnt<<", normMPS="<<normMPS<<", norm_N="<<norm_N<<", norm_Nm1="<<norm_Nm1<<endl;
  }
      
  // iteration
  for(int k=cnt;k<=M;k++){
    MPS* chebyT_Np1;double norm_Np1(0.);
    if(k==0){
      chebyT_Np1=chebyT_N;  
      norm_Np1=norm_N;
    }
    else{
      chebyT_Np1=new MPS(*chebyT_N);
      if(k==1){
	// for Type I, T_1 is simply x*T0=> only multiply by H
	contractor.optimize(hamil_,*chebyT_N,*chebyT_Np1,D);
	chebyT_Np1->gaugeCond('R',0);chebyT_Np1->gaugeCond('L',0);
	norm_Np1=log2(chebyT_Np1->getNormFact())+norm_N;chebyT_Np1->setNormFact(1.);
      }
      else{
	// TODO: Optimize the sum of ops acting on vecs!
	MPS aux(*chebyT_N);
	contractor.optimize(hamil_,*chebyT_N,aux,D);
	aux.gaugeCond('R',0);aux.gaugeCond('L',0);
	double auxN=log2(aux.getNormFact());aux.setNormFact(1.);
	vector<const MPS*> kets;
	kets.push_back(&aux);kets.push_back(chebyT_Nm1);
	vector<complex_t> coefs;coefs.push_back(2.*ONE_c);
	coefs.push_back(-pow(2,norm_Nm1-norm_N-auxN)*ONE_c);
	contractor.optimizeSum(kets,coefs,*chebyT_Np1,D);
	chebyT_Np1->gaugeCond('R',0);chebyT_Np1->gaugeCond('L',0);
	norm_Np1=log2(chebyT_Np1->getNormFact())+norm_N+auxN;chebyT_Np1->setNormFact(1.);
      }
      // Shift the pointers and delete the discarded one 
      delete chebyT_Nm1;
      chebyT_Nm1=chebyT_N;norm_Nm1=norm_N;
      chebyT_N=chebyT_Np1;norm_N=norm_Np1;      
    }
    // Now compute the estimation of the state and the expectation values 
    if(k>0&&k%2==0){ // o.w. already done, and odd terms have 0 coeff (only for the delta!)
      MPS aux(mps);
      vector<const MPS*> kets;
      kets.push_back(&aux);kets.push_back(chebyT_Np1);
      vector<complex_t> coefs;coefs.push_back(ONE_c);
      coeffK=getJacksonKPMcoeff(k,M);
      if((k/2)%2!=0) coeffK=-coeffK;
      coefs.push_back((2./M_PIl)*coeffK*pow(2,norm_Np1-normMPS)*ONE_c);
      contractor.optimizeSum(kets,coefs,mps,D);
      mps.gaugeCond('R',0);mps.gaugeCond('L',0);
      normMPS+=log2(mps.getNormFact());mps.setNormFact(1.);
    

      // Now record values (notice that I do not need the overall norm for this)
      valsE=contractor.contract(mps,hamil0,mps);
      valsE2=contractor.contract2(hamil0,mps);
      Shalf=contractor.getEntropy(mps);
      // Save also the Schmidt values
      vector<double> schmidt;
      double entroSch=getSchmidtVals(schmidt,mps,L/2,'L');

      // debugging
      //    cout<<"Entropy computed to be "<<Shalf<<" and from Schmidt coeffs "<<entroSch<<endl;
      outS=new ofstream(Sfname.data(),ios::app);
      *outS<<setprecision(15);
      *outS<<L<<"\t";
      if(ising)
	*outS<<J_<<"\t"<<g_<<"\t"<<h_<<"\t";
      else
	if(!isXY)
	  *outS<<Jx_<<"\t"<<Jy_<<"\t"<<Jz_<<"\t"<<h_<<"\t";
	else
	  *outS<<Jx_<<"\t"<<Jy_<<"\t"<<hx_<<"\t"<<hy_<<"\t"<<h_<<"\t";
      *outS<<D<<"\t"<<k<<"\t";
      for(int vS=0;vS<D;vS++)
	*outS<<schmidt[vS]<<"\t";
      *outS<<endl;
      outS->close(); delete outS;
    
      // if necessary, also RDM
      vector<double> trDistRDM;
      vector<int> sizeRDM;
      if(checkRDM){
	if(!avrg){
	  int minPosR=L/2; // assuming L even and smallest Lc is 2
	  for(int k=L-1;k>=minPosR+1;k--){
	    mps.gaugeCond(k,'L',true);
	  }
	  int Lc=Lcmax;
	  int posL=(L-Lc)/2;int posR=posL+Lc-1; // outside tracing them
	  while(Lc>=2){
	    mwArray oper;
	    expandRDM(mps,posL,posR,oper);
	    //cout<<"Created RDM for Lc="<<Lc<<endl;
	    double singvals=0.;
	    double idFc=1./pow(2,Lc);
	    mwArray U;vector<complex_t> Dv;
	    wrapper::eig(oper,Dv,U); // only eigenvalues
	    //	  cout<<"Diagonalized the rho(Lc)"<<endl;
	    // now sum the absolute values
	    for(int k=0;k<Dv.size();k++){
	      double aux=real(Dv[k]);
	      singvals+=abs(aux-idFc);
	    }
	    if(Dv.size()<pow(2,Lc)){
	      cout<<"WARNING! Seems the dimension of the operator for D="<<D
		  <<", Lc="<<Lc<<"is not full: dim(oper)="<<Dv.size()
		  <<", 2^Lc="<<pow(2,Lc)<<endl;
	      singvals+=idFc*(pow(2,Lc)-Dv.size());
	    }
	    sizeRDM.push_back(Lc);
	    trDistRDM.push_back(singvals);
	    posL++;posR--;Lc=(posR-posL+1);
	  }
	}
	else{ // averaging over different subsystems of the same size
	  for(int lc=1;lc<=Lcmax;lc++){sizeRDM.push_back(lc);} // position lc-1
	  trDistRDM=vector<double>(sizeRDM.size(),0.);
	  vector<int> avrgCnt(sizeRDM.size(),0.);
	  int posL=Lb;
	  while(posL<L-Lb-1){ // at least there is a 2-site RDM to compute
	    int lc=Lcmax;
	    int posR=posL+lc-1;
	    while(posR>=L){ // cannot compute this one, need to reduce the size
	      lc--;
	      posR=posL+lc-1;
	    }
	    int dimL=pow(d,lc);
	    int dimR=dimL;
	    mwArray oper=contractor.getRDM(mps,posL,posR);
	    while(lc>=1){ // Compute the distance and trace another site
	      if(posR<L-Lb){
		// cout<<"Computing distance for rmd("<<lc<<") sites, pos:"<<posL<<"-"<<posR<<endl;
		double dist_lc=computeRDMdistanceToId(oper,lc);
		trDistRDM[lc-1]+=dist_lc;
		avrgCnt[lc-1]++; // to take the average later
	      }
	      if(lc>=1){
		// trace out two sites
		oper.reshape(Indices(dimL/d,d,dimR/d,d));
		oper.permute(Indices(1,3,2,4));
		dimL=dimL/d;dimR=dimR/d;
		oper.reshape(Indices(dimL*dimR,d*d));
		oper.multiplyRight(reshape(identityMatrix(d),Indices(d*d,1)));
		oper.reshape(Indices(dimL,dimR));
		posR--;
	      }
	      lc--;
	    }	    
	    posL++;
	  }
	  for(int lc=1;lc<=Lcmax;lc++) trDistRDM[lc-1]=trDistRDM[lc-1]/avrgCnt[lc-1];
	}
      }    
      cout<<"State for k="<<k<<"("<<M<<") found with <H^2>="<<real(valsE2)
	  <<" <H>="<<real(valsE)<<", normMPS="<<normMPS<<", norm(N+1)="<<norm_Np1<<endl;
      out=new ofstream(outfname.data(),ios::app);
      if(!out->is_open()){
	cout<<"Error: impossible to open file "<<outfname<<
	  " for output"<<endl;
	exit(1);
      }
      *out<<setprecision(15);
      // OJO!!!!! This is a mistake, doing it wrong for Heisenberg!!!!
      if(ising)
	*out<<L<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t";
      else
	if(!isXY)
	  *out<<L<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t";
      //  *out<<L<<"\t"<<Jx_<<"\t"<<Jy_<<"\t"<<Jz_<<"\t"<<h_<<"\t"<<D<<"\t"; // should be this
	else
	  *out<<L<<"\t"<<Jx_<<"\t"<<Jy_<<"\t"<<hx_<<"\t"<<hy_<<"\t"<<h_<<"\t"<<D<<"\t";
      *out<<k<<"\t"<<real(valsE2)<<"\t"<<real(valsE)<<"\t";
      if(checkRDM){
	for(int ks=0;ks<sizeRDM.size();ks++){
	  *out<<sizeRDM[ks]<<"\t"<<trDistRDM[ks]<<"\t";
	}
      }
      *out<<Shalf<<endl;
      *out<<endl;
      out->close(); delete out;

      // Case of inhomogeneous initial distrib: get and save also energy distribution
      
      if(ising&&initSt==11&&Edistfname.size()>0){
	vector<double> locE;
	computeLocalEnergy(mps,J_,g_,h_,locE);
	outEd=new ofstream(Edistfname.data(),ios::app);
	*outEd<<setprecision(10);
	*outEd<<L<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t"<<k<<"\t";
	for(int lpos=0;lpos<locE.size();lpos++)
	  *outEd<<locE[lpos]<<"\t";
	*outEd<<endl;
	outEd->close();delete outEd;
      }

      
    }
    // If the nr of the step is multiple of the desired saving frequency, save tmp files
    if(saveFr>0&&k%saveFr==0&&k>0){
      cout<<"About to save tmp files for step "<<k<<", normMPS="<<normMPS<<", norm_N="<<norm_N<<", norm_Nm1="<<norm_Nm1<<endl;
      tmpSave(mps,normMPS,*chebyT_Nm1,norm_Nm1,*chebyT_N,norm_N,mpsdir,mpsfileroot,k,tmpFiles);
    }    
  }

  // at the very end, delete the two left
  delete chebyT_Nm1,chebyT_N;
  // and if tmp files were kept, clean them up // NO!!! May be interesting to look at them
  // if(tmpFiles.size()>0){
  //   for(int k=0;k<tmpFiles.size();k++){
  //     string tmpF=tmpFiles[k];
  //     if(tmpF.length()>0&&file_exists(tmpF)){
  // 	cout<<"Removing tmp file "<<tmpF<<endl;
  // 	remove(tmpF.data());
  //     }
  //   }
  // }

  
}

void sumMPO(MPO& result,const MPO& mpo1,const MPO& mpo2){
  int L=mpo1.getLength();
  if(mpo2.getLength()!=L){
    cout<<"Error: incompatible MPO lengths for sum"<<endl;
    exit(1);
  }
  if(L==1){
    cout<<"Error: cannot do sum for MPO of length 1!!"<<endl;
    exit(1);
  }
  result.initLength(L);

  cout<<"Summing MPOs "<<mpo1<<" and "<<mpo2<<endl;
  for(int k=0;k<L;k++){
    // take both Operators involved
    const mwArray& op1=mpo1.getOpData(k);
    const mwArray& op2=mpo2.getOpData(k);
    Indices dims1=op1.getDimensions();
    Indices dims2=op2.getDimensions();
    if((dims1[0]!=dims2[0])||(dims1[2]!=dims2[2])){
      cout<<"Error: incompatible physical dimensions at position "<<k<<endl;
      exit(1);
    }
    int Dl=dims1[1]+dims2[1];
    if(k==0) Dl=1;
    int Dr=dims1[3]+dims2[3];
    if(k==L-1) Dr=1;
    cout<<"sumMPO("<<k<<") op1:"<<dims1<<", op2:"<<dims2<<"; new Dl="<<Dl<<", Dr="<<Dr<<endl;
    // Now fill in element by element
    mwArray blockOp(Indices(dims1[0],Dl,dims1[2],Dr));
    for(int l=0;l<Dl;l++){
      for(int r=0;r<Dr;r++){
	bool copy1=0;bool copy2=0;
	if((l<dims1[1])&&(r<dims1[3])) copy1=1;
	if((l>=dims1[1])&&(r>=dims1[3])) copy2=1;
	// special: edges
	if((k==0)&&(r<dims1[3])) copy1=1;
	if((k==0)&&(r>=dims1[3])) copy2=1;
	if((k==L-1)&&(l<dims1[1])) copy1=1;
	if((k==L-1)&&(l>=dims1[1])) copy2=1;
	// (operator element Dl, Dr): if Dl, Dr<dims1[2,3]
	//cout<<"For element "<<l<<","<<r<<" copy1="<<copy1
	//  <<" copy2="<<copy2<<endl;
	if(copy1+copy2){ // sth to copy
	  for(int i1=0;i1<dims1[0];i1++){
	    for(int i2=0;i2<dims1[2];i2++){
	      if(copy1) blockOp.setElement(op1.getElement(Indices(i1,l,i2,r)),Indices(i1,l,i2,r));
	      if(copy2){
		int indl=k==0?0:l-dims1[1];
		int indr=k==L-1?0:r-dims1[3];
		blockOp.setElement(op2.getElement(Indices(i1,indl,i2,indr)),Indices(i1,l,i2,r));
	      }
	    }
	  }
	}
      }
    }
    result.setOp(k,new Operator(blockOp),true);
  }

}

double getJacksonKPMcoeff(int k,int M){
  return ((M+1.-k)*cos(M_PIl*k/(M+1))+sin(M_PIl*k/(M+1))*cos(M_PIl/(M+1))/sin(M_PIl/(M+1)))/(M+1);
}

void computeLocalEnergy(const MPS& mps,double J_,double g_,double h_,vector<double>& locE){
  locE.clear();
  int L=mps.getLength();
  Contractor& contractor=Contractor::theContractor();
  IsingHamiltonian hamLoc(2,d,J_,.5*g_,.5*h_);
  const MPO& mpo=hamLoc.getHMPO();
  MPO auxMPO(L);
  Operator idOp(reshape(identityMatrix(d),Indices(d,1,d,1)));
  for(int k=0;k<L;k++)
      auxMPO.setOp(k,&idOp,false);
  complex_t normMPS=contractor.contract(mps,mps);
  // the rest identity
  for(int k=0;k<L-1;k++){
    auxMPO.setOp(k,&mpo.getOp(0),false);
    auxMPO.setOp(k+1,&mpo.getOp(1),false);
    complex_t valE=contractor.contract(mps,auxMPO,mps);
    auxMPO.setOp(k,&idOp,false);
    auxMPO.setOp(k+1,&idOp,false);
    locE.push_back(real(valE/normMPS));
  }

}

void __computeLocalEnergy(const MPS& mps,double J_,double g_,double h_,vector<double>& locE){
  locE.clear();
  int L=mps.getLength();
  Contractor& contractor=Contractor::theContractor();
  mwArray tmpL(ONE_c);tmpL.reshape(Indices(1,1));
  MPS state(mps);
  state.gaugeCond('L',1); // so that contraction from right is identity
  IsingHamiltonian hamLoc(2,d,J_,.5*g_,.5*h_);
  const MPO& mpo=hamLoc.getHMPO();
  MPO h12(4);
  h12.setOp(1,&mpo.getOp(0),false);
  h12.setOp(2,&mpo.getOp(1),false);
  for(int kl=0;kl<L-1;kl++){
    int Dl=mps.getA(kl).getDl();
    int Dr=mps.getA(kl+1).getDr();
    h12.setOp(0,new Operator(reshape(tmpL,Indices(Dl,1,Dl,1))),true);
    h12.setOp(3,new Operator(reshape(identityMatrix(Dr),Indices(Dr,1,Dr,1))),true);

    MPS aux(4,Dl,d);
    aux.replaceSite(0,reshape(identityMatrix(Dl),Indices(Dl,1,Dl)),false);
    aux.replaceSite(1,state.getA(kl).getA(),false);
    aux.replaceSite(2,state.getA(kl+1).getA(),false);
    aux.replaceSite(3,reshape(identityMatrix(Dr),Indices(Dr,Dr,1)),false);

    complex_t valE=contractor.contract(aux,h12,aux);
    cout<<"Computed value of E for link "<<kl<<" "<<valE<<endl;
    locE.push_back(real(valE));
    
    //contract one more site to tmpL
    if(kl<L-2){ // do by hand what Contractor.contractL does      
      int d=mps.getA(kl).getd();int Dl=mps.getA(kl).getDl();
      int Dr=mps.getA(kl).getDr();
      int betab=tmpL.getDimension(0);int betak=tmpL.getDimension(1);
      mwArray result=permute(mps.getA(kl).getA(),Indices(2,1,3));
      result.reshape(Indices(Dl,d*Dr));
      result=tmpL*result; // betab x d*Dr
      result.reshape(Indices(betab*d,Dr));
      mwArray aux=permute(mps.getA(kl).getA(),Indices(3,2,1)); // Dr x Dl x d (of bra)
      aux.conjugate();aux.reshape(Indices(Dr,Dl*d));
      tmpL=aux*result;
    }
  }
  
  
}

void computeValues(const MPS& state,const MPO& hamil,vector<complex_t>& valH2,vector<complex_t>& valH){
  valH2.clear();valH.clear();
  int L=state.getLength();
  Contractor& contractor=Contractor::theContractor();
  // Id op
  mwArray Id=identityMatrix(d);Id.reshape(Indices(d,1,d,1));
  Operator opId(Id); // Since all TI, only need edges and middle
  // Init MPO with all Ids  
  MPO aux(L);
  // Now, to start, I set the left edge on 0 and the right one on 1
  aux.setOp(0,&hamil.getOp(0),false);
  aux.setOp(1,&hamil.getOp(L-1),false);
  for(int k=2;k<L;k++) aux.setOp(k,&opId,false); // and the rest to Id
  for(int k=2;k<=L;k++){
    // compute expectation values
    complex_t valE=contractor.contract(state,aux,state);
    complex_t valE2=contractor.contract2(aux,state); // since they are hermitian, I contract H^+ H
    valH.push_back(valE);
    valH2.push_back(valE2);
    if(k<L){
      // one by one, move the right edge one to the right, and substitute in the former position by the middle operator
      aux.setOp(k,&hamil.getOp(L-1),false);
      aux.setOp(k-1,&hamil.getOp(1),false); // if (L==2) it does not enter here
    }    
  }
  
}


const string mpsFileName(const string& mpsdir,const string& baseName,const string& label,int cnt){
  stringstream s;
  s<<mpsdir<<"/"<<baseName<<"_"<<label;
  if(cnt>0) s<<"_"<<cnt;
  return s.str();
}

void tmpSave(MPS& mps,double norm_MPS,MPS& chebyT_Nm1,double norm_Nm1,MPS& chebyT_N,
	     double norm_N,const string& baseDir,const string& baseName,int cnt,vector<string>& tmpFiles){
  // first create new tmpFile names
  vector<string> existing(tmpFiles); // copy to remove at the end
  tmpFiles.clear();
  // now save ony by one
  string newTmp=mpsFileName(baseDir,baseName,"mps",cnt);
  double exNorm=mps.getNormFact();  
  mps.setNormFact(norm_MPS+log2(exNorm));
  mps.exportMPS(newTmp.data()); // saved containing log of norm!
  mps.setNormFact(exNorm); // restore
  tmpFiles.push_back(newTmp);
  
  newTmp=mpsFileName(baseDir,baseName,"Nm1",cnt);
  exNorm=chebyT_Nm1.getNormFact();  
  chebyT_Nm1.setNormFact(norm_Nm1+log2(exNorm));
  chebyT_Nm1.exportMPS(newTmp.data()); // saved containing log of norm!
  chebyT_Nm1.setNormFact(exNorm); // restore
  tmpFiles.push_back(newTmp);

  newTmp=mpsFileName(baseDir,baseName,"N",cnt);
  exNorm=chebyT_N.getNormFact();  
  chebyT_N.setNormFact(norm_N+log2(exNorm));
  chebyT_N.exportMPS(newTmp.data()); // saved containing log of norm!
  chebyT_N.setNormFact(exNorm); // restore
  tmpFiles.push_back(newTmp);

  cout<<"Saved tmp files for step "<<cnt<<" to ";
  for(int p=0;p<tmpFiles.size();p++)
    cout<<tmpFiles[p]<<", ";
  cout<<endl;
  // now remove the old one
  if(existing.size()>0){
    for(int k=0;k<existing.size();k++){
      string tmpF=existing[k];
      if(tmpF.length()>0&&file_exists(tmpF)){
	remove(tmpF.data());
      }
    }
  }
}


int tmpRestore(MPS& mps,double& norm_MPS,MPS& chebyT_Nm1,double& norm_Nm1,MPS& chebyT_N,
	     double& norm_N,const string& mpsdir,const string& baseName,int M,vector<string>& tmpFiles){
  tmpFiles.clear();
  // first find if there are tmp files
  int cnt=M;
  bool found=0;
  while(cnt>0&&!found){
    string mpsfile=mpsFileName(mpsdir,baseName,"mps",cnt);
    if(file_exists(mpsfile)){
      cout<<"findTmpFile found mpsfile="<<mpsfile<<endl;
      // Check also for the others
      string Nm1file=mpsFileName(mpsdir,baseName,"Nm1",cnt);
      string Nfile=mpsFileName(mpsdir,baseName,"N",cnt);
      if(file_exists(Nm1file)&&file_exists(Nfile)){
	// Read the files and return
	mps.importMPS(mpsfile.data());
	chebyT_Nm1.importMPS(Nm1file.data());
	chebyT_N.importMPS(Nfile.data());
	norm_MPS=mps.getNormFact();mps.setNormFact(1.);
	norm_Nm1=chebyT_Nm1.getNormFact();chebyT_Nm1.setNormFact(1.);
	norm_N=chebyT_N.getNormFact();chebyT_N.setNormFact(1.);
	tmpFiles.push_back(mpsfile);
	tmpFiles.push_back(Nm1file);
	tmpFiles.push_back(Nfile);
	found=true;
	return cnt;
      }
      else{
	cout<<"ERROR: Found a file for MPS but not the corresponding ones for Nm1 and N=> no tmp file used"<<endl;
	return 0;
      }
    }
    else{
      cnt--;
    }
  }
  if(!found) return 0;
  cout<<"ERROR: Should not be here"<<endl;
  exit(1);
}


void expandRDM(const MPS& mps,int posL,int posR,mwArray& A){
  A=mps.getA(posL).getA();
  int d=A.getDimension(0);int Dl0=A.getDimension(1);int Dr=A.getDimension(2);
  A.permute(Indices(2,1,3));A.reshape(Indices(Dl0*d,Dr));
  int dphys=d; // temporary dimensions of the vector
  int Dv=Dr;
  for(int k=posL+1;k<=posR;k++){
    mwArray aux=mps.getA(k).getA();
    //cout<<"Multiplying tensor for site "<<k<<", A="<<aux.getDimensions()<<" with tmp result "<<A.getDimensions()<<endl;
    Indices dims=aux.getDimensions();
    int d=dims[0];int Dl=dims[1];int Dr=dims[2];
    if(Dl!=Dv){
      cout<<"Error: Dimensions do not agree in expandVec!"<<endl;
      exit(1);
    }
    aux.permute(Indices(2,1,3));
    aux.reshape(Indices(Dl,d*Dr));
    A.multiplyRight(aux);
    dphys=dphys*d;
    Dv=Dr;
    A.reshape(Indices(Dl0*dphys,Dv));
  }
  A.reshape(Indices(Dl0,dphys,Dv));A.permute(Indices(2,1,3));A.reshape(Indices(dphys,Dl0*Dv));
  // if the mps has a normalization factor, we need to include it!!
  A=mps.getNormFact()*A;
  A.multiplyRight(Hconjugate(A));
}  


double computeRDMdistanceToId(const mwArray& rdm,int Lc){
  double singvals=0.;
  double idFc=1./pow(2,Lc);
  mwArray U;vector<complex_t> Dv;
  wrapper::eig(rdm,Dv,U); // only eigenvalues
  // now sum the absolute values
  for(int k=0;k<Dv.size();k++){
    double aux=real(Dv[k]);
    singvals+=abs(aux-idFc);
  }
  if(Dv.size()<pow(2,Lc)){
    cout<<"WARNING! Seems the dimension of the operator for Lc="<<Lc
	<<"is not full: dim(oper)="<<Dv.size()
	<<", 2^Lc="<<pow(2,Lc)<<endl;
    singvals+=idFc*(pow(2,Lc)-Dv.size());
  }
  return singvals;
}


void initialState(MPS& mps,int initSt,int D){
  int L=mps.getLength();
  if(initSt==+1||abs(initSt)==4) mps.setProductState(p_xplus);
  else{
    if(initSt==-1) mps.setProductState(p_xminus);
    else{
      if(initSt==+2||abs(initSt)==5) mps.setProductState(p_yplus);
      else{
	if(initSt==-2) mps.setProductState(p_yminus);
	else{
	  if(initSt==+3||abs(initSt)==6) mps.setProductState(p_zero);
	  else{
	    if(initSt==-3) mps.setProductState(p_one);
	    else{ // any other number, I start random
	      cout<<"Initial state a non-translational random MPS with D="<<mps.getBond()<<endl;
	      mps.setRandomState(); // This creates it with complex values
	      // If real, discard imaginary parts
	      if(initSt==8){
		for(int pos=0;pos<L;pos++){
		  mwArray aux=mps.getA(pos).getA();
		  mps.replaceSite(pos,.5*(aux+conjugate(aux)));
		}
	      }
	      mps.gaugeCond('R',1);
	    }
	  }
	}
      }
    }
  }
  // Now for the staggered, change one every two by applying a single site op (sigma_z to x,y sigma_y to z)  
  if(abs(initSt)>=4&&abs(initSt)<=6){
    int firstPos=initSt<0?0:1; //+- changes first site 1
    complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
    mwArray oper;
    if(abs(initSt)<6) oper=mwArray(Indices(d,d),dataz);
    else oper=mwArray(Indices(d,d),datax);
    for(int k=firstPos;k<L;k+=2){
	mps.applyLocalOperator(k,oper);
    }
  }
  // And for the staggered with period 2 I apply the same but every second pair
  if(abs(initSt)==10){
    mps.setProductState(p_zero);
    complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    mwArray oper=mwArray(Indices(d,d),datax);
    for(int k=0;k<L;k++){
      if(initSt>0){
	if(k%4>1) mps.applyLocalOperator(k,oper);
      }
      else{
	if(k%4<=1) mps.applyLocalOperator(k,oper);
      }
    }
  }
  mps.increaseBondDimension(D);
}

void initialStepState(MPS& mps,double thetaL,double phiL,double thetaR,double phiR,int D){
  int L=mps.getLength();
  mwArray auxL(Indices(d,1));
  mwArray auxR(Indices(d,1));
  auxL.setElement(cos(thetaL)*ONE_c,Indices(0,0));
  auxL.setElement(sin(thetaL)*exp(phiL*I_c),Indices(1,0));
  auxR.setElement(cos(thetaR)*ONE_c,Indices(0,0));
  auxR.setElement(sin(thetaR)*exp(phiR*I_c),Indices(1,0));
  auxL.reshape(Indices(d,1,1));
  auxR.reshape(Indices(d,1,1));
  for(int pos=0;pos<L;pos++){
    if(pos<L/2)
      mps.replaceSite(pos,auxL,false);
    else
      mps.replaceSite(pos,auxR,false);
  }
  mps.gaugeCond('R',1);
  mps.increaseBondDimension(D);  
}

double getSchmidtVals(vector<double>& schmidt,const MPS& mps,int pos,char gaugeDir){
  Contractor& contractor=Contractor::theContractor();
  int len=mps.getLength();
  mwArray rdm=contractor.contract(mps,mps,pos,gaugeDir);
  // cout<<"getSchmidtVals of mps:"<<mps<<" at cut pos="<<pos<<"\nComputed the rdm: "<<rdm.getDimensions()<<endl;
  // Now diagonalize
  rdm.permute(Indices(2,1)); // here it does not matter
  vector<complex_t> Dval;mwArray U; //placeholder:ignored
  int D=rdm.getDimension(0); //mps.getA(pos).getDl();
  wrapper::eig(rdm,Dval,U);
  schmidt.clear();
  for(int k=0;k<D;k++) schmidt.push_back(real(Dval[k]));
  double entropy=0.;
  for(int k=0;k<D;k++){
    double tmp=real(Dval[k]);
    if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
  }
  //cout<<"getSchmidtVals computed vals="<<Dval<<"\n and Entropy="<<entropy<<endl;
  return entropy;
  
}

