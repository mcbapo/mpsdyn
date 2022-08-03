
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

double getJacksonKPMcoeff(int k,int M);

/** renyiCheby tests the approximation to the ensemble that optimizes
    the 2-Renyi entropy with Chebyshev polynomials for the step
    function.

    In renyChebiStep, the iteration constructs the polynomials and the
    theta function operator, but not the true rho2, which would be
    Theta(H)*(H+alpha), properly normalized.

    Notice that it is not necessary to construct the sum of
    polynomials for Theta(H), since the expectation values are linear
    in that, and can be computed by adding the contributions c_m g_m
    trace[ T_m(H) (H+alpha) Op]. In this version it is computed explicitly (rhoS2), 
    so that we can compare to exact and also compute the 2-Renyi entropy,
    but could be dropped for more efficiency!!
 */

int d=2;

int main(int argc,const char* argv[]){
  // Read input arguments
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
  double delta=props.getDoubleProperty("delta");
  if(delta<0) delta=0.02; 
  double alpha=props.getDoubleProperty("alpha"); // parameter of the ensemble
  bool positive=(props.getIntProperty("positiveBranch")>0); // which branch to take: positive means theta(E-alpha)
  int D=props.getIntProperty("D");
  // int dpur=props.getIntProperty("dpur");
  // if(dpur<0) dpur=d;
  string outfname=props.getProperty("outputfile");
  string outfnameFull=props.getProperty("outputfileFull");
  string mpsfname=props.getProperty("mpsfile");
  bool app=(props.getIntProperty("append")!=0);
  bool centerAlpha=(props.getIntProperty("center")==1); // then H+alpha is used as argument for the expansion, and theta finishes at 0
  int M=props.getIntProperty("M"); // LARGEST ORDER
  double Emin=props.getDoubleProperty("Emin");
  double Emax=props.getDoubleProperty("Emax");
  // TODO! Save energies to a file, so that I can reuse them!
  double scale,offset;
  // string mpsfname=props.getProperty("mpsfile");
  // string initmpsfname=props.getProperty("initmpsfile");

  cout<<"Initialized arguments: L="<<L
      <<", J="<<J_
      <<", g="<<g_
      <<", h="<<h_
      <<", alpha="<<alpha
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;


  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;

  // The original Hamiltonian to get the energy
  IsingHamiltonian hamH0(L,d,J_,g_,h_);
  const MPO& hamil0=hamH0.getHMPO();

  if(Emin==-1&&Emax==-1){
    cout<<"First need to estimate the energy band to rescale "
	<<"and shift the Hamiltonian"<<endl;
    // First rescale the Hamiltonian:
    //hamH.registerTimeDependence();
    cout<<"Created the non-rescaled Hamiltonian"<<endl;
    const MPO& hamil=hamH0.getHMPO();
    if(Emin==-1&&Emax==-1){// not provided
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
    }
    else{
      cout<<"Using for extreme energies the provided arguments Emin="<<Emin
	  <<" Emax="<<Emax<<endl;
    }
  }
  // Some cases are trivial and require no projection! Namely, if all
  // eigenvalues fullfill E+alpha>0 or E+alpha<0, the ensemble is just
  // H+alpha, with the proper normalization (which can include a
  // global sign). This is exactly a small MPO, but not computed here
  if(Emin>-alpha||Emax<-alpha){
    cout<<"ERROR: With this value of alpha, the ensemble is trivial!"<<endl;
    exit(1);    
  }
  double x0;  
  // // Now the scale has to be such that E+alpha fits in [-1+delta,1+delta]
  if(centerAlpha){
    scale=(1.-delta)/max(abs(Emax+alpha),abs(Emin+alpha));
    offset=alpha*scale;
    x0=0.; // the theta function is always wrt 0
  }
  else{
    // Alternative, we scale the H without offset, and apply a theta function with varying cut
    scale=(1.-delta)/max(abs(Emax),abs(Emin));
    offset=0.;
    x0=-scale*alpha; // where to stop the theta function (positive, then E needs to be >x0)
  }
  
  // Now the rescaled H for the expansion, such that the spec is in [-1+delta,1+delta]
  IsingHamiltonian hamH(L,d,J_*scale,g_*scale,h_*scale,offset);
  //hamH.registerTimeDependence();
  cout<<"Created the rescaled Hamiltonian: scale="<<scale<<", offset="<<offset<<endl;
  const MPO& hamil_=hamH.getHMPO();
  int Dh=3; // bond dimension of H, just in case

  
  // Prepare the observables to compute

  MPS idMPS(L,1,d*d);
  for(int k=0;k<L;k++)
    idMPS.setA(k,reshape(identityMatrix(d),Indices(d*d,1,1)));

  MPS hamilMPS(L,Dh,d*d);
  MPSfromMPO(hamil_,hamilMPS);

  MPS hamil0MPS(L,Dh,d*d);
  MPSfromMPO(hamil0,hamil0MPS);

  MPS SxCenter(idMPS);
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  sigX.reshape(Indices(d*d,1,1));
  SxCenter.setA(L/2,sigX);//sigX.reshape(Indices(d,d));

  MPS SyCenter(idMPS);
  complex_t dataY[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  sigY.reshape(Indices(d*d,1,1));
  SyCenter.setA(L/2,sigY);//sigY.reshape(Indices(d,d));

  MPS SzCenter(idMPS);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  sigZ.reshape(Indices(d*d,1,1));
  SzCenter.setA(L/2,sigZ);//sigZ.reshape(Indices(d,d));
  
  MPO hamilId(L); // needed for the recursion
  extendMPO(hamil_,hamilId,d);
  //  hamilId.exportForMatlab("hId.m");

  //
  MPO hamilAlpha(L);
  if(centerAlpha){
    hamilAlpha=hamilId;
  }
  else{
    IsingHamiltonian hamH(L,d,J_*scale,g_*scale,h_*scale,-x0);
    //hamH.registerTimeDependence();
    const MPO& hamilA_=hamH.getHMPO();
    extendMPO(hamilA_,hamilAlpha,d);
  }
  
  // I first compute the coefficients C_m of the Chebyshev expansion of the step function
  // TODO: Simply compute them as needed: no need to keep this vector
  vector<double> Cs(M+1,0.);
  if(centerAlpha){
    Cs[0]=.5*M_PIl;
    Cs[0]*=getJacksonKPMcoeff(0,M);
    for(int k=1;k<=M;k++){
      if(k%2==0){
	Cs[k]=0.;
      }
      else{
	double signS=(((k-1)/2)%2==0)?+1.:-1.;
	Cs[k]=positive?signS*2./k:-signS*2./k;
      }
      Cs[k]*=getJacksonKPMcoeff(k,M);
    }
  }
  else{
    if(positive)
      Cs[0]=acos(x0);
    else
      Cs[0]=M_PIl-acos(x0);
    Cs[0]*=getJacksonKPMcoeff(0,M);
    for(int k=1;k<=M;k++){
      // the factor 2 is for the normalization factor, leaving out a 1/Pi for every term
      if(positive)
	Cs[k]=2.*sin(k*acos(x0))/k;
      else
	Cs[k]=-2.*sin(k*acos(x0))/k;
      Cs[k]*=getJacksonKPMcoeff(k,M);
    }
  }

  // For output
  ofstream* out;
  if(!app||!file_exists(outfname)){
    out=new ofstream(outfname);
    *out<<"% Chebyshev expansion of the step, M="<<M<<endl;
    *out<<"% k\t J\t g\t h\t D\t tr(rho_k)\t Energy\t tr(rho_k^2)\t <X(L/2)>\t<Y(L/2)>\t<Z(L/2)>"<<endl;
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    out->close();delete out;
  }
  
  // Now start the Chebyshev expansion computing T0(=Id) and T1(=X)
  MPS* chebyT_Nm1=new MPS(idMPS); 
  MPS* chebyT_N=new MPS(hamilMPS); // for type I, no factor 2!
  chebyT_Nm1->gaugeCond('R',0);chebyT_Nm1->gaugeCond('L',0);
  chebyT_N->gaugeCond('R',0);chebyT_N->gaugeCond('L',0);
  double norm_Nm1=log2(chebyT_Nm1->getNormFact());chebyT_Nm1->setNormFact(1.);
  double norm_N=log2(chebyT_N->getNormFact());chebyT_N->setNormFact(1.);
  
  MPS rhoS2(idMPS); // times Cs[0], but I put it in the norm
  //  rhoS2.setA(0,Cs[0]*rhoS2.getA(0).getA()); // include the first coeff
  rhoS2.gaugeCond('R',0);rhoS2.gaugeCond('L',0);
  double normRho2=log2(Cs[0]*rhoS2.getNormFact());rhoS2.setNormFact(1.);
  
  // place for trace, energy, trace^2 and ev of X Y Z at the center
  double trRho(0.);
  double trRho2(0.);
  double trRhoH(0.);
  double valX(0.);
  double valY(0.);
  double valZ(0.);
      
  // iteration
  for(int k=0;k<=M;k++){
    MPS* chebyT_Np1;double norm_Np1(0.);
    if(k==0){
      chebyT_Np1=chebyT_Nm1; 
      norm_Np1=norm_Nm1;
    }
    else if(k==1){
      chebyT_Np1=chebyT_N;
      norm_Np1=norm_N;
    }
    else{
      chebyT_Np1=new MPS(*chebyT_N);
      // TODO: Optimize the sum of ops acting on vecs!
      MPS aux(*chebyT_N);
      contractor.optimize(hamilId,*chebyT_N,aux,D);
      aux.gaugeCond('R',0);aux.gaugeCond('L',0);
      double auxN=log2(aux.getNormFact());aux.setNormFact(1.);
      vector<const MPS*> kets;
      kets.push_back(&aux);kets.push_back(chebyT_Nm1);
      vector<complex_t> coefs;coefs.push_back(2.*ONE_c);
      coefs.push_back(-pow(2,norm_Nm1-norm_N-auxN)*ONE_c);
      contractor.optimizeSum(kets,coefs,*chebyT_Np1,D);
      chebyT_Np1->gaugeCond('R',0);chebyT_Np1->gaugeCond('L',0);
      norm_Np1=log2(chebyT_Np1->getNormFact())+norm_N+auxN;chebyT_Np1->setNormFact(1.);
      // Shift the pointers and delete the discarded one 
      delete chebyT_Nm1;
      chebyT_Nm1=chebyT_N;norm_Nm1=norm_N;
      chebyT_N=chebyT_Np1;norm_N=norm_Np1;      
    }
    // Now compute the estimation of rho and the expectation values 
    if(k>0){ // o.w. already done
      MPS aux(rhoS2);
      vector<const MPS*> kets;
      kets.push_back(&aux);kets.push_back(chebyT_Np1);
      vector<complex_t> coefs;coefs.push_back(ONE_c);
      coefs.push_back(Cs[k]*pow(2,norm_Np1-normRho2)*ONE_c);
      contractor.optimizeSum(kets,coefs,rhoS2,D);
      rhoS2.gaugeCond('R',0);rhoS2.gaugeCond('L',0);
      normRho2+=log2(rhoS2.getNormFact());rhoS2.setNormFact(1.);
    }
    double normF=pow(2,norm_Np1);
    trRho+=Cs[k]*normF*real(contractor.contract(*chebyT_Np1,hamilAlpha,idMPS));
    // tr(rho^2) cannot be computed accumulating unless I get overlaps between polys
    double traceFull=real(contractor.contract(rhoS2,hamilAlpha,idMPS));
    trRho2=real(contractor.contract2(hamilAlpha,rhoS2))/(traceFull*traceFull);
    // would all go with 1/scale in front
    trRhoH+=Cs[k]*normF*real(contractor.contract(*chebyT_Np1,hamilAlpha,hamil0MPS));
    valX+=Cs[k]*normF*real(contractor.contract(*chebyT_Np1,hamilAlpha,SxCenter));
    valY+=Cs[k]*normF*real(contractor.contract(*chebyT_Np1,hamilAlpha,SyCenter));
    valZ+=Cs[k]*normF*real(contractor.contract(*chebyT_Np1,hamilAlpha,SzCenter));
    
    // and save values to file
    out=new ofstream(outfname,ios::app);
    *out<<k<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t"
	<<trRho<<"\t"<<trRhoH<<"\t"<<trRho2<<"\t"
	<<valX<<"\t"<<valY<<"\t"<<valZ
	<<endl;
    out->close();
    delete out;
    // also values from Rho directly
    complex_t enerFull=contractor.contract(rhoS2,hamilAlpha,hamil0MPS);
    complex_t valXFull=contractor.contract(rhoS2,hamilAlpha,SxCenter);
    complex_t valYFull=contractor.contract(rhoS2,hamilAlpha,SyCenter);
    complex_t valZFull=contractor.contract(rhoS2,hamilAlpha,SzCenter);
    out=new ofstream(outfnameFull,ios::app);
    *out<<k<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t"
	<<traceFull<<"\t"<<real(enerFull)<<"\t"<<trRho2<<"\t"
	<<real(valXFull)<<"\t"<<real(valYFull)<<"\t"<<real(valZFull)
	<<"\t"<<normRho2<<"\t"<<norm_Np1
	<<endl;
    out->close();
    delete out;

    cout<<"After step k="<<k<<" tr(rho^2)="<<trRho2<<", tr(rhoH)="<<enerFull/traceFull
	<<", norm_Np1="<<norm_Np1<<", normR2="<<normRho2<<endl;
  }

  // at the very end, delete the two left
  delete chebyT_Nm1,chebyT_N;

  // And export the resulting rho, for comparison
  cout<<"Coming out of the iteration! with trR2="<<trRho2<<
    ", energy="<<trRhoH/trRho<<endl;
  if(mpsfname.length()>0){
    // Watch out: this is not rho2, but Theta(H) alone. Needs to be multiplied by H+alpha
    rhoS2.exportMPS(mpsfname.data());
  }
  // if(L<=12){
  //   MPO toExport(L);
  //       vector<int> dimS(L,d);
  //   MPOfromMPS(rhoS2,toExport,dimS);
  //   stringstream str;
  //   str<<"mpo2RenCh_M"<<M<<"_D"<<D<<".m";
  //   toExport.exportForMatlab(str.str().data()); //"mpo2Ren.m");
  // }
  if(L<=10){
    mwArray rhoExact;
    cout<<"Computing exact 2-Renyi ensemble for alpha="<<alpha<<endl;
    mwArray Uex;
    vector<complex_t> Dex;
    mwArray Hex;
    expandOper(hamil0,Hex);
    wrapper::eig(Hex,Dex,Uex,true);

    mwArray diagRhoExact(Hex.getDimensions());
    for(int k=0;k<Hex.getDimension(0);k++){
      if((positive&&(real(Dex[k])+alpha>=0))||(!positive&&(real(Dex[k])+alpha<0))){
	diagRhoExact.setElement(Dex[k]+alpha*ONE_c,Indices(k,k));
      }
    }
    complex_t trRex=diagRhoExact.trace();
    //    cout<<"Trace of rhoEx="<<trRex<<endl;
    diagRhoExact.multiplyLeft((1./real(trRex))*ONE_c);
    rhoExact=Uex*diagRhoExact*Hconjugate(Uex);

    cout<<"******** Comparing to exact solution ******* "<<endl;
    mwArray fullRho;
    MPO rhoMPO(L);
    MPOfromMPS(rhoS2,rhoMPO);
    expandOper(rhoMPO,fullRho);
    fullRho=fullRho*(Hex+alpha*identityMatrix(pow(d,L)));
    complex_t trRmps=fullRho.trace();
    fullRho=(ONE_c/trRmps)*fullRho;
    cout<<"Now trace of expanded rho: "<<fullRho.trace()<<endl;
    cout<<"Expanded MPS "<<rhoS2<<" to "<<fullRho.getDimensions()<<endl;
    complex_t trR2ex,trR2mps;
    mwArray aux(fullRho);
    trR2mps=(aux*Hconjugate(aux)).trace();
    trR2ex=(rhoExact*rhoExact).trace();
    complex_t Emps=(aux*Hex).trace();
    complex_t Eex=(rhoExact*Hex).trace();
    cout<<"Exact tr(rho^2)="<<trR2ex<<", tr(rho H)="<<Eex<<endl;
    cout<<"Computed MPS tr(rho^2)="<<trR2mps<<", tr(rho H)="<<Emps<<endl;
    
    // Distance to exact?
    mwArray U,S,Vd;
    wrapper::svd(aux-rhoExact,U,S,Vd);
    cout<<"(Trace) Distance between both: "<<S.trace()<<endl;
  }

  
}



double getJacksonKPMcoeff(int k,int M){
  return ((M+1.-k)*cos(M_PIl*k/(M+1))+sin(M_PIl*k/(M+1))*cos(M_PIl/(M+1))/sin(M_PIl/(M+1)))/(M+1);
}

