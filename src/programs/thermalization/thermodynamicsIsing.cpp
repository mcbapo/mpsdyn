
#include <math.h>
#include <iomanip>

#include "Properties.h"
#include "misc.h"
#include "MPS.h"
#include "FoldedOperator.h"
#include "DoubleOperator.h"
#include "Contractor.h"
#include "SpinMPO.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

/** Construct the MPO for the thermal state of a finite chain under
    the Ising Hamiltonian using the purification ansatz, i.e., an MPO
    for rho(beta/2) and rho(beta)=rho(beta/2)^dagger rho(beta/2).
    Compute total magnetization and energy as function of beta. Also
    ground state (and highest energy state) energy and magnetizations.
    
    Arguments:
    \param <N> (int) length of the chain
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian (transverse field)
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian (parallel field)
    \param <D> (int) bond dimension for the half MPO
    \param <beta> (double) inverse temperature
    \param <delta> (double) width of time step: this is the increment in beta
                   as imaginary time evolution is simulated (meaning that each
		   exponential has coefficient delta/2)
                   The time step turns out to be 2*delta, because the answer 
		   is rho^+ rho
    \param <stepBlk> (int) save results every stepBlk steps (means
            every delta*stepBlk*2 in beta)
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)
    \param <mpsfile> (char*) optional argument, with a file name where to store the 
                     final MPS.


*/

void computeGS(MPS& state,double* ener,const MPO& hamil,int Dgs,const string& mpsdir,const string modelstr,bool maxE=false);
const string mpsfilename(int N,int D,const string& mpsdir,const string& modelstr,bool maxE=false);

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

  int N = props.getIntProperty("L");  
  double J = props.getDoubleProperty("J");  
  double g = props.getDoubleProperty("g");  
  double h = props.getDoubleProperty("h");  
  int Dgs = props.getIntProperty("Dgs");  
  int D = props.getIntProperty("D");  
  double beta = props.getDoubleProperty("beta");  
  double delta = props.getDoubleProperty("delta");  
  int stepBlk = props.getIntProperty("stepBlock");
  if(stepBlk<=0) stepBlk=1;
  string outfnameGS=props.getProperty("outputFileGS");
  string outfname=props.getProperty("outputFile");
  bool app=(props.getIntProperty("append")!=0);
  string mpsDir=props.getProperty("mpsdir"); // for GS and max exctd states
  //  bool saveFinal=mpsFile.length()>0;


  cout<<"Initialized arguments: L="<<N
      <<", J="<<J
      <<", g="<<g
      <<", h="<<h
      <<", Tmin="<<1/beta
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  // String to identify states of this model
  stringstream str;
  str<<"J"<<J<<"_g"<<g<<"_h"<<g;
  
  ofstream* out;
  //  if(!app){
  out=new ofstream(outfname);
  *out<<"% N\t J\t g\t h\t D\t delta\t nrsteps\t beta\t T\t Energy/N\t <Sx>/N\t <Sz>/N\t <sigx(N/2)\t <sigz(N/2)>\t log(Z)"
      <<endl;
  //    *out<<setprecision(12);
  // }
  // else
  //   out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);

  int M=beta/(2*delta); // nr of steps I need in each beta/2 (positive and negative!)
  // Maybe I need an extra step, if beta/2 is not divisible by delta
  double remainder=beta/2-M*delta;
  // I will need a step of remainder on each half of the state


  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d0=2;
  MPS thS(N,1,d0*d0);
  thS.setProductState(p_maxent); // initial identity, with funny normalization
  cout<<"Initialized identity state. Check norm (unnormalized should be 2^N="<<pow(2,N)<<")="
      <<contractor.contract(thS,thS)<<endl;
  // I keep it for later (identity contraction)
  MPS idState(thS);

  // Now I normalize and keep the norm
  thS.gaugeCond('R',1);
  thS.gaugeCond('L',1);

  double logZ=.5*N*log(2); // place for the log of the norm (log of sqrt(Z))
  double logZneg=logZ; // same for negative beta
  
  MPS thSneg(thS); // for negative beta

  IsingHamiltonian hamI(N,d0,J,g,h); //
  const MPO& mpoH=hamI.getHMPO();
  cout<<"Initialized basic Hamiltonian"<<endl;

  
  // MPOs fot total magnetization
  MPO doubleSx(N);
  MPO doubleSz(N);

  // I now also want the individual spin operators on a single site!!
  mwArray opX,opZ;
  mwArray sig0=identityMatrix(d);//identity
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),datax);//sigmax
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz
  constructOperatorProduct(opX,sigX,sig0);
  constructOperatorProduct(opZ,sigZ,sig0);
  
  {
    MPS gs(N,Dgs,d0);
    MPS gsNeg(N,Dgs,d0); // max. excited state
    
    double Emin=0.;
    double Emax=0.;
    
    computeGS(gs,&Emin,mpoH,Dgs,mpsDir,str.str(),0);
    cout<<"Ground state found with eigenvalue "<<Emin
	<<" (Energy per particle="<<Emin/N<<")"<<endl;
    computeGS(gsNeg,&Emax,mpoH,Dgs,mpsDir,str.str(),1);
    cout<<"Highest energy state found with eigenvalue "<<Emax
	<<" (Energy per particle="<<Emax/N<<")"<<endl;
    MPO Sx(N);
    SpinMPO::getSxMPO(N,d0,Sx);
    extendMPO(Sx,doubleSx,d0);
    MPO Sz(N);
    SpinMPO::getSzMPO(N,d0,Sz);
    extendMPO(Sz,doubleSz,d0);

    complex_t Sxgs=contractor.contract(gs,Sx,gs);
    complex_t Szgs=contractor.contract(gs,Sz,gs);
    complex_t SxgsNeg=contractor.contract(gsNeg,Sx,gsNeg);
    complex_t SzgsNeg=contractor.contract(gsNeg,Sz,gsNeg);


    complex_t SxHalf,SzHalf,SxHalfNeg,SzHalfNeg;
    {
      MPS aux(gs);
      aux.applyLocalOperator(N/2,sigX);
      SxHalf=contractor.contract(aux,gs);
      aux.applyLocalOperator(N/2,sigX);aux.applyLocalOperator(N/2,sigZ);
      SzHalf=contractor.contract(aux,gs);
      aux=gsNeg;
      aux.applyLocalOperator(N/2,sigX);
      SxHalfNeg=contractor.contract(aux,gsNeg);
      aux.applyLocalOperator(N/2,sigX);aux.applyLocalOperator(N/2,sigZ);
      SzHalfNeg=contractor.contract(aux,gsNeg);
    }
    
    ofstream* outGS;
    if(!app){
      outGS=new ofstream(outfnameGS);
      *outGS<<"% N\t J\t g\t h\t D\t level(0/1)\t Energy\t <Sx>\t <Sz>\t <sigx(N/2)\t <sigz(N/2)>"
	  <<endl;
      *outGS<<setprecision(12);
    }
    else
      outGS=new ofstream(outfnameGS,ios::app);
    if(!outGS->is_open()){
      cout<<"Error: impossible to open file "<<outfnameGS<<
	" for output"<<endl;
      exit(1);
    }
    *outGS<<setprecision(12);
    *outGS<<N<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<Dgs<<"\t"<<0<<"\t"
	  <<Emin<<"\t"<<real(Sxgs)<<"\t"<<real(Szgs)
	  <<"\t"<<real(SxHalf)<<"\t"<<real(SzHalf)
	  <<endl;
    *outGS<<N<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<Dgs<<"\t"<<1<<"\t"
	  <<Emax<<"\t"<<real(SxgsNeg)<<"\t"<<real(SzgsNeg)
	  <<"\t"<<real(SxHalfNeg)<<"\t"<<real(SzHalfNeg)
	  <<endl;

    outGS->close();
    delete outGS;
  }

  // double MPO for the H to be applied on rho
  MPO doubleH(N);
  extendMPO(mpoH,doubleH,d0);

  MPO evolU(N); // MPO for imaginary time evolution
  MPO evolUneg(N); // MPO for imaginary time evolution with -H
  {
    MPO auxU(N);
    hamI.getUMPO(auxU,-delta*.5*ONE_c,0.);
    doubleMPO(auxU,evolU); // the folded e^(-delta/2 H)xe^(-delta/2 H^T)
    hamI.getUMPO(auxU,delta*.5*ONE_c,0.);
    doubleMPO(auxU,evolUneg); // the folded e^(delta/2 H)xe^(delta/2 H^T)
  }

  cout<<"Before applying any evolution step, energy is E(Id)="
      <<contractor.contract(thS,doubleH,thS)/contractor.contract(thS,thS)
      <<endl;

  // Now apply the evolution
  int leftSteps=M; // nr of steps left to apply
  double betai=0; // current value of beta
  int cnt=0;
  while(leftSteps>0){
    int nrToApply=min(stepBlk,leftSteps);
    for(int k=0;k<nrToApply;k++){
      MPS aux(thS);
      contractor.optimize(evolU,aux,thS,D);
      MPS aux2(thSneg);
      contractor.optimize(evolUneg,aux2,thSneg,D); 
    }
    thS.gaugeCond('R',false);
    thSneg.gaugeCond('R',false); 
    logZ+=log(thS.getNormFact());
    logZneg+=log(thSneg.getNormFact());
    thS.setNormFact(1.);
    thSneg.setNormFact(1.);
    leftSteps-=nrToApply;
    betai+=2*delta*nrToApply;
    cnt+=nrToApply;
    // compute/record observables now
    double norm2=real(contractor.contract(thS,thS)); //trace of rho^2
    complex_t energy2=contractor.contract(thS,doubleH,thS);
    complex_t magX=contractor.contract(thS,doubleSx,thS);
    complex_t magZ=contractor.contract(thS,doubleSz,thS);
    // Also compute the sigx and sigz on a single site in the center
    complex_t SxHalf,SzHalf;
    {
      MPS aux(thS);
      aux.applyLocalOperator(N/2,opX);
      SxHalf=contractor.contract(aux,thS);
      aux.applyLocalOperator(N/2,opX);aux.applyLocalOperator(N/2,opZ);
      SzHalf=contractor.contract(aux,thS);
    }
    
    
    *out<<N<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"
	<<delta<<"\t"<<2*cnt<<"\t"<<betai<<"\t"<<1./betai<<"\t"
	<<real(energy2)/norm2<<"\t"
	<<real(magX)/norm2<<"\t"
	<<real(magZ)/norm2<<"\t"
	<<real(SxHalf)/norm2<<"\t"
	<<real(SzHalf)/norm2<<"\t"
	<<2*logZ<<"\t"
	<<endl;
    cout<<"At the end of step "<<cnt<<", beta="<<betai<<", T="<<1./betai
	<<", energy="<<real(energy2)/norm2
	<<endl;

    norm2=real(contractor.contract(thSneg,thSneg)); //trace of rho^2
    energy2=contractor.contract(thSneg,doubleH,thSneg);
    magX=contractor.contract(thSneg,doubleSx,thSneg);
    magZ=contractor.contract(thSneg,doubleSz,thSneg);
    {
      MPS aux(thSneg);
      aux.applyLocalOperator(N/2,opX);
      SxHalf=contractor.contract(aux,thSneg);
      aux.applyLocalOperator(N/2,opX);aux.applyLocalOperator(N/2,opZ);
      SzHalf=contractor.contract(aux,thSneg);
    }

    *out<<N<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"
	<<-delta<<"\t"<<2*cnt<<"\t"<<-betai<<"\t"<<-1./betai<<"\t"
	<<real(energy2)/norm2<<"\t"
	<<real(magX)/norm2<<"\t"
	<<real(magZ)/norm2<<"\t"
	<<real(SxHalf)/norm2<<"\t"
	<<real(SzHalf)/norm2<<"\t"
	<<2*logZneg<<"\t"
	<<endl;
  }
  // If I am not exactly at the desired beta, I have to apply an extra step
  if(remainder!=0){
    MPO auxU(N);
    hamI.getUMPO(auxU,remainder/2,0.,true);
    doubleMPO(auxU,evolU); // the folded e^(-delta/2 H)xe^(-delta/2 H^T)
    hamI.getUMPO(auxU,-remainder/2,0.,true);
    doubleMPO(auxU,evolUneg); // the folded e^(-delta/2 H)xe^(-delta/2 H^T)
    MPS aux(thS);
    contractor.optimize(evolU,aux,thS,D); 
    MPS aux2(thSneg);
    contractor.optimize(evolUneg,aux2,thSneg,D); 
    betai+=remainder*2;cnt++;
    thS.gaugeCond('R',false);
    thSneg.gaugeCond('R',false);
    logZ+=log(thS.getNormFact());
    logZneg+=log(thSneg.getNormFact());
    thS.setNormFact(1.);
    thSneg.setNormFact(1.);
    // compute/record observables now
    double norm2=real(contractor.contract(thS,thS)); //trace of rho^2
    complex_t energy2=contractor.contract(thS,doubleH,thS);
    complex_t magX=contractor.contract(thS,doubleSx,thS);
    complex_t magZ=contractor.contract(thS,doubleSz,thS);
    *out<<N<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"
	<<delta<<"\t"<<2*cnt<<"\t"<<betai<<"\t"<<1./betai<<"\t"
	<<real(energy2)/norm2<<"\t"
	<<real(magX)/norm2<<"\t"
	<<real(magZ)/norm2<<"\t"
	<<2*logZ<<"\t"
	<<endl;
    cout<<"At the end of step "<<cnt<<", beta="<<betai<<", T="<<1./betai
	<<", energy="<<real(energy2)/norm2
	<<endl;
    norm2=real(contractor.contract(thSneg,thSneg)); //trace of rho^2
    energy2=contractor.contract(thSneg,doubleH,thSneg);
    magX=contractor.contract(thSneg,doubleSx,thSneg);
    magZ=contractor.contract(thSneg,doubleSz,thSneg);
    *out<<N<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"
	<<-delta<<"\t"<<2*cnt<<"\t"<<-betai<<"\t"<<-1./betai<<"\t"
	<<real(energy2)/norm2<<"\t"
	<<real(magX)/norm2<<"\t"
	<<real(magZ)/norm2<<"\t"
	<<2*logZneg<<"\t"
	<<endl;
  }
  // Finally, if needed, save the state
  // if(saveFinal)
  //   thS.exportMPS(mpsFile.data());

}

const string mpsfilename(int N,int D,const string& mpsdir,const string& modelstr,bool maxE){
  stringstream str;
  int level=maxE?1:0;
  str<<mpsdir<<"/MPS_L"<<N<<"_"<<modelstr<<"_l"<<level<<"_D"<<D;
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return str.str();
}

void computeGS(MPS& state,double* ener,const MPO& hamil,int Dgs,const string& mpsdir,const string modelstr,bool maxE){
  int N=hamil.getLength();
  Contractor& contractor=Contractor::theContractor();
  MPO hamilNeg(N);
  if(maxE){
    // need to change sign of H
    hamilNeg.setOp(0,new Operator(-1.*hamil.getOp(0).getFullData()),true);
    for(int k=1;k<N;k++)
      hamilNeg.setOp(k,&hamil.getOp(k),0);
  }
  bool found=0;
  // Now try to check for the MPS file
  const string mpsfile=mpsfilename(N,Dgs,mpsdir,modelstr,maxE);
  if(file_exists(mpsfile)){
    state.importMPS(mpsfile.data());
    cout<<"Recovered existing file for D="<<Dgs<<"=> no optimization"<<endl;
    state.gaugeCond('R',1);
    found=1;
    // compute the energy and return  
    if(maxE){
      *ener=real(contractor.contract(state,hamilNeg,state));
    }
    else{
      *ener=real(contractor.contract(state,hamil,state));
    }
    return;
  }
  else{
    // try any smaller D
    int D_=Dgs-1;
    while(!found&&D_>0){ 
      const string mpsfileD=mpsfilename(N,D_,mpsdir,modelstr,maxE);
      if(file_exists(mpsfileD)){
	state.importMPS(mpsfileD.data());
	cout<<"Imported file for D="<<D_<<endl;
	found=1;
      }
      else D_--;
    }
    
    if(!found){
      cout<<"No file found for "<<(!maxE)<<"; Start with random state."<<endl;
      state.setRandomState(); // the intial state, random
    }
    state.gaugeCond('R',1);
    state.gaugeCond('L',1);
    // Now try to find the ground state
    if(maxE){
      contractor.findGroundState(hamilNeg,Dgs,ener,state);
      *ener=-*ener;
    }
    else{
      contractor.findGroundState(hamil,Dgs,ener,state);      
    }
    // Now save it
    state.exportMPS(mpsfile.data());
  }
  
  
}
