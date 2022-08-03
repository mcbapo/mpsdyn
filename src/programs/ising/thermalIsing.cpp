
#include <math.h>
#include <iomanip>

#include "MPS.h"
//#include "FoldedOperator.h"
#include "DoubleOperator.h"
#include "Contractor.h"
#include "Properties.h"
#include "misc.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

/** Construct the MPO for the thermal state of a finite chain under
    the Ising Hamiltonian using the purification ansatz, i.e., an MPO
    for rho(beta/2) and rho(beta)=rho(beta/2)^dagger rho(beta/2).
    
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

void computeZZcorrelations(MPS& rho,ofstream* out);
void computeLocalObservables(MPS& rho,ofstream* out);
// Compute all strings of subL Pauli matrices from site n0
void computeNrdm(MPS& rho,int subL,int n0,ofstream* out);

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
  int D = props.getIntProperty("D");  
  double beta = props.getDoubleProperty("beta");  
  double delta = props.getDoubleProperty("delta");  
  int stepBlk = props.getIntProperty("stepBlock");  
  string outfname=props.getProperty("outputFile");
  bool app=(props.getIntProperty("append")!=0);
  string mpsFile=props.getProperty("mpsfile");
  bool saveFinal=mpsFile.length()>0;

    
  cout<<"Initialized arguments: L="<<N
      <<", g="<<g
      <<", h="<<h
      <<", Tmin="<<1/beta
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app||!file_exists(outfname.data())){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% Thermal Ising with J="<<J<<", g="<<g<<" h="<<h<<", L="<<N<<endl;

    *out<<"% N\t g\t h\t D\t delta\t nrsteps\t beta\t <H>\t <h(L/2)> "
	<<"<sigx(L/2)>\t <sigy(L/2)>\t <sigz(L/2)>\t "<<endl; 
    
    //    *out<<setprecision(12);
    out->close();delete out;
  }

  int M=abs(beta/(2*delta)); // nr of steps I need in each beta/2
  if(beta<0) delta=-abs(delta); // they have to have the same sign
  // Maybe I need an extra step, if beta/2 is not divisible by delta
  double remainder=beta/2-M*delta;
  // I will need a step of remainder on each half of the state

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d0=2;
  MPS thS(N,1,d0*d0);
  thS.setProductState(p_maxent); // initial identity, with funny normalization
  // cout<<"Initialized identity state. Check norm (unnormalized should be 2^N="<<pow(2,N)<<")="
  //     <<contractor.contract(thS,thS)<<endl;
  // I keep it for later (identity contraction)
  //  MPS idState(thS);
 
  IsingHamiltonian hamI(N,d0,J,g,h); 
  const MPO& mpoH=hamI.getHMPO();
  cout<<"Initialized basic Hamiltonian"<<endl;

  // double MPO for the H to be applied on rho
  MPO doubleH(N);
  extendMPO(mpoH,doubleH,d0);

  // To compute the local energy density in the center, I use a trick:
  // define a Hamiltonian on only two sites, with half the g and h
  // In gnral, it would be better to compute the rdm in the center and
  // compute all 3-site (for instance) products of Paulis
  MPO localH(N); Operator* idOp;
  {
    IsingHamiltonian ham12(2,d0,J,.5*g,.5*h);
    const MPO& h12=ham12.getHMPO();
    // fill in the long MPO with identities
    mwArray id=identityMatrix(d0*d0);id.reshape(Indices(d0*d0,1,d0*d0,1));
    idOp=new Operator(id);
    MPO auxH12(2);extendMPO(h12,auxH12,d0);
    for(int k=0;k<N;k++){ 
      localH.setOp(k,idOp,false);
      if(k==N/2-1)
	localH.setOp(N/2-1,new Operator(auxH12.getOp(0).getFullData()),true);
      if(k==N/2)
	localH.setOp(N/2,new Operator(auxH12.getOp(1).getFullData()),true);
    }
  }


  MPO evolU(N); // MPO for imaginary time evolution
  {
    MPO auxU(N);
    hamI.getUMPO(auxU,-delta*.5*ONE_c,0.);
    doubleMPO(auxU,evolU); // the folded e^(-delta/2 H)xe^(-delta/2 H^T)
  }
  
  cout<<"Before applying any evolution step, energy is E(Id)="
      <<contractor.contract(thS,doubleH,thS)/contractor.contract(thS,thS)
      <<" and local energy is "<<contractor.contract(thS,localH,thS)
      <<endl;

  // Now apply the evolution
  int leftSteps=M; // nr of steps left to apply
  double betai=0; // current value of beta
  int cnt=0;
  while(leftSteps>0){
    int nrToApply=min(stepBlk,leftSteps);
    //cout<<"stepBlk="<<stepBlk<<" leftSteps="<<leftSteps<<", nrToApply="<<nrToApply<<", M="<<M<<" delta="<<delta<<endl;
    for(int k=0;k<nrToApply;k++){
      MPS aux(thS);
      contractor.optimize(evolU,aux,thS,D); 
    }
    thS.gaugeCond('R',1);
    leftSteps-=nrToApply;
    betai+=2*delta*nrToApply;
    cnt+=nrToApply;
    // compute/record observables now
    double norm2=real(contractor.contract(thS,thS)); //trace of rho^2
    complex_t energy2=contractor.contract(thS,doubleH,thS);
    complex_t locEnergy=contractor.contract(thS,localH,thS);
    //MPS aux(thS);
    //contractor.optimize(evolU,aux,thS,D); 
    cout<<"After "<<nrToApply<<" steps, E="<<energy2<<endl;
    cout<<"At the end of step "<<cnt<<", beta="<<betai<<", T="<<1./betai
	<<", energy="<<real(energy2)/norm2<<", local energy="<<real(locEnergy)/norm2
	<<endl;
    out=new ofstream(outfname.data(),ios::app);
    *out<<N<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"
	<<delta<<"\t"<<cnt<<"\t"<<betai<<"\t"
	<<real(energy2)<<"\t"
	<<real(locEnergy)<<"\t";
    computeLocalObservables(thS,out);
    *out<<endl;
    out->close();delete out;
  }
  // If I am not exactly at the desired beta, I have to apply an extra step
  if(remainder!=0){
    MPO auxU(N);
    hamI.getUMPO(auxU,remainder/2,0.,true);
    doubleMPO(auxU,evolU); // the folded e^(-delta/2 H)xe^(-delta/2 H^T)
    MPS aux(thS);
    contractor.optimize(evolU,aux,thS,D); 
    betai+=remainder*2;cnt++;
    // compute/record observables now
    double norm2=real(contractor.contract(thS,thS)); //trace of rho^2
    complex_t energy2=contractor.contract(thS,doubleH,thS);
    complex_t locEnergy=contractor.contract(thS,localH,thS);
    cout<<"At the end of step "<<cnt<<", beta="<<betai<<", T="<<1./betai
	<<", energy="<<real(energy2)/norm2<<", local energy="<<real(locEnergy)/norm2
	<<endl;
    out=new ofstream(outfname.data(),ios::app);
    *out<<N<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"
	<<delta<<"\t"<<cnt<<"\t"<<betai<<"\t"
	<<real(energy2)<<"\t"
	<<real(locEnergy)<<"\t";
    computeLocalObservables(thS,out);
    *out<<endl;
    out->close();delete out;
  }
  // Finally, if needed, save the state
  if(saveFinal)
    thS.exportMPS(mpsFile.data());

  
}

void computeZZcorrelations(MPS& rho,ofstream* out){
  // Z operators
  static bool init(false);
  static mwArray sigZ_Id;
  if(!init){
    init=true;
    int d0=2;
    mwArray sig0=identityMatrix(d0);sig0.reshape(Indices(1,d0*d0));//identity
    complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
    sigZ_Id=mwArray(Indices(d0*d0,1),dataz);//sigmaz
    sigZ_Id.multiplyRight(sig0);
    sigZ_Id.reshape(Indices(d0,d0,d0,d0));
    sigZ_Id.permute(Indices(1,3,2,4));
    sigZ_Id.reshape(Indices(d0*d0,d0*d0));
  }
  int L=rho.getLength();
  // correls of L/2,L/2+i for i=1 to the end
  MPS aux(rho);
  aux.applyLocalOperator(L/2,sigZ_Id);
  Contractor& contractor=Contractor::theContractor();
  complex_t valZ=contractor.contract(aux,rho);
  *out<<real(valZ)<<"\t";
  // now, for distances i=1 to the largest one, compute <Z(L/2)Z(L/2+i)> and <Z(L/2+i)>
  for(int k=L/2+1;k<L;k++){
    // apply Z on k=L/2+i
    aux.applyLocalOperator(k,sigZ_Id);
    complex_t valZZ=contractor.contract(aux,rho);
    // undo the L/2
    aux.applyLocalOperator(L/2,sigZ_Id);
    valZ=contractor.contract(aux,rho);
    // undo the k
    aux.applyLocalOperator(k,sigZ_Id);
    // and redo the L/2
    if(k<L-1) aux.applyLocalOperator(L/2,sigZ_Id);
    *out<<k-L/2<<"\t"<<real(valZ)<<"\t"<<real(valZZ)<<"\t";
  }
}



void computeLocalObservables(MPS& rho,ofstream* out){
  // Z operators
  static bool init(false);
  static mwArray sigX0,sigY0,sigZ0;
  if(!init){
    init=true;
    int d0=2;
    mwArray sig0=identityMatrix(d0);
    complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    complex_t datay[]={ZERO_c,I_c,-I_c,ZERO_c};
    complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
    constructOperatorProduct(sigX0,mwArray(Indices(d0,d0),datax),sig0);
    constructOperatorProduct(sigY0,mwArray(Indices(d0,d0),datay),sig0);
    constructOperatorProduct(sigZ0,mwArray(Indices(d0,d0),dataz),sig0);
  }
  int L=rho.getLength();
  Contractor& contractor=Contractor::theContractor();
  MPS aux(rho);
  aux.applyLocalOperator(L/2,sigX0);
  complex_t valX=contractor.contract(aux,rho);
  aux.applyLocalOperator(L/2,sigY0*sigX0);
  complex_t valY=contractor.contract(aux,rho);
  aux.applyLocalOperator(L/2,sigZ0*sigY0);
  complex_t valZ=contractor.contract(aux,rho);
  *out<<real(valX)<<"\t"<<real(valY)<<"\t"<<real(valZ)<<"\t";
}
 
