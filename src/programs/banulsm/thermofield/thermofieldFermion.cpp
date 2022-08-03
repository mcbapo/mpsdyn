
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"

#include "ThermofieldFermionHamiltonian.h"

using namespace std;
using namespace shrt;

/** Global matrices */
complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
mwArray sigX(Indices(2,2),dataX);
complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
mwArray sigY(Indices(2,2),dataY);
complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
mwArray sigZ(Indices(2,2),dataZ);
complex_t dataN[]={ONE_c,ZERO_c,ZERO_c,ZERO_c};
mwArray spinN(Indices(2,2),dataN); // n_mode

/** 
    Program thermofield constructs an initial state and lets it evolve
    Arguments for the couplings and free terms of the boson modes
    should be read from a file, where a mwArray is stored in text mode.
*/

// To try fourth order Trotter!
//#define ORDER4 1
// order 4 complex (faster?)
//#define ORDER4_C 1

int main(int argc,const char* argv[]){
  int cntr=0;
  const char* infile=argv[++cntr];
  string directory="";
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  int L = props.getIntProperty("NumberOfBathModes");  
  int d = props.getIntProperty("d");  
  if(d!=2){
    cout<<"ERROR: Only d=2 supported for the central system"<<endl;
    exit(1);
  }
  double deltat=props.getDoubleProperty("TimeStep");
  int M = props.getIntProperty("NumberOfTimeSteps");  
  int D = props.getIntProperty("BondDimension");  
  int rate = props.getIntProperty("SavingRate");  

  bool fileCoeff = (props.getIntProperty("ReadCoefficientsFromFiles")==1);  
  vector<double> A1n,A2n,B1n,B2n;
  if(!fileCoeff){
    // Then I read from Properties
    B2n=props.getVectorDoubleProperty("B2n"); // coefs of the couplings (R)
    B1n=props.getVectorDoubleProperty("B1n"); // idem (L)
    A2n=props.getVectorDoubleProperty("A2n"); // coeffs diagonal boson
					      // terms (R)
    A1n=props.getVectorDoubleProperty("A1n"); // idem (L)
    if(A1n.size()<L||A2n.size()<L||B1n.size()<L||B2n.size()<L){
      cout<<"ERROR: Not enough coefficients given int he properties file"<<endl;
      exit(1);
    }
  }
  else{ // read from files
    string AB1file=props.getProperty("FileCoefficients1");
    cout<<"Reading coefficients from file "<<AB1file<<endl;
    // read A1n and B2n from here
    ifstream inlog(AB1file.data());
    int ptr=0;
    while(inlog.good()&&ptr<L){
      double alpha,beta;
      inlog>>alpha;
      inlog>>beta;
      //cout<<"Read coeffs ("<<ptr<<") "<<alpha<<", "<<beta<<endl;
      A1n.push_back(alpha);
      B1n.push_back(beta);
      ptr++;
    }
    string AB2file=props.getProperty("FileCoefficients2");
    cout<<"Reading coefficients from file "<<AB2file<<endl;
    ifstream inlog2(AB2file.data());
    ptr=0;
    while(inlog2.good()&&ptr<L){
      double alpha,beta;
      inlog2>>alpha;
      inlog2>>beta;
      //cout<<"Read coeffs ("<<ptr<<") "<<alpha<<", "<<beta<<endl;
      A2n.push_back(alpha);
      B2n.push_back(beta);
      ptr++;
    }
  }
  double U=props.getDoubleProperty("U");
  double V=props.getDoubleProperty("V");

  double couplingT=props.getDoubleProperty("t"); // coupling
  if(couplingT==-1.) couplingT=1.; // should multiply the first B terms
  B1n[0]*=couplingT;
  B2n[0]*=couplingT;
  
  string outfile=props.getProperty("OutputFile");

  // Initialize the Hamiltonian
  ThermofieldFermionHamiltonian thH(L,A1n,A2n,B1n,B2n,U,V);
  cout<<"Created the Hamiltonian"<<endl;


  // Now I should create the intial state, construct the evolution
  // operators, and start applying time steps

  // Read the initial state
  vector<double> Psi00=props.getVectorDoubleProperty("Psi00");
  vector<double> Psi01=props.getVectorDoubleProperty("Psi01");
  vector<double> Psi10=props.getVectorDoubleProperty("Psi10");
  vector<double> Psi11=props.getVectorDoubleProperty("Psi11");

  vector<int> dims(2*L+1,d);
  dims[L]=4;
  // The initial state I construct by hand, now
  MPS state(2*L+1,D,dims); 
  //  state.setProductState(p_zero); // everyone in 0
  state.setProductState(p_one); // everyone in 1, as this is equivalent to fermionic vacuum
  // Change the state of the spin
  mwArray auxA(Indices(d,d,D,D));auxA.fillWithZero();
  auxA.setElement(Psi00[0],Psi00[1],Indices(0,0,0,0));
  auxA.setElement(Psi01[0],Psi01[1],Indices(0,1,0,0));
  auxA.setElement(Psi10[0],Psi10[1],Indices(1,0,0,0));
  auxA.setElement(Psi11[0],Psi11[1],Indices(1,1,0,0));
  auxA.reshape(Indices(d*d,D,D));
  state.setA(L,auxA);
  cout<<"Constructed initial state "<<endl;

  // NOw get the exponentials I need
#ifdef ORDER4
  cout<<"Running fourth Trotter order!"<<endl;
  // For Trotter fourth order I need 5 different operators
  MPO expHe_s(2*L+1),expHe_2s(2*L+1),expHo_s(2*L+1),expHo_ms(2*L+1),expHo_s_2(2*L+1);
  double sconst=1.351207191959657;
  complex_t deltaC=-deltat*I_c;
  thH.getExponentialMPOeven(expHe_s,sconst*deltaC);
  thH.getExponentialMPOeven(expHe_2s,(1-2*sconst)*deltaC);
  thH.getExponentialMPOodd(expHo_s,sconst*deltaC);
  thH.getExponentialMPOodd(expHo_s_2,.5*sconst*deltaC);
  thH.getExponentialMPOodd(expHo_ms,.5*(1-sconst)*deltaC);
#else
#ifdef ORDER4_C
  cout<<"Running fourth Trotter order with complex parameters!"<<endl;
  complex_t Apar=(.5*ONE_c+(sqrt(3)/6.)*I_c);
  // For Trotter fourth order I need 5 different operators
  MPO expHe_a(2*L+1),expHe_aC(2*L+1),expHo_aC_2(2*L+1),expHo_a_2(2*L+1),expHo_2(2*L+1);
  complex_t deltaC=-deltat*I_c;
  thH.getExponentialMPOeven(expHe_a,Apar*deltaC);
  thH.getExponentialMPOeven(expHe_aC,conjugate(Apar)*deltaC);
  thH.getExponentialMPOodd(expHo_a_2,.5*Apar*deltaC);
  thH.getExponentialMPOodd(expHo_2,.5*deltaC);
  thH.getExponentialMPOodd(expHo_aC_2,.5*conjugate(Apar)*deltaC);
#else
  cout<<"Running order 2 Trotter"<<endl;
  MPO expHe(2*L+1),expHo(2*L+1),expHo_2(2*L+1);
  complex_t deltaC=-deltat*I_c;
  thH.getExponentialMPOeven(expHe,deltaC);
  cout<<"Constructed the exponential even MPO "<<endl;
  //cout<<" With result "<<expHe<<endl;
  thH.getExponentialMPOodd(expHo,deltaC);
  cout<<"Constructed the exponential odd MPO "<<endl;
  //cout<<" With result "<<expHo<<endl;
  thH.getExponentialMPOodd(expHo_2,.5*deltaC);
  cout<<"Constructed the exponential MPOs "<<endl;
  //expHo_2.exportForMatlab("expHo2.m");
  //  cout<<" With result "<<expHo<<endl;
#endif
#endif

  // Prepare the observables for the system (occupation of each fermionic mode)
  mwArray Ndown(spinN);
  Ndown.reshape(Indices(d*d,1));
  Ndown.multiplyRight(reshape(identityMatrix(2),Indices(1,4)));
  Ndown.reshape(Indices(d,d,d,d));
  mwArray Nup(Ndown);
  Ndown.permute(Indices(1,3,2,4));
  Ndown.reshape(Indices(d*d,d*d));
  Nup.permute(Indices(3,1,4,2));
  Nup.reshape(Indices(d*d,d*d));


  Contractor& contractor=Contractor::theContractor();
  int cnt=0;double time=0;

  mwArray rdm=contractor.getRDM(state,L);
  complex_t trRho=rdm.trace();
  mwArray calcuProd=Ndown*rdm;
  complex_t nd=calcuProd.trace();
  calcuProd=Nup*rdm;
  complex_t nu=calcuProd.trace();

  //  ofstream* out=new ofstream(outfile.data(),ios::app);
  ofstream* out=new ofstream(outfile.data());
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output"<<endl;
    exit(1);
  }
  cout<<"Writing output to files: "<<outfile<<endl;
  *out<<"% Measured the occupation number of both fermionic modes (called up and down)"<<endl;
  *out<<"% nstep\t delta\t time\t D\t Nup\t Ndown\t trace(rho)"<<endl;
  *out<<cnt<<"\t"<<deltat<<"\t"<<time<<"\t"<<D<<"\t"<<nd<<"\t"<<nu
      <<"\t"<<trRho<<endl;
  while(cnt<M){
    //cout<<"Step nr: "<<cnt<<endl;
    MPS aux(state);
#ifdef ORDER4
    contractor.optimize(expHo_s_2,aux,state,D);
#else
#ifdef ORDER4_C
    contractor.optimize(expHo_aC_2,aux,state,D);
#else
    contractor.optimize(expHo_2,aux,state,D);
#endif
#endif
    // Now apply r-1 times the pair Ho He
    int cntLoop=0;
    while(cntLoop<rate&&cnt<M){
#ifdef ORDER4
      contractor.optimize(expHe_s,state,aux,D);
      contractor.optimize(expHo_ms,aux,state,D);
      contractor.optimize(expHe_2s,state,aux,D);
      contractor.optimize(expHo_ms,aux,state,D);
      contractor.optimize(expHe_s,state,aux,D);
#else
#ifdef ORDER4_C
      contractor.optimize(expHe_aC,state,aux,D);
      contractor.optimize(expHo_2,aux,state,D);
      contractor.optimize(expHe_a,state,aux,D);
#else
      contractor.optimize(expHe,state,aux,D);
#endif
#endif
      if(cntLoop<rate-1){
#ifdef ORDER4
	contractor.optimize(expHo_s,aux,state,D);
#else
#ifdef ORDER4_C
	contractor.optimize(expHo_2,aux,state,D);
#else
	contractor.optimize(expHo,aux,state,D);
#endif
#endif
      }
      else{
#ifdef ORDER4
	contractor.optimize(expHo_s_2,aux,state,D);
#else
#ifdef ORDER4_C
	contractor.optimize(expHo_a_2,aux,state,D);
#else
	contractor.optimize(expHo_2,aux,state,D);
#endif
#endif
      }
      cnt++;cntLoop++;time+=deltat;
    }
    // Although norm should not change under unitary evolution, it can
    // do it due to truncation, so I normalize again
    state.gaugeCond('R',1);
    // Now compute and save observables
    // First, the RDM
    rdm=contractor.getRDM(state,L);
    // Its trace?
    trRho=rdm.trace();
    // Now occupation of each fermionic mode (1+sigz)/2 for the corresponding mode (left=down)

    // rho is i'i as |phi><phi|
    calcuProd=Ndown*rdm;
    nd=calcuProd.trace();
    calcuProd=Nup*rdm;
    nu=calcuProd.trace();
    *out<<cnt<<"\t"<<deltat<<"\t"<<time<<"\t"<<D<<"\t"<<nd<<"\t"<<nu
	<<"\t"<<trRho<<endl;
    cout<<"Step nr "<<cnt<<", t="<<time<<", nu="<<nu<<", trace="<<trRho<<endl;
  }
  out->close();
  delete out;

}
