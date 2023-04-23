
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"

#include "ThermofieldHamiltonian.h"

using namespace std;
using namespace shrt;

/** Global Pauli matrices */
int d=2;
mwArray sig0=identityMatrix(2);
complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
mwArray sigX(Indices(2,2),dataX);
complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
mwArray sigY(Indices(2,2),dataY);
complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
mwArray sigZ(Indices(2,2),dataZ);
mwArray sigP=.5*(sigX+I_c*sigY); //sigma^+
mwArray sigM=Hconjugate(sigP); //sigma^-

/** 
    Program thermofieldQubitFluctuator constructs an initial state and lets it evolve
    Arguments for the couplings and free terms of the boson modes
    should be read from a file, where a mwArray is stored in text mode.

    The system Hamiltonian (qubit+fluctuator) is of the form:
    \f[H_{\mathrm{QF}}=\frac{1}{2}\omega_{\mathrm{q}} \sigma_{\mathrm{q}}^z + \frac{1}{2}\epsilon \sigma_{\mathrm{f}}^z 
    +\frac{1}{2}\Delta \sigma_{\mathrm{f}}^x+
    \frac{1}{2}v \sigma_{\mathrm{q}}^z\otimes  \sigma_{\mathrm{f}}^z
    \f]
    and the parameters are read from the properties file (or command line).

    The interaction of the system with the bath is only through the fluctuator, with a term

    \f[
    \sigma_{\mathrm{f}}^{+} b_0+ \mathrm{h.c.}
    \f]
    where \f$b_0\f$ is the first mode of the chain.

    Additionally, DD pulses are applied at regular intervals of time, given by parameter \param<timeDD>
*/

const string mpsTmpFileNameRoot(const string baseDir,int D,int L,double deltat,vector<int> allNs);
const string mpsTmpFileName(const string rootName,int label,int cnt);

// Evolution according to different Trotter orders
void evolveSteps(int trotter,MPS& state,int D,const ThermofieldHamiltonian& thH,double deltat,int nr);

void evolveStepsTrotter2(MPS& state,int D,const ThermofieldHamiltonian& thH,double deltat,int nr);
// To try fourth order Trotter
void evolveStepsTrotter4(MPS& state,int D,const ThermofieldHamiltonian& thH,double deltat,int nr);
// order 4 complex (faster?)
void evolveStepsTrotter4C(MPS& state,int D,const ThermofieldHamiltonian& thH,double deltat,int nr);


// auxiliary to write results to text file
void writeComponents(ofstream* out,const mwArray& rdm);

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

  int L = props.getIntProperty("NumberOfBosonicModes");  
  double deltat=props.getDoubleProperty("TimeStep");
  int M = props.getIntProperty("NumberOfTimeSteps");  
  int D = props.getIntProperty("BondDimension");  
  int rate = props.getIntProperty("SavingRate");  
  int stepDD=props.getIntProperty("stepDD");
  bool applyingDD=stepDD>0;
  // I need tDD to be integer multiple of deltat, so I read an integer
  // rate cannot be larger than stepDD (it has to be a multiple), so I fix, if necessary
  if(applyingDD){
    if(rate>stepDD){
      rate=stepDD;
    }
    else{
      rate=gcd(rate,stepDD);
      cout<<"Modified SavingRate to "<<rate<<" for compatiibility with DD (every "<<stepDD<<" steps)"<<endl;
    }
  }

  
  string baseDir = props.getProperty("tmpMPSdir");
  int mpsrate=props.getIntProperty("SavingMPSRate");

  int trotter=1; //default is 2
  bool trotter2=(props.getIntProperty("Trotter2")!=0);
  if(trotter2){
    cout<<"Setting Trotter order to 2"<<endl;
    trotter=1;
  }
  bool trotter4=(props.getIntProperty("Trotter4")!=0);
  if(trotter4){
    cout<<"Setting Trotter order to 4"<<endl;
    trotter=2;
  }
  bool trotter4C=(props.getIntProperty("Trotter4C")!=0);
  if(trotter4C){
    cout<<"Setting Trotter order to 4C"<<endl;
    trotter=3;
  }
  
  vector<int> allNs=props.getVectorIntProperty("MaxOccupationOfBosonicModes");
  if(allNs.size()<L){
    int sz=allNs.size();
    int lastN=allNs.back();
    for(int k=sz;k<L;k++)
      allNs.push_back(lastN);
  }

  bool fileCoeff = (props.getIntProperty("ReadCoefficientsFromFiles")==1);  
  vector<double> A1n,A2n,B1n,B2n;
  if(!fileCoeff){
    cout<<"WARNING!!! No file for coefficients of the chain(s)??"<<endl;
    // Then I read from Properties
    B2n=props.getVectorDoubleProperty("B2n"); // coefs of the couplings (R)
    B1n=props.getVectorDoubleProperty("B1n"); // idem (L)
    A2n=props.getVectorDoubleProperty("A2n"); // coeffs diagonal boson
					      // terms (R)
    A1n=props.getVectorDoubleProperty("A1n"); // idem (L)
    if(A1n.size()<L||A2n.size()<L||B1n.size()<L||B2n.size()<L){
      cout<<"ERROR: Not enough coefficients given in the properties file"<<endl;
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

  // Coefficients for the qubit+fluctuator Hamiltonian
  double omegaq=props.getDoubleProperty("omegaq");  
  double epsilon=props.getDoubleProperty("epsilon");  
  double Delta=props.getDoubleProperty("Delta");  
  double v=props.getDoubleProperty("v");  

  string qubitfile=props.getProperty("FileQubitState"); // An array, in text mode, with the initial state of the qubit (as two complex numbers, written as the sequence 3 2 1 1 real,imag,real,imag

  
  // construct the local Hamiltonian and interaction operator for the q+f system

  mwArray Hspin;
  {
    constructOperatorProduct(Hspin,.5*omegaq*sigZ,sig0); // free qubit
    mwArray auxH;
    constructOperatorProduct(auxH,sig0,.5*epsilon*sigZ+.5*Delta*sigX); // free fluctuator
    Hspin=Hspin+auxH;
    constructOperatorProduct(auxH,.5*v*sigZ,sigZ); // qubit-fluctuator
    Hspin=Hspin+auxH;
  }

  // Same for interaction operator, defined as the term that goes with boson destruction
  // b0, so it has to be sigma^+ on the fluctuator
  mwArray OL;
  constructOperatorProduct(OL,sig0,sigP);

  //  cout<<"sigP="<<sigP<<endl;
  // cout<<"OL="<<OL<<endl;

  // if DD, need the corresponding operator (I assume sigmax on qubit)
  mwArray opDD;
  if(applyingDD){
    constructOperatorProduct(opDD,sigX,sig0);
  }
  
  
  string outfile=props.getProperty("OutputFile");

  // Now construct the list of dimensions
  vector<int> dims(2*L+1,d*d);
  for(int k=0;k<L;k++){
    dims[L+k+1]=allNs[k]+1;
    dims[L-k-1]=allNs[k]+1;
  }

  // And initialize the Hamiltonian
  ThermofieldHamiltonian thH(L,dims,A1n,A2n,B1n,B2n,Hspin,OL);
  cout<<"Created the Hamiltonian"<<endl;


  // Now I should create the initial state, construct the evolution
  // operators, and start applying time steps

  // Because the fluctuator is in a thermal state, I will evolve each of its eigenvectors
  // independently=> there will be two MPS

  MPS state1(2*L+1,D,dims); 
  MPS state2(2*L+1,D,dims); 

  // First, create the root name for backup files
  string backupFileNameBasis=mpsTmpFileNameRoot(baseDir,D,L,deltat,allNs);
  // Now find, if any, the latest backup file stored
  string backupFile1,backupFile2;
  bool found=false;int cnt=M;
  while((!found)&&cnt>0){
    backupFile1=mpsTmpFileName(backupFileNameBasis,1,cnt);
    if(file_exists(backupFile1)){
      found=true;
    }
    else{ backupFile1="";cnt--;}
  }

  // Now, if found, I initialize the state with it, and the step with the corresponding count
  bool readFile=0;
  if(backupFile1.length()!=0){
    cout<<"Found backup file for MPS 1 "<<backupFile1<<", step "<<cnt<<endl;
    state1.importMPS(backupFile1.data());
    // check for the second one too
    backupFile2=mpsTmpFileName(backupFileNameBasis,2,cnt);
    if(file_exists(backupFile2)){
      cout<<"Found backup file for MPS 2 "<<backupFile2<<", step "<<cnt<<endl;
      state2.importMPS(backupFile1.data());
      readFile=true;
    }
    else{
      cout<<"Backup files for MPS1 and 2 do not correspond to same step (or one not found)=> restart"<<endl;
      readFile=false;state1.clear();
    }
  }
  if(!readFile){
    // otherwise, start from scratch
    cout<<"Starting from the beginning"<<endl;
    // The initial state I construct by hand, now
    state1=MPS(2*L+1,1,dims);
    state1.setProductState(p_zero); // everyone in 0
    state2=MPS(2*L+1,1,dims);
    state2.setProductState(p_zero); // everyone in 0

    // Change the state of the system to each eigenstate of the fluctuator,
    // times the qubit in a fixed state, read from file:
    ifstream inqubit(qubitfile.data());
    mwArray psiQ;psiQ.loadtext(inqubit,false); // should be dims 2,1,1
    inqubit.close();
    psiQ.reshape(Indices(2,1));
    double normQ=real((Hconjugate(psiQ)*psiQ).getElement(Indices(0,0)));
    psiQ=(1./sqrt(normQ))*psiQ;
    cout<<"Read initial qubit state: "<<psiQ<<endl;
    
    // eigenstates of the fluctuator: psiF1 the negative one
    mwArray psiF1(Indices(2,1,1));
    mwArray psiF2(Indices(2,1,1));
    if(epsilon!=0.){
      double x=Delta/epsilon;
      if(x==0){
	psiF1.setElement(ZERO_c,Indices(0,0,0)); // comp 0
	psiF1.setElement(ONE_c,Indices(1,0,0)); // comp 1
	
	psiF2.setElement(ONE_c,Indices(0,0,0)); // comp 0
	psiF2.setElement(ZERO_c,Indices(1,0,0)); // comp 1

      }
      else{
	psiF1.setElement((1-sqrt(1+x*x))*ONE_c,Indices(0,0,0)); // comp 0
	psiF1.setElement(x*ONE_c,Indices(1,0,0)); // comp 1
	
	psiF2.setElement((1+sqrt(1+x*x))*ONE_c,Indices(0,0,0)); // comp 0
	psiF2.setElement(x*ONE_c,Indices(1,0,0)); // comp 1
	double norm=x*x+pow(1-sqrt(1+x*x),2);
	psiF1=(1./sqrt(norm))*psiF1;
	norm=x*x+pow(1+sqrt(1+x*x),2);
	psiF2=(1./sqrt(norm))*psiF2;
      }
    }
    else{ // special case epsilon=0
      cout<<"Special case: epsilon=0, thus eigenstates of F are x+-"<<endl;
      psiF1.setElement((1./sqrt(2))*ONE_c,Indices(0,0,0)); // comp 0
      psiF1.setElement(-(1./sqrt(2))*ONE_c,Indices(1,0,0)); // comp 1

      psiF2.setElement((1./sqrt(2))*ONE_c,Indices(0,0,0)); // comp 0
      psiF2.setElement((1./sqrt(2))*ONE_c,Indices(1,0,0)); // comp 1      
    }
    cout<<"States of the fluctuator psiF1="<<psiF1<<"\n\t psiF2="<<psiF2<<endl;

    // now construct the state of Qbit times fluctuator
    psiF1.reshape(Indices(1,2));
    psiF1.multiplyLeft(psiQ); // 2x2
    psiF2.reshape(Indices(1,2));
    psiF2.multiplyLeft(psiQ); // 2x2

    psiF1.reshape(Indices(2*2,1,1));
    psiF2.reshape(Indices(2*2,1,1));
    
    state1.setA(L,psiF1);
    state2.setA(L,psiF2);
    cout<<"Constructed initial states "<<endl;
  }

  //  output file prepared
  ofstream* out;
  if(!found){
    out=new ofstream(outfile.data()); // brand new file
    *out<<"% nstep\t delta\t time\t D\t rho1_11\t rho1_22\t real(rho1_12)\t imag(rho1_12) \t rho2_11\t rho2_22\t real(rho2_12)\t imag(rho2_12)"<<endl;
  }
  else{
    out=new ofstream(outfile.data(),ios::app); // appending to old file
  }
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output (check permissions?)"<<endl;
    exit(1);
  }
  out->close(); delete out;
  cout<<"Writing output to file: "<<outfile<<endl;
    
  Contractor& contractor=Contractor::theContractor();
  //int cnt=0;
  double time=cnt*deltat;

  mwArray rdm1=contractor.getRDM(state1,L);
  cout<<"Test: at the very beginning, e.v. of system H in state1 "<<(Hspin*rdm1).trace()<<endl;
  
  // trace out the fluctuator, too
  rdm1.reshape(Indices(2,2,2,2));
  rdm1.permute(Indices(1,3,2,4));
  rdm1.reshape(Indices(2*2,2*2));
  rdm1.multiplyRight(reshape(identityMatrix(2),Indices(2*2,1)));
  rdm1.reshape(Indices(2,2));
  
  mwArray rdm2=contractor.getRDM(state2,L);
  cout<<"Test: at the very beginning, e.v. of system H in state2: "<<(Hspin*rdm2).trace()<<endl;
  rdm2.reshape(Indices(2,2,2,2));
  rdm2.permute(Indices(1,3,2,4));
  rdm2.reshape(Indices(2*2,2*2));
  rdm2.multiplyRight(reshape(identityMatrix(2),Indices(2*2,1)));
  rdm2.reshape(Indices(2,2));

  // I will store the diagonal elements of each rho, and the off-diagonal (complex) one, total 4 real parameters

  out=new ofstream(outfile.data(),ios::app); // appending 
  *out<<setprecision(10);
  *out<<cnt<<"\t"<<deltat<<"\t"<<time<<"\t"<<D<<"\t";
  writeComponents(out,rdm1);
  writeComponents(out,rdm2);
  *out<<endl;
  out->close();delete out;
  
  MPS aux(state1);
  while(cnt<M){
    int nrToApply=min(rate,M-cnt);
    //    cout<<"Step nr: "<<cnt<<", applying "<<nrToApply<<endl;
    // Evolve state1 and state2 nrToApply steps
    evolveSteps(trotter,state1,D,thH,deltat,nrToApply);
    //cout<<"Testing after application of "<<nrToApply<<" steps, overlap of state1: ("<<time+nrToApply*deltat<<","<<contractor.contract(state1,aux)<<")"<<endl;
    evolveSteps(trotter,state2,D,thH,deltat,nrToApply);
    //    cout<<"\t overlap of state2: ("<<time+nrToApply*deltat<<","<<contractor.contract(state1,state2)<<")"<<endl;
    cnt+=nrToApply;time+=nrToApply*deltat;

    // If needed, apply DD pulse
    if(applyingDD){
      if(cnt>0&&cnt%stepDD==0){
	cout<<"Applying DD pulse after step "<<cnt<<endl;
	state1.applyLocalOperator(L,opDD,false);
	state2.applyLocalOperator(L,opDD,false); 
      }
    }
    
    // Although norm should not change under unitary evolution, it can
    // do it due to truncation, so I normalize again
    state1.gaugeCond('R',1);
    state2.gaugeCond('R',1);
    // At multiples of mpsrate, save the MPs as tmp file
    if(cnt%mpsrate==0||cnt>=M){
      cout<<"Saving backup at cnt="<<cnt<<endl;
      string oldBackupFile1(backupFile1);
      string oldBackupFile2(backupFile2);
      backupFile1=mpsTmpFileName(backupFileNameBasis,1,cnt);
      backupFile2=mpsTmpFileName(backupFileNameBasis,2,cnt);
      state1.exportMPS(backupFile1.data());
      state2.exportMPS(backupFile2.data());
      // removing the old one
      if(file_exists(oldBackupFile1.data())){
	cout<<"Removing old backup file "<<oldBackupFile1<<endl;
	remove(oldBackupFile1.data());
      }  
      if(file_exists(oldBackupFile2.data())){
	cout<<"Removing old backup file "<<oldBackupFile2<<endl;
	remove(oldBackupFile2.data());
      }  
    }
    // Now save observables
      mwArray rdm1=contractor.getRDM(state1,L);
      // trace out the fluctuator, too
      rdm1.reshape(Indices(2,2,2,2));
      rdm1.permute(Indices(1,3,2,4));
      rdm1.reshape(Indices(2*2,2*2));
      rdm1.multiplyRight(reshape(identityMatrix(2),Indices(2*2,1)));
      rdm1.reshape(Indices(2,2));
      
      mwArray rdm2=contractor.getRDM(state2,L);
      rdm2.reshape(Indices(2,2,2,2));
      rdm2.permute(Indices(1,3,2,4));
      rdm2.reshape(Indices(2*2,2*2));
      rdm2.multiplyRight(reshape(identityMatrix(2),Indices(2*2,1)));
      rdm2.reshape(Indices(2,2));

      out=new ofstream(outfile.data(),ios::app); // appending 
      *out<<setprecision(10);
      *out<<cnt<<"\t"<<deltat<<"\t"<<time<<"\t"<<D<<"\t";
      writeComponents(out,rdm1);
      writeComponents(out,rdm2);
      *out<<endl;
      out->close();
      delete out;
  }
}


const string mpsTmpFileNameRoot(const string baseDir,int D,int L,double deltat,vector<int> allNs){
  stringstream s;
  s<<baseDir<<"/tmpMPS_D"<<D<<"_L"<<L<<"_deltat"<<deltat;
  s<<"_Nb";
  for(int k=0;k<allNs.size();k++)
    s<<allNs[k];
  s<<"_";
  return s.str();
}

const string mpsTmpFileName(const string rootName,int label,int cnt){
  stringstream s;
  s<<rootName<<label<<"_"<<cnt<<".dat";
  return s.str();
}

void writeComponents(ofstream* out,const mwArray& rdm){
  complex_t trRho=rdm.trace();
  complex_t aux=rdm.getElement(Indices(0,0));
  *out<<real(aux)<<"\t";
  aux=rdm.getElement(Indices(1,1));
  *out<<real(aux)<<"\t";
  aux=rdm.getElement(Indices(0,1));
  *out<<real(aux)<<"\t"<<imag(aux)<<"\t";
  //*out<<trRho<<"\t";
}


void evolveStepsTrotter2(MPS& state,int D,const ThermofieldHamiltonian& thH,double deltat,int nr){
  static bool init=false;
  int L=state.getLength(); //2*L+1
  static MPO expHe(2*L+1);
  static MPO expHo(2*L+1);
  static MPO expHo_2(2*L+1);
  static complex_t deltaC=-deltat*I_c;
  if(!init||deltaC!=(-deltat*I_c)){
    cout<<"Initializing Order 2 Trotter!"<<endl;
    deltaC=-deltat*I_c;
    thH.getExponentialMPOeven(expHe,deltaC);
    //cout<<"Constructed the exponential even MPO "<<endl;
    thH.getExponentialMPOodd(expHo,deltaC);
    thH.getExponentialMPOodd(expHo_2,.5*deltaC);
    cout<<"Constructed the exponential MPOs "<<endl;
    // expHe.exportForMatlab("expHe.m");
    // expHo.exportForMatlab("expHo.m");
    // expHo_2.exportForMatlab("expHo_2.m");
    init=true;
  } // now ready to apply

  MPS aux(state);
  Contractor& contractor=Contractor::theContractor();
  contractor.optimize(expHo_2,aux,state,D);
  while(nr>0){
    contractor.optimize(expHe,state,aux,D);
    if(nr>1)
      contractor.optimize(expHo,aux,state,D);
    else
      contractor.optimize(expHo_2,aux,state,D);
    nr--;
    //cout<<"\tTo apply: "<<nr<<endl;
  }
}


void evolveStepsTrotter4(MPS& state,int D,const ThermofieldHamiltonian& thH,double deltat,int nr){
  static bool init=false;
  int L=state.getLength(); //2*L+1
  // For Trotter fourth order I need 5 different operators
  static MPO expHe_s(2*L+1);
  static MPO expHe_2s(2*L+1);
  static MPO expHo_s(2*L+1);
  static MPO expHo_ms(2*L+1);
  static MPO expHo_s_2(2*L+1);
  static double sconst=1.351207191959657;
  static complex_t deltaC=-deltat*I_c;
  if(!init||deltaC!=(-deltat*I_c)){
    cout<<"Initializing Order 4 Trotter!"<<endl;
    deltaC=-deltat*I_c;
    thH.getExponentialMPOeven(expHe_s,sconst*deltaC);
    thH.getExponentialMPOeven(expHe_2s,(1-2*sconst)*deltaC);
    thH.getExponentialMPOodd(expHo_s,sconst*deltaC);
    thH.getExponentialMPOodd(expHo_s_2,.5*sconst*deltaC);
    thH.getExponentialMPOodd(expHo_ms,.5*(1-sconst)*deltaC);
    init=true;
  } // now ready to apply
  MPS aux(state);
  Contractor& contractor=Contractor::theContractor();
  contractor.optimize(expHo_s_2,aux,state,D);
  while(nr>0){
    contractor.optimize(expHe_s,state,aux,D);
    contractor.optimize(expHo_ms,aux,state,D);
    contractor.optimize(expHe_2s,state,aux,D);
    contractor.optimize(expHo_ms,aux,state,D);
    contractor.optimize(expHe_s,state,aux,D);
    if(nr>1)
      contractor.optimize(expHo_s,aux,state,D);
    else
      contractor.optimize(expHo_s_2,aux,state,D);
    nr--;
  }
}

void evolveStepsTrotter4C(MPS& state,int D,const ThermofieldHamiltonian& thH,double deltat,int nr){
  static bool init=false;
  int L=state.getLength(); //2*L+1
  // For Trotter fourth order I need 5 different operators
  complex_t Apar=(.5*ONE_c+(sqrt(3)/6.)*I_c);
  static MPO expHe_a(2*L+1);
  static MPO expHe_aC(2*L+1);
  static MPO expHo_aC_2(2*L+1);
  static MPO expHo_a_2(2*L+1);
  static MPO expHo_2(2*L+1);
  static complex_t deltaC=-deltat*I_c;
  if(!init||deltaC!=(-deltat*I_c)){
    cout<<"Initializing Order 4 Trotter!"<<endl;
    deltaC=-deltat*I_c;
    thH.getExponentialMPOeven(expHe_a,Apar*deltaC);
    thH.getExponentialMPOeven(expHe_aC,conjugate(Apar)*deltaC);
    thH.getExponentialMPOodd(expHo_a_2,.5*Apar*deltaC);
    thH.getExponentialMPOodd(expHo_2,.5*deltaC);
    thH.getExponentialMPOodd(expHo_aC_2,.5*conjugate(Apar)*deltaC);
    init=true;
  } // now ready to apply
  MPS aux(state);
  Contractor& contractor=Contractor::theContractor();
  contractor.optimize(expHo_aC_2,aux,state,D);
  while(nr>0){
    contractor.optimize(expHe_aC,state,aux,D);
    contractor.optimize(expHo_2,aux,state,D);
    contractor.optimize(expHe_a,state,aux,D);
    if(nr>1)
      contractor.optimize(expHo_2,aux,state,D);
    else
      contractor.optimize(expHo_a_2,aux,state,D);
    nr--;
  }
}


void evolveSteps(int trotter,MPS& state,int D,const ThermofieldHamiltonian& thH,double deltat,int nr){
  switch(trotter){
  case 1:
    evolveStepsTrotter2(state,D,thH,deltat,nr);
    break;
  case 2:
    evolveStepsTrotter4(state,D,thH,deltat,nr);
    break;
  case 3:
    evolveStepsTrotter4C(state,D,thH,deltat,nr);
    break;
  default:
    cout<<"ERROR: Unknown Trotter order"<<endl;
    exit(1);
  }
  
}
