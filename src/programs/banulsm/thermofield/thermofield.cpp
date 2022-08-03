
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"

#include "ThermofieldHamiltonian.h"

using namespace std;
using namespace shrt;

/** Global matrices */
complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
mwArray sigX(Indices(2,2),dataX);
complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
mwArray sigY(Indices(2,2),dataY);
complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
mwArray sigZ(Indices(2,2),dataZ);

/** 
    Program thermofield constructs an initial state and lets it evolve
    Arguments for the couplings and free terms of the boson modes
    should be read from a file, where a mwArray is stored in text mode.
*/

// To try fourth order Trotter!
//#define ORDER4 1
// order 4 complex (faster?)
#define ORDER4_C 1

const string mpsTmpFileNameRoot(const string baseDir,int D,int L,double deltat,vector<int> allNs);
const string mpsTmpFileName(const string rootName,int cnt);


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
  int d = props.getIntProperty("d");  
  if(d!=2){
    cout<<"ERROR: Only d=2 supported for the central system"<<endl;
    exit(1);
  }
  double deltat=props.getDoubleProperty("TimeStep");
  int M = props.getIntProperty("NumberOfTimeSteps");  
  int D = props.getIntProperty("BondDimension");  
  int rate = props.getIntProperty("SavingRate");  

  string baseDir = props.getProperty("tmpMPSdir");
  int mpsrate=props.getIntProperty("SavingMPSRate");

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
  
  // Read the spin operators from lists of matrix elements indicated
  // in the Properties file
  vector<double> Hspin00=props.getVectorDoubleProperty("Hspin00");
  vector<double> Hspin01=props.getVectorDoubleProperty("Hspin01");
  vector<double> Hspin10=props.getVectorDoubleProperty("Hspin10");
  vector<double> Hspin11=props.getVectorDoubleProperty("Hspin11");

  vector<double> L00=props.getVectorDoubleProperty("L00");
  vector<double> L01=props.getVectorDoubleProperty("L01");
  vector<double> L10=props.getVectorDoubleProperty("L10");
  vector<double> L11=props.getVectorDoubleProperty("L11");

  string outfile=props.getProperty("OutputFile");

  // By hand I initialize the mwArrays now
  mwArray Hspin(Indices(d,d));
  mwArray OL(Indices(d,d));
  Hspin.setElement(Hspin00[0],Hspin00[1],Indices(0,0));
  Hspin.setElement(Hspin01[0],Hspin01[1],Indices(0,1));
  Hspin.setElement(Hspin10[0],Hspin10[1],Indices(1,0));
  Hspin.setElement(Hspin11[0],Hspin11[1],Indices(1,1));
  OL.setElement(L00[0],L00[1],Indices(0,0));
  OL.setElement(L01[0],L01[1],Indices(0,1));
  OL.setElement(L10[0],L10[1],Indices(1,0));
  OL.setElement(L11[0],L11[1],Indices(1,1));

  // Now construct the list of dimensions
  vector<int> dims(2*L+1,d);
  for(int k=0;k<L;k++){
    dims[L+k+1]=allNs[k]+1;
    dims[L-k-1]=allNs[k]+1;
  }

  // And initialize the Hamiltonian
  ThermofieldHamiltonian thH(L,dims,A1n,A2n,B1n,B2n,Hspin,OL);
  cout<<"Created the Hamiltonian"<<endl;


  // Now I should create the intial state, construct the evolution
  // operators, and start applying time steps
  MPS state(2*L+1,D,dims); 

  // First, create the root name for backup files
  string backupFileNameBasis=mpsTmpFileNameRoot(baseDir,D,L,deltat,allNs);
  // Now find, if any, the latest backup file stored
  string backupFile;
  bool found=false;int cnt=M;
  while((!found)&&cnt>0){
    backupFile=mpsTmpFileName(backupFileNameBasis,cnt);
    if(file_exists(backupFile)){
      found=true;
    }
    else{ backupFile="";cnt--;}
  }

  // Now, if found, I initialize the state with it, and the step with the corresponding count
  if(backupFile.length()!=0){
    cout<<"Using backup file "<<backupFile<<", step "<<cnt<<endl;
    state.importMPS(backupFile.data());
  }
  else{
    // otherwise, start from scratch
    cout<<"Starting from the beginning"<<endl;
    // Read the initial state
    vector<double> Psi0=props.getVectorDoubleProperty("Psi0");
    vector<double> Psi1=props.getVectorDoubleProperty("Psi1");

    // The initial state I construct by hand, now
    state.setProductState(p_zero); // everyone in 0
    // Change the state of the spin
    mwArray auxA(Indices(d,D,D));auxA.fillWithZero();
    auxA.setElement(Psi0[0],Psi0[1],Indices(0,0,0));
    auxA.setElement(Psi1[0],Psi1[1],Indices(1,0,0));
    state.setA(L,auxA);
    cout<<"Constructed initial state "<<endl;
  }

  //  output file prepared
  ofstream* out;
  if(found)
    out=new ofstream(outfile.data(),ios::app); // appending to old file
  else{
    out=new ofstream(outfile.data()); // brand new file
  }
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output"<<endl;
    exit(1);
  }
  cout<<"Writing output to files: "<<outfile<<endl;
  if(!found) // only if new file, write the header
    *out<<"% nstep\t delta\t time\t D\t <sigma_x>\t <sigma_y>\t <sigma_z>\t trace(rho)"<<endl;

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
  cout<<"Order 2 Trotter!"<<endl;
  MPO expHe(2*L+1),expHo(2*L+1),expHo_2(2*L+1);
  complex_t deltaC=-deltat*I_c;
  thH.getExponentialMPOeven(expHe,deltaC);
  cout<<"Constructed the exponential even MPO "<<endl;
  //  cout<<" With result "<<expHe<<endl;
  thH.getExponentialMPOodd(expHo,deltaC);
  thH.getExponentialMPOodd(expHo_2,.5*deltaC);
  cout<<"Constructed the exponential MPOs "<<endl;
  expHo_2.exportForMatlab("expHo2.m");
  //  cout<<" With result "<<expHo<<endl;
#endif
#endif

  Contractor& contractor=Contractor::theContractor();
  //int cnt=0;
  double time=cnt*deltat;

  mwArray rdm=contractor.getRDM(state,L);
  complex_t trRho=rdm.trace();
  mwArray calcuProd=sigX*rdm;
  complex_t sx=calcuProd.trace();
  calcuProd=sigY*rdm;
  complex_t sy=calcuProd.trace();
  calcuProd=sigZ*rdm; 
  complex_t sz=calcuProd.trace();

  *out<<cnt<<"\t"<<deltat<<"\t"<<time<<"\t"<<D<<"\t"<<sx<<"\t"<<sy
      <<"\t"<<sz<<"\t"<<trRho<<endl;
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
    // At multiples of mpsrate, save the MPs as tmp file
    if(cnt%mpsrate==0||cnt>=M){
      cout<<"Saving backup at cnt="<<cnt<<endl;
      string oldBackupFile(backupFile);
      backupFile=mpsTmpFileName(backupFileNameBasis,cnt);
      state.exportMPS(backupFile.data());
      // removing the old one
      if(file_exists(oldBackupFile.data())){
	cout<<"Removing old backup file "<<oldBackupFile<<endl;
	remove(oldBackupFile.data());
      }  
    }
    // Now compute and save observables
    // First, the RDM
    rdm=contractor.getRDM(state,L);
    // Its trace?
    trRho=rdm.trace();
    // Now each of the sigmas
    // rho is i'i as |phi><phi|
    calcuProd=sigX*rdm;
    sx=calcuProd.trace();
    calcuProd=sigY*rdm;
    sy=calcuProd.trace();
    calcuProd=sigZ*rdm; 
    sz=calcuProd.trace();
    *out<<cnt<<"\t"<<deltat<<"\t"<<time<<"\t"<<D<<"\t"<<sx<<"\t"<<sy
	<<"\t"<<sz<<"\t"<<trRho<<endl;
    cout<<"Step nr "<<cnt<<", t="<<time<<", sx="<<sx<<", trace="<<trRho<<endl;
  }
  out->close();
  delete out;

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

const string mpsTmpFileName(const string rootName,int cnt){
  stringstream s;
  s<<rootName<<cnt<<".dat";
  return s.str();
}
