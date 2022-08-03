/**
   \file SU2NoLinksSpectrum.cpp
   
   This program is designed for large production runs to ultimately determine the spectrum of the SU(2) model in the formulation where the links are no longer explicitly present. The aim is to compute MPS approximations for the vector and the scalar and to determine in an automated way if they have already been found such that I do not have to take care of that myself. Furthermore options can be overridden to manually take control over the number of excited states which should be computed. 
   
   Additionally there is quite some error checking going on to make automated productions runs safe and easy to restart.
   
   
  \param <inputfile> (string) File containing the parameters for the run
  
  \author Stefan Kühn
  \date 17/11/2015

*/


#include <cstdlib>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>
//Lib which is needed for sleep command on *nix Systems
#include "unistd.h"

#include "MPS.h"
#include "MPO.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include "Indices.h"
#include "misc.h"

#include "Properties.h"
#include "SU2HamiltonianNoLinks.h"
#include "SU2HamiltonianNoLinksQ.h"
#include "SU2HamiltonianNoLinksDimless.h"
#include "SU2HamiltonianNoLinksDimlessQ.h"

#define NUM_OF_PARAMS 1
#define FILENAME "SU2nolinks.txt"
#define CONFIGFILE "SU2nolinksconfiguration.txt"
#define TMPFILE "tmp.dat"
#define LOGFILE "log.txt"


//Maxtime (currently 2 days)
#define MAXTIME  172800
//Maximum number of excited states to be computed
#define MAXSTATES 100

using namespace std;
using namespace shrt;

/** Get spin and flux configuration*/
vector<complex_t> get_full_config(vector<MPO*> NocOps, vector<MPO*> FluxOps, Contractor &contractor, MPS &mps, int N);
/** Get charge configuration*/
vector<complex_t> get_config(vector<MPO*> Ops, Contractor &contractor, MPS &mps, int N);
/**Save the current configuration to a given file*/
int save_configuration(string filename, vector<complex_t> config);
/** Compute the sum of the absolute values of the vector entries*/
double sum_abs(vector<complex_t> input);
/**Given the parameters, construct the name of the state*/
string constructFileName(int N, int D,double x, double mg, int dlink,int level_nr);
/**Write to logfile such that I can already monitor during the run what is going on*/
void write_logfile(string message);
/**Function to determine the next state and what do with it*/
int getNextState(int N, int D, int Dprev, double x, double mg, int dlink, int dlinkprev ,int level_nr,string inputdir, string &filename);
/**Construct an effective Hamiltonian for the subspace of states I already have and diagonalize it to determine the coefficients for the states*/
int getCoefficients(vector<MPS*> states,vector<complex_t> &coefficients, const MPO &hamiltonian, Contractor & contractor);

/** Struct defining a state consisting of energy, level number and the file name it obtains*/
struct state
{
  int number;
  double energy;
  double phase;
  double momentum;
  string name;
  bool operator<(const state& ex) const{ return (energy<ex.energy);}  
  bool operator<=(const state& ex) const{ return (energy<=ex.energy);} 
  bool operator>=(const state& ex) const{ return (energy>=ex.energy);} 
  bool operator>(const state& ex) const{ return (energy>ex.energy);} 
};

/** Given a bunch of computed states, decide wether the scalar has already been found*/
bool foundScalar(vector<state> states);

int main(int argc,const char* argv[])
{
  //Code Revision (magic field, automatically filled in by subversion)
  string revision="$Rev: 1483 $";
  //Parse revision string
  revision=revision.substr(6,4);
  
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1)
  {
    cout << "Program compiled on " << __DATE__ << ", " << __TIME__ << " with code revision " << revision << endl;
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage:   " << argv[0] << " <inputfile>" << endl;
    return -1;
  }
  
  /**************************************************************************************
				 * Set up parameters*
   * ***********************************************************************************/
  
  const char* infile=argv[1];
  Properties props(infile);
  
  //System size
  int N = props.getIntProperty("N");
  //Bond dimension
  int D=props.getIntProperty("D");
  //Hamiltonian parameters
  double x=props.getDoubleProperty("x");
  double mg=props.getDoubleProperty("mg");  
  int dlink=props.getIntProperty("dlink");
  //Penalty for unphysical states
  double penalty=props.getDoubleProperty("penalty");  
  //Penalty for states with non-vanishing flux at the site N
  double szpenalty=props.getDoubleProperty("szpenalty");  
  //Penalty for states with non-vanishing U(1) charge
  double qpenalty=props.getDoubleProperty("qpenalty");  
  //Minimum number of excited states to compute, if automode is not active this is the number of excited states which are computed
  int nmin = props.getIntProperty("nmin");
  //Contractor tolerance
  double tol=props.getDoubleProperty("tol");
  //Whether or not the results should be appended to a (possibly) existing file
  int app=props.getIntProperty("append");
  //Whether I run the automatic mode which tries to detect the scalar itself or not
  bool automode=(props.getIntProperty("automode")==1);
  //Input directory where I read the initial guess from
  string input_directory=props.getProperty("inputdir");
  //Results of previous bond dimension which should be reused
  int Dinit=props.getIntProperty("Dinit");
  //Results of previous link dimension which should be reused (optional)
  int dlinkinit=props.getIntProperty("dlinkinit");
  //Flag if the scalar state has been found (if I am not running in automode, I simply pretend that the state has been found to be able to stop after nmin states)
  bool scalar_found = !automode;
  
  //Make sure that input directory ends with the token "/"
  if(!input_directory.empty() && *input_directory.end()!='/')
    input_directory+='/';
  
  cout << "------------------------------------------------------------------------"<<endl;
  cout << "Parameters: " << endl
       << " --> N           = " << N << endl
       << " --> D           = " << D << endl
       << " --> x           = " << x << endl
       << " --> mg          = " << mg << endl
       << " --> dlink       = " << dlink << endl
       << " --> pen         = " << penalty << endl
       << " --> szpen       = " << szpenalty << endl
       << " --> qpen        = " << qpenalty << endl
       << " --> nmin        = " << nmin << endl
       << " --> tol         = " << tol << endl
       << " --> app         = " << app << endl
       << " --> automode    = " << automode << endl
       << " --> input_dir   = " << input_directory << endl;      
  cout << "------------------------------------------------------------------------"<<endl;
  
  /**************************************************************************************
			      * Set up MPOs & Contractor*
   * ***********************************************************************************/
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  contractor.setConvTol(tol);

  //Get the Hamiltonian (dimensionless version with charge penalty)
  SU2HamiltonianNoLinksDimlessQ HNGI(N,x,mg,dlink,penalty,szpenalty,qpenalty);
  const MPO& hamil=HNGI.getHMPO();        
  cout << "Successfully retrieved SU2Hamiltonian" << endl;
  
  //Now prepare a MPS for the ground state
  vector<int> dims(N,4);
  MPS gs(N,D,dims);  
  
  //Build the MPOs for the observables
  vector<MPO*> FluxOps,NocOps;
  FluxOps.resize(N-1);
  NocOps.resize(N);
  MPO* tmp_mpo;
  #pragma omp parallel for
  for (int i=1;i<N; i++)
  {    
    FluxOps[i-1]=new MPO(1);
    HNGI.getJsquaredMPO(*FluxOps[i-1],i);
    
    NocOps[i-1]=new MPO(1);
    HNGI.getNocMPO(*NocOps[i-1],i);
  } 
  NocOps[N-1]=new MPO(1);
  HNGI.getNocMPO(*NocOps[N-1],N);
  
  
  MPO PenaltyMPO(1),CShift(1),CShiftsquare(1),Ptrafo(1),CConjugation(1),CConjugationP(1),Ntot(1),Projector(1);
  HNGI.getPenaltyMPO2(PenaltyMPO,true,1.0);
  HNGI.getNocMPO(Ntot);
  HNGI.getProjectorMPO(Projector);
    
  //Get the cyclic shift square operator
  HNGI.getCyclicShiftMPO(CShift);
  const MPO* oprsC[2]={&CShift,&CShift};
  CShiftsquare.initLength(N);
  MPO::join(2,oprsC,CShiftsquare);
  
  //Charge conjugation operator
  HNGI.getParticleTrafoMPO(Ptrafo);  
  const MPO* oprsCC[2]={&Ptrafo,&CShift};
  CConjugation.initLength(N);
  MPO::join(2,oprsCC,CConjugation);
  
  const MPO* oprsCCP[2]={&CConjugation,&Projector};
  CConjugationP.initLength(N);
  MPO::join(2,oprsCCP,CConjugationP);
  cout << "Setting up operators done " << endl;
  
  //Prepare the file for the results
  ofstream file; 
  if(file_exists(FILENAME) && app)
  {
    file.open(FILENAME, ios_base::app);
    file << setprecision(15) << "#Created with code revision " << revision << endl;
  }
  else
  {
    file.open(FILENAME);
    file << "#Created with code revision " << revision << endl;
    file << setprecision(16) << "#N\tD\tdlink\tx\tmg" << endl;
    file << "#E0\t<Penalty>\t<C^(2)>\t<CConj>\tphase(<CConj>)\t<Ntot>" << endl;
    file << "#E1\t<Penalty>\t<C^(2)>\t<CConj>\tphase(<CConj>)\t<Ntot>" << endl;
    file << "#..." << endl;
    file << setprecision(15) << N << "\t" << D << "\t" << dlink << "\t" << x << "\t" << mg << endl;
  }
  file.close();
  
  
  /**************************************************************************************
			      * Actual computation*
   * ***********************************************************************************/
  
  //Struct for the current level
  state curr_state;
  //Vectors to save the states
  vector <state> states;
  vector<MPS*> computed_states;
  //Pointer for the current state I am working on
  MPS* ex;  
  //Vector to save expectation values of operators
  vector<complex_t> result;
  //Variables for the energy
  double E0=0.0,E1=0.0;
  //Variables for observables
  complex_t penalty_val, cconjugation_val, cshift_val,particle_val;
  //Auxiliary variables for filenames and statuses
  stringstream sstm,tmpname;
  int status;
  string nextfile; 
  
  
  status=getNextState(N,D,Dinit,x,mg,dlink,dlinkinit,0,input_directory,nextfile);
  if(status==0)
  {
    //Nothing to compute
    gs.importMPS(nextfile.c_str());
    cout << "Found ground state ready to use in " << nextfile << endl;
    E0=real(contractor.contract(gs,hamil,gs));
    //cshift_val=contractor.contract(gs,CShiftsquare,gs);
    //cconjugation_val=contractor.contract(gs,CConjugation,gs);
    //Now rename the state to prevent being overwritten at the end while reordering
    status=rename(nextfile.c_str(),"GS.dat");
    if(status!=0)
	write_logfile("!Could not rename input state " + sstm.str() +", might be overwritten while reordering");
  }
  else 
  {
    if(status==1)
    {
      gs.importMPS(nextfile.c_str());
      cout << "Found initial guess for the ground state " << nextfile << endl;
    }
    else if(status==2)
    {
      vector <int> dummy;    
      HNGI.constructProductStateMPS(gs,dummy);
    }
    contractor.findGroundState (hamil,D,&E0,gs);
    gs.exportMPS("GS.dat");
    //Compute some observables
    //cshift_val=contractor.contract(gs,CShiftsquare,gs);
    //cconjugation_val=contractor.contract(gs,CConjugation,gs);    
    //cconjugation_val=contractor.contract(gs,CConjugationP,gs)/sqrt(contractor.contract2(Projector,gs));
    //cout << "Norm loss=" << ONE_c - sqrt(contractor.contract2(Projector,gs)) << endl;
    //penalty_val=contractor.contract(gs,PenaltyMPO,gs);
    //particle_val=contractor.contract(gs,Ntot,gs);
    
    file.open(FILENAME, ios_base::app);
    file << E0 << "\t"
         << penalty_val << "\t"
	 << cshift_val << "\t" 
	 << cconjugation_val << "\t"
	 << phase(cconjugation_val) << "\t"
	 << particle_val << endl;
    file.close();
    
    cout << "<mps|C^(2)|mps>=" << cshift_val << endl;
    cout << "<mps|CConjugation|mps>=" << cconjugation_val << ", phi="<< phase(cconjugation_val)  << endl;
    cout << "<mps|P|mps>=" << penalty_val << endl;
    cout << "<mps|Ntot|mps>=" << particle_val << endl<<endl;
    
    //Some error checking
    if(fabs(real(penalty_val))>100.0*tol)
      write_logfile("!Expectation value of penalty for ground state might is larger than tolerance");
    if(fabs(real(particle_val)-N)>100.0*tol)
      write_logfile("!Expectation value of total particle number might deviate than tolerance");
  }
  
  curr_state.number=0;
  curr_state.energy=E0;
  curr_state.phase=phase(cconjugation_val);
  curr_state.momentum=real(cshift_val);
  states.push_back(curr_state);
  computed_states.push_back(&gs);
  write_logfile("Finished computing ground state");
  
  
  for(int i=1; i<=MAXSTATES; i++)
  {
    //Still something to compute or can I break the loop
    if(i>nmin && scalar_found)
      break;
    
    ex=new MPS(N,D,4);
    
    sstm.str("");
    sstm << i << "EX.dat";
    tmpname.str("");
    tmpname << i << "EX.tmp";
    status=getNextState(N,D,Dinit,x,mg,dlink,dlinkinit,i,input_directory,nextfile);    
    if(status==0)
    {
      //Nothing to compute      
      ex->importMPS(nextfile.c_str());
      cout << "Found " << i <<". excited state ready to use in " << nextfile << endl;
      E1=real(contractor.contract(*ex,hamil,*ex));
      //cshift_val=contractor.contract(*ex,CShiftsquare,*ex);
      //cconjugation_val=contractor.contract(*ex,CConjugation,*ex);
      //Now rename the state to prevent being overwritten at the end while reordering
      status=rename(nextfile.c_str(),sstm.str().c_str());
      if(status!=0)
	write_logfile("!Could not rename input state " + sstm.str() +", might be overwritten while reordering");
    }
    else
    {
      if(status==1)	
      {
	ex->importMPS(nextfile.c_str());
	cout << "Found initial guess for " << i <<". excited state in " << nextfile << endl;	
      }
      else
	ex->setRandomState();
      
      contractor.findNextExcitedState(hamil,D,computed_states,&E1,*ex,0,1,tmpname.str().c_str(),MAXTIME);
      ex->exportMPS(sstm.str().c_str());
      //Computation finished, remove temporary file (if there was any)
      if(remove(tmpname.str().c_str())!=0)
	cout << "Warning could not remove temporary file" << endl;     
      //Compute some observables
      //cshift_val=contractor.contract(*ex,CShiftsquare,*ex);
      //cconjugation_val=contractor.contract(*ex,CConjugation,*ex);
      //cconjugation_val=contractor.contract(*ex,CConjugationP,*ex)/sqrt(contractor.contract2(Projector,*ex));
      //cout << "Norm loss=" << ONE_c - sqrt(contractor.contract2(Projector,*ex)) << endl;;
      //penalty_val=contractor.contract(*ex,PenaltyMPO,*ex);
      //particle_val=contractor.contract(*ex,Ntot,*ex);
      
      file.open(FILENAME, ios_base::app);
      file << E1 << "\t"
	  << penalty_val << "\t"
	  << cshift_val << "\t" 
	  << cconjugation_val << "\t"
	  << phase(cconjugation_val) << "\t"
	  << particle_val << endl;
      file.close();
      
      cout << "<mps|C^(2)|mps>=" << cshift_val << endl;
      cout << "<mps|CConjugation|mps>=" << cconjugation_val << ", phi="<< phase(cconjugation_val)  << endl;
      cout << "<mps|P|mps>=" << penalty_val << endl;
      cout << "<mps|Ntot|mps>=" << particle_val << endl<<endl;   
      
      
      //Some error checking
      if(fabs(real(penalty_val))>100.0*tol)
	write_logfile("!Expectation value of penalty for ground state might is larger than tolerance");
      if(fabs(real(particle_val)-N)>100.0*tol)
	write_logfile("!Expectation value of total particle number might deviate than tolerance");
    } 
    
    curr_state.number=i;
    curr_state.energy=E1;
    curr_state.phase=phase(cconjugation_val);
    curr_state.momentum=real(cshift_val);
    states.push_back(curr_state);
    computed_states.push_back(ex);
    sstm.str("");
    sstm << i;
    write_logfile("Finished computing " + sstm.str() +". excited state");
    
    if(automode && foundScalar(states))
    {
      //TODO Do something which allows to detect the scalar
      cout << "Detected candidate for the scalar state, will end computation" << endl;
      break;
    }
  }

  
  
  /**************************************************************************************
				 * Clean up and free memory *
   * ***********************************************************************************/
  file.close(); 
  
  //As the excited states my not come in order, I reorder the in the end according to their energy, therefore I sort them according to energy, generate the appropriate filenames and relabel the output states which are simply numbered in increasing order
  std::sort(states.begin(),states.end());
      
  for(int i=0; i<states.size(); i++)
  {    
    cout << "State " << i << " has energy " << states[i].energy << endl;
    sstm.str("");
    if(states[i].number==0)
      sstm << "GS.dat";
    else
      sstm << states[i].number << "EX.dat";
    states[i].name=constructFileName(N,D,x,mg,dlink,i);
    status=rename(sstm.str().c_str(),states[i].name.c_str());
    if(status!=0)
      write_logfile("!Could not rename level " + sstm.str());
    cout << states[i].number << ". state is now " << i << ". energy level" << endl;
  } 
  
  //Free the memory, do NOT delete the ground state!
  for(int i=1; i<computed_states.size(); i++)
  {
    computed_states[i]->clear();
    delete computed_states[i];
  }
    
  //Free memory
  for(int k=0; k<N-1; k++)
  {
    FluxOps[k]->clear();
    delete FluxOps[k];
    
    NocOps[k]->clear();
    delete NocOps[k];
  }
  
  NocOps[N-1]->clear();
  delete NocOps[N-1];
  
  cout << "End of program" << endl;
  return 0;
  
}
void print_MPS(MPS mps)
{
  cout << "MPS of length " << mps. getLength() << endl;
  for(int i=0; i<mps.getLength(); i++)
  {
    cout << "[" << i << "]=" <<  (mps.getA(i)).getA() << endl;
  }
}

vector<complex_t> get_full_config(vector<MPO*> NocOps, vector<MPO*> FluxOps, Contractor &contractor, MPS &mps, int N)
{
  vector <complex_t> res;
  res.resize(2*N-1);
  
#pragma omp parallel for
  for(int i=1; i<=N; i++)
  {
    res[2*i-2]=contractor.contract(mps,*NocOps[i-1],mps);
    if(i<N)
      res[2*i-1]=contractor.contract(mps,*FluxOps[i-1],mps);    
  }    
  return res;
}

vector<complex_t> get_config(vector<MPO*> Ops, Contractor &contractor, MPS &mps, int N)
{
  vector <complex_t> res;
  res.resize(Ops.size());
  
  #pragma omp parallel for
  for(int i=0; i<Ops.size(); i++)
    res[i]=contractor.contract(mps,*Ops[i],mps);
    
  return res;
}

//Once a configuration vector has been obtained with the previous function, this function saves it
int save_configuration(string filename, vector<complex_t> config)
{
  ofstream file;
  file.open(filename.c_str(),ios_base::app);
  for(int k=0; k<config.size()-1; k++)
    file << setprecision(16) << scientific << config[k] << "\t";
  file << config.back() << endl;
  file.close();
  return 0;
}

double sum_abs(vector<complex_t> input)
{
  double result=0.0;
  for(int i=0; i<input.size(); i++)
    result += + abs(input[i]);
  return result;
}

string constructFileName(int N, int D,double x, double mg, int dlink,int level_nr)
{
  stringstream sstm;    
  sstm.str("");
  if(level_nr==0)
    sstm << "GS_";
  else
    sstm<< level_nr << "EX_";
  sstm << "N" << N << "D" << D << "d" << dlink << "x" << (int) x << "mg" << (int) (mg*1000) << ".dat";
  return sstm.str();
}

int getNextState(int N, int D, int Dinit, double x, double mg, int dlink, int dlinkinit ,int level_nr,string inputdir, string &filename)
{
  int status;
  stringstream sstm; 
  
  /*Return status signals what to do:
   * 0: Nothing to compute, just read level
   * 1: Found valid initial guess, start computation with this level
   * 2: Nothing found, start with random state*/
  
   if(level_nr==0)
    sstm << "GS";
  else
    sstm<< level_nr << "EX";
  
  if(file_exists(constructFileName(N,D,x,mg,dlink,level_nr)))
  {
    //File has already been computed in a previously completed run, so simply read it
    filename=constructFileName(N,D,x,mg,dlink,level_nr);
    return 0;
  }
  else if(file_exists(sstm.str()+".dat"))
  {
    //The current run has not been completed yet, but I already have this state, just read it
    filename=sstm.str()+".dat";
    return 0;
  }
  else if(file_exists(sstm.str()+".tmp"))
  {
    //Temporary file from an incomplete sweep, use as initial guess
    filename=sstm.str()+".tmp";
    return 1;
  }
  else if(file_exists(inputdir+constructFileName(N,Dinit,x,mg,dlinkinit,level_nr)))
  {
    //Use initial guess from previous run with different bond dimension (and potentially different dlink)
    filename=inputdir+constructFileName(N,Dinit,x,mg,dlinkinit,level_nr);
    return 1;
  }
  else
  {
    //Nothing found, start from random state
    filename="";
    return 2;
  }
}

bool foundScalar(vector<state> states)
{
  bool scalar_found=false;
  //First reorder the states according to their energy (original vector is passed as copy and should stay unaffected)
  std::sort(states.begin(),states.end());
  //Then compute the momentum quantum number with the difference to the ground state
  vector<std::pair<double,int> > momentumdiffrence;
  for(int i=0; i<states.size(); i++)
    momentumdiffrence.push_back(make_pair(fabs(states[i].momentum - states[0].momentum),i));
  //Now sort the momentum differences
  std::sort(momentumdiffrence.begin(), momentumdiffrence.end());
  //Read through them and find which one is the scalar candidate
  for(int i=1; i<momentumdiffrence.size(); i++)
  {
    cout << "Delta P=" << momentumdiffrence[i].first << ", index=" << momentumdiffrence[i].second << endl;
    cout << "Phase=" << states[momentumdiffrence[i].second].phase << endl;
    if(fabs(states[momentumdiffrence[i].second].phase)<M_PIl/2.0)
    {
      scalar_found=true;
      break;
    }
    else 
      cout << "Does not qualify as scalar as phase is not correct" << endl;
  }
  
  return scalar_found;
}

void write_logfile(string message)
{  
  //Get date and time
  time_t t = time(0);   
  struct tm * now = localtime(&t);
  
  ofstream file;
  file.open(LOGFILE, ios_base::app);
  if(!file.is_open())
  {
    cout << "Error could not write to logfile on " <<  asctime(now) << endl;
    return;
  }
  file << "--> Entry from " <<  asctime(now);
  file << message << endl << endl;
  file.close();
}

int getCoefficients(vector<MPS*> states,vector<complex_t> &coefficients, const MPO &hamiltonian, Contractor & contractor)
{
  //Construct effective Hamiltonian
  int num_basis_states=states.size();
  complex_t entry;
  mwArray Heff(Indices(num_basis_states,num_basis_states)),U;
  Heff.fillWithZero();
  
  //Hamiltonian is Hermitan, so only compute lower part
  for(int i=0; i<num_basis_states; i++)
    for(int j=0; j<=i; j++)
    {
      entry=contractor.contract(*states[i],hamiltonian,*states[j]);
      if(i!=j)
      {
	Heff.setElement(entry,Indices(i,j));
	Heff.setElement(conjugate(entry),Indices(j,i));
      }
      else
	Heff.setElement(ONE_c*real(entry),Indices(i,j));
    }
    
  cout << "Heff=" << Heff << endl;
  wrapper::eig(Heff,coefficients,U);  
}
    