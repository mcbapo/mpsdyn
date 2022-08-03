/**
   \file SU2NoLinksSpectrum2.cpp
   
   This program is designed for large production runs to ultimately determine the spectrum of the SU(2) model in the formulation where the links are no longer explicitly present. The aim is to compute MPS approximations for the vector and the scalar and to determine in an automated way if they have already been found such that I do not have to take care of that myself. Furthermore options can be overridden to manually take control over the number of excited states which should be computed. 
   
   Additionally there is quite some error checking going on to make automated productions runs safe and easy to restart.   
   
  \param <inputfile> (string) File containing the parameters for the run
  
  \author Stefan Kühn
  \date 23/11/2015

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
#define FILENAME_TMP "SU2nolinks.tmp"
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
string constructFileName(int N, int D,double x, double mg, int dlink, double tol, int level_nr);
/**Write to logfile such that I can already monitor during the run what is going on*/
void write_logfile(string message);
/**Function to determine the next state and what do with it*/
int getNextState(int N, int D, int Dprev, double x, double mg, int dlink, int dlinkprev ,int level_nr, double tol, double tolinit, string inputdir, string &filename);
/**Construct an effective Hamiltonian for the subspace of states I already have and diagonalize it to determine the coefficients for the states*/
void getCoefficients(vector<MPS*> states,vector<complex_t> &coefficients, const MPO &hamiltonian, Contractor & contractor);
//Prepare initial state for refining
void getInitStateRefining(vector <MPS*> computed_states,vector <complex_t> coefficients,MPS &init_state_refine,Contractor & contractor);
//Refine computed states with the help of an effective Hamiltonian
void refineStates(vector <MPS*> &computed_states, const MPO &hamiltonian, Contractor & contractor,int D,int maxSweeps);

/** Struct defining a state consisting of energy, level number and the file name it obtains*/
struct state
{
  int number;
  double energy;
  double phase;
  double momentum;
  double penalty;
  double Ntot;
  complex_t cconjugation;
  string name;
  bool operator<(const state& ex) const{ return (energy<ex.energy);}  
  bool operator<=(const state& ex) const{ return (energy<=ex.energy);} 
  bool operator>=(const state& ex) const{ return (energy>=ex.energy);} 
  bool operator>(const state& ex) const{ return (energy>ex.energy);} 
};

/** Given a bunch of computed states, decide whether the vector/scalar has already been found*/
bool foundScalar(vector<state> states);
bool foundVector(vector<state> states);
bool finishedComputation(vector<state> states, int nmin, int nmax, string mode);

int main(int argc,const char* argv[])
{
  //Code Revision (magic field, automatically filled in by subversion)
  string revision="$Rev: 1496 $";
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
  //Minimum number of excited states to compute
  int nmin = props.getIntProperty("nmin");
  //Maximum number of excited states to compute
  int nmax = props.getIntProperty("nmax");
  //Contractor tolerance
  double tol=props.getDoubleProperty("tol");
  //Whether or not the results should be appended to a (possibly) existing file
  int app=props.getIntProperty("append");
  //Whether I run the automatic mode which tries to detect the scalar itself or not
  string automode=(props.getProperty("automode"));
  std::transform(automode.begin(), automode.end(), automode.begin(), ::tolower);
  if(automode.compare("1")==0)
  {
    //Ensure compatibility with older program version
    automode="scal";
  }
  else if(automode.compare("vec")!=0 && automode.compare("scal")!=0 && automode.compare("0")!=0)
  {
    //Neither scalar or vector given, invalid input, just go up to the scalar
    cout << "Warning, could not detect type of state for automode, going for the scalar" << endl;
    automode="scal";
  }
  //Input directory where I read the initial guess from
  string input_directory=props.getProperty("inputdir");
  //Results of previous bond dimension which should be reused
  int Dinit=props.getIntProperty("Dinit");
  //Results of previous link dimension which should be reused
  int dlinkinit=props.getIntProperty("dlinkinit");
  //Results of previous tolerance which should be reused
  double tolinit=props.getDoubleProperty("tolinit");
  //Maximum amount of sweeps allowed for the computation of the excited states
  int maxSweeps=props.getIntProperty("maxSweeps");
  //Should the results be refined if they did no converge properly within the given amount of sweeps
  int refine=props.getIntProperty("refine");
  
  //Make sure that input directory ends with the token "/", therefore trim at the end and add "/" if needed
  input_directory.erase(input_directory.find_last_not_of(" \n\r\t")+1);
  if(!input_directory.empty() && *input_directory.rbegin()!='/')
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
       << " --> input_dir   = " << input_directory << endl
       << " --> Dinit       = " << Dinit << endl
       << " --> dlinkinit   = " << dlinkinit << endl
       << " --> tolinit     = " << tolinit << endl
       << " --> maxSweeps   = " << maxSweeps << endl
       << " --> refine      = " << refine << endl;
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
  /*vector<MPO*> FluxOps,NocOps;
  FluxOps.resize(N-1);
  NocOps.resize(N);
  #pragma omp parallel for
  for (int i=1;i<N; i++)
  {    
    FluxOps[i-1]=new MPO(1);
    HNGI.getJsquaredMPO(*FluxOps[i-1],i);
    
    NocOps[i-1]=new MPO(1);
    HNGI.getNocMPO(*NocOps[i-1],i);
  } 
  NocOps[N-1]=new MPO(1);
  HNGI.getNocMPO(*NocOps[N-1],N);*/
  
  
  MPO PenaltyMPO(1),CShift(1),CShiftsquare(1),CConjugation(1),CConjugationP(1),Ntot(1),Shift(1);
  //HNGI.getShiftMPO(Shift,true);
  HNGI.getPenaltyMPO2(PenaltyMPO,true,1.0);
  HNGI.getNocMPO(Ntot);
  //HNGI.getChargeConjugationMPO(CConjugationP,true);
  HNGI.getChargeConjugationMPO(CConjugation,false);
    
  //Get the cyclic shift square operator
  HNGI.getCyclicShiftMPO(CShift);
  const MPO* oprsC[2]={&CShift,&CShift};
  CShiftsquare.initLength(N);
  MPO::join(2,oprsC,CShiftsquare);
  cout << "Setting up operators done " << endl;
  
  //Prepare the file for the results
  ofstream file; 
  if(file_exists(FILENAME_TMP) && app)
  {
    file.open(FILENAME_TMP, ios_base::app);
    file << setprecision(15) << "#Created with code revision " << revision << endl;
  }
  else
  {
    file.open(FILENAME_TMP);
    file << "#Created with code revision " << revision << endl;
    file << setprecision(16) << "#N\tD\tdlink\tx\tmg\ttol\tpen\tszpen\tqpen" << endl;
    file << "#E0\t<Penalty>\t<C^(2)>\t<CConj>\tphase(<CConj>)\t<Ntot>" << endl;
    file << "#E1\t<Penalty>\t<C^(2)>\t<CConj>\tphase(<CConj>)\t<Ntot>" << endl;
    file << "#..." << endl;
    file << setprecision(15) << N << "\t" << D << "\t" << dlink << "\t" << x << "\t" << mg << "\t" << tol << "\t"  << penalty << "\t" << szpenalty << "\t" << qpenalty << endl;
  }
  file.close();
  
  
  /**************************************************************************************
			      * Actual computation*
   * ***********************************************************************************/
  
  //Struct for the current level
  state curr_state;
  //Vectors to save the states one for the pointers to the MPSs and one for the observables
  vector <state> states;
  vector<MPS*> computed_states;
  //Pointer for the current state I am working on
  MPS* ex;  
  //Vector to save expectation values of operators
  vector<complex_t> result;
  //Variable for the energy
  double E1=0.0;
  //Variables for observables
  complex_t penalty_val, cconjugation_val, cshift_val,particle_val, shift_val;
  //Auxiliary variables for filenames and statuses
  stringstream sstm,tmpname,curr_level_nr;
  int status;
  string nextfile; 
  
  MPS init_state_refine;
  
  int level_nr=0;
  bool refining_needed=false;
  bool refining_done=false;
  while(1)
  {
    //Still something to compute or can I break the loop/start refining my results
    if(finishedComputation(states,nmin,nmax,automode))
      break;

    //Allocate new state and determine preliminary name                                                                                                        
    ex=new MPS(N,D,4);
    sstm.str("");
    tmpname.str("");
    curr_level_nr.str("");
    curr_level_nr << level_nr;
    if(level_nr>0)
    {
      sstm << level_nr << "EX.dat";      
      tmpname << level_nr << "EX.tmp";
    }
    else
    {
      sstm << "GS.dat";
      tmpname << "GS.tmp";
    }
    status=getNextState(N,D,Dinit,x,mg,dlink,dlinkinit,level_nr,tol,tolinit,input_directory,nextfile);
    if(status==0)
    {
      //Nothing to compute, just import appropriate state      
      ex->importMPS(nextfile.c_str());
      cout << "Found " << level_nr <<". state ready to use in " << nextfile << endl;
      write_logfile("Found "+curr_level_nr.str()+". state ready to use in "+nextfile);
      E1=real(contractor.contract(*ex,hamil,*ex));
      cshift_val=contractor.contract(*ex,CShiftsquare,*ex);
      cconjugation_val=contractor.contract(*ex,CConjugation,*ex);
      /*cconjugation_val=contractor.contract(*ex,CConjugationP,*ex)/sqrt(contractor.contract2(CConjugationP,*ex));
      //cout << "Norm loss=" << ONE_c - sqrt(contractor.contract2(CConjugationP,*ex)) << endl;*/
      penalty_val=contractor.contract(*ex,PenaltyMPO,*ex);
      particle_val=contractor.contract(*ex,Ntot,*ex);
      //shift_val=contractor.contract(*ex,Shift,*ex);
      //Now rename the state to prevent being overwritten at the end while reordering
      status=rename(nextfile.c_str(),sstm.str().c_str());
      if(status!=0)
	write_logfile("!Could not rename input state " + sstm.str() +", might be overwritten while reordering");
    }
    else
    {
      if(status==1)	
      {
	//Initial guess for input exists, use it for further computation
	ex->importMPS(nextfile.c_str());
	cout << "Found initial guess for " << level_nr << ". state in " << nextfile << " with energy " << contractor.contract(*ex,hamil,*ex) << endl;
	write_logfile("Found initial guess for " + curr_level_nr.str() + ". state in "+nextfile);
      }
      else if(status==2)
      {
	//No initial guess, use a strong coupling GS for ground state search or a random one for excited states
	if(level_nr==0)
	{
	  vector<int> dummy;
	  HNGI.constructProductStateMPS(*ex,dummy);
	}
	else
	    ex->setRandomState();
      }
      //Now compute the state
      if(level_nr==0)
      {
	contractor.findGroundState (hamil,D,&E1,*ex);
	status=0;
      }
      else
	status = contractor.findNextExcitedState(hamil,D,computed_states,&E1,*ex,0,1,tmpname.str().c_str(),MAXTIME,maxSweeps);
      
      //Check if there are problems with convergence, if there are some, set the flag for refining
      if(status==-1)
      {
	refining_needed=true;
	write_logfile("#Could not converge results for " + curr_level_nr.str() + ". state with maxSweeps during regular run, going for refinement run if refining flag is set");
      }
      
      //Computation finished, remove temporary file (if there is any) and save the final state
      ex->exportMPS(sstm.str().c_str());     
      if(remove(tmpname.str().c_str())!=0)
	cout << "Warning could not remove temporary file" << endl;     
      //Compute some observables
      cshift_val=contractor.contract(*ex,CShiftsquare,*ex);
      cconjugation_val=contractor.contract(*ex,CConjugation,*ex);
      /*cconjugation_val=contractor.contract(*ex,CConjugationP,*ex)/sqrt(contractor.contract2(CConjugationP,*ex));
      cout << "Norm loss=" << ONE_c - sqrt(contractor.contract2(CConjugationP,*ex)) << endl;*/
      penalty_val=contractor.contract(*ex,PenaltyMPO,*ex);
      particle_val=contractor.contract(*ex,Ntot,*ex);
      //shift_val=contractor.contract(*ex,Shift,*ex);
      
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
    //Write results to temporary file to allow for itnermediate evalutions
    file.open(FILENAME_TMP, ios_base::app);
    file << E1 << "\t"
	<< penalty_val << "\t"
	<< cshift_val << "\t" 
	<< cconjugation_val << "\t"
	<< phase(cconjugation_val) << "\t"
	<< particle_val << endl;
    file.close();
    
    curr_state.number=level_nr;
    curr_state.energy=E1;
    curr_state.phase=phase(cconjugation_val);
    curr_state.cconjugation=cconjugation_val;
    curr_state.momentum=real(cshift_val);
    curr_state.penalty=real(penalty_val);
    curr_state.Ntot=real(particle_val);
    states.push_back(curr_state);
    computed_states.push_back(ex);
    sstm.str("");
    sstm << level_nr;
    write_logfile("Finished computing " + sstm.str() +". excited state");
    
    level_nr++;
  }
  
  if(refining_needed && refine==1)
  {
    write_logfile("Going for refinement run");
    refineStates(computed_states,hamil,contractor,D,-1);
    //Update observables
    for(int i=1; i<computed_states.size(); i++)
    {
      states[i].energy=real(contractor.contract(*computed_states[i],hamil,*computed_states[i]));
      cconjugation_val=contractor.contract(*computed_states[i],CConjugationP,*computed_states[i])/sqrt(contractor.contract2(CConjugationP,*computed_states[i]));
      states[i].phase=phase(cconjugation_val);
      states[i].cconjugation=cconjugation_val;
      states[i].momentum=real(contractor.contract(*computed_states[i],CShiftsquare,*computed_states[i]));
      states[i].penalty=real(contractor.contract(*computed_states[i],PenaltyMPO,*computed_states[i]));
      states[i].Ntot=real(contractor.contract(*computed_states[i],Ntot,*computed_states[i]));
    }      
  } 
  
  /**************************************************************************************
				 * Clean up and free memory *
   * ***********************************************************************************/
  file.close(); 
  
  //As the excited states my not come in order, I reorder the in the end according to their energy, therefore I sort them according to energy, generate the appropriate filenames and relabel the output states which are simply numbered in increasing order
  std::sort(states.begin(),states.end());
  
  //Write to final file in order such that it is easier to read and evaluate and rename the states properly
  if(file_exists(FILENAME) && app)
  {
    file.open(FILENAME, ios_base::app);
    file << setprecision(15) << "#Created with code revision " << revision << endl;
    file << setprecision(15) << N << "\t" << D << "\t" << dlink << "\t" << x << "\t" << mg << "\t" << tol <<  endl;
  }
  else
  {
    file.open(FILENAME);
    file << "#Created with code revision " << revision << endl;
    file << setprecision(16) << "#N\tD\tdlink\tx\tmg\ttol\tpen\tszpen\tqpen" << endl;
    file << "#E0\t<Penalty>\t<C^(2)>\t<CConj>\tphase(<CConj>)\t<Ntot>" << endl;
    file << "#E1\t<Penalty>\t<C^(2)>\t<CConj>\tphase(<CConj>)\t<Ntot>" << endl;
    file << "#..." << endl;
    file << setprecision(15) << N << "\t" << D << "\t" << dlink << "\t" << x << "\t" << mg << "\t" << tol << "\t" << penalty << "\t" << szpenalty << "\t" << qpenalty << endl;
  }
      
  for(int i=0; i<states.size(); i++)
  {    
    cout << "State " << i << " has energy " << states[i].energy << endl;
    sstm.str("");
    if(states[i].number==0)
      sstm << "GS.dat";
    else
      sstm << states[i].number << "EX.dat";
    states[i].name=constructFileName(N,D,x,mg,dlink,tol,i);
    status=rename(sstm.str().c_str(),states[i].name.c_str());
    if(status!=0)
      write_logfile("!Could not rename level " + sstm.str());
    cout << states[i].number << ". state is now " << i << ". energy level" << endl;
    
    file  << states[i].energy << "\t"
	  << states[i].penalty << "\t"
	  << states[i].momentum << "\t" 
	  << states[i].cconjugation << "\t" 
	  << states[i].phase << "\t" 
	  << states[i].Ntot << endl;
  } 
  file.close();
  
  //Free the memory
  for(int i=0; i<computed_states.size(); i++)
  {
    computed_states[i]->clear();
    delete computed_states[i];
  }
    
  //Free memory
  /*for(int k=0; k<N-1; k++)
  {
    FluxOps[k]->clear();
    delete FluxOps[k];
    
    NocOps[k]->clear();
    delete NocOps[k];
  }
  
  NocOps[N-1]->clear();
  delete NocOps[N-1];*/
  
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

string constructFileName(int N, int D,double x, double mg, int dlink, double tol, int level_nr)
{
  stringstream sstm;    
  sstm.str("");
  if(level_nr==0)
    sstm << "GS_";
  else
    sstm<< level_nr << "EX_";
  sstm << "N" << N << "D" << D << "d" << dlink << "x" << (int) x << "mg" << (int) (mg*1000) << "tol" << tol << ".dat";
  return sstm.str();
}

int getNextState(int N, int D, int Dinit, double x, double mg, int dlink, int dlinkinit ,int level_nr, double tol, double tolinit, string inputdir, string &filename)
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
  
  if(file_exists(constructFileName(N,D,x,mg,dlink,tol,level_nr)))
  {
    //File has already been computed in a previously completed run, so simply read it
    filename=constructFileName(N,D,x,mg,dlink,tol,level_nr);
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
  else if(file_exists(constructFileName(N,D,x,mg,dlink,tolinit,level_nr)))
  {
    //File has already been computed in a previously completed run but with a different tolerance level, so I use it as initial guess
    filename=constructFileName(N,D,x,mg,dlink,tolinit,level_nr);
    return 1;
  }
  else if(file_exists(inputdir+constructFileName(N,Dinit,x,mg,dlinkinit,tolinit,level_nr)))
  {
    //Use initial guess from previous run with different bond dimension (and potentially different dlink)
    filename=inputdir+constructFileName(N,Dinit,x,mg,dlinkinit,tolinit,level_nr);
    return 1;
  }
  else
  {
    //Nothing found, start from scratch
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
    //cout << "Delta P=" << momentumdiffrence[i].first << ", index=" << momentumdiffrence[i].second << endl;
    //cout << "Phase=" << states[momentumdiffrence[i].second].phase << endl;
    if(fabs(states[momentumdiffrence[i].second].phase)<M_PIl/2.0)
    {
      scalar_found=true;
      break;
    }
    /*else 
      cout << "State number " << i <<" does not qualify as scalar as phase is not correct" << endl;*/
  }
  return scalar_found;
}

bool foundVector(vector<state> states)
{
  bool vector_found=false;
  //First reorder the states according to their energy (original vector is passed as copy and should stay unaffected)
  std::sort(states.begin(),states.end());
  //Now look if there is a state having a phase of pi
  for(int i=0; i<states.size(); i++)
  {    
    if(fabs(states[i].phase)>M_PIl/2.0)
    {
      vector_found=true;
      if(i==0)
	cout << "Warning ground state seems to be the vector" << endl;
      break;
    }
    /*else 
      cout << "State number " << i <<" does not qualify as vector as phase is not correct" << endl;*/
  }
  return vector_found;
}

bool finishedComputation(vector<state> states, int nmin, int nmax, string mode)
{
  bool finished=false;
  //If I have already more than nmax states, I simply stop
  if(states.size()>=nmax)
  {
    cout << "Computed (more than) the maximum number of states, will stop computation" << endl;
    finished=true;
  }
  else if(mode.compare("vec")==0 && states.size()>=nmin && foundVector(states))
  {
    //Found the vector and I have more than nmin states
    cout << "Found vector state, will finish computation" << endl;
    finished=true;
  }
  else if(mode.compare("scal")==0 && states.size()>=nmin && foundScalar(states))
  {
    //Found the scalar and I have more than nmin states
    cout << "Found scalar state, will finish computation" << endl;
    finished=true;
  }
  cout << endl;
  return finished;    
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

void getCoefficients(vector<MPS*> states,vector<complex_t> &coefficients, const MPO &hamiltonian, Contractor & contractor)
{
  //Construct effective Hamiltonian
  int num_basis_states=states.size();
  complex_t entry;
  mwArray Heff(Indices(num_basis_states-1,num_basis_states-1)),U;
  Heff.fillWithZero();
  
  //Hamiltonian is Hermitian, so only compute lower part
  for(int i=1; i<num_basis_states; i++)
    for(int j=1; j<=i; j++)
    {
      entry=contractor.contract(*states[i],hamiltonian,*states[j]);
      if(i!=j)
      {
	Heff.setElement(entry,Indices(i-1,j-1));
	Heff.setElement(conjugate(entry),Indices(j-1,i-1));
      }
      else
	Heff.setElement(ONE_c*real(entry),Indices(i-1,j-1));
    }
    
  wrapper::eig(Heff,coefficients,U);  
}

void getInitStateRefining(vector <MPS*> computed_states,vector <complex_t> coefficients,MPS &init_state_refine,Contractor & contractor)
{
  //Copy pointers in another vector to be compliant with the interface
  vector<const MPS*> kets;
  for (int i=1; i<computed_states.size(); i++)
    kets.push_back(computed_states[i]);
  
  contractor.optimizeSum (kets,coefficients,init_state_refine);
  init_state_refine.gaugeCond('L',true);
}

void refineStates(vector <MPS*> &computed_states, const MPO &hamiltonian, Contractor & contractor,int D,int maxSweeps)
{
  //First get the coefficients from the effective Hamiltonian
  vector<complex_t> coefficients;
  getCoefficients(computed_states,coefficients,hamiltonian,contractor);
  //Construct initial state
  MPS init_state(*computed_states[0]);
  init_state.setRandomState();
  getInitStateRefining(computed_states, coefficients,init_state, contractor);  
  //Delete the old states
  int num_states=computed_states.size();
  for(int i=1; i<num_states; i++)
  {
    computed_states[i]->clear();
    delete computed_states[i];
  }
  computed_states.resize(1);
  cout << "Computes states: " << *computed_states[0] << endl;
  //Now recompute the new states
  MPS *ex;
  stringstream sstm, tmpname;
  int status;
  double E1;
  for(int i=1; i<num_states; i++)
  {
    sstm.str("");
    tmpname.str("");
    sstm << i << "EX_refine.dat";      
    tmpname << i << "EX.tmp";
    ex=new MPS(init_state);
    status =contractor.findNextExcitedState(hamiltonian,D,computed_states,&E1,*ex,0,1,tmpname.str().c_str(),MAXTIME,maxSweeps);
    computed_states.push_back(ex);
    ex->exportMPS(sstm.str().c_str());
  }
  //After successfully computing all states give them the usual name again
  for(int i=1; i<num_states; i++)
  {
    sstm.str("");
    tmpname.str("");
    tmpname << i << "EX_refine.dat";      
    sstm << i << "EX.dat";
    if(rename(tmpname.str().c_str(),sstm.str().c_str())!=0)
      write_logfile("!Could not rename state " + tmpname.str());
  }
}

    