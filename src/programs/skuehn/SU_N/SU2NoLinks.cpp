/**
   \file SU2NoLinks.cpp
   
   
  \param <N> (int) Number of Fermions
  \param <D> (int) Bond dimension
  \param <µ> (double) Dimensionless constant of the Hamiltonian
  \param <epsilon> (double) Hopping parameter
  \param <mu> (double) Mass
  \param <g> (double) Coupling constant
  
  \author Stefan Kühn
  \date 09/10/2015

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

#include "SU2HamiltonianNoLinks.h"
#include "SU2HamiltonianNoLinksQ.h"
#include "SU2HamiltonianNoLinksDimless.h"
#include "SU2HamiltonianNoLinksDimlessQ.h"

#define EPS 1E-8
#define NUM_OF_PARAMS 10
#define FILENAME "SU2nolinks.txt"
#define CONFIGFILE "SU2nolinksconfiguration.txt"
#define TMPFILE "tmp.dat"

//Maxtime (currently 3 days)
#define MAXTIME 259200

using namespace std;
using namespace shrt;

//bool file_exists(string filename);
void print_MPS(MPS);
//Get spin and flux configuration
vector<complex_t> get_full_config(vector<MPO*> NocOps, vector<MPO*> FluxOps, Contractor &contractor, MPS &mps, int N);
//Get charge configuration
vector<complex_t> get_config(vector<MPO*> Ops, Contractor &contractor, MPS &mps, int N);
//Save the current configuration to a given file
int save_configuration(string filename, vector<complex_t> config);
//Time Evolution
void compute_evolution(MPO &exp_even,MPO &exp_odd,MPO &exp_odd_half,Contractor &contractor, MPS &mps, int steps,bool normalize);
void compute_evolution(MPO &Umpo,Contractor &contractor, MPS &mps, int steps,bool normalize);
void print_config(vector<complex_t> confi);
//Compute the sum of the absolute values of the vector entries
double sum_abs(vector<complex_t> input);
//Given the projectors, project back into the Gauss Law fulfilling subspace
void project(MPO &Podd, MPO &Peven, Contractor &contractor, MPS &mps, bool normalize=false);
//Given the projector, project back into the Gauss Law fulfilling subspace
void project(MPO &Projector, Contractor &contractor, MPS &mps, bool normalize=false);

int main(int argc,const char* argv[]){
  int cntr=0;
  
  //Case that number of input arguments is to small, display some instruction
  if(argc < NUM_OF_PARAMS + 1)
  {
    cout << "Program compiled on " << __DATE__ << ", " << __TIME__ << endl;
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage:   " << argv[0] << " N D µ g epsilon" << endl;
    cout << "N:       " << "Number of spins" << endl
	 << "D:       " << "Bond dimension of the MPS" << endl
	 << "epsilon: " << "Hopping parameter" << endl
	 << "µ:       " << "Mass" << endl
	 << "g:       " << "Coupling constant" << endl
	 << "n:       " << "Number of states on a link" << endl
	 << "lambda:  " << "Constant for the penalty" << endl
	 << "sz_pen:  " << "Strength of the penalty targeting the subspace of total charge zero" << endl	 
	 << "n_ex:    " << "Number of excited states to be computed, n_ex<=0 means none" << endl
	 << "qpen:    " << "Penalty for the balance of particles and anti particles" << endl
	 << "mode:    " << "If set to something different from zero, before computing a state, it is checked wether it is available from a file and used without further calculation" << endl
	 << "auto:    " << "If auto is set to one, the program automatically tries to compute all states up to the scalar (still experimental)" << endl;
    return -1;
  }
  
  //Number of spin sites
  int N=atoi(argv[++cntr]);
  //Bond dimension
  int D=atoi(argv[++cntr]);
  //Hopping parameter
  double epsilon=atof(argv[++cntr]);
  //Mass
  double mu=atof(argv[++cntr]);
  //Coupling constant
  double g=atof(argv[++cntr]);
  //Number of excited states to be computed
  int dlink=atoi(argv[++cntr]);
  //Number of time steps to be computed
  double lambda=atof(argv[++cntr]);
  //Number of time steps to be computed
  double sz_pen=atof(argv[++cntr]);
  //Number of excited states to be computed
  int n_ex=atoi(argv[++cntr]); 
  //Penalty strength for the balance of particles and antiparticles
  double qpen=atof(argv[++cntr]); 
  //Mode, wether states should be used as inital guess for a computation, or simply replace the computation
  int mode=0;
  //Whether I want to automatically try to compute all states up to the scalar
  bool automode=false;
  
  if(argc>NUM_OF_PARAMS+1)
    mode = atoi(argv[++cntr]);
  if(argc>NUM_OF_PARAMS+2)
    if(atoi(argv[++cntr])==1)
    {
      automode = true;
      //Just set n_ex reasonably high, such that I find the scalar, if I haven't found it until then, I might probably not find it anyway
      n_ex= 50;
      cout << "Switching to automode" << endl;
    }
  //Revision
  string revision="$Rev: 1483 $";
  //Parse revision string
  revision=revision.substr(6,4);
  
  cout << "------------------------------------------------------------------------"<<endl;
  cout << "Parameters: " << endl
       << " --> N        = " << N << endl
       << " --> D        = " << D << endl
       << " --> epsilon  = " << epsilon << endl
       << " --> mu       = " << mu << endl
       << " --> g        = " << g << endl
       << " --> dlink    = " << dlink << endl
       << " --> lambda   = " << lambda << endl
       << " --> sz_pen   = " << sz_pen << endl
       << " --> n_ex     = " << n_ex << endl
       << " --> qpen     = " << qpen << endl;       
  cout << "------------------------------------------------------------------------"<<endl;
       
  
  //File for results
  ofstream file; 
  //Vector to save expectation values of operators
  vector<complex_t> result;
  //Variables for the energy
  double E0=0.0,E1=0.0;
  //Variable for the penalty term
  complex_t pen_val;
  //Variable for charge conjugation
  complex_t cconjugation;
  //Values for the distance in the cyclic_shift by two sites ("pseudo momentum")
  complex_t cgs=ZERO_c,c2=ZERO_c;
  double cdiff=0.0;
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif
  contractor.setConvTol(EPS);

  //Get the Hamiltonian
  //SU2HamiltonianNoLinks HNGI(N,epsilon,mu,g,dlink,lambda,sz_pen);
  //Get the Hamiltonian (version with charge penalty)
  //SU2HamiltonianNoLinksQ HNGI(N,epsilon,mu,g,dlink,lambda,sz_pen,qpen);
  //Get the Hamiltonian (dimensionless version)
  //SU2HamiltonianNoLinksDimless HNGI(N,epsilon,mu,dlink,lambda,sz_pen);
  //Get the Hamiltonian (dimensionless version with charge penalty)
  SU2HamiltonianNoLinksDimlessQ HNGI(N,epsilon,mu,dlink,lambda,sz_pen,qpen);
  const MPO& hamil=HNGI.getHMPO();    
    
  cout << "Successfully retrieved SU2Hamiltonian" << endl;
  
  srand(time(NULL));
  
  //Now prepare and MPS
  vector<int> dims(N,4);
  MPS mps(N,D,dims);
  cout << "Building MPS done" << endl;
  
  
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
  
  
  MPO penalty(1),penalty2(1), penalty3(1), particlediff(1), Qsquaretot(1), Pmpo(1), Pmpo_tmp(1), particlebalancepenalty(1),CShift(1),CShiftsquare(1),Ptrafo(1),CConjugation(1);
  HNGI.getPenaltyMPO(penalty3,1.0);
  HNGI.getPenaltyMPO2(penalty,true,1.0);
  HNGI.getPenaltyMPO2(penalty2,false,1.0);
  HNGI.getParticleBalanceMPO(particlediff);
  HNGI.getChargeSquareMPO(Qsquaretot);
  HNGI.getParticleBalancePenaltyMPO(particlebalancepenalty);  
  
  //Get P² operator
  /*HNGI.getMomentumMPO3(Pmpo_tmp);
  const MPO* oprsP[2]={&Pmpo_tmp,&Pmpo_tmp};
  Pmpo.initLength(N);
  MPO::join(2,oprsP,Pmpo);*/  
    
  //Get the cyclic shift square operator
  HNGI.getCyclicShiftMPO(CShift);
  const MPO* oprsC[2]={&CShift,&CShift};
  CShiftsquare.initLength(N);
  MPO::join(2,oprsC,CShiftsquare);
  
  HNGI.getParticleTrafoMPO(Ptrafo);
  //Charge conjugation operator?
  const MPO* oprsCC[2]={&Ptrafo,&CShift};
  CConjugation.initLength(N);
  MPO::join(2,oprsCC,CConjugation);
  
  
#ifdef MATOUT
  //Pmpo.exportForMatlab("Psquare.m");
  CShiftsquare.exportForMatlab("CShiftsquare.m");
  CConjugation.exportForMatlab("CConjugation.m");
#endif
  
  cout << "Setting up operators done " << endl;
 
  vector <int> dummy;    
#ifdef CHECKDIMS
  MPO Shift(1);
  HNGI.getShiftMPO(Shift,true);
  cout << endl << "Warning, doing some mapping for debugging" << endl;
  //Testing initial states
  dummy.resize(N,3);
  for(int i=1; i<N; i+=2)
    dummy[i]=0;
  for(int i=2; i<N; i+=2)
    dummy[i]=0;
  cout << "Dummy: " << dummy << endl;
  MPS tmp(N,D*4,4);
  tmp.setRandomState();
  HNGI.constructProductStateMPS(mps,dummy);
  cout << "<mps|tmp>=" << contractor.contract(mps,tmp) << endl;
  cout << "<mps|Shift|tmp>=" << contractor.contract(mps,Shift,tmp) << endl;
  cout << "<mps|Shift|mps>=" << contractor.contract(mps,Shift,mps) << endl;
  contractor.optimize(Shift,mps,tmp);
  cout << "Configuration=(";
  for(int i=0; i<N; i++)
    cout << contractor.contract(tmp,*NocOps[i],tmp) << ",";
  cout << ")" << endl;
  /*MPO nocnew(1);
  MPS mappedmps(N,D,3);
  HNGI.constructProductStateMPS(mps,dummy);
  cout << "MPS=" << mps << endl;
  cout << "Mapping MPO=" << Pmpo << endl;
  //contractor.optimize(Pmpo,mps,mappedmps,500);
  mappedmps.approximate(Pmpo,mps,500);
  //mappedmps.setRandomState();
  cout << "... done" << endl;
  cout << "Mapped mps=" << mappedmps << endl;
  cout << "Noc = (";
  for(int i=1; i<=N; i++)
  {
    HNGI.getNocMPO_other_basis(nocnew,i);
    cout << contractor.contract(mappedmps,nocnew,mappedmps) << ", " ;
  }  
  cout << ")" << endl;*/
  //Test the "momentum MPO" 
  /*HNGI.constructProductStateMPS(mps,dummy);
  cout << "Momentum of dummy configuration " << contractor.contract(mps,Pmpo,mps) << endl << endl;*/
  
  //Test the cyclic shift operator
  //cout << "Cyclic shift: " << contractor.contract(mps,CShift,mps) << endl;
  
  //Check the penalties, therefore create a bunch of configurations which react differently
  for(int i=0; i<N; i++)
    dummy[i]=0;
  dummy[N-1]=2;
  cout << "Dummy: " << dummy << endl;
  HNGI.constructProductStateMPS(mps,dummy);
  cout << "Old penalty:              " << contractor.contract(mps,penalty,mps) << endl;
  cout << "New penalty:              " << contractor.contract(mps,penalty2,mps) << endl;
  cout << "New penalty no last site: " << contractor.contract(mps,penalty3,mps) << endl;
  
  for(int i=0; i<N; i++)
    dummy[i]=2;
  cout << "Dummy: " << dummy << endl;
  HNGI.constructProductStateMPS(mps,dummy);
  cout << "Old penalty: " << contractor.contract(mps,penalty,mps) << endl;
  cout << "New penalty: " << contractor.contract(mps,penalty2,mps) << endl;
  cout << "New penalty no last site: " << contractor.contract(mps,penalty3,mps) << endl;
  
  return 0;
#endif
  
  //Check if there is an initial guess
  if(file_exists("GS.dat") && mode!=0)
  {
    //State exists and should be simply used
    mps.importMPS("GS.dat");
    cout << "Found ground state" << endl;
    E0 = real(contractor.contract(mps,hamil,mps));
    cout << "Energy = " << E0 << endl;
  }
  else
  {
    if(file_exists("GS.dat") && mode==0)
    {
      //State exists but should be used as inital guess
      mps.importMPS("GS.dat");
      cout << "Found inital guess for ground state with energy " << contractor.contract(mps,hamil,mps) << endl;
      cgs = contractor.contract(mps,CShiftsquare,mps);
    }
    else
      HNGI.constructProductStateMPS(mps,dummy);
    
    //Now compute the ground state and save it to disk
    contractor.findGroundState (hamil,D,&E0,mps);
    mps.exportMPS("GS.dat");    
  
    //Get the penalty expectation value
    pen_val = contractor.contract(mps,penalty,mps);
    //Get the full configuration of the ground state
    result=get_full_config(NocOps,FluxOps,contractor,mps,N);
    save_configuration(CONFIGFILE,result);    
    //Get the cyclic shift values
    cgs = contractor.contract(mps,CShiftsquare,mps);
    //Get the charge conjugation value
    cconjugation=contractor.contract(mps,CConjugation,mps);
    
    cout << "Ground state:" << endl;
    cout  << "    - Energy =           " << E0 << endl
	  << "    - Penalty =          " << pen_val << endl
	  << "    - P^2 =              " << cgs << endl
	  << "    - phase(CConj) =     " <<  phase(cconjugation) << endl;
	  /*<< "    - Particle balance = " << contractor.contract(mps,particlebalance,mps) << endl
	  << "    - Particle diff =    " << contractor.contract(mps,particlediff,mps) <<endl
	  << "    - Q²_tot =           " << contractor.contract(mps,Qsquaretot,mps) << endl
	  << "    - P_balance =        " << contractor.contract(mps,particlebalancepenalty,mps) << endl
	  << "    - SR =               " << contractor.contract(mps,SR,mps) << endl
	  << "    - Pmpo =             " << contractor.contract(mps,Pmpo_tmp,mps) << endl
	  << "    - Pmpo^2 =           " << contractor.contract(mps,Pmpo,mps) << endl;*/
    
    
    if(file_exists(FILENAME))
    {
      file.open(FILENAME, ios_base::app);
      file << endl << "#Created with code revision " << revision;
      file << endl << setprecision(16) << N << "\t" << D << "\t" << dlink << "\t" << epsilon << "\t" << mu << "\t" << g << endl;
      file << E0 << "\t" << pen_val << "\t" << cgs << "\t" << cconjugation << "\t" << phase(cconjugation);
    }
    else
    {
      file.open(FILENAME, ios_base::app);
      file << "#Created with code revision " << revision << endl;
      file << setprecision(16) << "#N\tD\tdlink\tepsilon\tmu\tg" << endl;
      file << "#E0\t<Penalty>\t<C^(2)>\t<CConj>\tphase(<CConj>)" << endl;
      file << "#E1\t<Penalty>\t<C^(2)>\t<CConj>\tphase(<CConj>)" << endl;
      file << "#..." << endl;
      file << setprecision(16) << N << "\t" << D << "\t" << dlink << "\t" << epsilon << "\t" << mu << "\t" << g << endl;
      file << E0 << "\t" << pen_val << "\t" << cgs << "\t" << cconjugation << "\t" << phase(cconjugation);
    }
    file.close();
  }
  
  
  //Try to get an excited state
  MPS first_ex;
  first_ex = mps;
  stringstream sstm;
    
  vector<MPS*> computed_states;
  computed_states.push_back(&mps);
  
  vector<double> energies;
  vector<complex_t> Gx_vals_ex,Gy_vals_ex,Gz_vals_ex;
  
  //Save also the ground state
  energies.push_back(E0);
  
  for(int i=1; i<=n_ex; i++)
  {
    sstm.str("");
    sstm << i << "EX.dat";
    //Check if there is an initial guess
    if(file_exists(sstm.str().c_str()) && mode!=0)
    {
      //State exists and should be simply used
      first_ex.importMPS(sstm.str().c_str());
      cout << "Found "<< i <<". excited state" << endl;
      computed_states.push_back(new MPS(first_ex));
      c2 = contractor.contract(first_ex,CShiftsquare,first_ex);
      cdiff = abs(c2-cgs);
    }
    else
    {
      if(file_exists(sstm.str().c_str()) && mode==0)
      {
	//State exists but should be used as initial guess
	first_ex.importMPS(sstm.str().c_str());
	cout << "Found initial guess for " << i << ". excited state with energy " << contractor.contract(first_ex,hamil,first_ex) << endl;
	if(first_ex.getBond()<D)
	{
	  first_ex.increaseBondDimension(D);
	  cout << " - Increased bond dimension, new energy " << contractor.contract(first_ex,hamil,first_ex) << endl;
	}
      }
      else
      {
	if(file_exists(TMPFILE))
	{
	  //Temporary result exists, reuse it
	  first_ex.importMPS(TMPFILE);
	  cout << "Found temporary result, going to reuse it" << endl;
	}
	else
	{
	  first_ex.setRandomState();
	  first_ex.gaugeCond('R',true);
	}
      }
      
      //Now compute the ground state and save it to disk
      contractor.findNextExcitedState(hamil,D,computed_states,&E1,first_ex,0,1,TMPFILE,MAXTIME);
      
      first_ex.exportMPS(sstm.str().c_str());
      computed_states.push_back(new MPS(first_ex));
      //Get the penalty expectation value      
      pen_val = contractor.contract(first_ex,penalty,first_ex);
      //Get the full configuration
      result=get_full_config(NocOps,FluxOps,contractor,first_ex,N);
      save_configuration(CONFIGFILE,result);
      //Get the cyclic shift values
      c2 = contractor.contract(first_ex,CShiftsquare,first_ex);
      //Get the charge conjugation value
      cconjugation=contractor.contract(first_ex,CConjugation,first_ex);
      cout << i <<". excited state:" << endl;
      cout << "    - Energy =           " << E1 << endl
	  << "    - Penalty =          " << pen_val << endl
	  << "    - phase(CConj) =     " << phase(cconjugation) << endl;
	  /*<< "    - Particle balance = " << contractor.contract(first_ex,particlebalance,first_ex) << endl
	  << "    - Particle diff =    " << contractor.contract(first_ex,particlediff,first_ex) <<endl
	  << "    - Q²_tot =           " << contractor.contract(first_ex,Qsquaretot,first_ex) << endl
	  << "    - P_balance =        " << contractor.contract(first_ex,particlebalancepenalty,first_ex) << endl
	  << "    - SR =               " << contractor.contract(first_ex,SR,first_ex) << endl
	  << "    - Pmpo =            " << contractor.contract(first_ex,Pmpo_tmp,first_ex) << endl
	  << "    - Pmpo^2 =          " << contractor.contract(first_ex,Pmpo,first_ex) << endl;*/
      file.open(FILENAME, ios_base::app);
      //file << endl << setprecision(16) << E1 << "\t" << pen_val << "\t";
      file << endl << setprecision(16) << E1 << "\t" << pen_val << "\t" << c2 << "\t" << cconjugation << "\t" << phase(cconjugation);
      file.close();   
      
      //Remove the temporary output to avoid using a temporary file for excited state n-1 as initial guess for excited state n, as the two states are orthogonal
      if(remove(TMPFILE)!=0)
	cout << "Warning could not remove temporary file" << endl;     
      
      //Check if I am about to cross the zero line, then I can stop because thing are not converging anymore
      if(abs(E1)<1E-1)
      {
	cout << "Energy is very close to zero, further results will not converge, aborting program" << endl;
	break;
      }
      
      //Check if I can end the search because I found the scalar
      if(abs(c2-cgs)<cdiff && real(c2)>0.0 && automode)
      {
	cout << "Difference between "<< i << ". Ex and GS: " << abs(c2-cgs) << endl;
	cout << "--> Found scalar state, will end program" << endl << endl;
	break;
      }
      else
      {
	cdiff = abs(c2-cgs);
	cout << "Difference between "<< i << ". Ex and GS: " << cdiff << endl;
      }
    }  
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


//Little function to compute the trotter splitting, as it is quite complicated and I want to take two time steps together to improve performance
void compute_evolution(MPO &exp_even,MPO &exp_odd,MPO &exp_odd_half,Contractor &contractor, MPS &mps, int steps,bool normalize)
{
   MPS aux=MPS(mps);
   
  //Now with a bit more optimization to save time  
  //Inital half step
  contractor.optimize(exp_odd_half,mps,aux);
  for(int i=1; i<=(steps-1); i++)
  {
    //cout << "\rStep " << i <<"      " ;
    cout.flush();   
    mps.setRandomState();
    contractor.optimize(exp_even,aux,mps);
    aux.setRandomState();
    contractor.optimize(exp_odd,mps,aux);
    //cout << endl << "Norm=" << contractor.contract(mps,mps) << endl;
    if(normalize)
    {
      //cout << "Normalizing" << endl;
      mps.gaugeCond('R',true);      
    }
  }
  
  //Complete last step
  contractor.optimize(exp_even,aux,mps);
  contractor.optimize(exp_odd_half,mps,aux);  
  mps=aux;
  if(normalize)
  {
    //cout << "Normalizing" << endl;
    mps.gaugeCond('R',true);      
  }  
  //cout << endl;
}

void compute_evolution(MPO &Umpo,Contractor &contractor, MPS &mps, int steps,bool normalize)
{
  MPS aux=MPS(mps);
  for(int i=1; i<=steps; i++)
  {
    contractor.optimize(Umpo,mps,aux);
    mps=aux;
    if(normalize)
      mps.gaugeCond('R',true);
  }
}

void project(MPO &Podd, MPO &Peven, Contractor &contractor, MPS &mps,bool normalize)
{
  MPS aux=MPS(mps);
  contractor.optimize(Podd,mps,aux);
  contractor.optimize(Peven,aux,mps);
  if(normalize)
    mps.gaugeCond('R',true);
}

void project(MPO &Projector, Contractor &contractor, MPS &mps,bool normalize)
{
  MPS aux=MPS(mps);
  contractor.optimize(Projector,mps,aux);
  mps = aux;
  if(normalize)
    mps.gaugeCond('R',true);
}


void print_config(vector<complex_t> confi)
{
  cout << "Config=(";
  int count=0;
  for(int i=0; i<confi.size(); i++)
  {
    count ++;
    if(count !=2 )
      cout << confi[i] << ",";
    else
    {
      cout << confi[i] << "   ";
      count = 0;
    }
  }
  cout << ")" << endl;
}

double sum_abs(vector<complex_t> input)
{
  double result=0.0;
  for(int i=0; i<input.size(); i++)
    result += + abs(input[i]);
  return result;
}
    