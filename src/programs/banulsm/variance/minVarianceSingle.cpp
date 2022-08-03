
#include <math.h>
#include <iomanip>
#include <cstring>
#include <strstream>

#include "misc.h"
#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results
#define MAXSTEPS 1E6 // If so many steps, I stop

#include "IsingHamiltonian.h"

using namespace shrt;

/** To find the GS of H^2, instead of running findGroundState,
    it applies repeatedly cos(x H/N), which approximates the 
    effect of imaginary time evolution of H^2.
    The Hamiltonian is Ising Hamiltonian \ref <IsingHamiltonian> 
    After running, compute the exp value of H^2 for all cuts

    Receives arguments:
    \param <propsfile> (string) mandatory: name of config props file where to read the rest. The rest can be set in the command line
    \param <L> (int) length of the chain
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D> (int) bond dimension
    \param <tol> (double) relative change in the exp value of H^2 to decide convergence (df 1E-5)
    \param <dt> (double) time step (for Trotter)
    \param <x> (double) paramenter multiplying the H in the argument of the cosine
    \param <outfname> (char*) name of the output file for the energy
    \param <rate> (int) frequency to record results on file
    \param <mpsdir> (char*) directory where to store the tmp (or final)
                     state, so that I can continue the evolution
    \param <saveFreq> (int) frequency to save MPS tmp results
    \param <mpsfile> (char*) base of the filename where to save or from 
                     which to read (only if newInstance is false) the
		     tmp MPS. The full name will be changed to include 
		     the bond dimension it corresponds to. The basis of
		     thename should be such that we do not mix
		     evolutions with different parameters (e.g. changing step width).
    \param <initmpsfile> (char*) if given explicitly, this file is used to read an initial guess. Otherwise, the program searches for some tmp file that itself wrote in a previous run.
    \param <app> (int) whether the output file is to be kept (app==1)
*/

void computeValues(const MPS& state,const MPO& hamil,vector<complex_t>& valH2,vector<complex_t>& valH);
int findTmpFile(int M,const string& mpsdir,const string& baseName,string& mpsfile);
const string mpsFileName(const string& mpsdir,const string& baseName,int cnt);

int lcm(int a,int b);

int d=2;

int main(int argc,const char* argv[]){
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

  int L = props.getIntProperty("L");  
  double J = props.getDoubleProperty("J");  
  double g = props.getDoubleProperty("g");  
  double h = props.getDoubleProperty("h");  
  int D = props.getIntProperty("D");  
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir");
  string mpsfile=props.getProperty("mpsfile");
  string initmpsfile=props.getProperty("initmpsfile");
  bool useInitMPS=initmpsfile.size()>0;
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  double deltaX=props.getDoubleProperty("x");
  double dt=props.getDoubleProperty("dt");
  if(dt<0) dt=deltaX; // default both steps the same
  int mS=(int)(deltaX/dt); // if not really an integer divisor, adjust
  if(abs(deltaX-mS*dt)>1E-6){
    mS=ceil(deltaX/dt);dt=deltaX/mS;
    cout<<"Adjusted the Trotter step to dt="<<dt<<", i.e. "<<mS
	<<" steps for x="<<deltaX<<endl;
  }
  int rate=props.getIntProperty("rate");
  if(rate<0) rate=1;
  int saveFreq=props.getIntProperty("saveFreq");
  if(saveFreq<0) saveFreq=MAXSTEPS;
  saveFreq=lcm(rate,saveFreq);
  double tolC=props.getDoubleProperty("tol");
    
  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    *out<<"% N\t J\t g\t h\t D\t dt\t deltaX\t cnt\t time \t<H2>\t <H>"<<endl;
    out->close();
    delete out;
  }
  cout<<setprecision(10);

  int cnt=0;
  Contractor& contractor=Contractor::theContractor();
  //  contractor.setEigenSolver(primme);
  cout<<"Initialized Contractor"<<endl;
  MPS gs(L,D,d);
  // Do we read from somefile?
  if(useInitMPS){
    cout<<"Will try to read initial state from "<<initmpsfile<<endl;
    gs.importMPS(initmpsfile.data());
  }
  else{
    // try to find a previoius tmp file with same root or init from scratch
    string tmpMPSfile;
    double time=0.;
    cnt=findTmpFile(MAXSTEPS,mpsdir,mpsfile,tmpMPSfile);
    if(cnt>0){
      cout<<"Recovered state for step "<<cnt<<" at file "<<tmpMPSfile<<endl;
      gs.importMPS(tmpMPSfile.data());
    }
    else{
      gs.setRandomState(); // the intial state, random
      cout<<"Initialized random state from scratch, norm "<<contractor.contract(gs,gs)<<endl;
    }
  }
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  
  // First: H
  IsingHamiltonian hamI(L,d,J,g,h); 
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamI.getHMPO(0.);
  cout<<"Constructed the hamil MPO"<<endl;

  cout<<"Initial state has <H>="<<contractor.contract(gs,hamil,gs)<<", <H^2>="<<contractor.contract2(hamil,gs)<<endl;

  // Run the evol step with positive and negative exponentials separately
  MPO expP(L),expM(L); // double layer positive and negative exponentials of the commutator
  // get the idividual MPOs from the Hamiltonian
  hamI.getUMPO(expP,dt*I_c,0); // the single layer positive exponential
  hamI.getUMPO(expM,-dt*I_c,0); // the negative exponential

 
  bool done=0;
  double time=cnt*deltaX;
  double expH2=L*100.;
  // place for other params
  //  complex_t E,EL,varL,ER,varR;
  complex_t valsE,valsE2;
  while(!done&&cnt<MAXSTEPS){ // find GS of H^2 until convergence of valE2 is sufficient
    // instead of findGS, I run cos(xH) and set a limit of up to MAXSTEPS steps 
    for(int k=0;k<rate;k++){
      // If mS>1, need to apply several operators of each sign before summing
      MPS aux1(gs),aux2(gs);
      for(int k=0;k<mS;k++){
	// could be optimized by saving copies, but depends on whether mS even or odd
	MPS _aux1(aux1),_aux2(aux2);
	contractor.optimize(expP,_aux1,aux1,D);
	contractor.optimize(expM,_aux2,aux2,D);
      }
      aux1.gaugeCond('R',1);aux2.gaugeCond('R',1);
      vector<const MPS*> vecs;vecs.push_back(&aux1);vecs.push_back(&aux2);
      vector<complex_t> coeffs(2,.5*ONE_c);
      //      cout<<"About to optimize the sum"<<endl;
      contractor.optimizeSum(vecs,coeffs,gs,D);
      cnt++;time+=deltaX;
      gs.gaugeCond('R',1);
    }
    valsE=contractor.contract(gs,hamil,gs);
    valsE2=contractor.contract2(hamil,gs);
    cout<<"State for step "<<cnt<<" found with <H^2>="<<valsE2<<" <H>="<<valsE<<endl;
    out=new ofstream(outfname.data(),ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<setprecision(15);
    *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"<<dt<<"\t"<<deltaX<<"\t";
    *out<<cnt<<"\t"<<time<<"\t"<<real(valsE2)<<"\t"<<real(valsE);
    *out<<endl;
    out->close(); delete out;
    // check for convergence
    if(abs(ONE_c-valsE2/expH2)<tolC){
      done=true;
    }
    expH2=real(valsE2);
    // save the MPS if needed
    if(done||(cnt%saveFreq==0)){
      string newMPSfile=mpsFileName(mpsdir,mpsfile,cnt);
      gs.exportMPS(newMPSfile.data());
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


const string mpsFileName(const string& mpsdir,const string& baseName,int cnt){
  stringstream s;
  s<<mpsdir<<"/"<<baseName;
  if(cnt>0) s<<"_"<<cnt;
  return s.str();
}

void tmpSave(MPS& rho,const string& baseDir,const string& baseName,int cnt,string& tmpFile){
  // first create a new tmpFile name
  string newTmp=mpsFileName(baseDir,baseName,cnt);
  rho.exportMPS(newTmp.data());
  // now remove the old one
  if(tmpFile.length()>0&&file_exists(tmpFile)){
    remove(tmpFile.data());
  }
  // and replace the string
  tmpFile=newTmp;
}

int findTmpFile(int M,const string& mpsdir,const string& baseName,string& mpsfile){
  int cnt=M;
  bool found=0;
  while(cnt>0&&!found){
    mpsfile=mpsFileName(mpsdir,baseName,cnt);
    if(file_exists(mpsfile)){
      cout<<"findTmpFile found mpsfile="<<mpsfile<<endl;
      found=true;
      return cnt;
    }
    else{
      mpsfile="";
      cnt--;
    }
  }
  if(!found) return 0;
  cout<<"ERROR: Should not be here"<<endl;
  exit(1);
}

int gcd(int a,int b){
  cout<<"Computing gcd("<<a<<","<<b<<")";
  int p1=max(a,b);
  int p2=min(a,b);
  while(p2>0){
    int aux=(p1%p2);
    p1=max(aux,p2);
    p2=min(aux,p2);
  }
  cout<<p1<<endl;
  if(p2==0) return p1;
  else{
    cout<<"What is going on here??"<<endl;exit(1);
  }
}

int lcm(int a,int b){
  return (a*b)/gcd(a,b);
}
