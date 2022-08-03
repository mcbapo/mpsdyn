
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

#include "IsingHamiltonian.h"

using namespace shrt;

/** To find the GS of H^2, instead of running findGroundState,
    it applies repeatedly cos(x H/N), which approximates the 
    effect of imaginary time evolution of H^2.
    In this version, the final result for M steps is approximated
    as the sum of O(sqrt(M)) terms.
    The Hamiltonian is Ising Hamiltonian \ref <IsingHamiltonian> 
    After running, compute the exp value of H^2 for all cuts

    Receives arguments:
    \param <propsfile> (string) mandatory: name of config props file where to read the rest. The rest can be set in the command line
    \param <L> (int) length of the chain
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D1> (int) maximum bond dimension
    \param <D2> (int) maximum bond dimension
    \param <deltaD> (int) step change in bond dim
    \param <dt> (double) time step (for Trotter)
    \param <x> (double) paramenter multiplying the H in the argument of the cosine
    \param <outfname> (char*) name of the output file for the energy
    \param <mpsdir> (char*) directory where to store the tmp (or final)
                     state, so that I can continue the evolution
    \param <mpsfile> (char*) base of the filename where to save or from 
                     which to read (only if newInstance is false) the
		     tmp MPS. The full name will be changed to include 
		     the bond dimension it corresponds to. The basis of
		     thename should be such that we do not mix
		     evolutions with different parameters (e.g. changing step width).
    \param <app> (int) whether the output file is to be kept (app==1)
*/

void computeValues(const MPS& state,const MPO& hamil,vector<complex_t>& valH2,vector<complex_t>& valH);
int findTmpFile(int M,const string& mpsdir,const string& baseName,string& mpsfile);
const string mpsFileName(const string& mpsdir,const string& baseName,int cnt);
void initialState(MPS& mps,int intSt,int D);

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
  //string mpsdir=props.getProperty("mpsdir");
  string mpsfile=props.getProperty("mpsfile"); // save final MPS
  string initmpsfile=props.getProperty("initmpsfile"); // use for initial state
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  double deltaX=props.getDoubleProperty("x"); // cos(H x) is the operator (x=1/N*eps)
  double dt=props.getDoubleProperty("dt"); // Trotter step
  if(dt<0||dt>2*deltaX) dt=2*deltaX; // default both steps the same
  int mS=(int)(2*deltaX/dt); // if not really an integer divisor, adjust
  if(abs(2*deltaX-mS*dt)>1E-6){
    mS=ceil(2*deltaX/dt);dt=2*deltaX/mS;
    cout<<"Adjusted the Trotter step to dt="<<dt<<", i.e. "<<mS
	<<" steps for x="<<deltaX<<endl;
  }
  else{
    cout<<"x="<<deltaX<<" thus nr of trotter steps "<<mS<<endl;
  }
  int initSt=props.getIntProperty("initSt"); // initial state to use:
  //                               +1=X+; -1=X-; (2,3 for Y, Z)
  //                               +4,+5,+6 are the staggered(+-) for X, Y and Z and -4,-5,-6 the -+
  double a=props.getDoubleProperty("a"); // factor such that a sqrt(M) terms are used in the sum
  int M=props.getIntProperty("M");
  int nSteps=ceil(a*sqrt(M));
  int rate=props.getIntProperty("rate");
  if(rate<0) rate=1;
  // int saveFreq=props.getIntProperty("saveFreq");
  // if(saveFreq<0) saveFreq=M;
  // saveFreq=lcm(rate,saveFreq);
  
  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
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
  // Now initialize the MPS
  bool realState=0;
  if(!initmpsfile.empty()){
    gs.importMPS(initmpsfile.data());
  }
  else{
    initialState(gs,initSt,D);
    realState=(abs(initSt)!=2)&&(abs(initSt)!=5); // all except Ys
  }
  // In this version, if I save something it should be the accumulated
  // state, plus both + and - evolved. Then I would need to recover
  // all of them consistently. So, for the moment I do not save them.
  double time=0.;
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);


  MPS auxP(gs),auxM(gs); // for plus and minus evol
  
  // First: H
  IsingHamiltonian hamI(L,d,J,g,h); 
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamI.getHMPO(0.);
  cout<<"Constructed the hamil MPO"<<endl;

  cout<<"Initial state has <H>="<<contractor.contract(gs,hamil,gs)<<", <H^2>="<<contractor.contract2(hamil,gs)<<endl;

  // Run the evol step with positive and negatieve exponentials separately
  MPO expP(L),expM(L); // positive and negative exponentials of the commutator
  // get the idividual MPOs from the Hamiltonian
  hamI.getUMPO(expP,dt*I_c,0); // the single layer positive exponential
  hamI.getUMPO(expM,-dt*I_c,0); // the negative exponential

  bool done=0;
  double expH2=L*100.;
  // place for other params
  //  complex_t E,EL,varL,ER,varR;
  complex_t valsE,valsE2;
  valsE=contractor.contract(gs,hamil,gs);
  valsE2=contractor.contract2(hamil,gs);
  cout<<"For step "<<cnt<<", <H^2>="<<valsE2<<" <H>="<<valsE<<endl;
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
  double factor=0.; // factor with which we sum the contribs of higher terms
  // instead of findGS, I run cos(xH) for up to M steps and approximate it as:
  //cos(x H eps)^M~1/2^M sum_{}
  int stepsLeft=nSteps-cnt;
  double normLast=log(gs.getNormFact());
  while(stepsLeft>0){
    int nrToApply=min(stepsLeft,rate);
    for(int k=0;k<nrToApply;k++){
      // If mS>1, need to apply several operators of eqch sign before summing
      for(int k2=0;k2<mS;k2++){
	MPS _aux1(auxP),_aux2(auxM);
	contractor.optimize(expP,_aux1,auxP,D);
	contractor.optimize(expM,_aux2,auxM,D);
      }
      auxP.gaugeCond('R',1);auxM.gaugeCond('R',1); // these have norm 1 (if exact)
      cnt++;stepsLeft--;time+=deltaX;
      factor+=(log(.5*M-cnt+1.)-log(.5*M+cnt)-normLast);
      //      cout<<"At step "<<cnt<<", log(factor)="<<factor<<", exp()="<<exp(factor)<<endl;
      // cout<<"Check overlaps:\n\t <state|vecPlus>="<<contractor.contract(auxP,gs)<<endl;
      // cout<<"\t <state|vecMinus>="<<contractor.contract(auxM,gs)<<endl;
      // cout<<"\t <vecPlus|vecMinus>="<<contractor.contract(auxM,auxP)<<endl;
      //      cout<<"M="<<M<<", M/2="<<.5*M<<", cnt="<<cnt<<", (M/2-cnt+1)="<<.5*M-cnt+1.<<endl;
      //cout<<"log(eso)="<<log(.5*M-cnt+1.)<<", (M/2+cnt)="<<.5*M+cnt<<", log(eso)="<<log(.5*M+cnt)<<endl;
      MPS auxGS(gs);
      vector<const MPS*> vecs;vecs.push_back(&auxGS);vecs.push_back(&auxP);vecs.push_back(&auxM);
      vector<complex_t> coeffs;
      coeffs.push_back(ONE_c);coeffs.push_back(exp(factor)*ONE_c);coeffs.push_back(exp(factor)*ONE_c);
      contractor.optimizeSum(vecs,coeffs,gs,D);
      // cout<<"After optimizeSum, now <state|vecPlus>="<<contractor.contract(auxP,gs)<<endl;
      // cout<<"\t <state|vecMinus>="<<contractor.contract(auxM,gs)<<endl;
      // double normGS=real(contractor.contract(gs,gs));
      // double normP=real(contractor.contract(auxP,auxP));
      // double normM=real(contractor.contract(auxM,auxM));
      // cout<<"And the same, but normalizing: <state|vecPlus>="<<contractor.contract(auxP,gs)/sqrt(normGS*normP)<<endl;
      // cout<<"\t <state|vecMinus>="<<contractor.contract(auxM,gs)/sqrt(normGS*normM)<<endl;
      gs.gaugeCond('R',0);
      //normLast+=log(gs.getNormFact());gs.setNormFact(1.);
    }
    //    exit(1);
    // Now record values
    valsE=contractor.contract(gs,hamil,gs);
    valsE2=contractor.contract2(hamil,gs);
    complex_t normV=contractor.contract(gs,gs);
    cout<<"State for step "<<cnt<<" found with <H^2>="<<valsE2/normV<<" <H>="<<valsE/normV<<", norm now="<<normV<<endl;
    out=new ofstream(outfname.data(),ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<setprecision(15);
    *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"<<dt<<"\t"<<deltaX<<"\t";
    *out<<cnt<<"\t"<<time<<"\t"<<real(valsE2)/real(normV)<<"\t"<<real(valsE)/real(normV);
    *out<<endl;
    out->close(); delete out;
    // If the added weight is too small, stop early
    if(exp(factor)<1E-4){
      cout<<"Stoppping the iteration because added weight is small: "<<exp(factor)<<endl;
      stepsLeft=0;
    }
  }
  if(!mpsfile.empty()){
    gs.exportMPS(mpsfile.data());
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

void initialState(MPS& mps,int initSt,int D){
  int L=mps.getLength();
  if(initSt==+1||abs(initSt)==4) mps.setProductState(p_xplus);
  else{
    if(initSt==-1) mps.setProductState(p_xminus);
    else{
      if(initSt==+2||abs(initSt)==5) mps.setProductState(p_yplus);
      else{
	if(initSt==-2) mps.setProductState(p_yminus);
	else{
	  if(initSt==+3||abs(initSt)==6) mps.setProductState(p_zero);
	  else{
	    if(initSt==-3) mps.setProductState(p_one);
	    else{ // any other number, I start random
	      cout<<"Initial state a non-translational random MPS with D="<<mps.getBond()<<endl;
	      mps.setRandomState();mps.gaugeCond('R',1);
	    }
	  }
	}
      }
    }
  }
  mps.increaseBondDimension(D);
  // Now for the staggered, change one every two by applying a single site op (sigma_z to x,y sigma_y to z)  
  if(abs(initSt)>=4&&abs(initSt)<=6){
    int firstPos=initSt<0?0:1; //+- changes first site 1
    complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
    mwArray oper;
    if(abs(initSt)<6) oper=mwArray(Indices(d,d),dataz);
    else oper=mwArray(Indices(d,d),datax);
    for(int k=firstPos;k<L;k+=2){
	mps.applyLocalOperator(k,oper);
    }
  }
}
