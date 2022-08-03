
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

#define REALPSI 1

using namespace shrt;

/** To find the GS of H^2, instead of running findGroundState,
    it applies sequentially cos(x*H*m), for x=eps/N and m=1 to maxM,
    (and the step can be exponential in m)
    which approximates the effect of imaginary time evolution of H^2.
    The Hamiltonian is Ising Hamiltonian \ref <IsingHamiltonian> 
    After running, compute the exp value of H^2

    Receives arguments:
    \param <propsfile> (string) mandatory: name of config props file where to read the rest. The rest can be set in the command line
    \param <L> (int) length of the chain
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D> (int) bond dimension
    \param <dt> (double) time step (for Trotter)
    \param <x> (double) paramenter multiplying the H in the argument of the cosine
    \param <M> (int) maximum step to apply
    \param <exponential> (bool) whether the step for increasing the cosine is 2^m (default false)
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
    \param <saveFreq> (int) frequency to save MPS to file, if any
*/

void computeValues(const MPS& state,const MPO& hamil,vector<complex_t>& valH2,vector<complex_t>& valH);
int findTmpFile(int M,const string& mpsdir,const string& baseName,string& mpsfile);
const string mpsFileName(const string& mpsdir,const string& baseName,int cnt);

void expandRDM(const MPS& mps,int posL,int posR,mwArray& A);
double computeTraceDist(const MPS& gs,int posL,int posR);
int factorCos(int m,bool expon);

//int lcm(int a,int b);

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
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  double deltaX=props.getDoubleProperty("x");
  double dt=props.getDoubleProperty("dt");
  double tolC=props.getDoubleProperty("convTol"); // relative tolerance for convergence of <H^2>
  if(tolC<0) tolC=1E-5;
  int stepsTol=props.getIntProperty("stepsConvTol"); // nr of steps such
						  // taht is after
						  // these many, the
						  // relative change
						  // in H^2 is always
						  // smaller than
						  // tolC, it stops.
  if(stepsTol<0) stepsTol=1;
  if(dt<0||dt>deltaX) dt=deltaX; // default both steps the same
  int mS=(int)(deltaX/dt); // if not really an integer divisor, adjust
  if(abs(deltaX-mS*dt)>1E-6){
    mS=ceil(deltaX/dt);dt=deltaX/mS;
    cout<<"Adjusted the Trotter step to dt="<<dt<<", i.e. "<<mS
	<<" steps for x="<<deltaX<<endl;
  }
  int M=props.getIntProperty("M");
  bool expon=(props.getIntProperty("exponential")>0);
  int saveFreq=props.getIntProperty("saveFreq");
  if(saveFreq<0) saveFreq=M;
  
  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% N\t J\t g\t h\t D\t dt\t deltaX\t cnt\t<H2>\t <H>\t S(L/2)\t Lc \t dist(rho_Lc)..."<<endl;
    out->close();
    delete out;
  }
  cout<<setprecision(10);

  int cnt=0;
  Contractor& contractor=Contractor::theContractor();
  //  contractor.setEigenSolver(primme);
  cout<<"Initialized Contractor"<<endl;
  MPS gs(L,D,d);
  // Now if we start from scratch, initialize the MPS. Otherwise, read
  // the tmp file
  string tmpMPSfile;
  double time=0.;
  cnt=findTmpFile(M,mpsdir,mpsfile,tmpMPSfile);
  if(cnt>0)
    cout<<"Recovered state for step="<<cnt<<" at file "<<tmpMPSfile<<endl;
  if(cnt==0){ // init from scratch
    gs.setProductState(p_xplus);
    gs.gaugeCond('L',1);
    cout<<"Initialized initial X+ product state, norm "<<contractor.contract(gs,gs)<<endl;
    //gs.setRandomState(); // the intial state, random
    //cout<<"Initialized random state from scratch, norm "<<contractor.contract(gs,gs)<<endl;
  }
  else{ 
    gs.importMPS(tmpMPSfile.data());
    cout<<"Recovered state at step "<<cnt<<endl;
  }
  gs.gaugeCond('R',1);
  
  // First: H
  IsingHamiltonian hamI(L,d,J,g,h); 
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamI.getHMPO(0.);
  cout<<"Constructed the hamil MPO"<<endl;

  cout<<"Initial state has <H>="<<contractor.contract(gs,hamil,gs)<<", <H^2>="<<contractor.contract2(hamil,gs)<<endl;

  // Run the evol step with positive and negatieve exponentials separately
  // get the idividual MPOs from the Hamiltonian
  MPO expP(L);
  hamI.getUMPO(expP,dt*I_c,0); // the single layer positive exponential
#ifndef REALPSI
  MPO expM(L); 
  hamI.getUMPO(expM,-dt*I_c,0); // the negative exponential
#endif
  
  bool done=0;int cntTol=0;
  double expH2=L*100.;
  // place for other params
  //  complex_t E,EL,varL,ER,varR;
  complex_t valsE,valsE2;
  
  // instead of findGS, I run cos(xH m ) incrementing m from 1 to M (or 2^m)
  while(cnt<=M&&!done){
    // Now record values
    valsE=contractor.contract(gs,hamil,gs);
    valsE2=contractor.contract2(hamil,gs);
    //    cout<<"State for step "<<cnt<<" found with <H^2>="<<valsE2<<" <H>="<<valsE<<endl;
    // Prepare right gauge to compute trace distance!
    int minPosR=L/2; // assuming L even and smallest Lc is 2
    for(int k=L-1;k>=minPosR+1;k--){
      gs.gaugeCond(k,'L',true);
    }
    int Lc=8;
    int posL=(L-Lc)/2;int posR=posL+Lc-1; // outside tracing them
    vector<double> distances;vector<int> Lcs;
    while(Lc>=2){
      double trDist=computeTraceDist(gs,posL,posR);
      distances.push_back(trDist);Lcs.push_back(Lc);
      posL++;posR--;Lc=(posR-posL+1);
    }
    cout<<"State for step "<<cnt<<" found with <H^2>="<<valsE2<<" <H>="<<valsE<<", dist(Lc=2)="<<distances[distances.size()-1]<<endl;
    // Compute entropy too
    double entropy=contractor.getEntropy(gs);
    // and write everything
    out=new ofstream(outfname.data(),ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<setprecision(15);
    *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"<<dt<<"\t"<<deltaX<<"\t";
    *out<<factorCos(cnt,expon)<<"\t"<<real(valsE2)<<"\t"<<real(valsE)<<"\t"<<entropy<<"\t";
    for(int l=0;l<Lcs.size();l++)
      *out<<Lcs[l]<<"\t"<<distances[l]<<"\t";
    *out<<endl;
    out->close(); delete out;

    cnt++; // have to apply the step m=cnt
    if(cnt<=M){
      // Apply one step with m: decide if dt has to change, if so,cahnge evolution operators
      // TODO THIS!!!!^^^^
      int newF=factorCos(cnt,expon);
      int mS_=newF*mS;
      // If mS>1, need to apply several operators of each sign before summing
      MPS aux1(gs);
      for(int k=0;k<mS_;k++){
	MPS _aux1(aux1);
	contractor.optimize(expP,_aux1,aux1,D);
      }
      aux1.gaugeCond('R',1);
#ifndef REALPSI
      MPS aux2(gs);
      for(int k=0;k<mS_;k++){
	MPS _aux2(aux2);
	contractor.optimize(expM,_aux2,aux2,D);
      }
      aux2.gaugeCond('R',1);
#else
      // copy and conjugate
      MPS aux2(aux1);
      for(int k=0;k<L;k++)
	aux2.setRotatedA(k,aux2.getA(k).getA(),Indices(1,2,3),true); // simply conjugate
#endif
      vector<const MPS*> vecs;vecs.push_back(&aux1);vecs.push_back(&aux2);
      vector<complex_t> coeffs(2,.5*ONE_c);
      contractor.optimizeSum(vecs,coeffs,gs,D);
      gs.gaugeCond('R',1);

      // And if needed, save the MPS to file
      if(cnt%saveFreq==0){
	string newMPSfile=mpsFileName(mpsdir,mpsfile,cnt);
	gs.exportMPS(newMPSfile.data());
      }
    }
    // check for doneness: if variance didn't change more than tolC
    if(abs(1.-real(valsE2)/expH2)<=tolC) cntTol++;
    else cntTol=0;
    done=(cntTol>=stepsTol);
    expH2=real(valsE2);
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

// int gcd(int a,int b){
//   cout<<"Computing gcd("<<a<<","<<b<<")";
//   int p1=max(a,b);
//   int p2=min(a,b);
//   while(p2>0){
//     int aux=(p1%p2);
//     p1=max(aux,p2);
//     p2=min(aux,p2);
//   }
//   cout<<p1<<endl;
//   if(p2==0) return p1;
//   else{
//     cout<<"What is going on here??"<<endl;exit(1);
//   }
// }

// int lcm(int a,int b){
//   return (a*b)/gcd(a,b);
// }

void expandRDM(const MPS& mps,int posL,int posR,mwArray& A){
  A=mps.getA(posL).getA();
  int d=A.getDimension(0);int Dl0=A.getDimension(1);int Dr=A.getDimension(2);
  A.permute(Indices(2,1,3));A.reshape(Indices(Dl0*d,Dr));
  int dphys=d; // temporary dimensions of the vector
  int Dv=Dr;
  for(int k=posL+1;k<=posR;k++){
    mwArray aux=mps.getA(k).getA();
    //cout<<"Multiplying tensor for site "<<k<<", A="<<aux.getDimensions()<<" with tmp result "<<A.getDimensions()<<endl;
    Indices dims=aux.getDimensions();
    int d=dims[0];int Dl=dims[1];int Dr=dims[2];
    if(Dl!=Dv){
      cout<<"Error: Dimensions do not agree in expandRDM!"<<endl;
      exit(1);
    }
    aux.permute(Indices(2,1,3));
    aux.reshape(Indices(Dl,d*Dr));
    A.multiplyRight(aux);
    dphys=dphys*d;
    Dv=Dr;
    A.reshape(Indices(Dl0*dphys,Dv));
  }
  A.reshape(Indices(Dl0,dphys,Dv));A.permute(Indices(2,1,3));A.reshape(Indices(dphys,Dl0*Dv));
  // if the mps has a normalization factor, we need to include it!!
  A=mps.getNormFact()*A;
  A.multiplyRight(Hconjugate(A));
}  

double computeTraceDist(const MPS& gs,int posL,int posR){
  mwArray oper;
  expandRDM(gs,posL,posR,oper);
  int Lc=posR-posL+1;
  //  cout<<"Created RDM for Lc="<<Lc<<endl;
  double singvals=0.;
  double idFc=1./pow(2,Lc);
  mwArray U;vector<complex_t> Dv;
  wrapper::eig(oper,Dv,U); // only eigenvalues
  //	  cout<<"Diagonalized the rho(Lc)"<<endl;
  // now sum the absolute values
  for(int k=0;k<Dv.size();k++){
    double aux=real(Dv[k]);
    singvals+=abs(aux-idFc);
    if(Dv.size()<pow(2,Lc)){
      cout<<"WARNING! Seems the dimension of the operator for "
	  <<"Lc="<<Lc<<"is not full: dim(oper)="<<Dv.size()
	  <<", 2^Lc="<<pow(2,Lc)<<endl;
      singvals+=idFc*(pow(2,Lc)-Dv.size());
    }
  }
  return singvals;
}

int factorCos(int m,bool expon){
  //return 2;
  if(!expon) return m;
  else return pow(2,m); // TODO add a random part?
}
