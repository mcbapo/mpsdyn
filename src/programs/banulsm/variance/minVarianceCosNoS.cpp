
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
// #define TEST_RE 1

#include "IsingHamiltonian.h"

using namespace shrt;

/** To find the GS of H^2, instead of running findGroundState,
    it applies repeatedly cos(x H), which approximates the 
    effect of imaginary time evolution of (H/N)^2.
    The Hamiltonian is Ising Hamiltonian \ref <IsingHamiltonian> 
    After running, compute the exp value of H^2 for all cuts

    Receives arguments:
    \param <propsfile> (string) mandatory: name of config props file where to read the rest. The rest can be set in the command line
    \param <L> (int) length of the chain
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <dt> (double) time step (for Trotter)
    \param <x> (double) paramenter multiplying the H in the argument of the cosine
    \param <M> (int) maximum number of steps
    \param <rate> (int) rate to store results
    \param <outfname> (char*) name of the output file for the energy and variance
    \param <mpsdir> (char*) directory where to store the tmp (or final)
                     state, so that I can continue the evolution
    \param <mpsfile> (char*) base of the filename where to save or from 
                     which to read (only if newInstance is false) the
		     tmp MPS. The full name will be changed to include 
		     the bond dimension it corresponds to. The basis of
		     thename should be such that we do not mix
		     evolutions with different parameters (e.g. changing step width).
    \param <saveFreq> (int) when to store intermediate MPS data (default at the end only)
    \param <app> (int) whether the output file is to be kept (app==1)
*/

void computeValues(const MPS& state,const MPO& hamil,vector<complex_t>& valH2,vector<complex_t>& valH);
int findTmpFile(int M,const string& mpsdir,const string& baseName,string& mpsfile);
const string mpsFileName(const string& mpsdir,const string& baseName,int cnt);
void expandRDM(const MPS& mps,int posL,int posR,mwArray& A); // obsolete!
double computeRDMdistanceToId(const mwArray& rdm,int Lc);
void initialState(MPS& mps,int intSt,int D);



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
  // double noise = props.getDoubleProperty("noise");
  // if(noise<0) noise=0;
  int D = props.getIntProperty("D");  
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir");
  string mpsfile=props.getProperty("mpsfile");
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  double deltaX=props.getDoubleProperty("x");
  double dt=props.getDoubleProperty("dt");
  if(dt<0||dt>deltaX) dt=deltaX; // default both steps the same
  int mS=(int)(deltaX/dt); // if not really an integer divisor, adjust
  if(abs(deltaX-mS*dt)>1E-6){
    mS=ceil(deltaX/dt);dt=deltaX/mS;
    cout<<"Adjusted the Trotter step to dt="<<dt<<", i.e. "<<mS
	<<" steps for x="<<deltaX<<endl;
  }
  int initSt=props.getIntProperty("initSt"); // initial state to use:
  //                               +1=X+; -1=X-; (2,3 for Y, Z)
  //                               +4,+5,+6 are the staggered(+-) for X, Y and Z and -4,-5,-6 the -+
  //                               +7 is a (complex) random state and +8 a real random (in which cases, the random seed may be given)
  //                               +9 is the real part of the Yplus state, to have a real state that has zero energy
  // These three parameters are obsolete
  bool initStag=(abs(initSt)>=4&&abs(initSt)<=6); //(props.getIntProperty("staggered")>0);
  bool initRand=(initSt>=7); //(props.getIntProperty("random")>0);
  bool initY=(initSt==+2); //(props.getIntProperty("initY")>0);
  bool realInit=(abs(initSt)!=2&&abs(initSt)!=5&&initSt!=7);
  int randSeed=props.getIntProperty("seed");
  int M=props.getIntProperty("M");
  int rate=props.getIntProperty("rate");
  bool checkRDM=(props.getIntProperty("checkRDM")>0); // whether to check how close RDM are to thermal
  if(rate<0) rate=1;
  int saveFreq=props.getIntProperty("saveFreq");
  if(saveFreq<0) saveFreq=M;
  saveFreq=lcm(rate,saveFreq);
  int lastSavedCnt=-1;
  int Lcmax=props.getIntProperty("Lcmax");
  if(Lcmax<0||Lcmax>10) Lcmax=6;
  int Lb=props.getIntProperty("Lb");
  if(Lb<0) Lb=6;
  bool avrg=1;
  avrg=!(props.getIntProperty("average")==0); // whether to average for each Lc subsystem
  
  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% N\t J\t g\t h\t D\t dt\t deltaX\t cnt\t time \t<H2>\t <H>";
    if(checkRDM){
      for(int kl=Lcmax;kl>=2;kl-=2){
	*out<<"\t"<<kl<<"\t dist(rho_"<<kl<<")";
      }
    }
    *out<<endl;
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
  if(cnt==0){ // init from scratch
    initialState(gs,initSt,D);
    cout<<"Initialized state ("<<initSt<<") from scratch, norm "<<contractor.contract(gs,gs)<<endl;
  }
  else{ 
    cout<<"Recovering state for step="<<cnt<<" from file "<<tmpMPSfile<<endl;
    gs.importMPS(tmpMPSfile.data());
  }
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  
  // First: H
  IsingHamiltonian hamI(L,d,J,g,h); 
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamI.getHMPO(0.);
  cout<<"Constructed the hamil MPO"<<endl;

  cout<<"Initial state has <H>="<<contractor.contract(gs,hamil,gs)<<", <H^2>="<<contractor.contract2(hamil,gs)<<endl;

  // Run the evol step with positive and negatieve exponentials separately
  MPO expP(L),expM(L); // double layer positive and negative exponentials of the commutator
  // get the idividual MPOs from the Hamiltonian
  hamI.getUMPO(expP,dt*I_c,0); // the single layer positive exponential
  if(!realInit)
    hamI.getUMPO(expM,-dt*I_c,0); // the negative exponential

  bool done=0;
  double expH2=L*100.;
  // place for other params
  //  complex_t E,EL,varL,ER,varR;
  complex_t valsE,valsE2;
  // instead of findGS, I run cos(xH) for up to M steps
  int stepsLeft=M-cnt;
  clock_t start,finish,finishSum;
  while(stepsLeft>0){
    int nrToApply=min(stepsLeft,rate);
    start=clock();int stepsInTime=nrToApply;
    for(int k=0;k<nrToApply;k++){
      // If mS>1, need to apply several operators of eqch sign before summing
      MPS aux1(gs);
      MPS aux2(gs);
      for(int k=0;k<mS;k++){
	MPS _aux1(aux1);
	contractor.optimize(expP,_aux1,aux1,D);
	if(!realInit){
	  MPS _aux2(aux2);
	  contractor.optimize(expM,_aux2,aux2,D);
	}
      }
      aux1.gaugeCond('R',1);
      if(!realInit){
	aux2.gaugeCond('R',1);
      }
// #ifdef TEST_RE
//       else{
// 	cout<<"Taking for GS the real part only!!"<<endl;
// 	// Test: take just real part of everything! Cannot work!!
// 	for(int pos=0;pos<L;pos++){
// 	  mwArray auxTensor=aux1.getA(pos).getA();
// 	  gs.replaceSite(pos,.5*(auxTensor+conjugate(auxTensor)));
// 	}	
// 	// End test
//       }
// #else
      else{ // aux2 is just the complex conjugate of aux1
	//cout<<"Taking for GS the real part from aux1 and aux2"<<endl;
	for(int pos=0;pos<L;pos++){
	  aux2.replaceSite(pos,conjugate(aux1.getA(pos).getA()));
	}
	// I should not need gaugeCond again
      }
      finish=clock();
      vector<const MPS*> vecs;vecs.push_back(&aux1);vecs.push_back(&aux2);
      vector<complex_t> coeffs(2,.5*ONE_c);
      contractor.optimizeSum(vecs,coeffs,gs,D);
// #endif
      cnt++;stepsLeft--;time+=deltaX;
      gs.gaugeCond('R',1);
      finishSum=clock();	   
    }
    // Now record values
    valsE=contractor.contract(gs,hamil,gs);
    valsE2=contractor.contract2(hamil,gs);
    // if necessary, also RDM
    vector<double> trDistRDM;
    vector<int> sizeRDM;
    if(checkRDM){
      if(!avrg){
	int minPosR=L/2; // assuming L even and smallest Lc is 2
	for(int k=L-1;k>=minPosR+1;k--){
	  gs.gaugeCond(k,'L',true);
	}
	int Lc=Lcmax;
	int posL=(L-Lc)/2;int posR=posL+Lc-1; // outside tracing them
	while(Lc>=2){
	  mwArray oper;
	  expandRDM(gs,posL,posR,oper);
	  //cout<<"Created RDM for Lc="<<Lc<<endl;
	  double singvals=0.;
	  double idFc=1./pow(2,Lc);
	  mwArray U;vector<complex_t> Dv;
	  wrapper::eig(oper,Dv,U); // only eigenvalues
	  //	  cout<<"Diagonalized the rho(Lc)"<<endl;
	  // now sum the absolute values
	  for(int k=0;k<Dv.size();k++){
	    double aux=real(Dv[k]);
	    singvals+=abs(aux-idFc);
	  }
	  if(Dv.size()<pow(2,Lc)){
	    cout<<"WARNING! Seems the dimension of the operator for D="<<D
		<<", Lc="<<Lc<<"is not full: dim(oper)="<<Dv.size()
		<<", 2^Lc="<<pow(2,Lc)<<endl;
	    singvals+=idFc*(pow(2,Lc)-Dv.size());
	  }
	  sizeRDM.push_back(Lc);
	  trDistRDM.push_back(singvals);
	  posL++;posR--;Lc=(posR-posL+1);
	}
      }
      else{ // averaging over different subsystems of the same size
	for(int lc=1;lc<=Lcmax;lc++){sizeRDM.push_back(lc);} // position lc-1
	trDistRDM=vector<double>(sizeRDM.size(),0.);
	vector<int> avrgCnt(sizeRDM.size(),0.);
	int posL=Lb;
	while(posL<L-Lb-1){ // at least there is a 2-site RDM to compute
	  int lc=Lcmax;
	  int posR=posL+lc-1;
	  while(posR>=L){ // cannot compute this one, need to reduce the size
	    lc--;
	    posR=posL+lc-1;
	  }
	  int dimL=pow(d,lc);
	  int dimR=dimL;
	  mwArray oper=contractor.getRDM(gs,posL,posR);
	  while(lc>=1){ // Compute the distance and trace another site
	    if(posR<L-Lb){
	      // cout<<"Computing distance for rmd("<<lc<<") sites, pos:"<<posL<<"-"<<posR<<endl;
	      double dist_lc=computeRDMdistanceToId(oper,lc);
	      trDistRDM[lc-1]+=dist_lc;
	      avrgCnt[lc-1]++; // to take the average later
	    }
	    if(lc>=1){
	      // trace out two sites
	      oper.reshape(Indices(dimL/d,d,dimR/d,d));
	      oper.permute(Indices(1,3,2,4));
	      dimL=dimL/d;dimR=dimR/d;
	      oper.reshape(Indices(dimL*dimR,d*d));
	      oper.multiplyRight(reshape(identityMatrix(d),Indices(d*d,1)));
	      oper.reshape(Indices(dimL,dimR));
	      posR--;
	    }
	    lc--;
	  }	    
	  posL++;
	}
	for(int lc=1;lc<=Lcmax;lc++) trDistRDM[lc-1]=trDistRDM[lc-1]/avrgCnt[lc-1];
      }
    }    
    cout<<"State for step "<<cnt<<" found with <H^2>="<<valsE2<<" <H>="<<valsE<<endl;
    // cout<<" time "<<(finishSum-start)/(float)CLOCKS_PER_SEC<<"(per "<<stepsInTime<<" evol step "
    // 	<<(1./(stepsInTime))*(finish-start)/(float)CLOCKS_PER_SEC
    // 	<<"; sum "<<(finishSum-finish)/(float)CLOCKS_PER_SEC<<")"<<endl;
     // cout<<"DDD\t"<<(finishSum-start)/(float)CLOCKS_PER_SEC<<"\t"<<stepsInTime<<"\t"<<(finish-start)/(float)CLOCKS_PER_SEC
     // 	<<"\t"<<(finish-finishSum)/(float)CLOCKS_PER_SEC<<")"<<endl;
    out=new ofstream(outfname.data(),ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<setprecision(15);
    *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"<<dt<<"\t"<<deltaX<<"\t";
    *out<<cnt<<"\t"<<time<<"\t"<<real(valsE2)<<"\t"<<real(valsE)<<"\t";
    if(checkRDM){
      for(int ks=0;ks<sizeRDM.size();ks++){
	*out<<sizeRDM[ks]<<"\t"<<trDistRDM[ks]<<"\t";
      }
    }
    *out<<endl;
    out->close(); delete out;
    // And if needed, save the MPS to file
    if(cnt%saveFreq==0){
      string newMPSfile=mpsFileName(mpsdir,mpsfile,cnt);
      gs.exportMPS(newMPSfile.data());
      // remove the previous one
      if(lastSavedCnt>0&&lastSavedCnt!=cnt){
	string newMPSfile=mpsFileName(mpsdir,mpsfile,lastSavedCnt);
	remove(newMPSfile.data());
      }
      lastSavedCnt=cnt;
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
      cout<<"Error: Dimensions do not agree in expandVec!"<<endl;
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


double computeRDMdistanceToId(const mwArray& rdm,int Lc){
  double singvals=0.;
  double idFc=1./pow(2,Lc);
  mwArray U;vector<complex_t> Dv;
  wrapper::eig(rdm,Dv,U); // only eigenvalues
  // now sum the absolute values
  for(int k=0;k<Dv.size();k++){
    double aux=real(Dv[k]);
    singvals+=abs(aux-idFc);
  }
  if(Dv.size()<pow(2,Lc)){
    cout<<"WARNING! Seems the dimension of the operator for Lc="<<Lc
	<<"is not full: dim(oper)="<<Dv.size()
	<<", 2^Lc="<<pow(2,Lc)<<endl;
    singvals+=idFc*(pow(2,Lc)-Dv.size());
  }
  return singvals;
}


void initialState(MPS& mps,int initSt,int D){
  int L=mps.getLength();
  if(initSt==+1||abs(initSt)==4) mps.setProductState(p_xplus);
  else{
    if(initSt==-1) mps.setProductState(p_xminus);
    else{
      if(initSt==+2||abs(initSt)==5||initSt==9) mps.setProductState(p_yplus);
      else{
	if(initSt==-2) mps.setProductState(p_yminus);
	else{
	  if(initSt==+3||abs(initSt)==6) mps.setProductState(p_zero);
	  else{
	    if(initSt==-3) mps.setProductState(p_one);
	    else{ // any other number, I start random
	      cout<<"Initial state a non-translational random MPS with D="<<mps.getBond()<<endl;
	      mps.setRandomState(); // This creates it with complex values
	      // If real, discard imaginary parts
	      if(initSt==8){
		for(int pos=0;pos<L;pos++){
		  mwArray aux=mps.getA(pos).getA();
		  mps.replaceSite(pos,.5*(aux+conjugate(aux)));
		}
	      }
	      mps.gaugeCond('R',1);
	    }
	  }
	}
      }
    }
  }
  // For the initSt=9 case, Re(|Yplus>) now construct the real part
  if(initSt==9){
    MPS auxMPS(mps);
    auxMPS.setProductState(p_yminus); // a trick for the conjugate!
    // for(int pos=0;pos<L;pos++){
    //   auxMPS.replaceSite(pos,conjugate(auxMPS.getA(pos).getA()));
    // }
    Contractor& contractor=Contractor::theContractor();
    MPS mpsC(mps);
    vector<const MPS*> vecs;vecs.push_back(&auxMPS);vecs.push_back(&mpsC);
    vector<complex_t> coeffs(2,.5*ONE_c);
    contractor.optimizeSum(vecs,coeffs,mps,max(D,2));
    mps.gaugeCond('R',1);    
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

