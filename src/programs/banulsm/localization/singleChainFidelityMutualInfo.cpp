#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "HeisenbergHamiltonian.h"

using namespace shrt;

/** Based on qubitSingleChain and mutualInfoSingleChain, this program 
    starts or continues the evolution using for the
    Hamiltonian the random parameters contained in the mwArray file
    provided as argument.
    It thus simulates the time evolution of an MPO under a
    disordered Heisenberg Hamiltonian \ref <HeisenbergHamiltonian>.
    A number of sites in the center of the chain are used to encode a single qubit, so that 
    \f$| \tilde{0} \rangle \langle \tilde{0}|\f$ corresponds to 
    \f$\left ( | 0 \rangle \langle 0| \right )^{\otimes L_0}\f$
    tensored with identities on both sides.
    Starting from an MPO which encodes a pure state, Z+, Z-, X+ or X-,
    the program simulates the evolution and computes the trace norm of the resulting 
    simga_Z or sigma_X evolved state when tracing out the edges, and keeping \f$L_0,L0+2,\ldots 10\f$
    spins. This will allow us to compute the recovery fidelity.

    To be able to control the error in the trace, instead of evolving
    the sigmas, I evolve two pure states and correct the trace oe each
    after every step so that it is one. Then I compute the trace norm
    of the difference between the reduced density matrix for states +
    and -.

    Additionally, we compute the (2-Renyi) mutual information for the
    central spins (from L0 to N-2) wrt the rest and across all cuts of
    the chain. The program can compute the average and localized
    magnetization along the way.  To keep track of how good the
    evolution is, I will also compute the energy along the way.
    

    Receives arguments as Properties file, as config/rfmi.conf, 
    containing the required parameters 
    \param <L> (int) number of sites
    \param <L0> (int) number of sites (in the original chain)
    \param <isXY> (bool) whether the model to be considered is only XY (no ZZ term) 
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian (Jx=Jy=Jz or Jx=Jy)
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian (local values 
                        randomly distributed between -h and h)			
    \param <isCont> (bool) whether the model to be considered has continuous randomness
    \param <paramFile> (char*) file with the coefficients for the magnetic field
    \param <nrInst> (int) Nr of row in the file to be used (if rows
                       are longer than L, the first L values of the row will be used).
    \param <init> (char) initial state in the middle, can be 1 (X),
                        or 3 (Z) (" for Y is equivalent to X)
    \param <delta> (double) width of the time step used
    \param <M> (int) number of steps
    \param <rate> (int) frequency of recording data (nr of steps)
                       Results will be saved for 0,rate,2*rate...,M
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the results
                   (will be appended if newInstance is 0 and it exists)
    \param <mpsdir> (char*) name of the dir where to read/save the MPS files (tmp and final)
    \param <mpsfileP> (char*) name of the file where to read/save the MPS for + state
    \param <mpsfileM> (char*) name of the file where to read/save the MPS for - state
    \param <outfname_traces> (char*) name of the output file for the trace norms results
    \param <savingFreq> (int) frequency (steps) to save the tmp MPS file
    \param <newInstance> (bool) whether to start a new evolution from
                         the inhomogeneous polarization distribution
			 (if not, the MPS in mpsfile is used)

    To be able to resume evolution after an interruption
    (e.g. exhausted time) I will save a temporary file with the same
    name as mpsfile, but a suffix indicating the latest advance. Every
    certain number of steps (normFreq), the temporary file is saved
    with the current count, and the previous one is removed. when the
    evolution starts, if newInstance is set to 0, the directory is
    checked for the existence of such file.

*/

/** Vectorize an MPO to make an MPS out of it. If up==true (default)
    the first component of the double physical index will be the one
    which was going up.*/
void MPSfromMPO(const MPO& mpo,MPS& mps,bool up=true);
void MPOfromMPS(const MPS& mps,MPO& mpo,bool up=true);

/** Set initial states to initSt=|alpha + - > (alpha=x or z), encoded in L0 spins (from pos0) */
void setExtendedInitStates(MPS& rho0p,MPS& rho0m,int L,int d,int initSt,int L0,int pos0);

void computeLocalOp(const MPS& rho,int pos,int d,ofstream* out);
void computeTraceNorms(const MPS& rho0p,const MPS& rho0m,const char* outfname,int L,int L0,bool isXY,double J_,double h_,int initSt,
		       double delta,int M,int D);
void constructLocal(int L,int pos,int initSt,int d,MPS& result);
double traceNormDiff(const MPS& rhop,const MPS& rhom);

/** Compute the 2-Renyi entropies of two pieces of the chain, namely
    A, defined to be from pos1A until pos2A (both included and
    numbered from 0 to L) and the complementary, B */
double compute2RenyiS(const MPS& redRho,complex_t traceR);
void compute2RenyiEntropiesAB(const MPS& rho,complex_t trR,int pos1A,int pos2A,double* SA,double*SB);
void computeEntropies(const MPS& rho,complex_t trR,ofstream* out,int L0,bool header=false);

void normalizeTrace(MPS& rho,complex_t trR);

/** 
Generates a file name for the temporary MPS 
*/
const string mpsFileName(const string& baseDir,const string& baseName,int cnt=0);

/** Saves the state to a tmp file with new name and destroys the previous one */
void tmpSave(MPS& rho,const string& baseDir,const string& baseName,int cnt,string& tmpFile);

/** Check if there are tmp files with earlier stages of evolution */
int findTmpFiles(int M,const string& baseDir,const string& baseNameP,const string& baseNameM,string& mpsFileP,string& mpsFileM);


int d=2;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  int L=props.getIntProperty("L");
  int L0=props.getIntProperty("L0");
  int _aux_=props.getIntProperty("isXY");
  bool isXY=_aux_>0;
  double J_=props.getDoubleProperty("J");
  double h_=props.getDoubleProperty("h");
  _aux_=props.getIntProperty("isCont");
  bool isCont=_aux_>0;
  string paramfname=props.getProperty("paramFile");
  int nrInst=props.getIntProperty("nrInst");
  int initSt=props.getIntProperty("initSt");
  double delta=props.getDoubleProperty("delta");
  int M=props.getIntProperty("M");
  int rate=props.getIntProperty("rate");
  int D=props.getIntProperty("D");
  string outfname=props.getProperty("outputResults");
  string outfname_Traces=props.getProperty("outputTraces");
  string mpsdir=props.getProperty("outputDirMPS");
  string mpsfileP=props.getProperty("mpsFileP");
  string mpsfileM=props.getProperty("mpsFileM");
  int savingFreq=props.getIntProperty("savingFreq");
  _aux_=props.getIntProperty("newInstance");
  bool newInstance=_aux_>0;
    

  cout<<"Initialized arguments "
      <<", L="<<L
      <<", L0="<<L0
      <<", isXY?="<<isXY
      <<", isCont?="<<isCont
      <<", J="<<J_
      <<", h="<<h_
      <<", delta="<<delta
      <<", M="<<M
      <<", D="<<D
      <<", outfile="<<outfname
      <<", paramfile="<<paramfname
      <<", nrInst="<<nrInst
      <<", mpsfile(+)="<<mpsfileP
      <<", mpsfile(-)="<<mpsfileM
      <<", newInstance="<<newInstance
      <<", savingFreq="<<savingFreq
      <<", rate="<<rate
      <<endl;

  // First of all, I try to locate any existing calculation (tmp MPS files) for 
  // step nr smaller than M
  string tmpFileM,tmpFileP;
  int M0=findTmpFiles(M,mpsdir,mpsfileP,mpsfileM,tmpFileP,tmpFileM);

  if(M0>0){ 
    cout<<"Found tmp files for step M0="<<M0<<endl;
    newInstance=false;
  }
  else newInstance=true;

  ofstream* out;
  bool appendingFile=(!newInstance&&file_exists(outfname));
  if(appendingFile)
    out=new ofstream(outfname.data(),ios::app);
  else // created again
    out=new ofstream(outfname.data());
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  vector<double> J(L-1,J_);
  vector<double> Jz(L-1,isXY?0.:J_);
  vector<double> h(L,h_);
  
  if(!isCont)
    cout<<"WARNING!! Turning to discrete values of the random fields!"<<endl;

  // Now read the parameters
  {
    mwArray H;
    ifstream inData(paramfname.data());
    if(!inData.is_open()){
      cout<<"Error: Cannot open the mwArray file "<<paramfname<<" to read data "<<endl;
    }
    H.load(inData);
    double aux_(0.);
    for(int k=0;k<L;k++){
      double haux(0.);
      H.getElement(haux,aux_,Indices(nrInst,k));
      h[k]=h_*haux;
      if(!isCont)
	h[k]=haux>0?h_:-h_;
    }
    // TODO: Check for right dimensions
    inData.close();
  }

  cout<<"Read h="<<h<<endl;  
  if(!appendingFile){
    *out<<"% L="<<L<<", L0="<<L0<<", J="<<J_<<", h="<<h_
	<<", isXY="<<isXY<<", initSt="<<initSt<<", D="<<D<<endl;
    *out<<"% hs in "<<paramfname<<", row "<<nrInst<<endl;
    *out<<"% Random values h:"<<endl;
    *out<<setprecision(18);
    for(int k=0;k<L;k++) *out<<"% h("<<k+1<<")="<<h[k]<<endl;
    *out<<setprecision(15);
    *out<<"% MPS dir is "<<mpsdir<<endl;
    *out<<"% MPS(+) in "<<mpsfileP<<endl;
    *out<<"% MPS(-) in "<<mpsfileM<<endl;
    // Now what the file contains 
    *out<<"% M\t delta\t t\t D\t state\t Re[tr(rho H)]\t Im[tr(rho H)]\t ";
    *out<<"S2(rho)\t";
    computeEntropies(MPS(L,1,1),I_c,out,L0,true);
    *out<<endl;
  }
  int Mtot=M;

  // First: create the Hamiltonian
  HeisenbergHamiltonian hamH(L,J,J,Jz,h,d);
  cout<<"Created the Hamiltonian"<<endl;
  //const MPO& hamilMPO=hamH.getHMPO();
  //hamilMPO.exportForMatlab("hamilMPO.m");
  //exit(1);

  // The evolution operator for a rho will result from double operators
  // Which MPOs do I need for the indices in both sides
  MPO expHe2(L),expHe(L),expHo(L);
  // I construct directly the double operators
  hamH.getDoubleExponentialMPOeven(expHe2,-delta*.5*I_c);
  hamH.getDoubleExponentialMPOeven(expHe,-delta*I_c);
  hamH.getDoubleExponentialMPOodd(expHo,-delta*I_c);

  // Now construct the initial rho MPO, in this case, the identity,
  // except for L0 sites 
  MPS rho0p(L,1,d*d),rho0m(L,1,d*d);
  // Now read the initial rho MPO from the file or construct it for new instance
  if(!newInstance){
    rho0p.importMPS(tmpFileP.data());
    rho0m.importMPS(tmpFileM.data());
  }
  else{
    //rho0.setProductState(p_maxent);
    setExtendedInitStates(rho0p,rho0m,L,d,initSt,L0,(L-L0)/2);
  }
  //MPO auxRho(L);
  //MPOfromMPS(rho0,auxRho);
  //auxRho.exportForMatlab("rho0.m");
  //exit(1);

  // Also the identity, to compute the trace
  MPS Id(L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++)
    Id.setA(k,id2);
  
  MPS hmps(L,1,d*d);
  const MPO& hamil=hamH.getHMPO();
  MPSfromMPO(hamil,hmps);
  
  //  M1.exportForMatlab("M1.m");
  //Id.exportForMatlab("Id.m");
  
  // Now do the evolution step by step, and compute the expectation value of M1 every time
  int cnt=M0;
  int posCenter=L/2; // for local opers
  double time=delta*cnt;
  complex_t trRp=contractor.contract(rho0p,Id);
  complex_t trRm=contractor.contract(rho0m,Id);
  // Renormalize, only after evolving, as the first should be ok
  cout<<"Before starting, traces of "<<initSt<<"+ and "<<initSt<<"- are "<<trRp<<" and "<<trRm<<endl;
  if(trRp!=ONE_c) normalizeTrace(rho0p,trRp);
  if(trRm!=ONE_c) normalizeTrace(rho0m,trRm);
  complex_t Eevp=contractor.contract(rho0p,hmps);
  complex_t Eevm=contractor.contract(rho0m,hmps);
  double S2Rp=compute2RenyiS(rho0p,ONE_c);
  double S2Rm=compute2RenyiS(rho0m,ONE_c);
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<initSt<<"\t"
      <<real(Eevp)<<"\t"<<imag(Eevp)<<"\t";
  *out<<S2Rp<<"\t";
  computeEntropies(rho0p,ONE_c,out,L0);
  *out<<endl;
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<-initSt<<"\t"
      <<real(Eevm)<<"\t"<<imag(Eevm)<<"\t";
  *out<<S2Rm<<"\t";
  computeEntropies(rho0m,ONE_c,out,L0);
  *out<<endl;
  out->close();delete out;

  computeTraceNorms(rho0p,rho0m,outfname_Traces.data(),L,L0,isXY,J_,h_,initSt,delta,cnt,D);
  while(cnt<Mtot){
    //cout<<"Time evolution step nr "<<cnt<<", time="<<time<<endl;
    //", value "<<M1ev<<", trace "<<trR<<endl;
    MPS* vecs[2]={&rho0p,&rho0m};
    int cntL=cnt;
    for(int nrSt=0;nrSt<2;nrSt++){
      MPS& rho0=*vecs[nrSt];
      MPS aux(rho0); // temporary copy
      contractor.optimize(expHe2,aux,rho0,D);
      // cout<<"After expHe2, v'*H*v="<<contractor.contract(imaGS,hamil,imaGS)
      // 	   <<", v*v="<<contractor.contract(imaGS,imaGS)<<endl;
      // Now apply r-1 times the pair Ho He
      int cntLoop=0;
      while(cntLoop<rate&&cnt<Mtot){
	contractor.optimize(expHo,rho0,aux,D);
	// cout<<"After expHo, v'*H*v="<<contractor.contract(imaGS,hamil,imaGS)
	// 	   <<", v*v="<<contractor.contract(imaGS,imaGS)<<endl;
	if(cntLoop<rate-1&&cnt<Mtot-1){
	  contractor.optimize(expHe,aux,rho0,D);
	}
	else{
	  contractor.optimize(expHe2,aux,rho0,D);
	}
	cnt++;cntLoop++;
	if(nrSt==1) time+=delta; // time increased just once
      }
      rho0.setNormFact(1);
      // now renormalize them
      complex_t trR=contractor.contract(rho0,Id);
      normalizeTrace(rho0,trR);
      if(nrSt==0) cnt=cntL; // reset the counter for the second state
    }

    Eevp=contractor.contract(rho0p,hmps);
    Eevm=contractor.contract(rho0m,hmps);
    out=new ofstream(outfname.data(),ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    // trRp=contractor.contract(rho0p,Id);
    // trRm=contractor.contract(rho0m,Id);
    trRp=ONE_c;trRm=ONE_c;
    S2Rp=compute2RenyiS(rho0p,trRp);
    S2Rm=compute2RenyiS(rho0m,trRm);
    *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<initSt<<"\t"
	<<real(Eevp)<<"\t"<<imag(Eevp)<<"\t";
    *out<<S2Rp<<"\t";
    computeEntropies(rho0p,trRm,out,L0);
    *out<<endl;
    *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<-initSt<<"\t"
	<<real(Eevm)<<"\t"<<imag(Eevm)<<"\t";
    *out<<S2Rm<<"\t";
    computeEntropies(rho0m,trRp,out,L0);
    *out<<endl;
    out->close();delete out;
    computeTraceNorms(rho0p,rho0m,outfname_Traces.data(),L,L0,isXY,J_,h_,initSt,delta,cnt,D);
    
    // export the MPS if cnt is multiple of savingFreq
    if(savingFreq&&cnt%savingFreq==0){
      tmpSave(rho0p,mpsdir,mpsfileP,cnt,tmpFileP);
      tmpSave(rho0m,mpsdir,mpsfileM,cnt,tmpFileM);
    }
  }

  rho0p.exportMPS(mpsFileName(mpsdir,mpsfileP).data());
  rho0m.exportMPS(mpsFileName(mpsdir,mpsfileM).data());

}

const string mpsFileName(const string& baseDir,const string& baseName,int cnt){
  stringstream s;
  s<<baseDir<<"/"<<baseName;
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

int findTmpFiles(int M,const string& baseDir,const string& baseNameP,const string& baseNameM,string& mpsFileP,string& mpsFileM){
  int cnt=M;
  bool found=0;
  while(cnt>0&&!found){
    mpsFileP=mpsFileName(baseDir,baseNameP,cnt);
    mpsFileM=mpsFileName(baseDir,baseNameM,cnt);
    if(file_exists(mpsFileP)&&file_exists(mpsFileM)){
      // only if both are found!
      cout<<"findTmpFiles found mpsFileP="<<mpsFileP<<" and mpsFileM="<<mpsFileM<<endl;
      found=true;
      return cnt;
    }
    else{
      mpsFileP="";mpsFileM="";
      cnt--;
    }
  }
  if(!found) return 0;
}


void normalizeTrace(MPS& rho0,complex_t trR){
  //  return;
  int L=rho0.getLength();
  double absTr=abs(trR);
  complex_t phase=(1/absTr)*trR;
  double factor=pow(1./absTr,double(1./L));
  rho0.setA(0,factor*conjugate(phase)*rho0.getA(0).getA());  
  for(int k=1;k<L;k++){
    rho0.setA(k,factor*rho0.getA(k).getA());
  }
}


void MPSfromMPO(const MPO& mpo,MPS& mps,bool up){
  int L=mpo.getLength();
  // vector<int> du=mpo.getDimensions();
  // vector<int> dd=mpo.getOriginalDimensions();
  // int D=mpo.getOp(0).getDr();

  // vector<int> newDims(du);
  // for(int k=0;k<L;k++) newDims[k]*=dd[k];

  mps=MPS(L,1,1);

  for(int k=0;k<L;k++){
    mwArray aux=mpo.getOp(k).getFullData();
    Indices dims=aux.getDimensions();
    if(up)
      aux.permute(Indices(1,3,2,4));
    else
      aux.permute(Indices(3,1,2,4));
    aux.reshape(Indices(dims[0]*dims[2],dims[1],dims[3]));
    mps.replaceSite(k,aux,0); // brute force
  }

}


void computeLocalOp(const MPS& rho,int pos,int d,ofstream* out){
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d*d,1,1),dataX);
  complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigY(Indices(d*d,1,1),dataY);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d*d,1,1),dataZ);

  // cout<<"Computing local expectation values "<<endl;
  mwArray sig0=identityMatrix(d);
  sig0.reshape(Indices(d*d,1,1));
  int L=rho.getLength();
  MPS oper(L,1,d*d);
  for(int k=0;k<L;k++){
    oper.setA(k,sig0);
  }

  Contractor& contractor=Contractor::theContractor();
  complex_t resL;
  oper.setA(pos,sigX);
  resL=contractor.contract(rho,oper);
  *out<<real(resL)<<"\t"<<imag(resL)<<"\t";
  oper.setA(pos,sigY);
  resL=contractor.contract(rho,oper);
  *out<<real(resL)<<"\t"<<imag(resL)<<"\t";
  oper.setA(pos,sigZ);
  resL=contractor.contract(rho,oper);
  *out<<real(resL)<<"\t"<<imag(resL)<<"\t";
}

void constructLocal(int L,int pos,int initSt,int d,MPS& result){
  mwArray sig0=identityMatrix(d);
  sig0.reshape(Indices(d*d,1,1));
  complex_t Xop[]={ZERO_c,ONE_c,ONE_c,ZERO_c}; 
  complex_t Yop[]={ZERO_c,I_c,-1.*I_c,ZERO_c}; 
  complex_t Zop[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c}; 
  mwArray X(Indices(d*d,1,1),Xop);
  mwArray Y(Indices(d*d,1,1),Yop);
  mwArray Z(Indices(d*d,1,1),Zop);
  mwArray* op;
  switch(initSt){
  case 1:
    { op=&X;
      break;}
  case 2:
    { op=&Y;
      break;}
  case 3:
    { op=&Z;
      break;}
  default:{
    cout<<"Error! Unknown initSt "<<initSt<<endl;
    exit(1);
  }
  }
  op->reshape(Indices(d*d,1,1));
  result=MPS(L,1,d*d);
  for(int k=0;k<L;k++){
    result.setA(k,sig0);
  }
  result.setA(pos,*op);
}


void computeTraceNorms(const MPS& rho0p,const MPS& rho0m,const char* outfname,int L,int L0,bool isXY,double J_,double h_,int initSt,
		       double delta,int M,int D){
  ofstream* out;
  bool appendingFile=file_exists(outfname)&&(M!=0);
  if(appendingFile)
    out=new ofstream(outfname,ios::app);
  else // created again
    out=new ofstream(outfname);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  if(!appendingFile){
    *out<<"% L="<<L<<", J="<<J_<<", h="<<h_
	<<" D="<<D<<", init St="<<initSt<<endl;
    *out<<"% M\t delta\t t\t D\t tr(sigma"<<initSt<<")\t |rho(t)_2|_1"<<endl;
  }
  *out<<setprecision(15);

  Contractor& contractor=Contractor::theContractor();
  //cout<<"Initialized Contractor"<<endl;
  // Also the identity, to compute the trace
  MPS Id(L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++)
    Id.setA(k,id2);

  MPS rhop(rho0p); // copying
  MPS rhom(rho0m); // copying
  // Shouldn't normalize here!!!! The trace was already one! (but just in case)
  complex_t trRp=contractor.contract(rhop,Id);
  complex_t trRm=contractor.contract(rhom,Id);
  if(abs(trRp-ONE_c)>1E-13) normalizeTrace(rhop,trRp);
  if(abs(trRm-ONE_c)>1E-13) normalizeTrace(rhom,trRm);

  *out<<M<<"\t"<<delta<<"\t"<<M*delta<<"\t"<<D<<"\t"; // 4 columns for nothing (3 is time)

  // Now I need to start tracing out things, to construct reduced operators and to compute their trace norms.
  // I will keep 10, 8,... L0+2, L0
  // Column 4+2 is 10; 8 is 8; 10 is 6; 12 is 4; 14 is 2

  for(int cut=5;cut>=L0/2;cut--){
    vector<int> toTrace;
    // 1) trace all but the central block of 2*cut spins
    int lastTracedL=max(L/2-cut-1,-1);int firstTracedR=min(L/2+cut,L);
    for(int k=0;k<=lastTracedL;k++){
      toTrace.push_back(k);
    }
    for(int k=firstTracedR;k<L;k++){
      toTrace.push_back(k);
    }
    //cout<<"To get the central "<<cut*2<<" spins (L="<<L<<"), tracing out "<<toTrace<<endl;
    MPS redRhoP(rhop);
    redRhoP.traceOutSites(toTrace);
    MPS redRhoM(rhom);
    redRhoM.traceOutSites(toTrace);
    double tN10=traceNormDiff(redRhoP,redRhoM);
    // And finally write the results, preceeded by the number kept!
    *out<<2*cut<<"\t";
    *out<<tN10<<"\t";
  }
  //cout<<"Computation finished"<<endl;
  *out<<endl;
  out->close();
  delete out;
}



void MPOfromMPS(const MPS& mps,MPO& mpo,bool up){
  int L=mps.getLength();

  mpo.initLength(L);

  for(int k=0;k<L;k++){
    mwArray aux=mps.getA(k).getA();
    Indices dims=aux.getDimensions();
    int d0=sqrt(dims[0]);
    if(d0*d0!=dims[0]){
      cout<<"Error: Dimension of site "<<k<<" ("<<dims[0]
	  <<") does not seem a square=> don't know how t divide it"<<endl;
      exit(1);
    }
    aux.reshape(Indices(d0,d0,dims[1],dims[2]));
    if(up)
      aux.permute(Indices(1,3,2,4));
    else
      aux.permute(Indices(2,3,1,4));
    mpo.setOp(k,new Operator(aux),true);
  }

}

double traceNormDiff(const MPS& rhop,const MPS& rhom){
  mwArray OpP,OpM;
  MPO rhoMPO(rhop.getLength());
  MPOfromMPS(rhop,rhoMPO);
  // if(rho.getLength()==4){
  //   rhoMPO.exportForMatlab("rhoTest4.m");
  //   exit(1);
  // }
  //cout<<"Computing trace norm of mps, length "<<rho.getLength()<<endl;
  expandOper(rhoMPO,OpP);
  MPOfromMPS(rhom,rhoMPO);
  expandOper(rhoMPO,OpM);
  // cout<<"Operator expanded to dims "<<Op.getDimensions()<<endl;
  // trace norm, is the sum of singular values
  mwArray U,S,Vdag;
  wrapper::svd(OpP-OpM,U,S,Vdag);
  return real(S.trace());
}

// void constructM1(int L,int d,MPS& result){
void setExtendedInitStates(MPS& rho0p,MPS& rho0m,int L,int d,int initSt,int L0,int pos0){
   mwArray sig0=identityMatrix(d);
   sig0.reshape(Indices(d*d,1,1));
   complex_t data00[]={ONE_c,ZERO_c,ZERO_c,ZERO_c};
   mwArray op00(Indices(d*d,1,1),data00);
   complex_t data01[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
   mwArray op01(Indices(d*d,1,1),data01);
   complex_t data10[]={ZERO_c,ONE_c,ZERO_c,ZERO_c};
   mwArray op10(Indices(d*d,1,1),data10);
   complex_t data11[]={ZERO_c,ZERO_c,ZERO_c,ONE_c};
   mwArray op11(Indices(d*d,1,1),data11);
   mwArray Z(Indices(4,d*d));
   for(int j=0;j<d*d;j++){
     Z.setElement(op00.getElement(Indices(j,0,0)),Indices(0,j));
     Z.setElement(op01.getElement(Indices(j,0,0)),Indices(1,j));
     Z.setElement(op10.getElement(Indices(j,0,0)),Indices(2,j));
     Z.setElement(op11.getElement(Indices(j,0,0)),Indices(3,j));
   }
   // Set the product of identities on the edges
   for(int k=0;k<pos0;k++){
     rho0p.replaceSite(k,.5*sig0,false);
     rho0m.replaceSite(k,.5*sig0,false);
   }
   for(int k=pos0+L0;k<L;k++){
     rho0p.replaceSite(k,.5*sig0,false);
     rho0m.replaceSite(k,.5*sig0,false);
   }
   if(initSt==3){
     for(int k=pos0;k<pos0+L0;k++){
       rho0p.replaceSite(k,op00,false);
       rho0m.replaceSite(k,op11,false);
     }
   }
   else{
     int D=4;// Dimension of the MPS
     mwArray C(Indices(D,D,4)),C1(Indices(1,D,4)),CL(Indices(D,1,4));
     for(int j=0;j<D;j++){
       C.setElement(ONE_c,Indices(j,j,j));
       C1.setElement(.5*ONE_c,Indices(0,j,j));
       CL.setElement(ONE_c,Indices(j,0,j));
     }
     mwArray C1m(C1); // for the minus sign, I put the sign in the first coefficient
     C1m.setElement(-.5*ONE_c,Indices(0,1,1));
     C1m.setElement(-.5*ONE_c,Indices(0,2,2));
     C1.reshape(Indices(1*D,4));C1.multiplyRight(Z);C1.reshape(Indices(1,D,d*d));C1.permute(Indices(3,1,2));
     C1m.reshape(Indices(1*D,4));C1m.multiplyRight(Z);C1m.reshape(Indices(1,D,d*d));C1m.permute(Indices(3,1,2));
     rho0p.replaceSite(pos0,C1,false);
     rho0m.replaceSite(pos0,C1m,false);
     C.reshape(Indices(D*D,4));C.multiplyRight(Z);C.reshape(Indices(D,D,d*d));C.permute(Indices(3,1,2));
     for(int k=pos0+1;k<pos0+L0-1;k++){
       rho0p.replaceSite(k,C,false);
       rho0m.replaceSite(k,C,false);
     }
     CL.reshape(Indices(D*1,4));CL.multiplyRight(Z);CL.reshape(Indices(D,1,d*d));CL.permute(Indices(3,1,2));
     rho0p.replaceSite(pos0+L0-1,CL,false);
     rho0m.replaceSite(pos0+L0-1,CL,false);
   }
}


/** Compute the 2-Renyi entropies of two pieces of the chain, namely
    A, defined to be from pos1A until pos2A (both included and
    numbered from 0 to L) and the complementary, B */
void compute2RenyiEntropiesAB(const MPS& rho,complex_t trR,int pos1A,int pos2A,double* SA,double*SB){
  int L=rho.getLength();
  // First trace out the complementary of A
  MPS rhoA(rho);
  vector<int> toTrace;
  for(int k=0;k<L;k++){
    if(k<pos1A||k>pos2A)
      toTrace.push_back(k);
  }
  rhoA.traceOutSites(toTrace);
  // Now compute its 2-Renyi entropy
  *SA=compute2RenyiS(rhoA,trR);
  // And now the same for the complementary
  toTrace.clear();
  for(int k=pos1A;k<=pos2A;k++)
    toTrace.push_back(k);
  rhoA=rho;
  rhoA.traceOutSites(toTrace);
  // Now compute its 2-Renyi entropy
  *SB=compute2RenyiS(rhoA,trR);  
}

double compute2RenyiS(const MPS& redRho,complex_t trR){
  // I need trace(rho^2)
  
  // If rho was Hermitian, as it ideally should,
  // trace(rho^2)=trace(rho+ rho) We will have problably some
  // numerical error, so that rho has an antihermitian part (hopefully
  // small). Then trace(rho+ rho) is second order in that error, while
  // trace(rho^2) is first order. With a bit more work, I could get
  // trace(rho^2) too, and then the trace of the hermitian part
  // squared is exactly .5*(Re[tr(rho^2)]+tr(rho+ rho))
  Contractor& contractor=Contractor::theContractor();
  complex_t trR2_c=contractor.contract(redRho,redRho);
  // check for real character?
  double trR2=real(trR2_c);
  if(abs(imag(trR2_c))>1E-10&&abs(imag(trR2_c))>1E-8*abs(trR2)){
    cout<<"WARNING (or error): the trace of rho^+ rho does not seem to be real!"<<endl;
    exit(1);
  }
  double normRho=real(conjugate(trR)*trR);
  return -log2(trR2/normRho); // log2(trace(rho^2))/(1-2)
  
  // For trace(rho^2) I construct an MPO that swaps up and doen
  // indices in the rho MPS, and compute the contraction of this in
  // between rho and the bra for its complex conjugate (to compensate
  // the conjugation that happens within contract)
  // static Operator swap;
  // static bool init=false;
  // if(!init){
  //   mwArray id2=identityMatrix(d*d);
  //   id2.reshape(Indices(d*d,d,d));
  //   id2.permute(Indices(1,3,2));
  //   id2.reshape(Indices(d*d,1,d*d,1));
  //   swap=Operator(id2);
  //   init=true;
  // }
  // int L=redRho.getL();
  // MPO aux(L);
  // for(int k=0;k<L;k++)
  //   aux.setOp(k,&swap,false);

}

void computeEntropies(const MPS& rho,complex_t trR,ofstream* out,int L0,bool header){
  *out<<setprecision(15);
  int L=rho.getLength();
  int pos1,pos2;
  // Now compute a series of predetermined partitions A-B(pos,L-pos-1)-C,
  // from pos=1 until below L/2. For each, A:BC and B:AC are computed,
  // which means I write Sa, Sbc, Sb, Sac
  for(int posC=1;posC<=L/2;posC++){
    if(!header){
      double Sa,Sb,Sbc,Sac;
      compute2RenyiEntropiesAB(rho,trR,0,posC-1,&Sa,&Sbc);
      compute2RenyiEntropiesAB(rho,trR,posC,L-posC-1,&Sb,&Sac);
      *out<<Sa<<"\t"<<Sbc<<"\t"<<Sb<<"\t"<<Sac<<"\t";
    }
    else{
      *out<<"S2[A(0-"<<posC-1<<")]\t S2[BC(L-A)] \t";
      *out<<"S2[B("<<posC<<"-"<<L-posC-1<<")]\t S2[AC] \t";
    }
  }
}

