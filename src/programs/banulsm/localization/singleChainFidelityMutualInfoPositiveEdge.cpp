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
    The leftmost L0 sites in the chain are used to encode a single qubit, so that 
    \f$| \tilde{0} \rangle \langle \tilde{0}|\f$ corresponds to 
    \f$\left ( | 0 \rangle \langle 0| \right )^{\otimes L_0}\f$
    tensored with identities on the remaining L-L0 sites.
    Starting from an MPO which encodes a pure state, Z+, Z-, X+ or X-,
    the program simulates the evolution and computes the trace norm of the resulting 
    simga_Z or sigma_X evolved state when tracing out the rightmost sites, 
    keeping \f$L_0,L0+2,\ldots 10\f$ spins. This will allow us to compute 
    the recovery fidelity.

    To be able to control the error in the trace, instead of evolving
    the sigmas, I evolve two pure states and correct the trace of each
    one after every step so that it is one. Then I compute the trace
    norm of the difference between the reduced density matrix for
    states + and -.

    Additionally, we compute the (2-Renyi) mutual information for the
    leftmost L0 to N-1 spins wrt the rest and across all cuts of the
    chain. The program can compute the average and localized
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
    \param <outfname_pol> (char*) name of the output file for the local polarization results
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

    In this version (positive) we use the trick that our starting
    states are \f$\rho=\rho_{1/2}^{\dagger} \rho_{1/2)\f$, where
    \f$\rho_{1/2}= |\phi_0\rangle \langle \phi_0| \otimes
    \left ( \frac{Id}{\sqrt{d}}\right )^{\otimes L-L_0}\f$ 
    (so, basically \f$\rho\f$ with different normalization, and 
    the evolved state can be written as 
    \f$rho(t)=\rho_{1/2}(t)^{\dagger} \rho_{1/2}(t)\f$,
    with \f$\rho_{1/2}(t)=U(t) \rho_{1/2} U(t)^{\dagger}\f$.  In
    practice, this means that squaring the obtained state (and
    properly normalizing) we obtain a positive MPO ansatz for the
    evolved state.
    To achieve the results with this construction, I use specific 
    (MPO-origami) tricks. Some optimization could be needed for the 
    expandOper step.
*/

/** Set initial states to initSt=|alpha + - > (alpha=x or z), encoded in L0 spins (from pos0) */
void setExtendedInitStates(MPS& rho0p,MPS& rho0m,int L,int d,int initSt,int L0,int pos0);

void computeLocalOp(const MPS& rho,int pos,int d,ofstream* out);

void computeTraceNormsPurification(const MPS& rho0p,const MPS& rho0m,const char* outfname,int L,int L0,bool isXY,double J_,double h_,int initSt,
		       double delta,int M,int D);
void constructLocal(int L,int pos,int initSt,int d,MPS& result);

/** Compute the trace norm of the difference of purfications after
    tracing out the specified sites. If gauge==true, it is assumed
    that both MPS have been properly ggauge conditioned, that the left
    and right parts give identity. */
double traceNormDiffPurification(const MPS& rhop,const MPS& rhom,const vector<int>& toTrace,bool gauge=false);

/** Compute the 2-Renyi entropy of the purification when the specified
    sites are traced out (if the vector of sites is empty, the whole
    state is used). */
double compute2RenyiSpurification(const MPS& rho,const vector<int>& toTrace);

/** Prepare the MPO that together with the state MPS performs the
    partial trace over the specified sites (the vector may be empty). */
void partialTraceMPO(const MPS& mps,const vector<int>& toTrace,MPO& mpo);

/** Compute the 2-Renyi entropies of two pieces of the chain, namely
    A, defined to be from pos1A until pos2A (both included and
    numbered from 0 to L) and the complementary, B */
void compute2RenyiEntropiesPurificationAB(const MPS& rho,int pos1A,int pos2A,double* SA,double*SB);
void computeEntropiesPurification(const MPS& rho,ofstream* out,int L0,bool header=false);

/** Same, but not using the purification structure */
void computeEntropies(const MPS& rho,complex_t trR,ofstream* out,int L0,bool header=false);
double compute2RenyiS(const MPS& redRho,complex_t traceR);
void compute2RenyiEntropiesAB(const MPS& rho,complex_t trR,int pos1A,int pos2A,double* SA,double*SB);

void normalizeTrace(MPS& rho,complex_t trR);

/** 
Generates a file name for the temporary MPS 
*/
const string mpsFileName(const string& baseDir,const string& baseName,int cnt=0);

/** Saves the state to a tmp file with new name and destroys the previous one */
void tmpSave(MPS& rho,const string& baseDir,const string& baseName,int cnt,string& tmpFile);

/** Check if there are tmp files with earlier stages of evolution */
int findTmpFiles(int M,const string& baseDir,const string& baseNameP,const string& baseNameM,string& mpsFileP,string& mpsFileM);

void computeLocalPolarization(const MPS& rho0p,const MPS& rho0m,const char* outfname,int L,int L0,bool isXY,
			      double J_,double h_,int initSt,double delta,int M,int D);


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
  string outfname_Sigz=props.getProperty("outputPolarization");
  string mpsdir=props.getProperty("outputDirMPS");
  string mpsfileP=props.getProperty("mpsFileP");
  string mpsfileM=props.getProperty("mpsFileM");
  int savingFreq=props.getIntProperty("savingFreq");
  _aux_=props.getIntProperty("newInstance");
  bool newInstance=_aux_>0;

  int pos0=0;    

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

  // Now read the parameters So that instances with larger L but same
  // number coincide in the left edge, I now fill from the left
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
      h[k]=isCont?h_*haux:(haux>0?h_:-h_);
    }
    inData.close();
  }

  cout<<"Read h="<<h<<endl;  
  if(!appendingFile){
    *out<<"% L="<<L<<", LEFT L0="<<L0<<", J="<<J_<<", h="<<h_
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
    computeEntropiesPurification(MPS(L,1,1),out,L0,true);
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
  cout<<"Created the Hamiltonian exponentials"<<endl;

  // Now construct the initial rho MPO, in this case, the identity,
  // except for L0 sites 
  MPS rho0p(L,1,d*d),rho0m(L,1,d*d);
  // Now read the initial rho MPO from the file or construct it for new instance
  if(!newInstance){
    rho0p.importMPS(tmpFileP.data());
    rho0m.importMPS(tmpFileM.data());
  }
  else{
    setExtendedInitStates(rho0p,rho0m,L,d,initSt,L0,pos0);
  }
  cout<<"Created the initial states"<<endl;

  //MPO auxRho(L);
  //MPOfromMPS(rho0,auxRho);
  //auxRho.exportForMatlab("rho0.m");
  //exit(1);

  //MPS hmps(L,1,d*d);
  MPO hamil(L);
  {
    const MPO& hamil_=hamH.getHMPO();
    extendMPO(hamil_,hamil,d);
  }
  //  MPSfromMPO(hamil,hmps);
  
  // Also the identity, to compute the trace
  MPS Id(L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++)
    Id.setA(k,id2);
  //  M1.exportForMatlab("M1.m");
  //Id.exportForMatlab("Id.m");
  
  // Now do the evolution step by step
  int cnt=M0;
  //  int posCenter=L/2; // for local opers
  double time=delta*cnt;
  // normalize purification:
  rho0p.gaugeCond('r',1);
  rho0m.gaugeCond('r',1);
  complex_t trRp=contractor.contract(rho0p,Id);
  complex_t trRm=contractor.contract(rho0m,Id);
  complex_t Eevp=contractor.contract(rho0p,hamil,rho0p);
  complex_t Eevm=contractor.contract(rho0m,hamil,rho0m);
  vector<int> toTrace; // empty
  // double S2Rp=compute2RenyiSpurification(rho0p,toTrace);
  // double S2Rm=compute2RenyiSpurification(rho0m,toTrace);
  double S2Rp=compute2RenyiS(rho0p,trRp);
  double S2Rm=compute2RenyiS(rho0m,trRm);
  //cout<<"Computed the intial S2 entropies of the whole + and -states"<<endl;
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<initSt<<"\t"
      <<real(Eevp)<<"\t"<<imag(Eevp)<<"\t";
  *out<<S2Rp<<"\t";
  computeEntropies(rho0p,trRp,out,L0);
  // computeEntropiesPurification(rho0p,out,L0);
  *out<<endl;
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<-initSt<<"\t"
      <<real(Eevm)<<"\t"<<imag(Eevm)<<"\t";
  *out<<S2Rm<<"\t";
  computeEntropies(rho0p,trRp,out,L0);
  //  computeEntropiesPurification(rho0m,out,L0);
  *out<<endl;
  out->close();delete out;

  computeTraceNormsPurification(rho0p,rho0m,outfname_Traces.data(),L,L0,isXY,J_,h_,initSt,delta,cnt,D);
  computeLocalPolarization(rho0p,rho0m,outfname_Sigz.data(),L,L0,isXY,J_,h_,initSt,delta,cnt,D);

  //cout<<"Computed the initial trace norms "<<endl;
  while(cnt<Mtot){
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
	//cout<<"Time evolution step nr "<<cnt<<", time="<<time<<" for state "<<nrSt<<endl;
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
      rho0.gaugeCond('R',1);
      // now renormalize them
      if(nrSt==0) cnt=cntL; // reset the counter for the second state
    }
    Eevp=contractor.contract(rho0p,hamil,rho0p);
    Eevm=contractor.contract(rho0m,hamil,rho0m);
    out=new ofstream(outfname.data(),ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    // For the most costly Renyi entropies I use directly rho (square
    //of the one before) which should only differ from the computed
    //one in the normalization=> I compute the trace to keep track of
    //that.
    complex_t trRp=contractor.contract(rho0p,Id); 
    S2Rp=compute2RenyiS(rho0p,trRp);
    // S2Rp=compute2RenyiSpurification(rho0p,toTrace);
    //S2Rm=compute2RenyiSpurification(rho0m,toTrace);
    *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<initSt<<"\t"
	<<real(Eevp)<<"\t"<<imag(Eevp)<<"\t";
    *out<<S2Rp<<"\t";
    computeEntropies(rho0p,trRp,out,L0);
    //computeEntropiesPurification(rho0p,out,L0);
    *out<<endl;
    complex_t trRm=contractor.contract(rho0m,Id); 
    S2Rm=compute2RenyiS(rho0m,trRm);
    *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<-initSt<<"\t"
	<<real(Eevm)<<"\t"<<imag(Eevm)<<"\t";
    *out<<S2Rm<<"\t";
    computeEntropies(rho0m,trRm,out,L0);
    //computeEntropiesPurification(rho0m,out,L0);
    *out<<endl;
    out->close();delete out;
    //cout<<"About to compute trace norms for time evolution step nr "<<cnt<<", time="<<time<<endl;
    computeTraceNormsPurification(rho0p,rho0m,outfname_Traces.data(),L,L0,isXY,J_,h_,initSt,delta,cnt,D);
    //cout<<"Computed the trace norms "<<endl;
    computeLocalPolarization(rho0p,rho0m,outfname_Sigz.data(),L,L0,isXY,J_,h_,initSt,delta,cnt,D);
   
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


// void MPSfromMPO(const MPO& mpo,MPS& mps,bool up){
//   int L=mpo.getLength();
//   // vector<int> du=mpo.getDimensions();
//   // vector<int> dd=mpo.getOriginalDimensions();
//   // int D=mpo.getOp(0).getDr();

//   // vector<int> newDims(du);
//   // for(int k=0;k<L;k++) newDims[k]*=dd[k];

//   mps=MPS(L,1,1);

//   for(int k=0;k<L;k++){
//     mwArray aux=mpo.getOp(k).getFullData();
//     Indices dims=aux.getDimensions();
//     if(up)
//       aux.permute(Indices(1,3,2,4));
//     else
//       aux.permute(Indices(3,1,2,4));
//     aux.reshape(Indices(dims[0]*dims[2],dims[1],dims[3]));
//     mps.replaceSite(k,aux,0); // brute force
//   }

// }


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


void computeTraceNormsPurification(const MPS& rho0p,const MPS& rho0m,const char* outfname,int L,int L0,bool isXY,double J_,double h_,int initSt,
		       double delta,int M,int D){
  //cout<<"ComputeTraceNormsPurification "<<endl;
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

  MPS rhop(rho0p); // copying
  MPS rhom(rho0m); // copying
  //cout<<"Copied states"<<endl;
  // First I impose gauge condition on the right part (outside L0)
  for(int k=L-1;k>=L0-1;k--){
    rhop.gaugeCond(k,'L',true);
    rhom.gaugeCond(k,'L',true);
  }
  //cout<<"Imposed gauge conditions"<<endl;
  // Now I know the result of tracing out any part of the edges
  *out<<M<<"\t"<<delta<<"\t"<<M*delta<<"\t"<<D<<"\t"; // 4 columns for nothing (3 is time)

  // Now I need to start tracing out things, to construct reduced operators and to compute their trace norms.
  // I will keep 10, 9, 8, 7, 6, ... L0+1, L0 on the left
  // Column 4+2 is 10; 8 is 9; 10 is 8; 12 is 7; 14 is 6; ...

  for(int firstTracedR=min(10,L-1);firstTracedR>=L0;firstTracedR--){
    vector<int> toTrace;
    for(int k=firstTracedR;k<L;k++){
      toTrace.push_back(k);
    }
    cout<<"To get the leftmost "<<firstTracedR<<" spins (L="<<L<<"), tracing out "<<toTrace<<endl;
    double tN10=traceNormDiffPurification(rhop,rhom,toTrace,true);
    //cout<<"trace norm="<<tN10<<endl;
    // And finally write the results, preceeded by the number kept!
    *out<<firstTracedR<<"\t";
    *out<<tN10<<"\t";
  }
  //cout<<"Computation finished"<<endl;
  *out<<endl;
  out->close();
  delete out;
}


bool containsSite(int pos,const vector<int>& toTrace){
  bool found=0;
  vector<int>::const_iterator it=toTrace.begin();
  while(!found&&it!=toTrace.end()){
    if(*it==pos){
      found=true;
    }
    it++;
  }
  return found;
}

void complementaryMPO(const MPS& rho,const vector<int>&toTrace,MPO& mpo,bool dagger=false){
  int L=rho.getLength();
  vector<int> toKeep;
  for(int k=0;k<L;k++) if(!containsSite(k,toTrace)) toKeep.push_back(k);
  //  cout<<"complementaryMPO of rho="<<rho<<", toTrace="<<toTrace<<", toKeep="<<toKeep<<endl;
  mpo.initLength(toKeep.size()+1);
  int Dr=rho.getA(toKeep.back()).getDr();
  mwArray IdDr=identityMatrix(Dl);
  if(!dagger)
    IdDr.reshape(Indices(Dr,Dr,1,1));
  else
    IdDr.reshape(Indices(1,Dr,Dr,1));
  for(int k=0;k<toKeep.size();k++){
    mwArray aux=rho.getA(toKeep[k]).getA();
    Indices dims=aux.getDimensions();
    int d0=sqrt(dims[0]);
    if(d0*d0!=dims[0]){
      cout<<"Error: Dimension of site "<<k<<" ("<<dims[0]
	  <<") does not seem a square=> don't know how t divide it"<<endl;
      exit(1);
    }
    aux.reshape(Indices(d0,d0,dims[1],dims[2]));
    if(!dagger)
      aux.permute(Indices(1,3,2,4),false);
    else
      aux.permute(Indices(2,3,1,4),true);
    mpo.setOp(k,new Operator(aux),true);
  }
  mpo.setOp(toKeep.size(),new Operator(IdDr),true);
}

void expandDoubleOper(const MPO& rhop,mwArray& oper){
  //cout<<"expandDoubleOper "<<rhop<<endl;
  int L=rhop.getLength();
  oper=mwArray(ONE_c);
  mwArray operR=oper; // also start from the right
  int Lmid=floor(L/2); // where to break into left and right parts

  int dphys1=1; // temporary output dimensions 
  int dphys2=1; // temporary input (right) dimensions 
  int Do=1; // open bond dimension
  int Dor=1;
  int dphys1r=1; // temporary output dimensions 
  int dphys2r=1; // temporary input (right) dimensions 
  oper.reshape(Indices(dphys1*dphys2,Do*Do));
  operR.reshape(Indices(Dor*Dor,dphys1r*dphys2r));      
  for(int k=0;k<=Lmid;k++){
    oper.reshape(Indices(dphys1*dphys2*Do,Do));
    mwArray aux=rhop.getOp(k).getFullData();
    //cout<<"Now oper is "<<oper.getDimensions()<<". Multiplying oper{"<<k<<"}="<<aux.getDimensions()<<endl;
    Indices dims=aux.getDimensions();
    int d1=dims[0];int Dl=dims[1];int d2=dims[2];int Dr=dims[3];
    if(Dl!=Do){
      cout<<"Error: Dimensions do not agree in expandOper!"<<endl;
      exit(1);
    }
    aux.permute(Indices(2,1,3,4));
    aux.reshape(Indices(Dl,d1*d2*Dr));
    oper.multiplyRight(aux);
    oper.reshape(Indices(dphys1*dphys2,Do*d1,d2*Dr));
    oper.permute(Indices(1,3,2));oper.reshape(Indices(dphys1*dphys2*d2*Dr,Do*d1));
    aux=rhop.getOp(k).getFullData(); // recover the same operator for rho^dagger
    aux.permute(Indices(2,1,3,4),true);
    aux.reshape(Indices(Dl*d1,d2*Dr));
    oper.multiplyRight(aux);
    oper.reshape(Indices(dphys1,dphys2,d2,Dr,d2,Dr));
    oper.permute(Indices(1,5,2,3,6,4));
    dphys1=dphys1*d2;
    dphys2=dphys2*d2;
    Do=Dr;
    oper.reshape(Indices(dphys1*dphys2,Do*Do));
    //cout<<" oper is "<<oper.getDimensions()<<endl;
  }
  // cout<<"After multiplying the "<<Lmid+1<<" left-most sites, dphys1="<<dphys1
  //     <<", dphys2="<<dphys2<<", Do="<<Do<<", oper is "<<oper.getDimensions()<<endl;
  // Now, if Lmid is not the end, start from right.
  if(Lmid<(L-1)){
    for(int k=L-1;k>Lmid;k--){
      operR.reshape(Indices(Dor,Dor*dphys1r*dphys2r));
      mwArray aux=rhop.getOp(k).getFullData();
      //cout<<"Now operR is "<<operR.getDimensions()<<". Multiplying oper{"<<k<<"}="<<aux.getDimensions()<<endl;
      Indices dims=aux.getDimensions();
      int d1=dims[0];int Dl=dims[1];int d2=dims[2];int Dr=dims[3];
      if(Dr!=Dor){
	cout<<"Error: Dimensions do not agree in expandOper!"<<endl;
	exit(1);
      }
      //aux.permute(Indices(1,2,3,4),true); 
      aux.conjugate(); /// rho dagger
      aux.reshape(Indices(d1*Dl*d2,Dr));
      operR.multiplyLeft(aux);
      //cout<<"Multiplication done: operR is "<<operR.getDimensions()<<endl;
      operR.reshape(Indices(d1,Dl,d2,Dor,dphys1r*dphys2r));
      //cout<<"Reshaped to d1("<<d1<<"),Dl("<<Dl<<"),d2("<<d2<<"),Dor("<<Dor<<"),dphys1r("<<dphys1r<<")*dphys2r("<<dphys2r<<") operR is "<<operR.getDimensions()<<endl;
      operR.permute(Indices(1,4,2,3,5));operR.reshape(Indices(d1*Dor,Dl*d2*dphys1r*dphys2r));
      aux=rhop.getOp(k).getFullData(); // recover the same operator for rho
      //cout<<"Multiplying now by aux "<<aux.getDimensions()<<endl;
      aux.permute(Indices(2,3,1,4));
      aux.reshape(Indices(Dl*d2,d1*Dr));
      operR.multiplyLeft(aux);
      //cout<<"Multiplication done: operR is "<<operR.getDimensions()<<endl;
      operR.reshape(Indices(Dl,d2,Dl,d2,dphys1r,dphys2r));
      //cout<<"Reshaped to Dl("<<Dl<<"),d2("<<d2<<"),Dl,d2,dphys1r("<<dphys1r<<"),dphys2r("<<dphys2r<<"), operR is "<<operR.getDimensions()<<endl;
      operR.permute(Indices(3,1,4,5,2,6));
      dphys1r=dphys1r*d2;
      dphys2r=dphys2r*d2;
      Dor=Dl;
      //cout<<"Reshaping to new Dor("<<Dor<<")*Dor,dphys1r("<<dphys1r<<")*dphys2r("<<dphys2r<<")"<<endl;
      operR.reshape(Indices(Dor*Dor,dphys1r*dphys2r));      
    }
  }
  // Now contract both
  if(Do!=Dor){
    cout<<"Error: Dimensions do not agree in expandDoubleOper!: Do="<<Do<<", Dor="<<Dor<<endl;
    exit(1);
  }
  oper.multiplyRight(operR);
  oper.reshape(Indices(dphys1,dphys2,dphys1r,dphys2r));
  oper.permute(Indices(1,3,2,4));
  oper.reshape(Indices(dphys1*dphys1r,dphys2*dphys2r));
}


double traceNormDiffPurification(const MPS& rhop,const MPS& rhom,const vector<int>& toTrace,bool gauge){
  if(!gauge){
    cout<<"ERROR: The proper gauge condition has to be applied before calling traceNormDiffPurification!"<<endl;
    exit(1);
    // Actually, I should add this case
  }

  //  cout<<"traceNormDiffpurification of rhop="<<rhop<<", toTrace="<<toTrace<<endl;

  mwArray OpP,OpM;

  // Construct the basic MPO for the complementary to toTrace
  int L=rhop.getLength();
  MPO rhoMPO(L-toTrace.size());
  //MPO rhoMPOd(L-toTrace.size());
  complementaryMPO(rhop,toTrace,rhoMPO);
  //  cout<<"complementaryMPO of rhoP found: "<<rhoMPO<<endl; 
  expandDoubleOper(rhoMPO,OpP);
  // expandOper(rhoMPO,OpP);
  // cout<<"Operator OpP expanded to dims "<<OpP.getDimensions()<<endl;
  // OpP.multiplyLeft(Hconjugate(OpP));
  complementaryMPO(rhom,toTrace,rhoMPO);
  // cout<<"complementaryMPO of rhoP found: "<<rhoMPO<<endl; 
  expandDoubleOper(rhoMPO,OpM);
  // expandOper(rhoMPO,OpM);
  // cout<<"Operator OpM expanded to dims "<<OpM.getDimensions()<<endl;
  // OpM.multiplyLeft(Hconjugate(OpM));


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
  //if(trR!=ONE_c) normalizeTrace(rhoA,trR);
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
  //if(trR!=ONE_c) normalizeTrace(rhoA,trR);
  rhoA.traceOutSites(toTrace);
  // Now compute its 2-Renyi entropy
  *SB=compute2RenyiS(rhoA,trR);  
}

/** Compute the 2-Renyi entropies of two pieces of the chain, namely
    A, defined to be from pos1A until pos2A (both included and
    numbered from 0 to L-1) and the complementary, B */
void compute2RenyiEntropiesPurificationAB(const MPS& rho,int pos1A,int pos2A,double* SA,double*SB){
  int L=rho.getLength();
  // First trace out the complementary of A
  vector<int> toTrace;
  for(int k=0;k<L;k++){
    if(k<pos1A||k>pos2A)
      toTrace.push_back(k);
  }
  // Now compute its 2-Renyi entropy
  *SA=compute2RenyiSpurification(rho,toTrace);
  // And now the same for the complementary
  toTrace.clear();
  for(int k=pos1A;k<=pos2A;k++)
    toTrace.push_back(k);
  // Now compute its 2-Renyi entropy
  *SB=compute2RenyiSpurification(rho,toTrace);  
}

void partialTraceMPO(const MPS& rho,const vector<int>& toTrace,MPO& mpo){
  Contractor& contractor=Contractor::theContractor();
  int L=rho.getLength();
  mpo.initLength(L);
  // the tensors as MPO, for the adjoint of rho
  MPOfromMPS(rho,mpo,false,true);
  // now reshape the ones to be traced, and include the extra identity for the ones that are not traced
  for(int k=0;k<L;k++){
    mwArray A=mpo.getOp(k).getFullData(); 
    Indices dimsA=A.getDimensions(); // du,Dl,dd,Dr
    // is it traced?
    if(containsSite(k,toTrace)){
      A.permute(Indices(2,3,1,4)); // Dl,dd,du,Dr
      A.reshape(Indices(1,dimsA[1],dimsA[2]*dimsA[0],dimsA[3])); // 1,Dl,dd*du,Dr
      mpo.setOp(k,new Operator(A),true);
    }
    else{ // not traced => add an auxiliary identity for the ancilla
      // I could check that du and dd are the same, but it really does not matter
      // I add the identity that leaves the ancilla (index going up in MPO) unchanged
      mwArray Id=identityMatrix(dimsA[0]);Id.reshape(Indices(dimsA[0],1,dimsA[0],1));
      mpo.setOp(k,new DoubleOperator(A,Id),true);
    }
  }
}

double compute2RenyiSpurification(const MPS& rho,const vector<int>& toTrace){
  //  cout<<"Compute2RenyiSpurification of rho="<<rho<<", toTrace="<<toTrace<<endl;
  Contractor& contractor=Contractor::theContractor();
  if(abs(contractor.contract(rho,rho)-ONE_c)>1E-13){
    cout<<"ERROR: Not normalized state for compute2RenyiSpurification!!"<<endl;
    exit(1);
  }
  int L=rho.getLength();
  // I construct an auxiliary MPO with the rho^dagger, such that the
  // trace is over toTrace sites performed upon contraction, and the
  // other sites contain an extra identity on the ancilla
  MPO rhoDagtr(L);
  partialTraceMPO(rho,toTrace,rhoDagtr);
  //  cout<<"Created an MPO for the adjoint and partial trace over "<<toTrace<<", namely "<<rhoDagtr<<endl;
  // Now I contract this with the adjoint
  complex_t trR2_c=contractor.contract2(rhoDagtr,rho);
  // check for real character?
  double trR2=real(trR2_c);
  if(abs(imag(trR2_c))>1E-10&&abs(imag(trR2_c))>1E-8*abs(trR2)){
    cout<<"WARNING (or error): the trace of rho^+ rho does not seem to be real!"<<endl;
    exit(1);
  }
  return -log2(trR2); // log2(trace(rho^2))/(1-2)
}

double compute2RenyiS(const MPS& redRho,complex_t trR){
  // I need trace(rho^2)
  
  // If rho was Hermitian, as it ideally should,
  // trace(rho^2)=trace(rho+ rho) We will have probably some
  // numerical error, so that rho has an antihermitian part (hopefully
  // small). Then trace(rho+ rho) is second order in that error, while
  // trace(rho^2) is first order. With a bit more work, I could get
  // trace(rho^2) too, and then the trace of the hermitian part
  // squared is exactly .5*(Re[tr(rho^2)]+tr(rho+ rho))
  if(redRho.getLength()==0) return 0.;
  Contractor& contractor=Contractor::theContractor();
  complex_t trR2_c=contractor.contract(redRho,redRho);
  // check for real character?
  double trR2=real(trR2_c);
  if(abs(imag(trR2_c))>1E-10&&abs(imag(trR2_c))>1E-8*abs(trR2)){
    cout<<"WARNING (or error): the trace of rho^+ rho does not seem to be real!"<<endl;
    exit(1);
  }
  double normRho=real(conjugate(trR)*trR);
  //   cout<<"Computing 2Renyi of state with length "<<redRho.getLength()<<", results of tr(rho2)="
  //       <<trR2_c<<", normalization="<<normRho<<", RESULT="<<-log2(trR2/normRho)<<endl;
  return -log2(trR2/normRho); // log2(trace(rho^2))/(1-2)
}

void computeEntropiesPurification(const MPS& rho,ofstream* out,int L0,bool header){
  *out<<setprecision(15);
  int L=rho.getLength();
  int pos1,pos2;
  // Now compute a series of predetermined partitions A(0:pos-1)-B(pos,L-1),
  // from pos=1 until L-2. For each, A:B I write Sa and Sb
  for(int posC=1;posC<=L-2;posC++){
    if(!header){
      double Sa,Sb;
      compute2RenyiEntropiesPurificationAB(rho,0,posC-1,&Sa,&Sb);
      *out<<Sa<<"\t"<<Sb<<"\t";
    }
    else{
      *out<<"S2[A(0-"<<posC-1<<")]\t S2[B(L-A)] \t";
    }
  }
}

void computeEntropies(const MPS& rho,complex_t trR,ofstream* out,int L0,bool header){
  *out<<setprecision(15);
  int L=rho.getLength();
  int pos1,pos2;
  // Now compute a series of predetermined partitions A(0:pos-1)-B(pos,L-1),
  // from pos=1 until L-2. For each, A:B I write Sa and Sb
  for(int posC=1;posC<=L-2;posC++){
    if(!header){
      double Sa,Sb;
      compute2RenyiEntropiesAB(rho,trR,0,posC-1,&Sa,&Sb);
      *out<<Sa<<"\t"<<Sb<<"\t";
    }
    else{
      *out<<"S2[A(0-"<<posC-1<<")]\t S2[BC(L-A)] \t";
    }
  }
}

void computeLocalPolarization(const MPS& rho0p,const MPS& rho0m,const char* outfname,int L,int L0,bool isXY,double J_,double h_,int initSt,
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
    *out<<"% M\t delta\t t\t D\t <sigZ(l)>"<<endl;
  }
  *out<<setprecision(15);

  Contractor& contractor=Contractor::theContractor();

  // I need to prepare the MPo/MPS to compute the observable
  // Prepare the observables to compute
  // The identity, to compute the trace
  static MPS Id=MPS(L,1,d*d);
  // And I will also need the sigmaz operator on each site
  static MPS SzLocal=Id;
  static mwArray sigZ;
  static bool initialized=false;
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  if(!initialized){
    for(int k=0;k<L;k++)
      Id.setA(k,id2);
    SzLocal=Id;
    complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
    sigZ=mwArray(Indices(d,d),dataZ);
    sigZ.reshape(Indices(d*d,1,1));
  }

  // Site by site, compute values and store in vectors
  vector<double> sigZplus,sigZminus;
  for(int k=0;k<L;k++){
    SzLocal.setA(k,sigZ);
    // Compute and store both values
    complex_t val=contractor.contract(rho0p,SzLocal);
    complex_t trace=contractor.contract(rho0p,Id);
    sigZplus.push_back(real(val/trace));
    val=contractor.contract(rho0m,SzLocal);
    trace=contractor.contract(rho0m,Id);
    sigZminus.push_back(real(val/trace));
    // set the identity back in place
    SzLocal.setA(k,id2);
  }
  // Write numbers to the file  
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<initSt<<"\t"; // 4 columns for nothing (3 is time)
  for(int k=0;k<L;k++)
    *out<<sigZplus[k]<<"\t";
  *out<<endl;
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<-initSt<<"\t"; // 4 columns for nothing (3 is time)
  for(int k=0;k<L;k++)
    *out<<sigZminus[k]<<"\t";
  *out<<endl;
  out->close();
  delete out;
}
