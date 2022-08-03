
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"
#include "HeisenbergHamiltonian.h"

using namespace shrt;

void sumMPO(MPO& result,const MPO& mpo1,const MPO& mpo2);

void computeLocalEnergy(const MPS& mps,double J_,double g_,double h_,vector<double>& locE);

int computeBlockEntropies(MPS& mps,vector<double>& trDistRDM,
			   vector<double>& entropyRDM,vector<int>& sizeRDM);

int computeBlockEntropiesSlow(MPS& mps,vector<double>& trDistRDM,
			   vector<double>& entropyRDM,vector<int>& sizeRDM);

int computeBlockEntropiesExact(MPS& mps,vector<double>& trDistRDM,
			   vector<double>& entropyRDM,vector<int>& sizeRDM);

void computeTwoBodyTerm(mwArray& result,double Jx,double Jy,double Jz,
			double hx1,double hx2,double hy1,double hy2,double hz1,double hz2);

void splitTwoBodyTerms(const mwArray& h12,mwArray& Ol,mwArray& Or,mwArray& Omid12,mwArray& Omid21);

void initialState(MPS& mps,int intSt,int D);

void computeEnergyEnergy(const MPS& mps,mwArray& h12,vector<int>& pos,vector<double>& h0hm,vector<double>& hm,int offset=0);

/** Construct a file name to save (tmp) results */
const string mpsFileName(const string& mpsdir,const string& baseName,const string& label,int cnt);

/** Save the tmp files */
void tmpSave(MPS& mps,double norm_MPS,MPS& chebyT_Nm1,double norm_Nm1,MPS& chebyT_N,
	     double norm_N,const string& baseDir,const string& baseName,int cnt,vector<string>& tmpFiles);

/** Find tmp files and restore the status */
int tmpRestore(MPS& mps,double& norm_MPS,MPS& chebyT_Nm1,double& norm_Nm1,MPS& chebyT_N,
	       double& norm_N,const string& baseDir,const string& baseName,int M,vector<string>& tmpFiles);


/** analyzeMinVarChebyDelta reads the final states produced by minVarChebyDelta 
    with the same input parameters and (re)computes entropies of middle blocks, 
    and energy-density correlations.
 */

int d=2;

int main(int argc,const char* argv[]){
  // Read input arguments
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

  int L=props.getIntProperty("L");
  double J_=props.getDoubleProperty("J");
  double h_=props.getDoubleProperty("h");
  double g_=props.getDoubleProperty("g");
  bool ising=(props.getIntProperty("ising")!=0); // only if -ising=0 is given it will be Heisenberg (or XY)  
  // These are only used if Heisenberg model
  double Jx_ = props.getDoubleProperty("Jx");  
  double Jy_ = props.getDoubleProperty("Jy");  
  double Jz_ = props.getDoubleProperty("Jz");
  // It can also include x and y local fields Heisenberg model
  double hx_ = props.getDoubleProperty("hx");  
  double hy_ = props.getDoubleProperty("hy");  
  double delta=props.getDoubleProperty("delta");
  if(delta<0) delta=0.02; 
  int D = props.getIntProperty("D");  
  //string outfname=props.getProperty("output");
  //string Sfname=props.getProperty("outputS");
  string Sredfname=props.getProperty("outputSblock"); // for the entropy of the middle blocks
  string corrfname=props.getProperty("outputCorr"); // for the energy correlations
  string mpsdir=props.getProperty("mpsdir");
  string mpsfileroot=props.getProperty("mpsfileroot");
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  int M=props.getIntProperty("M"); // LARGEST ORDER: which result to look for
  int initSt=props.getIntProperty("initSt"); // initial state to use:
  //                               +1=X+; -1=X-; (2,3 for Y, Z)
  //                               +4,+5,+6 are the staggered(+-) for X, Y and Z and -4,-5,-6 the -+
  //                               +7 is a (complex) random state and +8 a real random (in which cases, the random seed may be given)
  //                               +10 is the staggered with period 4 (++--) in the Z direction, and -10 is (--++)
  //                               +11 is a step configuration (valid right now only for Ising)
  //                                   then requires Thleft and Thright to fix the excess of energy left-right
  //                                   and also a file to save the distribution of energy density

  //  bool checkRDM=(props.getIntProperty("checkRDM")>0); // whether to check how close RDM are to thermal
  // if(rate<0) rate=1;

  int lastSavedCnt=-1;
  //  int Lb=props.getIntProperty("Lb");
  //if(Lb<0) Lb=6;
  bool avrg=1;
  avrg=!(props.getIntProperty("average")==0); // whether to average for each Lc subsystem

  bool isXY=(!ising&&Jz_==0);
  cout<<"Initialized arguments: L="<<L;
  if(ising){
    cout<<"; ISING: J="<<J_
	<<", g="<<g_
	<<", h="<<h_;
  }
  else{
    if(!isXY){
      cout<<"; HEISENBERG: Jx="<<Jx_
	  <<", Jy="<<Jy_
	  <<", Jz="<<Jz_
	  <<", h="<<h_;
    }
    else{ // XY case
      cout<<"; XY: Jx="<<Jx_
	  <<", Jy="<<Jy_
	  <<", hx="<<hx_
      	  <<", hy="<<hy_
	  <<", hz="<<h_;
    }    
  }
  cout<<", app="<<app
      <<", D="<<D<<endl;

  if((!file_exists(Sredfname.data()))||(!file_exists(corrfname.data()))) app=0; // also if it does not exist (new)
  ofstream* outSred;
  ofstream* outEcorr;
  if(!app){
    outSred=new ofstream(Sredfname.data());
    if(!outSred->is_open()){
      cout<<"Error: impossible to open file "<<Sredfname<<
	" for output"<<endl;
      exit(1);
    }
    if(ising){
      *outSred<<"% N\t J\t g\t h\t D\t Dcut\t cnt\t Lc\t Sred(Lc)\t trDist"<<endl;
    }
    else{
      if(!isXY)
	*outSred<<"% N\t Jx\t Jy\t Jz\t h\t D\t cnt\t Lc\t Sred\t trDist"<<endl;
      else
	*outSred<<"% N\t Jx\t Jy\t hx\t hy\t hz\t D\t cnt\t Lc\t Sred(Lc)\t trDist"<<endl;
    }
    outSred->close();
    delete outSred;

    outEcorr=new ofstream(corrfname.data());
    if(!outEcorr->is_open()){
      cout<<"Error: impossible to open file "<<corrfname<<
	" for output of energy correls"<<endl;
      exit(1);
    }
    if(ising){
      *outEcorr<<"% N\t J\t g\t h\t D\t cnt\t E\t <H^2>\t x\t <h0 h(x)>\t <h(x)>"<<endl;
    }
    else{
      if(!isXY)
	*outEcorr<<"% N\t Jx\t Jy\t Jz\t h\t D\t cnt\t E\t <H^2>\t x\t <h0 h(x)>\t <h(x)>"<<endl;
      else
	*outEcorr<<"% N\t Jx\t Jy\t hx\t hy\t hz\t D\t cnt\t E\t <H^2>\t x\t <h0 h(x)>\t <h(x)>"<<endl;
    }
    outEcorr->close();
    delete outEcorr;        
  }
  
  
  cout<<setprecision(10);

  
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;

  MPS mps(L,D,d);

  // cout<<"Initialized state ("<<initSt<<") from scratch, norm "<<contractor.contract(mps,mps)<<endl;
  
  
  // The original Hamiltonian
  // First: H
  //  IsingHamiltonian hamI(L,d,J,g,h); 
  Hamiltonian* ham0;
  mwArray h12,hlast;
  if(ising){
    ham0=new IsingHamiltonian(L,d,J_,g_,h_);
    //computeTwoBodyTerm(h12,0.,0.,J_,.5*g_,.5*g_,0.,0.,.5*h_,.5*h_);
    computeTwoBodyTerm(h12,0.,0.,J_,g_,0.,0.,0.,h_,0.);
  }
  else{
    ham0=new HeisenbergHamiltonian(L,Jx_,Jy_,Jz_,hx_,hy_,h_);
    //computeTwoBodyTerm(h12,.25*Jx_,.25*Jy_,.25*Jz_,.25*hx_,.25*hx_,.25*hy_,.25*hy_,.25*h_,.25*h_);
    computeTwoBodyTerm(h12,.25*Jx_,.25*Jy_,.25*Jz_,.5*hx_,0.,.5*hy_,0.,.5*h_,0.);	
  }
  //IsingHamiltonian hamH0(L,d,J_,g_,h_);
  const MPO& hamil0=ham0->getHMPO();

  //cout<<"The two body term is h12="<<h12<<endl;
  
  // These are just placeholders for the tmp data
  MPS* chebyT_Nm1=new MPS(mps); // poly(-1)=0
  // Set the first tensor to 0=> ignore norm factor in this case
  // chebyT_Nm1->gaugeCond('R',0);chebyT_Nm1->gaugeCond('L',1);
  // chebyT_Nm1->setA(0,ZERO_c*chebyT_Nm1->getA(0).getA());
  double norm_Nm1=0.;
  MPS* chebyT_N=new MPS(mps); // normalized already
  double norm_N=0.;

  double normMPS=0;
  
  complex_t valsE,valsE2;
  double Shalf;

  vector<string> tmpFiles;
  int cnt=0;
  if(M==0){ // Special case, just to compute initial values (mostly a check, as I can compute them exactly)
    mps=MPS(L,1,d);
    initialState(mps,initSt,1);
    mps.gaugeCond('R',1);mps.gaugeCond('L',1); // not needed
  }
  else{
    cnt=tmpRestore(mps,normMPS,*chebyT_Nm1,norm_Nm1,*chebyT_N,norm_N,mpsdir,mpsfileroot,M,tmpFiles);

    if(cnt==0){
      cout<<"ERROR: No stored data found for this case!!!!"<<endl;
      exit(1);
    }
    cout<<"Just read tmp files for step "<<cnt<<", normMPS="<<normMPS<<", norm_N="<<norm_N<<", norm_Nm1="<<norm_Nm1<<endl;
    
    if(cnt<M-1){
      cout<<"Data for proper M is not stored. Please, run chebH2 program again"<<endl;
      exit(1);
    }
  }
  
  // Now compute all values 
  valsE=contractor.contract(mps,hamil0,mps); // <H>
  valsE2=contractor.contract2(hamil0,mps); // <H^2>
  Shalf=contractor.getEntropy(mps); // S(L/2)

  // entropies and distances to thermal of central blocks 
  
  vector<double> trDistRDM;
  vector<double> entropyRDM;
  vector<int> sizeRDM;
  
  // Also energy-energy correlations
  vector<double> h0hm(1,0.);vector<double> hm(1,0.);
  vector<int> posM(1,0);


  int offset=0;
  if(!ising&&abs(initSt)==10){
    if(L/2%4!=0)offset=2;
  }   
  computeEnergyEnergy(mps,h12,posM,h0hm,hm,offset);
  // Write in file
  outEcorr=new ofstream(corrfname.data(),ios::app);
  if(!outEcorr->is_open()){
    cout<<"Error: impossible to open file "<<corrfname<<
      " for output"<<endl;
    exit(1);
  }
  *outEcorr<<setprecision(15);
  for(int k=0;k<posM.size();k++){
    // OJO!!!!! This is a mistake, doing it wrong for Heisenberg!!!!
    if(ising)
      *outEcorr<<L<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t";
    else
      if(!isXY)
	*outEcorr<<L<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t";
      else
	*outEcorr<<L<<"\t"<<Jx_<<"\t"<<Jy_<<"\t"<<hx_<<"\t"<<hy_<<"\t"<<h_<<"\t"<<D<<"\t";
    *outEcorr<<cnt<<"\t"<<real(valsE)<<"\t"<<real(valsE2)<<"\t";
    *outEcorr<<posM[k]<<"\t"<<h0hm[k]<<"\t"<<hm[k]<<"\t";
    *outEcorr<<endl;
  }
  outEcorr->close(); delete outEcorr;

  //  int Dcut=computeBlockEntropiesSlow(mps,trDistRDM,entropyRDM,sizeRDM);
  int Dcut=D;computeBlockEntropiesExact(mps,trDistRDM,entropyRDM,sizeRDM);
  // Write in file
  outSred=new ofstream(Sredfname.data(),ios::app);
  if(!outSred->is_open()){
    cout<<"Error: impossible to open file "<<Sredfname<<
      " for output"<<endl;
    exit(1);
  }
  *outSred<<setprecision(15);
  for(int k=0;k<sizeRDM.size();k++){
    // OJO!!!!! This is a mistake, doing it wrong for Heisenberg!!!!
    if(ising)
      *outSred<<L<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t"<<Dcut<<"\t";
    else
      if(!isXY)
	*outSred<<L<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t"<<Dcut<<"\t";
      else
	*outSred<<L<<"\t"<<Jx_<<"\t"<<Jy_<<"\t"<<hx_<<"\t"<<hy_<<"\t"<<h_<<"\t"<<D<<"\t"<<Dcut<<"\t";
    *outSred<<cnt<<"\t"<<sizeRDM[k]<<"\t"<<entropyRDM[k]<<"\t"<<trDistRDM[k];
    *outSred<<endl;
  }
  outSred->close(); delete outSred;
  
}


void computeLocalEnergy(const MPS& mps,double J_,double g_,double h_,vector<double>& locE){
  locE.clear();
  int L=mps.getLength();
  Contractor& contractor=Contractor::theContractor();
  IsingHamiltonian hamLoc(2,d,J_,.5*g_,.5*h_);
  const MPO& mpo=hamLoc.getHMPO();
  MPO auxMPO(L);
  Operator idOp(reshape(identityMatrix(d),Indices(d,1,d,1)));
  for(int k=0;k<L;k++)
      auxMPO.setOp(k,&idOp,false);
  complex_t normMPS=contractor.contract(mps,mps);
  // the rest identity
  for(int k=0;k<L-1;k++){
    auxMPO.setOp(k,&mpo.getOp(0),false);
    auxMPO.setOp(k+1,&mpo.getOp(1),false);
    complex_t valE=contractor.contract(mps,auxMPO,mps);
    auxMPO.setOp(k,&idOp,false);
    auxMPO.setOp(k+1,&idOp,false);
    locE.push_back(real(valE/normMPS));
  }

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
  // And for the staggered with period 2 I apply the same but every second pair
  if(abs(initSt)==10){
    mps.setProductState(p_zero);
    complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    mwArray oper=mwArray(Indices(d,d),datax);
    for(int k=0;k<L;k++){
      if(initSt>0){
	if(k%4>1) mps.applyLocalOperator(k,oper);
      }
      else{
	if(k%4<=1) mps.applyLocalOperator(k,oper);
      }
    }
  }
  //  mps.increaseBondDimension(D);
}

void __computeLocalEnergy(const MPS& mps,double J_,double g_,double h_,vector<double>& locE){
  locE.clear();
  int L=mps.getLength();
  Contractor& contractor=Contractor::theContractor();
  mwArray tmpL(ONE_c);tmpL.reshape(Indices(1,1));
  MPS state(mps);
  state.gaugeCond('L',1); // so that contraction from right is identity
  IsingHamiltonian hamLoc(2,d,J_,.5*g_,.5*h_);
  const MPO& mpo=hamLoc.getHMPO();
  MPO h12(4);
  h12.setOp(1,&mpo.getOp(0),false);
  h12.setOp(2,&mpo.getOp(1),false);
  for(int kl=0;kl<L-1;kl++){
    int Dl=mps.getA(kl).getDl();
    int Dr=mps.getA(kl+1).getDr();
    h12.setOp(0,new Operator(reshape(tmpL,Indices(Dl,1,Dl,1))),true);
    h12.setOp(3,new Operator(reshape(identityMatrix(Dr),Indices(Dr,1,Dr,1))),true);

    MPS aux(4,Dl,d);
    aux.replaceSite(0,reshape(identityMatrix(Dl),Indices(Dl,1,Dl)),false);
    aux.replaceSite(1,state.getA(kl).getA(),false);
    aux.replaceSite(2,state.getA(kl+1).getA(),false);
    aux.replaceSite(3,reshape(identityMatrix(Dr),Indices(Dr,Dr,1)),false);

    complex_t valE=contractor.contract(aux,h12,aux);
    cout<<"Computed value of E for link "<<kl<<" "<<valE<<endl;
    locE.push_back(real(valE));
    
    //contract one more site to tmpL
    if(kl<L-2){ // do by hand what Contractor.contractL does      
      int d=mps.getA(kl).getd();int Dl=mps.getA(kl).getDl();
      int Dr=mps.getA(kl).getDr();
      int betab=tmpL.getDimension(0);int betak=tmpL.getDimension(1);
      mwArray result=permute(mps.getA(kl).getA(),Indices(2,1,3));
      result.reshape(Indices(Dl,d*Dr));
      result=tmpL*result; // betab x d*Dr
      result.reshape(Indices(betab*d,Dr));
      mwArray aux=permute(mps.getA(kl).getA(),Indices(3,2,1)); // Dr x Dl x d (of bra)
      aux.conjugate();aux.reshape(Indices(Dr,Dl*d));
      tmpL=aux*result;
    }
  }
  
  
}



const string mpsFileName(const string& mpsdir,const string& baseName,const string& label,int cnt){
  stringstream s;
  s<<mpsdir<<"/"<<baseName<<"_"<<label;
  if(cnt>0) s<<"_"<<cnt;
  return s.str();
}

void tmpSave(MPS& mps,double norm_MPS,MPS& chebyT_Nm1,double norm_Nm1,MPS& chebyT_N,
	     double norm_N,const string& baseDir,const string& baseName,int cnt,vector<string>& tmpFiles){
  // first create new tmpFile names
  vector<string> existing(tmpFiles); // copy to remove at the end
  tmpFiles.clear();
  // now save ony by one
  string newTmp=mpsFileName(baseDir,baseName,"mps",cnt);
  double exNorm=mps.getNormFact();  
  mps.setNormFact(norm_MPS+log2(exNorm));
  mps.exportMPS(newTmp.data()); // saved containing log of norm!
  mps.setNormFact(exNorm); // restore
  tmpFiles.push_back(newTmp);
  
  newTmp=mpsFileName(baseDir,baseName,"Nm1",cnt);
  exNorm=chebyT_Nm1.getNormFact();  
  chebyT_Nm1.setNormFact(norm_Nm1+log2(exNorm));
  chebyT_Nm1.exportMPS(newTmp.data()); // saved containing log of norm!
  chebyT_Nm1.setNormFact(exNorm); // restore
  tmpFiles.push_back(newTmp);

  newTmp=mpsFileName(baseDir,baseName,"N",cnt);
  exNorm=chebyT_N.getNormFact();  
  chebyT_N.setNormFact(norm_N+log2(exNorm));
  chebyT_N.exportMPS(newTmp.data()); // saved containing log of norm!
  chebyT_N.setNormFact(exNorm); // restore
  tmpFiles.push_back(newTmp);

  cout<<"Saved tmp files for step "<<cnt<<" to ";
  for(int p=0;p<tmpFiles.size();p++)
    cout<<tmpFiles[p]<<", ";
  cout<<endl;
  // now remove the old one
  if(existing.size()>0){
    for(int k=0;k<existing.size();k++){
      string tmpF=existing[k];
      if(tmpF.length()>0&&file_exists(tmpF)){
	remove(tmpF.data());
      }
    }
  }
}


int tmpRestore(MPS& mps,double& norm_MPS,MPS& chebyT_Nm1,double& norm_Nm1,MPS& chebyT_N,
	     double& norm_N,const string& mpsdir,const string& baseName,int M,vector<string>& tmpFiles){
  tmpFiles.clear();
  // first find if there are tmp files
  int cnt=M;
  bool found=0;
  while(cnt>0&&!found){
    string mpsfile=mpsFileName(mpsdir,baseName,"mps",cnt);
    if(file_exists(mpsfile)){
      cout<<"findTmpFile found mpsfile="<<mpsfile<<endl;
      // Check also for the others
      string Nm1file=mpsFileName(mpsdir,baseName,"Nm1",cnt);
      string Nfile=mpsFileName(mpsdir,baseName,"N",cnt);
      if(file_exists(Nm1file)&&file_exists(Nfile)){
	// Read the files and return
	mps.importMPS(mpsfile.data());
	chebyT_Nm1.importMPS(Nm1file.data());
	chebyT_N.importMPS(Nfile.data());
	norm_MPS=mps.getNormFact();mps.setNormFact(1.);
	norm_Nm1=chebyT_Nm1.getNormFact();chebyT_Nm1.setNormFact(1.);
	norm_N=chebyT_N.getNormFact();chebyT_N.setNormFact(1.);
	tmpFiles.push_back(mpsfile);
	tmpFiles.push_back(Nm1file);
	tmpFiles.push_back(Nfile);
	found=true;
	return cnt;
      }
      else{
	cout<<"ERROR: Found a file for MPS but not the corresponding ones for Nm1 and N=> no tmp file used"<<endl;
	return 0;
      }
    }
    else{
      cnt--;
    }
  }
  if(!found) return 0;
  cout<<"ERROR: Should not be here"<<endl;
  exit(1);
}




double computeRDMentropy(const mwArray& rdm,int Lc){
  double entropy=0.;
  mwArray U;vector<complex_t> Dv;
  wrapper::eig(rdm,Dv,U); // only eigenvalues
  for(int ik=0;ik<Dv.size();ik++){
    double tmp=real(Dv[ik]);
    if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
  }
  return entropy;
}

double computeRDMdistanceToId(const mwArray& rdm,int Lc,double& entropy){
  double singvals=0.;entropy=0.;
  double idFc=1./pow(2,Lc);
  mwArray U;vector<complex_t> Dv;
  wrapper::eig(rdm,Dv,U); // only eigenvalues
  // now sum the absolute values
  for(int k=0;k<Dv.size();k++){
    double aux=real(Dv[k]);
    singvals+=abs(aux-idFc);
    if(abs(aux)>1E-12) entropy+=-aux*log2(aux);
  }
  if(Dv.size()<pow(2,Lc)){
    cout<<"WARNING! Seems the dimension of the operator for Lc="<<Lc
	<<"is not full: dim(oper)="<<Dv.size()
	<<", 2^Lc="<<pow(2,Lc)<<endl;
    singvals+=idFc*(pow(2,Lc)-Dv.size());
  }
  return singvals;
}


// Copied from Contractor, as no access
void contractLmiddle(mwArray& result,const mwArray& tmpL,
		     const Site& ket,const Site& bra){
  int d=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  if(tmpL.isEmpty()){ // first contraction, only MPS
    //cout<<"First contraction, tmpL empty; bra:"<<bra.getA().getDimensions()<<", ket:"<<ket.getA().getDimensions()<<endl;
    result=bra.getA();
    result.reshape(Indices(d,Dl*Dr));
    result.conjugate();
    mwArray aux=ket.getA();
    aux.reshape(Indices(dp,Dlp*Drp));
    aux.permute(Indices(2,1));
    result.multiplyLeft(aux);
    result.reshape(Indices(Dlp,Drp,Dl,Dr));
    result.permute(Indices(3,1,4,2));
    result.reshape(Indices(Dl,Dlp,Dr,Drp));
  }
  else{
    int Dtmp=tmpL.getDimension(0);int Dtmpp=tmpL.getDimension(1);
    int chi2=Dtmp*Dtmpp;
    int betab=tmpL.getDimension(2);int betak=tmpL.getDimension(3);
    result=permute(ket.getA(),Indices(2,1,3));
    result.reshape(Indices(Dlp,dp*Drp));
    mwArray aux0(tmpL); aux0.reshape(Indices(chi2*betab,betak));
    result.multiplyLeft(aux0);
    //cout<<"After multiplying in ket, "<<result.getDimensions()<<endl;
    result.reshape(Indices(chi2,betab*dp,Drp));
    result.permute(Indices(2,3,1));
    result.reshape(Indices(betab*dp,Drp*chi2));
    mwArray aux=bra.getA();
    aux.permute(Indices(3,2,1),true);
    aux.reshape(Indices(Dr,Dl*d));
    result.multiplyLeft(aux);
    //cout<<"After multiplying in bra, "<<result.getDimensions()<<endl;
    result.reshape(Indices(Dr,Drp,chi2));
    result.permute(Indices(3,1,2));
    result.reshape(Indices(Dtmp,Dtmpp,Dr,Drp));
  }
}

void contractRmiddle(mwArray& result,mwArray& tmpR,
		     const Site& ket,const Site& bra){
  int d=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int alphab=tmpR.getDimension(0);int alphak=tmpR.getDimension(1);
  int betab=tmpR.getDimension(2);int betak=tmpR.getDimension(3);
  tmpR.permute(Indices(3,4,1,2));
  contractLmiddle(result,tmpR,Site(0,permute(ket.getA(),Indices(1,3,2))),Site(0,permute(bra.getA(),Indices(1,3,2))));
  result.permute(Indices(3,4,1,2));
}




int computeBlockEntropies(MPS& mps_,vector<double>& trDistRDM,
			   vector<double>& entropyRDM,vector<int>& sizeRDM){
  cout<<"computeBlockEntropies "<<endl;
  int L=mps_.getLength();
  int Dcut=mps_.getBond();
  if(Dcut>200) Dcut=100;
  MPS mps(mps_,Dcut);
  
  mps.gaugeCond('R',1); // should be already normalized
  for(int k=L-1;k>=L/2;k--){
    mps.gaugeCond(k,'L',true);
  }
  // Now contraction from left or from right until mid chain gives ID
  // I compute each block starting from the middle one (2 sites), which is cheaper
  
  int Lc=2;
  double Sred=0.;
  double thermalEV=0.;
  double distLc=0;

  mwArray midContr;
  for(int posL=L/2-1;posL>=1;posL--){
    int posR=posL+Lc-1;
    cout<<"posL="<<posL<<", posR="<<posR<<", Lc="<<Lc<<endl;
    mwArray tmp;
    contractLmiddle(tmp,midContr,mps.getA(posR),mps.getA(posR));
    cout<<"After contracting site "<<posR<<" tmp:"<<tmp.getDimensions()<<endl;
    contractRmiddle(midContr,tmp,mps.getA(posL),mps.getA(posL));
    cout<<"After contracting site "<<posL<<" tmp:"<<midContr.getDimensions()<<endl;

    tmp=midContr;
    // Reshape as due
    int Dl=mps.getA(posL).getDl();
    int Dr=mps.getA(posR).getDr();
    tmp.reshape(Indices(Dl,Dl,Dr,Dr));
    tmp.permute(Indices(1,3,2,4));
    tmp.reshape(Indices(Dl*Dr,Dl*Dr));
    // And perform the non orthogonal diagonalization
    vector<complex_t> Dval;mwArray U; //placeholder:ignored
    wrapper::eig(tmp,Dval,U); // only eigenvalues
    // And now compute distance and entropy, too
    Sred=0.;distLc=0.;
    int nrVals=Dval.size();
    for(int k=0;k<nrVals;k++){
      double aux=real(Dval[k]);
      if(abs(aux)>1E-12) Sred+=-aux*log2(aux);
      distLc+=abs(aux-thermalEV);
    }
    distLc+=thermalEV*(pow(2,L)-nrVals);
    sizeRDM.push_back(Lc);  
    entropyRDM.push_back(Sred);
    trDistRDM.push_back(distLc);
    cout<<"For middle "<<Lc<<" sites ("<<posL<<" to "<<posR<<"), entropy="<<Sred<<", dis "<<distLc<<endl;
    Lc+=2;
  }
  // Last block
  Lc=L;
  Sred=0.;
  thermalEV=1./pow(2,Lc);
  distLc=2*(1.-thermalEV);
  sizeRDM.push_back(Lc);  
  entropyRDM.push_back(Sred);
  trDistRDM.push_back(distLc);
  
  return Dcut;
}


//constructOperatorProduct(mwArray& result,const mwArray& opA,const mwArray& opB);


void computeTwoBodyTerm(mwArray& result,double Jx,double Jy,double Jz,
			double hx1,double hx2,double hy1,double hy2,double hz1,double hz2){
  mwArray sig0=identityMatrix(d);//identity
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),datax);//sigmax
  complex_t datay[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigY(Indices(d,d),datay);//sigmay
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz
  
  mwArray XX,YY,ZZ,XI,IX,YI,IY,ZI,IZ;
  constructOperatorProduct(XX,sigX,sigX);
  constructOperatorProduct(YY,sigY,sigY);
  constructOperatorProduct(ZZ,sigZ,sigZ);
  constructOperatorProduct(XI,sigX,sig0);
  constructOperatorProduct(YI,sigY,sig0);
  constructOperatorProduct(ZI,sigZ,sig0);
  constructOperatorProduct(IX,sig0,sigX);
  constructOperatorProduct(IY,sig0,sigY);
  constructOperatorProduct(IZ,sig0,sigZ);

  result=Jx*XX+Jy*YY+Jz*ZZ+hx1*XI+hx2*IX+hy1*YI+hy2*IY+hz1*ZI+hz2*IZ;
  result.reshape(Indices(d,d,d,d));
  result.permute(Indices(1,3,2,4));
  result.reshape(Indices(d*d,d*d));
}

void splitTwoBodyTerms(const mwArray& h12,mwArray& Ol,mwArray& Or,mwArray& Omid12,mwArray& Omid21){
  mwArray S; // place for the singular values
  int nr=0;
  double tol=0.;
  //cout<<"Before decomposition, expH="<<expH<<endl;
  wrapper::svd(h12,tol,nr,Ol,S,Or);
  //  cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<", S="<<S<<endl;
  // redistribute the singular values
  mwArray sqrS=sqrt(S);
  // TODO: check no NaN???
  int xi=S.getDimension(0);
  Ol.multiplyRight(sqrS); // d*d*Xi
  Or.multiplyLeft(sqrS); // Xi*d*d
  // reshape for operators (i'x 1 x i x alpha ))
  Ol.reshape(Indices(d,1,d,-1));
  Or.reshape(Indices(-1,d,d,1));
  Or.permute(Indices(2,1,3,4));  

  Omid12=mwArray(Or); // d' xi d 1
  Omid12.reshape(Indices(d*xi,d));
  Omid12.multiplyRight(reshape(Ol,Indices(d,d*xi))); // d' *xi x d * xi
  Omid12.reshape(Indices(d,xi,d,xi));

  Omid21=mwArray(Ol); // d'*1*d*xi
  Omid21.permute(Indices(1,4,3,2)); //d' xi d 1
  Omid21.reshape(Indices(d*xi,d));
  Omid21.multiplyRight(reshape(Or,Indices(d,xi*d))); // d'*xi x xi(l)*d
  Omid21.reshape(Indices(d,xi,xi,d));
  Omid21.permute(Indices(1,3,4,2));
  
}

void computeEnergyEnergy(const MPS& mps,mwArray& h12,vector<int>& pos,vector<double>& h0hm,vector<double>& hm,int offset){
  cout<<"computeEnergyEnergy with h12="<<h12<<endl;
  h0hm.clear();hm.clear();pos.clear();
  // split the two body op as MPO
  mwArray Ol,Or,Omid12,Omid21;
  splitTwoBodyTerms(h12,Ol,Or,Omid12,Omid21);

   // // cout<<"After splitTwoBodyTerms,"<<endl;
   //  putForMatlab(cout,Ol,"Ol");
   //  putForMatlab(cout,Or,"Or");
   //  putForMatlab(cout,Omid12,"Omid12");
   //  //exit(1);

  Operator* Olop=new Operator(Ol);
  Operator* Orop=new Operator(Or);
  Operator* Omid=new Operator(Omid12);
  
  int L=mps.getLength();
  MPO corrHH(L);
  MPO corrHH1(L); // just for the first term
  MPO singleH(L);
  Operator* idOp=new Operator(reshape(identityMatrix(d),Indices(d,1,d,1)));
  for(int k=0;k<L;k++){
    corrHH.setOp(k,idOp,0);
    corrHH1.setOp(k,idOp,0);
    singleH.setOp(k,idOp,0);
  }
  // Now all MPOs are just identity
  corrHH1.setOp(L/2-1-offset,Olop,0);
  corrHH1.setOp(L/2-offset,Omid,0);
  corrHH1.setOp(L/2+1-offset,Orop,0);

  Contractor& contractor=Contractor::theContractor();

  // first the e.v. of h12 on the middle pair (L/2-1 and L/2) -offset
  corrHH.setOp(L/2-1-offset,Olop,0);
  corrHH.setOp(L/2-offset,Orop,0);
  double hm_=real(contractor.contract(mps,corrHH,mps));
  double h0hm_=real(contractor.contract2(corrHH,mps));

  // cout<<"In this state, <h12>="<<hm_<<" and <h12^2>="<<h0hm_<<endl;
  // cout<<"State is "<<mps<<endl;
  
  pos.push_back(L/2-1-offset);h0hm.push_back(h0hm_);hm.push_back(hm_);
  // the site L/2 is special, too
  singleH.setOp(L/2-offset,Olop,0);
  singleH.setOp(L/2+1-offset,Orop,0);
  hm_=real(contractor.contract(mps,singleH,mps));
  singleH.setOp(L/2-offset,idOp,0);
  singleH.setOp(L/2+1-offset,idOp,0);

  h0hm_=real(contractor.contract(mps,corrHH1,mps));
  pos.push_back(L/2-offset);
  hm.push_back(hm_);
  h0hm.push_back(h0hm_);

  for(int k=L/2+1-offset;k<L-1;k++){
    singleH.setOp(k,Olop,0);
    singleH.setOp(k+1,Orop,0);
    hm_=real(contractor.contract(mps,singleH,mps));    
    singleH.setOp(k,idOp,0);
    singleH.setOp(k+1,idOp,0);
    corrHH.setOp(k,Olop,0);
    corrHH.setOp(k+1,Orop,0);
    h0hm_=real(contractor.contract(mps,corrHH,mps));
    corrHH.setOp(k,idOp,0);
    corrHH.setOp(k+1,idOp,0);
    pos.push_back(k);hm.push_back(hm_);h0hm.push_back(h0hm_);
 
  }
  cout<<"At the end of computeEnergyEnergy, pos="<<pos<<endl;
  cout<<"At the end of computeEnergyEnergy, hm="<<hm<<endl;
  cout<<"At the end of computeEnergyEnergy, hmho="<<h0hm<<endl;
  delete idOp,Olop,Orop,Omid;
}

// Works, but not efficient (contracts repeatedly each block without reusing)
int computeBlockEntropiesSlow(MPS& mps_,vector<double>& trDistRDM,
			   vector<double>& entropyRDM,vector<int>& sizeRDM){
  cout<<"computeBlockEntropies "<<endl;
  int L=mps_.getLength();
  mps_.gaugeCond('R',1);
  int Dcut=mps_.getBond();
  MPS mps(mps_,min(Dcut,100));
  if(Dcut>100){
    Contractor& contractor=Contractor::theContractor();
    contractor.optimizeMPS(mps_,mps,100);
    Dcut=100;
  }
  
  // *********************
  // I compute from the boundary, inwards, removing 2 (one per site) every time, until only 2 are left
  // first, I set gaugeCond to right

  mps.gaugeCond('R',1); // should be already normalized
  for(int k=L-1;k>=L/2;k--){
    mps.gaugeCond(k,'L',true);
  }
  // Now contraction from left or from right until mid chain gives ID
  
  int Lc=L;
  double Sred=0.;
  double thermalEV=1./pow(2,Lc);
  double distLc=2*(1.-thermalEV);

  for(int lb=1;lb<=L/2-1;lb++){
    Lc-=2;

    // Now contract the part of the state, from lb to L-lb-1 with the bra, to obtain the
    // D^2 x D^2 overlap matrix for the L-2*lb block
    mwArray overl;
    for(int k=lb;k<=L-lb-1;k++){
      mwArray tmp(overl);
      contractLmiddle(overl,tmp,mps.getA(k),mps.getA(k));
      //      cout<<"After contracting in site "<<k<<", tmp is "<<overl.getDimensions()<<endl;
    } 
    // Reshape as due
    int Dl=mps.getA(lb).getDl();
    int Dr=mps.getA(L-lb-1).getDr();
    overl.reshape(Indices(Dl,Dl,Dr,Dr));
    overl.permute(Indices(1,3,2,4));
    overl.reshape(Indices(Dl*Dr,Dl*Dr));
    // And perform the non orthogonal diagonalization
    vector<complex_t> Dval;mwArray U; //placeholder:ignored
    wrapper::eig(overl,Dval,U); // only eigenvalues
    // And now compute distance and entropy, too
    Sred=0.;distLc=0.;
    for(int k=0;k<Dval.size();k++){
      double tmp=real(Dval[k]);
      if(abs(tmp)>1E-12) Sred+=-tmp*log2(tmp);
      distLc+=abs(real(Dval[k])-thermalEV);
    }
    distLc+=thermalEV*(pow(2,L)-Dval.size());
    cout<<"For middle "<<Lc<<" sites ("<<lb<<" to "<<L-lb-1<<"), entropy="<<Sred<<", dis "<<distLc<<endl;
    sizeRDM.push_back(Lc);  
    entropyRDM.push_back(Sred);
    trDistRDM.push_back(distLc);
  }
  return Dcut;
}

// for at most Lc=12 (max I can do exactly, I think), compute the
// exact RDM of the central Lc sites, its entropy and distance to the maximally mixed state.
// It could be possible to use a numerical approximation to the entropy for something larger

// Assuming proper "central" gauge
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


int computeBlockEntropiesExact(MPS& mps_,vector<double>& trDistRDM,
			   vector<double>& entropyRDM,vector<int>& sizeRDM){
  int Lcmax=10;
  cout<<"computeBlockEntropiesExact "<<endl;
  int L=mps_.getLength();
  int minPosR=L/2; // assuming L even and smallest Lc is 2
  mps_.gaugeCond('R',1);
  for(int k=L-1;k>=minPosR+1;k--){
    mps_.gaugeCond(k,'L',true);
  } // now both sides are identity
  int Lc=Lcmax;
  int posL=(L-Lc)/2;int posR=posL+Lc-1; // outside tracing them
  while(Lc>=2){
    mwArray oper;
    expandRDM(mps_,posL,posR,oper);
    cout<<"Created RDM for Lc="<<Lc<<endl;
    double singvals=0.;
    double idFc=1./pow(2,Lc);
    mwArray U;vector<complex_t> Dv;
    wrapper::eig(oper,Dv,U); // only eigenvalues
    double entropy=0.;
    for(int ik=0;ik<Dv.size();ik++){
      double tmp=real(Dv[ik]);
      if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
    }
    cout<<"Diagonalized the rho(Lc)"<<endl;
    // now sum the absolute values
    for(int k=0;k<Dv.size();k++){
      double aux=real(Dv[k]);
      singvals+=abs(aux-idFc);
    }
    if(Dv.size()<pow(2,Lc)){
      cout<<"WARNING! Seems the dimension of the RDM operator for D="<<mps_.getA(L/2).getDl()
	  <<", Lc="<<Lc<<"is not full: dim(oper)="<<Dv.size()
	  <<", 2^Lc="<<pow(2,Lc)<<endl;
      singvals+=idFc*(pow(2,Lc)-Dv.size());
    }
    sizeRDM.push_back(Lc);
    trDistRDM.push_back(singvals);
    entropyRDM.push_back(entropy);
    posL++;posR--;Lc=(posR-posL+1);
  }
  
}

