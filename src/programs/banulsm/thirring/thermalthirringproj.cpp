
#include <math.h>
#include <iomanip>
#include <sstream>

#include "misc.h"
#include "Properties.h"
//#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "SpinMPO.h"
#include "quicksort.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "ThirringHamiltonian.h"

const string mpsTmpFileNameRoot(const string baseDir,int D,double g2,double gtilde2,double ma,double Starget,int L,double delta,double tol);
// Add a particular counter
const string mpsTmpFileName(const string rootName,int cnt);
// And a final name
const string mpsSaveFileName(const string rootName,double beta);

void getStepNumbers(vector<int>& stepNrs,const vector<double>& targetBs,double delta);

void evolveExtra(const MPS& mpsPos,ThirringHamiltonian& hamH,const MPO& projS0,const MPO& Hproj,double extraB,
		 const string& nameToSavePos,int D);

using namespace shrt;
using namespace std;

/** thermalthirringproj approximates the thermal state of
    Thirring Hamiltonian \ref <ThirringHamiltonian>,
    for a certain inverse temperature beta,
    via imaginary time evolution.

    In this case, because we are interested in a sector of constant
    Sz, we project only when measuring observables, instead of
    starting with the projector.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <ma> (double) adimensional mass parameter \f$ma\f$ of the Hamiltonian
    \param <g2> (double) adimensional coupling parameter \f$g^2\f$ of the Hamiltonian
    \param <gt2> (double) adimensional coupling parameter \f$\tilde{g}^2\f$
    of the Hamiltonian, to allow for a different coupling for even links
    \param <lambda> (double) penalty term \f$\lambda\f$ (\f$>0\f$)
    \param <mu> (double) chemical potential \f$\mu\f$ 
    \param <Starget> (int) targetted total \f$S_z\f$
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the energy
    \param <beta> (double) inverse temperature
    \param <delta> (double) step width for imaginary time evolution
    \param <rate> (int) frequency with with results are stored
    \param <outfname> (char*) name of the output file for the results
    \param <app> (int) whether the output file is to be kept (app==1)
*/

int d=2;
complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
double TOLPROJ=1E-4;

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

  int L = props.getIntProperty("L");  
  if(L%2!=0){
    cout<<"ERROR: Nr of sites L should be even"<<endl;
    exit(1);
  }
  double ma=props.getDoubleProperty("ma");  
  double gt2_=props.getDoubleProperty("gtilde2");
  double g2_=props.getDoubleProperty("g2"); // the non-alternating one
  double lambda_=props.getDoubleProperty("lambda");
  if(lambda_<0) lambda_=0.;
  double mu_=props.getDoubleProperty("mu");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  int Stot_=props.getIntProperty("Starget");
  int D=props.getIntProperty("D");
  double beta=props.getDoubleProperty("beta");
  double delta=props.getDoubleProperty("delta");
  int rate=props.getIntProperty("stepBlock");
  int freqBackup=props.getIntProperty("backupFreq");
  if(freqBackup<=0) freqBackup=10;
  freqBackup=lcm(freqBackup,rate);
  cout<<"The frequency of backup decided to be "<<freqBackup<<endl;
  bool bothSigns=props.getIntProperty("onlyPositive")!=1; // default is both
  
  vector<double> targetBs=props.getVectorDoubleProperty("targets");
  if(targetBs.size()>0){
    // order, just in case
    quicksort(targetBs,0,targetBs.size()-1);
    if(targetBs[0]<0) bothSigns=1; // independently of option
  }
  
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir");
  if(freqBackup>0&&mpsdir.length()==0){
    cout<<"Error: no option mpsdir specified for tmp backups!"<<endl;
    exit(1);
  }
  
  bool app=(props.getIntProperty("append")!=0);
  bool findGS=(props.getIntProperty("findGS")==1);

  cout<<"Initialized arguments: L="<<L
      <<", ma="<<ma
      <<", g2="<<g2_
      <<", gt2="<<gt2_
      <<", lambda="<<lambda_
      <<", mu="<<mu_
      <<", Starget="<<Stot_
      <<", outfile="<<outfname
      <<", beta="<<beta
      <<", delta="<<delta
      <<", target betas: "<<targetBs
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    *out<<"% L="<<L<<", ma="<<ma<<", g2="<<g2_<<", gt2="<<gt2_<<", Starget="<<Stot_<<", D="<<D<<endl;
  }
  else
    out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }

  *out<<setprecision(10);
  int M=beta/(2*delta); // each step applies delta twice (each side)
  vector<int> stepNrs;
  getStepNumbers(stepNrs,targetBs,delta);
  cout<<"M="<<M<<", I will need to stop at stepNrs: "<<stepNrs<<"; also bothSigns="<<bothSigns<<endl;
  
  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  cout<<"Initialized Contractor"<<endl;
  contractor.setConvTol(tol);

  int d=2;
  int cnt=0;
  MPS thS(L,D,d*d);
  thS.setProductState(p_maxent); // the initial state, maximally entangled
  MPS thSneg(thS); // same, to evolve in negative beta

  MPO _projS0(L); // projector onto total spin 0
  SpinMPO::getProjectorTotalSz(L,d,_projS0); // TODO: generalize for other total spins
  // Notice:this is a projector with exponentially large Frobenius norm choose(N,N/2)
  cout<<"After creating the proj (as MPO) its norm is "<<contractor.contract(thS,thS)<<endl;

  // First, create the root name for backup files
  string backupFileNameBasis=mpsTmpFileNameRoot(mpsdir,D,g2_,gt2_,ma,Stot_,L,delta,tol);
  // Now find, if any, the latest backup file stored
  string backupFileP,backupFileM;
  bool found=false;cnt=M;
  while((!found)&&cnt>0){
    backupFileP=mpsTmpFileName(backupFileNameBasis,cnt);
    backupFileM=mpsTmpFileName(backupFileNameBasis,-cnt);
    if(file_exists(backupFileP)&&(file_exists(backupFileM)||!bothSigns)){
      found=true;
    }
    else{ backupFileP="";backupFileM="";cnt--;}
  }
  if(backupFileP.length()!=0&&(backupFileM.length()!=0||!bothSigns)){
    cout<<"Using backup files "<<backupFileP<<" and "<<backupFileM<<", step "<<cnt<<endl;
    thS.importMPS(backupFileP.data());
    if(bothSigns) thSneg.importMPS(backupFileM.data());
  }

  // cout<<"Initialized identity state"<<endl;
  //   // TODO: Substitute by projector onto right value of Stot!!

  //MPSfromMPO(_projS0,thS);
  // cout<<"After creating the proj (as MPO) its norm is "<<contractor.contract(thS,thS)<<endl;
  // thS.gaugeCond('R',true);    
  // cout<<"After gauge L its norm is "<<contractor.contract(thS,thS)<<endl;
  // thS.gaugeCond('L',true);    
  // cout<<"After gauge R its norm is "<<contractor.contract(thS,thS)<<endl;

  // which ofthe desired betas is next
  int ptr=0;
  if(cnt>0&&targetBs.size()>0){
    bool found=0;
    while(!found){
      if(stepNrs[ptr]>=cnt) found=true;
      else ptr++;
      if(ptr>=stepNrs.size()) found=true;
    }
  }
  cout<<"About to start evolution with ptr of stepNrs in "<<ptr<<"( out of "<<stepNrs.size()<<"), cnt="<<cnt<<endl;
  
  MPO Szmpo(L);
  SpinMPO::getSzMPO(L,d,Szmpo);

 
  ThirringHamiltonian hamH(L,ma,g2_,lambda_,Stot_,d,gt2_,mu_);
  //cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamH.getHMPO();

  MPS gs(L,D,d);
  double E0=0;
  if(findGS){
    contractor.findGroundState(hamil,D,&E0,gs);
    cout<<"Ground state energy="<<E0<<endl;
    if(!app) *out<<"% Ground state energy="<<E0<<endl;
  }
  else{
    if(!app) *out<<"% Ground state calculation skipped"<<endl;
  }
  if(!app){
    //*out<<"%\n%delta\t cnt\t beta\t T\t <H>\t <ProjS0>\t <Z(L/2-1)Z(L/2)>\t <Z(L/2-1)>\t <Z(L/2)>\t <Sztot>";
    *out<<"%\n%delta\t cnt\t beta\t T\t <H>\t <ProjS0>";
    if(findGS) *out<<"\t <ProjGS>";
    *out<<endl;
  }
  
  // Now try evolution with the exponential, instead
  MPO expHe(L);
  MPO expHo(L);
  MPO expHe2(L);
  MPO expHeneg(L);
  MPO expHoneg(L);
  MPO expHe2neg(L);
  MPO Hdoub(L);
  MPO projS0(L);
  //  MPO projS0_2(L); // projector on both sides: unnecessary!
  MPO Szdoub(L);
  MPO Hproj(L); // H on one side, projector on the other
  
  // I construct the simple ones first, and then from them the DoubleOperators
  {
    // Regular ones
    MPO _expHe2(L),_expHe(L),_expHo(L);
    hamH.getExponentialMPOeven(_expHe2,-delta*.5*ONE_c);
    hamH.getExponentialMPOeven(_expHe,-delta*ONE_c);
    hamH.getExponentialMPOodd(_expHo,-delta*ONE_c);

    MPO _expHe2neg(L),_expHeneg(L),_expHoneg(L);
    hamH.getExponentialMPOeven(_expHe2neg,delta*.5*ONE_c);
    hamH.getExponentialMPOeven(_expHeneg,delta*ONE_c);
    hamH.getExponentialMPOodd(_expHoneg,delta*ONE_c);

    // _expHe.exportForMatlab("expHe.m");
    // _expHo.exportForMatlab("expHo.m");
    // hamil.exportForMatlab("Htirr.m");

    
    mwArray identPhys=identityMatrix(d);
    identPhys.reshape(Indices(d,1,d,1)); 
    for(int k=0;k<L;k++){
      mwArray aux=_expHe2.getOp(k).getFullData();
      expHe2.setOp(k,new DoubleOperator(_expHe2.getOp(k).getFullData(),
					permute(aux,Indices(3,2,1,4))),true);
      aux=_expHe.getOp(k).getFullData();
      expHe.setOp(k,new DoubleOperator(_expHe.getOp(k).getFullData(),
				       permute(aux,Indices(3,2,1,4))),true);
      aux=_expHo.getOp(k).getFullData();
      expHo.setOp(k,new DoubleOperator(_expHo.getOp(k).getFullData(),
				       permute(aux,Indices(3,2,1,4))),true);

      aux=_expHe2neg.getOp(k).getFullData();
      expHe2neg.setOp(k,new DoubleOperator(_expHe2neg.getOp(k).getFullData(),
					permute(aux,Indices(3,2,1,4))),true);
      aux=_expHeneg.getOp(k).getFullData();
      expHeneg.setOp(k,new DoubleOperator(_expHeneg.getOp(k).getFullData(),
				       permute(aux,Indices(3,2,1,4))),true);
      aux=_expHoneg.getOp(k).getFullData();
      expHoneg.setOp(k,new DoubleOperator(_expHoneg.getOp(k).getFullData(),
				       permute(aux,Indices(3,2,1,4))),true);
      
      Hdoub.setOp(k,new DoubleOperator(hamil.getOp(k).getFullData(),identPhys),true);
      projS0.setOp(k,new DoubleOperator(_projS0.getOp(k).getFullData(),identPhys),true);
      aux=_projS0.getOp(k).getFullData();
      //      projS0_2.setOp(k,new DoubleOperator(_projS0.getOp(k).getFullData(),
      //				  permute(aux,Indices(3,2,1,4))),true);
      Hproj.setOp(k,new DoubleOperator(hamil.getOp(k).getFullData(),
				       permute(aux,Indices(3,2,1,4))),true);
      Szdoub.setOp(k,new DoubleOperator(Szmpo.getOp(k).getFullData(),identPhys),true);
    }
    // The simple ones will be destroyed now. I already have them copied
    // what if I instead do a joined operator
    const MPO* toJoin[]={&Hdoub,&projS0};
    MPO::join(2,toJoin,Hproj);

    
  }
  cout<<"Created all the exponential operators"<<endl;
  // Also the Z operator
  mwArray sigZ(Indices(d*d,1),dataZ);
  mwArray idOp=identityMatrix(d);idOp.reshape(Indices(1,d*d));
  sigZ.multiplyRight(idOp);sigZ.reshape(Indices(d,d,d,d));
  sigZ.permute(Indices(1,3,2,4));sigZ.reshape(Indices(d*d,d*d));

  // // To compute the exp. values, I need to apply the projector on both sides, so I construct a joined MPO
  // const MPO* ptrsH[]={&projS0,&Hdoub}; // enough on one side?
  // MPO hamProj(L);
  // MPO::join(2,ptrsH,hamProj);// MPO for H^2
  
  while(cnt<=M){
    *out<<delta<<"\t";
    *out<<cnt<<"\t";
    *out<<delta*cnt*2.<<"\t";
    *out<<1./(2.*delta*cnt)<<"\t";
    complex_t energy=contractor.contract(thS,Hproj,thS);
    *out<<real(energy)<<"\t";
    complex_t weightS0=contractor.contract(thS,projS0,thS);
    *out<<real(weightS0)<<"\t";
    // // Trick to compute Sz Sz and each of them
    // MPS Z1(thS);
    // Z1.applyLocalOperator(L/2-1,sigZ);
    // MPS Z2(thS);
    // Z2.applyLocalOperator(L/2,sigZ);
    // complex_t zz=contractor.contract(Z1,Z2);
    // complex_t z1=contractor.contract(Z1,thS);
    // complex_t z2=contractor.contract(thS,Z2);
    // *out<<real(zz)<<"\t";
    // *out<<real(z1)<<"\t";
    // *out<<real(z2)<<"\t";
    // complex_t valSz=contractor.contract(thS,Szdoub,thS);
    // *out<<real(valSz)<<"\t";
    if(findGS){
    // Also, what was the overlap with the GS?
      MPO gibbs(L);
      MPOfromMPS(thS,gibbs);
      complex_t projGS=contractor.contract2(gibbs,gs);
      *out<<abs(projGS)<<"\t";
    }
    *out<<endl;
    
    // And also the negative one
    if(bothSigns&&cnt>0){
      *out<<delta<<"\t";
      *out<<cnt<<"\t";
      *out<<-delta*cnt*2.<<"\t";
      *out<<-1./(2.*delta*cnt)<<"\t";
      complex_t energyN=contractor.contract(thSneg,Hproj,thSneg);
      *out<<real(energyN)<<"\t";
      complex_t weightS0neg=contractor.contract(thSneg,projS0,thSneg);
      *out<<real(weightS0neg)<<"\t";
      if(findGS){
	// Also, what was the overlap with the GS?
	MPO gibbsN(L);
	MPOfromMPS(thSneg,gibbsN);
	complex_t projGSneg=contractor.contract2(gibbsN,gs);
	*out<<abs(projGSneg)<<"\t";
      }
      *out<<endl;
    }

    // in the special case I am at a target stepNr, now I need to
    // evolve for the remainder imaginary time, save the state and
    // compute whatever observables. Then increment the ptr
    if(ptr<stepNrs.size()){
      while(ptr<stepNrs.size()&&cnt==stepNrs[ptr]){      
	string nameMPSpos=mpsSaveFileName(backupFileNameBasis,targetBs[ptr]);
	//string nameMPSneg=mpsSaveFileName(backupFileNameBasis,-targetBs[ptr]);

	double extraB=targetBs[ptr]-delta*cnt*2;
	cout<<"At cnt "<<cnt<<", evolving extra to reach beta="<<targetBs[ptr]<<endl;
	evolveExtra(thS,hamH,projS0,Hproj,extraB,nameMPSpos,D); //nameMPSneg,D);
	ptr++;
      }
    }    
    
    // how many asteps left in total, until next desired step
    int nrLeft=ptr<stepNrs.size()?stepNrs[ptr]-cnt:M-cnt;
    int toapply=min(rate,nrLeft);

    cout<<"At cnt="<<cnt<<"(of "<<M<<"), beta="<<delta*cnt*2<<", next stepNr:"<<(ptr<stepNrs.size()?stepNrs[ptr]:0)<<", nrLeft="<<nrLeft<<" nr to apply "<<toapply<<endl;
    MPS aux(thS); // temporary copy
    MPS auxN(thSneg); // temporary copy
    // cout<<"Created a tmp copy of the state with norm "<<contractor.contract(aux,aux)<<endl;
    // cout<<"Orig:"<<contractor.contract(thS,thS)<<endl;
    // cout<<"Overlap:"<<contractor.contract(thS,aux)<<endl;
    contractor.optimize(expHe2,aux,thS,D);thS.setNormFact(1.);
    if(bothSigns) contractor.optimize(expHe2neg,auxN,thSneg,D);thSneg.setNormFact(1.);
    //    cout<<"Applied the first exponential (e2)"<<endl;
    // Now apply r-1 times the pair Ho He
    int cntLoop=0;
    while(cntLoop<toapply-1&&cnt<M-1){
      //cout<<"In the loop cntLoop="<<cntLoop<<" (cnt="<<cnt<<")"<<endl;
      contractor.optimize(expHo,thS,aux,D);aux.setNormFact(1.);
      contractor.optimize(expHe,aux,thS,D);thS.setNormFact(1.);
      //cout<<"Time evolution step nr "<<cnt+1<<", time="<<time<<", value "<<M1ev<<", trace "<<trR<<endl;
      if(bothSigns) contractor.optimize(expHoneg,thSneg,auxN,D);auxN.setNormFact(1.);
      if(bothSigns) contractor.optimize(expHeneg,auxN,thSneg,D);thSneg.setNormFact(1.);
      cnt++;cntLoop++;
    }
    //cout<<"Out of the loop, norms thS:"<<thS.getNormFact()<<", thSneg:"<<thSneg.getNormFact()<<endl;
    contractor.optimize(expHo,thS,aux,D);aux.setNormFact(1.);
    contractor.optimize(expHe2,aux,thS,D);thS.setNormFact(1.);
    thS.gaugeCond('R',true);
    if(bothSigns){
      contractor.optimize(expHoneg,thSneg,auxN,D);auxN.setNormFact(1.);
      contractor.optimize(expHe2neg,auxN,thSneg,D);thSneg.setNormFact(1.);
      thSneg.gaugeCond('R',true);    
    }
    cnt++;
    
    // // from time to time, I should project back to the right sector
    // BUT ONLY IF THe DEVIATION IS LARGE (o.w. extra cost and probably truncation error)
    // if((cnt%10==0)){
    //   cout<<"Applying projector at cnt="<<cnt<<" (rate="<<rate<<")"<<endl;
    //   MPS aux(thS); // temporary copy
    //   contractor.optimize(projS0,aux,thS,D);
    //   thS.gaugeCond('R',true);    
    //   MPS auxN(thSneg); // temporary copy
    //   contractor.optimize(projS0,auxN,thSneg,D);
    //   thSneg.gaugeCond('R',true);    
    // }
          // if too much errror in the projector:
    weightS0=contractor.contract(thS,projS0,thS);
    if(abs(weightS0-ONE_c)>=TOLPROJ){
      cout<<"Applying projector to thS at cnt="<<cnt<<" (weightS0="<<weightS0<<")"<<endl;
      MPS aux(thS); // temporary copy
      contractor.optimize(projS0,aux,thS,D);
      thS.gaugeCond('R',true);    
    }
    if(bothSigns){
      complex_t weightS0neg=contractor.contract(thSneg,projS0,thSneg);
      if(abs(weightS0neg-ONE_c)>=TOLPROJ){
	cout<<"Applying projector to thSneg at cnt="<<cnt<<endl;
	MPS auxN(thSneg); // temporary copy
	contractor.optimize(projS0,auxN,thSneg,D);
	thSneg.gaugeCond('R',true);
      }
    }
    // Now decide if I need to save the backup file
    if(cnt>0&&(cnt%freqBackup==0||cnt==M)){
      cout<<"Saving backup at cnt="<<cnt<<endl;
      string oldBackupFileP(backupFileP);
      string oldBackupFileM(backupFileM);
      backupFileP=mpsTmpFileName(backupFileNameBasis,cnt);
      backupFileM=mpsTmpFileName(backupFileNameBasis,-cnt);
      thS.exportMPS(backupFileP.data());
      if(bothSigns) thSneg.exportMPS(backupFileM.data());
      // removing the old one
      if(file_exists(oldBackupFileP.data())){
	cout<<"Removing old backup file "<<oldBackupFileP<<endl;
	remove(oldBackupFileP.data());
      }
      if(bothSigns&&file_exists(oldBackupFileM.data())){
	cout<<"Removing old backup file "<<oldBackupFileM<<endl;
	remove(oldBackupFileM.data());
      }
    }
   }

   out->close();
}


const string mpsTmpFileNameRoot(const string baseDir,int D,double g2,double gtilde2,double ma,double Starget,int L,double delta,double tol){
  stringstream s;
  s<<baseDir<<"/tmpMPS_D"<<D<<"_g2"<<g2<<"_gt2"<<gtilde2<<"_ma"<<ma
   <<"_S"<<Starget<<"_L"<<L<<"_delta"<<delta
   <<"_eps"<<int(-log10(tol));
  s<<"_";
  return s.str();
}

const string mpsTmpFileName(const string rootName,int cnt){
  stringstream s;
  s<<rootName<<cnt<<".dat";
  return s.str();
}

const string mpsSaveFileName(const string rootName,double beta){
  stringstream s;
  s<<rootName<<"beta"<<beta<<".dat";
  return s.str();
}

void getStepNumbers(vector<int>& stepNrs,const vector<double>& targetBs,double delta){
  // betas are ordered
  // M=betaM/(2*delta)
  if(targetBs.size()>0){
    stepNrs=vector<int>(targetBs.size());
    for(int k=0;k<targetBs.size();k++){
      int nr=floor(targetBs[k]/(2*delta));
      stepNrs[k]=nr; //.push_back(nr);
    }
  }
}

void evolveExtra(const MPS& mpsPos,ThirringHamiltonian& hamH,const MPO& projS0,const MPO& Hproj,double extraB,
		 const string& nameToSavePos,int D){
		 //const string& nameToSaveNeg,int D){
  int L=mpsPos.getLength();
  double delta=extraB*.5; // half on each side
  MPO expHo(L);
  MPO expHe2(L);
  // MPO expHeneg(L);
  // MPO expHoneg(L);
  // MPO expHe2neg(L);
  // I construct the simple ones first, and then from them the DoubleOperators
  {
    // Regular ones
    MPO _expHe2(L),_expHo(L);
    hamH.getExponentialMPOeven(_expHe2,-delta*.5*.5*ONE_c);
    //    hamH.getExponentialMPOeven(_expHe,-delta*.5*ONE_c);
    hamH.getExponentialMPOodd(_expHo,-delta*.5*ONE_c);

    // MPO _expHe2neg(L),_expHeneg(L),_expHoneg(L);
    // hamH.getExponentialMPOeven(_expHe2neg,delta*.5*.5*ONE_c);
    // hamH.getExponentialMPOeven(_expHeneg,delta*.5*ONE_c);
    // hamH.getExponentialMPOodd(_expHoneg,delta*.5*ONE_c);
    
    mwArray identPhys=identityMatrix(d);
    identPhys.reshape(Indices(d,1,d,1)); 
    for(int k=0;k<L;k++){
      mwArray aux=_expHe2.getOp(k).getFullData();
      expHe2.setOp(k,new DoubleOperator(_expHe2.getOp(k).getFullData(),
					permute(aux,Indices(3,2,1,4))),true);
      // aux=_expHe.getOp(k).getFullData();
      // expHe.setOp(k,new DoubleOperator(_expHe.getOp(k).getFullData(),
      // 				       permute(aux,Indices(3,2,1,4))),true);
      aux=_expHo.getOp(k).getFullData();
      expHo.setOp(k,new DoubleOperator(_expHo.getOp(k).getFullData(),
				       permute(aux,Indices(3,2,1,4))),true);

      // aux=_expHe2neg.getOp(k).getFullData();
      // expHe2neg.setOp(k,new DoubleOperator(_expHe2neg.getOp(k).getFullData(),
      // 					permute(aux,Indices(3,2,1,4))),true);
      // aux=_expHeneg.getOp(k).getFullData();
      // expHeneg.setOp(k,new DoubleOperator(_expHeneg.getOp(k).getFullData(),
      // 				       permute(aux,Indices(3,2,1,4))),true);
      // aux=_expHoneg.getOp(k).getFullData();
      // expHoneg.setOp(k,new DoubleOperator(_expHoneg.getOp(k).getFullData(),
      // 				       permute(aux,Indices(3,2,1,4))),true);
      
    }
  }

  Contractor& contractor=Contractor::theContractor();
  MPS thS(mpsPos);//MPS thSneg(mpsNeg);
  MPS aux(thS);//MPS auxN(thSneg);
  contractor.optimize(expHe2,mpsPos,thS,D);
  //contractor.optimize(expHe2neg,mpsNeg,thSneg,D);
  contractor.optimize(expHo,thS,aux,D);
  //contractor.optimize(expHoneg,thSneg,auxN,D);
  contractor.optimize(expHe2,aux,thS,D);
  //contractor.optimize(expHe2neg,auxN,thSneg,D);

  // I should also project=> they end up in aux
  contractor.optimize(projS0,thS,aux,D);
  aux.gaugeCond('R',true);    
  //contractor.optimize(projS0,thSneg,auxN,D);
  //auxN.gaugeCond('R',true);    

  aux.exportMPS(nameToSavePos.data());
  //auxN.exportMPS(nameToSaveNeg.data());

  cout<<"Computed and saved state for additional beta="<<extraB<<" E=";
  complex_t energy=contractor.contract(aux,Hproj,aux);
  cout<<real(energy)<<"\t P(S0)=";
  complex_t weightS0=contractor.contract(aux,projS0,aux);
  cout<<real(weightS0)<<endl;

}
