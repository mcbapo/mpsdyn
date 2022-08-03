#include <math.h>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string.h>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "SchwingerHamiltonian.h"
#include "Properties.h"

#include <vector>
using namespace std;
using namespace shrt;

#define LOCALDIR "./"

//const string mpsTmpFileName(const string baseDir,int D,double x,double mg,int L,double delta,int cnt,int deltaChg);
const string mpsTmpFileNameRoot(const string baseDir,int D,double x,double mg,int L,double delta,double tol,int Lmax);
// Add a particular counter
const string mpsTmpFileName(const string rootName,int cnt);

// wrapper for the contractor function
void Optimize(const  MPO& mpo,const MPS& initMps,MPS& destMPS,int D);

// Useful functions to find gcd and lcm of two integer numbers
int gcd(int a,int b);
int lcm(int a,int b);

/** thermalWithExpL runs imaginary time evolution to approximate the
    thermal equilibrium state of the Schwinger Hamiltonian \ref
    <SchwingerHamiltonian>, and saves values of imaginary time,
    crresponding g beta, energy and condensate.

    Receives as argument the name of a Properties file, where
    arguments are read.  Additional parameters can be given to
    override the values in the Properties file.

    This version keeps track of the norm, to be able to cmopute the partition function.


 */

int main(int argc,const char* argv[]){

  int cntr=0;
  const char* infile=argv[++cntr];
  string directory="";
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  double alpha(0.), L0(0.);
  int L = props.getIntProperty("L");
  double mg=props.getDoubleProperty("mg");
  double x=props.getDoubleProperty("x");
  directory=props.getProperty("outputdir");
  const string outfname=directory+"/"+props.getProperty("outfile");
  int app=props.getIntProperty("append");
  if(app<0) app=1; // default will be to append
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-5;
  int D=props.getIntProperty("D");
  int M=props.getIntProperty("M"); // nr of steps
  int Lmax=props.getIntProperty("Lmax"); // cutoff for max L value on the links when computing expHL
  //  double beta=props.getDoubleProperty("beta");
  //double stepBeta=props.getDoubleProp("stepBeta"); // interval between recorded values of g beta
  double delta=props.getDoubleProperty("delta");
  string mpsdir=props.getProperty("mpsdir"); // default current
  if(mpsdir.empty()) 
    mpsdir=".";
  int freqBackup=props.getIntProperty("freqBkup");
  int blockSteps=props.getIntProperty("blockSteps"); // how many steps to block before measuring values

  // The frequency of backups will have to be a multiple of the steps blocked=> least common multiple of both numbers
  freqBackup=lcm(freqBackup,blockSteps);
  cout<<"The fequanc of backup decided to be "<<freqBackup<<endl;

  double mu=2*mg*sqrt(x);
  int d=2;
  ofstream *out;  // output file

  Contractor& contractor = Contractor::theContractor();
  contractor.setConvTol(tol);
  MPS imaGS(L,D,d*d);
  double E0=0.; double Norm=1.;
  complex_t CondSch;

  double logZ=0; // place for the log of the norm
  
  SchwingerHamiltonian hamSch(L, mu, x, L0, alpha);

  MPO Hdoub(L), Conddoub(L);   // Build double operators for observables
  {
    const MPO& hamil=hamSch.getHMPO();
    MPO CondMPO(L);
    hamSch.constructCondensateMPO(CondMPO);
    extendMPO(hamil,Hdoub,d);
    extendMPO(CondMPO,Conddoub,d);
  }


  // First, create the root name for backup files
  string backupFileNameBasis=mpsTmpFileNameRoot(mpsdir,D,x,mg,L,delta,tol,Lmax);
  // Now find, if any, the latest backup file stored
  string backupFile;
  bool found=false;int cnt=M;
  while((!found)&&cnt>0){
    backupFile=mpsTmpFileName(backupFileNameBasis,cnt);
    if(file_exists(backupFile)){
      found=true;
    }
    else{ backupFile="";cnt--;}
  }

  double time=2*cnt*delta;
  double gbeta=2*sqrt(x)*time;

  // Now, if found, I initialize the state with it, and the step with the corresponding count
  if(backupFile.length()!=0){
    cout<<"Using backup file "<<backupFile<<", step "<<cnt<<endl;
    imaGS.importMPS(backupFile.data());
    logZ=imaGS.getNormFact();
    imaGS.gaugeCond('R',1);
    imaGS.gaugeCond('L',1);
  }
  else{
    // otherwise, start from scratch
    cout<<"Starting from the beginning"<<endl;
    imaGS.setProductState(p_maxent);
    //Norm=real(contractor.contract(imaGS, imaGS)); // should be 2^L
    imaGS.gaugeCond('R',1);
    imaGS.gaugeCond('L',1);
    // initial values
    Norm=real(contractor.contract(imaGS, imaGS)); // just set to 1
    E0=real(contractor.contract(imaGS, Hdoub, imaGS))/Norm;
    CondSch=contractor.contract(imaGS, Conddoub, imaGS)/Norm;
    logZ=.5*L*log(2); // sqrt of the norm
    
    out=new ofstream(outfname.data());
    *out<<"% L : "<<L
	<<", mg : "<<mg
	<<", x : "<<x
	<<", D : "<<D
	<<", delta : "<<delta<<endl;
    *out<<"% backup files : "<<backupFileNameBasis<<endl;
    *out<<"% count\t time\t g beta \t energy \t condensate "<<endl;
    *out<<setprecision(15)<<cnt<<"\t"<<time<<"\t"<<gbeta<<"\t"
	<<E0<<"\t"<<real(CondSch)<<"\t"<<2*logZ<<endl;
    out->close();
    delete out;
  }

  // Before starting the evolution, prepare the required exponential operators
  // Cutoff for exponential of HL!
  int maxXi=min(L-1,2*Lmax+1);double tolHL=1E-16;

  // operators
  MPO evolxo_half(L), evolxe(L), evolxe_half(L);
  MPO evolHLexS(L),evolHLexA(L); // to separate the e(-beta L^2) MPOs on system and ancilla in two MPOs to be applied sequentially
  
  { // I do not need the single MPOs later
    MPO auxExp(L); // auxiliary placeholder for single layer MPOs
    hamSch.getExponentialMPOx(auxExp, -delta*.5*.5*ONE_c,true); // He half
    doubleMPO(auxExp,evolxe_half);
    hamSch.getExponentialMPOx(auxExp, -delta*.5*ONE_c,true); // He
    doubleMPO(auxExp,evolxe);
    hamSch.getExponentialMPOx(auxExp, -delta*.5*.5*ONE_c,false); // Ho half
    doubleMPO(auxExp,evolxo_half);
    hamSch.getExactExponentialMPOL(auxExp,-.5*delta*ONE_c,maxXi,tolHL,Lmax); // HL
    extendMPO(auxExp,evolHLexS,d);
    extendTransposeMPO(auxExp,evolHLexA,d);
  }


  //======================
  // TIME EVOLUTION PART
  //======================
  cout<<"=== STARTING TIME EVOLUTION FROM STEP "<<cnt<<" ==="<<endl;

  int stepsLeft=M-cnt;
  while(stepsLeft>0){
    int nrToApply=min(stepsLeft,blockSteps);
    MPS aux(imaGS);
    Optimize(evolxe_half, aux, imaGS, D);
    // Apply nrToApply-1 steps, because for the last one I need to end in evolxe_half
    for(int k=1;k<nrToApply;k++){
      Optimize(evolxo_half, imaGS, aux, D);
      Optimize(evolHLexS, aux, imaGS, D);
      Optimize(evolHLexA, imaGS, aux, D);
      Optimize(evolxo_half, aux, imaGS, D);aux=imaGS;
      Optimize(evolxe, aux, imaGS, D);
      cnt++;stepsLeft--;
    }
    // And the last one of the block
    Optimize(evolxo_half, imaGS, aux, D);
    Optimize(evolHLexS, aux, imaGS, D);
    Optimize(evolHLexA, imaGS, aux, D);
    Optimize(evolxo_half, aux, imaGS, D);aux=imaGS;
    Optimize(evolxe_half, aux, imaGS, D);
    cnt++;stepsLeft--;
    time=2*delta*cnt;
    gbeta=time*2*sqrt(x);

    // Normalize again, but keep the norm
    imaGS.gaugeCond('R',false);
    // store the factor
    logZ+=log(imaGS.getNormFact());
    imaGS.setNormFact(1.);
    // And compute and record observables
    // This is not really necessary: after normalizing, it has to be 1.
    Norm=real(contractor.contract(imaGS, imaGS));
    E0=real(contractor.contract(imaGS, Hdoub, imaGS))/Norm;
    CondSch=contractor.contract(imaGS, Conddoub, imaGS)/Norm;

    out=new ofstream(outfname.data(),ios::app);
    *out<<setprecision(15)<<cnt<<"\t"<<time<<"\t"<<gbeta<<"\t"
	<<E0<<"\t"<<real(CondSch)<<"\t"<<2*logZ<<endl;
    out->close();    
    delete out;

    // Now decide if I need to save the backup file
    if(cnt>0&&(cnt%freqBackup==0||cnt==M)){
      cout<<"Saving backup at cnt="<<cnt<<endl;
      string oldBackupFile(backupFile);
      backupFile=mpsTmpFileName(backupFileNameBasis,cnt);
      imaGS.setNormFact(logZ);
      imaGS.exportMPS(backupFile.data());
      imaGS.setNormFact(1.);
      // removing the old one
      if(file_exists(oldBackupFile.data())){
	cout<<"Removing old backup file "<<oldBackupFile<<endl;
	remove(oldBackupFile.data());
      }
    }
  }

}


/** The name of the file for the tmp storage of the MPs every few
    steps should be generated automatically, but such that in is
    unique for the simulation at hand. I choose the following format:
    tmpMPS_D[D]_x[x]_L[L]_delta[delta]_[cnt] (if this is the first delta)

    WARNING!!! This is not unambiguous, because I do not keep record
    in the name of which were the earlier deltas and when did the change
    happen, but I hope it is enough for my set of runs.
*/
// const string mpsTmpFileName(const char* baseDir,int D,double x,double mg,int L,double delta,int cnt){
//   stringstream s;
//   s<<baseDir<<"/tmpMPS_D"<<D<<"_x"<<x<<"_L"<<L<<"_mg"<<mg<<<<"_delta"<<delta<<"_";
//   if(cnt>0) s<<"_"<<cnt;
//   return s.str();
// }

const string mpsTmpFileNameRoot(const string baseDir,int D,double x,double mg,int L,double delta,double tol,int Lmax){
  stringstream s;
  s<<baseDir<<"/tmpMPS_D"<<D<<"_x"<<x<<"_L"<<L<<"_mg"<<mg<<"_delta"<<delta
   <<"_Lcut"<<Lmax<<"_eps"<<int(-log10(tol));
  s<<"_";
  return s.str();
}

const string mpsTmpFileName(const string rootName,int cnt){
  stringstream s;
  s<<rootName<<cnt<<".dat";
  return s.str();
}

void Optimize(const  MPO& mpo,const MPS& initMps,MPS& destMPS,int D){
  Contractor& contractor = Contractor::theContractor();
  contractor.optimize(mpo,initMps,destMPS,D);
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
// }

// int lcm(int a,int b){
//   return (a*b)/gcd(a,b);
// }
