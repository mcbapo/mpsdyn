
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"
//#include "HermitianTensorMultiplier.h"

#include "SoftCoulombHamiltonian.h"
//#include "SpinMPO.h"
#include <cmath>

using namespace std;
using namespace shrt;

#define MAXLEN 160

/** 
Generates a file name for the MPS corresponding to a certain excitation.
*/
const string mpsfilename(int L,int D,int M,int Ze,int level,const string mpsdir);
const string observablefilename(int L,int D,int M,int Ze,const string outputdir);

/** Given a state, compute all the expectation values (on each
    position) of a given single body operator. The results are
    written to the ostream (could be a file, or cout). 
    Returns 0 if values are ok, -1 if some site is outside the limits
*/
int computeSingleBodyEV(const MPS& exc,const mwArray& op1,ostream& os,double low=0.,double high=0.);
int computeCrossedSingleBodyEV(const MPS& gs,const MPS& exc,const mwArray& op1,vector<complex_t>& vals,double low=0.,double high=0.);

int d=2;

/** Reads states corresponding to SoftCoulombHamiltonian, as specified
    by a Properties file, and computes densities (maybe to be extended
    with other observables).

    Receives as argument a Properties file, such as
    config/softCoulombProp.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.

    To replace some of the values in the file by a command line
    option, the additional arguments have to be of the form
    -propName=propValue
    (except for the vector properties)
*/

int main(int argc,const char* argv[]){
  // Read input arguments
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

  // Store the command line, for this is what should be relaunched if not all levels are done
  stringstream commandstr;
  for(int c=0;c<argc;c++) commandstr<<argv[c]<<" ";
  // Now, str contains the command line, that I need to relaunch, if this was not the last run

  int L = props.getIntProperty("L");
  double M=props.getDoubleProperty("M");
  double Z=props.getDoubleProperty("Ze");
  double a=props.getDoubleProperty("a");
  double Delta=props.getDoubleProperty("Delta");
  directory=props.getProperty("outputdir");
  int app=props.getIntProperty("append");
  if(app<0) app=1; // default will be to append
  double tol=props.getDoubleProperty("tol");
  // int unif=1; // All same D 
  // int unif=props.getIntProperty("unif");
  double Npen=props.getDoubleProperty("penaltyNe");
  int Ntot=props.getIntProperty("Ne");
  double Spen=props.getDoubleProperty("penaltyS");
  if(Spen<0) Spen=0;
  int D=props.getIntProperty("D");
  string mpsdir=props.getProperty("mpsdir"); // default current
  if(mpsdir.empty()) 
    mpsdir=".";
  int nLev=props.getIntProperty("nLev");
  if(nLev<0) nLev=0; // where to stop
  int stopAtGS=(nLev==0);

  double offset=props.getDoubleProperty("offset");
  if(offset==-1) offset=0;

  cout<<"Initialized arguments: L="<<L
      <<", Delta="<<Delta
      <<", Ze="<<Z
      <<", a="<<a
      <<", M="<<M
      <<", app="<<app
      <<", tol="<<tol
      <<", penaltyNe="<<Npen
      <<", Ne="<<Ntot
      <<", D0="<<D<<endl;


  const string outfile=observablefilename(L,D,M,Z,directory);
  ofstream* out;

  out=new ofstream(outfile.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output (append="<<app<<")"<<endl;
    exit(1);
  }
  *out<<setprecision(15);

  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);
  SoftCoulombHamiltonian H(L,Delta,Z,a,M,Npen,Ntot,Spen);
  const MPO& hamil=H.getHMPO();
  // Observables I will compute
  // Compute the expectation value of N, as it commutes with H
  MPO Nmpo(L);
  H.getNumberMPO(Nmpo);
  MPO Smpo(L);
  H.getSpinMPO(Smpo);
  mwArray localN;
  H.getLocalNumberOp(localN);

  vector<MPS*> levels;

  // Read the computed levels
  bool oldMPS=true;
  int lval=0;
  MPS* nextLev;
  while(oldMPS&&(lval<=nLev)){
    const string filename=mpsfilename(L,D,M,Z,lval,mpsdir);
    cout<<"Checking for file ("<<lval<<") "<<filename<<endl;
    // Try to open file for MPS level lval
    if(file_exists(filename)){
      nextLev=new MPS(L,D,d);
      nextLev->importMPS(filename.data());
      levels.push_back(nextLev);
      cout<<"Recovered MPS file for level "<<lval<<" from file "<<filename<<endl;
      lval++;
    }
    else{
      oldMPS=false; // break the loop 
    }
  }
  if(lval<nLev){
    cout<<"Apparently not enough levels computed for the desired nLev("<<nLev
	<<")=> computing for "<<lval<<endl;
  }

  // Now I already have read all existing MPS
  cout<<"Recovered "<<lval<<" already computed MPS levels from disk. "<<endl;

  for(int l=0;l<lval;l++){
    cout<<"Computing for level "<<l<<endl;
    MPS* exc=levels[l];
    // values for each level
    complex_t Hval,H2val,normN,Nval,Sval;
    Hval=contractor.contract(*exc,hamil,*exc);
    //    H2val=contractor.contract2(hamil,*exc);
    normN=contractor.contract(*exc,*exc);
    Nval=contractor.contract(*exc,Nmpo,*exc);
    Sval=contractor.contract(*exc,Smpo,*exc);
    *out<<lval<<"\t"<<exc->getBond()<<"\t";
    *out<<real(normN)<<"\t";
    //*out<<real(H2val)<<"\t"; //<<imag(H2)<<"\t";
    *out<<0<<"\t"; //<<imag(H2)<<"\t";
    *out<<real(Hval)<<"\t"; //<<imag(Hval)<<"\t";
    *out<<real(Nval)<<"\t";
    *out<<real(Sval)<<"\t"; // densities start on column 8
    computeSingleBodyEV(*exc,localN,*out);
    *out<<endl;
    // I also wanted the cross-elements: <Psi_g|n_l|Psi_e>
    if(l>0){
      MPS* gs=levels[0];
      vector<complex_t> vals;
      computeCrossedSingleBodyEV(*gs,*exc,localN,vals);
      *out<<lval<<"\t"<<exc->getBond()<<"\t";
      *out<<real(normN)<<"\t";
      *out<<real(H2val)<<"\t"; //<<imag(H2)<<"\t";
      *out<<real(Hval)<<"\t"; //<<imag(Hval)<<"\t";
      *out<<real(Nval)<<"\t";
      *out<<real(Sval)<<"\t"; // densities start on column 8
      for(int p=0;p<L;p++){
	*out<<real(vals[p]/normN)<<"\t";
      }
      *out<<endl;
      *out<<lval<<"\t"<<exc->getBond()<<"\t";
      *out<<real(normN)<<"\t";
      *out<<real(H2val)<<"\t"; //<<imag(H2)<<"\t";
      *out<<real(Hval)<<"\t"; //<<imag(Hval)<<"\t";
      *out<<real(Nval)<<"\t";
      *out<<real(Sval)<<"\t"; // densities start on column 8
      for(int p=0;p<L;p++){
	*out<<imag(vals[p]/normN)<<"\t";
      }
      *out<<endl;
    }


  }


  
  out->close();
  delete out;

  // Clear the MPS vector
  while(levels.size()>1){
    delete levels.back();
    levels.pop_back();
  }  


   // for(int z=1;z<L;z++){
   //   delete Skmpos[z];
   // }
   // delete []Skmpos;
}

#include <sstream>    

const string mpsfilename(int L,int D,int M,int Ze,int level,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS_L"<<L<<"_M"<<M<<"_Z"<<Ze<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

const string observablefilename(int L,int D,int M,int Ze,const string outputdir){
  stringstream s;
  s<<outputdir<<"/density_L"<<L<<"_M"<<M<<"_Z"<<Ze<<"_D"<<D<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}


int computeSingleBodyEV(const MPS& exc,const mwArray& op1,ostream& os,double low,double high){
  bool result=0;
  bool checking=(high-low>0);
  int L=exc.getLength();
  MPO aux(L);
  mwArray basic_=identityMatrix(d*d); // basic piece of the MPO
  basic_.reshape(Indices(d*d,1,d*d,1));
  Operator basic(basic_);
  for(int k=0;k<L;k++){
    aux.setOp(k,&basic,0);
  }

  Contractor& contractor=Contractor::theContractor();
  double norm=real(contractor.contract(exc,exc));
  if(abs(norm-1)>1E-10){
    cout<<"** WARNING: State not normalized"<<endl;
  }

  mwArray op1_(op1);
  op1_.reshape(Indices(d*d,1,d*d,1));
  Operator op(op1_);

  // Now, site by site, replace the identity by the operator,
  // contract, write the result to os, and replace the operator by the
  // identity again
  for(int k=0;k<L;k++){
    aux.setOp(k,&op,0); // operator on site k
    complex_t valK=contractor.contract(exc,aux,exc);
    os<<real(valK)/norm<<"\t";
    aux.setOp(k,&basic,0); // operator removed
    if(checking&&(abs(real(valK)/norm)>high||abs(real(valK)/norm)<low)){
      cout<<"** WARNING: Non physical state detected, site "<<k<<" shows EV="<<real(valK)/norm
	  <<"+i*"<<imag(valK)/norm<<", norm="<<norm<<endl;
      result=-1;
    }
  }
  return result;
}

int computeCrossedSingleBodyEV(const MPS& gs,const MPS& exc,const mwArray& op1,vector<complex_t>& vals,double low,double high){
  bool result=0;
  bool checking=(high-low>0);
  int L=exc.getLength();
  MPO aux(L);
  mwArray basic_=identityMatrix(d*d); // basic piece of the MPO
  basic_.reshape(Indices(d*d,1,d*d,1));
  Operator basic(basic_);
  for(int k=0;k<L;k++){
    aux.setOp(k,&basic,0);
  }

  Contractor& contractor=Contractor::theContractor();
  double norm=real(contractor.contract(exc,exc));
  if(abs(norm-1)>1E-10){
    cout<<"** WARNING: State not normalized"<<endl;
  }

  mwArray op1_(op1);
  op1_.reshape(Indices(d*d,1,d*d,1));
  Operator op(op1_);

  // Now, site by site, replace the identity by the operator,
  // contract, write the result to os, and replace the operator by the
  // identity again
  for(int k=0;k<L;k++){
    aux.setOp(k,&op,0); // operator on site k
    complex_t valK=contractor.contract(exc,aux,gs);
    vals.push_back(valK);
    aux.setOp(k,&basic,0); // operator removed
    if(checking&&(abs(real(valK)/norm)>high||abs(real(valK)/norm)<low)){
      cout<<"** WARNING: Non physical state detected, site "<<k<<" shows EV="<<real(valK)/norm
	  <<"+i*"<<imag(valK)/norm<<", norm="<<norm<<endl;
      result=-1;
    }
  }
}
