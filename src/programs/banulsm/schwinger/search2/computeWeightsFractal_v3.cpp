
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"
//#include "HermitianTensorMultiplier.h"

#include "SchwingerHamiltonianSz.h"
#include "SpinMPO.h"
#include <cmath>
#include <unistd.h>

#include "misc.h"
#include "quicksort.h"

using namespace std;
using namespace shrt;

#define MAXLEN 120

/** 
Generates a file name for the MPS corresponding to a certain excitation.
*/
const string mpsfilename(int L,int D,double x,double mg,int level,const string mpsdir);
/** 
Generate the name of the file with the command in it
*/
const string jobfilename(int L,int D,double x,double mg,double tol,const string mpsdir);

// /** Check if the filename exists */
// bool file_exists(const string filename);

int getNextD(int D);

/** Auxiliary object, level-energy, to sort */
class Level{
public:
  int lev;
  double energy;
  Level(int l, double E):lev(l),energy(E){};
  Level(const Level& l2):lev(l2.lev),energy(l2.energy){};
  bool operator<=(const Level& l2){return energy<=l2.energy;}
  bool operator<(const Level& l2){return energy<l2.energy;}
  bool operator==(const Level& l2){return energy==l2.energy;}
  friend bool operator<=(const Level& l1,const Level& l2){return l1.energy<=l2.energy;}
  friend bool operator<(const Level& l1,const Level& l2){return l1.energy<l2.energy;}
  friend bool operator==(const Level& l1,const Level& l2){return l1.energy==l2.energy;}
};


/** This program should read the GS computed by searchOrder and
    compute its energy and some weights to analyze its fractal
    structure. 

    In this case, project the end of the chain. But what I should do is to take into account that in my convention, sites are ordered as even-odd, not odd-even as in Hamer's

    Receives as argument a Properties file, such as
    config/excitProp.conf, from where the required parameters are
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

  int L = props.getIntProperty("L");
  double mg=props.getDoubleProperty("mg");
  double alpha=props.getDoubleProperty("alpha");
  double x=props.getDoubleProperty("x");
  int nLev = props.getIntProperty("nLev");
  directory=props.getProperty("outputdir");
  const string outfile=directory+"/"+props.getProperty("outfile");
  const string outfileHeff=directory+"/"+props.getProperty("outfileHeff");
  int app=props.getIntProperty("append");
  if(app<0) app=1; // default will be to append
  double tol=props.getDoubleProperty("tol");
  int D=props.getIntProperty("D");
  double zpen=0.; //props.getDoubleProperty("penalty");
  double offset=0.; //props.getDoubleProperty("offset");
  string mpsdir=props.getProperty("mpsdir"); // default current
  if(mpsdir.empty()) 
    mpsdir=".";

  bool contD=props.getIntProperty("singleD")!=1; // only if -singleD=1 given, it will stop at the given value of D

  
  // Parameters I do not use here: chemical potential and weight of the gauge term
  double nu=0;
  // I might want to switch off the gauge terms
  //double gweight=props.getDoubleProperty("gweight");
  //if(gweight==-1) gweight=1; // default: do not alter
  double gweight=1.; // Not touching the gauge terms (std)

  double mu=2*mg*sqrt(x);
  double L0=0.; // l0 parameter, which is irrelevant, as only appears
		// as alpha+l0

  //  double t0=.01; // todo: differently -> reference time

  cout<<"Initialized arguments: L="<<L
      <<", mu="<<mu
      <<", x="<<x
      <<", alpha="<<alpha
      <<", outfile="<<outfile
      <<", app="<<app
      <<", tol="<<tol
      <<", level="<<nLev
      <<", D="<<D<<endl;

  ofstream* out;

  // Now I will project the last spins
  if(!app){
      out=new ofstream(outfile.data());
      *out<<"%lev\tL\tx\tmu\tD\t<Psi|Psi>\t<H>\t<Sz>\t<H^2>\t";
      *out<<"P000\tP001\tP01\tP10\tP110\tP111\t";
      *out<<"P10_000\tP10_001\tP10_01\tP10_10\tP10_110\tP10_111\t";
      *out<<"P01_000\tP01_001\tP01_01\tP01_10\tP01_110\tP01_111"<<endl;
  }
  else{
    out=new ofstream(outfile.data(),ios::app);
  }
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output (append="<<app<<")"<<endl;
    exit(1);
  }
  *out<<setprecision(15);

  ofstream* outH;
  if(!app){
      outH=new ofstream(outfileHeff.data());
  }
  else{
    outH=new ofstream(outfileHeff.data(),ios::app);
  }
  if(!outH->is_open()){
    cout<<"Error: impossible to open file "<<outfileHeff<<
      " for output (append="<<app<<")"<<endl;
    exit(1);
  }
  *outH<<setprecision(10);

  
  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);
  contractor.setExcOrthTol(1E-8); // max overlap between excitations
  int d=2;
  // First: put H a
  int Dop=5; // bond dimension of Schwinger H
  double lambda=0.;
  SchwingerHamiltonianSz hSch(L,mu,x,L0,alpha,zpen,offset,nu,gweight);
  //  SchwingerHamiltonian hSch(L,mu,x,L0,alpha,offset,nu,gweight);
  const MPO& hamil=hSch.getHMPO();
  // Observables I will compute
  // Compute the expectation value of Sz, as it commutes with H
  MPO Szmpo(L);
  SpinMPO::getSzMPO(L,d,Szmpo);
  // Sz2 is the penalty term (we need to substract it) 
  MPO Sz2mpo(L);
  SpinMPO::getSz2MPO(L,d,Sz2mpo);

  // Prepare the file name where I have to find the MPS

  bool done=0;
  while(!done){
  
    MPS state(L,D,d);
    const string filename=mpsfilename(L,D,x,mg,nLev,mpsdir);
    // Try to open file for MPS level nLev
    if(file_exists(filename)){
      state.importMPS(filename.data());
      cout<<"Recovered MPS file for level "<<nLev<<" from file "<<filename<<endl;
      state.gaugeCond('R',1);
      state.gaugeCond('L',1);
    }
    else{
      cout<<"ERROR: Could not find file with MPS "<<filename<<" for these parameters"<<endl;
      exit(1);
    }

    complex_t Hval,H2val,normN,Szval,Sz2val,P000,P001,P01,P10,P110,P111;

    Hval=contractor.contract(state,hamil,state);
    H2val=contractor.contract2(hamil,state);
    normN=contractor.contract(state,state);

    // cout<<"Read state from "<<filename<<", norm: "<<normN<<" <H^2>-<H>^2="<<H2val-Hval*Hval<<endl;
    // if(real(H2val-Hval*Hval)<0){
    //   cout<<"WARNING: Variance negative:"<<real(H2val-Hval*Hval)<<" but relative value "<<real(H2val-Hval*Hval)/real(Hval*Hval)<<" and imaginary part of energy "<<imag(Hval)/real(Hval)<<endl;
    //   exit(1);
    // }
    
    Szval=contractor.contract(state,Szmpo,state);
    Sz2val=contractor.contract(state,Sz2mpo,state);
    
    
    // Basic operators
    mwArray sig0=identityMatrix(d);//identity
    // complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    // mwArray sigX(Indices(d,d),datax);//sigmax
    complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
    mwArray sigZ(Indices(d,d),dataz);//sigmaz
    
    Operator proj0(reshape(.5*(sig0+sigZ),Indices(d,1,d,1)));
    Operator proj1(reshape(.5*(sig0-sigZ),Indices(d,1,d,1)));
    Operator idOp(reshape(sig0,Indices(d,1,d,1)));
    
    MPO proj(L);
    for(int k=0;k<L;k++)
      proj.setOp(k,&idOp,false);
    int p0(L-1),p1(L-2),p2(L-3); // the last three, i.e. start at the end
    // 000
    proj.setOp(p0,&proj0,false); //0
    proj.setOp(p1,&proj0,false); //00
    proj.setOp(p2,&proj0,false); //000
    P000=contractor.contract(state,proj,state);
    proj.setOp(p2,&proj1,false); //001
    P001=contractor.contract(state,proj,state);
    proj.setOp(p2,&idOp,false); //00
    proj.setOp(p1,&proj1,false); //01
    P01=contractor.contract(state,proj,state);
    proj.setOp(p1,&proj0,false); //00
    proj.setOp(p0,&proj1,false); //10
    P10=contractor.contract(state,proj,state);
    proj.setOp(p1,&proj1,false); //11
    proj.setOp(p2,&proj0,false); //110
    P110=contractor.contract(state,proj,state);
    proj.setOp(p2,&proj1,false); //111
    P111=contractor.contract(state,proj,state);
    
    // But also inside the 3rd block:
    int p3(L-4),p4(L-5);
    proj.setOp(p0,&proj1,false); 
    proj.setOp(p1,&proj0,false); //10
    proj.setOp(p2,&proj0,false);proj.setOp(p3,&proj0,false);proj.setOp(p4,&proj0,false); // 000
    complex_t P10_000=contractor.contract(state,proj,state);
    proj.setOp(p2,&proj0,false);proj.setOp(p3,&proj0,false);proj.setOp(p4,&proj1,false); // 001
    complex_t P10_001=contractor.contract(state,proj,state);
    proj.setOp(p2,&proj0,false);proj.setOp(p3,&proj1,false);proj.setOp(p4,&idOp,false); // 01
    complex_t P10_01=contractor.contract(state,proj,state);
    proj.setOp(p2,&proj1,false);proj.setOp(p3,&proj0,false);proj.setOp(p4,&idOp,false); // 10
    complex_t P10_10=contractor.contract(state,proj,state);
    proj.setOp(p2,&proj1,false);proj.setOp(p3,&proj1,false);proj.setOp(p4,&proj0,false); // 110
    complex_t P10_110=contractor.contract(state,proj,state);
    proj.setOp(p2,&proj1,false);proj.setOp(p3,&proj1,false);proj.setOp(p4,&proj1,false); // 111
    complex_t P10_111=contractor.contract(state,proj,state);
  
    cout<<"Checking sumof projections:"<<P10_000+P10_001+P10_01+P10_10+P10_110+P10_111<<" vs P10="<<P10<<endl;

    // But also inside the 2nd block:
    proj.setOp(p0,&proj0,false); 
    proj.setOp(p1,&proj1,false); //01
    proj.setOp(p2,&proj0,false);proj.setOp(p3,&proj0,false);proj.setOp(p4,&proj0,false); // 000
    complex_t P01_000=contractor.contract(state,proj,state);
    proj.setOp(p2,&proj0,false);proj.setOp(p3,&proj0,false);proj.setOp(p4,&proj1,false); // 001
    complex_t P01_001=contractor.contract(state,proj,state);
    proj.setOp(p2,&proj0,false);proj.setOp(p3,&proj1,false);proj.setOp(p4,&idOp,false); // 01
    complex_t P01_01=contractor.contract(state,proj,state);
    proj.setOp(p2,&proj1,false);proj.setOp(p3,&proj0,false);proj.setOp(p4,&idOp,false); // 10
    complex_t P01_10=contractor.contract(state,proj,state);
    proj.setOp(p2,&proj1,false);proj.setOp(p3,&proj1,false);proj.setOp(p4,&proj0,false); // 110
    complex_t P01_110=contractor.contract(state,proj,state);
    proj.setOp(p2,&proj1,false);proj.setOp(p3,&proj1,false);proj.setOp(p4,&proj1,false); // 111
    complex_t P01_111=contractor.contract(state,proj,state);
  
    cout<<"Checking sum of projections:"<<P01_000+P01_001+P01_01+P01_10+P01_110+P01_111<<" vs P01="<<P01<<endl;

    // And the efffective Hamiltonian (not very efficient, as I create four states, but it is ok here)
    MPS psi1(state),psi2(state),psi3(state),psi4(state);
    mwArray pi0=.5*(sig0+sigZ);mwArray pi1=.5*(sig0-sigZ);
    // Old version, now all 00 together
    //    psi1.applyLocalOperator(p0,pi0,true);psi1.applyLocalOperator(p1,pi0,true);psi1.applyLocalOperator(p2,pi1,true);
    psi1.applyLocalOperator(p0,pi0,true);psi1.applyLocalOperator(p1,pi0,true);
    psi2.applyLocalOperator(p0,pi0,true);psi2.applyLocalOperator(p1,pi1,true);
    psi3.applyLocalOperator(p0,pi1,true);psi3.applyLocalOperator(p1,pi0,true);
    //    psi4.applyLocalOperator(p0,pi1,true);psi4.applyLocalOperator(p1,pi1,true);psi4.applyLocalOperator(p2,pi0,true);
    psi4.applyLocalOperator(p0,pi1,true);psi4.applyLocalOperator(p1,pi1,true);

    const MPS* mpss[]={&psi1,&psi2,&psi3,&psi4};
    mwArray Heff(Indices(4,4));
    *outH<<nLev<<"\t"<<L<<"\t"<<x<<"\t"<<mu<<"\t"<<D<<"\t";
    for(int k=0;k<4;k++)
      for(int p=0;p<4;p++){
	Heff.setElement(contractor.contract(*mpss[k],hamil,*mpss[p]),Indices(k,p));
	*outH<<Heff.getElement(Indices(k,p))<<"\t";
      }
    *outH<<endl;
    putForMatlab(cout,Heff,"Heff");
    //cout<<Heff<<endl;
    
    // write out results to the file
    *out<<nLev<<"\t"<<L<<"\t"<<x<<"\t"<<mu<<"\t"<<D<<"\t";
    *out<<real(normN)<<"\t";
    *out<<real(Hval)<<"\t"; //<<imag(Hval)<<"\t";
    *out<<real(Szval)<<"\t";
    *out<<real(H2val)<<"\t"; //<<imag(H2)<<"\t";
    *out<<real(P000)<<"\t"; 
    *out<<real(P001)<<"\t"; 
    *out<<real(P01)<<"\t"; 
    *out<<real(P10)<<"\t"; 
    *out<<real(P110)<<"\t"; 
    *out<<real(P111)<<"\t"; 
    *out<<real(P10_000)<<"\t"; 
    *out<<real(P10_001)<<"\t"; 
    *out<<real(P10_01)<<"\t"; 
    *out<<real(P10_10)<<"\t"; 
    *out<<real(P10_110)<<"\t"; 
    *out<<real(P10_111)<<"\t"; 
    *out<<real(P01_000)<<"\t"; 
    *out<<real(P01_001)<<"\t"; 
    *out<<real(P01_01)<<"\t"; 
    *out<<real(P01_10)<<"\t"; 
    *out<<real(P01_110)<<"\t"; 
    *out<<real(P01_111)<<"\t"; 
    *out<<endl;

    D=getNextD(D);
    if(D<0) done=true;
  }
    
  out->close();
  delete out;
  outH->close();delete outH;
  
} // End of main

#include <sstream>    

int getNextD(int D){
  if(D>=20){
    if(D<160)
      return D+20;
    else
      return -1;
  }
  else{
    if(D==10) return 20;
    if(D==4) return 10;
    if(D==2) return 4;  
    return -1;
  }
}

const string mpsfilename(int L,int D,double x,double mg,int level,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS_L"<<L<<"_m"<<mg<<"_x"<<x<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

const string jobfilename(int L,int D,double x,double mg,double tol,const string jobsdir){
  stringstream s;
  s<<jobsdir<<"/job_L"<<L<<"_m"<<mg<<"_x"<<x<<"_D"<<D;
  if(tol!=1E-7) s<<"_tol"<<-log10(tol);
  s<<".sh";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

// bool file_exists(const string filename)
// {
//   ifstream ifile(filename.data());
//   //cout<<"file_exists("<<filename<<")="<<ifile<<endl;
//   return ifile;
// }
