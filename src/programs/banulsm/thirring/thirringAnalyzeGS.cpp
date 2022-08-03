
#include <math.h>
#include <iomanip>

#include "Properties.h"
#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "SpinMPO.h"
#include "misc.h"


#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "ThirringHamiltonian.h"

using namespace shrt;

/** thirringGS reads a MPS produced by thirringP, and analyzes
    expectation values. Arguments refer to the parameters of the GS to
    be analyzed, except the file name for output.

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
    \param <app> (int) whether the output file is to be kept (app==1)
*/


int findGSFile(int D,const string& mpsdir,const string& baseName,string& mpsfile);//,int level=0);
const string mpsFileName(const string& mpsdir,const string& baseName,int D,int level=0);
const string baseFileName(const Properties& props);
void values(const MPS& MPS,ofstream* out);
complex_t valZ(const MPS& gs,int pos);

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
    cout<<"ERROR: Nr of site L should be even"<<endl;
    exit(1);
  }
  double ma=props.getDoubleProperty("ma");  
  double gt2_=props.getDoubleProperty("gtilde2");
  double g2_=props.getDoubleProperty("g2"); // the non-alternating one (Delta)
  double lambda_=props.getDoubleProperty("lambda");
  if(lambda_<0) lambda_=0.;
  double mu_=props.getDoubleProperty("mu");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  int Stot_=props.getIntProperty("Starget");
  int D=props.getIntProperty("D");

  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir");
  //  string mpsfile=props.getProperty("mpsfile");
  //  string mpstostretch=props.getProperty("stretchmps");
  bool app=(props.getIntProperty("append")!=0);

  cout<<"Initialized arguments: L="<<L
      <<", ma="<<ma
      <<", g2="<<g2_
      <<", gt2="<<gt2_
      <<", lambda="<<lambda_
      <<", mu="<<mu_
      <<", Starget="<<Stot_
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  cout<<"Initialized Contractor"<<endl;
  contractor.setConvTol(tol);
  int d=2;
  MPS gs(L,D,d);
    // Check if initial state was given
  string mpsfile=baseFileName(props); // base of file name containing these parameters

  string tmpMPSfile;
  if(findGSFile(D,mpsdir,mpsfile,tmpMPSfile)){// found EXACTLY this file
    cout<<"Recovered initial guess from file "<<tmpMPSfile<<endl;
    if(gs.getLength()!=L){
      cout<<"Error! File "<<tmpMPSfile<<" does not contain valid MPS "<<endl;
      exit(1);
    }
    gs.importMPS(tmpMPSfile.data());
  }
  else{
    cout<<"ERROR: No GS found"<<endl;
    exit(1);
    //    gs.exportMPStext("randomMPS.txt");
  }
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  cout<<setprecision(10);
  
  MPO Szmpo(L);
  SpinMPO::getSzMPO(L,d,Szmpo);
 
  ThirringHamiltonian hamH(L,ma,g2_,lambda_,Stot_,d,gt2_,mu_);
  //cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamH.getHMPO();

  complex_t E0=contractor.contract(gs,hamil,gs);
  complex_t Sztot=contractor.contract(gs,Szmpo,gs);
  cout<<"Energy contraction is:"<<E0<<" and Sz contraction "<<Sztot<<endl;
  

  cout<<"Ground state found for S_tot="<<Stot_<<" with energy "<<E0
      <<" (Energy per particle="<<E0/L<<", plus offset from David->"
      <<real(E0)/L+gt2_*.25+g2_*(L-1.)/(4.*L)<<")"
      <<" Total Sz="<<real(Sztot)<<endl;

  complex_t Sz2=contractor.contract2(Szmpo,gs);
  complex_t H2=contractor.contract2(hamil,gs);
  cout<<"DeltaH="<<real(H2)-real(E0*E0)<<endl;
  cout<<"Delta S_z="<<Sz2-Sztot*Sztot<<endl;

  // cout<<"Computing pieces of H "<<endl;
  // MPO Hterm(L);
  // SpinMPO::getSxSxMPO(L,d,Hterm);
  // complex_t SxSxtot=contractor.contract(gs,Hterm,gs);
  // cout<<"Obtained <SxSx>="<<SxSxtot<<endl;
  // SpinMPO::getSySyMPO(L,d,Hterm);
  // complex_t SySytot=contractor.contract(gs,Hterm,gs);
  // cout<<"Obtained <SySy>="<<SySytot<<endl;
  // SpinMPO::getSzSzMPO(L,d,Hterm);
  // complex_t SzSztot=contractor.contract(gs,Hterm,gs);
  // cout<<"Obtained <SzSz>="<<SzSztot<<endl;
  // SpinMPO::getSzMPO(L,d,Hterm,true);
  // complex_t SzStaggered=contractor.contract(gs,Hterm,gs);
  // cout<<"Obtained <(-1)^n Sz>="<<SzStaggered<<endl;

  mwArray sig0=identityMatrix(d);
  complex_t dataX[]={ZERO_c,.5*ONE_c,.5*ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  complex_t dataY[]={ZERO_c,.5*I_c,-.5*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);


  
  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    *out<<"% GS analysys for N="<<L<<", ma="<<ma<<", gt2="<<gt2_<<", g2="<<g2_<<", Stot="<<Stot_<<endl;
    *out<<"% With:\n%\t E="<<real(E0)<<", DeltaH^2="<<real(H2)-real(E0*E0)<<endl;
    *out<<"%\tSz="<<real(Sztot)<<", DeltaSz="<<real(Sz2-Sztot*Sztot)<<endl;
    *out<<"%pos1\tpos2\t<XX>\t<YY>\t<ZZ>\t<Z1>\t<Z2>\tS(pos1:rest)"<<endl;
  }
  else
    out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(15);
  values(gs,out);


  out->close();
  delete out;
}

const string baseFileName(const Properties& props){
  int L=props.getIntProperty("L");
  double ma=props.getDoubleProperty("ma");  
  double gt2_=props.getDoubleProperty("gtilde2");
  double g2_=props.getDoubleProperty("g2"); // the non-alternating one
  double lambda_=props.getDoubleProperty("lambda");
  if(lambda_<0) lambda_=0.;
  double mu_=props.getDoubleProperty("mu");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  int Stot_=props.getIntProperty("Starget");
  stringstream s;
  s<<"MPSthirring_L"<<L<<"_ma"<<ma<<"_gt2"<<gt2_<<"_g2"<<g2_<<"_mu"<<mu_<<"_Stot"<<Stot_;
  return s.str();
}

const string mpsFileName(const string& mpsdir,const string& baseName,int D,int level){
  stringstream s;
  s<<mpsdir<<"/"<<baseName;
  if(level>0) s<<"_"<<level;
  if(D>0) s<<"_"<<D;
  return s.str();
}

int findGSFile(int D,const string& mpsdir,const string& baseName,string& mpsfile){
  mpsfile=mpsFileName(mpsdir,baseName,D,0);
  bool found=0;
  if(file_exists(mpsfile)){
    cout<<"findGSFile found mpsfile="<<mpsfile<<endl;
    found=true;
    return 1;
  }
  else{
    cout<<"File for GS "<<mpsfile<<" not found!"<<endl;
    return 0;
  }
  cout<<"ERROR: Should not be here"<<endl;
  exit(1);
}

complex_t valZ(const MPS& gs,int pos){
  int d=2;
  complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  Contractor& contractor=Contractor::theContractor();
  MPS aux(gs);
  aux.applyLocalOperator(pos,sigZ,0);
  return contractor.contract(aux,gs);  
}


void values(const MPS& gs,ofstream* out){
  //void check(const MPS& gs,double gt2,double g2,double ma,double mu,const MPO& hamil,double lambda,double Stot){
  int d=2;
  mwArray sig0=identityMatrix(d);
  complex_t dataX[]={ZERO_c,.5*ONE_c,.5*ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  complex_t dataY[]={ZERO_c,.5*I_c,-.5*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);


  vector<complex_t> valsXX;
  vector<complex_t> valsYY;
  vector<complex_t> valsZZ;
  vector<complex_t> valsZ;

  // Notice this is very innefficient: should use canonical forms and reuse contractions as much as possible
  
  Contractor& contractor=Contractor::theContractor();
  int L=gs.getLength();
  MPS aux(gs);
  for(int i=0;i<L-1;i++){

    complex_t XX,YY,ZZ,Z1,Z2;
    int pos1=i;int pos2=i+1; // should loop the second one
    aux.applyLocalOperator(i,sigX,0);
    aux.applyLocalOperator(i+1,sigX,0);
    XX=contractor.contract(aux,gs);
    
    aux.applyLocalOperator(i,4*sigY*sigX,0);
    aux.applyLocalOperator(i+1,4*sigY*sigX,0);
    YY=contractor.contract(aux,gs);
    
    aux.applyLocalOperator(i,4*sigZ*sigY,0);
    aux.applyLocalOperator(i+1,4*sigZ*sigY,0);
    ZZ=contractor.contract(aux,gs);
    
    aux.applyLocalOperator(i+1,4*sigZ,0);
    Z1=contractor.contract(aux,gs);

    aux.applyLocalOperator(i,4*sigZ,0);
    aux.applyLocalOperator(i+1,sigZ,0);
    Z2=contractor.contract(aux,gs);
    aux.applyLocalOperator(i+1,4*sigZ,0);

    double Spos=0.;
    if(i>0)Spos=contractor.getEntropy(gs,i);

    *out<<pos1<<"\t"<<pos2<<"\t"<<real(XX)<<"\t"<<real(YY)<<"\t"<<real(ZZ)<<"\t"<<real(Z1)<<"\t"<<real(Z2)<<"\t"<<Spos<<endl;
  }
  
}
