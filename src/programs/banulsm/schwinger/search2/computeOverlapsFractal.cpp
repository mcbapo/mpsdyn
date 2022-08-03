
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
const string mpsfilename(int L,int D,double x,double mg,int level,const string& mpsdir);
/** 
Generate the name of the file with the command in it
*/
const string jobfilename(int L,int D,double x,double mg,double tol,const string& mpsdir);

void recoverMPS(MPS& state,int L,int D,double x,double mg,int nLev,const string& mpsdir,bool posSz=true);

complex_t computeOverlap(double& normP,const MPS& bra,const MPS& state,const vector<int>& pos,const vector<int>& lab,bool flipBra=false,bool flipKet=false);

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


/** This program should read the GSs computed by searchOrder and
    compute its energy and some overlaps to analyze its fractal
    structure. 

    In this version, I compute literally the elements they wrote in the
    paper. Without consideration for my even-odd exchange, since now it
    seems they start in 0!!!!

    Projecting the beginning of the chain.

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
      *out<<"% Notice Psi4 means Psi(N-4), Psi5+ means Psi(N-5,Sz=+1/2) etc"<<endl;
      *out<<"% The norm terms (header Nrhs) for each overlap are the norm of the ket of the next overlap, e.g. |P10 Psi2|^2"<<endl;
      *out<<"%lev\tL\tx\tmu\tD\t";
      *out<<"N1\treal(<Psi(4)|P01 Psi2>/sqrt(N1))\timag(<Psi(4)|P01 Psi2>/sqrt(N1))\t";
      *out<<"N2\treal(<Psi4|P10 Psi2>/sqrt(N2))\timag(<Psi4|P10 Psi2>/sqrt(N2))\t";
      *out<<"N3\treal(<Psi6|P0011 Psi2>/sqrt(N3))\timag(<Psi6|P0011 Psi2>/sqrt(N3))\t";
      *out<<"N4\treal(<Psi6|P01 Psi4>/sqrt(N4))\timag(<Psi6|P01 Psi4>/sqrt(N4))\t";
      *out<<"N5\treal(<Psi5+|P110 Psi2>/sqrt(N5))\timag(<Psi5+|P110 Psi2>/sqrt(N5))\t";
      *out<<"N6\treal(<Psi5+|P1 Psi4>/sqrt(N6))\timag(<Psi5+|P1 Psi4>/sqrt(N6))\t";
      *out<<"N7\treal(<Psi4|P0 Psi3+>/sqrt(N7))\timag(<Psi4|P0 Psi3+>/sqrt(N7))\t";
      *out<<"N8\treal(<Psi5+|P10 Psi3+>/sqrt(N8))\timag(<Psi5+|P10 Psi3+>/sqrt(N8))\t";
      *out<<"N9\treal(<Psi6|P1100 Psi2>/sqrt(N9))\timag(<Psi6|P1100 Psi2>/sqrt(N9))\t";
      *out<<"N10\treal(<T Psi5+|P001 Psi2>/sqrt(N10))\timag(<T Psi5+|P001 Psi2>/sqrt(N10))\t";
      *out<<"N11\treal(<Psi4|P1 T Psi3+>/sqrt(N11))\timag(<Psi4|P1 T Psi3+>/sqrt(N11))\t";
      *out<<"N12\treal(<Psi5+|P01 Psi3+>/sqrt(N12))\timag(<Psi5+|P01 Psi3+>/sqrt(N12))\t";
      *out<<"N13\treal(<Psi6|P10 Psi4>/sqrt(N13))\timag(<Psi6|P10 Psi4>/sqrt(N13))\t";
      *out<<"N14\treal(<T Psi5+|P0 Psi4>/sqrt(N14))\timag(<T Psi5+|P0 Psi4>/sqrt(N14))\t";
      *out<<endl;
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


  
  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);
  contractor.setExcOrthTol(1E-8); // max overlap between excitations
  int d=2;
  // First: put H a
  int Dop=5; // bond dimension of Schwinger H
  double lambda=0.;
  // SchwingerHamiltonianSz hSch(L,mu,x,L0,alpha,zpen,offset,nu,gweight);
  // //  SchwingerHamiltonian hSch(L,mu,x,L0,alpha,offset,nu,gweight);
  // const MPO& hamil=hSch.getHMPO();
  // // Observables I will compute
  // // Compute the expectation value of Sz, as it commutes with H
  // MPO Szmpo(L);
  // SpinMPO::getSzMPO(L,d,Szmpo);
  // // Sz2 is the penalty term (we need to substract it) 
  // MPO Sz2mpo(L);
  // SpinMPO::getSz2MPO(L,d,Sz2mpo);

  // Prepare the file name where I have to find the MPS

  bool done=0;
  while(!done){
    // read all the states
    MPS state2(L-2,D,d);
    MPS state4(L-4,D,d);
    MPS state6(L-6,D,d);
    MPS state3(L-3,D,d);MPS state3m(L-3,D,d);
    MPS state5(L-5,D,d);MPS state5m(L-5,D,d);
    recoverMPS(state2,L-2,D,x,mg,nLev,mpsdir);
    recoverMPS(state4,L-4,D,x,mg,nLev,mpsdir);
    recoverMPS(state6,L-6,D,x,mg,nLev,mpsdir);
    recoverMPS(state3,L-3,D,x,mg,nLev,mpsdir);
    recoverMPS(state5,L-5,D,x,mg,nLev,mpsdir);
    //recoverMPS(state3m,L-3,D,x,mg,nLev,mpsdir,false); // the ones with minus Sz
    //recoverMPS(state5m,L-5,D,x,mg,nLev,mpsdir,false);
    
    //    complex_t Ov42_10,Ov42_01,Ov62_0011,Ov64_10,Ov52_11I0,Ov54_I1,Ov43_I0,Ov53_01;
    *out<<nLev<<"\t"<<L<<"\t"<<x<<"\t"<<mu<<"\t"<<D<<"\t";
    // One by one the overlaps are computed and printed to file
    double normN;
    complex_t overl;
    vector<int> pos;vector<int> label;

    // 1:
    cout<<"<Psi(4)|P01 Psi2>"<<endl;
    pos.push_back(0);pos.push_back(1);
    label.push_back(0);label.push_back(1);
    overl=computeOverlap(normN,state4,state2,pos,label);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";

    // 2:
    cout<<"<Psi4|P10 Psi2>"<<endl;
    pos.clear();label.clear();
    pos.push_back(0);pos.push_back(1);
    label.push_back(1);label.push_back(0);
    overl=computeOverlap(normN,state4,state2,pos,label);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";

    // 3:
    cout<<"<Psi6|P0011 Psi2>"<<endl;
    pos.clear();label.clear();
    pos.push_back(0);pos.push_back(1);pos.push_back(2);pos.push_back(3);
    label.push_back(0);label.push_back(0);label.push_back(1);label.push_back(1);
    overl=computeOverlap(normN,state6,state2,pos,label);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";

    // 4:
    cout<<"<Psi6|P01 Psi4>"<<endl;
    pos.clear();label.clear();
    pos.push_back(0);pos.push_back(1);
    label.push_back(0);label.push_back(1);
    overl=computeOverlap(normN,state6,state4,pos,label);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";

    // 5:
    cout<<"<Psi5+|P110 Psi2>"<<endl;
    pos.clear();label.clear();
    pos.push_back(0);pos.push_back(1);pos.push_back(2);
    label.push_back(1);label.push_back(1);label.push_back(0);
    overl=computeOverlap(normN,state5,state2,pos,label);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";

    // 6:
    cout<<"<Psi5+|P1 Psi4>"<<endl;
    pos.clear();label.clear();
    pos.push_back(0);
    label.push_back(1);
    overl=computeOverlap(normN,state5,state4,pos,label);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";

    // 7:
    cout<<"<Psi4|P0 Psi3+>"<<endl;
    pos.clear();label.clear();
    pos.push_back(0);
    label.push_back(0);
    overl=computeOverlap(normN,state4,state3,pos,label);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";

    // 8:
    cout<<"<Psi5+|P10 Psi3+>"<<endl;
    pos.clear();label.clear();
    pos.push_back(0);pos.push_back(1);
    label.push_back(1);label.push_back(0);
    overl=computeOverlap(normN,state5,state3,pos,label);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";

    // 9:
    cout<<"<Psi6|1100 Psi2>"<<endl;
    pos.clear();label.clear();
    pos.push_back(0);pos.push_back(1);pos.push_back(2);pos.push_back(3);
    label.push_back(1);label.push_back(1);
    label.push_back(0);label.push_back(0);
    overl=computeOverlap(normN,state6,state2,pos,label);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";

    // 10:
    cout<<"<T Psi5|001 Psi2>"<<endl;
    pos.clear();label.clear();
    pos.push_back(0);pos.push_back(1);pos.push_back(2);
    label.push_back(0);label.push_back(0);label.push_back(1);
    //    pos.push_back(L-2-3);pos.push_back(L-2-2);pos.push_back(L-2-1);
    //label.push_back(1);label.push_back(0);label.push_back(0);
    overl=computeOverlap(normN,state5,state2,pos,label,true,false);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";

    // 11: 
    cout<<"<Psi4|1 T Psi3>"<<endl;
    pos.clear();label.clear();
    pos.push_back(0);
    label.push_back(1);
    overl=computeOverlap(normN,state4,state3,pos,label,false,true);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";

    // 12:
    cout<<"<Psi5|01 Psi3>"<<endl;
    pos.clear();label.clear();
    pos.push_back(0);pos.push_back(1);
    label.push_back(0);label.push_back(1);
    overl=computeOverlap(normN,state5,state3,pos,label);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";

    // 13: 
    cout<<"<Psi6|10 Psi4> "<<endl;
    pos.clear();label.clear();
    pos.push_back(0);pos.push_back(1);
    label.push_back(1);label.push_back(0);
    overl=computeOverlap(normN,state6,state4,pos,label);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";

    // 14:
    cout<<"<T Psi5|0 Psi4>"<<endl;
    pos.clear();label.clear();
    pos.push_back(0);
    label.push_back(0);
    overl=computeOverlap(normN,state5,state4,pos,label,true,false);
    *out<<normN<<"\t"<<real(overl)<<"\t"<<imag(overl)<<"\t";
    
    *out<<endl;
    D=getNextD(D);
    if(D<0) done=true;
  // just once to check!!
    if(done&&L==12){
      char* name=new char[100];
      sprintf(name,"mps_L%d_x%g_mu%g",L-6,x,mu);
      state6.exportForMatlab(name,15);
      cout<<"Written MPS to file "<<name<<endl;
      delete name;
    }
  }
    
  out->close();
  delete out;
  
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

const string mpsfilename(int L,int D,double x,double mg,int level,const string& mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS_L"<<L<<"_m"<<mg<<"_x"<<x<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

const string jobfilename(int L,int D,double x,double mg,double tol,const string& jobsdir){
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

void recoverMPS(MPS& state,int L,int D,double x,double mg,int nLev,const string& mpsdir,bool posSz){
  stringstream s;
  s<<mpsdir<<"/L"<<L<<"/mps";
  if(L%2!=0) s<<"/Stot"<<(posSz?.5:-.5);

  double mu=2*mg*sqrt(x);
  int d=2;
  
  const string filename=mpsfilename(L,D,x,mg,nLev,s.str());
  // Try to open file for MPS level nLev
  if(file_exists(filename)){
    state.importMPS(filename.data());
    cout<<"Recovered MPS file for level "<<nLev<<" from file "<<filename<<endl;
    Contractor& contractor=Contractor::theContractor();
    SchwingerHamiltonianSz hSch(L,mu,x,0.,0.,0.,0.,0.,1.);
    //  SchwingerHamiltonian hSch(L,mu,x,L0,alpha,offset,nu,gweight);
    const MPO& hamil=hSch.getHMPO();
    // Observables I will compute
    // Compute the expectation value of Sz, as it commutes with H
    MPO Szmpo(L);
    SpinMPO::getSzMPO(L,d,Szmpo);
    // Sz2 is the penalty term (we need to substract it) 
    // MPO Sz2mpo(L);
    //SpinMPO::getSz2MPO(L,d,Sz2mpo);
    cout<<"\tenergy: "<<contractor.contract(state,hamil,state)<<", Sz="<<contractor.contract(state,Szmpo,state)<<endl;
    state.gaugeCond('R',1);
    state.gaugeCond('L',1);
  }
  else{
    cout<<"ERROR: Could not find file with MPS ("<<L<<") "<<filename<<" for these parameters"<<endl;
    exit(1);
  }
}

complex_t computeOverlap(double& normP,const MPS& bra,const MPS& state,const vector<int>& pos,const vector<int>& lab,bool flipBra,bool flipKet){
  cout<<"computeOverlap  of bra:"<<bra.getLength()<<" and ket:"<<state.getLength()
      <<", with the latter projected as pos:"<<pos<<", label:"<<lab<<endl;
  complex_t data0[]={ONE_c,ZERO_c};
  complex_t data1[]={ZERO_c,ONE_c};
  mwArray p0(Indices(1,2),data0);//<v0|
  mwArray p1(Indices(1,2),data1);//<v1|
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(2,2),datax);
  // Make a copy of the long state (rhs)
  MPS aux(state);
  if(flipKet){
    for(int k=0;k<state.getLength();k++)
      aux.applyLocalOperator(k,sigX);
  }
  vector<int> projected(aux.getLength(),0);
  for(int k=0;k<pos.size();k++){
    if(lab[k]==0)
      aux.applyLocalOperator(pos[k],p0);
    else
      aux.applyLocalOperator(pos[k],p1);
    projected[pos[k]]=1; 
  }
  Contractor& contractor=Contractor::theContractor();
  normP=real(contractor.contract(aux,aux));
  cout<<"After projecting, norm of ket:"<<normP<<", projected:"<<projected<<endl;
  aux.gaugeCond('L',1); // now normalized
  if(flipBra){ // apply the flip now, after the projection, but only to sites not projected
    for(int k=0;k<aux.getLength();k++)
      if(aux.getA(k).getd()!=1)
	aux.applyLocalOperator(k,sigX);
  }
  // Now extend the length of the shorter MPS
  MPS aux2(aux.getLength(),aux.getBond(),2);
  int offset=0; // whenever finding a projected site, the positions are shifted by one
  for(int k=0;k<aux.getLength();k++){
    if(!projected[k]){
      aux2.replaceSite(k,bra.getA(k-offset),false);
      // cout<<"Replaced site "<<k<<" with tensor of pos "<<k-offset<<"in original bra"<<endl;
    }
    else{
      int Dbond=k==0?1:aux2.getA(k-1).getDr();
      aux2.replaceSite(k,reshape(identityMatrix(Dbond),Indices(1,Dbond,Dbond)),false);
      offset++;
    }
  }
  cout<<"Checking, after extending bra from "<<bra.getLength()<<" to "<<aux2<<"\t with norm "<<contractor.contract(aux2,aux2)<<" to be multiplied with "<<aux<<" of norm "<<contractor.contract(aux,aux)<<endl;
  cout<<"The overlap is: "<<contractor.contract(aux2,aux)<<endl;
  return contractor.contract(aux2,aux);
  
}


// void getProjectorMPO(MPO& proj,int L,const vector<int>& pos,const vector<int>& lab){
//   proj=MPO(L);
//   int d=2;
//   // Basic operators
//   mwArray sig0=identityMatrix(d);//identity
//   // complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
//   // mwArray sigX(Indices(d,d),datax);//sigmax
//   complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
//   mwArray sigZ(Indices(d,d),dataz);//sigmaz
    
//   Operator proj0(reshape(.5*(sig0+sigZ),Indices(d,1,d,1)));
//   Operator proj1(reshape(.5*(sig0-sigZ),Indices(d,1,d,1)));
  
//   MPO proj(L);
//   for(int k=0;k<L;k++)
//     proj.setOp(k,new Operator(reshape(sig0,Indices(d,1,d,1))),true);
//   for(int k=0;k<pos.size();k++){
//     if(lab[k]==0)
      
//     else
//       aux.applyLocalOperator(pos[k],p1);
//     projected[pos[k]]=1; 
//   }

// }
