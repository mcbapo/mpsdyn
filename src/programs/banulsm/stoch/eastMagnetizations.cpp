
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#include "SpinMPO.h"
#include "StochasticHamiltonian.h"
#include "Properties.h"
#include "misc.h"

#include "quicksort.h"

#define FREQSAVE 18000 // frequency of tmp saving: default 5 hrs 


/** 
Generates a file name for the MPS corresponding to a certain excitation.
*/
const string mpsfilename(int L,int D,double s,double c,int level,const string mpsdir);
const string jobfilename(int L,int D,double s,double c,double q,const string jobsdir);

/** Compute individual occupations of sites */
void computeOccupations(const MPS& gs,int L,int d,vector<double>& occ);

/** Same for polarizations */
void computePolarizations(const MPS& gs,int L,int d,vector<double>& occ);


/** eastMagnetizations takes the states computed by basicSpec
    and computes observables on them.  In particular, the modulated (with
    a phase k) magnetizations in Z and X

    Receives as argument a Properties file, such as
    config/stochspec.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.

    Receives, among others, arguments:
    \param <L> (int) length of the chain
    \param <s> (double) parameter \f$s\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)

*/

using namespace shrt;

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
  double s=props.getDoubleProperty("s");
  double q=props.getDoubleProperty("q");
  int D=props.getIntProperty("D");
  int maxD=props.getIntProperty("maxD");
  int incrD=props.getIntProperty("incrD");
  if(incrD<0) incrD=20;
  double c=props.getDoubleProperty("c");
  if(c<0||c>1){
    cout<<"Error: value of c not allowed (c in[0,1])"<<endl;
    exit(1);
  }
  int nlev=props.getIntProperty("nlevel");
  string outfname=props.getProperty("outputfile");
  string mpsdir=props.getProperty("mpsdir");
  string jobsdir=props.getProperty("jobsdir");
  int append=props.getIntProperty("append");
  bool app=(append!=0);
  bool sectorOne=(props.getIntProperty("firstN")!=0);
  cout<<"Occupation of first site fixed to "<<(sectorOne?1:0)<<endl;

  cout<<"Appending to output?"<<app<<endl;  

  bool origApp=app;
  ofstream* out;
  if(!app||!file_exists(outfname)){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% Modulated polarization of stochastic model L="<<L<<", D="<<D<<endl;
    *out<<"% s="<<s<<endl;
    *out<<"% q="<<q<<endl;
    *out<<"% nlev\t D\t q\t<H>\t<H^2>\t <Z>\t <Z^2>\t <X>\t<X^2>\t<Z(q)^dag Z(q)>\t <X(q)^dag X(q)>"<<endl;
    out->close();delete out;
  }
  
  //out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
  //out->close();delete out;
    
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;

  int d=2;

  double offset=0.;
  StochasticHamiltonian hSto(L-1,s,offset,sectorOne,c);
  const MPO& hamil=hSto.getHMPO();
  //hamil.exportForMatlab("testSto.m");//exit(1);

  // Operators I want to compute:
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray sig0=identityMatrix(d);
  
  MPO SzMPO(L);
  MPO SxMPO(L);
  MPO SzqMPO(L);
  MPO SxqMPO(L);
  SpinMPO::getSingleBodyMPO(L,sig0,sigX,SxMPO);
  SpinMPO::getSingleBodyMPO(L,sig0,sigZ,SzMPO);
  SpinMPO::getModulatedSingleBodyMPO(L,sig0,sigZ,q,SzqMPO);
  SpinMPO::getModulatedSingleBodyMPO(L,sig0,sigX,q,SxqMPO);
  
  // REad the existing computed levels
  bool oldMPS=true;
  int lval=0;

  while(oldMPS&&lval<=nlev-1){
    const string filename=mpsfilename(L,D,s,c,lval,mpsdir);
    // Try to open file for MPS level lval
    if(file_exists(filename)){
      MPS nextLev(L-1,D,d);
      nextLev.importMPS(filename.data());
      cout<<"Recovered MPS file for level "<<lval<<" from file "<<filename<<endl;
      nextLev.gaugeCond('R',1);
      nextLev.gaugeCond('L',1);
      complex_t Ek=contractor.contract(nextLev,hamil,nextLev);
      complex_t Ek2=contractor.contract2(hamil,nextLev);
      complex_t Sz,Sz2,Sx,Sx2,Szq,Sxq,Szq2,Sxq2;
      {// but since first site is fixed, I need a trick for the staggered magnetizations

	MPS nextLev_(L,1,d);
	if(!sectorOne) nextLev_.setProductState(p_zero);
	else nextLev_.setProductState(p_one);
	for(int pos=1;pos<L;pos++)
	  nextLev_.replaceSite(pos,nextLev.getA(pos-1).getA(),false);
	
	Sz=contractor.contract(nextLev_,SzMPO,nextLev_);
	Sz2=contractor.contract2(SzMPO,nextLev_);
	Sx=contractor.contract(nextLev_,SxMPO,nextLev_);
	Sx2=contractor.contract2(SxMPO,nextLev_);
	Szq2=contractor.contract2(SzqMPO,nextLev_);
	Sxq2=contractor.contract2(SxqMPO,nextLev_);
	Szq=contractor.contract(nextLev_,SzqMPO,nextLev_);
	Sxq=contractor.contract(nextLev_,SxqMPO,nextLev_);

      }

      
      // Compute magnetizations
      
      out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
      *out<<lval<<"\t"<<D<<"\t"<<q<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t";
      *out<<real(Sz)<<"\t"<<real(Sz2)<<"\t";
      *out<<real(Sx)<<"\t"<<real(Sx2)<<"\t";
      *out<<real(Szq2)<<"\t"<<real(Sxq2)<<"\t";
      *out<<real(Szq)<<"\t"<<imag(Szq)<<"\t";
      *out<<real(Sxq)<<"\t"<<imag(Sxq)<<"\t";
      *out<<endl;
      out->close();delete out;

      lval++;
    }
    else{
      oldMPS=false; // break the loop 
    }
  }
  

  // Prepare the script with the new job
  if(jobsdir.length()>0&&D<maxD){
    int newD=D+incrD;
    string newOutput(outfname);
    stringstream newSuff;newSuff<<"_D"<<newD; 
    newOutput.replace(newOutput.find("_D"),string::npos,newSuff.str());
    stringstream commandstr;
    commandstr<<argv[0]; // exec
    commandstr<<" "<<argv[1]; // config file

    commandstr<<" -L="<<L<<" -s="<<s<<" -q="<<q;
    commandstr<<" -D="<<newD;
    commandstr<<" -maxD="<<maxD;
    commandstr<<" -incrD="<<incrD;
    commandstr<<" -c="<<c<<" -nlevel="<<nlev;
    commandstr<<" -outputfile="<<newOutput;
    commandstr<<" -mpsdir="<<mpsdir<<" -initmpsdir="<<mpsdir;
    commandstr<<" -jobsdir="<<jobsdir;
    commandstr<<" -append="<<origApp;
    commandstr<<" -firstN="<<sectorOne;
    string jobfile=jobfilename(L,newD,s,c,q,jobsdir);

    ofstream* ojob=new ofstream(jobfile.data());
    if(!ojob->is_open()){
      cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
      cout<<commandstr.str()<<endl;
      cout<<endl;
    }
    else{
      int nPar=4;
      if(newD>=100) nPar=10;
      if(newD>=200) nPar=10;

      *ojob<<"msub_slurm -N magQ.c"<<c<<".N"<<L<<".D"<<newD<<".s"<<s<<".q"<<q;
      *ojob<<" -P "<<nPar;
      *ojob<<" -- ";
      *ojob<<commandstr.str()<<endl;
      ojob->close();
    }
    delete ojob;

  }

  
}


const string mpsfilename(int L,int D,double s,double c,int level,const string mpsdir){
  stringstream str;
  str<<mpsdir<<"/MPS_L"<<L<<"_s"<<s<<"_c"<<c<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return str.str();
}

const string jobfilename(int L,int D,double s,double c,double q,const string jobsdir){
  stringstream str;
  str<<jobsdir<<"/job_mag_L"<<L<<"_s"<<s<<"_c"<<c<<"_q"<<q<<"_D"<<D;
  str<<".dat";
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return str.str();
}

void computeOccupations(const MPS& gs,int L,int d,vector<double>& occ){
  // Construct a MPO for the particle number on a given site
  static MPO mpoId(L); // Basic MPO: all Identity
  static Operator idOp(reshape(identityMatrix(d),Indices(d,1,d,1)));
  static Operator* nOp; // individual operator number
  static bool init(false);
  if(!init){
    mwArray nSite(Indices(d,d));
    for(int k=0;k<d;k++){
      nSite.setElement(k*ONE_c,Indices(k,k));
    }
    nSite.reshape(Indices(d,1,d,1));
    nOp=new Operator(nSite);
    for(int k=0;k<L;k++)
      mpoId.setOp(k,&idOp,0);
    init=true;
  }
  Contractor& contractor=Contractor::theContractor();
  occ.clear();
  for(int pos=0;pos<L;pos++){
    mpoId.setOp(pos,nOp,false);
    complex_t nk=contractor.contract(gs,mpoId,gs);
    mpoId.setOp(pos,&idOp,false);
    occ.push_back(real(nk));
  }
}

void computePolarizations(const MPS& gs,int L,int d,vector<double>& polX){
  // Construct a MPO for the sigma_x on a given site
  static MPO mpoId(L); // Basic MPO: all Identity
  static Operator idOp(reshape(identityMatrix(d),Indices(d,1,d,1)));
  static Operator* xOp; // individual operator sigmax
  static bool init(false);
  if(!init){
    complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    mwArray sigX(Indices(d,1,d,1),dataX);
    xOp=new Operator(sigX);
    for(int k=0;k<L;k++)
      mpoId.setOp(k,&idOp,0);
    init=true;
  }
  Contractor& contractor=Contractor::theContractor();
  polX.clear();
  for(int pos=0;pos<L;pos++){
    mpoId.setOp(pos,xOp,false);
    complex_t nk=contractor.contract(gs,mpoId,gs);
    mpoId.setOp(pos,&idOp,false);
    polX.push_back(real(nk));
  }
}


