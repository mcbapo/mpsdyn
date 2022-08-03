
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#include "SpinMPO.h"
#include "StochasticHamiltonian.h"
#include "FAHamiltonian.h"
#include "Properties.h"
#include "misc.h"

#include "quicksort.h"

#define FREQSAVE 18000 // frequency of tmp saving: default 5 hrs 


/** 
Generates a file name for the MPS corresponding to a certain excitation.
*/
const string mpsfilename(int L,int D,double s,double c,int level,const string mpsdir);
const string jobfilename(int L,int D,double s,double c,const string jobsdir);

/** Compute individual occupations of sites */
void computeOccupations(const MPS& gs,int L,int d,vector<double>& occ);

/** Same for polarizations */
void computePolarizations(const MPS& gs,int L,int d,vector<double>& occ);


/** eastEntropy takes the states computed by basicSpec (FASpec)
    and computes their entropy.

    Receives as argument a Properties file, such as
    config/stochspec.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.

    Receives, among others, arguments:
    \param <L> (int) length of the chain
    \param <s> (double) parameter \f$s\f$ of the Hamiltonian
    \param <c> (double) parameter \f$c\f$ of the Hamiltonian
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
  double c=props.getDoubleProperty("c");
  if(c<0||c>1){
    cout<<"Error: value of c not allowed (c in[0,1])"<<endl;
    exit(1);
  }

  int D0=props.getIntProperty("D");
  int maxD=props.getIntProperty("maxD");
  int incrD=props.getIntProperty("incrD");
  if(incrD<0) incrD=20;
  int nlev=props.getIntProperty("nlevel");
  if(nlev<0) nlev=0;
  string outfname=props.getProperty("outputfile");
  string outfnameSch=props.getProperty("outputfileSchmidt");
  bool saveSchmidt=outfnameSch.length()>0; // whether to save also Schmidt values in the middle of the chain
  string mpsdir=props.getProperty("mpsdir");
  //  string jobsdir=props.getProperty("jobsdir");
  int append=props.getIntProperty("append");
  bool app=(append!=0);


  bool sectorOne=(props.getIntProperty("firstN")!=0);
  cout<<"Occupation of first site fixed to "<<(sectorOne?1:0)<<endl;

  cout<<"Appending to output?"<<app<<endl;  

  bool origApp=app;
  ofstream* out;ofstream* outS;
  if(!app||!file_exists(outfname)){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% Entropy of East model L="<<L<<endl;
    *out<<"% c="<<c<<endl;
    *out<<"% s="<<s<<endl;
    *out<<"% nlev\t D\t <H>\t<H^2>\t pos(L/2)\t S(L/2)"<<endl;
    out->close();delete out;
  }
  
  if(saveSchmidt){
    if(!app||!file_exists(outfnameSch)){
      outS=new ofstream(outfnameSch.data());
      if(!outS->is_open()){
	cout<<"Error: impossible to open file "<<outfnameSch<<
	  " for output"<<endl;
	exit(1);
      }
      *outS<<"% Schmidt values of East model L="<<L<<endl;
      *outS<<"% c="<<c<<endl;
      *outS<<"% s="<<s<<endl;
      *outS<<"% nlev\t D\t values at L/2"<<endl;
      outS->close();delete outS;
    }
  }

  
  //out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
  //out->close();delete out;
    
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;

  int d=2;

  double offset=0.;
  StochasticHamiltonian hSto(L-1,s,offset,sectorOne,c);
  const MPO& hamil=hSto.getHMPO();

    // Now loop over Ds
  int D=D0;
  while(D<=maxD){

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
	
	double entropy=contractor.getEntropy(nextLev,L/2-1);
	
	out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
	*out<<lval<<"\t"<<D<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t";
	*out<<L/2-1<<"\t"<<entropy<<"\t";      
	*out<<endl;
	out->close();delete out;

	if(saveSchmidt){
	  vector<complex_t> schmVals;
	  contractor.getSchmidtValues(nextLev,schmVals,L/2-1);
	  //while(schmVals.size()<maxD) schmVals.push_back(ZERO_c);
	  outS=new ofstream(outfnameSch.data(),ios::app);*outS<<setprecision(15);
	  *outS<<lval<<"\t"<<D<<"\t";
	  int nrVals=schmVals.size();
	  for(int id=0;id<maxD;id++){
	    if(id<nrVals)
	      *outS<<abs(schmVals[nrVals-id-1])<<"\t";
	    else
	      *outS<<0.<<"\t";
	  }
	  *outS<<endl;
	  outS->close();delete outS;
	}
	
	lval++;
      }
      else{
	oldMPS=false; // break the loop 
      }
    }
    D+=incrD;
  }
  
}


const string mpsfilename(int L,int D,double s,double c,int level,const string mpsdir){
  stringstream str;
  str<<mpsdir<<"/MPS_L"<<L<<"_s"<<s<<"_c"<<c<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return str.str();
}

const string jobfilename(int L,int D,double s,double c,const string jobsdir){
  stringstream str;
  str<<jobsdir<<"/job_entropy_L"<<L<<"_s"<<s<<"_c"<<c<<"_D"<<D;
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


