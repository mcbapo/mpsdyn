
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#include "BoseHubbardHamiltonian.h"
#include "BosonMPO.h"
#include "Properties.h"
#include "misc.h"

/** 
Generates a file name for the MPS corresponding to a certain excitation.
*/
const string mpsfilename(int L,int N,int D,double t,double U,double mu,
			 int level,const string mpsdir);

/** Compute individual occupations of sites */
void computeOccupations(const MPS& gs,int L,int d,vector<double>& occ);

/** Compute all b(pos0)b(pos0+l) correlations (which will be returned in bbcorr,
    and all bDagger(pos0)b(pos0+l) correlations (returned in bdagbcorr) */
void computeCorrelations(const MPS& gs,int L,int d,int pos0,vector<double>& bbcorr,
			 vector<double>& bdagbcorr);

/** Compute all b(pos0)b(pos0+l) correlations (which will be returned
    in bbcorr, and all bDagger(pos0)b(pos0+l) correlations (returned
    in bdagbcorr), but also individual values b(pos) and
    bdagger(pos). */
void computeCorrelationsAndSingleOps(const MPS& gs,int L,int d,int pos0,
				     vector<complex_t>& b,vector<complex_t>& bdag,
				     vector<complex_t>& bbcorr,vector<complex_t>& bdagbcorr);

/** Check whether there is a file with smaller bond dimension that can be used as initial guess. */
int findInitFile(string& str,int L,int N,int D,double t,double U,double mu,
		  int level,const string initmpsdir);

/** bosehub tries to find the ground state and excitations of the Bose-Hubbard Hamiltonian
    \ref <BoseHubbardHamiltonian>

    Receives as argument a Properties file, such as
    config/bhspec.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <N> (int) cutoff for the occupation number
    \param <t> (double) parameter \f$t\f$ of the Hamiltonian
    \param <U> (double) parameter \f$U\f$ of the Hamiltonian
    \param <mu> (double) parameter \f$\mu\f$ of the Hamiltonian
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
  int N=props.getIntProperty("N");
  double t=props.getDoubleProperty("t");
  double U=props.getDoubleProperty("U");
  double mu=props.getDoubleProperty("mu"); // chemical potential
  int D=props.getIntProperty("D");
  int nlev=props.getIntProperty("nlevel");
  string outfname=props.getProperty("outputfile");
  string outfnameC=props.getProperty("correlationsfile");
  int posC=props.getIntProperty("pos0");
  if(posC>=0&&outfnameC.size()==0){
    cout<<"ERROR: cannot compute correlations because I do not have a file where to save them"<<endl;
    exit(1);
  }
  cout<<"Read pos0="<<posC<<endl;
  if(posC<0) posC=0;
  string mpsdir=props.getProperty("mpsdir");
  string initmpsdir=props.getProperty("initmpsdir");
  bool app=(props.getIntProperty("append")!=0);
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  
  bool rewrite=false;
  ofstream* out;
  if(!app||!file_exists(outfname)){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% Spectrum of BH model L="<<L<<", N="<<N<<", D="<<D<<endl;
    *out<<"% t="<<t<<", U="<<U<<", mu="<<mu<<endl;
    *out<<"% nlev\t <H>\t<H^2>\t <N>\t <N^2>\t <n(0)>\t <n(1)>\t ..."<<endl;
    out->close();delete out;
    rewrite=true;
  }
  //out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
  //out->close();delete out;
    
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);
  int d=N+1;

  double penN=0.; // no penalty term
  int Ntot=0.; // no total nr of bosons imposed
  BoseHubbardHamiltonian hBH(L,N,t,U,Ntot,penN,mu);
  const MPO& hamil=hBH.getHMPO();
  //hamil.exportForMatlab("testBH.m");exit(1);
  
 // Observable(s) to compute
  MPO Nmpo(L);
  BosonMPO::getNumberMPO(L,N,Nmpo);
  
  vector<MPS*> levels;

  // First, check how many MPS I already have in mpsdir
  bool oldMPS=true;
  int lval=0;
  MPS* nextLev;
  while(oldMPS&&lval<=nlev-1){
    const string filename=mpsfilename(L,N,D,t,U,mu,lval,mpsdir);
    // Try to open file for MPS level lval
    if(file_exists(filename)){
      nextLev=new MPS(L,D,d);
      nextLev->importMPS(filename.data());
      nextLev->gaugeCond('R',1);
      nextLev->gaugeCond('L',1);
      levels.push_back(nextLev);
      cout<<"Recovered MPS file for level "<<lval<<" from file "<<filename<<endl;
      if(rewrite){
	complex_t Ek=contractor.contract(*nextLev,hamil,*nextLev);
	complex_t Ek2=contractor.contract2(hamil,*nextLev);
	complex_t Nk=contractor.contract(*nextLev,Nmpo,*nextLev);
	complex_t Nk2=contractor.contract2(Nmpo,*nextLev);
	vector<double> occ;
	computeOccupations(*nextLev,L,d,occ);
	out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
	*out<<lval<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t"<<real(Nk)<<"\t"<<real(Nk2)<<"\t";
	for(int pos=0;pos<L;pos++)
	  *out<<occ[pos]<<"\t";
	*out<<endl;
	out->close();delete out;
	if(outfnameC.size()>0){
	  // vector<double> bbcorr;
	  // vector<double> bdagbcorr;
	  //computeCorrelations(*nextLev,L,d,posC,bbcorr,bdagbcorr);
	  vector<complex_t> bvec;
	  vector<complex_t> bdagvec;
	  vector<complex_t> bbcorr;
	  vector<complex_t> bdagbcorr;
	  computeCorrelationsAndSingleOps(*nextLev,L,d,posC,bvec,bdagvec,bbcorr,bdagbcorr);
	  out=new ofstream(outfnameC.data(),ios::app);*out<<setprecision(15);
	  *out<<lval<<"\t"<<0<<"\t";
	  for(int k=0;k<bvec.size();k++){
	    *out<<bvec[k]<<"\t";
	  }
	  *out<<endl;
	  *out<<lval<<"\t"<<1<<"\t";
	  for(int k=0;k<bdagvec.size();k++){
	    *out<<bvec[k]<<"\t";
	  }
	  *out<<endl;
	  *out<<lval<<"\t"<<2<<"\t";
	  for(int k=0;k<bbcorr.size();k++){
	    *out<<bbcorr[k]<<"\t";
	  }
	  *out<<endl;
	  *out<<lval<<"\t"<<3<<"\t";
	  for(int k=0;k<bdagbcorr.size();k++){
	    *out<<bdagbcorr[k]<<"\t";
	  }
	  *out<<endl;
	  out->close();delete out;
	}
      }      
      lval++;
    }
    else{
      oldMPS=false; // break the loop 
    }
  }
  cout<<"Recovered "<<lval<<" already computed MPS levels from disk. "
      <<"Computing next level now."<<endl;

  while(lval<=nlev-1){
    MPS* exc=new MPS(L,D,d);
    // First check if a valid initial guess exists in the initmpsdir directory
    bool found=0;
    if(initmpsdir.size()>0){
      string filename;
      int Dinit=findInitFile(filename,L,N,D,t,U,mu,lval,initmpsdir);
      if(Dinit>0){
	cout<<"Using initial guess for D="<<D<<" from result for "<<Dinit<<" in file "
	    <<filename<<endl;
	exc->importMPS(filename.data());
	found=true;
      }
    }
    if(!found){
      cout<<"Initializing the search for level "<<lval<<" with a random MPS"<<endl;
      exc->setRandomState(); // the initial state, random -> not to the GS!!
    }
    exc->gaugeCond('R',1);
    exc->gaugeCond('L',1);
    double lambda=0.;

    if(lval==0){
      contractor.findGroundState(hamil,D,&lambda,*exc);
      cout<<"Ground state found with eigenvalue "<<lambda<<endl;
    }
    else{
      contractor.findNextExcitedState(hamil,D,levels,&lambda,*exc,0.,1);
      cout<<"Next("<<lval<<") excited state found with E="<<lambda<<endl;
    }
    exc->gaugeCond('L',1);
    // Save to file
    const string filename=mpsfilename(L,N,D,t,U,mu,lval,mpsdir);
    exc->exportMPS(filename.data());

    levels.push_back(exc);
    complex_t Ek=contractor.contract(*exc,hamil,*exc);
    complex_t Ek2=contractor.contract2(hamil,*exc);
    complex_t Nk=contractor.contract(*exc,Nmpo,*exc);
    complex_t Nk2=contractor.contract2(Nmpo,*exc);
    vector<double> occ;
    computeOccupations(*exc,L,d,occ);
    out=new ofstream(outfname.data(),ios::app);*out<<setprecision(15);
    *out<<lval<<"\t"<<real(Ek)<<"\t"<<real(Ek2)<<"\t"<<real(Nk)<<"\t"<<real(Nk2)<<"\t";
    for(int pos=0;pos<L;pos++)
      *out<<occ[pos]<<"\t";
    *out<<endl;
    out->close();delete out;
    cout<<"Computed and saved occupations"<<endl;
    if(outfnameC.size()>0){
      // vector<double> bbcorr;
      // vector<double> bdagbcorr;
      // computeCorrelations(*exc,L,d,posC,bbcorr,bdagbcorr);
      vector<complex_t> bvec;
      vector<complex_t> bdagvec;
      vector<complex_t> bbcorr;
      vector<complex_t> bdagbcorr;
      computeCorrelationsAndSingleOps(*exc,L,d,posC,bvec,bdagvec,bbcorr,bdagbcorr);
      out=new ofstream(outfnameC.data(),ios::app);*out<<setprecision(15);
      *out<<lval<<"\t"<<0<<"\t";
      for(int k=0;k<bvec.size();k++){
	*out<<bvec[k]<<"\t";
      }
      *out<<endl;
      *out<<lval<<"\t"<<1<<"\t";
      for(int k=0;k<bdagvec.size();k++){
	*out<<bdagvec[k]<<"\t";
      }
      *out<<endl;
      *out<<lval<<"\t"<<2<<"\t";
      for(int k=0;k<bbcorr.size();k++){
	*out<<bbcorr[k]<<"\t";
      }
      *out<<endl;
      *out<<lval<<"\t"<<3<<"\t";
      for(int k=0;k<bdagbcorr.size();k++){
	*out<<bdagbcorr[k]<<"\t";
      }
      *out<<endl;
      out->close();delete out;
    }
    lval++;
  }

}


int findInitFile(string& str,int L,int N,int D,double t,double U,double mu,
			 int level,const string initmpsdir){
  cout<<"Searching for init file for level "<<level<<" in dir "<<initmpsdir<<endl;
  for(int myD=D-1;myD>0;myD--){
    // TRy with one smaller D
    str=mpsfilename(L,N,myD,t,U,mu,level,initmpsdir);
    if(file_exists(str)){
      cout<<"Found file "<<str<<endl;
      return myD;
    }
    // else{
    //   cout<<"No file "<<str<<" found"<<endl;
    // }
  }
  // if nothing found
  str="";
  return 0; 
}

const string mpsfilename(int L,int N,int D,double t,double U,double mu,
			 int level,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS_L"<<L<<"_t"<<t<<"_U"<<U<<"_mu"<<mu<<"_nmax"<<N<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
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


void computeCorrelations(const MPS& gs,int L,int d,int pos0,vector<double>& bbcorr,
			 vector<double>& bdagbcorr){
  static MPO mpoId(L); // Basic MPO: all Identity
  static Operator idOp(reshape(identityMatrix(d),Indices(d,1,d,1)));
  static Operator* bOp; // individual annihilation operator 
  static Operator* bDagOp; // individual creation operator 
  static bool init(false);
  if(!init){
    mwArray bCreat(Indices(d,d)); // individual operator bDagger
    int N=d-1;
    for(int k=0;k<N;k++){
      bCreat.setElement(sqrt(k+1),0.,Indices(k+1,k));
    }
    mwArray bDestr(Hconjugate(bCreat));
    init=true;
    bOp=new Operator(bDestr);
    bDagOp=new Operator(bCreat);
    for(int k=0;k<L;k++)
      mpoId.setOp(k,&idOp,0);
    init=true;
  }
  Contractor& contractor=Contractor::theContractor();
  bbcorr.clear();bdagbcorr.clear();
  for(int pos=pos0+1;pos<L;pos++){
    mpoId.setOp(pos,bOp,false);
    mpoId.setOp(pos0,bOp,false);
    complex_t corrBB=contractor.contract(gs,mpoId,gs);
    mpoId.setOp(pos0,bDagOp,false);
    complex_t corrBdagB=contractor.contract(gs,mpoId,gs);
    mpoId.setOp(pos0,&idOp,false); // restore
    mpoId.setOp(pos,&idOp,false); // restore
    bbcorr.push_back(real(corrBB));
    bdagbcorr.push_back(real(corrBdagB));
  }

}

void computeCorrelationsAndSingleOps(const MPS& gs,int L,int d,int pos0,
				     vector<complex_t>& b,vector<complex_t>& bdag,
				     vector<complex_t>& bbcorr,vector<complex_t>& bdagbcorr){
  static MPO mpoId(L); // Basic MPO: all Identity
  static Operator idOp(reshape(identityMatrix(d),Indices(d,1,d,1)));
  static Operator* bOp; // individual annihilation operator 
  static Operator* bDagOp; // individual creation operator
  static Operator* numOp; // individual number operator (bdagb)
  static bool init(false);
  if(!init){
    mwArray bCreat(Indices(d,d)); // individual operator bDagger
    int N=d-1;
    for(int k=0;k<N;k++){
      bCreat.setElement(sqrt(k+1)*ONE_c,Indices(k+1,k));
    }
    mwArray bDestr(bCreat);bDestr.Hconjugate();
    init=true;
    bOp=new Operator(bDestr);
    bDagOp=new Operator(bCreat);
    numOp=new Operator(bCreat*bDestr);
    for(int k=0;k<L;k++)
      mpoId.setOp(k,&idOp,0); // all Id
    init=true;
  }
  Contractor& contractor=Contractor::theContractor();
  b.clear();bdag.clear();
  bbcorr.clear();bdagbcorr.clear();
  for(int pos=0;pos<L;pos++){
    mpoId.setOp(pos,bOp,false); // b in pos
    complex_t valB=contractor.contract(gs,mpoId,gs);
    complex_t corrBB,corrBdagB;
    if(pos0!=pos){
      mpoId.setOp(pos0,bOp,false); // b in pos0
      corrBB=contractor.contract(gs,mpoId,gs);
      mpoId.setOp(pos0,bDagOp,false); // bDag in pos0
      corrBdagB=contractor.contract(gs,mpoId,gs);
      mpoId.setOp(pos0,&idOp,false);   // restore Id in pos0
    }
    else{
      corrBB=ZERO_c; // b^2=0
      mpoId.setOp(pos0,numOp,false);   // bdag b on same pos=pos0
      corrBdagB=contractor.contract(gs,mpoId,gs);
      mpoId.setOp(pos0,&idOp,false);   // restore Id in pos0
    }
    mpoId.setOp(pos,bDagOp,false); // bDag in pos
    complex_t valBdag=contractor.contract(gs,mpoId,gs);
    mpoId.setOp(pos,&idOp,false); // restore Id in pos
    bbcorr.push_back(corrBB);
    bdagbcorr.push_back(corrBdagB);
    b.push_back(valB);
    bdag.push_back(valBdag);
  }

}
