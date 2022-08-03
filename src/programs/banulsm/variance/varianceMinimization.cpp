
#include <math.h>
#include <iomanip>

#include "mwArray.h"
#include "MPS.h"
#include "misc.h"
#include "Contractor.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"
#include "HeisenbergHamiltonian.h"

using namespace shrt;

/** Runs the minimizeVariance routine with the MPO for the Ising
    Hamiltonian \ref <IsingHamiltonian> at a given energy, and saves
    the results (MPS) in a file.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <E0> (double) energy value
    \param <penH> (double) penalty fort the energy term
    \param <outfname> (char*) name of the output file for the MPS
*/

const string mpsFileName(const string& mpsdir,const int N,const double E0,const int D);
int findInitFile(string& str,const string& initmpsdir,const int N,const double E0,const int D);
const string jobfilename(int L,int D,double E0,const string& strMod,const string jobsdir,int parNr);
double getSchmidtVals(vector<double>& schmidt,const MPS& mps,int pos,char gaugeDir);
void expandRDM(const MPS& mps,int posL,int posR,mwArray& A);
double computeRDMdistanceToId(const mwArray& rdm,int Lc);

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

  double E0=props.getDoubleProperty("E");
  double penH=props.getDoubleProperty("penH");
  double maxDevE=props.getDoubleProperty("toleranceE");
  if(maxDevE<0) maxDevE=max(1E-2,.1*E0);
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  int D=props.getIntProperty("D");
  int incrD=props.getIntProperty("incrD");
  int maxD=props.getIntProperty("maxD");
  string outfname=props.getProperty("output");
  string Sfname=props.getProperty("outputS");
  string mpsdir=props.getProperty("mpsdir");
  string jobsdir=props.getProperty("jobsdir");
  int Lcmax=props.getIntProperty("Lcmax");
  if(Lcmax<0||Lcmax>10) Lcmax=6;
  int Lb=props.getIntProperty("Lb");
  if(Lb<0) Lb=6;
  bool avrg=1;
  avrg=!(props.getIntProperty("average")==0); // whether to average for each Lc subsystem
  bool app=(props.getIntProperty("append")!=0);

  bool isXY=(!ising&&Jz_==0);
  int parNr=props.getIntProperty("parNr"); // Just for job identifier
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
  cout<<", outfile="<<outfname
      <<", E="<<E0
      <<", app="<<app
      <<", D="<<D<<endl;

  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    if(ising){
      *out<<"% N\t J\t g\t h\t D\t <H2>\t <H>";
    }
    else{
      if(!isXY)
	*out<<"% N\t Jx\t Jy\t Jz\t h\t D\t <H2>\t <H>";
      else
	*out<<"% N\t Jx\t Jy\t hx\t hy\t hz\t D\t <H2>\t <H>";
    }
    for(int kl=Lcmax;kl>=2;kl-=2){
      *out<<"\t"<<kl<<"\t dist(rho_"<<kl<<")";
    }
    *out<<"\t S(L/2)";
    *out<<endl;
    out->close();
    delete out;
  }

  ofstream* outS;
  if(!file_exists(Sfname.data())){
    outS=new ofstream(Sfname.data());
    if(!outS->is_open()){
      cout<<"Error: impossible to open file "<<Sfname<<
	" for output"<<endl;
      exit(1);
    }
    if(ising){
      *outS<<"% N\t J\t g\t h\t D\t lambda[0]^2\t lamdba[1]^2..."<<endl;
    }
    else{
      if(!isXY)
	*outS<<"% N\t Jx\t Jy\t Jz\t h\t D\t lambda[0]^2\t lamdba[1]^2..."<<endl;
      else
	*outS<<"% N\t Jx\t Jy\t hx\t hy\t hz\t D\t lambda[0]^2\t lamdba[1]^2..."<<endl;
    }
    *outS<<endl;
    outS->close();
    delete outS;
  }

  cout<<setprecision(10);

  string mpsfile=mpsFileName(mpsdir,L,E0,D);

  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  contractor.setConvTol(tol);
  cout<<"Initialized Contractor"<<endl;
  string tmpfile(mpsfile);tmpfile.append(".tmp");
  int d=2;
  MPS gs(L,D,d);
  bool readInit=false; // whether I success in reading an initial guess
  // FIRST: Check if the same file exisits (then refine)
  if(file_exists(mpsfile)){
    gs.importMPS(mpsfile.data());
    cout<<"Imported initial state from the same file "<<mpsfile<<endl;
    readInit=true;
  }
  else{ // check for tmpfile
    if(file_exists(tmpfile)){
      gs.importMPS(tmpfile.data());
      cout<<"Imported initial state from tmp file "<<tmpfile<<endl;
      readInit=true;
    }    
    else{ // check for a MPS with smaller D
      string initmpsfile;
      int Dinit=findInitFile(initmpsfile,mpsdir,L,E0,D);
      if(Dinit>0){
	gs.importMPS(initmpsfile.data());    readInit=true;
	cout<<"Imported initial state from file "<<initmpsfile<<endl;
      }
    }
  }
  if(!readInit){
    gs.setRandomState(); // the intial state, random
    cout<<"Initialized random state"<<endl;
  }
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  
  // First: put H and find the GS
  Hamiltonian* ham0;
  if(ising){
    ham0=new IsingHamiltonian(L,d,J_,g_,h_);
  }
  else{
    ham0=new HeisenbergHamiltonian(L,Jx_,Jy_,Jz_,hx_,hy_,h_);
  }
  //IsingHamiltonian hamH0(L,d,J_,g_,h_);
  const MPO& hamil=ham0->getHMPO();
  cout<<"Constructed the hamil MPO"<<endl;
  
  // place for other params
  complex_t E=contractor.contract(gs,hamil,gs);

  complex_t E2=contractor.contract2(hamil,gs);

  //hamil.exportForMatlab("isingHmpo.m");
  cout<<"Initial value of E, with initial state"<<endl;
  cout<<E<<endl;
  cout<<"and <H^2>="<<E2<<endl;

  // const MPO* ptrs[]={&hamil,&hamil};
  // MPO ham2(L);
  // MPO::join(2,ptrs,ham2);// MPO for H^2
  
  double varE=real(E2)-real(E)*real(E);
  cout<<"Initial variance: "<<varE<<endl;
  //  contractor.findGroundState(ham2,D,&expH2,gs);
  contractor.minimizeVariance(gs,D,hamil,E0,penH,&varE,tmpfile,500);
  gs.gaugeCond('R',1);gs.gaugeCond('L',1);

  E=contractor.contract(gs,hamil,gs);
  E2=contractor.contract2(hamil,gs);
  cout<<"After the optimization, cost function is "<<varE<<endl;
  cout<<"Variance="<<E2-E*E<<" <H^2>="<<E2<<", <H>="<<E<<endl;

  // Q: Am I close enough to target energy?; do I need to adjust penalty?
  bool repeat=0;
  double errorE=abs(E-E0*ONE_c);
  if(errorE>maxDevE){
    repeat=1;
    penH=10*penH;
    //penH=10*tol/errorE*errorE;
  }

  if(mpsfile.size()>0){
    gs.exportMPS(mpsfile.data());
    cout<<"Exported solution to file "<<mpsfile<<endl;
  }
  remove(tmpfile.data());
  
  // Compute and save same things as with Cheby

  double Shalf=contractor.getEntropy(gs);
  // Save also the Schmidt values
  vector<double> schmidt;
  double entroSch=getSchmidtVals(schmidt,gs,L/2,'L');
  // debugging
  //    cout<<"Entropy computed to be "<<Shalf<<" and from Schmidt coeffs "<<entroSch<<endl;
  outS=new ofstream(Sfname.data(),ios::app);
  *outS<<setprecision(15);
  *outS<<L<<"\t";
  if(ising)
    *outS<<J_<<"\t"<<g_<<"\t"<<h_<<"\t";
  else
    if(!isXY)
      *outS<<Jx_<<"\t"<<Jy_<<"\t"<<Jz_<<"\t"<<h_<<"\t";
    else
      *outS<<Jx_<<"\t"<<Jy_<<"\t"<<hx_<<"\t"<<hy_<<"\t"<<h_<<"\t";
  *outS<<D<<"\t";
  for(int vS=0;vS<D;vS++)
    *outS<<schmidt[vS]<<"\t";
  *outS<<endl;
  outS->close(); delete outS;
    
  // also RDM
  vector<double> trDistRDM;
  vector<int> sizeRDM;
  if(!avrg){
    int minPosR=L/2; // assuming L even and smallest Lc is 2
    for(int k=L-1;k>=minPosR+1;k--){
      gs.gaugeCond(k,'L',true);
    }
    int Lc=Lcmax;
    int posL=(L-Lc)/2;int posR=posL+Lc-1; // outside tracing them
    while(Lc>=2){
      mwArray oper;
      expandRDM(gs,posL,posR,oper);
      //cout<<"Created RDM for Lc="<<Lc<<endl;
      double singvals=0.;
      double idFc=1./pow(2,Lc);
      mwArray U;vector<complex_t> Dv;
      wrapper::eig(oper,Dv,U); // only eigenvalues
      //  cout<<"Diagonalized the rho(Lc)"<<endl;
      // now sum the absolute values
      for(int k=0;k<Dv.size();k++){
	double aux=real(Dv[k]);
	singvals+=abs(aux-idFc);
      }
      if(Dv.size()<pow(2,Lc)){
	cout<<"WARNING! Seems the dimension of the operator for D="<<D
	    <<", Lc="<<Lc<<"is not full: dim(oper)="<<Dv.size()
	    <<", 2^Lc="<<pow(2,Lc)<<endl;
	singvals+=idFc*(pow(2,Lc)-Dv.size());
      }
      sizeRDM.push_back(Lc);
      trDistRDM.push_back(singvals);
      posL++;posR--;Lc=(posR-posL+1);
    }
  }
  else{ // averaging over different subsystems of the same size
    for(int lc=1;lc<=Lcmax;lc++){sizeRDM.push_back(lc);} // position lc-1
    trDistRDM=vector<double>(sizeRDM.size(),0.);
    vector<int> avrgCnt(sizeRDM.size(),0.);
    int posL=Lb;
    while(posL<L-Lb-1){ // at least there is a 2-site RDM to compute
      int lc=Lcmax;
      int posR=posL+lc-1;
      while(posR>=L){ // cannot compute this one, need to reduce the size
	lc--;
	posR=posL+lc-1;
      }
      int dimL=pow(d,lc);
      int dimR=dimL;
      mwArray oper=contractor.getRDM(gs,posL,posR);
      while(lc>=1){ // Compute the distance and trace another site
	if(posR<L-Lb){
	  // cout<<"Computing distance for rmd("<<lc<<") sites, pos:"<<posL<<"-"<<posR<<endl;
	  double dist_lc=computeRDMdistanceToId(oper,lc);
	  trDistRDM[lc-1]+=dist_lc;
	  avrgCnt[lc-1]++; // to take the average later
	}
	if(lc>=1){
	  // trace out two sites
	  oper.reshape(Indices(dimL/d,d,dimR/d,d));
	  oper.permute(Indices(1,3,2,4));
	  dimL=dimL/d;dimR=dimR/d;
	  oper.reshape(Indices(dimL*dimR,d*d));
	  oper.multiplyRight(reshape(identityMatrix(d),Indices(d*d,1)));
	  oper.reshape(Indices(dimL,dimR));
	  posR--;
	}
	lc--;
      }    
      posL++;
    }
    for(int lc=1;lc<=Lcmax;lc++) trDistRDM[lc-1]=trDistRDM[lc-1]/avrgCnt[lc-1];
  }

  out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(15);
  // OJO!!!!! This is a mistake, doing it wrong for Heisenberg!!!!
  if(ising)
    *out<<L<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t";
  else
    if(!isXY)
      *out<<L<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t";
  //  *out<<L<<"\t"<<Jx_<<"\t"<<Jy_<<"\t"<<Jz_<<"\t"<<h_<<"\t"<<D<<"\t"; // should be this
    else
      *out<<L<<"\t"<<Jx_<<"\t"<<Jy_<<"\t"<<hx_<<"\t"<<hy_<<"\t"<<h_<<"\t"<<D<<"\t";
  *out<<real(E2)<<"\t"<<real(E)<<"\t";
  for(int ks=0;ks<sizeRDM.size();ks++){
    *out<<sizeRDM[ks]<<"\t"<<trDistRDM[ks]<<"\t";
  }
  *out<<Shalf<<endl;
  *out<<endl;
  out->close(); delete out;


  // ************* End of saving data

  // Prepare the script with the new job
  if(jobsdir.length()>0&&D<maxD){
    int newD=D+incrD;
    string newOutput(outfname);
    string newSfname(Sfname);
    //    stringstream newSuff;newSuff<<"_D"<<newD; 
    //newOutput.replace(newOutput.find("_D"),string::npos,newSuff.str());
    stringstream commandstr;
    commandstr<<argv[0]; // exec
    commandstr<<" "<<argv[1]; // config file
    commandstr<<" -L="<<L;
    commandstr<<" -J="<<J_<<" -h="<<h_<<" -g="<<g_;
    commandstr<<" -ising="<<ising;
    commandstr<<" -Jx="<<Jx_<<" -Jy="<<Jy_<<" -Jz="<<Jz_<<" -hx="<<hx_<<" -hy="<<hy_;
    commandstr<<" -D="<<newD;
    if(newD==100){maxD=1000;incrD=100;}
    commandstr<<" -maxD="<<maxD;
    commandstr<<" -incrD="<<incrD;
    commandstr<<" -E="<<E0;
    commandstr<<" -penH="<<penH;
    commandstr<<" -parNr="<<parNr;
    commandstr<<" -output="<<newOutput;
    commandstr<<" -outputS="<<newSfname;
    commandstr<<" -mpsdir="<<mpsdir; //<<" -initmpsfile="<<mpsfile; // use this as init
    commandstr<<" -jobsdir="<<jobsdir;
    commandstr<<" -append="<<1;
    commandstr<<" -tol="<<tol;
    commandstr<<" -Lcmax="<<Lcmax<<" -Lb="<<Lb<<" -average="<<avrg;

    string strMod=ising?"Ising":(isXY?"XY":"Heis");
    string jobfile=jobfilename(L,newD,E0,strMod,jobsdir,parNr);

    ofstream* ojob=new ofstream(jobfile.data());
    if(!ojob->is_open()){
      cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
      cout<<commandstr.str()<<endl;
      cout<<endl;
    }
    else{
      int nPar=4;
      if(newD>=80) nPar=8;
      if(newD>=100) nPar=10;
      if(newD>=200) nPar=12;
      if(newD>=500){ nPar=20; if(L>=100) nPar=40;}
      *ojob<<"msub_slurm -N minVar."<<strMod<<".N"<<L<<".D"<<newD<<".E0"<<E0<<".par"<<parNr;
      *ojob<<" -P "<<nPar;
      *ojob<<" -- ";
      *ojob<<commandstr.str()<<endl;
      ojob->close();
    }
    delete ojob;

  }

}


const string jobfilename(int L,int D,double E0,const string& strMod,const string jobsdir,int parNr){
  stringstream str;
  str<<jobsdir<<"/job_"<<strMod;
  str<<"_L"<<L<<"_E0"<<E0<<"_D"<<D<<"_par"<<parNr;
  str<<".dat";
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return str.str();
}

int findInitFile(string& str,const string& initmpsdir,const int N,const double E0,const int D){
  cout<<"Searching for init file for D "<<D<<" in dir "<<initmpsdir<<endl;
  for(int myD=D-1;myD>0;myD--){
    // TRy with one smaller D
    str=mpsFileName(initmpsdir,N,E0,myD);
    if(file_exists(str)){
      cout<<"Found file "<<str<<endl;
      return myD;
    }
  }
  // if nothing found
  str="";
  return 0; 
}

const string mpsFileName(const string& mpsdir,const int N,const double E0,const int D){
  stringstream s;
  s<<mpsdir<<"/MPS_N"<<N<<"_E0"<<E0<<"_D"<<D;
  return s.str();
}

double getSchmidtVals(vector<double>& schmidt,const MPS& mps,int pos,char gaugeDir){
  Contractor& contractor=Contractor::theContractor();
  int len=mps.getLength();
  mwArray rdm=contractor.contract(mps,mps,pos,gaugeDir);
  // cout<<"getSchmidtVals of mps:"<<mps<<" at cut pos="<<pos<<"\nComputed the rdm: "<<rdm.getDimensions()<<endl;
  // Now diagonalize
  rdm.permute(Indices(2,1)); // here it does not matter
  vector<complex_t> Dval;mwArray U; //placeholder:ignored
  int D=rdm.getDimension(0); //mps.getA(pos).getDl();
  wrapper::eig(rdm,Dval,U);
  schmidt.clear();
  for(int k=0;k<D;k++) schmidt.push_back(real(Dval[k]));
  double entropy=0.;
  for(int k=0;k<D;k++){
    double tmp=real(Dval[k]);
    if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
  }
  //cout<<"getSchmidtVals computed vals="<<Dval<<"\n and Entropy="<<entropy<<endl;
  return entropy;
  
}

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


 double computeRDMdistanceToId(const mwArray& rdm,int Lc){
   double singvals=0.;
   double idFc=1./pow(2,Lc);
   mwArray U;vector<complex_t> Dv;
   wrapper::eig(rdm,Dv,U); // only eigenvalues
   // now sum the absolute values
   for(int k=0;k<Dv.size();k++){
     double aux=real(Dv[k]);
     singvals+=abs(aux-idFc);
   }
   if(Dv.size()<pow(2,Lc)){
     cout<<"WARNING! Seems the dimension of the operator for Lc="<<Lc
	 <<"is not full: dim(oper)="<<Dv.size()
	 <<", 2^Lc="<<pow(2,Lc)<<endl;
     singvals+=idFc*(pow(2,Lc)-Dv.size());
   }
   return singvals;
 }
