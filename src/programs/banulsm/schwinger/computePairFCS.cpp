#include <math.h>
#include <iomanip>
//#include <cmath>

#include "misc.h"
#include "MPS.h"
#include "MPO.h"
#include "Contractor.h"
#include "Properties.h"

#include "SchwingerHamiltonianSz.h"
//#include "SpinMPO.h"


using namespace std;
using namespace shrt;

#define MAXLEN 120

/** If set to 0, I take fermion to be a 1 on even, 0 on odd sites; if set to 1, reverse */
int FERMEVEN=1;

/** 
Generates a file name for the MPS corresponding to a certain excitation. Copied from single_state.cpp, to match the name
*/
const string mpsfilename(int L,int D,double x,double mg,int level,const string mpsdir);

void prepareMPOfermions(MPO& mpo,mwArray& midOp,int L);
void prepareMPOantifermions(MPO& mpo,mwArray& midOp,int L);

/** Read a MPS file (generated by single_state, potentially reordered
    later) and compute the probability distribution of the number of
    fermion-antifermion pairs. */

int d=2;

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
  double mg=props.getDoubleProperty("mg");
  double alpha=props.getDoubleProperty("alpha");
  double x=props.getDoubleProperty("x");
  int lev=props.getIntProperty("level");
  int D0=props.getIntProperty("D");
  double zpen=props.getDoubleProperty("penalty");
  double offset=props.getDoubleProperty("offset");
  string mpsdir=props.getProperty("mpsdir"); // default current
  if(mpsdir.empty()) 
    mpsdir=".";
  // Output
  string directory=props.getProperty("outputdir");
  const string outfile=directory+"/"+props.getProperty("outfile");
  int app=props.getIntProperty("append");
  if(app<0) app=1; // default will be to append


  int d=2;
  // Read the MPS
  const string filename=mpsfilename(L,D0,x,mg,lev,mpsdir);
  MPS state(L,D0,d);
  if(!file_exists(filename)){
    cout<<"ERROR: No file found for the desired MPS ("<<filename<<endl;
    exit(1);
  }
  state.importMPS(filename.data());
  cout<<"Read MPS from file"<<endl;

  // Prepare the MPO(s) to compute the FCS
  MPO allProjFerm(L);mwArray connFerm;
  MPO allProjAntiFerm(L);mwArray connAntiFerm; // the connection is actually the same
  prepareMPOfermions(allProjFerm,connFerm,L);
  prepareMPOantifermions(allProjAntiFerm,connAntiFerm,L);

  // and now I use a trick to get all FCS
  Contractor& contractor=Contractor::theContractor();

  // contract the sandwich MPS-MPO-MPS until mid point on both sides
  mwArray contrFermL=contractor.contract(state,allProjFerm,state,L/2,'L'); // D x xiL x D
  mwArray contrFermR=contractor.contract(state,allProjFerm,state,L/2-1,'R'); // D x xiR x D

  Indices midDim=connFerm.getDimensions(); // xiL,L+1,xiR
  int D=contrFermL.getDimension(0);
  contrFermL.permute(Indices(2,1,3)); // xiL x D x D
  contrFermL.reshape(Indices(midDim[0],D*D));
  contrFermR.permute(Indices(1,3,2)); // D x D x xiR
  contrFermR.reshape(Indices(D*D,midDim[2]));
  contrFermL.multiplyRight(contrFermR); // xiL x xiR
  contrFermL.reshape(Indices(midDim[0]*midDim[2],1));
  // contract with the midOp and get a vector where wach component is the prob of a nr
  connFerm.permute(Indices(2,1,3));
  connFerm.reshape(Indices(midDim[1],midDim[0]*midDim[2]));
  connFerm.multiplyRight(contrFermL);

  // Same for the antifermions part, as double check
  // contract the sandwich MPS-MPO-MPS until mid point on both sides
  mwArray contrAntiFermL=contractor.contract(state,allProjAntiFerm,state,L/2,'L'); // D x xiL x D
  mwArray contrAntiFermR=contractor.contract(state,allProjAntiFerm,state,L/2-1,'R'); // D x xiR x D

  Indices midAntiDim=connAntiFerm.getDimensions(); // xiL,L/2+1,xiR
  D=contrAntiFermL.getDimension(0);
  contrAntiFermL.permute(Indices(2,1,3)); // xiL x D x D
  contrAntiFermL.reshape(Indices(midAntiDim[0],D*D));
  contrAntiFermR.permute(Indices(1,3,2)); // D x D x xiR
  contrAntiFermR.reshape(Indices(D*D,midAntiDim[2]));
  contrAntiFermL.multiplyRight(contrAntiFermR); // xiL x xiR
  contrAntiFermL.reshape(Indices(midAntiDim[0]*midAntiDim[2],1));
  // contract with the midOp and get a vector where wach component is the prob of a nr
  connAntiFerm.permute(Indices(2,1,3));
  connAntiFerm.reshape(Indices(midAntiDim[1],midAntiDim[0]*midAntiDim[2]));
  connAntiFerm.multiplyRight(contrAntiFermL);

  ofstream* out;

  if(!file_exists(outfile)||!app){
    out=new ofstream(outfile.data());
    if(!out->is_open()){
      cout<<"Error: impossible to (re)create file "<<outfile<<
	" for output "<<endl;
      exit(1);
    }
    // Write header
    *out<<"%lev \t L \t mg\t x\t D\t type(0 ferm, 1 antif) \t list of L/2+1 pairs of nr-probability"<<endl;
  }
  else{
    out=new ofstream(outfile.data(),ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfile<<
	" for output (append="<<app<<")"<<endl;
      exit(1);
    }
  }
  *out<<setprecision(15);
  *out<<lev<<"\t"<<L<<"\t"<<mg<<"\t"<<x<<"\t"<<D0<<"\t"<<0<<"\t";
  for(int s=0;s<connFerm.getDimension(0);s++){
    *out<<s<<"\t"<<real(connFerm.getElement(Indices(s,0)))<<"\t";
  }
  *out<<endl;
  *out<<lev<<"\t"<<L<<"\t"<<mg<<"\t"<<x<<"\t"<<D0<<"\t"<<1<<"\t";
  for(int s=0;s<connAntiFerm.getDimension(0);s++){
    *out<<s<<"\t"<<real(connAntiFerm.getElement(Indices(s,0)))<<"\t";
  }
  *out<<endl;
      
  
}

/** I use a trick to prepare the MPO which computes all projectors at
    once: I prepare the MPO for the projector onto the largest
    possible number (so no cut) and provide a "connecting" operator to
    be inserted in the middle link. This operator is a rank-3 tensor,
    ordered Dlx dN x Dr, where Dl and Dr match the bond dimension of
    the MPO, and dN=L/2+1 corresponds to each allowed value of the
    number of pairs */
void prepareMPOfermions(MPO& mpo,mwArray& midOp,int L){
  // L has to be even
  if(L%2!=0){
    cout<<"ERROR: The number of sites must be even"<<endl;
    exit(1);
  }
  mpo.initLength(L);
  // The even sites are simply identity in physical and in virtual
  // dimensions
  // TODO: check who is even and odd and which is the "occupied" state
  // in each, as I had different conventions. I assume now that 1 is
  // the occupied for even (and start at 0, even)
  int D=1;
  for(int k=0;k<L/2;k+=2){
    // even site
    mwArray evenOp(Indices(d,D,d,D+1));evenOp.fillWithZero(); // redundant, but in case
    // I count the FERMEVEN
    int oneFerm=FERMEVEN;
    int noFerm=1-FERMEVEN;
    for(int p=0;p<D;p++){
      evenOp.setElement(ONE_c,Indices(noFerm,p,noFerm,p));
      evenOp.setElement(ONE_c,Indices(oneFerm,p,oneFerm,p+1));
    }
    mpo.setOp(k,new Operator(evenOp),true);
    D++;
    // odd site
    mwArray oddOp=reshape(identityMatrix(D),Indices(D*D,1))*reshape(identityMatrix(d),Indices(1,d*d));
    oddOp.reshape(Indices(D,D,d,d));
    oddOp.permute(Indices(3,1,4,2));
    mpo.setOp(k+1,new Operator(oddOp),true);
  }
  int Dl=D; // the one in the middle
  // now the right part (actually, symmetric, could use reflections)
  // This is so that the bond dimension needed is only L/2 and not L,
  // but we could just go to the right end with the previous even and
  // odd operators.
  D=1;
  for(int k=L-1;k>=L/2;k-=2){
    // odd site
    mwArray oddOp=reshape(identityMatrix(D),Indices(D*D,1))*reshape(identityMatrix(d),Indices(1,d*d));
    oddOp.reshape(Indices(D,D,d,d));
    oddOp.permute(Indices(3,1,4,2));
    mpo.setOp(k,new Operator(oddOp),true);
    // even site
    mwArray evenOp(Indices(d,D+1,d,D));evenOp.fillWithZero(); // redundant, but in case
    int oneFerm=FERMEVEN;
    int noFerm=1-FERMEVEN;
    for(int p=0;p<D;p++){
      evenOp.setElement(ONE_c,Indices(noFerm,p,noFerm,p));
      evenOp.setElement(ONE_c,Indices(oneFerm,p+1,oneFerm,p));
    }
    mpo.setOp(k-1,new Operator(evenOp),true);
    D++;    
  }
  int Dr=D;
  midOp=mwArray(Indices(Dl,L/2+1,Dr));
  for(int kl=0;kl<Dl;kl++)
    for(int kr=0;kr<Dr;kr++)
      midOp.setElement(ONE_c,Indices(kl,kl+kr,kr));
}

// Same for odd (count antifermions)
void prepareMPOantifermions(MPO& mpo,mwArray& midOp,int L){
  if(L%2!=0){
    cout<<"ERROR: The number of sites must be even"<<endl;
    exit(1);
  }
  mpo.initLength(L);
  // The even sites are simply identity in physical and in virtual
  // dimensions
  // TODO: check who is even and odd and which is the "occupied" state
  // in each, as I had different conventions. I assume now that 0 is
  // the occupied for odd (and start at 0, even)
  int D=1;
  for(int k=0;k<L/2;k+=2){
    // even site
    mwArray evenOp=reshape(identityMatrix(D),Indices(D*D,1))*reshape(identityMatrix(d),Indices(1,d*d));
    evenOp.reshape(Indices(D,D,d,d));
    evenOp.permute(Indices(3,1,4,2));
    mpo.setOp(k,new Operator(evenOp),true);
    // odd site
    mwArray oddOp(Indices(d,D,d,D+1));oddOp.fillWithZero(); // redundant, but in case
    int oneFerm=1-FERMEVEN;
    int noFerm=FERMEVEN;
    for(int p=0;p<D;p++){
      oddOp.setElement(ONE_c,Indices(noFerm,p,noFerm,p));
      oddOp.setElement(ONE_c,Indices(oneFerm,p,oneFerm,p+1));
    }
    mpo.setOp(k+1,new Operator(oddOp),true);
    D++;
  }
  int Dl=D; // the one in the middle
  // now the right part (actually, symmetric, could use reflections)
  // This is so that the bond dimension needed is only L/2 and not L,
  // but we could just go to the right end with the previous even and
  // odd operators.
  D=1;
  for(int k=L-1;k>=L/2;k-=2){
    // odd site
    mwArray oddOp(Indices(d,D+1,d,D));oddOp.fillWithZero(); // redundant, but in case
    int oneFerm=1-FERMEVEN;
    int noFerm=FERMEVEN;
    for(int p=0;p<D;p++){
      oddOp.setElement(ONE_c,Indices(noFerm,p,noFerm,p));
      oddOp.setElement(ONE_c,Indices(oneFerm,p+1,oneFerm,p));
    }
    mpo.setOp(k,new Operator(oddOp),true);
    D++;    
    // even site
    mwArray evenOp=reshape(identityMatrix(D),Indices(D*D,1))*reshape(identityMatrix(d),Indices(1,d*d));
    evenOp.reshape(Indices(D,D,d,d));
    evenOp.permute(Indices(3,1,4,2));
    mpo.setOp(k-1,new Operator(evenOp),true);
  }
  int Dr=D;
  midOp=mwArray(Indices(Dl,L/2+1,Dr));
  for(int kl=0;kl<Dl;kl++)
    for(int kr=0;kr<Dr;kr++)
      midOp.setElement(ONE_c,Indices(kl,kl+kr,kr));
}



#include <sstream>    

const string mpsfilename(int L,int D,double x,double mg,int level,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS_L"<<L<<"_m"<<mg<<"_x"<<x<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}
