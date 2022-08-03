#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "JoinedOperator.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0


using namespace shrt;
using namespace std;

complex_t computeTrace(const MPS& Omps);

const string mpsfilename(int M,int D,const string mpsdir);
const string jobfilename(int M,int D,double gP,double hP,const string jobsdir);
void preparejobfile(Properties& props);
void constructMPOrange(MPO& mpoL,int range,int len);
int d=2;

int Ds[]={2,4,10,20,40,60,80,100,120,140,160,180,200};
int MAXPTR=13; // nr elems Ds

/** Compute numerically the weights of different range in the best MPO
    for the Hamiltonian problem.  Uses the results (MPS) saved by
    slowMPOham and computes contractions with appropriate projectors.    
*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  //  string directory="";
  // if(argc>2){
  //   const char* indir=argv[++cntr];
  //   directory=string(indir)+"/";
  // }
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


  int M=props.getIntProperty("M");
  int D=props.getIntProperty("D");
  bool allDs=false;
  if(D<0){
    // none provided: run all of them and replace file!
    cout<<"No valid D provided, running them all"<<endl;
    allDs=true;
  }
  else{
    cout<<"Running with single D="<<D<<endl;
  }
  //  const char* outfnameG=argv[++cntr];
  // const char* outfnameL=argv[++cntr];
  string mpsdir=props.getProperty("mpsdir");
  string jobsdir=props.getProperty("jobsdir"); 
  string outfile=props.getProperty("weightsfile");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-5;
  double gP=props.getDoubleProperty("gP");
  double hP=props.getDoubleProperty("hP");

  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  cout<<"Initialized Contractor"<<endl;
  //  srandom(rSeed);

  ofstream* out;
  if(!allDs){
    cout<<"Opening "<<outfile<<" to append"<<endl;
    out=new ofstream(outfile.data(),ios::app);
  }
  else{
    cout<<"Oppening "<<outfile<<" anew"<<endl;
    out=new ofstream(outfile.data());
  }
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output"<<endl;
    exit(1);
  }
  bool done=false;
  int ptr=0;
  while(!done){
    if(allDs){
      D=Ds[ptr];
    }
    string fileMPS=mpsfilename(M,D,mpsdir);
    // 1) Read the already computed state
    MPS Omps(M,D,d*d);
    if(file_exists(fileMPS)){
      Omps.importMPS(fileMPS.data());
      //cout<<"Imported MPS from file "<<fileMPS<<endl;
      //cout<<"Omps="<<Omps<<endl;
      mwArray R(Indices(M+1,1));
      // The 0-term is proportional to the trace and should then be zero (check normalization!)
      complex_t trO=computeTrace(Omps);
      R.setElement(1./pow(d,M)*conjugate(trO)*trO,Indices(0,0));
      *out<<0<<"\t"<<1./pow(d,M)*conjugate(trO)*trO<<"\t"<<D<<endl;
      cout<<0<<"\t"<<1./pow(d,M)*conjugate(trO)*trO<<"\t"<<D<<endl;
      // 2) For each range, from 1 to M, compute the weight squared
      for(int r=1;r<=M;r++){
	MPO mpoR(M);
	constructMPOrange(mpoR,r,M); // the projector
	//cout<<"Constructed MPO for range "<<r<<": "<<mpoR<<endl;
	complex_t value=contractor.contract(Omps,mpoR,Omps);
	R.setElement(value,Indices(r,0));	  
	*out<<r<<"\t"<<real(value)<<"\t"<<D<<endl;
	cout<<r<<"\t"<<real(value)<<"\t"<<D<<endl;
      }
    }
    else{  // no tmp file
      cout<<"There is no MPS file "<<fileMPS<<" for M="<<M<<", D="<<D<<endl;
      // Pass the required info
      // props.addProperty("command",argv[0]);
      // props.addProperty("propFile",argv[1]);
      // preparejobfile(jobsdir,props);
      if(~allDs)      
	exit(15);
    }
    if(allDs){
      ptr++;
      if(ptr>=MAXPTR) done=true;
    }
    else done=true; // only once
  }

  out->close();delete out;
  // if(D<260){
  //   // Check if result was converged
  //   int newD=D+20; // in general
  //   // other cases
  //   if(D<20){
  //     if(D==1) newD=2;
  //     if(D==2) newD=4;
  //     if(D==4) newD=10;
  //     if(D==10) newD=D+10;
  //   }
  //   // Once everything is done, launch the next job, for next D
  //   // First, manipulate the properties, to include the new M 
  //   props.addProperty("command",argv[0]);
  //   props.addProperty("propFile",argv[1]);
  //   props.addProperty("D",newD);
  //   preparejobfile(jobsdir,props);
  // }


}


void preparejobfile(Properties& props){
  int M=props.getIntProperty("M");
  int D=props.getIntProperty("D");
  //  const char* outfnameG=argv[++cntr];
  // const char* outfnameL=argv[++cntr];
  string mpsdir=props.getProperty("mpsdir");
  string jobsdir=props.getProperty("jobsdir"); 
  string outfile=props.getProperty("weightsfile");
  string command=props.getProperty("command");
  string propFile=props.getProperty("propFile");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-5;
  double gP=props.getDoubleProperty("gP");
  double hP=props.getDoubleProperty("hP");
  stringstream commandstr;
  commandstr<<command;
  commandstr<<" "<<propFile;
    commandstr<<" -gP="<<gP<<" -hP="<<hP;
    commandstr<<" -M="<<M<<" -weightsfile="<<outfile;
    commandstr<<" -jobsdir="<<jobsdir;
    commandstr<<" -tol="<<tol<<" -D="<<D;
    commandstr<<" -mpsdir="<<mpsdir;
    const string jobfile=jobfilename(M,D,gP,hP,jobsdir);
    ofstream* ojob=new ofstream(jobfile.data());
    if(!ojob->is_open()){
      cout<<"WARNING!!!! Could not open the JOB file "<<jobfile<<". Launch by hand with the following command:"<<endl;
      cout<<commandstr.str()<<endl;
      cout<<endl;
      exit(15);
    }
    else{
      *ojob<<"#!/bin/bash"<<endl;
      *ojob<<"msub_modified_tqo097 ";
      //if(D+20>60) *ojob<<" -P 4";
      *ojob<<" -N slowHloc."<<hP<<"."<<M<<"."<<D;
      *ojob<<" -raw ";
      *ojob<<commandstr.str()<<endl;
      ojob->close();
    }
    delete ojob;
}


const string mpsfilename(int M,int D,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS"<<"_M"<<M<<"_D"<<D<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

const string jobfilename(int M,int D,double gP,double hP,const string jobsdir){
  stringstream s;
  s<<jobsdir<<"/job_Hloc_M"<<M<<"_D"<<D<<"_g"<<gP<<"_h"<<hP;  
  s<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}


void constructMPOrange(MPO& mpoL,int range,int len){
  //cout<<"In MPO(L="<<len<<"), range="<<range<<endl;
  if(len!=mpoL.getLength())
    mpoL.initLength(len);
  if(range>len){
    cout<<"ERROR! range cannot be larger than MPO!"<<endl;
    exit(1);
  }
  // Basic operators, construct them once, then reuse
  static bool init=false;
  static mwArray Z(Indices(3,d*d,d*d));
  if(!init){
    mwArray O0=.5*identityMatrix(d*d);
    O0.reshape(Indices(d,d,d,d));
    O0.permute(Indices(2,4,1,3));
    O0.reshape(Indices(d*d,d*d));
    mwArray O1=identityMatrix(d*d);
    // construct the Z mwArray
    for(int i1=0;i1<d*d;i1++)
      for(int i2=0;i2<d*d;i2++){
	Z.setElement(O0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
	Z.setElement(O1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
	Z.setElement(O1.getElement(Indices(i1,i2))-O0.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      }
    Z.reshape(Indices(3,d*d*d*d));
    init=true;
    //cout<<"Initialized the Z operator for the MPO"<<endl;
  }

  // DEBUGGING ONLY: I build the MPS alone to check coeffs
  MPS aux(len,range+1,3);

  // MPO terms. Special case, k=1=L
  int D=range+1;int Dl(D),Dr(D);
  for(int p=0;p<len;p++){
    if(p==0) Dl=1; else Dl=D;
    if(p==len-1) Dr=1; else Dr=D;
    mwArray C(Indices(Dl,Dr,3));
    if(len==1){ // special case: just one term
      C.setElement(ONE_c,Indices(Dl,Dr,2));
    }
    else{
      if(p<len-range)
	C.setElement(ONE_c,Indices(0,0,0));
      if(p>range-1)
	C.setElement(ONE_c,Indices(Dl-1,Dr-1,0));
      if(p<=len-range)
	if(p<len-1)
	  C.setElement(ONE_c,Indices(0,1,2)); // left operator
	else
	  C.setElement(ONE_c,Indices(0,0,2)); // left operator is also last
      for(int ik=1;ik<range-1;ik++){
	if(p>=ik&&(len-1-p)>=range-ik-1) // place for the rest of the term on both sides
	  C.setElement(ONE_c,Indices(ik,ik+1,1)); // in between
      }
      if(p>=range-1)
	C.setElement(ONE_c,Indices(range-1,Dr-1,2)); // in between
    }
    aux.replaceSite(p,permute(C,Indices(3,1,2)));
    //cout<<"Coefficient of site "<<p<<", C="<<C<<endl;
    C.reshape(Indices(Dl*Dr,3));
    C.multiplyRight(Z);
    C.reshape(Indices(Dl,Dr,d*d,d*d));
    C.permute(Indices(3,1,4,2));
    mpoL.setOp(p,new Operator(C),true); // no reusing operators!
  }
  // stringstream nameMPS;
  // nameMPS<<"mpsForRange"<<range<<"_L"<<len<<".m";
  // aux.exportForMatlab(nameMPS.str().data());
}

 
complex_t computeTrace(const MPS& Omps){
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  int len=Omps.getLength();
  MPS aux(len,1,d*d);
  for(int k=0;k<len;k++){
    aux.setA(k,id2);
  }
  Contractor& contractor=Contractor::theContractor();
  return contractor.contract(Omps,aux);
}
