
//#include "mwArray.h"
//#include "MPO.h"
#include "CoherentDissipationLiouvillian.h"
//#include "Operator.h"
#include "Contractor.h"
#include "Properties.h"
//#include "JoinedOperator.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <stdio.h>

#define MAXLEN 200  // Max nr of characters in a filename


/** 
    Program coherentAdiab: Compute the steady state of a certain
    Liouvillian
    as in \class <CoherentDissipationLiouvillian>,
    with g, gamma given as input parameters.
    In this version, another MPS computed for the same D and slightly larger g 
    is used as input, and when the program finishes, it starts a new job for 
    a smaller g.

    To compile: make coherentAd

    Receives as argument a Properties file, such as
    config/excitProp.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.
    To replace some of the values in the file by a command line
    option, the additional arguments have to be of the form
    -propName=propValue

    The properties that will be required:
    + L (int) length of the chain
    + D (int) maximum bond dimension
    + g (double) parameter $g$ of the Hamiltonian
    + deltaG (double) change in the parameter $g$ of the Hamiltonian 
                      for the next job
    + gamma (double) parameter $gamma$ of the dissipation
    + scale (double) scale factor to multiply Lmpo
    + tol (double) tolerance asked from Contractor
    + mpsfile (string) file where the resulting MPS will be stored.
    + oldmpsfile (string) name of another (binary) file, containing 
                          an earlier MPS to be used as starting point
    + outdir (string) directory where to write the (automatically
                      named) file(s) containing observables.
*/

using namespace std;
using namespace shrt;

/** Given a state, compute all the expectation values (on each
    position) of a given single body operator. The results are
    written to the ostream (could be a file, or cout). 
    Returns 0 if values are ok, -1 if some site is non physical
*/
int computeSingleBodyEV(const MPS& rho,const mwArray& op1,ostream& os);

/** Same as before, but compute the expectation value of the product
    of op1 acting on one site and op2 acting on the site plus l.  If
    maxLen!=0 is provided, it writes up to maxLen elements in the
    line, padding with 0 when needed.  
*/
void computeTwoBodyEV(const MPS& rho,const mwArray& op1,
		      const mwArray& op2,int l,ostream& os,
		      int maxLen=0);

/** Prepare header for a file */
void writeHeader(ofstream& out,const CoherentDissipationLiouvillian& L,
		 double lambda0,
		 int D,double scal,double costT,const string oldMPSfile);

/** Force a site in the MPS to be Hermitian */
void makeHermitian(MPS& mps,int pos,int d);

bool file_exists(const string filename);

/** Global matrices */
int d=2;
complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
mwArray sigX(Indices(d,d),dataX);
complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
mwArray sigY(Indices(d,d),dataY);
complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
mwArray sigZ(Indices(d,d),dataZ);

//double noise=1E-60;

int main(int argc,const char* argv[]){
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


  int L=props.getIntProperty("L");
  int D=props.getIntProperty("D");
  double g=props.getDoubleProperty("g");
  double deltaG=props.getDoubleProperty("deltaG");
  double gamma=props.getDoubleProperty("gamma");
  double scal=props.getDoubleProperty("scale");
  double tol=props.getDoubleProperty("tol");
  const string newMPSdir=props.getProperty("mpsdir");
  string oldMPSfile=props.getProperty("oldmpsfile");
  bool usingOldMPS=(!oldMPSfile.empty()&&oldMPSfile.length()>1);
  const string outdir=props.getProperty("outdir");
  double noise=props.getDoubleProperty("noise");
  if(noise<0.) noise=0.;
  const string jobsdir=props.getProperty("jobsdir");
  int freqSv=props.getIntProperty("saveTmp");
  if(freqSv<0) freqSv=0;

  //  cout<<"oldMPSfile is "<<oldMPSfile<<endl;
  // cout<<"its size is "<<oldMPSfile.size()<<endl;
  // cout<<"its length is "<<oldMPSfile.length()<<endl;

  //  int initL=props.getIntProperty("initL");

  if(!usingOldMPS){
    cout<<"Error: the adiabatic version of coherentAdiab cannot run without reference"<<endl;
    exit(1);
  }
  if(oldMPSfile=="auto"){
    cout<<"Error: the adiabatic version of coherentAdiab needs explicit reference"<<endl;
    exit(1);
  }
  
  {
    ifstream test(oldMPSfile.data());
    if(test)
      cout<<"Found reference file!"<<endl;  
    else{
      cout<<"Error: Cannot open refernce file "<<oldMPSfile<<endl;
      exit(1);
    } 
  }

  // I construct here the name of the text output files, where I will
  // save the results
  char fileNameC[MAXLEN];
  // For the single body observables (x,y,z, site after site, each in
  // one line)
  sprintf(fileNameC,"coherentSingle_%d_%1.2f_%1.2f_%d.txt",L,g,gamma,D);
  const string out1bodyfile=outdir+"/"+fileNameC;
  // For the two body observables (xx,yy,zz, site after site, first
  // for distance l=1, until N/2, each (xx, yy, zz) for each l
  // in one line, whichstarts with an identifier and the value of l
  sprintf(fileNameC,"coherentDouble_%d_%1.2f_%1.2f_%d.txt",L,g,gamma,D);
  const string out2bodyfile=outdir+"/"+fileNameC;

  // Also the name of the MPS file
  sprintf(fileNameC,"MPS_%d_%1.2f_%1.2f_%d.mat",L,g,gamma,D);
  const string newMPSfile=newMPSdir+"/"+fileNameC;
  const string tmpMPSfile=newMPSfile+".tmp"; // temporary place

  // Construct the Liouvillian
  CoherentDissipationLiouvillian model(L,scal*g,scal*gamma);
  
  // Liouvillian
  const MPO& Lop=model.getLMPO();  
  // Adjoint Liouvillian
  const MPO& Ldagger=model.getLAdjointMPO();  

  // And now construct a Joined MPO for the L+ L
  const MPO* ptrs[2]={&Ldagger,&Lop};
  MPO LdaggerL(L);
  MPO::join(2,ptrs,LdaggerL);

  // And directly try for the GS!
  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  if(D<=12)
    contractor.setEigenSolver(fulleig);
  contractor.setConvTol(tol);

  MPS myGS(L,D,d*d);
  double lambda;
  double offset=0.; // Contractor will set it alone

  // Initialize the MPS with something already existing, maybe?
  myGS.importMPS(oldMPSfile.data());
  myGS.gaugeCond('R',1);
  myGS.gaugeCond('L',1);

  lambda=real(contractor.contract(myGS,LdaggerL,myGS));
  cout<<"Starting from an energy value "<<lambda<<endl;
  //if(abs(lambda)<1E-5)
  //offset=10*lambda;

  clock_t start=clock();
  if(freqSv)
    contractor.findGroundState(LdaggerL,D,&lambda,myGS,offset,1,tmpMPSfile,freqSv);
  else
    contractor.findGroundState(LdaggerL,D,&lambda,myGS,offset);
  clock_t end=clock();

  cout<<"Finished GS search with lambda="<<lambda<<endl;

  // Now I have the GS, I will write it to a file
  myGS.exportMPS(newMPSfile.data());

  // And delete the temporary file, if present
  if(file_exists(tmpMPSfile))
    remove(tmpMPSfile.data());

  // I write out some results to output file, and close it (instead, I
  // could start computing correaltions afterwards) Actually, I could
  // open out just now, because it does not need to be open all during
  // the calculation of Contractor
  double cost=(double)(end-start)/CLOCKS_PER_SEC;
  
  // Open here the output file: for normal output Repeat the same for
  // as many output files as required, but not necessarily here at the
  // beginning!
  ofstream* out1=new ofstream(out1bodyfile.data(),ios::app);
  if(!out1->is_open()){
    cout<<"Error: impossible to open file "<<out1bodyfile<<
      " for output"<<endl;
    exit(1);
  }
  cout<<"Writing output to files: \n\t"<<out1bodyfile
      <<"\n\t"<<out2bodyfile<<endl;
  // Now compute observables and write results to files
  writeHeader(*out1,model,lambda,D,scal,cost,oldMPSfile);

  int errors=0;
  // Compute now expectation values of sigmaX,Y,Z
  *out1<<0<<"\t"<<1<<"\t";
  errors=computeSingleBodyEV(myGS,sigX,*out1);
  *out1<<endl;
  *out1<<0<<"\t"<<2<<"\t";
  errors=computeSingleBodyEV(myGS,sigY,*out1);
  *out1<<endl;
  *out1<<0<<"\t"<<3<<"\t";
  errors=computeSingleBodyEV(myGS,sigZ,*out1);
  *out1<<endl;
  out1->close();
  delete out1;

  ofstream* out2=new ofstream(out2bodyfile.data(),ios::app);
  if(!out2->is_open()){
    cout<<"Error: impossible to open file "<<out2bodyfile<<
      " for output"<<endl;
    exit(1);
  }
  writeHeader(*out2,model,lambda,D,scal,cost,oldMPSfile);
  *out2<<setprecision(10);
  *out2<<"% Every line starts with an integer: 0 for the lowest eigenvalue "
       <<"(would be one for the next one). Then another integer: 1 for XX, "
       <<"2 for YY, 3 for ZZ correlators; the next value is the distance l "
       <<"(so correlators  <X(i)X(i+3)> will be in a line starting with 0 1"
       <<" 3 and then a series of double values from i=1 to the end"<<endl;
  // Compute now expectation values of sigmaXX,YY,ZZ at distance 1-3
  int maxLen=L-1; // the largest number of terms
  for(int l=1;l<=L-1;l++){
    *out2<<0<<"\t"<<1<<"\t"<<l<<"\t";
    computeTwoBodyEV(myGS,sigX,sigX,l,*out2,maxLen);
    *out2<<endl;
    *out2<<0<<"\t"<<2<<"\t"<<l<<"\t";
    computeTwoBodyEV(myGS,sigY,sigY,l,*out2,maxLen);
    *out2<<endl;
    *out2<<0<<"\t"<<3<<"\t"<<l<<"\t";
    computeTwoBodyEV(myGS,sigZ,sigZ,l,*out2,maxLen);
    *out2<<endl;
  }
  out2->close();
  delete(out2);

  if(g-deltaG>=0&&g-deltaG<=4)
    if(!jobsdir.empty()){
      if(errors){
	cout<<"Error happened: adiabatic didn't work at g="<<g<<", L="<<L<<endl;
	//return 111;
      }
      int mpsSize=ceil(16*L*(D*D*(d*d+2)+2*25*d*d*d*d+d*d)*1E-6)+200;
      if(L>=80&&mpsSize<500) mpsSize=500;
      stringstream jobName;
      jobName<<"ad."<<L<<"."<<g-deltaG<<"."<<D;
      stringstream commandstr;
      commandstr<<"msub_modified_tqo097 -N "<<jobName.str()<<" -l h_vmem="<<mpsSize<<"M -raw ";
      commandstr<<argv[0]<<" "<<argv[1]<<" -L="<<L<<" -D="<<D<<" -g="<<g-deltaG<<" -deltaG="<<deltaG<<" -gamma="<<gamma<<" -scale="<<scal;
      commandstr<<" -tol="<<tol;
      commandstr<<" -noise="<<noise<<" -outdir="<<outdir<<" -jobsdir="<<jobsdir;
      commandstr<<" -oldmpsfile="<<newMPSfile;
      commandstr<<" -mpsdir="<<newMPSdir;
      
      string jobfile=jobsdir+"/"+jobName.str()+".sh";
      ofstream* ojob=new ofstream(jobfile.data());
      if(!ojob->is_open()){
	cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
	cout<<commandstr.str()<<endl;
	cout<<endl;
      }
      else{
	*ojob<<commandstr.str()<<endl;
	ojob->close();
      }
      delete ojob;    
    }
    else{
      cout<<"g=0 reached: no further job prepared for this D, L"<<endl;
    }
}


int computeSingleBodyEV(const MPS& rho,const mwArray& op1,ostream& os){
  bool result=0;
  int L=rho.getLength();
  MPS aux(L,1,d*d);
  mwArray basic=identityMatrix(d); // basic piece of the MPO
  basic.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++){
    aux.setA(k,basic);
  }

  Contractor& contractor=Contractor::theContractor();
  double traceRho=real(contractor.contract(rho,aux));
  if(abs(traceRho)<1E-6){
    cout<<"** WARNING: Non physical state detected, trace(rho)="<<traceRho<<endl;
    result=-1;    
  }

  // WARNING: No check of dimensions done!!
  mwArray op1_=reshape(op1,Indices(d*d,1,1));
  op1_.conjugate(); // It gets conjugated in the scalar product! 

  // Now, site by site, replace the identity by the operator,
  // contract, write the result to os, and replace the operator by the
  // identity again
  for(int k=0;k<L;k++){
    aux.setA(k,op1_); // operator on site k
    complex_t valK=contractor.contract(rho,aux);
    os<<real(valK)/traceRho<<"\t";
    aux.setA(k,basic); // operator removed
    if(abs(real(valK)/traceRho)>1){
      cout<<"** WARNING: Non physical state detected, site "<<k<<" shows EV="<<real(valK)/traceRho
	  <<"+i*"<<imag(valK)/traceRho<<", trace(rho)="<<traceRho<<endl;
      //<<", and nasty tensor is "<<rho.getA(k).getA()<<endl;
      // If D is very small, I allow for non-physical values that arre not too bad:
      if(rho.getBond()<=2&&abs(abs(real(valK)/traceRho)-1)<1E-2) result=0;
      else result=-1;
    }
  }
  return result;
}

void computeTwoBodyEV(const MPS& rho,const mwArray& op1,
		      const mwArray& op2,int l,ostream& os,
		      int maxLen){
  int L=rho.getLength();
  MPS aux(L,1,d*d);
  mwArray basic=identityMatrix(d); // basic piece of the MPO
  basic.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++){
    aux.setA(k,basic);
  }

  Contractor& contractor=Contractor::theContractor();
  double traceRho=real(contractor.contract(rho,aux));

  // WARNING: No check of dimensions done!!
  mwArray op1_=reshape(op1,Indices(d*d,1,1));
  op1_.conjugate(); // It gets conjugated in the scalar product! 
  mwArray op2_=reshape(op2,Indices(d*d,1,1));
  op2_.conjugate(); // It gets conjugated in the scalar product! 

  // Now, site by site, replace the identity at k by the operator op1,
  // the identity in k+l by op2,
  // contract, write the result to os, and replace the operators by the
  // identity again
  int cnt=0; // nr of terms written
  for(int k=0;k+l<L;k++){
    aux.setA(k,op1_); // operator on site k
    aux.setA(k+l,op2_); // operator on site k+l
    complex_t valK=contractor.contract(rho,aux);
    os<<real(valK)/traceRho<<"\t";
    cnt++;
    aux.setA(k,basic); // operator1 removed
    aux.setA(k+l,basic); // operator2 removed
  }
  if(maxLen>0)
    while(cnt<maxLen){
      os<<0.<<"\t";
      cnt++;
    }

}

void writeHeader(ofstream& out,const CoherentDissipationLiouvillian& L,
		 double lambda0,
		 int D,double scal,double costT,const string oldMPSfile){
  out<<"% L="<<L.getLength()<<endl;
  out<<"% g="<<L.getG()<<endl;
  out<<"% gamma="<<L.getGamma()<<endl;
  out<<"% D="<<D<<endl;
  out<<"% scale="<<scal<<endl;
  Contractor& contr=Contractor::theContractor();
  out<<"% Contractor tol="<<contr.getConvTol()<<endl;
  out<<"% Contractor solver="<<contr.getSolverName()<<endl;
  out<<"% Time (s) needed = "<<costT<<endl;
  if(!oldMPSfile.empty())
    out<<"% Using initial state from "<<oldMPSfile<<endl;
  out<<"% Lowest eigenvalue found to be "<<lambda0<<endl;
}

void makeHermitian(MPS& mps,int pos,int d){
  mwArray aux=mps.getA(pos).getA();
  Indices dims=aux.getDimensions();
  if(dims[0]!=d*d){
    cout<<"ERROR: Cannot make Hermitian this tensor (pos "<<pos<<") "<<aux<<endl;
  }
  aux.reshape(Indices(d,d,dims[1]*dims[2]));
  aux=.5*(aux+conjugate(permute(aux,Indices(2,1,3))));
  aux.reshape(dims);
  mps.setA(pos,aux);
}

bool file_exists(const string filename){
  ifstream ifile(filename.data());
  return ifile;
}
