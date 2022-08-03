
//#include "mwArray.h"
//#include "MPO.h"
#include "CoherentDissipationLiouvillianLR.h"
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
#include "SpinMPO.h"

#define MAXLEN 200  // Max nr of characters in a filename


/** 
    Program coherentProps: Compute the steady state of a certain
    Liouvillian
    as in \class <CoherentDissipationLiouvillianLR>,
    with g, gamma given as input parameters.

    To compile: make coherent

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
    + g (double) parameter \f$g/( L/4)\f$ of the Hamiltonian
    + gamma (double) parameter $gamma$ of the dissipation
    + scale (double) scale factor to multiply Lmpo
    + tol (double) tolerance asked from Contractor
    + mpsfile (string) file where the resulting MPS will be stored.
    + oldmpsfile [OPTIONAL] (string) name of another (binary) file, containing 
                    an earlier MPS to be used as starting point
    + outdir (string) directory where to write the (automatically
                    named) file(s) containing observables.

*/

using namespace std;
using namespace shrt;

/** Given a state, compute all the expectation values (on each
    position) of a given single body operator. The results are
    written to the ostream (could be a file, or cout). */
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
		 int D,double scal,double costT,const string oldMPSfile,int mode);
/** Vectorize an MPO to make an MPS out of it. If up==true (default)
    the first component of the double physical index will be the one
    which was going up.*/
void MPSfromMPO(const MPO& mpo,MPS& mps,bool up=true);

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
  double gamma=props.getDoubleProperty("gamma");
  double scal=props.getDoubleProperty("scale");
  double tol=props.getDoubleProperty("tol");
  const string newMPSdir=props.getProperty("mpsdir");
  const string oldMPSfile=props.getProperty("oldmpsfile");
  bool usingOldMPS=(!oldMPSfile.empty()&&oldMPSfile.length()>1);
  const string outdir=props.getProperty("outdir");
  double noise=props.getDoubleProperty("noise");
  if(noise<0.) noise=0.;
  int wur=props.getIntProperty("wurounds");
  if(wur<0) wur=0;
  int mode=props.getIntProperty("mode");
  const string jobsdir=props.getProperty("jobsdir");

  cout<<"oldMPSfile is "<<oldMPSfile<<endl;
  cout<<"its size is "<<oldMPSfile.size()<<endl;
  cout<<"its length is "<<oldMPSfile.length()<<endl;

  vector<int> allDs=props.getVectorIntProperty("allDs");
  int initL=props.getIntProperty("initL");

  vector<double> allGs=props.getVectorDoubleProperty("allGs");

  // I construct here the name of the text output files, where I will
  // save the results
  char fileNameC[MAXLEN];
  // For the single body observables (x,y,z, site after site, each in
  // one line)
  sprintf(fileNameC,"coherentSingleLR_%d_%1.2f_%1.2f_%d.txt",L,g,gamma,D);
  const string out1bodyfile=outdir+"/"+fileNameC;
  // For the two body observables (xx,yy,zz, site after site, first
  // for distance l=1, until N/2, each (xx, yy, zz) for each l
  // in one line, whichstarts with an identifier and the value of l
  sprintf(fileNameC,"coherentDoubleLR_%d_%1.2f_%1.2f_%d.txt",L,g,gamma,D);
  const string out2bodyfile=outdir+"/"+fileNameC;

  // Also the name of the MPS file
  sprintf(fileNameC,"MPSLR_%d_%1.2f_%1.2f_%d.mat",L,g,gamma,D);
  const string newMPSfile=newMPSdir+"/"+fileNameC;

  // // Open here the output file: for normal output Repeat the same for
  // // as many output files as required, but not necessarily here at the
  // // beginning!
  // ofstream* out1=new ofstream(out1bodyfile.data(),ios::app);
  // ofstream* out2=new ofstream(out2bodyfile.data(),ios::app);
  // if(!out1->is_open()){
  //   cout<<"Error: impossible to open file "<<out1bodyfile<<
  //     " for output"<<endl;
  //   exit(1);
  // }
  // if(!out2->is_open()){
  //   cout<<"Error: impossible to open file "<<out2bodyfile<<
  //     " for output"<<endl;
  //   exit(1);
  // }
  // cout<<"Writing output to files: \n\t"<<out1bodyfile
  //     <<"\n\t"<<out2bodyfile<<endl;

  //srandom(time(NULL));
  // Construct the Liouvillian
  CoherentDissipationLiouvillianLR model(L,scal*g*L*.25,scal*gamma);
  
  // Liouvillian
  const MPO& Lop=model.getLMPO();  
  // Adjoint Liouvillian
  const MPO& Ldagger=model.getLAdjointMPO();  

  // // To debug
  // if(L<=5)
  //   Lop.exportForMatlab("myLmpo.m");


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
  if(usingOldMPS){
    if(initL==-1){
      myGS.importMPS(oldMPSfile.data());
      myGS.increaseBondDimensionWithNoise(D,noise);
    }
    else{
      MPS aux(initL,D,d*d);
      aux.importMPS(oldMPSfile.data());
      myGS.stretchMPS(L,aux,L/2);
      cout<<"MPS was stretched from "<<initL<<" to "<<L<<endl;
    }
    //myGS.increaseBondDimension(D);
    cout<<"Using initial MPS from file "<<oldMPSfile<<endl;
  }
  else{
    if(D>1){
      myGS.setRandomState();
      cout<<"Starting from random MPS=> sweep first from center!"<<endl;
    }
    else{
      // The D=1 state is initialized according to the mode argument
      if(mode!=-1){ // Default is |0> state
	ProductState stat=p_zero;
	switch(mode){
	case 1:{
	  stat=p_one;
	  break;
	}
	case 2:{
	  stat=p_special;
	  break;
	}
	case 3:{
	  stat=p_maxent;
	  break;
	}
	default:{
	  stat=p_zero;
	}
	}
	if(mode<4)
	  myGS.setProductState(stat);
	else{ // Random cases
	  myGS.setRandomState();
	  if(mode==4){ // Force symmetry
	    int posL=0;int posR=L-1;
	    while(posL<posR){
	      myGS.setA(posR--,myGS.getA(posL++).getA());
	    }
	  }
	}
      }
      else
	cout<<"Starting the D=1 optimization with default (|0>) product state";
    }
  }

  MPO S2(L),Sz2(L),Sx2(L),Sy2(L);
  SpinMPO::getS2MPO(L,d,S2,0.);
  MPS auxS2(L,1,d*d);
  MPSfromMPO(S2,auxS2,false);
  MPS auxSx(L,1,d*d),auxSy(L,1,d*d),auxSz(L,1,d*d);
  SpinMPO::getSx2MPO(L,d,Sx2);
  SpinMPO::getSy2MPO(L,d,Sy2);
  SpinMPO::getSz2MPO(L,d,Sz2);
  MPSfromMPO(Sx2,auxSx,false);
  MPSfromMPO(Sy2,auxSy,false);
  MPSfromMPO(Sz2,auxSz,false);
  complex_t _sx2_=.25*contractor.contract(auxSx,myGS);
  complex_t _sy2_=.25*contractor.contract(auxSy,myGS);
  complex_t _sz2_=.25*contractor.contract(auxSz,myGS);
  
  cout<<"The very first initial state has <S^2>="
      <<contractor.contract(auxS2,myGS)
      <<" (when max N/2(N/2+1)="<<.5*L*(.5*L+1)<<" and sqrt(N)="
      <<sqrt(L)*(sqrt(L)+1)<<")"<<endl;
  cout<<"The individual terms are: "
      <<"<Sx2>="<<_sx2_
      <<"<Sy2>="<<_sy2_
      <<"<Sz2>="<<_sz2_
      <<" => sum="<<_sx2_+_sy2_+_sz2_<<endl;

  if(wur>0){
    cout<<"Applying "<<wur<<" symmetric warm up rounds"<<endl;
    for(int pp=0;pp<wur;pp++){
      for(int pos=0;pos<(L-1)/2+1;pos++){
	contractor.sweepPart(pos,pos+1,LdaggerL,D,myGS,offset); // one right
	contractor.sweepPart(L-1-pos,L-1-pos-1,LdaggerL,D,myGS,offset); // one left
      }
      lambda=real(contractor.contract(myGS,LdaggerL,myGS));
      cout<<"After symmetric sweep nr "<<pp+1<<" energy value "<<lambda<<endl;
    }
  }

  myGS.gaugeCond('R',1);
  myGS.gaugeCond('L',1);

  lambda=real(contractor.contract(myGS,LdaggerL,myGS));
  cout<<"Starting from an energy value "<<lambda<<endl;
  //if(abs(lambda)<1E-5)
  //offset=10*lambda;

  clock_t start=clock();
  contractor.findGroundState(LdaggerL,D,&lambda,myGS,offset);
  clock_t end=clock();

  cout<<"Finished GS search with lambda="<<lambda<<endl;

  complex_t  _s2_=contractor.contract(myGS,auxS2);
  cout<<"At the end, expectation value of S^2?="
      <<_s2_<<endl;
  _sx2_=.25*contractor.contract(auxSx,myGS);
  _sy2_=.25*contractor.contract(auxSy,myGS);
  _sz2_=.25*contractor.contract(auxSz,myGS);
  
  cout<<"The individual terms are: "
      <<"<Sx2>="<<_sx2_
      <<"<Sy2>="<<_sy2_
      <<"<Sz2>="<<_sz2_
      <<" => sum="<<_sx2_+_sy2_+_sz2_<<endl;



  // Now I have the GS, I will write it to a file
  myGS.exportMPS(newMPSfile.data());


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
  writeHeader(*out1,model,lambda,D,scal,cost,oldMPSfile,mode);
  *out1<<"% S2="<<_s2_<<endl;
  *out1<<"% Sx2="<<_sx2_<<endl;
  *out1<<"% Sy2="<<_sy2_<<endl;
  *out1<<"% Sz2="<<_sz2_<<endl;

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
  writeHeader(*out2,model,lambda,D,scal,cost,oldMPSfile,mode);
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


  if(!jobsdir.empty()){
    if(allDs.size()==0){
      cout<<"Do not know the list of Ds=> cannot predict next job. Exiting"<<endl;
      return 1;
    }
    // I need to determine nextD
    int nextD=D;
    if(errors&&(D==1)){
      cout<<"Error state found for D=1 => Try with a different initial state (mode)"<<endl;
      if(mode==0) mode=3;
      else if(mode==3) mode=4;
      else mode=-1; // to indicate this is not a mode change
    }
    else{ // or no error (=> launch with nextD), or D>1  
      if(!errors){
	vector<int>::iterator ptr=allDs.begin();
	while(ptr!=allDs.end()){
	  if(*ptr!=D) ptr++;
	  else{
	    if(ptr+1!=allDs.end())
	      nextD=*(ptr+1);
	    break;
	  }
	}
	cout<<"Next D (from "<<D<<") must be "<<nextD<<endl;
	if(nextD==D){
	  cout<<"There does not seem to be a larger D=> exiting with no job written"<<endl;
	  return errors?1:0;
	}
      }
    }
    int mpsSize=ceil(16*L*nextD*nextD*(d*d+2)*1E-6)+200;
    stringstream jobName;
    jobName<<"lrlr."<<L<<"."<<g<<"."<<gamma<<"."<<nextD;
    if(nextD==1) jobName<<"."<<mode<<"."<<wur;
    stringstream commandstr;
    commandstr<<"msub_modified_tqo097 -N "<<jobName.str()<<" -l h_vmem="<<mpsSize<<"M -raw ";
    commandstr<<argv[0]<<" "<<argv[1]<<" -L="<<L<<" -D="<<nextD<<" -g="<<g<<" -gamma="<<gamma<<" -scale="<<scal;
    if(nextD<10)
      commandstr<<" -tol="<<tol;
    else
      commandstr<<" -tol=1E-8";
    commandstr<<" -noise="<<noise<<" -outdir="<<outdir<<" -jobsdir="<<jobsdir;
    if(!errors)
      commandstr<<" -oldmpsfile="<<newMPSfile;
    else{
      //	commandstr<<" -mode="<<mode<<" -wurounds="<<wur;
      // If errors happened, I should try again, with another seed
      if(nextD>1||mode<0){
	// Guess neighboring g
	double nextG=g;
	vector<double>::iterator ptr=allGs.begin();
	while(ptr!=allGs.end()){
	  if(*ptr!=g) ptr++;
	  else{ // ptr on g
	    if(ptr!=allGs.begin()&&(g<1||ptr==allGs.end()))
	      nextG=*(ptr-1);
	    else
	      nextG=*(ptr+1);
	  }
	}
	cout<<"Next g (from "<<g<<") must be "<<nextG<<endl;

	// build MPS reference filename, from g value before
	sprintf(fileNameC,"MPSLR_%d_%1.2f_%1.2f_%d.mat",L,nextG,gamma,nextD);
	commandstr<<" -oldmpsfile="<<newMPSdir<<"/"<<fileNameC;
      }

    }
    commandstr<<" -mpsdir="<<newMPSdir;
    
    string jobfile=jobsdir+"/"+jobName.str()+".sh";
    ofstream* ojob=new ofstream(jobfile.data());
    if(!ojob->is_open()){
      cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
      cout<<commandstr.str()<<endl;
      cout<<endl;
      // Do not exit, because I still need to save the MPS
    }
    else{
      *ojob<<commandstr.str()<<endl;
      ojob->close();
      if(errors&&D>1){
	cout<<"WARNING!!!! Next job launched ( "<<jobName<<" ) although errors in the intial state!!"<<endl;
	return 212;
      }
    }
    delete ojob;    
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
	  <<"+i*"<<imag(valK)/traceRho<<", trace(rho)="<<traceRho<<", and nasty tensor is "<<rho.getA(k).getA()<<endl;
      result=-1;
    }
  }
  os<<endl;
  os<<"% Trace="<<traceRho;
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
		 int D,double scal,double costT,const string oldMPSfile,int mode){
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
  else{
    if(D>1)
      out<<"% Using random initial state"<<endl;
    else
      out<<"% Using independent initial state, mode "<<mode<<endl;
  }
  out<<"% Lowest eigenvalue found to be "<<lambda0<<endl;
}



void MPSfromMPO(const MPO& mpo,MPS& mps,bool up){
  int L=mpo.getLength();
  // vector<int> du=mpo.getDimensions();
  // vector<int> dd=mpo.getOriginalDimensions();
  // int D=mpo.getOp(0).getDr();

  // vector<int> newDims(du);
  // for(int k=0;k<L;k++) newDims[k]*=dd[k];

  mps=MPS(L,1,1);

  for(int k=0;k<L;k++){
    mwArray aux=mpo.getOp(k).getFullData();
    Indices dims=aux.getDimensions();
    if(up)
      aux.permute(Indices(1,3,2,4));
    else
      aux.permute(Indices(3,1,2,4));
    aux.reshape(Indices(dims[0]*dims[2],dims[1],dims[3]));
    mps.replaceSite(k,aux,0); // brute force
  }

}
