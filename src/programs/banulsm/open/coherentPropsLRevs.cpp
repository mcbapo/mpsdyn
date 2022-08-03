
//#include "mwArray.h"
//#include "MPO.h"
#include "CoherentDissipationLiouvillianLR.h"
//#include "Operator.h"
#include "Contractor.h"
#include "Properties.h"
#include "DoubleOperator.h"
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
    Program coherentPropsLRevs: Compute some expectation values from
    the final result obtained by coherentPropsLR. It receives as
    oldmpsfile argument the corresponding mps file.

    To compile: make coherentLR2

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
    + gamma (double) parameter $gamma$ of the dissipation
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

complex_t getTraceRho(const MPS& rho);

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
  const string oldMPSfile=props.getProperty("oldmpsfile");
  const string outdir=props.getProperty("outdir");
  const string jobsdir=props.getProperty("jobsdir");

  // cout<<"oldMPSfile is "<<oldMPSfile<<endl;
  // cout<<"its size is "<<oldMPSfile.size()<<endl;
  // cout<<"its length is "<<oldMPSfile.length()<<endl;

  // I construct here the name of the text output files, where I will
  // save the results
  char fileNameC[MAXLEN];
  // For the single body observables (x,y,z, site after site, each in
  // one line)
  sprintf(fileNameC,"extraObservablesLR_%d_%1.2f_%1.2f_%d.txt",L,g,gamma,D);
  const string out1bodyfile=outdir+"/"+fileNameC;

  // And directly try for the GS!
  Contractor& contractor=Contractor::theContractor();

  MPS myGS(L,D,d*d);
  myGS.importMPS(oldMPSfile.data());

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

  // Open here the output file: for normal output Repeat the same for
  // as many output files as required, but not necessarily here at the
  // beginning!
  ofstream* out1=new ofstream(out1bodyfile.data(),ios::app);
  if(!out1->is_open()){
    cout<<"Error: impossible to open file "<<out1bodyfile<<
      " for output"<<endl;
    exit(1);
  }
  cout<<"Writing output to files: \n\t"<<out1bodyfile<<endl;

  // Now compute observables and write results to files
  *out1<<"% L="<<L<<endl;
  *out1<<"% g="<<g<<endl;
  *out1<<"% gamma="<<gamma<<endl;
  *out1<<"% D="<<D<<endl;
  *out1<<"% trace="<<getTraceRho(myGS)<<endl;
  // Normalize:
  myGS.setNormFact(myGS.getNormFact()/abs(getTraceRho(myGS)));
  *out1<<"% <S2>="<<contractor.contract(auxS2,myGS)<<endl;

  // Now project out
  // 1) Construct the projector onto pair of spins 1
  MPO P1e(L),P1o(L);
  SpinMPO::getProjectorTotalSPair(L,d,P1e,true);
  SpinMPO::getProjectorTotalSPair(L,d,P1o,false);
  MPS auxP1e(L,4,d*d),auxP1o(L,4,d*d);
  MPSfromMPO(P1e,auxP1e,false);
  MPSfromMPO(P1o,auxP1o,false);
  // Probabilities ofP1e and P1o
  *out1<<"% Before projecting, the probability should be (even)="
       <<contractor.contract(auxP1e,myGS)<<endl;
  *out1<<"% Before projecting, the probability should be (odd)="
       <<contractor.contract(auxP1o,myGS)<<endl;

  // 2) compose each of them to have the double one
  // The second one should be transposed, but in this case it remains
  // the same (each Sy term gets a sig, but coming in pairs, the
  // product remains constant)
  MPO ProjSe(L),ProjSo(L);
  // Should keep pointers and avoid new Ops, but I do not need so much
  // memory for this program
  for(int k=0;k<L;k++){
    ProjSe.setOp(k,new DoubleOperator(P1e.getOp(k).getFullData(),
				      P1e.getOp(k).getFullData()),true);
    ProjSo.setOp(k,new DoubleOperator(P1o.getOp(k).getFullData(),
				      P1o.getOp(k).getFullData()),true);
  }
  // Find best approximation of the ProjSe one
  MPS aux(myGS);
  contractor.optimize(ProjSe,myGS,aux);
  *out1<<"% After projecting with the even ones trace="<<
    getTraceRho(aux)<<endl;
  *out1<<"% S2="<<contractor.contract(auxS2,aux)<<endl;
  contractor.optimize(ProjSo,aux,myGS);
  *out1<<"% After projecting with the odd ones trace="<<
    getTraceRho(myGS)<<endl;
  *out1<<"% S2="<<contractor.contract(auxS2,myGS)<<endl;

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


}


complex_t getTraceRho(const MPS& rho){
  int L=rho.getLength();
  MPS aux(L,1,d*d);
  mwArray basic=identityMatrix(d); // basic piece of the MPO
  basic.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++){
    aux.setA(k,basic);
  }
  Contractor& contractor=Contractor::theContractor();
  return (contractor.contract(rho,aux));
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
