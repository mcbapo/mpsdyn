
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
    Program coherentEffectiveDS: Read the MPS computed by
    coherentDegeneracy
    and try to extract info about phyiscal steady states
    for Liouvillian 
    as in \class <CoherentDissipationLiouvillian>,
    with g, gamma given as input parameters,
    in the region where we expect degeneracy.

    To compile: make coherentDS

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
    + scale (double) scale factor to multiply Lmpo
    + tol (double) tolerance asked from Contractor
    + nrlev (int) number of levels to read (four)
    + mpsdir (string) name of the directory where 
                 MPS files are stored 
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
		 int D,double scal,double costT,const string oldMPSfile,int mode);

/** Force a site in the MPS to be Hermitian */
void makeHermitian(MPS& mps,int pos,int d);


/** Compute matrix form in the subspace */
void computeMatrix(const vector<MPS*>& states,const MPO& oper,mwArray& result);

/** Special overlap matrix */
void computeOverlapMatrix(const vector<MPS*>& states,mwArray& result,int d);


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
  const string oldMPSdir=props.getProperty("mpsdir");
  //const string oldMPSdir=props.getProperty("oldmpsdir");
  //bool usingOldMPS=(!oldMPSdir.empty()&&oldMPSdir.length()>1);
  const string outdir=props.getProperty("outdir");
  //double noise=props.getDoubleProperty("noise");
  //if(noise<0.) noise=0.;
  //int wur=props.getIntProperty("wurounds");
  //if(wur<0) wur=0;
  //int mode=props.getIntProperty("mode");
  ////  int makeH=1;
  //const string jobsdir=props.getProperty("jobsdir");
  int nrLevs=props.getIntProperty("numberoflevels");
  if(nrLevs<0) nrLevs=4;

  //  cout<<"oldMPSfile is "<<oldMPSfile<<endl;
  // cout<<"its size is "<<oldMPSfile.size()<<endl;
  // cout<<"its length is "<<oldMPSfile.length()<<endl;

  //vector<int> allDs=props.getVectorIntProperty("allDs");
  //cout<<"Read vector allDs="<<allDs<<endl;
  //int initL=props.getIntProperty("initL");

  //srandom(time(NULL));
  // Construct the Liouvillian
  CoherentDissipationLiouvillian model(L,scal*g,scal*gamma);
  // Liouvillian
  const MPO& Lop=model.getLMPO();  
  // Adjoint Liouvillian
  const MPO& Ldagger=model.getLAdjointMPO();  
  // if(L<=5) Lop.exportForMatlab("myLmpo.m"); // Just debugging

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


  vector<MPS*> levels;
  for(int k=0;k<nrLevs;k++){
      char oldFile[MAXLEN];
      sprintf(oldFile,"MPS_%d_%1.2f_%1.2f_%d_%d.mat",L,g,gamma,D,k);
      const string filename=oldMPSdir+"/"+oldFile;
      MPS* mps=new MPS(L,D,d*d);
      mps->importMPS(filename.data());
      cout<<"Read level "<<k<<" from file "<<oldFile<<endl;
      levels.push_back(mps);
  }

  // Now compute and write overlap, L, and L+L  matrices
  mwArray overl,Lmat,LdL;
  
  computeMatrix(levels,Lop,Lmat);
  cout<<"Matrix for L: "<<Lmat<<endl;
  computeMatrix(levels,LdaggerL,LdL);
  cout<<"Matrix for L+L: "<<LdL<<endl;

  computeOverlapMatrix(levels,overl,d);
  cout<<"Matrix for overlap: "<<overl<<endl;

}

void computeMatrix(const vector<MPS*>& vec,const MPO& oper,mwArray& result){
  Contractor& contractor=Contractor::theContractor();  
  int dim=vec.size();
  result=mwArray(Indices(dim,dim));
  for(int k1=0;k1<dim;k1++){
    for(int k2=0;k2<dim;k2++){
      complex_t val=contractor.contract(*vec[k2],oper,*vec[k1]);
      result.setElement(val,Indices(k1,k2));
    }
  }
}

void computeOverlapMatrix(const vector<MPS*>& states,mwArray& result,int d){
  Contractor& contractor=Contractor::theContractor();  
  int dim=states.size();
  result=mwArray(Indices(dim,dim));
  int L=states[0]->getLength();
  vector<MPS*> herm;
  for(int k=0;k<dim;k++){
    const MPS* ref=states[k];
    MPS* mpsH=new MPS(*ref);
    for(int l=0;l<L;l++){
      // The Hermitian conjugate 
      mwArray aux=ref->getA(l).getA();
      int Dl=aux.getDimension(1);
      int Dr=aux.getDimension(2);
      aux.reshape(Indices(d,d,Dl,Dr));
      aux.permute(Indices(2,1,3,4),true);
      aux.reshape(Indices(d*d,Dl,Dr));
      mpsH->setA(l,aux);
    }
    herm.push_back(mpsH);
  }
  for(int k1=0;k1<dim;k1++){
    for(int k2=0;k2<dim;k2++){
      complex_t val=contractor.contract(*herm[k2],*states[k1]);
      result.setElement(val,Indices(k1,k2));
    }
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

