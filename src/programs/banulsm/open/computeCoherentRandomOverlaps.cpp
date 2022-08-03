
//#include "mwArray.h"
#include "MPO.h"
#include "CoherentDissipationLiouvillian.h"
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

#define MAXLEN 200  // Max nr of characters in a filename


/** 
    Recompute the overlap matrix from the states computed by
    the program coherentPropsRand: 

    Receive as argument a Properties file, with all the range of 
    parameters to work with.

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
void writeHeader(ofstream& out,double Length,double g,double gamma,int D);

/** Force a site in the MPS to be Hermitian */
void makeHermitian(MPS& mps,int pos,int d);

bool file_exists(const string filename);

/** extendMPO takes a normal MPO, and transforms it into an MPO acting
    on a mixed state only on one side (the system indices). On the
    rest of the double index (whose dimension must be given in dimA,
    but is typically the same as for the physical) only the identity
    acts.*/
void extendMPO(const MPO& simpleMPO,MPO& doubleMPO,int dimA);

/** Compute purity */
void computePurity(const MPS& rho,double* traceRho2,double* traceRho);

/** Compute entropy of teh half chain */
void computeEntropy(const MPS& rho,double* traceRho2,double* traceRho);

/** Compute projector onto dark subspace */
void computeProjectorOntoSubspace(const MPS& rho,const MPO& proj,double* traceRhoProj);

/** Prepare the MPO that allows us to compute the overlap between
    reduced density matrices for n central sites */
void prepareOverlapMPO(MPO& mpo,int L,int d,int n);

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

  const int L=props.getIntProperty("L");
  const int D=props.getIntProperty("D");
  const double g=props.getDoubleProperty("g");
  const double gamma=props.getDoubleProperty("gamma");
  const string newMPSdir=props.getProperty("mpsdir");
  const string outdir=props.getProperty("outdir");
  const int rand0=props.getIntProperty("rand0");
  const int rand1=props.getIntProperty("rand1");
  const int n=props.getIntProperty("n"); // How many sites kept in the center

  char fileNameC[MAXLEN];
  Contractor& contractor=Contractor::theContractor();
  
  // Loop over the random labels from rand0 to rand1 and read MPSs
  vector<MPS> states(rand1-rand0+1,MPS(1,1,1));
  char fileName[MAXLEN];
  int nr=0;
  for(int r=rand0;r<=rand1;r++){
    sprintf(fileName,"MPS_%d_%1.2f_%1.2f_%d_r%d.mat",L,g,gamma,D,r);
    const string mpsfile=newMPSdir+"/"+fileName;
    if(file_exists(mpsfile))
      states[nr++].importMPS(mpsfile.data());
    // TODO: what if some is missing
  }

  // construct the MPO that implements partial trace and contractions
  MPO parTraceMPO(L);
  prepareOverlapMPO(parTraceMPO,L,d,n);

  // Now construct the overlap matrix, and save it 
  mwArray overlapM(Indices(nr,nr));
  for(int i=0;i<nr;i++)
    for(int j=0;j<nr;j++){
      complex_t Oij=contractor.contract(states[j],parTraceMPO,states[i]);
      overlapM.setElement(Oij,Indices(i,j));
    }

  // Write the overlap matrix to a binary file in outdir 
  sprintf(fileNameC,"overlapMatrix_%d_%1.2f_%1.2f_%d_n%d_r%d-%d.m",L,g,gamma,D,n,rand0,rand1);
  const string outOverlap=outdir+"/"+fileNameC;
  {
    ofstream ofmat(outOverlap.data());
    //overlapM.save(ofmat);
    putForMatlab(ofmat,overlapM,"overlap");
  }
  // SVD the matrix
  mwArray U,Vdagger,Dvals;
  wrapper::svd(overlapM,U,Dvals,Vdagger);

  // Also write the result of the rank calculation (singular values)
  sprintf(fileNameC,"overlapSVD_%d_%1.2f_%1.2f_%d_n%d_r%d-%d.txt",L,g,gamma,D,n,rand0,rand1);
  const string outRank=outdir+"/"+fileNameC;
  ofstream* out1=new ofstream(outRank.data());
  if(!out1->is_open()){
    cout<<"Error: impossible to open file "<<outRank<<
      " for output"<<endl;
    exit(1);
  }
  *out1<<"% singular values of the overlap matrix for n="<<n<<" and random states from "
       <<rand0<<" to "<<rand1<<endl;
  *out1<<setprecision(10);
  for(int k=0;k<nr;k++)
    *out1<<real(Dvals.getElement(Indices(k,k)))<<"\t";
  *out1<<endl;
  // Number above some tolerance (1E-6, 1E-4, 1E-2)
  *out1<<"% Nr of Singular Values above threshold: "<<endl;
  double tols[]={1E-8,1E-6,1E-4,1E-2};
  *out1<<"% ";
  for(int i=0;i<4;i++){
    double tol=tols[i];
    int cnt=0;
    for(int k=0;k<nr;k++){
      double elem=real(Dvals.getElement(Indices(k,k)));
      if(elem>=tol) cnt++;
    }
    *out1<<tol<<"\t"<<cnt<<"\t";
  }
  *out1<<endl;

  out1->close();
  delete out1;
}


void computeProjectorOntoSubspace(const MPS& rho,const MPO& proj,double* traceRhoProj){
  int L=rho.getLength();
  MPS aux(L,1,d*d);
  mwArray basic=identityMatrix(d); // basic piece of the MPO
  basic.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++){
    aux.setA(k,basic);
  }
  Contractor& contractor=Contractor::theContractor();
  *traceRhoProj=real(contractor.contract(rho,proj,aux));
}

void computePurity(const MPS& rho,double* traceRho2,double* traceRho){
  int L=rho.getLength();
  MPS aux(L,1,d*d);
  mwArray basic=identityMatrix(d); // basic piece of the MPO
  basic.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++){
    aux.setA(k,basic);
  }
  double result;
  Contractor& contractor=Contractor::theContractor();
  *traceRho=real(contractor.contract(rho,aux));
  if(abs(*traceRho)<1E-6){
    cout<<"** WARNING: Non physical state detected, trace(rho)="<<*traceRho<<endl;
    //result=-1;    
  }

  // to compute the trace of [(rho + rho^+)/2]^2, i.e. the Hermitian
  // part of rho, we need trace(rho^2). And to compute that in our
  // language, we need to construct an extra MPS, for rho^+, such that
  // in contract it gives exactly what we need.
  MPS rhoDagger(L,rho.getBond(),d*d);
  for(int k=0;k<L;k++){
    mwArray data=rho.getA(k).getA();
    int Dl=data.getDimension(1);
    int Dr=data.getDimension(2);
    data.reshape(Indices(d,d,Dl,Dr));
    data.permute(Indices(2,1,3,4),true);
    data.reshape(Indices(d*d,Dl,Dr));
    rhoDagger.replaceSite(k,data,0);
  }

  complex_t traceRho2C=contractor.contract(rho,rhoDagger);
  *traceRho2=.5*(real(traceRho2C)+real(contractor.contract(rho,rho)));

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

void writeHeader(ofstream& out,double Length,double g,double Gamma,
		 int D){
  out<<"% L="<<Length<<endl;
  out<<"% g="<<g<<endl;
  out<<"% gamma="<<Gamma<<endl;
  out<<"% D="<<D<<endl;
  Contractor& contr=Contractor::theContractor();
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


void extendMPO(const MPO& simpleMPO,MPO& doubleMPO,int dimA){
  int L=simpleMPO.getLength();
  doubleMPO.initLength(L); // just in case
  mwArray identPhys=identityMatrix(dimA);
  identPhys.reshape(Indices(dimA,1,dimA,1)); 
  for(int k=0;k<L;k++){  
    doubleMPO.setOp(k,new DoubleOperator(simpleMPO.getOp(k).getFullData(),identPhys),true);
  }
}

void prepareOverlapMPO(MPO& mpo,int L,int d,int n){
  mpo.initLength(L);
  // Two operators needed: the trace on both sides, and the true
  // identity in the central sites.
  mwArray basicId=identityMatrix(d*d);
  mwArray edgeOp(basicId);
  edgeOp.reshape(Indices(d,d,d,d));
  edgeOp.permute(Indices(1,3,2,4));
  edgeOp.reshape(Indices(d*d,d*d));
  Operator* edge=new Operator(edgeOp);
  Operator* center=new Operator(basicId);
  bool setEdge=0;int posEdge=0;
  bool setCenter=0;int posCenter=0;
  int left,right;// leftmost and rightmost sites in the n selected
  if(n%2==0)
    if(L%2==0) left=(L-n)/2;
    else left=(L-n+1)/2;
  else // n odd
    if(L%2==0) left=(L-n+1)/2;
    else left=(L-n)/2;
  right=left+n-1;
  for(int k=0;k<L;k++){
    if(k<left||k>right){
      if(!setEdge){
	mpo.setOp(k,edge,true);posEdge=k;setEdge=true;
      }
      else
	mpo.setOp(k,&mpo.getOp(posEdge),false);
    }
    else
      if(!setCenter){
	mpo.setOp(k,center,true);posCenter=k;setCenter=true;
      }
      else
	mpo.setOp(k,&mpo.getOp(posCenter),false);
  }

}


bool file_exists(const string filename){
  ifstream ifile(filename.data());
  return ifile;
}

