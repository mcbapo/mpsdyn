
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
    Recompute the observables from the states computed by
    the program coherentProps: 

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

  double gamma=props.getDoubleProperty("gamma");
  const string newMPSdir=props.getProperty("mpsdir");
  const string outdir=props.getProperty("outdir");
  const string suffix=props.getProperty("suffix");
  const char* suff=suffix.data();
  int onlyPurity=props.getIntProperty("onlypurity");
  if(onlyPurity<0) onlyPurity=0;

  vector<int> allDs=props.getVectorIntProperty("allDs");
  vector<int> allLs=props.getVectorIntProperty("allLs");
  vector<double> allGs=props.getVectorDoubleProperty("allGs1");
  vector<double> allGs2=props.getVectorDoubleProperty("allGs2");
  vector<double> allGs3=props.getVectorDoubleProperty("allGs3");
  vector<double> allGs4=props.getVectorDoubleProperty("allGs4");

  for(int k=0;k<allGs2.size();k++)
    allGs.push_back(allGs2[k]);
  for(int k=0;k<allGs3.size();k++)
    allGs.push_back(allGs3[k]);
  for(int k=0;k<allGs4.size();k++)
    allGs.push_back(allGs4[k]);

  cout<<"Read from Properties: \n allLs="<<allLs<<endl;
  cout<<" allDs="<<allDs<<endl;
  cout<<" allGammas="<<allGs<<endl;

  char fileNameC[MAXLEN];
  Contractor& contractor=Contractor::theContractor();
  
  // Loop over the values
  for(int il=0;il<allLs.size();il++){
    int L=allLs[il];
    MPO DSP(L); // the dark subspace projector as a MPO on the vectorized rho
    {
      CoherentDissipationLiouvillian model(L,1.,1.); // Just for the dark state projector
      const MPO& DSP_=model.getDarkSubspaceProjectorMPO();
      extendMPO(DSP_,DSP,d);
    }
    for(int id=0;id<allDs.size();id++){
      int D=allDs[id];
      // File for purities
      sprintf(fileNameC,"L%d/EV/purity%s_%d_%1.2f_%d.txt",
	      L,suff,L,gamma,D);
      const string outPurityfile=outdir+"/"+fileNameC;
      ofstream* outPur;
      if(onlyPurity){
	outPur=new ofstream(outPurityfile.data());
	if(!outPur->is_open()){
	  cout<<"Error: impossible to open file "<<outPurityfile<<
	    " for output"<<endl;
	  exit(1);
	}
      }

      for(double g=0.00;g<=4.00;g=g+0.01){
	//      for(int ig=0;ig<allGs.size();ig++){
	//double g=allGs[ig];	
	// I construct here the name of the text output files, where I will
	// save the results
	// For the single body observables (x,y,z, site after site, each in
	// one line)
	sprintf(fileNameC,"L%d/EV/coherentSingle%s_%d_%1.2f_%1.2f_%d.txt",
		L,suff,L,g,gamma,D);
	const string out1bodyfile=outdir+"/"+fileNameC;
	// For the two body observables (xx,yy,zz, site after site, first
	// for distance l=1, until N/2, each (xx, yy, zz) for each l
	// in one line, whichstarts with an identifier and the value of l
	sprintf(fileNameC,"L%d/EV/coherentDouble%s_%d_%1.2f_%1.2f_%d.txt",
		L,suff,L,g,gamma,D);
	const string out2bodyfile=outdir+"/"+fileNameC;
	// Name of the MPS file to read
	sprintf(fileNameC,"L%d/mps/MPS%s_%d_%1.2f_%1.2f_%d.mat",L,suff,L,g,gamma,D);
	const string mpsfile=newMPSdir+"/"+fileNameC;
	
	// check if the MPs file exists
	ifstream infile(mpsfile.data());
	if(!infile.good()){
	  //cout<<"Ignoring parameters L="<<L<<", g="<<g<<", gamma="<<gamma<<", D="<<D<<" because no MPS file "<<mpsfile<<" is found"<<endl;
	  continue;
	}

	// Read the mps
	MPS myGS(L,D,d*d);
	myGS.importMPS(mpsfile.data());

	if(!onlyPurity){
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
	  writeHeader(*out1,L,g,gamma,D);
	  
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
	  writeHeader(*out2,L,g,gamma,D);
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
	}
	else{
	  cout<<"Computing only purity of the MPS for L="<<L<<", g="<<g<<", D="<<D<<endl;
	  double traceRho2,traceRho,proj;
	  computePurity(myGS,&traceRho2,&traceRho);
	  computeProjectorOntoSubspace(myGS,DSP,&proj);
	  *outPur<<g<<"\t"<<traceRho2<<"\t"<<traceRho<<"\t"<<proj<<endl;
	}
      } // end Gamma loop
      if(onlyPurity){
	outPur->close();
	delete outPur;
      }
    } // end D loop
  } // end L loop

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
