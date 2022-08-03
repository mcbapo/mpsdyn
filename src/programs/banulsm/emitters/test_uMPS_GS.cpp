#include <math.h>
#include <iomanip>
#include <sstream>

#include "misc.h"
#include "uMPS.h"
//#include "Contractor.h"
//#include "DoubleOperator.h"
//#include "SpinMPO.h"
#include "Properties.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"

#include "IsingHamiltonian.h"

using namespace std;
using namespace shrt;

/** 
    Evolve an infinite Ising chain to demonstrate the usage of uMPS. In this case, imaginary time veolution for GS
*/

int d=2; // physical dimension: really silly to have it here
double tolSVD=1E-10; // default precision parameter: singular values below this are considered zero when truncating the MPS (can be redefined as argument)
double tolE=1E-5;

/** Compute the entropy of a rdm. */
double computeEntropy(const mwArray& rdm);


int main(int argc,const char* argv[]){
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

  double J_=props.getDoubleProperty("J");
  double h_=props.getDoubleProperty("h");
  double g_=props.getDoubleProperty("g");
  int nB=props.getIntProperty("nB"); // length of the unit cell in uMPS
  if(nB<0) nB=1; // default is 1

  //  int initSt=props.getIntProperty("initSt"); // initial state
  
  double delta=props.getDoubleProperty("delta"); // Trotter step
  //int M=props.getIntProperty("M");
  int Dcut=props.getIntProperty("Dcut"); // bond dimension
  string outfname=props.getProperty("outputfile");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=tolSVD;    

  bool app=0; // This program does not have the possibility to store
	      // intermediate results and be relaunched, but it can be
	      // added.


  // Get the basic piece of the evolution
  mwArray midU;
  {
    // Since the Ising Hamiltonian allows a translationally invariant
    // evolution operator, I use the trick of constructing a small
    // chain of 3 sites and use the middle one.    
    IsingHamiltonian hamH(3,d,J_,g_,h_);
    MPO _Uevol_1(3); // held just temporarily
    hamH.getUMPO(_Uevol_1,-delta,0); // the single layer normal evolution
    midU=_Uevol_1.getOp(1).getFullData();
  }

  // Now I construct a vector where all the elements are the
  // same. Notice that these are pointers to the same operator! (not
  // very good, either)
  vector<const mwArray*> Uops(nB,&midU);
  
  // and the operators to measure
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1*I_c,ZERO_c};
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c}; //Z
  mwArray sig0=identityMatrix(d);//sig0.reshape(Indices(1,d*d));
  mwArray sigX=mwArray(Indices(d,d),dataX);//sigX.reshape(Indices(1,d*d));
  mwArray sigY(Indices(d,d),dataY);//sigY.reshape(Indices(1,d*d));
  mwArray sigZ(Indices(d,d),dataZ);//sigZ.reshape(Indices(1,d*d));

  // time
  int cnt=0;
  double t=0.;

  // First of all, I initialize a uMPS and evolve it
  ofstream* out;
  if(app){
    // out=new ofstream(outfname.data(),ios::app);  *out<<setprecision(15);
  }  
  else{
    out=new ofstream(outfname.data());
    *out<<"% t\t real(tr)\t imag(tr)\t real(sigX)\t imag(sigX)\t real(sigY)\t imag(sigY)\t"
	<<" real(sigZ)\t imag(sigZ)\t real(ZZ)\t imag(ZZ)\t Dcut"<<endl;
    out->close();delete out;

  }

  int D0=1; // starting with product state
  uMPS state(nB,d,D0);
  cout<<"Created empty initial state "<<state<<endl;
  // switch(initSt){
  // case 1:{state.setProductState(p_xplus);break;}
  // case 2:{state.setProductState(p_yplus);break;}
  // case 3:{state.setProductState(p_zero);break;}
  // }
  // To test: initialize random
  srand(time(NULL));srandom(time(NULL));state.setRandomState();
  cout<<"Set random initial state "<<state<<endl;

  double E=1000;
  bool done=0;
  while(delta>=1E-4){
    {
      // Since the Ising Hamiltonian allows a translationally invariant
      // evolution operator, I use the trick of constructing a small
      // chain of 3 sites and use the middle one.    
      IsingHamiltonian hamH(3,d,J_,g_,h_);
      MPO _Uevol_1(3); // held just temporarily
      hamH.getUMPO(_Uevol_1,-delta,0); // the single layer normal evolution
      midU=_Uevol_1.getOp(1).getFullData();
    }
    done=0;
    while(!done){
      // evaluate expectation values and evolve one time step
      //cout<<"Step nr "<<cnt<<" uMPS: "<<state<<endl;
      
      // For single site exp. values, I get the rdm of one site,
      // averaged over tensors in the unit cell
      mwArray rdm(Indices(d,d));rdm.fillWithZero();
      for(int k=0;k<nB;k++){
	mwArray auxRDM;
	state.getSingleSiteRDM(k,auxRDM);      
	rdm=rdm+auxRDM;
      }
      rdm.multiplyLeft((1./nB)*ONE_c); // Avrg of nB sites
      //cout<<"Now RDM="<<rdm<<endl;
      rdm.reshape(Indices(d*d,1));
      mwArray res=reshape(sigX,Indices(1,d*d))*rdm;
      complex_t resultX=res.getElement(Indices(0,0));
      res=reshape(sigY,Indices(1,d*d))*rdm;
      complex_t resultY=res.getElement(Indices(0,0));
      res=reshape(sigZ,Indices(1,d*d))*rdm;
      complex_t resultZ=res.getElement(Indices(0,0));
      res=reshape(sig0,Indices(1,d*d))*rdm;
      complex_t trV=res.getElement(Indices(0,0));
      
      // For the ZZ value, with only one tensor I just need one calculation
      complex_t resultZZ(ZERO_c);
      vector<int> pos(2,0);pos[1]=1;
      vector<mwArray*> ops;ops.push_back(&sigZ);ops.push_back(&sigZ);
      if(nB==1){
	resultZZ=state.computeExpectationValue(2,pos,ops);
      }
      else{ // unit cell larger than 1, I need to average
	// should loop to get an average!
	for(int k=0;k<nB;k++){
	  pos[0]=k;
	  if(k<nB-1) pos[1]=k+1;
	  else pos[1]=0;
	  resultZZ+=state.computeExpectationValue(2,pos,ops);
	}
	resultZZ=(1./nB)*resultZZ;
      }
      //cout<<"resultZZ="<<resultZZ<<endl;
      // and write to screen and file
      double newE=real((J_*(resultZZ)+g_*(resultX)+h_*(resultZ))/trV);
      done=abs(1.-newE/E)<tolE;      
      if(newE>E){
	cout<<"ERROR: energy increasing from "<<E<<" to "<<newE<<endl;
      }
      E=newE;
      cout<<"t="<<t<<", <O>="<<(resultX)/trV<<", trace="
	  <<(trV)
	  <<" E="<<(J_*(resultZZ)+g_*(resultX)+h_*(resultZ))/trV
	  <<", D="<<state.getLeftVirtualDimension(0)<<endl;
      
      out=new ofstream(outfname.data(),ios::app);
      *out<<setprecision(15);
      *out<<t
	  <<"\t"<<E
	  <<"\t"<<real(trV)<<"\t"<<imag(trV)
	  <<"\t"<<real(resultX)<<"\t"<<imag(resultX)
	  <<"\t"<<real(resultY)<<"\t"<<imag(resultY)
	  <<"\t"<<real(resultZ)<<"\t"<<imag(resultZ)
	  <<"\t"<<real(resultZZ)<<"\t"<<imag(resultZZ)
	  <<"\t"<<E
	  <<"\t"<<state.getLeftVirtualDimension(0)
	  <<endl;
      out->close();delete out;
    
    
      // And now evolve (and truncate)
      state.applyMPO(Uops,tol,Dcut);
      cnt++;t+=delta;
    }
    if(done) delta=delta*.5;
  }
}
  
double computeEntropy(const mwArray& rdm){
  mwArray res=rdm;
  int d=rdm.getDimension(0);
  res=res+Hconjugate(res);  
  res=(ONE_c/res.trace())*res;// normalized rdm
  vector<complex_t> Dval;mwArray U; //placeholder:ignored
  wrapper::eig(res,Dval,U);
  double entropy=0.;
  for(int k=0;k<d;k++){
    double tmp=real(Dval[k]);
    if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
  }
  return entropy;
}


