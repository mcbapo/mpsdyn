
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "Properties.h"

#include "HeisenbergHamiltonian.h"
#include "MPOMultiplierHermitian.h"
#include "SpinMPO.h"
#include <cmath>
#include <unistd.h>

#include "misc.h"

//#include "quicksort.h"

using namespace std;
using namespace shrt;

#define MAXLEN 120

int d=2;
double tol=1E-8;

/** Read a GS MPS from a file (needed input)
    and compute the spectrum of the effective local Hamiltonian for a certain number of sites over the whole chain.
    In the XY case, I am also computing the total Sz of the resulting "excitations"
*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  string directory="";
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

  const string mpsfile=props.getProperty("mpsfile"); // full path expected
  // Hamiltonian parameters
  double Jx_=props.getDoubleProperty("Jx");
  double Jy_=props.getDoubleProperty("Jy");
  double Jz_=props.getDoubleProperty("Jz");
  double h_=props.getDoubleProperty("h");  
  int block = props.getIntProperty("block");
  if(block<=0) block=2;
  const string outfile=props.getProperty("outfile");  
  bool app=(props.getIntProperty("append")!=0); // default: append
  int nrEV=props.getIntProperty("numberEV"); // max nr of eigenvalues to record
  const string outfileS=props.getProperty("outfileS");  
  bool computeS=!outfileS.empty();
  const string outfileSz=props.getProperty("outfileSz");  
  bool computeSz=!outfileSz.empty();

  int pos1 = props.getIntProperty("initpos");
  int pos2 = props.getIntProperty("finalpos");
  
  // Read the MPS
  if(!file_exists(mpsfile)){
    cout<<"ERROR: File for MPS "<<mpsfile<<" not found!"<<endl;
    exit(1);
  }

  MPS state(1,1,1);
  state.importMPS(mpsfile.data());
  cout<<"Read GS file "<<endl;
  int L=state.getLength();

  state.gaugeCond('R',1);
  state.gaugeCond('L',1);

  if(pos1<=0) pos1=0;
  if(pos2<=0)pos2=L-block-1;

  
  // Prepare output file
  ofstream* out;
  if(!app||!file_exists(outfile))
    out=new ofstream(outfile.data());
  else{
    out=new ofstream(outfile.data(),ios::app);
    // Write the header
  }
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output (append="<<app<<")"<<endl;
    exit(1);
  }
  out->close();delete out;

  if(computeS){
    if(!app||!file_exists(outfileS))
      out=new ofstream(outfileS.data());
    else
      out=new ofstream(outfileS.data(),ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfileS<<
	" for output (append="<<app<<")"<<endl;
      exit(1);
    }
    out->close();delete out;
  }

  vector<double> Jx(L-1,Jx_);
  vector<double> Jy(L-1,Jy_);
  vector<double> Jz(L-1,Jz_);
  vector<double> h(L,h_);

  HeisenbergHamiltonian hamH(L,Jx,Jy,Jz,h,d);
  const MPO& mpoH=hamH.getHMPO();
  MPO Szmpo(L);
  SpinMPO::getSzMPO(L,d,Szmpo);
 
  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);

  // Now it may be that my state was a purification: check
  bool purif=0; MPO doubleH(L);MPO doubleSz(L);
  int dphys=state.getA(0).getd();
  if(dphys==d*d){
    cout<<"The state read was a purification: extending the MPO of H"<<endl;
    purif=true;
    extendMPO(mpoH,doubleH,d);
    extendMPO(Szmpo,doubleSz,d);
  }

  
  complex_t initE;
  if(!purif)
    initE=contractor.contract(state,mpoH,state);
  else
    initE=contractor.contract(state,doubleH,state);
  cout<<"The energy of the read state is "<<initE<<endl;
  out=new ofstream(outfile.data(),ios::app);
  *out<<setprecision(15);
  *out<<"% E="<<real(initE)<<endl;
  *out<<endl;
  out->close();delete out;

  // max dim of the SVD, if I compute entropies
  int maxSVD=state.getA(L/2).getDr()*pow(d,(int)block/2);
    
  for(int k=pos1;k<=pos2;k++){

    MPOMultiplierHermitian Heff(block+2); // maybe could use non-Hermitian, as it is in fact Hermitian (factor 2 cheaper)
    MPOMultiplier Szeff(block+2); 
    if(!purif){
      contractor.getEffectiveOperatorMPOMultiplier(state,mpoH,block,k,Heff);
      contractor.getEffectiveOperatorMPOMultiplier(state,Szmpo,block,k,Szeff);
    }
    else{
      contractor.getEffectiveOperatorMPOMultiplier(state,doubleH,block,k,Heff);
      contractor.getEffectiveOperatorMPOMultiplier(state,doubleSz,block,k,Szeff);
    }
    
    if(purif){
      // my Heff is a product of boundaries witht eh true H and
      // identities on ancillas. The identity is tensorized with the
      // rest and I want to remove it, as it only gives degeneracy
      int firstOpPos=1;
      int lastOpPos=block;
      for(int p=firstOpPos;p<=lastOpPos;p++){
	const FoldedOperator& op=(FoldedOperator&)Heff.getOp(p);
	//cout<<"Operator at p is "<<op<<endl;
	mwArray aux=op.getDataLeft();
	Heff.setOp(p,new Operator(aux),true);
	//cout<<"Substituted op "<<p<<" by simple"<<endl;
	const FoldedOperator& opZ=(FoldedOperator&)Szeff.getOp(p);
	//cout<<"Operator at p is "<<op<<endl;
	aux=opZ.getDataLeft();
	Szeff.setOp(p,new Operator(aux),true);
      }
    }

    
    int totDim=Heff.getSize();
    cout<<"Hamiltonian for site "<<k<<" range "<<block<<" dimensions "<<totDim<<endl;
    // Now extract the exact spectrum (no eigenvectors)
    vector<complex_t> Dval;
    mwArray U;
    //wrapper::eig(Heff,Dval,U,false);
    //    wrapper::eigs(Heff,nrEV,"SR",Dval,U,false);
    //bool computeVec=computeS?true:false;
    wrapper::eigs_primme(Heff,min(nrEV,totDim),primme_smallest,Dval,U);
    out=new ofstream(outfile.data(),ios::app);
    *out<<setprecision(15);
    *out<<k<<"\t";
    for(int p=0;p<min(nrEV,(int)Dval.size());p++){
      *out<<real(Dval[p])<<"\t";
    }
    if(Dval.size()<nrEV) // fill in with 0
      for(int p=Dval.size();p<nrEV;p++) *out<<0<<"\t";
    *out<<endl;
    out->close();delete out;

    if(computeS){
      int dL=Heff.getOp(0).getdorig()*pow(d,(int)block/2);
      int dR=totDim/dL;
      out=new ofstream(outfileS.data(),ios::app);
      *out<<setprecision(15);
      for(int p=0;p<min(nrEV,(int)Dval.size());p++){
	*out<<k<<"\t"<<p<<"\t";
	// Compute the corresponding entropy
	mwArray stateP=U.subArray(Indices(-1,p)); // column p of U
	stateP.reshape(Indices(dL,dR)); 
	int nr=0;
	mwArray U,S,Vd;
	wrapper::svd(stateP,U,S,Vd);
	for(int ip=0;ip<min(dL,dR);ip++){
	  *out<<real(S.getElement(Indices(ip,ip)))<<"\t";
	}
	for(int ip=min(dL,dR);ip<maxSVD;ip++) // fill in
	  *out<<0.<<"\t";
	*out<<endl;
      }
      out->close();delete out;
    }

    if(computeSz){
      int dL=Heff.getOp(0).getdorig()*pow(d,(int)block/2);
      int dR=totDim/dL;
      //cout<<"Will now try to apply Szeff:"<<Szeff<<endl;
      out=new ofstream(outfileSz.data(),ios::app);
      *out<<setprecision(15);
      for(int p=0;p<min(nrEV,(int)Dval.size());p++){
	*out<<k<<"\t"<<p<<"\t";
	// Compute the corresponding expectation value of Sz
	mwArray stateP=U.subArray(Indices(-1,p)); // column p of U
	mwArray result;
	Szeff.product(stateP,result);
	stateP.reshape(Indices(totDim,1));
	stateP.Hconjugate();
	stateP.multiplyRight(result);
	//cout<<"Computed Sz in state "<<p<<": "<<stateP<<endl;
	complex_t valSz=stateP.getElement(0);
	*out<<real(valSz)<<"\t"<<imag(valSz)<<endl;
      }
      out->close();delete out;
    }
    
  }

}
