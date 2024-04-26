
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"

#include "IsingHamiltonian.h"
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
    For Ising model and Heisenberg model (depending on input parameters).
    This one should be general, so that block of size 0 is supported (and default).
    Computing also entropies, total Sz, total S2 of the excitations (all optional: if output files specified).
    It supports enforcing a total magnetization by a penalty term (if Starget and penH are provided as parameters)
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
  double J=props.getDoubleProperty("J");
  double g=props.getDoubleProperty("g");
  double h=props.getDoubleProperty("h");
  bool ising=(props.getIntProperty("ising")!=0); // only if -ising=0 is given it will be Heisenberg (or XY)  
  // These are only used if Heisenberg model (h is the same)
  double Jx_ = props.getDoubleProperty("Jx");  
  double Jy_ = props.getDoubleProperty("Jy");  
  double Jz_ = props.getDoubleProperty("Jz");
  // It can also include x and y local fields Heisenberg model
  double hx_ = props.getDoubleProperty("hx");  
  double hy_ = props.getDoubleProperty("hy");  

  int block = props.getIntProperty("block");
  if(block<0) block=0;
  const string outfile=props.getProperty("outfile");  
  bool app=(props.getIntProperty("append")!=0); // default: append
  int nrEV=props.getIntProperty("numberEV"); // max nr of eigenvalues to record
  const string outfileS=props.getProperty("outfileS");  
  bool computeS=!outfileS.empty();
  const string outfileSz=props.getProperty("outfileSz");  
  bool computeSz=!outfileSz.empty();
  const string outfileS2=props.getProperty("outfileS2");  
  bool computeS2=!outfileS2.empty();
  const string outfileH2=props.getProperty("outfileVariance");  
  bool computeH2=!outfileH2.empty(); // whether to compute the variance of the MPSs for the excitations  

  int pos1 = props.getIntProperty("initpos");
  int pos2 = props.getIntProperty("finalpos");
  int Dmps = props.getIntProperty("D");

  // optional penalty term to enforce total Sz (only in Heisenberg case)
  double penH=props.getDoubleProperty("penaltySz");  
  if(penH<0) penH=0;
  int Starget=props.getIntProperty("Starget");
  
  // Read the MPS
  if(!file_exists(mpsfile)){
    cout<<"ERROR: File for MPS "<<mpsfile<<" not found!"<<endl;
    exit(1);
  }

  MPS state(1,1,1);
  state.importMPS(mpsfile.data());
  cout<<"Read GS file "<<endl;
  int L=state.getLength();

  //  state.gaugeCond('R',1);
  state.gaugeCond('L',1); // all sites gauge to L

  MPS stateR(state);
  stateR.gaugeCond('R',1); // all sites gauge to R in the copy
  
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

  // First: put H and find the GS
  Hamiltonian* ham0;
  if(ising){
    ham0=new IsingHamiltonian(L,d,J,g,h);
  }
  else{
    ham0=new HeisenbergHamiltonian(L,Jx_,Jy_,Jz_,hx_,hy_,h,d,0.,Starget,penH);
  }

  const MPO& mpoH=ham0->getHMPO();
  MPO Szmpo(L);
  SpinMPO::getSzMPO(L,d,Szmpo);
  MPO S2mpo(L);
  SpinMPO::getS2MPO(L,d,S2mpo);

  
  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);

  complex_t initE;
  initE=contractor.contract(state,mpoH,state);
  cout<<"The energy of the read state is "<<initE<<endl;
  out=new ofstream(outfile.data(),ios::app);
  *out<<setprecision(15);
  *out<<"% E="<<real(initE)<<endl;
  *out<<endl;
  out->close();delete out;

  // max dim of the SVD, if I compute entropies
  int maxSVD=state.getA(L/2).getDr()*pow(d,(int)block/2);
  
  // Use a copy of the state to include the gauge
  
  for(int k=pos1;k<=pos2;k++){
    // could be optimized by copying only the next one
    for(int ik=0;ik<=k;ik++){ // copy the left part of the MPS (gauge to R) onto the one with gauge to L
      state.replaceSite(ik,stateR.getA(ik).getA()); // dimensions should match
    } // notice state has changed, for gauge is not consistent
      // (central piece not included), but it does not matter for
      // Heff, which only needs the orthogonal basis
    bool gauged=true; // since I already have the proper gauge, Contractor does not need to do it for each operator
    MPOMultiplierHermitian Heff(block+2); // maybe could use non-Hermitian, as it is in fact Hermitian (factor 2 cheaper)
    MPOMultiplier Szeff(block+2); 
    MPOMultiplier S2eff(block+2); 
    MPOMultiplier H2eff(block+2); // for variance 
    contractor.getEffectiveOperatorMPOMultiplier(state,mpoH,block,k,Heff,gauged);
    contractor.getEffectiveOperatorMPOMultiplier(state,Szmpo,block,k,Szeff,gauged);
    contractor.getEffectiveOperatorMPOMultiplier(state,S2mpo,block,k,S2eff,gauged);
    {
      const MPO* ptrs[]={&mpoH,&mpoH};
      MPO ham2(L);
      MPO::join(2,ptrs,ham2);// MPO for H^2

      contractor.getEffectiveOperatorMPOMultiplier(state,ham2,block,k,H2eff,gauged);

    }
    
    // If block is 0-sized, the MPS is not the GS, as I didn't keep
    // the central tensor, so I skip initializing with the current
    // state and compute it instead    
    MPS gs(block+2,Dmps,d);
    // get the lowest
    double lambda=L;
    contractor.findGroundState(Heff,Dmps,&lambda,gs);
    complex_t initE=contractor.contract(gs,Heff,gs);
    complex_t initH2=contractor.contract(gs,H2eff,gs);
    complex_t valSz=contractor.contract(gs,Szeff,gs);
    complex_t valS2=contractor.contract(gs,S2eff,gs);
    
    vector<MPS*> levels; levels.push_back(&gs);
    vector<double> energies; energies.push_back(real(initE));
    vector<double> variances; variances.push_back(real(initH2-initE*conjugate(initE)));
    vector<double> valsSz; valsSz.push_back(real(valSz));
    vector<double> valsS2; valsS2.push_back(real(valS2));

    // // Check the energy
    // cout<<"Energy of the GS computed with Heff("<<k<<") "<<real(computedE0)<<endl;
    // if(abs(computedE0-initE)>1E-5*abs(initE)){
    //   cout<<"Apparently there is an error!!"<<endl;
    //   exit(1);
    // }
    
    int lval=0; // last one computed
    
    while(lval<nrEV-1){
      double lambda=L;
      // compute one more
      MPS* exc=new MPS(gs);
      exc->setRandomState();
      exc->gaugeCond('R',1);
      exc->gaugeCond('L',1);

      //cout<<"In newly created MPS energy is "<<contractor.contract(*exc,Heff,*exc)<<endl;
      // orthogonalize it to gs
      //contractor.orthogonalize(*exc,Dmps,levels,vector<double>(levels.size(),1.));
      contractor.findNextExcitedState(Heff,Dmps,levels,&lambda,*exc,-L);

      complex_t initH2=contractor.contract(gs,H2eff,gs);
      complex_t valSz=contractor.contract(*exc,Szeff,*exc);	
      complex_t valS2=contractor.contract(*exc,S2eff,*exc);	
      
      levels.push_back(exc);
      energies.push_back(lambda);
      variances.push_back(real(initH2)-lambda*lambda);
      valsSz.push_back(real(valSz));
      valsS2.push_back(real(valS2));
      
      lval++;
      
      if(computeS){
	vector<complex_t> sv_;
	contractor.getSchmidtValues(*exc,sv_);	
	out=new ofstream(outfileS.data(),ios::app);
	*out<<setprecision(15);
	*out<<k<<"\t"<<lval<<"\t";
	for(int ip=0;ip<sv_.size();ip++){
	  *out<<real(sv_[ip])<<"\t";
	}
	for(int ip=sv_.size();ip<maxSVD;ip++) // fill in
	  *out<<0.<<"\t";
	*out<<endl;
	out->close();delete out;
      }
    }
    
    // Now write all output to files
    out=new ofstream(outfile.data(),ios::app);
    *out<<setprecision(15);
    *out<<k<<"\t";
    for(int p=0;p<min(nrEV,(int)energies.size());p++){
	*out<<energies[p]<<"\t";
    }
    if(energies.size()<nrEV) // fill in with 0
      for(int p=energies.size();p<nrEV;p++) *out<<0<<"\t";
    *out<<endl;
    out->close();delete out;

    if(computeH2){
      out=new ofstream(outfileH2.data(),ios::app);
      *out<<setprecision(15);
      *out<<k<<"\t";
      for(int p=0;p<min(nrEV,(int)energies.size());p++){
	*out<<variances[p]<<"\t";
      }
      if(variances.size()<nrEV) // fill in with 0
	for(int p=variances.size();p<nrEV;p++) *out<<0<<"\t";
      *out<<endl;
      out->close();delete out;
    }

    if(computeSz){
      //cout<<"Will now try to apply Szeff:"<<Szeff<<endl;
      out=new ofstream(outfileSz.data(),ios::app);
      *out<<setprecision(15);
      *out<<k<<"\t";
      for(int p=0;p<min(nrEV,(int)valsSz.size());p++){
	*out<<real(valsSz[p])<<"\t";
      }
      if(valsSz.size()<nrEV) // fill in with 0
	for(int p=valsSz.size();p<nrEV;p++) *out<<0<<"\t";
      *out<<endl;
      out->close();delete out;
    }

    if(computeS2){
      out=new ofstream(outfileS2.data(),ios::app);
      *out<<setprecision(15);
      //      *out<<k<<"\t"<<lval<<"\t"<<real(valSz)<<"\t"<<imag(valSz)<<endl;
      *out<<k<<"\t";
      for(int p=0;p<min(nrEV,(int)valsS2.size());p++){
	*out<<real(valsS2[p])<<"\t";
      }
      if(valsS2.size()<nrEV) // fill in with 0
	for(int p=valsS2.size();p<nrEV;p++) *out<<0<<"\t";
      *out<<endl;
      out->close();delete out;
    }

  }
}
