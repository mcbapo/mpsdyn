
#include <math.h>
#include <iomanip>
#include <sstream>

#include "misc.h"
#include "Properties.h"
//#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "SpinMPO.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "ThirringHamiltonian.h"

using namespace shrt;
using namespace std;

/** thermalThirring approximates the thermal state of
    Thirring Hamiltonian \ref <ThirringHamiltonian>,
    for a certain inverse temperature beta,
    via imaginary time evolution.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <ma> (double) adimensional mass parameter \f$ma\f$ of the Hamiltonian
    \param <g2> (double) adimensional coupling parameter \f$g^2\f$ of the Hamiltonian
    \param <gt2> (double) adimensional coupling parameter \f$\tilde{g}^2\f$
    of the Hamiltonian, to allow for a different coupling for even links
    \param <lambda> (double) penalty term \f$\lambda\f$ (\f$>0\f$)
    \param <mu> (double) chemical potential \f$\mu\f$ 
    \param <Starget> (int) targetted total \f$S_z\f$
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the energy
    \param <beta> (double) inverse temperature
    \param <delta> (double) step width for imaginary time evolution
    \param <rate> (int) frequency with with results are stored
    \param <outfname> (char*) name of the output file for the results
    \param <app> (int) whether the output file is to be kept (app==1)
*/

complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};

int main(int argc,const char* argv[]){
  // Read input arguments
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

  int L = props.getIntProperty("L");  
  if(L%2!=0){
    cout<<"ERROR: Nr of sites L should be even"<<endl;
    exit(1);
  }
  double ma=props.getDoubleProperty("ma");  
  double gt2_=props.getDoubleProperty("gtilde2");
  double g2_=props.getDoubleProperty("g2"); // the non-alternating one
  double lambda_=props.getDoubleProperty("lambda");
  if(lambda_<0) lambda_=0.;
  double mu_=props.getDoubleProperty("mu");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  int Stot_=props.getIntProperty("Starget");
  int D=props.getIntProperty("D");
  double beta=props.getDoubleProperty("beta");
  double delta=props.getDoubleProperty("delta");
  int rate=props.getIntProperty("stepBlock");
  
  string outfname=props.getProperty("output");
  //  string mpsdir=props.getProperty("mpsdir");
  //  string mpsfile=props.getProperty("mpsfile");
  //  string mpstostretch=props.getProperty("stretchmps");
  // string outfnameS=props.getProperty("outputSz"); // for entropies and polarizations
  // bool saveS=(outfnameS.length()>0);
  bool app=(props.getIntProperty("append")!=0);
  bool findGS=(props.getIntProperty("findGS")==1);

  cout<<"Initialized arguments: L="<<L
      <<", ma="<<ma
      <<", g2="<<g2_
      <<", gt2="<<gt2_
      <<", lambda="<<lambda_
      <<", mu="<<mu_
      <<", Starget="<<Stot_
      <<", outfile="<<outfname
      <<", beta="<<beta
      <<", delta="<<delta
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    *out<<"% L="<<L<<", ma="<<ma<<", g2="<<g2_<<", gt2="<<gt2_<<", Starget="<<Stot_<<", D="<<D<<endl;
  }
  else
    out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }

  *out<<setprecision(10);
  int M=beta/(2*delta);

  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  cout<<"Initialized Contractor"<<endl;
  contractor.setConvTol(tol);

  int d=2;
  MPS thS(L,D,d*d);

  MPO _projS0(L); // projector onto total spin 0
  SpinMPO::getProjectorTotalSz(L,d,_projS0); // TODO: generalize for other total spins
  // Notice:this is a projector with exponentially large Frobenius norm choose(N,N/2)
  cout<<"After creating the proj (as MPO) its norm is "<<contractor.contract(thS,thS)<<endl;
  
  // thS.setProductState(p_maxent); // the initial state, maximally entangled
  // cout<<"Initialized identity state"<<endl;
  //   // TODO: Substitute by projector onto right value of Stot!!

  MPSfromMPO(_projS0,thS);
  cout<<"After creating the proj (as MPO) its norm is "<<contractor.contract(thS,thS)<<endl;
  thS.gaugeCond('R',true);    
  cout<<"After gauge L its norm is "<<contractor.contract(thS,thS)<<endl;
  thS.gaugeCond('L',true);    
  cout<<"After gauge R its norm is "<<contractor.contract(thS,thS)<<endl;
  
  MPO Szmpo(L);
  SpinMPO::getSzMPO(L,d,Szmpo);

 
  ThirringHamiltonian hamH(L,ma,g2_,lambda_,Stot_,d,gt2_,mu_);
  //cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamH.getHMPO();

  MPS gs(L,D,d);
  double E0=0;
  if(findGS){
    contractor.findGroundState(hamil,D,&E0,gs);
    cout<<"Ground state energy="<<E0<<endl;
    if(!app) *out<<"% Ground state energy="<<E0<<endl;
  }
  else{
    if(!app) *out<<"% Ground state calculation skipped"<<endl;
  }
  if(!app){
    *out<<"%\n%delta\t cnt\t beta\t T\t <H>\t <ProjS0>\t <Z(L/2-1)Z(L/2)>\t <Z(L/2-1)>\t <Z(L/2)>\t <Sztot>";
    if(findGS) *out<<"\t <ProjGS>";
    *out<<endl;
  }
  
  int cnt=0;
  double betai=0;
  // Now try evolution with the exponential, instead
  MPO expHe(L);
  MPO expHo(L);
  MPO expHe2(L);
  MPO Hdoub(L);
  MPO projS0(L);
  MPO projS0_2(L);
  MPO Szdoub(L);
  // I construct the simple ones first, and then from them the DoubleOperators
  {
    // Regular ones
    MPO _expHe2(L),_expHe(L),_expHo(L);
    hamH.getExponentialMPOeven(_expHe2,-delta*.5*.5*ONE_c);
    hamH.getExponentialMPOeven(_expHe,-delta*.5*ONE_c);
    hamH.getExponentialMPOodd(_expHo,-delta*.5*ONE_c);

    // _expHe.exportForMatlab("expHe.m");
    // _expHo.exportForMatlab("expHo.m");
    // hamil.exportForMatlab("Htirr.m");

    
    mwArray identPhys=identityMatrix(d);
    identPhys.reshape(Indices(d,1,d,1)); 
    for(int k=0;k<L;k++){
      mwArray aux=_expHe2.getOp(k).getFullData();
      expHe2.setOp(k,new DoubleOperator(_expHe2.getOp(k).getFullData(),
					permute(aux,Indices(3,2,1,4))),true);
      aux=_expHe.getOp(k).getFullData();
      expHe.setOp(k,new DoubleOperator(_expHe.getOp(k).getFullData(),
				       permute(aux,Indices(3,2,1,4))),true);
      aux=_expHo.getOp(k).getFullData();
      expHo.setOp(k,new DoubleOperator(_expHo.getOp(k).getFullData(),
				       permute(aux,Indices(3,2,1,4))),true);
      Hdoub.setOp(k,new DoubleOperator(hamil.getOp(k).getFullData(),identPhys),true);
      projS0.setOp(k,new DoubleOperator(_projS0.getOp(k).getFullData(),identPhys),true);
      aux=_projS0.getOp(k).getFullData();
      projS0_2.setOp(k,new DoubleOperator(_projS0.getOp(k).getFullData(),
					  permute(aux,Indices(3,2,1,4))),true);
      Szdoub.setOp(k,new DoubleOperator(Szmpo.getOp(k).getFullData(),identPhys),true);
    }
    // The simple ones will be destroyed now. I already have them copied
  }
  cout<<"Created all the exponential operators"<<endl;
  // Also the Z operator
  mwArray sigZ(Indices(d*d,1),dataZ);
  mwArray idOp=identityMatrix(d);idOp.reshape(Indices(1,d*d));
  sigZ.multiplyRight(idOp);sigZ.reshape(Indices(d,d,d,d));
  sigZ.permute(Indices(1,3,2,4));sigZ.reshape(Indices(d*d,d*d));

  while(cnt<=M){
    *out<<delta<<"\t";
    *out<<cnt<<"\t";
    *out<<delta*cnt*2.<<"\t";
    *out<<1./(2.*delta*cnt)<<"\t";
    complex_t energy=contractor.contract(thS,Hdoub,thS);
    *out<<real(energy)<<"\t";
    complex_t weightS0=contractor.contract(thS,projS0,thS);
    *out<<real(weightS0)<<"\t";
    // Trick to compute Sz Sz and each of them
    MPS Z1(thS);
    Z1.applyLocalOperator(L/2-1,sigZ);
    MPS Z2(thS);
    Z2.applyLocalOperator(L/2,sigZ);
    complex_t zz=contractor.contract(Z1,Z2);
    complex_t z1=contractor.contract(Z1,thS);
    complex_t z2=contractor.contract(thS,Z2);
    *out<<real(zz)<<"\t";
    *out<<real(z1)<<"\t";
    *out<<real(z2)<<"\t";
    complex_t valSz=contractor.contract(thS,Szdoub,thS);
    *out<<real(valSz)<<"\t";
    if(findGS){
    // Also, what was the overlap with the GS?
      MPO gibbs(L);
      MPOfromMPS(thS,gibbs);
      complex_t projGS=contractor.contract2(gibbs,gs);
      *out<<abs(projGS)<<"\t";
    }
    *out<<endl;


    MPS aux(thS); // temporary copy
    cout<<"Created a tmp copy of the state with norm "<<contractor.contract(aux,aux)<<endl;
    cout<<"Orig:"<<contractor.contract(thS,thS)<<endl;
    cout<<"Overlap:"<<contractor.contract(thS,aux)<<endl;
    contractor.optimize(expHe2,aux,thS,D);
    cout<<"Applied the first exponential (e2)"<<endl;
    // Now apply r-1 times the pair Ho He
    int cntLoop=0;
    while(cntLoop<rate-1&&cnt<M-1){
      contractor.optimize(expHo,thS,aux,D);
      //cout<<"Time evolution step nr "<<cnt+1<<", time="<<time<<", value "<<M1ev<<", trace "<<trR<<endl;
      contractor.optimize(expHe,aux,thS,D);
      cnt++;cntLoop++;
    }
    contractor.optimize(expHo,thS,aux,D);
    contractor.optimize(expHe2,aux,thS,D);
    thS.gaugeCond('R',true);    
    cnt++;
    
    // from time to time, I should project back to the right sector
    if((cnt%10<rate)){
      cout<<"Applying projector at cnt="<<cnt<<" (rate="<<rate<<")"<<endl;
      MPS aux(thS); // temporary copy
      contractor.optimize(projS0_2,aux,thS,D);
      thS.gaugeCond('R',true);    
    }
    
    
   }

   out->close();
}
