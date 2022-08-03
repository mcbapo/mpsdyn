
#include <math.h>
#include <iomanip>

#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "HeisenbergHamiltonian.h"

using namespace shrt;

/** thermalHeisenberg approximates the thermal state of
    Heisenberg Hamiltonian \ref <HeisenbergHamiltonian>,
    for a certain inverse temperature beta,
    via imaginary time evolution.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <Jx> (double) parameter \f$Jx\f$ of the Hamiltonian
    \param <Jy> (double) parameter \f$Jy\f$ of the Hamiltonian
    \param <Jz> (double) parameter \f$Jz\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <beta> (double) inverse temperature
    \param <delta> (double) step width for imaginary time evolution
    \param <rate> (int) frequency with with results are stored
    \param <outfname> (char*) name of the output file for the results
*/

complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double Jx_=atof(argv[++cntr]);
  double Jy_=atof(argv[++cntr]);
  double Jz_=atof(argv[++cntr]);
  double h_=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);

  double beta=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int rate=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=1;

  cout<<"Initialized arguments: L="<<L
      <<", Jx="<<Jx_
      <<", Jy="<<Jy_
      <<", Jz="<<Jz_
      <<", h="<<h_
      <<", beta="<<beta
      <<", delta="<<delta
      <<", outfile="<<outfname
      <<", D="<<D<<endl;

  ofstream* out;
  out=new ofstream(outfname,ios::app);
  *out<<"% L="<<L<<", J=("<<Jx_<<","<<Jy_<<","<<Jz_<<"), h="<<h_
      <<", D="<<D<<endl;
  *out<<"% delta\t nrsteps\t beta\t T\t Energy\t <Sz(N/2-1)Sz(N/2)>"
      <<"\t <Sz(N/2-1)>\t<Sz(N/2)>"<<endl;
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }

  *out<<setprecision(10);
  int M=beta/(2*delta);

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d=2;
  MPS thS(L,D,d*d);
  thS.setProductState(p_maxent); // the initial state, maximally entangled
  cout<<"Initialized identity state"<<endl;
   
   // First: put H and find the GS
   // Set coefficients
  vector<double> Jx(L-1,Jx_);
  vector<double> Jy(L-1,Jy_);
  vector<double> Jz(L-1,Jz_);
  vector<double> h(L,h_);
  HeisenbergHamiltonian hamH(L,Jx,Jy,Jz,h,d);
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamH.getHMPO();

  MPS gs(L,D,d);
  double E0=0;
  contractor.findGroundState(hamil,D,&E0,gs);
  cout<<"Ground state energy="<<E0<<endl;
  *out<<"% Ground state energy="<<E0<<endl;

  contractor.setConvTol(1E-8);
  int cnt=0;
  double betai=0;
  // Now try evolution with the exponential, instead
  MPO expHe(L);
  MPO expHo(L);
  MPO expHe2(L);
  MPO Hdoub(L);
  // I construct the simple ones first, and then from them the DoubleOperators
  {
    // Regular ones
    MPO _expHe2(L),_expHe(L),_expHo(L);
    hamH.getExponentialMPOeven(_expHe2,-delta*.5*.5*ONE_c);
    hamH.getExponentialMPOeven(_expHe,-delta*.5*ONE_c);
    hamH.getExponentialMPOodd(_expHo,-delta*.5*ONE_c);
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
    }
    // The simple ones will be destroyed now. I already have them copied
  }
  cout<<"Created all the exponential operators"<<endl;
  // Also the Z operator
  mwArray sigZ(Indices(d*d,1),dataZ);
  mwArray idOp=identityMatrix(d);idOp.reshape(Indices(1,d*d));
  sigZ.multiplyRight(idOp);sigZ.reshape(Indices(d,d,d,d));
  sigZ.permute(Indices(1,3,2,4));sigZ.reshape(Indices(d*d,d*d));

  while(cnt<M){
    MPS aux(thS); // temporary copy
    contractor.optimize(expHe2,aux,thS,D);
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

    *out<<delta<<"\t";
    *out<<++cnt<<"\t";
    *out<<delta*cnt*2.<<"\t";
    *out<<1./(2.*delta*cnt)<<"\t";
    complex_t energy=contractor.contract(thS,Hdoub,thS);
    *out<<real(energy)<<"\t";
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
    *out<<endl;
   }

   out->close();
}
