
#include <math.h>
#include <iomanip>

#include "mwArray.h"
#include "MPS.h"
#include "SpinMPO.h"
#include "Contractor.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define SPLIT 1 // The MPOs for the evolution are split 

#include "HeisenbergHamiltonian2D.h"

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
  int Lx=atoi(argv[++cntr]);
  int Ly=atoi(argv[++cntr]);
  int nL=atoi(argv[++cntr]);
  double Jx=atof(argv[++cntr]);
  double Jy=atof(argv[++cntr]);
  double Jz=atof(argv[++cntr]);
  double h=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  double beta=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int rate=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=1;

  cout<<"Initialized arguments: Lx="<<Lx<<" Ly="<<Ly<<" nL="<<nL
      <<", Jx="<<Jx
      <<", Jy="<<Jy
      <<", Jz="<<Jz
      <<", hz="<<h
      <<", beta="<<beta
      <<", delta="<<delta
      <<", outfile="<<outfname
      <<", D="<<D<<endl;

  ofstream* out;
  out=new ofstream(outfname,ios::app);
  *out<<"% Lx="<<Lx<<" Ly="<<Ly<<" nL="<<nL<<", J=("<<Jx<<","<<Jy<<","<<Jz<<"), h="<<h
      <<", D="<<D<<endl;
  *out<<"% delta\t nrsteps\t beta\t T\t log(Z)\t Energy\t <Sx^2>\t <Sz^2>\t <H^2>"<<endl;
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }

  *out<<setprecision(10);
  int M=beta/(2*delta);

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int dsp=2;
  int d=pow(dsp,nL);
   // First: put H and find the GS
  HeisenbergHamiltonian2D hamH(Lx,Ly,nL,Jx,Jy,Jz,0.,0.,h);
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamH.getHMPO();

  int L=hamH.getLegNr()*Lx;

  MPS thS(L,D,d*d);
  thS.setProductState(p_maxent); // the initial state, maximally entangled
  cout<<"Initialized identity state"<<endl;
  thS.gaugeCond('R',1);
  thS.gaugeCond('L',1);
  double logZ=.5*Lx*Ly*log(2);  

  MPS gs(L,D,d);gs.setRandomState();gs.gaugeCond('R',true);    
  double E0=0;
  contractor.findGroundState(hamil,D,&E0,gs);
  //hamil.exportForMatlab("Heis2D.m");
  cout<<"Ground state energy="<<E0<<endl;
  *out<<"% Ground state energy="<<E0<<endl;

  contractor.setConvTol(1E-8);
  int cnt=0;
  double betai=0;
  

  // Now try evolution with the exponential, instead
  // 1st order Trotter for simplicity
  MPO expHex(L);MPO expHe2(L);
  MPO expHox(L);
  MPO expHey(L);
  MPO expHoy(L); // only if nL>2
  MPO Hdoub(L);
  // Operators I want to measure: Sx^2+Sy^2
  MPO Sx2(L),Sy2(L);

  // If dims are large, I split the evolution ops, to be able to handle the tensors!
  MPO expHexB(L);MPO expHe2B(L);
  MPO expHoxB(L);
  MPO expHeyB(L);
  MPO expHoyB(L); // only if nL>2
  

  // hamH.getDoubleExponentialMPOevenH(expHex,-delta*.5*ONE_c);
  // hamH.getDoubleExponentialMPOoddH(expHox,-delta*.5*ONE_c);
  // if(hamH.getLegNr()>1) hamH.getDoubleExponentialMPOevenV(expHey,-delta*.5*ONE_c);
  // if(hamH.getLegNr()>2) hamH.getDoubleExponentialMPOoddV(expHoy,-delta*.5*ONE_c);
  
  // I construct the simple ones first, and then from them the DoubleOperators
  {
    MPO _expHex(L),_expHox(L),_expHey(L),_expHoy(L),_expHe2(L);
    hamH.getExponentialMPOevenH(_expHex,-delta*.5*ONE_c);
    hamH.getExponentialMPOoddH(_expHox,-delta*.5*ONE_c);
    if(hamH.getLegNr()>1){
      hamH.getExponentialMPOevenV(_expHey,-delta*.5*ONE_c);
      if(hamH.getLegNr()>2) hamH.getExponentialMPOoddV(_expHoy,-delta*.5*ONE_c);
    }
    else{
      _expHe2.initLength(L);
      hamH.getExponentialMPOevenH(_expHe2,-delta*.5*.5*ONE_c);
    }

    // For Sx^2 and Sy 2
    MPO _Sx2(L),_Sy2(L);
    {if(nL==1){
	SpinMPO::getSx2MPO(Lx*Ly,dsp,_Sx2);
	SpinMPO::getSy2MPO(Lx*Ly,dsp,_Sy2);
      }
      else{
	cout<<"TODO: Block Sx^2 and Sy^2 MPOs for multiple sites"<<endl;
	exit(1);
      }
    }

    
    mwArray identPhys=identityMatrix(d);
    identPhys.reshape(Indices(d,1,d,1)); 
    for(int k=0;k<L;k++){
      Hdoub.setOp(k,new DoubleOperator(hamil.getOp(k).getFullData(),identPhys),true);
      Sx2.setOp(k,new DoubleOperator(_Sx2.getOp(k).getFullData(),identPhys),true);
      Sy2.setOp(k,new DoubleOperator(_Sy2.getOp(k).getFullData(),identPhys),true);
      //cout<<"Hi, k="<<k<<" MPO expHex:"<<_expHex<<"; Op(k):"<<_expHex.getOp(k)<<endl;
      mwArray aux=_expHex.getOp(k).getFullData();
      if(!SPLIT)
	expHex.setOp(k,new DoubleOperator(_expHex.getOp(k).getFullData(),
					  permute(aux,Indices(3,2,1,4))),true);
      else{
	expHex.setOp(k,new DoubleOperator(_expHex.getOp(k).getFullData(),
					  identPhys),true);
	expHexB.setOp(k,new DoubleOperator(identPhys,
					   permute(aux,Indices(3,2,1,4))),true);
      }
      aux=_expHox.getOp(k).getFullData();
      if(!SPLIT)
	expHox.setOp(k,new DoubleOperator(_expHox.getOp(k).getFullData(),
					  permute(aux,Indices(3,2,1,4))),true);
      else{
	expHox.setOp(k,new DoubleOperator(_expHox.getOp(k).getFullData(),identPhys),true);
	expHoxB.setOp(k,new DoubleOperator(identPhys,permute(aux,Indices(3,2,1,4))),true);
      }
      if(hamH.getLegNr()>1){
	aux=_expHey.getOp(k).getFullData();
	if(!SPLIT)
	  expHey.setOp(k,new DoubleOperator(_expHey.getOp(k).getFullData(),
					    permute(aux,Indices(3,2,1,4))),true);
	else{
	  expHey.setOp(k,new DoubleOperator(_expHey.getOp(k).getFullData(),identPhys),true);
	  expHeyB.setOp(k,new DoubleOperator(identPhys,permute(aux,Indices(3,2,1,4))),true);
	}
	if(hamH.getLegNr()>2){
	  aux=_expHoy.getOp(k).getFullData();
	  if(!SPLIT)
	    expHoy.setOp(k,new DoubleOperator(_expHoy.getOp(k).getFullData(),
					      permute(aux,Indices(3,2,1,4))),true);
	  else{
	    expHoy.setOp(k,new DoubleOperator(_expHoy.getOp(k).getFullData(),identPhys),true);
	    expHoyB.setOp(k,new DoubleOperator(identPhys,permute(aux,Indices(3,2,1,4))),true);
	  }
	}
      }
      else{
	aux=_expHe2.getOp(k).getFullData();
	if(!SPLIT)
	  expHe2.setOp(k,new DoubleOperator(_expHe2.getOp(k).getFullData(),
					    permute(aux,Indices(3,2,1,4))),true);
	else{
	  expHe2.setOp(k,new DoubleOperator(_expHe2.getOp(k).getFullData(),identPhys),true);
	  expHe2B.setOp(k,new DoubleOperator(identPhys,permute(aux,Indices(3,2,1,4))),true);
	}
      }
      if(k==L/2){
	cout<<"Dimension of exponential operator (ex) at k="<<k<<" : "<<expHex.getOp(k).getDimensions()<<endl;
	cout<<"Dimension of exponential operator (ey) at k="<<k<<" : "<<expHey.getOp(k).getDimensions()<<endl;
      }
    }
  }
  cout<<"Created all the exponential operators"<<endl;
  // // Also the Z operator
  // mwArray sigZ(Indices(d*d,1),dataZ);
  // mwArray idOp=identityMatrix(d);idOp.reshape(Indices(1,d*d));
  // sigZ.multiplyRight(idOp);sigZ.reshape(Indices(d,d,d,d));
  // sigZ.permute(Indices(1,3,2,4));sigZ.reshape(Indices(d*d,d*d));

  while(cnt<M){
    MPS aux(thS); // temporary copy
    if(hamH.getLegNr()==1){
      // contractor.optimize(expHox,thS,aux,D);
      // contractor.optimize(expHex,aux,thS,D);
      // thS=aux;
      if(!SPLIT){
	contractor.optimize(expHe2,aux,thS,D);
	contractor.optimize(expHox,thS,aux,D);
	contractor.optimize(expHe2,aux,thS,D);
      }
      else{
	contractor.optimize(expHe2,thS,aux,D);
	contractor.optimize(expHe2B,aux,thS,D);
	contractor.optimize(expHox,thS,aux,D);
	contractor.optimize(expHoxB,aux,thS,D);
	contractor.optimize(expHe2,thS,aux,D);
	contractor.optimize(expHe2B,aux,thS,D);
      }
    }
    else{
      if(!SPLIT){
	contractor.optimize(expHex,aux,thS,D);
	contractor.optimize(expHox,thS,aux,D);
	contractor.optimize(expHey,aux,thS,D);
	if(hamH.getLegNr()>2){
	  aux=thS;
	  contractor.optimize(expHoy,aux,thS,D);
	}
      }
      else{
	contractor.optimize(expHex,thS,aux,D);
	contractor.optimize(expHexB,aux,thS,D);
	contractor.optimize(expHox,thS,aux,D);
	contractor.optimize(expHoxB,aux,thS,D);
	contractor.optimize(expHey,thS,aux,D);
	contractor.optimize(expHeyB,aux,thS,D);
	if(hamH.getLegNr()>2){
	  aux=thS;
	  contractor.optimize(expHoy,thS,aux,D);
	  contractor.optimize(expHoyB,aux,thS,D);
	}
      }
    }
    thS.gaugeCond('R',false);    
    logZ+=log(thS.getNormFact());
    thS.setNormFact(1.);
    
    *out<<delta<<"\t";
    *out<<++cnt<<"\t";
    *out<<delta*cnt*2.<<"\t";
    *out<<1./(2.*delta*cnt)<<"\t";
    *out<<logZ<<"\t";
    complex_t energy=contractor.contract(thS,Hdoub,thS);
    *out<<real(energy)<<"\t";
    complex_t sx2=contractor.contract(thS,Sx2,thS);
    *out<<real(sx2)<<"\t";
    complex_t sy2=contractor.contract(thS,Sy2,thS);
    *out<<real(sy2)<<"\t";
    complex_t energy2=contractor.contract2(Hdoub,thS);
    *out<<real(energy2)<<"\t";
    // // Trick to compute Sz Sz and each of them
    // MPS Z1(thS);
    // Z1.applyLocalOperator(L/2-1,sigZ);
    // MPS Z2(thS);
    // Z2.applyLocalOperator(L/2,sigZ);
    // complex_t zz=contractor.contract(Z1,Z2);
    // complex_t z1=contractor.contract(Z1,thS);
    // complex_t z2=contractor.contract(thS,Z2);
    // *out<<real(zz)<<"\t";
    // *out<<real(z1)<<"\t";
    // *out<<real(z2)<<"\t";
    *out<<endl;

    cout<<"Step "<<cnt<<" beta "<<delta*cnt*2.<<" E="<<energy<<" <DeltaM>="<<sx2+sy2<<endl;
   }

   out->close();
}
