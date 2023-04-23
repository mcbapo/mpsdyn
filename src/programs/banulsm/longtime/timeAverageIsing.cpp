#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "SpinMPO.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

#include "IsingHamiltonian.h"

using namespace shrt;

/** timeAverageIsing simulates the time averaged evolution of a
    (initially pure) state under the Ising Hamiltonian \ref
    <IsingHamiltonian>.  Starting from a certain MPO (product pure
    state), it simulates the averaged evolution by repeated
    application of a map to the MPO. Then local observables are
    computed.

    Receives arguments:
    \param <L> (int) number of sites
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian 
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian 
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <initState> (int) which initial rho to consider. Options are:
                        1 - |X+>
			2 - |Y+>
			3 - |Z+>
    \param <delta> (double) width of the time step used
    \param <M> (int) max number of steps
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the results
                   (will be replaced if it exists)
    \param <Dpsi> (int) [optional] bond dimension used for the pure 
                   state (default D)
    \param <T0> (double) [optional] initial time from which to take the 
                   average (default 0). Until this time, normal evolution 
		   is applied to rho.
*/

/** Construct the initial state */
void constructInitialRho(int L,int d,int initSt,MPS& rho0);

/** Product of all identities, but sigz on site pos */
void constructLocal(int L,int pos,int d,MPS& result);

void computeLocalSigZ(const MPS& rho,int d,ofstream* out);
void computeEntropies(const MPS& rho,ofstream* out);

/** Vectorize an MPO to make an MPS out of it. If up==true (default)
    the first component of the double physical index will be the one
    which was going up.*/
void MPSfromMPO(const MPO& mpo,MPS& mps,bool up=true);

/** Set the step dependent coefficients of the MPO */
void setStepMPOCoefficients(int step,int L,mwArray& C,mwArray& Cl,mwArray& Cr);
/** construct the step dependent MPO */
void assembleMPO(int step,MPO& mapMPO,IsingHamiltonian& hamH,const MPS& rho0,double delta);

int d=2;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double J_=atof(argv[++cntr]);
  double g_=atof(argv[++cntr]);
  double h_=atof(argv[++cntr]);
  int initSt=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int Dpsi=D;
  if(argc>10)
    Dpsi=atoi(argv[++cntr]);
  double T0=0;
  if(argc>11)
    T0=atof(argv[++cntr]);

  int M0=ceil(T0/delta);

  cout<<"Initialized arguments "
      <<", L="<<L
      <<", J="<<J_
      <<", g="<<g_
      <<", h="<<h_
      <<", initSt="<<initSt
      <<", delta="<<delta
      <<", M="<<M
      <<", D="<<D
      <<", outfile="<<outfname
      <<endl;

  ofstream* out;
  out=new ofstream(outfname);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<"% L="<<L<<", J="<<J_<<", g="<<g_<<", h="<<h_
      <<" D="<<D<<endl;

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;

  // First: create the Hamiltonian
  *out<<setprecision(15);
  *out<<"% M\t delta\t t\t D\t <M1>\t tr(rho)\t tr(rho H)"<<endl;

  IsingHamiltonian hamH(L,d,J_,g_,h_);
  hamH.registerTimeDependence();
  cout<<"Created the Hamiltonian"<<endl;
  //hamH.getHMPO().exportForMatlab("isingH.m");
  // Now I construct by hand the MPO. But since it will be time
  // dependent, here I construct only the basic pieces, and then, step
  // by step, I need to replace the proper coefficients and reassemble
  // the MPO.

  // Initial state rho0
  MPS rho0(L,1,d*d);
  constructInitialRho(L,d,initSt,rho0);
  MPS rhoInit(rho0); // need a copy to reinsert it!
  // Place for the MPO
  MPO mapMPO(L);

  // Also an initial pure state, to check
  MPS psi(L,1,d);
  switch(initSt){
  case 1:{psi.setProductState(p_xplus);break;}
  case 2:{psi.setProductState(p_yplus);break;}
  case 3:{psi.setProductState(p_zero);break;}
  }
  // And the evolution operator
  MPO Uevol(L);
  hamH.getUMPO(Uevol,-delta*I_c,0);
  

  // Also the identity, to compute the trace
  MPS Id(L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++)
    Id.setA(k,id2);
  
  MPS hmps(L,1,d*d);
  const MPO& hamil=hamH.getHMPO();
  MPSfromMPO(hamil,hmps);

  // Sum of sigmaz/sigmax operators
  MPS SxSum(L,2,d*d),SzSum(L,2,d*d);
  {
    MPO Saux(L);
    SpinMPO::getSxMPO(L,d,Saux);
    MPSfromMPO(Saux,SxSum);
    SpinMPO::getSzMPO(L,d,Saux);
    MPSfromMPO(Saux,SzSum);
  }
  // And I will also need the sigmaz/sigmax operator on central site, which I could compute ony by one, but I assemble it here
  MPS SxCenter(Id);
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  sigX.reshape(Indices(d*d,1,1));
  SxCenter.setA(L/2,sigX);sigX.reshape(Indices(d,d));
  
  // Now do the evolution step by step, and compute the expectation value of H, Sum(sigmax) Sum(sigmaz) and sigmax in center every time
  int cnt=0;
  double time=0.;
  complex_t trR=contractor.contract(rho0,Id);
  complex_t Eev=contractor.contract(rho0,hmps);
  complex_t Sxev=contractor.contract(rho0,SxSum);
  complex_t Szev=contractor.contract(rho0,SzSum);
  complex_t SxCev=contractor.contract(rho0,SxCenter);
  MPS auxP(psi);auxP.applyLocalOperator(L/2,sigX,false);
  complex_t vectorSx=contractor.contract(auxP,psi);

  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D
      <<"\t"<<real(trR)<<"\t"<<imag(trR)
      <<"\t"<<real(Eev)<<"\t"<<imag(Eev)
      <<"\t"<<real(Sxev)<<"\t"<<imag(Sxev)
      <<"\t"<<real(Szev)<<"\t"<<imag(Szev)
      <<"\t"<<real(SxCev)<<"\t"<<imag(SxCev)
      <<"\t"<<real(vectorSx)<<"\t"<<imag(vectorSx)
      <<endl;

  MPO expH(L);
  { // Just for the initial evolution
    MPO _expHA(L); // for the second indices (first is Uevol)
    hamH.getUMPO(_expHA,delta*I_c,0.); // not imaginary time
    for(int k=0;k<L;k++){
      expH.setOp(k,new DoubleOperator(Uevol.getOp(k).getFullData(),permute(_expHA.getOp(k).getFullData(),Indices(3,2,1,4))),true);
    }
  }
  for(int step=1;step<=M0;step++){
    MPS aux(rho0);
    contractor.optimize(expH,aux,rho0,D);
    trR=contractor.contract(rho0,Id);
    double absTr=abs(trR);
    complex_t phase=(1/absTr)*trR;
    double factor=pow(1./absTr,double(1./L));
    rho0.setA(0,factor*conjugate(phase)*rho0.getA(0).getA());
    for(int k=1;k<L;k++){
      rho0.setA(k,factor*rho0.getA(k).getA());
    }
    auxP=psi;
    contractor.optimize(Uevol,auxP,psi,Dpsi);psi.gaugeCond('R',1);
    //cout<<"Time evolution step nr "<<cnt+1<<", time="<<time<<", value "<<M1ev<<", trace "<<trR<<endl;
    time=step*delta;
    trR=contractor.contract(rho0,Id);
    Eev=contractor.contract(rho0,hmps);
    Sxev=contractor.contract(rho0,SxSum);
    Szev=contractor.contract(rho0,SzSum);
    SxCev=contractor.contract(rho0,SxCenter);
    auxP=psi;auxP.applyLocalOperator(L/2,sigX,false);
    vectorSx=contractor.contract(auxP,psi);

    *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D
	<<"\t"<<real(trR)<<"\t"<<imag(trR)
	<<"\t"<<real(Eev)<<"\t"<<imag(Eev)
	<<"\t"<<real(Sxev)<<"\t"<<imag(Sxev)
	<<"\t"<<real(Szev)<<"\t"<<imag(Szev)
	<<"\t"<<real(SxCev)<<"\t"<<imag(SxCev)
	<<"\t"<<real(vectorSx)<<"\t"<<imag(vectorSx)
	<<endl;
    cout<<"Time evolution step nr "<<step<<", time="<<time<<", value "<<SxCev/trR<<", trace "<<trR<<endl;
  }
  rhoInit=rho0;
  for(int step=M0+1;step<M;step++){
    assembleMPO(step-M0,mapMPO,hamH,rhoInit,delta);
    // if(step==5)
    //   mapMPO.exportForMatlab("myMapMPO.m");
    MPS aux(rho0);
    contractor.optimize(mapMPO,aux,rho0,D);
    // After each step, I need to renormalize:
    trR=contractor.contract(rho0,Id);
    double absTr=abs(trR);
    complex_t phase=(1/absTr)*trR;
    double factor=pow(1./absTr,double(1./L));
    rho0.setA(0,factor*conjugate(phase)*rho0.getA(0).getA());
    for(int k=1;k<L;k++){
      rho0.setA(k,factor*rho0.getA(k).getA());
    }

    auxP=psi;
    contractor.optimize(Uevol,auxP,psi,Dpsi);psi.gaugeCond('R',1);
    //cout<<"Time evolution step nr "<<cnt+1<<", time="<<time<<", value "<<M1ev<<", trace "<<trR<<endl;
    time=step*delta;

    trR=contractor.contract(rho0,Id);
    Eev=contractor.contract(rho0,hmps);
    Sxev=contractor.contract(rho0,SxSum);
    Szev=contractor.contract(rho0,SzSum);
    SxCev=contractor.contract(rho0,SxCenter);
    auxP=psi;auxP.applyLocalOperator(L/2,sigX,false);
    vectorSx=contractor.contract(auxP,psi);

    *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D
	<<"\t"<<real(trR)<<"\t"<<imag(trR)
	<<"\t"<<real(Eev)<<"\t"<<imag(Eev)
	<<"\t"<<real(Sxev)<<"\t"<<imag(Sxev)
	<<"\t"<<real(Szev)<<"\t"<<imag(Szev)
	<<"\t"<<real(SxCev)<<"\t"<<imag(SxCev)
	<<"\t"<<real(vectorSx)<<"\t"<<imag(vectorSx)
	<<endl;
    cout<<"Time evolution step nr "<<step<<", time="<<time<<", value "<<SxCev/trR<<", trace "<<trR<<endl;

  }

  out->close();
  delete out;
}


void constructLocal(int L,int pos,int d,MPS& result){
   mwArray sig0=identityMatrix(d);
   sig0.reshape(Indices(d*d,1,1));
   complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
   mwArray sigZ(Indices(d*d,1,1),dataZ);
   result=MPS(L,1,d*d);
   for(int k=0;k<L;k++){
     result.setA(k,sig0);
   }
   result.setA(pos,sigZ);
}


void computeLocalSigZ(const MPS& rho,int d,ofstream* out){
  // cout<<"Computing local expectation values "<<endl;
   mwArray sig0=identityMatrix(d);
   sig0.reshape(Indices(d*d,1,1));
   complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
   mwArray sigZ(Indices(d*d,1,1),dataZ);
   int L=rho.getLength();
   MPS oper(L,1,d*d);
   for(int k=0;k<L;k++){
     oper.setA(k,sig0);
   }

   Contractor& contractor=Contractor::theContractor();
   for(int l=0;l<L;l++){
     oper.setA(l,sigZ);
     complex_t resL=contractor.contract(rho,oper);
     *out<<real(resL)<<"\t"<<imag(resL)<<"\t";
     oper.setA(l,sig0);
   }


}

void computeEntropies(const MPS& rho,ofstream* out){
  //cout<<"Computing entropies "<<endl;
   int D=rho.getBond();
   if(D>40) 
     cout<<"Skipping entropy calculation due to too large D"<<endl;
   int L=rho.getLength();
   int LmaxL,LminR,cnt;
   if(L%2==0){ //even length
     LmaxL=L/2-1;LminR=L/2;
     if(L/2%2==0) cnt=L/4;
     else cnt=(L/2-1)/2;
   }
   else{ // odd length
     LmaxL=(L-1)/2;LminR=LmaxL;
     if((L-1)/2%2==0) cnt=(L-1)/4+1;
     else cnt=(L-3)/4+1;
   }
   Contractor& contractor=Contractor::theContractor();
   int pos1=LmaxL;int pos2=LminR;
   for(int l=0;l<cnt;l++){
     //cout<<"\t length "<<pos2-pos1+1<<endl;
     double entropy=0.;
     if(D<=40)
       entropy=contractor.getEntropy(rho,pos1,pos2);
     *out<<pos2-pos1+1<<"\t"<<entropy<<"\t";
     pos1--;pos2++;
   }
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


void setStepMPOCoefficients(int step,int L,mwArray& C,mwArray& Cl,mwArray& Cr){
  double coeff1=pow(1./step,1./double(L));
  double coeff2=pow(double(step-1.)/double(step),1./double(L));
  // coeff2=1.;coeff1=0.;
  cout<<"Setting coefficients for step "<<step-1<<" to "<<coeff1<<" and "<<coeff2<<endl;
  C.reshape(Indices(2,2,2));Cl.reshape(Indices(1,2,2));Cr.reshape(Indices(2,1,2));
  C.setElement(coeff1*ONE_c,Indices(0,0,0));
  C.setElement(coeff2*ONE_c,Indices(1,1,1));
  Cl.setElement(coeff1*ONE_c,Indices(0,0,0));
  Cl.setElement(coeff2*ONE_c,Indices(0,1,1));
  Cr.setElement(coeff1*ONE_c,Indices(0,0,0));
  Cr.setElement(coeff2*ONE_c,Indices(1,0,1));
  C.reshape(Indices(2*2,2));Cl.reshape(Indices(1*2,2));Cr.reshape(Indices(2*1,2));
}

void constructInitialRho(int L,int d,int initSt,MPS& rho0){
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sig0=identityMatrix(d);
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray sigX(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  
  mwArray state;
  switch(initSt){
  case 1: // X+
    {
      state=.5*(sig0+sigX);
      break;
    }
  case 2: // Y+
    {
      state=.5*(sig0+sigY);
      break;
    }
  case 3: // Z+
    {
      state=.5*(sig0+sigZ);
      break;
    }
  default:
    cout<<"Error: unknown type onf intSt="<<initSt<<endl;
    exit(1);
  }
  state.reshape(Indices(d*d,1,1));

  rho0=MPS(L,1,d*d);
  for(int k=0;k<L;k++)
    rho0.setA(k,state);
}

void assembleMPO(int step,MPO& mapMPO,IsingHamiltonian& hamH,const MPS& rho0,double delta){
  static bool initialized(false);
  // Basic operators I need to use (here I cannot use the structure inside, but 
  // have to build them as block mwArray
  // I will keep the two blocks that are then combined with different weights depending on the step
  int D=4;
  static mwArray Z(Indices(d*d,D+1,d*d,D+1));
  static mwArray Z0(Indices(d*d,D+1,d*d,D+1));
  static mwArray Zl(Indices(d*d,1,d*d,D+1));
  static mwArray Zl0(Indices(d*d,1,d*d,D+1));
  static mwArray Zr(Indices(d*d,D+1,d*d,1));
  static mwArray Zr0(Indices(d*d,D+1,d*d,1));
  int L=hamH.getHMPO().getLength();
  if(!initialized){
    mwArray UUd,UUdl,UUdr,rho0Phi;
    // I construct the simple ones first, and then from them the DoubleOperators
    // Regular ones
    MPO _expH(L);
    //MPO _expHA(L); // for the second indices
    hamH.getUMPO(_expH,-delta*I_c,0.); // not imaginary time
    // _expH.exportForMatlab("SoloU.m");
    //hamH.getUMPO(_expHA,delta*I_c,0.);
    // Instead of the full MPO for UUd, I just get the tensors at the edges and center
    //DoubleOperator _UUdl(_expH.getOp(0).getFullData(),permute(_expHA.getOp(0).getFullData(),Indices(3,2,1,4)));
    DoubleOperator _UUdl(_expH.getOp(0).getFullData(),conjugate(_expH.getOp(0).getFullData()));
    UUdl=_UUdl.getFullData();
    //putForMatlab(cout,UUdl,"UUdl");
    if(L>2){
      //DoubleOperator _UUd(_expH.getOp(1).getFullData(),permute(_expHA.getOp(1).getFullData(),Indices(3,2,1,4)));
      DoubleOperator _UUd(_expH.getOp(1).getFullData(),conjugate(_expH.getOp(1).getFullData()));
      UUd=_UUd.getFullData();
    }
    //DoubleOperator _UUdr(_expH.getOp(L-1).getFullData(),permute(_expHA.getOp(L-1).getFullData(),Indices(3,2,1,4)));
    DoubleOperator _UUdr(_expH.getOp(L-1).getFullData(),conjugate(_expH.getOp(L-1).getFullData()));
    UUdr=_UUdr.getFullData();
    //putForMatlab(cout,UUdr,"UUdr");
    // And the basic piece for the rho0 contribution
    rho0Phi=identityMatrix(d);rho0Phi.reshape(Indices(1,d*d));
    mwArray aux=rho0.getA(0).getA();aux.reshape(Indices(d*d,1));
    rho0Phi.multiplyLeft(aux);rho0Phi.reshape(Indices(d*d,d*d));
    // Now place them in the corresponding blocks of Z tensors
    for(int l=0;l<D;l++)
      for(int r=0;r<D;r++)
	for(int i=0;i<d*d;i++)
	  for(int j=0;j<d*d;j++){
	    Zl.setElement(UUdl.getElement(Indices(i,0,j,r)),Indices(i,0,j,r+1));
	    if(L>2) Z.setElement(UUd.getElement(Indices(i,l,j,r)),Indices(i,l+1,j,r+1));
	    Zr.setElement(UUdr.getElement(Indices(i,l,j,0)),Indices(i,l+1,j,0));
	  }
    // And the same for the term with rho0
    for(int i=0;i<d*d;i++)
      for(int j=0;j<d*d;j++){
	Zl0.setElement(rho0Phi.getElement(Indices(i,j)),Indices(i,0,j,0));
	if(L>2) Z0.setElement(rho0Phi.getElement(Indices(i,j)),Indices(i,0,j,0));
	Zr0.setElement(rho0Phi.getElement(Indices(i,j)),Indices(i,0,j,0));
      }
    cout<<"Created all the basic pieces of the operators"<<endl;
    // putForMatlab(cout,Zl,"Zl");
    // putForMatlab(cout,Zl0,"Zl0");
    // putForMatlab(cout,Zr,"Zr");
    // putForMatlab(cout,Zr0,"Zr0");
    initialized=true;
  }
  // Now compute the coefficients that multiply each block depending on the step:
  double coeff1=pow(1./(step+1.),1./double(L));
  double coeff2=pow(double(step)/double(step+1.),1./double(L));
  // Now assemble and link operators
  // Left one
  Operator* OpL=new Operator(coeff1*Zl0+coeff2*Zl);
  Operator* OpR=new Operator(coeff1*Zr0+coeff2*Zr);
  mapMPO.setOp(0,OpL,true);
  mapMPO.setOp(L-1,OpR,true);
  if(L>2){
    Operator* Op=new Operator(coeff1*Z0+coeff2*Z);
    mapMPO.setOp(1,Op,true);
    for(int k=2;k<L-1;k++) mapMPO.setOp(k,Op,false);
  }
}


// void assembleMPO(int step,MPO& mapMPO,IsingHamiltonian& hamH,const MPS& rho0,double delta){
//   static bool initialized(false);
//   // Basic operators I need to use (here I cannot use the structure inside, but 
//   // have to build them as block mwArray
//   // I will keep the two blocks that are then combined with different weights depending on the step
//   static mwArray Z(Indices(2,d*d,2,d*d,2));
//   static mwArray Zl(Indices(2,d*d,1,d*d,2));
//   static mwArray Zr(Indices(2,d*d,2,d*d,1));
//   // Coefficients of the MPO (middle, left and right edges)
//   static mwArray C(Indices(2,2,2));
//   static mwArray Cl(Indices(1,2,2));
//   static mwArray Cr(Indices(2,1,2));
//   int L=hamH.getHMPO().getLength();
//   if(!initialized){
//     mwArray UUd,UUdl,UUdr,rho0Phi;
//     // I construct the simple ones first, and then from them the DoubleOperators
//     // Regular ones
//     MPO _expH(L);
//     //MPO _expHA(L); // for the second indices
//     hamH.getUMPO(_expH,-delta*I_c,0.); // not imaginary time
//     _expH.exportForMatlab("SoloU.m");
//     //hamH.getUMPO(_expHA,delta*I_c,0.);
//     // Instead of the full MPO for UUd, I just get the tensors at the edges and center
//     //DoubleOperator _UUdl(_expH.getOp(0).getFullData(),permute(_expHA.getOp(0).getFullData(),Indices(3,2,1,4)));
//     DoubleOperator _UUdl(_expH.getOp(0).getFullData(),conjugate(_expH.getOp(0).getFullData()));
//     UUdl=_UUdl.getFullData();
//     //putForMatlab(cout,UUdl,"UUdl");
//     if(L>2){
//       //DoubleOperator _UUd(_expH.getOp(1).getFullData(),permute(_expHA.getOp(1).getFullData(),Indices(3,2,1,4)));
//       DoubleOperator _UUd(_expH.getOp(1).getFullData(),conjugate(_expH.getOp(1).getFullData()));
//       UUd=_UUd.getFullData();
//     }
//     //DoubleOperator _UUdr(_expH.getOp(L-1).getFullData(),permute(_expHA.getOp(L-1).getFullData(),Indices(3,2,1,4)));
//     DoubleOperator _UUdr(_expH.getOp(L-1).getFullData(),conjugate(_expH.getOp(L-1).getFullData()));
//     UUdr=_UUdr.getFullData();
//     //putForMatlab(cout,UUdr,"UUdr");
//     // And the basic piece for the rho0 contribution
//     rho0Phi=identityMatrix(d);rho0Phi.reshape(Indices(1,d*d));
//     mwArray aux=rho0.getA(0).getA();aux.reshape(Indices(d*d,1));
//     rho0Phi.multiplyLeft(aux);rho0Phi.reshape(Indices(d*d,d*d));
//     for(int j1=0;j1<d*d;j1++)
//       for(int j3=0;j3<d*d;j3++){
// 	//cout<<"Setting in Z element ("<<j1<<","<<j3<<") of rho0Phi="<<rho0Phi.getElement(Indices(j1,0,j3,0))<<endl;
// 	if(L>2)
// 	  Z.setElement(rho0Phi.getElement(Indices(j1,j3)),Indices(0,j1,0,j3,0));
// 	Zl.setElement(rho0Phi.getElement(Indices(j1,j3)),Indices(0,j1,0,j3,0));
// 	Zr.setElement(rho0Phi.getElement(Indices(j1,j3)),Indices(0,j1,0,j3,0));
// 	for(int j2=0;j2<2;j2++){
// 	  Zr.setElement(UUdr.getElement(Indices(j1,j2,j3,0)),Indices(1,j1,j2,j3,0));
// 	  //cout<<"Setting in Zr element ("<<j1<<","<<j3<<"), j2="<<j2<<" of UdUr="<<UUdr.getElement(Indices(j1,j2,j3,0))<<endl;
// 	  for(int j4=0;j4<2;j4++){
// 	    if(L>2)
// 	      Z.setElement(UUd.getElement(Indices(j1,j2,j3,j4)),Indices(1,j1,j2,j3,j4));
// 	    Zl.setElement(UUdl.getElement(Indices(j1,0,j3,j4)),Indices(1,j1,0,j3,j4));
// 	    //cout<<"Setting in Zl element ("<<j1<<","<<j3<<"), j4="<<j4<<" of UdUl="<<UUdl.getElement(Indices(j1,0,j3,j4))<<endl;
// 	  }
// 	}
//       }    
//     Z.reshape(Indices(2,d*d*2*d*d*2));Zl.reshape(Indices(2,d*d*1*d*d*2));Zr.reshape(Indices(2,d*d*2*d*d*1));
//     cout<<"Created all the basic pieces of the operators"<<endl;
//     initialized=true;
//   }
//   // For this step, set the coefficients of Cs
//   setStepMPOCoefficients(step+1,L,C,Cl,Cr);
//   // Now assemble and link operators
//   // Left one
//   mwArray OpL=Cl*Zl;OpL.reshape(Indices(1,2,d*d,1,d*d,2));OpL.permute(Indices(3,1,4,5,2,6));OpL.reshape(Indices(d*d,1*1,d*d,2*2));
//   mwArray OpR=Cr*Zr;OpR.reshape(Indices(2,1,d*d,2,d*d,1));OpR.permute(Indices(3,1,4,5,2,6));OpR.reshape(Indices(d*d,2*2,d*d,1*1));
//   mapMPO.setOp(0,new Operator(OpL),true);
//   mapMPO.setOp(L-1,new Operator(OpR),true);
//   if(L>2){
//     mwArray Op=C*Z;Op.reshape(Indices(2,2,d*d,2,d*d,2));Op.permute(Indices(3,1,4,5,2,6));Op.reshape(Indices(d*d,2*2,d*d,2*2));
//     Operator* Op_=new Operator(Op);
//     mapMPO.setOp(1,Op_,true);
//     for(int k=2;k<L-1;k++) mapMPO.setOp(k,Op_,false);
//   }
// }
