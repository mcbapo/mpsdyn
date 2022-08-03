#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "SpinMPO.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

#define PROJECTEDZ0 0

#include "HeisenbergAncillaryNoiseHamiltonian.h"

using namespace shrt;

/** localizationTI simulates the time evolution of an MPO under the
    special Hamiltonian \ref <HeisenbergAncillaryNoiseHamiltonian>, in
    which the effect of the random magnetic field is achieved by a
    second spin chain, interacting with the first via an Ising type
    term. Starting from a certain MPO, in whih the second chain is
    NOT maximally mixed 
    BUT in a pure state with equal prob of any z1...zL eigenstate
    and the first one has a non-homogeneous
    polarization, it simulates the evolution and computes the average
    and localized magnetization (of the first chain) along the way.
    To keep track of how good the evolution is, I will also compute
    the energy.

    Receives arguments:
    \param <L> (int) number of sites (in the original chain)
    \param <J> (double) parameter \f$J\f$ of the Heisenberg
                        Hamiltonian for the ancillary chain
    \param <B> (double) parameter \f$B\f$ of the interaction
                        between both chains
    \param <delta> (double) width of the time step used
    \param <M> (int) number of steps
    \param <rate> (int) frequency of recording data (nr of steps)
                       Results will be saved for 0,rate,2*rate...,M
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the results
                   (will be replaced if it exists)
    \param <mpsfile> (char*) name of the file where to save the MPS 
    \param <savingFreq> (int) [optional] frequency (steps) to save 
                    intermediate MPS files. If 0 or none, no intermediate 
		    states will be saved
    \param <newInstance> (bool) whether to start a new evolution from
                         the inhomogeneous polarization distribution
			 (if not, the MPS in mpsfile is used)
    \param <M0> (int) OPTIONAL, only read if newInstance is false, it
                         is the number of steps already done in the
			 existing file.
*/

/** Construct the initial state and observable, but vectorized*/
void constructM1(int L,int d,MPS& result);
/** Idem for the one with the exponential complex coefficients */
void constructM1exp(int L,int d,MPS& result);

void constructM1single(int L,int d,MPS& result,bool cosi=false);
/** Construct the initial state for the stm (M1) plus ancilla (X+) */
void constructInitState(int L,int d,MPS& result);

/** Construct a special MPO to compute the average variance of 
\f$\langle S_z^{(n)}\rangle\f$ over the first spin chain, the 
second or both (depending on whether the fourth argument is 0 (default), 1 or 2). */
void constructVarSzMPO(int L,int d,MPO& varSz,int chain=0);

/** Product of all identities, but sigz on site pos */
void constructLocal(int L,int pos,int d,MPS& result);

/** Construct the projector onto total spin 0 for one of the subchains
    (if ancilla==true, the ancillas, otherwise the first chain). 
    If both is true, then the projector is the product of both Sz
    projectors (and the value of ancilla is ignored). Each of them
    has bond dimension L/2, and on both, it is the square!! */
void constructProjZ0(int L,int d,MPO& result,bool ancilla=false,bool both=false);

/** Construct, as initial state, the EXACT projection of M1 onto the
    Sz=0 subspace, tensor the same projector for the ancillas. */
void constructM1P0(int L,int d, MPS& result);

/** construct also the observable (as a MPS which I can contract with
    rho) for total Sz in the primary or ancillary system. */
void constructSzMPS(int L,int d,MPS& result,bool ancilla=false);

void computeLocalSigZ(const MPS& rho,int d,ofstream* out);
void computeEntropies(const MPS& rho,ofstream* out);
void computeErrors(const MPS& rho,ofstream* out,const vector<int>& Ds);

/** Vectorize an MPO to make an MPS out of it. If up==true (default)
    the first component of the double physical index will be the one
    which was going up.*/
void MPSfromMPO(const MPO& mpo,MPS& mps,bool up=true);
/** doubleMPO taxes a normal MPO, acting on (system)
    physical dimensions, and constracts a double version, as in the
    evolution of a thermal state, for instance, where a second
    (transposed) copy of the MPO acts on the (ancillary) double indices. */
void doubleMPO(const MPO& simpleMPO,MPO& doubleMPO);
/** extendMPO takes a normal MPO, and transforms it into an MPO acting
    on a mixed state only on one side (the system indices). On the
    rest of the double index (whose dimension must be given in dimA,
    but is typically the same as for the physical) only the identity
    acts.*/
void extendMPO(const MPO& simpleMPO,MPO& doubleMPO,int dimA);
/** Analogous to the one above, but puts the transposed MPO onto the ancillary indices*/
void extendTransposeMPO(const MPO& simpleMPO,MPO& doubleMPO,int dimA);

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double J_=atof(argv[++cntr]);
  double B_=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int rate=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  const char* mpsfile=argv[++cntr];
  int savingFreq=0;
  if(argc>10)
    savingFreq=atoi(argv[++cntr]);
  int newInstance=atoi(argv[++cntr]);
  if(newInstance<0) newInstance=1;
  int M0=0;
  if(!newInstance&&argc>12) M0=atoi(argv[++cntr]);


  cout<<"Initialized arguments "
      <<", L="<<L
      <<", J="<<J_
      <<", B="<<B_
      <<", delta="<<delta
      <<", M="<<M
      <<", D="<<D
      <<", outfile="<<outfname
      <<", mpsfile="<<mpsfile
      <<", newInstance="<<newInstance
      <<", M0="<<M0
      <<endl;

  ofstream* out;
  bool appendingFile=(!newInstance&&file_exists(outfname));
  if(appendingFile)
    out=new ofstream(outfname,ios::app);
  else // created again
    out=new ofstream(outfname);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  if(!appendingFile){
    *out<<"% L="<<L<<", J="<<J_<<", B="<<B_
	<<" D="<<D<<endl;
    *out<<"% MPS in "<<mpsfile<<endl;
    *out<<"% M\t delta\t t\t D\t <M1>\t tr(rho)\t tr(rho H)"<<endl;
  }
  *out<<setprecision(15);
  int Mtot=M0+M;

  Contractor& contractor=Contractor::theContractor();
  //contractor.setConvTol(1E-10);
  cout<<"Initialized Contractor"<<endl;
  int d=2;

  // //  if(1){
  //   // srandom(time(NULL));
  //   srandom(12);
  //   // Set random coefficients
  //   vector<double> Bs(L,B_);
  //   for(int k=0;k<L;k++){
  //     Bs[k]=(1.*random()/RAND_MAX)*B_;
  //   }
  //   cout<<"Bs="<<Bs<<endl;
  //   //}

  // First: create the Hamiltonian
  HeisenbergAncillaryNoiseHamiltonian hamH(L,J_,B_);
  //HeisenbergAncillaryNoiseHamiltonian hamH(L,J_,Bs);
  //hamH.getHMPO().exportForMatlab("hamilHeisAnc.m");  
  cout<<"Created the Hamiltonian"<<endl;
  if(0&&L<10){
    // I could now compute the ground state to test
    MPS testGS(2*L,D,d);double lambda=0;
    contractor.setEigenSolver(primme);
    contractor.findGroundState(hamH.getHMPO(),D,&lambda,testGS);
    cout<<setprecision(15);
    cout<<"Found GS with energy "<<lambda<<endl;
    //    exit(1);
  }

  // The evolution operator for a rho will result from double operators
  // Which MPOs do I need for the indices in both sides
  MPO expHSe(2*L),expHSo(2*L),expHSo2(2*L),
    expHTe(2*L),expHTo(2*L),expHTo2(2*L),
    expHB(2*L),expHB2(2*L);
  hamH.getDoubleExponentialMPOeven(expHSe,-delta*I_c,true);
  hamH.getDoubleExponentialMPOodd(expHSo,-delta*I_c,true);
  hamH.getDoubleExponentialMPOodd(expHSo2,-delta*.5*I_c,true);
  hamH.getDoubleExponentialMPOeven(expHTe,-delta*I_c,false);
  hamH.getDoubleExponentialMPOodd(expHTo,-delta*I_c,false);
  hamH.getDoubleExponentialMPOodd(expHTo2,-delta*.5*I_c,false);
  hamH.getDoubleExponentialMPOmixed(expHB,-delta*I_c);
  hamH.getDoubleExponentialMPOmixed(expHB2,-delta*.5*I_c);
  //   cout<<"Created all the exponential operators:"<<endl;
  // cout<<"expHeSigma:"<<expHSe<<endl;
  // cout<<"expHeTau:"<<expHTe<<endl;  
  // cout<<"expHoSigma:"<<expHSo<<endl;  
  // cout<<"expHoTau:"<<expHTo<<endl;  
  // cout<<"expHoSigma(1/2):"<<expHSo2<<endl;  
  // cout<<"expHoTau(1/2):"<<expHTo2<<endl;  
  // cout<<"expHB:"<<expHB<<endl; 
  // { 
  // // Just to test placement of terms!
  // MPO _expHSe(L),_expHSo(L),_expHB(L);
  // hamH.getExponentialMPOeven(_expHSe,-delta*I_c,false);
  // cout<<"simple expHeSigma:"<<_expHSe<<endl;
  // hamH.getExponentialMPOodd(_expHSo,-delta*I_c,false);
  // hamH.getExponentialMPOmixed(_expHB,-delta*I_c);
  // _expHSe.exportForMatlab("expHeSigma.m");  
  // _expHSo.exportForMatlab("expHoSigma.m");  
  // _expHB.exportForMatlab("expHB.m");  
  // }

  // Now construct the initial rho MPO (term 1x...xsigmazx...x1)
  MPS rho0(2*L,2,d*d);
  MPS M1(2*L,2,d*d);
  //  constructInitState(L,d,rho0);
  constructM1exp(L,d,M1);
  if(!newInstance)
    rho0.importMPS(mpsfile);
  else
    rho0=M1;
  //constructM1(L,d,M1);
  //cout<<"Created the intial state M1 "<<rho0<<endl;

  // I would also like to compute the expectation values of Sz (total) and Sz^2
  MPS Sz(2*L,2,d*d);
  MPS Sz2(2*L,2,d*d);
  {
    MPO aux(2*L);
    SpinMPO::getSzMPO(2*L,d,aux);
    MPSfromMPO(aux,Sz);
    SpinMPO::getSz2MPO(2*L,d,aux);
    MPSfromMPO(aux,Sz2);
  }
  //cout<<"Created the total Sz "<<Sz<<endl;
  // And also the total Sz on ancilla and primary sites
  MPS SzSig(2*L,2,d*d);
  MPS SzTau(2*L,2,d*d);
  constructSzMPS(L,d,SzSig,false);
  constructSzMPS(L,d,SzTau,true);
  //  cout<<"Created the total Sz on sigmas (taus) "<<SzSig<<endl;

  // A special one to compute the average variance of the individual 
  // Sz values (in sigma sites!). this has a trick, as I will compute 
  // the expectation value of the MPO in vector rho, to obtain the sum 
  // of squared <S_z^{(n)}> for all n's
  MPO varSz(2*L);
  constructVarSzMPO(L,d,varSz,0);


  // I will now construct an initial state which is mixture of all
  // Sz=0 states
  //  if(0)
  // {
  //   MPO aux(2*L),aux2(2*L);
  //   constructProjZ0(L,d,aux,false,true);
  //   doubleMPO(aux,aux2);
  //   MPS tmp(rho0);
  //   contractor.optimize(aux2,tmp,rho0,2*(L/2+1)*(L/2+1)*(L/2+1));
    //    MPSfromMPO(aux,rho0);
  // }
  /** De hecho, valdria con proyectar solo las sigmas, por izqda y
      drcha, si el M1, en vez de con Id en taus, lo hiciera con el proyector*/
  MPS M1proj(2*L,2,d*d);
  constructM1P0(L,d,M1proj);

  /** Instead of starting with the projected state, I can just project
      fter evolution, actually, before (or at the same time as) each
      of the measurements. For this, I need the projector (prod of
      projectors) onto Sz=0 for each subchain. I can apply them in
      sequence, so I can just construct a (double) MPO for each and
      then join them. */
  MPO projZ0stm(2*L),projZ0anc(2*L);
  {
    MPO auxZ0(2*L);bool ancilla=false;
    constructProjZ0(L,d,auxZ0,ancilla); // for the stm
    doubleMPO(auxZ0,projZ0stm);
    ancilla=true;
    constructProjZ0(L,d,auxZ0,ancilla);
    doubleMPO(auxZ0,projZ0anc);
  }
  MPO joinedProj(2*L);
  MPO joinedProjVar(2*L);
  if(PROJECTEDZ0){
//       const MPO* projs[]={&projZ0stm,&projZ0anc,&varSz};
//       MPO::join(3,projs,joinedProjVar);

    const MPO* projs[]={&projZ0stm,&projZ0anc,&varSz};
    MPO::join(2,projs,joinedProj);
    MPO::join(3,projs,joinedProjVar);
  }

  cout<<"Constructed the variance MPO "<<varSz<<endl;
  // Also the identity, to compute the trace
  MPS Id(2*L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<2*L;k++)
    Id.setA(k,id2);
  
  MPS hmps(L,1,d*d);
  const MPO& hamil=hamH.getHMPO();
  MPSfromMPO(hamil,hmps);

  // And I will also need the sigmaz operator on every site, which I can compute ony by one  

  //  M1.exportForMatlab("M1.m");
  //Id.exportForMatlab("Id.m");

  double normF=pow(2,2*L); // trace of the identity rho
  
  // Now do the evolution step by step, and compute the expectation value of M1 every time
  int cnt=M0;
  double time=cnt*delta;
  complex_t M1ev,trR,Eev,SZev,SZ2ev,vari,SZSigev,SZTauev;
  if(PROJECTEDZ0){
    cout<<"Time evolution step nr "<<cnt<<", time="<<time<<", contracting with projectors in between"<<endl;
    //    M1ev=contractor.contract(rho0,joinedProj,M1); // takes conj(M1) and
    //                                 // thus effectively M'
    M1ev=contractor.contract(rho0,M1proj); // takes conj(M1) and
    // trR=contractor.contract(rho0,joinedProj,Id);
//     Eev=contractor.contract(rho0,joinedProj,hmps);
//     SZev=contractor.contract(rho0,joinedProj,Sz);
//     SZ2ev=contractor.contract(rho0,joinedProj,Sz2);
// //     MPO joinedProjVar(2*L);
// //     const MPO* projs[]={&projZ0stm,&projZ0anc,&varSz};
// //     MPO::join(3,projs,joinedProjVar);
//     vari=contractor.contract(rho0,joinedProjVar,rho0);
//     SZSigev=contractor.contract(rho0,joinedProj,SzSig);
    // SZTauev=contractor.contract(rho0,joinedProj,SzTau);
  }
  else{ // normal evolution
    M1ev=contractor.contract(rho0,M1); // takes conj(M1) and
                                       // thus effectively M'
    trR=contractor.contract(rho0,Id);
    Eev=contractor.contract(rho0,hmps);
    SZev=contractor.contract(rho0,Sz);
    SZ2ev=contractor.contract(rho0,Sz2);
    vari=contractor.contract(rho0,varSz,rho0);
    SZSigev=contractor.contract(rho0,SzSig);
    SZTauev=contractor.contract(rho0,SzTau);
  }
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<real(M1ev)<<"\t"<<imag(M1ev)
      <<"\t"<<real(trR)<<"\t"<<imag(trR)
      <<"\t"<<real(Eev)<<"\t"<<imag(Eev)<<"\t"
      <<"\t"<<real(SZev)<<"\t"<<imag(SZev)<<"\t"
      <<"\t"<<real(SZ2ev)<<"\t"<<imag(SZ2ev)<<"\t"
      <<"\t"<<real(vari)<<"\t"<<imag(vari)<<"\t"
      <<"\t"<<real(SZSigev)<<"\t"<<imag(SZSigev)<<"\t"
      <<"\t"<<real(SZTauev)<<"\t"<<imag(SZTauev)<<"\t";
  // vector<complex_t> lambda2;
  // contractor.getSchmidtValues(rho0,lambda2); // lambda^2 in the middle
  // for(int id=0;id<D;id++)
  //   *out<<real(lambda2[id])<<"\t";
  // //  computeEntropies(rho0,out);
  // computeLocalSigZ(rho0,d,out);
  // //  computeErrors(rho0,out,Ds);
  *out<<endl;
  out->close(); // I probably don't want to keep it open for days!
  delete out;
  while(cnt<Mtot){
    MPS aux(rho0); // temporary copy
    contractor.optimize(expHSo2,rho0,aux,D);
    if(J_!=0)
      contractor.optimize(expHTo2,aux,rho0,D);
    else
      rho0=aux;
    // Now apply r-1 times the pair Ho He
    int cntLoop=0;
    while(cntLoop<rate&&cnt<Mtot){
      aux=rho0;  // AQUI!!!!*****
      contractor.optimize(expHB2,aux,rho0,D);
      contractor.optimize(expHSe,rho0,aux,D);
      if(J_!=0)
	contractor.optimize(expHTe,aux,rho0,D);
      else
	rho0=aux;
      contractor.optimize(expHB2,rho0,aux,D);
      if(cntLoop<rate-1){
	// contractor.optimize(expHSo,rho0,aux,D);
	// contractor.optimize(expHTo,aux,rho0,D);
	contractor.optimize(expHSo,aux,rho0,D);
	if(J_!=0)
	  contractor.optimize(expHTo,rho0,aux,D);
	else
	  aux=rho0;
      }
      else{
	// contractor.optimize(expHSo2,rho0,aux,D);
	// contractor.optimize(expHTo2,aux,rho0,D);
	contractor.optimize(expHSo2,aux,rho0,D);
	if(J_!=0)
	  contractor.optimize(expHTo2,rho0,aux,D);
	else
	  aux=rho0;
      }
      rho0=aux;
      //cout<<"Time evolution step nr "<<cnt+1<<", time="<<time<<", value "<<M1ev<<", trace "<<trR<<endl;
      cnt++;cntLoop++;time+=delta;
    }

    if(PROJECTEDZ0){
      cout<<"Time evolution step nr "<<cnt+1<<", time="<<time<<", contracting with projectors in between"<<endl;
      M1ev=contractor.contract(rho0,M1proj); // takes conj(M1) and
//       M1ev=contractor.contract(rho0,joinedProj,M1); // takes conj(M1) and
//       // thus effectively M'
//       trR=contractor.contract(rho0,joinedProj,Id);
//       Eev=contractor.contract(rho0,joinedProj,hmps);
//       SZev=contractor.contract(rho0,joinedProj,Sz);
//       SZ2ev=contractor.contract(rho0,joinedProj,Sz2);
// //       MPO joinedProjVar(2*L);
// //       const MPO* projs[]={&projZ0stm,&projZ0anc,&varSz};
// //       MPO::join(3,projs,joinedProjVar);
//       vari=contractor.contract(rho0,joinedProjVar,rho0);
//       SZSigev=contractor.contract(rho0,joinedProj,SzSig);
//       SZTauev=contractor.contract(rho0,joinedProj,SzTau);
    }
    else{ // normal evolution
      M1ev=contractor.contract(rho0,M1);
      trR=contractor.contract(rho0,Id);
      Eev=contractor.contract(rho0,hmps);
      SZev=contractor.contract(rho0,Sz);
      SZ2ev=contractor.contract(rho0,Sz2);
      vari=contractor.contract(rho0,varSz,rho0);
      SZSigev=contractor.contract(rho0,SzSig);
      SZTauev=contractor.contract(rho0,SzTau);
    }
    cout<<"Time evolution step nr "<<cnt<<", time="<<time<<", value "<<M1ev<<", trace "<<trR<<endl;
    out=new ofstream(outfname,ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<real(M1ev)<<"\t"<<imag(M1ev)
	<<"\t"<<real(trR)<<"\t"<<imag(trR)
	<<"\t"<<real(Eev)<<"\t"<<imag(Eev)<<"\t"
	<<"\t"<<real(SZev)<<"\t"<<imag(SZev)<<"\t"
	<<"\t"<<real(SZ2ev)<<"\t"<<imag(SZ2ev)<<"\t"
	<<"\t"<<real(vari)<<"\t"<<imag(vari)<<"\t"
	<<"\t"<<real(SZSigev)<<"\t"<<imag(SZSigev)<<"\t"
	<<"\t"<<real(SZTauev)<<"\t"<<imag(SZTauev)<<"\t";

    // contractor.getSchmidtValues(rho0,lambda2); // lambda^2 in the middle
    // for(int id=0;id<D;id++)
    //   *out<<real(lambda2[id])<<"\t";
    // //    computeEntropies(rho0,out);
    // computeLocalSigZ(rho0,d,out);
    // //    computeErrors(rho0,out,Ds);
    *out<<endl;
    out->close();
    delete out;
    // export the MPS if cnt is multiple of 100
    if(savingFreq&&cnt%savingFreq==0){
      char name[120];
      cout<<"Saving step "<<cnt<<" to file"<<endl;
      sprintf(name,"%s_%d",mpsfile,cnt);
      rho0.exportMPS(name);
    }
  }

  rho0.exportMPS(mpsfile);
}

void constructM1(int L,int d,MPS& result){
   mwArray sig0=identityMatrix(d);
   sig0.reshape(Indices(d*d,1,1));
   complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
   mwArray sigZ(Indices(d*d,1,1),dataZ);
   mwArray Z(Indices(2,d*d));
   for(int j=0;j<d*d;j++){
     Z.setElement(sig0.getElement(Indices(j,0,0)),Indices(0,j));
     Z.setElement(sigZ.getElement(Indices(j,0,0)),Indices(1,j));
   }
   int D=2; // Dimension of the MPS
   // In this case, the even sites (in [0,2*(L-1)]) are like the
   // original ones, and on the even ones, I have the identity on the
   // tau spins, and the identity on the virtual bond, except for the
   // last one (2*L-1) which finishes
   // Original terms (for even sites)
   mwArray C(Indices(D,D,2)),C1(Indices(1,D,2));
   C.setElement(ONE_c,Indices(0,0,0));
   C1.setElement(ONE_c,Indices(0,0,0));
   C.setElement(ONE_c,Indices(D-1,D-1,0));
   // There is a site-dependent element C(0,D-1,1)
   C1.setElement(ONE_c,Indices(0,D-1,1));
   C1.reshape(Indices(1*D,2));
   C1.multiplyRight(Z);
   C1.reshape(Indices(1,D,d*d));
   C1.permute(Indices(Indices(3,1,2)));
   result.setA(0,C1);
   for(int k=1;k<L;k++){
     //C.setElement(exp((-2*M_PIl*k/L)*I_c),Indices(0,D-1,1));
     C.setElement(cos((2*M_PIl*k/L))*ONE_c,Indices(0,D-1,1));
     //C.setElement(cos((M_PIl*k/L))*ONE_c,Indices(0,D-1,1));
     if(k==L-1) // the last one HAS TO put the operator
       C.setElement(ZERO_c,Indices(0,0,0));
     mwArray C_(C);
     C_.reshape(Indices(D*D,2));
     C_.multiplyRight(Z);
     C_.reshape(Indices(D,D,d*d));
     C_.permute(Indices(Indices(3,1,2)));
     result.setA(2*k,C_);
   }
   // Now for the odd ones it is just the product of identities,
   // except for the last site
   mwArray idPhys=identityMatrix(d);
   mwArray idVirt=identityMatrix(D);idVirt.reshape(Indices(1,D*D));
   idVirt.multiplyLeft(reshape(idPhys,Indices(d*d,1)));
   idVirt.reshape(Indices(d*d,D,D));
   for(int k=0;k<L-1;k++)
     result.setA(2*k+1,idVirt);
   // The last one does not pass on D
   complex_t data[]={ONE_c,ONE_c};
   idVirt=mwArray(Indices(1,D),data);
   idVirt.multiplyLeft(reshape(idPhys,Indices(d*d,1)));
   idVirt.reshape(Indices(d*d,D,1));
   result.setA(2*L-1,idVirt);
}

void constructM1exp(int L,int d,MPS& result){
   mwArray sig0=identityMatrix(d);
   sig0.reshape(Indices(d*d,1,1));
   complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
   mwArray sigZ(Indices(d*d,1,1),dataZ);
   mwArray Z(Indices(2,d*d));
   for(int j=0;j<d*d;j++){
     Z.setElement(sig0.getElement(Indices(j,0,0)),Indices(0,j));
     Z.setElement(sigZ.getElement(Indices(j,0,0)),Indices(1,j));
   }
   int D=2; // Dimension of the MPS
   // In this case, the even sites (in [0,2*(L-1)]) are like the
   // original ones, and on the even ones, I have the identity on the
   // tau spins, and the identity on the virtual bond, except for the
   // last one (2*L-1) which finishes
   // Original terms (for even sites)
   mwArray C(Indices(D,D,2)),C1(Indices(1,D,2));
   C.setElement(ONE_c,Indices(0,0,0));
   C1.setElement(ONE_c,Indices(0,0,0));
   C.setElement(ONE_c,Indices(D-1,D-1,0));
   // There is a site-dependent element C(0,D-1,1)
   C1.setElement(ONE_c,Indices(0,D-1,1)); // Coeff is 1 for n=0
   C1.reshape(Indices(1*D,2));
   C1.multiplyRight(Z);
   C1.reshape(Indices(1,D,d*d));
   C1.permute(Indices(Indices(3,1,2)));
   result.setA(0,C1);
   for(int k=1;k<L;k++){
     C.setElement(exp((-2*M_PIl*k/L)*I_c),Indices(0,D-1,1));
     if(k==L-1) // the last one HAS TO put the operator
       C.setElement(ZERO_c,Indices(0,0,0));
     mwArray C_(C);
     C_.reshape(Indices(D*D,2));
     C_.multiplyRight(Z);
     C_.reshape(Indices(D,D,d*d));
     C_.permute(Indices(Indices(3,1,2)));
     result.setA(2*k,C_);
   }
   // Now for the odd ones it is just the product of identities,
   // except for the last site
   mwArray idPhys=identityMatrix(d);
   mwArray idVirt=identityMatrix(D);idVirt.reshape(Indices(1,D*D));
   idVirt.multiplyLeft(reshape(idPhys,Indices(d*d,1)));
   idVirt.reshape(Indices(d*d,D,D));
   for(int k=0;k<L-1;k++)
     result.setA(2*k+1,idVirt);
   // The last one does not pass on D
   complex_t data[]={ONE_c,ONE_c};
   idVirt=mwArray(Indices(1,D),data);
   idVirt.multiplyLeft(reshape(idPhys,Indices(d*d,1)));
   idVirt.reshape(Indices(d*d,D,1));
   result.setA(2*L-1,idVirt);
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

void computeErrors(const MPS& rho,ofstream* out,const vector<int>& Ds){
  int nr=Ds.size();
  int L=rho.getLength();
  int D=rho.getBond();
  Contractor& contractor=Contractor::theContractor();
  for(int k=0;k<nr;k++){
    int Dref=Ds[k];
    double error=0.;
    if(Dref<D){
    // Find best approximation with bond Dref to the given state
      MPS aux(rho,Dref);
      contractor.optimizeMPS(rho,aux,Dref,&error);
    }
    //    else Dref>=D // no error
    *out<<Dref<<"\t"<<error<<"\t";
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

void constructInitState(int L,int d,MPS& result){
  constructM1(L,d,result);
  // And exchange the odd ones by identity on the virutal level tensor the pure state in the other
  // Pure state X+:
  complex_t Xplus[]={ONE_c,ONE_c,ONE_c,ONE_c}; // already the tho |X+><X+|
  mwArray Xp(Indices(d*d,1),Xplus);
  int D=2; // Dimension of the MPS
  mwArray idVirt=identityMatrix(D);idVirt.reshape(Indices(1,D*D));
  idVirt.multiplyLeft(reshape(Xp,Indices(d*d,1)));
  idVirt.reshape(Indices(d*d,D,D));
  for(int k=0;k<L-1;k++)
    result.setA(2*k+1,idVirt);
  // The last one does not pass on D
  complex_t data[]={ONE_c,ONE_c};
  idVirt=mwArray(Indices(1,D),data);
  idVirt.multiplyLeft(reshape(Xp,Indices(d*d,1)));
  idVirt.reshape(Indices(d*d,D,1));
  result.setA(2*L-1,idVirt);
}

void constructVarSzMPO(int L,int d,MPO& varSz,int chain){
  varSz.initLength(2*L);
  // individual operators I need: Id and Sz
  mwArray sig0=identityMatrix(d);
  sig0.reshape(Indices(d*d,1));
  complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
  mwArray sigZ(Indices(d*d,1),dataZ);
  // Now my basic operators are the tensor product of Id with itself and Sz with itself,
  mwArray IdId=sig0*Hconjugate(sig0);
  mwArray SzSz=sigZ*Hconjugate(sigZ);

  int D=2; // MPO bond dimension
  // Now, depending on the case, construct a single body MPO
  // case chain==0 or 1, the MPO only acts on half the chain
  // case chain==2, directly let SpinMPO construct the sum
  if(chain==2){
    SpinMPO::getSingleBodyMPO(2*L,IdId,SzSz,varSz);
    return;
  }
  if(chain!=0&&chain!=1){
    cout<<"ERROR: constructVarSzMPO does not know chain value "<<chain<<endl;
    exit(1);
  }
  // For chain 0 or 1 I just insert identities (idVirtual tensor
  // idPhysical) in the sites I do not want to involve and use the other ones
  MPO aux(3);
  SpinMPO::getSingleBodyMPO(3,IdId,SzSz,aux); // Just need three
					      // tensors and the identities
  mwArray idBoth=identityMatrix(D);idBoth.reshape(Indices(D*D,1));
  mwArray idPhys=identityMatrix(d*d);idPhys.reshape(Indices(1,d*d*d*d));
  idBoth.multiplyRight(idPhys);idBoth.reshape(Indices(D,D,d*d,d*d));
  idBoth.permute(Indices(3,1,4,2));
  Operator* opId=new Operator(idBoth);
  idPhys.reshape(Indices(d*d,1,d*d,1));
  Operator* opIdPhys=new Operator(idPhys);
  if(chain==0){ // Then I fill in the odd ones with Id (last only phys)
    varSz.setRotatedOp(0,&aux.getOp(0),Indices(1,2,3,4)); // trick to copy
    cout<<"Set operator 0 to aux[0]->"<<aux.getOp(0).getDimensions()<<endl; 
    varSz.setOp(1,opId,true);
    cout<<"Set operator 1 to Id->"<<opId->getDimensions()<<endl; 
    if(L>2){
      varSz.setRotatedOp(2,&aux.getOp(1),Indices(1,2,3,4)); // trick to copy
      cout<<"Set operator 2 to aux[1]->"<<aux.getOp(1).getDimensions()<<endl; 
    }
    for(int k=1;k<L-1;k++){
      varSz.setOp(2*k+1,opId,false);
      cout<<"Set operator "<<2*k+1<<" to Id->"<<opId->getDimensions()<<endl; 
    if(k>1){
	varSz.setOp(2*k,&varSz.getOp(1),false);
	cout<<"Set operator "<<2*k<<" to Op[1]->"<<varSz.getOp(1).getDimensions()<<endl; 
    }
    }
    varSz.setRotatedOp(2*(L-1),&aux.getOp(2),Indices(1,2,3,4)); // trick to copy
    cout<<"Set operator "<<2*(L-1)<<" to aux[2]->"<<aux.getOp(2).getDimensions()<<endl; 
    varSz.setOp(2*L-1,opIdPhys,true);
  }
  else{ // Then I "remove" the even ones (0,2,...)
    varSz.setOp(0,opIdPhys,true);
    varSz.setRotatedOp(1,&aux.getOp(0),Indices(1,2,3,4)); // trick to copy
    cout<<"Set operator 1 to aux[0]->"<<aux.getOp(0).getDimensions()<<endl; 
    varSz.setOp(2,opId,true);
    cout<<"Set operator "<<2<<" to Id->"<<opId->getDimensions()<<endl; 
    if(L>2){
      varSz.setRotatedOp(3,&aux.getOp(1),Indices(1,2,3,4)); // trick to copy
      cout<<"Set operator 3 to aux[1]->"<<aux.getOp(1).getDimensions()<<endl; 
    }
    for(int k=2;k<L;k++){
      if(k<L-1){
	varSz.setOp(2*k+1,&varSz.getOp(3),false);
	cout<<"Set operator "<<2*k+1<<" to Op[3]->"<<varSz.getOp(3).getDimensions()<<endl; 
      }
      varSz.setOp(2*k,opId,false);
      cout<<"Set operator "<<2*k<<" to Id->"<<opId->getDimensions()<<endl; 
    }
    varSz.setRotatedOp(2*L-1,&aux.getOp(2),Indices(1,2,3,4)); // trick to copy
    cout<<"Set operator "<<2*L-1<<" to aux[2]->"<<aux.getOp(2).getDimensions()<<endl; 
  }
}

void constructProjZ0(int L,int d,MPO& result,bool ancilla,bool both){
  // First, recover the basic Z0 projector onto a single chain of
  // length L
  result.initLength(2*L);
  MPO aux(L);
  SpinMPO::getProjectorTotalSz(L,d,aux);
  // Now, if both sectors are projected, I need to compose the
  // MPOs. Otherwise, just fill in with identities
  if(!both){
    // Auxiliary identity matrices, I get from dimensions of the
    // neighbors
    for(int k=0;k<L;k++){
      // Set the original MPO in place
      int dest=ancilla?2*k+1:2*k;
      mwArray op=aux.getOp(k).getFullData();
      result.setRotatedOp(dest,op,Indices(1,2,3,4)); // Copy them
      // the neighbor is set to identity
      int destN=ancilla?dest-1:dest+1;
      int dimVirt=ancilla?op.getDimension(1):op.getDimension(3);
      mwArray idOp=identityMatrix(d);idOp.reshape(Indices(d*d,1));
      mwArray idVirt=identityMatrix(dimVirt);idVirt.reshape(Indices(1,dimVirt*dimVirt));
      idOp.multiplyRight(idVirt);idOp.reshape(Indices(d,d,dimVirt,dimVirt));
      idOp.permute(Indices(1,3,2,4));
      result.setOp(destN,new Operator(idOp),true);
    }
  }
  else{
    // I need to compose the intermediate operators one by one
    // except first and last, which are simple
    result.setRotatedOp(0,aux.getOp(0).getFullData(),Indices(1,2,3,4));
    result.setRotatedOp(2*L-1,aux.getOp(L-1).getFullData(),Indices(1,2,3,4));
    for(int k=1;k<2*L-1;k++){
      int orig=k%2==0?k/2:(k-1)/2;
      mwArray opData=aux.getOp(orig).getFullData();
      Indices dimOp=opData.getDimensions();
      // Dim of extra Id for the other system: in case of even site,
      // that is the dim of the last ancilla operator, which was for
      // the site before, i.e. it is the Dl now
      int extraDl=k%2==0?dimOp[1]:dimOp[3];
      mwArray auxId=identityMatrix(extraDl);auxId.reshape(Indices(1,extraDl*extraDl));
      opData.reshape(Indices(dimOp[0]*dimOp[1]*dimOp[2]*dimOp[3],1));
      opData.multiplyRight(auxId);
      opData.reshape(Indices(dimOp[0],dimOp[1],dimOp[2],dimOp[3],extraDl,extraDl));
      // The new order of the indices is such that ancillas are above
      // (i.e. first)
      Indices newOrder(1,5,2,3,6,4);
      if(k%2!=0) newOrder=Indices(1,2,5,3,4,6);
      opData.permute(newOrder);
      opData.reshape(Indices(dimOp[0],dimOp[1]*extraDl,dimOp[2],dimOp[3]*extraDl));
      result.setOp(k,new Operator(opData),true);
    }
  }
}

void constructSzMPS(int L,int d,MPS& result,bool ancilla){
  // First, construct the MPO from the single chain Sz one
  MPO aux(L);
  SpinMPO::getSzMPO(L,d,aux);
  // Now insert identiti in between, depending on whether the ancilla
  // or the other systems are being measured
  MPO aux2(2*L);
  for(int k=0;k<L;k++){
    // Set the original MPO in place
    int dest=ancilla?2*k+1:2*k;
    mwArray op=aux.getOp(k).getFullData();
    aux2.setRotatedOp(dest,op,Indices(1,2,3,4),true); // Copy them and
						      // conjugate to
						      // compensate
						      // for later!
    // the neighbor is set to identity
    int destN=ancilla?dest-1:dest+1;
    int dimVirt=ancilla?op.getDimension(1):op.getDimension(3);
    mwArray idOp=identityMatrix(d);idOp.reshape(Indices(d*d,1));
    mwArray idVirt=identityMatrix(dimVirt);idVirt.reshape(Indices(1,dimVirt*dimVirt));
    idOp.multiplyRight(idVirt);idOp.reshape(Indices(d,d,dimVirt,dimVirt));
    idOp.permute(Indices(1,3,2,4));
    aux2.setOp(destN,new Operator(idOp),true);
  }
  // Finally, transform in a MPS
  MPSfromMPO(aux2,result,false); // down index ends up first, so that
				 // it is contracted with outgoing
				 // index of rho
}


void doubleMPO(const MPO& simpleMPO,MPO& doubleMPO){
  int L=simpleMPO.getLength();
  doubleMPO.initLength(L); // just in case
  for(int k=0;k<L;k++){
    mwArray aux=simpleMPO.getOp(k).getFullData();
    doubleMPO.setOp(k,new DoubleOperator(aux,permute(aux,Indices(3,2,1,4))),true);
  }
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

void extendTransposeMPO(const MPO& simpleMPO,MPO& doubleMPO,int dimA){
  int L=simpleMPO.getLength();
  doubleMPO.initLength(L); // just in case
  mwArray identPhys=identityMatrix(dimA);
  identPhys.reshape(Indices(dimA,1,dimA,1)); 
  for(int k=0;k<L;k++){  
    doubleMPO.setOp(k,new DoubleOperator(identPhys,permute(simpleMPO.getOp(k).getFullData(),Indices(3,2,1,4))),true);
  }
}


void constructM1P0(int L,int d, MPS& result){
  // Get the pieces: projectors and M1 exp
  MPO P0(L);
  SpinMPO::getProjectorTotalSz(L,d,P0);
  // I will compose, by hand, the projector and M1
  MPS auxL(L,2,d*d);
  constructM1single(L,d,auxL);
  int D=auxL.getBond();
  //Now, replace all tensors adequately
  cout<<"Prepared auxL MPS "<<auxL<<endl;
  // Now apply exactly the projector
  Contractor& contractor=Contractor::theContractor();
  {
    MPS tmp(auxL);
    MPO projSig(L);
    extendMPO(P0,projSig,d);
    contractor.optimize(projSig,auxL,tmp,D*(L/2+1));
    extendTransposeMPO(P0,projSig,d);
    contractor.optimize(projSig,tmp,auxL,D*(L/2+1));
  }
  // And replace the tensors in the original MPS, while also
  // substituting the ancilla ones by the projector times proper Id
  for(int k=0;k<L;k++){
    // pair of tensors involved in the pair
    mwArray m1tensor=auxL.getA(k).getA();
    mwArray p0tensor=P0.getOp(k).getFullData();
    int Dr_sigma=m1tensor.getDimension(2);
    int Dl_tau=p0tensor.getDimension(1);
    mwArray id_sigma=identityMatrix(Dl_tau);id_sigma.reshape(Indices(1,Dl_tau*Dl_tau));    
    mwArray id_tau=identityMatrix(Dr_sigma);id_tau.reshape(Indices(1,Dr_sigma*Dr_sigma));    
    // Now the new site for even is tensor prod of m1tensor with id_sigma
    Indices oldDims=m1tensor.getDimensions();
    m1tensor.reshape(Indices(-1,1));m1tensor.multiplyRight(id_sigma);
    m1tensor.reshape(Indices(oldDims[0],oldDims[1],oldDims[2],Dl_tau,Dl_tau));
    m1tensor.permute(Indices(1,2,4,3,5));
    m1tensor.reshape(Indices(oldDims[0],oldDims[1]*Dl_tau,oldDims[2]*Dl_tau));
    oldDims=p0tensor.getDimensions();
    p0tensor.reshape(Indices(-1,1));p0tensor.multiplyRight(id_tau);
    p0tensor.reshape(Indices(oldDims[0],oldDims[1],oldDims[2],oldDims[3],Dr_sigma,Dr_sigma));
    p0tensor.permute(Indices(1,3,5,2,6,4));
    p0tensor.reshape(Indices(oldDims[0]*oldDims[2],oldDims[1]*Dr_sigma,oldDims[3]*Dr_sigma));
    result.replaceSite(2*k,m1tensor,false); // prevent checking with neighbours
    result.replaceSite(2*k+1,p0tensor,false); // prevent checking with neighbours
 }
}

void constructM1single(int L,int d,MPS& result,bool cosi){
   mwArray sig0=identityMatrix(d);
   sig0.reshape(Indices(d*d,1,1));
   complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
   mwArray sigZ(Indices(d*d,1,1),dataZ);
   mwArray Z(Indices(2,d*d));
   for(int j=0;j<d*d;j++){
     Z.setElement(sig0.getElement(Indices(j,0,0)),Indices(0,j));
     Z.setElement(sigZ.getElement(Indices(j,0,0)),Indices(1,j));
   }
   int D=2; // Dimension of the MPS
   mwArray C(Indices(D,D,2)),C1(Indices(1,D,2)),CL(Indices(D,1,2));
   C.setElement(ONE_c,Indices(0,0,0));
   C1.setElement(ONE_c,Indices(0,0,0));
   C.setElement(ONE_c,Indices(D-1,D-1,0));
   CL.setElement(ONE_c,Indices(D-1,0,0));
   // There is a site-dependent element C(0,D-1,1)
   C1.setElement(ONE_c,Indices(0,D-1,1));
   if(cosi)
     CL.setElement(cos(2*M_PIl*(L-1.)/L)*ONE_c,Indices(0,0,1));
   else
     CL.setElement(exp((-2*M_PIl*(L-1.)/L)*I_c),Indices(0,0,1));
   C1.reshape(Indices(1*D,2));
   C1.multiplyRight(Z);
   C1.reshape(Indices(1,D,d*d));
   C1.permute(Indices(Indices(3,1,2)));
   result.setA(0,C1);
   CL.reshape(Indices(D*1,2));
   CL.multiplyRight(Z);
   CL.reshape(Indices(D,1,d*d));
   CL.permute(Indices(Indices(3,1,2)));
   result.setA(L-1,CL);
   for(int k=1;k<L-1;k++){
     if(cosi)
       C.setElement(cos(2*M_PIl*k/L)*ONE_c,Indices(0,D-1,1));
     else
       C.setElement(exp((-2*M_PIl*k/L)*I_c),Indices(0,D-1,1));
     mwArray C_(C);
     C_.reshape(Indices(D*D,2));
     C_.multiplyRight(Z);
     C_.reshape(Indices(D,D,d*d));
     C_.permute(Indices(Indices(3,1,2)));
     result.setA(k,C_);
   }
}
