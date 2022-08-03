#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "SpinMPO.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "HeisenbergAncillaryNoiseHamiltonian.h"

using namespace shrt;

/** localizationTI_Single simulates the evolution of a mixed state,
    which initially is maximally mixed for every site except for one.
    In this one, we start with stma_x, sigma_y or sigma_z, evolve, and
    save the evolved MPS after a certain number of steps.

    Receives arguments:
    \param <L> (int) number of sites (in the original chain)
    \param <L0> (int) number of sites (in the original chain)
    \param <isXY> (bool) whether the model to be considered is only XY (no ZZ term) 
    \param <J> (double) parameter \f$J\f$ of the Heisenberg
                  Hamiltonian for the ancillary chain
    \param <B> (double) parameter \f$B\f$ of the interaction
                        between both chains
    \param <init> (char) initial state in the middle, can be 1 (X),
                        2(Y) or 3 (Z)
    \param <delta> (double) width of the time step used
    \param <M> (int) number of steps
    \param <rate> (int) frequency of recording data (nr of steps)
                        Results will be saved for 0,rate,2*rate...,M
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the results
                        (will be replaced if it exists, unless newInstance==0)
    \param <mpsfile> (char*) name of the file where to save the MPS 
    \param <jobsdir> (char*) name of the dir to write the new jobs
    \param <outfname_traces> (char*) name of the output file for the trace norms results
    \param <normFreq> (int) frequency (steps) to compute trace norms
                        If 0 or none, no traces will be computed.
    \param <newInstance> (bool) whether to start a new evolution 
			(if not, the MPS in mpsfile is used)
    \param <M0> (int)  only read if newInstance is false, it
                       is the number of steps already done in the
		       existing file.
*/


/** Product of all identities, but single site operator Op on site pos */
void constructLocal(int L,int pos,mwArray op,int d,MPS& result);
/** Vectorize an MPO to make an MPS out of it. If up==true (default)
    the first component of the double physical index will be the one
    which was going up.*/
void MPSfromMPO(const MPO& mpo,MPS& mps,bool up=true);
void MPOfromMPS(const MPS& mps,MPO& mpo,bool up=true);

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

/** Substitute a single site in the MPS by the one corresponding to
    the single site operator initSt indicates*/
void setSingleInitState(MPS& rho0,int L,int d,int initSt,int pos);

/** Substitute L0 sites in the MPS by the state that encodes initSt in L0 sites, 
    starting at pos0. */
void setExtendedInitState(MPS& rho0,int L,int d,int initSt,int L0,int pos0);

/** Compute and print to ofstrem the expectation values of X, Y, Z on
    pos */
void computeLocalOp(const MPS& rho,int pos,int d,ofstream* out);

//double traceNorm(const MPS& rho);
// Just to do here the same as mblS2 does
void computeTraceNorms(const MPS& rho0,const char* outfname,int L,int L0,bool isXY,double J_,double B_,int initSt,
		       double delta,int M,int D);

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  int L0=atoi(argv[++cntr]);
  bool isXY=false;
  isXY=atoi(argv[++cntr])!=0;
  double J_=atof(argv[++cntr]);
  double B_=atof(argv[++cntr]);
  int initSt=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int rate=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  const char* mpsfile=argv[++cntr];
  const char* jobsdir=argv[++cntr];
  const char* outfname_Traces=argv[++cntr];
  int savingFreq=0;
  if(argc>11)
    savingFreq=atoi(argv[++cntr]);
  int newInstance=atoi(argv[++cntr]);
  if(newInstance<0) newInstance=1;
  int M0=0;
  if(!newInstance&&argc>13) M0=atoi(argv[++cntr]);


  cout<<"Initialized arguments "
      <<", L="<<L
      <<", L0="<<L0
      <<", isXY?="<<isXY
      <<", J="<<J_
      <<", B="<<B_
      <<", delta="<<delta
      <<", M="<<M
      <<", D="<<D
      <<", outfile="<<outfname
      <<", mpsfile="<<mpsfile
      <<", newInstance="<<newInstance
      <<", M0="<<M0
      <<", savingFreq="<<savingFreq
      <<", rate="<<rate
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
	<<" D="<<D<<", init St="<<initSt<<endl;
    *out<<"% MPS in "<<mpsfile<<endl;
    *out<<"% M\t delta\t t\t D\t tr(rho H)\t tr(rho) \t <X>\t <Y>\t <Z>"<<endl;
  }
  *out<<setprecision(15);
  int Mtot=M0+M;

  Contractor& contractor=Contractor::theContractor();
  //contractor.setConvTol(1E-10);
  //cout<<"Initialized Contractor"<<endl;
  int d=2;

  // First: create the Hamiltonian
  double J1x(1.),J1y(1.),J1z(isXY?0.:1.);
  double J2x(J_),J2y(J_),J2z(isXY?0.:J_);
  HeisenbergAncillaryNoiseHamiltonian hamH(L,J1x,J1y,J1z,J2x,J2y,J2z,B_);
  //hamH.getHMPO().exportForMatlab("hamilHeisAnc.m");  
  //cout<<"Created the Hamiltonian"<<endl;
  // if(L<10){
  //   // I could now compute the ground state to test
  //   MPS testGS(2*L,D,d);double lambda=0;
  //   contractor.setEigenSolver(primme);
  //   contractor.findGroundState(hamH.getHMPO(),D,&lambda,testGS);
  //   cout<<setprecision(15);
  //   cout<<"Found GS with energy "<<lambda<<endl;
  //   //    exit(1);
  // }

  // The evolution operator for a rho will result from double operators
  // Which MPOs do I need for the evolution
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

  // cout<<"Created all the exponential operators:"<<endl;
  // cout<<"expHeSigma:"<<expHSe<<endl;
  // cout<<"expHeTau:"<<expHTe<<endl;  
  // cout<<"expHoSigma:"<<expHSo<<endl;  
  // cout<<"expHoTau:"<<expHTo<<endl;  
  // cout<<"expHoSigma(1/2):"<<expHSo2<<endl;  
  // cout<<"expHoTau(1/2):"<<expHTo2<<endl;  
  // cout<<"expHB:"<<expHB<<endl; 

  // Now construct the initial rho MPO, in this case, the identity,
  // except for one site 
  MPS rho0(2*L,1,d*d);
  rho0.setProductState(p_maxent);
  int posCenter=L; // position in the chain where I set the different operator
  //  setSingleInitState(rho0,L,d,initSt,L);
  setExtendedInitState(rho0,L,d,initSt,L0,(L-L0)/2);

  if(!newInstance)
    rho0.importMPS(mpsfile);
  rho0.gaugeCond('R',1); // normalize
  
  // Now do the evolution step by step, and compute the energy,
  // and the expectation value of the local operators where we started
  MPO hamil(2*L);
  const MPO& _hamil=hamH.getHMPO();
  extendMPO(_hamil,hamil,d);

  // Also the identity, to compute the trace
  MPS Id(2*L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<2*L;k++)
    Id.setA(k,id2);

  int cnt=M0;
  double time=cnt*delta;
  complex_t Eev,trRho;
  Eev=contractor.contract(rho0,hamil,Id);
  trRho=contractor.contract(rho0,Id);
  //cout<<"Time evolution step nr "<<cnt<<", time="<<time
  //  <<", trace="<<trRho<<", <H>="<<Eev
  //  <<endl;
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"
      <<real(Eev)<<"\t"<<imag(Eev)<<"\t"<<real(trRho)<<"\t"<<imag(trRho)<<"\t";
  computeLocalOp(rho0,posCenter,d,out);
  // vector<complex_t> lambda2;
  // contractor.getSchmidtValues(rho0,lambda2); // lambda^2 in the middle
  // for(int id=0;id<D;id++)
  //   *out<<real(lambda2[id])<<"\t";
  // //  computeEntropies(rho0,out);
  // computeLocalSigZ(rho0,d,out);
  // //  computeErrors(rho0,out,Ds);
  *out<<endl;
  out->close(); // I probably don't want to keep it open for days!
  if(D<100){ // I can just run the traces in the same process and skip saving
    computeTraceNorms(rho0,outfname_Traces,L,L0,isXY,J_,B_,initSt,delta,cnt,D);
  }
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
    //    rho0.gaugeCond('R',1);
    Eev=contractor.contract(rho0,hamil,Id);
    trRho=contractor.contract(rho0,Id);
    //cout<<"Time evolution step nr "<<cnt<<", time="<<time
    //	<<", trace="<<trRho<<", <H>="<<Eev
    //	<<endl;
    out=new ofstream(outfname,ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"
	<<real(Eev)<<"\t"<<imag(Eev)<<"\t"<<real(trRho)<<"\t"<<imag(trRho)<<"\t";
    computeLocalOp(rho0,posCenter,d,out);
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
      //      if(D<100){ // I can just run the traces in the same process and skip saving
      computeTraceNorms(rho0,outfname_Traces,L,L0,isXY,J_,B_,initSt,delta,cnt,D);
      // }
      // else{
      // 	char name[120];
      // 	//cout<<"Saving step "<<cnt<<" to file"<<endl;
      // 	sprintf(name,"%s_%d",mpsfile,cnt);
      // 	rho0.exportMPS(name);
      // 	// Also, after exporting, I will launch the trace norm calculation!
      // 	stringstream jobfile;
      // 	jobfile<<jobsdir<<"/job_"<<initSt<<"_L"<<L<<"_J"<<J_<<"_B"<<B_<<"_D"<<D<<"_M"<<cnt<<".dat";
      // 	ofstream* ojob=new ofstream(jobfile.str().data());
      // 	if(!ojob->is_open()){
      // 	  cout<<"WARNING!!!! Could not open the JOB file. Launch by hand!!!!"<<endl;
      // 	  exit(1);
      // 	}
      // 	//cout<<"Saving next job to file "<<jobfile.str()<<endl;
      // 	int mem=6000; //if(D==80) mem=6000;
      // 	*ojob<<"#!/bin/bash"<<endl;
      // 	*ojob<<"msub_modified_tqo097 -N "<<(isXY?"XY":"heis")<<".M"<<cnt<<".L"<<L<<".J"<<J_<<".B"<<B_<<".D"<<D<<" -l h_vmem="
      // 	     <<mem<<"M -raw ";
      // 	*ojob<<"./mblS2 "<<L<<" "<<isXY<<" "<<J_<<" "<<B_<<" "<<initSt<<" "
      // 	     <<delta<<" "<<cnt<<" "<<D<<" "
      // 	     <<outfname_Traces<<" "<<name<<endl;
      // 	ojob->close();
      // 	delete ojob;
      // }
    }
  }

  rho0.exportMPS(mpsfile);
}


void constructLocal(int L,int pos,mwArray op,int d,MPS& result){
   mwArray sig0=identityMatrix(d);
   sig0.reshape(Indices(d*d,1,1));
   mwArray sigZ(op);sigZ.reshape(Indices(d*d,1,1));
   result=MPS(L,1,d*d);
   for(int k=0;k<L;k++){
     result.setA(k,sig0);
   }
   result.setA(pos,sigZ);
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

void setSingleInitState(MPS& rho0,int L,int d,int initSt,int pos){
  complex_t Xop[]={ZERO_c,ONE_c,ONE_c,ZERO_c}; 
  complex_t Yop[]={ZERO_c,I_c,-1.*I_c,ZERO_c}; 
  complex_t Zop[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c}; 
  mwArray X(Indices(d*d,1,1),Xop);
  mwArray Y(Indices(d*d,1,1),Yop);
  mwArray Z(Indices(d*d,1,1),Zop);
  switch(initSt){
  case 1:
    { cout<<"Setting init State to X"<<endl;
      rho0.setA(pos,X);
      break;}
  case 2:
    {  cout<<"Setting init State to Y"<<endl;
      rho0.setA(pos,Y);
      break;}
  case 3:
    {  cout<<"Setting init State to Z"<<endl;
      rho0.setA(pos,Z);
      break;}
  default:{
    cout<<"Error! Unknown initSt "<<initSt<<endl;
    exit(1);
  }
  }
}

// void setExtendedInitState(MPS& rho0,int L,int d,int initSt,int L0,int pos0){
//   int ind01=2; //complex_t op01[]={ZERO_c,ZERO_c,ONE_c,ZERO_c}; 
//   int ind10=1; //complex_t op10[]={ZERO_c,ONE_c,ZERO_c,ZERO_c}; 
//   int ind00=0; //complex_t op00[]={ONE_c,ZERO_c,ZERO_c,ZERO_c}; 
//   int ind11=3; //complex_t op11[]={ZERO_c,ZERO_c,ZERO_c,ONE_c}; 
//   int ind0,ind1;complex_t coeff0,coeff1;
//   switch(initSt){
//   case 1:
//     { cout<<"Setting init State to X"<<endl;
//       ind0=ind01;ind1=ind10;
//       coeff0=ONE_c;coeff1=ONE_c;
//       break;}
//   case 2:
//     {  cout<<"Setting init State to Y"<<endl;
//       ind0=ind01;ind1=ind10;
//       coeff0=-1.*I_c;coeff1=I_c;
//       break;}
//   case 3:
//     {  cout<<"Setting init State to Z"<<endl;
//       ind0=ind00;ind1=ind11;
//       coeff0=ONE_c;coeff1=-1.*ONE_c;
//       break;}
//   default:{
//     cout<<"Error! Unknown initSt "<<initSt<<endl;
//     exit(1);
//   }
//   }
//   // The original MPS is changed s.t. a shorter MPS is inserted (from pos0 to pos0+L0-1) with the right state.
//   // But the state is stretched over the ancillary sites!
//   //left edge 
//   int D=L0==1?1:2;
//   mwArray Al(Indices(d*d,rho0.getA(2*pos0).getDl(),D));Al.fillWithZero();
//   mwArray Am(Indices(d*d,D,D));Am.fillWithZero();
//   mwArray Ar(Indices(d*d,D,rho0.getA(2*(pos0+L0-1)).getDr()));Ar.fillWithZero();
//   Al.setElement(coeff0,Indices(ind0,0,0));Al.setElement(coeff1,Indices(ind1,0,D-1));
//   Am.setElement(ONE_c,Indices(ind0,0,0));Am.setElement(ONE_c,Indices(ind1,D-1,D-1));
//   Ar.setElement(ONE_c,Indices(ind0,0,0));Ar.setElement(ONE_c,Indices(ind1,D-1,0));
//   // And the special identity for the intermediate ancillas:
//   mwArray idA=identityMatrix(d);idA.reshape(Indices(d*d,1));
//   mwArray idVirt=identityMatrix(D);idVirt.reshape(Indices(1,D*D));
//   idA.multiplyRight(idVirt);idA.reshape(Indices(d*d,D,D));
//   //  rho0.setA(2*pos0,Al);rho0.setA(2*pos0+1,idA);
//   rho0.replaceSite(2*pos0,Al,false);rho0.replaceSite(2*pos0+1,idA,false);
//   for(int k=1;k<L0-1;k++){
//     //    rho0.setA(2*(pos0+k),Am);
//     rho0.replaceSite(2*(pos0+k),Am,false);
//     //rho0.setA(2*(pos0+k)+1,idA);
//     rho0.replaceSite(2*(pos0+k)+1,idA,false);
//   }
//   //  rho0.setA(pos0+L0-1,Ar);
//   rho0.replaceSite(2*(pos0+L0-1),Ar,false);
// }

void setExtendedInitState(MPS& rho0,int L,int d,int initSt,int L0,int pos0){
   mwArray sig0=identityMatrix(d);
   sig0.reshape(Indices(d*d,1,1));
   complex_t data00[]={ONE_c,ZERO_c,ZERO_c,ZERO_c};
   mwArray op00(Indices(d*d,1,1),data00);
   complex_t data01[]={ZERO_c,ZERO_c,ONE_c,ZERO_c};
   mwArray op01(Indices(d*d,1,1),data01);
   complex_t data10[]={ZERO_c,ONE_c,ZERO_c,ZERO_c};
   mwArray op10(Indices(d*d,1,1),data10);
   complex_t data11[]={ZERO_c,ZERO_c,ZERO_c,ONE_c};
   mwArray op11(Indices(d*d,1,1),data11);
   mwArray Z(Indices(4,d*d));
   for(int j=0;j<d*d;j++){
     Z.setElement(op00.getElement(Indices(j,0,0)),Indices(0,j));
     Z.setElement(op01.getElement(Indices(j,0,0)),Indices(1,j));
     Z.setElement(op10.getElement(Indices(j,0,0)),Indices(2,j));
     Z.setElement(op11.getElement(Indices(j,0,0)),Indices(3,j));
   }
   // Set the product of identities on the edges
   for(int k=0;k<pos0;k++)
     rho0.replaceSite(k,sig0,false);
   for(int k=pos0+L0;k<L;k++)
     rho0.replaceSite(k,sig0,false);
   int ind1(0),ind2(3);
   complex_t coeff1(ONE_c),coeff2(-1*ONE_c);
   if(initSt<3){
     ind1=1;ind2=2;
     coeff1=ONE_c;coeff2=ONE_c;
     if(initSt==2){
       coeff1=-1*I_c;coeff2=I_c;
     }
   }
   int D=2; // Dimension of the MPS
   mwArray C(Indices(D,D,4)),C1(Indices(1,D,4)),CL(Indices(D,1,4));
   C.setElement(ONE_c,Indices(0,0,ind1));
   C.setElement(ONE_c,Indices(D-1,D-1,ind2));
   C1.setElement(coeff1,Indices(0,0,ind1));
   C1.setElement(coeff2,Indices(0,D-1,ind2));
   CL.setElement(ONE_c,Indices(0,0,ind1));
   CL.setElement(ONE_c,Indices(D-1,0,ind2));

   // And the special identity for the intermediate ancillas:
   mwArray idA=identityMatrix(d);idA.reshape(Indices(d*d,1));
   mwArray idVirt=identityMatrix(D);idVirt.reshape(Indices(1,D*D));
   idA.multiplyRight(idVirt);idA.reshape(Indices(d*d,D,D));
   
   C1.reshape(Indices(1*D,4));C1.multiplyRight(Z);C1.reshape(Indices(1,D,d*d));C1.permute(Indices(3,1,2));
   rho0.replaceSite(2*pos0,C1,false);rho0.replaceSite(2*pos0+1,idA,false);
   C.reshape(Indices(D*D,4));C.multiplyRight(Z);C.reshape(Indices(D,D,d*d));C.permute(Indices(3,1,2));
   for(int k=pos0+1;k<pos0+L0-1;k++){
     rho0.replaceSite(2*k,C,false);rho0.replaceSite(2*k+1,idA,false);
   }
   CL.reshape(Indices(D*1,4));CL.multiplyRight(Z);CL.reshape(Indices(D,1,d*d));CL.permute(Indices(3,1,2));
   rho0.replaceSite(2*(pos0+L0-1),CL,false);
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

void computeLocalOp(const MPS& rho,int pos,int d,ofstream* out){
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d*d,1,1),dataX);
  complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigY(Indices(d*d,1,1),dataY);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d*d,1,1),dataZ);

  // cout<<"Computing local expectation values "<<endl;
  mwArray sig0=identityMatrix(d);
  sig0.reshape(Indices(d*d,1,1));
  int L=rho.getLength();
  MPS oper(L,1,d*d);
  for(int k=0;k<L;k++){
    oper.setA(k,sig0);
  }

  Contractor& contractor=Contractor::theContractor();
  complex_t resL;
  oper.setA(pos,sigX);
  resL=contractor.contract(rho,oper);
  *out<<real(resL)<<"\t"<<imag(resL)<<"\t";
  oper.setA(pos,sigY);
  resL=contractor.contract(rho,oper);
  *out<<real(resL)<<"\t"<<imag(resL)<<"\t";
  oper.setA(pos,sigZ);
  resL=contractor.contract(rho,oper);
  *out<<real(resL)<<"\t"<<imag(resL)<<"\t";
}






void constructLocal(int L,int pos,int initSt,int d,MPS& result){
   mwArray sig0=identityMatrix(d);
   sig0.reshape(Indices(d*d,1,1));
  complex_t Xop[]={ZERO_c,ONE_c,ONE_c,ZERO_c}; 
  complex_t Yop[]={ZERO_c,I_c,-1.*I_c,ZERO_c}; 
  complex_t Zop[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c}; 
  mwArray X(Indices(d*d,1,1),Xop);
  mwArray Y(Indices(d*d,1,1),Yop);
  mwArray Z(Indices(d*d,1,1),Zop);
  mwArray* op;
  switch(initSt){
  case 1:
    { op=&X;
      break;}
  case 2:
    { op=&Y;
      break;}
  case 3:
    { op=&Z;
      break;}
  default:{
    cout<<"Error! Unknown initSt "<<initSt<<endl;
    exit(1);
  }
  }
  op->reshape(Indices(d*d,1,1));
  result=MPS(L,1,d*d);
  for(int k=0;k<L;k++){
    result.setA(k,sig0);
  }
  result.setA(pos,*op);
}



void MPOfromMPS(const MPS& mps,MPO& mpo,bool up){
  int L=mps.getLength();

  mpo.initLength(L);

  for(int k=0;k<L;k++){
    mwArray aux=mps.getA(k).getA();
    Indices dims=aux.getDimensions();
    int d0=sqrt(dims[0]);
    if(d0*d0!=dims[0]){
      cout<<"Error: Dimension of site "<<k<<" ("<<dims[0]
	  <<") does not seem a square=> don't know how t divide it"<<endl;
      exit(1);
    }
    aux.reshape(Indices(d0,d0,dims[1],dims[2]));
    if(up)
      aux.permute(Indices(1,3,2,4));
    else
      aux.permute(Indices(2,3,1,4));
    mpo.setOp(k,new Operator(aux),true);
  }

}

double traceNorm(const MPS& rho){
  MPO rhoMPO(rho.getLength());
  MPOfromMPS(rho,rhoMPO);
  // if(rho.getLength()==4){
  //   rhoMPO.exportForMatlab("rhoTest4.m");
  //   exit(1);
  // }
  mwArray Op;
  //cout<<"Computing trace norm of mps, length "<<rho.getLength()<<endl;
  expandOper(rhoMPO,Op);
  // cout<<"Operator expanded to dims "<<Op.getDimensions()<<endl;
  // trace norm, is the sum of singular values
  mwArray U,S,Vdag;
  wrapper::svd(Op,U,S,Vdag);
  return real(S.trace());
}

int d=2;
void computeTraceNorms(const MPS& rho0,const char* outfname,int L,int L0,bool isXY,double J_,double B_,int initSt,
		       double delta,int M,int D){

  ofstream* out;
  bool appendingFile=file_exists(outfname)&&(M!=0);
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
	<<" D="<<D<<", init St="<<initSt<<endl;
    *out<<"% M\t delta\t t\t D\t tr(rho H)\t tr(rho(t)) \t tr(sigma"<<initSt<<")\t |rho(t)_2|_1"<<endl;
  }
  *out<<setprecision(15);

  Contractor& contractor=Contractor::theContractor();
  //cout<<"Initialized Contractor"<<endl;

  // Also the identity, to compute the trace
  MPS Id(2*L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<2*L;k++)
    Id.setA(k,id2);

  //  cout<<"Analyzing state in file "<<mpsfile<<endl;
  //  MPS rho(L,1,1);
  // rho.importMPS(mpsfile);
  MPS rho(rho0); // copying
  rho.gaugeCond('R',1);

  complex_t trRho,Eev,sigev;
  {
    double J1x(1.),J1y(1.),J1z(isXY?0.:1.);
    double J2x(J_),J2y(J_),J2z(isXY?0.:J_);
    HeisenbergAncillaryNoiseHamiltonian hamH(L,J1x,J1y,J1z,J2x,J2y,J2z,B_);
    //    HeisenbergAncillaryNoiseHamiltonian hamH(L,J_,B_); 
    const MPO& hamil=hamH.getHMPO();
    MPS hmps(2*L,1,d*d);
    MPSfromMPO(hamil,hmps);
    Eev=contractor.contract(rho,hmps);
    trRho=contractor.contract(rho,Id);
    constructLocal(2*L,L,initSt,d,hmps);
    sigev=contractor.contract(rho,hmps);
  }

  *out<<M<<"\t"<<delta<<"\t"<<M*delta<<"\t"<<D<<"\t"<<real(Eev)<<"\t"<<imag(Eev)<<"\t"
      <<real(trRho)<<"\t"<<imag(trRho)<<"\t"
      <<real(sigev)<<"\t"<<imag(sigev)<<"\t"; // 10 columns for nothing (3 is time)


  // Now I need to start tracing out things, to construct reduced operators and to compute their trace norms.
  // I will keep 10,8,6,4,L0 sites of the primary chain (cols 11,12,13,14,15)
  // and 4(8), 2(4) of both (cols 16,17)


  for(int cut=5;cut>=L0/2;cut--){
    vector<int> toTrace;
    // 1) trace all but the central block of 2*cut spins, in the first chain
    int lastTracedL=max(L/2-cut-1,-1);int firstTracedR=min(L/2+cut,L);
    for(int k=0;k<=lastTracedL;k++){
      toTrace.push_back(2*k);
      toTrace.push_back(2*k+1);
    }
    for(int k=lastTracedL+1;k<firstTracedR;k++){
      toTrace.push_back(2*k+1); // all the ancillas
    }
    for(int k=firstTracedR;k<L;k++){
      toTrace.push_back(2*k);
      toTrace.push_back(2*k+1);
    }
    //cout<<"To get the central "<<cut*2<<" spins (L="<<L<<"), tracing out "<<toTrace<<endl;
    MPS redRho(rho);
    redRho.traceOutSites(toTrace);
    double tN12=traceNorm(redRho);
    // And finally write the results
    *out<<(firstTracedR-1-lastTracedL)<<"\t";
    *out<<tN12<<"\t";
  }
  for(int cut=2;cut>=L0/2;cut--){
    vector<int> toTrace;
    // 1) trace all but the central block of 2*cut spins, in both chains
    int lastTracedL=max(L/2-cut-1,-1);int firstTracedR=min(L/2+cut,L);
    for(int k=0;k<=lastTracedL;k++){
      toTrace.push_back(2*k);
      toTrace.push_back(2*k+1);
    }
    for(int k=firstTracedR;k<L;k++){
      toTrace.push_back(2*k);
      toTrace.push_back(2*k+1);
    }
    //    cout<<"To get the central "<<cut*2<<" spins (L="<<L<<"), tracing out "<<toTrace<<endl;
    MPS redRho(rho);
    redRho.traceOutSites(toTrace);
    double tN12=traceNorm(redRho);
    // And finally write the results
    *out<<(firstTracedR-1-lastTracedL)*2<<"\t";
    *out<<tN12<<"\t";
  }

  //cout<<"Computation finished"<<endl;
  *out<<endl;
  out->close();
  delete out;


}
