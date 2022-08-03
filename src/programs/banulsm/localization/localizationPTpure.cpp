#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

#include "HeisenbergHamiltonian.h"

using namespace shrt;

/** localizationPTpure simulates the time evolution of an MPO under a
    disordered Heisenberg Hamiltonian \ref <HeisenbergHamiltonian>.
    Starting from a certain MPO, in this case a pure state, 
    \f$\rho_0=|\Psi_0 \rangle|langle \Psi_0 |\f$, it
    simulates the evolution and computes the average and localized
    magnetization along the way.  To keep track of how good the
    evolution is, I will also compute the energy along the way.

    Receives arguments:
    \param <L> (int) number of sites
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian (Jx=Jy=Jz)
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian (local values 
                        randomly distributed between -h and h)			
    \param <initState> (int) which initial tho to consider. Options are:
                        0 - \f$|\Psi_0\rangle=|0\rangle^{\otimes L}\f$;
                        1 - \f$|\Psi_0\rangle=|1\rangle^{\otimes L}\f$;			
                        5 - \f$|\Psi_0\rangle=|\alpha\rangle^{\otimes L}\f$, 
			    random, but the same on every site;
			9 - \f$|GS_{J=0}\rangle\f$, GS of the single particle term.
    \param <delta> (double) width of the time step used
    \param <M> (int) number of steps
    \param <rate> (int) frequency of recording data (nr of steps)
                       Results will be saved for 0,rate,2*rate...,M
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the results
                   (will be replaced if it exists)
    \param <paramfname> (char*) name of the output file for the 
                            parameters and the name of the MPS file
    \param <mpsfile> (char*) name of the file where to save the MPS 
    \param <seed> (int) a seed for the random number generator (if it is 
                            zero, TIME will be used).
*/

/** Construct the initial state and observable, but vectorized*/
void constructM1(int L,int d,MPS& result);

/** Idem for the one with the exponential complex coefficients */
void constructM1exp(int L,int d,MPS& result);

/** Product of all identities, but sigz on site pos */
void constructLocal(int L,int pos,int d,MPS& result);

/** Construct the MPO version (vectorized) of the pure product initial state 
    determined by input parameter initSt.
 */
void constructPureM1(int L,int d,MPS& result,int initSt,const vector<double>& h);

void computeLocalSigZ(const MPS& rho,int d,ofstream* out);
void computeEntropies(const MPS& rho,ofstream* out);
void computeErrors(const MPS& rho,ofstream* out,const vector<int>& Ds);

/** Vectorize an MPO to make an MPS out of it. If up==true (default)
    the first component of the double physical index will be the one
    which was going up.*/
void MPSfromMPO(const MPO& mpo,MPS& mps,bool up=true);

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double J_=atof(argv[++cntr]);
  double h_=atof(argv[++cntr]);
  int initSt=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int rate=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  const char* paramfname=argv[++cntr];
  const char* mpsfile=argv[++cntr];
  int seed=0;
  if(argc>11) seed=atoi(argv[++cntr]);

  cout<<"Initialized arguments "
      <<", L="<<L
      <<", J="<<J_
      <<", h="<<h_
      <<", initSt="<<initSt
      <<", delta="<<delta
      <<", M="<<M
      <<", D="<<D
      <<", outfile="<<outfname
      <<", paramfile="<<paramfname
      <<", mpsfile="<<mpsfile
      <<", seed="<<seed
      <<endl;

  ofstream* out;
  out=new ofstream(outfname);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<"% L="<<L<<", J="<<J_<<", h="<<h_
      <<" D="<<D<<endl;
  *out<<"% hs in "<<paramfname<<endl;
  *out<<"% MPS in "<<mpsfile<<endl;

  vector<int> Ds; // vector of reference Ds to compare to and determine errors
  for(int Dref=10;Dref<D;Dref+=10){
    Ds.push_back(Dref);
  }

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d=2;

  // First: create the Hamiltonian
  // Set the random coefficients for the magnetic field
  if(seed==0)
    srandom(time(NULL));
  else
    srandom(seed);
  // Set random coefficients
  vector<double> J(L-1,J_);
  vector<double> h(L,h_);
  for(int k=0;k<L;k++){
    h[k]=(2.*random()/RAND_MAX-1.)*h_;
  }
  cout<<"h="<<h<<endl;

  *out<<"% Random values (seed "<<seed<<") h:"<<endl;
  *out<<setprecision(18);
  for(int k=0;k<L;k++) *out<<"% h("<<k+1<<")="<<h[k]<<endl;
  *out<<setprecision(15);
  *out<<"% M\t delta\t t\t D\t <M1>\t tr(rho)\t tr(rho H)"<<endl;


  ofstream* outP;
  outP=new ofstream(paramfname);
  if(!outP->is_open()){
    cout<<"Error: impossible to open file "<<paramfname<<
      " for parameters"<<endl;
    exit(1);
  }
  // Basic configuration:
  outP->write((char*)&L,sizeof(int));
  outP->write((char*)&D,sizeof(int));
  outP->write((char*)&J_,sizeof(double));
  outP->write((char*)&h_,sizeof(double));
  outP->write((char*)&delta,sizeof(double));
  outP->write((char*)&M,sizeof(int));
  for(int k=0;k<L;k++) outP->write((char*)&h[k],sizeof(double));
  // Name of the MPS file
  int lenName=strlen(mpsfile);
  outP->write((char*)&lenName,sizeof(int));
  outP->write(mpsfile,lenName);
  outP->close();
  delete outP;

  HeisenbergHamiltonian hamH(L,J,J,J,h,d);
  cout<<"Created the Hamiltonian"<<endl;
  // The evolution operator for a rho will result from double operators
  // Which MPOs do I need for the indices in both sides
  MPO expHe2(L),expHe(L),expHo(L);
  // I construct the simple ones first, and then from them the DoubleOperators
  {
    // Regular ones
    MPO _expHe2(L),_expHe(L),_expHo(L);
    MPO _expHe2A(L),_expHeA(L),_expHoA(L); // for the second indices
    hamH.getExponentialMPOeven(_expHe2,-delta*.5*I_c);
    hamH.getExponentialMPOeven(_expHe,-delta*I_c);
    hamH.getExponentialMPOodd(_expHo,-delta*I_c);
    hamH.getExponentialMPOeven(_expHe2A,delta*.5*I_c);
    hamH.getExponentialMPOeven(_expHeA,delta*I_c);
    hamH.getExponentialMPOodd(_expHoA,delta*I_c);
    for(int k=0;k<L;k++){
      expHe2.setOp(k,new DoubleOperator(_expHe2.getOp(k).getFullData(),_expHe2A.getOp(k).getFullData()),true);
      expHe.setOp(k,new DoubleOperator(_expHe.getOp(k).getFullData(),_expHeA.getOp(k).getFullData()),true);
      expHo.setOp(k,new DoubleOperator(_expHo.getOp(k).getFullData(),_expHoA.getOp(k).getFullData()),true);
    }
    // The simple ones will be destroyed now. I already have them copied
  }
  cout<<"Created all the exponential operators"<<endl;
  
  // Now construct the initial rho MPO (term 1x...xsigmazx...x1)
  MPS rho0(L,2,d*d);
  MPS M1(L,2,d*d);
  constructM1(L,d,M1);
  constructPureM1(L,d,rho0,initSt,h);

  // Also the identity, to compute the trace
  MPS Id(L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++)
    Id.setA(k,id2);
  
  MPS hmps(L,1,d*d);
  const MPO& hamil=hamH.getHMPO();
  MPSfromMPO(hamil,hmps);

  // And I will also need the sigmaz operator on every site, which I can compute ony by one  

  //  M1.exportForMatlab("M1.m");
  //Id.exportForMatlab("Id.m");
  
  // Now do the evolution step by step, and compute the expectation value of M1 every time
  int cnt=0;
  double time=0.;
  complex_t M1ev=contractor.contract(rho0,M1);
  complex_t trR=contractor.contract(rho0,Id);
  complex_t Eev=contractor.contract(rho0,hmps);
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<real(M1ev)<<"\t"<<imag(M1ev)
      <<"\t"<<real(trR)<<"\t"<<imag(trR)
      <<"\t"<<real(Eev)<<"\t"<<imag(Eev)<<"\t";
  vector<complex_t> lambda2;
  contractor.getSchmidtValues(rho0,lambda2); // lambda^2 in the middle
  for(int id=0;id<D;id++)
    *out<<real(lambda2[id])<<"\t";
  //  computeEntropies(rho0,out);
  computeLocalSigZ(rho0,d,out);
  computeErrors(rho0,out,Ds);
  *out<<endl;
  while(cnt<M){
    MPS aux(rho0); // temporary copy
    contractor.optimize(expHe2,aux,rho0,D);
    // Now apply r-1 times the pair Ho He
    int cntLoop=0;
    while(cntLoop<rate-1&&cnt<M-1){
      contractor.optimize(expHo,rho0,aux,D);
      //cout<<"Time evolution step nr "<<cnt+1<<", time="<<time<<", value "<<M1ev<<", trace "<<trR<<endl;
      contractor.optimize(expHe,aux,rho0,D);
      cnt++;cntLoop++;
    }
    contractor.optimize(expHo,rho0,aux,D);
    contractor.optimize(expHe2,aux,rho0,D);
    //cout<<"Time evolution step nr "<<cnt<<", time="<<time<<", value "<<M1ev<<", trace "<<trR<<endl;

    M1ev=contractor.contract(rho0,M1);
    trR=contractor.contract(rho0,Id);
    Eev=contractor.contract(rho0,hmps);
    time+=delta;
    *out<<++cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<real(M1ev)<<"\t"<<imag(M1ev)
	<<"\t"<<real(trR)<<"\t"<<imag(trR)
	<<"\t"<<real(Eev)<<"\t"<<imag(Eev)<<"\t";
    contractor.getSchmidtValues(rho0,lambda2); // lambda^2 in the middle
    for(int id=0;id<D;id++)
      *out<<real(lambda2[id])<<"\t";
    //    computeEntropies(rho0,out);
    computeLocalSigZ(rho0,d,out);
    computeErrors(rho0,out,Ds);
    *out<<endl;
    // export the MPS if cnt is multiple of 100
    if(cnt%100==0){
      char name[120];
      sprintf(name,"%s_%d",mpsfile,cnt);
      rho0.exportMPS(name);
    }
  }

  rho0.exportMPS(mpsfile);

  out->close();
  delete out;
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
   mwArray C(Indices(D,D,2)),C1(Indices(1,D,2)),CL(Indices(D,1,2));
   C.setElement(ONE_c,Indices(0,0,0));
   C1.setElement(ONE_c,Indices(0,0,0));
   C.setElement(ONE_c,Indices(D-1,D-1,0));
   CL.setElement(ONE_c,Indices(D-1,0,0));
   // There is a site-dependent element C(0,D-1,1)
   C1.setElement(ONE_c,Indices(0,D-1,1));
   //CL.setElement(exp((-2*M_PIl*(L-1.)/L)*I_c),Indices(0,0,1));
   CL.setElement(cos((M_PIl*(L-1.)/L))*ONE_c,Indices(0,0,1));
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
     //C.setElement(exp((-2*M_PIl*k/L)*I_c),Indices(0,D-1,1));
     C.setElement(cos((M_PIl*k/L))*ONE_c,Indices(0,D-1,1));
     mwArray C_(C);
     C_.reshape(Indices(D*D,2));
     C_.multiplyRight(Z);
     C_.reshape(Indices(D,D,d*d));
     C_.permute(Indices(Indices(3,1,2)));
     result.setA(k,C_);
   }

}

void constructM1exp(int L,int d,MPS& result){
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
   mwArray C(Indices(D,D,2)),C1(Indices(1,D,2)),CL(Indices(D,1,2));
   C.setElement(ONE_c,Indices(0,0,0));
   C1.setElement(ONE_c,Indices(0,0,0));
   C.setElement(ONE_c,Indices(D-1,D-1,0));
   CL.setElement(ONE_c,Indices(D-1,0,0));
   // There is a site-dependent element C(0,D-1,1)
   C1.setElement(ONE_c,Indices(0,D-1,1));
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
     C.setElement(exp((-2*M_PIl*k/L)*I_c),Indices(0,D-1,1));
     mwArray C_(C);
     C_.reshape(Indices(D*D,2));
     C_.multiplyRight(Z);
     C_.reshape(Indices(D,D,d*d));
     C_.permute(Indices(Indices(3,1,2)));
     result.setA(k,C_);
   }

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

void constructPureM1(int L,int d,MPS& result,int initSt,const vector<double>& h){
  result=MPS(L,1,d*d);
  if(initSt!=9){
    mwArray basicA(Indices(d,1)); // local vector
    switch(initSt){
    case 1:{
      cout<<"Pure state |Z->"<<endl;
      basicA.setElement(ONE_c,Indices(1,0));
      break;
    }
    case 5:{
      cout<<"Pure (TI) random state "<<endl;
      basicA.fillRandom();
      basicA=(1./norm(basicA))*basicA;
      break;
    }
    default:
      cout<<"Default pure state |Z+>"<<endl;
      basicA.setElement(ONE_c,Indices(0,0));
    }
    mwArray aux=conjugate(basicA);
    aux.transpose();
    basicA.multiplyRight(aux); // this is my basic piece now
    cout<<"Basic tensor now is "<<basicA<<endl;
    basicA.reshape(Indices(d*d,1,1));
    for(int k=0;k<L;k++){
      result.setA(k,basicA);
    }
  }
  else{
    mwArray basicA0(Indices(d,d));
    basicA0.setElement(ONE_c,Indices(0,0));
    basicA0.reshape(Indices(d*d,1,1));
    mwArray basicA1(Indices(d,d));
    basicA1.setElement(ONE_c,Indices(1,1));
    basicA1.reshape(Indices(d*d,1,1));
    for(int k=0;k<L;k++){
      if(h[k]<0)
	result.setA(k,basicA0);
      else
	result.setA(k,basicA1);
    }
  }
}
