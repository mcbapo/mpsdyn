#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

#include "IsingHamiltonian.h"

using namespace shrt;

/** modulatedEnergyIsing simulates the time evolution of an MPO under the
    Ising Hamiltonian \ref <IsingHamiltonian>.
    In this case, we start from the modulated energy density as initial state,
    we simulate the evolution of the mixed state, and compute tr(rho rho(t))
    Additionally, we compute also the energy.

    Receives arguments:
    \param <L> (int) number of sites
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian 
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian 
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <delta> (double) width of the time step used
    \param <M> (int) number of steps
    \param <rate> (int) frequency of recording data (nr of steps)
                       Results will be saved for 0,rate,2*rate...,M
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the results
                   (will be replaced if it exists)
    \param <mpsfile> (char*) name of the file where to save the MPS 
*/

/** Construct the initial state, vectorized*/
void constructInitialState(int L,int d,double J,double g,double h,MPS& result);

/** Idem for the one with the exponential complex coefficients */
void constructM1exp(int L,int d,MPS& result);

/** Product of all identities, but sigz on site pos */
void constructLocal(int L,int pos,int d,MPS& result);

void computeLocalSigZ(const MPS& rho,int d,ofstream* out);
void computeEntropies(const MPS& rho,ofstream* out);
void computeErrors(const MPS& rho,ofstream* out,const vector<int>& Ds);

/** Vectorize an MPO to make an MPS out of it. If up==true (default)
    the first component of the double physical index will be the one
    which was going up.*/
void MPSfromMPO(const MPO& mpo,MPS& mps,bool up=true);

/** I could try with higher Trotter order, just using the recursive formula from second one */
//#define ORDER4 1
#define ORDER4_C 1

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double J_=atof(argv[++cntr]);
  double g_=atof(argv[++cntr]);
  double h_=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int rate=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  const char* mpsfile=argv[++cntr];

  cout<<"Initialized arguments "
      <<", L="<<L
      <<", J="<<J_
      <<", g="<<g_
      <<", h="<<h_
      <<", delta="<<delta
      <<", M="<<M
      <<", D="<<D
      <<", outfile="<<outfname
      <<", mpsfile="<<mpsfile
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
  *out<<"% MPS in "<<mpsfile<<endl;
  *out<<"% Output organized in "<<10+2*L<<" columns as follows "<<endl;
  *out<<"% 1) M step nr "<<endl;
  *out<<"% 2) delta time step (const)"<<endl;
  *out<<"% 3) time (delta*M) "<<endl;
  *out<<"% 4) D bond dimension (const) "<<endl;
  *out<<"% 5) re[tr(rho0 rho(t))]"<<endl;
  *out<<"% 6) im[tr(rho0 rho(t))]"<<endl;
  *out<<"% 7) re[tr(rho(t))]"<<endl;
  *out<<"% 8) im[tr(rho(t))]"<<endl;
  *out<<"% 9) re[tr(H rho(t))]"<<endl;
  *out<<"% 10) im[tr(H rho(t))]"<<endl;
  *out<<"% 11) re[tr(sigma_z[0] rho(t))]"<<endl;
  *out<<"% 12) im[tr(sigma_z[0] rho(t))]"<<endl;
  *out<<"% And similar (site by site) until column "<<10+2*L<<endl;
  vector<int> Ds; // vector of reference Ds to compare to and determine errors
  for(int Dref=10;Dref<D;Dref+=10){
    Ds.push_back(Dref);
  }

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d=2;

  // First: create the Hamiltonian
  *out<<setprecision(15);
  *out<<"% M\t delta\t t\t D\t <M1>\t tr(rho)\t tr(rho H)"<<endl;

  IsingHamiltonian hamH(L,d,J_,g_,h_);
  hamH.registerTimeDependence();
  cout<<"Created the Hamiltonian"<<endl;
  // The evolution operator for a rho will result from double operators
  // Which MPOs do I need for the indices in both sides
#ifdef ORDER4  // 4th order Trotter S2(sx) S2((1-2s)x) S2(sx) with s=1/(2-2^(1/3))
  cout<<"Order Trotter 4"<<endl;
  double sconst=1.351207191959657;
  MPO expHs(L); // the one with sx
  MPO expH2s(L); // the one with (1-2s)x
  // I construct the simple ones first, and then from them the DoubleOperators
  {
    // Regular ones
    MPO _expH(L);
    MPO _expHA(L); // for the second indices
    hamH.getUMPO(_expH,delta*sconst,0.,false); // not imaginary time
    hamH.getUMPO(_expHA,-delta*sconst,0.,false);
    for(int k=0;k<L;k++){
      expHs.setOp(k,new DoubleOperator(_expH.getOp(k).getFullData(),permute(_expHA.getOp(k).getFullData(),Indices(3,2,1,4))),true);
    }
    // and also for the second one
    hamH.getUMPO(_expH,delta*(1-2*sconst),0.,false); // not imaginary time
    hamH.getUMPO(_expHA,-delta*(1-2*sconst),0.,false);
    for(int k=0;k<L;k++){
      expH2s.setOp(k,new DoubleOperator(_expH.getOp(k).getFullData(),permute(_expHA.getOp(k).getFullData(),Indices(3,2,1,4))),true);
    }
    // The simple ones will be destroyed now. I already have them copied
  }
#else
#ifdef ORDER4_C
  cout<<"Order Trotter 4 complex"<<endl;
  complex_t Apar=.5*ONE_c+(sqrt(3)/6.)*I_c;
  MPO expHa(L); // the one with ax
  MPO expHaC(L); // the one with conj(a)x
  // I construct the simple ones first, and then from them the DoubleOperators
  {
    complex_t deltaC=-delta*I_c; // 'x' in Suzuki-trotter
    // Regular ones
    MPO _expH(L);
    MPO _expHA(L); // for the second indices
    hamH.getUMPO(_expH,Apar*deltaC,0.);
    hamH.getUMPO(_expHA,-conjugate(Apar)*deltaC,0.);
    for(int k=0;k<L;k++){
      expHa.setOp(k,new DoubleOperator(_expH.getOp(k).getFullData(),permute(_expHA.getOp(k).getFullData(),Indices(3,2,1,4))),true);
    }
    // and also for the second one
    hamH.getUMPO(_expH,conjugate(Apar)*deltaC,0.);
    hamH.getUMPO(_expHA,-Apar*deltaC,0.);
    for(int k=0;k<L;k++){
      expHaC.setOp(k,new DoubleOperator(_expH.getOp(k).getFullData(),permute(_expHA.getOp(k).getFullData(),Indices(3,2,1,4))),true);
    }
    // The simple ones will be destroyed now. I already have them copied
  }

#else
  cout<<"Order Trotter 2"<<endl;
  MPO expH(L);
  // I construct the simple ones first, and then from them the DoubleOperators
  {
    // Regular ones
    MPO _expH(L);
    MPO _expHA(L); // for the second indices
    hamH.getUMPO(_expH,delta,0.,false); // not imaginary time
    hamH.getUMPO(_expHA,-delta,0.,false);
    for(int k=0;k<L;k++){
      expH.setOp(k,new DoubleOperator(_expH.getOp(k).getFullData(),permute(_expHA.getOp(k).getFullData(),Indices(3,2,1,4))),true);
    }
    // The simple ones will be destroyed now. I already have them copied
  }
#endif
#endif

  cout<<"Created all the exponential operators"<<endl;
  
  // Now construct the initial rho MPO (term 1x...xcos() h_12x...x1)
  MPS rho(L,2,d*d);
  constructInitialState(L,d,J_,g_,h_,rho);

  MPS A(rho); // for measuring
    
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
  double normF=pow(2,L);
  complex_t Aev=contractor.contract(rho,A);
  complex_t trR=contractor.contract(rho,Id);
  complex_t Eev=contractor.contract(rho,hmps);
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<real(Aev)/normF<<"\t"<<imag(Aev)/normF
      <<"\t"<<real(trR)<<"\t"<<imag(trR)
      <<"\t"<<real(Eev)<<"\t"<<imag(Eev)<<"\t";
  // vector<complex_t> lambda2;
  // contractor.getSchmidtValues(rho,lambda2); // lambda^2 in the middle
  // for(int id=0;id<D;id++){
  //   if(id<lambda2.size())
  //     *out<<real(lambda2[id])<<"\t";
  //   else
  //     *out<<"0.\t";
  // }

  //  computeEntropies(rho0,out);
  computeLocalSigZ(rho,d,out);
  //computeErrors(rho,out,Ds);
  *out<<endl;
  while(cnt<M){
    // Now apply r times the evolution operator
    int cntLoop=0;
    while(cntLoop<rate&&cnt<M){
      cout<<"Inside the loop for "<<rate<<" iterations, countLoop="<<cntLoop<<", total cnt="<<cnt<<" (M="<<M<<")"<<endl;
      MPS aux(rho); // temporary copy
#ifdef ORDER4
      contractor.optimize(expHs,aux,rho,D);
      contractor.optimize(expH2s,rho,aux,D);
      contractor.optimize(expHs,aux,rho,D);
#else
#ifdef ORDER4_C
      contractor.optimize(expHaC,rho,aux,D);
      contractor.optimize(expHa,aux,rho,D);
#else
      contractor.optimize(expH,aux,rho,D);
#endif
#endif
      //cout<<"Time evolution step nr "<<cnt+1<<", time="<<time<<", value "<<M1ev<<", trace "<<trR<<endl;
      cnt++;cntLoop++;
      time+=delta;
    }

    Aev=contractor.contract(rho,A);
    trR=contractor.contract(rho,Id);
    Eev=contractor.contract(rho,hmps);
    *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<real(Aev)/normF<<"\t"<<imag(Aev)/normF
	<<"\t"<<real(trR)<<"\t"<<imag(trR)
	<<"\t"<<real(Eev)<<"\t"<<imag(Eev)<<"\t";
    // contractor.getSchmidtValues(rho,lambda2); // lambda^2 in the middle
    // for(int id=0;id<D;id++)
    //   if(id<lambda2.size())
    // 	*out<<real(lambda2[id])<<"\t";
    //   else *out<<"0.\t";
    //    computeEntropies(rho0,out);
    computeLocalSigZ(rho,d,out);
    //computeErrors(rho,out,Ds);
    *out<<endl;
    cout<<"Time evolution step nr "<<cnt<<", time="<<time<<", value "<<Aev/normF<<", trace "<<trR<<endl;
    // // export the MPS if cnt is multiple of 100
    // if(cnt%100==0){
    //   char name[120];
    //   sprintf(name,"%s_%d",mpsfile,cnt);
    //   rho.exportMPS(name);
    // }
  }

  rho.exportMPS(mpsfile);

  out->close();
  delete out;
}

void constructInitialState(int N,int d,double J,double g,double h,MPS& result){
  // basic spin operators appearing
  mwArray sig0=identityMatrix(d);
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sig3(Indices(d,d),dataz);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig1(Indices(d,d),datax);
  mwArray Z(Indices(3,d,d));Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      //Z.setElement(sig1.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      Z.setElement(sig3.getElement(Indices(i1,i2)),Indices(1,i1,i2));
      //Z.setElement(sig3.getElement(Indices(i1,i2)),Indices(2,i1,i2));
      Z.setElement(sig1.getElement(Indices(i1,i2)),Indices(2,i1,i2));
    }
  Z.reshape(Indices(3,d*d));
  // Now prepare the MPO but instead of TI as in the IsingHamiltonian, make it site dependent
  // and set it as an MPS before returning
  //MPO mpo(N);
  int D=3; // Dimension of the MPS
  for(int k=0;k<N;k++){
    int Dl=D; int Dr=D;
    if(k==0) Dl=1;
    if(k==N-1) Dr=1;
    mwArray C(Indices(Dl,Dr,3));C.fillWithZero();
    complex_t effJ=(J<0)?sqrt(abs(J))*I_c:sqrt(J)*ONE_c;
    //    complex_t effG=g*cos(k*M_PIl/N)*ONE_c;
    complex_t effG=g*cos(k*M_PIl/(N-1))*ONE_c;
    //    complex_t effH=h*cos(k*M_PIl/N)*ONE_c;
    complex_t effH=h*cos(k*M_PIl/(N-1))*ONE_c;
    //cout<<"Coefficients: effJ="<<effJ<<", effG="<<effG<<", effH="<<effH
    //  <<endl;
    // Identity
    if(k!=N-1)
      C.setElement(ONE_c,Indices(0,0,0)); // just pass 1 (before anything)
    if(k!=0)
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // just pass 1 (after)
    // sigmax terms (in all of them)
    C.setElement(effG,Indices(0,Dr-1,2));
    // sigmaz terms (in all of them)
    C.setElement(effH,Indices(0,Dr-1,1));
    // left sigz term (all but last)
    if(k!=N-1)
      C.setElement(cos((2*k-1)*M_PIl/(2*(N-1)))*effJ,Indices(0,1,1));
    //      C.setElement(cos(k*M_PIl/N)*effJ,Indices(0,1,1));
    // right sigz term (all but first)
    if(k!=0)
      C.setElement(effJ,Indices(1,Dr-1,1));
    
    // Now contract and give the operator
    // Now contract MPS and operators and set proper index ordering
    C.reshape(Indices(Dl*Dr,3));
    //cout<<"C for site "<<k<<"="<<C<<endl;
    mwArray res=C*Z;
    //cout<<"Computed res "<<res<<endl;
    res.reshape(Indices(Dl,Dr,d,d));
    //    res.permute(Indices(3,1,4,2));
    //cout<<"  with size: "<<res.getDimensions()<<endl;
    // Create and store an Operator
    //    mpo.setOp(k,new Operator(res),true);
    res.permute(Indices(3,4,1,2));
    res.reshape(Indices(d*d,Dl,Dr));
    result.replaceSite(k,res,false);
  }

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

