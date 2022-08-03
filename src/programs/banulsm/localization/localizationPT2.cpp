#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0


#include "HeisenbergHamiltonian.h"

using namespace shrt;

/** localizationPT2 simulates the time evolution of an MPO under a
    disordered Heisenberg Hamiltonian \ref
    <HeisenbergHamiltonian>. Instead of generating the random fields,
    it reads it from the parameter file, to be able to check the same
    realization with different bond dimensions, etc.

    Starting from a certain MPO, it simulates the evolution and
    computes the average and localized magnetization along the way.
    To keep track of how good the evolution is, I will also compute
    the energy along the way.

    Receives arguments:
    \param <L> (int) number of sites
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian (Jx=Jy=Jz)
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian (local values 
                        randomly distributed between -h and h)
    \param <delta> (double) width of the time step used
    \param <M> (int) number of steps
    \param <rate> (int) frequency of recording data (nr of steps)
                       Results will be saved for 0,rate,2*rate...,M
    \param <D> (int) maximum bond dimension
    \param <inputparamfname> (char*) name of the file where the 
                            parameters of H were stored
    \param <outfname> (char*) name of the output file for the results
                   (will be replaced if it exists)
    \param <paramfname> (char*) name of the file where the 
                            parameters and the name of the MPS file
			    will be stored (to eventually continue with 
			    localizationPTcont)
    \param <mpsfile> (char*) name of the file where to save the MPS 
*/

/** Construct the initial state and observable, but vectorized*/
void constructM1(int L,int d,MPS& result);

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
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int rate=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* inputparamfname=argv[++cntr];
  const char* outfname=argv[++cntr];
  const char* paramfname=argv[++cntr];
  const char* mpsfile=argv[++cntr];

  cout<<"Initialized arguments "
      <<", L="<<L
      <<", J="<<J_
      <<", h="<<h_
      <<", delta="<<delta
      <<", M="<<M
      <<", D="<<D
      <<", outfile="<<outfname
      <<", paramfile="<<paramfname
      <<", mpsfile="<<mpsfile
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
  *out<<"% hs in "<<inputparamfname<<endl;
  *out<<"% MPS in "<<mpsfile<<endl;
  *out<<"% M\t delta\t t\t D\t <M1>\t tr(rho)\t tr(rho H)"<<endl;

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d=2;

  // Place for coefficients
  vector<double> J(L-1,J_);
  vector<double> h(L,h_);

  // First, read from the parameters file 
  ifstream* inP;
  inP=new ifstream(inputparamfname);
  if(!inP->is_open()){
    cout<<"Error: impossible to open file "<<inputparamfname<<
      " to read the parameters"<<endl;
    exit(1);
  }
  {
    int _L,_D,_M0;
    double _J,_h,_delta;
    char _mpsfile[200];
    // Basic configuration:
    inP->read((char*)&_L,sizeof(int));
    inP->read((char*)&_D,sizeof(int));
    inP->read((char*)&_J,sizeof(double));
    inP->read((char*)&_h,sizeof(double));
    inP->read((char*)&_delta,sizeof(double));
    inP->read((char*)&_M0,sizeof(int));
    for(int k=0;k<L;k++) inP->read((char*)&h[k],sizeof(double));
    // Name of the MPS file
    int lenName=0;
    inP->read((char*)&lenName,sizeof(int));
    inP->read(_mpsfile,lenName);
    inP->close();
    delete inP;
    // Check consistency
    if(_L!=L||_J!=J_||_h!=h_){
      cout<<"Error: parameters are inconsistent with contents of file"<<endl;
      exit(1);
    }
  }

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

  // Create the Hamiltonian
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
  
  // Now construct the initial rho MPO (term 1x...xsigmazx...x1)
  MPS M1(L,2,d*d);
  constructM1(L,d,M1);
  MPS rho0(M1);
  
  // Also the identity, to compute the trace
  MPS Id(L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++)
    Id.setA(k,id2);
  
  MPS hmps(L,1,d*d);
  const MPO& hamil=hamH.getHMPO();
  MPSfromMPO(hamil,hmps);
  
  //  M1.exportForMatlab("M1.m");
  //Id.exportForMatlab("Id.m");
  
  // Now do the evolution step by step, and compute the expectation value of M1 every time
  int cnt=0;
  double time=0.;
  complex_t M1ev=contractor.contract(rho0,M1);
  complex_t trR=contractor.contract(rho0,Id);
  complex_t Eev=contractor.contract(rho0,hmps);
  *out<<cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<M1ev<<"\t"<<trR<<"\t"<<Eev<<endl;
  while(cnt<M){
    //     cout<<"Time evolution step nr "<<cnt<<", time="<<time<<", value "<<M1ev<<", trace "<<trR<<endl;
    MPS aux(rho0); // temporary copy
    contractor.optimize(expHe2,aux,rho0,D);
    // cout<<"After expHe2, v'*H*v="<<contractor.contract(imaGS,hamil,imaGS)
    // 	   <<", v*v="<<contractor.contract(imaGS,imaGS)<<endl;
    contractor.optimize(expHo,rho0,aux,D);
    // cout<<"After expHo, v'*H*v="<<contractor.contract(imaGS,hamil,imaGS)
    // 	   <<", v*v="<<contractor.contract(imaGS,imaGS)<<endl;
    contractor.optimize(expHe2,aux,rho0,D);
    // cout<<"After expHe2, v'*H*v="<<contractor.contract(imaGS,hamil,imaGS)
    // 	   <<", v*v="<<contractor.contract(imaGS,imaGS)<<endl;
    //rho0.gaugeCond('R',true);
    M1ev=contractor.contract(rho0,M1);
    trR=contractor.contract(rho0,Id);
    Eev=contractor.contract(rho0,hmps);
    time+=delta;
    *out<<++cnt<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<M1ev<<"\t"<<trR<<"\t"<<Eev<<endl;
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

