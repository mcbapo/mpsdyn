
#include <math.h>
#include <iomanip>

#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

/** Runs the findGroundState routine with the MPO for the
    Ising Hamiltonian \ref <IsingHamiltonian> squared, to check 
    states with min variance for fixed bond dimension.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D1> (int) maximum bond dimension
    \param <D2> (int) maximum bond dimension
    \param <deltaD> (int) step change in bond dim
    \param <outfname> (char*) name of the output file for the energy
    \param <app> (int) whether the output file is to be kept (app==1)
*/


int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double J=atof(argv[++cntr]);
  double g=atof(argv[++cntr]);
  double h=atof(argv[++cntr]);
  int D1=atoi(argv[++cntr]);
  int D2=atoi(argv[++cntr]);
  int incrD=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  cout<<"Initialized arguments: L="<<L
      <<", J="<<J
      <<", g="<<g
      <<", h="<<h
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D1="<<D1
      <<", D2="<<D2
      <<", deltaD="<<incrD<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% N\t J\t g\t h\t D\t <H^2>\t <H_L^2>\t <H_R^2> \t <H>\t <HL>\t <HR>"<<endl;
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  cout<<setprecision(10);

  int D=D1; // to start
  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  cout<<"Initialized Contractor"<<endl;
  int d=2;
  MPS gs(L,D,d);
  gs.setRandomState(); // the intial state, random
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  cout<<"Initialized random state, norm "<<contractor.contract(gs,gs)<<endl;
  
  // First: put H and find the GS
  IsingHamiltonian hamI(L,d,J,g,h); 
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamI.getHMPO(0.);
  cout<<"Constructed the hamil MPO"<<endl;
  // And now two MPOs with half the Hamil on each and the rest identities
  MPO HL(L),HR(L);
  mwArray Id=identityMatrix(d);Id.reshape(Indices(d,1,d,1));
  Operator opId(Id); // Since all TI, only need edges and middle
  for(int k=0;k<L;k++){
    if(k<L/2){
      if(k==L/2-1)
	HL.setOp(k,&hamil.getOp(L-1),false); // last one
      else
	HL.setOp(k,&hamil.getOp(k),false);
      HR.setOp(k,&opId,false);
    }
    else{
      HL.setOp(k,&opId,false);
      if(k==L/2)
	HR.setOp(k,&hamil.getOp(0),false); // first one
      else
	HR.setOp(k,&hamil.getOp(k),false);
    }
  }
  
  //hamil.exportForMatlab("isingHmpo.m");
  cout<<"Initial value of E, with initial state"<<endl;
  cout<<contractor.contract(gs,hamil,gs)<<endl;

  
  const MPO* ptrs[]={&hamil,&hamil};
  MPO ham2(L);
  MPO::join(2,ptrs,ham2);// MPO for H^2
  //  MPO ham2(2,ptrs); 
  
  bool done=0;
  double expH2=L*100.;
  // place for other params
  complex_t E,EL,varL,ER,varR;
  while(!done){ // find GS of H^2 until D is large enough
    contractor.findGroundState(ham2,D,&expH2,gs);
    // compute also half-chain variances
    E=contractor.contract(gs,hamil,gs);
    EL=contractor.contract(gs,HL,gs);
    ER=contractor.contract(gs,HR,gs);
    varL=contractor.contract2(gs,HL,gs); // since they are hermitian, I contract H^+ H
    varR=contractor.contract2(gs,HR,gs);
    cout<<"State for D="<<D<<" found with variance "<<expH2
	<<" energy="<<E<<", EL="<<EL<<", ER="<<ER<<", varL="<<varL
	<<", varR="<<varR<<endl;
     
    //    *out<<"% N\t J\t g\t h\t D\t <H^2>\t <H_L^2>\t <H_R^2> \t <H>\t <HL>\t <HR>"<<endl;
    *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t"
	<<expH2<<"\t"<<varL<<"\t"<<varR<<"\t"
	<<E<<"\t"<<EL<<"\t"<<ER<<"\t";
    *out<<endl;
     
    D+=incrD;
    if(D>D2) done=true;
    else gs.increaseBondDimension(D);
  }
  
  out->close();
  delete out;
}
