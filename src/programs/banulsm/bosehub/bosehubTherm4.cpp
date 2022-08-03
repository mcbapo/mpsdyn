
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "DoubleOperator.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "BoseHubbardHamiltonian.h"

/** bosehubTherm4 tries to find the thermal state of the Bose-Hubbard
    Hamiltonian \ref <BoseHubbardHamiltonian> using purification and
    imaginary time evolution with an MPO expression for the exponential
    operator.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <N> (int) cutoff for the occupation number
    \param <t> (double) parameter \f$t\f$ of the Hamiltonian
    \param <U> (double) parameter \f$U\f$ of the Hamiltonian
    \param <mu> (double) parameter \f$\mu\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <beta> (double) inverse temperature, \f$\beta\f$
    \param <delta> (double) one half of the width of time step (actually, 
    this is the number in the exponent of Hamiltonian terms, but appears 
    in both sides, so the "step" is 2*delta
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)

*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  int N=atoi(argv[++cntr]);
  double t=atof(argv[++cntr]);
  double U=atof(argv[++cntr]);
  double mu=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  double beta=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  cout<<"Initialized arguments: L="<<L
      <<", t="<<t
      <<", U="<<U
      <<", mu="<<mu
      <<", Tmin="<<1/beta
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% L="<<L<<", t="<<t<<", U="<<U<<", mu="<<mu
	<<", Nmax="<<N<<", D="<<D<<endl;
    *out<<"% delta\t nrsteps\t beta\t T\t Energy"
	<<endl;
    *out<<setprecision(10);
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);

  int M=beta/(2*delta); // Nr steps

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d0=N+1;
  int d=d0*d0;
  MPS thS(L,D,d);
  thS.setProductState(p_maxent); // the initial state, maximally entangled
  cout<<"Initialized identity state"<<endl;

  //int Dop=5; // bond dimension of Bose-Hubbard H

  //  double offset=0.; // done by Contractor, if needed
  BoseHubbardHamiltonian hBH(L,N,t,U,mu);
  //SchwingerHamiltonian hSchM(L,mu,x,L0M,alpha);
  const MPO& hamil=hBH.getHMPO();
  cout<<"Initialized basic Hamiltonian"<<endl;

  complex_t deltaC=(complex_t){-delta,0.}; // delta  as complex
  complex_t deltaC_2=(complex_t){-delta*.5,0.}; // delta/2
  // The required DoubleOperators
  DoubleOperator *HfirstL,*HfirstL_2,*HfirstR,*HfirstR_2;
  DoubleOperator *HL,*HL_2,*HR,*HR_2;
  DoubleOperator *HlastL,*HlastL_2,*HlastR,*HlastR_2;
  DoubleOperator* doubleId;

  { // I do not need to keep the mwArrays, for they are copied in the
    // FoldedOperators.
    // Now get the required terms for the MPO of the exponential
    mwArray Op1first,Op2first,Op1last,Op2last,Op1,Op2;
    // The term starting on pos 0
    hBH.getTwoBodyTermExponential(Op1first,Op2first,deltaC,0);
    // term starting on first to last
    hBH.getTwoBodyTermExponential(Op1last,Op2last,deltaC,L-2);
    // Basic term
    hBH.getTwoBodyTermExponential(Op1,Op2,deltaC,1);
    // Also the ones for the exponential with delta/2
    // These ones I only need for one of the terms in H, choose H_even 
    // (i.e., the one starting with 01)
    mwArray Op1first_2,Op2first_2,Op1last_2,Op2last_2,Op1_2,Op2_2;
    // The term starting on pos 0
    hBH.getTwoBodyTermExponential(Op1first_2,Op2first_2,deltaC_2,0);
    // term starting on first to last
    hBH.getTwoBodyTermExponential(Op1last_2,Op2last_2,deltaC_2,L-2);
    // Basic term
    hBH.getTwoBodyTermExponential(Op1_2,Op2_2,deltaC_2,L-2);
    
    // Construct the corresponding DoubleOperators
    HfirstL=new DoubleOperator(Op1first,conjugate(Op1first));
    HfirstR=new DoubleOperator(Op2first,conjugate(Op2first));
    HfirstL_2=new DoubleOperator(Op1first_2,conjugate(Op1first_2));
    HfirstR_2=new DoubleOperator(Op2first_2,conjugate(Op2first_2));

    HL=new DoubleOperator(Op1,conjugate(Op1));
    HR=new DoubleOperator(Op2,conjugate(Op2));
    HL_2=new DoubleOperator(Op1_2,conjugate(Op1_2));
    HR_2=new DoubleOperator(Op2_2,conjugate(Op2_2));

    HlastL=new DoubleOperator(Op1last,conjugate(Op1last));
    HlastR=new DoubleOperator(Op2last,conjugate(Op2last));
    HlastL_2=new DoubleOperator(Op1last_2,conjugate(Op1last_2));
    HlastR_2=new DoubleOperator(Op2last_2,conjugate(Op2last_2));

  }
  // now an auxiliary MPO, just identities (for the H only)
  mwArray ident=identityMatrix(d0);
  ident.reshape(shrt::Indices(d0,1,d0,1)); // just dx1xdx1
  doubleId=new DoubleOperator(ident,ident);
  // Simple MPO: H
  MPO foldedH(L); // For H (easy one)
  for(int k=0;k<L;k++){
    // Remember left one is reflected
    foldedH.setOp(k,new FoldedOperator(ident,hamil.getOpData(k)),true);
  }
  // And now assemble the MPOs
  MPO expHe(L),expHe_2(L),expHo(L); // exponential terms
  // set the first and second by hand
  expHo.setOp(0,doubleId);expHo.setOp(1,HL);
  expHe.setOp(0,HfirstL);expHe.setOp(1,HfirstR);
  expHe_2.setOp(0,HfirstL_2);expHe_2.setOp(1,HfirstR_2);
  // For the last and first-to-last there are two possibilities:
  if((L%2)==0){ // last one (L-1) is odd
    expHo.setOp(L-1,doubleId);expHo.setOp(L-2,HR);
    expHe.setOp(L-1,HlastR);expHe.setOp(L-2,HlastL);
    expHe_2.setOp(L-1,HlastR_2);expHe_2.setOp(L-2,HlastL_2);    
  }
  else{ // last one is even
    expHo.setOp(L-1,HlastR);expHo.setOp(L-2,HlastL);
    expHe.setOp(L-1,doubleId);expHe.setOp(L-2,HR);
    expHe_2.setOp(L-1,doubleId);expHe_2.setOp(L-2,HR_2);
  }
  // for the rest, loop
  for(int k=2;k<L-2;k++){
    bool even=(k%2)==0;
    if(even){ 
      expHe.setOp(k,HL);expHe_2.setOp(k,HL_2);
      expHo.setOp(k,HR);
    }
    else{ 
      expHe.setOp(k,HR);expHe_2.setOp(k,HR_2);
      expHo.setOp(k,HL);
    }
    cout<<"Set DoubleOperators at position "<<k<<endl;
  }

  // for(int k=0;k<L;k++){
  //cout<<"In expHe, op "<<k<<" has dims "<<expHe.getOp(k).getDimensions()<<endl;
  //}

  //expHe.exportMPOtext("/home/banulsm/transfer/mpoExpHe.txt");
  //expHe_2.exportMPOtext("/home/banulsm/transfer/mpoExpHe_2.txt");
  //expHo.exportMPOtext("/home/banulsm/transfer/mpoExpHo.txt");


  // Check the ground state?
  // double lambda=0.;
  // MPS gS(L,D,d);gS.setRandomState();
  // gS.gaugeCond('R',1);
  // contractor.findGroundState(foldedH,D,&lambda,gS);
  //cout<<"Ground state found with eigenvalue "<<lambda<<endl;

  // exit(1);

  // Now try an optimization
  int cnt=0;
  double betai=0;

  int stepBlk=5; //blocks of steps

  while(cnt<M){
    MPS aux(thS);
    contractor.optimize(expHe_2,aux,thS,D);
    int k=0;
    while((k<stepBlk-1)&&(cnt<M-1)){
      aux=thS;
      contractor.optimize(expHo,aux,thS,D);
      aux=thS;
      contractor.optimize(expHe,aux,thS,D);
      betai+=2*delta;
      k++;cnt++;
    }
    aux=thS;
    contractor.optimize(expHo,aux,thS,D);
    aux=thS;
    contractor.optimize(expHe_2,aux,thS,D);
    betai+=2*delta;
    cnt++;
    // Now I need to normalize
    thS.gaugeCond('R',1);
    complex_t energy=contractor.contract(thS,foldedH,thS);
    *out<<delta<<"\t"<<cnt<<"\t"<<betai<<"\t"<<1/betai<<"\t"
	<<energy<<endl;
    cout<<"At the end of step "<<cnt<<", beta="<<betai<<", T="<<1/betai
	<<", energy="<<energy<<endl;
  }
  
  // At the very end, delete everything
  delete HL,HR,HL_2,HR_2,HfirstL,HfirstR,HfirstL_2,HfirstR_2,HlastL,HlastR,HlastL_2,HlastR_2;
  delete doubleId;
}
