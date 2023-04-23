/**
   \file lightcone2.cpp Find the time evolution of a local excitation,
   starting in the ground state of some Hamiltonian, and using a
   light-cone + folding approach.  In this case, I use the lightcone
   only in that I apply the Hamiltonian to the center and leave the
   rest as at the beginning, as in flatLightCone
   
   \author Mari-Carmen Banuls
   \date 8/05/2012

*/


#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "IsingHamiltonian.h"
#include "FoldedOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

using namespace std;
using namespace shrt;

void extendFoldedState(MPS& result,const MPS& foldedS,
		       const mwArray& extraSite1,const mwArray& extraSite2,
		       int d,int Dl);
void computeExpectation(MPS& foldedS,const mwArray& oper,ofstream& out,
			int d);

/** 
    lightcone2 tries to approximate the evolution after a local
    perturbation using a "light-cone" approach.
    Arguments:
    \param <L> (int) length of the chain used to find the ground state
    \param <Lb> (int) initial number of sites to the sides of the
    center used in the evolution
    \param <J> (double) Parameter of the Ising Hamiltonian
    \param <g> (double) Parameter of the Ising Hamiltonian
    \param <h> (double) Parameter of the Ising Hamiltonian
    \param <D2> (int) maximum bond dimension in the (foldd) evolved state
    \param <D0> (int) maximum bond dimension used for the ground state
    \param <delta> (double) step width
    \param <M> (int) number of time steps
    \param <outfname> (char*) name of the output file for the magnetization
    \param <app> (int) whether the output file is to be kept (app==1)
    \param <outfname2> (char*) name of the output file for the entropies
    \param <app2> (int) whether the output file is to be kept (app2==1)
*/
 
int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  int Lb=atoi(argv[++cntr]);
  double J=atof(argv[++cntr]);
  double g=atof(argv[++cntr]);
  double h=atof(argv[++cntr]);
  int D2=atoi(argv[++cntr]);
  int D0=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);
  const char* outfname2=argv[++cntr];
  int app2=atoi(argv[++cntr]);
  
  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% N\t J\t g\t h\t D\t Energy"<<endl;
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  ofstream* out2;
  if(!app2){
    out2=new ofstream(outfname2);
  }
  else
    out2=new ofstream(outfname2,ios::app);
  if(!out2->is_open()){
    cout<<"Error: impossible to open file "<<outfname2<<
      " for output"<<endl;
    exit(1);
  }

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d=2;
  //int D0=D2/2; // D for GS
  int Dl=D0; // D of the single term

  int len=1+2*Lb; // length of the Hamiltonian applied
  MPS foldSt0(L,D2,d);
  {
    if(L%2==0) L++; // I want an odd number of sites, to use the
		    // middle one
    MPS gs(L,D0,d);
    gs.setRandomState(); // the intial state, random
    gs.gaugeCond('R',1);
    gs.gaugeCond('L',1);
    cout<<"Initialized random state "<<endl;
  
    // First, I compute the GS of the Ising H in a "long" system, as a
    // substitute for iTEBD
    double E0=0.;
    IsingHamiltonian hamI(L,d,J,g,h);
    cout<<"Created the Hamiltonian"<<endl;
    const MPO& hamil=hamI.getHMPO(0.);
    cout<<"Constructed the hamil MPO"<<endl;
    //   hamil.exportMPOtext("flatHmpo.txt");
    cout<<"Initial value, with initial state"<<endl;
    cout<<contractor.contract(gs,hamil,gs)<<endl;
    
    contractor.findGroundState(hamil,D0,&E0,gs);
    cout<<"Ground state found with eigenvalue "<<E0<<endl;
    *out<<"% "<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"
	<<D0<<"\t"<<E0<<"\t";
    *out<<endl;

    gs.fold(foldSt0); // folded initial state
  }
  //  cout<<"Constructed center MPS: "<<center<<endl;
  // cout<<"Constructed center2 MPS: "<<center2<<endl;

  MPS foldSt(foldSt0,D2);

  // Apply the local perturbation (now, strictly local)
  complex_t sigmaX_[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),sigmaX_);
  // Since I apply sigmax, I shouldn't worry about the norm (only about
  // the gauge condition, but this will be fixed after folding)
  foldSt.applyLocalOperator(0,sigX,false);
  //  cout<<"Constructed perturbed foldSt"<<endl;

  // I need an Identity operator on the edges
  mwArray id_=identityMatrix(d);
  id_.reshape(Indices(d,1,d,1));
  FoldedOperator idOp(id_,id_);
    
  // First, get the central site
  int L0=(L-1)/2;
  int count=0; // nr of step
  double t=0.;// H is constant
  int oT=2; // Trotter order
  bool imag=false; // no imaginary time evolution
  int Ltot=foldSt.getLength();
  while(count<M){
    // Get the Hamiltonian for a finite chain
    IsingHamiltonian shortH(len,d,J,g,h);
    // and the evolution MPOs
    MPO evol(len);
    shortH.getUMPO(evol,delta,t,imag,oT);
    double time=count*(delta+1);
    cout<<"Time evolution step nr "<<count+1<<", time="<<delta*(count+1)<<endl;
    // *****************    
    // fold the evolution
    MPO fEvol0((len+1)/2);
    evol.fold(fEvol0);
    MPS aux(foldSt); // temporary copy
    // Now I do the trick and enlarge it with Identity on the edges
    MPO fEvol(Ltot);
    if(fEvol0.getLength()<Ltot){
      for(int k=0;k<fEvol0.getLength();k++)
	fEvol.setOp(k,&fEvol0.getOp(k),false);
      for(int k=fEvol0.getLength();k<Ltot;k++){
	//cout<<"Setting fold MPO op site "<<k<<" to "<<idOp<<endl;
	fEvol.setOp(k,&idOp,false);
      }
    // Apply (folded) evolution
      cout<<"Applying cut evolution to "<<fEvol0.getLength()<<endl;
      contractor.optimize(fEvol,aux,foldSt,D2);
    }
    else{
      cout<<"Applying full evolution"<<endl;
      contractor.optimize(fEvol0,aux,foldSt,D2);
    }
    // Compute something
    cout<<"Computing magnetization"<<endl;
    computeExpectation(foldSt,sigX,*out,d);
    // Compute also the entropies (pair by pair)
    cout<<"Computing entropies"<<endl;
    for(int k=0;k<Ltot;k++){
      *out2<<contractor.getEntropy(foldSt,k,k)<<"\t";
    }
    *out2<<endl;

    if((len+1)/2<Ltot) len+=2;
    count++;
  }

  out->close();
  out2->close();
  delete out,out2;
}

void extendFoldedState(MPS& result,const MPS& foldedS,
		       const mwArray& extraSite1,const mwArray& extraSite2,
		       int d,int Dl){
  cout<<"extendFoldedState"<<endl;
  int len=foldedS.getLength();
  mwArray last=foldedS.getA(len-1).getA(); // this has both sqrt(rhos)
  Indices dims=last.getDimensions();
  last.reshape(Indices(Dl*Dl,Dl*Dl));
  mwArray extraSite1_(extraSite1);
  extraSite1_.permute(Indices(2,1,3));
  extraSite1_.reshape(Indices(Dl*Dl,d*d*Dl*Dl));
  last.multiplyRight(extraSite1_);
  last.reshape(Indices(Dl*Dl,d*d,Dl*Dl));
  last.permute(Indices(2,1,3));
  result=MPS(len+1,foldedS.getBond(),d*d);
  for(int k=0;k<len-1;k++){
    result.replaceSite(k,foldedS.getA(k).getA(),false);
  }
  result.replaceSite(len-1,last,false);
  result.replaceSite(len,extraSite2,false);
}

void computeExpectation(MPS& state,const mwArray& oper,ofstream& out,int d){
  int len=state.getLength(); // of the folded one!
  MPO ops(len);
  mwArray id_=identityMatrix(d);
  mwArray oper_(oper); oper_.reshape(Indices(d,1,d,1));
  id_.reshape(Indices(d,1,d,1));
  FoldedOperator identAll(id_,id_);
  for(int k=1;k<len;k++) ops.setOp(k,&identAll,false);
  // Center is special!
  // Center?
  int dcent=state.getA(0).getd();
  Operator* idCenter;
  if(dcent==d){ // folded with one in the middle
    idCenter=new Operator(id_);
  }
  else{ // folded with nothing in the middle
    mwArray fakeId(Indices(1,1,1,1));
    fakeId.setElement(ONE_c,0);
    idCenter=new Operator(fakeId);
  }
  ops.setOp(0,idCenter,false);

  //cout<<"Before starting sweeping, state "<<state<<", idMPO "<<ops<<endl;
  // Now sweep computing expectation values
  Contractor& contractor=Contractor::theContractor();
  double norm=real(contractor.contract(state,state));
  FoldedOperator leftOp(oper_,id_);
  FoldedOperator rightOp(id_,oper_);
  for(int pos=len-1;pos>0;pos--){
    ops.setOp(pos,&leftOp,false);
    // Compute and write result
    complex_t valk=contractor.contract(state,ops,state);
    out<<real(valk)/norm<<"\t";
    ops.setOp(pos,&identAll,false);
  }
  // Central one: different if it is double or not
  if(dcent==d){ // otherwise, there is nothing there
    // although I should fold the oper, as side dims are 1, it is ok
    ops.setOp(0,new Operator(oper_),true);
    complex_t valk=contractor.contract(state,ops,state);
    out<<real(valk)/norm<<"\t";
    ops.setOp(0,idCenter,false);
  }
  for(int pos=1;pos<len;pos++){
    ops.setOp(pos,&rightOp,false);
    // Compute and write result
    complex_t valk=contractor.contract(state,ops,state);
    out<<real(valk)/norm<<"\t";
    ops.setOp(pos,&identAll,false);
  }
  out<<endl;
  // delete trash
  delete idCenter;
}
