/**
   \file lightcone.cpp
   Find the time evolution of a local excitation, starting in the
   ground state of some Hamiltonian, and using a light-cone + folding 
   approach.
   
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
			int d,int Dl);

/** 
    lightcone tries to approximate the evolution after a local
    perturbation using a "light-cone" approach.
    Arguments:
    \param <L> (int) length of the chain used to find the ground state
    \param <Lb> (int) initial number of sites to the sides of the
    center used in the evolution
    \param <J> (double) Parameter of the Ising Hamiltonian
    \param <g> (double) Parameter of the Ising Hamiltonian
    \param <h> (double) Parameter of the Ising Hamiltonian
    \param <D2> (int) maximum bond dimension in the (foldd) evolved state
    \param <D0> (int) maximum bond dimension in the (foldd) evolved state
    \param <delta> (double) step width
    \param <M> (int) number of time steps
    \param <outfname> (char*) name of the output file for the magnetization
    \param <app> (int) whether the output file is to be kept (app==1)
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

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d=2;
  //  int D0=D2/2; // D for GS
  int Dl=D0; // D of the single term

  int len=1+2*Lb+2;
  MPS center(len,D0,d);
  MPS center2(len+2,D0,d); // to hold the extension
  // and place for the RDM of the rest of chain
  mwArray rhoL,rhoR;

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
    gs.gaugeCond('R');
    cout<<"Ground state found with eigenvalue "<<E0<<endl;
    *out<<"% "<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"
	<<D0<<"\t"<<E0<<"\t";
    *out<<endl;
    
    // Now the tricks start
    // First, get the central site
    int L0=(L-1)/2;
    for(int k=1;k<len-1;k++)
      center.setA(k,gs.getA(L0).getA());
    for(int k=1;k<len+1;k++)
      center2.setA(k,gs.getA(L0).getA());
    // The trick now is to use an extra site on both edges, with sqrt(rho)    
    Dl=gs.getA(L0).getDl();
    //I need the square root and the inverse of rho
    mwArray sqrhoL,sqrhoR,sqrhoLinv,sqrhoRinv;
    {
      // Maybe better approach is to contract the environment we have
      // in the long system (in the TI case should be replaced by true
      // EV). A bit stupid, here, since one of them (rhoL) will be 1
      //      rhoL=contractor.contract(gs,gs,L0,'L');
      rhoL=identityMatrix(Dl);
      sqrhoL=identityMatrix(Dl);
      sqrhoLinv=identityMatrix(Dl);
      //cout<<"Computed rhoL="<<rhoL<<endl;
      rhoR=contractor.contract(gs,gs,L0,'R');
      //cout<<"Computed rhoR="<<rhoR<<endl;
      vector<complex_t> vals; 
      mwArray U;
      wrapper::eig(rhoR,vals,U,true);
      mwArray aux=diag(vals);
      aux=sqrt(aux);
      //sqrhoL=U*aux*Hconjugate(U);
      wrapper::eig(rhoR,vals,U,true);
      aux=diag(vals);
      aux=sqrt(aux);
      int nr=0;
      mwArray auxInv=inverse(aux,nr);
      // the one multiplying on the right
      sqrhoR=conjugate(U)*aux;
      //sqrhoLinv=inverse(sqrhoL,nr);
      U.transpose();
      sqrhoRinv=auxInv*U;
    }
    // With this I compute the edge tensors
    mwArray Al=sqrhoL;Al.reshape(Indices(Dl,1,Dl));
    mwArray Ar=sqrhoR;Ar.reshape(Indices(Dl,Dl,1));
    center.replaceSite(0,Al);
    center2.replaceSite(0,Al);
    //cout<<"Set tensor 0 for center and center2"<<endl;
    center.replaceSite(len-1,Ar);
    center2.replaceSite(len+1,Ar);
    //cout<<"Set tensor "<<len-1<<"("<<len+1<<") for center(2)"<<endl;
    // But the 1 for center2 has to have sqrtinv!!!
    //Al=gs.getA(L0).getA();Ar=Al;
    //Indices dims=Al.getDimensions();
    //Al.reshape(Indices(dims[0]*dims[1],dims[2]));
    //Al.multiplyRight(sqrhoLinv);
    //Al.reshape(Indices(dims[0],dims[1],dims[2]));
    //center2.replaceSite(1,Al);
    Ar=gs.getA(L0).getA();
    Indices dims=Ar.getDimensions();
    Ar.permute(Indices(2,1,3));
    Ar.reshape(Indices(dims[1],dims[0]*dims[2]));
    Ar.multiplyLeft(sqrhoRinv);
    Ar.reshape(Indices(dims[1],dims[0],dims[2]));
    Ar.permute(Indices(2,1,3));
    center2.replaceSite(len,Ar);
    // Now I forget about the Hamiltonian, and everything

  }
  //  cout<<"Constructed center MPS: "<<center<<endl;
  // cout<<"Constructed center2 MPS: "<<center2<<endl;

  int middC=Lb+1; // center of the "short" chain
  // Apply the local perturbation (now, strictly local)
  complex_t sigmaX_[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),sigmaX_);
  // Since I apply sigmax, I shouldn't worry about the norm (only about
  // the gauge condition, but this will be fixed after folding)
  //center.applyLocalOperator(Lb,sigX,false);
  // And fold the stuff
  MPS foldedS;
  center.fold(foldedS);
  cout<<"Length of center="<<center.getLength()<<", length of foldedS="
      <<foldedS.getLength()<<endl;
  MPS foldedS2;
  center2.fold(foldedS2);
  const mwArray& extraSite1=foldedS2.getA(foldedS2.getLength()-2).getA();
  const mwArray& extraSite2=foldedS2.getA(foldedS2.getLength()-1).getA();
  center2.clear();

  foldedS.applyLocalOperator(0,sigX,false);

  // I need an Identity operator on the edges
  mwArray iDl=identityMatrix(Dl);
  iDl.reshape(Indices(Dl,1,Dl,1));
  FoldedOperator idEdge(iDl,iDl);

  //  cout<<"Constructed perturbed foldedS"<<endl;
  int count=0; // nr of step
  double t=0.;// H is constant
  int oT=2; // Trotter order
  bool imag=false; // no imaginary time evolution
  while(count<M){
    // Get the Hamiltonian (edges are just the bond: do not count
    IsingHamiltonian shortH(len-2,d,J,g,h);
    // and the evolution MPOs
    MPO evol(len-2);
    shortH.getUMPO(evol,delta,t,imag,oT);
    double time=count*(delta+1);
    cout<<"Time evolution step nr "<<count+1<<", time="<<delta*(count+1)<<endl;
    cout<<"len="<<len<<", MPs length="<<foldedS.getLength()<<", MPO length"
	<<evol.getLength()<<endl;
    // fold
    MPO fEvol0((len-1)/2);
    evol.fold(fEvol0);
    // Now I do the trick and enlarge it with Identity on the edges
    MPO fEvol(fEvol0.getLength()+1);
    for(int k=0;k<fEvol0.getLength();k++)
      fEvol.setOp(k,&fEvol0.getOp(k),false);
    fEvol.setOp(fEvol0.getLength(),&idEdge,false);
    // Apply (folded) evolution
    MPS aux(foldedS); // temporary copy
    contractor.optimize(fEvol,aux,foldedS,Dl);
    // Compute something
    computeExpectation(foldedS,sigX,*out,d,Dl);
    // Enlarge the state
    MPS aux2(foldedS); //(foldedS.getLength()+1,Dl,d);
    // the extra tensor
    extendFoldedState(foldedS,aux2,extraSite1,extraSite2,d,Dl);
    foldedS=aux2;
    len+=2;
    count++;
  }


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

void computeExpectation(MPS& state,const mwArray& oper,ofstream& out,int d,
			int Dl){
  int len=state.getLength(); // of the folded one!
  MPO ops(len);
  mwArray id_=identityMatrix(d);
  mwArray oper_(oper); oper_.reshape(Indices(d,1,d,1));
  id_.reshape(Indices(d,1,d,1));
  FoldedOperator identAll(id_,id_);
  for(int k=1;k<len-1;k++) ops.setOp(k,&identAll,false);
  // Center and edges are special!
  mwArray idEd_=identityMatrix(Dl);idEd_.reshape(Indices(Dl,1,Dl,1));
  FoldedOperator identEdge(idEd_,idEd_);
  ops.setOp(len-1,&identEdge,false);
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
  for(int pos=len-2;pos>0;pos--){
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
  for(int pos=1;pos<len-1;pos++){
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
