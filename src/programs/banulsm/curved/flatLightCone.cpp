/**
   \file flatLightCone.cpp
   Evolve Ising ground state after a local perturbation, using only a
   few sites of the Hamiltonian
   
   \author Mari-Carmen Banuls
   \date 11/05/2012

*/

#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

void applyPerturbation(const MPS& current,MPS& result,int center);
void computeExpectation(MPS& state,const Operator& oper,ofstream& out);

/** flatLightCone runs the findGroundState routine with the MPO for the
    Ising Hamiltonian \ref <IsingHamiltonian>, then applies a local
    perturbation and evolves in real time for M steps using width delta.
    In the output file, we write the magnetization site by site 
    as a function of time.
    In the second file, entropies are written across every cut. 

    Receives arguments:
    \param <L> (int) length of the chain
    \param <Lb> (int) how many sites around the center, from start
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <M> (int) number of time steps
    \param <stepBlock> (int) frequency (number of time steps) of writing
    \param <delta> (double) step width
    \param <outfname> (char*) name of the output file for the energy
    \param <app> (int) whether the output file is to be kept (app==1)
    \param <outfname2> (char*) name of the output file for the energy
    \param <app2> (int) whether the output file is to be kept (app==1)
*/


int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  int Lb=atoi(argv[++cntr]);
  double J=atof(argv[++cntr]);
  double g=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int stpBlk=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);
  const char* outfname2=argv[++cntr];
  int app2=atoi(argv[++cntr]);

  double h=0.;
  if(L%2==0)L++;
  int center=L%2==0?L/2:(L-1)/2;

  cout<<"Initialized arguments: L="<<L
      <<", J="<<J
      <<", g="<<g
      <<", outfile="<<outfname
      <<", app="<<app
      <<", M="<<M
      <<", delta="<<delta
      <<", D="<<D<<endl;

  ofstream* out,*out2;
  if(!app){
    out=new ofstream(outfname);
    *out<<"# N\t J\t g\t h\t D\t Energy"<<endl;
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
   MPS gs(L,D,d);
   gs.setRandomState(); // the intial state, random
   gs.gaugeCond('R',1);
   gs.gaugeCond('L',1);
   cout<<"Initialized random state, norm "<<contractor.contract(gs,gs)<<endl;
   
   // First: put H and find the GS
   int Dop=3; // bond dimension of Ising H
   double E0=0.;
   IsingHamiltonian hamI(L,d,J,g,h);
   cout<<"Created the Hamiltonian"<<endl;
   const MPO& hamil=hamI.getHMPO(0.);
   cout<<"Constructed the hamil MPO"<<endl;
   //   hamil.exportMPOtext("flatHmpo.txt");
   cout<<"Initial value, with initial state"<<endl;
   cout<<contractor.contract(gs,hamil,gs)<<endl;

   contractor.findGroundState(hamil,D,&E0,gs);
   cout<<"Ground state found with eigenvalue "<<E0<<endl;
   *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"
       <<D<<"\t"<<E0<<"\t";
   *out<<endl;

   // ///// DEBUG
   // gs.exportMPS("gsflat.dat");

   // Now apply perturbation (and normalize)
   // Create an MPO for it
   MPS state(gs);
   complex_t sigmaX_[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
   mwArray sigX(Indices(d,1,d,1),sigmaX_);
   Operator sigmaX(sigX);
   applyPerturbation(gs,state,center);

   // Normalize
   state.gaugeCond('R',true);

   // // DEBUG
   // state.exportMPS("perflat.dat");

   *out<<"% Evolution in time, with "<<M<<" steps of width "<<delta<<endl;
   *out<<"% M\t t\t <sigx(pos)>"<<endl;
   out->close();

  if(!app2){
    out2=new ofstream(outfname2);
   *out2<<"% Evolution in time, with "<<M<<" steps of width "<<delta<<endl;
   *out2<<"% M\t t\t E \t <S(pos)>"<<endl;
  }
  else
    out2=new ofstream(outfname2,ios::app);
  if(!out2->is_open()){
    cout<<"Error: impossible to open file "<<outfname2<<
      " for output"<<endl;
    exit(1);
  }


   // Now run evolution with the exponential
   MPO evol(L);
   double t=0.;// H is constant
   int oT=2; // Trotter order
   bool imag=false; // no imaginary time evolution
   hamI.getUMPO(evol,delta,t,imag,oT);
   //   evol.exportMPOtext("operFlat.mat");
   int cnt=0;
   //int Lb=1; // oper applied from center-Lb till center+Lb
   MPO auxMPO(L);
   mwArray id=identityMatrix(d);id.reshape(Indices(d,1,d,1));
   Operator idOp(id);
   for(int k=0;k<L;k++)
     auxMPO.setOp(k,&idOp,false);
   auxMPO.setOp(center,&evol.getOp(center),false);
   for(int k=center-Lb+1;k<center+Lb;k++)
     auxMPO.setOp(k,&evol.getOp(k),false);
   while(cnt<M){
     double time=cnt*(delta+1);
     cout<<"Time evolution step nr "<<cnt+1<<", time="<<delta*(cnt+1)<<endl;
     // ///// DEBUG Save it now
     //   if(cnt==2)   state.exportMPS("evolflat.dat");

     auxMPO.setOp(center-Lb,&evol.getOp(0),false);
     auxMPO.setOp(center-Lb+1,&evol.getOp(center),false);
     auxMPO.setOp(center+Lb,&evol.getOp(evol.getLength()-1),false);
     auxMPO.setOp(center+Lb-1,&evol.getOp(center),false);

     MPS aux(state); // temporary copy
     contractor.optimize(auxMPO,aux,state,D);
     if(cnt%stpBlk==0||cnt==M||cnt==0){
       out->open(outfname,ios::app);
       *out<<cnt<<"\t"<<time<<"\t";
       // normalize?
       computeExpectation(state,sigmaX,*out);
       *out<<endl;
       out->close();
       state.gaugeCond('R',1); // normalize
       double energy=real(contractor.contract(state,hamil,state));
       out2->open(outfname2,ios::app);
       *out2<<cnt<<"\t"<<time<<"\t"<<energy<<"\t";
       // Now I would also like to compute the entropy between the right 
       // part and the left one across every cut in the chain!
       for(int pl=1;pl<=L-1;pl++){
	 double entr=contractor.getEntropy(state,pl);
	 *out2<<entr<<"\t";
       }
       *out2<<endl;
       out2->close();
     }
     cnt++;
     if(Lb<center) Lb++;
   }
   if(out->is_open())
     out->close();
   delete out;
}

void applyPerturbation(const MPS& current,MPS& result,int center){
  cout<<"Inside applyPerturbation"<<endl;
  int L=current.getLength();
  int d=current.getA(0).getd();
  // the MPO
  MPO perturbation(L);
  mwArray ident=identityMatrix(d);ident.reshape(Indices(d,1,d,1));
  Operator id(ident);
  for(int k=0;k<L;k++)
    if(k!=center)
      perturbation.setOp(k,&id,false);
  // The operator in the middle I do by hand, using a trick
  //IsingHamiltonian auxH(2,2,0.,.5,0.);
  //const MPO& hamil=auxH.getHMPO(0.);
  //perturbation.setOp(center,&(hamil.getOp(0)),false);
  //perturbation.setOp(center+1,&(hamil.getOp(1)),false);
   complex_t sigmaX_[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
   mwArray sigX(Indices(d,1,d,1),sigmaX_);
   Operator sigmaX(sigX);
   perturbation.setOp(center,&sigmaX,false);
  Contractor& contractor=Contractor::theContractor();
  // // The operators on sites center and center+1 by hand
  //  mwArray delta0(Indices(2,1),sigmaX_);  
  //  mwArray delta1(Indices(2,1),&sigmaX_[2]);
  //  sigX.reshape(Indices(1,d*d));ident.reshape(Indices(1,d*d));
  //  mwArray pertL=.5*delta0*sigX+delta1*ident;
  //  mwArray pertR=.5*delta1*sigX+delta0*ident;
  //  pertL.reshape(Indices(2,d,1,d));pertL.permute(Indices(2,3,4,1));
  //  pertR.reshape(Indices(2,d,d,1));pertR.permute(Indices(2,1,3,4));
  //  perturbation.setOp(center,new Operator(pertL),true);
  //  perturbation.setOp(center+1,new Operator(pertL),true);
  contractor.optimize(perturbation,current,result);
   
  cout<<"Leaving applyPerturbation"<<endl;

}

void computeExpectation(MPS& state,const Operator& oper,ofstream& out){
  int len=state.getLength();
  int d=2;
  static MPO ops(len);
  static bool init(false);
  static Operator ident(Indices(d,1,d,1));
  if(!init){
    mwArray id_=identityMatrix(d);
    id_.reshape(Indices(d,1,d,1));
    ident.setData(&id_);
    // A trick to store copies of the mwArrays in the opers
    ident.permuteOp(Indices(1,2,3,4));
    init=true;
    for(int k=0;k<len;k++) ops.setOp(k,&ident,false);
  }
  Contractor& contractor=Contractor::theContractor();
  double norm=real(contractor.contract(state,state));
  for(int pos=0;pos<len;pos++){
    ops.setOp(pos,&oper,false);
    // Compute and write result
    complex_t valk=contractor.contract(state,ops,state);
    out<<real(valk)/norm<<"\t";
    ops.setOp(pos,&ident,false);
  }
}
