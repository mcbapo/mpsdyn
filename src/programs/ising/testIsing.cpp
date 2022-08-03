
#include <math.h>
#include <iomanip>

#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

/** \file testIsing.cpp
 
    testIsing runs the findGroundState routine with the MPO for the
    Ising Hamiltonian \ref <IsingHamiltonian.cpp>, to check convergence
    and then runs imaginary time evolution to check also the
    implementation of the unitary evolution.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <M> (int) maximum number of imaginary time steps 
           (if 0, no imaginary time phase is run)
    \param <delta> (double) step width for imaginary time evolution
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
  int D=atoi(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  cout<<"Initialized arguments: L="<<L
      <<", J="<<J
      <<", g="<<g
      <<", h="<<h
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% N\t J\t g\t h\t D\t Energy\t Energy/N"<<endl;
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  cout<<setprecision(10);

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
   //hamil.exportForMatlab("isingHmpo.m");
   cout<<"Initial value, with initial state"<<endl;
   cout<<contractor.contract(gs,hamil,gs)<<endl;

   contractor.findGroundState(hamil,D,&E0,gs);
   cout<<"Ground state found with eigenvalue "<<E0
       <<" (Energy per particle="<<E0/L<<")"<<endl;
   // mwArray Amid=gs.getA(L/2).getA();
   // ofstream Amid_("Amid.dat");
   // Amid.savetext(Amid_);
   // Amid_.close();
   // //   putForMatlab(cout,Amid,"Amid");
   // cout<<endl;
   // gs.exportForMatlab("gsIsingNoG.m");

   *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"
       <<D<<"\t"<<E0<<"\t"<<E0/L<<"\t";
   *out<<endl;

   if(M>0){
     // Now try evolution with the exponential, instead
     MPO evol(L);
     double t=0.;
     int oT=2; // Trotter order
     bool imag=true; 
     hamI.getUMPO(evol,delta,t,imag,oT);
     *out<<endl;
     *out<<"% Imaginary time evolution!! "<<endl;
     *out<<"% N\t J\t g\t h\t t\t Energy\t Energy/N"<<endl;

     MPS imaGS(L,D,d);
     double oldE=1000;
     double time=0.;
     int cnt=0;
     bool conver=0;
     while(cnt<M&&!conver){
       time+=delta;
       cout<<"Time evolution step nr "<<cnt+1<<", time="<<time
	   <<", E0="<<E0<<endl;
       MPS aux(imaGS); // temporary copy
       contractor.optimize(evol,aux,imaGS,D);
       
       imaGS.gaugeCond('R',true);
       E0=real(contractor.contract(imaGS,hamil,imaGS));
       
       *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"
	   <<time<<"\t"<<E0<<"\t"<<E0/L<<"\t";
       *out<<endl;
       if(cnt>0&&abs((E0-oldE)/oldE)<1E-8) conver=true;
       cnt++;
       oldE=E0;
     }

   //  evol.exportForMatlab("isingUmpo.m");
   // hamil.exportForMatlab("isingHmpo.m");
   
   cout << "Overlap: " << abs(contractor.contract(imaGS,gs)) << endl;
   }

   out->close();
}
