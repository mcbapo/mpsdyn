/**
   \file hypIsing.cpp
   Find ground staate for a hyperbolically deformed Ising model.
   
   \author Mari-Carmen Banuls
   \date 30/09/2011

*/

#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "HyperbolicIsingHamiltonian.h"

using namespace shrt;

/** hypIsing runs the findGroundState routine with the MPO for the
    hyperbolically deformed Ising Hamiltonian 
    \ref <HyperbolicIsingHamiltonian>, to check convergence and
    energies for different values of lambda.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <lambda> (double) parameter \f$\lambda\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <outfnameCorr> (char*) name of the output file for correlations 
                      (always overwritten)
    \param <outfname> (char*) name of the output file for teh energy
    \param <app> (int) whether the output file is to be kept (app==1)

*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double J=atof(argv[++cntr]);
  double g=atof(argv[++cntr]);
  double h=atof(argv[++cntr]);
  double lambda=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfnameC=argv[++cntr];
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  cout<<"Initialized arguments: L="<<L
      <<", J="<<J
      <<", g="<<g
      <<", h="<<h
      <<", lambda="<<lambda
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"# N\t J\t g\t h\t lambda\t D\t Energy"<<endl;
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
   cout<<"Initialized random state "<<endl;
   
   // First: put H and find the GS
   int Dop=3; // bond dimension of Ising H
   double E0=0.;
   HyperbolicIsingHamiltonian hypI(L,J,g,h,lambda);
   cout<<"Created the Hamiltonian"<<endl;
   const MPO& hamil=hypI.getHMPO();
   cout<<"Constructed the hamil MPO"<<endl;
   contractor.findGroundState(hamil,D,&E0,gs);
   cout<<"Ground state found with eigenvalue "<<E0<<endl;

   *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<lambda<<"\t"
       <<D<<"\t"<<E0<<"\t";
   *out<<endl;
   out->close();
   // Now measure all terms
   double* valZZ=new double[L-1];
   double* valZ=new double[L];
   double* valX=new double[L];
   mwArray identPhys=identityMatrix(d);
   identPhys.reshape(Indices(d,1,d,1));
   Operator* idOp=new Operator(identPhys);
   complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
   mwArray sigZ(Indices(d,d),dataz);
   sigZ.reshape(Indices(d,1,d,1));
   Operator* opZ=new Operator(sigZ);
   complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
   mwArray sigX(Indices(d,d),datax);
   sigX.reshape(Indices(d,1,d,1));
   Operator* opX=new Operator(sigX);
   MPO idMPO(L);
   for(int k=0;k<L;k++)
     idMPO.setOp(k,idOp);
   // Compute ZZ terms
   for(int k=0;k<L;k++){
     idMPO.setOp(k,opZ);
     if(k<L-1){
       idMPO.setOp(k+1,opZ);
       valZZ[k]=real(contractor.contract(gs,idMPO,gs));     
       idMPO.setOp(k+1,idOp);
     }
     valZ[k]=real(contractor.contract(gs,idMPO,gs));     
     idMPO.setOp(k,opX);
     valX[k]=real(contractor.contract(gs,idMPO,gs));     
     idMPO.setOp(k,idOp);
   }
   // Now print them all??
   // 
  ofstream* out2=new ofstream(outfnameC);
  if(!out2->is_open()){
    cout<<"Error: impossible to open file "<<outfnameC<<
      " for output"<<endl;
    exit(1);
  }
  *out2<<"% N="<<L<<", J="<<J<<", g="<<g<<", h="<<h<<", lambda="<<lambda
       <<", D="<<D<<", E0="<<E0<<endl;
  *out2<<"% N\t <ZZ>\t <X>(N)\t <Z>(N)\t <h_12(sym)>\t <h_12(right)>\t coefZZ\t coefX\t coefZ\t entropy"<<endl;
  for(int k=0;k<L-1;k++){
    *out2<<k<<"\t"<<valZZ[k]<<"\t"<<valX[k]<<"\t"<<valZ[k]<<"\t";
    double energyLink=-hypI.getCoefZZ(k)*valZZ[k];
    double energyLinkR=energyLink; // el de la derecha
    double factL=k==0?1.:.5;
    double factR=k==L-1?1.:.5;
    energyLink+=-factL*hypI.getCoefX(k)*valX[k];
    energyLink+=-factR*hypI.getCoefX(k+1)*valX[k+1];
    energyLink+=-factL*hypI.getCoefZ(k)*valZ[k];
    energyLink+=-factR*hypI.getCoefZ(k+1)*valZ[k+1];
    energyLinkR+=-hypI.getCoefX(k+1)*valX[k+1];
    energyLinkR+=-hypI.getCoefZ(k+1)*valZ[k+1];
    if(k==0){
      energyLinkR+=-hypI.getCoefX(k)*valX[k];
      energyLinkR+=-hypI.getCoefZ(k)*valZ[k];
    }
    *out2<<energyLink<<"\t"<<energyLinkR<<"\t";
    // To end, I copy also the values of the coeffs
    *out2<<hypI.getCoefZZ(k)<<"\t"<<hypI.getCoefX(k)<<"\t"
	 <<hypI.getCoefZ(k)<<"\t";
    *out2<<contractor.getEntropy(gs,k);
    *out2<<endl;
  }
  // And for the lastlink, write the corresponding values (NaN for not def)
  *out2<<L-1<<"\t NaN"<<"\t"<<valX[L-1]<<"\t"<<valZ[L-1]<<"\t";
  *out2<<"NaN \t Nan\t NaN\t"<<hypI.getCoefX(L-1)<<"\t"
       <<hypI.getCoefZ(L-1)<<"\t";
  *out2<<"NaN";
  *out2<<endl;
  out2->close();

  // Now compute the evolution MPO
  MPO expU(L);double delta=0.1;
  cout<<"About to construct the MPO for the exponential"<<endl;
  hypI.getUMPO(expU,delta,2,false);
  expU.exportMPOtext("testMPOU.txt");

  cout<<"Saved? MPO"<<endl;
}
