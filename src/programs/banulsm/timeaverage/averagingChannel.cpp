#include <math.h>
#include <iomanip>

#include "TIMPS.h"
#include "Contractor.h"
#include "VariableOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

/** averagingChannel uses the same operator of averageIsing to
    represent the time average over some number of evolution steps.
    In this case, the Ising evolution is applied to a TIMPS, by just
    truncating the evolved TIMPS after each application of the
    channel.  The procedure is repeated for a maximum number of
    applications, and as observable we record the local
    magnetization. 

    Receives arguments:
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <delta> (double) width of the time step used
    \param <M> (int) number of steps in each channel application
    \param <repe> (int) maximum number of iterations of the channel application
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the energy
    \param <app> (int) whether the output file is to be kept (app==1)


*/

int main(int argc,const char* argv[]){
  // Read input arguments

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  double J=atof(argv[++cntr]);
  double g=atof(argv[++cntr]);
  double h=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int repe=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  cout<<"averagingChannel -> Initialized arguments "
      <<", J="<<J
      <<", g="<<g
      <<", h="<<h
      <<", delta="<<delta
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% J="<<J<<", g="<<g<<", h="<<h
	<<" D="<<D<<endl;
    *out<<"% M="<<M<<", delta="<<delta<<endl;
    *out<<"% k\t delta\t t/2 \t D\t avrg sigx"<<endl;
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

   // First: create the Hamiltonian (only three sites, as I am interested 
   // in the middle, TI one)
   int Dop=2; // bond dimension of Ising exponential mpo
   IsingHamiltonian hamI(3,d,J,g,h);
   cout<<"Created the Hamiltonian"<<endl;
   mwArray basicU=hamI.getUOperator(delta,0.,false,2);
   // This is the du+Dl+dd+Dr basic piece
   // I want to increase it with an identity
   // and combine with the adjoints for the ancillary side
   mwArray identOp=identityMatrix(d*d); // double identity
   mwArray basicUad(basicU);
   basicUad.conjugate(); // the adjoint for the ancillas

   mwArray doubleU(basicU);
   doubleU.reshape(Indices(-1,1));
   basicUad.reshape(Indices(1,-1));
   doubleU.multiplyRight(basicUad);
   doubleU.reshape(Indices(d,Dop,d,Dop,d,Dop,d,Dop));
   doubleU.permute(Indices(1,5,2,6,3,7,4,8));
   doubleU.reshape(Indices(d*d,Dop*Dop,d*d,Dop*Dop));

   // The coefficients of the middle "mpo"
   mwArray C(Indices(2,2,2));
   C.setElement(ONE_c,Indices(0,0,0));
   C.setElement(ONE_c,Indices(0,1,1));
   C.setElement(ONE_c,Indices(1,1,1));
   //   C.setElement(ONE_c,Indices(1,0,0));

   // And for the initial step (0)
   mwArray Cl(Indices(2,2));
   Cl.setElement(ONE_c,Indices(0,0));
   Cl.setElement(ONE_c,Indices(1,1));

   // The "operator basis"
   mwArray Z(Indices(2,d*d,Dop*Dop+1,d*d,Dop*Dop+1));
   for(int i=0;i<d*d;i++)
     for(int j=0;j<d*d;j++){
       Z.setElement(identOp.getElement(Indices(i,j)),
		    Indices(0,i,0,j,0));
       for(int l=0;l<Dop*Dop;l++)
	 for(int r=0;r<Dop*Dop;r++){
	   Z.setElement(doubleU.getElement(Indices(i,l,j,r)),
			Indices(1,i,l+1,j,r+1));
	 }
     }

   // Now construct the Operator itselt
   C.reshape(Indices(-1,2));
   Z.reshape(Indices(2,-1));
   mwArray data=C*Z;
   data.reshape(Indices(2,2,d*d,Dop*Dop+1,d*d,Dop*Dop+1));
   data.permute(Indices(2,3,4,1,5,6));
   data.reshape(Indices(2*d*d,Dop*Dop+1,2*d*d,Dop*Dop+1));
   cout<<"Created operator data "<<data.getDimensions()<<endl;
   Operator* oper=new Operator(data);

   mwArray datal=Cl*Z;
   datal.reshape(Indices(1,2,d*d,Dop*Dop+1,d*d,Dop*Dop+1));
   datal.permute(Indices(2,3,4,1,5,6));
   datal.reshape(Indices(2*d*d,Dop*Dop+1,d*d,Dop*Dop+1));
   cout<<"Created operator datal "<<datal.getDimensions()<<endl;
   Operator* operl=new Operator(datal);

   // Initial state (X+)
   mwArray stateX(Indices(2,1));
   stateX.setElement(1./sqrt(2)*ONE_c,Indices(0,0));
   stateX.setElement(1./sqrt(2)*ONE_c,Indices(1,0));
   mwArray adjSt(Hconjugate(stateX));
   mwArray dataX=stateX*adjSt; // the density matrix
   TIMPS rhoX;
   rhoX.setA(reshape(dataX(d*d,1,1)));

   // Observables
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigz(Indices(d,1,d,1),dataz);
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigx(Indices(d,1,d,1),datax);

  Operator sZ(sigz);
  Operator sX(sigx);

   int cntR=0;
   while(cntR<repe){
     cntR++;
     cout<<"Starting application nr. "<<cntR<<" of the channel "<<endl;
     // first application

     rhoX.evolve(*operl,D,false,false);
     for(int k=0;k<M-1;k++)
       rhoX.evolve(*oper,D,false,false);
     

     //TODO!!!! Me falta calcular valores esperados en RHO!!!!!
     double Ex=rhoX.expectationValueSingle(sX);
     double Ez=rhoZ.expectationValueSingle(sZ);



   }

}
