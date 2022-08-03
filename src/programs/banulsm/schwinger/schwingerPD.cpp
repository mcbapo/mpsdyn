
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "SchwingerHamiltonian.h"


#ifdef MACOSX
//#define LOCALDIR "/Users/banulsm/Desktop/SVN/projects/mpsdyn/src/bin/"
#define LOCALDIR "/Users/banuls/SVN/peps/projects/mpsdyn/src/bin/"
#else
#define LOCALDIR "/home/banulsm/ALaPorra/bin/"
#endif
#define MAXLEN 120

/** SchwingerPD runs the findGroundState routine with the MPO for the
    Schwinger Hamiltonian \ref <SchwingerHamiltonian>, and saves
    values of order parameters, to construct a phase diagram. Since
    the order parameters \f$\Gamma^{\alpha}\f$ and \f$\Gamma^5\f$ have
    to be computed as expectation values between two ground states,
    for sectors \f$L_0+\alpha=-1/2\f$ and $L_o+\alpha=1/2$, it computes the 
    ground states for \f$L_0=0\f$ and \f$L_0=-1\f$ with \f$\alpha=1/2\f$.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <mg> (double) parameter \f$m/g\f$ of the Hamiltonian
    \param <x> (double) parameter \f$x\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)

*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double mg=atof(argv[++cntr]);
  double x=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  char filename[MAXLEN];

  double mu=2*mg*sqrt(x);

  double alpha=.5;
  double L0P=0.;
  double L0M=-1.;

  cout<<"Initialized arguments: L="<<L
      <<", mu="<<mu
      <<", x="<<x
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;


  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"# N\t x\t mg\t D\t Energy(+1/2)\t Energy(-1/2)\t w0(+)/(2Nx)  \t "
	<<"w0(-)/(2Nx)\t overlap\t GammaA(+)-L0-alpha\t GammaA(-) \t GammaA(+-)\t"
	<<"(L/sqrt(x))Gamma5(+)\t Gamma5(-)\t Gamma5(+-)"<<endl;
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }

   Contractor& contractor=Contractor::theContractor();
   cout<<"Initialized Contractor?"<<endl;
   int d=2;
   MPS gsP(L,D,d);
   gsP.setRandomState(); // the intial state, random
   gsP.gaugeCond('R',1);
   gsP.gaugeCond('L',1);
   MPS gsM(gsP);
   
   // First: put H and find the GS
   int Dop=5; // bond dimension of Schwinger H
   double lambdaP=0.;
   double lambdaM=0.;
   SchwingerHamiltonian hSchP(L,mu,x,L0P,alpha);
   SchwingerHamiltonian hSchM(L,mu,x,L0M,alpha);
   const MPO& hamilP=hSchP.getHMPO();
   contractor.findGroundState(hamilP,D,&lambdaP,gsP);
   cout<<"Ground state (sector +1/2) found with eigenvalue "<<lambdaP<<endl;
   strcpy(filename,LOCALDIR);
   gsP.exportMPS(strcat(filename,"mpsGS.mat"));

   const MPO& hamilM=hSchM.getHMPO();
   contractor.findGroundState(hamilM,D,&lambdaM,gsM);
   cout<<"Ground state (sector -1/2) found with eigenvalue "<<lambdaM<<endl;

   MPO GammaA(L);
   MPO Gamma5(L);
   hSchP.constructGammaAlphaMPO(GammaA);
   hSchP.constructGamma5MPO(Gamma5);

   // Now compute and save all the expectation values
   complex_t gammaAlpha=contractor.contract(gsP,GammaA,gsM);
   complex_t gammaAlphaP=contractor.contract(gsP,GammaA,gsP);
   complex_t gammaAlphaM=contractor.contract(gsM,GammaA,gsM);
   complex_t gammaAlphaCt=contractor.contract(gsP,gsM);
   complex_t gamma5=contractor.contract(gsP,Gamma5,gsM);
   complex_t gamma5P=contractor.contract(gsP,Gamma5,gsP);
   complex_t gamma5M=contractor.contract(gsM,Gamma5,gsM);

   *out<<L<<"\t"<<x<<"\t"<<mg<<"\t"
       <<D<<"\t"<<lambdaP<<"\t"<<lambdaM<<"\t"
       <<lambdaP/(2*L*x)<<"\t"
       <<lambdaM/(2*L*x)<<"\t";
   // overlap
   *out<<gammaAlphaCt<<"\t";
   // *out<<gammaAlphaCt.re<<"\t";
   // *out<<gammaAlphaCt.im<<"\t";
   //Gamma alpha (+)
   *out<<gammaAlphaP.re<<"\t";
   // *out<<gammaAlphaP.im<<"\t";
   //Gamma alpha (-)
   *out<<gammaAlphaM.re<<"\t";
   // *out<<gammaAlphaM.im<<"\t";
   //Gamma alpha (+-)
   *out<<.5*(gammaAlphaP-gammaAlphaM-gammaAlpha+contractor.contract(gsM,GammaA,gsP))<<"\t";
   //   *out<<gammaAlpha.re<<"\t";
   // *out<<gammaAlpha.im<<"\t";
   //Gamma5(+)
   *out<<gamma5P.re<<"\t";
   // *out<<gamma5P.im<<"\t";
   //Gamma5(-)
   *out<<gamma5M.re<<"\t";
   //*out<<gamma5M.im<<"\t";
   //Gamma5(+-)
   *out<<.5*(gamma5P-gamma5M-gamma5+contractor.contract(gsM,Gamma5,gsP))<<"\t";
   // *out<<gamma5.re<<"\t";
   // *out<<gamma5.im<<endl;
   out->close();
}
