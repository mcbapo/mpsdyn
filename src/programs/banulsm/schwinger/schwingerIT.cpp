
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "SchwingerHamiltonian.h"

#ifdef TESTMULTI
#define BLOCK 4
#endif

#ifdef MACOSX
#define LOCALDIR "/Users/banulsm/Desktop/SVN/projects/mpsdyn/src/bin/"
#else
#define LOCALDIR "/home/banulsm/SVN/peps/projects/folding/src/bin/"
#endif

/** SchwingerIT tries to find the ground state of the Schwinger Hamiltonian
    \ref <SchwingerHamiltonian> using imaginary time evolution.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <mg> (double) parameter \f$m/g\f$ of the Hamiltonian
    \param <x> (double) parameter \f$x\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <infname> (char*) name of the input file containing the 
    time configuration
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
  const char* infname=argv[++cntr];
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  char filename[120];

  double mu=2*mg*sqrt(x);

  double alpha=.5;
  double L0P=0.;
  //double L0M=-1.;

  cout<<"Initialized arguments: L="<<L
      <<", mu="<<mu
      <<", x="<<x
      <<", infile="<<infname
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"# N\t x\t mg\t D\t delta\t M\t Energy(+1/2)\t w0(+)/(2Nx)"
	<<endl;
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);
  // Read step widths
  int M;
  ifstream inlog(infname);
  double t1=0;
  inlog>>M; // true final number of steps
  double* deltas=new double[M];
  // Initialize the step widths
  double tmp;int cnt;double normFact;
  int k=0; // actual step
  while(k<M-1){ // read all the file
    inlog>>cnt;
    inlog>>tmp;
    for(int l=0;l<cnt;l++){
      deltas[k++]=tmp;
      //cout<<"Step "<<k<<" is "<<tmp<<endl;
    }
  }
  if(inlog.good())
    inlog>>normFact;
  else
    normFact=1.;
  inlog.close();

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
   int d=2;
   MPS gsP(L,D,d);
   gsP.setRandomState(); // the initial state, random -> not to the GS!!
   //gsP.setProductState(p_xplus);
   gsP.gaugeCond('R',1);
   gsP.gaugeCond('L',1);
   //gsP.importMPS("/home/banulsm/SVN/peps/projects/mpsdyn/src/bin/mpsGS.mat");
   // I import the true GS!
   MPS gsPtrue(gsP);
#ifdef TESTMULTI
   strcpy(filename,LOCALDIR);
   gsPtrue.importMPS(strcat(filename,"mpsGS_N10.mat")); //   MPS gsM(gsP);
#else
   strcpy(filename,LOCALDIR);
   gsPtrue.importMPS(strcat(filename,"mpsGS.mat")); //   MPS gsM(gsP);
#endif
   
   // First: put H and find the GS
   int Dop=5; // bond dimension of Schwinger H
   double lambdaP=0.;
   //double lambdaM=0.;
   double offset=0.;
   SchwingerHamiltonian hSchP(L,mu,x,L0P,alpha,offset);
   //SchwingerHamiltonian hSchM(L,mu,x,L0M,alpha);
   const MPO& hamilP=hSchP.getHMPO();
   cout<<"Before starting the evolution, expectation value is="
       <<contractor.contract(gsP,hamilP,gsP)<<", and norm="
       <<contractor.contract(gsP,gsP)<<endl;
   cout<<"Overlap with true GS is="<<contractor.contract(gsP,gsPtrue)<<endl;
   // Perform the imaginary time evolution
   k=0;
   bool imag=true;
   //bool imag=false;
   while(k<M){
     double delta=deltas[k];
     MPS aux(gsP);
#ifdef TESTMULTI
     contractor.approximateMultiExponential(hamilP,aux,gsP,delta,BLOCK,D,imag);
     k=k+BLOCK-1;
#else
#ifdef RESCALE
     contractor.approximateExponential(hamilP,aux,gsP,delta,D,imag,-1175,2500);
#else
     contractor.approximateExponential(hamilP,aux,gsP,delta,D,imag,-1175);
#endif
#endif
     // Normalize the state
     gsP.gaugeCond('R',1);
     // Compute the energy
     complex_t lambdaP=contractor.contract(gsP,hamilP,gsP);
     cout<<"After step "<<k<<" of imaginary time evolution, "
	 <<"energy="<<lambdaP<<endl;
     complex_t overlap=contractor.contract(gsP,gsPtrue);
     cout<<"Overlap with true GS is="<<overlap<<endl;
     // How close am I from eigenstate? energy variance
     MPS aux2(gsP);
     contractor.optimize(hamilP,gsP,aux2);
     complex_t Evar=contractor.contract(aux2,aux2)-lambdaP*lambdaP;

     *out<<L<<"\t"<<x<<"\t"<<mg<<"\t"
	 <<D<<"\t" <<delta<<"\t" <<k+1<<"\t"<<lambdaP<<"\t"
	 <<(1./(2*L*x))*lambdaP<<"\t"<<abs(overlap)<<"\t"<<Evar<<endl;

     k++;
   }


   out->close();
}
