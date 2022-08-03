
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "SchwingerHamiltonian.h"

#ifdef TESTMULTI
#define BLOCK 4
#endif

/** SchwingerIT2 tries to find the ground state of the Schwinger
    Hamiltonian \ref <SchwingerHamiltonian> using imaginary time
    evolution. In this case, instead of the Cayley like
    approximations, a Trotter decomposition is used, such that the x
    term is approximated by a product of even/odd exponentials, while
    the long-range term still uses the Cayley approach (maybe with
    offset required).

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

  double mu=2*mg*sqrt(x);

  double alpha=.5*0.;
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
    *out<<"% N\t x\t mg\t D\t delta\t M\t Energy\t w0/(2Nx)"
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
  MPS aux(gsP); // I will be testing my gauge here
  //gsP.setProductState(p_xplus);
  gsP.gaugeCond('R',1);
  gsP.gaugeCond('L',1);
  //gsP.importMPS("/home/banulsm/SVN/peps/projects/folding/src/bin/mpsGS.mat");
  // I import the true GS!
  MPS gsPtrue(gsP);
  //#ifdef TESTMULTI
  //gsPtrue.importMPS("/home/banulsm/SVN/peps/projects/mpsdyn/src/bin/mpsGS_N10.mat");   //   MPS gsM(gsP);
  //#else
  //gsPtrue.importMPS("/home/banulsm/SVN/peps/projects/folding/src/bin/mpsGS.mat");   //   MPS gsM(gsP);
  //#endif
  // First: put H and find the GS
  int Dop=5; // bond dimension of Schwinger H
  double lambdaP=0.;
  //double lambdaM=0.;
  double offset=0.;
  SchwingerHamiltonian hSchP(L,mu,x,L0P,alpha,offset);
  // A trick to construct H_L alone
  //SchwingerHamiltonian hSchL(L,mu,0.,L0P,alpha,offset);
  //SchwingerHamiltonian hSchM(L,mu,x,L0M,alpha);
  const MPO& hamilP=hSchP.getHMPO();

  // aux.exportMPS("aRandomMPS.mat");
  // clock_t start=clock();
  // aux.gaugeCond('R',0);
  // clock_t timeQR=clock()-start;
  // MPS aux2(aux);
  // start=clock();
  // aux2.gaugeCond('R',0,D);
  // clock_t timeSVD=clock()-start;
  // cout<<"With QR gauge (time "<<timeQR<<"), aux H aux="<<contractor.contract(aux,hamilP,aux)
  //     <<" and with SVD gauge ("<<timeSVD<<"), "<<contractor.contract(aux2,hamilP,aux2)
  //     <<", overlap "<<contractor.contract(aux,aux2)<<endl;

  // aux.importMPS("aRandomMPS.mat");
  // aux2.importMPS("aRandomMPS.mat");
  // start=clock();
  // aux.gaugeCond('L',0);timeQR=clock()-start;
  // start=clock();
  // aux2.gaugeCond('L',0,D);timeSVD=clock()-start;
  // cout<<"With LQ gauge (time "<<timeQR<<"), aux H aux="<<contractor.contract(aux,hamilP,aux)
  //     <<" and with SVD gauge (L) (time "<<timeSVD<<") "<<contractor.contract(aux2,hamilP,aux2)
  //     <<", overlap "<<contractor.contract(aux,aux2)<<endl;
  //  exit(1);


  //const MPO& hamilL=hSchL.getHMPO();
  double E0=0;
  contractor.findGroundState(hamilP,D,&E0,gsPtrue);
   cout<<"True GS (for bond dimensions D="<<D<<") has energy E0="<<E0<<endl;
   
   cout<<"Before starting the evolution, expectation value is="
       <<contractor.contract(gsP,hamilP,gsP)<<", and norm="
       <<contractor.contract(gsP,gsP)<<endl;
   cout<<"Overlap with true GS is="<<abs(contractor.contract(gsP,gsPtrue))<<endl;
   
   double delta=deltas[0];
   complex_t deltaC=(complex_t){-delta,0.}; // delta  as complex
   complex_t deltaC_2=(complex_t){-delta*.5,0.}; // delta/2
   
   // Get the exponential MPOs
   MPO expHxe_2(L),expHxo(L),expHL(L),expHL_2(L);
   hSchP.getExponentialMPOx(expHxe_2,deltaC_2,true);  
   hSchP.getExponentialMPOx(expHxo,deltaC,false);
   hSchP.getExponentialMPOL(expHL,deltaC);  
   hSchP.getExponentialMPOL(expHL_2,deltaC_2);
   cout<<"Constructed all exponential MPOs"<<endl;

   //cout<<"Some checking now?"<<endl;
   //gsP.exportMPStext("initMPS.txt");

   //expHL_2.exportMPOtext("expHL2.txt");
   //expHL.exportMPOtext("expHL.txt");
   //hamilL.exportMPOtext("mpoHL.txt");
   *out<<setprecision(10);

   // Perform the imaginary time evolution
   k=0;double time=0;
   bool imag=true;
   //bool imag=false;
   while(k<M){
     double delta_=deltas[k];
     if(delta_!=delta){
       cout<<"Recomputing exponential operators"<<endl;
       delta=delta_;
       deltaC=(complex_t){-delta,0.};
       deltaC_2=(complex_t){-delta*.5,0.};
       hSchP.getExponentialMPOx(expHxe_2,deltaC_2,true);  
       hSchP.getExponentialMPOx(expHxo,deltaC,false);
       //hSchP.getExponentialMPOL(expHL,deltaC);  
       hSchP.getExponentialMPOL(expHL_2,deltaC_2);
     }
     MPS aux(gsP);
     if(1){
       // Not very optimized: apply step by step
       contractor.optimize(expHL_2,aux,gsP,D);
       aux=gsP;
       contractor.optimize(expHxe_2,aux,gsP,D);
       aux=gsP;
       contractor.optimize(expHxo,aux,gsP,D);
       aux=gsP;
       contractor.optimize(expHxe_2,aux,gsP,D);
       aux=gsP;
       contractor.optimize(expHL_2,aux,gsP,D);
     }
     else{
       contractor.approximateMultiExponential(hamilP,aux,gsP,delta,2,D,true);
       //      contractor.approximateExponential(hamilP,aux,gsP,delta,D,true);
       // Cutre-option for higher order Taylor
       //       contractor.approximateTaylorExponential(hamilP,aux,gsP,delta,1,D,true);
       //contractor.approximateTaylorExponential(hamilP,aux,gsP,delta,1,D,true);

     }
     // Normalize the state
     gsP.gaugeCond('R',1);
     // Compute the energy
     complex_t lambdaP=contractor.contract(gsP,hamilP,gsP);
     cout<<"After step "<<k<<" of imaginary time evolution, "
	 <<"energy="<<lambdaP<<" true GS "<<E0<<endl;
     complex_t overlap=contractor.contract(gsP,gsPtrue);
     cout<<"Overlap with true GS is="<<abs(overlap)<<endl;
     // How close am I from eigenstate? energy variance
     MPS aux2(gsP);
     //     contractor.optimize(hamilP,gsP,aux2);
     //     complex_t Evar=contractor.contract(aux2,aux2)-lambdaP*lambdaP;

     time+=delta;
     *out<<L<<"\t"<<x<<"\t"<<mg<<"\t"
	 <<D<<"\t" <<delta<<"\t" <<time<<"\t"<<real(lambdaP)<<"\t"
	 <<(1./(2*L*x))*real(lambdaP)<<"\t"<<abs(overlap)
	 <<"\t"<<abs((real(lambdaP)-E0)/E0)
	 <<endl; //<<"\t"<<Evar<<endl;

     k++;
   }


   out->close();
}
