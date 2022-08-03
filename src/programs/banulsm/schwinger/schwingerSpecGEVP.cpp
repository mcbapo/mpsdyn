
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "SchwingerHamiltonian.h"
#include <cmath>
#include "mergesort.h"

using namespace std;
using namespace shrt;

#ifdef MACOSX
#define LOCALDIR "/Users/banulsm/Desktop/SVN/projects/mpsdyn/src/bin/"
#else
#define LOCALDIR "/home/banulsm/ALaPorra/bin/"
#endif
#define MAXLEN 120

void constructCorrelatorMatrix(int nrCorr,MPS** states,mwArray& Mmatr);
void constructReferenceMatrix(const mwArray& M0matr,mwArray& C0);

/** SchwingerSpecGEVP tries to find the scalar mass gap by combining
    DMRG and GEVP techniques. First we find the MPS approximation to
    the ground state and to some low lying excited states. The latter
    part is achieved by diagonalizing the effective Hamiltonian in the
    middle of the chain, a la DMRG. This provides an MPS ansatz for
    some low energy excitations. The DMRG idea is to use these as
    estimations of excited levels. Here we construct a correlation
    matrix from precisely these levels. We set t0=0, as we consider
    our initial state to be already close enough to the desired
    levels. They are then evolved in imaginary time, and the GEVP is
    solved at later times.

    In the output file we write the parameters used, the ground state
    found, and the lowest part of the spectrum of the effective
    Hamiltonian in the middle of the chain (with one site, for the
    moment), all as comments. Then the evolution is registered as for
    schwingerMS, i.e. every (uncommented) line of the file contains
    the step number, the time (k*delta) and the list of generalized
    eigenvalues.    

    \TODO Use some projector to try to ensure good quantum numbers?

    Receives arguments:
    \param <L> (int) length of the chain (if it is even, the site 
             to the right of the center will be taken for the effective H)
    \param <mg> (double) parameter \f$m/g\f$ of the Hamiltonian
    \param <alpha> (double) parameter \f$\alpha\f$ of the Hamiltonian
                            (actually, \f$\alpha+\ell_0\f$)
    \param <x> (double) parameter \f$x\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <nrCorr> (int) how many states to use for the correlation matrix
    \param <delta> (double) time step width
    \param <M> (int) maximum number of time steps to apply (i.e., evolve until t=M*delta)
    \param <outdirname> (char*) name of the output directory
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)
*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double mg=atof(argv[++cntr]);
  double alpha=atof(argv[++cntr]);
  double x=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  int nrCorr=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  const char* outdirname=argv[++cntr];
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  char filename[MAXLEN];
  sprintf(filename,"%s/%s",outdirname,outfname);

  double mu=2*mg*sqrt(x);
  double L0=0.; // l0 parameter, which is irrelevant, as only appears
		// as alpha+l0

  int midSites=1; // I just let one site change for the variational eigenstates

  //  double t0=.01; // todo: differently -> reference time

  cout<<"Initialized arguments: L="<<L
      <<", mu="<<mu
      <<", x="<<x
      <<", alpha="<<alpha
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;

  if(!app){
    out=new ofstream(filename);
  }
  else{
    out=new ofstream(filename,ios::app);
    *out<<"%%%%%%%%%%%%%%%%%%%"<<endl;
  }
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<filename<<
      " for output (append="<<app<<")"<<endl;
    exit(1);
  }

   Contractor& contractor=Contractor::theContractor();
   int d=2;
   MPS gs(L,D,d);
   gs.setRandomState(); // the intial state, random
   gs.gaugeCond('R',1);
   gs.gaugeCond('L',1);
   
   // First: put H and find the GS
   int Dop=5; // bond dimension of Schwinger H
   double lambda=0.;
   SchwingerHamiltonian hSch(L,mu,x,L0,alpha);
   const MPO& hamil=hSch.getHMPO();
   // To check!!
   //   hamil.exportForMatlab("schwHamil.m");
   contractor.findGroundState(hamil,D,&lambda,gs);
   cout<<"Ground state found with eigenvalue "<<lambda<<endl;
   *out<<"% L="<<L<<endl;
   *out<<"% mg="<<mg<<endl;
   *out<<"% x="<<x<<endl;
   *out<<"% alpha="<<alpha<<endl;
   *out<<"% D="<<D<<endl;
   *out<<setprecision(10);
   *out<<"% E0="<<lambda<<endl;

   
   // Place for the results of the diagonalization
   vector<complex_t> diagVals;
   mwArray U;

   int Lmid=L%2==0?L/2:(L-1)/2; 

   //   mwArray Heff=contractor.getEffectiveOperatorSingleSite(gs,hamil,Lmid);
   // cout<<"Effective H="<<endl;
   // cout<<Heff<<endl;
   // // Force Hermiticity:
   // Heff=.5*(Heff+Hconjugate(Heff));
   // wrapper::eigs(Heff,10,"SR",diagVals,U,false);


   int errEigs=0;
   if(midSites==1){
     // Now come the spectrum
     *out<<"%%% Spectrum of effective H on site "<<Lmid<<", dimension="<<d*D*D<<endl;
     TensorMultiplier Heff=contractor.getEffectiveOperatorMultiplierSingleSite(gs,hamil,Lmid);
     cout<<"Got the Heff, computing eigs!!"<<endl;
     errEigs=wrapper::eigs(Heff,nrCorr,"SR",diagVals,U,true);
   }
   else{ // not yet used, but maybe in the future?
     int Llow=midSites%2!=0?Lmid-(midSites-1)/2:Lmid-midSites/2;
     if(L%2!=0&&midSites%2==0) Llow+=1;
     *out<<"%%% Spectrum of effective H on "<<midSites<<" sites, starting on "<<Llow<<", dimension="<<pow(d,midSites)*D*D<<endl;
     TensorMultiplier Heff=contractor.getEffectiveOperatorMultiplierMultiSite(gs,hamil,midSites,Lmid);
     errEigs=wrapper::eigs(Heff,nrCorr,"SR",diagVals,U,true);
   }

   // Since from eigs I might get them in the reverse order, I have to sort them in ascending order
   cout<<"Computed values "<<diagVals<<endl;
   // cout<<"and matrix U="<<endl;  //cout<<U<<endl;
   cout<<setprecision(10);
   // cout<<"With error value "<<errEigs<<endl;
   vector<double> diagValsR(diagVals.size());
   for(int k=0;k<diagVals.size();k++){
     diagValsR[k]=real(diagVals[k]);
   }
   cout<<"REAL values "<<diagValsR<<endl;
   vector<int> indices(diagVals.size());
   for(int k=0;k<diagVals.size();k++) indices[k]=k;
   mergesort<double>(diagValsR,indices,false);
   cout<<"Sorted values "<<diagValsR<<" as "<<indices<<", i.e. "<<endl;
   for(int k=0;k<diagVals.size();k++)
     cout<<diagValsR[indices[k]]<<",";
   cout<<endl;

   int nrVals=diagVals.size();
   for(int k=0;k<nrVals;k++){
     *out<<"%"<<real(diagVals[k])<<"\t"<<imag(diagVals[k])<<endl;
   }

   // Construct the initial MPSs, as the GS and the excited states from the previous step
   MPS** states=new MPS*[nrCorr]; // place for pointers to the states
   states[0]=&gs; // the first one is the ground state, but I apply a certain (convenient) gauge
   gs.gaugeCond('R',true);
   for(int pos=L-1;pos>Lmid;pos--)
     gs.gaugeCond(pos,'L',true);
   for(int k=0;k<nrCorr;k++){ // prepare, one by one
     MPS* st=new MPS(gs); // first, initialize with gs
     // And then replace the middle site by the one computed in the diagonalization of the effective Hamiltonian,
     // but in the proper gauge, i.e., with the right tensors fulfilling gauge to the left!

     mwArray A=U.subArray(Indices(-1,indices[k]));A.reshape(Indices(-1,1));
     //cout<<"Extracted column "<<indices[k]<<" of U, "<<endl;
     //     cout<<A<<endl;
     //cout<<" with norm "<<Hconjugate(A)*A<<" and <Heff>="<<Hconjugate(A)*Heff*A<<endl;
     st->setA(Lmid,A);
     //st->gaugeCond('R',1);
     //cout<<"Computed state "<<k<<", with norm "<<contractor.contract(*st,*st)<<" and energy "<<contractor.contract(*st,hamil,*st)<<endl;
     states[k]=st;
  }

   cout<<"Computed now the states"<<endl;
   // At this moment, the correlation matrix, by construction, is the identity!! This is my perfect C0, because C0^(-1/2) is also the identity.

   // JUST TO CHECK!!!
   mwArray C0;
   constructCorrelatorMatrix(nrCorr,states,C0);
   cout<<"And the obtained correlation matrix is..."<<endl;
   cout<<C0<<endl;
   // JUST TO CHECK!!!

   // And now, evolve in imaginary time, and see

   double t1=0; // current time
   mwArray Mmatr(C0);
   
   // Perform the imaginary time evolution
   int k=0;
   complex_t deltaC=(complex_t){-delta,0.}; // delta  as complex
   complex_t deltaC_2=(complex_t){-delta*.5,0.}; // delta/2

   // Get the exponential MPOs
   MPO expHxe_2(L),expHxo(L),expHL(L),expHL_2(L);
   hSch.getExponentialMPOx(expHxe_2,deltaC_2,true);  
   hSch.getExponentialMPOx(expHxo,deltaC,false);
   hSch.getExponentialMPOL(expHL,deltaC);  
   hSch.getExponentialMPOL(expHL_2,deltaC_2);
   cout<<"Constructed all exponential MPOs"<<endl;

   bool imag=true;
   //bool imag=false;
   while(k<=M){
    mwArray U;vector<complex_t> eigVs;
    //      matr=.5*(matr+Hconjugate(matr));
    //wrapper::eig(Mmatr,eigVs,U,false); // not computing eigenvalues
    wrapper::eig(Mmatr,eigVs,U,true); // computing eigenvalues
    bool isNaN=false;
    *out<<k<<"\t"<<t1<<"\t";
    for(int rr=0;rr<nrCorr;rr++){
      if(std::isnan(real(eigVs[rr]))){
	isNaN=true;
	break; // do not write the NaN
      }
      else
	*out<<real(eigVs[rr])<<"\t";
    }
    *out<<endl;
    if(isNaN){
      cout<<"Not a Number in the diagonalization! Aborting (try smaller delta?)"<<endl;
      exit(1);
    }
     //    double delta_=deltas[k];
    double delta_=delta;
    if(delta_!=delta){ // It won't happen now
      cout<<"Recomputing exponential operators"<<endl;
      delta=delta_;
      deltaC=(complex_t){-delta,0.};
      deltaC_2=(complex_t){-delta*.5,0.};
      hSch.getExponentialMPOx(expHxe_2,deltaC_2,true);  
      hSch.getExponentialMPOx(expHxo,deltaC,false);
      //hSchP.getExponentialMPOL(expHL,deltaC);  
      hSch.getExponentialMPOL(expHL_2,deltaC_2);
    }     
    // Not very optimized: apply step by step
    for(int st=0;st<nrCorr;st++){
      // Evolve the desired states, one after the other
      MPS aux(*states[st]);MPS& st1=*states[st];
      contractor.optimize(expHL_2,aux,st1,D);
      aux=st1;
      contractor.optimize(expHxe_2,aux,st1,D);
      aux=st1;
      contractor.optimize(expHxo,aux,st1,D);
      aux=st1;
      contractor.optimize(expHxe_2,aux,st1,D);
      aux=st1;
      contractor.optimize(expHL_2,aux,st1,D);
      // Normalize the state ?
      // st1.gaugeCond('r',1);
      st1.setNormFact(st1.getNormFact()*exp(-abs(lambda)*delta)); // keep
								   // norms
								   // under control
    }
    k++;t1+=delta;
    constructCorrelatorMatrix(nrCorr,states,Mmatr);
   }


   out->close();
   delete out;

}

void constructCorrelatorMatrix(int nrCorr,MPS** states,mwArray& Mmatr){
  Contractor& contractor=Contractor::theContractor();
  Mmatr=mwArray(Indices(nrCorr,nrCorr));
  for(int k1=0;k1<nrCorr;k1++){
    complex_t diagK=contractor.contract(*states[k1],*states[k1]);
    Mmatr.setElement(.5*(diagK+conjugate(diagK)),Indices(k1,k1));
    for(int k2=k1+1;k2<nrCorr;k2++){
      complex_t K1K2=contractor.contract(*states[k2],*states[k1]);
      Mmatr.setElement(K1K2,Indices(k1,k2));
      Mmatr.setElement(conjugate(K1K2),Indices(k2,k1));
    }
  }
  //  cout<<"Constructed Mmatr="<<Mmatr<<endl;
  // should be hermitian, now
}

void constructReferenceMatrix(const mwArray& M0matr,mwArray& C0){
  mwArray U;vector<complex_t> eigVs;
  //cout<<"constructReferenceMatrix from M0="<<M0matr<<endl;

  wrapper::eig(M0matr,eigVs,U,true); // computing eigenvalues
  //cout<<"Reference matrix is "<<M0matr<<endl; 
  cout<<"It has eigenvalues "<<eigVs<<endl;

  //cout<<"Checking diagonalization of C0!"<<endl;
  //cout<<M0matr-U*diag(eigVs)*Hconjugate(U);
  // C0=U*1/sqrt(eigVs)*U+ 

  vector<complex_t> invSq;
  for(int k=0;k<eigVs.size();k++){
    invSq.push_back(ONE_c*1./sqrt(real(eigVs[k])));
    cout<<"In the diagonal of C0^(-1/2), "<<ONE_c/sqrt(eigVs[k])<<endl;
  }
  C0=U*diag(invSq)*Hconjugate(U);
}
