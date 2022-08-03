
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


/** Auxiliary function which constructs a blocked version of a MPO */
void blockMPO(const MPO& orig,MPO& final,int nrSites,int Llow);

/** Construct a block version given a list of Operator pointers. Returns a new Operator. */
Operator* blockOps(Operator** list,int nrsites);


/** SchwingerSpecGEVP2 tries to find the scalar mass gap by combining
    DMRG and GEVP techniques. First we find the MPS approximation to
    the ground state and to some low lying excited states. The latter
    part is achieved by diagonalizing the effective Hamiltonian in a
    few sites in the middle of the chain, a la DMRG. This provides an
    MPS ansatz for some low energy excitations. The DMRG idea is to
    use these as estimations of excited levels. Here we construct a
    correlation matrix from precisely these levels. But since the DMRG
    calculation includes all possible momenta, we select the levels to
    be used in the GEVP as the ones with the lowest value of
    \f$P^2\f$.  We set t0=0, as we consider our initial state to be
    already close enough to the desired levels. They are then evolved
    in imaginary time, and the GEVP is solved at later times.

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
    \param <nSt> (int) number of sites in the middle
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
  int nSt=atoi(argv[++cntr]);
  const char* outdirname=argv[++cntr];
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  char filename[MAXLEN];
  sprintf(filename,"%s/%s",outdirname,outfname);

  double mu=2*mg*sqrt(x);
  double L0=0.; // l0 parameter, which is irrelevant, as only appears
		// as alpha+l0

  int midSites=nSt; // I just let nSt sites change for the variational eigenstates

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
   int Llow=Lmid;
   if(midSites==1){
     // Now come the spectrum
     *out<<"%%% Spectrum of effective H on site "<<Lmid<<", dimension="<<d*D*D<<endl;
     TensorMultiplier Heff=contractor.getEffectiveOperatorMultiplierSingleSite(gs,hamil,Lmid);
     cout<<"Got the Heff, computing eigs!!"<<endl;
     errEigs=wrapper::eigs(Heff,nrCorr*10,"SR",diagVals,U,true);
   }
   else{ 
     Llow=midSites%2!=0?Lmid-(midSites-1)/2:Lmid-midSites/2;
     if(L%2!=0&&midSites%2==0) Llow+=1;
     *out<<"%%% Spectrum of effective H on "<<midSites<<" sites, starting on "<<Llow<<", dimension="<<pow(d,midSites)*D*D<<endl;
     TensorMultiplier Heff=contractor.getEffectiveOperatorMultiplierMultiSite(gs,hamil,midSites,Llow);
     errEigs=wrapper::eigs(Heff,nrCorr*10,"SR",diagVals,U,true);
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
   //   cout<<"REAL values "<<diagValsR<<endl;
   vector<int> indices(diagVals.size());
   for(int k=0;k<diagVals.size();k++) indices[k]=k;
   mergesort<double>(diagValsR,indices,false);
   // cout<<"Sorted values "<<diagValsR<<" as "<<indices<<", i.e. "<<endl;
   // for(int k=0;k<diagVals.size();k++)
   //   cout<<diagValsR[indices[k]]<<",";
   // cout<<endl;

   // The momentum operator
   MPO Pmpo(L),Pblocked(L-midSites+1);
   hSch.constructMomentumMPO(Pmpo);
   blockMPO(Pmpo,Pblocked,midSites,Llow);
   const MPO* oprsP[2]={&Pblocked,&Pblocked};
   MPO PblockedSq(L-midSites+1);
   MPO::join(2,oprsP,PblockedSq);


   // 2) Do the same with the ground state, but do not block sites
   // explicitely, as I will substitute them by the thing computed by
   // eigs
   gs.gaugeCond('R',true);
   for(int pos=L-1;pos>Llow+midSites-1;pos--)
     gs.gaugeCond(pos,'L',true);
   MPS auxSt(L-midSites+1,D,d);
   auxSt.gaugeCond('R',1);
   auxSt.gaugeCond('L',1);
   for(int k=0;k<L-midSites+1;k++){
     //cout<<"Copying old gs tensor to site "<<k<<" of the new one "<<endl;
     if(k<Llow)
       auxSt.setA(k,gs.getA(k).getA());
     else if(k>Llow)
       auxSt.setA(k,gs.getA(k+midSites-1).getA());
   }

   int nrVals=diagVals.size();
   vector<double> pvals(nrVals);
   complex_t PSqval0; //reference value(I want to choose the ones closest to the GS mom) 
   complex_t* vari=new complex_t[nrVals]; // place for the variance
   for(int k=0;k<nrVals;k++){
     *out<<"%"<<real(diagVals[k])<<"\t"<<imag(diagVals[k])<<"\t";
     mwArray A=U.subArray(Indices(-1,k));
     A.reshape(Indices(pow(d,midSites),auxSt.getA(Llow).getDl(),auxSt.getA(Llow).getDr()));
     auxSt.replaceSite(Llow,A);
     complex_t PSqval=contractor.contract(auxSt,PblockedSq,auxSt);
     if(k==0) PSqval0=PSqval;
     *out<<real(PSqval)<<"\t"<<imag(PSqval)<<endl;
     pvals[k]=abs(real(PSqval-PSqval0));
   }

   // Now select the ones with the lowest momentum
   // First, order them
   vector<int> indicesP(nrVals);
   for(int k=0;k<nrVals;k++) indicesP[k]=k;
   mergesort<double>(pvals,indicesP,false);

  // Construct the initial MPSs, as the GS and the excited states from the previous step
   MPS** states=new MPS*[nrCorr]; // place for pointers to the states
   states[0]=&gs; // the first one is the ground state, oalready in the proper gauge
   for(int k=0;k<nrCorr;k++){ // prepare, one by one
     MPS* st=new MPS(gs); // first, initialize with gs
     // And then replace the middle sites by the one computed in the diagonalization of the effective Hamiltonian,
     // but in the proper gauge, i.e., with the right tensors fulfilling gauge to the left!
     mwArray A=U.subArray(Indices(-1,indicesP[k]));
     A.reshape(Indices(-1,1));
     //cout<<"Looking at column "<<indices[k]<<" of the diagonalization, size "<<A.getDimensions()<<endl;
     if(nSt==1)
       st->setA(Lmid,A);
     else{
       // Divide the tensor (SVD? how to enlarge D?)
       mwArray W,S,Vdagger; // A=W*(S*Vdagger)
       int Dl=gs.getA(Llow).getDl(); // leftmost bond
       int Dr=gs.getA(Llow+nSt-1).getDr(); // rightmost bond
       int restdim=pow(d,nSt-1); //composite physical dimension
       A.reshape(Indices(d*restdim,Dl,Dr));
       //cout<<"The vector found is now "<<A.getDimensions()<<endl;
       A.permute(Indices(2,1,3)); // Dl x dphys x Dr
       for(int il=0;il<nSt-1;il++){
	 A.reshape(Indices(Dl*d,restdim*Dr));
	 wrapper::svd(A,W,S,Vdagger);
	 // cout<<"Result of svd of A="<<A<<endl;
	 // cout<<"W="<<W<<endl;
	 // cout<<"S="<<S<<endl;
	 // cout<<"Vdagger="<<Vdagger<<endl;
	 W.multiplyRight(S);
	 int nr=S.getDimension(1);
	 W.reshape(Indices(Dl,d,nr));
	 st->setRotatedA(Llow+il,W,Indices(2,1,3),false,false); // do not check dimensions
	 if(il==nSt-2){ // last one, so I set the last tensor to the rest
	   Vdagger.reshape(Indices(nr,d,Dr));
	   st->setRotatedA(Llow+il+1,Vdagger,Indices(2,1,3),false,false);
	 }
	 else{
	   A=Vdagger;
	   restdim=restdim/d;
	   Dl=nr;
	 }
       }
     }
     //st->gaugeCond('R',1);
     //cout<<"MPS st"<<*st<<endl;
     complex_t normK=contractor.contract(*st,*st);
     complex_t expecH=contractor.contract(*st,hamil,*st);
     complex_t expecH2=contractor.contract2(hamil,*st);
     vari[k]=expecH2/normK-expecH*expecH/(normK*normK);
     cout<<"Computed state "<<k<<", with norm "<<normK<<", energy "<<expecH
	 <<" and <H^2>="<<expecH2<<" => Delta H2="<<vari[k]<<endl;
     states[k]=st;
  }

   cout<<"Computed now the states"<<endl;
   // At this moment, the correlation matrix, by construction, is the identity!! This is my perfect C0, because C0^(-1/2) is also the identity.

   // JUST TO CHECK!!!
   mwArray C0;
   constructCorrelatorMatrix(nrCorr,states,C0);
   // cout<<"And the obtained correlation matrix is..."<<endl;
   //cout<<C0<<endl;
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
   delete [] vari;
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


void blockMPO(const MPO& orig,MPO& final,int midSites,int Llow){
  int L=orig.getLength();
  final.initLength(L-midSites+1);
   for(int k=0;k<L-midSites+1;k++){
     if(k<Llow){
       //final.setOp(k,&orig.getOp(k),false);
       final.setOp(k,new Operator(orig.getOp(k).getFullData()),true);
     }
     else if(k>Llow){
       // final.setOp(k,&orig.getOp(k+midSites-1),false);
       final.setOp(k,new Operator(orig.getOp(k+midSites-1).getFullData()),true);
     }
     else{ // the middle site!
       mwArray effSite=orig.getOpData(Llow); // du dl dd dr
       int du=effSite.getDimension(0);
       int dd=effSite.getDimension(2);
       int Dl=effSite.getDimension(1);int Dr=effSite.getDimension(3);
       for(int p=1;p<midSites;p++){
	 effSite.reshape(Indices(du*Dl*dd,Dr));
	 mwArray auxOp=orig.getOpData(Llow+p); // du Dr dd Drp
	 Indices newdims=auxOp.getDimensions();
	 // cout<<"Adding site "<<p+Llow<<" to effective operator, dims now are"
	 //     <<effSite.getDimensions()<<" and the new piece has "<<newdims<<endl;
	 auxOp.permute(Indices(2,1,3,4));
	 auxOp.reshape(Indices(newdims[1],newdims[0]*newdims[2]*newdims[3]));
	 effSite.multiplyRight(auxOp);
	 effSite.reshape(Indices(du,Dl,dd,newdims[0],newdims[2],newdims[3]));
	 effSite.permute(Indices(1,4,2,3,5,6));
	 du=du*newdims[0];dd=dd*newdims[2];
	 effSite.reshape(Indices(du,Dl,dd,Dr));
       }
       final.setOp(k,new Operator(effSite),true);
     }
   }
}


 Operator* blockOps(Operator** list,int nrsites){
   mwArray aux=list[0]->getFullData();
   int dimU=aux.getDimension(0);
   int Dl=aux.getDimension(1);
   int dimD=aux.getDimension(2);
   int Dr=aux.getDimension(3);
   for(int s=1;s<nrsites;s++){
     aux.reshape(Indices(dimU*Dl*dimD,Dr));
     mwArray aux2=list[s]->getFullData();
     Indices dim2=aux2.getDimensions();
     aux2.permute(Indices(2,1,3,4));
     aux2.reshape(Indices(dim2[1],dim2[0]*dim2[2]*dim2[3]));
     aux.multiplyRight(aux2);
     aux.reshape(Indices(dimU,Dl,dimD,dim2[0],dim2[2],dim2[3]));
     aux.permute(Indices(1,4,2,3,5,6));
     dimU*=dim2[0];
     dimD*=dim2[2];
     Dr=dim2[3];
     aux.reshape(Indices(dimU,Dl,dimD,Dr));
   }
   return new Operator(aux);
}
