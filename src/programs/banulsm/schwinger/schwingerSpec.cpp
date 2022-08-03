
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
//#include "HermitianTensorMultiplier.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "SchwingerHamiltonian.h"
#include "SpinMPO.h"
#include <cmath>

using namespace std;
using namespace shrt;

#ifdef MACOSX
#define LOCALDIR "/Users/banulsm/Desktop/SVN/projects/mpsdyn/src/bin/"
#else
#define LOCALDIR "/home/banulsm/ALaPorra/bin/"
#endif
#define MAXLEN 120


/** Auxiliary function which constructs a blocked version of a MPO */
void blockMPO(const MPO& orig,MPO& final,int nrSites,int Llow);

/** Construct a block version given a list of Operator pointers. Returns a new Operator. */
Operator* blockOps(Operator** list,int nrsites);

/** SchwingerSpec tries to find the scalar mass gap by first finding
    the MPS approximation to the ground state, then diagonalizing the
    effective Hamiltonian in the middle of the chain, a la DMRG.  In
    the output file I will write the parameters used, the ground state
    found, and the full spectrum of the effective Hamiltonian in the
    middle of the chain.

    Receives arguments:
    \param <L> (int) length of the chain (if it is even, the site 
             to the right of the center will be taken for the effective H)
    \param <mg> (double) parameter \f$m/g\f$ of the Hamiltonian
    \param <alpha> (double) parameter \f$\alpha\f$ of the Hamiltonian
                            (actually, \f$\alpha+\ell_0\f$)
    \param <x> (double) parameter \f$x\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <midSites> (int) how many sites to use for the effective Hamiltonian (one, two) 
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
  int midSites=atoi(argv[++cntr]);
  const char* outdirname=argv[++cntr];
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  int nrCorr=50;

  char filename[MAXLEN];
  sprintf(filename,"%s/%s",outdirname,outfname);

  double mu=2*mg*sqrt(x);
  double L0=0.; // l0 parameter, which is irrelevant, as only appears
		// as alpha+l0

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

   { // try to compute the expectation value of the condensate, too
     MPO cond(L);
     hSch.constructCondensateMPO(cond);
     complex_t fC=contractor.contract(gs,cond,gs);
     *out<<"% <condensate>="<<real(fC)<<endl;
   }
   //mwArray Heff=contractor.getEffectiveOperatorSingleSite(gs,hamil,Lmid);
   // // Force Hermiticity:
   // Heff=.5*(Heff+Hconjugate(Heff));
   // wrapper::eigs(Heff,10,"SR",diagVals,U,false);

   // // I could also get eigenvectors, and compute the <H^2>, to check how close to an eigenstate I am!
   
   // Place for the results of the diagonalization
   vector<complex_t> diagVals;
   mwArray U;

   int Lmid=L%2==0?L/2:(L-1)/2; 
   int Llow=Lmid;
   if(midSites==1){
     // Now come the spectrum
     *out<<"%%% Spectrum of effective H on site "<<Lmid<<", dimension="<<d*D*D<<endl;
     TensorMultiplier Heff=contractor.getEffectiveOperatorMultiplierSingleSite(gs,hamil,Lmid);
     cout<<"Got the Heff, computing eigs!!"<<endl;
     //wrapper::eigs(Heff,10,"SR",diagVals,U,false);
     wrapper::eigs(Heff,nrCorr,"SR",diagVals,U,true);
   }
   else{
     Llow=midSites%2!=0?Lmid-(midSites-1)/2:Lmid-midSites/2;
     if(L%2!=0&&midSites%2==0) Llow+=1;
     *out<<"%%% Spectrum of effective H on "<<midSites<<" sites, starting on "<<Llow<<", dimension="<<pow(d,midSites)*D*D<<endl;
     TensorMultiplier Heff=contractor.getEffectiveOperatorMultiplierMultiSite(gs,hamil,midSites,Llow);
     //HermitianTensorMultiplier Heff_(Heff);
     //wrapper::eigs(Heff,10,"SR",diagVals,U,false);
     wrapper::eigs(Heff,nrCorr,"SR",diagVals,U,true);
   }

   // I want to check how close to eigenvalues these are
   // 1) construct a "fake" MPO with the middle sites grouped
   MPO Hblocked(L-midSites+1);
   blockMPO(hamil,Hblocked,midSites,Llow);
   // for(int k=0;k<L-midSites+1;k++){
   //   if(k<Llow)
   //     Hblocked.setOp(k,&hamil.getOp(k),false);
   //   else if(k>Llow)
   //     Hblocked.setOp(k,&hamil.getOp(k+midSites-1),false);
   //   else{ // the middle site!
   //     mwArray effSite=hamil.getOpData(Llow); // du dl dd dr
   //     int du=d;int dd=d;int Dl=effSite.getDimension(1);int Dr=effSite.getDimension(3);
   //     for(int p=1;p<midSites;p++){
   // 	 effSite.reshape(Indices(du*Dl*dd,Dr));
   // 	 mwArray auxOp=hamil.getOpData(Llow+p); // du Dr dd Drp
   // 	 Indices newdims=auxOp.getDimensions();
   // 	 cout<<"Adding site "<<p+Llow<<" to effective H, dims now are "<<effSite.getDimensions()<<" and the new piece has "<<newdims<<endl;
   // 	 auxOp.permute(Indices(2,1,3,4));
   // 	 auxOp.reshape(Indices(newdims[1],newdims[0]*newdims[2]*newdims[3]));
   // 	 effSite.multiplyRight(auxOp);
   // 	 effSite.reshape(Indices(du,Dl,dd,newdims[0],newdims[2],newdims[3]));
   // 	 effSite.permute(Indices(1,4,2,3,5,6));
   // 	 du=du*newdims[0];dd=dd*newdims[2];
   // 	 effSite.reshape(Indices(du,Dl,dd,Dr));
   //     }
   //     Hblocked.setOp(k,new Operator(effSite),true);
   //   }
   // }

   // Compute the expectation value of Sz, as it commutes with H
   MPO Szmpo(L),Szblocked(L-midSites+1);
   SpinMPO::getSzMPO(L,d,Szmpo);
   blockMPO(Szmpo,Szblocked,midSites,Llow);
   Szmpo.exportForMatlab("Szmpo.m");

   // The momentum operator
   MPO Pmpo(L),Pblocked(L-midSites+1);
   hSch.constructMomentumMPO(Pmpo);
   blockMPO(Pmpo,Pblocked,midSites,Llow);
   Pmpo.exportForMatlab("Pmpo.m");
   const MPO* oprsP[2]={&Pblocked,&Pblocked};
   MPO PblockedSq(L-midSites+1);
   MPO::join(2,oprsP,PblockedSq);

   // The momentum squared operator
   MPO P2mpo(L),P2blocked(L-midSites+1);
   hSch.constructMomentumSquaredMPO(P2mpo);
   blockMPO(P2mpo,P2blocked,midSites,Llow);
   P2mpo.exportForMatlab("P2mpo.m");

   // And the momentum distribution
   //MPO Nkmpo(L),Nkblocked(L-midSites+1);
   // cout<<"Initializing the relevant momentum occupation operators"<<endl; 
   // MPO** Skmpos=new MPO*[L-1]; // the blocked ones!
   // MPO Skmpo(L); // the unblocked (reuse it)
   // for(int z=0;z<L;z++){
   //   SpinMPO::getOBCSkMPO(L,d,z+1,Skmpo,false);
   //   Skmpos[z]=new MPO(L-midSites+1);
   //   blockMPO(Skmpo,*Skmpos[z],midSites,Llow);       
   // }

   cout<<"Now computing properties of the eigenvectors"<<endl;

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

   // Do that also for the translation-by-two operator
   MPO Tmpo(L),TLmpo(L), Tblocked(L-midSites+1),TLblocked(L-midSites+1);
   hSch.constructCyclicTranslationMPO(Tmpo); // T to right
   hSch.constructCyclicTranslationMPO(TLmpo,false); // T' (to left)
   blockMPO(Tmpo,Tblocked,midSites,Llow);   
   blockMPO(TLmpo,TLblocked,midSites,Llow);   
   MPO T2blocked(L-midSites+1);
   MPO T4blocked(L-midSites+1);
   MPO T2LRblocked(L-midSites+1);
   MPO T2RLblocked(L-midSites+1);
   {   
     const MPO* oprs[2]={&Tblocked,&Tblocked};
     const MPO* oprs4[4]={&Tblocked,&Tblocked,&Tblocked,&Tblocked};
     const MPO* oprsLR[4]={&TLblocked,&TLblocked,&Tblocked,&TLblocked};
     const MPO* oprsRL[4]={&Tblocked,&TLblocked,&TLblocked,&TLblocked};
     MPO::join(2,oprs,T2blocked);
     MPO::join(4,oprs4,T4blocked);
     MPO::join(4,oprsLR,T2LRblocked);
     MPO::join(4,oprsRL,T2RLblocked);
   }
   // I also want the probability distribution of fermions, which should be something like
   // (1+sigma_z)_e/2+(1+sigma_z)_o/2
   // For this, I need a local operator
   complex_t locO[]={ONE_c,ZERO_c,ZERO_c,ZERO_c};
   mwArray locOp(Indices(d,1,d,1),locO);
   mwArray locI=identityMatrix(d);locI.reshape(Indices(d,1,d,1));
   Operator* Nop=new Operator(locOp);
   Operator* Idop=new Operator(locI);
   // the identity for the blocked sites
   mwArray blockI=identityMatrix(pow(d,midSites));blockI.reshape(Indices(pow(d,midSites),1,pow(d,midSites),1));
   Operator* Idblock=new Operator(blockI);;
   // Construct by hand the blocked operators
   Operator** Nblocks=new Operator*[midSites];
   for(int s=0;s<midSites;s++){
     Operator** listOps=new Operator*[midSites];
     for(int ls=0;ls<midSites;ls++){
       if(ls==s)
	 listOps[ls]=Nop;
       else
	 listOps[ls]=Idop;
     }
     Nblocks[s]=blockOps(listOps,midSites);
     delete [] listOps;
   }
   // Prepare the MPO with identities
   MPO probMPO(L-midSites+1);
   for(int l=0;l<L-midSites+1;l++)
     if(l!=Llow)
       probMPO.setOp(l,Idop,false);
     else
       probMPO.setOp(l,Idblock,false);


   int nrVals=diagVals.size();
   for(int k=0;k<nrVals;k++){
     cout<<"Nr "<<k<<" (E="<<diagVals[k]<<")"<<endl;
     *out<<real(diagVals[k])<<"\t"<<imag(diagVals[k])<<"\t";
     // compute the expectation value of H and H^2
     mwArray A=U.subArray(Indices(-1,k));
     A.reshape(Indices(pow(d,midSites),auxSt.getA(Llow).getDl(),auxSt.getA(Llow).getDr()));
     auxSt.replaceSite(Llow,A);
     complex_t H2=contractor.contract2(Hblocked,auxSt);
     complex_t normN=contractor.contract(auxSt,auxSt);
     complex_t Hval=contractor.contract(auxSt,Hblocked,auxSt);
     complex_t Szval=contractor.contract(auxSt,Szblocked,auxSt);
     complex_t Pval=contractor.contract(auxSt,Pblocked,auxSt);
     complex_t P2val=contractor.contract(auxSt,P2blocked,auxSt);
     complex_t PSqval=contractor.contract(auxSt,PblockedSq,auxSt);
     complex_t T2val=contractor.contract(auxSt,T2blocked,auxSt);
     // I'd also like to compute [i*(T2-T2')]^2, which should be sth like the (2nd) derivative
     // complex_t T4LRval=contractor.contract(auxSt,T2LRblocked,auxSt);
     // complex_t T4RLval=contractor.contract(auxSt,T2RLblocked,auxSt);
     complex_t T4LRval=contractor.contract2(T2blocked,auxSt);
     complex_t T4RLval=contractor.contract2(auxSt,T2blocked);
     complex_t T4val=contractor.contract(auxSt,T4blocked,auxSt);
     double secDer=2*real(T4val)-real(T4LRval)-real(T4RLval);

     *out<<real(normN)<<"\t";
     *out<<real(H2)<<"\t"; //<<imag(H2)<<"\t";
     *out<<real(Hval)<<"\t"; //<<imag(Hval)<<"\t";
     //*out<<real(H2/normN-Hval*Hval/(normN*normN))<<"\t";
     *out<<real(Szval)<<"\t";
     *out<<real(Pval)<<"\t";
     *out<<real(P2val)<<"\t";
     *out<<real(PSqval)<<"\t";
     *out<<real(T2val)<<"\t"<<imag(T2val)<<"\t";
     *out<<secDer<<"\t";
     // compute the probability distribution
     if(1){
       for(int s=0;s<L;s++){
	 int ptr=s;
	 if(s<Llow)
	   probMPO.setOp(s,Nop,false);
	 else if(s>=Llow+midSites){
	   ptr=s-midSites+1;
	   probMPO.setOp(ptr,Nop,false);
	 }
	 else{
	   ptr=Llow;
	   probMPO.setOp(Llow,Nblocks[s-Llow],false);
	 }
	 complex_t probVal=contractor.contract(auxSt,probMPO,auxSt);
	 *out<<real(probVal)<<"\t";
	 if(ptr==Llow)
	   probMPO.setOp(ptr,Idblock,false);
	 else
	   probMPO.setOp(ptr,Idop,false);
       }
     }
     //*out<<real(Tval)<<"\t"<<imag(Tval)<<"\t";
     // // And the momentum distribution
     // for(int kval=0;kval<L;kval++){
     //   //cout<<"About to retrieve Nk mpo for k="<<kval<<endl;
     //   //SpinMPO::getMomentumDistributionMPO(L,d,kval,Nkmpo,false);
     //   //blockMPO(Nkmpo,Nkblocked,midSites,Llow);       
     //   //Skmpo.exportForMatlab("Skmpo.m");exit(1);
     //   complex_t nkval=contractor.contract2(*Skmpos[kval],auxSt);
     //   *out<<real(nkval)<<"\t";
     //   if(abs(imag(nkval))>1E-10){
     // 	 cout<<"ERROR: Non-real expectation value for n_k with k="<<kval
     // 	     <<" ("<<nkval<<")"<<endl;
     // 	 exit(1);
     //   }
     // }
     *out<<endl;
   }
   out->close();
   delete out;

   // for(int z=1;z<L;z++){
   //   delete Skmpos[z];
   // }
   // delete []Skmpos;
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
