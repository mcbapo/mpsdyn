
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "SchwingerHamiltonian.h"
#include "SpinMPO.h"

using namespace std;
using namespace shrt;

#ifdef MACOSX
#define LOCALDIR "/Users/banulsm/Desktop/SVN/projects/mpsdyn/src/bin/"
#else
#define LOCALDIR "/home/banulsm/ALaPorra/bin/"
#endif
#define MAXLEN 120

/** Initialize a set of states by  */
MPS** initializeVStates(int nr,const MPO& Vmpo,const MPS& gs,int D);
/** Initialize the vectors by applying just two-body terms  */
MPS** initializeVStatesAlt(int nr,const MPO& Vmpo,const MPS& gs,int D,const MPO& theOp);
/** Trampa para inicializar los estados, unos con Vec y otros con Scal */
MPS** initializeVStatesTR(int nr,const MPO& Vmpo,const MPS& gs,int D,SchwingerHamiltonian& hSch);
MPS** initializeRandomStates(int nr,const MPS& gs,int D);
void deleteStates(MPS** array,int nr);
void constructCorrelatorMatrix(int nrCorr,MPS** states,mwArray& Mmatr);
void printMatrixElements(ofstream& out,const mwArray& Mmatr);
void constructReferenceMatrix(const mwArray& M0Matr,mwArray& C0);
void constructOverlapMatrices(mwArray& expecV,mwArray& overl,const MPO& vecMPO,MPS** states,int nrCorr);
void constructSpinMatrices(mwArray& expecSp,mwArray& overl,const MPO& sxMPO,
			   const MPO& syMPO,const MPO& szMPO,MPS** states,int nrCorr);


/** SchwingerMS tries to find the scalar mass gap by first finding the
    MPS approximation to the ground state, then applying the operator
    V one or two times and evolving in imaginary time the resulting
    states. 

    Receives arguments:
    \param <L> (int) length of the chain
    \param <mg> (double) parameter \f$m/g\f$ of the Hamiltonian
    \param <alpha> (double) parameter \f$\alpha\f$ of the Hamiltonian
                            (actually, \f$\alpha+\ell_0\f$)
    \param <x> (double) parameter \f$x\f$ of the Hamiltonian
    \param <t0> (double) reference time
    \param <delta> (double) time step (it will run until t=2*t0)
    \param <D> (int) maximum bond dimension
    \param <nrCorr> (int) how many "states" are to be kept (dim of corr. 
                             matrix, thus cost increases as nCorr)
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
  double t0=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  int nrCorr=atoi(argv[++cntr]);
  const char* outdirname=argv[++cntr];
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  char filename[MAXLEN];
  char filenameE[MAXLEN];
  sprintf(filename,"%s/%s",outdirname,outfname);
  sprintf(filenameE,"%s/Eigs_%s",outdirname,outfname);

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
  ofstream *out2;

  if(!app){
    //    out=new ofstream(filename);
    //*out<<"% step\t t\t Upper diagonal of the Corr Matr"<<endl;
    char prefix[MAXLEN];
    out2=new ofstream(filenameE);
    *out2<<"% step\t t\t eigenvalues"<<endl;
  }
  else{
    //out=new ofstream(filename,ios::app);
    out2=new ofstream(filenameE,ios::app);
  }
  if(!out2->is_open()){
    cout<<"Error: impossible to open file "<<outdirname<<"/Eigs_"<<outfname<<
      " for output"<<endl;
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
   if(!app)
     *out2<<"% E0="<<lambda<<endl;
   // I want to know the spin(z) of the ground state
   MPO sx(L),sy(L),sz(L);
   //MPO sz(L);
   SpinMPO::getSxMPO(L,d,sx);
   SpinMPO::getSyMPO(L,d,sy);
   SpinMPO::getSzMPO(L,d,sz);
   complex_t sz_=contractor.contract(gs,sz,gs);
   //   complex_t s2=contractor.contract2(sx,gs)+contractor.contract2(sy,gs)
   //+contractor.contract2(sz,gs);
   cout<<"Ground state has spin Z "<<sz_/contractor.contract(gs,gs)<<endl;
   //sx.exportForMatlab("sx.m");
   //sy.exportForMatlab("sy.m");
   //sz.exportForMatlab("sz.m");
   //gs.exportForMatlab("schwGS.m");
   //   MPO par(L);
   // hSch.constructParityMPO(par);
   // cout<<"Ground state has parity "<<contractor.contract(gs,par,gs)/contractor.contract(gs,gs)<<endl;
   //MPS aux;
   //complex_t parGS=hSch.computeParity(gs);
   //cout<<"Ground state has parity "<<parGS/contractor.contract(gs,gs)<<endl;

   // If applied twice, should recover the original!
   // MPS auxP;
   // hSch.applyParity(aux,auxP);
   // cout<<"After applying it twice, the norm is "<<contractor.contract(auxP,auxP)<<", and the overlap with the original, which should be one, is "<<contractor.contract(gs,auxP)<<endl;

   // //  But what is the parity symm. part? (1+P)|Psi>
   // MPS aux2(aux);
   // vector<const MPS*> terms;terms.push_back(&gs);terms.push_back(&aux);
   // vector<complex_t> beta(2,.5*ONE_c);
   // contractor.optimizeSum(terms,beta,aux2,4*D);
   // MPS aux2a(aux);beta[1]=-1.*beta[1]; // antisymmetric part
   // contractor.optimizeSum(terms,beta,aux2a,4*D);
   // cout<<"The norm of the symmetric part is "<<contractor.contract(aux2,aux2)<<endl;
   // cout<<"The norm of the antisymmetric part is "<<contractor.contract(aux2a,aux2a)<<endl;
   // cout<<"And their overlap (should be zero) "<<contractor.contract(aux2,aux2a)<<endl;
   // MPS aux3;   hSch.applyParity(aux2,aux3);
   // cout<<"And its parity (sym) is "<<contractor.contract(aux2,aux3)/contractor.contract(aux2,aux2)<<endl;
   // hSch.applyParity(aux2a,aux3);
   // cout<<"And its parity (asym) is "<<contractor.contract(aux2a,aux3)/contractor.contract(aux2a,aux2a)<<endl;

   // exit(1);

   //   strcpy(filename,LOCALDIR);
   //   gsP.exportMPS(strcat(filename,"mpsGS.mat"));

   // Now get the V operator and apply a number of times
   MPS** states;
   //int nrCorr=5+1; // how many do I want
   { // I just don't need the MPOs after applying them)
     //states=initializeRandomStates(nrCorr,gs,D);
     // // I first apply sx (flip one) and then construct interpolants on that
     // MPS aux(gs);
     // contractor.optimize(sx,gs,aux);
     // MPO Vmpo(L);hSch.constructVMPO(Vmpo);
     // states=initializeVStates(nrCorr,Vmpo,aux,D);
     // // Interpolators from powers of Sz
     // states=initializeVStates(nrCorr,sz,gs,D);
     // // Interpolators from powers of Sx
     //  states=initializeVStates(nrCorr,sx,gs,D);
     MPO Vmpo(L);
     //  hSch.constructScalMPO(Vmpo);
     //  MPO theMpo(2);
     //  hSch.constructVecTwoSiteMPO(theMpo);
     //  MPO Vmpo(L);hSch.constructVecMPO(Vmpo);
     //  Vmpo.exportForMatlab("Ovec.m");
     //  states=initializeVStatesTR(nrCorr,Vmpo,gs,D,hSch);
     //  states=initializeVStatesAlt(nrCorr,Vmpo,gs,D,theMpo);
     hSch.constructVecMPOobc(Vmpo);
     states=initializeVStates(nrCorr,Vmpo,gs,D);
   }
   cout<<"All the states prepared "<<endl;

   // And now, evolve in imaginary time, and see

   int M;
   double t1=0; // current time
   M=ceil(2*t0/delta);
   mwArray Mmatr(Indices(nrCorr,nrCorr));
   constructCorrelatorMatrix(nrCorr,states,Mmatr);
   ofstream filetest("./primeraCorr.txt");
   Mmatr.savetext(filetest);
   filetest.close();
   //   *out<<k<<"\t"<<t1<<"\t";
   //printMatrixElements(*out,Mmatr);
   //*out<<endl;

   mwArray C0; // place for C0^-(1/2)
   
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

   char* name=new char[120];
   sprintf(name,"mpo_expHe2.m");
   expHxe_2.exportForMatlab(name);
   sprintf(name,"mpo_expHo.m");
   expHxo.exportForMatlab(name);
   //sprintf(name,"mpo_expHL.m");
   //expHL->exportForMatlab(name);
   sprintf(name,"mpo_expHL2.m");
   expHL_2.exportForMatlab(name);
   delete [] name;



#ifdef EXPEC // Computing also some expectation values
    // The operator I want to compute
   MPO vecMPO(L),scalMPO(L);
   hSch.constructVecMPO(vecMPO);
   hSch.constructScalMPO(scalMPO);
   // Other observables (spin?)
#endif // EXPEC // Computing also some expectation values

   bool imag=true;
   bool done=false;
   //bool imag=false;
   while(k<M&&!done){
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
    
    // cout<<"Not normalizing, so norms now: "<<contractor.contract(st1,st1)
    // 	<<", "<<contractor.contract(st2,st2)
    // 	<<", "<<contractor.contract(st3,st3)<<endl;
    // exit(1);

    k++;t1+=delta;

    constructCorrelatorMatrix(nrCorr,states,Mmatr);
    //*out<<k<<"\t"<<t1<<"\t";
    //cout<<"Printing matrix elements for time t1="<<t1<<endl;
    //printMatrixElements(*out,Mmatr);
    // cout<<"Printed matrix elements for time t1="<<t1<<endl;
    //*out<<endl;

    mwArray expecV,expecS; // place for the expectation matrix 
    mwArray overl; // place for the overlap matrix 
    mwArray expecSp; // place for <S^2>
    //    if(abs(t1-t0/2)<1E-10){
    if(abs(t1-t0)<1E-10){
      // cout<<"Exactly on reference time: t="<<t1<<" (t0="<<t0<<")"<<endl;
      // cout<<"Mmatr="<<Mmatr<<endl;
      constructReferenceMatrix(Mmatr,C0);
      // cout<<"And supposed reference matrix"<<endl;
      // cout<<C0<<endl;
      // exit(1);
    }
    //    else if(t1>t0/2){
    else if(t1>t0){
      if(t1>=2.*t0) done=true; // stop the simulation
      mwArray U;vector<complex_t> eigVs;
      // cout<<"Checking Hermiticity!!"<<endl;
      // cout<<"C0^(-1/2): "<<setprecision(10)<<C0<<endl;
      // cout<<"Mmatr: "<<Mmatr<<endl;
      // cout<<"C0*M="<<C0*Mmatr<<endl;
      // cout<<"M*C0="<<Mmatr*C0<<endl;
      mwArray matr(C0);matr.multiplyLeft(Mmatr);matr.multiplyLeft(C0);
      //      cout<<"C0*M*C0="<<setprecision(10)<<matr<<endl;
      matr=.5*(matr+Hconjugate(matr));
      //wrapper::eig(matr,eigVs,U,false); // not computing eigenvalues
      wrapper::eig(matr,eigVs,U,true); // computing eigenvalues
      // And now, check the expectation value of the vector operator
      // in each of them
#ifdef EXPEC // Computing also some expectation values
      constructOverlapMatrices(expecV,overl,vecMPO,states,nrCorr);
      constructOverlapMatrices(expecS,overl,scalMPO,states,nrCorr);
      constructSpinMatrices(expecSp,overl,sx,sy,sz,states,nrCorr);
      mwArray values=Hconjugate(U)*expecV*U;
      mwArray valuesS=Hconjugate(U)*expecS*U;
      mwArray valuesSp=Hconjugate(U)*expecSp*U;
      mwArray norms=Hconjugate(U)*overl*U;
      vector<double> theVals;values.getRealDiagonal(theVals);
      vector<double> theValsS;valuesS.getRealDiagonal(theValsS);
      vector<double> theValsSp;valuesSp.getRealDiagonal(theValsSp);
      vector<double> theNorms;norms.getRealDiagonal(theNorms);
#endif // EXPEC // Computing also some expectation values
      //cout<<k<<"\t"<<t1*2-t0<<"\t"<<eigVs<<endl;
      bool isNaN=false;
      *out2<<k<<"\t"<<t1<<"\t";
      for(int rr=0;rr<nrCorr;rr++){
	*out2<<real(eigVs[rr])<<"\t";
	if(std::isnan(real(eigVs[rr]))) isNaN=true;
#ifdef EXPEC // Computing also some expectation values
	*out2<<theVals[rr]/theNorms[rr]<<"\t";
	*out2<<theValsS[rr]/theNorms[rr]<<"\t";
	*out2<<theValsSp[rr]/theNorms[rr]<<"\t";
#endif // EXPEC // Computing also some expectation values
      }
      *out2<<endl;
      if(isNaN){
	cout<<"Not a Number in the diagonalization! Aborting (try shorter t0)"<<endl;
	exit(1);
      }
    }
   }

   deleteStates(states,nrCorr);
   //   out->close();
   //delete out;
   out2->close();
   delete out2;

}


  
void deleteStates(MPS** array,int nr){
  for(int k=0;k<nr;k++)
    delete array[k];
  delete [] array;
}

MPS** initializeVStates(int nr,const MPO& Vmpo,const MPS& gs,int D){
  MPS** array=new MPS*[nr]; // place for pointers to the states
  Contractor& contractor=Contractor::theContractor();
  for(int k=0;k<nr;k++){ // prepare, one by one
    cout<<"Preparing state V^"<<k<<"|GS>"<<endl;
    MPS* st=new MPS(gs,D); // first, initialize with gs
    // I need to apply Vmpo k+1 times
    const MPO** mpoPtr=new const MPO*[k];
    if(k>0){
      for(int l=0;l<k;l++){
	mpoPtr[l]=&Vmpo;
      }
      MPO Vkmpo(Vmpo.getLength());
      MPO::join(k,mpoPtr,Vkmpo); // applying twice
      contractor.optimize(Vkmpo,gs,*st,D);
      st->gaugeCond('R',1);
    }
    array[k]=st;
    delete [] mpoPtr;
    
  }
  return array;
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

void printMatrixElements(ofstream& out,const mwArray& Mmatr){
  int nrCorr=Mmatr.getDimension(0);
  for(int k1=0;k1<nrCorr;k1++)
    for(int k2=k1;k2<nrCorr;k2++)
      out<<Mmatr.getElement(Indices(k1,k2))<<"\t";
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
   
void constructOverlapMatrices(mwArray& expecV,mwArray& overl,const MPO& vecMPO,
			    MPS** states,int nrCorr){
  expecV=mwArray(Indices(nrCorr,nrCorr));
  overl=mwArray(Indices(nrCorr,nrCorr));
  Contractor& contractor=Contractor::theContractor();
  for(int k=0;k<nrCorr;k++){
    complex_t diagV=contractor.contract(*states[k],vecMPO,*states[k]);
    expecV.setElement(ONE_c*real(diagV),Indices(k,k));
    complex_t diag=contractor.contract(*states[k],*states[k]);
    overl.setElement(ONE_c*real(diag),Indices(k,k));
    for(int j=k+1;j<nrCorr;j++){
      complex_t offdiagV=contractor.contract(*states[j],vecMPO,*states[k]);
      expecV.setElement(offdiagV,Indices(k,j));
      expecV.setElement(conjugate(offdiagV),Indices(j,k));
      complex_t offdiag=contractor.contract(*states[j],*states[k]);
      overl.setElement(offdiag,Indices(k,j));
      overl.setElement(conjugate(offdiag),Indices(j,k));
    }
  }
}

void constructSpinMatrices(mwArray& expecSp,mwArray& overl,const MPO& sxMPO,
			   const MPO& syMPO,const MPO& szMPO,MPS** states,int nrCorr){
  expecSp=mwArray(Indices(nrCorr,nrCorr));
  overl=mwArray(Indices(nrCorr,nrCorr));
  Contractor& contractor=Contractor::theContractor();
  // Prepare three joined operators for sx^2, etc
  // const MPO* sxptr[2]={&sxMPO,&sxMPO};
  // const MPO* syptr[2]={&syMPO,&syMPO};
  // const MPO* szptr[2]={&szMPO,&szMPO};
  // MPO sx2(sxMPO.getLength()),sy2(syMPO.getLength()),sz2(szMPO.getLength());
  // MPO::join(2,sxptr,sx2);
  // MPO::join(2,syptr,sy2);
  // MPO::join(2,szptr,sz2);
  for(int k=0;k<nrCorr;k++){
    //    complex_t diagVx=contractor.contract(*states[k],sx2,*states[k]);
    // complex_t diagVy=contractor.contract(*states[k],sy2,*states[k]);
    //    complex_t
    // diagVz=contractor.contract(*states[k],sz2,*states[k]);
    //complex_t diagV=diagVx+diagVy+diagVz;
    complex_t diagV=contractor.contract(*states[k],szMPO,*states[k]);
    expecSp.setElement(ONE_c*real(diagV),Indices(k,k));
    complex_t diag=contractor.contract(*states[k],*states[k]);
    overl.setElement(ONE_c*real(diag),Indices(k,k));
    for(int j=k+1;j<nrCorr;j++){
      // complex_t offdiagVx=contractor.contract(*states[j],sx2,*states[k]);
      // complex_t offdiagVy=contractor.contract(*states[j],sy2,*states[k]);
      // complex_t
      // offdiagV=contractor.contract(*states[j],sz2,*states[k]);
      // complex_t value=offdiagVx+offdiagVy+offdiagVz;
      complex_t offdiagV=contractor.contract(*states[j],szMPO,*states[k]);
      expecSp.setElement(offdiagV,Indices(k,j));
      expecSp.setElement(conjugate(offdiagV),Indices(j,k));
      complex_t offdiag=contractor.contract(*states[j],*states[k]);
      overl.setElement(offdiag,Indices(k,j));
      overl.setElement(conjugate(offdiag),Indices(j,k));
    }
  }
}

MPS** initializeRandomStates(int nr,const MPS& gs,int D){
  MPS** array=new MPS*[nr]; // place for pointers to the states
  for(int k=0;k<nr;k++){ // prepare, one by one
    MPS* st=new MPS(gs,D); // first, initialize with gs
    if(k>0)
      st->setRandomState();
    st->gaugeCond('R',1);
    array[k]=st;
    // Now, instead of normalizing each state to one, I will normalize such that the overlap is larger
    if(k>0){
      //      MPS::gaugeCond(*array[k],*array[k-1],'R');
      MPS::gaugeCond(*array[k],*array[0],'R');
    }
    // char* name=new char[120];
    // sprintf(name,"mps_n%d.m",k+1);
    // st->exportForMatlab(name);
    // delete [] name;
  }
  return array;
}


MPS** initializeVStatesTR(int nr,const MPO& Vmpo,const MPS& gs,int D,
			  SchwingerHamiltonian& hSch){
  MPS** array=new MPS*[nr]; // place for pointers to the states
  MPO trampa(Vmpo.getLength()); //TR
  hSch.constructVecMPO(trampa); //TR
  MPO trampa2(Vmpo.getLength()); //TR
  hSch.constructScalMPO(trampa2); //TR
  Contractor& contractor=Contractor::theContractor();
  for(int k=0;k<nr;k++){ // prepare, one by one
    cout<<"Preparing state V^"<<k+1<<"|GS>"<<endl;
    MPS* st=new MPS(gs,D); // first, initialize with gs
    // I need to apply Vmpo k+1 times
    const MPO** mpoPtr=new const MPO*[k+1];
    for(int l=0;l<k+1;l++){
      mpoPtr[l]=&Vmpo;
      if(l%3==0)mpoPtr[l]=&trampa; //TR
      if(l%3==1)mpoPtr[l]=&trampa2; //TR
    }
    MPO Vkmpo(Vmpo.getLength());
    MPO::join(k+1,mpoPtr,Vkmpo); // applying twice
    contractor.optimize(Vkmpo,gs,*st,D);
    st->gaugeCond('R',1);
    array[k]=st;
    delete [] mpoPtr;

  }
  return array;
}

/** Applies a two body operator somewhere in the chain!
NO TI! */

MPS** initializeVStatesAlt(int nr,const MPO& Vmpo,const MPS& gs,int D,const MPO& theMpo){
  int L=gs.getLength();
  MPS** array=new MPS*[nr]; // place for pointers to the states
  MPO trampa(Vmpo.getLength()); //TR
  mwArray id=identityMatrix(2); //TR
  id.reshape(Indices(2,1,2,1));//TR
  Operator idOp(id);//TR
  for(int k=0;k<Vmpo.getLength();k++)//TR
    trampa.setOp(k,&idOp,false);//TR
  Contractor& contractor=Contractor::theContractor();
  for(int k=0;k<nr;k++){ // prepare, one by one
    cout<<"Preparing state V^"<<k<<"|GS>"<<endl;
    MPS* st=new MPS(gs,D); // first, initialize with gs
    if(k>0){
      if(2*k+1<gs.getLength()){
	// I need to apply some special Vmpo 
	trampa.setOp(2*k,&theMpo.getOp(0),false);
	trampa.setOp(2*k+1,&theMpo.getOp(1),false);
	contractor.optimize(trampa,gs,*st,D);
	st->gaugeCond('R',1);
	trampa.setOp(2*k,&idOp,false);
	trampa.setOp(2*k+1,&idOp,false);
      }
      else{ // I have to do something more
	int kp=2*k%L; // pos mod L
	trampa.setOp(kp,&theMpo.getOp(0),false);
	trampa.setOp(kp+1,&theMpo.getOp(1),false);
	trampa.setOp((kp+L/2)%L,&theMpo.getOp(0),false);
	trampa.setOp((kp+1+L/2)%L,&theMpo.getOp(1),false);
	contractor.optimize(trampa,gs,*st,D);
	st->gaugeCond('R',1);
	trampa.setOp(kp,&idOp,false);
	trampa.setOp(kp+1,&idOp,false);
	trampa.setOp((kp+L/2)%L,&idOp,false);
	trampa.setOp((kp+1+L/2)%L,&idOp,false);      
      }
    }
    array[k]=st;
  }
  return array;
}
