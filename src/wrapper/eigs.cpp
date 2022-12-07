#include "wrapper.h"
#include "mwArray.h"
#include "Indices.h"
#include <cstdio>
#include <vector>
#include "BasicMultiplier.h"

using namespace std;
using namespace shrt;

#ifdef PRIMMEv2
#define PRIMME_COMPLEX PRIMME_COMPLEX_DOUBLE
#else
#define PRIMME_COMPLEX Complex_Z
#endif


#define NUMPREC 2.22E-16

extern "C" {
  /** ARPACK routine. Used to compute a few eigenvalues of a certain
      linear operator. */ 
  void znaupd_(int* ido,const char* bmat,const int* n,const char* which,
	       int* nev,double* tol,double* RESID,int* ncv,
	       double* V,int* ldv,int* IPARAM,int* IPNTR,
	       double* WORKD,double* WORKL,int* lworkl,
	       double* RWORK,int* info);
  /** ARPACK routine which extracts eigenvalues and eigenvectors from
      the results of znaupd. See documentation:
      \url<http://www.mathkeisan.com/UsersGuide/man/zneupd.html> */
  void zneupd_(int* rvec,char* howmny,int* SELECT,double* D,
	       double* Z,int* ldz,double* sigma,double* workev,
	       const char* bmat,const int *N, const char *which, 
	       int *nev, double *tol,
	       double* RESID,int* ncv,double* V,int* ldv,int* IPARAM,
	       int* IPNTR,double* WORKD,double* WORKL,int* lworkl,
	       double* RWORK,int* info);
}


int wrapper::eigs(const mwArray& A,const int k,const char* which,
		   std::vector<complex_t>& D,mwArray& U,bool eigv,double tol){

  // The basic callback would be the actual product, which I will
  // implement as a "nested" function and then use to call the other
  // flavour of eigs.
  BasicMultiplier aux(A);
  return wrapper::eigs(aux,k,which,D,U,eigv,tol);

}


/*
#include "time.h"
timespec diff(timespec start, timespec end)
{
  timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}
*/

#define MAXREC 10

int wrapper::eigs(Multiplier& multi,const int k,
		   const char* which,std::vector<complex_t>& Dval,mwArray& U,
		   bool eigv,double tol){
  int N=multi.getSize();
  int ido=0; // reverse communication flag, init 0
  char bmat[2]="I"; // Only regular eigenvalue problem (not generalized)
  int nev=k;
  static int recursive=0; // count of recursive calls
  //double tol=0.;
  //double tol=1E-10;
  if(tol==0.){
#ifndef LOWTOL
    tol=NUMPREC;
#else
    tol=1E-6;
#endif
  }
  double* resid=new double[2*N];
  //int ncv=min(N,nev+2); //4*nev; 
  int ncv=min(N,max(20,3*nev)); //4*nev; 
  int ldv=N;
  double* V=new double[2*(ldv*ncv)];
  int* iparam=new int[11];
  iparam[0]=1; // not sure what that means
  iparam[2]=max(300,(int)(2*N/ncv)+1);//3*N;
  iparam[6]=1;
  int* ipntr=new int[14];
  double* workd=new double[2*3*N];
  int lworkl=3*ncv*ncv+5*ncv;
  double* workl=new double[2*lworkl];
  double* rwork=new double[ncv];
  int info=0;
  // Repeatedly call the Arpack routine 
  bool done=0;
  if(!U.isEmpty()){ // use an initial vector passed by the caller
    info=1;
    if(U.getDimension(0)!=N){
      cout<<"Error: wrong dimension of initial vector for eigs!"<<endl;
      exit(212);
    }
    for(int ik=0;ik<N;ik++)
      U.getElement(resid[2*ik],resid[2*ik+1],shrt::Indices(ik,0));
  }
  while(!done){
    //clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time1);//start=clock();
    znaupd_(&ido, bmat, &N, which, &nev, &tol, resid, 
	    &ncv, V, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, rwork, &info);
    //countCalls++;clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time2);//finish=clock();
    //cout<<"Called znaupd with args: ido="<<ido<<",bmat="<<bmat<<",N="<<N
    //	<<",which="<<which<<",nev="<<nev<<",tol="<<tol<<",ncv="
    //	<<ncv<<",ldv="<<ldv<<",lworkl="<<lworkl<<",time="
    //	<<diff(time1,time2).tv_nsec*1.E-6;
    ////<<(finish-start)*1./(countCalls*CLOCKS_PER_SEC)<<endl;
    if ((ido==1)||(ido==-1)){ // need to multiply
      // improve efficiency by avoiding copy
      mwArray X;X.setPointer(Indices(N,1),&workd[2*(ipntr[0]-1)]);
      mwArray Y;Y.setPointer(Indices(N,1),&workd[2*(ipntr[1]-1)]);
      //X.setData(N,(complex_t*)&workd[2*(ipntr[0]-1)]); // copied
      //Y=X;
      //clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time2);
      multi.product(X,Y);
      //clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time3);
      //cout<<", time in product="<<diff(time2,time3).tv_nsec*1.E-6<<endl;
      // copy the result CHECK!!!
      //memcpy((void*)&workd[2*(ipntr[1]-1)],(void*)Y.components,
      //     sizeof(double)*2*N);
    }
    done=((ido!=1)&&(ido!=-1));
  }
  //finish=clock();
  // Finished iteration: rearrange results, if ok
  if(info<0){
    cout<<"Error within znaupd routine, info ="<<info<<endl;
    cout<<"See documentation in "
	<<"http://www.mathkeisan.com/UsersGuide/man/znaupd.html"<<endl;
    exit(212);
  } 
  else{
    //cout<<"Out of znaupd, after "<<countCalls<<" calls, avg. time= "
    //	<<(finish-start)*1./(countCalls*CLOCKS_PER_SEC)<<endl;
    int revec=eigv?1:0; //whether we have to compute eigenvalues
    char howmny[2]="A";
    int* select = new int[ncv];    
    double* D=new double[2*(nev+1)]; // For eigenvalues
    double* Z=revec?new double[2*N*nev]:0;
    int ldz=revec?N:1;
    double sigma;
    double* workev=new double[2*3*ncv];
    zneupd_(&revec,howmny, select, D, Z, &ldz, &sigma, workev, 
	    bmat, &N, which, &nev, &tol, resid, &ncv, V, &ldv,
	    iparam, ipntr, workd, workl, &lworkl, rwork, &info);
    if (info!=0) {
      if(info==-14){ //REMOVE THIS BLOCK
	if(recursive>MAXREC||k>=N){ // Abort
	  cout<<"Got the -14 error, more than max("<<MAXREC<<") times, with tol="<<tol
	      <<"but meanwhile, the ("<<k<<")computed values were"
	      <<endl;
	  for(int kn=0;kn<nev;kn++){
	    cout<<(complex_t){D[2*kn],D[2*kn+1]}<<",";
	  }
	  U.clear();
	  if(revec){
	    U=mwArray(2,Indices(N,nev));// Q: Which constructor is this line calling??? (int rank,const int* dims)???
	    U.setData(N*nev,(complex_t*)&Z[0]);
	    //cout<<"), vectors "<<U;
	  }
	  cout<<endl;
	  return info; //exit(212);
	}
	else{
	  cout<<"zneupd error:-14. I am trying again with larger tolerance, as, according to MAtlab technical note, "
	      <<"the error comes from bad conditioned matrix"<<endl;
	  // try again!
	  // // cout<<"zneupd error:-14. I am trying again for twice as many eigenvalues. Now: ";
	  // // for(int kn=0;kn<nev;kn++){
	  // //   cout<<(complex_t){D[2*kn],D[2*kn+1]}<<",";
	  // // }
	  // // cout<<endl;
	  // Clear up!
	  delete[] select;
	  delete[] D;
	  if(revec) delete[] Z;
	  delete[] workev;    
	  delete[] V;
	  recursive++;
	  // //	  return eigs(multi,2*k,which,Dval,U,eigv);
	  tol=1E-6;
	  return eigs(multi,2*k,which,Dval,U,eigv,tol);
	  //	exit(212);
	} //REMOVE THIS BLOCK
      }
      cout << "Error with zneupd, info = " << info <<endl;
      cout<<"See documentation in "
	  <<"[http://www.mathkeisan.com/UsersGuide/man/zneupd.html]"
	  <<endl;
      exit(212); 
    } 
    else if(info ==1){
      cout << "Maximum number of iterations reached.\n\n";
    } 
    else if(info ==3){
      cout << "No shifts could be applied during implicit\n";
      cout << "Arnoldi update, try increasing NCV.\n\n";
    }
    // Now rearrange results as expected and exit
    Dval.clear();
    for(int k=0;k<nev;k++){
      Dval.push_back((complex_t){D[2*k],D[2*k+1]});
    }
    U.clear();
    if(revec){
      U=mwArray(2,Indices(N,nev));
      U.setData(N*nev,(complex_t*)&Z[0]);
    }
    // clean up here
    delete[] select;
    delete[] D;
    if(revec) delete[] Z;
    delete[] workev;    
    delete[] V;
  }

  // Clean up
  delete[] resid;
  delete[] iparam;
  delete[] ipntr;
  delete[] workd;
  delete[] workl;
  delete[] rwork;

  recursive=0; // clear count of recursion
  return info;
}

#ifdef USING_PRIMME

/** To adapt to the interface of PRIMME eigensolver, we need to
    specify a callback for the matrix-vector product. So, we need to
    wrap the ususal Multiplier member function in the type that PRIMME
    wants. Moreover, since we would need to call a member function
    (for every different Multiplier implementation), and we cannot
    pass a pointer to that, we "hide" the pointer to the Multiplier as
    one member of an "extended" primme_params struct, that is only
    visible here and used only by eigs and the matrix-vector product
    function. */

// The extended primme_params structure

struct primme_params_ext {
  primme_params primme;
  Multiplier* multi;
  Multiplier* multiB;
};

#ifndef PRIMMEv2
// The method with the interface that primme wanted
void matrixMatvec(void* x,void* y,int* blockSize,primme_params* primme){
  // I receive a primme_params_ext and call the product function of the enclosed Multiplier
  // A bit of an ugly hack, here :(
  *blockSize=1;
  //primme_params_ext* primmePl=(primme_params_ext*)primme;
  //(primmePl->multi)->*product_primme(x,y);
  // fPtr thePtr=&Multiplier::product_primme;
  // (primmePl->multi->*thePtr)(x,y);
  //cout<<"In matrixMatvec(x,y)"<<endl;
  Multiplier* multi=((primme_params_ext*)primme)->multi;
  multi->product_primme(x,y);
}

// The method with the interface that primme wanted
void massMatrixMatvec(void* x,void* y,int* blockSize,primme_params* primme){
  // I receive a primme_params_ext and call the product function of the enclosed Multiplier
  // A bit of an ugly hack, here :(
  *blockSize=1;
  //primme_params_ext* primmePl=(primme_params_ext*)primme;
  //(primmePl->multi)->*product_primme(x,y);
  // fPtr thePtr=&Multiplier::product_primme;
  // (primmePl->multi->*thePtr)(x,y);
  Multiplier* multiB=((primme_params_ext*)primme)->multiB;
  multiB->product_primme(x,y);
}
#else
// Primme v2 requires different interface for matrixMatvec!!

void matrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
		     primme_params *primme, int *ierr){
  // I receive a primme_params_ext and call the product function of the enclosed Multiplier
  *blockSize=1;
  Multiplier* multi=((primme_params_ext*)primme)->multi;
  multi->product_primme(x,y);
}

void massMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
			 primme_params *primme, int *ierr){
  *blockSize=1;
  Multiplier* multiB=((primme_params_ext*)primme)->multiB;
  multiB->product_primme(x,y);
}

#endif //PRIMMEv2



int wrapper::eigs_primme(const mwArray& A,const int k,primme_target which,
			 std::vector<complex_t>& D,mwArray& U,double tol,primme_preset_method method){

  // The basic callback would be the actual product, which I will
  // implement as a "nested" function and then use to call the other
  // flavour of eigs.
  BasicMultiplier aux(A);
  return wrapper::eigs_primme(aux,k,which,D,U,tol,method);
}

int wrapper::eigs_primme(Multiplier& multi,const int nev,
			 primme_target which,std::vector<complex_t>& Dval,mwArray& U,
			 double tol,primme_preset_method method,int numShifts,double* targetShifts){

  int N=multi.getSize();
  if(nev>N){
    cout<<"Error: required number of eigenvalues eigenvalues "<<nev
	<<" larger than dimension "<<N<<endl;
    exit(1);
  }
  primme_params_ext primmePlus;
  primme_params& primme=primmePlus.primme;
  primme_initialize(&primme); // default initialization
  primme.n=N;
  primme.nLocal=N;
  primme.matrixMatvec=matrixMatvec;
  primmePlus.multi=&multi;
  primme.numEvals=nev;


  if(tol!=0) primme.eps=tol; // Q: Is this the right number here?
  primme.printLevel=0;
  //primme.printLevel=0; // For no output
  //primme.printLevel=2; // (3,4,5) For (increasingly) more information along the process 

  // Which eigenvalues to compute
  if((which!=primme_smallest)&&(which!=primme_largest)){
    primme.numTargetShifts=numShifts;
    primme.targetShifts=targetShifts;
  }
  primme.target=which;

  primme_set_method(method,&primme); 
  //  primme_set_method(DYNAMIC,&primme); // Q: Is this the best thing?
  //  primme_set_method(JDQMR,&primme); // Q: Is this the best thing? => Seems better than DYNAMIC for (H-E)^2
  // primme_set_method(RQI,&primme); // This seems quite bad

  int ret=wrapper::call_eigs_primme(&primme,Dval,U);
#ifdef PRIMMEv2
  primme_free(&primme);
#else
  primme_Free(&primme);
#endif
  return ret;
}

int wrapper::call_eigs_primme(primme_params* primme,std::vector<complex_t>& Dval,mwArray& U){
  int nev=primme->numEvals;
  int N=primme->n;
  // Place for results, as PRIMME likes it
  double *evals=new double[2*nev]; // real eigenvalues
  PRIMME_COMPLEX* evecs=new PRIMME_COMPLEX[N*2*nev];
  double* rnorms=new double[2*nev];
  //cout<<"In call_eigs_primme(U="<<U<<")"<<endl;
  //primme_display_params(*primme);


  // If U was passed with something inside, it is used as intial vector(s ?)  
  if(!U.isEmpty()){
    for(int ik=0;ik<N;ik++){
#ifndef PRIMMEv2
      U.getElement(evecs[ik].r,evecs[ik].i,shrt::Indices(ik,0));
#else
      double re(0.),im(0.);
      U.getElement(re,im,shrt::Indices(ik,0));
      evecs[ik].real(re);evecs[ik].imag(im);
      primme->initSize=1;
#endif
    }
  }
  // cout<<"Calling eigs_primme with init U="<<U<<endl;

  // Call the solver itself
  int ret = zprimme(evals, evecs, rnorms, primme);
  
  // If some error occurred, do something // TODO!!
  if (ret != 0) {
    cout<<"Error: zprimme returned with nonzero exit status "<<ret<<endl;
    primme_display_params(*primme);
    return ret;
  }
  
  // Reformat/copy the results for the output of eigs
  // 1) Copy eigenvalues
  Dval.clear();
  for(int k=0;k<nev;k++){
    Dval.push_back(evals[k]*ONE_c);
  }
  // 2) Copy eigenvectors
  U.clear();
  U=mwArray(Indices(N,nev));
  U.setData(N*nev,(complex_t*)&evecs[0]);
  
  delete []evals;
  delete []evecs;
  delete []rnorms;
  return 0; // Successful exit
}

int wrapper::eigsgen_primme(Multiplier& multiA,Multiplier& multiB,const int nev,
			    primme_target which,std::vector<complex_t>& Dval,
			    mwArray& U,
			    double tol,
			    primme_preset_method method,
			    int numShifts,double* targetShifts){
  int N=multiA.getSize();
  if(nev>N){
    cout<<"Error: required number of eigenvalues eigenvalues "<<nev
	<<" larger than dimension "<<N<<endl;
    exit(1);
  }
  primme_params_ext primmePlus;
  primme_params& primme=primmePlus.primme;
  primme_initialize(&primme); // default initialization
  primme.n=N;
  primme.nLocal=N;
  primme.matrixMatvec=matrixMatvec;
  primmePlus.multi=&multiA;
  primmePlus.multiB=&multiB;
  primme.numEvals=nev;


  if(tol!=0) primme.eps=tol; // Q: Is this the right number here?
  //primme.printLevel=3;
  //primme.printLevel=0; // For no output
  //primme.printLevel=2; // (3,4,5) For (increasingly) more information along the process 

  // Which eigenvalues to compute
  if((which!=primme_smallest)&&(which!=primme_largest)){
    primme.numTargetShifts=numShifts;
    primme.targetShifts=targetShifts;
  }
  primme.target=which;

  primme_set_method(method,&primme); 
  //  primme_set_method(DYNAMIC,&primme); // Q: Is this the best thing?
  //  primme_set_method(JDQMR,&primme); // Q: Is this the best thing? => Seems better than DYNAMIC for (H-E)^2
  // primme_set_method(RQI,&primme); // This seems quite bad

  int ret=wrapper::call_eigs_primme(&primme,Dval,U);
#ifndef PRIMMEv2
  primme_Free(&primme);
#else
  primme_free(&primme);
#endif

  return ret;
}


#endif
