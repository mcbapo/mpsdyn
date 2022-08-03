#include "SchwingerHamiltonianSzSR.h"
#include "Indices.h"
#include "Contractor.h"

using namespace std;
using namespace shrt;

SchwingerHamiltonianSzSR::SchwingerHamiltonianSzSR(int N_,double mu_,double x_,
					       double l0_,double alpha_,
						   double zpen_,double lambda_,
					       double offset_,double nu_,
						   double y_):
  zpen(zpen_),lambda(lambda_){
  N=N_;mu=mu_;x=x_;l0=l0_;alpha=alpha_;nu=nu_;gweight=y_;
  offset=offset_;
  hamil.initLength(N);
  // Construct the MPO
  initZ();
  initHMPO(); // my own
}


void SchwingerHamiltonianSzSR::initHMPO(){
  int D=5; // bond dimension of the MPO without the extra term
  int d=2; // physical dimension
  double alpha_= alpha + l0;
  // Auxiliary: the values that depend on the site
  // even-odd alpha
  mwArray alphas(Indices(1,N));
  for(int k=0;k<N;k++){
    double signo=(k+1)%2==0?+1.:-1.; // since this k starts in 0
    alphas.setElement((complex_t){alpha_-.25*(1.-signo)},Indices(0,k));
  }

  double signZpen=zpen>0?1:-1;
  // sums
  mwArray sumA(Indices(1,N));sumA.fillWithZero();
  complex_t sumTot=ZERO_c;
  for(int k=N-2;k>=0;k--){
    sumTot=sumTot+alphas.getElement(Indices(0,k));
    sumA.setElement(sumTot,Indices(0,k));
  }
  // The terms for SR and SR+ are constructed by SchwingerHamiltonian::constructShiftRotateMPO
  MPO Sr(N),SrD(N);
  constructShiftRotateMPO(Sr);
  //  constructShiftRotateMPO(SrD,false);
  // Distribute the lambda factor among all terms in the MPO
  double lambdaF=pow(.5*lambda,1./N);
  cout<<"In initHMPO, lambdaF="<<lambdaF<<", lambda="<<lambda<<" and lambda/N="<<lambda/N<<", while offset="<<offset<<endl;
  // in H, I see -(lambda/2)( Sr + Sr+ ), so on each term I can put (lambda/2)^(1/N) and on the first the minus sign

  // \TODO To optimize the performance I should only change the C elements
  // when calling the function several times
  for(int k=0;k<N;k++){
    // the MPS containing all the coefficients. But they will depend on the
    // position
    int Dl=D,Dr=D;
    if(k==0) Dl=1;
    if(k==N-1) Dr=1;
    mwArray C(Indices(Dl,Dr,4));C.fillWithZero();
    // identities
    if(k!=N-1) // all but the last one can pass on a D=0
      C.setElement(ONE_c,Indices(0,0,0));

    if(k!=0) // in the first one I canot put 1 and finish
      C.setElement(ONE_c,Indices(Dl-1,Dr-1,0));

    //if((k>0)&&(k<N-2)) // all but the last and first to last can be intermediate
    if((k>0)&&(k<N-1))   // with the penalty term, the first to last can also pass 3
      C.setElement(ONE_c,Indices(3,3,0));

    complex_t alphak=alphas.getElement(Indices(0,k));
    double signok=(k+1)%2==0?+1.:-1.;
    if(k!=N-1){
      C.setElement(gweight*alphak*alphak+(complex_t){gweight*(k+1)/4.+mu/2.+nu/2.+(double)offset/(double)N+zpen+(double)lambda/(double)N,0.},
		   Indices(0,Dr-1,0)); // fk' 1
      C.setElement((complex_t){signok*mu/2.+nu/2.,0.}+
		   gweight*sumA.getElement(Indices(0,k)),Indices(0,Dr-1,3));
    }
    else{
      C.setElement((complex_t){mu/2.+nu/2.+(double)offset/(double)N+zpen+(double)lambda/(double)N,0.},Indices(0,Dr-1,0)); // fk' 1
      C.setElement((complex_t){signok*mu/2.+nu/2.,0.},Indices(0,Dr-1,3));
    }
    //C(1,Dr,4)=(mu/2)*(-1)^k + sum(alphas(k:N)); // fk Z

    if(k!=(N-1)){
      C.setElement((complex_t){sqrt(x/2),0.},Indices(0,1,1)); // XX
      C.setElement((complex_t){sqrt(x/2),0.},Indices(0,2,2)); // YY
      //      if(k<N-2)
      C.setElement(ONE_c,Indices(0,3,3)); // k ZZ
      if(k!=0){
	C.setElement((complex_t){sqrt(x/2),0.},Indices(1,4,1)); // XX
	C.setElement((complex_t){sqrt(x/2),0.},Indices(2,4,2)); // YY
	C.setElement((complex_t){gweight*(N-k-1)/2.+2*zpen,0.},Indices(3,4,3)); // N-k ZZ
	//if(k<(N-2)) C.setElement(ONE_c,Indices(3,3,0)); // 
      }
    }
    else{ // in the last site I can only have the second half
	C.setElement((complex_t){sqrt(x/2),0.},Indices(1,0,1)); // XX
	C.setElement((complex_t){sqrt(x/2),0.},Indices(2,0,2)); // YY
	C.setElement((complex_t){2*zpen,0.},Indices(3,0,3)); //  ZZ only from penalty term
	//C.setElement((complex_t){(N-k)/2.,0.},Indices(3,0,3)); // N-k+1 ZZ
    }

    // Now contract MPS and operators and set proper index ordering
    mwArray res=reshape(reshape(C,Indices(Dl*Dr,4))*
			reshape(Z,Indices(4,d*d)),Indices(Dl,Dr,d,d));
    res.permute(Indices(3,1,4,2));
    // Now I expand it to include the corresponding lambda(SR + SR+) term
    mwArray extraSr=Sr.getOp(k).getFullData();
    //mwArray extraSrD=SrD.getOp(k).getFullData();
#ifdef CHECKDIMS
    if(//(extraSr.getDimensions()!=extraSrD.getDimensions())||
       (extraSr.getDimension(0)!=d)
       ||(extraSr.getDimension(2)!=d)){
      cout<<"Error when trying to patch the MPO for H with Sr+Sr+"<<endl;
      exit(212);
    }
#endif
    int extraDimL(2), extraDimR(2);int signL=1;
    if(k==0){extraDimL=0;signL=-1;}
    if(k==N-1) extraDimR=0;
    res.resize(Indices(d,Dl+2*extraDimL,d,Dr+2*extraDimR));
    // fill in the extra elements
    for(int i=0;i<d;i++)
      for(int j=0;j<d;j++){
	if(k!=0&&k!=N-1){
	  for(int l=0;l<extraDimL;l++)
	    for(int r=0;r<extraDimR;r++){
	      res.setElement(signL*lambdaF*extraSr.getElement(Indices(i,l,j,r)),Indices(i,Dl+l,j,Dr+r));
	      res.setElement(conjugate(signL*lambdaF*extraSr.getElement(Indices(j,l,i,r))),Indices(i,Dl+extraDimL+l,j,Dr+extraDimR+r));
	      //	    res.setElement(signL*lambdaF*extraSrD.getElement(Indices(i,l,j,r)),Indices(i,Dl+extraDimL+l,j,Dr+extraDimR+r));
	    }
	}
	else if(k==0){
	  int l=0;
	  for(int r=0;r<extraDimR;r++){
	    res.setElement(signL*lambdaF*extraSr.getElement(Indices(i,l,j,r)),Indices(i,Dl-1,j,Dr+r));
	    res.setElement(conjugate(signL*lambdaF*extraSr.getElement(Indices(j,l,i,r))),Indices(i,Dl-1,j,Dr+extraDimR+r));
	  }
	}
	else {//k==N-1
	  int r=0;
	  for(int l=0;l<extraDimL;l++){
	    res.setElement(signL*lambdaF*extraSr.getElement(Indices(i,l,j,r)),Indices(i,Dl+l,j,Dr-1));
	    res.setElement(conjugate(signL*lambdaF*extraSr.getElement(Indices(j,l,i,r))),Indices(i,Dl+extraDimL+l,j,Dr-1));
	  }
	}
      }
    // Create and store an Operator
    hamil.setOp(k,new Operator(res),true);
  }

}
