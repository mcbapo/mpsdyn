#include "SchwingerHamiltonianSz.h"
#include "Indices.h"
#include "Contractor.h"

using namespace std;
using namespace shrt;

SchwingerHamiltonianSz::SchwingerHamiltonianSz(int N_,double mu_,double x_,
					       double l0_,double alpha_,
					       double zpen_,
					       double offset_,double nu_,
					       double y_,double scale_,double Stot_):
  zpen(zpen_),scale(scale_),Stot(Stot_){
  N=N_;mu=mu_;x=x_;l0=l0_;alpha=alpha_;nu=nu_;gweight=y_;zpen=zpen_;
  offset=offset_;
  scale=scale_;
  if(N%2!=(int(2*Stot+.5-(Stot<0)))%2){
    cout<<"WARNING: Invalid value of total spin "<<Stot<<" for a chain of length "<<N<<", needs to be fixed before optimizing"<<endl;
  }
  hamil.initLength(N);
  // Construct the MPO
  initZ();
  initHMPO(); // my own
}

void SchwingerHamiltonianSz::setStarget(double Stot_){
  Stot=Stot_;
  if(N%2!=abs((int(2*Stot+.5-(Stot<0)))%2)){
    cout<<"ERROR: Invalid value of total spin "<<Stot<<" for a chain of length "<<N<<"(N%2="<<N%2<<", Stot"<<Stot<<" int(2*Stot+.5-(Stot<0))="<<(int(2*Stot+.5-(Stot<0)))<<")"<<endl;
    exit(1);
  }
  initHMPO();
}

void SchwingerHamiltonianSz::initHMPO(){
  int D=5; // bond dimension of the MPO
  int d=2; // physical dimension
  double alpha_= alpha + l0;
  // If a scale factor is applied, the factor in from of every
  // Operator
  double factor=scale!=1?pow(scale,1/((double)N)):1.;

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
      C.setElement(gweight*alphak*alphak+(complex_t){gweight*(k+1)/4.+mu/2.+nu/2.+(double)offset/(double)N+zpen*(1.+4*Stot*Stot/N),0.},
		   Indices(0,Dr-1,0)); // fk' 1
      C.setElement((complex_t){signok*mu/2.+nu/2.-4*zpen*Stot,0.}+
		   gweight*sumA.getElement(Indices(0,k)),Indices(0,Dr-1,3));
    }
    else{
      C.setElement((complex_t){mu/2.+nu/2.+(double)offset/(double)N+zpen*(1+4*Stot*Stot/N),0.},Indices(0,Dr-1,0)); // fk' 1
      C.setElement((complex_t){signok*mu/2.+nu/2.-4*zpen*Stot,0.},Indices(0,Dr-1,3));
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
    // Create and store an Operator
    hamil.setOp(k,new Operator(permute(factor*res,Indices(3,1,4,2))),true);
  }

}
