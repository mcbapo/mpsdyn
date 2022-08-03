
#include "CoherentDissipationLiouvillianLR.h"
#include "DoubleOperator.h"

using namespace std;
using namespace shrt;


CoherentDissipationLiouvillianLR::CoherentDissipationLiouvillianLR(int N_,double g_,double gamma_):
  CoherentDissipationLiouvillian(N_,g_,gamma_){
  isInit=false;
  initL();
  initLAdjoint();
  initTrL();
}

CoherentDissipationLiouvillianLR::~CoherentDissipationLiouvillianLR(){
  clear();
}

// void CoherentDissipationLiouvillianLR::clear(){
//   Lmpo.clear();
//   isInit=false;
//   LAdjmpo.clear();
// }

// void CoherentDissipationLiouvillianLR::initLAdjoint(){
//   // I can only do this if the Lmpo is already prepared
//   initL();

//   // I will create copies of the edge and middle operators with the
//   // adequate rotation of indices (up<->down) and complex conjugation
//   // and for the rest of the chain, I just point to the operator in
//   // the second site, as L itself does

//   LAdjmpo.setRotatedOp(0,&Lmpo.getOp(0),Indices(3,2,1,4),true);
//   LAdjmpo.setRotatedOp(N-1,&Lmpo.getOp(N-1),Indices(3,2,1,4),true);
//   if(N>2){
//     LAdjmpo.setRotatedOp(1,&Lmpo.getOp(1),Indices(3,2,1,4),true);
//     for(int k=2;k<N-1;k++)
//       LAdjmpo.setOp(k,&LAdjmpo.getOp(1),false);            
//   }
// }



void CoherentDissipationLiouvillianLR::initL(){
  if(isInit) return;
  cout<<"CoherentDissipationLiouvillianLR::initL"<<endl;
  int D=6; // bond dimension of the MPO
  mwArray Z; // The operators
  initZ(Z);
  int nrOps=Z.getDimension(0); // Nr of operators (10)
  // Now construct the coefficients for every term. I need different
  // ones for the first and last sites, but the rest are all the same
  mwArray C1(Indices(1,D,nrOps)); // the first site
  mwArray C(Indices(D,D,nrOps)); // the middle site
  mwArray CN(Indices(D,1,nrOps)); // the last site

  // Identity terms when nothing has happened yet
  C.setElement(ONE_c,Indices(0,0,0));
  C1.setElement(ONE_c,Indices(0,0,0));
  // Identity terms after the end
  C.setElement(ONE_c,Indices(D-1,D-1,0));
  CN.setElement(ONE_c,Indices(D-1,0,0));
  // Identity terms in between operators (all but edges)
  C.setElement(ONE_c,Indices(1,1,0));
  C.setElement(ONE_c,Indices(2,2,0));
  C.setElement(ONE_c,Indices(3,3,0));
  C.setElement(ONE_c,Indices(4,4,0));
  // Single body terms (one for each operator, with proper weights)
  C.setElement(-g*I_c,Indices(0,D-1,1));
  C.setElement(g*I_c,Indices(0,D-1,2));
  C.setElement(-gamma*.5*ONE_c,Indices(0,D-1,3));
  C.setElement(-gamma*.5*ONE_c,Indices(0,D-1,4));
  C.setElement(gamma*ONE_c,Indices(0,D-1,5));
  // The same for first site 
  C1.setElement(-g*I_c,Indices(0,D-1,1));
  C1.setElement(g*I_c,Indices(0,D-1,2));
  C1.setElement(-.5*gamma*ONE_c,Indices(0,D-1,3));
  C1.setElement(-.5*gamma*ONE_c,Indices(0,D-1,4));
  C1.setElement(gamma*ONE_c,Indices(0,D-1,5));
  // And also for the last site
  CN.setElement(-g*I_c,Indices(0,0,1));
  CN.setElement(g*I_c,Indices(0,0,2));
  CN.setElement(-.5*gamma*ONE_c,Indices(0,0,3));
  CN.setElement(-.5*gamma*ONE_c,Indices(0,0,4));
  CN.setElement(gamma*ONE_c,Indices(0,0,5));
  // And now the two-body terms
  // gamma ( sigma_minus x Id ) (x) ( Id x sigma_minus ) 
  C.setElement(ONE_c,Indices(0,1,6));
  C.setElement(gamma*ONE_c,Indices(1,D-1,8));
  C1.setElement(ONE_c,Indices(0,1,6));
  CN.setElement(gamma*ONE_c,Indices(1,0,8));
  // gamma ( Id x sigma_minus ) (x) ( sigma_minus x Id ) 
  C.setElement(ONE_c,Indices(0,3,8));
  C.setElement(gamma*ONE_c,Indices(3,D-1,6));
  C1.setElement(ONE_c,Indices(0,3,8));
  CN.setElement(gamma*ONE_c,Indices(3,0,6));
  // -0.5*gamma ( sigma_plus x Id ) (x) ( sigma_minus x Id ) 
  C.setElement(ONE_c,Indices(0,2,7));
  C.setElement(-.5*gamma*ONE_c,Indices(2,D-1,6));
  C1.setElement(ONE_c,Indices(0,2,7));
  CN.setElement(-.5*gamma*ONE_c,Indices(2,0,6));
  // -.5*gamma ( sigma_minus x Id ) (x) ( sigma_plus x Id ) 
  C.setElement(-.5*gamma*ONE_c,Indices(1,D-1,7));
  CN.setElement(-.5*gamma*ONE_c,Indices(1,0,7));
  // -0.5*gamma ( Id x sigma_plus ) (x) ( Id x sigma_minus ) 
  C.setElement(ONE_c,Indices(0,4,9));
  C.setElement(-.5*gamma*ONE_c,Indices(4,D-1,8));
  C1.setElement(ONE_c,Indices(0,4,9));
  CN.setElement(-.5*gamma*ONE_c,Indices(4,0,8));
  // -0.5*gamma ( Id x sigma_minus ) (x) ( Id x sigma_plus ) 
  C.setElement(-.5*gamma*ONE_c,Indices(3,D-1,9));
  CN.setElement(-.5*gamma*ONE_c,Indices(3,0,9));

  // Now I reshape, contract with Z and set the indices in proper order
  C.reshape(Indices(D*D,nrOps));
  C1.reshape(Indices(1*D,nrOps));
  CN.reshape(Indices(D*1,nrOps));
  Z.reshape(Indices(nrOps,d*d*d*d));
  C.multiplyRight(Z);C.reshape(Indices(D,D,d*d,d*d));C.permute(Indices(3,1,4,2));
  C1.multiplyRight(Z);C1.reshape(Indices(1,D,d*d,d*d));C1.permute(Indices(3,1,4,2));
  CN.multiplyRight(Z);CN.reshape(Indices(D,1,d*d,d*d));CN.permute(Indices(3,1,4,2));

  // Construct the three operators
  Operator* Op1=new Operator(C1);
  Operator* Op=new Operator(C);
  Operator* OpN=new Operator(CN);

  // And set them to the adequate positions in the MPO!
  Lmpo.setOp(0,Op1,true); // Lmpo is now resposible of this pointer
  Lmpo.setOp(N-1,OpN,true);
  if(N>2){
    Lmpo.setOp(1,Op,true);
    for(int k=2;k<N-1;k++){
      // For all the remaining middle sites, I do not need to create
      // new Operators, but just point to the same one!
      Lmpo.setOp(k,Op,false);      
    }
  }
  else{
    // delete the extra Op, to avoid memory leaks, in this very silly case of only two sites
    delete Op;
  }

  // to avoid repeating all this, set a flag that indicates that L is already prepared
  isInit=true;
}

// void CoherentDissipationLiouvillianLR::initTrL(){
//   // I can only do this if the Lmpo is already prepared
//   initL();

//   // I will create copies of the edge and middle operators with the
//   // adequate contraction of indices (upi,upj) and complex conjugation
//   // and for the rest of the chain, I just point to the operator in
//   // the second site, as L itself does

//   mwArray auxI=identityMatrix(d);
//   auxI.reshape(Indices(1,d*d));

//   // First of the chain
//   mwArray op=Lmpo.getOp(0).getFullData();
//   Indices dimOp=op.getDimensions();
//   op.reshape(Indices(dimOp[0],-1));
//   op.multiplyLeft(auxI);
//   op.reshape(Indices(dimOp[1],d,d,dimOp[3]));
//   op.permute(Indices(3,1,2,4));
//   trLmpo.setOp(0,new Operator(op),true);

//   // Last of the chain
//   op=Lmpo.getOp(N-1).getFullData();
//   dimOp=op.getDimensions();
//   op.reshape(Indices(dimOp[0],-1));
//   op.multiplyLeft(auxI);
//   op.reshape(Indices(dimOp[1],d,d,dimOp[3]));
//   op.permute(Indices(3,1,2,4));
//   trLmpo.setOp(N-1,new Operator(op),true);

//   if(N>2){
//     op=Lmpo.getOp(1).getFullData();
//     dimOp=op.getDimensions();
//     op.reshape(Indices(dimOp[0],-1));
//     op.multiplyLeft(auxI);
//     op.reshape(Indices(dimOp[1],d,d,dimOp[3]));
//     op.permute(Indices(3,1,2,4));
//     trLmpo.setOp(1,new Operator(op),true);
//     for(int k=2;k<N-1;k++)
//       trLmpo.setOp(k,&trLmpo.getOp(1),false);            
//   }
// }

// void CoherentDissipationLiouvillianLR::getTrLMPO(int dpur,MPO& result){
//   result.initLength(N);
//   mwArray auxId=identityMatrix(dpur);
//   auxId.reshape(Indices(dpur,1,dpur,1));

//   // Set operator 0
//   mwArray op=trLmpo.getOp(0).getFullData();
//   result.setOp(0,new DoubleOperator(op,auxId),true);
//   // Set operator N-1
//   op=trLmpo.getOp(N-1).getFullData();
//   result.setOp(N-1,new DoubleOperator(op,auxId),true);
//   if(N>2){
//     op=trLmpo.getOp(1).getFullData();
//     result.setOp(1,new DoubleOperator(op,auxId),true);
//     for(int k=2;k<N-1;k++)
//       result.setOp(k,&result.getOp(1),false);            
//   }
// }
