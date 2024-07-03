/**
   \file Contractor.cpp
   Implementation of the basic functionality with MPS-MPO operations,
   and the main (finite) MPS algorithms.
   
   \author Mari-Carmen Banuls
   \date 16/03/2011
*/


#include "Contractor.h"
#include "BasicMultiplier.h"
// // #include "TensorMultiplier.h"
#include "TensorMultiplierProj.h"
#include "TensorMultiplierOffset.h"
#include "TensorMultiplierHermitianOffset.h"
#include "TensorMultiplierProjMulti.h"
//#include "HermitianTensorMultiplier.h"
#include "JoinedOperator.h"
#include "FoldedOperator.h"

#include <iomanip>

#define MAXITER 20000 // max nr of trials for eigenvector
#define MAXPERT 1 // max nr of perturbations and retrials for eigenvector
#define MAX_ROUND_EXCTD 1000 // max nr of rounds for excited states

#define NUMPREC 2.2E-16 // numerical precision

#define LOWEST_PRIMME 801 // min dim to call Primme eigensolver

//#define SOLVEN 1  // version in which A=M/N at each step
//#undef SOLVEN

using namespace std;
using namespace shrt;

// The private constructor only sets the default tolerances
Contractor::Contractor():
  svdtol(1E-8),convtol(1E-5),eigtol(0.),excorthtol(0.),solver(arpack){
  cout<<"Constructed default Contractor"<<endl;
}


Contractor::~Contractor(){
  // Destroy temporary storage, if any?
}

Contractor& Contractor::theContractor(){
  static Contractor theContr;
  // TODO: Do I need more things in here?
  return theContr;
}

const string Contractor::getSolverName() const{
  switch (solver){
  case arpack:
    return "arpack";
    break;
  case fulleig:
    return "fulleig";
    break;
  case primme:
    return "primme";
    break;
  case primme_JDQR:
    return "primme_JDQR";
    break;
  case primme_arnoldi:
    return "primme_arnoldi";
    break;
  default:
    return "Unknown solver "+solver;
  }
}

void Contractor::optimize(const MPO& ops,const MPS& orig,
			  MPS& init,int D,double* errN,double* errT){
  /** I approximate by an MPS satisfying gauge condition, times a 
      normalization factor. This makes the behaviour more stable
      and allows a simpler optimization. Moreover, I avoid the calculation
      of the (expensive) N matrix, as this is guaranteed to be the 
      identity. 
  */
  // Check some dimensions
  if(ops.getLength()!=orig.getLength()){
    cout<<"Error: incompatible sizes in Contractor::optimize(MPO-"
	<<ops.getLength()<<",MPS-"<<orig.getLength()<<")"<<endl;
    exit(2);
  }
  // cout<<"Contractor::optimize starts with orig ket "<<orig;
  // cout<<"\nStarting with init ket "<<init;
  // cout<<"\nUsing Dmax="<<D<<endl;
  if(init.isEmpty()) init=orig;
  vector<int> findims=ops.getDimensions();
  init.adjustPhysDimensions(findims);
  //cout<<"*Adjusted physical dimensions of init "<<findims<<", now init "<<init<<endl;
  if(D>0)
    if(D>init.getBond())
      init.increaseBondDimension(D);
    else if(D<init.getBond()){
      cout<<"WARNING: Contractor::optimize called with D("<<D<<") smaller "
	  <<"than in the starting MPS (D_init="<<init.getBond()<<"). Will truncate "
	  <<"the initial MPS!"<<endl;
      init=MPS(init,D);
    }
  //  if(D>0){
  //init.increaseBondDimension(D);
  //}
  init.gaugeCond('R',1); // this first gauge ensures right Dr dimensions
  //cout<<"*Applied ritht gauge condition"<<endl;
  init.gaugeCond('L',1); /* this ensures proper Dl and allows starting 
			    from left.*/
  //cout<<"*Applied left gauge condition"<<endl;
  //  else cout<<"Keeping my old bond dimension of init MPS: "<<init.getBond()<<endl;
 
  // cout<<"Contractor::optimize from orig of norm "<<contract(orig,orig);
  // cout<<", norm fact="<<orig.getNormFact()<<", and starting with init of norm="
  //   <<contract(init,init)<<", norm fact="<<init.getNormFact()<<endl;
  // cout<<"Overlap INICIAL <init|ops|orig>="<<contract(orig,ops,init)
  //   <<", <orig|init>="<<contract(init,orig)<<endl;

  int nrsites=orig.getLength();
  //cout<<"Contractor:optimize, using nrsites="<<nrsites<<endl;

  // Initialize temporary storage with gauge condition=true
  bool gauge=true;
  bool ketN=0; // For this optimization, the norm terms should be done with bra (init)
  TmpManager tmpMgr(nrsites,gauge);

  // Start optimization: sweep to right and left until convergence
  bool done=0;
  bool right=1;
  int pos=0;
  double distance=abs(contract(orig,orig))+abs(contract(init,init)); // some absurdly large initial value (it would be better to take the true contraction)
  int round=1;
  mwArray M; //mwArray N;
  mwArray contr1,contr2,newA;
  //  double start=clock();double finish(0.);
  while(!done){
    /* Optimize matrix for position pos */
    //int d=init.getA(pos).getd();
    // Calculate M(pos) matrix (N not needed)
    // 1. Obtain left and right terms for everything
    //if(pos==0){
    // finish=clock();
    // cout<<"Contractor::optimize for pos "<<pos<<", round "<<round
    // 	<<", distance="<<distance<<" time (previous round)"<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
    // start=clock();
    //}
    calculateL(pos,tmpMgr,ops,orig,init,gauge,ketN);
    calculateR(pos,tmpMgr,ops,orig,init,gauge,ketN);
    // 2. Contract the terms appropriately for M matrix
    //cout<<"After calculateL->"<<tmpMgr.operL[pos]<<", calculateR->"
    //	<<tmpMgr.operR[pos]<<", with op->"<<ops.getOp(pos)<<endl;
    calculateMket(M,tmpMgr.operL[pos],orig.getA(pos),ops.getOp(pos),
		  tmpMgr.operR[pos]);
    if(pos==0&&round==1){
      if(M==mwArray(ZERO_c)){ // \TODO check if this test is enough!!
	init.perturb(1.);
	for(int l=0;l<nrsites;l++){
	  tmpMgr.outofdateL(l);
	  tmpMgr.outofdateR(l);
	}
	continue;
      }
    }
    // 3. Solve the linear equations
    solveSiteGauge(newA,contr1,contr2,M,init.getNormFact()); // Solo el vector! 
    //  cout<<"For pos "<<pos<<", where A="<<init.getA(pos).getA()<<
    // 	"the solution is "<<newA<<endl;
    init.setA(pos,newA); 

    //cout<<"After solveSiteGauge contr1="<<contr1<<", contr2="<<contr2<<endl;
    //  cout<<"For pos "<<pos<<", M matrix calculated of size "<<M.GetDimensions()
    // 	   <<" yields new site "<<init.getA(pos).getA().GetDimensions()<<endl;
    //        cout<<"Namely A="<<init.getA(pos).getA()<<endl;
    if(right)
      init.gaugeCond(pos,'R');
    else
      init.gaugeCond(pos,'L');
    // This has changed contr1 AND contr2!!!
    mwArray A=init.getA(pos).getA();
    A.reshape(Indices(A.getNrComponents(),1));// vector form
    contr1=Hconjugate(A)*A;
    contr2=Hconjugate(A)*M; // Calculates <A|M>     
    //cout<<"After changing the tensor, and applying gauge, contr1="
    //<<contr1<<", contr2="<<contr2<<endl;

    // when the edge is reached, find the norm
    //       if((right&&pos==nrsites-1)||((!right)&&pos==0)){
    // 	//cout<<"At pos "<<pos<<", round "<<round<<", going right="<<right
    // 	// <<", when norm was "<<setprecision(12)<<(double)init.getNormFact();
    // 	init.setNormFact(((double)contr2.Real())/(double)contr1);
    // 	//cout<<" setting norm to "<<init.getNormFact()<<endl;
    // 	if(init.getNormFact()<1E-12){
    // 	  cout<<"WARNING!!! Norm="<<init.getNormFact()
    // 	      <<"<1E-12, forcing end of optimize!!"<<endl;
    // 	  done=1;break;
    // 	}
    //       }
    
    // notify temporary storage as appropriate
    for(int l=pos+1;l<nrsites;l++)
      tmpMgr.outofdateL(l);
    for(int l=pos-1;l>=0;l--)
      tmpMgr.outofdateR(l);
    
    // Implementar la lógica de avanzar uno, a izqda o drcha
    if(right){ //moving right
      if(pos==nrsites-1){ //last reached: check convergence and change sense
	// 	cout<<"Last site, round "<<round<<"=> now move left, and check convergence, "
	// 	    <<"when contr2="<<contr2<<endl;
	// 	//pos--; // Do this one again for proper gauge! (I could also apply gauge here!)
	right=0;
	done=checkDistanceConvergence(distance,contr1,contr2,
				      init.getNormFact(),
				      orig.getNormFact(),round);
	if(!done){
	  init.gaugeCond(pos,'L');if(pos>0)pos--;
	}
	round++;
      }
      else{
	pos++;
      }
    }
    else{// moving left
      if(pos==0){ //first reached: check convergence and change sense
	//cout<<"optimize at pos 0, round "<<round<<", distance "<<distance<<endl;
	// 	cout<<"First site, round "<<round<<"=> now move right, and check convergence, "
	// 	    <<"when contr2="<<contr2<<endl;
	// 	//pos++;
	right=1;
	done=checkDistanceConvergence(distance,contr1,contr2,
				      init.getNormFact(),orig.getNormFact(),
				      round);
	if(!done){
	  init.gaugeCond(pos,'R');if(pos<nrsites-1)pos++;
	}
	round++;
      }
      else{
	pos--;
      }
    }    
  }
  
  // Normalize only here!
  //init.setNormFact(((double)contr2.Real())/(double)contr1);
  double reMA=real(contr2.getElement(0)); // first element
  double ANA=real(contr1.getElement(0)); // first element
  init.setNormFact(reMA/ANA);
  
  // When the state has converged, I have something which is not normalized.
  // I could normalize it as
  // gaugeCond
  // But this depends on what I am doing=> do outside!
  // In any case, if the original vector had a normFactor, it was ignored 
  // in the above calculation, hence it has to be added by hand now!
  //cout<<"At the end of optimize, changing norm "<<init.getNormFact();
  init.setNormFact(init.getNormFact()*orig.getNormFact());
  //  cout<<" -> "<<init.getNormFact()<<endl;
  // If error has to be computed...
  if(errN!=0){
    double norm0=real(contract2(ops,orig));
    double norm1=real(contract(init,init));
    complex_t crossR=contract(orig,ops,init);
    *errN=norm0+norm1-2*real(crossR);
    if(*errN<-1E-12){
      cout<<"WARNING: Contractor finds err^2<0!!!"<<endl;
    }
    if(errT!=0){
      *errT=norm0;
    }
  }
  tmpMgr.clear();
}

void Contractor::optimizeL(const MPO& ops,const MPS& orig,
			   MPS& init,int D,double* errN,double* errT){

  // Check some dimensions
  if(ops.getLength()!=orig.getLength()){
    cout<<"Error: incompatible sizes in Contractor::optimizeL(MPO-"
	<<ops.getLength()<<",MPS-"<<orig.getLength()<<")"<<endl;
    exit(2);
  }
  if(init.isEmpty()) init=orig;
  // I approximate by an MPS satisfying gauge condition, times a 
  // normalization factor. This makes the behaviour more stable
  // and allows a simpler optimization. Moreover, I avoid the calculation
  // of the (expensive) N matrix, as this is guaranteed to be the 
  // identity.

  // WARNING: The dimensions of the target MPS must agree with the ORIGINAL
  // physical dimensions of ops. 
  vector<int> findims=ops.getOriginalDimensions();
  init.adjustPhysDimensions(findims);
  if(D>0)
    if(D>init.getBond())
      init.increaseBondDimension(D);
    else if(D<init.getBond()){
      cout<<"WARNING: Contractor::optimizeL called with D("<<D<<") smaller "
	  <<"than in the starting MPS (D_init="<<init.getBond()<<"). Will truncate "
	  <<"the initial MPS!"<<endl;
      init=MPS(init,D);
    }
  init.gaugeCond('R',1); // this first gauge ensures right Dr dimensions
  init.gaugeCond('L',1); /* this ensures proper Dl and allows starting 
			    from left.*/
  //  if(D>0){
  //init.increaseBondDimension(D);
    //cout<<"Increased Bond Dimension to "<<D<<endl;
  //}
  //else cout<<"Keeping my old bond dimension of init MPS: "<<init.getBond()<<endl;
  //   cout<<"Contractor::optimizeL from orig of norm "<<contract(orig,orig);
  //   cout<<", norm fact="<<orig.getNormFact()<<", and starting with init of norm="
  //     <<contract(init,init)<<", norm fact="<<init.getNormFact()<<endl;
  int nrsites=orig.getLength();
  //cout<<"Contractor:optimizeL, using nrsites="<<nrsites<<endl;
  // Initialize temporary storage with gauge condition=true
  bool gauge=true;
  bool ketN=true; // Only for Norm terms=>with gauge this has no effect
  TmpManager tmpMgr(nrsites,gauge);

  // Start optimization: sweep to right and left until convergence
  bool done=0;
  bool right=1;
  int pos=0;
  // bool right=0;
  // int pos=nrsites-1;
  double distance=1E10;
  int round=1;

  mwArray M; 
  mwArray contr1,contr2,newA;
  mwArray auxA;
  while(!done){
    /* Optimize matrix for position pos */
    //int d=init.getA(pos).getd();      
    // Calculate M(pos) matrix (N not needed)
    // 1. Obtain left and right terms for everything
    //cout<<"Contractor::optimizeL for pos "<<pos<<", round "<<round<<endl;
    calculateL(pos,tmpMgr,ops,init,orig,gauge,ketN);
    //cout<<"Contractor::optimizeL for pos "<<pos<<", after calculateL "
    //	<<"tmpL term: "<<tmpMgr.operL[pos]<<endl;
    calculateR(pos,tmpMgr,ops,init,orig,gauge,ketN);
    //cout<<"Contractor::optimizeL for pos "<<pos<<", after calculateR "
    //	<<"tmpR term: "<<tmpMgr.operR[pos]<<endl;
    // 2. Contract the terms appropriately for M matrix
    //calculateMket(M,tmpMgr.operL[pos],orig.getA(pos),ops.getOp(pos),
    //		  tmpMgr.operR[pos]);
    calculateMbra(M,tmpMgr.operL[pos],orig.getA(pos),ops.getOp(pos),
		  tmpMgr.operR[pos]);
    //cout<<"termL="<<tmpMgr.operL[pos]<<", termR="<<tmpMgr.operR[pos]
    //	  <<", site="<<orig.getA(pos)<<endl;

    if(pos==0&&round==1){
      if(M==mwArray(ZERO_c)){
	init.perturb(1.);
	for(int l=0;l<nrsites;l++){
	  tmpMgr.outofdateL(l);
	  tmpMgr.outofdateR(l);
	}
	continue;
      }
    }
    // 3. Solve the linear equations
    solveSiteGauge(newA,contr1,contr2,M,init.getNormFact()); // Solo el vector! (normalizado)  
    // 4. Set the site in the MPS, but take into account that I have solved the equation for 
    // the conjugate
    //cout<<"After solveSite@"<<pos<<", contr1="<<contr1<<", contr2="<<contr2<<endl;
    //cout<<"newA="<<newA<<endl;
    init.setA(pos,conjugate(newA)); 
    if(right) init.gaugeCond(pos,'R');
    else init.gaugeCond(pos,'L');
    // This has changed contr1 AND contr2!!!
 
      // TODO: Update more efficiently?
    mwArray A=init.getA(pos).getA();
    A.reshape(Indices(A.getNrComponents(),1));// vector form
    contr1=Hconjugate(A)*A;
    contr2=permute(M,Indices(2,1))*A;
    //cout<<"After updating@"<<pos<<", contr1="<<contr1<<", contr2="<<contr2<<endl;
    //       // When the edge is reached, find the norm
    //       double oldN=init.getNormFact();
    //       if((right&&pos==nrsites-1)||(!right&&pos==0)){
    // 	init.setNormFact(((double)contr2.Real())/(double)contr1);
    // 	if(init.getNormFact()<1E-12){
    // 	  cout<<"WARNING!!! N<1E-12, forcing end of optimize!!"<<endl;
    // 	  done=1;break;
    // 	}
    // 	// init.setNormFact((double)contr2.Real());
    //        // cout<<"Setting norm factor of state to "<<init.getNormFact()<<endl;
    //       }
    
    // notify temporary storage as appropriate
    for(int l=pos+1;l<nrsites;l++)
      tmpMgr.outofdateL(l);
    for(int l=pos-1;l>=0;l--)
      tmpMgr.outofdateR(l);
    
    // Implementar la lógica de avanzar uno, a izqda o drcha
    if(right){ //moving right
      if(pos==nrsites-1){ //last reached: check convergence and change sense
	// 	cout<<"Last site, round "<<round<<"=> now move left, and check convergence, "
	// 	    <<"when contr2="<<contr2<<endl;
	right=0;
	done=checkDistanceConvergence(distance,contr1,contr2,
				      init.getNormFact(),
				      orig.getNormFact(),round);
	//cout<<"Current distance "<<distance<<endl;
	round++;
      }
      else{
	pos++;
      }
    }
    else{// moving left
      if(pos==0){ //first reached: check convergence and change sense
	//cout<<"optimizeL at pos 0, round "<<round<<", distance "
	//  <<distance<<endl;	
	// 	cout<<"First site, round "<<round<<"=> now move right, and check convergence, "
	// 	    <<"when contr2="<<contr2<<endl;
	right=1;
	done=checkDistanceConvergence(distance,contr1,contr2,
				      init.getNormFact(),orig.getNormFact(),
				      round);
	//cout<<"Current distance "<<distance<<endl;
	round++;
      }
      else{
	pos--;
      }
    }    
  }
    // TEST! Set NOW the norm factor
    // When the edge is reached, find the norm
  double reMA=real(contr2.getElement(0)); // first element
  double ANA=real(contr1.getElement(0)); // first element
  init.setNormFact(reMA/ANA);
  // When the state has converged, I have something which is not normalized.
  // I could normalize it as
  // gaugeCond
  // But this depends on what I am doing=> do outside!
  // In any case, if the original vector had a normFactor, it was ignored 
  // in the above calculation, hence it has to be added by hand now!
  //cout<<"At the end of optimizeL, changing norm "<<init.getNormFact();
  init.setNormFact(init.getNormFact()*orig.getNormFact());
  //cout<<" -> "<<init.getNormFact()<<endl;
  // If error has to be computed...
  if(errN!=0){
    double norm0=real(contract2(orig,ops));
    double norm1=real(contract(init,init));
    complex_t crossL=contract(init,ops,orig);
    *errN=norm0+norm1-2*real(crossL);
    if(*errN<-1E-12){
      cout<<"WARNING: Contractor finds err^2<0!!!"<<endl;
    }
    if(errT!=0){
      *errT=norm0;
    }
  }
}

complex_t Contractor::contract(const MPS& ket,const MPS& bra,char dir){
  // TODO!!! Use the lower method
  //  cout<<"Contractor:contract("<<bra<<","<<ket<<")"<<endl;
  int nrsites=ket.getLength();
  if(nrsites!=bra.getLength()){
    cout<<"Error in contract: dimensions do not coincide: bra ";
    cout<<bra.getLength()<<", ket "<<nrsites<<endl;
    exit(212);
  }
  mwArray result(ONE_c); // here the result is to be saved
  //  mwArray& tmp=init;
  switch(dir){
  case 'R':
    for(int k=nrsites-1;k>=0;k--){
      mwArray tmp(result);
      contractR(result,tmp,ket.getA(k),bra.getA(k));
    }
    break;
  case 'L':
  default:
    for(int k=0;k<nrsites;k++){
      mwArray tmp(result);
      contractL(result,tmp,ket.getA(k),bra.getA(k));
    }
    break;
  }
  // At the end: check size (1x1) and return value!
  if(!result.isScalar()){
    cout<<"Error: Contractor::contract(MPS,MPS) produces result "<<result<<endl;
    exit(212);
  }
  return (ket.getNormFact()*bra.getNormFact())*result.getElement(0);
}

mwArray Contractor::contract(const MPS& ket,const MPS& bra,int pos,char dir){
  int nrsites=ket.getLength();
  if(nrsites!=bra.getLength()){
    cout<<"Error in contract: dimensions do not coincide: bra ";
    cout<<bra.getLength()<<", ket "<<nrsites<<endl;
    exit(212);
  }
  mwArray result(ONE_c); // here the result is to be saved
  //  mwArray& tmp=init;
  switch(dir){
  case 'R':
    for(int k=nrsites-1;k>pos;k--){
      mwArray tmp(result);
      contractR(result,tmp,ket.getA(k),bra.getA(k));
    }
    break;
  case 'L':
  default:
    for(int k=0;k<pos;k++){
      mwArray tmp(result);
      contractL(result,tmp,ket.getA(k),bra.getA(k));
    }
    break;
  }
  return result;
}

complex_t Contractor::contract(const MPS& ket,const MPO& ops,
			     const MPS& bra){
  //cout<<"Contractor:contract("<<bra<<","<<ops<<","<<ket<<")"<<endl;
  int nrsites=ket.getLength();
  if((nrsites!=bra.getLength())||(nrsites!=ops.getLength())){
    cout<<"Error in contract: dimensions do not coincide: bra ";
    cout<<bra.getLength()<<", ket "<<nrsites<<", ops "<<
      ops.getLength()<<endl;
    exit(212);
  }
  // TODO: Check lengths and dimensions!
  mwArray result(ONE_c); // here the result is to be saved
  //mwArray& tmp=init;
  for(int k=0;k<nrsites;k++){
    //cout<<"Before step "<<k<<", result="<<result.getDimensions()<<endl;
    //cout<<"Before step "<<k<<", result="<<result<<endl;
    mwArray tmp(result);
    contractOperL(result,tmp,ket.getA(k),ops.getOp(k),bra.getA(k));
    //tmp=result; // CHECK!! Funciona esto? Le paso igual input y output!!
  }
  // At the end: check size (1x1) and return value!
  // if the given MPSs have normalization factors, I must multiply them now
  if(!result.isScalar()){
    cout<<"Error: after Contractor::contract(MPS,MPO,MPS), result is "
	<<result<<endl;
    exit(212);
  }
  return (ket.getNormFact()*bra.getNormFact())*result.getElement(0);
}

mwArray Contractor::contract(const MPS& ket,const MPO& ops,
			       const MPS& bra,int pos,char dir){
  //cout<<"Contractor:contract("<<bra<<","<<ops<<","<<ket<<")"<<endl;
  int nrsites=ket.getLength();
  if((nrsites!=bra.getLength())||(nrsites!=ops.getLength())){
    cout<<"Error in contract: dimensions do not coincide: bra ";
    cout<<bra.getLength()<<", ket "<<nrsites<<", ops "<<
      ops.getLength()<<endl;
    exit(212);
  }
  mwArray result(ONE_c); // here the result is to be saved
  switch(dir){
  case 'R':
    for(int k=nrsites-1;k>pos;k--){
      //cout<<"Before step "<<k<<", result="<<result<<endl;
      mwArray tmp(result);
      contractOperR(result,tmp,ket.getA(k),ops.getOp(k),bra.getA(k));
    }
    break;
  case 'L':
  default:
    for(int k=0;k<pos;k++){
      //cout<<"Before step "<<k<<", result="<<result.getDimensions()<<endl;
      //cout<<"Before step "<<k<<", result="<<result<<endl;
      mwArray tmp(result);
      contractOperL(result,tmp,ket.getA(k),ops.getOp(k),bra.getA(k));
      //tmp=result; // CHECK!! Funciona esto? Le paso igual input y output!!
    }
    break;
  }
    return result;
}

mwArray Contractor::contractPart(const MPS& ket,const MPO& ops,const MPS& bra,int posL,int posR,char dir){
  //cout<<"contractPart("<<posL<<","<<posR<<")"<<endl;
  int nrsites=ket.getLength();
  if((nrsites!=bra.getLength())||(nrsites!=ops.getLength())){
    cout<<"Error in contract: dimensions do not coincide: bra ";
    cout<<bra.getLength()<<", ket "<<nrsites<<", ops "<<
      ops.getLength()<<endl;
    exit(212);
  }
  mwArray result; // here the result is to be saved
  // for simplicity, I construct a first term with the identity of the
  // right dimensions
  int dimL=ket.getA(posL).getDl()*bra.getA(posL).getDl()*ops.getOp(posL).getDl();
  int dimR=ket.getA(posR).getDr()*bra.getA(posR).getDr()*ops.getOp(posR).getDr();
  switch(dir){
  case 'R':
    cout<<"Not yet supported, contracted from left"<<endl;
    // result=identityMatrix(dimR);
    // result.reshape(Indices(ket.getA(posR).getDr(),bra.getA(posR).getDr(),ops.getOp(posR).getDr(),dmR));
    // for(int k=posR;k>=posL;k--){
    //   //cout<<"Before step "<<k<<", result="<<result<<endl;
    //   mwArray tmp(result);
    //   contractOperRmiddle(result,tmp,ket.getA(k),ops.getOp(k),bra.getA(k));
    // }
    // break;
  case 'L':
  default:
    result=identityMatrix(dimL);
    result.reshape(Indices(dimL,ket.getA(posL).getDl(),ops.getOp(posL).getDl(),bra.getA(posL).getDl()));
    for(int k=posL;k<=posR;k++){
      //cout<<"Before step "<<k<<", result="<<result.getDimensions()<<endl;
      //cout<<"Before step "<<k<<", result="<<result<<endl;
      mwArray tmp(result);
      contractOperLmiddle(result,tmp,ket.getA(k),ops.getOp(k),bra.getA(k));
    }
    result.reshape(Indices(dimL,-1));
    break;
  }
  return result;

}


complex_t Contractor::contract2(const MPO& ops,const MPS& ket){
  //cout<<"Contractor:contract("<<bra<<","<<ops<<","<<ket<<")"<<endl;
  int nrsites=ket.getLength();
  if(nrsites!=ops.getLength()){
    cout<<"Error in contract: dimensions do not coincide: ket ";
    cout<<nrsites<<", ops "<<ops.getLength()<<endl;
    exit(212);
  }
  mwArray result(ONE_c); // here the result is to be saved

  for(int k=0;k<nrsites;k++){
    //cout<<"contract2:: Before step "<<k<<", result="<<result<<endl;
    mwArray tmp(result);
    ops.getOp(k).contractL2(result,tmp,ket.getA(k),ket.getA(k));
  }
  // At the end: check size (1x1) and return value
  // if the given MPSs have normalization factors, I must multiply them now
  if(!result.isScalar()){
    cout<<"Error: after Contractor::contract2(MPO,MPS), result is "<<result<<endl;
    exit(212);
  }
  return (ket.getNormFact()*ket.getNormFact())*result.getElement(0);
  //  cout<<"And now a TEST!!!"<<endl;
  //   TmpMgr tmpMgr2(nrsites,0);int pos=nrsites-1;
  //   calculateL2(pos,tmpMgr2,ops,ket,ket);
  //   calculateR2(pos,tmpMgr2,ops,ket,ket);
  //   mwArray N;
  //   ops.getOp(pos).contractN2(N,tmpMgr2.normL[pos],tmpMgr2.normR[pos]);
  //   m_times(1,Nn,N,mwArray(ket.getNormFact()*ket.getNormFact()));
  // 	m_times(1,Mn,N,mwArray(init.getNormFact()));
  //        solveSite(newA,contr1,contr2,M,N); // Solve Matrix Equation (N not Id)
  //       double newdist=constantN+(double)contr1-2*(double)contr2.Real();
  //       cout<<"Distance by tmpmgr2 "<<newdist<<", simple contraction is "<<(double)contr2.Real()*orig.getNormFact();
  //       if(contr2.IsComplex())
  // 	cout<<"+ "<<(double)contr2.Imag()*orig.getNormFact()<<" i ";
  //       cout<<", and double contraction is "<<(double)contr1<<endl;
  //       cout<<"newA this time is "<<newA<<endl;  /** EEEEHHHH*/
}

complex_t Contractor::contract2(const MPS& bra,const MPO& ops,const MPS& ket){
  //cout<<"Contractor:contract("<<bra<<","<<ops<<","<<ket<<")"<<endl;
  int nrsites=ket.getLength();
  if(nrsites!=ops.getLength()){
    cout<<"Error in contract2(bra,Op,ket): dimensions do not coincide: ket ";
    cout<<nrsites<<", ops "<<ops.getLength()<<endl;
    exit(212);
  }
  mwArray result(ONE_c); // here the result is to be saved
  for(int k=0;k<nrsites;k++){
    //cout<<"Before step "<<k<<", result="<<result.GetDimensions()<<endl;
    mwArray tmp(result);
    ops.getOp(k).contractL2(result,tmp,ket.getA(k),bra.getA(k));
  }
  // At the end: check size (1x1) and return value!
  // if the given MPSs have normalization factors, I must multiply them now
  if(!result.isScalar()){
    cout<<"Error: after Contractor::contract2(MPS,MPO,MPS), result is "<<result<<endl;
    exit(212);
  }
  return (bra.getNormFact()*ket.getNormFact())*result.getElement(0);
}


complex_t Contractor::contract2(const MPS& bra,const MPO& ops){
  //cout<<"Contractor:contract("<<bra<<","<<ops<<","<<ket<<")"<<endl;
  int nrsites=bra.getLength();
  if(nrsites!=ops.getLength()){
    cout<<"Error in contract2: dimensions do not coincide: bra ";
    cout<<nrsites<<", ops "<<ops.getLength()<<endl;
    exit(212);
  }
  mwArray result(ONE_c); // here the result is to be saved
  for(int k=0;k<nrsites;k++){
    //cout<<"Before step "<<k<<", result="<<result.GetDimensions()<<endl;
    mwArray tmp(result);
    ops.getOp(k).contractL2(result,tmp,bra.getA(k),bra.getA(k),true);
  }
  // At the end: check size (1x1) and return value!
  // if the given MPSs have normalization factors, I must multiply them now
  if(!result.isScalar()){
    cout<<"Error: after Contractor::contract2(MPS,MPO), result is "<<result<<endl;
    exit(212);
  }
  return (bra.getNormFact()*bra.getNormFact())*result.getElement(0);
}




void Contractor::calculateL(int pos,TmpManager& theMgr,const MPO& ops,
			    const MPS& ket,const MPS& bra,
			    bool gauge,bool ketN){
  //cout<<"Contractor::calculateL(pos="<<pos<<", tmpMgr="<<&theMgr<<")"<<endl;
  if(!theMgr.isuptodateL(pos)){ //else, nothing to do, just return!
    calculateL(pos-1,theMgr,ops,ket,bra,gauge,ketN); // take the one before
    //cout<<"Using tmpL of pos "<<pos-1<<" uptodate"<<endl;
    contractOperL(theMgr.operL[pos],theMgr.operL[pos-1],ket.getA(pos-1),
 		  ops.getOp(pos-1),bra.getA(pos-1));
     //contractOperL(theMgr.getOperL(pos),theMgr.getOperL(pos-1),ket.getA(pos-1),
     //	  ops.getOp(pos-1),bra.getA(pos-1));
    if(!gauge){
      if(!ketN)
	contractNormL(theMgr.normL[pos],theMgr.normL[pos-1],bra.getA(pos-1));
      else
	contractNormL(theMgr.normL[pos],theMgr.normL[pos-1],ket.getA(pos-1));
    }
    theMgr.storedL(pos);    
  }
}

void Contractor::calculateR(int pos,TmpManager& theMgr,const MPO& ops,
			    const MPS& ket,const MPS& bra,
			    bool gauge,bool ketN){
  //cout<<"Contractor::calculateR(pos="<<pos<<", tmpMgr="<<&theMgr<<")"<<endl;
  if(!theMgr.isuptodateR(pos)){ //else, nothing to do, just return!
    calculateR(pos+1,theMgr,ops,ket,bra,gauge,ketN); // take next one and use it!
    // cout<<"Using tmpR of pos "<<pos+1<<" uptodate, with ket "
    //   //	<<endl;
    // 	<<ket.getA(pos+1).getA().getDimensions()<<", op "<<ops.getOp(pos+1).getFullData().getDimensions()<<" and bra "<<bra.getA(pos+1).getA().getDimensions()<<endl;
    contractOperR(theMgr.operR[pos],theMgr.operR[pos+1],ket.getA(pos+1),
		  ops.getOp(pos+1),bra.getA(pos+1));
    if(!gauge){
      if(!ketN)
	contractNormR(theMgr.normR[pos],theMgr.normR[pos+1],bra.getA(pos+1));
      else
	contractNormR(theMgr.normR[pos],theMgr.normR[pos+1],ket.getA(pos+1));
    }
    //cout<<"Computed tmpMgr.operR["<<pos<<"]"<<theMgr.operR[pos]<<endl;
    theMgr.storedR(pos);
  }
}



void Contractor::contractOperL(mwArray& result,const mwArray& oper,
			       const Site& ket,const Operator& op,
			       const Site& bra){
  op.contractL(result,oper,ket,bra);
}

void Contractor::contractOperR(mwArray& result,const mwArray& oper,
			       const Site& ket,const Operator& op,
			       const Site& bra){
  op.contractR(result,oper,ket,bra);
}

void Contractor::contractNormL(mwArray& result,const mwArray& norm,
			       const Site& ket){
  contractL(result,norm,ket,ket);
}

void Contractor::contractNormR(mwArray& result,const mwArray& norm,
			       const Site& ket){
  contractR(result,norm,ket,ket);
}

void Contractor::calculateMket(mwArray& result,const mwArray& tmpL,
			       const Site& ket,const Operator& op,
			       const mwArray& tmpR){
  //cout<<"Contractor::calculateMket calling contractMket on op "<<op<<endl;
  op.contractMket(result,tmpL,ket,tmpR);
}

void Contractor::calculateMbra(mwArray& result,const mwArray& tmpL,
			    const Site& bra,const Operator& op,
			    const mwArray& tmpR){
  //cout<<"Contractor::calculateMbra calling contractMbra on op "<<op<<endl;
  //const FoldedOperator& op2=dynamic_cast<const FoldedOperator&>(op);
  //op2.contractMbra(result,tmpL,bra,tmpR);
  op.contractMbra(result,tmpL,bra,tmpR);
}

void Contractor::calculateN(mwArray& result,const mwArray& tmpL,
			    const mwArray& tmpR,int d){
  int betab=tmpL.getDimension(0);int betak=tmpL.getDimension(1);
  int betabr=tmpR.getDimension(0);int betakr=tmpR.getDimension(1);
  result=reshape(permute(reshape(reshape(identityMatrix(d),Indices(d*d,1))*
				 reshape(reshape(tmpL,Indices(betab*betak,1))*
					 reshape(tmpR,Indices(1,betabr*betakr)),
					 Indices(1,betab*betak*betabr*betakr)),
				 Indices(d,d,betab,betak,betabr,betakr)),
			 Indices(1,3,5,2,4,6)),Indices(d*betab*betabr,d*betak*betakr));
}


void Contractor::contractL(mwArray& result,const mwArray& tmpL,
			   const Site& ket,const Site& bra){
  int d=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int betab=tmpL.getDimension(0);int betak=tmpL.getDimension(1);
#ifdef CHECKDIMS
  if((Dl!=betab)||(Dlp!=betak)||(d!=dp)){
    cout<<"Error in Contractor::contractL dimensions: bra="<<bra.getDimensions()
	<<", ket="<<ket.getDimensions()<<", tmpL="<<tmpL.getDimensions()<<endl;
    exit(212);
  }
#endif
  result=permute(ket.getA(),Indices(2,1,3));
  result.reshape(Indices(Dlp,dp*Drp));
  result=tmpL*result;
  result.reshape(Indices(betab*d,Drp));
  mwArray aux=permute(bra.getA(),Indices(3,2,1));
  aux.conjugate();aux.reshape(Indices(Dr,Dl*d));
  result=aux*result;
}

void Contractor::contractR(mwArray& result,const mwArray& tmpR,
			   const Site& ket,const Site& bra){
  int d=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int betab=tmpR.getDimension(0);int betak=tmpR.getDimension(1);
#ifdef CHECKDIMS
  if((Dr!=betab)||(Drp!=betak)||(d!=dp)){
    cout<<"Error in Contractor::contractR dimensions: bra="<<bra
	<<", ket="<<ket<<", tmpR="<<tmpR<<endl;
    exit(212);
  }
#endif
  result=permute(bra.getA(),Indices(2,1,3));
  result.conjugate();
  result.reshape(Indices(Dl*d,Dr));
  result=result*tmpR;
  result.reshape(Indices(Dl,d*betak));
  result=result*reshape(permute(ket.getA(),Indices(1,3,2)),Indices(dp*Drp,Dlp));
}

void Contractor::contractMketId(mwArray& result,const mwArray& tmpL,
				const Site& ket,const mwArray& tmpR){
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int betabr=tmpR.getDimension(0);int betakr=tmpR.getDimension(1);
  int betabl=tmpL.getDimension(0);int betakl=tmpL.getDimension(1);
#ifdef CHECKDIMS
  if((betakl!=Dlp)||(betakr!=Drp)){
    cout<<"Error in Contractor::contractMketId dimensions: tmpL"<<tmpL
	<<", ket="<<ket<<", tmpR="<<tmpR<<endl;
    exit(212);
}
#endif
  // First step: contract ket with tmpL
  mwArray tmpres=tmpL*reshape(permute(ket.getA(),Indices(2,1,3)),Indices(Dlp,dp*Drp));
  // Now add termR
  tmpres=reshape(tmpres,Indices(betabl*dp,Drp))*permute(tmpR,Indices(2,1));
  result=reshape(permute(reshape(tmpres,Indices(betabl,dp,betabr)),Indices(2,1,3)),
		 Indices(dp*betabl*betabr,1));
}

void Contractor::contractMbraId(mwArray& result,const mwArray& tmpL,
				const Site& bra,const mwArray& tmpR){
  int dp=bra.getd();int Dlp=bra.getDl();int Drp=bra.getDr();
  int betabr=tmpR.getDimension(0);int betakr=tmpR.getDimension(1);
  int betabl=tmpL.getDimension(0);int betakl=tmpL.getDimension(1);
#ifdef CHECKDIMS
  if((betabl!=Dlp)||(betabr!=Drp)){
    cout<<"Error in Contractor::contractMbraId dimensions: tmpL"<<tmpL
	<<", bra="<<bra<<", tmpR="<<tmpR<<endl;
    exit(212);
}
#endif
  // First step: contract bra (conj) with tmpR
  result=bra.getA();result.conjugate();
  result.reshape(Indices(dp*Dlp,Drp));
  result.multiplyRight(tmpR);
  result.reshape(Indices(dp,Dlp,betakr));
  result.permute(Indices(1,3,2));
  result.reshape(Indices(dp*betakr,Dlp));
  result.multiplyRight(tmpL);
  result.reshape(Indices(dp,betakr,betakl));
  result.permute(Indices(1,3,2));
  result.reshape(Indices(1,dp*betakl*betakr));
}

void Contractor::solveSite(const mwArray& N,const mwArray& M,mwArray& result,
			   complex_t& contraction1,
			   complex_t& contraction2,double svdtol_){

  //if(svdtol_==0) svdtol_=svdtol;
  if(N.isSquare())
    wrapper::lslu(N,M,result);
  else
    wrapper::lsd(N,M,result,svdtol_);
  // TODO!!! With the other lapack routine when convenient
  mwArray contr1=Hconjugate(result)*N*result;
  contraction1=contr1.getElement(0);
  mwArray contr2=Hconjugate(result)*M;
  contraction2=contr2.getElement(0);
}

void Contractor::solveSiteGauge(mwArray& result,mwArray& contraction1,
				mwArray& contraction2,
				const mwArray& M,double N){
  // Solves the matrix equation N*vector=M, when N is taken to be the
  // identity (times the constant N). It returns simply a vector
  // proportional to M, M/N and in the contractions,
  // contraction1=sol*sol and contraction2=sol*M
  //result=(1./N)*M;
  //cout<<"Contractor::solveSite receives M"<<M<<" and N "<<N<<endl;
  double normM=norm(M);
  result=(1./normM)*M;
  //contraction1=Hconjugate(result)*result;
  //contraction2=Hconjugate(result)*M;
  contraction1=mwArray(ONE_c);
  contraction2=mwArray((complex_t){normM,0.});
}

bool isnan_double(const double& g){return g!=g;}

bool Contractor::checkDistanceConvergence(double& distance,
					  const mwArray& contA,
					  const mwArray& contAB,
					  double normA,
					  double normB,int round){
  bool done;
  double newdistance;
  // Check convergence and set done=1, if needed
  //  cout<<"TEST: Checking convergence round "<<round<<", old distance "
  //  <<distance<<endl;
  // cout<<"Received : contr1="<<contA<<", contr2="<<contAB<<
  //   ", fact1="<<normA<<", normB ="<<normB<<endl;
  newdistance=normA*(normA*real(contA.getElement(0))-
		     2*real(contAB.getElement(0)))+normB*normB; // was ignoring normB
  if(isnan_double(newdistance)){
    cout<<"checkDistanceConvergence (round "<<round<<") received : distance="<<distance<<", contr1="
	<<contA<<", contr2="<<contAB<<
      ", fact1="<<normA<<", fact0 ="<<normB<<")"<<endl;
    cout<<"newdistance="<<newdistance<<endl;
    exit(1);
  }
  //cout<<"newdistance="<<newdistance<<", difference="<<newdistance-distance
  //     <<endl;
  if((newdistance-distance)>convtol*abs(distance)){
    if(abs(newdistance)>convtol&&abs(newdistance-distance)>1E-12){
      cout<<"Error, the norm distance (round "<<round<<") is increasing from the last step, and value is larger than "<<convtol<<"!"<<endl;
      cout<<"Received : distance="<<distance<<",  contr1="<<contA<<", contr2="<<contAB<<
	", fact1="<<normA<<", fact0 ="<<normB<<")"<<endl;
      cout<<" last="<<distance<<", new="<<newdistance<<
	", absolute difference "<<(newdistance-distance)<<
	", rel difference="<<(newdistance-distance)/distance
	  <<", tolerance="<<convtol<<endl;
      //cout<<"The norm of the state was "<<normA<<endl;
      exit(212);
    }
    else{
      // cout<<"WARNING: Assuming convergence because distance is smaller than "<<convtol<<"!"
      // 	  <<" Or difference is too small "<<endl;
      done=1;
    }
  }
  if(abs((distance-newdistance)/distance)<convtol){
    //cout<<"Contractor reached convergence after round "<<round<<
    //     "; rel. difference in distance="<<(distance-newdistance)/distance<<endl;
    done=1;
  }
  else{
    //     cout<<"Contractor did not reach convergence at end of round "<<round<<
    //       "; rel. difference in distance="<<(distance-newdistance)/distance<<endl;
    done=0;
  }
  // special case: if the distance I am looking at is small compared to convtol, then the criterion above may not work.
  if((abs(distance)+abs(newdistance))*convtol<1E-14){
    // cout<<"WARNING: Contractor::checkDistanceConvergence is trying to compare small numbers=>turning to absolute values"<<endl;
    if(abs(distance-newdistance)<convtol)
      done=1;
  }
  distance=newdistance;
  return done;
}

complex_t Contractor::contractR(const MPS& ket,const MPO& ops,
			      const MPS& bra){
  //  cout<<"Contractor:contractR("<<bra<<","<<ops<<","<<ket<<")"<<endl;
  int nrsites=ket.getLength();
  if((nrsites!=bra.getLength())||(nrsites!=ops.getLength())){
    cout<<"Error in contractR: dimensions do not coincide: bra ";
    cout<<bra.getLength()<<", ket "<<nrsites<<", ops "<<
      ops.getLength()<<endl;
    exit(212);
  }
  // TODO: Check dimensions!
  mwArray result(ONE_c); //initial value
  for(int k=nrsites-1;k>=0;k--){
    //cout<<"Before step "<<k<<", result="<<result<<endl;
    mwArray tmp(result);
    contractOperR(result,tmp,ket.getA(k),ops.getOp(k),bra.getA(k));
  }
  // At the end: check size (1x1) and return value!
  if(!result.isScalar()){
    cout<<"Error: after Contractor::contractR(MPS,MPO,MPS), result "
	<<result<<endl;
    exit(212);
  }
  // if the given MPSs have normalization factors, I must multiply them now
  return (ket.getNormFact()*bra.getNormFact())*result.getElement(0);
}


// TODO: Set max number of iterations for these ops?
void Contractor::findRightEigenvector(const MPO& ops,int D,
				      double* lambda,MPS& resultMPS){
  return findMaxEigenvector(ops,D,lambda,'R',resultMPS);
}
void Contractor::findLeftEigenvector(const MPO& ops,int D,
				     double* lambda,MPS& resultMPS){
  return findMaxEigenvector(ops,D,lambda,'L',resultMPS);
}

void Contractor::findMaxEigenvector(const MPO& ops,int D,
				    double* lambda,char dir,MPS& lastMPS){
  if(lastMPS.isEmpty()){
    lastMPS=MPS(ops.getLength(),D,ops.getDimensions());
    lastMPS.setProductState(p_zero);
  }
  else{ // Check dimensions
    if(lastMPS.getLength()!=ops.getLength()){
      cout<<"Error: invalid sizes in Contractor::findMaxEigenvector, "
	  <<"ops "<<ops<<", lastMPS "<<lastMPS<<endl;
      exit(212);
    }
    else if(D!=lastMPS.getBond())
      lastMPS.increaseBondDimension(D);
    // I need to check also the PHYSICAL dimensions!
    vector<int> findims=ops.getDimensions();
    lastMPS.adjustPhysDimensions(findims);
  }
  lastMPS.gaugeCond('R',1);lastMPS.gaugeCond('L',1);
  // Repeat the following procedure:
  // (1) Apply Operator -> find new (unnormalized) MPS
  // (2) Check if normalization factor converges
  // (3) If not, normalize obtained state and repeat from (1)
  bool done=0;int count=0;
  *lambda=0.; // initialize
  int nrkicks=0,perturbcnt=0; // how many times I have retried and perturbed (todo: set a limit?)
  // cout<<"Contractor find MaxEigenvector("<<dir<<") with initial MPS "
  //     <<lastMPS<<", and operator "<<ops<<endl;
  while(!done){
    count++; // count number of attempts => set a limit
    //cout<<"Contractor::findMaxEigenvector"<<dir<<"["<<count<<"] lamb="
       // 	<<*lambda<<"; "<<endl;
    // 6-10-2008 Try with approximation!!!
    //cout<<"Copying the MPS "<<lastMPS<<" into aux"<<endl;
    MPS aux(lastMPS); // on which ops is acting (constant after optimize)
    if(count<2){
      switch(dir){
      case 'R':{
	lastMPS.approximate(ops,aux,D);
	break;
      }
      case 'L':{
	lastMPS.approximate(ops,aux,D,'U'); // op acting upwards
	break;
      }
      }
      //cout<<"First round: approximation for initial guess"<<endl;
    }
    switch(dir){
    case 'R':
      optimize(ops,aux,lastMPS);
      break;
    case 'L':
      optimizeL(ops,aux,lastMPS);
      break;
    }
    
    double newLamb=lastMPS.getNormFact();
    // At the end, check convergence
    if(abs(newLamb-*lambda)<abs((*lambda)*convtol)){
      cout<<"Contractor::find "<<dir<<" MaxEigenvector converged after "<<count
      	  <<" attempts with lambda="<<newLamb<<" and tol="<<convtol<<endl;
      done=1;
    }
    else{
      //cout<<"last var="<<abs(newLamb-*lambda)/abs(*lambda)<<"; ";
      if(count>MAXITER){
	cout<<"No convergence of findMaxEigenvector after "<<count<<" trials"<<endl;
	lastMPS.perturb(.01); // Perturbation on 1 out of 100 sites
	count=0;
	perturbcnt++;
	if(perturbcnt>=MAXPERT){
	  cout<<"Contract could not find an eigenvector after "<<perturbcnt
	      <<" trials (each of them = perturbation-"<<MAXITER<<" applications"
	      <<" of the operator=> FAILED"<<endl;
	  exit(1); // TODO: Seguir de todas formas??
	}
      }    
    }
    *lambda=newLamb;
    lastMPS.gaugeCond('R',1);lastMPS.gaugeCond('L',1); 
    // Check again?
     // cout<<"Now <vec "<<dir<<"|O|old vec>="<<contract(aux,ops,lastMPS)<<", and <vec|vec>="
//  	<<contract(lastMPS,lastMPS)<<", <oldV|oldV>="<<contract(aux,aux)<<", <vec|oldV>="
//  	<<contract(aux,lastMPS)<<endl;
  }

}

#include "time.h"

mwArray Contractor::getEffectiveOperatorSingleSite(const MPS& ket,const MPO& ops,int pos){
  int nrsites=ops.getLength();  
  if(pos<0||pos>=nrsites){
    cout<<"Position "<<pos<<" out of range in Contractor::getEffectiveOperatorSingleSite "<<endl;
    exit(1);
  }
  // Make a copy of the MPS to ensure the proper gauge conditions (TODO!!!! CHECK!!!)
  MPS auxMPS(ket);
  //  double origNorm=contract(auxMPS,auxMPS); // keep the original normalization, just in case?
  if(!auxMPS.isGaugeR()) auxMPS.gaugeCond('R',1);
  // Now left part is normalized => contract right part until pos
  for(int k=nrsites-1;k>=pos+1;k--){
    auxMPS.gaugeCond(k,'L',true);
  }
  // Now N matrix is identity
  //##################### HERE!!!!
  bool gauge=true;
  bool ketN=false; 
  TmpManager tmpMgr(nrsites,gauge);
  // left and right terms
  calculateL(pos,tmpMgr,ops,auxMPS,auxMPS,gauge,ketN);
  calculateR(pos,tmpMgr,ops,auxMPS,auxMPS,gauge,ketN);
  mwArray M; // the effective operator
  ops.getOp(pos).contractN(M,tmpMgr.operL[pos],tmpMgr.operR[pos]);
  return M;
}


mwArray Contractor::getEffectiveOperatorMultiSite(const MPS& ket,const MPO& ops,int k,int pos){
  mwArray M; // the effective operator
  // I reuse the work done for the TensorMultiplier
  TensorMultiplier aux=getEffectiveOperatorMultiplierMultiSite(ket,ops,k,pos);
  // Just use the full tensor from TensorMultiplier
  aux.getFullTensor(M);
  return M;
  // A trick: I create an operator with the contraction of all intermediate ones, 
  // to finish this contraction
  // Operator auxOp(aux.getMiddleTerm());
  // auxOp.contractN(M,aux.getLeftTerm(),aux.getRightTerm());
  // return M;
}

TensorMultiplier Contractor::getEffectiveOperatorMultiplierSingleSite(const MPS& ket,const MPO& ops,int pos){
  return getEffectiveOperatorMultiplierMultiSite(ket,ops,1,pos);
}

void Contractor::getEffectiveOperatorMPOMultiplier(const MPS& ket,const MPO& ops,int block,int pos,MPOMultiplier& mpoMulti,bool gauged){
  //cout<<"Contractor::getEffectiveOperatorMPOMultiplier from mps("<<ket.getLength()<<") and mpo("<<ops.getLength()<<") sites "
  //   <<pos<<" to "<<pos+block-1<<endl;
  int nrsites=ops.getLength();  
  if(pos<0||pos+block-1>=nrsites){
    cout<<"Position "<<pos<<" out of range in Contractor::getEffectiveOperatorSingleSite "<<endl;
    exit(1);
  }
  mpoMulti.initLength(block+2);

  bool gauge=true; // I don't want to compute Neff
  bool ketN=false; // neither some other term
  TmpManager tmpMgr(nrsites,gauge);
  // First, compute left and right terms

  // To ensure an orthogonal basis, I would need to impose gauge
  // conditions up to the edges of the operator, but if gauged is
  // true, I assume it is done outside
  if(!gauged){
     MPS auxMPS(ket);
     auxMPS.gaugeCond('L',0); // norm factor should be included
     for(int p=0;p<pos-1;p++){
       auxMPS.gaugeCond(p,'R',true); // apply gauge to right until pos-1
     }
     if(block>0){ // then I can apply the gauge properly, such taht the state does not change (it does not really matter, though)
       auxMPS.gaugeCond(pos-1,'R',true); // apply gauge to right until pos-1
     }
     else{
       auxMPS.gaugeCond(pos-1,'R',false); // If block length is 0, I cannot push the remianing tensor, or I will spoil the gauge to left of pos
     }
     calculateL(pos,tmpMgr,ops,auxMPS,auxMPS,gauge,ketN);
     calculateR(pos+block-1,tmpMgr,ops,auxMPS,auxMPS,gauge,ketN);
  }
  else{
    calculateL(pos,tmpMgr,ops,ket,ket,gauge,ketN);
    calculateR(pos+block-1,tmpMgr,ops,ket,ket,gauge,ketN);
  }

  mwArray& opL=tmpMgr.operL[pos];
  Indices dimsL=opL.getDimensions(); // bra,xi,ket
  mwArray& opR=tmpMgr.operR[pos+block-1];
  Indices dimsR=opR.getDimensions(); // bra,xi,ket
  mpoMulti.setOp(0,new Operator(permute(reshape(opL,Indices(dimsL,1)),Indices(1,4,3,2))),true);
  mpoMulti.setOp(block+1,new Operator(reshape(opR,Indices(dimsR,1))),true);
  for(int k=0;k<block;k++){
    //    cout<<"Contractor::getEffectiveOperatorMPOMultiplier setting op "<<k+1<<endl;
    mpoMulti.setOp(k+1,&ops.getOp(pos+k),false);
  }
  //cout<<"Constructed the MPOMultiplier "<<mpoMulti<<endl;
}


TensorMultiplier Contractor::getEffectiveOperatorMultiplierMultiSite(const MPS& ket,const MPO& ops,int k,int pos){
  cout<<"Contractor::getEffectiveOperatorMultiplierMultiSite from mps("<<ket.getLength()<<") and mpo("<<ops.getLength()<<") sites "<<pos<<" to "<<pos+k-1<<endl;
  int nrsites=ops.getLength();  
  if(pos<0||pos+k-1>=nrsites){
    cout<<"Position "<<pos<<" out of range in Contractor::getEffectiveOperatorSingleSite "<<endl;
    exit(1);
  }
  // I am ensuring that the operator is in orthogonal basis, so need to impose a couple of gauge conditions (only on the edges)
  MPS auxMPS(ket);
  // //  double origNorm=contract(auxMPS,auxMPS); // keep the original normalization, just in case?
  // if(!auxMPS.isGaugeR()) auxMPS.gaugeCond('R',1);
  // // Now left part is normalized => contract right part until pos+k-1
  // for(int ik=nrsites-1;ik>=pos+k;ik--){
  //   auxMPS.gaugeCond(ik,'L',true);
  // }
  auxMPS.gaugeCond(pos-1,'R',true);
  auxMPS.gaugeCond(pos+k,'L',true);

  // Now N matrix is identity
  //##################### HERE!!!!
  bool gauge=true;
  bool ketN=false; 
  TmpManager tmpMgr(nrsites,gauge);
  // left and right terms
  calculateL(pos,tmpMgr,ops,auxMPS,auxMPS,gauge,ketN);
  calculateR(pos+k-1,tmpMgr,ops,auxMPS,auxMPS,gauge,ketN);
  mwArray tmpOp;
  blockOperatorMultiSite(tmpOp,ops,k,pos);
  cout<<"Will construct a multiplier out of tmpOp:"<<tmpOp.getDimensions()<<", left "<<tmpMgr.operL[pos].getDimensions()<<"and right "<<tmpMgr.operR[pos+k-1].getDimensions()<<endl;
  return TensorMultiplier(tmpMgr.operL[pos],tmpOp,tmpMgr.operR[pos+k-1]);
}

void Contractor::sweepPart(int pos0,int lim,const MPO& ops,int D,
			   MPS& init,double offset,int knr){
  int nrsites=ops.getLength();  
  vector<int> findims=ops.getDimensions();
  init.adjustPhysDimensions(findims);
  bool right=lim>=pos0;
  const char dir=right?'R':'L';
  // const char dir2=right?'L':'R'; // the other direction
  // init.gaugeCond(dir,1); // this first gauge ensures right D dimensions
  // init.gaugeCond(dir2,1); /* this ensures proper D and allows starting 
  // 			     from appropriate edge.*/

  if(D>0)
    if(D>init.getBond())
      init.increaseBondDimension(D);
    else if(D<init.getBond()){
      cout<<"WARNING: Contractor::sweepPart called with D("<<D<<") smaller "
	  <<"than in the starting MPS (D_init="<<init.getBond()<<"). Will truncate "
	  <<"the initial MPS!"<<endl;
      init=MPS(init,D);
    }
  // I will assume that dimensions are ok from the beginning
  int pos=pos0;
  for(int k=0;k<pos;k++)
    init.gaugeCond(k,'R',true);
  for(int k=nrsites-1;k>pos;k--)
    init.gaugeCond(k,'L',true);


  double start=0.,finish=0.;
  bool gauge=true;bool ketN=false;
  bool done=false; 
  TmpManager tmpMgr(nrsites,gauge);
  double energy=real(contract(init,ops,init)); // initial value
  double val=energy;
  while(!done){
    // finish=clock();
    // cout<<"Sweep part, pos "<<pos<<", val="<<val<<", going in dir "<<dir;
    // cout<<" time (previous round)"<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
    // start=clock();
    // 1. Obtain left and right terms for everything
    calculateL(pos,tmpMgr,ops,init,init,gauge,ketN);
    calculateR(pos,tmpMgr,ops,init,init,gauge,ketN);
 
    // Place for the results of the diagonalization
    vector<complex_t> diagVals;
    mwArray U(init.getA(pos).getA()); 
    U.reshape(Indices(-1,1));
    mwArray tmpOp=ops.getOp(pos).getFullData();
    //TensorMultiplierHermitianOffset multi(tmpMgr.operL[pos],tmpOp,tmpMgr.operR[pos],offset);
    TensorMultiplier* multi;
    if(solver==fulleig)
      multi=new TensorMultiplierHermitianOffset(tmpMgr.operL[pos],tmpOp,tmpMgr.operR[pos],offset);
    else
      multi=new TensorMultiplierOffset(tmpMgr.operL[pos],tmpOp,tmpMgr.operR[pos],offset);
    int knr_=min(knr,multi->getSize());
    // if(abs(energy+offset)<convtol){
    //   // Estimate largest eigenvalue
    //   solveEigs(*multi,knr_,"LM",diagVals,U,false); 
    //   cout<<"Since I hold a small value of E("<<energy<<")+offset("<<offset
    // 	  <<"), and I suspect the condition "
    // 	  <<"of the matrix could be bad, will use an offset ";
    //   offset+=-2.5*real(*(complex_t*)&diagVals[0]); // first eigenvalue
    //   cout<<offset<<endl;
    //   TensorMultiplierHermitianOffset multiOff(tmpMgr.operL[pos],tmpOp,tmpMgr.operR[pos],offset);
    //   solveEigs(multiOff,knr_,"SR",diagVals,U,true); // Compute one eigenvalue-vector
    // }
    // else{    
    solveEigs(*multi,knr_,"SR",diagVals,U,true); // Compute one eigenvalue-vector
    // }
    delete multi;
    val=real(*(complex_t*)&diagVals[0])-offset; // first eigenvalue
    //cout<<"After solving Eigenvalue Eq, lambda="<<val<<endl;
    U.resize(Indices(U.getDimension(0),1)); // first column only!
    init.setA(pos,U); 
    bool passOnX=false; // whether to multiply next tensor by the transformation
    if((right&&pos==lim-1)||(!right&&pos==lim+1)||(pos==lim))
      passOnX=true;
    //    cout<<"About to apply gauge "<<dir<<" to "<<pos<<" and pass TErm="<<passOnX<<endl;
    init.gaugeCond(pos,dir,passOnX);
    // notify temporary storage as appropriate
    for(int l=pos+1;l<nrsites;l++)
      tmpMgr.outofdateL(l);
    for(int l=pos-1;l>=0;l--)
      tmpMgr.outofdateR(l);
    if(right) pos++;
    else pos--;
    done=right?pos>=lim:pos<=lim;
  }    

}

void Contractor::findGroundState(const MPO& ops,int D,double* lambda,
				 MPS& init,double offset,int knr,
				 const string& tmpfile,int freqSv,bool useLM){
#ifndef USING_PRIMME
  cout<<"Compiled without PRIMME support => eigensolver = ARPACK"<<endl;
  solver=arpack;
#endif
  int nrsites=ops.getLength();  
  vector<int> findims=ops.getDimensions();
  init.adjustPhysDimensions(findims);
  if(D>0)
    if(D>init.getBond())
      init.increaseBondDimension(D);
    else if(D<init.getBond()){
      cout<<"WARNING: Contractor::findGroundState called with D("<<D<<") smaller "
	  <<"than in the starting MPS (D_init="<<init.getBond()<<"). Will truncate "
	  <<"the initial MPS!"<<endl;
      init=MPS(init,D);
    }
  init.gaugeCond('R',1); // this first gauge ensures right Dr dimensions
  init.gaugeCond('L',1); /* this ensures proper Dl and allows starting 
			    from left.*/
  bool done=0;int round=0;
  bool gauge=true;
  bool right=1; int pos=0; //nrsites-1;
  *lambda=0.; // initialize
  bool ketN=false; 
  TmpManager tmpMgr(nrsites,gauge);
  double energy=real(contract(init,ops,init)); // initial value
  clock_t start=clock();clock_t finish=clock();
  double val=energy;
  cout<<"Starting findGroundState with initial value E="<<energy<<endl;
  while(!done){
    if(pos==0){
      finish=clock();
      cout<<"Round number "<<round<<", pos "<<pos<<", energy "<<val<<" (time "<<
	(finish-start)*1./CLOCKS_PER_SEC<<")"<<endl; //energy;
      start=clock();
    }
    // 1. Obtain left and right terms for everything
    calculateL(pos,tmpMgr,ops,init,init,gauge,ketN);
    calculateR(pos,tmpMgr,ops,init,init,gauge,ketN);
    // Place for the results of the diagonalization
    vector<complex_t> diagVals;
    // mwArray U;
    mwArray U(init.getA(pos).getA()); 
    // todo: include gauge from previous!
    U.reshape(Indices(-1,1));//U.fillWithOne();
    mwArray tmpOp=ops.getOp(pos).getFullData();
    //TensorMultiplierHermitianOffset multi(tmpMgr.operL[pos],tmpOp,tmpMgr.operR[pos],offset);
    TensorMultiplier* multi;
    if(solver==fulleig)
      multi=new TensorMultiplierHermitianOffset(tmpMgr.operL[pos],tmpOp,tmpMgr.operR[pos],offset);
    else
      multi=new TensorMultiplierOffset(tmpMgr.operL[pos],tmpOp,tmpMgr.operR[pos],offset);
    int knr_=min(knr,multi->getSize());
    const char* which="SR";
    if(useLM) which="LM";
    solveEigs(*multi,knr_,which,diagVals,U,true); // Compute one eigenvalue-vector
    delete multi;
    //if(useLM)      cout<<"Found the LM ev of this thing: "<<diagVals<<endl;
    val=real(*(complex_t*)&diagVals[0])-offset; // first eigenvalue
    //cout<<"After solving Eigenvalue Eq, ("<<diagVals[0]<<") lambda="<<val<<endl;
    U.resize(Indices(U.getDimension(0),1)); // first column only!
    init.setA(pos,U); 
    //if(val-energy>abs(energy)*1E-10){
    if((val-energy>1E-8)&&(val-energy>abs(energy+offset)*1E-10)){
      if((val-energy>abs(energy+offset)*convtol)){
	cout<<"Error: energy increasing in this step!! oldenergy="<<energy
	    <<", new one="<<val<<", difference="<<(val-energy)
	    <<", relative increase (incl offset="<<offset<<")="<<(val-energy)/(energy+offset)
	    <<", round="<<round<<", pos="<<pos
	    <<endl;
	exit(212);
      }
      else{
	cout<<"WARNING: energy increasing in this step!! oldenergy="<<energy
	    <<", new one="<<val<<", difference="<<(val-energy)
	    <<", relative increase="<<(val-energy)/(energy+offset)
	    <<", round="<<round<<", pos="<<pos
	    <<" I CONTINUE!!"<<endl;
      }
    }
    if(right){
      //if(pos<nrsites-1)
      init.gaugeCond(pos,'R',true); // I transform the next site to be able to use the tensor as starting point
    }
    else{
      //if(pos>0)
      init.gaugeCond(pos,'L',true);
    }
    // notify temporary storage as appropriate
    for(int l=pos+1;l<nrsites;l++)
      tmpMgr.outofdateL(l);
    for(int l=pos-1;l>=0;l--)
      tmpMgr.outofdateR(l);
    
    bool doneRound=(right&&pos==nrsites-1)||((!right)&&pos==0); // was this an endo of round?
    //if (pos==1) done=1;
    // Now move one pos and check convergence if we are at the end f
    // one sweep
    if(right){ //moving right
      if(pos==nrsites-1){ //last reached: check convergence and change sense
	right=0;
	done=abs(energy-val)<convtol*abs(energy+offset);
	//done=abs(energy-val)<convtol*abs(energy);
	if(!done){
	  //	  cout<<"Done test failed (R) with relative change "
	  //  <<abs((energy-val)/(energy+offset))<<" including offset "
	  //  <<offset<<", with val="<<val<<endl;
	  // Check if the value is too close to zero!
	  if(abs(energy-val)<1E-13){
	    cout<<"WARNING: Absolute value of the difference below 1E-13"<<endl;
	    //cout<<"Assuming convergence because difference and RELATIVE(!!) difference is small!"<<endl;
	    //done=true;
	  }
	  if((abs(val)<convtol)&&(abs(energy)<convtol)){
	    cout<<"WARNING: Absolute value of the energy, in this step and previous one, below tol"<<endl;
	    //cout<<"Assuming convergence because values are smaller than tolerance!"<<endl;
	    //done=true;
	  }
	  if(abs(val)<1E-13){
	    //cout<<"WARNING: Absolute value of the energy below 1E-13"<<endl;
	    //cout<<"WARNING***** Assuming convergence because E is small!"<<endl;
	    cout<<"WARNING***** Stopping because E is small! (<1E-13)"<<endl;
	    done=true;
	  }
	}
	energy=val;
	if(!done){
	  init.gaugeCond(pos,'L',true);pos--;
	}
	round++;
      }
      else{
	pos++;
      }
    }
    else{// moving left
      if(pos==0){ //first reached: check convergence and change sense
	right=1;
	done=abs(energy-val)<convtol*abs(energy+offset);
	//done=abs(energy-val)<convtol*abs(energy);
	if(!done){
	  //	  cout<<"Done test failed (L) with relative change "
	  //  <<abs((energy-val)/(energy+offset))<<" including offset "
	  //  <<offset<<", with val="<<val<<endl;
	  // Check if the value is too close to zero!
	  if(abs(energy-val)<1E-13){
	    cout<<"WARNING: Absolute value of the difference below 1E-13"<<endl;
	    //cout<<"Assuming convergence because difference is small!"<<endl;
	    //done=true;
	  }
	  if((abs(val)<convtol)&&(abs(energy)<convtol)){
	    cout<<"WARNING: Absolute value of the energy, in this step and previous one, below tol"<<endl;
	    //cout<<"Assuming convergence because values are smaller than tolerance!"<<endl;
	    //done=true;
	  }
	  if(abs(val)<1E-13){
	    //cout<<"WARNING: Absolute value of the energy below 1E-13"<<endl;
	    cout<<"WARNING***** Stopping because E is small! (<1E-13)"<<endl;
	    done=true;
	  }
	}
	energy=val;
	if(!done){
	  init.gaugeCond(pos,'R',true);pos++;
	}
	round++;
      }
      else{
	pos--;
      }
    }
    if(!done&&doneRound&&tmpfile!="") // Saving intermediate results (only at the end)
      if(round%freqSv==0){
	init.exportMPS(tmpfile.data());
	cout<<"Saved Temporary resuts to "<<tmpfile<<endl;
      }
  }
  *lambda=energy;
}



/*************** TEST FUNCTION ************************* */
#ifdef TESTINGCNTR
void Contractor::test(const MPO& ops,const MPS& psi1,
			  const MPS& psi2){

  int nrsites=psi1.getLength();
  int nrsites2=psi2.getLength();
  int nrsiteso=ops.getLength();
  if((nrsites!=nrsites2)||(nrsites!=nrsiteso)){
    cout<<"Error en las dimensiones!"<<endl;
    exit(2);
  }

  mwArray refSc,refOp; // reference values
  refSc=contract(psi2,psi2);
  refOp=contract(psi1,ops,psi2);


  /// HERE: Initialize temporary storage
  TmpManager* tmpMgptr=new TmpManager(nrsites);
  TmpManager& tmpMgr=*tmpMgptr;
  bool gauge=0;bool ketN=0;
  for(int k=0;k<nrsites;k++){
    int pos=k;
    Indices dims=psi2.getA(pos).getDimensions();
    cout<<"En el test, pos="<<pos<<", dims(A)="<<dims<<endl;
    int d=dims[0];
    calculateL(k,tmpMgr,ops,psi1,psi2,gauge,ketN);
    calculateR(k,tmpMgr,ops,psi1,psi2,gauge,ketN);
    //cout<<"Now checking tmpMgr "<<tmpMgr<<endl;
    mwArray M; mwArray N;
    calculateMket(M,tmpMgr.operL[pos],psi1.getA(pos),ops.getOp(pos),
		  tmpMgr.operR[pos]);
    //calculateMbra(M,tmpMgr.getOperL(pos),psi2.getA(pos),ops.getOp(pos),
    //	  tmpMgr.getOperR(pos));
    calculateN(N,tmpMgr.normL[pos],tmpMgr.normR[pos],d);

    // For calculateMket
    //m_vector(1,Avec2,psi2.getA(pos).getA());
    //m_adjoint(1,Atransp2,Avec2);
    //m_times(1,resultM,Atransp2,M);
    //       m_vector(1,Avec1,psi1.getA(pos).getA());
    //       m_adjoint(1,Atransp1,Avec1);
    //       m_times(1,resultN,Atransp2,N);
    //       m_times(1,resultN,resultN,Avec2);
    // For calculateMbra:
    mwArray Avec2=psi1.getA(pos).getA();
    Avec2.reshape(Avec2.getNrComponents(),1);
    mwArray resultM=Avec2*M;
    mwArray Avec1=psi2.getA(pos).getA();
    Avec1.reshape(Avec1.getNrComponents(),1);
    mwArray resultN=conjugate(Avec1)*N*Avec1;
    
    cout<<"In position "<<k<<", result from M ="<<resultM<
      <", reference="<<refOp<<endl;
    cout<<"In position "<<k<<", result from N ="<<resultN
	<<", reference="<<refSc<<endl;

  }
  delete tmpMgptr;
}
#endif // TESTINGCNTR


/** */

void Contractor::optimizeInv(const MPO& ops,const MPS& orig,
			     MPS& init,int D,double* errN,double* errT){
  // I approximate by an MPS satisfying gauge condition, times a 
  // normalization factor. 
  
  vector<int> findims=ops.getOriginalDimensions();
  init.adjustPhysDimensions(findims);
  if(D>0)
    if(D>init.getBond())
      init.increaseBondDimension(D);
    else if(D<init.getBond()){
      cout<<"WARNING: Contractor::optimizeInv called with D("<<D<<") smaller "
	  <<"than in the starting MPS (D_init="<<init.getBond()<<"). Will truncate "
	  <<"the initial MPS!"<<endl;
      init=MPS(init,D);
    }
  init.gaugeCond('R',true); /* this first gauge ensures right Dr dimensions */
  init.gaugeCond('L',true); /* this ensures proper Dl and allows starting from left.*/
  //  if(D>0)
  //init.increaseBondDimension(D);

  double origNorm=real(contract(orig,orig));
  double constantN=origNorm/(orig.getNormFact()*orig.getNormFact()); // to be sure it works even without gauge
  int nrsites=orig.getLength();

  // Initialize temporary storage with gauge condition=false (trick 
  // because we use norm now for other quantity)
  TmpManager tmpMgr(nrsites,0);

  // Start optimization: sweep to right and left until convergence
  bool done=0;
  bool right=1;
  int pos=0;
  double distance=1000*orig.getNormFact()*orig.getNormFact(); // some crazily large initial value
  double olddistance=distance; // from iteration to iteration
  int round=1;

  double mytol=1E-10;

  mwArray M;
  mwArray N,newA;
  complex_t contr1,contr2;
  while(!done){
    /* Optimize matrix for position pos */
    //    cout<<"optimizeInv round "<<round<<", pos="<<pos<<endl;
    //int d=init.getA(pos).getd(); // ignored?
    // Calculate M(pos) and N(pos) matrices 
    // 1. Obtain left and right terms for everything
    calculateL2(pos,tmpMgr,ops,orig,init);
    calculateR2(pos,tmpMgr,ops,orig,init);
    // 2. Contract the terms appropriately for M and N matrices
    bool dagger=true;
    ops.getOp(pos).contractMket(M,tmpMgr.operL[pos],orig.getA(pos),
				tmpMgr.operR[pos],dagger);
    // And the N term!!
    ops.getOp(pos).contractN2(N,tmpMgr.normL[pos],tmpMgr.normR[pos]);
    if(pos==0&&round==1){
      if(M==mwArray(ZERO_c)){
	cout<<"Warning: perturbing the initial MPS "<<endl;
	init.perturb(1.);
	for(int l=0;l<nrsites;l++){
	  tmpMgr.outofdateL(l);
	  tmpMgr.outofdateR(l);
	}
	continue;
      }
    }
    // Since my MPS contains a norm factor, which is not included 
    // by the contractions, I have to multiply it by hand
    N=(init.getNormFact()*init.getNormFact())*N;
    M=init.getNormFact()*M;
    // 3. Solve the linear equations
    //solveSite(N,M,newA,contr1,contr2,mytol); // Solve Matrix Equation (N not Id)
    solveSite(N,M,newA,contr1,contr2,0); // Solve Matrix Equation (N not Id)
    //if(round==1)
    //cout<<"Linear eq receives N="<<N<<", M="<<M<<", solution "<<newA<<endl;   
    //cout<<"Linear eq for pos "<<pos<<" gives C1="<<contr1<<", C2="<<contr2<<endl;
    // Computed with A before imposing gauge condition => the one solving the equation
    double newdist=constantN+real(contr1)-2*real(contr2);
    cout<<"Contractor::optimizeInv pos "<<pos<<", round "<<round<<", distance=="<<newdist<<endl;

    /** Now we have to substitute the solved A, but since I do not trust NA=M solution (numerical errors), if I get a distance 
	larger than before, I do not change the state. */
    if(newdist<distance){
      init.setA(pos,newA);
      distance=newdist;
    }
    else{
      cout<<"NOT Changing tensor at site "<<pos<<" as distance was "<<distance<<" and now is "<<newdist<<endl;
      newdist=distance;
    }
    // now impose gauge condition (according to the direction) 
    if(right){
      // the option 1 is required only in the edge, to keep norm factor, but
      // I multiply the neighbour, in case it does not change next time
      init.gaugeCond(pos,'R',true);  
      if(pos+2<nrsites)
	tmpMgr.outofdateL(pos+2);     
    }
    else{
      init.gaugeCond(pos,'L',true);
      if(pos-2>=0)
	tmpMgr.outofdateR(pos-2);     
    }
    // notify temporary storage as appropriate
    for(int l=pos+1;l<nrsites;l++)  // if gauge without changing next site
      tmpMgr.outofdateL(l);     
    for(int l=pos-1;l>=0;l--) // if gauge without changing next site
      tmpMgr.outofdateR(l);

    //    cout<<"After changing tmpMgr at pos "<<pos<<endl;

    // Decide the convergence only at the end of sweep
    if(((right)&&(pos==nrsites-1))||
       ((!right)&&(pos==0))){
      done=(abs(olddistance-distance)<max(convtol*abs(olddistance),NUMPREC))?1:0;
      if(!done){  // turn around and continue
	if(right){
	  // init.gaugeCond(pos,'R',1); // this sets the norm factor but is already applied
	  init.gaugeCond(pos,'L',true);pos--;
	  tmpMgr.outofdateR(pos);
	  if(pos-1>=0)
	    tmpMgr.outofdateR(pos-1);
	  if(pos+1<nrsites)
	    tmpMgr.outofdateL(pos+1);
	}
	else{
	  // init.gaugeCond(pos,'L',1); // this sets the norm factor but is already applied
	  init.gaugeCond(pos,'R',true);pos++;
	  tmpMgr.outofdateL(pos);
	  if(pos+1<nrsites)
	    tmpMgr.outofdateL(pos+1);
	  if(pos-1>=0)
	    tmpMgr.outofdateR(pos-1);
	}
	olddistance=distance;
	right=!right;
	round++;
      }
      }
    else{ // not in the edge=> just advance in the same direction
      pos=right?pos+1:pos-1;
    }       
  } // end of while(!done)
  // When the state has converged, I have something which is not normalized.
  // In any case, if the original vector had a normFactor, it was ignored 
  // in the above calculation, hence it has to be added by hand now!
  //cout<<"At the end of optimize, changing norm "<<init.getNormFact();
  init.setNormFact(init.getNormFact()*orig.getNormFact());
  // If error has to be computed...
  if(errN!=0){
    double norm0=real(contract2(ops,init));
    double norm1=real(contract(orig,orig));
    complex_t crossR=contract(init,ops,orig);
    *errN=norm0+norm1-2*real(crossR);
    if(*errN<-1E-12){
      cout<<"WARNING: Contractor finds err^2<0!!!"<<endl;
    }
    if(errT!=0){
      *errT=norm0;
    }
  }
  tmpMgr.clear();
}



void Contractor::calculateL2(int pos,TmpManager& theMgr,const MPO& ops,
			     const MPS& ket,const MPS& bra,bool resolvent){
  //cout<<"Contractor::calculateL2(pos="<<pos<<", tmpMgr="<<&theMgr<<")"<<endl;
  if(!theMgr.isuptodateL(pos)){ //else, nothing to do, just return!
    calculateL2(pos-1,theMgr,ops,ket,bra,resolvent); // take the one before
    //cout<<"Using tmpL of pos "<<pos-1<<" uptodate"<<endl;
    ops.getOp(pos-1).contractL2(theMgr.normL[pos],theMgr.normL[pos-1],
				bra.getA(pos-1),bra.getA(pos-1));
    bool dagger=true;
    ops.getOp(pos-1).contractL(theMgr.operL[pos],theMgr.operL[pos-1],
			       ket.getA(pos-1),bra.getA(pos-1),dagger);
    //TmpManagerResolvent* ptr=dynamic_cast<TmpManagerResolvent*>(&theMgr);
    if(resolvent){
      //ops.getOp(pos-1).contractL(ptr->normSL[pos],ptr->normSL[pos-1],
      //		       bra.getA(pos-1),bra.getA(pos-1),dagger);
      contractL(theMgr.operRL[pos],theMgr.operRL[pos-1],
		ket.getA(pos-1),bra.getA(pos-1));
    }
    theMgr.storedL(pos);    
  }
}

void Contractor::calculateR2(int pos,TmpManager& theMgr,const MPO& ops,
			     const MPS& ket,const MPS& bra,bool resolvent){
  //cout<<"Contractor::calculateR(pos="<<pos<<", tmpMgr="<<&theMgr<<")"<<endl;
  if(!theMgr.isuptodateR(pos)){ //else, nothing to do, just return!
    calculateR2(pos+1,theMgr,ops,ket,bra,resolvent); // take next one and use it!
    //cout<<"Using tmpR of pos "<<pos+1<<" uptodate"<<endl;
    ops.getOp(pos+1).contractR2(theMgr.normR[pos],theMgr.normR[pos+1],bra.getA(pos+1),bra.getA(pos+1));
    bool dagger=true; // For the other term, contract with adjoint of ops
    ops.getOp(pos+1).contractR(theMgr.operR[pos],theMgr.operR[pos+1],
			       ket.getA(pos+1),bra.getA(pos+1),dagger);
    // TmpManagerResolvent* ptr=dynamic_cast<TmpManagerResolvent*>(&theMgr);
    if(resolvent){
      //ops.getOp(pos+1).contractR(ptr->normSR[pos],ptr->normSR[pos+1],
      //		       bra.getA(pos+1),bra.getA(pos+1),dagger);
	// cout<<"contractR2, for operRR("<<pos<<"), contracting the next "
// 	    <<"one ("<<theMgr.operRR[pos+1].GetDimensions()<<" with "
// 	    <<"ket "<<ket.getA(pos+1).getA().GetDimensions()<<" and "
// 	    <<"bra "<<bra.getA(pos+1).getA().GetDimensions()<<endl;
      contractR(theMgr.operRR[pos],theMgr.operRR[pos+1],
		ket.getA(pos+1),bra.getA(pos+1));
    }
    theMgr.storedR(pos);
  }
}

void Contractor::optimizeResolvent(const MPO& ops,const MPS& orig,
				   MPS& init,double eta,int D,double* errN,
				   double* errT){
  // I approximate by an MPS satisfying gauge condition, times a 
  // normalization factor. 
  
  vector<int> findims=ops.getOriginalDimensions();
  init.adjustPhysDimensions(findims);
  if(D>0)
    if(D>init.getBond())
      init.increaseBondDimension(D);
    else if(D<init.getBond()){
      cout<<"WARNING: Contractor::optimizeResolvent called with D("<<D<<") smaller "
	  <<"than in the starting MPS (D_init="<<init.getBond()<<"). Will truncate "
	  <<"the initial MPS!"<<endl;
      init=MPS(init,D);
    }
  init.gaugeCond('R',1); /* this first gauge ensures right Dr dimensions */
  init.gaugeCond('L',1); /* this ensures proper Dl and allows starting from left.*/
  //  if(D>0)
  //init.increaseBondDimension(D);

  double origNorm=real(contract(orig,orig));
  double constantN=origNorm/(orig.getNormFact()*orig.getNormFact()); // to be sure it works even without gauge

  int nrsites=orig.getLength();

  // Initialize temporary storage with gauge condition=false (trick 
  // because we use norm now for other quantity)
  //  TmpManagerResolvent tmpMgr(nrsites);
  bool resolvent=true;
  TmpManager tmpMgr(nrsites,0,resolvent);

  // Start optimization: sweep to right and left until convergence
  bool done=0;
  bool right=1;
  int pos=0;
  double distance=1000*orig.getNormFact()*orig.getNormFact(); // some crazily large initial value
  double olddistance=distance; // from iteration to iteration
  int round=1;

  double mytol=0.;

  mwArray M,M1,M2; // two terms to compose M
  mwArray N,N2;
  complex_t contr1,contr2;
    
  while(!done){
    /* Optimize matrix for position pos */
    cout<<"Contractor::optimizeResolvent pos "<<pos<<" round "<<round
	<<" distance "<<distance<<endl;
    //int d=init.getA(pos).getd();
    //cout<<"The local MPS is "<<dims<<endl;
    // Calculate M(pos) and N(pos) matrices 
    // 1. Obtain left and right terms for everything
    calculateL2(pos,tmpMgr,ops,orig,init,resolvent);
    calculateR2(pos,tmpMgr,ops,orig,init,resolvent);
    // 2. Contract the terms appropriately for M and N matrices
    bool dagger=true;
    ops.getOp(pos).contractMket(M1,tmpMgr.operL[pos],orig.getA(pos),
				tmpMgr.operR[pos],dagger);
    // extra M term now
    //cout<<"Calling contractMketId with A="<<orig.getA(pos).getA().GetDimensions()<<", L="<<tmpMgr.operRL[pos].GetDimensions()<<", R="<<tmpMgr.operRR[pos].GetDimensions()<<endl;
    contractMketId(M2,tmpMgr.operRL[pos],orig.getA(pos),tmpMgr.operRR[pos]);
    // And the N term!!
    ops.getOp(pos).contractN2(N2,tmpMgr.normL[pos],tmpMgr.normR[pos]);

    // Now compose the real ones
    N=N2+eta*eta*identityMatrix(N2.getDimension(0));
    M=M1+I_c*eta*M2;
// 	mwArray aux; // temporary, for -i eta M2
// 	m_times(1,aux,M2,mwArray(0.,eta));
// 	cout<<"M=M1+i eta M2; M1:"<<M1.GetDimensions()<<", M2:"
// 	    <<M2.GetDimensions()
// 	    <<", aux:"<<aux.GetDimensions()<<endl;
// 	m_add(1,M,M1,aux);  // M=M1+i*eta*M2
//	cout<<"M: "<<M.GetDimensions()<<endl;
// 	mwArray dimN=N.GetDimensions();int dim1=dimN(1);
// 	//cout<<"N has dimensions "<<dimN<<" of which the first one is "<<dim1<<endl;
// 	mwArray auxN; // temporary, to store the identity
// 	m_eye(1,auxN,mwArray(dim1));
// 	m_times(1,aux,auxN,mwArray(eta*eta)); // eta^2 Id
// 	cout<<"N=N1+eta^2 Id; N1:"<<N.GetDimensions()<<", Id:"
// 	    <<auxN.GetDimensions()
// 	    <<", aux:"<<aux.GetDimensions()<<endl;
// 	mwArray N2(N); // temporary copy
// 	m_add(1,N,N2,auxN); // N=N+eta^2 Id
//	cout<<"N: "<<N.GetDimensions()<<endl;

      if(pos==0&&round==1){
	if(M==mwArray(ZERO_c)){
	  cout<<"Warning: perturbing the initial MPS "<<endl;
	  init.perturb(1.);
	  for(int l=0;l<nrsites;l++){
	    tmpMgr.outofdateL(l);
	    tmpMgr.outofdateR(l);
	  }
	  continue;
	}
      }
      // Since my MPS contains a norm factor, which is not included by the contractions, I have to multiply it by hand
      N=(init.getNormFact()*init.getNormFact())*N;
      M=init.getNormFact()*M;
      // 3. Solve the linear equations
      mwArray newA;
      solveSite(N,M,newA,contr1,contr2,mytol); // Solve Matrix Equation (N not Id)
      
      // Computed with A before imposing gauge condition => the one solving the equation
      double newdist=constantN+real(contr1)-2*real(contr2);
      cout<<"Contractor::optimize pos "<<pos<<", round "<<round<<", distance"<<newdist<<endl;

      /** Now we have to substitute the solved A, but since I do not trust NA=M solution (numerical errors), if I get a distance 
	  larger than before, I do not change the state. */
      if(newdist<distance){
	init.setA(pos,newA);
	distance=newdist;
      }
      else{
	//cout<<"NOT Changing tensor at site "<<pos<<" as distance was "<<distance<<" and now is "<<newdist<<endl;
	newdist=distance;
      }
      // now impose gauge condition (according to the direction) 
      if(right){
	// the option 1 is required only in the edge, to keep norm factor, but
	// I multiply the neighbour, in case it does not change next time
	  init.gaugeCond(pos,'R',true); 
	  if(pos+2<nrsites)
	    tmpMgr.outofdateL(pos+2);     
      }
      else{
	  init.gaugeCond(pos,'L',true);
	  if(pos-2>=0)
	    tmpMgr.outofdateR(pos-2);     
      }
      // notify temporary storage as appropriate
      for(int l=pos+1;l<nrsites;l++)  // if gauge without changing next site
	tmpMgr.outofdateL(l);     
      for(int l=pos-1;l>=0;l--) // if gauge without changing next site
	tmpMgr.outofdateR(l);

      // Decide the convergence only at the end of sweep
      if(((right)&&(pos==nrsites-1))||
	 ((!right)&&(pos==0))){
	done=(abs(olddistance-distance)<max(convtol*abs(olddistance),NUMPREC))?1:0;
	if(!done){  // turn around and continue
	  if(right){
	    // init.gaugeCond(pos,'R',1); // this sets the norm factor but is already applied
 	    init.gaugeCond(pos,'L',true);pos--;
	    tmpMgr.outofdateR(pos);
	    if(pos-1>=0)
	      tmpMgr.outofdateR(pos-1);
	    if(pos+1<nrsites)
	      tmpMgr.outofdateL(pos+1);
	  }
	  else{
 	    // init.gaugeCond(pos,'L',1); // this sets the norm factor but is already applied
	    init.gaugeCond(pos,'R',true);pos++;
	    tmpMgr.outofdateL(pos);
	    if(pos+1<nrsites)
	      tmpMgr.outofdateL(pos+1);
	    if(pos-1>=0)
	      tmpMgr.outofdateR(pos-1);
	  }
	  olddistance=distance;
	  right=!right;
	  round++;
	}
      }
      else{ // not in the edge=> just advance in the same direction
	pos=right?pos+1:pos-1;
      }       
    } // end of while(!done)
  // When the state has converged, I have something which is not normalized.
  // In any case, if the original vector had a normFactor, it was ignored 
  // in the above calculation, hence it has to be added by hand now!
  //cout<<"At the end of optimize, changing norm "<<init.getNormFact();
  init.setNormFact(init.getNormFact()*orig.getNormFact());
  // If error has to be computed...
  if(errN!=0){
    double norm0=real(contract2(ops,init));
    double norm1=real(contract(orig,orig));
    complex_t crossR=contract(init,ops,orig);
    *errN=norm0+norm1-2*real(crossR);
    if(*errN<-1E-12){
      cout<<"WARNING: Contractor finds err^2<0!!!"<<endl;
    }
    if(errT!=0){
      *errT=norm0;
    }
  }
  tmpMgr.clear();
  cout<<"At the end of optimizeResolvent, "<<endl;
  double norm0=real(contract2(ops,init));
  double norm1=real(contract(orig,orig));
  complex_t crossR=contract(init,ops,orig);
  double differ=norm0+norm1-2*real(crossR);
  cout<<"   <Orig|Orig>="<<norm1;
  cout<<"   |Ops|Psi>|^2="<<norm0;
  cout<<"   ||Orig>-Ops|Psi>|^2="<<differ<<endl;
}



/** Auxiliary function for debugging: I should remove it completely */

void checkDistanceNow(double delta,double constantN,const MPS& init,const MPO& ops,const MPS& orig){
  Contractor& thec=Contractor::theContractor();
  double termN2=real(thec.contract2(ops,init));
  double term1=real(thec.contract(init,init));
  complex_t termM2=thec.contract(orig,init);
  complex_t termM1=thec.contract(orig,ops,init);
  complex_t termM3=thec.contract2(init,ops,orig);
  double origDist=(delta*delta/4)*termN2+term1+constantN-
    2*(real(termM2)+delta*imag(termM1)-(delta*delta/4)*real(termM3));
  cout<<"distance="<<origDist<<endl;
  if(1){
    cout<<"AND NOW CHECK CONTRACTIONS WITH FULL MPSs"<<endl;
    cout<<"  <init|init>="<<term1<<endl;
    cout<<", <init|H^2|init>="<<termN2<<endl;
    cout<<", <init|orig>="<<termM2<<endl;
    cout<<", <init|H|orig>="<<termM1<<endl;
    cout<<", <init|H^2|orig>="<<termM3<<endl;
    cout<<" => A+NA="<<(delta*delta/4)*termN2+term1<<endl;
    cout<<"    MA="<<delta*imag(termM1)-(delta*delta/4)*real(termM3)<<endl;
      }
}


double Contractor::getEntropy(const MPS& state,int pos){
  int length=state.getLength();
  if(pos==length) return 0.;
  if(pos<0||pos>length){
    cout<<"Error: trying to compute entropy for subchain size "<<pos
	<<" when length is "<<length<<endl;
    exit(212);
  }
  MPS auxMPS(state);
  // I copy the MPS and compute S of the normalized, with gauge cond
  // to the R one
  if(!auxMPS.isGaugeR()) auxMPS.gaugeCond('R',1);
  if(pos==0)
    pos=length/2;
  // Now left part is normalized => contract right part for
  //Lambda(L/2) MPS auxMPS2(*this); // trick to have an on basis for
  //right auxMPS2.gaugeCond('L',1);
  mwArray result(ONE_c);
  for(int k=length-1;k>=pos;k--){
    mwArray aux(result);
    contractR(result,aux,auxMPS.getA(k),auxMPS.getA(k));
  }
  vector<complex_t> Dval;mwArray U; //placeholder:ignored
  int D=auxMPS.getA(pos).getDl();
  wrapper::eig(reshape(result,Indices(D,D)),Dval,U);
  double entropy=0.;
  for(int k=0;k<D;k++){
    double tmp=real(Dval[k]);
    if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
  }
  return entropy;
}

mwArray Contractor::getRDM(const MPS& state,int pos){
  int length=state.getLength();
  if(pos<0||pos>=length){
    cout<<"Error: trying to compute RDM for site "<<pos
	<<" when length is "<<length<<endl;
    exit(212);
  }
  MPS auxMPS(state);
  //cout<<"Copied MPS "<<auxMPS<<endl;
  // I copy the MPS and compute rdm of the copy with gauge cond
  // to the R one (so contraction to the right is Id
  if(!auxMPS.isGaugeR()) auxMPS.gaugeCond('R');
  // Now left part is normalized => contract right part for
  //Lambda(L/2) MPS auxMPS2(*this); // trick to have an on basis for
  //right auxMPS2.gaugeCond('L',1);
  //cout<<"Imposed gaugeR "<<endl;
  mwArray result(ONE_c);
  for(int k=length-1;k>pos;k--){
    mwArray aux(result);
    contractR(result,aux,auxMPS.getA(k),auxMPS.getA(k));
    //cout<<"after contractR for k="<<k<<", results:"<<result.getDimensions()<<endl;
  }
  // Contract the site tensors
  mwArray site=auxMPS.getA(pos).getA(); // d Dl Dr
  //cout<<"Accessed site at pos "<<pos<<endl;
  int d=site.getDimension(0);
  int Dl=site.getDimension(1);
  int Dr=site.getDimension(2);
  site.reshape(Indices(d*Dl,Dr));
  result.reshape(Indices(Dr,Dr));result.transpose();
  result.multiplyLeft(site); // d*Dl,Dr
  result.reshape(Indices(d,Dl*Dr));
  site.reshape(Indices(d,Dl*Dr));site.Hconjugate();
  result.multiplyRight(site); // d x d
  double normF=auxMPS.getNormFact();
  if(abs(normF-1.)>1E-3){
    cout<<"WARNING: in getRDM() for site "<<pos<<" the original MPS was not normalized."
	<<endl;
  }
  result=normF*normF*result;
  return result;

  // // site.permute(Indices(2,1,3));site.reshape(Indices(Dl,d*Dr));
  // // //cout<<"site tensor reshaped: "<<site.getDimensions()<<endl;
  // // mwArray siteC(site);siteC.Hconjugate();
  // // //cout<<"siteC:"<<siteC.getDimensions()<<endl;
  // // siteC.multiplyRight(site); //d*Dr(bra) x d*Dr (ket)
  // // siteC.reshape(Indices(d,Dr,d,Dr));
  // // siteC.permute(Indices(1,3,2,4));
  // // siteC.reshape(Indices(d*d,Dr*Dr));
  // // result.reshape(Indices(-1,1));
  // // siteC.multiplyRight(result);
  // // siteC.reshape(Indices(d,d));
  // if necessary, add the normalization of the state
  // double normF=auxMPS.getNormFact();
  // if(abs(normF-1.)>1E-3){
  //   cout<<"WARNING: in getRDM() for site "<<pos<<" the original MPS was not normalized."
  // 	<<endl;
  // }
  // //  siteC=normF*normF*siteC;
  // //  return permute(siteC,Indices(2,1));
}

mwArray Contractor::getRDM(const MPS& state,int posL,int posR){
  int length=state.getLength();
  if(posL<0||posR>=length){
    cout<<"Error: trying to compute RDM for sites "<<posL<<"-"<<posR
	<<" when length is "<<length<<endl;
    exit(212);
  }
  MPS auxMPS(state);
  if(!auxMPS.isGaugeR()) auxMPS.gaugeCond('R');
  // Contraction from the left is now identity
  mwArray result(ONE_c);
  for(int k=length-1;k>posR;k--){
    mwArray aux(result);
    contractR(result,aux,auxMPS.getA(k),auxMPS.getA(k));
  }
  // Notice that, as prepared by contractR, the first index is that of
  // the bra, so we need to transpose the result now
  result.permute(Indices(2,1));
  // And contraction from the right is result (Dr x Dr)  
  mwArray A=auxMPS.getA(posL).getA();
  int d=A.getDimension(0);int Dl0=A.getDimension(1);int Dr=A.getDimension(2);
  A.permute(Indices(2,1,3));A.reshape(Indices(Dl0*d,Dr));
  int dphys=d; // temporary dimensions of the vector
  int Dv=Dr;
  for(int k=posL+1;k<=posR;k++){
    mwArray aux=auxMPS.getA(k).getA();
    //cout<<"Multiplying tensor for site "<<k<<", A="<<aux.getDimensions()<<" with tmp result "<<A.getDimensions()<<endl;
    Indices dims=aux.getDimensions();
    int d=dims[0];int Dl=dims[1];int Dr=dims[2];
    if(Dl!=Dv){
      cout<<"Error: Dimensions do not agree in getRDM!"<<endl;
      exit(1);
    }
    aux.permute(Indices(2,1,3));
    aux.reshape(Indices(Dl,d*Dr));
    A.multiplyRight(aux);
    dphys=dphys*d;
    Dv=Dr;
    A.reshape(Indices(Dl0*dphys,Dv));
  }
  // if the mps has a normalization factor, we need to include it!!
  A=auxMPS.getNormFact()*A;
  // And at the end
  result.multiplyLeft(A); // Dl0*dphys x Dr
  result.reshape(Indices(Dl0,dphys,Dv));
  result.permute(Indices(2,1,3)); // dphys x Dl0 x Dv
  result.reshape(Indices(dphys,Dl0*Dv));
  A.reshape(Indices(Dl0,dphys,Dv));
  A.permute(Indices(2,1,3));A.reshape(Indices(dphys,Dl0*Dv));
  result.multiplyRight(Hconjugate(A)); // dphys x dphys
  return result;
}  


double Contractor::getEntropy(const MPS& state,int posL,int posR){
  int length=state.getLength();
  if(posL<0||posR<0||posL>=length||posR>=length||posL>posR){
    cout<<"Error: trying to compute entropy between sites "<<posL
	<<" and "<<posR<<" when length is "<<length<<endl;
    exit(212);
  }
  // simpler case: one of them the edge
  if(posL==0)
    return getEntropy(state,posR+1);
  if(posR==length-1)
    return getEntropy(state,posL);
  MPS auxMPS(state);
  // I copy the MPS and compute S of the normalized, with gauge cond
  // to the R one
  // TODO: impose gauge to simplify longest part to trace? 
  if(!auxMPS.isGaugeR()) auxMPS.gaugeCond('R',1);
  // Now left part is normalized => contract right part until posR
  for(int k=length-1;k>=posR+1;k--){
    auxMPS.gaugeCond(k,'L',true);
  }
  // Now the right part is also nomalized
  // Now contract the part of the state with the bra, to obtain the
  // D^2 x D^2 overlap matrix
  mwArray overl;
  for(int k=posL;k<=posR;k++){
    mwArray tmp(overl);
    contractLmiddle(overl,tmp,auxMPS.getA(k),auxMPS.getA(k));
    //cout<<"After contracting in site "<<k<<", tmp is "<<overl<<endl;
  } 
  // Reshape as due
  int Dl=auxMPS.getA(posL).getDl();
  int Dr=auxMPS.getA(posR).getDr();
  //cout<<"Reshaping the overlap matrix with Dl="<<Dl<<", Dr="<<Dr<<endl;
  overl.reshape(Indices(Dl,Dl,Dr,Dr));
  overl.permute(Indices(1,3,2,4));
  overl.reshape(Indices(Dl*Dr,Dl*Dr));
  // And perform the non orthogonal diagonalization
  vector<complex_t> Dval;mwArray U; //placeholder:ignored
  wrapper::eig(overl,Dval,U); // only eigenvalues
  double entropy=0.;
  for(int k=0;k<Dl*Dr;k++){
    double tmp=real(Dval[k]);
    if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
  }
  return entropy;
}

void Contractor::getSchmidtValues(const MPS& state,vector<complex_t>& lambdas,
				  int pos){
  int length=state.getLength();
  lambdas.clear();
  if(pos==length){lambdas.push_back(ONE_c); return;}
  if(pos<0||pos>length){
    cout<<"Error: trying to compute entropy for subchain size "<<pos
  	<<" when length is "<<length<<endl;
    exit(212);
  }
  MPS auxMPS(state);
  // I copy the MPS and compute S of the normalized, with gauge cond
  // to the R one
  if(!auxMPS.isGaugeR()) auxMPS.gaugeCond('R',1);
  if(pos==0)
    pos=length/2;
  // Now left part is normalized => contract right part for
  //Lambda(L/2) MPS auxMPS2(*this); // trick to have an on basis for
  //right auxMPS2.gaugeCond('L',1);
  mwArray result(ONE_c);
  for(int k=length-1;k>=pos;k--){
    mwArray aux(result);
    contractR(result,aux,auxMPS.getA(k),auxMPS.getA(k));
  }
  mwArray U; //placeholder:ignored
  int D=auxMPS.getA(pos).getDl();
  wrapper::eig(reshape(result,Indices(D,D)),lambdas,U);
}

void Contractor::calculateL0(int pos,TmpManager& theMgr,
			    const MPS& ket,const MPS& bra,
			    bool gauge){
  //cout<<"Contractor::calculateL(pos="<<pos<<", tmpMgr="<<&theMgr<<")"<<endl;
  if(!theMgr.isuptodateL(pos)){ //else, nothing to do, just return!
    calculateL0(pos-1,theMgr,ket,bra,gauge); // take the one before
    //cout<<"Using tmpL of pos "<<pos-1<<" uptodate"<<endl;
    contractL(theMgr.operL[pos],theMgr.operL[pos-1],ket.getA(pos-1),bra.getA(pos-1));
    if(!gauge)
      contractNormL(theMgr.normL[pos],theMgr.normL[pos-1],bra.getA(pos-1));
    theMgr.storedL(pos);    
  }
}

void Contractor::calculateR0(int pos,TmpManager& theMgr,
			    const MPS& ket,const MPS& bra,
			    bool gauge){
  //cout<<"Contractor::calculateR(pos="<<pos<<", tmpMgr="<<&theMgr<<")"<<endl;
  if(!theMgr.isuptodateR(pos)){ //else, nothing to do, just return!
    calculateR0(pos+1,theMgr,ket,bra,gauge); // take next one and use it!
    //cout<<"Using tmpR of pos "<<pos+1<<" uptodate"<<endl;
    contractR(theMgr.operR[pos],theMgr.operR[pos+1],ket.getA(pos+1),bra.getA(pos+1));
    if(!gauge)
      contractNormR(theMgr.normR[pos],theMgr.normR[pos+1],bra.getA(pos+1));
    theMgr.storedR(pos);
  }
}


/** Auxiliary debugging function */
double Contractor::distanceVectors(const vector<const MPS*>& kets,
				   const vector<complex_t>& beta,
				   const MPS& init){
  //cout<<"Distance from init to sum of vectors=";
  int nrterms=kets.size();  
  complex_t sum={0.,0.};
  double sumNorms=0.; //the constant term!
  for(int nr_=0;nr_<nrterms;nr_++){
    sum=sum+beta[nr_]*contract(*kets[nr_],init);
    complex_t betaK=beta[nr_];
    double modBeta=betaK.re*betaK.re+betaK.im+betaK.im;
    sumNorms+=modBeta*real(contract(*kets[nr_],*kets[nr_]));
  }
  return real(contract(init,init))-2*real(sum)+sumNorms;
}


void Contractor::optimizeMPS(const MPS& ket,MPS& init,int D,double* err){
  vector<const MPS*> kets(1,&ket);
  vector<complex_t> beta(1,ONE_c);
  optimizeSum(kets,beta,init,D);
  if(err!=0){ // After return, estimate error, i required
    complex_t norm0=contract(ket,ket);
    complex_t norm1=contract(init,init);
    complex_t overl=contract(init,ket);
    *err=real(norm0)+real(norm1)-2*real(overl);
  }
}


void  Contractor::optimizeSum(vector<const MPS*> kets,vector<complex_t> beta,
			   MPS& init,int D){
  // check length of vectors
  int nrterms=kets.size(); 
  if(nrterms!=beta.size()){
    cout<<"Error: incompatible length of arguments passed to optimize"<<endl;
    exit(1);
  }
  int nrsites=kets[0]->getLength();
  for(int nr_=1;nr_<nrterms;nr_++)
    if(nrsites!=kets[nr_]->getLength()){
      cout<<"ERROR: Not all MPS in optimize have the same number of sites"
	  <<endl;
      exit(1);
    }
  //if(D>0){
  //init.increaseBondDimension(D);
  //}
  if(D>0){
    if(D>init.getBond())
      init.increaseBondDimension(D);
    else if(D<init.getBond()){
      cout<<"WARNING: Contractor::optimizeSum called with D("<<D<<") smaller "
	  <<"than in the starting MPS (D_init="<<init.getBond()<<"). Will truncate "
	  <<"the initial MPS!"<<endl;
      init=MPS(init,D);
    }
  }
  init.gaugeCond('R',1); // this first gauge ensures right Dr dimensions
  //cout<<"*Applied ritht gauge condition"<<endl;
  init.gaugeCond('L',1); /* this ensures proper Dl and allows starting 
			    from left.*/
  //cout<<"*Applied left gauge condition"<<endl;
  // Initialize temporary storage with gauge condition=true
  bool gauge=true;
  vector<TmpManager*> tmpMgr;
  for(int nr_=0;nr_<nrterms;nr_++)
    tmpMgr.push_back(new TmpManager(nrsites,gauge));

  // Start optimization: sweep to right and left until convergence
  bool done=0;
  bool right=1;
  int pos=0;
  double distance=1000; // initially large value
  distance=distanceVectors(kets,beta,init);
  double olddistance=distance;
  //cout<<"Distance before starting: "<<distance<<", init.normFactor="
  //<<init.getNormFact()<<endl;
  int round=1;

  mwArray Mloc;
  mwArray M; //mwArray N;
  mwArray contr1,contr2,newA;
    
  while(!done){
    /* Optimize matrix for position pos */
    //cout<<"Contractor::optimizeSum for pos "<<pos<<", round "<<round
    //	<<", distance(stored)="<<distance<<", and computed now="
    //	<<distanceVectors(kets,beta,init)<<endl;
    // Calculate M(pos) matrices (N not needed)
    for(int nr_=0;nr_<nrterms;nr_++){
      // 1. Obtain left and right terms for everything
      calculateL0(pos,*tmpMgr.at(nr_),*kets[nr_],init,gauge);
      calculateR0(pos,*tmpMgr.at(nr_),*kets[nr_],init,gauge);
      // 2. Contract the terms appropriately for M matrix
      // cout<<"After calculateL->"<<tmpMgr.operL[pos]<<", calculateR->"
      // <<tmpMgr.operR[pos]<<", with op->"<<ops.getOp(pos)<<endl;
      contractMketId(Mloc,tmpMgr.at(nr_)->operL[pos],kets[nr_]->getA(pos),
		     tmpMgr.at(nr_)->operR[pos]);
      if(nr_==0) // M not initialized
	M=kets[nr_]->getNormFact()*beta[nr_]*Mloc;
      else
	M=M+kets[nr_]->getNormFact()*beta[nr_]*Mloc;
    }
    if(pos==0&&round==1){
      if(M.isNull()){
	init.perturb(1.);
	for(int l=0;l<nrsites;l++){
	  for(int nr_=0;nr_<nrterms;nr_++){
	    tmpMgr.at(nr_)->outofdateL(l);
	    tmpMgr.at(nr_)->outofdateR(l);
	  }
	}
	continue;
      }
    }
    // 3. Solve the linear equations => since I have gauge condition, I only set A=M/normfact
    solveSiteGauge(newA,contr1,contr2,M,init.getNormFact());
    newA=(norm(M)/init.getNormFact())*newA;
    
    init.setA(pos,newA); 
    // double auxDist=distanceVectors(kets,beta,init);
    //    cout<<"After setting A distance shoud be "<<newdist
    // 	  <<" and it is (brute force)"<<auxDist<<endl;
    if(right){
      if(pos<nrsites-1)
	init.gaugeCond(pos,'R');
      else
	init.gaugeCond(pos,'R',true); // at the end, keep norm factor
    }
    else{
      if(pos>0)
	init.gaugeCond(pos,'L');
      else
	init.gaugeCond(pos,'L',true);// at the end, keep norm factor
    }
    // notify temporary storage as appropriate
    for(int nr_=0;nr_<nrterms;nr_++){
      for(int l=pos+1;l<nrsites;l++) tmpMgr.at(nr_)->outofdateL(l);
      for(int l=pos-1;l>=0;l--) tmpMgr.at(nr_)->outofdateR(l);
    }
    // This has changed contr1 AND contr2!!!
    mwArray A=init.getA(pos).getA();
    A.reshape(Indices(A.getNrComponents(),1));// vector form
    contr1=Hconjugate(A)*A;
    contr2=Hconjugate(A)*M; // Calculates <A|M>        
    double normM=real((Hconjugate(M)*M).getElement(0));
   // Decide the convergence only at the end of sweep
    if(((right)&&(pos==nrsites-1))||
       ((!right)&&(pos==0))){
      done=checkDistanceConvergence(distance,contr1,contr2,
				    init.getNormFact(),sqrt(normM),round);
      //done=(abs(olddistance-distance)<max(convtol*abs(olddistance),NUMPREC))?1:0;
      if(!done){  // turn around and continue
	if(right){
	  init.gaugeCond(pos,'L');pos--;
	}
	else{
	  init.gaugeCond(pos,'R');pos++;
	}
	olddistance=distance;
	right=!right;
	round++;
      }
    }
    else{ // not in the edge=> just advance in the same direction
	pos=right?pos+1:pos-1;
      }    
    } // while(!done)
  for(int nr_=0;nr_<nrterms;nr_++){
    tmpMgr.at(nr_)->clear();
    delete tmpMgr[nr_];
  }
  tmpMgr.clear();
  //cout<<"Distance going out: "<<distance<<endl;
}


void Contractor::contractLmiddle(mwArray& result,const mwArray& tmpL,
			   const Site& ket,const Site& bra){
  int d=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  //  cout<<"Contractor::contractLmiddle with tmp "<<tmpL.getDimensions()
  //  <<", ket "<<ket.getDimensions()<<", bra "<<bra.getDimensions()<<endl;
  if(tmpL.isEmpty()){ // first contraction, only MPS
    result=bra.getA();
    result.reshape(Indices(d,Dl*Dr));
    result.conjugate();
    mwArray aux=ket.getA();
    aux.reshape(Indices(dp,Dlp*Drp));
    aux.permute(Indices(2,1));
    result.multiplyLeft(aux);
    result.reshape(Indices(Dlp,Drp,Dl,Dr));
    result.permute(Indices(3,1,4,2));
    result.reshape(Indices(Dl*Dlp,Dr,Drp));
    //cout<<"First contraction: result "<<result<<endl;
  }
  else{
    int chi2=tmpL.getDimension(0);
    int betab=tmpL.getDimension(1);int betak=tmpL.getDimension(2);
#ifdef CHECKDIMS
    if((Dl!=betab)||(Dlp!=betak)||(d!=dp)){
      cout<<"Error in Contractor::contractLmiddle dimensions: bra="
	  <<bra.getDimensions()
	  <<", ket="<<ket.getDimensions()<<", tmpL="
	  <<tmpL.getDimensions()<<endl;
      exit(212);
    }
#endif
    result=permute(ket.getA(),Indices(2,1,3));
    result.reshape(Indices(Dlp,dp*Drp));
    mwArray aux0(tmpL); aux0.reshape(Indices(chi2*betab,betak));
    result.multiplyLeft(aux0);
    //cout<<"After multiplying in ket, "<<result.getDimensions()<<endl;
    result.reshape(Indices(chi2,betab*dp,Drp));
    result.permute(Indices(2,3,1));
    result.reshape(Indices(betab*dp,Drp*chi2));
    mwArray aux=bra.getA();
    aux.permute(Indices(3,2,1),true);
    aux.reshape(Indices(Dr,Dl*d));
    result.multiplyLeft(aux);
    //cout<<"After multiplying in bra, "<<result.getDimensions()<<endl;
    result.reshape(Indices(Dr,Drp,chi2));
    result.permute(Indices(3,1,2));
  }
}

void Contractor::contractOperLmiddle(mwArray& result,const mwArray& tmpL,
				     const Site& ket,const Operator& op,const Site& bra){
  // I can cheat the operator by giving a modified bra that holds the extra leg
  //Indices dimsTmp=result.getDimensions(); // dimL,Dlbra,xi,Dlket
  //cout<<"Here tmpL="<<tmpL.getDimensions()<<", bra="<<bra.getDimensions()<<" ket="<<ket.getDimensions()<<", oper:"<<op.getDimensions()<<endl;
  int k=bra.getPos();
  int dimL=result.getDimension(0);
  int Dlb=result.getDimension(1);int xi=result.getDimension(2);
  int Dlk=result.getDimension(3);
  int Drb=bra.getDr();int d=bra.getd();
  //cout<<"k="<<k<<" dimL="<<dimL<<" Dlb="<<Dlb<<" Dlk="<<Dlk<<" xi="<<xi<<endl;
  mwArray auxBra(bra.getA());//
  result.reshape(Indices(dimL*Dlb,xi,Dlk));
  auxBra.reshape(Indices(d*Dlb*Drb,1));
  auxBra.multiplyRight(reshape(identityMatrix(dimL),Indices(1,dimL*dimL)));
  auxBra.reshape(Indices(d,Dlb,Drb,dimL,dimL));
  auxBra.permute(Indices(1,4,2,5,3));
  auxBra.reshape(Indices(d,dimL*Dlb,dimL*Drb));
  Site auxBra_(k,auxBra);
  mwArray tmp(result);
  op.contractL(result,tmp,ket,auxBra_);
  //cout<<"After contractL, I get resultL:"<<result<<" but I cannot read its dimensions!"<<endl;
  // now result will be thick Dbr, xir, Drk
  Indices dimsTmp(result.getDimensions());
  //cout<<"result.dimsTmp="<<dimsTmp<<endl;
  //cout<<result.getDimensions()<<endl;
  // dimsTmp=result.getDimensions(); // dimL*Drb,xir,Drk
  result.reshape(Indices(dimL,Drb,dimsTmp[1],dimsTmp[2]));
  //cout<<"After reshape, "<<result.getDimensions()<<endl;
  //result=result_;
}

void Contractor::sweepPartNextExcitedState(int pos0,int lim,const MPO& ops,int D,
					   const vector<MPS*>& computedLevels,
					   MPS& init,double offset,int knr){
  vector<double> energies(computedLevels.size(),0);
  int nextK=computedLevels.size();
  for(int k=0;k<nextK;k++){
    //cout<<"Contractor::findNextExcitedState "<<nextK<<" computing energy of "
    //	<<"computed state nr "<<k<<", "<<*computedLevels[k]<<endl;
    energies[k]=real(contract(*computedLevels[k],ops,*computedLevels[k]))+offset; 
  }
  sweepPartWithProjectorPenalty(pos0,lim,ops,D,energies,computedLevels,
				init,offset,knr);
}

void Contractor::sweepPartWithProjectorPenalty(int pos0,int lim,const MPO& ops,int D,
					       const vector<double>& energies,
					       const vector<MPS*>& computedLevels,
					       MPS& init,double offset,int knr){
  int nrsites=ops.getLength();  
  vector<int> findims=ops.getDimensions();
  init.adjustPhysDimensions(findims);
  bool right=lim>=pos0;
  const char dir=right?'R':'L';

  if(D>0)
    if(D>init.getBond())
      init.increaseBondDimension(D);
    else if(D<init.getBond()){
      cout<<"WARNING: Contractor::sweepPartWithProjPenalty called with D("<<D<<") smaller "
	  <<"than in the starting MPS (D_init="<<init.getBond()<<"). Will truncate "
	  <<"the initial MPS!"<<endl;
      init=MPS(init,D);
    }
  // I will assume that dimensions are ok from the beginning
  int pos=pos0;
  for(int k=0;k<pos;k++)
    init.gaugeCond(k,'R',true);
  for(int k=nrsites-1;k>pos;k--)
    init.gaugeCond(k,'L',true);

  int nextK=computedLevels.size();
  double energy=real(contract(init,ops,init));
  for(int k=0;k<nextK;k++){
    complex_t overlapK=contract(init,*computedLevels[k]);
    energy+=energies[k]*real(overlapK*conjugate(overlapK));
  }
  // The norm term in the TmpManager will hold the contraction with
  // GS, so I have to "trick" it into thinking gauge=false;
  TmpManagerExc tmpMgr(nrsites,nextK);
  clock_t start=clock();clock_t finish=clock();
  double val=energy;
  bool done=0;
  while(!done){
    // // finish=clock();
    // // cout<<"SweepPartWithProjectorPenalty["<<nextK<<"], pos "<<pos<<", energy "<<val; // energy;
    // // cout<<" time "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
    // // start=clock();
    // 1. Obtain left and right terms for everything. In this case,
    // this includes the partial contractions of the result MPS with
    // the MPO, and also its contraction with the GS, to implement the
    // projector term.
    calculateLexck(pos,tmpMgr,ops,init,computedLevels);
    calculateRexck(pos,tmpMgr,ops,init,computedLevels);
    // Place for the results of the diagonalization
    vector<complex_t> diagVals;
    mwArray U(init.getA(pos).getA()); 
    // todo: include gauge from previous!
    U.reshape(Indices(-1,1));//U.fillWithOne();
    // 2. contract the Heff matrix, but with the projection out of
    // each of the computed eigenstates, which here amounts to
    // substracting Ek times |MPS_k><MPS_k| for each of them, so we
    // need to compute the effective MPS_k vector (on the site, after contracting
    // the rest of the chain, and that's our local projector.
    vector<mwArray> effectiveMPS(nextK,mwArray::empty);
    for(int l=0;l<nextK;l++){
      //      cout<<"Computing effective vector for level "<<l<<endl;
      contractMbraId(effectiveMPS[l],tmpMgr.normL[l][pos],computedLevels[l]->getA(pos),
		     tmpMgr.normR[l][pos]);
    }
    mwArray tmpOp=ops.getOp(pos).getFullData();
    TensorMultiplierProjMulti multi(tmpMgr.operL[pos],tmpOp,tmpMgr.operR[pos],energies,
				    effectiveMPS,true,offset*ONE_c);  // forcing Hermiticity
    //				    effectiveMPS,false); // not forcing Hermiticity
    solveEigs(multi,min(knr,multi.getSize()),"SR",diagVals,U,true); // Compute one eigenvalue-vector
    val=real(*(complex_t*)&diagVals[0])-offset; // first eigenvalue
    U.resize(Indices(U.getDimension(0),1)); // first column only!
    init.setA(pos,U); 
    bool passOnX=false; // whether to multiply next tensor by the transformation
    if((right&&pos==lim-1)||(!right&&pos==lim+1)||(pos==lim))
      passOnX=true;
    //    cout<<"About to apply gauge "<<dir<<" to "<<pos<<" and pass TErm="<<passOnX<<endl;
    init.gaugeCond(pos,dir,passOnX);
    // notify temporary storage as appropriate
    for(int l=pos+1;l<nrsites;l++)
      tmpMgr.outofdateL(l);
    for(int l=pos-1;l>=0;l--)
      tmpMgr.outofdateR(l);
    if(right) pos++;
    else pos--;
    // Check if sweep is done
    done=right?pos>=lim:pos<=lim;
    energy=val;
  } // while(!done)


}

void Contractor::findNextExcitedStateLM(const MPO& ops,int D,
				      const vector<MPS*>& computedLevels,
				      double* lambdak,
				      MPS& init,double offset,int knr,
				      const string& tmpfile,int maxTime){
  findNextExcitedState(ops,D,computedLevels,lambdak,init,offset,knr,tmpfile,maxTime,true);
}

void Contractor::findNextExcitedState(const MPO& ops,int D,
				      const vector<MPS*>& computedLevels,
				      double* lambdak,
				      MPS& init,double offset,int knr,
				      const string& tmpfile,int maxTime,bool useLM){
  vector<double> energies(computedLevels.size(),0);
  int nextK=computedLevels.size();
  // I can estimate the offset I need to guarantee (more or less) that the next state is the
  // lowest of the new H
  // double minEk(0.),maxEk(0.);
  double minEk(10*ops.getLength()),maxEk(-10*ops.getLength());
  for(int k=0;k<nextK;k++){
    //cout<<"Contractor::findNextExcitedState "<<nextK<<" computing energy of "
    //<<"computed state nr "<<k<<", "<<contract(*computedLevels[k],ops,*computedLevels[k])<<endl;
    energies[k]=(real(contract(*computedLevels[k],ops,*computedLevels[k]))+offset);
    if(k==0||energies[k]<minEk) minEk=energies[k];
    if(k==0||energies[k]>maxEk) maxEk=energies[k];
  }
  // Now I decide the penalty for the computed levels
  double newOffset=max(20*(maxEk-minEk),abs(offset));
  if(useLM){ // in the special case that I am using the largest magnitude eigensolver, can set those energies to 0
    // TODO!!!! + offset???? YES!!: it is later substracted in TensorMultiplier
    newOffset=0.-offset;
  }
  else{
    // A special case is if I only had one level (or almost degenerate ones), then I have to guess, so I add L/2 (big)
    if(nextK==1||newOffset<1E-6) newOffset=computedLevels[0]->getLength()*.5;
  }
  //  cout<<"Setting new penalties for findNextExcitedState: "<<newOffset<<endl;//" (old was "<<offset<<")"<<endl;
  for(int k=0;k<nextK;k++){
    energies[k]=newOffset;
  }
  bool done=0;int cnt=0;
  while(!done){
    findGroundStateWithProjectorPenalty(ops,D,energies,computedLevels,
				      lambdak,init,offset,knr,tmpfile,maxTime,useLM);
    done=1;cnt++;
    if(excorthtol>0){ // only if it is set, this check is run
      //If the result is not orthogonal enough to existing levels, the offset was probably not enough
      // check orthogonality
      for(int k=0;k<nextK;k++){
	complex_t overl=contract(init,*computedLevels[k]);
	if(abs(overl)>excorthtol){
	  cout<<"Contractor::findNextExcitedState found state not orthogonal to existing level "
	      <<k<<", with overlap "<<overl<<". Reapeating with increased penalty ("<<energies[k]
	      <<"->"<<2*energies[k]<<")"<<endl;
	  done=0;energies[k]=2*energies[k];
	}
      }      
    }
  }  
}

void Contractor::findGroundStateWithProjectorPenaltyLM(const MPO& ops,int D,
						     const vector<double>& energies,
						     const vector<MPS*>& computedLevels,
						     double* lambdak,
						     MPS& init,double offset,int knr,
						     const string& tmpfile,
						     int maxTime){
  findGroundStateWithProjectorPenalty(ops,D,energies,computedLevels,lambdak,init,offset,knr,tmpfile,maxTime,true);
}
  
void Contractor::findGroundStateWithProjectorPenalty(const MPO& ops,int D,
						     const vector<double>& energies,
						     const vector<MPS*>& computedLevels,
						     double* lambdak,
						     MPS& init,double offset,int knr,
						     const string& tmpfile,
						     int maxTime,bool useLM){
  // cout<<"findGroundStateWithProjectorPenalty called with D="<<D<<", init.D= "
  //     <<init.getBond()<<", offset="<<offset<<" and penalties (-)"<<energies<<endl;
#ifndef USING_PRIMME
  cout<<"Compiled without PRIMME support => eigensolver = ARPACK"<<endl;
  solver=arpack;
#endif
  int nrsites=ops.getLength();  
  vector<int> findims=ops.getDimensions();
  init.adjustPhysDimensions(findims);
  // For saving intermediate results, if tmpfile and maxTime are given
  clock_t lastSaved=clock();
  clock_t now=clock();
  if(D>0)
    if(D>init.getBond())
      init.increaseBondDimension(D);
    else if(D<init.getBond()){
      cout<<"WARNING: Contractor::findGroundStateWithProjPenalty called with D("<<D<<") smaller "
	  <<"than in the starting MPS (D_init="<<init.getBond()<<"). Will truncate "
	  <<"the initial MPS!"<<endl;
      init=MPS(init,D);
    }
  init.gaugeCond('R',1); // this first gauge ensures right Dr dimensions
  init.gaugeCond('L',1); // this ensures proper Dl and allows starting 
			    // from left.
  bool done=0;int round=0;
  bool right=1; int pos=0; //nrsites-1;
  *lambdak=0.; // initialize
  if(energies.size()!=computedLevels.size()){
    cout<<"ERROR: findGroundStateWithProjectorPenalty received different number of projectors and penalties"<<endl;
    exit(1);
  }
  int nextK=computedLevels.size();
  double energy=real(contract(init,ops,init))+offset;
  //cout<<"exp. value of H in init = "<<energy-offset<<endl;
  // TensorMultiplierProj *substracts* sth, so I need to pass minus the penalties
  vector<double> energies_;
  for(int k=0;k<nextK;k++){
    complex_t overlapK=contract(init,*computedLevels[k]);
    //cout<<"overlap of init with level "<<k<<" = "<<overlapK<<endl;
    energy+=energies[k]*real(overlapK*conjugate(overlapK));
    energies_.push_back(-energies[k]);
  }
  energy-=offset;
  // The norm term in the TmpManager will hold the contraction with
  // GS, so I have to "trick" it into thinking gauge=false;
  bool gauge=false;bool ketN=false; 
  TmpManagerExc tmpMgr(nrsites,nextK);
  clock_t start=clock();clock_t finish=clock();
  double val=energy;
  while(!done){
    //if(pos==0){
    //finish=clock();
     //cout<<"Round number "<<round<<", pos "<<pos<<", energy "<<val<<" time "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
     //start=clock();
    // }
    //if(round>2*nrsites) convtol=1E-3;
    // 1. Obtain left and right terms for everything. In this case,
    // this includes the partial contractions of the result MPS with
    // the MPO, and also its contraction with the GS, to implement the
    // projector term.
    calculateLexck(pos,tmpMgr,ops,init,computedLevels);
    calculateRexck(pos,tmpMgr,ops,init,computedLevels);
    // Place for the results of the diagonalization
    vector<complex_t> diagVals;
    //mwArray U;
    mwArray U(init.getA(pos).getA()); 
    // todo: include gauge from previous!
    U.reshape(Indices(-1,1));//U.fillWithOne();
    // 2. contract the Heff matrix, but with the projection out of
    // each of the computed eigenstates, which here amounts to
    // substracting Ek times |MPS_k><MPS_k| for each of them, so we
    // need to compute the effective MPS_k vector (on the site, after contracting
    // the rest of the chain, and that's our local projector.
    vector<mwArray> effectiveMPS(nextK,mwArray::empty);
    for(int l=0;l<nextK;l++){
      //cout<<"Computing effective vector @pos "<<pos<<" for level "<<l<<endl;
      contractMbraId(effectiveMPS[l],tmpMgr.normL[l][pos],computedLevels[l]->getA(pos),
		     tmpMgr.normR[l][pos]);
    }
    mwArray tmpOp=ops.getOp(pos).getFullData();
    TensorMultiplierProjMulti multi(tmpMgr.operL[pos],tmpOp,tmpMgr.operR[pos],energies_,
				    effectiveMPS,true,offset*ONE_c);  // forcing Hermiticity
    //				    effectiveMPS,false); // not forcing Hermiticity
    const char* which="SR";
    if(useLM) which="LM";
    solveEigs(multi,min(knr,multi.getSize()),which,diagVals,U,true); // Compute one eigenvalue-vector
    val=real(*(complex_t*)&diagVals[0])-offset; // first eigenvalue
    U.resize(Indices(U.getDimension(0),1)); // first column only!
    init.setA(pos,U); 
    //if(val-energy>abs(energy)*1E-10){
    if((val-energy>1E-8)&&(val-energy>abs(energy+offset)*1E-10)){
      if((val-energy>abs(energy+offset)*1E-6)){
	cout<<"Error: energy increasing in this step!! oldenergy="<<energy
	    <<", new one="<<val<<", difference="<<(val-energy)
	    <<", relative increase="<<(val-energy)/(energy+offset)<<endl;
	// ofstream output("badM.m");
	// mwArray M; // the effective Hamiltonian 
	// ops.getOp(pos).contractN(M,tmpMgr.operL[pos],tmpMgr.operR[pos]);
	// M=.5*(M+Hconjugate(M));
	// putForMatlab(output,M);
	// cout<<endl;output.close();
	exit(212);
      }
      else{
	cout<<"WARNING: energy increasing in this step!! oldenergy="<<energy
	    <<", new one="<<val<<", difference="<<(val-energy)
	    <<", relative increase="<<(val-energy)/(energy+offset)
	    <<" I CONTINUE!!"<<endl;
      }
    }
    if(right){
      init.gaugeCond(pos,'R',true);
    }
    else{
      init.gaugeCond(pos,'L',true);
    }
    // notify temporary storage as appropriate
    for(int l=pos+1;l<nrsites;l++)
      tmpMgr.outofdateL(l);
    for(int l=pos-1;l>=0;l--)
      tmpMgr.outofdateR(l);
  
    if(tmpfile!=""){ // Saving intermediate results (only at the end)
      bool doneRound=(right&&pos==nrsites-1)||((!right)&&pos==0);
      bool saveNow=false;
      if(doneRound){
	now=clock();
	if(maxTime>0){
	  saveNow=((now-lastSaved)*1./CLOCKS_PER_SEC>=maxTime);
	}
	else saveNow=true;
	if(saveNow){
	  init.exportMPS(tmpfile.data());
	  //cout<<"Saved Temporary resuts to "<<tmpfile<<endl;
	  lastSaved=clock();
	}
      }
    }
  
    //if (pos==1) done=1;
    // Now move one pos and check convergence if we are at the end f
    // one sweep
    if(right){ //moving right
      if(pos==nrsites-1){ //last reached: check convergence and change sense
	right=0;
	done=abs(energy-val)<convtol*abs(energy+offset);
	//done=abs(energy-val)<convtol*abs(energy);
	if(!done){
	  //	  cout<<"Done test failed with relative change "
	  //  <<abs((energy-val)/(energy+offset))<<" including offset "<<offset<<endl;
	  // Check if the value is too close to zero!
	  if(abs(energy-val)<1E-12){
	    cout<<"Assuming convergence because difference is small!"<<endl;
	    done=true;
	  }
	  if(round>MAX_ROUND_EXCTD){
	    cout<<"Assuming convergence because too many rounds!!!"<<endl;
	    done=true;	    
	  }
	  // if(abs(val)<1E-10){
	  //   cout<<"WARNING***** Assuming convergence because E is small!"<<endl;
	  //   done=true;
	  // }
	}
	energy=val;
	if(!done){
	  init.gaugeCond(pos,'L',true);pos--;
	}
	round++;
      }
      else{
	pos++;
      }
    }
    else{// moving left
      if(pos==0){ //first reached: check convergence and change sense
	right=1;
	done=abs(energy-val)<convtol*abs(energy+offset);
	//done=abs(energy-val)<convtol*abs(energy);
	if(!done){
	  //	  cout<<"Done test failed with relative change "
	  //  <<abs((energy-val)/(energy+offset))<<" including offset "<<offset<<endl;
	  // Check if the value is too close to zero!
	  if(abs(energy-val)<1E-12){
	    cout<<"Assuming convergence because difference is small!"<<endl;
	    done=true;
	  }
	  if(round>MAX_ROUND_EXCTD){
	    cout<<"Assuming convergence because too many rounds!!!"<<endl;
	    done=true;	    
	  }
	  // if(abs(val)<1E-10){
	  //   cout<<"WARNING***** Assuming convergence because E is small!"<<endl;
	  //   done=true;
	  // }
	  
	}
	energy=val;
	if(!done){
	  init.gaugeCond(pos,'R',true);pos++;
	}
	round++;
      }
      else{
	pos--;
      }
    }
  }
  *lambdak=energy;
  // cout<<"Leaving Contractor::findNextExcitedState after "<<round<<" rounds"<<endl;
}

void Contractor::calculateLexck(int pos,TmpManagerExc& theMgr,const MPO& ops,
				const MPS& ket,const vector<MPS*>& braLevels){
  //cout<<"Contractor::calculateLexck(pos="<<pos<<", tmpMgr="<<&theMgr<<")"<<endl;
  if(!theMgr.isuptodateL(pos)){ //else, nothing to do, just return!
    calculateLexck(pos-1,theMgr,ops,ket,braLevels); // take the one before
    //cout<<"Using tmpL of pos "<<pos-1<<" uptodate"<<endl;
    contractOperL(theMgr.operL[pos],theMgr.operL[pos-1],ket.getA(pos-1),
		  ops.getOp(pos-1),ket.getA(pos-1));
    int nextK=theMgr.getNextK();
    for(int k=0;k<nextK;k++)
      contractL(theMgr.normL[k][pos],theMgr.normL[k][pos-1],ket.getA(pos-1),
		braLevels[k]->getA(pos-1));
    theMgr.storedL(pos);    
  }
}

void Contractor::calculateRexck(int pos,TmpManagerExc& theMgr,const MPO& ops,
				const MPS& ket,const vector<MPS*>& braLevels){
  //cout<<"Contractor::calculateRexck(pos="<<pos<<", tmpMgr="<<&theMgr<<")"<<endl;
  if(!theMgr.isuptodateR(pos)){ //else, nothing to do, just return!
    calculateRexck(pos+1,theMgr,ops,ket,braLevels); // take the one before
    //cout<<"Using tmpR of pos "<<pos+1<<" uptodate"<<endl;
    contractOperR(theMgr.operR[pos],theMgr.operR[pos+1],ket.getA(pos+1),
 		  ops.getOp(pos+1),ket.getA(pos+1));
    int nextK=theMgr.getNextK();
    for(int k=0;k<nextK;k++)
      contractR(theMgr.normR[k][pos],theMgr.normR[k][pos+1],ket.getA(pos+1),
		braLevels[k]->getA(pos+1));
    theMgr.storedR(pos);    
  }
}

  // For the GEVP, oprL(R) will hold the effective Hamiltonian (opsA) and normL(R) the "mass" term (opsB)
void Contractor::calculateLgev(int pos,TmpManager& theMgr,const MPO& opsA,const MPO& opsB,
			       const MPS& ket,const MPS& bra){
  //cout<<"Contractor::calculateLgev(pos="<<pos<<", tmpMgr="<<&theMgr<<")"<<endl;
  if(!theMgr.isuptodateL(pos)){ //else, nothing to do, just return!
    calculateLgev(pos-1,theMgr,opsA,opsB,ket,bra); // take the one before
    //cout<<"Using tmpL of pos "<<pos-1<<" uptodate"<<endl;
    contractOperL(theMgr.operL[pos],theMgr.operL[pos-1],ket.getA(pos-1),
 		  opsA.getOp(pos-1),bra.getA(pos-1));
     //contractOperL(theMgr.getOperL(pos),theMgr.getOperL(pos-1),ket.getA(pos-1),
     //	  ops.getOp(pos-1),bra.getA(pos-1));
    contractOperL(theMgr.normL[pos],theMgr.normL[pos-1],ket.getA(pos-1),
		  opsB.getOp(pos-1),bra.getA(pos-1));
    theMgr.storedL(pos);    
  }
}

void Contractor::calculateRgev(int pos,TmpManager& theMgr,const MPO& opsA,const MPO& opsB,
			       const MPS& ket,const MPS& bra){
  //cout<<"Contractor::calculateRgev(pos="<<pos<<", tmpMgr="<<&theMgr<<")"<<endl;
  if(!theMgr.isuptodateR(pos)){ //else, nothing to do, just return!
    calculateRgev(pos+1,theMgr,opsA,opsB,ket,bra); // take next one and use it!
    //cout<<"Using tmpR of pos "<<pos+1<<" uptodate, with ket "
    //	<<endl;
     //<<ket.getA(pos+1).getA()<<" and bra "<<bra.getA(pos+1).getA()<<endl;
    contractOperR(theMgr.operR[pos],theMgr.operR[pos+1],ket.getA(pos+1),
		  opsA.getOp(pos+1),bra.getA(pos+1));
    contractOperR(theMgr.normR[pos],theMgr.normR[pos+1],ket.getA(pos+1),
		  opsB.getOp(pos+1),bra.getA(pos+1));
    //cout<<"Computed tmpMgr.operR["<<pos<<"]"<<theMgr.operR[pos]<<endl;
    theMgr.storedR(pos);
  }
}

void Contractor::calculateLgevexck(int pos,TmpManagerExc& theMgr,const MPO& opsA,const MPO& opsB,
				const MPS& ket,const vector<MPS*>& braLevels){
  //cout<<"Contractor::calculateLexck(pos="<<pos<<", tmpMgr="<<&theMgr<<")"<<endl;
  if(!theMgr.isuptodateL(pos)){ //else, nothing to do, just return!
    calculateLgevexck(pos-1,theMgr,opsA,opsB,ket,braLevels); // take the one before
    //cout<<"Using tmpL of pos "<<pos-1<<" uptodate"<<endl;
    contractOperL(theMgr.operL[pos],theMgr.operL[pos-1],ket.getA(pos-1),
		  opsA.getOp(pos-1),ket.getA(pos-1));
    contractOperL(theMgr.normL0[pos],theMgr.normL0[pos-1],ket.getA(pos-1),
		  opsB.getOp(pos-1),ket.getA(pos-1));
    int nextK=theMgr.getNextK();
    for(int k=0;k<nextK;k++)
      contractL(theMgr.normL[k][pos],theMgr.normL[k][pos-1],ket.getA(pos-1),
		braLevels[k]->getA(pos-1));
    theMgr.storedL(pos);    
  }
}

void Contractor::calculateRgevexck(int pos,TmpManagerExc& theMgr,const MPO& opsA,const MPO& opsB,
				   const MPS& ket,const vector<MPS*>& braLevels){
  //cout<<"Contractor::calculateRexck(pos="<<pos<<", tmpMgr="<<&theMgr<<")"<<endl;
  if(!theMgr.isuptodateR(pos)){ //else, nothing to do, just return!
    calculateRgevexck(pos+1,theMgr,opsA,opsB,ket,braLevels); // take the one before
    //cout<<"Using tmpR of pos "<<pos+1<<" uptodate"<<endl;
    contractOperR(theMgr.operR[pos],theMgr.operR[pos+1],ket.getA(pos+1),
 		  opsA.getOp(pos+1),ket.getA(pos+1));
    contractOperR(theMgr.normR0[pos],theMgr.normR0[pos+1],ket.getA(pos+1),
 		  opsB.getOp(pos+1),ket.getA(pos+1));
    int nextK=theMgr.getNextK();
    for(int k=0;k<nextK;k++)
      contractR(theMgr.normR[k][pos],theMgr.normR[k][pos+1],ket.getA(pos+1),
		braLevels[k]->getA(pos+1));
    theMgr.storedR(pos);    
  }
}


#ifdef USING_PRIMME
extern "C" {
#include "primme.h"
}
#endif

void Contractor::solveEigs(Multiplier& multi,const int k,const char* which,
			   std::vector<complex_t>& D,mwArray& U,bool eigv,double tol){
  //  cout<<"solveEigs with solver ("<<solver<<")"<<endl;
  if(tol==0) tol=eigtol;
#ifdef USING_PRIMME
  primme_preset_method method;
  primme_target target; // which eigenvalues (largest or smaller?)
  int numTargetShifts=1; // only for LM
  double targetShifts[1]={0.};
  if(!strcmp("SM",which)||(!strcmp("SR",which))){ // strcmp gives 0 if they are equal
    target=primme_smallest;
  }
  else if(!strcmp("LM",which)){
    target=primme_largest_abs;
  }
  else{
    cout<<"WARNING: Targeting eigenvalues "<<which<<" not (yet) supported for PRIMME=> going back to default (SM)"
	<<endl;
    target=primme_smallest;
  }
#endif
  switch(solver){
  case fulleig:
    { // solving the full matrix!
      mwArray fullM;
      multi.getFullTensor(fullM);
      // rescale by largest element (? is this worse?)
      complex_t largest;Indices posL;fullM.getMaxElement(largest,posL);
      fullM=(1./abs(largest))*fullM;
      //cout<<"Calling full eig solver now"<<endl;
      wrapper::eig(fullM,D,U,eigv);
      // scale back the eigenvalues
      for(int i=0;i<D.size();i++) D[i]=D[i]*abs(largest);
      return;
    }
    break;
  case arpack:
    //cout<<"Calling arpack solver now"<<endl;
    wrapper::eigs(multi,k,which,D,U,eigv,tol);
    return;
    break;
#ifdef USING_PRIMME
  case primme:
    method=DYNAMIC_ver;
    break;
  case primme_JDQR:
    method=JDQR_ver;
    break;
  case primme_arnoldi:
    method=Arnoldi_ver;
    break;
#endif
  default:
    cout<<"Error: unknown type of eigensolver ("<<solver<<") in Contractor"<<endl;
    exit(1);
  }
#ifdef USING_PRIMME
  // WARNING: Primme not good with very small matrices, it seems. So, if the full dimension is below 
  // LOWEST_PRIMME, I switch to full solution!
  if(multi.getSize()<LOWEST_PRIMME){
    //    cout<<"WARNING: Switching from eigs to full solution for small matrix!"<<endl;
    mwArray fullM;
    multi.getFullTensor(fullM);
    fullM=.5*(fullM+Hconjugate(fullM));
    //cout<<"Calling full eig solver now"<<endl;
    wrapper::eig(fullM,D,U,eigv);
  }
  else
    //cout<<"Calling now eigs_primme(k="<<k<<",target="<<target<<",D="<<D<<",U="<<U<<",tol="<<tol<<")"<<endl;
    wrapper::eigs_primme(multi,k,target,D,U,tol,method,numTargetShifts,targetShifts);
#endif
}
 
void Contractor::solveEigsClosest(Multiplier& multi,const int k,const char* which,double targetVal,
				  std::vector<complex_t>& D,mwArray& U,double tol){
  if(solver==arpack){
    cout<<"Error: solveEigsClosest can only be called for PRIMME solver. For ARPACK, "
	<<"call solveEigs, which=SM and set an appropriate offset=-target."<<endl;
    exit(1);
  }
  if(solver==fulleig){
    cout<<"Error: calling solveEigsClosest with full diagonalization. Select values from the full spectrum, instead."<<endl;
    exit(1);
  }
#ifdef USING_PRIMME
  primme_preset_method method;
  primme_target target; 
  switch(solver){
  case primme:
    method=DYNAMIC_ver;
    break;
  case primme_arnoldi:
    method=Arnoldi_ver;
    break;
  case primme_JDQR:
    method=JDQR_ver;
    break;
  default:
    cout<<"Error, solver="<<solver<<" unknown or not supported in selveEigsClosest"<<endl;
    exit(1);
  }
  if(!strcmp("LEQ",which)||!strcmp("leq",which))
    target=primme_closest_leq;
  else if(!strcmp("GEQ",which)||!strcmp("geq",which))
    target=primme_closest_geq;
  else if(!strcmp("ABS",which)||!strcmp("abs",which))
    target=primme_closest_abs;
  else{
    cout<<"Unknown mode "<<which<<" for solveEigsClosest"<<endl;
    exit(1);
  }
  int numShift=1;
  wrapper::eigs_primme(multi,k,target,D,U,tol,method,numShift,&targetVal);
#endif
}

void Contractor::solveEigsGen(Multiplier& multiA,Multiplier& multiB,const int k,const char* which,
			      std::vector<complex_t>& D,mwArray& U,bool eigv,double tol){
  if(tol==0) tol=eigtol;
#ifdef USING_PRIMME
  primme_preset_method method;
  primme_target target; // which eigenvalues (largest or smaller?)
  if(!strcmp("SM",which)||(!strcmp("SR",which))){ // strcmp gives 0 if they are equal
    target=primme_smallest;
  }
  else if(!strcmp("LM",which)){
    target=primme_largest;
  }
  else{
    cout<<"WARNING: Targeting eigenvalues "<<which<<" not (yet) supported for PRIMME=> going back to default (SM)"
	<<endl;
    target=primme_smallest;
  }
#endif
  switch(solver){
  case fulleig:
    { // solving the full matrix!
      // mwArray fullM;
      // multi.getFullTensor(fullM);
      // wrapper::eig(fullM,D,U,eigv);
      cout<<"Not supported eig for generalized case"<<endl;
      exit(1);
    }
    break;
  case arpack:
    //    wrapper::eigs(multi,k,which,D,U,eigv,tol);
    cout<<"Not supported eigs for generalized case"<<endl;
    exit(1);
    break;
#ifdef USING_PRIMME
  case primme:
    method=DYNAMIC_ver;
    break;
  case primme_JDQR:
    method=JDQR_ver;
    break;
  case primme_arnoldi:
    method=Arnoldi_ver;
    break;
#endif
  default:
    cout<<"Error: unknown type of eigensolver ("<<solver<<") in Contractor"<<endl;
    exit(1);
  }
#ifdef USING_PRIMME
  wrapper::eigsgen_primme(multiA,multiB,k,target,D,U,tol,method);
#endif
}

#include "PowerMultiplier.h"
#include "LinearCombinationMultiplier.h"

// Version 3: Generalized EV (H-e)^2x=(H-E) x global
// ops should be H, and be Hermitian. I construct (H-e)^2
void Contractor::findClosestEigenstate(const MPO& ops,int D,double* lambda, //const char* which,
				       MPS& init,double target,int knr){
#ifndef USING_PRIMME
  cout<<" WARNING: Compiled without PRIMME support => eigensolver = ARPACK => only support for targetting"
      <<" interior eigenstates is SM (closer to zero) => call instead SM with an offset"<<endl;
  solver=arpack;
#endif
  if((solver==arpack)||(solver==fulleig)){
    cout<<"ERROR: findClosestEigenstate only supported for solvers in the PRIMME package"<<endl;
    exit(1);
  }
  double offset=-target;
  int nrsites=ops.getLength();  
  vector<int> findims=ops.getDimensions();
  init.adjustPhysDimensions(findims);
  //  if(D>0)
  // init.increaseBondDimension(D);
  if(D>0)
    if(D>init.getBond())
      init.increaseBondDimension(D);
    else if(D<init.getBond()){
      cout<<"WARNING: Contractor::findClosestEigenstate called with D("<<D<<") smaller "
	  <<"than in the starting MPS (D_init="<<init.getBond()<<"). Will truncate "
	  <<"the initial MPS!"<<endl;
      init=MPS(init,D);
    }
  init.gaugeCond('R',1); // this first gauge ensures right Dr dimensions
  init.gaugeCond('L',1); /* this ensures proper Dl and allows starting 
			    from left.*/
  
  // An auxiliary MPO for the square
  MPO ops2(ops.getLength());
  const MPO* auxPtr[2]={&ops,&ops};
  MPO::join(2,auxPtr,ops2);

  bool done=0;int round=0;
  bool gauge=false; // Need to keep two terms: H^2 and H (which plays the role of N term)
  bool right=1; int pos=0; //nrsites-1;
  *lambda=0.; // initialize
  bool ketN=false; 
  TmpManager tmpMgr(nrsites,gauge);
  double energy;
  {double auxV=real(contract(init,ops,init));
    energy=real(contract(init,ops2,init)-2*target*auxV*ONE_c+target*target*ONE_c)/(auxV-target);} // initial value 
  clock_t start=clock();clock_t finish=clock();
  double val=energy;
  double valL; //from last iteration
  while(!done){
    finish=clock();
    cout<<"Round number "<<round<<", pos "<<pos<<", energy "<<val; //energy;
    cout<<" time (previous round)"<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
    if(round>20) convtol=1E-3;
    start=clock();
    // 1. Obtain left and right terms for everything
    calculateLgev(pos,tmpMgr,ops2,ops,init,init);
    calculateRgev(pos,tmpMgr,ops2,ops,init,init);
 
    // Place for the results of the diagonalization
    vector<complex_t> diagVals;
    // mwArray U;
    mwArray U(init.getA(pos).getA()); 
    // todo: include gauge from previous!
    U.reshape(Indices(-1,1));//U.fillWithOne();
    mwArray tmpOp2=ops2.getOp(pos).getFullData();
    mwArray tmpOp=ops.getOp(pos).getFullData();
    TensorMultiplierHermitian multi2(tmpMgr.operL[pos],tmpOp2,tmpMgr.operR[pos]);    
    TensorMultiplierHermitianOffset multi(tmpMgr.normL[pos],tmpOp,tmpMgr.normR[pos],offset);    
    int knr_=min(knr,multi.getSize());
    // Now the appropriate LC of multipliers for (H-E)^2
    vector<complex_t> coeffs(2);coeffs[0]=ONE_c;coeffs[1]=-2*target*ONE_c;
    vector<Multiplier*> multis(2);multis[0]=&multi2;multis[1]=&multi;
    LinearCombinationMultiplier HE2(coeffs,multis,-ONE_c*target*target);
    //    solveEigs(multi2,1,"SR",diagVals,U,true);
    solveEigsGen(HE2,multi,knr,"SR",diagVals,U,true);
    //solveEigs(HE2,knr,"SR",diagVals,U,true);

    //    val=real(*(complex_t*)&diagVals[0])-offset; // first eigenvalue
    U.resize(Indices(U.getDimension(0),1)); // first column only!
    init.setA(pos,U); 
    // TODO: More efficient calculation here
    {
      mwArray num(U.getDimensions()),den(U.getDimensions());HE2.product(U,num);multi.product(U,den);
      valL=val;
      //      val=real((Hconjugate(U)*num).getElement(0))/real((Hconjugate(U)*den).getElement(0)); 
      //double auxV=real(contract(init,ops,init));
      //val=real(contract(init,ops2,init)-2*target*auxV*ONE_c+target*target*ONE_c)/(auxV-target);
      val=real((Hconjugate(U)*den).getElement(0));
      //cout<<"After solveEigenvalueEq, lambda="<<val<<", and first ev "<<(*(complex_t*)&diagVals[0])
      //<<" and norm(U)="<<Hconjugate(U)*U<<endl;
    }
    //val=real(*(complex_t*)&diagVals[0]); // less accurate?
    //if(val-energy>abs(energy)*1E-10){
    if((abs(val)-abs(valL)>1E-8)){  //&&(abs(val)-abs(energy)>abs(energy)*1E-10)){
      // if((val-energy>abs(energy+offset)*1E-6)){
      // 	cout<<"Error: energy increasing in this step!! oldenergy="<<energy
      // 	    <<", new one="<<val<<", difference="<<(val-energy)
      // 	    <<", relative increase="<<(val-energy)/(energy+offset)<<endl;
      // 	exit(212);
      // }
      //else{
      cout<<"WARNING: \"eigenvalue\" of H-E increasing in absolute value in this step!! old one="<<valL
	  <<", new one="<<val<<", difference="<<(abs(val)-abs(valL))
	  <<", relative increase (incl offset)="<<(abs(val)-abs(valL))/(target)
	  <<" I CONTINUE"<<endl; //", but skip this change!!"<<endl;
      //      val=valL;
    }
      //else
    init.setA(pos,U); 
      //}
    //init.setA(pos,U); 
    if(right){
      //if(pos<nrsites-1)
      init.gaugeCond(pos,'R',true);
    }
    else{
      //if(pos>0)
      init.gaugeCond(pos,'L',true);
    }
    // notify temporary storage as appropriate
    for(int l=pos+1;l<nrsites;l++)
      tmpMgr.outofdateL(l);
    for(int l=pos-1;l>=0;l--)
      tmpMgr.outofdateR(l);
    
    //if (pos==1) done=1;
    // Now move one pos and check convergence if we are at the end f
    // one sweep
    if(right){ //moving right
      if(pos==nrsites-1){ //last reached: check convergence and change sense
	right=0;
	done=abs(abs(energy)-abs(val))<convtol*abs(target);
	if(!done){
	  //	  cout<<"Done test failed (R) ("<<energy<<"->"<<val<<") with relative change "
	  //  <<abs((abs(energy)-abs(val))/(target))<<" wrt offset "<<target<<", with val="<<val<<endl;
	  // Check if the value is too close to zero!
	  if(abs(abs(energy)-abs(val))<1E-12){
	    cout<<"Assuming convergence because difference is small!"<<endl;
	    done=true;
	  }
	  if(abs(val)<1E-10){
	    cout<<"WARNING***** Assuming convergence because E is small!"<<endl;
	    done=true;
	  }
	}
	energy=val;
	if(!done){
	  init.gaugeCond(pos,'L',true);pos--;
	}
	round++;
      }
      else{
	pos++;
      }
    }
    else{// moving left
      if(pos==0){ //first reached: check convergence and change sense
	right=1;
	done=abs(abs(energy)-abs(val))<convtol*abs(target);
	if(!done){
	  //	  cout<<"Done test failed (L) ("<<energy<<"->"<<val<<") with relative change "
	  //  <<abs((abs(energy)-abs(val))/(target))<<" wrt offset "<<offset<<endl;
	  // Check if the value is too close to zero!
	  if(abs(abs(energy)-abs(val))<1E-12){
	    cout<<"Assuming convergence because difference is small!"<<endl;
	    done=true;
	  }
	}
	energy=val;
	if(!done){
	  init.gaugeCond(pos,'R',true);pos++;
	}
	round++;
      }
      else{
	pos--;
      }
    }
  }
  *lambda=energy;
}

void Contractor::findNextClosestEigenstate(const MPO& ops,int D,
					   const vector<MPS*>& computedLevels,
					   double* lambdak,
					   MPS& init,double target,int knr){
  if(target==0) // back to findNextExcitedState
    return findNextExcitedState(ops,D,computedLevels,lambdak,init,0.,knr);
#ifndef USING_PRIMME
  cout<<"Compiled without PRIMME support => eigensolver = ARPACK"<<endl;
  solver=arpack;
#endif
  int nrsites=ops.getLength();  
  vector<int> findims=ops.getDimensions();
  init.adjustPhysDimensions(findims);
  if(D>0)
    if(D>init.getBond())
      init.increaseBondDimension(D);
    else if(D<init.getBond()){
      cout<<"WARNING: Contractor::findNextClosestEigenstate called with D("<<D<<") smaller "
	  <<"than in the starting MPS (D_init="<<init.getBond()<<"). Will truncate "
	  <<"the initial MPS!"<<endl;
      init=MPS(init,D);
    }
  init.gaugeCond('R',1); // this first gauge ensures right Dr dimensions
  init.gaugeCond('L',1); // this ensures proper Dl and allows starting 
			    // from left.

  // An auxiliary MPO for the square
  MPO ops2(nrsites);
  const MPO* auxPtr[2]={&ops,&ops};
  MPO::join(2,auxPtr,ops2);

  bool done=0;int round=0;
  bool right=1; int pos=0; //nrsites-1;
  *lambdak=0.; // initialize
  vector<double> energies(computedLevels.size(),0);
  vector<double> energies2(computedLevels.size(),0);
  int nextK=computedLevels.size();
  double energyDen=real(contract(init,ops,init)); //<H>
  double energyNum=real(contract(init,ops2,init))-2*target*energyDen;//<H^2>-2E<H>
  for(int k=0;k<nextK;k++){
    //cout<<"Contractor::findNextExcitedState "<<nextK<<" computing energy of "
    //	<<"computed state nr "<<k<<", "<<*computedLevels[k]<<endl;
    energies[k]=real(contract(*computedLevels[k],ops,*computedLevels[k])); 
    energies2[k]=-(2*target-energies[k])*energies[k]; // weights for the H op in numerator op
    //cout<<"Result: "<<energies[k]<<endl;
    complex_t overlapK=contract(init,*computedLevels[k]);
    energyNum+=(2*target-energies[k])*energies[k]*real(overlapK*conjugate(overlapK));
    energyDen-=energies[k]*real(overlapK*conjugate(overlapK));
  }
  double energy=(energyNum+target*target)/energyDen;
  // The norm term in the TmpManager will hold the contraction with
  // GS, so I have to "trick" it into thinking gauge=false;
  bool gev=true; 
  TmpManagerExc tmpMgr(nrsites,nextK,gev);
  clock_t start=clock();clock_t finish=clock();
  double val=energy;
  double valL;
  while(!done){
    finish=clock();
    cout<<"Round number "<<round<<", pos "<<pos<<", energy "<<val; // energy;
    cout<<" time "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
    if(round>2*nrsites) convtol=1E-3;
    start=clock();
    // 1. Obtain left and right terms for everything. In this case,
    // this includes the partial contractions of the result MPS with
    // the MPO, and also its contraction with the GS, to implement the
    // projector term.
    calculateLgevexck(pos,tmpMgr,ops2,ops,init,computedLevels);
    calculateRgevexck(pos,tmpMgr,ops2,ops,init,computedLevels);
    // Place for the results of the diagonalization
    vector<complex_t> diagVals;
    //mwArray U;
    mwArray U(init.getA(pos).getA()); 
    // todo: include gauge from previous!
    U.reshape(Indices(-1,1));//U.fillWithOne();
    // 2. contract the Heff matrix, but with the projection out of
    // each of the computed eigenstates, which here amounts to
    // substracting Ek times |MPS_k><MPS_k| for each of them, so we
    // need to compute the effective MPS_k vector (on the site, after contracting
    // the rest of the chain, and that's our local projector.
    vector<mwArray> effectiveMPS(nextK,mwArray::empty);
    for(int l=0;l<nextK;l++){
      //      cout<<"Computing effective vector for level "<<l<<endl;
      contractMbraId(effectiveMPS[l],tmpMgr.normL[l][pos],computedLevels[l]->getA(pos),
		     tmpMgr.normR[l][pos]);
    }
    mwArray tmpOp=ops.getOp(pos).getFullData();
    mwArray tmpOp2=ops2.getOp(pos).getFullData();
    // The mass matrix term (<H-E-sum E_i>)
    TensorMultiplierProjMulti multiDen(tmpMgr.normL0[pos],tmpOp,tmpMgr.normR0[pos],energies,
				       effectiveMPS,true,-target*ONE_c); //Hermitian
    int knr_=min(knr,multiDen.getSize());
    // The numerator (H^2-2EH+sum(2E-E_i)E_i +E^2)
    TensorMultiplierHermitian multi2(tmpMgr.operL[pos],tmpOp2,tmpMgr.operR[pos]);    
    TensorMultiplierHermitian multi1(tmpMgr.normL0[pos],tmpOp,tmpMgr.normR0[pos]);    
    // Now the appropriate LC of multipliers for (H-E)^2
    vector<complex_t> coeffs(2);coeffs[0]=ONE_c;coeffs[1]=-2*target*ONE_c;
    vector<Multiplier*> multis(2);multis[0]=&multi2;multis[1]=&multi1;
    LinearCombinationMultiplier multiNum1(coeffs,multis); // H^2 -2 E H
    TensorMultiplierProjMulti multiNum(multiNum1,energies2,effectiveMPS,ONE_c*target*target);
    solveEigsGen(multiNum,multiDen,knr,"SR",diagVals,U,true);
    effectiveMPS.clear(); multis.clear();coeffs.clear();
    //solveEigs(multi,min(knr,multi.getSize()),"SR",diagVals,U,true); // Compute one eigenvalue-vector
    //cout<<"After solveEigenvalueEq, lambda="<<val<<", and first ev "<<(*(complex_t*)&diagVals[0])<<endl;
    U.resize(Indices(U.getDimension(0),1)); // first column only!

    //val=real(*(complex_t*)&diagVals[0]); // first eigenvalue
    {
      mwArray num;multiNum.product(U,num);
      //mwArray den;multiDen.product(U,den);
      valL=val; // the last one
      val=real((Hconjugate(U)*num).getElement(0));///real((Hconjugate(U)*den).getElement(0));
    }

    // The quotient computed here is now val
    //if(val-energy>abs(energy)*1E-10){
    if((abs(val)-abs(valL)>1E-8)){ //&&(val-energy>abs(energy)*1E-10)){
      // if((val-energy>abs(energy)*1E-6)){
      // 	cout<<"Error: energy increasing in this step!! oldenergy="<<energy
      // 	    <<", new one="<<val<<", difference="<<(val-energy)
      // 	    <<", relative increase="<<(val-energy)/(energy)<<endl;
	// ofstream output("badM.m");
	// mwArray M; // the effective Hamiltonian 
	// ops.getOp(pos).contractN(M,tmpMgr.operL[pos],tmpMgr.operR[pos]);
	// M=.5*(M+Hconjugate(M));
	// putForMatlab(output,M);
	// cout<<endl;output.close();
      // 	exit(212);
      // }
      // else{
      cout<<"WARNING: \"eigenvalue\" increasing in this step!! oldenergy="<<valL
	  <<", new one="<<val<<", difference="<<abs(val)-abs(valL)
	  <<", relative increase="<<(abs(val)-abs(valL))/(target)
	  <<" I CONTINUE!!"<<endl;
      //energy=val;
    }
    // else
    //   init.setA(pos,U); 
    init.setA(pos,U); 

    if(right){
      //if(pos<nrsites-1)
      init.gaugeCond(pos,'R',true);
    }
    else{
      //if(pos>0)
      init.gaugeCond(pos,'L',true);
    }
    // notify temporary storage as appropriate
    for(int l=pos+1;l<nrsites;l++)
      tmpMgr.outofdateL(l);
    for(int l=pos-1;l>=0;l--)
      tmpMgr.outofdateR(l);
    
    //if (pos==1) done=1;
    // Now move one pos and check convergence if we are at the end f
    // one sweep
    if(right){ //moving right
      if(pos==nrsites-1){ //last reached: check convergence and change sense
	right=0;
	done=abs(abs(energy)-abs(val))<convtol*abs(energy+target);
	if(!done){
	  //	  cout<<"Done test failed with relative change "
	  //  <<abs((abs(energy)-abs(val))/(energy+target))<<endl;
	  // Check if the value is too close to zero!
	  if(abs(energy-val)<1E-12){
	    cout<<"Assuming convergence because difference is small!"<<endl;
	    done=true;
	  }
	  // if(abs(val)<1E-10){
	  //   cout<<"WARNING***** Assuming convergence because E is small!"<<endl;
	  //   done=true;
	  // }
	}
	energy=val;
	if(!done){
	  init.gaugeCond(pos,'L',true);pos--;
	}
	round++;
      }
      else{
	pos++;
      }
    }
    else{// moving left
      if(pos==0){ //first reached: check convergence and change sense
	right=1;
	done=abs(energy-val)<convtol*abs(energy);
	if(!done){
	  //	  cout<<"Done test failed with relative change "
	  //  <<abs((energy-val)/(energy))<<endl;
	  // Check if the value is too close to zero!
	  if(abs(energy-val)<1E-12){
	    cout<<"Assuming convergence because difference is small!"<<endl;
	    done=true;
	  }
	  // if(abs(val)<1E-10){
	  //   cout<<"WARNING***** Assuming convergence because E is small!"<<endl;
	  //   done=true;
	  // }
	  
	}
	energy=val;
	if(!done){
	  init.gaugeCond(pos,'R',true);pos++;
	}
	round++;
      }
      else{
	pos--;
      }
    }
  }
  *lambdak=energy;

}


void Contractor::sweepPartWithProjPenalty(int pos0,int lim,const MPO& ops,int D,
					  double mu,MPS& projected,
					  MPS& init,double offset,int knr){
#ifndef USING_PRIMME
  cout<<"Compiled without PRIMME support => eigensolver = ARPACK"<<endl;
  solver=arpack;
#endif
  int nrsites=ops.getLength();  
  vector<int> findims=ops.getDimensions();
  init.adjustPhysDimensions(findims);
  bool right=lim>=pos0;
  const char dir=right?'R':'L';
  if(D>0)
    if(D>init.getBond())
      init.increaseBondDimension(D);
    else if(D<init.getBond()){
      cout<<"WARNING: Contractor::sweepPartWithProjPenalty called with D("<<D<<") smaller "
	  <<"than in the starting MPS (D_init="<<init.getBond()<<"). Will truncate "
	  <<"the initial MPS!"<<endl;
      init=MPS(init,D);
    }
  // I will assume that dimensions are ok from the beginning
  int pos=pos0;
  for(int k=0;k<pos;k++)
    init.gaugeCond(k,'R',true);
  for(int k=nrsites-1;k>pos;k--)
    init.gaugeCond(k,'L',true);

  vector<MPS*> computedLevels(1,&projected);
  vector<double> energies(1,-mu);
  // For saving intermediate results, if tmpfile and maxTime are given
  clock_t lastSaved=clock();
  clock_t now=clock();

  bool done=0;
  double energy=real(contract(init,ops,init));
  complex_t overlapK=contract(init,projected);
  energy+=mu*real(overlapK*conjugate(overlapK));
  // The norm term in the TmpManager will hold the contraction with
  // GS, so I have to "trick" it into thinking gauge=false;
  bool gauge=false;bool ketN=false; 
  TmpManagerExc tmpMgr(nrsites,1);
  clock_t start=clock();clock_t finish=clock();
  double val=energy;
  while(!done){
    //finish=clock();
    //    cout<<"Sweep part with penalty, pos "<<pos<<", energy "<<val; // energy;
    //cout<<" time "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
    //start=clock();
    // 1. Obtain left and right terms for everything. In this case,
    // this includes the partial contractions of the result MPS with
    // the MPO, and also its contraction with the GS, to implement the
    // projector term.
    calculateLexck(pos,tmpMgr,ops,init,computedLevels);
    calculateRexck(pos,tmpMgr,ops,init,computedLevels);
    // Place for the results of the diagonalization
    vector<complex_t> diagVals;
    mwArray U(init.getA(pos).getA()); 
    // todo: include gauge from previous!
    U.reshape(Indices(-1,1));//U.fillWithOne();
    // 2. contract the Heff matrix, but with the contribution of
    // the projector, which here amounts to adding mu times
    // |MPS_k><MPS_k|, so we need to compute the effective MPS_k 
    // vector (on the site, after contracting
    // the rest of the chain, and that's our local projector.
    vector<mwArray> effectiveMPS(1,mwArray::empty);
    contractMbraId(effectiveMPS[0],tmpMgr.normL[0][pos],
		   computedLevels[0]->getA(pos),
		   tmpMgr.normR[0][pos]);
    mwArray tmpOp=ops.getOp(pos).getFullData();
    TensorMultiplierProjMulti multi(tmpMgr.operL[pos],tmpOp,
				    tmpMgr.operR[pos],energies,
				    effectiveMPS,true,offset*ONE_c);  // forcing Hermiticity
    //				    effectiveMPS,false); // not forcing Hermiticity
    solveEigs(multi,min(knr,multi.getSize()),"SR",diagVals,U,true); // Compute one eigenvalue-vector
    val=real(*(complex_t*)&diagVals[0])-offset; // first eigenvalue
    U.resize(Indices(U.getDimension(0),1)); // first column only!
    // Now I project out the effectiveLevels and normalize
    init.setA(pos,U); 
    bool passOnX=true; // whether to multiply next tensor by the transformation
    //if((right&&pos==lim-1)||(!right&&pos==lim+1)||(pos==lim))
    //passOnX=true;
    //    cout<<"About to apply gauge "<<dir<<" to "<<pos<<" and pass TErm="<<passOnX<<endl;
    init.gaugeCond(pos,dir,passOnX);
    // notify temporary storage as appropriate
    for(int l=pos+1;l<nrsites;l++)
      tmpMgr.outofdateL(l);
    for(int l=pos-1;l>=0;l--)
      tmpMgr.outofdateR(l);
    if(right) pos++;
    else pos--;
    // Check if sweep is done
    done=right?pos>=lim:pos<=lim;
    energy=val;
  } // while(!done)
}

// Auxiliary funtion to see what's going on!

double costOrthogonalize(MPS& current,const MPS& init,
			     const vector<const MPS*>& kets,
			     const vector<double>& penalties){
  Contractor& contractor=Contractor::theContractor(); 
  int nrkets=kets.size();
  vector<complex_t> overlaps;
  for(int k=0;k<nrkets;k++){
    overlaps.push_back(contractor.contract(current,*kets[k]));
  }
  double result=real(contractor.contract(current,current)+contractor.contract(init,init))
    -2*real(contractor.contract(current,init));
  for(int l=0;l<nrkets;l++){
    result+=penalties[l]*real(conjugate(overlaps[l])*overlaps[l]);
  }
  return result;
}

void  Contractor::orthogonalize(MPS& init,int D,
				const vector<const MPS*>& kets,
				const vector<double>& penalties){
  // check length of vectors
  int nrterms=kets.size(); 
  int nrsites=init.getLength();
  for(int nr_=0;nr_<nrterms;nr_++)
    if(nrsites!=kets[nr_]->getLength()){
      cout<<"ERROR: Not all MPS in orthogonalize have the same number of sites"
	  <<endl;
      exit(1);
    }
  MPS aux(init); // to use in the algorithm (kept const)
  aux.gaugeCond('R',0);aux.gaugeCond('L',0); // keep orogonal norm
  double normFact=aux.getNormFact(); // I will need to renormalize at the end
  double constantN=real(contract(aux,aux)); // contrib to cost function (should be the sqare of prev)
  //  cout<<"normFact of init vec is "<<normFact<<", squared "<<normFact*normFact<<" and constantN "<<constantN<<endl; // just checking
  if(D<aux.getBond()){
    // cout<<"WARNING: Contractor::orthogonalize called with D("<<D<<") smaller "
    // 	<<"than in the starting MPS (D_init="<<aux.getBond()<<"). Will truncate "
    // 	<<"the initial MPS!"<<endl;
    init=MPS(aux,D);
  }
  init.gaugeCond('R',1); // this first gauge ensures right Dr dimensions
  //cout<<"*Applied ritht gauge condition"<<endl;
  init.gaugeCond('L',1); /* this ensures proper Dl and allows starting 
			    from left.*/
  //cout<<"*Applied left gauge condition"<<endl;
  // Initialize temporary storage with gauge condition=true
  bool gauge=true;
  vector<TmpManager*> tmpMgr;
  for(int nr_=0;nr_<nrterms+1;nr_++)
    tmpMgr.push_back(new TmpManager(nrsites,gauge));

  // Start optimization: sweep to right and left until convergence
  bool done=0;
  bool right=1;
  int pos=0;
  double distance=1000*nrterms*2*constantN; // initially large value for cosat function (not truely a distance here!)
  //  distance=costOrthogonalize(init,aux,kets,penalties); // can be used to debug, but a bit more costly 
  double olddistance=distance; // to save from round to round
    //  cout<<"Cost function before starting orthogonalize: "<<distance<<endl;
  int round=1;

  vector<mwArray> Mkets(nrterms); // space for each of the effective projectors
  mwArray M,newA; // the one for the true vector
  complex_t contr1(ONE_c),contr2(ONE_c);
    
  while(!done){
    /* Optimize matrix for position pos */
    //cout<<"Contractor::orthogonalize for pos "<<pos<<", round "<<round<<endl;
    // Calculate M(pos) matrices (N not needed)
    calculateL0(pos,*tmpMgr.at(0),aux,init,gauge);
    calculateR0(pos,*tmpMgr.at(0),aux,init,gauge);    
    contractMketId(M,tmpMgr.at(0)->operL[pos],aux.getA(pos),
		   tmpMgr.at(0)->operR[pos]);
    M.multiplyLeft(normFact*ONE_c);
    mwArray matV(Indices(M.getDimension(0),nrterms)); // space for the effective vectors as columns
    for(int nr_=0;nr_<nrterms;nr_++){
      // 1. Obtain left and right terms for everything
      calculateL0(pos,*tmpMgr.at(nr_+1),*kets[nr_],init,gauge);
      calculateR0(pos,*tmpMgr.at(nr_+1),*kets[nr_],init,gauge);
      // 2. Contract the terms appropriately for M matrix
      //mwArray vecKet;
      //      contractMketId(vecKet,tmpMgr.at(nr_+1)->operL[pos],kets[nr_]->getA(pos),
      contractMketId(Mkets[nr_],tmpMgr.at(nr_+1)->operL[pos],kets[nr_]->getA(pos),
		     tmpMgr.at(nr_+1)->operR[pos]);
      Mkets[nr_]=kets[nr_]->getNormFact()*Mkets[nr_];
      // save the column (is there sth more efficient??)
      for(int k=0;k<Mkets[nr_].getDimension(0);k++)
	matV.setElement(sqrt(penalties[nr_])*Mkets[nr_].getElement(Indices(k,0)),Indices(k,nr_)); // to column nr_
      //cout<<"Found vector "<<nr_<<":"<<Mkets[nr_]<<endl;
    }
    if(pos==0&&round==1){
      if(M.isNull()){ // no overlap with target vector!
	init.perturb(1.);
	for(int l=0;l<nrsites;l++){
	  for(int nr_=0;nr_<nrterms+1;nr_++){
	    tmpMgr.at(nr_)->outofdateL(l);
	    tmpMgr.at(nr_)->outofdateR(l);
	  }
	}
	continue;
      }
    }
    // 3. Solve the local problem. In this case, I have to solve
    // (N_eff +\sum \lambda_i |V_eff_i><V_eff_i|)A =M
    // with gauge, N_eff=Id, but I still need to invert the lhs
    // instead of calling solveSite, which inverts the full matrix Neff, I will try to
    //diagonalize exactly the matrix of projectors, to exploit the fact that it has low rank
    // and I know its expression in the non-orthogonal basis Mkets

    // cout<<"matV now:"<<matV<<endl;
    // Compute overlap matrix, and using the loop also the overlaps with M
    mwArray Over=Hconjugate(matV)*matV;
    mwArray effM=Hconjugate(matV)*M;
    // Diagonalize the overlap matrix
    mwArray V;vector<complex_t> eigenOv;
    wrapper::eig(Over,eigenOv,V,true); // Over=V*eig*V^dagger
    vector<complex_t> diagVals(eigenOv);
    // Now construct the inverse matrix I need
    for(int k=0;k<nrterms;k++){
      if(eigenOv[k]!=ZERO_c)
	eigenOv[k]=ONE_c*(1./(1.+real(eigenOv[k])));
      else eigenOv[k]=ONE_c; 
    }
    mwArray result=M-matV*(V*(diag(eigenOv)*(Hconjugate(V)*effM))); 
    newA=(1./init.getNormFact())*result;
    contr1=(Hconjugate(newA)*newA+(Hconjugate(newA)*matV)*(Hconjugate(matV)*newA)).getElement(0);
    contr2=(Hconjugate(newA)*M).getElement(0);
    //    cout<<"contr1="<<contr1<<", contr2="<<contr2<<endl;
    // This is a trick to use the checkDistance function!
     double newdist=constantN+real(contr1)*init.getNormFact()*init.getNormFact()-2*real(contr2)*init.getNormFact();
     // cout<<"Contractor::orthogonalize pos "<<pos<<", round "<<round<<", distance=="<<distance
     //  	<<", newdistance="<<newdist<<endl;

    /** Now we have to substitute the solved A, but I check first that
	the distance is not increasing (only allowed within numerical
	precision) */
    // NEW
    if((newdist-distance>1E-12)&&(newdist-distance>abs(distance)*1E-10)){
      if((newdist-distance>abs(distance)*convtol)){
	cout<<"ERROR! Distance is increasing in round "<<round<<", pos "<<pos
	    <<" from "<<distance<<" to "<<newdist<<", rel. change "
	    <<(newdist-distance)/distance<<endl;
	      exit(212);
      }
      else{
	cout<<"WARNING!!! at site "<<pos<<" as distance was "<<distance
	    <<" and now is "<<newdist<<", rel chg is "<<(newdist-distance)/distance<<endl;
      }
    }
    init.setA(pos,newA);
    distance=newdist;

    double oldNormFact=init.getNormFact();
    bool keepF=false;
    // only if at the edge, I keep the normalization in the norm factor
    if((right&&pos==nrsites-1)||((!right)&&pos==0)) keepF=true;
    if(right){
      init.gaugeCond(pos,'R',keepF);
    }
    else{
      init.gaugeCond(pos,'L',keepF); 
    }
    int posL=pos+1; // first L term altered
    int posR=pos-1; // first R term altered
    for(int nr_=0;nr_<nrterms+1;nr_++){
      for(int l=posL;l<nrsites;l++) tmpMgr.at(nr_)->outofdateL(l);
      for(int l=posR;l>=0;l--) tmpMgr.at(nr_)->outofdateR(l);
    }

    // Decide the convergence only at the end of sweep
    if(((right)&&(pos==nrsites-1))||
       ((!right)&&(pos==0))){
      done=checkDistanceConvergence(olddistance,contr1,contr2,
				    oldNormFact,sqrt(constantN),round);
      if(!done){  // turn around and continue
	if(right){
	  init.gaugeCond(pos,'L',true);pos--;
	}
	else{
	  init.gaugeCond(pos,'R',true);pos++;
	}
	olddistance=distance;
	right=!right;
	round++;
      }
    }
    else{ // not in the edge=> just advance in the same direction
	pos=right?pos+1:pos-1;
      }    
    } // while(!done)
  for(int nr_=0;nr_<nrterms+1;nr_++){
    tmpMgr.at(nr_)->clear();
    delete tmpMgr[nr_];
  }
  tmpMgr.clear();
  //cout<<"Distance going out: "<<distance<<endl;
}

double costOrthogonalize(MPS& current,const MPO& mpo,const MPS& orig,
			     const vector<const MPS*>& kets,
			     const vector<double>& penalties){
  Contractor& contractor=Contractor::theContractor(); 
  int nrkets=kets.size();
  vector<complex_t> overlaps;
  for(int k=0;k<nrkets;k++){
    overlaps.push_back(contractor.contract(current,*kets[k]));
  }
  double result=real(contractor.contract(current,current))
    -2*real(contractor.contract(orig,mpo,current));
  //result+=real(contractor.contract2(mpo,orig));
  for(int l=0;l<nrkets;l++){
    result+=penalties[l]*real(conjugate(overlaps[l])*overlaps[l]);
  }
  return result;
}
void Contractor::orthogonalize(MPS& init,int D,const MPO& ops,const MPS& orig,
			       const vector<const MPS*>& kets,
			       const vector<double>& penalties){
  // check length of vectors
  int nrterms=kets.size(); 
  int nrsites=init.getLength();
  for(int nr_=0;nr_<nrterms;nr_++)
    if(nrsites!=kets[nr_]->getLength()){
      cout<<"ERROR: Not all MPS in orthogonalize have the same number of sites"
	  <<endl;
      exit(1);
    }
  double constantN=0.; //real(contract2(ops,orig)); // contrib to cost function
  double normFact=orig.getNormFact(); // I will need to renormalize at the end
  if(D<init.getBond()){
    // cout<<"WARNING: Contractor::orthogonalize called with D("<<D<<") smaller "
    // 	<<"than in the starting MPS (D_init="<<aux.getBond()<<"). Will truncate "
    // 	<<"the initial MPS!"<<endl;
    init=MPS(init,D);
  }
  init.gaugeCond('R',1); // this first gauge ensures right Dr dimensions
  //cout<<"*Applied ritht gauge condition"<<endl;
  init.gaugeCond('L',1); /* this ensures proper Dl and allows starting 
			    from left.*/
  //cout<<"*Applied left gauge condition"<<endl;
  // Initialize temporary storage with gauge condition=true
  bool gauge=true;
  vector<TmpManager*> tmpMgr;
  for(int nr_=0;nr_<nrterms+1;nr_++)
    tmpMgr.push_back(new TmpManager(nrsites,gauge));

  // Start optimization: sweep to right and left until convergence
  bool done=0;
  bool right=1;
  int pos=0;
  double distance=1000*nrterms*2*constantN; // initially large value
  //  distance=costOrthogonalize(init,ops,orig,kets,penalties);
  double olddistance=distance; // to save from round to round
  //  cout<<"Distance before starting orthogonalize with MPO: "<<distance<<endl;
  int round=1;

  vector<mwArray> Mkets(nrterms); // space for each of the effective projectors
  mwArray M,newA; // the one for the true Op*vector
  complex_t contr1(ONE_c),contr2(ONE_c);
    
  while(!done){
    /* Optimize matrix for position pos */
    // cout<<"Contractor::orthogonalize(MPO) for pos "<<pos
    //  	<<", round "<<round<<", distance "<<distance<<endl;
    // Calculate M(pos) matrices (N not needed)
    calculateL(pos,*tmpMgr.at(0),ops,orig,init,gauge);
    calculateR(pos,*tmpMgr.at(0),ops,orig,init,gauge);    
    calculateMket(M,tmpMgr.at(0)->operL[pos],orig.getA(pos),ops.getOp(pos),
		  tmpMgr.at(0)->operR[pos]);
    M.multiplyLeft(normFact*ONE_c);
    mwArray matV(Indices(M.getDimension(0),nrterms)); // space for the effective vectors as columns
    for(int nr_=0;nr_<nrterms;nr_++){
      // 1. Obtain left and right terms for everything
      calculateL0(pos,*tmpMgr.at(nr_+1),*kets[nr_],init,gauge);
      calculateR0(pos,*tmpMgr.at(nr_+1),*kets[nr_],init,gauge);
      // 2. Contract the terms appropriately for M matrix
      //mwArray vecKet;
      //      contractMketId(vecKet,tmpMgr.at(nr_+1)->operL[pos],kets[nr_]->getA(pos),
      contractMketId(Mkets[nr_],tmpMgr.at(nr_+1)->operL[pos],kets[nr_]->getA(pos),
		     tmpMgr.at(nr_+1)->operR[pos]);
      Mkets[nr_]=kets[nr_]->getNormFact()*Mkets[nr_];
      // save the column (is there sth more efficient??)
      for(int k=0;k<Mkets[nr_].getDimension(0);k++)
	matV.setElement(sqrt(penalties[nr_])*Mkets[nr_].getElement(Indices(k,0)),Indices(k,nr_)); // to column nr_
      //cout<<"Found vector "<<nr_<<":"<<Mkets[nr_]<<endl;
    }
    if(pos==0&&round==1){
      if(M.isNull()){ // no overlap with target vector!
	init.perturb(1.);
	for(int l=0;l<nrsites;l++){
	  for(int nr_=0;nr_<nrterms+1;nr_++){
	    tmpMgr.at(nr_)->outofdateL(l);
	    tmpMgr.at(nr_)->outofdateR(l);
	  }
	}
	continue;
      }
    }
    // 3. Solve the local problem. In this case, I have to solve
    // (N_eff +\sum \lambda_i |V_eff_i><V_eff_i|)A =M
    // with gauge, N_eff=Id, but I still need to invert the lhs
    // instead of calling solveSite, which inverts the full matrix Neff, I will try to
    //diagonalize exactly the matrix of projectors, to exploit the fact that it has low rank
    // and I know its expression in the non-orthogonal basis Mkets

    // Compute overlap matrix, and using the loop also the overlaps with M
    mwArray Over=Hconjugate(matV)*matV;
    //Over=.5*(Over+Hconjugate(Over)); // Unnecessary, since by construction it is HErmitian!
    // take Hermitian part if numerical errors have spoiled this
    mwArray effM=Hconjugate(matV)*M;
    // Diagonalize the overlap matrix
    mwArray V;vector<complex_t> eigenOv;
    wrapper::eig(Over,eigenOv,V,true); // Over=V*eig*V^dagger
    vector<complex_t> diagVals(eigenOv);
    //    cout<<"Diagonalized the overlap matrix, with size "<<Over.getDimension(0)<<" and eigenvals "<<eigenOv<<endl;
    // Now construct the inverse matrix I need
    for(int k=0;k<nrterms;k++){
      if(eigenOv[k]!=ZERO_c)
	eigenOv[k]=ONE_c*(1./(1.+real(eigenOv[k])));
      else eigenOv[k]=ONE_c; 
    }
    //mwArray result=M-matV*(V*(diag(eigenOv)*(Hconjugate(V)*effM))); 
    //newA=(1./init.getNormFact())*result;
    newA=(-1./init.getNormFact())*(matV*(V*(diag(eigenOv)*(Hconjugate(V)*effM)))-M); 
    contr1=(Hconjugate(newA)*newA+(Hconjugate(newA)*matV)*(Hconjugate(matV)*newA)).getElement(0);
    contr2=(Hconjugate(newA)*M).getElement(0);
    //    cout<<"contr1="<<contr1<<", contr2="<<contr2<<endl;

    double newdist=constantN+real(contr1)*init.getNormFact()*init.getNormFact()-2*real(contr2)*init.getNormFact();
    //cout<<"Contractor::orthogonalize pos "<<pos<<", round "<<round<<", distance=="<<distance
    //	<<", newdistance="<<newdist<<endl;

    /** Now we have to substitute the solved A, but I check first that
	the distance is not increasing (only allowed within numerical
	precision) */
    // NEW
    if((newdist-distance>1E-12)&&(newdist-distance>abs(distance)*1E-10)){
      if((newdist-distance>abs(distance)*convtol)){
	cout<<"ERROR! Distance is increasing in round "<<round<<", pos "<<pos
	    <<" from "<<distance<<" to "<<newdist<<", rel. change "
	    <<(newdist-distance)/distance<<endl;
	      exit(212);
      }
      else{
	cout<<"WARNING!!! at site "<<pos<<" as distance was "<<distance
	    <<" and now is "<<newdist<<", rel chg is "<<(newdist-distance)/distance<<endl;
      }
    }
    init.setA(pos,newA);
    distance=newdist;
    double oldNormFact=init.getNormFact();
    bool keepF=false;
    // only if at thee dge, I keep the normalization in the norm factor
    if((right&&pos==nrsites-1)||(!right&&pos==0)) keepF=true;
    if(right){
      init.gaugeCond(pos,'R',keepF);
    }
    else{
      init.gaugeCond(pos,'L',keepF); 
    }
    int posL=pos+1; // first L term altered
    int posR=pos-1; // first R term altered
    for(int nr_=0;nr_<nrterms+1;nr_++){
      for(int l=posL;l<nrsites;l++) tmpMgr.at(nr_)->outofdateL(l);
      for(int l=posR;l>=0;l--) tmpMgr.at(nr_)->outofdateR(l);
    }

    // Decide the convergence only at the end of sweep
    if(((right)&&(pos==nrsites-1))||
       ((!right)&&(pos==0))){
      //done=(abs(olddistance-distance)<max(convtol*abs(olddistance),NUMPREC))?1:0;
      done=checkDistanceConvergence(olddistance,contr1,contr2,
				    oldNormFact,sqrt(constantN),round);
      if(!done){  // turn around and continue
	if(right){
	  init.gaugeCond(pos,'L',true);pos--;
	}
	else{
	  init.gaugeCond(pos,'R',true);pos++;
	}
	olddistance=distance;
	right=!right;
	round++;
      }
    }
    else{ // not in the edge=> just advance in the same direction
	pos=right?pos+1:pos-1;
      }    
    } // while(!done)
  for(int nr_=0;nr_<nrterms+1;nr_++){
    tmpMgr.at(nr_)->clear();
    delete tmpMgr[nr_];
  }
  tmpMgr.clear();
  //cout<<"Distance going out: "<<distance<<endl;

}

double costOrthogonalizeSum(const MPS& current,
			    const vector<const MPS*>& vecs,
			    const vector<double>& betas,
			    const vector<const MPS*>& kets,
			    const vector<double>& penalties){
  Contractor& contractor=Contractor::theContractor(); 
  int nrvecs=vecs.size();
  int nrkets=kets.size();
  vector<complex_t> overlaps;
  for(int k=0;k<nrkets;k++){
    overlaps.push_back(contractor.contract(current,*kets[k]));
  }
  double result=real(contractor.contract(current,current));
  for(int k=0;k<nrvecs;k++){
    complex_t overlapK=contractor.contract(current,*vecs[k]);
    result+=-2*betas[k]*real(overlapK);
  }
  for(int l=0;l<nrkets;l++){
    result+=penalties[l]*real(conjugate(overlaps[l])*overlaps[l]);
  }
  return result;
}

void  Contractor::orthogonalizeSum(MPS& init,int D,
				   const vector<const MPS*>& vecs,const vector<complex_t>& beta,
				   const vector<const MPS*>& kets,
				   const vector<double>& penalties){
  // check length of vectors
  int nrterms=kets.size(); 
  int nrvecs=vecs.size(); 
  int nrsites=init.getLength();
  for(int nr_=0;nr_<nrterms;nr_++)
    if(nrsites!=kets[nr_]->getLength()){
      cout<<"ERROR: Not all MPS in orthogonalizeSum have the same number of sites"
	  <<endl;
      exit(1);
    }
  double constantN=0.; //should be the contrib of |sum beta*vecs(k)|^2  to the cost function, but I am ignoring the constant!
  if(D<init.getBond()){
    // cout<<"WARNING: Contractor::orthogonalize called with D("<<D<<") smaller "
    // 	<<"than in the starting MPS (D_init="<<aux.getBond()<<"). Will truncate "
    // 	<<"the initial MPS!"<<endl;
    init=MPS(init,D);
  }
  init.gaugeCond('R',1); // this first gauge ensures right Dr dimensions
  //cout<<"*Applied right gauge condition"<<endl;
  init.gaugeCond('L',1); /* this ensures proper Dl and allows starting 
			    from left.*/
  //cout<<"*Applied left gauge condition"<<endl;
  // Initialize temporary storage with gauge condition=true
  bool gauge=true;
  vector<TmpManager*> tmpMgrV; // for the vectors in the sum
  for(int nr_=0;nr_<nrvecs+1;nr_++)
    tmpMgrV.push_back(new TmpManager(nrsites,gauge));
  vector<TmpManager*> tmpMgr; // for the penalty terms
  for(int nr_=0;nr_<nrterms+1;nr_++)
    tmpMgr.push_back(new TmpManager(nrsites,gauge));

  // Start optimization: sweep to right and left until convergence
  bool done=0;
  bool right=1;
  int pos=0;
  double distance=1000*(nrterms*2+nrvecs*2); // initially large value
  //  distance=costOrthogonalize(init,aux,kets,penalties); // can be used to debug, but a bit more costly 
  double olddistance=distance; // to save from round to round
  //cout<<"Distance before starting orthogonalize: "<<distance<<endl;
  int round=1;

  vector<mwArray> Mkets(nrterms); // space for each of the effective projectors
  mwArray Mloc,M,newA; // the one for the true vector
  complex_t contr1(ONE_c),contr2(ONE_c);
    
  while(!done){
    /* Optimize matrix for position pos */
    //cout<<"Contractor::orthogonalize for pos "<<pos<<", round "<<round<<endl;
    // Calculate M(pos) matrices (N not needed)
    for(int nr_=0;nr_<nrvecs;nr_++){
      // 1. Obtain left and right terms for everything
      calculateL0(pos,*tmpMgrV.at(nr_),*vecs[nr_],init,gauge);
      calculateR0(pos,*tmpMgrV.at(nr_),*vecs[nr_],init,gauge);
      // 2. Contract the terms appropriately for M matrix
      // cout<<"After calculateL->"<<tmpMgr.operL[pos]<<", calculateR->"
      // <<tmpMgr.operR[pos]<<", with op->"<<ops.getOp(pos)<<endl;
      contractMketId(Mloc,tmpMgrV.at(nr_)->operL[pos],vecs[nr_]->getA(pos),
		     tmpMgrV.at(nr_)->operR[pos]);
      if(nr_==0) // M not initialized
	M=vecs[nr_]->getNormFact()*beta[nr_]*Mloc;
      else
	M=M+vecs[nr_]->getNormFact()*beta[nr_]*Mloc;
    }
    mwArray matV(Indices(M.getDimension(0),nrterms)); // space for the effective vectors as columns
    for(int nr_=0;nr_<nrterms;nr_++){
      // 1. Obtain left and right terms for everything
      calculateL0(pos,*tmpMgr.at(nr_+1),*kets[nr_],init,gauge);
      calculateR0(pos,*tmpMgr.at(nr_+1),*kets[nr_],init,gauge);
      // 2. Contract the terms appropriately for M matrix
      //mwArray vecKet;
      //      contractMketId(vecKet,tmpMgr.at(nr_+1)->operL[pos],kets[nr_]->getA(pos),
      contractMketId(Mkets[nr_],tmpMgr.at(nr_+1)->operL[pos],kets[nr_]->getA(pos),
		     tmpMgr.at(nr_+1)->operR[pos]);
      Mkets[nr_]=kets[nr_]->getNormFact()*Mkets[nr_];
      // save the column (is there sth more efficient??)
      for(int k=0;k<Mkets[nr_].getDimension(0);k++)
	matV.setElement(sqrt(penalties[nr_])*Mkets[nr_].getElement(Indices(k,0)),Indices(k,nr_)); // to column nr_
      //cout<<"Found vector "<<nr_<<":"<<Mkets[nr_]<<endl;
    }
    if(pos==0&&round==1){
      if(M.isNull()){ // no overlap with target vector!
	init.perturb(1.);
	for(int l=0;l<nrsites;l++){
	  for(int nr_=0;nr_<nrterms+1;nr_++){
	    tmpMgr.at(nr_)->outofdateL(l);
	    tmpMgr.at(nr_)->outofdateR(l);
	  }
	}
	continue;
      }
    }
    // 3. Solve the local problem. In this case, I have to solve
    // (N_eff +\sum \lambda_i |V_eff_i><V_eff_i|)A =M
    // with gauge, N_eff=Id, but I still need to invert the lhs
    // instead of calling solveSite, which inverts the full matrix Neff, I will try to
    //diagonalize exactly the matrix of projectors, to exploit the fact that it has low rank
    // and I know its expression in the non-orthogonal basis Mkets

    // cout<<"matV now:"<<matV<<endl;
    // Compute overlap matrix, and using the loop also the overlaps with M
    mwArray Over=Hconjugate(matV)*matV;
    mwArray effM=Hconjugate(matV)*M;
    // Diagonalize the overlap matrix
    mwArray V;vector<complex_t> eigenOv;
    wrapper::eig(Over,eigenOv,V,true); // Over=V*eig*V^dagger
    vector<complex_t> diagVals(eigenOv);
    // Now construct the inverse matrix I need
    for(int k=0;k<nrterms;k++){
      if(eigenOv[k]!=ZERO_c)
	eigenOv[k]=ONE_c*(1./(1.+real(eigenOv[k])));
      else eigenOv[k]=ONE_c; 
    }
    mwArray result=M-matV*(V*(diag(eigenOv)*(Hconjugate(V)*effM))); 
    newA=(1./init.getNormFact())*result;
    contr1=(Hconjugate(newA)*newA+(Hconjugate(newA)*matV)*(Hconjugate(matV)*newA)).getElement(0);
    contr2=(Hconjugate(newA)*M).getElement(0);
    //    cout<<"contr1="<<contr1<<", contr2="<<contr2<<endl;

    double newdist=constantN+real(contr1)*init.getNormFact()*init.getNormFact()-2*real(contr2)*init.getNormFact();
    // cout<<"Contractor::orthogonalize pos "<<pos<<", round "<<round<<", distance=="<<distance
    //  	<<", newdistance="<<newdist<<endl;

    /** Now we have to substitute the solved A, but I check first that
	the distance is not increasing (only allowed within numerical
	precision) */
    // NEW
    if((newdist-distance>1E-12)&&(newdist-distance>abs(distance)*1E-10)){
      if((newdist-distance>abs(distance)*convtol)){
	cout<<"ERROR! Distance is increasing in round "<<round<<", pos "<<pos
	    <<" from "<<distance<<" to "<<newdist<<", rel. change "
	    <<(newdist-distance)/distance<<endl;
	      exit(212);
      }
      else{
	cout<<"WARNING!!! at site "<<pos<<" as distance was "<<distance
	    <<" and now is "<<newdist<<", rel chg is "<<(newdist-distance)/distance<<endl;
      }
    }
    init.setA(pos,newA);
    distance=newdist;
    double oldNormFact=init.getNormFact();
    bool keepF=false;
    // only if at thee dge, I keep the normalization in the norm factor
    if((right&&pos==nrsites-1)||(!right&&pos==0)) keepF=true;
    if(right){
      init.gaugeCond(pos,'R',keepF);
    }
    else{
      init.gaugeCond(pos,'L',keepF); 
    }
    int posL=pos+1; // first L term altered
    int posR=pos-1; // first R term altered
    for(int nr_=0;nr_<nrterms+1;nr_++){
      for(int l=posL;l<nrsites;l++) tmpMgr.at(nr_)->outofdateL(l);
      for(int l=posR;l>=0;l--) tmpMgr.at(nr_)->outofdateR(l);
    }
    for(int nr_=0;nr_<nrvecs+1;nr_++){
      for(int l=posL;l<nrsites;l++) tmpMgrV.at(nr_)->outofdateL(l);
      for(int l=posR;l>=0;l--) tmpMgrV.at(nr_)->outofdateR(l);
    }

    // Decide the convergence only at the end of sweep
    if(((right)&&(pos==nrsites-1))||
       ((!right)&&(pos==0))){
      done=checkDistanceConvergence(olddistance,contr1,contr2,
				    init.getNormFact(),sqrt(constantN),round);
      //init.getNormFact(),1.,round);
      //done=(abs(olddistance-distance)<max(convtol*abs(olddistance),NUMPREC))?1:0;
      if(!done){  // turn around and continue
	if(right){
	  init.gaugeCond(pos,'L',true);pos--;
	}
	else{
	  init.gaugeCond(pos,'R',true);pos++;
	}
	olddistance=distance;
	right=!right;
	round++;
      }
    }
    else{ // not in the edge=> just advance in the same direction
	pos=right?pos+1:pos-1;
      }    
    } // while(!done)
  for(int nr_=0;nr_<nrterms+1;nr_++){
    tmpMgr.at(nr_)->clear();
    delete tmpMgr[nr_];
  }
  for(int nr_=0;nr_<nrvecs+1;nr_++){
    tmpMgrV.at(nr_)->clear();
    delete tmpMgrV[nr_];
  }
  tmpMgr.clear();
  //cout<<"Distance going out: "<<distance<<endl;
}

void Contractor::blockOperatorMultiSite(mwArray& blockOp,const MPO& ops,int k,int pos){
  cout<<"Contractor::blockOperatorMultiSite from sites "<<pos<<" to "<<k+pos-1<<endl;
  blockOp=ops.getOp(pos).getFullData();
  //  cout<<"blockOp:"<<blockOp.getDimensions()<<endl;
  for(int l=1;l<k;l++){
    mwArray nextOp=ops.getOp(pos+l).getFullData();
    // contract together the previous one and the local term
    Indices dimsOld=blockOp.getDimensions();
    Indices dimsNew=nextOp.getDimensions();
    blockOp.reshape(Indices(-1,dimsOld[3]));
    nextOp.permute(Indices(2,1,3,4));
    nextOp.reshape(Indices(dimsNew[1],-1));
    blockOp.multiplyRight(nextOp);
    blockOp.reshape(Indices(dimsOld[0],dimsOld[1],dimsOld[2],dimsNew[0],dimsNew[2],dimsNew[3]));
    blockOp.permute(Indices(1,4,2,3,5,6));
    blockOp.reshape(Indices(dimsOld[0]*dimsNew[0],dimsOld[1],dimsOld[2]*dimsNew[2],dimsNew[3]));
    //cout<<"After reshaping:"<<blockOp.getDimensions()<<endl;
  }
}

double Contractor::gradientMinVariance(mwArray& A,TensorMultiplier* multiH2,
				       TensorMultiplier* multiH,double E0,double penH){
  mwArray P; //,Plast(-1.*Glast); // p_k and p_{k-1}
  double value,lastvalue(1E3);
  double alpha0=5;// TODO: fix!!!!
  int k=0;
  bool stop=false;
  while(!stop){
    double alpha=alpha0;
    // Compute the value of the gradient now
    mwArray gk;
    // 1) Term from MPO^2
    multiH2->product(reshape(A,Indices(-1,1)),gk);
    gk.reshape(Indices(1,-1)); // a row
    mwArray valH(A);
    valH.reshape(Indices(1,-1));valH.Hconjugate();
    value=real((gk*valH).getElement(0));
    // 2) Derivative of the term with the Hamiltonian
    mwArray Gh;
    multiH->product(reshape(A,Indices(-1,1)),Gh);
    Gh.reshape(Indices(1,-1)); // a row
    valH.multiplyLeft(Gh); // value of the contr. with H
    mwArray valH0(valH); // just current E, if I want to include -<H>^2 explicitly in the cost
    valH=valH-mwArray(E0*ONE_c);
    //cout<<"Now contr with H="<<valH<<"=>+"<<penH*real((valH*valH).getElement(0))<<endl;
    gk=gk+2*penH*valH*Gh; // current gradient
    gk=gk-2*valH0*Gh; //if I want to include -<H>^2 explicitly in the cost
    value+=penH*real((valH*valH).getElement(0));
    value-=real((valH0*valH0).getElement(0)); // if I want to include -<H>^2 explicitly in the cost
    stop=(abs(1.-value/lastvalue)<1E-10)||k>1E6;    
    if(isnan_double(value)){
      cout<<"Contractor::gradientMinVariance NAN!!! iter "<<k<<" value="<<value<<" chg="<<abs(1.-value/lastvalue)<<endl;
      exit(1);
    }
    if(!stop){
      lastvalue=value;
      P=-1.*gk;
      // Now I should decide on alpha, by choosing the min along direction P****
      bool decreasing;
      do{
	mwArray aux(A+alpha*reshape(P,A.getDimensions()));
	mwArray normA(aux);normA.reshape(Indices(-1,1));
	normA.multiplyLeft(Hconjugate(normA));
	aux.multiplyLeft(1./sqrt(real(normA.getElement(0)))*ONE_c);

	// Need to evaluate again the cost function
	double newval=0.;
	{
	  mwArray gk2,gk1;
	  multiH2->product(reshape(aux,Indices(-1,1)),gk2);
	  multiH->product(reshape(aux,Indices(-1,1)),gk1);
	  gk2.reshape(Indices(1,-1)); // a row
	  gk1.reshape(Indices(1,-1)); // a row
	  mwArray valH(aux);
	  valH.reshape(Indices(1,-1));valH.Hconjugate();
	  complex_t val2=(gk2*valH).getElement(0);
	  complex_t val1=(gk1*valH).getElement(0);
	  newval+=real(val2)+penH*real((val1-E0*ONE_c)*(val1-E0*ONE_c));
	  newval-=real(val1*val1);// if I want to include -<H>^2 explicitly in the cost
	}
	decreasing=newval<lastvalue;
	if(decreasing){
	  A=aux; //A+alpha*reshape(P,A.getDimensions());
	  lastvalue=newval;
	}
	else{
	  alpha=0.5*alpha;
	}
      } while(decreasing||alpha>1E-5);
      k++;
    }
  }
  return value;
}



void Contractor::minimizeVariance(MPS& init,const int D,const MPO& ops,
				  const double E0,double penH,double* lambda,
				  const string& tmpfile,int freqSv,int stopRound){

  cout<<"minimizeVariance called with D="<<D<<", init.D= "
      <<init.getBond()<<" with <ops>="<<E0<<endl;
#ifndef USING_PRIMME
  cout<<"Compiled without PRIMME support => eigensolver = ARPACK"<<endl;
  solver=arpack;
#endif
  int nrsites=ops.getLength();  
  vector<int> findims=ops.getDimensions();
  init.adjustPhysDimensions(findims);
  if(D>0){
    if(D>init.getBond())
      init.increaseBondDimension(D);
    else if(D<init.getBond()){
      cout<<"WARNING: Contractor::minimizeVariance called with D("<<D<<") smaller "
	  <<"than in the starting MPS (D_init="<<init.getBond()<<"). Will truncate "
	  <<"the initial MPS!"<<endl;
      init=MPS(init,D);
    }
  }
  init.gaugeCond('R',1); // this first gauge ensures right Dr dimensions
  init.gaugeCond('L',1); // this ensures proper Dl and allows starting 
			    // from left.
  bool done=0;int round=0;
  bool right=1; int pos=0; //nrsites-1;
  // Initial value of the cost function *****////
  double valE=real(contract(init,ops,init));
  double energy=real(contract2(ops,init))+penH*(valE-E0)*(valE-E0);
  energy-=valE*valE; // if including -H^2
  // I make a joined operator for the MPO squared
  const MPO* ptrs[]={&ops,&ops};
  MPO ops2(nrsites);
  MPO::join(2,ptrs,ops2); // Joined
  
  // Two TmpManager instances to keep the (regular) contraction with H,
  // and the tmp terms for H^2 
  bool gauge=true;
  TmpManager tmpMgrH(nrsites,gauge);
  TmpManager tmpMgrH2(nrsites,gauge);
  clock_t start=clock();clock_t finish=clock();
  double val=energy;
  while(!done){
    if(pos==0&&round%20==1){
      // finish=clock();
       cout<<"Round number "<<round<<", pos "<<pos
       	  <<", cost function="<<val; // energy;
       cout<<endl;
	    // cout<<" time "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
	      // start=clock();
    }
    // 1. Obtain left and right terms for everything. In this case,
    // this includes the partial contractions of the result MPS with
    // the MPO, and also its contraction with MPO^2
    calculateL(pos,tmpMgrH,ops,init,init,gauge);
    calculateR(pos,tmpMgrH,ops,init,init,gauge);
    calculateL(pos,tmpMgrH2,ops2,init,init,gauge);
    calculateR(pos,tmpMgrH2,ops2,init,init,gauge);
    // Now I need to solve the local problem for A. I try a gradient
    // search, which I program adhoc for this method
    // As help, I pass the TensorMultipliers to the method
    TensorMultiplier* multi=new TensorMultiplier(tmpMgrH.operL[pos],ops.getOp(pos).getFullData(),
						 tmpMgrH.operR[pos]);
    TensorMultiplier* multi2=new TensorMultiplier(tmpMgrH2.operL[pos],ops2.getOp(pos).getFullData(),
						 tmpMgrH2.operR[pos]);
    mwArray A(init.getA(pos).getA()); 
    val=gradientMinVariance(A,multi2,multi,E0,penH);

    delete multi;delete multi2;
    //    cout<<"After the gradient min, value="<<val<<" and tr(rho^2)="<<*trR2<<endl;
     //val=linearMinRho2WithE(A,pos,d,dpur,multi,tmpMgrRho2,penH,E,penN,trR2);
     //cout<<"After the linearized min, value="<<val<<" and tr(rho^2)="<<*trR2<<endl;

    // Normalize and set the tensor
    double normA=real((Hconjugate(reshape(A,Indices(-1,1)))*reshape(A,Indices(-1,1))).getElement(0));
    init.setA(pos,(1./sqrt(normA))*A); 

    // cout<<"After setting the tensor, computed exactly: trR2="<<contractRho2Puri(init,d,dpur)
    // 	<<", <H>="<<contract(init,ops,init)<<", tr(rho)="<<contract(init,init)<<endl;
    // exit(1);
    if((val-energy>1E-8)&&(val-energy>abs(energy)*1E-10)){
      if((val-energy>abs(energy)*1E-6)){
	cout<<"Error: energy increasing in this step!! oldcost="<<energy
	    <<", new one="<<val<<", difference="<<(val-energy)
	    <<", relative increase="<<(val-energy)/(energy)<<endl;
      }
      else{
	cout<<"WARNING: energy increasing in this step!! oldcost="<<energy
	    <<", new one="<<val<<", difference="<<(val-energy)
	    <<", relative increase="<<(val-energy)/(energy)
	    <<" I CONTINUE!!"<<endl;
      }
    }
    if(right){
      //if(pos<nrsites-1)
      init.gaugeCond(pos,'R',true);
    }
    else{
      //if(pos>0)
      init.gaugeCond(pos,'L',true);
    }
    // notify temporary storage as appropriate
    for(int l=pos+1;l<nrsites;l++){
      tmpMgrH.outofdateL(l);
      tmpMgrH2.outofdateL(l);
    }
    for(int l=pos-1;l>=0;l--){
      tmpMgrH.outofdateR(l);
      tmpMgrH2.outofdateR(l);
    }  
    bool doneRound=(right&&pos==nrsites-1)||((!right)&&pos==0);
    //if (pos==1) done=1;
    // Now move one pos and check convergence if we are at the end f
    // one sweep
    if(right){ //moving right
      if(pos==nrsites-1){ //last reached: check convergence and change sense
	right=0;
	done=abs(energy-val)<convtol*abs(energy);
	if(!done){
	  //	  cout<<"Done test failed with relative change "
	  //  <<abs((energy-val)/(energy))<<endl;
	  // Check if the value is too close to zero!
	  if(abs(energy-val)<1E-16){
	    cout<<"Assuming convergence because difference is small!"<<endl;
	    done=true;
	  }
	  // if(abs(val)<1E-10){
	  //   cout<<"WARNING***** Assuming convergence because E is small!"<<endl;
	  //   done=true;
	  // }
	}
	else{
	  cout<<"done=true bc abs(energy-val)="<<abs(energy-val)<<", energy="<<energy<<" convtol="<<convtol<<endl;
	}
	energy=val;
	if(!done){
	  init.gaugeCond(pos,'L',true);pos--;
	}
	round++;
      }
      else{
	pos++;
      }
    }
    else{// moving left
      if(pos==0){ //first reached: check convergence and change sense
	right=1;
	done=abs(energy-val)<convtol*abs(energy);
	if(!done){
	  //	  cout<<"Done test failed with relative change "
	  //  <<abs((energy-val)/(energy+offset))<<" including offset "<<offset<<endl;
	  // Check if the value is too close to zero!
	  if(abs(energy-val)<1E-12){
	    cout<<"Assuming convergence because difference is small!"<<endl;
	    done=true;
	  }
	  // if(abs(val)<1E-10){
	  //   cout<<"WARNING***** Assuming convergence because E is small!"<<endl;
	  //   done=true;
	  // }
	  
	}
	else{
	  cout<<"done=true bc abs(energy-val)="<<abs(energy-val)<<", energy="<<energy<<" convtol="<<convtol<<endl;
	}
	energy=val;
	if(!done){
	  init.gaugeCond(pos,'R',true);pos++;
	}
	round++;
      }
      else{
	pos--;
      }
    }
    if(doneRound&&round>=stopRound-1){
      cout<<"Stopping iteration after round "<<round+1<<" independent of convergence"<<endl;
      done=true;
    }  
    //      cout<<"Saving Temporary resuts to "<<tmpfile<<endl;
    if(tmpfile!=""){ // Saving intermediate results (only at the end)
      bool saveNow=doneRound&&(round%freqSv==0);
      if(saveNow){
	init.exportMPS(tmpfile.data());
	//cout<<"Saved Temporary resuts to "<<tmpfile<<endl;
      }
    }
  }
}



