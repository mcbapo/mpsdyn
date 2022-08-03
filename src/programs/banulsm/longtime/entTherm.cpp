#include <math.h>
#include <iomanip>
#include <sstream>

#include "misc.h"
#include "uMPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "SpinMPO.h"
#include "Properties.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

//#define SYMM 1 // If SYMM=1, only left part is computed explicitly, and right is just the reflection

#include "IsingHamiltonian.h"

using namespace std;
using namespace shrt;

/** Evolve an infinite chain in order to test the entanglement
    dynamics. So far,it just evolves with uMPS and computes some decompositions.
*/

int d=2;
double tolSVD=1E-10;      


int maxIter=100;
double precCut=1E-8; // where to cut the eigenvalues

double computeEntropy(const mwArray& rdm);
/** Construct the effective state for the nB block, with dimensions Dlxd^nxDr
    It starts on site k0 (default 0) */
void buildBlockState(mwArray& res,uMPS& state,int k0);
void buildBlockState(mwArray& res,uMPS& state,int k0,int L);
void checkBlockState(uMPS& state,int k0,ofstream* outSVD,int lB=0);
void entropiesBlockState(uMPS& state,int k0,ofstream* outS,int lB=0);
void negativitiesBlockState(uMPS& state,int k0,ofstream* outS,int lB=0);
void channelBlock(uMPS& state,int k0,ofstream* outS,int lB);
// Auxiliary:
void buildIndexList(vector<int>& ind,int k0,int L,int nB);


/** Analyze the entanglement between a pair of spins at k0 af beyond a
    block of size lB to the right. See if a maximally entangled pair
    is popssible.  */
void distantPair(uMPS& state,int k0,int lB,ofstream* outS);
// Same between pair of blocks
void distantBlock(uMPS& state,int k0,int nS,int lB,ofstream* outS);

/** Same but for several pairs of lengths between lBmin and lBmax */
void distantPairs(uMPS& state,int k0,int lBmin,int lBmax,ofstream* outS);


int main(int argc,const char* argv[]){
  int cntr=0;
  const char* infile=argv[++cntr];
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  double J_=props.getDoubleProperty("J");
  double h_=props.getDoubleProperty("h");
  double g_=props.getDoubleProperty("g");
  int nB=props.getIntProperty("nB"); // length of the block in uMPS
  if(nB<0) nB=1; // default is 1
  int lB=props.getIntProperty("lB"); // length of the block studied
  if(lB<0) lB=nB; // default is 1
  int initSt=props.getIntProperty("initSt");
  
  double delta=props.getDoubleProperty("delta");
  int M=props.getIntProperty("M");
  int Mcut=props.getIntProperty("Mcut");
  int Dcut=props.getIntProperty("Dcut"); // the truncation 
  string outfname=props.getProperty("outputfile");
  string outSfname=props.getProperty("outputS");
  string outSVDfname=props.getProperty("outputSVD");
  //double tol=props.getDoubleProperty("tol");
  //if(tol<0) tol=1E-5;    
  int savingfreq=props.getIntProperty("savingfreq"); // frequency with
						     // which
						     // intermediate
						     // results are
						     // saved to disk
  bool savingTmp=savingfreq>0;
  string tmpdir=props.getProperty("tmpdir"); // directory to save
					     // temporary tensors
					     // (big!)
  if(savingTmp&&tmpdir.empty()){
    cout<<"ERROR! No tmpdir specified to save temporary data"<<endl;
    exit(1);
  }
  bool app=0;
  double eps=props.getDoubleProperty("eps");
  // int maxIter_=props.getIntProperty("maxIter");
  // if(maxIter_>0) maxIter=maxIter_;
  // double precCut_=props.getDoubleProperty("precCut");
  // if(precCut_>0) precCut=precCut_;


  // Get the basic piece of the evolution
  mwArray midU;
  {
    IsingHamiltonian hamH(3,d,J_,g_,h_);
    //hamH.registerTimeDependence();
    MPO _Uevol_1(3); // held just temporarily
    hamH.getUMPO(_Uevol_1,-delta*I_c,0); // the single layer normal evolution
    midU=_Uevol_1.getOp(1).getFullData();
  }
  vector<const mwArray*> Uops(nB,&midU);
  
  // and the opsto measure
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1*I_c,ZERO_c};
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c}; //Z
  mwArray sig0=identityMatrix(d);//sig0.reshape(Indices(1,d*d));
  mwArray sigX=mwArray(Indices(d,d),dataX);//sigX.reshape(Indices(1,d*d));
  mwArray sigY(Indices(d,d),dataY);//sigY.reshape(Indices(1,d*d));
  mwArray sigZ(Indices(d,d),dataZ);//sigZ.reshape(Indices(1,d*d));

  int cnt=0;
  double t=0.;

  // First of all, I initialize a uMPS and evolve it
  ofstream* out;
  if(app){
    // out=new ofstream(outfname.data(),ios::app);  *out<<setprecision(15);
  }  
  else{
    out=new ofstream(outfname.data());
    *out<<"% t\t real(tr)\t imag(tr)\t real(sigX)\t imag(sigX)\t real(sigY)\t imag(sigY)\t"
	<<" real(sigZ)\t imag(sigZ)\t real(ZZ)\t imag(ZZ)\t Dcut"<<endl;
    out->close();delete out;

    ofstream* outS=new ofstream(outSfname.data());
    //*outS<<"% t\t S(nB)\t S(L/2)"<<endl;
    *outS<<"%t\t SL\t Sph\t SR"<<endl;
    outS->close();delete outS;

    ofstream* outSVD=new ofstream(outSVDfname.data());
    *outSVD<<"% t\t svd"<<endl;
    outSVD->close();delete outSVD;

  }

  int D0=1; // starting with product state
  uMPS state(nB,d,D0);
  cout<<"Created empty initial state "<<state<<endl;
  switch(initSt){
  case 1:{state.setProductState(p_xplus);break;}
  case 2:{state.setProductState(p_yplus);break;}
  case 3:{state.setProductState(p_zero);break;}
  }
  //  srand(time(NULL));srandom(time(NULL));state.setRandomState();//TEST
  cout<<"Set initial state "<<state<<endl;

  // // Save state and operator
  // {ofstream* debu=new ofstream("state0.m");
  //   putForMatlab(*debu,state.getA(0),"A0");
  //   debu->close();delete debu;
  //   debu=new ofstream("operatorEvol.m");
  //   putForMatlab(*debu,midU,"mpoU");
  //   debu->close();delete debu;
  // }
  



  
  while(cnt<=M){
    // // if(cnt==2){
    // //   ofstream* debu=new ofstream("state2.m");
    // //   putForMatlab(*debu,state.getA(0),"A2");
    // //   debu->close();delete debu;
    // // }
    // evaluate EV and evolve one step
    // single site EV are averaged over tensors
    //cout<<"Step nr "<<cnt<<" uMPS: "<<state<<endl;
    mwArray rdm(Indices(d,d));rdm.fillWithZero();
    for(int k=0;k<nB;k++){
      mwArray auxRDM;
      state.getSingleSiteRDM(k,auxRDM);      
      rdm=rdm+auxRDM;
      //cout<<"Now RDM="<<rdm<<endl;
      // if(k==0){
      // 	ofstream* outS=new ofstream(outSfname.data(),ios::app);
      // 	*outS<<setprecision(15);
      // 	double Sb=computeEntropy(rdm);
      // 	*outS<<t<<"\t"<<Sb<<"\t"<<state.getEntropy(k)<<"\t";
      // 	outS<<endl;outS->close();delete outS;
      // }
    }
    rdm.multiplyLeft((1./nB)*ONE_c); // Avrg of nB sites
    //cout<<"Now RDM="<<rdm<<endl;
    rdm.reshape(Indices(d*d,1));
    mwArray res=reshape(sigX,Indices(1,d*d))*rdm;
    complex_t resultX=res.getElement(Indices(0,0));
    res=reshape(sigY,Indices(1,d*d))*rdm;
    complex_t resultY=res.getElement(Indices(0,0));
    res=reshape(sigZ,Indices(1,d*d))*rdm;
    complex_t resultZ=res.getElement(Indices(0,0));
    res=reshape(sig0,Indices(1,d*d))*rdm;
    complex_t trV=res.getElement(Indices(0,0));
    
    // For the ZZ value, with only one tensor I just need one calculation
    complex_t resultZZ(ZERO_c);
    vector<int> pos(2,0);pos[1]=1;
    vector<mwArray*> ops;ops.push_back(&sigZ);ops.push_back(&sigZ);
    if(nB==1){
      resultZZ=state.computeExpectationValue(2,pos,ops);
    }
    else{
      // should loop to get an average!
      for(int k=0;k<nB;k++){
	  pos[0]=k;
	  if(k<nB-1) pos[1]=k+1;
	  else pos[1]=0;
	  resultZZ+=state.computeExpectationValue(2,pos,ops);
      }
      resultZZ=(1./nB)*resultZZ;
    }
    //cout<<"resultZZ="<<resultZZ<<endl;
    // and write to file and screen
    cout<<"t="<<t<<", <O>="<<(resultX)/trV<<", trace="
	<<(trV)
	<<" E="<<(J_*(resultZZ)+g_*(resultX)+h_*(resultZ))/trV
      //<<" S(block)="<<Sblock
	<<", D="<<state.getLeftVirtualDimension(0)<<endl;
    out=new ofstream(outfname.data(),ios::app);  *out<<setprecision(15);
    *out<<t<<"\t"<<real(trV)<<"\t"<<imag(trV)
	<<"\t"<<real(resultX)<<"\t"<<imag(resultX)
	<<"\t"<<real(resultY)<<"\t"<<imag(resultY)
	<<"\t"<<real(resultZ)<<"\t"<<imag(resultZ)
	<<"\t"<<real(resultZZ)<<"\t"<<imag(resultZZ)
	<<"\t"<<state.getLeftVirtualDimension(0)
	<<endl;
    out->close();delete out;
    // Identify disentangled/thermalized dof
    // Stage0: just test what is going on.

    // Reconstruct the state for the block, and analyze its entanglement 
    //mwArray Lambda=state.getLambda(0);
    //vector<double> lVals;Lambda.getRealDiagonal(lVals);
    //cout<<"Lambda is "<<lVals<<endl;

    if(0){
    ofstream* outSVD=new ofstream(outSVDfname.data(),ios::app);
    *outSVD<<setprecision(15);
    *outSVD<<t<<"\t";
    checkBlockState(state,0,outSVD,lB);
    *outSVD<<endl;
    outSVD->close();delete outSVD;
    }

    ofstream* outS=new ofstream(outSfname.data(),ios::app);
    *outS<<setprecision(15);
    *outS<<t<<"\t";
    if(0){
    entropiesBlockState(state,0,outS,lB);
    }
    //if(0){
    //    distantPair(state,0,lB,outS);
    distantBlock(state,0,3,lB,outS);
    //}
    if(0)    channelBlock(state,0,outS,lB);
    *outS<<endl;
    outS->close();delete outS;

    // // if(cnt==M||cnt==Mcut){
    // //   if(cnt==Mcut){
    // // 	// Read a modified tensor to do some tests
    // // 	mwArray aux;
    // // 	ifstream data("newTensor.dat");
    // // 	aux.loadtext(data,0);
    // // 	cout<<"Read tensor from file: "<<aux.getDimensions()<<endl;
    // // 	mwArray aux2=state.getA(0); // previous state
    // // 	cout<<"Replacing current tensor "<<aux2.getDimensions()<<endl;
    // // 	// state.setRandomState();//
    // // 	// mwArray aux3=state.getA(0); 
    // // 	// cout<<"Read tensor now (random) "<<aux3.getDimensions()<<endl;
    // // 	state.replaceSite(0,aux);//aux3=state.getA(0); 
    // // 	//state.replaceSite(0,aux2);
    // // 	//	cout<<"Read tensor now "<<aux3.getDimensions()<<endl;
    // // 	//cout<<"Difference? "<<aux-aux2<<endl;exit(1);
    // // 	// Apply canonical form 
    // //   }	
    // //    }
    // And now evolve (and truncate, for test); change Dcut->0 for gilt;
    if(cnt<M){
      //state.applyMPOnoCanonical(Uops);
      state.applyMPO(Uops,tolSVD,Dcut);
    }
    cnt++;t+=delta;
  }
}

double computeEntropy(const mwArray& rdm){
  mwArray res=rdm;
  int d=rdm.getDimension(0);
  res=res+Hconjugate(res);  
  res=(ONE_c/res.trace())*res;// normalized rdm
  vector<complex_t> Dval;mwArray U; //placeholder:ignored
  wrapper::eig(res,Dval,U);
  double entropy=0.;
  for(int k=0;k<d;k++){
    double tmp=real(Dval[k]);
    if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
  }
  return entropy;
}


/** Construct the effective state for the nB block, with dimensions Dlxd^nxDr
    It starts on site k0 (default 0). It needs to include Lambda! */
void buildBlockState(mwArray& res,uMPS& state,int k0){
  int nB=state.getN();
  res=mwArray(ONE_c);
  int dres=1;int Dl=1;int Dr=1;
  vector<int> ind;
  for(int k=k0;k<nB;k++) ind.push_back(k);
  for(int k=0;k<k0;k++) ind.push_back(k);
      
  for(int p=0;p<nB;p++){
    int k=ind[p];
    int dk=state.getPhysicalDimension(k);
    int Dlk=state.getLeftVirtualDimension(k);
    int Drk=state.getRightVirtualDimension(k);
    mwArray aux=state.getA(k); // dk Dlk Drk
    aux.permute(Indices(2,1,3));aux.reshape(Indices(Dlk,dk*Drk));
    if(k==k0){
      Dl=Dlk;res=aux;
    }
    else{ res.multiplyRight(aux);}
    res.reshape(Indices(Dl*dres*dk,Drk));
    Dr=Drk;dres*=dk;
  }
  res.multiplyRight(state.getLambda(k0)); // Lambda included: nowit is the full state
  res.reshape(Indices(Dl,dres,Dr));
}

/** Construct the effective state for a block of length l, with
    dimensions Dlxd^lxDr It starts on site k0 within the nB of the
    state (default 0), but does not need to have exactly nB sites. It
    needs to include Lambda! */
void buildBlockState(mwArray& res,uMPS& state,int k0,int L){
  int nB=state.getN();
  if(L==nB){
    buildBlockState(res,state,k0);
    return;
  }
  res=mwArray(ONE_c);
  int dres=1;int Dl=1;int Dr=1;
  vector<int> ind;
  //int lastK=(k0+L-1)%nB; // pos in which I stop in the last loop
  int lambdaK=(k0+L)%nB; // pos from which to recover the Lambda
  int nrLoops=(L-(L%nB))/nB;
  int extraSites=L%nB;
  //  cout<<"L="<<L<<", nB="<<nB<<", Nr_loops="<<nrLoops<<", extra="<<extraSites<<" , last="<<lastK<<", lambda from "<<lambdaK<<endl;
  for(int p=0;p<nrLoops;p++){
    for(int k=k0;k<nB;k++) ind.push_back(k);
    for(int k=0;k<k0;k++) ind.push_back(k);
  }
  if(extraSites>0){ // continue adding
    for(int k=k0;k<min(nB,extraSites);k++) ind.push_back(k);
    for(int k=0;k<min(k0,extraSites-nB);k++) ind.push_back(k);    
  }
  //  cout<<"After filling the list: "<<ind<<endl;
  
  for(int p=0;p<L;p++){
    int k=ind[p];
     int dk=state.getPhysicalDimension(k);
    int Dlk=state.getLeftVirtualDimension(k);
    int Drk=state.getRightVirtualDimension(k);
    mwArray aux=state.getA(k); // dk Dlk Drk
    // cout<<"Tensor "<<p<<" is the "<<k<<" in uMPS: "<<aux.getDimensions()<<endl;
    aux.permute(Indices(2,1,3));aux.reshape(Indices(Dlk,dk*Drk));
    if(p==0){//(k==k0){
      Dl=Dlk;res=aux;
    }
    else{ res.multiplyRight(aux);}
    res.reshape(Indices(Dl*dres*dk,Drk));
    Dr=Drk;dres*=dk;
  }
  res.multiplyRight(state.getLambda(lambdaK)); // Lambda included: now it is the full state
  res.reshape(Indices(Dl,dres,Dr));
}



void checkBlockState(uMPS& state,int k0,ofstream* outSVD,int lB){
  mwArray psiLphR;
  if(lB==0)
    buildBlockState(psiLphR,state,0);
  else
    buildBlockState(psiLphR,state,0,lB);
  int Dl=psiLphR.getDimension(0);
  int Dr=psiLphR.getDimension(2);
  int dph=psiLphR.getDimension(1);
  psiLphR.permute(Indices(2,1,3));
  psiLphR.reshape(Indices(dph,Dl*Dr));
  mwArray U,S,Vd;  
  wrapper::svd(psiLphR,U,S,Vd);
  vector<double> sVals;S.getRealDiagonal(sVals);
  cout<<"After SVD, S="<<sVals<<endl;
  for(int p=0;p<dph;p++){
    if(p<sVals.size())
      *outSVD<<sVals[p]<<"\t";
    else
      *outSVD<<0<<"\t";
  }
  if(dph<Dl*Dr){ // need the extra dim
    mwArray projSlow(Vd);projSlow.Hconjugate();
    projSlow.multiplyRight(Vd);
    mwArray projFast=identityMatrix(Dl*Dr);
    projFast=projFast-projSlow;

    // Assuming rank of projSlow is max (=dph) => this constrains Dfast
    cout<<"Dfast should be "<<sqrt(Dl*Dr/dph)<<endl;
    
    // mwArray auxS(psiLphR);auxS.multiplyRight(projSlow);
    // mwArray auxF(psiLphR);auxF.multiplyRight(projFast);
    // cout<<"Norms: psiABC: "<<norm(reshape(psiLphR,Indices(-1,1)))
    // 	<<" projSlow: "<<norm(reshape(auxS,Indices(-1,1)))
    // 	<<" projFast: "<<norm(reshape(auxF,Indices(-1,1)))
    // 	<<endl;
    //cout<<"projSlow:"<<projSlow<<endl;
    //cout<<"projFast:"<<projFast<<endl;
    
    if(Dl>1){
      // Another test
      mwArray block(Vd); // dphys x Dl Dr
      // block.reshape(Indices(dph,Dl,Dr));
      // block.permute(Indices(2,1,3));
      // block.reshape(Indices(Dl,dph*Dr));
      // block.multiplyRight(Hconjugate(block));
      // mwArray Ul;vector<complex_t> Ll;
      // wrapper::eig(block,Ll,Ul,true);
      // cout<<"Diagonalized the block on left to give "<<Ll<<endl;
      // //cout<<"CHECK:"<<Ul*diag(Ll)*Hconjugate(Ul)-block<<endl;
            
	mwArray Ui,Si,Vdi;
	for(int i=0;i<Vd.getDimension(0);i++){
	mwArray blockI=Vd.subArray(Indices(i,-1));
	blockI.reshape(Indices(Dl,Dr));
	// if(i>0){
	//   // apply the previously found SVD
	//   blockI.multiplyLeft(Hconjugate(Ui));
	//   blockI.multiplyRight(Hconjugate(Vdi));
	//   cout<<"Vector "<<i<<" now "<<blockI<<endl;
	// }
	
	// Now try its SVD
	wrapper::svd(blockI,Ui,Si,Vdi);
	vector<double> sValsi;Si.getRealDiagonal(sValsi);
	cout<<"After SVD of vector "<<i<<", Si="<<sValsi<<endl;
	//cout<<"while Ui="<<Ui<<", Vi="<<Vi<<endl;
	
	// for(int j=0;j<sVals.size();j++){
	//   mwArray blockJ=Vd.subArray(Indices(j,-1));
	//   blockJ.reshape(Indices(Dl,Dr));
	//   blockJ.Hconjugate();
	//   blockJ.multiplyLeft(blockI); // Dl*Dl
	//   // Now check what happens when multilpying the Ul
	//   blockJ.multiplyLeft(Hconjugate(Ul));
	//   blockJ.multiplyRight(Ul);
	//   cout<<"For component ("<<i<<","<<j<<")->"<<blockJ<<endl;
	//   if(i==j&&i==0){
	//     wrapper::eig(block,Ll,Ul,true);
	//     cout<<"Diagonalized the block "<<i<<","<<j<<" on left to give "<<Ll<<endl;
	//     //cout<<"CHECK:"<<Ul*diag(Ll)*Hconjugate(Ul)-block<<endl
	//   }

	// }
      }

    }
  }
}

 double computeNegativity(mwArray rho,int d1,int d2){  
  rho.reshape(Indices(d1,d2,d1,d2));
  rho.permute(Indices(1,4,3,2));
  rho.reshape(Indices(d1*d2,d1*d2));
  mwArray U,S,Vd; 
  wrapper::svd(rho,U,S,Vd);
  return real(S.trace());
}

void negativitiesBlockState(mwArray& psiLphR,ofstream* outS){
  int Dl=psiLphR.getDimension(0);
  int Dr=psiLphR.getDimension(2);
  int dph=psiLphR.getDimension(1);
  mwArray rho_phR(psiLphR);  
  rho_phR.reshape(Indices(Dl,dph*Dr));
  rho_phR.permute(Indices(2,1));
  rho_phR.multiplyRight(Hconjugate(rho_phR)); // dph*Dr,dph*Dr
  double NphR=computeNegativity(rho_phR,dph,Dr);

  mwArray rho_phL(psiLphR);  
  rho_phL.reshape(Indices(Dl*dph,Dr));
  rho_phL.multiplyRight(Hconjugate(rho_phL)); // Dl*dph,Dl*dph
  double NphL=computeNegativity(rho_phL,Dl,dph);

  mwArray rho_LR(psiLphR);  
  rho_LR.permute(Indices(1,3,2));
  rho_LR.reshape(Indices(Dl*Dr,dph));
  rho_LR.multiplyRight(Hconjugate(rho_LR));
  double NLR=computeNegativity(rho_LR,Dl,Dr);



  *outS<<NphL<<"\t"<<NphR<<"\t"<<NLR<<"\t";
}

void entropiesBlockState(uMPS& state,int k0,ofstream* outS,int lB){
  mwArray psiLphR;
  buildBlockState(psiLphR,state,0,lB);
  int Dl=psiLphR.getDimension(0);
  int Dr=psiLphR.getDimension(2);
  int dph=psiLphR.getDimension(1);
  //  cout<<"Constructed the block state size: "<<psiLphR.getDimensions()
  //  <<", norm "<<norm(reshape(psiLphR,Indices(Dl*dph*Dr)))<<endl;

  double Sph(0.),SL(0.),SR(0.);
  
  psiLphR.reshape(Indices(Dl,dph*Dr));
  mwArray U,S,Vd;  
  wrapper::svd(psiLphR,U,S,Vd);
  SL=computeEntropy(S*S);
  double NL=real(S.trace()); NL*=NL;

  
  psiLphR.reshape(Indices(Dl,dph,Dr));
  psiLphR.permute(Indices(2,1,3));
  psiLphR.reshape(Indices(dph,Dl*Dr));
  wrapper::svd(psiLphR,U,S,Vd);
  Sph=computeEntropy(S*S);
  double Nph=real(S.trace()); Nph*=Nph;
  //  vector<double> sValsPh;S.getRealDiagonal(sValsPh);

  psiLphR.reshape(Indices(dph*Dl,Dr));
  wrapper::svd(psiLphR,U,S,Vd);
  SR=computeEntropy(S*S);
  double NR=real(S.trace()); NR*=NR;
  *outS<<SL<<"\t"<<Sph<<"\t"<<SR<<"\t";

  psiLphR.reshape(Indices(dph,Dl,Dr));
  psiLphR.permute(Indices(2,1,3)); // back to original
  negativitiesBlockState(psiLphR,outS);
  *outS<<NL<<"\t"<<Nph<<"\t"<<NR<<"\t";
}

void negativitiesBlockState(uMPS& state,int k0,ofstream* outS,int lB){
  mwArray psiLphR;
  buildBlockState(psiLphR,state,0,lB);
  negativitiesBlockState(psiLphR,outS);
}


int applyTransferOperator(mwArray& res,uMPS& state,int k1,int L,bool withLambda=0){
  int dphys=res.getDimension(0); // dphys(dxd'), Dr' Dr
  int Dr=res.getDimension(1);
  int nB=state.getN();
  if(k1>=nB){ cout<<"Invalid initial point "<<k1<<endl;exit(1);}
  // list of sites to add
  vector<int> ind;
  int nrLoops=(L-(L%nB))/nB;
  int extraSites=L%nB;
  //    cout<<"L="<<L<<", nB="<<nB<<", Nr_loops="<<nrLoops<<", extra="<<extraSites<<endl;
  for(int p=0;p<nrLoops;p++){
    for(int k=k1;k<nB;k++) ind.push_back(k);
    for(int k=0;k<k1;k++) ind.push_back(k);
  }
  if(extraSites>0){ // continue adding
    for(int k=k1;k<min(nB,extraSites);k++) ind.push_back(k);
    for(int k=0;k<min(k1,extraSites-nB);k++) ind.push_back(k);    
  }
  //cout<<"applyTransferOp from k1="<<k1<<", L="<<L<<" ind:"<<ind<<endl;
  // cout<<"received res:"<<res.getDimensions()<<endl;
  for(int p=0;p<L;p++){
    int k=ind[p];
    int dk=state.getPhysicalDimension(k);
    int Dlk=state.getLeftVirtualDimension(k);
    int Drk=state.getRightVirtualDimension(k);
    mwArray aux=state.getA(k); // dk Dlk Drk
    if(withLambda&&p==L-1){
      aux.reshape(Indices(dk*Dlk,Drk));
      aux.multiplyRight(state.getLambda((k+1)%nB));
      aux.reshape(Indices(dk,Dlk,Drk));
    }
    //cout<<"Tensor "<<p<<" is the "<<k<<" in uMPS: "<<aux.getDimensions()<<endl;
    aux.permute(Indices(2,1,3));aux.reshape(Indices(Dlk,dk*Drk));
    res.reshape(Indices(dphys*Dr,Dr));
    res.multiplyRight(aux); // dphys*Dr' x dk*Drk        
    res.reshape(Indices(dphys,Dr*dk,Drk));
    res.permute(Indices(1,3,2));res.reshape(Indices(dphys*Drk,Dr*dk)); // dphys*Dr, Dr'*dk'
    aux.reshape(Indices(Dlk*dk,Drk));
    aux.conjugate();
    res.multiplyRight(aux); //dphys*Drk,Drk'
    res.reshape(Indices(dphys,Drk,Drk));
    res.permute(Indices(1,3,2));
    Dr=Drk;
  }
  return ind.back();
}


void distantPair(uMPS& state,int k0,int lB,ofstream* outS){
  // RM for sites k0 and
  int nB=state.getN();
  if(k0>=nB){ cout<<"Invalid initial point "<<k0<<endl;exit(1);}
  int d=state.getPhysicalDimension(k0);
  int Dl=state.getLeftVirtualDimension(k0);
  int Dr=state.getRightVirtualDimension(k0);
  mwArray aux=state.getA(k0); // d Dl Dr
  // cout<<"aux:"<<aux<<endl;
  // cout<<"Lambda:"<<state.getLambda(k0)<<endl;
  aux.permute(Indices(2,1,3));aux.reshape(Indices(Dl,d*Dr));
  aux.multiplyLeft(Hconjugate(aux)); // d'*Dr', d*Dr
  aux.reshape(Indices(d,Dr,d,Dr));
  aux.permute(Indices(3,1,2,4)); // d, d', Dr', Dr
  int dphys=d*d;
  aux.reshape(Indices(dphys,Dr,Dr)); // dphysL Dr' Dr
  int lastK=lB==0?k0:applyTransferOperator(aux,state,(k0+1)%nB,lB);
  //cout<<"After applyTransferOp, aux:"<<aux.getDimensions()<<" last: "<<lastK<<endl;
  Dr=aux.getDimension(1); // now aux dphys x Dr' x Dr
  aux.reshape(Indices(dphys,Dr*Dr)); // dphys Dr' Dr
  int k=(lastK+1)%nB;
  int dk=state.getPhysicalDimension(k);
  int Dlk=state.getLeftVirtualDimension(k);
  int Drk=state.getRightVirtualDimension(k);
  mwArray auxA=state.getA(k); // dk Dlk Drk
  auxA.reshape(Indices(dk*Dlk,Drk));
  auxA.multiplyRight(state.getLambda((lastK+2)%nB));
  //cout<<auxA<<endl;
  auxA.multiplyRight(Hconjugate(auxA)); // dk*Dlk,dk'*Dlk'
  auxA.reshape(Indices(dk,Dlk,dk,Dlk));
  auxA.permute(Indices(4,2,1,3));
  auxA.reshape(Indices(Dlk*Dlk,dk*dk));
  aux.multiplyRight(auxA);
  aux.reshape(Indices(d,d,dk,dk));
  aux.permute(Indices(1,3,2,4));
  aux.reshape(Indices(d*dk,d*dk));
  // Get negativity
  double neg=computeNegativity(aux,d,dk);  
  *outS<<neg<<"\t";
}

void buildIndexList(vector<int>& ind,int k0,int L,int nB){
  ind.clear();
  if(k0>=nB){ cout<<"Invalid initial point "<<k0<<endl;exit(1);}
  int nrLoops=(L-(L%nB))/nB;
  int extraSites=L%nB;
  //  cout<<"L="<<L<<", nB="<<nB<<", Nr_loops="<<nrLoops<<", extra="<<extraSites<<" , last="<<lastK<<", lambda from "<<lambdaK<<endl;
  for(int p=0;p<nrLoops;p++){
    for(int k=k0;k<nB;k++) ind.push_back(k);
    for(int k=0;k<k0;k++) ind.push_back(k);
  }
  if(extraSites>0){ // continue adding
    for(int k=k0;k<min(nB,extraSites);k++) ind.push_back(k);
    for(int k=0;k<min(k0,extraSites-nB);k++) ind.push_back(k);    
  }
}

void distantBlock(uMPS& state,int k0,int nS,int lB,ofstream* outS){
  // RM for sites k0 and
  int nB=state.getN();
  vector<int> ind;
  buildIndexList(ind,k0,nS,nB);
  int dphys=1;int Dr=state.getLeftVirtualDimension(ind[0]);
  mwArray rdm=identityMatrix(Dr);
  rdm.reshape(Indices(dphys*dphys,Dr,Dr)); //d^l x d^l' x Dr' x Dr 
  for(int p=0;p<ind.size();p++){
    int pos=ind[p];
    int dp=state.getPhysicalDimension(pos);
    int Dlp=state.getLeftVirtualDimension(pos);
    int Drp=state.getRightVirtualDimension(pos);
    mwArray aux=state.getA(pos); // d Dl Dr
    aux.permute(Indices(2,1,3));aux.reshape(Indices(Dlp,dp*Drp));
    rdm.reshape(Indices(dphys*dphys*Dr,Dr));
    rdm.multiplyRight(aux); // dphys^2 *Dr' , dp*Drp
    rdm.reshape(Indices(dphys*dphys,Dr,dp*Drp));
    rdm.permute(Indices(2,1,3));
    rdm.reshape(Indices(Dr,dphys*dphys*dp*Drp));
    rdm.multiplyLeft(Hconjugate(aux)); // dp*Drp', dphys*dphys * dp*Drp
    rdm.reshape(Indices(dp,Drp,dphys,dphys,dp,Drp));
    rdm.permute(Indices(3,5,4,1,2,6));
    dphys*=dp;Dr=Drp;
    rdm.reshape(Indices(dphys*dphys,Dr,Dr));
  }
  int lastK=lB==0?ind.back():applyTransferOperator(rdm,state,(ind.back()+1)%nB,lB);
  Dr=rdm.getDimension(1); // now aux dphys^2 x Dr' x Dr
  rdm.reshape(Indices(dphys*dphys,Dr*Dr)); // dphys^2 Dr' Dr
  int k=(lastK+1)%nB;
  ind.clear();  buildIndexList(ind,k,nS,nB);
  int dphysR=1;int DlR=state.getRightVirtualDimension(ind.back());
  mwArray rdmR=identityMatrix(DlR);
  rdmR.reshape(Indices(DlR,DlR,dphysR*dphysR)); // Dl' x Dl x d^l*d^l'
  for(int p=ind.size()-1;p>=0;p--){
    int pos=ind[p];
    int dp=state.getPhysicalDimension(pos);
    int Dlp=state.getLeftVirtualDimension(pos);
    int Drp=state.getRightVirtualDimension(pos);
    mwArray auxA=state.getA(pos); // dk Dlk Drk
    auxA.reshape(Indices(dp*Dlp,Drp));
    if(p==ind.size()-1){
      auxA.multiplyRight(state.getLambda((pos+1)%nB));
    }
    rdmR.reshape(Indices(DlR,DlR*dphysR*dphysR));
    rdmR.multiplyLeft(conjugate(auxA)); //dp'*Dlp' Dl  dphys^2
    rdmR.reshape(Indices(dp*Dlp,DlR,dphysR*dphysR));
    rdmR.permute(Indices(2,1,3)); // DlR dp*Dlp' dphys^2;
    rdmR.reshape(Indices(DlR,dp*Dlp*dphysR*dphysR));
    rdmR.multiplyLeft(auxA); // dp*Dlp dp'*Dlp' dphys^2
    rdmR.reshape(Indices(dp,Dlp,dp,Dlp,dphysR,dphysR));
    rdmR.permute(Indices(4,2,1,5,3,6));
    dphysR*=dp;DlR=Dlp;
    rdmR.reshape(Indices(DlR,DlR,dphysR*dphysR));
  }
  rdmR.reshape(Indices(DlR*DlR,dphysR*dphysR));
  rdm.multiplyRight(rdmR); // dphys*dphys x dphysR*dphysR
  rdm.reshape(Indices(dphys,dphys,dphysR,dphysR));
  rdm.permute(Indices(1,3,2,4));
  rdm.reshape(Indices(dphys*dphysR,dphys*dphysR));
  // putForMatlab(cout,rdm,"rdm");
  // Get negativity
  double neg=computeNegativity(rdm,dphys,dphysR);  
  *outS<<neg<<"\t";  
}

void channelBlock(uMPS& state,int k0,ofstream* outS,int lB){
  int nB=state.getN();
  int dphys=1;int Dr=state.getLeftVirtualDimension(k0);
  mwArray channel=identityMatrix(Dr*Dr);
  // trick
  channel.reshape(Indices(Dr*Dr,Dr,Dr));
  int lastK=applyTransferOperator(channel,state,k0,lB,true);
  int Dr2=channel.getDimension(1); // now aux Dr^2 x Dr2' x Dr2
  channel.reshape(Indices(Dr*Dr,Dr2*Dr2));
  //Problem: not Hermitian! Would need sth else to figure out spectrum
   mwArray U;vector<complex_t> Lamb;
   wrapper::eigNH(channel,Lamb,U,true);
  for(int p=0;p<Dr*Dr;p++){
     if(p<Lamb.size()) *outS<<real(Lamb[p])<<"\t"<<imag(Lamb[p])<<"\t";
     else *outS<<"0\t 0\t";   
   }
  if(Dr*Dr<10) putForMatlab(cout,channel,"T");
}

// void distantPairs(uMPS& state,int k0,int lBmin,int lBmax,ofstream* outS){
//   // RM for sites k0 and
//   int nB=state.getN();
//   if(k0>=nB){ cout<<"Invalid initial point "<<k0<<endl;exit(1);}
//   int d=state.getPhysicalDimension(k0);
//   int Dl=state.getLeftVirtualDimension(k0);
//   int Dr=state.getRightVirtualDimension(k0);
//   mwArray aux=state.getA(k0); // d Dl Dr
//   // cout<<"aux:"<<aux<<endl;
//   // cout<<"Lambda:"<<state.getLambda(k0)<<endl;
//   aux.permute(Indices(2,1,3));aux.reshape(Indices(Dl,d*Dr));
//   aux.multiplyLeft(Hconjugate(aux)); // d'*Dr', d*Dr
//   aux.reshape(Indices(d,Dr,d,Dr));
//   aux.permute(Indices(3,1,2,4)); // d, d', Dr', Dr
//   int dphys=d*d;
//   aux.reshape(Indices(dphys,Dr,Dr)); // dphys Dr' Dr
//   int lB=lBmin;
//   while(lB<=lBmax){
//     int toApply=(lB==lBmin)?lB:1;
//     int lastK=applyTransferOperator(aux,state,(k0+1)%nB,toApply);
//     //cout<<"After applyTransferOp, aux:"<<aux.getDimensions()<<" last: "<<lastK<<endl;
//     Dr=aux.getDimension(1); // now aux dphys x Dr' x Dr
//     aux.reshape(Indices(dphys,Dr*Dr)); // dphys Dr' Dr
//   int k=(lastK+1)%nB;
//   int dk=state.getPhysicalDimension(k);
//   int Dlk=state.getLeftVirtualDimension(k);
//   int Drk=state.getRightVirtualDimension(k);
//   mwArray auxA=state.getA(k); // dk Dlk Drk
//   auxA.reshape(Indices(dk*Dlk,Drk));
//   auxA.multiplyRight(state.getLambda((lastK+2)%nB));
//   //cout<<auxA<<endl;
//   auxA.multiplyRight(Hconjugate(auxA)); // dk*Dlk,dk'*Dlk'
//   auxA.reshape(Indices(dk,Dlk,dk,Dlk));
//   auxA.permute(Indices(4,2,1,3));
//   auxA.reshape(Indices(Dlk*Dlk,dk*dk));
//   aux.multiplyRight(auxA);
//   aux.reshape(Indices(d,d,dk,dk));
//   aux.permute(Indices(1,3,2,4));
//   aux.reshape(Indices(d*dk,d*dk));




//   // auxA.permute(Indices(2,1,3));auxA.reshape(Indices(Dlk,dk*Drk));
//   // aux.reshape(Indices(dphys*Dr,Dr));
//   // aux.multiplyRight(auxA);
//   // aux.reshape(Indices(dphys,Dr,dk,Drk));
//   // aux.permute(Indices(1,3,2,4));
//   // aux.reshape(Indices(dphys*dk,Dr*Drk));
//   // auxA.reshape(Indices(Dlk,dk,Drk));
//   // auxA.permute(Indices(1,3,2),true); // already conj
//   // auxA.reshape(Indices(Dlk*Drk,dk));
//   // aux.multiplyRight(auxA); //dphys*dk,dk'
//   // aux.reshape(Indices(d,d,dk,dk));
//   // aux.permute(Indices(1,3,2,4));
//   // aux.reshape(Indices(d*dk,d*dk));
//   //  cout<<"rdm:"<<aux<<endl;

//   // Get negativity
//   double neg=computeNegativity(aux,d,dk);  

//   *outS<<neg<<"\t";
// }




// #include "quicksort.h"

// double factorVector(const vector<double>& V_,int D1,int D2,double tol,vector<double>& v1,vector<double>& v2){
//   double error=-1;
//   vector<double> V(V_);
//   // sort the vector
//   quicksort(V,0,V.size()-1);
//   if(D1*D2>V.size()){
//     cout<<"Adding zeros to the vector"<<endl;
//     while(V.size()<D1*D2) V.push_back(0);    
//   }
//   if(D1*D2<V.size()){
//     double firstExtra=V[D1*D2];
//     cout<<"WARNING: vector too long, check if discarding is possible: largest remaining value "<<firstExtra<<endl;
//     if(abs(firstExtra)>tol){ // make relative
//       cout<<"Discarding"<<endl;
//     }
//   }

  
// }





