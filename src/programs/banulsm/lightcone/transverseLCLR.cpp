#include <math.h>
#include <iomanip>
#include <sstream>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "Properties.h"
#include "time.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

//#define SYMM 1 // If SYMM=1, only left part is computed explicitly, and right is just the reflection

#include "IsingHamiltonian.h"

using namespace shrt;

/** Try the folded transverse contraction strategy but instead of
    looking for the fixed point of the column, iterate as time
    increases using the light cone structure. (31.10.2021)

    Add the possibility to renormalize the old times, to maintain a
    constant length of the MPS and control the resources needed in a
    better way (allow for larger D). (2.11.2021)

    This version allows a LC with different slope, to try to use LR 
    velocity (23.12.2021)

    PROBLEM: No tmp storage implemented here!

    (17.1.2022) Computing the observables at earlier times

  */

int d=2;
double tolSVD=1E-12;      

void constructInitialProductTensors(int d,int D,int initSt,mwArray& W,mwArray& leftB);
void saveObservables(const mwArray& rdm2,int cnt,double delta,const string& outfname,
		     double J,double g,double h);
void computeAndSaveObservables(const MPS& mpsB,const Operator& middleUL,const Operator& rightmostOpL,
			       const mwArray& edgeUL, const mwArray& edgeUR,
			       int cnt,double delta,const string& outfname,
			       double J,double g,double h);


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
  double v_=props.getDoubleProperty("vLR");
  int initSt=props.getIntProperty("initSt");
  double delta=props.getDoubleProperty("delta");
  int M=props.getIntProperty("M"); // max nr of time steps to be simulated
  //double alpha=props.getDoubleProperty("alpha");
  //if(alpha<0) alpha=1.;
  int rate=props.getIntProperty("rate");
  if(rate<=0) rate=1; // default frequency
  int Dcut=props.getIntProperty("D"); // the truncation  
  //  int Dcontr=props.getIntProperty("Drho"); // the dimension used for the effective rho+boundary as it evolves
  string outfname=props.getProperty("outputfile");
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
  bool app=props.getIntProperty("append")!=0; // default is appending

  int deepF=props.getIntProperty("bufferLR"); // how deep we go to avoid the edge of the LC
  if(deepF<0) deepF=2; // default is leaving one site on each side
  
  // Preliminaries
  int d2=d*d;

 // The operator I want to compute on L sites (for now, just a sigma_x in the center)
  // MPS mpsOp(L+2,1,d2);
  // MPS mpsTr(L+2,1,d2); // Another one for the trace
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d2,1,1));
  // for(int k=1;k<L+1;k++){
  //   mpsOp.setA(k,id2);
  //   mpsTr.setA(k,id2);
  // }
  // and the middle op I set to sigmax
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1*I_c,ZERO_c};
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c}; //Z
  mwArray sig0=identityMatrix(d);
  mwArray sigX=mwArray(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  mwArray sigZ(Indices(d,d),dataZ);

  // Decide proportion of time steps vs space ones: how many t steps
  // with one of x (1 to 1 is the usual program)
  int nT=ceil(1./(v_*delta));
  cout<<"I should be doing "<<nT<<" per space site"<<endl;
  
  
  // Prepare the initial tensors (initial state, boundary, edges, system)

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;

  // Basic terms for the exponential: I assume TI case, infinite
  // chain, but I only need the central (TI MPO) and the edges for a
  // finite case, which appear when the light-cone ends.
  MPO Uevol(3); // this is the double one
  {
    IsingHamiltonian hamH(3,d,J_,g_,h_);
    //hamH.registerTimeDependence();
    MPO _Uevol_1(3); // held just temprarily
    hamH.getUMPO(_Uevol_1,-delta*I_c,0); // the single layer normal evolution
    // But I need to double it
    doubleMPO(_Uevol_1,Uevol,true);
    cout<<"Created the double Uevol"<<endl;
  }

  // The most repeated tensor in the network (doubled U)  
  mwArray midU(Uevol.getOp(1).getFullData()); // d2 x xi x d2 x xi
  //  putForMatlab(cout,midU,"midU");
  int xi=midU.getDimension(1);
  //mwArray midUl(midU);midUl.permute(Indices(1,4,2,3));midUl.reshape(Indices(d2*xi,xi*d2));
  //  midU.reshape(Indices(d2*xi,d2*xi)); // to contract with right boundary

  // Boundary (without tracing)
  mwArray boundaryUL=Uevol.getOp(0).getFullData(); //d2 1 d2 xi
  mwArray boundaryUR=Uevol.getOp(2).getFullData(); //d2 xi d2 1

  // Initial edges
  mwArray edgeUL=Uevol.getOp(0).getFullData(); //d2 1 d2 xi
  //cout<<"Recovered edgeUL from Uevol: "<<edgeUL<<endl;
  mwArray edgeUR=Uevol.getOp(2).getFullData(); //d2 xi d2 1
  //cout<<"Recovered edgeUR from Uevol: "<<edgeUR<<endl;
  // For the edges, I only need the contraction with id2 above
  id2.reshape(Indices(1,d2));
  edgeUL.reshape(Indices(d2,d2*xi));edgeUL.multiplyLeft(id2);edgeUL.reshape(Indices(d2,xi)); //d2 x xi2
  //edgeUL.transpose(); // xi2 x d2
  edgeUR.reshape(Indices(d2,xi*d2));edgeUR.multiplyLeft(id2);edgeUR.reshape(Indices(xi,d2)); // xi2 x d2

  mwArray leftB,W;
  int D0=1; // initial bond dim
  constructInitialProductTensors(d,D0,initSt,W,leftB);
  // Now leftB is (1xD0^2) and W (d2xD0^2xD0^2) 

  MPS mpsB(2,d2,D0*D0);
  {
    mpsB.replaceSite(0,reshape(identityMatrix(d2),Indices(d2,1,d2)),false);
    mwArray aux(W); // d2xD0xD0
    aux.permute(Indices(2,1,3));
    aux.reshape(Indices(D0,d2*D0));
    aux.multiplyLeft(leftB); // 1 x d2*D0;
    mpsB.replaceSite(1,permute(reshape(aux,Indices(d2,D0,1)),Indices(2,1,3)),false);
    cout<<"Initialized boundary MPS: "<<mpsB<<endl;
  }
  int lenB=mpsB.getLength(); // 1+ nr of steps in the next step

  // Common pieces for the evolution operator are the edge tensor(s) and the state one

  Operator leftmostOpL(reshape(identityMatrix(d2),Indices(d2,1,1,d2)));
  Operator secondOpL(reshape(edgeUL*reshape(permute(midU,Indices(2,1,3,4)),
					    Indices(xi,d2*d2*xi)), // d2 x d2*d2*xi
			     Indices(d2,d2,d2,xi)),Indices(4,2,1,3)); 
  Operator rightmostOpL(reshape(W,Indices(d2,D0*D0,D0*D0,1)),Indices(3,1,2,4));
  Operator middleUL(midU,Indices(4,1,2,3)); // xi(r) x d2(u) x xi x d2
  // if more than one step each time, I also need a reshaped boundary
  Operator extraUL(boundaryUL,Indices(4,1,2,3)); // xi x d2 x 1 x d2
  
  Operator upperU(reshape(id2*reshape(midU,Indices(d2,xi*d2*xi)),Indices(1,xi,d2,xi)),Indices(4,1,2,3)); // tracing out te uppermost ones!


  int cnt=0;double time=0.;

  // For the first stp, the rdm is special
  mwArray rdm0(W);// d2xD02xD02
  rdm0.permute(Indices(2,1,3));rdm0.reshape(Indices(D0*D0,d2*D0*D0));
  rdm0.multiplyLeft(leftB);rdm0.reshape(Indices(d2,D0*D0));
  rdm0.multiplyRight(permute(rdm0,Indices(2,1))); // d2xd2
  rdm0.permute(Indices(2,1));rdm0.reshape(Indices(d2*d2,1)); // as will appear later, with index of rightmost index first
    
  saveObservables(rdm0,cnt,delta,outfname,J_,g_,h_);

  //cout<<"At the beginning cnt="<<cnt<<" L="<<mpsB.getLength()<<endl;;
  while(cnt<M){
  // step by step:
  // 1) evolve the boundary (and grow it) 
  // 2) increase the MPO(s) for the system
  // 3) contract expectation values

    // how many time steps to add now?
    int toadd=1; // first I add a time step
    lenB=mpsB.getLength(); // 1+ nr of steps in the next step
    int Lt=lenB+toadd;
    MPO evol(Lt);
    // initialize the MPS
    MPS aux(Lt,Dcut,xi);
    for(int k=0;k<Lt;k++){
      if(k<lenB)
	aux.replaceSite(Lt-k-1,mpsB.getA(lenB-k-1),false);
      else
	aux.replaceSite(Lt-k-1,reshape(mwArray(ONE_c),Indices(1,1,1)),false);
    }
    //cout<<"Step"<<cnt<<" Extended the MPS "<<mpsB<<" to "<<aux<<endl;
    evol.setOp(0,&leftmostOpL,false);
    evol.setOp(1,&secondOpL,false); // which should be at pos lenB from R
    evol.setOp(Lt-1,&rightmostOpL,false);
    for(int k=2;k<Lt-1;k++)
      evol.setOp(k,&middleUL,false);

    //cout<<"Prepared evolution MPO: "<<evol<<endl;
    mpsB.approximate(evol,aux,Dcut); // first estimate by SVDs
    contractor.optimize(evol,aux,mpsB,Dcut);
    //    cout<<"After optimize,mpsB:"<<mpsB<<endl;
    //    MPS mpsIso(toadd+2,Dcut,xi);

    // //mpsB.gaugeCond('R',1);
    //    mpsB.gaugeCond('R',1);
    cnt+=toadd;
    time+=toadd*delta;
    //cout<<"Applied one step now L="<<mpsB.getLength()<<", cnt="<<cnt<<endl;

    // and then stretch the MPS with the nT-1 extra sites (one already added)
    if(nT>1){
      mwArray extraL(boundaryUL); // d2 1 d2 xi
      extraL.permute(Indices(4,1,3,2)); // xi d2(u) d2(d) 1
      extraL.reshape(Indices(xi*d2,d2)); // ready to multiply
      mwArray leftId=identityMatrix(d2);leftId.reshape(Indices(d2,1,d2)); // 
      mwArray auxA=mpsB.getA(0).getA(); // d2 x 1 x d2
      auxA.reshape(Indices(d2,d2));
      auxA.multiplyLeft(extraL); //xi*d2(u) x d2
      auxA.reshape(Indices(xi,d2,d2));
      mpsB.replaceSite(0,auxA,false);
      extraL.reshape(Indices(xi,d2,d2)); // ready to use as MPS tensor
      for(int p=1;p<nT-1;p++){
	mpsB.insertSite(0,leftId,false,false);
	mpsB.replaceSite(0,extraL,false); // no checkdim, no norm
      }
      mpsB.insertSite(0,leftId,false,false);
      cnt+=(nT-1);
      time+=delta*(nT-1);
      //cout<<"Inserted nT-1 now L="<<mpsB.getLength()<<", cnt="<<cnt<<endl;
    }
    mpsB.gaugeCond('R',1);

    // now compute observables inside (would be sightly better to do before stretching)
    //int deepF=4; // how deep
    if(cnt>=deepF*nT){
      Lt=mpsB.getLength();
      MPO midMPO(Lt);midMPO.setOp(Lt-1,&rightmostOpL,false);
      for(int k=1;k<Lt-1;k++)midMPO.setOp(k,&middleUL,false);
      midMPO.setOp(0,&upperU,false);
      // First one special: I trace out phys indices
      MPO mpoCenter(Lt);
      const MPO* ptrs[]={&midMPO,&midMPO};
      MPO::join(2,ptrs,mpoCenter);

      aux=mpsB;aux.applyLocalOperator(0,permute(edgeUL,Indices(2,1)));
      // And (only because this case is symmetric) copy the left boundary, but complex, to use as bra
      MPS bra(aux);bra.conjugateMPS();

      for(int p=deepF*nT-1;p>=(deepF-1)*nT;p--){
	int cntP=cnt-p; // time to which this contraction corresponds
	mwArray auxR=contractor.contract(aux,mpoCenter,bra,p,'R'); // D x d2(r)*d2(l) x D
	mwArray auxL=contractor.contract(aux,mpoCenter,bra,p+1,'L'); // D x d2(r)*d2(l) x D
	int Dp=aux.getA(p).getDr();
	auxR.reshape(Indices(Dp,d,d,d,d,Dp));
	auxR.permute(Indices(2,4,1,3,5,6));
	auxR.reshape(Indices(d*d,Dp*d*d*Dp));
	auxL.reshape(Indices(Dp,d,d,d,d,Dp));
	auxL.permute(Indices(1,3,5,6,2,4));
	auxL.reshape(Indices(Dp*d*d*Dp,d*d));
	mwArray rdm2=auxR*auxL;
	rdm2.reshape(Indices(d,d,d,d));rdm2.permute(Indices(1,3,2,4));
	rdm2.reshape(Indices(d2*d2,1));
	saveObservables(rdm2,cntP,delta,outfname,J_,g_,h_);
	//cout<<"For contraction (R) up to (excluded) "<<p<<", cnt="<<cntP<<endl;
      }
    }    
	
    // if(nT>1&&cnt>=21){
    //   // extend the boundary by nT sites
    //   mwArray extraL(boundaryUL); // d2 1 d2 xi
    //   extraL.permute(Indices(4,1,3,2)); // xi d2(u) d2(d) 1
    //   extraL.reshape(Indices(xi*d2,d2)); // ready to multiply
    //   mwArray leftId=identityMatrix(d2);leftId.reshape(Indices(d2,1,d2)); // 
    //   for(int p=0;p<nT-1;p++){
    // 	// 1) insert an extra site and change the (up to now) first site of the mps
    // 	mwArray auxA=mpsB.getA(0).getA(); // d2 x 1 x d2
    // 	auxA.reshape(Indices(d2,d2));
    // 	auxA.multiplyLeft(extraL); //xi*d2(u) x d2
    // 	auxA.reshape(Indices(xi,d2,d2));
    // 	mpsB.insertSite(0,leftId,false,false); // no checkdim, no norm
    // 	mpsB.replaceSite(1,auxA,false); // no check dim
    // 	// do the same on the bra (for simplicity, I just copy and
    // 	// gauge, as before, but could be made more efficient keeping
    // 	// // // // the rdm part from before and from one step to the next)
    // 	// // // bra=mpsB;bra.conjugateMPS();
    // 	// // // //bra.gaugeCond('R',1);
    // 	// // // // try gauge again? 
    // 	// // // MPS::gaugeCond(mpsB,bra,'L');
    // 	// // // mwArray rdm2=contractor.contract(mpsB,bra,0,'R'); // d2 x 1 x d2
    // 	// // // rdm2.reshape(Indices(d2,d2));
    // 	// // // //mwArray auxRDM(rdm2); // so that I can recompute with various time steps
    // 	// // // rdm2.multiplyLeft(reshape(conjugate(bra.getA(0).getA()),Indices(d2,d2)));
    // 	// // // rdm2.multiplyRight(permute(reshape(mpsB.getA(0).getA(),Indices(d2,d2)),Indices(2,1)));
    // 	// // // //cout<<"Here rdm2="<<rdm2<<endl;if(cnt>10)exit(1);
    // 	// // // rdm2.reshape(Indices(d2*d2,1));
    // 	cnt++;
    // 	time+=delta;

    // 	if(1){
    // 	  // // To compute expectation values, prepare the middle MPO
    // 	  Lt=mpsB.getLength();
    // 	  MPO midMPO(Lt);midMPO.setOp(Lt-1,&rightmostOpL,false);
    // 	  for(int k=0;k<Lt-1;k++)midMPO.setOp(k,&middleUL,false);
    // 	  MPO mpoCenter(Lt);
    // 	  const MPO* ptrs[]={&midMPO,&midMPO};
    // 	  MPO::join(2,ptrs,mpoCenter);
    
    // 	  aux=mpsB;aux.applyLocalOperator(0,permute(edgeUL,Indices(2,1)));
    // 	  // And (only because this case is symmetric) copy the left boundary, but complex, to use as bra
    // 	  MPS bra(aux);bra.conjugateMPS();
    // 	  mwArray rdm2=contractor.contract(aux,mpoCenter,bra,-1,'R'); // 1 x d2*d2 x 1
    // 	  rdm2.reshape(Indices(d2*d2,1));
	  
    // 	  saveObservables(rdm2,cnt+1,delta,outfname,J_,g_,h_);
    // 	}
    //   }      
    // }    
  }

}

void constructInitialProductTensors(int d,int D,int initSt,mwArray& W,mwArray& leftB){
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sig0=identityMatrix(d);
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray sigX(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  
  switch(initSt){
  case 1: // X+
    {
      W=.5*(sig0+sigX);
      break;
    }
  case 2: // Y+
    {
      W=.5*(sig0+sigY);
      break;
    }
  case 3: // Z+
    {
      W=.5*(sig0+sigZ);
      break;
    }
  case 4: // mixed state at infinite temperature (Id)
    {
      W=.5*(sig0);
      break;
    }
  default:
    cout<<"Error: unknown type of intSt="<<initSt<<endl;
    exit(1);
  }
  W.reshape(Indices(d*d,1,1));
  W.resize(Indices(d*d,D,D));
  leftB=mwArray(ONE_c);leftB.reshape(Indices(1,1));leftB.resize(Indices(1,D));
}

/* ***************************************************************** */
/* A direct product of single site operators, returned as a  vector  */
/* to be used with rdm (indices ud ud)                               */
/* ***************************************************************** */

void prepareOperatorProduct(mwArray& result,const mwArray& opA,const mwArray& opB){
  result=opA; 
  result.reshape(Indices(d*d,1));
  result.multiplyRight(reshape(opB,Indices(1,d*d)));
  result.reshape(Indices(d,d,d,d)); // order here: ud ud
  result.reshape(Indices(1,d*d*d*d));
}

void saveObservables(const mwArray& rdm2,int cnt,double delta,const string& outfname,double J_,double g_,double h_){
  static bool first=true;
  // Operators needed (only local)
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1*I_c,ZERO_c};
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c}; //Z
  mwArray sig0=identityMatrix(d);
  mwArray sigX=mwArray(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  mwArray sigZ(Indices(d,d),dataZ);
  int d2=d*d;

   // for two sites in the system, I will instead contract product
  // operators (e.g. sigX->.5 sigX x Id + .5 Id x sigX
  mwArray sig0v,sigXv,sigYv,sigZv,sigZsigZv;

  prepareOperatorProduct(sigZsigZv,sigZ,sigZ);
  prepareOperatorProduct(sig0v,sig0,sig0);
  prepareOperatorProduct(sigXv,sigX,sig0);
  {mwArray aux;prepareOperatorProduct(aux,sig0,sigX);
    sigXv=.5*(sigXv+aux);
  }
  prepareOperatorProduct(sigYv,sigY,sig0);
  {mwArray aux;prepareOperatorProduct(aux,sig0,sigY);
    sigYv=.5*(sigYv+aux);
  }
  prepareOperatorProduct(sigZv,sigZ,sig0);
  {mwArray aux;prepareOperatorProduct(aux,sig0,sigZ);
    sigZv=.5*(sigZv+aux);
  }
  
  ofstream* out;

  if(first){
    // Only the first time, if needed, write the header of the file
    if(cnt==0||!file_exists(outfname.data())){
      out=new ofstream(outfname.data());
      if(!out->is_open()){
	cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
	exit(1);
      }
      *out<<"%t\t real(tr(rho2))\timag(tr(rho2))\treal(<sigX>)\t imag(<sigX>)";
      *out<<"\t real(<sigY>)\t imag(<sigY>)";
      *out<<"\t real(<sigZ>)\t imag(<sigZ>)";
      *out<<"\t real(<sigZ sigZ>)\t imag(<sigZ sigZ>)"<<endl;
      out->close();delete out;
    }    
    first =false;
  }
  // Compute and write expectation values
  // Open the file to write them
  out=new ofstream(outfname.data(),ios::app);
      if(!out->is_open()){
	cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
	exit(1);
      }
  *out<<setprecision(15);

  mwArray res=sigXv*rdm2;
  complex_t resultX=res.getElement(Indices(0,0));
  res=sigYv*rdm2;
  complex_t resultY=res.getElement(Indices(0,0));
  res=sigZv*rdm2;
  complex_t resultZ=res.getElement(Indices(0,0));
  res=sig0v*rdm2;
  complex_t trV=res.getElement(Indices(0,0));
  res=sigZsigZv*rdm2;
  complex_t resultZZ=res.getElement(Indices(0,0));
  //complex_t resultZZ=ZERO_c;

  // and write to file and screen
  cout<<"t="<<cnt*delta<<", <O>="<<(resultX/trV)
      <<" E="<<(J_*resultZZ+g_*resultX+h_*resultZ)/trV<<endl;
  *out<<cnt*delta<<"\t"<<real(trV)<<"\t"<<imag(trV)
      <<"\t"<<real(resultX)<<"\t"<<imag(resultX)
      <<"\t"<<real(resultY)<<"\t"<<imag(resultY)
      <<"\t"<<real(resultZ)<<"\t"<<imag(resultZ)
      <<"\t"<<real(resultZZ)<<"\t"<<imag(resultZZ)
      <<endl;
  
  out->close();delete out;
}

// NOT GOOD!!! Would need to debug this carefully!
void computeAndSaveObservables(const MPS& mpsB,const Operator& middleUL,const Operator& rightmostOpL,
			       const mwArray& edgeUL, const mwArray& edgeUR,
			       int cnt,double delta,const string& outfname,
			       double J_,double g_,double h_){
  static bool first=true;
  // Operators needed (only local)
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1*I_c,ZERO_c};
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c}; //Z
  mwArray sig0=identityMatrix(d);
  mwArray sigX=mwArray(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  mwArray sigZ(Indices(d,d),dataZ);
  int d2=d*d;

  // for two sites in the system, I will instead contract product
  // operators (e.g. sigX->.5 sigX x Id + .5 Id x sigX
  mwArray sig0v,sigXv,sigYv,sigZv,sigZsigZv;

  prepareOperatorProduct(sigZsigZv,sigZ,sigZ);
  prepareOperatorProduct(sig0v,sig0,sig0);
  prepareOperatorProduct(sigXv,sigX,sig0);
  {mwArray aux;prepareOperatorProduct(aux,sig0,sigX);
    sigXv=.5*(sigXv+aux);
  }
  prepareOperatorProduct(sigYv,sigY,sig0);
  {mwArray aux;prepareOperatorProduct(aux,sig0,sigY);
    sigYv=.5*(sigYv+aux);
  }
  prepareOperatorProduct(sigZv,sigZ,sig0);
  {mwArray aux;prepareOperatorProduct(aux,sig0,sigZ);
    sigZv=.5*(sigZv+aux);
  }

  ofstream* out;
  
  if(first){
    // Only the first time, if needed, write the header of the file
    if(cnt==0||!file_exists(outfname.data())){
      out=new ofstream(outfname.data());
      if(!out->is_open()){
	cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
	exit(1);
      }
      *out<<"%t\t real(tr(rho2))\timag(tr(rho2))\treal(<sigX>)\t imag(<sigX>)";
      *out<<"\t real(<sigY>)\t imag(<sigY>)";
      *out<<"\t real(<sigZ>)\t imag(<sigZ>)";
      *out<<"\t real(<sigZ sigZ>)\t imag(<sigZ sigZ>)"<<endl;
      out->close();delete out;
    }    
    first =false;
  }
  // Compute and write expectation values
  // Open the file to write them
  out=new ofstream(outfname.data(),ios::app);
      if(!out->is_open()){
	cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
	exit(1);
      }
  *out<<setprecision(15);

  int Lt=mpsB.getLength(); // now I include the edge at the end of this
  int D0=mpsB.getA(Lt-1).getd();
  int xi=edgeUL.getDimension(1); // d2 x xi2
  Operator idOp(reshape(identityMatrix(xi),Indices(xi,1,xi,1)));
  Operator idOpD0(reshape(identityMatrix(D0),Indices(D0,1,D0,1)));
  Operator edgeLop(reshape(permute(edgeUL,Indices(2,1)),Indices(xi,1,d2,1)));
  Operator edgeRop(reshape(edgeUR,Indices(d2,1,xi,1)));

  MPO midMPO(Lt);
  MPO leftMPO(Lt),rightMPO(Lt); // to hold the edges
  for(int k=0;k<Lt-1;k++){
    midMPO.setOp(k,&middleUL,false);
    leftMPO.setOp(k,&idOp,false);
    rightMPO.setOp(k,&idOp,false);
  }
  midMPO.setOp(Lt-1,&rightmostOpL,false);
  leftMPO.setOp(Lt-1,&idOpD0,false);
  rightMPO.setOp(Lt-1,&idOpD0,false);
  leftMPO.setOp(0,&edgeLop,false);
  rightMPO.setOp(0,&edgeRop,false);

  
  MPO mpoCenter(Lt);
  const MPO* ptrs[]={&rightMPO,&midMPO,&midMPO,&leftMPO};
  MPO::join(4,ptrs,mpoCenter);
  
  // And (only because this case is symmetric) copy the left boundary, but complex, to use as bra
  MPS bra(mpsB);bra.conjugateMPS();
  //bra.gaugeCond('R',1);
  // try gauge again? 
  //MPS::gaugeCond(mpsB,bra,'L'); // maybe skip?

  Contractor& contractor=Contractor::theContractor();
  mwArray rdm2=contractor.contract(mpsB,mpoCenter,bra,-1,'R'); // 1 x d2*d2 x 1
  rdm2.reshape(Indices(d2*d2,1));

  mwArray res=sigXv*rdm2;
  complex_t resultX=res.getElement(Indices(0,0));
  res=sigYv*rdm2;
  complex_t resultY=res.getElement(Indices(0,0));
  res=sigZv*rdm2;
  complex_t resultZ=res.getElement(Indices(0,0));
  res=sig0v*rdm2;
  complex_t trV=res.getElement(Indices(0,0));
  res=sigZsigZv*rdm2;
  complex_t resultZZ=res.getElement(Indices(0,0));
  //complex_t resultZZ=ZERO_c;

  // and write to file and screen
  cout<<"t="<<cnt*delta<<", <O>="<<(resultX/trV)
      <<" E="<<(J_*resultZZ+g_*resultX+h_*resultZ)/trV<<endl;
  *out<<cnt*delta<<"\t"<<real(trV)<<"\t"<<imag(trV)
      <<"\t"<<real(resultX)<<"\t"<<imag(resultX)
      <<"\t"<<real(resultY)<<"\t"<<imag(resultY)
      <<"\t"<<real(resultZ)<<"\t"<<imag(resultZ)
      <<"\t"<<real(resultZZ)<<"\t"<<imag(resultZZ)
      <<endl;
  
  out->close();delete out;
  
}
