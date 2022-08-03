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

    Compute also XX correlators as we go

  */

int d=2;
double tolSVD=1E-12;      

void constructInitialProductTensors(int d,int D,int initSt,mwArray& W,mwArray& leftB);
void saveObservables(const mwArray& rdm2,int cnt,double delta,const string& outfname,
		     double J,double g,double h);

void saveCorrelations(const MPS& mpsB,const MPS& bra,const mwArray& W,const Operator& rightmostOpL,const Operator& middleUL,int cnt,double time,const string& outfnameCorr);

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
  int initSt=props.getIntProperty("initSt");
  double delta=props.getDoubleProperty("delta");
  int M=props.getIntProperty("M"); // max nr of time steps to be simulated
  //double alpha=props.getDoubleProperty("alpha");
  //if(alpha<0) alpha=1.;
  bool cutL=true;
  int Lmax=props.getIntProperty("Lmax"); // max nr of tensors to keep in the MPS
  if(Lmax<=0){Lmax=0;cutL=false;} // absolute min
  int rate=props.getIntProperty("rate");
  if(rate<=0) rate=1; // default frequency
  int Dcut=props.getIntProperty("D"); // the truncation  
  //  int Dcontr=props.getIntProperty("Drho"); // the dimension used for the effective rho+boundary as it evolves
  string outfname=props.getProperty("outputfile");
  string outfnameCorr=props.getProperty("outputfileCorr");
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

  
  // Prepare the initial tensors (initial state, boundary, edges, system)

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;

  // Basic terms for the exponential: I assume TI case, infinite
  // chain, but I only need the central (TI MPO) and the edges for a
  // finite case, which appear when the light-cone ends.
  MPO Uevol(3); // this is the double one
  //  mwArray midU_1;
  {
    IsingHamiltonian hamH(3,d,J_,g_,h_);
    //hamH.registerTimeDependence();
    MPO _Uevol_1(3); // held just temprarily
    hamH.getUMPO(_Uevol_1,-delta*I_c,0); // the single layer normal evolution
    //midU_1=_Uevol_1.getOp(1).getFullData();
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
  //mwArray edgeUR=Uevol.getOp(2).getFullData(); //d2 xi d2 1
  //cout<<"Recovered edgeUR from Uevol: "<<edgeUR<<endl;
  // For the edges, I only need the contraction with id2 above
  id2.reshape(Indices(1,d2));
  edgeUL.reshape(Indices(d2,d2*xi));edgeUL.multiplyLeft(id2);edgeUL.reshape(Indices(d2,xi)); //d2 x xi2
  //edgeUL.transpose(); // xi2 x d2
  //edgeUR.reshape(Indices(d2,xi*d2));edgeUR.multiplyLeft(id2);edgeUR.reshape(Indices(xi,d2)); // xi2 x d2

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
  //  DoubleOperator middleUL(permute(midU_1,Indices(4,1,2,3)),permute(conjugate(midU_1),Indices(4,1,2,3)));
  // if more than one step each time, I also need a reshaped boundary
  Operator extraUL(boundaryUL,Indices(4,1,3,2)); // xi x d2 x d2 x 1

  // if I renormalize the edge, I will need a link operator from left
  // to right, to convert the basis (iso*isoT). When not renormalized,
  // this is just the identity of the oriinal bond dimension
  mwArray linkOp=identityMatrix(D0*D0);
  
  int cnt=0;double time=0.;

  // For the first stp, the rdm is special
  mwArray rdm0(W);// d2xD02xD02
  rdm0.permute(Indices(2,1,3));rdm0.reshape(Indices(D0*D0,d2*D0*D0));
  rdm0.multiplyLeft(leftB);rdm0.reshape(Indices(d2,D0*D0));
  rdm0.multiplyRight(permute(rdm0,Indices(2,1))); // d2xd2
  rdm0.permute(Indices(2,1));rdm0.reshape(Indices(d2*d2,1)); // as will appear later, with index of rightmost index first
    
  saveObservables(rdm0,cnt,delta,outfname,J_,g_,h_);

  
  clock_t start,finish; start=clock();
  while(cnt<M){
  // step by step:
  // 1) evolve the boundary (and grow it) 
  // 2) increase the MPO(s) for the system
  // 3) contract expectation values

    // how many time steps to add now?
    int toadd=1; // to be redone for quadratic LC
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
    if(toadd>1){ // TODO!!
      cout<<"Not yet implemented!!"<<endl;exit(1);
    }
    else{
      evol.setOp(0,&leftmostOpL,false);
      evol.setOp(1,&secondOpL,false); // which should be at pos lenB from R
      evol.setOp(Lt-1,&rightmostOpL,false);
      for(int k=2;k<Lt-1;k++)
	evol.setOp(k,&middleUL,false);
    }
    //cout<<"Prepared evolution MPO: "<<evol<<endl;
    mpsB.approximate(evol,aux,Dcut); // first estimate by SVDs
    contractor.optimize(evol,aux,mpsB,Dcut);
    //    cout<<"After optimize,mpsB:"<<mpsB<<endl;
    //    MPS mpsIso(toadd+2,Dcut,xi);
    if(cutL){
      if(Lt>Lmax){
	// Start cutting only when computational advantage
	int toblock=Lt-Lmax+1;
	if(mpsB.getA(Lt-1-toblock).getDr()>=Dcut){
	  //cout<<"Cutting the MPS (len Lt "<<Lt<<", Lmax="<<Lmax<<"), toblock="<<toblock<<endl;
	  MPS mpsIso(toblock+1,Dcut,xi);
	  // should decide the isometry. Here from the MPS
	  mpsB.gaugeCond('L',1); // I could just apply the gauge to the last toadd+1 sites on the right	
	  //cout<<"Checking gauge to pos "<<Lt-1-toadd<<": "<<contractor.contract(mpsB,mpsB,Lt-1-toadd,'R')<<endl;
	  // Now my isometry is just a piece of the MPS, which I keep as a MPS with toadd+2 sites, where the leftmost tensor is just an identity xi
	  // for(int k=0;k<toadd+1;k++){
	    //	    mpsIso.replaceSite(toadd+1-k,mpsB.getA(Lt-1-k).getA(),false);
	  for(int k=0;k<toblock;k++){
	    mpsIso.replaceSite(toblock-k,mpsB.getA(Lt-1-k).getA(),false);
	  }
	  //int Diso=mpsB.getA(Lt-1-toadd).getDl();
	  int Diso=mpsB.getA(Lt-toblock).getDl();
	  mpsIso.replaceSite(0,reshape(identityMatrix(Diso),Indices(Diso,1,Diso)),false);
	  //cout<<"MPS for isometry:"<<mpsIso<<endl;
	  //cout<<"Checking contraction: "<<contractor.contract(mpsIso,mpsIso,0,'R')<<endl;
	  // midU // d2 x xi x d2 x xi
	  // use it to transform the rightmostOp
	  //MPO mpoCol(toadd+2);
	  //mpoCol.setOp(toadd+1,&rightmostOpL,false);
	  //for(int p=0;p<toadd;p++) mpoCol.setOp(toadd-p,&middleUL,false);
	  MPO mpoCol(toblock+1);
	  mpoCol.setOp(toblock,&rightmostOpL,false); 
	  for(int p=0;p<toblock-1;p++) mpoCol.setOp(toblock-1-p,&middleUL,false);
	  mwArray aux=contractor.contract(mpsIso,mpoCol,mpsIso,0,'R'); // Diso (r) x d2 x Diso(l)
	  rightmostOpL.setRotatedData(reshape(aux,Indices(Diso,d2,Diso,1)),Indices(1,2,3,4)); // trick to replace the data in Op
	  // And I also need to replace my mpsB by the cut version!
	  // MPS auxMPS(Lt-toadd,Dcut,xi);
	  // for(int p=0;p<Lt-toadd-1;p++) auxMPS.replaceSite(p,mpsB.getA(p).getA(),false); // there must be omething more efficient
	  // auxMPS.replaceSite(Lt-toadd-1,reshape(identityMatrix(Diso),Indices(Diso,Diso,1)),false);
	  MPS auxMPS(Lt-toblock+1,Dcut,xi);
	  for(int p=0;p<Lt-toblock;p++) auxMPS.replaceSite(p,mpsB.getA(p).getA(),false); // there must be omething more efficient
	  auxMPS.replaceSite(Lt-toblock,reshape(identityMatrix(Diso),Indices(Diso,Diso,1)),false);
	  mpsB=auxMPS;
	  // And prepare the link operator that converts between left and right part
	  aux=linkOp;
	  MPS mpsIsoT(mpsIso);mpsIsoT.conjugateMPS();
	  //	  mpsIsoT.applyLocalOperator(toadd+1,aux,false);
	  mpsIsoT.applyLocalOperator(toblock,aux,false);
	  linkOp=contractor.contract(mpsIso,mpsIsoT,0,'R'); // Dcut (r) x Dcut
	  linkOp.Hconjugate(); // As I need to apply it on bra always
	  //exit(1);
	}
	// else{
	//   cout<<"Not cutting because although Lt="<<Lt<<">Lmax="<<Lmax<<" the bond dim of link "<<Lt-1-toblock<<" (toblock="<<toblock<<") is "<<mpsB.getA(Lt-1-toblock).getDr()<<"<"<<Dcut<<endl;
	// }
      }
    }

    // //mpsB.gaugeCond('R',1);
    mpsB.gaugeCond('R',1);
    cnt+=toadd;
    time+=toadd*delta;

    // to compute expectation values on two sites, I only need this
    // mps. If I want a larger subsystem, I should used the commented
    // part or a modification of it (as this one produces the rdm one
    // step forward!!)
    MPS bra(mpsB);bra.conjugateMPS();
    bra.applyLocalOperator(bra.getLength()-1,linkOp,false);
    //bra.gaugeCond('R',1);
    // try gauge again? 
    MPS::gaugeCond(mpsB,bra,'L');
    mwArray rdm2=contractor.contract(mpsB,bra,0,'R'); // d2 x 1 x d2
    rdm2.reshape(Indices(d2,d2));
    rdm2.multiplyLeft(reshape(conjugate(bra.getA(0).getA()),Indices(d2,d2)));
    rdm2.multiplyRight(permute(reshape(mpsB.getA(0).getA(),Indices(d2,d2)),Indices(2,1)));
    //cout<<"Here rdm2="<<rdm2<<endl;if(cnt>10)exit(1);
    rdm2.reshape(Indices(d2*d2,1));
    saveObservables(rdm2,cnt,delta,outfname,J_,g_,h_);
        
    saveCorrelations(mpsB,bra,W,rightmostOpL,middleUL,cnt,time,outfnameCorr);
    finish=clock();
    cout<<"\t time for step "<<cnt<<": "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
    start=clock();

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
  case 4: // Id!
    {
      W=.5*sig0;
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


void saveCorrelations(const MPS& mpsB,const MPS& bra,const mwArray& W,const Operator& rightmostOpL,const Operator& middleUL,int cnt,double time,const string& outfname){
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
  
  ofstream* out;
  Contractor& contractor=Contractor::theContractor();
  // static Operator* opXW;
  // static Operator* opYW;
  // static Operator* opZW;

  if(first){
    // Only the first time, if needed, write the header of the file
    if(cnt==0||!file_exists(outfname.data())){
      out=new ofstream(outfname.data());
      if(!out->is_open()){
	cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
	exit(1);
      }
      *out<<"%t\t real(tr(rho2))\timag(tr(rho2))\treal(<sigX(t) sigX>)\t imag(<sigX(t) sigX>)";
      *out<<"\t real(<sigY(t) sigZ>)\t imag(<sigY(t) sigY>)";
      *out<<"\t real(<sigZ(t) sigZ>)\t imag(<sigZ(t) sigZ>)";
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

  int Lt=mpsB.getLength();
  int D0=W.getDimension(1);
  // To compute correlations at different times, I need a special  MPO
  MPO midMPO(Lt);midMPO.setOp(Lt-1,&rightmostOpL,false);
  midMPO.setOp(0,new Operator(reshape(permute(reshape(identityMatrix(d*d*d),Indices(d,d,d,d,d,d)),Indices(1,4,2,5,3,6)),Indices(d*d,1,d*d,d*d))),true);
  for(int k=1;k<Lt-1;k++)midMPO.setOp(k,&middleUL,false);
  complex_t norm=contractor.contract(mpsB,midMPO,bra);
  // XX
  midMPO.setOp(0,new Operator(reshape(reshape(reshape(identityMatrix(d),Indices(d*d,1))
					      *reshape(identityMatrix(d),Indices(1,d*d)),Indices(d*d*d*d,1))
				      *reshape(sigX,Indices(1,d*d)),Indices(d*d,1,d*d,d*d))),true);
  midMPO.setOp(Lt-1,new Operator(reshape(sigX*reshape(W,Indices(d,d*D0*D0*D0*D0)),Indices(d2,D0*D0,D0*D0,1)),Indices(3,1,2,4)),true);
  complex_t valXX=contractor.contract(mpsB,midMPO,bra);
  // YY
  midMPO.setOp(0,new Operator(reshape(reshape(reshape(identityMatrix(d),Indices(d*d,1))
					      *reshape(identityMatrix(d),Indices(1,d*d)),Indices(d*d*d*d,1))
				      *reshape(sigY,Indices(1,d*d)),Indices(d*d,1,d*d,d*d))),true);
  midMPO.setOp(Lt-1,new Operator(reshape(conjugate(sigY)*reshape(W,Indices(d,d*D0*D0*D0*D0)),Indices(d2,D0*D0,D0*D0,1)),Indices(3,1,2,4)),true);
  complex_t valYY=contractor.contract(mpsB,midMPO,bra);
  // ZZ
  midMPO.setOp(0,new Operator(reshape(reshape(reshape(identityMatrix(d),Indices(d*d,1))
					      *reshape(identityMatrix(d),Indices(1,d*d)),Indices(d*d*d*d,1))
				      *reshape(sigZ,Indices(1,d*d)),Indices(d*d,1,d*d,d*d))),true);
  midMPO.setOp(Lt-1,new Operator(reshape(sigZ*reshape(W,Indices(d,d*D0*D0*D0*D0)),Indices(d2,D0*D0,D0*D0,1)),Indices(3,1,2,4)),true);
  complex_t valZZ=contractor.contract(mpsB,midMPO,bra);
  *out<<time<<"\t";
  *out<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(valXX)<<"\t"<<imag(valXX)<<"\t";
  *out<<real(valYY)<<"\t"<<imag(valYY)<<"\t";
  *out<<real(valZZ)<<"\t"<<imag(valZZ)<<"\t";
  *out<<endl;
  out->close();delete out;
  cout<<"t="<<time<<", <X(t)X>="<<valXX/norm<<", <Y(t)Y>="<<valYY/norm<<", <Z(t)Z>="<<valZZ/norm<<endl;
  
}
