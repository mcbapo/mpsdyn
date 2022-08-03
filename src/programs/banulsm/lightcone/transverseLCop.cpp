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
    increases using the light cone structure.  This version
    concentrates in correlations at different times, averagd in a
    (potentially thermal, but now infinite temperature) state, namely
    <O_x O_0(t)>=tr(O_x U O_0 U^dagger).  (25.11.2021)

    After checking, the most stable option to comppute correlators (at
    same position, different times) is to use a MPO with proper
    insertion of the operato (rather than minimal contraciton of the
    MPSs,for which error seems larger. (28.11.2021)

    Now also the possibility to compute correlations at distance
    larger than 1 (29.11.2021)
    And correlations between a pair of MPOs (of length exactly two, for the moment)

    TODO: A better implementation could encapsulate all the Operators,
    which need to be computed only once, in a class, and then offer
    methods to compute the correlations and evolve the boundaries.
  */

int d=2;
double tolSVD=1E-12;      
int labelE=17; // the label with which I indicate energy correlators


void constructInitialProductTensors(int d,int D,int initSt,mwArray& W,mwArray& leftB);
void computeCorrelators(MPS& mpsB,int cnt,double delta,const string& outfname,
			const Operator& middleUL,
			const mwArray& specialOp,const mwArray& specialOpb);
void computeCorrelatorsSameStep(const MPS& mpsB,const MPS& bra,
				     const Operator& middleUL,const mwArray& specialOp,const mwArray& specialOpb,
				     int cnt,double delta,const string& outfname);
void computeCorrelatorsStepAfter(const MPS& mpsB,const MPS& bra,const Operator& middleUL,
				  const mwArray& specialOp,const mwArray& specialOpb,
				  int cnt,double delta,const string& outfname);
void computeCorrelatorsTwoStepsAfter(const MPS& mpsB,const MPS& bra,
				     const Operator& middleUL,const mwArray& specialOp,const mwArray& specialOpb,
				     int cnt,double delta,const string& outfname);
void computeDistantCorrelators(const MPS& mpsB,int maxM,
			       const Operator& middleUL,const mwArray& specialOp,const mwArray& specialOpb,
			       int cnt,double delta,const string& outfname);

void computeDistantCorrelators(const MPS& mpsB,
			       const mwArray& midU,
			       const mwArray& edgeUL,const mwArray& edgeULb,
			       const mwArray& edgeUR,const mwArray& edgeURb,
			       int cnt,double delta,const string& outfname,
			       int Dcut,int maxM,int whichOp);

// Compute two body correlators at distance >=1, for products of two Pauli matrices
// NOT PROPERLY TESTED!!
void computeDistantCorrelatorsTwoBody(const MPS& mpsB,
				      const mwArray& midU,
				      const mwArray& edgeUL,const mwArray& edgeULb,
				      const mwArray& edgeUR,const mwArray& edgeURb,
				      int cnt,double delta,const string& outfname,
				      int Dcut,int maxM,int whichOp1u,int whichOp2u,int whichOp1b,int whichOp2b);

// As above, but for two MPOs, and the last integer is the label that will appear in the file
void computeDistantCorrelatorsTwoBody(const MPS& mpsB,
				      const mwArray& midU,
				      const mwArray& edgeUL,const mwArray& edgeULb,
				      const mwArray& edgeUR,const mwArray& edgeURb,
				      int cnt,double delta,const string& outfname,
				      int Dcut,int maxM,const MPO& mpoLow,const MPO& mpoUp,int label);
// As above, but with the optimization that the right vector is not optimized independently but computed from the left one (under test)
void computeDistantCorrelatorsTwoBodyOpt(const MPS& mpsB,
				      const mwArray& midU,
				      const mwArray& edgeUL,const mwArray& edgeULb,
				      const mwArray& edgeUR,const mwArray& edgeURb,
				      int cnt,double delta,const string& outfname,
				      int Dcut,int maxM,const MPO& mpoLow,const MPO& mpoUp,int label);
void computeMPOCorrelatorsSameSite(MPS& mpsB,const mwArray& midU,
				   const mwArray& edgeUL,const mwArray& edgeULb,
				   const mwArray& edgeUR,const mwArray& edgeURb,
				   int cnt,double delta,const string& outfname,
				   int Dcut,int maxM,const MPO& mpoLow,const MPO& mpoUp,int label);

// The special terms where t=l, which have 1D structure in the TN
void computeDiagonalCorrelators(const mwArray& edgeUL,const mwArray& edgeURb,
				double delta,const string& outfname,int maxM);
void computeMPODiagonalCorrelators(const mwArray& edgeUL,const mwArray& edgeURb,const mwArray& midU,
				   double delta,const string& outfname,int maxM,const MPO& mpoLow,
				   const MPO& mpoUp,int label,int Dcut);

// auxiliary, to dot he manipulations required to contract with the links shifted and edges traced out
complex_t contractWithTracedEdges(MPS& left,MPS& right,const MPO& mpoTraces);
// auxiliary functions to aply one "empty" column to the left or right (but flipped) vectors
void applyOneStepLeft(MPS& mps,const MPO& oplow,int Dcut);
void applyOneStepRight(MPS& mps,const MPO& ophigh,int Dcut);

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
  //int rate=props.getIntProperty("rate");
  //if(rate<=0) rate=1; // default frequency
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

  bool allCs=props.getIntProperty("computeAll")==1; // all correlators between Pauli matrices (diff distances) computed
  bool energyCs=props.getIntProperty("computeEnergyCorr")==1; // compute correlators of the energy density
  //  string outfnameCorr=props.getProperty("outputCorr"); // for two-body ones
  // if(allCs&outfnameCorr.empty()){
  //   cout<<"ERROR! No file specified to save correlators at different distances"<<endl;
  //   exit(1);
  // }
  
  
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
  //  mwArray midU_1; // for a single layer (I need it to get the inverse and apply mid steps)
  {
    IsingHamiltonian hamH(3,d,J_,g_,h_);
    //hamH.registerTimeDependence();
    MPO _Uevol_1(3); // held just temprarily
    hamH.getUMPO(_Uevol_1,-delta*I_c,0); // the single layer normal evolution
    //midU_1=_Uevol_1.getOp(1).getFullData();
    //putForMatlab(cout,_Uevol_1.getOp(1).getFullData(),"Umid");
    // But I need to double it
    doubleMPO(_Uevol_1,Uevol,true);
    cout<<"Created the double Uevol"<<endl;
  }
  
  // The most repeated tensor in the network (doubled U)
  // Notice it is not using the double structure! May be optimized
  mwArray midU(Uevol.getOp(1).getFullData()); // d2 x xi x d2 x xi
  //  putForMatlab(cout,midU,"midU");
  int xi=midU.getDimension(1);
  //mwArray midUl(midU);midUl.permute(Indices(1,4,2,3));midUl.reshape(Indices(d2*xi,xi*d2));
  //  midU.reshape(Indices(d2*xi,d2*xi)); // to contract with right boundary

  // /*/ // // Boundary (without tracing)
  // /*/ // mwArray boundaryUL=Uevol.getOp(0).getFullData(); //d2 1 d2 xi
  // /*/ // mwArray boundaryUR=Uevol.getOp(2).getFullData(); //d2 xi d2 1

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

  // and the edge "downwards"
  mwArray edgeULb=Uevol.getOp(0).getFullData(); //d2 1 d2 xi
  edgeULb.permute(Indices(3,1,2,4)); // d2(b) d2(u) 1 xi 
  edgeULb.reshape(Indices(d2,d2*xi));edgeULb.multiplyLeft(id2);edgeULb.reshape(Indices(d2,xi));

  id2.reshape(Indices(d2,1));
  mwArray edgeURb=Uevol.getOp(2).getFullData(); //d2 xi d2 1
  edgeURb.reshape(Indices(d2*xi,d2));edgeURb.multiplyRight(id2);
  edgeURb.reshape(Indices(d2,xi));edgeURb.permute(Indices(2,1)); // xi (l) x d2 (u)
  
  int D0=d2; // initial bond dim

  MPS mpsB(2,d2,D0);
  {
    mpsB.replaceSite(0,reshape(identityMatrix(d2),Indices(d2,1,d2)),false);
    mpsB.replaceSite(1,reshape(identityMatrix(d2),Indices(d2,d2,1)),false);
    cout<<"Initialized boundary MPS: "<<mpsB<<endl;
    // In this case, it is trivial, as it is only the identity
  }
  int lenB=mpsB.getLength(); // 1+ nr of steps 

  // Common pieces for the evolution operator are the edge tensor(s) 

  Operator leftmostOpL(reshape(identityMatrix(d2),Indices(d2,1,1,d2)));
  Operator secondOpL(reshape(edgeUL*reshape(permute(midU,Indices(2,1,3,4)),
					    Indices(xi,d2*d2*xi)), // d2 x d2*d2*xi
			     Indices(d2,d2,d2,xi)),Indices(4,2,1,3)); // xi x d2 x d2 (incl edge) x d2
  Operator rightmostOpL(reshape(identityMatrix(d2),Indices(d2,d2,1,1)));
  // next one identical to secondOpL in this case (I could skip the copy!)
  Operator rightsecondOpL(reshape(edgeULb*reshape(permute(midU,Indices(2,1,3,4)),
						  Indices(xi,d2*d2*xi)), // d2 x d2*d2*xi
				  Indices(d2,d2,d2,xi)),Indices(4,2,1,3)); // xi x d2 x d2 (incl edge) x d2
  Operator middleUL(midU,Indices(4,1,2,3)); // xi(r) x d2(u) x xi x d2 // This is actually faster than the double
  //DoubleOperator middleUL(permute(midU_1,Indices(4,1,2,3)),permute(conjugate(midU_1),Indices(4,1,2,3))); // less efficient!
  // Special op with two edges and midU, to be able to compute intermediate steps
  mwArray specialOp=permute(reshape(reshape(edgeUL*reshape(permute(midU,Indices(2,1,3,4)),
							   Indices(xi,d2*d2*xi)),Indices(d2*d2*d2,xi)) // d2 x d2*d2*xi
				    *edgeUR,
				    Indices(d2,d2,d2,d2)),Indices(4,2,1,3)); // d2 (incl edgeR) x d2 x d2 (incl edge) x d2
  mwArray specialOpb=permute(reshape(reshape(edgeULb*reshape(permute(midU,Indices(2,1,3,4)),
							     Indices(xi,d2*d2*xi)),Indices(d2*d2*d2,xi)) // d2 x d2*d2*xi
				     *edgeURb,
				     Indices(d2,d2,d2,d2)),Indices(4,2,1,3)); // d2 (incl edgeR) x d2 x d2 (incl edge) x d2

  // I also need 
  
  // /*/ // if more than one step each time, I also need a reshaped boundary
  // /*/ Operator extraUL(boundaryUL,Indices(4,1,3,2)); // xi x d2 x d2 x 1

  // if I renormalize the edge, I will need a link operator from left
  // to right, to convert the basis (iso*isoT). When not renormalized,
  // this is just the identity of the oriinal bond dimension
  mwArray linkOp=identityMatrix(D0*D0);
  
  int cnt=0;double time=0.;

  // To compute the correlator, I prepare a special local operator
  // that I apply to the first and last sites of the MPS before
  // contracting with the transpose
  computeCorrelators(mpsB,cnt,delta,outfname,middleUL,specialOp,specialOpb);
  if(allCs){
    computeDiagonalCorrelators(edgeUL,edgeURb,delta,outfname,M);
    computeDistantCorrelators(mpsB,midU,edgeUL,edgeULb,edgeUR,edgeURb,
			      cnt,delta,outfname,
			      Dcut,M,1);
    computeDistantCorrelators(mpsB,midU,edgeUL,edgeULb,edgeUR,edgeURb,
			      cnt,delta,outfname,
			      Dcut,M,2);
    computeDistantCorrelators(mpsB,midU,edgeUL,edgeULb,edgeUR,edgeURb,
			      cnt,delta,outfname,
			      Dcut,M,3);
  }
  IsingHamiltonian hamH(2,d,J_,.5*g_,.5*h_);
  const MPO& eDens=hamH.getHMPO(); // special MPO to be able to compute the energy density correlator, constructed as the H for 2 sites with single site terms divided by 2
  if(energyCs){
    computeMPODiagonalCorrelators(edgeUL,edgeURb,midU,delta,outfname,M,eDens,eDens,labelE,Dcut);
    computeMPOCorrelatorsSameSite(mpsB,midU,edgeUL,edgeULb,edgeUR,edgeURb,
				  cnt,delta,outfname,
				  Dcut,M,eDens,eDens,labelE);
    computeDistantCorrelatorsTwoBodyOpt(mpsB,midU,edgeUL,edgeULb,edgeUR,edgeURb,
				     cnt,delta,outfname,
				     Dcut,M,eDens,eDens,labelE);
  }

  //  saveObservables(rdm0,cnt,delta,outfname,J_,g_,h_);

  clock_t start,finish; start=clock();
  while(cnt<M){
  // step by step:
  // 1) evolve the boundary (and grow it) 
  // 2) increase the MPO(s) for the system
  // 3) contract expectation values

    // how many time steps to add now?
    int toadd=1; // to be redone for quadratic LC
    lenB=mpsB.getLength(); // 2+ current nr of steps 
    // initialize the MPS, expanding on both edges
    int Lt=lenB+2*toadd;
    MPS aux(Lt,Dcut,xi);
    for(int k=0;k<lenB;k++){
      aux.replaceSite(toadd+k,mpsB.getA(k),false);
    }
    for(int k=0;k<toadd;k++){
      aux.replaceSite(k,reshape(mwArray(ONE_c),Indices(1,1,1)),false);
      aux.replaceSite(Lt-1-k,reshape(mwArray(ONE_c),Indices(1,1,1)),false);
    }
    //cout<<"Step"<<cnt<<" Extended the MPS "<<mpsB<<" to "<<aux<<endl;
    MPO evol(Lt);
    if(toadd>1){ // TODO!!
      cout<<"Not yet implemented!!"<<endl;exit(1);
    }
    else{
      evol.setOp(0,&leftmostOpL,false);
      evol.setOp(1,&secondOpL,false); // which should be at pos lenB from R
      evol.setOp(Lt-1,&rightmostOpL,false);
      //evol.setOp(Lt-2,&secondOpL,false);
      evol.setOp(Lt-2,&rightsecondOpL,false);
      for(int k=2;k<Lt-2;k++)
	evol.setOp(k,&middleUL,false);
    }
    //cout<<"Prepared evolution MPO: "<<evol<<endl;
    mpsB.approximate(evol,aux,Dcut); // first estimate by SVDs (might be better starting somewhere else?)
    contractor.optimize(evol,aux,mpsB,Dcut);
    //    cout<<"After optimize,mpsB:"<<mpsB<<endl;

    // //mpsB.gaugeCond('R',1);
    mpsB.gaugeCond('R',1);
    cnt+=2*toadd;
    time+=2*toadd*delta;

    finish=clock();
    cout<<"\t time for step "<<cnt<<": "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;

    computeCorrelators(mpsB,cnt,delta,outfname,middleUL,specialOp,specialOpb);
    if(allCs){
      start=clock();
      computeDistantCorrelators(mpsB,midU,edgeUL,edgeULb,edgeUR,edgeURb,
				cnt,delta,outfname,
				Dcut,M,1);
      computeDistantCorrelators(mpsB,midU,edgeUL,edgeULb,edgeUR,edgeURb,
				cnt,delta,outfname,
				Dcut,M,2);
      computeDistantCorrelators(mpsB,midU,edgeUL,edgeULb,edgeUR,edgeURb,
				cnt,delta,outfname,
				Dcut,M,3);
      finish=clock();
      cout<<"\t time for correlators from step "<<cnt<<": "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
    }
    if(energyCs){
      start=clock();
      computeMPOCorrelatorsSameSite(mpsB,midU,edgeUL,edgeULb,edgeUR,edgeURb,
				    cnt,delta,outfname,
				    Dcut,M,eDens,eDens,labelE);
      computeDistantCorrelatorsTwoBodyOpt(mpsB,midU,edgeUL,edgeULb,edgeUR,edgeURb,
				       cnt,delta,outfname,
				       Dcut,M,eDens,eDens,labelE);
      finish=clock();
      cout<<"\t time for energy correlators from step "<<cnt<<": "<<(finish-start)*1./CLOCKS_PER_SEC<<endl;
    }
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

void computeCorrelators(MPS& mpsB,int cnt,double delta,const string& outfname,
			const Operator& middleUL,const mwArray& specialOp,const mwArray& specialOpb){
  static bool first=true;
  // Operators needed (only local)
  // complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  // complex_t dataY[]={ZERO_c,I_c,-1*I_c,ZERO_c};
  // complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c}; //Z
  // mwArray sig0=identityMatrix(d);
  // mwArray sigX=mwArray(Indices(d,d),dataX);
  // mwArray sigY(Indices(d,d),dataY);
  // mwArray sigZ(Indices(d,d),dataZ);
  // int d2=d*d;

  // To compute the correlator, I prepare a special local operator
  // that I apply to the first and last sites of the MPS before
  // contracting with the transpose
  // mwArray op0=reshape(sig0,Indices(d*d,1))*reshape(sig0,Indices(1,d*d));
  // mwArray opX=reshape(sig0,Indices(d*d,1))*reshape(sigX,Indices(1,d*d));
  // mwArray opY=reshape(sig0,Indices(d*d,1))*reshape(sigY,Indices(1,d*d));
  // mwArray opZ=reshape(sig0,Indices(d*d,1))*reshape(sigZ,Indices(1,d*d));
  // the one in 0 needs to be transposed, which for all is the same except Y, which needs to be conj
  
  ofstream* out;

  if(first){
    // Only the first time, if needed, write the header of the file
    if(cnt==0||!file_exists(outfname.data())){
      out=new ofstream(outfname.data());
      if(!out->is_open()){
	cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
	exit(1);
      }
      *out<<"%t\t l\t alpha(1,2,3)\t real(tr(Id))\timag(tr(Id))";
      *out<<"\t real(tr(sig_alpha(t,l)sig_alpha(0,0))) \t imag(tr(sig_alpha(t,l)sig_alpha(0,0))";
      //*out<<"\t real(tr(sigy(t)sigy)) \t imag(tr(sigy(t)sigy))";
      //*out<<"\t real(tr(sigz(t)sigz)) \t imag(tr(sigz(t)sigz))";
      *out<<endl;
      out->close();delete out;
    }    
    first =false;
  }
  // Compute and write expectation values
  // Open the file to write them
  // out=new ofstream(outfname.data(),ios::app);
  // if(!out->is_open()){
  //   cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
  //   exit(1);
  // }
  // *out<<setprecision(15);
  MPS bra(mpsB);bra.conjugateMPS();
  MPS::gaugeCond(mpsB,bra,'L');
  computeCorrelatorsStepAfter(mpsB,bra,middleUL,specialOp,specialOpb,cnt,delta,outfname);
  computeCorrelatorsTwoStepsAfter(mpsB,bra,middleUL,specialOp,specialOpb,cnt,delta,outfname);

  // Contractor& contractor=Contractor::theContractor();

  // MPO link(mpsB.getLength());
  // Operator* idXi;
  // if(mpsB.getLength()>2){
  //   int xi=mpsB.getA(1).getd();
  //   idXi=new Operator(reshape(identityMatrix(xi),Indices(xi,1,xi,1)));
  // }
  // for(int k=1;k<mpsB.getLength()-1;k++)
  //   link.setOp(k,idXi,false); 
  // link.setOp(0,new Operator(op0),true);
  // link.setOp(mpsB.getLength()-1,new Operator(op0),true);

  // complex_t norm=contractor.contract(mpsB,link,bra);

  // link.setOp(0,new Operator(opX),true);
  // link.setOp(mpsB.getLength()-1,new Operator(opX),true);
  // complex_t valX=contractor.contract(mpsB,link,bra);

  // link.setOp(0,new Operator(conjugate(opY)),true);
  // link.setOp(mpsB.getLength()-1,new Operator(opY),true);
  // complex_t valY=contractor.contract(mpsB,link,bra);

  // link.setOp(0,new Operator(opZ),true);
  // link.setOp(mpsB.getLength()-1,new Operator(opZ),true);
  // complex_t valZ=contractor.contract(mpsB,link,bra);

  // *out<<cnt*delta<<"\t";
  // *out<<real(norm)<<"\t"<<imag(norm)<<"\t";
  // *out<<real(valX)<<"\t"<<imag(valX)<<"\t";
  // *out<<real(valY)<<"\t"<<imag(valY)<<"\t";
  // *out<<real(valZ)<<"\t"<<imag(valZ)<<"\t";
  // *out<<endl;

  // // I can also compute one step before, 
  
  // out->close();delete out;

  // cout<<"t="<<cnt*delta<<", <X(t)X>="<<valX/norm<<", <Y(t)Y>="<<valY/norm<<", <Z(t)Z>="<<valZ/norm<<endl;
  // if(mpsB.getLength()>2){ delete idXi;}
  //  computeCorrelatorsSameStep(mpsB,bra,middleUL,specialOp,specialOpb,cnt,delta,outfname);
    

}



void computeCorrelatorsStepAfter(const MPS& mpsB,const MPS& bra,
				  const Operator& middleUL,const mwArray& specialOp,const mwArray& specialOpb,
				  int cnt,double delta,const string& outfname){
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

  // To compute the correlator, I construct a MPO

  static Operator *lastOpX,*lastOpY,*lastOpZ,*lastOp0;
  static Operator *firstOpX,*firstOpY,*firstOpZ,*firstOp0;
  if(first){
    // The ones at the beginning and end of the chain are different, in this case
    mwArray aux(specialOp); // d2(r) d2(u) d2 (l) d2 (d)
    cout<<"specialOp:"<<specialOp.getDimensions()<<endl;
    aux.reshape(Indices(d2,d2,d2*d2));
    aux.permute(Indices(2,1,3)); // d2(u) x d2(r) d2(l) d2(d)
    aux.reshape(Indices(d2,d2*d2*d2));
    aux.multiplyLeft(reshape(identityMatrix(d),Indices(1,d2))); // 1 x d2(r) d2(l) d2(r)
    aux.reshape(Indices(d2*d2,d,d));
    aux.permute(Indices(1,3,2)); // d2(r)*d2(l)*d'(d) x d(d)
    aux.reshape(Indices(d2*d2*d,d)); 
    firstOpX=new Operator(reshape(permute(reshape(aux*sigX,Indices(d2*d2,d,d)),Indices(1,3,2)),Indices(d2,1,d2,d2)));
    firstOpY=new Operator(reshape(permute(reshape(aux*sigY,Indices(d2*d2,d,d)),Indices(1,3,2)),Indices(d2,1,d2,d2)));
    firstOpZ=new Operator(reshape(permute(reshape(aux*sigZ,Indices(d2*d2,d,d)),Indices(1,3,2)),Indices(d2,1,d2,d2)));
    firstOp0=new Operator(reshape(permute(reshape(aux*sig0,Indices(d2*d2,d,d)),Indices(1,3,2)),Indices(d2,1,d2,d2)));

    aux=specialOpb; // d2 (incl edgeR) x d2(u) x d2 (incl edgeL) x d2
    aux.reshape(Indices(d2*d2*d2,d2));
    lastOpX=new Operator(reshape(aux*reshape(sigX,Indices(d2,1)),Indices(d2,d2,d2,1)));
    lastOpY=new Operator(reshape(aux*reshape(sigY,Indices(d2,1)),Indices(d2,d2,d2,1)));
    lastOpZ=new Operator(reshape(aux*reshape(sigZ,Indices(d2,1)),Indices(d2,d2,d2,1)));
    lastOp0=new Operator(reshape(aux*reshape(sig0,Indices(d2,1)),Indices(d2,d2,d2,1)));
    
    first=false;
  }
  
  // Open the file to write them
  ofstream* out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
    exit(1);
  }
  *out<<setprecision(15);

  Contractor& contractor=Contractor::theContractor();

  //  cout<<"Step after "<<cnt<<" (length of mps:"<<mpsB.getLength()<<")"<<endl;
  MPO link(mpsB.getLength());
  for(int k=1;k<mpsB.getLength()-1;k++)
    link.setOp(k,&middleUL,false); 
  link.setOp(0,firstOp0,false);
  link.setOp(mpsB.getLength()-1,lastOp0,false);

  complex_t norm=contractor.contract(mpsB,link,bra);
  link.setOp(0,firstOpX,false);
  link.setOp(mpsB.getLength()-1,lastOpX,false);
  complex_t valX=contractor.contract(mpsB,link,bra);

  link.setOp(0,firstOpY,false);
  link.setOp(mpsB.getLength()-1,lastOpY,false);
  complex_t valY=contractor.contract(mpsB,link,bra);

  link.setOp(0,firstOpZ,false);
  link.setOp(mpsB.getLength()-1,lastOpZ,false);
  complex_t valZ=contractor.contract(mpsB,link,bra);

  // I can also compute ZZ for the energy density:
  const MPO* ptrs[]={&link,&link};
  MPO link2(mpsB.getLength());
  MPO::join(2,ptrs,link2);
  complex_t valZZ=contractor.contract(mpsB,link2,bra);
  
  *out<<(cnt+1)*delta<<"\t"<<0<<"\t"<<1<<"\t"<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(valX)<<"\t"<<imag(valX)<<endl;
  *out<<(cnt+1)*delta<<"\t"<<0<<"\t"<<2<<"\t"<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(valY)<<"\t"<<imag(valY)<<endl;
  *out<<(cnt+1)*delta<<"\t"<<0<<"\t"<<3<<"\t"<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(valZ)<<"\t"<<imag(valZ)<<endl;
  *out<<(cnt+1)*delta<<"\t"<<0<<"\t"<<4<<"\t"<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(valZZ)<<"\t"<<imag(valZZ)<<endl;

  out->close();delete out;
  cout<<"t="<<(cnt+1)*delta<<", <X(t)X>="<<valX/norm<<", <Y(t)Y>="<<valY/norm<<", <Z(t)Z>="<<valZ/norm<<endl;

}

void computeCorrelatorsSameStep(const MPS& mpsB,const MPS& bra,
				const Operator& middleUL,const mwArray& specialOp,const mwArray& specialOpb,
				int cnt,double delta,const string& outfname){
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

  // To compute the correlator, I construct a MPO

  static Operator *lastOpX,*lastOpY,*lastOpZ,*lastOp0;
  static Operator *firstOpX,*firstOpY,*firstOpZ,*firstOp0;
  if(first){
    // The ones at the beginning and end of the chain are different, in this case
    mwArray aux(specialOp); // d2(r) d2(u) d2 (l) d2 (d)
    cout<<"specialOp:"<<specialOp.getDimensions()<<endl;
    aux.reshape(Indices(d2,d2,d2*d2));
    aux.permute(Indices(2,1,3)); // d2(u) x d2(r) d2(l) d2(d)
    aux.reshape(Indices(d2,d2*d2*d2));
    aux.multiplyLeft(reshape(identityMatrix(d),Indices(1,d2))); // 1 x d2(r) d2(l) d2(r)
    aux.reshape(Indices(d2*d2,d,d));
    aux.permute(Indices(1,3,2)); // d2(r)*d2(l)*d'(d) x d(d)
    aux.reshape(Indices(d2*d2*d,d)); 
    firstOpX=new Operator(reshape(permute(reshape(aux*sigX,Indices(d2*d2,d,d)),Indices(1,3,2)),Indices(d2,1,d2,d2)));
    firstOpY=new Operator(reshape(permute(reshape(aux*sigY,Indices(d2*d2,d,d)),Indices(1,3,2)),Indices(d2,1,d2,d2)));
    firstOpZ=new Operator(reshape(permute(reshape(aux*sigZ,Indices(d2*d2,d,d)),Indices(1,3,2)),Indices(d2,1,d2,d2)));
    firstOp0=new Operator(reshape(permute(reshape(aux*sig0,Indices(d2*d2,d,d)),Indices(1,3,2)),Indices(d2,1,d2,d2)));

    aux=specialOpb; // d2 (incl edgeR) x d2(u) x d2 (incl edgeL) x d2
    aux.reshape(Indices(d2*d2*d2,d2));
    aux.multiplyRight(reshape(identityMatrix(d),Indices(d2,1))); // d2(incl edgeR) d2(u) d2(incl edgeL)
    aux.reshape(Indices(d2,d,d*d2));
    aux.permute(Indices(2,1,3)); // d d2(R) d'(u) d2(L)
    aux.reshape(Indices(d,d2*d*d2));
    lastOpX=new Operator(reshape(permute(reshape(sigX*aux,Indices(d,d2,d*d2)),Indices(2,1,3)),Indices(d2,d2,d2,1)));
    lastOpY=new Operator(reshape(permute(reshape(sigY*aux,Indices(d,d2,d*d2)),Indices(2,1,3)),Indices(d2,d2,d2,1)));
    lastOpZ=new Operator(reshape(permute(reshape(sigZ*aux,Indices(d,d2,d*d2)),Indices(2,1,3)),Indices(d2,d2,d2,1)));
    lastOp0=new Operator(reshape(permute(reshape(sig0*aux,Indices(d,d2,d*d2)),Indices(2,1,3)),Indices(d2,d2,d2,1)));
   
    first=false;
  }
  
  // Open the file to write them
  ofstream* out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
    exit(1);
  }
  *out<<setprecision(15);

  Contractor& contractor=Contractor::theContractor();

  //  cout<<"Step after "<<cnt<<" (length of mps:"<<mpsB.getLength()<<")"<<endl;
  MPO link(mpsB.getLength());
  for(int k=1;k<mpsB.getLength()-1;k++)
    link.setOp(k,&middleUL,false); 
  link.setOp(0,firstOp0,false);
  link.setOp(mpsB.getLength()-1,lastOp0,false);

  complex_t norm=contractor.contract(mpsB,link,bra);
  link.setOp(0,firstOpX,false);
  link.setOp(mpsB.getLength()-1,lastOpX,false);
  complex_t valX=contractor.contract(mpsB,link,bra);

  link.setOp(0,firstOpY,false);
  link.setOp(mpsB.getLength()-1,lastOpY,false);
  complex_t valY=contractor.contract(mpsB,link,bra);

  link.setOp(0,firstOpZ,false);
  link.setOp(mpsB.getLength()-1,lastOpZ,false);
  complex_t valZ=contractor.contract(mpsB,link,bra);

  *out<<(cnt)*delta<<"\t"<<0<<"\t"<<1<<"\t"<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(valX)<<"\t"<<imag(valX)<<endl;
  *out<<(cnt)*delta<<"\t"<<0<<"\t"<<2<<"\t"<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(valY)<<"\t"<<imag(valY)<<endl;
  *out<<(cnt)*delta<<"\t"<<0<<"\t"<<3<<"\t"<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(valZ)<<"\t"<<imag(valZ)<<endl;
  // *out<<(cnt)*delta<<"\t"<<0<<"\t";
  // *out<<real(norm)<<"\t"<<imag(norm)<<"\t";
  // *out<<real(valX)<<"\t"<<imag(valX)<<"\t";
  // *out<<real(valY)<<"\t"<<imag(valY)<<"\t";
  // *out<<real(valZ)<<"\t"<<imag(valZ)<<"\t";
  // *out<<endl;

  out->close();delete out;
  cout<<"t="<<(cnt)*delta<<", <X(t)X>="<<valX/norm<<", <Y(t)Y>="<<valY/norm<<", <Z(t)Z>="<<valZ/norm<<endl;

}

void computeCorrelatorsTwoStepsAfter(const MPS& mpsB,const MPS& bra,
				  const Operator& middleUL,const mwArray& specialOp,const mwArray& specialOpb,
				  int cnt,double delta,const string& outfname){
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

  // To compute the correlator, I construct a MPO

  static Operator *lastOpX,*lastOpY,*lastOpZ,*lastOp0;
  static Operator *firstOpX,*firstOpY,*firstOpZ,*firstOp0;
  if(first){
    // The ones at the beginning and end of the chain are different, in this case
    mwArray aux(specialOp); // d2(r) d2(u) d2 (l) d2 (d)
    //cout<<"specialOp:"<<specialOp.getDimensions()<<endl;
    aux.reshape(Indices(d2,d2,d2*d2));
    aux.permute(Indices(2,1,3)); // d2(u) x d2(r) d2(l) d2(d)
    aux.reshape(Indices(d2,d2*d2*d2));

    firstOpX=new Operator(reshape(reshape(permute(sigX,Indices(2,1)),Indices(1,d2))*aux,Indices(d2,1,d2,d2)));
    firstOpY=new Operator(reshape(reshape(permute(sigY,Indices(2,1)),Indices(1,d2))*aux,Indices(d2,1,d2,d2)));
    firstOpZ=new Operator(reshape(reshape(permute(sigZ,Indices(2,1)),Indices(1,d2))*aux,Indices(d2,1,d2,d2)));
    firstOp0=new Operator(reshape(reshape(permute(sig0,Indices(2,1)),Indices(1,d2))*aux,Indices(d2,1,d2,d2)));


    aux=specialOpb; // d2 (incl edgeR) x d2(u) x d2 (incl edgeL) x d2
    aux.reshape(Indices(d2*d2*d2,d2));
    lastOpX=new Operator(reshape(aux*reshape(sigX,Indices(d2,1)),Indices(d2,d2,d2,1)));
    lastOpY=new Operator(reshape(aux*reshape(sigY,Indices(d2,1)),Indices(d2,d2,d2,1)));
    lastOpZ=new Operator(reshape(aux*reshape(sigZ,Indices(d2,1)),Indices(d2,d2,d2,1)));
    lastOp0=new Operator(reshape(aux*reshape(sig0,Indices(d2,1)),Indices(d2,d2,d2,1)));
    
    first=false;
  }
  
  // Open the file to write them
  ofstream* out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
    exit(1);
  }
  *out<<setprecision(15);

  Contractor& contractor=Contractor::theContractor();

  //  cout<<"Step after "<<cnt<<" (length of mps:"<<mpsB.getLength()<<")"<<endl;
  MPO link(mpsB.getLength());
  for(int k=1;k<mpsB.getLength()-1;k++)
    link.setOp(k,&middleUL,false); 
  link.setOp(0,firstOp0,false);
  link.setOp(mpsB.getLength()-1,lastOp0,false);

  complex_t norm=contractor.contract(mpsB,link,bra);
  link.setOp(0,firstOpX,false);
  link.setOp(mpsB.getLength()-1,lastOpX,false);
  complex_t valX=contractor.contract(mpsB,link,bra);

  link.setOp(0,firstOpY,false);
  link.setOp(mpsB.getLength()-1,lastOpY,false);
  complex_t valY=contractor.contract(mpsB,link,bra);

  link.setOp(0,firstOpZ,false);
  link.setOp(mpsB.getLength()-1,lastOpZ,false);
  complex_t valZ=contractor.contract(mpsB,link,bra);

  // I can also compute ZZ for the energy density:
  const MPO* ptrs[]={&link,&link};
  MPO link2(mpsB.getLength());
  MPO::join(2,ptrs,link2);
  complex_t valZZ=contractor.contract(mpsB,link2,bra);

  *out<<(cnt+2)*delta<<"\t"<<0<<"\t"<<1<<"\t"<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(valX)<<"\t"<<imag(valX)<<endl;
  *out<<(cnt+2)*delta<<"\t"<<0<<"\t"<<2<<"\t"<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(valY)<<"\t"<<imag(valY)<<endl;
  *out<<(cnt+2)*delta<<"\t"<<0<<"\t"<<3<<"\t"<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(valZ)<<"\t"<<imag(valZ)<<endl;
  *out<<(cnt+2)*delta<<"\t"<<0<<"\t"<<4<<"\t"<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(valZZ)<<"\t"<<imag(valZZ)<<endl;

  out->close();delete out;
  cout<<"t="<<(cnt+2)*delta<<", <X(t)X>="<<valX/norm<<", <Y(t)Y>="<<valY/norm<<", <Z(t)Z>="<<valZ/norm<<endl;

}


 
void computeDistantCorrelators(const MPS& mpsB,
			       const mwArray& midU,
			       const mwArray& edgeUL,const mwArray& edgeULb,
			       const mwArray& edgeUR,const mwArray& edgeURb,
			       int cnt,double delta,const string& outfname,
			       int Dcut,int maxM,int whichOp){
  //cout<<"computeDistantCorrelators("<<whichOp<<") cnt="<<cnt<<", mps length "<<mpsB.getLength()<<endl;
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
  int xi=midU.getDimension(1);
  
  // To compute the correlator, I construct a MPO
  static Operator *middleU=new Operator(midU,Indices(4,1,2,3)); // xi(r) d(u) xi(l) d(d)
  static Operator *middleUtr=new Operator(conjugate(midU),Indices(2,1,4,3)); // xi(l) d(u) xi(r) d(d)
  static Operator *firstOpX,*firstOpY,*firstOpZ; //*firstOp0; // mid with op an edgeUL to start the contractions (with hook down right) 
  static Operator *lastOpX,*lastOpY,*lastOpZ; //*lastOp0; // mid with op an edgeULb to start the contractions (with hook up right)
  static Operator *firstOp0,*lastOp0; // only for debugging (I think)
  Operator leftmostHook(reshape(identityMatrix(d2),Indices(d2,1,1,d2))); // "hooks" (just identities)
  Operator rightmostHook(reshape(identityMatrix(d2),Indices(d2,d2,1,1))); // would conjugate, no need as it is Id
  Operator leftHookWithEdge(reshape(permute(edgeURb,Indices(2,1)),Indices(1,d2,xi,1)));
  Operator rightHookWithEdge(reshape(conjugate(edgeUR),Indices(1,1,xi,d2))); // for subsequent steps the hook carries the "closing" edge
  Operator edgeOpL(reshape(permute(edgeUL,Indices(2,1)),Indices(xi,1,d2,1)));
  Operator edgeOpR(reshape(permute(edgeURb,Indices(2,1)),Indices(d2,1,xi,1)));
  //cout<<"edgeUL*edgeUR="<<edgeUL*edgeUR<<endl;
  //cout<<" in  edgeOpL:"<<edgeOpL.getFullData()<<endl;exit(1);
  //cout<<"edgeULb"<<edgeULb<<endl;exit(1);
  Operator idOpXi(reshape(identityMatrix(xi),Indices(xi,1,xi,1)));
  Operator traceIdL(reshape(identityMatrix(d),Indices(1,1,d*d,1)));
  Operator traceIdR(reshape(identityMatrix(d),Indices(d*d,1,1,1)));
  Operator secondOpL(reshape(edgeUL*reshape(permute(midU,Indices(2,1,3,4)),
					    Indices(xi,d2*d2*xi)), // d2 x d2*d2*xi
			     Indices(d2,d2,d2,xi)),Indices(4,2,1,3)); // xi x d2 x d2 (incl edge) x d2
  Operator rightsecondOpL(conjugate(permute(reshape(reshape(midU,Indices(d2*xi*d2,d2))*edgeURb,
						    Indices(d2,xi,d2,d2)),Indices(2,1,4,3))));

  //cout<<"Created the constant Ops"<<endl;
  if(first){
    // The ones at the beginning and end of the chain are different, in this case
    mwArray aux(midU); // du xil dd xir
    aux.reshape(Indices(d2*xi*d2,xi));
    aux.multiplyRight(edgeUR);// du xil dd dr
    aux.reshape(Indices(d2,xi*d2*d2));

    mwArray aux0(aux);
    aux0.multiplyLeft(reshape(sig0,Indices(1,d2))); //xil dd dr
    aux0.reshape(Indices(1,xi,d2,d2)); // 1 xil dd dr
    aux0.permute(Indices(2,1,4,3));aux0.conjugate();
    firstOp0=new Operator(aux0);
    aux0=aux;
    aux0.multiplyLeft(reshape(permute(sigX,Indices(2,1)),Indices(1,d2))); //xil dd dr
    aux0.reshape(Indices(1,xi,d2,d2)); // 1 xil dd dr
    aux0.permute(Indices(2,1,4,3));aux0.conjugate();
    firstOpX=new Operator(aux0);
    aux0=aux;
    aux0.multiplyLeft(reshape(permute(sigY,Indices(2,1)),Indices(1,d2))); //xil dd dr
    aux0.reshape(Indices(1,xi,d2,d2)); // 1 xil dd dr
    aux0.permute(Indices(2,1,4,3));aux0.conjugate();
    firstOpY=new Operator(aux0);
    aux0=aux;
    aux0.multiplyLeft(reshape(permute(sigZ,Indices(2,1)),Indices(1,d2))); //xil dd dr
    aux0.reshape(Indices(1,xi,d2,d2)); // 1 xil dd dr
    aux0.permute(Indices(2,1,4,3));aux0.conjugate();
    firstOpZ=new Operator(aux0);

    aux=midU; // d2(u) xi(l) d2(d) xi(r)
    aux.permute(Indices(2,1,3,4)); // xi(l) d2(u) d2(d) xi(r)
    aux.reshape(Indices(xi,d2*d2*xi));
    aux.multiplyLeft(edgeULb); // d2(l) x d2(u)*d2(d)*xi(r)
    aux.reshape(Indices(d2,d2,d2,xi)); //d2(l) d2(u) d2(d) xi(r)
    aux.permute(Indices(4,2,1,3)); // xi(r) d2(u) d2(l) d2(d)
    aux.reshape(Indices(xi*d2*d2,d2));
    lastOpX=new Operator(reshape(aux*reshape(sigX,Indices(d2,1)),Indices(xi,d2,d2,1)));
    lastOpY=new Operator(reshape(aux*reshape(sigY,Indices(d2,1)),Indices(xi,d2,d2,1)));
    lastOpZ=new Operator(reshape(aux*reshape(sigZ,Indices(d2,1)),Indices(xi,d2,d2,1)));
    lastOp0=new Operator(reshape(aux*reshape(sig0,Indices(d2,1)),Indices(xi,d2,d2,1)));

    first=false;
    //cout<<"Initialized the static Ops"<<endl;
  }

  Contractor& contractor=Contractor::theContractor();

  int len=mpsB.getLength();

  // Initialize the MPS for left and right boundaries
  MPS mpsLop(mpsB); // operator at the bottom (left boundary)=> extra site at the beginning
  MPS mpsRop(mpsB); // operator at the top (right boundary)=> extra site at the end
  mpsRop.conjugateMPS();
  MPS::gaugeCond(mpsLop,mpsRop,'L');
  // I am not sure that this helps, anyway!
  double normFL=0.;
  double normFR=0.; // to keep normalization under control, as I expect it to go as d^(l+1)
  // First of all, get the (common) norm factor (then for each term tehre is also a d^(l+1) factor
  complex_t norm0;
  {// the MPO is shorter than later, and it needs both edges (R and L) because this case is special
    MPO mpo(len);
    mpo.setOp(0,new Operator(reshape(permute(edgeUL*edgeUR,Indices(2,1)),Indices(d2,1,d2,1))),true);
    mwArray id2=identityMatrix(d);id2.reshape(Indices(d2,1));
    mpo.setOp(len-1,new Operator(reshape(permute(edgeULb*edgeURb,Indices(2,1)),Indices(d2,1,d2,1))),true);
    for(int k=1;k<len-1;k++)
      mpo.setOp(k,&idOpXi,false);
    norm0=contractor.contract(mpsLop,mpo,mpsRop);
  }  
  //cout<<"norm0="<<norm0<<endl;

  // Now apply the operators at the edges
  mpsLop.insertSite(0,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);
  mpsRop.insertSite(len,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);

    // prepare the MPOs for the contraction
  MPO mpoEdges(len+1);
  mpoEdges.setOp(0,&edgeOpL,false);
  mpoEdges.setOp(len,&edgeOpR,false);
  for(int k=1;k<len;k++) mpoEdges.setOp(k,&idOpXi,false);
  MPO mpoTraces(len+2);
  mpoTraces.setOp(0,&traceIdL,false);
  mpoTraces.setOp(len+1,&traceIdR,false);
  for(int k=1;k<len+1;k++) mpoTraces.setOp(k,&idOpXi,false);
  //cout<<"Created MPOs for contraction MPOedges: "<<mpoEdges<<" and mpoTraces:"<<mpoTraces<<endl;

  { // since the MPOs with the op are only needed once, I close a scope here
    // Modify those MPSs to include the operators at the edges  
    MPO oplow(len+1); // the one with the operator below
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    oplow.setOp(0,&leftmostHook,false);
    oplow.setOp(1,&secondOpL,false);
    for(int k=2;k<len;k++)
      oplow.setOp(k,middleU,false);
    MPO ophigh(len+1); // the one with the operator above
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    ophigh.setOp(len,&rightmostHook,false);
    ophigh.setOp(len-1,&rightsecondOpL,false);
    for(int k=1;k<len-1;k++)
      ophigh.setOp(k,middleUtr,false);
    switch(whichOp){
    case 0:
      ophigh.setOp(0,firstOp0,false); 
      oplow.setOp(len,lastOp0,false); 
      break;
    case 1:
      ophigh.setOp(0,firstOpX,false); 
      oplow.setOp(len,lastOpX,false); 
      break;
    case 2:
      ophigh.setOp(0,firstOpY,false); 
      oplow.setOp(len,lastOpY,false); 
      break;
    case 3:
      ophigh.setOp(0,firstOpZ,false); 
      oplow.setOp(len,lastOpZ,false); 
      break;
    }
    MPS aux(mpsLop);
    contractor.optimize(oplow,aux,mpsLop,Dcut);
    mpsLop.gaugeCond('L',0);
    normFL+=log2(mpsLop.getNormFact());
    mpsLop.setNormFact(1.);
    aux=mpsRop;
    contractor.optimize(ophigh,aux,mpsRop,Dcut);
    mpsRop.gaugeCond('L',0);
    normFR+=log2(mpsRop.getNormFact());
    mpsRop.setNormFact(1.);
  }
  // I need to keep the normalization factors!
  //cout<<"After applying the operator columns, length of mps:"<<mpsLop.getLength()<<endl;
  
  int cntT=cnt+2;
  int l=1; // distance
  // Now iterate, producing contractions and evolving once more

  // Prepare MPOs for the "evolution"
  MPO oplow(len+2); // the one for the left term with the operator below(now no operator, though)
  oplow.setOp(0,&leftmostHook,false);
  oplow.setOp(1,&secondOpL,false);
  for(int k=2;k<len+1;k++)
    oplow.setOp(k,middleU,false);
  oplow.setOp(len+1,&leftHookWithEdge,false);
  MPO ophigh(len+2); // the one with the operator above
  // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
  ophigh.setOp(len+1,&rightmostHook,false);
  ophigh.setOp(len,&rightsecondOpL,false);
  for(int k=1;k<len;k++)
    ophigh.setOp(k,middleUtr,false);
  ophigh.setOp(0,&rightHookWithEdge,false);

  //cout<<"Created MPOs for evol oplow: "<<oplow<<" and ophigh:"<<ophigh<<endl;
  // Open the file to write them
  ofstream* out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
    exit(1);
  }
  
  while(cntT<=maxM){
    //cout<<"Starting contractions within step "<<cnt<<" now cntT "<<cntT<<", length of mps "<<mpsRop.getLength()<<endl;
    // first, contract the left with the right with traced out edges (t+2, i.e. cntT, l=1)
    complex_t val0_1=contractWithTracedEdges(mpsLop,mpsRop,mpoTraces);
    val0_1*=pow(2,normFL+normFR-(l+1)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<whichOp<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val0_1)<<"\t"<<imag(val0_1)<<endl;
    //cout<<"CC("<<whichOp<<")[t="<<cntT<<",l="<<l<<"]="<<val0_1/norm0<<endl;
    // then contract both with the operator with edges (t+3,l=1)
    cntT++;
    complex_t val1_1=contractor.contract(mpsLop,mpoEdges,mpsRop);
    val1_1*=pow(2,normFL+normFR-(l+1)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<whichOp<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val1_1)<<"\t"<<imag(val1_1)<<endl;
    //cout<<"CC("<<whichOp<<")[t="<<cntT<<",l="<<l<<"]="<<val1_1/norm0<<endl;
    // now evolve the left (op below)
    applyOneStepLeft(mpsLop,oplow,Dcut);
    mpsLop.gaugeCond('L',0);
    normFL+=log2(mpsLop.getNormFact());
    mpsLop.setNormFact(1.);
    // contract with previous right with traced out (t+3, l=2)
    l++;
    complex_t val1_2=contractWithTracedEdges(mpsLop,mpsRop,mpoTraces);    
    val1_2*=pow(2,normFL+normFR-(l+1)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<whichOp<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val1_2)<<"\t"<<imag(val1_2)<<endl;
    // contract with previous right with edges (t+4, l=2)
    cntT++;
    complex_t val2_2=contractor.contract(mpsLop,mpoEdges,mpsRop);
    val2_2*=pow(2,normFL+normFR-(l+1)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<whichOp<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val2_2)<<"\t"<<imag(val2_2)<<endl;
    // evolve right (and I will be at cnt+=2, l+=2)
    applyOneStepRight(mpsRop,ophigh,Dcut);
    mpsRop.gaugeCond('L',0);
    normFR+=log2(mpsRop.getNormFact());
    mpsRop.setNormFact(1.);
    l++;
  }
  out->close();
  delete out;
}

// auxiliary to apply one step of evolution to the left vector
// (requires creating an auxiliary site of dim 1 at the beginning and
// at the end absorbing the last site (which should also be dim 1
void applyOneStepLeft(MPS& mps,const MPO& oplow,int Dcut){
  int len=mps.getLength();
  Contractor& contractor=Contractor::theContractor();
  MPS aux(mps);
  aux.insertSite(0,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);
  mps.insertSite(len,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);
  contractor.optimize(oplow,aux,mps,Dcut);
  // now absorb the last site in the previous one and copy everything to the original mps
  // I use the trick of "tracing out" the site
  vector<int> pos(1,len);
  mps.traceOutSites(pos,false); // do not normalize again
}

void applyOneStepRight(MPS& mps,const MPO& ophigh,int Dcut){
  int len=mps.getLength();
  Contractor& contractor=Contractor::theContractor();
  MPS aux(mps);
  aux.insertSite(len,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);// empty site at the end  
  mps.insertSite(0,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);//  For a better starting point, I "move" the sites down
  contractor.optimize(ophigh,aux,mps,Dcut);
  // now absorb the last site in the previous one and copy everything to the original mps
  // I use the trick of "tracing out" the site
  vector<int> pos(1,0);
  mps.traceOutSites(pos,false); // do not normalize again
}


// auxiliary, to dot he manipulations required to contract with the links shifted and edges traced out
complex_t contractWithTracedEdges(MPS& left,MPS& right,const MPO& mpoTraces){
  // for the left I create a new MPS with a1 site at the right edge, for the right the extra site is at the beginning
  Contractor& contractor=Contractor::theContractor();
  int len=left.getLength();  
  left.insertSite(len,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);// empty site at the end
  right.insertSite(0,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);// empty site at the beginning
  complex_t value=contractor.contract(left,mpoTraces,right);  
  vector<int> pos(1,0);
  right.traceOutSites(pos,false); // do not normalize again
  pos[0]=len;
  left.traceOutSites(pos,false); // do not normalize again
  return value;
}


void computeDiagonalCorrelators(const mwArray& edgeUL,const mwArray& edgeURb,
				double delta,const string& outfname,int maxM){
  int d2=d*d;
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1*I_c,ZERO_c};
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c}; //Z
  mwArray sig0=identityMatrix(d);
  mwArray sigX=mwArray(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray step(edgeUL); //d2 x xi
  step.multiplyRight(.5*edgeURb); // d2 x d2
  //  putForMatlab(cout,step,"stepM");
  // cancelling one factor d per l
  mwArray stepL(step);
  int cnt=1;
  ofstream* out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
    exit(1);
  }
  while(cnt<maxM){
    complex_t val0=(reshape(sig0,Indices(1,d2))*stepL
		    *reshape(permute(sig0,Indices(2,1)),Indices(d2,1))).getElement(0);
    complex_t valX=(reshape(sigX,Indices(1,d2))*stepL
		    *reshape(permute(sigX,Indices(2,1)),Indices(d2,1))).getElement(0);
    complex_t valY=(reshape(sigY,Indices(1,d2))*stepL
		    *reshape(permute(sigY,Indices(2,1)),Indices(d2,1))).getElement(0);
    complex_t valZ=(reshape(sigZ,Indices(1,d2))*stepL
		    *reshape(permute(sigZ,Indices(2,1)),Indices(d2,1))).getElement(0);
    *out<<cnt*delta<<"\t"<<cnt<<"\t"<<1<<"\t"<<real(val0)<<"\t"<<imag(val0)<<"\t"<<real(valX)<<"\t"<<imag(valX)<<endl;
    *out<<cnt*delta<<"\t"<<cnt<<"\t"<<2<<"\t"<<real(val0)<<"\t"<<imag(val0)<<"\t"<<real(valY)<<"\t"<<imag(valY)<<endl;
    *out<<cnt*delta<<"\t"<<cnt<<"\t"<<3<<"\t"<<real(val0)<<"\t"<<imag(val0)<<"\t"<<real(valZ)<<"\t"<<imag(valZ)<<endl;
    stepL.multiplyRight(step);
    cnt++;
  }
  out->close();delete out;
}

void computeDistantCorrelatorsTwoBody(const MPS& mpsB,
				      const mwArray& midU,
				      const mwArray& edgeUL,const mwArray& edgeULb,
				      const mwArray& edgeUR,const mwArray& edgeURb,
				      int cnt,double delta,const string& outfname,
				      int Dcut,int maxM,int whichOp1u,int whichOp2u,int whichOp1b,int whichOp2b){
  cout<<"computeDistantCorrelatorsTwoBody WARNING!!! Not yet properly tested!!!"<<endl;
  //cout<<"computeDistantCorrelators("<<whichOp<<") cnt="<<cnt<<", mps length "<<mpsB.getLength()<<endl;
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
  int xi=midU.getDimension(1);
  
  // To compute the correlator, I construct a MPO
  static Operator *middleU=new Operator(midU,Indices(4,1,2,3)); // xi(r) d(u) xi(l) d(d)
  static Operator *middleUtr=new Operator(conjugate(midU),Indices(2,1,4,3)); // xi(l) d(u) xi(r) d(d)
  static Operator *firstOpX,*firstOpY,*firstOpZ; //*firstOp0; // mid with op an edgeUL to start the contractions (with hook down right) 
  static Operator *lastOpX,*lastOpY,*lastOpZ; //*lastOp0; // mid with op an edgeULb to start the contractions (with hook up right)
  static Operator *firstOp0,*lastOp0; // only for debugging (I think)
  static Operator *firstOpXnoE,*firstOpYnoE,*firstOpZnoE; // mid with op but no edge (for the second in the pair
  static Operator *lastOpXnoE,*lastOpYnoE,*lastOpZnoE; // mid with op but no edge (for the second in the pair
  Operator leftmostHook(reshape(identityMatrix(d2),Indices(d2,1,1,d2))); // "hooks" (just identities)
  Operator rightmostHook(reshape(identityMatrix(d2),Indices(d2,d2,1,1))); // would conjugate, no need as it is Id
  Operator leftHookWithEdge(reshape(permute(edgeURb,Indices(2,1)),Indices(1,d2,xi,1)));
  Operator rightHookWithEdge(reshape(conjugate(edgeUR),Indices(1,1,xi,d2))); // for subsequent steps the hook carries the "closing" edge
  Operator edgeOpL(reshape(permute(edgeUL,Indices(2,1)),Indices(xi,1,d2,1)));
  Operator edgeOpR(reshape(permute(edgeURb,Indices(2,1)),Indices(d2,1,xi,1)));
  //cout<<"edgeUL*edgeUR="<<edgeUL*edgeUR<<endl;
  //cout<<" in  edgeOpL:"<<edgeOpL.getFullData()<<endl;exit(1);
  //cout<<"edgeULb"<<edgeULb<<endl;exit(1);
  Operator idOpXi(reshape(identityMatrix(xi),Indices(xi,1,xi,1)));
  Operator traceIdL(reshape(identityMatrix(d),Indices(1,1,d*d,1)));
  Operator traceIdR(reshape(identityMatrix(d),Indices(d*d,1,1,1)));
  Operator secondOpL(reshape(edgeUL*reshape(permute(midU,Indices(2,1,3,4)),
					    Indices(xi,d2*d2*xi)), // d2 x d2*d2*xi
			     Indices(d2,d2,d2,xi)),Indices(4,2,1,3)); // xi x d2 x d2 (incl edge) x d2
  Operator rightsecondOpL(conjugate(permute(reshape(reshape(midU,Indices(d2*xi*d2,d2))*edgeURb,
						    Indices(d2,xi,d2,d2)),Indices(2,1,4,3))));

  //cout<<"Created the constant Ops"<<endl;
  if(first){
    // The ones at the beginning and end of the chain are different, in this case
    mwArray aux(midU); // du xil dd xir
    aux.reshape(Indices(d2*xi*d2,xi));
    aux.multiplyRight(edgeUR);// du xil dd dr
    aux.reshape(Indices(d2,xi*d2*d2));

    mwArray aux0(aux);
    aux0.multiplyLeft(reshape(sig0,Indices(1,d2))); //xil dd dr
    aux0.reshape(Indices(1,xi,d2,d2)); // 1 xil dd dr
    aux0.permute(Indices(2,1,4,3));aux0.conjugate();
    firstOp0=new Operator(aux0);
    aux0=aux;
    aux0.multiplyLeft(reshape(permute(sigX,Indices(2,1)),Indices(1,d2))); //xil dd dr
    aux0.reshape(Indices(1,xi,d2,d2)); // 1 xil dd dr
    aux0.permute(Indices(2,1,4,3));aux0.conjugate();
    firstOpX=new Operator(aux0);
    aux0=aux;
    aux0.multiplyLeft(reshape(permute(sigY,Indices(2,1)),Indices(1,d2))); //xil dd dr
    aux0.reshape(Indices(1,xi,d2,d2)); // 1 xil dd dr
    aux0.permute(Indices(2,1,4,3));aux0.conjugate();
    firstOpY=new Operator(aux0);
    aux0=aux;
    aux0.multiplyLeft(reshape(permute(sigZ,Indices(2,1)),Indices(1,d2))); //xil dd dr
    aux0.reshape(Indices(1,xi,d2,d2)); // 1 xil dd dr
    aux0.permute(Indices(2,1,4,3));aux0.conjugate();
    firstOpZ=new Operator(aux0);

    aux=midU; // d2(u) xi(l) d2(d) xi(r)
    aux.permute(Indices(2,1,3,4)); // xi(l) d2(u) d2(d) xi(r)
    aux.reshape(Indices(xi,d2*d2*xi));
    aux.multiplyLeft(edgeULb); // d2(l) x d2(u)*d2(d)*xi(r)
    aux.reshape(Indices(d2,d2,d2,xi)); //d2(l) d2(u) d2(d) xi(r)
    aux.permute(Indices(4,2,1,3)); // xi(r) d2(u) d2(l) d2(d)
    aux.reshape(Indices(xi*d2*d2,d2));
    lastOpX=new Operator(reshape(aux*reshape(sigX,Indices(d2,1)),Indices(xi,d2,d2,1)));
    lastOpY=new Operator(reshape(aux*reshape(sigY,Indices(d2,1)),Indices(xi,d2,d2,1)));
    lastOpZ=new Operator(reshape(aux*reshape(sigZ,Indices(d2,1)),Indices(xi,d2,d2,1)));
    lastOp0=new Operator(reshape(aux*reshape(sig0,Indices(d2,1)),Indices(xi,d2,d2,1)));

    // also special edges for the second in the pair
    aux=midU; // d2(u) xi(l) d2(d) xi(r)
    aux.reshape(Indices(d2,xi*d2*xi));
    firstOpXnoE=new Operator(reshape(reshape(permute(sigX,Indices(2,1)),Indices(1,d2))*aux,Indices(1,xi,d2,xi)),Indices(4,1,2,3));
    firstOpYnoE=new Operator(reshape(reshape(permute(sigY,Indices(2,1)),Indices(1,d2))*aux,Indices(1,xi,d2,xi)),Indices(4,1,2,3));
    firstOpZnoE=new Operator(reshape(reshape(permute(sigZ,Indices(2,1)),Indices(1,d2))*aux,Indices(1,xi,d2,xi)),Indices(4,1,2,3));

    aux.reshape(Indices(d2*xi,d2,xi)); // d2(u)*xi(l) d2(d) xi(r)
    aux.permute(Indices(1,3,2)); // d2(u)*xi(l)*xi(r) d2(d) 
    aux.reshape(Indices(d2*xi*xi,d2));
    lastOpXnoE=new Operator(reshape(aux*reshape(sigX,Indices(d2,1)),Indices(d2,xi,xi,1)),Indices(3,1,2,4));
    lastOpYnoE=new Operator(reshape(aux*reshape(sigY,Indices(d2,1)),Indices(d2,xi,xi,1)),Indices(3,1,2,4));
    lastOpZnoE=new Operator(reshape(aux*reshape(sigZ,Indices(d2,1)),Indices(d2,xi,xi,1)),Indices(3,1,2,4));
    
    first=false;
    //cout<<"Initialized the static Ops"<<endl;
  }

  // For two sites I need a couple of special operators
  int cntT=cnt+2;
  int l=1; // distance
  // Open the file to write them
  ofstream* out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
    exit(1);
  }


  Contractor& contractor=Contractor::theContractor();

  int len=mpsB.getLength();

  // Initialize the MPS for left and right boundaries
  MPS mpsLop(mpsB); // operator at the bottom (left boundary)=> extra site at the beginning
  MPS mpsRop(mpsB); // operator at the top (right boundary)=> extra site at the end
  mpsRop.conjugateMPS();
  MPS::gaugeCond(mpsLop,mpsRop,'L');
  // I am not sure that this helps, anyway!
  double normFL=0.;
  double normFR=0.; // to keep normalization under control, as I expect it to go as d^(l+1)
  // First of all, get the (common) norm factor (then for each term there is also a d^(l+1) factor
  complex_t norm0;
  {// the MPO is shorter than later, and it needs both edges (R and L) because this case is special
    MPO mpo(len);
    mpo.setOp(0,new Operator(reshape(permute(edgeUL*edgeUR,Indices(2,1)),Indices(d2,1,d2,1))),true);
    mwArray id2=identityMatrix(d);id2.reshape(Indices(d2,1));
    mpo.setOp(len-1,new Operator(reshape(permute(edgeULb*edgeURb,Indices(2,1)),Indices(d2,1,d2,1))),true);
    for(int k=1;k<len-1;k++)
      mpo.setOp(k,&idOpXi,false);
    norm0=contractor.contract(mpsLop,mpo,mpsRop);
  }  
  //cout<<"norm0="<<norm0<<endl;

  // Now apply the operators at the edges
  mpsLop.insertSite(0,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);
  mpsRop.insertSite(len,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);

    // prepare the MPOs for the contraction; in this case, tehy are one site longer each
  MPO mpoEdges(len+2);
  mpoEdges.setOp(0,&edgeOpL,false);
  mpoEdges.setOp(len+1,&edgeOpR,false);
  for(int k=1;k<len+1;k++) mpoEdges.setOp(k,&idOpXi,false);
  MPO mpoTraces(len+3);
  mpoTraces.setOp(0,&traceIdL,false);
  mpoTraces.setOp(len+2,&traceIdR,false);
  for(int k=1;k<len+2;k++) mpoTraces.setOp(k,&idOpXi,false);
  //cout<<"Created MPOs for contraction MPOedges: "<<mpoEdges<<" and mpoTraces:"<<mpoTraces<<endl;

  { // since the MPOs with the op are only needed once, I close a scope here
    // Modify those MPSs to include the operators at the edges  
    MPO oplow(len+1); // the one with the operator below
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    oplow.setOp(0,&leftmostHook,false);
    oplow.setOp(1,&secondOpL,false);
    for(int k=2;k<len;k++)
      oplow.setOp(k,middleU,false);
    MPO ophigh(len+1); // the one with the operator above
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    ophigh.setOp(len,&rightmostHook,false);
    ophigh.setOp(len-1,&rightsecondOpL,false);
    for(int k=1;k<len-1;k++)
      ophigh.setOp(k,middleUtr,false);
    switch(whichOp1b){
    case 0:
      oplow.setOp(len,lastOp0,false); 
      break;
    case 1:
      oplow.setOp(len,lastOpX,false); 
      break;
    case 2:
      oplow.setOp(len,lastOpY,false); 
      break;
    case 3:
      oplow.setOp(len,lastOpZ,false); 
      break;
    }
    switch(whichOp2u){
    case 0:
      ophigh.setOp(0,firstOp0,false); 
      break;
    case 1:
      ophigh.setOp(0,firstOpX,false); 
      break;
    case 2:
      ophigh.setOp(0,firstOpY,false); 
      break;
    case 3:
      ophigh.setOp(0,firstOpZ,false); 
      break;
    }
    MPS aux(mpsLop);
    contractor.optimize(oplow,aux,mpsLop,Dcut);
    mpsLop.gaugeCond('L',0);
    normFL+=log2(mpsLop.getNormFact());
    mpsLop.setNormFact(1.);
    aux=mpsRop;
    contractor.optimize(ophigh,aux,mpsRop,Dcut);
    mpsRop.gaugeCond('L',0);
    normFR+=log2(mpsRop.getNormFact());
    mpsRop.setNormFact(1.);
  }
  {
    // In this case, the first term, l=1, corresponds to inserting a
    // special MPO, so I should do this, then evolve again with
    // special (second op) and then the systematic part

    // First the one with the traces (t+2,l=1)
    MPO step1a(len+2);    
    switch(whichOp1u){
    case 1:
      step1a.setOp(1,firstOpXnoE,false);
      break;
    case 2:
      step1a.setOp(1,firstOpYnoE,false);
      break;
    case 3:
      step1a.setOp(1,firstOpZnoE,false);
      break;
    }
    switch(whichOp2b){
    case 1:
      step1a.setOp(len,lastOpXnoE,false);
      break;
    case 2:
      step1a.setOp(len,lastOpYnoE,false);
      break;
    case 3:
      step1a.setOp(len,lastOpZnoE,false);
      break;
    }
    for(int k=2;k<len;k++)step1a.setOp(k,&idOpXi,false);
    step1a.setOp(0,&traceIdL,false);
    step1a.setOp(len+1,&traceIdR,false);
    complex_t val0_1=contractWithTracedEdges(mpsLop,mpsRop,step1a);
    *out<<cntT*delta<<"\t"<<l<<"\t"<<whichOp1b<<"\t"<<whichOp2b<<"\t"<<whichOp1u<<"\t"<<whichOp2u<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val0_1)<<"\t"<<imag(val0_1)<<endl;
    cntT++;
    // Then the one with edges (t+3, l=1)
    MPO step1(len+1);    
    switch(whichOp1u){
    case 1:
      step1.setOp(0,firstOpXnoE,false);
      break;
    case 2:
      step1.setOp(0,firstOpYnoE,false);
      break;
    case 3:
      step1.setOp(0,firstOpZnoE,false);
      break;
    }
    switch(whichOp2b){
    case 1:
      step1.setOp(len,lastOpXnoE,false);
      break;
    case 2:
      step1.setOp(len,lastOpYnoE,false);
      break;
    case 3:
      step1.setOp(len,lastOpZnoE,false);
      break;
    }
    for(int k=1;k<len;k++)step1.setOp(k,&idOpXi,false);

    // And i will use extra mpos for the contraction
    MPO aux1(len+1);
    for(int k=1;k<len+1;k++) aux1.setOp(k,&idOpXi,false);
    aux1.setOp(0,&edgeOpL,false);
    MPO aux2(len+1);
    for(int k=0;k<len;k++) aux2.setOp(k,&idOpXi,false);
    aux2.setOp(len,&edgeOpR,false);

    const MPO* ptrs[]={&aux2,&step1,&aux1};
    MPO _mpo(len+1);
    MPO::join(3,ptrs,_mpo);
    complex_t val1_1=contractor.contract(mpsLop,_mpo,mpsRop);
    *out<<cntT*delta<<"\t"<<l<<"\t"<<whichOp1b<<"\t"<<whichOp2b<<"\t"<<whichOp1u<<"\t"<<whichOp2u<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val1_1)<<"\t"<<imag(val1_1)<<endl;
    l++;
  }
  mpsLop.insertSite(0,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);
  mpsRop.insertSite(len,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);
  {
    // And now evolve them with the secoond op in the pair ****
    //** construct the left and right evolution ops that increase the length by one.
    MPO oplow(len+2); // the one with the operator below
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    oplow.setOp(0,&leftmostHook,false);
    oplow.setOp(1,&secondOpL,false);
    for(int k=2;k<len+1;k++)
      oplow.setOp(k,middleU,false);
    MPO ophigh(len+2); // the one with the operator above
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    ophigh.setOp(len+1,&rightmostHook,false);
    ophigh.setOp(len,&rightsecondOpL,false);
    for(int k=1;k<len;k++)
      ophigh.setOp(k,middleUtr,false);
    switch(whichOp2b){
    case 1:
      oplow.setOp(len,lastOpXnoE,false); 
      break;
    case 2:
      oplow.setOp(len,lastOpYnoE,false); 
      break;
    case 3:
      oplow.setOp(len,lastOpZnoE,false); 
      break;
    }
    switch(whichOp1u){
    case 1:
      ophigh.setOp(0,firstOpXnoE,false); 
      break;
    case 2:
      ophigh.setOp(0,firstOpYnoE,false); 
      break;
    case 3:
      ophigh.setOp(0,firstOpZnoE,false); 
      break;
    }
    MPS aux(mpsLop);
    contractor.optimize(oplow,aux,mpsLop,Dcut);
    mpsLop.gaugeCond('L',0);
    normFL+=log2(mpsLop.getNormFact());
    mpsLop.setNormFact(1.);
    aux=mpsRop;
    contractor.optimize(ophigh,aux,mpsRop,Dcut);
    mpsRop.gaugeCond('L',0);
    normFR+=log2(mpsRop.getNormFact());
    mpsRop.setNormFact(1.);
    
  }

  //cout<<"After applying the operator columns, length of mps:"<<mpsLop.getLength()<<endl;
  
  // Now iterate, producing contractions and evolving once more

  // Prepare MPOs for the "evolution"
  MPO oplow(len+2); // the one for the left term with the operator below(now no operator, though)
  oplow.setOp(0,&leftmostHook,false);
  oplow.setOp(1,&secondOpL,false);
  for(int k=2;k<len+1;k++)
    oplow.setOp(k,middleU,false);
  oplow.setOp(len+1,&leftHookWithEdge,false);
  MPO ophigh(len+2); // the one with the operator above
  // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
  ophigh.setOp(len+1,&rightmostHook,false);
  ophigh.setOp(len,&rightsecondOpL,false);
  for(int k=1;k<len;k++)
    ophigh.setOp(k,middleUtr,false);
  ophigh.setOp(0,&rightHookWithEdge,false);

  //cout<<"Created MPOs for evol oplow: "<<oplow<<" and ophigh:"<<ophigh<<endl;
  
  while(cntT<=maxM){
    //cout<<"Starting contractions within step "<<cnt<<" now cntT "<<cntT<<", length of mps "<<mpsRop.getLength()<<endl;
    // first, contract the left with the right with traced out edges (t+2, i.e. cntT, l=1)
    complex_t val0_1=contractWithTracedEdges(mpsLop,mpsRop,mpoTraces);
    val0_1*=pow(2,normFL+normFR-(l+1)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<whichOp1b<<"\t"<<whichOp2b<<"\t"<<whichOp1u<<"\t"<<whichOp2u<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val0_1)<<"\t"<<imag(val0_1)<<endl;
    //cout<<"CC("<<whichOp<<")[t="<<cntT<<",l="<<l<<"]="<<val0_1/norm0<<endl;
    // then contract both with the operator with edges (t+3,l=1)
    cntT++;
    complex_t val1_1=contractor.contract(mpsLop,mpoEdges,mpsRop);
    val1_1*=pow(2,normFL+normFR-(l+1)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<whichOp1b<<"\t"<<whichOp2b<<"\t"<<whichOp1u<<"\t"<<whichOp2u<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val1_1)<<"\t"<<imag(val1_1)<<endl;
    //cout<<"CC("<<whichOp<<")[t="<<cntT<<",l="<<l<<"]="<<val1_1/norm0<<endl;
    // now evolve the left (op below)
    applyOneStepLeft(mpsLop,oplow,Dcut);
    mpsLop.gaugeCond('L',0);
    normFL+=log2(mpsLop.getNormFact());
    mpsLop.setNormFact(1.);
    // contract with previous right with traced out (t+3, l=2)
    l++;
    complex_t val1_2=contractWithTracedEdges(mpsLop,mpsRop,mpoTraces);    
    val1_2*=pow(2,normFL+normFR-(l+1)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<whichOp1b<<"\t"<<whichOp2b<<"\t"<<whichOp1u<<"\t"<<whichOp2u<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val1_2)<<"\t"<<imag(val1_2)<<endl;
    // contract with previous right with edges (t+4, l=2)
    cntT++;
    complex_t val2_2=contractor.contract(mpsLop,mpoEdges,mpsRop);
    val2_2*=pow(2,normFL+normFR-(l+1)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<whichOp1b<<"\t"<<whichOp2b<<"\t"<<whichOp1u<<"\t"<<whichOp2u<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val2_2)<<"\t"<<imag(val2_2)<<endl;
    // evolve right (and I will be at cnt+=2, l+=2)
    applyOneStepRight(mpsRop,ophigh,Dcut);
    mpsRop.gaugeCond('L',0);
    normFR+=log2(mpsRop.getNormFact());
    mpsRop.setNormFact(1.);
    l++;
  }
  out->close();
  delete out;
}

void computeDistantCorrelatorsTwoBody(const MPS& mpsB,
				      const mwArray& midU,
				      const mwArray& edgeUL,const mwArray& edgeULb,
				      const mwArray& edgeUR,const mwArray& edgeURb,
				      int cnt,double delta,const string& outfname,
				      int Dcut,int maxM,const MPO& mpoLow,const MPO& mpoUp,int label){
  int nops=mpoLow.getLength();
  if(nops!=2||mpoUp.getLength()!=nops){
    cout<<"ERROR: For now, correlators of MPO can only be computed for 2-site case on both sides"<<endl;
    exit(1);
  }

  int d2=d*d;
  int xi=midU.getDimension(1);
  
  // To compute the correlator, I construct a MPO
  Operator middleU(midU,Indices(4,1,2,3)); // xi(r) d(u) xi(l) d(d)
  Operator middleUtr(conjugate(midU),Indices(2,1,4,3)); // xi(l) d(u) xi(r) d(d)
  Operator leftmostHook(reshape(identityMatrix(d2),Indices(d2,1,1,d2))); // "hooks" (just identities)
  Operator rightmostHook(reshape(identityMatrix(d2),Indices(d2,d2,1,1))); // would conjugate, no need as it is Id
  Operator leftHookWithEdge(reshape(permute(edgeURb,Indices(2,1)),Indices(1,d2,xi,1)));
  Operator rightHookWithEdge(reshape(conjugate(edgeUR),Indices(1,1,xi,d2))); // for subsequent steps the hook carries the "closing" edge
  Operator edgeOpL(reshape(permute(edgeUL,Indices(2,1)),Indices(xi,1,d2,1)));
  Operator edgeOpR(reshape(permute(edgeURb,Indices(2,1)),Indices(d2,1,xi,1)));
  //cout<<"edgeUL*edgeUR="<<edgeUL*edgeUR<<endl;
  //cout<<" in  edgeOpL:"<<edgeOpL.getFullData()<<endl;exit(1);
  //cout<<"edgeULb"<<edgeULb<<endl;exit(1);
  Operator idOpXi(reshape(identityMatrix(xi),Indices(xi,1,xi,1)));
  Operator traceIdL(reshape(identityMatrix(d),Indices(1,1,d*d,1)));
  Operator traceIdR(reshape(identityMatrix(d),Indices(d*d,1,1,1)));
  Operator secondOpL(reshape(edgeUL*reshape(permute(midU,Indices(2,1,3,4)),
					    Indices(xi,d2*d2*xi)), // d2 x d2*d2*xi
			     Indices(d2,d2,d2,xi)),Indices(4,2,1,3)); // xi x d2 x d2 (incl edge) x d2
  Operator rightsecondOpL(conjugate(permute(reshape(reshape(midU,Indices(d2*xi*d2,d2))*edgeURb,
						    Indices(d2,xi,d2,d2)),Indices(2,1,4,3))));

  // Now the operators that depend on the MPOs
  // 1) for the upper part (from mpoUp), which get contracted witht eh right side, but in conjugate way
  mwArray opUp1=mpoUp.getOp(0).getFullData(); // d 1 d alpha
  int alphaU=opUp1.getDimension(3);
  opUp1.permute(Indices(4,3,1,2)); // alpha d(d)d(u) 1
  opUp1.reshape(Indices(alphaU,d2));
  mwArray opUp2=mpoUp.getOp(1).getFullData(); // d alpha d 1
  opUp2.permute(Indices(2,3,1,4)); // alpha d(d) d(u) 1
  opUp2.reshape(Indices(alphaU,d2));

  mwArray aux(midU);// du xil dd xir
  aux.reshape(Indices(d2,xi*d2*xi));
  aux.multiplyLeft(opUp1);
  aux.reshape(Indices(alphaU,xi,d2,xi)); // alpha xil dd xir
  aux.permute(Indices(2,1,4,3));
  aux.reshape(Indices(xi,1,alphaU*xi,d2));aux.conjugate();
  Operator firstOp1(aux); // no edge

  aux=midU; // du xil dd xir
  aux.reshape(Indices(d2*xi*d2,xi));
  aux.multiplyRight(edgeUR);// du xil dd dr
  aux.reshape(Indices(d2,xi*d2*d2));
  aux.multiplyLeft(opUp2); // alphaU x xi*dd*dr
  aux.reshape(Indices(1,alphaU*xi,d2,d2)); // check?
  aux.permute(Indices(2,1,4,3));aux.conjugate();
  Operator firstOp2(aux); // includes an edgeUR

  // 2) from mpoLow for the lower part, which gets contracted with the left
  mwArray opLow1=mpoLow.getOp(0).getFullData(); // d 1 d alpha
  mwArray opLow2=mpoLow.getOp(1).getFullData(); // d alpha d 1
  int alphaL=opLow1.getDimension(3);
  opLow1.permute(Indices(1,3,4,2));
  opLow1.reshape(Indices(d2,alphaL));
  opLow2.permute(Indices(1,3,2,4));
  opLow2.reshape(Indices(d2,alphaL));

  aux=midU; // d2(u) xi(l) d2(d) xi(r)
  aux.permute(Indices(2,1,3,4)); // xi(l) d2(u) d2(d) xi(r)
  aux.reshape(Indices(xi,d2*d2*xi));
  aux.multiplyLeft(edgeULb); // d2(l) x d2(u)*d2(d)*xi(r)
  aux.reshape(Indices(d2,d2,d2,xi)); //d2(l) d2(u) d2(d) xi(r)
  aux.permute(Indices(4,2,1,3)); // xi(r) d2(u) d2(l) d2(d)
  aux.reshape(Indices(xi*d2*d2,d2));
  aux.multiplyRight(opLow1); // xi(r) d2(u) d2(l) alpha
  aux.reshape(Indices(xi,d2,d2,alphaL));
  aux.permute(Indices(4,1,2,3)); // alpha*xi(r) d2(u) d2(l)
  aux.reshape(Indices(alphaL*xi,d2,d2,1));
  Operator lastOp1(aux); // includes an edgeULb

  aux=midU;
  aux.reshape(Indices(d2*xi,d2,xi)); // d2(u)*xi(l) d2(d) xi(r)
  aux.permute(Indices(1,3,2)); // d2(u)*xi(l)*xi(r) d2(d) 
  aux.reshape(Indices(d2*xi*xi,d2));
  aux.multiplyRight(opLow2); // d2(u)xi(l)xi(r) alpha
  aux.reshape(Indices(d2,xi,xi,alphaL));
  aux.permute(Indices(3,1,4,2));
  aux.reshape(Indices(xi,d2,alphaL*xi,1));
  Operator lastOp2(aux); // no edge
  
  // For two sites I need a couple of special operators
  int cntT=cnt+2;
  int l=1; // distance
  // Open the file to write them
  ofstream* out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
    exit(1);
  }


  Contractor& contractor=Contractor::theContractor();

  int len=mpsB.getLength();

  // Initialize the MPS for left and right boundaries
  MPS mpsLop(mpsB); // operator at the bottom (left boundary)=> extra site at the beginning
  MPS mpsRop(mpsB); // operator at the top (right boundary)=> extra site at the end
  mpsRop.conjugateMPS();
  MPS::gaugeCond(mpsLop,mpsRop,'L');
  // I am not sure that this helps, anyway!
  double normFL=0.;
  double normFR=0.; // to keep normalization under control, as I expect it to go as d^(l+1)
  // First of all, get the (common) norm factor (then for each term there is also a d^(l+1) factor
  complex_t norm0;
  {// the MPO is shorter than later, and it needs both edges (R and L) because this case is special
    MPO mpo(len);
    mpo.setOp(0,new Operator(reshape(permute(edgeUL*edgeUR,Indices(2,1)),Indices(d2,1,d2,1))),true);
    mwArray id2=identityMatrix(d);id2.reshape(Indices(d2,1));
    mpo.setOp(len-1,new Operator(reshape(permute(edgeULb*edgeURb,Indices(2,1)),Indices(d2,1,d2,1))),true);
    for(int k=1;k<len-1;k++)
      mpo.setOp(k,&idOpXi,false);
    norm0=contractor.contract(mpsLop,mpo,mpsRop);
  }  
  //cout<<"norm0="<<norm0<<endl;

  // Now apply the operators at the edges
  mpsLop.insertSite(0,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);
  mpsRop.insertSite(len,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);

    // prepare the MPOs for the contraction; in this case, they are one site longer each
  MPO mpoEdges(len+2);
  mpoEdges.setOp(0,&edgeOpL,false);
  mpoEdges.setOp(len+1,&edgeOpR,false);
  for(int k=1;k<len+1;k++) mpoEdges.setOp(k,&idOpXi,false);
  MPO mpoTraces(len+3);
  mpoTraces.setOp(0,&traceIdL,false);
  mpoTraces.setOp(len+2,&traceIdR,false);
  for(int k=1;k<len+2;k++) mpoTraces.setOp(k,&idOpXi,false);
  //cout<<"Created MPOs for contraction MPOedges: "<<mpoEdges<<" and mpoTraces:"<<mpoTraces<<endl;

  { // since the MPOs with the op are only needed once, I close a scope here
    // Modify those MPSs to include the operators at the edges  
    MPO oplow(len+1); // the one with the operator below
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    oplow.setOp(0,&leftmostHook,false);
    oplow.setOp(1,&secondOpL,false);
    for(int k=2;k<len;k++)
      oplow.setOp(k,&middleU,false);
    oplow.setOp(len,&lastOp1,false); 
    MPO ophigh(len+1); // the one with the operator above
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    ophigh.setOp(len,&rightmostHook,false);
    ophigh.setOp(len-1,&rightsecondOpL,false);
    for(int k=1;k<len-1;k++)
      ophigh.setOp(k,&middleUtr,false);
    ophigh.setOp(0,&firstOp2,false); 
    //    cout<<"First step, apply to L oplow "<<oplow<<"\n\t to R ophigh "<<ophigh<<endl;
    MPS aux(mpsLop);
    contractor.optimize(oplow,aux,mpsLop,Dcut);
    mpsLop.gaugeCond('L',0);
    normFL+=log2(mpsLop.getNormFact());
    mpsLop.setNormFact(1.);
    aux=mpsRop;
    contractor.optimize(ophigh,aux,mpsRop,Dcut);
    mpsRop.gaugeCond('L',0);
    normFR+=log2(mpsRop.getNormFact());
    mpsRop.setNormFact(1.);
  }
  //cout<<"After applying the first operator: \n\tmpsLop:"<<mpsLop<<"\n\tmpsRop:"<<mpsRop<<endl;
    //      <<"\n\tauxL"<<auxL<<"\n\tmpo1"<<mpo1
    //<<"\n\tmpo2"<<mpo2<<"\n\tauxR"<<auxR<<endl;
  {
    // In this case, the first term, l=1, corresponds to inserting a
    // special MPO, so I should do this, then evolve again with
    // special (second op) and then the systematic part

    // First the one with the traces (t+2,l=1)
    MPO step1a(len+2);    
    for(int k=2;k<len;k++)
      step1a.setOp(k,&middleU,false);
    step1a.setOp(0,&traceIdL,false);
    step1a.setOp(len+1,&traceIdR,false);
    // Since in firstOp1 everything is transposed and conj, I need to change
    Operator firstOp1_b(permute(conjugate(firstOp1.getFullData()),Indices(3,2,1,4)));
    step1a.setOp(1,&firstOp1_b,false);
    step1a.setOp(len,&lastOp2,false);
    //cout<<"For t+2, l=1, will contract with step1a "<<step1a<<endl;
    complex_t val0_1=contractWithTracedEdges(mpsLop,mpsRop,step1a);
    val0_1*=pow(2,normFL+normFR-(l+2)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<label<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val0_1)<<"\t"<<imag(val0_1)<<endl;
    cntT++;
    // Then the one with edges (t+3, l=1)
    MPO step1(len+1);    
    for(int k=1;k<len;k++)
      step1.setOp(k,&middleU,false);
    step1.setOp(0,&firstOp1_b,false);
    step1.setOp(len,&lastOp2,false);

    // And i will use extra mpos for the contraction, to include the missing edges
    MPO aux1(len+1);
    for(int k=1;k<len;k++) aux1.setOp(k,&idOpXi,false);
    aux1.setOp(0,&edgeOpL,false);
    Operator idAlXiL(reshape(identityMatrix(alphaL*xi),Indices(alphaL*xi,1,alphaL*xi,1)));
    aux1.setOp(len,&idAlXiL,false);
    MPO aux2(len+1);
    for(int k=1;k<len;k++) aux2.setOp(k,&idOpXi,false);
    aux2.setOp(len,&edgeOpR,false);
    Operator idAlXiU(reshape(identityMatrix(alphaU*xi),Indices(alphaU*xi,1,alphaU*xi,1)));
    aux2.setOp(0,&idAlXiU,false);
    
    //cout<<"For t+3, l=1, will contract with aux1:"<<aux1<<" step1 "<<step1<<" aux2 "<<aux2<<endl;
    const MPO* ptrs[]={&aux2,&step1,&aux1};
    MPO _mpo(len+1);
    MPO::join(3,ptrs,_mpo);
    complex_t val1_1=contractor.contract(mpsLop,_mpo,mpsRop);
    val1_1*=pow(2,normFL+normFR-(l+2)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<label<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val1_1)<<"\t"<<imag(val1_1)<<endl;
    l++;
  }
  mpsLop.insertSite(0,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);
  mpsRop.insertSite(len+1,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);
  //  cout<<"Before applying the second operator: \n\tmpsLop:"<<mpsLop<<"\n\tmpsRop:"<<mpsRop<<endl;
  {
    // And now evolve them with the secoond op in the pair ****
    //** construct the left and right evolution ops that increase the length by one.
    MPO oplow(len+2); // the one with the operator below
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    oplow.setOp(0,&leftmostHook,false);
    oplow.setOp(1,&secondOpL,false);
    for(int k=2;k<len+1;k++)
      oplow.setOp(k,&middleU,false);
    oplow.setOp(len+1,&lastOp2,false); 
    MPO ophigh(len+2); // the one with the operator above
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    ophigh.setOp(len+1,&rightmostHook,false);
    ophigh.setOp(len,&rightsecondOpL,false);
    for(int k=1;k<len;k++)
      ophigh.setOp(k,&middleUtr,false);
    ophigh.setOp(0,&firstOp1,false); 

    //cout<<"Second step, apply to L oplow "<<oplow<<"\n\t to R ophigh "<<ophigh<<endl;
    MPS aux(mpsLop);
    contractor.optimize(oplow,aux,mpsLop,Dcut);
    mpsLop.gaugeCond('L',0);
    normFL+=log2(mpsLop.getNormFact());
    mpsLop.setNormFact(1.);
    aux=mpsRop;
    contractor.optimize(ophigh,aux,mpsRop,Dcut);
    mpsRop.gaugeCond('L',0);
    normFR+=log2(mpsRop.getNormFact());
    mpsRop.setNormFact(1.);
    
  }

  //cout<<"After applying the operator columns, length of mps:"<<mpsLop.getLength()<<endl;
  
  // Now iterate, producing contractions and evolving once more
  // Prepare MPOs for the "evolution"
  MPO oplow(len+3); // the one for the left term with the operator below(now no operator, though)
  oplow.setOp(0,&leftmostHook,false);
  oplow.setOp(1,&secondOpL,false);
  for(int k=2;k<len+2;k++)
    oplow.setOp(k,&middleU,false);
  oplow.setOp(len+2,&leftHookWithEdge,false);
  MPO ophigh(len+3); // the one with the operator above
  // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
  ophigh.setOp(len+2,&rightmostHook,false);
  ophigh.setOp(len+1,&rightsecondOpL,false);
  for(int k=1;k<len+1;k++)
    ophigh.setOp(k,&middleUtr,false);
  ophigh.setOp(0,&rightHookWithEdge,false);

  //cout<<"Created MPOs for evol oplow: "<<oplow<<" and ophigh:"<<ophigh<<endl;
  
  while(cntT<=maxM){
    //cout<<"Starting contractions within step "<<cnt<<" now cntT "<<cntT<<", length of mps "<<mpsRop.getLength()<<endl;
    // first, contract the left with the right with traced out edges (t+2, i.e. cntT, l=1)
    complex_t val0_1=contractWithTracedEdges(mpsLop,mpsRop,mpoTraces);
    val0_1*=pow(2,normFL+normFR-(l+2)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<label<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val0_1)<<"\t"<<imag(val0_1)<<endl;
    //cout<<"CC("<<whichOp<<")[t="<<cntT<<",l="<<l<<"]="<<val0_1/norm0<<endl;
    // then contract both with the operator with edges (t+3,l=1)
    cntT++;
    complex_t val1_1=contractor.contract(mpsLop,mpoEdges,mpsRop);
    val1_1*=pow(2,normFL+normFR-(l+2)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<label<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val1_1)<<"\t"<<imag(val1_1)<<endl;
    //cout<<"CC("<<whichOp<<")[t="<<cntT<<",l="<<l<<"]="<<val1_1/norm0<<endl;
    // now evolve the left (op below)
    applyOneStepLeft(mpsLop,oplow,Dcut);
    mpsLop.gaugeCond('L',0);
    normFL+=log2(mpsLop.getNormFact());
    mpsLop.setNormFact(1.);
    // contract with previous right with traced out (t+3, l=2)
    l++;
    complex_t val1_2=contractWithTracedEdges(mpsLop,mpsRop,mpoTraces);    
    val1_2*=pow(2,normFL+normFR-(l+2)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<label<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val1_2)<<"\t"<<imag(val1_2)<<endl;
    // contract with previous right with edges (t+4, l=2)
    cntT++;
    complex_t val2_2=contractor.contract(mpsLop,mpoEdges,mpsRop);
    val2_2*=pow(2,normFL+normFR-(l+2)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<label<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val2_2)<<"\t"<<imag(val2_2)<<endl;
    // evolve right (and I will be at cnt+=2, l+=2)
    applyOneStepRight(mpsRop,ophigh,Dcut);
    mpsRop.gaugeCond('L',0);
    normFR+=log2(mpsRop.getNormFact());
    mpsRop.setNormFact(1.);
    l++;
  }
  out->close();
  delete out;
}

 void computeMPOCorrelatorsSameSite(MPS& mpsB,const mwArray& midU,
				    const mwArray& edgeUL,const mwArray& edgeULb,
				    const mwArray& edgeUR,const mwArray& edgeURb,
				    int cnt,double delta,const string& outfname,
				    int Dcut,int maxM,const MPO& mpoLow,const MPO& mpoUp,int label){
   int nops=mpoLow.getLength();
   if(nops!=2||mpoUp.getLength()!=nops){
     cout<<"ERROR: For now, correlators of MPO can only be computed for 2-site case on both sides"<<endl;
     exit(1);
   }

   int d2=d*d;
   int xi=midU.getDimension(1);
  
   // To compute the correlator, I construct a MPO
   Operator middleU(midU,Indices(4,1,2,3)); // xi(r) d(u) xi(l) d(d)
   Operator edgeOpL(reshape(permute(edgeUL,Indices(2,1)),Indices(xi,1,d2,1)));
   Operator edgeOpR(reshape(permute(edgeUR,Indices(2,1)),Indices(d2,1,xi,1)));
   Operator edgeOpLb(reshape(permute(edgeULb,Indices(2,1)),Indices(xi,1,d2,1)));
   Operator edgeOpRb(reshape(permute(edgeURb,Indices(2,1)),Indices(d2,1,xi,1)));
   Operator idOpXi(reshape(identityMatrix(xi),Indices(xi,1,xi,1)));
   Operator traceIdL(reshape(identityMatrix(d),Indices(1,1,d*d,1)));
   Operator traceIdR(reshape(identityMatrix(d),Indices(d*d,1,1,1)));

   // Now the operators that depend on the MPOs (no edges included) 1)
   // for the upper part (from mpoUp), which get contracted with the
   // right side, (here no conj, as there is a single contraction!)
   mwArray opUp1=mpoUp.getOp(0).getFullData(); // d 1 d alpha
   int alphaU=opUp1.getDimension(3);
   opUp1.permute(Indices(4,3,1,2)); // alpha d(d)d(u) 1
   opUp1.reshape(Indices(alphaU,d2));
   mwArray opUp2=mpoUp.getOp(1).getFullData(); // d alpha d 1
   opUp2.permute(Indices(2,3,1,4)); // alpha d(d) d(u) 1
   opUp2.reshape(Indices(alphaU,d2));

   mwArray aux(midU);// du xil dd xir
   aux.reshape(Indices(d2,xi*d2*xi));
   aux.multiplyLeft(opUp1);
   aux.reshape(Indices(alphaU,xi,d2,xi)); // alpha xil dd xir
   aux.permute(Indices(1,4,2,3));
   aux.reshape(Indices(alphaU*xi,1,xi,d2));
   Operator firstOp1(aux); // no edge

   aux=midU; // du xil dd xir
   aux.reshape(Indices(d2,xi*d2*xi));
   aux.multiplyLeft(opUp2);
   aux.reshape(Indices(alphaU,xi,d2,xi)); // alpha xil dd xir
   aux.permute(Indices(4,1,2,3)); // xir alpha*xil, dd
   aux.reshape(Indices(xi,1,alphaU*xi,d2)); // alpha xil dd xir
   Operator firstOp2(aux); // no edge

   // 2) from mpoLow for the lower part, which gets contracted with the left
   mwArray opLow1=mpoLow.getOp(0).getFullData(); // d 1 d alpha
   mwArray opLow2=mpoLow.getOp(1).getFullData(); // d alpha d 1
   int alphaL=opLow1.getDimension(3);
   opLow1.permute(Indices(1,3,4,2));
   opLow1.reshape(Indices(d2,alphaL));
   opLow2.permute(Indices(1,3,2,4));
   opLow2.reshape(Indices(d2,alphaL));

   aux=midU;  // d2(u) xi(l) d2(d) xi(r)
   aux.reshape(Indices(d2*xi,d2,xi)); // d2(u)*xi(l) d2(d) xi(r)
   aux.permute(Indices(1,3,2)); // d2(u)*xi(l)*xi(r) d2(d) 
   aux.reshape(Indices(d2*xi*xi,d2));
   aux.multiplyRight(opLow1); // d2(u)xi(l)xi(r) alpha
   aux.reshape(Indices(d2,xi,xi,alphaL));
   aux.permute(Indices(4,3,1,2)); // alpha*xi(r) d2(u) xil
   aux.reshape(Indices(alphaL*xi,d2,xi,1));
   Operator lastOp1(aux); // no edge
   
   aux=midU;
   aux.reshape(Indices(d2*xi,d2,xi)); // d2(u)*xi(l) d2(d) xi(r)
   aux.permute(Indices(1,3,2)); // d2(u)*xi(l)*xi(r) d2(d) 
   aux.reshape(Indices(d2*xi*xi,d2));
   aux.multiplyRight(opLow2); // d2(u)xi(l)xi(r) alpha
   aux.reshape(Indices(d2,xi,xi,alphaL));
   aux.permute(Indices(3,1,4,2));
   aux.reshape(Indices(xi,d2,alphaL*xi,1));
   Operator lastOp2(aux); // no edge

   // Open the file to write them
   ofstream* out=new ofstream(outfname.data(),ios::app);
   if(!out->is_open()){
     cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
     exit(1);
   }


   Contractor& contractor=Contractor::theContractor();
   int len=mpsB.getLength();

   // Suboptimal!! I may be copying (and gauging) too many times
   MPS bra(mpsB);bra.conjugateMPS();
   //MPS::gaugeCond(mpsB,bra,'L');

   //cout<<"norm="<<norm<<endl;  
   MPO mpo1(len);
   mpo1.setOp(0,&firstOp1,false);
   mpo1.setOp(len-1,&lastOp1,false);
   MPO mpo2(len);
   mpo2.setOp(0,&firstOp2,false);
   mpo2.setOp(len-1,&lastOp2,false);
   for(int k=1;k<len-1;k++){
     mpo1.setOp(k,&middleU,false);
     mpo2.setOp(k,&middleU,false);
   }
   // and two auxiliary MPOs with just the missing edges
   MPO auxL(len);MPO auxR(len);
   auxL.setOp(0,&edgeOpL,false);
   auxR.setOp(0,&edgeOpR,false);
   auxL.setOp(len-1,&edgeOpLb,false);
   auxR.setOp(len-1,&edgeOpRb,false);
   for(int k=1;k<len-1;k++){
     auxL.setOp(k,&idOpXi,false);
     auxR.setOp(k,&idOpXi,false);
   }
   
   // cout<<"For step+2 I am going to use: \n\tmpsB:"<<mpsB<<"\n\tbra:"<<bra
   //     <<"\n\tauxL"<<auxL<<"\n\tmpo1"<<mpo1
   //     <<"\n\tmpo2"<<mpo2<<"\n\tauxR"<<auxR<<endl;
   const MPO* ptrs[]={&auxR,&mpo2,&mpo1,&auxL};
   MPO twoSteps(len);
   MPO::join(4,ptrs,twoSteps);

   complex_t norm; // The nor factor wrt the terms needs to include the edges that are missing here (i.e. auxL and auxR)
   {const MPO* ptrs2[]={&auxR,&auxL};
     MPO normMPO(len);
     MPO::join(2,ptrs2,normMPO);
     norm=contractor.contract(mpsB,normMPO,bra);
   }
  
  complex_t val2=contractor.contract(mpsB,twoSteps,bra);
  val2*=pow(d,-2)*ONE_c;
  *out<<(cnt+2)*delta<<"\t"<<0<<"\t"<<label<<"\t"<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(val2)<<"\t"<<imag(val2)<<endl;


  // I compute also the corr one step after, as I can reuse operators

  // I reuse mpo1 and 2 replacing the first operator by the mpos without midU
  mpo1.setOp(0,new Operator(reshape(opUp1,Indices(alphaU,1,1,d2))),true);
  mpo2.setOp(0,new Operator(reshape(opUp2,Indices(1,1,alphaU,d2))),true);
  // And the extra ones with identities are also different (edges below but above traces)
  auxL.setOp(0,&traceIdL,false);
  auxR.setOp(0,&traceIdR,false);
  MPO::join(4,ptrs,twoSteps);

  // cout<<"For step+1 I am going to use: \n\tmpsB:"<<mpsB<<"\n\tbra:"<<bra
  //     <<"\n\tauxL"<<auxL<<"\n\tmpo1"<<mpo1
  //     <<"\n\tmpo2"<<mpo2<<"\n\tauxR"<<auxR<<endl;

  complex_t val1=contractor.contract(mpsB,twoSteps,bra);
  val1*=pow(d,-2)*ONE_c;
  *out<<(cnt+1)*delta<<"\t"<<0<<"\t"<<label<<"\t"<<real(norm)<<"\t"<<imag(norm)<<"\t";
  *out<<real(val1)<<"\t"<<imag(val1)<<endl;

  
  out->close();delete out;

 }

void computeMPODiagonalCorrelators(const mwArray& edgeUL,const mwArray& edgeURb,const mwArray& midU,
				   double delta,const string& outfname,int maxM,const MPO& mpoLow,
				   const MPO& mpoUp,int label,int Dcut){
   int nops=mpoLow.getLength();
   // if(nops!=2||mpoUp.getLength()!=nops){
   //   cout<<"ERROR: For now, correlators of MPO can only be computed for 2-site case on both sides"<<endl;
   //   exit(1);
   // }

   int d2=d*d;
   int xi=midU.getDimension(1);
  
   MPO step(nops+1);
   step.setOp(0,new Operator(reshape(edgeUL,Indices(1,1,d2,xi))),true);
   for(int k=1;k<nops;k++)
     step.setOp(k,new Operator(midU),true); // could include a factor 2 in here
   step.setOp(nops,new Operator(reshape(permute(edgeURb,Indices(2,1)),Indices(d2,xi,1,1))),true);

   MPS mpsLow;
   MPSfromMPO(mpoLow,mpsLow,true);
   MPS mpsUp;
   MPSfromMPO(mpoUp,mpsUp,false);
   mpsUp.insertSite(0,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);//  For a better starting point, I "move" the sites down
   mpsUp.conjugateMPS(); // to compensate for later

   double normF=0.;
   
   int cnt=1;
   ofstream* out=new ofstream(outfname.data(),ios::app);
   if(!out->is_open()){
     cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
     exit(1);
   }
   Contractor& contractor=Contractor::theContractor();
   while(cnt<maxM){
     //     applyOneStepRight(mpsLow,step,Dcut); // could use this, but I want to skip the tracing out
     MPS aux(mpsLow);
     aux.insertSite(nops,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);// empty site at the end  
     mpsLow.insertSite(0,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);//  For a better starting point, I "move" the sites down
     contractor.optimize(step,aux,mpsLow,Dcut);

     mpsLow.gaugeCond('L',0);
     normF+=log2(mpsLow.getNormFact());
     mpsLow.setNormFact(1.);
     
     complex_t val=contractor.contract(mpsLow,mpsUp);
     val*=pow(2,normF-(cnt+2)*log2(d))*ONE_c;
     // now absorb the last site in the previous one and copy everything to the original mps
     // I use the trick of "tracing out" the site
     vector<int> pos(1,0);
     mpsLow.traceOutSites(pos,false); // do not normalize again

     *out<<cnt*delta<<"\t"<<cnt<<"\t"<<label<<"\t"<<1.<<"\t"<<0.<<"\t"<<real(val)<<"\t"<<imag(val)<<endl;
     cnt++;
   }
   out->close();delete out;
   
}


///
void computeDistantCorrelatorsTwoBodyOpt(const MPS& mpsB,
					 const mwArray& midU,
					 const mwArray& edgeUL,const mwArray& edgeULb,
					 const mwArray& edgeUR,const mwArray& edgeURb,
					 int cnt,double delta,const string& outfname,
					 int Dcut,int maxM,const MPO& mpoLow,const MPO& mpoUp,int label){
  int nops=mpoLow.getLength();
  if(nops!=2||mpoUp.getLength()!=nops){
    cout<<"ERROR: For now, correlators of MPO can only be computed for 2-site case on both sides"<<endl;
    exit(1);
  }

  int d2=d*d;
  int xi=midU.getDimension(1);
  
  // To compute the correlator, I construct a MPO
  Operator middleU(midU,Indices(4,1,2,3)); // xi(r) d(u) xi(l) d(d)
  Operator middleUtr(conjugate(midU),Indices(2,1,4,3)); // xi(l) d(u) xi(r) d(d)
  Operator leftmostHook(reshape(identityMatrix(d2),Indices(d2,1,1,d2))); // "hooks" (just identities)
  Operator rightmostHook(reshape(identityMatrix(d2),Indices(d2,d2,1,1))); // would conjugate, no need as it is Id
  Operator leftHookWithEdge(reshape(permute(edgeURb,Indices(2,1)),Indices(1,d2,xi,1)));
  Operator rightHookWithEdge(reshape(conjugate(edgeUR),Indices(1,1,xi,d2))); // for subsequent steps the hook carries the "closing" edge
  Operator edgeOpL(reshape(permute(edgeUL,Indices(2,1)),Indices(xi,1,d2,1)));
  Operator edgeOpR(reshape(permute(edgeURb,Indices(2,1)),Indices(d2,1,xi,1)));
  //cout<<"edgeUL*edgeUR="<<edgeUL*edgeUR<<endl;
  //cout<<" in  edgeOpL:"<<edgeOpL.getFullData()<<endl;exit(1);
  //cout<<"edgeULb"<<edgeULb<<endl;exit(1);
  Operator idOpXi(reshape(identityMatrix(xi),Indices(xi,1,xi,1)));
  Operator traceIdL(reshape(identityMatrix(d),Indices(1,1,d*d,1)));
  Operator traceIdR(reshape(identityMatrix(d),Indices(d*d,1,1,1)));
  Operator secondOpL(reshape(edgeUL*reshape(permute(midU,Indices(2,1,3,4)),
					    Indices(xi,d2*d2*xi)), // d2 x d2*d2*xi
			     Indices(d2,d2,d2,xi)),Indices(4,2,1,3)); // xi x d2 x d2 (incl edge) x d2
  Operator rightsecondOpL(conjugate(permute(reshape(reshape(midU,Indices(d2*xi*d2,d2))*edgeURb,
						    Indices(d2,xi,d2,d2)),Indices(2,1,4,3))));

  // Now the operators that depend on the MPOs
  // 1) for the upper part (from mpoUp), which get contracted witht eh right side, but in conjugate way
  mwArray opUp1=mpoUp.getOp(0).getFullData(); // d 1 d alpha
  int alphaU=opUp1.getDimension(3);
  opUp1.permute(Indices(4,3,1,2)); // alpha d(d)d(u) 1
  opUp1.reshape(Indices(alphaU,d2));
  mwArray opUp2=mpoUp.getOp(1).getFullData(); // d alpha d 1
  opUp2.permute(Indices(2,3,1,4)); // alpha d(d) d(u) 1
  opUp2.reshape(Indices(alphaU,d2));

  mwArray aux(midU);// du xil dd xir
  aux.reshape(Indices(d2,xi*d2*xi));
  aux.multiplyLeft(opUp1);
  aux.reshape(Indices(alphaU,xi,d2,xi)); // alpha xil dd xir
  aux.permute(Indices(2,1,4,3));
  aux.reshape(Indices(xi,1,alphaU*xi,d2));aux.conjugate();
  Operator firstOp1(aux); // no edge

  aux=midU; // du xil dd xir
  aux.reshape(Indices(d2*xi*d2,xi));
  aux.multiplyRight(edgeUR);// du xil dd dr
  aux.reshape(Indices(d2,xi*d2*d2));
  aux.multiplyLeft(opUp2); // alphaU x xi*dd*dr
  aux.reshape(Indices(1,alphaU*xi,d2,d2)); // check?
  aux.permute(Indices(2,1,4,3));aux.conjugate();
  Operator firstOp2(aux); // includes an edgeUR

  // 2) from mpoLow for the lower part, which gets contracted with the left
  mwArray opLow1=mpoLow.getOp(0).getFullData(); // d 1 d alpha
  mwArray opLow2=mpoLow.getOp(1).getFullData(); // d alpha d 1
  int alphaL=opLow1.getDimension(3);
  opLow1.permute(Indices(1,3,4,2));
  opLow1.reshape(Indices(d2,alphaL));
  opLow2.permute(Indices(1,3,2,4));
  opLow2.reshape(Indices(d2,alphaL));

  aux=midU; // d2(u) xi(l) d2(d) xi(r)
  aux.permute(Indices(2,1,3,4)); // xi(l) d2(u) d2(d) xi(r)
  aux.reshape(Indices(xi,d2*d2*xi));
  aux.multiplyLeft(edgeULb); // d2(l) x d2(u)*d2(d)*xi(r)
  aux.reshape(Indices(d2,d2,d2,xi)); //d2(l) d2(u) d2(d) xi(r)
  aux.permute(Indices(4,2,1,3)); // xi(r) d2(u) d2(l) d2(d)
  aux.reshape(Indices(xi*d2*d2,d2));
  aux.multiplyRight(opLow1); // xi(r) d2(u) d2(l) alpha
  aux.reshape(Indices(xi,d2,d2,alphaL));
  aux.permute(Indices(4,1,2,3)); // alpha*xi(r) d2(u) d2(l)
  aux.reshape(Indices(alphaL*xi,d2,d2,1));
  Operator lastOp1(aux); // includes an edgeULb

  aux=midU;
  aux.reshape(Indices(d2*xi,d2,xi)); // d2(u)*xi(l) d2(d) xi(r)
  aux.permute(Indices(1,3,2)); // d2(u)*xi(l)*xi(r) d2(d) 
  aux.reshape(Indices(d2*xi*xi,d2));
  aux.multiplyRight(opLow2); // d2(u)xi(l)xi(r) alpha
  aux.reshape(Indices(d2,xi,xi,alphaL));
  aux.permute(Indices(3,1,4,2));
  aux.reshape(Indices(xi,d2,alphaL*xi,1));
  Operator lastOp2(aux); // no edge
  
  // For two sites I need a couple of special operators
  int cntT=cnt+2;
  int l=1; // distance
  // Open the file to write them
  ofstream* out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"ERROR! Cannot open file "<<outfname<<" for output!"<<endl;
    exit(1);
  }


  Contractor& contractor=Contractor::theContractor();

  int len=mpsB.getLength();

  // Initialize the MPS for left and right boundaries
  MPS mpsLop(mpsB); // operator at the bottom (left boundary)=> extra site at the beginning
  MPS mpsRop(mpsB); // operator at the top (right boundary)=> extra site at the end
  mpsRop.conjugateMPS();
  MPS::gaugeCond(mpsLop,mpsRop,'L');
  // I am not sure that this helps, anyway!
  double normFL=0.;
  double normFR=0.; // to keep normalization under control, as I expect it to go as d^(l+1)
  // First of all, get the (common) norm factor (then for each term there is also a d^(l+1) factor
  complex_t norm0;
  {// the MPO is shorter than later, and it needs both edges (R and L) because this case is special
    MPO mpo(len);
    mpo.setOp(0,new Operator(reshape(permute(edgeUL*edgeUR,Indices(2,1)),Indices(d2,1,d2,1))),true);
    mwArray id2=identityMatrix(d);id2.reshape(Indices(d2,1));
    mpo.setOp(len-1,new Operator(reshape(permute(edgeULb*edgeURb,Indices(2,1)),Indices(d2,1,d2,1))),true);
    for(int k=1;k<len-1;k++)
      mpo.setOp(k,&idOpXi,false);
    norm0=contractor.contract(mpsLop,mpo,mpsRop);
  }  
  //cout<<"norm0="<<norm0<<endl;

  // Now apply the operators at the edges
  mpsLop.insertSite(0,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);
  mpsRop.insertSite(len,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);

    // prepare the MPOs for the contraction; in this case, they are one site longer each
  MPO mpoEdges(len+2);
  mpoEdges.setOp(0,&edgeOpL,false);
  mpoEdges.setOp(len+1,&edgeOpR,false);
  for(int k=1;k<len+1;k++) mpoEdges.setOp(k,&idOpXi,false);
  MPO mpoTraces(len+3);
  mpoTraces.setOp(0,&traceIdL,false);
  mpoTraces.setOp(len+2,&traceIdR,false);
  for(int k=1;k<len+2;k++) mpoTraces.setOp(k,&idOpXi,false);
  //cout<<"Created MPOs for contraction MPOedges: "<<mpoEdges<<" and mpoTraces:"<<mpoTraces<<endl;

  { // since the MPOs with the op are only needed once, I close a scope here
    // Modify those MPSs to include the operators at the edges  
    MPO oplow(len+1); // the one with the operator below
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    oplow.setOp(0,&leftmostHook,false);
    oplow.setOp(1,&secondOpL,false);
    for(int k=2;k<len;k++)
      oplow.setOp(k,&middleU,false);
    oplow.setOp(len,&lastOp1,false); 
    MPO ophigh(len+1); // the one with the operator above
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    ophigh.setOp(len,&rightmostHook,false);
    ophigh.setOp(len-1,&rightsecondOpL,false);
    for(int k=1;k<len-1;k++)
      ophigh.setOp(k,&middleUtr,false);
    ophigh.setOp(0,&firstOp2,false); 
    //    cout<<"First step, apply to L oplow "<<oplow<<"\n\t to R ophigh "<<ophigh<<endl;
    MPS aux(mpsLop);
    contractor.optimize(oplow,aux,mpsLop,Dcut);
    mpsLop.gaugeCond('L',0);
    normFL+=log2(mpsLop.getNormFact());
    mpsLop.setNormFact(1.);
    aux=mpsRop;
    contractor.optimize(ophigh,aux,mpsRop,Dcut);
    mpsRop.gaugeCond('L',0);
    normFR+=log2(mpsRop.getNormFact());
    mpsRop.setNormFact(1.);
  }
  //cout<<"After applying the first operator: \n\tmpsLop:"<<mpsLop<<"\n\tmpsRop:"<<mpsRop<<endl;
    //      <<"\n\tauxL"<<auxL<<"\n\tmpo1"<<mpo1
    //<<"\n\tmpo2"<<mpo2<<"\n\tauxR"<<auxR<<endl;
  {
    // In this case, the first term, l=1, corresponds to inserting a
    // special MPO, so I should do this, then evolve again with
    // special (second op) and then the systematic part

    // First the one with the traces (t+2,l=1)
    MPO step1a(len+2);    
    for(int k=2;k<len;k++)
      step1a.setOp(k,&middleU,false);
    step1a.setOp(0,&traceIdL,false);
    step1a.setOp(len+1,&traceIdR,false);
    // Since in firstOp1 everything is transposed and conj, I need to change
    Operator firstOp1_b(permute(conjugate(firstOp1.getFullData()),Indices(3,2,1,4)));
    step1a.setOp(1,&firstOp1_b,false);
    step1a.setOp(len,&lastOp2,false);
    //cout<<"For t+2, l=1, will contract with step1a "<<step1a<<endl;
    complex_t val0_1=contractWithTracedEdges(mpsLop,mpsRop,step1a);
    val0_1*=pow(2,normFL+normFR-(l+2)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<label<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val0_1)<<"\t"<<imag(val0_1)<<endl;
    cntT++;
    // Then the one with edges (t+3, l=1)
    MPO step1(len+1);    
    for(int k=1;k<len;k++)
      step1.setOp(k,&middleU,false);
    step1.setOp(0,&firstOp1_b,false);
    step1.setOp(len,&lastOp2,false);

    // And i will use extra mpos for the contraction, to include the missing edges
    MPO aux1(len+1);
    for(int k=1;k<len;k++) aux1.setOp(k,&idOpXi,false);
    aux1.setOp(0,&edgeOpL,false);
    Operator idAlXiL(reshape(identityMatrix(alphaL*xi),Indices(alphaL*xi,1,alphaL*xi,1)));
    aux1.setOp(len,&idAlXiL,false);
    MPO aux2(len+1);
    for(int k=1;k<len;k++) aux2.setOp(k,&idOpXi,false);
    aux2.setOp(len,&edgeOpR,false);
    Operator idAlXiU(reshape(identityMatrix(alphaU*xi),Indices(alphaU*xi,1,alphaU*xi,1)));
    aux2.setOp(0,&idAlXiU,false);
    
    //cout<<"For t+3, l=1, will contract with aux1:"<<aux1<<" step1 "<<step1<<" aux2 "<<aux2<<endl;
    const MPO* ptrs[]={&aux2,&step1,&aux1};
    MPO _mpo(len+1);
    MPO::join(3,ptrs,_mpo);
    complex_t val1_1=contractor.contract(mpsLop,_mpo,mpsRop);
    val1_1*=pow(2,normFL+normFR-(l+2)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<label<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val1_1)<<"\t"<<imag(val1_1)<<endl;
    l++;
  }
  mpsLop.insertSite(0,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);
  mpsRop.insertSite(len+1,reshape(mwArray(ONE_c),Indices(1,1,1)),true,false);
  //  cout<<"Before applying the second operator: \n\tmpsLop:"<<mpsLop<<"\n\tmpsRop:"<<mpsRop<<endl;
  {
    // And now evolve them with the secoond op in the pair ****
    //** construct the left and right evolution ops that increase the length by one.
    MPO oplow(len+2); // the one with the operator below
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    oplow.setOp(0,&leftmostHook,false);
    oplow.setOp(1,&secondOpL,false);
    for(int k=2;k<len+1;k++)
      oplow.setOp(k,&middleU,false);
    oplow.setOp(len+1,&lastOp2,false); 
    MPO ophigh(len+2); // the one with the operator above
    // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
    ophigh.setOp(len+1,&rightmostHook,false);
    ophigh.setOp(len,&rightsecondOpL,false);
    for(int k=1;k<len;k++)
      ophigh.setOp(k,&middleUtr,false);
    ophigh.setOp(0,&firstOp1,false); 

    //cout<<"Second step, apply to L oplow "<<oplow<<"\n\t to R ophigh "<<ophigh<<endl;
    MPS aux(mpsLop);
    contractor.optimize(oplow,aux,mpsLop,Dcut);
    mpsLop.gaugeCond('L',0);
    normFL+=log2(mpsLop.getNormFact());
    mpsLop.setNormFact(1.);
    aux=mpsRop;
    contractor.optimize(ophigh,aux,mpsRop,Dcut);
    mpsRop.gaugeCond('L',0);
    normFR+=log2(mpsRop.getNormFact());
    mpsRop.setNormFact(1.);
    
  }

  //cout<<"After applying the operator columns, length of mps:"<<mpsLop.getLength()<<endl;
  
  // Now iterate, producing contractions and evolving once more
  // Prepare MPOs for the "evolution"
  MPO oplow(len+3); // the one for the left term with the operator below(now no operator, though)
  oplow.setOp(0,&leftmostHook,false);
  oplow.setOp(1,&secondOpL,false);
  for(int k=2;k<len+2;k++)
    oplow.setOp(k,&middleU,false);
  oplow.setOp(len+2,&leftHookWithEdge,false);
  MPO ophigh(len+3); // the one with the operator above
  // for the first column, an op is inserted below, but the part above is exactly as in the evolution MPO
  ophigh.setOp(len+2,&rightmostHook,false);
  ophigh.setOp(len+1,&rightsecondOpL,false);
  for(int k=1;k<len+1;k++)
    ophigh.setOp(k,&middleUtr,false);
  ophigh.setOp(0,&rightHookWithEdge,false);

  //cout<<"Created MPOs for evol oplow: "<<oplow<<" and ophigh:"<<ophigh<<endl;
  
  while(cntT<=maxM){
    //cout<<"Starting contractions within step "<<cnt<<" now cntT "<<cntT<<", length of mps "<<mpsRop.getLength()<<endl;
    // first, contract the left with the right with traced out edges (t+2, i.e. cntT, l=1)
    complex_t val0_1=contractWithTracedEdges(mpsLop,mpsRop,mpoTraces);
    val0_1*=pow(2,normFL+normFR-(l+2)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<label<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val0_1)<<"\t"<<imag(val0_1)<<endl;
    //cout<<"CC("<<whichOp<<")[t="<<cntT<<",l="<<l<<"]="<<val0_1/norm0<<endl;
    // then contract both with the operator with edges (t+3,l=1)
    cntT++;
    complex_t val1_1=contractor.contract(mpsLop,mpoEdges,mpsRop);
    val1_1*=pow(2,normFL+normFR-(l+2)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<label<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val1_1)<<"\t"<<imag(val1_1)<<endl;
    //cout<<"CC("<<whichOp<<")[t="<<cntT<<",l="<<l<<"]="<<val1_1/norm0<<endl;
    // now evolve the left (op below)
    applyOneStepLeft(mpsLop,oplow,Dcut);
    mpsLop.gaugeCond('L',0);
    normFL+=log2(mpsLop.getNormFact());
    mpsLop.setNormFact(1.);
    // contract with previous right with traced out (t+3, l=2)
    l++;
    complex_t val1_2=contractWithTracedEdges(mpsLop,mpsRop,mpoTraces);    
    val1_2*=pow(2,normFL+normFR-(l+2)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<label<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val1_2)<<"\t"<<imag(val1_2)<<endl;
    // contract with previous right with edges (t+4, l=2)
    cntT++;
    complex_t val2_2=contractor.contract(mpsLop,mpoEdges,mpsRop);
    val2_2*=pow(2,normFL+normFR-(l+2)*log2(d))*ONE_c;
    *out<<cntT*delta<<"\t"<<l<<"\t"<<label<<"\t"<<real(norm0)<<"\t"<<imag(norm0)<<"\t"<<real(val2_2)<<"\t"<<imag(val2_2)<<endl;
    // evolve right (and I will be at cnt+=2, l+=2)
    //applyOneStepRight(mpsRop,ophigh,Dcut);
    //mpsRop.gaugeCond('L',0);
    // instead, just copy the left one in the opposite order and conjugate
    int lenL=mpsLop.getLength();
    mpsRop=MPS(lenL,1,1);
    for(int pos=0;pos<mpsLop.getLength();pos++){
      //mpsRop.setRotatedA(pos,mpsLop.getA(lenL-1-pos).getA(),Indices(1,3,2),true,false); // flip and conjugate, do not check Ds
      mpsRop.setRotatedA(pos,mpsLop.getA(lenL-1-pos).getA(),Indices(1,3,2),false,false); // flip and not conjugate (done later), do not check Ds
    }
    //mpsRop.gaugeCond('L',0); // not needed
    normFR=normFL;//+=log2(mpsRop.getNormFact());
    mpsRop.setNormFact(1.);
    l++;
  }
  out->close();
  delete out;
}
