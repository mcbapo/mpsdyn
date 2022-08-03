#include <math.h>
#include <iomanip>
#include <sstream>

#include "misc.h"
#include "uMPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "SpinMPO.h"
#include "Properties.h"
#include "time.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

//#define SYMM 1 // If SYMM=1, only left part is computed explicitly, and right is just the reflection

#include "IsingHamiltonian.h"

using namespace std;
using namespace shrt;

/** Instead of simulating real time evolution, this one does the same
    for the thermal state of the infinite chain, with uMPS, as newly implemented.
 */


int d=2;
double tolSVD=1E-8;      


/** Although uMPS offers a method to compute the purity of the RDM for
    fixed length, if I want several sizes, it is more efficient to
    reuse contractions, so I repeat the code here */
void getTwoRenyiEntropies(const vector<int>& nReps,uMPS& state,vector<complex_t>& results);


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
  double scale_=props.getDoubleProperty("scale");
  if(scale_<0) scale_=1.;
  else{
    J_=J_/scale_;
    g_=g_/scale_;
    h_=h_/scale_;
  }
  int nB=props.getIntProperty("nB"); // length of the block
  if(nB<0) nB=1; // default is 1
  //int nB=1; // this MPS is TI
  double beta=props.getDoubleProperty("beta");
  double delta=props.getDoubleProperty("delta");
  int M=beta/(2*delta);
  double remainder=beta/2-M*delta;
  int Dcut=props.getIntProperty("D"); // the truncation 
  int stepBlk=props.getIntProperty("stepBlock");
  if(stepBlk<0)stepBlk=1;
  string outfname=props.getProperty("outputfile");
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
  //double eps=props.getDoubleProperty("eps");

  int maxL=props.getIntProperty("maxL"); // max size of the block
  if(maxL<0) maxL=10*nB;
  vector<int> nReps;
  for(int p=0;p<maxL/nB;p++) nReps.push_back(p+1);

  // // The initial state here is the identity (beta=0)
  // mwArray W=(1./sqrt(2.))*identityMatrix(d); // the state
  // W.reshape(Indices(d*d,1,1));
  // mwArray lambdaR(ONE_c); 

  //constructInitialState(initSt,W);

  // Get the basic piece of the evolution
  mwArray midU;
  {
    IsingHamiltonian hamH(3,d,J_,g_,h_);
    //hamH.registerTimeDependence();
    MPO _Uevol_1(3),Uevol(3); // held just temporarily
    hamH.getUMPO(_Uevol_1,-delta*.5*ONE_c,0); // the single layer normal evolution
    doubleMPO(_Uevol_1,Uevol,true); // the double layer
    midU=Uevol.getOp(1).getFullData(); // d2 x xi2_l x d2 x xi2_r 
    //ofstream offs("midU.m");
    //putForMatlab(offs,_Uevol_1.getOp(1).getFullData(),"midUop");
    //putForMatlab(offs,midU,"midUdouble");
    //offs.close();
  }
  vector<const mwArray*> Uops(nB,&midU);

  // and the ops to measure
  int d2=d*d;

  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1*I_c,ZERO_c};
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c}; //Z
  mwArray sig0_=identityMatrix(d);
  mwArray sigX_(Indices(d,d),dataX);
  mwArray sigY_(Indices(d,d),dataY);
  mwArray sigZ_(Indices(d,d),dataZ);
  mwArray sigX,sigY,sigZ,sig0; // the double ones
  constructOperatorProduct(sig0,sig0_,sig0_);
  constructOperatorProduct(sigX,sigX_,sig0_);
  constructOperatorProduct(sigY,sigY_,sig0_);
  constructOperatorProduct(sigZ,sigZ_,sig0_);
  sig0.reshape(Indices(1,d2*d2));
  sigX.reshape(Indices(1,d2*d2));
  sigY.reshape(Indices(1,d2*d2));
  sigZ.reshape(Indices(1,d2*d2));

  int cnt=0;
  double iTime=0.; // tmp beta

  ofstream* out;
  //    out=new ofstream(outfname.data(),ios::app);  
  if(!app||!file_exists(outfname.data())){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"ERROR: Cannot open file "<<outfname<<" for output "<<endl;
      exit(1);
    }
    *out<<"% t\t real(tr)\t imag(tr)\t real(sigX)\t imag(sigX)\t real(sigY)\t imag(sigY)\t"
	<<" real(sigZ)\t imag(sigZ)\t real(ZZ)\t imag(ZZ)\t <h12>"<<endl;
    out->close();delete out;
  }
  //*out<<setprecision(15);

  // Initial state is a product, the identity
  int D0=1; // starting with product state
  uMPS state(nB,d2,D0);
  state.setProductState(p_maxent); // identity state
  
  int leftSteps=M;
  while(leftSteps>0){
    // evaluate EV
    mwArray rdm(Indices(d2,d2));rdm.fillWithZero();
    for(int k=0;k<nB;k++){
      mwArray auxRDM;
      state.getSingleSiteRDM(k,auxRDM);      
      rdm=rdm+auxRDM;
      // //cout<<"Now RDM="<<rdm<<endl;
      // if(k==0){
      // 	ofstream* outS=new ofstream(outSfname.data(),ios::app);
      // 	*outS<<setprecision(15);
      // 	double Sb=computeEntropy(rdm);
      // 	*outS<<t<<"\t"<<Sb<<"\t"<<state.getEntropy(k)<<endl;
      // 	outS->close();delete outS;
    }
    rdm.multiplyLeft((1./nB)*ONE_c);
    //cout<<"Now RDM="<<rdm<<endl;
    rdm.reshape(Indices(d2*d2,1));
    mwArray res=reshape(sigX,Indices(1,d2*d2))*rdm;
    complex_t resultX=res.getElement(Indices(0,0));
    res=reshape(sigY,Indices(1,d2*d2))*rdm;
    complex_t resultY=res.getElement(Indices(0,0));
    res=reshape(sigZ,Indices(1,d2*d2))*rdm;
    complex_t resultZ=res.getElement(Indices(0,0));
    res=reshape(sig0,Indices(1,d2*d2))*rdm;
    complex_t trV=res.getElement(Indices(0,0));
    
    // For the ZZ value, with only one tensor I just need one calculation
    complex_t resultZZ(ZERO_c),resultXZ(ZERO_c);
    vector<int> pos(2,0);pos[1]=1;
    vector<mwArray*> ops;ops.push_back(&sigZ);ops.push_back(&sigZ);
    vector<mwArray*> opsXZ;opsXZ.push_back(&sigX);opsXZ.push_back(&sigZ);
    if(nB==1){
      resultZZ=state.computeExpectationValue(2,pos,ops);
      resultXZ=state.computeExpectationValue(2,pos,opsXZ);
    }
    else{
      // should loop to get an average!
      for(int k=0;k<nB;k++){
	  pos[0]=k;
	  if(k<nB-1) pos[1]=k+1;
	  else pos[1]=0;
	  resultZZ+=state.computeExpectationValue(2,pos,ops);
	  resultXZ+=state.computeExpectationValue(2,pos,opsXZ);
      }
      resultZZ=(1./nB)*resultZZ;
      resultXZ=(1./nB)*resultXZ;
    }
    // Get also purity (2-Renyi entropy)

    complex_t eta=ZERO_c;
    state.getPurity(eta);
    
    //and also the oones for blocks
    vector<complex_t> results;
    getTwoRenyiEntropies(nReps,state,results);
    
    //cout<<"resultZZ="<<resultZZ<<endl;
    // and write to file and screen
    cout<<"beta="<<iTime<<", <O>="<<(resultX)/trV<<", trace="
	<<(trV)
	<<" E="<<(J_*(resultZZ)+g_*(resultX)+h_*(resultZ))/trV
	<<" S2="<<-log(real(eta))
      //<<" S(block)="<<Sblock
	<<", D="<<state.getLeftVirtualDimension(0)<<endl;
    out=new ofstream(outfname.data(),ios::app);  *out<<setprecision(15);
    *out<<iTime<<"\t"<<(J_*(resultZZ)+g_*(resultX)+h_*(resultZ))/trV
	<<"\t"<<real(trV)<<"\t"<<imag(trV)
	<<"\t"<<real(resultX)<<"\t"<<imag(resultX)
	<<"\t"<<real(resultY)<<"\t"<<imag(resultY)
	<<"\t"<<real(resultZ)<<"\t"<<imag(resultZ)
	<<"\t"<<real(resultZZ)<<"\t"<<imag(resultZZ)
	<<"\t"<<real(resultXZ)<<"\t"<<imag(resultXZ)
	<<"\t"<<state.getLeftVirtualDimension(0)
	<<"\t"<<real(eta)<<"\t"<<imag(eta);
    for(int nr=0;nr<results.size();nr++){
      *out<<"\t"<<real(results[nr])<<"\t"<<imag(results[nr]);
    }
    *out<<endl;
    out->close();delete out;

    // apply as many steps as required
    int nrToApply=min(stepBlk,leftSteps);
    for(int k=0;k<nrToApply;k++)
      state.applyMPO(Uops,tolSVD,Dcut);

    cnt+=nrToApply;
    iTime+=2*delta*nrToApply;
    leftSteps-=nrToApply;
  }

}

#include "MPOMultiplier.h"

void getTwoRenyiEntropies(const vector<int>& nReps,uMPS& state,vector<complex_t>& results){
  results.clear();
  int nA=state.getN();
  // Construct the MPO(s) for the transfer operator in rho^2
  // Each tensor gives a MPOMultiplier
  vector<Multiplier*> E;
  for(int k=0;k<nA;k++){
    mwArray W_=state.getA(k); // dkxDlxDr
    int dk=W_.getDimension(0);
    int dk0=floor(sqrt(dk));
    if(dk0*dk0!=dk){
      cout<<"ERROR: The "<<k<<"-th tensor in this uMPO does not correspond to a purification"
	  <<", with d="<<dk<<endl;
      exit(1);
    }
    int Dl=W_.getDimension(1);
    int Dr=W_.getDimension(2);
    MPOMultiplier* Ek=new MPOMultiplier(4);
    Operator* A0=new Operator(reshape(W_,Indices(dk,Dl,Dr,1)),Indices(2,4,3,1),0);
    Ek->setOp(0,A0,true); // all mine
    Operator* A3=new Operator(reshape(W_,Indices(dk,Dl,Dr,1)),Indices(2,1,3,4),true);
    Ek->setOp(3,A3,true); // all mine
    // the middle ones are doubled with identities
    mwArray idPhys=identityMatrix(dk0);idPhys.reshape(Indices(1,dk0,1,dk0));
    W_.reshape(Indices(dk0,dk0,Dl,Dr));
    Ek->setOp(1,new DoubleOperator(permute(conjugate(W_),Indices(3,1,4,2)),idPhys),true);
    Ek->setOp(2,new DoubleOperator(permute(W_,Indices(3,2,4,1)),idPhys),true);
    E.push_back(Ek);
  }
  // Now I need the boundary vectors
  mwArray V(state.getLambda(0)); // I'll take std block
  V=V*V; // Dr x Dr
  int Dr=V.getDimension(0);
  V.reshape(Indices(Dr*Dr,1));
  V.multiplyRight(permute(V,Indices(2,1)));
  V.reshape(Indices(Dr,Dr,Dr,Dr));
  V.permute(Indices(3,2,1,4));
  V.reshape(Indices(Dr*Dr*Dr*Dr,1));

  mwArray Vl=identityMatrix(Dr*Dr);
  Vl.reshape(Indices(Dr,Dr,Dr,Dr));
  Vl.permute(Indices(1,2,4,3));
  Vl.reshape(Indices(1,Dr*Dr*Dr*Dr));

  // And now iterate block contractions until required nReps, which are saved
  // I assume they come sorted!
  int nrLs=nReps.size();int trgt=0;
  for(int kR=0;kR<nReps[nrLs-1];kR++){
    for(int k=nA-1;k>=0;k--){
      mwArray aux(V);
      E[k]->product(aux,V);
    }
    if(kR+1==nReps[trgt]){
      mwArray aux=Vl*V;
      results.push_back(aux.getElement(0));
      trgt++;
    }
  }  
}
