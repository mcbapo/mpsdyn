#include <math.h>
#include <iomanip>
#include <sstream>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "SpinMPO.h"
#include "Properties.h"
#include "time.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

//#define SYMM 1 // If SYMM=1, only left part is computed explicitly, and right is just the reflection

#include "IsingHamiltonian.h"

using namespace shrt;

// Initial product state
void constructInitialState(int initSt,mwArray& W);

void computeSingleSiteRDM(const mwArray& W,const mwArray& lambda,mwArray &rdm);

complex_t computeZZEV(const mwArray& W,const mwArray& lambda_);

void applyTruncateMPO(mwArray& W,mwArray& lambdaR,const mwArray& midU,int Dcut,double tolSVD);

double computeBlockEntanglement(const mwArray& W,const mwArray& lambda,int nB,
				vector<double>& schm,vector<double>& schm2,double& negVal);
void computeBlockDecompositions(const mwArray& W,const mwArray& lambda,int nB,
				vector<complex_t>& valsPh,vector<double>& schm);
void computeBlockPartTrace(const mwArray& W,const mwArray& lambda,int nB,vector<complex_t>& eigvals);

void truncDoubleBlockLeft(mwArray& W,mwArray& lambda,int nB,int Dcut,mwArray& S,double eps);
void truncDoubleBlock(mwArray& W,mwArray& lambda,int nB,int Dcut,mwArray& S,double eps);
void truncSingleBlock(mwArray& W,mwArray& lambda,int nB,int Dcut,mwArray& S,double eps); //,bool hardCut=false);
void specDoubleBlock(mwArray& W,mwArray& lambda,int nB,int Dcut,mwArray& S,mwArray& U,mwArray& V);

void applyMPO(mwArray& W,mwArray& lambdaR,const mwArray& midU,double tolSVD);
void applyCanonical(mwArray& W,mwArray& lambdaR,double tolSVD); // No truncation, just canonical form

//void applyMPO(mwArray& W,mwArray& lambda,const mwArray& midU);

//void truncateMPS(mwArray& W,mwArray& lambda,int Dcut,double tolSVD=0.);

int d=2;
double tolSVD=1E-8;      

/** Evolve with normal iTEBD (TI) and explore the entanglement
    spectrum of the effective block tensor for length m.  New version:
    look for CDL in the RDM. */


int maxIter=100;
double precCut=1E-8; // where to cut the eigenvalues


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
  int nB=props.getIntProperty("nB"); // length of the block
  int initSt=props.getIntProperty("initSt");
  
  double delta=props.getDoubleProperty("delta");
  int M=props.getIntProperty("M");
  int Dcut=props.getIntProperty("Dcut"); // the truncation 
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
  string outfnameSchm=props.getProperty("outputfileSchm");
  //  string outfnameSchm2=props.getProperty("outputfileSchmVert");
  string outfnameEig=props.getProperty("outputfileEig");
  double eps=props.getDoubleProperty("eps");
  int maxIter_=props.getIntProperty("maxIter");
  if(maxIter_>0) maxIter=maxIter_;
  double precCut_=props.getDoubleProperty("precCut");
  if(precCut_>0) precCut=precCut_;
  
  mwArray W; // the state
  mwArray lambdaR(ONE_c); 
  constructInitialState(initSt,W);

  // Get the basic piece of the evolution
  mwArray midU;
  {
    IsingHamiltonian hamH(3,d,J_,g_,h_);
    //hamH.registerTimeDependence();
    MPO _Uevol_1(3); // held just temporarily
    hamH.getUMPO(_Uevol_1,-delta*I_c,0); // the single layer normal evolution
    midU=_Uevol_1.getOp(1).getFullData();
  }

  // and the opsto measure
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1*I_c,ZERO_c};
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c}; //Z
  mwArray sig0=identityMatrix(d);sig0.reshape(Indices(1,d*d));
  mwArray sigX=mwArray(Indices(d,d),dataX);sigX.reshape(Indices(1,d*d));
  mwArray sigY(Indices(d,d),dataY);sigY.reshape(Indices(1,d*d));
  mwArray sigZ(Indices(d,d),dataZ);sigZ.reshape(Indices(1,d*d));

  int cnt=0;
  double time=0.;



  ofstream* out;
  ofstream* outS;
  ofstream* outE;
  if(app){
    // out=new ofstream(outfname.data(),ios::app);  *out<<setprecision(15);
    // outS=new ofstream(outfnameSchm.data(),ios::app);  *outS<<setprecision(15);
    // outE=new ofstream(outfnameEig.data(),ios::app);  *outE<<setprecision(15);
  }  
  else{
    out=new ofstream(outfname.data());
    *out<<"% t\t real(tr)\t imag(tr)\t real(sigX)\t imag(sigX)\t real(sigY)\t imag(sigY)\t"
	<<" real(sigZ)\t imag(sigZ)\t real(ZZ)\t imag(ZZ)\t Dcut"<<endl;
    out->close();delete out;
    outS=new ofstream(outfnameSchm.data());
    *outS<<"% t\t Schmidt vals " <<endl;
    outS->close();delete outS;
    //outS=new ofstream(outfnameSchm2.data());
    //*outS<<"% t\t Schmidt vals vertical" <<endl;
    //outS->close();delete outS;
    outE=new ofstream(outfnameEig.data());
    *outE<<"% t\t Eigenvalues vals " <<endl;
    outE->close();delete outE;
  }



  while(cnt<=M){
    // evaluate EV
    mwArray rdm;
    computeSingleSiteRDM(W,lambdaR,rdm);
    mwArray res=sigX*rdm;
    complex_t resultX=res.getElement(Indices(0,0));
    res=sigY*rdm;
    complex_t resultY=res.getElement(Indices(0,0));
    res=sigZ*rdm;
    complex_t resultZ=res.getElement(Indices(0,0));
    res=sig0*rdm;
    complex_t trV=res.getElement(Indices(0,0));
    complex_t resultZZ=computeZZEV(W,lambdaR);
    //complex_t resultZZ=ZERO_c;
    

    // Evaluate Schmidt decomp of block transfer op
    vector<double> schm,schm2;vector<complex_t> eigvals;
    double negVal;
    //    double Sblock=0;
    //Sblock=computeBlockEntanglement(W,lambdaR,nB,schm,schm2,negVal);


    //    computeBlockPartTrace(W,lambdaR,nB,eigvals);
    mwArray S(ZERO_c);
    //specDoubleBlock(W,lambdaR,nB,Dcut,S);
    if(W.getDimension(1)>Dcut)
      //      truncDoubleBlock(W,lambdaR,nB,Dcut,S,eps);    
      truncSingleBlock(W,lambdaR,nB,Dcut,S,eps);    
    
    // int rnk=0;
    // for(int kk=0;kk<schm.size();kk++) if(schm[kk]>1E-6) rnk++;
    // cout<<"Got Schmidt vals "<<schm.size()<<" of which "<<rnk<<" are larger than 1E-6"<<endl;
				     
    // and write to file and screen
    cout<<"t="<<time<<", <O>="<<(resultX/trV)<<", trace="
	<<(trV)
	<<" E="<<(J_*resultZZ+g_*resultX+h_*resultZ)/trV
      //<<" S(block)="<<Sblock
	<<", D="<<W.getDimension(1)<<endl;
    out=new ofstream(outfname.data(),ios::app);  *out<<setprecision(15);
    *out<<time<<"\t"<<real(trV)<<"\t"<<imag(trV)
	<<"\t"<<real(resultX)<<"\t"<<imag(resultX)
	<<"\t"<<real(resultY)<<"\t"<<imag(resultY)
	<<"\t"<<real(resultZ)<<"\t"<<imag(resultZ)
	<<"\t"<<real(resultZZ)<<"\t"<<imag(resultZZ)
      //	<<"\t"<<Sblock<<"\t"<<negVal
	<<"\t"<<W.getDimension(1)<<"\t"<<negVal
	<<endl;
    out->close();delete out;
    outS=new ofstream(outfnameSchm.data(),ios::app);  *outS<<setprecision(10);
    *outS<<time<<"\t";
    for(int k=0;k<Dcut;k++)
      if(k<schm.size()) *outS<<schm[k]<<"\t";else *outS<<0<<"\t";

    //    for(int k=0;k<schm.size();k++) *outS<<schm[k]<<"\t";	  
    *outS<<endl;
    outS->close();delete outS;

    //   outS=new ofstream(outfnameSchm2.data(),ios::app);  *outS<<setprecision(10);
    //     *outS<<time<<"\t";
    //     for(int k=0;k<Dcut;k++)
    //       if(k<schm2.size()) *outS<<schm2[k]<<"\t";else *outS<<0<<"\t";

    outE=new ofstream(outfnameEig.data(),ios::app);  //*outE<<setprecision(10);
    *outE<<time<<"\t";
    for(int k=0;k<Dcut*Dcut;k++){
      //if(k<eigvals.size()) *outE<<abs(eigvals[k])<<"\t";else *outE<<0<<"\t";
      if(k<S.getDimension(0)) *outE<<abs(S.getElement(Indices(k,k)))<<"\t";else *outE<<0<<"\t";
    }
    *outE<<endl;
      outE->close();delete outE;
      // apply one step 
      //    applyTruncateMPO(W,lambdaR,midU,Dcut,tolSVD);
      int D=W.getDimension(1);
      int xi=midU.getDimension(1); // d xi d xi
      
      applyMPO(W,lambdaR,midU,tolSVD); // just apply the evolution, no cut
      // // applyTruncateMPO(W,lambdaR,midU,D*xi,tolSVD); // no cut
      //applyTruncateMPO(W,lambdaR,midU,Dcut,tolSVD); // no cut
      //    applyMPO(W,midU); // changes the MPS
      //truncateMPS(W,lambdaR,Dcut); // truncates with canonical form
    cnt++;
    time+=delta;
  }

  //out->close();outS->close();delete out;delete outS;
}

// Inefficient version
void getLargestEigenvalue(const mwArray& M,mwArray& V,complex_t& lambda){
  int D=M.getDimension(0);
  if(M.getDimension(1)!=D){
    cout<<"Error trying power method"<<endl;
    exit(1);
  }
  //cout<<"getLargesteigenvalue M:"<<M.getDimensions()<<endl;
  bool conv=0;
  V=mwArray(Indices(D,1));V.fillRandom();
  V=(1./norm(V))*V;
  int numIt=1;
  while(!conv&&numIt<10000){
    mwArray newV=M*V;
    mwArray aux=Hconjugate(newV)*V;
    complex_t new_lambda=aux.getElement(0);
    if(abs(new_lambda-lambda)/abs(lambda)<1E-8) conv=true;    
    V=(1./norm(newV))*newV;lambda=new_lambda;
    numIt++;
  }
  //  cout<<"Found largest eigenvalue of M"<<M.getDimensions()<<" to be "<<lambda<<endl;
}

void getLargestEigenvalue(TensorMultiplier& M,mwArray& V,complex_t& lambda){
  int D=M.getSize();
  //cout<<"getLargesteigenvalue M:"<<M.getSize()<<endl;
  bool conv=0;
  V=mwArray(Indices(D,1));V.fillRandom();
  V=(1./norm(V))*V;
  int numIt=1;
  while(!conv&&numIt<10000){
    mwArray newV;
    M.product(V,newV);
    mwArray aux=Hconjugate(newV)*V;
    complex_t new_lambda=aux.getElement(0);
    if(abs(new_lambda-lambda)/abs(lambda)<1E-13) conv=true;    
    V=(1./norm(newV))*newV;lambda=new_lambda;
    numIt++;
  }
  // cout<<"It nr "<<numIt<<". Found largest eigenvalue to be "<<lambda<<endl;
  // mwArray aux;M.product(V,aux);
  // cout<<"CHECK!! M*V-lambda*V"<<norm(aux-lambda*V)<<endl;
}

complex_t getFirstNonZeroDiag(const mwArray& X){
  bool done=false;
  int D=X.getDimension(0);
  int cnt=0;
  //cout<<"getFirstNonZeroDiag of X:"<<X.getDimensions()<<endl;
  while(!done&&cnt<D){
    complex_t res=X.getElement(Indices(cnt,cnt));
    //cout<<"Looked al element ("<<cnt<<")="<<res<<endl;
    if(abs(res)>1E-3) return res;
    //cout<<"Looked al element ("<<cnt<<") was not big enough"<<endl;
    cnt++;
  }
  if(!done){
    cout<<"ERROR! All diagonal elements of X are too small!"<<endl;
    putForMatlab(cout,X,"X");
    exit(1);
  }
}

void applyTruncateMPO(mwArray& W,mwArray& lambdaR,const mwArray& midU,int Dcut,double tolSVD){
  //cout<<"Apply MPO "<<midU.getDimensions()<<" to W:"<<W.getDimensions()<<" and truncate to "<<Dcut<<endl;
  int D=W.getDimension(1);
  int xi=midU.getDimension(1); // d xi d xi
  mwArray auxOp(midU); auxOp.reshape(Indices(d,xi*d*xi)); // d_u x
							  // xi_l d_d xi_r
  auxOp.multiplyLeft(Hconjugate(auxOp)); // xi*d*xi x xi*d*xi
  auxOp.reshape(Indices(xi,d,xi,xi,d,xi));
  auxOp.permute(Indices(4,1,5,6,3,2));
  auxOp.reshape(Indices(xi*xi,d,xi*xi,d));
  TensorMultiplier E(permute(W,Indices(2,1,3)),auxOp,permute(conjugate(W),Indices(2,1,3)));
  mwArray Xt,Yt,X,Y;
  complex_t etaL,etaR;
  getLargestEigenvalue(E,Xt,etaR);
  //cout<<"Found largest R eigenvalue "<<etaR<<endl;
  TensorMultiplier El(permute(conjugate(W),Indices(3,1,2)),permute(conjugate(auxOp),Indices(3,2,1,4)),permute(W,Indices(3,1,2)));
  getLargestEigenvalue(El,Yt,etaL);Yt.Hconjugate();  // 1 x xi*xi'*D*D'
  W=(ONE_c/sqrt(etaR))*W;
  //  cout<<"Found largest L eigenvalue "<<etaL<<", and <Yt|Xt>="<<Yt*Xt<<endl;
  Xt.reshape(Indices(xi,xi,D,D));
  Xt.permute(Indices(1,3,2,4));
  Xt.reshape(Indices(xi*D,xi*D));
  // This guy should be Hermitian, but the global phase may be
  // arbitrary from the algorithm=> I find it looking at the diagonal,
  // and remove it.
  complex_t elX=getFirstNonZeroDiag(Xt);
  complex_t phase=conjugate(elX)/abs(elX);
  Xt=phase*Xt;
  if(!Xt.isHermitian(1E-10)){
    cout<<"% ERROR!??! Xt not Hermitian (should be Hermitian times a phase!): "<<endl;
    // exp(i 2* phi) X = X^dagger
    /// gt a value which is not zero 
    //putForMatlab(cout,Xt,"Xt");
  }
  Xt=.5*(Xt+Hconjugate(Xt));
  //cout<<" Xt="<<Xt<<endl;
  //  cout<<"Empezando con Yt:"<<Yt.getDimensions()<<endl;
  Yt.reshape(Indices(xi,xi,D,D));
  //cout<<"reshaped: "<<Yt.getDimensions()<<endl;
  Yt.permute(Indices(2,4,1,3));
  //cout<<"permuted: "<<Yt.getDimensions()<<endl;
  Yt.reshape(Indices(xi*D,xi*D));
  //cout<<"reshaped: "<<Yt.getDimensions()<<endl;
  complex_t elY=getFirstNonZeroDiag(Yt);
  //cout<<"Found element in the diag "<<elY<<endl;
  phase=conjugate(elY)/abs(elY);
  //cout<<"phase to remove is "<<phase<<endl;
  Yt.multiplyLeft(phase);
  //cout<<"rephased: "<<Yt.getDimensions()<<endl;
  if(!Yt.isHermitian(1E-10)){
    cout<<"% ERROR!??! Yt not Hermitian: "<<endl;
    //putForMatlab(cout,Yt,"Yt");
  }
  Yt=.5*(Yt+Hconjugate(Yt));
//cout<<" Yt="<<Yt<<endl;
  complex_t trYX;
  {
    mwArray auxTr=Yt*Xt;
    trYX=auxTr.trace(); // norhtonormalize Yt*Xt
  }
  if(abs(imag(trYX))>1E-10){
    cout<<"ERROR!Why is the product of Y and X not a positive real number??"<<endl;
    exit(1);
  }
  Yt=1/sqrt(abs(trYX))*Yt;
  Xt=1/sqrt(abs(trYX))*Xt;
  vector<complex_t> Dx,Dy;
  mwArray Ux,Uy;
  wrapper::eig(Xt,Dx,Ux,true);
  wrapper::eig(Yt,Dy,Uy,true);
  // cout<<"Diagonalization of Xt:"<<Dx<<endl;
  //cout<<"Diagonalization of Yt:"<<Dy<<endl;
  if(real(Dx[0])<-1E12){
    cout<<"X seems to be all negative=> changing sign"<<endl;
    X=Ux*sqrt(-ONE_c*diag(Dx))*Hconjugate(Ux);
  }
  else X=Ux*sqrt(diag(Dx))*Hconjugate(Ux);
  if(real(Dy[0])<-1E12){
    cout<<"Y seems to be all negative=> changing sign"<<endl;
    Y=Uy*sqrt(-ONE_c*diag(Dy))*Hconjugate(Uy);
  }
  else Y=Uy*sqrt(diag(Dy))*Hconjugate(Uy);
  //cout<<"Created X:"<<X.getDimensions()<<" and Y:"<<Y.getDimensions()<<endl;
  mwArray U,S,V;
  wrapper::svd(Y*X,tolSVD,Dcut,U,S,V);
  lambdaR=S;
  //  mwArray valsS;S.getDiagonal(valsS);
  //cout<<"After svd of YX, S="<<valsS<<endl;
  W.reshape(Indices(d,D*D));
  W.multiplyLeft(reshape(permute(midU,Indices(1,2,4,3)),Indices(d*xi*xi,d))); //d*xi*xiXD*D
  W.reshape(Indices(d,xi,xi,D,D));
  W.permute(Indices(1,2,4,3,5));
  W.reshape(Indices(d*xi*D,xi*D));
  int nrS=0;
  mwArray auxS=invertDiag(S,nrS,1E-10);
  //  cout<<"nrS="<<nrS<<endl;
  W.multiplyRight(X*Hconjugate(V)*auxS);//invertDiag(S,nrS));
  //W.multiplyRight(X*Hconjugate(V)*invertDiag(S,nrS));
  W.reshape(Indices(d,xi*D,Dcut));
  W.permute(Indices(2,1,3));
  W.reshape(Indices(xi*D,d*Dcut));
  W.multiplyLeft(Hconjugate(U)*Y);
  W.reshape(Indices(Dcut,d,Dcut));
  W.permute(Indices(2,1,3));
}

void applyMPO_onlyW(mwArray& W,const mwArray& midU){
  //cout<<"Apply MPO "<<midU.getDimensions()<<" to W:"<<W.getDimensions()<<endl;
  int D=W.getDimension(1);
  int xi=midU.getDimension(1); // d xi d xi
  W.reshape(Indices(d,D*D));
  W.multiplyLeft(reshape(permute(midU,Indices(1,2,4,3)),Indices(d*xi*xi,d))); //d*xi*xiXD*D
  W.reshape(Indices(d,xi,xi,D,D));
  W.permute(Indices(1,2,4,3,5));
  W.reshape(Indices(d,xi*D,xi*D));
}

void applyMPO_WLambda(mwArray& W,mwArray& lambda,const mwArray& midU){
  cout<<"Apply MPO "<<midU.getDimensions()<<" to W:"<<W.getDimensions()<<", lambda="<<lambda.getDimensions()<<endl;
  int D=W.getDimension(1);
  int xi=midU.getDimension(1); // d xi d xi
  W.reshape(Indices(d,D*D));
  W.multiplyLeft(reshape(permute(midU,Indices(1,2,4,3)),Indices(d*xi*xi,d))); //d*xi*xiXD*D
  W.reshape(Indices(d,xi,xi,D,D));
  W.permute(Indices(1,2,4,3,5));
  W.reshape(Indices(d,xi*D,xi*D));
  lambda.reshape(Indices(1,D*D));
  lambda.multiplyLeft(reshape(identityMatrix(xi),Indices(xi*xi,1)));
  lambda.reshape(Indices(xi,xi,D,D));
  lambda.permute(Indices(1,3,2,4));
  lambda.reshape(Indices(xi*D,xi*D));
}



void computeSingleSiteRDM(const mwArray& W,const mwArray& lambda,mwArray &rdm){
  // assuming canonical form
  // W is d x D x D, lambda is Dx D
  //cout<<"computeSingleSiteRDM, W:"<<W.getDimensions()<<", lambda:"<<lambda.getDimensions()<<endl;
  rdm=W;
  int D=W.getDimension(1);
  rdm.reshape(Indices(d*D,D));
  rdm.multiplyRight(lambda);
  rdm.reshape(Indices(d,D*D));
  rdm.multiplyRight(Hconjugate(rdm)); // d x d
  rdm.transpose();
  rdm.reshape(Indices(d*d,1));
}


complex_t computeZZEV(const mwArray& W,const mwArray& lambda){
  //  cout<<"computeZZEV"<<endl;
  int D=W.getDimension(1);
  mwArray aux(W); // d x D x D
  aux.permute(Indices(2,1,3));
  aux.reshape(Indices(D*d,D));
  aux.multiplyRight(reshape(aux,Indices(D,d*D))); //D_l*d x d*D_r
  aux.reshape(Indices(D*d*d,D));
  aux.multiplyRight(lambda);
  aux.reshape(Indices(D,d*d,D));
  aux.permute(Indices(2,1,3));
  aux.reshape(Indices(d*d,D*D));
  aux.multiplyRight(Hconjugate(aux)); // d*d x d*d
  aux.reshape(Indices(d,d,d,d));
  aux.permute(Indices(3,1,4,2));
  aux.reshape(Indices(d*d,d*d));
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  aux.multiplyLeft(reshape(sigZ,Indices(1,d*d))); // should be
						  // transposed, but
						  // this one is sym!
  aux.multiplyRight(reshape(sigZ,Indices(d*d,1))); 
  return aux.getElement(0);
}

void computeSingleSiteRDM_WLambda(const mwArray& W,const mwArray& lambda,mwArray &rdm){
  // assuming canonical form
  // W is d x D x D, lambda is Dx D
  //cout<<"computeSingleSiteRDM, W:"<<W.getDimensions()<<", lambda:"<<lambda.getDimensions()<<endl;
  rdm=W;
  int D=W.getDimension(1);
  rdm.reshape(Indices(d*D,D));
  rdm.multiplyRight(lambda);
  rdm.multiplyRight(Hconjugate(lambda));
  rdm.reshape(Indices(d,D,D));
  rdm.permute(Indices(2,1,3));
  rdm.reshape(Indices(D,d*D));
  rdm.multiplyLeft(lambda);
  rdm.multiplyLeft(Hconjugate(lambda)); // D x d*D
  rdm.reshape(Indices(D,d,D));
  rdm.permute(Indices(2,1,3));
  rdm.reshape(Indices(d,D*D));
  rdm.multiplyRight(Hconjugate(reshape(W,Indices(d,D*D))));
  rdm.reshape(Indices(d*d,1));
}

double computeBlockEntanglement(const mwArray& W,const mwArray& lambda,int nB,vector<double>& schm,vector<double>& schm2,double& negVal){
  // Not very efficient, right now: get the transfer op and diagonalize the nB-th power exactly
  // 

  mwArray E=W;
  int D=W.getDimension(1);
  E.reshape(Indices(d*D,D));
  E.multiplyRight(lambda);
  E.reshape(Indices(d,D*D));
  E.multiplyLeft(Hconjugate(E)); // Dl'xDr' DlxDr
  E.reshape(Indices(D,D,D,D));
  E.permute(Indices(1,2,4,3)); // Dl'xDr' DrxDl (just for convenience)
  //  E.permute(Indices(1,3,2,4));
  //E.reshape(Indices(D*D,D*D)); // Dl' Dl x Dr' Dr
  for(int k=1;k<nB;k++){ // could group terms to save time
    // multiply in two phases W and W^dagger, both on the left
    E.reshape(Indices(D,D*D*D)); // Dl'x Dr' Dr Dl
    mwArray Eb(W);Eb.reshape(Indices(d*D,D));
    E.multiplyLeft(conjugate(Eb)); // d*Dl'x Dr' Dr Dl
    E.reshape(Indices(d,D*D*D,D));
    E.permute(Indices(1,3,2));
    E.reshape(Indices(d*D,D*D*D));
    Eb.reshape(Indices(d,D,D));Eb.permute(Indices(2,1,3));
    Eb.reshape(Indices(D,d*D));
    E.multiplyLeft(Eb); // Dl x Dl' Dr' Dr
    E.reshape(Indices(D,D,D,D));
    E.permute(Indices(2,3,4,1)); // Dl'xDr' DrxDl (just for convenience)
  }
  E.permute(Indices(1,4,2,3)); // Dl' Dl Dr' Dr
  E.reshape(Indices(D*D,D*D)); // Dl' Dl x Dr' Dr
 
  // // // Insert the lambda from the left (NO! from left they are orthogonal)
  // // mwArray auxLamb(conjugate(lambda));auxLamb.reshape(Indices(D*D,1));
  // // auxLamb.multiplyRight(Hconjugate(auxLamb));
  // // auxLamb.reshape(Indices(D,D,D,D));
  // // auxLamb.permute(Indices(1,3,2,4));
  // // auxLamb.reshape(Indices(D*D,D*D));
  // // Eb.multiplyLeft(auxLamb);
  // Now Eb has the tensor for the block: diagonalize
  mwArray U,V,S;
  int Dcut;
  wrapper::svd(E,U,S,V); //tolSVD,Dcut,U,S,V);
  // cout<<"After SVD S="<<endl;
  // putForMatlab(cout,S,"S");

  // If computing entropy
  double Sval=0;
  //  for(int k=0;k<Dcut;k++){
  for(int k=0;k<S.getDimension(0);k++){
    double valK=real(S.getElement(Indices(k,k)));
    schm.push_back(valK);
    if(abs(valK)>1E-10)
      Sval-=valK*log2(valK);
  }

  // If computing negativity
  mwArray rhoPT(E);rhoPT.reshape(Indices(D,D,D,D)); // Dl' Dl Dr' Dr
  rhoPT.permute(Indices(1,4,2,3)); // PT and regroup Dl' Dr x Dl Dr'
  rhoPT.reshape(Indices(D*D,D*D));
  wrapper::svd(rhoPT,U,S,V);
  negVal=log2(real(S.trace()));

  E.reshape(Indices(D,D,D,D));
  E.permute(Indices(1,3,2,4)); // Dl' Dr' x Dl Dr (vertical)
  E.reshape(Indices(D*D,D*D));
  wrapper::svd(E,U,S,V); //tolSVD,Dcut,U,S,V);
  for(int k=0;k<S.getDimension(0);k++){
    double valK=real(S.getElement(Indices(k,k)));
    schm2.push_back(valK);
  }  
  return Sval;

  
}

void computeBlockDecompositions(const mwArray& W,const mwArray& lambda,int nB,
				vector<complex_t>& valsPh,
				vector<double>& schm){
  // Not very efficient, right now: get the transfer op and diagonalize the nB-th power exactly

  mwArray E=W;
  int D=W.getDimension(1);
  E.reshape(Indices(d*D,D));
  E.multiplyRight(lambda);
  E.reshape(Indices(d,D*D));
  E.multiplyLeft(Hconjugate(E)); // Dl'xDr' DlxDr
  E.reshape(Indices(D,D,D,D));
  E.permute(Indices(1,3,2,4));
  E.reshape(Indices(D*D,D*D)); // Dl' Dl x Dr' Dr
  mwArray Eb=E;
  for(int k=1;k<nB;k++){ // could group terms to save time
    Eb.multiplyLeft(E);
  }
  Eb.reshape(Indices(D,D,D,D));
  Eb.permute(Indices(1,3,2,4));
  Eb.reshape(Indices(D*D,D*D));
  //vector<complex_t> Deb;
  mwArray Ueb;valsPh.clear();
  wrapper::eig(Eb,valsPh,Ueb,true);
  cout<<"Diagonalized, eigenvals "<<valsPh<<endl;  
  // Now Eb has the tensor for the block: diagonalize
  mwArray U,V,S;
  int Dcut;
  wrapper::svd(Eb,U,S,V); //tolSVD,Dcut,U,S,V);
  // cout<<"After SVD S="<<endl;
  // putForMatlab(cout,S,"S");
  double Sval=0;
  //  for(int k=0;k<Dcut;k++){
  for(int k=0;k<S.getDimension(0);k++){
    double valK=real(S.getElement(Indices(k,k)));
    schm.push_back(valK);
    if(abs(valK)>1E-10)
      Sval-=valK*log2(valK);
  }
  //return Sval;
}


complex_t computeZZEV_WLambda(const mwArray& W,const mwArray& lambda){
  //  cout<<"computeZZEV"<<endl;
  int D=W.getDimension(1);
  mwArray aux(W); // d x D x D
  aux.permute(Indices(2,1,3));
  aux.reshape(Indices(D,d*D));
  aux.multiplyLeft(lambda);
  aux.reshape(Indices(D,d,D));
  aux.permute(Indices(2,1,3));
  aux.reshape(Indices(d*D,D));
  aux.multiplyRight(reshape(permute(W,Indices(2,1,3)),Indices(D,d*D))); // d_l*D_lx d_r*D_r
  aux.reshape(Indices(d*D*d,D));
  aux.multiplyRight(lambda);
  aux.reshape(Indices(d,D,d,D));
  aux.permute(Indices(1,3,2,4));
  aux.reshape(Indices(d*d,D*D));
  aux.multiplyRight(Hconjugate(aux)); // d*d_s x d*d_a
  aux.reshape(Indices(d,d,d,d));
  aux.permute(Indices(1,3,2,4));
  aux.reshape(Indices(d*d,d*d));
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  aux.multiplyLeft(reshape(sigZ,Indices(1,d*d))); // should be
						  // transposed, but
						  // this one is sym!
  aux.multiplyRight(reshape(sigZ,Indices(d*d,1))); 
  return aux.getElement(0);
}


void constructInitialState(int initSt,mwArray& W){
  complex_t data0[]={ONE_c,ZERO_c};
  complex_t data1[]={ZERO_c,ONE_c};
  mwArray v0(Indices(d,1),data0);
  mwArray v1(Indices(d,1),data1);
  
  switch(initSt){
  case 1: // X+
    {
      W=(1./sqrt(2))*(v0+v1);
      break;
    }
  case 2: // Y+
    {
      W=(1./sqrt(2))*(v0+I_c*v1);
      break;
    }
  case 3: // Z+
    {
      W=v0;
      break;
    }
  default:
    cout<<"Error: unknown type of intSt="<<initSt<<endl;
    exit(1);
  }
  W.reshape(Indices(d,1,1));
}

void truncateMPS_WLambda(mwArray& W,mwArray& lambda,int Dcut,double tolSVD){
  int D=W.getDimension(1);
  // Find left and right dominant eigenvectors
  mwArray X,Y;
  complex_t etaL,etaR;
  // Could be more efficient, but right now, brute force
  // TODO: Implement the product as TensorMultiplier, and run the
  // power method more efficiently (Lanczos?)
  mwArray E(W); // transfer matrix
  E.reshape(Indices(d,D*D));
  E.multiplyLeft(Hconjugate(E)); // D*D,D*D
  E.reshape(Indices(D,D,D,D));
  E.permute(Indices(1,3,2,4));
  E.reshape(Indices(D*D,D*D));
  mwArray lambdaE(lambda);
  lambdaE.reshape(Indices(1,D*D));
  lambdaE.multiplyLeft(Hconjugate(lambdaE)); // D*D, D*D
  lambdaE.reshape(Indices(D,D,D,D));
  lambdaE.permute(Indices(1,3,2,4));
  lambdaE.reshape(Indices(D*D,D*D));
  getLargestEigenvalue(E*lambdaE,X,etaR);
  cout<<"Found largest R eigenvalue "<<etaR<<endl;
  getLargestEigenvalue(Hconjugate(lambdaE*E),Y,etaL);
  cout<<"Found largest L eigenvalue "<<etaL<<endl;
  W=(ONE_c/sqrt(etaR))*W;
  //  if(D<=Dcut) return;
  cout<<"truncateMPS W:"<<W.getDimensions()<<", lambda:"<<lambda.getDimensions()<<" to Dcut="<<Dcut<<endl;
  Y.Hconjugate(); // 1 x D*D
  X.reshape(Indices(D,D));X.transpose();
  X=X+Hconjugate(X);
  //  cout<<" with X-Xdagger="<<X-Hconjugate(X)<<endl;
  Y.reshape(Indices(D,D));
  Y=Y+Hconjugate(Y);
  //cout<<" with Y-Ydagger="<<Y-Hconjugate(Y)<<endl;
  // Need the sqr of X and Y => ful diagonalization
  vector<complex_t> Dx,Dy;
  mwArray Ux,Uy;
  wrapper::eig(X,Dx,Ux,true);
  wrapper::eig(Y,Dy,Uy,true);
  //cout<<"X-Ux*Dx*Uxdagger="<<X-Ux*diag(Dx)*Hconjugate(Ux)<<endl;
  mwArray sqY=Uy*sqrt(diag(Dy));
  mwArray sqX=Ux*sqrt(diag(Dx));
  mwArray U,V;
  mwArray aux=Hconjugate(sqY*lambda*sqX);
  wrapper::svd(aux,tolSVD,Dcut,U,lambda,V);
  int nrX=0; // nr of eigs of invX
  mwArray invX=inverse(X,nrX);
  int nrY=0; // nr of eigs of invYdagger
  mwArray invY=inverse(Hconjugate(Y),nrY);
  W.reshape(Indices(d*D,D));
  W.multiplyRight(invY);
  W.multiplyRight(U);
  W.reshape(Indices(d,D,D));
  W.permute(Indices(2,1,3));
  W.reshape(Indices(D,d*D));
  W.multiplyLeft(invX);
  W.multiplyLeft(V);
  W.reshape(Indices(D,d,D));
  W.permute(Indices(2,1,3));
  mwArray norm=lambda*Hconjugate(lambda);
  lambda=(ONE_c/sqrt(norm.trace()))*lambda;
}

void applyMPO(mwArray& W,mwArray& lambdaR,const mwArray& midU,double tolSVD){
  applyMPO_onlyW(W,midU);
  applyCanonical(W,lambdaR,tolSVD);
}


void applyCanonical(mwArray& W,mwArray& lambdaR,double tolSVD){
  int d=W.getDimension(0); // effective phys dim, here d*d
  int D=W.getDimension(1);
  int Dcut=D; // no truncation
  mwArray auxOp=identityMatrix(d);auxOp.reshape(Indices(1,d,1,d)); // trick for TensorMultiplier
  TensorMultiplier E(permute(W,Indices(2,1,3)),auxOp,permute(conjugate(W),Indices(2,1,3)));
  mwArray Xt,Yt,X,Y;
  complex_t etaL,etaR;
  getLargestEigenvalue(E,Xt,etaR);
  // cout<<"Found largest R eigenvalue "<<etaR<<endl;
  TensorMultiplier El(permute(conjugate(W),Indices(3,1,2)),permute(conjugate(auxOp),Indices(3,2,1,4)),permute(W,Indices(3,1,2)));
  getLargestEigenvalue(El,Yt,etaL);Yt.Hconjugate();  // 1 x xi*xi'*D*D'
  W=(ONE_c/sqrt(etaR))*W;
  // cout<<"Found largest L eigenvalue "<<etaL<<", and <Yt|Xt>="<<Yt*Xt<<endl;
  Xt.reshape(Indices(D,D));
  // This guy should be Hermitian, but the global phase may be
  // arbitrary from the algorithm=> I find it looking at the diagonal,
  // and remove it.
  complex_t elX=getFirstNonZeroDiag(Xt);
  complex_t phase=conjugate(elX)/abs(elX);
  Xt=phase*Xt;
  if(!Xt.isHermitian(1E-10)){
    cout<<"% ERROR!??! Xt not Hermitian (should be Hermitian times a phase!): "<<endl;
    // exp(i 2* phi) X = X^dagger
    /// gt a value which is not zero 
    //putForMatlab(cout,Xt,"Xt");
  }
  Xt=.5*(Xt+Hconjugate(Xt));
  //cout<<" Xt="<<Xt<<endl;
  // cout<<"Empezando con Yt:"<<Yt.getDimensions()<<endl;
  Yt.reshape(Indices(D,D));
  Yt.permute(Indices(2,1));
  complex_t elY=getFirstNonZeroDiag(Yt);
  phase=conjugate(elY)/abs(elY);
  Yt.multiplyLeft(phase);
  if(!Yt.isHermitian(1E-10)){
    cout<<"% ERROR!??! Yt not Hermitian: "<<endl;
    //putForMatlab(cout,Yt,"Yt");
  }
  Yt=.5*(Yt+Hconjugate(Yt));
  complex_t trYX;
  {
    mwArray auxTr=Yt*Xt;
    trYX=auxTr.trace(); // orhtonormalize Yt*Xt
  }
  if(abs(imag(trYX))>1E-10){
    cout<<"ERROR!Why is the product of Y and X not a positive real number??"<<endl;
    exit(1);
  }
  Yt=1/sqrt(abs(trYX))*Yt;
  Xt=1/sqrt(abs(trYX))*Xt;
  vector<complex_t> Dx,Dy;
  mwArray Ux,Uy;
  wrapper::eig(Xt,Dx,Ux,true);
  wrapper::eig(Yt,Dy,Uy,true);
  if(real(Dx[0])<-1E12){
    cout<<"X seems to be all negative=> changing sign"<<endl;
    X=Ux*sqrt(-ONE_c*diag(Dx))*Hconjugate(Ux);
  }
  else X=Ux*sqrt(diag(Dx))*Hconjugate(Ux);
  if(real(Dy[0])<-1E12){
    cout<<"Y seems to be all negative=> changing sign"<<endl;
    Y=Uy*sqrt(-ONE_c*diag(Dy))*Hconjugate(Uy);
  }
  else Y=Uy*sqrt(diag(Dy))*Hconjugate(Uy);
  //cout<<"Created X:"<<X.getDimensions()<<" and Y:"<<Y.getDimensions()<<endl;
  mwArray U,S,V;
  wrapper::svd(Y*X,tolSVD,Dcut,U,S,V);
  lambdaR=S;
  mwArray valsS;S.getDiagonal(valsS);
  //cout<<"After svd of YX, S="<<valsS<<endl;
  W.reshape(Indices(d*D,D));
  int nrS=0;
  mwArray auxS=invertDiag(S,nrS,1E-10);
  //  cout<<"nrS="<<nrS<<endl;
  W.multiplyRight(X*Hconjugate(V)*auxS);//invertDiag(S,nrS));
  //W.multiplyRight(X*Hconjugate(V)*invertDiag(S,nrS));
  W.reshape(Indices(d,D,Dcut));
  W.permute(Indices(2,1,3));
  W.reshape(Indices(D,d*Dcut));
  W.multiplyLeft(Hconjugate(U)*Y);
  W.reshape(Indices(Dcut,d,Dcut));
  W.permute(Indices(2,1,3));
}




void computeBlockPartTrace(const mwArray& W,const mwArray& lambda,int nB,vector<complex_t>& eigvals){
  // Now trace out the nB sites and boundary and try to guess the product structure?
  // Assuming canonical with Id on left and lambda on right, I'll look at the trace from the right

  int d=W.getDimension(0); // phys dim, here d 
  int D=W.getDimension(1);
  int Dcut=D; // no truncation
  mwArray auxOp=identityMatrix(d);auxOp.reshape(Indices(1,d,1,d)); // trick for TensorMultiplier
  TensorMultiplier E(permute(W,Indices(2,1,3)),auxOp,permute(conjugate(W),Indices(2,1,3)));

  mwArray rightPart=identityMatrix(D);rightPart.reshape(Indices(D,1,D));
  for(int k=0;k<nB;k++){
    mwArray aux(rightPart);
    E.product(aux,rightPart);
  }
  
  rightPart.reshape(Indices(D,D));
  rightPart=.5*(rightPart+Hconjugate(rightPart));
  
  mwArray Ux;
  wrapper::eig(rightPart,eigvals,Ux,false);
  

}

void specDoubleBlockLeft(mwArray& W,mwArray& lambda,int nB,int Dcut,mwArray& S,mwArray& U,mwArray& V){
  // Try GILD operations on the rdm to determine the projector that kills a CDL ???
  // First try if there is one, with the SVD
  // I need the RDM for 2*nB, with left side open => lambda squared on the right, and then copies of W
  if(nB>2){
    cout<<"WARNING!!! The tmp matrix for this step may be too big!!!"<<endl;    
  }
  // First put together 2*nB copies of W
  int d=W.getDimension(0); // phys dim, here d 
  int D=W.getDimension(1);

  mwArray aux=W;
  aux.permute(Indices(2,1,3));
  mwArray tmpW(aux);tmpW.reshape(Indices(D,d*D));
  aux.reshape(Indices(D*d,D));
  int dphys=d;
  for(int k=1;k<2*nB;k++){
    aux.multiplyRight(tmpW); // D*d^(k),d*D
    dphys*=d;
    aux.reshape(Indices(D*dphys,D));
  }
  aux.multiplyRight(lambda); // D*dphys,D
  aux.multiplyRight(Hconjugate(aux)); // D*dphys x D*dphys
  aux.reshape(Indices(D,dphys,D,dphys));
  aux.permute(Indices(1,3,2,4));
  aux.reshape(Indices(D*D,dphys*dphys));
  wrapper::svd(aux,U,S,V);    
}


void findProjectorRprimeLeft(mwArray& U,mwArray&S,int D,double eps,mwArray& Ur,vector<complex_t>& eigVals){
  mwArray trId=identityMatrix(sqrt(U.getDimension(0)));
  trId.reshape(Indices(1,U.getDimension(0)));
  trId.multiplyRight(U); // Elems are t_i
  // Now reweight the elements with those of S^2/(eps+S^2)
  mwArray auxS2(S);
  auxS2.multiplyRight(S);
  mwArray auxDen=auxS2+eps*identityMatrix(U.getDimension(1));
  int nr;
  auxDen=invertDiag(auxDen,nr,eps);
  auxS2.multiplyRight(auxDen);
  // for(int k=0;k<D*D;k++){
  //   double S2=real(S.getElement(Indices(k,k)));
  //   S2=S2*S2;
  //   S2=S2/(eps+S2);
  //   auxS2.setElement(S2*ONE_c,Indices(k,k));
  // }
  trId.multiplyRight(auxS2); // -> the t'_i
  //S=trId;//return;
  trId.multiplyRight(Hconjugate(U)); //1x D*D -> R'
  trId.reshape(Indices(sqrt(U.getDimension(0)),sqrt(U.getDimension(0))));
  //cout<<"R'="<<trId<<endl;
  // In this case, R' is hermitian, so I will correct for numerics
  trId=.5*(trId+Hconjugate(trId));
  wrapper::eig(trId,eigVals,Ur,true);
  //  cout<<"R'->"<<eigVals<<endl;
}


// This one does nor work
void truncDoubleBlockLeft(mwArray& W,mwArray& lambda,int nB,int Dcut,mwArray& S,double eps){
  int d=W.getDimension(0);int D=W.getDimension(1);
  mwArray U,V,Ur;
  mwArray Wcut(W);Wcut.permute(Indices(2,1,3)); // D*d*D
  specDoubleBlockLeft(W,lambda,nB,D,S,U,V);  //the SVD of environment ;
  if(D>Dcut){
    cout<<"Attempting a truncation, from D="<<D<<" dim to Dcut="<<Dcut<<endl;
  }
  int cnt=0;
  vector<complex_t> eigVals;
  while(cnt<maxIter){
    findProjectorRprimeLeft(U,S,D,eps,Ur,eigVals);
    Ur.multiplyRight(sqrt(diag(eigVals)));
    if(cnt==maxIter-1){// cut eigs below precCut (if below Dcut)
      bool* keep=new bool[D];int howmany=0;
      for(int p=0;p<D;p++){
	keep[p]=abs(eigVals[p])>=precCut?1:0;
	if(keep[p]) howmany++;
      }
      mwArray projMat(Indices(D,howmany));int pos=0;
      for(int p=0;p<D;p++){
	if(keep[p]){
	  projMat.setElement(ONE_c,Indices(p,pos));
	  pos++;
	}
      }
      delete[] keep;
      Ur.multiplyRight(projMat);
    } 


    // Apply the projector to U (on both legs) and repeat (unless it is the last one)
    int Dn=Ur.getDimension(1);
    if(cnt<maxIter-1){
      U.reshape(Indices(D,D*D*D));U.multiplyLeft(permute(Ur,Indices(2,1)));
      U.reshape(Indices(D,D,D*D));U.permute(Indices(1,3,2));U.reshape(Indices(D*D*D,D));
      U.multiplyRight(conjugate(Ur));
      U.reshape(Indices(D,D*D,D));U.permute(Indices(1,3,2));U.reshape(Indices(D*D,D*D));
    }
    // The same in W!
    //cout<<"Iter "<<cnt<<" Wcut="<<Wcut<<", Ur="<<Ur<<endl;
    Wcut.reshape(Indices(D,d*D));Wcut.multiplyLeft(Hconjugate(Ur));
    Wcut.reshape(Indices(Dn*d,D));Wcut.multiplyRight(Ur);Wcut.reshape(Indices(Dn,d,Dn));    
    cnt++;
  }
  S=diag(eigVals);
  //  cout<<" After the proc, W was "<<W<<", and now Wcut="<<Wcut<<", Ur="<<Ur<<endl;
  Wcut.permute(Indices(2,1,3));
  W=Wcut;
  applyCanonical(W,lambda,tolSVD);
  //cout<<" After applyCanonical W="<<W<<endl;
  
  return;
  //mwArray Ur,Sr,Vr;
  //wrapper::svd(trId,Ur,Sr,Vr);    
} 

void specDoubleBlock(mwArray& W,mwArray& lambda,int nB,int Dcut,mwArray& S,mwArray& U,mwArray& V){
  // Try GILD operations on the rdm to determine the projector that kills a CDL ???
  // First try if there is one, with the SVD
  // I need the RDM for 2*nB, open in the middle => lambda squared on the right, and then copies of W
  if(nB>2){
    cout<<"WARNING!!! The tmp matrix for this step may be too big!!!"<<endl;    
  }
  // First put together 2*nB copies of W
  int d=W.getDimension(0); // phys dim, here d 
  int D=W.getDimension(1);

  mwArray aux=W;
  aux.permute(Indices(2,1,3)); // D x d x D
  mwArray tmpW(aux);tmpW.reshape(Indices(D,d*D));
  aux.reshape(Indices(D*d,D));
  int dphys=d;
  for(int k=1;k<nB;k++){
    aux.multiplyRight(tmpW); // D*d^(k),d*D
    dphys*=d;
    aux.reshape(Indices(D*dphys,D));
  }
  // For the left half, contract with Id
  mwArray auxL(aux); auxL.reshape(Indices(D,dphys*D));
  auxL.multiplyLeft(Hconjugate(auxL)); // dphys*D x dphys*D
  auxL.reshape(Indices(dphys,D,dphys*D));
  auxL.permute(Indices(3,1,2));
  auxL.reshape(Indices(dphys*D*dphys,D)); // dphys_u*Du*dphys_d x D_d
    
  aux.multiplyRight(lambda); // D*dphys,D
  aux.multiplyRight(Hconjugate(aux)); // D*dphys x D*dphys
  aux.reshape(Indices(D*dphys,D,dphys)); //Du*dphys_u,Dd,dphys_d
  aux.permute(Indices(2,1,3)); // Dd x Du*dphys_u * dphys_d
  aux.reshape(Indices(D,D*dphys*dphys));
  aux.multiplyLeft(auxL); // dphys_u*Du*dphys_d x Du*dphys_u*dphys_d
  aux.reshape(Indices(dphys,D,dphys,D,dphys,dphys));
  aux.permute(Indices(2,4,1,3,5,6)); // Dl Dr x dlu dld dru drd
  aux.reshape(Indices(D*D,-1));
  
  wrapper::svd(aux,U,S,V);
  //mwArray enV;S.getDiagonal(enV);cout<<"After environment SVD, Svals:"<<enV<<endl;
}

void findProjectorRprime(mwArray& U,mwArray&S,int D,double eps,mwArray& Ur,mwArray& Sr,mwArray& Vr){
  mwArray trId=identityMatrix(D);
  trId.reshape(Indices(1,D*D));
  trId.multiplyRight(U); // Elems are t_i
  // Now reweight the elements with those of S^2/(eps+S^2)
  //cout<<"S="<<S<<endl;
  mwArray auxS2(S);
  //cout<<"S="<<auxS2<<endl;
  auxS2.multiplyRight(S);
  //cout<<"S2="<<auxS2<<endl;
  mwArray auxDen=auxS2+eps*identityMatrix(D*D);
  int nr;
  auxDen=invertDiag(auxDen,nr,eps*eps);
  //cout<<"1/(eps+S2)="<<auxDen<<endl;
  auxS2.multiplyRight(auxDen);
  //mwArray weights;auxS2.getDiagonal(weights);cout<<"S2/(eps+S2)="<<weights<<endl;
  // for(int k=0;k<D*D;k++){
  //   double S2=real(S.getElement(Indices(k,k)));
  //   S2=S2*S2;
  //   S2=S2/(eps+S2);
  //   auxS2.setElement(S2*ONE_c,Indices(k,k));
  // }
  trId.multiplyRight(auxS2); // -> the t'_i
  //cout<<"tprime="<<trId<<endl;
  //S=trId;//return;
  trId.multiplyRight(Hconjugate(U)); //1x D*D -> R'
  trId.reshape(Indices(D,D));
  //  cout<<"error (R-R')*S="<<reshape(identityMatrix(D)-trId,Indices(1,D*D))*S<<endl;
  // In this case, R' is not necessarily hermitian
  wrapper::svd(trId,Ur,Sr,Vr);    
  //mwArray auxD;Sr.getDiagonal(auxD);cout<<"R'->"<<auxD<<endl;
}


void truncDoubleBlock(mwArray& W,mwArray& lambda,int nB,int Dcut,mwArray& S,double eps){
  int d=W.getDimension(0);int D=W.getDimension(1);
  mwArray U,V,Ur;
  mwArray Wcut(W);Wcut.permute(Indices(2,1,3)); // D*d*D
  specDoubleBlock(W,lambda,nB,D,S,U,V);  //the SVD of environment ;
  cout<<"Attempting a truncation, from D="<<D<<endl;
  int cnt=0;
  //vector<complex_t> eigVals;
  mwArray Sr,Vr;
  while(cnt<maxIter){
    findProjectorRprime(U,S,D,eps,Ur,Sr,Vr);
    Ur.multiplyRight(sqrt(Sr));
    Vr.multiplyLeft(sqrt(Sr));
    if(cnt==maxIter-1){// cut eigs below precCut (if below Dcut)
      bool* keep=new bool[D];int howmany=0;
      for(int p=0;p<D;p++){
	keep[p]=abs(Sr.getElement(Indices(p,p)))>=precCut?1:0;
	if(keep[p]) howmany++;
      }
      mwArray projMat(Indices(D,howmany));int pos=0;
      for(int p=0;p<D;p++){
	if(keep[p]){
	  projMat.setElement(ONE_c,Indices(p,pos));
	  pos++;
	}
      }
      delete[] keep;
      if(howmany<D){
	cout<<"Reducing the dimension in the link to "<<howmany<<endl;
      }
      Ur.multiplyRight(projMat);
      Vr.multiplyLeft(Hconjugate(projMat));
    } 


    // Apply the projector to U (on both legs) and repeat (unless it is the last one)
    int Dn=Ur.getDimension(1);
    if(cnt<maxIter-1){
      U.reshape(Indices(D,D*D*D));U.multiplyLeft(permute(Ur,Indices(2,1)));
      U.reshape(Indices(D,D,D*D));U.permute(Indices(1,3,2));U.reshape(Indices(D*D*D,D));
      U.multiplyRight(conjugate(Ur));
      U.reshape(Indices(D,D*D,D));U.permute(Indices(1,3,2));U.reshape(Indices(D*D,D*D));
    }
    // The same in W!
    //cout<<"Iter "<<cnt<<" Wcut="<<Wcut<<", Ur="<<Ur<<endl;
    Wcut.reshape(Indices(D,d*D));Wcut.multiplyLeft(Vr); //Hconjugate(Vr));
    Wcut.reshape(Indices(Dn*d,D));Wcut.multiplyRight(Ur);Wcut.reshape(Indices(Dn,d,Dn));    
    cnt++;
  }
  //S=diag(eigVals);
  //  cout<<" After the proc, W was "<<W<<", and now Wcut="<<Wcut<<", Ur="<<Ur<<endl;
  Wcut.permute(Indices(2,1,3));
  W=Wcut;
  applyCanonical(W,lambda,tolSVD);
  //cout<<" After applyCanonical W="<<W<<endl;
  
  return;
  //mwArray Ur,Sr,Vr;
  //wrapper::svd(trId,Ur,Sr,Vr);    
} 


void specBlockLeft(mwArray& W,mwArray& lambda,int nB,int Dcut,mwArray& S,mwArray& U,mwArray& V){
  // Try GILD operations on the rdm to determine the projector that kills a CDL ???
  // First try if there is one, with the SVD
  // I need the RDM for nB, with left side open => lambda squared on the right, and then copies of W
  if(nB>4){
    cout<<"WARNING!!! The tmp matrix for this step may be too big!!!"<<endl;    
  }
  // First put together nB copies of W
  int d=W.getDimension(0); // phys dim, here d 
  int D=W.getDimension(1);

  mwArray aux=W;
  aux.permute(Indices(2,1,3));
  mwArray tmpW(aux);tmpW.reshape(Indices(D,d*D));
  aux.reshape(Indices(D*d,D));
  int dphys=d;
  for(int k=1;k<nB;k++){
    aux.multiplyRight(tmpW); // D*d^(k),d*D
    dphys*=d;
    aux.reshape(Indices(D*dphys,D));
  }
  aux.multiplyRight(lambda); // D*dphys,D
  aux.multiplyRight(Hconjugate(aux)); // D*dphys x D*dphys
  aux.reshape(Indices(D,dphys,D,dphys));
  aux.permute(Indices(1,3,2,4));
  aux.reshape(Indices(D*D,dphys*dphys));
  wrapper::svd(aux,U,S,V);    
}

void truncSingleBlock(mwArray& W,mwArray& lambda,int nB,int Dcut,mwArray& S,double eps){ //,bool hardCut){
  //cout<<"truncSingleBlock: W:"<<W.getDimensions()<<endl;
  int d=W.getDimension(0);int D=W.getDimension(1);
  mwArray U,V,Ur;
  mwArray Wcut(W);Wcut.permute(Indices(2,1,3)); // D*d*D
  specBlockLeft(W,lambda,nB,D,S,U,V);  //the SVD of environment ;
  //mwArray auxD;S.getDiagonal(auxD);cout<<"D="<<D<<", SVD env ->"<<auxD<<endl; 
  //cout<<"Attempting a truncation, from D="<<D<<endl;

  mwArray trId=identityMatrix(D);
  trId.reshape(Indices(1,D*D));
  trId.multiplyRight(U); // Elems are t_i
  // Now reweight the elements with those of (S^2/(eps+S^2))^maxIter
  //mwArray auxS2(S);
  for(int k=0;k<S.getDimension(0);k++){
    double val=real(S.getElement(Indices(k,k)));
    //auxS2.setElement(pow(val*val/(val*val+eps),maxIter)*ONE_c,Indices(k,k));
    trId.setElement(pow(val*val/(val*val+eps),maxIter)*trId.getElement(Indices(0,k)),Indices(0,k));
  }
    
  //trId.multiplyRight(auxS2*Hconjugate(U)); // This is R'
  trId.multiplyRight(Hconjugate(U)); // This is R'
  trId.reshape(Indices(D,D)); // Du x Dd
  trId=.5*(trId+Hconjugate(trId)); // TODO: check hermiticity
  vector<complex_t> eigVals;
  wrapper::eig(trId,eigVals,Ur,true); // R'= Ur eigs Ur' : Du x Dd
  //cout<<"R'->"<<eigVals<<endl;
  //cout<<"Ur="<<Ur<<endl;
  //Ur.multiplyRight(sqrt(diag(eigVals)));
  bool* keep=new bool[D];int howmany=0;
  for(int p=0;p<D;p++){
    //if(!hardCut)
    keep[p]=abs(eigVals[p])>=precCut?1:0;
    //    else keep[p]=(abs(eigVals[p])>=precCut&&howmany<Dcut)?1:0;// cut at max Dcut
    if(keep[p]) howmany++;
  }
  mwArray projMat(Indices(D,howmany));int pos=0;
  for(int p=0;p<D;p++){
    if(keep[p]){
      projMat.setElement(ONE_c,Indices(p,pos));
      pos++;
    }
  }
  delete[] keep;
  if(howmany<D){
    if(2-D/howmany>.1) return;
    cout<<"Reducing the dimension in the link from "<<D<<" to "<<howmany<<" ************ "<<D*1./howmany<<endl;    
    //    cout<<"Reducing the dimension in the link from "<<D<<" to "<<howmany<<endl;
  }
  Ur.multiplyRight(projMat);
  //cout<<"New Ur:"<<Ur<<endl;

  // Apply the projector to U (on both legs) and repeat (unless it is the last one)
  int Dn=Ur.getDimension(1);
  // The same in W!
  //cout<<"Iter "<<cnt<<" Wcut="<<Wcut<<", Ur="<<Ur<<endl;
  //  Wcut.reshape(Indices(D,d*D));Wcut.multiplyLeft(Hconjugate(Ur)); //Hconjugate(Vr));
  // Wcut.reshape(Indices(Dn*d,D));Wcut.multiplyRight(Ur);Wcut.reshape(Indices(Dn,d,Dn));    
  Wcut.reshape(Indices(D,d*D));Wcut.multiplyLeft(permute(Ur,Indices(2,1))); //Hconjugate(Vr));
  Wcut.reshape(Indices(Dn*d,D));Wcut.multiplyRight(conjugate(Ur));Wcut.reshape(Indices(Dn,d,Dn));    

  //S=diag(eigVals);
  //  cout<<" After the proc, W was "<<W<<", and now Wcut="<<Wcut<<", Ur="<<Ur<<endl;
  Wcut.permute(Indices(2,1,3));
  W=Wcut;
  applyCanonical(W,lambda,tolSVD);
  //cout<<" After applyCanonical W="<<W<<endl;
  
  return;
  //mwArray Ur,Sr,Vr;
  //wrapper::svd(trId,Ur,Sr,Vr);    
} 





void disentangleLegs(mwArray& W,mwArray& lambda,int nB,int Dcut,mwArray& Ufast,mwArray& Vfast){
  // En realidad me gustaría determinar numéricamente qué dof (y de
  // qué dimensión) se han desacoplado. No tiene por qué ser siempre
  // D=2, ni ocurrirá al ppio.
  // Si sé las dimensiones, puedo iterativamente tratar de optimizar el disentangler left y right.
  
  
  

}

// Optimize a unitry on legs leg1 and leg2=leg1+1, such that the rank of the cut between both is minimized
void optimizeDisentangler(const mwArray& state,int leg1,mwArray& Udis){
  Indices dims=state.getDimensions();
  if(leg1<0||leg1+1>=dims.size()){
    cout<<"Impossible number of the leg to be disentangled: "<<leg1<<" when dims="<<dims<<endl;
    exit(1);
  }
  if(Udis.getDimension(0)!=dims[leg1]*dims[leg1+1])
    Udis=identityMatrix(dims[leg1]*dims[leg1+1]); // start from identity? (or whatever came!)


}

