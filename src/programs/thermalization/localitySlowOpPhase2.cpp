#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0


using namespace shrt;
using namespace std;

/** After slowlocal(localitySlowOp) has computed the weights of terms
 with different range (in norm 2) in the individual O[n] terms and the
 crossed terms, compute the weight in the full A_M operator, as
 defined by a given set of coefficients, c_m. */ 

int d=2;
//bool traceless=true;
bool traceless=false;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int M=atoi(argv[++cntr]); // max for which R is computed
  const char* infnameR=argv[++cntr]; // input: matrix R(r)
  const char* infnameG=argv[++cntr]; // input: coefficients g (filtered op) => recomputing the Cs
  const char* outfnameX=argv[++cntr]; // output: weights of various ranges as a function of m<=M 
  //  int rSeed=atoi(argv[++cntr]);

  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(1E-8);
  cout<<"Initialized Contractor"<<endl;
  ifstream inR(infnameR);
  if(!inR.is_open()){
    cout<<"Error: impossible to open file "<<infnameR<<
      " for input"<<endl;
    exit(1);
  }
  // Read the R matrix line by line
  // discard first line
  string s;
  getline(inR,s);
  mwArray R(Indices(4*(M+1),2*M+1,2*M+1));
  int r,n,m;double val;
  int Mr=1;
  while(inR){
    inR>>r>>n>>m>>val;
    R.setElement(val*ONE_c,Indices(r,M+n,M+m));
    R.setElement(val*ONE_c,Indices(r,M+m,M+n));
    if(inR&&abs(m)>Mr&&m<0&&r==1&&n==m){ 
      // It might be the last element of step m
      Mr=abs(m);
    }
  }
  inR.close();
  cout<<"Reconstructed R "<<R.getDimensions()<<" with highest seen m="<<Mr<<endl;
  // But maybe this one was not complete!
  // Check if the last to be written was actually saved, or discard the last one
  //Mr=Mr-1;


  // Now read the g's that allowed us to compute the optimal Cs
  ifstream in(infnameG);
  if(!in.is_open()){
    cout<<"Error: impossible to open file "<<infnameG<<
      " for input"<<endl;
    exit(1);
  }
  // discard first line
  getline(in,s);

  vector<double> g; // results

  int p;
  double reG,imG;
  int Mg=M;
  for(int k=0;k<=2*M+1;k++){
    in>>p;
    in>>reG>>imG;
    g.push_back(reG);
    if(!in){
      // Largest read is k => max M is floor((k-1)/2)
      Mg=floor((k-1)/2);
      cout<<"# Error happened while reading from "<<infnameG<<"!"<<endl;
      cout<<"# Gs only computed for M="<<Mg<<endl;
      // cout<<"# Rerun: "<<endl;
      // cout<<"./slowLocal2 "<<M<<" "<<infnameR<<" "<<infnameG<<" "<< outfnameX<<endl;
      //exit(1);
      break;
    }
  }
  in.close();

  int Meff=min(Mr,Mg);
  if(Meff<M){ // cut the R accordingly
    cout<<"Cutting the R to Mg="<<Meff<<endl;
    mwArray aux(R);
    R=mwArray(Indices(4*(Meff+1)+1,2*Meff+1,2*Meff+1));
    for(int r=1;r<=4*(Meff+1);r++)
      for(int m=-Meff;m<=Meff;m++)
	for(int n=-Meff;n<=Meff;n++)
	  R.setElement(aux.getElement(Indices(r,M+m,M+n)),Indices(r,Meff+m,Meff+n));
    M=Meff;
    //cout<<"R cut to "<<R<<endl;
  }
  if(Meff<Mg){
    cout<<"Cutting the G to Mr="<<Meff<<endl; // Not needed: I just won't use it
  }

  M=Meff;

  // Now I can construct the matrices for the optimization problem, for m=1,..M, compute the Cs 
  // and then the weights
  ofstream* out=new ofstream(outfnameX);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfnameX<<
      " for output"<<endl;
    exit(1);
  }
  *out<<"% M="<<M<<", inputfileR="<<infnameR<<", inputfileG="<<infnameG<<endl;
  *out<<setprecision(10);
  *out<<"% M\t G(r): 1, 2, 3..."<<endl;
  for(int m=1;m<=M;m++){
    cout<<"Computing for m="<<m<<endl;
    mwArray X(Indices(2*m+1,2*m+1)),Y(Indices(2*m+1,2*m+1));
    for(int i=0;i<2*m+1;i++){
      for(int j=0;j<2*m+1;j++){
	X.setElement(ONE_c*(g[abs(j-i)]-g[abs(j-i+1)]),Indices(i,j));
	if(!traceless)
	  Y.setElement(ONE_c*(g[abs(j-i)]-.5),Indices(i,j));
	else
	  Y.setElement(ONE_c*(g[abs(j-i)]),Indices(i,j));
      }      
    }
    X=X+Hconjugate(X);
    //    cout<<"X="<<X<<endl;
    // cout<<"Y="<<Y<<endl;
    //Could be done more efficiently by proper gEV solver, butI hope the denominator is not singular, so I will invert it
    mwArray Uy;vector<complex_t> Dy;
    wrapper::eig(Y,Dy,Uy,true);
    // Now the generalizes eigenvalues, assuming Y is invertible, are the eigenvalues of Dy^-1/2 Uy^+ X Uy Dy^-1/2
    int nr=0;double tol=1E-10;
    mwArray D2=invertDiag(sqrt(diag(Dy)),nr,tol);
    D2.multiplyLeft(Uy);
    //    exit(1);
    X.multiplyRight(D2);
    X.multiplyLeft(Hconjugate(D2));
    mwArray Ux;vector<complex_t> Dx;
    wrapper::eig(X,Dx,Ux,true);
    // cout<<"Eigenvalues for m="<<m<<" Dx="<<Dx<<endl;
    //*out<<m<<"\t"<<Dx[0]<<endl;

    // Now I want to compute the Cs
    mwArray Ux0=Ux.subArray(Indices(-1,0)); // first column
    Ux0.multiplyLeft(D2);
    // I need to fill in with zeros up to the 2*M+1 dimension of R
    mwArray Cs(Indices(R.getDimension(1),1));
    for(int p=0;p<Ux0.getDimension(0);p++){
      Cs.setElement(Ux0.getElement(Indices(p,0)),Indices(p+M-m,0));
    }

    mwArray Ch(Cs);Ch.Hconjugate(); // Hermitian conjugate

    // To normalize, tr(A+A)
    mwArray trAA=Hconjugate(Ux0)*Y*Ux0;
    if(!trAA.isScalar()){
      cout<<"Error: didn't get a scalar when multiplying coefficients "<<Ux0<<" denominator="<<Y<<endl;
      exit(1);
    }
    double trA2=real(trAA.getElement(0));
    // and now, for each range, I need the contraction C* R(r) C
    *out<<m<<"\t"<<trA2<<"\t";
    for(int r=1;r<=4*(M+1);r++){
      if(r<=4*m+1){
	mwArray Rr=R.subArray(Indices(r,-1,-1));
	Rr.multiplyRight(Cs);
	Rr.multiplyLeft(Ch);
	if(!Rr.isScalar()){
	  cout<<"Error: didn't get a scalar when multiplying coefficients "<<Ux0<<" with Rr (m="<<m<<",r="<<r<<") "<<Rr<<endl;
	  exit(1);
	}
	*out<<real(Rr.getElement(0))<<"\t";
      }
      else
	*out<<0<<"\t";
    }
    *out<<endl;

  }
  out->close();
  delete out;
}
