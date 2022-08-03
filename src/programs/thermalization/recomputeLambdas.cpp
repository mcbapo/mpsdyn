#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0


using namespace shrt;
using namespace std;

/** Compute numerically the terms appearing in the coefficient of the
    slow relaxing operator, after the Gs have been computed by slowop,
    recover the file and compute the lambdas again. Also compute the cn's */

int d=2;
//bool traceless=true;
bool traceless=false;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int M=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* infnameG=argv[++cntr];
  const char* outfnameL=argv[++cntr];
  const char* outfnameC=argv[++cntr];
  //  int rSeed=atoi(argv[++cntr]);

  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(1E-8);
  cout<<"Initialized Contractor"<<endl;
  ifstream in(infnameG);
  if(!in.is_open()){
    cout<<"Error: impossible to open file "<<infnameG<<
      " for input"<<endl;
    exit(1);
  }
  // discard first line
  string s;
  getline(in,s);

  vector<double> g; // results

  int p;
  double reG,imG;
  for(int k=0;k<=2*M+1;k++){
    in>>p;
    in>>reG>>imG;
    g.push_back(reG);
    if(!in){
      // Largest read is k => max M is floor((k-1)/2)
      M=floor((k-1)/2);
      cout<<"# Error happened while reading from "<<infnameG<<"!"<<endl;
      cout<<"# Only computing for M="<<M<<endl;
      cout<<"# Rerun: "<<endl;
      cout<<"./slow2 "<<M<<" "<<D<<" "<<infnameG<<" "<< outfnameL<<" "<< outfnameC<<endl;
      //exit(1);
      break;
    }
  }

  in.close();


  // Now I can construct the matrices for the optimization problem, for m=1,..M
  ofstream* out=new ofstream(outfnameL);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfnameL<<
      " for output"<<endl;
    exit(1);
  }
  *out<<"% D="<<D<<", M="<<M<<", inputfile="<<infnameG<<endl;
  *out<<setprecision(10);
  *out<<"% M\t min(lambda)"<<endl;
  for(int m=1;m<=M;m++){
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
    //    cout<<"Eigenvalues for m="<<m<<" Dx="<<Dx<<endl;
    *out<<m<<"\t"<<Dx[0]<<endl;

    // Now I want to compute the Cs and write them to a file
    mwArray Ux0=Ux.subArray(Indices(-1,0)); // first column
    Ux0.multiplyLeft(D2);
    stringstream s;
    string nameC(outfnameC);
    s<<nameC<<"_M"<<m;
    ofstream* out2=new ofstream(s.str().data());
    if(!out2->is_open()){
      cout<<"Error: impossible to open file "<<outfnameC<<"_M"<<m<<
	" for output"<<endl;
      exit(1);
    }
    for(int k=0;k<Ux0.getDimension(0);k++)
      *out2<<real(Ux0.getElement(Indices(k,0)))<<endl;
    out2->close();
    delete out2;

  }
  out->close();
  delete out;
}
