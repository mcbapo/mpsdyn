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
    slow relaxing operator, afte the Gs have been computed by slowop,
    recover the file and compute the lambdas again. Also compute the cn's */

int d=2;
//bool traceless=true;
//bool traceless=false;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int M=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* infnameG=argv[++cntr];
  const char* outfnameL=argv[++cntr];
  const char* outfnameC=argv[++cntr];
  double contribE=1.;
  double contribO=1.;
  if(argc>6){
    contribE=atof(argv[++cntr]);
    contribO=atof(argv[++cntr]);
    cout<<"Using weights for even and odd contributions "<<contribE<<"/"<<contribO<<endl;
  }
  bool correct0=0;
  if(argc>8){
    cout<<"Special case, correcting for error in the first entry!!!"<<endl;
    correct0=true;
  }
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

  vector<double>  gEE,gEO,gOE,gOO; // results

  int p;
  double reGee,imGee;
  double reGeo,imGeo;
  double reGoe,imGoe;
  double reGoo,imGoo;
  for(int k=0;k<=2*M+1;k++){
    in>>p;
    in>>reGee>>imGee;
    in>>reGeo>>imGeo;
    in>>reGoe>>imGoe;
    in>>reGoo>>imGoo;
    if(correct0&&k==0){
      gEE.push_back(reGee-.5);
      gEO.push_back(reGeo-.5);
      gOE.push_back(reGoe-.5);
      gOO.push_back(reGoo-.5);
    }
    else{    
      gEE.push_back(reGee);
      gEO.push_back(reGeo);
      gOE.push_back(reGoe);
      gOO.push_back(reGoo);
    }
    if(!in){
      // Largest read is k => max M is floor((k-1)/2)
      int Mold=M;
      M=floor((k-1)/2);
      cout<<"# Error happened while reading from "<<infnameG<<"!"<<endl;
      cout<<"# Only computing for M="<<M<<endl;
      cout<<"# Rerun: "<<endl;
      cout<<"./slowTI2 "<<Mold<<" "<<D<<" "<<infnameG<<" "<< outfnameL<<" "<< outfnameC<<" "<<contribE<<" "<<contribO;
      if(correct0) cout<<" 1";
      cout<<endl;
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
  *out<<"% contribE="<<contribE<<", contribO="<<contribO<<endl;
  *out<<"% M\t min(lambda)"<<endl;
  for(int m=1;m<=M;m++){
    mwArray X(Indices(2*m+1,2*m+1)),Y(Indices(2*m+1,2*m+1));
    for(int i=0;i<2*m+1;i++){
      for(int j=0;j<2*m+1;j++){
	X.setElement(ONE_c*(contribE*(contribE*(gEE[abs(j-i)]-gEE[abs(j-i+1)])+contribO*(gEO[abs(j-i)]-gEO[abs(j-i+1)]))+
			    contribO*(contribE*(gOE[abs(j-i)]-gOE[abs(j-i+1)])+contribO*(gOO[abs(j-i)]-gOO[abs(j-i+1)]))),Indices(i,j));
// 	if(!traceless)
// 	  Y.setElement(ONE_c*(gEE[abs(j-i)]-.5+gEO[abs(j-i)]-.5+
// 			      gOE[abs(j-i)]-.5+gOO[abs(j-i)]-.5),Indices(i,j));
// 	else
	Y.setElement(ONE_c*(contribE*(contribE*gEE[abs(j-i)]+contribO*gEO[abs(j-i)])+
			    contribO*(contribE*gOE[abs(j-i)]+contribO*gOO[abs(j-i)])),Indices(i,j));
       }      
    }
    X=X+Hconjugate(X);
    //    cout<<"X="<<X<<endl;
    // cout<<"Y="<<Y<<endl;
    //Could be done more efficiently by proper gEV solver, butI hope the denominator is not singular, so I will invert it
    mwArray Uy;vector<complex_t> Dy;
    wrapper::eig(Y,Dy,Uy,true);
    // Now the generalizes eigenvalues, assuming Y is invertible, are the eigenvalues of Dy^-1/2 Uy^+ X Uy Dy^-1/2
    // I need now to do sth different if Y is singular!
    int nr=0;double tol=1E-8;
    vector<int> nonZeroEVy;
    complex_t largestEV=Dy.back();
    for(int k=0;k<Dy.size();k++){
      if(abs(Dy[k]/largestEV)>=tol){
	nonZeroEVy.push_back(k);
	// I'll also check for negative or non-real alues!
	if(abs(imag(Dy[k])/real(Dy[k]))>=tol){
	  cout<<"WARNING: Step "<<m<<" non-real eigenvalue("<<k<<") of Y "<<Dy[k]<<endl;
	}
	if(real(Dy[k])<0){
	  cout<<"WARNING: Step "<<m<<" negative eigenvalue("<<k<<") of Y "<<Dy[k]<<endl;
	  ofstream* outLog=new ofstream("logRecomputeLambdas.m",ios::app);
	  *outLog<<"# Error in step "<<m<<" from ="<<infnameG<<endl;
	  *outLog<<"# Negative eigenvalue("<<k<<") of Y "<<Dy[k]<<endl;
	  putForMatlab(*outLog,Y,"Y");
	  putForMatlab(*outLog,X,"X");
	  *outLog<<"### #"<<endl;
	  outLog->close();
	  delete outLog;
	  exit(1);
	}

      }
    }
    // if(m==1){
    //   putForMatlab(cout,X,"X");
    //   cout<<"Eigenvalues of the denominator for m="<<m<<" are "<<Dy<<endl;
    // }
    mwArray D2(Indices(nonZeroEVy.size(),nonZeroEVy.size())),Uyeff(Indices(Dy.size(),nonZeroEVy.size()));
    if(nonZeroEVy.size()<Dy.size()){
      cout<<"WARNING: singular denominator. Only "<<nonZeroEVy.size()<<" (out of "<<Dy.size()<<") non-zero values found at "<<nonZeroEVy<<" with values ";
      for(int p=0;p<nonZeroEVy.size();p++)
	cout<<Dy[nonZeroEVy[p]]<<",";
      cout<<endl;
      // Construct the isometry removing columns with 0
      for(int p=0;p<nonZeroEVy.size();p++){
	//cout<<"Taking column "<<nonZeroEVy[p]<<" of Uy "<<endl;
	D2.setElement(ONE_c/sqrt(Dy[nonZeroEVy[p]]),Indices(p,p));      
	for(int q=0;q<Dy.size();q++)
	  Uyeff.setElement(Uy.getElement(Indices(q,nonZeroEVy[p])),Indices(q,p));      
      }
    }
    else{
      D2=invertDiag(sqrt(diag(Dy)),nr,sqrt(tol));
      Uyeff=Uy;
    }
    //    D2.multiplyLeft(Uy);
    D2.multiplyLeft(Uyeff);
    //    exit(1);
    X.multiplyRight(D2);
    X.multiplyLeft(Hconjugate(D2));
    mwArray Ux;vector<complex_t> Dx;
    wrapper::eig(X,Dx,Ux,true);
    //cout<<"Eigenvalues for m="<<m<<" Dx="<<Dx<<", Ux="<<Ux<<endl;
    *out<<m<<"\t"<<Dx[0]<<endl;
    if(real(Dx[0])<0&&abs(Dx[0])>=tol){
      cout<<"ERROR! Negative minimum value, probably from numerical errors!"<<endl;
      exit(1);
    }

    // Now I want to compute the Cs and write them to a file
    mwArray Ux0;
    if(Ux.isScalar()) Ux0=Ux; // special case, when only one
    else
      Ux0=Ux.subArray(Indices(-1,0)); // first column
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
