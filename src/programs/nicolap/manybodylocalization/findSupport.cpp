#include <math.h>
#include <iomanip>
#include <fstream>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"
#include "mwArray.h"

#include "HeisenbergHamiltonian.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

using namespace shrt;
using namespace std;

/* This program compute the overlap of an input MPS with combinations of Pauli Matrices of different sizes */

const string mpsfilename(int M,int D,const string mpsdir);

#define MAXITER 100
#define MAXPERT 10

//#define ISINGCASE 1

int d=2;
double gP=0.853;
double tau=0.8;
bool traceless=true;

int main(int argc,const char* argv[]){

  int cntr=0;
  const char* infile=argv[++cntr];
 
  Properties props(infile);
  


  int M=props.getIntProperty("M");
  int D=props.getIntProperty("D");
  int corr=props.getIntProperty("corr");
  string mpsdir= props.getProperty("mpsdir");  
  string outfi=props.getProperty("outfile");
  string outavefi=props.getProperty("avefile");
  string outEntfi=props.getProperty("entfile");
  string outSxfi=props.getProperty("sx");
  string outSyfi=props.getProperty("sy");
  string outSzfi=props.getProperty("sz");
  string outSchmidtfi=props.getProperty("Schmidt");

  string binD=props.getProperty("binD");
  
  string outfile = binD+outfi;
  string outavefile = binD+outavefi;
  string outEntfile = binD+outEntfi;
  string outSxfile = binD+outSxfi;
  string outSyfile = binD+outSyfi;
  string outSzfile = binD+outSzfi;  
  string outSchmidtfile = binD+outSchmidtfi;  


  cout<<"Append output in: "<<outfile<<endl<<endl;
  cout<<"Append average output in: "<<outavefile<<endl;
  string inputMPS=mpsfilename(M,D,mpsdir);
  const string tmpinfile = inputMPS+".tmp";
  

  MPS Omps(M,D,d*d);
  if(file_exists(inputMPS)){
    
    if(inputMPS.empty()){ // no init file indicated in arguments
      cout<<"There is no input file for initD "<<endl;
    }else{ // initMPS file read
      Omps.importMPS(inputMPS.data());
      cout<<"Imported initMPS file for initial state"<<inputMPS<<endl;
    }
  }else { // no tmp file
    if(file_exists(tmpinfile)){
      cout<<"There is only tmp file for same D: "<<tmpinfile<<endl;
      Omps.importMPS(tmpinfile.data());
      cout<<"Imported tmp file for same D!"<<tmpinfile<<endl;
    } else {
      cout<<"File to import not found"<<endl;
      exit(1);
    }
  }

  Contractor& contractor=Contractor::theContractor(); //inizialize contractor

  mwArray id2=1/sqrt(2.)*identityMatrix(d);id2.reshape(Indices(4,1));
  mwArray id2T=1/sqrt(2.)*identityMatrix(d);id2T.reshape(Indices(1,4));
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),datax);sigX.reshape(Indices(4,1));
  mwArray sigXT(Indices(d,d),datax);sigXT.reshape(Indices(1,4));
  complex_t datay[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),datay);sigY.reshape(Indices(4,1));
  mwArray sigYT(Indices(d,d),datay);sigYT.reshape(Indices(1,4));
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);sigZ.reshape(Indices(4,1));
  mwArray sigZT(Indices(d,d),dataz);sigZT.reshape(Indices(1,4));
  

  mwArray ID4=identityMatrix(4);
  mwArray id4 = id2*id2T;
  mwArray Sx = sigX*sigXT;
  mwArray Sy = sigY*sigYT;
  mwArray Sz = sigZ*sigZT;
  mwArray sumPauli = Sx+Sy+Sz;

  Operator id4Op(id4), sPaulOp(sumPauli), ID4Op(ID4);
  Operator sigmaX(Sx), sigmaY(Sy), sigmaZ(Sz);

  vector<const mwArray*> tensors;
  for(int k=0; k<M; k++)
    tensors.push_back(&ID4);

  MPO findSupport(M);
  cout<<"Empty operator created"<<endl;

  int lengthSup,k,i;
  double sumPaulis, sum,sx,sy,sz;
  
  // double norm;
  // norm=real(contractor.contract(Omps,Omps));
  // cout<<"The norm is equal to: "<<norm<<endl;
  ofstream* out;  
  out=new ofstream(outfile.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<"for output"<<endl;
    exit(1);
  }
  *out<<setprecision(15);
  
  ofstream* outAve;
  outAve=new ofstream(outavefile.data(),ios::app);
  if(!outAve->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<"for output"<<endl;
    exit(1);
  }
  *outAve<<setprecision(15);

  ofstream* outEnt;
  outEnt=new ofstream(outEntfile.data(),ios::app);
  if(!outEnt->is_open()){
    cout<<"Error: impossible to open file "<<outEnt<<"for output"<<endl;
    exit(1);
  }
  *outEnt<<setprecision(15);
  
  ofstream* outSx;
  outSx=new ofstream(outSxfile.data(),ios::app);
  if(!outSx->is_open()){
    cout<<"Error: impossible to open file "<<outSx<<"for output"<<endl;
    exit(1);
  }
  *outSx<<setprecision(15);

  ofstream* outSy;
  outSy=new ofstream(outSyfile.data(),ios::app);
  if(!outSy->is_open()){
    cout<<"Error: impossible to open file "<<outSy<<"for output"<<endl;
    exit(1);
  }
  *outSy<<setprecision(15);

  ofstream* outSz;
  outSz=new ofstream(outSzfile.data(),ios::app);
  if(!outSz->is_open()){
    cout<<"Error: impossible to open file "<<outSz<<"for output"<<endl;
    exit(1);
  }
  *outSz<<setprecision(15);

  ofstream* outSchmidt;
  outSchmidt=new ofstream(outSchmidtfile.data(),ios::app);
  if(!outSchmidt->is_open()){
    cout<<"Error: impossible to open file "<<outSchmidt<<"for output"<<endl;
    exit(1);
  }
  *outSchmidt<<setprecision(15);



  const Operator *Pid4 = &id4Op;
  const Operator *PsumPauli = &sPaulOp;
  const Operator *SIGMAz = &sigmaZ;
  const Operator *SIGMAx = &sigmaX;
  const Operator *SIGMAy = &sigmaY;
  const Operator *PID4 = &ID4Op;

  vector<complex_t> SchmidtCoef(D);

  for(lengthSup=1; lengthSup<=M; lengthSup++){
    *out<<lengthSup<<"\t";
    *outSx<<lengthSup<<"\t";
    *outSy<<lengthSup<<"\t";
    *outSz<<lengthSup<<"\t";
    *outAve<<lengthSup<<"\t";

    sum=0;
    for(k=0; k<M-lengthSup+1; k++){
      findSupport.setOperatorArray(tensors);

      //Check sum of Paulis 
      findSupport.setOp(k,PsumPauli,false);
      if ( lengthSup > 1)
	findSupport.setOp(k+lengthSup-1,PsumPauli,false);
      sumPaulis=real(contractor.contract(Omps,findSupport,Omps));
      //Check Sx
      findSupport.setOp(k,SIGMAx,false);
      if ( lengthSup > 1)
	findSupport.setOp(k+lengthSup-1,SIGMAx,false);
      sx=real(contractor.contract(Omps,findSupport,Omps));
      //Check Sy
      findSupport.setOp(k,SIGMAy,false);
      if ( lengthSup > 1)
	findSupport.setOp(k+lengthSup-1,SIGMAy,false);
      sy=real(contractor.contract(Omps,findSupport,Omps));
      //Check Sz
      findSupport.setOp(k,SIGMAz,false);
      if ( lengthSup > 1)
	findSupport.setOp(k+lengthSup-1,SIGMAz,false);
      sz=real(contractor.contract(Omps,findSupport,Omps));

      if ( corr ) {
	for(i=k+1; i<k+lengthSup-1; i++)
	  findSupport.setOp(i,Pid4,false);
      }


      *out<<sumPaulis<<"  \t";
      *outSx<<sx<<"  \t";
      *outSy<<sy<<"  \t";
      *outSz<<sz<<"  \t";
      
      sum += sumPaulis;
      
    }
    *out<<endl;
    cout<<endl;
    *outAve<<sum/(M-lengthSup)<<endl;
    if (lengthSup > 0 && lengthSup < M){
      *outEnt<<lengthSup<<"\t";
      *outEnt<<contractor.getEntropy(Omps,lengthSup)<<endl;
      
      *outSchmidt<<lengthSup<<"\t";
      contractor.getSchmidtValues(Omps,SchmidtCoef, lengthSup);
      for (int i=0; i<SchmidtCoef.size();i++)
	*outSchmidt<<real(SchmidtCoef[i])<<"\t";
      *outSchmidt << endl;
	
	
    }


    
  }

  

  

  out->close();
  delete out;
  outAve->close();
  delete outAve;
  outEnt->close();
  delete outEnt;
  outSx->close();
  delete outSx;
  outSy->close();
  delete outSy;
  outSz->close();
  delete outSz;
  outSchmidt->close();
  delete outSchmidt;

}

const string mpsfilename(int M,int D,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS"<<"_M"<<M<<"_D"<<D<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}
