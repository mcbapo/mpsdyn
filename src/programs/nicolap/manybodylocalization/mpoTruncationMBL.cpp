#include <math.h>
#include <iomanip>
#include <fstream>
#include <deque>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "HeisenbergHamiltonian.h"
#include "Properties.h"
#include "mwArray.h"


#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0


using namespace shrt;
using namespace std;



string mpsfilename(int step,const string mpsdir);

int d=2;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  Properties props(infile);

  int M=props.getIntProperty("M");
  int D=props.getIntProperty("D");
  int truncD=props.getIntProperty("truncD");
  int steptruncation=props.getIntProperty("steptruncation");
  int stepMeasure=props.getIntProperty("stepMeasure");
  int rowInstances=props.getIntProperty("rowInstances");
  int importTime=props.getIntProperty("importTime");
  int lengthEvolution=props.getIntProperty("lengthEvolution");
  double h=props.getDoubleProperty("h");
  double tau=props.getDoubleProperty("tau");
  string outfname=props.getProperty("outfile");
  string infname=props.getProperty("instances");
  string mpsdir=props.getProperty("mpsdir");
  string outputMPS="";
  
  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(1E-8);
  cout<<"Initialized Contractor"<<endl;

  
  mwArray instances;
  ifstream in(infname.data());
  if(!in.is_open()){
    cout<<"Error: couldn't open file "<<infname<<" for output"<<endl;
    exit(1);
  }
  instances.load(in);
  in.close();

  ofstream* out;
  out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(15);

  ofstream fileInst("instFile.txt");
  instances.savetext(fileInst, true);
  fileInst.close(); 
   
  vector<double> randField;
  int c=0,newc=0;

    
  while(randField.size() < M ) {
    if(c%96==0 && c>0) {rowInstances+=1; newc+=96;}
    randField.insert(randField.begin(), h*real(instances.getElement(Indices(rowInstances,c-newc))));c+=1;
    randField.push_back(h*real(instances.getElement(Indices(rowInstances,c-newc))));c+=1;    
  }
  
  HeisenbergHamiltonian hamH(M,1.,1.,1.,randField,d);
  const MPO& mpoH=hamH.getHMPO();
  cout<<"Initialized basic Hamiltonian"<<endl;

  MPO doubleH(M);
  extendMPO(mpoH,doubleH,d);

  MPO mpoEven(M),mpoHalfEven(M), mpoOdd(M);
  {
    MPO auxU(M);
    hamH.getExponentialMPOeven(auxU,-tau*I_c);
    doubleMPO(auxU,mpoEven,false);
  }
  {
    MPO auxU(M);
    hamH.getExponentialMPOodd(auxU,-tau*I_c);
    doubleMPO(auxU,mpoOdd,false);
  }
  {
    MPO auxU(M);
    hamH.getExponentialMPOodd(auxU,-tau*I_c);
    doubleMPO(auxU,mpoHalfEven,false);
  }
  
  string inputMPS=mpsfilename(importTime,mpsdir);       
  MPS thS(M,1,d*d);
  if ( file_exists(inputMPS) ) {
    thS.importMPS(inputMPS.data());
  } else {
    thS.setProductState(p_maxent); // initial identity
  }

  double norm2;
  complex_t energy2,supInit;
  double err;

  const MPS initialState(thS);


  *out<<"#Truncation error for an MPO of initial bond dimension "<<D
      <<", length "<<M<<" and time step"<<tau<<endl;
  *out<<"#Time";

  int i = D;
  while (i >= truncD) {
    *out<<"\t bd"<<i;
    i -= steptruncation;
  }
  *out<<endl;

  {
  MPS aux(thS);
  contractor.optimize(mpoHalfEven,aux,thS,D);
  }

  for( int k=0 ; k<lengthEvolution ; k++ ) {
    
    
    MPS aux(thS);
    contractor.optimize(mpoOdd,aux,thS,D);    
    aux = thS;
    contractor.optimize(mpoEven,aux,thS,D);    

    thS.gaugeCond('R',true);
    norm2=real(contractor.contract(thS,thS)); 
    // energy2=contractor.contract(thS,doubleH,thS);
    // supInit=contractor.contract(thS,initialState);
    

    if (k % 100 == 0 ){
      outputMPS=mpsfilename(k,mpsdir);   
      thS.exportMPS(outputMPS.data());
    }

    if ( k % stepMeasure == 0 && k > 0) {
      i = D;
      *out<<2.*k*tau<<"\t"<<norm2;
      cout<<2.*k*tau<<"\t"<<norm2;
      MPS aux(thS);
      contractor.optimize(mpoHalfEven,aux,thS,D);
      
      while (i > truncD) {
	
	i -= steptruncation;
	MPS truncated(M,i,d*d);
	contractor.optimizeMPS(thS,truncated,i,&err);
	*out<<"\t"<<err;
	cout<<"\t"<<err;
      }
      *out<<endl;
      cout<<endl;
    }
  }
  out->close();delete out;
  
}



string mpsfilename(int step,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS"<<"_k"<<step<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

