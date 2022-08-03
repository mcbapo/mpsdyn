#include "mwArray.h"
#include "misc.h"
#include "Contractor.h"
#include "Properties.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>

using namespace std;
using namespace shrt;

/** Read a MPS (binary) from a file and export it to text or matlab format */


const string mpsFileName(const string& baseDir,const string& name,const string& ext="");

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  string directory="";
  // if(argc>2){
  //   const char* indir=argv[++cntr];
  //   directory=string(indir)+"/";
  // }
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  const string mpsdir = props.getProperty("inputMPSdir");
  const string mpsdirOut = props.getProperty("outputMPSdir");
  const string mpsname = props.getProperty("inputFile");
  const string mpsnameOut = props.getProperty("outputFile");
  bool matlabFormat=props.getIntProperty("matlab")>0;
  bool txtFormat=props.getIntProperty("txt")>0;

  if(!matlabFormat&&!txtFormat){
    cout<<"ERROR: No output format selected!"<<endl;
    exit(1);
  }
  
  string mpsfile=mpsFileName(mpsdir,mpsname);
  cout<<"Will read MPs from "<<mpsfile<<endl;

  MPS state(1,1,1);
  if(file_exists(mpsfile))
    state.importMPS(mpsfile.data());
  else{
    cout<<"ERROR: Cannot find input MPS in file "<<mpsfile<<endl;
    exit(1);
  }

  cout<<"Read MPS "<<state<<"from file"<<endl;
  int L=state.getLength();
  int d=state.getA(0).getd();
  int D=state.getBond();
  // for(int k=0;k<L;k++){
  //   cout<<"Showing tensor nr "<<k+1<<endl;
  //   cout<<state.getA(k).getA()<<endl;
  // }

  if(txtFormat){
    string mpsfileOut=mpsFileName(mpsdirOut,mpsnameOut,".txt");
    state.exportMPStext(mpsfileOut.data());
    cout<<"Exported to text format in file "<<mpsfileOut<<endl;
  }
  if(matlabFormat){
    string mpsfileOut=mpsFileName(mpsdirOut,mpsnameOut,".m");
    state.exportForMatlab(mpsfileOut.data());
    cout<<"Exported to matlab format in file "<<mpsfileOut<<endl;
  }

}

const string mpsFileName(const string& baseDir,const string& name,const string& ext){
  stringstream s;
  s<<baseDir<<"/"<<name<<ext;
  //s<<".txt";
  return s.str();
}

