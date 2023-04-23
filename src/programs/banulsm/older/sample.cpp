#include <math.h>
#include <iomanip>
#include <typeinfo>

#include "PEPS.h"
#include "Contractor.h"
#include "Sampler.h"
#include <vector>

#define MAXSMPL 1000

using namespace std;
using namespace shrt;

/** sample computes the magnetization in a position (x,y) of a lattice
    in the state given by a PEPS Receives arguments: 
    @param inputfile where the PEPS is stored as txt 
    @param x (0 -- Nrows-1) horizontal position of the operator
    @param y (0 -- Ncolumns-1) vertical position of the operator; if x
    or y are negative, the average magnetization is computed, instead
    of the local one;
    @param mode (0=mps/1=site) to decide which sampling to run 
    @param Meq nr of configurations before considering a valid one 
    @param Mtot nr stored for averaging
    @param Mb every Mb configs the average is saved 
    @param Dprime bond dimension used to approximate the environment
    @param oufile where to write the results
 */

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* inputfname=argv[++cntr];
  int x=atoi(argv[++cntr]);
  int y=atoi(argv[++cntr]);
  int mode=atoi(argv[++cntr]);
  int Meq=atoi(argv[++cntr]);
  int Mtot=atoi(argv[++cntr]);
  int Mb=atoi(argv[++cntr]);
  int Dprime=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  
  srandom(time(NULL));

  PEPS C;
  C.importPEPStext(inputfname);    //"../input/PEPSL3D2B10DExp50.txt"
  int Nr=C.getNrRows();
  int Nc=C.getNrColumns();
  cout<<"Read a PEPS "<<Nr<<"x"<<Nc<<"from file "<<inputfname<<endl;

  bool average=false; 
  Indices pos(x,y);
  if( (x>=Nr)||
      (y>=Nc) ){
    cout<<"Error, position "<<pos<<" exceeds the boundaries of this lattice!!"<<endl;
    exit(1);
  }
  if((x<0)||(y<0)){
    cout<<"Computing average operator"<<endl;
    average=true;
  }
     
  
  Sampler& sampl=Sampler::theSampler();
  
  int d=C.getd();
  if(d!=2){
    cout<<"Error, I cannot compute magnetizatioon in a PEPS of dimension "<<d<<endl;
    exit(1);
  }
  mwArray oper=identityMatrix(d); //  identity
  oper.setElement(-1*ONE_c,Indices(1,1)); // sigma_z
  
  ofstream* out=new ofstream(outfname);
  if(!out->is_open()){
    cout<<"Error: could not open file to write!"<<endl;
    exit(1);
  }
  *out<<setprecision(15);
  double avrg;
  if(!average)
  // For a single site magnetization
    avrg=sampl.sampleLocalOperator(C,oper,pos,Dprime,*out,SampleMode(mode),Meq,Mtot,Mb);
  else
  // For the average magnetization
    avrg=sampl.sampleAverageLocalOperator(C,oper,Dprime,*out,SampleMode(mode),Meq,Mtot,Mb);
  *out<<setprecision(15)<<"# Average="<<avrg<<endl;
  out->close();
  cout<<"Obtained "<<setprecision(15)<<avrg<<" with Mtot="<<Mtot<<endl;
  delete out;
}
