
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "Properties.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "HeisenbergHamiltonian.h"

using namespace shrt;

/** Auxiliary, to compute input required by effectiveH
*/


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

  const string mpsfile=props.getProperty("mpsfile"); // full path expected
  const string initmpsfile=props.getProperty("initmpsfile"); // full path expected
  // Hamiltonian parameters
  int L=props.getIntProperty("L");
  int D=props.getIntProperty("D");
  double Jx_=props.getDoubleProperty("Jx");
  double Jy_=props.getDoubleProperty("Jy");
  double Jz_=props.getDoubleProperty("Jz");
  double h_=props.getDoubleProperty("h");  
  const string outfname=props.getProperty("outfile");  
  bool app=(props.getIntProperty("append")!=0); // default: append
  
  cout<<"Initialized arguments: L="<<L
      <<", Jx="<<Jx_
      <<", Jy="<<Jy_
      <<", Jz="<<Jz_
      <<", h="<<h_
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  vector<double> Jx(L-1,Jx_);
  vector<double> Jy(L-1,Jy_);
  vector<double> Jz(L-1,Jz_);
  vector<double> h(L,h_);
  
  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% N\t Jx\t Jy\t Jz\t h\t D\t Energy\t Energy/N"<<endl;
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);
  
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d=2;
  MPS gs(L,D,d);
  if(file_exists(initmpsfile)){
    gs.importMPS(initmpsfile.data());
    cout<<"REad initial state from file"<<endl;
  }
  else{
    gs.setRandomState(); // the intial state, random
    gs.setProductState(p_zero);
    gs.gaugeCond('R',1);
    gs.gaugeCond('L',1);
    cout<<"Initialized random state, norm "<<contractor.contract(gs,gs)<<endl;
  }
  
  // First: put H and find the GS
  int Dop=3; // bond dimension of Ising H
  double E0=0.;
  HeisenbergHamiltonian hamH(L,Jx,Jy,Jz,h,d);
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamH.getHMPO();
  cout<<"Constructed the hamil MPO"<<endl;
  //hamil.exportForMatlab("isingHmpo.m");
  cout<<"Initial value, with initial state"<<endl;
  cout<<contractor.contract(gs,hamil,gs)<<endl;

  contractor.findGroundState(hamil,D,&E0,gs);
  // And save to file
  gs.exportMPS(mpsfile.data());  

  cout<<"Ground state found with eigenvalue "<<E0
      <<" (Energy per particle="<<E0/L<<")"<<endl;

  *out<<L<<"\t"<<Jx_<<"\t"<<Jy_<<"\t"<<Jz_<<"\t"<<h_<<"\t"
      <<D<<"\t"<<E0<<"\t"<<E0/L<<"\t";
  *out<<endl;
  out->close();delete out;

  
}
