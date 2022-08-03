
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "mwArray.h"
#include "Properties.h"
#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

/** Runs the findGroundState routine with the MPO for the
    Ising Hamiltonian \ref <IsingHamiltonian> squared, to check 
    states with min variance for fixed bond dimension.

    It stops after a certain convergence criterion or as soon as a
    given threshold value of the variance is attained.  If a continue
    argument is given, it searches for the largest D for which a MPS
    file is available, and starts from that one, and not from D1.

*/

int findTmpFile(int D1,int Dmax,int incrD,const string& mpsdir,const string& baseName,string& mpsfile);
const string mpsFileName(const string& mpsdir,const string& baseName,int D);

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  srandom(time(NULL));
  const char* infile=argv[++cntr];
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  int L = props.getIntProperty("L");  
  double J = props.getDoubleProperty("J");  
  double g = props.getDoubleProperty("g");  
  double h = props.getDoubleProperty("h");  
  int D1 = props.getIntProperty("D1");  
  int D2 = props.getIntProperty("D2");
  int incrD = props.getIntProperty("incrD");  
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir"); // to save final MPS
  string mpsfile=props.getProperty("mpsfile"); // basename to save final MPS
  bool saveMPS=(mpsfile.size()>0);
  string initmpsfile=props.getProperty("initmpsfile"); // use for initial state
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  bool continued=(props.getIntProperty("continue")>0);
  if(continued) app=true; // overrides the append option
  double varTh = props.getDoubleProperty("targetD");   // threshold
  double tolC = props.getDoubleProperty("convTol");  // if threshold not met, tolerance to stop iteration
  if(tolC<0) tolC=1E-2; // set quite high
  
  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% N\t J\t g\t h\t D\t dt\t deltaX\t cnt\t time \t<H2>\t <H>"<<endl;
    out->close();
    delete out;
  }
  cout<<setprecision(10);

  
  int D=D1; // to start
  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  cout<<"Initialized Contractor"<<endl;
  int d=2;
  MPS gs(L,D,d);
  // First, if continue, search for an initial state
  if(continued){
    string tmpMPSfile;
    int Dtmp=findTmpFile(D1,D2,incrD,mpsdir,mpsfile,tmpMPSfile);
    if(Dtmp>D1){
      cout<<"Found intermediate result with D="<<Dtmp<<endl;
      D=Dtmp;
      gs.importMPS(tmpMPSfile.data());
    }
    else{
      cout<<"Cannot continue previous results, because no tmp file was found"<<endl;
      continued=false; // deactivate this option
    }
  }
  if(!continued){   // happens both if no continued option was given or no file was found
    if(!initmpsfile.empty()){
      gs.importMPS(initmpsfile.data());
      cout<<"Read initial state from "<<initmpsfile<<endl;
    }
    else{
      gs.setRandomState(); // the intial state, random
      cout<<"Initialized random state, norm "<<contractor.contract(gs,gs)<<endl;
    }
  }  

  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  
  // First: put H and find the GS
  IsingHamiltonian hamI(L,d,J,g,h); 
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamI.getHMPO(0.);
  complex_t valH=contractor.contract(gs,hamil,gs);
  double valH2=real(contractor.contract2(hamil,gs));
  double S=contractor.getEntropy(gs);
  cout<<"Constructed the hamil MPO. ";
  //hamil.exportForMatlab("isingHmpo.m");
  cout<<"Initial value of E, with initial state ";
  cout<<"<H>="<<valH;
  cout<<", and <H^2>="<<valH2<<endl;

  bool done=0;

  
  if(valH2<varTh){
    cout<<"!!!!! Initial state is already below the threshold, so I am writing current"
	<<" values to the file (but a smaller D may also get it!"<<endl;
    // I also compute the entropy    
    out=new ofstream(outfname.data(),ios::app);*out<<setprecision(10);
    *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"
	<<varTh<<"\t"<<D<<"\t"<<valH2<<"\t"<<real(valH)
	<<"\t"<<S<<endl;
    out->close();
    delete out;
    done=true;
  }

  
  const MPO* ptrs[]={&hamil,&hamil};
  MPO ham2(L);
  MPO::join(2,ptrs,ham2);// MPO for H^2
  //  MPO ham2(2,ptrs); 

  contractor.setConvTol(tolC);
  while(!done){ // find GS of H^2 until D is large enough
    //    contractor.findGroundState(ham2,D,&expH2,gs);
    contractor.findLowestEnergyThreshold(ham2,D,valH2,gs,varTh);
    // compute values
    valH=contractor.contract(gs,hamil,gs);
    S=contractor.getEntropy(gs);
    cout<<"State for D="<<D<<" found with variance "<<valH2
	<<" energy="<<valH<<", S(L/2)="<<S<<endl;

    if(valH2<0){
      cout<<"ERROR!!!!!! Variance cannot be negative! Value held: "<<valH2<<
	" and computed "<<contractor.contract2(hamil,gs)<<endl;
      exit(1);
    }

    out=new ofstream(outfname.data(),ios::app);*out<<setprecision(10);
    *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"
	<<varTh<<"\t"<<D<<"\t"<<valH2<<"\t"<<valH
	<<"\t"<<S<<endl;
    out->close();
    delete out;
    
    if(saveMPS){
      string newMPSfile=mpsFileName(mpsdir,mpsfile,D);
      gs.exportMPS(newMPSfile.data());
    }
    
    if(valH2<varTh)
      done=true;
    else{
      D+=incrD;
      if(D>D2) done=true;
      else gs.increaseBondDimension(D);
    }
  }
  
}

const string mpsFileName(const string& mpsdir,const string& baseName,int cnt){
  stringstream s;
  s<<mpsdir<<"/"<<baseName;
  if(cnt>0) s<<"_"<<cnt;
  return s.str();
}

int findTmpFile(int D1,int Dmax,int incrD,const string& mpsdir,const string& baseName,string& mpsfile){
  int D=Dmax;
  bool found=0;
  while(D>=D1&&!found){
    mpsfile=mpsFileName(mpsdir,baseName,D);
    if(file_exists(mpsfile)){
      cout<<"findTmpFile found mpsfile="<<mpsfile<<endl;
      found=true;
      return D;
    }
    else{
      mpsfile="";
      D-=incrD;
    }
  }
  if(!found) return 0;
  cout<<"ERROR: Should not be here"<<endl;
  exit(1);
}
