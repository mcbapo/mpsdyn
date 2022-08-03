
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"


#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

/** Runs the findGroundState routine with the MPO for the
    Ising Hamiltonian \ref <IsingHamiltonian> squared, to check 
    states with min variance for fixed bond dimension.

In this version, find it with the largest D first, and then go cutting down!
Save the largest D result

    Receives arguments:
    \param <L> (int) length of the chain
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D1> (int) maximum bond dimension
    \param <D2> (int) maximum bond dimension
    \param <deltaD> (int) step change in bond dim
    \param <outfname> (char*) name of the output file for the energy
    \param <MPSname> (char*) name of the output file for the MPS
    \param <app> (int) whether the output file is to be kept (app==1)
*/

void computeValues(const MPS& state,const MPO& hamil,vector<complex_t>& valH2,vector<complex_t>& valH);

int d=2;

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
  string inMPS=props.getProperty("inputMPS");
  string outMPS=props.getProperty("outputMPS");
  bool readMPS=!(inMPS.empty());
  bool saveMPS=!(outMPS.empty());
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  bool optim=(props.getIntProperty("optimize")>0);

  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  ofstream* out;

  if(!app){
    out=new ofstream(outfname.data());
    *out<<"% N\t J\t g\t h\t D\t <H(2)^2> \t <H(2)>\t ...<H(L-1)^2> \t <H(L-1)>\t <H^2>\t <H> "<<endl;
  }
  else
    out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  cout<<setprecision(10);

  int D=D2; // to start
  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  cout<<"Initialized Contractor"<<endl;
  MPS gs(L,D,d);
  if(readMPS){ // import
    gs.importMPS(inMPS.data());
    cout<<"Read initial MPS from file "<<endl;
    saveMPS=0;
    if(D!=gs.getBond()){
      cout<<"The bond dim of the read MPS is "<<gs.getBond()<<" and D2="<<D2<<endl;
      exit(1);
    }
  }
  else{
    gs.setRandomState(); // the intial state, random
    gs.gaugeCond('R',1);
    gs.gaugeCond('L',1);
    cout<<"Initialized random state, norm "<<contractor.contract(gs,gs)<<endl;
  }
  
  // First: put H and find the GS
  IsingHamiltonian hamI(L,d,J,g,h); 
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamI.getHMPO(0.);
  cout<<"Constructed the hamil MPO"<<endl;
  
  //hamil.exportForMatlab("isingHmpo.m");
  cout<<"Initial value of E, with initial state"<<endl;
  cout<<contractor.contract(gs,hamil,gs)<<endl;

  const MPO* ptrs[]={&hamil,&hamil};
  MPO ham2(L);
  MPO::join(2,ptrs,ham2);// MPO for H^2
  //  MPO ham2(2,ptrs); 
  
  bool done=0;
  double expH2=L*100.;
  // place for other params
  vector<complex_t> valsE,valsE2;
  while(!done){ // find GS of H^2 until D is large enough
    if(optim){ // try to sweep every time, only using the trncation as initial guess
      contractor.findGroundState(ham2,D,&expH2,gs);
    }
    else{// just truncate
    if(D==D2){ // even if I read it, I sweep
      contractor.findGroundState(ham2,D,&expH2,gs);
      if(saveMPS){
    	stringstream nameMPSfile;
    	nameMPSfile<<outMPS;
    	//cout<<"MPS file is "<<nameMPSfile.str()<<endl;exit(1);
    	nameMPSfile<<"_D"<<D;
    	//cout<<"MPS file is "<<nameMPSfile.str()<<endl;exit(1);
    	gs.exportMPS(nameMPSfile.str().data());
      }
    }
    else{
      const MPS aux(gs);gs=MPS(aux,D); // simple truncation
      cout<<"About to optimize the truncation from "<<D+incrD<<" to "<<D<<"; norm orig="
    	  <<contractor.contract(aux,aux)<<", norm truncated="<<contractor.contract(gs,gs)
    	  <<"; overlap="<<contractor.contract(aux,gs)<<endl;
      contractor.optimizeMPS(aux,gs,D);
      gs.gaugeCond('R',1);
    }
    }
    // now compute all expectation values <H(x)^2> and <H(x)> for cuts on the left
    computeValues(gs,hamil,valsE2,valsE);
    cout<<"State for D="<<D<<" found with <H^2> "<<*valsE2.end()<<endl;
    //    *out<<"% N\t J\t g\t h\t D\t <H^2>\t <H_L^2>\t <H_R^2> \t <H>\t <HL>\t <HR>"<<endl;
    *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"<<D<<"\t";
    for(int k=0;k<valsE.size();k++){
      *out<<real(valsE2[k])<<"\t";
      *out<<real(valsE[k])<<"\t";
    }
    *out<<endl;
    D-=incrD;
    if(D<D1) done=true;
    //else gs.increaseBondDimension(D);
  }
  
  out->close();
  delete out;
}


void computeValues(const MPS& state,const MPO& hamil,vector<complex_t>& valH2,vector<complex_t>& valH){
  valH2.clear();valH.clear();
  int L=state.getLength();
  Contractor& contractor=Contractor::theContractor();
  // Id op
  mwArray Id=identityMatrix(d);Id.reshape(Indices(d,1,d,1));
  Operator opId(Id); // Since all TI, only need edges and middle
  // Init MPO with all Ids  
  MPO aux(L);
  // Now, to start, I set the left edge on 0 and the right one on 1
  aux.setOp(0,&hamil.getOp(0),false);
  aux.setOp(1,&hamil.getOp(L-1),false);
  for(int k=2;k<L;k++) aux.setOp(k,&opId,false); // and the rest to Id
  for(int k=2;k<=L;k++){
    // compute expectation values
    complex_t valE=contractor.contract(state,aux,state);
    complex_t valE2=contractor.contract2(aux,state); // since they are hermitian, I contract H^+ H
    valH.push_back(valE);
    valH2.push_back(valE2);
    if(k<L){
      // one by one, move the right edge one to the right, and substitute in the former position by the middle operator
      aux.setOp(k,&hamil.getOp(L-1),false);
      aux.setOp(k-1,&hamil.getOp(1),false); // if (L==2) it does not enter here
    }    
  }
  
}
