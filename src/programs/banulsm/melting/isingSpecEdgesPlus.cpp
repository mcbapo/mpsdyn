
#include <math.h>
#include <iomanip>

#include "Properties.h"
#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "misc.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

int findTmpFile(int D,const string& mpsdir,const string& baseName,string& mpsfile);
const string mpsFileName(const string& mpsdir,const string& baseName,int D);
//const string getMPSfilename(int L,double J,double g,double h,bool GS);
const string getMPSfilename(int L,double J,double g,double h,int levelN);

/** \file isingSpecEdgesPlus.cpp
 
    runs the findGroundState routine with the MPO for the
    Ising Hamiltonian \ref <IsingHamiltonian.cpp>, to check convergence
    and the same for the largest energy eigenstate. And for the first excited state.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the energy
    \param <app> (int) whether the output file is to be kept (app==1)
*/


int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  Properties props;
  if(argc>1){
    cntr++;
    cout<<"Reading properties from command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }
  int L = props.getIntProperty("L");  
  double J = props.getDoubleProperty("J");  
  double g = props.getDoubleProperty("g");  
  double h = props.getDoubleProperty("h");  
  int D = props.getIntProperty("D");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir"); // where are the MPS
  bool app=(props.getIntProperty("append")!=0);
 
  cout<<"Initialized arguments: L="<<L
      <<", J="<<J
      <<", g="<<g
      <<", h="<<h
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% N\t J\t g\t h\t D\t level\t Energy\t Energy/N"<<endl;
    out->close();delete out;
  }
  


  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);
  cout<<"Initialized Contractor"<<endl;
  int d=2;
  MPS gs(L,D,d);

  string mpsfile=getMPSfilename(L,J,g,h,0); // basename of MPS files
  string tmpMPSfile;
  if(findTmpFile(D,mpsdir,mpsfile,tmpMPSfile)){// found some tmp or smaller D file
    cout<<"Recovered initial guess from file "<<tmpMPSfile<<endl;
    gs.importMPS(tmpMPSfile.data());
    if(gs.getLength()!=L){
      cout<<"Error! File "<<tmpMPSfile<<" does not contain valid MPS "<<endl;
      exit(1);
    }
  }
  else{
    gs.setRandomState(); // the intial state, random
    //    cout<<"Initialized random state from scratch "<<endl;
  }
 
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
   
  // First: put H and find the GS
  double E0=0.;
  IsingHamiltonian hamI(L,d,J,g,h);
  cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamI.getHMPO(0.);
  cout<<"Constructed the hamil MPO"<<endl;
  cout<<"Initial value, with initial state"<<endl;
  cout<<contractor.contract(gs,hamil,gs)<<endl;

  contractor.findGroundState(hamil,D,&E0,gs);
  cout<<"Ground state found with eigenvalue "<<E0
      <<" (Energy per particle="<<E0/L<<")"<<endl;
  tmpMPSfile=mpsFileName(mpsdir,mpsfile,D);
  gs.exportMPS(tmpMPSfile.data());

  out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);
  *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"
      <<D<<"\t"<<0<<"\t"<<E0<<"\t"<<E0/L<<"\t";
  *out<<endl;
  out->close();delete out;

  // The first excited state
  vector<MPS*> levels;
  levels.push_back(&gs);
  double E1=0;
  MPS lev(L,D,d);
  mpsfile=getMPSfilename(L,J,g,h,1); // basename of MPS files
  if(findTmpFile(D,mpsdir,mpsfile,tmpMPSfile)){// found some tmp or smaller D file
    //    cout<<"Recovered initial guess for exc from file "<<tmpMPSfile1<<endl;
    lev.importMPS(tmpMPSfile.data());
    if(lev.getLength()!=L){
      cout<<"Error! File "<<tmpMPSfile<<" does not contain valid MPS "<<endl;
      exit(1);
    }
  }
  else{
    lev.setRandomState(); // the intial state, random
    //cout<<"Initialized random state for exc from scratch "<<endl;
  }
  lev.gaugeCond('R');lev.gaugeCond('L',1);
  contractor.findNextExcitedState(hamil,D,levels,&E1,lev);
  cout<<"First excited state found with eigenvalue "<<E1
      <<" (gap="<<E1-E0<<")"<<endl;
  tmpMPSfile=mpsFileName(mpsdir,mpsfile,D);
  lev.exportMPS(tmpMPSfile.data());

  out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);
  *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"
      <<D<<"\t"<<1<<"\t"<<E1<<"\t"<<E1/L<<"\t";
  *out<<endl;
  out->close();delete out;



  // Now the same for the largest eigenvalue

  // Same for largest eigenvalue:
  double Emax=0;
  MPS exc(L,D,d);
  string mpsfile1=getMPSfilename(L,J,g,h,pow(2,L)-1); // basename of MPS files
  string tmpMPSfile1;
  if(findTmpFile(D,mpsdir,mpsfile1,tmpMPSfile1)){// found some tmp or smaller D file
    //    cout<<"Recovered initial guess for exc from file "<<tmpMPSfile1<<endl;
    exc.importMPS(tmpMPSfile1.data());
    if(exc.getLength()!=L){
      cout<<"Error! File "<<tmpMPSfile1<<" does not contain valid MPS "<<endl;
      exit(1);
    }
  }
  else{
    exc.setRandomState(); // the intial state, random
    //cout<<"Initialized random state for exc from scratch "<<endl;
  }
  exc.gaugeCond('R');exc.gaugeCond('L',1);
  {
    MPO hamilMinus(L);
    hamilMinus.setOp(0,new Operator(-1.*hamil.getOp(0).getFullData()),true);
    for(int k=1;k<L;k++){
      hamilMinus.setOp(k,&hamil.getOp(k),false);
    }
    contractor.findGroundState(hamilMinus,D,&Emax,exc);
    Emax=-Emax;
  }
  cout<<"Found max Exc of H with E="<<Emax<<endl;
  tmpMPSfile1=mpsFileName(mpsdir,mpsfile1,D);
  exc.exportMPS(tmpMPSfile1.data());

  out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);
  *out<<L<<"\t"<<J<<"\t"<<g<<"\t"<<h<<"\t"
      <<D<<"\t"<<-1<<"\t"<<Emax<<"\t"<<Emax/L<<"\t";
  *out<<endl;
  out->close();delete out;
}



int findTmpFile(int D,const string& mpsdir,const string& baseName,string& mpsfile){
  mpsfile=mpsFileName(mpsdir,baseName,D)+"_tmp";
  bool found=0;
  if(file_exists(mpsfile)){
    cout<<"findTmpFile found mpsfile="<<mpsfile<<endl;
    found=true;
    return 1;
  }
  else{// search for smaller Ds
    int cnt=D;
    bool found=0;
    while(cnt>0&&!found){
      mpsfile=mpsFileName(mpsdir,baseName,cnt);
      if(file_exists(mpsfile)){
	cout<<"findTmpFile found mpsfile="<<mpsfile<<endl;
	found=true;
	return cnt;
      }
      else{
	mpsfile="";
	cnt--;
      }
    }
  }
  if(!found) return 0;
  cout<<"ERROR: Should not be here"<<endl;
  exit(1);
}
  
const string mpsFileName(const string& mpsdir,const string& baseName,int D){
  stringstream s;
  s<<mpsdir<<"/"<<baseName;
  if(D>0) s<<"_"<<D;
  return s.str();
}

const string getMPSfilename(int L,double J,double g,double h,int level){
  stringstream s;
  if(level==0)
    s<<"gs_";
  else
    if(level==pow(2,L)-1)
      s<<"exc_";
    else
      s<<"lev_"<<level;
  s<<"_L"<<L<<"_J"<<J<<"_g"<<g<<"_h"<<h;
  return s.str();
}


