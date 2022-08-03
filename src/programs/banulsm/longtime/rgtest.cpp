/**
   \file rgtest.cpp
   Try the stupid contraction of the tranverse time-evolved TN
   
   \author Mari-Carmen Banuls
   \date 28/01/2014
*/

#include <math.h>
#include <iomanip>
#include <typeinfo>

#include "MPS.h"
#include "Contractor.h"
#include "IsingHamiltonian.h"
#include "LongTimeRenormalizer.h"


using namespace shrt;

/** Program rgtest 

    \param <mode> (int) initial state to use, using the following key
                 0 = product state Z+
                 1 = product state Z-
                 2 = product state X+
                 3 = product state X-
                 4 = product state Y+
                 5 = product state Y-
                 (todo others)
    
*/

int main(int argc,const char* argv[]){

  // Read input arguments
  int cntr=0;
  double J=atof(argv[++cntr]);
  double g=atof(argv[++cntr]);
  double h=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int Dx=atoi(argv[++cntr]);
  int Dt=atoi(argv[++cntr]);
  int mode=atoi(argv[++cntr]); // initial state
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  cout<<"Initialized arguments "
      <<", J="<<J
      <<", g="<<g
      <<", h="<<h
      <<", delta="<<delta
      <<", outfile="<<outfname
      <<", app="<<app
      <<", Dx="<<Dx
      <<", Dt="<<Dt<<endl;

#ifdef DISKSTRG
  //  int errD=system("rm -rf /ptmp/mpq/banulsm/*");
  char tmpdir[140];
  sprintf(tmpdir,"tmpF/SitesXXXXXX");
  // if(jobname!=0)
  //   sprintf(tmpdir,"/ptmp/mpq/banulsm/Sites%sXXXXXX",jobname);
  // else
  //   sprintf(tmpdir,"/ptmp/mpq/banulsm/SitesXXXXXX");
  char* dirid=mkdtemp(tmpdir);
  if(dirid==0){
    cout<<"Error: couldn't apply mkdtemp "<<tmpdir<<endl;
    exit(1);
  }
  FileSite::setDir(tmpdir);
#endif

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% J="<<J<<", g="<<g<<", h="<<h
	<<" Dx="<<Dx<<" Dt="<<Dt<<endl;
    *out<<"% M\t delta\t t\t Dx\t Dt\t avrg sigx"<<endl;
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d=2;

  // First: create the basic evolution operator (actually, one piece
  // of the MPO) using an auxiliary, finite Hamiltonian
  // idem for the TI initial state (for now, only product)
  int Dop=2; // bond dimension of Ising exponential mpo
  mwArray basicU,A;
  {
    IsingHamiltonian hamI(3,d,J,g,h);
    cout<<"Created the Hamiltonian"<<endl;
    basicU=hamI.getUOperator(delta,0.,false,2);
    MPS initSt(3,1,d);
    if(mode<=5)
      initSt.setProductState(ProductState(mode));
    else{
      cout<<"Other modes not yet supported"<<endl;
      exit(1);
    }
    A=initSt.getA(1).getA();
    cout<<"Site mwArray is "<<A<<endl;
  }

  // Also the operator to be measured and the identity
  mwArray idOp=identityMatrix(d);idOp.reshape(Indices(d,1,d,1));
  mwArray idOpDt=identityMatrix(Dt);idOpDt.reshape(Indices(Dt,1,Dt,1));

  complex_t sigmaX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),sigmaX);sigX.reshape(Indices(d,1,d,1));

  // Now, for M steps, apply the renormalizing step, and compute the EV
  LongTimeRenormalizer& algo=LongTimeRenormalizer::theLongTimeRenormalizer();
  int cnt=0;
  while(cnt<M){
    cout<<"Step nr "<<cnt+1<<endl;
    complex_t num,den;
    mwArray Wt;
    double lambda; // dominant eigenvector of the transfer matrix. I
		   // want it to be one, so I will normalize the site
		   // by sqrt(lambda)
    // Compute
    mwArray auxId(idOpDt);
    if(A.getDimension(0)<Dt){
      auxId=identityMatrix(A.getDimension(0));
      auxId.reshape(Indices(A.getDimension(0),1,A.getDimension(0),1));
    }
    //    algo.computeContractions(A,idOp,sigX,num,den,lambda);
    algo.computeContractions(A,auxId,sigX,num,den,lambda);
    // Write
    *out<<cnt<<"\t"<<delta<<"\t"<<cnt*delta<<"\t"<<Dx<<"\t"
	<<Dt<<"\t"<<real(num/den)<<endl;
    // Update operators
    A=1/sqrt(lambda)*A;
    //    algo.unfoldedStep2(Dx,Dt,A,basicU,sigX,idOp,Wt);
    algo.unfoldedStep2(Dx,Dt,A,basicU,sigX,auxId,Wt);
    // Increment
    cnt++;
}

}
