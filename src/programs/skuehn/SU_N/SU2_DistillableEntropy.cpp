/**
   \file SU2_DistillableEntropy.cpp
   
   Little program to compute the half-chain entropies and sub-blocks up to a certain size
   
  \param <inputfile> (string) File containing the state for which I want to compute the entropy
  
  \author Stefan Kühn
  \date 05/04/2017

*/


#include <cstdlib>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>
//Lib which is needed for sleep command on *nix Systems
#include "unistd.h"

#include "MPS.h"
#include "MPO.h"
#include "FoldedOperator.h"
#include "Contractor.h"
#include "Indices.h"
#include "misc.h"

#include "Properties.h"

#define NUM_OF_PARAMS 1

using namespace std;
using namespace shrt;

enum Mode {single_cut,block};
enum Basis {full_basis,no_gauge_fields};

vector<double> getEntropyAllCuts(MPS& state, int dlink);
void contractLmiddle(mwArray& result,const Site& ket);
/**Given the parameters, construct the name of the state*/
string constructFileName(int N, int D,double x, double mg, int dlink, double tol, int level_nr);
int encodeState(double j, double m1, double m2,double jmax);
int encodeState(double j, double m, double jmax);
void MapToHamerBasisMPO(MPO &mpo, int N,double jmax);

int main(int argc,const char* argv[])
{
  //Code Revision (magic field, automatically filled in by subversion)
  string revision="$Rev: 1561 $";
  //Parse revision string
  revision=revision.substr(6,4);
  
  //Case that number of input arguments is too small, display some instruction
  if(argc < NUM_OF_PARAMS + 1)
  {
#ifdef DISKSTRG
    cout << "Program compiled on " << __DATE__ << ", " << __TIME__ << " using disk support with code revision " << revision << endl;
#else
    cout << "Program compiled on " << __DATE__ << ", " << __TIME__ << " with code revision " << revision << endl;
#endif  
    cout << "Number of parameters incorrect, program will be aborted..." << endl;
    cout << "Usage:   " << argv[0] << " <inputfile>" << endl;
    return -1;
  }
  
#ifdef DISKSTRG
  char tmpdir[140];
  sprintf(tmpdir,"./SitesXXXXXX");
  char* dirid=mkdtemp(tmpdir);
  if(dirid==0)
  {
    cout<<"Error: could not create temporary directory "<<tmpdir<<endl;
    exit(1);
  }
  FileSite::setDir(tmpdir);
#endif  

  /**************************************************************************************
				 * Set up parameters*
   * ***********************************************************************************/  
  const char* infile=argv[1];
  if(!file_exists(infile))
  {
    cout << "Input file '" << infile << "' does not exist, program will be aborted" << endl;
    exit(666);
  }  
  Properties props(infile);  
  //System size
  int N = props.getIntProperty("N");
  //Bond dimension
  int D=props.getIntProperty("D");
  //Parameters, I need to construct input file name and to run the simulation
  double x=props.getDoubleProperty("x");
  double mg=props.getDoubleProperty("mg");  
  int dlink=props.getIntProperty("dlink");  
  double tol=props.getDoubleProperty("tol");
  double jmax = (dlink-1.0)/2.0;
  int gauge_pos=-1;  
  
  MPS mps;
  string mps_file;
  mps_file = constructFileName(N,D,x,mg,dlink,tol,0);  
  if(!file_exists(mps_file))
  {
    cout << "MPS file '" << mps_file << "' does not exist, program will be aborted" << endl;
    exit(666);
  }  
  mps.importMPS(mps_file.c_str());
  N=mps.getLength();  
  cout << "Successfully read MPS" << endl;  
  
  //Generate a contractor
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl; 
  
  //In case PRIMME-support is possible set eigensolver to PRIMME
#ifdef USING_PRIMME
  cout << "Setting eigensolver to PRIMME" << endl;
  contractor.setEigenSolver(primme);
#endif      
  vector<double> entropies;    
  
  //Little error check, right now, I rely on the fact that the input state has an even length
  if((N%2)!=0)
  {
    cout << "Error, received an input state of length " << N << ", however, up to now only odd input sizes are supported" << endl;
    exit(666);
  }
    
  //First insert dummy sites up to the half of the original system
  MPS mps_dummy_sites;
  int Nnew = N+(int)round(N/2.0);
  vector<int> dims(Nnew,4);
  for(int i=1;i<N; i+=2)
    dims[i]=1;  
    
  mps_dummy_sites=MPS(Nnew,1,dims);
  int Dr;
  for (int i=0;i<N; i++)
  {
    if(i<(N/2))
    {
      mps_dummy_sites.replaceSite(2*i,mps.getA(i),0);
      Dr=mps.getA(i).getDr();
      mwArray dummySite=identityMatrix(Dr);
      dummySite.reshape(Indices(1,Dr,Dr));
      mps_dummy_sites.replaceSite(2*i+1,dummySite,0);      
    }
    else
      mps_dummy_sites.replaceSite(i+N/2,mps.getA(i),0);      
  }  
  
  //Now build the isometries and map back to hamers basis
  MPO isometries(1);    
  MapToHamerBasisMPO(isometries,N,jmax);
  contractor.contract(isometries, mps_dummy_sites, mps);  
  
  cout << "Computing entropies of single cut" << endl;
  entropies = getEntropyAllCuts(mps,dlink);
  
  cout << "End of program..." << endl;
  return 0;
}

//Memorywise more efficient for the single cut case
void contractLmiddle(mwArray& result,const Site& ket)
{  
  int d=ket.getd();
  int Dl=ket.getDl();
  int Dr=ket.getDr();
  mwArray aux;
  if(result.isEmpty())
  { 
    //First contraction, contract the MPS tensor    
    result=ket.getA();
    result.reshape(Indices(d,Dl*Dr));
    result.conjugate();
    mwArray aux=ket.getA();
    aux.reshape(Indices(d,Dl*Dr));
    aux.permute(Indices(2,1));
    result.multiplyLeft(aux);
    result.reshape(Indices(Dl,Dr,Dl,Dr));
    result.permute(Indices(3,1,4,2));    
  }
  else
  {
    int Dl1p=result.getDimension(0);
    int Dl2p=result.getDimension(1);
    int Dr1p=result.getDimension(2);
    int Dr2p=result.getDimension(3);    
    //Ket Tensor
    aux=ket.getA();
    aux.permute(Indices(2,1,3));
    aux.reshape(Indices(Dl,d*Dr));
    result.reshape(Indices(Dl1p*Dl2p*Dr1p,Dr2p));
    result.multiplyRight(aux);
    result.reshape(Indices(Dl1p,Dl2p,Dr1p,d,Dr));
    result.permute(Indices(1,2,5,4,3));
    result.reshape(Indices(Dl1p*Dl2p*Dr,d*Dr1p));
    //Bra Tensor
    aux=ket.getA();    
    aux.reshape(Indices(d*Dl,Dr));
    aux.conjugate();
    result.multiplyRight(aux);
    result.reshape(Indices(Dl1p,Dl2p,Dr,Dr));
    result.permute(Indices(1,2,4,3));
  }
}

vector<double> getEntropyAllCuts(MPS& state, int dlink)
{
  vector<double> res;
  int length=state.getLength();
  double entropy=0.;
  int Nnew = state.getLength();
  int N = (int) round(2.0/3.0*Nnew);
  double p_k;
  
  //Ensure proper gauge
  if(!state.isGaugeL())
    state.gaugeCond('L',true);
  //D x D overlap matrix
  mwArray overl,overl_tmp,tensor;     
  
  //Now compute the entropy and alternately add a site to the left/right  
  string filename="entropy_distillable.txt";
  
  ofstream f;
  f.open(filename.c_str());
  f << "#L\tHalf chain entropy, cut at L" << endl;
  f.close();
  
  
  int Dr;
  vector<double> pk(dlink,0.0), Sk(dlink,0.0);
  for(int pos=0; pos<N-1; pos++)
  {    
    /*double pkterm=0.0;
    double repterm=0.0;
    double entterm=0.0;*/    
    if((pos%2)!=0)
    {
      //I add a gauge site, thus I can try to compute the distillable entanglement
      pk.clear();
      Sk.clear();
      for(int i=0; i<dlink; i++)
      {
	overl_tmp = overl;
	tensor = state.getA(pos).getA();	
	int left=tensor.getDimension(1);
	int right=tensor.getDimension(2);
	tensor = tensor.subArray(Indices(i,-1,-1));
	tensor.reshape(Indices(1,left,right));
	Site curr_site(pos,tensor);	
	//Generate a site again and compute the reduced density matrix
	contractLmiddle(overl_tmp,curr_site);
	Dr=state.getA(pos).getDr();
	//Reshaping the overlap matrix to (Dl*Dr)x(Dl*Dr);
	overl_tmp.permute(Indices(1,3,2,4));
	overl_tmp.reshape(Indices(1*Dr,1*Dr));
	// And perform the non orthogonal diagonalization
	vector<complex_t> Dval;mwArray U; 
	wrapper::eig(overl_tmp,Dval,U); // only eigenvalues
	entropy=0.0;
	
	//Get the trace, which is the "p_k"
	p_k=0.0;
	for(int k=0;k<Dr;k++)
	  p_k += real(Dval[k]);
	
	
	entropy=0.0; 
	for(int k=0;k<Dr;k++)
	{
	  //Check if the values are larger than machine precision, otherwise I always have these values causing a nan
	  if(real(Dval[k])>1.0E-16 && p_k>1.0E-16)
	  {
	    double tmp=real(Dval[k])/p_k;	    
	    if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
	  }
	}	
	pk[i] = p_k;
	Sk[i] = entropy;
	//cout << "In sector j=" << i*0.5 << ": p_k=" << p_k << ", S=" << entropy << endl;
	
// 	if(p_k>1E-12)
// 	{
// 	  pkterm -= p_k*log2(p_k);
// 	  repterm += p_k*log2(i+1.0);
// 	  entterm += p_k*entropy;
// 	}
      }
      //cout << "Total: " << pkterm + repterm + entterm << endl;
      //cout << "Total (no representation term): " << pkterm + entterm << endl;
    }
    else
    {
      for(int i=0; i<dlink; i++)
      {
	pk[i] = -1.0;
	Sk[i] = -1.0;
      }
    }
    //Add next tensor    
    contractLmiddle(overl,state.getA(pos));    
    Dr=state.getA(pos).getDr();
    //Reshaping the overlap matrix to (Dl*Dr)x(Dl*Dr);
    overl.permute(Indices(1,3,2,4));
    overl.reshape(Indices(1*Dr,1*Dr));
    // And perform the non orthogonal diagonalization
    vector<complex_t> Dval;mwArray U; 
    wrapper::eig(overl,Dval,U); // only eigenvalues
    entropy=0.0;
    for(int k=0;k<Dr;k++)
    {
      double tmp=real(Dval[k]);
      if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
    }    
    f.open(filename.c_str(),ios_base::app);
    f << setprecision(15) << pos+1 << "\t" << entropy ;
    for(int i=0; i<dlink; i++)
    {
      f << "\t" << pk[i] << "\t" << Sk[i];
    }
    f << endl;
    f.close();      
    res.push_back(entropy);
    
    /*cout << "Stotal: " << entropy << endl ;
    cout << "Difference: " << pkterm + repterm + entterm - entropy << endl;
    cout << "Difference without rep term: " << pkterm +  entterm - entropy << endl << endl;*/
      
    //Now undo reshaping and permuting, as I reuse it to save some memory
    overl.reshape(Indices(1,Dr,1,Dr));    
    overl.permute(Indices(1,3,2,4));
  }
  return res;
}

void MapToHamerBasisMPO(MPO &mpo, int N,double jmax)
{
  int Nnew = N + (int)round(N/2.0);
  mpo.initLength(Nnew);
  
  int D=(int)round(2.0*jmax+1.0);
  int d=D;
    
  mwArray U,S,Sfirst,V,dummy_flux,dummy_fermi,T,Slast;
  
  U=mwArray(Indices(d,1,d,D));		U.fillWithZero();
  S=mwArray(Indices(3,D,4,D));		S.fillWithZero();
  Sfirst=mwArray(Indices(3,1,4,D));	Sfirst.fillWithZero();
  Slast=mwArray(Indices(3,D,4,1));	Slast.fillWithZero();
  V=mwArray(Indices(d,D,1,1));		V.fillWithZero();
  T=mwArray(Indices(d,D,1,D));		T.fillWithZero();
  dummy_fermi=identityMatrix(4);	dummy_fermi.reshape(Indices(4,1,4,1));
  
  for (int i=0; i<d; i++)
  {
    //Pass on the value of j on the link
    U.setElement(ONE_c,Indices(i,0,i,i));
    //Set the value of j according to the input
    V.setElement(ONE_c,Indices(i,i,0,0));
  }
  
  for (int i=0; i<D; i++)
  {
    //Empty site, just pass on the value of j
    S.setElement(ONE_c,Indices(0,i,0,i));
    if(i>0)
      //Singly occupied site, case that flux decreases |1->
      S.setElement(ONE_c,Indices(1,i,1,i-1));
    if(i<(D-1))
      //Singly occupied site, case that flux increases |1+>
      S.setElement(ONE_c,Indices(1,i,2,i+1));
    //Doubly occupied case, just pass the flux to the right
    S.setElement(ONE_c,Indices(2,i,3,i));
   
    T.setElement(ONE_c,Indices(i,i,0,i));
  }
    
  Sfirst.setElement(ONE_c,Indices(0,0,0,0));
  Sfirst.setElement(ONE_c,Indices(1,0,2,1));
  Sfirst.setElement(ONE_c,Indices(2,0,3,0));
    
  mpo.setOp(0,new Operator(Sfirst),true);
  bool saved=false;
  int pos_saved;
  for (int i=1; i<(N-1); i+=2)
  {
    if(saved)
    {
      mpo.setOp(i,&mpo.getOp(pos_saved),false);
      mpo.setOp(i+1,&mpo.getOp(pos_saved+1),false);
    }
    else
    {
      mpo.setOp(i,new Operator(T),true);
      mpo.setOp(i+1,new Operator(S),true);
      saved=true;
      pos_saved=i;
    }   
  }
  mpo.setOp(N-1,new Operator(V),true);  
  saved=false;
  pos_saved=-1;
  for(int i=N; i<Nnew; i++)
  {    
    if(saved)
      mpo.setOp(i,&mpo.getOp(pos_saved),false);
    else
    {
      mpo.setOp(i,new Operator(dummy_fermi),true);
      pos_saved=i;
      saved=true;
    }
  }  
}
  
string constructFileName(int N, int D,double x, double mg, int dlink, double tol, int level_nr)
{
  stringstream sstm;    
  sstm.str("");
  if(level_nr==0)
    sstm << "GS_";
  else
    sstm<< level_nr << "EX_";
  sstm << "N" << N << "D" << D << "d" << dlink << "x" << (int) x << "mg" << (int) (mg*1000) << "tol" << tol << ".dat";
  return sstm.str();
}