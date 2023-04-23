/**
   \file SU2_EntropyFull.cpp
   
   Little program to compute the half-chain entropies and sub-blocks up to a certain size
   
  \param <inputfile> (string) File containing the state for which I want to compute the entropy
  
  \author Stefan Kühn
  \date 17/03/2017

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

vector<double> getEntropyAllBlocks(const MPS& state,Basis which_basis);
vector<double> getEntropyAllCuts(const MPS& state,Basis which_basis,int gauge_pos=-1,int lmin=-1, int lmax=-1);
void contractLmiddle(mwArray& result, const Site& ket,char dir);
void contractLmiddle(mwArray& result, const Site& ket);
/**Given the parameters, construct the name of the state*/
string constructFileName(int N, int D,double x, double mg, int dlink, double tol, int level_nr);
double lgf(double x);
double ksum(double j1,double m1,double j2,double m2,double J,double M);
double ClebschGordan(double j1, double m1, double j2, double m2, double J, double M);
int encodeState(double j, double m1, double m2,double jmax);
int encodeState(double j, double m, double jmax);
void constructMappingMPOnew(MPO &mpo, int N,double jmax);
//Memory wise a slightly less demanding version
void constructMappingMPOsparse(MPO &mpo, int N,double jmax, int l);
void constructIsometricMPO(MPO &mpo, int N,double jmax);

int main(int argc,const char* argv[])
{
  //Code Revision (magic field, automatically filled in by subversion)
  string revision="$Rev: 1559 $";
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
  Mode computationmode = single_cut;
  Basis which_basis = no_gauge_fields;  
  int lmin=-1;
  int lmax=-1;
  if(argc > NUM_OF_PARAMS + 1)
  {
    if(atoi(argv[2])>0)
      which_basis=full_basis;
  }  
  if(argc > NUM_OF_PARAMS + 2)
  {
    if(atoi(argv[3])>0)
      computationmode=block;
  }
  if(argc > NUM_OF_PARAMS + 3)
    lmin=atoi(argv[4]);  
  if(argc > NUM_OF_PARAMS + 4)
    lmax=atoi(argv[5]);  
  
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
  
  if(which_basis==full_basis)
  {
    cout << "Unfolding state without gauge fields to full basis" << endl;
    
    MPS mps_dummy_sites,mps_hamer_basis,mps_hamer_basis_lowerD;
    
    //First insert dummy sites    
    vector<int> dims(2*N-1,4);
    for(int i=1;i<2*N-2; i+=2)
      dims[i]=1;  
    
    mps_dummy_sites=MPS(2*N-1,1,dims);
    int Dr;
    int Dmax=0;
    for (int i=0;i<N; i++)
    {
      mps_dummy_sites.replaceSite(2*i,mps.getA(i),0);
      if(i<N-1)
      {
	Dr=mps.getA(i).getDr();
	mwArray dummySite=identityMatrix(Dr);
	dummySite.reshape(Indices(1,Dr,Dr));
	mps_dummy_sites.replaceSite(2*i+1,dummySite,0);
	if(Dmax<Dr)
	  Dmax=Dr;
      }
    }
    
    //Now build the isometries and map back to hamers basis
    MPO isometries(2*N-1);
    constructIsometricMPO(isometries,N,jmax);  
    contractor.contract(isometries, mps_dummy_sites, mps_hamer_basis);  
    //Now try to cut back D
    mps_hamer_basis_lowerD.approximate(mps_hamer_basis,Dmax);
    contractor.optimize(isometries, mps_dummy_sites, mps_hamer_basis_lowerD,Dmax);
    mps_dummy_sites.clear();
    mps_hamer_basis.clear();
    mps_hamer_basis_lowerD.gaugeCond('L',true);    
    
    //Now bring it back to the full basis
    MPO map(2*N-1);
    MPS mps_full;      
    if(lmax<=0)
      constructMappingMPOnew(map,N,jmax);
    else
      constructMappingMPOsparse(map,N,jmax,lmax);
    contractor.contract(map, mps_hamer_basis_lowerD, mps_full);
    gauge_pos=2*lmax;
    lmax=2*lmax;
    lmin=2*lmin-1;
    
    //Save result in MPS again 
    mps = mps_full;
    mps_full.clear();
    
    complex_t norm=contractor.contract(mps,mps);
    cout << "Norm after unfolding: " << norm << endl;
    if(abs(norm-ONE_c)>1.0E-4)
    {
      ofstream f;
      f.open("norm_warning.txt");
      f << "#Norm after unfolding" << endl << setprecision(15) <<norm << endl;
      f.close();
    }
    
  }
  if(computationmode==block)
  {
    cout << "Computing entropies of center block" << endl;
    entropies = getEntropyAllBlocks(mps,which_basis);  
  }
  else
  {
    cout << "Computing entropies of single cut" << endl;
    entropies = getEntropyAllCuts(mps,which_basis,gauge_pos,lmin,lmax);
  }
  
  cout << "End of program..." << endl;
  return 0;
}


//Memory wise more efficient for my case as the the existing one in contractor for a block in the middle with my specific problem, as I contract differently and avoid copies
void contractLmiddle(mwArray& result,const Site& ket,char dir)
{  
  int d=ket.getd();
  int Dl=ket.getDl();
  int Dr=ket.getDr();
  mwArray aux1,aux2;
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
    switch(dir)
    {
      //Append on the right hand side
      case 'L':
      case 'l':
	aux1=ket.getA();
	aux2=ket.getA();
	//Contract MPS tensors first
	aux1.reshape(Indices(d,Dl*Dr));
	aux2.reshape(Indices(d,Dl*Dr));
	aux2.permute(Indices(2,1),true);
	aux2.multiplyRight(aux1);
	aux2.reshape(Indices(Dl,Dr,Dl,Dr));
	aux2.permute(Indices(1,3,2,4));
	aux2.reshape(Indices(Dl*Dl,Dr*Dr));
	//Now multiply to center part
	result.reshape(Indices(Dl1p*Dl2p,Dr1p*Dr2p));
	result.multiplyLeft(aux2);
	result.reshape(Indices(Dl,Dl,Dr1p,Dr2p));
	break;
      //Append on the right hand side
      case 'R':
      case 'r':
	aux1=ket.getA();
	aux2=ket.getA();
	//Contract MPS tensors first
	aux1.reshape(Indices(d,Dl*Dr));
	aux2.reshape(Indices(d,Dl*Dr));
	aux2.permute(Indices(2,1),true);
	aux2.multiplyRight(aux1);
	aux2.reshape(Indices(Dl,Dr,Dl,Dr));
	aux2.permute(Indices(1,3,2,4));
	aux2.reshape(Indices(Dl*Dl,Dr*Dr));
	//Now multiply to center part
	result.reshape(Indices(Dl1p*Dl2p,Dr1p*Dr2p));
	result.multiplyRight(aux2);
	result.reshape(Indices(Dl1p,Dl2p,Dr,Dr));
	break;
    }
  }
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

vector<double> getEntropyAllBlocks(const MPS& state,Basis which_basis)
{
  vector<double> res;
  int length=state.getLength();
  double entropy=0.;
  
  //Copy the MPS and 
  MPS auxMPS(state);
  
  //Determine the center of the chain and the positions where I start to compute the entropy from
  int posL,posR;
  int middle=(int) round(length/2.0);
  
  posL=middle;
  posR=middle+1;
  
  //Now normalize the left part up to the two center sites
  if(!auxMPS.isGaugeR()) 
  {
    for(int k=0; k<posL; k++)
    {
      //cout << "Right gauging site " << k << endl;
      auxMPS.gaugeCond(k,'R',true);
    }
  }    
  //Now the same for the right part
  for(int k=length-1;k>posR ;k--)
  {
    //cout << "Left gauging site " << k << endl;
    auxMPS.gaugeCond(k,'L',true);
  } 
    
  // D^2 x D^2 overlap matrix
  mwArray overl;   
  
  //Compute the initial center part
  contractLmiddle(overl,auxMPS.getA(posL),'l');  
  contractLmiddle(overl,auxMPS.getA(posR),'r');  
  
  
  //Now compute the entropy and alternately add a site to the left/right
  int flag=0;  
  string filename="entropy_block.txt";
  if(which_basis==full_basis)
    filename="entropy_block_full.txt";
  ofstream f;  
  f.open(filename.c_str());
  f << "#n_l\tn_r\tl\tBlock entropy of block from n_l to n_r" << endl;
  f.close();
  
  while(1)
  {
    
    int Dl=auxMPS.getA(posL).getDl();
    int Dr=auxMPS.getA(posR).getDr();
    //Reshaping the overlap matrix to (Dl*Dr)x(Dl*Dr);
    overl.permute(Indices(1,3,2,4));
    overl.reshape(Indices(Dl*Dr,Dl*Dr));
    // And perform the non orthogonal diagonalization
    vector<complex_t> Dval;mwArray U;
    wrapper::eig(overl,Dval,U); // only eigenvalues
    entropy=0.0;
    for(int k=0;k<Dl*Dr;k++)
    {
      double tmp=real(Dval[k]);
      if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
    }    
    f.open(filename.c_str(),ios_base::app);
    f << setprecision(15) << posL << "\t" << posR << "\t" << posR-posL+1 << "\t" << entropy << endl;
    f.close();
    
    cout << posL << "\t" << posR << "\t" << posR-posL+1 << "\t" << entropy << endl;
    
    res.push_back(entropy);
    
    //Now undo reshaping and permuting, as I reuse it to save some memory
    overl.reshape(Indices(Dl,Dr,Dl,Dr));    
    overl.permute(Indices(1,3,2,4));      
    
    if(flag==0 && posL>1)
    {
      posL--;
      flag=1;
      contractLmiddle(overl,auxMPS.getA(posL),'l');
    }
    else if(flag==1 && posR<(length-1))
    {
      posR++;
      flag=0;
      contractLmiddle(overl,auxMPS.getA(posR),'r');      
    }
    else
      break;   
    
  }
  return res;
}

vector<double> getEntropyAllCuts(const MPS& state,Basis which_basis,int gauge_pos,int lmin, int lmax)
{
  vector<double> res;
  int length=state.getLength();
  double entropy=0.;
   int pos=0;
  
  //Copy the MPS and 
  MPS auxMPS(state);
  
  //Now normalize the left part up to the two center sites
  if(lmin<0 && gauge_pos<0 && !auxMPS.isGaugeL())
  {
    auxMPS.gaugeCond('L',true);
    lmax = length;
  }      
  else
  {      
    if(lmin<0)
      lmin=0;
    if(lmax<0)
      lmax=length;
    if(gauge_pos<0)
      gauge_pos=length;        
    for(int k=gauge_pos-1;k>lmin;k--)
      auxMPS.gaugeCond(k,'L',true);    
  }
    
  // D^2 x D^2 overlap matrix
  mwArray overl;   
  
  //Compute the initial center part
  contractLmiddle(overl,auxMPS.getA(0),'l');  
  
  //Now compute the entropy and alternately add a site to the left/right  
  string filename="entropy_half_chain.txt";
  if(which_basis==full_basis)
    filename="entropy_half_chain_full.txt"; 
  
  ofstream f;
  f.open(filename.c_str());
  f << "#L\tHalf chain entropy, cut at L" << endl;
  f.close();
  
  int Dr;
  while(1)
  { 
    if(lmin<=pos && pos<lmax)
    {
      int Dr=auxMPS.getA(pos).getDr();
      //Reshaping the overlap matrix to (Dl*Dr)x(Dl*Dr);
      overl.permute(Indices(1,3,2,4));
      overl.reshape(Indices(1*Dr,1*Dr));
      // And perform the non orthogonal diagonalization
      vector<complex_t> Dval;mwArray U; 
      cout.flush();
      wrapper::eig(overl,Dval,U); // only eigenvalues
      entropy=0.0;
      for(int k=0;k<Dr;k++)
      {
	double tmp=real(Dval[k]);
	if(abs(tmp)>1E-12) entropy+=-tmp*log2(tmp);
      }    
      f.open(filename.c_str(),ios_base::app);
      f << setprecision(15) << pos+1 << "\t" << entropy << endl;
      f.close();
      
      res.push_back(entropy);
      
      //Now undo reshaping and permuting, as I reuse it to save some memory
      overl.reshape(Indices(1,Dr,1,Dr));    
      overl.permute(Indices(1,3,2,4));
    }
    
    if(pos<(lmax-2))
    {
      pos++;
      contractLmiddle(overl,auxMPS.getA(pos));
    }    
    else
      break;   
    
  }
  return res;
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

void constructIsometricMPO(MPO &mpo, int N,double jmax)
{
  mpo.initLength(2*N-1);
  
  int D=(int)round(2.0*jmax+1.0);
  int d=D;
    
  mwArray U,S,Sfirst,V,dummy_flux,dummy_fermi,T,Slast;
  
  U=mwArray(Indices(d,1,d,D));		U.fillWithZero();
  S=mwArray(Indices(3,D,4,D));		S.fillWithZero();
  Sfirst=mwArray(Indices(3,1,4,D));	Sfirst.fillWithZero();
  Slast=mwArray(Indices(3,D,4,1));	Slast.fillWithZero();
  V=mwArray(Indices(d,D,1,1));		V.fillWithZero();
  T=mwArray(Indices(d,D,1,D));		T.fillWithZero();
  dummy_flux=mwArray(Indices(1,1,1,1));	dummy_flux.fillWithZero();
  dummy_fermi=mwArray(Indices(4,1,4,1));dummy_fermi.fillWithZero();
  
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
  
  Slast.setElement(ONE_c,Indices(0,0,0,0));
  Slast.setElement(ONE_c,Indices(1,1,1,0));
  Slast.setElement(ONE_c,Indices(2,0,3,0));  
  
  mpo.setOp(0,new Operator(Sfirst),true);
  bool saved=false;
  int pos_saved;
  for (int i=1; i<2*N-3; i+=2)
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
  if(saved)
    mpo.setOp(2*N-3,&mpo.getOp(pos_saved),false);
  else
    mpo.setOp(2*N-3,new Operator(T),true);
  mpo.setOp(2*N-2,new Operator(Slast),true);
  
}
  

double ClebschGordan(double j1,double m1,double j2,double m2,double J,double M)
{
  //Calculates CG coefficient <j1m1j2m2|JM> for the angular momentum state 

  //Check if input is sound, otherwise return 0
  if((fabs(m1)>j1)||(fabs(m2)>j2))
  {
    //cout << "Killing coefficient, because one of the z-components is larger than the individual spins" << endl;
    return 0.0;
  }
  else if(((J<fabs(j1-j2))||(J>(j1+j2))))
  {
    //cout << "Killing coefficient, because total spin has fulfill |j1-j2| <= J <= (j1+j2)" << endl;
    return 0.0;
  }
  else if(fabs(m1+m2-M)>1.0E-5)
  {
    //cout << "Killing coefficient, because M has to be m1+m2" << endl;
    return 0.0;
  }
  else if(fabs(M)>J)
  {
    //cout << "Killing coefficient, because |M|<=J has to fulfilled" << endl;
    return 0.0;
  }
  else if(fabs(fabs(remainder(J,1.0))-fabs(remainder(M,1.0)))>1.0E-5)
  {
    //cout << "Killing coefficient, because spin seems to be half-integer/integer, but z-component belongs to integer/half-integer" << endl;
    return 0.0;
  }
  else
  {
    //Some special cases and finally the most general one    
    if(fabs(j2)<1.01E-5)
        return 1.0;
    else if(fabs(J)<1.0E-5)
    {
      return pow(-1.0,(double)(j1-m1))/sqrt(2.0*j1+1.0);
    }
    else
      return sqrt(2.0*J+1.0)*exp(0.5*(lgf(J+j1-j2)+lgf(J-j1+j2)+lgf(j1+j2-J)+lgf(J+M)+lgf(J-M)-lgf(j1+j2+J+1.0)-lgf(j1-m1)-lgf(j1+m1)-lgf(j2-m2)-lgf(j2+m2)))*ksum(j1,m1,j2,m2,J,M);
  }
}

double ksum(double j1,double m1,double j2,double m2,double J,double M)
{
  double Ck=0.0;
  double kmin=max(m1-j1,max(0.0,-j1+j2+M));
  double kmax=min(J-j1+j2,min(J+M,j2+m1+J));
  for(double k=kmin; k<=kmax; k++)
    Ck =Ck+pow(-1.0,k+j2+m2)*exp(lgf(j2+J+m1-k)+lgf(j1-m1+k)-lgf(k)-lgf(J-j1+j2-k)-lgf(J+M-k)-lgf(k+j1-j2-M));
  
  return Ck;
}

//Stirlings approximation ln(n!) = nln(n)-n+0.5ln(2pin)
double lgf(double x)
{
  if(x<170.0)
      return log(tgamma(x+1.0));
  else
      return x*log(x)-x+0.5*log(2.0*M_PIl*x)+1/(12.0*x)-1.0/360.0/pow(x,3.0)+1.0/1260.0/pow(x,5.0)-1.0/1680.0/pow(x,7.0)+1.0/1188.0/pow(x,9.0);
}

int encodeState(double j, double m1, double m2, double jmax)
{
  int count=0;
  int flag=1;
  for(double jl=0; jl<=jmax; jl=jl+0.5)
  {
    for(double m1l=-jl; m1l<=jl; m1l++)
    {
      for(double m2l=-jl; m2l<=jl; m2l++)
      {
	if((fabs(jl-j)<1E-5) && (fabs(m1l-m1)<1E-5) && (fabs(m2l-m2)<1E-5))
	{
	  flag=0;
	  return count;
	}
	else
	  count++;
      }
    }
  }
  if(flag)
    return -1;
}

int encodeState(double j, double m, double jmax)
{
  int count=0;
  int flag=1;
  for(double jl=0; jl<=jmax; jl=jl+0.5)
  {
    for(double ml=-jl; ml<=jl; ml++)
    {
	if((fabs(jl-j)<1E-5) && (fabs(ml-m)<1E-5))
	  return count;
	else
	  count++;
    }
  }
  if(flag)
    return -1;
}

vector<double> decodeState_jmmp(int num, double jmax)
{
  int count=0;  
  vector<double> res;
  for(double jl=0.0; jl<=jmax; jl=jl+0.5)
    for(double m1l=-jl; m1l<=jl; m1l++)
      for(double m2l=-jl; m2l<=jl; m2l++)
      {
	if(count == num)
	{
	  res.push_back(jl);
	  res.push_back(m1l);
	  res.push_back(m2l);
	  return res;
	}
	count++;
      }
      
  return res;  
}

vector<double> decodeState_jm(int num, double jmax)
{
  int count=0;  
  vector<double> res;
  for(double jl=0.0; jl<=jmax; jl=jl+0.5)
    for(double m1l=-jl; m1l<=jl; m1l++)      
    {
      if(count == num)
      {
	res.push_back(jl);
	res.push_back(m1l);  
	return res;
      }
      count++;
    }
	
  return res;  
}

void constructMappingMPOnew(MPO &mpo, int N, double jmax)
{
  cout << "Constructing unfolding map with jmax=" << jmax << " and N=" << N << endl;
  
  int D = encodeState(jmax,jmax,jmax)+1;
  int d_j = (int) round(2.0*jmax+1.0);
  int d_jm = encodeState(jmax,jmax,jmax)+1;
  int d_jmmp = encodeState(jmax,jmax,jmax,jmax)+1; 
  int d_fermi_in = 3;
  int d_fermi_out = 4;  
  int ind_in,ind_out,ind_left,ind_right;
  vector<double> state_left;
  
  //Set up tensors
  mwArray U,S,V,T,Sfirst,Slast,M; 
  
  U=mwArray(Indices(d_jmmp,1,d_jm,D));				U.fillWithZero();
  S=mwArray(Indices(d_fermi_out,D,d_fermi_in,3*D));		S.fillWithZero();
  V=mwArray(Indices(d_jm,3*D,d_j,1));				V.fillWithZero();
  
  Sfirst=mwArray(Indices(d_fermi_out,1,d_fermi_in,3*D));	Sfirst.fillWithZero();
  Slast=mwArray(Indices(d_fermi_out,D,d_fermi_in,1));		Slast.fillWithZero();  
  
  //First tensor, just encode the mp value and pass it to the right
  for(double j=0; j<=jmax; j+=0.5)
    for(double m=-j; m<=j; m++)
    {
      ind_in = encodeState(j,m,jmax);
      for(double mp=-j;mp<=j;mp++)
      {
	ind_out = encodeState(j,m,mp,jmax);
	//ind_right=mp+j;
	ind_right=encodeState(j,mp,jmax);
	U.setElement(ONE_c,Indices(ind_out,0,ind_in,ind_right));
      }
    }
    
  //S tensor, which simply passes on the input value from the left to the right in the right block, where
  // - block 0: value of m is unchanged
  // - block 1: value of m is decreasing by 1/2
  // - block 2: value of m is increasing by 1/2
  for(int i=0; i<D; i++)
  {    
    //Block 0
    S.setElement(ONE_c,Indices(0,i,0,i));
    S.setElement(ONE_c,Indices(3,i,2,i));
    //Block 1
    S.setElement(ONE_c,Indices(1,i,1,D+i));
    //Block 2
    S.setElement(ONE_c,Indices(2,i,1,2*D+i));
  }
  
  //Block 0
  Sfirst.setElement(ONE_c,Indices(0,0,0,0));
  Sfirst.setElement(ONE_c,Indices(3,0,2,0));
  //Block 1
  Sfirst.setElement(ONE_c,Indices(1,0,1,D+0));
  //Block 2
  Sfirst.setElement(ONE_c,Indices(2,0,1,2*D+0));
  
  //Last tensor, set correction factors and Clebsch Gordan coefficients if necessary
  double mp_left,mp,j_left,cg;
  for(double j=0.0; j<=jmax; j+=0.5)
  {
    ind_in = (int)round(2.0*j);
    for(ind_left=0; ind_left<D; ind_left++)
    {     
      //Decode the z-component to which the input index corresponds
      state_left =  decodeState_jm(ind_left,jmax);
      j_left = state_left[0];
      mp_left = state_left[1];
      
      //Block 0, j to the left is exactly the same value, now get my      
      ind_out = encodeState(j,mp_left,jmax);
      if(fabs(mp_left)<=j_left && ind_out>=0)
	V.setElement(1.0/sqrt(2.0*j+1)*ONE_c,Indices(ind_out,ind_left,ind_in,0));
            
      //Block 1, j to the left can differ by 1/2, my m-value has to be the one from the left -1/2            
      ind_out = encodeState(j,mp_left-0.5,jmax);
      if(j_left<=jmax && j_left>=0.0 && fabs(mp_left)<=j_left && ind_out>=0)
      {	
	cg = ClebschGordan(0.5,-0.5,j_left,mp_left,j,mp_left-0.5);	
	V.setElement(cg/sqrt(2.0*j+1.0)*ONE_c,Indices(ind_out,ind_left+D,ind_in,0));	
      }
      
      ind_out = encodeState(j,mp_left-0.5,jmax);
      if(j_left<=jmax && j_left<=jmax && fabs(mp_left)<=j_left && ind_out>=0)
      {	
	cg = ClebschGordan(0.5,-0.5,j_left,mp_left,j,mp_left-0.5);	
	V.setElement(cg/sqrt(2.0*j+1.0)*ONE_c,Indices(ind_out,ind_left+D,ind_in,0));	
      }
      
      //Block 2, j to the left can differ by 1/2, my m-value has to be the one from the left +1/2
      ind_out = encodeState(j,mp_left+0.5,jmax);
      if(j_left<=jmax && j_left>=0.0 && fabs(mp_left)<=j_left && ind_out>=0)
      {	
	cg = ClebschGordan(0.5,0.5,j_left,mp_left,j,mp_left+0.5);	
	V.setElement(cg/sqrt(2.0*j+1.0)*ONE_c,Indices(ind_out,ind_left+2*D,ind_in,0));	
      }
      
      ind_out = encodeState(j,mp_left+0.5,jmax);
      if(j_left<=jmax && j_left<=jmax && fabs(mp_left)<=j_left && ind_out>=0)
      {
	cg = ClebschGordan(0.5,0.5,j_left,mp_left,j,mp_left+0.5);	
	V.setElement(cg/sqrt(2.0*j+1.0)*ONE_c,Indices(ind_out,ind_left+2*D,ind_in,0));	
      }
    }
  }
  
  //Block 0
  Slast.setElement(ONE_c,Indices(0,0,0,0));
  Slast.setElement(ONE_c,Indices(3,0,2,0));
  //Block 1
  Slast.setElement(-1.0/sqrt(2.0)*ONE_c,Indices(1,encodeState(0.5,0.5,jmax),1,0));  
  //Block 2
  Slast.setElement(1.0/sqrt(2.0)*ONE_c,Indices(2,encodeState(0.5,-0.5,jmax),1,0));
  
  //To make it a single MPO, contract the first and last tensor of two local maps acting on different sites  
  U.permute(Indices(1,2,4,3));
  U.reshape(Indices(d_jmmp*D,d_jm));
  V.reshape(Indices(d_jm,d_j*3*D));  
  
  M = U*V;
  M.reshape(Indices(d_jmmp,D,3*D,d_j));
  M.permute(Indices(1,3,4,2));

  //Now allocate the MPO and fill in the tensors
  mpo.initLength(2*N-1);
    
  mpo.setOp(0,new Operator(Sfirst),true); 
  bool saved=false;
  int pos_saved;
  for(int i=1; i<2*N-3; i+=2)
  {
    if(saved)
    {
      mpo.setOp(i,&mpo.getOp(pos_saved),false);
      mpo.setOp(i+1,&mpo.getOp(pos_saved+1),false);
    }
    else
    {
      mpo.setOp(i,new Operator(M),true);    
      mpo.setOp(i+1,new Operator(S),true);
      pos_saved=i;
      saved=true;      
    }        
  }
  mpo.setOp(2*N-3,&mpo.getOp(pos_saved),false);
  mpo.setOp(2*N-2,new Operator(Slast),true);
}


void constructMappingMPOsparse(MPO &mpo, int N, double jmax,int l)
{
  cout << "Constructing sparse unfolding map with jmax=" << jmax << " and N=" << N << " up to site " << l <<  endl;
  //Check if the given l is reasonable
  if(l<1 || l>=N)
  {
    cout << "Given value of l=" << l << " is outside the reasonable range" << endl;
    exit(666);
  }
  int ind = 2*l-1;
  
  int D = encodeState(jmax,jmax,jmax)+1;
  int d_j = (int) round(2.0*jmax+1.0);
  int d_jm = encodeState(jmax,jmax,jmax)+1;
  int d_jmmp = encodeState(jmax,jmax,jmax,jmax)+1; 
  int d_fermi_in = 3;
  int d_fermi_out = 4;  
  int ind_in,ind_out,ind_left,ind_right;
  vector<double> state_left;
  
  //Set up tensors
  mwArray U,S,V,T,Sfirst,Slast,M; 
  
  U=mwArray(Indices(d_jmmp,1,d_jm,D));				U.fillWithZero();
  S=mwArray(Indices(d_fermi_out,D,d_fermi_in,3*D));		S.fillWithZero();
  V=mwArray(Indices(d_jm,3*D,d_j,1));				V.fillWithZero();
  
  Sfirst=mwArray(Indices(d_fermi_out,1,d_fermi_in,3*D));	Sfirst.fillWithZero();
  Slast=mwArray(Indices(d_fermi_out,D,d_fermi_in,1));		Slast.fillWithZero();  
  
  //First tensor, just encode the mp value and pass it to the right
  for(double j=0; j<=jmax; j+=0.5)
    for(double m=-j; m<=j; m++)
    {
      ind_in = encodeState(j,m,jmax);
      for(double mp=-j;mp<=j;mp++)
      {
	ind_out = encodeState(j,m,mp,jmax);
	//ind_right=mp+j;
	ind_right=encodeState(j,mp,jmax);
	U.setElement(ONE_c,Indices(ind_out,0,ind_in,ind_right));
      }
    }
    
  //S tensor, which simply passes on the input value from the left to the right in the right block, where
  // - block 0: value of m is unchanged
  // - block 1: value of m is decreasing by 1/2
  // - block 2: value of m is increasing by 1/2
  for(int i=0; i<D; i++)
  {    
    //Block 0
    S.setElement(ONE_c,Indices(0,i,0,i));
    S.setElement(ONE_c,Indices(3,i,2,i));
    //Block 1
    S.setElement(ONE_c,Indices(1,i,1,D+i));
    //Block 2
    S.setElement(ONE_c,Indices(2,i,1,2*D+i));
  }
  
  //Block 0
  Sfirst.setElement(ONE_c,Indices(0,0,0,0));
  Sfirst.setElement(ONE_c,Indices(3,0,2,0));
  //Block 1
  Sfirst.setElement(ONE_c,Indices(1,0,1,D+0));
  //Block 2
  Sfirst.setElement(ONE_c,Indices(2,0,1,2*D+0));
  
  //Last tensor, set correction factors and Clebsch Gordan coefficients if necessary
  double mp_left,mp,j_left,cg;
  for(double j=0.0; j<=jmax; j+=0.5)
  {
    ind_in = (int)round(2.0*j);
    for(ind_left=0; ind_left<D; ind_left++)
    {     
      //Decode the z-component to which the input index corresponds
      state_left =  decodeState_jm(ind_left,jmax);
      j_left = state_left[0];
      mp_left = state_left[1];
      
      //Block 0, j to the left is exactly the same value, now get my      
      ind_out = encodeState(j,mp_left,jmax);
      if(fabs(mp_left)<=j_left && ind_out>=0)
	V.setElement(1.0/sqrt(2.0*j+1)*ONE_c,Indices(ind_out,ind_left,ind_in,0));
            
      //Block 1, j to the left can differ by 1/2, my m-value has to be the one from the left -1/2            
      ind_out = encodeState(j,mp_left-0.5,jmax);
      if(j_left<=jmax && j_left>=0.0 && fabs(mp_left)<=j_left && ind_out>=0)
      {	
	cg = ClebschGordan(0.5,-0.5,j_left,mp_left,j,mp_left-0.5);	
	V.setElement(cg/sqrt(2.0*j+1.0)*ONE_c,Indices(ind_out,ind_left+D,ind_in,0));	
      }
      
      ind_out = encodeState(j,mp_left-0.5,jmax);
      if(j_left<=jmax && j_left<=jmax && fabs(mp_left)<=j_left && ind_out>=0)
      {	
	cg = ClebschGordan(0.5,-0.5,j_left,mp_left,j,mp_left-0.5);	
	V.setElement(cg/sqrt(2.0*j+1.0)*ONE_c,Indices(ind_out,ind_left+D,ind_in,0));	
      }
      
      //Block 2, j to the left can differ by 1/2, my m-value has to be the one from the left +1/2
      ind_out = encodeState(j,mp_left+0.5,jmax);
      if(j_left<=jmax && j_left>=0.0 && fabs(mp_left)<=j_left && ind_out>=0)
      {	
	cg = ClebschGordan(0.5,0.5,j_left,mp_left,j,mp_left+0.5);	
	V.setElement(cg/sqrt(2.0*j+1.0)*ONE_c,Indices(ind_out,ind_left+2*D,ind_in,0));	
      }
      
      ind_out = encodeState(j,mp_left+0.5,jmax);
      if(j_left<=jmax && j_left<=jmax && fabs(mp_left)<=j_left && ind_out>=0)
      {
	cg = ClebschGordan(0.5,0.5,j_left,mp_left,j,mp_left+0.5);	
	V.setElement(cg/sqrt(2.0*j+1.0)*ONE_c,Indices(ind_out,ind_left+2*D,ind_in,0));	
      }
    }
  }
  
  //Block 0
  Slast.setElement(ONE_c,Indices(0,0,0,0));
  Slast.setElement(ONE_c,Indices(3,0,2,0));
  //Block 1
  Slast.setElement(-1.0/sqrt(2.0)*ONE_c,Indices(1,encodeState(0.5,0.5,jmax),1,0));  
  //Block 2
  Slast.setElement(1.0/sqrt(2.0)*ONE_c,Indices(2,encodeState(0.5,-0.5,jmax),1,0));
  
  //To make it a single MPO, contract the first and last tensor of two local maps acting on different sites 
  mwArray Vtmp;
  Vtmp = V;
  U.permute(Indices(1,2,4,3));
  U.reshape(Indices(d_jmmp*D,d_jm));
  V.reshape(Indices(d_jm,d_j*3*D));  
  
  M = U*V;
  M.reshape(Indices(d_jmmp,D,3*D,d_j));
  M.permute(Indices(1,3,4,2));

  //Now allocate the MPO and fill in the tensors
  mpo.initLength(2*N-1);
    
  mpo.setOp(0,new Operator(Sfirst),true); 
  bool saved=false;
  int pos_saved;
  for(int i=1; i<ind; i+=2)
  {    
    if(saved)
    {
      mpo.setOp(i,&mpo.getOp(pos_saved),false);
      mpo.setOp(i+1,&mpo.getOp(pos_saved+1),false);
    }
    else
    {
      mpo.setOp(i,new Operator(M),true);    
      mpo.setOp(i+1,new Operator(S),true);
      pos_saved=i;
      saved=true;      
    }        
  }
  mpo.setOp(ind,new Operator(Vtmp),true);  
  
  //Prepare the dummy tensors and set them
  mwArray id_link,id_site;
  id_link = identityMatrix(d_j);	id_link.reshape(Indices(d_j,1,d_j,1));
  id_site = identityMatrix(3);		id_site.reshape(Indices(3,1,3,1));
  
  saved=false;
  for(int i=ind+1; i<2*N-1; i+=2)
  {
    if(saved)
    {
      mpo.setOp(i,&mpo.getOp(pos_saved),false);
      if(i<2*N-2)
	mpo.setOp(i+1,&mpo.getOp(pos_saved+1),false);
    }
    else
    {
      mpo.setOp(i,new Operator(id_site),true);
      if(i<2*N-2)
	mpo.setOp(i+1,new Operator(id_link),true);
      saved = true; 
      pos_saved=i;
    } 
  }  
}