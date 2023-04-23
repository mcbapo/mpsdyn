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

//Directly altering the state, thus memory-wise more efficient
vector<double> getEntropyAllCuts(MPS& state, int dlink, int npoints=-1);
void contractLmiddle(mwArray& result, const Site& ket);
/**Given the parameters, construct the name of the state*/
string constructFileName(int N, int D,double x, double mg, int dlink, double tol, int level_nr);
double lgf(double x);
double ksum(double j1,double m1,double j2,double m2,double J,double M);
double ClebschGordan(double j1, double m1, double j2, double m2, double J, double M);
int encodeState(double j, double m1, double m2,double jmax);
int encodeState(double j, double m, double jmax);
vector<double> decodeState_jmmp(int num, double jmax);
vector<double> decodeState_jm(int num, double jmax);
void constructMappingMPOsparse(MPO &mpo, int N, double jmax);
void MapToHamerBasisMPO(MPO &mpo, int N,double jmax);

int main(int argc,const char* argv[])
{
  //Code Revision (magic field, automatically filled in by subversion)
  string revision="$Rev: 1562 $";
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
  int npoints=-1;
  if(argc > NUM_OF_PARAMS + 1)
    npoints = atoi(argv[2]);
  
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
    
  MPS mps_dummy_sites,mps_hamer_basis,mps_hamer_basis_lowerD;
  //First insert dummy sites up to the half of the original system    
  int Nnew = N+(int)round(N/2.0),Dmax=0,Dr;
  vector<int> dims(Nnew,4);
  for(int i=1;i<N; i+=2)
    dims[i]=1;  
    
  mps_dummy_sites=MPS(Nnew,1,dims);
  for (int i=0;i<N; i++)
  {
    if(i<(N/2))
    {
      mps_dummy_sites.replaceSite(2*i,mps.getA(i),0);
      Dr=mps.getA(i).getDr();
      mwArray dummySite=identityMatrix(Dr);
      dummySite.reshape(Indices(1,Dr,Dr));
      mps_dummy_sites.replaceSite(2*i+1,dummySite,0);      
      if(Dr>Dmax)
	Dmax=Dr;
    }
    else
      mps_dummy_sites.replaceSite(i+N/2,mps.getA(i),0);      
  }    
  
  //Now build the isometries and map back to hamers basis
  MPO isometries(Nnew);    
  MapToHamerBasisMPO(isometries,N,jmax);
  contractor.contract(isometries, mps_dummy_sites, mps_hamer_basis);  
  //Now cut back D
  mps_hamer_basis_lowerD.approximate(mps_hamer_basis,Dmax);
  contractor.optimize(isometries, mps_dummy_sites, mps_hamer_basis_lowerD,Dmax);  
  mps_hamer_basis_lowerD.gaugeCond('L',true);
  mps_dummy_sites.clear();
  mps_hamer_basis.clear();  
  
  //Now bring it back to the full basis
  MPO map(1);
  MPS mps_full;    
  constructMappingMPOsparse(map,N,jmax);    
    
  contractor.contract(map, mps_hamer_basis_lowerD, mps);  
  complex_t norm=contractor.contract(mps,mps);
  cout << "Norm after unfolding: " << norm << endl;
  if(abs(norm-ONE_c)>1.0E-4)
  {
    ofstream f;
    f.open("norm_warning.txt");
    f << "#Norm after unfolding" << endl << setprecision(15) <<norm << endl;
    f.close();
  }    
  cout << "Computing entropies of single cut" << endl;  
  entropies = getEntropyAllCuts(mps,dlink,npoints);  
  
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

void constructMappingMPOsparse(MPO &mpo, int N, double jmax)
{
  //Check if the given l is reasonable
  int Nnew = N + (int)round(N/2.0);
  mpo.initLength(Nnew);
  
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
  
  mpo.setOp(0,new Operator(Sfirst),true); 
  bool saved=false;
  int pos_saved;
  for(int i=1; i<(N-1); i+=2)
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
  mpo.setOp(N-1,new Operator(Vtmp),true);  
  
  //Prepare the dummy tensors and set them
  mwArray dummy_fermi;
  dummy_fermi=identityMatrix(4);	dummy_fermi.reshape(Indices(4,1,4,1));
  
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

vector<double> getEntropyAllCuts(MPS& state, int dlink, int npoints)
{
  vector<double> res;
  int length=state.getLength();
  double entropy=0.;
  int Nnew = state.getLength();
  int N = (int) round(2.0/3.0*Nnew);
  int startpos;
  double jmax = (dlink-1.0)/2.0;
  double p_k;
  
  //Ensure proper gauge, but to the last site where I apply the map, it is anyway still gauged
  if(npoints<=0 || npoints>(N-1))
    npoints=N-1;  
  startpos = N-npoints;
  
  for(int i=N-1; i>=startpos; i--)  
    state.gaugeCond(i,'L',true);  
  
  //D x D overlap matrix
  mwArray overl,overl_tmp,tensor,subtensor;     
  
  //Now compute the entropy and alternately add a site to the left/right  
  string filename="entropy_distillable_full.txt";
  
  ofstream f;
  f.open(filename.c_str());
  f << "#L\tHalf chain entropy, cut at L" << endl;
  f.close();
  
  int Dr;
  vector<double> pk(dlink,0.0), Sk(dlink,0.0);
  vector<mwArray> rhok;
  rhok.resize(dlink);
  for(int pos=0; pos<N-1; pos++)
  {    
    /*double pkterm=0.0;
    double repterm=0.0;
    double entterm=0.0;*/    
    if((pos%2)!=0 && pos>=(startpos-1))
    {      
      //I add a gauge site, thus I can try to compute the distillable entanglement
      pk.clear();
      Sk.clear();
      for(int i=0; i<dlink; i++)
      {
	overl_tmp = overl;
	tensor = state.getA(pos).getA();
	//Extract indices belonging to the current value of, so basically |j -j -j> to |j j j>
	double jval = i/2.0;	
	int ind_lower = encodeState(jval,-jval,-jval,jmax);
	int ind_upper = encodeState(jval,jval,jval,jmax);	
	int left=tensor.getDimension(1);
	int right=tensor.getDimension(2);
	//Now extract the right subtensor
	subtensor = mwArray(Indices(ind_upper-ind_lower+1,left,right));
	for(int k=ind_lower; k<=ind_upper; k++)
	  for(int inl=0; inl<left; inl++)
	    for(int inr=0; inr<right; inr++)
	      subtensor.setElement(tensor.getElement(Indices(k,inl,inr)),Indices(k-ind_lower,inl,inr));	
	Site curr_site(pos,subtensor);
	//Generate a site again and compute the reduced density matrix
	contractLmiddle(overl_tmp,curr_site);
	Dr=state.getA(pos).getDr();
	//Reshaping the overlap matrix to (Dl*Dr)x(Dl*Dr);
	overl_tmp.permute(Indices(1,3,2,4));
	overl_tmp.reshape(Indices(1*Dr,1*Dr));
	//Save the rho_k 
	rhok[i] = overl_tmp;
	// And perform the non orthogonal diagonalization
	vector<complex_t> Dval;mwArray U; 
	wrapper::eig(overl_tmp,Dval,U); // only eigenvalues	
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
	cout << "In sector j=" << i*0.5 << ": p_k=" << p_k << ", S=" << entropy << endl;
      }
      mwArray tmp;
      for (int i=0; i<dlink; i++)
      {
	for(int k=i+1; k<dlink; k++)
	{
	  tmp = Hconjugate(rhok[i])*rhok[k];
	  if(abs(tmp.trace())>1E-10)
	    cout << "|tr(rho_" << i << "^dagger rho_" << k << ")|=" << abs(tmp.trace()) << endl;
	}
      }
      cout << endl;
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
    //Did I reach the range for which I want to compute the entropies?
    if(pos>=(startpos-1))
    {
      Dr=state.getA(pos).getDr();
      //Reshaping the overlap matrix to DrxDr
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
    
      //Now undo reshaping and permuting, as I reuse it to save some memory
      overl.reshape(Indices(1,Dr,1,Dr));    
      overl.permute(Indices(1,3,2,4));
    }
  }
  return res;
}