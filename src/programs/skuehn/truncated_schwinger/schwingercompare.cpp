
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results


typedef OperatorRow MPO;
#define CHECKDIMS 1

#include "SchwingerHamiltonian.h"


#include <vector>
using namespace std;
using namespace shrt;


/** SchwingerPD runs the findGroundState routine with the MPO for the
    Schwinger Hamiltonian \ref <SchwingerHamiltonian>, and saves
    values of order parameters, to construct a phase diagram. 

    To be used as test in the new implementation, compile make test9

    Receives arguments:
    \param <L> (int) length of the chain
    \param <mg> (double) parameter \f$m/g\f$ of the Hamiltonian
    \param <x> (double) parameter \f$x\f$ of the Hamiltonian
    \param <L0> (double) parameter \f$l_0\f$ of the Hamiltonian
    \param <alpha> (double) parameter \f$\alpha\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)

    For it to work in Mac Os X, a thread has to be created which
    receives the main function as argument.
*/


#include "time.h"

complex_t computeCoefficient(const MPS& state,long int elem);
complex_t computeCoefficient(const MPS& state,const vector<int>& binZ);
vector<int> intToBin(long int z,int L);
long int binToLong(const vector<int>& refStr);
void getSzMPO(int L,MPO& Smpo);
void getLMPO(int L,MPO& Ldiff, int n);
void constructsigmazMPO(int N_total, MPO &mpo, int n);


int main(int argc,const char* argv[]){
  char filename[120];
  int error=0;
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double mg=atof(argv[++cntr]);
  double x=atof(argv[++cntr]);
  double L0=atof(argv[++cntr]);
  double alpha=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  double mu=2*mg*sqrt(x);
  // WARNING!!! If x=0, this won't work, so then take mu=mg, but this is not right either!
  if(x==0.) {
    cout<<"WARNING!!!! The coefficients in case x=0 can be a bit funny"<<endl;
    mu=mg;
  }

  cout<<"Initialized arguments: L="<<L
      <<", mg="<<mg
      <<", mu="<<mu
      <<", x="<<x
      <<", l0="<<L0
      <<", alpha="<<alpha
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<setprecision(10)<<"# N\t x\t mg\t alpha\t l0\t D\t Energy \t Energy/(2Nx)\t <Psi-bar Psi>/g\t GammaA \t Gamma5"<<endl;
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  Contractor& contractor=Contractor::theContractor();
  //    contractor.setConvTol(1E-8);
  cout<<"Initialized Contractor?"<<endl;
  int d=2;
  MPS init(L,D,d);
  MPS gs(L,D,d);
  //    cout<<"Initialized empty MPS: "<<init<<endl;
  gs.setRandomState(); // the intial state, random
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  //cout<<"Set random MPS: "<<init<<endl;
  
  // First: put H and find the GS
  int Dop=5; // bond dimension of Schwinger 
  SchwingerHamiltonian hSch(L,mu,x,L0,alpha);
  const MPO& hamil=hSch.getHMPO(); 
    
    
  double lambda=0.;
  clock_t start=clock();
  contractor.findGroundState(hamil,D,&lambda,gs);
  cout<<"Ground state found with eigenvalue "<<lambda<<endl;


  clock_t final=clock();
  cout<<"Total time in findGroundState="<<(final-start)*1./CLOCKS_PER_SEC<<endl;
  cout<<"And expectation value of H="<<contractor.contract(gs,hamil,gs)<<endl;
  cout<<"Norm="<<contractor.contract(gs,gs)<<endl;

  bool printing=true;
    if(printing){
      // Use the MPOs for Gamma^alpha and Gamma^5
      MPO GammaA(L);
      MPO Gamma5(L);
      MPO fCond(L);
      hSch.constructGammaAlphaMPO(GammaA);
      hSch.constructGamma5MPO(Gamma5);
      hSch.constructCondensateMPO(fCond);
      complex_t gammaAlpha=contractor.contract(gs,GammaA,gs);
      complex_t gamma5=contractor.contract(gs,Gamma5,gs);
      complex_t condensate=contractor.contract(gs,fCond,gs);
      cout<<"Condensate expectation value found to be "<<condensate<<endl;
      complex_t gA=gammaAlpha+(complex_t){L0+alpha,0.};
      complex_t g5=(sqrt(x)/L)*gamma5;
      *out<<setprecision(10)<<L<<"\t"<<x<<"\t"<<mg<<"\t"<<alpha<<"\t"<<L0<<"\t"<<D<<"\t"<<lambda<<"\t"<<lambda/(2*L*x)<<"\t";
      *out<<condensate<<"\t";
      *out<<gA<<"\t";
      *out<<g5<<"\t";
      *out<<endl;
      
      
      //Get the L(n)
      complex_t Lold;
      Lold.re = 0;
      Lold.im = 0;
      vector <complex_t> Lvalues;
      for(int k=0; k<L;++k)
      {
	MPO Ldiff(L);
	getLMPO(L,Ldiff,k);
	Lvalues.push_back(contractor.contract(gs,Ldiff,gs)+Lold);	
	Lold = Lvalues[k];
      }
      
      ofstream electro;
      electro.open("Ln.txt",ios::app);
      for(int k=0; k<Lvalues.size(); ++k)
      {
	  //cout << Lvalues[k] << endl;
	  electro << Lvalues[k].re <<"\t";
      }
      electro << endl;
      electro.close();
      
      
      ofstream spins;
      spins.open("Spin.txt",ios::app);
      for(int k=0;k<L;++k)
      {
	MPO Sz(L);
	constructsigmazMPO(L,Sz,k);
	spins << contractor.contract(gs,Sz,gs).re << "\t";
      }
      spins << endl;
      spins.close();
	
	
	
    

    // Now also check the expectation value of Sz and how close it is to an eigenvector
    /*MPO Sz(L);
    getSzMPO(L,Sz);
    //hamil.exportForMatlab("hamil.m");
    //Sz.exportForMatlab("Sz.m");
    complex_t sz=contractor.contract(gs,Sz,gs);
    MPO Sz2(L);
    const MPO* ptrs[2]={&Sz,&Sz};
    MPO::join(2,ptrs,Sz2);
    complex_t sz2=contractor.contract(gs,Sz2,gs);
    cout<<"<Sz>="<<sz<<", <Sz^2>="<<sz2<<", <Sz^2>-<Sz>^2="<<sz2-sz*sz<<endl;
    *out<<"% <Sz>\t <Sz^2>-<Sz>^2"<<endl;
    *out<<sz<<"\t"<<sz2-sz*sz<<endl<<endl;

    // In this new version, the coefficients of the GS in the computational basis are computed one by one,
    // if the chain is up to 20 sites long
    *out<<"% Basis element \t Coefficient"<<endl;
    // I want to use as reference the |0> state from the SC expansion, which in my case is 01010101....
    complex_t reference=ONE_c; // to check signs, I rephase them all by dividing by the first non-zero one
    vector<int> scGS(L,0); for(int k=1;k<L;k=k+2) scGS[k]=1;
    reference=computeCoefficient(gs,scGS);
    reference=(1/abs(reference))*conjugate(reference); // the phase I want to remove
    if(L<=20){
      for(long int z=0;z<pow(2,L);z++){
	vector<int> binZ=intToBin(z,L);
	for(int k=0;k<L;k++) *out<<binZ[k];
	complex_t coeff=computeCoefficient(gs,z);
	if(abs(coeff)>1E-12){
	  coeff=coeff*reference;
	}
	else coeff=ZERO_c;
	int sumZ=0; for(int k=0;k<L;k++) sumZ+=pow(-1,binZ[k]);
	if(sumZ==0) *out<<"*";
	*out<<"\t"<<coeff<<endl;
      }
    }*/



    } // end of if(printing)
  out->close();
  delete out;
return error;

}

complex_t computeCoefficient(const MPS& state,long int elem){
  int L=state.getLength();
  vector<int> binZ=intToBin(elem,L);
  return computeCoefficient(state,binZ);
}
complex_t computeCoefficient(const MPS& state,const vector<int>& binZ){
  int L=state.getLength();
  complex_t data0[]={ONE_c,ZERO_c};
  complex_t data1[]={ZERO_c,ONE_c};
  static mwArray A0(Indices(2,1,1),data0);
  static mwArray A1(Indices(2,1,1),data1);
  MPS basic(L,1,2);
  for(int k=0;k<L;k++){
    if(binZ[k]==0)
      basic.setA(k,A0);
    else if(binZ[k]==1)
      basic.setA(k,A1);
    else{
      cout<<"Error: binary sequence "<<binZ<<" not valid at pos "<<k+1<<endl;
      exit(1);
    }
  }
  Contractor& contractor=Contractor::theContractor();
  return contractor.contract(state,basic);
}



vector<int> intToBin(long int z,int L){
  long int tmp=z;
  vector<int> result(L,0);
  for(int k=L-1;k>=0;k--){
    result[k]=(int)(tmp%2);
    tmp=(tmp-result[k])/2;
  }
  return result;
}


long int binToLong(const vector<int>& refStr){
  int L=refStr.size();
  long int z=0;
  for(int k=0;k<L;k++)
    z+=refStr[k]*pow(2,L-k-1);
  return z;
}


void getSzMPO(int L,MPO& Smpo){
  Smpo.initLength(L);
  int d=2;
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  mwArray idOp=identityMatrix(d);
  mwArray spinOp(Indices(d,d),dataZ);

  mwArray Z=mwArray(Indices(2,d,d));
  Z.fillWithZero();
  for(int i1=0;i1<d;i1++)
    for(int i2=0;i2<d;i2++){
      Z.setElement(idOp.getElement(Indices(i1,i2)),Indices(0,i1,i2));
      Z.setElement(spinOp.getElement(Indices(i1,i2)),Indices(1,i1,i2));
    }
  Z.reshape(Indices(2,d*d));

  int D=2;
  mwArray C(Indices(D,D,2));
  mwArray Cl(Indices(1,D,2)),Cr(Indices(D,1,2));
  Cl.setElement(ONE_c,Indices(0,0,0));
  Cr.setElement(ONE_c,Indices(D-1,0,0));
  C.setElement(ONE_c,Indices(0,0,0));
  C.setElement(ONE_c,Indices(D-1,D-1,0));
  C.setElement(ONE_c,Indices(0,D-1,1));
  Cl.setElement(ONE_c,Indices(0,D-1,1));
  Cr.setElement(ONE_c,Indices(0,0,1));
  C.reshape(Indices(D*D,2));
  Cl.reshape(Indices(D,2));  Cr.reshape(Indices(D,2));
  mwArray term=C*Z;
  term.reshape(Indices(D,D,d,d));term.permute(Indices(3,1,4,2));
  mwArray termL=Cl*Z;
  termL.reshape(Indices(1,D,d,d));termL.permute(Indices(3,1,4,2));
  mwArray termR=Cr*Z;
  termR.reshape(Indices(D,1,d,d));termR.permute(Indices(3,1,4,2));

  Smpo.setOp(0,new Operator(termL),true);
  Smpo.setOp(L-1,new Operator(termR),true);
  if(L>2){
    Smpo.setOp(1,new Operator(term),true);
    for(int k=2;k<L-1;k++){
      Smpo.setOp(k,&Smpo.getOp(1),false);
    }
  }
}

//Does not directly calculate the L(n) but L(n)-L(n-1)=1/2[sigma_z(n)+(-1)^n] which is easier to obtain, the L(n)s can than be subsequently calculated by using L(0)=0=L(N)
void getLMPO(int L,MPO& mpo, int n){
  //cout << "Constructing GaussMPO for site " << n << endl;
 
  //Make sure MPO is empty
  mpo.initLength(L);
  
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  int sign=pow(-1,n+1);  
  sigma_z = 0.5*(sigma_z+sign*id_fermionic);
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,2,2));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<2; ++i)
  {
    for(int j=0; j<2; ++j)
    {
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
      Fermionic_operators.setElement(sigma_z.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(2,2*2));
  
  //Bond dimension for the MPO Matrices
  int D_bond_left = 2;
  int D_bond_right = 2;
  
  //Different matrices appearing
  mwArray First_fermi(Indices(1,D_bond_right,2));
  mwArray Last_fermi(Indices(D_bond_left,1,2));
  mwArray Fermi_left(Indices(D_bond_left,D_bond_right,2));
  mwArray Fermi_right(Indices(D_bond_left,D_bond_right,2));
  mwArray Fermi_Op(Indices(D_bond_left,D_bond_right,2));
  
  
  //Fill with zeros
  First_fermi.fillWithZero();
  Last_fermi.fillWithZero();
  Fermi_left.fillWithZero();
  Fermi_right.fillWithZero();
  Fermi_Op.fillWithZero();
  
  //Set elements
  Fermi_left.setElement(ONE_c,Indices(0,0,0));
  Fermi_right.setElement(ONE_c,Indices(1,1,0));

  //reshape  
  Fermi_left.reshape(Indices(D_bond_left*D_bond_right,2));
  Fermi_right.reshape(Indices(D_bond_left*D_bond_right,2));
  
  //Array for the results
  mwArray res;
  
  //First, last and operator fermionic site
  if(n==0)
  {
    //Case where sigma_z is at the first site
    First_fermi.setElement(ONE_c,Indices(0,1,1));
    Last_fermi.setElement(ONE_c,Indices(1,0,0));
    
    First_fermi.reshape(Indices(D_bond_right,2));
    Last_fermi.reshape(Indices(D_bond_left,2));
    
    res=reshape(First_fermi*Fermionic_operators,Indices(1,D_bond_right,2,2));
    mpo.setOp(0,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Last_fermi*Fermionic_operators,Indices(D_bond_left,1,2,2));
    mpo.setOp(L-1,new Operator(permute(res,Indices(3,1,4,2))),true); 
  }
  else if(n==(L-1))
  {
    //Case where sigma_z is at the last site
    First_fermi.setElement(ONE_c,Indices(0,0,0));
    Last_fermi.setElement(ONE_c,Indices(0,0,1));
    
    First_fermi.reshape(Indices(D_bond_right,2));
    Last_fermi.reshape(Indices(D_bond_left,2));
    
    res=reshape(First_fermi*Fermionic_operators,Indices(1,D_bond_right,2,2));
    mpo.setOp(0,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Last_fermi*Fermionic_operators,Indices(D_bond_left,1,2,2));
    mpo.setOp(L-1,new Operator(permute(res,Indices(3,1,4,2))),true); 
  }
  else
  {
    //Case where sigma_z is not at the edge
    First_fermi.setElement(ONE_c,Indices(0,0,0));
    Last_fermi.setElement(ONE_c,Indices(1,0,0));
    Fermi_Op.setElement(ONE_c,Indices(0,1,1));
    
    First_fermi.reshape(Indices(D_bond_right,2));
    Last_fermi.reshape(Indices(D_bond_left,2));
    Fermi_Op.reshape(Indices(D_bond_left*D_bond_right,2));
    
    res=reshape(First_fermi*Fermionic_operators,Indices(1,D_bond_right,2,2));
    mpo.setOp(0,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Last_fermi*Fermionic_operators,Indices(D_bond_left,1,2,2));
    mpo.setOp(L-1,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Fermi_Op*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
    mpo.setOp(n,new Operator(permute(res,Indices(3,1,4,2))),true);     
  }
   
  //Set fermionic identities in between
  for(int k=1;k<L-1; k++)
  {
    if(k<n)
    {
      res=reshape(Fermi_left*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
    else if(k>n)
    {
      res=reshape(Fermi_right*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
  }
//   string Result;  
//   ostringstream convert; 
//   convert << n+1; 
//   Result = convert.str();  
//   string name="Lndiff"+Result+".m";
//   mpo.exportForMatlab(name.c_str());
}

void constructsigmazMPO(int N_total, MPO &mpo, int n)
{
  cout << "Constructing sigma_z("<< n << ")" << endl;  
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,2,2));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<2; ++i)
  {
    for(int j=0; j<2; ++j)
    {
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
      Fermionic_operators.setElement(sigma_z.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(2,2*2));
  
  //Bond dimension for the MPO Matrices
  int D_bond_left = 2;
  int D_bond_right = 2;
  
  //Different matrices appearing
  mwArray First_fermi(Indices(1,D_bond_right,2));
  mwArray Last_fermi(Indices(D_bond_left,1,2));
  mwArray Fermi_left(Indices(D_bond_left,D_bond_right,2));
  mwArray Fermi_right(Indices(D_bond_left,D_bond_right,2));
  mwArray Fermi_Op(Indices(D_bond_left,D_bond_right,2));

  //Fill with zeros
  First_fermi.fillWithZero();
  Last_fermi.fillWithZero();
  Fermi_left.fillWithZero();
  Fermi_right.fillWithZero();
  Fermi_Op.fillWithZero();
  
  //Set elements
  Fermi_left.setElement(ONE_c,Indices(0,0,0));
  Fermi_right.setElement(ONE_c,Indices(1,1,0));  
  
  //reshape  
  Fermi_left.reshape(Indices(D_bond_left*D_bond_right,2));
  Fermi_right.reshape(Indices(D_bond_left*D_bond_right,2));
  
  //Array for the results
  mwArray res;
  
  //First, last and operator fermionic site
  if(n==0)
  {
    //Case where sigma_z is at the first site
    First_fermi.setElement(ONE_c,Indices(0,1,1));
    Last_fermi.setElement(ONE_c,Indices(1,0,0));
    
    First_fermi.reshape(Indices(D_bond_right,2));
    Last_fermi.reshape(Indices(D_bond_left,2));
    
    res=reshape(First_fermi*Fermionic_operators,Indices(1,D_bond_right,2,2));
    mpo.setOp(0,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Last_fermi*Fermionic_operators,Indices(D_bond_left,1,2,2));
    mpo.setOp(N_total-1,new Operator(permute(res,Indices(3,1,4,2))),true); 
  }
  else if(n==(N_total-1))
  {
    //Case where sigma_z is at the last site
    First_fermi.setElement(ONE_c,Indices(0,0,0));
    Last_fermi.setElement(ONE_c,Indices(0,0,1));
    
    First_fermi.reshape(Indices(D_bond_right,2));
    Last_fermi.reshape(Indices(D_bond_left,2));
    
    res=reshape(First_fermi*Fermionic_operators,Indices(1,D_bond_right,2,2));
    mpo.setOp(0,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Last_fermi*Fermionic_operators,Indices(D_bond_left,1,2,2));
    mpo.setOp(N_total-1,new Operator(permute(res,Indices(3,1,4,2))),true); 
  }
  else
  {
    //Case where sigma_z is not at the edge
    First_fermi.setElement(ONE_c,Indices(0,0,0));
    Last_fermi.setElement(ONE_c,Indices(1,0,0));
    Fermi_Op.setElement(ONE_c,Indices(0,1,1));
    
    First_fermi.reshape(Indices(D_bond_right,2));
    Last_fermi.reshape(Indices(D_bond_left,2));
    Fermi_Op.reshape(Indices(D_bond_left*D_bond_right,2));
    
    res=reshape(First_fermi*Fermionic_operators,Indices(1,D_bond_right,2,2));
    mpo.setOp(0,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Last_fermi*Fermionic_operators,Indices(D_bond_left,1,2,2));
    mpo.setOp(N_total-1,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Fermi_Op*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
    mpo.setOp(n,new Operator(permute(res,Indices(3,1,4,2))),true);     
  }
  
  
  
  //Set fermionic identities in between
  for(int k=1;k<N_total-1; k++)
  {
    if(k<n)
    {
      res=reshape(Fermi_left*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
    else if(k>n)
    {
      res=reshape(Fermi_right*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
  }
}
