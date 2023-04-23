#include "SHGIGaussLawAngular.h"
#include "Indices.h"
#include "Contractor.h"

using namespace std;
using namespace shrt;


SHGIGaussLawAngular::SHGIGaussLawAngular(int N_, int d_,int D_,double mu_,double x_,double e_,double lambda_):hamil(2*N_-1),charge(2*N_-1)
{
  N_fermions = N_;
  //As a bosonic link is between two fermions, the total number of sites in the chain is given by
  N_total = 2*N_-1;
  D = D_;
  mu = mu_;
  x = x_;
  d_bose = d_;
  d_fermi = 2;
  e = e_;
  lambda = lambda_;
  
  
  cout << "Constructing gauge invariant Hamiltonian with " << endl
       << "N_fermions = " << N_fermions << endl
       << "N_total = " << N_total << endl
       << "D = " << D << endl
       << "mu = " << mu << endl
       << "x = " << x << endl
       << "d = " << d_bose << endl
       << "e = " << e << endl
       << "la = " << lambda << endl;
  
  //Construct an MPO for the hamiltonian
  initOperators();
  //initMPO();
  initChargeoperator();
}

//Destructor
SHGIGaussLawAngular::~SHGIGaussLawAngular(){
    hamil.clear();
    charge.clear();    
}

void SHGIGaussLawAngular::initOperators()
{
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  
  complex_t data_x[] = {ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigma_x(Indices(2,2),data_x);
  
  complex_t data_y[] = {ZERO_c,I_c,-1.0*I_c,ZERO_c};
  mwArray sigma_y(Indices(2,2),data_y);
  
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
  mwArray sigma_minus(Indices(2,2),data_sm);
  
  complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
  mwArray sigma_plus(Indices(2,2),data_sp);
  
  complex_t data_idpsz[] = {2.0*ONE_c,ZERO_c,ZERO_c,ZERO_c};  
  mwArray idpsz(Indices(2,2),data_idpsz); 
  
  //Provide Operators for gauge field
  //mwArray Op1 = identityMatrix(this->d_bose);
  //mwArray Op2 = identityMatrix(this->d_bose);
  mwArray Op1(Indices(this->d_bose,this->d_bose));
  mwArray Op2(Indices(this->d_bose,this->d_bose));
  mwArray Theta(Indices(this->d_bose,this->d_bose));
  mwArray L_operator(Indices(this->d_bose,this->d_bose));
  mwArray Lsq_operator(Indices(this->d_bose,this->d_bose));
  mwArray id_bosonic = identityMatrix(this->d_bose);
  
  L_operator.fillWithZero();
  Lsq_operator.fillWithZero();
  Op1.fillWithZero();
  Op2.fillWithZero();
  Theta.fillWithZero();
  
  int offset = (int) floor(this->d_bose/2.0);  
  cout << "offset = " << offset << endl;
  
  
  
  for(int i=0; i<this->d_bose; ++i){
    for(int j=0; j<this->d_bose; ++j){
      if(i==j)
      {
	L_operator.setElement((double) (i-offset), 0.0,Indices(i,j));
      }
//       else if(i==(j+1))
//       {
// 	// "Like a^\dagger"
// 	Op2.setElement(1.0,0.0,Indices(i,j));
//       }
//       else if((i+1)==j)
//       {
// 	// "Like a"
// 	Op1.setElement(1.0,0.0,Indices(i,j));
//       }
    }
  }
  
  //Construct the theta Operator
  double xval,yval,k_val;
  complex_t tmp;
  for(int i=0; i<this->d_bose; ++i)
  {
    xval = -1.0*(this->d_bose-1)/2.0 + (double) i;
    for(int j=0; j<this->d_bose; ++j)
    {
      yval = -1.0*(this->d_bose-1)/2.0 + (double) j;
      tmp=ZERO_c;
      for(int k=0; k<this->d_bose; ++k)
      {
	k_val= -1.0*(this->d_bose-1)/2.0 + (double) k;
	//cout << "x=" << xval << " y=" << y << " k=" << k_val << endl;
	tmp = tmp + 2.0*M_PI/((double) this->d_bose)*k_val*exp(-2.0*I_c*M_PI/((double) this->d_bose)*k_val*(xval-yval));
      }
      Theta.setElement(tmp/((double) this->d_bose),Indices(i,j));
    }
  }
  //Get the gauge invariant version of the operators  
  wrapper::expm(Theta,Op1,I_c);
  wrapper::expm(Theta,Op2,-1.0*I_c);
  
  //cout << "L = " << L_operator << endl;
  //cout << "Op1 = " << Op1 << endl;
  //cout << "Op2 = " << Op2 << endl;
  
  //Now get get the square of the L-operator
  Lsq_operator = L_operator*L_operator;
  
  /*cout << "L²=" << Lsq_operator << endl;
  cout << "Theta = " << Theta << endl;
  cout << "Op1 = " << Op1 << endl;
  cout << "Op2 = " << Op2 << endl;*/
  
   //Array to store the bosonic operators
   mwArray Bosonic_operators(Indices(8,this->d_bose,this->d_bose));
   Bosonic_operators.fillWithZero();
   
   for(int i=0; i<this->d_bose; ++i)
   {
     for(int j=0; j<this->d_bose; ++j)
     {
       //Identity
       Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));
       //01
       Bosonic_operators.setElement(Op1.getElement(Indices(i,j)),Indices(4,i,j));
       //02
       Bosonic_operators.setElement(Op2.getElement(Indices(i,j)),Indices(5,i,j));
       //L
       Bosonic_operators.setElement(L_operator.getElement(Indices(i,j)),Indices(6,i,j));
       //Lsquared
       Bosonic_operators.setElement(Lsq_operator.getElement(Indices(i,j)),Indices(7,i,j));
     }
   }
   
   Bosonic_operators.reshape(Indices(8,this->d_bose*this->d_bose));
  
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(8,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();
  
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma_plus
      Fermionic_operators.setElement(sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
      //sigma_minus
      Fermionic_operators.setElement(sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
      //sigma_z
      Fermionic_operators.setElement(sigma_z.getElement(Indices(i,j)),Indices(3,i,j));
    }
  }
  
  Fermionic_operators.reshape(Indices(8,this->d_fermi*this->d_fermi));
  
  /***************************************************************************
   *                   Actual construction of the MPO
   * ************************************************************************/
  
  mwArray Bmatrices;
  mwArray res;
  Operator* Local_op;
  int D_bond = 9;
  int fermion_counter = 0;
  double sign;
  
  for (int k=0; k<this->N_total; ++k)
  {
    //cout << "Constructing Hamiltonian MPO for site k = " << k << ":" << endl;       
    Bmatrices.clear();
    res.clear();
    
    if((k%2)==1)
    {
      //Bosonic site(odd sites)
      Bmatrices.resize(Indices(D_bond,D_bond,8));
      Bmatrices.fillWithZero();
      
      sign = pow(-1.0,fermion_counter);
      //cout << "Sign bosonic site: " << sign << endl;
      
      //Identiy
      Bmatrices.setElement(ONE_c,Indices(0,0,0));
      Bmatrices.setElement(ONE_c,Indices(8,8,0));
      //Op1
      Bmatrices.setElement(this->x,0,Indices(2,5,4));
      //Op2
      Bmatrices.setElement(this->x,0,Indices(3,6,5));
      //L
      Bmatrices.setElement(ONE_c,Indices(0,7,6));
      Bmatrices.setElement(-2.0*this->lambda*sign,0,Indices(0,8,6));
      Bmatrices.setElement(ONE_c,Indices(1,8,6));
      Bmatrices.setElement(ONE_c,Indices(4,8,6));     
      //Lsquared
      Bmatrices.setElement(1+2*this->lambda,0,Indices(0,8,7));  
      
      //Contract and set the operator
      Bmatrices.reshape(Indices(D_bond*D_bond,8));
      res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose,this->d_bose));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true); 
    }
    else
    {
      //Set the number of the current fermionic site and calculate sign
      fermion_counter++;
      sign = pow(-1.0,fermion_counter);
      
      //cout << "Sign fermionic site " << sign << endl;
      
      if(k==0)
      {
	//First fermionic site
	//cout<< "First fermionic site " << endl;
	Bmatrices.resize(Indices(1,D_bond,8));
	Bmatrices.fillWithZero();
	
	//Identiy
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	Bmatrices.setElement(0.5*(this->mu*sign+this->lambda),0,Indices(0,8,0));
	//sigma_plus
	Bmatrices.setElement(ONE_c,Indices(0,2,1));
	//sigma_minus
	Bmatrices.setElement(ONE_c,Indices(0,3,2));
	//sigma_z
	Bmatrices.setElement(-this->lambda,0,Indices(0,4,3));
	Bmatrices.setElement(0.5*sign*(this->mu+this->lambda),0,Indices(0,8,3));
	
	//Contract and set the operator
	Bmatrices.reshape(Indices(D_bond,8));
	res=reshape(Bmatrices*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true); 
      }
      else if(k==(N_total-1))
      {
	//Last fermionic site
	//cout << "Last fermionic site " << endl;
	Bmatrices.resize(Indices(D_bond,1,8));
	Bmatrices.fillWithZero();
	
	//Identiy
	Bmatrices.setElement(0.5*(this->mu*sign+this->lambda),0,Indices(0,0,0));
	Bmatrices.setElement(ONE_c,Indices(8,0,0));
	//sigma_plus
	Bmatrices.setElement(ONE_c,Indices(6,0,1));
	//sigma_minus
	Bmatrices.setElement(ONE_c,Indices(5,0,2));
	//sigma_z
	Bmatrices.setElement(0.5*sign*(this->mu+this->lambda),0,Indices(0,0,3));
	Bmatrices.setElement(this->lambda,0,Indices(7,0,3));
	
	//Contract and set the operator
	Bmatrices.reshape(Indices(D_bond,8));
	res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true); 
      }
      else
      {
	//Fermionic sites in between
	Bmatrices.resize(Indices(D_bond,D_bond,8));
	Bmatrices.fillWithZero();
	
	//Identiy
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	Bmatrices.setElement(0.5*(this->mu*sign+this->lambda),0,Indices(0,8,0));
	Bmatrices.setElement(-2.0*this->lambda,0,Indices(7,1,0));
	Bmatrices.setElement(ONE_c,Indices(8,8,0));
	//sigma_plus
	Bmatrices.setElement(ONE_c,Indices(0,2,1));
	Bmatrices.setElement(ONE_c,Indices(6,8,1));
	//sigma_minus
	Bmatrices.setElement(ONE_c,Indices(0,3,2));
	Bmatrices.setElement(ONE_c,Indices(5,8,2));
	//sigma_z
	Bmatrices.setElement(-this->lambda,0,Indices(0,4,3));
	Bmatrices.setElement(0.5*sign*(this->mu+this->lambda),0,Indices(0,8,3));
	Bmatrices.setElement(this->lambda,0,Indices(7,8,3));
	
	//Contract and set the operator
	Bmatrices.reshape(Indices(D_bond*D_bond,8));
	res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true); 
      }
    }      
  }
  //cout << hamil << endl;
  //hamil.exportForMatlab("HSHGIGaussLawAngular.m");
}// initOperators

void SHGIGaussLawAngular::initChargeoperator()
{
  cout << "Constructing the MPO representation of the Charge Operator" << endl;
  
  //Sigma_z for the fermionic sites
  mwArray id_fermionic = identityMatrix(2);
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1.0*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  //Identity for the bosonic sites
  mwArray id_bosonic = identityMatrix(this->d_bose);
  
  //Reshape operators for contraction later 
  
   //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,2,2));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi; ++j)
    {
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
      Fermionic_operators.setElement(sigma_z.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(2,2*2));
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(2,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<this->d_bose; ++i)
  {
    for(int j=0; j<this->d_bose; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
    }
  }
  Bosonic_operators.reshape(Indices(2,this->d_bose*this->d_bose));
  
  //Bond dimension for the MPO Matrices
  int D_bond = 2;
  
  for(int k=0; k<N_total;++k)
  {
    if((k%2)==1)
    {
      //Bosonic site(odd sites)
      mwArray Bmatrices(Indices(D_bond,D_bond,2));
      
      //Identity
      Bmatrices.setElement(ONE_c,Indices(0,0,0));
      Bmatrices.setElement(ONE_c,Indices(1,1,0));
      
      Bmatrices.reshape(Indices(D_bond*D_bond,2));
      mwArray res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose,this->d_bose));
      charge.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);      
    }
    else
    {
      //Fermionic site(even sites)
      if(k==0)
      {
	//First site
	mwArray Bmatrices(Indices(1,D_bond,2));
	Bmatrices.fillWithZero();
	
	//Identiy
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	
	//Sigma_z
	Bmatrices.setElement(ONE_c,Indices(0,1,1));
	
	Bmatrices.reshape(Indices(D_bond,2));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(1,D_bond,2,2));
	charge.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	
      }
      else if(k==(N_total-1))
      {
	//Last site
	mwArray Bmatrices(Indices(D_bond,1,2));
	
	//Identiy
	Bmatrices.setElement(ONE_c,Indices(1,0,0));
	
	//Sigma_z
	Bmatrices.setElement(ONE_c,Indices(0,0,1));
	
	Bmatrices.reshape(Indices(D_bond,2));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,1,2,2));
	charge.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else
      {
	//Sites in between
	mwArray Bmatrices(Indices(D_bond,D_bond,2));
	
	//Identiy
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	Bmatrices.setElement(ONE_c,Indices(1,1,0));
	
	//Sigma_z
	Bmatrices.setElement(ONE_c,Indices(0,1,1));
	
	Bmatrices.reshape(Indices(D_bond*D_bond,2));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,D_bond,2,2));
	charge.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      
    }
  }
}

/*Old, deprecated version for the L operator
void SHGIGaussLawAngular::constructLopMPO(MPO &mpo, int n)
{
  cout << "Constructing L("<<n<<")" << endl;
  
  //As n is a DMRG site and L(n) acts on bosonic sites only, check if n is a bosonic index
  if((n%2)!=1)
  {
    cerr << "Error in SHGIGaussLawAngular::constructLopMPO(MPO &mpo, int n): Index " << n << " is not a valid bosonic site" << endl;
    exit(1);
  }
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
  //Provide Fermionic Identiy
  mwArray id_fermionic = identityMatrix(2);
  
  //Provide Bosonic Operators
  mwArray L_operator(Indices(this->d_bose,this->d_bose));
  mwArray id_bosonic = identityMatrix(this->d_bose);
  
  L_operator.fillWithZero();
  
  int offset = (int) floor(this->d_bose/2.0);  
  
  for(int i=0; i<this->d_bose; ++i){
    for(int j=0; j<this->d_bose; ++j){
      if(i==j){
	L_operator.setElement((double) (i-offset),0.0,Indices(i,j));
      }
    }
  }
  
  //cout << "L = " << L_operator << endl;
  //Reshape operators for contraction later 
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(1,2,2));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi; ++j)
    {
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(1,this->d_fermi*this->d_fermi));
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(2,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<this->d_bose; ++i)
  {
    for(int j=0; j<this->d_bose; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
      Bosonic_operators.setElement(L_operator.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Bosonic_operators.reshape(Indices(2,this->d_bose*this->d_bose));
  
  //Bond dimension for the MPO Matrices
  int D_bond_left = 2;
  int D_bond_right = 2;
  
  //Different matrices appearing
  mwArray First_fermi(Indices(1,D_bond_right,1));
  mwArray Last_fermi(Indices(D_bond_left,1,1));
  mwArray Fermi_left(Indices(D_bond_left,D_bond_right,1));
  mwArray Fermi_right(Indices(D_bond_left,D_bond_right,1));
  mwArray Bose_left(Indices(D_bond_left,D_bond_right,2));
  mwArray Bose_right(Indices(D_bond_left,D_bond_right,2));
  mwArray Bose_Op(Indices(D_bond_left,D_bond_right,2));
  
  //Fill with zeros
  First_fermi.fillWithZero();
  Last_fermi.fillWithZero();
  Fermi_left.fillWithZero();
  Fermi_right.fillWithZero();
  Bose_left.fillWithZero();
  Bose_right.fillWithZero();
  Bose_Op.fillWithZero();
  
  //Set elements
  First_fermi.setElement(ONE_c,Indices(0,0,0));
  Last_fermi.setElement(ONE_c,Indices(1,0,0));
  Fermi_left.setElement(ONE_c,Indices(0,0,0));
  Fermi_right.setElement(ONE_c,Indices(1,1,0));
  
  Bose_left.setElement(ONE_c,Indices(0,0,0));
  Bose_right.setElement(ONE_c,Indices(1,1,0));
  Bose_Op.setElement(ONE_c,Indices(0,1,1));
  
  //reshape
  First_fermi.reshape(Indices(D_bond_right,1));
  Last_fermi.reshape(Indices(D_bond_left,1));
  Fermi_left.reshape(Indices(D_bond_left*D_bond_right,1));
  Fermi_right.reshape(Indices(D_bond_left*D_bond_right,1));
  
  Bose_left.reshape(Indices(D_bond_left*D_bond_right,2));
  Bose_right.reshape(Indices(D_bond_left*D_bond_right,2));
  Bose_Op.reshape(Indices(D_bond_left*D_bond_right,2));
  
  //Set first and last entry and the operator
  //cout << "Setting first fermi " << endl;
  mwArray res;  
  res=reshape(First_fermi*Fermionic_operators,Indices(1,D_bond_right,this->d_fermi,this->d_fermi));
  mpo.setOp(0,new Operator(permute(res,Indices(3,1,4,2))),true);
  
  //cout << "Setting last fermi " << endl;
  res.clear();
  res=reshape(Last_fermi*Fermionic_operators,Indices(D_bond_left,1,this->d_fermi,this->d_fermi));
  mpo.setOp(N_total-1,new Operator(permute(res,Indices(3,1,4,2))),true);
  
  //cout << "Setting Operator " << endl;
  res.clear();
  res=reshape(Bose_Op*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
  mpo.setOp(n,new Operator(permute(res,Indices(3,1,4,2))),true);
  
  //Set bosonic rest
  for(int k=1; k<N_total-1; k+=2)
  {
    res.clear();
    if(k<n)
    {
      //cout << "Bosonic id left" << endl;
      res=reshape(Bose_left*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
    else if(k>n)
    {
      //cout << "Bosonic id right" << endl;
      res=reshape(Bose_right*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
  }
  
  //Set fermionic rest
  for(int k=2; k<N_total-1; k+=2)
  {
    res.clear();
    if(k<n)
    {
      //cout << "Fermionic id left" << endl;
      res=reshape(Fermi_left*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true); 
    }
    else if(k>n)
    {
      //cout << "Fermionic id right" << endl;
      res=reshape(Fermi_right*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true); 
    }
  }  
}*/

//New version for the L(n)-operator (more efficient)
void SHGIGaussLawAngular::constructLopMPO(MPO &mpo, int n)
{
  cout << "Constructing L("<<n<<")" << endl;
  
  //As n is a DMRG site and L(n) acts on bosonic sites only, check if n is a bosonic index
  if((n%2)!=1)
  {
    cerr << "Error in SHGIGaussLawAngular::constructLopMPO(MPO &mpo, int n): Index " << n << " is not a valid bosonic site" << endl;
    exit(1);
  }
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
  //Provide Fermionic Identiy
  mwArray id_fermionic = identityMatrix(2);
  
  //Provide Bosonic Operators
  mwArray L_operator(Indices(this->d_bose,this->d_bose));
  mwArray id_bosonic = identityMatrix(this->d_bose);
  
  L_operator.fillWithZero();
  
  int offset = (int) floor(this->d_bose/2.0);  
  
  for(int i=0; i<this->d_bose; ++i){
    for(int j=0; j<this->d_bose; ++j){
      if(i==j){
	L_operator.setElement((double) (i-offset),0.0,Indices(i,j));
      }
    }
  }
   
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(1,2,2));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi; ++j)
    {
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(1,this->d_fermi*this->d_fermi));
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(2,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<this->d_bose; ++i)
  {
    for(int j=0; j<this->d_bose; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
      Bosonic_operators.setElement(L_operator.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Bosonic_operators.reshape(Indices(2,this->d_bose*this->d_bose));
  
  //Bond dimension for the MPO Matrices
  int D_bond_left = 1;
  int D_bond_right = 1;
  
  //Different matrices appearing
  mwArray Fermi_site(Indices(D_bond_left,D_bond_right,1));
  mwArray Bose_id_site(Indices(D_bond_left,D_bond_right,2));
  mwArray Bose_Op_site(Indices(D_bond_left,D_bond_right,2));
  
  //Fill with zeros
  Fermi_site.fillWithZero();
  Bose_id_site.fillWithZero();
  Bose_Op_site.fillWithZero();
  
  //Set elements
  Fermi_site.setElement(ONE_c,Indices(0,0,0));
  
  Bose_id_site.setElement(ONE_c,Indices(0,0,0));
  Bose_Op_site.setElement(ONE_c,Indices(0,0,1));
  
  //reshape
  Fermi_site.reshape(Indices(D_bond_left*D_bond_right,1));
  
  Bose_id_site.reshape(Indices(D_bond_left*D_bond_right,2));
  Bose_Op_site.reshape(Indices(D_bond_left*D_bond_right,2));
  
  //Array for temporary results
  mwArray res;  

  //Flags if operator is already saved
  bool bose_saved=false;
  bool fermi_saved=false;
  int first_saved_site;
  
  
  //Set operators
  for(int k=0; k<N_total; k++)
  {
    res.clear();
    if((k%2)==1)
    {
      //Bosonic site wher operator is located
      if(k==n)
      {
	res=reshape(Bose_Op_site*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else
      {
	if(!bose_saved)
	{
	  res=reshape(Bose_id_site*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
	  mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	  bose_saved = true;
	  first_saved_site = k;
	}
	else
	  mpo.setOp(k,&mpo.getOp(first_saved_site),false);
      }
    }
    else
    {
      //Fermionic site
      if(!fermi_saved)
      {
	res=reshape(Fermi_site*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	fermi_saved = true;
      }
      else
	mpo.setOp(k,&mpo.getOp(k-2),false);
    }
  } 
}//constructLopMPO

void SHGIGaussLawAngular::constructGaussMPO(MPO &mpo, int n)
{
  cout << "Constructing GaussMPO for site " << n << endl;
  
  //As n is a DMRG site and sigma_z acts on fermionic sites only, check if n is a bosonic index
  if((n%2)==1)
  {
    cerr << "Error in SHGIGaussLawAngular::constructGaussMPO(MPO &mpo, int n): Index " << n << " is not a valid fermionic site" << endl;
    exit(1);
  }
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  int sign=pow(-1,n/2+1);  
  sigma_z = 0.5*(sigma_z+sign*id_fermionic);
  
  //Provide bosonic identiy
  mwArray id_bosonic = identityMatrix(this->d_bose);
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi; ++j)
    {
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
      Fermionic_operators.setElement(sigma_z.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(2,this->d_fermi*this->d_fermi));
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(1,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<this->d_bose; ++i)
  {
    for(int j=0; j<this->d_bose; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
    }
  }
  Bosonic_operators.reshape(Indices(1,this->d_bose*this->d_bose));
  
  //Bond dimension for the MPO Matrices
  int D_bond_left = 2;
  int D_bond_right = 2;
  
  //Different matrices appearing
  mwArray First_fermi(Indices(1,D_bond_right,2));
  mwArray Last_fermi(Indices(D_bond_left,1,2));
  mwArray Fermi_left(Indices(D_bond_left,D_bond_right,2));
  mwArray Fermi_right(Indices(D_bond_left,D_bond_right,2));
  mwArray Fermi_Op(Indices(D_bond_left,D_bond_right,2));
  mwArray Bose_left(Indices(D_bond_left,D_bond_right,1));
  mwArray Bose_right(Indices(D_bond_left,D_bond_right,1));
  
  //Fill with zeros
  First_fermi.fillWithZero();
  Last_fermi.fillWithZero();
  Fermi_left.fillWithZero();
  Fermi_right.fillWithZero();
  Fermi_Op.fillWithZero();
  Bose_left.fillWithZero();
  Bose_right.fillWithZero();
  
  //Set elements
  Fermi_left.setElement(ONE_c,Indices(0,0,0));
  Fermi_right.setElement(ONE_c,Indices(1,1,0));
  
  Bose_left.setElement(ONE_c,Indices(0,0,0));
  Bose_right.setElement(ONE_c,Indices(1,1,0));
  
  //reshape  
  Bose_left.reshape(Indices(D_bond_left*D_bond_right,1));
  Bose_right.reshape(Indices(D_bond_left*D_bond_right,1));
  
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
    
    res=reshape(First_fermi*Fermionic_operators,Indices(1,D_bond_right,this->d_fermi,this->d_fermi));
    mpo.setOp(0,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Last_fermi*Fermionic_operators,Indices(D_bond_left,1,this->d_fermi,this->d_fermi));
    mpo.setOp(N_total-1,new Operator(permute(res,Indices(3,1,4,2))),true); 
  }
  else if(n==(N_total-1))
  {
    //Case where sigma_z is at the last site
    First_fermi.setElement(ONE_c,Indices(0,0,0));
    Last_fermi.setElement(ONE_c,Indices(0,0,1));
    
    First_fermi.reshape(Indices(D_bond_right,2));
    Last_fermi.reshape(Indices(D_bond_left,2));
    
    res=reshape(First_fermi*Fermionic_operators,Indices(1,D_bond_right,this->d_fermi,this->d_fermi));
    mpo.setOp(0,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Last_fermi*Fermionic_operators,Indices(D_bond_left,1,this->d_fermi,this->d_fermi));
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
    
    res=reshape(First_fermi*Fermionic_operators,Indices(1,D_bond_right,this->d_fermi,this->d_fermi));
    mpo.setOp(0,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Last_fermi*Fermionic_operators,Indices(D_bond_left,1,this->d_fermi,this->d_fermi));
    mpo.setOp(N_total-1,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Fermi_Op*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
    mpo.setOp(n,new Operator(permute(res,Indices(3,1,4,2))),true);     
  }
  
  //Set bosonic identities
  for(int k=1;k<N_total-1; k+=2)
  {
    if(k<n)
    {
      res=reshape(Bose_left*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
    else if(k>n)
    {
      res=reshape(Bose_right*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
  }
  
  //Set fermionic identities in between
  for(int k=2;k<N_total-2; k+=2)
  {
    if(k<n)
    {
      res=reshape(Fermi_left*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
    else if(k>n)
    {
      res=reshape(Fermi_right*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
  }
}

void SHGIGaussLawAngular::constructsigmazMPO(MPO &mpo, int n)
{
  cout << "Constructing sigma_z("<< n << ")" << endl;
  
  //As n is a DMRG site and sigma_z acts on fermionic sites only, check if n is a bosonic index
  if((n%2)==1)
  {
    cerr << "Error in SHGIGaussLawAngular::constructsigmazMPO(MPO &mpo, int n): Index " << n << " is not a valid fermionic site" << endl;
    exit(1);
  }
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  //Provide bosonic identiy
  mwArray id_bosonic = identityMatrix(this->d_bose);
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi; ++j)
    {
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
      Fermionic_operators.setElement(sigma_z.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(2,this->d_fermi*this->d_fermi));
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(1,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<this->d_bose; ++i)
  {
    for(int j=0; j<this->d_bose; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
    }
  }
  Bosonic_operators.reshape(Indices(1,this->d_bose*this->d_bose));
  
  //Bond dimension for the MPO Matrices
  int D_bond_left = 2;
  int D_bond_right = 2;
  
  //Different matrices appearing
  mwArray First_fermi(Indices(1,D_bond_right,2));
  mwArray Last_fermi(Indices(D_bond_left,1,2));
  mwArray Fermi_left(Indices(D_bond_left,D_bond_right,2));
  mwArray Fermi_right(Indices(D_bond_left,D_bond_right,2));
  mwArray Fermi_Op(Indices(D_bond_left,D_bond_right,2));
  mwArray Bose_left(Indices(D_bond_left,D_bond_right,1));
  mwArray Bose_right(Indices(D_bond_left,D_bond_right,1));
  
  //Fill with zeros
  First_fermi.fillWithZero();
  Last_fermi.fillWithZero();
  Fermi_left.fillWithZero();
  Fermi_right.fillWithZero();
  Fermi_Op.fillWithZero();
  Bose_left.fillWithZero();
  Bose_right.fillWithZero();
  
  //Set elements
  Fermi_left.setElement(ONE_c,Indices(0,0,0));
  Fermi_right.setElement(ONE_c,Indices(1,1,0));
  
  Bose_left.setElement(ONE_c,Indices(0,0,0));
  Bose_right.setElement(ONE_c,Indices(1,1,0));
  
  //reshape  
  Bose_left.reshape(Indices(D_bond_left*D_bond_right,1));
  Bose_right.reshape(Indices(D_bond_left*D_bond_right,1));
  
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
    
    res=reshape(First_fermi*Fermionic_operators,Indices(1,D_bond_right,this->d_fermi,this->d_fermi));
    mpo.setOp(0,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Last_fermi*Fermionic_operators,Indices(D_bond_left,1,this->d_fermi,this->d_fermi));
    mpo.setOp(N_total-1,new Operator(permute(res,Indices(3,1,4,2))),true); 
  }
  else if(n==(N_total-1))
  {
    //Case where sigma_z is at the last site
    First_fermi.setElement(ONE_c,Indices(0,0,0));
    Last_fermi.setElement(ONE_c,Indices(0,0,1));
    
    First_fermi.reshape(Indices(D_bond_right,2));
    Last_fermi.reshape(Indices(D_bond_left,2));
    
    res=reshape(First_fermi*Fermionic_operators,Indices(1,D_bond_right,this->d_fermi,this->d_fermi));
    mpo.setOp(0,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Last_fermi*Fermionic_operators,Indices(D_bond_left,1,this->d_fermi,this->d_fermi));
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
    
    res=reshape(First_fermi*Fermionic_operators,Indices(1,D_bond_right,this->d_fermi,this->d_fermi));
    mpo.setOp(0,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Last_fermi*Fermionic_operators,Indices(D_bond_left,1,this->d_fermi,this->d_fermi));
    mpo.setOp(N_total-1,new Operator(permute(res,Indices(3,1,4,2))),true); 
    res=reshape(Fermi_Op*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
    mpo.setOp(n,new Operator(permute(res,Indices(3,1,4,2))),true);     
  }
  
  
  //Set bosonic identities
  for(int k=1;k<N_total-1; k+=2)
  {
    if(k<n)
    {
      res=reshape(Bose_left*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
    else if(k>n)
    {
      res=reshape(Bose_right*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
  }
  
  //Set fermionic identities in between
  for(int k=2;k<N_total-2; k+=2)
  {
    if(k<n)
    {
      res=reshape(Fermi_left*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
    else if(k>n)
    {
      res=reshape(Fermi_right*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
  }
}

void SHGIGaussLawAngular::constructCondensateMPO(MPO &mpo)
{
  cout << "Constructing condensate MPO"<< endl;
  
  //Make sure MPO is empty
  mpo.initLength(this->N_total);
  
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  mwArray condesateOp = id_fermionic + sigma_z;
  //cout << "Condensate Op" << condesateOp << endl;
  
  //Provide bosonic identiy
  mwArray id_bosonic = identityMatrix(this->d_bose);
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi; ++j)
    {
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
      Fermionic_operators.setElement(condesateOp.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(2,this->d_fermi*this->d_fermi));
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(2,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<this->d_bose; ++i)
  {
    for(int j=0; j<this->d_bose; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
    }
  }
  Bosonic_operators.reshape(Indices(2,this->d_bose*this->d_bose));
  
  //Bond dimension for the MPO Matrices
  int D_bond_left = 2;
  int D_bond_right = 2;
  
  //Count the number of Fermions
  int count_fermi=0;
  //Sign of the current site
  double sign;
  
  mwArray res;
  Operator* Local_op;
  
  for(int k=0; k<N_total; ++k)
  {
    if((k%2)==0)
    {
      count_fermi++;
      sign = pow(-1.0,count_fermi);
      //Fermionic site
      if(k==0)
      {
	//cout << "Site " << k << " x=" << this->x << " N_fermi=" << this->N_fermions <<" sign=" << sign << " Prefactor=" << 0.5*sign*sqrt(this->x)/this->N_fermions <<endl;
	//First fermionic site
	mwArray Bmatrices(Indices(1,D_bond_right,2));
	Bmatrices.fillWithZero();
	
	//Identity
	Bmatrices.setElement(ONE_c,Indices(0,0,0));	
	//Condensate operator
	Bmatrices.setElement(0.5*sign*sqrt(this->x)/this->N_fermions,0.0,Indices(0,1,1));
	
	//Contract and set the operator
	Bmatrices.reshape(Indices(D_bond_right,2));
	res=reshape(Bmatrices*Fermionic_operators,Indices(1,D_bond_right,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	mpo.setOp(k,Local_op,true); 
      }
      else if(k==(N_total-1))
      {
	//cout << "Site " << k << " x=" << this->x << " N_fermi=" << this->N_fermions <<" sign=" << sign << " Prefactor=" << 0.5*sign*sqrt(this->x)/this->N_fermions <<endl;
	//Last fermionic site
	mwArray Bmatrices(Indices(D_bond_left,1,2));
	Bmatrices.fillWithZero();
	
	//Identity
	Bmatrices.setElement(ONE_c,Indices(1,0,0));	
	//Condensate operator
	Bmatrices.setElement(0.5*sign*sqrt(this->x)/this->N_fermions,0.0,Indices(0,0,1));
	
	//Contract and set the operator
	Bmatrices.reshape(Indices(D_bond_left,2));
	res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond_left,1,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	mpo.setOp(k,Local_op,true); 
      }
      else
      {
	//cout << "Site " << k << " x=" << this->x << " N_fermi=" << this->N_fermions <<" sign=" << sign <<" Prefactor=" << 0.5*sign*sqrt(this->x)/this->N_fermions << endl;
	//Sites in between
	mwArray Bmatrices(Indices(D_bond_left,D_bond_right,2));
	Bmatrices.fillWithZero();
	
	//Identity
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	Bmatrices.setElement(ONE_c,Indices(1,1,0));	
	//Condensate operator
	Bmatrices.setElement(0.5*sign*sqrt(this->x)/this->N_fermions,0.0,Indices(0,1,1));
	
	//Contract and set the operator
	Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));	
	res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	mpo.setOp(k,Local_op,true); 
      }
    }
    else
    {
      //cout << "Site " << k << " x=" << this->x << " N_fermi=" << this->N_fermions <<" sign=" << sign <<" Prefactor=" << 0.5*sign*sqrt(this->x)/this->N_fermions << endl;
      //Bosonic site
      mwArray Bmatrices(Indices(D_bond_left,D_bond_right,2));
      Bmatrices.fillWithZero();
      
      //Identity
      Bmatrices.setElement(ONE_c,Indices(0,0,0));
      Bmatrices.setElement(ONE_c,Indices(1,1,0));
      
      //Contract and set the operator
      Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));
      res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(k,Local_op,true); 
    }
  }  
  //cout << "Writing MPO" << endl;
  //cout << mpo << endl;
  //mpo.exportForMatlab("CondensateMPO.m");
}//CondensateMPO


//Construct an explicit MPO for the penalty term to be able to see its contribution 
void SHGIGaussLawAngular::constructPenaltyMPO(MPO &mpo, double strength)
{
  if(strength<0.0)
    strength = this->lambda;
  
  cout << "Constructing Penalty MPO with strength " << strength << endl;
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
 
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  //Provide Operators for gauge field
  mwArray L_operator(Indices(this->d_bose,this->d_bose));
  mwArray Lsq_operator(Indices(this->d_bose,this->d_bose));
  mwArray id_bosonic = identityMatrix(this->d_bose);
  
  Lsq_operator.fillWithZero();
  L_operator.fillWithZero();
  
  int offset = (int) floor(this->d_bose/2.0);  
  
  for(int i=0; i<this->d_bose; ++i){
    for(int j=0; j<this->d_bose; ++j){
      if(i==j){
	L_operator.setElement((double) (i-offset),0.0,Indices(i,j));
      }
    }
  }
  
  //Now get get the square of the L-operator
  Lsq_operator = L_operator*L_operator;
  //cout << "Lop" << L_operator << endl;
  
   //Array to store the bosonic operators
   mwArray Bosonic_operators(Indices(4,this->d_bose,this->d_bose));
   Bosonic_operators.fillWithZero();
   
   for(int i=0; i<this->d_bose; ++i)
   {
     for(int j=0; j<this->d_bose; ++j)
     {
       //Identity
       Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));
       //L
       Bosonic_operators.setElement(L_operator.getElement(Indices(i,j)),Indices(2,i,j));
       //Lsquared
       Bosonic_operators.setElement(Lsq_operator.getElement(Indices(i,j)),Indices(3,i,j));
     }
   }
   
   Bosonic_operators.reshape(Indices(4,this->d_bose*this->d_bose));  
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(4,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();
  
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma_z
      Fermionic_operators.setElement(sigma_z.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  
  Fermionic_operators.reshape(Indices(4,this->d_fermi*this->d_fermi));
  
  /***************************************************************************
   *                   Actual construction of the MPO
   * ************************************************************************/
  
  mwArray Bmatrices;
  mwArray res;
  Operator* Local_op;
  int D_bond = 5;
  int fermion_counter = 0;
  double sign;
  
  for (int k=0; k<this->N_total; ++k)
  {
    //cout << "Constructing Penalty MPO for site k = " << k << ":" << endl;       
    Bmatrices.clear();
    res.clear();
    
    if((k%2)==1)
    {
      //Bosonic site(odd sites)
      Bmatrices.resize(Indices(D_bond,D_bond,4));
      Bmatrices.fillWithZero();
      
      sign = pow(-1.0,fermion_counter);
      //cout << "Sign bosonic site: " << sign << endl;
      
      //Identiy
      Bmatrices.setElement(ONE_c,Indices(0,0,0));
      Bmatrices.setElement(ONE_c,Indices(4,4,0));
      //L
      Bmatrices.setElement(ONE_c,Indices(0,2,2));
      Bmatrices.setElement(-2.0*sign*strength*ONE_c,Indices(0,4,2));
      Bmatrices.setElement(-1.0*ONE_c,Indices(1,4,2));
      Bmatrices.setElement(-2.0*ONE_c,Indices(3,4,2));
      //Lsquared
      Bmatrices.setElement(2.0*strength*ONE_c,Indices(0,4,3));  
      
      //Contract and set the operator
      Bmatrices.reshape(Indices(D_bond*D_bond,4));
      res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose,this->d_bose));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(k,Local_op,true); 
    }
    else
    {
      //Set the number of the current fermionic site and calculate sign
      fermion_counter++;
      sign = pow(-1.0,fermion_counter);
      
      //cout << "Sign fermionic site " << sign << endl;
      
      if(k==0)
      {
	//First fermionic site
	//cout<< "First fermionic site " << endl;
	Bmatrices.resize(Indices(1,D_bond,4));
	Bmatrices.fillWithZero();
	
	//Identiy
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	Bmatrices.setElement(0.5*strength*ONE_c,Indices(0,4,0));
	
	//sigma_z
	Bmatrices.setElement(strength*ONE_c,Indices(0,1,1));
	Bmatrices.setElement(0.5*strength*sign*ONE_c,Indices(0,4,1));
	
	//Contract and set the operator
	Bmatrices.reshape(Indices(D_bond,4));
	res=reshape(Bmatrices*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	mpo.setOp(k,Local_op,true); 
      }
      else if(k==(N_total-1))
      {
	//Last fermionic site
	//cout << "Last fermionic site " << endl;
	Bmatrices.resize(Indices(D_bond,1,4));
	Bmatrices.fillWithZero();
	
	//Identiy
	Bmatrices.setElement(0.5*strength*ONE_c,Indices(0,0,0));
	Bmatrices.setElement(ONE_c,Indices(4,0,0));
	//sigma_z
	Bmatrices.setElement(0.5*sign*strength*ONE_c,Indices(0,0,1));
	Bmatrices.setElement(strength*ONE_c,Indices(2,0,1));
	
	//Contract and set the operator
	Bmatrices.reshape(Indices(D_bond,4));
	res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	mpo.setOp(k,Local_op,true); 
      }
      else
      {
	//Fermionic sites in between
	Bmatrices.resize(Indices(D_bond,D_bond,4));
	Bmatrices.fillWithZero();
	
	//Identiy
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	Bmatrices.setElement(0.5*strength*ONE_c,Indices(0,4,0));
	Bmatrices.setElement(strength*ONE_c,Indices(2,3,0));
	Bmatrices.setElement(ONE_c,Indices(4,4,0));
	//sigma_z
	Bmatrices.setElement(strength*ONE_c,Indices(0,1,1));
	Bmatrices.setElement(0.5*strength*sign*ONE_c,Indices(0,4,1));
	Bmatrices.setElement(strength*ONE_c,Indices(2,4,1));
	
	//Contract and set the operator
	Bmatrices.reshape(Indices(D_bond*D_bond,4));
	res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	mpo.setOp(k,Local_op,true); 
      }
    }      
  }
  //cout << mpo << endl;
  mpo.exportForMatlab("Penalty.m");
}//Penalty

void SHGIGaussLawAngular::constructEvolutionMPO(MPO &Uodd, MPO &Ueven, double dt)
{
  //Provide the operators
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
  mwArray sigma_minus(Indices(2,2),data_sm);
  
  complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
  mwArray sigma_plus(Indices(2,2),data_sp);
  
  complex_t data_idpsz[] = {2.0*ONE_c,ZERO_c,ZERO_c,ZERO_c};  
  mwArray idpsz(Indices(2,2),data_idpsz); 
  
  //Provide Operators for gauge field
  mwArray Op1(Indices(this->d_bose,this->d_bose));
  mwArray Op2(Indices(this->d_bose,this->d_bose));
  mwArray Theta(Indices(this->d_bose,this->d_bose));
  mwArray L_operator(Indices(this->d_bose,this->d_bose));
  mwArray Lsq_operator(Indices(this->d_bose,this->d_bose));
  mwArray id_bosonic = identityMatrix(this->d_bose);
  
  L_operator.fillWithZero();
  Lsq_operator.fillWithZero();
  Theta.fillWithZero();
  Op1.fillWithZero();
  Op2.fillWithZero();
  
  int offset = (int) floor(this->d_bose/2.0);  
  
  for(int i=0; i<this->d_bose; ++i){
    for(int j=0; j<this->d_bose; ++j){
      if(i==j)
      {
	L_operator.setElement((double) (i-offset), 0.0,Indices(i,j));
      }
    }
  }
  
  //Construct the theta Operator
  double xval,yval,k_val;
  complex_t tmp;
  for(int i=0; i<this->d_bose; ++i)
  {
    xval = -1.0*(this->d_bose-1)/2.0 + (double) i;
    for(int j=0; j<this->d_bose; ++j)
    {
      yval = -1.0*(this->d_bose-1)/2.0 + (double) j;
      tmp=ZERO_c;
      for(int k=0; k<this->d_bose; ++k)
      {
	k_val= -1.0*(this->d_bose-1)/2.0 + (double) k;
	//cout << "x=" << x << " y=" << y << " k=" << k_val << endl;
	tmp = tmp + 2.0*M_PI/((double) this->d_bose)*k_val*exp(-2.0*I_c*M_PI/((double) this->d_bose)*k_val*(xval-yval));
      }
      Theta.setElement(tmp/((double) this->d_bose),Indices(i,j));
    }
  }
  //Get the gauge invariant version of the operators  
  wrapper::expm(Theta,Op1,I_c);
  wrapper::expm(Theta,Op2,-1.0*I_c);
  
  //cout << "L = " << L_operator << endl;
  //cout << "Op1 = " << Op1 << endl;
  //cout << "Op2 = " << Op2 << endl;
  
  //Now get get the square of the L-operator
  Lsq_operator = L_operator*L_operator;
  
  //   cout << "Simga plus: " << sigma_plus << endl;
  //   cout << "Simga minus: " << sigma_minus << endl;
  //   cout << "I + s3: " << idpsz << endl;
  //   cout << "Op1: " << Op1 << endl;
  //   cout << "Op2: " << Op2 << endl;
  //   cout << "id_f: " << id_fermionic << endl;
  //   cout << "id_b: " << id_bosonic << endl;
  //   cout << "Lsq: " << Lsq_operator << endl;
  //   
  //   cout << "kron(sp,sm)=" << kron(sigma_plus,sigma_minus) << endl;
  //   cout << "kron(sm,sp)=" << kron(sigma_minus,sigma_plus) << endl;
  
  
  //Since I want to split id+s3 on every site in 0.5*(id+s3) + 0.5*(id+s3) and take one term to the even MPO and one to the odd MPO I get a problem at the first site and the last site, because one term is missing. To avoid this I have to calculate an extra MPO for the first site and an extra MPO for the site before the last site. Graphically the problem looks like this
  //x--o--x--o--x--o--x
  //|  |  |  |  |  |  |
  //-------  |  -------
  //|  |  -------  |  |
  //Note: That this scheme works, the chain has to have an even number of spin sites, otherwise it could happen, that this hast to be done for an even side MPO
  
  if((this->N_fermions%2)!=0)
    cout << "!!! WARINING, TIME EVOLUTION MPO WILL BE WRONG FOR AN ODD NUMBER OF FERMIONS !!!" << endl;
  
  //Print the kronecker products for debugging
  /*cout << "Kron1 = " << this->x*kron(sigma_plus,kron(Op1,sigma_minus)) << endl;
  cout << "Kron2 = " << this->x*kron(sigma_minus,kron(Op2,sigma_plus)) << endl;*/  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Calculate the tensorproducts for the exponentials, for FIRST ODD
  mwArray arg_odd,exp_odd;
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_odd = this->x*kron(sigma_plus,kron(Op1,sigma_minus)) + this->x*kron(sigma_minus,kron(Op2,sigma_plus)) + 0.5*this->mu*kron(-1.0*idpsz,kron(id_bosonic,id_fermionic)) + 0.5*this->mu*kron(id_fermionic,kron(id_bosonic,0.5*idpsz)) + kron(id_fermionic,kron(Lsq_operator,id_fermionic));
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  //cout << "Argument exp-odd " << arg_odd << endl;
  //cout << "Exp-odd mit dt "<< dt << ": " << exp_odd << endl;
  
  //Three matrices for the odd MPO
  mwArray W1_odd_first_site,W2_odd_first_site,W3_odd_first_site;
 
  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  mwArray U,S,V;
  wrapper::svd(exp_odd,U,S,V);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_odd_first_site = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  Indices dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  mwArray U2,S2,V2;
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_odd_first_site = U2;
  //Reshape V2 and get the last matrix of the MPO
  Indices dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_odd_first_site = V2;

  //Calculate the tensorproducts for the exponentials, for the site BEFORE THE LAST SITE
  arg_odd.clear(); exp_odd.clear();
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_odd = this->x*kron(sigma_plus,kron(Op1,sigma_minus)) + this->x*kron(sigma_minus,kron(Op2,sigma_plus)) + 0.5*this->mu*kron(-0.5*idpsz,kron(id_bosonic,id_fermionic)) + 0.5*this->mu*kron(id_fermionic,kron(id_bosonic,idpsz)) + kron(id_fermionic,kron(Lsq_operator,id_fermionic));
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  //cout << "Argument exp-odd " << arg_odd << endl;
  //cout << "Exp-odd mit dt "<< dt << ": " << exp_odd << endl;
  
  //Three matrices for the odd MPO
  mwArray W1_odd_last_site,W2_odd_last_site,W3_odd_last_site;
 
  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  U.clear(); S.clear(); V.clear();
  wrapper::svd(exp_odd,U,S,V);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_odd_last_site = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_odd_last_site = U2;
  //Reshape V2 and get the last matrix of the MPO
  dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_odd_last_site = V2;    
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  //Calculate the tensorproducts for the exponentials, now for the rest of the odd part
  arg_odd.clear(); exp_odd.clear();
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_odd = this->x*kron(sigma_plus,kron(Op1,sigma_minus)) + this->x*kron(sigma_minus,kron(Op2,sigma_plus)) + 0.5*this->mu*kron(-0.5*idpsz,kron(id_bosonic,id_fermionic)) + 0.5*this->mu*kron(id_fermionic,kron(id_bosonic,0.5*idpsz)) + kron(id_fermionic,kron(Lsq_operator,id_fermionic));
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  //cout << "Argument exp-odd " << arg_odd << endl;
  //cout << "Exp-odd mit dt "<< dt << ": " << exp_odd << endl;
  
  //Three matrices for the odd MPO
  mwArray W1_odd,W2_odd,W3_odd;
 
  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  U.clear(); S.clear(); V.clear();
  wrapper::svd(exp_odd,U,S,V);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_odd = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_odd = U2;
  //Reshape V2 and get the last matrix of the MPO
  dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_odd = V2;  
  
   
  //Three matrices for the even MPO
  mwArray W1_even,W2_even,W3_even;
  
  //Calculate the tensorproducts for the exponentials, now for even part
  mwArray arg_even,exp_even;
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_even = this->x*kron(sigma_plus,kron(Op1,sigma_minus)) + this->x*kron(sigma_minus,kron(Op2,sigma_plus)) + 0.5*this->mu*kron(0.5*idpsz,kron(id_bosonic,id_fermionic)) + 0.5*this->mu*kron(id_fermionic,kron(id_bosonic,-0.5*idpsz)) + kron(id_fermionic,kron(Lsq_operator,id_fermionic));
  wrapper::expm(arg_even,exp_even,-I_c*dt);
  //cout << "Argument exp-even " << arg_even << endl;
  //cout << "Exp-even mit dt "<< dt << ": " << exp_even << endl;
  

  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_even.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_even.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_even.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  U.clear(); S.clear(); V.clear();
  wrapper::svd(exp_even,U,S,V);
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_even = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_even = U2;
  //Reshape V2 and get the last matrix of the MPO
  dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_even = V2;  
  
  
  //Now fill the odd MPO
  Uodd.clear();
  Uodd.initLength(this->N_total);
  
  bool first=true;
  
  mwArray id_bose_mpo_odd=id_bosonic;
  mwArray id_bose_mpo_even=id_bosonic;
  mwArray id_fermi_mpo_even=id_fermionic;
  id_bose_mpo_odd.reshape(Indices(this->d_bose,1,this->d_bose,1));
  id_bose_mpo_even.reshape(Indices(this->d_bose,1,this->d_bose,1));
  id_fermi_mpo_even.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
  
  for(int i=0; i<N_total-2; i +=4)
  {
    //cout << "i=" << i << endl;
    if(i==0)
    {
      //First site of the chain (which is equal to the first fermionic site
      //cout << "Odd MPO at first site" << endl;
      Uodd.setOp(i,new Operator(W1_odd_first_site),true);
      Uodd.setOp(i+1,new Operator(W2_odd_first_site),true);
      Uodd.setOp(i+2,new Operator(W3_odd_first_site),true);
    }
    else if(i==(N_total-3))
    {
      //Site before the last fermionic site
      //cout << "Odd MPO at site before last site" << endl;
      Uodd.setOp(i,new Operator(W1_odd_last_site),true);
      Uodd.setOp(i+1,new Operator(W2_odd_last_site),true);
      Uodd.setOp(i+2,new Operator(W3_odd_last_site),true);
    }      
    else if(first && (i!=0) && i!=(N_total-3))
    {
      //cout << "First odd MPO in between" << endl;
      Uodd.setOp(i,new Operator(W1_odd),true);
      Uodd.setOp(i+1,new Operator(W2_odd),true);
      Uodd.setOp(i+2,new Operator(W3_odd),true);
      first = false;
    }
    else
    {
      //cout << "Odd MPO in between" << endl;
      Uodd.setOp(i,&Uodd.getOp(i-4));
      Uodd.setOp(i+1,&Uodd.getOp(i-3));
      Uodd.setOp(i+2,&Uodd.getOp(i-2));
    }
  }
  //Set Identities
  first = true;
  for(int i=3; i<N_total; i+=4)
  {
    //cout << "i=" << i << endl;
    if(first)
      Uodd.setOp(i,new Operator(id_bose_mpo_odd),true);
    else
      Uodd.setOp(i,&Uodd.getOp(i-4));
  }
  //cout <<"Uodd " << Uodd << endl;
  //Uodd.exportForMatlab("Uodd.m");
  
  /*MPO dummyMPO(3);
  
  dummyMPO.setOp(0,new Operator(W1_odd)); 
  dummyMPO.setOp(1,new Operator(W2_odd)); 
  dummyMPO.setOp(2,new Operator(W3_odd)); 
  
  dummyMPO.exportForMatlab("Uoddsmall.m");*/
    
  
  
  //Now fill the even MPO
  Ueven.clear();
  Ueven.initLength(this->N_total);
  
   //First and last fermonic identity
  Ueven.setOp(0,new Operator(id_fermi_mpo_even),true);
  Ueven.setOp(N_total-1,&Ueven.getOp(0));
  
  //Rest
  first = true;
  for(int i=2; i<N_total-2; i +=4)
  {
    //cout << "i=" << i << endl;
    if(first)
    {
      Ueven.setOp(i,new Operator(W1_even),true);
      Ueven.setOp(i+1,new Operator(W2_even),true);
      Ueven.setOp(i+2,new Operator(W3_even),true);
      first = false;
    }
    else
    {
      Ueven.setOp(i,&Ueven.getOp(i-4));
      Ueven.setOp(i+1,&Ueven.getOp(i-3));
      Ueven.setOp(i+2,&Ueven.getOp(i-2));
    }
  }
  //Set bosonic Identities
  first = true;
  for(int i=1; i<N_total; i+=4)
  {
    //cout << "i=" << i << endl;
    if(first)
      Ueven.setOp(i,new Operator(id_bose_mpo_even),true);
    else
      Ueven.setOp(i,&Ueven.getOp(i-4));
  }
  //Ueven.exportForMatlab("Ueven.m");

  //cout <<"Ueven " << Ueven << endl;
  
  /*dummyMPO.setOp(0,new Operator(W1_even)); 
  dummyMPO.setOp(1,new Operator(W2_even)); 
  dummyMPO.setOp(2,new Operator(W3_even)); 
  
  dummyMPO.exportForMatlab("Uevensmall.m");
  dummyMPO.clear();*/
  
  //cout << "Passt oder so " << endl;
  //Ueven.exportForMatlab("Ueven.m");
  
}//constructEvolutionMPO

MPS* SHGIGaussLawAngular::constructInitalMPS(InitState state, int extent)
{  
  //Up to know the MPS has to have the right dimensions and so on
  //initmps.clear();
  mwArray spinsiteup(Indices(2,1,1));
  mwArray spinsitedown(Indices(2,1,1));
  mwArray bosesite(Indices(this->d_bose,1,1));
  mwArray bosesiteplus(Indices(this->d_bose,1,1));
  mwArray bosesiteminus(Indices(this->d_bose,1,1));
  
  spinsiteup.fillWithZero();
  spinsitedown.fillWithZero();
  bosesite.fillWithZero();
  bosesiteplus.fillWithZero();
  bosesiteminus.fillWithZero();
  
  spinsiteup.setElement(ONE_c,Indices(0,0,0));
  spinsitedown.setElement(ONE_c,Indices(1,0,0));
  bosesite.setElement(ONE_c,Indices((this->d_bose-1)/2,0,0));
  bosesiteplus.setElement(ONE_c,Indices((this->d_bose-1)/2+1,0,0));
  bosesiteminus.setElement(ONE_c,Indices((this->d_bose-1)/2-1,0,0));
  
  MPS *res_mps;
  
  res_mps = new MPS(this->N_total,1,1);
  
  vector<int> phys_dims;
  
  for(int i=0; i<this->N_total; i++)
  {
    if((i%2)==0)
      {
	phys_dims.push_back(2);
      }
      else
      {
	phys_dims.push_back(this->d_bose);
      }
    
  }
  res_mps->adjustPhysDimensions(phys_dims);
  
  
  
  //Counter to count the number of fermionic sites
  int count=0;
  
  if(state == is_downup)
  {
    cout << "Preparing state |d>|0>|u>..." << endl;
    for(int i=0; i<this->N_total; i++)
    {
      if((i%2)==0)
      {
	//Fermionic sites
	if((count%2)==0)
	{
	  //initmps.setA(i,spinsiteup);
	  res_mps->setA(i,spinsitedown);
	}
	else
	{
	  //initmps.setA(i,spinsitedown);
	  res_mps->setA(i,spinsiteup);
	}
	count ++;
      }
      else
      {
	//Bosonic sites
	//initmps.setA(i,bosesite);
	res_mps->setA(i,bosesite);
      }
    }
  }
  else if(state == is_updown)
  {
    cout << "Preparing state |u>|0>|d>..." << endl;
    for(int i=0; i<this->N_total; i++)
    {
      if((i%2)==0)
      {
	//Fermionic sites
	if((count%2)==0)
	{
	  //initmps.setA(i,spinsitedown);
	  res_mps->setA(i,spinsiteup);
	}
	else
	{
	  //initmps.setA(i,spinsiteup);
	  res_mps->setA(i,spinsitedown);
	}
	count ++;
      }
      else
      {
	//Bosonic sites
	//initmps.setA(i,bosesite);
	res_mps->setA(i,bosesite);
      }
    }
  }
  else if(state == is_upup)
  {
    cout << "Preparing state |u>|0>|u>..." << endl;
    for(int i=0; i<this->N_total; i++)
    {
      if((i%2)==0)
      {
	//Fermionic sites
	res_mps->setA(i,spinsiteup);
      }
      else
      {
	//Bosonic sites
	//initmps.setA(i,bosesite);
	res_mps->setA(i,bosesite);
      }
    }
  }
  else if(state == is_downdown)
  {
    cout << "Preparing state |d>|0>|d>..." << endl;
    for(int i=0; i<this->N_total; i++)
    {
      if((i%2)==0)
      {
	//Fermionic sites
	res_mps->setA(i,spinsitedown);
      }
      else
      {
	//Bosonic sites
	//initmps.setA(i,bosesite);
	res_mps->setA(i,bosesite);
      }
    }
  }
  else if(state == is_string)
  {
    cout << "Preparing string approximately in the middle of the chain" << endl;
    //First, set everything |u>|0>|d>...
    for(int i=0; i<this->N_total; i++)
    {
      if((i%2)==0)
      {
	//Fermionic sites
	if((count%2)==0)
	{
	  //initmps.setA(i,spinsitedown);
	  res_mps->setA(i,spinsiteup);
	}
	else
	{
	  //initmps.setA(i,spinsiteup);
	  res_mps->setA(i,spinsitedown);
	}
	count ++;
      }
      else
      {
	//Bosonic sites
	//initmps.setA(i,bosesite);
	res_mps->setA(i,bosesite);
      }
    }   
    //Check if user specified extent and it is suitable for the chosen chain-length, otherwise just set it to 1
    if((extent<=0) || (extent >= (this->N_fermions-3)) || ((extent%2)==0))
      extent = 1;
    //Find the middle of the chain where the string should reside
    int middle = (int) (this->N_fermions/2.0);
    cout << "Middle of the chain is " << middle << endl;
    cout << "Extent is " << extent << endl;
    //Find out if it is a positive or a negative site
    int sign = pow(-1,middle- (extent-1)/2);
    //Translate into absolute site number
    int startsite = 2*(middle-(extent-1)/2)-2;
    int endsite = 2*(middle+(extent-1)/2+1)-2;    
    
    cout<< "Startsite " << startsite << ", Endsite " << endsite << endl;
    
    for(int i=startsite; i<=endsite; i++)
    {
      if((i%2)==0)
      {
	if(i==startsite)
	{
	  //Fermionic sites
	  if(sign == -1)//Case of a former spin-up-site
	    res_mps->setA(i,spinsitedown);
	  else
	    res_mps->setA(i,spinsiteup);
	}
	else if(i==endsite)
	{
	  if(sign == -1)//Startsite was spin up, so endsite is spin down and has to be flipped up
	    res_mps->setA(i,spinsiteup);
	  else
	    res_mps->setA(i,spinsitedown);
	}
      }
      else
      {
	//Bosonic sites
	if(sign==-1)
	  res_mps->setA(i,bosesiteminus);
	else
	  res_mps->setA(i,bosesiteplus);
      }
    }
  }
  else if(state == is_pzoller_string)
  {
    cout << "Preparing state |u>|-1>|d>|-1>|u>|-1>|d>..." << endl;
    //First, set everything |u>|0>|d>...
    for(int i=0; i<this->N_total; i++)
    {
      if((i%2)==0)
      {
	//Fermionic sites
	if((count%2)==0)
	{
	  //initmps.setA(i,spinsitedown);
	  res_mps->setA(i,spinsiteup);
	}
	else
	{
	  //initmps.setA(i,spinsiteup);
	  res_mps->setA(i,spinsitedown);
	}
	count ++;
      }
      else
      {
	//Bosonic sites
	res_mps->setA(i,bosesiteminus);
      }
    }
  }
  else
    cout << "Unknown state, state cannot be prepared, return empty MPS" << endl;
  
  /*cout << "Res mps:" << *res_mps << endl;
  
  for(int i=0; i<this->N_total; i++)
  {
    cout << (res_mps->getA(i)).getA() << endl;
  }*/
  
  return res_mps;
}//constructInitalMPS

void SHGIGaussLawAngular::constructEvolutionMPOefficient(MPO &Uodd, MPO &Ueven, double dt)
{
  if(!this->terms_saved)
  {
    //Case in which the single terms are not yet built and have to be generated
    //Provide the operators
    //Provide Pauli-Matrices
    mwArray id_fermionic = identityMatrix(2);
    
    complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
    mwArray sigma_z(Indices(2,2),data_z);
    
    complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
    mwArray sigma_minus(Indices(2,2),data_sm);
    
    complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
    mwArray sigma_plus(Indices(2,2),data_sp);
    
    complex_t data_idpsz[] = {2.0*ONE_c,ZERO_c,ZERO_c,ZERO_c};  
    mwArray idpsz(Indices(2,2),data_idpsz); 
    
    //Provide Operators for gauge field
    mwArray Op1(Indices(this->d_bose,this->d_bose));
    mwArray Op2(Indices(this->d_bose,this->d_bose));
    mwArray Theta(Indices(this->d_bose,this->d_bose));
    mwArray L_operator(Indices(this->d_bose,this->d_bose));
    mwArray Lsq_operator(Indices(this->d_bose,this->d_bose));
    mwArray id_bosonic = identityMatrix(this->d_bose);
    
    L_operator.fillWithZero();
    Lsq_operator.fillWithZero();
    Op1.fillWithZero();
    Op2.fillWithZero();
    Theta.fillWithZero();
    
    int offset = (int) floor(this->d_bose/2.0);  
    cout << "offset = " << offset << endl;  
  
    for(int i=0; i<this->d_bose; ++i){
      for(int j=0; j<this->d_bose; ++j){
	if(i==j)
	{
	  L_operator.setElement((double) (i-offset), 0.0,Indices(i,j));
	}
      }
    }
    
    //Construct the theta Operator
    double xval,yval,k_val;
    complex_t tmp;
    for(int i=0; i<this->d_bose; ++i)
    {
      xval = -1.0*(this->d_bose-1)/2.0 + (double) i;
      for(int j=0; j<this->d_bose; ++j)
      {
	yval = -1.0*(this->d_bose-1)/2.0 + (double) j;
	tmp=ZERO_c;
	for(int k=0; k<this->d_bose; ++k)
	{
	  k_val= -1.0*(this->d_bose-1)/2.0 + (double) k;
	  //cout << "x=" << x << " y=" << y << " k=" << k_val << endl;
	  tmp = tmp + 2.0*M_PI/((double) this->d_bose)*k_val*exp(-2.0*I_c*M_PI/((double) this->d_bose)*k_val*(xval-yval));
	}
	Theta.setElement(tmp/((double) this->d_bose),Indices(i,j));
      }
    }
    //Get the gauge invariant version of the operators  
    wrapper::expm(Theta,Op1,I_c);
    wrapper::expm(Theta,Op2,-1.0*I_c);
    
    //I have problems with small detrimental entries in the operators Op1 and Op2 because I calculate them numerically. Since they have a rather clear structure, I try to clean them here in order to avoid numerical instabilities
    complex_t entry;
    double CLEAN_TOL=1.0E-14;
    cout << "WARNING in SHGIGaussLawAngular::constructEvolutionMPOefficient(MPO &Uodd, MPO &Ueven, double dt): I will clean the two operators Op1 and Op2 with tolerance " << CLEAN_TOL << endl;
    
    for(int i=0; i<this->d_bose; i++)
    {
      for(int j=0; j<this->d_bose; j++)
      {
	entry = Op1.getElement(Indices(i,j));
	if(fabs(entry.re)<CLEAN_TOL)
	  entry.re = 0.0;
	if(fabs(entry.im)<CLEAN_TOL)
	  entry.im = 0.0;
	Op1.setElement(entry,Indices(i,j));
	entry = Op2.getElement(Indices(i,j));
      if(fabs(entry.re)<CLEAN_TOL)
	  entry.re = 0.0;
	if(fabs(entry.im)<CLEAN_TOL)
	  entry.im = 0.0;
	Op2.setElement(entry,Indices(i,j));
      }
    }
	
    
    //cout << "L = " << L_operator << endl;
    //cout << "Op1 = " << Op1 << endl;
    //cout << "Op2 = " << Op2 << endl;
    
    //Now get get the square of the L-operator
    Lsq_operator = L_operator*L_operator;
  
    //   cout << "Simga plus: " << sigma_plus << endl;
    //   cout << "Simga minus: " << sigma_minus << endl;
    //   cout << "I + s3: " << idpsz << endl;
    //   cout << "Op1: " << Op1 << endl;
    //   cout << "Op2: " << Op2 << endl;
    //   cout << "id_f: " << id_fermionic << endl;
    //   cout << "id_b: " << id_bosonic << endl;
    //   cout << "Lsq: " << Lsq_operator << endl;
    //   
    //   cout << "kron(sp,sm)=" << kron(sigma_plus,sigma_minus) << endl;
    //   cout << "kron(sm,sp)=" << kron(sigma_minus,sigma_plus) << endl;
    
    
    //Since I want to split id+s3 on every site in 0.5*(id+s3) + 0.5*(id+s3) and take one term to the even MPO and one to the odd MPO I get a problem at the first site and the last site, because one term is missing. To avoid this I have to calculate an extra MPO for the first site and an extra MPO for the site before the last site. Graphically the problem looks like this
    //x--o--x--o--x--o--x
    //|  |  |  |  |  |  |
    //-------  |  -------
    //|  |  -------  |  |
    //Note: That this scheme works, the chain has to have an even number of spin sites, otherwise it could happen, that this hast to be done for an even side MPO
    
    cout << "Generating terms and storing them" << endl;
    
    /*cout << endl << "RAW Terms" << endl;
    cout << "sigma_plus " << sigma_plus << endl;
    cout << "sigma_minus " << sigma_minus << endl;
    cout << "I+sz " << idpsz << endl;
    cout << "Op1 " << Op1 << endl;
    cout << "Op2 " << Op2 << endl << endl;*/
    
    
    //sigma_plus Op1 sigma_minus
    this->hopping1=kron(sigma_plus,kron(Op1,sigma_minus));
    //sigma_minus Op2 sigma_plus
    this->hopping2=kron(sigma_minus,kron(Op2,sigma_plus));
    //mwArray id L^2  id
    this->lterm=kron(id_fermionic,kron(Lsq_operator,id_fermionic));
    // 0.5*idpsz id id
    this->idpszleft=kron(0.5*idpsz,kron(id_bosonic,id_fermionic));
    // id id 0.5*idpsz
    this->idpszright=kron(id_fermionic,kron(id_bosonic,0.5*idpsz));;
    // idpsz id id
    this->firstidpszleft=kron(-1.0*idpsz,kron(id_bosonic,id_fermionic));;  
    // id id idpsz
    this->lastidpszright=kron(id_fermionic,kron(id_bosonic,idpsz));
    
    this->id_bose_mpo=id_bosonic;
    this->id_fermi_mpo=id_fermionic;
    this->id_bose_mpo.reshape(Indices(this->d_bose,1,this->d_bose,1));
    this->id_fermi_mpo.reshape(Indices(this->d_fermi,1,this->d_fermi,1));    
    
    //Terms are now saved
    this->terms_saved=true;
  }
  
  if((this->N_fermions%2)!=0)
    cout << "!!! WARINING, TIME EVOLUTION MPO WILL BE WRONG FOR AN ODD NUMBER OF FERMIONS !!!" << endl;
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Calculate the tensorproducts for the exponentials, for FIRST ODD
  mwArray arg_odd,exp_odd;
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_odd = this->x*hopping1 + this->x*hopping2 + 0.5*this->mu*firstidpszleft + 0.5*this->mu*idpszright + lterm;  
  
  //If it is not hermitian enough I have to ensure that by hand
  if(!arg_odd.isHermitian(1.0E-14))
  {
    mwArray temp;
    temp = arg_odd;
    cout << "Warning arg_odd is not hermitian, will clean it " << endl;
    temp.Hconjugate();
    arg_odd = 0.5*(arg_odd+temp);
  }
      
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  //cout << "Argument exp-odd " << arg_odd << endl;
  //cout << "Exp-odd mit dt "<< dt << ": " << exp_odd << endl;
  
  //Three matrices for the odd MPO
  mwArray W1_odd_first_site,W2_odd_first_site,W3_odd_first_site;
 
  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  mwArray U,S,V;
  wrapper::svd(exp_odd,U,S,V);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_odd_first_site = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  Indices dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  mwArray U2,S2,V2;
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_odd_first_site = U2;
  //Reshape V2 and get the last matrix of the MPO
  Indices dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_odd_first_site = V2;

  //Calculate the tensorproducts for the exponentials, for the site BEFORE THE LAST SITE
  arg_odd.clear(); exp_odd.clear();
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_odd = this->x*hopping1 + this->x*hopping2 - 0.5*this->mu*idpszleft + 0.5*this->mu*lastidpszright + lterm;
  
  //If it is not hermitian enough I have to ensure that by hand
  if(!arg_odd.isHermitian(1.0E-14))
  {
    mwArray temp;
    temp = arg_odd;
    cout << "Warning arg_odd is not hermitian, will clean it " << endl;
    temp.Hconjugate();
    arg_odd = 0.5*(arg_odd+temp);
  }
  
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  //cout << "Argument exp-odd " << arg_odd << endl;
  //cout << "Exp-odd mit dt "<< dt << ": " << exp_odd << endl;
  
  //Three matrices for the odd MPO
  mwArray W1_odd_last_site,W2_odd_last_site,W3_odd_last_site;
 
  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  U.clear(); S.clear(); V.clear();
  wrapper::svd(exp_odd,U,S,V);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_odd_last_site = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_odd_last_site = U2;
  //Reshape V2 and get the last matrix of the MPO
  dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_odd_last_site = V2;    
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  //Calculate the tensorproducts for the exponentials, now for the rest of the odd part
  arg_odd.clear(); exp_odd.clear();
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_odd = this->x*hopping1 + this->x*hopping2 - 0.5*this->mu*idpszleft + 0.5*this->mu*idpszright + lterm;
  
  //If it is not hermitian enough i have to ensure that by hand
  if(!arg_odd.isHermitian(1.0E-14))
  {
    mwArray temp;
    temp = arg_odd;
    cout << "Warning arg_odd is not hermitian, will clean it " << endl;
    temp.Hconjugate();
    arg_odd = 0.5*(arg_odd+temp);
  }
  
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  //cout << "Argument exp-odd " << arg_odd << endl;
  //cout << "Exp-odd mit dt "<< dt << ": " << exp_odd << endl;
  
  //Three matrices for the odd MPO
  mwArray W1_odd,W2_odd,W3_odd;
 
  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  U.clear(); S.clear(); V.clear();
  wrapper::svd(exp_odd,U,S,V);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_odd = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_odd = U2;
  //Reshape V2 and get the last matrix of the MPO
  dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_odd = V2;  
  
   
  //Three matrices for the even MPO
  mwArray W1_even,W2_even,W3_even;
  
  //Calculate the tensorproducts for the exponentials, now for even part
  mwArray arg_even,exp_even;
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_even = this->x*hopping1 + this->x*hopping2 + 0.5*this->mu*idpszleft - 0.5*this->mu*idpszright + lterm;
  
  //If it is not hermitian enough I have to ensure that by hand
  if(!arg_even.isHermitian(1.0E-14))
  {
    mwArray temp;
    temp = arg_even;
    cout << "Warning arg_even is not hermitian, will clean it " << endl;
    temp.Hconjugate();
    arg_even = 0.5*(arg_even+temp);
  }
  
  wrapper::expm(arg_even,exp_even,-I_c*dt);
  //cout << "Argument exp-even " << arg_even << endl;
  //cout << "Exp-even mit dt "<< dt << ": " << exp_even << endl;
  

  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_even.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_even.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_even.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  U.clear(); S.clear(); V.clear();
  wrapper::svd(exp_even,U,S,V);
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_even = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_even = U2;
  //Reshape V2 and get the last matrix of the MPO
  dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_even = V2;  
  
  
  //Now fill the odd MPO
  Uodd.clear();
  Uodd.initLength(this->N_total);
  
  bool first=true;
  
  for(int i=0; i<N_total-2; i +=4)
  {
    //cout << "i=" << i << endl;
    if(i==0)
    {
      //First site of the chain (which is equal to the first fermionic site
      //cout << "Odd MPO at first site" << endl;
      Uodd.setOp(i,new Operator(W1_odd_first_site),true);
      Uodd.setOp(i+1,new Operator(W2_odd_first_site),true);
      Uodd.setOp(i+2,new Operator(W3_odd_first_site),true);
    }
    else if(i==(N_total-3))
    {
      //Site before the last fermionic site
      //cout << "Odd MPO at site before last site" << endl;
      Uodd.setOp(i,new Operator(W1_odd_last_site),true);
      Uodd.setOp(i+1,new Operator(W2_odd_last_site),true);
      Uodd.setOp(i+2,new Operator(W3_odd_last_site),true);
    }      
    else if(first && (i!=0) && i!=(N_total-3))
    {
      //cout << "First odd MPO in between" << endl;
      Uodd.setOp(i,new Operator(W1_odd),true);
      Uodd.setOp(i+1,new Operator(W2_odd),true);
      Uodd.setOp(i+2,new Operator(W3_odd),true);
      first = false;
    }
    else
    {
      //cout << "Odd MPO in between" << endl;
      Uodd.setOp(i,&Uodd.getOp(i-4));
      Uodd.setOp(i+1,&Uodd.getOp(i-3));
      Uodd.setOp(i+2,&Uodd.getOp(i-2));
    }
  }
  //Set Identities
  first = true;
  for(int i=3; i<N_total; i+=4)
  {
    //cout << "i=" << i << endl;
    if(first)
      Uodd.setOp(i,new Operator(id_bose_mpo),true);
    else
      Uodd.setOp(i,&Uodd.getOp(i-4));
  }
  //cout <<"Uodd " << Uodd << endl;
  //Uodd.exportForMatlab("Uodd.m");
  
  /*MPO dummyMPO(3);
  
  dummyMPO.setOp(0,new Operator(W1_odd),false); 
  dummyMPO.setOp(1,new Operator(W2_odd),false); 
  dummyMPO.setOp(2,new Operator(W3_odd),false); 
  
  dummyMPO.exportForMatlab("Uoddsmall.m");*/
    
  
  
  //Now fill the even MPO
  Ueven.clear();
  Ueven.initLength(this->N_total);
  
  first = true;
  for(int i=2; i<N_total-2; i +=4)
  {
    //cout << "i=" << i << endl;
    if(first)
    {
      Ueven.setOp(i,new Operator(W1_even),true);
      Ueven.setOp(i+1,new Operator(W2_even),true);
      Ueven.setOp(i+2,new Operator(W3_even),true);
      first = false;
    }
    else
    {
      Ueven.setOp(i,&Ueven.getOp(i-4));
      Ueven.setOp(i+1,&Ueven.getOp(i-3));
      Ueven.setOp(i+2,&Ueven.getOp(i-2));
    }
  }
  //Set bosonic Identities
  first = true;
  for(int i=1; i<N_total; i+=4)
  {
    //cout << "i=" << i << endl;
    if(first)
      Ueven.setOp(i,new Operator(id_bose_mpo),true);
    else
      Ueven.setOp(i,&Ueven.getOp(i-4));
  }
  
  //First and last fermonic identity
  //cout << "Ueven  " << Ueven << endl;
  Ueven.setOp(0,new Operator(id_fermi_mpo),true);
  Ueven.setOp(N_total-1,&Ueven.getOp(0));
  //cout << "Ueven after setting " << Ueven << endl;
  //cout <<"Ueven " << Ueven << endl;
  
  /*dummyMPO.setOp(0,new Operator(W1_even),false); 
  dummyMPO.setOp(1,new Operator(W2_even),false); 
  dummyMPO.setOp(2,new Operator(W3_even),false); 
  
  dummyMPO.exportForMatlab("Uevensmall.m");*/
  
  //cout << "Passt oder so " << endl;
  //Ueven.exportForMatlab("Ueven.m");
  
}//constructEvolutionMPOefficient


void SHGIGaussLawAngular::constructUoddMPO(MPO &Uodd, double dt)
{
  if(!this->terms_saved)
  {
    //Case in which the single terms are not yet built and have to be generated (also I construct only Uodd, I provide the terms for Ueven, since it is very likely that I will use this one afterwards)
    
    //Provide the operators
    //Provide Pauli-Matrices
    mwArray id_fermionic = identityMatrix(2);
    
    complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
    mwArray sigma_z(Indices(2,2),data_z);
    
    complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
    mwArray sigma_minus(Indices(2,2),data_sm);
    
    complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
    mwArray sigma_plus(Indices(2,2),data_sp);
    
    complex_t data_idpsz[] = {2.0*ONE_c,ZERO_c,ZERO_c,ZERO_c};  
    mwArray idpsz(Indices(2,2),data_idpsz); 
    
    //Provide Operators for gauge field
    mwArray Op1(Indices(this->d_bose,this->d_bose));
    mwArray Op2(Indices(this->d_bose,this->d_bose));
    mwArray Theta(Indices(this->d_bose,this->d_bose));
    mwArray L_operator(Indices(this->d_bose,this->d_bose));
    mwArray Lsq_operator(Indices(this->d_bose,this->d_bose));
    mwArray id_bosonic = identityMatrix(this->d_bose);
    
    L_operator.fillWithZero();
    Lsq_operator.fillWithZero();
    Op1.fillWithZero();
    Op2.fillWithZero();
    Theta.fillWithZero();
    
    int offset = (int) floor(this->d_bose/2.0);  
    cout << "offset = " << offset << endl;
  
  
  
    for(int i=0; i<this->d_bose; ++i){
      for(int j=0; j<this->d_bose; ++j){
	if(i==j)
	{
	  L_operator.setElement((double) (i-offset), 0.0,Indices(i,j));
	}
      }
    }
  
  //Construct the theta Operator
  double xval,yval,k_val;
  complex_t tmp;
  for(int i=0; i<this->d_bose; ++i)
  {
    xval = -1.0*(this->d_bose-1)/2.0 + (double) i;
    for(int j=0; j<this->d_bose; ++j)
    {
      yval = -1.0*(this->d_bose-1)/2.0 + (double) j;
      tmp=ZERO_c;
      for(int k=0; k<this->d_bose; ++k)
      {
	k_val= -1.0*(this->d_bose-1)/2.0 + (double) k;
	//cout << "x=" << x << " y=" << y << " k=" << k_val << endl;
	tmp = tmp + 2.0*M_PI/((double) this->d_bose)*k_val*exp(-2.0*I_c*M_PI/((double) this->d_bose)*k_val*(xval-yval));
      }
      Theta.setElement(tmp/((double) this->d_bose),Indices(i,j));
    }
  }
  //Get the gauge invariant version of the operators  
  wrapper::expm(Theta,Op1,I_c);
  wrapper::expm(Theta,Op2,-1.0*I_c);
  
  //cout << "L = " << L_operator << endl;
  //cout << "Op1 = " << Op1 << endl;
  //cout << "Op2 = " << Op2 << endl;
  
  //Now get get the square of the L-operator
  Lsq_operator = L_operator*L_operator;
  
    //   cout << "Simga plus: " << sigma_plus << endl;
    //   cout << "Simga minus: " << sigma_minus << endl;
    //   cout << "I + s3: " << idpsz << endl;
    //   cout << "Op1: " << Op1 << endl;
    //   cout << "Op2: " << Op2 << endl;
    //   cout << "id_f: " << id_fermionic << endl;
    //   cout << "id_b: " << id_bosonic << endl;
    //   cout << "Lsq: " << Lsq_operator << endl;
    //   
    //   cout << "kron(sp,sm)=" << kron(sigma_plus,sigma_minus) << endl;
    //   cout << "kron(sm,sp)=" << kron(sigma_minus,sigma_plus) << endl;
    
    
    //Since I want to split id+s3 on every site in 0.5*(id+s3) + 0.5*(id+s3) and take one term to the even MPO and one to the odd MPO I get a problem at the first site and the last site, because one term is missing. To avoid this I have to calculate an extra MPO for the first site and an extra MPO for the site before the last site. Graphically the problem looks like this
    //x--o--x--o--x--o--x
    //|  |  |  |  |  |  |
    //-------  |  -------
    //|  |  -------  |  |
    //Note: That this scheme works, the chain has to have an even number of spin sites, otherwise it could happen, that this hast to be done for an even side MPO
    
    cout << "Generating terms and storing them" << endl;
    
    //sigma_plus Op1 sigma_minus
    this->hopping1=kron(sigma_plus,kron(Op1,sigma_minus));
    //sigma_minus Op2 sigma_plus
    this->hopping2=kron(sigma_minus,kron(Op2,sigma_plus));
    //mwArray id L^2  id
    this->lterm=kron(id_fermionic,kron(Lsq_operator,id_fermionic));
    // 0.5*idpsz id id
    this->idpszleft=kron(0.5*idpsz,kron(id_bosonic,id_fermionic));
    // id id 0.5*idpsz
    this->idpszright=kron(id_fermionic,kron(id_bosonic,0.5*idpsz));;
    // idpsz id id
    this->firstidpszleft=kron(-1.0*idpsz,kron(id_bosonic,id_fermionic));;  
    // id id idpsz
    this->lastidpszright=kron(id_fermionic,kron(id_bosonic,idpsz));
    
    this->id_bose_mpo=id_bosonic;
    this->id_fermi_mpo=id_fermionic;
    this->id_bose_mpo.reshape(Indices(this->d_bose,1,this->d_bose,1));
    this->id_fermi_mpo.reshape(Indices(this->d_fermi,1,this->d_fermi,1));    
    
    //Terms are now saved
    this->terms_saved=true;
  }
  
  if((this->N_fermions%2)!=0)
    cout << "!!! WARINING, TIME EVOLUTION MPO WILL BE WRONG FOR AN ODD NUMBER OF FERMIONS !!!" << endl;
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Calculate the tensorproducts for the exponentials, for FIRST ODD
  mwArray arg_odd,exp_odd;
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_odd = this->x*hopping1 + this->x*hopping2 + 0.5*this->mu*firstidpszleft + 0.5*this->mu*idpszright + lterm;
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  //cout << "Argument exp-odd " << arg_odd << endl;
  //cout << "Exp-odd mit dt "<< dt << ": " << exp_odd << endl;
  
  //Three matrices for the odd MPO
  mwArray W1_odd_first_site,W2_odd_first_site,W3_odd_first_site;
 
  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  mwArray U,S,V;
  wrapper::svd(exp_odd,U,S,V);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_odd_first_site = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  Indices dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  mwArray U2,S2,V2;
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_odd_first_site = U2;
  //Reshape V2 and get the last matrix of the MPO
  Indices dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_odd_first_site = V2;

  //Calculate the tensorproducts for the exponentials, for the site BEFORE THE LAST SITE
  arg_odd.clear(); exp_odd.clear();
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_odd = this->x*hopping1 + this->x*hopping2 - 0.5*this->mu*idpszleft + 0.5*this->mu*lastidpszright + lterm;
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  //cout << "Argument exp-odd " << arg_odd << endl;
  //cout << "Exp-odd mit dt "<< dt << ": " << exp_odd << endl;
  
  //Three matrices for the odd MPO
  mwArray W1_odd_last_site,W2_odd_last_site,W3_odd_last_site;
 
  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  U.clear(); S.clear(); V.clear();
  wrapper::svd(exp_odd,U,S,V);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_odd_last_site = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_odd_last_site = U2;
  //Reshape V2 and get the last matrix of the MPO
  dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_odd_last_site = V2;    
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  //Calculate the tensorproducts for the exponentials, now for the rest of the odd part
  arg_odd.clear(); exp_odd.clear();
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_odd = this->x*hopping1 + this->x*hopping2 - 0.5*this->mu*idpszleft + 0.5*this->mu*idpszright + lterm;
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  //cout << "Argument exp-odd " << arg_odd << endl;
  //cout << "Exp-odd mit dt "<< dt << ": " << exp_odd << endl;
  
  //Three matrices for the odd MPO
  mwArray W1_odd,W2_odd,W3_odd;
 
  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  U.clear(); S.clear(); V.clear();
  wrapper::svd(exp_odd,U,S,V);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_odd = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_odd = U2;
  //Reshape V2 and get the last matrix of the MPO
  dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_odd = V2;  
  
  
  //Now fill the odd MPO
  Uodd.clear();
  Uodd.initLength(this->N_total);
  
  bool first=true;
  
  for(int i=0; i<N_total-2; i +=4)
  {
    //cout << "i=" << i << endl;
    if(i==0)
    {
      //First site of the chain (which is equal to the first fermionic site
      //cout << "Odd MPO at first site" << endl;
      Uodd.setOp(i,new Operator(W1_odd_first_site),true);
      Uodd.setOp(i+1,new Operator(W2_odd_first_site),true);
      Uodd.setOp(i+2,new Operator(W3_odd_first_site),true);
    }
    else if(i==(N_total-3))
    {
      //Site before the last fermionic site
      //cout << "Odd MPO at site before last site" << endl;
      Uodd.setOp(i,new Operator(W1_odd_last_site),true);
      Uodd.setOp(i+1,new Operator(W2_odd_last_site),true);
      Uodd.setOp(i+2,new Operator(W3_odd_last_site),true);
    }      
    else if(first && (i!=0) && i!=(N_total-3))
    {
      //cout << "First odd MPO in between" << endl;
      Uodd.setOp(i,new Operator(W1_odd),true);
      Uodd.setOp(i+1,new Operator(W2_odd),true);
      Uodd.setOp(i+2,new Operator(W3_odd),true);
      first = false;
    }
    else
    {
      //cout << "Odd MPO in between" << endl;
      Uodd.setOp(i,&Uodd.getOp(i-4));
      Uodd.setOp(i+1,&Uodd.getOp(i-3));
      Uodd.setOp(i+2,&Uodd.getOp(i-2));
    }
  }
  //Set Identities
  first = true;
  for(int i=3; i<N_total; i+=4)
  {
    //cout << "i=" << i << endl;
    if(first)
      Uodd.setOp(i,new Operator(id_bose_mpo),true);
    else
      Uodd.setOp(i,&Uodd.getOp(i-4));
  }
  //cout <<"Uodd " << Uodd << endl;
  //Uodd.exportForMatlab("Uodd.m");
  
  /*MPO dummyMPO(3);
  
  dummyMPO.setOp(0,new Operator(W1_odd),false); 
  dummyMPO.setOp(1,new Operator(W2_odd),false); 
  dummyMPO.setOp(2,new Operator(W3_odd),false); 
  
  dummyMPO.exportForMatlab("Uoddsmall.m");*/  
}//constructUoddMPO

void SHGIGaussLawAngular::constructUoddMPO(MPO &Uodd, double dt, double x1, double x2)
{
  if(!this->terms_saved)
  {
    //Case in which the single terms are not yet built and have to be generated (also I construct only Uodd, I provide the terms for Ueven, since it is very likely that I will use this one afterwards)
    
    //Provide the operators
    //Provide Pauli-Matrices
    mwArray id_fermionic = identityMatrix(2);
    
    complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
    mwArray sigma_z(Indices(2,2),data_z);
    
    complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
    mwArray sigma_minus(Indices(2,2),data_sm);
    
    complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
    mwArray sigma_plus(Indices(2,2),data_sp);
    
    complex_t data_idpsz[] = {2.0*ONE_c,ZERO_c,ZERO_c,ZERO_c};  
    mwArray idpsz(Indices(2,2),data_idpsz); 
    
    //Provide Operators for gauge field
    mwArray Op1(Indices(this->d_bose,this->d_bose));
    mwArray Op2(Indices(this->d_bose,this->d_bose));
    mwArray Theta(Indices(this->d_bose,this->d_bose));
    mwArray L_operator(Indices(this->d_bose,this->d_bose));
    mwArray Lsq_operator(Indices(this->d_bose,this->d_bose));
    mwArray id_bosonic = identityMatrix(this->d_bose);
    
    L_operator.fillWithZero();
    Lsq_operator.fillWithZero();
    Op1.fillWithZero();
    Op2.fillWithZero();
    Theta.fillWithZero(); 
    
    int offset = (int) floor(this->d_bose/2.0);  
    cout << "offset = " << offset << endl; 
  
    for(int i=0; i<this->d_bose; ++i){
      for(int j=0; j<this->d_bose; ++j){
	if(i==j)
	{
	  L_operator.setElement((double) (i-offset), 0.0,Indices(i,j));
	}
      }
    }
  
    //Construct the theta Operator
    double xval,yval,k_val;
    complex_t tmp;
    for(int i=0; i<this->d_bose; ++i)
    {
      xval = -1.0*(this->d_bose-1)/2.0 + (double) i;
      for(int j=0; j<this->d_bose; ++j)
      {
	yval = -1.0*(this->d_bose-1)/2.0 + (double) j;
	tmp=ZERO_c;
	for(int k=0; k<this->d_bose; ++k)
	{
	  k_val= -1.0*(this->d_bose-1)/2.0 + (double) k;
	  //cout << "x=" << x << " y=" << y << " k=" << k_val << endl;
	  tmp = tmp + 2.0*M_PI/((double) this->d_bose)*k_val*exp(-2.0*I_c*M_PI/((double) this->d_bose)*k_val*(xval-yval));
	}
	Theta.setElement(tmp/((double) this->d_bose),Indices(i,j));
      }
    }
    //Get the gauge invariant version of the operators  
    wrapper::expm(Theta,Op1,I_c);
    wrapper::expm(Theta,Op2,-1.0*I_c);
    
    //cout << "L = " << L_operator << endl;
    //cout << "Op1 = " << Op1 << endl;
    //cout << "Op2 = " << Op2 << endl;
    
    //Now get get the square of the L-operator
    Lsq_operator = L_operator*L_operator;
  
    //   cout << "Simga plus: " << sigma_plus << endl;
    //   cout << "Simga minus: " << sigma_minus << endl;
    //   cout << "I + s3: " << idpsz << endl;
    //   cout << "Op1: " << Op1 << endl;
    //   cout << "Op2: " << Op2 << endl;
    //   cout << "id_f: " << id_fermionic << endl;
    //   cout << "id_b: " << id_bosonic << endl;
    //   cout << "Lsq: " << Lsq_operator << endl;
    //   
    //   cout << "kron(sp,sm)=" << kron(sigma_plus,sigma_minus) << endl;
    //   cout << "kron(sm,sp)=" << kron(sigma_minus,sigma_plus) << endl;
    
    
    //Since I want to split id+s3 on every site in 0.5*(id+s3) + 0.5*(id+s3) and take one term to the even MPO and one to the odd MPO I get a problem at the first site and the last site, because one term is missing. To avoid this I have to calculate an extra MPO for the first site and an extra MPO for the site before the last site. Graphically the problem looks like this
    //x--o--x--o--x--o--x
    //|  |  |  |  |  |  |
    //-------  |  -------
    //|  |  -------  |  |
    //Note: That this scheme works, the chain has to have an even number of spin sites, otherwise it could happen, that this hast to be done for an even side MPO
    
    cout << "Generating terms and storing them" << endl;
    
    //sigma_plus Op1 sigma_minus
    this->hopping1=kron(sigma_plus,kron(Op1,sigma_minus));
    //sigma_minus Op2 sigma_plus
    this->hopping2=kron(sigma_minus,kron(Op2,sigma_plus));
    //mwArray id L^2  id
    this->lterm=kron(id_fermionic,kron(Lsq_operator,id_fermionic));
    // 0.5*idpsz id id
    this->idpszleft=kron(0.5*idpsz,kron(id_bosonic,id_fermionic));
    // id id 0.5*idpsz
    this->idpszright=kron(id_fermionic,kron(id_bosonic,0.5*idpsz));;
    // idpsz id id
    this->firstidpszleft=kron(-1.0*idpsz,kron(id_bosonic,id_fermionic));;  
    // id id idpsz
    this->lastidpszright=kron(id_fermionic,kron(id_bosonic,idpsz));
    
    this->id_bose_mpo=id_bosonic;
    this->id_fermi_mpo=id_fermionic;
    this->id_bose_mpo.reshape(Indices(this->d_bose,1,this->d_bose,1));
    this->id_fermi_mpo.reshape(Indices(this->d_fermi,1,this->d_fermi,1));    
    
    //Terms are now saved
    this->terms_saved=true;
  }
  
  if((this->N_fermions%2)!=0)
    cout << "!!! WARINING, TIME EVOLUTION MPO WILL BE WRONG FOR AN ODD NUMBER OF FERMIONS !!!" << endl;
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Calculate the tensorproducts for the exponentials, for FIRST ODD
  mwArray arg_odd,exp_odd;
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_odd = 0.5*(x1+x2)*hopping1 + 0.5*(x1+x2)*hopping2 + 0.5*this->mu*firstidpszleft + 0.5*this->mu*idpszright + lterm;
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  //cout << "Argument exp-odd " << arg_odd << endl;
  //cout << "Exp-odd mit dt "<< dt << ": " << exp_odd << endl;
  
  //Three matrices for the odd MPO
  mwArray W1_odd_first_site,W2_odd_first_site,W3_odd_first_site;
 
  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  mwArray U,S,V;
  wrapper::svd(exp_odd,U,S,V);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_odd_first_site = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  Indices dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  mwArray U2,S2,V2;
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_odd_first_site = U2;
  //Reshape V2 and get the last matrix of the MPO
  Indices dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_odd_first_site = V2;

  //Calculate the tensorproducts for the exponentials, for the site BEFORE THE LAST SITE
  arg_odd.clear(); exp_odd.clear();
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_odd = 0.5*(x1+x2)*hopping1 + 0.5*(x1+x2)*hopping2 - 0.5*this->mu*idpszleft + 0.5*this->mu*lastidpszright + lterm;
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  //cout << "Argument exp-odd " << arg_odd << endl;
  //cout << "Exp-odd mit dt "<< dt << ": " << exp_odd << endl;
  
  //Three matrices for the odd MPO
  mwArray W1_odd_last_site,W2_odd_last_site,W3_odd_last_site;
 
  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  U.clear(); S.clear(); V.clear();
  wrapper::svd(exp_odd,U,S,V);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_odd_last_site = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_odd_last_site = U2;
  //Reshape V2 and get the last matrix of the MPO
  dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_odd_last_site = V2;    
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  //Calculate the tensorproducts for the exponentials, now for the rest of the odd part
  arg_odd.clear(); exp_odd.clear();
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_odd = 0.5*(x1+x2)*hopping1 + 0.5*(x1+x2)*hopping2 - 0.5*this->mu*idpszleft + 0.5*this->mu*idpszright + lterm;
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  //cout << "Argument exp-odd " << arg_odd << endl;
  //cout << "Exp-odd mit dt "<< dt << ": " << exp_odd << endl;
  
  //Three matrices for the odd MPO
  mwArray W1_odd,W2_odd,W3_odd;
 
  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  U.clear(); S.clear(); V.clear();
  wrapper::svd(exp_odd,U,S,V);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_odd = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_odd = U2;
  //Reshape V2 and get the last matrix of the MPO
  dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_odd = V2;   
  
  //Now fill the odd MPO
  Uodd.clear();
  Uodd.initLength(this->N_total);
  
  bool first=true;
  
  for(int i=0; i<N_total-2; i +=4)
  {
    //cout << "i=" << i << endl;
    if(i==0)
    {
      //First site of the chain (which is equal to the first fermionic site
      //cout << "Odd MPO at first site" << endl;
      Uodd.setOp(i,new Operator(W1_odd_first_site),true);
      Uodd.setOp(i+1,new Operator(W2_odd_first_site),true);
      Uodd.setOp(i+2,new Operator(W3_odd_first_site),true);
    }
    else if(i==(N_total-3))
    {
      //Site before the last fermionic site
      //cout << "Odd MPO at site before last site" << endl;
      Uodd.setOp(i,new Operator(W1_odd_last_site),true);
      Uodd.setOp(i+1,new Operator(W2_odd_last_site),true);
      Uodd.setOp(i+2,new Operator(W3_odd_last_site),true);
    }      
    else if(first && (i!=0) && i!=(N_total-3))
    {
      //cout << "First odd MPO in between" << endl;
      Uodd.setOp(i,new Operator(W1_odd),true);
      Uodd.setOp(i+1,new Operator(W2_odd),true);
      Uodd.setOp(i+2,new Operator(W3_odd),true);
      first = false;
    }
    else
    {
      //cout << "Odd MPO in between" << endl;
      Uodd.setOp(i,&Uodd.getOp(i-4));
      Uodd.setOp(i+1,&Uodd.getOp(i-3));
      Uodd.setOp(i+2,&Uodd.getOp(i-2));
    }
  }
  //Set Identities
  first = true;
  for(int i=3; i<N_total; i+=4)
  {
    //cout << "i=" << i << endl;
    if(first)
      Uodd.setOp(i,new Operator(id_bose_mpo),true);
    else
      Uodd.setOp(i,&Uodd.getOp(i-4));
  }
  //cout <<"Uodd " << Uodd << endl;
  //Uodd.exportForMatlab("Uodd.m");
  
  /*MPO dummyMPO(3);
  
  dummyMPO.setOp(0,new Operator(W1_odd),false); 
  dummyMPO.setOp(1,new Operator(W2_odd),false); 
  dummyMPO.setOp(2,new Operator(W3_odd),false); 
  
  dummyMPO.exportForMatlab("Uoddsmall.m");*/  
}//constructUoddMPO (version with x1 and x2 given)

void SHGIGaussLawAngular::constructUevenMPO(MPO &Ueven, double dt)
{
  if(!this->terms_saved)
  {
    //Case in which the single terms are not yet built and have to be generated
    //Provide the operators
    //Provide Pauli-Matrices
    mwArray id_fermionic = identityMatrix(2);
    
    complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
    mwArray sigma_z(Indices(2,2),data_z);
    
    complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
    mwArray sigma_minus(Indices(2,2),data_sm);
    
    complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
    mwArray sigma_plus(Indices(2,2),data_sp);
    
    complex_t data_idpsz[] = {2.0*ONE_c,ZERO_c,ZERO_c,ZERO_c};  
    mwArray idpsz(Indices(2,2),data_idpsz); 
    
    //Provide Operators for gauge field
    mwArray Op1(Indices(this->d_bose,this->d_bose));
    mwArray Op2(Indices(this->d_bose,this->d_bose));
    mwArray Theta(Indices(this->d_bose,this->d_bose));
    mwArray L_operator(Indices(this->d_bose,this->d_bose));
    mwArray Lsq_operator(Indices(this->d_bose,this->d_bose));
    mwArray id_bosonic = identityMatrix(this->d_bose);
    
    L_operator.fillWithZero();
    Lsq_operator.fillWithZero();
    Op1.fillWithZero();
    Op2.fillWithZero();
    Theta.fillWithZero();
  
    int offset = (int) floor(this->d_bose/2.0);  
    cout << "offset = " << offset << endl;   
    
    for(int i=0; i<this->d_bose; ++i){
      for(int j=0; j<this->d_bose; ++j){
	if(i==j)
	{
	  L_operator.setElement((double) (i-offset), 0.0,Indices(i,j));
	}
      }
    }
    
    //Construct the theta Operator
    double xval,yval,k_val;
    complex_t tmp;
    for(int i=0; i<this->d_bose; ++i)
    {
      xval = -1.0*(this->d_bose-1)/2.0 + (double) i;
      for(int j=0; j<this->d_bose; ++j)
      {
	yval = -1.0*(this->d_bose-1)/2.0 + (double) j;
	tmp=ZERO_c;
	for(int k=0; k<this->d_bose; ++k)
	{
	  k_val= -1.0*(this->d_bose-1)/2.0 + (double) k;
	  //cout << "x=" << x << " y=" << y << " k=" << k_val << endl;
	  tmp = tmp + 2.0*M_PI/((double) this->d_bose)*k_val*exp(-2.0*I_c*M_PI/((double) this->d_bose)*k_val*(xval-yval));
	}
	Theta.setElement(tmp/((double) this->d_bose),Indices(i,j));
      }
    }
    //Get the gauge invariant version of the operators  
    wrapper::expm(Theta,Op1,I_c);
    wrapper::expm(Theta,Op2,-1.0*I_c);
    
    //cout << "L = " << L_operator << endl;
    //cout << "Op1 = " << Op1 << endl;
    //cout << "Op2 = " << Op2 << endl;
    
    //Now get get the square of the L-operator
    Lsq_operator = L_operator*L_operator;

    //   cout << "Simga plus: " << sigma_plus << endl;
    //   cout << "Simga minus: " << sigma_minus << endl;
    //   cout << "I + s3: " << idpsz << endl;
    //   cout << "Op1: " << Op1 << endl;
    //   cout << "Op2: " << Op2 << endl;
    //   cout << "id_f: " << id_fermionic << endl;
    //   cout << "id_b: " << id_bosonic << endl;
    //   cout << "Lsq: " << Lsq_operator << endl;
    //   
    //   cout << "kron(sp,sm)=" << kron(sigma_plus,sigma_minus) << endl;
    //   cout << "kron(sm,sp)=" << kron(sigma_minus,sigma_plus) << endl;
    
    
    //Since I want to split id+s3 on every site in 0.5*(id+s3) + 0.5*(id+s3) and take one term to the even MPO and one to the odd MPO I get a problem at the first site and the last site, because one term is missing. To avoid this I have to calculate an extra MPO for the first site and an extra MPO for the site before the last site. Graphically the problem looks like this
    //x--o--x--o--x--o--x
    //|  |  |  |  |  |  |
    //-------  |  -------
    //|  |  -------  |  |
    //Note: That this scheme works, the chain has to have an even number of spin sites, otherwise it could happen, that this hast to be done for an even side MPO
    
    cout << "Generating terms and storing them" << endl;
    
    //sigma_plus Op1 sigma_minus
    this->hopping1=kron(sigma_plus,kron(Op1,sigma_minus));
    //sigma_minus Op2 sigma_plus
    this->hopping2=kron(sigma_minus,kron(Op2,sigma_plus));
    //mwArray id L^2  id
    this->lterm=kron(id_fermionic,kron(Lsq_operator,id_fermionic));
    // 0.5*idpsz id id
    this->idpszleft=kron(0.5*idpsz,kron(id_bosonic,id_fermionic));
    // id id 0.5*idpsz
    this->idpszright=kron(id_fermionic,kron(id_bosonic,0.5*idpsz));;
    // idpsz id id
    this->firstidpszleft=kron(-1.0*idpsz,kron(id_bosonic,id_fermionic));;  
    // id id idpsz
    this->lastidpszright=kron(id_fermionic,kron(id_bosonic,idpsz));
    
    this->id_bose_mpo=id_bosonic;
    this->id_fermi_mpo=id_fermionic;
    this->id_bose_mpo.reshape(Indices(this->d_bose,1,this->d_bose,1));
    this->id_fermi_mpo.reshape(Indices(this->d_fermi,1,this->d_fermi,1));    
    
    //Terms are now saved
    this->terms_saved=true;
  }
  
  if((this->N_fermions%2)!=0)
    cout << "!!! WARINING, TIME EVOLUTION MPO WILL BE WRONG FOR AN ODD NUMBER OF FERMIONS !!!" << endl;
  
  
  //Three matrices for the even MPO
  mwArray W1_even,W2_even,W3_even;
  
  //Calculate the tensorproducts for the exponentials, now for even part
  mwArray arg_even,exp_even;
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_even = this->x*hopping1 + this->x*hopping2 + 0.5*this->mu*idpszleft - 0.5*this->mu*idpszright + lterm;
  wrapper::expm(arg_even,exp_even,-I_c*dt);
  //cout << "Argument exp-even " << arg_even << endl;
  //cout << "Exp-even mit dt "<< dt << ": " << exp_even << endl;
  

  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_even.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_fermi,this->d_bose,this->d_fermi));
  exp_even.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_even.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose*this->d_fermi*this->d_fermi));
  mwArray U,S,V;
  wrapper::svd(exp_even,U,S,V);
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_even = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  Indices dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  mwArray U2,S2,V2;
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose,this->d_bose,-1));
  U2.permute(Indices(2,1,3,4));
  W2_even = U2;
  //Reshape V2 and get the last matrix of the MPO
  Indices dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V2.permute(Indices(2,1,3,4));
  W3_even = V2;   
  
  //Now fill the even MPO
  Ueven.clear();
  Ueven.initLength(this->N_total);
  
  bool first = true;
  for(int i=2; i<N_total-2; i +=4)
  {
    //cout << "i=" << i << endl;
    if(first)
    {
      Ueven.setOp(i,new Operator(W1_even),true);
      Ueven.setOp(i+1,new Operator(W2_even),true);
      Ueven.setOp(i+2,new Operator(W3_even),true);
      first = false;
    }
    else
    {
      Ueven.setOp(i,&Ueven.getOp(i-4));
      Ueven.setOp(i+1,&Ueven.getOp(i-3));
      Ueven.setOp(i+2,&Ueven.getOp(i-2));
    }
  }
  //Set bosonic Identities
  first = true;
  for(int i=1; i<N_total; i+=4)
  {
    //cout << "i=" << i << endl;
    if(first)
      Ueven.setOp(i,new Operator(id_bose_mpo),true);
    else
      Ueven.setOp(i,&Ueven.getOp(i-4));
  }
  
  //First and last fermonic identity
  //cout << "Ueven  " << Ueven << endl;
  Ueven.setOp(0,new Operator(id_fermi_mpo),true);
  Ueven.setOp(N_total-1,&Ueven.getOp(0));
  //cout << "Ueven after setting " << Ueven << endl;
  //cout <<"Ueven " << Ueven << endl;
  
  /*dummyMPO.setOp(0,new Operator(W1_even),false); 
  dummyMPO.setOp(1,new Operator(W2_even),false); 
  dummyMPO.setOp(2,new Operator(W3_even),false); 
  
  dummyMPO.exportForMatlab("Uevensmall.m");*/
  
  //cout << "Passt oder so " << endl;
  //Ueven.exportForMatlab("Ueven.m");
  
}//constructUevenMPO

void SHGIGaussLawAngular::constructUnitarySubspaceMPO(MPO &Uodd, MPO &Ueven, double dt, double strength)
{
  
  //First get the operators
  mwArray Uspin(Indices(this->d_fermi,this->d_fermi));
  mwArray Uspindagger(Indices(this->d_fermi,this->d_fermi));
  mwArray Uleft(Indices(this->d_bose,this->d_bose));
  mwArray Uright(Indices(this->d_bose,this->d_bose));
  mwArray idbose = identityMatrix(this->d_bose);
  mwArray idfermi = identityMatrix(this->d_fermi);
  
  //In case there are not yet saved
  this->id_bose_mpo=idbose;
  this->id_fermi_mpo=idfermi;
  this->id_bose_mpo.reshape(Indices(this->d_bose,1,this->d_bose,1));
  this->id_fermi_mpo.reshape(Indices(this->d_fermi,1,this->d_fermi,1));    
  
  Uspin.fillWithZero();
  Uleft.fillWithZero();
  Uright.fillWithZero();
  
  //Factor which I choose 1 for convenience 
  int m=1;
  
  //Constant for the phase (since only a couple of phases are possible)
  double alpha = 2.0*M_PI/this->d_bose*m;
  //cout << "Phase constant alpha: " << alpha << endl;
  
  Uspin.setElement(exp(-I_c*0.5*alpha),Indices(0,0));
  Uspin.setElement(exp(I_c*0.5*alpha),Indices(1,1));
  
  Uspindagger = Uspin;
  Uspindagger.transpose(true);
  
  //cout << "Uspin " << Uspin << endl;
  //cout << "Uspindagger " << Uspindagger << endl;
  
  int l = (int) ((this->d_bose-1)/2);
  cout << "l=" << l << endl;
  
  for(int i=0; i<this->d_bose; i++)
  {
    Uleft.setElement(exp(-I_c*(-l+i)*alpha),Indices(i,i));
    Uright.setElement(exp(I_c*(-l+i)*alpha),Indices(i,i));
  }
  
  cout << "Uleft " << Uleft << endl;
  //cout << "Uright " << Uright << endl;
  
  //Possible phases
  complex_t phase_up_site = exp(-I_c*alpha/2.0);
  complex_t phase_down_site = exp(I_c*alpha/2.0);
  
  //cout << "Phase up " << phase_up_site << endl;
  //cout << "Phase down " << phase_down_site << endl;
  
  //Various kronecker products appearing
  mwArray UrUsUl, UlUsUr, identities;
  
  UrUsUl = kron(Uright,kron(Uspin,Uleft));
  UlUsUr = kron(Uleft,kron(Uspindagger,Uright));
  identities = kron(idbose,kron(idfermi,idbose));
  
  //Get exponentials (attention: the three body terms start at a bose site, those which start at a even bose site have a odd spin site in between an the other way round)
  //Term starting at odd bosonic site (so have an even spin in between)
  mwArray arg_odd,exp_odd;
  //Term starting at even bosonic site (so have an odd spin in between)
  mwArray arg_even,exp_even;
  
  arg_odd = strength*(2*identities - phase_down_site*UrUsUl - conjugate(phase_down_site)*UlUsUr);  
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  
  arg_even = strength*(2*identities - phase_up_site*UrUsUl -  conjugate(phase_up_site)*UlUsUr);
  wrapper::expm(arg_even,exp_even,-I_c*dt);
  
  
  //Now proceed with SVD    
  
  //Three matrices for the odd MPO
  mwArray W1_odd,W2_odd,W3_odd;
 
  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_odd.reshape(Indices(this->d_bose,this->d_fermi,this->d_bose,this->d_bose,this->d_fermi,this->d_bose));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_bose*this->d_bose,this->d_fermi*this->d_fermi*this->d_bose*this->d_bose));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  mwArray U,S,V;
  wrapper::svd(exp_odd,U,S,V);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_bose,1,this->d_bose,-1));
  W1_odd = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  Indices dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_fermi*this->d_fermi,this->d_bose*this->d_bose));
  
  //No proceed with another SVD
  mwArray U2,S2,V2;
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_fermi,this->d_fermi,-1));
  U2.permute(Indices(2,1,3,4));
  W2_odd = U2;
  //Reshape V2 and get the last matrix of the MPO
  Indices dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_bose,this->d_bose,1));
  V2.permute(Indices(2,1,3,4));
  W3_odd = V2;  
  
   
  //Three matrices for the even MPO
  mwArray W1_even,W2_even,W3_even;

  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_even.reshape(Indices(this->d_bose,this->d_fermi,this->d_bose,this->d_bose,this->d_fermi,this->d_bose));
  exp_even.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_even.reshape(Indices(this->d_bose*this->d_bose,this->d_fermi*this->d_fermi*this->d_bose*this->d_bose));
  U.clear(); S.clear(); V.clear();
  wrapper::svd(exp_even,U,S,V);
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  U.reshape(Indices(this->d_bose,1,this->d_bose,-1));
  W1_even = U;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  dimV=V.getDimensions();
  V.reshape(Indices(dimV[0]*this->d_fermi*this->d_fermi,this->d_bose*this->d_bose));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_fermi,this->d_fermi,-1));
  U2.permute(Indices(2,1,3,4));
  W2_even = U2;
  //Reshape V2 and get the last matrix of the MPO
  dimV2=V2.getDimensions();
  V2.reshape(Indices(-1,this->d_bose,this->d_bose,1));
  V2.permute(Indices(2,1,3,4));
  W3_even = V2;   
  
  //Now get the left edge term concerning the Ueven MPO
  arg_even.clear(); exp_even.clear();
  mwArray W1_even_edge,W2_even_edge;
  
  arg_even = strength*(2.0*kron(idfermi,idbose) - phase_up_site*kron(Uspindagger,Uleft) - conjugate(phase_up_site)*kron(Uspin,Uright));  
  mwArray exp_even_edge;
  wrapper::expm(arg_even,exp_even_edge,-I_c*dt);
  
  //cout << "arg_even for edges " << arg_even << endl;
  
  exp_even_edge.reshape(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_bose));
  exp_even_edge.permute(Indices(1,3,2,4));
  
  //Now take the first two indices together and perform a SVD
  exp_even_edge.reshape(Indices(this->d_fermi*this->d_fermi,this->d_bose*this->d_bose));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  mwArray Uedge,Sedge,Vedge;  
  wrapper::svd(exp_even_edge,Uedge,Sedge,Vedge);
  //Redistribute the singular values  
  Sedge=sqrt(Sedge);
  Uedge.multiplyRight(Sedge);
  Vedge.multiplyLeft(Sedge);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  Uedge.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  W1_even_edge = Uedge;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  Vedge.reshape(Indices(-1,this->d_bose,this->d_bose,1));
  Vedge.permute(Indices(2,1,3,4));
  W2_even_edge = Vedge; 
  
  //Now get the right edge term concerning the Uodd MPO
  arg_odd.clear(); exp_odd.clear();
  mwArray W1_odd_edge,W2_odd_edge;
  
  arg_odd = strength*(2*kron(idbose,idfermi) - phase_down_site*kron(Uright,Uspindagger) - conjugate(phase_down_site)*kron(Uleft,Uspin));  
  wrapper::expm(arg_odd,exp_odd,-I_c*dt);
  
  exp_odd.reshape(Indices(this->d_bose,this->d_fermi,this->d_bose,this->d_fermi));
  exp_odd.permute(Indices(1,3,2,4));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_bose*this->d_bose,this->d_fermi*this->d_fermi));
  //cout << "Exp-odd gereshaped und permutiert " << exp_odd << endl;
  mwArray Uedge2,Sedge2,Vedge2;
  wrapper::svd(exp_odd,Uedge2,Sedge2,Vedge2);
  /*cout << "U=" << U << endl;
  cout << "S=" << S << endl;
  cout << "V=" << V << endl;*/
  //Redistribute the singular values
  Sedge2=sqrt(Sedge2);
  Uedge2.multiplyRight(Sedge2);
  Vedge2.multiplyLeft(Sedge2);
  //Reshape U and get the first matrix of the MPO (add dummy-index 1 on the left side)
  Uedge2.reshape(Indices(this->d_bose,1,this->d_bose,-1));
  W1_odd_edge = Uedge2;
  //Prepare V for the next SVD by taking the next two indices to the row-index introduced by the SVD
  Vedge2.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  Vedge2.permute(Indices(2,1,3,4));

  W2_odd_edge = Vedge2; 
  
  /*cout << "Edgematrices W1_even_edge " << W1_even_edge << endl;
  cout << "Edgematrices W2_even_edge " << W2_even_edge << endl;
  
  cout << "Edgematrices W1_odd_edge " << W1_odd_edge << endl;
  cout << "Edgematrices W2_odd_edge " << W2_odd_edge << endl;*/
  
  
  //Now fill the odd MPO
  Uodd.clear();
  Uodd.initLength(this->N_total);
  
  bool first=true;
  
  for(int i=1; i<N_total-3; i +=4)
  {
   if(first)
    {
      //cout << "First odd MPO in between" << endl;
      Uodd.setOp(i,new Operator(W1_odd),true);
      Uodd.setOp(i+1,new Operator(W2_odd),true);
      Uodd.setOp(i+2,new Operator(W3_odd),true);
      first = false;
    }
    else
    {
      //cout << "Odd MPO in between" << endl;
      Uodd.setOp(i,&Uodd.getOp(i-4));
      Uodd.setOp(i+1,&Uodd.getOp(i-3));
      Uodd.setOp(i+2,&Uodd.getOp(i-2));
    }
  }
  //Set fermionic Identities
  first = true;
  for(int i=0; i<N_total; i+=4)
  {
    //cout << "i=" << i << endl;
    if(first)
      Uodd.setOp(i,new Operator(id_fermi_mpo),true);
    else
      Uodd.setOp(i,&Uodd.getOp(i-4));
  }
  //Edge terms at the end
  Uodd.setOp(this->N_total-2,new Operator(W1_odd_edge),true);
  Uodd.setOp(this->N_total-1,new Operator(W2_odd_edge),true);
  
  
  //Set edge terms
  
  //cout <<"Uodd " << Uodd << endl;
  //Uodd.exportForMatlab("Uodd.m");
  
  /*MPO dummyMPO(3);
  
  dummyMPO.setOp(0,new Operator(W1_odd),false); 
  dummyMPO.setOp(1,new Operator(W2_odd),false); 
  dummyMPO.setOp(2,new Operator(W3_odd),false); 
  
  dummyMPO.exportForMatlab("Uoddsmall.m");*/
    
  
  
  //Now fill the even MPO
  Ueven.clear();
  Ueven.initLength(this->N_total);
  
  //Edge terms at the beginning
  Ueven.setOp(0,new Operator(W1_even_edge),true);
  Ueven.setOp(1,new Operator(W2_even_edge),true);
  
  first = true;
  for(int i=3; i<this->N_total-3; i +=4)
  {
    //cout << "i=" << i << endl;
    if(first)
    {
      Ueven.setOp(i,new Operator(W1_even),true);
      Ueven.setOp(i+1,new Operator(W2_even),true);
      Ueven.setOp(i+2,new Operator(W3_even),true);
      first = false;
    }
    else
    {
      Ueven.setOp(i,&Ueven.getOp(i-4));
      Ueven.setOp(i+1,&Ueven.getOp(i-3));
      Ueven.setOp(i+2,&Ueven.getOp(i-2));
    }
  }
  //Set fermionic Identities
  first = true;
  for(int i=2; i<this->N_total; i+=4)
  {
    //cout << "i=" << i << endl;
    if(first)
      Ueven.setOp(i,new Operator(id_fermi_mpo),true);
    else
      Ueven.setOp(i,&Ueven.getOp(i-4));
  }
 
  
  //First and last fermonic identity
  //cout << "Ueven  " << Ueven << endl;
  //Ueven.setOp(0,new Operator(id_fermi_mpo),true);
  //Ueven.setOp(N_total-1,&Ueven.getOp(0));
  //cout << "Ueven after setting " << Ueven << endl;
  //cout <<"Ueven " << Ueven << endl;
  
  /*dummyMPO.setOp(0,new Operator(W1_even),false); 
  dummyMPO.setOp(1,new Operator(W2_even),false); 
  dummyMPO.setOp(2,new Operator(W3_even),false); 
  
  dummyMPO.exportForMatlab("Uevensmall.m");*/
  
  //cout << "Passt oder so " << endl;
  //Ueven.exportForMatlab("Ueven.m");
  
}//constructUnitarySubspaceMPO

void SHGIGaussLawAngular::GaussProjector(MPO &mpo, int n)
{
 
  //Clear the MPO and initialize the right length
  mpo.clear();
  mpo.initLength(this->N_total);
  
  mwArray LeftMatrix(Indices(this->d_bose,1,this->d_bose,this->d_bose));
  mwArray MiddleMatrix(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_bose));
  mwArray RightMatrix(Indices(this->d_bose,this->d_bose,this->d_bose,1));
  
  mwArray bosonic_identity = identityMatrix(this->d_bose);
  mwArray fermionic_identity = identityMatrix(this->d_fermi);
  
  bosonic_identity.reshape(Indices(this->d_bose,1,this->d_bose,1));
  fermionic_identity.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
  
  LeftMatrix.fillWithZero();
  MiddleMatrix.fillWithZero();
  RightMatrix.fillWithZero();
  
  //Fill left Matrix (basically just diagonal, on the right bond side put the value of the physical index
  for(int i=0; i<this->d_bose; i++)
    LeftMatrix.setElement(ONE_c,Indices(i,0,i,i));
  
  //Fill the Middle matrix
  for(int i=0; i<this->d_fermi; i++)
  {
    for(int left_bond_index=0; left_bond_index<this->d_bose; left_bond_index++)
    {
      for(int right_bond_index=0; right_bond_index<this->d_bose; right_bond_index++)
      {
	int factor = (i==0) ? 1 : -1;
	if(left_bond_index+(factor+pow(-1,n))/2==right_bond_index)
	  MiddleMatrix.setElement(ONE_c,Indices(i,left_bond_index,i,right_bond_index));
      }
    }
  }
  
  //Fill right Matrix (basically just diagonal, on the right bond side put the value of the physical index
  for(int i=0; i<this->d_bose; i++)
    RightMatrix.setElement(ONE_c,Indices(i,i,i,0));
  
  cout << "Left Matrix" << LeftMatrix << endl;
  cout << "Middle Matrix" << MiddleMatrix << endl;
  cout << "Right Matrix" << RightMatrix << endl;
  
  
  //Fill the MPO
  for(int i=0; i<this->N_total; i++)
  {
    if(i==((2*n-2)-1))
      mpo.setOp(i,new Operator(LeftMatrix),true);
    else if(i==(2*n-2))
      mpo.setOp(i,new Operator(MiddleMatrix),true);
    else if(i==((2*n-2)+1))
      mpo.setOp(i,new Operator(RightMatrix),true);
    else
      if((i%2)==0)
	mpo.setOp(i,new Operator(fermionic_identity),true);
	//cout << ":-)" << endl;
      else
	mpo.setOp(i,new Operator(bosonic_identity),true);
	//cout << ":-(" << endl;
  }  
  
  cout << "Warning, not yet working" << endl;
  
}

//Test of Gauss Projector, alternative which I came up with while going to Wuerzburg
void SHGIGaussLawAngular::GaussProjector2(MPO &mpo, int n)
{
  //Clear the MPO and initialize the right length
  mpo.clear();
  mpo.initLength(this->N_total);
  
  mwArray bosonic_identity = identityMatrix(this->d_bose);
  mwArray fermionic_identity = identityMatrix(this->d_fermi);
  bosonic_identity.reshape(Indices(this->d_bose,1,this->d_bose,1));
  fermionic_identity.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
  
  if(n>1 && n<(this->N_fermions))
  {
	  //Sites in between
  
	  //Prepare matrices which are used in the MPO, these are the ones where the acutal operator is located
	  mwArray LeftMatrix(Indices(this->d_bose,1,this->d_bose,this->d_bose));
	  mwArray MiddleMatrix(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_bose));
	  mwArray RightMatrix(Indices(this->d_bose,this->d_bose,this->d_bose,1));
	  
	  LeftMatrix.fillWithZero();
	  MiddleMatrix.fillWithZero();
	  RightMatrix.fillWithZero();
	  
	  //Variable for the field value
	  int field_value_left,field_value_right;
	  
	  //Fill left Matrix (basically just diagonal, on the right bond side put the value of the physical index
	  for(int i=0; i<this->d_bose; i++)
		LeftMatrix.setElement(ONE_c,Indices(i,0,i,i));
	  
	  //Fill the Middle matrix
	  for(int i=0; i<this->d_fermi; i++)
	  {
		for(int left_bond_index=0; left_bond_index<this->d_bose; left_bond_index++)
		{
		  for(int right_bond_index=0; right_bond_index<this->d_bose; right_bond_index++)
		  {
			  //Factor to determine whether one is at a spin up or spin down site (basically the value of sigma_z)
			  int factor = (i==0) ? 1 : -1;
			  field_value_left = left_bond_index - (this->d_bose-1)/2;
			  field_value_right = right_bond_index - (this->d_bose-1)/2;
			  //cout << "Field value " << field_value_left << endl;
			  if((field_value_left+(factor+pow(-1,n))/2)==field_value_right)
			  MiddleMatrix.setElement(ONE_c,Indices(i,left_bond_index,i,right_bond_index));
		  }
		}
	  }
	  
	  //Fill right Matrix (basically just diagonal, on the right bond side put the value of the physical index
	  for(int i=0; i<this->d_bose; i++)
		RightMatrix.setElement(ONE_c,Indices(i,i,i,0));
	  
	  /*cout << "Left Matrix" << LeftMatrix << endl;
	  cout << "Middle Matrix" << MiddleMatrix << endl;
	  cout << "Right Matrix" << RightMatrix << endl;*/
	  
	  
	  //Fill the MPO
	  for(int i=0; i<this->N_total; i++)
	  {
		if(i==((2*n-2)-1))
		  mpo.setOp(i,new Operator(LeftMatrix),true);
		else if(i==(2*n-2))
		  mpo.setOp(i,new Operator(MiddleMatrix),true);
		else if(i==((2*n-2)+1))
		  mpo.setOp(i,new Operator(RightMatrix),true);
		else
		  if((i%2)==0)
		    mpo.setOp(i,new Operator(fermionic_identity),true);
		    //cout << ":-)" << endl;
		  else
		    mpo.setOp(i,new Operator(bosonic_identity),true);
		    //cout << ":-(" << endl;
	  } 
	}
  else if(n==1)
  {
	  //First site
	  
	  //Prepare the matrices
	  mwArray MiddleMatrix(Indices(this->d_fermi,1,this->d_fermi,this->d_bose));
	  mwArray RightMatrix(Indices(this->d_bose,this->d_bose,this->d_bose,1));
	  
	  MiddleMatrix.fillWithZero();
	  RightMatrix.fillWithZero();
	  
	  //Variable for the field value
	  int field_value_left=0;
	  int field_value_right;
	  
	  //Fill the Middle matrix
	  for(int i=0; i<this->d_fermi; i++)
	  {
		  for(int right_bond_index=0; right_bond_index<this->d_bose; right_bond_index++)
		  {
			  //Factor to determine whether one is at a spin up or spin down site (basically the value of sigma_z)
			  int factor = (i==0) ? 1 : -1;
			  field_value_right = right_bond_index - (this->d_bose-1)/2;
			  //cout << "Field value " << field_value_left << endl;
			  if((factor+pow(-1,n))/2==field_value_right)
				MiddleMatrix.setElement(ONE_c,Indices(i,0,i,right_bond_index));
		  }
	  }
	  
	  //Fill right Matrix (basically just diagonal, on the right bond side put the value of the physical index
	  for(int i=0; i<this->d_bose; i++)
		RightMatrix.setElement(ONE_c,Indices(i,i,i,0));
	  
	  /*cout << "Left Matrix" << LeftMatrix << endl;
	  cout << "Middle Matrix" << MiddleMatrix << endl;
	  cout << "Right Matrix" << RightMatrix << endl;*/
	  
	  
	  //Fill the MPO
	  for(int i=0; i<this->N_total; i++)
	  {
		if(i==0)
		  mpo.setOp(i,new Operator(MiddleMatrix),true);
		else if(i==1)
		  mpo.setOp(i,new Operator(RightMatrix),true);
		else
		{
		  if((i%2)==0)
			mpo.setOp(i,new Operator(fermionic_identity),true);
		//cout << ":-)" << endl;
		  else
			mpo.setOp(i,new Operator(bosonic_identity),true);
		//cout << ":-(" << endl;
		}
	  }  
	  
  }
  else
  {
	  //Last site
	  mwArray LeftMatrix(Indices(this->d_bose,1,this->d_bose,this->d_bose));
	  mwArray MiddleMatrix(Indices(this->d_fermi,this->d_bose,this->d_fermi,1));
	  
	  LeftMatrix.fillWithZero();
	  MiddleMatrix.fillWithZero();
	  
	  //Fill left Matrix (basically just diagonal, on the right bond side put the value of the physical index
	  for(int i=0; i<this->d_bose; i++)
		LeftMatrix.setElement(ONE_c,Indices(i,0,i,i));
	  
	  //Variable for the field value
	  int field_value_left;
	  int field_value_right=0;
	  
	  //Fill the Middle matrix
	  for(int i=0; i<this->d_fermi; i++)
	  {
		for(int left_bond_index=0; left_bond_index<this->d_bose; left_bond_index++)
		{
			  //Factor to determine whether one is at a spin up or spin down site (basically the value of sigma_z)
			  int factor = (i==0) ? 1 : -1;
			  field_value_left = left_bond_index - (this->d_bose-1)/2;
			  if((field_value_left+(factor+pow(-1,n))/2)==field_value_right)
			  MiddleMatrix.setElement(ONE_c,Indices(i,left_bond_index,i,0));
		}
	  }
	  
	  //Fill the MPO
	  for(int i=0; i<this->N_total; i++)
	  {
		if(i==(this->N_total-2))
		  mpo.setOp(i,new Operator(LeftMatrix),true);
		else if(i==(this->N_total-1))
		  mpo.setOp(i,new Operator(MiddleMatrix),true);
		else
		  if((i%2)==0)
		    mpo.setOp(i,new Operator(fermionic_identity),true);
		    //cout << ":-)" << endl;
		  else
		    mpo.setOp(i,new Operator(bosonic_identity),true);
		    //cout << ":-(" << endl;
	  }  
  }
  
  //cout << "Warning, not yet tested" << endl;
  
}

//Get the MPO to ensure the unitary penalty thing
void SHGIGaussLawAngular::UnitaryPenaltyMPO(MPO &mpo, int n)
{
  //Clear the MPO and initialize the right length
  mpo.clear();
  mpo.initLength(this->N_total);
  
  double FIELD_TOL = 1E-12; 
  
  int l = (int) ((this->d_bose-1)/2);
  
  //Factor which I choose 1 for convenience 
  int m=1;  
  //Constant for the phase (since only a couple of phases are possible)
  double alpha = 2.0*M_PI/this->d_bose*m;
  
  mwArray bosonic_identity = identityMatrix(this->d_bose);
  mwArray fermionic_identity = identityMatrix(this->d_fermi);
  bosonic_identity.reshape(Indices(this->d_bose,1,this->d_bose,1));
  fermionic_identity.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
  
  //Variable for the field value
  complex_t field_value_left,field_value_right,field_value_middle;
  
  if(n>1 && n<this->N_fermions)
  {
    //Prepare matrices which are used in the MPO, these are the ones where the acutal operator is located
    mwArray LeftMatrix(Indices(this->d_bose,1,this->d_bose,this->d_bose));
    mwArray MiddleMatrix(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_bose));
    mwArray RightMatrix(Indices(this->d_bose,this->d_bose,this->d_bose,1));
    
    LeftMatrix.fillWithZero();
    MiddleMatrix.fillWithZero();
    RightMatrix.fillWithZero();
    
    //Possible phases (index n is from 1 to N for fermionic sites, so 1,3,5,... are spin-up sites, the rest is down sites
    complex_t phase;
    if((n%2)!=0)
      phase = exp(-I_c*alpha/2.0);
    else
      phase= exp(I_c*alpha/2.0);
    
    for(int i=0; i<this->d_bose; i++)
    {
      LeftMatrix.setElement(ONE_c,Indices(i,0,i,i));
      RightMatrix.setElement(ONE_c,Indices(i,i,i,0));
    }
    
    for(int i=0; i<this->d_fermi; i++)
    {
      for(int left_bond_index=0; left_bond_index<this->d_bose; left_bond_index++)
      {
	for(int right_bond_index=0; right_bond_index<this->d_bose; right_bond_index++)
	{
	  //Factor to determine whether one is at a spin up or spin down site (basically the value of sigma_z)
	  field_value_left = exp(-I_c*(-l+left_bond_index)*alpha);
	  field_value_right = exp(I_c*(-l+right_bond_index)*alpha);
	  field_value_middle = (i==0) ? exp(-I_c*alpha/2.0) : exp(I_c*alpha/2.0);
	  //cout << "Field value " << field_value_left << endl;
	  if(abs(field_value_left*field_value_middle*field_value_right-phase)< FIELD_TOL)
	  {
	    //cout << "i_l " << left_bond_index << ", i_r " << right_bond_index << ", i_m " << i << ", Feld_l " << field_value_left <<", Feld_r " << field_value_right << ", Feld_m " << field_value_middle << ", product " << field_value_left*field_value_middle*field_value_right << ", Phase " << phase << ", Differenz " << abs(field_value_left*field_value_middle*field_value_right-phase) << endl;	 
	    MiddleMatrix.setElement(ONE_c,Indices(i,left_bond_index,i,right_bond_index));
	  }
	}
      }
    }
    //cout << "Middle Matrix " << MiddleMatrix << endl;
    
    //Fill the MPO
    for(int i=0; i<this->N_total; i++)
    {
      if(i==((2*n-2)-1))
	mpo.setOp(i,new Operator(LeftMatrix),true);
      else if(i==(2*n-2))
	mpo.setOp(i,new Operator(MiddleMatrix),true);
      else if(i==((2*n-2)+1))
	mpo.setOp(i,new Operator(RightMatrix),true);
      else
	if((i%2)==0)
	  mpo.setOp(i,new Operator(fermionic_identity),true);
	  //cout << ":-)" << endl;
	else
	  mpo.setOp(i,new Operator(bosonic_identity),true);
	    //cout << ":-(" << endl;
    }
  }
  else if(n==1)
  {
    //Prepare matrices which are used in the MPO, these are the ones where the acutal operator is located
    mwArray MiddleMatrix(Indices(this->d_fermi,1,this->d_fermi,this->d_bose));
    mwArray RightMatrix(Indices(this->d_bose,this->d_bose,this->d_bose,1));
    
    MiddleMatrix.fillWithZero();
    RightMatrix.fillWithZero();
    
    //Possible phases (index n is from 1 to N for fermionic sites, so 1,3,5,... are spin-up sites, the rest is down sites
    complex_t phase;
    if((n%2)!=0)
      phase = exp(-I_c*alpha/2.0);
    else
      phase= exp(I_c*alpha/2.0);
    
    for(int i=0; i<this->d_bose; i++)
    {
      RightMatrix.setElement(ONE_c,Indices(i,i,i,0));
    }
    
    for(int i=0; i<this->d_fermi; i++)
    {
      for(int right_bond_index=0; right_bond_index<this->d_bose; right_bond_index++)
      {
	//Factor to determine whether one is at a spin up or spin down site (basically the value of sigma_z)
	field_value_right = exp(I_c*(-l+right_bond_index)*alpha);
	field_value_middle = (i==0) ? exp(-I_c*alpha/2.0) : exp(I_c*alpha/2.0);
	//cout << "Field value " << field_value_left << endl;
	if(abs(field_value_middle*field_value_right-phase)< FIELD_TOL)
	{
	  //cout << "i_r " << right_bond_index << ", i_m " << i << ", Feld_r " << field_value_right << ", Feld_m " << field_value_middle << ", product " << field_value_middle*field_value_right << ", Phase " << phase << ", Differenz " << abs(field_value_middle*field_value_right-phase) << endl;	  
	  MiddleMatrix.setElement(ONE_c,Indices(i,0,i,right_bond_index));
	}
      }
    }
    
    //Fill the MPO
    for(int i=0; i<this->N_total; i++)
    {
      if(i==0)
	mpo.setOp(i,new Operator(MiddleMatrix),true);
      else if(i==1)
	mpo.setOp(i,new Operator(RightMatrix),true);
      else
	if((i%2)==0)
	  mpo.setOp(i,new Operator(fermionic_identity),true);
	  //cout << ":-)" << endl;
	else
	  mpo.setOp(i,new Operator(bosonic_identity),true);
	    //cout << ":-(" << endl;
    }
  }
  else
  {
    //Prepare matrices which are used in the MPO, these are the ones where the acutal operator is located
    mwArray LeftMatrix(Indices(this->d_bose,1,this->d_bose,this->d_bose));
    mwArray MiddleMatrix(Indices(this->d_fermi,this->d_bose,this->d_fermi,1));
    
    LeftMatrix.fillWithZero();
    MiddleMatrix.fillWithZero();
    
    //Possible phases (index n is from 1 to N for fermionic sites, so 1,3,5,... are spin-up sites, the rest is down sites
    complex_t phase;
    if((n%2)!=0)
      phase = exp(-I_c*alpha/2.0);
    else
      phase= exp(I_c*alpha/2.0);
    
    for(int i=0; i<this->d_bose; i++)
    {
      LeftMatrix.setElement(ONE_c,Indices(i,0,i,i));
    }
    
    for(int i=0; i<this->d_fermi; i++)
    {
      for(int left_bond_index=0; left_bond_index<this->d_bose; left_bond_index++)
      {
	//Factor to determine whether one is at a spin up or spin down site (basically the value of sigma_z)
	field_value_left = exp(-I_c*(-l+left_bond_index)*alpha);
	field_value_middle = (i==0) ? exp(-I_c*alpha/2.0) : exp(I_c*alpha/2.0);
	//cout << "Field value " << field_value_left << endl;
	if(abs(field_value_left*field_value_middle-phase)< FIELD_TOL)
	{
	  //cout << "i_l " << left_bond_index << ", i_m " << i << ", Feld_l " << field_value_left << ", Feld_m " << field_value_middle << ", product " << field_value_left*field_value_middle << ", Phase " << phase << ", Differenz " << abs(field_value_left*field_value_middle-phase) << endl;
	  MiddleMatrix.setElement(ONE_c,Indices(i,left_bond_index,i,0));
	}
      }
    }
    //cout << "Middle Matrix " << MiddleMatrix << endl;
    
    //Fill the MPO
    for(int i=0; i<this->N_total; i++)
    {
      if(i==(this->N_total-2))
	mpo.setOp(i,new Operator(LeftMatrix),true);
      else if(i==(this->N_total-1))
	mpo.setOp(i,new Operator(MiddleMatrix),true);
      else
	if((i%2)==0)
	  mpo.setOp(i,new Operator(fermionic_identity),true);
	  //cout << ":-)" << endl;
	else
	  mpo.setOp(i,new Operator(bosonic_identity),true);
	    //cout << ":-(" << endl;
    }
  } 
}

//Build a projector for the charge zero subspace
//Idea: I start counting from the left and the right and encode the current charge in the bond dimension. Starting from left and right limits the bond dimension to N_fermion/2. In the middle I have a site where the left and right part must sum up to a zero.
void SHGIGaussLawAngular::ChargeProjectorMPO(MPO &mpo)
{
  //A vector to store all the different matrices
  vector<mwArray> Matrices;
  Matrices.resize(this->N_total);
  
  //Determine the splitting spin site (total number of sites should be odd)
  if((this->N_total%2)==0)
  {
    cout << "Error in SHGIGaussLawAngular::ChargeProjectorMPO(MPO &mpo), total number of sites must be odd, instead it is " << this->N_total << endl;
    exit(666);
  }
  int split_site = (this->N_total-1)/2;
  
  //Variable to store the current bond dimension
  int bond_count;
  
  //Spin up is encodes as +1, spin down as -1
  
  //Left half of the chain
  bond_count=1;
  for(int i=0; i<split_site; i++)
  {
    if((i%2)==0)
    {
      //Fermionic site
      mwArray Matrix(Indices(this->d_fermi,bond_count,this->d_fermi,bond_count+1));
      Matrix.fillWithZero();
      
      Matrices[i]=Matrix;
      bond_count ++;
    }
    else
    {
      //Bosonic site
      mwArray Matrix(Indices(this->d_bose,bond_count,this->d_bose,bond_count));
      Matrix.fillWithZero();
      for(int k=0; k<this->d_bose; k++)
	for(int l=0; l<bond_count; l++)
	  Matrix.setElement(ONE_c,Indices(k,l,k,l));
      
      Matrices[i]=Matrix;
    }    
  }
  
  //Right half of the chain
  bond_count=1;
  for(int i=(this->N_total-1); i>split_site; i--)
  {
    if((i%2)==0)
    {
      //Fermionic site
      mwArray Matrix(Indices(this->d_fermi,bond_count+1,this->d_fermi,bond_count));
      Matrix.fillWithZero();
      
      Matrices[i]=Matrix;
      bond_count ++;
    }
    else
    {
      //Bosonic site
      mwArray Matrix(Indices(this->d_bose,bond_count,this->d_bose,bond_count));
      Matrix.fillWithZero();
      for(int k=0; k<this->d_bose; k++)
	for(int l=0; l<bond_count; l++)
	  Matrix.setElement(ONE_c,Indices(k,l,k,l));
      
      Matrices[i]=Matrix;
    }
  }
  
  //Connecting site
  mwArray Matrix(Indices(this->d_bose,bond_count,this->d_bose,bond_count));
  Matrix.fillWithZero();
  
  Matrices[split_site]=Matrix;
  
  //Fill the MPO
  mpo.clear();
  mpo.initLength(this->N_total);
  for(int i=0; i<this->N_total;i++)
    mpo.setOp(i,new Operator(Matrices[i]),true); 
}

void SHGIGaussLawAngular::EvolutionMPO(MPO &mpo,double dt)
{
  /************************************************************************************************
   * Provide the operators
  ************************************************************************************************/
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  mwArray id_fermionic_modified = 1.0/((double)this->N_total)*id_fermionic;
  
  complex_t data_x[] = {ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigma_x(Indices(2,2),data_x);
  
  complex_t data_y[] = {ZERO_c,I_c,-1.0*I_c,ZERO_c};
  mwArray sigma_y(Indices(2,2),data_y);
  
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
  mwArray sigma_minus(Indices(2,2),data_sm);
  
  complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
  mwArray sigma_plus(Indices(2,2),data_sp);
  
  complex_t data_idpsz[] = {2.0*ONE_c,ZERO_c,ZERO_c,ZERO_c};  
  mwArray idpsz(Indices(2,2),data_idpsz); 
  
  //Provide Operators for gauge field
  mwArray Op1(Indices(this->d_bose,this->d_bose));
  mwArray Op2(Indices(this->d_bose,this->d_bose));
  mwArray Theta(Indices(this->d_bose,this->d_bose));
  mwArray L_operator(Indices(this->d_bose,this->d_bose));
  mwArray Lsq_operator(Indices(this->d_bose,this->d_bose));
  mwArray id_bosonic = identityMatrix(this->d_bose);
  mwArray  Lsqpmodid;
  
  L_operator.fillWithZero();
  Lsq_operator.fillWithZero();
  Op1.fillWithZero();
  Op2.fillWithZero();
  Theta.fillWithZero();
  
  int offset = (int) floor(this->d_bose/2.0);   
  
  for(int i=0; i<this->d_bose; ++i){
    for(int j=0; j<this->d_bose; ++j){
      if(i==j)
	L_operator.setElement((double) (i-offset), 0.0,Indices(i,j));
    }
  }
  
  //Construct the theta Operator
  double xval,yval,k_val;
  complex_t tmp;
  for(int i=0; i<this->d_bose; ++i)
  {
    xval = -1.0*(this->d_bose-1)/2.0 + (double) i;
    for(int j=0; j<this->d_bose; ++j)
    {
      yval = -1.0*(this->d_bose-1)/2.0 + (double) j;
      tmp=ZERO_c;
      for(int k=0; k<this->d_bose; ++k)
      {
	k_val= -1.0*(this->d_bose-1)/2.0 + (double) k;
	//cout << "x=" << xval << " y=" << yval << " k=" << k_val << endl;
	tmp = tmp + 2.0*M_PI/((double) this->d_bose)*k_val*exp(-2.0*I_c*M_PI/((double) this->d_bose)*k_val*(xval-yval));
      }
      Theta.setElement(tmp/((double) this->d_bose),Indices(i,j));
    }
  }
  //Get the gauge invariant version of the operators  
  wrapper::expm(Theta,Op1,I_c);
  wrapper::expm(Theta,Op2,-1.0*I_c);
  
  //Now get get the square of the L-operator
  Lsq_operator = L_operator*L_operator;
 
  
  
  //I need -i*dt all the time, that's why I save it WITH MINUS SIGN
  complex_t taylor_coefficient;
  taylor_coefficient = -1.0*dt*I_c;
  
  Lsqpmodid = 1.0/((double)this->N_total)*id_bosonic + taylor_coefficient*Lsq_operator;
  
  /************************************************************************************************
   * Prepare the operator arrays
  ************************************************************************************************/
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(4,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<this->d_bose; ++i)
  {
    for(int j=0; j<this->d_bose; ++j)
    {
      //Identity
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));
      //01
      Bosonic_operators.setElement(Op1.getElement(Indices(i,j)),Indices(1,i,j));
      //02
      Bosonic_operators.setElement(Op2.getElement(Indices(i,j)),Indices(2,i,j));
      //Lsquared+modified identity
      Bosonic_operators.setElement(Lsqpmodid.getElement(Indices(i,j)),Indices(3,i,j));
    }
  }
  
  Bosonic_operators.reshape(Indices(4,this->d_bose*this->d_bose));

  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(5,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();


  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
      //Modified Identity
      Fermionic_operators.setElement(id_fermionic_modified.getElement(Indices(i,j)),Indices(1,i,j));
      //id+sigma_z
      Fermionic_operators.setElement(idpsz.getElement(Indices(i,j)),Indices(2,i,j));
      //sigma_plus
      Fermionic_operators.setElement(sigma_plus.getElement(Indices(i,j)),Indices(3,i,j));
      //sigma_minus
      Fermionic_operators.setElement(sigma_minus.getElement(Indices(i,j)),Indices(4,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(5,this->d_fermi*this->d_fermi));
  
  /************************************************************************************************
   * Now prepare the matrices
  ************************************************************************************************/
  //Bond dimension of my MPO
  int D_bond = 4;
  double sign;
  
  //cout << "Filling matrices" << endl;
  //First fermionic site
  mwArray fermi_first(Indices(1,D_bond,5)); fermi_first.fillWithZero();
  sign=-1.0;
  //Identity
  fermi_first.setElement(ONE_c,Indices(0,0,0));  
  //Modified identity
  fermi_first.setElement(ONE_c,Indices(0,3,1));  
  //I+sigma_z(n)
  fermi_first.setElement(this->mu/2.0*sign*ONE_c*taylor_coefficient,Indices(0,3,2));  
  //sigma_plus
  fermi_first.setElement(ONE_c,Indices(0,1,3));  
  //sigma_minus
  fermi_first.setElement(ONE_c,Indices(0,2,4));   
  
  //Even fermionic sites in between
  mwArray fermi_even(Indices(D_bond,D_bond,5)); fermi_even.fillWithZero();
  sign=1.0;
  //Identity
  fermi_even.setElement(ONE_c,Indices(0,0,0));
  fermi_even.setElement(ONE_c,Indices(3,3,0));  
  //Modified identity
  fermi_even.setElement(ONE_c,Indices(0,3,1));  
  //I+sigma_z(n)
  fermi_even.setElement(this->mu/2.0*sign*ONE_c*taylor_coefficient,Indices(0,3,2));  
  //sigma_plus
  fermi_even.setElement(ONE_c,Indices(0,1,3));
  fermi_even.setElement(ONE_c,Indices(2,3,3));  
  //sigma_minus
  fermi_even.setElement(ONE_c,Indices(0,2,4)); 
  fermi_even.setElement(ONE_c,Indices(1,3,4)); 
  
  //Odd fermionic sites in between
  mwArray fermi_odd(Indices(D_bond,D_bond,5)); fermi_odd.fillWithZero();
  sign=-1.0;
  //Identity
  fermi_odd.setElement(ONE_c,Indices(0,0,0));
  fermi_odd.setElement(ONE_c,Indices(3,3,0));  
  //Modified identity
  fermi_even.setElement(ONE_c,Indices(0,3,1));  
  //I+sigma_z(n)
  fermi_odd.setElement(this->mu/2.0*sign*ONE_c*taylor_coefficient,Indices(0,3,2));  
  //sigma_plus
  fermi_odd.setElement(ONE_c,Indices(0,1,3));
  fermi_odd.setElement(ONE_c,Indices(2,3,3));  
  //sigma_minus
  fermi_odd.setElement(ONE_c,Indices(0,2,4)); 
  fermi_odd.setElement(ONE_c,Indices(1,3,4));
  
  //Bosonic sites in between
  mwArray bose(Indices(D_bond,D_bond,4)); bose.fillWithZero();
  //Identity
  bose.setElement(ONE_c,Indices(0,0,0));
  bose.setElement(ONE_c,Indices(3,3,0));  
  //Op1
  bose.setElement(this->x*taylor_coefficient*ONE_c,Indices(1,1,1));  
  //Op2
  bose.setElement(this->x*taylor_coefficient*ONE_c,Indices(2,2,2));  
  //L^2+ modified identity
  bose.setElement(ONE_c,Indices(0,3,3));  
  
  //Last fermionic site
  mwArray fermi_last(Indices(D_bond,1,5)); fermi_last.fillWithZero();
  sign = pow(-1.0,((double)this->N_fermions));
  //Identity
  fermi_last.setElement(ONE_c,Indices(3,0,0));
  //Modified identity
  fermi_last.setElement(ONE_c,Indices(0,0,1));
  //I+sigma_z(n)
  fermi_last.setElement(this->mu/2.0*sign*ONE_c*taylor_coefficient,Indices(0,0,2));
  //sigma_plus
  fermi_last.setElement(ONE_c,Indices(2,0,3));
  //sigma_minus
  fermi_last.setElement(ONE_c,Indices(1,0,4)); 
 
  
  /*cout << "fermi_first " << fermi_first << endl;
  cout << "fermi_even " << fermi_even << endl;
  cout << "fermi_odd " << fermi_odd << endl;
  cout << "fermi_last " << fermi_last << endl;
  cout << "bose " << bose << endl;*/
  
  
  //Reshape operators
  fermi_first.reshape(Indices(1*D_bond,5));
  fermi_even.reshape(Indices(D_bond*D_bond,5));
  fermi_odd.reshape(Indices(D_bond*D_bond,5));
  fermi_last.reshape(Indices(1*D_bond,5));
  bose.reshape(Indices(D_bond*D_bond,4));
  
  

  
  /************************************************************************************************
  * Contract and fill the MPO
  ************************************************************************************************/
  //cout << "Filling MPO" << endl;
  mpo.clear();
  mpo.initLength(this->N_total);  
  int fermi_count=1;
  bool fermi_odd_set=false;
  bool fermi_even_set=false;
  bool bose_set=false;
  mwArray res;
  for(int i=0; i<this->N_total; i++)
  {
    if(i==0)
    {
      //First fermionic site
      res=reshape(fermi_first*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
      mpo.setOp(i,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
    else if(i==(this->N_total-1))
    {
      //Last fermionic site
      res=reshape(fermi_last*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
      mpo.setOp(i,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
    else if((i%2)==0)
    {
      if((fermi_count%2)==0)
      {
	//Even Fermionic sites in between
	if(fermi_even_set)
	{
	  //Operator already set, just get pointer (keep in mind to go to the next even fermionic site)
	  mpo.setOp(i,&mpo.getOp(i-4),false);
	}
	else
	{
	  //Set the operator and the flag
	  res=reshape(fermi_even*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	  mpo.setOp(i,new Operator(permute(res,Indices(3,1,4,2))),true);
	  fermi_even_set = true;
	}
      }
      else
      {
	//Odd Fermionic sites in between
	if(fermi_odd_set)
	{
	  //Operator already set, just get pointer (keep in mind to go to the next odd fermionic site)
	  mpo.setOp(i,&mpo.getOp(i-4),false);
	}
	else
	{
	  //Set the operator and the flag
	  res=reshape(fermi_odd*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	  mpo.setOp(i,new Operator(permute(res,Indices(3,1,4,2))),true);
	  fermi_odd_set = true;
	}
      }
      fermi_count++;
    }
    else
    {
      //Bosonic sites in between
      if(bose_set)
      {
	//Operator already set, just get pointer
	mpo.setOp(i,&mpo.getOp(i-2),false);
      }
      else
      {
	//Set the operator and the flag
	res=reshape(bose*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose,this->d_bose));
	mpo.setOp(i,new Operator(permute(res,Indices(3,1,4,2))),true);
	bose_set = true;
      }
    }
  }
}