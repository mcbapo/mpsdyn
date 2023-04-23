#include "SchwingerHamiltonianAGIGaussLaw.h"
#include "Indices.h"
#include "Contractor.h"

using namespace std;
using namespace shrt;


SchwingerHamiltonianAGIGaussLaw::SchwingerHamiltonianAGIGaussLaw(int N_, int d_,int D_,double mu_,double x_,double e_,double lambda_,unsigned int exp_order_):hamil(2*N_-1),charge(2*N_-1)
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
  if(exp_order_<170)
    exp_order = exp_order_;
  else
  {
    cout << "!!! Warning in SchwingerHamiltonianAGIGaussLaw::SchwingerHamiltonianAGIGaussLaw: value for the order of the exponential approximation is too large, setting to maximum possible value !!! " << endl;
    exp_order = 169;
  }
    
  
  
  cout << "Constructing Hamiltonian with " << endl
       << "N_fermions = " << N_fermions << endl
       << "N_total = " << N_total << endl
       << "D = " << D << endl
       << "mu = " << mu << endl
       << "x = " << x << endl
       << "d = " << d_bose << endl
       << "e = " << e << endl
       << "la = " << lambda << endl 
       << "exp-order = " << exp_order << endl;
  
  //Construct an MPO for the hamiltonian
  initOperators();
  //initMPO();
  initChargeoperator();
}

//Destructor
SchwingerHamiltonianAGIGaussLaw::~SchwingerHamiltonianAGIGaussLaw(){
    hamil.clear();
    charge.clear();    
}

////Construct an explicit MPO for the Hamiltonian
void SchwingerHamiltonianAGIGaussLaw::initOperators()
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
  mwArray Op1(Indices(this->d_bose,this->d_bose));
  mwArray Op2(Indices(this->d_bose,this->d_bose));
  mwArray argOp1(Indices(this->d_bose,this->d_bose));
  mwArray argOp2(Indices(this->d_bose,this->d_bose));
  mwArray L_operator(Indices(this->d_bose,this->d_bose));
  mwArray Lsq_operator(Indices(this->d_bose,this->d_bose));
  mwArray id_bosonic = identityMatrix(this->d_bose);
  
  Op1.fillWithZero();
  Op2.fillWithZero();
  L_operator.fillWithZero();
  Lsq_operator.fillWithZero();
  
  for(int i=0; i<this->d_bose; ++i){
    for(int j=0; j<this->d_bose; ++j){
      if((i+1)==j){
	argOp1.setElement(0.0,e*sqrt(j),Indices(i,j));
	argOp2.setElement(0.0,-e*sqrt(j),Indices(i,j));
	L_operator.setElement(0.0,-sqrt(j)/2.0/e,Indices(i,j));
      }
      else if(i==(j+1)){
	argOp1.setElement(0.0,e*sqrt(i),Indices(i,j));
	argOp2.setElement(0.0,-e*sqrt(i),Indices(i,j));
	L_operator.setElement(0.0,sqrt(i)/2.0/e,Indices(i,j));
      }
    }
  }
  
  //Approximate epx(argOp1) and exp(argOp2) up to the desired order
  mwArray TempOp1 = identityMatrix(this->d_bose);
  mwArray TempOp2 = identityMatrix(this->d_bose);
  double faculty = 1;
  for(int i=0; i<=this->exp_order; i++)
  {
    Op1 = Op1 + 1/faculty*TempOp1;
    Op2 = Op2 + 1/faculty*TempOp2;
    
    TempOp1 = TempOp1*argOp1;
    TempOp2 = TempOp2*argOp2;
    
    faculty = faculty*(i+1.0);
  }
  
  
  //Now get get the square of the L-operator
  Lsq_operator = L_operator*L_operator;
  
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
    cout << "Constructing Hamiltonian MPO for site k = " << k << ":" << endl;       
    Bmatrices.clear();
    res.clear();
    
    if((k%2)==1)
    {
      //Bosonic site(odd sites)
      Bmatrices.resize(Indices(D_bond,D_bond,8));
      Bmatrices.fillWithZero();
      
      sign = pow(-1.0,fermion_counter);
      cout << "Sign bosonic site: " << sign << endl;
      
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
      
      cout << "Sign fermionic site " << sign << endl;
      
      if(k==0)
      {
	//First fermionic site
	cout<< "First fermionic site " << endl;
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
	cout << "Last fermionic site " << endl;
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
  cout << hamil << endl;
  hamil.exportForMatlab("HPenaltyAGI.m");
}// initOperators

////Construct an explicit MPO for the total charge
void SchwingerHamiltonianAGIGaussLaw::initChargeoperator()
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

//Construct an explicit MPO to be able to get the electric field and check the Gauss Law
void SchwingerHamiltonianAGIGaussLaw::constructLopMPO(MPO &mpo, int n)
{
  cout << "Constructing L("<<n<<")" << endl;
  
  //As n is a DMRG site and L(n) acts on bosonic sites only, check if n is a bosonic index
  if((n%2)!=1)
  {
    cerr << "Error in SchwingerHamiltonianAGIGaussLaw::constructLopMPO(MPO &mpo, int n): Index " << n << " is not a valid bosonic site" << endl;
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
  
  for(int i=0; i<this->d_bose; ++i){
    for(int j=0; j<this->d_bose; ++j){
      if((i+1)==j){
	L_operator.setElement(0.0,-sqrt(j)/2.0/e,Indices(i,j));
      }
      else if(i==(j+1)){
	L_operator.setElement(0.0,sqrt(i)/2.0/e,Indices(i,j));
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
}

//Construct an explicit MPO to be able to check the Gauss Law
void SchwingerHamiltonianAGIGaussLaw::constructGaussMPO(MPO &mpo, int n)
{
  cout << "Constructing GaussMPO for site " << n << endl;
  
  //As n is a DMRG site and sigma_z acts on fermionic sites only, check if n is a bosonic index
  if((n%2)==1)
  {
    cerr << "Error in SSchwingerHamiltonianAGIGaussLaw::constructGaussMPO(MPO &mpo, int n): Index " << n << " is not a valid fermionic site" << endl;
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

//Construct an explicit MPO to check the local spin
void SchwingerHamiltonianAGIGaussLaw::constructsigmazMPO(MPO &mpo, int n)
{
  cout << "Constructing sigma_z("<< n << ")" << endl;
  
  //As n is a DMRG site and sigma_z acts on fermionic sites only, check if n is a bosonic index
  if((n%2)==1)
  {
    cerr << "Error in SSchwingerHamiltonianAGIGaussLaw::constructsigmazMPO(MPO &mpo, int n): Index " << n << " is not a valid fermionic site" << endl;
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

//Construct an explicit MPO for the Condensate
void SchwingerHamiltonianAGIGaussLaw::constructCondensateMPO(MPO &mpo)
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
	Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));	res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
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
void SchwingerHamiltonianAGIGaussLaw::constructPenaltyMPO(MPO &mpo)
{
  cout << "Constructing Penalty MPO" << endl;
  
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
  
  for(int i=0; i<this->d_bose; ++i){
    for(int j=0; j<this->d_bose; ++j){
      if((i+1)==j){
	L_operator.setElement(0.0,-sqrt(j)/2.0/e,Indices(i,j));
      }
      else if(i==(j+1)){
	L_operator.setElement(0.0,sqrt(i)/2.0/e,Indices(i,j));
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
      Bmatrices.setElement(-2.0*sign*this->lambda*ONE_c,Indices(0,4,2));
      Bmatrices.setElement(-1.0*ONE_c,Indices(1,4,2));
      Bmatrices.setElement(-2.0*ONE_c,Indices(3,4,2));
      //Lsquared
      Bmatrices.setElement(2.0*this->lambda*ONE_c,Indices(0,4,3));  
      
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
	Bmatrices.setElement(0.5*this->lambda*ONE_c,Indices(0,4,0));
	
	//sigma_z
	Bmatrices.setElement(this->lambda*ONE_c,Indices(0,1,1));
	Bmatrices.setElement(0.5*this->lambda*sign*ONE_c,Indices(0,4,1));
	
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
	Bmatrices.setElement(0.5*this->lambda*ONE_c,Indices(0,0,0));
	Bmatrices.setElement(ONE_c,Indices(4,0,0));
	//sigma_z
	Bmatrices.setElement(0.5*sign*this->lambda*ONE_c,Indices(0,0,1));
	Bmatrices.setElement(this->lambda*ONE_c,Indices(2,0,1));
	
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
	Bmatrices.setElement(0.5*this->lambda*ONE_c,Indices(0,4,0));
	Bmatrices.setElement(this->lambda*ONE_c,Indices(2,3,0));
	Bmatrices.setElement(ONE_c,Indices(4,4,0));
	//sigma_z
	Bmatrices.setElement(this->lambda*ONE_c,Indices(0,1,1));
	Bmatrices.setElement(0.5*this->lambda*sign*ONE_c,Indices(0,4,1));
	Bmatrices.setElement(this->lambda*ONE_c,Indices(2,4,1));
	
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
