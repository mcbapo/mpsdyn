#include "SchwingerHamiltonianNonGaugeinvariant.h"
#include "Indices.h"
#include "Contractor.h"


using namespace std;
using namespace shrt;


SchwingerHamiltonianNonGaugeinvariant::SchwingerHamiltonianNonGaugeinvariant(int N_, int d_,int D_,double mu_,double x_,double e_):hamil(2*N_-1),charge(2*N_-1)
{
  N_fermions = N_;
  //As a bosonic link is between two fermions, the total number of sites in the chain is given by
  N_total = 2*N_-1;
  D = D_;
  mu = mu_;
  x = x_;
  d = d_;
  e = e_;
  
  cout << "Constructing Hamiltonian with " << endl
       << "N_fermions = " << N_fermions << endl
       << "N_total = " << N_total << endl
       << "D = " << D << endl
       << "mu = " << mu << endl
       << "x = " << x << endl
       << "d = " << d << endl
       << "e = " << e << endl;
  
  //Construct an MPO for the hamiltonian
  initMPO();
  initChargeoperator();
}

//Destructor
SchwingerHamiltonianNonGaugeinvariant::~SchwingerHamiltonianNonGaugeinvariant(){
    hamil.clear();
    charge.clear();
}

//Construct the MPO representation of the Schwinger Hamiltonian
void SchwingerHamiltonianNonGaugeinvariant::initMPO(void)
{
  cout << "Constructing the MPO representation of the Hamiltonian" << endl;
  
  /*************************************************
   * 		Preliminary setup
   * **********************************************/  
  
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
  
//   cout << " sigma_x=" << sigma_x << endl;
//   cout << " sigma_y=" << sigma_y << endl;
//   cout << " sigma_z=" << sigma_z << endl;
//   
//   cout << " Identitaet= " << id_fermionic << endl;
//   cout << " sigma_+=" << sigma_plus << endl;
//   cout << " sigma_-=" << sigma_minus << endl;
//   
//   cout << " (I+sigma_z)=" << idpsz << endl;
  
  //Provide Operators for gauge field
  mwArray Op1(Indices(this->d,this->d));
  mwArray Op2(Indices(this->d,this->d));
  mwArray Lsq_operator(Indices(this->d,this->d));
  mwArray id_bosonic = identityMatrix(this->d);
  
  Op1.fillWithZero();
  Op2.fillWithZero();
  Lsq_operator.fillWithZero();
  
  for(int i=0; i<this->d; ++i){
    for(int j=0; j<this->d; ++j){
      if(i==j){
	Op1.setElement(1.0,0.0,Indices(i,j));
	Op2.setElement(1.0,0.0,Indices(i,j));
      }
      else if((i+1)==j){
	Op1.setElement(0.0,e*sqrt(j),Indices(i,j));
	Op2.setElement(0.0,-e*sqrt(j),Indices(i,j));
	Lsq_operator.setElement(0.0,-sqrt(j)/2.0/e,Indices(i,j));
      }
      else if(i==(j+1)){
	Op1.setElement(0.0,e*sqrt(i),Indices(i,j));
	Op2.setElement(0.0,-e*sqrt(i),Indices(i,j));
	Lsq_operator.setElement(0.0,sqrt(i)/2.0/e,Indices(i,j));
      }
    }
  }
  
  //Now get get the square of the L-operator
  Lsq_operator = Lsq_operator*Lsq_operator;
  
//  cout << endl << endl << "!!! ACHTUNG, L^2, O_1 und O_2 auf 0 gesetzt !!!" << endl << endl;
//  Op1.fillWithZero();
//   Op2.fillWithZero();
//   Lsq_operator.fillWithZero();
//   id_bosonic.fillWithZero();
  
//   cout << "Op1=" << Op1 << endl;
//   cout << "Op2=" << Op2 << endl;
//   cout << "Lsq_operator=" << Lsq_operator << endl;
//   cout << "Identitaet=" << id_bosonic << endl;
  
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(7,d,d));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<d; ++i)
  {
    for(int j=0; j<d; ++j)
    { 
      //Identiy
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));
      //Operator O_1
      Bosonic_operators.setElement(Op1.getElement(Indices(i,j)),Indices(3,i,j));
      //Operator O_2
      Bosonic_operators.setElement(Op2.getElement(Indices(i,j)),Indices(4,i,j));
      //Operator L^2
      Bosonic_operators.setElement(Lsq_operator.getElement(Indices(i,j)),Indices(5,i,j));
    }
  }
  
  Bosonic_operators.reshape(Indices(7,d*d));
  
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(7,2,2));
  Fermionic_operators.fillWithZero();
  
//  idpsz.fillWithZero();
//  sigma_plus.fillWithZero();
//  sigma_minus.fillWithZero();
//  id_fermionic.fillWithZero();
  
//  cout << "Fermionische Identitaet: " << id_fermionic << endl;
  
  for(int i=0; i<2; ++i)
  {
    for(int j=0; j<2; j++)
    { 
      //Identiy
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
      //Operator sigma_+
      Fermionic_operators.setElement(sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
      //Operator sigma_-
      Fermionic_operators.setElement(sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
      //Operator id+sigma_z
      Fermionic_operators.setElement(idpsz.getElement(Indices(i,j)),Indices(6,i,j));
    }
  } 
  
  Fermionic_operators.reshape(Indices(7,2*2));
  
  /*************************************************
  * 		Actual MPO
  * **********************************************/
  
  //Bond dimension for the MPO Matrices
  int D_bond = 6;
  //Counter to keep track, whether one is at an even fermionic site or an odd fermionic site
  int fermion_counter=0;
  double sign;
  
  
  for(int k=0; k<N_total;++k)
  {
    cout << "k = " << k << ":" << endl;
       
    //Bosonic site(odd sites)
    if((k%2)==1)
    {
      cout << "Bosonic site" << endl << endl;
      //Array to store the matrices which contain the information about the operators
      mwArray Bmatrices(Indices(D_bond,D_bond,7));
      Bmatrices.fillWithZero();
      
      //Identity
      Bmatrices.setElement(ONE_c,Indices(0,0,0));
      Bmatrices.setElement(ONE_c,Indices(5,5,0));
      
      //Operator O_1
      //Bmatrices.setElement(ONE_c,Indices(1,3,3));
      Bmatrices.setElement(x,0.0,Indices(1,3,3));
      
      //Operator O_2
      //Bmatrices.setElement(ONE_c,Indices(2,4,4));
      Bmatrices.setElement(x,0.0,Indices(2,4,4));
      
      //Operator L^2
      Bmatrices.setElement(ONE_c,Indices(0,5,5));  
      
      //Contract and set the MPO
      Bmatrices.reshape(Indices(D_bond*D_bond,7));
      mwArray res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond,D_bond,d,d));
      hamil.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      
    }
    //Fermonic site (even sites)
    if((k%2)==0)
    {
      //Start with 1
      fermion_counter ++;
      //Set sign
      sign = ((fermion_counter%2)==0) ? 1.0 : -1.0;
      
      //First and last fermionic site have to be treated differently
      if(k==0)
      {
	cout << "First Fermionic site, sign = " << sign << endl << endl;
	mwArray Bmatrices(Indices(1,D_bond,7));
	Bmatrices.fillWithZero();
	
	//Identity
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	
	//sigma_+
	Bmatrices.setElement(ONE_c,Indices(0,1,1));
	
	//sigma_-
	Bmatrices.setElement(ONE_c,Indices(0,2,2));	
	
	//0.5*(I+sigma_z)
	//Bmatrices.setElement(sign*ONE_c,Indices(0,5,6));
	Bmatrices.setElement(sign*mu/2.0,0.0,Indices(0,5,6));
	
	//Contract and set the MPO
	Bmatrices.reshape(Indices(D_bond,7));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(1,D_bond,2,2));
	hamil.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	
      }
      else if(k==(N_total-1))
      {
	cout << "Last Fermionic site, sign = " << sign << endl << endl;
	mwArray Bmatrices(Indices(D_bond,1,7));
	Bmatrices.fillWithZero();
	
	//Identity
	Bmatrices.setElement(ONE_c,Indices(5,0,0));
	
	//sigma_+
	Bmatrices.setElement(ONE_c,Indices(4,0,1));
	
	//sigma_-
	Bmatrices.setElement(ONE_c,Indices(3,0,2));	
	
	//0.5*(I+sigma_z)
	//Bmatrices.setElement(sign*ONE_c,Indices(0,0,6));
	Bmatrices.setElement(sign*mu/2.0,0.0,Indices(0,0,6));
	
	//Contract and set the MPO
	Bmatrices.reshape(Indices(D_bond,7));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,1,2,2));
	hamil.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else 
      {
	cout << "Fermonic site, sign = " << sign << endl << endl;
	mwArray Bmatrices(Indices(D_bond,D_bond,7));
	Bmatrices.fillWithZero();
	
	//Identity
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	Bmatrices.setElement(ONE_c,Indices(5,5,0));
	
	//sigma_+
	Bmatrices.setElement(ONE_c,Indices(0,1,1));
	Bmatrices.setElement(ONE_c,Indices(4,5,1));
	
	//sigma_-
	Bmatrices.setElement(ONE_c,Indices(0,2,2));
	Bmatrices.setElement(ONE_c,Indices(3,5,2));	
	
	//0.5*(I+sigma_z)
	//Bmatrices.setElement(sign*ONE_c,Indices(0,5,6));
	Bmatrices.setElement(sign*mu/2.0,0.0,Indices(0,5,6));
	
	//Contract and set the MPO
	Bmatrices.reshape(Indices(D_bond*D_bond,7));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,D_bond,2,2));
	hamil.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }      
    }
  }
  //hamil.exportForMatlab("HSngi.m");
}

void SchwingerHamiltonianNonGaugeinvariant::initChargeoperator(void)
{
  cout << "Constructing the MPO representation of the Charge Operator" << endl;
  
  //Sigma_z for the fermionic sites
  mwArray id_fermionic = identityMatrix(2);
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1.0*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  //Identity for the bosonic sites
  mwArray id_bosonic = identityMatrix(this->d);
  
  //Reshape operators for contraction later 
  
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
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(2,d,d));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<d; ++i)
  {
    for(int j=0; j<d; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
    }
  }
  Bosonic_operators.reshape(Indices(2,d*d));
  
  //Bond dimension for the MPO Matrices
  int D_bond = 2;
  
  for(int k=0; k<N_total;++k)
  {
    cout << "k= " << k << endl;
    if((k%2)==1)
    {
      //Bosonic site(odd sites)
      mwArray Bmatrices(Indices(D_bond,D_bond,2));
      
      //Identity
      Bmatrices.setElement(ONE_c,Indices(0,0,0));
      Bmatrices.setElement(ONE_c,Indices(1,1,0));
      
      Bmatrices.reshape(Indices(D_bond*D_bond,2));
      mwArray res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond,D_bond,d,d));
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
  //charge.exportForMatlab("ChargeOp.m");
}

void SchwingerHamiltonianNonGaugeinvariant::constructLopMPO(MPO &mpo, int n)
{
  cout << endl <<"Constructing MPO for L("<<n<<")" << endl;
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
  
  //Provide Fermionic Identiy
  mwArray id_fermionic = identityMatrix(2);
  
  //Provide Bosonic Operators
  mwArray L_operator(Indices(this->d,this->d));
  mwArray id_bosonic = identityMatrix(this->d);
  
  L_operator.fillWithZero();
  
  for(int i=0; i<this->d; ++i){
    for(int j=0; j<this->d; ++j){
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
  mwArray Fermionic_operators(Indices(2,2,2));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<2; ++i)
  {
    for(int j=0; j<2; ++j)
    {
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(2,2*2));
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(2,d,d));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<d; ++i)
  {
    for(int j=0; j<d; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
      Bosonic_operators.setElement(L_operator.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Bosonic_operators.reshape(Indices(2,d*d));
  
  //Bond dimension for the MPO Matrices
  int D_bond = 2;
  
  for(int k=0; k<N_total;++k)
  {
    if((k%2)==1)
    {
      //Bosonic site(odd sites)       
      if(k<n)
      {
	 //Sites left of site n
	 //cout << "Bosonic site " << k << " left of " << n << endl;
	 mwArray Bmatrices(Indices(D_bond,D_bond,2));
	 Bmatrices.fillWithZero();
	
	//Identity
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	
	Bmatrices.reshape(Indices(D_bond*D_bond,2));
	mwArray res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond,D_bond,d,d));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else if(k>n)
      {
	//Sites right of site n
	//cout << "Bosonic site " << k << " right of " << n << endl;
	mwArray Bmatrices(Indices(D_bond,D_bond,2));
	Bmatrices.fillWithZero();
	
	//Identity
	Bmatrices.setElement(ONE_c,Indices(1,1,0));
	
	Bmatrices.reshape(Indices(D_bond*D_bond,2));
	mwArray res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond,D_bond,d,d));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else if(k==n)
      {
	//Site whre the operator has to be placed
	//cout << "Bosonic site " << k << " equal to " << n << endl;
	mwArray Bmatrices(Indices(D_bond,D_bond,2));
	Bmatrices.fillWithZero();
	
	//Put the operator L
	Bmatrices.setElement(ONE_c,Indices(0,1,1));
	
	Bmatrices.reshape(Indices(D_bond*D_bond,2));
	mwArray res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond,D_bond,d,d));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }	
    }
    else
    {
      //Fermionic site(even sites)
      if(k==0)
      {
	//First site
	//cout << "First Fermionic site " << k << " left of " << n << endl;
	mwArray Bmatrices(Indices(1,D_bond,2));
	Bmatrices.fillWithZero();
	
	//Identiy
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	
	Bmatrices.reshape(Indices(D_bond,2));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(1,D_bond,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	
      }
      else if(k==(N_total-1))
      {
	//Last site
	//cout << "Last Fermionic site " << k << " right of " << n << endl;
	mwArray Bmatrices(Indices(D_bond,1,2));
	Bmatrices.fillWithZero();
	
	//Identiy
	Bmatrices.setElement(ONE_c,Indices(1,0,0));	
	
	Bmatrices.reshape(Indices(D_bond,2));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,1,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else if(k<n)
      {
	//Sites left of n
	//cout << "Fermionic site " << k << " left of " << n << endl;
	mwArray Bmatrices(Indices(D_bond,D_bond,2));
	Bmatrices.fillWithZero();
	
	//Identiy
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	
	Bmatrices.reshape(Indices(D_bond*D_bond,2));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,D_bond,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else if(k>n)
      {
	//Sites left of n
	//cout << "Fermionic site " << k << " right of " << n << endl;
	mwArray Bmatrices(Indices(D_bond,D_bond,2));
	Bmatrices.fillWithZero();
	
	//Identiy
	Bmatrices.setElement(ONE_c,Indices(1,1,0));
	
	Bmatrices.reshape(Indices(D_bond*D_bond,2));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,D_bond,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }     
    }
    
  }
  string Result;  
  ostringstream convert; 
  //Add one to be compliant with Matlab which starts counting by 1
//   convert << n+1; 
//   Result = convert.str();  
//   string name="LnOp"+Result+".m";
//   mpo.exportForMatlab(name.c_str());  
}


void SchwingerHamiltonianNonGaugeinvariant::constructGaussMPO(MPO &mpo, int n)
{
  cout << endl <<"Constructing MPO for 0.5*(sigma_z+(-1)^"<<n<<"]" << endl;
  
  //Make sure MPO is empty
  mpo.initLength(N_total);  
  
  //Provide Fermionic Identiy
  mwArray id_fermionic = identityMatrix(2);
  
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  //Provide Bosonic Operators
  mwArray id_bosonic = identityMatrix(this->d);
  
  if((n%2)!=0)
    cout << "!!! Error, no fermionic site !!!"<< endl;
  
  int sign=pow(-1,n/2+1);
  
  sigma_z = 0.5*(sigma_z+sign*id_fermionic);
  
  
  //Reshape operators for contraction later 
  
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
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(2,d,d));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<d; ++i)
  {
    for(int j=0; j<d; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
    }
  }
  Bosonic_operators.reshape(Indices(2,d*d));
  
  //Bond dimension for the MPO Matrices
  int D_bond = 2;
  
  for(int k=0; k<N_total;++k)
  {
    if((k%2)==1)
    {
      if(k<n)
      {
	//Bosonic sites
	//cout << "Bosonic site " << k <<" left of " << n << endl;
	mwArray Bmatrices(Indices(D_bond,D_bond,2));
	Bmatrices.fillWithZero();
	  
	//Identity
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	  
	Bmatrices.reshape(Indices(D_bond*D_bond,2));
	mwArray res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond,D_bond,d,d));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else
      {
	//Bosonic sites
	//cout << "Bosonic site " << k << " right of " << endl;
	mwArray Bmatrices(Indices(D_bond,D_bond,2));
	Bmatrices.fillWithZero();
	  
	//Identity
	Bmatrices.setElement(ONE_c,Indices(1,1,0));
	  
	Bmatrices.reshape(Indices(D_bond*D_bond,2));
	mwArray res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond,D_bond,d,d));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
    }
    else
    {
      //Fermionic site(even sites)
      if(k==n)
      {
	//Site where the operator is located
	if((k!=0) && (k!=(N_total-1)))
	{
	  //Site in between
	  //cout << "Fermionic site " << k << " equal to " << n << endl;
	  mwArray Bmatrices(Indices(D_bond,D_bond,2));
	  Bmatrices.fillWithZero();
	  
	  //0.5*(sigma_z+(-1)^n)
	  Bmatrices.setElement(ONE_c,Indices(0,1,1));
	  
	  Bmatrices.reshape(Indices(D_bond*D_bond,2));
	  mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,D_bond,2,2));
	  mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	}
	else if(k==0)
	{
	  //First site
	  //cout << "First Fermionic site " << k << " equal to " << n << endl;
	  mwArray Bmatrices(Indices(1,D_bond,2));
	  Bmatrices.fillWithZero();
	  
	  //0.5*(sigma_z+(-1)^n)
	  Bmatrices.setElement(ONE_c,Indices(0,1,1));
	  
	  Bmatrices.reshape(Indices(D_bond,2));
	  mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(1,D_bond,2,2));
	  mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	}
	else if(k==(N_total-1))
	{
	  //Last site
	  //cout << "Last Fermionic site " << k << " equal to " << n << endl;
	  mwArray Bmatrices(Indices(D_bond,1,2));
	  Bmatrices.fillWithZero();
	  
	  //0.5*(sigma_z+(-1)^n)
	  Bmatrices.setElement(ONE_c,Indices(0,0,1));
	  
	  Bmatrices.reshape(Indices(D_bond,2));
	  mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,1,2,2));
	  mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	}	
      }
      else
      {
	//Sites whith no operator
	if(k<n)
	{
	  if(k==0)
	  {
	    //Site in between
	    //cout << "Fermionic site " << k << " left of " << n << endl;
	    mwArray Bmatrices(Indices(1,D_bond,2));
	    Bmatrices.fillWithZero();
	    
	    //0.5*(sigma_z+(-1)^n)
	    Bmatrices.setElement(ONE_c,Indices(0,0,0));
	    
	    Bmatrices.reshape(Indices(D_bond,2));
	    mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(1,D_bond,2,2));
	    mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	  }
	  else
	  {
	    //Site in between
	    //cout << "Fermionic site " << k << " left of " << n << endl;
	    mwArray Bmatrices(Indices(D_bond,D_bond,2));
	    Bmatrices.fillWithZero();
	    
	    //0.5*(sigma_z+(-1)^n)
	    Bmatrices.setElement(ONE_c,Indices(0,0,0));
	    
	    Bmatrices.reshape(Indices(D_bond*D_bond,2));
	    mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,D_bond,2,2));
	    mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	  }
	}
	else
	{
	  if(k==(N_total-1))
	  {
	    //Site in between
	    //cout << "Last Fermionic site " << k << " right of " << n << endl;
	    mwArray Bmatrices(Indices(D_bond,1,2));
	    Bmatrices.fillWithZero();
	    
	    //0.5*(sigma_z+(-1)^n)
	    Bmatrices.setElement(ONE_c,Indices(1,0,0));
	    
	    Bmatrices.reshape(Indices(D_bond,2));
	    mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,1,2,2));
	    mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	  }
	  else
	  {
	    //Site in between
	    //cout << "Fermionic site " << k << " right of " << n << endl;
	    mwArray Bmatrices(Indices(D_bond,D_bond,2));
	    Bmatrices.fillWithZero();
	    
	    //Identitaet
	    Bmatrices.setElement(ONE_c,Indices(1,1,0));
	    
	    Bmatrices.reshape(Indices(D_bond*D_bond,2));
	    mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond,D_bond,2,2));
	    mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	  }
	}	
      }
    }    
  }
/*  string Result;  
  ostringstream convert; 
  convert << n+1; 
  Result = convert.str();  
  string name="GaussOp"+Result+".m";
  mpo.exportForMatlab(name.c_str()); */ 
}

void SchwingerHamiltonianNonGaugeinvariant::constructsigmazMPO(MPO &mpo, int n)
{
  cout << endl <<  "Constructing MPO for sigma_z("<<n<<")" << endl;
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
  //Sigma_z for the fermionic sites
  mwArray id_fermionic = identityMatrix(2);
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1.0*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  //Identity for the bosonic sites
  mwArray id_bosonic = identityMatrix(this->d);
  
  //Reshape operators for contraction later 
  
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
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(2,d,d));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<d; ++i)
  {
    for(int j=0; j<d; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
    }
  }
  Bosonic_operators.reshape(Indices(2,d*d));
  
  //Bond dimension for the MPO Matrices
  int D_bond_left;
  int D_bond_right;
  
  for(int k=0; k<N_total;++k)
  {
    //cout << "k= " << k << endl;
    D_bond_left = 2;
    D_bond_right = 2;
    if(k<n)
    {
      if(k==0)
	D_bond_left = 1;
      
      mwArray Bmatrices(Indices(D_bond_left,D_bond_right,2));      
      Bmatrices.fillWithZero();
      //Just put the identiy
      Bmatrices.setElement(ONE_c,Indices(0,0,0));
      
      if((k%2)==0)
      {
	//Fermionic site
	Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true); 
      }
      else if((k%2)==1)
      {
	//Bosonic site
	Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));
	mwArray res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond_left,D_bond_right,d,d));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true); 
      }
      else 
	cout << "Error constructing the sigma_z("<<n<<") MPO at site "<< k << endl;     
    }    
    else if(k>n)
    {
      //Here last site needs extra treatment, since the identities are placed at index (1,1) but the last site is a vector and there it has to be (0,1)
      if(k==(N_total-1))
      {
	D_bond_right = 1;
	mwArray Bmatrices(Indices(D_bond_left,D_bond_right,2));      
	Bmatrices.fillWithZero();
	//Just put the identiy
	Bmatrices.setElement(ONE_c,Indices(1,0,0));
	if((k%2)==0)
	{
	  //Fermionic site
	  Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));
	  mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	  mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true); 
	}
	else if((k%2)==1)
	{
	  //Bosonic site
	  Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));
	  mwArray res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond_left,D_bond_right,d,d));
	  mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true); 
	}
	else 
	  cout << "Error constructing the sigma_z("<<n<<") MPO at site "<< k << endl;     
      }
      else
      {
	mwArray Bmatrices(Indices(D_bond_left,D_bond_right,2));      
	Bmatrices.fillWithZero();
	//Just put the identiy
	Bmatrices.setElement(ONE_c,Indices(1,1,0));
	if((k%2)==0)
	{
	  //Fermionic site
	  Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));
	  mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	  mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true); 
	}
	else if((k%2)==1)
	{
	  //Bosonic site
	  Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));
	  mwArray res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond_left,D_bond_right,d,d));
	  mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true); 
	}
	else 
	  cout << "Error constructing the sigma_z("<<n<<") MPO at site "<< k << endl;     
      }
    }
    else
    {
      if(k==0)
      {
	//cout << "Put operator at start" << endl;
	D_bond_left = 1;
	mwArray Bmatrices(Indices(D_bond_left,D_bond_right,2));      
	Bmatrices.fillWithZero();
	//Just put the sigma_z
	Bmatrices.setElement(ONE_c,Indices(0,1,1));
	Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else if(k==(N_total-1))
      {
	//cout << "Sigmaz at end" << endl;
	D_bond_right = 1;
	mwArray Bmatrices(Indices(D_bond_left,D_bond_right,2));      
	Bmatrices.fillWithZero();
	//Just put the sigma_z
	Bmatrices.setElement(ONE_c,Indices(0,0,1));
	Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else
      {
	mwArray Bmatrices(Indices(D_bond_left,D_bond_right,2));      
	Bmatrices.fillWithZero();
	//Just put the sigma_z
	Bmatrices.setElement(ONE_c,Indices(0,1,1));
	Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));
	mwArray res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      
    }
  }
  //mpo.exportForMatlab("sigma_zn.m");
}

void SchwingerHamiltonianNonGaugeinvariant::constructCondensateMPO(MPO &mpo)
{
  cout << "Constructing condensate MPO"<< endl;
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  mwArray condesateOp = id_fermionic + sigma_z;
  
  //Provide bosonic identiy
  mwArray id_bosonic = identityMatrix(this->d);
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,2,2));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<2; ++i)
  {
    for(int j=0; j<2; ++j)
    {
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
      Fermionic_operators.setElement(condesateOp.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(2,2*2));
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(2,this->d,this->d));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<this->d; ++i)
  {
    for(int j=0; j<this->d; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
    }
  }
  Bosonic_operators.reshape(Indices(2,this->d*this->d));
  
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
	//First fermionic site
	mwArray Bmatrices(Indices(1,D_bond_right,2));
	Bmatrices.fillWithZero();
	
	//Identity
	Bmatrices.setElement(ONE_c,Indices(0,0,0));	
	//Condensate operator
	Bmatrices.setElement(0.5*sign*sqrt(this->x)/this->N_fermions,0.0,Indices(0,1,1));
	
	//Contract and set the operator
	Bmatrices.reshape(Indices(D_bond_right,2));
	res=reshape(Bmatrices*Fermionic_operators,Indices(1,D_bond_right,2,2));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	mpo.setOp(k,Local_op,true); 
      }
      else if(k==(N_total-1))
      {
	//Last fermionic site
	mwArray Bmatrices(Indices(D_bond_left,1,2));
	Bmatrices.fillWithZero();
	
	//Identity
	Bmatrices.setElement(ONE_c,Indices(1,0,0));	
	//Condensate operator
	Bmatrices.setElement(0.5*sign*sqrt(this->x)/this->N_fermions,0.0,Indices(0,0,1));
	
	//Contract and set the operator
	Bmatrices.reshape(Indices(D_bond_left,2));
	res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond_left,1,2,2));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	mpo.setOp(k,Local_op,true); 
      }
      else
      {
	//Sites in between
	mwArray Bmatrices(Indices(D_bond_left,D_bond_right,2));
	Bmatrices.fillWithZero();
	
	//Identity
	Bmatrices.setElement(ONE_c,Indices(0,0,0));
	Bmatrices.setElement(ONE_c,Indices(1,1,0));	
	//Condensate operator
	Bmatrices.setElement(0.5*sign*sqrt(this->x)/this->N_fermions,0.0,Indices(0,1,1));
	
	//Contract and set the operator
	Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));	res=reshape(Bmatrices*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	mpo.setOp(k,Local_op,true); 
      }
    }
    else
    {
      //Bosonic site
      mwArray Bmatrices(Indices(D_bond_left,D_bond_right,2));
      Bmatrices.fillWithZero();
      
      //Identity
      Bmatrices.setElement(ONE_c,Indices(0,0,0));
      Bmatrices.setElement(ONE_c,Indices(1,1,0));
      
      //Contract and set the operator
      Bmatrices.reshape(Indices(D_bond_left*D_bond_right,2));
      res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d,this->d));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(k,Local_op,true); 
    }
  }
  
}