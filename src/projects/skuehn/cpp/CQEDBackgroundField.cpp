#include "CQEDBackgroundField.h"
#include "Indices.h"
#include "Contractor.h"

using namespace std;
using namespace shrt;


CQEDBackgroundField::CQEDBackgroundField(int N_, int d_,int D_,double mu_,double x_,double e_,double lambda_):hamil(2*N_-1),charge(2*N_-1)
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
  terms_saved = false;
  
  
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
CQEDBackgroundField::~CQEDBackgroundField(){
    hamil.clear();
    charge.clear();    
}

void CQEDBackgroundField::initOperators()
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
  mwArray Op1(Indices(this->d_bose+1,this->d_bose+1));
  mwArray Op2(Indices(this->d_bose+1,this->d_bose+1));
  mwArray L_operator(Indices(this->d_bose+1,this->d_bose+1));
  mwArray Lsq_operator(Indices(this->d_bose+1,this->d_bose+1));
  mwArray id_bosonic = identityMatrix(this->d_bose+1);
  
  L_operator.fillWithZero();
  Lsq_operator.fillWithZero();
  Op1.fillWithZero();
  Op2.fillWithZero();
  
  for(int i=0; i<=this->d_bose; ++i)
  {
    L_operator.setElement(-0.5*this->d_bose+i,0.0,Indices(i,i));
  }
  
  //Construct Op1 und Op2
  //Since I want to reuse my old implementation which has hopping terms like sp Op1 sm + sm op2 sp I absorb i in Op1 and -i in Op2
  double norm_factor = sqrt(this->d_bose/2*(this->d_bose/2+1));
  complex_t tmp;
  for(int i=0; i<=this->d_bose; ++i)
  {
    for(int j=0; j<=this->d_bose; ++j)
    {
      //Lp
      if(i==j+1)
      {
	Op1.setElement(0.0,sqrt(i*(this->d_bose-j))/norm_factor,Indices(i,j));
      }
//       //-1.0*Lm
//       else if(i==j-1)
//       {
// 	Op2.setElement(0.0,-sqrt(i*(this->d_bose-j-2))/norm_factor,Indices(i,j));
//       }
    }
  }
  
  //Get Op2 as the hermitian conjugate of Op1
  Op2=Op1;
  Op2.transpose(true);
  //Op2.transpose();
  
  /*cout << "L = " << L_operator << endl;
  cout << "Op1 = " << Op1 << endl;
  cout << "Op2 = " << Op2 << endl;*/
  
  //Now get get the square of the L-operator
  Lsq_operator = L_operator*L_operator;
  
  /*cout << "L²=" << Lsq_operator << endl;
  cout << "Theta = " << Theta << endl;
  cout << "Op1 = " << Op1 << endl;
  cout << "Op2 = " << Op2 << endl;*/
  
   //Array to store the bosonic operators
   mwArray Bosonic_operators(Indices(8,this->d_bose+1,this->d_bose+1));
   Bosonic_operators.fillWithZero();
   
   for(int i=0; i<=this->d_bose; ++i)
   {
     for(int j=0; j<=this->d_bose; ++j)
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
   
  Bosonic_operators.reshape(Indices(8,(this->d_bose+1)*(this->d_bose+1)));  
  
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
  
  //cout << "Constructing Hamiltonian MPO" << endl;
  
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
      res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose+1,this->d_bose+1));
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
  //hamil.exportForMatlab("HCQEDBackgroundField.m");
}// initOperators

void CQEDBackgroundField::initChargeoperator()
{
  cout << "Constructing the MPO representation of the Charge Operator" << endl;
  
  //Sigma_z for the fermionic sites
  mwArray id_fermionic = identityMatrix(2);
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1.0*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  //Identity for the bosonic sites
  mwArray id_bosonic = identityMatrix(this->d_bose+1);
  
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
  mwArray Bosonic_operators(Indices(2,this->d_bose+1,this->d_bose+1));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<=this->d_bose; ++i)
  {
    for(int j=0; j<=this->d_bose; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
    }
  }
  Bosonic_operators.reshape(Indices(2,(this->d_bose+1)*(this->d_bose+1)));
  
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
      mwArray res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose+1,this->d_bose+1));
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

void CQEDBackgroundField::constructLopMPO(MPO &mpo, int n)
{
  //cout << "Constructing L("<<n<<")" << endl;
  
  //As n is a DMRG site and L(n) acts on bosonic sites only, check if n is a bosonic index
  if((n%2)!=1)
  {
    cerr << "Error in CQEDBackgroundField::constructLopMPO(MPO &mpo, int n): Index " << n << " is not a valid bosonic site" << endl;
    exit(1);
  }
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
  //Provide Fermionic Identiy
  mwArray id_fermionic = identityMatrix(2);
  
  //Provide Bosonic Operators
  mwArray L_operator(Indices(this->d_bose+1,this->d_bose+1));
  mwArray id_bosonic = identityMatrix(this->d_bose+1);
  
  L_operator.fillWithZero();
  
  for(int i=0; i<=this->d_bose; ++i)
  {
    L_operator.setElement(-0.5*this->d_bose+i,0.0,Indices(i,i));
  } 
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(1,2,2));
  Fermionic_operators.fillWithZero();  
 
  for(int i=0; i<this->d_fermi; ++i)
    for(int j=0; j<this->d_fermi; ++j)
      Fermionic_operators.setElement(id_fermionic.getElement(Indices(i,j)),Indices(0,i,j));
      
  Fermionic_operators.reshape(Indices(1,this->d_fermi*this->d_fermi));
 
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(2,this->d_bose+1,this->d_bose+1));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<=this->d_bose; ++i)
  {
    for(int j=0; j<=this->d_bose; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
      Bosonic_operators.setElement(L_operator.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Bosonic_operators.reshape(Indices(2,(this->d_bose+1)*(this->d_bose+1)));
  
  //Bond dimension for the MPO Matrices
  int D_bond_left = 1;
  int D_bond_right = 1;
  
  //Different matrices appearing
  mwArray Fermi_identity(Indices(D_bond_left,D_bond_right,1));
  mwArray Bose_identity(Indices(D_bond_left,D_bond_right,2));
  mwArray Bose_Op(Indices(D_bond_left,D_bond_right,2));
  
  //Fill with zeros
  Fermi_identity.fillWithZero();
  Bose_identity.fillWithZero();
  Bose_Op.fillWithZero();
  
  //Set elements
  Fermi_identity.setElement(ONE_c,Indices(0,0,0));
  
  Bose_identity.setElement(ONE_c,Indices(0,0,0));
  Bose_Op.setElement(ONE_c,Indices(0,0,1));
  
  //reshape
  Fermi_identity.reshape(Indices(D_bond_right,1));
    
  Bose_identity.reshape(Indices(D_bond_left*D_bond_right,2));
  Bose_Op.reshape(Indices(D_bond_left*D_bond_right,2));
  
  //Array to store temporary terms
  mwArray res;
  
  //Flags for memory
  bool bose_saved=false;
  bool fermi_saved=false;
  int first_saved_site;
  
  //Set elements
  for(int k=0; k<this->N_total; k++)
  {
    res.clear();
    if((k%2)==1)
    {
		if(k==n)
		{
			//Site where the operator is located
			res=reshape(Bose_Op*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose+1,this->d_bose+1));
			mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true); 
		}
		else
		{
			//Just put the bosonic identity
			if(!bose_saved)
			{
				res=reshape(Bose_identity*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose+1,this->d_bose+1));
				mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true); 
				bose_saved=true;
				first_saved_site=k;
			}
			else
				mpo.setOp(k,&mpo.getOp(first_saved_site),false); 	
		}
	}
	else
	{
		//Fermionic site, just set identity
		if(!fermi_saved)
		{
			res=reshape(Fermi_identity*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
			mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true); 
		}
		else
			mpo.setOp(k,&mpo.getOp(k-2),false); 
	}	
  }
}

void CQEDBackgroundField::constructGaussMPO(MPO &mpo, int n)
{
  //cout << "Constructing GaussMPO for site " << n << endl;
  
  //As n is a DMRG site and sigma_z acts on fermionic sites only, check if n is a bosonic index
  if((n%2)==1)
  {
    cerr << "Error in CQEDBackgroundField::constructGaussMPO(MPO &mpo, int n): Index " << n << " is not a valid fermionic site" << endl;
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
  mwArray id_bosonic = identityMatrix(this->d_bose+1);
  
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
  mwArray Bosonic_operators(Indices(1,this->d_bose+1,this->d_bose+1));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<=this->d_bose; ++i)
    for(int j=0; j<=this->d_bose; ++j)
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
      
  Bosonic_operators.reshape(Indices(1,(this->d_bose+1)*(this->d_bose+1)));
  
  //Bond dimension for the MPO Matrices
  int D_bond_left = 1;
  int D_bond_right = 1;
  
  //Different matrices appearing
  mwArray Fermi_identity(Indices(D_bond_left,D_bond_right,2));
  mwArray Fermi_sigma_z(Indices(D_bond_left,D_bond_right,2));
  mwArray Bose_identity(Indices(D_bond_left,D_bond_right,1));
  
  //Fill with zeros
  Fermi_identity.fillWithZero();
  Fermi_sigma_z.fillWithZero();
  Bose_identity.fillWithZero();
  
  //Set elements
  Fermi_identity.setElement(ONE_c,Indices(0,0,0));
  Fermi_sigma_z.setElement(ONE_c,Indices(0,0,1));
  
  Bose_identity.setElement(ONE_c,Indices(0,0,0));
  
  //reshape  
  Bose_identity.reshape(Indices(D_bond_left*D_bond_right,1));
  
  Fermi_identity.reshape(Indices(D_bond_left*D_bond_right,2));
  Fermi_sigma_z.reshape(Indices(D_bond_left*D_bond_right,2));
  
  //Array for the results
  mwArray res;  
  
  //Flags for memory
  bool bose_saved=false;
  bool fermi_saved=false;
  int first_saved_site;
  
  //Set the operators
  for(int k=0;k<this->N_total; k++)
  {
    if((k%2)==1)
    {
	  //Bosonic site, just put identity
	  if(!bose_saved)
	  {
		res=reshape(Bose_identity*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose+1,this->d_bose+1));
		mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
		bose_saved=true;
	  }
	  else
		mpo.setOp(k,&mpo.getOp(k-2),false);
	}
    else 
    {
	  if(k==n)
	  {
		  //Site where the operator is located
		  res=reshape(Fermi_sigma_z*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
		  mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	  }
	  else
	  {
		  //Fermionic sites with identity
		  if(!fermi_saved)
		  {
			res=reshape(Fermi_identity*Fermionic_operators,Indices(D_bond_left,D_bond_right,this->d_fermi,this->d_fermi));
			mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
			fermi_saved=true;
			first_saved_site=k;
		 }
		 else
			mpo.setOp(k,&mpo.getOp(first_saved_site),false);		
	  }
    }
  }
}

void CQEDBackgroundField::constructsigmazMPO(MPO &mpo, int n)
{
  //cout << "Constructing sigma_z(" << n << ")" << endl;
   
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
 
  //Provide bosonic operators
  mwArray id_bosonic = identityMatrix(this->d_bose+1);
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(1,this->d_bose+1,this->d_bose+1));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<=this->d_bose; ++i)
    for(int j=0; j<=this->d_bose; ++j)
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));
  
  Bosonic_operators.reshape(Indices(1,(this->d_bose+1)*(this->d_bose+1)));  
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,this->d_fermi,this->d_fermi));
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
  Fermionic_operators.reshape(Indices(2,this->d_fermi*this->d_fermi));
  
  //Bond dimension for the MPO Matrices
  int D_bond_left = 1;
  int D_bond_right = 1;
  
  //Different matrices appearing
  mwArray Fermi_identity(Indices(D_bond_left,D_bond_right,2));
  mwArray Fermi_sigma_z(Indices(D_bond_left,D_bond_right,2));
  
  mwArray Bose_identity(Indices(D_bond_left,D_bond_right,1));
  
  //Fill with zeros
  Fermi_identity.fillWithZero();
  Fermi_sigma_z.fillWithZero();
  Bose_identity.fillWithZero();
  
  //Set entries
  Fermi_identity.setElement(ONE_c,Indices(0,0,0));
  Fermi_sigma_z.setElement(ONE_c,Indices(0,0,1));
  
  Bose_identity.setElement(ONE_c,Indices(0,0,0));
  
  //reshape
  Fermi_identity.reshape(Indices(D_bond_left*D_bond_right,2));
  Fermi_sigma_z.reshape(Indices(D_bond_left*D_bond_right,2));
    
  Bose_identity.reshape(Indices(D_bond_left*D_bond_right,1));
  
  //Make sure MPO is empty before filling
  mpo.initLength(this->N_total);
  
  //Temporary array for MPO and stuff
  mwArray res;
  
  //Flags for memory
  bool bose_saved=false;
  bool fermi_saved=false;
  int first_saved_site;
 
  
  for(int k=0; k<this->N_total; k++)
  {
    if((k%2)==1)
    {
      //Bosonic site where only identity can reside
      if(!bose_saved)
      {
		  res=reshape(Bose_identity*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose+1,this->d_bose+1));
		  mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
		  bose_saved = true;
	  }
	  else
		  mpo.setOp(k,&mpo.getOp(k-2),false);
    }
    else
    {
      //Fermionic site
      if(k==n)
      {
		//Site at which the operator starts
		res=reshape(Fermi_sigma_z*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
		mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else
      {
		//Just set fermionic identity
		if(!fermi_saved)
		{
			res=reshape(Fermi_identity*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
			mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
			fermi_saved = true;
			first_saved_site = k;
		}
		else
			mpo.setOp(k,&mpo.getOp(first_saved_site),false);
			
      }
    }	
  } 
}


void CQEDBackgroundField::constructCondensateMPO(MPO &mpo)
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
  mwArray id_bosonic = identityMatrix(this->d_bose+1);
  
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
  mwArray Bosonic_operators(Indices(2,this->d_bose+1,this->d_bose+1));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<this->d_bose; ++i)
  {
    for(int j=0; j<this->d_bose; ++j)
    {
      Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));      
    }
  }
  Bosonic_operators.reshape(Indices(2,(this->d_bose+1)*(this->d_bose+1)));
  
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
      res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose+1,this->d_bose+1));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(k,Local_op,true); 
    }
  }  
  //cout << "Writing MPO" << endl;
  //cout << mpo << endl;
  //mpo.exportForMatlab("CondensateMPO.m");
}//CondensateMPO


//Construct an explicit MPO for the penalty term to be able to see its contribution 
void CQEDBackgroundField::constructPenaltyMPO(MPO &mpo)
{
  cout << "Constructing Penalty MPO" << endl;
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
 
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sigma_z(Indices(2,2),data_z);
  
  //Provide Operators for gauge field
  mwArray L_operator(Indices(this->d_bose+1,this->d_bose+1));
  mwArray Lsq_operator(Indices(this->d_bose+1,this->d_bose+1));
  mwArray id_bosonic = identityMatrix(this->d_bose+1);
  
  Lsq_operator.fillWithZero();
  L_operator.fillWithZero();
  
  for(int i=0; i<=this->d_bose; ++i)
  {
    L_operator.setElement(-0.5*this->d_bose+i,0.0,Indices(i,i));
  }
  
  //Now get get the square of the L-operator
  Lsq_operator = L_operator*L_operator;
  //cout << "Lop" << L_operator << endl;
  
   //Array to store the bosonic operators
   mwArray Bosonic_operators(Indices(4,this->d_bose+1,this->d_bose+1));
   Bosonic_operators.fillWithZero();
   
   for(int i=0; i<=this->d_bose; ++i)
   {
     for(int j=0; j<=this->d_bose; ++j)
     {
       //Identity
       Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));
       //L
       Bosonic_operators.setElement(L_operator.getElement(Indices(i,j)),Indices(2,i,j));
       //Lsquared
       Bosonic_operators.setElement(Lsq_operator.getElement(Indices(i,j)),Indices(3,i,j));
     }
   }
   
   Bosonic_operators.reshape(Indices(4,(this->d_bose+1)*(this->d_bose+1)));  
  
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
      res=reshape(Bmatrices*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose+1,this->d_bose+1));
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

void CQEDBackgroundField::constructEvolutionMPO(MPO &Uodd, MPO &Ueven, double dt)
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
  mwArray Op1(Indices(this->d_bose+1,this->d_bose+1));
  mwArray Op2(Indices(this->d_bose+1,this->d_bose+1));
  mwArray L_operator(Indices(this->d_bose+1,this->d_bose+1));
  mwArray Lsq_operator(Indices(this->d_bose+1,this->d_bose+1));
  mwArray id_bosonic = identityMatrix(this->d_bose+1);
  
  L_operator.fillWithZero();
  Lsq_operator.fillWithZero();
  Op1.fillWithZero();
  Op2.fillWithZero();
  
  for(int i=0; i<=this->d_bose; ++i)
  {
    L_operator.setElement(-0.5*this->d_bose+i,0.0,Indices(i,i));
  }
  
  //Now get L_squared
  Lsq_operator = L_operator*L_operator;
  
  //Construct Op1 und Op2 (imaginary unit including its sign is absorbed in the operators)
  double norm_factor = sqrt(this->d_bose/2*(this->d_bose/2+1));
  complex_t tmp;
  for(int i=0; i<=this->d_bose; ++i)
  {
    for(int j=0; j<=this->d_bose; ++j)
    {
      //Lp
      if(i==j+1)
      {
	Op1.setElement(0.0,sqrt(i*(this->d_bose-j))/norm_factor,Indices(i,j));
      }
    }
  }
  
  //Get Lm aut of Lp by simply taking the hermitian conjugate (and therefore getting the right sign!)
  Op2=Op1;  
  Op2.transpose(true);
  
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
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  mwArray U2,S2,V2;
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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
  exp_even.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_even.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_even.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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
  id_bose_mpo_odd.reshape(Indices(this->d_bose+1,1,this->d_bose+1,1));
  id_bose_mpo_even.reshape(Indices(this->d_bose+1,1,this->d_bose+1,1));
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

MPS* CQEDBackgroundField::constructInitalMPS(InitState state, int extent)
{  
  //Up to know the MPS has to have the right dimensions and so on
  //initmps.clear();
  mwArray spinsiteup(Indices(2,1,1));
  mwArray spinsitedown(Indices(2,1,1));
  mwArray bosesite(Indices(this->d_bose+1,1,1));
  mwArray bosesiteplus(Indices(this->d_bose+1,1,1));
  mwArray bosesiteminus(Indices(this->d_bose+1,1,1));
  
  spinsiteup.fillWithZero();
  spinsitedown.fillWithZero();
  bosesite.fillWithZero();
  bosesiteplus.fillWithZero();
  bosesiteminus.fillWithZero();
  
  spinsiteup.setElement(ONE_c,Indices(0,0,0));
  spinsitedown.setElement(ONE_c,Indices(1,0,0));
  bosesite.setElement(ONE_c,Indices(this->d_bose/2,0,0));
  bosesiteplus.setElement(ONE_c,Indices(this->d_bose/2+1,0,0));
  bosesiteminus.setElement(ONE_c,Indices(this->d_bose/2-1,0,0));
  
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
	phys_dims.push_back(this->d_bose+1);
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

void CQEDBackgroundField::constructEvolutionMPOefficient(MPO &Uodd, MPO &Ueven, double dt)
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
    mwArray Op1(Indices(this->d_bose+1,this->d_bose+1));
    mwArray Op2(Indices(this->d_bose+1,this->d_bose+1));
    mwArray L_operator(Indices(this->d_bose+1,this->d_bose+1));
    mwArray Lsq_operator(Indices(this->d_bose+1,this->d_bose+1));
    mwArray id_bosonic = identityMatrix(this->d_bose+1);
    
    L_operator.fillWithZero();
    Lsq_operator.fillWithZero();
    Op1.fillWithZero();
    Op2.fillWithZero();
    
    for(int i=0; i<=this->d_bose; ++i)
    {
      L_operator.setElement(-0.5*this->d_bose+i,0.0,Indices(i,i));
    }
    
    //Now get L_squared
    Lsq_operator = L_operator*L_operator;
    
    //Construct Op1 und Op2 (imaginary unit including its sign is absorbed in the operators)
    double norm_factor = sqrt(this->d_bose/2*(this->d_bose/2+1));
    complex_t tmp;
    for(int i=0; i<=this->d_bose; ++i)
    {
      for(int j=0; j<=this->d_bose; ++j)
      {
	//Lp
	if(i==j+1)
	{
	  Op1.setElement(0.0,sqrt(i*(this->d_bose-j))/norm_factor,Indices(i,j));
	}
      }      
    }
    
    //Get Lm aut of Lp by simply taking the hermitian conjugate (and therefore getting the right sign!)
    Op2=Op1;  
    Op2.transpose(true);
  
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
    this->id_bose_mpo.reshape(Indices(this->d_bose+1,1,this->d_bose+1,1));
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
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  mwArray U2,S2,V2;
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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
  wrapper::expm(arg_even,exp_even,-I_c*dt);
  //cout << "Argument exp-even " << arg_even << endl;
  //cout << "Exp-even mit dt "<< dt << ": " << exp_even << endl;
  

  //Since reshape is acting columnwise the indices have to be rearranged in the order (ij)(kl)(mn). This is done by reshaping and permuting
  exp_even.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_even.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_even.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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


void CQEDBackgroundField::constructUoddMPO(MPO &Uodd, double dt)
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
    mwArray Op1(Indices(this->d_bose+1,this->d_bose+1));
    mwArray Op2(Indices(this->d_bose+1,this->d_bose+1));
    mwArray L_operator(Indices(this->d_bose+1,this->d_bose+1));
    mwArray Lsq_operator(Indices(this->d_bose+1,this->d_bose+1));
    mwArray id_bosonic = identityMatrix(this->d_bose+1);
    
    L_operator.fillWithZero();
    Lsq_operator.fillWithZero();
    Op1.fillWithZero();
    Op2.fillWithZero();
    
    for(int i=0; i<=this->d_bose; ++i)
    {
      L_operator.setElement(-0.5*this->d_bose+i,0.0,Indices(i,i));
    }
    
    //Now get L_squared
    Lsq_operator = L_operator*L_operator;
    
    //Construct Op1 und Op2 (imaginary unit including its sign is absorbed in the operators)
    double norm_factor = sqrt(this->d_bose/2*(this->d_bose/2+1));
    complex_t tmp;
    for(int i=0; i<=this->d_bose; ++i)
    {
      for(int j=0; j<=this->d_bose; ++j)
      {
	//Lp
	if(i==j+1)
	{
	  Op1.setElement(0.0,sqrt(i*(this->d_bose-j))/norm_factor,Indices(i,j));
	}
      }      
    }
    
    //Get Lm aut of Lp by simply taking the hermitian conjugate (and therefore getting the right sign!)
    Op2=Op1;  
    Op2.transpose(true);
  
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
    this->id_bose_mpo.reshape(Indices(this->d_bose+1,1,this->d_bose+1,1));
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
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  mwArray U2,S2,V2;
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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

void CQEDBackgroundField::constructUoddMPO(MPO &Uodd, double dt, double x1, double x2)
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
    mwArray Op1(Indices(this->d_bose+1,this->d_bose+1));
    mwArray Op2(Indices(this->d_bose+1,this->d_bose+1));
    mwArray L_operator(Indices(this->d_bose+1,this->d_bose+1));
    mwArray Lsq_operator(Indices(this->d_bose+1,this->d_bose+1));
    mwArray id_bosonic = identityMatrix(this->d_bose+1);
    
    L_operator.fillWithZero();
    Lsq_operator.fillWithZero();
    Op1.fillWithZero();
    Op2.fillWithZero();
    
    for(int i=0; i<=this->d_bose; ++i)
    {
      L_operator.setElement(-0.5*this->d_bose+i,0.0,Indices(i,i));
    }
    
    //Now get L_squared
    Lsq_operator = L_operator*L_operator;
    
    //Construct Op1 und Op2 (imaginary unit including its sign is absorbed in the operators)
    double norm_factor = sqrt(this->d_bose/2*(this->d_bose/2+1));
    complex_t tmp;
    for(int i=0; i<=this->d_bose; ++i)
    {
      for(int j=0; j<=this->d_bose; ++j)
      {
	//Lp
	if(i==j+1)
	{
	  Op1.setElement(0.0,sqrt(i*(this->d_bose-j))/norm_factor,Indices(i,j));
	}
      }      
    }
    
    //Get Lm aut of Lp by simply taking the hermitian conjugate (and therefore getting the right sign!)
    Op2=Op1;  
    Op2.transpose(true);
  
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
    this->id_bose_mpo.reshape(Indices(this->d_bose+1,1,this->d_bose+1,1));
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
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  mwArray U2,S2,V2;
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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
  exp_odd.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_odd.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_odd.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  U2.clear(); S2.clear(); V2.clear();
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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


void CQEDBackgroundField::constructUevenMPO(MPO &Ueven, double dt)
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
    mwArray Op1(Indices(this->d_bose+1,this->d_bose+1));
    mwArray Op2(Indices(this->d_bose+1,this->d_bose+1));
    mwArray L_operator(Indices(this->d_bose+1,this->d_bose+1));
    mwArray Lsq_operator(Indices(this->d_bose+1,this->d_bose+1));
    mwArray id_bosonic = identityMatrix(this->d_bose+1);
    
    L_operator.fillWithZero();
    Lsq_operator.fillWithZero();
    Op1.fillWithZero();
    Op2.fillWithZero();
    
    for(int i=0; i<=this->d_bose; ++i)
    {
      L_operator.setElement(-0.5*this->d_bose+i,0.0,Indices(i,i));
    }
    
    //Now get L_squared
    Lsq_operator = L_operator*L_operator;
    
    //Construct Op1 und Op2 (imaginary unit including its sign is absorbed in the operators)
    double norm_factor = sqrt(this->d_bose/2*(this->d_bose/2+1));
    complex_t tmp;
    for(int i=0; i<=this->d_bose; ++i)
    {
      for(int j=0; j<=this->d_bose; ++j)
      {
	//Lp
	if(i==j+1)
	{
	  Op1.setElement(0.0,sqrt(i*(this->d_bose-j))/norm_factor,Indices(i,j));
	}
      }      
    }
    
    //Get Lm aut of Lp by simply taking the hermitian conjugate (and therefore getting the right sign!)
    Op2=Op1;  
    Op2.transpose(true);
  
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
    this->id_bose_mpo.reshape(Indices(this->d_bose+1,1,this->d_bose+1,1));
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
  exp_even.reshape(Indices(this->d_fermi,this->d_bose+1,this->d_fermi,this->d_fermi,this->d_bose+1,this->d_fermi));
  exp_even.permute(Indices(1,4,2,5,3,6));
  
  //Now take the first two indices together and perform a SVD
  exp_even.reshape(Indices(this->d_fermi*this->d_fermi,(this->d_bose+1)*(this->d_bose+1)*this->d_fermi*this->d_fermi));
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
  V.reshape(Indices(dimV[0]*(this->d_bose+1)*(this->d_bose+1),this->d_fermi*this->d_fermi));
  
  //No proceed with another SVD
  mwArray U2,S2,V2;
  wrapper::svd(V,U2,S2,V2);
  //Redistribute the singular values
  S2=sqrt(S2);
  U2.multiplyRight(S2);
  V2.multiplyLeft(S2);
  //Reshape U2 and get the second matrix of the MPO
  U2.reshape(Indices(dimV[0],this->d_bose+1,this->d_bose+1,-1));
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

void CQEDBackgroundField::constructFlipduMPO(MPO &mpo, int n)
{
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  
  complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
  mwArray sigma_minus(Indices(2,2),data_sm);
  
  complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
  mwArray sigma_plus(Indices(2,2),data_sp);
 
  //Provide bosonic operators
  mwArray Op1(Indices(this->d_bose+1,this->d_bose+1));
  mwArray id_bosonic = identityMatrix(this->d_bose+1);
  
  Op1.fillWithZero();
  
  //Construct Op1 
  //Since I want to reuse my old implementation which has hopping terms like sp Op1 sm + sm op2 sp I absorb i in Op1 and -i in Op2
  double norm_factor = sqrt(this->d_bose/2*(this->d_bose/2+1));
  complex_t tmp;
  for(int i=0; i<=this->d_bose; ++i)
  {
    for(int j=0; j<=this->d_bose; ++j)
    {
      //Lp
      if(i==j+1)
      {
	Op1.setElement(0.0,sqrt(i*(this->d_bose-j))/norm_factor,Indices(i,j));
      }
    }
  }
  
  /*cout << "Op1 = " << Op1 << endl;
  cout << "sp = " << sigma_plus << endl;
  cout << "sm = " << sigma_minus << endl;*/
  
   //Array to store the bosonic operators
   mwArray Bosonic_operators(Indices(4,this->d_bose+1,this->d_bose+1));
   Bosonic_operators.fillWithZero();
   
   for(int i=0; i<=this->d_bose; ++i)
   {
     for(int j=0; j<=this->d_bose; ++j)
     {
       //Identity
       Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));
       //01
       Bosonic_operators.setElement(Op1.getElement(Indices(i,j)),Indices(3,i,j));
     }
   }
   
  Bosonic_operators.reshape(Indices(4,(this->d_bose+1)*(this->d_bose+1)));  
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(4,this->d_fermi,this->d_fermi));
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
    }
  }
  
  Fermionic_operators.reshape(Indices(4,this->d_fermi*this->d_fermi));
  ///////////////////////////////////////////////////////////////////////////////////
  //Bond dimension for the MPO Matrices
  int D_bond_left = 1;
  int D_bond_right = 1;
  
  //Different matrices appearing
  mwArray Fermi_identity(Indices(D_bond_left,D_bond_right,4));
  mwArray Fermi_sigma_plus(Indices(D_bond_left,D_bond_right,4));
  mwArray Fermi_sigma_minus(Indices(D_bond_left,D_bond_right,4));
  
  mwArray Bose_identity(Indices(D_bond_left,D_bond_right,4));
  mwArray Bose_Op(Indices(D_bond_left,D_bond_right,4));
  
  //Fill with zeros
  Fermi_identity.fillWithZero();
  Fermi_sigma_plus.fillWithZero();
  Fermi_sigma_minus.fillWithZero();
  
  Bose_identity.fillWithZero();
  Bose_Op.fillWithZero();
  
  //Set elements
  Fermi_identity.setElement(ONE_c,Indices(0,0,0));
  Fermi_sigma_plus.setElement(ONE_c,Indices(0,0,1));
  Fermi_sigma_minus.setElement(ONE_c,Indices(0,0,2));
  
  Bose_identity.setElement(ONE_c,Indices(0,0,0));
  Bose_Op.setElement(ONE_c,Indices(0,0,3));
  
  
  //reshape
  Fermi_identity.reshape(Indices(D_bond_left*D_bond_right,4));
  Fermi_sigma_plus.reshape(Indices(D_bond_left*D_bond_right,4));
  Fermi_sigma_minus.reshape(Indices(D_bond_left*D_bond_right,4));
  
  Bose_identity.reshape(Indices(D_bond_left*D_bond_right,4));
  Bose_Op.reshape(Indices(D_bond_left*D_bond_right,4));
  
  //Make sure MPO is empty before filling
  mpo.initLength(this->N_total);
  
  //Temporary array for MPO and stuff
  mwArray res;
  
  for(int k=0; k<this->N_total; k++)
  {
    if((k%2)==1)
    {
      //Bosonic site where only identity can site
      res=reshape(Bose_identity*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose+1,this->d_bose+1));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
    else
    {
      //Fermionic site
      if(k==n)
      {
	//Site at which the operator starts
	res=reshape(Fermi_sigma_plus*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	//Operator for the flux
	res=reshape(Bose_Op*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose+1,this->d_bose+1));
	mpo.setOp(k+1,new Operator(permute(res,Indices(3,1,4,2))),true);
	//Finally put sigma_minus
	res=reshape(Fermi_sigma_minus*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	mpo.setOp(k+2,new Operator(permute(res,Indices(3,1,4,2))),true);
	//Correct k for next iteration
	k += 2;
      }
      else
      {
	//Just set fermionic identity
	res=reshape(Fermi_identity*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
    }	
  }
  
  //mpo.exportForMatlab("flipdu.m");  
}//constructFlipduMPO

void CQEDBackgroundField::constructFlipudMPO(MPO &mpo, int n)
{
   //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  
  complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
  mwArray sigma_minus(Indices(2,2),data_sm);
  
  complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
  mwArray sigma_plus(Indices(2,2),data_sp);
 
  //Provide bosonic operators
  mwArray Op2(Indices(this->d_bose+1,this->d_bose+1));
  mwArray id_bosonic = identityMatrix(this->d_bose+1);
  
  Op2.fillWithZero();
  
  //Construct Op2
  //Since I want to reuse my old implementation which has hopping terms like sp Op1 sm + sm op2 sp I absorb i in Op1 and -i in Op2
  double norm_factor = sqrt(this->d_bose/2*(this->d_bose/2+1));
  complex_t tmp;
  for(int i=0; i<=this->d_bose; ++i)
  {
    for(int j=0; j<=this->d_bose; ++j)
    {
      //Lp
      if(i+1==j)
      {
	Op2.setElement(0.0,-sqrt(i*(this->d_bose-j))/norm_factor,Indices(i,j));
      }
    }
  }
  
  /*cout << "Op2 = " << Op2 << endl;
  cout << "sp = " << sigma_plus << endl;
  cout << "sm = " << sigma_minus << endl;*/
  
   //Array to store the bosonic operators
   mwArray Bosonic_operators(Indices(4,this->d_bose+1,this->d_bose+1));
   Bosonic_operators.fillWithZero();
   
   for(int i=0; i<=this->d_bose; ++i)
   {
     for(int j=0; j<=this->d_bose; ++j)
     {
       //Identity
       Bosonic_operators.setElement(id_bosonic.getElement(Indices(i,j)),Indices(0,i,j));
       //02
       Bosonic_operators.setElement(Op2.getElement(Indices(i,j)),Indices(3,i,j));
     }
   }
   
  Bosonic_operators.reshape(Indices(4,(this->d_bose+1)*(this->d_bose+1)));  
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(4,this->d_fermi,this->d_fermi));
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
    }
  }
  
  Fermionic_operators.reshape(Indices(4,this->d_fermi*this->d_fermi));
  ///////////////////////////////////////////////////////////////////////////////////
  //Bond dimension for the MPO Matrices
  int D_bond_left = 1;
  int D_bond_right = 1;
  
  //Different matrices appearing
  mwArray Fermi_identity(Indices(D_bond_left,D_bond_right,4));
  mwArray Fermi_sigma_plus(Indices(D_bond_left,D_bond_right,4));
  mwArray Fermi_sigma_minus(Indices(D_bond_left,D_bond_right,4));
  
  mwArray Bose_identity(Indices(D_bond_left,D_bond_right,4));
  mwArray Bose_Op(Indices(D_bond_left,D_bond_right,4));
  
  //Fill with zeros
  Fermi_identity.fillWithZero();
  Fermi_sigma_plus.fillWithZero();
  Fermi_sigma_minus.fillWithZero();
  
  Bose_identity.fillWithZero();
  Bose_Op.fillWithZero();
  
  //Set elements
  Fermi_identity.setElement(ONE_c,Indices(0,0,0));
  Fermi_sigma_plus.setElement(ONE_c,Indices(0,0,1));
  Fermi_sigma_minus.setElement(ONE_c,Indices(0,0,2));
  
  Bose_identity.setElement(ONE_c,Indices(0,0,0));
  Bose_Op.setElement(ONE_c,Indices(0,0,3));
  
  
  //reshape
  Fermi_identity.reshape(Indices(D_bond_left*D_bond_right,4));
  Fermi_sigma_plus.reshape(Indices(D_bond_left*D_bond_right,4));
  Fermi_sigma_minus.reshape(Indices(D_bond_left*D_bond_right,4));
  
  Bose_identity.reshape(Indices(D_bond_left*D_bond_right,4));
  Bose_Op.reshape(Indices(D_bond_left*D_bond_right,4));
  
  //Make sure MPO is empty before filling
  mpo.initLength(this->N_total);
  
  //Temporary array for MPO and stuff
  mwArray res;
  
  for(int k=0; k<this->N_total; k++)
  {
    if((k%2)==1)
    {
      //Bosonic site where only identity can site
      res=reshape(Bose_identity*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose+1,this->d_bose+1));
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
    }
    else
    {
      //Fermionic site
      if(k==n)
      {
	//Site at which the operator starts
	res=reshape(Fermi_sigma_minus*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	//Operator for the flux
	res=reshape(Bose_Op*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose+1,this->d_bose+1));
	mpo.setOp(k+1,new Operator(permute(res,Indices(3,1,4,2))),true);
	//Finally put sigma_minus
	res=reshape(Fermi_sigma_plus*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	mpo.setOp(k+2,new Operator(permute(res,Indices(3,1,4,2))),true);
	//Correct k for next iteration
	k += 2;
      }
      else
      {
	//Just set fermionic identity
	res=reshape(Fermi_identity*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
    }	
  }
  
  //mpo.exportForMatlab("flipud.m");  
}//constructFlipudMPO

void CQEDBackgroundField::makeString(MPS &mps, int pos_start, int length)
{
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  
  complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
  mwArray sigma_minus(Indices(2,2),data_sm);
  
  complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
  mwArray sigma_plus(Indices(2,2),data_sp);
 
  //Provide bosonic operators
  mwArray Op1(Indices(this->d_bose+1,this->d_bose+1));
  mwArray id_bosonic = identityMatrix(this->d_bose+1);
  mwArray Op2;
  
  Op1.fillWithZero();
  
  //Construct Op1 
  //Since I want to reuse my old implementation which has hopping terms like sp Op1 sm + sm op2 sp I absorb i in Op1 and -i in Op2
  double norm_factor = sqrt(this->d_bose/2*(this->d_bose/2+1));
  complex_t tmp;
  for(int i=0; i<=this->d_bose; ++i)
  {
    for(int j=0; j<=this->d_bose; ++j)
    {
      //Lp
      if(i==j+1)
      {
	Op1.setElement(0.0,sqrt(i*(this->d_bose-j))/norm_factor,Indices(i,j));
      }
    }
  }
  
  //Get Op2 as the hermitian conjugate of Op1
  Op2=Op1;
  Op2.transpose(true);
  
  //In case pos_start is a negative number just prepare a string approximately in the middle
  if(pos_start<0)
  {
    //Find the middle of the chain where the string should reside
    int middle = (int) (this->N_fermions/2.0);
    //Translate into absolute site number
    pos_start = 2*(middle-(length-1)/2)-2;    
    cout << "Construction string starting at site " << pos_start << endl;
  }
  
  //Check if string fits on chain in principle (not paying attention to Gauss Law and stuff which might be an additionl limitation)
  if((mps.getLength()-pos_start)< 2*length)
  {
	  cerr << "MPS is not long enough to generate the desires string, program will be aborted" << endl;
	  exit(1);
  }
  
  //Check if startposition is a valid spin site and if it is a spin up or a spin down site
  if((pos_start%2)==1)
  {
	  cerr << "Cannot generate a string starting at a bosonic site, program will be aborted..." << endl;
	  exit(1);
  }
  int sign = pow(-1,pos_start/2+1);
  cout << "sign=" << sign << endl;
  cout << "length=" << length << endl;
  
  if(sign==1)
  {
	 mps.applyLocalOperator(pos_start,sigma_plus,true); 
	 for(int i=pos_start+1;i<pos_start+2*length;i+=2)
		mps.applyLocalOperator(i,Op1,true);
	 mps.applyLocalOperator(pos_start+2*length,sigma_minus,true);
  }
  else
  {
	 mps.applyLocalOperator(pos_start,sigma_minus,true); 
	 for(int i=pos_start+1;i<pos_start+2*length;i+=2)
		mps.applyLocalOperator(i,Op2,true);
	 mps.applyLocalOperator(pos_start+2*length,sigma_plus,true);	
  }
}//makeString

