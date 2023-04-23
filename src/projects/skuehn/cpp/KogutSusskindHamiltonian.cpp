#include "KogutSusskindHamiltonian.h"
#include "Indices.h"
#include "Contractor.h"

using namespace std;
using namespace shrt;

KogutSusskindHamiltonian::KogutSusskindHamiltonian(int N_, int d_,int D_,double mu_,double x_,double lambda_, Model model_,double noise_strength_):hamil(2*N_-1)
{
  N_fermions = N_;
  //As a bosonic link is between two fermions, the total number of sites in the chain is given by
  N_total = 2*N_-1;
  D = D_;
  mu = mu_;
  x = x_;
  d_bose = d_;
  d_fermi = 2;
  lambda = lambda_;
  model = model_;
  noise_strength = noise_strength_;
  terms_saved = false;
  lxprefactor = 1.0;
  
  
  cout << "Constructing Kogut-Susskind-Hamiltonian with " << endl
       << "N_fermions = " << N_fermions << endl
       << "N_total    = " << N_total << endl
       << "D          = " << D << endl
       << "mu         = " << mu << endl
       << "x          = " << x << endl
       << "d          = " << d_bose << endl
       << "la         = " << lambda << endl;
  
  //Construct the basic operators
  initOperators();
  //Construct the actual Hamiltonian
  if((this->model != cQEDeffective)&&(this->model != Zncgl)&&(this->model != Zncglnoise)&&(this->model != realZn))
    initHamiltonian();
  else if((this->model==Zncgl)||(this->model==Zncglnoise)||(this->model==realZn))
    initZncglHamiltonian();
  else
    initeffectiveHamiltonian();
}

//Destructor
KogutSusskindHamiltonian::~KogutSusskindHamiltonian(){
  //Clear Hamiltonian MPO
  hamil.clear();  
  //Clear all operators
  id_fermi.clear();
  sigma_x.clear();
  sigma_y.clear();
  sigma_z.clear();
  sigma_plus.clear();
  sigma_minus.clear();
  idpsz.clear();
  id_bose.clear();
  Lp.clear();
  Lm.clear();
  Lz.clear();
  Lzsquare.clear();
  //Clear temporary saved tensor products for time evolution operators in case the exists
  if(terms_saved)
  {
    hopping1.clear();
    hopping2.clear();
    lterm.clear();
    idpszleft.clear();
    idpszright.clear();
    firstidpszleft.clear();
    lastidpszright.clear();
    id_bose_mpo.clear();
    id_fermi_mpo.clear();
    if(this->model==cQEDnoise)
    {
      Lpnoise.clear();
      Lmnoise.clear();
    }    
  }  
}

//Update the Hamiltonian MPO, once it has been constructed constructed
void KogutSusskindHamiltonian::updateHMPO()
{
  if((this->model==cQED)||(this->model==cQEDnoise)||(this->model==Zn))
  {
    cout << "Updating cQED/Zn Hamiltonian" << endl;
    initHamiltonian();
  }
  else if((this->model==Zncgl)||(this->model==Zncglnoise))
  {
    cout << "Updating Zncgl Hamiltonian" << endl;
    initZncglHamiltonian();
  }
  else if(this->model==cQEDeffective)
  {
    cout << "Updating effective Hamiltonian" << endl;
    initeffectiveHamiltonian();
  }
  else
  {
    cout << "Error, unknown model, cannot update Hamiltonian" << endl;
    cout << "Program will be aborted..." << endl;
    exit(666);
  }
}

//Provide the basic operators according to the model
void KogutSusskindHamiltonian::initOperators()
{
  //First the fermionic operators (basically all kinds of Pauli Matrices)
  this->id_fermi=identityMatrix(this->d_fermi);
  
  complex_t data_x[] = {ZERO_c,ONE_c,ONE_c,ZERO_c};
  this->sigma_x=mwArray(Indices(2,2),data_x);
  
  complex_t data_y[] = {ZERO_c,I_c,-1.0*I_c,ZERO_c};
  this->sigma_y=mwArray(Indices(2,2),data_y);
  
  complex_t data_z[] = {ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  this->sigma_z=mwArray(Indices(2,2),data_z);
  
  complex_t data_sm[] = {ZERO_c,ONE_c,ZERO_c,ZERO_c};  
  this->sigma_minus=mwArray(Indices(2,2),data_sm);
  
  complex_t data_sp[] = {ZERO_c,ZERO_c,ONE_c,ZERO_c};  
  this->sigma_plus=mwArray(Indices(2,2),data_sp);
  
  complex_t data_idpsz[] = {2.0*ONE_c,ZERO_c,ZERO_c,ZERO_c};  
  this->idpsz=mwArray(Indices(2,2),data_idpsz); 
  
  //Generate empty matrices for the links
  this->Lp=mwArray(Indices(this->d_bose,this->d_bose));
  this->Lm=mwArray(Indices(this->d_bose,this->d_bose));
  this->Lz=mwArray(Indices(this->d_bose,this->d_bose));
  this->Lzsquare=mwArray(Indices(this->d_bose,this->d_bose));
  this->id_bose=identityMatrix(this->d_bose);

  Lz.fillWithZero();
  Lzsquare.fillWithZero();
  Lp.fillWithZero();
  Lm.fillWithZero();
  
  //Identities matrices already reshaped as MPO
  this->id_bose_mpo=this->id_bose;
  this->id_fermi_mpo=this->id_fermi;
  this->id_bose_mpo.reshape(Indices(this->d_bose,1,this->d_bose,1));
  this->id_fermi_mpo.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
  
  //Provide the operators for the links, according to the model
  if((this->model==cQED) || (this->model==cQEDnoise))
  {  
    cout << "Preparing links for cQED model" << endl;
    //This model works only for a fixed, even number of bosons on the links, which means an odd bosonic dimension
    if(((this->d_bose)%2)==0)
    {
      cout << "Error, cQED model works only for odd bosonic Hilbertspace dimension, program will be aborted..." << endl;
      exit(666);
    }
    for(int i=0; i<this->d_bose; ++i)
    {
	    Lz.setElement(-0.5*(this->d_bose-1)+i,0.0,Indices(i,i));
    }

    //Construct Lp and Lm
    //Since I want to reuse my old implementation which has hopping terms like sp Lp sm + sm Lm sp I absorb i in Op1 and -i in Op2
    double norm_factor = sqrt((this->d_bose-1)/2*((this->d_bose-1)/2+1));
    complex_t tmp;
    for(int i=0; i<this->d_bose; ++i)
    {
      for(int j=0; j<this->d_bose; ++j)
      {
	//Lp
	if(i==j+1)
	  this->Lp.setElement(0.0,sqrt(i*(this->d_bose-j-1))/norm_factor,Indices(i,j));
      }
    }

    //Get Op2 as the hermitian conjugate of Op1
    this->Lm=this->Lp;
    this->Lm.transpose(true);
    //Now get get the square of the L-operator
    Lzsquare = Lz*Lz;

  }
  else if(this->model==cQEDeffective)
  {
    cout << "Preparing links for effective cQED model" << endl;
    //This model works only for a fixed, even number of bosons on the links, which means an odd bosonic dimension
    if(((this->d_bose)%2)==0)
    {
      cout << "Error, cQED model works only for odd bosonic Hilbertspace dimension, program will be aborted..." << endl;
      exit(666);
    }
    for(int i=0; i<this->d_bose; ++i)
    {
	    Lz.setElement(-0.5*(this->d_bose-1)+i,0.0,Indices(i,i));
    }

    //Construct Lp and Lm
    //Here I need Lp and Lm without the imaginary unit absorbed in them, because I want to build Lx, which is (Lp + Lm)/2
    double norm_factor = sqrt((this->d_bose-1)/2*((this->d_bose-1)/2+1));
    complex_t tmp;
    for(int i=0; i<this->d_bose; ++i)
    {
      for(int j=0; j<this->d_bose; ++j)
      {
	//Lp
	if(i==j+1)
	  this->Lp.setElement(sqrt(i*(this->d_bose-j-1))/norm_factor,0.0,Indices(i,j));
      }
    }

    //Get Op2 as the hermitian conjugate of Op1
    this->Lm=this->Lp;
    this->Lm.transpose(true);
    
    //Now get Lx
    //double lxprefactor=1.0;
    //mwArray Lx = lxprefactor*0.5/(1.0+2*this->lambda)*(Lp+Lm);
    
    //Now, since I already have Lx, I can set Lp and Lm to the bosonic identity
    //this->Lm = this->id_bose;
    //this->Lp = this->id_bose;
    
    //Now get get the square of the L-operator plus the other stuff I need
    this->Lzsquare = this->Lz*this->Lz;
    
    //cout << "Lx=" << Lx << endl;
    
  }
  else if ((this->model==Zn) || (this->model==Zncgl) || (this->model==Zncglnoise) || (this->model==realZn))
  {
    cout << "Preparing links for Zn model" << endl;
    double edge_val = ((double)this->d_bose-1.0)/2.0;
    for(int i=0; i<this->d_bose; i++)
      for(int j=0; j<this->d_bose; j++)
	if(i==j)
	  this->Lz.setElement((-1.0*edge_val+i)*ONE_c,Indices(i,j));
	else if(((i+1)%this->d_bose)==j)
	  this->Lm.setElement(ONE_c,Indices(i,j));

    this->Lp=this->Lm;
    this->Lp.transpose(true);
    this->Lzsquare = this->Lz*this->Lz;
    
    //Case that I want to simulate a real Zn theory, where I have \sum_n(1- cos(2*Pi/d_bose*L_z,n))
    if(this->model==realZn)
    {
      //As the code allows me to compute the Matrix exponentials, I use cos(A) = 1/2*(exp(-iA)+exp(iA))
      cout << "Preparing modular electric term for Z_d model" << endl;
      mwArray expPA;
      wrapper::expm(this->Lz,expPA,2.0*M_PIl/this->d_bose*I_c);
      mwArray expMA;
      wrapper::expm(this->Lz,expMA,-2.0*M_PIl/this->d_bose*I_c);
      this->Lzsquare=this->id_bose-0.5*(expPA+expMA);
    }
  }
  else
  {
	cout << "Given model is not supported, program will be aborted"<< endl;
	exit(666);
  }
  //In case I need to see the operators for debugging
  /*cout << "Lp=" << this->Lp << endl;
  cout << "Lm=" << this->Lm << endl;
  cout << "Lz=" << this->Lz << endl;
  cout << "Lz²=" << this->Lzsquare << endl;*/
}

void KogutSusskindHamiltonian::initHamiltonian()
{
   //Array to store the bosonic operators
   mwArray Bosonic_operators(Indices(8,this->d_bose,this->d_bose));
   Bosonic_operators.fillWithZero();
   
   for(int i=0; i<this->d_bose; ++i)
   {
     for(int j=0; j<this->d_bose; ++j)
     {
       //Identity
       Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));
       //Lp
       Bosonic_operators.setElement(this->Lp.getElement(Indices(i,j)),Indices(4,i,j));
       //Lm
       Bosonic_operators.setElement(this->Lm.getElement(Indices(i,j)),Indices(5,i,j));
       //Lz
       Bosonic_operators.setElement(this->Lz.getElement(Indices(i,j)),Indices(6,i,j));
       //Lsquared
       Bosonic_operators.setElement(this->Lzsquare.getElement(Indices(i,j)),Indices(7,i,j));
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
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma_plus
      Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
      //sigma_minus
      Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
      //sigma_z
      Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(3,i,j));
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
      //Lp
      Bmatrices.setElement(this->x,0,Indices(2,5,4));
      if(this->model==cQEDnoise)
      {
	//cout << "Adding noise term for L+ with strength "<< this->noise_strength << endl;
	double norm_factor = sqrt((this->d_bose-1)/2*((this->d_bose-1)/2+1));
	Bmatrices.setElement(this->noise_strength*(-1.0)*I_c*norm_factor,Indices(0,8,4));
      }
      //Lm
      Bmatrices.setElement(this->x,0,Indices(3,6,5));
      if(this->model==cQEDnoise)
      {
	//cout << "Adding noise term for L- with strength "<< this->noise_strength << endl;
	double norm_factor = sqrt((this->d_bose-1)/2*((this->d_bose-1)/2+1));
	Bmatrices.setElement(this->noise_strength*I_c*norm_factor,Indices(0,8,5));
      }
      //Lz
      Bmatrices.setElement(ONE_c,Indices(0,7,6));
      Bmatrices.setElement(-2.0*this->lambda*sign,0,Indices(0,8,6));
      Bmatrices.setElement(ONE_c,Indices(1,8,6));
      Bmatrices.setElement(ONE_c,Indices(4,8,6));     
      //Lzsquare
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
	
	//Identity
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
}// initHamiltonian

//Initalize the Hamiltonian for the effective cQED Model
void KogutSusskindHamiltonian::initeffectiveHamiltonian()
{
  //First provide the modified single body term for bosonic sites, the factor 0.5 arises from the fact that Lx=1/2*(L_+ + L_-)
  mwArray L_mod = this->Lzsquare + this->lxprefactor*0.5/(1.0+2.0*this->lambda)*(Lp+Lm);
  
  //cout << "Lx=" << lxprefactor*0.5/(1.0+2.0*this->lambda)*(Lp+Lm) << endl;
  
  cout << "Trying to run effective model " << endl;
  //cout << "L_mod = " << L_mod << endl;
  
  //Array to store the bosonic operators
   mwArray Bosonic_operators(Indices(3,this->d_bose,this->d_bose));
   Bosonic_operators.fillWithZero();
   
   for(int i=0; i<this->d_bose; ++i)
   {
     for(int j=0; j<this->d_bose; ++j)
     {
       //Identity
       Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));
       //Lz
       Bosonic_operators.setElement(this->Lz.getElement(Indices(i,j)),Indices(1,i,j));
       //L_mod
       Bosonic_operators.setElement(L_mod.getElement(Indices(i,j)),Indices(2,i,j));
     }
   }
   
  Bosonic_operators.reshape(Indices(3,this->d_bose*this->d_bose));  
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(4,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();  
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma_plus
      Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
      //sigma_minus
      Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
      //sigma_z
      Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(3,i,j));
    }
  }
  
  Fermionic_operators.reshape(Indices(4,this->d_fermi*this->d_fermi));  
  
  cout << "Finished setting up operators" << endl;
  
  //Actual construciton of the MPO
  mwArray FirstFermi,LastFermi,EvenFermi,OddFermi,EvenBose,OddBose,res;
  int D_bond = 6;
  
  //Get the sizes right and fill with zeros, and reshape 
  FirstFermi.resize(Indices(1,D_bond,4));	FirstFermi.fillWithZero();
  LastFermi.resize(Indices(D_bond,1,4));	LastFermi.fillWithZero();
  OddFermi.resize(Indices(D_bond,D_bond,4));	OddFermi.fillWithZero();
  EvenFermi.resize(Indices(D_bond,D_bond,4));	EvenFermi.fillWithZero();
  OddBose.resize(Indices(D_bond,D_bond,3));	OddBose.fillWithZero();
  EvenBose.resize(Indices(D_bond,D_bond,3));	EvenBose.fillWithZero();
  
  //Now fill the MPO matrices
  //Determine the sign of the last site (I start counting by 1)
  double sign = pow(-1.0,this->N_fermions);
  
  //Identity, first fermi site (which is always odd)
  FirstFermi.setElement(ONE_c,Indices(0,0,0));
  FirstFermi.setElement((this->lambda/2.0-this->mu/2.0)*ONE_c,Indices(0,5,0));
  //sigma_plus, first fermi site
  FirstFermi.setElement(this->x*ONE_c,Indices(0,1,1));
  //sigma_minus, first fermi site
  FirstFermi.setElement(this->x*ONE_c,Indices(0,2,2));
  //sigma_z, first fermi site
  FirstFermi.setElement(-1.0*this->lambda*ONE_c,Indices(0,4,3));
  FirstFermi.setElement(-1.0*(this->lambda/2.0+this->mu/2.0)*ONE_c,Indices(0,5,3));
  
  //Identity, last fermi site (which is always odd)
  //cout << "sign = " << sign << endl;
  LastFermi.setElement((this->lambda/2.0+sign*this->mu/2.0)*ONE_c,Indices(0,0,0));
  LastFermi.setElement(ONE_c,Indices(5,0,0));
  //sigma_plus, last fermi site
  LastFermi.setElement(ONE_c,Indices(2,0,1));
  //sigma_minus, last fermi site
  LastFermi.setElement(ONE_c,Indices(1,0,2));
  //sigma_z, last fermi site
  LastFermi.setElement(sign*(this->lambda/2.0+this->mu/2.0)*ONE_c,Indices(0,0,3));
  LastFermi.setElement(this->lambda*ONE_c,Indices(3,0,3));
  
  //Identity, odd fermi site
  OddFermi.setElement(ONE_c,Indices(0,0,0));
  OddFermi.setElement((this->lambda/2.0-this->mu/2.0)*ONE_c,Indices(0,5,0));
  OddFermi.setElement(-2.0*this->lambda*ONE_c,Indices(3,3,0));
  OddFermi.setElement(ONE_c,Indices(5,5,0));
  //sigma_plus, odd fermi site
  OddFermi.setElement(this->x*ONE_c,Indices(0,1,1));
  OddFermi.setElement(ONE_c,Indices(2,5,1));
  //sigma_minus, odd fermi site
  OddFermi.setElement(this->x*ONE_c,Indices(0,2,2));
  OddFermi.setElement(ONE_c,Indices(1,5,2));
  //sigma_z, odd fermi site
  OddFermi.setElement(-1.0*this->lambda*ONE_c,Indices(0,4,3));
  OddFermi.setElement(-1.0*(this->lambda/2.0+this->mu/2.0)*ONE_c,Indices(0,5,3));
  OddFermi.setElement(this->lambda*ONE_c,Indices(3,5,3));
  
  //Identity, even fermi site
  EvenFermi.setElement(ONE_c,Indices(0,0,0));
  EvenFermi.setElement((this->lambda/2.0+this->mu/2.0)*ONE_c,Indices(0,5,0));
  EvenFermi.setElement(-2.0*this->lambda*ONE_c,Indices(3,3,0));
  EvenFermi.setElement(ONE_c,Indices(5,5,0));
  //sigma_plus, even fermi site
  EvenFermi.setElement(this->x*ONE_c,Indices(0,1,1));
  EvenFermi.setElement(ONE_c,Indices(2,5,1));
  //sigma_minus, even fermi site
  EvenFermi.setElement(this->x*ONE_c,Indices(0,2,2));
  EvenFermi.setElement(ONE_c,Indices(1,5,2));
  //sigma_z, even fermi site
  EvenFermi.setElement(-1.0*this->lambda*ONE_c,Indices(0,4,3));
  EvenFermi.setElement((this->lambda/2.0+this->mu/2.0)*ONE_c,Indices(0,5,3));
  EvenFermi.setElement(this->lambda*ONE_c,Indices(3,5,3));
  
  //Identity, odd bose site
  OddBose.setElement(ONE_c,Indices(0,0,0));
  OddBose.setElement(ONE_c,Indices(1,1,0));
  OddBose.setElement(ONE_c,Indices(2,2,0));
  OddBose.setElement(ONE_c,Indices(5,5,0));
  //Lz, odd bose site
  OddBose.setElement(ONE_c,Indices(0,3,1));
  OddBose.setElement(2.0*this->lambda*ONE_c,Indices(0,5,1));
  OddBose.setElement(ONE_c,Indices(3,5,1));
  OddBose.setElement(ONE_c,Indices(4,5,1));
  //Modified L^2, odd bose site
  OddBose.setElement((1.0+2.0*this->lambda)*ONE_c,Indices(0,5,2));

  //Identity, even bose site
  EvenBose.setElement(ONE_c,Indices(0,0,0));
  EvenBose.setElement(ONE_c,Indices(1,1,0));
  EvenBose.setElement(ONE_c,Indices(2,2,0));
  EvenBose.setElement(ONE_c,Indices(5,5,0));
  //Lz, even bose site
  EvenBose.setElement(ONE_c,Indices(0,3,1));
  EvenBose.setElement(-2.0*this->lambda*ONE_c,Indices(0,5,1));
  EvenBose.setElement(ONE_c,Indices(3,5,1));
  EvenBose.setElement(ONE_c,Indices(4,5,1));
  //Modified L^2, even bose site
  EvenBose.setElement((1.0+2.0*this->lambda)*ONE_c,Indices(0,5,2));  
  
  //Now reshape
  FirstFermi.reshape(Indices(D_bond,4));
  LastFermi.reshape(Indices(D_bond,4));
  OddFermi.reshape(Indices(D_bond*D_bond,4));
  EvenFermi.reshape(Indices(D_bond*D_bond,4));
  OddBose.reshape(Indices(D_bond*D_bond,3));
  EvenBose.reshape(Indices(D_bond*D_bond,3));
  
  cout << "Finished to set up matrices" << endl;
  
  //Now contract and fill the MPO
  bool odd_bose_set=false,even_bose_set=false,odd_fermi_set=false,even_fermi_set=false; 
  int fermi_counter=0;
  Operator* Local_op;
  
  for(int k=0; k<this->N_total; k++)
  {
    if((k%2)==0)
    {
      //Fermionic site
      fermi_counter++;
      if(k==0)
      {
	//First fermionic site
	//cout << "First fermi site " << endl;
	res=reshape(FirstFermi*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true); 
      }
      else if(k==(this->N_total-1))
      {
	//Last fermionic site
	//cout << "Last fermi site " << endl;
	res=reshape(LastFermi*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true); 
      }
      else if(((fermi_counter%2)==0) && (!even_fermi_set))
      {
	//First even site
	//cout << "Fermi site " << fermi_counter << ", setting even fermi operator" << endl;
	res=reshape(EvenFermi*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true); 
	even_fermi_set=true;
      }
      else if(((fermi_counter%2)==1) && (!odd_fermi_set))
      {
	//First odd site
	//cout << "Fermi site " << fermi_counter <<", setting odd fermi operator" << endl;
	res=reshape(OddFermi*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true); 
	odd_fermi_set=true;
      }
      else
      {
	//Previous site was already set
	//cout << "Fermi site " << fermi_counter << ", taking operator form site " << k-4 << endl;
	hamil.setOp(k,&hamil.getOp(k-4),false);
      }
    }
    else
    {
      //Bosonic site
      //cout << "Sign = " << pow(-1.0,fermi_counter) << endl;
      if(((fermi_counter%2)==0) && (!even_bose_set))
      {
	//First even site
	//cout << "Bose site " << fermi_counter <<", setting even bose operator" << endl;
	res=reshape(EvenBose*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose,this->d_bose));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true); 
	even_bose_set=true;
      }
      else if(((fermi_counter%2)==1) && (!odd_bose_set))
      {
	//First odd site
	//cout << "Bose site " << fermi_counter <<", setting odd bose operator" << endl;
	res=reshape(OddBose*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose,this->d_bose));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true); 
	odd_bose_set=true;
      }
      else
      {
	//Previous site was already set
	//cout << "Bose site " << fermi_counter << ", taking operator form site " << k-4 << endl;
	hamil.setOp(k,&hamil.getOp(k-4),false);
      }
    }
  }
  //cout << hamil << endl;
  //hamil.exportForMatlab("cQEDeffective.m");
}//initeffectiveHamiltonian()

void KogutSusskindHamiltonian::initZncglHamiltonian()
{
    //Arrays for the additional operators I need
    mwArray Op,Om,Kp,Km;
    
    wrapper::expm(this->Lz,Op,2.0*M_PIl*I_c/((double) this->d_bose));
    wrapper::expm(this->Lz,Om,-2.0*M_PIl*I_c/((double) this->d_bose));
    wrapper::expm(this->sigma_z,Kp,M_PIl*I_c/((double) this->d_bose));
    wrapper::expm(this->sigma_z,Km,-1.0*M_PIl*I_c/((double) this->d_bose));
    
    /*cout << "Op " << Op << endl;
    cout << "Om " << Om << endl;
    cout << "Kp " << Kp << endl;
    cout << "Km " << Km << endl;*/    
  
    //Array to store the bosonic operators
    mwArray Bosonic_operators(Indices(6,this->d_bose,this->d_bose));
    Bosonic_operators.fillWithZero();
    
    for(int i=0; i<this->d_bose; ++i)
    {
      for(int j=0; j<this->d_bose; ++j)
      {
	//Identity
	Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));
	//Lp
	Bosonic_operators.setElement(this->Lp.getElement(Indices(i,j)),Indices(1,i,j));
	//Lm
	Bosonic_operators.setElement(this->Lm.getElement(Indices(i,j)),Indices(2,i,j));
	//Lsquared
	Bosonic_operators.setElement(this->Lzsquare.getElement(Indices(i,j)),Indices(3,i,j));
	//Op
	Bosonic_operators.setElement(Op.getElement(Indices(i,j)),Indices(4,i,j));
	//Om
	Bosonic_operators.setElement(Om.getElement(Indices(i,j)),Indices(5,i,j));
      }
    }
    
  Bosonic_operators.reshape(Indices(6,this->d_bose*this->d_bose));  

  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(6,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();    

  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma_plus
      Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
      //sigma_minus
      Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
      //id+sigma_z
      Fermionic_operators.setElement(this->idpsz.getElement(Indices(i,j)),Indices(3,i,j));
      //Kp-operator
      Fermionic_operators.setElement(Kp.getElement(Indices(i,j)),Indices(4,i,j));
      //Km-operator
      Fermionic_operators.setElement(Km.getElement(Indices(i,j)),Indices(5,i,j));
    }
  }

  Fermionic_operators.reshape(Indices(6,this->d_fermi*this->d_fermi));
  
  //Site dependend constant (sites are counted as spin sites 1,2,...,N_fermi)
  complex_t c_odd,c_even;
  c_odd = exp(-1.0*I_c*M_PIl/((double) this->d_bose));
  c_even = exp(1.0*I_c*M_PIl/((double) this->d_bose));
  
  //Construct the matrices
  mwArray FermiFirst,FermiLast,FermiOdd,FermiEven,Bose;
  
  //Bond dimension of my MPO
  int D_bond = 6;
  
  FermiFirst.resize(Indices(1,D_bond,6));	FermiFirst.fillWithZero();
  FermiLast.resize(Indices(D_bond,1,6));	FermiLast.fillWithZero();
  FermiOdd.resize(Indices(D_bond,D_bond,6));	FermiOdd.fillWithZero();
  FermiEven.resize(Indices(D_bond,D_bond,6));	FermiEven.fillWithZero();
  Bose.resize(Indices(D_bond,D_bond,6));	Bose.fillWithZero();
  
  //Identity first Fermi site
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(2.0*this->lambda*ONE_c,Indices(0,5,0));
  //sigma_plus firt Fermi site
  FermiFirst.setElement(ONE_c,Indices(0,1,1));
  //sigma_plus firt Fermi site
  FermiFirst.setElement(ONE_c,Indices(0,2,2));
  //id+sigma_z firt Fermi site
  FermiFirst.setElement(-1.0*this->mu/2.0*ONE_c,Indices(0,5,3));
  //Kp first Fermi site
  FermiFirst.setElement(-1.0*this->lambda*c_odd,Indices(0,3,4));
  //Km first Fermi site
  FermiFirst.setElement(-1.0*this->lambda*conjugate(c_odd),Indices(0,4,5));
  
  
  //Determine the constant for the last site
  complex_t c_last;
  if((this->N_fermions%2)==0)
    c_last=c_even;
  else
    c_last=c_odd;
 
  //Identity last Fermi site
  FermiLast.setElement(2.0*this->lambda*ONE_c,Indices(0,0,0));
  FermiLast.setElement(ONE_c,Indices(5,0,0));
  //sigma_plus last Fermi site
  FermiLast.setElement(ONE_c,Indices(2,0,1));
  //sigma_minus last Fermi site
  FermiLast.setElement(ONE_c,Indices(1,0,2));
  //id+sigma_z last Fermi site
  FermiLast.setElement(pow(-1.0,this->N_fermions)*this->mu/2.0*ONE_c,Indices(0,0,3));
  //Kp last Fermi site
  FermiLast.setElement(-1.0*this->lambda*c_last,Indices(3,0,4));
  //Km last Fermi site
  FermiLast.setElement(-1.0*this->lambda*conjugate(c_last),Indices(4,0,5));
  
  //Identity odd Fermi site
  FermiOdd.setElement(ONE_c,Indices(0,0,0));
  FermiOdd.setElement(2.0*this->lambda*ONE_c,Indices(0,5,0));
  FermiOdd.setElement(ONE_c,Indices(5,5,0));
  //sigma_plus odd Fermi site
  FermiOdd.setElement(ONE_c,Indices(0,1,1));
  FermiOdd.setElement(ONE_c,Indices(2,5,1));
  //sigma_minus odd Fermi site
  FermiOdd.setElement(ONE_c,Indices(0,2,2));
  FermiOdd.setElement(ONE_c,Indices(1,5,2));
  //id+sigma_z odd Fermi site
  FermiOdd.setElement(-1.0*this->mu/2.0*ONE_c,Indices(0,5,3));
  //Kp odd Fermi site
  FermiOdd.setElement(-1.0*this->lambda*c_odd,Indices(3,3,4));
  //Km odd Fermi site
  FermiOdd.setElement(-1.0*this->lambda*conjugate(c_odd),Indices(4,4,5));
  
  //Identity even Fermi site
  FermiEven.setElement(ONE_c,Indices(0,0,0));
  FermiEven.setElement(2.0*this->lambda*ONE_c,Indices(0,5,0));
  FermiEven.setElement(ONE_c,Indices(5,5,0));
  //sigma_plus even Fermi site
  FermiEven.setElement(ONE_c,Indices(0,1,1));
  FermiEven.setElement(ONE_c,Indices(2,5,1));
  //sigma_minus even Fermi site
  FermiEven.setElement(ONE_c,Indices(0,2,2));
  FermiEven.setElement(ONE_c,Indices(1,5,2));
  //id+sigma_z even Fermi site
  FermiEven.setElement(1.0*this->mu/2.0*ONE_c,Indices(0,5,3));
  //Kp even Fermi site
  FermiEven.setElement(-1.0*this->lambda*c_even,Indices(3,3,4));
  //Km even Fermi site
  FermiEven.setElement(-1.0*this->lambda*conjugate(c_even),Indices(4,4,5));
  
  //Identity Bose site
  Bose.setElement(ONE_c,Indices(0,0,0));
  Bose.setElement(ONE_c,Indices(5,5,0));
  //Lp Bose site
  Bose.setElement(this->x*ONE_c,Indices(1,1,1));
  if(this->model==Zncglnoise)
  {
    cout << "Adding noise term for L+ with strength "<< this->noise_strength << endl;
    Bose.setElement(this->noise_strength*ONE_c,Indices(0,5,1));
  }
  //Lm Bose site
  Bose.setElement(this->x*ONE_c,Indices(2,2,2));
  if(this->model==Zncglnoise)
  {
    cout << "Adding noise term for L- with strength "<< this->noise_strength << endl;
    Bose.setElement(this->noise_strength*ONE_c,Indices(0,5,2));
  }
  //Lz^2 Bose site
  Bose.setElement(ONE_c,Indices(0,5,3));
  //Op Bose site
  Bose.setElement(ONE_c,Indices(0,3,4));
  Bose.setElement(ONE_c,Indices(4,5,4));
  //Om Bose site
  Bose.setElement(ONE_c,Indices(0,4,5));
  Bose.setElement(ONE_c,Indices(3,5,5));
  
  //Now reshape the matrices built previously
  FermiFirst.reshape(Indices(D_bond,6));
  FermiLast.reshape(Indices(D_bond,6));
  FermiOdd.reshape(Indices(D_bond*D_bond,6));
  FermiEven.reshape(Indices(D_bond*D_bond,6));
  Bose.reshape(Indices(D_bond*D_bond,6));
  
  //Now contract and fill the MPO
  bool bose_set=false,even_bose_set=false,odd_fermi_set=false,even_fermi_set=false; 
  mwArray res;
  Operator* Local_op; 
  
  for(int k=0; k<this->N_total; k++)
  {
    if(k==0)
    {
      //First fermionic site
      //cout << "First fermi site" << endl;
      res=reshape(FermiFirst*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true); 
    }
    else if(k==(this->N_total-1))
    {
      //Last fermionic site
      //cout << "Last fermi site" << endl;
      res=reshape(FermiLast*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true); 
    }
    else if(((k%2)==0) && ((((k+2)/2)%2)!=0))
    {
      //Odd fermionic site
      //cout << "Odd fermi site" << endl;
      if(odd_fermi_set)
	hamil.setOp(k,&hamil.getOp(k-4),false);
      else
      {
	res=reshape(FermiOdd*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true); 
	odd_fermi_set=true;
      }
    }
    else if (((k%2)==0) && ((((k+2)/2)%2)==0))
    {
      //Even Fermionic site
      //cout << "Even fermi site" << endl;
      if(even_fermi_set)
	hamil.setOp(k,&hamil.getOp(k-4),false);
      else
      {
	res=reshape(FermiEven*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true); 
	even_fermi_set=true;
      }
    }
    else
    {
      //Bosonic site
      //cout << "Bose site" << endl;
      if(bose_set)
	hamil.setOp(k,&hamil.getOp(k-2),false);
      else
      {
	res=reshape(Bose*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose,this->d_bose));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true); 
	bose_set=true;
      }
    }
  }  
  //hamil.exportForMatlab("Zncgl.m"); 
}//initZncglHamiltonian()

void KogutSusskindHamiltonian::constructTotalSzMPO(MPO &mpo) const
{
  cout << "Constructing the MPO total Sz-Operator" << endl;
  
  //Make sure given MPO is empty
  mpo.initLength(this->N_total);
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,2,2));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi; ++j)
    {
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(1,i,j));
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
      Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));      
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
      mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);      
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
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	
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
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
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
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      
    }
  }
}

void KogutSusskindHamiltonian::constructChargeMPO(MPO &mpo) const
{
  cout << "Constructing the MPO representation of the Charge Operator" << endl;
  
  //Make sure given MPO is empty
  mpo.initLength(this->N_total);
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi; ++j)
    {
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(1,i,j));
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
      Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));      
    }
  }
  Bosonic_operators.reshape(Indices(1,this->d_bose*this->d_bose));
  
  //Bond dimension for the MPO Matrices
  int D_bond = 2;
  
  
  //Matrices
  mwArray Bose_site(Indices(D_bond,D_bond,1));		Bose_site.fillWithZero();
  mwArray Fermi_first(Indices(1,D_bond,2));		Fermi_first.fillWithZero();
  mwArray Fermi_last(Indices(D_bond,1,2));		Fermi_last.fillWithZero();
  mwArray Fermi_even(Indices(D_bond,D_bond,2));		Fermi_even.fillWithZero();
  mwArray Fermi_odd(Indices(D_bond,D_bond,2));		Fermi_odd.fillWithZero();
  
  //Set non-zero entries
  Bose_site.setElement(ONE_c,Indices(0,0,0));
  Bose_site.setElement(ONE_c,Indices(1,1,0));
  
  Fermi_first.setElement(ONE_c,Indices(0,0,0));
  Fermi_first.setElement(-0.5*ONE_c,Indices(0,1,0));
  Fermi_first.setElement(0.5*ONE_c,Indices(0,1,1));  
  
  Fermi_last.setElement(0.5*pow(-1.0,(double) this->N_fermions)*ONE_c,Indices(0,0,0));
  Fermi_last.setElement(ONE_c,Indices(1,0,0));
  Fermi_last.setElement(0.5*ONE_c,Indices(0,0,1));
  
  Fermi_even.setElement(ONE_c,Indices(0,0,0));
  Fermi_even.setElement(ONE_c,Indices(1,1,0));
  Fermi_even.setElement(0.5*ONE_c,Indices(0,1,0));
  Fermi_even.setElement(0.5*ONE_c,Indices(0,1,1));
  
  Fermi_odd.setElement(ONE_c,Indices(0,0,0));
  Fermi_odd.setElement(ONE_c,Indices(1,1,0));
  Fermi_odd.setElement(-0.5*ONE_c,Indices(0,1,0));
  Fermi_odd.setElement(0.5*ONE_c,Indices(0,1,1));
  
  //Reshape
  Bose_site.reshape(Indices(D_bond*D_bond,1));
  Fermi_first.reshape(Indices(1*D_bond,2));
  Fermi_last.reshape(Indices(D_bond*1,2));
  Fermi_even.reshape(Indices(D_bond*D_bond,2));
  Fermi_odd.reshape(Indices(D_bond*D_bond,2));
  
  //Flags, if operators were already saved
  bool bose_saved=false, fermi_odd_saved=false, fermi_even_saved=false;
  
  //Now fill the MPO
  mwArray res;
  int fermi_site_number;
  //cout << "Starting to fill the MPO" << endl;
  for(int k=0; k<N_total;++k)
  {
    if((k%2)==1)
    {
      if(bose_saved)
      {
	//cout << "Bose site "<<k<<", taking operator from site " << k-2 << endl;
	mpo.setOp(k,&mpo.getOp(k-2),false);
      }
      else
      {
	//cout << "Bose site "<<k<<", setting operator " << endl;
	res=reshape(Bose_site*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose,this->d_bose));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);      
	bose_saved=true;
      }
    }
    else
    {
      //Fermionic site(even sites)
      if(k==0)
      {
	//First site
	res=reshape(Fermi_first*Fermionic_operators,Indices(1,D_bond,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	
      }
      else if(k==(N_total-1))
      {
	//Last site
	res=reshape(Fermi_last*Fermionic_operators,Indices(D_bond,1,2,2));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else
      {
	fermi_site_number = (int) (0.5*k+1);
	//cout << "Fermi Site Number " << fermi_site_number << endl;
	//Determine if site is even or odd 
	if((fermi_site_number%2)==0)
	{
	  //Even site
	  if(fermi_even_saved)
	  {
	    //cout << "Even fermi site "<<k<<", taking operator from "<< k-4 << endl;
	    mpo.setOp(k,&mpo.getOp(k-4),false);
	  }
	  else
	  {
	    //cout << "Even fermi site "<<k<<", setting operator " << endl;
	    res=reshape(Fermi_even*Fermionic_operators,Indices(D_bond,D_bond,2,2));
	    mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	    fermi_even_saved=true;
	  }
	}
	else
	{
	  //Odd site
	  if(fermi_odd_saved)
	  {
	    //cout << "Odd fermi site "<<k<<", taking operator from "<< k-4 << endl;
	    mpo.setOp(k,&mpo.getOp(k-4),false);
	  }
	  else
	  {
	    //cout << "Odd fermi site "<<k<<", setting operator " << endl;
	    res=reshape(Fermi_odd*Fermionic_operators,Indices(D_bond,D_bond,2,2));
	    mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	    fermi_odd_saved=true;
	  }
	}
      }      
    }
  }
  //cout << "Charge MPO " << mpo << endl;
  //mpo.exportForMatlab("NewChargeMPO.m");
}

void KogutSusskindHamiltonian::constructLopMPO(MPO &mpo, int n) const
{
  //cout << "Constructing L("<<n<<")" << endl;
  
  //As n is a DMRG site and L(n) acts on bosonic sites only, check if n is a bosonic index
  if((n%2)!=1)
  {
    cerr << "Error in KogutSusskindHamiltonian::constructLopMPO(MPO &mpo, int n): Index " << n << " is not a valid bosonic site" << endl;
    exit(1);
  }
  
  //Make sure MPO is empty
  mpo.initLength(this->N_total);

  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(1,2,2));
  Fermionic_operators.fillWithZero();  
 
  for(int i=0; i<this->d_fermi; ++i)
    for(int j=0; j<this->d_fermi; ++j)
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      
  Fermionic_operators.reshape(Indices(1,this->d_fermi*this->d_fermi));
 
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(2,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<this->d_bose; ++i)
  {
    for(int j=0; j<this->d_bose; ++j)
    {
      Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));      
      Bosonic_operators.setElement(this->Lz.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Bosonic_operators.reshape(Indices(2,(this->d_bose*this->d_bose)));
  
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
	res=reshape(Bose_Op*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true); 
      }
      else
      {
	//Just put the bosonic identity
	if(!bose_saved)
	{
	  res=reshape(Bose_identity*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
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

void KogutSusskindHamiltonian::constructLopsquareMPO(MPO &mpo, int n) const
{
  mpo.initLength(this->N_total);
  
  //As n is a DMRG site and L(n) acts on bosonic sites only, check if n is a bosonic index
  if((n%2)!=1 || n<1 || n>=this->N_total)
  {
    cerr << "Error in KogutSusskindHamiltonian::constructLopsquareMPO(MPO &mpo, int n): Index " << n << " is not a valid bosonic site" << endl;
    exit(1);
  }
  
  //Set the operator
  mpo.setOp(n,new Operator(reshape(this->Lzsquare,Indices(this->d_bose,1,this->d_bose,1))),true);
  
  //Fermionic identities
  mpo.setOp(0,new Operator(this->id_fermi_mpo),true);
  for(int i=2; i<this->N_total; i+=2)
    mpo.setOp(i,&mpo.getOp(0),false);
  
  //Remaining bosonic identities
  bool op_saved=false; 
  int pos_saved;
  for(int i=1; i<this->N_total; i++)
  {
    if(op_saved && mpo.isEmpty(i))
      mpo.setOp(i,&mpo.getOp(pos_saved),false);
    else if(mpo.isEmpty(i))
    {
      mpo.setOp(i,new Operator(this->id_bose_mpo),true);
      op_saved=true;
      pos_saved=i;
    }
  }
}

void KogutSusskindHamiltonian::constructLtotMPO(MPO &mpo) const
{
  //cout << "Constructing L("<<n<<")" << endl;
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(1,2,2));
  Fermionic_operators.fillWithZero();  
 
  for(int i=0; i<this->d_fermi; ++i)
    for(int j=0; j<this->d_fermi; ++j)
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      
  Fermionic_operators.reshape(Indices(1,this->d_fermi*this->d_fermi));
 
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(2,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<this->d_bose; ++i)
  {
    for(int j=0; j<this->d_bose; ++j)
    {
      Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));      
      Bosonic_operators.setElement(this->Lz.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Bosonic_operators.reshape(Indices(2,(this->d_bose)*(this->d_bose)));
  
  //Bond dimension for the MPO Matrices
  int D_bond_left = 2;
  int D_bond_right = 2;
  
  //Different matrices appearing
  mwArray Fermi_site(Indices(D_bond_left,D_bond_right,1));
  mwArray Fermi_first(Indices(1,D_bond_right,1));
  mwArray Fermi_last(Indices(D_bond_left,1,1));
  mwArray Bose_site(Indices(D_bond_left,D_bond_right,2));
  
  //Fill with zeros
  Fermi_site.fillWithZero();
  Fermi_first.fillWithZero();
  Fermi_last.fillWithZero();
  Bose_site.fillWithZero();
  
  //Set elements
  Fermi_site.setElement(ONE_c,Indices(0,0,0));
  Fermi_site.setElement(ONE_c,Indices(1,1,0));
  Fermi_first.setElement(ONE_c,Indices(0,0,0));
  Fermi_last.setElement(ONE_c,Indices(1,0,0));
  
  Bose_site.setElement(ONE_c,Indices(0,0,0));
  Bose_site.setElement(ONE_c,Indices(1,1,0));
  Bose_site.setElement(ONE_c,Indices(0,1,1));
  
  //reshape
  Fermi_site.reshape(Indices(D_bond_left*D_bond_right,1));
  Fermi_first.reshape(Indices(D_bond_right,1));
  Fermi_last.reshape(Indices(D_bond_left,1));
  
  Bose_site.reshape(Indices(D_bond_left*D_bond_right,2));
  
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
      //Bosonic site
      if(!bose_saved)
      {
	res=reshape(Bose_site*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
	bose_saved = true;
      }
      else
	mpo.setOp(k,&mpo.getOp(k-2),false);
    }
    else
    {
      if(k==0)
      {
	//First fermionic site
	res=reshape(Fermi_first*Fermionic_operators,Indices(1,D_bond_right,this->d_fermi,this->d_fermi));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else if(k==(this->N_total-1))
      {
	//Last fermionic site
	res=reshape(Fermi_last*Fermionic_operators,Indices(D_bond_left,1,this->d_fermi,this->d_fermi));
	mpo.setOp(k,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else
      {
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
  }
}//ConstructLtotMPO


void KogutSusskindHamiltonian::constructGaussMPO(MPO &mpo, int n) const
{
  //cout << "Constructing GaussMPO for site " << n << endl;
  
  //As n is a DMRG site and sigma_z acts on fermionic sites only, check if n is a bosonic index
  if((n%2)==1)
  {
    cerr << "Error in KogutSusskindHamiltonian::constructGaussMPO(MPO &mpo, int n): Index " << n << " is not a valid fermionic site" << endl;
    exit(1);
  }
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
  //Stagered sum of identiy and sigma_z
  mwArray stagg_idpsz;
  int sign=pow(-1,n/2+1);  
  stagg_idpsz = 0.5*(this->sigma_z+sign*this->id_fermi);
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi; ++j)
    {
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      Fermionic_operators.setElement(stagg_idpsz.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(2,this->d_fermi*this->d_fermi));
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(1,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<this->d_bose; ++i)
    for(int j=0; j<this->d_bose; ++j)
      Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));      
      
  Bosonic_operators.reshape(Indices(1,this->d_bose*this->d_bose));
  
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
	res=reshape(Bose_identity*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
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

void KogutSusskindHamiltonian::constructsigmazMPO(MPO &mpo, int n) const
{
  //cout << "Constructing sigma_z(" << n << ")" << endl;
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(1,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<this->d_bose; ++i)
    for(int j=0; j<this->d_bose; ++j)
      Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));
  
  Bosonic_operators.reshape(Indices(1,this->d_bose*this->d_bose));  
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();  
  
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma_z
      Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(1,i,j));
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
		  res=reshape(Bose_identity*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
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

void KogutSusskindHamiltonian::constructCondensateMPO(MPO &mpo) const
{
  //Make sure MPO is empty
  mpo.initLength(this->N_total);
  
  mwArray Op_fermi = sqrt(this->x)/this->N_fermions*0.5*this->idpsz;
  
  mwArray Bosonic_operators(Indices(1,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();
  
  mwArray Fermionic_operators(Indices(2,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();
  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //Operator
      Fermionic_operators.setElement(Op_fermi.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }
 
  for(int i=0; i<this->d_bose; ++i)
  {
    for(int j=0; j<this->d_bose;++j)
    {
      //Identity
      Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));
    }
  }
  
  //reshape
  Fermionic_operators.reshape(Indices(2,this->d_fermi*this->d_fermi));
  Bosonic_operators.reshape(Indices(1,this->d_bose*this->d_bose));
  
  //Bond dimension for the MPO Matrices
  int D_bond_left = 2;
  int D_bond_right = 2;
  
  //Different matrices appearing
  mwArray Fermi_Op_odd(Indices(D_bond_left,D_bond_right,2));
  mwArray Fermi_Op_even(Indices(D_bond_left,D_bond_right,2));
  mwArray Fermi_first(Indices(1,D_bond_right,2));
  mwArray Fermi_last(Indices(D_bond_left,1,2));
  
  mwArray Bose_identity(Indices(D_bond_left,D_bond_right,1));
  
  //Fill with zeros
  Fermi_Op_odd.fillWithZero();
  Fermi_Op_even.fillWithZero();
  Fermi_first.fillWithZero();
  Fermi_last.fillWithZero();
  
  Bose_identity.fillWithZero();
  
  //Determine if last site is an even or odd fermionic site
  double factor=1.0;
  if((this->N_fermions%2)!=0)
    factor = -1.0;
  
  //Set elements
  Fermi_Op_odd.setElement(ONE_c,Indices(0,0,0));
  Fermi_Op_odd.setElement(ONE_c,Indices(1,1,0));
  Fermi_Op_odd.setElement(-1.0*ONE_c,Indices(0,1,1));
  Fermi_Op_even.setElement(ONE_c,Indices(0,0,0));
  Fermi_Op_even.setElement(ONE_c,Indices(1,1,0));
  Fermi_Op_even.setElement(ONE_c,Indices(0,1,1));
  Fermi_first.setElement(ONE_c,Indices(0,0,0));
  Fermi_first.setElement(-1.0*ONE_c,Indices(0,1,1));
  Fermi_last.setElement(ONE_c,Indices(1,0,0));
  Fermi_last.setElement(factor*ONE_c,Indices(0,0,1));
  
  Bose_identity.setElement(ONE_c,Indices(0,0,0));
  Bose_identity.setElement(ONE_c,Indices(1,1,0));
  
  //Reshape matrices
  Fermi_Op_odd.reshape(Indices(D_bond_left*D_bond_right,2));
  Fermi_Op_even.reshape(Indices(D_bond_left*D_bond_right,2));
  Fermi_first.reshape(Indices(D_bond_right,2));
  Fermi_last.reshape(Indices(D_bond_left,2));
  Bose_identity.reshape(Indices(D_bond_left*D_bond_right,1));
  
  //Temporary storage
  mwArray res;
  
  //Sign which tells if actual is an even or odd fermionic site
  int sign=-1;
  
  //Flags if the operator was saved already
  bool fermi_odd_saved=false;
  bool fermi_even_saved=false;
  bool bose_saved=false;
  
  for(int i=0; i<this->N_total; i++)
  {
    
    if((i%2)==0)
    {
      //cout << "Sign of site "<< i << " is "  << sign <<endl;
      if(i==0)
      {
	//First fermionic site
	res=reshape(Fermi_first*Fermionic_operators,Indices(1,D_bond_right,2,2));
	mpo.setOp(i,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else if(i==(this->N_total-1))
      {
	//Last fermionic site
	res=reshape(Fermi_last*Fermionic_operators,Indices(D_bond_left,1,2,2));
	mpo.setOp(i,new Operator(permute(res,Indices(3,1,4,2))),true);
      }
      else
      {
	//Fermionic sites in between
	if(fermi_even_saved && fermi_odd_saved)
	  mpo.setOp(i,&mpo.getOp(i-4),false);
	
	else if(sign==-1)
	{
	    res=reshape(Fermi_Op_odd*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	    mpo.setOp(i,new Operator(permute(res,Indices(3,1,4,2))),true);
	    fermi_odd_saved = true;
	}
	else if(sign == 1)
	{
	  res=reshape(Fermi_Op_even*Fermionic_operators,Indices(D_bond_left,D_bond_right,2,2));
	  mpo.setOp(i,new Operator(permute(res,Indices(3,1,4,2))),true);
	  fermi_even_saved = true;
	}
      }      
      sign *= -1;      
    }
    else
    {
      if(!bose_saved)
      {
	res=reshape(Bose_identity*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
	mpo.setOp(i,new Operator(permute(res,Indices(3,1,4,2))),true);
	bose_saved = true;
      }
      else
	mpo.setOp(i,&mpo.getOp(i-2),false);
    }
  }
//mpo.exportForMatlab("CondensateMPO.m");
}//CondensateMPO

//Construct an explicit MPO for the penalty term to be able to see its contribution 
void KogutSusskindHamiltonian::constructPenaltyMPO(MPO &mpo, double strength) const
{
  if(strength<0.0)
    strength = this->lambda;
    
  cout << "Constructing Penalty MPO with strength " << strength << endl;  
  
  //Make sure MPO is empty
  mpo.initLength(N_total);
  
   //Array to store the bosonic operators
   mwArray Bosonic_operators(Indices(4,this->d_bose,this->d_bose));
   Bosonic_operators.fillWithZero();
   
   for(int i=0; i<this->d_bose; ++i)
   {
     for(int j=0; j<this->d_bose; ++j)
     {
       //Identity
       Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));
       //L
       Bosonic_operators.setElement(this->Lz.getElement(Indices(i,j)),Indices(2,i,j));
       //Lsquared
       Bosonic_operators.setElement(this->Lzsquare.getElement(Indices(i,j)),Indices(3,i,j));
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
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma_z
      Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(1,i,j));
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
  //mpo.exportForMatlab("PenaltyMPO.m");
}//Penalty

void KogutSusskindHamiltonian::constructCglPenaltyMPO(MPO &mpo, double strength) const
{
  if(strength<0.0)
    strength = this->lambda;
  
  cout << "Constructing cyclic Penalty MPO with strength " << strength << endl; 
  
  //Provide the additional operators needed
  //TODO These operators are duplicated from initZncglHamiltonian(), so I should maybe refactor this part
  mwArray Op,Om,Kp,Km;
  
  wrapper::expm(this->Lz,Op,2.0*M_PIl*I_c/((double) this->d_bose));
  wrapper::expm(this->Lz,Om,-2.0*M_PIl*I_c/((double) this->d_bose));
  wrapper::expm(this->sigma_z,Kp,M_PIl*I_c/((double) this->d_bose));
  wrapper::expm(this->sigma_z,Km,-1.0*M_PIl*I_c/((double) this->d_bose));
  
  //Make sure MPO is empty
  mpo.initLength(this->N_total);
  
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(3,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();   
  for(int i=0; i<this->d_bose; ++i)
  {
    for(int j=0; j<this->d_bose; ++j)
    {
      //Identity
      Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));
      //Op
      Bosonic_operators.setElement(Op.getElement(Indices(i,j)),Indices(1,i,j));
      //Om
      Bosonic_operators.setElement(Om.getElement(Indices(i,j)),Indices(2,i,j));
    }
  }   
  Bosonic_operators.reshape(Indices(3,this->d_bose*this->d_bose));  
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(3,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();  
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //Kp
      Fermionic_operators.setElement(Kp.getElement(Indices(i,j)),Indices(1,i,j));
      //Km
      Fermionic_operators.setElement(Km.getElement(Indices(i,j)),Indices(2,i,j));
    }
  } 
  Fermionic_operators.reshape(Indices(3,this->d_fermi*this->d_fermi));
  
  /***************************************************************************
   *                   Actual construction of the MPO
   * ************************************************************************/
  complex_t c_odd,c_even;
  c_odd = exp(-1.0*I_c*M_PIl/((double) this->d_bose));
  c_even = exp(1.0*I_c*M_PIl/((double) this->d_bose));
  
  //Construct the matrices
  mwArray FermiFirst,FermiLast,FermiOdd,FermiEven,Bose;
  
  //Bond dimension of my MPO
  int D_bond = 4;
  
  FermiFirst.resize(Indices(1,D_bond,3));	FermiFirst.fillWithZero();
  FermiLast.resize(Indices(D_bond,1,3));	FermiLast.fillWithZero();
  FermiOdd.resize(Indices(D_bond,D_bond,3));	FermiOdd.fillWithZero();
  FermiEven.resize(Indices(D_bond,D_bond,3));	FermiEven.fillWithZero();
  Bose.resize(Indices(D_bond,D_bond,3));	Bose.fillWithZero();
  
  //Identity first Fermi site
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(2.0*strength*ONE_c,Indices(0,3,0));
  //Kp first Fermi site
  FermiFirst.setElement(-1.0*strength*c_odd,Indices(0,1,1));
  //Km first Fermi site
  FermiFirst.setElement(-1.0*strength*conjugate(c_odd),Indices(0,2,2));  
  
  //Determine the constant for the last site
  complex_t c_last;
  if((this->N_fermions%2)==0)
    c_last=c_even;
  else
    c_last=c_odd;
  
  //Identity last Fermi site
  FermiLast.setElement(2.0*strength*ONE_c,Indices(0,0,0));
  FermiLast.setElement(ONE_c,Indices(3,0,0));
  //Kp last Fermi site
  FermiLast.setElement(-1.0*strength*c_last,Indices(1,0,1));
  //Km last Fermi site
  FermiLast.setElement(-1.0*strength*conjugate(c_last),Indices(2,0,2));
  
  //Identity odd Fermi site
  FermiOdd.setElement(ONE_c,Indices(0,0,0));
  FermiOdd.setElement(2.0*strength*ONE_c,Indices(0,3,0));
  FermiOdd.setElement(ONE_c,Indices(3,3,0));
  //Kp odd Fermi site
  FermiOdd.setElement(-1.0*strength*c_odd*ONE_c,Indices(1,1,1));
  //Km odd Fermi site
  FermiOdd.setElement(-1.0*strength*conjugate(c_odd)*ONE_c,Indices(2,2,2));
  
  //Identity even Fermi site
  FermiEven.setElement(ONE_c,Indices(0,0,0));
  FermiEven.setElement(2.0*strength*ONE_c,Indices(0,3,0));
  FermiEven.setElement(ONE_c,Indices(3,3,0));
  //Kp even Fermi site
  FermiEven.setElement(-1.0*strength*c_even*ONE_c,Indices(1,1,1));
  //Km even Fermi site
  FermiEven.setElement(-1.0*strength*conjugate(c_even)*ONE_c,Indices(2,2,2));
  
  //Identity Bose site
  Bose.setElement(ONE_c,Indices(0,0,0));
  Bose.setElement(ONE_c,Indices(3,3,0));
  //Op Bose site
  Bose.setElement(ONE_c,Indices(0,1,1));
  Bose.setElement(ONE_c,Indices(2,3,1));
  //Om Bose site
  Bose.setElement(ONE_c,Indices(0,2,2));
  Bose.setElement(ONE_c,Indices(1,3,2));
  
  //Now reshape the matrices built previously
  FermiFirst.reshape(Indices(D_bond,3));
  FermiLast.reshape(Indices(D_bond,3));
  FermiOdd.reshape(Indices(D_bond*D_bond,3));
  FermiEven.reshape(Indices(D_bond*D_bond,3));
  Bose.reshape(Indices(D_bond*D_bond,3));
  
  //Now contract and fill the MPO
  bool bose_set=false,even_bose_set=false,odd_fermi_set=false,even_fermi_set=false; 
  mwArray res;
  Operator* Local_op; 
  
  for(int k=0; k<this->N_total; k++)
  {
    if(k==0)
    {
      //First fermionic site
      //cout << "First fermi site" << endl;
      res=reshape(FermiFirst*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(k,Local_op,true); 
    }
    else if(k==(this->N_total-1))
    {
      //Last fermionic site
      //cout << "Last fermi site" << endl;
      res=reshape(FermiLast*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(k,Local_op,true); 
    }
    else if(((k%2)==0) && ((((k+2)/2)%2)!=0))
    {
      //Odd fermionic site
      //cout << "Odd fermi site" << endl;
      if(odd_fermi_set)
	mpo.setOp(k,&mpo.getOp(k-4),false);
      else
      {
	res=reshape(FermiOdd*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	mpo.setOp(k,Local_op,true); 
	odd_fermi_set=true;
      }
    }
    else if (((k%2)==0) && ((((k+2)/2)%2)==0))
    {
      //Even Fermionic site
      //cout << "Even fermi site" << endl;
      if(even_fermi_set)
	mpo.setOp(k,&mpo.getOp(k-4),false);
      else
      {
	res=reshape(FermiEven*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	mpo.setOp(k,Local_op,true); 
	even_fermi_set=true;
      }
    }
    else
    {
      //Bosonic site
      //cout << "Bose site" << endl;
      if(bose_set)
	mpo.setOp(k,&mpo.getOp(k-2),false);
      else
      {
	res=reshape(Bose*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose,this->d_bose));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	mpo.setOp(k,Local_op,true); 
	bose_set=true;
      }
    }
  }
  //mpo.exportForMatlab("cglPenalty.m");
}//constructCglPenaltyMPO(MPO &mpo, double strength)

//Basically deprecated, because I have operators for Uodd and Ueven seperately, but I keep it as a wrapper
void KogutSusskindHamiltonian::constructEvolutionMPO(MPO &Uodd, MPO &Ueven, double dt)
{    
  this->constructUoddMPO(Uodd,dt);
  this->constructUevenMPO(Ueven,dt);
}//constructEvolutionMPO

MPS* KogutSusskindHamiltonian::constructInitalMPS(InitState state, int extent)
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



void KogutSusskindHamiltonian::constructUoddMPO(MPO &Uodd, double dt)
{
  if(!this->terms_saved)
  {
    //Case in which the single terms are not yet built and have to be generated (also I construct only Uodd, I provide the terms for Ueven, since it is very likely that I will use this one afterwards)
    
    //Since I want to split id+s3 on every site in 0.5*(id+s3) + 0.5*(id+s3) and take one term to the even MPO and one to the odd MPO I get a problem at the first site and the last site, because one term is missing. To avoid this I have to calculate an extra MPO for the first site and an extra MPO for the site before the last site. Graphically the problem looks like this
    //x--o--x--o--x--o--x
    //|  |  |  |  |  |  |
    //-------  |  -------
    //|  |  -------  |  |
    //Note: That this scheme works, the chain has to have an even number of spin sites, otherwise it could happen, that this hast to be done for an even side MPO
    
    cout << "Generating terms and storing them" << endl;
    
    //sigma_plus Op1 sigma_minus
    this->hopping1=kron(this->sigma_plus,kron(this->Lp,this->sigma_minus));
    //sigma_minus Op2 sigma_plus
    this->hopping2=kron(this->sigma_minus,kron(this->Lm,this->sigma_plus));
    //mwArray id L^2  id
    this->lterm=kron(this->id_fermi,kron(this->Lzsquare,this->id_fermi));
    // 0.5*idpsz id id
    this->idpszleft=kron(0.5*this->idpsz,kron(this->id_bose,this->id_fermi));
    // id id 0.5*idpsz
    this->idpszright=kron(this->id_fermi,kron(this->id_bose,0.5*this->idpsz));;
    // idpsz id id
    this->firstidpszleft=kron(-1.0*this->idpsz,kron(this->id_bose,this->id_fermi));;  
    // id id idpsz
    this->lastidpszright=kron(this->id_fermi,kron(this->id_bose,this->idpsz));  
    
    //Terms for the noise
    if(this->model==cQEDnoise)
    {
      double norm_factor = sqrt((this->d_bose-1)/2*((this->d_bose-1)/2+1));
      this->Lpnoise=kron(this->id_fermi,kron((-1.0)*I_c*norm_factor*this->Lp,this->id_fermi));
      this->Lmnoise=kron(this->id_fermi,kron(I_c*norm_factor*this->Lm,this->id_fermi));
    }
    else if(this->model==Zncglnoise)
    {
      cout << "Preparing noise terms" << endl;
      this->Lpnoise=kron(this->id_fermi,kron(this->Lp,this->id_fermi));
      this->Lmnoise=kron(this->id_fermi,kron(this->Lm,this->id_fermi));
    }
    
    /*cout << "Lp-noise: " << this->Lpnoise << endl;
    cout << "Lm-noise: " << this->Lmnoise << endl;*/
    
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
  if((this->model==cQEDnoise)||(this->model==Zncglnoise))
    arg_odd = arg_odd + this->noise_strength*Lpnoise + this->noise_strength*Lmnoise;
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
  if((this->model==cQEDnoise)||(this->model==Zncglnoise))
    arg_odd = arg_odd + this->noise_strength*Lpnoise + this->noise_strength*Lmnoise;
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
  if((this->model==cQEDnoise)||(this->model==Zncglnoise))
    arg_odd = arg_odd + this->noise_strength*Lpnoise + this->noise_strength*Lmnoise;
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
}//constructUoddMPO

void KogutSusskindHamiltonian::constructUoddMPO(MPO &Uodd, double dt, double x1, double x2)
{
  //WARNING: This routine is deprecated and not compliant with noisy evolutions, therefore I should not be used anymore
  cout << "WARNING: This routine is deprecated and not compliant with noisy evolutions, use with care!" << endl;
  if(!this->terms_saved)
  {
    //Case in which the single terms are not yet built and have to be generated (also I construct only Uodd, I provide the terms for Ueven, since it is very likely that I will use this one afterwards)
     
    //Since I want to split id+s3 on every site in 0.5*(id+s3) + 0.5*(id+s3) and take one term to the even MPO and one to the odd MPO I get a problem at the first site and the last site, because one term is missing. To avoid this I have to calculate an extra MPO for the first site and an extra MPO for the site before the last site. Graphically the problem looks like this
    //x--o--x--o--x--o--x
    //|  |  |  |  |  |  |
    //-------  |  -------
    //|  |  -------  |  |
    //Note: That this scheme works, the chain has to have an even number of spin sites, otherwise it could happen, that this hast to be done for an even side MPO
    
    cout << "Generating terms and storing them" << endl;
    
    //sigma_plus Op1 sigma_minus
    this->hopping1=kron(this->sigma_plus,kron(this->Lp,this->sigma_minus));
    //sigma_minus Op2 sigma_plus
    this->hopping2=kron(this->sigma_minus,kron(this->Lm,this->sigma_plus));
    //mwArray id L^2  id
    this->lterm=kron(this->id_fermi,kron(this->Lzsquare,this->id_fermi));
    // 0.5*idpsz id id
    this->idpszleft=kron(0.5*this->idpsz,kron(this->id_bose,this->id_fermi));
    // id id 0.5*idpsz
    this->idpszright=kron(this->id_fermi,kron(this->id_bose,0.5*this->idpsz));;
    // idpsz id id
    this->firstidpszleft=kron(-1.0*this->idpsz,kron(this->id_bose,this->id_fermi));;  
    // id id idpsz
    this->lastidpszright=kron(this->id_fermi,kron(this->id_bose,this->idpsz));
    
    //Terms for the noise
    if(this->model==cQEDnoise)
    {
      double norm_factor = sqrt((this->d_bose-1)/2*((this->d_bose-1)/2+1));
      this->Lpnoise=kron(this->id_fermi,kron((-1.0)*I_c*norm_factor*this->Lp,this->id_fermi));
      this->Lmnoise=kron(this->id_fermi,kron(I_c*norm_factor*this->Lm,this->id_fermi));
    }
    else if(this->model==Zncglnoise)
    {
      this->Lpnoise=kron(this->id_fermi,kron(this->Lp,this->id_fermi));
      this->Lmnoise=kron(this->id_fermi,kron(this->Lm,this->id_fermi));
    }
    
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
}//constructUoddMPO (version with x1 and x2 given)


void KogutSusskindHamiltonian::constructUevenMPO(MPO &Ueven, double dt)
{
  if(!this->terms_saved)
  {
    
    //Since I want to split id+s3 on every site in 0.5*(id+s3) + 0.5*(id+s3) and take one term to the even MPO and one to the odd MPO I get a problem at the first site and the last site, because one term is missing. To avoid this I have to calculate an extra MPO for the first site and an extra MPO for the site before the last site. Graphically the problem looks like this
    //x--o--x--o--x--o--x
    //|  |  |  |  |  |  |
    //-------  |  -------
    //|  |  -------  |  |
    //Note: That this scheme works, the chain has to have an even number of spin sites, otherwise it could happen, that this hast to be done for an even side MPO
    
    cout << "Generating terms and storing them" << endl;
    
    //sigma_plus Op1 sigma_minus
    this->hopping1=kron(this->sigma_plus,kron(this->Lp,this->sigma_minus));
    //sigma_minus Op2 sigma_plus
    this->hopping2=kron(this->sigma_minus,kron(this->Lm,this->sigma_plus));
    //mwArray id L^2  id
    this->lterm=kron(this->id_fermi,kron(this->Lzsquare,this->id_fermi));
    // 0.5*idpsz id id
    this->idpszleft=kron(0.5*this->idpsz,kron(this->id_bose,this->id_fermi));
    // id id 0.5*idpsz
    this->idpszright=kron(this->id_fermi,kron(this->id_bose,0.5*this->idpsz));
    // idpsz id id
    this->firstidpszleft=kron(-1.0*this->idpsz,kron(this->id_bose,this->id_fermi));  
    // id id idpsz
    this->lastidpszright=kron(this->id_fermi,kron(this->id_bose,this->idpsz)); 
    
    //Terms for the noise
    if(this->model==cQEDnoise)
    {
      double norm_factor = sqrt((this->d_bose-1)/2*((this->d_bose-1)/2+1));
      this->Lpnoise=kron(this->id_fermi,kron((-1.0)*I_c*norm_factor*this->Lp,this->id_fermi));
      this->Lmnoise=kron(this->id_fermi,kron(I_c*norm_factor*this->Lm,this->id_fermi));
    }
    else if(this->model==Zncglnoise)
    {
      cout << "Preparing noise terms" << endl;
      this->Lpnoise=kron(this->id_fermi,kron(this->Lp,this->id_fermi));
      this->Lmnoise=kron(this->id_fermi,kron(this->Lm,this->id_fermi));
    }
    /*cout << "Lp-noise: " << this->Lpnoise << endl;
    cout << "Lm-noise: " << this->Lmnoise << endl;*/
    
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
  if((this->model==cQEDnoise)||(this->model==Zncglnoise))
    arg_even = arg_even + this->noise_strength*Lpnoise + this->noise_strength*Lmnoise;
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
}//constructUevenMPO

void KogutSusskindHamiltonian::constructFlipduMPO(MPO &mpo, int n)
{
   //Array to store the bosonic operators
   mwArray Bosonic_operators(Indices(4,this->d_bose,this->d_bose));
   Bosonic_operators.fillWithZero();
   
   for(int i=0; i<this->d_bose; ++i)
   {
     for(int j=0; j<this->d_bose; ++j)
     {
       //Identity
       Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));
       //Lp
       Bosonic_operators.setElement(this->Lp.getElement(Indices(i,j)),Indices(3,i,j));
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
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma_plus
      Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
      //sigma_minus
      Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
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
      res=reshape(Bose_identity*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
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
	res=reshape(Bose_Op*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
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

void KogutSusskindHamiltonian::constructFlipudMPO(MPO &mpo, int n)
{
   //Array to store the bosonic operators
   mwArray Bosonic_operators(Indices(4,this->d_bose,this->d_bose));
   Bosonic_operators.fillWithZero();
   
   for(int i=0; i<this->d_bose; ++i)
   {
     for(int j=0; j<this->d_bose; ++j)
     {
       //Identity
       Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));
       //Lm
       Bosonic_operators.setElement(this->Lm.getElement(Indices(i,j)),Indices(3,i,j));
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
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma_plus
      Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
      //sigma_minus
      Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
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
      res=reshape(Bose_identity*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
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
	res=reshape(Bose_Op*Bosonic_operators,Indices(D_bond_left,D_bond_right,this->d_bose,this->d_bose));
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

void KogutSusskindHamiltonian::makeString(MPS &mps, int pos_start, int length)
{  
  //In case pos_start is a negative number just prepare a string approximately in the middle
  if(pos_start<0)
  {
    //Find the middle of the chain where the string should reside
    int middle = (int) (this->N_fermions/2.0);
    //Translate into absolute site number
    pos_start = 2*(middle-(length-1)/2)-2;    
    cout << "Construction string starting at site " << pos_start << endl;
  }
  
  //Check if string fits on chain in principle (not paying attention to Gauss Law and stuff which might be an additional limitation)
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
  /*cout << "sign=" << sign << endl;
  cout << "length=" << length << endl;*/
  
  if(sign==1)
  {
	 mps.applyLocalOperator(pos_start,this->sigma_plus,true); 
	 for(int i=pos_start+1;i<pos_start+2*length;i+=2)
		mps.applyLocalOperator(i,Lp,true);
	 mps.applyLocalOperator(pos_start+2*length,this->sigma_minus,true);
  }
  else
  {
	 mps.applyLocalOperator(pos_start,this->sigma_minus,true); 
	 for(int i=pos_start+1;i<pos_start+2*length;i+=2)
		mps.applyLocalOperator(i,Lm,true);
	 mps.applyLocalOperator(pos_start+2*length,this->sigma_plus,true);	
  }
}//makeString

void KogutSusskindHamiltonian::EvolutionMPO(MPO &mpo,double dt)
{
  /************************************************************************************************
   * Provide the operators
  ************************************************************************************************/
  //Provide Pauli-Matrices
  mwArray id_fermionic = identityMatrix(2);
  cout << "WARNING: changed modified identity, correction is missing" << endl;
  //mwArray id_fermionic_modified = 1.0/((double)this->N_total)*id_fermionic;
  mwArray id_fermionic_modified = id_fermionic;
  
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
  mwArray Op1(Indices(this->d_bose+1,this->d_bose+1));
  mwArray Op2(Indices(this->d_bose+1,this->d_bose+1));
  mwArray Theta(Indices(this->d_bose+1,this->d_bose+1));
  mwArray L_operator(Indices(this->d_bose+1,this->d_bose+1));
  mwArray Lsq_operator(Indices(this->d_bose+1,this->d_bose+1));
  mwArray id_bosonic = identityMatrix(this->d_bose+1);
  mwArray  Lsqpmodid;
  
  L_operator.fillWithZero();
  Lsq_operator.fillWithZero();
  Op1.fillWithZero();
  Op2.fillWithZero();
  Theta.fillWithZero();
  
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
    }
  }
  
  //Get Op2 as the hermitian conjugate of Op1
  Op2=Op1;
  Op2.transpose(true);

  
  //Now get get the square of the L-operator
  Lsq_operator = L_operator*L_operator;
 
  //I need -i*dt all the time, that's why I save it WITH MINUS SIGN
  complex_t taylor_coefficient;
  taylor_coefficient = -1.0*dt*I_c;
  
  cout << "Korrekturfaktor: " << 1.0/((double)this->N_total) << endl;
  
  //Lsqpmodid = 1.0/((double)this->N_total)*id_bosonic + taylor_coefficient*Lsq_operator;
  cout<< "WARNING: changed L^2 operator, correction factor is missing" << endl;
  Lsqpmodid = id_bosonic + taylor_coefficient*Lsq_operator;
  
  /************************************************************************************************
   * Prepare the operator arrays
  ************************************************************************************************/
  //Array to store the bosonic operators
  mwArray Bosonic_operators(Indices(4,this->d_bose+1,this->d_bose+1));
  Bosonic_operators.fillWithZero();
  
  for(int i=0; i<(this->d_bose+1); ++i)
  {
    for(int j=0; j<(this->d_bose+1); ++j)
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
  
  Bosonic_operators.reshape(Indices(4,(this->d_bose+1)*(this->d_bose+1)));

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
	res=reshape(bose*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose+1,this->d_bose+1));
	mpo.setOp(i,new Operator(permute(res,Indices(3,1,4,2))),true);
	bose_set = true;
      }
    }
  }
}

//More optimized version, which avoids having the identy seperately, maybe this circumvents the bug I'm dealing with right now
void KogutSusskindHamiltonian::EvolutionMPO2(MPO &mpo,double dt)
{
  //I need -i*dt all the time, that's why I save it WITH MINUS SIGN
  complex_t taylor_coefficient;
  taylor_coefficient = -1.0*dt*I_c; 
  
  if((this->model==cQED) || (this->model==Zn) || (this->model==Zncgl))
  {    
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
	Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));
	//Lp
	Bosonic_operators.setElement(this->Lp.getElement(Indices(i,j)),Indices(1,i,j));
	//Lm
	Bosonic_operators.setElement(this->Lm.getElement(Indices(i,j)),Indices(2,i,j));
	//Lzsquare
	Bosonic_operators.setElement(this->Lzsquare.getElement(Indices(i,j)),Indices(3,i,j));
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
	Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
	//id+sigma_z
	Fermionic_operators.setElement(this->idpsz.getElement(Indices(i,j)),Indices(1,i,j));
	//sigma_plus
	Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(2,i,j));
	//sigma_minus
	Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(3,i,j));
      }
    }
    Fermionic_operators.reshape(Indices(4,this->d_fermi*this->d_fermi));
    
    /************************************************************************************************
    * Now prepare the matrices
    ************************************************************************************************/
    //Bond dimension of my MPO
    int D_bond = 4;
    double sign;
    
    //cout << "Filling matrices" << endl;
    //First fermionic site
    mwArray fermi_first(Indices(1,D_bond,4)); fermi_first.fillWithZero();
    sign=-1.0;
    //Identity
    fermi_first.setElement(ONE_c,Indices(0,0,0));    
    fermi_first.setElement(1.0/((double)this->N_total)*ONE_c,Indices(0,3,0));
    //I+sigma_z(n)
    fermi_first.setElement(this->mu/2.0*sign*ONE_c*taylor_coefficient,Indices(0,3,1));  
    //sigma_plus
    fermi_first.setElement(ONE_c,Indices(0,1,2));  
    //sigma_minus
    fermi_first.setElement(ONE_c,Indices(0,2,3));   
    
    //Even fermionic sites in between
    mwArray fermi_even(Indices(D_bond,D_bond,4)); fermi_even.fillWithZero();
    sign=1.0;
    //Identity
    fermi_even.setElement(ONE_c,Indices(0,0,0));
    fermi_even.setElement(ONE_c,Indices(3,3,0));  
    fermi_even.setElement(1.0/((double)this->N_total)*ONE_c,Indices(0,3,0)); 
    //I+sigma_z(n)
    fermi_even.setElement(this->mu/2.0*sign*ONE_c*taylor_coefficient,Indices(0,3,1));  
    //sigma_plus
    fermi_even.setElement(ONE_c,Indices(0,1,2));
    fermi_even.setElement(ONE_c,Indices(2,3,2));  
    //sigma_minus
    fermi_even.setElement(ONE_c,Indices(0,2,3)); 
    fermi_even.setElement(ONE_c,Indices(1,3,3)); 
    
    //Odd fermionic sites in between
    mwArray fermi_odd(Indices(D_bond,D_bond,4)); fermi_odd.fillWithZero();
    sign=-1.0;
    //Identity
    fermi_odd.setElement(ONE_c,Indices(0,0,0));
    fermi_odd.setElement(ONE_c,Indices(3,3,0));  
    fermi_odd.setElement(1.0/((double)this->N_total)*ONE_c,Indices(0,3,0));  
    //I+sigma_z(n)
    fermi_odd.setElement(this->mu/2.0*sign*ONE_c*taylor_coefficient,Indices(0,3,1));  
    //sigma_plus
    fermi_odd.setElement(ONE_c,Indices(0,1,2));
    fermi_odd.setElement(ONE_c,Indices(2,3,2));  
    //sigma_minus
    fermi_odd.setElement(ONE_c,Indices(0,2,3)); 
    fermi_odd.setElement(ONE_c,Indices(1,3,3));
    
    //Bosonic sites in between
    mwArray bose(Indices(D_bond,D_bond,4)); bose.fillWithZero();
    //Identity
    bose.setElement(ONE_c,Indices(0,0,0));
    bose.setElement(ONE_c,Indices(3,3,0));  
    bose.setElement(1.0/((double)this->N_total)*ONE_c,Indices(0,3,0));  
    //Op1
    bose.setElement(this->x*taylor_coefficient*ONE_c,Indices(1,1,1));  
    //Op2
    bose.setElement(this->x*taylor_coefficient*ONE_c,Indices(2,2,2));  
    //L^2
    bose.setElement(taylor_coefficient*ONE_c,Indices(0,3,3));  
    
    //Last fermionic site
    mwArray fermi_last(Indices(D_bond,1,4)); fermi_last.fillWithZero();
    sign = pow(-1.0,((double)this->N_fermions));
    //Identity
    fermi_last.setElement(ONE_c,Indices(3,0,0));
    fermi_last.setElement(1.0/((double)this->N_total)*ONE_c,Indices(0,0,0));
    //I+sigma_z(n)
    fermi_last.setElement(this->mu/2.0*sign*ONE_c*taylor_coefficient,Indices(0,0,1));
    //sigma_plus
    fermi_last.setElement(ONE_c,Indices(2,0,2));
    //sigma_minus
    fermi_last.setElement(ONE_c,Indices(1,0,3)); 
  
    //Reshape operators
    fermi_first.reshape(Indices(1*D_bond,4));
    fermi_even.reshape(Indices(D_bond*D_bond,4));
    fermi_odd.reshape(Indices(D_bond*D_bond,4));
    fermi_last.reshape(Indices(1*D_bond,4));
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
  else if(this->model==cQEDeffective)
  {
    cout << "Time evolution for effective cQED model is still experimental, use with care" << endl;
    //First provide the modified single body term for bosonic sites, the factor 0.5 arises from the fact that Lx=1/2*(L_+ + L_-)
    mwArray L_mod = this->Lzsquare + this->lxprefactor*0.5/(1.0+2.0*this->lambda)*(Lp+Lm);
    
    //cout << "Lx=" << lxprefactor*0.5/(1.0+2.0*this->lambda)*(Lp+Lm) << endl;
   
    //cout << "L_mod = " << L_mod << endl;
    
    //Array to store the bosonic operators
    mwArray Bosonic_operators(Indices(3,this->d_bose,this->d_bose));
    Bosonic_operators.fillWithZero();
    
    for(int i=0; i<this->d_bose; ++i)
    {
      for(int j=0; j<this->d_bose; ++j)
      {
	//Identity
	Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));
	//Lz
	Bosonic_operators.setElement(this->Lz.getElement(Indices(i,j)),Indices(1,i,j));
	//L_mod
	Bosonic_operators.setElement(L_mod.getElement(Indices(i,j)),Indices(2,i,j));
      }
    }
    
    Bosonic_operators.reshape(Indices(3,this->d_bose*this->d_bose));  
    
    //Array to store the fermionic operators
    mwArray Fermionic_operators(Indices(4,this->d_fermi,this->d_fermi));
    Fermionic_operators.fillWithZero();  
    
    for(int i=0; i<this->d_fermi; ++i)
    {
      for(int j=0; j<this->d_fermi;++j)
      {
	//Identity
	Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
	//sigma_plus
	Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
	//sigma_minus
	Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
	//sigma_z
	Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(3,i,j));
      }
    }
    
    Fermionic_operators.reshape(Indices(4,this->d_fermi*this->d_fermi));  
    
    cout << "Finished setting up operators" << endl;
    
    //Actual construciton of the MPO
    mwArray FirstFermi,LastFermi,EvenFermi,OddFermi,EvenBose,OddBose,res;
    int D_bond = 6;
    
    //Get the sizes right and fill with zeros, and reshape 
    FirstFermi.resize(Indices(1,D_bond,4));	FirstFermi.fillWithZero();
    LastFermi.resize(Indices(D_bond,1,4));	LastFermi.fillWithZero();
    OddFermi.resize(Indices(D_bond,D_bond,4));	OddFermi.fillWithZero();
    EvenFermi.resize(Indices(D_bond,D_bond,4));	EvenFermi.fillWithZero();
    OddBose.resize(Indices(D_bond,D_bond,3));	OddBose.fillWithZero();
    EvenBose.resize(Indices(D_bond,D_bond,3));	EvenBose.fillWithZero();
    
    //Now fill the MPO matrices
    //Determine the sign of the last site (I start counting by 1)
    double sign = pow(-1.0,this->N_fermions);
    
    //Identity, first fermi site (which is always odd)
    FirstFermi.setElement(ONE_c,Indices(0,0,0));
    FirstFermi.setElement(ONE_c/this->N_fermions+(this->lambda/2.0-this->mu/2.0)*taylor_coefficient,Indices(0,5,0));
    //sigma_plus, first fermi site
    FirstFermi.setElement(this->x*taylor_coefficient,Indices(0,1,1));
    //sigma_minus, first fermi site
    FirstFermi.setElement(this->x*taylor_coefficient,Indices(0,2,2));
    //sigma_z, first fermi site
    FirstFermi.setElement(-1.0*this->lambda*taylor_coefficient,Indices(0,4,3));
    FirstFermi.setElement(-1.0*(this->lambda/2.0+this->mu/2.0)*taylor_coefficient,Indices(0,5,3));
    
    //Identity, last fermi site (which is always odd)
    //cout << "sign = " << sign << endl;
    LastFermi.setElement(ONE_c/this->N_fermions+(this->lambda/2.0+sign*this->mu/2.0)*taylor_coefficient,Indices(0,0,0));
    LastFermi.setElement(ONE_c,Indices(5,0,0));
    //sigma_plus, last fermi site
    LastFermi.setElement(ONE_c,Indices(2,0,1));
    //sigma_minus, last fermi site
    LastFermi.setElement(ONE_c,Indices(1,0,2));
    //sigma_z, last fermi site
    LastFermi.setElement(sign*(this->lambda/2.0+this->mu/2.0)*taylor_coefficient,Indices(0,0,3));
    LastFermi.setElement(this->lambda*taylor_coefficient,Indices(3,0,3));
    
    //Identity, odd fermi site
    OddFermi.setElement(ONE_c,Indices(0,0,0));
    OddFermi.setElement(ONE_c/this->N_fermions+(this->lambda/2.0-this->mu/2.0)*taylor_coefficient,Indices(0,5,0));
    OddFermi.setElement(-2.0*this->lambda*taylor_coefficient,Indices(3,3,0));
    OddFermi.setElement(ONE_c,Indices(5,5,0));
    //sigma_plus, odd fermi site
    OddFermi.setElement(this->x*taylor_coefficient,Indices(0,1,1));
    OddFermi.setElement(ONE_c,Indices(2,5,1));
    //sigma_minus, odd fermi site
    OddFermi.setElement(this->x*taylor_coefficient,Indices(0,2,2));
    OddFermi.setElement(ONE_c,Indices(1,5,2));
    //sigma_z, odd fermi site
    OddFermi.setElement(-1.0*this->lambda*taylor_coefficient,Indices(0,4,3));
    OddFermi.setElement(-1.0*(this->lambda/2.0+this->mu/2.0)*taylor_coefficient,Indices(0,5,3));
    OddFermi.setElement(this->lambda*taylor_coefficient,Indices(3,5,3));
    
    //Identity, even fermi site
    EvenFermi.setElement(ONE_c,Indices(0,0,0));
    EvenFermi.setElement(ONE_c/this->N_fermions+(this->lambda/2.0+this->mu/2.0)*taylor_coefficient,Indices(0,5,0));
    EvenFermi.setElement(-2.0*this->lambda*taylor_coefficient,Indices(3,3,0));
    EvenFermi.setElement(ONE_c,Indices(5,5,0));
    //sigma_plus, even fermi site
    EvenFermi.setElement(this->x*taylor_coefficient,Indices(0,1,1));
    EvenFermi.setElement(ONE_c,Indices(2,5,1));
    //sigma_minus, even fermi site
    EvenFermi.setElement(this->x*taylor_coefficient,Indices(0,2,2));
    EvenFermi.setElement(ONE_c,Indices(1,5,2));
    //sigma_z, even fermi site
    EvenFermi.setElement(-1.0*this->lambda*taylor_coefficient,Indices(0,4,3));
    EvenFermi.setElement((this->lambda/2.0+this->mu/2.0)*taylor_coefficient,Indices(0,5,3));
    EvenFermi.setElement(this->lambda*taylor_coefficient,Indices(3,5,3));
    
    //Identity, odd bose site
    OddBose.setElement(ONE_c,Indices(0,0,0));
    OddBose.setElement(ONE_c,Indices(1,1,0));
    OddBose.setElement(ONE_c,Indices(2,2,0));
    OddBose.setElement(ONE_c,Indices(5,5,0));
    //Lz, odd bose site
    OddBose.setElement(ONE_c,Indices(0,3,1));
    OddBose.setElement(2.0*this->lambda*taylor_coefficient,Indices(0,5,1));
    OddBose.setElement(ONE_c,Indices(3,5,1));
    OddBose.setElement(ONE_c,Indices(4,5,1));
    //Modified L^2, odd bose site
    OddBose.setElement((1.0+2.0*this->lambda)*taylor_coefficient,Indices(0,5,2));

    //Identity, even bose site
    EvenBose.setElement(ONE_c,Indices(0,0,0));
    EvenBose.setElement(ONE_c,Indices(1,1,0));
    EvenBose.setElement(ONE_c,Indices(2,2,0));
    EvenBose.setElement(ONE_c,Indices(5,5,0));
    //Lz, even bose site
    EvenBose.setElement(ONE_c,Indices(0,3,1));
    EvenBose.setElement(-2.0*this->lambda*taylor_coefficient,Indices(0,5,1));
    EvenBose.setElement(ONE_c,Indices(3,5,1));
    EvenBose.setElement(ONE_c,Indices(4,5,1));
    //Modified L^2, even bose site
    EvenBose.setElement((1.0+2.0*this->lambda)*taylor_coefficient,Indices(0,5,2));  
    
    //Now reshape
    FirstFermi.reshape(Indices(D_bond,4));
    LastFermi.reshape(Indices(D_bond,4));
    OddFermi.reshape(Indices(D_bond*D_bond,4));
    EvenFermi.reshape(Indices(D_bond*D_bond,4));
    OddBose.reshape(Indices(D_bond*D_bond,3));
    EvenBose.reshape(Indices(D_bond*D_bond,3));
    
    cout << "Finished to set up matrices" << endl;
    
    //Now contract and fill the MPO
    bool odd_bose_set=false,even_bose_set=false,odd_fermi_set=false,even_fermi_set=false; 
    int fermi_counter=0;
    Operator* Local_op;
    
    for(int k=0; k<this->N_total; k++)
    {
      if((k%2)==0)
      {
	//Fermionic site
	fermi_counter++;
	if(k==0)
	{
	  //First fermionic site
	  //cout << "First fermi site " << endl;
	  res=reshape(FirstFermi*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
	  Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	  hamil.setOp(k,Local_op,true); 
	}
	else if(k==(this->N_total-1))
	{
	  //Last fermionic site
	  //cout << "Last fermi site " << endl;
	  res=reshape(LastFermi*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
	  Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	  hamil.setOp(k,Local_op,true); 
	}
	else if(((fermi_counter%2)==0) && (!even_fermi_set))
	{
	  //First even site
	  //cout << "Fermi site " << fermi_counter << ", setting even fermi operator" << endl;
	  res=reshape(EvenFermi*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	  Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	  hamil.setOp(k,Local_op,true); 
	  even_fermi_set=true;
	}
	else if(((fermi_counter%2)==1) && (!odd_fermi_set))
	{
	  //First odd site
	  //cout << "Fermi site " << fermi_counter <<", setting odd fermi operator" << endl;
	  res=reshape(OddFermi*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	  Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	  hamil.setOp(k,Local_op,true); 
	  odd_fermi_set=true;
	}
	else
	{
	  //Previous site was already set
	  //cout << "Fermi site " << fermi_counter << ", taking operator form site " << k-4 << endl;
	  hamil.setOp(k,&hamil.getOp(k-4),false);
	}
      }
      else
      {
	//Bosonic site
	//cout << "Sign = " << pow(-1.0,fermi_counter) << endl;
	if(((fermi_counter%2)==0) && (!even_bose_set))
	{
	  //First even site
	  //cout << "Bose site " << fermi_counter <<", setting even bose operator" << endl;
	  res=reshape(EvenBose*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose,this->d_bose));
	  Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	  hamil.setOp(k,Local_op,true); 
	  even_bose_set=true;
	}
	else if(((fermi_counter%2)==1) && (!odd_bose_set))
	{
	  //First odd site
	  //cout << "Bose site " << fermi_counter <<", setting odd bose operator" << endl;
	  res=reshape(OddBose*Bosonic_operators,Indices(D_bond,D_bond,this->d_bose,this->d_bose));
	  Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	  hamil.setOp(k,Local_op,true); 
	  odd_bose_set=true;
	}
	else
	{
	  //Previous site was already set
	  //cout << "Bose site " << fermi_counter << ", taking operator form site " << k-4 << endl;
	  hamil.setOp(k,&hamil.getOp(k-4),false);
	}
      }
    }    
  }
  else
    cout << "Given model not supported, did NOT create an MPO" << endl;
}

//Test of Gauss Projector, alternative which I came up with while going to Wuerzburg
void KogutSusskindHamiltonian::GaussProjector(MPO &mpo, int n)
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

//Projector, project all sites at once
void KogutSusskindHamiltonian::GaussProjector(MPO &mpo)
{
  //Initalize the MPO
  mpo.clear();
  mpo.initLength(this->N_total);
  //cout << "Constructing global projector" << endl;
  
  //Prepare matrices
  mwArray LeftSpinMatrix(Indices(this->d_fermi,1,this->d_fermi,this->d_bose));
  mwArray RightSpinMatrix(Indices(this->d_fermi,this->d_bose,this->d_fermi,1));
  mwArray SpinMatrixOdd(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_bose));
  mwArray SpinMatrixEven(Indices(this->d_fermi,this->d_bose,this->d_fermi,this->d_bose));
  mwArray FieldMatrix(Indices(this->d_bose,this->d_bose,this->d_bose,this->d_bose));
  
  LeftSpinMatrix.fillWithZero();
  RightSpinMatrix.fillWithZero();
  SpinMatrixOdd.fillWithZero();
  SpinMatrixEven.fillWithZero();
  FieldMatrix.fillWithZero();
  
  //Fill the field matrices (nothing but an identity)
  for(int i=0; i<this->d_bose; i++)
    FieldMatrix.setElement(ONE_c,Indices(i,i,i,i));
  
  //Generate the matrix for the spin sites, this is where the actual projection takes place
  int spin,field_value_left,field_value_right;  
  for(int i=0; i<this->d_fermi; i++)
  {
    for(int bond_index_left=0; bond_index_left<this->d_bose; bond_index_left++)
    {
      for(int bond_index_right=0;bond_index_right<this->d_bose; bond_index_right++)
      {
	//Determine spin value
	spin = (i==0) ? 1 : -1;
	field_value_left = bond_index_left - (int)((this->d_bose-1.0)/2.0);
	field_value_right = bond_index_right - (int)((this->d_bose-1.0)/2.0);
	
	//Odd site (-1)^n=-1 in the middle
	if((field_value_left+(int)((-1+spin)/2.0))==field_value_right)
	  SpinMatrixOdd.setElement(ONE_c,Indices(i,bond_index_left,i,bond_index_right));	
	
	//Even site (-1)^n=1 in the middle
	if((field_value_left+(int)((1+spin)/2.0))==field_value_right)
	  SpinMatrixEven.setElement(ONE_c,Indices(i,bond_index_left,i,bond_index_right));	
      }
    }
  }
  
  //Do the same for the first and last spin matrix, keep in mind that L(0)=0=L(N)
  //First spin matrix
  for(int i=0; i<this->d_fermi; i++)
  {
      for(int bond_index_right=0;bond_index_right<this->d_bose; bond_index_right++)
      {
	//Determine spin value
	spin = (i==0) ? 1 : -1;
	field_value_left = 0;
	field_value_right = bond_index_right - (int)((this->d_bose-1.0)/2.0);
	
	//First site (always odd)
	if((field_value_left+(int)((-1+spin)/2.0))==field_value_right)
	  LeftSpinMatrix.setElement(ONE_c,Indices(i,0,i,bond_index_right));	
      }
   }
  //Last spin matrix
  for(int i=0; i<this->d_fermi; i++)
  {
    for(int bond_index_left=0; bond_index_left<this->d_bose; bond_index_left++)
    {
	//Determine spin values and field values
	spin = (i==0) ? 1 : -1;
	field_value_left = bond_index_left - (int)((this->d_bose-1.0)/2.0);
	field_value_right = 0;
	
	//Last site (should in princial be even, because I work with a even number of spins to ensure total charge zero, but for the sake of generatlity I calculate it explicitly)
	if((field_value_left+(int)((pow(-1,this->N_fermions)+spin)/2.0))==field_value_right)
	  RightSpinMatrix.setElement(ONE_c,Indices(i,bond_index_left,i,0));	
      }
   }
  
  //Generate the MPO
  int spin_even_set = 0;
  int spin_odd_set = 0;
  int field_set = 0;
  int spin_counter = 0;
  
  for(int i=0; i< this->N_total; i++)
  {
    if((i%2)==0)
    {
      //Spin site
      spin_counter++;
      if(spin_counter==1)
      {
	//First spin site
	mpo.setOp(i,new Operator(LeftSpinMatrix),true);
      }
      else if(spin_counter==this->N_fermions)
      {
	//Last spin site
	mpo.setOp(i,new Operator(RightSpinMatrix),true);
      }
      else
      {
	if((spin_counter%2)==0)
	{
	  //Even spin site
	  if(spin_even_set!=0)
	  {
	    //I already set the matrix and can work with a pointer
	    mpo.setOp(i,&mpo.getOp(spin_even_set),false);
	  }
	  else
	  {
	    //I have to set the matrix
	    mpo.setOp(i,new Operator(SpinMatrixEven),true);
	    spin_even_set = i;
	  }
	}
	else
	{
	  //Odd spin site
	  if(spin_odd_set!=0)
	  {
	    //I already set the matrix and can work with a pointer
	    mpo.setOp(i,&mpo.getOp(spin_odd_set),false);
	  }
	  else
	  {
	    //I have to set the matrix
	    mpo.setOp(i,new Operator(SpinMatrixOdd),true);
	    spin_odd_set = i;
	  }
	}
      }      
    }
    else
    {
      //Bosonic site 
      if(field_set!=0)
      {
	//I already set the matrix and can work with a pointer
	mpo.setOp(i,&mpo.getOp(field_set),false);
      }
      else
      { 
	//I have to set the matrix
	mpo.setOp(i,new Operator(FieldMatrix),true);
	field_set = i;
      }
    }
  } 
}


void KogutSusskindHamiltonian::getStringOperator(MPO &mpo, int l, int n, bool adjoint) const
{
  mpo.initLength(this->N_total);

  if(n<=0)
    n=round(this->N_fermions/2.0 - l/2.0);
  
  //Some error checking
  if(l<1 || (n+l)>this->N_fermions || n>this->N_fermions)
  {
    cout << "Error in KogutSusskindHamiltonian::getStringOperator(), received l=" << l << " and n=" << n << " which does not allow to construct a valid string" << endl;
    exit(666);
  }
  if((l%2)==0)    
  {
    cout << "Error in KogutSusskindHamiltonian::getStringOperator(), received l=" << l << ", only odd length are supported" << endl;
    exit(666);
  }
  
  //Determine if starting site is odd or even
  bool is_odd = ((n%2)==1);
  //If I want the adjoint operator, I simply pretend the starting site is the opposite
  if(adjoint)
    is_odd = !is_odd;
  
  //Compute starting and end index
  int pos_start = 2*n-2;
  int pos_end = 2*(n+l)-2;
  
  
  //Set the operators for start and end
  if(is_odd)
  {
    mpo.setOp(pos_start,new Operator(reshape(this->sigma_minus,Indices(this->d_fermi,1,this->d_fermi,1))),true);
    mpo.setOp(pos_end,new Operator(reshape(this->sigma_plus,Indices(this->d_fermi,1,this->d_fermi,1))),true);
  }
  else
  {
    mpo.setOp(pos_start,new Operator(reshape(this->sigma_plus,Indices(this->d_fermi,1,this->d_fermi,1))),true);
    mpo.setOp(pos_end,new Operator(reshape(this->sigma_minus,Indices(this->d_fermi,1,this->d_fermi,1))),true);
  }
  
  //Set the L+ or L- operators
  bool op_saved=false;
  int pos_op_saved;
  
  for(int i=pos_start+1; i<pos_end; i+=2)
  {
    if(op_saved)
      mpo.setOp(i,&mpo.getOp(pos_op_saved),false);
    else if(is_odd)
    {
      mpo.setOp(i,new Operator(reshape(this->Lm,Indices(this->d_bose,1,this->d_bose,1))),true);
      pos_op_saved=i;
      op_saved=true;
    }
    else
    {
      mpo.setOp(i,new Operator(reshape(this->Lp,Indices(this->d_bose,1,this->d_bose,1))),true);
      pos_op_saved=i;
      op_saved=true;
    }
  } 
  
  //Now the identities
  bool id_bose_saved=false,id_fermi_saved=false;
  int pos_id_bose,pos_id_fermi;
  for(int i=0; i<this->N_total; i++)
  {
    //Set operators
    if((i%2)==0 && id_fermi_saved && mpo.isEmpty(i))
      mpo.setOp(i,&mpo.getOp(pos_id_fermi),false);
    else if((i%2)==0 && !id_fermi_saved && mpo.isEmpty(i))
    {
      mpo.setOp(i,new Operator(this->id_fermi_mpo),true);
      id_fermi_saved=true;
      pos_id_fermi=i;
    }
    else if((i%2)==1 && id_bose_saved && mpo.isEmpty(i))
      mpo.setOp(i,&mpo.getOp(pos_id_bose),false);
    else if((i%2)==1 && !id_bose_saved && mpo.isEmpty(i))
    {
      mpo.setOp(i,new Operator(this->id_bose_mpo),true);
      id_bose_saved=true;
      pos_id_bose=i;
    }
  } 
}

void KogutSusskindHamiltonian::getMesonOperator(MPO &mpo, int n, double k, double width) const
{
  mpo.initLength(this->N_total);
  
  //Some error checking
  if(n<1 || n>(this->N_fermions-1))
  {
    cout << "Error in KogutSusskindHamiltonian::getMesonOperator(), invalid value of n, received " << n << endl;
    exit(666);
  }
  //Check if N is even, otherwise it won't as I assume for the last matrix that the site is even
  if((this->N_fermions%2)!=0)    
  {
    cout << "Error in KogutSusskindHamiltonian::getMesonOperator(), only working for an odd number of sithes, received " << this->N_fermions << endl;
    exit(666);
  }
  
  if(width<=0.0)
    width=2.0;
  
  cout << "Building Meson distribution peaked around " << n << " with width " << width << " and k-vector " << k << endl;
    
  
  mwArray Bosonic_operators(Indices(2,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();   
   for(int i=0; i<this->d_bose; ++i)
   {
     for(int j=0; j<this->d_bose; ++j)
     {
       //Identity
       Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));
       //Lm
       Bosonic_operators.setElement(this->Lm.getElement(Indices(i,j)),Indices(1,i,j));
     }
   }   
  Bosonic_operators.reshape(Indices(2,this->d_bose*this->d_bose));  
  
  mwArray Fermionic_operators(Indices(3,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();    
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma_plus
      Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
      //sigma_minus
      Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
    }
  }  
  Fermionic_operators.reshape(Indices(3,this->d_fermi*this->d_fermi));
  
  //Bond dimension
  int D_bond=3;
  
  //Matrices for the finite automata
  mwArray Fermiodd,Fermieven,Gauge,Fermifirst,Fermilast;
  Fermifirst=mwArray(Indices(1,D_bond,3));	Fermifirst.fillWithZero();
  Fermiodd=mwArray(Indices(D_bond,D_bond,3));	Fermiodd.fillWithZero();
  Fermieven=mwArray(Indices(D_bond,D_bond,3));	Fermieven.fillWithZero();
  Fermilast=mwArray(Indices(D_bond,1,3));	Fermilast.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,2));	Gauge.fillWithZero();
  
  //Variables for gauss distribution and momentum
  double gauss;
  complex_t momentum,phase;
  phase = -I_c*k*M_PIl;
  
  //Set site independent entries
  Fermifirst.setElement(ONE_c,Indices(0,0,0));
  Fermifirst.setElement(ONE_c,Indices(0,1,2));
  
  Fermiodd.setElement(ONE_c,Indices(0,0,0));
  Fermiodd.setElement(ONE_c,Indices(2,2,0));
  Fermiodd.setElement(ONE_c,Indices(0,1,2));
  
  Fermieven.setElement(ONE_c,Indices(0,0,0));
  Fermieven.setElement(ONE_c,Indices(2,2,0));
  Fermieven.setElement(ONE_c,Indices(1,2,1));
  
  Fermilast.setElement(ONE_c,Indices(2,0,0));
  Fermilast.setElement(ONE_c,Indices(1,0,1));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(2,2,0));
  
  //Contract, reshape and set operators
  mwArray res;
  Operator *Local_op;
  bool op_saved=false;
  int pos_op_saved;
  
  res = reshape(Fermifirst,Indices(1*D_bond,3))*Fermionic_operators;
  res.reshape(Indices(1,D_bond,this->d_fermi,this->d_fermi));
  Local_op=new Operator(permute(res,Indices(3,1,4,2)));
  mpo.setOp(0,Local_op,true);
  
  res = reshape(Fermilast,Indices(D_bond*1,3))*Fermionic_operators;
  res.reshape(Indices(D_bond,1,this->d_fermi,this->d_fermi));
  Local_op=new Operator(permute(res,Indices(3,1,4,2)));
  mpo.setOp(this->N_total-1,Local_op,true);
  
  res = reshape(Fermiodd,Indices(D_bond*D_bond,3))*Fermionic_operators;
  res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
  Local_op=new Operator(permute(res,Indices(3,1,4,2)));
  for(int i=3; i<=this->N_fermions; i+=2)
  {
    if(op_saved)
      mpo.setOp(2*i-2,&mpo.getOp(pos_op_saved),false);
    else
    {
      mpo.setOp(2*i-2,Local_op,true);
      op_saved=true;
      pos_op_saved=2*i-2;
    }
  }  
  
  res = reshape(Fermieven,Indices(D_bond*D_bond,3))*Fermionic_operators;
  res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
  Local_op=new Operator(permute(res,Indices(3,1,4,2)));
  op_saved=false;
  for(int i=2; i<this->N_fermions; i+=2)
  {
    if(op_saved)
      mpo.setOp(2*i-2,&mpo.getOp(pos_op_saved),false);
    else
    {
      mpo.setOp(2*i-2,Local_op,true);
      op_saved=true;
      pos_op_saved=2*i-2;
    }
  }
  
  for(int i=1; i<this->N_fermions; i++)
  {
    gauss=exp(-0.5/width/width*pow(i-n,2.0));
    momentum=exp(phase*((double) i));
    //cout << "Factor: " <<  gauss*momentum << endl;
    Gauge.setElement(gauss*momentum,Indices(1,1,1));
    res = reshape(Gauge,Indices(D_bond*D_bond,2))*Bosonic_operators;
    res.reshape(Indices(D_bond,D_bond,this->d_bose,this->d_bose));
    Local_op=new Operator(permute(res,Indices(3,1,4,2)));
    mpo.setOp(2*i-1,Local_op,true);
  } 
}

void KogutSusskindHamiltonian::getMesonAntiMesonOperator(MPO &mpo, int n, double k, double width) const
{
  mpo.initLength(this->N_total);
  
  //Some error checking
  if(n<1 || n>(this->N_fermions-1))
  {
    cout << "Error in KogutSusskindHamiltonian::getMesonAntiMesonOperator(), invalid value of n, received " << n << endl;
    exit(666);
  }
  //Check if N is even, otherwise it won't as I assume for the last matrix that the site is even
  if((this->N_fermions%2)!=0)    
  {
    cout << "Error in KogutSusskindHamiltonian::getMesonAntiMesonOperator(), only working for an odd number of sithes, received " << this->N_fermions << endl;
    exit(666);
  }
  
  if(width<=0.0)
    width=2.0;
  
  cout << "Building Meson Antimeson distribution peaked around " << n << " with width " << width << " and k-vector " << k << endl;
    
  
  mwArray Bosonic_operators(Indices(3,this->d_bose,this->d_bose));
  Bosonic_operators.fillWithZero();   
   for(int i=0; i<this->d_bose; ++i)
   {
     for(int j=0; j<this->d_bose; ++j)
     {
       //Identity
       Bosonic_operators.setElement(this->id_bose.getElement(Indices(i,j)),Indices(0,i,j));
       //Lm
       Bosonic_operators.setElement(this->Lm.getElement(Indices(i,j)),Indices(1,i,j));
       //Lm
       Bosonic_operators.setElement(this->Lp.getElement(Indices(i,j)),Indices(2,i,j));
     }
   }   
  Bosonic_operators.reshape(Indices(3,this->d_bose*this->d_bose));  
  
  mwArray Fermionic_operators(Indices(3,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();    
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma_plus
      Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
      //sigma_minus
      Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
    }
  }  
  Fermionic_operators.reshape(Indices(3,this->d_fermi*this->d_fermi));
  
  //Bond dimension
  int D_bond=4;
  
  //Matrices for the finite automata
  mwArray Fermiodd,Fermieven,Gauge,Fermifirst,Fermilast;
  Fermifirst=mwArray(Indices(1,D_bond,3));	Fermifirst.fillWithZero();
  Fermiodd=mwArray(Indices(D_bond,D_bond,3));	Fermiodd.fillWithZero();
  Fermieven=mwArray(Indices(D_bond,D_bond,3));	Fermieven.fillWithZero();
  Fermilast=mwArray(Indices(D_bond,1,3));	Fermilast.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,3));	Gauge.fillWithZero();
  
  //Variables for gauss distribution and momentum
  double gauss;
  complex_t momentum,phase;
  phase = -I_c*k*M_PIl;
  
  //Set site independent entries
  Fermifirst.setElement(ONE_c,Indices(0,0,0));
  Fermifirst.setElement(ONE_c,Indices(0,1,2));
  
  Fermiodd.setElement(ONE_c,Indices(0,0,0));
  Fermiodd.setElement(ONE_c,Indices(3,3,0));
  Fermiodd.setElement(ONE_c,Indices(0,1,2));
  Fermiodd.setElement(ONE_c,Indices(2,3,2));
  
  Fermieven.setElement(ONE_c,Indices(0,0,0));
  Fermieven.setElement(ONE_c,Indices(3,3,0));
  Fermieven.setElement(ONE_c,Indices(0,2,1));
  Fermieven.setElement(ONE_c,Indices(1,3,1));
  
  Fermilast.setElement(ONE_c,Indices(3,0,0));
  Fermilast.setElement(ONE_c,Indices(1,0,1));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(3,3,0));
  
  //Contract, reshape and set operators
  mwArray res;
  Operator *Local_op;
  bool op_saved=false;
  int pos_op_saved;
  
  res = reshape(Fermifirst,Indices(1*D_bond,3))*Fermionic_operators;
  res.reshape(Indices(1,D_bond,this->d_fermi,this->d_fermi));
  Local_op=new Operator(permute(res,Indices(3,1,4,2)));
  mpo.setOp(0,Local_op,true);
  
  res = reshape(Fermilast,Indices(D_bond*1,3))*Fermionic_operators;
  res.reshape(Indices(D_bond,1,this->d_fermi,this->d_fermi));
  Local_op=new Operator(permute(res,Indices(3,1,4,2)));
  mpo.setOp(this->N_total-1,Local_op,true);
  
  res = reshape(Fermiodd,Indices(D_bond*D_bond,3))*Fermionic_operators;
  res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
  Local_op=new Operator(permute(res,Indices(3,1,4,2)));
  for(int i=3; i<=this->N_fermions; i+=2)
  {
    if(op_saved)
      mpo.setOp(2*i-2,&mpo.getOp(pos_op_saved),false);
    else
    {
      mpo.setOp(2*i-2,Local_op,true);
      op_saved=true;
      pos_op_saved=2*i-2;
    }
  }  
  
  res = reshape(Fermieven,Indices(D_bond*D_bond,3))*Fermionic_operators;
  res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
  Local_op=new Operator(permute(res,Indices(3,1,4,2)));
  op_saved=false;
  for(int i=2; i<this->N_fermions; i+=2)
  {
    if(op_saved)
      mpo.setOp(2*i-2,&mpo.getOp(pos_op_saved),false);
    else
    {
      mpo.setOp(2*i-2,Local_op,true);
      op_saved=true;
      pos_op_saved=2*i-2;
    }
  }
  
  for(int i=1; i<this->N_fermions; i++)
  {
    gauss=exp(-0.5/width/width*pow(i-n,2.0));
    momentum=exp(phase*((double) i));
    //cout << "Factor: " <<  gauss*momentum << endl;
    Gauge.setElement(gauss*momentum,Indices(1,1,1));
    Gauge.setElement(gauss*momentum,Indices(2,2,2));
    res = reshape(Gauge,Indices(D_bond*D_bond,3))*Bosonic_operators;
    res.reshape(Indices(D_bond,D_bond,this->d_bose,this->d_bose));
    Local_op=new Operator(permute(res,Indices(3,1,4,2)));
    mpo.setOp(2*i-1,Local_op,true);
  } 
}

void KogutSusskindHamiltonian::constructMomentumMPO(MPO& mpo) const
{ 
  mpo.initLength(this->N_total);
  
  //Some error checking
  if(this->N_fermions<3)
  {
    cout << "Error in KogutSusskindHamiltonian::constructMomentumMPO(), cannot fit a three fermion operator on a system of size " << this->N_fermions << endl;
    exit(666);
  }
  
  mwArray Fermionic_operators(Indices(4,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();    
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma_plus
      Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
      //sigma_minus
      Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
      //sigma_z
      Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(3,i,j));
    }
  }  
  Fermionic_operators.reshape(Indices(4,this->d_fermi*this->d_fermi));
  
  //Bond dimension
  int D_bond=6;
  
  //Matrices for the finite automata
  mwArray Fermi,Gauge,Fermifirst,Fermilast;
  Fermifirst=mwArray(Indices(1,D_bond,4));	Fermifirst.fillWithZero();
  Fermi=mwArray(Indices(D_bond,D_bond,4));	Fermi.fillWithZero();
  Fermilast=mwArray(Indices(D_bond,1,4));	Fermilast.fillWithZero();
  Gauge=mwArray(Indices(this->d_bose,D_bond,this->d_bose,D_bond));	Gauge.fillWithZero();
  
  Fermifirst.setElement(ONE_c,Indices(0,0,0));
  Fermifirst.setElement(I_c,Indices(0,1,1));
  Fermifirst.setElement(-I_c,Indices(0,2,2));
  
  Fermi.setElement(ONE_c,Indices(0,0,0));
  Fermi.setElement(ONE_c,Indices(5,5,0));
  Fermi.setElement(I_c,Indices(0,1,1));
  Fermi.setElement(ONE_c,Indices(4,5,1));
  Fermi.setElement(-I_c,Indices(0,2,2));
  Fermi.setElement(ONE_c,Indices(3,5,2));
  Fermi.setElement(ONE_c,Indices(1,3,3));
  Fermi.setElement(ONE_c,Indices(2,4,3));
  
  Fermilast.setElement(ONE_c,Indices(5,0,0));
  Fermilast.setElement(ONE_c,Indices(4,0,1));
  Fermilast.setElement(ONE_c,Indices(3,0,2));
    
  for(int i=0; i<this->d_bose; i++)
    for(int k=0; k<D_bond; k++)
      Gauge.setElement(ONE_c,Indices(i,k,i,k));
      
  
  Fermifirst.reshape(Indices(1*D_bond,4));
  Fermi.reshape(Indices(D_bond*D_bond,4));
  Fermilast.reshape(Indices(D_bond*1,4));
  
  Operator *op1=new Operator(permute(reshape(Fermifirst*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi)),Indices(3,1,4,2)));
  Operator *op2=new Operator(permute(reshape(Fermi*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi)),Indices(3,1,4,2)));
  Operator *op3=new Operator(permute(reshape(Fermilast*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi)),Indices(3,1,4,2)));
  Operator *op4=new Operator(Gauge);

  mpo.setOp(0,op1,true);
  mpo.setOp(this->N_total-1,op3,true);
  mpo.setOp(2,op2,true);
  for(int i=4; i<(this->N_total-1); i+=2)
    mpo.setOp(i,&mpo.getOp(2),false);
  mpo.setOp(1,op4,true);
  for(int i=3; i<(this->N_total-1); i+=2)
    mpo.setOp(i,&mpo.getOp(1),false); 
    
  //mpo.exportForMatlab("Pop.m");
}

