#include "KSHamiltonian.h"
#include "Indices.h"
#include "Contractor.h"
#include "splitLocalOp.h"

using namespace std;
using namespace shrt;

KSHamiltonian::KSHamiltonian(int N_, int d_,double mu_,double x_,double alpha_,double lambda_):hamiltonian(2*N_-1)
{
  this->N_sites = N_;  
  this->mu = mu_;
  this->x = x_;
  this->d_link = d_;
  this->d_fermi = 2;
  this->alpha = alpha_;
  this->lambda = lambda_;  
  
  //This model works only for a fixed, even number of bosons on the links, which means an odd link dimension
  if(((this->d_link)%2)==0)
  {
    cout << "Error, model works only for odd link Hilbert space dimension, program will be aborted..." << endl;
    exit(666);
  }
  
  cout << "Constructing Kogut-Susskind-Hamiltonian with " << endl
       << "N_sites   = " << this->N_sites << endl       
       << "mu        = " << this->mu << endl
       << "x         = " << this->x << endl
       << "d         = " << this->d_link << endl
       << "alpha     = " << this->alpha << endl
       << "la        = " << this->lambda << endl;
  
  //Construct the basic operators
  initOperators();
  //Construct hamiltonian MPO
  initHamiltonian();
}

//Destructor
KSHamiltonian::~KSHamiltonian(){
  //Clear Hamiltonian MPO
  hamiltonian.clear();
}

//Update the Hamiltonian MPO, once it has been constructed constructed
void KSHamiltonian::updateHMPO()
{
  cout << "Updating Hamiltonian" << endl;
  this->initHamiltonian();
}

//Provide the basic operators according to the model
void KSHamiltonian::initOperators()
{
  
  //Identity Matrix
  mwArray id_fermi=identityMatrix(this->d_fermi);  
  //sigma_x
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig_x(Indices(this->d_fermi,this->d_fermi),datax);
  //sigma_y
  complex_t datay[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sig_y(Indices(this->d_fermi,this->d_fermi),datay);
  //sigma_z
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sig_z(Indices(this->d_fermi,this->d_fermi),dataz);
  //sigma_+
  mwArray sig_p=0.5*(sig_x+I_c*sig_y);
  //sigma_-
  mwArray sig_m=0.5*(sig_x-I_c*sig_y);
  
  this->FermiOperators=mwArray(Indices(6,this->d_fermi,this->d_fermi));
  for(int i=0; i<this->d_fermi; i++)
    for(int k=0; k<this->d_fermi; k++)
    {
      this->FermiOperators.setElement(id_fermi.getElement(Indices(i,k)),Indices(0,i,k));
      this->FermiOperators.setElement(sig_p.getElement(Indices(i,k)),Indices(1,i,k));
      this->FermiOperators.setElement(sig_m.getElement(Indices(i,k)),Indices(2,i,k));
      this->FermiOperators.setElement(sig_z.getElement(Indices(i,k)),Indices(3,i,k));
      this->FermiOperators.setElement(sig_y.getElement(Indices(i,k)),Indices(4,i,k));
      this->FermiOperators.setElement(sig_x.getElement(Indices(i,k)),Indices(5,i,k));      
    }   
  
  //Generate empty matrices for the links
  mwArray Lp(Indices(this->d_link,this->d_link));		Lp.fillWithZero();
  mwArray Lm(Indices(this->d_link,this->d_link));		Lm.fillWithZero();
  mwArray Lz(Indices(this->d_link,this->d_link));		Lz.fillWithZero();
  mwArray Lzsquare(Indices(this->d_link,this->d_link));		Lzsquare.fillWithZero();
  mwArray id_link=identityMatrix(this->d_link);
   
  //Provide the operators for the links, according to the model
  for(int i=0; i<this->d_link; ++i)
    Lz.setElement(-0.5*(this->d_link-1)+i,0.0,Indices(i,i));

  //Construct Lp and Lm
  double norm_factor = sqrt((this->d_link-1)/2*((this->d_link-1)/2+1));
  complex_t tmp;
  for(int i=0; i<this->d_link; ++i)
    for(int j=0; j<this->d_link; ++j)
      if(i==j+1)
	Lp.setElement(0.0,sqrt(i*(this->d_link-j-1))/norm_factor,Indices(i,j));

  //Get Lm as the hermitian conjugate of Lp
  Lm=Lp;
  Lm.transpose(true);
  //Now get get the square of the L-operator
  Lzsquare = Lz*Lz;
  
  this->LinkOperators=mwArray(Indices(5,this->d_link,this->d_link));
  for(int i=0; i<this->d_link; i++)
    for(int k=0; k<this->d_link; k++)
    {
      this->LinkOperators.setElement(id_link.getElement(Indices(i,k)),Indices(0,i,k));
      this->LinkOperators.setElement(Lp.getElement(Indices(i,k)),Indices(1,i,k));
      this->LinkOperators.setElement(Lm.getElement(Indices(i,k)),Indices(2,i,k));
      this->LinkOperators.setElement(Lz.getElement(Indices(i,k)),Indices(3,i,k));
      this->LinkOperators.setElement(Lzsquare.getElement(Indices(i,k)),Indices(4,i,k));   
    }    
  
  //Identities matrices already reshaped as MPO
  this->id_link_mpo=id_link;
  this->id_fermi_mpo=id_fermi;
  this->id_link_mpo.reshape(Indices(this->d_link,1,this->d_link,1));
  this->id_fermi_mpo.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
  
  //In case I need to see the operators for debugging
  /*cout << "Lp=" << Lp << endl;
  cout << "Lm=" << Lm << endl;
  cout << "Lz=" << Lz << endl;
  cout << "Lz²=" << Lzsquare << endl;
  cout << "sig_x = " << sig_x << endl;
  cout << "sig_y = " << sig_y << endl;
  cout << "sig_z = " << sig_z << endl;
  cout << "sig_p = " << sig_p << endl;
  cout << "sig_m = " << sig_m << endl;*/ 
}

void KSHamiltonian::initHamiltonian()
{
  int D=7;
  
  mwArray FermiFirst(Indices(1,D,6));	FermiFirst.fillWithZero();
  mwArray FermiOdd(Indices(D,D,6));	FermiOdd.fillWithZero();
  mwArray FermiEven(Indices(D,D,6));	FermiEven.fillWithZero();
  mwArray FermiLast(Indices(D,1,6));	FermiLast.fillWithZero();
  
  mwArray BoseOdd(Indices(D,D,5));	BoseOdd.fillWithZero();
  mwArray BoseEven(Indices(D,D,5));	BoseEven.fillWithZero();
  
  //Sign of the last site (which is positive for an even number of sites, because I count form 1 to N)
  double sign=(this->N_sites%2)==0 ? 1.0:-1.0;

  //Identity
  FermiOdd.setElement(ONE_c,Indices(0,0,0));
  FermiOdd.setElement(ONE_c,Indices(D-1,D-1,0));
  FermiOdd.setElement(ONE_c,Indices(3,3,0));
  FermiOdd.setElement(ONE_c*0.5*(this->lambda-this->mu),Indices(0,D-1,0));
  
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(ONE_c*0.5*(this->lambda-this->mu),Indices(0,D-1,0));
  
  FermiLast.setElement(ONE_c*0.5*(this->lambda+sign*this->mu),Indices(0,0,0));
  FermiLast.setElement(ONE_c,Indices(D-1,0,0));
  
  //sp
  FermiOdd.setElement(ONE_c,Indices(0,1,1));
  FermiOdd.setElement(ONE_c,Indices(2,D-1,1));
  
  FermiFirst.setElement(ONE_c,Indices(0,1,1));
  
  FermiLast.setElement(ONE_c,Indices(2,0,1));
  
  //sm
  FermiOdd.setElement(ONE_c,Indices(0,2,2));
  FermiOdd.setElement(ONE_c,Indices(1,D-1,2));
  
  FermiFirst.setElement(ONE_c,Indices(0,2,2));
  
  FermiLast.setElement(ONE_c,Indices(1,0,2));
  
  //sz
  FermiOdd.setElement(-0.5*ONE_c*(this->mu+this->lambda),Indices(0,D-1,3));
  FermiOdd.setElement(-this->lambda*ONE_c,Indices(0,4,3));
  FermiOdd.setElement(ONE_c,Indices(5,D-1,3));
  
  FermiFirst.setElement(-0.5*ONE_c*(this->mu+this->lambda),Indices(0,D-1,3));
  FermiFirst.setElement(-this->lambda*ONE_c,Indices(0,4,3));
  
  FermiLast.setElement(sign*0.5*ONE_c*(this->mu+this->lambda),Indices(0,0,3));
  FermiLast.setElement(ONE_c,Indices(5,0,3));
  
  //As only two entries of even version are different get set it equal to FermiOdd and replace them
  FermiEven = FermiOdd;
  FermiEven.setElement(ONE_c*0.5*(this->lambda+this->mu),Indices(0,D-1,0));
  FermiEven.setElement(0.5*ONE_c*(this->mu+this->lambda),Indices(0,D-1,3));
  
  //Identity
  BoseOdd.setElement(ONE_c,Indices(0,0,0));
  BoseOdd.setElement(ONE_c,Indices(D-1,D-1,0));
  BoseOdd.setElement(this->alpha*this->alpha*ONE_c,Indices(0,D-1,0));
  
  //Lp and Lm
  BoseOdd.setElement(this->x*ONE_c,Indices(1,1,1));
  BoseOdd.setElement(this->x*ONE_c,Indices(2,2,2));
  
  //Lz
  BoseOdd.setElement(2.0*ONE_c*(this->alpha+this->lambda),Indices(0,D-1,3));
  BoseOdd.setElement(-2.0*ONE_c*this->lambda,Indices(0,3,3));
  BoseOdd.setElement(ONE_c,Indices(3,D-1,3));
  BoseOdd.setElement(ONE_c,Indices(4,D-1,3));
  BoseOdd.setElement(ONE_c*this->lambda,Indices(0,5,3));
  
  //Lz^2
  BoseOdd.setElement(ONE_c*(1.0+2.0*this->lambda),Indices(0,D-1,4));
  
  //Same as for the fermionic sites
  BoseEven=BoseOdd;
  BoseEven.setElement(2.0*ONE_c*(this->alpha-this->lambda),Indices(0,D-1,3));
  
  FermiFirst.reshape(Indices(1*D,6)); FermiFirst.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiFirst.reshape(Indices(1,D,this->d_fermi,this->d_fermi));	FermiFirst.permute(Indices(3,1,4,2));
  
  FermiOdd.reshape(Indices(D*D,6)); FermiOdd.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiOdd.reshape(Indices(D,D,this->d_fermi,this->d_fermi));	FermiOdd.permute(Indices(3,1,4,2));
  
  FermiEven.reshape(Indices(D*D,6)); FermiEven.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiEven.reshape(Indices(D,D,this->d_fermi,this->d_fermi));	FermiEven.permute(Indices(3,1,4,2));
  
  FermiLast.reshape(Indices(D*1,6)); FermiLast.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiLast.reshape(Indices(D,1,this->d_fermi,this->d_fermi));	FermiLast.permute(Indices(3,1,4,2));
  
  BoseOdd.reshape(Indices(D*D,5)); BoseOdd.multiplyRight(reshape(this->LinkOperators,Indices(5,this->d_link*this->d_link)));
  BoseOdd.reshape(Indices(D,D,this->d_link,this->d_link));	BoseOdd.permute(Indices(3,1,4,2));
  
  BoseEven.reshape(Indices(D*D,5)); BoseEven.multiplyRight(reshape(this->LinkOperators,Indices(5,this->d_link*this->d_link)));
  BoseEven.reshape(Indices(D,D,this->d_link,this->d_link));	BoseEven.permute(Indices(3,1,4,2));
  
  //Now set the operators
  this->hamiltonian.setOp(0,new Operator(FermiFirst),true);
  this->hamiltonian.setOp(1,new Operator(BoseOdd),true);
  this->hamiltonian.setOp(2*this->N_sites-2,new Operator(FermiLast),true);
  
  bool odd_saved=false, even_saved=false;
  int pos_odd_saved, pos_even_saved;
  
  for(int i=2; i<this->N_sites; i++)
  {
    if((i%2)==0)
    {
      //Even Site
      if(even_saved)
      {
	hamiltonian.setOp(2*i-2,&hamiltonian.getOp(pos_even_saved),false);
	hamiltonian.setOp(2*i-1,&hamiltonian.getOp(pos_even_saved+1),false);
      }
      else
      {
	hamiltonian.setOp(2*i-2,new Operator(FermiEven),true);
	hamiltonian.setOp(2*i-1,new Operator(BoseEven),true);
	even_saved=true;
	pos_even_saved=2*i-2;
      }
    }
    else
    {
      //Odd site
      if(odd_saved)
      {
	hamiltonian.setOp(2*i-2,&hamiltonian.getOp(pos_odd_saved),false);
	hamiltonian.setOp(2*i-1,&hamiltonian.getOp(1),false);
      }
      else
      {
	hamiltonian.setOp(2*i-2,new Operator(FermiOdd),true);
	hamiltonian.setOp(2*i-1,&hamiltonian.getOp(1),false);
	odd_saved=true;
	pos_odd_saved=2*i-2;
      }
    }
  }
}// initHamiltonian


void KSHamiltonian::constructPenaltyMPO(MPO &mpo, double strength) const
{
  int D=5;
  mpo.initLength(2*this->N_sites-1);
  
  if(strength<0.0)
    strength = this->lambda;
  
  //Sign of the last site (which is positive for an even number of sites, because I count form 1 to N)
  double sign=(this->N_sites%2)==0 ? 1.0:-1.0;
  
  mwArray FermiFirst(Indices(1,D,6));	FermiFirst.fillWithZero();
  mwArray FermiOdd(Indices(D,D,6));	FermiOdd.fillWithZero();
  mwArray FermiEven(Indices(D,D,6));	FermiEven.fillWithZero();
  mwArray FermiLast(Indices(D,1,6));	FermiLast.fillWithZero();
  
  mwArray BoseOdd(Indices(D,D,5));	BoseOdd.fillWithZero();
  mwArray BoseEven(Indices(D,D,5));	BoseEven.fillWithZero();
  
  //Identity 
  FermiOdd.setElement(ONE_c,Indices(0,0,0));
  FermiOdd.setElement(ONE_c,Indices(1,1,0));
  FermiOdd.setElement(ONE_c,Indices(D-1,D-1,0));
  FermiOdd.setElement(ONE_c*0.5*strength,Indices(0,D-1,0));
  
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(ONE_c*0.5*strength,Indices(0,D-1,0));
  
  FermiLast.setElement(ONE_c*0.5*strength,Indices(0,0,0));
  FermiLast.setElement(ONE_c,Indices(D-1,0,0));  
  
  //sigma_z
  FermiOdd.setElement(-0.5*strength*ONE_c,Indices(0,D-1,3));
  FermiOdd.setElement(-strength*ONE_c,Indices(0,2,3));
  FermiOdd.setElement(ONE_c,Indices(3,D-1,3));
  
  FermiFirst.setElement(-0.5*strength*ONE_c,Indices(0,D-1,3));
  FermiFirst.setElement(-strength*ONE_c,Indices(0,2,3));
  
  FermiLast.setElement(sign*0.5*strength*ONE_c,Indices(0,0,3));
  FermiLast.setElement(ONE_c,Indices(3,0,3));
  
  //As only a single entry of even version is different get set it equal to FermiOdd and replace them
  FermiEven = FermiOdd;
  FermiEven.setElement(0.5*strength*ONE_c,Indices(0,D-1,3));
  
  //Identity
  BoseOdd.setElement(ONE_c,Indices(0,0,0));
  BoseOdd.setElement(ONE_c,Indices(D-1,D-1,0));
  
  //Lz
  BoseOdd.setElement(-2.0*strength*ONE_c,Indices(0,1,3));
  BoseOdd.setElement(ONE_c,Indices(1,D-1,3));
  BoseOdd.setElement(ONE_c,Indices(2,D-1,3));
  BoseOdd.setElement(strength*ONE_c,Indices(0,3,3));
  BoseOdd.setElement(2.0*strength*ONE_c,Indices(0,D-1,3));
  
  //Lz^2
  BoseOdd.setElement(2.0*strength*ONE_c,Indices(0,D-1,4));
  
  //Same as for the fermionic part
  BoseEven=BoseOdd;
  BoseEven.setElement(-2.0*strength*ONE_c,Indices(0,D-1,3));
  
  FermiFirst.reshape(Indices(1*D,6)); FermiFirst.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiFirst.reshape(Indices(1,D,this->d_fermi,this->d_fermi));	FermiFirst.permute(Indices(3,1,4,2));
  
  FermiOdd.reshape(Indices(D*D,6)); FermiOdd.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiOdd.reshape(Indices(D,D,this->d_fermi,this->d_fermi));	FermiOdd.permute(Indices(3,1,4,2));
  
  FermiEven.reshape(Indices(D*D,6)); FermiEven.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiEven.reshape(Indices(D,D,this->d_fermi,this->d_fermi));	FermiEven.permute(Indices(3,1,4,2));
  
  FermiLast.reshape(Indices(D*1,6)); FermiLast.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiLast.reshape(Indices(D,1,this->d_fermi,this->d_fermi));	FermiLast.permute(Indices(3,1,4,2));
  
  BoseOdd.reshape(Indices(D*D,5)); BoseOdd.multiplyRight(reshape(this->LinkOperators,Indices(5,this->d_link*this->d_link)));
  BoseOdd.reshape(Indices(D,D,this->d_link,this->d_link));	BoseOdd.permute(Indices(3,1,4,2));
  
  BoseEven.reshape(Indices(D*D,5)); BoseEven.multiplyRight(reshape(this->LinkOperators,Indices(5,this->d_link*this->d_link)));
  BoseEven.reshape(Indices(D,D,this->d_link,this->d_link));	BoseEven.permute(Indices(3,1,4,2));
  
  //Now set the operators
  mpo.setOp(0,new Operator(FermiFirst),true);
  mpo.setOp(1,new Operator(BoseOdd),true);
  mpo.setOp(2*this->N_sites-2,new Operator(FermiLast),true);
    
  bool odd_saved=false, even_saved=false;
  int pos_odd_saved, pos_even_saved;
  
  for(int i=2; i<this->N_sites; i++)
  {
    if((i%2)==0)
    {
      //Even Site
      if(even_saved)
      {
	mpo.setOp(2*i-2,&mpo.getOp(pos_even_saved),false);
	mpo.setOp(2*i-1,&mpo.getOp(pos_even_saved+1),false);
      }
      else
      {
	mpo.setOp(2*i-2,new Operator(FermiEven),true);
	mpo.setOp(2*i-1,new Operator(BoseEven),true);
	even_saved=true;
	pos_even_saved=2*i-2;
      }
    }
    else
    {
      //Odd site
      if(odd_saved)
      {
	mpo.setOp(2*i-2,&mpo.getOp(pos_odd_saved),false);
	mpo.setOp(2*i-1,&mpo.getOp(1),false);
      }
      else
      {
	mpo.setOp(2*i-2,new Operator(FermiOdd),true);
	mpo.setOp(2*i-1,&mpo.getOp(1),false);
	odd_saved=true;
	pos_odd_saved=2*i-2;
      }
    }
  }
}//constructPenaltyMPO

void KSHamiltonian::constructGamma5MPO(MPO &mpo, bool prefactor_off) const
{
  int D=4;
  mpo.initLength(2*this->N_sites-1);
  
  double factor=sqrt(this->x)/this->N_sites;
  if(prefactor_off)
    factor=1.0;
  
  mwArray FermiFirst(Indices(1,D,6));	FermiFirst.fillWithZero();
  mwArray FermiOdd(Indices(D,D,6));	FermiOdd.fillWithZero();
  mwArray FermiEven(Indices(D,D,6));	FermiEven.fillWithZero();
  mwArray FermiLast(Indices(D,1,6));	FermiLast.fillWithZero();
  
  mwArray Bose(Indices(D,D,5));		Bose.fillWithZero();  
  
  //Identity
  FermiOdd.setElement(ONE_c,Indices(0,0,0));
  FermiOdd.setElement(ONE_c,Indices(D-1,D-1,0));
  
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  
  FermiLast.setElement(ONE_c,Indices(D-1,0,0));
  
  //sigma_+ and sigma_-
  FermiOdd.setElement(-factor*ONE_c,Indices(0,1,1));
  FermiOdd.setElement(ONE_c,Indices(2,D-1,1));
  
  FermiOdd.setElement(-factor*ONE_c,Indices(0,2,2));
  FermiOdd.setElement(ONE_c,Indices(1,D-1,2));
  
  FermiFirst.setElement(factor*ONE_c,Indices(0,1,1));  
  FermiFirst.setElement(factor*ONE_c,Indices(0,2,2));
  
  FermiLast.setElement(ONE_c,Indices(2,0,1));
  FermiLast.setElement(ONE_c,Indices(1,0,2));  
  
  //As only two entries of even version are different get set it equal to FermiOdd and replace them
  FermiEven=FermiOdd;
  FermiEven.setElement(factor*ONE_c,Indices(0,1,1)); 
  FermiEven.setElement(factor*ONE_c,Indices(0,2,2));
  
  //Identity
  Bose.setElement(ONE_c,Indices(0,0,0));
  Bose.setElement(ONE_c,Indices(D-1,D-1,0));
  
  //Lp and Lm
  Bose.setElement(ONE_c,Indices(1,1,1));
  Bose.setElement(ONE_c,Indices(2,2,2));
  
  //Reshape and contract
  FermiFirst.reshape(Indices(1*D,6)); FermiFirst.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiFirst.reshape(Indices(1,D,this->d_fermi,this->d_fermi));	FermiFirst.permute(Indices(3,1,4,2));
  
  FermiOdd.reshape(Indices(D*D,6)); FermiOdd.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiOdd.reshape(Indices(D,D,this->d_fermi,this->d_fermi));	FermiOdd.permute(Indices(3,1,4,2));
  
  FermiEven.reshape(Indices(D*D,6)); FermiEven.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiEven.reshape(Indices(D,D,this->d_fermi,this->d_fermi));	FermiEven.permute(Indices(3,1,4,2));
  
  FermiLast.reshape(Indices(D*1,6)); FermiLast.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiLast.reshape(Indices(D,1,this->d_fermi,this->d_fermi));	FermiLast.permute(Indices(3,1,4,2));
  
  Bose.reshape(Indices(D*D,5)); Bose.multiplyRight(reshape(this->LinkOperators,Indices(5,this->d_link*this->d_link)));
  Bose.reshape(Indices(D,D,this->d_link,this->d_link));	Bose.permute(Indices(3,1,4,2));
  
  
  //Now set the operators
  mpo.setOp(0,new Operator(FermiFirst),true);
  mpo.setOp(1,new Operator(Bose),true);
  mpo.setOp(2*this->N_sites-2,new Operator(FermiLast),true);
  
  //cout << "Penalty MPO" << endl << mpo << endl;
  
  bool odd_saved=false, even_saved=false;
  int pos_odd_saved, pos_even_saved;
  
  for(int i=1; i<this->N_sites-1; i++)
  {
    if((i%2)==0)
    {
      //Even Site
      if(even_saved)
      {
	mpo.setOp(2*i,&mpo.getOp(pos_even_saved),false);
	mpo.setOp(2*i+1,&mpo.getOp(1),false);
      }
      else
      {
	mpo.setOp(2*i,new Operator(FermiEven),true);
	mpo.setOp(2*i+1,&mpo.getOp(1),false);
	even_saved=true;
	pos_even_saved=2*i;
      }
    }
    else
    {
      //Odd site
      if(odd_saved)
      {
	mpo.setOp(2*i,&mpo.getOp(pos_odd_saved),false);
	mpo.setOp(2*i+1,&mpo.getOp(1),false);
      }
      else
      {
	mpo.setOp(2*i,new Operator(FermiOdd),true);
	mpo.setOp(2*i+1,&mpo.getOp(1),false);
	odd_saved=true;
	pos_odd_saved=2*i;
      }
    }
  }  
}//constructGamma5MPO

void KSHamiltonian::constructGammaAlphaMPO(MPO &mpo, bool prefactor_off) const
{
  int D=2;
  mpo.initLength(2*this->N_sites-1);
  
  double factor=1.0/this->N_sites;
  if(prefactor_off)
    factor=1.0;
  
  mwArray FermiFirst(Indices(1,D,6));	FermiFirst.fillWithZero();
  mwArray Fermi(Indices(D,D,6));	Fermi.fillWithZero();
  mwArray FermiLast(Indices(D,1,6));	FermiLast.fillWithZero();
  
  mwArray Bose(Indices(D,D,5));		Bose.fillWithZero();  
  
  //Identity
  Fermi.setElement(ONE_c,Indices(0,0,0));
  Fermi.setElement(ONE_c,Indices(D-1,D-1,0));  
  
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(factor*this->alpha*ONE_c,Indices(0,D-1,0)); //Additional alpha/N to be compliant with Mari Carmen's conventions
  FermiLast.setElement(ONE_c,Indices(D-1,0,0));
  
  //Identity
  Bose.setElement(ONE_c,Indices(0,0,0));
  Bose.setElement(ONE_c,Indices(D-1,D-1,0));
  Bose.setElement(factor*this->alpha*ONE_c,Indices(0,D-1,0));
  
  //Lz  
  Bose.setElement(factor*ONE_c,Indices(0,D-1,3));
  
  //Reshape and contract
  FermiFirst.reshape(Indices(1*D,6)); FermiFirst.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiFirst.reshape(Indices(1,D,this->d_fermi,this->d_fermi));	FermiFirst.permute(Indices(3,1,4,2));
  
  Fermi.reshape(Indices(D*D,6)); Fermi.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  Fermi.reshape(Indices(D,D,this->d_fermi,this->d_fermi));	Fermi.permute(Indices(3,1,4,2));
  
  FermiLast.reshape(Indices(D*1,6)); FermiLast.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiLast.reshape(Indices(D,1,this->d_fermi,this->d_fermi));	FermiLast.permute(Indices(3,1,4,2));
  
  Bose.reshape(Indices(D*D,5)); Bose.multiplyRight(reshape(this->LinkOperators,Indices(5,this->d_link*this->d_link)));
  Bose.reshape(Indices(D,D,this->d_link,this->d_link));	Bose.permute(Indices(3,1,4,2));
  
  
  //Now set the operators
  mpo.setOp(0,new Operator(FermiFirst),true);
  mpo.setOp(1,new Operator(Bose),true);
  mpo.setOp(2*this->N_sites-2,new Operator(FermiLast),true);
  
  //cout << "Penalty MPO" << endl << mpo << endl;
  
  bool saved=false;
  int pos_saved;
  
  for(int i=1; i<this->N_sites-1; i++)
  {
    if(saved)
    {
      mpo.setOp(2*i,&mpo.getOp(pos_saved),false);
      mpo.setOp(2*i+1,&mpo.getOp(1),false);
    }
    else
    {
      mpo.setOp(2*i,new Operator(Fermi),true);
      mpo.setOp(2*i+1,&mpo.getOp(1),false);
      saved=true;
      pos_saved=2*i;
    }
  }
  
}//constructGammaAlphaMPO

void KSHamiltonian::constructCondensateMPO(MPO &mpo, bool prefactor_off) const
{
  int D=2;
  mpo.initLength(2*this->N_sites-1);
  
  //Sign of the last site (which is negative for an even number of sites, because I count form 0 to N-1)
  double sign=(this->N_sites%2)==0 ? -1.0:1.0;
  
  double factor=sqrt(this->x)/this->N_sites;
  if(prefactor_off)
    factor=1.0;
  
  mwArray FermiFirst(Indices(1,D,6));	FermiFirst.fillWithZero();
  mwArray FermiOdd(Indices(D,D,6));	FermiOdd.fillWithZero();
  mwArray FermiEven(Indices(D,D,6));	FermiEven.fillWithZero();
  mwArray FermiLast(Indices(D,1,6));	FermiLast.fillWithZero();
  
  mwArray Bose(Indices(D,D,5));		Bose.fillWithZero();
  
  //Identity
  FermiOdd.setElement(ONE_c,Indices(0,0,0));
  FermiOdd.setElement(ONE_c,Indices(D-1,D-1,0));
  FermiOdd.setElement(-0.5*factor*ONE_c,Indices(0,D-1,0));
  
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(0.5*factor*ONE_c,Indices(0,D-1,0));
  
  FermiLast.setElement(ONE_c,Indices(D-1,0,0));
  FermiLast.setElement(sign*0.5*factor*ONE_c,Indices(0,0,0));
  
  //sigma_z
  FermiOdd.setElement(-0.5*factor*ONE_c,Indices(0,D-1,3));
  
  FermiFirst.setElement(0.5*factor*ONE_c,Indices(0,D-1,3));
  
  FermiLast.setElement(sign*0.5*factor*ONE_c,Indices(0,0,3));
  
  //As only two entries of even version are different get set it equal to FermiOdd and replace them
  FermiEven=FermiOdd;
  FermiEven.setElement(0.5*factor*ONE_c,Indices(0,D-1,0));
  FermiEven.setElement(0.5*factor*ONE_c,Indices(0,D-1,3));
  
  //Identity
  Bose.setElement(ONE_c,Indices(0,0,0));
  Bose.setElement(ONE_c,Indices(D-1,D-1,0));
  
  //Reshape and contract
  FermiFirst.reshape(Indices(1*D,6)); FermiFirst.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiFirst.reshape(Indices(1,D,this->d_fermi,this->d_fermi));	FermiFirst.permute(Indices(3,1,4,2));
  
  FermiOdd.reshape(Indices(D*D,6)); FermiOdd.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiOdd.reshape(Indices(D,D,this->d_fermi,this->d_fermi));	FermiOdd.permute(Indices(3,1,4,2));
  
  FermiEven.reshape(Indices(D*D,6)); FermiEven.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiEven.reshape(Indices(D,D,this->d_fermi,this->d_fermi));	FermiEven.permute(Indices(3,1,4,2));
  
  FermiLast.reshape(Indices(D*1,6)); FermiLast.multiplyRight(reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)));
  FermiLast.reshape(Indices(D,1,this->d_fermi,this->d_fermi));	FermiLast.permute(Indices(3,1,4,2));
  
  Bose.reshape(Indices(D*D,5)); Bose.multiplyRight(reshape(this->LinkOperators,Indices(5,this->d_link*this->d_link)));
  Bose.reshape(Indices(D,D,this->d_link,this->d_link));	Bose.permute(Indices(3,1,4,2));
  
  
  //Now set the operators
  mpo.setOp(0,new Operator(FermiFirst),true);
  mpo.setOp(1,new Operator(Bose),true);
  mpo.setOp(2*this->N_sites-2,new Operator(FermiLast),true);
  
  
  bool odd_saved=false, even_saved=false;
  int pos_odd_saved, pos_even_saved;
  
  for(int i=1; i<this->N_sites-1; i++)
  {
    if((i%2)==0)
    {
      //Even Site
      if(even_saved)
      {
	mpo.setOp(2*i,&mpo.getOp(pos_even_saved),false);
	mpo.setOp(2*i+1,&mpo.getOp(1),false);
      }
      else
      {
	mpo.setOp(2*i,new Operator(FermiEven),true);
	mpo.setOp(2*i+1,&mpo.getOp(1),false);
	even_saved=true;
	pos_even_saved=2*i;
      }
    }
    else
    {
      //Odd site
      if(odd_saved)
      {
	mpo.setOp(2*i,&mpo.getOp(pos_odd_saved),false);
	mpo.setOp(2*i+1,&mpo.getOp(1),false);
      }
      else
      {
	mpo.setOp(2*i,new Operator(FermiOdd),true);
	mpo.setOp(2*i+1,&mpo.getOp(1),false);
	odd_saved=true;
	pos_odd_saved=2*i;
      }
    }
  }   
}//constructCondensateMPO


void KSHamiltonian::constructMomentumMPO(MPO& mpo) const
{ 
  mpo.initLength(2*this->N_sites-1);
  
  //Some error checking
  if(this->N_sites<3)
  {
    cout << "Error in KSHamiltonian::constructMomentumMPO(), cannot fit a three fermion operator on a system of size " << this->N_sites << endl;
    exit(666);
  }
  //Bond dimension
  int D_bond=6;
  
  //Matrices for the finite automata
  mwArray Fermi,Gauge,Fermifirst,Fermilast;
  Fermifirst=mwArray(Indices(1,D_bond,6));	Fermifirst.fillWithZero();
  Fermi=mwArray(Indices(D_bond,D_bond,6));	Fermi.fillWithZero();
  Fermilast=mwArray(Indices(D_bond,1,6));	Fermilast.fillWithZero();
  Gauge=mwArray(Indices(this->d_link,D_bond,this->d_link,D_bond));	Gauge.fillWithZero();
  
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
    
  for(int i=0; i<this->d_link; i++)
    for(int k=0; k<D_bond; k++)
      Gauge.setElement(ONE_c,Indices(i,k,i,k));
      
  
  Fermifirst.reshape(Indices(1*D_bond,6));
  Fermi.reshape(Indices(D_bond*D_bond,6));
  Fermilast.reshape(Indices(D_bond*1,6));
  
  Operator *op1=new Operator(permute(reshape(Fermifirst*reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)),Indices(1,D_bond,this->d_fermi,this->d_fermi)),Indices(3,1,4,2)));
  Operator *op2=new Operator(permute(reshape(Fermi*reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)),Indices(D_bond,D_bond,this->d_fermi,this->d_fermi)),Indices(3,1,4,2)));
  Operator *op3=new Operator(permute(reshape(Fermilast*reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)),Indices(D_bond,1,this->d_fermi,this->d_fermi)),Indices(3,1,4,2)));
  Operator *op4=new Operator(Gauge);
  
  mpo.setOp(0,op1,true);
  mpo.setOp(2*this->N_sites-2,op3,true);
  mpo.setOp(2,op2,true);
  for(int i=4; i<(2*this->N_sites-2); i+=2)
    mpo.setOp(i,&mpo.getOp(2),false);
  mpo.setOp(1,op4,true);
  for(int i=3; i<(2*this->N_sites-2); i+=2)
    mpo.setOp(i,&mpo.getOp(1),false); 
}//constructMomentumMPO


void KSHamiltonian::constructInitialMPS(MPS &mps, InitState state) const
{ 
  //Different kinds of matrices for different sites
  mwArray spinsiteup(Indices(2,1,1));
  mwArray spinsitedown(Indices(2,1,1));
  mwArray bosesite(Indices(this->d_link,1,1));
  mwArray bosesiteplus(Indices(this->d_link,1,1));
  mwArray bosesiteminus(Indices(this->d_link,1,1));
  
  spinsiteup.fillWithZero();
  spinsitedown.fillWithZero();
  bosesite.fillWithZero();
  bosesiteplus.fillWithZero();
  bosesiteminus.fillWithZero();
  
  spinsiteup.setElement(ONE_c,Indices(0,0,0));
  spinsitedown.setElement(ONE_c,Indices(1,0,0));
  bosesite.setElement(ONE_c,Indices(this->d_link/2,0,0));
  bosesiteplus.setElement(ONE_c,Indices(this->d_link/2+1,0,0));
  bosesiteminus.setElement(ONE_c,Indices(this->d_link/2-1,0,0));
  
  //Adjust the physical dimensions of the MPS, therefore overwrite the old one and create soemthing with matching dimensions
  vector<int> phys_dims(2*this->N_sites-1,this->d_fermi);  
  for(int k=1; k<(2*this->N_sites-1); k+=2)
    phys_dims[k]=this->d_link;   
  
  mps.clear();
  mps=MPS(2*this->N_sites-1,1,phys_dims);
  
  int count=0;
  if(state == is_downup)
  {
    cout << "Preparing state |d>|0>|u>..." << endl;
    for(int i=0; i<2*this->N_sites-1; i++)
    {
      if((i%2)==0)
      {
	//Fermionic sites
	if((count%2)==0)
	{
	  //initmps.setA(i,spinsiteup);
	  mps.setA(i,spinsitedown);
	}
	else
	{
	  //initmps.setA(i,spinsitedown);
	  mps.setA(i,spinsiteup);
	}
	count ++;
      }
      else
      {
	//Bosonic sites
	//initmps.setA(i,bosesite);
	mps.setA(i,bosesite);
      }
    }
  }
  else if(state == is_updown)
  {
    cout << "Preparing state |u>|0>|d>..." << endl;
    for(int i=0; i<2*this->N_sites-1; i++)
    {
      if((i%2)==0)
      {
	//Fermionic sites
	if((count%2)==0)
	{
	  //initmps.setA(i,spinsitedown);
	  mps.setA(i,spinsiteup);
	}
	else
	{
	  //initmps.setA(i,spinsiteup);
	  mps.setA(i,spinsitedown);
	}
	count ++;
      }
      else
      {
	//Bosonic sites
	//initmps.setA(i,bosesite);
	mps.setA(i,bosesite);
      }
    }
  }
  else if(state == is_upup)
  {
    cout << "Preparing state |u>|0>|u>..." << endl;
    for(int i=0; i<2*this->N_sites-1; i++)
    {
      if((i%2)==0)
      {
	//Fermionic sites
	mps.setA(i,spinsiteup);
      }
      else
      {
	//Bosonic sites
	//initmps.setA(i,bosesite);
	mps.setA(i,bosesite);
      }
    }
  }
  else if(state == is_downdown)
  {
    cout << "Preparing state |d>|0>|d>..." << endl;
    for(int i=0; i<2*this->N_sites-1; i++)
    {
      if((i%2)==0)
      {
	//Fermionic sites
	mps.setA(i,spinsitedown);
      }
      else
      {
	//Bosonic sites
	//initmps.setA(i,bosesite);
	mps.setA(i,bosesite);
      }
    }
  }
  else if(state == is_bgfield)
  {
    //Little error check
    if((this->N_sites%2)!=0)
    {
      cout << "State |d>|-1>|d>|-1>|u>|-1>|d>|-1>|u>...|-1>|u> can only be constructed for an even number of sitest" << endl;
      cout << "Program will be aborted..." << endl;
      exit(666);
    }
    cout << "Constructing state |d>|-1>|d>|-1>|u>|-1>|d>|-1>|u>...|-1>|u>" << endl;
    
    mps.setA(0,spinsitedown);
    mps.setA(1,bosesiteminus);
    mps.setA(2*this->N_sites-2,spinsiteup);
    
    for(int i=2; i<this->N_sites; i++)
    {
      if((i%2)==0)
	mps.setA(2*i-2,spinsitedown);
      else
	mps.setA(2*i-2,spinsiteup);
      mps.setA(2*i-1,bosesiteminus);
    }
  }
  else if(state == is_bgfield_superpos)
  {
    //Little error check
    if((this->N_sites%2)!=0)
    {
      cout << "State can only be constructed for an even number of sitest" << endl;
      cout << "Program will be aborted..." << endl;
      exit(666);
    }    
    cout << "Constructing superposition" << endl;
    cout.flush();
    
    //Now I have bond dimension two, so adjust MPS accordingly
    mps.clear();
    mps=MPS(2*this->N_sites-1,2,phys_dims);
      
    mwArray First,Odd,Even,Last;
    mwArray Gauge(Indices(this->d_link,2,2));
    First=mwArray(Indices(this->d_fermi,1,2));	First.fillWithZero();
    Odd=mwArray(Indices(this->d_fermi,2,2));	Odd.fillWithZero();
    Even=mwArray(Indices(this->d_fermi,2,2));	Even.fillWithZero();
    Last=mwArray(Indices(this->d_fermi,2,1));	Last.fillWithZero();
    
    //Set elements
    First.setElement(ONE_c,Indices(0,0,0));
    First.setElement(ONE_c,Indices(1,0,1));
    
    Odd.setElement(ONE_c,Indices(0,0,0));
    Odd.setElement(ONE_c,Indices(0,1,1));
    
    Even.setElement(ONE_c,Indices(1,0,0));
    Even.setElement(ONE_c,Indices(1,1,1));
    
    Last.setElement(ONE_c,Indices(0,1,0));
    Last.setElement(ONE_c,Indices(1,0,0));
    
    Gauge.setElement(ONE_c,Indices(this->d_link/2,0,0));
    Gauge.setElement(ONE_c,Indices(this->d_link/2-1,1,1));
    
    mps.setA(0,First);
    mps.setA(1,Gauge);
    mps.setA(2*this->N_sites-2,Last);
    
    for(int i=2; i<this->N_sites; i++)
    {
      if((i%2)==0)
	mps.setA(2*i-2,Even);
      else
	mps.setA(2*i-2,Odd);
      mps.setA(2*i-1,Gauge);
    }
    
  }
  else
    cout << "State not (yet) supported" << endl;  
}

void KSHamiltonian::constructChargeMPO(MPO &mpo) const
{
  //Make sure given MPO is empty
  mpo.initLength(2*this->N_sites-1);
  
  cout << "Constructing charge MPO" << endl;
  
  //Bond dimension for the MPO Matrices
  int D_bond = 2;  
  
  //Matrices
  mwArray Link(Indices(D_bond,D_bond,5));		Link.fillWithZero();
  mwArray Fermi_first(Indices(1,D_bond,6));		Fermi_first.fillWithZero();
  mwArray Fermi_last(Indices(D_bond,1,6));		Fermi_last.fillWithZero();
  mwArray Fermi_even(Indices(D_bond,D_bond,6));		Fermi_even.fillWithZero();
  mwArray Fermi_odd(Indices(D_bond,D_bond,6));		Fermi_odd.fillWithZero();
  
  //Set non-zero entries
  Link.setElement(ONE_c,Indices(0,0,0));
  Link.setElement(ONE_c,Indices(1,1,0));
  
  Fermi_first.setElement(ONE_c,Indices(0,0,0));
  Fermi_first.setElement(-0.5*ONE_c,Indices(0,1,0));
  Fermi_first.setElement(0.5*ONE_c,Indices(0,1,3));  
  
  Fermi_last.setElement(0.5*pow(-1.0,this->N_sites)*ONE_c,Indices(0,0,0));
  Fermi_last.setElement(ONE_c,Indices(1,0,0));
  Fermi_last.setElement(0.5*ONE_c,Indices(0,0,3));
  
  Fermi_even.setElement(ONE_c,Indices(0,0,0));
  Fermi_even.setElement(ONE_c,Indices(1,1,0));
  Fermi_even.setElement(0.5*ONE_c,Indices(0,1,0));
  Fermi_even.setElement(0.5*ONE_c,Indices(0,1,3));
  
  Fermi_odd.setElement(ONE_c,Indices(0,0,0));
  Fermi_odd.setElement(ONE_c,Indices(1,1,0));
  Fermi_odd.setElement(-0.5*ONE_c,Indices(0,1,0));
  Fermi_odd.setElement(0.5*ONE_c,Indices(0,1,3));
  
  //Reshape
  Link.reshape(Indices(D_bond*D_bond,5));
  Fermi_first.reshape(Indices(1*D_bond,6));
  Fermi_last.reshape(Indices(D_bond*1,6));
  Fermi_even.reshape(Indices(D_bond*D_bond,6));
  Fermi_odd.reshape(Indices(D_bond*D_bond,6));
  
  //Contract and permute
  Fermi_first=permute(reshape(Fermi_first*reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)),Indices(1,D_bond,this->d_fermi,this->d_fermi)),Indices(3,1,4,2));
  Fermi_last=permute(reshape(Fermi_last*reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)),Indices(D_bond,1,this->d_fermi,this->d_fermi)),Indices(3,1,4,2));
  Fermi_odd=permute(reshape(Fermi_odd*reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)),Indices(D_bond,D_bond,this->d_fermi,this->d_fermi)),Indices(3,1,4,2));  
  Fermi_even=permute(reshape(Fermi_even*reshape(this->FermiOperators,Indices(6,this->d_fermi*this->d_fermi)),Indices(D_bond,D_bond,this->d_fermi,this->d_fermi)),Indices(3,1,4,2));   
  Link=permute(reshape(Link*reshape(this->LinkOperators,Indices(5,this->d_link*this->d_link)),Indices(D_bond,D_bond,this->d_link,this->d_link)),Indices(3,1,4,2));
  
  //Flags, if operators were already saved
  bool odd_saved=false, even_saved=false;
  int pos_odd_saved,pos_even_saved;
  
  //Now fill the MPO
  mpo.setOp(0,new Operator(Fermi_first),true);
  mpo.setOp(1,new Operator(Link),true);
  mpo.setOp(2*this->N_sites-2,new Operator(Fermi_last),true);
  for(int i=2; i<this->N_sites; i++)
  {
    if((i%2)==0)
    {
      //Even Site
      if(even_saved)
	mpo.setOp(2*i-2,&mpo.getOp(pos_even_saved),false);
      else
      {
	mpo.setOp(2*i-2,new Operator(Fermi_even),true);
	even_saved=true;
	pos_even_saved=2*i-2;
      }
    }
    else
    {
      //Odd site
      if(odd_saved)
	mpo.setOp(2*i-2,&mpo.getOp(pos_odd_saved),false);
      else
      {
	mpo.setOp(2*i-2,new Operator(Fermi_odd),true);
	odd_saved=true;
	pos_odd_saved=2*i-2;
      }
    }
    //Link
    mpo.setOp(2*i-1,&mpo.getOp(1),false);
  }  
}//constructChargeMPO
  
void KSHamiltonian::constructUoddMPO(MPO &Uodd, double dt, bool imag_time)
{  
  complex_t prefactor=-1.0*dt*I_c;
  if(imag_time)
    prefactor = -I_c*prefactor;
    
  //Since I want to split id+s3 on every site in 0.5*(id+s3) + 0.5*(id+s3) and take one term to the even MPO and one to the odd MPO I get a problem at the first site and the last site, because one term is missing. To avoid this I have to calculate an extra kronecker product for the first site and for the site before the last site. Graphically the problem looks like this
  //x--o--x--o--x--o--x
  //|  |  |  |  |  |  |
  //-------  |  -------
  //|  |  -------  |  |
  //Note: That this scheme works, the chain has to have an even number of spin sites, otherwise it could happen, that this has to be done for an even side MPO

  if((this->N_sites%2)!=0)
    cout << "!!! Warning in KSHamiltonian::constructUoddMPO(), time evolution MPO will be wrong for an odd number of fermions !!!" << endl;
      
  //sigma_plus Op1 sigma_minus
  mwArray hopping1=kron(this->FermiOperators.subArray(Indices(1,-1,-1)),kron(this->LinkOperators.subArray(Indices(1,-1,-1)),this->FermiOperators.subArray(Indices(2,-1,-1))));
  //sigma_minus Op2 sigma_plus
  mwArray hopping2=kron(this->FermiOperators.subArray(Indices(2,-1,-1)),kron(this->LinkOperators.subArray(Indices(2,-1,-1)),this->FermiOperators.subArray(Indices(1,-1,-1))));
  //(L+alpha)^2
  mwArray lz_a_square = this->LinkOperators.subArray(Indices(3,-1,-1)) + this->alpha*identityMatrix(this->d_link);
  lz_a_square = lz_a_square*lz_a_square;
  mwArray lterm=kron(this->FermiOperators.subArray(Indices(0,-1,-1)),kron(lz_a_square,this->FermiOperators.subArray(Indices(0,-1,-1))));
  // 0.5*idpsz id id
  mwArray idpsz=this->FermiOperators.subArray(Indices(0,-1,-1))+this->FermiOperators.subArray(Indices(3,-1,-1));
  mwArray idpszleft=kron(0.5*idpsz,kron(this->LinkOperators.subArray(Indices(0,-1,-1)),this->FermiOperators.subArray(Indices(0,-1,-1))));
  // id id 0.5*idpsz
  mwArray idpszright=kron(this->FermiOperators.subArray(Indices(0,-1,-1)),kron(this->LinkOperators.subArray(Indices(0,-1,-1)),0.5*idpsz));
  // idpsz id id
  mwArray firstidpszleft=kron(-1.0*idpsz,kron(this->LinkOperators.subArray(Indices(0,-1,-1)),this->FermiOperators.subArray(Indices(0,-1,-1))));;  
  // id id idpsz
  mwArray lastidpszright=kron(this->FermiOperators.subArray(Indices(0,-1,-1)),kron(this->LinkOperators.subArray(Indices(0,-1,-1)),idpsz));  
  
  vector<mwArray> matrices;
  vector<int> dims;  
  dims.push_back(this->d_fermi);
  dims.push_back(this->d_link);
  dims.push_back(this->d_fermi);
  
  //Calculate the tensor products for the exponentials, for FIRST ODD
  mwArray arg_odd,exp_odd;
  //Attention: The function for the Kronecker product does not the same as the matlab routine. If A=(a_ij) and B=(b_kl) then kron(A,B) gives B \otimes A (exactly the other way round than the matlab Kronecker product
  arg_odd = this->x*hopping1 + this->x*hopping2 + 0.5*this->mu*firstidpszleft + 0.5*this->mu*idpszright + lterm;
  wrapper::expm(arg_odd,exp_odd,prefactor);  
  splitLocalOp (exp_odd,dims,matrices);
  
  //Three matrices for the odd MPO
  mwArray W1_odd_first_site,W2_odd_first_site,W3_odd_first_site;
  
  W1_odd_first_site=matrices[0];
  W2_odd_first_site=matrices[1];
  W3_odd_first_site=matrices[2];
  
  //Calculate the tensor products for the exponentials, for the site BEFORE THE LAST SITE
  arg_odd.clear(); exp_odd.clear(); matrices.clear();
  arg_odd = this->x*hopping1 + this->x*hopping2 - 0.5*this->mu*idpszleft + 0.5*this->mu*lastidpszright + lterm;  
  wrapper::expm(arg_odd,exp_odd,prefactor);
  splitLocalOp (exp_odd,dims,matrices);
  
  //Three matrices for the odd MPO
  mwArray W1_odd_last_site,W2_odd_last_site,W3_odd_last_site;
  
  W1_odd_last_site=matrices[0];
  W2_odd_last_site=matrices[1];
  W3_odd_last_site=matrices[2]; 
 
  //Calculate the tensor products for the exponentials, now for the rest of the odd part
  arg_odd.clear(); exp_odd.clear(); matrices.clear();
  arg_odd = this->x*hopping1 + this->x*hopping2 - 0.5*this->mu*idpszleft + 0.5*this->mu*idpszright + lterm;  
  wrapper::expm(arg_odd,exp_odd,prefactor);
  splitLocalOp (exp_odd,dims,matrices);
  
  //Three matrices for the odd MPO
  mwArray W1_odd,W2_odd,W3_odd;
  
  W1_odd=matrices[0];
  W2_odd=matrices[1];
  W3_odd=matrices[2];
  
  //Now fill the odd MPO
  Uodd.initLength(2*this->N_sites-1);
  
  bool first=true;  
  for(int i=0; i<2*this->N_sites-1; i +=4)
  {
    if(i==0)
    {
      //First site of the chain (which is equal to the first fermionic site)
      Uodd.setOp(i,new Operator(W1_odd_first_site),true);
      Uodd.setOp(i+1,new Operator(W2_odd_first_site),true);
      Uodd.setOp(i+2,new Operator(W3_odd_first_site),true);
    }
    else if(i==(2*this->N_sites-1-3))
    {
      //Site before the last fermionic site
      Uodd.setOp(i,new Operator(W1_odd_last_site),true);
      Uodd.setOp(i+1,new Operator(W2_odd_last_site),true);
      Uodd.setOp(i+2,new Operator(W3_odd_last_site),true);
    }      
    else if(first && (i!=0) && i!=(2*this->N_sites-1-3))
    {
      Uodd.setOp(i,new Operator(W1_odd),true);
      Uodd.setOp(i+1,new Operator(W2_odd),true);
      Uodd.setOp(i+2,new Operator(W3_odd),true);
      first = false;
    }
    else
    {
      Uodd.setOp(i,&Uodd.getOp(i-4),false);
      Uodd.setOp(i+1,&Uodd.getOp(i-3),false);
      Uodd.setOp(i+2,&Uodd.getOp(i-2),false);
    }
  }
  //Set Identities
  first = true;
  for(int i=3; i<2*this->N_sites-1; i+=4)
  {
    if(first)
    {
      Uodd.setOp(i,new Operator(this->id_link_mpo),true);
      first=false;
    }
    else
      Uodd.setOp(i,&Uodd.getOp(i-4),false);
  }
}//constructUoddMPO

void KSHamiltonian::constructUevenMPO(MPO &Ueven, double dt, bool imag_time)
{
  complex_t prefactor=-1.0*dt*I_c;
  if(imag_time)
    prefactor = -I_c*prefactor;
  
  //Since I want to split id+s3 on every site in 0.5*(id+s3) + 0.5*(id+s3) and take one term to the even MPO and one to the odd MPO I get a problem at the first site and the last site, because one term is missing. To avoid this I have to calculate an extra MPO for the first site and an extra MPO for the site before the last site. Graphically the problem looks like this
  //x--o--x--o--x--o--x
  //|  |  |  |  |  |  |
  //-------  |  -------
  //|  |  -------  |  |
  //Note: That this scheme works, the chain has to have an even number of spin sites, otherwise it could happen, that this hast to be done for an even side MPO
  
  //sigma_plus Op1 sigma_minus
  mwArray hopping1=kron(this->FermiOperators.subArray(Indices(1,-1,-1)),kron(this->LinkOperators.subArray(Indices(1,-1,-1)),this->FermiOperators.subArray(Indices(2,-1,-1))));
  //sigma_minus Op2 sigma_plus
  mwArray hopping2=kron(this->FermiOperators.subArray(Indices(2,-1,-1)),kron(this->LinkOperators.subArray(Indices(2,-1,-1)),this->FermiOperators.subArray(Indices(1,-1,-1))));
  //mwArray id L^2  id
  mwArray lz_a_square = this->LinkOperators.subArray(Indices(3,-1,-1)) + this->alpha*identityMatrix(this->d_link);
  lz_a_square = lz_a_square*lz_a_square;
  mwArray lterm=kron(this->FermiOperators.subArray(Indices(0,-1,-1)),kron(lz_a_square,this->FermiOperators.subArray(Indices(0,-1,-1))));
  // 0.5*idpsz id id
  mwArray idpsz=this->FermiOperators.subArray(Indices(0,-1,-1))+this->FermiOperators.subArray(Indices(3,-1,-1));
  mwArray idpszleft=kron(0.5*idpsz,kron(this->LinkOperators.subArray(Indices(0,-1,-1)),this->FermiOperators.subArray(Indices(0,-1,-1))));
  // id id 0.5*idpsz
  mwArray idpszright=kron(this->FermiOperators.subArray(Indices(0,-1,-1)),kron(this->LinkOperators.subArray(Indices(0,-1,-1)),0.5*idpsz));
  
  if((this->N_sites%2)!=0)
    cout << "!!! Warning in KSHamiltonian::constructUevenMPO(), time evolution MPO will be wrong for an odd number of fermions !!!" << endl;
    
  //Three matrices for the even MPO
  mwArray W1_even,W2_even,W3_even;
  
  vector<mwArray> matrices;
  vector<int> dims;  
  dims.push_back(this->d_fermi);
  dims.push_back(this->d_link);
  dims.push_back(this->d_fermi);
  
  //Calculate the tensorproducts for the exponentials, now for even part
  mwArray arg_even,exp_even;
  arg_even = this->x*hopping1 + this->x*hopping2 + 0.5*this->mu*idpszleft - 0.5*this->mu*idpszright + lterm;  
  wrapper::expm(arg_even,exp_even,prefactor);
  splitLocalOp (exp_even,dims,matrices);
  
  W1_even=matrices[0];
  W2_even=matrices[1];
  W3_even=matrices[2];
  
  //Now fill the even MPO
  Ueven.initLength(2*this->N_sites-1);
  
  bool first = true;
  for(int i=2; i<2*this->N_sites-1-2; i +=4)
  {   
    if(first)
    {
      Ueven.setOp(i,new Operator(W1_even),true);
      Ueven.setOp(i+1,new Operator(W2_even),true);
      Ueven.setOp(i+2,new Operator(W3_even),true);
      first = false;
    }
    else
    {
      Ueven.setOp(i,&Ueven.getOp(i-4),false);
      Ueven.setOp(i+1,&Ueven.getOp(i-3),false);
      Ueven.setOp(i+2,&Ueven.getOp(i-2),false);
    }
  }
  //Set link identities
  first = true;
  for(int i=1; i<2*this->N_sites-1; i+=4)
  {
    if(first)
    {
      Ueven.setOp(i,new Operator(this->id_link_mpo),true);
      first=false;
    }
    else
      Ueven.setOp(i,&Ueven.getOp(i-4),false);
  }
  
  //First and last fermonic identity
  Ueven.setOp(0,new Operator(this->id_fermi_mpo),true);
  Ueven.setOp(2*this->N_sites-2,&Ueven.getOp(0));
}//constructUevenMPO
