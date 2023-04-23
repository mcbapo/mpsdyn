/**
   \file ReducedSU2Hamiltonian.cpp
   Implementation of the truncated version of the SU(2) Hamiltonian from Erez, Benni and Ignacio
   
   \author Stefan Kühn
   \date 15/08/2015
*/

#include "ReducedSU2Hamiltonian.h"
#include "Indices.h"
#include "Contractor.h"
#include "splitLocalOp.h"
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>

using namespace std;
using namespace shrt;


//Standard Constructor
ReducedSU2Hamiltonian::ReducedSU2Hamiltonian(int N_,double epsilon_,double mu_,double g_,int Nmax_):hamil(2*N_-1)
{
  //Number of fermions
  this->N_fermions=N_;
  //Total number of sites (2*N_fermions-1)
  this->N_total=2*N_-1;  
  //Hopping constant
  this->epsilon=epsilon_;  
  //Mass
  this->mu=mu_;
  //Coupling
  this->g=g_;
  //Dimension of the link Hilbert spaces
  this->d_link=Nmax_;
  //Dimension of the fermionic Hilbertspaces
  this->d_fermi=3;

  //Construct the basic operators
  initOperators();
  //Construct the actual hamiltonian MPO
  initHamiltonian();
}

//Constructor with penalty
ReducedSU2Hamiltonian::ReducedSU2Hamiltonian(int N_,double epsilon_,double mu_,double g_,int Nmax_, double lambda_):hamil(2*N_-1)
{
  //Number of fermions
  this->N_fermions=N_;
  //Total number of sites (2*N_fermions-1)
  this->N_total=2*N_-1;  
  //Hopping constant
  this->epsilon=epsilon_;  
  //Mass
  this->mu=mu_;
  //Coupling
  this->g=g_;
  //Dimension of the link Hilbert spaces
  this->d_link=Nmax_;
  //Dimension of the fermionic Hilbertspaces
  this->d_fermi=3;
  //Penalty strength
  this->lambda=lambda_;

  //Construct the basic operators
  initOperators();
  //Construct the actual hamiltonian MPO
  initHamiltonianPenalty();
}



//Destructor
ReducedSU2Hamiltonian::~ReducedSU2Hamiltonian(){
  //Clear Hamiltonian MPO
  hamil.clear();  
  //Clear all operators

}

//Provide the basic operator matrices used for MPOs
void ReducedSU2Hamiltonian::initOperators(void)
{    

  //Fermionic operators
  this->O1=mwArray(Indices(this->d_fermi,this->d_fermi));	this->O1.fillWithZero();	this->O1dagger=mwArray(Indices(this->d_fermi,this->d_fermi));	this->O1dagger.fillWithZero();
  this->O2=mwArray(Indices(this->d_fermi,this->d_fermi));	this->O2.fillWithZero();	this->O2dagger=mwArray(Indices(this->d_fermi,this->d_fermi));	this->O2dagger.fillWithZero();
  this->Noc=mwArray(Indices(this->d_fermi,this->d_fermi));	this->Noc.fillWithZero();
  
  this->O1.setElement(ONE_c,Indices(0,1));	this->O1dagger=this->O1; this->O1dagger.transpose(true);
  this->O2.setElement(ONE_c,Indices(1,2));	this->O2dagger=this->O2; this->O2dagger.transpose(true);
  this->Noc.setElement(ONE_c,Indices(1,1));	this->Noc.setElement(2.0*ONE_c,Indices(2,2));  
  
  //Gauge operators
  this->Jsquared=mwArray(Indices(this->d_link,this->d_link));	this->Jsquared.fillWithZero();
  this->Lvalop=mwArray(Indices(this->d_link,this->d_link));	this->Lvalop.fillWithZero();
  this->Lpjdep=mwArray(Indices(this->d_link,this->d_link));	this->Lpjdep.fillWithZero();
  this->Lmjdep=mwArray(Indices(this->d_link,this->d_link));	this->Lmjdep.fillWithZero();
  this->Lm=mwArray(Indices(this->d_link,this->d_link));		this->Lm.fillWithZero();
  this->Lp=mwArray(Indices(this->d_link,this->d_link));		this->Lp.fillWithZero();
  
  double lval;
  for (int l=0;l<this->d_link; l++)
  {
    lval=((double)l)/2.0;
    this->Jsquared.setElement(lval*(lval+1)*ONE_c,Indices(l,l));
    this->Lvalop.setElement(lval*ONE_c,Indices(l,l));
    if((l-1)>=0)
    {
      this->Lmjdep.setElement(-sqrt(2.0*lval/(1.0+2.0*lval))*ONE_c,Indices(l-1,l));
      this->Lm.setElement(ONE_c,Indices(l-1,l));
    }
    if((l+1)<this->d_link)
    {
      this->Lpjdep.setElement(sqrt(2.0*(1.0+lval)/(1.0+2.0*lval))*ONE_c,Indices(l+1,l));
      this->Lp.setElement(ONE_c,Indices(l+1,l));
    }
  }  
  
  //Get the adjoints
  this->Lpjdepdagger=this->Lpjdep;	this->Lpjdepdagger.transpose(true);
  this->Lmjdepdagger=this->Lmjdep;	this->Lmjdepdagger.transpose(true);
  this->Lmdagger=this->Lm;		this->Lmdagger.transpose(true);
  this->Lpdagger=this->Lp;		this->Lpdagger.transpose(true); 
  
  //Identity for fermionic site
  this->id_fermi=identityMatrix(this->d_fermi);  
  //Identity for a link
  this->id_link=identityMatrix(this->d_link);
  
  //I need often identities to construct a single body operator, therefore I prepare those and save them
  this->id_link_mpo=this->id_link;
  this->id_fermi_mpo=this->id_fermi;
  this->id_link_mpo.reshape(Indices(this->d_link,1,this->d_link,1));
  this->id_fermi_mpo.reshape(Indices(this->d_fermi,1,this->d_fermi,1));    
  
  //For debugging
  /*cout << "O1      = " << this->O1 << endl;
  cout << "O2      = " << this->O2 << endl;
  cout << "O1d     = " << this->O1dagger << endl;
  cout << "O2d     = " << this->O2dagger << endl;
  cout << "Noc     = " << this->Noc << endl;
  cout << "Lpjdep  = " << this->Lpjdep << endl;
  cout << "Lmjdep  = " << this->Lmjdep << endl;
  cout << "Lp      = " << this->Lp << endl;
  cout << "Lm      = " << this->Lm << endl;
  cout << "Lpjdepd = " << this->Lpjdepdagger << endl;
  cout << "Lmjdepd = " << this->Lmjdepdagger << endl;
  cout << "Lpd     = " << this->Lpdagger << endl;
  cout << "Lmd     = " << this->Lmdagger << endl;
  cout << "J²      = " << this->Jsquared << endl;*/
}

//Build the Hamiltonian MPO
void ReducedSU2Hamiltonian::initHamiltonian()
{
  cout << "Building Hamiltonian MPO with" << endl
       << "epsilon: " << this->epsilon << endl
       << "mu:      " << this->mu << endl
       << "g:       " << this->g << endl
       << "dlink:   " << this->d_link << endl;
  
  //Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(10,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //Lp/Lm site dependent
      Gauge_operators.setElement(this->Lmjdep.getElement(Indices(i,j)),Indices(1,i,j));
      Gauge_operators.setElement(this->Lpjdep.getElement(Indices(i,j)),Indices(2,i,j));
      Gauge_operators.setElement(this->Lmjdepdagger.getElement(Indices(i,j)),Indices(3,i,j));
      Gauge_operators.setElement(this->Lpjdepdagger.getElement(Indices(i,j)),Indices(4,i,j));
      //Lp/Lm site independent
      Gauge_operators.setElement(this->Lm.getElement(Indices(i,j)),Indices(5,i,j));
      Gauge_operators.setElement(this->Lp.getElement(Indices(i,j)),Indices(6,i,j));
      Gauge_operators.setElement(this->Lmdagger.getElement(Indices(i,j)),Indices(7,i,j));
      Gauge_operators.setElement(this->Lpdagger.getElement(Indices(i,j)),Indices(8,i,j));
      //Jsquared
      Gauge_operators.setElement(this->Jsquared.getElement(Indices(i,j)),Indices(9,i,j));
     }
  }   
  //cout << "Gauge_operators" << Gauge_operators <<endl;
  Gauge_operators.reshape(Indices(10,this->d_link*this->d_link));  
  //cout << "Bosonic operators done" << endl;
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(6,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();    
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //O1/O2
      Fermionic_operators.setElement(this->O1.getElement(Indices(i,j)),Indices(1,i,j));
      Fermionic_operators.setElement(this->O2.getElement(Indices(i,j)),Indices(2,i,j));
      Fermionic_operators.setElement(this->O1dagger.getElement(Indices(i,j)),Indices(3,i,j));
      Fermionic_operators.setElement(this->O2dagger.getElement(Indices(i,j)),Indices(4,i,j));
      //Occupation number operator
      Fermionic_operators.setElement(this->Noc.getElement(Indices(i,j)),Indices(5,i,j));
      
    }
  }  
  Fermionic_operators.reshape(Indices(6,this->d_fermi*this->d_fermi));
  //cout << "Fermionic operators done" << endl;
  
  //Now prepare the matrices
  int D_bond=10;
  mwArray FermiFirst,FermiEven,FermiOdd,FermiLast,Gauge;
  FermiFirst=mwArray(Indices(1,D_bond,6));	FermiFirst.fillWithZero();
  FermiEven=mwArray(Indices(D_bond,D_bond,6));	FermiEven.fillWithZero();
  FermiOdd=mwArray(Indices(D_bond,D_bond,6));	FermiOdd.fillWithZero();
  FermiLast=mwArray(Indices(D_bond,1,6));	FermiLast.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,10));	Gauge.fillWithZero();
  
  //FermiOdd
  FermiOdd.setElement(ONE_c,Indices(0,0,0));
  FermiOdd.setElement(ONE_c,Indices(D_bond-1,D_bond-1,0));
  
  FermiOdd.setElement(this->epsilon*ONE_c,Indices(0,1,1));
  FermiOdd.setElement(ONE_c,Indices(4,D_bond-1,1));
  FermiOdd.setElement(ONE_c,Indices(6,D_bond-1,1));
  
  FermiOdd.setElement(this->epsilon*ONE_c,Indices(0,2,2));
  FermiOdd.setElement(ONE_c,Indices(3,D_bond-1,2));
  FermiOdd.setElement(ONE_c,Indices(8,D_bond-1,2));
  
  FermiOdd.setElement(this->epsilon*ONE_c,Indices(0,3,3));
  FermiOdd.setElement(ONE_c,Indices(2,D_bond-1,3));
  FermiOdd.setElement(ONE_c,Indices(5,D_bond-1,3));
  
  FermiOdd.setElement(this->epsilon*ONE_c,Indices(0,4,4));
  FermiOdd.setElement(ONE_c,Indices(1,D_bond-1,4));
  FermiOdd.setElement(ONE_c,Indices(7,D_bond-1,4));
  
  FermiOdd.setElement(-this->mu*ONE_c,Indices(0,D_bond-1,5));
  
  //FermiEven
  FermiEven = FermiOdd;
  FermiEven.setElement(this->mu*ONE_c,Indices(0,D_bond-1,5));
  
  //Fermi1
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(this->epsilon*ONE_c,Indices(0,1,1));
  FermiFirst.setElement(this->epsilon*ONE_c,Indices(0,2,2));
  FermiFirst.setElement(this->epsilon*ONE_c,Indices(0,3,3));
  FermiFirst.setElement(this->epsilon*ONE_c,Indices(0,4,4));
  FermiFirst.setElement(-1.0*this->mu*ONE_c,Indices(0,D_bond-1,5));
  
  //FermiLast
  FermiLast.setElement(ONE_c,Indices(D_bond-1,0,0));
  
  FermiLast.setElement(ONE_c,Indices(4,0,1));
  FermiLast.setElement(ONE_c,Indices(6,0,1));
  
  FermiLast.setElement(ONE_c,Indices(3,0,2));
  FermiLast.setElement(ONE_c,Indices(8,0,2));
  
  FermiLast.setElement(ONE_c,Indices(2,0,3));
  FermiLast.setElement(ONE_c,Indices(5,0,3));
  
  FermiLast.setElement(ONE_c,Indices(1,0,4));
  FermiLast.setElement(ONE_c,Indices(7,0,4));
  
  FermiLast.setElement(pow(-1.0,((double) this->N_fermions))*this->mu*ONE_c,Indices(0,0,5)); 
  
  //Gauge
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(D_bond-1,D_bond-1,0));
  
  Gauge.setElement(ONE_c,Indices(2,2,1));
  Gauge.setElement(ONE_c,Indices(3,3,1));
  
  Gauge.setElement(ONE_c,Indices(2,2,2));
  Gauge.setElement(ONE_c,Indices(3,3,2));
  
  Gauge.setElement(ONE_c,Indices(1,1,3));
  Gauge.setElement(ONE_c,Indices(4,4,3));
  
  Gauge.setElement(ONE_c,Indices(1,1,4));
  Gauge.setElement(ONE_c,Indices(4,4,4));
  
  Gauge.setElement(ONE_c,Indices(3,6,5));
  Gauge.setElement(-1.0*ONE_c,Indices(2,7,5));
  
  Gauge.setElement(ONE_c,Indices(3,6,6));
  Gauge.setElement(-1.0*ONE_c,Indices(2,7,6));
  
  Gauge.setElement(ONE_c,Indices(1,5,7));
  Gauge.setElement(-1.0*ONE_c,Indices(4,8,7));
  
  Gauge.setElement(ONE_c,Indices(1,5,8));
  Gauge.setElement(-1.0*ONE_c,Indices(4,8,8));
  
  Gauge.setElement(0.5*this->g*this->g*ONE_c,Indices(0,D_bond-1,9));
  
  //Reshape and contract
  FermiFirst.reshape(Indices(1*D_bond,6));	FermiFirst.multiplyRight(Fermionic_operators);	FermiFirst.reshape(Indices(1,D_bond,this->d_fermi,this->d_fermi));	FermiFirst.permute(Indices(3,1,4,2));
  FermiOdd.reshape(Indices(D_bond*D_bond,6));	FermiOdd.multiplyRight(Fermionic_operators);	FermiOdd.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));	FermiOdd.permute(Indices(3,1,4,2));
  FermiEven.reshape(Indices(D_bond*D_bond,6));	FermiEven.multiplyRight(Fermionic_operators);	FermiEven.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));	FermiEven.permute(Indices(3,1,4,2));
  FermiLast.reshape(Indices(D_bond*1,6));	FermiLast.multiplyRight(Fermionic_operators);	FermiLast.reshape(Indices(D_bond,1,this->d_fermi,this->d_fermi));	FermiLast.permute(Indices(3,1,4,2));
  Gauge.reshape(Indices(D_bond*D_bond,10));	Gauge.multiplyRight(Gauge_operators);		Gauge.reshape(Indices(D_bond,D_bond,this->d_link,this->d_link));	Gauge.permute(Indices(3,1,4,2));
  
  //First and last entry
  hamil.setOp(0,new Operator(FermiFirst),true);
  hamil.setOp(this->N_total-1,new Operator(FermiLast),true);
  
  int pos_odd_saved,pos_even_saved,pos_gauge_saved;
  bool odd_saved=false,even_saved=false,gauge_saved=false;
  
  for(int i=1; i<(this->N_total-1); i+=2)
  {
    if(gauge_saved)
      hamil.setOp(i,&hamil.getOp(pos_gauge_saved),false);
    else
    {
      hamil.setOp(i,new Operator(Gauge),true);
      pos_gauge_saved=i;
      gauge_saved=true;
    }
  }
  
  for(int i=2; i<(this->N_total-1); i+=4)
  {
    if(even_saved)
      hamil.setOp(i,&hamil.getOp(pos_even_saved),false);
    else
    {
      hamil.setOp(i,new Operator(FermiEven),true);
      pos_even_saved=i;
      even_saved=true;
    }
  }
  
  for(int i=4; i<(this->N_total-1); i+=4)
  {
    if(odd_saved)
      hamil.setOp(i,&hamil.getOp(pos_odd_saved),false);
    else
    {
      hamil.setOp(i,new Operator(FermiOdd),true);
      pos_odd_saved=i;
      odd_saved=true;
    }
  }
  
  //cout << "Hamiltonian " << hamil << endl; 
#ifdef MATOUT
  hamil.exportForMatlab("ReducedSU2Hamiltonian.m");
#endif
}

//Build the Hamiltonian MPO
void ReducedSU2Hamiltonian::initHamiltonianPenalty()
{
  cout << "Building Hamiltonian MPO with" << endl
       << "epsilon: " << this->epsilon << endl
       << "mu:      " << this->mu << endl
       << "g:       " << this->g << endl
       << "dlink:   " << this->d_link << endl
       << "lambda:  " << this->lambda << endl;
       
       
  mwArray J,Jsq,Jcube,Jquarter,Qsquare,Qquarter;
  J=mwArray(Indices(this->d_link,this->d_link));		J.fillWithZero();
  Qsquare=mwArray(Indices(this->d_fermi,this->d_fermi));	Qsquare.fillWithZero();
  
  Qsquare.setElement(0.25*ONE_c,Indices(1,1));  
  Qquarter=Qsquare*Qsquare;
  for(int i=0; i<this->d_link; i++)
    J.setElement(((double)i)/2.0*ONE_c,Indices(i,i));  
  Jsq=J*J;
  Jcube=Jsq*J;
  Jquarter=Jcube*J;
  
  //Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(14,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //Lp/Lm site dependent
      Gauge_operators.setElement(this->Lmjdep.getElement(Indices(i,j)),Indices(1,i,j));
      Gauge_operators.setElement(this->Lpjdep.getElement(Indices(i,j)),Indices(2,i,j));
      Gauge_operators.setElement(this->Lmjdepdagger.getElement(Indices(i,j)),Indices(3,i,j));
      Gauge_operators.setElement(this->Lpjdepdagger.getElement(Indices(i,j)),Indices(4,i,j));
      //Lp/Lm site independent
      Gauge_operators.setElement(this->Lm.getElement(Indices(i,j)),Indices(5,i,j));
      Gauge_operators.setElement(this->Lp.getElement(Indices(i,j)),Indices(6,i,j));
      Gauge_operators.setElement(this->Lmdagger.getElement(Indices(i,j)),Indices(7,i,j));
      Gauge_operators.setElement(this->Lpdagger.getElement(Indices(i,j)),Indices(8,i,j));
      //Jsquared
      Gauge_operators.setElement(this->Jsquared.getElement(Indices(i,j)),Indices(9,i,j));
      //J,J^2,J^3,J^4
      Gauge_operators.setElement(J.getElement(Indices(i,j)),Indices(10,i,j));
      Gauge_operators.setElement(Jsq.getElement(Indices(i,j)),Indices(11,i,j));
      Gauge_operators.setElement(Jcube.getElement(Indices(i,j)),Indices(12,i,j));
      Gauge_operators.setElement(Jquarter.getElement(Indices(i,j)),Indices(13,i,j));
     }
  }   
  //cout << "Gauge_operators" << Gauge_operators <<endl;
  Gauge_operators.reshape(Indices(14,this->d_link*this->d_link));  
  //cout << "Bosonic operators done" << endl;
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(8,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();    
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //O1/O2
      Fermionic_operators.setElement(this->O1.getElement(Indices(i,j)),Indices(1,i,j));
      Fermionic_operators.setElement(this->O2.getElement(Indices(i,j)),Indices(2,i,j));
      Fermionic_operators.setElement(this->O1dagger.getElement(Indices(i,j)),Indices(3,i,j));
      Fermionic_operators.setElement(this->O2dagger.getElement(Indices(i,j)),Indices(4,i,j));
      //Occupation number operator
      Fermionic_operators.setElement(this->Noc.getElement(Indices(i,j)),Indices(5,i,j));
      //Q² and Q^4
      Fermionic_operators.setElement(Qsquare.getElement(Indices(i,j)),Indices(6,i,j));
      Fermionic_operators.setElement(Qquarter.getElement(Indices(i,j)),Indices(7,i,j));            
      
    }
  }  
  Fermionic_operators.reshape(Indices(8,this->d_fermi*this->d_fermi));
  //cout << "Fermionic operators done" << endl;
  
  //Now prepare the matrices
  int D_bond=15;
  mwArray FermiFirst,FermiEven,FermiOdd,FermiLast,Gauge;
  FermiFirst=mwArray(Indices(1,D_bond,8));	FermiFirst.fillWithZero();
  FermiEven=mwArray(Indices(D_bond,D_bond,8));	FermiEven.fillWithZero();
  FermiOdd=mwArray(Indices(D_bond,D_bond,8));	FermiOdd.fillWithZero();
  FermiLast=mwArray(Indices(D_bond,1,8));	FermiLast.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,14));	Gauge.fillWithZero();
  
  //FermiOdd
  FermiOdd.setElement(ONE_c,Indices(0,0,0));
  FermiOdd.setElement(ONE_c,Indices(D_bond-1,D_bond-1,0));
  FermiOdd.setElement(ONE_c,Indices(9,9,0));  
  FermiOdd.setElement(ONE_c,Indices(11,11,0));  
  FermiOdd.setElement(ONE_c,Indices(12,12,0));  
  
  FermiOdd.setElement(this->epsilon*ONE_c,Indices(0,1,1));
  FermiOdd.setElement(ONE_c,Indices(4,D_bond-1,1));
  FermiOdd.setElement(ONE_c,Indices(6,D_bond-1,1));
  
  FermiOdd.setElement(this->epsilon*ONE_c,Indices(0,2,2));
  FermiOdd.setElement(ONE_c,Indices(3,D_bond-1,2));
  FermiOdd.setElement(ONE_c,Indices(8,D_bond-1,2));
  
  FermiOdd.setElement(this->epsilon*ONE_c,Indices(0,3,3));
  FermiOdd.setElement(ONE_c,Indices(2,D_bond-1,3));
  FermiOdd.setElement(ONE_c,Indices(5,D_bond-1,3));
  
  FermiOdd.setElement(this->epsilon*ONE_c,Indices(0,4,4));
  FermiOdd.setElement(ONE_c,Indices(1,D_bond-1,4));
  FermiOdd.setElement(ONE_c,Indices(7,D_bond-1,4));
  
  FermiOdd.setElement(-this->mu*ONE_c,Indices(0,D_bond-1,5));
  
  FermiOdd.setElement(ONE_c,Indices(9,10,6));
  FermiOdd.setElement(-2.0*ONE_c,Indices(11,D_bond-1,6));
  FermiOdd.setElement(this->lambda*ONE_c,Indices(0,13,6));
  
  FermiOdd.setElement(this->lambda*ONE_c,Indices(0,D_bond-1,7));
  
  //FermiEven
  FermiEven = FermiOdd;
  FermiEven.setElement(this->mu*ONE_c,Indices(0,D_bond-1,5));
  
  //Fermi1
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(this->epsilon*ONE_c,Indices(0,1,1));
  FermiFirst.setElement(this->epsilon*ONE_c,Indices(0,2,2));
  FermiFirst.setElement(this->epsilon*ONE_c,Indices(0,3,3));
  FermiFirst.setElement(this->epsilon*ONE_c,Indices(0,4,4));
  FermiFirst.setElement(-1.0*this->mu*ONE_c,Indices(0,D_bond-1,5));
  FermiFirst.setElement(this->lambda*ONE_c,Indices(0,13,6));
  FermiFirst.setElement(this->lambda*ONE_c,Indices(0,D_bond-1,7));
  
  //FermiLast
  FermiLast.setElement(ONE_c,Indices(D_bond-1,0,0));
  
  FermiLast.setElement(ONE_c,Indices(4,0,1));
  FermiLast.setElement(ONE_c,Indices(6,0,1));
  
  FermiLast.setElement(ONE_c,Indices(3,0,2));
  FermiLast.setElement(ONE_c,Indices(8,0,2));
  
  FermiLast.setElement(ONE_c,Indices(2,0,3));
  FermiLast.setElement(ONE_c,Indices(5,0,3));
  
  FermiLast.setElement(ONE_c,Indices(1,0,4));
  FermiLast.setElement(ONE_c,Indices(7,0,4));
  
  FermiLast.setElement(pow(-1.0,((double) this->N_fermions))*this->mu*ONE_c,Indices(0,0,5)); 
  
  FermiLast.setElement(-2.0*ONE_c,Indices(11,0,6));
  
  FermiLast.setElement(this->lambda*ONE_c,Indices(0,0,7));
  
  //Gauge
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(D_bond-1,D_bond-1,0));
  
  Gauge.setElement(ONE_c,Indices(2,2,1));
  Gauge.setElement(ONE_c,Indices(3,3,1));
  
  Gauge.setElement(ONE_c,Indices(2,2,2));
  Gauge.setElement(ONE_c,Indices(3,3,2));
  
  Gauge.setElement(ONE_c,Indices(1,1,3));
  Gauge.setElement(ONE_c,Indices(4,4,3));
  
  Gauge.setElement(ONE_c,Indices(1,1,4));
  Gauge.setElement(ONE_c,Indices(4,4,4));
  
  Gauge.setElement(ONE_c,Indices(3,6,5));
  Gauge.setElement(-1.0*ONE_c,Indices(2,7,5));
  
  Gauge.setElement(ONE_c,Indices(3,6,6));
  Gauge.setElement(-1.0*ONE_c,Indices(2,7,6));
  
  Gauge.setElement(ONE_c,Indices(1,5,7));
  Gauge.setElement(-1.0*ONE_c,Indices(4,8,7));
  
  Gauge.setElement(ONE_c,Indices(1,5,8));
  Gauge.setElement(-1.0*ONE_c,Indices(4,8,8));
  
  Gauge.setElement(0.5*this->g*this->g*ONE_c,Indices(0,D_bond-1,9));
  
  Gauge.setElement(this->lambda*ONE_c,Indices(0,9,10));
  Gauge.setElement(4.0*ONE_c,Indices(10,D_bond-1,10));
  Gauge.setElement(-4.0*ONE_c,Indices(12,D_bond-1,10));
  
  Gauge.setElement(this->lambda*ONE_c,Indices(0,11,11));
  Gauge.setElement(6.0*ONE_c,Indices(11,D_bond-1,11));
  Gauge.setElement(-2.0*ONE_c,Indices(13,D_bond-1,11));
  
  Gauge.setElement(-4.0*ONE_c,Indices(9,D_bond-1,12));
  Gauge.setElement(this->lambda*ONE_c,Indices(0,12,12));
  
  Gauge.setElement(this->lambda*2.0*ONE_c,Indices(0,D_bond-1,13));
  
  //Reshape and contract
  FermiFirst.reshape(Indices(1*D_bond,8));	FermiFirst.multiplyRight(Fermionic_operators);	FermiFirst.reshape(Indices(1,D_bond,this->d_fermi,this->d_fermi));	FermiFirst.permute(Indices(3,1,4,2));
  FermiOdd.reshape(Indices(D_bond*D_bond,8));	FermiOdd.multiplyRight(Fermionic_operators);	FermiOdd.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));	FermiOdd.permute(Indices(3,1,4,2));
  FermiEven.reshape(Indices(D_bond*D_bond,8));	FermiEven.multiplyRight(Fermionic_operators);	FermiEven.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));	FermiEven.permute(Indices(3,1,4,2));
  FermiLast.reshape(Indices(D_bond*1,8));	FermiLast.multiplyRight(Fermionic_operators);	FermiLast.reshape(Indices(D_bond,1,this->d_fermi,this->d_fermi));	FermiLast.permute(Indices(3,1,4,2));
  Gauge.reshape(Indices(D_bond*D_bond,14));	Gauge.multiplyRight(Gauge_operators);		Gauge.reshape(Indices(D_bond,D_bond,this->d_link,this->d_link));	Gauge.permute(Indices(3,1,4,2));
  
  //First and last entry
  hamil.setOp(0,new Operator(FermiFirst),true);
  hamil.setOp(this->N_total-1,new Operator(FermiLast),true);
  
  int pos_odd_saved,pos_even_saved,pos_gauge_saved;
  bool odd_saved=false,even_saved=false,gauge_saved=false;
  
  for(int i=1; i<(this->N_total-1); i+=2)
  {
    if(gauge_saved)
      hamil.setOp(i,&hamil.getOp(pos_gauge_saved),false);
    else
    {
      hamil.setOp(i,new Operator(Gauge),true);
      pos_gauge_saved=i;
      gauge_saved=true;
    }
  }
  
  for(int i=2; i<(this->N_total-1); i+=4)
  {
    if(even_saved)
      hamil.setOp(i,&hamil.getOp(pos_even_saved),false);
    else
    {
      hamil.setOp(i,new Operator(FermiEven),true);
      pos_even_saved=i;
      even_saved=true;
    }
  }
  
  for(int i=4; i<(this->N_total-1); i+=4)
  {
    if(odd_saved)
      hamil.setOp(i,&hamil.getOp(pos_odd_saved),false);
    else
    {
      hamil.setOp(i,new Operator(FermiOdd),true);
      pos_odd_saved=i;
      odd_saved=true;
    }
  }
  
  //cout << "Hamiltonian " << hamil << endl; 
#ifdef MATOUT
  hamil.exportForMatlab("ReducedSU2HamiltonianPenalty.m");
#endif
}

void ReducedSU2Hamiltonian::removeGaugeInteraction()
{
  cout << "WARNING, will remove the gauge interaction, the Hamiltonian MPO has to be updated accordingly!" << endl;
  this->Lpjdep=identityMatrix(this->d_link);
  this->Lmjdep=identityMatrix(this->d_link);
  this->Lpjdepdagger=identityMatrix(this->d_link);
  this->Lmjdepdagger=identityMatrix(this->d_link);
  this->Lm=identityMatrix(this->d_link);
  this->Lp =identityMatrix(this->d_link);
  this->Lmdagger=identityMatrix(this->d_link);
  this->Lpdagger=identityMatrix(this->d_link);
}

  
void ReducedSU2Hamiltonian::getJsquaredMPO(MPO &Spin, int index) const
{
  //Some error checking
  if(index<1 || index>(this->N_fermions-1))
  {
    cout << "Error in ReducedSU2Hamiltonian::getJsquaredMPO, index is invalid" << endl;
    exit(666);
  }
  //Call the more general routine with appropriately transformed index
  getSingleBodyMPO(Spin,2*index-1,JsqOp);
}

void ReducedSU2Hamiltonian::getLMPO(MPO &mpo, int index) const
{
  //Some error checking
  if(index<1 || index>(this->N_fermions-1))
  {
    cout << "Error in ReducedSU2Hamiltonian::getLMPO, index is invalid" << endl;
    exit(666);
  }
  //Call the more general routine with appropriately transformed index
  getSingleBodyMPO(mpo,2*index-1,LOp);
}

void ReducedSU2Hamiltonian::getNocMPO(MPO &Spin, int index) const
{
  //Some error checking
  if(index<1 || index>this->N_fermions)
  {
    cout << "Error in ReducedSU2Hamiltonian::getJsquaredMPO, index is invalid" << endl;
    exit(666);
  }
  //Call the more general routine with appropriately transformed index
  getSingleBodyMPO(Spin,2*index-2,NfOp);
}


void ReducedSU2Hamiltonian::getSingleBodyMPO(MPO &Spin, int index, SingleBodyOperator Op) const
{
  //Frist clear the MPO and resize it
  Spin.clear();
  Spin.initLength(this->N_total);
  
  cout.flush();
  
  //Prepare the matrices
  mwArray Link_id=this->id_link_mpo;
  mwArray Link_op;
  mwArray Fermi_id=this->id_fermi_mpo;
  mwArray Fermi_op;
  
  switch(Op)
  {
    case NfOp:
      Fermi_op=this->Noc; Fermi_op.reshape(Indices(this->d_fermi,1,this->d_fermi,1)); break;
    case JsqOp:
      Link_op=this->Jsquared; Link_op.reshape(Indices(this->d_link,1,this->d_link,1)); break;
    case LOp:
      Link_op=this->Lvalop; Link_op.reshape(Indices(this->d_link,1,this->d_link,1)); break;
    default:
      cout << "Unknown operator" << endl;
      exit(666);
  }
  
  //Variable for the Operators
  Operator* Local_op; 
  
  //Now set operators for fermionic sites
  bool fermi_saved=false;
  int saved_site;
  for(int k=0; k<this->N_total; k+=2)
  {
    if(fermi_saved && (k!=index))
      Spin.setOp(k,&Spin.getOp(saved_site),false);
    else
    {
      if(k==index)
      {
	Local_op=new Operator(Fermi_op);
	Spin.setOp(k,Local_op,true);
      }
      else
      {
	Local_op=new Operator(Fermi_id);
	Spin.setOp(k,Local_op,true);
	fermi_saved=true;
	saved_site=k;
      }
    }
  }
  
  //Now set the link operators
  bool link_saved=false;
  for(int k=1; k<this->N_total; k+=2)
  {
    if(link_saved && (k!=index))
      Spin.setOp(k,&Spin.getOp(saved_site),false);
    else
    {
      if(k==index)
      {
	Local_op=new Operator(Link_op);
	Spin.setOp(k,Local_op,true);
      }
      else
      {
	Local_op=new Operator(Link_id);
	Spin.setOp(k,Local_op,true);
	link_saved=true;
	saved_site=k;
      }
    }
  }
  //cout << "Spin: " << Spin << endl;
}

void ReducedSU2Hamiltonian::getUMPO(MPO &mpo, double dt, bool odd, bool imag_time) const
{
  mpo.initLength(this->N_total);
  complex_t delta_t=-I_c*dt;
  
  //In case i want imaginary time
  if(imag_time)
    delta_t=-I_c*delta_t;
  
  cout << "delta_t=" << delta_t<< endl;
  
  //Construct the necessary terms and sum them
  mwArray Local_op,exp_op,Local_op_first,Local_op_last;
  //Since I have two body terms and N_f flavors, the dimension of the final operator is
  
  Local_op = this->epsilon*kron(this->O1dagger,kron(this->Lmjdep,this->O2))
           + this->epsilon*kron(this->O1dagger,kron(this->Lpjdep,this->O2))
	   + this->epsilon*kron(this->O1,kron(this->Lmjdepdagger,this->O2dagger))
	   + this->epsilon*kron(this->O1,kron(this->Lpjdepdagger,this->O2dagger))
	   + this->epsilon*kron(this->O2dagger,kron(this->Lmjdepdagger,this->O1))
	   + this->epsilon*kron(this->O2dagger,kron(this->Lpjdepdagger,this->O1))
	   + this->epsilon*kron(this->O2,kron(this->Lmjdep,this->O1dagger))
	   + this->epsilon*kron(this->O2,kron(this->Lpjdep,this->O1dagger))
	   + this->epsilon*kron(this->O1,kron(this->Lmdagger,this->O1dagger))
	   + this->epsilon*kron(this->O1,kron(this->Lpdagger,this->O1dagger))
	   + this->epsilon*kron(this->O1dagger,kron(this->Lm,this->O1))
	   + this->epsilon*kron(this->O1dagger,kron(this->Lp,this->O1))
	   + this->epsilon*kron(this->O2,kron(-1.0*this->Lm,this->O2dagger))
	   + this->epsilon*kron(this->O2,kron(-1.0*this->Lp,this->O2dagger))
	   + this->epsilon*kron(this->O2dagger,kron(-1.0*this->Lmdagger,this->O2))
	   + this->epsilon*kron(this->O2dagger,kron(-1.0*this->Lpdagger,this->O2))
	   + 0.5*this->g*this->g*kron(this->id_fermi,kron(this->Jsquared,this->id_fermi));
  if(odd)
  {
    //Before adding the mass term, construct special first and last operators
    Local_op_first = Local_op + -1.0*this->mu*kron(this->Noc,kron(this->id_link,this->id_fermi)) + 0.5*this->mu*kron(this->id_fermi,kron(this->id_link,this->Noc));
    Local_op_last = Local_op + -0.5*this->mu*kron(this->Noc,kron(this->id_link,this->id_fermi)) + 1.0*this->mu*kron(this->id_fermi,kron(this->id_link,this->Noc));
    //Now add the mass term for sites in between
    Local_op  = Local_op + -0.5*this->mu*kron(this->Noc,kron(this->id_link,this->id_fermi)) + 0.5*this->mu*kron(this->id_fermi,kron(this->id_link,this->Noc));
    
  }
  else
  {
    //Since I start counting with 1, the first site is always odd, if I look at the even operator, I do not have to care about the first site    
    Local_op_last  = Local_op + 0.5*this->mu*kron(this->Noc,kron(this->id_link,this->id_fermi)) + -1.0*this->mu*kron(this->id_fermi,kron(this->id_link,this->Noc));
    Local_op  = Local_op + 0.5*this->mu*kron(this->Noc,kron(this->id_link,this->id_fermi)) + -0.5*this->mu*kron(this->id_fermi,kron(this->id_link,this->Noc));
  }
  
  wrapper::expm(Local_op,exp_op,delta_t);
  
  //Now split into single matrices
  vector<int> phys_dims;
  phys_dims.push_back(this->d_fermi);	phys_dims.push_back(this->d_link);	phys_dims.push_back(this->d_fermi);
  vector<mwArray> matrices,matrices_last,matrices_first;
  splitLocalOp(exp_op,phys_dims,matrices);
  
  //Now put the matrices
  int start_site=odd?1:2;
  int site_index;
  
  for(int i=start_site; i<this->N_fermions; i+=2)
  {
    if(i!=1 && i!=(this->N_fermions-1))
    {
      mpo.setOp(2*i-2,new Operator(matrices[0]),true);
      mpo.setOp(2*i-1,new Operator(matrices[1]),true);
      mpo.setOp(2*(i+1)-2,new Operator(matrices[2]),true);
    }
    else if(i==1)
    {
      cout << "Setting special first entry" << endl;
      //Compute decomposition on the fly as it might not be needed and therefore pre-computation would be a waste of time
      wrapper::expm(Local_op_first,exp_op,delta_t);
      splitLocalOp(exp_op,phys_dims,matrices_first);
      mpo.setOp(2*i-2,new Operator(matrices_first[0]),true);
      mpo.setOp(2*i-1,new Operator(matrices_first[1]),true);
      mpo.setOp(2*(i+1)-2,new Operator(matrices_first[2]),true);
    }
    else if(i==(this->N_fermions-1))
    {
      cout << "Setting special last entry" << endl;
      //Same here
      wrapper::expm(Local_op_last,exp_op,delta_t);
      splitLocalOp(exp_op,phys_dims,matrices_last);
      mpo.setOp(2*i-2,new Operator(matrices_last[0]),true);
      mpo.setOp(2*i-1,new Operator(matrices_last[1]),true);
      mpo.setOp(2*(i+1)-2,new Operator(matrices_last[2]),true);
    }
  }  
  
  //Remaining identities
  for(int i=0; i<this->N_total; i++)
  {
    if(mpo.isEmpty(i) && ((i%2)==0))
      mpo.setOp(i,new Operator(this->id_fermi_mpo),true);
    else if (mpo.isEmpty(i) && ((i%2)==1))
      mpo.setOp(i,new Operator(this->id_link_mpo),true);
  }
#ifdef MATOUT
  if(odd)
    mpo.exportForMatlab("Uodd.m");
  else
    mpo.exportForMatlab("Ueven.m");
#endif
}

void ReducedSU2Hamiltonian::getUoddMPO(MPO &mpo, double dt, bool imag_time) const
{
  this->getUMPO(mpo,dt,true,imag_time);
}
void ReducedSU2Hamiltonian::getUevenMPO(MPO &mpo, double dt, bool imag_time) const
{
  this->getUMPO(mpo,dt,false,imag_time);
}

void ReducedSU2Hamiltonian::getUMPO(MPO &mpo, double dt, bool imag_time) const
{
  mpo.initLength(this->N_total);
  complex_t delta_t=-I_c*dt;
  
  //In case i want imaginary time
  if(imag_time)
    delta_t=-I_c*delta_t;
  
  cout << "delta_t=" << delta_t<< endl;
  
  //Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(10,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //Lp/Lm site dependent
      Gauge_operators.setElement(this->Lmjdep.getElement(Indices(i,j)),Indices(1,i,j));
      Gauge_operators.setElement(this->Lpjdep.getElement(Indices(i,j)),Indices(2,i,j));
      Gauge_operators.setElement(this->Lmjdepdagger.getElement(Indices(i,j)),Indices(3,i,j));
      Gauge_operators.setElement(this->Lpjdepdagger.getElement(Indices(i,j)),Indices(4,i,j));
      //Lp/Lm site independent
      Gauge_operators.setElement(this->Lm.getElement(Indices(i,j)),Indices(5,i,j));
      Gauge_operators.setElement(this->Lp.getElement(Indices(i,j)),Indices(6,i,j));
      Gauge_operators.setElement(this->Lmdagger.getElement(Indices(i,j)),Indices(7,i,j));
      Gauge_operators.setElement(this->Lpdagger.getElement(Indices(i,j)),Indices(8,i,j));
      //Jsquared
      Gauge_operators.setElement(this->Jsquared.getElement(Indices(i,j)),Indices(9,i,j));
     }
  }   
  //cout << "Gauge_operators" << Gauge_operators <<endl;
  Gauge_operators.reshape(Indices(10,this->d_link*this->d_link));  
  //cout << "Bosonic operators done" << endl;
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(6,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();    
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //O1/O2
      Fermionic_operators.setElement(this->O1.getElement(Indices(i,j)),Indices(1,i,j));
      Fermionic_operators.setElement(this->O2.getElement(Indices(i,j)),Indices(2,i,j));
      Fermionic_operators.setElement(this->O1dagger.getElement(Indices(i,j)),Indices(3,i,j));
      Fermionic_operators.setElement(this->O2dagger.getElement(Indices(i,j)),Indices(4,i,j));
      //Occupation number operator
      Fermionic_operators.setElement(this->Noc.getElement(Indices(i,j)),Indices(5,i,j));
      
    }
  }  
  Fermionic_operators.reshape(Indices(6,this->d_fermi*this->d_fermi));
  //cout << "Fermionic operators done" << endl;
  
  //Now prepare the matrices
  int D_bond=10;
  mwArray FermiFirst,FermiEven,FermiOdd,FermiLast,Gauge;
  FermiFirst=mwArray(Indices(1,D_bond,6));	FermiFirst.fillWithZero();
  FermiEven=mwArray(Indices(D_bond,D_bond,6));	FermiEven.fillWithZero();
  FermiOdd=mwArray(Indices(D_bond,D_bond,6));	FermiOdd.fillWithZero();
  FermiLast=mwArray(Indices(D_bond,1,6));	FermiLast.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,10));	Gauge.fillWithZero();
  
  //FermiOdd
  FermiOdd.setElement(ONE_c,Indices(0,0,0));
  FermiOdd.setElement(ONE_c,Indices(D_bond-1,D_bond-1,0));
  FermiOdd.setElement(ONE_c/this->N_fermions,Indices(0,D_bond-1,0));
  
  FermiOdd.setElement(delta_t*this->epsilon*ONE_c,Indices(0,1,1));
  FermiOdd.setElement(ONE_c,Indices(4,D_bond-1,1));
  FermiOdd.setElement(ONE_c,Indices(6,D_bond-1,1));
  
  FermiOdd.setElement(delta_t*this->epsilon*ONE_c,Indices(0,2,2));
  FermiOdd.setElement(ONE_c,Indices(3,D_bond-1,2));
  FermiOdd.setElement(ONE_c,Indices(8,D_bond-1,2));
  
  FermiOdd.setElement(delta_t*this->epsilon*ONE_c,Indices(0,3,3));
  FermiOdd.setElement(ONE_c,Indices(2,D_bond-1,3));
  FermiOdd.setElement(ONE_c,Indices(5,D_bond-1,3));
  
  FermiOdd.setElement(delta_t*this->epsilon*ONE_c,Indices(0,4,4));
  FermiOdd.setElement(ONE_c,Indices(1,D_bond-1,4));
  FermiOdd.setElement(ONE_c,Indices(7,D_bond-1,4));
  
  FermiOdd.setElement(delta_t*(-this->mu*ONE_c),Indices(0,D_bond-1,5));
  
  //FermiEven
  FermiEven = FermiOdd;
  FermiEven.setElement(delta_t*this->mu*ONE_c,Indices(0,D_bond-1,5));
  
  //Fermi1
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(ONE_c/this->N_fermions,Indices(0,D_bond-1,0));
  FermiFirst.setElement(delta_t*this->epsilon*ONE_c,Indices(0,1,1));
  FermiFirst.setElement(delta_t*this->epsilon*ONE_c,Indices(0,2,2));
  FermiFirst.setElement(delta_t*this->epsilon*ONE_c,Indices(0,3,3));
  FermiFirst.setElement(delta_t*this->epsilon*ONE_c,Indices(0,4,4));
  FermiFirst.setElement(delta_t*(-1.0*this->mu*ONE_c),Indices(0,D_bond-1,5));
  
  //FermiLast
  FermiLast.setElement(ONE_c,Indices(D_bond-1,0,0));
  FermiLast.setElement(ONE_c/this->N_fermions,Indices(0,0,0));
  
  FermiLast.setElement(ONE_c,Indices(4,0,1));
  FermiLast.setElement(ONE_c,Indices(6,0,1));
  
  FermiLast.setElement(ONE_c,Indices(3,0,2));
  FermiLast.setElement(ONE_c,Indices(8,0,2));
  
  FermiLast.setElement(ONE_c,Indices(2,0,3));
  FermiLast.setElement(ONE_c,Indices(5,0,3));
  
  FermiLast.setElement(ONE_c,Indices(1,0,4));
  FermiLast.setElement(ONE_c,Indices(7,0,4));
  
  FermiLast.setElement(delta_t*pow(-1.0,((double) this->N_fermions))*this->mu*ONE_c,Indices(0,0,5)); 
  
  //Gauge
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(D_bond-1,D_bond-1,0));
  
  Gauge.setElement(ONE_c,Indices(2,2,1));
  Gauge.setElement(ONE_c,Indices(3,3,1));
  
  Gauge.setElement(ONE_c,Indices(2,2,2));
  Gauge.setElement(ONE_c,Indices(3,3,2));
  
  Gauge.setElement(ONE_c,Indices(1,1,3));
  Gauge.setElement(ONE_c,Indices(4,4,3));
  
  Gauge.setElement(ONE_c,Indices(1,1,4));
  Gauge.setElement(ONE_c,Indices(4,4,4));
  
  Gauge.setElement(ONE_c,Indices(3,6,5));
  Gauge.setElement(-1.0*ONE_c,Indices(2,7,5));
  
  Gauge.setElement(ONE_c,Indices(3,6,6));
  Gauge.setElement(-1.0*ONE_c,Indices(2,7,6));
  
  Gauge.setElement(ONE_c,Indices(1,5,7));
  Gauge.setElement(-1.0*ONE_c,Indices(4,8,7));
  
  Gauge.setElement(ONE_c,Indices(1,5,8));
  Gauge.setElement(-1.0*ONE_c,Indices(4,8,8));
  
  Gauge.setElement(delta_t*0.5*this->g*this->g*ONE_c,Indices(0,D_bond-1,9));
  
  //Reshape and contract
  FermiFirst.reshape(Indices(1*D_bond,6));	FermiFirst.multiplyRight(Fermionic_operators);	FermiFirst.reshape(Indices(1,D_bond,this->d_fermi,this->d_fermi));	FermiFirst.permute(Indices(3,1,4,2));
  FermiOdd.reshape(Indices(D_bond*D_bond,6));	FermiOdd.multiplyRight(Fermionic_operators);	FermiOdd.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));	FermiOdd.permute(Indices(3,1,4,2));
  FermiEven.reshape(Indices(D_bond*D_bond,6));	FermiEven.multiplyRight(Fermionic_operators);	FermiEven.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));	FermiEven.permute(Indices(3,1,4,2));
  FermiLast.reshape(Indices(D_bond*1,6));	FermiLast.multiplyRight(Fermionic_operators);	FermiLast.reshape(Indices(D_bond,1,this->d_fermi,this->d_fermi));	FermiLast.permute(Indices(3,1,4,2));
  Gauge.reshape(Indices(D_bond*D_bond,10));	Gauge.multiplyRight(Gauge_operators);		Gauge.reshape(Indices(D_bond,D_bond,this->d_link,this->d_link));	Gauge.permute(Indices(3,1,4,2));
  
  //First and last entry
  mpo.setOp(0,new Operator(FermiFirst),true);
  mpo.setOp(this->N_total-1,new Operator(FermiLast),true);
  
  int pos_odd_saved,pos_even_saved,pos_gauge_saved;
  bool odd_saved=false,even_saved=false,gauge_saved=false;
  
  for(int i=1; i<(this->N_total-1); i+=2)
  {
    if(gauge_saved)
      mpo.setOp(i,&mpo.getOp(pos_gauge_saved),false);
    else
    {
      mpo.setOp(i,new Operator(Gauge),true);
      pos_gauge_saved=i;
      gauge_saved=true;
    }
  }
  
  for(int i=2; i<(this->N_total-1); i+=4)
  {
    if(even_saved)
      mpo.setOp(i,&mpo.getOp(pos_even_saved),false);
    else
    {
      mpo.setOp(i,new Operator(FermiEven),true);
      pos_even_saved=i;
      even_saved=true;
    }
  }
  
  for(int i=4; i<(this->N_total-1); i+=4)
  {
    if(odd_saved)
      mpo.setOp(i,&mpo.getOp(pos_odd_saved),false);
    else
    {
      mpo.setOp(i,new Operator(FermiOdd),true);
      pos_odd_saved=i;
      odd_saved=true;
    }
  }
  //cout << "UMPO " << mpo << endl;
#ifdef MATOUT
  mpo.exportForMatlab("Umpo.m");
#endif
}


void ReducedSU2Hamiltonian::constructInitialMPS(MPS &mps) const
{
  //Adjust the physical dimensions of the MPS, therefore overwrite the old one and create soemthing with matching dimensions
//   vector<int> phys_dims(this->N_total,this->d_fermi);  
//   for(int k=1; k<(this->N_total-1); k+=2)
//     phys_dims[k]=this->d_link;    
//   
//   mps.clear();
//   mps=MPS(this->N_total,1,phys_dims);
//   
//   
//   mwArray fermi_empty(Indices(this->d_fermi,1,1));
//   mwArray fermi_doubly_occupied(Indices(this->d_fermi,1,1));
//   mwArray link(Indices(this->d_link,1,1));
//   
//   fermi_empty.fillWithZero();
//   fermi_doubly_occupied.fillWithZero();
//   link.fillWithZero();
//   
//   fermi_empty.setElement(ONE_c,Indices(0,0,0));
//   fermi_doubly_occupied.setElement(ONE_c,Indices(2,0,0));
//   link.setElement(ONE_c,Indices(0,0,0));
//   
//   
//   //Set the matrices 
//   for(int k=1; k<=this->N_fermions; k++)
//   {    
//     if((k%2)!=0)
//       mps.setA(2*k-2,fermi_doubly_occupied);
//     else
//       mps.setA(2*k-2,fermi_empty);
//     if(k<this->N_fermions)
//       mps.setA(2*k-1,link);
//   }    
//   mps.gaugeCond('R',true);
  vector<int> sites(this->N_total,2);
  
  for(int i=2; i<=this->N_fermions; i+=2)
  {
    //Set fermionic site and previous link
    sites[2*(i-1)-1]=0;
    sites[2*i-2]=0;
    //If it is still inside the system, also take the link to the right
    if(i<this->N_fermions)
      sites[2*i-1]=0;
  }
  this->constructProductStateMPS(mps,sites); 
}

void ReducedSU2Hamiltonian::constructProductStateMPS(MPS &mps,vector<int> sites) const
{
  //Check if the amount of sites is ok
  if(sites.size()!= this->N_total)
  {
    cout << "Error in ReducedSU2Hamiltonian::constructProductStateMPS(), input vector has size " << sites.size() << " which is not compliant with the total number of sites " << this->N_total << endl;
    exit(666);
  }
  //Check if the fermionic occupation numbers and the link states are within a valid range
  for(int i=0; i<sites.size(); i++)
  {
    if((i%2)==0)
    {
      //Fermionic site, occupation number can be 0,1,2
      if(sites[i]<0 || sites[i]>2)
      {
	cout << "Error in ReducedSU2Hamiltonian::constructProductStateMPS(), input vector has invalid entries for fermionic sites" << endl;
	exit(666);
      }
    }
    else
    {
      //Link, state number can go from 0 .. d_link-1
      if(sites[i]<0 || sites[i]>= d_link)
      {
	cout << "Error in ReducedSU2Hamiltonian::constructProductStateMPS(), input vector has invalid entries for gauge links" << endl;
	exit(666);
      }
    }     
  } 
  //Adjust the physical dimensions of the MPS, therefore overwrite the old one and create soemthing with matching dimensions
  vector<int> phys_dims(this->N_total,this->d_fermi);  
  for(int k=1; k<(this->N_total-1); k+=2)
    phys_dims[k]=this->d_link;    
  
  mps.clear();
  mps=MPS(this->N_total,1,phys_dims);
  
  mwArray fermi(Indices(this->d_fermi,1,1));
  mwArray link(Indices(this->d_link,1,1));
  //Set the matrices 
  for(int k=0; k<this->N_total; k++)
  {    
    if((k%2)==0)
    {
      //Fermionic site
      fermi.fillWithZero();
      fermi.setElement(ONE_c,Indices(sites[k],0,0));
      mps.setA(k,fermi);
    }
    else
    {
      //Link
      link.fillWithZero();
      link.setElement(ONE_c,Indices(sites[k],0,0));
      mps.setA(k,link);
    }
  }    
  mps.gaugeCond('R',true); 
}

void ReducedSU2Hamiltonian::constructStringMPS(MPS &mps,vector<string_t> str_config,int quanta) const
{
  //First get a strong coupling state, the start replacing the tensors accordingly
  this->constructInitialMPS(mps);
  
  mwArray link(Indices(this->d_link,1,1));  
  link.fillWithZero();
  
  link.setElement(ONE_c,Indices(quanta,0,0));  
  
  //Set the matrices
  int curr_pos,curr_leng;
  for(int i=0; i<str_config.size(); i++)
  {
    curr_leng=str_config[i].leng;
    curr_pos=str_config[i].pos;
    for(int k=(2*curr_pos-1); k<=(2*(curr_pos+curr_leng-1)-1); k+=2)
      mps.setA(k,link);
  }
  mps.gaugeCond('R',true);
}

void ReducedSU2Hamiltonian::getHeavyString(MPS &mps, unsigned int length, unsigned int flux_quanta, unsigned int pos) const
{
  //Adjust the physical dimensions of the MPS, therefore overwrite the old one and create soemthing with matching dimensions
  vector<int> phys_dims(this->N_total,this->d_fermi);  
  for(int k=1; k<(this->N_total-1); k+=2)
    phys_dims[k]=this->d_link;    
  
  //Check if desired flux_quanta are possible with number of states allowed
  if(flux_quanta>=this->d_link || flux_quanta==0)
  {
    cout << "Error in ReducedSU2Hamiltonian::getHeavyString(), number of flux quanta is not compliant with size of link Hilbert space, will use largest flux quantum possible" << endl;
    flux_quanta=this->d_link-1;
  }
  
  mps.clear();
  mps=MPS(this->N_total,1,phys_dims);  
  
  mwArray fermi_empty(Indices(this->d_fermi,1,1));
  mwArray fermi_doubly_occupied(Indices(this->d_fermi,1,1));
  mwArray link(Indices(this->d_link,1,1));
  mwArray link_flux(Indices(this->d_link,1,1));
  
  fermi_empty.fillWithZero();
  fermi_doubly_occupied.fillWithZero();
  link.fillWithZero();
  link_flux.fillWithZero();
  
  fermi_empty.setElement(ONE_c,Indices(0,0,0));
  fermi_doubly_occupied.setElement(ONE_c,Indices(2,0,0));
  link.setElement(ONE_c,Indices(0,0,0));  
  link_flux.setElement(ONE_c,Indices(flux_quanta,0,0));
  
  //Determine the starting position of the string
  if(pos<1)
    pos = round(this->N_fermions/2)-round(length/2);
  //Some error checking
  if(pos>this->N_total || pos<1 || (pos+length)>this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getHeavyString, desired string of length " << length<< " starting on position " << pos << " is outside the the chain of " << this->N_fermions << " spins" << endl;
    exit(666);
  }
  
  cout << "Creating strong coupling string between heavy external charges with " << endl
       << "pos:    " << pos << endl
       << "length: " << length << endl  
       << "quanta: " << flux_quanta << endl;
  
  
  //Set the matrices
  int index; 
  for(int k=1; k<=this->N_fermions; k++)
  {    
    if((k%2)!=0)
      mps.setA(2*k-2,fermi_doubly_occupied);
    else
      mps.setA(2*k-2,fermi_empty);
    if(k<this->N_fermions && k>=pos && k<(pos+length))
      mps.setA(2*k-1,link_flux);
    else if(k<this->N_fermions)
      mps.setA(2*k-1,link);
  }    
  mps.gaugeCond('R',true);
}

void ReducedSU2Hamiltonian::getChargeSquareMPO(MPO &mpo, int pos) const
{
  mpo.initLength(this->N_total);
  
  //Check if given position is valid
  if(pos<1 || pos>this->N_fermions)
  {
    cout << "Error in ReducedSU2Hamiltonian::getChargeSquareMPO(), given position is not valid" << endl;
    exit(666);
  }
  
  mwArray Local_op(Indices(this->d_fermi,this->d_fermi));	Local_op.fillWithZero();	Local_op.setElement(0.25*ONE_c,Indices(1,1));
  Local_op.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
  
  bool id_fermi_saved=false,id_link_saved=false;
  int pos_id_fermi_saved,pos_id_link_saved;
  
  for(int i=1; i<=this->N_fermions; i++)
  {
    if(i==pos)
    {
      mpo.setOp(2*i-2,new Operator(Local_op),true);
      if(id_link_saved && i<this->N_fermions)
	mpo.setOp(2*i-1,&mpo.getOp(pos_id_link_saved),false);
      else if(i<this->N_fermions)
      {
	mpo.setOp(2*i-1,new Operator(this->id_link_mpo),true);
	id_link_saved=true;
	pos_id_link_saved=2*i-1;
      }
    }
    else
    {
      if(id_fermi_saved)
	mpo.setOp(2*i-2,&mpo.getOp(pos_id_fermi_saved),false);
      else
      {
	mpo.setOp(2*i-2,new Operator(this->id_fermi_mpo),true);
	id_fermi_saved=true;
	pos_id_fermi_saved=2*i-2;
      }
      if(id_link_saved && i<this->N_fermions)
	mpo.setOp(2*i-1,&mpo.getOp(pos_id_link_saved),false);
      else if(i<this->N_fermions)
      {
	mpo.setOp(2*i-1,new Operator(this->id_link_mpo),true);
	id_link_saved=true;
	pos_id_link_saved=2*i-1;
      }
    }
  }
}

void ReducedSU2Hamiltonian::getCondensateMPO(MPO &mpo, bool factor_off) const
{
  mpo.initLength(this->N_total);  
  
  double factor = sqrt(this->epsilon)/this->N_fermions;
  if(factor_off)
    factor=1.0;
  
  int D_bond = 2; 
  
  cout << "Constructing condensate MPO" << endl;
  
  //Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(1,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
     }
  }   
  //cout << "Gauge_operators" << Gauge_operators <<endl;
  Gauge_operators.reshape(Indices(1,this->d_link*this->d_link));  
  //cout << "Bosonic operators done" << endl;
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();    
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));      
      //Occupation number operator
      Fermionic_operators.setElement(this->Noc.getElement(Indices(i,j)),Indices(1,i,j));
      
    }
  }  
  Fermionic_operators.reshape(Indices(2,this->d_fermi*this->d_fermi));
  
  mwArray FermiFirst, FermiLast, FermiOdd, FermiEven, Gauge;
  FermiFirst=mwArray(Indices(1,D_bond,2));	FermiFirst.fillWithZero();
  FermiOdd=mwArray(Indices(D_bond,D_bond,2));	FermiOdd.fillWithZero();
  FermiEven=mwArray(Indices(D_bond,D_bond,2));	FermiEven.fillWithZero();
  FermiLast=mwArray(Indices(D_bond,1,2));		FermiLast.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,1));	Gauge.fillWithZero();
  
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(-1.0*factor*ONE_c,Indices(0,D_bond-1,1));
  
  FermiOdd.setElement(ONE_c,Indices(0,0,0));
  FermiOdd.setElement(ONE_c,Indices(D_bond-1,D_bond-1,0));
  FermiOdd.setElement(-1.0*factor*ONE_c,Indices(0,D_bond-1,1));
  
  FermiEven.setElement(ONE_c,Indices(0,0,0));
  FermiEven.setElement(ONE_c,Indices(D_bond-1,D_bond-1,0));
  FermiEven.setElement(factor*ONE_c,Indices(0,D_bond-1,1));
  
  FermiLast.setElement(ONE_c,Indices(D_bond-1,0,0));
  FermiLast.setElement(pow(-1.0,this->N_fermions)*factor*ONE_c,Indices(0,0,1));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(D_bond-1,D_bond-1,0));
  
  //Reshape and contract
  FermiFirst.reshape(Indices(1*D_bond,2));	FermiFirst.multiplyRight(Fermionic_operators);	FermiFirst.reshape(Indices(1,D_bond,this->d_fermi,this->d_fermi));	FermiFirst.permute(Indices(3,1,4,2));
  FermiOdd.reshape(Indices(D_bond*D_bond,2));	FermiOdd.multiplyRight(Fermionic_operators);	FermiOdd.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));	FermiOdd.permute(Indices(3,1,4,2));
  FermiEven.reshape(Indices(D_bond*D_bond,2));	FermiEven.multiplyRight(Fermionic_operators);	FermiEven.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));	FermiEven.permute(Indices(3,1,4,2));
  FermiLast.reshape(Indices(D_bond*1,2));	FermiLast.multiplyRight(Fermionic_operators);	FermiLast.reshape(Indices(D_bond,1,this->d_fermi,this->d_fermi));	FermiLast.permute(Indices(3,1,4,2));
  Gauge.reshape(Indices(D_bond*D_bond,1));	Gauge.multiplyRight(Gauge_operators);		Gauge.reshape(Indices(D_bond,D_bond,this->d_link,this->d_link));	Gauge.permute(Indices(3,1,4,2));
  
  int pos_odd_saved,pos_even_saved;
  bool odd_saved=false,even_saved=false;
  
  mpo.setOp(0,new Operator(FermiFirst),true);
  mpo.setOp(this->N_total-1,new Operator(FermiLast),true);
  mpo.setOp(1,new Operator(Gauge),true);
  
  for(int i=2; i<this->N_fermions; i++)
  {
    if((i%2)==0)
    {
      if(even_saved)
	mpo.setOp(2*i-2,&mpo.getOp(pos_even_saved),false);
      else
      {
	mpo.setOp(2*i-2,new Operator(FermiEven),true);
	pos_even_saved=2*i-2;
	even_saved=true;
      }
    }
    else
    {
      if(odd_saved)
	mpo.setOp(2*i-2,&mpo.getOp(pos_odd_saved),false);
      else
      {
	mpo.setOp(2*i-2,new Operator(FermiOdd),true);
	pos_odd_saved=2*i-2;
	odd_saved=true;
      }
    }
    mpo.setOp(2*i-1,&mpo.getOp(1),false);
  } 
#ifdef MATOUT
  mpo.exportForMatlab("CondensateMPO.m");
#endif
}

void ReducedSU2Hamiltonian::getGaussLawMPO(MPO &mpo, int pos) const
{
  cout << "Constructing Gauss Law MPO for position " << pos << endl;
  
  //Some error checking
  if(pos<1 || pos>this->N_fermions)
  {
    cout << "Error in ReducedSU2Hamiltonian::getGaussLawMPO(), the given position " << pos <<" is not valid" << endl;
    exit(666);
  }  
  mpo.initLength(this->N_total);
  
  //Prepare the operators (J here is simply a diagonal operator giving the l-value and J² just contains the squares (so not like the angular momentum operators)
  mwArray Jsq,J,Qsquare;
  J=mwArray(Indices(this->d_link,this->d_link));	J.fillWithZero();
  Qsquare=mwArray(Indices(this->d_fermi,this->d_fermi));	Qsquare.fillWithZero();
  
  Qsquare.setElement(0.25*ONE_c,Indices(1,1));  
  for(int i=0; i<this->d_link; i++)
    J.setElement(((double)i)/2.0*ONE_c,Indices(i,i));
  
  Jsq=J*J;
  
  ///Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(3,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //J
      Gauge_operators.setElement(J.getElement(Indices(i,j)),Indices(1,i,j));
      //J^2
      Gauge_operators.setElement(Jsq.getElement(Indices(i,j)),Indices(2,i,j));
     }
  }   
  Gauge_operators.reshape(Indices(3,this->d_link*this->d_link));  
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();    
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));      
      //Charge square
      Fermionic_operators.setElement(Qsquare.getElement(Indices(i,j)),Indices(1,i,j));      
    }
  }  
  Fermionic_operators.reshape(Indices(2,this->d_fermi*this->d_fermi));
  
  int Dbond;
  mwArray GaugeLeft,Fermi,GaugeRight;    
  //Special case for the edge
  if(pos==1)
  {
    Dbond=2;
    Fermi=mwArray(Indices(1,Dbond,2));		Fermi.fillWithZero();
    GaugeRight=mwArray(Indices(Dbond,1,3));	GaugeRight.fillWithZero();
    
    Fermi.setElement(ONE_c,Indices(0,0,0));
    Fermi.setElement(-1.0*ONE_c,Indices(0,Dbond-1,1));
    
    GaugeRight.setElement(ONE_c,Indices(Dbond-1,0,0));
    GaugeRight.setElement(ONE_c,Indices(0,0,2));
    
    Fermi.reshape(Indices(1*Dbond,2));		Fermi.multiplyRight(Fermionic_operators);	Fermi.reshape(Indices(1,Dbond,this->d_fermi,this->d_fermi));		Fermi.permute(Indices(3,1,4,2));
    GaugeRight.reshape(Indices(Dbond*1,3));	GaugeRight.multiplyRight(Gauge_operators);	GaugeRight.reshape(Indices(Dbond,1,this->d_link,this->d_link));		GaugeRight.permute(Indices(3,1,4,2));
    
    //Set operators
    mpo.setOp(2*pos-2,new Operator(Fermi),true);
    mpo.setOp(2*pos-1,new Operator(GaugeRight),true);   
  }
  else if(pos==this->N_fermions)
  {
    Dbond=2;
    GaugeLeft=mwArray(Indices(1,Dbond,3));	GaugeLeft.fillWithZero();
    Fermi=mwArray(Indices(Dbond,1,2));		Fermi.fillWithZero();
    
    GaugeLeft.setElement(ONE_c,Indices(0,0,0));
    GaugeLeft.setElement(ONE_c,Indices(0,Dbond-1,2));
    
    Fermi.setElement(ONE_c,Indices(Dbond-1,0,0));
    Fermi.setElement(-1.0*ONE_c,Indices(Dbond-1,0,1));
    
    GaugeLeft.reshape(Indices(1*Dbond,3));	GaugeLeft.multiplyRight(Gauge_operators);	GaugeLeft.reshape(Indices(1,Dbond,this->d_link,this->d_link));		GaugeLeft.permute(Indices(3,1,4,2));
    Fermi.reshape(Indices(Dbond*1,2));		Fermi.multiplyRight(Fermionic_operators);	Fermi.reshape(Indices(Dbond,1,this->d_fermi,this->d_fermi));		Fermi.permute(Indices(3,1,4,2));
    
    //Set operators
    mpo.setOp(2*(pos-1)-1,new Operator(GaugeLeft),true);
    mpo.setOp(2*pos-2,new Operator(Fermi),true);            
  }
  else
  {
    Dbond=3;
    
    GaugeLeft=mwArray(Indices(1,Dbond,3));	GaugeLeft.fillWithZero();
    GaugeRight=mwArray(Indices(Dbond,1,3));	GaugeRight.fillWithZero();
    Fermi=mwArray(Indices(Dbond,Dbond,2));	Fermi.fillWithZero();
    
    GaugeLeft.setElement(ONE_c,Indices(0,0,0));
    GaugeLeft.setElement(-2.0*ONE_c,Indices(0,1,1));
    GaugeLeft.setElement(ONE_c,Indices(0,Dbond-1,2));
    
    Fermi.setElement(ONE_c,Indices(0,0,0));
    Fermi.setElement(ONE_c,Indices(1,1,0));
    Fermi.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
    Fermi.setElement(-1.0*ONE_c,Indices(0,Dbond-1,1));
    
    GaugeRight.setElement(ONE_c,Indices(Dbond-1,0,0));
    GaugeRight.setElement(ONE_c,Indices(1,0,1));
    GaugeRight.setElement(ONE_c,Indices(0,0,2));
    
    GaugeLeft.reshape(Indices(1*Dbond,3));	GaugeLeft.multiplyRight(Gauge_operators);	GaugeLeft.reshape(Indices(1,Dbond,this->d_link,this->d_link));		GaugeLeft.permute(Indices(3,1,4,2));
    Fermi.reshape(Indices(Dbond*Dbond,2));	Fermi.multiplyRight(Fermionic_operators);	Fermi.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Fermi.permute(Indices(3,1,4,2));
    GaugeRight.reshape(Indices(Dbond*1,3));	GaugeRight.multiplyRight(Gauge_operators);	GaugeRight.reshape(Indices(Dbond,1,this->d_link,this->d_link));		GaugeRight.permute(Indices(3,1,4,2));
    
    //Set operators
    mpo.setOp(2*(pos-1)-1,new Operator(GaugeLeft),true);
    mpo.setOp(2*pos-2,new Operator(Fermi),true);
    mpo.setOp(2*pos-1,new Operator(GaugeRight),true);
  }
  
  //Add the missing identities
  int pos_id_fermi_saved,pos_id_link_saved;
  bool id_fermi_saved=false,id_link_saved=false;
  
  for(int i=1; i<=this->N_fermions; i++)
  {
    if(mpo.isEmpty(2*i-2) && id_fermi_saved)
      mpo.setOp(2*i-2,&mpo.getOp(pos_id_fermi_saved),false);
    else if(mpo.isEmpty(2*i-2))
    {
      mpo.setOp(2*i-2,new Operator(this->id_fermi_mpo),true);
      pos_id_fermi_saved=2*i-2;
      id_fermi_saved=true;
    }
    if(i<this->N_fermions && mpo.isEmpty(2*i-1) && id_link_saved)
      mpo.setOp(2*i-1,&mpo.getOp(pos_id_link_saved),false);
    else if(i<this->N_fermions && mpo.isEmpty(2*i-1))
    {
      mpo.setOp(2*i-1,new Operator(this->id_link_mpo),true);
      pos_id_link_saved=2*i-1;
      id_link_saved=true;
    }      
  }
}

void ReducedSU2Hamiltonian::getGaussLawPenaltyMPO(MPO &mpo, double strength) const
{
  mpo.initLength(this->N_total);
  if(strength==0.0)
  {
    cout << "Taking the strength from the Hamiltonian: "<< this->lambda << endl;
    strength=this->lambda;    
  }
  
  //Prepare all kinds of operators
  mwArray J,Jsq,Jcube,Jquarter,Qsquare,Qquarter;
  J=mwArray(Indices(this->d_link,this->d_link));	J.fillWithZero();
  Qsquare=mwArray(Indices(this->d_fermi,this->d_fermi));	Qsquare.fillWithZero();
  
  Qsquare.setElement(0.25*ONE_c,Indices(1,1));  
  Qquarter=Qsquare*Qsquare;
  for(int i=0; i<this->d_link; i++)
    J.setElement(((double)i)/2.0*ONE_c,Indices(i,i));  
  Jsq=J*J;
  Jcube=Jsq*J;
  Jquarter=Jcube*J;
  
  
  ///Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(5,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      Gauge_operators.setElement(J.getElement(Indices(i,j)),Indices(1,i,j));
      Gauge_operators.setElement(Jsq.getElement(Indices(i,j)),Indices(2,i,j));
      Gauge_operators.setElement(Jcube.getElement(Indices(i,j)),Indices(3,i,j));
      Gauge_operators.setElement(Jquarter.getElement(Indices(i,j)),Indices(4,i,j));
     }
  }   
  Gauge_operators.reshape(Indices(5,this->d_link*this->d_link));  
  
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(3,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();    
  for(int i=0; i<this->d_fermi; ++i)
  {
    for(int j=0; j<this->d_fermi;++j)
    {
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));      
      Fermionic_operators.setElement(Qsquare.getElement(Indices(i,j)),Indices(1,i,j));
      Fermionic_operators.setElement(Qquarter.getElement(Indices(i,j)),Indices(2,i,j));            
    }
  }  
  Fermionic_operators.reshape(Indices(3,this->d_fermi*this->d_fermi));
  
  int Dbond=7;
  mwArray FermiFirst,FermiLast,Fermi,Gauge;
  
  FermiFirst=mwArray(Indices(1,Dbond,3));	FermiFirst.fillWithZero();
  FermiLast=mwArray(Indices(Dbond,1,3));	FermiLast.fillWithZero();
  Fermi=mwArray(Indices(Dbond,Dbond,3));	Fermi.fillWithZero();
  Gauge=mwArray(Indices(Dbond,Dbond,5));	Gauge.fillWithZero();
  
  Fermi.setElement(ONE_c,Indices(0,0,0));
  Fermi.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  Fermi.setElement(ONE_c,Indices(1,1,0));
  Fermi.setElement(ONE_c,Indices(3,3,0));
  Fermi.setElement(ONE_c,Indices(4,4,0));
  
  Fermi.setElement(ONE_c,Indices(1,2,1));
  Fermi.setElement(-2.0*ONE_c,Indices(3,Dbond-1,1));
  Fermi.setElement(strength*ONE_c,Indices(0,5,1));
  
  Fermi.setElement(strength*ONE_c,Indices(0,Dbond-1,2));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  
  Gauge.setElement(strength*ONE_c,Indices(0,1,1));
  Gauge.setElement(4.0*ONE_c,Indices(2,Dbond-1,1));
  Gauge.setElement(-4.0*ONE_c,Indices(4,Dbond-1,1));
  
  Gauge.setElement(strength*ONE_c,Indices(0,3,2));
  Gauge.setElement(6.0*ONE_c,Indices(3,Dbond-1,2));
  Gauge.setElement(-2.0*ONE_c,Indices(5,Dbond-1,2));
  
  Gauge.setElement(-4.0*ONE_c,Indices(1,Dbond-1,3));
  Gauge.setElement(strength*ONE_c,Indices(0,4,3));
  
  Gauge.setElement(2.0*strength*ONE_c,Indices(0,Dbond-1,4));
  
  FermiFirst=Fermi.subArray(Indices(0,-1,-1));
  FermiLast=Fermi.subArray(Indices(-1,Dbond-1,-1));
  
  //Reshape and contract
  FermiFirst.reshape(Indices(1*Dbond,3));	FermiFirst.multiplyRight(Fermionic_operators);	FermiFirst.reshape(Indices(1,Dbond,this->d_fermi,this->d_fermi));	FermiFirst.permute(Indices(3,1,4,2));
  Fermi.reshape(Indices(Dbond*Dbond,3));	Fermi.multiplyRight(Fermionic_operators);	Fermi.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Fermi.permute(Indices(3,1,4,2));
  FermiLast.reshape(Indices(Dbond*1,3));	FermiLast.multiplyRight(Fermionic_operators);	FermiLast.reshape(Indices(Dbond,1,this->d_fermi,this->d_fermi));	FermiLast.permute(Indices(3,1,4,2));
  Gauge.reshape(Indices(Dbond*Dbond,5));	Gauge.multiplyRight(Gauge_operators);		Gauge.reshape(Indices(Dbond,Dbond,this->d_link,this->d_link));		Gauge.permute(Indices(3,1,4,2));
  
  mpo.setOp(0,new Operator(FermiFirst),true);
  mpo.setOp(1,new Operator(Gauge),true);
  mpo.setOp(this->N_total-1,new Operator(FermiLast),true);
  
  int pos_fermi_saved;
  bool fermi_saved=false;
  
  for(int i=2; i<this->N_fermions; i++)
  {
    if(fermi_saved)
      mpo.setOp(2*i-2,&mpo.getOp(pos_fermi_saved),false);
    else
    {
      mpo.setOp(2*i-2,new Operator(Fermi),true);
      pos_fermi_saved=2*i-2;
      fermi_saved=true;
    }
    mpo.setOp(2*i-1,&mpo.getOp(1),false);
  }  
  
#ifdef MATOUT
  mpo.exportForMatlab("PenaltyMPO.m");
#endif
  //cout << "Gauss Law Penalty MPO " << mpo << endl;
}

void ReducedSU2Hamiltonian::getProjectorMPO(MPO &Podd, MPO &Peven) const
{
  Podd.initLength(this->N_total);
  Peven.initLength(this->N_total);
  
  mwArray LeftEdge,RightEdge,RegularSite;
  mwArray FermiEmpty,FermiSingle,FermiDouble,FluxLeft,FluxRightlower,FluxRightupper;
  
  FermiEmpty=mwArray(Indices(this->d_fermi,this->d_fermi));	FermiEmpty.fillWithZero();	FermiEmpty.setElement(ONE_c,Indices(0,0));
  FermiSingle=mwArray(Indices(this->d_fermi,this->d_fermi));	FermiSingle.fillWithZero();	FermiSingle.setElement(ONE_c,Indices(1,1));
  FermiDouble=mwArray(Indices(this->d_fermi,this->d_fermi));	FermiDouble.fillWithZero();	FermiDouble.setElement(ONE_c,Indices(2,2));
  FluxLeft=mwArray(Indices(this->d_link,this->d_link));
  FluxRightlower=mwArray(Indices(this->d_link,this->d_link));
  FluxRightupper=mwArray(Indices(this->d_link,this->d_link));
  
  LeftEdge=mwArray(Indices(this->d_fermi*this->d_link,this->d_fermi*this->d_link));					LeftEdge.fillWithZero();
  RegularSite=mwArray(Indices(this->d_link*this->d_fermi*this->d_link,this->d_link*this->d_fermi*this->d_link));	RegularSite.fillWithZero();
  RightEdge=mwArray(Indices(this->d_link*this->d_fermi,this->d_link*this->d_fermi));					RightEdge.fillWithZero();
  
  
  for(int k=0; k<this->d_link; k++)
  {
    FluxLeft.fillWithZero();
    FluxRightupper.fillWithZero();
    FluxRightlower.fillWithZero();
    
    FluxLeft.setElement(ONE_c,Indices(k,k));
    if(k>0)
      FluxRightlower.setElement(ONE_c,Indices(k-1,k-1));
    if(k<(this->d_link-1))
      FluxRightupper.setElement(ONE_c,Indices(k+1,k+1));
    
    RegularSite = RegularSite + kron(FluxLeft,kron(FermiEmpty,FluxLeft))
                              + kron(FluxLeft,kron(FermiDouble,FluxLeft))
			      + kron(FluxLeft,kron(FermiSingle,FluxRightlower))
			      + kron(FluxLeft,kron(FermiSingle,FluxRightupper));   
  }
  //Edge terms, as the flux outside the chain is zero, there are only a few possibilities
  FluxLeft.fillWithZero();		FluxLeft.setElement(ONE_c,Indices(0,0));
  FluxRightupper.fillWithZero();	FluxRightupper.setElement(ONE_c,Indices(1,1));
  LeftEdge = LeftEdge + kron(FermiEmpty,FluxLeft) + kron(FermiDouble,FluxLeft) + kron(FermiSingle,FluxRightupper);
  RightEdge = RightEdge + kron(FluxLeft,FermiEmpty) + kron(FluxLeft,FermiDouble) + kron(FluxRightupper,FermiSingle);
  
  vector<mwArray> MatricesLeft, MatricesSite, MatricesRight;
  vector<int> dimsleft, dims, dimsright;
  dimsleft.push_back(this->d_fermi);	dimsleft.push_back(this->d_link);
  dims.push_back(this->d_link);		dims.push_back(this->d_fermi);		dims.push_back(this->d_link);
  dimsright.push_back(this->d_link);	dimsright.push_back(this->d_fermi);
  
  splitLocalOp(LeftEdge,dimsleft,MatricesLeft);
  splitLocalOp(RegularSite,dims,MatricesSite);
  splitLocalOp(RightEdge,dimsright,MatricesRight);
  
  //Now fill the MPOs
  Podd.setOp(0,new Operator(MatricesLeft[0]),true);
  Podd.setOp(1,new Operator(MatricesLeft[1]),true);

  int pos_odd_saved,pos_even_saved;
  bool odd_saved=false,even_saved=false;
  
  for(int i=3; i<=this->N_fermions; i+=2)
  {
    if(i!=this->N_fermions)
    {
      if(odd_saved)
      {
	Podd.setOp(2*(i-1)-1,&Podd.getOp(pos_odd_saved),false);
	Podd.setOp(2*i-2,&Podd.getOp(pos_odd_saved+1),false);
	Podd.setOp(2*i-1,&Podd.getOp(pos_odd_saved+2),false);
      }
      else
      {
	Podd.setOp(2*(i-1)-1,new Operator(MatricesSite[0]),true);
	Podd.setOp(2*i-2,new Operator(MatricesSite[1]),true);
	Podd.setOp(2*i-1,new Operator(MatricesSite[2]),true);
	pos_odd_saved = 2*(i-1)-1;
	odd_saved=true;
      }
    }
    else
    {
      Podd.setOp(2*(i-1)-1,new Operator(MatricesRight[0]),true);
      Podd.setOp(2*i-2,new Operator(MatricesRight[1]),true);
    }
  }
  
  for(int i=2; i<=this->N_fermions; i+=2)
  {
    if(i!=this->N_fermions)
    {
      if(even_saved)
      {
	Peven.setOp(2*(i-1)-1,&Peven.getOp(pos_even_saved),false);
	Peven.setOp(2*i-2,&Peven.getOp(pos_even_saved+1),false);
	Peven.setOp(2*i-1,&Peven.getOp(pos_even_saved+2),false);
      }
      else
      {
	Peven.setOp(2*(i-1)-1,new Operator(MatricesSite[0]),true);
	Peven.setOp(2*i-2,new Operator(MatricesSite[1]),true);
	Peven.setOp(2*i-1,new Operator(MatricesSite[2]),true);
	pos_even_saved=2*(i-1)-1;
	even_saved=true;
      }
    }
    else
    {
      Peven.setOp(2*(i-1)-1,new Operator(MatricesRight[0]),true);
      Peven.setOp(2*i-2,new Operator(MatricesRight[1]),true);
    }
  }
  
  //Remaining identities (only fermionic sites, as I already set all gauge links 
  odd_saved=false;
  even_saved=false;
  for(int i=1; i<=this->N_fermions; i++)
  {
    if(Podd.isEmpty(2*i-2) && odd_saved)
      Podd.setOp(2*i-2,&Podd.getOp(pos_odd_saved),false);
    else if(Podd.isEmpty(2*i-2))
    {
      Podd.setOp(2*i-2,new Operator(this->id_fermi_mpo),true);
      pos_odd_saved=2*i-2;
      odd_saved=true;      
    }
    if(Peven.isEmpty(2*i-2) && even_saved)
      Peven.setOp(2*i-2,&Peven.getOp(pos_even_saved),false);
    else if(Peven.isEmpty(2*i-2))
    {
      Peven.setOp(2*i-2,new Operator(this->id_fermi_mpo),true);
      pos_even_saved=2*i-2;
      even_saved=true;
    }
  }   
#ifdef MATOUT
  Podd.exportForMatlab("Podd_kron.m");
  Peven.exportForMatlab("Peven_kron.m");
#endif
}

void ReducedSU2Hamiltonian::getProjectorMPO2(MPO &Podd, MPO &Peven) const
{
  Podd.initLength(this->N_total);
  Peven.initLength(this->N_total);
  
  
  int Dbond=this->d_link;
  
  mwArray Gauge_operators,Fermionic_operators;
  Gauge_operators=mwArray(Indices(this->d_link,this->d_link,this->d_link));		Gauge_operators.fillWithZero();
  Fermionic_operators=mwArray(Indices(this->d_fermi,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  
  //Projectors on state i (d_link different ones)
  for(int i=0; i<this->d_link; i++)
    Gauge_operators.setElement(ONE_c,Indices(i,i,i));
  
  //Projectors on state i (d_fermi different ones)
  for(int i=0; i<this->d_fermi; i++)
    Fermionic_operators.setElement(ONE_c,Indices(i,i,i));
  
  //Projectors on state i (d_link different ones)
  for(int i=0; i<this->d_link; i++)
    cout << "Gauge operator " << i  << ": " << Gauge_operators.subArray(Indices(i,-1,-1)) << endl;
  
  //Projectors on state i (d_link different ones)
  for(int i=0; i<this->d_fermi; i++)
    cout << "Fermi operator " << i  << ": " << Fermionic_operators.subArray(Indices(i,-1,-1)) << endl; 
  
  cout << "Reshaping operators " << endl;
  Gauge_operators.reshape(Indices(this->d_link,this->d_link*this->d_link));
  Fermionic_operators.reshape(Indices(this->d_fermi,this->d_fermi*this->d_fermi));
  
  mwArray GaugeLeft,GaugeRight,Fermi,FermiLeft,FermiRight;
  
  GaugeLeft=mwArray(Indices(1,Dbond,this->d_link));	GaugeLeft.fillWithZero();
  Fermi=mwArray(Indices(Dbond,Dbond,this->d_fermi));	Fermi.fillWithZero();
  GaugeRight=mwArray(Indices(Dbond,1,this->d_link));	GaugeRight.fillWithZero();
  FermiLeft=mwArray(Indices(1,Dbond,this->d_fermi));	FermiLeft.fillWithZero();
  FermiRight=mwArray(Indices(Dbond,1,this->d_fermi));	FermiRight.fillWithZero();
  
  for(int i=0; i<Dbond; i++)
  {
    GaugeLeft.setElement(ONE_c,Indices(0,i,i));
    GaugeRight.setElement(ONE_c,Indices(i,0,i));
    Fermi.setElement(ONE_c,Indices(i,i,0));
    Fermi.setElement(ONE_c,Indices(i,i,2));
    if((i-1)>=0)
      Fermi.setElement(ONE_c,Indices(i,i-1,1));
    if((i+1)<Dbond)
      Fermi.setElement(ONE_c,Indices(i,i+1,1));
  }
  
  FermiLeft.setElement(ONE_c,Indices(0,0,0));
  FermiLeft.setElement(ONE_c,Indices(0,1,1));
  FermiLeft.setElement(ONE_c,Indices(0,0,2));
  
  FermiRight.setElement(ONE_c,Indices(0,0,0));
  FermiRight.setElement(ONE_c,Indices(1,0,1));
  FermiRight.setElement(ONE_c,Indices(0,0,2));
  
  for(int i=0; i<this->d_fermi; i++)
    cout << "Fermi Left " << i << ": " << FermiLeft.subArray(Indices(-1,-1,i));
  for(int i=0; i<this->d_fermi; i++)
    cout << "Fermi " << i << ": " << Fermi.subArray(Indices(-1,-1,i));
  for(int i=0; i<this->d_fermi; i++)
    cout << "Fermi Right " << i << ": " << FermiRight.subArray(Indices(-1,-1,i));
  
  for(int i=0; i<this->d_link; i++)
    cout << "Gauge Left " << i << ": " << GaugeLeft.subArray(Indices(-1,-1,i));
  for(int i=0; i<this->d_link; i++)
    cout << "Gauge Right " << i << ": " << GaugeRight.subArray(Indices(-1,-1,i));
  
  //Reshape and contract operators
  FermiLeft.reshape(Indices(1*Dbond,this->d_fermi));	FermiLeft.multiplyRight(Fermionic_operators); 	FermiLeft.reshape(Indices(1,Dbond,this->d_fermi,this->d_fermi));	FermiLeft.permute(Indices(3,1,4,2));
  FermiRight.reshape(Indices(Dbond*1,this->d_fermi));	FermiRight.multiplyRight(Fermionic_operators);	FermiRight.reshape(Indices(Dbond,1,this->d_fermi,this->d_fermi));	FermiRight.permute(Indices(3,1,4,2));
  Fermi.reshape(Indices(Dbond*Dbond,this->d_fermi));	Fermi.multiplyRight(Fermionic_operators);	Fermi.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Fermi.permute(Indices(3,1,4,2));
  GaugeLeft.reshape(Indices(1*Dbond,this->d_link));	GaugeLeft.multiplyRight(Gauge_operators);	GaugeLeft.reshape(Indices(1,Dbond,this->d_link,this->d_link));		GaugeLeft.permute(Indices(3,1,4,2));
  GaugeRight.reshape(Indices(Dbond*1,this->d_link));	GaugeRight.multiplyRight(Gauge_operators);	GaugeRight.reshape(Indices(Dbond,1,this->d_link,this->d_link));		GaugeRight.permute(Indices(3,1,4,2));  
  
  
  //Now fill the MPOs
  Podd.setOp(0,new Operator(FermiLeft),true);
  Podd.setOp(1,new Operator(GaugeRight),true);

  int pos_odd_saved,pos_even_saved;
  bool odd_saved=false,even_saved=false;
  
  for(int i=3; i<=this->N_fermions; i+=2)
  {
    if(i!=this->N_fermions)
    {
      if(odd_saved)
      {
	Podd.setOp(2*(i-1)-1,&Podd.getOp(pos_odd_saved),false);
	Podd.setOp(2*i-2,&Podd.getOp(pos_odd_saved+1),false);
	Podd.setOp(2*i-1,&Podd.getOp(pos_odd_saved+2),false);
      }
      else
      {
	Podd.setOp(2*(i-1)-1,new Operator(GaugeLeft),true);
	Podd.setOp(2*i-2,new Operator(Fermi),true);
	Podd.setOp(2*i-1,new Operator(GaugeRight),true);
	pos_odd_saved = 2*(i-1)-1;
	odd_saved=true;
      }
    }
    else
    {
      Podd.setOp(2*(i-1)-1,new Operator(GaugeLeft),true);
      Podd.setOp(2*i-2,new Operator(FermiRight),true);
    }
  }
  
  for(int i=2; i<=this->N_fermions; i+=2)
  {
    if(i!=this->N_fermions)
    {
      if(even_saved)
      {
	Peven.setOp(2*(i-1)-1,&Peven.getOp(pos_even_saved),false);
	Peven.setOp(2*i-2,&Peven.getOp(pos_even_saved+1),false);
	Peven.setOp(2*i-1,&Peven.getOp(pos_even_saved+2),false);
      }
      else
      {
	Peven.setOp(2*(i-1)-1,new Operator(GaugeLeft),true);
	Peven.setOp(2*i-2,new Operator(Fermi),true);
	Peven.setOp(2*i-1,new Operator(GaugeRight),true);
	pos_even_saved=2*(i-1)-1;
	even_saved=true;
      }
    }
    else
    {
      Peven.setOp(2*(i-1)-1,new Operator(GaugeLeft),true);
      Peven.setOp(2*i-2,new Operator(FermiRight),true);
    }
  }
  
  //Remaining identities (only fermionic sites, as I already set all gauge links 
  odd_saved=false;
  even_saved=false;
  for(int i=1; i<=this->N_fermions; i++)
  {
    if(Podd.isEmpty(2*i-2) && odd_saved)
      Podd.setOp(2*i-2,&Podd.getOp(pos_odd_saved),false);
    else if(Podd.isEmpty(2*i-2))
    {
      Podd.setOp(2*i-2,new Operator(this->id_fermi_mpo),true);
      pos_odd_saved=2*i-2;
      odd_saved=true;      
    }
    if(Peven.isEmpty(2*i-2) && even_saved)
      Peven.setOp(2*i-2,&Peven.getOp(pos_even_saved),false);
    else if(Peven.isEmpty(2*i-2))
    {
      Peven.setOp(2*i-2,new Operator(this->id_fermi_mpo),true);
      pos_even_saved=2*i-2;
      even_saved=true;
    }
  }   
#ifdef MATOUT
  Podd.exportForMatlab("Podd_mpo.m");
  Peven.exportForMatlab("Peven_mpo.m");
#endif
}


void ReducedSU2Hamiltonian::getProjectorSumMPO(MPO &Pgauss) const
{
  Pgauss.initLength(this->N_total);
  int Dbond=this->d_link+2;
  
  mwArray Gauge_operators,Fermionic_operators;
  Gauge_operators=mwArray(Indices(this->d_link+1,this->d_link,this->d_link));
  Fermionic_operators=mwArray(Indices(this->d_fermi+1,this->d_fermi,this->d_fermi));
  
  for(int i=0; i<this->d_link; i++)
  {
    //Identity
    Gauge_operators.setElement(ONE_c,Indices(0,i,i));
    //Projectors on state i (d_link different ones)
    Gauge_operators.setElement(ONE_c,Indices(i+1,i,i));
    //cout << "Gauge operator " << i+1 << ": " << Gauge_operators.subArray(Indices(i+1,-1,-1));
  }
  
  for(int i=0; i<this->d_fermi; i++)
  {
    //Identity
    Fermionic_operators.setElement(ONE_c,Indices(0,i,i));
    //Projectors on state i (d_link different ones)
    Fermionic_operators.setElement(ONE_c,Indices(i+1,i,i));
    //cout << "Fermionic operator " << i+1 << ": " << Fermionic_operators.subArray(Indices(i+1,-1,-1));
  }
  
  //For testing
  vector <mwArray> krons;
  for(int i=1; i<=this->d_link; i++)
    for(int j=1; j<=this->d_fermi; j++)
      for(int k=1;k<=this->d_link;k++)
      {
	krons.push_back(kron(Gauge_operators.subArray(Indices(i,-1,-1)),kron(Fermionic_operators.subArray(Indices(j,-1,-1)),Gauge_operators.subArray(Indices(k,-1,-1)))));
      }
      
  mwArray A,B;
  for(int i=0; i<krons.size(); i++)
  {
    A=krons[i];
    for(int k=i+1; k<krons.size();k++)
    {
      B=krons[k];
      if(A==B)
	cout << "Eintrag " << i << " ist identisch zu Eintrag " << k << endl;
    }
  }
  
  Gauge_operators.reshape(Indices(this->d_link+1,this->d_link*this->d_link));
  Fermionic_operators.reshape(Indices(this->d_fermi+1,this->d_fermi*this->d_fermi));
  
  mwArray FermiFirst,FermiLast,Fermi,Gauge;
  
  FermiFirst=mwArray(Indices(1,Dbond,this->d_fermi+1));	FermiFirst.fillWithZero();
  Fermi=mwArray(Indices(Dbond,Dbond,this->d_fermi+1));	Fermi.fillWithZero();
  FermiLast=mwArray(Indices(Dbond,1,this->d_fermi+1));	FermiLast.fillWithZero();
  Gauge=mwArray(Indices(Dbond,Dbond,this->d_link+1));		Gauge.fillWithZero();
  
  Fermi.setElement(ONE_c,Indices(0,0,0));
  Fermi.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  for(int i=1; i<(Dbond-1); i++)
  {
    Fermi.setElement(ONE_c,Indices(i,i,1));
    Fermi.setElement(ONE_c,Indices(i,i,3));
    if((i-1)>=1)
      Fermi.setElement(ONE_c,Indices(i,i-1,2));
    if((i+1)<=this->d_link)
      Fermi.setElement(ONE_c,Indices(i,i+1,2));
  }
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  for(int i=1; i<(Dbond-1); i++)
  {
    Gauge.setElement(ONE_c,Indices(0,i,i));
    Gauge.setElement(ONE_c,Indices(i,Dbond-1,i));
  }
  
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(ONE_c,Indices(0,1,1));
  FermiFirst.setElement(ONE_c,Indices(0,2,2));
  FermiFirst.setElement(ONE_c,Indices(0,1,3));
  
  FermiLast.setElement(ONE_c,Indices(Dbond-1,0,0));
  FermiLast.setElement(ONE_c,Indices(1,0,1));
  FermiLast.setElement(ONE_c,Indices(2,0,2));
  FermiLast.setElement(ONE_c,Indices(1,0,3));
  
  FermiFirst.reshape(Indices(1*Dbond,this->d_fermi+1));	FermiFirst.multiplyRight(Fermionic_operators);	FermiFirst.reshape(Indices(1,Dbond,this->d_fermi,this->d_fermi));	FermiFirst.permute(Indices(3,1,4,2));
  Fermi.reshape(Indices(Dbond*Dbond,this->d_fermi+1));	Fermi.multiplyRight(Fermionic_operators);	Fermi.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Fermi.permute(Indices(3,1,4,2));
  FermiLast.reshape(Indices(Dbond*1,this->d_fermi+1));	FermiLast.multiplyRight(Fermionic_operators);	FermiLast.reshape(Indices(Dbond,1,this->d_fermi,this->d_fermi));	FermiLast.permute(Indices(3,1,4,2));
  Gauge.reshape(Indices(Dbond*Dbond,this->d_link+1));	Gauge.multiplyRight(Gauge_operators);		Gauge.reshape(Indices(Dbond,Dbond,this->d_link,this->d_link));		Gauge.permute(Indices(3,1,4,2));
  
  Pgauss.setOp(0,new Operator(FermiFirst),true);
  Pgauss.setOp(1,new Operator(Gauge),true);
  Pgauss.setOp(this->N_total-1,new Operator(FermiLast),true);
  
  int pos_fermi_saved;
  bool fermi_saved=false;
  for(int i=2; i<(this->N_total-1); i+=2)
  {
    if(fermi_saved)
      Pgauss.setOp(i,&Pgauss.getOp(pos_fermi_saved),false);
    else      
    {
      Pgauss.setOp(i,new Operator(Fermi),true);
      pos_fermi_saved=i;
      fermi_saved=true;
    }
    Pgauss.setOp(i+1,&Pgauss.getOp(1),false);
  }
  
  
#ifdef MATOUT
  Pgauss.exportForMatlab("Pgauss.m");
#endif
  
}

void ReducedSU2Hamiltonian::getProjectorMPO(MPO &Podd, MPO &Peven,int length, int flux_quanta, int pos) const
{
  
  //Check if desired flux_quanta are possible with number of states allowed
  if(flux_quanta>=this->d_link || flux_quanta==0)
  {
    cout << "Error in ReducedSU2Hamiltonian::getProjectorMPO(), number of flux quanta is not compliant with size of link Hilbert space, will use largest flux quantum possible" << endl;
    flux_quanta=this->d_link-1;
  }  
  //Determine the starting position of the string
  if(pos<1)
    pos = round(this->N_fermions/2)-round(length/2);
  //Some error checking
  if(pos>this->N_total || pos<1 || (pos+length)>this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getProjectorMPO(), desired string of length " << length<< " starting on position " << pos << " is outside the the chain of " << this->N_fermions << " spins" << endl;
    exit(666);
  }
  
  cout << "Warning, right now only strings inside the chain work!!!!" << endl;
  cout << "Allowing strings with " << flux_quanta << " flux quanta starting at " << pos << " and " << pos+length <<  endl;
  
  Podd.initLength(this->N_total);
  Peven.initLength(this->N_total);
  
  mwArray LeftEdge,RightEdge,RegularSite,StringStart,StringEnd;
  mwArray FermiEmpty,FermiSingle,FermiDouble,FluxLeft,FluxRight,FluxRightlower,FluxRightupper;
  
  FermiEmpty=mwArray(Indices(this->d_fermi,this->d_fermi));	FermiEmpty.fillWithZero();	FermiEmpty.setElement(ONE_c,Indices(0,0));
  FermiSingle=mwArray(Indices(this->d_fermi,this->d_fermi));	FermiSingle.fillWithZero();	FermiSingle.setElement(ONE_c,Indices(1,1));
  FermiDouble=mwArray(Indices(this->d_fermi,this->d_fermi));	FermiDouble.fillWithZero();	FermiDouble.setElement(ONE_c,Indices(2,2));
  FluxLeft=mwArray(Indices(this->d_link,this->d_link));
  FluxRight=mwArray(Indices(this->d_link,this->d_link));
  FluxRightlower=mwArray(Indices(this->d_link,this->d_link));
  FluxRightupper=mwArray(Indices(this->d_link,this->d_link));
  
  LeftEdge=mwArray(Indices(this->d_fermi*this->d_link,this->d_fermi*this->d_link));					LeftEdge.fillWithZero();
  RegularSite=mwArray(Indices(this->d_link*this->d_fermi*this->d_link,this->d_link*this->d_fermi*this->d_link));	RegularSite.fillWithZero();
  StringStart=mwArray(Indices(this->d_link*this->d_fermi*this->d_link,this->d_link*this->d_fermi*this->d_link));	StringStart.fillWithZero();
  StringEnd=mwArray(Indices(this->d_link*this->d_fermi*this->d_link,this->d_link*this->d_fermi*this->d_link));		StringEnd.fillWithZero();
  RightEdge=mwArray(Indices(this->d_link*this->d_fermi,this->d_link*this->d_fermi));					RightEdge.fillWithZero();
  
  for(int k=0; k<this->d_link; k++)
  {
    //Regular Site
    FluxLeft.fillWithZero();
    FluxRight.fillWithZero();
    FluxRightupper.fillWithZero();
    FluxRightlower.fillWithZero();

    FluxLeft.setElement(ONE_c,Indices(k,k));
    FluxRight=FluxLeft;
    if(k>0)
      FluxRightlower.setElement(ONE_c,Indices(k-1,k-1));
    if(k<(this->d_link-1))
      FluxRightupper.setElement(ONE_c,Indices(k+1,k+1));
    
    RegularSite = RegularSite + kron(FluxLeft,kron(FermiEmpty,FluxRight))
                              + kron(FluxLeft,kron(FermiDouble,FluxRight))
			      + kron(FluxLeft,kron(FermiSingle,FluxRightlower))
			      + kron(FluxLeft,kron(FermiSingle,FluxRightupper)); 
			      
    //String starting point, flux to the right should be higher
    FluxRight.fillWithZero();
    FluxRightupper.fillWithZero();
    FluxRightlower.fillWithZero();
    
    //For q=2 it can happen that l+q-1 (l-1+1) are equal to l+1 (l-1) which actually means no Gauss Law violation. Therefore, besides checking if the desired flux change is possible, I also have to check if Gauss Law is really violated
    //\TODO Right now I enforce the flux on the right to be higher, but in general I think also lower is possible
    if((k+flux_quanta)<this->d_link)
      FluxRight.setElement(ONE_c,Indices(k+flux_quanta,k+flux_quanta));
    if((k-1+flux_quanta)>0 && (k-1+flux_quanta)<this->d_link  && abs(-1+flux_quanta)!=1)
      FluxRightlower.setElement(ONE_c,Indices(k-1+flux_quanta,k-1+flux_quanta));
    if((k+1+flux_quanta)>0 && (k+1+flux_quanta)<this->d_link)
      FluxRightupper.setElement(ONE_c,Indices(k+1+flux_quanta,k+1+flux_quanta));
    
    StringStart = StringStart + kron(FluxLeft,kron(FermiEmpty,FluxRight))
                              + kron(FluxLeft,kron(FermiDouble,FluxRight))
			      + kron(FluxLeft,kron(FermiSingle,FluxRightlower))
			      + kron(FluxLeft,kron(FermiSingle,FluxRightupper)); 
			      
    //String end point, flux to the right should be lower
    /*FluxRight.fillWithZero();
    FluxRightupper.fillWithZero();
    FluxRightlower.fillWithZero();
    
    if((k-flux_quanta)>0)
      FluxRight.setElement(ONE_c,Indices(k-flux_quanta,k-flux_quanta));
    if((k-1-flux_quanta)>0 && (k-1-flux_quanta)<this->d_link)
      FluxRightlower.setElement(ONE_c,Indices(k-1-flux_quanta,k-1-flux_quanta));
    if((k+1-flux_quanta)>0 && (k+1-flux_quanta)<this->d_link)
      FluxRightupper.setElement(ONE_c,Indices(k+1-flux_quanta,k+1-flux_quanta));
    
    StringEnd = StringEnd + kron(FluxLeft,kron(FermiEmpty,FluxRight))
                          + kron(FluxLeft,kron(FermiDouble,FluxRight))
			  + kron(FluxLeft,kron(FermiSingle,FluxRightlower))
			  + kron(FluxLeft,kron(FermiSingle,FluxRightupper)); */ 
    StringEnd = StringEnd + kron(FluxRight,kron(FermiEmpty,FluxLeft))
                          + kron(FluxRight,kron(FermiDouble,FluxLeft))
			  + kron(FluxRightlower,kron(FermiSingle,FluxLeft))
			  + kron(FluxRightupper,kron(FermiSingle,FluxLeft));
    
  }
  cout << "Warning, due to testing will replace StringEnd with String Start" << endl;
  //StringStart=StringEnd;
  //Edge terms, as the flux outside the chain is zero, there are only a few possibilities
  FluxLeft.fillWithZero();		FluxLeft.setElement(ONE_c,Indices(0,0));
  FluxRightupper.fillWithZero();	FluxRightupper.setElement(ONE_c,Indices(1,1));
  LeftEdge =  kron(FermiEmpty,FluxLeft) + kron(FermiDouble,FluxLeft) + kron(FermiSingle,FluxRightupper);
  RightEdge = kron(FluxLeft,FermiEmpty) + kron(FluxLeft,FermiDouble) + kron(FluxRightupper,FermiSingle);
  
  vector<mwArray> MatricesLeft, MatricesSite, MatricesRight, MatricesStart, MatricesEnd;
  vector<int> dimsleft, dims, dimsright;
  dimsleft.push_back(this->d_fermi);	dimsleft.push_back(this->d_link);
  dims.push_back(this->d_link);		dims.push_back(this->d_fermi);		dims.push_back(this->d_link);
  dimsright.push_back(this->d_link);	dimsright.push_back(this->d_fermi);
  
  splitLocalOp(LeftEdge,dimsleft,MatricesLeft);
  splitLocalOp(RegularSite,dims,MatricesSite);
  splitLocalOp(StringStart,dims,MatricesStart);
  splitLocalOp(StringEnd,dims,MatricesEnd);
  splitLocalOp(RightEdge,dimsright,MatricesRight);
  
  //Now fill the MPOs
  Podd.setOp(0,new Operator(MatricesLeft[0]),true);
  Podd.setOp(1,new Operator(MatricesLeft[1]),true);

  int pos_odd_saved,pos_even_saved;
  bool odd_saved=false,even_saved=false;
  
  for(int i=3; i<=this->N_fermions; i+=2)
  {
    if(i==pos)
    {
      cout << "Setting odd string projector at site " << i << endl;
      Podd.setOp(2*(i-1)-1,new Operator(MatricesStart[0]),true);
      Podd.setOp(2*i-2,new Operator(MatricesStart[1]),true);
      Podd.setOp(2*i-1,new Operator(MatricesStart[2]),true);
    }
    else if( i==(pos+length))
    {
      cout << "Setting odd string projector at site " << i << endl;
      Podd.setOp(2*(i-1)-1,new Operator(MatricesEnd[0]),true);
      Podd.setOp(2*i-2,new Operator(MatricesEnd[1]),true);
      Podd.setOp(2*i-1,new Operator(MatricesEnd[2]),true);
    }
    else if(i!=this->N_fermions)
    {
      if(odd_saved)
      {
	Podd.setOp(2*(i-1)-1,&Podd.getOp(pos_odd_saved),false);
	Podd.setOp(2*i-2,&Podd.getOp(pos_odd_saved+1),false);
	Podd.setOp(2*i-1,&Podd.getOp(pos_odd_saved+2),false);
      }
      else
      {
	Podd.setOp(2*(i-1)-1,new Operator(MatricesSite[0]),true);
	Podd.setOp(2*i-2,new Operator(MatricesSite[1]),true);
	Podd.setOp(2*i-1,new Operator(MatricesSite[2]),true);
	pos_odd_saved = 2*(i-1)-1;
	odd_saved=true;
      }
    }
    else
    {
      Podd.setOp(2*(i-1)-1,new Operator(MatricesRight[0]),true);
      Podd.setOp(2*i-2,new Operator(MatricesRight[1]),true);
    }
  }
  
  for(int i=2; i<=this->N_fermions; i+=2)
  {
    if(i==pos)
    {
      cout << "Setting even string projector at site " << i << endl;      
      Peven.setOp(2*(i-1)-1,new Operator(MatricesStart[0]),true);
      Peven.setOp(2*i-2,new Operator(MatricesStart[1]),true);
      Peven.setOp(2*i-1,new Operator(MatricesStart[2]),true);
    }
    else if(i==(pos+length))
    {
      cout << "Setting even string projector at site " << i << endl;      
      Peven.setOp(2*(i-1)-1,new Operator(MatricesEnd[0]),true);
      Peven.setOp(2*i-2,new Operator(MatricesEnd[1]),true);
      Peven.setOp(2*i-1,new Operator(MatricesEnd[2]),true);
    }
    else if(i!=this->N_fermions)
    {
      if(even_saved)
      {
	Peven.setOp(2*(i-1)-1,&Peven.getOp(pos_even_saved),false);
	Peven.setOp(2*i-2,&Peven.getOp(pos_even_saved+1),false);
	Peven.setOp(2*i-1,&Peven.getOp(pos_even_saved+2),false);
      }
      else
      {
	Peven.setOp(2*(i-1)-1,new Operator(MatricesSite[0]),true);
	Peven.setOp(2*i-2,new Operator(MatricesSite[1]),true);
	Peven.setOp(2*i-1,new Operator(MatricesSite[2]),true);
	pos_even_saved=2*(i-1)-1;
	even_saved=true;
      }
    }
    else
    {
      Peven.setOp(2*(i-1)-1,new Operator(MatricesRight[0]),true);
      Peven.setOp(2*i-2,new Operator(MatricesRight[1]),true);
    }
  }
  
  //Remaining identities (only fermionic sites, as I already set all gauge links)
  odd_saved=false;
  even_saved=false;
  for(int i=1; i<=this->N_fermions; i++)
  {
    if(Podd.isEmpty(2*i-2) && odd_saved)
      Podd.setOp(2*i-2,&Podd.getOp(pos_odd_saved),false);
    else if(Podd.isEmpty(2*i-2))
    {
      Podd.setOp(2*i-2,new Operator(this->id_fermi_mpo),true);
      pos_odd_saved=2*i-2;
      odd_saved=true;      
    }
    if(Peven.isEmpty(2*i-2) && even_saved)
      Peven.setOp(2*i-2,&Peven.getOp(pos_even_saved),false);
    else if(Peven.isEmpty(2*i-2))
    {
      Peven.setOp(2*i-2,new Operator(this->id_fermi_mpo),true);
      pos_even_saved=2*i-2;
      even_saved=true;
    }
  }   
#ifdef MATOUT
  MPO dummy(3);
  dummy.setOp(0,new Operator(MatricesStart[0]),true);
  dummy.setOp(1,new Operator(MatricesStart[1]),true);
  dummy.setOp(2,new Operator(MatricesStart[2]),true);
  dummy.exportForMatlab("Pstringsmall_start.m");
  
  dummy.clear();
  dummy.initLength(3);
  dummy.setOp(0,new Operator(MatricesEnd[0]),true);
  dummy.setOp(1,new Operator(MatricesEnd[1]),true);
  dummy.setOp(2,new Operator(MatricesEnd[2]),true);
  dummy.exportForMatlab("Pstringsmall_end.m");
  
  Podd.exportForMatlab("Podd_string_kron.m");
  Peven.exportForMatlab("Peven_string_kron.m");
  
#endif
}

void ReducedSU2Hamiltonian::getProjectorMPO2(MPO &Podd, MPO &Peven,int length, int flux_quanta) const
{ 
  int pos = round(this->N_fermions/2)-round(length/2);
  //Some error checking
  if((pos>this->N_fermions || pos<1 || (pos+length)>this->N_fermions) && flux_quanta>0)
  {
    cout << "Error in SU2Hamiltonian::getProjectorMPO(), desired string of length " << length<< " starting on position " << pos << " is outside the the chain of " << this->N_fermions << " spins" << endl;
    exit(666);
  }
  vector <string_t> config;
  //Generate string configuartion if one is given, if not, simply take empty vector
  if(flux_quanta>0 && pos>0)    
    config.push_back((string_t) {pos,length}); 
  
  this->getProjectorMPO2(Podd,Peven,config,flux_quanta);  
}

void ReducedSU2Hamiltonian::getProjectorMPO2(MPO &Podd, MPO &Peven, vector<string_t> str_config, int flux_quanta) const
{
  //cout << "Warning, mixed terms deactivated" << endl;
  
  //Check if desired amount of flux_quanta are possible with number of states allowed
  if(flux_quanta>=this->d_link || flux_quanta<0)
  {
    cout << "Error in ReducedSU2Hamiltonian::getProjectorMPO(), number of flux quanta is not compliant with size of link Hilbert space, will use largest flux quantum possible" << endl;
    flux_quanta=this->d_link-1;
  }  
  
  if(flux_quanta>0 && str_config.size()>0)
    cout << "Constructing projectors on the gauge invariant subspace with external charges having the value " << flux_quanta <<  endl;
  else
  {
    //First make sure that pos is really 0 if flux_quanta is 0 such that later on while setting the operators I do not set an uninitialized one
    cout << "Constructing projectors on the gauge invariant subspace without external charges" << endl;
  }
  
  //TODO Check that strings are not overlapping
  
  Podd.initLength(this->N_total);
  Peven.initLength(this->N_total);
  
  mwArray Fermionic_operators,Gauge_operators;
  
  Fermionic_operators=mwArray(Indices(this->d_fermi,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(this->d_link,this->d_link,this->d_link));		Gauge_operators.fillWithZero();
  
  for(int i=0;i<this->d_fermi; i++)
    Fermionic_operators.setElement(ONE_c,Indices(i,i,i));
  for(int i=0; i<this->d_link; i++)
    Gauge_operators.setElement(ONE_c,Indices(i,i,i));
  
  Fermionic_operators.reshape(Indices(this->d_fermi,this->d_fermi*this->d_fermi));
  Gauge_operators.reshape(Indices(this->d_link,this->d_link*this->d_link));
  
  int Dbond = this->d_link;
  mwArray FermiFirst,FermiLast,Fermi,GaugeLeft,GaugeRight,FermiFirstString,FermiLastString,FermiString;
  
  FermiFirst=mwArray(Indices(1,Dbond,this->d_fermi));		FermiFirst.fillWithZero();
  FermiFirstString=mwArray(Indices(1,Dbond,this->d_fermi));	FermiFirstString.fillWithZero();
  Fermi=mwArray(Indices(Dbond,Dbond,this->d_fermi));		Fermi.fillWithZero();
  FermiLastString=mwArray(Indices(Dbond,1,this->d_fermi));	FermiLastString.fillWithZero();
  FermiLast=mwArray(Indices(Dbond,1,this->d_fermi));		FermiLast.fillWithZero();
  GaugeLeft=mwArray(Indices(1,Dbond,this->d_link));	GaugeLeft.fillWithZero();
  GaugeRight=mwArray(Indices(Dbond,1,this->d_link));GaugeRight.fillWithZero();
  FermiString=mwArray(Indices(Dbond,Dbond,this->d_fermi));	FermiString.fillWithZero();
  
  
  //String is starting on the first site, FermiLeft has to be set accordingly
  FermiFirstString.setElement(ONE_c,Indices(0,flux_quanta,0));
  if((flux_quanta+1)>=0 && (flux_quanta+1)<Dbond)
    FermiFirstString.setElement(ONE_c,Indices(0,flux_quanta+1,1));
  if((flux_quanta-1)>=0 && (flux_quanta-1)<Dbond && abs(flux_quanta-1)!=1)
    FermiFirstString.setElement(ONE_c,Indices(0,flux_quanta-1,1));
  if((-flux_quanta+1)>=0 && (-flux_quanta+1)<Dbond && abs(-flux_quanta+1)!=1)
    FermiFirstString.setElement(ONE_c,Indices(0,-flux_quanta+1,1));
  if((-flux_quanta-1)>=0 && (-flux_quanta-1)<Dbond)
    FermiFirstString.setElement(ONE_c,Indices(0,-flux_quanta-1,1));
  FermiFirstString.setElement(ONE_c,Indices(0,flux_quanta,2));

  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(ONE_c,Indices(0,1,1));
  FermiFirst.setElement(ONE_c,Indices(0,0,2));
  
  for(int i=0; i<Dbond; i++)
  {
    GaugeLeft.setElement(ONE_c,Indices(0,i,i));
    GaugeRight.setElement(ONE_c,Indices(i,0,i));
    Fermi.setElement(ONE_c,Indices(i,i,0));
    Fermi.setElement(ONE_c,Indices(i,i,2));
    if((i+1)<Dbond)
      Fermi.setElement(ONE_c,Indices(i,i+1,1));
    if((i-1)>=0)
      Fermi.setElement(ONE_c,Indices(i,i-1,1));      
    if(flux_quanta>0)
    {
      //Fermionic site where the string starts/ends
      if((i+flux_quanta)<Dbond)
      {
	FermiString.setElement(ONE_c,Indices(i,i+flux_quanta,0));
	FermiString.setElement(ONE_c,Indices(i,i+flux_quanta,2));
      }
      if((i-flux_quanta)>=0)
      {
	FermiString.setElement(ONE_c,Indices(i,i-flux_quanta,0));
	FermiString.setElement(ONE_c,Indices(i,i-flux_quanta,2));
      }
      //For q=2 can happen that l+q-1 and l-q+1 result in l+1 and l+1 or l-1 which actually means no Gauss Law violation, so besides checking if it is possible to raise (lower) the flux by the desired amount of quanta, I also have to make sure, that it is really Gauss Law violating 
      if((i+flux_quanta+1)>=0 && (i+flux_quanta+1)<Dbond)
	FermiString.setElement(ONE_c,Indices(i,i+flux_quanta+1,1));
      if((i+flux_quanta-1)>=0 && (i+flux_quanta-1)<Dbond && abs(flux_quanta-1)!=1)
	FermiString.setElement(ONE_c,Indices(i,i+flux_quanta-1,1));
      if((i-flux_quanta+1)>=0 && (i-flux_quanta+1)<Dbond  && abs(-flux_quanta+1)!=1)
	FermiString.setElement(ONE_c,Indices(i,i-flux_quanta+1,1));
      if((i-flux_quanta-1)>=0 && (i-flux_quanta-1)<Dbond)
	FermiString.setElement(ONE_c,Indices(i,i-flux_quanta-1,1));
    }
  }
  
  //String at the last site, FermiLast has to be set accordingly
  FermiLastString.setElement(ONE_c,Indices(flux_quanta,0,0));
  if((flux_quanta+1)>=0 && (flux_quanta+1)<Dbond)
    FermiLastString.setElement(ONE_c,Indices(flux_quanta+1,0,1));
  if((flux_quanta-1)>=0 && (flux_quanta-1)<Dbond && abs(flux_quanta-1)!=1)
    FermiLastString.setElement(ONE_c,Indices(flux_quanta-1,0,1));
  if((-flux_quanta+1)>=0 && (-flux_quanta+1)<Dbond  && abs(-flux_quanta+1)!=1)
    FermiLastString.setElement(ONE_c,Indices(-flux_quanta+1,0,1));
  if((-flux_quanta-1)>=0 && (-flux_quanta-1)<Dbond)
    FermiLastString.setElement(ONE_c,Indices(-flux_quanta-1,0,1));
  FermiLast.setElement(ONE_c,Indices(flux_quanta,0,2));
  
  //No string at last site
  FermiLast.setElement(ONE_c,Indices(0,0,0));
  FermiLast.setElement(ONE_c,Indices(1,0,1));
  FermiLast.setElement(ONE_c,Indices(0,0,2));
  
  //Reshape and contract
  FermiFirst.reshape(Indices(1*Dbond,this->d_fermi));	FermiFirst.multiplyRight(Fermionic_operators);	FermiFirst.reshape(Indices(1,Dbond,this->d_fermi,this->d_fermi));	FermiFirst.permute(Indices(3,1,4,2));
  Fermi.reshape(Indices(Dbond*Dbond,this->d_fermi));	Fermi.multiplyRight(Fermionic_operators);	Fermi.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Fermi.permute(Indices(3,1,4,2));
  FermiLast.reshape(Indices(Dbond*1,this->d_fermi));	FermiLast.multiplyRight(Fermionic_operators);	FermiLast.reshape(Indices(Dbond,1,this->d_fermi,this->d_fermi));	FermiLast.permute(Indices(3,1,4,2));
  GaugeLeft.reshape(Indices(1*Dbond,this->d_link));	GaugeLeft.multiplyRight(Gauge_operators);	GaugeLeft.reshape(Indices(1,Dbond,this->d_link,this->d_link));		GaugeLeft.permute(Indices(3,1,4,2));
  GaugeRight.reshape(Indices(Dbond*1,this->d_link));	GaugeRight.multiplyRight(Gauge_operators);	GaugeRight.reshape(Indices(Dbond,1,this->d_link,this->d_link));		GaugeRight.permute(Indices(3,1,4,2));
  
  FermiFirstString.reshape(Indices(1*Dbond,this->d_fermi));
  FermiFirstString.multiplyRight(Fermionic_operators);
  FermiFirstString.reshape(Indices(1,Dbond,this->d_fermi,this->d_fermi));
  FermiFirstString.permute(Indices(3,1,4,2));
  
  FermiString.reshape(Indices(Dbond*Dbond,this->d_fermi));
  FermiString.multiplyRight(Fermionic_operators);
  FermiString.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));
  FermiString.permute(Indices(3,1,4,2));
  
  FermiLastString.reshape(Indices(Dbond*1,this->d_fermi));
  FermiLastString.multiplyRight(Fermionic_operators);
  FermiLastString.reshape(Indices(Dbond,1,this->d_fermi,this->d_fermi));
  FermiLastString.permute(Indices(3,1,4,2));  
  
  //First ignore the string and simply save all matrices projecting on the Gauss Law fulfilling subspace, later on replace the ones
  Podd.setOp(0,new Operator(FermiFirst),true);
  Podd.setOp(1,new Operator(GaugeRight),true);  
  //Rest in between
  for(int i=3; i<this->N_fermions; i+=2)
  {    
    Podd.setOp(2*(i-1)-1,new Operator(GaugeLeft),true);
    Podd.setOp(2*i-2,new Operator(Fermi),true);
    Podd.setOp(2*i-1,new Operator(GaugeRight),true);
  }  
  //Last entry if necessary
  if((this->N_fermions%2)!=0)
  {
    Podd.setOp(this->N_total-2,new Operator(GaugeLeft),true);
    Podd.setOp(this->N_total-1,new Operator(FermiLast),true);
  }
  
  for(int i=2; i<this->N_fermions; i+=2)
  {
    Peven.setOp(2*(i-1)-1,new Operator(GaugeLeft),true);
    Peven.setOp(2*i-2,new Operator(Fermi),true);
    Peven.setOp(2*i-1,new Operator(GaugeRight),true);
  }  
  //Last entry if necessary
  if((this->N_fermions%2)==0)
  {
    Peven.setOp(this->N_total-2,new Operator(GaugeLeft),true);
    Peven.setOp(this->N_total-1,new Operator(FermiLast),true);
  }
  
  //Now fill in the missing fermionic identities
  bool id_even_saved=false,id_odd_saved=false;
  int pos_even_saved,pos_odd_saved;
  for(int i=0; i<this->N_total; i+=2)
  {
    if(Podd.isEmpty(i) && id_odd_saved)
      Podd.setOp(i,&Podd.getOp(pos_odd_saved),false);
    else if(Podd.isEmpty(i))
    {
      Podd.setOp(i,new Operator(this->id_fermi_mpo),true);
      pos_odd_saved=i;
      id_odd_saved=true;
    }
    if(Peven.isEmpty(i) && id_even_saved)
      Peven.setOp(i,&Peven.getOp(pos_even_saved),false);
    else if(Peven.isEmpty(i))
    {
      Peven.setOp(i,new Operator(this->id_fermi_mpo),true);
      pos_even_saved=i;
      id_even_saved=true;
    }
  }
  
  //Now go over the string configurations and replace the entries accordingly
  int start,end;
  for(int i=0; i<str_config.size(); i++)
  {
    //Extract pos and leng
    start=str_config[i].pos;
    end=start + str_config[i].leng;
    if((start%2)==0)
    {
      //Put it in Peven
      Peven.setOp(2*start-2,new Operator(FermiString),true);
    }
    else
    {
      //Put it in Podd, now look if it is the first site or not
      if(start==1)
	Podd.setOp(2*start-2,new Operator(FermiFirstString),true);
      else
	Podd.setOp(2*start-2,new Operator(FermiString),true);
    }
    if((end%2)==0)
    {
      //Put it in Peven
      if(end==this->N_fermions)
	Peven.setOp(2*end-2,new Operator(FermiLastString),true);
      else
	Peven.setOp(2*end-2,new Operator(FermiString),true);
    }
    else
    {
      //Put it in Podd
      if(end==this->N_fermions)
	Podd.setOp(2*end-2,new Operator(FermiLastString),true);
      else
	Podd.setOp(2*end-2,new Operator(FermiString),true);
    }    
  }
#ifdef MATOUT
   Podd.exportForMatlab("Podd_string_mpo.m");
   Peven.exportForMatlab("Peven_string_mpo.m");  
#endif
}

