/**
   \file Z2Hamiltonian.cpp Definition of the class that implements the Z2 Hamiltonian in spin formulation
   
   \author Stefan Kühn
   \date 14/12/2017
*/

#include "Z2Hamiltonian.h"
#include "Indices.h"
#include "splitLocalOp.h"
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>

using namespace std;
using namespace shrt;


//Standard Constructors
Z2Hamiltonian::Z2Hamiltonian(int L_,double x_,double mu_, double lambda_):hamiltonian(2*L_-1)
{
  ///Number of fermions
  this->L=L_;
  //Hopping constant
  this->x=x_;  
  //Mass
  this->mu=mu_;
  //Constant for the penalty term
  this->lambda=lambda_;
  //Set values for boundary spins
  this->c_left=1.0;
  this->c_right=1.0;
  //No background field, simply set to zero
  this->alpha=0.0;
  
  //Construct the basic operators
  initOperators();
  //Construct the actual Hamiltonian MPO (now with the new routine overwriting the parent one)
  initHamiltonian();
}

Z2Hamiltonian::Z2Hamiltonian(int L_,double x_,double mu_, double alpha_, double lambda_):hamiltonian(2*L_-1)
{
  ///Number of fermions
  this->L=L_;
  //Hopping constant
  this->x=x_;  
  //Mass
  this->mu=mu_;
  //Constant for the penalty term
  this->lambda=lambda_;
  //Set values for boundary spins
  this->c_left=1.0;
  this->c_right=1.0;
  //With background field
  this->alpha=alpha_;
  
  //Construct the basic operators
  initOperators();
  //Construct the actual Hamiltonian MPO (now with the new routine overwriting the parent one)
  initHamiltonian();
}

Z2Hamiltonian::Z2Hamiltonian(int L_,double x_,double mu_):hamiltonian(2*L_-1)
{
  //Number of fermions
  this->L=L_;
  //Hopping constant
  this->x=x_;  
  //Mass
  this->mu=mu_;
  
  //Construct the basic operators
  initOperators();
  //Construct the actual Hamiltonian MPO (now with the new routine overwriting the parent one)
  initHamiltonianNoPenalty();
}

//Destructor
Z2Hamiltonian::~Z2Hamiltonian(){}

void Z2Hamiltonian::initOperators()
{
  Z=mwArray(Indices(7,2,2));
  
  //Identity Matrix
  mwArray id=identityMatrix(2);  
  //sigma_x
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sig_x(Indices(2,2),datax);
  //sigma_y
  complex_t datay[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sig_y(Indices(2,2),datay);
  //sigma_z
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
  mwArray sig_z(Indices(2,2),dataz);
  //sigma_+
  mwArray sig_p=0.5*(sig_x+I_c*sig_y);
  //sigma_-
  mwArray sig_m=0.5*(sig_x-I_c*sig_y);
  //P+P^\dagger for the background field case
  complex_t datael[]={2.0*ONE_c*cos(M_PIl*this->alpha),ZERO_c,ZERO_c,-2.0*ONE_c*cos(M_PIl*this->alpha)};
  mwArray el(Indices(2,2),datael);
  
  cout << "el = " << el << endl;
  
  for(int i=0; i<2; i++)
    for(int k=0; k<2; k++)
    {
      Z.setElement(id.getElement(Indices(i,k)),Indices(0,i,k));
      Z.setElement(sig_p.getElement(Indices(i,k)),Indices(1,i,k));
      Z.setElement(sig_m.getElement(Indices(i,k)),Indices(2,i,k));
      Z.setElement(sig_z.getElement(Indices(i,k)),Indices(3,i,k));
      Z.setElement(sig_y.getElement(Indices(i,k)),Indices(4,i,k));
      Z.setElement(sig_x.getElement(Indices(i,k)),Indices(5,i,k)); 
      Z.setElement(el.getElement(Indices(i,k)),Indices(6,i,k)); 
    }   

  //I need often identities to construct a single body operator, therefore I prepare those and save them
  this->id_mpo=id;
  this->id_mpo.reshape(Indices(2,1,2,1));
}

//Build the Hamiltonian MPO
void Z2Hamiltonian::initHamiltonianNoPenalty()
{
  cout << "Still to come" << endl;  
}

//Build the Hamiltonian MPO
void Z2Hamiltonian::initHamiltonian()
{
  
  //Set up Hamiltonian
  this->hamiltonian.initLength(this->L*2-1);
  //Compared to the model without charge penalty, my bond dimension is increased by one
  int Dbond=5;
  
  mwArray site_first, site_odd, site_even, site_last, link;
  
  site_first=mwArray(Indices(1,Dbond,7));	site_first.fillWithZero();
  site_last=mwArray(Indices(Dbond,1,7));	site_last.fillWithZero();
  site_odd=mwArray(Indices(Dbond,Dbond,7));	site_odd.fillWithZero();
  site_even=mwArray(Indices(Dbond,Dbond,7));	site_even.fillWithZero();
  link=mwArray(Indices(Dbond,Dbond,7));		link.fillWithZero();
  
  //identity
  site_first.setElement(ONE_c,Indices(0,0,0));
  site_first.setElement(ONE_c*(this->mu/2.0+2.0*this->lambda),Indices(0,Dbond-1,0));
  
  site_odd.setElement(ONE_c,Indices(0,0,0));
  site_odd.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  site_odd.setElement(ONE_c*(-1.0*this->mu/2.0+2.0*this->lambda),Indices(0,Dbond-1,0));
  
  site_even.setElement(ONE_c,Indices(0,0,0));
  site_even.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  site_even.setElement(ONE_c*(this->mu/2.0+2.0*this->lambda),Indices(0,Dbond-1,0));
  
  site_last.setElement(2.0*ONE_c*this->lambda+pow(-1.0,this->L-1)*ONE_c*this->mu/2.0,Indices(0,0,0));
  site_last.setElement(ONE_c,Indices(Dbond-1,0,0));
  
  link.setElement(ONE_c,Indices(0,0,0));
  link.setElement(2.0*ONE_c,Indices(0,Dbond-1,0));
  link.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  
  //sigma_plus
  site_first.setElement(-this->x*ONE_c,Indices(0,1,1));
  
  site_odd.setElement(-this->x*ONE_c,Indices(0,1,1));
  site_odd.setElement(ONE_c,Indices(2,Dbond-1,1));
  
  site_even.setElement(-this->x*ONE_c,Indices(0,1,1));
  site_even.setElement(ONE_c,Indices(2,Dbond-1,1));
  
  site_last.setElement(ONE_c,Indices(2,0,1));
  
  //sigma_minus
  site_first.setElement(-this->x*ONE_c,Indices(0,2,2));
  
  site_odd.setElement(-this->x*ONE_c,Indices(0,2,2));
  site_odd.setElement(ONE_c,Indices(1,Dbond-1,2));
  
  site_even.setElement(-this->x*ONE_c,Indices(0,2,2));
  site_even.setElement(ONE_c,Indices(1,Dbond-1,2));
  
  site_last.setElement(ONE_c,Indices(1,0,2));
  
  //sigma_z
  site_first.setElement(ONE_c*this->lambda*2.0*this->c_left,Indices(0,3,3));
  site_first.setElement(ONE_c*this->mu/2.0,Indices(0,Dbond-1,3));
  
  site_odd.setElement(-1.0*ONE_c*this->mu/2.0,Indices(0,Dbond-1,3));
  site_odd.setElement(-2.0*ONE_c*this->lambda,Indices(3,3,3));
  
  site_even.setElement(1.0*ONE_c*this->mu/2.0,Indices(0,Dbond-1,3));
  site_even.setElement(2.0*ONE_c*this->lambda,Indices(3,3,3));
  
  site_last.setElement(pow(-1.0,this->L-1)*ONE_c*this->mu/2.0,Indices(0,0,3));
  site_last.setElement(pow(-1.0,this->L-1)*2.0*ONE_c*this->lambda*this->c_right,Indices(3,0,3));
  
  link.setElement(ONE_c,Indices(0,3,3));
  link.setElement(ONE_c,Indices(3,Dbond-1,3));
  
  //sigma_x
  link.setElement(ONE_c,Indices(1,1,5));
  link.setElement(ONE_c,Indices(2,2,5));  
  
  //"Electric field" operator
  link.setElement(-1.0*ONE_c,Indices(0,Dbond-1,6));
  
  //Reshape and contract
  site_first.reshape(Indices(1*Dbond,7));
  site_odd.reshape(Indices(Dbond*Dbond,7));
  site_even.reshape(Indices(Dbond*Dbond,7));
  site_last.reshape(Indices(Dbond*1,7));
  link.reshape(Indices(Dbond*Dbond,7));  
  
  site_first.multiplyRight(reshape(Z,Indices(7,2*2)));
  site_odd.multiplyRight(reshape(Z,Indices(7,2*2)));
  site_even.multiplyRight(reshape(Z,Indices(7,2*2)));
  site_last.multiplyRight(reshape(Z,Indices(7,2*2)));
  link.multiplyRight(reshape(Z,Indices(7,2*2)));
  
  site_first.reshape(Indices(1,Dbond,2,2));
  site_odd.reshape(Indices(Dbond,Dbond,2,2));
  site_even.reshape(Indices(Dbond,Dbond,2,2));
  site_last.reshape(Indices(Dbond,1,2,2));
  link.reshape(Indices(Dbond,Dbond,2,2));  
  
  site_first.permute(Indices(3,1,4,2));
  site_odd.permute(Indices(3,1,4,2));
  site_even.permute(Indices(3,1,4,2));
  site_last.permute(Indices(3,1,4,2));
  link.permute(Indices(3,1,4,2));
  
  
  //Now fill the MPO
  this->hamiltonian.setOp(0,new Operator(site_first),true);
  this->hamiltonian.setOp(1,new Operator(link),true);
  this->hamiltonian.setOp(2*this->L-2,new Operator(site_last),true);
  
  //Links first
  for(int i=3; i<2*this->L-2; i+=2)
    this->hamiltonian.setOp(i,&this->hamiltonian.getOp(1),false);
   
  //Sites in between the edges
  bool odd_saved=false, even_saved=false;
  int ind_odd_saved, ind_even_saved;
  for(int i=1; i<this->L-1; i++)
  {
   if((i%2)==0)
   {
     //Even site 
     if(even_saved)
       this->hamiltonian.setOp(2*i,&this->hamiltonian.getOp(ind_even_saved),false);
     else
     {
       this->hamiltonian.setOp(2*i,new Operator(site_even),true);
       ind_even_saved=2*i;
       even_saved=true;
     }
   }
   else
   {
     //Even site 
     if(odd_saved)
       this->hamiltonian.setOp(2*i,&this->hamiltonian.getOp(ind_odd_saved),false);
     else
     {
       this->hamiltonian.setOp(2*i,new Operator(site_odd),true);
       ind_odd_saved=2*i;
       odd_saved=true;
     }
   }
  }
  
#ifdef MATOUT
  this->hamiltonian.exportForMatlab("Z2Hamiltonian.m",11);
#endif
}

void Z2Hamiltonian::getSingleBodyMPO(MPO &mpo, int index, SingleBodyOperator Op) const
{
  //Frist clear the MPO and resize it
  mpo.clear();
  mpo.initLength(2*this->L-1);
  
  mwArray LocalOperator;
  
     
   switch(Op)
  {
    case spOp:
      LocalOperator=this->Z.subArray(Indices(-1,-1,1)); break;
    case smOp:
      LocalOperator=this->Z.subArray(Indices(-1,-1,2)); break;
    case szOp:
      LocalOperator=this->Z.subArray(Indices(-1,-1,3)); break;
    case syOp:
      LocalOperator=this->Z.subArray(Indices(-1,-1,4)); break;
    case sxOp:
      LocalOperator=this->Z.subArray(Indices(-1,-1,5)); break;
    default:
      cout << "Unknown operator" << endl;
      exit(666);
  }
   LocalOperator.reshape(Indices(2,1,2,1));

    //Now set identities
   bool saved=false;
   int saved_site;
   for(int k=0; k<2*this->L-1; k++)
   {
     if(saved && (k!=index))
      mpo.setOp(k,&mpo.getOp(saved_site),false);
    else
    {
      if(k==index)
	mpo.setOp(k,new Operator(LocalOperator),true);      
      else
      {
	mpo.setOp(k,new Operator(this->id_mpo),true);
	saved=true;
	saved_site=k;
      }
    }
  }
}

//Get the MPO for the penalty enforching the gauge invariant sector
void Z2Hamiltonian::getPenaltyMPO(MPO &mpo) const
{
  //Set up MPO
  mpo.initLength(this->L*2-1);
  int Dbond=3;
  
  mwArray site_first, site_odd, site_even, site_last, link;
  
  site_first=mwArray(Indices(1,Dbond,7));	site_first.fillWithZero();
  site_last=mwArray(Indices(Dbond,1,7));	site_last.fillWithZero();
  site_odd=mwArray(Indices(Dbond,Dbond,7));	site_odd.fillWithZero();
  site_even=mwArray(Indices(Dbond,Dbond,7));	site_even.fillWithZero();
  link=mwArray(Indices(Dbond,Dbond,7));		link.fillWithZero();
  
  //identity
  site_first.setElement(ONE_c,Indices(0,0,0));
  site_first.setElement(2.0*ONE_c,Indices(0,Dbond-1,0));
  
  site_odd.setElement(ONE_c,Indices(0,0,0));
  site_odd.setElement(2.0*ONE_c,Indices(0,Dbond-1,0));
  site_odd.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  
  site_even=site_odd;
  
  site_last.setElement(2.0*ONE_c,Indices(0,0,0));
  site_last.setElement(ONE_c,Indices(Dbond-1,0,0));
  
  link.setElement(ONE_c,Indices(0,0,0));
  link.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  
  //sigma_z
  site_first.setElement(2.0*ONE_c*this->c_left,Indices(0,1,3));
  
  site_odd.setElement(-2.0*ONE_c,Indices(1,1,3));
  
  site_even.setElement(2.0*ONE_c,Indices(1,1,3));
  
  site_last.setElement(2.0*pow(-1.0,this->L-1)*ONE_c,Indices(1,0,3));
  
  link.setElement(ONE_c,Indices(0,1,3));
  link.setElement(ONE_c,Indices(1,Dbond-1,3)); 
  
  site_first.reshape(Indices(1*Dbond,7));
  site_odd.reshape(Indices(Dbond*Dbond,7));
  site_even.reshape(Indices(Dbond*Dbond,7));
  site_last.reshape(Indices(Dbond*1,7));
  link.reshape(Indices(Dbond*Dbond,7));  
  
  site_first.multiplyRight(reshape(Z,Indices(7,2*2)));
  site_odd.multiplyRight(reshape(Z,Indices(7,2*2)));
  site_even.multiplyRight(reshape(Z,Indices(7,2*2)));
  site_last.multiplyRight(reshape(Z,Indices(7,2*2)));
  link.multiplyRight(reshape(Z,Indices(7,2*2)));
  
  site_first.reshape(Indices(1,Dbond,2,2));
  site_odd.reshape(Indices(Dbond,Dbond,2,2));
  site_even.reshape(Indices(Dbond,Dbond,2,2));
  site_last.reshape(Indices(Dbond,1,2,2));
  link.reshape(Indices(Dbond,Dbond,2,2));  
  
  site_first.permute(Indices(3,1,4,2));
  site_odd.permute(Indices(3,1,4,2));
  site_even.permute(Indices(3,1,4,2));
  site_last.permute(Indices(3,1,4,2));
  link.permute(Indices(3,1,4,2));
  
  //Now fill the MPO
  mpo.setOp(0,new Operator(site_first),true);
  mpo.setOp(1,new Operator(link),true);
  mpo.setOp(2*this->L-2,new Operator(site_last),true);
  
  //Links first
  for(int i=3; i<2*this->L-2; i+=2)
    mpo.setOp(i,&mpo.getOp(1),false);
   
  //Sites in between the edges
  bool odd_saved=false, even_saved=false;
  int ind_odd_saved, ind_even_saved;
  for(int i=1; i<this->L-1; i++)
  {
   if((i%2)==0)
   {
     //Even site 
     if(even_saved)
       mpo.setOp(2*i,&mpo.getOp(ind_even_saved),false);
     else
     {
       mpo.setOp(2*i,new Operator(site_even),true);
       ind_even_saved=2*i;
       even_saved=true;
     }
   }
   else
   {
     //Even site 
     if(odd_saved)
       mpo.setOp(2*i,&mpo.getOp(ind_odd_saved),false);
     else
     {
       mpo.setOp(2*i,new Operator(site_odd),true);
       ind_odd_saved=2*i;
       odd_saved=true;
     }
   }
  }
  
#ifdef MATOUT
  mpo.exportForMatlab("Penalty.m",11);
#endif
  
}