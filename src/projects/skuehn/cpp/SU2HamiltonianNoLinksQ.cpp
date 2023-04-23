/**
   \file SU2HamiltonianNoLinksQ.cpp
   Implementation of the truncated version of the SU(2) Hamiltonian from Hamer, C. Nuclear Physics B , 1982, 195, 503 - 521
   
   \author Stefan Kühn
   \date 20/10/2015
*/

#include "SU2HamiltonianNoLinksQ.h"
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
SU2HamiltonianNoLinksQ::SU2HamiltonianNoLinksQ(int N_,double epsilon_,double mu_,double g_,int dlink_, double lambda_, double sz_penalty_, double charge_penalty_)
{
  ///Number of fermions
  this->N_fermions=N_;
  //Hopping constant
  this->epsilon=epsilon_;  
  //Mass
  this->mu=mu_;
  //Coupling
  this->g=g_;
  //Dimension of the link Hilbert spaces
  this->d_link=dlink_;
  //Dimension of the fermionic Hilbertspaces
  this->d_fermi=4;
  //Constant for the penalty term
  this->lambda=lambda_;
  //Target the sector with total charge equal to zero
  this->sz_penalty=sz_penalty_;
  //Set the charge penalty
  this->charge_penalty=charge_penalty_;
  
  //Set up Hamiltonian
  this->hamil.initLength(this->N_fermions);
  
  //For the code to work, d_link>=2 is necessary (otherwise the electric field is zero everywhere and the model is anyway trivial)
  if(this->d_link<2)
  {
    cout << "Error in SU2HamiltonianNoLinksQ::SU2HamiltonianNoLinksQ(), link dimension has to be greater or equal to 2, instead a value of " << this->d_link << " was given" << endl;
    exit(666);
  }

  //Construct the basic operators
  initOperators();
  //Construct the actual Hamiltonian MPO (now with the new routine overwriting the parent one)
  initHamiltonian();
}

//Destructor
SU2HamiltonianNoLinksQ::~SU2HamiltonianNoLinksQ(){}

//Build the Hamiltonian MPO
void SU2HamiltonianNoLinksQ::initHamiltonian()
{
  //Compared to the model without charge penalty, my bond dimension is increased by one
  int Dbond=this->d_link+16+1+1+1,Dl,Dr;
  
  cout << "Building Hamiltonian MPO containing charge penalty with" << endl
       << "epsilon: " << this->epsilon << endl
       << "mu:      " << this->mu << endl
       << "g:       " << this->g << endl
       << "dlink:   " << this->d_link << endl
       << "Dbond:   " << Dbond << endl;
  if(this->sz_penalty!=0.0)
    cout << "Imposing subsector with total charge zero with strength " << this->sz_penalty << endl;
  if(this->charge_penalty!=0.0)
    cout << "Imposing balance between particles and antiparticles with strength " << this->sz_penalty << endl;
  
  //Two additional operators I need for the charge penalty
  mwArray OpOdd,OpEven;
  OpOdd=mwArray(Indices(this->d_fermi,this->d_fermi));		OpOdd.fillWithZero();
  OpEven=mwArray(Indices(this->d_fermi,this->d_fermi));		OpEven.fillWithZero();
  OpOdd = this->Noc-2.0*identityMatrix(this->d_fermi);
  OpEven = this->Noc;
  
  /*cout << "OpEven = " << OpEven << endl;
  cout << "OpOdd = " << OpOdd << endl;*/
  
  mwArray counting_tensor, mass_tensor, hopping_tensor, electric_tensor, identity_tensor, charge_penalty_tensor;
  int offset=this->d_link-1,l1;
  
  //cout << endl << "16+offset=" << 16+offset << endl;
  
  /****************************************************************************************
   * 					First site
   * *************************************************************************************/
  counting_tensor=mwArray(Indices(this->d_fermi,1,this->d_fermi,Dbond));	counting_tensor.fillWithZero();
  mass_tensor=mwArray(Indices(this->d_fermi,1,this->d_fermi,Dbond));		mass_tensor.fillWithZero();
  hopping_tensor=mwArray(Indices(this->d_fermi,1,this->d_fermi,Dbond));		hopping_tensor.fillWithZero();
  electric_tensor=mwArray(Indices(this->d_fermi,1,this->d_fermi,Dbond));	electric_tensor.fillWithZero();
  charge_penalty_tensor=mwArray(Indices(this->d_fermi,1,this->d_fermi,Dbond));	charge_penalty_tensor.fillWithZero();
  identity_tensor=mwArray(Indices(this->d_fermi,1,this->d_fermi,Dbond));	identity_tensor.fillWithZero();
  
  for(int d=0; d<this->d_fermi; d++)
  {
    //Counting tensor and electric tensor
    if(d==0 || d==(this->d_fermi-1))
      counting_tensor.setElement(ONE_c,Indices(d,0,d,0));
    else if(d==1)
      counting_tensor.setElement(this->lambda*ONE_c,Indices(d,0,d,Dbond-1));
    else if(d==2)
    {
      counting_tensor.setElement(ONE_c,Indices(d,0,d,1));
      electric_tensor.setElement(this->g*this->g/2.0*1.0/2.0*(1.0/2.0+1.0)*ONE_c,Indices(d,0,d,Dbond-1));
    }
    //Mass tensor
    mass_tensor.setElement(-1.0*this->mu*this->Noc.getElement(Indices(d,d)),Indices(d,0,d,Dbond-1));
    
    //Charge penalty (odd site, as first site is always odd)
    //cout << "Taking odd for site 1" << endl;
    charge_penalty_tensor.setElement(this->charge_penalty*OpOdd.getElement(Indices(d,d))*OpOdd.getElement(Indices(d,d)),Indices(d,0,d,Dbond-1)); //Set the Q² and finish
    charge_penalty_tensor.setElement(2.0*this->charge_penalty*OpOdd.getElement(Indices(d,d)),Indices(d,0,d,Dbond-2)); //Start long range two body term
    
    //Identity tensor
    identity_tensor.setElement(ONE_c,Indices(d,0,d,Dbond-3));
  }
  //cout << "Penalty tensor: " << endl;
  //print_tensor(charge_penalty_tensor,false);
  //cout << "Counting tensor: " << endl;
  //print_tensor(counting_tensor,false);
  //Input is simply zero flux
  l1=0;
  //O1
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1-1,0)*ONE_c,Indices(0,l1,1,1+offset));
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1-1,1)*ONE_c,Indices(0,l1,1,9+offset));
  //O2
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1+1,0)*ONE_c,Indices(0,l1,2,2+offset));
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1+1,1)*ONE_c,Indices(0,l1,2,10+offset));
  //O3
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1-1,0)*ONE_c,Indices(1,l1,3,3+offset));
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1-1,2)*ONE_c,Indices(1,l1,3,11+offset));
  //O4
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1+1,0)*ONE_c,Indices(2,l1,3,4+offset));
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1+1,2)*ONE_c,Indices(2,l1,3,12+offset));
  //O1^\dagger
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1-1,0)*ONE_c,Indices(1,l1,0,5+offset));
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1-1,1)*ONE_c,Indices(1,l1,0,13+offset));
  //O2^\dagger
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1+1,0)*ONE_c,Indices(2,l1,0,6+offset));
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1+1,1)*ONE_c,Indices(2,l1,0,14+offset));
  //O3^\dagger
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1-1,0)*ONE_c,Indices(3,l1,1,7+offset));
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1-1,2)*ONE_c,Indices(3,l1,1,15+offset));
  //O4^\dagger
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1+1,0)*ONE_c,Indices(3,l1,2,8+offset));
  hopping_tensor.setElement(this->epsilon*get_coefficient(l1,l1+1,2)*ONE_c,Indices(3,l1,2,16+offset));

  
  //Set operator
  this->hamil.setOp(0,new Operator(counting_tensor+mass_tensor+hopping_tensor+electric_tensor+charge_penalty_tensor+identity_tensor),true);  
  
  /****************************************************************************************
   * 					Sites in between
   * *************************************************************************************/
  counting_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));
  mass_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));
  hopping_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));
  electric_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));
  identity_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));
  charge_penalty_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));
  
  for(int i=1; i<(this->N_fermions-1); i++)
  {
    counting_tensor.fillWithZero();
    mass_tensor.fillWithZero();
    hopping_tensor.fillWithZero();
    electric_tensor.fillWithZero();
    identity_tensor.fillWithZero();
    charge_penalty_tensor.fillWithZero();
    
    //Count the value for the electric field
    for(int ind_left=0; ind_left<this->d_link; ind_left++)
    {
      for(int d=0; d<this->d_fermi; d++)
      {
	//Counting tensor & electric tensor
	if(d==0 || d==(this->d_fermi-1))
	{
	  counting_tensor.setElement(ONE_c,Indices(d,ind_left,d,ind_left));
	  electric_tensor.setElement(this->g*this->g/2.0*ind_left/2.0*(ind_left/2.0+1)*ONE_c,Indices(d,ind_left,d,Dbond-1));
	}
	else if(d==1 && ind_left-1>=0)
	{
	  counting_tensor.setElement(ONE_c,Indices(d,ind_left,d,ind_left-1));
	  electric_tensor.setElement(this->g*this->g/2.0*(ind_left-1)/2.0*((ind_left-1)/2.0+1)*ONE_c,Indices(d,ind_left,d,Dbond-1));
	}
	else if(d==1)
	  counting_tensor.setElement(this->lambda*ONE_c,Indices(d,ind_left,d,Dbond-1)); //State with negative flux, put penalty and finish
	else if(d==2 && ind_left+1<this->d_link)
	{
	  counting_tensor.setElement(ONE_c,Indices(d,ind_left,d,ind_left+1));
	  electric_tensor.setElement(this->g*this->g/2.0*(ind_left+1)/2.0*((ind_left+1)/2.0+1)*ONE_c,Indices(d,ind_left,d,Dbond-1));
	}
	else if(d==2)
	  counting_tensor.setElement(this->lambda*ONE_c,Indices(d,ind_left,d,Dbond-1));  //State with flux exceeding the maximum value, put penalty and finish
	
	//Mass tensor
	mass_tensor.setElement(pow(-1.0,i+1.0)*this->mu*this->Noc.getElement(Indices(d,d)),Indices(d,ind_left,d,Dbond-1));
	
	//Identities to the right
	identity_tensor.setElement(ONE_c,Indices(d,Dbond-1,d,Dbond-1));
	//Identities between the long range two body term
	identity_tensor.setElement(ONE_c,Indices(d,Dbond-2,d,Dbond-2));
	identity_tensor.setElement(ONE_c,Indices(d,Dbond-3,d,Dbond-3));	
	
	//O1
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left-1,0)*ONE_c,Indices(0,ind_left,1,1+offset));
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left-1,1)*ONE_c,Indices(0,ind_left,1,9+offset));
	hopping_tensor.setElement(ONE_c,Indices(0,8+offset,1,Dbond-1));
	hopping_tensor.setElement(ONE_c,Indices(0,13+offset,1,Dbond-1));
	//O2
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left+1,0)*ONE_c,Indices(0,ind_left,2,2+offset));
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left+1,1)*ONE_c,Indices(0,ind_left,2,10+offset));
	hopping_tensor.setElement(ONE_c,Indices(0,7+offset,2,Dbond-1));
	hopping_tensor.setElement(ONE_c,Indices(0,14+offset,2,Dbond-1));
	//O3
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left-1,0)*ONE_c,Indices(1,ind_left,3,3+offset));
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left-1,2)*ONE_c,Indices(1,ind_left,3,11+offset));
	hopping_tensor.setElement(ONE_c,Indices(1,6+offset,3,Dbond-1));
	hopping_tensor.setElement(ONE_c,Indices(1,15+offset,3,Dbond-1));
	//O4
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left+1,0)*ONE_c,Indices(2,ind_left,3,4+offset));
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left+1,2)*ONE_c,Indices(2,ind_left,3,12+offset));
	hopping_tensor.setElement(ONE_c,Indices(2,5+offset,3,Dbond-1));
	hopping_tensor.setElement(ONE_c,Indices(2,16+offset,3,Dbond-1));
	//O1^\dagger
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left-1,0)*ONE_c,Indices(1,ind_left,0,5+offset));
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left-1,1)*ONE_c,Indices(1,ind_left,0,13+offset));
	hopping_tensor.setElement(ONE_c,Indices(1,4+offset,0,Dbond-1));
	hopping_tensor.setElement(ONE_c,Indices(1,9+offset,0,Dbond-1));
	//O2^\dagger
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left+1,0)*ONE_c,Indices(2,ind_left,0,6+offset));
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left+1,1)*ONE_c,Indices(2,ind_left,0,14+offset));
	hopping_tensor.setElement(ONE_c,Indices(2,3+offset,0,Dbond-1));
	hopping_tensor.setElement(ONE_c,Indices(2,10+offset,0,Dbond-1));
	//O3^\dagger
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left-1,0)*ONE_c,Indices(3,ind_left,1,7+offset));
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left-1,2)*ONE_c,Indices(3,ind_left,1,15+offset));
	hopping_tensor.setElement(ONE_c,Indices(3,2+offset,1,Dbond-1));
	hopping_tensor.setElement(ONE_c,Indices(3,11+offset,1,Dbond-1));
	//O4^\dagger
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left+1,0)*ONE_c,Indices(3,ind_left,2,8+offset));
	hopping_tensor.setElement(this->epsilon*get_coefficient(ind_left,ind_left+1,2)*ONE_c,Indices(3,ind_left,2,16+offset));
	hopping_tensor.setElement(ONE_c,Indices(3,1+offset,2,Dbond-1));
	hopping_tensor.setElement(ONE_c,Indices(3,12+offset,2,Dbond-1));
	
	//Charge penalty tensor
	if((i%2)==1)
	{
	  //cout << "Taking even for site "<< i+1 << endl;
	  //Even site (take into account that I start counting my indices by zero but sites range form 1 to this->N_fermions)
	  charge_penalty_tensor.setElement(this->charge_penalty*OpEven.getElement(Indices(d,d))*OpEven.getElement(Indices(d,d)),Indices(d,Dbond-3,d,Dbond-1));
	  charge_penalty_tensor.setElement(2.0*this->charge_penalty*OpEven.getElement(Indices(d,d)),Indices(d,Dbond-3,d,Dbond-2));
	  charge_penalty_tensor.setElement(OpEven.getElement(Indices(d,d)),Indices(d,Dbond-2,d,Dbond-1));
	}
	else
	{
	  //cout << "Taking odd for site " << i+1 << endl;
	  //Odd site
	  charge_penalty_tensor.setElement(this->charge_penalty*OpOdd.getElement(Indices(d,d))*OpOdd.getElement(Indices(d,d)),Indices(d,Dbond-3,d,Dbond-1));
	  charge_penalty_tensor.setElement(2.0*this->charge_penalty*OpOdd.getElement(Indices(d,d)),Indices(d,Dbond-3,d,Dbond-2));
	  charge_penalty_tensor.setElement(OpOdd.getElement(Indices(d,d)),Indices(d,Dbond-2,d,Dbond-1));
	}
      }
    } 
    //cout << "Penalty tensor: " << endl;
    //print_tensor(charge_penalty_tensor,false);
    //cout << "Counting tensor: " << endl;
    //print_tensor(counting_tensor,false);
    //Set operator
    this->hamil.setOp(i,new Operator(counting_tensor+mass_tensor+hopping_tensor+electric_tensor+identity_tensor+charge_penalty_tensor),true);
  }
  
  /****************************************************************************************
   * 					Last site
   * *************************************************************************************/
  counting_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,1));	counting_tensor.fillWithZero();
  mass_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,1));		mass_tensor.fillWithZero();
  hopping_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,1));		hopping_tensor.fillWithZero();
  electric_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,1));	electric_tensor.fillWithZero();
  identity_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,1));	identity_tensor.fillWithZero();
  charge_penalty_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,1));	charge_penalty_tensor.fillWithZero();
  
  for(int ind_left=0; ind_left<this->d_link; ind_left++)
  {
    for(int d=0; d<this->d_fermi; d++)
    { 
      //Mass tensor
      mass_tensor.setElement(pow(-1.0,this->d_fermi)*this->mu*this->Noc.getElement(Indices(d,d)),Indices(d,ind_left,d,0));
      
      //Identities to the right
      identity_tensor.setElement(ONE_c,Indices(d,Dbond-1,d,0));
      
      //If I want to project on the subspace with total charge equal to zero, I do it with the counting tensor
      //\TODO Optimize a bit, it is anyway filled with zeros, so I don't have to care about that
      if(d==1 && ind_left==1)
      {
	//cout << "Having case that d=1 and ind_left=1" << endl;
	counting_tensor.setElement(ZERO_c,Indices(d,ind_left,d,0));
      }
      else if((d==0 || d==(this->d_fermi-1)) && ind_left==0)
      {
	//cout << "Having case that d=0 or d=3 and ind_left=0" << endl;
	counting_tensor.setElement(ZERO_c,Indices(d,ind_left,d,0));
      }
      else
      {
	//cout << "Having case which should be projected out" << endl;
	counting_tensor.setElement(this->sz_penalty*ONE_c,Indices(d,ind_left,d,0));	//Penalize all states that lead to a non vanishing electric flux on site n
      }
      
      //Charge charge penalty tensor
      if((this->N_fermions%2)==0)
      {
	//cout << "Taking even for site "<< this->N_fermions << endl;
	//Even site
	charge_penalty_tensor.setElement(this->charge_penalty*OpEven.getElement(Indices(d,d))*OpEven.getElement(Indices(d,d)),Indices(d,Dbond-3,d,0));
	charge_penalty_tensor.setElement(OpEven.getElement(Indices(d,d)),Indices(d,Dbond-2,d,0));
      }
      else
      {
	//cout << "Taking odd for site "<< this->N_fermions << endl;
	//Odd site
	charge_penalty_tensor.setElement(this->charge_penalty*OpOdd.getElement(Indices(d,d))*OpOdd.getElement(Indices(d,d)),Indices(d,Dbond-3,d,0));
	charge_penalty_tensor.setElement(OpOdd.getElement(Indices(d,d)),Indices(d,Dbond-2,d,0));
      }
    }
  }
  //Remove the state which would imply a negative flux for J_n (happens if J_n-1=0 and d=1), this entry could already be set from the penalty for targeting the total charge zero subsector, so I overwrite it here
  counting_tensor.setElement(this->lambda*ONE_c,Indices(1,0,1,0));
  
  //cout << "Counting tensor last site" << endl;
  ////print_tensor(counting_tensor);
  
  //cout << "Penalty tensor: " << endl;
  //print_tensor(charge_penalty_tensor,false);
  //cout << "Counting tensor: " << endl;
  //print_tensor(counting_tensor,false);
  
  //O1
  hopping_tensor.setElement(ONE_c,Indices(0,8+offset,1,0));
  hopping_tensor.setElement(ONE_c,Indices(0,13+offset,1,0));
  //O2
  hopping_tensor.setElement(ONE_c,Indices(0,7+offset,2,0));
  hopping_tensor.setElement(ONE_c,Indices(0,14+offset,2,0));
  //O3
  hopping_tensor.setElement(ONE_c,Indices(1,6+offset,3,0));
  hopping_tensor.setElement(ONE_c,Indices(1,15+offset,3,0));
  //O4
  hopping_tensor.setElement(ONE_c,Indices(2,5+offset,3,0));
  hopping_tensor.setElement(ONE_c,Indices(2,16+offset,3,0));
  //O1^\dagger
  hopping_tensor.setElement(ONE_c,Indices(1,4+offset,0,0));
  hopping_tensor.setElement(ONE_c,Indices(1,9+offset,0,0));
  //O2^\dagger
  hopping_tensor.setElement(ONE_c,Indices(2,3+offset,0,0));
  hopping_tensor.setElement(ONE_c,Indices(2,10+offset,0,0));
  //O3^\dagger
  hopping_tensor.setElement(ONE_c,Indices(3,2+offset,1,0));
  hopping_tensor.setElement(ONE_c,Indices(3,11+offset,1,0));
  //O4^\dagger
  hopping_tensor.setElement(ONE_c,Indices(3,1+offset,2,0));
  hopping_tensor.setElement(ONE_c,Indices(3,12+offset,2,0)); 
  
  //Set operator
  this->hamil.setOp(this->N_fermions-1,new Operator(counting_tensor+mass_tensor+hopping_tensor+electric_tensor+identity_tensor+charge_penalty_tensor),true);
  
  //cout << "Hamiltonian " << this->hamil << endl; 
#ifdef MATOUT
  hamil.exportForMatlab("SU2HamiltonianNoLinksQ.m");
#endif
}
