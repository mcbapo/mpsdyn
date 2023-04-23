/**
   \file SU2HamiltonianNoLinksDimless.cpp
   Implementation of the truncated version of the SU(2) Hamiltonian from Hamer, C. Nuclear Physics B , 1982, 195, 503 - 521, this version is dimensionless and thus suitable for doing calculations with the aim to extrapolate to the continuum
   
   \author Stefan Kühn
   \date 23/10/2015
*/

#include "SU2HamiltonianNoLinksDimless.h"
#include "Indices.h"
#include "Contractor.h"
#include "splitLocalOp.h"
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>

using namespace std;
using namespace shrt;

//Constructor
SU2HamiltonianNoLinksDimless::SU2HamiltonianNoLinksDimless(int N_,double x_,double mg, int dlink_, double lambda_, double sz_penalty_)
{
  //Number of fermions
  this->N_fermions=N_;
  //Hopping constant
  this->epsilon=x_;  
  //Mass
  this->mu=2.0*mg*sqrt(this->epsilon);
  //Dimension of the link Hilbert spaces
  this->d_link=dlink_;
  //Dimension of the fermionic Hilbertspaces
  this->d_fermi=4;
  //Constant for the penalty term
  this->lambda=lambda_;
  //Target the sector with total charge equal to zero
  this->sz_penalty=sz_penalty_;
  //Initialize Hamiltonian
  this->hamil.initLength(this->N_fermions);
  
  //For the code to work, d_link>=2 is necessary (otherwise the electric field is zero everywhere and the model is anyway trivial)
  if(this->d_link<2)
  {
    cout << "Error in SU2HamiltonianNoLinks::SU2HamiltonianNoLinks(), link dimension has to be greater or equal to 2, instead a value of " << this->d_link << " was given" << endl;
    exit(666);
  }

  //Construct the basic operators
  initOperators();
  //Construct the actual hamiltonian MPO
  initHamiltonian();
}

//Destructor
SU2HamiltonianNoLinksDimless::~SU2HamiltonianNoLinksDimless(){}


//Build the Hamiltonian MPO
void SU2HamiltonianNoLinksDimless::initHamiltonian()
{
  int Dbond=this->d_link+16+1,Dl,Dr;
  
  cout << "Building Hamiltonian MPO with" << endl
       << "x:       " << this->epsilon << endl
       << "mu:      " << this->mu << endl
       << "dlink:   " << this->d_link << endl
       << "Dbond:   " << Dbond << endl;
  if(this->sz_penalty!=0.0)
    cout << "Imposing subsector with total charge zero with strength " << this->sz_penalty << endl;
  
  mwArray counting_tensor, mass_tensor, hopping_tensor, electric_tensor, identity_tensor;
  int offset=this->d_link-1,l1;
  
  /****************************************************************************************
   * 					First site
   * *************************************************************************************/
  counting_tensor=mwArray(Indices(this->d_fermi,1,this->d_fermi,Dbond));	counting_tensor.fillWithZero();
  mass_tensor=mwArray(Indices(this->d_fermi,1,this->d_fermi,Dbond));		mass_tensor.fillWithZero();
  hopping_tensor=mwArray(Indices(this->d_fermi,1,this->d_fermi,Dbond));		hopping_tensor.fillWithZero();
  electric_tensor=mwArray(Indices(this->d_fermi,1,this->d_fermi,Dbond));	electric_tensor.fillWithZero();
  
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
      electric_tensor.setElement(1.0/2.0*(1.0/2.0+1)*ONE_c,Indices(d,0,d,Dbond-1));
    }
    //Mass tensor
    mass_tensor.setElement(-1.0*this->mu*this->Noc.getElement(Indices(d,d)),Indices(d,0,d,Dbond-1));
  }
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
  this->hamil.setOp(0,new Operator(counting_tensor+mass_tensor+hopping_tensor+electric_tensor),true);
  
  /****************************************************************************************
   * 					Sites in between
   * *************************************************************************************/
  counting_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));
  mass_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));
  hopping_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));
  electric_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));
  identity_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));
  
  for(int i=1; i<(this->N_fermions-1); i++)
  {
    counting_tensor.fillWithZero();
    mass_tensor.fillWithZero();
    hopping_tensor.fillWithZero();
    electric_tensor.fillWithZero();
    identity_tensor.fillWithZero();
    
    //Count the value for the electric field
    for(int ind_left=0; ind_left<this->d_link; ind_left++)
    {
      for(int d=0; d<this->d_fermi; d++)
      {
	//Counting tensor & electric tensor
	if(d==0 || d==(this->d_fermi-1))
	{
	  counting_tensor.setElement(ONE_c,Indices(d,ind_left,d,ind_left));
	  electric_tensor.setElement(ind_left/2.0*(ind_left/2.0+1)*ONE_c,Indices(d,ind_left,d,Dbond-1));
	}
	else if(d==1 && ind_left-1>=0)
	{
	  counting_tensor.setElement(ONE_c,Indices(d,ind_left,d,ind_left-1));
	  electric_tensor.setElement((ind_left-1)/2.0*((ind_left-1)/2.0+1)*ONE_c,Indices(d,ind_left,d,Dbond-1));
	}
	else if(d==1)
	  counting_tensor.setElement(this->lambda*ONE_c,Indices(d,ind_left,d,Dbond-1)); //State with negative flux, put penalty and finish
	else if(d==2 && ind_left+1<this->d_link)
	{
	  counting_tensor.setElement(ONE_c,Indices(d,ind_left,d,ind_left+1));
	  electric_tensor.setElement((ind_left+1)/2.0*((ind_left+1)/2.0+1)*ONE_c,Indices(d,ind_left,d,Dbond-1));
	}
	else if(d==2)
	  counting_tensor.setElement(this->lambda*ONE_c,Indices(d,ind_left,d,Dbond-1));  //State with flux exceeding the maximum value, put penalty and finish
	
	//Mass tensor
	mass_tensor.setElement(pow(-1.0,i+1.0)*this->mu*this->Noc.getElement(Indices(d,d)),Indices(d,ind_left,d,Dbond-1));
	
	//Identities to the right
	identity_tensor.setElement(ONE_c,Indices(d,Dbond-1,d,Dbond-1));
	
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
      }
    }    
    
    //Set operator
    this->hamil.setOp(i,new Operator(counting_tensor+mass_tensor+hopping_tensor+electric_tensor+identity_tensor),true);
  }
  
  /****************************************************************************************
   * 					Last site
   * *************************************************************************************/
  counting_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,1));	counting_tensor.fillWithZero();
  mass_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,1));		mass_tensor.fillWithZero();
  hopping_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,1));		hopping_tensor.fillWithZero();
  electric_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,1));	electric_tensor.fillWithZero();
  identity_tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,1));	identity_tensor.fillWithZero();
  
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
    }
  }
  //Remove the state which would imply a negative flux for J_n (happens if J_n-1=0 and d=1), this entry could already be set from the penalty for targeting the total charge zero subsector, so I overwrite it here
  counting_tensor.setElement(this->lambda*ONE_c,Indices(1,0,1,0));
  
  //cout << "Coun;ing tensor last site" << endl;
  //print_tensor(counting_tensor);
  
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
  this->hamil.setOp(this->N_fermions-1,new Operator(counting_tensor+mass_tensor+hopping_tensor+electric_tensor+identity_tensor),true);
   
  //cout << "Hamiltonian " << this->hamil << endl; 
#ifdef MATOUT
  hamil.exportForMatlab("SU2HamiltonianNoLinksDimless.m");
#endif
}