/**
   \file SU2HamiltonianNoLinks.cpp
   Implementation of the truncated version of the SU(2) Hamiltonian from Hamer, C. Nuclear Physics B , 1982, 195, 503 - 521
   
   \author Stefan Kühn
   \date 09/10/2015
*/

#include "SU2HamiltonianNoLinks.h"
#include "Indices.h"
#include "Contractor.h"
#include "splitLocalOp.h"
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>

using namespace std;
using namespace shrt;

//Empty constructor for deriving classes
SU2HamiltonianNoLinks::SU2HamiltonianNoLinks():hamil(1){}


//Standard Constructor
SU2HamiltonianNoLinks::SU2HamiltonianNoLinks(int N_,double epsilon_,double mu_,double g_,int dlink_, double lambda_, double sz_penalty_):hamil(N_)
{
  //Number of fermions
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
SU2HamiltonianNoLinks::~SU2HamiltonianNoLinks(){
  //Clear Hamiltonian MPO
  hamil.clear();  
  //Clear other operators?
}

//Provide the basic operator matrices used for MPOs
void SU2HamiltonianNoLinks::initOperators(void)
{    

  //Fermionic operators
  this->O1=mwArray(Indices(this->d_fermi,this->d_fermi));	this->O1.fillWithZero();	
  this->O1dagger=mwArray(Indices(this->d_fermi,this->d_fermi));	this->O1dagger.fillWithZero();
  this->O2=mwArray(Indices(this->d_fermi,this->d_fermi));	this->O2.fillWithZero();	
  this->O2dagger=mwArray(Indices(this->d_fermi,this->d_fermi));	this->O2dagger.fillWithZero();
  this->O3=mwArray(Indices(this->d_fermi,this->d_fermi));	this->O3.fillWithZero();	
  this->O3dagger=mwArray(Indices(this->d_fermi,this->d_fermi));	this->O3dagger.fillWithZero();
  this->O4=mwArray(Indices(this->d_fermi,this->d_fermi));	this->O4.fillWithZero();	
  this->O4dagger=mwArray(Indices(this->d_fermi,this->d_fermi));	this->O4dagger.fillWithZero();
  
  this->Noc=mwArray(Indices(this->d_fermi,this->d_fermi));	this->Noc.fillWithZero();
  
  this->O1.setElement(ONE_c,Indices(0,1));	this->O1dagger=this->O1; 			this->O1dagger.transpose(true);
  this->O2.setElement(ONE_c,Indices(0,2));	this->O2dagger=this->O2; 			this->O2dagger.transpose(true);
  this->O3.setElement(ONE_c,Indices(1,3));	this->O3dagger=this->O3; 			this->O3dagger.transpose(true);
  this->O4.setElement(ONE_c,Indices(2,3));	this->O4dagger=this->O4; 			this->O4dagger.transpose(true);
  this->Noc.setElement(ONE_c,Indices(1,1));	this->Noc.setElement(ONE_c,Indices(2,2));	this->Noc.setElement(2.0*ONE_c,Indices(3,3));  
  
  //Identity for fermionic site
  this->id_fermi=identityMatrix(this->d_fermi);  
  
  //I often need identities to construct a single body operator, therefore I prepare those and save them
  this->id_fermi_mpo=this->id_fermi;
  this->id_fermi_mpo.reshape(Indices(this->d_fermi,1,this->d_fermi,1));    
  
  //For debugging
  /*cout << "O1      = " << this->O1 << endl;
  cout << "O2      = " << this->O2 << endl;
  cout << "O3      = " << this->O3 << endl;
  cout << "O4      = " << this->O4 << endl;
  cout << "O1d     = " << this->O1dagger << endl;
  cout << "O2d     = " << this->O2dagger << endl;
  cout << "O3d     = " << this->O3dagger << endl;
  cout << "O4d     = " << this->O4dagger << endl;
  cout << "Noc     = " << this->Noc << endl;*/
}


void SU2HamiltonianNoLinks::print_tensor(const mwArray &tensor,bool fix_physical) const
{
  int dl,dr;
  //Check the dimensions
  if(fix_physical)
  {
    if(tensor.getRank()==3)
    {
      dl=tensor.getDimension(0);
      for(int i=0; i<dl; i++)
	if(!tensor.subArray(Indices(i,-1,-1)).isNull(1.0E-6))
	  cout << "Array for dl=" << i << ": " << tensor.subArray(Indices(i,-1,-1)) << endl;
    }
    else if(tensor.getRank()==4)
    {
      dl=tensor.getDimension(0);
      dr=tensor.getDimension(2);
      for(int i=0; i<dl; i++)
	for(int j=0; j<dr; j++)
	{
	  if(!tensor.subArray(Indices(i,-1,j,-1)).isNull(1.0E-6))
	    cout << "Array for (dl,dr)=(" << i << "," << j << "): " << tensor.subArray(Indices(i,-1,j,-1)) << endl;
	}
    }
    else
      cout << "Only tensors of rank 3 and 4 are supported, given tensor has rank " << tensor.getRank() << endl;
  }
  else
  {
    if(tensor.getRank()==3)
    {
      dl=tensor.getDimension(2);
      dr=tensor.getDimension(3);
      for(int i=0; i<dl; i++)
	for(int j=0; j<dr; j++)
	{
	  if(!tensor.subArray(Indices(-1,i,j)).isNull(1.0E-6))
	    cout << "Array for (Dl,Dr)=(" << i << "," << j << "): " << tensor.subArray(Indices(-1,i,j)) << endl;
	}
    }
    else if(tensor.getRank()==4)
    {
      dl=tensor.getDimension(1);
      dr=tensor.getDimension(3);
      for(int i=0; i<dl; i++)
	for(int j=0; j<dr; j++)
	{
	  if(!tensor.subArray(Indices(-1,i,-1,j)).isNull(1.0E-6))
	    cout << "Array for (Dl,Dr)=(" << i << "," << j << "): " << tensor.subArray(Indices(-1,i,-1,j)) << endl;
	}
    }
    else
      cout << "Only tensors of rank 3 and 4 are supported, given tensor has rank " << tensor.getRank() << endl;
  }
}

//Little function to compute the coefficients
double SU2HamiltonianNoLinks::get_coefficient(int l1, int l2, int flag) const
{
  //First check, if I obtained a case which is physically impossible, because the flux would be negative, in that case, I simply return 0 meaning that the corresponding operator is simply a zero
  double coefficient;
  if(l2<0 || l2>= this->d_link)
    coefficient=0.0;
  else if(flag==0)
  {    
    //Convert l1 and l2 to flux values
    double flux1=((double)l1)/2.0;
    double flux2=((double)l2)/2.0;
    coefficient = pow(-1.0,flux2-flux1-0.5)*sqrt((2.0*flux2+1.0)/(2.0*flux1+1.0));
  }
  else if(flag==1)
    coefficient = 1.0;
  else
    coefficient = -1.0;
  /*if(flag=0 &&coefficient<=0.0)
    cout << "Got value " << coefficient << " for coefficient "<< endl;
  if(coefficient>=999.0)
    cout << "Penalizing transition for l1=" << l1 << ", l2=" << l2 << ", c=" << coefficient << endl;*/
  return coefficient;    
}


//Build the Hamiltonian MPO
void SU2HamiltonianNoLinks::initHamiltonian()
{
  int Dbond=this->d_link+16+1,Dl,Dr;
  
  cout << "Building Hamiltonian MPO with" << endl
       << "epsilon: " << this->epsilon << endl
       << "mu:      " << this->mu << endl
       << "g:       " << this->g << endl
       << "dlink:   " << this->d_link << endl
       << "Dbond:   " << Dbond << endl;
  if(this->sz_penalty!=0.0)
    cout << "Imposing subsector with total charge zero with strength " << this->sz_penalty << endl;
  
  mwArray tensor;
  
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
      electric_tensor.setElement(this->g*this->g/2.0*1.0/2.0*(1.0/2.0+1)*ONE_c,Indices(d,0,d,Dbond-1));
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
  hamil.exportForMatlab("SU2HamiltonianNoLinks.m");
#endif
}

void SU2HamiltonianNoLinks::constructProductStateMPS(MPS &mps,vector<int> sites) const
{
  mps.clear();
  vector<int> dims(this->N_fermions,this->d_fermi);
  mps=MPS(this->N_fermions,1,dims);  
  
  mwArray fermi(Indices(this->d_fermi,1,1));
  if(sites.empty())
  {
    for(int i=0; i<this->N_fermions; i++)
    {
      fermi.fillWithZero();
      if((i%2)==0)
	fermi.setElement(ONE_c,Indices(this->d_fermi-1,0,0));
      else
	fermi.setElement(ONE_c,Indices(0,0,0));
      mps.setA(i,fermi);
    }
  }
  else
  {
    //Check if given values are valid
    int min = *min_element(sites.begin(), sites.end());
    int max = *max_element(sites.begin(), sites.end());
    if(min<0 || max>=this->d_fermi)
    {
      cout << "Error in SU2HamiltonianNoLinks::constructProductStateMPS(), invalid range of values" << endl;
      exit(666);
    }
    //Check if enough values are given
    if(sites.size()!=this->N_fermions)
    {
      cout << "Error in SU2HamiltonianNoLinks::constructProductStateMPS(), number of sites specified is not compliant with existing number of sites" << endl;
      exit(666);
    }
    //Another little check, if the specified state is gauge invariant and compliant with the given maximum flux. As I rely on the fact, that I only deal with valid states for my Implementation to work, I give a warning if this is not the case
    int l_left=0,l_right;
    for (int i=0; i<sites.size()-1; i++)
    {
      if(sites[i]==0 || sites[i]==(this->d_fermi-1))
	l_right=l_left;
      else if(sites[i]==1)
	l_right=l_left-1;
      else if(sites[i]==2)
	l_right=l_left+1;
      if(l_right>=this->d_link)
	cout << "Warning in SU2HamiltonianNoLinks::constructProductStateMPS(), electric field on site " << i+1 << " has value " << l_right/2.0 << " and exceeds maximum possible value " << (this->d_link-1.0)/2.0 << ", state might cause undefined behaviour" << endl;
      else if(l_right<0)
	cout << "Warning in SU2HamiltonianNoLinks::constructProductStateMPS(), electric field on site " << i+1 << " has value " << l_right/2.0 << " and exceeds minimum possible value 0, state might cause undefined behaviour" << endl;
      l_left=l_right;
    }
    for(int i=0; i<sites.size(); i++)
    {
      fermi.fillWithZero();
      fermi.setElement(ONE_c,Indices(sites[i],0,0));
      mps.setA(i,fermi);
    }
  }
  //Normalize and gauge
  mps.gaugeCond('R',true);  
}

void SU2HamiltonianNoLinks::get_J_Jsquared_MPO(MPO &JsquaredMPO, int index, bool squared) const
{
  if(index<1 || index>this->N_fermions)
  {
   cout << "Error in SU2HamiltonianNoLinks::getJsquaredMPO(), given site index " << index << " is out of range" << endl; 
   exit(666);
  }
  JsquaredMPO.initLength(this->N_fermions);
  
  cout << "--> Creating MPO for site " << index << endl;
  
  mwArray tensor;
  int Dl,Dr;
  int ind_left,ind_right;
  complex_t val;
  
  //Count until one reaches the site in question
  for(int i=0; i<index-1; i++)
  {
    Dr = (i+2<=this->d_link) ? i+2:this->d_link;
    Dl = (i+1<=this->d_link) ? i+1:this->d_link;
    tensor = mwArray(Indices(this->d_fermi,Dl,this->d_fermi,Dr));
    tensor.fillWithZero();
    for(int ind_left=0; ind_left<Dl; ind_left++)
    {
      for(int d=0; d<this->d_fermi; d++)
      {
	if(d==0 || d==(this->d_fermi-1))
	  tensor.setElement(ONE_c,Indices(d,ind_left,d,ind_left));
	else if(d==1 && ind_left-1>=0)
	  tensor.setElement(ONE_c,Indices(d,ind_left,d,ind_left-1));
	else if(d==2 && ind_left+1<Dr)
	  tensor.setElement(ONE_c,Indices(d,ind_left,d,ind_left+1));
      }
    }
    JsquaredMPO.setOp(i,new Operator(tensor),true);
  }
  //Site in question
  Dl=(index<=this->d_link) ? index:this->d_link;
  tensor = mwArray(Indices(this->d_fermi,Dl,this->d_fermi,1));
  tensor.fillWithZero();
  for(int ind_left=0; ind_left<Dl; ind_left++)
  {
    for(int d=0; d<this->d_fermi; d++)
    {
      if(d==0 || d==(this->d_fermi-1))
      {
	val = squared ? ind_left/2.0*(ind_left/2.0+1.0)*ONE_c : ind_left/2.0*ONE_c;
	tensor.setElement(val,Indices(d,ind_left,d,0));
      }
      else if(d==1 && ind_left-1>=0)
      {
	val = squared ? (ind_left-1.0)/2.0*((ind_left-1.0)/2.0+1.0)*ONE_c : (ind_left-1.0)/2.0*ONE_c;
	tensor.setElement(val,Indices(d,ind_left,d,0));
      }
      else if(d==2 && ind_left+1<this->d_link)
      {
	val = squared ? (ind_left+1.0)/2.0*((ind_left+1.0)/2.0+1.0)*ONE_c : (ind_left+1.0)/2.0*ONE_c;
	tensor.setElement(val,Indices(d,ind_left,d,0));
      }
    }
  }
  JsquaredMPO.setOp(index-1,new Operator(tensor),true);
  
  //Simply identity to the right of the site
  bool id_saved=false;
  int pos_saved;
  for(int i=index; i<this->N_fermions; i++)
  {
    if(id_saved)
      JsquaredMPO.setOp(i,&JsquaredMPO.getOp(pos_saved),false);
    else
    {
      JsquaredMPO.setOp(i,new Operator(this->id_fermi_mpo),true);
      id_saved=true;
      pos_saved=i;
    }
  }  
}

void SU2HamiltonianNoLinks::getJMPO(MPO &JMPO, int index) const
{
  this->get_J_Jsquared_MPO(JMPO, index, false);
}

void SU2HamiltonianNoLinks::getJsquaredMPO(MPO &JsquaredMPO, int index) const
{
  this->get_J_Jsquared_MPO(JsquaredMPO, index, true);
}

int SU2HamiltonianNoLinks::getSingleBodyMPO(MPO &mpo, int index, const mwArray &Op) const
{
  //Some error checking
  if(index<1 || index>this->N_fermions)
    return -1;
  //Frist clear the MPO and resize it
  mpo.clear();
  mpo.initLength(this->N_fermions);  
  
  //Now set operators for fermionic sites
  bool fermi_saved=false;
  int saved_site;
  for(int k=1; k<=this->N_fermions; k++)
  {
    if(fermi_saved && (k!=index))
      mpo.setOp(k-1,&mpo.getOp(saved_site),false);
    else
    {
      if(k==index)
	mpo.setOp(k-1,new Operator(reshape(Op,Indices(this->d_fermi,1,this->d_fermi,1))),true);
      else
      {
	mpo.setOp(k-1,new Operator(this->id_fermi_mpo),true);
	fermi_saved=true;
	saved_site=k-1;
      }
    }
  }
  return 0;
}

void SU2HamiltonianNoLinks::getNocMPO(MPO &mpo, int index) const
{
  int status;
  status = this->getSingleBodyMPO(mpo,index,this->Noc);
  if(status<0)
  {
    cout << "Error in SU2HamiltonianNoLinks::getNocMPO()" << endl;
    exit(666);
  }  
}

void SU2HamiltonianNoLinks::getChargeSquareMPO(MPO &mpo, int pos) const
{
  int status;
  mwArray Operator(Indices(this->d_fermi,this->d_fermi));
  Operator.fillWithZero();
  Operator.setElement(0.25*ONE_c,Indices(1,1));
  Operator.setElement(0.25*ONE_c,Indices(2,2));
  status = this->getSingleBodyMPO(mpo,pos,Operator);
  if(status<0)
  {
    cout << "Error in SU2HamiltonianNoLinks::getChargeSquareMPO()" << endl;
    exit(666);
  }
}

void SU2HamiltonianNoLinks::getNocMPO(MPO &mpo) const
{
  mpo.initLength(this->N_fermions);
  
  mwArray Operators(Indices(2,this->d_fermi,this->d_fermi));
  Operators.fillWithZero();
  
  for(int i=0; i<this->d_fermi; i++)
  {
    Operators.setElement(this->id_fermi.getElement(Indices(i,i)),Indices(0,i,i));
    Operators.setElement(this->Noc.getElement(Indices(i,i)),Indices(1,i,i));
  }
  Operators.reshape(Indices(2,this->d_fermi*this->d_fermi));
  
  int Dbond=2;
  mwArray Left,Middle,Right;
  Left=mwArray(Indices(1,Dbond,2));		Left.fillWithZero();
  Middle=mwArray(Indices(Dbond,Dbond,2));	Middle.fillWithZero();
  Right=mwArray(Indices(Dbond,1,2));		Right.fillWithZero();
  
  Left.setElement(ONE_c,Indices(0,0,0));
  Left.setElement(ONE_c,Indices(0,1,1));
  
  Middle.setElement(ONE_c,Indices(0,0,0));
  Middle.setElement(ONE_c,Indices(1,1,0));
  Middle.setElement(ONE_c,Indices(0,1,1));
  
  Right.setElement(ONE_c,Indices(1,0,0));
  Right.setElement(ONE_c,Indices(0,0,1));
  
  
  Left.reshape(Indices(1*Dbond,2));		Left.multiplyRight(Operators);		Left.reshape(Indices(1,Dbond,this->d_fermi,this->d_fermi));		Left.permute(Indices(3,1,4,2));
  Middle.reshape(Indices(Dbond*Dbond,2));	Middle.multiplyRight(Operators);	Middle.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Middle.permute(Indices(3,1,4,2));
  Right.reshape(Indices(Dbond*1,2));		Right.multiplyRight(Operators);		Right.reshape(Indices(Dbond,1,this->d_fermi,this->d_fermi));		Right.permute(Indices(3,1,4,2));
  
  mpo.setOp(0,new Operator(Left),true);
  mpo.setOp(this->N_fermions-1,new Operator(Right),true);
  
  if(this->N_fermions>2)
    mpo.setOp(1,new Operator(Middle),true);
  for(int i=2; i<this->N_fermions-1; i++)
    mpo.setOp(i,&mpo.getOp(1),false);
}

void SU2HamiltonianNoLinks::getChargeSquareMPO(MPO &mpo) const
{
  mpo.initLength(this->N_fermions);
  
  mwArray QsquareOp(Indices(this->d_fermi,this->d_fermi));
  QsquareOp.fillWithZero();
  QsquareOp.setElement(0.25*ONE_c,Indices(1,1));
  QsquareOp.setElement(0.25*ONE_c,Indices(2,2));
  
  int Dbond=2;
  mwArray Operators,First,Last,Middle;
  Operators = mwArray(Indices(2,this->d_fermi,this->d_fermi));	Operators.fillWithZero();
  First = mwArray(Indices(1,Dbond,2));				First.fillWithZero();
  Middle = mwArray(Indices(Dbond,Dbond,2));			Middle.fillWithZero();
  Last = mwArray(Indices(Dbond,1,2));				Last.fillWithZero();
  
  
  //Operators are anyway diagonal
  for(int i=0; i<this->d_fermi; i++)
  {
    Operators.setElement(this->id_fermi.getElement(Indices(i,i)),Indices(0,i,i));
    Operators.setElement(QsquareOp.getElement(Indices(i,i)),Indices(1,i,i));
  }
  Operators.reshape(Indices(2,this->d_fermi*this->d_fermi));  
  
  First.setElement(ONE_c,Indices(0,0,0));
  First.setElement(ONE_c,Indices(0,1,1));
  
  Middle.setElement(ONE_c,Indices(0,0,0));
  Middle.setElement(ONE_c,Indices(1,1,0));
  Middle.setElement(ONE_c,Indices(0,1,1));
  
  Last.setElement(ONE_c,Indices(1,0,0));
  Last.setElement(ONE_c,Indices(0,0,1));
  
  First.reshape(Indices(1*Dbond,2));		First.multiplyRight(Operators);		First.reshape(Indices(1,Dbond,this->d_fermi,this->d_fermi));		First.permute(Indices(3,1,4,2));
  Middle.reshape(Indices(Dbond*Dbond,2));	Middle.multiplyRight(Operators);	Middle.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Middle.permute(Indices(3,1,4,2));
  Last.reshape(Indices(Dbond*1,2));		Last.multiplyRight(Operators);		Last.reshape(Indices(Dbond,1,this->d_fermi,this->d_fermi));		Last.permute(Indices(3,1,4,2));
  
  int pos_saved;
  bool saved=false;
  
  mpo.setOp(0,new Operator(First),true);
  mpo.setOp(this->N_fermions-1,new Operator(Last),true);
  
  for(int i=1; i<(this->N_fermions-1); i++)
  {
    if(saved)
      mpo.setOp(i,&mpo.getOp(pos_saved),false);
    else
    {
      mpo.setOp(i,new Operator(Middle),true);
      pos_saved=i;
      saved=true;
    }
  }  
}

void SU2HamiltonianNoLinks::getPenaltyMPO(MPO &mpo, double strength) const
{
  if(strength==0.0)
    strength=this->lambda;
  
  mpo.initLength(this->N_fermions);
  
  mwArray tensor;
  int Dbond=this->d_link+1;
  
  //Frist site
  tensor = mwArray(Indices(this->d_fermi,1,this->d_fermi,Dbond));
  tensor.fillWithZero();
  for(int d=0; d<this->d_fermi; d++)
  {
      if(d==0 || d==(this->d_fermi-1))
	tensor.setElement(ONE_c,Indices(d,0,d,0));
      else if(d==1)
	tensor.setElement(strength*ONE_c,Indices(d,0,d,Dbond-1));
      else if(d==2)
	tensor.setElement(ONE_c,Indices(d,0,d,1));
  }
  mpo.setOp(0,new Operator(tensor),true);
  
  //Sites in between, count until a reach a violating site
  for(int i=1; i<(this->N_fermions-1); i++)
  {
    tensor = mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));
    tensor.fillWithZero();
    for(int d=0; d<this->d_fermi; d++)
    {
      for(int ind_left=0; ind_left<(Dbond-1); ind_left++)
      {
	if(d==0 || d==(this->d_fermi-1))
	  tensor.setElement(ONE_c,Indices(d,ind_left,d,ind_left));
	else if(d==1 && ind_left-1>=0)
	  tensor.setElement(ONE_c,Indices(d,ind_left,d,ind_left-1));
	else if(d==1)
	  tensor.setElement(strength*ONE_c,Indices(d,ind_left,d,Dbond-1));
	else if(d==2 && ind_left+1<this->d_link)
	  tensor.setElement(ONE_c,Indices(d,ind_left,d,ind_left+1));
	else if(d==2)
	  tensor.setElement(strength*ONE_c,Indices(d,ind_left,d,Dbond-1));
      }
      //For each d, pass on in case I the signal of an over-/underflow
      tensor.setElement(ONE_c,Indices(d,Dbond-1,d,Dbond-1));
    }
    mpo.setOp(i,new Operator(tensor),true);
  }
  
  tensor = mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,1));
  tensor.fillWithZero();
  for(int d=0; d<this->d_fermi; d++)
  {
    for(int ind_left=0; ind_left<(Dbond-1); ind_left++)
    {
      if(d==0 || d==(this->d_fermi-1))
	tensor.setElement(ZERO_c,Indices(d,ind_left,d,0));
      else if(d==1 && ind_left-1>=0)
	tensor.setElement(ZERO_c,Indices(d,ind_left,d,0));
      else if(d==1)
	tensor.setElement(strength*ONE_c,Indices(d,ind_left,d,0));
      else if(d==2 && ind_left+1<this->d_link)
	tensor.setElement(ZERO_c,Indices(d,ind_left,d,0));
      else if(d==2)
	tensor.setElement(strength*ONE_c,Indices(d,ind_left,d,0));
    }
    //Case I already get an underflow signaled
    tensor.setElement(ONE_c,Indices(d,Dbond-1,d,0));
  }
  mpo.setOp(this->N_fermions-1,new Operator(tensor),true); 
}


void SU2HamiltonianNoLinks::getPenaltyMPO2(MPO &mpo, bool check_ln,  double strength)
const
{
  if(strength==0.0)
    strength=this->lambda;
  
  mpo.initLength(this->N_fermions);
  
  mwArray tensor;  
  int Dbond=this->d_link+1;
  
  //Essentially I simply count until I find a violating site, if I have found one, I simply finish; so I have a delta in the physical indices and encode in the bond dimension what is going on 
  
  //First site, l_0 is zero
  tensor=mwArray(Indices(this->d_fermi,1,this->d_fermi,Dbond));	tensor.fillWithZero();
  //i=0
  tensor.setElement(ONE_c,Indices(0,0,0,0));
  //i=1 (flux has to go down and would yield unphysical state
  tensor.setElement(strength*ONE_c,Indices(1,0,1,Dbond-1));
  //i=2 (flux has to go up and as dlink>=2 this is always possible)
  tensor.setElement(ONE_c,Indices(2,0,2,1));
  //i=3
  tensor.setElement(ONE_c,Indices(3,0,3,0));
  
  mpo.setOp(0,new Operator(tensor),true);
  
  //Sites in between
  tensor=mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));	tensor.fillWithZero();
  for(int i=0; i<this->d_link; i++)
  {
    //i=0, 2 (just pass on)
    tensor.setElement(ONE_c,Indices(0,i,0,i));
    tensor.setElement(ONE_c,Indices(3,i,3,i));
    //i=1, check if the flux can still be lowered
    if((i-1)>=0)
      tensor.setElement(ONE_c,Indices(1,i,1,i-1));
    else
      tensor.setElement(ONE_c,Indices(1,i,1,Dbond-1));
    //i=2  check if the flux can still be increased
    if((i+1)<this->d_link)
      tensor.setElement(ONE_c,Indices(2,i,2,i+1));
    else
      tensor.setElement(ONE_c,Indices(2,i,2,Dbond-1));    
  }
  //For all values of the physical index pass on, if there is already a violation
  tensor.setElement(ONE_c,Indices(0,Dbond-1,0,Dbond-1));
  tensor.setElement(ONE_c,Indices(1,Dbond-1,1,Dbond-1));
  tensor.setElement(ONE_c,Indices(2,Dbond-1,2,Dbond-1));
  tensor.setElement(ONE_c,Indices(3,Dbond-1,3,Dbond-1));
  
  if(this->N_fermions>2)
    mpo.setOp(1,new Operator(tensor),true);
  for(int i=2; i<this->N_fermions-1; i++)
    mpo.setOp(i,&mpo.getOp(1),false);
  
  //Last site
  tensor = mwArray(Indices(this->d_fermi,Dbond,this->d_fermi,1));	tensor.fillWithZero();
  for(int i=0; i<this->d_link; i++)
  {
    //I only care about the fact that the given state falls in my physical Hilbert space
    if(!check_ln)
    {
      //d=1
      if((i-1)<0)
	tensor.setElement(strength*ONE_c,Indices(1,i,1,0));
      //d=2
      if((i+1)>=this->d_link)
	tensor.setElement(strength*ONE_c,Indices(2,i,2,0));
    }
    else
    {
      //Check if the configuration results in l_n=0, otherwise give a signal
      if(i!=0)
      {
	tensor.setElement(strength*ONE_c,Indices(0,i,0,0));
	tensor.setElement(strength*ONE_c,Indices(3,i,3,0));
      }
      if(i!=1)
	tensor.setElement(strength*ONE_c,Indices(1,i,1,0));
      //No matter which value of l_n-1 there is, l_n will never be zero if i=2, therefore penalize all    
      tensor.setElement(strength*ONE_c,Indices(2,i,2,0));
    }     
  }
  //Now take into account if I receive a signal already from somewhere before the last site
  tensor.setElement(ONE_c,Indices(0,Dbond-1,0,0));
  tensor.setElement(ONE_c,Indices(1,Dbond-1,1,0));
  tensor.setElement(ONE_c,Indices(2,Dbond-1,2,0));
  tensor.setElement(ONE_c,Indices(3,Dbond-1,3,0));
  
  mpo.setOp(this->N_fermions-1,new Operator(tensor),true);
#ifdef MATOUT
  mpo.exportForMatlab("Penalty.m");
#endif
}


void SU2HamiltonianNoLinks::getParticleBalanceMPO(MPO &mpo) const
{
  mpo.initLength(this->N_fermions);
  
  int Dbond=2;
  mwArray First,Last,Odd,Even,Ops;
  First=mwArray(Indices(1,Dbond,3));	First.fillWithZero();
  Odd=mwArray(Indices(Dbond,Dbond,3));	Odd.fillWithZero();
  Even=mwArray(Indices(Dbond,Dbond,3));	Even.fillWithZero();
  Last=mwArray(Indices(Dbond,1,3));	Last.fillWithZero();
  Ops=mwArray(Indices(3,this->d_fermi,this->d_fermi));	Ops.fillWithZero();
  
  //Get the operator N_{f,n} - \left(1-(-1)^n\right) for odd and even n
  for(int i=0; i<this->d_fermi; i++)
    Ops.setElement(this->id_fermi.getElement(Indices(i,i)),Indices(0,i,i));
  
  Ops.setElement(this->Noc.getElement(Indices(1,1)),Indices(1,1,1));
  Ops.setElement(this->Noc.getElement(Indices(2,2)),Indices(1,2,2));
  Ops.setElement(this->Noc.getElement(Indices(3,3)),Indices(1,3,3));
  
  Ops.setElement(-1.0*this->Noc.getElement(Indices(3,3)),Indices(2,0,0));
  Ops.setElement(-1.0*this->Noc.getElement(Indices(1,1)),Indices(2,1,1));
  Ops.setElement(-1.0*this->Noc.getElement(Indices(2,2)),Indices(2,2,2));

  /*cout << "Op1=" << Ops.subArray(Indices(1,-1,-1)) << endl;
  cout << "Op2=" << Ops.subArray(Indices(2,-1,-1)) << endl;*/
  
  Ops.reshape(Indices(3,this->d_fermi*this->d_fermi));
  
  First.setElement(ONE_c,Indices(0,0,0));
  First.setElement(ONE_c,Indices(0,1,2));
  
  Odd.setElement(ONE_c,Indices(0,0,0));
  Odd.setElement(ONE_c,Indices(1,1,0));
  Odd.setElement(ONE_c,Indices(0,1,2));
  
  Even.setElement(ONE_c,Indices(0,0,0));
  Even.setElement(ONE_c,Indices(1,1,0));
  Even.setElement(ONE_c,Indices(0,1,1));
  
  Last.setElement(ONE_c,Indices(1,0,0));
  if((this->N_fermions%2)==0)
    Last.setElement(ONE_c,Indices(0,0,1));
  else
    Last.setElement(ONE_c,Indices(0,0,2));
  
  First.reshape(Indices(1*Dbond,3));	First.multiplyRight(Ops);	First.reshape(Indices(1,Dbond,this->d_fermi,this->d_fermi));	First.permute(Indices(3,1,4,2));
  Odd.reshape(Indices(Dbond*Dbond,3));	Odd.multiplyRight(Ops);		Odd.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Odd.permute(Indices(3,1,4,2));
  Even.reshape(Indices(Dbond*Dbond,3));	Even.multiplyRight(Ops);	Even.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Even.permute(Indices(3,1,4,2));
  Last.reshape(Indices(Dbond*1,3));	Last.multiplyRight(Ops);	Last.reshape(Indices(Dbond,1,this->d_fermi,this->d_fermi));	Last.permute(Indices(3,1,4,2));
  
  mpo.setOp(0,new Operator(First),true);
  mpo.setOp(this->N_fermions-1,new Operator(Last),true);
  
  bool odd_saved=false,even_saved=false;
  int pos_odd_saved,pos_even_saved;
  for(int i=1; i<(this->N_fermions-1); i++)
  {
    if((i+1)%2==0 && even_saved)
      mpo.setOp(i,&mpo.getOp(pos_even_saved),false);
    else if((i+1)%2==0)
    {
      mpo.setOp(i,new Operator(Even),true);
      pos_even_saved=i;
      even_saved=true;
    }
    else if((i+1)%2!=0 && odd_saved)
      mpo.setOp(i,&mpo.getOp(pos_odd_saved),false);
    else if((i+1)%2!=0)
    {
      mpo.setOp(i,new Operator(Odd),true);
      pos_odd_saved=i;
      odd_saved=true;
    }
    else
      cout << "Error, unexpected case" << endl;      
  }
  
#ifdef MATOUT
  mpo.exportForMatlab("particle_diff.m");
#endif
}

void SU2HamiltonianNoLinks::getParticleBalancePenaltyMPO(MPO &mpo) const
{
  mpo.initLength(this->N_fermions);
  
  int Dbond=3;  
  mwArray Operators, First, Last, Odd, Even;
  Operators=mwArray(Indices(5,this->d_fermi,this->d_fermi));
  First=mwArray(Indices(1,Dbond,5));	First.fillWithZero();
  Odd=mwArray(Indices(Dbond,Dbond,5));	Odd.fillWithZero();
  Even=mwArray(Indices(Dbond,Dbond,5));	Even.fillWithZero();
  Last=mwArray(Indices(Dbond,1,5));	Last.fillWithZero();
  
  mwArray OpOdd,OpEven,OpOddsquare, OpEvensquare;
  OpOdd = this->Noc-2.0*identityMatrix(this->d_fermi);
  OpEven = this->Noc;
  
  OpOddsquare = OpOdd*OpOdd;
  OpEvensquare = OpEven*OpEven;
  
  /*cout << "OpOdd=" << OpOdd << endl;
  cout << "OpOddsquare=" << OpOddsquare << endl;
  cout << "OpEven=" << OpEven << endl;
  cout << "OpEvensquare=" << OpEvensquare << endl;  */
  
  for (int i=0; i<this->d_fermi; i++)
  {
    Operators.setElement(this->id_fermi.getElement(Indices(i,i)),Indices(0,i,i));
    Operators.setElement(OpOdd.getElement(Indices(i,i)),Indices(1,i,i));
    Operators.setElement(OpEven.getElement(Indices(i,i)),Indices(2,i,i));
    Operators.setElement(OpOddsquare.getElement(Indices(i,i)),Indices(3,i,i));
    Operators.setElement(OpEvensquare.getElement(Indices(i,i)),Indices(4,i,i));
  }
  
  First.setElement(ONE_c,Indices(0,0,0));
  First.setElement(2.0*ONE_c,Indices(0,1,1));
  First.setElement(ONE_c,Indices(0,Dbond-1,3));
  
  Odd.setElement(ONE_c,Indices(0,0,0));
  Odd.setElement(ONE_c,Indices(1,1,0));
  Odd.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  Odd.setElement(2.0*ONE_c,Indices(0,1,1));
  Odd.setElement(ONE_c,Indices(1,Dbond-1,1));
  Odd.setElement(ONE_c,Indices(0,Dbond-1,3));
  
  Even.setElement(ONE_c,Indices(0,0,0));
  Even.setElement(ONE_c,Indices(1,1,0));
  Even.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  Even.setElement(2.0*ONE_c,Indices(0,1,2));
  Even.setElement(ONE_c,Indices(1,Dbond-1,2));
  Even.setElement(ONE_c,Indices(0,Dbond-1,4));
  
  Last.setElement(ONE_c,Indices(Dbond-1,0,0));
  if((this->N_fermions%2)==0)
  {
    Last.setElement(ONE_c,Indices(1,0,2));
    Last.setElement(ONE_c,Indices(0,0,4));
  }
  else
  {
    Last.setElement(ONE_c,Indices(1,0,1));
    Last.setElement(ONE_c,Indices(0,0,3));
  }
  
  Operators.reshape(Indices(5,this->d_fermi*this->d_fermi));
  
  First.reshape(Indices(1*Dbond,5));	First.multiplyRight(Operators);	First.reshape(Indices(1,Dbond,this->d_fermi,this->d_fermi));	First.permute(Indices(3,1,4,2));
  Odd.reshape(Indices(Dbond*Dbond,5));	Odd.multiplyRight(Operators);	Odd.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Odd.permute(Indices(3,1,4,2));
  Even.reshape(Indices(Dbond*Dbond,5));	Even.multiplyRight(Operators);	Even.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Even.permute(Indices(3,1,4,2));
  Last.reshape(Indices(Dbond*1,5));	Last.multiplyRight(Operators);	Last.reshape(Indices(Dbond,1,this->d_fermi,this->d_fermi));	Last.permute(Indices(3,1,4,2));
  
  /*cout << "First=" << endl;
  print_tensor(First);
  cout << "Odd=" << endl;
  print_tensor(Odd);
  cout << "Even" << endl;
  print_tensor(Even);
  cout << "Last=" << endl;
  print_tensor(Last);*/
  
  mpo.setOp(0,new Operator(First),true);
  mpo.setOp(this->N_fermions-1,new Operator(Last),true);
  
  for(int i=1; i<(this->N_fermions-1); i++)
  {
    if((i%2)==0)
    {
      //Odd site
      mpo.setOp(i,new Operator(Odd),true);
    }
    else
      mpo.setOp(i,new Operator(Even),true);
  }
#ifdef MATOUT
  mpo.exportForMatlab("particle_balance_penalty.m");
#endif
}



void SU2HamiltonianNoLinks::getMomentumMPO(MPO &mpo) const
{
  /***********************************************************
   * Step 1: construct the momentum MPO in a different basis
   * *********************************************************/
  int Dbond=6;
  mwArray Operators;  
  mwArray Tleft,Tright,Tmiddle;  
  Tleft=mwArray(Indices(1,Dbond,3));		Tleft.fillWithZero();
  Tmiddle=mwArray(Indices(Dbond,Dbond,3));	Tmiddle.fillWithZero();
  Tright=mwArray(Indices(Dbond,1,3));		Tright.fillWithZero();
  Operators=mwArray(Indices(3,3,3));		Operators.fillWithZero();
  
  Operators.setElement(ONE_c,Indices(0,0,0));
  Operators.setElement(ONE_c,Indices(0,1,1));
  Operators.setElement(ONE_c,Indices(0,2,2));
  
  Operators.setElement(ONE_c,Indices(1,1,0));
  Operators.setElement(ONE_c,Indices(1,2,1));
  
  Operators.setElement(ONE_c,Indices(2,0,1));
  Operators.setElement(ONE_c,Indices(2,1,2));
  
  cout << "Id=" << Operators.subArray(Indices(0,-1,-1)) << endl;
  cout << "Phi_plus=" << Operators.subArray(Indices(1,-1,-1)) << endl;
  cout << "Phi_minus=" << Operators.subArray(Indices(2,-1,-1)) << endl;
  
  Operators.reshape(Indices(3,3*3));
  
  Tleft.setElement(ONE_c,Indices(0,0,0));
  Tleft.setElement(ONE_c,Indices(0,1,1));
  Tleft.setElement(ONE_c,Indices(0,2,2));
  
  Tmiddle.setElement(ONE_c,Indices(0,0,0));
  Tmiddle.setElement(ONE_c,Indices(1,3,0));
  Tmiddle.setElement(ONE_c,Indices(2,4,0));
  Tmiddle.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  
  Tmiddle.setElement(ONE_c,Indices(0,1,1));
  Tmiddle.setElement(ONE_c,Indices(4,Dbond-1,1));
  
  Tmiddle.setElement(ONE_c,Indices(0,2,2));
  Tmiddle.setElement(ONE_c,Indices(3,Dbond-1,2));
  
  Tright.setElement(ONE_c,Indices(Dbond-1,0,0));
  Tright.setElement(ONE_c,Indices(4,0,1));
  Tright.setElement(ONE_c,Indices(3,0,2));
  
  //Reshape and contract to obtain the three different tensors
  Tleft.reshape(Indices(1*Dbond,3));		Tleft.multiplyRight(Operators);		Tleft.reshape(Indices(1,Dbond,3,3));		Tleft.permute(Indices(3,1,4,2));
  Tmiddle.reshape(Indices(Dbond*Dbond,3));	Tmiddle.multiplyRight(Operators);	Tmiddle.reshape(Indices(Dbond,Dbond,3,3));	Tmiddle.permute(Indices(3,1,4,2));
  Tright.reshape(Indices(Dbond*1,3));		Tright.multiplyRight(Operators);	Tright.reshape(Indices(Dbond,1,3,3));		Tright.permute(Indices(3,1,4,2));
  
  //For debugging
  cout << "Warning replacing momentum MPO with identity for debugging" << endl;
  Dbond=1;
  Tleft=identityMatrix(3);	Tleft.reshape(Indices(3,1,3,1));
  Tmiddle=identityMatrix(3);	Tmiddle.reshape(Indices(3,1,3,1));
  Tright=identityMatrix(3);	Tright.reshape(Indices(3,1,3,1));  

  
#ifdef MATOUT
  MPO tmp_mpo(this->N_fermions);
  tmp_mpo.setOp(0,new Operator(Tleft),true);
  tmp_mpo.setOp(this->N_fermions-1,new Operator(Tright),true);
  if(this->N_fermions>2)
    tmp_mpo.setOp(1,new Operator(Tmiddle),true);
  for(int i=2; i<this->N_fermions-1; i++)
    tmp_mpo.setOp(i,&tmp_mpo.getOp(1),false);
  tmp_mpo.exportForMatlab("Pmpo_nocbasis.m");
#endif
  
  /***********************************************************
   * Step 2: Get the A and B tensors
   * *********************************************************/
  cout << "Step 2: Get the A and B tensors" << endl;
  
  mwArray Aleft,Amiddle,Aright,Bmiddle;
  
  Aleft=mwArray(Indices(3,1,this->d_fermi,this->d_link));			Aleft.fillWithZero();
  Amiddle=mwArray(Indices(3,this->d_link,this->d_fermi,this->d_link));	Amiddle.fillWithZero();
  Aright=mwArray(Indices(3,this->d_link,this->d_fermi,1));			Aright.fillWithZero();
  Bmiddle=mwArray(Indices(this->d_link,this->d_link,this->d_link));		Bmiddle.fillWithZero();
  
  for(int i=0; i<this->d_link; i++)
  {
    Amiddle.setElement(ONE_c,Indices(0,i,0,i));
    Amiddle.setElement(ONE_c,Indices(2,i,3,i));
    
    Aright.setElement(ONE_c,Indices(0,i,0,0));
    Aright.setElement(ONE_c,Indices(2,i,3,0));
    Aright.setElement(ONE_c,Indices(1,i,1,0));
    Aright.setElement(ONE_c,Indices(1,i,2,0));
    
    if((i-1)>=0)
      Amiddle.setElement(ONE_c,Indices(1,i,1,i-1));
    if((i+1)<this->d_link)
      Amiddle.setElement(ONE_c,Indices(1,i,2,i+1));
    
   Bmiddle.setElement(ONE_c,Indices(i,i,i));
  }  
  Aleft.setElement(ONE_c,Indices(0,0,0,0));
  Aleft.setElement(ONE_c,Indices(1,0,2,1));
  Aleft.setElement(ONE_c,Indices(2,0,3,0));  
  
#ifdef CHECKDIMS
  cout << "Warning, giving only the conversion MPO!" << endl;
  mpo.initLength(this->N_fermions);
  
  cout << "Aleft: " << Aleft.getDimensions() << endl;
  cout << "Amiddle: " << Amiddle.getDimensions() << endl;
  cout << "Aright: " << Aright.getDimensions() << endl; 
  
  /*Aleft.permute(Indices(3,2,1,4));
  Aright.permute(Indices(3,2,1,4));
  Amiddle.permute(Indices(3,2,1,4));
  
  cout << "Aleft nach permute: " << Aleft.getDimensions() << endl;
  cout << "Amiddle nach permute: " << Amiddle.getDimensions() << endl;
  cout << "Aright nach permute: " << Aright.getDimensions() << endl;*/
  
  mpo.setOp(0,new Operator(Aleft),true);
  mpo.setOp(this->N_fermions-1,new Operator(Aright),true);
  if(this->N_fermions>2)
    mpo.setOp(1,new Operator(Amiddle),true);
  for(int i=2;i<this->N_fermions-1; i++)
    mpo.setOp(i,&mpo.getOp(1),false);
  
#else
  /***********************************************************
   * Step 3: Contract everything together and build the MPO
   * *********************************************************/
  cout << "Step 3: Contract everything together and build the MPO" << endl;
  mpo.initLength(this->N_fermions);
  
  mwArray A,B,S,T;  
  /////////////
  //First site
  /////////////
  //Get copy of A, B and T as contract leg destroys shape
  A=Aleft;	B=Bmiddle;	T=Tleft;
  //Now contract A with B and get a copy of the object
  A.contractLeg(4,B,2);
  S=A;
  //Contract the result with the MPO tensor of the momentum operator
  A.contractLeg(1,T,3);
  //Now contract the rest of the legs
  A.contractLegs(Indices(3,5),S,Indices(4,1));
  //Reshape and permute to have the final MPO tensor
  cout << "A=" << A.getDimensions() << endl;
  A.permute(Indices(7,6,4,1,2,8,5,3));
  A.reshape(Indices(this->d_fermi,1,this->d_fermi,Dbond*this->d_link*this->d_link));
  
  mpo.setOp(0,new Operator(A),true);
  
  /////////////
  //Last site
  /////////////  
  //Get copy of A, T as contract leg destroys shape
  A=Aright;	T=Tright;
  //Contract A with T
  A.contractLeg(1,T,3);
  //Contract with the hermitian conjugate (take Aright directly, as it is not needed anymore afterwards)
  A.contractLeg(4,Aright,1);
  //Reshape and permute to have the final MPO tensor
  cout << "A=" << A.getDimensions() << endl;
  A.permute(Indices(7,6,4,1,2,8,5,3));
  A.reshape(Indices(this->d_fermi,Dbond*this->d_link*this->d_link,this->d_fermi,1));
  
  mpo.setOp(this->N_fermions-1,new Operator(A),true);  
  
  /////////////////////
  //Site in the middle
  /////////////////////
  //Get copy of A, B and T as contract leg destroys shape
  A=Amiddle;	B=Bmiddle;	T=Tmiddle;
  //Now contract A with B and get a copy of the object
  A.contractLeg(4,B,2);
  S=A;
  //Contract the result with the MPO tensor of the momentum operator
  A.contractLeg(1,T,3);
  //Now contract the rest of the legs
  A.contractLegs(Indices(3,5),S,Indices(4,1));
  //Reshape and permute to have the final MPO tensor
  cout << "A=" << A.getDimensions() << endl;
  A.permute(Indices(7,6,4,1,2,8,5,3));
  A.reshape(Indices(this->d_fermi,Dbond*this->d_link*this->d_link,this->d_fermi,Dbond*this->d_link*this->d_link));
  
  if(this->N_fermions>2)
    mpo.setOp(1,new Operator(A),true);
  for (int i=2; i<(this->N_fermions-1); i++)
    mpo.setOp(i,&mpo.getOp(1), false);
  
  
  
  cout << "Mega MPO: " << mpo << endl;
#ifdef MATOUT
  mpo.exportForMatlab("Pmpo.m");
#endif
#endif
}

void SU2HamiltonianNoLinks::getNocMPO_other_basis(MPO &mpo, int n) const
{
  mwArray NocMPO=mwArray(Indices(3,1,3,1));
  NocMPO.fillWithZero();
  NocMPO.setElement(ONE_c,Indices(1,0,1,0));
  NocMPO.setElement(2.0*ONE_c,Indices(2,0,2,0));
  
  //print_tensor(NocMPO);
  
  mwArray idmatrix = identityMatrix(3);
  idmatrix.reshape(Indices(3,1,3,1));
  
  mpo.initLength(this->N_fermions);
  mpo.setOp(n-1, new Operator(NocMPO), true);
  for(int i=0; i<this->N_fermions; i++)
    if(i!=(n-1))
      mpo.setOp(i,new Operator(idmatrix),true);
    
}

void SU2HamiltonianNoLinks::getMomentumMPO2(MPO &mpo) const
{
  /***********************************************************
   * Step 1: construct the momentum MPO in a different basis
   * *********************************************************/
  int Dbond=6;
  mwArray Operators;  
  mwArray Tleft,Tright,Tmiddle;  
  Tleft=mwArray(Indices(1,Dbond,3));		Tleft.fillWithZero();
  Tmiddle=mwArray(Indices(Dbond,Dbond,3));	Tmiddle.fillWithZero();
  Tright=mwArray(Indices(Dbond,1,3));		Tright.fillWithZero();
  Operators=mwArray(Indices(3,3,3));		Operators.fillWithZero();
  
  Operators.setElement(ONE_c,Indices(0,0,0));
  Operators.setElement(ONE_c,Indices(0,1,1));
  Operators.setElement(ONE_c,Indices(0,2,2));
  
  Operators.setElement(ONE_c,Indices(1,1,0));
  Operators.setElement(ONE_c,Indices(1,2,1));
  
  Operators.setElement(ONE_c,Indices(2,0,1));
  Operators.setElement(ONE_c,Indices(2,1,2));
  
  Operators.reshape(Indices(3,3*3));
  
  complex_t c1=ONE_c;
  complex_t c2=-ONE_c;  
  //Uncomment for making something -i\Phi^\dagger(n)\Phi(n+2)+i\Phi(n)\Phi(n+2)^\dagger
  c1 = I_c;
  c2 = -I_c;
  
  Tleft.setElement(ONE_c,Indices(0,0,0));
  Tleft.setElement(c1*ONE_c,Indices(0,1,1));
  Tleft.setElement(c2*ONE_c,Indices(0,2,2));
  
  Tmiddle.setElement(ONE_c,Indices(0,0,0));
  Tmiddle.setElement(ONE_c,Indices(1,3,0));
  Tmiddle.setElement(ONE_c,Indices(2,4,0));
  Tmiddle.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  
  Tmiddle.setElement(c1*ONE_c,Indices(0,1,1));
  Tmiddle.setElement(ONE_c,Indices(4,Dbond-1,1));
  
  Tmiddle.setElement(c2*ONE_c,Indices(0,2,2));
  Tmiddle.setElement(ONE_c,Indices(3,Dbond-1,2));
  
  Tright.setElement(ONE_c,Indices(Dbond-1,0,0));
  Tright.setElement(ONE_c,Indices(4,0,1));
  Tright.setElement(ONE_c,Indices(3,0,2));
  
  //Reshape and contract to obtain the three different tensors
  Tleft.reshape(Indices(1*Dbond,3));		Tleft.multiplyRight(Operators);		Tleft.reshape(Indices(1,Dbond,3,3));		Tleft.permute(Indices(3,1,4,2));
  Tmiddle.reshape(Indices(Dbond*Dbond,3));	Tmiddle.multiplyRight(Operators);	Tmiddle.reshape(Indices(Dbond,Dbond,3,3));	Tmiddle.permute(Indices(3,1,4,2));
  Tright.reshape(Indices(Dbond*1,3));		Tright.multiplyRight(Operators);	Tright.reshape(Indices(Dbond,1,3,3));		Tright.permute(Indices(3,1,4,2));
  
  //Construct the square of the MPO
  /*mwArray dummyT=Tleft;
  Tleft.contractLeg(3,dummyT,1);
  Tleft.permute(Indices(1,2,4,5,3,6));
  Tleft.reshape(Indices(3,1,3,Dbond*Dbond));
  
  dummyT=Tmiddle;
  Tmiddle.contractLeg(3,dummyT,1);
  Tmiddle.permute(Indices(1,2,4,5,3,6));
  Tmiddle.reshape(Indices(3,Dbond*Dbond,3,Dbond*Dbond));
  
  dummyT=Tright;
  Tright.contractLeg(3,dummyT,1);
  Tright.permute(Indices(1,2,4,5,3,6));
  Tright.reshape(Indices(3,Dbond*Dbond,3,1));
  
  Dbond=Dbond*Dbond;*/
  
  //For debugging, put only identities
  /*Tleft=identityMatrix(3);	Tleft.reshape(Indices(3,1,3,1));
  Tmiddle=Tleft;
  Tright=Tleft;
  Dbond=1;*/

  
#ifdef MATOUT
  MPO tmp_mpo(this->N_fermions);
  tmp_mpo.setOp(0,new Operator(Tleft),true);
  tmp_mpo.setOp(this->N_fermions-1,new Operator(Tright),true);
  if(this->N_fermions>2)
    tmp_mpo.setOp(1,new Operator(Tmiddle),true);
  for(int i=2; i<this->N_fermions-1; i++)
    tmp_mpo.setOp(i,&tmp_mpo.getOp(1),false);
  tmp_mpo.exportForMatlab("Pmpo2_nocbasis.m");
#endif
  
  /***********************************************************
   * Step 2: Get the A and B tensors
   * *********************************************************/
  //Here I don't care about the gauge field, I simply map the fermionic content
  cout << "Step 2: Get the A " << endl;
  
  mwArray Aleft,Amiddle,Aright,Bmiddle;
  
  Aleft=mwArray(Indices(3,1,this->d_fermi,1));		Aleft.fillWithZero();
  Amiddle=mwArray(Indices(3,1,this->d_fermi,1));	Amiddle.fillWithZero();
  Aright=mwArray(Indices(3,1,this->d_fermi,1));		Aright.fillWithZero();
  
  Aleft.setElement(ONE_c,Indices(0,0,0,0));
  Aleft.setElement(ONE_c,Indices(1,0,1,0));
  Aleft.setElement(ONE_c,Indices(1,0,2,0));
  Aleft.setElement(ONE_c,Indices(2,0,3,0));
  
  Amiddle=Aleft;
  Aright=Aleft;
#ifdef CHECKDIMS
  cout << "Warning, giving only the conversion MPO!" << endl;
  mpo.initLength(this->N_fermions);
  
  cout << "Aleft: " << Aleft.getDimensions() << endl;
  cout << "Amiddle: " << Amiddle.getDimensions() << endl;
  cout << "Aright: " << Aright.getDimensions() << endl; 
  
  /*Aleft.permute(Indices(3,2,1,4));
  Aright.permute(Indices(3,2,1,4));
  Amiddle.permute(Indices(3,2,1,4));
  
  cout << "Aleft nach permute: " << Aleft.getDimensions() << endl;
  cout << "Amiddle nach permute: " << Amiddle.getDimensions() << endl;
  cout << "Aright nach permute: " << Aright.getDimensions() << endl;*/
  
  mpo.setOp(0,new Operator(Aleft),true);
  mpo.setOp(this->N_fermions-1,new Operator(Aright),true);
  if(this->N_fermions>2)
    mpo.setOp(1,new Operator(Amiddle),true);
  for(int i=2;i<this->N_fermions-1; i++)
    mpo.setOp(i,&mpo.getOp(1),false);
  
#else
  /***********************************************************
   * Step 3: Contract everything together and build the MPO
   * *********************************************************/
  cout << "Step 3: Contract everything together and build the MPO" << endl;
  mpo.initLength(this->N_fermions);
  
  mwArray A,B,S,T;  
  /////////////
  //First site
  /////////////
  //Get copy of A, B and T as contract leg destroys shape
  A=Aleft;	T=Tleft;
  //Contract A with T
  A.contractLeg(1,T,3);
  //Contract with the hermitian conjugate (take Aleft directly, as it is not needed anymore afterwards)
  A.contractLeg(4,Aleft,1);
  //Reshape and permute to have the final MPO tensor
  cout << "A=" << A.getDimensions() << endl;
  A.permute(Indices(7,6,4,1,2,8,5,3));
  A.reshape(Indices(this->d_fermi,1,this->d_fermi,Dbond));  

  mpo.setOp(0,new Operator(A),true);
  
  /////////////
  //Last site
  /////////////  
  //Get copy of A, T as contract leg destroys shape
  A=Aright;	T=Tright;
  //Contract A with T
  A.contractLeg(1,T,3);
  //Contract with the hermitian conjugate (take Aright directly, as it is not needed anymore afterwards)
  A.contractLeg(4,Aright,1);
  //Reshape and permute to have the final MPO tensor
  cout << "A=" << A.getDimensions() << endl;
  A.permute(Indices(7,6,4,1,2,8,5,3));
  A.reshape(Indices(this->d_fermi,Dbond,this->d_fermi,1));
  
  mpo.setOp(this->N_fermions-1,new Operator(A),true);  
  
  /////////////////////
  //Site in the middle
  /////////////////////
  //Get copy of A, B and T as contract leg destroys shape
  A=Amiddle;	T=Tmiddle;
  //Contract A with T
  A.contractLeg(1,T,3);
  //Contract with the hermitian conjugate (take Amiddle directly, as it is not needed anymore afterwards)
  A.contractLeg(4,Amiddle,1);
  //Reshape and permute to have the final MPO tensor
  cout << "A=" << A.getDimensions() << endl;
  A.permute(Indices(7,6,4,1,2,8,5,3));
  A.reshape(Indices(this->d_fermi,Dbond,this->d_fermi,Dbond));  
  
  if(this->N_fermions>2)
    mpo.setOp(1,new Operator(A),true);
  for (int i=2; i<(this->N_fermions-1); i++)
    mpo.setOp(i,&mpo.getOp(1), false); 
  
  //cout << "Mega MPO2: " << mpo << endl;
#ifdef MATOUT
  mpo.exportForMatlab("Pmpo2.m");
#endif
#endif
}

void SU2HamiltonianNoLinks::getMomentumMPO3(MPO &mpo) const
{
  //Set up operators
  mwArray Op1,Op2,Op3,Op4,Op1dag,Op2dag,Op3dag,Op4dag,D1,D2;
  
  Op1=mwArray(Indices(this->d_fermi,this->d_fermi));	Op1.fillWithZero();
  Op2=mwArray(Indices(this->d_fermi,this->d_fermi));	Op2.fillWithZero();
  Op3=mwArray(Indices(this->d_fermi,this->d_fermi));	Op3.fillWithZero();
  Op4=mwArray(Indices(this->d_fermi,this->d_fermi));	Op4.fillWithZero();
  D1=mwArray(Indices(this->d_fermi,this->d_fermi));	D1.fillWithZero();
  
  Op1.setElement(ONE_c,Indices(0,1));
  Op2.setElement(ONE_c,Indices(0,2));
  Op3.setElement(ONE_c,Indices(1,3));
  Op4.setElement(ONE_c,Indices(2,3));
  D1.setElement(ONE_c,Indices(1,2));
  
  Op1dag = Op1;	Op1dag.Hconjugate();
  Op2dag = Op2;	Op2dag.Hconjugate();
  Op3dag = Op3;	Op3dag.Hconjugate();
  Op4dag = Op4;	Op4dag.Hconjugate();
  D2 = D1;	D2.Hconjugate();
  
  //For debugging 
  /*cout << "Op1=" << Op1 << endl;
  cout << "Op2=" << Op2 << endl;
  cout << "Op3=" << Op3 << endl;
  cout << "Op4=" << Op4 << endl;
  
  cout << "Op1dag=" << Op1dag << endl;
  cout << "Op2dag=" << Op2dag << endl;
  cout << "Op3dag=" << Op3dag << endl;
  cout << "Op4dag=" << Op4dag << endl;
  
  cout << "D1=" << D1 << endl;
  cout << "D2=" << D2 << endl;*/
  
  mwArray Ops(Indices(11,this->d_fermi,this->d_fermi));  
  for(int i=0; i<this->d_fermi;i++)
  {
    for(int j=0; j<this->d_fermi;j++)
    {
      Ops.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      Ops.setElement(Op1.getElement(Indices(i,j)),Indices(1,i,j));
      Ops.setElement(Op2.getElement(Indices(i,j)),Indices(2,i,j));
      Ops.setElement(Op3.getElement(Indices(i,j)),Indices(3,i,j));
      Ops.setElement(Op4.getElement(Indices(i,j)),Indices(4,i,j));
      Ops.setElement(Op1dag.getElement(Indices(i,j)),Indices(5,i,j));
      Ops.setElement(Op2dag.getElement(Indices(i,j)),Indices(6,i,j));
      Ops.setElement(Op3dag.getElement(Indices(i,j)),Indices(7,i,j));
      Ops.setElement(Op4dag.getElement(Indices(i,j)),Indices(8,i,j));
      Ops.setElement(D1.getElement(Indices(i,j)),Indices(9,i,j));
      Ops.setElement(D2.getElement(Indices(i,j)),Indices(10,i,j));
    }
  }
  Ops.reshape(Indices(11,this->d_fermi*this->d_fermi));
  int Dbond=22;
  complex_t c1,c2;
  mwArray Mleft,M,Mright;
  M=mwArray(Indices(Dbond,Dbond,11));	M.fillWithZero();
  
  c1= I_c;
  c2= -I_c;
  
  M.setElement(ONE_c,Indices(0,0,0));
  M.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  M.setElement(ONE_c,Indices(1,9,0));
  M.setElement(ONE_c,Indices(2,10,0));
  M.setElement(ONE_c,Indices(3,11,0));
  M.setElement(ONE_c,Indices(4,12,0));
  M.setElement(ONE_c,Indices(5,13,0));
  M.setElement(ONE_c,Indices(6,14,0));
  M.setElement(ONE_c,Indices(7,15,0));
  M.setElement(ONE_c,Indices(8,16,0));
  
  M.setElement(c2*ONE_c,Indices(0,1,1));
  M.setElement(ONE_c,Indices(13,Dbond-1,1));
  M.setElement(ONE_c,Indices(16,Dbond-1,1));
  M.setElement(ONE_c,Indices(18,Dbond-1,1));
  
  M.setElement(c2*ONE_c,Indices(0,2,2));
  M.setElement(ONE_c,Indices(14,Dbond-1,2));
  M.setElement(ONE_c,Indices(15,Dbond-1,2));
  M.setElement(ONE_c,Indices(20,Dbond-1,2));
  
  M.setElement(c2*ONE_c,Indices(0,3,3));
  M.setElement(ONE_c,Indices(15,Dbond-1,3));
  M.setElement(ONE_c,Indices(14,Dbond-1,3));
  M.setElement(ONE_c,Indices(20,Dbond-1,3));
  
  M.setElement(c2*ONE_c,Indices(0,4,4));
  M.setElement(ONE_c,Indices(16,Dbond-1,4));
  M.setElement(ONE_c,Indices(13,Dbond-1,4));
  M.setElement(ONE_c,Indices(18,Dbond-1,4));
  
  M.setElement(c1*ONE_c,Indices(0,5,5));
  M.setElement(ONE_c,Indices(9,Dbond-1,5));
  M.setElement(ONE_c,Indices(12,Dbond-1,5));
  M.setElement(ONE_c,Indices(19,Dbond-1,5));
  
  M.setElement(c1*ONE_c,Indices(0,6,6));
  M.setElement(ONE_c,Indices(10,Dbond-1,6));
  M.setElement(ONE_c,Indices(11,Dbond-1,6));
  M.setElement(ONE_c,Indices(17,Dbond-1,6));
  
  M.setElement(c1*ONE_c,Indices(0,7,7));
  M.setElement(ONE_c,Indices(11,Dbond-1,7));
  M.setElement(ONE_c,Indices(10,Dbond-1,7));
  M.setElement(ONE_c,Indices(17,Dbond-1,7));
  
  M.setElement(c1*ONE_c,Indices(0,8,8));
  M.setElement(ONE_c,Indices(12,Dbond-1,8));
  M.setElement(ONE_c,Indices(9,Dbond-1,8));
  M.setElement(ONE_c,Indices(19,Dbond-1,8));
  
  M.setElement(ONE_c,Indices(1,17,9));
  M.setElement(ONE_c,Indices(4,17,9));
  M.setElement(ONE_c,Indices(6,18,9));
  M.setElement(ONE_c,Indices(7,18,9));
  
  M.setElement(ONE_c,Indices(2,19,10));
  M.setElement(ONE_c,Indices(3,19,10));
  M.setElement(ONE_c,Indices(5,20,10));
  M.setElement(ONE_c,Indices(8,20,10));
  
  Mleft=M.subArray(Indices(0,-1,-1));
  Mright=M.subArray(Indices(-1,Dbond-1,-1));
  
  
  Mleft.reshape(Indices(1*Dbond,11));	Mleft.multiplyRight(Ops);	Mleft.reshape(Indices(1,Dbond,this->d_fermi,this->d_fermi));	Mleft.permute(Indices(3,1,4,2));
  M.reshape(Indices(Dbond*Dbond,11));	M.multiplyRight(Ops);		M.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	M.permute(Indices(3,1,4,2));
  Mright.reshape(Indices(Dbond*1,11));	Mright.multiplyRight(Ops);	Mright.reshape(Indices(Dbond,1,this->d_fermi,this->d_fermi));	Mright.permute(Indices(3,1,4,2));
  
  mpo.initLength(this->N_fermions);
  mpo.setOp(0,new Operator(Mleft),true);
  mpo.setOp(this->N_fermions-1,new Operator(Mright),true);
  if(this->N_fermions>2)
    mpo.setOp(1,new Operator(M),true);
  for(int i=2; i<(this->N_fermions-1); i++)
    mpo.setOp(i,&mpo.getOp(1),false);
  
#ifdef MATOUT
  mpo.exportForMatlab("Pmpo_new.m");
#endif
  
}


void SU2HamiltonianNoLinks::getCyclicShiftMPO(MPO &mpo) const
{
  mpo.initLength(this->N_fermions);
  int d=this->d_fermi;
  int D=d*d;
  // The basic operator on any middle site
  mwArray basicOp=identityMatrix(d*d*d);
  basicOp.reshape(Indices(d,d,d,d,d,d));
  basicOp.permute(Indices(5,1,2,3,4,6));
  basicOp.reshape(Indices(d,d*d,d,d*d));
  // The operator on the left
  mwArray leftOp=identityMatrix(d*d);
  leftOp.reshape(Indices(d,1,d,d*d));
  mwArray rightOp=identityMatrix(d*d);
  rightOp.reshape(Indices(d,d,d,d));
  rightOp.permute(Indices(4,1,2,3));
  rightOp.reshape(Indices(d,d*d,d,1));

  // Now set the operators in the chain
  mpo.setOp(0,new Operator(leftOp),true);
  mpo.setOp(this->N_fermions-1,new Operator(rightOp),true);
  if(this->N_fermions>2)
    mpo.setOp(1,new Operator(basicOp),true);
  for(int k=2;k<this->N_fermions-1;k++)
    mpo.setOp(k,&mpo.getOp(1),false);
#ifdef MATOUT
  mpo.exportForMatlab("CyclicShiftMPO.m");
#endif
}

void SU2HamiltonianNoLinks::getParticleTrafoMPO(MPO &mpo) const
{
  mpo.initLength(this->N_fermions);
  mwArray Op(Indices(this->d_fermi,1,this->d_fermi,1));
  Op.fillWithZero();
  
  Op.setElement(ONE_c,Indices(0,0,3,0));
  Op.setElement(ONE_c,Indices(3,0,0,0));
  Op.setElement(ONE_c,Indices(1,0,1,0));
  Op.setElement(ONE_c,Indices(2,0,2,0));
  
  mpo.setOp(0,new Operator(Op), true);
  for(int i=1; i<this->N_fermions; i++)
    mpo.setOp(i,&mpo.getOp(0),false);
#ifdef MATOUT
  mpo.exportForMatlab("ParticleTrafo.m");
#endif
}

void SU2HamiltonianNoLinks::getProjectorMPO(MPO &mpo) const
{
  mpo.initLength(this->N_fermions);
  mwArray Op(Indices(this->d_fermi,this->d_fermi));
  Op.fillWithZero();
  Op.setElement(ONE_c,Indices(0,0));
  Op.setElement(ONE_c,Indices(2,2));
  Op.setElement(ONE_c,Indices(3,3));
  Op.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
  
  mpo.setOp(0,new Operator(this->id_fermi_mpo),true);
  for(int i=1; i<this->N_fermions-1; i++)
    mpo.setOp(i,&mpo.getOp(0),false);
  mpo.setOp(this->N_fermions-1,new Operator(Op),true);
}

void SU2HamiltonianNoLinks::getShiftMPO(MPO &mpo, bool flip_particles) const
{
  mpo.initLength(this->N_fermions);
  
  mwArray First,Middle,Last;
  First=mwArray(Indices(this->d_fermi,1,this->d_fermi,this->d_fermi));			First.fillWithZero();
  Middle=mwArray(Indices(this->d_fermi,this->d_fermi,this->d_fermi,this->d_fermi));	Middle.fillWithZero();
  Last=mwArray(Indices(this->d_fermi,this->d_fermi,this->d_fermi,1));			Last.fillWithZero();
  
  for(int i=0; i<this->d_fermi; i++)
    for(int j=0; j<this->d_fermi; j++)
    {
      Middle.setElement(ONE_c,Indices(i,i,j,j));
      //Variant, just set the remaining index to a |0> at the beginning
      //First.setElement(ONE_c,Indices(0,0,j,j));
      //Variant, just copy the physical input of the old state
      First.setElement(ONE_c,Indices(j,0,j,j));
      Last.setElement(ONE_c,Indices(i,i,0,0));
    }
  
  if(flip_particles)
  {
    mwArray Op(Indices(this->d_fermi,1,this->d_fermi,1));	Op.fillWithZero();
    
    Op.setElement(ONE_c,Indices(0,0,3,0));
    Op.setElement(ONE_c,Indices(3,0,0,0));
    Op.setElement(ONE_c,Indices(1,0,1,0));
    Op.setElement(ONE_c,Indices(2,0,2,0));
    
    Op.permute(Indices(1,2,4,3));
    Op.reshape(Indices(this->d_fermi*1*1,this->d_fermi));
    Middle.reshape(Indices(this->d_fermi,this->d_fermi*this->d_fermi*this->d_fermi));
    Last.reshape(Indices(this->d_fermi,this->d_fermi*this->d_fermi*1));
    
    Middle.multiplyLeft(Op);
    Last.multiplyLeft(Op);
    
    Middle.permute(Indices(1,2,4,5,3,6));
    Last.permute(Indices(1,2,4,5,3,6));
    
    Middle.reshape(Indices(this->d_fermi,this->d_fermi,this->d_fermi,this->d_fermi));
    Last.reshape(Indices(this->d_fermi,this->d_fermi,this->d_fermi,1));
  }  
  
  mpo.setOp(0,new Operator(First),true);
  mpo.setOp(this->N_fermions-1,new Operator(Last),true);
  if(this->N_fermions>2)
    mpo.setOp(1,new Operator(Middle),true);
  for(int i=2; i<this->N_fermions-1; i++)
    mpo.setOp(i,&mpo.getOp(1),false);  
#ifdef MATOUT
  mpo.exportForMatlab("ShiftMPO.m");
#endif
}

//Basically a wrapper for the more general function visible from the outside
void SU2HamiltonianNoLinks::getChargeConjugationMPO(MPO &mpo, bool project) const
{
  this->getChargeConjugationMPO(mpo,project,this->N_fermions);
}

void SU2HamiltonianNoLinks::getChargeConjugationMPO(MPO &mpo, bool project, int N) const
{
  mpo.initLength(N);
  
  //Cyclic Shift
  int d=this->d_fermi;
  int D=d*d;
  // The basic operator on any middle site
  mwArray basicOp=identityMatrix(d*d*d);
  basicOp.reshape(Indices(d,d,d,d,d,d));
  basicOp.permute(Indices(5,1,2,3,4,6));
  basicOp.reshape(Indices(d,d*d,d,d*d));
  // The operator on the left
  mwArray leftOp=identityMatrix(d*d);
  leftOp.reshape(Indices(d,1,d,d*d));
  mwArray rightOp=identityMatrix(d*d);
  rightOp.reshape(Indices(d,d,d,d));
  rightOp.permute(Indices(4,1,2,3));
  rightOp.reshape(Indices(d,d*d,d,1));  
  
  //Add projector if desired
  if(project)
  {
    mwArray Projector(Indices(d,d));
    Projector.fillWithZero();
    Projector.setElement(ONE_c,Indices(0,0));
    Projector.setElement(ONE_c,Indices(2,2));
    Projector.setElement(ONE_c,Indices(3,3));
    //Contract it with last tensor (first apply projector, then shift)
    rightOp.permute(Indices(1,2,4,3));
    rightOp.reshape(Indices(d*D*1,d));
    rightOp.multiplyRight(Projector);
    rightOp.reshape(Indices(d,D,1,d));
    rightOp.permute(Indices(1,2,4,3));    
  }
  
  //Add transformation of particles to antiparticles and vice versa on top
  mwArray Op(Indices(d,d));	Op.fillWithZero();    
  Op.setElement(ONE_c,Indices(0,3));
  Op.setElement(ONE_c,Indices(3,0));
  Op.setElement(ONE_c,Indices(1,1));
  Op.setElement(ONE_c,Indices(2,2));  
  
  leftOp.reshape(Indices(d,D*d));
  leftOp.multiplyLeft(Op);
  leftOp.reshape(Indices(d,1,d,D));
    
  basicOp.reshape(Indices(d,D*d*D));
  basicOp.multiplyLeft(Op);
  basicOp.reshape(Indices(d,D,d,D));
  
  rightOp.reshape(Indices(d,D*d));
  rightOp.multiplyLeft(Op);
  rightOp.reshape(Indices(d,D,d,1));
  
  mpo.setOp(0,new Operator(leftOp),true);
  mpo.setOp(N-1,new Operator(rightOp),true);
  if(N>2)
    mpo.setOp(1,new Operator(basicOp),true);
  for(int i=2; i<N-1; i++)
    mpo.setOp(i,&mpo.getOp(1),false); 
}

void SU2HamiltonianNoLinks::getCyclicShiftSquareMPO(MPO &mpo) const
{
  mpo.initLength(this->N_fermions);
//   int d=this->d_fermi;
//   int D=2*d*d;
//   // The basic operator on any middle site
//   mwArray OpOdd,OpEven,OpLeft,OpRight;
//   
//   OpOdd = mwArray(Indices(this->d_fermi,D,this->d_fermi,D));	OpOdd.fillWithZero();
//   OpEven = mwArray(Indices(this->d_fermi,D,this->d_fermi,D));	OpEven.fillWithZero();
//   
//   
//   for(int i=0; i<this->d_fermi; i++)
//     for(int j=0; j<this->d_fermi; j++)
//     {
//       for(int k=0; k<this->d_fermi; k++)
//       {
// 	OpOdd.setElement(ONE_c,Indices(i,k+i*this->d_fermi,j,k+j*this->d_fermi));
// 	
//     }
//     
// 
//   // Now set the operators in the chain
//   mpo.setOp(0,new Operator(leftOp),true);
//   mpo.setOp(this->N_fermions-1,new Operator(rightOp),true);
//   if(this->N_fermions>2)
//     mpo.setOp(1,new Operator(basicOp),true);
//   for(int k=2;k<this->N_fermions-1;k++)
//     mpo.setOp(k,&mpo.getOp(1),false);
#ifdef MATOUT
  mpo.exportForMatlab("CyclicShiftMPO.m");
#endif
}

complex_t SU2HamiltonianNoLinks::ChargeConjugationValue(const MPS &mps, Contractor &contractor, int offset) const
{
  complex_t result=ZERO_c;
  //Some error checking, do I have enough sites left, that after tracing out 2*offset sites still at least two are left?
  if(this->N_fermions-2*offset<2)
  {
    cout << "Error in SU2HamiltonianNoLinks::ChargeConjugationValue(), cannot trace away " << offset << " sites at each end for a system of size " << this->N_fermions << " because the resulting number of sites would be smaller than 2" << endl;
    exit(666);
  }
  //First get a copy of the MPS which I can manipulate and two states bra and ket for later
  MPS state(mps),ket1(this->N_fermions-2*offset,1,4),ket2(this->N_fermions-2*offset,1,4);
  
  //Properly gauge from each side, such that I don't have to compute the trace explicitly, it is simply an identity 
  for(int i=0; i<offset; i++)
  {
    state.gaugeCond(i,'R',true);
    state.gaugeCond(this->N_fermions-1-i,'L',true);
  }
  //cout << "Norm after gauging: " << contractor.contract(state,state) << endl;
  
  //Now prepare the states later used as bra and ket
  for(int i=offset; i<this->N_fermions-offset; i++)
  {
    //cout << "Reading site " << i << endl;
    ket1.replaceSite(i-offset,state.getA(i),false);
    ket2.replaceSite(i-offset,state.getA(i),false);
  }
  
  //Now contract the sites at the edges explicitly to check if they are really giving identities
#ifdef CHECKDIMS
  mwArray resLeft,resRight;
  resLeft=contractor.contract(state,state,offset-1,'L');
  resRight=contractor.contract(state,state,this->N_fermions-offset-1,'R');
  cout << "resLeft=" << resLeft << endl;
  cout << "resRight=" << resRight << endl;
  
  cout << "|ket1>=" << ket1 << endl;
  cout << "|ket2>=" << ket2 << endl;
#endif
  
  //Now get the charge conjugation operator and apply it in the middle of the chain
  MPO CC(1);
  this->getChargeConjugationMPO(CC,false,this->N_fermions-2*offset);  
  
  //contractor.contract(CC,ket1,ket2);
  
  
  //Now |ket1> is in the original state and |ket2>=CC|ket1>
  //cout << "|ket1>=" << ket1 << endl;
  //cout << "|ket2>=" << ket2 << endl;  
  
  //Now trace once again, but this time the remaining sites
  /*mwArray tmp,ket,bra,res;
  res = identityMatrix(ket1.getA(0).getDl());
  int d1,Dl1,Dr1,d2,Dl2,Dr2;
  res.reshape(Indices(1,-1));
  for(int i=0; i<this->N_fermions-2*offset; i++)
  {
    ket=ket1.getA(i).getA();
    bra=ket2.getA(i).getA();	bra.conjugate();
    
    d1=ket1.getA(i).getd();	Dl1=ket1.getA(i).getDl();	Dr1=ket1.getA(i).getDr();
    d2=ket2.getA(i).getd();	Dl2=ket2.getA(i).getDl();	Dr2=ket2.getA(i).getDr();
    //Now contract the physical index
    ket.permute(Indices(2,3,1));
    ket.reshape(Indices(Dl1*Dr1,d1));
    bra.reshape(Indices(d2,Dl2*Dr2));
    tmp = ket*bra;
    tmp.reshape(Indices(Dl1,Dr1,Dl2,Dr2));
    tmp.permute(Indices(1,3,2,4));
    tmp.reshape(Indices(Dl1*Dl2,Dr1*Dr2));
    res.multiplyRight(tmp);    
  }
  //Final contraction
  res.reshape(Indices(Dr1,Dr2));
  result = res.trace(); */
  
  //Now trace once again, but this time the remaining sites but this time hopefully more efficient
  /*mwArray ket,bra,res,op;
  res = identityMatrix(ket1.getA(0).getDl());
  int d1,Dl1,Dr1,d2,Dl2,Dr2,od1,od2,oDl,oDr;
  res.reshape(Indices(1,-1));
  for(int i=0; i<this->N_fermions-2*offset; i++)
  {
    ket=ket1.getA(i).getA();
    bra=ket2.getA(i).getA();	bra.conjugate();
    op=CC.getOpData(i);
    
    d1=ket.getDimension(0);	Dl1=ket.getDimension(1);	Dr1=ket.getDimension(2);
    d2=bra.getDimension(0);	Dl2=bra.getDimension(1);	Dr2=bra.getDimension(2);
    od1=op.getDimension(0);	oDl=op.getDimension(1);		od2=op.getDimension(2);		oDr=op.getDimension(3);
    
    //1. Contract bra with operator
    bra.permute(Indices(2,3,1));
    bra.reshape(Indices(Dl2*Dr2,d2));
    op.reshape(Indices(od1,oDl*od2*oDr));
    op.multiplyLeft(bra);
    cout << "Op after first stage: " << op.getDimensions() << endl;
    
    //2. Reshape and contract with ket tensor
    op.reshape(Indices(Dl2,Dr2,oDl,od2,oDr));
    op.permute(Indices(1,2,3,5,4));
    op.reshape(Indices(Dl2*Dr2*oDl*oDr,od2));
    ket.reshape(Indices(d1,Dl1*Dr1));
    op.multiplyRight(ket);
    
    cout << "Op after second stage: " << op.getDimensions() << endl;
    
    //3. Now take all bond indices from the left and the right together
    op.reshape(Indices(Dl2,Dr2,oDl,oDr,Dl1,Dr1));
    op.permute(Indices(1,3,5,2,4,6));
    op.reshape(Indices(Dl2*oDl*Dl1,Dr2*oDr*Dr1));
    
    cout << "Op after third stage: " << op.getDimensions() << endl;
    
    //4. Multiply with result
    res = res*op;   
    
    cout << "Result at end of step " << i << ": " << res.getDimensions() << endl << endl;
  }
  res.reshape(Indices(Dr2,Dr1));
  result=res.trace();*/
  
  //Now trace once again, but this time the remaining sites but this time hopefully more memory efficient
  mwArray ket,bra,res,op;
  res = identityMatrix(ket1.getA(0).getDl());
  int d1,Dl1,Dr1,d2,Dl2,Dr2,od1,od2,oDl,oDr;
  for(int i=0; i<this->N_fermions-2*offset; i++)
  {
    ket=ket1.getA(i).getA();
    bra=ket2.getA(i).getA();	bra.conjugate();
    op=CC.getOpData(i);
    
    d1=ket.getDimension(0);	Dl1=ket.getDimension(1);	Dr1=ket.getDimension(2);
    d2=bra.getDimension(0);	Dl2=bra.getDimension(1);	Dr2=bra.getDimension(2);
    od1=op.getDimension(0);	oDl=op.getDimension(1);		od2=op.getDimension(2);		oDr=op.getDimension(3);
    
    //1. Contract Dl of bra with result
    bra.permute(Indices(1,3,2));
    bra.reshape(Indices(d2*Dr2,Dl2));
    res.multiplyLeft(bra);    
    //cout << "res after first stage: " << res.getDimensions() << endl;
    
    //2. Reshape and contract with ket tensor
    res.reshape(Indices(d2*Dr2*oDl,Dl1));
    ket.permute(Indices(2,1,3));
    ket.reshape(Indices(Dl1,d1*Dr1));
    res.multiplyRight(ket);
    //cout << "res after second stage: " << res.getDimensions() << endl;
    
    //3. Now with the operator
    op.reshape(Indices(od1*oDl*od2,oDr));
    res.reshape(Indices(d2,Dr2,oDl,d1,Dr1));
    res.permute(Indices(2,5,1,3,4));
    res.reshape(Indices(Dr2*Dr1,d2*oDl*d1));
    res.multiplyRight(op);    
    //cout << "res after third stage: " << res.getDimensions() << endl;
    
    //5. Reshape and prepare for next round
    res.reshape(Indices(Dr2,Dr1,oDr));
    res.permute(Indices(1,3,2));
    res.reshape(Indices(Dr1,oDr*Dr1));    
    //cout << "Result at end of step " << i << ": " << res.getDimensions() << endl << endl;
  }
  res.reshape(Indices(Dr2,Dr1));
  result=res.trace();
  
  return result;
}

