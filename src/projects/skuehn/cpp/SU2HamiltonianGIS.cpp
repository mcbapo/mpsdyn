/**
   \file SU2HamiltonianGIS.cpp
   Implementation of the SU(2) Hamiltonian on the gauge invariant subspace after integrating out the gauge field as developed by Tao and Pablo

   \author Stefan Kühn
   \date 02/23/2018
*/

#include "SU2HamiltonianGIS.h"
#include "Indices.h"
#include "Contractor.h"
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>

using namespace std;
using namespace shrt;


//Standard Constructor
SU2HamiltonianGIS::SU2HamiltonianGIS(int N_,double epsilon_,double mu_,double g_, double lambda_):hamil(N_)
{
  //Number of sites
  this->N_sites=N_;
  //Hopping constant
  this->epsilon=epsilon_;  
  //Mass
  this->mu=mu_;
  //Coupling
  this->g=g_;
  //Penalty strength
  this->lambda = lambda_;
  //Dimension of a site
  this->d_site=4; 
  //Set the unused variables to  dummy values
  this->cext_left = -1;
  this->cext_right = -1;
  this->eta = -1.0;
  
  //Construct the basic operators
  initOperators();
  //Construct the actual hamiltonian MPO
  initHamiltonian();
}

//Constructor for external charge case
SU2HamiltonianGIS::SU2HamiltonianGIS(int N_,double epsilon_,double mu_,double g_,int k_, int l_, double lambda_,double eta_):hamil(N_+2)
{
  //Number of sites
  this->N_sites=N_;
  //Hopping constant
  this->epsilon=epsilon_;  
  //Mass
  this->mu=mu_;
  //Coupling
  this->g=g_;
  //Penalty strength
  this->lambda = lambda_;
  //Dimension of a site
  this->d_site=4;  
  //Penalty strength for single occupancy of external charges
  this->eta = eta_;
  //indices of the external charges
  if(k_<l_ && k_>0 && l_<N_-1)
  {
    this->cext_left = k_;
    this->cext_right = l_;
  }
  else if(k_>l_  && l_>0 && k_<N_-1)
  {
    this->cext_left = l_;
    this->cext_right = k_;
  }
  else
  {
    cout << "Error in SU2HamiltonianGIS::SU2HamiltonianGIS(), the external charges have to sit inside the chain at different sites" << endl;
    exit(1);
  }
  
  //Construct the basic operators
  initOperators();
  //Construct the actual hamiltonian MPO
  initHamiltonianExtCharges();
}


//Destructor
SU2HamiltonianGIS::~SU2HamiltonianGIS(){
  //Clear Hamiltonian MPO
  hamil.clear();  
  //Clear all operators
  id_fermi.clear();
  sigma_x.clear();
  sigma_y.clear();
  sigma_z.clear();
  sigma_plus.clear();
  sigma_minus.clear();

}

//Provide the basic operator matrices used for MPOs
void SU2HamiltonianGIS::initOperators(void)
{
  //First the fermionic operators (basically all kinds of Pauli Matrices)
  this->id_fermi = identityMatrix(2);
  
  this->id_site_mpo = identityMatrix(this->d_site);
  this->id_site_mpo.reshape(Indices(this->d_site,1,this->d_site,1));
  
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
   
}

//Build the Hamiltonian MPO
void SU2HamiltonianGIS::initHamiltonian()
{
  cout << "Building Hamiltonian MPO with" << endl
       << "epsilon: " << this->epsilon << endl
       << "mu:      " << this->mu << endl
       << "g:       " << this->g << endl
       << "lambda:  " << this->lambda << endl;
       
  //Initialize the operators acting on a site 
  mwArray Id,O1,O2,O3,O4,O1d,O2d,O3d,O4d,M1,M2,Qsquare,Qx,Qy,Qz;
  
  Id = identityMatrix(this->d_site);
  
  O1 = kron(this->sigma_plus,this->sigma_z);
  O2 = kron(this->sigma_minus,this->id_fermi);
  O3 = kron(this->id_fermi,this->sigma_plus);
  O4 = kron(this->sigma_z,this->sigma_minus);
  
  O1d = Hconjugate(O1);
  O2d = Hconjugate(O2);
  O3d = Hconjugate(O3);
  O4d = Hconjugate(O4);
  
  M1 = kron(this->sigma_z,this->id_fermi) + kron(this->id_fermi,this->id_fermi);
  M2 = kron(this->id_fermi,this->sigma_z) + kron(this->id_fermi,this->id_fermi);
  
  Qx = 0.5*I_c*(kron(this->sigma_minus,this->sigma_plus) - kron(this->sigma_plus,this->sigma_minus));
  Qy = -0.5*(kron(this->sigma_minus,this->sigma_plus) + kron(this->sigma_plus,this->sigma_minus));
  Qz = 0.25*(kron(this->sigma_z,this->id_fermi) - kron(this->id_fermi,this->sigma_z));
  
  /*Check if Q^alpha squared is the same for all components
  cout << "Qx^2 = " << Qx*Qx  << endl;
  cout << "Qy^2 = " << Qy*Qy  << endl;
  cout << "Qz^2 = " << Qz*Qz  << endl;*/
  
  Qsquare = Qx*Qx;
  
  
  //Array to store the fermionic operators
  int N_ops=15;
  mwArray Fermionic_operators(Indices(N_ops,this->d_site,this->d_site));
  Fermionic_operators.fillWithZero();
  complex_t constant;
  for(int i=0; i<this->d_site; ++i)
  {
    for(int j=0; j<this->d_site;++j)
    {
      //Identity
      Fermionic_operators.setElement(Id.getElement(Indices(i,j)),Indices(0,i,j));
      //O1 to O4 and hermitian conjugates thereof
      Fermionic_operators.setElement(O1.getElement(Indices(i,j)),Indices(1,i,j));
      Fermionic_operators.setElement(O2.getElement(Indices(i,j)),Indices(2,i,j));
      Fermionic_operators.setElement(O3.getElement(Indices(i,j)),Indices(3,i,j));
      Fermionic_operators.setElement(O4.getElement(Indices(i,j)),Indices(4,i,j));
      Fermionic_operators.setElement(O1d.getElement(Indices(i,j)),Indices(5,i,j));
      Fermionic_operators.setElement(O2d.getElement(Indices(i,j)),Indices(6,i,j));
      Fermionic_operators.setElement(O3d.getElement(Indices(i,j)),Indices(7,i,j));
      Fermionic_operators.setElement(O4d.getElement(Indices(i,j)),Indices(8,i,j));
      //M1 and M2
      Fermionic_operators.setElement(M1.getElement(Indices(i,j)),Indices(9,i,j));
      Fermionic_operators.setElement(M2.getElement(Indices(i,j)),Indices(10,i,j));
      //Q^2 and components thereof
      Fermionic_operators.setElement(Qsquare.getElement(Indices(i,j)),Indices(11,i,j));
      Fermionic_operators.setElement(Qx.getElement(Indices(i,j)),Indices(12,i,j));
      Fermionic_operators.setElement(Qy.getElement(Indices(i,j)),Indices(13,i,j));
      Fermionic_operators.setElement(Qz.getElement(Indices(i,j)),Indices(14,i,j));
      
    }
  }  
  Fermionic_operators.reshape(Indices(N_ops,this->d_site*this->d_site));
  //cout << "Fermionic operators done" << endl;
  
  //Now prepare the matrices
  int D_bond=9;
  mwArray Site;
  
  for(int n=0; n<this->N_sites; n++)
  {
    if(n==0)
    {
      //First site
      Site = mwArray(Indices(1,D_bond,N_ops));		Site.fillWithZero();
      
      //Identity
      Site.setElement(ONE_c,Indices(0,0,0));
      
      //Hopping part
      Site.setElement(this->epsilon*ONE_c,Indices(0,1,1));
      Site.setElement(this->epsilon*ONE_c,Indices(0,2,3));
      Site.setElement(this->epsilon*ONE_c,Indices(0,3,5));
      Site.setElement(this->epsilon*ONE_c,Indices(0,4,7));
      //Mass terms
      constant = pow(-1.0,(double) n)*this->mu/2.0*ONE_c;
      Site.setElement(constant,Indices(0,D_bond-1,9));
      Site.setElement(constant,Indices(0,D_bond-1,10));
      //Electric part and penalty 
      constant = 3.0*(this->g*this->g/2.0*(this->N_sites-1.0-n)+this->lambda)*ONE_c;
      Site.setElement(constant,Indices(0,D_bond-1,11));
      Site.setElement(ONE_c,Indices(0,5,12));
      Site.setElement(ONE_c,Indices(0,6,13));
      Site.setElement(ONE_c,Indices(0,7,14));
            
      Site.reshape(Indices(1*D_bond,N_ops));
      Site.multiplyRight(Fermionic_operators);
      Site.reshape(Indices(1,D_bond,this->d_site,this->d_site));
      Site.permute(Indices(3,1,4,2));
    }
    else if(n==this->N_sites-1)
    {
      //Last site
      Site = mwArray(Indices(D_bond,1,N_ops));		Site.fillWithZero();
      
      //Identity
      Site.setElement(ONE_c,Indices(D_bond-1,0,0));
      
      //Hopping part
      Site.setElement(ONE_c,Indices(1,0,2));
      Site.setElement(ONE_c,Indices(2,0,4));
      Site.setElement(ONE_c,Indices(3,0,6));
      Site.setElement(ONE_c,Indices(4,0,8));
      //Mass terms
      constant = pow(-1.0,(double) n)*this->mu/2.0*ONE_c;
      Site.setElement(constant,Indices(0,0,9));
      Site.setElement(constant,Indices(0,0,10));
      //Electric part and penalty 
      constant = 3.0*this->lambda*ONE_c;
      Site.setElement(constant,Indices(0,0,11));
      constant = 2.0*this->lambda*ONE_c;
      Site.setElement(constant,Indices(5,0,12));
      Site.setElement(constant,Indices(6,0,13));
      Site.setElement(constant,Indices(7,0,14));
      
      Site.reshape(Indices(D_bond*1,N_ops));
      Site.multiplyRight(Fermionic_operators);
      Site.reshape(Indices(D_bond,1,this->d_site,this->d_site));
      Site.permute(Indices(3,1,4,2));
    }
    else
    {
      //Sites in the center
      Site = mwArray(Indices(D_bond,D_bond,N_ops));	Site.fillWithZero();
      
      //Identity
      Site.setElement(ONE_c,Indices(0,0,0));
      Site.setElement(ONE_c,Indices(D_bond-1,D_bond-1,0));
      Site.setElement(ONE_c,Indices(5,5,0));
      Site.setElement(ONE_c,Indices(6,6,0));
      Site.setElement(ONE_c,Indices(7,7,0));
      
      //Hopping part
      Site.setElement(this->epsilon*ONE_c,Indices(0,1,1));
      Site.setElement(ONE_c,Indices(1,D_bond-1,2));
      Site.setElement(this->epsilon*ONE_c,Indices(0,2,3));
      Site.setElement(ONE_c,Indices(2,D_bond-1,4));
      Site.setElement(this->epsilon*ONE_c,Indices(0,3,5));
      Site.setElement(ONE_c,Indices(3,D_bond-1,6));
      Site.setElement(this->epsilon*ONE_c,Indices(0,4,7));
      Site.setElement(ONE_c,Indices(4,D_bond-1,8));
      //Mass terms
      constant = pow(-1.0,(double) n)*this->mu/2.0*ONE_c;
      Site.setElement(constant,Indices(0,D_bond-1,9));
      Site.setElement(constant,Indices(0,D_bond-1,10));
      //Electric part and penalty 
      constant = 3.0*(this->g*this->g/2.0*(this->N_sites-1.0-n)+this->lambda)*ONE_c;
      Site.setElement(constant,Indices(0,D_bond-1,11));
      constant = (this->g*this->g*(this->N_sites-1.0-n)+2.0*this->lambda)*ONE_c;
      Site.setElement(ONE_c,Indices(0,5,12));
      Site.setElement(constant,Indices(5,D_bond-1,12));
      Site.setElement(ONE_c,Indices(0,6,13));
      Site.setElement(constant,Indices(6,D_bond-1,13));
      Site.setElement(ONE_c,Indices(0,7,14));
      Site.setElement(constant,Indices(7,D_bond-1,14));    
      
      Site.reshape(Indices(D_bond*D_bond,N_ops));
      Site.multiplyRight(Fermionic_operators);
      Site.reshape(Indices(D_bond,D_bond,this->d_site,this->d_site));
      Site.permute(Indices(3,1,4,2));      
    }
    
    hamil.setOp(n,new Operator(Site),true);
  }
#ifdef MATOUT
  ostringstream convert;  
  string mponame;
  //Build the filename
  convert << this->N_sites;
  mponame = "HSU2_"+mponame+"N"+convert.str();
  convert.str("");
  convert.clear();
  convert << this->epsilon;
  mponame = mponame+"_t"+convert.str();
  convert.str("");
  convert.clear();
  convert << this->mu;
  mponame = mponame+"_m"+convert.str();
  convert.str("");
  convert.clear();
  convert << this->g;
  mponame = mponame+"_g"+convert.str();
  convert.str("");
  convert.clear();
  convert << this->lambda;
  mponame = mponame+"_la"+convert.str()+".m";
  convert.str("");
  convert.clear();
  hamil.exportForMatlab(mponame.c_str(),16);
#endif
}

//Build the Hamiltonian MPO
void SU2HamiltonianGIS::initHamiltonianExtCharges()
{
  cout << "Building Hamiltonian MPO with external charges and" << endl
       << "epsilon:    " << this->epsilon << endl
       << "mu:         " << this->mu << endl
       << "g:          " << this->g << endl
       << "lambda:     " << this->lambda << endl
       << "eta:        " << this->eta << endl
       << "cext_left:  " << this->cext_left << endl
       << "cext_right: " << this->cext_right << endl;
       
  //Initialize the operators acting on a site 
  mwArray Id,O1,O2,O3,O4,O1d,O2d,O3d,O4d,M1,M2,Qsquare,Qx,Qy,Qz,P;
  
  Id = identityMatrix(this->d_site);
  
  O1 = kron(this->sigma_plus,this->sigma_z);
  O2 = kron(this->sigma_minus,this->id_fermi);
  O3 = kron(this->id_fermi,this->sigma_plus);
  O4 = kron(this->sigma_z,this->sigma_minus);
  
  O1d = Hconjugate(O1);
  O2d = Hconjugate(O2);
  O3d = Hconjugate(O3);
  O4d = Hconjugate(O4);
  
  M1 = kron(this->sigma_z,this->id_fermi) + kron(this->id_fermi,this->id_fermi);
  M2 = kron(this->id_fermi,this->sigma_z) + kron(this->id_fermi,this->id_fermi);
  
  Qx = 0.5*I_c*(kron(this->sigma_minus,this->sigma_plus) - kron(this->sigma_plus,this->sigma_minus));
  Qy = -0.5*(kron(this->sigma_minus,this->sigma_plus) + kron(this->sigma_plus,this->sigma_minus));
  Qz = 0.25*(kron(this->sigma_z,this->id_fermi) - kron(this->id_fermi,this->sigma_z));

  Qsquare = Qx*Qx;
  
  P = identityMatrix(this->d_site) + kron(this->sigma_z,this->sigma_z);
  
  //Array to store the fermionic operators
  int N_ops=16;
  mwArray Fermionic_operators(Indices(N_ops,this->d_site,this->d_site));
  Fermionic_operators.fillWithZero();
  complex_t constant;
  for(int i=0; i<this->d_site; ++i)
  {
    for(int j=0; j<this->d_site;++j)
    {
      //Identity
      Fermionic_operators.setElement(Id.getElement(Indices(i,j)),Indices(0,i,j));
      //O1 to O4 and hermitian conjugates thereof
      Fermionic_operators.setElement(O1.getElement(Indices(i,j)),Indices(1,i,j));
      Fermionic_operators.setElement(O2.getElement(Indices(i,j)),Indices(2,i,j));
      Fermionic_operators.setElement(O3.getElement(Indices(i,j)),Indices(3,i,j));
      Fermionic_operators.setElement(O4.getElement(Indices(i,j)),Indices(4,i,j));
      Fermionic_operators.setElement(O1d.getElement(Indices(i,j)),Indices(5,i,j));
      Fermionic_operators.setElement(O2d.getElement(Indices(i,j)),Indices(6,i,j));
      Fermionic_operators.setElement(O3d.getElement(Indices(i,j)),Indices(7,i,j));
      Fermionic_operators.setElement(O4d.getElement(Indices(i,j)),Indices(8,i,j));
      //M1 and M2
      Fermionic_operators.setElement(M1.getElement(Indices(i,j)),Indices(9,i,j));
      Fermionic_operators.setElement(M2.getElement(Indices(i,j)),Indices(10,i,j));
      //Q^2 and components thereof
      Fermionic_operators.setElement(Qsquare.getElement(Indices(i,j)),Indices(11,i,j));
      Fermionic_operators.setElement(Qx.getElement(Indices(i,j)),Indices(12,i,j));
      Fermionic_operators.setElement(Qy.getElement(Indices(i,j)),Indices(13,i,j));
      Fermionic_operators.setElement(Qz.getElement(Indices(i,j)),Indices(14,i,j));
      //Penalty
      Fermionic_operators.setElement(P.getElement(Indices(i,j)),Indices(15,i,j));      
    }
  }  
  Fermionic_operators.reshape(Indices(N_ops,this->d_site*this->d_site));
  //cout << "Fermionic operators done" << endl;
  
  //Now prepare the matrices
  int D_bond=18;
  mwArray Site;
  
  for(int n=0; n<this->N_sites; n++)
  {
    if(n==0)
    {
      //First site
      Site = mwArray(Indices(1,D_bond,N_ops));		Site.fillWithZero();
      
      //Identity
      Site.setElement(ONE_c,Indices(0,0,0));
      
      //Hopping part
      Site.setElement(this->epsilon*ONE_c,Indices(0,1,1));
      Site.setElement(this->epsilon*ONE_c,Indices(0,2,3));
      Site.setElement(this->epsilon*ONE_c,Indices(0,3,5));
      Site.setElement(this->epsilon*ONE_c,Indices(0,4,7));
      //Mass terms
      constant = pow(-1.0,(double) n)*this->mu/2.0*ONE_c;
      Site.setElement(constant,Indices(0,D_bond-1,9));
      Site.setElement(constant,Indices(0,D_bond-1,10));
      //Electric part and penalty 
      constant = 3.0*(this->g*this->g/2.0*(this->N_sites-1.0-n)+this->lambda)*ONE_c;
      Site.setElement(constant,Indices(0,D_bond-1,11));
      Site.setElement(ONE_c,Indices(0,5,12));
      Site.setElement(ONE_c,Indices(0,6,13));
      Site.setElement(ONE_c,Indices(0,7,14));
      //Electric part and penalty corresponding to external charge part
      constant = (this->g*this->g*(this->N_sites-1.0-max(n,this->cext_left))+2.0*this->lambda)*ONE_c;
      Site.setElement(constant,Indices(0,8,12));
      Site.setElement(constant,Indices(0,9,13));
      Site.setElement(constant,Indices(0,10,14));
      constant = (this->g*this->g*(this->N_sites-1.0-max(n,this->cext_right))+2.0*this->lambda)*ONE_c;
      Site.setElement(constant,Indices(0,11,12));
      Site.setElement(constant,Indices(0,12,13));
      Site.setElement(constant,Indices(0,13,14));      
            
      Site.reshape(Indices(1*D_bond,N_ops));
      Site.multiplyRight(Fermionic_operators);
      Site.reshape(Indices(1,D_bond,this->d_site,this->d_site));
      Site.permute(Indices(3,1,4,2));
    }
    else if(n==this->N_sites-1)
    {
      //Last site
      Site = mwArray(Indices(D_bond,D_bond,N_ops));	Site.fillWithZero();
      
      //Identity
      Site.setElement(ONE_c,Indices(0,0,0));
      Site.setElement(ONE_c,Indices(D_bond-1,D_bond-1,0));
      Site.setElement(ONE_c,Indices(8,8,0));
      Site.setElement(ONE_c,Indices(9,9,0));
      Site.setElement(ONE_c,Indices(10,10,0));
      Site.setElement(ONE_c,Indices(11,11,0));
      Site.setElement(ONE_c,Indices(12,12,0));
      Site.setElement(ONE_c,Indices(13,13,0)); 
      
      //Hopping part
      Site.setElement(ONE_c,Indices(1,D_bond-1,2));
      Site.setElement(ONE_c,Indices(2,D_bond-1,4));
      Site.setElement(ONE_c,Indices(3,D_bond-1,6));
      Site.setElement(ONE_c,Indices(4,D_bond-1,8));
      //Mass terms
      constant = pow(-1.0,(double) n)*this->mu/2.0*ONE_c;
      Site.setElement(constant,Indices(0,D_bond-1,9));
      Site.setElement(constant,Indices(0,D_bond-1,10));
      //Electric part and penalty 
      constant = 3.0*this->lambda*ONE_c;
      Site.setElement(constant,Indices(0,D_bond-1,11));
      constant = 2.0*this->lambda*ONE_c;
      Site.setElement(constant,Indices(5,D_bond-1,12));
      Site.setElement(constant,Indices(6,D_bond-1,13));
      Site.setElement(constant,Indices(7,D_bond-1,14));    
      //Electric part and penalty corresponding to external charge part
      constant = 2.0*this->lambda*ONE_c;
      Site.setElement(constant,Indices(0,8,12));
      Site.setElement(constant,Indices(0,9,13));
      Site.setElement(constant,Indices(0,10,14));
      constant = 2.0*this->lambda*ONE_c;
      Site.setElement(constant,Indices(0,11,12));
      Site.setElement(constant,Indices(0,12,13));
      Site.setElement(constant,Indices(0,13,14));
      
      Site.reshape(Indices(D_bond*D_bond,N_ops));
      Site.multiplyRight(Fermionic_operators);
      Site.reshape(Indices(D_bond,D_bond,this->d_site,this->d_site));
      Site.permute(Indices(3,1,4,2));
    }
    else
    {
      //Sites in the center
      Site = mwArray(Indices(D_bond,D_bond,N_ops));	Site.fillWithZero();
      
      //Identity
      Site.setElement(ONE_c,Indices(0,0,0));
      Site.setElement(ONE_c,Indices(D_bond-1,D_bond-1,0));
      Site.setElement(ONE_c,Indices(5,5,0));
      Site.setElement(ONE_c,Indices(6,6,0));
      Site.setElement(ONE_c,Indices(7,7,0));
      Site.setElement(ONE_c,Indices(8,8,0));
      Site.setElement(ONE_c,Indices(9,9,0));
      Site.setElement(ONE_c,Indices(10,10,0));
      Site.setElement(ONE_c,Indices(11,11,0));
      Site.setElement(ONE_c,Indices(12,12,0));
      Site.setElement(ONE_c,Indices(13,13,0)); 
      
      //Hopping part
      Site.setElement(this->epsilon*ONE_c,Indices(0,1,1));
      Site.setElement(ONE_c,Indices(1,D_bond-1,2));
      Site.setElement(this->epsilon*ONE_c,Indices(0,2,3));
      Site.setElement(ONE_c,Indices(2,D_bond-1,4));
      Site.setElement(this->epsilon*ONE_c,Indices(0,3,5));
      Site.setElement(ONE_c,Indices(3,D_bond-1,6));
      Site.setElement(this->epsilon*ONE_c,Indices(0,4,7));
      Site.setElement(ONE_c,Indices(4,D_bond-1,8));
      //Mass terms
      constant = pow(-1.0,(double) n)*this->mu/2.0*ONE_c;
      Site.setElement(constant,Indices(0,D_bond-1,9));
      Site.setElement(constant,Indices(0,D_bond-1,10));
      //Electric part and penalty 
      constant = 3.0*(this->g*this->g/2.0*(this->N_sites-1.0-n)+this->lambda)*ONE_c;
      Site.setElement(constant,Indices(0,D_bond-1,11));
      constant = (this->g*this->g*(this->N_sites-1.0-n)+2.0*this->lambda)*ONE_c;
      Site.setElement(ONE_c,Indices(0,5,12));
      Site.setElement(constant,Indices(5,D_bond-1,12));
      Site.setElement(ONE_c,Indices(0,6,13));
      Site.setElement(constant,Indices(6,D_bond-1,13));
      Site.setElement(ONE_c,Indices(0,7,14));
      Site.setElement(constant,Indices(7,D_bond-1,14));    
      //Electric part and penalty corresponding to external charge part
      constant = (this->g*this->g*(this->N_sites-1.0-max(n,this->cext_left))+2.0*this->lambda)*ONE_c;
      Site.setElement(constant,Indices(0,8,12));
      Site.setElement(constant,Indices(0,9,13));
      Site.setElement(constant,Indices(0,10,14));
      constant = (this->g*this->g*(this->N_sites-1.0-max(n,this->cext_right))+2.0*this->lambda)*ONE_c;
      Site.setElement(constant,Indices(0,11,12));
      Site.setElement(constant,Indices(0,12,13));
      Site.setElement(constant,Indices(0,13,14));            
      
      Site.reshape(Indices(D_bond*D_bond,N_ops));
      Site.multiplyRight(Fermionic_operators);
      Site.reshape(Indices(D_bond,D_bond,this->d_site,this->d_site));
      Site.permute(Indices(3,1,4,2));      
    }
    
    hamil.setOp(n,new Operator(Site),true);
  }
  
  //Now the second last site
  Site = mwArray(Indices(D_bond,D_bond,N_ops));	Site.fillWithZero();
      
  //Identity
  Site.setElement(ONE_c,Indices(0,0,0));
  Site.setElement(ONE_c,Indices(D_bond-1,D_bond-1,0));
  Site.setElement(ONE_c,Indices(11,11,0));
  Site.setElement(ONE_c,Indices(12,12,0));
  Site.setElement(ONE_c,Indices(13,13,0)); 
  
  //Electric part and penalty corresponding to external charge part
  constant = 3.0*(this->g*this->g/2.0*(this->N_sites-1.0-this->cext_left)+this->lambda)*ONE_c;
  Site.setElement(constant,Indices(0,D_bond-1,11));
  constant = (this->g*this->g*(this->N_sites-1.0-this->cext_right)+2.0*this->lambda)*ONE_c;
  Site.setElement(constant*ONE_c,Indices(0,14,12));
  Site.setElement(ONE_c,Indices(8,D_bond-1,12));
  Site.setElement(constant*ONE_c,Indices(0,15,13));
  Site.setElement(ONE_c,Indices(9,D_bond-1,13));
  Site.setElement(constant*ONE_c,Indices(0,16,14));
  Site.setElement(ONE_c,Indices(10,D_bond-1,14));
  
  //Penalty for single occupancy
  Site.setElement(this->eta*ONE_c,Indices(0,D_bond-1,15));
  
  Site.reshape(Indices(D_bond*D_bond,N_ops));
  Site.multiplyRight(Fermionic_operators);
  Site.reshape(Indices(D_bond,D_bond,this->d_site,this->d_site));
  Site.permute(Indices(3,1,4,2));
  
  hamil.setOp(this->N_sites,new Operator(Site),true);
  
  //Now the last site
  Site = mwArray(Indices(D_bond,1,N_ops));	Site.fillWithZero();
      
  //Identity
  Site.setElement(ONE_c,Indices(D_bond-1,0,0));
  
  //Electric part and penalty corresponding to external charge part
  constant = 3.0*(this->g*this->g/2.0*(this->N_sites-1.0-this->cext_right)+this->lambda)*ONE_c;
  Site.setElement(constant,Indices(0,0,11));
  Site.setElement(ONE_c,Indices(11,0,12));
  Site.setElement(ONE_c,Indices(14,0,12));
  Site.setElement(ONE_c,Indices(12,0,13));
  Site.setElement(ONE_c,Indices(15,0,13));
  Site.setElement(ONE_c,Indices(13,0,14));
  Site.setElement(ONE_c,Indices(16,0,14));
  
  //Penalty for single occupancy
  Site.setElement(this->eta*ONE_c,Indices(0,0,15));
  
  Site.reshape(Indices(D_bond*1,N_ops));
  Site.multiplyRight(Fermionic_operators);
  Site.reshape(Indices(D_bond,1,this->d_site,this->d_site));
  Site.permute(Indices(3,1,4,2));
  
  hamil.setOp(this->N_sites+1,new Operator(Site),true);
 
#ifdef MATOUT
  ostringstream convert;  
  string mponame;
  //Build the filename
  convert << this->N_sites;
  mponame = "HSU2ExtCharges_"+mponame+"N"+convert.str();
  convert.str("");
  convert.clear();
  convert << this->epsilon;
  mponame = mponame+"_t"+convert.str();
  convert.str("");
  convert.clear();
  convert << this->mu;
  mponame = mponame+"_m"+convert.str();
  convert.str("");
  convert.clear();  
  convert << this->g;
  mponame = mponame+"_g"+convert.str();
  convert.str("");
  convert.clear();  
  convert << this->cext_left;
  mponame = mponame+"_cleft"+convert.str();
  convert.str("");
  convert.clear();  
  convert << this->cext_right;
  mponame = mponame+"_cright"+convert.str();
  convert.str("");
  convert.clear();  
  convert << this->lambda;
  mponame = mponame+"_la"+convert.str();
  convert.str("");
  convert.clear();
  convert << this->eta;
  mponame = mponame+"_eta"+convert.str()+".m";
  hamil.exportForMatlab(mponame.c_str(),16);
#endif
}

void SU2HamiltonianGIS::getQsquareMPO(MPO &mpo, int site) const
{   
  this->getSingleBodyMPO(mpo,Q2Op,site);
}

void SU2HamiltonianGIS::getQalphaMPO(MPO &mpo, int site,Component alpha) const
{
  if(alpha==x_comp)
    this->getSingleBodyMPO(mpo,QxOp,site);
  else if(alpha==y_comp)
    this->getSingleBodyMPO(mpo,QyOp,site);
  else if(alpha==z_comp)
    this->getSingleBodyMPO(mpo,QzOp,site);
  else
  {
    cout << "Error in SU2HamiltonianGIS::getQalphaMPO, unknown component" << endl;
    exit(-1);
  }
}

void SU2HamiltonianGIS::getSingleBodyMPO(MPO &mpo,SingleBodyOperator Op, int site) const
{
  //Check if the site index is valid
  if(site<0 || site>=this->N_sites)
  {
    cout << "Error in SU2HamiltonianGIS::getSingleBodyMPO, received an invalid site" << endl;
    exit(-1);
  }
  
  //Initialize the MPO 
  mpo.initLength(this->N_sites);
  
  //Retrieve the operator
  mwArray M;
  if(Op == QxOp)
    M = 0.5*I_c*(kron(this->sigma_minus,this->sigma_plus) - kron(this->sigma_plus,this->sigma_minus));
  else if(Op == QyOp)
    M = -0.5*(kron(this->sigma_minus,this->sigma_plus) + kron(this->sigma_plus,this->sigma_minus));
  else if(Op == QzOp)
    M = 0.25*(kron(this->sigma_z,this->id_fermi) - kron(this->id_fermi,this->sigma_z));
  else if(Op == Q2Op)
  {
    M = 0.25*(kron(this->sigma_z,this->id_fermi) - kron(this->id_fermi,this->sigma_z));
    M = M*M;
  }
  cout << "Operator = " << M << endl;
  M.reshape(Indices(this->d_site,1,this->d_site,1));
  
  mpo.setOp(site,new Operator(M),true);
  
  int pos_saved;
  bool saved=false;
  for(int i=0;i<this->N_sites; i++)
  {
    if(saved && i!=site)
      mpo.setOp(i,&mpo.getOp(pos_saved),false);
    else if(!saved && i!=site)
    {
      mpo.setOp(i,new Operator(this->id_site_mpo),true);
      pos_saved = i;
      saved = true;
    }
  }    
}

void SU2HamiltonianGIS::getExtChargeCorrelationMPO(MPO &mpo,Component alpha) const
{
  mwArray cOp;
  
  //Select appropriate operator
  if(alpha==x_comp)
    cOp = 0.5*I_c*(kron(this->sigma_minus,this->sigma_plus) - kron(this->sigma_plus,this->sigma_minus));
  else if(alpha==y_comp)
    cOp = -0.5*(kron(this->sigma_minus,this->sigma_plus) + kron(this->sigma_plus,this->sigma_minus));
  else if(alpha==z_comp)
    cOp = 0.25*(kron(this->sigma_z,this->id_fermi) - kron(this->id_fermi,this->sigma_z));
  else
  {
    cout << "Warning in SU2HamiltonianGIS::getExtChargeCorrelationMPO(), received unknown component, will take x-component" << endl;
    cOp = 0.5*I_c*(kron(this->sigma_minus,this->sigma_plus) - kron(this->sigma_plus,this->sigma_minus));
  }
  cOp.reshape(Indices(this->d_site,1,this->d_site,1));
  
  mpo.initLength(this->N_sites+2);
  mpo.setOp(0,new Operator(this->id_site_mpo),true);
  for (int i=1; i<this->N_sites; i++)
    mpo.setOp(i,&mpo.getOp(0),false);
  
  mpo.setOp(this->N_sites,new Operator(cOp),true);
  mpo.setOp(this->N_sites+1,&mpo.getOp(this->N_sites),false);
}