/**
   \file SU2Hamiltonian.cpp
   Implementation of the truncated version of the SU(2) Hamiltonian from Erez, Benni and Ignacio
   
   \author Stefan Kühn
   \date 12/08/2014
*/

#include "SU2Hamiltonian.h"
#include "Indices.h"
#include "Contractor.h"
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>

using namespace std;
using namespace shrt;

#include "buildProjectorMatrices.h"

//Standard Constructor
SU2Hamiltonian::SU2Hamiltonian(int N_,double epsilon_,double mu_,double g_):hamil(2*N_+N_-1)
{
  //Number of fermions
  this->N_fermions=N_;
  //Total number of sites (2*N_fermions-1)
  this->N_total=2*N_+N_-1;  
  //Hopping constant
  this->epsilon=epsilon_;  
  //Mass
  this->mu=mu_;
  //Coupling
  this->g=g_;
  //Dimension of the link Hilbert spaces
  this->d_link=5;
  //Dimension of the fermionic Hilbertspaces
  this->d_fermi=2;
  //Correlation MPO matrices are not yet saved
  this->correlation_matrices_saved=false;
  
  //Construct the basic operators
  initOperators();
  //Construct the actual hamiltonian MPO
  initHamiltonian();
}

//Constructor for site depended mass
SU2Hamiltonian::SU2Hamiltonian(int N_,double epsilon_,vector<double> mu_,double g_):hamil(2*N_+N_-1)
{
  //First check, of probably only a single mass is given, in this case I assume, that probably the other version of the constructor was meant to be called (as it is more efficient). Therefore I call it manually and put a warning
  if(mu_.size()==1)
  {
    cout << "WARNING, as only a single mass is given, I assume a constant mass of " << endl;
    //Number of fermions
    this->N_fermions=N_;
    //Total number of sites (2*N_fermions-1)
    this->N_total=2*N_+N_-1;  
    //Hopping constant
    this->epsilon=epsilon_;  
    //Mass
    this->mu=mu_[0];
    //Coupling
    this->g=g_;
    //Dimension of the link Hilbert spaces
    this->d_link=5;
    //Dimension of the fermionic Hilbertspaces
    this->d_fermi=2;
    //Correlation MPO matrices are not yet saved
    this->correlation_matrices_saved=false;
    
    //Construct the basic operators
    initOperators();
    //Construct the actual hamiltonian MPO
    initHamiltonian();
  }
  else
  {
    //Number of fermions
    this->N_fermions=N_;
    //Total number of sites (2*N_fermions-1)
    this->N_total=2*N_+N_-1;  
    //Hopping constant
    this->epsilon=epsilon_;  
    //Mass
    if(mu_.size()==N_)
      this->mu_site=mu_;
    else
    {
      cout << "Error in SU2Hamiltonian::SU2Hamiltonian(), number of masses does not correspond to number of sites, program will be aborted" << endl;
      exit(666);
    }
    //Coupling
    this->g=g_;
    //Dimension of the link Hilbert spaces
    this->d_link=5;
    //Dimension of the fermionic Hilbertspaces
    this->d_fermi=2;
    //Correlation MPO matrices are not yet saved
    this->correlation_matrices_saved=false;
    
    //Construct the basic operators
    initOperators();
    //Construct the actual hamiltonian MPO with site dependend masses
    initHamiltonian_sitedependendmass();
  }
}

//Constructor for site depended mass and charges at the specified sites
SU2Hamiltonian::SU2Hamiltonian(int N_,double epsilon_,vector<double> mu_,double g_, vector<int> charge_pos_, double lambda):hamil(2*N_+N_-1)
{
  //Number of fermions
  this->N_fermions=N_;
  //Total number of sites (2*N_fermions-1)
  this->N_total=2*N_+N_-1;  
  //Hopping constant
  this->epsilon=epsilon_;  
  //Mass
  if(mu_.size()==N_)
    this->mu_site=mu_;
  else
  {
    cout << "Error in SU2Hamiltonian::SU2Hamiltonian(), number of masses does not correspond to number of sites, program will be aborted" << endl;
    exit(666);
  }
  //Coupling
  this->g=g_;
  //Dimension of the link Hilbert spaces
  this->d_link=5;
  //Dimension of the fermionic Hilbertspaces
  this->d_fermi=2;
  
  //Construct the basic operators
  initOperators();
  //Construct the actual hamiltonian MPO with site depended masses and charges given at the sites specified in charge_pos_, where lambda is the strength of the penalty
  initHamiltonian_sitedependendmass_charge(charge_pos_,lambda);
}

//Destructor
SU2Hamiltonian::~SU2Hamiltonian(){
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

}

//Provide the basic operator matrices used for MPOs
void SU2Hamiltonian::initOperators(void)
{
  //First the fermionic operators (basically all kinds of Pauli Matrices)
  this->id_fermi=identityMatrix(2);
  
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
  
  //Left generators
  this->Lx=mwArray(Indices(5,5));	this->Lx.fillWithZero();
  this->Ly=mwArray(Indices(5,5));	this->Ly.fillWithZero();
  this->Lz=mwArray(Indices(5,5));	this->Lz.fillWithZero();
  
  this->Lx.setElement(0.5*ONE_c,Indices(0,3));	this->Ly.setElement(-0.5*I_c,Indices(0,3));	this->Lz.setElement(-0.5*ONE_c,Indices(0,0));
  this->Lx.setElement(0.5*ONE_c,Indices(1,4));	this->Ly.setElement(-0.5*I_c,Indices(1,4));	this->Lz.setElement(-0.5*ONE_c,Indices(1,1));
  this->Lx.setElement(0.5*ONE_c,Indices(3,0));	this->Ly.setElement(0.5*I_c,Indices(3,0));	this->Lz.setElement(0.5*ONE_c,Indices(3,3));
  this->Lx.setElement(0.5*ONE_c,Indices(4,1));	this->Ly.setElement(0.5*I_c,Indices(4,1));	this->Lz.setElement(0.5*ONE_c,Indices(4,4));
  
  //Right generators
  this->Rx=mwArray(Indices(5,5));	this->Rx.fillWithZero();
  this->Ry=mwArray(Indices(5,5));	this->Ry.fillWithZero();
  this->Rz=mwArray(Indices(5,5));	this->Rz.fillWithZero();
  
  this->Rx.setElement(0.5*ONE_c,Indices(0,1));	this->Ry.setElement(0.5*I_c,Indices(0,1));	this->Rz.setElement(-0.5*ONE_c,Indices(0,0));
  this->Rx.setElement(0.5*ONE_c,Indices(1,0));	this->Ry.setElement(-0.5*I_c,Indices(1,0));	this->Rz.setElement(0.5*ONE_c,Indices(1,1));
  this->Rx.setElement(0.5*ONE_c,Indices(3,4));	this->Ry.setElement(0.5*I_c,Indices(3,4));	this->Rz.setElement(-0.5*ONE_c,Indices(3,3));
  this->Rx.setElement(0.5*ONE_c,Indices(4,3));	this->Ry.setElement(-0.5*I_c,Indices(4,3));	this->Rz.setElement(0.5*ONE_c,Indices(4,4));
  
  //Get Jsquare as Casimir Operator
  this->Jsquared=this->Lx*this->Lx+this->Ly*this->Ly+this->Lz*this->Lz;
  //this->Jsquared=this->Rx*this->Rx+this->Ry*this->Ry+this->Rz*this->Rz;
  
  //Gauge operators
  this->U00=mwArray(Indices(5,5));	this->U00.fillWithZero();
  this->U01=mwArray(Indices(5,5));	this->U01.fillWithZero();
  this->U10=mwArray(Indices(5,5));	this->U10.fillWithZero();
  this->U11=mwArray(Indices(5,5));	this->U11.fillWithZero();
  
  this->U00.setElement(1.0/sqrt(2.0)*ONE_c,Indices(2,0));	this->U00.setElement(1.0/sqrt(2.0)*ONE_c,Indices(4,2));
  this->U01.setElement(-1.0/sqrt(2.0)*ONE_c,Indices(2,1));	this->U01.setElement(1.0/sqrt(2.0)*ONE_c,Indices(3,2));
  this->U10.setElement(1.0/sqrt(2.0)*ONE_c,Indices(1,2));	this->U10.setElement(-1.0/sqrt(2.0)*ONE_c,Indices(2,3));
  this->U11.setElement(1.0/sqrt(2.0)*ONE_c,Indices(0,2));	this->U11.setElement(1.0/sqrt(2.0)*ONE_c,Indices(2,4));
  
  //Adjoints of the Gauge operators
  this->U00ad=this->U00;	this->U00ad.transpose(true);
  this->U01ad=this->U01;	this->U01ad.transpose(true);
  this->U10ad=this->U10;	this->U10ad.transpose(true);
  this->U11ad=this->U11;	this->U11ad.transpose(true);
  
  //Identity for a link
  this->id_link=identityMatrix(this->d_link);
  
  //I need often identites to construct a single body operator, therefore I prepare those and save them
  this->id_link_mpo=this->id_link;
  this->id_fermi_mpo=this->id_fermi;
  this->id_link_mpo.reshape(Indices(this->d_link,1,this->d_link,1));
  this->id_fermi_mpo.reshape(Indices(this->d_fermi,1,this->d_fermi,1));    
  
  
  //For debugging
  /*cout << "Lx=" << this->Lx << endl;
  cout << "Ly=" << this->Ly << endl;
  cout << "Lz=" << this->Lz << endl;
  
  cout << "Rx=" << this->Rx << endl;
  cout << "Ry=" << this->Ry << endl;
  cout << "Rz=" << this->Rz << endl;
  
  cout << "U00=" << this->U00 << endl;
  cout << "U01=" << this->U01 << endl;
  cout << "U10=" << this->U10 << endl;
  cout << "U11=" << this->U11 << endl;
  
  cout << "U00ad=" << this->U00ad << endl;
  cout << "U01ad=" << this->U01ad << endl;
  cout << "U10ad=" << this->U10ad << endl;
  cout << "U11ad=" << this->U11ad << endl;
  
  cout << "J²=" << this->Jsquared << endl;
  cout << "J² aus rechten Generatoren=" << Rx*this->Rx+this->Ry*this->Ry+this->Rz*this->Rz << endl;*/
  
}

//Build the Hamiltonian MPO
void SU2Hamiltonian::initHamiltonian()
{
  cout << "Building Hamiltonian MPO with" << endl
       << "epsilon: " << this->epsilon << endl
       << "mu:      " << this->mu << endl
       << "g:       " << this->g << endl;
  
  //Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(10,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //Jsquared
      Gauge_operators.setElement(this->Jsquared.getElement(Indices(i,j)),Indices(1,i,j));
      //U00
      Gauge_operators.setElement(this->U00.getElement(Indices(i,j)),Indices(2,i,j));
      //U01
      Gauge_operators.setElement(this->U01.getElement(Indices(i,j)),Indices(3,i,j));
      //U10
      Gauge_operators.setElement(this->U10.getElement(Indices(i,j)),Indices(4,i,j));
      //U11
      Gauge_operators.setElement(this->U11.getElement(Indices(i,j)),Indices(5,i,j));
      //U00 dagger
      Gauge_operators.setElement(this->U00ad.getElement(Indices(i,j)),Indices(6,i,j));
      //U01 dagger
      Gauge_operators.setElement(this->U01ad.getElement(Indices(i,j)),Indices(7,i,j));
      //U10 dagger
      Gauge_operators.setElement(this->U10ad.getElement(Indices(i,j)),Indices(8,i,j));
      //U11 dagger
      Gauge_operators.setElement(this->U11ad.getElement(Indices(i,j)),Indices(9,i,j));       
     }
  }   
  //cout << "Gauge_operators" << Gauge_operators <<endl;
  Gauge_operators.reshape(Indices(10,this->d_link*this->d_link));  
  //cout << "Bosonic operators done" << endl;
  
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
  //cout << "Fermionic operators done" << endl;
  
  //Now prepare the matrices
  int D_bond=10;
  mwArray FermiFirst,Fermi1Odd,Fermi2Odd,Fermi1Even,Fermi2Even,Fermi1Last,Fermi2Last,Gauge;
  
  FermiFirst=mwArray(Indices(1,D_bond,4));	FermiFirst.fillWithZero();
  Fermi1Odd=mwArray(Indices(D_bond,D_bond,4));	Fermi1Odd.fillWithZero();
  Fermi2Odd=mwArray(Indices(D_bond,D_bond,4));	Fermi2Odd.fillWithZero();
  Fermi1Even=mwArray(Indices(D_bond,D_bond,4));	Fermi1Even.fillWithZero();
  Fermi2Even=mwArray(Indices(D_bond,D_bond,4));	Fermi2Even.fillWithZero();
  Fermi1Last=mwArray(Indices(D_bond,D_bond,4));	Fermi1Last.fillWithZero();
  Fermi2Last=mwArray(Indices(D_bond,1,4));	Fermi2Last.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,10));	Gauge.fillWithZero();
  
  //Determine the sign of the last fermionic site(s)
  double sign = (double) pow(-1,this->N_fermions);
  //cout << "Sign=" << sign << endl;
  
  //Identity
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(-0.5*this->mu*ONE_c,Indices(0,9,0));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi1Odd.setElement(-0.5*this->mu*ONE_c,Indices(0,9,0));
  Fermi1Odd.setElement(ONE_c,Indices(9,9,0));
  
  Fermi1Even.setElement(ONE_c,Indices(0,0,0));
  Fermi1Even.setElement(0.5*this->mu*ONE_c,Indices(0,9,0));
  Fermi1Even.setElement(ONE_c,Indices(9,9,0));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi2Odd.setElement(-0.5*this->mu*ONE_c,Indices(0,9,0));
  Fermi2Odd.setElement(ONE_c,Indices(9,9,0));
  
  Fermi2Even.setElement(ONE_c,Indices(0,0,0));
  Fermi2Even.setElement(0.5*this->mu*ONE_c,Indices(0,9,0));
  Fermi2Even.setElement(ONE_c,Indices(9,9,0));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(9,9,0));
  
  Fermi1Last.setElement(ONE_c,Indices(0,0,0));
  Fermi1Last.setElement(sign*0.5*this->mu*ONE_c,Indices(0,9,0));
  Fermi1Last.setElement(ONE_c,Indices(9,9,0));
  
  Fermi2Last.setElement(sign*0.5*this->mu*ONE_c,Indices(0,0,0));
  Fermi2Last.setElement(ONE_c,Indices(9,0,0));
  
  //sigma_plus
  FermiFirst.setElement(ONE_c,Indices(0,1,1));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,1,1));
  Fermi1Odd.setElement(ONE_c,Indices(7,9,1));
  
  Fermi1Even.setElement(ONE_c,Indices(0,1,1));
  Fermi1Even.setElement(ONE_c,Indices(7,9,1));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,3,1));
  Fermi2Odd.setElement(ONE_c,Indices(8,9,1));
  
  Fermi2Even.setElement(ONE_c,Indices(0,3,1));
  Fermi2Even.setElement(ONE_c,Indices(8,9,1));
  
  Fermi1Last.setElement(ONE_c,Indices(7,9,1));
  
  Fermi2Last.setElement(ONE_c,Indices(8,0,1));
  
  //sigma_minus
  FermiFirst.setElement(ONE_c,Indices(0,2,2));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,2,2));
  Fermi1Odd.setElement(ONE_c,Indices(5,9,2));
  
  Fermi1Even.setElement(ONE_c,Indices(0,2,2));
  Fermi1Even.setElement(ONE_c,Indices(5,9,2));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,4,2));
  Fermi2Odd.setElement(ONE_c,Indices(6,9,2));
  
  Fermi2Even.setElement(ONE_c,Indices(0,4,2));
  Fermi2Even.setElement(ONE_c,Indices(6,9,2));
  
  Fermi1Last.setElement(ONE_c,Indices(5,9,2));
  
  Fermi2Last.setElement(ONE_c,Indices(6,0,2));
  
  //sigma_z
  FermiFirst.setElement(-0.5*this->mu*ONE_c,Indices(0,9,3));
  
  Fermi1Odd.setElement(-0.5*this->mu*ONE_c,Indices(0,9,3));
  Fermi1Odd.setElement(ONE_c,Indices(6,6,3));
  Fermi1Odd.setElement(ONE_c,Indices(8,8,3));
  
  Fermi1Even.setElement(0.5*this->mu*ONE_c,Indices(0,9,3));
  Fermi1Even.setElement(ONE_c,Indices(6,6,3));
  Fermi1Even.setElement(ONE_c,Indices(8,8,3));
  
  Fermi2Odd.setElement(-0.5*this->mu*ONE_c,Indices(0,9,3));
  Fermi2Odd.setElement(ONE_c,Indices(1,1,3));
  Fermi2Odd.setElement(ONE_c,Indices(2,2,3));
  
  Fermi2Even.setElement(0.5*this->mu*ONE_c,Indices(0,9,3));
  Fermi2Even.setElement(ONE_c,Indices(1,1,3));
  Fermi2Even.setElement(ONE_c,Indices(2,2,3));
  
  Fermi1Last.setElement(sign*0.5*this->mu*ONE_c,Indices(0,9,3));
  Fermi1Last.setElement(ONE_c,Indices(6,6,3));
  Fermi1Last.setElement(ONE_c,Indices(8,8,3));
  
  Fermi2Last.setElement(sign*0.5*this->mu*ONE_c,Indices(0,0,3));
  
  //Jsquared
  Gauge.setElement(0.5*this->g*this->g*ONE_c,Indices(0,9,1));
  
  //U00
  Gauge.setElement(this->epsilon*ONE_c,Indices(1,5,2));
  
  //U01
  Gauge.setElement(this->epsilon*I_c,Indices(1,6,3));
  
  //U10
  Gauge.setElement(-this->epsilon*I_c,Indices(3,5,4));
  
  //U11
  Gauge.setElement(this->epsilon*ONE_c,Indices(3,6,5));
  
  //U00ad
  Gauge.setElement(this->epsilon*ONE_c,Indices(2,7,6));
  
  //U01ad
  Gauge.setElement(-this->epsilon*I_c,Indices(2,8,7));
  
  //U10ad
  Gauge.setElement(this->epsilon*I_c,Indices(4,7,8));
  
  //U11ad
  Gauge.setElement(this->epsilon*ONE_c,Indices(4,8,9));
  
  //Reshape the Matrices
  FermiFirst.reshape(Indices(1*D_bond,4));
  Fermi1Odd.reshape(Indices(D_bond*D_bond,4));
  Fermi1Even.reshape(Indices(D_bond*D_bond,4));
  Fermi2Odd.reshape(Indices(D_bond*D_bond,4));
  Fermi2Even.reshape(Indices(D_bond*D_bond,4));
  Fermi1Last.reshape(Indices(D_bond*D_bond,4));
  Fermi2Last.reshape(Indices(D_bond*1,4));
  Gauge.reshape(Indices(D_bond*D_bond,10));
  
  /*cout << "FermiFirst " << FermiFirst 
       << "Fermi1Odd  " << Fermi1Odd 
       << "Fermi1Even " << Fermi1Even 
       << "Fermi2Odd  " << Fermi2Odd 
       << "Fermi2Even " << Fermi2Even 
       << "Fermi1Last " << Fermi1Last 
       << "Fermi2Last " << Fermi2Last 
       << "Gauge      " << Gauge << endl;*/
  
  //Start filling the matrices in the MPO
  mwArray res;
  bool saved,odd_saved,even_saved;
  Operator* Local_op; 
  int site_index;
  
  //Fill the MPO with operators for fermionic species 1
  odd_saved=false;
  even_saved=false;
  site_index=0;
  for(int k=0; k<this->N_total-1; k+=3)
  {
    //Get the site index, counter is safer than direct computation
    site_index++;
    //cout << "Site index (constant mass): " << site_index << endl;
    if(k==0)
    {
      //First site (always fermionic)
      res=reshape(FermiFirst*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true); 
    }
    else if(k==(this->N_total-2))
    {
      //Last site of fermionic species 1
      res=reshape(Fermi1Last*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true);       
    }
    else
    {
      //Rest of sites
      if((site_index%2)==0 && !even_saved)
      {
	res=reshape(Fermi1Even*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true);       
	even_saved=true;
      }
      else if((site_index%2)==1 && !odd_saved)
      {
	res=reshape(Fermi1Odd*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true);       
	odd_saved=true;
      }
      else
      {
	hamil.setOp(k,&hamil.getOp(k-6),false);
      }
    }    
  }
  
  //Fill the MPO with operators for fermionic species 2
  odd_saved=false;
  even_saved=false;
  site_index=0;
  for(int k=1; k<this->N_total; k+=3)
  {
    //Get the site index, counter is safer than direct computation
    site_index++;
    //cout << "Site index (constant mass): " << site_index << endl;
    
    if(k==(this->N_total-1))
    {
      //Last site of fermionic species 1
      res=reshape(Fermi2Last*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true);       
    }
    else
    {
      //Rest of sites
      if((site_index%2)==0 && !even_saved)
      {
	res=reshape(Fermi2Even*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true);       
	even_saved=true;
      }
      else if((site_index%2)==1 && !odd_saved)
      {
	res=reshape(Fermi2Odd*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true);       
	odd_saved=true;
      }
      else
      {
	hamil.setOp(k,&hamil.getOp(k-6),false);
      }
    }    
  }
  
  //Fill the MPO with operators for the links
  saved=false;
  for(int k=2; k<this->N_total-2; k+=3)
  {
    //cout << "k=" << k << endl;
    if(saved)
    {
      hamil.setOp(k,&hamil.getOp(k-3),false);
    }
    else
    {
      //Save the gauge operator once
      res=reshape(Gauge*Gauge_operators,Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true);       
      saved=true;
    }    
  } 
  
  cout << "Setting Matrix Elements done" << endl;
  
  //cout << "Hamiltonian " << hamil << endl;
  //char name[100];
  //sprintf (name, "SU2Hamiltonian_N%deps%fM%fg%f.m",this->N_fermions,this->epsilon,this->mu,this->g);
  //hamil.exportForMatlab("SU2Hamiltonian.m");
  
}

void SU2Hamiltonian::initHamiltonian_sitedependendmass()
{
  //TODO somehow bad style because of duplicated code. Could be optimized such that I always deal with this structure
  
  
  cout << "Building Hamiltonian MPO with" << endl
       << "epsilon: " << this->epsilon << endl
       << "mu:      " << this->mu_site << endl
       << "g:       " << this->g << endl;
  
  //Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(10,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //Jsquared
      Gauge_operators.setElement(this->Jsquared.getElement(Indices(i,j)),Indices(1,i,j));
      //U00
      Gauge_operators.setElement(this->U00.getElement(Indices(i,j)),Indices(2,i,j));
      //U01
      Gauge_operators.setElement(this->U01.getElement(Indices(i,j)),Indices(3,i,j));
      //U10
      Gauge_operators.setElement(this->U10.getElement(Indices(i,j)),Indices(4,i,j));
      //U11
      Gauge_operators.setElement(this->U11.getElement(Indices(i,j)),Indices(5,i,j));
      //U00 dagger
      Gauge_operators.setElement(this->U00ad.getElement(Indices(i,j)),Indices(6,i,j));
      //U01 dagger
      Gauge_operators.setElement(this->U01ad.getElement(Indices(i,j)),Indices(7,i,j));
      //U10 dagger
      Gauge_operators.setElement(this->U10ad.getElement(Indices(i,j)),Indices(8,i,j));
      //U11 dagger
      Gauge_operators.setElement(this->U11ad.getElement(Indices(i,j)),Indices(9,i,j));       
     }
  }   
  //cout << "Gauge_operators" << Gauge_operators <<endl;
  Gauge_operators.reshape(Indices(10,this->d_link*this->d_link));  
  //cout << "Bosonic operators done" << endl;
  
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
  //cout << "Fermionic operators done" << endl;
  
  //Now prepare the matrices
  int D_bond=10;
  mwArray FermiFirst,Fermi1Odd,Fermi2Odd,Fermi1Even,Fermi2Even,Fermi1Last,Fermi2Last,Gauge;
  
  FermiFirst=mwArray(Indices(1,D_bond,4));	FermiFirst.fillWithZero();
  Fermi1Odd=mwArray(Indices(D_bond,D_bond,4));	Fermi1Odd.fillWithZero();
  Fermi2Odd=mwArray(Indices(D_bond,D_bond,4));	Fermi2Odd.fillWithZero();
  Fermi1Even=mwArray(Indices(D_bond,D_bond,4));	Fermi1Even.fillWithZero();
  Fermi2Even=mwArray(Indices(D_bond,D_bond,4));	Fermi2Even.fillWithZero();
  Fermi1Last=mwArray(Indices(D_bond,D_bond,4));	Fermi1Last.fillWithZero();
  Fermi2Last=mwArray(Indices(D_bond,1,4));	Fermi2Last.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,10));	Gauge.fillWithZero();
  
  //Determine the sign of the last fermionic site(s)
  double sign = (double) pow(-1,this->N_fermions);
  //cout << "Sign=" << sign << endl;
  
  //Identity
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(-0.5*this->mu*ONE_c,Indices(0,9,0));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi1Odd.setElement(-0.5*this->mu*ONE_c,Indices(0,9,0));
  Fermi1Odd.setElement(ONE_c,Indices(9,9,0));
  
  Fermi1Even.setElement(ONE_c,Indices(0,0,0));
  Fermi1Even.setElement(0.5*this->mu*ONE_c,Indices(0,9,0));
  Fermi1Even.setElement(ONE_c,Indices(9,9,0));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi2Odd.setElement(-0.5*this->mu*ONE_c,Indices(0,9,0));
  Fermi2Odd.setElement(ONE_c,Indices(9,9,0));
  
  Fermi2Even.setElement(ONE_c,Indices(0,0,0));
  Fermi2Even.setElement(0.5*this->mu*ONE_c,Indices(0,9,0));
  Fermi2Even.setElement(ONE_c,Indices(9,9,0));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(9,9,0));
  
  Fermi1Last.setElement(ONE_c,Indices(0,0,0));
  Fermi1Last.setElement(sign*0.5*this->mu*ONE_c,Indices(0,9,0));
  Fermi1Last.setElement(ONE_c,Indices(9,9,0));
  
  Fermi2Last.setElement(sign*0.5*this->mu*ONE_c,Indices(0,0,0));
  Fermi2Last.setElement(ONE_c,Indices(9,0,0));
  
  //sigma_plus
  FermiFirst.setElement(ONE_c,Indices(0,1,1));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,1,1));
  Fermi1Odd.setElement(ONE_c,Indices(7,9,1));
  
  Fermi1Even.setElement(ONE_c,Indices(0,1,1));
  Fermi1Even.setElement(ONE_c,Indices(7,9,1));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,3,1));
  Fermi2Odd.setElement(ONE_c,Indices(8,9,1));
  
  Fermi2Even.setElement(ONE_c,Indices(0,3,1));
  Fermi2Even.setElement(ONE_c,Indices(8,9,1));
  
  Fermi1Last.setElement(ONE_c,Indices(7,9,1));
  
  Fermi2Last.setElement(ONE_c,Indices(8,0,1));
  
  //sigma_minus
  FermiFirst.setElement(ONE_c,Indices(0,2,2));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,2,2));
  Fermi1Odd.setElement(ONE_c,Indices(5,9,2));
  
  Fermi1Even.setElement(ONE_c,Indices(0,2,2));
  Fermi1Even.setElement(ONE_c,Indices(5,9,2));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,4,2));
  Fermi2Odd.setElement(ONE_c,Indices(6,9,2));
  
  Fermi2Even.setElement(ONE_c,Indices(0,4,2));
  Fermi2Even.setElement(ONE_c,Indices(6,9,2));
  
  Fermi1Last.setElement(ONE_c,Indices(5,9,2));
  
  Fermi2Last.setElement(ONE_c,Indices(6,0,2));
  
  //sigma_z
  FermiFirst.setElement(-0.5*this->mu*ONE_c,Indices(0,9,3));
  
  Fermi1Odd.setElement(-0.5*this->mu*ONE_c,Indices(0,9,3));
  Fermi1Odd.setElement(ONE_c,Indices(6,6,3));
  Fermi1Odd.setElement(ONE_c,Indices(8,8,3));
  
  Fermi1Even.setElement(0.5*this->mu*ONE_c,Indices(0,9,3));
  Fermi1Even.setElement(ONE_c,Indices(6,6,3));
  Fermi1Even.setElement(ONE_c,Indices(8,8,3));
  
  Fermi2Odd.setElement(-0.5*this->mu*ONE_c,Indices(0,9,3));
  Fermi2Odd.setElement(ONE_c,Indices(1,1,3));
  Fermi2Odd.setElement(ONE_c,Indices(2,2,3));
  
  Fermi2Even.setElement(0.5*this->mu*ONE_c,Indices(0,9,3));
  Fermi2Even.setElement(ONE_c,Indices(1,1,3));
  Fermi2Even.setElement(ONE_c,Indices(2,2,3));
  
  Fermi1Last.setElement(sign*0.5*this->mu*ONE_c,Indices(0,9,3));
  Fermi1Last.setElement(ONE_c,Indices(6,6,3));
  Fermi1Last.setElement(ONE_c,Indices(8,8,3));
  
  Fermi2Last.setElement(sign*0.5*this->mu*ONE_c,Indices(0,0,3));
  
  //Jsquared
  Gauge.setElement(0.5*this->g*this->g*ONE_c,Indices(0,9,1));
  
  //U00
  Gauge.setElement(this->epsilon*ONE_c,Indices(1,5,2));
  
  //U01
  Gauge.setElement(this->epsilon*I_c,Indices(1,6,3));
  
  //U10
  Gauge.setElement(-this->epsilon*I_c,Indices(3,5,4));
  
  //U11
  Gauge.setElement(this->epsilon*ONE_c,Indices(3,6,5));
  
  //U00ad
  Gauge.setElement(this->epsilon*ONE_c,Indices(2,7,6));
  
  //U01ad
  Gauge.setElement(-this->epsilon*I_c,Indices(2,8,7));
  
  //U10ad
  Gauge.setElement(this->epsilon*I_c,Indices(4,7,8));
  
  //U11ad
  Gauge.setElement(this->epsilon*ONE_c,Indices(4,8,9));
  
  //Reshape the gauge Matrices
  Gauge.reshape(Indices(D_bond*D_bond,10));
  
  /*cout << "FermiFirst " << FermiFirst 
       << "Fermi1Odd  " << Fermi1Odd 
       << "Fermi1Even " << Fermi1Even 
       << "Fermi2Odd  " << Fermi2Odd 
       << "Fermi2Even " << Fermi2Even 
       << "Fermi1Last " << Fermi1Last 
       << "Fermi2Last " << Fermi2Last 
       << "Gauge      " << Gauge << endl;*/
  
  //Start filling the matrices in the MPO
  mwArray res,fermi_matrix;
  bool saved,odd_saved,even_saved;
  Operator* Local_op; 
  int site_index;
  
  //Fill the MPO with operators for fermionic species 1
  site_index=0;
  for(int k=0; k<this->N_total-1; k+=3)
  {
    //Get the site index (counter is safer than direct computation)
    site_index++;    
    //cout << "Site index (site dependend mass): " << site_index << endl;
    
    //(Re)Set the mass terms
    //Identity
    FermiFirst.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,9,0));    
    Fermi1Odd.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,9,0));    
    Fermi1Even.setElement(0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,9,0));    
    Fermi1Last.setElement(sign*0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,9,0));
    
    //sigma_z
    FermiFirst.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,9,3));    
    Fermi1Odd.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,9,3));    
    Fermi1Even.setElement(0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,9,3));    
    Fermi1Last.setElement(sign*0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,9,3));
    if(k==0)
    {
      //First site (always fermionic)
      fermi_matrix=FermiFirst;
      fermi_matrix.reshape(Indices(1*D_bond,4));
      res=reshape(fermi_matrix*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true); 
    }
    else if(k==(this->N_total-2))
    {
      //Last site of fermionic species 1
      fermi_matrix=Fermi1Last;
      fermi_matrix.reshape(Indices(D_bond*D_bond,4));
      res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true);       
    }
    else
    {
      //Rest of sites
      if((site_index%2)==0)
      {
	fermi_matrix=Fermi1Even;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true);       
	even_saved=true;
      }
      else
      {
	fermi_matrix=Fermi1Odd;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true);       
	odd_saved=true;
      }
    }    
  }
  
  //Fill the MPO with operators for fermionic species 2
  site_index=0;
  for(int k=1; k<this->N_total; k+=3)
  {
    //Get the site index (counter is safer than direct computation)
    site_index++;
    //cout << "Site index (site dependend mass): " << site_index << endl;
    
    //Identity    
    Fermi2Odd.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,9,0));
    Fermi2Even.setElement(0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,9,0));    
    Fermi2Last.setElement(sign*0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,0,0));
    
    //sigma_z
    Fermi2Odd.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,9,3));
    Fermi2Even.setElement(0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,9,3));    
    Fermi2Last.setElement(sign*0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,0,3));
    
    if(k==(this->N_total-1))
    {
      //Last site of fermionic species 2
      fermi_matrix=Fermi2Last;
      fermi_matrix.reshape(Indices(D_bond*1,4));
      res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true);       
    }
    else
    {
      //Rest of sites
      if((site_index%2)==0)
      {
	fermi_matrix=Fermi2Even;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true);       
      }
      else
      {
	fermi_matrix=Fermi2Odd;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true);       
      }
    }    
  }
  
  //Fill the MPO with operators for the links
  saved=false;
  for(int k=2; k<this->N_total-2; k+=3)
  {
    //cout << "k=" << k << endl;
    if(saved)
    {
      hamil.setOp(k,&hamil.getOp(k-3),false);
    }
    else
    {
      //Save the gauge operator once
      res=reshape(Gauge*Gauge_operators,Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true);       
      saved=true;
    }    
  } 
  
  cout << "Setting Matrix Elements done" << endl;
  
  //cout << "Hamiltonian " << hamil << endl;
  //char name[100];
  //sprintf (name, "SU2Hamiltonian_N%deps%fM%fg%f.m",this->N_fermions,this->epsilon,this->mu,this->g);
  //hamil.exportForMatlab("SU2HamiltonianlocalSiteDependendMass.m");
  
}

void SU2Hamiltonian::initHamiltonian_sitedependendmass_charge(vector<int> charge_pos, double lambda)
{
  //TODO somehow bad style because of duplicated code. Could be optimized such that I always deal with this structure
  
  
  cout << "Building Hamiltonian MPO with" << endl
       << "epsilon: " << this->epsilon << endl
       << "mu:      " << this->mu_site << endl
       << "g:       " << this->g << endl
       << "charges: " << charge_pos << endl
       << "lambda:  " << lambda << endl;
  
  //Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(10,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //Jsquared
      Gauge_operators.setElement(this->Jsquared.getElement(Indices(i,j)),Indices(1,i,j));
      //U00
      Gauge_operators.setElement(this->U00.getElement(Indices(i,j)),Indices(2,i,j));
      //U01
      Gauge_operators.setElement(this->U01.getElement(Indices(i,j)),Indices(3,i,j));
      //U10
      Gauge_operators.setElement(this->U10.getElement(Indices(i,j)),Indices(4,i,j));
      //U11
      Gauge_operators.setElement(this->U11.getElement(Indices(i,j)),Indices(5,i,j));
      //U00 dagger
      Gauge_operators.setElement(this->U00ad.getElement(Indices(i,j)),Indices(6,i,j));
      //U01 dagger
      Gauge_operators.setElement(this->U01ad.getElement(Indices(i,j)),Indices(7,i,j));
      //U10 dagger
      Gauge_operators.setElement(this->U10ad.getElement(Indices(i,j)),Indices(8,i,j));
      //U11 dagger
      Gauge_operators.setElement(this->U11ad.getElement(Indices(i,j)),Indices(9,i,j));       
     }
  }   
  //cout << "Gauge_operators" << Gauge_operators <<endl;
  Gauge_operators.reshape(Indices(10,this->d_link*this->d_link));  
  //cout << "Bosonic operators done" << endl;
  
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
  //cout << "Fermionic operators done" << endl;
  
  //Now prepare the matrices
  int D_bond=11;
  mwArray FermiFirst,Fermi1Odd,Fermi2Odd,Fermi1Even,Fermi2Even,Fermi1Last,Fermi2Last,Gauge,FermiFirstCharge,Fermi1OddCharge,Fermi2OddCharge,Fermi1EvenCharge,Fermi2EvenCharge,Fermi1LastCharge,Fermi2LastCharge;
  
  FermiFirst=mwArray(Indices(1,D_bond,4));	FermiFirst.fillWithZero();
  Fermi1Odd=mwArray(Indices(D_bond,D_bond,4));	Fermi1Odd.fillWithZero();
  Fermi2Odd=mwArray(Indices(D_bond,D_bond,4));	Fermi2Odd.fillWithZero();
  Fermi1Even=mwArray(Indices(D_bond,D_bond,4));	Fermi1Even.fillWithZero();
  Fermi2Even=mwArray(Indices(D_bond,D_bond,4));	Fermi2Even.fillWithZero();
  Fermi1Last=mwArray(Indices(D_bond,D_bond,4));	Fermi1Last.fillWithZero();
  Fermi2Last=mwArray(Indices(D_bond,1,4));	Fermi2Last.fillWithZero();
  
  FermiFirstCharge=mwArray(Indices(1,D_bond,4));	FermiFirstCharge.fillWithZero();
  Fermi1OddCharge=mwArray(Indices(D_bond,D_bond,4));	Fermi1OddCharge.fillWithZero();
  Fermi2OddCharge=mwArray(Indices(D_bond,D_bond,4));	Fermi2OddCharge.fillWithZero();
  Fermi1EvenCharge=mwArray(Indices(D_bond,D_bond,4));	Fermi1EvenCharge.fillWithZero();
  Fermi2EvenCharge=mwArray(Indices(D_bond,D_bond,4));	Fermi2EvenCharge.fillWithZero();
  Fermi1LastCharge=mwArray(Indices(D_bond,D_bond,4));	Fermi1LastCharge.fillWithZero();
  Fermi2LastCharge=mwArray(Indices(D_bond,1,4));	Fermi2LastCharge.fillWithZero();
  
  Gauge=mwArray(Indices(D_bond,D_bond,10));	Gauge.fillWithZero();
  
  //Determine the sign of the last fermionic site(s)
  double sign = (double) pow(-1,this->N_fermions);
  //cout << "Sign=" << sign << endl;
  
  //Identity
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(-0.5*this->mu*ONE_c,Indices(0,10,0));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi1Odd.setElement(-0.5*this->mu*ONE_c,Indices(0,10,0));
  Fermi1Odd.setElement(ONE_c,Indices(10,10,0));
  
  Fermi1Even.setElement(ONE_c,Indices(0,0,0));
  Fermi1Even.setElement(0.5*this->mu*ONE_c,Indices(0,10,0));
  Fermi1Even.setElement(ONE_c,Indices(10,10,0));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi2Odd.setElement(-0.5*this->mu*ONE_c,Indices(0,10,0));
  Fermi2Odd.setElement(ONE_c,Indices(10,10,0));
  
  Fermi2Even.setElement(ONE_c,Indices(0,0,0));
  Fermi2Even.setElement(0.5*this->mu*ONE_c,Indices(0,10,0));
  Fermi2Even.setElement(ONE_c,Indices(10,10,0));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(10,10,0));
  
  Fermi1Last.setElement(ONE_c,Indices(0,0,0));
  Fermi1Last.setElement(sign*0.5*this->mu*ONE_c,Indices(0,10,0));
  Fermi1Last.setElement(ONE_c,Indices(10,10,0));
  
  Fermi2Last.setElement(sign*0.5*this->mu*ONE_c,Indices(0,0,0));
  Fermi2Last.setElement(ONE_c,Indices(10,0,0));
  
  //sigma_plus
  FermiFirst.setElement(ONE_c,Indices(0,1,1));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,1,1));
  Fermi1Odd.setElement(ONE_c,Indices(7,10,1));
  
  Fermi1Even.setElement(ONE_c,Indices(0,1,1));
  Fermi1Even.setElement(ONE_c,Indices(7,10,1));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,3,1));
  Fermi2Odd.setElement(ONE_c,Indices(8,10,1));
  
  Fermi2Even.setElement(ONE_c,Indices(0,3,1));
  Fermi2Even.setElement(ONE_c,Indices(8,10,1));
  
  Fermi1Last.setElement(ONE_c,Indices(7,10,1));
  
  Fermi2Last.setElement(ONE_c,Indices(8,0,1));
  
  //sigma_minus
  FermiFirst.setElement(ONE_c,Indices(0,2,2));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,2,2));
  Fermi1Odd.setElement(ONE_c,Indices(5,10,2));
  
  Fermi1Even.setElement(ONE_c,Indices(0,2,2));
  Fermi1Even.setElement(ONE_c,Indices(5,10,2));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,4,2));
  Fermi2Odd.setElement(ONE_c,Indices(6,10,2));
  
  Fermi2Even.setElement(ONE_c,Indices(0,4,2));
  Fermi2Even.setElement(ONE_c,Indices(6,10,2));
  
  Fermi1Last.setElement(ONE_c,Indices(5,10,2));
  
  Fermi2Last.setElement(ONE_c,Indices(6,0,2));
  
  //sigma_z
  FermiFirst.setElement(-0.5*this->mu*ONE_c,Indices(0,10,3));
  
  Fermi1Odd.setElement(-0.5*this->mu*ONE_c,Indices(0,10,3));
  Fermi1Odd.setElement(ONE_c,Indices(6,6,3));
  Fermi1Odd.setElement(ONE_c,Indices(8,8,3));
  
  Fermi1Even.setElement(0.5*this->mu*ONE_c,Indices(0,10,3));
  Fermi1Even.setElement(ONE_c,Indices(6,6,3));
  Fermi1Even.setElement(ONE_c,Indices(8,8,3));
  
  Fermi2Odd.setElement(-0.5*this->mu*ONE_c,Indices(0,10,3));
  Fermi2Odd.setElement(ONE_c,Indices(1,1,3));
  Fermi2Odd.setElement(ONE_c,Indices(2,2,3));
  
  Fermi2Even.setElement(0.5*this->mu*ONE_c,Indices(0,10,3));
  Fermi2Even.setElement(ONE_c,Indices(1,1,3));
  Fermi2Even.setElement(ONE_c,Indices(2,2,3));
  
  Fermi1Last.setElement(sign*0.5*this->mu*ONE_c,Indices(0,10,3));
  Fermi1Last.setElement(ONE_c,Indices(6,6,3));
  Fermi1Last.setElement(ONE_c,Indices(8,8,3));
  
  Fermi2Last.setElement(sign*0.5*this->mu*ONE_c,Indices(0,0,3));
  
  //Copy the same thing to the charge sites and set the additional entries
  FermiFirstCharge = FermiFirst;
  Fermi1OddCharge = Fermi1Odd;
  Fermi1EvenCharge = Fermi1Even;
  Fermi2OddCharge = Fermi2Odd;
  Fermi2EvenCharge = Fermi2Even;
  Fermi1LastCharge = Fermi1Last;
  Fermi2LastCharge = Fermi2Last;  
  
  //Identities
  FermiFirstCharge.setElement((-0.5*this->mu+lambda/2.0)*ONE_c,Indices(0,10,0));  
  Fermi1OddCharge.setElement((-0.5*this->mu+lambda/2.0)*ONE_c,Indices(0,10,0));  
  Fermi1EvenCharge.setElement((0.5*this->mu+lambda/2.0)*ONE_c,Indices(0,10,0));  
  Fermi2OddCharge.setElement((-0.5*this->mu+lambda/2.0)*ONE_c,Indices(0,10,0));  
  Fermi2EvenCharge.setElement((0.5*this->mu+lambda/2.0)*ONE_c,Indices(0,10,0));  
  Fermi1LastCharge.setElement((sign*0.5*this->mu+lambda/2.0)*ONE_c,Indices(0,10,0));  
  Fermi2LastCharge.setElement((sign*0.5*this->mu+lambda/2.0)*ONE_c,Indices(0,0,0));
  
  FermiFirstCharge.setElement(lambda*ONE_c,Indices(0,9,3));  
  Fermi1OddCharge.setElement(lambda*ONE_c,Indices(0,9,3));  
  Fermi2OddCharge.setElement(ONE_c,Indices(9,10,3));  
  Fermi1EvenCharge.setElement(lambda*ONE_c,Indices(0,9,3));  
  Fermi2EvenCharge.setElement(ONE_c,Indices(9,10,3));
  Fermi1LastCharge.setElement(lambda*ONE_c,Indices(0,9,3));
  Fermi2LastCharge.setElement(ONE_c,Indices(9,0,3));
  
  //Jsquared
  Gauge.setElement(0.5*this->g*this->g*ONE_c,Indices(0,10,1));  
  //U00
  Gauge.setElement(this->epsilon*ONE_c,Indices(1,5,2));  
  //U01
  Gauge.setElement(this->epsilon*I_c,Indices(1,6,3));  
  //U10
  Gauge.setElement(-this->epsilon*I_c,Indices(3,5,4));  
  //U11
  Gauge.setElement(this->epsilon*ONE_c,Indices(3,6,5));  
  //U00ad
  Gauge.setElement(this->epsilon*ONE_c,Indices(2,7,6));  
  //U01ad
  Gauge.setElement(-this->epsilon*I_c,Indices(2,8,7));  
  //U10ad
  Gauge.setElement(this->epsilon*I_c,Indices(4,7,8));  
  //U11ad
  Gauge.setElement(this->epsilon*ONE_c,Indices(4,8,9));
  
  //Reshape the gauge Matrices
  Gauge.reshape(Indices(D_bond*D_bond,10));
  
  /*cout << "FermiFirst " << FermiFirst 
       << "Fermi1Odd  " << Fermi1Odd 
       << "Fermi1Even " << Fermi1Even 
       << "Fermi2Odd  " << Fermi2Odd 
       << "Fermi2Even " << Fermi2Even 
       << "Fermi1Last " << Fermi1Last 
       << "Fermi2Last " << Fermi2Last 
       << "Gauge      " << Gauge << endl;*/
  
  //Start filling the matrices in the MPO
  mwArray res,fermi_matrix,First,Last,Odd,Even;
  bool saved,odd_saved,even_saved;
  Operator* Local_op; 
  int site_index;
  
  //Fill the MPO with operators for fermionic species 1
  site_index=0;
  for(int k=0; k<this->N_total-1; k+=3)
  {
    //Get the site index (counter is safer than direct computation)
    site_index++;    
    //cout << "Site index (site depended mass): " << site_index << endl;
    
    //(Re)Set the mass terms
    if(std::find(charge_pos.begin(), charge_pos.end(), site_index)!=charge_pos.end())
    {
      //Case that I encountered a charge site
      cout << "Charge site at " << site_index << endl;
      //Identity
      FermiFirstCharge.setElement((-0.5*this->mu_site[site_index-1]+lambda/2.0)*ONE_c,Indices(0,10,0));    
      Fermi1OddCharge.setElement((-0.5*this->mu_site[site_index-1]+lambda/2.0)*ONE_c,Indices(0,10,0));    
      Fermi1EvenCharge.setElement((0.5*this->mu_site[site_index-1]+lambda/2.0)*ONE_c,Indices(0,10,0));    
      Fermi1LastCharge.setElement((sign*0.5*this->mu_site[site_index-1]+lambda/2.0)*ONE_c,Indices(0,10,0));
      
      //sigma_z
      FermiFirstCharge.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,3));    
      Fermi1OddCharge.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,3));    
      Fermi1EvenCharge.setElement(0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,3));    
      Fermi1LastCharge.setElement(sign*0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,3));
      
      First=FermiFirstCharge;
      Odd=Fermi1OddCharge;
      Even=Fermi1EvenCharge;
      Last=Fermi1LastCharge;
    }
    else
    {    
      //Identity
      FermiFirst.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,0));    
      Fermi1Odd.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,0));    
      Fermi1Even.setElement(0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,0));    
      Fermi1Last.setElement(sign*0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,0));
      
      //sigma_z
      FermiFirst.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,3));    
      Fermi1Odd.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,3));    
      Fermi1Even.setElement(0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,3));    
      Fermi1Last.setElement(sign*0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,3));
      
      First=FermiFirst;
      Odd=Fermi1Odd;
      Even=Fermi1Even;
      Last = Fermi1Last;
    }
    if(k==0)
    {
      //First site (always fermionic)
      fermi_matrix=First;
      fermi_matrix.reshape(Indices(1*D_bond,4));
      res=reshape(fermi_matrix*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true); 
    }
    else if(k==(this->N_total-2))
    {
      //Last site of fermionic species 1
      fermi_matrix=Last;
      fermi_matrix.reshape(Indices(D_bond*D_bond,4));
      res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true);       
    }
    else
    {
      //Rest of sites
      if((site_index%2)==0)
      {
	fermi_matrix=Even;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true);       
      }
      else
      {
	fermi_matrix=Odd;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true);       
      }
    }    
  }
  
  //Fill the MPO with operators for fermionic species 2
  site_index=0;
  for(int k=1; k<this->N_total; k+=3)
  {
    //Get the site index (counter is safer than direct computation)
    site_index++;
    //cout << "Site index (site depended mass): " << site_index << endl;
    
    if(std::find(charge_pos.begin(), charge_pos.end(), site_index)!=charge_pos.end())
    {
      //Case that I encounter a charge site
      cout << "Charge site at " << site_index << endl;
      //Identity    
      Fermi2OddCharge.setElement((-0.5*this->mu_site[site_index-1]+lambda/2.0)*ONE_c,Indices(0,10,0));
      Fermi2EvenCharge.setElement((0.5*this->mu_site[site_index-1]+lambda/2.0)*ONE_c,Indices(0,10,0));    
      Fermi2LastCharge.setElement((sign*0.5*this->mu_site[site_index-1]+lambda/2.0)*ONE_c,Indices(0,0,0));
      
      //sigma_z
      Fermi2OddCharge.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,3));
      Fermi2EvenCharge.setElement(0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,3));    
      Fermi2LastCharge.setElement(sign*0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,0,3));
      
      Odd=Fermi2OddCharge;
      Even=Fermi2EvenCharge;
      Last=Fermi2LastCharge;
    }
    else
    {
      //Identity    
      Fermi2Odd.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,0));
      Fermi2Even.setElement(0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,0));    
      Fermi2Last.setElement(sign*0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,0,0));
      
      //sigma_z
      Fermi2Odd.setElement(-0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,3));
      Fermi2Even.setElement(0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,10,3));    
      Fermi2Last.setElement(sign*0.5*this->mu_site[site_index-1]*ONE_c,Indices(0,0,3));
      
      Odd=Fermi2Odd;
      Even=Fermi2Even;
      Last=Fermi2Last;
    }
    
    if(k==(this->N_total-1))
    {
      //Last site of fermionic species 2
      fermi_matrix=Last;
      fermi_matrix.reshape(Indices(D_bond*1,4));
      res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true);       
    }
    else
    {
      //Rest of sites
      if((site_index%2)==0)
      {
	fermi_matrix=Even;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true);       
      }
      else
      {
	fermi_matrix=Odd;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	hamil.setOp(k,Local_op,true);       
      }
    }    
  }
  
  //Fill the MPO with operators for the links
  saved=false;
  for(int k=2; k<this->N_total-2; k+=3)
  {
    //cout << "k=" << k << endl;
    if(saved)
    {
      hamil.setOp(k,&hamil.getOp(k-3),false);
    }
    else
    {
      //Save the gauge operator once
      res=reshape(Gauge*Gauge_operators,Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      hamil.setOp(k,Local_op,true);       
      saved=true;
    }    
  } 
  
  cout << "Setting Matrix Elements done" << endl;
  
  //cout << "Hamiltonian " << hamil << endl;
  //char name[100];
  //sprintf (name, "SU2Hamiltonian_N%deps%fM%fg%f.m",this->N_fermions,this->epsilon,this->mu,this->g);
  //hamil.exportForMatlab("SU2HamiltonianlocalSiteDependendMassCharge.m");
  
}

//Construct the time evolution MPO
void SU2Hamiltonian::constructEvolutionMPO(MPO &Umpo, double dt, bool imag_time) const
{
  //Factor from the taylor expansion of the time evolution operator
  complex_t taylor_factor=-1.0*dt*I_c;
  //cout << endl << "Taylor factor " << taylor_factor << endl;
  
  if(imag_time)
  {
    cout << "Constructing imaginary time evolution MPO" << endl;
    taylor_factor = -I_c*taylor_factor;
  }
  
  //Prepare the MPO
  Umpo.initLength(this->N_total);
  
  //Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(10,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //Jsquared
      Gauge_operators.setElement(this->Jsquared.getElement(Indices(i,j)),Indices(1,i,j));
      //U00
      Gauge_operators.setElement(this->U00.getElement(Indices(i,j)),Indices(2,i,j));
      //U01
      Gauge_operators.setElement(this->U01.getElement(Indices(i,j)),Indices(3,i,j));
      //U10
      Gauge_operators.setElement(this->U10.getElement(Indices(i,j)),Indices(4,i,j));
      //U11
      Gauge_operators.setElement(this->U11.getElement(Indices(i,j)),Indices(5,i,j));
      //U00 dagger
      Gauge_operators.setElement(this->U00ad.getElement(Indices(i,j)),Indices(6,i,j));
      //U01 dagger
      Gauge_operators.setElement(this->U01ad.getElement(Indices(i,j)),Indices(7,i,j));
      //U10 dagger
      Gauge_operators.setElement(this->U10ad.getElement(Indices(i,j)),Indices(8,i,j));
      //U11 dagger
      Gauge_operators.setElement(this->U11ad.getElement(Indices(i,j)),Indices(9,i,j));       
     }
  }   
  //cout << "Gauge_operators" << Gauge_operators <<endl;
  Gauge_operators.reshape(Indices(10,this->d_link*this->d_link));  
  
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
  
  //Now prepare the matrices
  int D_bond=10;
  mwArray FermiFirst,Fermi1Odd,Fermi2Odd,Fermi1Even,Fermi2Even,Fermi1Last,Fermi2Last,Gauge;
  
  FermiFirst=mwArray(Indices(1,D_bond,4));	FermiFirst.fillWithZero();
  Fermi1Odd=mwArray(Indices(D_bond,D_bond,4));	Fermi1Odd.fillWithZero();
  Fermi2Odd=mwArray(Indices(D_bond,D_bond,4));	Fermi2Odd.fillWithZero();
  Fermi1Even=mwArray(Indices(D_bond,D_bond,4));	Fermi1Even.fillWithZero();
  Fermi2Even=mwArray(Indices(D_bond,D_bond,4));	Fermi2Even.fillWithZero();
  Fermi1Last=mwArray(Indices(D_bond,D_bond,4));	Fermi1Last.fillWithZero();
  Fermi2Last=mwArray(Indices(D_bond,1,4));	Fermi2Last.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,10));	Gauge.fillWithZero();
  
  //Determine the sign of the last fermionic site(s)
  double sign = (double) pow(-1,this->N_fermions);
  
  //Identity
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu*taylor_factor,Indices(0,9,0));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi1Odd.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu*taylor_factor,Indices(0,9,0));
  Fermi1Odd.setElement(ONE_c,Indices(9,9,0));
  
  Fermi1Even.setElement(ONE_c,Indices(0,0,0));
  Fermi1Even.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu*taylor_factor,Indices(0,9,0));
  Fermi1Even.setElement(ONE_c,Indices(9,9,0));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi2Odd.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu*taylor_factor,Indices(0,9,0));
  Fermi2Odd.setElement(ONE_c,Indices(9,9,0));
  
  Fermi2Even.setElement(ONE_c,Indices(0,0,0));
  Fermi2Even.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu*taylor_factor,Indices(0,9,0));
  Fermi2Even.setElement(ONE_c,Indices(9,9,0));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(9,9,0));
  
  Fermi1Last.setElement(ONE_c,Indices(0,0,0));
  Fermi1Last.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu*taylor_factor,Indices(0,9,0));
  Fermi1Last.setElement(ONE_c,Indices(9,9,0));
  
  Fermi2Last.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu*taylor_factor,Indices(0,0,0));
  Fermi2Last.setElement(ONE_c,Indices(9,0,0));
  
  //sigma_plus
  FermiFirst.setElement(ONE_c,Indices(0,1,1));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,1,1));
  Fermi1Odd.setElement(ONE_c,Indices(7,9,1));
  
  Fermi1Even.setElement(ONE_c,Indices(0,1,1));
  Fermi1Even.setElement(ONE_c,Indices(7,9,1));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,3,1));
  Fermi2Odd.setElement(ONE_c,Indices(8,9,1));
  
  Fermi2Even.setElement(ONE_c,Indices(0,3,1));
  Fermi2Even.setElement(ONE_c,Indices(8,9,1));
  
  Fermi1Last.setElement(ONE_c,Indices(7,9,1));
  
  Fermi2Last.setElement(ONE_c,Indices(8,0,1));
  
  //sigma_minus
  FermiFirst.setElement(ONE_c,Indices(0,2,2));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,2,2));
  Fermi1Odd.setElement(ONE_c,Indices(5,9,2));
  
  Fermi1Even.setElement(ONE_c,Indices(0,2,2));
  Fermi1Even.setElement(ONE_c,Indices(5,9,2));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,4,2));
  Fermi2Odd.setElement(ONE_c,Indices(6,9,2));
  
  Fermi2Even.setElement(ONE_c,Indices(0,4,2));
  Fermi2Even.setElement(ONE_c,Indices(6,9,2));
  
  Fermi1Last.setElement(ONE_c,Indices(5,9,2));
  
  Fermi2Last.setElement(ONE_c,Indices(6,0,2));
  
  //sigma_z
  FermiFirst.setElement(-0.5*this->mu*taylor_factor,Indices(0,9,3));
  
  Fermi1Odd.setElement(-0.5*this->mu*taylor_factor,Indices(0,9,3));
  Fermi1Odd.setElement(ONE_c,Indices(6,6,3));
  Fermi1Odd.setElement(ONE_c,Indices(8,8,3));
  
  Fermi1Even.setElement(0.5*this->mu*taylor_factor,Indices(0,9,3));
  Fermi1Even.setElement(ONE_c,Indices(6,6,3));
  Fermi1Even.setElement(ONE_c,Indices(8,8,3));
  
  Fermi2Odd.setElement(-0.5*this->mu*taylor_factor,Indices(0,9,3));
  Fermi2Odd.setElement(ONE_c,Indices(1,1,3));
  Fermi2Odd.setElement(ONE_c,Indices(2,2,3));
  
  Fermi2Even.setElement(0.5*this->mu*taylor_factor,Indices(0,9,3));
  Fermi2Even.setElement(ONE_c,Indices(1,1,3));
  Fermi2Even.setElement(ONE_c,Indices(2,2,3));
  
  Fermi1Last.setElement(sign*0.5*this->mu*taylor_factor,Indices(0,9,3));
  Fermi1Last.setElement(ONE_c,Indices(6,6,3));
  Fermi1Last.setElement(ONE_c,Indices(8,8,3));
  
  Fermi2Last.setElement(sign*0.5*this->mu*taylor_factor,Indices(0,0,3));
  
  //Jsquared
  Gauge.setElement(0.5*this->g*this->g*taylor_factor,Indices(0,9,1));
  
  //U00
  Gauge.setElement(this->epsilon*taylor_factor,Indices(1,5,2));
  
  //U01
  Gauge.setElement(this->epsilon*I_c*taylor_factor,Indices(1,6,3));
  
  //U10
  Gauge.setElement(-this->epsilon*I_c*taylor_factor,Indices(3,5,4));
  
  //U11
  Gauge.setElement(this->epsilon*taylor_factor,Indices(3,6,5));
  
  //U00ad
  Gauge.setElement(this->epsilon*taylor_factor,Indices(2,7,6));
  
  //U01ad
  Gauge.setElement(-this->epsilon*I_c*taylor_factor,Indices(2,8,7));
  
  //U10ad
  Gauge.setElement(this->epsilon*I_c*taylor_factor,Indices(4,7,8));
  
  //U11ad
  Gauge.setElement(this->epsilon*taylor_factor,Indices(4,8,9));
  
  //Reshape the Matrices
  FermiFirst.reshape(Indices(1*D_bond,4));
  Fermi1Odd.reshape(Indices(D_bond*D_bond,4));
  Fermi1Even.reshape(Indices(D_bond*D_bond,4));
  Fermi2Odd.reshape(Indices(D_bond*D_bond,4));
  Fermi2Even.reshape(Indices(D_bond*D_bond,4));
  Fermi1Last.reshape(Indices(D_bond*D_bond,4));
  Fermi2Last.reshape(Indices(D_bond*1,4));
  Gauge.reshape(Indices(D_bond*D_bond,10));
  
  /*cout << "FermiFirst " << FermiFirst 
       << "Fermi1Odd  " << Fermi1Odd 
       << "Fermi1Even " << Fermi1Even 
       << "Fermi2Odd  " << Fermi2Odd 
       << "Fermi2Even " << Fermi2Even 
       << "Fermi1Last " << Fermi1Last 
       << "Fermi2Last " << Fermi2Last 
       << "Gauge      " << Gauge << endl;*/
  
  //Start filling the matrices in the MPO
  mwArray res;
  bool saved,odd_saved,even_saved;
  Operator* Local_op; 
  int site_index;
  
  //Fill the MPO with operators for fermionic species 1
  odd_saved=false;
  even_saved=false;
  site_index=0;
  for(int k=0; k<this->N_total-1; k+=3)
  {
    //Get the site index, counter is safer than direct computation
    site_index++;
    //cout << "Site index (constant mass): " << site_index << endl;
    if(k==0)
    {
      //First site (always fermionic)
      res=reshape(FermiFirst*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      Umpo.setOp(k,Local_op,true); 
    }
    else if(k==(this->N_total-2))
    {
      //Last site of fermionic species 1
      res=reshape(Fermi1Last*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      Umpo.setOp(k,Local_op,true);       
    }
    else
    {
      //Rest of sites
      if((site_index%2)==0 && !even_saved)
      {
	res=reshape(Fermi1Even*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	Umpo.setOp(k,Local_op,true);       
	even_saved=true;
      }
      else if((site_index%2)==1 && !odd_saved)
      {
	res=reshape(Fermi1Odd*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	Umpo.setOp(k,Local_op,true);       
	odd_saved=true;
      }
      else
      {
	Umpo.setOp(k,&Umpo.getOp(k-6),false);
      }
    }    
  }
  
  //Fill the MPO with operators for fermionic species 2
  odd_saved=false;
  even_saved=false;
  site_index=0;
  for(int k=1; k<this->N_total; k+=3)
  {
    //Get the site index, counter is safer than direct computation
    site_index++;
    //cout << "Site index (constant mass): " << site_index << endl;
    
    if(k==(this->N_total-1))
    {
      //Last site of fermionic species 1
      res=reshape(Fermi2Last*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      Umpo.setOp(k,Local_op,true);       
    }
    else
    {
      //Rest of sites
      if((site_index%2)==0 && !even_saved)
      {
	res=reshape(Fermi2Even*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	Umpo.setOp(k,Local_op,true);       
	even_saved=true;
      }
      else if((site_index%2)==1 && !odd_saved)
      {
	res=reshape(Fermi2Odd*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	Umpo.setOp(k,Local_op,true);       
	odd_saved=true;
      }
      else
      {
	Umpo.setOp(k,&Umpo.getOp(k-6),false);
      }
    }    
  }
  
  //Fill the MPO with operators for the links
  saved=false;
  for(int k=2; k<this->N_total-2; k+=3)
  {
    //cout << "k=" << k << endl;
    if(saved)
    {
      Umpo.setOp(k,&Umpo.getOp(k-3),false);
    }
    else
    {
      //Save the gauge operator once
      res=reshape(Gauge*Gauge_operators,Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      Umpo.setOp(k,Local_op,true);       
      saved=true;
    }    
  } 
  //Umpo.exportForMatlab("SU2Umpo.m");
}

//Construct the time evolution MPO for the case of site dependend mass
void SU2Hamiltonian::constructEvolutionMPOSiteDepMass(MPO &Umpo, double dt, bool imag_time) const
{
  //Factor from the taylor expansion of the time evolution operator
  complex_t taylor_factor=-1.0*dt*I_c;
  //cout << endl << "Taylor factor " << taylor_factor << endl;
  
  if(imag_time)
  {
    cout << "Constructing imaginary time evolution MPO" << endl;
    taylor_factor = -I_c*taylor_factor;
  }
  
  //Prepare the MPO
  Umpo.initLength(this->N_total);
  
  //Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(10,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //Jsquared
      Gauge_operators.setElement(this->Jsquared.getElement(Indices(i,j)),Indices(1,i,j));
      //U00
      Gauge_operators.setElement(this->U00.getElement(Indices(i,j)),Indices(2,i,j));
      //U01
      Gauge_operators.setElement(this->U01.getElement(Indices(i,j)),Indices(3,i,j));
      //U10
      Gauge_operators.setElement(this->U10.getElement(Indices(i,j)),Indices(4,i,j));
      //U11
      Gauge_operators.setElement(this->U11.getElement(Indices(i,j)),Indices(5,i,j));
      //U00 dagger
      Gauge_operators.setElement(this->U00ad.getElement(Indices(i,j)),Indices(6,i,j));
      //U01 dagger
      Gauge_operators.setElement(this->U01ad.getElement(Indices(i,j)),Indices(7,i,j));
      //U10 dagger
      Gauge_operators.setElement(this->U10ad.getElement(Indices(i,j)),Indices(8,i,j));
      //U11 dagger
      Gauge_operators.setElement(this->U11ad.getElement(Indices(i,j)),Indices(9,i,j));       
     }
  }   
  //cout << "Gauge_operators" << Gauge_operators <<endl;
  Gauge_operators.reshape(Indices(10,this->d_link*this->d_link));  
  
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
  
  //Now prepare the matrices
  int D_bond=10;
  mwArray FermiFirst,Fermi1Odd,Fermi2Odd,Fermi1Even,Fermi2Even,Fermi1Last,Fermi2Last,Gauge;
  
  FermiFirst=mwArray(Indices(1,D_bond,4));	FermiFirst.fillWithZero();
  Fermi1Odd=mwArray(Indices(D_bond,D_bond,4));	Fermi1Odd.fillWithZero();
  Fermi2Odd=mwArray(Indices(D_bond,D_bond,4));	Fermi2Odd.fillWithZero();
  Fermi1Even=mwArray(Indices(D_bond,D_bond,4));	Fermi1Even.fillWithZero();
  Fermi2Even=mwArray(Indices(D_bond,D_bond,4));	Fermi2Even.fillWithZero();
  Fermi1Last=mwArray(Indices(D_bond,D_bond,4));	Fermi1Last.fillWithZero();
  Fermi2Last=mwArray(Indices(D_bond,1,4));	Fermi2Last.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,10));	Gauge.fillWithZero();
  
  //Determine the sign of the last fermionic site(s)
  double sign = (double) pow(-1,this->N_fermions);
  
  //Identity
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu*taylor_factor,Indices(0,9,0));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi1Odd.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu*taylor_factor,Indices(0,9,0));
  Fermi1Odd.setElement(ONE_c,Indices(9,9,0));
  
  Fermi1Even.setElement(ONE_c,Indices(0,0,0));
  Fermi1Even.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu*taylor_factor,Indices(0,9,0));
  Fermi1Even.setElement(ONE_c,Indices(9,9,0));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi2Odd.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu*taylor_factor,Indices(0,9,0));
  Fermi2Odd.setElement(ONE_c,Indices(9,9,0));
  
  Fermi2Even.setElement(ONE_c,Indices(0,0,0));
  Fermi2Even.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu*taylor_factor,Indices(0,9,0));
  Fermi2Even.setElement(ONE_c,Indices(9,9,0));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(9,9,0));
  
  Fermi1Last.setElement(ONE_c,Indices(0,0,0));
  Fermi1Last.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu*taylor_factor,Indices(0,9,0));
  Fermi1Last.setElement(ONE_c,Indices(9,9,0));
  
  Fermi2Last.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu*taylor_factor,Indices(0,0,0));
  Fermi2Last.setElement(ONE_c,Indices(9,0,0));
  
  //sigma_plus
  FermiFirst.setElement(ONE_c,Indices(0,1,1));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,1,1));
  Fermi1Odd.setElement(ONE_c,Indices(7,9,1));
  
  Fermi1Even.setElement(ONE_c,Indices(0,1,1));
  Fermi1Even.setElement(ONE_c,Indices(7,9,1));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,3,1));
  Fermi2Odd.setElement(ONE_c,Indices(8,9,1));
  
  Fermi2Even.setElement(ONE_c,Indices(0,3,1));
  Fermi2Even.setElement(ONE_c,Indices(8,9,1));
  
  Fermi1Last.setElement(ONE_c,Indices(7,9,1));
  
  Fermi2Last.setElement(ONE_c,Indices(8,0,1));
  
  //sigma_minus
  FermiFirst.setElement(ONE_c,Indices(0,2,2));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,2,2));
  Fermi1Odd.setElement(ONE_c,Indices(5,9,2));
  
  Fermi1Even.setElement(ONE_c,Indices(0,2,2));
  Fermi1Even.setElement(ONE_c,Indices(5,9,2));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,4,2));
  Fermi2Odd.setElement(ONE_c,Indices(6,9,2));
  
  Fermi2Even.setElement(ONE_c,Indices(0,4,2));
  Fermi2Even.setElement(ONE_c,Indices(6,9,2));
  
  Fermi1Last.setElement(ONE_c,Indices(5,9,2));
  
  Fermi2Last.setElement(ONE_c,Indices(6,0,2));
  
  //sigma_z
  FermiFirst.setElement(-0.5*this->mu*taylor_factor,Indices(0,9,3));
  
  Fermi1Odd.setElement(-0.5*this->mu*taylor_factor,Indices(0,9,3));
  Fermi1Odd.setElement(ONE_c,Indices(6,6,3));
  Fermi1Odd.setElement(ONE_c,Indices(8,8,3));
  
  Fermi1Even.setElement(0.5*this->mu*taylor_factor,Indices(0,9,3));
  Fermi1Even.setElement(ONE_c,Indices(6,6,3));
  Fermi1Even.setElement(ONE_c,Indices(8,8,3));
  
  Fermi2Odd.setElement(-0.5*this->mu*taylor_factor,Indices(0,9,3));
  Fermi2Odd.setElement(ONE_c,Indices(1,1,3));
  Fermi2Odd.setElement(ONE_c,Indices(2,2,3));
  
  Fermi2Even.setElement(0.5*this->mu*taylor_factor,Indices(0,9,3));
  Fermi2Even.setElement(ONE_c,Indices(1,1,3));
  Fermi2Even.setElement(ONE_c,Indices(2,2,3));
  
  Fermi1Last.setElement(sign*0.5*this->mu*taylor_factor,Indices(0,9,3));
  Fermi1Last.setElement(ONE_c,Indices(6,6,3));
  Fermi1Last.setElement(ONE_c,Indices(8,8,3));
  
  Fermi2Last.setElement(sign*0.5*this->mu*taylor_factor,Indices(0,0,3));
  
  //Jsquared
  Gauge.setElement(0.5*this->g*this->g*taylor_factor,Indices(0,9,1));
  
  //U00
  Gauge.setElement(this->epsilon*taylor_factor,Indices(1,5,2));
  
  //U01
  Gauge.setElement(this->epsilon*I_c*taylor_factor,Indices(1,6,3));
  
  //U10
  Gauge.setElement(-this->epsilon*I_c*taylor_factor,Indices(3,5,4));
  
  //U11
  Gauge.setElement(this->epsilon*taylor_factor,Indices(3,6,5));
  
  //U00ad
  Gauge.setElement(this->epsilon*taylor_factor,Indices(2,7,6));
  
  //U01ad
  Gauge.setElement(-this->epsilon*I_c*taylor_factor,Indices(2,8,7));
  
  //U10ad
  Gauge.setElement(this->epsilon*I_c*taylor_factor,Indices(4,7,8));
  
  //U11ad
  Gauge.setElement(this->epsilon*taylor_factor,Indices(4,8,9));
  
  //Reshape the gauge Matrices
  Gauge.reshape(Indices(D_bond*D_bond,10));
  
  /*cout << "FermiFirst " << FermiFirst 
       << "Fermi1Odd  " << Fermi1Odd 
       << "Fermi1Even " << Fermi1Even 
       << "Fermi2Odd  " << Fermi2Odd 
       << "Fermi2Even " << Fermi2Even 
       << "Fermi1Last " << Fermi1Last 
       << "Fermi2Last " << Fermi2Last 
       << "Gauge      " << Gauge << endl;*/
  
  //Start filling the matrices in the MPO
  mwArray res,fermi_matrix;
  bool saved,odd_saved,even_saved;
  Operator* Local_op; 
  int site_index;
  
  //Fill the MPO with operators for fermionic species 1
  site_index=0;
  for(int k=0; k<this->N_total-1; k+=3)
  {
    //Get the site index (counter is safer than direct computation)
    site_index ++;
    //cout << "Site index (site dependend mass): " << site_index << endl;
    
    //Identity
    FermiFirst.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,9,0));
    Fermi1Odd.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,9,0));
    Fermi1Even.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,9,0));
    Fermi1Last.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,9,0));
    
    //sigma_z
    FermiFirst.setElement(-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,9,3));    
    Fermi1Odd.setElement(-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,9,3));    
    Fermi1Even.setElement(0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,9,3));    
    Fermi1Last.setElement(sign*0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,9,3));
    
    if(k==0)
    {
      //First site (always fermionic)
      fermi_matrix=FermiFirst;
      fermi_matrix.reshape(Indices(1*D_bond,4));
      res=reshape(fermi_matrix*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      Umpo.setOp(k,Local_op,true); 
    }
    else if(k==(this->N_total-2))
    {
      //Last site of fermionic species 1
      fermi_matrix=Fermi1Last;
      fermi_matrix.reshape(Indices(D_bond*D_bond,4));
      res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      Umpo.setOp(k,Local_op,true);       
    }
    else
    {
      //Rest of sites
      if((site_index%2)==0)
      {
	fermi_matrix=Fermi1Even;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	Umpo.setOp(k,Local_op,true);       
      }
      else if((site_index%2)==1 && !odd_saved)
      {
	fermi_matrix=Fermi1Odd;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	Umpo.setOp(k,Local_op,true);       
      }
    }    
  }
  
  //Fill the MPO with operators for fermionic species 2
  site_index=0;
  for(int k=1; k<this->N_total; k+=3)
  {
    //Get the site index (counter is safer than direct computation)
    site_index++;
    //cout << "Site index (site dependend mass): " << site_index << endl;
    
    //Identity    
    Fermi2Odd.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,9,0));
    Fermi2Even.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,9,0));
    Fermi2Last.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,0,0));
    
    //sigma_z
    Fermi2Odd.setElement(-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,9,3));
    Fermi2Even.setElement(0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,9,3));
    Fermi1Last.setElement(sign*0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,9,3));
    Fermi2Last.setElement(sign*0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,0,3));//Identity
  
    if(k==(this->N_total-1))
    {
      //Last site of fermionic species 1
      fermi_matrix=Fermi2Last;
      fermi_matrix.reshape(Indices(D_bond*1,4));
      res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      Umpo.setOp(k,Local_op,true);       
    }
    else
    {
      //Rest of sites
      if((site_index%2)==0)
      {
	fermi_matrix=Fermi2Even;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	Umpo.setOp(k,Local_op,true);       
	even_saved=true;
      }
      else
      {
	fermi_matrix=Fermi2Odd;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	Umpo.setOp(k,Local_op,true);       
	odd_saved=true;
      }
    }    
  }
  
  //Fill the MPO with operators for the links
  saved=false;
  for(int k=2; k<this->N_total-2; k+=3)
  {
    //cout << "k=" << k << endl;
    if(saved)
    {
      Umpo.setOp(k,&Umpo.getOp(k-3),false);
    }
    else
    {
      //Save the gauge operator once
      res=reshape(Gauge*Gauge_operators,Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      Umpo.setOp(k,Local_op,true);       
      saved=true;
    }    
  } 
  //Umpo.exportForMatlab("SU2Umpo.m");
}

void SU2Hamiltonian::constructEvolutionMPOSiteDepMassCharge(MPO &Umpo, double dt, vector<int> charge_pos, double lambda, bool imag_time) const
{
  //TODO somehow bad style because of duplicated code. Could be optimized such that I always deal with this structure
  
  //Factor from the taylor expansion of the time evolution operator
  complex_t taylor_factor=-1.0*dt*I_c;
  //cout << endl << "Taylor factor " << taylor_factor << endl;
  
  if(imag_time)
  {
    cout << "Constructing imaginary time evolution MPO" << endl;
    taylor_factor = -I_c*taylor_factor;
  }
  
  //Prepare the MPO
  Umpo.initLength(this->N_total);  
  
  //Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(10,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //Jsquared
      Gauge_operators.setElement(this->Jsquared.getElement(Indices(i,j)),Indices(1,i,j));
      //U00
      Gauge_operators.setElement(this->U00.getElement(Indices(i,j)),Indices(2,i,j));
      //U01
      Gauge_operators.setElement(this->U01.getElement(Indices(i,j)),Indices(3,i,j));
      //U10
      Gauge_operators.setElement(this->U10.getElement(Indices(i,j)),Indices(4,i,j));
      //U11
      Gauge_operators.setElement(this->U11.getElement(Indices(i,j)),Indices(5,i,j));
      //U00 dagger
      Gauge_operators.setElement(this->U00ad.getElement(Indices(i,j)),Indices(6,i,j));
      //U01 dagger
      Gauge_operators.setElement(this->U01ad.getElement(Indices(i,j)),Indices(7,i,j));
      //U10 dagger
      Gauge_operators.setElement(this->U10ad.getElement(Indices(i,j)),Indices(8,i,j));
      //U11 dagger
      Gauge_operators.setElement(this->U11ad.getElement(Indices(i,j)),Indices(9,i,j));       
     }
  }   
  //cout << "Gauge_operators" << Gauge_operators <<endl;
  Gauge_operators.reshape(Indices(10,this->d_link*this->d_link));  
  //cout << "Bosonic operators done" << endl;
  
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
  //cout << "Fermionic operators done" << endl;
  
  //Now prepare the matrices
  int D_bond=11;
  mwArray FermiFirst,Fermi1Odd,Fermi2Odd,Fermi1Even,Fermi2Even,Fermi1Last,Fermi2Last,Gauge,FermiFirstCharge,Fermi1OddCharge,Fermi2OddCharge,Fermi1EvenCharge,Fermi2EvenCharge,Fermi1LastCharge,Fermi2LastCharge;
  
  FermiFirst=mwArray(Indices(1,D_bond,4));	FermiFirst.fillWithZero();
  Fermi1Odd=mwArray(Indices(D_bond,D_bond,4));	Fermi1Odd.fillWithZero();
  Fermi2Odd=mwArray(Indices(D_bond,D_bond,4));	Fermi2Odd.fillWithZero();
  Fermi1Even=mwArray(Indices(D_bond,D_bond,4));	Fermi1Even.fillWithZero();
  Fermi2Even=mwArray(Indices(D_bond,D_bond,4));	Fermi2Even.fillWithZero();
  Fermi1Last=mwArray(Indices(D_bond,D_bond,4));	Fermi1Last.fillWithZero();
  Fermi2Last=mwArray(Indices(D_bond,1,4));	Fermi2Last.fillWithZero();
  
  FermiFirstCharge=mwArray(Indices(1,D_bond,4));	FermiFirstCharge.fillWithZero();
  Fermi1OddCharge=mwArray(Indices(D_bond,D_bond,4));	Fermi1OddCharge.fillWithZero();
  Fermi2OddCharge=mwArray(Indices(D_bond,D_bond,4));	Fermi2OddCharge.fillWithZero();
  Fermi1EvenCharge=mwArray(Indices(D_bond,D_bond,4));	Fermi1EvenCharge.fillWithZero();
  Fermi2EvenCharge=mwArray(Indices(D_bond,D_bond,4));	Fermi2EvenCharge.fillWithZero();
  Fermi1LastCharge=mwArray(Indices(D_bond,D_bond,4));	Fermi1LastCharge.fillWithZero();
  Fermi2LastCharge=mwArray(Indices(D_bond,1,4));	Fermi2LastCharge.fillWithZero();
  
  Gauge=mwArray(Indices(D_bond,D_bond,10));	Gauge.fillWithZero();
  
  //Determine the sign of the last fermionic site(s)
  double sign = (double) pow(-1,this->N_fermions);
  //cout << "Sign=" << sign << endl;
  
  //Identity
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu*ONE_c*taylor_factor,Indices(0,10,0));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi1Odd.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu*ONE_c*taylor_factor,Indices(0,10,0));
  Fermi1Odd.setElement(ONE_c,Indices(10,10,0));
  
  Fermi1Even.setElement(ONE_c,Indices(0,0,0));
  Fermi1Even.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu*ONE_c*taylor_factor,Indices(0,10,0));
  Fermi1Even.setElement(ONE_c,Indices(10,10,0));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi2Odd.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu*ONE_c*taylor_factor,Indices(0,10,0));
  Fermi2Odd.setElement(ONE_c,Indices(10,10,0));
  
  Fermi2Even.setElement(ONE_c,Indices(0,0,0));
  Fermi2Even.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu*ONE_c*taylor_factor,Indices(0,10,0));
  Fermi2Even.setElement(ONE_c,Indices(10,10,0));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(10,10,0));
  
  Fermi1Last.setElement(ONE_c,Indices(0,0,0));
  Fermi1Last.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu*ONE_c*taylor_factor,Indices(0,10,0));
  Fermi1Last.setElement(ONE_c,Indices(10,10,0));
  
  Fermi2Last.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu*ONE_c*taylor_factor,Indices(0,0,0));
  Fermi2Last.setElement(ONE_c,Indices(10,0,0));
  
  //sigma_plus
  FermiFirst.setElement(ONE_c,Indices(0,1,1));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,1,1));
  Fermi1Odd.setElement(ONE_c,Indices(7,10,1));
  
  Fermi1Even.setElement(ONE_c,Indices(0,1,1));
  Fermi1Even.setElement(ONE_c,Indices(7,10,1));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,3,1));
  Fermi2Odd.setElement(ONE_c,Indices(8,10,1));
  
  Fermi2Even.setElement(ONE_c,Indices(0,3,1));
  Fermi2Even.setElement(ONE_c,Indices(8,10,1));
  
  Fermi1Last.setElement(ONE_c,Indices(7,10,1));
  
  Fermi2Last.setElement(ONE_c,Indices(8,0,1));
  
  //sigma_minus
  FermiFirst.setElement(ONE_c,Indices(0,2,2));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,2,2));
  Fermi1Odd.setElement(ONE_c,Indices(5,10,2));
  
  Fermi1Even.setElement(ONE_c,Indices(0,2,2));
  Fermi1Even.setElement(ONE_c,Indices(5,10,2));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,4,2));
  Fermi2Odd.setElement(ONE_c,Indices(6,10,2));
  
  Fermi2Even.setElement(ONE_c,Indices(0,4,2));
  Fermi2Even.setElement(ONE_c,Indices(6,10,2));
  
  Fermi1Last.setElement(ONE_c,Indices(5,10,2));
  
  Fermi2Last.setElement(ONE_c,Indices(6,0,2));
  
  //sigma_z
  FermiFirst.setElement(-0.5*this->mu*taylor_factor,Indices(0,10,3));
  
  Fermi1Odd.setElement(-0.5*this->mu*taylor_factor,Indices(0,10,3));
  Fermi1Odd.setElement(ONE_c,Indices(6,6,3));
  Fermi1Odd.setElement(ONE_c,Indices(8,8,3));
  
  Fermi1Even.setElement(0.5*this->mu*taylor_factor,Indices(0,10,3));
  Fermi1Even.setElement(ONE_c,Indices(6,6,3));
  Fermi1Even.setElement(ONE_c,Indices(8,8,3));
  
  Fermi2Odd.setElement(-0.5*this->mu*taylor_factor,Indices(0,10,3));
  Fermi2Odd.setElement(ONE_c,Indices(1,1,3));
  Fermi2Odd.setElement(ONE_c,Indices(2,2,3));
  
  Fermi2Even.setElement(0.5*this->mu*taylor_factor,Indices(0,10,3));
  Fermi2Even.setElement(ONE_c,Indices(1,1,3));
  Fermi2Even.setElement(ONE_c,Indices(2,2,3));
  
  Fermi1Last.setElement(sign*0.5*this->mu*taylor_factor,Indices(0,10,3));
  Fermi1Last.setElement(ONE_c,Indices(6,6,3));
  Fermi1Last.setElement(ONE_c,Indices(8,8,3));
  
  Fermi2Last.setElement(sign*0.5*this->mu*taylor_factor,Indices(0,0,3));
  
  //Copy the same thing to the charge sites and set the additional entries
  FermiFirstCharge = FermiFirst;
  Fermi1OddCharge = Fermi1Odd;
  Fermi1EvenCharge = Fermi1Even;
  Fermi2OddCharge = Fermi2Odd;
  Fermi2EvenCharge = Fermi2Even;
  Fermi1LastCharge = Fermi1Last;
  Fermi2LastCharge = Fermi2Last;  
  
  //Identities
  FermiFirstCharge.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu*taylor_factor+lambda/2.0*taylor_factor,Indices(0,10,0));  
  Fermi1OddCharge.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu*taylor_factor+lambda/2.0*taylor_factor,Indices(0,10,0));  
  Fermi1EvenCharge.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu*taylor_factor+lambda/2.0*taylor_factor,Indices(0,10,0));  
  Fermi2OddCharge.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu*taylor_factor+lambda/2.0*taylor_factor,Indices(0,10,0));  
  Fermi2EvenCharge.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu*taylor_factor+lambda/2.0*taylor_factor,Indices(0,10,0));  
  Fermi1LastCharge.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu*taylor_factor+lambda/2.0*taylor_factor,Indices(0,10,0));  
  Fermi2LastCharge.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu*taylor_factor+lambda/2.0*taylor_factor,Indices(0,0,0));
  
  //sigma_z
  FermiFirstCharge.setElement(lambda*taylor_factor,Indices(0,9,3));  
  Fermi1OddCharge.setElement(lambda*taylor_factor,Indices(0,9,3));  
  Fermi2OddCharge.setElement(ONE_c,Indices(9,10,3));  
  Fermi1EvenCharge.setElement(lambda*taylor_factor,Indices(0,9,3));  
  Fermi2EvenCharge.setElement(ONE_c,Indices(9,10,3));
  Fermi1LastCharge.setElement(lambda*taylor_factor,Indices(0,9,3));
  Fermi2LastCharge.setElement(ONE_c,Indices(9,0,3));
  
  //Jsquared
  Gauge.setElement(0.5*this->g*this->g*taylor_factor,Indices(0,10,1));  
  //U00
  Gauge.setElement(this->epsilon*taylor_factor,Indices(1,5,2));  
  //U01
  Gauge.setElement(this->epsilon*I_c*taylor_factor,Indices(1,6,3));  
  //U10
  Gauge.setElement(-this->epsilon*I_c*taylor_factor,Indices(3,5,4));  
  //U11
  Gauge.setElement(this->epsilon*taylor_factor,Indices(3,6,5));  
  //U00ad
  Gauge.setElement(this->epsilon*taylor_factor,Indices(2,7,6));  
  //U01ad
  Gauge.setElement(-this->epsilon*I_c*taylor_factor,Indices(2,8,7));  
  //U10ad
  Gauge.setElement(this->epsilon*I_c*taylor_factor,Indices(4,7,8));  
  //U11ad
  Gauge.setElement(this->epsilon*taylor_factor,Indices(4,8,9));
  
  //Reshape the gauge Matrices
  Gauge.reshape(Indices(D_bond*D_bond,10));
  
  /*cout << "FermiFirst " << FermiFirst 
       << "Fermi1Odd  " << Fermi1Odd 
       << "Fermi1Even " << Fermi1Even 
       << "Fermi2Odd  " << Fermi2Odd 
       << "Fermi2Even " << Fermi2Even 
       << "Fermi1Last " << Fermi1Last 
       << "Fermi2Last " << Fermi2Last 
       << "Gauge      " << Gauge << endl;*/
  
  //Start filling the matrices in the MPO
  mwArray res,fermi_matrix,First,Last,Odd,Even;
  bool saved;
  Operator* Local_op; 
  int site_index;
  
  //Fill the MPO with operators for fermionic species 1
  site_index=0;
  for(int k=0; k<this->N_total-1; k+=3)
  {
    //Get the site index (counter is safer than direct computation)
    site_index++;    
    //cout << "Site index (site depended mass): " << site_index << endl;
    
    //(Re)Set the mass terms
    if(std::find(charge_pos.begin(), charge_pos.end(), site_index)!=charge_pos.end())
    {
      //Case that I encountered a charge site
      cout << "Charge site at " << site_index << endl;
      //Identity
      FermiFirstCharge.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu_site[site_index-1]*taylor_factor+lambda/2.0*taylor_factor,Indices(0,10,0));    
      Fermi1OddCharge.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu_site[site_index-1]*taylor_factor+lambda/2.0*taylor_factor,Indices(0,10,0));    
      Fermi1EvenCharge.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu_site[site_index-1]*taylor_factor+lambda/2.0*taylor_factor,Indices(0,10,0));    
      Fermi1LastCharge.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu_site[site_index-1]*taylor_factor+lambda/2.0*taylor_factor,Indices(0,10,0));
      
      //sigma_z
      FermiFirstCharge.setElement(-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,3));    
      Fermi1OddCharge.setElement(-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,3));    
      Fermi1EvenCharge.setElement(0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,3));    
      Fermi1LastCharge.setElement(sign*0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,3));
      
      First=FermiFirstCharge;
      Odd=Fermi1OddCharge;
      Even=Fermi1EvenCharge;
      Last=Fermi1LastCharge;
    }
    else
    {    
      //Identity
      FermiFirst.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,0));    
      Fermi1Odd.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,0));    
      Fermi1Even.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,0));    
      Fermi1Last.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,0));
      
      //sigma_z
      FermiFirst.setElement(-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,3));    
      Fermi1Odd.setElement(-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,3));    
      Fermi1Even.setElement(0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,3));    
      Fermi1Last.setElement(sign*0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,3));
      
      First=FermiFirst;
      Odd=Fermi1Odd;
      Even=Fermi1Even;
      Last = Fermi1Last;
    }
    if(k==0)
    {
      //First site (always fermionic)
      fermi_matrix=First;
      fermi_matrix.reshape(Indices(1*D_bond,4));
      res=reshape(fermi_matrix*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      Umpo.setOp(k,Local_op,true); 
    }
    else if(k==(this->N_total-2))
    {
      //Last site of fermionic species 1
      fermi_matrix=Last;
      fermi_matrix.reshape(Indices(D_bond*D_bond,4));
      res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      Umpo.setOp(k,Local_op,true);       
    }
    else
    {
      //Rest of sites
      if((site_index%2)==0)
      {
	fermi_matrix=Even;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	Umpo.setOp(k,Local_op,true);       
      }
      else
      {
	fermi_matrix=Odd;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	Umpo.setOp(k,Local_op,true);       
      }
    }    
  }
  
  //Fill the MPO with operators for fermionic species 2
  site_index=0;
  for(int k=1; k<this->N_total; k+=3)
  {
    //Get the site index (counter is safer than direct computation)
    site_index++;
    //cout << "Site index (site depended mass): " << site_index << endl;
    
    if(std::find(charge_pos.begin(), charge_pos.end(), site_index)!=charge_pos.end())
    {
      //Case that I encounter a charge site
      cout << "Charge site at " << site_index << endl;
      //Identity    
      Fermi2OddCharge.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu_site[site_index-1]*taylor_factor+lambda/2.0*taylor_factor,Indices(0,10,0));
      Fermi2EvenCharge.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu_site[site_index-1]*taylor_factor+lambda/2.0*taylor_factor,Indices(0,10,0));    
      Fermi2LastCharge.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu_site[site_index-1]*taylor_factor+lambda/2.0*taylor_factor,Indices(0,0,0));
      
      //sigma_z
      Fermi2OddCharge.setElement(-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,3));
      Fermi2EvenCharge.setElement(0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,3));    
      Fermi2LastCharge.setElement(sign*0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,0,3));
      
      Odd=Fermi2OddCharge;
      Even=Fermi2EvenCharge;
      Last=Fermi2LastCharge;
    }
    else
    {
      //Identity    
      Fermi2Odd.setElement(1.0/2.0/this->N_fermions*ONE_c-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,0));
      Fermi2Even.setElement(1.0/2.0/this->N_fermions*ONE_c+0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,0));    
      Fermi2Last.setElement(1.0/2.0/this->N_fermions*ONE_c+sign*0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,0,0));
      
      //sigma_z
      Fermi2Odd.setElement(-0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,3));
      Fermi2Even.setElement(0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,10,3));    
      Fermi2Last.setElement(sign*0.5*this->mu_site[site_index-1]*taylor_factor,Indices(0,0,3));
      
      Odd=Fermi2Odd;
      Even=Fermi2Even;
      Last=Fermi2Last;
    }
    
    if(k==(this->N_total-1))
    {
      //Last site of fermionic species 2
      fermi_matrix=Last;
      fermi_matrix.reshape(Indices(D_bond*1,4));
      res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      Umpo.setOp(k,Local_op,true);       
    }
    else
    {
      //Rest of sites
      if((site_index%2)==0)
      {
	fermi_matrix=Even;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	Umpo.setOp(k,Local_op,true);       
      }
      else
      {
	fermi_matrix=Odd;
	fermi_matrix.reshape(Indices(D_bond*D_bond,4));
	res=reshape(fermi_matrix*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
	Local_op = new Operator(permute(res,Indices(3,1,4,2)));
	Umpo.setOp(k,Local_op,true);       
      }
    }    
  }
  
  //Fill the MPO with operators for the links
  saved=false;
  for(int k=2; k<this->N_total-2; k+=3)
  {
    //cout << "k=" << k << endl;
    if(saved)
    {
      Umpo.setOp(k,&Umpo.getOp(k-3),false);
    }
    else
    {
      //Save the gauge operator once
      res=reshape(Gauge*Gauge_operators,Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      Umpo.setOp(k,Local_op,true);       
      saved=true;
    }    
  } 
  
  cout << "Setting Matrix Elements done" << endl; 
  //Umpo.exportForMatlab("SU2UmpoStaticCharge.m");
}

void SU2Hamiltonian::removeGaugeInteraction()
{
  cout << "WARNING, will remove the gauge interaction, the Hamiltonian MPO has to be updated accordingly!" << endl;
  this->U00=identityMatrix(this->d_link);
  this->U01=identityMatrix(this->d_link);
  this->U10=identityMatrix(this->d_link);
  this->U11=identityMatrix(this->d_link);
  
  this->U00ad=identityMatrix(this->d_link);
  this->U01ad=identityMatrix(this->d_link);
  this->U10ad=identityMatrix(this->d_link);
  this->U11ad=identityMatrix(this->d_link);  
}


void SU2Hamiltonian::getSpin1MPO(MPO &Spin, int index, Component comp) const
{
  //Some error checking
  if(index<1 || index>this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getSpin1MPO, index is invalid" << endl;
    exit(666);
  }
  //Call the more general routine with appropriately transformed index
  switch(comp)
  {
    case x_comp:
      getSingleBodyMPO(Spin,3*index-3,sxOp);	break;
    case y_comp:
      getSingleBodyMPO(Spin,3*index-3,syOp);	break;
    case z_comp:
      getSingleBodyMPO(Spin,3*index-3,szOp);	break;
  }
}

void SU2Hamiltonian::getSpin2MPO(MPO &Spin, int index, Component comp) const
{
  //Some error checking
  if(index<1 || index>this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getSpin2MPO, index is invalid" << endl;
    exit(666);
  }
  //Call the more general routine with appropriately transformed index
  switch(comp)
  {
    case x_comp:
      getSingleBodyMPO(Spin,3*index-2,sxOp);	break;
    case y_comp:
      getSingleBodyMPO(Spin,3*index-2,syOp);	break;
    case z_comp:
      getSingleBodyMPO(Spin,3*index-2,szOp);	break;
  }
      
}

void SU2Hamiltonian::getLMPO(MPO &Spin, int index, Component comp) const
{
  //Some error checking
  if(index<1 || index>(this->N_fermions-1))
  {
    cout << "Error in SU2Hamiltonian::getLMPO, index is invalid" << endl;
    exit(666);
  }
  //Call the more general routine with appropriately transformed index
  switch(comp)
  {
    case x_comp:
      getSingleBodyMPO(Spin,3*index-1,LxOp);	break;
    case y_comp:
      getSingleBodyMPO(Spin,3*index-1,LyOp);	break;
    case z_comp:
      getSingleBodyMPO(Spin,3*index-1,LzOp);	break;
  }
}
  
void SU2Hamiltonian::getRMPO(MPO &Spin, int index, Component comp) const
{
  //Some error checking
  if(index<1 || index>(this->N_fermions-1))
  {
    cout << "Error in SU2Hamiltonian::getRMPO, index is invalid" << endl;
    exit(666);
  }
  //Call the more general routine with appropriately transformed index
  switch(comp)
  {
    case x_comp:
      getSingleBodyMPO(Spin,3*index-1,RxOp);	break;
    case y_comp:
      getSingleBodyMPO(Spin,3*index-1,RyOp);	break;
    case z_comp:
      getSingleBodyMPO(Spin,3*index-1,RzOp);	break;
  }
}
  
void SU2Hamiltonian::getJsquaredMPO(MPO &Spin, int index) const
{
  //Some error checking
  if(index<1 || index>(this->N_fermions-1))
  {
    cout << "Error in SU2Hamiltonian::getJsquaredMPO, index is invalid" << endl;
    exit(666);
  }
  //Call the more general routine with appropriately transformed index
  getSingleBodyMPO(Spin,3*index-1,JsqOp);
}

void SU2Hamiltonian::getSingleBodyMPO(MPO &Spin, int index, SingleBodyOperator Op) const
{
  //Frist clear the MPO and resize it
  Spin.clear();
  Spin.initLength(this->N_total);
  
  //Prepare the matrices
  mwArray Link_id=this->id_link_mpo;
  mwArray Link_op;
  mwArray Fermi_id=this->id_fermi_mpo;
  mwArray Fermi_op;
  
  switch(Op)
  {
    case sxOp:
      Fermi_op=this->sigma_x; Fermi_op.reshape(Indices(this->d_fermi,1,this->d_fermi,1)); break;
    case syOp:
      Fermi_op=this->sigma_y; Fermi_op.reshape(Indices(this->d_fermi,1,this->d_fermi,1)); break;
    case szOp:
      Fermi_op=this->sigma_z; Fermi_op.reshape(Indices(this->d_fermi,1,this->d_fermi,1)); break;
    case spOp:
      Fermi_op=this->sigma_plus; Fermi_op.reshape(Indices(this->d_fermi,1,this->d_fermi,1)); break;
    case smOp:
      Fermi_op=this->sigma_minus; Fermi_op.reshape(Indices(this->d_fermi,1,this->d_fermi,1)); break;
    case LxOp:
      Link_op=this->Lx; Link_op.reshape(Indices(this->d_link,1,this->d_link,1)); break;
    case LyOp:
      Link_op=this->Ly; Link_op.reshape(Indices(this->d_link,1,this->d_link,1)); break;
    case LzOp:
      Link_op=this->Lz; Link_op.reshape(Indices(this->d_link,1,this->d_link,1)); break;
    case RxOp:
      Link_op=this->Rx; Link_op.reshape(Indices(this->d_link,1,this->d_link,1)); break;
    case RyOp:
      Link_op=this->Ry; Link_op.reshape(Indices(this->d_link,1,this->d_link,1)); break;
    case RzOp:
      Link_op=this->Rz; Link_op.reshape(Indices(this->d_link,1,this->d_link,1)); break;
    case JsqOp:
      Link_op=this->Jsquared; Link_op.reshape(Indices(this->d_link,1,this->d_link,1)); break;
    default:
      cout << "Unknown operator" << endl;
      exit(666);
  }
  
  //Variable for the Operators
  Operator* Local_op; 
  
  //Now set operators for fermionic species 1
  bool fermi_saved=false;
  int saved_site;
  for(int k=0; k<(this->N_total-1); k+=3)
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
  
  //Now set operators for fermionic species 2
  fermi_saved=false;
  for(int k=1; k<this->N_total; k+=3)
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
  for(int k=2; k<(this->N_total-2); k+=3)
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

void SU2Hamiltonian::getGaussLawMPO(MPO &GaussLaw, Component comp) const
{
  switch(comp){
    case z_comp:
      getGaussLawMPOz(GaussLaw); break;
    default:
      getGaussLawMPOxy(GaussLaw,comp);
  }
}

void SU2Hamiltonian::getGaussLawMPOxy(MPO &GaussLaw, Component comp) const
{
  //Prepare the MPO
  GaussLaw.initLength(this->N_total);
  
  //Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(4,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      if(comp==x_comp)
      {
	//Lx
	Gauge_operators.setElement(this->Lx.getElement(Indices(i,j)),Indices(1,i,j));
	//Rx
	Gauge_operators.setElement(this->Rx.getElement(Indices(i,j)),Indices(2,i,j));
      }
      else if(comp==y_comp)
      {
	//Ly
	Gauge_operators.setElement(this->Ly.getElement(Indices(i,j)),Indices(1,i,j));
	//Ry
	Gauge_operators.setElement(this->Ry.getElement(Indices(i,j)),Indices(2,i,j));
      }
      //Jsquared
      Gauge_operators.setElement(this->Jsquared.getElement(Indices(i,j)),Indices(3,i,j));
     }
  }   
  //cout << "Gauge_operators" << Gauge_operators <<endl;
  Gauge_operators.reshape(Indices(4,this->d_link*this->d_link));  
  
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
  
  //Now prepare the matrices
  int D_bond=8;
  mwArray FermiFirst,Fermi1,Fermi2,Fermi1Last,Fermi2Last,Gauge;
  
  FermiFirst=mwArray(Indices(1,D_bond,4));	FermiFirst.fillWithZero();
  Fermi1=mwArray(Indices(D_bond,D_bond,4));	Fermi1.fillWithZero();
  Fermi2=mwArray(Indices(D_bond,D_bond,4));	Fermi2.fillWithZero();
  Fermi1Last=mwArray(Indices(D_bond,D_bond,4));	Fermi1Last.fillWithZero();
  Fermi2Last=mwArray(Indices(D_bond,1,4));	Fermi2Last.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,4));	Gauge.fillWithZero();
  
  //Identity
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(1.0/16.0*ONE_c,Indices(0,7,0));
  
  Fermi1.setElement(ONE_c,Indices(0,0,0));
  Fermi1.setElement(1.0/16.0*ONE_c,Indices(0,7,0));
  Fermi1.setElement(ONE_c,Indices(4,4,0));
  Fermi1.setElement(ONE_c,Indices(7,7,0));
  
  Fermi2.setElement(ONE_c,Indices(0,0,0));
  Fermi2.setElement(1.0/16.0*ONE_c,Indices(0,7,0));
  Fermi2.setElement(ONE_c,Indices(4,4,0));
  Fermi2.setElement(ONE_c,Indices(7,7,0));
  
  Fermi1Last.setElement(ONE_c,Indices(0,0,0));
  Fermi1Last.setElement(1.0/16.0*ONE_c,Indices(0,7,0));
  Fermi1Last.setElement(ONE_c,Indices(4,4,0));
  Fermi1Last.setElement(ONE_c,Indices(7,7,0));
  
  Fermi2Last.setElement(1.0/16.0*ONE_c,Indices(0,0,0));
  Fermi2Last.setElement(ONE_c,Indices(7,0,0));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(7,7,0));
  
  //sigma_plus
  FermiFirst.setElement(ONE_c,Indices(0,1,1));
  
  Fermi1.setElement(ONE_c,Indices(0,1,1));
  Fermi1.setElement(ONE_c,Indices(4,5,1));  
  
  Fermi2.setElement(ONE_c,Indices(2,2,1));
  if(comp==x_comp)
    Fermi2.setElement(I_c,Indices(6,7,1));
  else if(comp==y_comp)
    Fermi2.setElement(-1.0*ONE_c,Indices(6,7,1));
  
  Fermi1Last.setElement(ONE_c,Indices(4,5,1));  
  
  if(comp==x_comp)
    Fermi2Last.setElement(I_c,Indices(6,0,1)); 
  else if(comp==y_comp)
    Fermi2Last.setElement(-1.0*ONE_c,Indices(6,0,1)); 
  
  //sigma_minus
  FermiFirst.setElement(ONE_c,Indices(0,2,2));
  
  Fermi1.setElement(ONE_c,Indices(0,2,2));
  Fermi1.setElement(ONE_c,Indices(4,6,2));  
  
  Fermi2.setElement(ONE_c,Indices(1,1,2));
  if(comp==x_comp)
    Fermi2.setElement(-1.0*I_c,Indices(5,7,2));
  else if(comp==y_comp)
    Fermi2.setElement(-1.0*ONE_c,Indices(5,7,2));
  
  Fermi1Last.setElement(ONE_c,Indices(4,6,2));  
  
  if(comp==x_comp)
    Fermi2Last.setElement(-1.0*I_c,Indices(5,0,2));
  else if(comp==y_comp)
    Fermi2Last.setElement(-1.0*ONE_c,Indices(5,0,2));
  
  //sigma_z
  FermiFirst.setElement(-1.0/8.0*ONE_c,Indices(0,3,3));
  
  Fermi1.setElement(-1.0/8.0*ONE_c,Indices(0,3,3));
  
  Fermi2.setElement(ONE_c,Indices(3,7,3));
  
  Fermi1Last.setElement(-1.0/8.0*ONE_c,Indices(0,3,3));
  
  Fermi2Last.setElement(ONE_c,Indices(3,0,3));
  
  //Lx
  if(comp==x_comp)
  {
    Gauge.setElement(I_c,Indices(1,7,1));
    Gauge.setElement(-1.0*I_c,Indices(2,7,1));
  }
  else if(comp==y_comp)
  {
    Gauge.setElement(ONE_c,Indices(1,7,1));
    Gauge.setElement(ONE_c,Indices(2,7,1));
  }
  Gauge.setElement(-2.0*ONE_c,Indices(4,7,1));
  
  //Rx
  Gauge.setElement(ONE_c,Indices(0,4,2));
  
  //Jsquared
  Gauge.setElement(2.0/3.0*ONE_c,Indices(0,7,3));  
  
  FermiFirst.reshape(Indices(1*D_bond,4));
  Fermi1.reshape(Indices(D_bond*D_bond,4));
  Fermi2.reshape(Indices(D_bond*D_bond,4));
  Fermi1Last.reshape(Indices(D_bond*D_bond,4));
  Fermi2Last.reshape(Indices(D_bond*1,4));
  Gauge.reshape(Indices(D_bond*D_bond,4));
  
  //Variable for the Operators
  Operator* Local_op; 
  mwArray res;
  
  //Now set operators for fermionic species 1
  bool fermi_saved=false;
  int saved_site;
  for(int k=0; k<(this->N_total-1); k+=3)
  {
    if(k==0)
    {
      res=reshape(FermiFirst*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      GaussLaw.setOp(k,Local_op,true);
    }
    else if(k==(this->N_total-2))
    {
      res=reshape(Fermi1Last*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      GaussLaw.setOp(k,Local_op,true);
    }
    else if(fermi_saved)
      GaussLaw.setOp(k,&GaussLaw.getOp(k-3),false);      
    else
    {
      res=reshape(Fermi1*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      GaussLaw.setOp(k,Local_op,true);
      fermi_saved=true;
    }
  }
  
  //Now set operators for fermionic species 2
  fermi_saved=false;
  for(int k=1; k<this->N_total; k+=3)
  {
    if(k==(this->N_total-1))
    {
      res=reshape(Fermi2Last*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      GaussLaw.setOp(k,Local_op,true);
    }
    else if(fermi_saved)
      GaussLaw.setOp(k,&GaussLaw.getOp(k-3),false);
    else
    {
      res=reshape(Fermi2*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      GaussLaw.setOp(k,Local_op,true);      
    }
  } 
  
  //Now set the link operators
  bool link_saved=false;
  for(int k=2; k<(this->N_total-2); k+=3)
  {
    if(link_saved)
      GaussLaw.setOp(k,&GaussLaw.getOp(k-3),false);
    else
    {
      res=reshape(Gauge*Gauge_operators,Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      GaussLaw.setOp(k,Local_op,true);
      link_saved=true;
    }
  }
  
  /*if(comp==x_comp)
    GaussLaw.exportForMatlab("Gx.m");
  else
    GaussLaw.exportForMatlab("Gy.m");*/
}

void SU2Hamiltonian::getGaussLawMPOz(MPO &GaussLaw) const
{
  //Prepare the MPO
  GaussLaw.initLength(this->N_total);
  
  //Array to store the operators acting on a link
  mwArray Gauge_operators(Indices(4,this->d_link,this->d_link));
  Gauge_operators.fillWithZero();   
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //Lz
      Gauge_operators.setElement(this->Lz.getElement(Indices(i,j)),Indices(1,i,j));
      //Rz
      Gauge_operators.setElement(this->Rz.getElement(Indices(i,j)),Indices(2,i,j));
      //Jsquared
      Gauge_operators.setElement(this->Jsquared.getElement(Indices(i,j)),Indices(3,i,j));
     }
  }   
  //cout << "Gauge_operators" << Gauge_operators <<endl;
  Gauge_operators.reshape(Indices(4,this->d_link*this->d_link));  
  
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
  
  //Now prepare the matrices
  int D_bond=5;
  mwArray FermiFirst,Fermi1,Fermi2,Fermi2Last,Gauge;
  
  FermiFirst=mwArray(Indices(1,D_bond,2));	FermiFirst.fillWithZero();
  Fermi1=mwArray(Indices(D_bond,D_bond,2));	Fermi1.fillWithZero();
  Fermi2=mwArray(Indices(D_bond,D_bond,2));	Fermi2.fillWithZero();
  Fermi2Last=mwArray(Indices(D_bond,1,2));	Fermi2Last.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,4));	Gauge.fillWithZero();
  
  //Identity
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(1.0/16.0*ONE_c,Indices(0,4,0));
  
  Fermi1.setElement(ONE_c,Indices(0,0,0));
  Fermi1.setElement(1.0/16.0*ONE_c,Indices(0,4,0));
  Fermi1.setElement(ONE_c,Indices(3,3,0));
  Fermi1.setElement(ONE_c,Indices(4,4,0));
  
  Fermi2.setElement(ONE_c,Indices(0,0,0));
  Fermi2.setElement(1.0/16.0*ONE_c,Indices(0,4,0));
  Fermi2.setElement(ONE_c,Indices(1,1,0));
  Fermi2.setElement(ONE_c,Indices(3,3,0));
  Fermi2.setElement(ONE_c,Indices(4,4,0));
  
  Fermi2Last.setElement(1.0/16.0*ONE_c,Indices(0,0,0));
  Fermi2Last.setElement(ONE_c,Indices(4,0,0));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(4,4,0));
    
  //sigma_z
  FermiFirst.setElement(ONE_c,Indices(0,1,1));
  
  Fermi1.setElement(ONE_c,Indices(0,1,1));
  Fermi1.setElement(0.5*ONE_c,Indices(3,4,1));
  
  Fermi2.setElement(ONE_c,Indices(0,2,1));
  Fermi2.setElement(-1.0/8.0*ONE_c,Indices(1,4,1));
  Fermi2.setElement(-1.0/2.0*ONE_c,Indices(3,4,1));
  
  Fermi2Last.setElement(-1.0/8.0*ONE_c,Indices(1,0,1));
  Fermi2Last.setElement(-1.0/2.0*ONE_c,Indices(3,0,1));
  
  //Lz
  Gauge.setElement(-0.5*ONE_c,Indices(1,4,1));
  Gauge.setElement(0.5*ONE_c,Indices(2,4,1));
  Gauge.setElement(-2.0*ONE_c,Indices(3,4,1));
  
  //Rz
  Gauge.setElement(ONE_c,Indices(0,3,2));
  
  //Jsquared
  Gauge.setElement(2.0/3.0*ONE_c,Indices(0,4,3));  
  
  FermiFirst.reshape(Indices(1*D_bond,2));
  Fermi1.reshape(Indices(D_bond*D_bond,2));
  Fermi2.reshape(Indices(D_bond*D_bond,2));
  Fermi2Last.reshape(Indices(D_bond*1,2));
  Gauge.reshape(Indices(D_bond*D_bond,4));
  
  //Variable for the Operators
  Operator* Local_op; 
  mwArray res;
  
  //Now set operators for fermionic species 1
  bool fermi_saved=false;
  int saved_site;
  for(int k=0; k<(this->N_total-1); k+=3)
  {
    if(k==0)
    {
      res=reshape(FermiFirst*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      GaussLaw.setOp(k,Local_op,true);
    }
    else if(fermi_saved)
      GaussLaw.setOp(k,&GaussLaw.getOp(k-3),false);      
    else
    {
      res=reshape(Fermi1*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      GaussLaw.setOp(k,Local_op,true);
      fermi_saved=true;
    }
  }
  
  //Now set operators for fermionic species 2
  fermi_saved=false;
  for(int k=1; k<this->N_total; k+=3)
  {
    if(k==(this->N_total-1))
    {
      res=reshape(Fermi2Last*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      GaussLaw.setOp(k,Local_op,true);
    }
    else if(fermi_saved)
      GaussLaw.setOp(k,&GaussLaw.getOp(k-3),false);
    else
    {
      res=reshape(Fermi2*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      GaussLaw.setOp(k,Local_op,true);      
    }
  } 
  
  //Now set the link operators
  bool link_saved=false;
  for(int k=2; k<(this->N_total-2); k+=3)
  {
    if(link_saved)
      GaussLaw.setOp(k,&GaussLaw.getOp(k-3),false);
    else
    {
      res=reshape(Gauge*Gauge_operators,Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      GaussLaw.setOp(k,Local_op,true);
      link_saved=true;
    }
  }
  //GaussLaw.exportForMatlab("Gz.m");
}


void SU2Hamiltonian::constructInitialMPS(MPS &mps, InitState state) const
{
  //Adjust the physical dimensions of the MPS, therefore overwrite the old one and create soemthing with matching dimensions
  vector<int> phys_dims(this->N_total,this->d_fermi);  
  for(int k=2; k<(this->N_total-2); k+=3)
    phys_dims[k]=this->d_link;    
  
  mps.clear();
  mps=MPS(this->N_total,1,phys_dims);
  
  
  mwArray spinsiteup(Indices(2,1,1));
  mwArray spinsitedown(Indices(2,1,1));
  mwArray link(Indices(this->d_link,1,1));
  mwArray link_flux(Indices(this->d_link,1,1));
  
  spinsiteup.fillWithZero();
  spinsitedown.fillWithZero();
  link.fillWithZero();
  
  spinsiteup.setElement(ONE_c,Indices(0,0,0));
  spinsitedown.setElement(ONE_c,Indices(1,0,0));
  link.setElement(ONE_c,Indices((this->d_link-1)/2,0,0));
  link_flux.setElement(ONE_c,Indices(0,0,0));
  //link.setElement(ONE_c,Indices(1,0,0));
  
  
  //Set the matrices
  int index;
  if(state==is_updown)
  {
    for(int k=0; k<(this->N_total-1); k+=3)
    {
      //Fermionic sites
      index=k/3+1;    
      if((index%2)!=0)
      {
	mps.setA(k,spinsiteup);
	mps.setA(k+1,spinsiteup);
      }
      else
      {
	mps.setA(k,spinsitedown);
	mps.setA(k+1,spinsitedown);
      }
      //Links
      if(k<(this->N_total-2))
	mps.setA(k+2,link);
    } 
  }
  else if(state==is_downup)
  {
    for(int k=0; k<(this->N_total-1); k+=3)
    {
      //Fermionic sites
      index=k/3+1;    
      if((index%2)==0)
      {
	mps.setA(k,spinsiteup);
	mps.setA(k+1,spinsiteup);
      }
      else
      {
	mps.setA(k,spinsitedown);
	mps.setA(k+1,spinsitedown);
      }
      //Links
      if(k<(this->N_total-2))
	mps.setA(k+2,link);
    } 
  }
  else if(state==is_down)
  {
    for(int k=0; k<(this->N_total-1); k+=3)
    {
      //Fermionic sites
      mps.setA(k,spinsitedown);
      mps.setA(k+1,spinsitedown);
      //Links
      if(k<(this->N_total-2))
	mps.setA(k+2,link);
    } 
  }
  else if(state==is_up)
  {
    for(int k=0; k<(this->N_total-1); k+=3)
    {
      //Fermionic sites
      mps.setA(k,spinsiteup);
      mps.setA(k+1,spinsiteup);
      //Links
      if(k<(this->N_total-2))
	mps.setA(k+2,link);
    } 
  }
  else if(state==is_updown_onsite)
  {
    for(int k=0; k<(this->N_total-1); k+=3)
    {
      //Fermionic sites
      mps.setA(k,spinsiteup);
      mps.setA(k+1,spinsitedown);
      //Links
      if(k<(this->N_total-2))
	mps.setA(k+2,link);
    } 
  }
  else if(state==is_downup_onsite)
  {
    for(int k=0; k<(this->N_total-1); k+=3)
    {
      //Fermionic sites
      mps.setA(k,spinsitedown);
      mps.setA(k+1,spinsiteup);
      //Links
      if(k<(this->N_total-2))
	mps.setA(k+2,link);
    } 
  }
  else if(state==str_teststate_1)
  {
    for(int k=0; k<(this->N_total-1); k+=3)
    {
      //Fermionic sites
      mps.setA(k,spinsitedown);
      mps.setA(k+1,spinsiteup);
      //Links
      if(k<(this->N_total-2))
      {
	if((k+2)>=5 && (k+2)<=9)
	  mps.setA(k+2,link_flux);
	else
	  mps.setA(k+2,link);
      }
    } 
  }
  else if(state==str_teststate_2)
  {
    for(int k=0; k<(this->N_total-1); k+=3)
    {
      //Fermionic sites
      mps.setA(k,spinsitedown);
      mps.setA(k+1,spinsitedown);
      //Links
      if(k<(this->N_total-2))
      {
	if((k+2)>=5 && (k+2)<=9)
	  mps.setA(k+2,link_flux);
	else
	  mps.setA(k+2,link);
      }
    }
  }
    
  mps.gaugeCond('R',true);
}

void SU2Hamiltonian::constructProjectors(MPO& Podd, MPO& Peven) const
{
  //Resize the MPOs
  Podd.initLength(this->N_total);
  Peven.initLength(this->N_total);
  
  //Set up mwArrays
  mwArray LeftEdge(Indices(this->d_fermi*this->d_fermi*this->d_link,this->d_fermi*this->d_fermi*this->d_link));
  mwArray RightEdge(Indices(this->d_link*this->d_fermi*this->d_fermi,this->d_link*this->d_fermi*this->d_fermi));
  mwArray RegularSite(Indices(this->d_link*this->d_fermi*this->d_fermi*this->d_link,this->d_link*this->d_fermi*this->d_fermi*this->d_link));
  
  mwArray Preg1,Preg2,Preg3,Preg4,Pleft1,Pleft2,Pleft3,Pright1,Pright2,Pright3;
  
  buildProjectorMatrices(LeftEdge,RightEdge,RegularSite);
  
  /*cout << "LeftEdge: " << LeftEdge << endl;
  cout << "RightEdge: " << RightEdge << endl;
  cout << "RegularSite: " << RegularSite << endl;*/
  
  //Matrices for the SVDs
  mwArray U,S,V;
  Indices dim;
  
  
  /*****************************************************************************
   * Warning
   * ***************************************************************************
   * In the MPS code, the Kronecker product works like this
   *  1 2 3
   *  | | |
   * -x-x-x-
   *  | | |
   *  5 6 7
   * and the new indices are formed as i=(1 2 3) and j=(5 6 7). As I computed analytically with Maple the common subspace of the three Gauss Law operators, the ordering of the Kronecker Product is differen, as Maple has the same convetion as Matlab. The new indices are therefore i=(3 2 1) and j=(7 6 5). Therefore in this part of code the indexing is a bit different then in all the others
      
  
  /*****************************************************************************
   * Regular sites
   * ***************************************************************************/  
  //First permute the indices in the right order  
  RegularSite.reshape(Indices(this->d_link,this->d_fermi,this->d_fermi,this->d_link,this->d_link,this->d_fermi,this->d_fermi,this->d_link));
  RegularSite.permute(Indices(1,5,2,6,3,7,4,8));
  //Now take the first two indices together and do a SVD
  RegularSite.reshape(Indices(this->d_link*this->d_link,pow(this->d_fermi,4)*this->d_link*this->d_link));
  wrapper::svd(RegularSite,U,S,V);  
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape the matrix and save results
  U.reshape(Indices(this->d_link,1,this->d_link,-1));
  Preg1=U;
  RegularSite=V;
  
  //Proceed with the next SVD
  dim=RegularSite.getDimensions();
  RegularSite.reshape(Indices(dim[0]*this->d_fermi*this->d_fermi,this->d_fermi*this->d_fermi*this->d_link*this->d_link));
  wrapper::svd(RegularSite,U,S,V);  
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape the matrix and save results
  U.reshape(Indices(dim[0],this->d_fermi,this->d_fermi,-1));
  U.permute(Indices(2,1,3,4));
  Preg2=U;
  RegularSite=V;
  
  //Proceed with the next SVD
  dim=RegularSite.getDimensions();
  RegularSite.reshape(Indices(dim[0]*this->d_fermi*this->d_fermi,this->d_link*this->d_link));
  wrapper::svd(RegularSite,U,S,V);  
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape the matrix and save results
  U.reshape(Indices(dim[0],this->d_fermi,this->d_fermi,-1));
  U.permute(Indices(2,1,3,4));
  Preg3=U;
  V.reshape(Indices(-1,this->d_link,this->d_link,1));
  V.permute(Indices(2,1,3,4));
  Preg4=V;
  
  /*****************************************************************************
   * Left Edge
   * ***************************************************************************/
  //Now take the first two indices together and do a SVD  
  LeftEdge.reshape(Indices(this->d_fermi,this->d_fermi,this->d_link,this->d_fermi,this->d_fermi,this->d_link));
  LeftEdge.permute(Indices(1,4,2,5,3,6));
  LeftEdge.reshape(Indices(this->d_fermi*this->d_fermi,this->d_fermi*this->d_fermi*this->d_link*this->d_link));
  wrapper::svd(LeftEdge,U,S,V);  
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape the matrix and save results
  U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
  Pleft1=U;
  LeftEdge=V;
  
  //Proceed with the next SVD
  dim=LeftEdge.getDimensions();
  LeftEdge.reshape(Indices(dim[0]*this->d_fermi*this->d_fermi,this->d_link*this->d_link));
  wrapper::svd(LeftEdge,U,S,V);  
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape the matrix and save results
  U.reshape(Indices(dim[0],this->d_fermi,this->d_fermi,-1));
  U.permute(Indices(2,1,3,4));
  Pleft2=U;
  V.reshape(Indices(-1,this->d_link,this->d_link,1));
  V.permute(Indices(2,1,3,4));
  Pleft3=V;
  
  /*****************************************************************************
   * Right Edge
   * ***************************************************************************/
  //Now take the first two indices together and do a SVD  
  RightEdge.reshape(Indices(this->d_link,this->d_fermi,this->d_fermi,this->d_link,this->d_fermi,this->d_fermi));
  RightEdge.permute(Indices(1,4,2,5,3,6));
  RightEdge.reshape(Indices(this->d_link*this->d_link,pow(this->d_fermi,4)));
  wrapper::svd(RightEdge,U,S,V);  
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape the matrix and save results
  U.reshape(Indices(this->d_link,1,this->d_link,-1));
  Pright1=U;
  RightEdge=V;
  
  //Proceed with the next SVD
  dim=RightEdge.getDimensions();
  RightEdge.reshape(Indices(dim[0]*this->d_fermi*this->d_fermi,this->d_fermi*this->d_fermi));
  wrapper::svd(RightEdge,U,S,V);  
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape the matrix and save results
  U.reshape(Indices(dim[0],this->d_fermi,this->d_fermi,-1));
  U.permute(Indices(2,1,3,4));
  Pright2=U;
  V.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
  V.permute(Indices(2,1,3,4));
  Pright3=V;
  
    
  //Add the missing identities for Podd on the left edge
  Podd.setOp(0,new Operator(this->id_fermi_mpo),true);
  Podd.setOp(1,new Operator(this->id_fermi_mpo),true);
  
  //Add the missing operators for Peven on the left edge
  Peven.setOp(0,new Operator(Pleft1),true);
  Peven.setOp(1,new Operator(Pleft2),true);
  Peven.setOp(2,new Operator(Pleft3),true);  
    
  //Add the missing opertors for Podd on the right edge
  Podd.setOp(this->N_total-3,new Operator(Pright1),true);
  Podd.setOp(this->N_total-2,new Operator(Pright2),true);
  Podd.setOp(this->N_total-1,new Operator(Pright3),true);
  
  //Add the missing opertors for Peven on the left edge
  Peven.setOp(this->N_total-2,new Operator(this->id_fermi_mpo),true);
  Peven.setOp(this->N_total-1,new Operator(this->id_fermi_mpo),true);
  
  //For testing put only the operators in the center, and padd with identities
  /*Podd.setOp(0,new Operator(this->id_fermi_mpo),true);
  Podd.setOp(1,new Operator(this->id_fermi_mpo),true);
  
  //Add the missing operators for Peven on the left edge
  Peven.setOp(0,new Operator(this->id_fermi_mpo),true);
  Peven.setOp(1,new Operator(this->id_fermi_mpo),true);
  Peven.setOp(2,new Operator(this->id_link_mpo),true);  
    
  //Add the missing opertors for Podd on the right edge
  Podd.setOp(this->N_total-3,new Operator(this->id_link_mpo),true);
  Podd.setOp(this->N_total-2,new Operator(this->id_fermi_mpo),true);
  Podd.setOp(this->N_total-1,new Operator(this->id_fermi_mpo),true);
  
  //Add the missing opertors for Peven on the left edge
  Peven.setOp(this->N_total-2,new Operator(this->id_fermi_mpo),true);
  Peven.setOp(this->N_total-1,new Operator(this->id_fermi_mpo),true);*/
  
  //Little dummy MPO for error checking
  /*MPO dummy(4);
  dummy.setOp(0,new Operator(Preg1),true);
  dummy.setOp(1,new Operator(Preg2),true);
  dummy.setOp(2,new Operator(Preg3),true);
  dummy.setOp(3,new Operator(Preg4),true);
  
  dummy.exportForMatlab("DummyProjector.m");  */
  
  //Fill the rest of the projectors
  int index;
  for(int k=2; k<(this->N_total-5); k+=3)
  {
    //Transform k into link index
    index=(int)(((double) k)/3.0+1.0/3.0);
    if((index%2)==0)
    {
      Peven.setOp(k,new Operator(Preg1),true);
      Peven.setOp(k+1,new Operator(Preg2),true);
      Peven.setOp(k+2,new Operator(Preg3),true);
      Peven.setOp(k+3,new Operator(Preg4),true);
      
      Podd.setOp(k+1,new Operator(id_fermi_mpo),true);
      Podd.setOp(k+2,new Operator(id_fermi_mpo),true);
    }
    else
    {
      Podd.setOp(k,new Operator(Preg1),true);
      Podd.setOp(k+1,new Operator(Preg2),true);
      Podd.setOp(k+2,new Operator(Preg3),true);
      Podd.setOp(k+3,new Operator(Preg4),true);
      
      Peven.setOp(k+1,new Operator(id_fermi_mpo),true);
      Peven.setOp(k+2,new Operator(id_fermi_mpo),true);
    }
  }
  
  /*cout << "Podd " << Podd << endl;
  cout << "Peven " << Peven << endl;*/
  
  /*Podd.exportForMatlab("Podd.m");
  Peven.exportForMatlab("Peven.m");*/
}

void SU2Hamiltonian::constructProjectors(MPO& Podd, MPO& Peven, int pos, int leng, Component comp) const
{
  cout << "Constructing Projectors" << endl;
  
  //Resize the MPOs
  Podd.initLength(this->N_total);
  Peven.initLength(this->N_total);
  
  if(pos<0)
    pos = round(this->N_fermions/2)-round(leng/2);
  
  //Set up mwArrays
  mwArray PL_m_half,P_m_half,PR_m_half,PL_zero,P_zero,PR_zero,PL_p_half,P_p_half,PR_p_half;  
  mwArray Preg1,Preg2,Preg3,Preg4,Pleft1,Pleft2,Pleft3,Pright1,Pright2,Pright3,ChargeL1,ChargeL2,ChargeL3,ChargeL4,ChargeR1,ChargeR2,ChargeR3,ChargeR4;
  mwArray Left,Right,Regular,Charge1,Charge2;
  
  buildProjectorMatrices(PL_m_half,P_m_half,PR_m_half,PL_zero,P_zero,PR_zero,PL_p_half,P_p_half,PR_p_half,comp);
  
  //Check if the matrices are really projectors, meaning that the spectrum consists of 1 and 0 only
  /*vector<complex_t> eval;
  mwArray evec;
  wrapper::eig(PL_m_half,eval,evec);
  cout << "Eigenvalues PL_m_half: " << eval << endl;
  wrapper::eig(P_m_half,eval,evec);
  cout << "Eigenvalues P_m_half: " << eval << endl;
  wrapper::eig(PR_m_half,eval,evec);
  cout << "Eigenvalues PR_m_half: " << eval << endl;
  
  wrapper::eig(PL_p_half,eval,evec);
  cout << "Eigenvalues PL_p_half: " << eval << endl;
  wrapper::eig(P_p_half,eval,evec);
  cout << "Eigenvalues P_p_half: " << eval << endl;
  wrapper::eig(PR_p_half,eval,evec);
  cout << "Eigenvalues PR_p_half: " << eval << endl;
  
  wrapper::eig(PL_zero,eval,evec);
  cout << "Eigenvalues PL_zero: " << eval << endl;
  wrapper::eig(P_zero,eval,evec);
  cout << "Eigenvalues P_zero: " << eval << endl;
  wrapper::eig(PR_zero,eval,evec);
  cout << "Eigenvalues PR_zero: " << eval << endl;*/
  
  Regular=P_zero;
  if((pos%2)==0)
  {
    cout << "Even start pos " << pos << endl;
    //Even start site, meaning charges come in + - configuration
    Left=PL_zero;
    Charge1=P_p_half;
    if((pos+leng)==this->N_fermions)
    {
      //Case that the right charge is on the end
      Right=PR_m_half;
      Charge2=P_zero;
    }
    else
    {
      Right=PR_zero;
      Charge2=P_m_half;
    }
  }
  else
  {
    cout << "Odd start pos " << pos << endl;
    //Odd start site, meaning charges come in - +  configuration
    if(pos==1)
    {
      //Case that the left charge is on first position
      Left=PL_m_half;
      Charge1=P_zero;
    }
    else
    {
      Left=PL_zero;
      Charge1=P_m_half;
    }
    if((pos+leng)==this->N_fermions)
    {
      //Case that the right charge is on the end
      Right=PR_p_half;
      Charge2=P_zero;
    }
    else
    {
      Right=PR_zero;
      Charge2=P_p_half;
    }
  }
  
  /*cout << "LeftEdge: " << LeftEdge << endl;
  cout << "RightEdge: " << RightEdge << endl;
  cout << "RegularSite: " << RegularSite << endl;*/ 
  
  //Prepare vectors for output
  vector <mwArray> LeftM, RightM, RegularM, Charge1M, Charge2M;
  vector <int> dimsleft, dimsright, dims;
  
  dimsleft.push_back(this->d_fermi);	dimsleft.push_back(this->d_fermi);	dimsleft.push_back(this->d_link);
  dimsright.push_back(this->d_link);	dimsright.push_back(this->d_fermi);	dimsright.push_back(this->d_fermi);
  dims.push_back(this->d_link);		dims.push_back(this->d_fermi);		dims.push_back(this->d_fermi);		dims.push_back(this->d_link);  
  
  /*****************************************************************************
   * Warning
   * ***************************************************************************
   * In the MPS code, the Kronecker product works like this
   *  1 2 3
   *  | | |
   * -x-x-x-
   *  | | |
   *  5 6 7
   * and the new indices are formed as i=(1 2 3) and j=(5 6 7).      
   */
  
  //Regular sites
  this->splitLocalOp(Regular,dims,RegularM);    
  // Left Edge
  this->splitLocalOp(Left,dimsleft,LeftM);  
  //Right Edge
  this->splitLocalOp(Right,dimsright,RightM);  
  //Charge 1
  this->splitLocalOp(Charge1,dims,Charge1M);   
  //Charge 2
  this->splitLocalOp(Charge2,dims,Charge2M); 
  
  /*MPO dummy(3);  
  dummy.setOp(0,new Operator(LeftM[0]),true);
  dummy.setOp(1,new Operator(LeftM[1]),true);
  dummy.setOp(2,new Operator(LeftM[2]),true);
  
  dummy.exportForMatlab("Left.m");
  
  dummy.clear();
  dummy.initLength(3);
  dummy.setOp(0,new Operator(RightM[0]),true);
  dummy.setOp(1,new Operator(RightM[1]),true);
  dummy.setOp(2,new Operator(RightM[2]),true);
  
  dummy.exportForMatlab("Right.m");
  
  dummy.clear();
  dummy.initLength(4);
  dummy.setOp(0,new Operator(RegularM[0]),true);
  dummy.setOp(1,new Operator(RegularM[1]),true);
  dummy.setOp(2,new Operator(RegularM[2]),true);
  dummy.setOp(3,new Operator(RegularM[3]),true);
  
  dummy.exportForMatlab("Regular.m");*/
  
  //Set left and right edges as well as the charge operators
  Podd.setOp(0,new Operator(this->id_fermi_mpo), true);
  Podd.setOp(1,new Operator(this->id_fermi_mpo), true);
  
  Peven.setOp(0,new Operator(LeftM[0]), true);
  Peven.setOp(1,new Operator(LeftM[1]), true);
  Peven.setOp(2,new Operator(LeftM[2]), true);
  
  Podd.setOp(this->N_total-3,new Operator(RightM[0]),true);
  Podd.setOp(this->N_total-2,new Operator(RightM[1]),true);
  Podd.setOp(this->N_total-1,new Operator(RightM[2]),true);  
  
  Peven.setOp(this->N_total-2,new Operator(this->id_fermi_mpo),true);
  Peven.setOp(this->N_total-1,new Operator(this->id_fermi_mpo),true);
  
  //Fill the rest of the projectors
  int index;
  for(int k=2; k<(this->N_total-5); k+=3)
  {
    //Transform k into link index
    index=(int)(((double) k)/3.0+1.0/3.0);
    if((index%2)==0)
    {
      if(index==(pos-1))
      {
	Peven.setOp(k,new Operator(Charge1M[0]),true);
	Peven.setOp(k+1,new Operator(Charge1M[1]),true);
	Peven.setOp(k+2,new Operator(Charge1M[2]),true);
	Peven.setOp(k+3,new Operator(Charge1M[3]),true);
      }
      else if(index==(pos+leng-1))
      {
	Peven.setOp(k,new Operator(Charge2M[0]),true);
	Peven.setOp(k+1,new Operator(Charge2M[1]),true);
	Peven.setOp(k+2,new Operator(Charge2M[2]),true);
	Peven.setOp(k+3,new Operator(Charge2M[3]),true);
      }
      else
      {
	Peven.setOp(k,new Operator(RegularM[0]),true);
	Peven.setOp(k+1,new Operator(RegularM[1]),true);
	Peven.setOp(k+2,new Operator(RegularM[2]),true);
	Peven.setOp(k+3,new Operator(RegularM[3]),true);
      }
      
      Podd.setOp(k+1,new Operator(id_fermi_mpo),true);
      Podd.setOp(k+2,new Operator(id_fermi_mpo),true);
    }
    else
    {
      if(index==(pos-1))
      {
	Podd.setOp(k,new Operator(Charge1M[0]),true);
	Podd.setOp(k+1,new Operator(Charge1M[1]),true);
	Podd.setOp(k+2,new Operator(Charge1M[2]),true);
	Podd.setOp(k+3,new Operator(Charge1M[3]),true);
      }
      else if(index==(pos+leng-1))
      {
	Podd.setOp(k,new Operator(Charge2M[0]),true);
	Podd.setOp(k+1,new Operator(Charge2M[1]),true);
	Podd.setOp(k+2,new Operator(Charge2M[2]),true);
	Podd.setOp(k+3,new Operator(Charge2M[3]),true);
      }
      else
      {
	Podd.setOp(k,new Operator(RegularM[0]),true);
	Podd.setOp(k+1,new Operator(RegularM[1]),true);
	Podd.setOp(k+2,new Operator(RegularM[2]),true);
	Podd.setOp(k+3,new Operator(RegularM[3]),true);
      }
      
      Peven.setOp(k+1,new Operator(id_fermi_mpo),true);
      Peven.setOp(k+2,new Operator(id_fermi_mpo),true);
    }
  }
  
  /*cout << "Podd " << Podd << endl;
  cout << "Peven " << Peven << endl;*/
  
  /*Podd.exportForMatlab("Podd.m");
  Peven.exportForMatlab("Peven.m");*/
}

void SU2Hamiltonian::createString(MPS& mps, Contractor &contractor,int D, unsigned int length, int pos)
{
  //Determine the starting position of the string
  if(pos<0)
    pos = round(this->N_fermions/2)-round(length/2);
  
  //TODO: Some more error checking?
  if(mps.getLength()!= this->N_total)
  {
    cout << "Error in SU2Hamiltonian::createString, total number of sites is " << this->N_total << ", length of MPS is "<< mps.getLength() << endl;
    exit(666);
  }
  if(pos>this->N_total || pos<0 || (pos+length)>this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::createString, desired string of length " << length<< " starting on position " << pos << " is outside the the chain of " << this->N_fermions << " spins" << endl;
    exit(666);
  }
  if((length%2)==0)
  {
    cout << "Error in SU2Hamiltonian::createString, only odd values for length are supported, received " << length << ", which is even" << endl;
    exit(666);
  }
 
  //MPOs to create the strings
  MPO StringOdd(1),StringEven(1);
  getStringMPO(StringOdd,StringEven,pos,pos+length);
  
  MPS aux=MPS(mps);
  aux.increaseBondDimension(D);
  aux.setRandomState();  
  
  if((pos%2)==0)
  {
    cout << "Applying String Even" << endl;    
    aux.approximate(StringEven,mps,D);
    contractor.optimize(StringEven,mps,aux,D);      
    aux.setNormFact(1.0);
    cout << "Applying String Odd" << endl;
    mps.approximate(StringOdd,aux,D);
    contractor.optimize(StringOdd,aux,mps,D);
  }
  else
  {
    //Starting with the even sites, therefore apply StringOdd first
    cout << "Applying String Odd" << endl;
    aux.approximate(StringOdd,mps,D);
    contractor.optimize(StringOdd,mps,aux,D);
    cout << "Applying String Even" << endl;
    mps.approximate(StringEven,aux,D);
    contractor.optimize(StringEven,aux,mps,D);
  } 
  //Restore Norm
  //mps.setNormFact(1.0);
  mps.gaugeCond('R',true);
}

void SU2Hamiltonian::getStringMPO(MPO &Podd, MPO &Peven, int startpos, int endpos) const
{  
  //Initialize the MPOs
  Podd.initLength(this->N_total);
  Peven.initLength(this->N_total);
 
  //Determine whether my startposition is a even or odd fermionic index
  bool even_start_site=((startpos%2)==0);
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(4,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(5,this->d_link,this->d_link));	Gauge_operators.fillWithZero();  
  
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
  
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      if(even_start_site)
      {
	//U00
	Gauge_operators.setElement(this->U00.getElement(Indices(i,j)),Indices(1,i,j));
	//U01
	Gauge_operators.setElement(this->U01.getElement(Indices(i,j)),Indices(2,i,j));
	//U10
	Gauge_operators.setElement(this->U10.getElement(Indices(i,j)),Indices(3,i,j));
	//U11
	Gauge_operators.setElement(this->U11.getElement(Indices(i,j)),Indices(4,i,j));
      }
      else
      {      
	//U00 dagger
	Gauge_operators.setElement(this->U00ad.getElement(Indices(i,j)),Indices(1,i,j));
	//U01 dagger
	Gauge_operators.setElement(this->U01ad.getElement(Indices(i,j)),Indices(2,i,j));
	//U10 dagger
	Gauge_operators.setElement(this->U10ad.getElement(Indices(i,j)),Indices(3,i,j));
	//U11 dagger
	Gauge_operators.setElement(this->U11ad.getElement(Indices(i,j)),Indices(4,i,j));     
      }
     }
  }   
  Gauge_operators.reshape(Indices(5,this->d_link*this->d_link));  
  
  //Determine wether I start at an odd or even position, depending on that prepare the matrices  
  mwArray Fermi1left,Fermi1,Fermi2,Fermi2right,Gauge;
  int D_bond=8;
  
  Fermi1left=mwArray(Indices(1,D_bond,4));		Fermi1left.fillWithZero();
  Fermi1=mwArray(Indices(D_bond,D_bond,4));		Fermi1.fillWithZero();
  Fermi2=mwArray(Indices(D_bond,D_bond,4));		Fermi2.fillWithZero();
  Fermi2right=mwArray(Indices(D_bond,1,4));		Fermi2right.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,5));		Gauge.fillWithZero();  
  
  if(even_start_site)
  {
    cout << "Starting at even site" << endl;
    //Fermi1left
    Fermi1left.setElement(ONE_c,Indices(0,0,0));
    
    Fermi1left.setElement(ONE_c,Indices(0,1,1));
    
    //Fermi1
    Fermi1.setElement(ONE_c,Indices(0,0,0));
    Fermi1.setElement(ONE_c,Indices(7,7,0));
    
    Fermi1.setElement(ONE_c,Indices(0,1,1));
    
    Fermi1.setElement(ONE_c,Indices(2,7,2));
    Fermi1.setElement(ONE_c,Indices(5,7,2));
    
    Fermi1.setElement(ONE_c,Indices(3,3,3));
    Fermi1.setElement(ONE_c,Indices(6,6,3));
    
    //Fermi2
    Fermi2.setElement(ONE_c,Indices(0,0,0));
    Fermi2.setElement(ONE_c,Indices(7,7,0));
    
    Fermi2.setElement(ONE_c,Indices(0,4,1));
    
    Fermi2.setElement(ONE_c,Indices(3,7,2));
    Fermi2.setElement(ONE_c,Indices(6,7,2));
    
    Fermi2.setElement(ONE_c,Indices(1,1,3));
    
    //Fermi2right
    Fermi2right.setElement(ONE_c,Indices(7,0,0));
      
    Fermi2right.setElement(ONE_c,Indices(3,0,2));
    Fermi2right.setElement(ONE_c,Indices(6,0,2));   
    
    //Gauge
    Gauge.setElement(ONE_c,Indices(0,0,0));
    Gauge.setElement(ONE_c,Indices(7,7,0));
    
    Gauge.setElement(ONE_c,Indices(1,2,1));
    
    Gauge.setElement(I_c,Indices(1,3,2));
    
    Gauge.setElement(-1.0*I_c,Indices(4,5,3));
    
    Gauge.setElement(ONE_c,Indices(4,6,4));
   
  }
  else
  {    
    cout << "Starting at even site" << endl;
    //Fermi1left
    Fermi1left.setElement(ONE_c,Indices(0,0,0));
    
    Fermi1left.setElement(ONE_c,Indices(0,1,2));
    
    //Fermi1
    Fermi1.setElement(ONE_c,Indices(0,0,0));
    Fermi1.setElement(ONE_c,Indices(7,7,0));
    
    Fermi1.setElement(ONE_c,Indices(2,7,1));
    Fermi1.setElement(ONE_c,Indices(5,7,1));
    
    Fermi1.setElement(ONE_c,Indices(0,1,2));
    
    Fermi1.setElement(ONE_c,Indices(3,3,3));
    Fermi1.setElement(ONE_c,Indices(6,6,3));
    
    //Fermi2
    Fermi2.setElement(ONE_c,Indices(0,0,0));
    Fermi2.setElement(ONE_c,Indices(7,7,0));
    
    Fermi2.setElement(ONE_c,Indices(3,7,1));
    Fermi2.setElement(ONE_c,Indices(6,7,1));
    
    Fermi2.setElement(ONE_c,Indices(0,4,2));
    
    Fermi2.setElement(ONE_c,Indices(1,1,3));
    
    //Fermi2right
    Fermi2right.setElement(ONE_c,Indices(7,0,0));
      
    Fermi2right.setElement(ONE_c,Indices(3,0,1));
    Fermi2right.setElement(ONE_c,Indices(6,0,1));   
    
    //Gauge
    Gauge.setElement(ONE_c,Indices(0,0,0));
    Gauge.setElement(ONE_c,Indices(7,7,0));
    
    Gauge.setElement(ONE_c,Indices(1,2,1));
    
    Gauge.setElement(-1.0*I_c,Indices(1,3,2));
    
    Gauge.setElement(I_c,Indices(4,5,3));
    
    Gauge.setElement(ONE_c,Indices(4,6,4));
  }
  
  /*cout << "Fermi1: " << Fermi1 << endl;
  cout << "Fermi2: " << Fermi2 << endl;
  cout << "Gauge: " << Gauge << endl;
  cout << "GaugeFill: " << GaugeFiller << endl;
  cout << "Fermi1start: " << Fermi1start << endl;
  cout << "Fermi2end: " << Fermi2end << endl;
  cout << "FermiLeft: " << LeftFermi << endl;
  cout << "FermiRight: " << RightFermi << endl;
  cout << "GaugeLeft: " << LeftGauge << endl;
  cout << "GaugeRight: " << RightGauge << endl;*/
  
  //Now reshape and set the matrices
  Fermi1left.reshape(Indices(1*D_bond,4));
  Fermi1.reshape(Indices(D_bond*D_bond,4));
  Fermi2.reshape(Indices(D_bond*D_bond,4));
  Fermi2right.reshape(Indices(D_bond*1,4));
  
  Gauge.reshape(Indices(D_bond*D_bond,5)); 
  
  //Now contract the matrices
  mwArray M_Fermi1left,M_Fermi1,M_Fermi2,M_Fermi2right,M_Gauge;  
  
  M_Fermi1left=reshape(Fermi1left*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
  M_Fermi1=reshape(Fermi1*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
  M_Fermi2=reshape(Fermi2*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
  M_Fermi2right=reshape(Fermi2right*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
  M_Gauge=reshape(Gauge*Gauge_operators,Indices(D_bond,D_bond,this->d_link,this->d_link));
  
  //For testing:
  /*MPO test(5);
  test.setOp(0,new Operator(permute(M_Fermi1left,Indices(3,1,4,2))),true);
  test.setOp(1,new Operator(permute(M_Fermi2,Indices(3,1,4,2))),true);
  test.setOp(2,new Operator(permute(M_Gauge,Indices(3,1,4,2))),true);
  test.setOp(3,new Operator(permute(M_Fermi1,Indices(3,1,4,2))),true);
  test.setOp(4,new Operator(permute(M_Fermi2right,Indices(3,1,4,2))),true);
  
  test.exportForMatlab("SingleStringPos.m");
  test.clear();*/
  
  //Operator and flags, if the Operator was already saved
  bool even_site_saved=false,odd_site_saved=false,id_odd_saved,id_even_saved;
  int pos_even_saved,pos_odd_saved,pos_id_even_saved,pos_id_odd_saved;
  Operator *Local_op;
  
  int index;
  for(int k=startpos; k<endpos; k++)
  {
    //Transform k into a fermionic index for a fermion of species 1
    index=3*k-3;
    //Even site
    if((k%2)==0)
    {
      if(even_site_saved)
      {
	Peven.setOp(index,&Peven.getOp(pos_even_saved),false);
	Peven.setOp(index+1,&Peven.getOp(pos_even_saved+1),false);
	Peven.setOp(index+2,&Peven.getOp(pos_even_saved+2),false);
	Peven.setOp(index+3,&Peven.getOp(pos_even_saved+3),false);
	Peven.setOp(index+4,&Peven.getOp(pos_even_saved+4),false);
      }
      else
      {
	Local_op = new Operator(permute(M_Fermi1left,Indices(3,1,4,2)));
	Peven.setOp(index,Local_op,true);	
	Local_op = new Operator(permute(M_Fermi2,Indices(3,1,4,2)));
	Peven.setOp(index+1,Local_op,true);	
	Local_op = new Operator(permute(M_Gauge,Indices(3,1,4,2)));
	Peven.setOp(index+2,Local_op,true);	
	Local_op = new Operator(permute(M_Fermi1,Indices(3,1,4,2)));
	Peven.setOp(index+3,Local_op,true);	
	Local_op = new Operator(permute(M_Fermi2right,Indices(3,1,4,2)));
	Peven.setOp(index+4,Local_op,true);	
	even_site_saved=true;
	pos_even_saved=index;
      }
    }
    else
    {
      if(odd_site_saved)
      {
	Podd.setOp(index,&Podd.getOp(pos_odd_saved),false);
	Podd.setOp(index+1,&Podd.getOp(pos_odd_saved+1),false);
	Podd.setOp(index+2,&Podd.getOp(pos_odd_saved+2),false);
	Podd.setOp(index+3,&Podd.getOp(pos_odd_saved+3),false);
	Podd.setOp(index+4,&Podd.getOp(pos_odd_saved+4),false);
      }
      else
      {
	Local_op = new Operator(permute(M_Fermi1left,Indices(3,1,4,2)));
	Podd.setOp(index,Local_op,true);	
	Local_op = new Operator(permute(M_Fermi2,Indices(3,1,4,2)));
	Podd.setOp(index+1,Local_op,true);	
	Local_op = new Operator(permute(M_Gauge,Indices(3,1,4,2)));
	Podd.setOp(index+2,Local_op,true);	
	Local_op = new Operator(permute(M_Fermi1,Indices(3,1,4,2)));
	Podd.setOp(index+3,Local_op,true);	
	Local_op = new Operator(permute(M_Fermi2right,Indices(3,1,4,2)));
	Podd.setOp(index+4,Local_op,true);	
	odd_site_saved=true;
	pos_odd_saved=index;
      }
    }
  }    
  /*cout << "Podd neu: " << Podd << endl;
  cout << "Peven neu: " << Peven << endl;*/
  
  //Remaining identities for fermionic species 1 and 2
  id_odd_saved=false; 
  id_even_saved=false;
  for(int k=0; k<(this->N_total-1); k+=3)
  {
    if(Peven.isEmpty(k))
    {
      if(id_even_saved)
      {
	Peven.setOp(k,&Peven.getOp(pos_id_even_saved),false);
	Peven.setOp(k+1,&Peven.getOp(pos_id_even_saved+1),false);
      }
      else
      {
	Local_op = new Operator(this->id_fermi_mpo);
	Peven.setOp(k,Local_op,true);
	Peven.setOp(k+1,&Peven.getOp(k),false);
	id_even_saved=true;
	pos_id_even_saved=k;
      }
    }
    if(Podd.isEmpty(k))
    {
      if(id_odd_saved)
      {
	Podd.setOp(k,&Podd.getOp(pos_id_odd_saved),false);
	Podd.setOp(k+1,&Podd.getOp(pos_id_odd_saved+1),false);
      }
      else
      {
	Local_op = new Operator(this->id_fermi_mpo);
	Podd.setOp(k,Local_op,true);
	Podd.setOp(k+1,&Podd.getOp(k),false);
	id_odd_saved=true;
	pos_id_odd_saved=k;
      }
   }
  } 
  
  //Remaining identities for links
  id_odd_saved=false; 
  id_even_saved=false;
  for(int k=2; k<(this->N_total-2); k+=3)
  {
    if(Peven.isEmpty(k))
    {
      if(id_even_saved)
	Peven.setOp(k,&Peven.getOp(pos_id_even_saved),false);
      else
      {
	Local_op = new Operator(this->id_link_mpo);
	Peven.setOp(k,Local_op,true);
	id_even_saved=true;
	pos_id_even_saved=k;
      }
    }
    if(Podd.isEmpty(k))
    {
      if(id_odd_saved)
	Podd.setOp(k,&Podd.getOp(pos_id_odd_saved),false);
      else
      {
	Local_op = new Operator(this->id_link_mpo);
	Podd.setOp(k,Local_op,true);
	id_odd_saved=true;
	pos_id_odd_saved=k;
      }
    }
  }
  /*cout << "Podd neu: " << Podd << endl;
  cout << "Peven neu: " << Peven << endl;*/
}

void SU2Hamiltonian::getStringMPO2(MPO &mpo, unsigned int length, unsigned int pos,bool conjugate) const
{
  //Determine the starting position of the string
  if(pos<1)
    pos = round(this->N_fermions/2)-round(length/2);
  //Some error checking
  if(pos>this->N_total || pos<1 || (pos+length)>this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getStringMPO2, desired string of length " << length<< " starting on position " << pos << " is outside the the chain of " << this->N_fermions << " spins" << endl;
    exit(666);
  }
  if((length%2)==0)
  {
    cout << "Error in SU2Hamiltonian::getStringMPO2, only odd values for length are supported, received " << length << ", which is even" << endl;
    exit(666);
  }
  
  //cout << "Creating a string starting at position " << pos << " with length " << length << endl;
  
  //Set up the MPO accordingly
  mpo.initLength(this->N_total);
  
  //Determine whether my startposition is a even or odd fermionic index
  bool even_start_site;
  if(!conjugate)
  {
    //If I do not want the conjugate, I check if the startsite is even
    //cout << "Not taking the conjugate" << endl;
    even_start_site=((pos%2)==0);
  }
  else
  {
    //I simply pretend the startsite is opposite, therefore getting the conjugate operators
    //cout << "Taking the conjugate" << endl;
    even_start_site=!((pos%2)==0);    
  }
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(4,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(5,this->d_link,this->d_link));	Gauge_operators.fillWithZero();  
  
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
  
  for(int i=0; i<this->d_link; ++i)
  {
    for(int j=0; j<this->d_link; ++j)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      if(even_start_site)
      {
	//U00
	Gauge_operators.setElement(this->U00.getElement(Indices(i,j)),Indices(1,i,j));
	//U01
	Gauge_operators.setElement(this->U01.getElement(Indices(i,j)),Indices(2,i,j));
	//U10
	Gauge_operators.setElement(this->U10.getElement(Indices(i,j)),Indices(3,i,j));
	//U11
	Gauge_operators.setElement(this->U11.getElement(Indices(i,j)),Indices(4,i,j));
      }
      else
      {      
	//U00 dagger
	Gauge_operators.setElement(this->U00ad.getElement(Indices(i,j)),Indices(1,i,j));
	//U01 dagger
	Gauge_operators.setElement(this->U01ad.getElement(Indices(i,j)),Indices(2,i,j));
	//U10 dagger
	Gauge_operators.setElement(this->U10ad.getElement(Indices(i,j)),Indices(3,i,j));
	//U11 dagger
	Gauge_operators.setElement(this->U11ad.getElement(Indices(i,j)),Indices(4,i,j));     
      }
     }
  }   
  Gauge_operators.reshape(Indices(5,this->d_link*this->d_link)); 
  
  //Determine wether I start at an odd or even position, depending on that prepare the matrices  
  mwArray Fermi1left,Fermi1right,Fermi2left,Fermi2right,FermiId,GaugeFirst,Gauge;
  int D_bond=2;
  
  Fermi1left=mwArray(Indices(1,D_bond,4));		Fermi1left.fillWithZero();
  Fermi1right=mwArray(Indices(D_bond,D_bond,4));	Fermi1right.fillWithZero();
  Fermi2left=mwArray(Indices(D_bond,D_bond,4));		Fermi2left.fillWithZero();
  Fermi2right=mwArray(Indices(D_bond,1,4));		Fermi2right.fillWithZero();
  FermiId=mwArray(Indices(D_bond,D_bond,4));		FermiId.fillWithZero();
  GaugeFirst=mwArray(Indices(D_bond,D_bond,5));		GaugeFirst.fillWithZero();  
  Gauge=mwArray(Indices(D_bond,D_bond,5));		Gauge.fillWithZero();
  
  if(even_start_site)
  {
    Fermi1left.setElement(ONE_c,Indices(0,1,0));
    Fermi1left.setElement(ONE_c,Indices(0,0,1));
    
    Fermi2left.setElement(ONE_c,Indices(1,1,1));
    Fermi2left.setElement(ONE_c,Indices(0,0,3));
    
    Fermi1right.setElement(ONE_c,Indices(0,0,2));
    Fermi1right.setElement(ONE_c,Indices(1,1,3));
    
    Fermi2right.setElement(ONE_c,Indices(0,0,0));
    Fermi2right.setElement(ONE_c,Indices(1,0,2));
    
    GaugeFirst.setElement(ONE_c,Indices(0,0,1));
    GaugeFirst.setElement(I_c,Indices(0,1,2));
    GaugeFirst.setElement(-1.0*I_c,Indices(1,0,3));
    GaugeFirst.setElement(ONE_c,Indices(1,1,4));
  }
  else
  {
    Fermi1left.setElement(ONE_c,Indices(0,1,0));
    Fermi1left.setElement(ONE_c,Indices(0,0,2));
    
    Fermi2left.setElement(ONE_c,Indices(1,1,2));
    Fermi2left.setElement(ONE_c,Indices(0,0,3));
    
    Fermi1right.setElement(ONE_c,Indices(0,0,1));
    Fermi1right.setElement(ONE_c,Indices(1,1,3));
    
    Fermi2right.setElement(ONE_c,Indices(0,0,0));
    Fermi2right.setElement(ONE_c,Indices(1,0,1));
    
    GaugeFirst.setElement(ONE_c,Indices(0,0,1));
    GaugeFirst.setElement(-1.0*I_c,Indices(0,1,2));
    GaugeFirst.setElement(I_c,Indices(1,0,3));
    GaugeFirst.setElement(ONE_c,Indices(1,1,4));
  }
  FermiId.setElement(ONE_c,Indices(0,0,0));
  FermiId.setElement(ONE_c,Indices(1,1,0));
  
  //I have to take care of the factor i in front of the fermionic operators coming from the Jordan Wigner transformation, therefore I set the Gauge Matrix this way
  Gauge=GaugeFirst;
  
  Fermi1left.reshape(Indices(1*D_bond,4));
  Fermi1right.reshape(Indices(D_bond*D_bond,4));
  Fermi2left.reshape(Indices(D_bond*D_bond,4));
  Fermi2right.reshape(Indices(D_bond*1,4));
  FermiId.reshape(Indices(D_bond*D_bond,4));
  
  GaugeFirst.reshape(Indices(D_bond*D_bond,5)); 
  Gauge.reshape(Indices(D_bond*D_bond,5)); 
  
  //Now fill the MPO, first the starting fermionic sites and the ending fermionic sites
  mwArray res;
  Operator *Local_op;  
  
  //cout << "Setting start sites " << endl;
  
  res=reshape(Fermi1left*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
  Local_op = new Operator(permute(res,Indices(3,1,4,2)));
  mpo.setOp(3*pos-3,Local_op,true);
  res=reshape(Fermi2left*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
  Local_op = new Operator(permute(res,Indices(3,1,4,2)));
  mpo.setOp(3*pos-2,Local_op,true);
  
  res=reshape(GaugeFirst*Gauge_operators,Indices(D_bond,D_bond,this->d_link,this->d_link));
  Local_op = new Operator(permute(res,Indices(3,1,4,2)));
  mpo.setOp(3*pos-1,Local_op,true);
  
  res=reshape(Fermi1right*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
  Local_op = new Operator(permute(res,Indices(3,1,4,2)));
  mpo.setOp(3*(pos+length)-3,Local_op,true);
  res=reshape(Fermi2right*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
  Local_op = new Operator(permute(res,Indices(3,1,4,2)));
  mpo.setOp(3*(pos+length)-2,Local_op,true); 
  
  //cout << " ...done " << endl;
  
  //cout << "StringOp: " << mpo << endl;
  
    
  //No fill the sites inside the string 
  bool gauge_saved=false;
  bool fermi_saved=false;
  int ind_gauge,ind_fermi1, ind_fermi2,pos_gauge_saved,pos_fermi_saved;
  //cout << "Setting operator sites " << endl;
  for (int k=pos+1; k<pos+length;k++)
  {
    ind_fermi1=3*k-3;
    ind_fermi2=3*k-2;
    ind_gauge=3*k-1;
    
    /*cout << "Working on site " << k << endl;
    cout << "- Setting Op at site " << ind_fermi1 << endl;
    cout << "- Setting Op at site " << ind_fermi2 << endl;
    cout << "- Setting Op at site " << ind_gauge << endl;*/
    
    if(fermi_saved)
    {
      mpo.setOp(ind_fermi1,&mpo.getOp(pos_fermi_saved),false);
      mpo.setOp(ind_fermi2,&mpo.getOp(pos_fermi_saved),false);
    }
    else
    {
      res=reshape(FermiId*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(ind_fermi1,Local_op,true);
      fermi_saved=true;
      pos_fermi_saved=ind_fermi1;
      mpo.setOp(ind_fermi2,&mpo.getOp(pos_fermi_saved),false);
    }
       
    if(gauge_saved)
      mpo.setOp(ind_gauge,&mpo.getOp(pos_gauge_saved),false);
    else
    {
      res=reshape(Gauge*Gauge_operators,Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op = new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(ind_gauge,Local_op,true);
      gauge_saved=true;
      pos_gauge_saved=ind_gauge;
    }
  }
  //cout << " ...done" << endl;
  
  //cout << "StringOp: " << mpo << endl;
  
  //Remaining identities for fermionic species 1 and 2
  fermi_saved=false; 
  //cout << "Setting fermionic identities " << endl;
  for(int k=0; k<(this->N_total-1); k+=3)
  {
    if(mpo.isEmpty(k))
    {
      if(fermi_saved)
      {
	mpo.setOp(k,&mpo.getOp(pos_fermi_saved),false);
	mpo.setOp(k+1,&mpo.getOp(pos_fermi_saved),false);
      }
      else
      {
	Local_op = new Operator(this->id_fermi_mpo);
	mpo.setOp(k,Local_op,true);
	mpo.setOp(k+1,&mpo.getOp(k),false);
	fermi_saved=true;
	pos_fermi_saved=k;
      }
    }
  } 
  //cout << " ...done" << endl;
  
  //Remaining identities for links
  gauge_saved=false;
  //cout << "Setting bosonic identitites " << endl;
  for(int k=2; k<(this->N_total-2); k+=3)
  {
    if(mpo.isEmpty(k))
    {
      if(gauge_saved)
	mpo.setOp(k,&mpo.getOp(pos_gauge_saved),false);
      else
      {
	Local_op = new Operator(this->id_link_mpo);
	mpo.setOp(k,Local_op,true);
	gauge_saved=true;
	pos_gauge_saved=k;
      }
    }
  } 
  //cout << " ...done" << endl;
}

void SU2Hamiltonian::getChargeMPO(MPO &mpo, unsigned int pos, Component comp) const
{
  //Before doing anything, check if position is valid
  if(pos<1 || pos >this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getChargeMPO(), trying to retrieve charge operator for site "<< pos << " while total number of sites ranges form 1 to " << this->N_fermions << endl;
    exit(666);
  }
  //Initialize the MPO
  mpo.initLength(this->N_total);
  
  //Compute the two site inidces, where I have to put operators
  unsigned int fermi1_index=3*pos-3;
  unsigned int fermi2_index=3*pos-2;
  
  bool saved=false;
  int pos_saved;
  
  //Set the fermionic sites
  for (int i=0; i<this->N_total; i+=3)
  {
    if(i!=fermi1_index && !saved)
    {
      mpo.setOp(i,new Operator(this->id_fermi_mpo),true);
      saved=true;
      pos_saved=i;
      mpo.setOp(i+1,&mpo.getOp(i),false);
    }
    else if(i!=fermi1_index)
    {
      mpo.setOp(i,&mpo.getOp(pos_saved),false);
      mpo.setOp(i+1,&mpo.getOp(pos_saved),false);
    }
    
  }
  //Set the link sites
  saved=false;
  for (int i=2; i<this->N_total; i+=3)
  {
    if(!saved)
      mpo.setOp(i,new Operator(this->id_link_mpo),true);
    else
      mpo.setOp(i,&mpo.getOp(i-3),true);
  }
  
  //Prepare the matrices for the link
  int D_bond=2;
  //Array to store the fermionic operators
  mwArray Fermionic_operators(Indices(2,this->d_fermi,this->d_fermi));
  Fermionic_operators.fillWithZero();    
  if(comp==z_comp)
  {
    for(int i=0; i<this->d_fermi; ++i)
    {
      for(int j=0; j<this->d_fermi;++j)
      {
	//sigma_plus
	Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
	//sigma_minus
	Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(1,i,j));
      }
    }
  }
  else
  {
    for(int i=0; i<this->d_fermi; ++i)
    {
      for(int j=0; j<this->d_fermi;++j)
      {
	//sigma_plus
	Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(0,i,j));
	//sigma_minus
	Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(1,i,j));
      }
    }
  }
  Fermionic_operators.reshape(Indices(2,this->d_fermi*this->d_fermi));
  //cout << "Fermionic operators done" << endl;
  mwArray Fermi1, Fermi2;
  Fermi1=mwArray(Indices(1,D_bond,2));	Fermi1.fillWithZero();
  Fermi2=mwArray(Indices(D_bond,1,2));	Fermi2.fillWithZero();
  
  if(comp==x_comp)
  {
    //Case that I have the x-component
    Fermi1.setElement(0.5*ONE_c,Indices(0,0,0));
    Fermi1.setElement(0.5*ONE_c,Indices(0,1,1));
    
    Fermi2.setElement(I_c,Indices(1,0,0));
    Fermi2.setElement(-I_c,Indices(0,0,1));    
  }
  else if(comp==y_comp)
  {
    //Case that I have the y-component
    Fermi1.setElement(-0.5*ONE_c,Indices(0,0,0));
    Fermi1.setElement(-0.5*ONE_c,Indices(0,1,1));
    
    Fermi2.setElement(ONE_c,Indices(1,0,0));
    Fermi2.setElement(ONE_c,Indices(0,0,1));    
  }
  else
  {
    Fermi1.setElement(0.25*ONE_c,Indices(0,0,0));
    Fermi1.setElement(-0.25*ONE_c,Indices(0,1,1));
    
    Fermi2.setElement(ONE_c,Indices(1,0,0));
    Fermi2.setElement(ONE_c,Indices(0,0,1));    
  }
  
  Fermi1.reshape(Indices(1*D_bond,2));
  Fermi2.reshape(Indices(D_bond*1,2));
  
  //As I do not need the matrices Fermi1 and Fermi2 any further, I contract them with the operators and reshape them appropriately
  Fermi1=reshape(Fermi1*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
  Fermi2=reshape(Fermi2*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
  
  //Plug the operators in the MPO
  mpo.setOp(fermi1_index,new Operator(permute(Fermi1,Indices(3,1,4,2))),true);
  mpo.setOp(fermi2_index,new Operator(permute(Fermi2,Indices(3,1,4,2))),true);
  
  //cout << "Charge MPO: " << mpo << endl;
  
  /*char buffer [50];
  if(comp==x_comp)
    sprintf (buffer, "ChargeOp_xcomp_n%d.m",pos);
  else if(comp==y_comp)
    sprintf (buffer, "ChargeOp_ycomp_n%d.m",pos);
  else
    sprintf (buffer, "ChargeOp_zcomp_n%d.m",pos);
  mpo.exportForMatlab(buffer);*/
}

void SU2Hamiltonian::getChargeSquareMPO(MPO &mpo) const
{
  mpo.initLength(this->N_total);
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(2,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(1,this->d_link,this->d_link));	Gauge_operators.fillWithZero();  
  
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
  
  for(int i=0; i<this->d_link; ++i)
    for(int j=0; j<this->d_link; ++j)
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
    
  Gauge_operators.reshape(Indices(1,this->d_link*this->d_link)); 
  
  
  //Prepare the matrices  
  mwArray Fermi1,Fermi2,Fermifirst,Fermilast,Gauge;
  int D_bond=3;
  
  Fermifirst=mwArray(Indices(1,D_bond,2));	Fermifirst.fillWithZero();
  Fermi1=mwArray(Indices(D_bond,D_bond,2));	Fermi1.fillWithZero();
  Fermi2=mwArray(Indices(D_bond,D_bond,2));	Fermi2.fillWithZero();
  Fermilast=mwArray(Indices(D_bond,1,2));	Fermilast.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,1));	Gauge.fillWithZero();
  
  //Identity
  Fermifirst.setElement(ONE_c,Indices(0,0,0));
  Fermifirst.setElement(1.0/8.0*ONE_c,Indices(0,2,0));
  
  Fermi1.setElement(ONE_c,Indices(0,0,0));
  Fermi1.setElement(1.0/8.0*ONE_c,Indices(0,2,0));
  Fermi1.setElement(ONE_c,Indices(2,2,0));
  
  Fermi2.setElement(ONE_c,Indices(0,0,0));
  Fermi2.setElement(ONE_c,Indices(2,2,0));
  
  Fermilast.setElement(ONE_c,Indices(2,0,0));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(2,2,0));
    
  //sigma_z
  Fermifirst.setElement(-1.0/8.0*ONE_c,Indices(0,1,1));
    
  Fermi1.setElement(-1.0/8.0*ONE_c,Indices(0,1,1));
  
  Fermi2.setElement(ONE_c,Indices(1,2,1));
  
  Fermilast.setElement(ONE_c,Indices(1,0,1));
  
  //Reshape the matrices accordingly
  Fermifirst.reshape(Indices(1*D_bond,2));
  Fermi1.reshape(Indices(D_bond*D_bond,2));
  Fermi2.reshape(Indices(D_bond*D_bond,2));
  Fermilast.reshape(Indices(D_bond*1,2));
  Gauge.reshape(Indices(D_bond*D_bond,1));

  mwArray res;
  Operator* Local_op;
  
  bool fermi1_saved=false,fermi2_saved=false,gauge_saved=false;
  int pos_fermi1,pos_fermi2;
  
  //Set the first and the last operator
  res=reshape(Fermifirst*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
  res=permute(res,Indices(3,1,4,2));
  Local_op=new Operator(res);
  mpo.setOp(0,Local_op,true);
  
  res=reshape(Fermilast*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
  res=permute(res,Indices(3,1,4,2));
  Local_op=new Operator(res);
  mpo.setOp(this->N_total-1,Local_op,true);
  
  
  //Rest of the fermionic sites of species 1
  for(int i=3;i<this->N_total; i+=3)
  {
    if(fermi1_saved)
      mpo.setOp(i,&mpo.getOp(pos_fermi1),false);
    else
    {
      res=reshape(Fermi1*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      res=permute(res,Indices(3,1,4,2));
      Local_op=new Operator(res);
      mpo.setOp(i,Local_op,true);
      fermi1_saved=true;
      pos_fermi1=i;
    }
  }
  
  //Rest of the fermionic sites of species 1
  for(int i=1;i<this->N_total-1; i+=3)
  {
    if(fermi2_saved)
      mpo.setOp(i,&mpo.getOp(pos_fermi2),false);
    else
    {
      res=reshape(Fermi2*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      res=permute(res,Indices(3,1,4,2));
      Local_op=new Operator(res);
      mpo.setOp(i,Local_op,true);
      fermi2_saved=true;
      pos_fermi2=i;
    }
  }
  
  //Rest of the fermionic sites of species 1
  for(int i=2;i<this->N_total; i+=3)
  {
    if(gauge_saved)
      mpo.setOp(i,&mpo.getOp(i-3),false);
    else
    {
      res=reshape(Gauge*Gauge_operators,Indices(D_bond,D_bond,this->d_link,this->d_link));
      res=permute(res,Indices(3,1,4,2));
      Local_op=new Operator(res);
      mpo.setOp(i,Local_op,true);
      gauge_saved=true;
    }
  }
}
  
  
void SU2Hamiltonian::getChargeSquareMPO(MPO &mpo, unsigned int pos) const
{
  //Before doing anything, check if position is valid
  if(pos<1 || pos >this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getChargeSquareMPO(), trying to retrieve charge operator for site "<< pos << " while total number of sites ranges form 1 to " << this->N_fermions << endl;
    exit(666);
  }
  mpo.initLength(this->N_total);
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(2,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  
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
  
  //Prepare the matrices  
  mwArray Fermi1,Fermi2;
  int D_bond=2;
  
  Fermi1=mwArray(Indices(1,D_bond,2));	Fermi1.fillWithZero();
  Fermi2=mwArray(Indices(D_bond,1,2));	Fermi2.fillWithZero();
  
  //Identity
  Fermi1.setElement(1.0/8.0*ONE_c,Indices(0,0,0));  
  Fermi2.setElement(ONE_c,Indices(0,0,0));
  
  //sigma_z    
  Fermi1.setElement(-1.0/8.0*ONE_c,Indices(0,1,1));  
  Fermi2.setElement(ONE_c,Indices(1,0,1));
    
  //Reshape the matrices accordingly
  Fermi1.reshape(Indices(1*D_bond,2));
  Fermi2.reshape(Indices(D_bond*1,2));

  mwArray res;
  Operator* Local_op;
  
  bool saved=false;
  int pos_saved;
  
  //Set the operators
  res=reshape(Fermi1*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
  res=permute(res,Indices(3,1,4,2));
  Local_op=new Operator(res);
  mpo.setOp(3*pos-3,Local_op,true);
  
  res=reshape(Fermi2*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
  res=permute(res,Indices(3,1,4,2));
  Local_op=new Operator(res);
  mpo.setOp(3*pos-2,Local_op,true);
  
  
  //Rest of the fermionic sites
  for(int i=0;i<this->N_total; i+=3)
  {
    if(saved && (3*pos-3)!=i)
    {
      mpo.setOp(i,&mpo.getOp(pos_saved),false);
      mpo.setOp(i+1,&mpo.getOp(pos_saved),false);
    }
    else if((3*pos-3)!=i)
    {
      mpo.setOp(i,new Operator(this->id_fermi_mpo),true);
      saved=true;
      pos_saved=i;
      mpo.setOp(i+1,&mpo.getOp(pos_saved),false);
    }
  }
  
  //Gauge variables
  saved=false;
  for(int i=2;i<this->N_total; i+=3)
  {
    if(saved)
      mpo.setOp(i,&mpo.getOp(i-3),false);
    else
    {
      mpo.setOp(i,new Operator(this->id_link_mpo),true);
      saved=true;
    }
  }
}

void SU2Hamiltonian::getFermionNumberMPO(MPO &mpo) const
{
  mpo.initLength(this->N_total);
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(2,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(1,this->d_link,this->d_link));	Gauge_operators.fillWithZero();  
  
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
  
  for(int i=0; i<this->d_link; ++i)
    for(int j=0; j<this->d_link; ++j)
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
    
  Gauge_operators.reshape(Indices(1,this->d_link*this->d_link)); 
  
  
  //Prepare the matrices  
  mwArray Fermi1,Fermi2,Fermifirst,Fermilast,Gauge;
  int D_bond=2;
  
  Fermifirst=mwArray(Indices(1,D_bond,2));	Fermifirst.fillWithZero();
  Fermi1=mwArray(Indices(D_bond,D_bond,2));	Fermi1.fillWithZero();
  Fermi2=mwArray(Indices(D_bond,D_bond,2));	Fermi2.fillWithZero();
  Fermilast=mwArray(Indices(D_bond,1,2));	Fermilast.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,1));	Gauge.fillWithZero();
  
  //Identity
  Fermifirst.setElement(ONE_c,Indices(0,0,0));
  Fermifirst.setElement(ONE_c,Indices(0,1,0));
  
  Fermi1.setElement(ONE_c,Indices(0,0,0));
  Fermi1.setElement(ONE_c,Indices(0,1,0));
  Fermi1.setElement(ONE_c,Indices(1,1,0));
  
  Fermi2.setElement(ONE_c,Indices(0,0,0));
  Fermi2.setElement(ONE_c,Indices(1,1,0));
  
  Fermilast.setElement(ONE_c,Indices(1,0,0));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(1,1,0));
    
  //sigma_z
  Fermifirst.setElement(1.0/2.0*ONE_c,Indices(0,1,1));    
  Fermi1.setElement(1.0/2.0*ONE_c,Indices(0,1,1));  
  Fermi2.setElement(1.0/2.0*ONE_c,Indices(0,1,1));  
  Fermilast.setElement(1.0/2.0*ONE_c,Indices(0,0,1));
  
  //Reshape the matrices accordingly
  Fermifirst.reshape(Indices(1*D_bond,2));
  Fermi1.reshape(Indices(D_bond*D_bond,2));
  Fermi2.reshape(Indices(D_bond*D_bond,2));
  Fermilast.reshape(Indices(D_bond*1,2));
  Gauge.reshape(Indices(D_bond*D_bond,1));

  mwArray res;
  Operator* Local_op;
  
  bool fermi1_saved=false,fermi2_saved=false,gauge_saved=false;
  int pos_fermi1,pos_fermi2;
  
  //Set the first and the last operator
  res=reshape(Fermifirst*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
  res=permute(res,Indices(3,1,4,2));
  Local_op=new Operator(res);
  mpo.setOp(0,Local_op,true);
  
  res=reshape(Fermilast*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
  res=permute(res,Indices(3,1,4,2));
  Local_op=new Operator(res);
  mpo.setOp(this->N_total-1,Local_op,true);
  
  
  //Rest of the fermionic sites of species 1
  for(int i=3;i<this->N_total; i+=3)
  {
    if(fermi1_saved)
      mpo.setOp(i,&mpo.getOp(pos_fermi1),false);
    else
    {
      res=reshape(Fermi1*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      res=permute(res,Indices(3,1,4,2));
      Local_op=new Operator(res);
      mpo.setOp(i,Local_op,true);
      fermi1_saved=true;
      pos_fermi1=i;
    }
  }
  
  //Rest of the fermionic sites of species 2
  for(int i=1;i<this->N_total-1; i+=3)
  {
    if(fermi2_saved)
      mpo.setOp(i,&mpo.getOp(pos_fermi2),false);
    else
    {
      res=reshape(Fermi2*Fermionic_operators,Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      res=permute(res,Indices(3,1,4,2));
      Local_op=new Operator(res);
      mpo.setOp(i,Local_op,true);
      fermi2_saved=true;
      pos_fermi2=i;
    }
  }
  
  //Gauge links
  for(int i=2;i<this->N_total; i+=3)
  {
    if(gauge_saved)
      mpo.setOp(i,&mpo.getOp(i-3),false);
    else
    {
      res=reshape(Gauge*Gauge_operators,Indices(D_bond,D_bond,this->d_link,this->d_link));
      res=permute(res,Indices(3,1,4,2));
      Local_op=new Operator(res);
      mpo.setOp(i,Local_op,true);
      gauge_saved=true;
    }
  }
}
  
void SU2Hamiltonian::getFermionNumberMPO(MPO &mpo, unsigned int pos) const
{
  //Before doing anything, check if position is valid
  if(pos<1 || pos >this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getChargeSquareMPO(), trying to retrieve charge operator for site "<< pos << " while total number of sites ranges form 1 to " << this->N_fermions << endl;
    exit(666);
  }
  mpo.initLength(this->N_total);
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(2,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  
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
  
  //Prepare the matrices  
  mwArray Fermi1,Fermi2;
  int D_bond=2;
  
  Fermi1=mwArray(Indices(1,D_bond,2));	Fermi1.fillWithZero();
  Fermi2=mwArray(Indices(D_bond,1,2));	Fermi2.fillWithZero();
  
  //Identity  
  Fermi1.setElement(ONE_c,Indices(0,0,0));
  Fermi1.setElement(ONE_c,Indices(0,1,0));
  
  Fermi2.setElement(ONE_c,Indices(1,0,0));
 
  //sigma_z
  Fermi1.setElement(1.0/2.0*ONE_c,Indices(0,1,1));  
  Fermi2.setElement(1.0/2.0*ONE_c,Indices(0,0,1));  
  
  //Reshape the matrices accordingly
  Fermi1.reshape(Indices(1*D_bond,2));
  Fermi2.reshape(Indices(D_bond*1,2));

  bool saved=false;
  int pos_saved;
  
  //Set the operators
  Fermi1=reshape(Fermi1*Fermionic_operators,Indices(1,D_bond,this->d_fermi,this->d_fermi));
  Fermi1=permute(Fermi1,Indices(3,1,4,2));
  mpo.setOp(3*pos-3,new Operator(Fermi1),true);
  
  Fermi2=reshape(Fermi2*Fermionic_operators,Indices(D_bond,1,this->d_fermi,this->d_fermi));
  Fermi2=permute(Fermi2,Indices(3,1,4,2));
  mpo.setOp(3*pos-2,new Operator(Fermi2),true);
  
  
  //Rest of the fermionic sites
  for(int i=0;i<this->N_total; i+=3)
  {
    if(saved && i!=(3*pos-3))
    {
      mpo.setOp(i,&mpo.getOp(pos_saved),false);
      mpo.setOp(i+1,&mpo.getOp(pos_saved),false);
    }
    else if(i!=(3*pos-3))
    {
      mpo.setOp(i,new Operator(this->id_fermi_mpo),true);
      saved=true;
      pos_saved=i;
      mpo.setOp(i+1,&mpo.getOp(pos_saved),false);
    }
  }
  
  //The gauge links
  saved=false;
  for(int i=2;i<this->N_total; i+=3)
  {
    if(saved)
      mpo.setOp(i,&mpo.getOp(i-3),false);
    else
    {
      mpo.setOp(i,new Operator(this->id_link_mpo),true);
      saved=true;
    }
  }
}


//Build the reduced MPO to detect strings
void SU2Hamiltonian::buildReducedStringMPO(MPO &mpo,const MPS &mps,int index_left_site,int index_right_site) const
{
  
  //Make sure the MPO has the right length
  mpo.initLength(this->N_total);
  
  //Get a contractor, for testing
  //Contractor& contractor=Contractor::theContractor();
  
  //I rely on the fact, that the contractions to the left and those to the right are simply identites in case the MPS is properly gauged. As I do not want to destroy the input, I get a copy of it and gauge it accordingly
  MPS local_mps(mps);
  
  local_mps.gaugeCond('R',true);
  
  for(int i=(this->N_total-1); i>(3*index_right_site-2); i--)
    local_mps.gaugeCond(i,'L',true);
  
  //Now set the identities
  bool fermi_id_saved=false,link_id_saved=false;
  int pos_fermi_id_saved,pos_link_id_saved;
  for(int i=1; i<=this->N_fermions; i++)
  {
    if(fermi_id_saved && (i<index_left_site || i>index_right_site))
    {
      mpo.setOp(3*i-3,&mpo.getOp(pos_fermi_id_saved),false);
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }
    else if(i<index_left_site || i>index_right_site)
    {
      mpo.setOp(3*i-3,new Operator(this->id_fermi_mpo),true);
      fermi_id_saved=true;
      pos_fermi_id_saved=3*i-3;
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }
    if(link_id_saved && (i<index_left_site || i>=index_right_site) && i<this->N_fermions)
      mpo.setOp(3*i-1,&mpo.getOp(pos_link_id_saved),false);
    else if((i<index_left_site || i>=index_right_site) && i<this->N_fermions)
    {
      mpo.setOp(3*i-1,new Operator(this->id_link_mpo),true);
      link_id_saved=true;
      pos_link_id_saved=3*i-1;
    }
  }
  //cout << "MPO after identites: " << mpo << endl;
  
  //Now proceed with computing the left and right edge tensor
  mwArray BraMatrix,KetMatrix,LeftTensor,RightTensor;
  int Dl,Dr,d;
  //LeftTensor
  Dl=local_mps.getA(3*index_left_site-3).getDl();
  Dr=local_mps.getA(3*index_left_site-3).getDr();
  d=local_mps.getA(3*index_left_site-3).getd();
  KetMatrix=local_mps.getA(3*index_left_site-3).getA();
  BraMatrix=KetMatrix; BraMatrix.conjugate();
  
  KetMatrix.permute(Indices(2,1,3));
  KetMatrix.reshape(Indices(Dl,d*Dr));
  BraMatrix.permute(Indices(1,3,2));
  BraMatrix.reshape(Indices(d*Dr,Dl));
  LeftTensor=BraMatrix*KetMatrix;
  LeftTensor.reshape(Indices(d,Dr,d,Dr));
  LeftTensor.permute(Indices(3,1,4,2));
  LeftTensor.reshape(Indices(d,1,d,Dr*Dr));
  
  //RightTensor
  Dl=local_mps.getA(3*index_right_site-2).getDl();
  Dr=local_mps.getA(3*index_right_site-2).getDr();
  d=local_mps.getA(3*index_right_site-2).getd();
  KetMatrix=local_mps.getA(3*index_right_site-2).getA();
  BraMatrix=KetMatrix; BraMatrix.conjugate();
  
  KetMatrix.reshape(Indices(d*Dl,Dr));
  BraMatrix.permute(Indices(3,1,2));
  BraMatrix.reshape(Indices(Dr,d*Dl));
  RightTensor=KetMatrix*BraMatrix;
  RightTensor.reshape(Indices(d,Dl,d,Dl));
  RightTensor.permute(Indices(1,2,4,3));
  RightTensor.reshape(Indices(d,Dl*Dl,d,1));
  
  //Set the edge tensors
  mpo.setOp(3*index_left_site-3,new Operator(LeftTensor),true);
  mpo.setOp(3*index_right_site-2,new Operator(RightTensor),true);
  
  //cout << "MPO after setting edge Tensors: " << mpo << endl;
  
  //Build a testmpo containing only the interesting part in order to be able to test this thing
  /*MPO testmpo((3*index_right_site-2)-(3*index_left_site-3)+1);  
  testmpo.setOp(0,new Operator(LeftTensor), true);
  testmpo.setOp((3*index_right_site-2)-(3*index_left_site-3),new Operator(RightTensor), true);*/
  
  //Now get the rest of the tensors in between the edge tensors
  mwArray localMatrix,mpsMatrix;
  int count=1;
  for (int i=index_left_site; i<index_right_site; i++)
  {
    //Site of fermi species 2
    Dl=local_mps.getA(3*i-2).getDl();
    Dr=local_mps.getA(3*i-2).getDr();
    d=local_mps.getA(3*i-2).getd();
    mpsMatrix=local_mps.getA(3*i-2).getA();
    localMatrix=mwArray(Indices(d,Dl,Dl,d,Dr,Dr));
    for(int Dlk=0; Dlk<Dl; Dlk++)
      for(int Drk=0; Drk<Dr; Drk++)
	for(int dk=0; dk<d; dk++)
	  for(int Dlb=0; Dlb<Dl; Dlb++)
	    for(int Drb=0; Drb<Dr; Drb++)
	      for(int db=0; db<d; db++)
		localMatrix.setElement(mpsMatrix.getElement(Indices(dk,Dlk,Drk))*conjugate(mpsMatrix.getElement(Indices(db,Dlb,Drb))),Indices(dk,Dlk,Dlb,db,Drk,Drb));
    
    localMatrix.reshape(Indices(d,Dl*Dl,d,Dr*Dr));
    //cout << "Setting site " << 3*i-2 << endl;
    mpo.setOp(3*i-2,new Operator(localMatrix),true);
    //testmpo.setOp(count,new Operator(localMatrix),true);
    //count++;
    
    //Link
    Dl=local_mps.getA(3*i-1).getDl();
    Dr=local_mps.getA(3*i-1).getDr();
    d=local_mps.getA(3*i-1).getd();
    mpsMatrix=local_mps.getA(3*i-1).getA();
    localMatrix=mwArray(Indices(d,Dl,Dl,d,Dr,Dr));    
    for(int Dlk=0; Dlk<Dl; Dlk++)
      for(int Drk=0; Drk<Dr; Drk++)
	for(int dk=0; dk<d; dk++)
	  for(int Dlb=0; Dlb<Dl; Dlb++)
	    for(int Drb=0; Drb<Dr; Drb++)
	      for(int db=0; db<d; db++)
		localMatrix.setElement(mpsMatrix.getElement(Indices(dk,Dlk,Drk))*conjugate(mpsMatrix.getElement(Indices(db,Dlb,Drb))),Indices(dk,Dlk,Dlb,db,Drk,Drb));
    
    localMatrix.reshape(Indices(d,Dl*Dl,d,Dr*Dr));
    //cout << "Setting site " << 3*i-1 << endl;
    mpo.setOp(3*i-1,new Operator(localMatrix),true);
    //testmpo.setOp(count,new Operator(localMatrix),true);
    //count++;
    
    //Next site with fermi species 1
    Dl=local_mps.getA(3*(i+1)-3).getDl();
    Dr=local_mps.getA(3*(i+1)-3).getDr();
    d=local_mps.getA(3*(i+1)-3).getd();
    mpsMatrix=local_mps.getA(3*(i+1)-3).getA();
    localMatrix=mwArray(Indices(d,Dl,Dl,d,Dr,Dr));
    for(int Dlk=0; Dlk<Dl; Dlk++)
      for(int Drk=0; Drk<Dr; Drk++)
	for(int dk=0; dk<d; dk++)
	  for(int Dlb=0; Dlb<Dl; Dlb++)
	    for(int Drb=0; Drb<Dr; Drb++)
	      for(int db=0; db<d; db++)
		localMatrix.setElement(mpsMatrix.getElement(Indices(dk,Dlk,Drk))*conjugate(mpsMatrix.getElement(Indices(db,Dlb,Drb))),Indices(dk,Dlk,Dlb,db,Drk,Drb));
    
    localMatrix.reshape(Indices(d,Dl*Dl,d,Dr*Dr));
    //cout << "Setting site " << 3*(i+1)-3 << endl;
    mpo.setOp(3*(i+1)-3,new Operator(localMatrix),true);
    //testmpo.setOp(count,new Operator(localMatrix),true);
    //count++;
  } 
  /*char name[30];
  sprintf(name,"Trace_start%i_end%i_leng%i.m",index_left_site,index_right_site,(3*index_right_site-2)-(3*index_left_site-3)+1);
  testmpo.exportForMatlab(name);*/
}

//Needed for searching, sorting and finding elements in vectors (for debugging)
#include <algorithm>
void SU2Hamiltonian::getStringProjector(MPO &mpo,int length,int pos, unsigned int fermions) const
{  
  //Some error checking
  if(pos>this->N_total || pos<1 || (pos+length)>this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getStringProjector(), desired string of length " << length<< " starting on position " << pos << " is outside the the chain of " << this->N_fermions << " spins" << endl;
    exit(666);
  }
  if(fermions > 2)
  {
    cout << "Warning in SU2Hamiltonian::getStringProjector(), desired mode for the fermions is not supported, negelecting the fermionic configuration" << endl;
    fermions=0;
  }
  
  //cout << "Creating string with length " << length << " at position " << pos << endl;
  
  vector <int> setting_pos;    
  int item;  
  
  
  //Initialize the MP0
  mpo.initLength(this->N_total);
  
  //Fermionic identities outside the string and first site inside the string
  bool fermi_id_saved=false,link_id_saved=false;
  int pos_fermi_id_saved=-1,pos_link_id_saved=-1;
  //To the left
  //cout << "- Left identities" << endl;
  for(int i=1; i<=pos; i++)
  {
    if(!fermi_id_saved)
    {
      //I have to save the operator
      //cout << "Setting fermi_id at " << i << " and saving it" << endl;
      //cout << "Setting operator at " << 3*i-3 << endl;
      item = 3*i-3;
      if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
      setting_pos.push_back(3*i-3);
      mpo.setOp(3*i-3,new Operator(this->id_fermi_mpo),true);
      fermi_id_saved=true;
      pos_fermi_id_saved=3*i-3;
      //cout << "Setting pointer at " << 3*i-2 << endl;
      item = 3*i-2;
      if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
      setting_pos.push_back(3*i-2);
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }
    else
    {
      //Just get the pointer
      //cout << "Setting fermi_id at " << i << endl;
      //cout << "Setting pointer at " << 3*i-3 << endl;
      item = 3*i-3;
      if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
      setting_pos.push_back(3*i-3);
      mpo.setOp(3*i-3,&mpo.getOp(pos_fermi_id_saved),false);
      //cout << "Setting pointer at " << 3*i-2 << endl;
      item = 3*i-2;
      if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
      setting_pos.push_back(3*i-2);
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }
    if(!link_id_saved && i<(pos-1))
    {
      //I have to save the operator
      //cout << "Setting link_id at " << i << " and saving it" << endl;
      //cout << "Setting operator at " << 3*i-1 << endl;
      item = 3*i-1;
      if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
      setting_pos.push_back(3*i-1);
      mpo.setOp(3*i-1,new Operator(this->id_link_mpo),true);
      link_id_saved=true;
      pos_link_id_saved=3*i-1;
    }
    else if(i<(pos-1))
    {
      //cout << "Setting link_id at " << i << endl;
      //cout << "Setting pointer at " << 3*i-1 << endl;
      item = 3*i-1;
      if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
      setting_pos.push_back(3*i-1);
      mpo.setOp(3*i-1,&mpo.getOp(pos_link_id_saved),false);
    }
  }
  //To the right
  //cout << "- Right identities" << endl;
  for(int i=(pos+length); i<=this->N_fermions; i++)
  {
    if(!fermi_id_saved)
    {
      //I have to save the operator
      //cout << "Setting fermi_id at " << i << " and saving it" <<endl;
      //cout << "Setting operator at " << 3*i-3 << endl;
      item = 3*i-3;
      if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
      setting_pos.push_back(3*i-3);
      mpo.setOp(3*i-3,new Operator(this->id_fermi_mpo),true);
      fermi_id_saved=true;
      pos_fermi_id_saved=3*i-3;
      //cout << "Setting pointer at " << 3*i-2 << endl;
      item = 3*i-2;
      if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
      setting_pos.push_back(3*i-2);
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }
    else
    {
      //Just get the pointer
      //cout << "Setting fermi_id at " << i << endl;
      //cout << "Setting pointer at " << 3*i-3 << endl;
      item = 3*i-3;
      if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
      setting_pos.push_back(3*i-3);
      mpo.setOp(3*i-3,&mpo.getOp(pos_fermi_id_saved),false);
      //cout << "Setting pointer at " << 3*i-2 << endl;
      item = 3*i-2;
      if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
      setting_pos.push_back(3*i-2);
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }
    if(!link_id_saved && i>(pos+length) && i<this->N_fermions)
    {
      //I have to save the operator
      //cout << "Setting link_id at " << i << " and saving it" << endl;
      //cout << "Setting operator at " << 3*i-1 << endl;
      item = 3*i-1;
      if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
      setting_pos.push_back(3*i-1);
      mpo.setOp(3*i-1,new Operator(this->id_link_mpo),true);
      link_id_saved=true;
      pos_link_id_saved=3*i-1;
    }
    else if(i>(pos+length) && i<this->N_fermions)
    {
      //cout << "Setting link_id at " << i << endl;
      //cout << "Setting pointer at " << 3*i-1 << endl;
      item = 3*i-1;
      if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
      setting_pos.push_back(3*i-1);
      mpo.setOp(3*i-1,&mpo.getOp(pos_link_id_saved),false);
    }
      
  } 
  
  //Matrix sitting at the ends of the string detecting whether there is a zero flux state
  mwArray endmatrix(Indices(this->d_link,1,this->d_link,1));
  endmatrix.fillWithZero();
  endmatrix.setElement(ONE_c,Indices(2,0,2,0));
  //Matrix for the region in between detecting only non-zero flux states
  mwArray innermatrix(Indices(this->d_link,1,this->d_link,1));
  innermatrix.fillWithZero();
  for(int i=0; i<this->d_link; i++)
  {
    if(i!=2)
      innermatrix.setElement(ONE_c,Indices(i,0,i,0));
  }
  /*cout << "Endmatrix: " << endmatrix << endl;
  cout << "Innermatrix: " << innermatrix << endl;*/
    
  //cout << "MPO: " << mpo << endl;
  //cout.flush();
  
  //Set initial matrix left of the start position (in case this is not already the end of the chain)
  //cout << "- Initial and end Matrix" << endl;
  if((3*(pos-1)-1)>0)
  {
    //cout << "Setting operator at " << 3*(pos-1)-1 << endl;
    item = 3*(pos-1)-1;
    if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
    setting_pos.push_back(3*(pos-1)-1);
    mpo.setOp(3*(pos-1)-1,new Operator(endmatrix),true);
  }
  if((3*(pos+length)-1)<this->N_total)
  {
    //cout << "Setting operator at " << 3*(pos+length)-1 << endl;
    item = 3*(pos+length)-1;
    if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
    setting_pos.push_back(3*(pos+length)-1);
    mpo.setOp(3*(pos+length)-1,new Operator(endmatrix),true);
  }
  //cout << "Matrices in the string" << endl;
  //Set the matrices in between
  for(int i=pos; i<(pos+length); i++)
  {
    //cout << "Setting operator at " << 3*i-1 << endl;
    item = 3*i-1;
    if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	cout << "----> Error" << endl;
    setting_pos.push_back(3*i-1);
    mpo.setOp(3*i-1,new Operator(innermatrix),true);
  }

  //Now fill the rest of the fermionic sites accordingly
  //cout << "- Fermionic sites in the string" << endl;
  if(fermions==0)
  {
    //In this case I do not care about the fermionic configuration, just put identities
    for(int i=(pos+1); i<(pos+length); i++)
    {
      if(fermi_id_saved) 
      {
	//I have to save the operator
	//cout << "Setting fermi_id at " << i << " and saving it" << endl;
	//cout << "Setting operator at " << 3*i-3 << endl;
	item = 3*i-3;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-3);
	mpo.setOp(3*i-3,new Operator(this->id_fermi_mpo),true);
	fermi_id_saved=true;
	pos_fermi_id_saved=3*i-3;
	//cout << "Setting pointer at " << 3*i-2 << endl;
	item = 3*i-2;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-2);
	mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
      }
      else
      {
	//Just get the pointer
	//cout << "Setting fermi_id at " << i << endl;
	//cout << "Setting pointer at " << 3*i-3 << endl;
	item = 3*i-3;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-3);
	mpo.setOp(3*i-3,&mpo.getOp(pos_fermi_id_saved),false);
	//cout << "Setting pointer at " << 3*i-2 << endl;
	item = 3*i-2;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-2);
	mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
      }
    }
  }
  else if(fermions==1)
  {
    //Project out antiparallel configurations
    mwArray fermi1(Indices(this->d_fermi,1,this->d_fermi,2));
    mwArray fermi2(Indices(this->d_fermi,2,this->d_fermi,1));
    fermi1.fillWithZero();	fermi1.setElement(ONE_c,Indices(0,0,0,0));	fermi1.setElement(ONE_c,Indices(1,0,1,1));
    fermi2.fillWithZero();	fermi2.setElement(ONE_c,Indices(0,0,0,0));	fermi2.setElement(ONE_c,Indices(1,1,1,0));
    bool fermi1_saved=false,fermi2_saved=false;
    int pos_fermi1_saved=-1, pos_fermi2_saved=-1;
    for(int i=(pos+1); i<(pos+length); i++)
    {
      if(fermi1_saved) 
      {
	//cout << "Setting pointer at " << 3*i-3 << endl;
	item = 3*i-3;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-3);
	mpo.setOp(3*i-3,&mpo.getOp(pos_fermi1_saved),false);
      }
      else
      {
	//cout << "Setting operator at " << 3*i-3 << endl;
	item=3*i-3;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-3);
	mpo.setOp(3*i-3,new Operator(fermi1),true);
	fermi1_saved=true;
	pos_fermi1_saved=3*i-3;
      }
	
      if(fermi2_saved) 
      {
	//cout << "Setting pointer at " << 3*i-2 << endl;
	item = 3*i-2;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-2);
	mpo.setOp(3*i-2,&mpo.getOp(pos_fermi2_saved),false);
      }
      else
      {
	//cout << "Setting operator at " << 3*i-2 << endl;
	item = 3*i-2;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-2);
	mpo.setOp(3*i-2,new Operator(fermi2),true);
	fermi2_saved=true;
	pos_fermi2_saved=3*i-2;
      }
    }
  }
  else
  {
    //Project out antiparallel and parallel configurations with the wrong polarisation
    mwArray fermiodd(Indices(this->d_fermi,1,this->d_fermi,1));
    mwArray fermieven(Indices(this->d_fermi,1,this->d_fermi,1));
    fermiodd.fillWithZero();	fermiodd.setElement(ONE_c,Indices(0,0,0,0));
    fermieven.fillWithZero();	fermieven.setElement(ONE_c,Indices(1,0,1,0));
    bool fermi_odd_saved=false,fermi_even_saved=false;
    int pos_fermi_odd_saved=-1,pos_fermi_even_saved=-1;
    for(int i=(pos+1); i<(pos+length); i++)
    {
      if((i%2)==0 && fermi_even_saved) 
      {
	//Even site, already saved
	//cout << "Setting pointer at " << 3*i-3 << endl;
	item = 3*i-3;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-3);
	mpo.setOp(3*i-3,&mpo.getOp(pos_fermi_even_saved),false);
	//cout << "Setting pointer at " << 3*i-2 << endl;
	item = 3*i-2;
	setting_pos.push_back(3*i-2);
	mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_even_saved),false);
      }
      else if((i%2)==0)
      {
	//cout << "Setting operator at " << 3*i-3 << endl;
	item = 3*i-3;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-3);
	mpo.setOp(3*i-3,new Operator(fermieven),true);
	fermi_even_saved=true;
	pos_fermi_even_saved=3*i-3;
	//cout << "Setting pointer at " << 3*i-2 << endl;
	item = 3*i-2;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-2);
	mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_even_saved),false);
      }
	
      if((i%2)!=0 && fermi_odd_saved) 
      {
	//Odd site
	//cout << "Setting pointer at " << 3*i-3 << endl;
	item = 3*i-3;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-3);
	mpo.setOp(3*i-3,&mpo.getOp(pos_fermi_odd_saved),false);
	//cout << "Setting pointer at " << 3*i-2 << endl;
	item = 3*i-2;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-2);
	mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_odd_saved),false);
      }
      else if((i%2)!=0)
      {
	//cout << "Setting operator at " << 3*i-3 << endl;
	item = 3*i-3;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-3);
	mpo.setOp(3*i-3,new Operator(fermiodd),true);
	fermi_odd_saved=true;
	pos_fermi_odd_saved=3*i-3;
	//cout << "Setting pointer at " << 3*i-2 << endl;
	item = 3*i-2;
	if(find(setting_pos.begin(), setting_pos.end(), item)!=setting_pos.end())
	  cout << "----> Error" << endl;
	setting_pos.push_back(3*i-2);
	mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_odd_saved),false);
      }
    }    
  }
  
  sort(setting_pos.begin(), setting_pos.end());    
  //cout << "Setting_pos: " << setting_pos << endl;
  bool duplicat_flag=false;
  for(int i=0; i<(setting_pos.size()-1); i++)
  {
    if(setting_pos[i]==setting_pos[i+1])
      duplicat_flag=true;
  }
  if(duplicat_flag)
    cout << "--> WARNING FOUND DUPLICATE ENTRIES, method was " << fermions << ", length " << length << ", pos " << pos << endl;
  
  //cout <<"---------------------------------------------------------------" << endl;
  /*char filename[50];  
  sprintf(filename,"%s%i%s%i%s","StringProjector_leng",length,"_pos",pos,".m");  
  mpo.exportForMatlab(filename);*/
}

void SU2Hamiltonian::getMovingString_particles_antiparticles(MPS &mps,int pos_left, int pos_right, double width_left, double width_right, double k_vec, bool factor_on) const
{
  //Bond dimension
  int Dbond=15;
  
  //k-vector 
  double kvector=2.0*M_PIl/(this->N_fermions+1)*k_vec;
  
  vector<int> phys_dims(this->N_total,this->d_fermi);  
  for(int k=2; k<(this->N_total-2); k+=3)
    phys_dims[k]=this->d_link;    
  
  mps.clear();
  mps=MPS(this->N_total,Dbond,phys_dims);
  
  //Determine width and position (if not set externally or given values are out of range)
  if(pos_left<=0 || pos_left>this->N_fermions)
    pos_left=round(this->N_fermions/3.0);
  if(pos_right<=0 || pos_left>this->N_fermions)
    pos_right=round(2.0*this->N_fermions/3.0);
  if(width_left<=0)
    width_left=2.0;
  if(width_right<=0)
    width_right=2.0;
  
  cout << "Creating particle antiparticle string with" << endl
       << "pos_left:    " << pos_left << endl
       << "pos_right:   " << pos_right << endl
       << "width_left:  " << width_left << endl
       << "width_right: " << width_right << endl
       << "k_vec:       " << k_vec << endl
       << "kvector:     " << kvector << endl;
       
  //I assume that I have an even number of fermionic sites, if this is not the case I give a warning
  if((this->N_fermions%2)!=0)
    cout << "Warning in SU2Hamiltonian::getMovingString(), expected an even number of sites, instead got " << this->N_fermions << endl;
  
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(2,this->d_fermi,1));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(5,this->d_link,1));	Gauge_operators.fillWithZero();  
  
  //Spin up
  Fermionic_operators.setElement(ONE_c,Indices(0,0,0));
  //Spin down
  Fermionic_operators.setElement(ONE_c,Indices(1,1,0));
  
  Fermionic_operators.reshape(Indices(2,this->d_fermi));
  
  //State e_0
  Gauge_operators.setElement(ONE_c,Indices(0,0,0));
  //State e_1
  Gauge_operators.setElement(ONE_c,Indices(1,1,0));
  //State e_2
  Gauge_operators.setElement(ONE_c,Indices(2,2,0));
  //State e_3
  Gauge_operators.setElement(ONE_c,Indices(3,3,0));
  //State e_4
  Gauge_operators.setElement(ONE_c,Indices(4,4,0));
  
  Gauge_operators.reshape(Indices(5,this->d_link)); 
  
  //Factor for the gauge matrices (as they are not unitary, I would in principle have to put a factor of 1/sqrt(2) everytime I apply one, but this would give different weights to strings of different length in the superposition, so I have this factor exteranlly to have more control)
  double factor = 1.0/sqrt(2.0);
  if (!factor_on)
  {
    cout << "Will disable supression factor for strings!" << endl;
    factor=1.0;
  } 
  
  mwArray Fermi1First(Indices(1,Dbond,2)),Fermi1Odd(Indices(Dbond,Dbond,2)),Fermi2Odd(Indices(Dbond,Dbond,2)),Gauge(Indices(Dbond,Dbond,5)),Fermi1Even(Indices(Dbond,Dbond,2)),Fermi2Even(Indices(Dbond,Dbond,2)),Fermi2End(Indices(Dbond,1,2));
  
  Fermi1First.fillWithZero();
  Fermi1Odd.fillWithZero();
  Fermi2Odd.fillWithZero();
  Gauge.fillWithZero();
  Fermi1Even.fillWithZero();
  Fermi2Even.fillWithZero();
  Fermi2End.fillWithZero();
  
  //Set the site independent entries of the matrices
  Fermi1First.setElement(ONE_c,Indices(0,0,0));
  
  Fermi1First.setElement(ONE_c,Indices(0,2,1));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi1Odd.setElement(ONE_c,Indices(5,5,0));
  Fermi1Odd.setElement(ONE_c,Indices(6,6,0));
  Fermi1Odd.setElement(ONE_c,Indices(7,7,0));
  Fermi1Odd.setElement(ONE_c,Indices(8,8,0));
  Fermi1Odd.setElement(ONE_c,Indices(9,9,0));
  Fermi1Odd.setElement(ONE_c,Indices(10,10,0));
  Fermi1Odd.setElement(ONE_c,Indices(11,11,0));
  Fermi1Odd.setElement(ONE_c,Indices(12,12,0));
  Fermi1Odd.setElement(ONE_c,Indices(14,14,0));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,2,1));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi2Odd.setElement(ONE_c,Indices(5,5,0));
  Fermi2Odd.setElement(ONE_c,Indices(6,6,0));
  Fermi2Odd.setElement(ONE_c,Indices(7,7,0));
  Fermi2Odd.setElement(ONE_c,Indices(8,8,0));
  Fermi2Odd.setElement(ONE_c,Indices(9,9,0));
  Fermi2Odd.setElement(ONE_c,Indices(10,10,0));
  Fermi2Odd.setElement(ONE_c,Indices(11,11,0));
  Fermi2Odd.setElement(ONE_c,Indices(12,12,0));
  Fermi2Odd.setElement(ONE_c,Indices(14,14,0));
  
  Fermi2Odd.setElement(ONE_c,Indices(13,14,1));
  
  Gauge.setElement(factor*ONE_c,Indices(2,5,0));
  Gauge.setElement(factor*ONE_c,Indices(4,9,0));
  Gauge.setElement(factor*ONE_c,Indices(5,5,0));
  Gauge.setElement(factor*ONE_c,Indices(9,9,0));
  Gauge.setElement(factor*ONE_c,Indices(7,5,0));
  Gauge.setElement(factor*ONE_c,Indices(11,9,0));
  
  Gauge.setElement(1.0*factor*I_c,Indices(2,6,1));
  Gauge.setElement(-1.0*factor*I_c,Indices(4,10,1));
  Gauge.setElement(1.0*factor*I_c,Indices(5,6,1));
  Gauge.setElement(-1.0*factor*I_c,Indices(9,10,1));
  Gauge.setElement(1.0*factor*I_c,Indices(7,6,1));
  Gauge.setElement(-1.0*factor*I_c,Indices(11,10,1));
  
  Gauge.setElement(ONE_c,Indices(0,0,2));
  Gauge.setElement(ONE_c,Indices(14,14,2));
  
  Gauge.setElement(-1.0*factor*I_c,Indices(1,7,3));
  Gauge.setElement(1.0*factor*I_c,Indices(3,11,3));
  Gauge.setElement(-1.0*factor*I_c,Indices(6,7,3));
  Gauge.setElement(1.0*factor*I_c,Indices(10,11,3));
  Gauge.setElement(-1.0*factor*I_c,Indices(8,7,3));
  Gauge.setElement(1.0*factor*I_c,Indices(12,11,3));
  
  Gauge.setElement(factor*ONE_c,Indices(1,8,4));
  Gauge.setElement(factor*ONE_c,Indices(3,12,4));
  Gauge.setElement(factor*ONE_c,Indices(8,8,4));
  Gauge.setElement(factor*ONE_c,Indices(12,12,4));
  Gauge.setElement(factor*ONE_c,Indices(6,8,4));
  Gauge.setElement(factor*ONE_c,Indices(10,12,4));
  
  Fermi1Even.setElement(-1.0*ONE_c,Indices(0,3,0));
  
  Fermi1Even.setElement(ONE_c,Indices(0,0,1));
  Fermi1Even.setElement(ONE_c,Indices(5,5,1));
  Fermi1Even.setElement(ONE_c,Indices(6,6,1));
  Fermi1Even.setElement(ONE_c,Indices(7,7,1));
  Fermi1Even.setElement(ONE_c,Indices(8,8,1));
  Fermi1Even.setElement(ONE_c,Indices(9,9,1));
  Fermi1Even.setElement(ONE_c,Indices(10,10,1));
  Fermi1Even.setElement(ONE_c,Indices(11,11,1));
  Fermi1Even.setElement(ONE_c,Indices(12,12,1));
  Fermi1Even.setElement(ONE_c,Indices(14,14,1));
  
  Fermi2Even.setElement(-1.0*ONE_c,Indices(13,14,0));
  
  Fermi2Even.setElement(ONE_c,Indices(0,0,1));
  Fermi2Even.setElement(-1.0*ONE_c,Indices(5,5,1));
  Fermi2Even.setElement(-1.0*ONE_c,Indices(6,6,1));
  Fermi2Even.setElement(-1.0*ONE_c,Indices(7,7,1));
  Fermi2Even.setElement(-1.0*ONE_c,Indices(8,8,1));
  Fermi2Even.setElement(-1.0*ONE_c,Indices(9,9,1));
  Fermi2Even.setElement(-1.0*ONE_c,Indices(10,10,1));
  Fermi2Even.setElement(-1.0*ONE_c,Indices(11,11,1));
  Fermi2Even.setElement(-1.0*ONE_c,Indices(12,12,1));
  Fermi2Even.setElement(ONE_c,Indices(14,14,1));  
  
  Fermi2End.setElement(-1.0*ONE_c,Indices(13,0,0));
  
  Fermi2End.setElement(ONE_c,Indices(14,0,1));
  
  /*cout << "Fermi1First: " << Fermi1First << endl
       << "Fermi1Odd: " << Fermi1Odd << endl
       << "Fermi2Odd: " << Fermi2Odd << endl
       << "Gauge: " << Gauge << endl
       << "Fermi1Even: " << Fermi1Even << endl
       << "Fermi2Even: " << Fermi2Even << endl
       << "Fermi2End: " << Fermi2End << endl;*/
  
  //mwArray to generate the local MPS matrix
  mwArray Matrix;
  //Variable for Momentum + Gaussian distribution
  complex_t prefactor;
  
  for(int i=1; i<=this->N_fermions; i++)
  {
    if(i==1)
    {
      //Put the special vector at the beginning
      Matrix = reshape(Fermi1First,Indices(1*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(1,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-3,Matrix);
      
      //Site dependend entries
      prefactor=exp(-I_c*kvector*i-ONE_c*(i-pos_left)*(i-pos_left)/2.0/width_left/width_left);
      Fermi2Odd.setElement(prefactor,Indices(2,2,0));
      Fermi2Odd.setElement(prefactor,Indices(0,1,1));
      //Contract and set matrix
      Matrix = reshape(Fermi2Odd,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-2,Matrix);
      
      //Contract and set matrix
      Matrix = reshape(Gauge,Indices(Dbond*Dbond,5))*Gauge_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_link));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-1,Matrix);      
    }
    else if(i==this->N_fermions)
    { 
      //Site dependend entries
      prefactor=exp(I_c*kvector*i-ONE_c*(i-pos_right)*(i-pos_right)/2.0/width_right/width_right);
      Fermi1Even.setElement(prefactor,Indices(5,14,0));
      Fermi1Even.setElement(prefactor,Indices(7,14,0));
      Fermi1Even.setElement(prefactor,Indices(6,13,1));
      Fermi1Even.setElement(prefactor,Indices(8,13,1));
      //Contract and set matrix
      Matrix = reshape(Fermi1Even,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-3,Matrix);
     
      //Put the special vector at the end
      Matrix = reshape(Fermi2End,Indices(Dbond*1,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,1,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-2,Matrix);
    }
    else if((i%2)==0)
    {
      //Even site
      
      //Site dependend entries
      prefactor=exp(I_c*kvector*i-ONE_c*(i-pos_right)*(i-pos_right)/2.0/width_right/width_right);
      Fermi1Even.setElement(prefactor,Indices(5,14,0));
      Fermi1Even.setElement(prefactor,Indices(7,14,0));
      Fermi1Even.setElement(prefactor,Indices(6,13,1));
      Fermi1Even.setElement(prefactor,Indices(8,13,1));
      //Contract and set matrix
      Matrix = reshape(Fermi1Even,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-3,Matrix);
      
      //Site dependend entries
      prefactor=exp(-I_c*kvector*i-ONE_c*(i-pos_left)*(i-pos_left)/2.0/width_left/width_left);
      Fermi2Even.setElement(prefactor,Indices(0,4,0));
      Fermi2Even.setElement(prefactor,Indices(3,3,1));
      //Contract and set matrix
      Matrix = reshape(Fermi2Even,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-2,Matrix);
      
      Matrix = reshape(Gauge,Indices(Dbond*Dbond,5))*Gauge_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_link));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-1,Matrix); 
    }
    else
    {
      //Odd site
      
      //Site dependend entries
      prefactor=exp(I_c*kvector*i-ONE_c*(i-pos_right)*(i-pos_right)/2.0/width_right/width_right);
      Fermi1Odd.setElement(prefactor,Indices(9,13,0));
      Fermi1Odd.setElement(prefactor,Indices(11,13,0));
      Fermi1Odd.setElement(prefactor,Indices(10,14,1));
      Fermi1Odd.setElement(prefactor,Indices(12,14,1));
      //Contract and set matrix
      Matrix = reshape(Fermi1Odd,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-3,Matrix);
      
      //Site dependend entries
      prefactor=exp(-I_c*kvector*i-ONE_c*(i-pos_left)*(i-pos_left)/2.0/width_left/width_left);
      Fermi2Odd.setElement(prefactor,Indices(2,2,0));
      Fermi2Odd.setElement(prefactor,Indices(0,1,1));
      //Contract and set matrix
      Matrix = reshape(Fermi2Odd,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-2,Matrix);
      
      //Site dependend entries
      //Contract and set matrix
      Matrix = reshape(Gauge,Indices(Dbond*Dbond,5))*Gauge_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_link));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-1,Matrix);      
    } 
  }
  mps.gaugeCond('R',true);
}

void SU2Hamiltonian::getMovingString_particles(MPS &mps,int pos_left, int pos_right, double width_left, double width_right, double k_vec ,bool factor_on) const
{
  //Bond dimension
  int Dbond=8;
  
  //k-vector 
  double kvector=2.0*M_PIl/(this->N_fermions+1)*k_vec;
  
  vector<int> phys_dims(this->N_total,this->d_fermi);  
  for(int k=2; k<(this->N_total-2); k+=3)
    phys_dims[k]=this->d_link;    
  
  mps.clear();
  mps=MPS(this->N_total,Dbond,phys_dims);
  
  //Determine width and position (if not set externally or given values are out of range)
  if(pos_left<=0 || pos_left>this->N_fermions)
    pos_left=round(this->N_fermions/4.0);
  if(pos_right<=0 || pos_left>this->N_fermions)
    pos_right=round(3.0*this->N_fermions/4.0);
  if(width_left<=0)
    width_left=2.0;
  if(width_right<=0)
    width_right=2.0;
  
  cout << "Creating particle string with " << endl
       << "pos_left:    " << pos_left << endl
       << "pos_right:   " << pos_right << endl
       << "width_left:  " << width_left << endl
       << "width_right: " << width_right << endl
       << "k_vec:       " << k_vec << endl
       << "kvector:     " << kvector << endl;
       
  cout << "WARNING: Using randomized kvector for every site!" << endl;
  //Seed the random number generator 
  srand(time(NULL));
       
  //I assume that I have an even number of fermionic sites, if this is not the case I give a warning
  if((this->N_fermions%2)!=0)
    cout << "Warning in SU2Hamiltonian::getMovingString(), expected an even number of sites, instead got " << this->N_fermions << endl;
  
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(2,this->d_fermi,1));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(5,this->d_link,1));	Gauge_operators.fillWithZero();  
  
  //Spin up
  Fermionic_operators.setElement(ONE_c,Indices(0,0,0));
  //Spin down
  Fermionic_operators.setElement(ONE_c,Indices(1,1,0));
  
  Fermionic_operators.reshape(Indices(2,this->d_fermi));
  
  //State e_0
  Gauge_operators.setElement(ONE_c,Indices(0,0,0));
  //State e_1
  Gauge_operators.setElement(ONE_c,Indices(1,1,0));
  //State e_2
  Gauge_operators.setElement(ONE_c,Indices(2,2,0));
  //State e_3
  Gauge_operators.setElement(ONE_c,Indices(3,3,0));
  //State e_4
  Gauge_operators.setElement(ONE_c,Indices(4,4,0));
  
  Gauge_operators.reshape(Indices(5,this->d_link)); 
  
  //Factor for the gauge matrices (as they are not unitary, I would in principle have to put a factor of 1/sqrt(2) everytime I apply one, but this would give different weights to strings of different length in the superposition, so I have this factor exteranlly to have more control)
  double factor = 1.0/sqrt(2.0);
  if (!factor_on)
  {
    cout << "Will disable supression factor for strings!" << endl;
    factor=1.0;
  }
  
  mwArray Fermi1First(Indices(1,Dbond,2)),Fermi1Odd(Indices(Dbond,Dbond,2)),Fermi2Odd(Indices(Dbond,Dbond,2)),Gauge(Indices(Dbond,Dbond,5)),Fermi1Even(Indices(Dbond,Dbond,2)),Fermi2Even(Indices(Dbond,Dbond,2)),Fermi2End(Indices(Dbond,1,2));
  
  Fermi1First.fillWithZero();
  Fermi1Odd.fillWithZero();
  Fermi2Odd.fillWithZero();
  Gauge.fillWithZero();
  Fermi1Even.fillWithZero();
  Fermi2Even.fillWithZero();
  Fermi2End.fillWithZero();
  
  //Set the site independent entries of the matrices
  Fermi1First.setElement(ONE_c,Indices(0,0,0));
  
  Fermi1First.setElement(ONE_c,Indices(0,1,1));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi1Odd.setElement(ONE_c,Indices(3,3,0));
  Fermi1Odd.setElement(ONE_c,Indices(4,4,0));
  Fermi1Odd.setElement(ONE_c,Indices(5,5,0));
  Fermi1Odd.setElement(ONE_c,Indices(6,6,0));
  Fermi1Odd.setElement(ONE_c,Indices(7,7,0));
  
  Fermi1Odd.setElement(ONE_c,Indices(0,1,1));
  
  Fermi2Odd.setElement(ONE_c,Indices(0,0,0));
  Fermi2Odd.setElement(ONE_c,Indices(3,3,0));
  Fermi2Odd.setElement(ONE_c,Indices(4,4,0));
  Fermi2Odd.setElement(ONE_c,Indices(5,5,0));
  Fermi2Odd.setElement(ONE_c,Indices(6,6,0));
  Fermi2Odd.setElement(ONE_c,Indices(7,7,0));
  
  Gauge.setElement(factor*ONE_c,Indices(1,3,0)); 
  Gauge.setElement(factor*ONE_c,Indices(3,3,0));
  Gauge.setElement(factor*ONE_c,Indices(5,3,0));
  
  Gauge.setElement(1.0*factor*I_c,Indices(1,4,1));
  Gauge.setElement(1.0*factor*I_c,Indices(3,4,1));
  Gauge.setElement(1.0*factor*I_c,Indices(5,4,1));
  
  Gauge.setElement(ONE_c,Indices(0,0,2));
  Gauge.setElement(ONE_c,Indices(7,7,2));
 
  Gauge.setElement(-1.0*factor*I_c,Indices(2,5,3));
  Gauge.setElement(-1.0*factor*I_c,Indices(4,5,3));
  Gauge.setElement(-1.0*factor*I_c,Indices(6,5,3));
  
  Gauge.setElement(factor*ONE_c,Indices(2,6,4));
  Gauge.setElement(factor*ONE_c,Indices(4,6,4));
  Gauge.setElement(factor*ONE_c,Indices(6,6,4));  
  
  Fermi1Even.setElement(ONE_c,Indices(0,0,1));
  Fermi1Even.setElement(ONE_c,Indices(3,3,1));
  Fermi1Even.setElement(ONE_c,Indices(4,4,1));
  Fermi1Even.setElement(ONE_c,Indices(5,5,1));
  Fermi1Even.setElement(ONE_c,Indices(6,6,1));
  Fermi1Even.setElement(ONE_c,Indices(7,7,1));
  
  Fermi2Even.setElement(ONE_c,Indices(0,0,1));
  Fermi2Even.setElement(-1.0*ONE_c,Indices(3,3,1));
  Fermi2Even.setElement(-1.0*ONE_c,Indices(4,4,1));
  Fermi2Even.setElement(-1.0*ONE_c,Indices(5,5,1));
  Fermi2Even.setElement(-1.0*ONE_c,Indices(6,6,1));
  Fermi2Even.setElement(ONE_c,Indices(7,7,1));
  
  Fermi2End.setElement(ONE_c,Indices(7,0,1));   
  
  /*cout << "Fermi1First: " << Fermi1First << endl
       << "Fermi1Odd: " << Fermi1Odd << endl
       << "Fermi2Odd: " << Fermi2Odd << endl
       << "Gauge: " << Gauge << endl
       << "Fermi1Even: " << Fermi1Even << endl
       << "Fermi2Even: " << Fermi2Even << endl
       << "Fermi2End: " << Fermi2End << endl;*/
  
  //mwArray to generate the local MPS matrix
  mwArray Matrix;
  //Variable for Momentum + Gaussian distribution
  complex_t prefactor;
  
  for(int i=1; i<=this->N_fermions; i++)
  {
    k_vec=rand();
    kvector=2.0*M_PIl/(this->N_fermions+1)*k_vec;    
    if(i==1)
    {
      //Put the special vector at the beginning
      Matrix = reshape(Fermi1First,Indices(1*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(1,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-3,Matrix);
      
      //Site dependend entries
      prefactor=exp(-I_c*kvector*i-ONE_c*(i-pos_left)*(i-pos_left)/2.0/width_left/width_left);
      //In case I want both sorts moving in the same direction
      //prefactor=exp(I_c*kvector*i-ONE_c*(i-pos_left)*(i-pos_left)/2.0/width_left/width_left);
      Fermi2Odd.setElement(prefactor,Indices(1,1,0));
      Fermi2Odd.setElement(prefactor,Indices(0,2,1));
      //Contract and set matrix
      Matrix = reshape(Fermi2Odd,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-2,Matrix);
      
      //Contract and set matrix
      Matrix = reshape(Gauge,Indices(Dbond*Dbond,5))*Gauge_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_link));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-1,Matrix);      
    }
    else if(i==this->N_fermions)
    { 
      //Site dependend entries
      prefactor=exp(I_c*kvector*i-ONE_c*(i-pos_right)*(i-pos_right)/2.0/width_right/width_right);
      Fermi1Even.setElement(prefactor,Indices(3,7,0));
      Fermi1Even.setElement(prefactor,Indices(5,7,0));
      //Contract and set matrix
      Matrix = reshape(Fermi1Even,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-3,Matrix);
     
      //Site dependend entries
      prefactor=exp(I_c*kvector*i-ONE_c*(i-pos_right)*(i-pos_right)/2.0/width_right/width_right);
      Fermi2End.setElement(-1.0*prefactor,Indices(4,0,0));
      Fermi2End.setElement(-1.0*prefactor,Indices(6,0,0));
      //Contract and set special matrix
      Matrix = reshape(Fermi2End,Indices(Dbond*1,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,1,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-2,Matrix);
    }
    else if((i%2)==0)
    {
      //Even site
      
      //Site dependend entries
      prefactor=exp(I_c*kvector*i-ONE_c*(i-pos_right)*(i-pos_right)/2.0/width_right/width_right);
      Fermi1Even.setElement(prefactor,Indices(3,7,0));
      Fermi1Even.setElement(prefactor,Indices(5,7,0));
      //Contract and set matrix
      Matrix = reshape(Fermi1Even,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-3,Matrix);
      
      //Site dependend entries
      Fermi2Even.setElement(-1.0*prefactor,Indices(4,7,0));
      Fermi2Even.setElement(-1.0*prefactor,Indices(6,7,0));
      //Contract and set matrix
      Matrix = reshape(Fermi2Even,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-2,Matrix);
      
      Matrix = reshape(Gauge,Indices(Dbond*Dbond,5))*Gauge_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_link));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-1,Matrix); 
    }
    else
    {
      //Odd site
      
      //Contract and set matrix
      Matrix = reshape(Fermi1Odd,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-3,Matrix);
      
      //Site dependend entries
      prefactor=exp(-I_c*kvector*i-ONE_c*(i-pos_left)*(i-pos_left)/2.0/width_left/width_left);
      //In case I want both sorts moving in the same direction
      //prefactor=exp(I_c*kvector*i-ONE_c*(i-pos_left)*(i-pos_left)/2.0/width_left/width_left);
      Fermi2Odd.setElement(prefactor,Indices(1,1,0));
      Fermi2Odd.setElement(prefactor,Indices(0,2,1));
      //Contract and set matrix
      Matrix = reshape(Fermi2Odd,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-2,Matrix);
      
      //Site dependend entries
      //Contract and set matrix
      Matrix = reshape(Gauge,Indices(Dbond*Dbond,5))*Gauge_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_link));
      Matrix.permute(Indices(3,1,2));
      mps.setA(3*i-1,Matrix);      
    } 
  }
  //Normalize the MPS
  mps.gaugeCond('R',true);
}

void SU2Hamiltonian::getMovingString(MPS &mps,bool antiparticles,int pos_left, int pos_right, double width_left, double width_right,double k_vec,bool factor_on) const
{
  if(antiparticles)
    this->getMovingString_particles_antiparticles(mps,pos_left,pos_right,width_left,width_right,k_vec,factor_on);
  else
    this->getMovingString_particles(mps,pos_left,pos_right,width_left,width_right,k_vec,factor_on);
}

void SU2Hamiltonian::getMovingMeson(MPO &mpo,int pos, double width, double k_vec) const
{
  //Bond dimension
  int D_bond=7;
  
  //k-vector 
  double kvector=M_PIl*k_vec;
  
  //Initialize the MPO
  mpo.initLength(this->N_total);

  
  //Determine width and position (if not set externally or given values are out of range)
  if(pos<=0 || pos>this->N_fermions)
    pos=round(this->N_fermions/2.0);
  if(width<=0)
    width=2.0;
  
  cout << "Creating MPO for meson superposition with " << endl
       << " pos:     " << pos << endl       
       << " width:   " << width << endl       
       << " k_vec:   " << k_vec << endl
       << " kvector: " << kvector << endl;
       
       
  //I assume that I have an even number of fermionic sites, if this is not the case I give a warning
  if((this->N_fermions%2)!=0)
    cout << "Warning in SU2Hamiltonian::getMovingMeson(), expected an even number of sites, instead got " << this->N_fermions << endl;
  
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(4,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(5,this->d_link,this->d_link));	Gauge_operators.fillWithZero();  
  
  for(int i=0; i<this->d_fermi; i++)
  {
    for(int j=0; j<this->d_fermi; j++)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma plus
      Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
      //sigma minus
      Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
      //sigma z
      Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(3,i,j));
    }
  }
  for(int i=0; i<this->d_link; i++)
  {
    for(int j=0; j<this->d_link; j++)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //U00^\dagger
      Gauge_operators.setElement(this->U00ad.getElement(Indices(i,j)),Indices(1,i,j));
      //U01^\dagger
      Gauge_operators.setElement(this->U01ad.getElement(Indices(i,j)),Indices(2,i,j));
      //U10^\dagger
      Gauge_operators.setElement(this->U10ad.getElement(Indices(i,j)),Indices(3,i,j));
      //U11^\dagger
      Gauge_operators.setElement(this->U11ad.getElement(Indices(i,j)),Indices(4,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(4,this->d_fermi*this->d_fermi));
  Gauge_operators.reshape(Indices(5,this->d_link*this->d_link));
  
  
  //Matrices for the finite automata
  mwArray Fermi1odd,Fermi2odd,Fermi1even,Fermi2even,Gauge,Fermi1first,Fermi2last;
  Fermi1first=mwArray(Indices(1,D_bond,4));	Fermi1first.fillWithZero();
  Fermi1odd=mwArray(Indices(D_bond,D_bond,4));	Fermi1odd.fillWithZero();
  Fermi2odd=mwArray(Indices(D_bond,D_bond,4));	Fermi2odd.fillWithZero();
  Fermi1even=mwArray(Indices(D_bond,D_bond,4));	Fermi1even.fillWithZero();
  Fermi2even=mwArray(Indices(D_bond,D_bond,4));	Fermi2even.fillWithZero();
  Fermi2last=mwArray(Indices(D_bond,1,4));	Fermi2last.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,5));	Gauge.fillWithZero();
  
  //Set the (site independent) entries
  Fermi1first.setElement(ONE_c,Indices(0,0,0));
  Fermi1first.setElement(ONE_c,Indices(0,1,2));
  
  Fermi1odd.setElement(ONE_c,Indices(0,0,0));
  Fermi1odd.setElement(ONE_c,Indices(6,6,0));  
  
  Fermi1odd.setElement(ONE_c,Indices(0,1,2));
  
  Fermi2odd.setElement(ONE_c,Indices(0,0,0));
  Fermi2odd.setElement(ONE_c,Indices(6,6,0));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(6,6,0));
  
  Gauge.setElement(ONE_c,Indices(1,3,1));
  Gauge.setElement(ONE_c,Indices(3,3,1));
  
  Gauge.setElement(-I_c,Indices(1,4,2));
  Gauge.setElement(-I_c,Indices(3,4,2));
  
  Gauge.setElement(I_c,Indices(2,3,3));
  Gauge.setElement(I_c,Indices(4,3,3));
  
  Gauge.setElement(ONE_c,Indices(2,4,4));
  Gauge.setElement(ONE_c,Indices(4,4,4));
  
  Fermi1even.setElement(ONE_c,Indices(0,0,0));
  Fermi1even.setElement(ONE_c,Indices(6,6,0));
  
  Fermi1even.setElement(ONE_c,Indices(3,6,1));
  
  Fermi1even.setElement(ONE_c,Indices(4,5,3));
  
  Fermi2even.setElement(ONE_c,Indices(0,0,0));
  Fermi2even.setElement(ONE_c,Indices(6,6,0));
  
  Fermi2even.setElement(ONE_c,Indices(5,6,1));
  
  Fermi2last.setElement(ONE_c,Indices(6,0,0));
  Fermi2last.setElement(ONE_c,Indices(5,0,1));
    
  //Contract and fill the MPO
  mwArray res;
  Operator* Local_op;
  double gauss;
  complex_t momentum;
  //cout << "Starting to fill MPO" << endl;
  for(int i=1; i<=this->N_fermions; i++)
  {
    //cout << "--> Site " << i << endl;
    if(i==1)
    {
      //Set site depenend entries for first site (always odd)
      gauss = exp(-1.0/2.0/width/width*pow(i-pos,2.0));
      momentum=exp(-I_c*kvector*i);
      Fermi2odd.setElement(gauss*momentum,Indices(0,2,2));
      Fermi2odd.setElement(gauss*momentum,Indices(1,1,3));
      //Contract and set entries
      res = reshape(Fermi1first,Indices(1*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(1,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-3,Local_op,true);
      
      res = reshape(Fermi2odd,Indices(D_bond*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-2,Local_op,true);
      
      res = reshape(Gauge,Indices(D_bond*D_bond,5))*Gauge_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-1,Local_op,true);     
    }
    else if(i==this->N_fermions)
    {
      //Contract and set entries
      res = reshape(Fermi1even,Indices(D_bond*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-3,Local_op,true);
      
      res = reshape(Fermi2last,Indices(D_bond*1,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,1,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-2,Local_op,true);
    }
    else if((i%2)==0)
    {
      //Contract and set entries
      res = reshape(Fermi1even,Indices(D_bond*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-3,Local_op,true);
      
      res = reshape(Fermi2even,Indices(D_bond*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-2,Local_op,true);
      
      res = reshape(Gauge,Indices(D_bond*D_bond,5))*Gauge_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-1,Local_op,true);
    }
    else
    {
      //Set site dependend entries for odd site 
      gauss = exp(-1.0/2.0/width/width*pow(i-pos,2.0));
      momentum=exp(-I_c*kvector*i);
      Fermi2odd.setElement(gauss*momentum,Indices(0,2,2));
      Fermi2odd.setElement(gauss*momentum,Indices(1,1,3));
      //Contract and set entries
      res = reshape(Fermi1odd,Indices(D_bond*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-3,Local_op,true);
      
      res = reshape(Fermi2odd,Indices(D_bond*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-2,Local_op,true);
      
      res = reshape(Gauge,Indices(D_bond*D_bond,5))*Gauge_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-1,Local_op,true); 
    }
  }
}

void SU2Hamiltonian::getMovingAntimesonMeson(MPO &mpo,int pos, double width, double k_vec, bool antimesons, bool mesons) const
{
  //Bond dimension
  int Dbond=6;
  
  if(!(antimesons || mesons))
  {
    //Cast that none of them is selected
    cout << "Warning in SU2Hamiltonian::getMovingAntimesonMeson(), must at least have one type of mesons, will take antimesons" << endl;
    antimesons=true;
  }
   
  //k-vector 
  complex_t phase=-I_c*M_PIl*k_vec;
  
  //Initialize the MPO
  mpo.initLength(this->N_total);
  
  //Determine width and position (if not set externally or given values are out of range)
  if(pos<=0 || pos>this->N_fermions)
    pos=round(this->N_fermions/2.0);
  if(width<=0)
    width=2.0;
  
  cout << "Creating MPO for antimeson meson superposition with " << endl
       << " pos:     " << pos << endl       
       << " width:   " << width << endl       
       << " k_vec:   " << k_vec << endl
       << " kvector: " << k_vec*M_PIl << endl;      
  if(antimesons)
    cout << "Taking into account antimesons" << endl;
  if(mesons)
    cout << "Taking into account mesons" << endl;
       
  //I assume that I have an even number of fermionic sites, if this is not the case I give a warning
  if((this->N_fermions%2)!=0)
    cout << "Warning in SU2Hamiltonian::getMovingAntimesonMeson(), expected an even number of sites, instead got " << this->N_fermions << endl;
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(4,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(9,this->d_link,this->d_link));	Gauge_operators.fillWithZero();  
  
  for(int i=0; i<this->d_fermi; i++)
  {
    for(int j=0; j<this->d_fermi; j++)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma plus
      Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
      //sigma minus
      Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
      //sigma z
      Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(3,i,j));
    }
  }
  for(int i=0; i<this->d_link; i++)
  {
    for(int j=0; j<this->d_link; j++)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //U00
      Gauge_operators.setElement(this->U00.getElement(Indices(i,j)),Indices(1,i,j));
      //U01
      Gauge_operators.setElement(this->U01.getElement(Indices(i,j)),Indices(2,i,j));
      //U10
      Gauge_operators.setElement(this->U10.getElement(Indices(i,j)),Indices(3,i,j));
      //U11
      Gauge_operators.setElement(this->U11.getElement(Indices(i,j)),Indices(4,i,j));
      //U00^\dagger
      Gauge_operators.setElement(this->U00ad.getElement(Indices(i,j)),Indices(5,i,j));
      //U01^\dagger
      Gauge_operators.setElement(this->U01ad.getElement(Indices(i,j)),Indices(6,i,j));
      //U10^\dagger
      Gauge_operators.setElement(this->U10ad.getElement(Indices(i,j)),Indices(7,i,j));
      //U11^\dagger
      Gauge_operators.setElement(this->U11ad.getElement(Indices(i,j)),Indices(8,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(4,this->d_fermi*this->d_fermi));
  Gauge_operators.reshape(Indices(9,this->d_link*this->d_link));
  
  //Matrices 
  mwArray Fermi1first, Fermi1odd, Fermi2odd, Fermi1even, Fermi2even, Fermi2last, Gaugeodd, Gaugeeven;
  
  Fermi1first=mwArray(Indices(1,Dbond,4));	Fermi1first.fillWithZero();
  Fermi1odd=mwArray(Indices(Dbond,Dbond,4));	Fermi1odd.fillWithZero();
  Fermi2odd=mwArray(Indices(Dbond,Dbond,4));	Fermi2odd.fillWithZero();
  Fermi1even=mwArray(Indices(Dbond,Dbond,4));	Fermi1even.fillWithZero();
  Fermi2even=mwArray(Indices(Dbond,Dbond,4));	Fermi2even.fillWithZero();
  Fermi2last=mwArray(Indices(Dbond,1,4));	Fermi2last.fillWithZero();
  Gaugeodd=mwArray(Indices(Dbond,Dbond,9));	Gaugeodd.fillWithZero();
  Gaugeeven=mwArray(Indices(Dbond,Dbond,9));	Gaugeeven.fillWithZero();
  
  //Fermi1first
  Fermi1first.setElement(ONE_c,Indices(0,0,0));
  Fermi1first.setElement(ONE_c,Indices(0,1,2));
  //Fermi1Odd
  Fermi1odd.setElement(ONE_c,Indices(0,0,0));
  Fermi1odd.setElement(ONE_c,Indices(5,5,0));
  Fermi1odd.setElement(ONE_c,Indices(0,1,2));
  Fermi1odd.setElement(ONE_c,Indices(3,5,2));
  Fermi1odd.setElement(ONE_c,Indices(4,4,3));
  //Fermi2Odd
  Fermi2odd.setElement(ONE_c,Indices(0,0,0));
  Fermi2odd.setElement(ONE_c,Indices(5,5,0));
  Fermi2odd.setElement(ONE_c,Indices(0,2,2));
  Fermi2odd.setElement(ONE_c,Indices(4,5,2));
  Fermi2odd.setElement(ONE_c,Indices(1,1,3));
  //Fermi1Even
  Fermi1even.setElement(ONE_c,Indices(0,0,0));
  Fermi1even.setElement(ONE_c,Indices(5,5,0));
  Fermi1even.setElement(ONE_c,Indices(0,1,1));
  Fermi1even.setElement(ONE_c,Indices(3,5,1));
  Fermi1even.setElement(ONE_c,Indices(4,4,3));
  //Fermi2Even
  Fermi2even.setElement(ONE_c,Indices(0,0,0));
  Fermi2even.setElement(ONE_c,Indices(5,5,0));
  Fermi2even.setElement(ONE_c,Indices(0,2,1));
  Fermi2even.setElement(ONE_c,Indices(4,5,1));
  Fermi2even.setElement(ONE_c,Indices(1,1,3));
  //Fermi2last
  Fermi2last.setElement(ONE_c,Indices(5,0,0));
  Fermi2last.setElement(ONE_c,Indices(4,0,1));
  //Site independent entries gauge odd and even
  Gaugeodd.setElement(ONE_c,Indices(0,0,0));
  Gaugeodd.setElement(ONE_c,Indices(5,5,0));
  Gaugeeven.setElement(ONE_c,Indices(0,0,0));
  Gaugeeven.setElement(ONE_c,Indices(5,5,0));
  
  //Reshape the fermionic matrices as they are site independent
  Fermi1first.reshape(Indices(1*Dbond,4));	Fermi1first=reshape(Fermi1first*Fermionic_operators,Indices(1,Dbond,this->d_fermi,this->d_fermi));	Fermi1first.permute(Indices(3,1,4,2));
  Fermi1odd.reshape(Indices(Dbond*Dbond,4));	Fermi1odd=reshape(Fermi1odd*Fermionic_operators,Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Fermi1odd.permute(Indices(3,1,4,2));
  Fermi2odd.reshape(Indices(Dbond*Dbond,4));	Fermi2odd=reshape(Fermi2odd*Fermionic_operators,Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Fermi2odd.permute(Indices(3,1,4,2));
  Fermi1even.reshape(Indices(Dbond*Dbond,4));	Fermi1even=reshape(Fermi1even*Fermionic_operators,Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Fermi1even.permute(Indices(3,1,4,2));
  Fermi2even.reshape(Indices(Dbond*Dbond,4));	Fermi2even=reshape(Fermi2even*Fermionic_operators,Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	Fermi2even.permute(Indices(3,1,4,2));
  Fermi2last.reshape(Indices(Dbond*1,4));	Fermi2last=reshape(Fermi2last*Fermionic_operators,Indices(Dbond,1,this->d_fermi,this->d_fermi));	Fermi2last.permute(Indices(3,1,4,2));
  
  //Set entries of MPO
  double gauss;
  complex_t momentum;
  mwArray res;
  bool even_saved=false,odd_saved=false;
  int pos_even_saved,pos_odd_saved;
  
  //First and last one
  mpo.setOp(0,new Operator(Fermi1first),true);
  mpo.setOp(this->N_total-1,new Operator(Fermi2last),true);
  
  //Fermi1
  for(int i=2; i<=this->N_fermions; i++)
  {
    if((i%2)==0)
    {
      //Even site
      if(even_saved)
	mpo.setOp(3*i-3,&mpo.getOp(pos_even_saved),false);
      else
      {
	mpo.setOp(3*i-3,new Operator(Fermi1even),true);
	even_saved=true;
	pos_even_saved=3*i-3;
      }
    }
    else
    {
      //Odd site
      if(odd_saved)
	mpo.setOp(3*i-3,&mpo.getOp(pos_odd_saved),false);
      else
      {
	mpo.setOp(3*i-3,new Operator(Fermi1odd),true);
	odd_saved=true;
	pos_odd_saved=3*i-3;
      }      
    }
  }
  //Fermi2
  even_saved=false; 
  odd_saved=false;
  for(int i=1; i<this->N_fermions; i++)
  {
    if((i%2)==0)
    {
      //Even site
      if(even_saved)
	mpo.setOp(3*i-2,&mpo.getOp(pos_even_saved),false);
      else
      {
	mpo.setOp(3*i-2,new Operator(Fermi2even),true);
	even_saved=true;
	pos_even_saved=3*i-2;
      }
    }
    else
    {
      //Odd site
      if(odd_saved)
	mpo.setOp(3*i-2,&mpo.getOp(pos_odd_saved),false);
      else
      {
	mpo.setOp(3*i-2,new Operator(Fermi2odd),true);
	odd_saved=true;
	pos_odd_saved=3*i-2;
      }   
    }
  }
  //Gauge links
  for(int i=1; i<this->N_fermions; i++)
  {
    gauss=exp(-0.5/width/width*pow(i-pos,2.0));
    momentum=exp(phase*((double) i));
    //cout << "Factor: " << gauss*momentum << ", abs=" << abs(gauss*momentum)<<  endl;
    if((i%2)==0)
    {
      //Even site
      if(antimesons)
      {
	Gaugeeven.setElement(gauss*momentum,Indices(1,3,1));
	Gaugeeven.setElement(I_c*gauss*momentum,Indices(1,4,2));
	Gaugeeven.setElement(-I_c*gauss*momentum,Indices(2,3,3));
	Gaugeeven.setElement(gauss*momentum,Indices(2,4,4));
      }
      res = reshape(Gaugeeven,Indices(Dbond*Dbond,9))*Gauge_operators;
      res.reshape(Indices(Dbond,Dbond,this->d_link,this->d_link));
      res.permute(Indices(3,1,4,2));
    }
    else      
    {
      //Odd site
      if(mesons)
      {
	Gaugeodd.setElement(gauss*momentum,Indices(1,3,5));
	Gaugeodd.setElement(-I_c*gauss*momentum,Indices(1,4,6));
	Gaugeodd.setElement(I_c*gauss*momentum,Indices(2,3,7));
	Gaugeodd.setElement(gauss*momentum,Indices(2,4,8));
      }
      res = reshape(Gaugeodd,Indices(Dbond*Dbond,9))*Gauge_operators;
      res.reshape(Indices(Dbond,Dbond,this->d_link,this->d_link));
      res.permute(Indices(3,1,4,2));
    }
    mpo.setOp(3*i-1,new Operator(res),true);
  }  
}

void SU2Hamiltonian::getMovingString(MPO &mpo, int pos_left, int pos_right, double width_left, double width_right,double k_vec,bool factor_on) const
{
  //Bond dimension
  int D_bond=7;
  
  //k-vector 
  double kvector=2.0*M_PIl/this->N_fermions*k_vec;
  
  //Initialize the MPO
  mpo.initLength(this->N_total);

  
  //Determine width and position (if not set externally or given values are out of range)
  if(pos_left<=0 || pos_left>this->N_fermions)
    pos_left=round(this->N_fermions/4.0);
  if(pos_right<=0 || pos_left>this->N_fermions)
    pos_right=round(3.0*this->N_fermions/4.0);
  if(width_left<=0)
    width_left=2.0;
  if(width_right<=0)
    width_right=2.0;
  
  cout << "Creating MPO for string superposition with " << endl
       << " pos_left:    " << pos_left << endl
       << " pos_right:   " << pos_right << endl
       << " width_left:  " << width_left << endl
       << " width_right: " << width_right << endl
       << " k_vec:       " << k_vec << endl
       << " kvector:     " << kvector << endl;
       
       
  //I assume that I have an even number of fermionic sites, if this is not the case I give a warning
  if((this->N_fermions%2)!=0)
    cout << "Warning in SU2Hamiltonian::getMovingString(), expected an even number of sites, instead got " << this->N_fermions << endl;
  
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(4,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(5,this->d_link,this->d_link));	Gauge_operators.fillWithZero();  
  
  for(int i=0; i<this->d_fermi; i++)
  {
    for(int j=0; j<this->d_fermi; j++)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma plus
      Fermionic_operators.setElement(this->sigma_plus.getElement(Indices(i,j)),Indices(1,i,j));
      //sigma minus
      Fermionic_operators.setElement(this->sigma_minus.getElement(Indices(i,j)),Indices(2,i,j));
      //sigma z
      Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(3,i,j));
    }
  }
  for(int i=0; i<this->d_link; i++)
  {
    for(int j=0; j<this->d_link; j++)
    {
      //Identity
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
      //U00^\dagger
      Gauge_operators.setElement(this->U00ad.getElement(Indices(i,j)),Indices(1,i,j));
      //U01^\dagger
      Gauge_operators.setElement(this->U01ad.getElement(Indices(i,j)),Indices(2,i,j));
      //U10^\dagger
      Gauge_operators.setElement(this->U10ad.getElement(Indices(i,j)),Indices(3,i,j));
      //U11^\dagger
      Gauge_operators.setElement(this->U11ad.getElement(Indices(i,j)),Indices(4,i,j));
    }
  }
  Fermionic_operators.reshape(Indices(4,this->d_fermi*this->d_fermi));
  Gauge_operators.reshape(Indices(5,this->d_link*this->d_link));
  
  //Factor for the gauge matrices I can put in front of the gauge matrices to compensate the 1/sqrt(2) in them such that longer strings wont be suppressed
  double factor = 1.0;
  if (!factor_on)
  {
    cout << "Will disable supression factor for strings!" << endl;
    factor=sqrt(2.0);
  }
  
  //Matrices for the finite automata
  mwArray Fermi1odd,Fermi2odd,Fermi1even,Fermi2even,Gauge,Fermi1first,Fermi2last;
  Fermi1first=mwArray(Indices(1,D_bond,4));	Fermi1first.fillWithZero();
  Fermi1odd=mwArray(Indices(D_bond,D_bond,4));	Fermi1odd.fillWithZero();
  Fermi2odd=mwArray(Indices(D_bond,D_bond,4));	Fermi2odd.fillWithZero();
  Fermi1even=mwArray(Indices(D_bond,D_bond,4));	Fermi1even.fillWithZero();
  Fermi2even=mwArray(Indices(D_bond,D_bond,4));	Fermi2even.fillWithZero();
  Fermi2last=mwArray(Indices(D_bond,1,4));	Fermi2last.fillWithZero();
  Gauge=mwArray(Indices(D_bond,D_bond,5));	Gauge.fillWithZero();
  
  //Set the (site independent) entries
  Fermi1first.setElement(ONE_c,Indices(0,0,0));
  Fermi1first.setElement(ONE_c,Indices(0,1,2));
  
  Fermi1odd.setElement(ONE_c,Indices(0,0,0));
  Fermi1odd.setElement(ONE_c,Indices(3,3,0));
  Fermi1odd.setElement(ONE_c,Indices(4,4,0));
  Fermi1odd.setElement(ONE_c,Indices(6,6,0));  
  
  Fermi1odd.setElement(ONE_c,Indices(0,1,2));
  
  Fermi2odd.setElement(ONE_c,Indices(0,0,0));
  Fermi2odd.setElement(ONE_c,Indices(3,3,0));
  Fermi2odd.setElement(ONE_c,Indices(4,4,0));
  Fermi2odd.setElement(ONE_c,Indices(6,6,0));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(6,6,0));
  
  Gauge.setElement(factor*ONE_c,Indices(1,3,1));
  Gauge.setElement(factor*ONE_c,Indices(3,3,1));
  
  Gauge.setElement(-I_c*factor,Indices(1,4,2));
  Gauge.setElement(-I_c*factor,Indices(3,4,2));
  
  Gauge.setElement(I_c*factor,Indices(2,3,3));
  Gauge.setElement(I_c*factor,Indices(4,3,3));
  
  Gauge.setElement(factor*ONE_c,Indices(2,4,4));
  Gauge.setElement(factor*ONE_c,Indices(4,4,4));
  
  Fermi1even.setElement(ONE_c,Indices(0,0,0));
  Fermi1even.setElement(ONE_c,Indices(3,3,0));
  Fermi1even.setElement(ONE_c,Indices(4,4,0));
  Fermi1even.setElement(ONE_c,Indices(6,6,0));
  
  Fermi2even.setElement(ONE_c,Indices(0,0,0));
  Fermi2even.setElement(ONE_c,Indices(3,3,0));
  Fermi2even.setElement(ONE_c,Indices(4,4,0));
  Fermi2even.setElement(ONE_c,Indices(6,6,0));
  
  Fermi2even.setElement(ONE_c,Indices(5,6,1));
  
  Fermi2last.setElement(ONE_c,Indices(6,0,0));
  Fermi2last.setElement(ONE_c,Indices(5,0,1));
    
  //Contract and fill the MPO
  mwArray res;
  Operator* Local_op;
  double gauss;
  complex_t momentum;
  //cout << "Starting to fill MPO" << endl;
  for(int i=1; i<=this->N_fermions; i++)
  {
    //cout << "--> Site " << i << endl;
    if(i==1)
    {
      //Set site depenend entries for first site (always odd)
      gauss = exp(-1.0/2.0/width_left/width_left*pow(i-pos_left,2.0));
      momentum=exp(-I_c*kvector*i);
      Fermi2odd.setElement(gauss*momentum,Indices(0,2,2));
      Fermi2odd.setElement(gauss*momentum,Indices(1,1,3));
      //Contract and set entries
      res = reshape(Fermi1first,Indices(1*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(1,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-3,Local_op,true);
      
      res = reshape(Fermi2odd,Indices(D_bond*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-2,Local_op,true);
      
      res = reshape(Gauge,Indices(D_bond*D_bond,5))*Gauge_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-1,Local_op,true);     
    }
    else if(i==this->N_fermions)
    {
      //Set site dependend entries for last site (always even)
      gauss = exp(-1.0/2.0/width_right/width_right*pow(i-pos_right,2.0));
      momentum=exp(I_c*kvector*i);
      Fermi1even.setElement(gauss*momentum,Indices(3,6,1));
      Fermi1even.setElement(gauss*momentum,Indices(4,5,3));
      //Contract and set entries
      res = reshape(Fermi1even,Indices(D_bond*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-3,Local_op,true);
      
      res = reshape(Fermi2last,Indices(D_bond*1,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,1,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-2,Local_op,true);
    }
    else if((i%2)==0)
    {
      //Set site dependend entries for even site 
      gauss = exp(-1.0/2.0/width_right/width_right*pow(i-pos_right,2.0));
      momentum=exp(I_c*kvector*i);
      Fermi1even.setElement(gauss*momentum,Indices(3,6,1));
      Fermi1even.setElement(gauss*momentum,Indices(4,5,3));
      //Contract and set entries
      res = reshape(Fermi1even,Indices(D_bond*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-3,Local_op,true);
      
      res = reshape(Fermi2even,Indices(D_bond*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-2,Local_op,true);
      
      res = reshape(Gauge,Indices(D_bond*D_bond,5))*Gauge_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-1,Local_op,true);
    }
    else
    {
      //Set site dependend entries for odd site 
      gauss = exp(-1.0/2.0/width_left/width_left*pow(i-pos_left,2.0));
      momentum=exp(-I_c*kvector*i);
      Fermi2odd.setElement(gauss*momentum,Indices(0,2,2));
      Fermi2odd.setElement(gauss*momentum,Indices(1,1,3));
      //Contract and set entries
      res = reshape(Fermi1odd,Indices(D_bond*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-3,Local_op,true);
      
      res = reshape(Fermi2odd,Indices(D_bond*D_bond,4))*Fermionic_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_fermi,this->d_fermi));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-2,Local_op,true);
      
      res = reshape(Gauge,Indices(D_bond*D_bond,5))*Gauge_operators;
      res.reshape(Indices(D_bond,D_bond,this->d_link,this->d_link));
      Local_op=new Operator(permute(res,Indices(3,1,4,2)));
      mpo.setOp(3*i-1,Local_op,true); 
    }
  }
  
  
}



void SU2Hamiltonian::getMovingQuarks(MPO &mpo, double kvector) const
{
  mpo.initLength(this->N_total);
  
  
  //Construct the Matrices for even and odd site
  mwArray spin_1,spin_2;
  spin_1=mwArray(Indices(2,1,2,2));	spin_1.fillWithZero();
  spin_2=mwArray(Indices(2,2,2,1));	spin_2.fillWithZero();
  
  spin_1.setElement(ONE_c,Indices(0,0,0,0));	spin_1.setElement(ONE_c,Indices(1,0,1,1));
  spin_2.setElement(ONE_c,Indices(0,0,0,0));	spin_2.setElement(ONE_c,Indices(1,1,1,0));
  
  //Momentum
  complex_t momentum;
  
  bool spin_1_saved=false;
  for(int i=1; i<=this->N_fermions; i++)
  {
    
    if((i%2)==0)
    {
      //Anti-Quark site, momentum to the right
      momentum = exp(I_c*2.0*M_PIl/(this->N_fermions+1)*kvector*((double) i));      
      spin_2.setElement(momentum,Indices(0,1,0,0));
      spin_2.setElement(momentum,Indices(1,0,1,0));
    }
    else
    {
      //Quark site, momentum to the left
      momentum = exp(-I_c*2.0*M_PIl/(this->N_fermions+1)*kvector*((double) i));
      spin_2.setElement(momentum,Indices(0,1,0,0));
      spin_2.setElement(momentum,Indices(1,0,1,0));
    }
    //cout << "Spin 2: " << spin_2 << endl;
    if(spin_1_saved)
      mpo.setOp(3*i-3,&mpo.getOp(0),false);
    else
    {
      mpo.setOp(3*i-3,new Operator(spin_1),true);
      spin_1_saved=true;
    }
    mpo.setOp(3*i-2,new Operator(spin_2),true);
  }  
  
  //Put identities in the links
  for(int i=1; i<this->N_fermions; i++)
  {
    if(i==1)
    {
      //First link, I have to save the operator
      mpo.setOp(3*i-1,new Operator(this->id_link_mpo),true);
    }
    else
    {
      //Now I can just take the pointer of the first link
      mpo.setOp(3*i-1,&mpo.getOp(3*1-1),false);
    }
  }  
  //cout << "MPO: " << mpo << endl;
}

void SU2Hamiltonian::getMovingQuarks2(MPO &mpo, double kvector, int species) const
{
  mpo.initLength(this->N_total);
  
  if(species<=0 || species>2)
  {
    cout << "Error in getMovingQuarks2(), species has to be 1 or 2, instead recieved " << species << endl;
    exit(666);
  }
  int offset_op_site, offset_id_site;
  if(species==1)
  {
    offset_op_site=3;
    offset_id_site=2;
  }
  else
  {
    offset_id_site=3;
    offset_op_site=2;
  }
    
  
  //Construct the Matrices for even and odd site
  mwArray spin_matrix;
  spin_matrix=mwArray(Indices(2,1,2,1));	spin_matrix.fillWithZero();
  
  //Momentum
  complex_t momentum;
  
  //Put identities on the species  
  bool spin_id_saved=false;
  int pos_spin_id_saved;
  for(int i=1; i<=this->N_fermions; i++)
  {
    if(spin_id_saved)
      mpo.setOp(3*i-offset_id_site,&mpo.getOp(pos_spin_id_saved),false);
    else
    {
      pos_spin_id_saved=3*i-offset_id_site;
      mpo.setOp(3*i-offset_id_site,new Operator(this->id_fermi_mpo),true);
    }      
  }
  
  for(int i=1; i<=this->N_fermions; i++)
  {
    
    if((i%2)==0)
    {
      //Anti-Quark site, momentum to the right
      momentum = exp(I_c*2.0*M_PIl/(this->N_fermions+1)*kvector*((double) i));      
      spin_matrix.setElement(momentum,Indices(0,0,0,0));
      spin_matrix.setElement(ONE_c,Indices(1,0,1,0));
    }
    else
    {
      //Quark site, momentum to the left
      momentum = exp(-I_c*2.0*M_PIl/(this->N_fermions+1)*kvector*((double) i));
      spin_matrix.setElement(ONE_c,Indices(0,0,0,0));
      spin_matrix.setElement(momentum,Indices(1,0,1,0));
    }
    mpo.setOp(3*i-offset_op_site,new Operator(spin_matrix),true);
  }  
  
  //Put identities in the links
  for(int i=1; i<this->N_fermions; i++)
  {
    if(i==1)
    {
      //First link, I have to save the operator
      mpo.setOp(3*i-1,new Operator(this->id_link_mpo),true);
    }
    else
    {
      //Now I can just take the pointer of the first link
      mpo.setOp(3*i-1,&mpo.getOp(3*1-1),false);
    }
  }
}

void SU2Hamiltonian::getMovingParticles(MPS &mps,int pos_left, int pos_right, double width_left, double width_right, double k_vec)
{
  //Bond dimension
  int Dbond=3;
  
  //k-vector 
  double kvector=2.0*M_PIl/(this->N_fermions/2.0)*k_vec;
  
  vector<int> phys_dims(this->N_total,this->d_fermi);
  vector<int> bond_dims(this->N_total-1,Dbond*Dbond);
  for(int k=2; k<(this->N_total-2); k+=3)
    phys_dims[k]=this->d_link;    
  
  bond_dims[0]=Dbond;
  bond_dims[1]=Dbond;
  bond_dims[2]=Dbond;
  
  bond_dims[this->N_total-4]=Dbond;
  bond_dims[this->N_total-3]=Dbond;
  bond_dims[this->N_total-2]=Dbond;
  
  
  mps.clear();
  mps=MPS(this->N_total,bond_dims,phys_dims);
  
  //Determine width and position (if not set externally or given values are out of range)
  if(pos_left<=0 || pos_left>this->N_fermions)
    pos_left=round(this->N_fermions/4.0);
  if(pos_right<=0 || pos_left>this->N_fermions)
    pos_right=round(3.0*this->N_fermions/4.0);
  if(width_left<=0)
    width_left=2.0;
  if(width_right<=0)
    width_right=2.0;
  
  cout << "Creating superposition with " << endl
       << "pos_left:    " << pos_left << endl
       << "pos_right:   " << pos_right << endl
       << "width_left:  " << width_left << endl
       << "width_right: " << width_right << endl
       << "k_vec:       " << k_vec << endl
       << "kvector:     " << kvector << endl;
       
  //Seed the random number generator 
  srand(time(NULL));
       
  //I assume that I have an even number of fermionic sites, if this is not the case I give a warning
  if((this->N_fermions%2)!=0)
    cout << "Warning in SU2Hamiltonian::getMovingParticles(), expected an even number of sites, instead got " << this->N_fermions << endl;
  
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(2,this->d_fermi,1));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(5,this->d_link,1));	Gauge_operators.fillWithZero();  
  
  //Spin up
  Fermionic_operators.setElement(ONE_c,Indices(0,0,0));
  //Spin down
  Fermionic_operators.setElement(ONE_c,Indices(1,1,0));
  
  Fermionic_operators.reshape(Indices(2,this->d_fermi));
  
  //State e_0
  Gauge_operators.setElement(ONE_c,Indices(0,0,0));
  //State e_1
  Gauge_operators.setElement(ONE_c,Indices(1,1,0));
  //State e_2
  Gauge_operators.setElement(ONE_c,Indices(2,2,0));
  //State e_3
  Gauge_operators.setElement(ONE_c,Indices(3,3,0));
  //State e_4
  Gauge_operators.setElement(ONE_c,Indices(4,4,0));
  
  Gauge_operators.reshape(Indices(5,this->d_link)); 
  
  mwArray Fermi1First(Indices(1,Dbond,2)),Fermi1Spin1(Indices(Dbond,Dbond,2)),Fermi1Spin2(Indices(Dbond,Dbond,2)),Gauge(Indices(1,1,5)),Fermi1Last(Indices(Dbond,1,2)),Fermi2First(Indices(1,Dbond,2)),Fermi2Last(Indices(Dbond,1,2)),Fermi2Spin1(Indices(Dbond,Dbond,2)),Fermi2Spin2(Indices(Dbond,Dbond,2));
  
  Fermi1First.fillWithZero();
  Fermi1Last.fillWithZero();
  Fermi1Spin1.fillWithZero();
  Fermi1Spin2.fillWithZero();  
  Fermi2First.fillWithZero();
  Fermi2Last.fillWithZero();
  Fermi2Spin1.fillWithZero();
  Fermi2Spin2.fillWithZero();
  Gauge.fillWithZero();
  
  //Set the site independent entries of the matrices 
  Fermi1First.setElement(ONE_c,Indices(0,0,1));
  
  Fermi1Last.setElement(ONE_c,Indices(1,0,0));
  Fermi1Last.setElement(ONE_c,Indices(2,0,1));
  
  Fermi1Spin1.setElement(ONE_c,Indices(0,0,1));
  Fermi1Spin1.setElement(ONE_c,Indices(2,2,1));
  
  Fermi1Spin2.setElement(ONE_c,Indices(1,2,0));
  Fermi1Spin2.setElement(ONE_c,Indices(0,0,1));
  Fermi1Spin2.setElement(ONE_c,Indices(2,2,1));
  
  Fermi2First.setElement(ONE_c,Indices(0,0,0));
  
  Fermi2Last.setElement(ONE_c,Indices(2,0,0));
  Fermi2Last.setElement(ONE_c,Indices(1,0,1));
  
  Fermi2Spin1.setElement(ONE_c,Indices(0,0,0));
  Fermi2Spin1.setElement(ONE_c,Indices(2,2,0));
  
  Fermi2Spin2.setElement(ONE_c,Indices(0,0,0));
  Fermi2Spin2.setElement(ONE_c,Indices(2,2,0));
  Fermi2Spin2.setElement(ONE_c,Indices(1,2,1));
  
  Gauge.setElement(ONE_c,Indices(0,0,2));
  
  //Matrices to store the contraction of the finite automata matrices with the states
  mwArray GaugeMatrix,Matrix,tmp,MPSMatrix;
  
  GaugeMatrix = reshape(Gauge,Indices(1*1,5))*Gauge_operators;
  GaugeMatrix.reshape(Indices(1,1,this->d_link));
  GaugeMatrix.permute(Indices(3,1,2));
  
  complex_t phase_factor;
  
  int count=0;
  
  //Particles
  cout << endl << "Setting particles " << endl;
  for(int i=1; i<=(this->N_fermions/2); i++)
  {
    phase_factor=exp(-I_c*kvector*i);
    if(i==1)
    {
      //First site
      Fermi1First.setElement(exp(-1.0/2.0/width_left*(i-pos_left)*(i-pos_left))*phase_factor,Indices(0,1,0));
      
      Matrix = reshape(Fermi1First,Indices(1*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(1,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));    
      cout << "Setting site " << count << endl;
      mps.setA(count,Matrix);
      count++;
      
      Matrix = reshape(Fermi1Spin2,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      cout << "Setting site " << count << endl;      
      mps.setA(count,Matrix);
      //Jump to the next particle site
      count += 5;
    }
    else if(i==(this->N_fermions/2))
    {
      //Site before the last one
      Fermi1Spin1.setElement(exp(-1.0/2.0/width_left*(i-pos_left)*(i-pos_left))*phase_factor,Indices(0,1,0));
      
      Matrix = reshape(Fermi1Spin1,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      //Now compute the Kronecker product
      MPSMatrix = mwArray(Indices(this->d_fermi,Dbond*Dbond,Dbond*Dbond));
      for(int k=0; k<this->d_fermi; k++)
      {
	tmp = Matrix.subArray(Indices(k,-1,-1));
	tmp = kron(tmp,identityMatrix(Dbond));
	for(int ind1=0; ind1<Dbond*Dbond; ind1++)
	  for(int ind2=0; ind2<Dbond*Dbond; ind2++)
	    MPSMatrix.setElement(tmp.getElement(Indices(ind1,ind2)),Indices(k,ind1,ind2));
      }
      cout << "Setting site " << count << endl;
      mps.setA(count,MPSMatrix);
      count++;
      //Last site (antiparticle site as I assume N_fermions to be even)
      Matrix = reshape(Fermi1Last,Indices(Dbond*1,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,1,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      
      //Now compute the Kronecker product
      MPSMatrix = mwArray(Indices(this->d_fermi,Dbond*Dbond,Dbond));
      for(int k=0; k<this->d_fermi; k++)
      {
	tmp = Matrix.subArray(Indices(k,-1,-1));
	tmp = kron(tmp,identityMatrix(Dbond));
	for(int ind1=0; ind1<Dbond*Dbond; ind1++)
	  for(int ind2=0; ind2<Dbond; ind2++)
	    MPSMatrix.setElement(tmp.getElement(Indices(ind1,ind2)),Indices(k,ind1,ind2));
      }
      cout << "Setting site " << count << endl;
      mps.setA(count,MPSMatrix);
    }
    else
    {
      //Odd site (particle site)
      Fermi1Spin1.setElement(exp(-1.0/2.0/width_left*(i-pos_left)*(i-pos_left))*phase_factor,Indices(0,1,0));
      
      Matrix = reshape(Fermi1Spin1,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      //Now compute the Kronecker product
      MPSMatrix = mwArray(Indices(this->d_fermi,Dbond*Dbond,Dbond*Dbond));
      for(int k=0; k<this->d_fermi; k++)
      {
	tmp = Matrix.subArray(Indices(k,-1,-1));
	tmp = kron(tmp,identityMatrix(Dbond));
	for(int ind1=0; ind1<Dbond*Dbond; ind1++)
	  for(int ind2=0; ind2<Dbond*Dbond; ind2++)
	    MPSMatrix.setElement(tmp.getElement(Indices(ind1,ind2)),Indices(k,ind1,ind2));
      }
      cout << "Setting site " << count << endl;
      mps.setA(count,MPSMatrix);
      count++;
      
      Matrix = reshape(Fermi1Spin2,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      //Now compute the Kronecker product
      MPSMatrix = mwArray(Indices(this->d_fermi,Dbond*Dbond,Dbond*Dbond));
      for(int k=0; k<this->d_fermi; k++)
      {
	tmp = Matrix.subArray(Indices(k,-1,-1));
	tmp = kron(tmp,identityMatrix(Dbond));
	for(int ind1=0; ind1<Dbond*Dbond; ind1++)
	  for(int ind2=0; ind2<Dbond*Dbond; ind2++)
	    MPSMatrix.setElement(tmp.getElement(Indices(ind1,ind2)),Indices(k,ind1,ind2));
      }
      cout << "Setting site " << count << endl;
      mps.setA(count,MPSMatrix);
      //Jump to the next particle site
      count += 5;
    }
  } 
  
  //Antiparticles
  cout << endl << "Setting antiparticles" << endl;
  count = 3;
  for(int i=1; i<=(this->N_fermions/2); i++)
  {
    phase_factor=exp(I_c*kvector*i);
    if(i==1)
    {
      //First site
      Fermi2First.setElement(exp(-1.0/2.0/width_right*(i-pos_right)*(i-pos_right))*phase_factor,Indices(0,1,1));
      
      Matrix = reshape(Fermi2First,Indices(1*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(1,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      //Now compute the Kronecker product
      MPSMatrix = mwArray(Indices(this->d_fermi,Dbond,Dbond*Dbond));
      for(int k=0; k<this->d_fermi; k++)
      {
	tmp = Matrix.subArray(Indices(k,-1,-1));
	tmp = kron(identityMatrix(Dbond),tmp);
	for(int ind1=0; ind1<Dbond; ind1++)
	  for(int ind2=0; ind2<Dbond*Dbond; ind2++)
	    MPSMatrix.setElement(tmp.getElement(Indices(ind1,ind2)),Indices(k,ind1,ind2));
      }
      cout << "Setting site " << count << endl;
      mps.setA(count,MPSMatrix);
      count++;
      
      Matrix = reshape(Fermi2Spin2,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      //Now compute the Kronecker product
      MPSMatrix = mwArray(Indices(this->d_fermi,Dbond*Dbond,Dbond*Dbond));
      for(int k=0; k<this->d_fermi; k++)
      {
	tmp = Matrix.subArray(Indices(k,-1,-1));
	tmp = kron(identityMatrix(Dbond),tmp);
	for(int ind1=0; ind1<Dbond*Dbond; ind1++)
	  for(int ind2=0; ind2<Dbond*Dbond; ind2++)
	    MPSMatrix.setElement(tmp.getElement(Indices(ind1,ind2)),Indices(k,ind1,ind2));
      }
      cout << "Setting site " << count << endl;
      mps.setA(count,MPSMatrix);
      //Jump to the next antiparticle site
      count+=5;
    }
    else if(i==(this->N_fermions/2))
    {
      //Last site
      Fermi2Spin1.setElement(exp(-1.0/2.0/width_right*(i-pos_right)*(i-pos_right))*phase_factor,Indices(0,1,0));
      
      Matrix = reshape(Fermi2Spin1,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      cout << "Setting site " << count << endl;
      mps.setA(count,Matrix);
      count++;
      
      Matrix = reshape(Fermi2Last,Indices(Dbond*1,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,1,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      cout << "Setting site " << count << endl;
      mps.setA(count,Matrix);      
    }
    else
    {
      //Sites in between
      Fermi2Spin1.setElement(exp(-1.0/2.0/width_right*(i-pos_right)*(i-pos_right))*phase_factor,Indices(0,1,0));
      
      Matrix = reshape(Fermi2Spin1,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      //Now compute the Kronecker product
      MPSMatrix = mwArray(Indices(this->d_fermi,Dbond*Dbond,Dbond*Dbond));
      for(int k=0; k<this->d_fermi; k++)
      {
	tmp = Matrix.subArray(Indices(k,-1,-1));
	tmp = kron(identityMatrix(Dbond),tmp);
	for(int ind1=0; ind1<Dbond*Dbond; ind1++)
	  for(int ind2=0; ind2<Dbond*Dbond; ind2++)
	    MPSMatrix.setElement(tmp.getElement(Indices(ind1,ind2)),Indices(k,ind1,ind2));
      }
      cout << "Setting site " << count << endl;
      mps.setA(count,MPSMatrix);
      count++;
      
      Matrix = reshape(Fermi2Spin2,Indices(Dbond*Dbond,2))*Fermionic_operators;
      Matrix.reshape(Indices(Dbond,Dbond,this->d_fermi));
      Matrix.permute(Indices(3,1,2));
      //Now compute the Kronecker product
      MPSMatrix = mwArray(Indices(this->d_fermi,Dbond*Dbond,Dbond*Dbond));
      for(int k=0; k<this->d_fermi; k++)
      {
	tmp = Matrix.subArray(Indices(k,-1,-1));
	tmp = kron(identityMatrix(Dbond),tmp);
	for(int ind1=0; ind1<Dbond*Dbond; ind1++)
	  for(int ind2=0; ind2<Dbond*Dbond; ind2++)
	    MPSMatrix.setElement(tmp.getElement(Indices(ind1,ind2)),Indices(k,ind1,ind2));
      }
      cout << "Setting site " << count << endl;
      mps.setA(count,MPSMatrix);
      //Jump to the next antiparticle site
      count+=5;
    }
  }  
  //Gauge Links
  cout << endl <<"Setting Gauge Links " <<  endl;
  for(int i=2; i<=(this->N_total-3); i+=3)
  {
    if(i==2)
    {
      MPSMatrix = mwArray(Indices(this->d_link,Dbond,Dbond));
      for(int k=0; k<this->d_link; k++)
      {
	tmp = GaugeMatrix.subArray(Indices(k,-1,-1));
	tmp = kron(identityMatrix(Dbond),tmp);
	for(int ind1=0; ind1<Dbond; ind1++)
	  for(int ind2=0; ind2<Dbond; ind2++)
	    MPSMatrix.setElement(tmp.getElement(Indices(ind1,ind2)),Indices(k,ind1,ind2));
      }
      cout << "Setting site " << i << endl;
      mps.setA(i,MPSMatrix);
    }
    else if(i==(this->N_total-3))
    {
      MPSMatrix = mwArray(Indices(this->d_link,Dbond,Dbond));
      for(int k=0; k<this->d_link; k++)
      {
	tmp = GaugeMatrix.subArray(Indices(k,-1,-1));
	tmp = kron(tmp,identityMatrix(Dbond));
	for(int ind1=0; ind1<Dbond; ind1++)
	  for(int ind2=0; ind2<Dbond; ind2++)
	    MPSMatrix.setElement(tmp.getElement(Indices(ind1,ind2)),Indices(k,ind1,ind2));
      }
      cout << "Setting site " << i << endl;
      mps.setA(i,MPSMatrix);
    }
    else
    {
      MPSMatrix = mwArray(Indices(this->d_link,Dbond*Dbond,Dbond*Dbond));
      for(int k=0; k<this->d_link; k++)
      {
	tmp = GaugeMatrix.subArray(Indices(k,-1,-1));
	tmp = kron(identityMatrix(Dbond),kron(tmp,identityMatrix(Dbond)));
	for(int ind1=0; ind1<Dbond*Dbond; ind1++)
	  for(int ind2=0; ind2<Dbond*Dbond; ind2++)
	    MPSMatrix.setElement(tmp.getElement(Indices(ind1,ind2)),Indices(k,ind1,ind2));
      }
      cout << "Setting site " << i << endl;
      mps.setA(i,MPSMatrix);      
    } 
  }
  mps.gaugeCond('R',true);
}


void SU2Hamiltonian::getHeavyString(MPS &mps, unsigned int length) const
{
  int pos = round(this->N_fermions/2)-round(length/2);
  //Some error checking
  if(pos>this->N_total || pos<1 || (pos+length)>this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getHeavyString, desired string of length " << length<< " starting on position " << pos << " is outside the the chain of " << this->N_fermions << " spins" << endl;
    exit(666);
  }
  
  cout << "Creating heavy strong coupling string with " << endl
       << "pos:    " << pos << endl
       << "length: " << length << endl;  
  vector <string_t> config;
  config.push_back((string_t) {pos,length});
  
  this->constructStringMPS(mps,config);
}

void SU2Hamiltonian::constructStringMPS(MPS &mps, vector<string_t> string_config) const
{  
  if(string_config.empty())
  {
    cout << "Warning in SU2Hamiltonian::getStringMPS(), no string configurations given, will return a strong coupling state" << endl;
    this->constructInitialMPS(mps);
    return;
  }
  
  //Bond dimension
  int Dbond=2; 
  
  vector<int> phys_dims(this->N_total,this->d_fermi);
  vector<int> bond_dims(this->N_total-1,1);
  for(int k=2; k<(this->N_total-2); k+=3)
    phys_dims[k]=this->d_link;    
  
   int pos,length;
//   for (int k=0; k<string_config.size(); k++)
//   {
//     pos = string_config[k].pos;
//     length = string_config[k].leng;
//     for (int i=pos; i<(pos+length-1); i++)
//     {
//       bond_dims[3*i-1]=Dbond;
//       if((3*(i+1)-2)<(3*(pos+length)-3))
//       {
// 	bond_dims[3*(i+1)-3]=Dbond;
// 	bond_dims[3*(i+1)-2]=Dbond;
//       }
//     }   
//   }
  
  mps.clear();
  mps=MPS(this->N_total,bond_dims,phys_dims);
  //cout << "MPS: " << mps << endl;

  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(2,this->d_fermi,1));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(5,this->d_link,1));	Gauge_operators.fillWithZero();  
  
  //Spin up
  Fermionic_operators.setElement(ONE_c,Indices(0,0,0));
  //Spin down
  Fermionic_operators.setElement(ONE_c,Indices(1,1,0));
  
  Fermionic_operators.reshape(Indices(2,this->d_fermi));
  
  //State e_0
  Gauge_operators.setElement(ONE_c,Indices(0,0,0));
  //State e_1
  Gauge_operators.setElement(ONE_c,Indices(1,1,0));
  //State e_2
  Gauge_operators.setElement(ONE_c,Indices(2,2,0));
  //State e_3
  Gauge_operators.setElement(ONE_c,Indices(3,3,0));
  //State e_4
  Gauge_operators.setElement(ONE_c,Indices(4,4,0));
  
  Gauge_operators.reshape(Indices(5,this->d_link)); 
  
  //Matrices outside the string region (simply bond dimension one)
  mwArray FermiOddout,FermiEvenout,Gaugeout;
  FermiOddout = mwArray(Indices(1,1,2));	FermiOddout.fillWithZero();
  FermiEvenout = mwArray(Indices(1,1,2));	FermiEvenout.fillWithZero();
  Gaugeout = mwArray(Indices(1,1,5));		Gaugeout.fillWithZero();
  
  FermiOddout.setElement(ONE_c,Indices(0,0,0));  
  FermiEvenout.setElement(ONE_c,Indices(0,0,1));  
  Gaugeout.setElement(ONE_c,Indices(0,0,2));
  
  //Matrices inside the string region
  mwArray FermiOdd,FermiEven,GaugeOdd,GaugeFirstOdd,GaugeLastOdd,GaugeEven,GaugeFirstEven,GaugeLastEven,GaugeSingleOdd,GaugeSingleEven;
  FermiOdd = mwArray(Indices(Dbond,Dbond,2));	FermiOdd.fillWithZero();;
  FermiEven = mwArray(Indices(Dbond,Dbond,2));	FermiEven.fillWithZero();;
  
  GaugeOdd = mwArray(Indices(Dbond,Dbond,5));		GaugeOdd.fillWithZero();
  GaugeFirstOdd = mwArray(Indices(1,Dbond,5));		GaugeFirstOdd.fillWithZero();
  GaugeLastOdd = mwArray(Indices(Dbond,1,5));		GaugeLastOdd.fillWithZero();
  GaugeSingleOdd = mwArray(Indices(1,1,5));		GaugeSingleOdd.fillWithZero();
  GaugeSingleEven = mwArray(Indices(1,1,5));		GaugeSingleEven.fillWithZero();
  
  GaugeEven = mwArray(Indices(Dbond,Dbond,5));		GaugeEven.fillWithZero();
  GaugeFirstEven = mwArray(Indices(1,Dbond,5));		GaugeFirstEven.fillWithZero();
  GaugeLastEven = mwArray(Indices(Dbond,1,5));		GaugeLastEven.fillWithZero();
  
  GaugeFirstOdd.setElement(ONE_c,Indices(0,0,0));
  GaugeFirstOdd.setElement(I_c,Indices(0,1,1));
  GaugeFirstOdd.setElement(-I_c,Indices(0,0,3));
  GaugeFirstOdd.setElement(ONE_c,Indices(0,1,4));
  
  GaugeOdd.setElement(ONE_c,Indices(0,0,0));
  GaugeOdd.setElement(I_c,Indices(0,1,1));
  GaugeOdd.setElement(-I_c,Indices(1,0,3));
  GaugeOdd.setElement(ONE_c,Indices(1,1,4));
  
  GaugeLastOdd.setElement(ONE_c,Indices(0,0,0));
  GaugeLastOdd.setElement(I_c,Indices(0,0,1));
  GaugeLastOdd.setElement(-I_c,Indices(1,0,3));
  GaugeLastOdd.setElement(ONE_c,Indices(1,0,4));
  
  GaugeFirstEven.setElement(ONE_c,Indices(0,0,0));
  GaugeFirstEven.setElement(-I_c,Indices(0,1,1));
  GaugeFirstEven.setElement(I_c,Indices(0,0,3));
  GaugeFirstEven.setElement(ONE_c,Indices(0,1,4));
  
  GaugeEven.setElement(ONE_c,Indices(0,0,0));
  GaugeEven.setElement(-I_c,Indices(0,1,1));
  GaugeEven.setElement(I_c,Indices(1,0,3));
  GaugeEven.setElement(ONE_c,Indices(1,1,4));
  
  GaugeLastEven.setElement(ONE_c,Indices(0,0,0));
  GaugeLastEven.setElement(-I_c,Indices(0,0,1));
  GaugeLastEven.setElement(I_c,Indices(1,0,3));
  GaugeLastEven.setElement(ONE_c,Indices(1,0,4));
  
  GaugeSingleOdd.setElement(ONE_c,Indices(0,0,0));
  GaugeSingleOdd.setElement(I_c,Indices(0,0,1));
  GaugeSingleOdd.setElement(-I_c,Indices(0,0,3));
  GaugeSingleOdd.setElement(ONE_c,Indices(0,0,4));
  
  GaugeSingleEven.setElement(ONE_c,Indices(0,0,0));
  GaugeSingleEven.setElement(-I_c,Indices(0,0,1));
  GaugeSingleEven.setElement(I_c,Indices(0,0,3));
  GaugeSingleEven.setElement(ONE_c,Indices(0,0,4));
  
  FermiOdd.setElement(ONE_c,Indices(0,0,0));
  FermiOdd.setElement(ONE_c,Indices(1,1,0));
  
  FermiEven.setElement(ONE_c,Indices(0,0,1));
  FermiEven.setElement(ONE_c,Indices(1,1,1));
  
  //Provide the matrices
  mwArray MFermiOddout, MFermiEvenout, MGaugeout, MFermiEven, MFermiOdd, MGaugeOdd, MGaugeLastOdd, MGaugeFirstOdd, MGaugeEven, MGaugeLastEven, MGaugeFirstEven, MGaugeSingleOdd, MGaugeSingleEven;
  
  MFermiOddout = reshape(FermiOddout,Indices(1,2))*Fermionic_operators;
  MFermiOddout.reshape(Indices(1,1,this->d_fermi));
  MFermiOddout.permute(Indices(3,1,2));
  
  MFermiEvenout = reshape(FermiEvenout,Indices(1,2))*Fermionic_operators;
  MFermiEvenout.reshape(Indices(1,1,this->d_fermi));
  MFermiEvenout.permute(Indices(3,1,2));
  
  MGaugeout = reshape(Gaugeout,Indices(1,5))*Gauge_operators;
  MGaugeout.reshape(Indices(1,1,this->d_link));
  MGaugeout.permute(Indices(3,1,2));
  
  MFermiOdd = reshape(FermiOdd,Indices(Dbond*Dbond,2))*Fermionic_operators;
  MFermiOdd.reshape(Indices(Dbond,Dbond,this->d_fermi));
  MFermiOdd.permute(Indices(3,1,2));
  
  MFermiEven = reshape(FermiEven,Indices(Dbond*Dbond,2))*Fermionic_operators;
  MFermiEven.reshape(Indices(Dbond,Dbond,this->d_fermi));
  MFermiEven.permute(Indices(3,1,2));
  
  MGaugeOdd = reshape(GaugeOdd,Indices(Dbond*Dbond,5))*Gauge_operators;
  MGaugeOdd.reshape(Indices(Dbond,Dbond,this->d_link));
  MGaugeOdd.permute(Indices(3,1,2));
  
  MGaugeFirstOdd = reshape(GaugeFirstOdd,Indices(1*Dbond,5))*Gauge_operators;
  MGaugeFirstOdd.reshape(Indices(1,Dbond,this->d_link));
  MGaugeFirstOdd.permute(Indices(3,1,2));
  
  MGaugeLastOdd = reshape(GaugeLastOdd,Indices(Dbond*1,5))*Gauge_operators;
  MGaugeLastOdd.reshape(Indices(Dbond,1,this->d_link));
  MGaugeLastOdd.permute(Indices(3,1,2));
  
  MGaugeEven = reshape(GaugeEven,Indices(Dbond*Dbond,5))*Gauge_operators;
  MGaugeEven.reshape(Indices(Dbond,Dbond,this->d_link));
  MGaugeEven.permute(Indices(3,1,2));
  
  MGaugeFirstEven = reshape(GaugeFirstEven,Indices(1*Dbond,5))*Gauge_operators;
  MGaugeFirstEven.reshape(Indices(1,Dbond,this->d_link));
  MGaugeFirstEven.permute(Indices(3,1,2));
  
  MGaugeLastEven = reshape(GaugeLastEven,Indices(Dbond*1,5))*Gauge_operators;
  MGaugeLastEven.reshape(Indices(Dbond,1,this->d_link));
  MGaugeLastEven.permute(Indices(3,1,2));
  
  MGaugeSingleOdd = reshape(GaugeSingleOdd,Indices(1,5))*Gauge_operators;
  MGaugeSingleOdd.reshape(Indices(1,1,this->d_link));
  MGaugeSingleOdd.permute(Indices(3,1,2)); 
  
  MGaugeSingleEven = reshape(GaugeSingleEven,Indices(1,5))*Gauge_operators;
  MGaugeSingleEven.reshape(Indices(1,1,this->d_link));
  MGaugeSingleEven.permute(Indices(3,1,2));  
  
  /*cout << "MFermiOddout: " << MFermiOddout << endl;
  cout << "MFermiEvenout: " << MFermiEvenout << endl;
  cout << "MGaugeout: " << MGaugeout << endl;
  
  cout << "MFermiOdd: " << MFermiOdd << endl;
  cout << "MFermiEven: " << MFermiEven << endl;
  cout << "MGauge: " << MGauge << endl;
  cout << "MGaugeFirst: " << MGaugeFirst << endl;
  cout << "MGaugeLast: " << MGaugeLast << endl;*/
  
  //Fill completely with matrices for non string region
  for(int i=1; i<=this->N_fermions; i++)
  {
    if((i%2)==0)
    {
      mps.setA(3*i-3,MFermiEvenout);
      mps.setA(3*i-2,MFermiEvenout);
    }
    else
    {
      mps.setA(3*i-3,MFermiOddout);
      mps.setA(3*i-2,MFermiOddout);
    }
    if(i<this->N_fermions)
      mps.setA(3*i-1,MGaugeout); 
  } 
  
  //Replace the string parts
  for(int k=0; k<string_config.size(); k++)
  {
    pos = string_config[k].pos;
    length = string_config[k].leng;
    
    if((pos%2)==0 && length>1)
      mps.replaceSite(3*pos-1,MGaugeFirstEven,0);
    else if((pos%2)!=0 && length>1)
      mps.replaceSite(3*pos-1,MGaugeFirstOdd,0);
    else if((pos%2)==0)
      mps.replaceSite(3*pos-1,MGaugeSingleEven,0);
    else
      mps.replaceSite(3*pos-1,MGaugeSingleOdd,0);
      
    for (int i=pos+1; i<pos+length; i++)
    {
      if((i%2)==0)      
      {
	mps.replaceSite(3*i-3,MFermiEven,0);
	mps.replaceSite(3*i-2,MFermiEven,0);
      }
      else
      {
	mps.replaceSite(3*i-3,MFermiOdd,0);
	mps.replaceSite(3*i-2,MFermiOdd,0);
      }
      if(i==(pos+length-1) && ((pos%2)==0))
	mps.replaceSite(3*i-1,MGaugeLastEven,0);
      else if(i==(pos+length-1) && ((pos%2)!=0))
	mps.replaceSite(3*i-1,MGaugeLastOdd,0);
      else if((pos%2)==0)
	mps.replaceSite(3*i-1,MGaugeEven,0);
      else
	mps.replaceSite(3*i-1,MGaugeOdd,0);
    }
  } 
  mps.gaugeCond('R',true);
}

void SU2Hamiltonian::getHeavyString(MPO &mpo, unsigned int length, unsigned int pos, bool conjugate) const
{
  //Bond dimension
  int Dbond=2;
  bool pos_is_odd;
    
  //Determine the starting position of the string
  if(pos<1)
    pos = round(this->N_fermions/2)-round(length/2);
  //Some error checking
  if(pos>this->N_total || pos<1 || (pos+length)>this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getHeavyString(), desired string of length " << length<< " starting on position " << pos << " is outside the the chain of " << this->N_fermions << " spins" << endl;
    exit(666);
  }
  //Should also work for even length
  /*if((length%2)==0 || length==1)
  {
    cout << "Error in SU2Hamiltonian::getHeavyString(), only odd values greater than one for length are supported, received " << length << endl;
    exit(666);
  }*/
  if(length==1)
  {
    cout << "Error in SU2Hamiltonian::getHeavyString(), only values greater than one for length are supported, received " << length << endl;
    exit(666);
  }
  //Initialize MPO
  mpo.initLength(this->N_total);
  //Determine if startposition is odd 
  pos_is_odd=(pos%2);
  
  //If I want the conjugate operators, I simply pretend that the startposition is odd (even) when I start at an even (odd) site
  if(conjugate)
    pos_is_odd = (!pos_is_odd);
  
  //Set up the operators
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(1,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(4,this->d_link,this->d_link));	Gauge_operators.fillWithZero();  
  
  for(int i=0; i<this->d_fermi; i++)
    for(int j=0; j<this->d_fermi; j++)
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
  
  for(int i=0; i<this->d_link; i++)
  {
    for(int j=0; j<this->d_link; j++)
    {
      if(pos_is_odd)
      {
	Gauge_operators.setElement(this->U00ad.getElement(Indices(i,j)),Indices(0,i,j));
	Gauge_operators.setElement(-1.0*this->U01ad.getElement(Indices(i,j))*I_c,Indices(1,i,j));
	Gauge_operators.setElement(this->U10ad.getElement(Indices(i,j))*I_c,Indices(2,i,j));
	Gauge_operators.setElement(this->U11ad.getElement(Indices(i,j)),Indices(3,i,j));
      }
      else
      {
	Gauge_operators.setElement(this->U00.getElement(Indices(i,j)),Indices(0,i,j));
	Gauge_operators.setElement(this->U01.getElement(Indices(i,j))*I_c,Indices(1,i,j));
	Gauge_operators.setElement(-1.0*this->U10.getElement(Indices(i,j))*I_c,Indices(2,i,j));
	Gauge_operators.setElement(this->U11.getElement(Indices(i,j)),Indices(3,i,j));
      }
    }
  }
    
  Fermionic_operators.reshape(Indices(1,this->d_fermi*this->d_fermi));
  Gauge_operators.reshape(Indices(4,this->d_link*this->d_link));  
  
  //Matrices
  mwArray GaugeFirst,GaugeEnd,Gauge,Fermi;
  
  GaugeFirst=mwArray(Indices(1,Dbond,4));	GaugeFirst.fillWithZero();
  Gauge=mwArray(Indices(Dbond,Dbond,4));	Gauge.fillWithZero();
  GaugeEnd=mwArray(Indices(Dbond,1,4));		GaugeEnd.fillWithZero();
  Fermi=mwArray(Indices(Dbond,Dbond,1));	Fermi.fillWithZero();
  
  Fermi.setElement(ONE_c,Indices(0,0,0));
  Fermi.setElement(ONE_c,Indices(1,1,0));
  
  GaugeFirst.setElement(ONE_c,Indices(0,0,0));
  GaugeFirst.setElement(ONE_c,Indices(0,1,1));
  GaugeFirst.setElement(ONE_c,Indices(0,0,2));
  GaugeFirst.setElement(ONE_c,Indices(0,1,3));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(0,1,1));
  Gauge.setElement(ONE_c,Indices(1,0,2));
  Gauge.setElement(ONE_c,Indices(1,1,3));
  
  GaugeEnd.setElement(ONE_c,Indices(0,0,0));
  GaugeEnd.setElement(ONE_c,Indices(0,0,1));
  GaugeEnd.setElement(ONE_c,Indices(1,0,2));
  GaugeEnd.setElement(ONE_c,Indices(1,0,3));
  
  Fermi.reshape(Indices(Dbond*Dbond,1));
  GaugeFirst.reshape(Indices(1*Dbond,4));
  Gauge.reshape(Indices(Dbond*Dbond,4));
  GaugeEnd.reshape(Indices(Dbond*1,4));
  
  //Contract and reshape accordingly
  Fermi=reshape(Fermi*Fermionic_operators,Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));		Fermi.permute(Indices(3,1,4,2));
  GaugeFirst=reshape(GaugeFirst*Gauge_operators,Indices(1,Dbond,this->d_link,this->d_link));	GaugeFirst.permute(Indices(3,1,4,2));
  Gauge=reshape(Gauge*Gauge_operators,Indices(Dbond,Dbond,this->d_link,this->d_link));			Gauge.permute(Indices(3,1,4,2));
  GaugeEnd=reshape(GaugeEnd*Gauge_operators,Indices(Dbond,1,this->d_link,this->d_link));		GaugeEnd.permute(Indices(3,1,4,2));
  
  //Save start and end operators
  mpo.setOp(3*pos-1, new Operator(GaugeFirst),true);
  mpo.setOp(3*(pos+length-1)-1,new Operator(GaugeEnd),true);
  
  bool fermi_saved=false,link_saved=false;
  int pos_fermi_saved,pos_link_saved;
  //Remaining sites in between
  for(int i=(pos+1); i<=(pos+length-1); i++)
  {
    if(fermi_saved)
    {
      mpo.setOp(3*i-3,&mpo.getOp(pos_fermi_saved),false);
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_saved),false);
    }
    else
    {
      mpo.setOp(3*i-3,new Operator(Fermi),true);
      pos_fermi_saved=3*i-3;
      fermi_saved=true;      
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_saved),false);
    }
    //Can I still put the gauge operator?
    if(i<(pos+length-1))
    {
      if(link_saved)
	mpo.setOp(3*i-1,&mpo.getOp(pos_link_saved),false);
      else
      {
	mpo.setOp(3*i-1,new Operator(Gauge),true);
	pos_link_saved=3*i-1;
	link_saved=true;
      }
    } 
  }

  //Remaining identities
  fermi_saved=false;
  link_saved=false;
  for(int i=1; i<=pos; i++)
  {
    if(fermi_saved)
      {
	mpo.setOp(3*i-3,&mpo.getOp(pos_fermi_saved),false);
	mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_saved),false);
      }
      else
      {
	mpo.setOp(3*i-3,new Operator(this->id_fermi_mpo),true);
	pos_fermi_saved=3*i-3;
	fermi_saved=true;      
	mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_saved),false);
    }
    //Can I still put the gauge operator?
    if(i<pos)
    {
      if(link_saved)
	mpo.setOp(3*i-1,&mpo.getOp(pos_link_saved),false);
      else
      {
	mpo.setOp(3*i-1,new Operator(this->id_link_mpo),true);
	pos_link_saved=3*i-1;
	link_saved=true;
      }
    } 
  }

  for(int i=(pos+length); i<=this->N_fermions; i++)
  {
    if(fermi_saved)
    {
      mpo.setOp(3*i-3,&mpo.getOp(pos_fermi_saved),false);
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_saved),false);
    }
    else
    {
      mpo.setOp(3*i-3,new Operator(this->id_fermi_mpo),true);
      pos_fermi_saved=3*i-3;
      fermi_saved=true;      
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_saved),false);
    }
    //Can I still put the gauge operator?
    if(i<this->N_fermions)
    {
      if(link_saved)
	mpo.setOp(3*i-1,&mpo.getOp(pos_link_saved),false);
      else
      {
	mpo.setOp(3*i-1,new Operator(this->id_link_mpo),true);
	pos_link_saved=3*i-1;
	link_saved=true;
      }
    } 
  }
} 

void SU2Hamiltonian::createSingleCharge(MPO &mpo, int n, int spin) const
{
  //Determine the site index in terms of an MPO index
  int index;
  if(spin==1)
    index = 3*n-3;
  else
    index = 3*n-2;
  //Determine whether I use sigma_plus or sigma_minus and generate the MPO
  if((n%2)==0)
    this->getSingleBodyMPO(mpo,index,spOp);
  else
    this->getSingleBodyMPO(mpo,index,smOp);    
}

void SU2Hamiltonian::getCountStringsMPO(MPO &mpo, int length, int pos) const
{
   //Some error checking
  if(pos>this->N_total || pos<1 || (pos+length)>this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getCountStringsMPO(), desired string of length " << length<< " starting on position " << pos << " is outside the the chain of " << this->N_fermions << " spins" << endl;
    exit(666);
  }
  //Intialize the MPO
  mpo.initLength(this->N_total);
  
  //MPO matrices
  mwArray fermi1antiparallel,fermi2antiparallel,fermi1parallel,fermi2parallel;
  
  fermi1antiparallel=mwArray(Indices(this->d_fermi,1,this->d_fermi,2));	fermi1antiparallel.fillWithZero();
  fermi2antiparallel=mwArray(Indices(this->d_fermi,2,this->d_fermi,1));	fermi2antiparallel.fillWithZero();
  fermi1parallel=mwArray(Indices(this->d_fermi,1,this->d_fermi,2));	fermi1parallel.fillWithZero();
  fermi2parallel=mwArray(Indices(this->d_fermi,2,this->d_fermi,1));	fermi2parallel.fillWithZero();
  
  //Only conserve parallel configurations, project out the antiparallel ones
  fermi1parallel.setElement(ONE_c,Indices(0,0,0,0));	fermi1parallel.setElement(ONE_c,Indices(1,0,1,1));
  fermi2parallel.setElement(ONE_c,Indices(0,0,0,0));	fermi2parallel.setElement(ONE_c,Indices(1,1,1,0));
  
  //Only conserve the antiparallel configurations, project out the parallel ones
  fermi1antiparallel.setElement(ONE_c,Indices(0,0,0,0));	fermi1antiparallel.setElement(ONE_c,Indices(1,0,1,1));
  fermi2antiparallel.setElement(ONE_c,Indices(0,1,0,0));	fermi2antiparallel.setElement(ONE_c,Indices(1,0,1,0));  
  
  //Set the operators
  mpo.setOp(3*pos-3,new Operator(fermi1antiparallel),true);
  mpo.setOp(3*pos-2,new Operator(fermi2antiparallel),true);
  mpo.setOp(3*(pos+length)-3,&mpo.getOp(3*pos-3),false);
  mpo.setOp(3*(pos+length)-2,&mpo.getOp(3*pos-2),false);  
  
  //Operators in between
  bool saved=false;
  int pos_saved;
  for(int i=(pos+1); i<(pos+length); i++)
  {
    if(saved)
    {
      mpo.setOp(3*i-3,&mpo.getOp(3*pos_saved-3),false);
      mpo.setOp(3*i-2,&mpo.getOp(3*pos_saved-2),false);
    }
    else
    {
      mpo.setOp(3*i-3,new Operator(fermi1parallel),true);
      mpo.setOp(3*i-2,new Operator(fermi2parallel),true);
      pos_saved=i;
      saved=true;
    }
  }
  //Set the gauge identites
  saved=false;
  for(int i=2; i<this->N_total; i+=3)
  {
    if(saved)
      mpo.setOp(i,&mpo.getOp(pos_saved),false);
    else
    {
      mpo.setOp(i,new Operator(this->id_link_mpo),true);
      pos_saved=i;
      saved=true;
    }
  }
  //Set the remaining fermionic identities
  saved=false;
  for(int i=1; i<=this->N_fermions; i++)
  {
    //Check if I'm outside the string area
    if(i<pos || i>(pos+length))
    {
      if(saved)
      {
	mpo.setOp(3*i-3,&mpo.getOp(pos_saved),false);
	mpo.setOp(3*i-2,&mpo.getOp(pos_saved),false);
      }
      else
      {
	mpo.setOp(3*i-3,new Operator(this->id_fermi_mpo),true);
	pos_saved=3*i-3;
	saved=true;
	mpo.setOp(3*i-2,&mpo.getOp(pos_saved),false);
      }
    }
  }  
}


void SU2Hamiltonian::getChargeCorrelatorMPO(MPO &mpo, int length, int pos, bool projector_on) const
{
  //Some error checking
  if(pos>this->N_total || pos<1 || (pos+length)>this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getChargeCorrelatorMPO(), desired correlator of length " << length<< " starting on position " << pos << " is outside the the chain of " << this->N_fermions << " spins" << endl;
    exit(666);
  }
  //Intialize the MPO
  mpo.initLength(this->N_total);
  
  int Dbond=3;
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(2,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(1,this->d_link,this->d_link));	Gauge_operators.fillWithZero();  
  
  for(int i=0; i<this->d_fermi; i++)
  {
    for(int j=0; j<this->d_fermi; j++)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma z
      Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }  
  for(int i=0; i<this->d_link; i++)
    for(int j=0; j<this->d_link; j++)
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
    
  Fermionic_operators.reshape(Indices(2,this->d_fermi*this->d_fermi));
  Gauge_operators.reshape(Indices(1,this->d_link*this->d_link));  
  
  mwArray Fermi1left,Fermi2left,Fermi1right,Fermi2right,Gauge,FermiInside1,FermiInside2;
  
  Fermi1left=mwArray(Indices(1,Dbond,2));		Fermi1left.fillWithZero();
  Fermi2left=mwArray(Indices(Dbond,Dbond,2));		Fermi2left.fillWithZero();
  Fermi1right=mwArray(Indices(Dbond,Dbond,2));		Fermi1right.fillWithZero();
  Fermi2right=mwArray(Indices(Dbond,1,2));		Fermi2right.fillWithZero();
  if(projector_on)
  {
    FermiInside1=mwArray(Indices(2,Dbond,2,2*Dbond));	FermiInside1.fillWithZero();
    FermiInside2=mwArray(Indices(2,2*Dbond,2,Dbond));	FermiInside2.fillWithZero();
  }
  else
  {
    FermiInside1=mwArray(Indices(Dbond,Dbond,2));	FermiInside1.fillWithZero();
    FermiInside2=mwArray(Indices(Dbond,Dbond,2));	FermiInside2.fillWithZero();
  }
  Gauge=mwArray(Indices(Dbond,Dbond,1));		Gauge.fillWithZero();  
  
  Fermi1left.setElement(ONE_c,Indices(0,0,0));
  Fermi1left.setElement(-1.0*ONE_c,Indices(0,1,1));
  
  Fermi2left.setElement(ONE_c,Indices(0,0,0));
  Fermi2left.setElement(ONE_c,Indices(1,1,1));
  
  Fermi1right.setElement(ONE_c,Indices(0,0,0));
  Fermi1right.setElement(ONE_c,Indices(1,0,0));
  Fermi1right.setElement(-1.0*ONE_c,Indices(0,2,1));
  Fermi1right.setElement(-1.0*ONE_c,Indices(1,2,1));
  
  Fermi2right.setElement(ONE_c,Indices(0,0,0));
  Fermi2right.setElement(ONE_c,Indices(2,0,1));

  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(1,1,0));  
  
  Fermi1left.reshape(Indices(1*Dbond,2));
  Fermi2left.reshape(Indices(Dbond*Dbond,2));
  Fermi1right.reshape(Indices(Dbond*Dbond,2));
  Fermi2right.reshape(Indices(Dbond*1,2));
  Gauge.reshape(Indices(Dbond*Dbond,1));
  
  if(projector_on)
  {
    for(int i=0; i<Dbond; i++)
    {
      FermiInside1.setElement(ONE_c,Indices(0,i,0,i));
      FermiInside1.setElement(ONE_c,Indices(1,i,1,i+Dbond));
      FermiInside2.setElement(ONE_c,Indices(0,i,0,i));
      FermiInside2.setElement(ONE_c,Indices(1,i+Dbond,1,i));
    }
  }
  else
  {
    FermiInside1.setElement(ONE_c,Indices(0,0,0));
    FermiInside1.setElement(ONE_c,Indices(1,1,0));
    
    FermiInside2.setElement(ONE_c,Indices(0,0,0));
    FermiInside2.setElement(ONE_c,Indices(1,1,0));
    
    FermiInside1.reshape(Indices(Dbond*Dbond,2));
    FermiInside2.reshape(Indices(Dbond*Dbond,2));
    
    FermiInside1 = FermiInside1*Fermionic_operators;
    FermiInside2 = FermiInside2*Fermionic_operators;
    
    FermiInside1.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));
    FermiInside2.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));
    
    FermiInside1.permute(Indices(3,1,4,2));
    FermiInside2.permute(Indices(3,1,4,2));
  }
  
  //Contract and set operators
  mwArray res;
  
  res=reshape(Fermi1left*Fermionic_operators,Indices(1,Dbond,this->d_fermi,this->d_fermi));
  res.permute(Indices(3,1,4,2));
  mpo.setOp(3*pos-3,new Operator(res),true);
  
  res=reshape(Fermi2left*Fermionic_operators,Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));
  res.permute(Indices(3,1,4,2));
  mpo.setOp(3*pos-2,new Operator(res),true);
  
  res=reshape(Gauge*Gauge_operators,Indices(Dbond,Dbond,this->d_link,this->d_link));
  res.permute(Indices(3,1,4,2));
  mpo.setOp(3*pos-1,new Operator(res),true);
  
  res=reshape(Fermi1right*Fermionic_operators,Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));
  res.permute(Indices(3,1,4,2));
  mpo.setOp(3*(pos+length)-3,new Operator(res),true);
  
  res=reshape(Fermi2right*Fermionic_operators,Indices(Dbond,1,this->d_fermi,this->d_fermi));
  res.permute(Indices(3,1,4,2));
  mpo.setOp(3*(pos+length)-2,new Operator(res),true);
  
  bool fermi_inside_saved=false;
  int pos_fermi_inside_saved;
  

  for(int i=(pos+1); i<(pos+length); i++)
  {
    if(fermi_inside_saved)
    {
      mpo.setOp(3*i-3,&mpo.getOp(3*pos_fermi_inside_saved-3),false);
      mpo.setOp(3*i-2,&mpo.getOp(3*pos_fermi_inside_saved-2),false);
    }
    else
    {
      mpo.setOp(3*i-3,new Operator(FermiInside1),true);
      fermi_inside_saved=true;
      pos_fermi_inside_saved=i;
      mpo.setOp(3*i-2,new Operator(FermiInside2),true);
    }
    mpo.setOp(3*i-1,&mpo.getOp(3*pos-1),false);      
  }
  //Identities to the left
  bool fermi_id_saved=false,link_id_saved=false;
  int pos_fermi_id_saved, pos_link_id_saved;

  for(int i=1; i<pos; i++)
  {
    if(fermi_id_saved)
    {
      mpo.setOp(3*i-3,&mpo.getOp(pos_fermi_id_saved),false);
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }
    else
    {
      mpo.setOp(3*i-3,new Operator(this->id_fermi_mpo),true);
      fermi_id_saved=true;
      pos_fermi_id_saved=3*i-3;
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }
    if(link_id_saved)
      mpo.setOp(3*i-1,&mpo.getOp(pos_link_id_saved),false);
    else
    {
      mpo.setOp(3*i-1,new Operator(this->id_link_mpo),true);
      link_id_saved=true;
      pos_link_id_saved=3*i-1;
    }    
  }

  for(int i=(pos+length+1); i<=this->N_fermions; i++)
  {
    if(fermi_id_saved)
    {
      mpo.setOp(3*i-3,&mpo.getOp(pos_fermi_id_saved),false);
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }
    else
    {
      mpo.setOp(3*i-3,new Operator(this->id_fermi_mpo),true);
      fermi_id_saved=true;
      pos_fermi_id_saved=3*i-3;
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }
    if(link_id_saved)
      mpo.setOp(3*(i-1)-1,&mpo.getOp(pos_link_id_saved),false);
    else
    {
      mpo.setOp(3*(i-1)-1,new Operator(this->id_link_mpo),false);
      link_id_saved=true;
      pos_link_id_saved=3*(i-1)-1;
    }    
  }
}

void SU2Hamiltonian::getChargeCorrelatorMPO2(MPO &mpo, int length, int pos, double factor) 
{
  //Some error checking
  if(pos>this->N_total || pos<1 || (pos+length)>this->N_fermions)
  {
    cout << "Error in SU2Hamiltonian::getChargeCorrelatorMPO2(), desired correlator of length " << length<< " starting on position " << pos << " is outside the the chain of " << this->N_fermions << " spins" << endl;
    exit(666);
  }
  //Intialize the MPO
  mpo.initLength(this->N_total);
  
  if(!this->correlation_matrices_saved)
  {
    //I have two different two site operators which act locally. I split these via SVD and then simply put them where I need them
    mwArray qsquare,unitary;
    //Matrices for the SVD
    mwArray U,S,V;
    
    //Get two copies of sigma_z for the kronecker product (ist working in place, so for saftey I copy it)
    mwArray copy_1_s3,copy_2_s3;
    copy_1_s3=mwArray(this->sigma_z);
    copy_2_s3=mwArray(this->sigma_z);
    
    cout << "Using factor " << factor << endl;
    cout.flush();
    //Compute the two site operators
    qsquare = identityMatrix(4)-kron(copy_1_s3,copy_2_s3);
    //cout << "Sigma_z after kronecker product " << this->sigma_z << endl;
    wrapper::expm(qsquare,unitary,factor*M_PIl*I_c);
    //Reshape them accordingly
    qsquare.reshape(Indices(this->d_fermi,this->d_fermi,this->d_fermi,this->d_fermi));
    unitary.reshape(Indices(this->d_fermi,this->d_fermi,this->d_fermi,this->d_fermi));  
    //Permute the indices, such that the indices of a single body operator are next to each other
    qsquare.permute(Indices(1,3,2,4));
    unitary.permute(Indices(1,3,2,4));
    //Now take the first two indices together and perform a SVD
    qsquare.reshape(Indices(this->d_fermi*this->d_fermi,this->d_fermi*this->d_fermi));
    unitary.reshape(Indices(this->d_fermi*this->d_fermi,this->d_fermi*this->d_fermi));
    
    wrapper::svd(qsquare,U,S,V);
    S=sqrt(S);
    U.multiplyRight(S);	V.multiplyLeft(S);
    U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
    V.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
    V.permute(Indices(2,1,3,4));
    
    this->charge_left=mwArray(U);
    this->charge_right=mwArray(V);
    
    U.clear();	S.clear();	V.clear();
    
    wrapper::svd(unitary,U,S,V);
    S=sqrt(S);
    U.multiplyRight(S);	V.multiplyLeft(S);
    U.reshape(Indices(this->d_fermi,1,this->d_fermi,-1));
    V.reshape(Indices(-1,this->d_fermi,this->d_fermi,1));
    V.permute(Indices(2,1,3,4));
    
    this->unitary_left=mwArray(U);
    this->unitary_right=mwArray(V);
    
    /*cout << "charge_left: " << charge_left << endl;
    cout << "charge_right: " << charge_right << endl;
    cout << "unitary_left: " << unitary_left << endl;
    cout << "unitary_right: " << unitary_right << endl;*/
    //Set the flag that the matrices are now saved
    this->correlation_matrices_saved=true;
  }
  
  //Now fill the MPO
  mpo.setOp(3*pos-3,new Operator(charge_left),true);
  mpo.setOp(3*pos-2,new Operator(charge_right),true);
  mpo.setOp(3*pos-1,new Operator(this->id_link_mpo),true);
  mpo.setOp(3*(pos+length)-3,&mpo.getOp(3*pos-3),false);
  mpo.setOp(3*(pos+length)-2,&mpo.getOp(3*pos-2),false);
  
  //Sites in between the two charge correlations sites
  bool fermi_inside_saved=false;
  int pos_fermi_inside_saved;
  for(int i=(pos+1); i<(pos+length); i++)
  {
    if(fermi_inside_saved)
    {
      mpo.setOp(3*i-3,&mpo.getOp(3*pos_fermi_inside_saved-3),false);
      mpo.setOp(3*i-2,&mpo.getOp(3*pos_fermi_inside_saved-2),false);
    }
    else
    {
      mpo.setOp(3*i-3,new Operator(unitary_left),true);
      mpo.setOp(3*i-2,new Operator(unitary_right),true);
      fermi_inside_saved=true;
      pos_fermi_inside_saved=i;
    }
    mpo.setOp(3*i-1,&mpo.getOp(3*pos-1),false);      
  }
  
  //Identities to the left
  bool fermi_id_saved=false;
  int pos_fermi_id_saved;
  for(int i=1; i<pos; i++)
  {
    if(fermi_id_saved)
    {
      mpo.setOp(3*i-3,&mpo.getOp(pos_fermi_id_saved),false);
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }
    else
    {
      mpo.setOp(3*i-3,new Operator(this->id_fermi_mpo),true);
      fermi_id_saved=true;
      pos_fermi_id_saved=3*i-3;
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }
  }
  //Identities to the right
  for(int i=(pos+length+1); i<=this->N_fermions; i++)
  {
    if(fermi_id_saved)
    {
      mpo.setOp(3*i-3,&mpo.getOp(pos_fermi_id_saved),false);
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }
    else
    {
      mpo.setOp(3*i-3,new Operator(this->id_fermi_mpo),true);
      fermi_id_saved=true;
      pos_fermi_id_saved=3*i-3;
      mpo.setOp(3*i-2,&mpo.getOp(pos_fermi_id_saved),false);
    }  
  }
  //Remaining links
  for(int i=2; i<this->N_total; i+=3)
  {
    if(mpo.isEmpty(i))
      mpo.setOp(i,&mpo.getOp(3*pos-1),false);
  }   
  /*char buffer [50];
  sprintf (buffer, "Correlation_n%dl%d.m",pos,length);  
  mpo.exportForMatlab(buffer);*/
}

void SU2Hamiltonian::getCondensateMPO(MPO &mpo, bool factor_off) const
{
  cout << "Constructing condensate MPO" << endl;
  
  mpo.initLength(this->N_total);  
  
  double factor = sqrt(this->epsilon)/this->N_fermions;
  if(factor_off)
    factor=1.0;
  
  int Dbond = 2; 
  
  //Intialize the operators 
  mwArray Fermionic_operators,Gauge_operators;  
  Fermionic_operators=mwArray(Indices(2,this->d_fermi,this->d_fermi));	Fermionic_operators.fillWithZero();
  Gauge_operators=mwArray(Indices(1,this->d_link,this->d_link));	Gauge_operators.fillWithZero();  
  
  for(int i=0; i<this->d_fermi; i++)
  {
    for(int j=0; j<this->d_fermi; j++)
    {
      //Identity
      Fermionic_operators.setElement(this->id_fermi.getElement(Indices(i,j)),Indices(0,i,j));
      //sigma z
      Fermionic_operators.setElement(this->sigma_z.getElement(Indices(i,j)),Indices(1,i,j));
    }
  }  
  for(int i=0; i<this->d_link; i++)
    for(int j=0; j<this->d_link; j++)
      Gauge_operators.setElement(this->id_link.getElement(Indices(i,j)),Indices(0,i,j));
    
  Fermionic_operators.reshape(Indices(2,this->d_fermi*this->d_fermi));
  Gauge_operators.reshape(Indices(1,this->d_link*this->d_link));  
  
  mwArray FermiOdd, FermiEven, FermiFirst, FermiLast, Gauge;
  
  FermiFirst=mwArray(Indices(1,Dbond,2));	FermiFirst.fillWithZero();
  FermiOdd=mwArray(Indices(Dbond,Dbond,2));	FermiOdd.fillWithZero();
  FermiEven=mwArray(Indices(Dbond,Dbond,2));	FermiEven.fillWithZero();
  FermiLast=mwArray(Indices(Dbond,1,2));	FermiLast.fillWithZero();
  Gauge=mwArray(Indices(Dbond,Dbond,1));	Gauge.fillWithZero();
  
  FermiFirst.setElement(ONE_c,Indices(0,0,0));
  FermiFirst.setElement(-0.5*factor*ONE_c,Indices(0,1,1));
  
  FermiOdd.setElement(ONE_c,Indices(0,0,0));
  FermiOdd.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  FermiOdd.setElement(-0.5*factor*ONE_c,Indices(0,Dbond-1,0));
  FermiOdd.setElement(-0.5*factor*ONE_c,Indices(0,Dbond-1,1));
  
  FermiEven.setElement(ONE_c,Indices(0,0,0));
  FermiEven.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  FermiEven.setElement(0.5*factor*ONE_c,Indices(0,Dbond-1,0));
  FermiEven.setElement(0.5*factor*ONE_c,Indices(0,Dbond-1,1));
  
  FermiLast.setElement(ONE_c,Indices(1,0,0));
  FermiLast.setElement(pow(-1.0,this->N_fermions)*0.5*factor*ONE_c,Indices(0,0,1));
  
  Gauge.setElement(ONE_c,Indices(0,0,0));
  Gauge.setElement(ONE_c,Indices(Dbond-1,Dbond-1,0));
  
  //Reshape and contract
  FermiFirst.reshape(Indices(1*Dbond,2));	FermiFirst.multiplyRight(Fermionic_operators);	FermiFirst.reshape(Indices(1,Dbond,this->d_fermi,this->d_fermi));	FermiFirst.permute(Indices(3,1,4,2));
  FermiOdd.reshape(Indices(Dbond*Dbond,2));	FermiOdd.multiplyRight(Fermionic_operators);	FermiOdd.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	FermiOdd.permute(Indices(3,1,4,2));
  FermiEven.reshape(Indices(Dbond*Dbond,2));	FermiEven.multiplyRight(Fermionic_operators);	FermiEven.reshape(Indices(Dbond,Dbond,this->d_fermi,this->d_fermi));	FermiEven.permute(Indices(3,1,4,2));
  FermiLast.reshape(Indices(Dbond*1,2));	FermiLast.multiplyRight(Fermionic_operators);	FermiLast.reshape(Indices(Dbond,1,this->d_fermi,this->d_fermi));	FermiLast.permute(Indices(3,1,4,2));
  Gauge.reshape(Indices(Dbond*Dbond,1));	Gauge.multiplyRight(Gauge_operators);		Gauge.reshape(Indices(Dbond,Dbond,this->d_link,this->d_link));		Gauge.permute(Indices(3,1,4,2));
  
  mpo.setOp(0,new Operator(FermiFirst),true);
  mpo.setOp(1,new Operator(FermiOdd),true);
  mpo.setOp(2,new Operator(Gauge),true);
  mpo.setOp(this->N_total-1,new Operator(FermiLast),true);
  
  int pos_even_saved;
  bool even_saved=false;
  
  for(int i=2; i<this->N_fermions; i++)
  {
    if((i%2)==0)
    {
      //Even site
      if(even_saved)
      {
	mpo.setOp(3*i-3,&mpo.getOp(pos_even_saved),false);
	mpo.setOp(3*i-2,&mpo.getOp(pos_even_saved),false);
      }
      else
      {
	mpo.setOp(3*i-3,new Operator(FermiEven),true);
	pos_even_saved=3*i-3;
	mpo.setOp(3*i-2,&mpo.getOp(pos_even_saved),false);
      }
    }
    else
    {
      //Odd site, already saved outise the loop
      mpo.setOp(3*i-3,&mpo.getOp(1),false);
      mpo.setOp(3*i-2,&mpo.getOp(1),false);
    }
    mpo.setOp(3*i-1,&mpo.getOp(2),false);
  }
  
  //Last site
  if((this->N_fermions%2)==0)
  {
    //Even site
    if(even_saved)
      mpo.setOp(3*this->N_fermions-3,&mpo.getOp(pos_even_saved),false);
    else
      mpo.setOp(3*this->N_fermions-3,new Operator(FermiEven),true);
  }
  else
    mpo.setOp(3*this->N_fermions-3,&mpo.getOp(1),false); 
#ifdef MATOUT
  mpo.exportForMatlab("SU2Condensate.m");
#endif
}


void SU2Hamiltonian::splitLocalOp(mwArray Op,vector<int> dims,vector<mwArray> &matrices) const
{
  //Delete potential content of the vector
  matrices.clear();
  
  //Matrices for the SVDs
  mwArray U,S,V;  
  
  ////////////////////////////////////////////////
  //Some error checking
  ////////////////////////////////////////////////
  
  //Is the given Operator a matrix?
  if(!Op.isMatrix())
  {
    cout << "Error in splitLocalOp(), the given operator is not a matrix" << endl;
    exit(666);
  }
  
  //Check if the product of the given physical dimensions is compliant with the given mwArray
  int proddim=accumulate (dims.begin(),dims.end(),1,multiplies<int>());
  Indices matrix_dim=Op.getDimensions();
  if(matrix_dim!=Indices(proddim,proddim))
  {
    cout << "Error in splitLocalOp(), the given physical dimensions are not compliant with the given operator" << endl;
    cout << "Dimension of matrix " << matrix_dim << ", product of dimensions " << proddim << endl;
    exit(666);
  }
  
  ////////////////////////////////////////////////
  //Check for special cases
  ////////////////////////////////////////////////
  if(dims.size()==1)
  {
    //Single site operator given, all I have to do is to reshape it accordingly to get an object with bond dimension one
    cout << "Single site operator given, will only reshape it" << endl;
    matrices.push_back(reshape(Op,Indices(dims[0],1,dims[0],1)));
    return;
  }
  if(dims.size()==2)
  {
    //I only need a single SVD
    cout << "Two site operator given, only computing a single SVD" << endl;
    Op.reshape(Indices(dims[0],dims[1],dims[0],dims[1]));
    Op.permute(Indices(1,3,2,4));
    Op.reshape(Indices(dims[0]*dims[0],dims[1]*dims[1]));
    wrapper::svd(Op,U,S,V);  
    //Redistribute the singular values
    S=sqrt(S);
    U.multiplyRight(S);
    V.multiplyLeft(S);
    //Reshape the matrix and save results
    U.reshape(Indices(dims[0],1,dims[0],-1));
    matrices.push_back(U);
    V.reshape(Indices(-1,dims[1],dims[1],1));
    V.permute(Indices(2,1,3,4));
    matrices.push_back(V);  
    return;
  }  
  ////////////////////////////////////////////////
  //First permute the indices in the right order  
  ////////////////////////////////////////////////  
  //Dimensions of the input matrix
  vector <int> indices;
  //As an operator has two physical legs and the kronecker product in the code orders the indices as described above, I simply double the physical dimensions
  indices = dims;	indices.insert(indices.end(), dims.begin(), dims.end());  
  //Reshape
  Op.reshape(Indices(indices));  
  //Now compute the order of the indices
  indices.clear();
  vector<int> opdims;
  for(int i=1; i<=dims.size(); i++)
  {
    indices.push_back(i);
    indices.push_back(i+dims.size());
    opdims.push_back(dims[i-1]);
    opdims.push_back(dims[i-1]);
  }  
  //Permute the indices accordingly, such that two physical legs of a single site operator are successive
  Op.permute(Indices(indices));
  
  ////////////////////////////////////////////////
  //Compute the SVDs  
  ////////////////////////////////////////////////   
  //Now take the first two indices together and do a SVD  
  Op.reshape(Indices(accumulate(opdims.begin(),opdims.begin()+2,1,multiplies<int>()),accumulate(opdims.begin()+2,opdims.end(),1,multiplies<int>())));
  wrapper::svd(Op,U,S,V);  
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape the matrix and save results
  U.reshape(Indices(opdims[0],1,opdims[1],-1));
  matrices.push_back(U);
  Op=V;  
  //Proceed with the SVDs up to the last one
  int count = 2;
  Indices dim;
  for (int i=2; i<(dims.size()-1); i++)
  {
    dim=Op.getDimensions();
    Op.reshape(Indices(dim[0]*accumulate(opdims.begin()+count,opdims.begin()+count+2,1,multiplies<int>()),accumulate(opdims.begin()+count+2,opdims.end(),1,multiplies<int>())));
    wrapper::svd(Op,U,S,V);  
    //Redistribute the singular values
    S=sqrt(S);
    U.multiplyRight(S);
    V.multiplyLeft(S);
    //Reshape the matrix and save results
    U.reshape(Indices(dim[0],opdims[count],opdims[count+1],-1));
    U.permute(Indices(2,1,3,4));
    matrices.push_back(U);
    Op=V;
    count+=2;
  }  
  //Last SVD
  dim=Op.getDimensions();
  Op.reshape(Indices(dim[0]*accumulate(opdims.begin()+count,opdims.begin()+count+2,1,multiplies<int>()),accumulate(opdims.begin()+count+2,opdims.end(),1,multiplies<int>())));
  wrapper::svd(Op,U,S,V);  
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape the matrix and save results
  U.reshape(Indices(dim[0],opdims[count],opdims[count+1],-1));
  U.permute(Indices(2,1,3,4));
  matrices.push_back(U);
  V.reshape(Indices(-1,opdims[count+2],opdims[count+3],1));
  V.permute(Indices(2,1,3,4));
  matrices.push_back(V);  
  cout << "Successfully computed SVDs" << endl;
}