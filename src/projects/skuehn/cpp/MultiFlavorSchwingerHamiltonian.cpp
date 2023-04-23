#include "MultiFlavorSchwingerHamiltonian.h"
#include "Indices.h"
#include "Contractor.h"
#include "splitLocalOp.h"
#include <numeric>
#include <algorithm>

using namespace std;
using namespace shrt;


MultiFlavorSchwingerHamiltonian::MultiFlavorSchwingerHamiltonian(int N_sites_,int N_flavor_,double x_,vector<double> mu_,vector<double> nu_,double lambda_,double l0_):hamiltonian(N_flavor_*N_sites_)
{
  //Physical dimension of a fermionic site
  this->d_fermi=2;
  //Number of sites
  this->N_sites = N_sites_;
  //Number of flavors
  this->N_flavor=N_flavor_;
  //Total number of sites
  this->N_total=this->N_flavor*this->N_sites;
  //Mass
  if(mu_.size()==this->N_flavor)
    this->mu=mu_;
  else
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian(): Number of masses given is not compatible with the number of flavors" << endl;
    exit(1);
  }
  //Chemical potential
  if(nu_.size()==this->N_flavor)
    this->nu=nu_;
  else
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian(): Number of chemical potentials given is not compatible with the number of flavors" << endl;
    exit(1);
  }
  //Hopping constant
  this->x=x_;
  //Background field
  this->l0=l0_;
  //Charge penalty strength
  this->lambda=lambda_;
  //Penalty strength (unused in this case)
  this->eta=0.0;
  //Flavor which should be fixed to a particle number (unused in this case thus set to negative value)
  this->penalized_flavor=-1;
  this->Npenalized_flavor=-1;
  
  initOperators();
  initHamiltonian();
}

MultiFlavorSchwingerHamiltonian::MultiFlavorSchwingerHamiltonian(int N_sites_,int N_flavor_,double x_,vector<double> mu_,vector<double> nu_,double lambda_,double l0_,int penalized_flavor_, int Npenalized_flavor_, double eta_):hamiltonian(N_flavor_*N_sites_)
{
  //Physical dimension of a fermionic site
  this->d_fermi=2;
  //Number of sites
  this->N_sites = N_sites_;
  //Number of flavors
  this->N_flavor=N_flavor_;
  //Total number of sites
  this->N_total=this->N_flavor*this->N_sites;
  //Mass
  if(mu_.size()==this->N_flavor)
    this->mu=mu_;
  else
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian(): Number of masses given is not compatible with the number of flavors" << endl;
    exit(1);
  }
  //Chemical potential
  if(nu_.size()==this->N_flavor)
    this->nu=nu_;
  else
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian(): Number of chemical potentials given is not compatible with the number of flavors" << endl;
    exit(1);
  }
  //Hopping constant
  this->x=x_;
  //Background field
  this->l0=l0_;
  //Charge penalty strength
  this->lambda=lambda_;
  //Penalty strength
  this->eta=eta_;
  if(this->eta==0)
    cout << "Warning in MultiFlavorSchwingerHamiltonian::MultiFlavorSchwingerHamiltonian(), penalty value is zero, more efficient penalty free Hamiltonian version could be used" << endl;
  //Flavor which should be fixed to a particle number (unused in this case thus set to negative value)
  this->penalized_flavor=penalized_flavor_;
  this->Npenalized_flavor=Npenalized_flavor_;
  if(this->penalized_flavor<0 || this->penalized_flavor>=this->N_flavor)
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::MultiFlavorSchwingerHamiltonian(), flavor which should be fixed to a particle number is not valid, program will be aborted" << endl;
    exit(666);
  }
  if(this->Npenalized_flavor<0 || this->penalized_flavor>this->N_sites)
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::MultiFlavorSchwingerHamiltonian(), flavor " << this->penalized_flavor << " can only be set to values in between 0 and " << this->N_sites << ", program will be aborted" << endl;
    exit(666);
  }
  
  initOperators();
  initHamiltonianPenalty();
}

MultiFlavorSchwingerHamiltonian::~MultiFlavorSchwingerHamiltonian()
{
  this->Z.clear();
}

void MultiFlavorSchwingerHamiltonian::initOperators()
{
  Z=mwArray(Indices(6,this->d_fermi,this->d_fermi));
  
  //Identity Matrix
  mwArray id=identityMatrix(this->d_fermi);  
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
  
  for(int i=0; i<2; i++)
    for(int k=0; k<2; k++)
    {
      Z.setElement(id.getElement(Indices(i,k)),Indices(0,i,k));
      Z.setElement(sig_p.getElement(Indices(i,k)),Indices(1,i,k));
      Z.setElement(sig_m.getElement(Indices(i,k)),Indices(2,i,k));
      Z.setElement(sig_z.getElement(Indices(i,k)),Indices(3,i,k));
      Z.setElement(sig_y.getElement(Indices(i,k)),Indices(4,i,k));
      Z.setElement(sig_x.getElement(Indices(i,k)),Indices(5,i,k));      
    }   

  //I need often identities to construct a single body operator, therefore I prepare those and save them
  this->id_fermi_mpo=id;
  this->id_fermi_mpo.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
}
void MultiFlavorSchwingerHamiltonian::initHamiltonian()
{
  //Get the bond dimension
  int D=3+2*this->N_flavor;
  mwArray Site;
  
  //Site dependent constants alpha
  vector <double> alpha;
  
  alpha.resize(this->N_sites,this->l0);
  for(int i=0; i<alpha.size(); i+=2)
    alpha[i]= this->l0+this->N_flavor/2.0;
  
  double constant;
  int site_index,index1,index2,index3;   
  index3 = 2*this->N_flavor+1; 
  
  //cout << "Alphas: " << alpha << endl;
  cout << "Constructing Hamiltonian MPO with bond dimension " << D << endl;
  //cout << "Index3=" << index3 << endl;
  
  //Right now the penalty is only working for an even number of sites, therefore I give a warning if the total number of sites is not even
  if(this->lambda!=0 && (this->N_sites%2)!=0)
    cout << "Warning in MultiFlavorSchwingerHamiltonian::initHamiltonian(), penalty only supports an even number of sites" << endl;
  
  for(int n=0; n<this->N_sites; n++)
    for(int f=0; f<this->N_flavor; f++)
    {
      //Compute the site index (the "p")
      site_index = n*this->N_flavor+f;
      index1=2*f+1;
      index2=index1+1;;
      //cout << "Site index: " << site_index << ", site=" << n << ", flavor=" << f <<", index1=" << index1 << ", index2=" << index2 << ", total number of sites " << this->N_total <<  endl;
      
      if(site_index==0)
      {
	//First site, has to be treated separately
	Site = mwArray(Indices(1,D,6));	Site.fillWithZero();
	
	//Identity
	constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f])+1.0/this->N_flavor*(alpha[n]*alpha[n]+(n+1.0)*this->N_flavor/4.0)+this->lambda;    
	Site.setElement(ONE_c,Indices(0,0,0));
	Site.setElement(constant*ONE_c,Indices(0,D-1,0));
	
	//sigma_p
	Site.setElement(-this->x*ONE_c,Indices(0,index1,1));
	
	//sigma_m
	Site.setElement(-this->x*ONE_c,Indices(0,index2,2));
	
	//sigma_z
	Site.setElement(ONE_c,Indices(0,index3,3));
	double tilde_alpha=std::accumulate(alpha.begin()+n, alpha.end()-1, 0.0);
	constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f])+tilde_alpha;
	Site.setElement(constant*ONE_c,Indices(0,D-1,3));
	
	Site.reshape(Indices(1*D,6));	Site.multiplyRight(reshape(Z,Indices(6,this->d_fermi*this->d_fermi)));
	Site.reshape(Indices(1,D,this->d_fermi,this->d_fermi));	Site.permute(Indices(3,1,4,2));
      }
      else if(site_index==(this->N_total-1))
      {
	//Last site, has to be treated separately
	Site = mwArray(Indices(D,1,6));	Site.fillWithZero();
	
	//Identity
	constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f])+this->lambda;
	Site.setElement(ONE_c,Indices(D-1,0,0));
	Site.setElement(constant*ONE_c,Indices(0,0,0));
	
	//sigma_p
	Site.setElement(ONE_c,Indices(index2,0,1));
	
	//sigma_m
	Site.setElement(ONE_c,Indices(index1,0,2));
	
	//sigma_z
	//Single body term
	constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f]);
	Site.setElement(constant*ONE_c,Indices(0,0,3));
	//Long range spin penalty
	Site.setElement(2.0*this->lambda*ONE_c,Indices(index3,0,3));
	
	Site.reshape(Indices(D*1,6));	Site.multiplyRight(reshape(Z,Indices(6,this->d_fermi*this->d_fermi)));
	Site.reshape(Indices(D,1,this->d_fermi,this->d_fermi));	Site.permute(Indices(3,1,4,2));
      }
      else
      {
	//Sites in between
	Site = mwArray(Indices(D,D,6));	Site.fillWithZero();
	
	//Identity
	if(n==(this->N_sites-1))
	  constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f])+this->lambda;
	else
	  constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f])+1.0/this->N_flavor*(alpha[n]*alpha[n]+(n+1.0)*this->N_flavor/4.0)+this->lambda;
	Site.setElement(ONE_c,Indices(0,0,0));
	Site.setElement(ONE_c,Indices(D-1,D-1,0));
	Site.setElement(constant*ONE_c,Indices(0,D-1,0));
	Site.setElement(ONE_c,Indices(index3,index3,0));
	
	//sigma_p
	Site.setElement(-this->x*ONE_c,Indices(0,index1,1));
	Site.setElement(ONE_c,Indices(index2,D-1,1));
	
	//sigma_m
	Site.setElement(-this->x*ONE_c,Indices(0,index2,2));
	Site.setElement(ONE_c,Indices(index1,D-1,2));
	
	//sigma_z
	for(int k=0; k<this->N_flavor; k++)
	{
	  //Pass on sigma_z for the other flavors
	  if(((2*k+1)!=index1) && (2*k+2)!=index2)
	  {
	    //cout << "Passing on for flavor " << k << endl;
	    Site.setElement(I_c,Indices(2*k+1,2*k+1,3));
	    Site.setElement(-I_c,Indices(2*k+2,2*k+2,3));
	  }
	}
	//Long range term (has to end before the last site!)
	if(n==(this->N_sites-1))
	  constant=2.0*this->lambda;
	else
	  constant=0.5*(this->N_sites-1-n)+2.0*this->lambda;
	Site.setElement(ONE_c,Indices(0,index3,3));
	Site.setElement(constant*ONE_c,Indices(index3,D-1,3));
	//Single body term
	if(n==(this->N_sites-1))
	  constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f]);
	else
	{
	  double tilde_alpha=std::accumulate(alpha.begin()+n, alpha.end()-1, 0.0);
	  //cout << "alpha=" << alpha << endl;
	  //cout << "tilde alpha=" << tilde_alpha << endl;
	  constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f])+tilde_alpha;
	}
	Site.setElement(constant*ONE_c,Indices(0,D-1,3));
	
	Site.reshape(Indices(D*D,6));	Site.multiplyRight(reshape(Z,Indices(6,this->d_fermi*this->d_fermi)));
	Site.reshape(Indices(D,D,this->d_fermi,this->d_fermi));	Site.permute(Indices(3,1,4,2));
      }     
      
      hamiltonian.setOp(site_index,new Operator(Site),true);
    }
  //cout << "Hamiltonian: " << this->hamiltonian << endl;
#ifdef MATOUT
  stringstream name;    
  name.str("");
  name << "MFSH_N" << this->N_sites << "F" << this->N_flavor << "x" << (int) (x*1000);
  for (int i=0; i<this->N_flavor; i++)
    name << "mu" << (int) (this->mu[i]*1000);
  for (int i=0; i<this->N_flavor; i++)
    name << "nu" << (int) (this->nu[i]*1000);
  name << "la" << (int) (this->lambda*1000) ;
  name << "lbg" << (int) (this->l0*1000) << ".m"; 
  
  hamiltonian.exportForMatlab(name.str().c_str(),11);
#endif
}

void MultiFlavorSchwingerHamiltonian::initHamiltonianPenalty()
{
  //Get the bond dimension
  int D=3+2*this->N_flavor+1;
  mwArray Site;
  
  //Site dependent constants alpha
  vector <double> alpha;
  
  alpha.resize(this->N_sites,this->l0);
  for(int i=0; i<alpha.size(); i+=2)
    alpha[i]= this->l0+this->N_flavor/2.0;
  
  double constant;
  int site_index,index1,index2,index3;   
  index3 = 2*this->N_flavor+1; 
  
  //cout << "Alphas: " << alpha << endl;
  cout << "Constructing Hamiltonian MPO with additional penalty imposing the particle number " << this->Npenalized_flavor << " for flavor " << this-> penalized_flavor << " with bond dimension " << D << endl;
  //cout << "Index3=" << index3 << endl;
  
  //Right now the penalty is only working for an even number of sites, therefore I give a warning if the total number of sites is not even
  if(this->lambda!=0 && (this->N_sites%2)!=0)
    cout << "Warning in MultiFlavorSchwingerHamiltonian::initHamiltonianPenalty(), penalty only supports an even number of sites" << endl;
  
  cout << "Putting penalty on flavor " << this->penalized_flavor << " with strength " << this->eta << endl;
  
  for(int n=0; n<this->N_sites; n++)
    for(int f=0; f<this->N_flavor; f++)
    {
      //Compute the site index (the "p")
      site_index = n*this->N_flavor+f;
      index1=2*f+1;
      index2=index1+1;;
      //cout << "Site index: " << site_index << ", site=" << n << ", flavor=" << f <<", index1=" << index1 << ", index2=" << index2 << ", total number of sites " << this->N_total <<  endl;
      
      if(site_index==0)
      {
	//First site, has to be treated separately
	Site = mwArray(Indices(1,D,6));	Site.fillWithZero();
	
	//Identity
	constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f])+1.0/this->N_flavor*(alpha[n]*alpha[n]+(n+1.0)*this->N_flavor/4.0)+this->lambda;
	if(f==penalized_flavor)
	  constant += this->eta*(0.25*this->N_sites + pow(this->Npenalized_flavor,2.0)/this->N_sites-this->Npenalized_flavor + 0.25);
	Site.setElement(ONE_c,Indices(0,0,0));
	Site.setElement(constant*ONE_c,Indices(0,D-1,0));
	
	//sigma_p
	Site.setElement(-this->x*ONE_c,Indices(0,index1,1));
	
	//sigma_m
	Site.setElement(-this->x*ONE_c,Indices(0,index2,2));
	
	//sigma_z
	Site.setElement(ONE_c,Indices(0,index3,3));
	double tilde_alpha=std::accumulate(alpha.begin()+n, alpha.end()-1, 0.0);
	constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f])+tilde_alpha;
	if(f==penalized_flavor)
	{
	  constant += this->eta*(0.5*this->N_sites-this->Npenalized_flavor);
	  //Start long-range term from particle number penalty on flavor 
	  Site.setElement(ONE_c,Indices(0,D-2,3));
	}
	Site.setElement(constant*ONE_c,Indices(0,D-1,3));	
	
	Site.reshape(Indices(1*D,6));	Site.multiplyRight(reshape(Z,Indices(6,this->d_fermi*this->d_fermi)));
	Site.reshape(Indices(1,D,this->d_fermi,this->d_fermi));	Site.permute(Indices(3,1,4,2));
      }
      else if(site_index==(this->N_total-1))
      {
	//Last site, has to be treated separately
	Site = mwArray(Indices(D,1,6));	Site.fillWithZero();
	
	//Identity
	constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f])+this->lambda;
	if(f==penalized_flavor)
	  constant += this->eta*(0.25*this->N_sites + pow(this->Npenalized_flavor,2.0)/this->N_sites-this->Npenalized_flavor + 0.25);
	Site.setElement(ONE_c,Indices(D-1,0,0));
	Site.setElement(constant*ONE_c,Indices(0,0,0));
	
	//sigma_p
	Site.setElement(ONE_c,Indices(index2,0,1));
	
	//sigma_m
	Site.setElement(ONE_c,Indices(index1,0,2));
	
	//sigma_z
	//Single body term
	constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f]);
	if(f==penalized_flavor)
	{
	  constant += this->eta*(0.5*this->N_sites-this->Npenalized_flavor);
	  //Finish long-range term from particle number penalty on flavor 
	  Site.setElement(0.5*this->eta*ONE_c,Indices(D-2,0,3));
	}
	Site.setElement(constant*ONE_c,Indices(0,0,3));
	//Long range spin penalty
	Site.setElement(2.0*this->lambda*ONE_c,Indices(index3,0,3));
	
	Site.reshape(Indices(D*1,6));	Site.multiplyRight(reshape(Z,Indices(6,this->d_fermi*this->d_fermi)));
	Site.reshape(Indices(D,1,this->d_fermi,this->d_fermi));	Site.permute(Indices(3,1,4,2));
      }
      else
      {
	//Sites in between
	Site = mwArray(Indices(D,D,6));	Site.fillWithZero();
	
	//Identity
	if(n==(this->N_sites-1))
	  constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f])+this->lambda;
	else
	  constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f])+1.0/this->N_flavor*(alpha[n]*alpha[n]+(n+1.0)*this->N_flavor/4.0)+this->lambda;
	if(f==penalized_flavor)
	  constant += this->eta*(0.25*this->N_sites + pow(this->Npenalized_flavor,2.0)/this->N_sites-this->Npenalized_flavor + 0.25);
	Site.setElement(ONE_c,Indices(0,0,0));
	Site.setElement(ONE_c,Indices(D-1,D-1,0));
	Site.setElement(constant*ONE_c,Indices(0,D-1,0));
	Site.setElement(ONE_c,Indices(index3,index3,0));
	//Pass on long range-interaction for particle penalty on flavor
	Site.setElement(ONE_c,Indices(D-2,D-2,0));
	
	
	//sigma_p
	Site.setElement(-this->x*ONE_c,Indices(0,index1,1));
	Site.setElement(ONE_c,Indices(index2,D-1,1));
	
	//sigma_m
	Site.setElement(-this->x*ONE_c,Indices(0,index2,2));
	Site.setElement(ONE_c,Indices(index1,D-1,2));
	
	//sigma_z
	for(int k=0; k<this->N_flavor; k++)
	{
	  //Pass on sigma_z for the other flavors
	  if(((2*k+1)!=index1) && (2*k+2)!=index2)
	  {
	    //cout << "Passing on for flavor " << k << endl;
	    Site.setElement(I_c,Indices(2*k+1,2*k+1,3));
	    Site.setElement(-I_c,Indices(2*k+2,2*k+2,3));
	  }
	}
	//Long range term (has to end before the last site!)
	if(n==(this->N_sites-1))
	  constant=2.0*this->lambda;
	else
	  constant=0.5*(this->N_sites-1-n)+2.0*this->lambda;
	Site.setElement(ONE_c,Indices(0,index3,3));
	Site.setElement(constant*ONE_c,Indices(index3,D-1,3));
	//Single body term
	if(n==(this->N_sites-1))
	  constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f]);
	else
	{
	  double tilde_alpha=std::accumulate(alpha.begin()+n, alpha.end()-1, 0.0);
	  //cout << "alpha=" << alpha << endl;
	  //cout << "tilde alpha=" << tilde_alpha << endl;
	  constant=0.5*(this->mu[f]*pow(-1.0,n)+this->nu[f])+tilde_alpha;
	}
	if(f==penalized_flavor)
	{
	  constant += this->eta*(0.5*this->N_sites-this->Npenalized_flavor);
	  //Long-range interaction for particle number penalty on flavor
	  Site.setElement(ONE_c,Indices(0,D-2,3));
	  Site.setElement(0.5*this->eta*ONE_c,Indices(D-2,D-1,3));
	}
	Site.setElement(constant*ONE_c,Indices(0,D-1,3));
	
	Site.reshape(Indices(D*D,6));	Site.multiplyRight(reshape(Z,Indices(6,this->d_fermi*this->d_fermi)));
	Site.reshape(Indices(D,D,this->d_fermi,this->d_fermi));	Site.permute(Indices(3,1,4,2));
      }     
      
      hamiltonian.setOp(site_index,new Operator(Site),true);
    }
  //cout << "Hamiltonian: " << this->hamiltonian << endl;
#ifdef MATOUT
  stringstream name;    
  name.str("");
  name << "MFSH_N" << this->N_sites << "F" << this->N_flavor << "x" << (int) (x*1000);
  for (int i=0; i<this->N_flavor; i++)
    name << "mu" << (int) (this->mu[i]*1000);
  for (int i=0; i<this->N_flavor; i++)
    name << "nu" << (int) (this->nu[i]*1000);
  name << "la" << (int) (this->lambda*1000) ;
  name << "lbg" << (int) (this->l0*1000) << ".m"; 
  
  hamiltonian.exportForMatlab(name.str().c_str(),11);
#endif
}


void MultiFlavorSchwingerHamiltonian::getSingleBodyMPO(MPO &mpo, int site, int flavor, SingleBodyOperator Op) const
{
  //First clear the MPO and resize it
  mpo.initLength(this->N_total);
  
  //Prepare the matrices
  mwArray Fermi_op;
  
  int index = this->N_flavor*site+flavor;
  
  switch(Op)
  {
    case sxOp:
      Fermi_op=this->Z.subArray(Indices(5,-1,-1)); Fermi_op.reshape(Indices(this->d_fermi,1,this->d_fermi,1)); break;
    case syOp:
      Fermi_op=this->Z.subArray(Indices(4,-1,-1)); Fermi_op.reshape(Indices(this->d_fermi,1,this->d_fermi,1)); break;
    case szOp:
      Fermi_op=this->Z.subArray(Indices(3,-1,-1)); Fermi_op.reshape(Indices(this->d_fermi,1,this->d_fermi,1)); break;
    case spOp:
      Fermi_op=this->Z.subArray(Indices(1,-1,-1)); Fermi_op.reshape(Indices(this->d_fermi,1,this->d_fermi,1)); break;
    case smOp:
      Fermi_op=this->Z.subArray(Indices(2,-1,-1)); Fermi_op.reshape(Indices(this->d_fermi,1,this->d_fermi,1)); break;
    case chargeOp:
      Fermi_op = 0.5*(pow(-1.0,site)*this->Z.subArray(Indices(0,-1,-1)) + this->Z.subArray(Indices(3,-1,-1)));      
      Fermi_op.reshape(Indices(Indices(this->d_fermi,1,this->d_fermi,1)));
      break;
    case condensateOp:
      Fermi_op = 0.5*sqrt(this->x)/this->N_sites*pow(-1.0,site)*(this->Z.subArray(Indices(0,-1,-1)) + this->Z.subArray(Indices(3,-1,-1)));      
      Fermi_op.reshape(Indices(Indices(this->d_fermi,1,this->d_fermi,1)));
      break;
    case condensate_nofactor_Op:
      Fermi_op = 0.5*pow(-1.0,site)*(this->Z.subArray(Indices(0,-1,-1)) + this->Z.subArray(Indices(3,-1,-1)));      
      Fermi_op.reshape(Indices(Indices(this->d_fermi,1,this->d_fermi,1)));
      break;
    default:
      cout << "Unknown operator" << endl;
      exit(666);
  } 
  
  //Now set operators
  bool saved=false;
  int saved_site;
  for(int k=0; k<this->N_total; k++)
  {
    if(saved && (k!=index))
      mpo.setOp(k,&mpo.getOp(saved_site),false);
    else
    {
      if(k==index)
	mpo.setOp(k,new Operator(Fermi_op),true);
      else
      {
	mpo.setOp(k,new Operator(this->id_fermi_mpo),true);
	saved=true;
	saved_site=k;
      }
    }
  }  
}

void MultiFlavorSchwingerHamiltonian::getSpinMPO(MPO &mpo, int site, int flavor, Component comp) const
{
  //Some error checking
  if(site<0 || site>=this->N_sites)
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::getSpinMPO(), given site index is invalid" << endl;
    exit(666);
  }
  if(flavor<0 || flavor>=this->N_flavor)
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::getSpinMPO(), given flavor is invalid" << endl;
    exit(666);
  }
  //Call the more general routine
  switch(comp)
  {
    case x_comp:
      getSingleBodyMPO(mpo,site,flavor,sxOp); break;
    case y_comp:
      getSingleBodyMPO(mpo,site,flavor,syOp); break;
    case z_comp:
      getSingleBodyMPO(mpo,site,flavor,szOp); break;
  }
}

void MultiFlavorSchwingerHamiltonian::getCondensateMPO(MPO &mpo, int site, int flavor, bool factor_off) const
{
  //Some error checking
  if(site<0 || site>=this->N_sites)
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::getCondensateMPO(), given site index is invalid" << endl;
    exit(666);
  }
  if(flavor<0 || flavor>=this->N_flavor)
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::getCondensateMPO(), given flavor is invalid" << endl;
    exit(666);
  }
  if(factor_off)
    this->getSingleBodyMPO(mpo,site,flavor,condensate_nofactor_Op);  
  else
    this->getSingleBodyMPO(mpo,site,flavor,condensateOp);  
}

void MultiFlavorSchwingerHamiltonian::getCondensateMPO(MPO &mpo,bool factor_off) const
{
  //Simply a wrapper
  this->MultiFlavorSchwingerHamiltonian::getCondensateMPO(mpo,-1,factor_off);
}

void MultiFlavorSchwingerHamiltonian::getCondensateMPO(MPO &mpo,int flavor, bool factor_off) const
{
  //Some error checking
  if(flavor>=this->N_flavor)
  {
    cout << "Warning in MultiFlavorSchwingerHamiltonian::getCondensateMPO(), specified flavor " << flavor << " larger than the number of flavors, will sum over all flavors" << endl;
    flavor=-1;
  }
  
  mpo.initLength(this->N_total);
  
  int D=2;
  double factor=factor_off?0.5:0.5*sqrt(this->x)/((double) this->N_sites);
  double sign=(((this->N_sites-1)%2)==0)?1.0:-1.0;
    
  mwArray First,Last,Odd,Even,Middle;
  First=mwArray(Indices(1,D,6));	First.fillWithZero();
  Odd=mwArray(Indices(D,D,6));		Odd.fillWithZero();
  Even=mwArray(Indices(D,D,6));		Even.fillWithZero();
  Middle=mwArray(Indices(D,D,6));	Middle.fillWithZero();
  Last=mwArray(Indices(D,1,6));		Last.fillWithZero();
  
  if(flavor==0 || flavor<0)
  {
    First.setElement(ONE_c,Indices(0,0,0));
    First.setElement(factor*ONE_c,Indices(0,1,0));
    First.setElement(factor*ONE_c,Indices(0,1,3));
  }
  else
    First.setElement(ONE_c,Indices(0,0,0));
  
  Odd.setElement(ONE_c,Indices(0,0,0));
  Odd.setElement(ONE_c,Indices(D-1,D-1,0));
  Odd.setElement(-factor*ONE_c,Indices(0,1,0));
  Odd.setElement(-factor*ONE_c,Indices(0,1,3)); 
  
  Even.setElement(ONE_c,Indices(0,0,0));
  Even.setElement(ONE_c,Indices(D-1,D-1,0));
  Even.setElement(factor*ONE_c,Indices(0,1,0));
  Even.setElement(factor*ONE_c,Indices(0,1,3));
  
  Middle.setElement(ONE_c,Indices(0,0,0));
  Middle.setElement(ONE_c,Indices(D-1,D-1,0));
  
  if(flavor==(this->N_flavor-1) || flavor<0)
  {
    Last.setElement(sign*factor*ONE_c,Indices(0,0,0));
    Last.setElement(ONE_c,Indices(1,0,0));
    Last.setElement(sign*factor*ONE_c,Indices(0,0,3));
  }
  else
    Last.setElement(ONE_c,Indices(1,0,0));
  
  //Reshape and contract
  First.reshape(Indices(1*D,6)); First.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); First.reshape(Indices(1,D,this->d_fermi,this->d_fermi)); First.permute(Indices(3,1,4,2));
  Odd.reshape(Indices(D*D,6)); Odd.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Odd.reshape(Indices(D,D,this->d_fermi,this->d_fermi)); Odd.permute(Indices(3,1,4,2));
  Even.reshape(Indices(D*D,6)); Even.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Even.reshape(Indices(D,D,this->d_fermi,this->d_fermi)); Even.permute(Indices(3,1,4,2));
  Middle.reshape(Indices(D*D,6)); Middle.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Middle.reshape(Indices(D,D,this->d_fermi,this->d_fermi)); Middle.permute(Indices(3,1,4,2));
  Last.reshape(Indices(D*1,6)); Last.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Last.reshape(Indices(D,1,this->d_fermi,this->d_fermi)); Last.permute(Indices(3,1,4,2));
  
  bool odd_saved=false, even_saved=false, middle_saved=false;
  int pos_odd_saved, pos_even_saved, pos_middle_saved, site_index;
  
  for(int n=0; n<this->N_sites; n++)
    for(int f=0; f<this->N_flavor; f++)
    {
      site_index=n*this->N_flavor+f;
      if(site_index==0)
	mpo.setOp(site_index,new Operator(First),true);
      else if(site_index==(this->N_total-1))
	mpo.setOp(site_index,new Operator(Last),true);
      else if(f==flavor || flavor<0)
      {
	//Site of the kind of flavor, where I would like to compute the condensate thereof
	if((n%2)==0)
	{
	  //Even site
	  if(even_saved)
	    mpo.setOp(site_index,&mpo.getOp(pos_even_saved),false);
	  else
	  {
	    mpo.setOp(site_index,new Operator(Even),true);
	    even_saved=true;
	    pos_even_saved=site_index;
	  }
	}
	else
	{
	  //Odd site
	  if(odd_saved)
	    mpo.setOp(site_index,&mpo.getOp(pos_odd_saved),false);
	  else
	  {
	    mpo.setOp(site_index,new Operator(Odd),true);
	    odd_saved=true;
	    pos_odd_saved=site_index;
	  }
	}
      }
      else
      {
	//Just some site in between (in principle to the left a matrix [1 0; 0 0] and to the right a matrix [0 0; 1 0] would be sufficient, but as there cannot be a one on the additional "channels" the middle matrix has, that is not a problem and I can spare two matrices
	if(middle_saved)
	  mpo.setOp(site_index,&mpo.getOp(pos_middle_saved),false);
	else
	{
	  mpo.setOp(site_index,new Operator(Middle),true);
	  middle_saved=true;
	  pos_middle_saved=site_index;
	}
      }	
    }  
#ifdef MATOUT
  stringstream name;    
  name.str("");
  if(flavor>=0)
    name << "CondensateMPO_f" << flavor << ".m";
  else
    name << "CondensateMPO_all_flavor.m";
  mpo.exportForMatlab(name.str().c_str());
#endif
}

void MultiFlavorSchwingerHamiltonian::getParticleNumberMPO(MPO &mpo, int flavor) const
{
  //Some error checking
  if(flavor>=this->N_flavor)
  {
    cout << "Warning in MultiFlavorSchwingerHamiltonian::getParticleNumberMPO(), specified flavor " << flavor << " larger than the number of flavors, will sum over all flavors" << endl;
    flavor=-1;
  }
  
  mpo.initLength(this->N_total);
  
  int D=2;
  double factor=0.5;
    
  mwArray First,Last,Op,Middle;
  First=mwArray(Indices(1,D,6));	First.fillWithZero();
  Op=mwArray(Indices(D,D,6));		Op.fillWithZero();
  Middle=mwArray(Indices(D,D,6));	Middle.fillWithZero();
  Last=mwArray(Indices(D,1,6));		Last.fillWithZero();
  
  if(flavor==0 || flavor<0)
  {
    First.setElement(ONE_c,Indices(0,0,0));
    First.setElement(factor*ONE_c,Indices(0,1,0));
    First.setElement(factor*ONE_c,Indices(0,1,3));
  }
  else
    First.setElement(ONE_c,Indices(0,0,0));
  
  Op.setElement(ONE_c,Indices(0,0,0));
  Op.setElement(ONE_c,Indices(D-1,D-1,0));
  Op.setElement(factor*ONE_c,Indices(0,1,0));
  Op.setElement(factor*ONE_c,Indices(0,1,3)); 
    
  Middle.setElement(ONE_c,Indices(0,0,0));
  Middle.setElement(ONE_c,Indices(D-1,D-1,0));
  
  if(flavor==(this->N_flavor-1) || flavor<0)
  {
    Last.setElement(factor*ONE_c,Indices(0,0,0));
    Last.setElement(ONE_c,Indices(1,0,0));
    Last.setElement(factor*ONE_c,Indices(0,0,3));
  }
  else
    Last.setElement(ONE_c,Indices(1,0,0));
  
  //Reshape and contract
  First.reshape(Indices(1*D,6)); First.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); First.reshape(Indices(1,D,this->d_fermi,this->d_fermi)); First.permute(Indices(3,1,4,2));
  Op.reshape(Indices(D*D,6)); Op.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Op.reshape(Indices(D,D,this->d_fermi,this->d_fermi)); Op.permute(Indices(3,1,4,2));
  Middle.reshape(Indices(D*D,6)); Middle.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Middle.reshape(Indices(D,D,this->d_fermi,this->d_fermi)); Middle.permute(Indices(3,1,4,2));
  Last.reshape(Indices(D*1,6)); Last.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Last.reshape(Indices(D,1,this->d_fermi,this->d_fermi)); Last.permute(Indices(3,1,4,2));
  
  bool op_saved=false, middle_saved=false;
  int pos_op_saved, pos_middle_saved, site_index;
  
  for(int n=0; n<this->N_sites; n++)
    for(int f=0; f<this->N_flavor; f++)
    {
      site_index=n*this->N_flavor+f;
      if(site_index==0)
	mpo.setOp(site_index,new Operator(First),true);
      else if(site_index==(this->N_total-1))
	mpo.setOp(site_index,new Operator(Last),true);
      else if(f==flavor || flavor<0)
      {
	//Site of the kind of flavor, where I would like to compute the condensate thereof
	if(op_saved)
	  mpo.setOp(site_index,&mpo.getOp(pos_op_saved),false);
	else
	{
	  mpo.setOp(site_index,new Operator(Op),true);
	  op_saved=true;
	  pos_op_saved=site_index;
	}
      }
      else
      {
	//Just some site in between (in principle to the left a matrix [1 0; 0 0] and to the right a matrix [0 0; 1 0] would be sufficient, but as there cannot be a one on the additional "channels" the middle matrix has, that is not a problem and I can spare two matrices
	if(middle_saved)
	  mpo.setOp(site_index,&mpo.getOp(pos_middle_saved),false);
	else
	{
	  mpo.setOp(site_index,new Operator(Middle),true);
	  middle_saved=true;
	  pos_middle_saved=site_index;
	}
      }	
    }  
#ifdef MATOUT
  stringstream name;    
  name.str("");
  if(flavor>=0)
    name << "ParcticleNumberMPO_f" << flavor << ".m";
  else
    name << "getParticleNumberMPO_all_flavor.m";
  mpo.exportForMatlab(name.str().c_str());
#endif
}

void MultiFlavorSchwingerHamiltonian::getChargeMPO(MPO &mpo, int site, int flavor) const
{
  //Some error checking
  if(site<0 || site>=this->N_sites)
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::getChargeMPO(), given site index is invalid" << endl;
    exit(666);
  }
  if(flavor<0 || flavor>=this->N_flavor)
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::getChargeMPO(), given flavor is invalid" << endl;
    exit(666);
  }
  this->getSingleBodyMPO(mpo,site,flavor,chargeOp);  
}

void MultiFlavorSchwingerHamiltonian::getChargeMPO(MPO &mpo) const
{
  //Simply a wrapper
  this->getChargeMPO(mpo,-1);
}

void MultiFlavorSchwingerHamiltonian::getChargeMPO(MPO &mpo,int flavor) const
{
  //Some error checking
  if(flavor>=this->N_flavor)
  {
    cout << "Warning in MultiFlavorSchwingerHamiltonian::getChargeMPO(), specified flavor " << flavor << " larger than the number of flavors, will sum over all flavors" << endl;
    flavor=-1;
  }
  
  mpo.initLength(this->N_total);
  
  int D=2;
  double sign=(((this->N_sites-1)%2)==0)?1.0:-1.0;
    
  mwArray First,Last,Odd,Even,Middle;
  First=mwArray(Indices(1,D,6));	First.fillWithZero();
  Odd=mwArray(Indices(D,D,6));		Odd.fillWithZero();
  Even=mwArray(Indices(D,D,6));		Even.fillWithZero();
  Middle=mwArray(Indices(D,D,6));	Middle.fillWithZero();
  Last=mwArray(Indices(D,1,6));		Last.fillWithZero();
  
  if(flavor==0 || flavor<0)
  {
    First.setElement(ONE_c,Indices(0,0,0));
    First.setElement(0.5*ONE_c,Indices(0,1,0));
    First.setElement(0.5*ONE_c,Indices(0,1,3));
  }
  else
    First.setElement(ONE_c,Indices(0,0,0));
  
  Odd.setElement(ONE_c,Indices(0,0,0));
  Odd.setElement(ONE_c,Indices(D-1,D-1,0));
  Odd.setElement(-0.5*ONE_c,Indices(0,1,0));
  Odd.setElement(0.5*ONE_c,Indices(0,1,3)); 
  
  Even.setElement(ONE_c,Indices(0,0,0));
  Even.setElement(ONE_c,Indices(D-1,D-1,0));
  Even.setElement(0.5*ONE_c,Indices(0,1,0));
  Even.setElement(0.5*ONE_c,Indices(0,1,3));
  
  Middle.setElement(ONE_c,Indices(0,0,0));
  Middle.setElement(ONE_c,Indices(D-1,D-1,0));
  
  if(flavor==(this->N_flavor-1) || flavor<0)
  {
    Last.setElement(sign*0.5*ONE_c,Indices(0,0,0));
    Last.setElement(ONE_c,Indices(1,0,0));
    Last.setElement(0.5*ONE_c,Indices(0,0,3));
  }
  else
    Last.setElement(ONE_c,Indices(1,0,0));
  
  //Reshape and contract
  First.reshape(Indices(1*D,6)); First.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); First.reshape(Indices(1,D,this->d_fermi,this->d_fermi)); First.permute(Indices(3,1,4,2));
  Odd.reshape(Indices(D*D,6)); Odd.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Odd.reshape(Indices(D,D,this->d_fermi,this->d_fermi)); Odd.permute(Indices(3,1,4,2));
  Even.reshape(Indices(D*D,6)); Even.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Even.reshape(Indices(D,D,this->d_fermi,this->d_fermi)); Even.permute(Indices(3,1,4,2));
  Middle.reshape(Indices(D*D,6)); Middle.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Middle.reshape(Indices(D,D,this->d_fermi,this->d_fermi)); Middle.permute(Indices(3,1,4,2));
  Last.reshape(Indices(D*1,6)); Last.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Last.reshape(Indices(D,1,this->d_fermi,this->d_fermi)); Last.permute(Indices(3,1,4,2));
  
  bool odd_saved=false, even_saved=false, middle_saved=false;
  int pos_odd_saved, pos_even_saved, pos_middle_saved, site_index;
  
  for(int n=0; n<this->N_sites; n++)
    for(int f=0; f<this->N_flavor; f++)
    {
      site_index=n*this->N_flavor+f;
      if(site_index==0)
	mpo.setOp(site_index,new Operator(First),true);
      else if(site_index==(this->N_total-1))
	mpo.setOp(site_index,new Operator(Last),true);
      else if(f==flavor || flavor<0)
      {
	//Site of the kind of flavor, where I would like to compute the condensate thereof
	if((n%2)==0)
	{
	  //Even site
	  if(even_saved)
	    mpo.setOp(site_index,&mpo.getOp(pos_even_saved),false);
	  else
	  {
	    mpo.setOp(site_index,new Operator(Even),true);
	    even_saved=true;
	    pos_even_saved=site_index;
	  }
	}
	else
	{
	  //Odd site
	  if(odd_saved)
	    mpo.setOp(site_index,&mpo.getOp(pos_odd_saved),false);
	  else
	  {
	    mpo.setOp(site_index,new Operator(Odd),true);
	    odd_saved=true;
	    pos_odd_saved=site_index;
	  }
	}
      }
      else
      {
	//Just some site in between (in principle to the left a matrix [1 0; 0 0] and to the right a matrix [0 0; 1 0] would be sufficient, but as there cannot be a one on the additional "channels" the middle matrix has, that is not a problem and I can spare two matrices
	if(middle_saved)
	  mpo.setOp(site_index,&mpo.getOp(pos_middle_saved),false);
	else
	{
	  mpo.setOp(site_index,new Operator(Middle),true);
	  middle_saved=true;
	  pos_middle_saved=site_index;
	}
      }	
    }  
#ifdef MATOUT
  stringstream name;    
  name.str("");
  if(flavor>=0)
    name << "ChargeMPO_f" << flavor << ".m";
  else
    name << "ChargeMPO_all_flavor.m";
  mpo.exportForMatlab(name.str().c_str());
#endif
}

void MultiFlavorSchwingerHamiltonian::getChargePenaltyMPO(MPO &mpo, bool const_off) const
{
  mpo.initLength(this->N_total);
  
  int D=4;
  mwArray First,Odd,Even,Last;
  First=mwArray(Indices(1,D,6));	First.fillWithZero();
  Odd=mwArray(Indices(D,D,6));		Odd.fillWithZero();
  Even=mwArray(Indices(D,D,6));		Even.fillWithZero();
  Last=mwArray(Indices(D,1,6));		Last.fillWithZero();
  
  double constant=this->lambda;
  if(const_off)
    constant=1.0;
  
  First.setElement(ONE_c,Indices(0,0,0));
  First.setElement(ONE_c*constant/2.0,Indices(0,D-1,0));  
  First.setElement(ONE_c*constant/2.0,Indices(0,1,0));  
  
  First.setElement(ONE_c*constant/2.0,Indices(0,D-1,3));
  First.setElement(ONE_c*constant/2.0,Indices(0,2,3));
  
  Odd.setElement(ONE_c,Indices(0,0,0));
  Odd.setElement(ONE_c,Indices(1,1,0));
  Odd.setElement(ONE_c,Indices(2,2,0));
  Odd.setElement(ONE_c,Indices(D-1,D-1,0));
  Odd.setElement(ONE_c*constant/2.0,Indices(0,D-1,0));
  Odd.setElement(-ONE_c*constant/2.0,Indices(0,1,0));
  Odd.setElement(-ONE_c,Indices(1,D-1,0));
  Odd.setElement(-ONE_c,Indices(2,D-1,0));
  
  Odd.setElement(ONE_c,Indices(1,D-1,3));
  Odd.setElement(-ONE_c*constant/2.0,Indices(0,D-1,3));
  Odd.setElement(ONE_c*constant/2.0,Indices(0,2,3));
  Odd.setElement(ONE_c,Indices(2,D-1,3));
  
  Even.setElement(ONE_c,Indices(0,0,0));
  Even.setElement(ONE_c,Indices(1,1,0));
  Even.setElement(ONE_c,Indices(2,2,0));
  Even.setElement(ONE_c,Indices(D-1,D-1,0));
  Even.setElement(ONE_c*constant/2.0,Indices(0,D-1,0));
  Even.setElement(ONE_c*constant/2.0,Indices(0,1,0));
  Even.setElement(ONE_c,Indices(1,D-1,0));
  Even.setElement(ONE_c,Indices(2,D-1,0));
  
  Even.setElement(ONE_c,Indices(1,D-1,3));
  Even.setElement(ONE_c*constant/2.0,Indices(0,D-1,3));
  Even.setElement(ONE_c*constant/2.0,Indices(0,2,3));
  Even.setElement(ONE_c,Indices(2,D-1,3));
  
  Last.setElement(ONE_c*constant/2.0,Indices(0,0,0));
  Last.setElement(ONE_c*pow(-1.0,this->N_total-1),Indices(1,0,0));
  Last.setElement(ONE_c*pow(-1.0,this->N_total-1),Indices(2,0,0));
  Last.setElement(ONE_c,Indices(D-1,0,0));
  
  Last.setElement(ONE_c*pow(-1.0,this->N_total-1)*constant/2.0,Indices(0,0,3));
  Last.setElement(ONE_c,Indices(1,0,3));
  Last.setElement(ONE_c,Indices(2,0,3));
  
  First.reshape(Indices(1*D,6)); First.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); First.reshape(Indices(1,D,this->d_fermi,this->d_fermi)); First.permute(Indices(3,1,4,2));
  Odd.reshape(Indices(D*D,6)); Odd.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Odd.reshape(Indices(D,D,this->d_fermi,this->d_fermi)); Odd.permute(Indices(3,1,4,2));
  Even.reshape(Indices(D*D,6)); Even.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Even.reshape(Indices(D,D,this->d_fermi,this->d_fermi)); Even.permute(Indices(3,1,4,2));
  Last.reshape(Indices(D*1,6)); Last.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Last.reshape(Indices(D,1,this->d_fermi,this->d_fermi)); Last.permute(Indices(3,1,4,2));
  
  bool odd_saved = false, even_saved=false;
  int pos_odd_saved, pos_even_saved,site_index;
  
  for(int n=0; n<this->N_sites; n++)
    for(int f=0; f<this->N_flavor; f++)
    {
      site_index=n*this->N_flavor+f;
       if(site_index==0)
	mpo.setOp(site_index,new Operator(First),true);
      else if(site_index==(this->N_total-1))
	mpo.setOp(site_index,new Operator(Last),true);
      else if((n%2)==0)
      {
	//Even site
	if(even_saved)
	  mpo.setOp(site_index,&mpo.getOp(pos_even_saved),false);
	else
	{
	  mpo.setOp(site_index,new Operator(Even),true);
	  even_saved=true;
	  pos_even_saved=site_index;
	}
      }
      else
      {
	//Odd site
	if(odd_saved)
	  mpo.setOp(site_index,&mpo.getOp(pos_odd_saved),false);
	else
	{
	  mpo.setOp(site_index,new Operator(Odd),true);
	  odd_saved=true;
	  pos_odd_saved=site_index;
	}
      }
    }
#ifdef MATOUT
  mpo.exportForMatlab("ChargePenalty.m");
#endif
}

void MultiFlavorSchwingerHamiltonian::getUMPO(MPO &mpo, double dt, bool odd, bool imag_time) const
{
  mpo.initLength(this->N_total);
  complex_t delta_t=-I_c*dt;
  
  //In case i want imaginary time
  if(imag_time)
    delta_t=-I_c*delta_t;
  
  cout << "delta_t=" << delta_t<< endl;
  
  //Construct the necessary terms and sum them
  mwArray Local_op,exp_op;
  //Since I have two body terms and N_f flavors, the dimension of the final operator is
  Local_op = mwArray(Indices(pow(this->d_fermi,2*this->N_flavor),pow(this->d_fermi,2*this->N_flavor)));
  vector<string> term1,term2;
  for(int f=0; f<this->N_flavor; f++)
  {
    mwArray hopping1 = identityMatrix(1);
    mwArray hopping2 = identityMatrix(1);
    term1.clear(); term2.clear();
    for(int k=0; k<2*this->N_flavor; k++)
    {
      if(k==f)
      {
	//Start the hopping term
	hopping1 = kron(hopping1,this->Z.subArray(Indices(1,-1,-1)));
	hopping2 = kron(hopping2,this->Z.subArray(Indices(2,-1,-1)));
	//term1.push_back("sp ");
	//term2.push_back("sm ");
      }
      else if(k==(f+this->N_flavor))
      {
	//Finish the hopping term
	hopping1 = kron(hopping1,this->Z.subArray(Indices(2,-1,-1)));
	hopping2 = kron(hopping2,this->Z.subArray(Indices(1,-1,-1)));
	//term1.push_back("sm ");
	//term2.push_back("sp ");
      }
      else if(k>f && k<(f+this->N_flavor))
      {
	//Sigma_z
	hopping1 = kron(hopping1,I_c*this->Z.subArray(Indices(3,-1,-1)));
	hopping2 = kron(hopping2,-I_c*this->Z.subArray(Indices(3,-1,-1)));
	//term1.push_back("i*s3 ");
	//term2.push_back("-i*s3 ");
      }
      else
      {
	//Identity
	hopping1 = kron(hopping1,this->Z.subArray(Indices(0,-1,-1)));
	hopping2 = kron(hopping2,this->Z.subArray(Indices(0,-1,-1)));
	//term1.push_back("1 ");
	//term2.push_back("1 ");
      }
    }
    /*cout << "Term1: " ;
    for (int bla=0; bla<term1.size(); bla++)
      cout << term1[bla].c_str() ;
    cout << endl << "Term2: " ;
    for (int bla=0; bla<term2.size(); bla++)
      cout << term2[bla].c_str() ;
    cout << endl;*/
    
    Local_op = Local_op - this->x*hopping1 - this->x*hopping2;
  } 
  wrapper::expm(Local_op,exp_op,delta_t);
  //cout << "Local_op " << Local_op << endl;
  //cout << "exp(Local_op) " << exp_op << endl;  
  
  //Now split into single matrices
  vector<int> phys_dims;	phys_dims.resize(2*this->N_flavor,this->d_fermi);
  vector<mwArray> matrices;
  splitLocalOp(exp_op,phys_dims,matrices);
  
  //Now put the matrices
  int start_site=odd?1:0;
  int site_index;
  
  for(int n=start_site; n<(this->N_sites-1); n+=2)
    for(int f=0; f<this->N_flavor; f++)
    {
      site_index=n*this->N_flavor+f;
      mpo.setOp(site_index,new Operator(matrices[f]),true);
      mpo.setOp(site_index+this->N_flavor,new Operator(matrices[f+this->N_flavor]),true);
    }
  //cout << "UMPO " << mpo << endl;
  //Remaining identities
  for(int i=0; i<this->N_total; i++)
  {
    if(mpo.isEmpty(i))
      mpo.setOp(i,new Operator(this->id_fermi_mpo),true);
  }
#ifdef MATOUT
  if(odd)
    mpo.exportForMatlab("Uodd.m");
  else
    mpo.exportForMatlab("Ueven.m");
#endif
}

void MultiFlavorSchwingerHamiltonian::getUoddMPO(MPO &mpo, double dt, bool imag_time) const
{
  this->getUMPO(mpo,dt,true,imag_time);
}
void MultiFlavorSchwingerHamiltonian::getUevenMPO(MPO &mpo, double dt, bool imag_time) const
{
  this->getUMPO(mpo,dt,false,imag_time);
}

void MultiFlavorSchwingerHamiltonian::getUzMPO(MPO &mpo, double dt, bool imag_time,int lmax) const
{
  mpo.initLength(this->N_total);
  complex_t delta_t=-I_c*dt;
  
  //See if I have a cutoff
  if(lmax!=0)
    cout << "Cutoff at lmax=" << lmax << endl;
  
  //In case i want imaginary time
  if(imag_time)
    delta_t=-I_c*delta_t;
  
  cout << "delta_t=" << delta_t<< endl;
  
  
  mwArray M;
  int Dl=0,Dr=1;
  double munuterm;
  int site_index,charge,lval,llmax=this->l0,llmin=this->l0,lprev,index;
  
  vector <int> field_values, charge_values,lprev_values;
  field_values.clear();
  vector<int>::iterator it;
  double constant;
  
  for(int n=0; n<this->N_sites; n++)
    for(int f=0; f<this->N_flavor; f++)
    {
      site_index=n*this->N_flavor+f;
      //cout << "Working on n=" << n << ", f="<< f << ", site=" << site_index << endl;
      //Adjust bond dimension according to site where I look at
      Dl=Dr;
      Dr++;      
      //Once I am at the last flavor at the second last site, I do not have to pass on information anymore
      if(site_index>=(this->N_total-this->N_flavor-1))
	Dr=1;
      //If I am at the last site, I only have 1x1 matrices, so I can set Dl to 1
      if(site_index>(this->N_total-this->N_flavor-1))
	Dl=1;
      M=mwArray(Indices(this->d_fermi,Dl,this->d_fermi,Dr));
      M.fillWithZero();
      
      charge_values.clear();
      field_values.clear();
      lprev_values.clear();
      
      for(int dl=0; dl<Dl; dl++)
	for(int d=0; d<this->d_fermi; d++)
	{
	  //Mass and chemical potential
	  munuterm = d==0?(this->mu[f]*pow(-1.0,n)+this->nu[f]):0.0;
	  //Charge
	  charge=(d==0)?(1+pow(-1,n))/2:(-1+pow(-1,n))/2;
	  //Previous L-value
	  lprev=dl+llmin;
	  //New L-value (previous one plus charge)
	  lval=lprev+charge;
	  field_values.push_back(lval);
	  charge_values.push_back(charge);
	  lprev_values.push_back(lprev);
	  //Transform the new L-value in an index such that the fluxes are orderd increasingly from lmin to lmax as dr goes from 1 to Dr
	  if((n%2)==1)
	    index = lval - (llmin-1);
	  else
	    index = lval - llmin;	  
	  //cout << "Index=" << index << endl;
	  //cout << "dl=" << dl << " d=" << d << " lprev=" << lprev << " lval=" << lval << " charge=" << charge << " llmin=" << llmin << " n=" << n << endl;
	  //Set entries accordingly, I put the contribution from the electric field in the last spin on a site corresponding to flavor N_flavor
	  if(f==(this->N_flavor-1) && n<(this->N_sites-2))
	    M.setElement(exp(delta_t*munuterm)*exp(delta_t*lval*lval),Indices(d,dl,d,index));
	  else if(f==(this->N_flavor-1) && n==(this->N_sites-2))
	    M.setElement(exp(delta_t*munuterm)*exp(delta_t*lval*lval),Indices(d,dl,d,0));
	  else if(n<(this->N_sites-1))
	    M.setElement(exp(delta_t*munuterm),Indices(d,dl,d,index));
	  else
	    M.setElement(exp(delta_t*munuterm),Indices(d,0,d,0));  
	}  
      sort(field_values.begin(),field_values.end());
      sort(charge_values.begin(),charge_values.end());
      sort(lprev_values.begin(),lprev_values.end());
      field_values.erase(unique(field_values.begin(),field_values.end()),field_values.end());
      charge_values.erase(unique(charge_values.begin(),charge_values.end()),charge_values.end());
      lprev_values.erase(unique(lprev_values.begin(),lprev_values.end()),lprev_values.end());
      /*cout << "Previous field values n=" << n << ", f=" <<f<<": " << lprev_values << endl;
      cout << "Charge values n=" << n << ", f=" <<f<<": " << charge_values << endl ;
      cout << "Field values n=" << n << ", f=" <<f<<": " << field_values << endl<< endl;      */
      //Adjust current minimal and maximal value for the flux (odd sites only increase the flux, even sites only decrease)
      if((n%2)!=0)
	llmin--;
      else
	llmax++;
      
      mpo.setOp(site_index,new Operator(M),true);
    }
  //cout << "MPO: " << mpo << endl;
  
  //Now apply compression, if desired
  //TODO: Fix this part, compression is not yet working
  if(lmax!=0)
  {
    vector <int> bond_dims, phys_dims;
    int dl,dr;
    mwArray M;
    int Dmax = 2*lmax+1;
    for(int i=0; i<this->N_total; i++)
    {
      M = mpo.getOp(i).getFullData();
      phys_dims.push_back(M.getDimension(0)*M.getDimension(2));
      if(i<(this->N_total-1))
	bond_dims.push_back(M.getDimension(3));
    }
    
    MPS mps(this->N_total,bond_dims,phys_dims);
    //cout << "MPS = " << mps << endl;
    for(int i=0; i<this->N_total; i++)
    {
      M=mpo.getOp(i).getFullData();
      //cout << "M=(" << M.getDimension(0) << "," << M.getDimension(1) << "," << M.getDimension(2) << "," << M.getDimension(3) << ")" << endl;
      dl = M.getDimension(1);
      dr = M.getDimension(3);
      M.permute(Indices(1,3,2,4));
      //cout << "M=(" << M.getDimension(0) << "," << M.getDimension(1) << "," << M.getDimension(2) << "," << M.getDimension(3) << ")" << endl;
      M.reshape(Indices(this->d_fermi*this->d_fermi,dl,dr));
      mps.setA(i,M); 
      //cout << "M=(" << M.getDimension(0) << "," << M.getDimension(1) << "," << M.getDimension(2) << ")" << endl << endl;
    }
    //Generate a contractor
    Contractor& contractor=Contractor::theContractor();
    cout<<"Initialized Contractor"<<endl; 
    
    //In case PRIMME-support is possible set eigensolver to PRIMME
  #ifdef USING_PRIMME
    cout << "Setting eigensolver to PRIMME" << endl;
    contractor.setEigenSolver(primme);
  #endif
    contractor.setConvTol(1E-8);
    MPS res(mps);
    //contractor.optimizeMPS(mps,res,Dmax);
    
    //cout << "Overlap <res|mps>=" << contractor.contract(res,mps)/sqrt(contractor.contract(mps,mps)*contractor.contract(res,res)) << endl;
    
    for(int i=0; i<this->N_total; i++)
    {
      M=res.getA(i).getA();
      //cout << "M=" << M << endl;
      dl = M.getDimension(1);
      dr = M.getDimension(2);
      M.reshape(Indices(this->d_fermi,this->d_fermi,dl,dr));
      M.permute(Indices(1,3,2,4));
      mpo.setOp(i,new Operator(M),true);
    }
  }
#ifdef MATOUT
  //cout << "MPO after compression: " << mpo <<endl;
  mpo.exportForMatlab("Uzmpo.m");
#endif
}

void MultiFlavorSchwingerHamiltonian::constructInitialMPS(MPS &mps, InitState state) const
{
  //Adjust the mps  
  mps.clear();
  mps=MPS(this->N_total,1,2);
  
  int site_index;
  
  mwArray spinsiteup(Indices(2,1,1));
  mwArray spinsitedown(Indices(2,1,1));
  
  spinsiteup.fillWithZero();
  spinsitedown.fillWithZero();
  
  spinsiteup.setElement(ONE_c,Indices(0,0,0));
  spinsitedown.setElement(ONE_c,Indices(1,0,0));
  
  //Set the matrices
  int index;
  if(state==is_strong_coupling)
  {
    for(int n=0; n<this->N_sites; n++)
      for(int f=0; f<this->N_flavor; f++)
      {
	site_index=n*this->N_flavor+f;
	if((n%2)==0)
	  mps.setA(site_index,spinsitedown);
	else
	  mps.setA(site_index,spinsiteup);
      } 
  }
  else if(state==is_anti_strong_coupling)
  {
    for(int n=0; n<this->N_sites; n++)
      for(int f=0; f<this->N_flavor; f++)
      {
	site_index=n*this->N_flavor+f;
	if((n%2)==0)
	  mps.setA(site_index,spinsiteup);
	else
	  mps.setA(site_index,spinsitedown);
      } 
  }
  else if(state==is_downup)
  {
    for(int k=0; k<this->N_total; k++)
    {
      if((k%2)==0)
	mps.setA(k,spinsitedown);
      else
	mps.setA(k,spinsiteup);      
    } 
  }
  else if(state==is_updown)
  {
    for(int k=0; k<this->N_total; k++)
    {
      if((k%2)==0)
	mps.setA(k,spinsiteup);
      else
	mps.setA(k,spinsitedown);      
    } 
  }
  else if(state==is_down)
  {
    for(int k=0; k<this->N_total; k++)
      mps.setA(k,spinsitedown);
  }
  else if(state==is_up)
  {
    for(int k=0; k<this->N_total; k++)
      mps.setA(k,spinsiteup);
  }
  //Apply gauge condition and normalize
  mps.gaugeCond('R',true);
}

void MultiFlavorSchwingerHamiltonian::getParticleNumberProjector(MPO &mpo, int flavor, int Nf) const
{
  //Some error checking
  if(flavor<0 || flavor>=this->N_flavor)
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::getParticleNumberProjector(), given flavor must be in between 0 and " << this->N_flavor-1 << endl;
    exit(666);
  }
  if(Nf<0 || Nf>=this->N_sites)
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::getParticleNumberProjector(), given particle number must be in between 0 and " << this->N_sites << endl;
    exit(666);
  }
  mpo.initLength(this->N_total);
  
  //Determine the central site
  int Nmiddle=round(this->N_sites/2.0);
  int middle_index = Nmiddle*this->N_flavor+flavor;
  //cout << "Middle site: " << Nmiddle << endl;
  //cout << "Middle index: " << middle_index << endl;
  
  mwArray M;
  int Dl=1,Dr=1,site_index;
  //Now put the matrices to the left  
  for(int n=0; n<=Nmiddle; n++)
  {
    for(int f=0; f<this->N_flavor; f++)
    {
      //Prevent from going to the flavors on the same site which are right of the middle site
      site_index=n*this->N_flavor+f;
      if(site_index>=middle_index)
	break;
      
      if(f==flavor)
      {
	//Here I have to do something
	M=mwArray(Indices(this->d_fermi,Dl,this->d_fermi,Dl+1));
	M.fillWithZero();
	for(int l=0; l<Dl; l++)
	{
	  M.setElement(ONE_c,Indices(0,l,0,l+1));
	  M.setElement(ONE_c,Indices(1,l,1,l));
	}
	Dl++;
      }
      else
      {
	//Just pass on
	M=identityMatrix(this->d_fermi*Dl);
	M.reshape(Indices(this->d_fermi,Dl,this->d_fermi,Dl));
	//for(int l=0; l<this->d_fermi; l++)
	  //cout << "M(" << l << ",:," << l << ",:)=" << M.subArray(Indices(l,-1,l,-1)) << endl;
      }      
      mpo.setOp(site_index,new Operator(M),true);
    }
  }
  //cout << "mpo = " << mpo << endl;
  
  //Now put the matrices to the left  
  for(int n=this->N_sites-1; n>=Nmiddle; n--)
  {
    for(int f=this->N_flavor-1; f>=0; f--)
    {
      //Prevent from going to the flavors on the same site which are left of the middle site
      site_index=n*this->N_flavor+f;
      if(site_index<=middle_index)	
	break;
      
      if(f==flavor)      
      {
	//Here I have to do something
	M=mwArray(Indices(this->d_fermi,Dr+1,this->d_fermi,Dr));
	M.fillWithZero();
	for(int l=0; l<Dr; l++)
	{
	  M.setElement(ONE_c,Indices(0,l+1,0,l));
	  M.setElement(ONE_c,Indices(1,l,1,l));
	}
	Dr++;
      }
      else
      {
	//Just pass on
	M=identityMatrix(this->d_fermi*Dr);
	M.reshape(Indices(this->d_fermi,Dr,this->d_fermi,Dr));
	//for(int l=0; l<this->d_fermi; l++)
	  //cout << "M(" << l << ",:," << l << ",:)=" << M.subArray(Indices(l,-1,l,-1)) << endl;
      }
      //cout << "Setting site: " << site_index << endl;
      mpo.setOp(site_index,new Operator(M),true);
    }
  }
  
  //cout << "mpo: " << mpo << endl;
  
  
  //Center size  
  M=mwArray(Indices(this->d_fermi,Dl,this->d_fermi,Dr));
  M.fillWithZero();
  for (int i=0; i<Dl; i++)
    for (int k=0; k<Dr; k++)
    {
      if((i+k+1)==Nf)
	M.setElement(ONE_c,Indices(0,i,0,k));
      else if((i+k)==Nf)
	M.setElement(ONE_c,Indices(1,i,1,k));
    }  
  mpo.setOp(middle_index,new Operator(M),true);
  
#ifdef MATOUT
  stringstream name;    
  name.str("");
  name << "Projector_f" << flavor << "N" << Nf << ".m";
  mpo.exportForMatlab(name.str().c_str());
#endif
}

void MultiFlavorSchwingerHamiltonian::getParticlePenaltyMPO(MPO &mpo, int flavor, int Nf,bool const_off) const
{
  //Some error checking
  if(flavor<0 || flavor>=this->N_flavor)
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::getParticleNumberProjector(), given flavor must be in between 0 and " << this->N_flavor-1 << endl;
    exit(666);
  }
  if(Nf<0 || Nf>this->N_sites)
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::getParticleNumberProjector(), given particle number must be in between 0 and " << this->N_sites << endl;
    exit(666);
  }
  mpo.initLength(this->N_total);
  
  double constant=this->eta;
  if(const_off)
    constant=1.0;
  
  mwArray First,Last,Site_f,Site_non_f;
  int D=3;
  
  First = mwArray(Indices(1,D,6));	First.fillWithZero();
  Site_f = mwArray(Indices(D,D,6));	Site_f.fillWithZero();
  Site_non_f = mwArray(Indices(D,D,6));	Site_non_f.fillWithZero();
  Last = mwArray(Indices(D,1,6));	Last.fillWithZero();
  
  First.setElement(ONE_c,Indices(0,0,0));
  if(flavor==0)
  {
    //id
    First.setElement(constant*(0.25*this->N_sites + pow(Nf,2.0)/this->N_sites-Nf+0.25)*ONE_c,Indices(0,D-1,0));
    
    //s3
    First.setElement(ONE_c,Indices(0,1,3));
    First.setElement(constant*(0.5*this->N_sites-Nf)*ONE_c,Indices(0,D-1,3));
  }
  
  Site_non_f.setElement(ONE_c,Indices(0,0,0));
  Site_non_f.setElement(ONE_c,Indices(1,1,0));
  Site_non_f.setElement(ONE_c,Indices(D-1,D-1,0));
  
  Site_f = Site_non_f;
  Site_f.setElement(constant*(0.25*this->N_sites + pow(Nf,2.0)/this->N_sites-Nf+0.25)*ONE_c,Indices(0,D-1,0));
  Site_f.setElement(constant*(0.5*this->N_sites-Nf)*ONE_c,Indices(0,D-1,3));
  Site_f.setElement(ONE_c,Indices(0,1,3));
  Site_f.setElement(0.5*constant*ONE_c,Indices(1,D-1,3));
  
  Last.setElement(ONE_c,Indices(D-1,0,0));
  if(flavor==(this->N_flavor-1))
  {
    Last.setElement(constant*(0.25*this->N_sites + pow(Nf,2.0)/this->N_sites-Nf+0.25)*ONE_c,Indices(0,0,0));
    Last.setElement(constant*(0.5*this->N_sites-Nf)*ONE_c,Indices(0,0,3));    
    Last.setElement(0.5*constant*ONE_c,Indices(1,0,3));    
  }
  
  //Contract and reshape
  First.reshape(Indices(1*D,6)); First.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); First.reshape(Indices(1,D,this->d_fermi,this->d_fermi)); First.permute(Indices(3,1,4,2));
  Site_f.reshape(Indices(D*D,6)); Site_f.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Site_f.reshape(Indices(D,D,this->d_fermi,this->d_fermi)); Site_f.permute(Indices(3,1,4,2));
  Site_non_f.reshape(Indices(D*D,6)); Site_non_f.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Site_non_f.reshape(Indices(D,D,this->d_fermi,this->d_fermi)); Site_non_f.permute(Indices(3,1,4,2));
  Last.reshape(Indices(D*1,6)); Last.multiplyRight(reshape(this->Z,Indices(6,this->d_fermi*this->d_fermi))); Last.reshape(Indices(D,1,this->d_fermi,this->d_fermi)); Last.permute(Indices(3,1,4,2));
  
  mpo.setOp(0,new Operator(First),true);
  mpo.setOp(this->N_total-1,new Operator(Last),true);
  
  int saved=false;
  int pos_saved;
  for(int i=flavor; i<this->N_total-1; i+=this->N_flavor)
  {
    if(saved && i!=0)
      mpo.setOp(i,&mpo.getOp(pos_saved),false);
    else if(i!=0)
    {
      mpo.setOp(i,new Operator(Site_f),true);
      pos_saved=i;
      saved=true;
    }
  }
  
  saved=false;  
  for(int i=0; i<this->N_total-1; i++)
  {
    if(mpo.isEmpty(i))
    {
      if(saved)
	mpo.setOp(i,&mpo.getOp(pos_saved),false);
      else
      {
	mpo.setOp(i,new Operator(Site_non_f),true);
	pos_saved=i;
	saved=true;
      }
    }
  }  
}

void MultiFlavorSchwingerHamiltonian::getLocalUMPO(MPO &mpo, double dt, bool imag_time, int site) const
{
  mpo.initLength(this->N_total);
  complex_t delta_t=-I_c*dt;
  
  //In case i want imaginary time
  if(imag_time)
    delta_t=-I_c*delta_t;
  
  if(site<0)
    site = (int) (round(this->N_sites/2.0)-1.0);
  
  if(site<0 || site>=(this->N_sites-1))
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::getLocalUMPO(), specified site has to be in between 0 and " << this->N_sites-2 << " instead received " << site << endl;
  } 
  
  cout << "delta_t=" << delta_t<< endl;
  
  //Construct the necessary terms and sum them
  mwArray Local_op,exp_op;
  //Since I have two body terms and N_f flavors, the dimension of the final operator is
  Local_op = mwArray(Indices(pow(this->d_fermi,2*this->N_flavor),pow(this->d_fermi,2*this->N_flavor)));
  vector<string> term1,term2;
  mwArray munuterm_left,munuterm_right;
  double prefactor_left=(site%2)==0?1.0:-1.0;
  double prefactor_right=-1.0*prefactor_left;
  int ind1,ind2;
  for(int f=0; f<this->N_flavor; f++)
  {
    mwArray hopping1 = identityMatrix(1);
    mwArray hopping2 = identityMatrix(1);    
    term1.clear(); term2.clear();
    for(int k=0; k<2*this->N_flavor; k++)
    {
      if(k==f)
      {
	//Start the hopping term
	hopping1 = kron(hopping1,this->Z.subArray(Indices(1,-1,-1)));
	hopping2 = kron(hopping2,this->Z.subArray(Indices(2,-1,-1)));
	//term1.push_back("sp ");
	//term2.push_back("sm ");
      }
      else if(k==(f+this->N_flavor))
      {
	//Finish the hopping term
	hopping1 = kron(hopping1,this->Z.subArray(Indices(2,-1,-1)));
	hopping2 = kron(hopping2,this->Z.subArray(Indices(1,-1,-1)));
	//term1.push_back("sm ");
	//term2.push_back("sp ");
      }
      else if(k>f && k<(f+this->N_flavor))
      {
	//Sigma_z
	hopping1 = kron(hopping1,I_c*this->Z.subArray(Indices(3,-1,-1)));
	hopping2 = kron(hopping2,-I_c*this->Z.subArray(Indices(3,-1,-1)));
	//term1.push_back("i*s3 ");
	//term2.push_back("-i*s3 ");
      }
      else
      {
	//Identity
	hopping1 = kron(hopping1,this->Z.subArray(Indices(0,-1,-1)));
	hopping2 = kron(hopping2,this->Z.subArray(Indices(0,-1,-1)));
	//term1.push_back("1 ");
	//term2.push_back("1 ");
      }
    }
    /*cout << "Term1: " ;
    for (int bla=0; bla<term1.size(); bla++)
      cout << term1[bla].c_str() ;
    cout << endl << "Term2: " ;
    for (int bla=0; bla<term2.size(); bla++)
      cout << term2[bla].c_str() ;
    cout << endl;*/
    //Now add mass and chemical potential term
    
    ind1 = pow(this->d_fermi,f);
    ind2 = pow(this->d_fermi,2*this->N_flavor-f-1);
    //cout << "ind1= " << ind1 << ", ind2=" << ind2 << endl;
    
    munuterm_left = kron(identityMatrix(ind1),kron(0.5*(prefactor_left*this->mu[f]+this->nu[f])*(this->Z.subArray(Indices(0,-1,-1))+this->Z.subArray(Indices(3,-1,-1))),identityMatrix(ind2)));
    munuterm_right = kron(identityMatrix(ind2),kron(0.5*(prefactor_right*this->mu[f]+this->nu[f])*(this->Z.subArray(Indices(0,-1,-1))+this->Z.subArray(Indices(3,-1,-1))),identityMatrix(ind1)));    
    
    Local_op = Local_op - this->x*hopping1 - this->x*hopping2 + munuterm_left + munuterm_right;
  } 
  
  wrapper::expm(Local_op,exp_op,delta_t);
  //cout << "Local_op " << Local_op << endl;
  //cout << "exp(Local_op) " << exp_op << endl;  
  
  //Now split into single matrices
  vector<int> phys_dims;	phys_dims.resize(2*this->N_flavor,this->d_fermi);
  vector<mwArray> matrices;
  splitLocalOp(exp_op,phys_dims,matrices);
  
  //Now put the matrices in the mpo
  int site_index;
#ifdef MATOUT
  MPO corempo(2*this->N_flavor);
#endif 
  for(int f=0; f<this->N_flavor; f++)
  {
    site_index=site*this->N_flavor+f;
    mpo.setOp(site_index,new Operator(matrices[f]),true);
    mpo.setOp(site_index+this->N_flavor,new Operator(matrices[f+this->N_flavor]),true);
#ifdef MATOUT
    corempo.setOp(f,new Operator(matrices[f]),true);
    corempo.setOp(f+this->N_flavor,new Operator(matrices[f+this->N_flavor]),true);
#endif
  }
  //cout << "ULocalMPO " << mpo << endl;
  //Remaining identities
  bool saved=false;
  int pos_saved;
  for(int i=0; i<this->N_total; i++)
  {
    if(mpo.isEmpty(i))
    {
      if(saved)
	mpo.setOp(i,&mpo.getOp(pos_saved),false);
      else
      {
	mpo.setOp(i,new Operator(this->id_fermi_mpo),true);
	pos_saved=i;
	saved=true;
      }
    }
  }
#ifdef MATOUT
  corempo.exportForMatlab("Ulocal.m",15);
#endif
}

void MultiFlavorSchwingerHamiltonian::getCorrelationMPO(MPO &mpo, int n1, int f1, int n2, int f2, bool condensate) const
{
  //Some error checking
  if(n1<0 || n1>=this->N_sites || f1<0 || f1>=this->N_flavor || n2<0 || n2>=this->N_sites || f2<0 || f2>=this->N_flavor)
  {
    cout << "Error in MultiFlavorSchwingerHamiltonian::getCorrelationMPO(), invalid flavor/site number" << endl;
    exit(666);    
  }
  
  mpo.initLength(this->N_total);
  
  int ind1 = n1*this->N_flavor+f1;
  int ind2 = n2*this->N_flavor+f2;
  
  mwArray Op1,Op2;
  
  if(condensate)
  {
    cout << "Constructing condensate correlator <O_" << ind1 << " O_" << ind2 << ">" << endl;
    Op1 = pow(-1.0,n1)*0.5*(this->Z.subArray(Indices(0,-1,-1))+this->Z.subArray(Indices(3,-1,-1)));
    if(ind1==ind2)
    {
      Op1 = Op1*Op1;
      Op1.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
    }
    else
    {
      Op2 = pow(-1.0,n2)*0.5*(this->Z.subArray(Indices(0,-1,-1))+this->Z.subArray(Indices(3,-1,-1)));
      Op1.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
      Op2.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
    } 
  }
  else
  {
    cout << "Constructing spin correlator <O_" << ind1 << " O_" << ind2 << ">" << endl;
    Op1 = this->Z.subArray(Indices(3,-1,-1));
    if(ind1==ind2)
    {
      Op1 = Op1*Op1;
      Op1.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
    }
    else
    {
      Op2 = Op1;
      Op1.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
      Op2.reshape(Indices(this->d_fermi,1,this->d_fermi,1));
    } 
  }
  
  mpo.setOp(ind1,new Operator(Op1),true);
  if(ind1!=ind2)
    mpo.setOp(ind2,new Operator(Op2),true);
  
  bool id_saved=false;
  int pos_id_saved;
  for(int i=0; i<this->N_total; i++)
    if(mpo.isEmpty(i))
    {
      //cout << "Setting identity at site " << i << endl;
      if(id_saved)
	mpo.setOp(i,&mpo.getOp(pos_id_saved),false);
      else
      {
	mpo.setOp(i,new Operator(this->id_fermi_mpo),true);
	id_saved=true;
	pos_id_saved=i;
      }
    }
}

