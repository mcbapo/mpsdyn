#include "O3Hamiltonian.h"
#include "Indices.h"
#include "Contractor.h"
#include <numeric>
#include <algorithm>

using namespace std;
using namespace shrt;


O3Hamiltonian::O3Hamiltonian(int N_sites_,int lmax_, int mmax_,double c1_,double c2_,double c3_,double lambda_):hamiltonian(N_sites_)
{
  //Number of sites
  this->N_sites = N_sites_;  
  //Constants in the Hamiltonian
  this->c1=c1_;
  this->c2=c2_;
  this->c3=c3_;
  this->lambda=lambda_;
  //Truncation parameters
  this->lmax=lmax_;
  this->mmax=mmax_;
  
  initOperators();
  initHamiltonian(false);
}

O3Hamiltonian::O3Hamiltonian(int N_sites_,int lmax_, int mmax_,double c1_,double c2_,double c3_,double lambda_,double eta_,int q_):hamiltonian(N_sites_)
{
  //Number of sites
  this->N_sites = N_sites_;  
  //Constants in the Hamiltonian
  this->c1=c1_;
  this->c2=c2_;
  this->c3=c3_;
  this->lambda=lambda_;
  this->eta=eta_;
  //Truncation parameters
  this->lmax=lmax_;
  this->mmax=mmax_;
  //Charge sector
  this->q = q_;
  
  initOperators();
  initHamiltonian(true);
}

O3Hamiltonian::~O3Hamiltonian()
{
  this->Operators.clear();
}

void O3Hamiltonian::initOperators()
{  
  int count=0;
  vector<int> lvals,mvals;
  //Determine dimension and
  for(int l=0; l<=this->lmax; l++)
  {
    for(int m=max(-l,-this->mmax); m<=min(l,this->mmax); m++)
    {
      count++;
      lvals.push_back(l);
      mvals.push_back(m);
    }
  }
  this->d_site = count;
  //Current number of operators I am dealing with
  this->nops = 8;
  this->Operators=mwArray(Indices(this->nops,this->d_site,this->d_site));  
  
  //Identity Matrix
  mwArray id=identityMatrix(this->d_site);  
  mwArray nx,ny,nz,Jsq,Jz,Jpen,np,nm;
  
  nx = mwArray(Indices(this->d_site,this->d_site)); nx.fillWithZero();
  ny = mwArray(Indices(this->d_site,this->d_site)); ny.fillWithZero();
  nz = mwArray(Indices(this->d_site,this->d_site)); nz.fillWithZero();
  np = mwArray(Indices(this->d_site,this->d_site)); np.fillWithZero();
  nm = mwArray(Indices(this->d_site,this->d_site)); nm.fillWithZero();
  Jsq = mwArray(Indices(this->d_site,this->d_site)); Jsq.fillWithZero();
  Jz = mwArray(Indices(this->d_site,this->d_site)); Jz.fillWithZero();
  Jpen = mwArray(Indices(this->d_site,this->d_site)); Jpen.fillWithZero();
  
  int l,m,ll,mm;
  for(int i=0; i<this->d_site; i++)
  {
    for(int k=0; k<this->d_site; k++)
    {
      
      l = lvals[i];
      m = mvals[i];
      ll = lvals[k];
      mm = mvals[k];
      
      nz.setElement(ONE_c*pow(-1.,m)*sqrt((2.*l+1.)*(2.*ll+1.))*Wigner3j(l,0.,1.,0.,ll,0.)*Wigner3j(l,-m,1.,0.,ll,mm),Indices(i,k));
//       if(abs(nz.getElement(Indices(i,k)))>1E-3)
// 	cout << "nz: setting element <l=" << l <<", m=" << m << "|nz|l'="<< ll <<", m'" << mm << "> = " << nz.getElement(Indices(i,k)) << endl;
      np.setElement(-1.*ONE_c*pow(-1.,m)*sqrt(2.)*sqrt((2.*l+1.)*(2.*ll+1.))*Wigner3j(l,0.,1.,0.,ll,0.)*Wigner3j(l,-m,1.,1.,ll,mm),Indices(i,k));
//       if(abs(np.getElement(Indices(i,k)))>1E-3)
// 	cout << "np: setting element <l=" << l <<", m=" << m << "|np|l'="<< ll <<", m'" << mm << "> = " << np.getElement(Indices(i,k)) << endl;
      nm.setElement(ONE_c*pow(-1.,m)*sqrt(2.)*sqrt((2.*l+1.)*(2.*ll+1.))*Wigner3j(l,0.,1.,0.,ll,0.)*Wigner3j(l,-m,1.,-1.,ll,mm),Indices(i,k));
//       if(abs(nm.getElement(Indices(i,k)))>1E-3)
// 	cout << "nm: setting element <l=" << l <<", m=" << m << "|nm|l'="<< ll <<", m'" << mm << "> = " << nm.getElement(Indices(i,k)) << endl;
      if(i==k)
      {
	  Jz.setElement(ONE_c*m,Indices(i,k));
	  Jsq.setElement(ONE_c*l*(l+1.),Indices(i,k));
	  if(l==this->lmax)
	    Jpen.setElement(ONE_c,Indices(i,k));
      }
         
    }
  }
    
  nx = 0.5*(np+nm);
  ny = -0.5*I_c*(np-nm);
    
//   cout << "nx=" << nx << endl;
//   cout << "ny=" << ny << endl;
//   cout << "nz=" << nz << endl;
//   cout << "Jsq=" << Jsq << endl;
//   cout << "Jz=" << Jz << endl;
//   cout << "Jpen=" << Jpen << endl;
  
  for(int i=0; i<this->d_site; i++)
  {
    for(int k=0; k<this->d_site; k++)
    {    
      this->Operators.setElement(id.getElement(Indices(i,k)),Indices(0,i,k));
      this->Operators.setElement(Jsq.getElement(Indices(i,k)),Indices(1,i,k));
      this->Operators.setElement(Jz.getElement(Indices(i,k)),Indices(2,i,k));
      this->Operators.setElement(nx.getElement(Indices(i,k)),Indices(3,i,k));
      this->Operators.setElement(ny.getElement(Indices(i,k)),Indices(4,i,k));
      this->Operators.setElement(nz.getElement(Indices(i,k)),Indices(5,i,k));      
      this->Operators.setElement(Jpen.getElement(Indices(i,k)),Indices(6,i,k));         
      this->Operators.setElement(Jz.getElement(Indices(i,k))*Jz.getElement(Indices(i,k)),Indices(7,i,k));
    }
  }

  //I need often identities to construct a single body operator, therefore I prepare those and save them
  this->id_mpo=id;
  this->id_mpo.reshape(Indices(this->d_site,1,this->d_site,1));
}

double O3Hamiltonian::factorial(double n)
{
  return (n == 1. || n == 0.) ? 1. : factorial(n - 1.) * n;
}

double O3Hamiltonian::Wigner3j(double j1, double m1, double j2, double m2, double j3, double m3)
{
//   if ( 2*j1 ~= floor(2*j1) || 2*j2 ~= floor(2*j2) || 2*j3 ~= floor(2*j3)|| 2*m1 ~= floor(2*m1) || 2*m2 ~= floor(2*m2) || 2*m3 ~= floor(2*m3) )
//     return = 0.;
//   else if ( j1 - m1 ~= floor ( j1 - m1 ) )
//     return = 0.;
//   else if ( j2 - m2 ~= floor ( j2 - m2 ) )
//     return = 0.;
//   else if ( j3 - m3 ~= floor ( j3 - m3 ) )
//     return = 0.;
  if (j3 > (j1 + j2) || j3 < fabs(j1 - j2))
    return 0.;
  else if (fabs(m1) > j1)
    return 0.;
  else if (fabs(m2) > j2)
    return 0.;
  else if (fabs(m3) > j3)
    return 0.;
  else if (fabs(m1+m2+m3)>1E-3)
    return 0.;
  else
  {
    //Reasonable input, so I have something to compute
    double t1 = j2 - m1 - j3;
    double t2 = j1 + m2 - j3;
    double t3 = j1 + j2 - j3;
    double t4 = j1 - m1;
    double t5 = j2 + m2;
  
    double tmin = max(0.,max( t1, t2));
    double tmax = min(t3,min(t4,t5));

    double wigner = 0.;

    for(int t=tmin; t<=tmax; t++)
        wigner = wigner + pow(-1.,t)/(factorial(t)*factorial(t-t1)*factorial(t-t2)*factorial(t3-t)*factorial(t4-t)*factorial(t5-t));

    wigner *= pow(-1.,j1-j2-m3)* sqrt(factorial(j1+j2-j3)*factorial(j1-j2+j3)*factorial(-j1+j2+j3)/factorial(j1+j2+j3+1.0)*factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2)*factorial(j3+m3)*factorial(j3-m3));
    return wigner;
  }
}

void O3Hamiltonian::initHamiltonian(bool target_charge)
{
  
  cout << "Constructing Hamiltonian with:\n- lmax = " << this->lmax << "\n- mmax = " << this->mmax << "\n- d = " << this->d_site << "\n- c1 = " << this->c1 << "\n- c2 = " << this->c2<< "\n- c3 = " << this->c3 << "\n- la = " << this->lambda <<endl;
  if(target_charge)
    cout << "- eta = " << this->eta << endl <<"- q = " << this->q << endl;
           
  
  //Get the bond dimension
  int D=5;
  if(target_charge)
    D=6;
  
  mwArray Mfirst,M,Mlast;
  
  Mfirst = mwArray(Indices(1,D,nops));		Mfirst.fillWithZero();
  M = mwArray(Indices(D,D,nops));		M.fillWithZero();
  Mlast = mwArray(Indices(D,1,nops));		Mlast.fillWithZero();  
  
  //Identity
  Mfirst.setElement(ONE_c,Indices(0,0,0));
  M.setElement(ONE_c,Indices(0,0,0));
  M.setElement(ONE_c,Indices(D-1,D-1,0));
  Mlast.setElement(ONE_c,Indices(D-1,0,0));
  if(target_charge)
  {
    Mfirst.setElement(this->eta*this->q*this->q/((double)this->N_sites)*ONE_c,Indices(0,D-1,0));
    M.setElement(this->eta*this->q*this->q/((double)this->N_sites)*ONE_c,Indices(0,D-1,0));
    M.setElement(ONE_c,Indices(4,4,0));
    Mlast.setElement(this->eta*this->q*this->q/((double)this->N_sites)*ONE_c,Indices(0,0,0));
  }
  //J^2
  Mfirst.setElement(this->c1*ONE_c,Indices(0,D-1,1));
  M.setElement(this->c1*ONE_c,Indices(0,D-1,1));
  Mlast.setElement(this->c1*ONE_c,Indices(0,0,1));
  //J^z
  Mfirst.setElement(this->c2*ONE_c,Indices(0,D-1,2));
  M.setElement(this->c2*ONE_c,Indices(0,D-1,2));
  Mlast.setElement(this->c2*ONE_c,Indices(0,0,2));
  if(target_charge)
  {
    //I have to take into account that I overwrite the entries, so all what I have set above has to be reset here again
    Mfirst.setElement(this->c2*ONE_c-2.0*this->eta*this->q*ONE_c,Indices(0,D-1,2));
    Mfirst.setElement(2.0*this->eta*ONE_c,Indices(0,4,2));
    M.setElement(this->c2*ONE_c-2.0*this->eta*this->q*ONE_c,Indices(0,D-1,2));
    M.setElement(2.0*this->eta*ONE_c,Indices(0,4,2));
    M.setElement(ONE_c,Indices(4,D-1,2));
    Mlast.setElement(this->c2*ONE_c-2.0*this->eta*this->q*ONE_c,Indices(0,0,2));
    Mlast.setElement(ONE_c,Indices(4,0,2));
  }
  //nx
  Mfirst.setElement(this->c3*ONE_c,Indices(0,1,3));
  M.setElement(this->c3*ONE_c,Indices(0,1,3));
  M.setElement(ONE_c,Indices(1,D-1,3));
  Mlast.setElement(ONE_c,Indices(1,0,3));
  //ny
  Mfirst.setElement(this->c3*ONE_c,Indices(0,2,4));
  M.setElement(this->c3*ONE_c,Indices(0,2,4));
  M.setElement(ONE_c,Indices(2,D-1,4));
  Mlast.setElement(ONE_c,Indices(2,0,4));
  //nz
  Mfirst.setElement(this->c3*ONE_c,Indices(0,3,5));
  M.setElement(this->c3*ONE_c,Indices(0,3,5));
  M.setElement(ONE_c,Indices(3,D-1,5));
  Mlast.setElement(ONE_c,Indices(3,0,5));
  //Jpen
  Mfirst.setElement(this->lambda*ONE_c,Indices(0,D-1,6));
  M.setElement(this->lambda*ONE_c,Indices(0,D-1,6));
  Mlast.setElement(this->lambda*ONE_c,Indices(0,0,6));
  //(J^z)^2
  if(target_charge)
  {
    Mfirst.setElement(this->eta*ONE_c,Indices(0,D-1,7));
    M.setElement(this->eta*ONE_c,Indices(0,D-1,7));
    Mlast.setElement(this->eta*ONE_c,Indices(0,0,7));
  }
  
  Mfirst.reshape(Indices(1*D,nops)); Mfirst.multiplyRight(reshape(this->Operators,Indices(nops,this->d_site*this->d_site))); Mfirst.reshape(Indices(1,D,this->d_site,this->d_site)); Mfirst.permute(Indices(3,1,4,2));
  M.reshape(Indices(D*D,nops)); M.multiplyRight(reshape(this->Operators,Indices(nops,this->d_site*this->d_site)));  M.reshape(Indices(D,D,this->d_site,this->d_site)); M.permute(Indices(3,1,4,2));
  Mlast.reshape(Indices(D*1,nops)); Mlast.multiplyRight(reshape(this->Operators,Indices(nops,this->d_site*this->d_site)));  Mlast.reshape(Indices(D,1,this->d_site,this->d_site)); Mlast.permute(Indices(3,1,4,2));
    
  //First and last matrix
  this->hamiltonian.setOp(0,new Operator(Mfirst),true);
  this->hamiltonian.setOp(this->N_sites-1,new Operator(Mlast),true);
  for(int i=1; i<this->N_sites-1; i++)
  {
    if(i>1)
      this->hamiltonian.setOp(i,&this->hamiltonian.getOp(1),false);
    else
      this->hamiltonian.setOp(i,new Operator(M),true);
  }
  
  //cout << "Hamiltonian" << this->hamiltonian << endl;;
    
#ifdef MATOUT
  stringstream name;    
  name.str("");
  name << "HO3_N" << this->N_sites << "_c1" << this->c1 << "_c2" << this->c2 << "_c3" << this->c3 << "_la" << this->lambda;
  if(target_charge)
    name << "_eta" << this->eta  << "_q" << this->q;
  name << ".m";
    
//   for (int i=0; i<this->N_flavor; i++)
//     name << "mu" << (int) (this->mu[i]*1000);
//   for (int i=0; i<this->N_flavor; i++)
//     name << "nu" << (int) (this->nu[i]*1000);
//   name << "la" << (int) (this->lambda*1000) ;
//   name << "lbg" << (int) (this->l0*1000) << ".m";   
  this->hamiltonian.exportForMatlab(name.str().c_str(),11);
#endif
}


void O3Hamiltonian::getSingleSiteMPO(MPO &mpo, const mwArray &Op, int n) const
{
  //Check if the given site is valid
  if(n<0 || n>=this->N_sites)
  {
    cout << "Error in O3::getSingleSiteMPO(...), given site " << n << " is not a valid site" << endl;
    exit(666);
  }
  //Prepare the MPO 
  mpo.initLength(this->N_sites);
  mpo.setOp(n,new Operator(reshape(Op,Indices(this->d_site,1,this->d_site,1))),true);
  bool saved=false;
  int pos_saved=-1;
  for(int i=0; i<this->N_sites; i++)
  {
    if(i!=n)
    {
      if(saved)
	mpo.setOp(i,&mpo.getOp(pos_saved),false);
      else
      {
	mpo.setOp(i,new Operator(this->id_mpo),true);
	saved=true;
	pos_saved=i;
      }
    }
  }  
}

void O3Hamiltonian::getJsqMPO(MPO &mpo, int n) const
{
  this->getSingleSiteMPO(mpo,this->Operators.subArray(Indices(1,-1,-1)),n);
}
  
void O3Hamiltonian::getJzMPO(MPO &mpo, int n) const
{
  this->getSingleSiteMPO(mpo,this->Operators.subArray(Indices(2,-1,-1)),n);
}

void O3Hamiltonian::getPotentialMPO(MPO &mpo) const
{
  int D=5;
  mpo.initLength(this->N_sites);
  mwArray Mfirst,M,Mlast;
  
  Mfirst = mwArray(Indices(1,D,this->nops));	Mfirst.fillWithZero();
  M = mwArray(Indices(D,D,this->nops));		M.fillWithZero();
  Mlast = mwArray(Indices(D,1,this->nops));	Mlast.fillWithZero();  
  
  //Identity
  Mfirst.setElement(ONE_c,Indices(0,0,0));
  M.setElement(ONE_c,Indices(0,0,0));
  M.setElement(ONE_c,Indices(D-1,D-1,0));
  Mlast.setElement(ONE_c,Indices(D-1,0,0));
  //nx
  Mfirst.setElement(ONE_c,Indices(0,1,3));
  M.setElement(ONE_c,Indices(0,1,3));
  M.setElement(ONE_c,Indices(1,D-1,3));
  Mlast.setElement(ONE_c,Indices(1,0,3));
  //ny
  Mfirst.setElement(ONE_c,Indices(0,2,4));
  M.setElement(ONE_c,Indices(0,2,4));
  M.setElement(ONE_c,Indices(2,D-1,4));
  Mlast.setElement(ONE_c,Indices(2,0,4));
  //nz
  Mfirst.setElement(ONE_c,Indices(0,3,5));
  M.setElement(ONE_c,Indices(0,3,5));
  M.setElement(ONE_c,Indices(3,D-1,5));
  Mlast.setElement(ONE_c,Indices(3,0,5));
  
  Mfirst.reshape(Indices(1*D,this->nops)); Mfirst.multiplyRight(reshape(this->Operators,Indices(this->nops,this->d_site*this->d_site))); Mfirst.reshape(Indices(1,D,this->d_site,this->d_site)); Mfirst.permute(Indices(3,1,4,2));
  M.reshape(Indices(D*D,this->nops)); M.multiplyRight(reshape(this->Operators,Indices(this->nops,this->d_site*this->d_site)));  M.reshape(Indices(D,D,this->d_site,this->d_site)); M.permute(Indices(3,1,4,2));
  Mlast.reshape(Indices(D*1,this->nops)); Mlast.multiplyRight(reshape(this->Operators,Indices(this->nops,this->d_site*this->d_site)));  Mlast.reshape(Indices(D,1,this->d_site,this->d_site)); Mlast.permute(Indices(3,1,4,2));
    
  //First and last matrix
  mpo.setOp(0,new Operator(Mfirst),true);
  mpo.setOp(this->N_sites-1,new Operator(Mlast),true);
  for(int i=1; i<this->N_sites-1; i++)
  {
    if(i>1)
      mpo.setOp(i,&mpo.getOp(1),false);
    else
      mpo.setOp(i,new Operator(M),true);
  } 
}

void O3Hamiltonian::getChargeMPO(MPO &mpo) const
{
  int D=2;
  mpo.initLength(this->N_sites);
  mwArray Mfirst,M,Mlast;
  
  Mfirst = mwArray(Indices(1,D,this->nops));	Mfirst.fillWithZero();
  M = mwArray(Indices(D,D,this->nops));		M.fillWithZero();
  Mlast = mwArray(Indices(D,1,this->nops));	Mlast.fillWithZero();  
  
  //Identity
  Mfirst.setElement(ONE_c,Indices(0,0,0));
  M.setElement(ONE_c,Indices(0,0,0));
  M.setElement(ONE_c,Indices(1,1,0));
  Mlast.setElement(ONE_c,Indices(1,0,0));
  
  //J^z
  Mfirst.setElement(ONE_c,Indices(0,1,2));
  M.setElement(ONE_c,Indices(0,1,2));
  Mlast.setElement(ONE_c,Indices(0,0,2));
  
  Mfirst.reshape(Indices(1*D,this->nops)); Mfirst.multiplyRight(reshape(this->Operators,Indices(this->nops,this->d_site*this->d_site))); Mfirst.reshape(Indices(1,D,this->d_site,this->d_site)); Mfirst.permute(Indices(3,1,4,2));
  M.reshape(Indices(D*D,this->nops)); M.multiplyRight(reshape(this->Operators,Indices(this->nops,this->d_site*this->d_site)));  M.reshape(Indices(D,D,this->d_site,this->d_site)); M.permute(Indices(3,1,4,2));
  Mlast.reshape(Indices(D*1,this->nops)); Mlast.multiplyRight(reshape(this->Operators,Indices(this->nops,this->d_site*this->d_site)));  Mlast.reshape(Indices(D,1,this->d_site,this->d_site)); Mlast.permute(Indices(3,1,4,2));
  
  //First and last matrix
  mpo.setOp(0,new Operator(Mfirst),true);
  mpo.setOp(this->N_sites-1,new Operator(Mlast),true);
  for(int i=1; i<this->N_sites-1; i++)
  {
    if(i>1)
      mpo.setOp(i,&mpo.getOp(1),false);
    else
      mpo.setOp(i,new Operator(M),true);
  } 
#ifdef MATOUT
  stringstream name;    
  name.str("");
  name << "Q_N" << this->N_sites << "_lmax"<< this->lmax << "_mmax" << this->mmax << ".m";
  mpo.exportForMatlab(name.str().c_str(),11);
#endif
}
  
void O3Hamiltonian::getChargeSquareMPO(MPO &mpo) const
{
  int D=3;
  mpo.initLength(this->N_sites);
  mwArray Mfirst,M,Mlast;
  
  Mfirst = mwArray(Indices(1,D,this->nops));	Mfirst.fillWithZero();
  M = mwArray(Indices(D,D,this->nops));		M.fillWithZero();
  Mlast = mwArray(Indices(D,1,this->nops));	Mlast.fillWithZero();  
  
  //Identity
  Mfirst.setElement(ONE_c,Indices(0,0,0));
  M.setElement(ONE_c,Indices(0,0,0));
  M.setElement(ONE_c,Indices(1,1,0));
  M.setElement(ONE_c,Indices(2,2,0));
  Mlast.setElement(ONE_c,Indices(2,0,0));
  
  //J^z
  Mfirst.setElement(2.0*ONE_c,Indices(0,1,2));
  M.setElement(2.0*ONE_c,Indices(0,1,2));
  M.setElement(ONE_c,Indices(1,2,2));
  Mlast.setElement(ONE_c,Indices(1,0,2));
  
  //(J^z)^2
  Mfirst.setElement(ONE_c,Indices(0,2,7));
  M.setElement(ONE_c,Indices(0,D-1,7));
  Mlast.setElement(ONE_c,Indices(0,0,7));
  
  Mfirst.reshape(Indices(1*D,this->nops)); Mfirst.multiplyRight(reshape(this->Operators,Indices(this->nops,this->d_site*this->d_site))); Mfirst.reshape(Indices(1,D,this->d_site,this->d_site)); Mfirst.permute(Indices(3,1,4,2));
  M.reshape(Indices(D*D,this->nops)); M.multiplyRight(reshape(this->Operators,Indices(this->nops,this->d_site*this->d_site)));  M.reshape(Indices(D,D,this->d_site,this->d_site)); M.permute(Indices(3,1,4,2));
  Mlast.reshape(Indices(D*1,this->nops)); Mlast.multiplyRight(reshape(this->Operators,Indices(this->nops,this->d_site*this->d_site)));  Mlast.reshape(Indices(D,1,this->d_site,this->d_site)); Mlast.permute(Indices(3,1,4,2));
  
  //First and last matrix
  mpo.setOp(0,new Operator(Mfirst),true);
  mpo.setOp(this->N_sites-1,new Operator(Mlast),true);
  for(int i=1; i<this->N_sites-1; i++)
  {
    if(i>1)
      mpo.setOp(i,&mpo.getOp(1),false);
    else
      mpo.setOp(i,new Operator(M),true);
  }
#ifdef MATOUT
  stringstream name;    
  name.str("");
  name << "Qsq_N" << this->N_sites << "_lmax"<< this->lmax << "_mmax" << this->mmax << ".m";
  mpo.exportForMatlab(name.str().c_str(),11);
#endif
}

