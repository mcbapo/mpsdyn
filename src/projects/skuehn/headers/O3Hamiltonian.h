#ifndef O3_HAMILTONIAN_H
#define O3_HAMILTONIAN_H

#include "MPS.h"
#include "Hamiltonian.h"

/**
   \file O3Hamiltonian.h Definition of the class that implements the truncated O(3) Hamiltonian
   
   \author Stefan Kühn
   \date 03/30/2018
*/

class O3Hamiltonian: 
public Hamiltonian{
private:
  //Array to store all the operators
  mwArray Operators; 
  //Array for an identity MPO
  mwArray id_mpo;
  //Number of different operators I am storing at the init stage
  int nops;
  
  double Wigner3j(double j1, double m1, double j2, double m2, double j3, double m3);
  double factorial(double n);
  
protected:
  //Physical dimension of a site
  int d_site;
  //Number of sites
  int N_sites;
  //Parameters for truncation 
  int lmax,mmax;
  //Constants in the Hamiltonian
  double c1,c2,c3,lambda,eta;
  //Hamiltonian MPO
  MPO hamiltonian;
  //Charge sector 
  int q;
  
   /** Generate some matrices used for single and few body MPOs to be saved in the object, therefore in the MPO the pointers to the matrices can be stored and do not consume memory*/
  void initOperators(void);
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously)*/
  void initHamiltonian(bool target_charge); 
  
  /** Basic functionality to get a single site MPO acting on site n*/
  void getSingleSiteMPO(MPO &mpo,const mwArray &Op, int n) const;
  
 public:
  /** Creates an instance of the truncated O(3) Hamiltonian with chemical potential for a chain of length N_sites sites
  \f[ H = c_1 \sum_k \mathbf{J}^2_k + c_2\sum_k J^z_k + c_3\sum_k \mathbf{n}_k \mathbf{n}_{k+1}. \f]
  If \f$ \lambda\neq 0 \f$ is given, a penalty for the last sector is added which energetically penalizes states with an angular momentum of \f$ l_\mathrm{max} \f$.
   */
  O3Hamiltonian(int N_sites,int lmax, int mmax,double c1,double c2,double c3,double lambda=0.);
  /** Creates an instance of the truncated O(3) Hamiltonian with chemical potential for a chain of length N_sites sites
  \f[ H = c_1 \sum_k \mathbf{J}^2_k + c_2\sum_k J^z_k + c_3\sum_k \mathbf{n}_k \mathbf{n}_{k+1}. \f]
  If \f$ \lambda\neq 0 \f$ is given, a penalty for the last sector is added which energetically penalizes states with an angular momentum of \f$ l_\mathrm{max} \f$. Additionally the penalty term 
  \f[ \eta (\sum_k J^z_k - q)^2 \f]
  which allows for targeting a sector of a certain total charge \f$ q \f$.
   */
  O3Hamiltonian(int N_sites,int lmax, int mmax,double c1,double c2,double c3,double lambda,double eta,int q);
  
  ~O3Hamiltonian();

  /** Get the number of sites */
  int getN() const{return N_sites;}
  
  bool hasHMPO() const {return true;}
  /** Return a reference to the MPO for the Hamiltonian. */
  const MPO& getHMPO() const{return hamiltonian;}
  
  /** Get the MPO to compute J^2 locally at site n*/
  void getJsqMPO(MPO &mpo, int n) const;
  
  /** Get the MPO to compute J^z locally at site n*/
  void getJzMPO(MPO &mpo, int n) const;
  
  /** Get the MPO for the potential term in the Hamiltonian corresponding to \f$ \mathbf{n}_k \mathbf{n}_{k+1} \f$ */
  void getPotentialMPO(MPO &mpo) const;
  
  /** Get the MPO for the total charge \f$ Q=\sum_n J^z_n\f$ */
  void getChargeMPO(MPO &mpo) const;
  
  /** Get the MPO for the total charge squared \f$ Q^2=\left(\sum_n J^z_n\right)^2 = \sum_n (J^z_n)^2 + 2\sum_n \sum_{k>n} J^z_n J^z_k \f$ */
  void getChargeSquareMPO(MPO &mpo) const; 
  
};


#endif // O3_HAMILTONIAN_H
