#ifndef SU2GIS_H
#define SU2GIS_H

#include "MPS.h"
#include "Hamiltonian.h"
#include "Contractor.h"
#include "LGTstring.h"

/**
   \file SU2HamiltonianGIS.h Definition of the class that implements the SU(2) Hamiltonian on the gauge invariant subspace after integrating out the gauge field as developed by Tao and Pablo

   \author Stefan Kühn
   \date 02/23/2018
*/


/** List of possible single body operators*/
enum SingleBodyOperator{
  QxOp,
  QyOp,
  QzOp,
  Q2Op
};

/** List of possible components*/
enum Component{
  x_comp,y_comp,z_comp
};


class SU2HamiltonianGIS: 
public Hamiltonian{
private:
  //Pauli Matrices
  mwArray sigma_x,sigma_y,sigma_z,sigma_plus,sigma_minus;
  //Fermionic idenitity
  mwArray id_fermi;
  //Fermionic identity in MPO form
  mwArray id_site_mpo;
  //Number of sites
  int N_sites;
  //Mass
  double mu;
  //Coupling constant
  double g;
  //Hopping constant 
  double epsilon;
  //Penalty for subspace of vanishing total charge
  double lambda;
  //Penalty for single occupancy of the external charges
  double eta;
  //Indices where the external charges are sitting
  int cext_left,cext_right;
  //Dimension of a site
  int d_site;  
  //Hamiltonian MPO
  MPO hamil;
  
   /** Generate some matrices used for single and few body MPOs to be saved in the object, therefore in the MPO the pointers to the matrices can be stored and do not consume memory*/
  void initOperators(void);
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously)*/
  void initHamiltonian(void);
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously)*/
  void initHamiltonianExtCharges(void);  
  /** Construct the MPO representation of a single body MPO*/
  void getSingleBodyMPO(MPO &mpo,SingleBodyOperator Op, int site) const;
  
  
  public:
  /** Creates an instance of the truncated SU(2) Hamiltonian where the gauge fields are integrated out
 */
  SU2HamiltonianGIS(int N,double epsilon,double mu,double g,double lambda);
  
   /** Creates an instance of the truncated SU(2) Hamiltonian where the gauge fields are integrated out and external charges are placed at sites k, l
 */
  SU2HamiltonianGIS(int N,double epsilon,double mu,double g,int k,int l,double lambda,double eta);
  
  ~SU2HamiltonianGIS();
  
  /** Return the MPO representation of the Hamiltonian*/
  const MPO& getHMPO() const {return hamil;}
  /** Return whether the MPO representation exists*/
  bool hasHMPO() const {return true;}
  
  /** Getter for the value of the interaction \f$ \varepsilon \f$*/
  double get_epsilon() const {return this->epsilon;};
  
  /** Setter for the value of the interaction \f$ \varepsilon \f$*/
  void set_t(double epsilon_value) {this->epsilon=epsilon_value;}; 
  
  /** Getter for the value of the mass \f$ \mu \f$*/
  double get_mu() const {return this->mu;};
  
  /** Setter for the value of the interaction \f$ \mu \f$*/
  void set_mu(double mu_value) {this->mu=mu_value;}; 
  
  /** Getter for the value of the penalty strength \f$ \lambda \f$*/
  double get_lambda() const {return this->mu;};
  
  /** Setter for the value of the penalty strength \f$ \lambda \f$*/
  void set_lambda(double la_value) {this->lambda=la_value;}; 
   
  /** As the Hamiltonian after changing some parameters with setter functions, it might not be immediately needed, the user has to update the Hamiltonian MPO manually with this function */
  void updateHMPO(){this->initHamiltonian();};
  
  /** Get the MPO for \f$ Q^{\alpha}_n Q^{\alpha}_n \f$, the \f$ \alpha \f$ is not summed since it is the same for every component and would only add a factor of three (if the summed version is desired, one has to add the factor by hand). The site has to be between 0 and N_sites-1 */
  void getQsquareMPO(MPO &mpo, int site) const;
  
  /** Get the MPO for \f$ Q^{\alpha}_n\f$  The site has to be between 0 and N_sites-1 */
  void getQalphaMPO(MPO &mpo, int site, Component alpha) const;
  
  /** Get the MPO for the external charge correlation function \f$ q^{\alpha}_kq^{\alpha}_l \f$ */
  void getExtChargeCorrelationMPO(MPO &mpo,Component alpha) const; 
  
};

#endif //SU2GIS_H
