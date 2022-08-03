#ifndef Z2_H
#define Z2_H

#include "MPS.h"
#include "Hamiltonian.h"
#include "Contractor.h"

/**
   \file Z2Hamiltonian.h Definition of the class that implements the Z2 Hamiltonian in spin formulation
   
   \author Stefan Kühn
   \date 14/12/2017
*/


/** List of possible single body operators*/
enum SingleBodyOperator{
  sxOp,	//sigma_x
  syOp,	//sigma_y
  szOp,	//sigma_z
  spOp,	//sigma_plus
  smOp 	//sigma_minus
};

/** Enum for operator components*/
enum Component{
  x_comp, y_comp, z_comp
};


class Z2Hamiltonian: 
public Hamiltonian{
private:
  //Pauli Matrices
  mwArray sigma_x,sigma_y,sigma_z,sigma_plus,sigma_minus;
  //Idenitity
  mwArray id;
  //Properly reshaped identity which can be used for single body MPOs
  mwArray id_mpo;    
  //Array to hold all the relevant basic spin operators
  mwArray Z;
  
  /**Get the spin on site \f$ n \f$, where n is now an index ranging from \f$ 0 \f$ to \f$ 2L-1 \f$. */
  void getSingleBodyMPO(MPO &Spin, int index,SingleBodyOperator Op) const; 

  
  protected:
  //Number of sites
  int L;
  //Parameters of the Hamiltonian
  double x,mu,lambda,alpha;
  //Boundary spins
  double c_left,c_right;
  //Hamiltonian MPO
  MPO hamiltonian;
  
   /** Generate some matrices used for single and few body MPOs to be saved in the object, therefore in the MPO the pointers to the matrices can be stored and do not consume memory*/
  void initOperators(void);
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously)*/
  void initHamiltonian(void);
  /** Construct the actual Hamiltonian MPO (which uses the operators that have to be initialized previously)*/
  void initHamiltonianNoPenalty(void);
  
  public:
  /** Creates an instance of the Z_2-Hamiltonian in spin formulation which is given by 
   \f[ W= -x\sum_{l=0}^{L-2}\left(\sigma^{+}_l\sigma^{x}_{l,l+1}\sigma^{-}_{l+1} +\mathrm{h.c}\right)+\frac{\mu}{2}\sum_{l=0}^{L-1}(-1)^l  (1+\sigma^{z}_{l}) + \sum_{l=1}^{L-2} (2-2\sigma^{z}_{l,l+1}) \f]   
 */
  Z2Hamiltonian(int L,double x,double mu);
  
  /** Creates an instance of the Z_2-Hamiltonian in spin formulation
   \f[ W= -x\sum_{l=0}^{L-2}\left(\sigma^{+}_l\sigma^{x}_{l,l+1}\sigma^{-}_{l+1} +\mathrm{h.c}\right)+\frac{\mu}{2}\sum_{l=0}^{L-1}(-1)^l  (1+\sigma^{z}_{l}) + \sum_{l=1}^{L-2} (2-2\sigma^{z}_{l,l+1}) + \lambda P\f]
   including a penalty term 
   \f[ P = 2\sum_{l=0}^{L-1} \left(1+ (-1)^{l}\sigma^{z}_{l-1,l}\sigma^{z}_{l}\sigma^{z}_{l,l+1}\right) \f]
   to enforce the gauge invariant subspace */   
  Z2Hamiltonian(int L, double x, double mu, double lambda);
  
  /** Creates an instance of the Z_2-Hamiltonian in spin formulation whith background field alpha
   \f[ W= -x\sum_{l=0}^{L-2}\left(\sigma^{+}_l\sigma^{x}_{l,l+1}\sigma^{-}_{l+1} +\mathrm{h.c}\right)+\frac{\mu}{2}\sum_{l=0}^{L-1}(-1)^l  (1+\sigma^{z}_{l}) + \sum_{l=1}^{L-2} (2-L_{l,l+1}) + \lambda P\f]
   including a penalty term 
   \f[ P = 2\sum_{l=0}^{L-1} \left(1+ (-1)^{l}\sigma^{z}_{l-1,l}\sigma^{z}_{l}\sigma^{z}_{l,l+1}\right) \f]
   to enforce the gauge invariant subspace, the operator \f$ L=P^\dagger + P = \left(\begin{array}{cc} 2\cos(\pi\alpha) & 0 \\ 0 & -2\cos(\pi\alpha) \end{array}\right) \f$ reduces to the standard \f$ 2\sigma^z \f$ for \f$ \alpha=0 \f$.*/   
  Z2Hamiltonian(int L, double x, double mu, double alpha, double lambda);
  
  ~Z2Hamiltonian();
  
  /** Return the MPO representation of the Hamiltonian */
  const MPO& getHMPO() const {return hamiltonian;}
  
  /** Return whether the MPO representation exists */
  bool hasHMPO() const {return true;}
  
  /** Update the Hamiltonian after changing some parameters */
  void updateHMPO(){this->initHamiltonian();};
  
  /** Getter for the value of the interaction \f$ x \f$ */
  double get_x() const {return this->x;};
  
  /** Setter for the value of the interaction \f$ x \f$*/
  void set_x(double x_value) {this->x=x_value;}; 
  
  /** Getter for the value of the mass \f$ \mu \f$ */
  double get_mu() const {return this->mu;};
  
  /** Setter for the value of the interaction \f$ \mu \f$ */
  void set_mu(double mu_value) {this->mu=mu_value;};   
  
  /** Get the MPO for penalty enfocing the gauge invariant sector \f$ P=2\sum_l(1-G_l)^\dagger(1-G_l) \f$ */
  void getPenaltyMPO(MPO &mpo) const;
  
  /** Get the Number of sites */
  int getSites(){return this->L;};
};

#endif //Z2_H
