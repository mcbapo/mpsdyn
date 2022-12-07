/**
   \file IsingLocalDissipLiouvillian.h 
   Implementation of the Liouvillian for a spin chain of \f$N\f$ sites,
   governed by the following model.
   \f[
    \frac{d\rho}{d t}=i[\rho,H]+\sum_{\alpha}\left( L_{\alpha}
    \rho L_{\alpha}^{\dagger}-\frac{1}{2} \rho L_{\alpha}^{\dagger} L_{\alpha}
    -\rho L_{\alpha}^{\dagger} L_{\alpha} \rho\right),
    \f]
    where the system Hamiltonian is the Ising chain
   \f[
   H=\frac{V}{4}\sum_{i=0}^{N-2} {\sigma_{z}}^{[i]}{\sigma_{z}}^{[i+1]}
   +\frac{\Omega}{2}\sum_{i=0}^{N-1} {\sigma_{x}}^{[i]}-\frac{V-\Delta}{2}\sum_{i=0}^{N-1} {\sigma_{z}}^{[i]}
   +\frac{V}{4}{\sigma_{z}}^{[0]}+\frac{V}{4}{\sigma_{z}}^{[N-1]}
   \f]   
   and the dissipation is local, given by the term
   \f[
   L_{i}=\sqrt{\gamma}{\sigma^{+}}^{[i]}
   \f]
   The model is constructed with parameters \f$N\f$, \f$V\f$, \f$\Omega\f$, \f$\Delta\f$,
   and \f$\gamma\f$.

   \author Mari-Carmen Banuls
   \date 6/12/2022
*/

#ifndef ISINGOPENLIOUV_H
#define ISINGOPENLIOUV_H

#include "mwArray.h"
#include "Liouvillian.h"

class IsingLocalDissipLiouvillian:
public Liouvillian{
 protected:
  int d; // physical dimension (fixed to two in this model)
  int N; // system size
  double V,Omega,Delta; // Coefficients of system Hamiltonian
  double gamma; // Coefficient of dissipative terms

  double scale; // scale factor that multiplies the whole Liouvillian

  mwArray Z; // auxiliary place for operators
  
  MPO Lmpo; // The MPO for the Liouvillian
  MPO LAdjmpo; // The MPO for the adjoint of the Liouvillian

  MPO trLmpo; // Special MPO for L on the purification, after trace

  //bool isInit; // A flag to know if the mpo was already prepared

 public: 
  /** Constructor. Besides the parameters in the description, it includes a (last) optional parameter \param scale
      to rescale the whole operator.  */
  IsingLocalDissipLiouvillian(int N,double V,double Omega,double Delta,
			      double gamma,double scale=1.);
  ~IsingLocalDissipLiouvillian();
  /** Recover the MPO for the Liouvillian superoperator */
  const MPO& getLMPO(){return Lmpo;}
  const MPO& getLAdjointMPO(){return LAdjmpo;}

  /* /\** Return the MPO for trace(L rho) tensor the identity on the */
  /*     dimension of the purification *\/ */
  /* void getTrLMPO(int dpur,MPO& result); */

  /** Return the properties of the model */
  //  int getLength()const{ return N;}

 protected:

  /** This is the one doing the work, completely analogous to initH in Hamiltonians */
  virtual void initL();
  /** Auxiliary function to construct the operators */
  void initZ();

  /** This just takes the result from the previous method and changes
      some tensors. */
  void initLAdjoint();

  /** Idem, for the trace after the Liouvillian. This is only to
      check, as by construction, it should be zero.  */
  //  void initTrL();

  void clear();
};

#endif // ISINGOPENLIOUV_H
