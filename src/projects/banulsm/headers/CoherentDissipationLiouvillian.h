/**
   \file CoherentDissipationLiouvillian.h 
   Implementation of the Liouvillian for a spin chain governed by the
   following model.
    \f[
    \frac{d\rho}{d t}=i[\rho,H]+\sum_{\alpha}\gamma_{\alpha}\left( L_{\alpha}
    \rho L_{\alpha}^{\dagger}-\frac{1}{2} \rho L_{\alpha}^{\dagger} L_{\alpha}
    -\rho L_{\alpha}^{\dagger} L_{\alpha} \rho\right),
    \f]
    where the system Hamiltonian is
   \f[
   H=g\sum_i \sigma_x^{[i]},
   \f]   
   and the dissipation is given by
   \f[
   L_{\alpha}={\sigma^{-}}^{[i]}+\epsilon {\sigma^{-}}^{[i+1]}
   \f]
   for \f$\alpha=1,\ldots N-1\f$.
   The parameter \f$\epsilon\f$ can be decided upon construction (default +1).
   The model is constructed with parameters \f$N\f$, \f$g\f$ and
   \f$\gamma\f$.

   Notice that if $\g$ and $\gamma$ are multiplied by a common factor,
   the effect is that of rescaling the full superoperator by the same
   global factor.

   \author Mari-Carmen Banuls
   \date 7/02/2013
*/

#ifndef COHDISLIOUV_H
#define COHDISLIOUV_H

#include "mwArray.h"
#include "Liouvillian.h"

class CoherentDissipationLiouvillian:
public Liouvillian{
 protected:
  int d; // physical dimension (fixed to two in this model)
  int N; // system size
  double g; // Coefficient of system Hamiltonian
  double gamma; // Coefficient in front of dissipative terms
  complex_t epsilon; // relative phase of the sigma- terms in the dissipation
  MPO Lmpo; // The MPO for the Liouvillian
  MPO LAdjmpo; // The MPO for the adjoint of the Liouvillian

  MPO trLmpo; // Special MPO for L on the purification, after trace

  MPO darkSubspace; // projector onto the subspace of eigenvalue 0 of the purely dissipative model.

  bool isInit; // A flag to know if the mpo was already prepared

 public: 
  /** Constructor */
  CoherentDissipationLiouvillian(int N,double g,double gamma,complex_t epsilon=ONE_c);
  ~CoherentDissipationLiouvillian();
  /** Recover the MPO for the Liouvillian superoperator */
  const MPO& getLMPO(){return Lmpo;}
  const MPO& getLAdjointMPO(){return LAdjmpo;}
  const MPO& getDarkSubspaceProjectorMPO(){return darkSubspace;}


  /** Return the MPO for trace(L rho) tensor the identity on the
      dimension of the purification */
  void getTrLMPO(int dpur,MPO& result);

  /** Return the properties of the model */
  int getLength()const{ return N;}
  double getG() const{return g;}
  double getGamma() const{return gamma;}

 protected:
  /** Empty constructor for the daughter */
  CoherentDissipationLiouvillian(int N,double g,double gamma,bool priv);


  /** This is the one doing the work */
  virtual void initL();
  /** Auxiliary function to construct the operators */
  void initZ(mwArray& Z);

  /** This just takes the result from the previous method and changes
      some tensors. */
  void initLAdjoint();

  /** Idem, for the trace after the Liouvillian. This is only to
      check, as by construction, it should be zero.  */
  void initTrL();

  /** Prepare the projector onto the dark subspace of L alone (no H)*/
  void initProjectorDS();

  void clear();

  /** Auxiliary function that converts two individual single-site
      operators in their tensor product and sets the indices in the
      correct order for the MPO, as repeatedly needed by initZ(). */
  void constructOperatorProduct(mwArray& result,const mwArray& opA,
				const mwArray& opB);
};

#endif // COHDISLIOUV_H
