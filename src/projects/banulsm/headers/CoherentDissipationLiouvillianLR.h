/**
   \file CoherentDissipationLiouvillianLR.h 
   Implementation of the Liouvillian for a spin chain governed by the
   following model.
    \f[
    \frac{d\rho}{d t}=i[\rho,H]+\gamma \left( L
    \rho L^{\dagger}-\frac{1}{2} \rho L^{\dagger} L
    -\rho L^{\dagger} L \rho\right),
    \f]
    where the system Hamiltonian is
   \f[
   H=g\sum_i \sigma_x^{[i]},
   \f]   
   and the dissipation is given by
   \f[
   L=\sum_i {\sigma^{-}}^{[i]}
   \f]
   The model is constructed with parameters \f$N\f$, \f$g\f$ and
   \f$\gamma\f$.

   Notice that if $\g$ and $\gamma$ are multiplied by a common factor,
   the effect is that of rescaling the full superoperator by the same
   global factor.

   \author Mari-Carmen Banuls
   \date 25/03/2013
*/

#ifndef COHDISLIOUVLR_H
#define COHDISLIOUVLR_H

#include "mwArray.h"
#include "CoherentDissipationLiouvillian.h"

class CoherentDissipationLiouvillianLR:
public CoherentDissipationLiouvillian{
  /* int d; // physical dimension (fixed to two in this model) */
  /* int N; // system size */
  /* double g; // Coefficient of system Hamiltonian */
  /* double gamma; // Coefficient in front of dissipative terms */
  /* MPO Lmpo; // The MPO for the Liouvillian */
  /* MPO LAdjmpo; // The MPO for the adjoint of the Liouvillian */

  /* MPO trLmpo; // Special MPO for L on the purification, after trace */

  /* bool isInit; // A flag to know if the mpo was already prepared */

 public: 
  /** Constructor */
  CoherentDissipationLiouvillianLR(int N,double g,double gamma);
  ~CoherentDissipationLiouvillianLR();
  /** Recover the MPO for the Liouvillian superoperator */
  //  const MPO& getLMPO(){return Lmpo;}
  // const MPO& getLAdjointMPO(){return LAdjmpo;}

  /** Return the MPO for trace(L rho) tensor the identity on the
      dimension of the purification */
  /* void getTrLMPO(int dpur,MPO& result); */

  /* /\** Return the properties of the model *\/ */
  /* int getLength()const{ return N;} */
  /* double getG() const{return g;} */
  /* double getGamma() const{return gamma;} */

 protected:
  /** This is the one doing the work */
  void initL();
  /* /\** Auxiliary function to construct the operators *\/ */
  /* void initZ(mwArray& Z); */

  /* /\** This just takes the result from the previous method and changes */
  /*     some tensors. *\/ */
  /* void initLAdjoint(); */

  /* /\** Idem, for the trace after the Liouvillian. This is only to */
  /*     check, as by construction, it should be zero.  *\/ */
  /* void initTrL(); */

  /* void clear(); */

  /* /\** Auxiliary function that converts two individual single-site */
  /*     operators in their tensor product and sets the indices in the */
  /*     correct order for the MPO, as repeatedly needed by initZ(). *\/ */
  /* void constructOperatorProduct(mwArray& result,const mwArray& opA, */
  /* 				const mwArray& opB); */
};

#endif // COHDISLIOUVLR_H
