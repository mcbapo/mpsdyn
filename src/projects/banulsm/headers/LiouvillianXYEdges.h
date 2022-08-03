/**
   \file LiouvillianXYEdges.h 
   Implementation of the Liouvillian for a spin chain of \f$N\f$ sites,
   governed by the following model.
    \f[
    \frac{d\rho}{d t}=i[\rho,H]+\sum_{\alpha}\left( L_{\alpha}
    \rho L_{\alpha}^{\dagger}-\frac{1}{2} \rho L_{\alpha}^{\dagger} L_{\alpha}
    -\rho L_{\alpha}^{\dagger} L_{\alpha} \rho\right),
    \f]
    where the system Hamiltonian is the XY model
   \f[
   H=\frac{1+\gamma}{2}\sum {\sigma_{x}}^{[i]}{\sigma_{x}}^{[i+1]}
   +\frac{1-\gamma}{2}\sum {\sigma_{y}}^{[i]}{\sigma_{y}}^{[i+1]}
   + h\sum_i \sigma_z^{[i]},
   \f]   
   and the dissipation is given by
   \f[
   L_{1,2}=\sqrt{\Gamma_{1,2}^{L}}{\sigma^{\mp}}^{[1]}
   \f]
   \f[
   L_{3,4}=\sqrt{\Gamma_{3,4}^{L}}{\sigma^{\mp}}^{[N]}.
   \f]
   The model is constructed with parameters \f$N\f$, \f$h\f$,
   \f$\gamma\f$, \f$\Gamma_{1,2}^L\f$ and \f$\Gamma_{3,4}^R\f$ .

   \author Mari-Carmen Banuls
   \date 25/09/2013
*/

#ifndef LIOUVXY_H
#define LIOUVXY_H

#include "mwArray.h"
#include "Liouvillian.h"

class LiouvillianXYEdges:
public Liouvillian{
 protected:
  int d; // physical dimension (fixed to two in this model)
  int N; // system size
  double h; // Coefficient of system Hamiltonian
  double gamma; // Coefficient of system Hamiltonian
  double GammaL1,GammaL2; // Coefficient in front of dissipative terms on
                          // first site 
  double GammaR3,GammaR4; // Coefficient in front of dissipative terms on
                          // last site 

  double scale; // scale factor that multiplies the whole Liouvillian
  MPO Lmpo; // The MPO for the Liouvillian
  MPO LAdjmpo; // The MPO for the adjoint of the Liouvillian

  MPO trLmpo; // Special MPO for L on the purification, after trace

  bool isInit; // A flag to know if the mpo was already prepared

 public: 
  /** Constructor */
  LiouvillianXYEdges(int N,double gamma,double h,double GammaL1,
		     double GammaL2,double GammaR3,double GammaR4,double scale=1.);
  ~LiouvillianXYEdges();
  /** Recover the MPO for the Liouvillian superoperator */
  const MPO& getLMPO(){return Lmpo;}
  const MPO& getLAdjointMPO(){return LAdjmpo;}

  /** Return the MPO for trace(L rho) tensor the identity on the
      dimension of the purification */
  void getTrLMPO(int dpur,MPO& result);

  /** Return the properties of the model */
  int getLength()const{ return N;}
  double geth() const{return h;}
  double getgamma() const{return gamma;}
  double getGammaL1() const{return GammaL1;}
  double getGammaL2() const{return GammaL2;}
  double getGammaR3() const{return GammaR3;}
  double getGammaR4() const{return GammaR4;}

 protected:
  /** Empty constructor for the daughter */
//  LiouvillianXYEdges(int N,double g,double gamma,bool priv);


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

  void clear();

  /** Auxiliary function that converts two individual single-site
      operators in their tensor product and sets the indices in the
      correct order for the MPO, as repeatedly needed by initZ().
      TODO: Should be implemented in base class Liouvillian!
  */
  void constructOperatorProduct(mwArray& result,const mwArray& opA,
				const mwArray& opB);
};

#endif // LIOUVXY_H
