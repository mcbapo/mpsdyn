#ifndef SCHWINGER_H
#define SCHWINGER_H

#include "MPS.h"
#include "Hamiltonian.h"

/** 
    Implementation of the Schwinger Hamiltonian in the spin
    representation

*/

class SchwingerHamiltonian: 
public Hamiltonian{
protected:
  int N;
  double mu,x,l0,alpha,nu;
  MPO hamil;
  mwArray Z; 
  double offset;
  double gweight;

 public:
  /** Creates an instance of Schwinger Hamiltonian for a chain of
      length N, and parameters mu, x, l0 and alpha.  If offset is
      given, it is used to ofset the ground state energy, by
      constructing, instead H+offset*N If nu is explicit, it includes
      a chemical potential term (nu=mu/g) Finally, the parameter y can
      be set to change the relative weight of the gauge terms, and
      eventually check the non-interacting case, by setting it to
      zero. If not set, or set to 1, the normal Schwinger Hamiltonian
      is implemented.  WARNING: The y parameter is not included in the
      exponential of HL yet!!
  */
  SchwingerHamiltonian(int N,double mu,double x,double l0,double alpha,
		       double offset=0,double nu=0.,double y=1.);

  ~SchwingerHamiltonian();

  /** Basic functions to return the values of the parameters */
  int getN() const{return N;}
  double getMu() const{return mu;}
  double getX() const{return x;}
  double getNu() const{return nu;}
  double getOffset() const{return offset;}
  double getGweight() const{return gweight;}

  /** Return a reference to the MPO for the Hamiltonian. */
  const MPO& getHMPO() const{return hamil;}

  bool hasHMPO() const {return true;}

  /** Operators, as MPOs */
  /** Construct the MPO for the operator:
      \f[
      \frac{\sqrt{x}}{L}\sum_n (-1)^n \[ \sigma_n^{+} \sigma_{n+1}^{-} +\sigma_{n+1}^{+} \sigma_n^{-}\]
      \f]
  */
  void constructGamma5MPO(MPO& gamma5);

  /** Construct the MPO for the operator:
      \f[
      \ell_0-\frac{1}{4}+\frac{1}{2L}\sum_{n=1}^{L-1} (-1)^n (L-n) \sigma_z^{[n]}
      \f]
  */
  void constructGammaAlphaMPO(MPO& gammaA);

  /** The operator is:
      \f[
      \frac{\sqrt{x}}{L}\sum_n (-1)^n \frac{1+\sigma_z^{[n]}}{2}
      \f]
  */
  void constructCondensateMPO(MPO& gammaC);

  /** Construct the MPO that corresponds to the momentum in the
      continuum (hopefully) \f[-i x \sum_n [\sigma_n^- \sigma_{n+1}^z
      \sigma_{n+2}^+ - h.c.].\f] This corresponds to the
      discretization of the continuum operator \f[P=\int dx
      \Psi^{\dagger}(x)(i \partial_x+g A_1(x))\Psi(x))\f], rescaled by
      \f$2/(a g^2)\f$, as the Hamiltonian.  The MPO (with bond
      dimension 6) is returned in the argument, except for the factor
      \f$x\f$, which is not included.  */
  void constructMomentumMPO(MPO& Pmpo) const;

  /** Construct the MPO that corresponds to the momentum squared in the continuum (hopefully)
      \f[x \sum_n [(1+\sigma_n^z) +\sigma_{n}^-\sigma_{n+1}^z \sigma_{n+2}^+ +
      \sigma_{n}^+\sigma_{n+1}^z \sigma_{n+2}^- ].\f]
      In the continuum limit this is \f$\frac{(i\partial)^2}{g^2}\f$. 
      The MPO (with bond dimension 6) is returned in the argument.  */
  void constructMomentumSquaredMPO(MPO& Pmpo) const;

  /** 
      Construct an MPO (bond dimension \f$d^2\f$) that represents a
      cyclic translation in the chain (i.e., everyone is moved one to
      the right, and the last site is translated to the first one).
      The result is returned in T.
      If \param<right> is false, the adjoint translation is constructed
   */
  void constructCyclicTranslationMPO(MPO& T,bool right=true) const;

  /** This one, on top, rotates sigmax the odd sites */
  void constructCyclicTranslationRotateMPO(MPO& T,bool right=true) const;

  /** Construct a translation operator which exchanges also up and
      down, so it would be the translation in the basis where the x=0
      vacuum is all 0. For the boundaries, the first site is left in
      0, and the rightmost one is not shifted, if it was a 1, and it
      is overriden by the one to the left if it was a 0, so that |001>
      -> |001>, |010>->|001>, |011>->|001> 
      If right=false, the translation is done to the left. */
  void constructShiftRotateMPO(MPO& T,bool right=true) const;

/** Construct the MPO for the kinetic (XY) term, and return it in the
    argument.  */
  void constructVMPO(MPO& Vmpo) const;

/** Construct (and return in the argument) an MPO for the operator
\f[\sum (1+\sigma_z)_k (1+\sigma_z)_{k+1}\f], where the sum runs over even
sites if even=true (default) and over odd ones otherwise. */
  void constructZZMPO(MPO& Vmpo,bool even=true) const;

  /**
Construct (and return in the argument) an MPO for the operator
\f[\sum \sigma^+_k (1+\sigma^z)_{k+1} \sigma^-_{k+2} \f]
 */
  void constructShiftPairMPO(MPO& Vmpo) const;

/** Construct the vector current operator i(XY-YX) as a MPO and return
 * it in the argument. */
  void constructVecMPO(MPO& Vmpo) const;
  /** Idem for an operator corresponding to "quasimomentum" k0
      (between 1 and N) in OBC, i.e., ansatz \f[\sum sin\frac{k_0 \pi
      m}{N+1} \hat{O}_m \f]*/
  void constructVecMPOobc(MPO& Vmpo,int k0=1) const;
/** Construct the vector current operator i(XY-YX)_(2k,2k+1) as a MPO
 * and return it in the argument. Notice that in this version, the
 * two-body operator acts only on half the links (the even ones). */
  void constructVecMPOe(MPO& Vmpo) const;
/** Idem, for just the two body term */
  void constructVecTwoSiteMPO(MPO& Vmpo) const;

/** Construct the scalar operator as a MPO and return it in the
 * argument. */
  void constructScalMPO(MPO& Smpo) const;
  void constructScalMPOobc(MPO& Smpo,int k0=1) const;
/** Idem, for just the two body term */
  void constructScalTwoSiteMPO(MPO& Vmpo) const;

/** Return the parity operator, which in this case changes sign of
 * only lower components (as gamma_0) 
 TODO: Exchange position??*/
  void constructParityMPO(MPO& Pmpo) const;



/** Apply parity transformation to the given MPS */
  void applyParity(const MPS& orig,MPS& result) const;

/** Apply parity transformation to the given MPS (site n<->N-n-1), and also 
    a sigmax rotation  */
  void applyReflectionRotation(const MPS& orig,MPS& result) const;

  complex_t computeParity(const MPS& orig) const;

  // Implement now the unitaries as Dyson series
  bool hasUMPO() const{return true;}
  void getBasicUOperators(int pos,double initT,double delta,int orderT,
			  vector<Operator&>& opSet,
			  bool imag=false,bool transverse=true);

  int getUOpersPerStep(int orderT) const;
  shrt::Indices getUOperDimensions(int orderT,int opernr=-1,
				   bool transverse=true) const;

  /** Construct the MPO for the unitaries in the Trotter expansion
      with the splitting even-odd in the Hamiltonian. Take as argument
      the complex factor in the exponential, delta.  For real time
      evolution, delta is -i*tau, for imaginary time, it should be
      -tau.  The Hamiltonian is splitted in the \f$H_x\f$ and
      \f$H_L\f$ terms, the first one containing only the XY (two-body)
      term, proportional to the parameter x.  The exponential of the
      second is approximated by \f$1+\delta H_L\f$ (with \f$\delta\f$
      appropriately real or complex) This approximates the exponential
      of the \f$x\f$ term as a product of two exponentials, only, for
      the even adn odd terms, which are grouped together in a single
      MPO with bond dimensions D=4 everywhere. It would be a better
      approximation to use a Trotter-style decomposition for this
      term.
  */

  /** Construct the MPO corresponding to the evolution operator with
      the splitting even-odd in the Hamiltonian. Take as argument the
      complex factor in the exponential, delta.  For real time
      evolution, delta is -i*tau, for imaginary time, it should be
      -tau. */
  void getExponentialMPOx(MPO& expHx,complex_t delta) const;
  /** This version implements the exponential of either the even or
      the odd part of the two-body term */
  void getExponentialMPOx(MPO& expHx,complex_t delta,bool even) const;
  /** Returns an MPO, with bond dimension D=3, which approximates 
      the exponential of delta H_L by Id+delta H_L.
      It also contains the mass term. */
  void getExponentialMPOL(MPO& expHL,complex_t delta) const;

  /** To try higher order Taylor, I need HL alone */
  void getMPOHL(MPO& HL) const;

  /** Returns an MPO which computes the "exact" exponential of the H_l
      term by recovering the value of the electric flux on every link
      from the spin content before. The bond dimension will be up to
      N-1, because we are not imposing that total charge is zero. In
      that case (TODO) as SchwingerHamiltonianSz, it could be reduced
      to N/2. 
      If cutoff!=0 is given, the bond dimension is truncated to the 
      specified value using a SVD of the exact operator. 
      If lmax!=0 is also given, the maximum electric field allowed per 
      link is truncated to lmax before the SVD, and 2*lmax+1 is used as cutoff.
      This is equivalent to truncating the maximum number of bosons
      If normFactor is provided, each individual operator is normalized by 
      this number, so that the overall MPO is normalized by normFact^L
*/
  void getExactExponentialMPOL(MPO& expHL,complex_t delta,int cutoff=0,double tol=0.,int lmax=0,double normFactor=1.) const;

/** Get a series of exponentials for the long range terms, specifying
 * the last site involved (between 1 and N-2). Notice that, if this
 * option is used, the single body terms need to be implemented and
 * applied on its own (although they do not have an extra cost,
 * because they have D=1). */
//  void getPartExponentialMPOL(MPO& expHL,complex_t delta,int second) const;
// I also need the single body terms!!

 protected:
/** Empty constructor only for the daughter */
  SchwingerHamiltonian();
  void initZ();
  virtual void initHMPO();

  /** Prepare the even-odd decomposition of the evolution operator,
      for the \f$H_x\f$ term, which can be expressed in this way as an
      MPO. This method returns the two terms to be inserted in the MPO
      for the exp(delta H_x) (delta a complex number).  TI is assumed 
      here.
  */
  void getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta) const;

  /** Compute the two-body term of the Hamiltonian: \f[ h_{i,i+1}= 
      \frac{x}{2}
      (\sigma_{x}^{[i]} \sigma_{x}^{[i+1]} + \sigma_{y}^{[i]}
      \sigma_{y}^{[i+1]} )\f] */
  void computeTwoBodyTerm(mwArray& result) const;

  /** Auxiliaries for the construction of operators U. If
      transverse=true, the order of indices is that of a transverse
      MPO.  Returns just an mwArray, with the proper shape, to use in
      the creation of a new Operator as needed.
   */
  void constructOperatorUxy(mwArray& result,double delta,int pos,
			    bool imag=false,bool transverse=false) const;
  void constructDforUxy(mwArray& D,double delta,int pos,bool imag) const;

  /** Auxiliary function that constructs the MPO for one term in the
      parity mess. */
  void constructParityMPO(MPO& Pmpo,int k) const;

  void constructSignParityMPO(MPO& Pmpo) const;
};


#endif // SCHWINGER_H
