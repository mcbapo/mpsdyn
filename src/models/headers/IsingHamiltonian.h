/**
   \file IsingHamiltonian.h
   Hamiltonian for the Ising model (with parallel and transverse magnetic field)
   
   \author Mari-Carmen Banuls
   \date 20/10/2011

*/


#ifndef ISINGHAMIL_H
#define ISINGHAMIL_H

#include "Hamiltonian.h"
#include "MPO.h"

typedef double (*tdepfun)(double);

/** 
    Implement the (time dependent) Ising Hamiltonian
    \f[
    H=  \sum_i[J(t) sigmaz_i sigmaz_(i+1) + g(t) sigmax_i+ h(t) sigmaz_i] + K
    \f]
    with a possible offset K.
*/
class IsingHamiltonian: 
public Hamiltonian {
 protected:
  int N; // length (if infinite, <0?)
  int d; // physical dimension
  tdepfun J,g,h;
  bool timedep; // to save calculations if no time dependency ??
  double J0,g0,h0; // initial values of the constants
  double K; // offset
  mwArray Z; 
  double t;
  MPO hamil; // MPO for the Hamiltonian

 public:
  IsingHamiltonian(int N,int d,double J0,double g0,double h0=0.,double K=0.);
  virtual ~IsingHamiltonian();
  /* Register function pointers that calculate the time dependency of the
     coefficients J and g. */
  void registerTimeDependence(tdepfun J=0,tdepfun g=0,tdepfun h=0);
  bool hasHMPO() const {return true;}

  /* First change the value of t, for a tme-dependent case, and then
     return the MPO for the corresponding Hamiltonian. */
  const MPO& getHMPO(double t);
  const MPO& getHMPO() const {return hamil;}  //getHMPO(0.);}
  // Implement now the unitaries 
  bool hasUMPO() const{return true;}

  /** Get the MPO for the exponential of the Hamiltonian times a certain
   * coefficient, which can be real (-delta, negative, for imaginary
   * time evolution) if imag=true, or imaginary (-i delta), if
   * imag=false. It returns the MPO corresponding to second order Trotter
   * with the decomposition H2, H1. If orderT!=2, not supported. 
   * \param t is the time at which the Hamiltonian is computed
   * (relevant only if time dependent).
   */
  void getUMPO(MPO& mpo,double delta,double t,bool imag=false,
	       int orderT=2,int nrOp=0);
  /** Same, with complex argument delta, to allow recursive
      construction of higher order Trotter decompositions. */
  void getUMPO(MPO& mpo,complex_t deltaC,double t,
	       int orderT=2,int nrOp=0);

  /** Get only the mwArray corresponding to one of the operators in the 
      MPO, as required to construct a transverse column, for instance.      
   */
  mwArray getUOperator(double delta,double t,bool imag=false,
		       int orderT=2,int nrOp=0);

  /** Construct an MPo for the commutator with an operator, and return in in the argument.  */
  void getCommutatorMPO(MPO& mpo,double t=0.);

 private:
  void initZ();
  void setHMPO(MPO& mpo,double J,double g, double h);

  /** Auxiliary functions to construct the operators of the exponential 
      MPOs */
  /* ******* Inherited from Matlab ********** */

  /** Computes the MPO tensor (TI) for the evolution operator arising from
   * the two body term in the Ising Hamiltonian, 
   * \f[exp[\epsilon \sum \sigma_z\otimes\sigma_z]\f] 
   * Order Dr x du x Dl x dd
   * if pos==0 or N-1, it returns the corresponding ege term, otherwise, 
   * it returns the one for the middle sites (also the infinite)
   */
  mwArray isingU2MPO(complex_t epsilon,int pos=1) const;

  /** Computes the MPO tensor (TI) for the evolution operator arising from
   * the one body term in the Ising Hamiltonian, 
   * \f[exp[\epsilon\sum(g \sigma_x+h\sigma_z)]\f]
   */
  mwArray isingU1MPO(complex_t epsilon,double g,double h) const;

  /** Return the MPO (element at one site) for Ising Hamiltonian evolution,
   * second order of the Trotter expansion, using width \param<delta>, i.e.
   * \f[e^{\frac{\delta}{2} H_1(t) }e^{\delta H_2(t) }e^{\frac{\delta}{2} H_1(t) }\f] 
   * For real time evolution, \param<delta>=\f$-i\delta\f$ should be used, 
   * while for imaginary time evolution, \param<delta>=\f$-\delta\f$, for 
   * some real \f$\delta\f$.
   * The returned mwArray has dimensions du x Dl x dd x Dr
   */
  mwArray constructU2T(double t,complex_t delta,int pos=1);

  /** 
      Auxiliary function for the construction of the commutator
      MPO. If this should be more general (not only for this
      IsingHamiltonian case), these functions might need to be
      somewhere else. Specially the double
   */
  void initZdoub(mwArray& Z);
  void constructOperatorProduct(mwArray& result,const mwArray& opA,
				const mwArray& opB);


};


#endif // ISINGHAMIL_H
