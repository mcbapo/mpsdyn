#ifndef SCHWINGERZR_H
#define SCHWINGERZR_H

#include "MPS.h"
#include "SchwingerHamiltonian.h"

/** 
    This is the Schwinger Hamiltonian in the spin representation, with
    an additional parameter, z, which penalizes magnetization
    \f$<S_z>\neq 0\f$, and an extra term, with parameter
    \f$\lambda\f$, that penalizes states with a negative sign for the
    expectation value of the \f$S_R\f$ operator (approx, translation
    by one times \f$S_x\f$).  The later does not commute with the
    original Hamiltonian, except in the thermodynmic limit so that we
    expect that it introduces important finite size effects in the
    spectrum.  */

class SchwingerHamiltonianSzSR: 
public SchwingerHamiltonian{
	double zpen;
	double lambda;
public:
  SchwingerHamiltonianSzSR(int N,double mu,double x,double l0,double alpha,
			   double z,double lambda,
			   double offset=0,double nu=0.,double y=1.);
  ~SchwingerHamiltonianSzSR(){};
protected:
  void initHMPO();

};

#endif // SCHWINGERZR_H
