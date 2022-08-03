#ifndef SCHWINGERZ_H
#define SCHWINGERZ_H

#include "MPS.h"
#include "SchwingerHamiltonian.h"

/** 
    This is the Schwinger Hamiltonian in the spin
    representation, with an additional parameter, z, which penalizes
    magnetization \f[<S_z>\neq Stot\f].

*/

class SchwingerHamiltonianSz: 
public SchwingerHamiltonian{
  double zpen;double scale;double Stot;
public:
  SchwingerHamiltonianSz(int N,double mu,double x,double l0,double alpha,
			 double z,
			 double offset=0,double nu=0.,double y=1.,
			 double scale=1.,double Stot=0.);
  ~SchwingerHamiltonianSz(){};
  void setStarget(double Stot_);
protected:
  void initHMPO();

};

#endif // SCHWINGERZ_H
