/**
   \file Liouvillian.h 
   Definition of the basic interface of a Liouvillian.
   
   \author Mari-Carmen Banuls
   \date 7/02/2013
*/

#ifndef LIOUV_H
#define LIOUV_H

#include "Operator.h"
#include <map>
#include <string>
#include <iostream>

#include "MPO.h"

/** 
    Basic interface for a Liouvillian superoperator, coming from a
    Lindblad master equation
    \f[
    \frac{d\rho}{d t}=i[\rho,H]+\sum_{\alpha}\gamma_{\alpha}\left( L_{\alpha}
    \rho L_{\alpha}^{\dagger}-\frac{1}{2} \rho L_{\alpha} L_{\alpha}^{\dagger}
    -\rho L_{\alpha} L_{\alpha}^{\dagger} \rho\right)
    \f]
    for some short-range Hamiltonian and dissipation terms
    \f$L_{\alpha}\f$.  This interface encapsulates methods to
    construct an MPO for the superoperator, that can be used to find
    the steady state. It could in principle also implement time
    evolution.
    There is no generic implementation, and any particular model will
    need to implement its own MPO for L, etc.
 */

class Liouvillian{

 public: 

  /** Get the Liouvillian superoperator as an MPO, the physical
      dimensions of which will correspond to double sites, with the
      first one appearing corresponding to the output index of
      \f$\rho\f$ (terms acting on the left in the master equation). */
  virtual const MPO& getLMPO()=0;

  /** Get the adjoint of the Liouvillian superoperator */
  virtual const MPO& getLAdjointMPO()=0;

  
  /** 
     Whether the Liouvillian knows how to construct the evolution 
     MPO (default is false) 
   */
  virtual bool hasEvolutionMPO() const {return false;}

  /** Calculate the unitary operators corresponding to a certain 
      step of the time evolution, starting at initT and lasting 
      delta, using Trotter expansion to orderT, for position \param<pos>
      \todo different in the infinite case
      The different operators that result are inserted in the vector
      \param<opSet>, with the first one (opSet[0]) corresponding to
      the operator that acts first on the state. The vector must
      already contain Operator reference in the appropriate number.
      \param<initT> is the time at the beginning of the evolution step.
      \param<imag> is true when the imaginary time evolution is to be 
      computed.
      \param<transverse> indicates whether the standard ordering (uldr)
      or the one for the transverse operators (ruld) is returned.
*/
  virtual void getBasicUOperators(int pos,double initT,double delta,int orderT,
				  vector<Operator&>& opSet,
				  bool imag=false,bool transverse=true){
    std::cout<<"ERROR: getBasicUOperators called from basic Liouvillian"
	     <<std::endl;
    exit(212);
  };

  /** When the MPO for the unitaries is available, return the number
      of different terms per time step. This depends on the
      Liouvillian itself and on the Trotter order. */
  virtual int getUOpersPerStep(int orderT) const{
    std::cout<<"ERROR: getUOpersPerStep called from basic Liouvillian"
	     <<std::endl;
    exit(212);
  }

  /** If unitary MPO are available, return the dimensions of 
      a given Operator in the set. If \param<opernr>=-1, the last 
      one is assumed. This will determine how TransverseManager
      can pack operators from different steps to make a more 
      efficient MPO.*/
  virtual shrt::Indices getUOperDimensions(int orderT,int opernr=-1,
					   bool transverse=true) const{
    std::cout<<"ERROR: getUOperDimensions called from basic Liouvillian"
	     <<std::endl;
    exit(212);
  }

  virtual ~Liouvillian(){};
};



#endif //LIOUV_H
