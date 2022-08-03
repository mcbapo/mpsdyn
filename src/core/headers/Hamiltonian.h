/**
   \file Hamiltonian.h 
   Definition of the basic interface of a Hamiltonian.
   
   \author Mari-Carmen Banuls
   \date 16/03/2011
*/

#ifndef HAMIL_H
#define HAMIL_H

#include "Operator.h"
#include <map>
#include <string>
#include <iostream>

#include "MPO.h"

/** 
    Basic interface for a Hamiltonian which encapsulates methods 
    to construct an MPO for energy minimization or time evolution. 
*/

class Hamiltonian{

 public: 
  virtual ~Hamiltonian(){};
  virtual bool hasHMPO() const =0;

  /** This method should be declared const!!! */
  virtual const MPO& getHMPO() const=0;
  
  /** 
     Whether the Hamiltonian knows how to construct the evolution 
     MPO (default is false) 
   */
  virtual bool hasUMPO() const {return false;}

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
    std::cout<<"ERROR: getBasicUOperators called from basic Hamiltonian"
	     <<std::endl;
    exit(212);
  };

  /** When the MPO for the unitaries is available, return the number
      of different terms per time step. This depends on the
      Hamiltonian itself and on the Trotter order. */
  virtual int getUOpersPerStep(int orderT) const{
    std::cout<<"ERROR: getUOpersPerStep called from basic Hamiltonian"
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
    std::cout<<"ERROR: getUOperDimensions called from basic Hamiltonian"
	     <<std::endl;
    exit(212);
  }
};



#endif //HAMIL_H
