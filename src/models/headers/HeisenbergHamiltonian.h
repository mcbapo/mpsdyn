#ifndef HEISENBERG_H
#define HEISENBERG_H

#include "Hamiltonian.h"
#include <vector>

/** 
    Implementation of Heisenberg Hamiltonian with Z magnetic field and
    site-dependent coefficients.
    \f[
    H=\sum_{i=1}^{L-1} \left ( J_x^{[i]} S_x^{[i]}S_x^{[i+1]} +
    J_y^{i} S_y^{[i]}S_y^{[i+1]} +
    J_z^{i} S_z^{[i]}S_z^{[i+1]} \right )+
    \sum_{i=1}^{L} h_x^{i} S_x^{[i]}
    \sum_{i=1}^{L} h_y^{i} S_y^{[i]}
    \sum_{i=1}^{L} h_z^{i} S_z^{[i]}
    \f]
    Notice it is defined in terms of spin matrices, not sigmas.
*/

class HeisenbergHamiltonian:
public Hamiltonian{
protected:
  int L; // chain length
  int d; // dimension (only 2 supported)
  MPO hamil;
  mwArray Z;
  std::vector<double> Jxi,Jyi,Jzi,hi,hxi,hyi;

 public:
  /** Create a Heisenberg Hamiltonian with site-dependent
      coefficients, and dimension of the individual sites d (only d=2,
      i.e. S=1/2, supported). With only one vector for the magnetic
      field, Z direction is assumed.  An offset can be given, so that
      the MPO corresponds to H+offset.  */
  HeisenbergHamiltonian(int L,std::vector<double> Jx,std::vector<double> Jy,
			std::vector<double> Jz,std::vector<double> h,int d=2,double offset=0.);
  /** Create with constant Js and only hi changing with the site */
  HeisenbergHamiltonian(int L,double Jx,double Jy,double Jz,
			std::vector<double> h,int d=2,double offset=0.);
  /** Create with constant coefficients */
  HeisenbergHamiltonian(int L,double Jx,double Jy,double Jz,double h,int d=2,double offset=0.);

  /** More general version, admitting local X and Y fields, too. */
  HeisenbergHamiltonian(int L,std::vector<double> Jx,std::vector<double> Jy,
			std::vector<double> Jz,std::vector<double> hx,
			std::vector<double> hy,std::vector<double> hz,int d=2,double offset=0.);
  /** Create with constant Js and only hi changing with the site */
  HeisenbergHamiltonian(int L,double Jx,double Jy,double Jz,
			std::vector<double> hx,std::vector<double> hy,
			std::vector<double> hz,int d=2,double offset=0.);
  /** Create with constant coefficients */
  HeisenbergHamiltonian(int L,double Jx,double Jy,double Jz,double hx,double hy,double hz,
			int d=2,double offset=0.);

  
  ~HeisenbergHamiltonian();

  const MPO& getHMPO() const {return hamil;}

  bool hasHMPO() const {return true;}

  /** Construct the MPO corresponding to the evolution operator with
      the splitting even-odd in the Hamiltonian. Take as argument the
      complex factor in the exponential, delta.  For real time
      evolution, delta is -i*tau, for imaginary time, it should be
      -tau. */
  virtual void getExponentialMPOeven(MPO& expHe,complex_t delta) const ;
  virtual void getExponentialMPOodd(MPO& expHo,complex_t delta) const ;

  /** Idem for the evolution of mixed states */
  virtual void getDoubleExponentialMPOeven(MPO& expHe,complex_t delta) const ;
  virtual void getDoubleExponentialMPOodd(MPO& expHo,complex_t delta) const ;

  virtual void getExtendedExponentialMPOeven(MPO& expHe,complex_t delta) const ;
  virtual void getExtendedExponentialMPOodd(MPO& expHo,complex_t delta) const ;

  /** Construct an MPo for the commutator with an operator, and return in in the argument.  */
  void getCommutatorMPO(MPO& mpo);

 protected:
  /** Special contructor for the TI daughter with ancillary chain */
  HeisenbergHamiltonian(int L,int d=2);

  /** Prepare the even-odd decomposition of the evolution operator:
      returns the two terms to be inserted in the MPO for the 
      exp(delta H12) (delta a complex number).
      TI is not assumed and the exponential of the term
      starting in pos is returned.  
      TODO: why is this method (and the following one) public??
  */
  void getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,
				 int pos=0) const;
  void getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,
				 const mwArray& H12,complex_t delta) const;
  /** Compute the two-body term on positions (\param<pos>,\param<pos>+1):
      \f[
      h_{pos,pos+1}=J_x({pos}) S_x^{pos} S_x^{pos+1} +
      J_y({pos}) S_y^{pos} S_y^{pos+1}+
      J_z({pos}) S_z^{pos} S_z^{pos+1}+
      \frac{h({pos})}{2} S_z^{pos} +
      \frac{h({pos+1})}{2} S_z^{pos+1}
      \f]
      (except for edges, where the \f$h(0),\ h(L-1)\f$ terms are not divided by two).
*/
  void computeTwoBodyTerm(mwArray& result,int pos) const;
 protected:
  /** As the one before, but passing directly the values of the
      required parameters, to be used by a daughter class which is
      TI. */
  void computeTwoBodyTerm(mwArray& result,double Jx,double Jy,double Jz,double h1,double h2) const;

  void computeTwoBodyTerm(mwArray& result,double Jx,double Jy,double Jz,double hx1,double hx2,
			  double hy1,double hy2,double hz1,double hz2) const;

 private:
  void initZ();
  void initHMPO(double offset=0.);
  //void initHMPO(double* tis=0,double* Uis=0,double* muis=0,
  //	 double* Vis=0);
  /** Common computation for recovering the even and odd parts of the exp */
  void getExponentialMPO(MPO& expHo,complex_t delta,bool even) const;
  void getDoubleExponentialMPO(MPO& expHo,complex_t delta,bool even) const;
  void getExtendedExponentialMPO(MPO& expHo,complex_t delta,bool even) const;

  // This is copied from IsingHamiltonian, but it would be better to
  // have it somewhere else (common single place)
  /** 
      Auxiliary function for the construction of the commutator
      MPO. If this should be more general (not only for this
      IsingHamiltonian case), these functions might need to be
      somewhere else. Specially the double
   */
  void initZdoub(mwArray& Z);


};


#endif // HEISENBERG_H
