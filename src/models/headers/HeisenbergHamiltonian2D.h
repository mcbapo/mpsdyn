#ifndef HEISENBERG2D_H
#define HEISENBERG2D_H

#include "Hamiltonian.h"
#include <vector>

/** 
    Implementation of Heisenberg Hamiltonian with local magnetic field and
    translationally invariant coefficients, on a 2D system with a ladder geometry.
    \f[
    H=\sum_{i=1}^{L-1} \left ( J_x S_x^{[i]}S_x^{[i+1]} +
    J_y S_y^{[i]}S_y^{[i+1]} +
    J_z S_z^{[i]}S_z^{[i+1]} \right )+
    \sum_{i=1}^{L} h_x S_x^{[i]}
    \sum_{i=1}^{L} h_y S_y^{[i]}
    \sum_{i=1}^{L} h_z S_z^{[i]}
    \f]
    Notice it is defined in terms of spin matrices, not sigmas.
*/

class HeisenbergHamiltonian2D:
public Hamiltonian{
protected:
  int Lx,Ly; // total dimensions of the system (assume Ly<=Lx)
  int nL; // nr of spin systems grouped per site (at most)
  int d; // dimension (2)
  int legNr; // number of legs in the ladder (determines the range of the interactions)
  MPO hamil;
  mwArray Z;
  double Jx,Jy,Jz,hx,hy,hz;
  
 public:
  /** Create a Heisenberg Hamiltonian on a ladder geometry.  */
  HeisenbergHamiltonian2D(int Lx,int Ly,int nL,double Jx,double Jy,double Jz,double hx,double hy,
			  double hz,double offset=0.);
  
  ~HeisenbergHamiltonian2D();

  int getLegNr(){return legNr;}

  const MPO& getHMPO() const {return hamil;}

  bool hasHMPO() const {return true;}

  /** Because in this geometry, the Hamiltonian is no longer
      nearest-neighbour, there are more than two terms in the Trotter
      deccomposition of the evolution operator. There will be two
      (even-odd) horizontal terms and one (if there are only two legs)
      or two vertical ones. The following functions construct them.
      They take as argument the complex factor in the exponential,
      delta.  For real time evolution, delta is -i*tau, for imaginary
      time, it should be -tau. */
  virtual void getExponentialMPOevenH(MPO& expHe,complex_t delta) const ;
  virtual void getExponentialMPOoddH(MPO& expHo,complex_t delta) const ;
  virtual void getExponentialMPOevenV(MPO& expHe,complex_t delta) const ;
  virtual void getExponentialMPOoddV(MPO& expHo,complex_t delta) const ;

  /** Idem for the evolution of mixed states */
  virtual void getDoubleExponentialMPOevenH(MPO& expHe,complex_t delta) const ;
  virtual void getDoubleExponentialMPOoddH(MPO& expHo,complex_t delta) const ;
  virtual void getDoubleExponentialMPOevenV(MPO& expHe,complex_t delta) const ;
  virtual void getDoubleExponentialMPOoddV(MPO& expHo,complex_t delta) const ;

  virtual void getExtendedExponentialMPOevenH(MPO& expHe,complex_t delta) const ;
  virtual void getExtendedExponentialMPOoddH(MPO& expHo,complex_t delta) const ;
  virtual void getExtendedExponentialMPOevenV(MPO& expHe,complex_t delta) const ;
  virtual void getExtendedExponentialMPOoddV(MPO& expHo,complex_t delta) const ;


 protected:

  /** Prepare the even-odd decomposition of the evolution operator:
      returns the two terms to be inserted in the MPO for the 
      exp(delta H12) (delta a complex number).
      TI is not assumed and the exponential of the term
      starting in pos is returned.  
  */
  void getTwoBodyTermExponentialx(mwArray& Ol,mwArray& Or,complex_t delta,
				 int posX=0) const;
  void getTwoBodyTermExponentialy(mwArray& Ol,mwArray& Or,complex_t delta,
				 int posY=0) const;

  /* General method to expoonentiate and split the two body operator. */
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
  void computeTwoBodyTermHorizontal(mwArray& result,int pos) const;
  void computeTwoBodyTermVertical(mwArray& result,int pos) const;

 private:
  void initZ();
  void initHMPO(double offset=0.);
  //void initHMPO(double* tis=0,double* Uis=0,double* muis=0,
  //	 double* Vis=0);
  /** Common computation for recovering the even and odd parts of the exp in horizontal and vertical directions */
  void getExponentialMPOx(MPO& expHo,complex_t delta,bool even) const;
  void getExponentialMPOy(MPO& expHo,complex_t delta,bool even) const;
  void getDoubleExponentialMPOx(MPO& expHo,complex_t delta,bool even) const;
  void getDoubleExponentialMPOy(MPO& expHo,complex_t delta,bool even) const;
  void getExtendedExponentialMPOx(MPO& expHo,complex_t delta,bool even) const;
  void getExtendedExponentialMPOy(MPO& expHo,complex_t delta,bool even) const;


  void computeTwoBodyTerm(mwArray& result,double Jx,double Jy,double Jz,double hx1,double hx2,
			  double hy1,double hy2,double hz1,double hz2) const;

  
};


#endif // HEISENBERG2D_H
