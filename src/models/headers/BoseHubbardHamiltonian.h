#ifndef BOSEHUB_H
#define BOSEHUB_H

#include "Hamiltonian.h"
#include <vector>

const std::vector<double> NULL_VEC; // default argument

  /** Bose-Hubbard Hamiltonian for a chain of 
      length L, with maximum occupation number N (this means the local 
      physical dimension is N+1)
      and parameters t, U, V  such that
      \f[
      H=- \sum_i t_i (b_i^{\dagger} b_{i+1} + h.c.) + 
      \sum_i \frac{U_i}{2} n_i (n_i-1) 
      \sum_i V_i n_i
      +\lambda_N \left(\sum_i n_i-N_{\mathrm{target}}\right)^2
      \f]
      It includes a penalty term to force a certain total number of particles.
  */

class BoseHubbardHamiltonian: 
public Hamiltonian{
 protected:
  int L,N;
  //double t,U,mu,V;
  MPO hamil;
  mwArray Z; 
  // Notice: the value U kept in here is actually U/2
  std::vector<double> tis,Uis,Vis;
  int Ntarg;double penN;

 public:
  BoseHubbardHamiltonian(int L,int N,double t,double U,int Ntarget,double penaltyN,
			 double V=0.);

  /** Idem but for position dependent parameters t,U and V.
      TODO: This does not work if they are not all the same length!! Allow
      only one to be varying with the site!*/
  BoseHubbardHamiltonian(int L,int N,double* tis,double* Uis,int Ntarget,double penN,
			 double* Vis=0);
  /** Idem with vector arguments. If some argument is shorter, the
      values are used until the last one and this is repeated to fill
      the chain (value zero has to be given explicitly!). */
  BoseHubbardHamiltonian(int L,int N,const std::vector<double>& tis,
			 const std::vector<double>& Uis,int Ntarget,double penN,
			 const std::vector<double>& Vis=NULL_VEC);

  ~BoseHubbardHamiltonian(){};
  const MPO& getHMPO() const {return hamil;}

  bool hasHMPO() const {return true;}

  int getBondDimension() const {
	  return 5;
  }

  /** Construct the MPO corresponding to the evolution operator with
      the splitting even-odd in the Hamiltonian. Take as argument the
      complex factor in the exponential, delta.  For real time
      evolution, delta is -i*tau, for imaginary time, it should be
      -tau. */
  void getExponentialMPOeven(MPO& expHe,complex_t delta);
  void getExponentialMPOodd(MPO& expHo,complex_t delta);

 private:
  /** Prepare the even-odd decomposition of the evolution operator:
      returns the two terms to be inserted in the MPO for the 
      exp(delta H12) (delta a complex number).
      TI is not assumed and the exponential of the term
      starting in pos is returned.  
  */
  void getTwoBodyTermExponential(mwArray& Ol,mwArray& Or,complex_t delta,
				 int pos=0);
  /** Compute the two-body term on positions (\param<pos>,\param<pos>+1):
      \f[
      h_{pos,pos+1}=-t_{pos} (b_{pos}^{\dagger} b_{pos+1} + h.c.) + 
      \frac{U_{pos}}{2} n_{pos} (n_{pos}-1) + 
      \frac{U_{pos+1}}{2} n_{pos+1} (n_{pos+1}-1) 
      \frac{V_{pos}-\mu_{pos}}{2} n_{pos} + 
      \frac{V_{pos+1}-\mu_{pos+1}}{2} n_{pos+1} 
      \f]
*/
  void computeTwoBodyTerm(mwArray& result,int pos) const;
  
  void initZ();
  void initHMPO();
  //void initHMPO(double* tis=0,double* Uis=0,double* muis=0,
  //	 double* Vis=0);
  /** Common computation for recovering the even and odd parts of the exp */
  void getExponentialMPO(MPO& expHo,complex_t delta,bool even);

};


#endif // BOSEHUB_H
