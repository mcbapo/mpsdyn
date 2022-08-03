#ifndef SINGLEBATHHAMIL_H
#define SINGLEBATHHAMIL_H

#include <vector>
#include "ThermofieldHamiltonian.h"

/** 
    A Hamiltonian containing only nearest-neighbor terms, for the
    particular case of the bath-mapped-to-chain simulations
    (asthermofield, with T=0 bath). This means we have
    a spin site in the leftmost of the chain, and L bosonic modes 
    coupled to it on the left side.
    We assume couplings of the form
    \f[b_{1,0} L\left(d_0^{\dagger}+d_0\right)\f]
    for the coupling between the spin and the first bosonic mode on 
    the left and
    \f[b_{1 i}
    \left(d_{i+1}^{\dagger}d_i+d_{i}^{\dagger}d_{i+1}\right)\f].
    The spin system will also have a free Hamiltonian
    \f$H_{\mathrm{spin}}\f$. 
    The physical dimensions can be different for each site, if we
    allow different maximum occupation number for the bosons.
*/

class SingleBathHamiltonian:
public ThermofieldHamiltonian{
 protected:
	//int L;
	//std::vector<int> dims; // Physical dimensions (N+1 for bosons, 2
			    // for spin)
	//double g; // Coupling constant
	//std::vector<double> Bn; // Coefficients of the bath
				  // couplings (B1)
	//std::vector<double> An; // Coefficients of the bath local
				   // terms (A1)
	//mwArray OL; // Operator for the coupling with the 0 mode bosons
	//mwArray Hfree; // Free Hamiltonian of the spin 
	//mwArray Zbos;

	//MPO hamil;

public:
/** Create with all the required fields as arguments.  */
	SingleBathHamiltonian(int L,std::vector<int> dims,
			       const std::vector<double>& an,
			       const std::vector<double>& bn,
			       mwArray Hfree,mwArray OL);
	/** TODO: Create from file */
	//	SingleBathHamiltonian(const char* filename);
	~SingleBathHamiltonian();
	
	/** Return the length */
	int getLength() const { cout<<"SingleBathHamiltonian::getLength()"<<endl;return L+1;}
	int getSpinPosition() const {return 0;}	
	
	/** Return the vector of dimensions in the given container. */
	void getDimensions(std::vector<int>& dims_) const{dims_=dims;}
	
	bool hasHMPO() const{return false;}
	const MPO& getHMPO() const{return hamil;}
	bool hasUMPO() const {return true;}
	/** Construct the MPO corresponding to the evolution operator with
	    the splitting even-odd in the Hamiltonian. Take as argument the
	    complex factor in the exponential, delta.  For real time
	    evolution, delta is -i*tau, for imaginary time, it should be
	    -tau.
	    The even term includes boson pairs 01, 23,... AND the free
	    Hamiltonian of the central spin. The odd term includes the
	    bosonic pairs 12,34... AND both coupling terms between spin and
	    zero-modes, which we assume commute.*/
	void getExponentialMPOeven(MPO& expHe,complex_t delta);
	void getExponentialMPOodd(MPO& expHo,complex_t delta);
	
 protected:
	
	virtual void initHMPO();
	virtual void initZ();
	
 protected:

	/** Return the positions occupied by a two-body bosonic term,
	 * identified by the boson number k and k+1, meaning k+
	 * either on the left or right chain.*/
	shrt::Indices getTwoBodyTermPositions(int k) const;

	/** Compute the two body term in two bosonic modes k and k+1 on the
	 * right, or on the left (if left==true). Return it as a tensor with
	 * dimensions (d1xd2)^2, where d1 is the physical dimension
	 * (i.e. N-1) of the boson (k or k+1) that is to the left in
	 * the chain. */
	virtual void computeTwoBodyTermBath(mwArray& result,int k) const;
	
	/** */
	void split2term(mwArray& Ol,mwArray& Or,mwArray& expH,int dimL,int dimR,int& nr) 
	  const;
	/** Prepare the even-odd decomposition of the evolution operator:
	    returns the two terms to be inserted in the MPO for the 
	    exp(delta H12) (delta a complex number).
	    The exponential of the term starting in pos is returned.  
	*/
	virtual void getTwoBodyTermBathExponential(mwArray& Ol,mwArray& Or,complex_t delta,
						   int k) const;
	
	/** Construct the two-body MPO operator for the exponential of the
	 * coupling term between central spin and zero boosonic mode
	 * onits right. */
	virtual void getCouplingExponential(mwArray& Ospin,mwArray& Or,
					    complex_t delta) const;
	/** Compute the exponential of the single body term Hs */
	void getFreeSpinExponential(mwArray& Op,complex_t delta) const;
	
	virtual void setSingleModeTerm(MPO& expH,int pos,int dim,
				       complex_t delta) const;
};

#endif // SINGLEBATHHAMIL_H
