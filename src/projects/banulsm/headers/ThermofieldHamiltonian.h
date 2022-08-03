#ifndef THERMOFIELDHAMIL_H
#define THERMOFIELDHAMIL_H

#include <vector>
#include "Hamiltonian.h"

/** 
    A Hamiltonian containing only nearest-neighbor terms, for the
    particular case of the thermofield simulations. This means we have
    a spin site in the center of the chain, and L bosonic modes 
    coupled to it on each side.
    We assume couplings of the form
    \f[b_{2,0} L\left(c_0^{\dagger}+c_0\right)\f]
    for the coupling between the spin and the first bosonic mode on the
    right (similar for the left, with \f$b_{1,0}\f$
    and modes \f$d_0,\,d_0^{\dagger}\f$) and
    \f[b_{1 i}
    \left(d_{i+1}^{\dagger}d_i+d_{i}^{\dagger}d_{i+1}\right)\f]
    \f[- b_{2 i}
    \left(c_{i+1}^{\dagger}c_i+c_{i}^{\dagger}c_{i+1}\right)\f]
    for the right bosonic chain.
    The spin system will also have a free Hamiltonian
    \f$H_{\mathrm{spin}}\f$. 
    The physical dimensions can be different for each site, if we
    allow different maximum occupation number for the bosons.
*/

class ThermofieldHamiltonian:
public Hamiltonian{
 protected:
	int L;
	std::vector<int> dims; // Physical dimensions (N+1 for bosons, 2
			    // for spin)
	//double g; // Coupling constant
	std::vector<double> B2n; // Coefficients of the right couplings (B2)
	std::vector<double> B1n; // Coefficients of the left
				  // couplings (B1)
	std::vector<double> A2n; // Coefficients of the right local
				   // terms (A2)
	std::vector<double> A1n; // Coefficients of the left local
				   // terms (A1)
	mwArray OL; // Operator for the coupling with the 0 mode bosons
	mwArray Hfree; // Free Hamiltonian of the spin 
	mwArray Zbos;

	MPO hamil;

public:
/** Create with all the required fields as arguments.  */
	ThermofieldHamiltonian(int L,std::vector<int> dims,
			       const std::vector<double>& a1n,
			       const std::vector<double>& a2n,
			       const std::vector<double>& b1n,
			       const std::vector<double>& b2n,
			       mwArray Hfree,mwArray OL);
	/** TODO: Create from file */
	//	ThermofieldHamiltonian(const char* filename);
	~ThermofieldHamiltonian();
	
	/** Return the length */
	virtual int getLength() const {cout<<"ThermofieldHamiltonian::getLength()"<<endl;return 2*L+1;}

	/** Return the position where the system is */
	virtual int getSpinPosition() const {return L;}	

	/** Return the vector of dimensions in the given container. */
	virtual void getDimensions(std::vector<int>& dims_) const{dims_=dims;}
	
	virtual bool hasHMPO() const{return false;}
	virtual const MPO& getHMPO() const{return hamil;}
	virtual bool hasUMPO() const {return true;}
	/** Construct the MPO corresponding to the evolution operator with
	    the splitting even-odd in the Hamiltonian. Take as argument the
	    complex factor in the exponential, delta.  For real time
	    evolution, delta is -i*tau, for imaginary time, it should be
	    -tau.
	    The even term includes boson pairs 01, 23,... AND the free
	    Hamiltonian of the central spin. The odd term includes the
	    bosonic pairs 12,34... AND both coupling terms between spin and
	    zero-modes, which we assume commute.*/
	virtual void getExponentialMPOeven(MPO& expHe,complex_t delta);
	virtual void getExponentialMPOodd(MPO& expHo,complex_t delta);
	
 protected:
	
	virtual void initHMPO();
	virtual void initZ();
	
 protected:
	/** Constructor for daughter class */
	ThermofieldHamiltonian(int L_,std::vector<int> dims,const vector<double>& A1n_,
			       const vector<double>& A2n_,const vector<double>& B1n_,
			       const vector<double>& B2n_);

	/** Return the positions occupied by a two-body bosonic term,
	 * identified by the boson number k and k+1, either on the left or
	 * right chain.*/
	shrt::Indices getTwoBodyTermPositions(int k,bool left=false) const;

	/** Compute the two body term in two bosonic modes k and k+1 on the
	 * right, or on the left (if left==true). Return it as a tensor with
	 * dimensions (d1xd2)^2, where d1 is the physical dimension
	 * (i.e. N-1) of the boson (k or k+1) that is to the left in the chain. */
	virtual void computeTwoBodyTermBath(mwArray& result,int k,bool left=false) const;
	
	/** */
	void split2term(mwArray& Ol,mwArray& Or,mwArray& expH,int dimL,int dimR,int& nr) 
	  const;
	/** Prepare the even-odd decomposition of the evolution operator:
	    returns the two terms to be inserted in the MPO for the 
	    exp(delta H12) (delta a complex number).
	    The exponential of the term starting in pos is returned.  
	*/
	virtual void getTwoBodyTermBathExponential(mwArray& Ol,mwArray& Or,complex_t delta,
						   int k,bool left=false) const;
	
	/** Construct the three-body MPO operators for the exponential of the
	 * coupling terms between central spin and zero boosonic modes. */
	virtual void getCouplingExponential(mwArray& Ol,mwArray& Ospin,mwArray& Or,
					    complex_t delta) const;
	/** Compute the exponential of the single body term Hs */
	void getFreeSpinExponential(mwArray& Op,complex_t delta) const;
	
	virtual void setSingleModeTerm(MPO& expH,int pos,int dim,
				       complex_t delta,bool left) const;
};

#endif // THERMOFIELDHAMIL_H
