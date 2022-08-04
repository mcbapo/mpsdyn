/**
   \file Contractor.h Definition of the class that provides the basic
   functionality with MPS-MPO operations, and the main (finite) MPS
   algorithms.
   
   \author Mari-Carmen Banuls
   \date 16/03/2011
*/

#ifndef CONTRACTOR_H
#define CONTRACTOR_H

//#include "libdiminf.h"

#include "MPS.h"
#include "Operator.h"
#include "MPO.h"
#include "TmpManager.h"
#include "TmpManagerExc.h"
#include "TensorMultiplier.h"

#include <cmath>
#include <vector>
#include <string>

/** 
    Type of eigensolver to be used. 
    \todo Include the possibility to select the method in the case of PRIMME
*/
enum EigenSolver {
  arpack=0, // default
  fulleig, // Use eig and find full spectrum
  primme, // primme DYNAMIC (default)
  primme_JDQR, // primme with method JDQR
  primme_arnoldi, // primme with method Arnoldi
};


/** 
    Fundamental class that organizes the contractions to optimize the MPS 
    that best approximates the effect of applying an operator row onto 
    a given MPS.
*/

class Contractor{
  double svdtol;
  double convtol;
  double eigtol;
  double excorthtol;// tolerance for lack of orthogonality of excitations 
  EigenSolver solver;
  
 public:
  ~Contractor();

  static Contractor& theContractor();

  /** Finds an optimal MPS approximation to the action of the MPO on
      the original MPS (ket). The bond dimension used will be that of
      the initial MPS(default, equal to orig).  The result is returned
      inside init. This could be empty, in which case, a copy of the
      original ket is used as a starting point (after adjusting
      dimensions and gauge).
      If a value D!=0 is given, it is the value used for the maximum
      bond of the final MPS, provided it is larger or equal than the
      max bond of orig.
      If errN!=0 is given, this methods computes also the 
      (unnormalized) truncation error,
      \f$*errN=\| result - Ops ket \|^2\f$
      If errT!=0 is also given, then it returns
      *errT=norm1=norm(ops ket)^2. 
      WARNING: This is a costly calculation, and should be spared,
      if possible.
  */ 
  void optimize(const MPO& ops,const MPS& ket,MPS& init,int D=0,
		double* errN=0,double* errT=0);

  /** Find an optimal (bra) MPS approximation to the LEFT action of the
      MPO on the original (bra) MPS. As before, if no initial
      MPS is given, a copy of bra is taken as starting point. 
      IMPORTANT: the state bra on which the operator acts has to be 
      given as ket (i.e., it will be conjugated within the method) 
      and the result is also given as a ket (the return value is thus 
      |Psi> s.t. <Psi|~<bra| Ops)
      If a value D!=0 is given, it is the value used for the maximum
      bond of the final MPS, provided it is larger or equal than the max 
      bond of orig.
      If errN!=0 is given, this method computes also the truncation error 
      *errN=|bra ops - result|^2
      If errT!=0 is also given, then it returns
      *errT=norm1=norm(bra ops)^2. 
      See optimize() for details.
 */
  void optimizeL(const MPO& ops,const MPS& bra,MPS& init,int D=0,
		 double* errN=0,double* errT=0);

  /** Finds an MPS approximation to the action of the inverse of the
      given MPO on the original MPS (ket). To this end, 
      /f$ \| ops |init\rangle - |ket \rangle \| /f$
      is optimized. The bond dimension used will be 
      that of the initial MPS(default, equal to orig).
      MPS init is taken as a starting point for the optimization, 
      and at the end contains the result.
      If a value D!=0 is given, it is the value used for the maximum
      bond of the final MPS, provided it is larger or equal than the max 
      bond of orig.
      If errN!=0 is given, this methods computes also the 
      (unnormalized) truncation error,
      \f$*errN=\| ops result - ket \|^2\f$
      If errT!=0 is also given, then it returns
      *errT=norm1=norm(ops ket)^2. 
      WARNING: This is a costly calculation, and should be spared,
      if possible.
  */ 

  void optimizeInv(const MPO& ops,const MPS& ket,MPS& init,int D=0,
		double* errN=0,double* errT=0);

  /** 
      As the former one, but optimizes MPS init such that
      /f$ \| (ops-i \eta )|init \rangle - | ket \rangle \| /f$
      is minimal.
  */
  void optimizeResolvent(const MPO& ops,const MPS& ket,MPS& init,
			 double eta,int D=0,
			 double* errN=0,double* errT=0);
  /** 
      Optimize the MPS init such that it is the best MPS approximation
      (with bond dimension D) to a sum of MPS (in vector kets) with
      complex coefficients (in vector beta), /f$\sum_k
      \beta_k|ket_k\rangle /f$
  */
  void optimizeSum(vector<const MPS*> kets,std::vector<complex_t> beta,
		MPS& init,int D=0);

  /** A particular case of optimizeSum is the optimization of a single
      MPS as a MPS of smaller bond dimension. This method is just a
      simpler interface to the one above, with the option to estimate
      the error, if the double pointer err is provided. */
  void optimizeMPS(const MPS& ket,MPS& init,int D=0,double* err=0);

  /** Orthogonalize the initial MPS (init, containing the initial
      vector and chenged afterwards) wrt a list of vectors, by trying
      to minimize the function
      \f[\||\Psi \rangle -| \Phi_0\rangle \|^2+\sum_i \lambda_i \langle \Psi | v_i \rangle \langle v_i |\Psi \rangle \f] 
      It requires values for the penalty terms associated to each of
      the vectors in the list.
*/
  void orthogonalize(MPS& init,int D,
		     const vector<const MPS*>& vectors,
		     const vector<double>& penalties);

  /** Similar to \method<orthogonalize>, but the vector to be
      approximated is the result of an MPO acting on a MPS. */
  void orthogonalize(MPS& init,int D,const MPO& ops,const MPS& orig,
		     const vector<const MPS*>& vectors,
		     const vector<double>& penalties);

  /** As \method<optimizeSum>, but at the same time that the distance
      to a sum of vectors is minimized, a constraint for orthogonality
      w.r.t. other vectors is enforced by means of penalty terms, in
      the same way as in \method<orthogonalize>. */
  void orthogonalizeSum(MPS& init,int D,
			const vector<const MPS*>& kets,const vector<complex_t>& beta,
			const vector<const MPS*>& vectors,
			const vector<double>& penalties);

  /** Contracts two MPS and returns the scalar product (a complex
      value).  Optional argument dir indicates the
      direction to apply the contraction (from left to right, if 'L',
      or right to left if 'R').  Both ket and bra are actually to be
      given as kets, and bra will be appropriately conjugated for the
      computation. */
  complex_t contract(const MPS& ket,const MPS& bra,char dir='L');
  /** Same as before, but contract only part of the chain, until site
      pos (not included in the contraction). It returns a mwArray,
      with dimensions DuxDd */
  mwArray contract(const MPS& ket,const MPS& bra,int pos,char dir='L');

  /** Contracts two MPS with an MPO in between and returns 
      the product (a complex value). */
  complex_t contract(const MPS& ket,const MPO& ops,const MPS& bra);
  // Just to check: do the same starting from right
  complex_t contractR(const MPS& ket,const MPO& ops,const MPS& bra);
  /** As in the contraction of two MPS, this version contracts only
      part of the chain from left (default) or right, until site pos
      (excluded), and returns a tensor with dimensions
      D(bra)*Xi*D(ket) */
  mwArray contract(const MPS& ket,const MPO& ops,const MPS& bra,int pos,char dir='L');

  
  /** Contracts the MPO acting on the ket with itselt, to get 
      <ket|O+ O|ket>. This is in general more expensive than 
      any other contraction in the program, therefore is to be used just 
      for error calculations, but not in a regular way. */
  complex_t contract2(const MPO& ops,const MPS& ket);
  /** Contracts the MPO acting on the bra to the left 
      with itselt, to get 
      <bra|O O+|bra>. This is in general more expensive than 
      any other contraction in the program, therefore is to be used just 
      for error calculations, but not in a regular way. */
  complex_t contract2(const MPS& bra,const MPO& ops);
  /** Contracts <bra|O+ O|ket>. This is in general more expensive than 
      any other contraction in the program, therefore is to be used just 
      for error calculations, but not in a regular way. */
  complex_t contract2(const MPS& bra,const MPO& ops,const MPS& ket);
  
  /** Set the cutoff for the inversion of N in the system of equations */
  void setSVDTol(double value){svdtol=value;}

  /** Set the criterion of convergence (relative variation lower than value) */
  void setConvTol(double value){convtol=value;}

  void setExcOrthTol(double value){excorthtol=value;}
  
  /** Return the convergence criteria in use */
  double getConvTol() const {return convtol;}
  double getEigTol() const {return eigtol;}
  double getSVDTol() const {return svdtol;}
  double getExcOrthTol() const {return excorthtol;}

  /** Set the tolerance for the eigensolver */
  void setEigTol(double value){eigtol=value;}

  /** Set the eigensolver to be used in findGroundState and
      findNextExcitedState methods. */
  void setEigenSolver(EigenSolver solver_=arpack){solver=solver_;}

  /** Get the name of the eigensolver in use */
  const string getSolverName() const;

  /** Find the best MPS approximation to the dominant right(left) 
      eigenvector of the given operator row, with maximum bond 
      dimension D. The pointer lambda is filled with the (absolute 
      value of) corresponding eigenvalue. The result is returned in 
      resultMPS, whose initial value is used as starting point. */
  void findRightEigenvector(const MPO& ops,int D,double* lambda,
			    MPS& resultMPS);
  void findLeftEigenvector(const MPO& ops,int D,double* lambda,
			   MPS& resultMPS);

/** Find the lowest eigenstate by Rayleigh-Ritz variational method:
    \f[
    min \frac{\langle \Psi| H|\Psi \rangle}{\langle \Psi| \Psi \rangle}.
    \f]
    To be used on an MPO which represents a Hamiltonian. 
    If nrLocalEig is given andit is not 1, it is the number of
    eigenvalues targeted by the underlying local optimizations.
If tmpfile is given, every greqSv rounds the temporary result is saved
*/
  void findGroundState(const MPO& ops,int D,double* lambda,
		       MPS& resultMPS,double offset=0.,int nrLocalEig=1,
		       const string& tmpfile="",int freqSv=200,bool useLM=false);

  /** As findGroundState, but at the level of the local eigensolver,
      use the option for largest magnitude eigenvalue, instead of the
      smallest, as this seems more stable. To get this working, we
      need to shift the H, such that the smallest eigenvalue has also
      the largest magnitude. This can be done by adding a negative
      offset, but it has to be large enough. Typically, for NN
      Hamiltonians, the order of the system size should work.

      So far, in most cases analyzed, this was not a problem, except
      for StochasticHamiltonian, so, although in pple it could be the
      best option in general, I leave it for now as an optional
      method, with std calls still using the same. */
  void findGroundStateLM(const MPO& ops,int D,double* lambda,
		       MPS& resultMPS,double offset=0.,int nrLocalEig=1,
		       const string& tmpfile="",int freqSv=200);

  /** Apply part os a sweep of GS search between sites pos0 and pos1.*/
  void sweepPart(int pos0,int pos1,const MPO& ops,int D,
		 MPS& init,double offset=0.,int knr=1);


  /** Try to find the eigenstate that is closest to the given target value. 
      \param<ops> (MPO) represents the Hamiltonian
      \param<D> (int) Maximum bond dimension to be used
      \param<lambda> (double*) place to store the computed eigenvalue
      \param<resultMPS> (MPS&) place to store the computed eigenvector
      \param<target> (double) Value around which we want to find eigenpairs
                       (default value 0.)
      \param<nrLocalEig> (int) Number of eigenvalues targetted by every local 
                         eigenvalue problem, during sweeps (default 1).
      \warning Notice that this method is much slower than the regular 
      \method<findGroundState> or \method<findNectExcitedState>, due to the 
      slow convergence of the underlying routines to find interior 
      (generalized) eigenvalues.
      \todo Add support for ARPACK solver!
 */
  void findClosestEigenstate(const MPO& ops,int D,double* lambda,  //const char* which,
			     MPS& resultMPS,double target=0.,int nrLocalEig=1);


  /** After determining the GS as MPS, we can look for a first excited
      state, as the minimum energy of H-Eo|GS><GS|.  This can be
      repeated, substracting more eigenstates, to find higher
      excitations.  Takes as arguments a vector containing the already
      computed excited states.
      \warning It will not work if the absolute energy of the excited
      state that is to be computed lies above zero. Thus, if this
      could be the case, it is advisable to specify an offset, so that
      the Hamiltonian handled will be H+offset*1 
  */
  void findNextExcitedState(const MPO& ops,int D,
			    const std::vector<MPS*>& computedLevels,
			    double* lambdak,
			    MPS& init,double offset=0.,int nrLocalEig=1,
			    const string& tmpfile="",
			    int maxTime=0,bool useLM=false);
  /** Same, using the larget magnitude value, instead of SM */
  void findNextExcitedStateLM(const MPO& ops,int D,
			    const std::vector<MPS*>& computedLevels,
			    double* lambdak,
			    MPS& init,double offset=0.,int nrLocalEig=1,
			    const string& tmpfile="",
			    int maxTime=0);

  /** The basic method to find excited states is to find a ground
      state of a given MPO plus a list of projectors with their
      respective penalties. In the case of the excited state, those
      are minus the energies, but it could be more general. */
  void findGroundStateWithProjectorPenalty(const MPO& ops,int D,
					   const std::vector<double>& penalties,
					   const std::vector<MPS*>& projected,
					   double* lambdak,
					   MPS& init,double offset=0.,int nrLocalEig=1,
					   const string& tmpfile="",
					   int maxTime=0,bool useLM=false);

  void findGroundStateWithProjectorPenaltyLM(const MPO& ops,int D,
					   const std::vector<double>& penalties,
					   const std::vector<MPS*>& projected,
					   double* lambdak,
					   MPS& init,double offset=0.,int nrLocalEig=1,
					   const string& tmpfile="",
					   int maxTime=0);

  /** As \method<sweepPart>, but including the penalty projector. */
  void sweepPartWithProjPenalty(int pos0,int pos1,const MPO& ops,int D,
				double mu,MPS& projected,
				MPS& init,double offset=0.,int nrLocalEig=1);

  /** As \method<sweepPart>, but for the next excited state iteration. */
  void sweepPartNextExcitedState(int pos0,int pos1,const MPO& ops,int D,
				 const std::vector<MPS*>& computedLevels,
				 MPS& init,double offset=0.,int nrLocalEig=1);

  /** As \method<sweepPart>, but including the penalty projector. */
  void sweepPartWithProjectorPenalty(int pos0,int pos1,const MPO& ops,int D,
				     const std::vector<double>& penalties,
				     const std::vector<MPS*>& projected,
				     MPS& init,double offset=0.,int nrLocalEig=1);

  /** As \method<findNextExcitedstate> but for the GEVP
      (H-target)(H-target)x=lambda(H-target)x */
  void findNextClosestEigenstate(const MPO& ops,int D,
				 const std::vector<MPS*>& computedLevels,
				 double* lambdak,
				 MPS& init,double target=0.,int nrLocalEig=1);


  /** Given a MPS, calculate the entropy of part of the chain, for
      the first \param <len> sites. If len==0, the entropy of
      half of the chain is computed.
  */
  double getEntropy(const MPS& state,int len=0);

  /** Given an MPS, compute the entropy of an internal subchain,
      containing sites from \param <posL> to \param <posR>, both 
      included (sites numbered from 0 to L-1). 
  */
  double getEntropy(const MPS& state,int posL,int posR);

/** Given an MPS, construct the reduced density matrix for site pos
 * (starting from 0) */
  mwArray getRDM(const MPS& state,int pos);

  /** Given a MPS, construct the RDM for sites posL-posR (both included) */
  mwArray getRDM(const MPS& state,int posL, int posR);

  /** Given a MPS, calculate the Schmidt values when the chain is cut 
      after the first \param<len> sites. If len==0, the middle cut is used.
  */
  void getSchmidtValues(const MPS& state,std::vector<complex_t>& lambdas,
			int len=0);

/**
   Given an MPS and an MPO, compute the effective operator for one
   site, by contracting MPS-MPO-MPS to the left and the right, and
   then contracting with the local operator term. Returns a mwArray
   of dimension (d x Dl x Dr )^2, as would be acting on the ket
 */
  mwArray getEffectiveOperatorSingleSite(const MPS& ket,const MPO& ops,int pos);
  
  /** As the previous method, but gets the effective operator for k
      adjacent sites starting on position pos */
  mwArray getEffectiveOperatorMultiSite(const MPS& ket,const MPO& ops,int k,int pos);

  /** Given an MPS and an MPO, return the effective operator for one
   site as a TensorMultiplier (which can then be used by eigs). The
   Multiplier will consist of the contracted MPS-MPO-MPS to the left
   and the right, and the local operator term in the center.  In
   general, if the dimensions are large, this will be more convenient
   than getting the full operator as a mwArray. */
  TensorMultiplier getEffectiveOperatorMultiplierSingleSite(const MPS& ket,const MPO& ops,int pos);

  /** As the previous method, but gets the effective operator for k
      adjacent sites starting on position pos */
  TensorMultiplier getEffectiveOperatorMultiplierMultiSite(const MPS& ket,const MPO& ops,int k,int pos);

  /** Variational miniization of the variance of the given MPO,
      with a non-quadratic cost function. */
  void minimizeVariance(MPS& init,const int D,const MPO& ops,const double E0,
			double penH,double* var,
			const string& tmpfile="",int freqSv=200,int stopRound=1E6);
  
 private:
  Contractor();
  Contractor(const Contractor& c); // I only want a unique copy of this

  /** Auxiliary functions for the contractions: calculate left and 
     right temporary terms that are stored in TmpManager.
     bool gauge set to 1 indicates that the vectors are assumed to satisfy
     gauge condition (so that norm terms do not need to be calculated)
     bool ketN set to 1 indicates that the norm term is to be calculated
     with the ket MPS (as needed for optimizeL). This has no effect if
     gauge is not set. */
  void calculateL(int pos,TmpManager& theMgr,const MPO& ops,
		  const MPS& ket,const MPS& bra,bool gauge=0,bool ketN=0);
  void calculateR(int pos,TmpManager& theMgr,const MPO& ops,
		  const MPS& ket,const MPS& bra,bool gauge=0,bool ketN=0);
  /** Same, when no operator in between */
  void calculateL0(int pos,TmpManager& theMgr,
		  const MPS& ket,const MPS& bra,bool gauge=0);
  void calculateR0(int pos,TmpManager& theMgr,
		  const MPS& ket,const MPS& bra,bool gauge=0);
  /** Same for inverse optimization */
  void calculateL2(int pos,TmpManager& theMgr,const MPO& ops,
		   const MPS& ket,const MPS& bra,bool resolvent=0);
  void calculateR2(int pos,TmpManager& theMgr,const MPO& ops,
		  const MPS& ket,const MPS& bra,bool resolvent=0);

  /** Similar to calculateL(R), but with a trick to keep the partial
      contractions of the computed state with the already computed
      eigenvectors, in order to implement the search for the next
      excited state. */
  void calculateLexck(int pos,TmpManagerExc& theMgr,const MPO& ops,
		      const MPS& ket,const std::vector<MPS*>& braLevels);
  void calculateRexck(int pos,TmpManagerExc& theMgr,const MPO& ops,
		      const MPS& ket,const std::vector<MPS*>& braLevels);

  /** Similar to calculateL(R), but keep in operL(R) terms the
      effective operator for opsA, and in normL(R) the one for opsB,
      so that it can be used to solve a GEVP. */
  void calculateLgev(int pos,TmpManager& theMgr,const MPO& opsA,const MPO& opsB,
		     const MPS& ket,const MPS& bra);
  void calculateRgev(int pos,TmpManager& theMgr,const MPO& opsA,const MPO& opsB,
		     const MPS& ket,const MPS& bra);

  /** Similar to calculateL(R), but with a trick to keep the partial
      contractions of the computed state with the already computed
      eigenvectors, in order to implement the search for the next
      excited state, and in the case of GEVP, so that a second
      operator is used. */
  void calculateLgevexck(int pos,TmpManagerExc& theMgr,const MPO& opsA,const MPO& opsB,
			 const MPS& ket,const std::vector<MPS*>& braLevels);
  void calculateRgevexck(int pos,TmpManagerExc& theMgr,const MPO& opsA,const MPO& opsB,
			 const MPS& ket,const std::vector<MPS*>& braLevels);

  /** Perform one step in the contraction of <bra|Ops|ket> from the left, 
     given the temporary result and the tensors corresponding to one site. */
  void contractOperL(mwArray& result,const mwArray& tmp,
		     const Site& ket,const Operator& op,
		     const Site& bra);
  /** Perform one step in the contraction of <bra|Ops|ket> from the right, 
     given the temporary result and the tensors corresponding to one site. */
  void contractOperR(mwArray& result,const mwArray& tmp,
		     const Site& ket,const Operator& op,
		     const Site& bra);

  // WARNING: THis calculates norm: uses only one vector!
  void contractNormL(mwArray& result,const mwArray& norm,const Site& ket);
  void contractNormR(mwArray& result,const mwArray& norm,const Site& ket);

  void calculateMket(mwArray& result,const mwArray& tmpL,const Site& ket,
		     const Operator& op,const mwArray& tmpR);
  void calculateMbra(mwArray& result,const mwArray& tmpL,const Site& bra,
		  const Operator& op,const mwArray& tmpR);
  /** Contract letf (tmpL) and right (tmpR) terms relative to 
      one site to obtain the N matrix
      tmpL(tmpR) is the temporary result, from terms to the left(right) 
      (dimensions [beta,beta]) */
  void calculateN(mwArray& result,const mwArray& tmpL,
		  const mwArray& tmpR,int d);

  void contractL(mwArray& result,const mwArray& tmpL,const Site& ket,
		 const Site& bra);
  void contractR(mwArray& result,const mwArray& tmpR,const Site& ket,
		 const Site& bra);

  /** For a PBC case, or for computing the negativity, I may need to
      contract transfer matrices in the middle of the chain. This
      method takes a temporary calculation tmpL, with indices
      alpha (which could be betab_ x betak_),betab,betak (in this
      order) and contracts the next site (bra tensor has left virtual 
      index betab) returning thre result in the first array.*/
  void contractLmiddle(mwArray& result,const mwArray& tmpL,const Site& ket,
  		       const Site& bra);

  /** Contract the term from left, the term from right, and the ket,
      assuming operator is the identity and ket is an MPS site (from
      the ket), dimension [d,Dl,Dr], tmpL is the temporary result,
      from terms to the left (dimensions [betab,betak]) and tmpR the
      analogous term coming from right.  It returns a (column) vector
      of dimension d*Dl*Dr in result */
  void contractMketId(mwArray& result,const mwArray& tmpL,
		      const Site& ket,const mwArray& tmpR);
  /** The same, but contracting the bra */
  void contractMbraId(mwArray& result,const mwArray& tmpL,
		      const Site& bra,const mwArray& tmpR);


  /** Compute temporary terms to right and left for the inverse optimization */
  void contractL2(mwArray& result,const mwArray& tmpL,const Site& ket,const Site& bra);
  void contractR2(mwArray& result,const mwArray& tmpR,const Site& ket,const Site& bra);

  /** Solve the matrix equation that results from the optimization of
      a particular A tensor in optimize, N A= M. result returns the
      solution A, and contraction1=A+NA and contraction2=A+M  */
  void solveSite(const mwArray& N,const mwArray& M,mwArray& result,
		 complex_t& contraction1,
		 complex_t& contraction2,double svdtol=0);
  /** If gauge conditions are imposed at each step, the N matrix in the 
     equation is the identity, and the solution is just the rhs vector.*/
  void solveSiteGauge(mwArray& result,mwArray& contraction1,mwArray& contraction2,const mwArray& M,double N);

  bool checkDistanceConvergence(double& distance,const mwArray& contA,
				const mwArray& contAB,double normA,
				double normB,int pos);

  /** Common logic of finding Left and Right Max Eigenvectors, but starting
     with a given state.*/
  void findMaxEigenvector(const MPO& ops,int D,double* lambda,
			  char dir,MPS& initialMPS);

  /** Auxiliary to decide distance between MPS and sum of vectors */
  double distanceVectors(const std::vector<const MPS*>& kets,
			 const std::vector<complex_t>& beta,
			 const MPS& init);

  /** Private eigs method that will dispatch to the eigensolver specified by this->solver */
  void solveEigs(Multiplier& multi,const int k,const char* which,
		 std::vector<complex_t>& D,mwArray& U,bool eigv=0,double tol=0.);

  /** Idem, for eigenvalue closest to a particular target */ 
  void solveEigsClosest(Multiplier& multi,const int k,const char* which,double targetVal,
			std::vector<complex_t>& D,mwArray& U,double tol=0.);

  /** Attack a generalized eigenvalue problem */
  void solveEigsGen(Multiplier& multiA,Multiplier& multiB,const int k,const char* which,
		 std::vector<complex_t>& D,mwArray& U,bool eigv=0,double tol=0.);

  /** Auxiliary for findGSTwoSites and EffectiveOperatorMultiplier...*/
  void blockOperatorMultiSite(mwArray& blockOp,const MPO& ops,int k,int pos);

  /** Run a gradient descent to find the min variance of an operator, at
      fixed energy (and norm). */
  double gradientMinVariance(mwArray& A,
			     TensorMultiplier* multiH2,
			     TensorMultiplier* multiH,
			     double E0,double penH);

 public:

  /** Test functionality */
#ifdef TESTINGCNTR
  void test(const MPO& ops,const MPS& psi1,
	    const MPS& psi2);
#endif

};





#endif // CONTRACTOR_H
