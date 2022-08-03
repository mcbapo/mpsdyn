#ifndef PXP_H
#define PXP_H

#include "Hamiltonian.h"
#include <vector>

/** 
  Implementation of a PXP Hamiltonian,
  \f[
  H=-J (\sum_{i=1}^{L-2} P^{[i]} \sigma_x^{[i+1]} P^{[i+2]} + \sigma_x^{[1]} P^{[2]} + P^{[L-1]} \sigma_x^{[L]} ), 
  \f]
  where \f$P=(1+\sigma_z)/2\f$ is the projector onto \f$|0\rangle\f$.
  In the case of periodic boundary conditions, instead of the above, we have
  \f[
  H=-\sum_{i=1}^{L} P^{[i]} \sigma_x^{[i+1]} P^{[i+2]}.
  \f]
  The class also implements a projector onto the constrained
  subspace (annihilating states with two excitations next to each
  other).
  If a penalty lambda is provided, the Hamiltonian includes a
  penalty term associating that value to any pair of neighbouring
  excitations. This can be used instead of the projector for some
  purposes.

*/

class PXPHamiltonian:
    public Hamiltonian{
        protected:
            int L; // chain length
            int d; // dimension (only 2 supported)
            MPO hamil;
            mwArray Z;
            double J;
            double offset;
            double lambda; // penalty
            bool pbc;
            MPO proj;

        public:
            PXPHamiltonian(int L, double J=1., double lambda=0.,bool pbc=0,double offset=0.);

            ~PXPHamiltonian();

            const MPO& getHMPO() const {return hamil;}

            bool hasHMPO() const {return true;}

            /** Return a projector onto the constrained subspace (no two 1s together) */
            const MPO& getProjectorConstr() const {return proj;}
            
            /* Return e^(delta H)
             * Only for obc at present!
             * */
            void getUMPO(MPO &expH, complex_t delta) const;
            void getUMPO(MPO &expH, double delta) const {getUMPO(expH, ONE_c*delta);};


            /** Construct the MPO corresponding to the evolution operator with
              the splitting even-odd in the Hamiltonian. Take as argument the
              complex factor in the exponential, delta.  For real time
              evolution, delta is -i*tau, for imaginary time, it should be
              -tau.  NOT YET READY! In this case, three body terms, needs to
              be handled differently.*/
            virtual void getExponentialMPOeven(MPO& expHe,complex_t delta) const {std::cout<<"ERROR: Exponential for PXP HAmiltonian not ready!"<<endl;exit(1);};
            virtual void getExponentialMPOodd(MPO& expHo,complex_t delta) const  {std::cout<<"ERROR: Exponential for PXP HAmiltonian not ready!"<<endl;exit(1);};


            /** Construct an MPO for the commutator with an operator, and return in in the argument.  */
            void getCommutatorMPO(MPO& mpo);

        protected:

            /** Prepare the decomposition of the evolution operator:
              returns the two (three) terms to be inserted in the MPO for the 
              exp(delta H12) (exp(delta H123), respectively) (delta a complex number).
              TI is not assumed and the exponential of the term
              starting in pos is returned.  
              */
            void getThreeBodyTermExponential(mwArray &Ol, mwArray &Om, mwArray &Or, complex_t delta, int pos)  const;
            void getThreeBodyTermExponential(mwArray &Ol, mwArray &Om, mwArray &Or, const mwArray &H123, complex_t delta) const;
            void getTwoBodyTermExponential(mwArray &Ol, mwArray &Or, complex_t delta, int pos) const;
            void getTwoBodyTermExponential(mwArray &Ol, mwArray &Or, const mwArray &H12, complex_t delta) const;

            /** Compute the two(three)-body term on positions (\param<pos>,\param<pos>+1(, \param<pos>+2)):
            */
            void computeThreeBodyTerm(mwArray &result, int pos) const;
            void computeThreeBodyTerm(mwArray &result) const;
            void computeTwoBodyTerm(mwArray &result, int pos) const;
            void computeTwoBodyTerm(int lr, mwArray &result) const;

        private:
            void initZ();
            void initHMPO();
            void initHMPOpbc();
            void initProjMPO();
            void initProjMPOpbc();
            //void initHMPO(double* tis=0,double* Uis=0,double* muis=0,
            //	 double* Vis=0);
            /** Common computation for recovering the even and odd parts of the exp */
            void getExponentialMPO(MPO& expHo,complex_t delta,bool even) const;
            //void getDoubleExponentialMPO(MPO& expHo,complex_t delta,bool even) const;
            //void getExtendedExponentialMPO(MPO& expHo,complex_t delta,bool even) const;

            /* // This is copied from IsingHamiltonian, but it would be better to */
            /* // have it somewhere else (common single place) */
            /* /\**  */
            /*     Auxiliary function for the construction of the commutator */
            /*     MPO. If this should be more general (not only for this */
            /*     IsingHamiltonian case), these functions might need to be */
            /*     somewhere else. Specially the double */
            /*  *\/ */
            /* void initZdoub(mwArray& Z); */
            /* void constructOperatorProduct(mwArray& result,const mwArray& opA, */
            /* 				const mwArray& opB) const; */



    };

#endif // PXP_H
