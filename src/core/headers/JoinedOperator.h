#ifndef JOPERATOR_H
#define JOPERATOR_H

#include "Operator.h"
#include <vector>
#include <iostream>

/** 
    Implements a specific operator coming from a sequence of
    operators, so that it contains several operators, and knows how to
    contract them efficiently with proper left and right terms.  */


class JoinedOperator: public Operator {

 protected:
  int n; // number of operators inside
  //Jtype me;
  std::vector<mwArray> data; // A dynamic array with all data (copies!)

 public:
  /** 
      From a list of operators (passed as an array of pointers)
      compose the JoinedOperator. The first operator appearing in the
      list is the last to act on the ket (top one). In this first
      implementation, I will always keep a copy of the arrays
      involved.  d,D are the physical and bond dimensions of the
      original operators.
      \warning If called on composite Operator, it will save copies of
      the rearranged data tensors (losing one level of structure).
      \todo Write the ostream<< 
      \todo Implement keeping pointers, as Operator!!
  */
  JoinedOperator(int n,const Operator* ops[]);
  /** 
      A safer option: pass the pointers in a vector
      \todo Keep the pointers, not the copy!
  */
  JoinedOperator(int n,const std::vector<Operator*>& ops);

  /** Destructor */
  ~JoinedOperator();

  /** Copy constructor */
  JoinedOperator(const JoinedOperator& op);
  JoinedOperator& operator=(const JoinedOperator& op);

  // When I just want to store the pointers
  void setData(int nr,const mwArray* opers[]);
  // Safe version of the same
  void setData(int nr,const std::vector<mwArray*>& opers);

  /** Particular implementations of contraction methods */
  void contractL(mwArray& result,const mwArray& termL,const Site& ket,
		 const Site& bra,bool dagger=0) const;
  void contractR(mwArray& result,const mwArray& termR,const Site& ket,
		 const Site& bra,bool dagger=0) const;

  void contractMket(mwArray& result,const mwArray& termL,const Site& ket,
		    const mwArray& termR,bool dagger=0) const;
  void contractMbra(mwArray& result,const mwArray& termL,const Site& bra,
		    const mwArray& termR) const;
  void contractN(mwArray& result,const mwArray& termL,const mwArray& termR,
		 bool dagger=0) const ;
  /* Get a copy of the expanded data! */
  const mwArray getFullData() const;
  /** Get a copy of the idx-th data array */
  const mwArray getDataN(int idx) const{return data[idx];}

  /** Return the number of components */
  int getNumber() const{return n;}

  /** Save and retrieve from text files */
  void savetext(std::ofstream& outfile) const;
  void loadtext(std::ifstream& outfile);

  void conjugateOp();
  void permuteOp(shrt::Indices perm);
  // TODO: Overload the basic one!!
  void contractRight(const JoinedOperator* rightOp);

  void put(std::ostream& os) const;
  //friend std::ostream& operator<<(std::ostream& os,const JoinedOperator& oper);
 private:

  // Private contractions
  void contractleftjop(mwArray& result,const mwArray& termL,const Site& ket,
		       const Site& bra,bool dagger) const;
  void contractrightjop(mwArray& result,const mwArray& termR,const Site& ket,
		       const Site& bra,bool dagger) const;
  void contractMketjop(mwArray& result,const mwArray& termL,const Site& ket,
		       const mwArray& termR,bool dagger) const;
  void contractMbrajop(mwArray& result,const mwArray& termL,const Site& bra,
		       const mwArray& termR,bool dagger) const;
  void contractNjop(mwArray& result,const mwArray& termL,
		    const mwArray& termR,bool dagger) const;
 protected:

  friend class MPO;
  friend class TimeEvolutionRow;
  friend class TimeDependentEvolutionRow;
  friend class TimeEvolutionFactory;


};

#endif // JOPERATOR_H
