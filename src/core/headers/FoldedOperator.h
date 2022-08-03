#ifndef FOPERATOR_H
#define FOPERATOR_H

#include "Indices.h"
#include "Operator.h"

/** 
    Implements a specific operator coming from folding of a certain
    row, so that it contains two operators, and knows how to contract
    them efficiently with proper left and right terms.
*/



class FoldedOperator: public Operator {

  // mwArray* data_; already in Operator, corresponds to R part
 protected:
  const mwArray* dataL; // the second set
  bool mydataL;
  //shrt::Indices dims; // The effective dimensions
 
 public:
  /** Receive two arrays, corresponding to the left and right
      operators that are gruped together now. They will be stored like
      that, and then used efficiently in the contractions.  The
      mwArrays must come with indices sorted as usual, up, left, down,
      right. But the one corresponding to the left term will be
      flipped (left <-> right) as would correspond to a folded MPO.
      If \param<norev> is true, it creates a usual FoldedOperator, but
      without flipping the left and right dimensions of the left
      mwArray.
  */
   FoldedOperator(const mwArray& dataL,const mwArray& dataR,bool norev=false);

  /** Destructor */
  ~FoldedOperator();

  /** Copy constructor */
  FoldedOperator(const FoldedOperator& op);
  FoldedOperator& operator=(const FoldedOperator& op);

  void conjugateOp();
  void permuteOp(shrt::Indices perm);
  // TODO: Implement the basic one, with argument Operator, so that it
  // tries a dynamic cast and fails if the type does not match.
  void contractRight(const FoldedOperator* rightOp);

  // When I just want to store the pointers
  void setData(const mwArray* operL,const mwArray* operR);

  void contractL(mwArray& result,const mwArray& termL,const Site& ket,
		 const Site& bra,bool dagger=0) const;
  void contractR(mwArray& result,const mwArray& termR,const Site& ket,
		 const Site& bra,bool dagger=0) const;

  void contractMket(mwArray& result,const mwArray& termL,const Site& ket,
		    const mwArray& termR,bool dagger=0) const;

  void contractMbra(mwArray& result,const mwArray& termL,const Site& bra,
		    const mwArray& termR,bool dagger=0) const;

  void contractN(mwArray& result,const mwArray& termL,
			       const mwArray& termR,bool dagger=0) const;
  /** Returns the expanded operator resulting from contracting both parts 
      together. */
  const mwArray getFullData() const;
  /** Return copies of the left and right data*/
  const mwArray getDataLeft() const{ return *dataL;}
  const mwArray getDataRight() const{ return *data;}

  /** Function that actually prints the Operator, and gets overloaded
      by daughter classes. */
  void put(std::ostream& os) const;


 private:

  void contractleftfop(mwArray& result,const mwArray& termL,const Site& ket,
		       const Site& bra,bool dagger) const;
  void contractrightfop(mwArray& result,const mwArray& termR,const Site& ket,
		       const Site& bra,bool dagger) const;
  void contractMketfop(mwArray& result,const mwArray& termL,const Site& ket,
		       const mwArray& termR,bool dagger) const;
  void contractMbrafop(mwArray& result,const mwArray& termL,const Site& bra,
		       const mwArray& termR,bool dagger) const;
  void contractNfop(mwArray& result,const mwArray& termL,
  		    const mwArray& termR,bool dagger) const;

 protected:
  void clear();
  void updateDims();

  friend class MPO;
  friend class TimeEvolutionRow;
  friend class TimeDependentEvolutionRow;
  friend class TimeEvolutionFactory;
  //  friend ostream& operator<<(ostream& os,const FoldedOperator& sigma);
};

#endif // FOPERATOR_H
