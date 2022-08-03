#ifndef TMPMGR_H
#define TMPMGR_H

#include "mwArray.h"
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>

#include <vector>

using namespace std;
#ifndef DISKTMPMGR  
typedef vector<mwArray> vector_t ;
#else
#include "FileVector.h"
typedef  FileVector<mwArray> vector_t ;
#endif // DISKTMPMGR  

/** 
    Class that keeps temporary partial calculations, for efficiency.
    This is very specific for MPS (1D) optimization, with an operator
    row and a norm as the only things to be calculated (plus some
    optional terms for related calculations as resolvent or
    inverse). For other applications, it might have different
    aspect. It is not implemented as the most general base class...
    IMPORTANT: This class does not calculate anything, it merely
    stores results (and dispose them when outdated) \warning {WARNING:
    All indices run from 0 to len-1 This should be private to
    Contractor!!}  */


class TmpManager{
  
  //const mwArray empty;  

 protected:
  int len;
  vector_t operL;
  vector_t normL;
  vector_t operRL;
  vector_t operEL;
  vector_t operR;
  vector_t normR;
  vector_t operRR;
  vector_t operER;

  // Case of exponential and offset, need linear terms on both sides
  vector_t operLoff;
  vector_t operRoff;

  bool* uptodateL;
  bool* uptodateR;

  bool gauge; // whether gauge condition is granted (so that no norm terms are needed)
  bool resolvent; // whether operRL and operRR are needed
  bool exponential; // whether operEL and operER are kept
  bool exponentialOff; // whether operLoff and operRoff are kept

  friend class Contractor;

 public:
/** 
    If gauge=1, the norm term is not kept because the MPS is assumed
    to fulfill the gauge condition. However, this term is used in
    other optimizations (resolvent or exponential) to save more
    complex terms. In those cases, TmpManager has to be contructed
    with gauge=0; 
    If resolvent=1, extra terms are mantained for <Psi|Phi>. This is
    also needed in the exponential.
    If exponential=1, yet one extra term is mantained.
    \todo Rewrite this so that it is initialized with a label that
    specifies which optimization is run.
 */
  TmpManager(int len,bool gauge=0,bool resolvent=0,bool exponential=0,
	     bool exponentialOff=0);
  virtual ~TmpManager();

  bool isuptodateL(int site) const;
  bool isuptodateR(int site) const;

  void outofdateL(int site);
  void outofdateR(int site);

  /* ******************************** 
  // Fills the pointers with the location of stored data
  bool retrieveL(int site,mwArray** oper,mwArray** norm=0);
  bool retrieveR(int site,mwArray** oper,mwArray** norm=0);


  // To be used only if gauge=0
  void storeL(int site,mwArray& oper,mwArray& norm);
  void storeR(int site,mwArray& oper,mwArray& norm);
  
  ******************************** */

 protected:
  // Contractor will register matrices as updated!
  void storedL(int site);
  void storedR(int site);

  // For Contractor to retrieve data!
/*   mwArray& getOperL(int pos); */
/*   mwArray& getOperR(int pos); */
/*   mwArray& getNormL(int pos); */
/*   mwArray& getNormR(int pos); */
 
 // WARNING! I am deleting whatever was still stored!
  TmpManager();
  void clear();

 private:
  TmpManager(const TmpManager& tmp); // no copy wanted

// Assignment also forbidden
  TmpManager& operator=(const TmpManager& tmp);

 public:
#ifdef TESTFUN
  friend ostream& operator<<(ostream& os,const TmpManager& tmp);
#endif

};


#endif // TMPMGR_H
