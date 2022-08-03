#ifndef TMPMGREX_H
#define TMPMGREX_H

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
    Class that keeps temporary partial calculations need by \class
    <Contractor> in findNextExcitedState, for efficiency.
    IMPORTANT: This class does not calculate anything, it merely
    stores results (and dispose them when outdated) \warning {WARNING:
    All indices run from 0 to len-1 This should be private to
    Contractor!!}  */


class TmpManagerExc{
  
  //const mwArray empty;  

 protected:
  int len;
  int nextK;
  vector_t operL;
  vector_t operR;
  vector<vector_t > normL;
  vector<vector_t > normR;

  bool gev;
  vector_t normL0;
  vector_t normR0;

  bool* uptodateL;
  bool* uptodateR;

  friend class Contractor;

 public:

  TmpManagerExc(int len,int nextK,bool gev=0);
  virtual ~TmpManagerExc();

  bool isuptodateL(int site) const;
  bool isuptodateR(int site) const;

  int getNextK() const{return nextK;}

  void outofdateL(int site);
  void outofdateR(int site);

 protected:
  // Contractor will register matrices as updated!
  void storedL(int site);
  void storedR(int site);

 // WARNING! I am deleting whatever was still stored!
  TmpManagerExc();
  void clear();

 private:
  TmpManagerExc(const TmpManagerExc& tmp); // no copy wanted

  TmpManagerExc& operator=(const TmpManagerExc& tmp);

 public:
#ifdef TESTFUN
  friend ostream& operator<<(ostream& os,const TmpManagerExc& tmp);
#endif

};


#endif // TMPMGREX_H
