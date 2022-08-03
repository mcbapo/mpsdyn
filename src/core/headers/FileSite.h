/**
   \file FileSite.h
   Definition of the class that stores site tensors in disk
   
   \author Mari-Carmen Banuls
   \date 31/01/2012
*/

#ifndef FILESITE_H
#define FILESITE_H

#include "Site.h"

using namespace std;


/** 
    Class that contains the A tensor for a particular site in the
    chain, but stored always in the disk.  The mwArray will be read
    when the data are needed and written when modified.  The only
    difference is that a directory can be given as constructor argument, to
    know where top store the file. If not, a default directory is used.
*/

class FileSite:
public Site {

  // A common directory to store all Sites
  static char* dir;
  static bool dirset;

  char* filename; // My particular file
  bool inmem; // control 
  friend class MPS;
  shrt::Indices dims; // I keep the dimensions, for fast access

 public:
  /** Create a FileSite (default is state |0>) */
  FileSite(int pos,int d,int Dl,int Dr,const char* dir=0);
  /** Create a FileSite passing an array with dimensions */
  FileSite(int pos,int* dimensions,const char* dir=0);
  /** Copy constructor. Allow a different position to be set */
  FileSite(const FileSite& s,int pos=0);
  /** 
      Create a FileSite for the given position, with a copy of this array as
      data. The dimensions of the array are taken as dxDlxDr */
  FileSite(int pos,const mwArray& A);
  ~FileSite();

// Return values
  int getd() const{return dims[0];}
  int getDl() const{return dims[1];}
  int getDr() const{return dims[2];}
  shrt::Indices getDimensions() const {return dims;}

  static void setDir(const char* dir);

  const mwArray getA() const;

  // Copy this data.
  void setSite(const mwArray& newA);
  // Copy this data but permuting its dimensions as newdims indicates
  void setRotatedSite(const mwArray& newA,const shrt::Indices& newdims,
		      bool conjugate=0);
  // Set all tensor elements to zero
  void setToZero();

  /** Debugging utilities */
  //friend ostream& operator<<(ostream& os, const FileSite& site);

  /** Set the array so that it represents the state of one particle given 
     as argument. */
  void setState(SingleState state);
  /** Overload assignment */
  FileSite& operator=(const FileSite& s);

  /** Read from a binary file */
  void load(ifstream& data);
  /** Save to a binary file */
  void save(ofstream& data) const;
  /** Read from a text file */
  void loadtext(ifstream& data);
  /** Save to a text file */
  void savetext(ofstream& data) const;

#ifndef TESTINGSITE
 protected:
#endif
  /** Auxiliary method for output*/
  void put(ostream& os) const;

  /** Apply the gauge condition to the right. This should only be done 
      by MPS, as it has no sense on an isolated tensor. 
      Argument rightterm is a reference that upon return contains the
      tensor to be multiplied on the site to right, and leftterm is the 
      same term received from the site to the left. */
  void gaugeR(mwArray& rightterm,const mwArray& leftterm,int cutD=0);

  // The same functionality for gauge condition to the left
  void gaugeL(mwArray& leftterm,const mwArray& rightterm,int cutD=0);

  static void gaugeR(FileSite* ket,FileSite* bra,mwArray& righttermK,
		     mwArray& righttermB,const mwArray& lefttermK,
		     const mwArray& lefttermB);

  // Apply the corresponding dimension reduction if the gauge condition 
  // requires it.
  void reduceDimR(const mwArray& leftterm);
  void reduceDimL(const mwArray& rightterm);

  /** Apply gauge condition to the left, at the same time an operator
      (or its adjoint, if dagger==true) is applied and the resulting
      bond dimension is truncated according to cutD.  If cutD is
      received, it is used to cut off the SVD, so that the resulting
      Site will have Dl<=D.
  */
  void gaugeLOp(mwArray& leftterm,const mwArray& op,
		const mwArray& rightterm,bool dagger=false,
		int cutD=0);
  // Associated contraction of a remaining term
  void reduceDimLOp(const mwArray& op,const mwArray& rightterm,
		    bool dagger=false);

  // To keep track of a pi phase when converting from an extern set of 
  // tensors that is then subject to gauge condition. This is only useful 
  // for test programs.
  void changeSign();

  /** Allow increasing/decreasing the bond dimension, but only from MPS 
      (so that it is consistent along the whole state) */
  void increaseBondDim(int Dl,int Dr);
  void decreaseBondDim(int Dl,int Dr);
  /** Idem for the physical dimension */
  void increasePhysDim(int d);
  /** Cut the higher indices of d, in order to adjust the Site to a given 
     physical dimension. This should be used ONLY before any optimization 
     routine, in order to produce an initial MPS of the proper dimensions. */
  void decreasePhysDim(int d);

  /** The static method gaugeR should be able to change the inmem status 
      for operating and then leave it as initially. */
  bool isInMem(){return inmem;}

  /** Just to make it easier to read and write from file. read
      recovers the mwArray stored in the file, and puts it in the
      local variable, write replaces the one in the file and
      discard clears the stored one. 
  */
  void read();
  void write();

 private:
  void discard();
  /** Initialize the filename for each constructor */
  void initFile(const char* dir=0);
};


#endif //SITE_H
