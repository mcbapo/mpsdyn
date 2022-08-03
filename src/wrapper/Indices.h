#ifndef INDICES_H
#define INDICES_H

#include <vector>


namespace shrt{

/** Some of the most used vectors for dimensions and indices */

  class Indices:public std::vector<int> {
  public:
    // Create the most used lists of 1,2,3,4,6 or 8 elements
  Indices(int d1);
  Indices(int d1,int d2);
  Indices(int d1,int d2,int d3);
  Indices(int d1,int d2,int d3,int d4);
  Indices(int d1,int d2,int d3,int d4,int d5);
  Indices(int d1,int d2,int d3,int d4,int d5,int d6);
  Indices(int d1,int d2,int d3,int d4,int d5,int d6,int d7);
  Indices(int d1,int d2,int d3,int d4,int d5,int d6,int d7,int d8);
  Indices(const std::vector<int>& ind);
  Indices(const Indices& ind);
  Indices(){};
  ~Indices();
  // Concatenate two vectors
  Indices(const Indices& ind1,const Indices& ind2);
  // Concatenate just an integer
  Indices(const Indices& ind1,int dlast);
  Indices& operator=(const Indices& ind1);
  friend bool compare(const Indices& ind1,const Indices& ind2);
  bool operator<(const Indices& ind2) const;
  bool operator==(const Indices& ind2) const;
  bool operator!=(const Indices& ind2) const;
  //friend bool operator==(const Indices& ind1,const Indices& ind2);
  //friend bool operator!=(const Indices& ind1,const Indices& ind2);
// Increment each number by the argument
  void increment(int k=1);
  private:
  void append(const std::vector<int>& orig);
};

  bool compare(const Indices& ind1,const Indices& ind2);
  //bool operator==(const Indices& ind1,const Indices& ind2);
  //bool operator!=(const Indices& ind1,const Indices& ind2);

};

#endif // INDICES_H
