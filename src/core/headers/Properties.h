#ifndef PROP_H
#define PROP_H

#include <cstdlib>
#include <map>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

/** 
    This class is designed to be initialized at the beginning of the
    program from one or more files, and to store the values of the
    different arguments/parameters passed to the program as a
    configuration file. Then the main program will recover the needed
    parameters (using the standard names also included in the config
    file) from the Properties table, instead of requiring dozens of
    inputs in the command line.
*/

class Properties{
  map<string,string> table_;

  // separator key-value ('=')
  static char keysep; 

 private:

  void parseString(string tmp);

 public:
  /** 
      Default constructor: no property loaded
  */
 Properties():table_(){};

  /** 
      Construction from a text file, where pairs of property-value
      are written.
  */
  Properties(const char* filename);

  ~Properties();

  /**
     Recover the property stored as name. If it does not exist, an empty
     string is returned
  */
  const string getProperty(const string& name) const;

  /** 
      Recover the (integer) value associated to name in the table of
      properties. If it does not exist, an exception should be thrown
      (TODO). Now, it returns -1.
  */
  int getIntProperty(const string& name) const;

  /** 
      Recover the (double) value associated to name in the table of
      properties. If it does not exist, an exception should be thrown
      (TODO). Now, it returns -1.
  */
  double getDoubleProperty(const string& name) const;

  /** 
      Recover the vector (double) value associated to name in the table of
      properties. If it does not exist, an exception should be thrown
      (TODO). Now, it returns empty vector.
  */
  std::vector<double> getVectorDoubleProperty(const string& name) const;

  /** 
      Recover the vector (int) value associated to name in the table of
      properties. If it does not exist, an exception should be thrown
      (TODO). Now, it returns empty vector.
  */
  std::vector<int> getVectorIntProperty(const string& name) const;

  /** 
      Add a new pair to the table
  */
  void addProperty(const string& name,const string& value);
  
  /** 
      Read all the pairs contained in a file
  */
  void loadFromFile(const char* filename);

  /** 
      Parse properties from the arguments received by the program, in
      the form:  -name=value
  */
  void loadProperties(int nr,const char* args[]);
};

#endif // PROP_H
