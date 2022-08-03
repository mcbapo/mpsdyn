#include <math.h>
#include <iomanip>
#include "Properties.h"

using namespace std;

/** Test program, just to show how to use the Properties class.
    Properties can be given in a text file, as the one copied at the
    end of this file, and/or as command line options. In that case,
    the name of the property has to be preceded by a minus, followed
    by '=' and the value (without space), and the value will supersede
    the one in the file (if there was one).  In this particular
    example, if there is a configuration file including a list of
    properties, its name must be the first argument in the list.
    Properties names are strings (if containing white space, in the
    command line have to be enclosed in quotes).  This particular case
    has no default arguments, so the program will exit with error if
    not all of the expected properties are set.  */

int main(int argc,const char* argv[]){

  Properties props;
  // Initialize parameters to defaults
  if(argc>1){
    const char* fProp=argv[1];
    int cnt=1;
    if(fProp[0]!='-'){
      // A first argument (default properties)
      props.loadFromFile(argv[1]);
      if((argc>2))
	cout<<"Some properties may be replaced by command line arguments"
	    <<endl;
      cnt++;
    }
    props.loadProperties(argc-cnt,&argv[cnt]);
  }

  int wrongparameters=0;
  const string fileAOs=props.getProperty("fileAOs");
  if(fileAOs.empty()) wrongparameters++;
  double J=props.getDoubleProperty("J");
  if(J==-1) wrongparameters++;
  double mu=props.getDoubleProperty("mu");
  if(mu==-1) wrongparameters++;
  int M1=props.getIntProperty("timesteps");
  if(M1==-1) wrongparameters++;
  const string timefile=props.getProperty("timefile");
  if(timefile.empty()) wrongparameters++;
  const string outfname=props.getProperty("out");
  if(outfname.empty()) wrongparameters++;
  int app=props.getIntProperty("append");
  if(app==-1) wrongparameters++;
  int D=props.getIntProperty("D");
  if(D==-1) wrongparameters++;
  int K=props.getIntProperty("K");
  if(K==-1) wrongparameters++;
  int per=props.getIntProperty("per");
  if(per==-1) wrongparameters++;

  if(wrongparameters>0){
    cout<<"Wrong properties, or not enough values set "<<endl;
    cout<<endl;
    cout<<"Usage:"<<endl;
    cout<<endl;
    cout<<"     periodJ [props.conf] [-fileAOs=fname] [-J=double] "
	<<"[-mu=double] [-timesteps=int]"<<endl;	  
    cout<<"              [-timefile=fname] [-out=fname] [-append] [-D=int]"<<endl;	  
    cout<<"              [-K=int] [-per=int]"<<endl;
    exit(1);
  }

}

/**
 *  EXAMPLE OF A PROPERTIES FILE
 *  Names of the properties can be changed arbitrarily (string)
 *

# File with all properties required by the program
# individual values can be changed here or by command-line arguments
#
# Usage:
#
#	program [thisfile.conf] [-name1=value1] [-name2=value2] ...
#
# 
# file that keeps As, Os

fileAOs=input/datosPeriod.mat

# file with time data

timefile = times.dat

# chemical potential

mu=0.0

# J

J=.5

# number of time steps

timesteps	=	20

# bond dimension

D = 20

# output file

out=prueba.dat

# append to output file? (if not, file is recreated)

append=0

# range

K=2

# periodicity: one spin up every per sites (should be K+1)

per=3

*/
