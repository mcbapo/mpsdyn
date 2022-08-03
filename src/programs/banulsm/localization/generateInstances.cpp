#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream> 

#include "mwArray.h"

using namespace shrt;
using namespace std;

/** Generate a number of random instances for the localization
    problem. In particular, it just generates lists of Lmax uniformly
    distributed real numbers in [-1,1] that will be used by
    localizationPTcont.cpp as the values of the magnetic field (set to
    h*random values)

    Receives as arguments: 
    \param <Lmax> (int) Maximum length of the chain of interest 
    \param <N> (int) nr of instances to be generated
    \param <filename> (char*) Name of the file where the numbers are
                      stored. It will be saved as a mwArray (complex!) 
		      in binary format, with dimensions N x Lmax 
*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int Lmax=atoi(argv[++cntr]);
  int N=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];

  mwArray H(Indices(N,Lmax));
  H.fillRandom();
  H=H+conjugate(H); // now real between 0 and 2
  mwArray aux(Indices(N,Lmax));aux.fillWithOne();
  H=H-aux;
  ofstream out(outfname);
  if(!out.is_open()){
    cout<<"Error: couldn't open file "<<outfname<<" for output"<<endl;
    exit(1);
  }
  H.save(out);
  out.close();
}
