
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "FoldedOperator.h"
#include "DoubleOperator.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

using namespace shrt;

/** Read a thermal state for certain beta/2 from a MPS file (as
    produced by thermalStateIsing.cpp, for instance), and a certain
    local operator M from a text file, written in the basis of Pauli
    matrices, to then compute the expectation value of this M.  The
    MPS file should contain an MPO for rho(beta/2) and
    rho(beta)=rho(beta/2)^dagger rho(beta/2).
    
    Arguments:
    \param <mpsFile> (char*) file containing the MPS for rho(beta/2)
    \param <opFile> (char*) file containing the operator
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)
*/

// Construct the unitary that cnages basis from |i><j| to Pauli
mwArray basisChange();


int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* mpsfname=argv[++cntr];
  const char* operfname=argv[++cntr];
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  MPS rho2;
  rho2.importMPS(mpsfname);

  // This operator is in |i><j| basis, and the other one in the Pauli basis




}

mwArray basisChange(){
  complex_t* dataU={.5*ONE_c,ZERO_c_ZERO_c,.5*ONE_c, // first column
		    ZERO_c,.5*ONE_c,.5*ONE_c,ZERO_c,
		    ZERO_c,.5*I_c,-.5*I_c,ZERO_c,
		    .5*ONE_c,ZERO_c_ZERO_c,-.5*ONE_c} // last column
  mwArray U(Indices(4,4),dataU);

}
