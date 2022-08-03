
#include "LinearCombinationMultiplier.h"

using namespace std;
using namespace shrt;

LinearCombinationMultiplier::LinearCombinationMultiplier(const vector<complex_t>& coeffs_,			      
							 const vector<Multiplier*>& multis_,
							 complex_t offset_):
  coeffs(coeffs_),multis(multis_),offset(offset_){
  // check that sizes are ok!
  bool wrong=false;
  if(coeffs.size()!=multis.size()) wrong=true;
  else{
    N=multis[0]->getSize();
    for(int k=1;k<multis.size();k++){
      if(multis[k]->getSize()!=N){
	wrong=true; break;
      }
    }
  }
  if(wrong){
    cout<<"ERROR: Trying to create LinearCombinationMultiplier with incompatible arguments"<<endl;
    exit(1);
  }
}



LinearCombinationMultiplier::~LinearCombinationMultiplier(){
  coeffs.clear();
  multis.clear();
}

void LinearCombinationMultiplier::product(const mwArray& input,mwArray& result){
  mwArray aux(result.getDimensions()); // initialize to zeros
  for(int k=0;k<coeffs.size();k++){
    multis[k]->product(input,result);
    aux=aux+coeffs[k]*result;
  }    
  result=aux+offset*input;
}
