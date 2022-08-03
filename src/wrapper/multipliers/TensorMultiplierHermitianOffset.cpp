
#include "TensorMultiplierHermitianOffset.h"
#include "Indices.h"

using namespace std;
using namespace shrt;

TensorMultiplierHermitianOffset::TensorMultiplierHermitianOffset(const mwArray& AL_,
					       const mwArray& Aop_,
					       const mwArray& AR_,
					       const double offset_):
  TensorMultiplierHermitian(AL_,Aop_,AR_),offset(offset_){};
  
void TensorMultiplierHermitianOffset::product(const mwArray& input,mwArray& result){
  TensorMultiplierHermitian::product(input,result);
  result=result+offset*input;
}

TensorMultiplierHermitianOffset::~TensorMultiplierHermitianOffset(){};
