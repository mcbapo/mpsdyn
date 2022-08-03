
#include "TensorMultiplierOffset.h"
#include "Indices.h"

using namespace std;
using namespace shrt;

TensorMultiplierOffset::TensorMultiplierOffset(const mwArray& AL_,
					       const mwArray& Aop_,
					       const mwArray& AR_,
					       const double offset_):
  TensorMultiplier(AL_,Aop_,AR_),offset(offset_){};
  
void TensorMultiplierOffset::product(const mwArray& input,mwArray& result){
  TensorMultiplier::product(input,result);
  result=result+offset*input;
}

TensorMultiplierOffset::~TensorMultiplierOffset(){};

void TensorMultiplierOffset::getFullTensor(mwArray& result){
  //  std::cout<<"In TensorMultiplierOffset::getFullTensor"<<std::endl;
  TensorMultiplier::getFullTensor(result);
  result=result+offset*identityMatrix(getSize());
}
