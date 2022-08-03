
#include "BasicMultiplier.h"

BasicMultiplier::BasicMultiplier(const mwArray& A_):A(A_){};

void BasicMultiplier::product(const mwArray& vector,mwArray& result){
  result=A*vector;
}
