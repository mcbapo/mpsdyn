#include "DoubleOperator.h"

using namespace std;
using namespace shrt;

DoubleOperator::DoubleOperator(const mwArray& dataL_,const mwArray& dataR):
  FoldedOperator(dataL_,dataR,true){};

DoubleOperator::~DoubleOperator(){
  clear();
}

DoubleOperator::DoubleOperator(const DoubleOperator& op):
  FoldedOperator(op){};

DoubleOperator& DoubleOperator::operator=(const DoubleOperator& op){
  if(this!=&op){
    clear();
    this->FoldedOperator::operator=(op);
  }
  return *this;
}

void DoubleOperator::put(ostream& os) const{
  os<<"DoubleOperator: ";
  FoldedOperator::put(os);
}
