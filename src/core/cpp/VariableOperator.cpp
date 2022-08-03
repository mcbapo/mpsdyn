#include "math.h"
#include "VariableOperator.h"


using namespace shrt;
using namespace std;

VariableOperator::VariableOperator(int d,int Dl,int Dr,int d2):
  Operator(d,Dl,Dr,d2){};

VariableOperator::VariableOperator(const mwArray& data_):
  Operator(data_.getDimensions()){
  data=new mwArray(reshape(data_,dims));
  mydata=1;
}

VariableOperator::~VariableOperator(){
  clear();
}

void VariableOperator::setData(const mwArray* oper){
  clear();
  if(oper->getRank()!=4){
    cout<<"Error: Trying to setData in VariableOperator from mwArray pointer "
	<<" to "<<*oper<<endl;
    exit(2);
  }   
  dims=oper->getDimensions();
  data=new mwArray(*oper);
  mydata=1;
}

void VariableOperator::setElement(complex_t value,const Indices& idx){
  if(idx.size()!=dims.size()){
    cout<<"Error: Trying to setElement("<<idx<<") in VariableOperator "
	<<" which has dimensions "<<dims<<endl;
    exit(2);
  }
  mwArray* aux=(mwArray*) data;
  aux->setElement(value,idx);
}

void VariableOperator::put(std::ostream& os) const{
  os<<"VariableOperator, ";
  Operator::put(os);
}
