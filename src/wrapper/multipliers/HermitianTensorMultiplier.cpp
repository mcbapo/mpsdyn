
#include "HermitianTensorMultiplier.h"
#include "Indices.h"

using namespace std;
using namespace shrt;

HermitianTensorMultiplier::HermitianTensorMultiplier(const mwArray& AL_,
						     const mwArray& Aop_,
						     const mwArray& AR_):
  TensorMultiplier(AL_,Aop_,AR_),
  adjT(permute(conjugate(AL_),Indices(3,2,1)),
       permute(conjugate(Aop_),Indices(3,2,1,4)),
       permute(conjugate(AR_),Indices(3,2,1))){};

void HermitianTensorMultiplier::product(const mwArray& input,mwArray& result){
  mwArray tmp;
  TensorMultiplier::product(input,result);
  adjT.product(input,tmp);
  result=.5*(result+tmp);
  //cout<<"Using HermitianTM"<<endl;
}

HermitianTensorMultiplier::~HermitianTensorMultiplier(){
}

HermitianTensorMultiplier::HermitianTensorMultiplier(const TensorMultiplier& orig){
  AL_=orig.get
}:
  TensorMultiplier(orig.getLeftTerm(),orig.getMiddleTerm(),orig.getRightTerm()),  adjT(permute(conjugate(orig.getLeftTerm()),Indices(3,2,1)),
       permute(conjugate(orig.getMiddleTerm()),Indices(3,2,1,4)),
       permute(conjugate(orig.getRightTerm()),Indices(3,2,1))){};
