#include "TensorMultiplierProjMulti.h"

#include "Indices.h"

using namespace std;
using namespace shrt;

TensorMultiplierProjMulti::TensorMultiplierProjMulti(const TensorMultiplierProjMulti&):
  basic(NULL),vectors(),weights(),offset(ZERO_c),mine(0){};

TensorMultiplierProjMulti& TensorMultiplierProjMulti::operator=(const TensorMultiplierProjMulti&){return *this;}

TensorMultiplierProjMulti::TensorMultiplierProjMulti(const mwArray& AL_,
						     const mwArray& Aop_,
						     const mwArray& AR_,
						     vector<double> weights_,
						     const vector<mwArray>& vect_,			  
  						     bool hermitian,
						     complex_t offset_):
  vectors(vect_.size()),weights(weights_),offset(offset_),mine(true),
  basic(NULL){
#ifdef TESTFUN
  cout<<"Constructing a TensorMultiplierProjMultiple with "<<vect_.size()<<" terms"
      <<" from tensors "<<endl;
#endif
  basic=hermitian?new TensorMultiplierHermitian(AL_,Aop_,AR_):new TensorMultiplier(AL_,Aop_,AR_);
  for(int l=0;l<vect_.size();l++){
    vectors[l]=&vect_[l];
#ifdef TESTFUN
    cout<<" ["<<l<<"]="<<vect_[l];
#endif
  }
#ifdef TESTFUN
  cout<<endl;
#endif
  if(hermitian) offset.im=0.;
}

TensorMultiplierProjMulti::TensorMultiplierProjMulti(Multiplier& multi_,
						     vector<double> weights_,
						     const vector<mwArray>& vect_,			  
  						     complex_t offset_):
  vectors(vect_.size()),weights(weights_),offset(offset_),mine(false),
  basic(&multi_){
#ifdef TESTFUN
  cout<<"Constructing a TensorMultiplierProjMultiple with "<<vect_.size()<<" terms"
      <<" from tensors "<<endl;
#endif
  for(int l=0;l<vect_.size();l++){
    vectors[l]=&vect_[l];
#ifdef TESTFUN
    cout<<" ["<<l<<"]="<<vect_[l];
#endif
  }
#ifdef TESTFUN
  cout<<endl;
#endif
}

TensorMultiplierProjMulti::~TensorMultiplierProjMulti(){
  if(mine){ 
//     TensorMultiplierHermitian* left=dynamic_cast<TensorMultiplierHermitian*>(basic);
//     if(left!=0)
//       delete left;
//     else{
//       TensorMultiplier* left=dynamic_cast<TensorMultiplier*>(basic);
//       if(left!=0)
// 	delete left;
//       else{
// 	cout<<"ERROR: The Multiplier stored in TensorMultiplierProjMulti cannot be casted to Her or non Herm TM"<<endl;
// 	exit(1);
//       }
//     }
    delete basic;
  } 
  vectors.clear();
  weights.clear();
}


void TensorMultiplierProjMulti::product(const mwArray& input,mwArray& result){
  // The first contribution comes from the inner TensorMultiplier action
  basic->product(input,result);
  // and now we substract something (one term per vector)
  for(int k=0;k<vectors.size();k++){
    const mwArray& vect=*vectors[k];
    result=result-weights[k]*(vect*input)*Hconjugate(vect);
  }
  if(offset!=ZERO_c)
    result=result+offset*input;
}
  
void TensorMultiplierProjMulti::getFullTensor(mwArray& result){
  basic->getFullTensor(result);
  for(int k=0;k<vectors.size();k++){
    const mwArray& vect=*vectors[k];
    result=result-weights[k]*Hconjugate(vect)*vect;
  } 
  result=result+offset*identityMatrix(getSize());
  result=.5*(result+Hconjugate(result)); // should be done only if hermitian!
}
