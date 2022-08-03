#include "TransposerLong.h"

using namespace std;
using namespace shrt;

TransposerLong::TransposerLong(const Indices& dimensions,const Indices& perm)
  :dims(dimensions),perm(permutation){
  cout<<"Creating TransposerLong for permutation "<<perm<<" on "<<dims<<endl;
  int nrcomponents=1;
  for(int k=0;k<dims.size();k++) nrcomponents*=dims[k];
  int moved=0;
  for(int initk=0;moved<nrcomponents;initk++){
    // compute the cycle starting in initk
    //cout<<"Starting cycle in "<<initk<<endl;
    vector<int>* cyclek=new vector<int>;
    // current indices are i1=mod(k,d1), i2=k/d1 (integer div.)
    // new position is i2+d2*i1 = k/d1+d2*(k%d1)
    cyclek->push_back(2*initk);
    int next=initk;bool seen=false;
    do{
      next=d2*(next%d1)+next/d1;
      if(next<initk){
	seen=true;
	//cout<<"Cycle contains "<<next<<" already seen"<<endl;
      }
      else{
	cyclek->push_back(2*next);
	//moved++;
      }
    }
    while(!seen&&next!=initk);
    //    if(!seen&&cyclek->size()>2){
    if(!seen){
      cycles.push_back(cyclek);moved+=cyclek->size()-1;
      //cout<<"Added cycle "<<*cyclek<<", moved="<<moved<<endl;
    }
    else{
      delete cyclek;
    }
  }
  nrcycles=cycles.size();
}

ostream& operator<<(ostream& os,const TransposerLong& transp){
  os<<"TransposerLong "<<transp.dims<<" contains "<<transp.nrcycles<<" cycles:"<<endl;
  for(int k=0;k<transp.nrcycles;k++){
    os<<*transp.cycles[k]<<endl;
  }
  return os;
}

TransposerLong::~TransposerLong(){
  while(cycles.size()>0){
    delete cycles.back();
    cycles.pop_back();
  }
}

