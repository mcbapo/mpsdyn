#include "wrapper.h"
#include "Indices.h"

using namespace std;

using namespace shrt;

typedef vector<int> list_t;

Indices::Indices(int d1){
  push_back(d1);
}
Indices::Indices(int d1,int d2){
  push_back(d1);  push_back(d2);
}
Indices::Indices(int d1,int d2,int d3){
  push_back(d1);  push_back(d2);
  push_back(d3);  
}
Indices::Indices(int d1,int d2,int d3,int d4){
  push_back(d1);  push_back(d2);
  push_back(d3);  push_back(d4);
}

Indices::Indices(int d1,int d2,int d3,int d4,int d5){
  push_back(d1);  push_back(d2);
  push_back(d3);  push_back(d4);
  push_back(d5);  
}

Indices::Indices(int d1,int d2,int d3,int d4,int d5,int d6){
  push_back(d1);  push_back(d2);
  push_back(d3);  push_back(d4);
  push_back(d5);  push_back(d6);
}

Indices::Indices(int d1,int d2,int d3,int d4,int d5,int d6,int d7){
  push_back(d1);  push_back(d2);
  push_back(d3);  push_back(d4);
  push_back(d5);  push_back(d6);
  push_back(d7); 
}

Indices::Indices(int d1,int d2,int d3,int d4,int d5,int d6,int d7,int d8){
  push_back(d1);  push_back(d2);
  push_back(d3);  push_back(d4);
  push_back(d5);  push_back(d6);
  push_back(d7);  push_back(d8);
}

Indices::~Indices(){
  list_t::clear();
}

void Indices::append(const vector<int>& orig){
  for(int k=0;k<orig.size();k++)
    push_back(orig[k]);
}

Indices::Indices(const Indices& ind1,const Indices& ind2):list_t(ind1){
  append(ind2);
}

Indices::Indices(const vector<int>& ind):list_t(ind){};

Indices::Indices(const Indices& ind):list_t(ind){};

Indices::Indices(const Indices& ind,int dlast):list_t(ind){
  push_back(dlast);
}

Indices& Indices::operator=(const Indices& ind){
  if(this!=&ind){
    int sz=size();
    if(ind.size()!=sz)
      list_t::operator=(ind);
    else{
      for(int k=0;k<sz;k++)
	list_t::operator[](k)=ind[k];
    }
  }
  return *this;
}


bool shrt::compare(const Indices& ind1,const Indices& ind2){
  //cout<<"Indices comparing "<<*this<<" < "<<ind2<<endl;
  if(ind1.size()<ind2.size()) return true;
  if(ind1.size()>ind2.size()) return false;
  for(int k=0;k<ind1.size();k++){
    if(ind1[k]<ind2[k]) return true;
    if(ind1[k]>ind2[k]) return false;
  }
  return false;
}
  

bool Indices::operator<(const Indices& ind2) const{
  //cout<<"Indices comparing "<<*this<<" < "<<ind2<<endl;
  if(size()!=ind2.size()) return false;
  for(int k=0;k<size();k++){
    if(this->operator[](k)>=ind2[k]) return false;
  }
  return true;
}

void Indices::increment(int val){
  for(int k=0;k<size();k++) this->operator[](k)+=val;
}  

bool Indices::operator==(const Indices& ind2) const {
  //cout<<"Indices comparing "<<*this<<" and "<<ind2<<endl;
  int sz=size();
  if(sz!=ind2.size()) return false;
  for(int k=0;k<sz;k++)
    if(this->operator[](k)!=ind2[k]) return false;
  return true;
}


bool Indices::operator!=(const Indices& ind2) const {
  return !(this->operator==(ind2));
}
