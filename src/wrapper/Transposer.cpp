#include "Transposer.h"

using namespace std;
using namespace shrt;

Transposer::Transposer(const Indices& dimensions,const Indices& perm)
  :dims(dimensions),count(0),nrcycles(0),
   d1(dimensions[0]),d2(dimensions[1]),cycles(){
  if((dimensions.size()!=2)||(perm.size()!=2)||(perm[0]!=2)||perm[1]!=1){
    cout<<"Right now Transposer only supports transposing matrices! "<<endl;
    exit(212);
  }
  int nrcomponents=d1*d2;
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

ostream& operator<<(ostream& os,const Transposer& transp){
  os<<"Transposer "<<transp.dims<<" contains "<<transp.nrcycles<<" cycles:"<<endl;
  for(int k=0;k<transp.nrcycles;k++){
    os<<*transp.cycles[k]<<endl;
  }
  return os;
}

Transposer::~Transposer(){
  //  cout<<"Destroying Transposer("<<dims<<") after "<<count<<" uses"<<endl;
  while(cycles.size()>0){
    delete cycles.back();
    cycles.pop_back();
  }
}

/**
#define CACHE 500

void Transposer::transpose(mwArray& mwA,bool conj){
  // Assume the vector of cycles has already been created
  vector<int>* cyclesk[CACHE]; int read=0;
  complex_t hold[CACHE];
  int initk[CACHE];
  int last[CACHE];
  int ptr[CACHE];
  bool done[CACHE];
  while(read<nrcycles){
    int loaded=0;
    while(loaded<CACHE&&read<nrcycles){
      vector<int>* cycle_=cycles[read++];
      if(cycle_->size()==2){ // just conjugate, if needed and cont
	if(conj){
	  int single=cycle_->operator[](0);
	  //*(complex_t*)&mwA.components[single]=
	  //conjugate(*(complex_t*)&mwA.components[single]);
	  mwA.components[single+1]=-mwA.components[single+1];
	}
	continue;
      }
      cyclesk[loaded]=cycle_;
      initk[loaded]=cycle_->operator[](0);
      last[loaded]=cycle_->size()-1;
      done[loaded]=false;
      hold[loaded]=*(complex_t*)&mwA.components[initk[loaded]];
      ptr[loaded++]=1;
    }
    //cout<<"Loaded "<<loaded<<", read "<<read<<" of "<<nrcycles<<endl;
    int finished=0;
    while(finished<loaded){
      for(int k=0;k<loaded;k++){
	// if this cycle not yet finished
	if(!done[k]){
	// exchange hold and the one in the pos indicated by cyclek
	  complex_t aux=*(complex_t*)&mwA.components[cyclesk[k]->operator[](ptr[k])];
	  *(complex_t*)&mwA.components[cyclesk[k]->operator[](ptr[k]++)]=
	    conj?conjugate(hold[k]):hold[k];
	  hold[k]=aux;
	  if(ptr[k]==last[k]+1){
	    done[k]=true;finished++;
	  }
	}
      }
    }
  }

}
*/


void Transposer::transpose(mwArray& mwA,bool conj){
  count++;
  for(int k=0;k<nrcycles;k++){
    vector<int>& cyclek=*cycles[k];
    int initk=cyclek[0];int last=cyclek.size()-1;
    complex_t hold=*(complex_t*)&mwA.components[initk];
    int ptr=1;
    do{
      // exchange hold and the one in the pos indicated by cyclek
      complex_t aux=*(complex_t*)&mwA.components[cyclek[ptr]];
      *(complex_t*)&mwA.components[cyclek[ptr++]]=conj?conjugate(hold):hold;
      hold=aux;
    }
    while(ptr<=last);
    //cout<<"Transposer used cyle["<<k<<"]="<<cyclek<<endl;
  }
}

