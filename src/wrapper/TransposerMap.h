#ifndef TRANSPMAP_H
#define TRANSPMAP_H

#include <map>
#include "Indices.h"
#include "Transposer.h"
#include <iostream>

typedef std::map<shrt::Indices,Transposer*,bool(*)(const shrt::Indices&,const shrt::Indices&)> map_t;

typedef map_t::iterator iter_t;

typedef bool (*comp_fn)(const shrt::Indices&,const shrt::Indices&);

class TransposerMap:
public map_t{
public:
  static TransposerMap& getTransposerMap(){
    //    bool(*)(const shrt::Indices&,const shrt::Indices&)=shrt::compare;
    static TransposerMap one(shrt::compare);
    return one;
  }

// I replace the normal archiving function by one which uses the
// Indices inside the Transposer itself as key.
// TODO: maybe return bool, s.t. false if failed ?
  void add(Transposer* transp){
    //std::cout<<"TransposerMap added a Transposer with key "<<*transp->getDims()
    //	     <<", "<<*transp<<std::endl;
    insert(std::pair<shrt::Indices,Transposer*>(shrt::Indices(*transp->getDims()),transp)); 
  }


  void clear(){
    while(!empty()){
      iter_t iter=begin();
      delete iter->second;
      erase(iter);
    }

  }
  // When the program finishes, everything is cleared
  ~TransposerMap(){
    std::cout<<"Destroying TransposerMap, with "<<size()<<" Transposers"<<std::endl;
clear();}

private:
 TransposerMap(comp_fn pointer):map_t(pointer){};

};

#endif // TRANSPMAP_H
