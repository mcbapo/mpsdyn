#include "misc.h"
#include <fstream>

using namespace std;

bool file_exists(const string filename){
  ifstream ifile(filename.data());
  return !ifile.fail();
  //  return ifile.goodbit;
}

bool file_exists(const char* filename){
  ifstream ifile(filename);
  return !ifile.fail();
  //  return ifile.goodbit;
}

int gcd(int a,int b){
  //  cout<<"Computing gcd("<<a<<","<<b<<")";
  int p1=max(a,b);
  int p2=min(a,b);
  while(p2>0){
    int aux=(p1%p2);
    p1=max(aux,p2);
    p2=min(aux,p2);
  }
  //cout<<p1<<endl;
  if(p2==0) return p1;
  else{ exit(1);
    //   std::cout<<"What is going on here??"<<endl;exit(1);
  }
}

int lcm(int a,int b){
  return (a*b)/gcd(a,b);
}
