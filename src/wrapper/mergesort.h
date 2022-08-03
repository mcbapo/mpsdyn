#ifndef MERGESORT_H
#define MERGESORT_H

#include <cstdlib>
#include <cstdio>
#include <vector>

template <class T>
void split(const std::vector<T>& orig,const std::vector<int>& indices,
	   std::vector<int>& left,
	   std::vector<int>& right){
  int len=indices.size();
  int lenL=indices.size()%2==0?indices.size()/2:(indices.size()-1)/2;
  left=vector<int>(lenL,0);
  right=vector<int>(len-lenL,0);
  for(int k=0;k<len;k++)
    if(k<lenL) left[k]=indices[k];
    else right[k-lenL]=indices[k];
}

template <class T>
void merge(const std::vector<T>& orig,const std::vector<int>& left,
	   const std::vector<int>& right,
	   std::vector<int>& result,bool descend=true){
  int lenL=left.size();
  int lenR=right.size();
  result=vector<int>();
  bool doneL=false;
  bool doneR=false;
  int cntL=0;int cntR=0;
  while(!doneL||!doneR){
    if(doneL){
      result.push_back(right[cntR++]);
    }
    else if(doneR){
      result.push_back(left[cntL++]);
    }
    else{
      if(!descend){ //ascending order
	if(orig[left[cntL]]<=orig[right[cntR]])
	  result.push_back(left[cntL++]);
	else
	  result.push_back(right[cntR++]);
      }
      else{ //descending order
	if(orig[left[cntL]]>=orig[right[cntR]])
	  result.push_back(left[cntL++]);
	else
	  result.push_back(right[cntR++]);
      }
    }
    if(cntR>=lenR) doneR=true;
    if(cntL>=lenL) doneL=true;
  }
}

// Run mergesort on a vector containing elements of class T.
/** Because the original vector is constant, this method sorts the
    indices, so that the second argument, upon return, contains the
    right indices, sorted according to the criterion (if
    descend==true, which is the deault, it will be in descending
    order)
 */
template <class T>
void mergesort(const std::vector<T>& orig,std::vector<int>& indices,
	       bool descend=true){
  /* cout<<"mergesort received to sort "; */
  /* if(descend) cout<<"in descending order "; */
  /* else cout<<"in ascending order "; */
  /* for(int k=0;k<indices.size();k++) */
  /*   cout<<orig[indices[k]]<<","; */
  /* cout<<endl; */
  if(indices.size()==1) return;
  vector<int> left;  
  vector<int> right;
  split<T>(orig,indices,left,right);
  mergesort<T>(orig,left,descend);
  mergesort<T>(orig,right,descend);
  merge<T>(orig,left,right,indices,descend);
}



#endif // MERGESORT_H
