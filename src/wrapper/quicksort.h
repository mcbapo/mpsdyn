#ifndef QUICKSORT_H
#define QUICKSORT_H

#include <cstdlib>
#include <cstdio>
#include <vector>

// Partition the vector orig, from position l to r (pivot in l) 
// and return the new position of the pivot (pos)
template <class T>
void partition(std::vector<T>& orig,int l,int r,int& pos){
  // pivot is the first
  T p=orig[l]; // pivot
  int i=l+1; // in the partitioned part, first larger than p 
  for(int j=l+1;j<=r;j++){
    if(orig[j]<=p){
      T aux=orig[j];
      orig[j]=orig[i];
      orig[i++]=aux;
    }
  }
  orig[l]=orig[i-1];
  orig[i-1]=p;
  pos=i-1;
}

// Choose random pivot
template <class T>
void choosepivotR(std::vector<T>& orig,int l,int r){
  bool static initialized=false;
  if(!initialized){
    srandom(time(NULL));
    initialized=true;
  }
  int pivot=random()%(r-l+1)+l;
  T aux=orig[pivot];
  orig[pivot]=orig[l];
  orig[l]=aux;
  return;
}

// Choose an index for the pivot in the vector, between l and r
// and place it in position l
template <class T>
void choosepivot(std::vector<T>& orig,int l,int r){
  choosepivotR(orig,l,r);
}

// Run quicksort on a vector and return it sorted in increasing order
template <class T>
void quicksort(std::vector<T>& orig,int l,int r){
  if(r-l<=0) return;
  choosepivot(orig,l,r);
  int pos=l;
  partition(orig,l,r,pos);
  quicksort(orig,l,pos-1);
  quicksort(orig,pos+1,r);
}

// Run quicksort, but return in count the required number of comparisons
template <class T>
void quicksort(std::vector<T>& orig,int l,int r,long int& count){
  if(r-l<=0) return;
  choosepivot(orig,l,r);
  int pos=l;
  partition(orig,l,r,pos); // r-l comparisons
  quicksort(orig,l,pos-1,count);
  quicksort(orig,pos+1,r,count);
  count+=r-l;
}



  // 1) always the first => nothing to do
template <class T>
void choosepivot1(std::vector<T>& orig,int l,int r){
	return;
}
  // 2) always the last : exchange it with the first
template <class T>
void choosepivot2(std::vector<T>& orig,int l,int r){
  // 2) always the last : exchange it with the first
  T aux=orig[r];
  orig[r]=orig[l];
  orig[l]=aux;
  return;
}

  // 3) median of three
template <class T>
void choosepivot3(std::vector<T>& orig,int l,int r){
  // 3) median of three
  //cout<<"choosepivot3 called with vector["<<l<<"-"<<r<<"]"<<endl;
  int len=r-l+1;
  int med=(len%2==0)?l+len/2-1:l+(len-1)/2;
  //cout<<"middle element "<<med<<endl;
  //cout<<"values "<<orig[l]<<","<<orig[med]<<","<<orig[r]<<endl;
  int pivot=l;
  if(orig[l]<orig[r]){
    if(orig[med]<orig[l])
      pivot=l;
    else if(orig[med]<orig[r])
      pivot=med;
    else
      pivot=r;
  }
  else{
    if(orig[med]<orig[r])
      pivot=r;
    else if(orig[med]<orig[l])
      pivot=med;
    else
      pivot=l;
  }
  //cout<<"chosen "<<orig[pivot]<<endl;
  T aux=orig[pivot];
  orig[pivot]=orig[l];
  orig[l]=aux;
  return;
}


#endif // QUICKSORT_H
