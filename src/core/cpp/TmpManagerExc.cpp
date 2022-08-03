#include "TmpManagerExc.h"


TmpManagerExc::TmpManagerExc():
  len(0),nextK(0),gev(0),
  operL(),normL(),normL0(),
  operR(),normR(),normR0(),
  uptodateL(NULL),uptodateR(NULL){};

TmpManagerExc::TmpManagerExc(const TmpManagerExc& tmp):
  len(0),nextK(0),gev(0),
  operL(),normL(),normL0(),
  operR(),normR(),normR0(),
  uptodateL(NULL),uptodateR(NULL){};

TmpManagerExc& TmpManagerExc::operator=(const TmpManagerExc& tmp){
  return *this;
}


TmpManagerExc::TmpManagerExc(int length,int nextK_,bool gev_):
  len(length),nextK(nextK_),gev(gev_),
  uptodateL(NULL),uptodateR(NULL),
  operL(length,mwArray::empty),operR(length,mwArray::empty),
  normL0(length,mwArray::empty),normR0(length,mwArray::empty),
  normL(nextK_,vector_t(length,mwArray::empty)),
  normR(nextK_,vector_t(length,mwArray::empty)) {
  //  cout<<"Creating a TmpManagerExc!"<<endl;
  uptodateL=new bool[length];
  uptodateR=new bool[length];
  for(int k=0;k<length;k++){
    uptodateL[k]=uptodateR[k]=0;
  }
  // initialize: left and rightmost sites are uptodate and equal to 1(+0 i)
  uptodateL[0]=1;
  uptodateR[length-1]=1;
  
  mwArray edge(ONE_c);
  operL[0]=edge;
  operR[length-1]=edge;

  if(gev){
    normL0[0]=edge;
    normR0[length-1]=edge;
  }
  else{
    normL0.clear();normR0.clear();
  }

  for(int k=0;k<nextK;k++){
    normL[k][0]=edge;
    normR[k][length-1]=edge;
  }

#ifdef TESTFUN
    cout<<"Test(just constructed TMpMgrExc): "<<*this<<endl;
#endif

}


TmpManagerExc::~TmpManagerExc(){
  clear();
  delete []uptodateL;
  delete []uptodateR;
  //  cout<<"Destructing TmpManager"<<endl;
}

void TmpManagerExc::clear(){
  // TODO: Delete things that might be kept uptodate till now
  operL.clear();
  operR.clear();
  for(int k=0;k<nextK;k++){
    normL[k].clear();
    normR[k].clear();
  }
  normL.clear();
  normR.clear();
  normL0.clear();
  normR0.clear();
}

bool TmpManagerExc::isuptodateL(int site) const{
  return uptodateL[site];
}

bool TmpManagerExc::isuptodateR(int site) const{
  if((site<0)||(site>=len)){
    cout<<"Error? Asking TmpManager::isuptodateR("<<site<<"), when "
	<<"length="<<len<<endl;
#ifdef TESTFUN
    cout<<*this<<endl;
#endif
    exit(1);
  }
  return uptodateR[site];
}

void TmpManagerExc::outofdateL(int site){
  if((site<0)||(site>=len)){
    cout<<"Error? TmpManager::outofdateL("<<site<<"), when "
	<<"length="<<len<<endl;
    exit(1);
  }
  if(site>0){
    //cout<<"TmpMgr:outofdateL "<<site<<endl;
    if(uptodateL[site]){
      //delete [] operL[site];
      //delete [] normL[site];
    }
    uptodateL[site]=0;
#ifdef DISKTMPMGR
    operL.discard(site);
    if(gev) normL0.discard(site);
    for(int k=0;k<nextK;k++)
      normL[k].discard(site);
#endif //DISKTMPMGR
    // Recursively outofdate Left terms to the right of this one  
    // is done from Contractor, that knows the logic and the dependencies
  }
}

void TmpManagerExc::outofdateR(int site){
  if((site<0)||(site>=len)){
    cout<<"Error? TmpManager::outofdateR("<<site<<"), when "
	<<"length="<<len<<endl;
    exit(1);
  }
  if(site<len-1){
    //cout<<"TmpMgr:outofdateR "<<site<<endl;
    if(uptodateR[site]){
      //delete [] operR[site];
      //delete [] normR[site];
    }
    uptodateR[site]=0;
#ifdef DISKTMPMGR
    operR.discard(site);
    if(gev) normR0.discard(site);
    for(int k=0;k<nextK;k++)
      normR[k].discard(site);
#endif //DISKTMPMGR
  }
}



// For Contractor!

void TmpManagerExc::storedL(int site){
  if(uptodateL[site]){
    printf("WARNING: I'm storing new Left terms for site %d, which was already uptodate!\n",site);
    exit(212);
    // TODO: EXIT?
  }
  uptodateL[site]=1;
#ifdef DISKTMPMGR
  operL.store(site,0);
  if(gev) normL0.store(site,0);
  for(int k=0;k<nextK;k++)
    normL[k].store(site,0);
  if(site>0){
    operL.discard(site-1); // the one to the left, I'm probably not using
    if(gev) normL0.discard(site-1);
    for(int k=0;k<nextK;k++)
      normL[k].discard(site-1);
  }
#endif //DISKTMPMGR
  // TODO: Recursively outofdate Left terms to the right of this one
}

// For Contractor!

void TmpManagerExc::storedR(int site){
  if(uptodateR[site]){
    printf("WARNING: I'm storing new Right terms for site %d, which was already uptodate!\n",site);
    exit(212);
    // TODO: EXIT?
  }
  uptodateR[site]=1;
#ifdef DISKTMPMGR
  operR.store(site,0); // try not to remove too early
  if(gev) normR0.store(site,0); 
  for(int k=0;k<nextK;k++)
    normR[k].store(site,0); 
  if(site<len-1){
    operR.discard(site+1); 
    if(gev) normR0.discard(site+1); 
    for(int k=0;k<nextK;k++)
      normR[k].discard(site);
  }
#endif //DISKTMPMGR
  // TODO: Recursively outofdate Right terms to the left of this one
}


#ifdef TESTFUN
ostream& operator<<(ostream& os,const TmpManagerExc& tmp){
  os<<"TmpManagerExc ("<<tmp.len<<","<<tmp.getNextK()<<") uptodate:"<<endl;
  os<<"Left: uptodate = ";
  for(int k=0;k<tmp.len;k++){
    os<<tmp.isuptodateL(k)<<" ";
  }
  os<<endl;
  os<<"Right: uptodate = ";
  for(int k=0;k<tmp.len;k++){
    os<<tmp.isuptodateR(k)<<" ";
  }
  os<<endl;
  for(int pos=0;pos<tmp.len;pos++){
    os<<"Testing Operators (pos="<<pos<<")"<<endl;  
    os<<"operL="<<tmp.operL[pos]<<endl;
    os<<"operR="<<tmp.operR[pos]<<endl;
    if(gev){
      os<<"normL0="<<tmp.normL0[pos]<<endl;
      os<<"normR0="<<tmp.normR0[pos]<<endl;
    }
    os<<"Testing overlap terms"<<endl;
    for(int l=0;l<tmp.getNextK();l++)
      os<<"normL term "<<l<<tmp.normL[l][pos]<<endl;
    for(int l=0;l<tmp.getNextK();l++)
      os<<"normR term "<<l<<tmp.normR[l][pos]<<endl;
  }
    return os;
}
#endif
