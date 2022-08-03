#include "TmpManager.h"

TmpManager::TmpManager():
  len(0),gauge(0),resolvent(0),exponential(0),
  exponentialOff(0),operL(),normL(),
  operRL(),operEL(),operR(),normR(),
  operRR(),operER(),operLoff(),operRoff(),
  uptodateL(NULL),uptodateR(NULL){};

TmpManager::TmpManager(const TmpManager& tmp):
  len(0),gauge(0),resolvent(0),exponential(0),
  exponentialOff(0),operL(),normL(),
  operRL(),operEL(),operR(),normR(),
  operRR(),operER(),operLoff(),operRoff(),
  uptodateL(NULL),uptodateR(NULL){};

TmpManager& TmpManager::operator=(const TmpManager& tmp){
  return *this;
}


TmpManager::TmpManager(int length,bool gauge_,bool resolvent_,
		       bool exponential_,bool exponentialOff_):
  len(length),gauge(gauge_),resolvent(resolvent_),exponential(exponential_),
  exponentialOff(exponentialOff_),
  uptodateL(NULL),uptodateR(NULL),
  operL(length,mwArray::empty),normL(length,mwArray::empty),
  operRL(length,mwArray::empty),operEL(length,mwArray::empty),
  operR(length,mwArray::empty),normR(length,mwArray::empty),
  operRR(length,mwArray::empty),operER(length,mwArray::empty),
  operLoff(length,mwArray::empty),operRoff(length,mwArray::empty){
  if(exponentialOff)exponential=1;
  if(exponential)resolvent=1;
  if(resolvent)gauge=0;
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
  normL[0]=edge;
  operRL[0]=edge;
  operEL[0]=edge;
  operLoff[0]=edge;
  operR[length-1]=edge;
  normR[length-1]=edge;
  operRR[length-1]=edge;
  operER[length-1]=edge;
  operRoff[length-1]=edge;
  if(!exponentialOff){
    operLoff.clear();
    operRoff.clear();
    if(!exponential){
      operEL.clear();
      operER.clear();
      if(!resolvent){
	operRL.clear();
	operRR.clear();
      }
    }
  }
#ifdef TESTFUN
  //  cout<<"Test(just constructed TMpMgr): "<<*this<<endl;
#endif

}


TmpManager::~TmpManager(){
  clear();
  delete []uptodateL;
  delete []uptodateR;
  //  cout<<"Destructing TmpManager"<<endl;
}

void TmpManager::clear(){
  // TODO: Delete things that might be kept uptodate till now
  operL.clear();
  normL.clear();
  operR.clear();
  normR.clear();
  operRL.clear();
  operRR.clear();
  operEL.clear();
  operER.clear();
  operLoff.clear();
  operRoff.clear();
}

bool TmpManager::isuptodateL(int site) const{
  return uptodateL[site];
}

bool TmpManager::isuptodateR(int site) const{
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

void TmpManager::outofdateL(int site){
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
    if(exponentialOff) operLoff.discard(site);
    if(exponential) operEL.discard(site);
    if(resolvent) operRL.discard(site);
    if(!gauge) normL.discard(site);
#endif //DISKTMPMGR
    // Recursively outofdate Left terms to the right of this one  
    // is done from Contractor, that knows the logic and the dependencies
  }
}

void TmpManager::outofdateR(int site){
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
    if(exponentialOff) operRoff.discard(site);
    if(exponential) operER.discard(site);
    if(resolvent) operRR.discard(site);
    if(!gauge) normR.discard(site);
#endif //DISKTMPMGR
  }
}



// For Contractor!

void TmpManager::storedL(int site){
  if(uptodateL[site]){
    printf("WARNING: I'm storing new Left terms for site %d, which was already uptodate!\n",site);
    exit(212);
    // TODO: EXIT?
  }
  uptodateL[site]=1;
#ifdef DISKTMPMGR
  operL.store(site,0);
  if(exponentialOff) operLoff.store(site,0);
  if(exponential) operEL.store(site,0);
  if(resolvent) operRL.store(site,0);
  if(!gauge) normL.store(site,0);
  if(site>0){
    operL.discard(site-1); // the one to the left, I'm probably not using
    if(exponentialOff) operLoff.discard(site-1);
    if(exponential) operEL.discard(site-1);
    if(resolvent) operRL.discard(site-1);
    if(!gauge) normL.discard(site-1);
  }
#endif //DISKTMPMGR
  // TODO: Recursively outofdate Left terms to the right of this one
}

// For Contractor!

void TmpManager::storedR(int site){
  if(uptodateR[site]){
    printf("WARNING: I'm storing new Right terms for site %d, which was already uptodate!\n",site);
    exit(212);
    // TODO: EXIT?
  }
  uptodateR[site]=1;
#ifdef DISKTMPMGR
  operR.store(site,0); // try not to remove too early
  if(exponentialOff) operRoff.store(site,0); 
  if(exponential) operER.store(site,0); 
  if(resolvent) operRR.store(site,0); 
  if(!gauge) normR.store(site,0); 
  if(site<len-1){
    operR.discard(site+1);
    if(exponentialOff) operRoff.discard(site);
    if(exponential) operER.discard(site);
    if(resolvent) operRR.discard(site);
    if(!gauge) normR.discard(site);
  }
#endif //DISKTMPMGR
  // TODO: Recursively outofdate Right terms to the left of this one
}

/* 
mwArray& TmpManager::getOperL(int pos){
  try{
    return operL.Get(1,pos+1);
  }
  catch(mwException& e){
    cout<<"Error in TmpManager getOperL: "<<e.what()<<endl;
    exit(212);
  }
}
mwArray& TmpManager::getOperR(int pos){
  try{
    return operR.Get(1,pos+1);
  }
  catch(mwException& e){
    cout<<"Error in TmpManager getOperR: "<<e.what()<<endl;
    exit(212);
  }
}
mwArray& TmpManager::getNormL(int pos){
  try{
    return normL.Get(1,pos+1);
  }
  catch(mwException& e){
    cout<<"Error in TmpManager getNormL: "<<e.what()<<endl;
    exit(212);
  }
}

mwArray& TmpManager::getNormR(int pos){
  try{
    return normR.Get(1,pos+1);
  }
  catch(mwException& e){
    cout<<"Error in TmpManager getNormR: "<<e.what()<<endl;
    exit(212);
  }
}
*/
#ifdef TESTFUN
ostream& operator<<(ostream& os,const TmpManager& tmp){
  os<<"TmpManager ("<<tmp.len<<") uptodate:"<<endl;
  os<<"Left:";
  for(int k=0;k<tmp.len;k++){
    os<<tmp.isuptodateL(k)<<" ";
  }
  os<<endl;
  os<<"Right:";
  for(int k=0;k<tmp.len;k++){
    os<<tmp.isuptodateR(k)<<" ";
  }
  os<<endl;
  int pos=0;
  os<<"Test Operators! (pos="<<pos<<")"<<endl;  
  os<<tmp.operL[pos]<<endl;
  if(!tmp.gauge)
    os<<tmp.normL[pos]<<endl;
  os<<tmp.operR[pos]<<endl;
  if(!tmp.gauge)
    return os;
  else
    return os<<tmp.normR[pos]<<endl;
}
#endif
