#include <math.h>
#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include <iomanip>
using namespace shrt;
using namespace std;


int d=2;
double tau=0.8;
bool traceless=true;

int main(int argc,const char* argv[]){


  int M = 6;
  int D = 10;
  int dD = 5;
  double noise = 3.5;

  //cout<<setprecision(15);
  
  
  Contractor& contractor=Contractor::theContractor(); //inizialize contractor

  double lambda,lambda1,lambda2;

  for(M=4;M<24;M++){  
    MPS Omps(M,D,4),OmpsBigD(M,D,4), OmpsBigDnoise(M,D,4);
    Omps.setRandomState();
    Omps.gaugeCond('R',1);Omps.gaugeCond('L',1); 
    OmpsBigD=Omps;
    OmpsBigDnoise=Omps;
    
 
    OmpsBigD.increaseBondDimension(D+dD);
    
    OmpsBigDnoise.increaseBondDimensionWithNoise(D+dD, noise);
    OmpsBigDnoise.gaugeCond('R',1);OmpsBigDnoise.gaugeCond('L',1);
    lambda=abs(contractor.contract(Omps,OmpsBigD));
    lambda1=abs(contractor.contract(Omps,OmpsBigDnoise));
    lambda2=abs(contractor.contract(OmpsBigD,OmpsBigDnoise));
    cout<<"<O(D)|O(D+dD)>: "<<lambda<<"\t <O(D)|noiseO(D+dD)>: "<<lambda1<<"\t <O(D+dD)|noiseO(D+dD)>: "<<lambda2<<endl;
  }
}

