#include "Indices.h"
#include "mwArray.h"
#ifdef DISKSTRG
#include "FileSite.h"
typedef FileSite Site_t ;
#else
#include "Site.h"
typedef Site Site_t ;
#endif
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>

using namespace std;
using namespace shrt;

/** 
    Test program Nr. 2: Check the correct behaviour of Site class
    WARNING: To try the gauge methods, they should be unprotected in Site.h
*/

int main(){


#ifdef DISKSTRG
  //  int errD=system("rm -rf /ptmp/mpq/banulsm/*");
  char tmpdir[140];
  sprintf(tmpdir,"tmpF/SitesXXXXXX");
  // if(jobname!=0)
  //   sprintf(tmpdir,"/ptmp/mpq/banulsm/Sites%sXXXXXX",jobname);
  // else
  //   sprintf(tmpdir,"/ptmp/mpq/banulsm/SitesXXXXXX");
  char* dirid=mkdtemp(tmpdir);
  if(dirid==0){
    cout<<"Error: couldn't apply mkdtemp "<<tmpdir<<endl;
    exit(1);
  }
  FileSite::setDir(tmpdir);
#endif

  //Create a Site_t
  int d=2;int Dl=2;int Dr=3;
  int pos=10;
  Site_t S(pos,d,Dl,Dr);
  cout<<"Constructed site S="<<S<<endl;

  // setting to product
  S.setState(randmps);
  cout<<"Site_t S set to random ="<<S<<endl;

  // rotating
  vector<int> dims(3);dims[0]=3;dims[1]=2;dims[2]=3;
  mwArray A(dims);A.fillRandom();
  
  Site_t S2(pos,2,3,3);
  Indices order(2,3,1);
  S2.setRotatedSite(A,order);
  cout<<"From mwArray A="<<A<<endl;
  cout<<"Set rotated site S2="<<S2<<endl;

  // Test gauge conditions
  // Left
  mwArray tensor,leftterm1(ONE_c);
  mwArray leftterm3(Indices(3,3));
  leftterm3.fillRandom();
  S2.gaugeR(tensor,leftterm3);

  cout<<"GaugeR condition with leftterm="<<leftterm3<<endl;
  cout<<"After applying gauge condition S2="<<S2<<endl;
  cout<<"Extra tensor="<<tensor<<endl;
  // And right
  Site_t S3(pos,2,3,3); S3.setRotatedSite(A,order);
  S3.gaugeL(tensor,leftterm3);
  cout<<"GaugeL condition with rightterm="<<leftterm3<<endl;
  cout<<"After applying gauge condition S3="<<S3<<endl;
  cout<<"Extra tensor="<<tensor<<endl;

  // Test mixed gauge condition: need two sites satisfying right one
  int Dl4(3),Dr4(3);
  Site_t S4(pos,2,Dl4,Dr4);S4.setState(randmps);
  S4.gaugeR(tensor,leftterm1);
  cout<<"After applying gaugeR condition S4="<<S4<<endl;

  mwArray leftterm4(Indices(Dl4,Dl4));leftterm4.fillRandom();
  mwArray rightterm2,rightterm4;
  Site_t::gaugeR(&S2,&S4,rightterm2,rightterm4,leftterm3,leftterm4);

  cout<<"Used for gauge condition letterm3="<<leftterm3<<endl;
  cout<<"Used for gauge condition letterm4="<<leftterm4<<endl;

  cout<<"After applying gauge condition S2(ket)="<<S2<<endl;
  cout<<"After applying gauge condition S4(bra)="<<S4<<endl;
  cout<<"After applying gauge condition rT2(ket)="<<rightterm2<<endl;
  cout<<"After applying gauge condition rT4(bra)="<<rightterm4<<endl;

  // Finally, check the silly increase and decrease functions
  S4.increasePhysDim(3);
  cout<<"After increasing d->3, S4="<<S4<<endl;
  S4.decreaseBondDim(2,2);
  cout<<"After decreasing Dl,Dr->2, S4="<<S4<<endl;
  S4.increaseBondDim(2,3);
  cout<<"After increasing Dr->3, S4="<<S4<<endl;

  // Testing save and load
  char filename[]="dataSite_t.dat";
  ofstream out(filename,ios::binary);
  if(!out.is_open()){
    cout<<"Error: could not open to write"<<endl;
  }
  else{
    S4.save(out);
    out.close();
  }
  cout<<"S4 saved to file "<<filename<<endl;
  ifstream inf(filename,ios::binary);
  if(!inf.is_open()){
    cout<<"Error: could not open to read"<<endl;
  }
  else{
    Site_t S5(0,1,1,1);
    S5.load(inf);
    inf.close();
    cout<<"Read from file S5 "<<S5<<endl;
    cout<<"Difference with the original: "<<S4.getA()-S5.getA()<<endl;
  }

  cout<<"Now testing gauge condition in a special case"<<endl;
  Site_t S6(0,1,1,1);
  S6.setState(randmps);
  mwArray auxT;
  mwArray leftterm(ONE_c);
  cout<<"Site is "<<S6<<", A="<<S6.getA()<<endl;
  S6.gaugeR(auxT,leftterm);
  cout<<"After appling gaugeR, Site_t is "<<S6<<", A="<<S6.getA()
      <<", and tensor for next "<<auxT<<endl;
  mwArray aux2;
  S6.gaugeL(aux2,auxT);
  cout<<"After appling gaugeL, Site_t is "<<S6
      <<", and tensor "<<aux2<<endl;
    
}
