
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"
//#include "HermitianTensorMultiplier.h"

#include "SchwingerHamiltonianSz.h"
#include "SpinMPO.h"
#include <cmath>
#include <unistd.h>

#include "quicksort.h"

using namespace std;
using namespace shrt;

#define MAXLEN 120

/** 
Generates a file name for the MPS corresponding to a certain excitation.
*/
const string mpsfilename(int L,int D,double x,double mg,int level,const string mpsdir);
/** 
Generate the name of the file with the command in it
*/
const string jobfilename(int L,int D,double x,double mg,double tol,const string mpsdir);

/** Check if the filename exists */
bool file_exists(const string filename);

int getNextD(int D);

/** Auxiliary object, level-overlap, to sort */
class Level{
public:
  int lev;
  double energy;
  Level(int l, double E):lev(l),energy(E){};
  Level(const Level& l2):lev(l2.lev),energy(l2.energy){};
  bool operator<=(const Level& l2){return energy<=l2.energy;}
  bool operator<(const Level& l2){return energy<l2.energy;}
  bool operator==(const Level& l2){return energy==l2.energy;}
  friend bool operator<=(const Level& l1,const Level& l2){return l1.energy<=l2.energy;}
  friend bool operator<(const Level& l1,const Level& l2){return l1.energy<l2.energy;}
  friend bool operator==(const Level& l1,const Level& l2){return l1.energy==l2.energy;}
};


/** After something has been computed, this program looks at smaller D
    states and the computed ones to determine which ones to use as
    seed. Then prepares a searchO run to fill in the missing levels.

    Receives as argument a Properties file, such as
    config/excitProp.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.

    Property minnumber determines how many states are to be computed,
    so how many links will be prepared.

    Existing MPS files with smaller D are read from another directory,
    if they do not exist in the destination one, to be used as
    initial states. For this, specify the name of the old directory
    initmpsdir in Properties. If the files for a smaller D are to be
    used, specify the value in initD (in that case, initmpsdir could
    be the same as mpsdir). 

    Property jobsdir will indicate where to write the file with the 
    command that has to be relaunched 

    To replace some of the values in the file by a command line
    option, the additional arguments have to be of the form
    -propName=propValue
    (except for the vector properties)
*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  string directory="";
  // if(argc>2){
  //   const char* indir=argv[++cntr];
  //   directory=string(indir)+"/";
  // }
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  // Store the command line, for this is what should be relaunched for next D
  stringstream commandstr;
  commandstr<<"#!/bin/bash"<<endl;
  commandstr<<"# ";
  for(int c=0;c<argc;c++) commandstr<<argv[c]<<" ";
  // Now, str contains the command line, that I need to relaunch, if this was not the last run



  int L = props.getIntProperty("L");
  double mg=props.getDoubleProperty("mg");
  double alpha=props.getDoubleProperty("alpha");
  double x=props.getDoubleProperty("x");
  //  int D0 = props.getIntProperty("D0");
  //  int nLev = props.getIntProperty("nLev");
  directory=props.getProperty("outputdir");
  //const string outfile=directory+"/"+props.getProperty("outfile");
  int app=props.getIntProperty("append");
  if(app<0) app=1; // default will be to append
  double tol=props.getDoubleProperty("tol");
  // int unif=1; // All same D 
  // int unif=props.getIntProperty("unif");
  int D=props.getIntProperty("D");
  double zpen=props.getDoubleProperty("penalty");
  double offset=props.getDoubleProperty("offset");
  string mpsdir=props.getProperty("mpsdir"); // default current
  if(mpsdir.empty()) 
    mpsdir=".";
  int minNumber=props.getIntProperty("minnumber");

  string initmpsdir=props.getProperty("initmpsdir"); // default none
  if(initmpsdir.empty()){
    initmpsdir=props.getProperty("oldmpsdir"); // default none
  }
  int initD=props.getIntProperty("initD");
  if(!initmpsdir.empty()){
    if(initD<0){
      cout<<"ERROR: No initD specified"<<endl;
      exit(1);
    }
  }
  else{
    cout<<"ERROR: No initmpsdir specified!"<<endl;
    exit(1);
  }


  string jobsdir=props.getProperty("jobsdir"); 
  int maxTime=props.getIntProperty("maxTime");
  if(maxTime==-1) maxTime=0;


  // int noCopy=props.getIntProperty("nocopy");
  // bool sortAndCopy=true;
  // if(noCopy!=-1) sortAndCopy=false;

  // Parameters I do not use here: chemical potential and weight of the gauge term
  double nu=0;
  // I might want to switch off the gauge terms
  //double gweight=props.getDoubleProperty("gweight");
  //if(gweight==-1) gweight=1; // default: do not alter
  double gweight=1.; // Not touching the gauge terms (std)

  double mu=2*mg*sqrt(x);
  double L0=0.; // l0 parameter, which is irrelevant, as only appears
		// as alpha+l0

  //  double t0=.01; // todo: differently -> reference time

  cout<<"Initialized arguments: L="<<L
      <<", mu="<<mu
      <<", x="<<x
      <<", alpha="<<alpha
    //<<", outfile="<<outfile
      <<", app="<<app
      <<", tol="<<tol
      <<", zpen="<<zpen
      <<", gweight="<<gweight
      <<", D0="<<D<<endl;

#ifdef DISKSTRG
  //  int errD=system("rm -rf /ptmp/mpq/banulsm/*");
  char tmpdir[150];
  sprintf(tmpdir,"/ptmp/mpq/banulsm/Schwinger/Sitesx%dL%dD%dmg%dXXXXXX",x,L,D,(int)(mg*1000));
  char* dirid=mkdtemp(tmpdir);
  if(dirid==0){
    cout<<"Error: couldn't apply mkdtemp "<<tmpdir<<endl;
    exit(1);
  }
  FileSite::setDir(tmpdir);
#endif

  // ofstream* out;

  // out=new ofstream(outfile.data(),ios::app);
  // if(!out->is_open()){
  //   cout<<"Error: impossible to open file "<<outfile<<
  //     " for output (append="<<app<<")"<<endl;
  //   exit(1);
  // }
  // *out<<setprecision(15);

  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);
  int d=2;
  // // First: put H a
  // int Dop=5; // bond dimension of Schwinger H
  // double lambda=0.;
  //SchwingerHamiltonianSz hSch(L,mu,x,L0,alpha,zpen,offset,nu,gweight);
  //const MPO& hamil=hSch.getHMPO();
  // Observables I will compute
  // Compute the expectation value of Sz, as it commutes with H
  // MPO Szmpo(L);
  // SpinMPO::getSzMPO(L,d,Szmpo);
  // // Sz2 is the penalty term (we need to substract it) 
  // MPO Sz2mpo(L);
  // SpinMPO::getSz2MPO(L,d,Sz2mpo);
  // // The momentum operator <P2> is what we want
  // MPO Pmpo(L);
  // hSch.constructMomentumMPO(Pmpo);
  // const MPO* oprsP[2]={&Pmpo,&Pmpo};
  // MPO P2mpo(L);
  // MPO::join(2,oprsP,P2mpo);
  // // // The momentum squared operator
  // // MPO PSqmpo(L);
  // // hSch.constructMomentumSquaredMPO(PSqmpo);

  // // the shift by one and rotate sigmax
  // MPO SRmpo(L);
  // hSch.constructShiftRotateMPO(SRmpo);
  // //  SRmpo.exportForMatlab("SRmpo.m");

  // // The fermionic condensate and the Gamma_alpha, Gamma_5 parameters
  // MPO condMPO(L),GammaA(L),Gamma5(L);
  // hSch.constructCondensateMPO(condMPO);
  // hSch.constructGammaAlphaMPO(GammaA);
  // hSch.constructGamma5MPO(Gamma5);

  vector<MPS*> oldlevels;
  vector<MPS*> levels;

  // Read the already existing levels in mpsdir
  bool oldMPS=true;
  int lval=0;
  MPS* nextLev;
  while(oldMPS){
    const string filename=mpsfilename(L,D,x,mg,lval,mpsdir); // Try the definitive one
    if(file_exists(filename)){
      nextLev=new MPS(L,D,d);
      nextLev->importMPS(filename.data());
      levels.push_back(nextLev);
      cout<<"Recovered MPS file for level "<<lval<<" from file "<<filename<<endl;
      lval++;
    }
    else
      oldMPS=false; // break the loop 
  }

  // Now I already have read all existing MPS with bond D
  cout<<"Recovered "<<lval<<" already computed MPS levels with D="<<D<<" from disk. "<<endl;

  // Now do the same for the specified number of smaller D levels
  oldMPS=true;
  int l=0;
  while(l<=minNumber&&oldMPS){
    const string filename=mpsfilename(L,initD,x,mg,l,initmpsdir); // Try the definitive one
    if(file_exists(filename)){
      nextLev=new MPS(L,initD,d);
      nextLev->importMPS(filename.data());
      oldlevels.push_back(nextLev);
      cout<<"Recovered MPS file for level "<<l<<" from file "<<filename<<endl;
      l++;
    }
    else
      oldMPS=false; // break the loop 
  }
  cout<<"Recovered "<<l<<" already computed MPS levels with D="<<initD<<" from disk. "<<endl;

  // Now compute the overlaps (I write them (abs) in a matrix)
  mwArray overlaps(Indices(lval,l));
  for(int i1=0;i1<lval;i1++){
    for(int i2=0;i2<l;i2++){
      complex_t res=contractor.contract(*levels[i1],*oldlevels[i2]);
      overlaps.setElement(res,Indices(i1,i2));
    }
  }

  cout<<"Overlap matrix is "<<overlaps<<endl;

  // Since I have lval levels computed with D, and l>lval with Dref, I
  // will set l-lval links to use as the next seeds. I will select the
  // most orthogonal levels to the ones computed, excluding the scalar
  // candidate, though.



  // Now from the overlap matrix, identify the mos t likely matches for the computed levels
  vector<Level> matches;
  for(int i2=2;i2<l-1;i2++){
    double tmp=0.;
    for(int i1=1;i1<lval;i1++){
      complex_t aux=overlaps.getElement(Indices(i1,i2));
      tmp+=abs(aux)*abs(aux);
    }
    matches.push_back(Level(i2,tmp));
  }  
  int remain=l-lval;
  // pick the smallest l-lval
  quicksort(matches,0,matches.size()-1);
  
  if(l-lval>0)
    cout<<"The minimum overlap (excluded 0,1, and scalar) corresponds to levels: "<<endl;
  for(int k=0;k<l-lval;k++){
    cout<<"l="<<matches[k].lev<<"-ov="<<matches[k].energy<<","<<endl;
    // Suggested link
    const string filenameOld=mpsfilename(L,initD,x,mg,matches[k].lev,initmpsdir);
    const string filenameNew=mpsfilename(L,D,x,mg,lval+k,mpsdir);
    cout<<"ln -s "<<filenameOld<<" "<<filenameNew<<".tmp"<<endl;
    commandstr<<"ln -s "<<filenameOld<<" "<<filenameNew<<".tmp"<<endl;
  }

  // And finally launch the searchO program to refine the states
  const string jobfile=jobfilename(L,D,x,mg,tol,jobsdir);
  int mem=ceil(L*D*D*(12.+4.*minNumber)*16E-6)+200;
  commandstr<<"msub_modified_tqo097 -N D"<<D<<".M"<<props.getProperty("mg")<<"."<<x<<"."<<L<<" -l h_vmem="<<mem<<"M -raw ";
  commandstr<<"./single"<<" "<<argv[1];
  commandstr<<" -L="<<L<<" -mg="<<mg<<" -alpha="<<alpha;
  commandstr<<" -x="<<x<<" -D="<<D;
  commandstr<<" -append="<<app<<" -outfile=L"<<L<<"x"<<x<<"D"<<D;
  commandstr<<" -outputdir="<<directory<<" -tol="<<tol;
  commandstr<<" -penalty="<<zpen;
  commandstr<<" -offset="<<offset;
  commandstr<<" -mpsdir="<<mpsdir;
  commandstr<<" -jobsdir="<<jobsdir;
  commandstr<<" -maxTime="<<maxTime;
  commandstr<<" -minnumber="<<minNumber;
  commandstr<<" -msubmem="<<mem<<"M";
  ofstream* ojob=new ofstream(jobfile.data());  
  if(!ojob->is_open()){
    cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
    cout<<commandstr.str()<<endl;
    cout<<endl;
  }
  else{
    *ojob<<commandstr.str()<<endl;
    ojob->close();
  }
  delete ojob;


  // Clear the MPS vectors
  while(levels.size()>1){
    delete levels.back();
    levels.pop_back();
  }  

  while(oldlevels.size()>1){
    delete oldlevels.back();
    oldlevels.pop_back();
  }  

}

// void nada(){
//   // Write the next one to a file, so that it is relaunched
//   int nextD=getNextD(D);
//   if(nextD>0){
//     const string jobfile=jobfilename(L,D,x,mg,tol,jobsdir);
//     stringstream commandstr;
//     commandstr<<argv[0]<<" "<<argv[1];
//     commandstr<<" -L="<<L<<" -mg="<<mg<<" -alpha="<<alpha;
//     commandstr<<" -x="<<x<<" -D="<<nextD;
//     commandstr<<" -append="<<app<<" -outfile=L"<<L<<"x"<<x<<"D"<<nextD;
//     commandstr<<" -outputdir="<<directory<<" -tol="<<tol;
//     commandstr<<" -penalty="<<zpen;
//     commandstr<<" -offset="<<offset;
//     commandstr<<" -mpsdir="<<mpsdir;
//     commandstr<<" -initmpsdir="<<mpsdir; // what I just computed is the new initmpsdir
//     commandstr<<" -initD="<<D<<" ";
//     commandstr<<" -noise="<<noise;
//     commandstr<<" -jobsdir="<<jobsdir;
//     commandstr<<" -maxTime="<<maxTime;
//     commandstr<<" -minnumber="<<minNumber;
//     commandstr<<" -msubmem="<<ceil(L*nextD*nextD*(12.+4.*minNumber)*16E-6)+200<<"M";
//     ofstream* ojob=new ofstream(jobfile.data());
//     if(!ojob->is_open()){
//       cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
//       cout<<commandstr.str()<<endl;
//       cout<<endl;
//       // Do not exit, because I still need to save the MPS
//     }
//     else{
//       *ojob<<commandstr.str()<<endl;
//       ojob->close();
//     }
//     delete ojob;
//   }
//   else{ // no more launching of searchO
//     if(D==20){
//       nextD=40;
//       // Prepare the true simulation with larger Ds and program single
//       const string jobfile=jobfilename(L,D,x,mg,tol,jobsdir);
//       stringstream commandstr;
//       commandstr<<" ./single "<<argv[1];
//       commandstr<<" -L="<<L<<" -mg="<<mg<<" -alpha="<<alpha;
//       commandstr<<" -x="<<x<<" -D="<<nextD;
//       commandstr<<" -append="<<app<<" -outfile="<<props.getProperty("outfile");
//       commandstr<<" -outputdir="<<directory<<" -tol="<<tol;
//       commandstr<<" -penalty="<<zpen;
//       commandstr<<" -offset="<<offset;
//       commandstr<<" -mpsdir="<<mpsdir;
//       commandstr<<" -initmpsdir="<<mpsdir; // what I just computed is the new initmpsdir
//       commandstr<<" -initD="<<D<<" ";
//       commandstr<<" -noise="<<noise;
//       commandstr<<" -jobsdir="<<jobsdir;
//       commandstr<<" -maxTime="<<maxTime;
//       commandstr<<" -msubmem="<<ceil((double)L*nextD*nextD*(12.+4.*2)*16.0E-6+200.0)<<"M";
//       ofstream* ojob=new ofstream(jobfile.data());
//       if(!ojob->is_open()){
// 	cout<<"WARNING!!!! Could not open the JOB file. Launch by hand with the following command:"<<endl;
// 	cout<<commandstr.str()<<endl;
// 	cout<<endl;
// 	// Do not exit, because I still need to save the MPS
//       }
//       else{
// 	*ojob<<commandstr.str()<<endl;
// 	ojob->close();
//       }
//       delete ojob;
//     }
//   }
//    // for(int z=1;z<L;z++){
//    //   delete Skmpos[z];
//    // }
//    // delete []Skmpos;
// }

#include <sstream>    

int getNextD(int D){
  if(D==10) return 20;
  if(D==4) return 10;
  if(D==2) return 4;
  return -1;
}

const string mpsfilename(int L,int D,double x,double mg,int level,const string mpsdir){
  stringstream s;
  s<<mpsdir<<"/MPS_L"<<L<<"_m"<<mg<<"_x"<<x<<"_D"<<D<<"_l"<<level<<".dat";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

const string jobfilename(int L,int D,double x,double mg,double tol,const string jobsdir){
  stringstream s;
  s<<jobsdir<<"/job_L"<<L<<"_m"<<mg<<"_x"<<x<<"_D"<<D;
  if(tol!=1E-7) s<<"_tol"<<-log10(tol);
  s<<".sh";  
  //  cout<<"mpsfilename returns "<<s.str()<<endl;
  return s.str();
}

bool file_exists(const string filename)
{
  ifstream ifile(filename.data());
  //cout<<"file_exists("<<filename<<")="<<ifile<<endl;
  return ifile;
}
