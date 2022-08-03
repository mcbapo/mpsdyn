
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

void sumMPO(MPO& result,const MPO& mpo1,const MPO& mpo2);

//void computeVars(vector<int>& Ds,vector<double>& energies,vector<double>& vars,const string mpsprobelist,int L,int d,const MPO& hamil);
void computeVars(ofstream* out,const string mpsprobelist,const MPO& hamil,const string& header,int Lcmax,int Lb);
double computeRDMdistanceToId(const mwArray& rdm,int Lc);


/** Given a list of MPS in a file, evaluate and write the energies and variances. 
    I also compute distances of RDM in the center wrt RDM of the thermal state at 0 energy!!!
    (TODO: Should compare to thermal state at the FOUND energy)
 */

int d=2;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  int L=props.getIntProperty("L");
  double J_=props.getDoubleProperty("J");
  double h_=props.getDoubleProperty("h");
  double g_=props.getDoubleProperty("g");
  string outfname=props.getProperty("outputfile");
  bool app=(props.getIntProperty("append")==1); // has to be asked for explicitly
  string mpsprobelist=props.getProperty("mpslist"); // file with a list of MPS to use to probe the DoS
  int Lcmax=props.getIntProperty("Lcmax"); // max size of middle 
  if(Lcmax<0||Lcmax>10) Lcmax=6;
  int Lb=props.getIntProperty("Lb");
  if(Lb<0) Lb=6;
  
  cout<<"Initialized arguments: L="<<L
      <<", J="<<J_
      <<", g="<<g_
      <<", h="<<h_
      <<", outfile="<<outfname<<endl;


  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;

  // The original Hamiltonian to get the energy
  IsingHamiltonian hamH0(L,d,J_,g_,h_);
  const MPO& hamil0=hamH0.getHMPO();
 
  
  // For output
  ofstream* out;
  if(mpsprobelist.size()>0){
    if(!file_exists(outfname)||!app){
      out=new ofstream(outfname.data());
      *out<<"% J\t g\t h\t D\t E(k)\t <H^2(k)> for k=1...."<<endl;
      if(!out->is_open()){
	cout<<"Error: impossible to open file "<<outfname<<
	  " for output"<<endl;
	exit(1);
      }
      out->close();delete out;
    }
  }
    
  // If appropriate, probe the distribution
  if(mpsprobelist.size()){
    // vector<double> energies;
    // vector<double> vars;
    // vector<int> Ds;
    // computeVars(Ds,energies,vars,mpsprobelist,L,d,hamil0);
    // save values to file
    // out=new ofstream(outfname.data(),ios::app);
    // for(int p=0;p<vars.size();p++){
    //   *out<<J_<<"\t"<<g_<<"\t"<<h_<<"\t";
    //   *out<<Ds[p]<<"\t";
    //   *out<<energies[p]<<"\t";
    //   *out<<vars[p]<<"\t";
    //   *out<<endl;
    // }
    stringstream strstr;
    strstr<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"; // beginning of each line    
    // compute and save values to file
    out=new ofstream(outfname.data(),ios::app);
    computeVars(out,mpsprobelist,hamil0,strstr.str(),Lcmax,Lb);
    out->close();
    delete out;
  }

  
}

void computeVars(ofstream* out,const string mpsprobelist,const MPO& hamil,const string& header,int Lcmax,int Lb){
  vector<string> list;
  int nr=0;
  ifstream input(mpsprobelist.data());
  if(!input.is_open()){
    cout<<"Error trying to read list of files from "<<mpsprobelist<<endl;
    exit(1);
  }
  string s;
  while(input.is_open()){
    getline(input,s);
    if(input.good()){
      if(s.length()>0&&s[0]!='%'){ // ignore comments and empty lines
	list.push_back(s);
	nr++;
	cout<<"Read file name: "<<s<<endl;
      }
    }
    else{ // no more lines
      input.close();
      cout<<"Closed the list file"<<endl;
    }
  }
  Contractor& contractor=Contractor::theContractor();
  for(int k=0;k<nr;k++){
    // Read and evaluate state nr
    MPS mps(1,1,1);
    mps.importMPS(list[k].data());
    int L=mps.getLength();
    // The D of this guy
    int Dk=mps.getA(floor(L/2)).getDl();
    complex_t norm=contractor.contract(mps,mps);
    complex_t Ek=contractor.contract(mps,hamil,mps);
    complex_t E2k=contractor.contract2(hamil,mps);
    cout<<"Probe state "<<k<<" has D="<<Dk<<", norm="<<norm<<", energy "<<Ek<<" and <H^2>="<<E2k
	<<", thus variance="
	<<E2k-Ek*Ek<<endl;
    *out<<header<<Dk<<"\t"<<real(Ek/norm)<<"\t"<<real(E2k/norm)<<"\t";

    vector<double> trDistRDM;
    vector<int> sizeRDM;

    // averaging over different subsystems of the same size
    for(int lc=1;lc<=Lcmax;lc++){sizeRDM.push_back(lc);} // position lc-1
    trDistRDM=vector<double>(sizeRDM.size(),0.);
    vector<int> avrgCnt(sizeRDM.size(),0.);
    int posL=Lb;
    while(posL<L-Lb-1){ // at least there is a 2-site RDM to compute
      int lc=Lcmax;
      int posR=posL+lc-1;
      while(posR>=L){ // cannot compute this one, need to reduce the size
	lc--;
	posR=posL+lc-1;
      }
      int dimL=pow(d,lc);
      int dimR=dimL;
      mwArray oper=contractor.getRDM(mps,posL,posR);
      while(lc>=1){ // Compute the distance and trace another site
	if(posR<L-Lb){
	  // cout<<"Computing distance for rmd("<<lc<<") sites, pos:"<<posL<<"-"<<posR<<endl;
	  double dist_lc=computeRDMdistanceToId(oper,lc);
	  trDistRDM[lc-1]+=dist_lc;
	  avrgCnt[lc-1]++; // to take the average later
	}
	if(lc>=1){
	  // trace out two sites
	  oper.reshape(Indices(dimL/d,d,dimR/d,d));
	  oper.permute(Indices(1,3,2,4));
	  dimL=dimL/d;dimR=dimR/d;
	  oper.reshape(Indices(dimL*dimR,d*d));
	  oper.multiplyRight(reshape(identityMatrix(d),Indices(d*d,1)));
	  oper.reshape(Indices(dimL,dimR));
	  posR--;
	}
	lc--;
      }	    
      posL++;
    }
    for(int lc=1;lc<=Lcmax;lc++) trDistRDM[lc-1]=trDistRDM[lc-1]/avrgCnt[lc-1];
   
    for(int ks=0;ks<sizeRDM.size();ks++){
      *out<<sizeRDM[ks]<<"\t"<<trDistRDM[ks]<<"\t";
    }
    *out<<endl;
  }
  
}

double computeRDMdistanceToId(const mwArray& rdm,int Lc){
  double singvals=0.;
  double idFc=1./pow(2,Lc);
  mwArray U;vector<complex_t> Dv;
  wrapper::eig(rdm,Dv,U); // only eigenvalues
  // now sum the absolute values
  for(int k=0;k<Dv.size();k++){
    double aux=real(Dv[k]);
    singvals+=abs(aux-idFc);
  }
  if(Dv.size()<pow(2,Lc)){
    cout<<"WARNING! Seems the dimension of the operator for Lc="<<Lc
	<<"is not full: dim(oper)="<<Dv.size()
	<<", 2^Lc="<<pow(2,Lc)<<endl;
    singvals+=idFc*(pow(2,Lc)-Dv.size());
  }
  return singvals;
}


// void computeVars(vector<int>& Ds,vector<double>& energies,vector<double>& vars,const string mpsprobelist,int L,int d,const MPO& hamil){
//   vector<string> list;
//   Ds.clear();
//   energies.clear();
//   vars.clear();
//   int nr=0;
//   ifstream input(mpsprobelist.data());
//   if(!input.is_open()){
//     cout<<"Error trying to read list of files from "<<mpsprobelist<<endl;
//     exit(1);
//   }
//   string s;
//   while(input.is_open()){
//     getline(input,s);
//     if(input.good()){
//       if(s.length()>0&&s[0]!='%'){ // ignore comments and empty lines
// 	list.push_back(s);
// 	nr++;
// 	cout<<"Read file name: "<<s<<endl;
//       }
//     }
//     else{ // no more lines
//       input.close();
//       cout<<"Closed the list file"<<endl;
//     }
//   }
//   Contractor& contractor=Contractor::theContractor();
//   for(int k=0;k<nr;k++){
//     // Read and evaluate state nr
//     MPS mps(L,1,d);
//     mps.importMPS(list[k].data());
//     // The D of this guy
//     int Dk=mps.getA(floor(L/2)).getDl();
//     Ds.push_back(Dk);
//     complex_t norm=contractor.contract(mps,mps);
//     complex_t Ek=contractor.contract(mps,hamil,mps);
//     energies.push_back(real(Ek/norm));
//     complex_t E2k=contractor.contract2(hamil,mps);
//     vars.push_back(real(E2k/norm));
//     cout<<"Probe state "<<k<<" has D="<<Dk<<", energy "<<Ek<<" and <H^2>="<<E2k
// 	<<", thus variance="
// 	<<E2k-Ek*Ek<<endl;
//   }
// }

