
#include <math.h>
#include <iomanip>
#include <sstream>

#include "misc.h"
#include "Properties.h"
//#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "SpinMPO.h"
#include "quicksort.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "ThirringHamiltonian.h"

const string mpsTmpFileNameRoot(const string baseDir,int D,double g2,double gtilde2,double ma,double Starget,int L,double delta,double tol);
// Add a particular counter
const string mpsTmpFileName(const string rootName,int cnt);
// And a final name
const string mpsSaveFileName(const string rootName,double beta);

const string resultsFileName(const string rootName,double beta);

void computeCorrelators(const MPS& thS,int minX,int maxX,double beta,complex_t energy,ofstream* out);
void computeCorrelators_long(const MPS& thS,int minX,int maxX,double beta,complex_t energy,ofstream* out);

void contractL(mwArray& result,const mwArray& tmpL,
	       const Site& ket,const Site& bra,const mwArray& oper);
void contractR(mwArray& result,const mwArray& tmpR,
	       const Site& ket,const Site& bra,const mwArray& oper);

using namespace shrt;
using namespace std;

/** Complementary program to thermalthirringproj, this one reads the
    MPs files created by the other at the indicated betas and computes
    some observables in the thermal states. */




int d=2;
/** Global matrices */
complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
mwArray sigX(Indices(d,d),dataX);
complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
mwArray sigY(Indices(d,d),dataY);
complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
mwArray sigZ(Indices(d,d),dataZ);
mwArray sig0=identityMatrix(d);

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

  int L = props.getIntProperty("L");  
  if(L%2!=0){
    cout<<"ERROR: Nr of sites L should be even"<<endl;
    exit(1);
  }
  double ma=props.getDoubleProperty("ma");  
  double gt2_=props.getDoubleProperty("gtilde2");
  double g2_=props.getDoubleProperty("g2"); // the non-alternating one
  double lambda_=props.getDoubleProperty("lambda");
  if(lambda_<0) lambda_=0.;
  double mu_=props.getDoubleProperty("mu");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  int Stot_=props.getIntProperty("Starget");
  int D=props.getIntProperty("D");
  double beta=props.getDoubleProperty("beta");
  double delta=props.getDoubleProperty("delta");
  vector<double> targetBs=props.getVectorDoubleProperty("targets");
  // order, just in case
  quicksort(targetBs,0,targetBs.size()-1);
  
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir");

  bool app=(props.getIntProperty("append")!=0);
  //bool findGS=(props.getIntProperty("findGS")==1);

  int minX=props.getIntProperty("minX"); // min pos for a correlation 
  int maxX=props.getIntProperty("maxX"); // max pos for a correlation 
  if(minX<0) minX=0;
  if(maxX>0) maxX=L-1;
  
  cout<<"Initialized arguments: L="<<L
      <<", ma="<<ma
      <<", g2="<<g2_
      <<", gt2="<<gt2_
      <<", lambda="<<lambda_
      <<", mu="<<mu_
      <<", Starget="<<Stot_
      <<", outfile="<<outfname
      <<", beta="<<beta
      <<", delta="<<delta
      <<", target betas: "<<targetBs
      <<", app="<<app
      <<", D="<<D<<endl;

  int M=beta/(2*delta);

  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  cout<<"Initialized Contractor"<<endl;
  contractor.setConvTol(tol);

  ThirringHamiltonian hamH(L,ma,g2_,lambda_,Stot_,d,gt2_,mu_);
  //cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamH.getHMPO();
  
  MPO Hdoub(L); // Hamil on one side, nothing on the other
  MPO projS0(L); // projector onto total Sz=0 (on one side)
  //  MPO projS0_2(L); // projector on both sides: unnecessary!
  MPO Szdoub(L);
  MPO Hproj(L); // H on one side, projector on the other

  {
    MPO _projS0(L); // projector onto total spin 0 (on one layer)
    // Notice:this is a projector with exponentially large Frobenius norm choose(N,N/2)
    SpinMPO::getProjectorTotalSz(L,d,_projS0); // TODO: generalize for other total spins
    MPO Szmpo(L);
    SpinMPO::getSzMPO(L,d,Szmpo);
    mwArray identPhys=identityMatrix(d);
    identPhys.reshape(Indices(d,1,d,1)); 
    for(int k=0;k<L;k++){
      Hdoub.setOp(k,new DoubleOperator(hamil.getOp(k).getFullData(),identPhys),true);
      projS0.setOp(k,new DoubleOperator(_projS0.getOp(k).getFullData(),identPhys),true);
      mwArray aux=_projS0.getOp(k).getFullData();
      //      projS0_2.setOp(k,new DoubleOperator(_projS0.getOp(k).getFullData(),
      //				  permute(aux,Indices(3,2,1,4))),true);
      Hproj.setOp(k,new DoubleOperator(hamil.getOp(k).getFullData(),
				       permute(aux,Indices(3,2,1,4))),true);
      Szdoub.setOp(k,new DoubleOperator(Szmpo.getOp(k).getFullData(),identPhys),true);
    }
  }

  MPS thS(L,D,d*d);
  MPS thSneg(thS); // place to the states that I will be reading

  // First, create the root name for backup files
  string fileNameBasis=mpsTmpFileNameRoot(mpsdir,D,g2_,gt2_,ma,Stot_,L,delta,tol);

  for(int ptr=0;ptr<targetBs.size();ptr++){
    string nameMPSpos=mpsSaveFileName(fileNameBasis,targetBs[ptr]);
    //string nameMPSneg=mpsSaveFileName(fileNameBasis,-targetBs[ptr]);
    if(file_exists(nameMPSpos.data())){ //&&file_exists(nameMPSneg.data())){
      cout<<"Reading state for beta="<<targetBs[ptr]<<" from "<<nameMPSpos<<endl; //" and "<<nameMPSneg<<endl;
      thS.importMPS(nameMPSpos.data());
      thS.gaugeCond('L',1);
      //thSneg.importMPS(nameMPSneg.data());
      // Now compute stuff
      complex_t energy=contractor.contract(thS,Hproj,thS);
      cout<<"State has energy"<<energy<<endl;
      //complex_t weightS0=contractor.contract(thS,projS0,thS);
      //complex_t energyN=contractor.contract(thSneg,Hproj,thSneg);
      //complex_t weightS0N=contractor.contract(thSneg,projS0,thSneg);
      // file for this beta:
      string resFname=resultsFileName(outfname,targetBs[ptr]);
      ofstream* out;
      if(!app||!file_exists(outfname.data())){
	out=new ofstream(resFname.data());
	*out<<"% L="<<L<<", ma="<<ma<<", g2="<<g2_<<", gt2="<<gt2_<<", Starget="<<Stot_<<", D="<<D<<endl;
	*out<<"%beta\tE\tpos1\tpos2\t<Z(pos1)>\t<Z(pos2)>\t<Z1Z2>\t<Z1..Z..Z2>\t<sigmaPlus(pos1)..Z..sigmaMinus(pos2)>"<<endl;
      }
      else
	out=new ofstream(resFname.data(),ios::app);
      if(!out->is_open()){
	cout<<"Error: impossible to open file "<<resFname<<
	  " for output"<<endl;
	exit(1);
      }

      *out<<setprecision(10);
      
      // Correlators and single site expectation values
      computeCorrelators(thS,minX,maxX,targetBs[ptr],energy,out);
      //computeCorrelators(thS,minX,maxX,-targetBs[ptr],energyN,out);
      out->close();delete out;
    }
    else{
      cout<<"File(s) for beta="<<targetBs[ptr]<<" not found!!"<<endl;
    }
  }
}

void computeCorrelators_long(const MPS& thS,int minX,int maxX,double beta,complex_t energy,ofstream* out){
  // Not very efficient, but should work for moderate sizes
    // Single site ops have to be doubled with id on the other side
  mwArray ZI,PlusI,MinI; // sigz, sigPlus,sigMinus
  constructOperatorProduct(ZI,sigZ,sig0);
  constructOperatorProduct(PlusI,.5*(sigX+I_c*sigY),sig0);
  constructOperatorProduct(MinI,.5*(sigX-I_c*sigY),sig0);
  
  Contractor& contractor=Contractor::theContractor();
  complex_t norm=contractor.contract(thS,thS);
  for(int pos1=minX;pos1<maxX;pos1++){
    cout<<"pos1="<<pos1<<endl;
    MPS mpsZz(thS);
    MPS mpsPzM(thS);
    MPS mps(thS);
    mps.applyLocalOperator(pos1,ZI,false);
    complex_t valZ1=contractor.contract(mps,thS)/norm;
    mpsZz.applyLocalOperator(pos1,ZI,false);
    mpsPzM.applyLocalOperator(pos1,PlusI,false);
    for(int pos2=pos1+1;pos2<=maxX;pos2++){
      MPS aux(thS);
      aux.applyLocalOperator(pos2,ZI,false); // I could skip this one
      complex_t valZ2=contractor.contract(aux,thS)/norm;
      mps.applyLocalOperator(pos2,ZI,false); 
      complex_t valZ1Z2=contractor.contract(mps,thS)/norm;
      mps.applyLocalOperator(pos2,ZI,false); // remove Z2 for next pos2
      mpsZz.applyLocalOperator(pos2,ZI,false); // this is just the string
      complex_t valZstr=contractor.contract(mpsZz,thS)/norm;
      aux=mpsPzM;
      aux.applyLocalOperator(pos2,MinI,false);
      complex_t valPzM=contractor.contract(aux,thS)/norm;
      mpsPzM.applyLocalOperator(pos2,ZI,false); // for next position
      // write to file
      *out<<beta<<"\t"<<real(energy/norm)<<"\t"<<pos1<<"\t"<<pos2<<"\t";
      *out<<real(valZ1)<<"\t"; //<<imag(valZ1)<<"\t";
      *out<<real(valZ2)<<"\t"; //<<imag(valZ2)<<"\t";
      *out<<real(valZ1Z2)<<"\t"; //<<imag(valZ1Z2)<<"\t";
      *out<<real(valZstr)<<"\t"; //<<imag(valZstr)<<"\t";
      *out<<real(valPzM)<<"\t"<<imag(valPzM)<<"\t";
      *out<<endl;
    }    
  }
  
  
}



const string mpsTmpFileNameRoot(const string baseDir,int D,double g2,double gtilde2,double ma,double Starget,int L,double delta,double tol){
  stringstream s;
  s<<baseDir<<"/tmpMPS_D"<<D<<"_g2"<<g2<<"_gt2"<<gtilde2<<"_ma"<<ma
   <<"_S"<<Starget<<"_L"<<L<<"_delta"<<delta
   <<"_eps"<<int(-log10(tol));
  s<<"_";
  return s.str();
}

const string mpsSaveFileName(const string rootName,double beta){
  stringstream s;
  s<<rootName<<"beta"<<beta<<".dat";
  return s.str();
}

const string resultsFileName(const string rootName,double beta){
  stringstream s;
  s<<rootName<<"_beta"<<beta;
  return s.str();
}


// A dirty trick: copy Contractor methods here!
void contractL(mwArray& result,const mwArray& tmpL,
	       const Site& ket,const Site& bra,const mwArray& oper){
  int d=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int betab=tmpL.getDimension(0);int betak=tmpL.getDimension(1);
// #ifdef CHECKDIMS
//   if((Dl!=betab)||(Dlp!=betak)||(d!=dp)){
//     cout<<"Error in Contractor::contractL dimensions: bra="<<bra.getDimensions()
// 	<<", ket="<<ket.getDimensions()<<", tmpL="<<tmpL.getDimensions()<<endl;
//     exit(212);
//   }
// #endif
  result=ket.getA();
  result.reshape(Indices(dp,Dlp*Drp));
  result.multiplyLeft(oper);
  result.reshape(Indices(d,Dlp,Drp));
  result.permute(Indices(2,1,3));
  result.reshape(Indices(Dlp,d*Drp));
  result.multiplyLeft(tmpL);
  result.reshape(Indices(betab*d,Drp));
  mwArray aux=bra.getA();
  aux.permute(Indices(3,2,1));
  aux.conjugate();
  aux.reshape(Indices(Dr,Dl*d));
  result.multiplyLeft(aux);
}

void contractR(mwArray& result,const mwArray& tmpR,
	       const Site& ket,const Site& bra,const mwArray& oper){
  int d=bra.getd();int Dl=bra.getDl();int Dr=bra.getDr();
  int dp=ket.getd();int Dlp=ket.getDl();int Drp=ket.getDr();
  int betab=tmpR.getDimension(0);int betak=tmpR.getDimension(1);
// #ifdef CHECKDIMS
//   if((Dr!=betab)||(Drp!=betak)||(d!=dp)){
//     cout<<"Error in Contractor::contractR dimensions: bra="<<bra
// 	<<", ket="<<ket<<", tmpR="<<tmpR<<endl;
//     exit(212);
//   }
// #endif
  result=bra.getA();result.conjugate();
  result.reshape(Indices(d,Dl*Dr));
  result.multiplyLeft(permute(oper,Indices(2,1)));
  result.reshape(Indices(dp,Dl,Dr));
  result.permute(Indices(2,1,3));
  result.reshape(Indices(Dl*dp,Dr));
  result.multiplyRight(tmpR);
  result.reshape(Indices(Dl,d*betak));
  result.multiplyRight(reshape(permute(ket.getA(),Indices(1,3,2)),Indices(dp*Drp,Dlp)));
}

void computeCorrelators(const MPS& thS,int minX,int maxX,double beta,complex_t energy,ofstream* out){
  // The efficient version
  // Single site ops have to be doubled with id on the other side
  mwArray ZI,PlusI,MinI; // sigz, sigPlus,sigMinus
  constructOperatorProduct(ZI,sigZ,sig0);
  constructOperatorProduct(PlusI,.5*(sigX+I_c*sigY),sig0);
  constructOperatorProduct(MinI,.5*(sigX-I_c*sigY),sig0);
  
  Contractor& contractor=Contractor::theContractor();
  MPS aux(thS);aux.gaugeCond('L',1); // normalized gauge cond to left
  mwArray leftId=contractor.contract(aux,aux,minX,'L');
  //  mwArray rightId=contractor.contract(aux,aux,maxX,'R'); //(should be Id)
  for(int pos1=minX;pos1<maxX;pos1++){
    cout<<"pos1="<<pos1<<endl;
    mwArray tmpZ1,tmpId,tmpZz,tmpPz;
    contractL(tmpZ1,leftId,aux.getA(pos1),aux.getA(pos1),ZI);
    contractL(tmpId,leftId,aux.getA(pos1),aux.getA(pos1),identityMatrix(d*d));
    tmpZz=tmpZ1;
    contractL(tmpPz,leftId,aux.getA(pos1),aux.getA(pos1),PlusI);
    complex_t valZ1=tmpZ1.trace();
    // update leftId with one more Id before moving pos1
    leftId=tmpId;
    for(int pos2=pos1+1;pos2<=maxX;pos2++){
      mwArray tmpZ2;
      contractL(tmpZ2,tmpId,aux.getA(pos2),aux.getA(pos2),ZI);
      complex_t valZ2=tmpZ2.trace();

      contractL(tmpZ2,tmpZ1,aux.getA(pos2),aux.getA(pos2),ZI);
      complex_t valZ1Z2=tmpZ2.trace();
      tmpZ2=tmpZz;
      contractL(tmpZz,tmpZ2,aux.getA(pos2),aux.getA(pos2),ZI);
      complex_t valZstr=tmpZz.trace();
      
      contractL(tmpZ2,tmpPz,aux.getA(pos2),aux.getA(pos2),MinI);
      complex_t valPzM=tmpZ2.trace();

      // write to file
      *out<<beta<<"\t"<<real(energy)<<"\t"<<pos1<<"\t"<<pos2<<"\t";
      *out<<real(valZ1)<<"\t"; //<<imag(valZ1)<<"\t";
      *out<<real(valZ2)<<"\t"; //<<imag(valZ2)<<"\t";
      *out<<real(valZ1Z2)<<"\t"; //<<imag(valZ1Z2)<<"\t";
      *out<<real(valZstr)<<"\t"; //<<imag(valZstr)<<"\t";
      *out<<real(valPzM)<<"\t"<<imag(valPzM)<<"\t";
      *out<<endl;
      // now update the terms, except Zz, which is already done
      if(pos2<maxX){
	tmpZ2=tmpZ1;
	contractL(tmpZ1,tmpZ2,aux.getA(pos2),aux.getA(pos2),identityMatrix(d*d));
	tmpZ2=tmpId;
	contractL(tmpId,tmpZ2,aux.getA(pos2),aux.getA(pos2),identityMatrix(d*d));
	tmpZ2=tmpPz;
	contractL(tmpPz,tmpZ2,aux.getA(pos2),aux.getA(pos2),ZI);
      }
    }
  }
  
  
}
