#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

// Make the random values discrete (-1 or 1)
#define DISCRETE 0

#include "HeisenbergHamiltonian.h"

using namespace shrt;

/** 
    Takes the MPS files produced by locaPT2 (localizationPTcont) at different time steps, and 
    compute the distance between the rho and the corresponding truncated one, when the truncation 
    is done exactly in the middle, in terms of trace norm. 
    I can also compute the difference in <M1> from the full rho to its truncations.

    Receives arguments:
    \param <L> (int) number of sites
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian (Jx=Jy=Jz)
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian (local values 
                        randomly distributed between -h and h)			
    \param <nrInst> (int) Nr of row in the file to be used (if rows
                       are longer than L, the first L values of the row will be used).
    \param <delta> (double) width of the time step used
    \param <Mmin> (int) min number of steps
    \param <Mmax> (int) max number of steps
    \param <deltaM> (int) interval in number of steps
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the results
    \param <mpsfile> (char*) basename of the files where to read the MPS 
*/



/** Construct the initial state and observable, but vectorized*/
void constructM1(int L,int d,MPS& result);

/** Vectorize an MPO to make an MPS out of it. If up==true (default)
    the first component of the double physical index will be the one
    which was going up.*/
void MPSfromMPO(const MPO& mpo,MPS& mps,bool up=true);

/** Get the whole list of MPS files with common basename for time steps from Mmin to Mmax, in steps of DeltaM */
void getListFiles(const string& basename,vector<string>& results,vector<int>& cnts,int Mmin,int Mmax,int DeltaM=100);

/** A function that prepared the Schmidt decomposition of an MPS
    across the cut after l sites. It should move to contractor or
    MPS, once it is properly tested! */
void schmidt(MPS& mps,int pos,vector<double>& lambdas);
/** Dirty trick to finish the schmidt decomp */
void siteGaugeL(mwArray& A,mwArray& leftterm);

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double J_=atof(argv[++cntr]);
  double h_=atof(argv[++cntr]);
  int nrInst=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int Mmin=atoi(argv[++cntr]);
  int Mmax=atoi(argv[++cntr]);
  int DeltaM=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  const char* mpsfile=argv[++cntr];

  int d=2;
  bool initFile=!file_exists(outfname);
  ofstream* out;
  out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  if(initFile){
    *out<<"% L="<<L<<", J="<<J_<<", h="<<h_
	<<" D="<<D<<", nrInst="<<nrInst<<endl;
    *out<<"% MPS in "<<mpsfile<<endl;
    *out<<"% cnt\t delta\t time\t k(nr SchVal)\t lambda(k)\t norm1L\t norm1R\t norm2L\t norm2R\t M1 (Re,Im)\t M1_truncK(Re,Im)\t overlap(truncMid) (Re,Im)\t M1_truncK(Re,Im)"<<endl;
  }

  MPS M1(L,2,d*d);
  constructM1(L,d,M1);

  // List of files to process
  vector<string> files;
  vector<int> cnts;
  getListFiles(mpsfile,files,cnts,Mmin,Mmax,DeltaM);

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;

  // I need to iterate reading the MPS file, computing 
  // 1) EV of M1
  // 2) truncation(s) of the MPS to smaller Ds and their M1s
  // 3) In particular, ALL left and right MPOs when truncating in the center, their trace 
  // norms (also operator norms, as I am at it)

  vector<string>::iterator it=files.begin();
  int ptr=0;
  while(it!=files.end()){
    cout<<"Analyzing state in file "<<*it<<", cnt="<<cnts[ptr]<<endl;
    MPS rho(L,1,1);
    rho.importMPS(it->data());
    rho.gaugeCond('R',1);
    complex_t M1ev=contractor.contract(rho,M1);

    vector<double> SchVal; // Schmidt values of rho seen as MPS (center)
    //    contractor.getSchmidtValues(rho,SchVal,L/2);

    // Now, for every term I compute the norms of left and right
    // "Schmidt vectors" (the normalization is not correct, at least
    // for the left ones, but I could recover from that knowing the
    // lambdas) 
    MPS auxRho(rho);
    schmidt(auxRho,L/2-1,SchVal);    

    cout<<"Found Schmidt decomposition with lambdas "<<SchVal<<endl;

    // Construct two MPOs for left and right part, with an open bond
    // at the end (the right one is constructed from right to left to
    // make it easier for expandOper)
    MPO lMPO(L/2),rMPO(L/2);
    for(int pos=0;pos<L;pos++){
      mwArray aux=auxRho.getA(pos).getA();
      int Dl=aux.getDimension(1);
      int Dr=aux.getDimension(2);
      aux.reshape(Indices(d,d,Dl,Dr));
      if(pos<L/2)
	lMPO.setRotatedOp(pos,aux,Indices(1,3,2,4));
      else // right ones are reflected
	rMPO.setRotatedOp(L-1-pos,aux,Indices(1,4,2,3));
    }
    // Could expand now, but will give me Dx(2^10)^2 big objects
    mwArray OpL,OpR;
    expandOper(lMPO,OpL);
    cout<<"Expanded left term to dims "<<OpL.getDimensions()<<endl;
    expandOper(rMPO,OpR);
    cout<<"Expanded right term to dims "<<OpR.getDimensions()<<endl;
    // They will have open bond dim D on the right each
    for(int k=0;k<D;k++){
      cout<<"*** k="<<k<<endl;
      // Trick: Cut the edges
      mwArray auxL=OpL.subArray(Indices(-1,-1,k));
      cout<<"Cut left term to dims "<<auxL.getDimensions()<<endl;
      mwArray auxR=OpR.subArray(Indices(-1,-1,k));
      cout<<"Cut right term to dims "<<auxR.getDimensions()<<endl;
      // And compute norms, etc
      // trace norm, is the sum of singular values
      mwArray U,S,Vdag;
      wrapper::svd(auxL,U,S,Vdag);
      double norm1L=real(S.trace());
      // the operator norm is the max sing value (first one,a s they come in order)
      double norm2L=real(S.getElement(Indices(0,0)));
      cout<<"Found SVD of left term, so norm1="<<norm1L<<", norm2="<<norm2L<<endl;
      wrapper::svd(auxR,U,S,Vdag);
      double norm1R=real(S.trace());
      double norm2R=real(S.getElement(Indices(0,0)));
      cout<<"Found SVD of right term, so norm1="<<norm1R<<", norm2="<<norm2R<<endl;
      // Now write down results
      *out<<cnts[ptr]<<"\t"<<delta<<"\t"<<cnts[ptr]*delta<<"\t";
      *out<<k+1<<"\t"<<SchVal[k]<<"\t";
      *out<<norm1L<<"\t"<<norm1R<<"\t";
      *out<<norm2L<<"\t"<<norm2R<<"\t";
      // and also the <M1>, overlap and norm of rho truncated only in the
      // middle bond
      mwArray projK(Indices(D,D));
      for(int l=0;l<D-k;l++) projK.setElement(SchVal[k]*ONE_c,Indices(l,l));
      MPS mpsTruncMid(auxRho);
      auxL=mpsTruncMid.getA(L/2-1).getA();
      Indices dims=auxL.getDimensions();
      auxL.reshape(Indices(-1,dims[1]));
      auxL.multiplyRight(projK);
      auxL.reshape(dims);
      mpsTruncMid.setA(L/2-1,auxL);
      // Now compute overlap and value of M1
      complex_t M1evMidT=contractor.contract(mpsTruncMid,M1);
      complex_t normMidT=contractor.contract(mpsTruncMid,mpsTruncMid);
      complex_t overl=contractor.contract(mpsTruncMid,rho);
      *out<<real(M1ev)<<"\t"<<imag(M1ev)<<"\t";
      *out<<real(M1evMidT)<<"\t"<<imag(M1evMidT)<<"\t";
      *out<<real(overl)<<"\t"<<imag(overl)<<"\t";
      *out<<real(normMidT)<<"\t"<<imag(normMidT)<<"\t";

      // also if I truncate the whole MPS at (D-k), what is the error in M1
      MPS trunc(rho,D-k);
      contractor.optimizeMPS(rho,trunc,D-k);
      complex_t M1evTr=contractor.contract(trunc,M1);
      *out<<real(M1evTr)<<"\t"<<imag(M1evTr)<<"\t";
      *out<<endl;
    }

    it++;ptr++;

    rho.clear();
    auxRho.clear();
    SchVal.clear();
    lMPO.clear();
    rMPO.clear();
  }





}


void constructM1(int L,int d,MPS& result){
   mwArray sig0=identityMatrix(d);
   sig0.reshape(Indices(d*d,1,1));
   complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
   mwArray sigZ(Indices(d*d,1,1),dataZ);
   mwArray Z(Indices(2,d*d));
   for(int j=0;j<d*d;j++){
     Z.setElement(sig0.getElement(Indices(j,0,0)),Indices(0,j));
     Z.setElement(sigZ.getElement(Indices(j,0,0)),Indices(1,j));
   }
   int D=2; // Dimension of the MPS
   mwArray C(Indices(D,D,2)),C1(Indices(1,D,2)),CL(Indices(D,1,2));
   C.setElement(ONE_c,Indices(0,0,0));
   C1.setElement(ONE_c,Indices(0,0,0));
   C.setElement(ONE_c,Indices(D-1,D-1,0));
   CL.setElement(ONE_c,Indices(D-1,0,0));
   // There is a site-dependent element C(0,D-1,1)
   C1.setElement(ONE_c,Indices(0,D-1,1));
   CL.setElement(exp((-2*M_PIl*(L-1.)/L)*I_c),Indices(0,0,1));
   C1.reshape(Indices(1*D,2));
   C1.multiplyRight(Z);
   C1.reshape(Indices(1,D,d*d));
   C1.permute(Indices(Indices(3,1,2)));
   result.setA(0,C1);
   CL.reshape(Indices(D*1,2));
   CL.multiplyRight(Z);
   CL.reshape(Indices(D,1,d*d));
   CL.permute(Indices(Indices(3,1,2)));
   result.setA(L-1,CL);
   for(int k=1;k<L-1;k++){
     C.setElement(exp((-2*M_PIl*k/L)*I_c),Indices(0,D-1,1));
     mwArray C_(C);
     C_.reshape(Indices(D*D,2));
     C_.multiplyRight(Z);
     C_.reshape(Indices(D,D,d*d));
     C_.permute(Indices(Indices(3,1,2)));
     result.setA(k,C_);
   }

}

void MPSfromMPO(const MPO& mpo,MPS& mps,bool up){
  int L=mpo.getLength();
  // vector<int> du=mpo.getDimensions();
  // vector<int> dd=mpo.getOriginalDimensions();
  // int D=mpo.getOp(0).getDr();

  // vector<int> newDims(du);
  // for(int k=0;k<L;k++) newDims[k]*=dd[k];

  mps=MPS(L,1,1);

  for(int k=0;k<L;k++){
    mwArray aux=mpo.getOp(k).getFullData();
    Indices dims=aux.getDimensions();
    if(up)
      aux.permute(Indices(1,3,2,4));
    else
      aux.permute(Indices(3,1,2,4));
    aux.reshape(Indices(dims[0]*dims[2],dims[1],dims[3]));
    mps.replaceSite(k,aux,0); // brute force
  }

}


void getListFiles(const string& basename,vector<string>& results,vector<int>& cnts,int Mmin,int Mmax,int DeltaM){
  int M=Mmin;
  while(M<=Mmax){
    // see if file with extension _M exisits and return name
    stringstream filename;
    filename<<basename<<"_"<<M;
    if(file_exists(filename.str())){
      results.push_back(filename.str());
      cnts.push_back(M);
      cout<<"Added "<<filename.str()<<" to the file list"<<endl;
    }
    M+=DeltaM;
  }
}

void schmidt(MPS& mps,int pos,vector<double>& lambdas){
  Contractor& contractor=Contractor::theContractor();
  // First, normalize, and impose gauge to R, but keep normalization
  // factor
  mps.gaugeCond('R',0);
  double norm=mps.getNormFact();
  mps.gaugeCond('R',1);
  int L=mps.getLength();
  // Now impose gauge to left, but only until site pos (pos-1 is
  // already on the left of the cut)
  for(int k=L-1;k>pos;k--){
    // except for the last one, I absorb the transformation in the next site
    mps.gaugeCond(k,'L',true);
  }
  // Now the next one has to be done by hand, because Site does not
  // let me use it from here
  mwArray A=mps.getA(pos).getA();
  mwArray Lamb;
  siteGaugeL(A,Lamb);
  // Now I need to diagonalize Lamb
  mwArray U,S,Vdag;
  wrapper::svd(Lamb,U,S,Vdag);
  A=mps.getA(pos-1).getA();
  int newD=U.getDimension(1); // new right dim of this tensor
  Indices dims=A.getDimensions();
  Indices packedDims(dims[0]*dims[1],dims[2]);
  A.reshape(packedDims);
  A=A*U;
  if(newD!=dims[2]) dims[2]=newD;
  A.reshape(dims);
  mps.setA(pos-1,A);
  // Noqw the right one
  A=mps.getA(pos).getA();
  newD=Vdag.getDimension(0);
  Indices neworder(2,1,3);
  dims=A.getDimensions(); // d Dl Dr
  packedDims=Indices(dims[1],dims[0]*dims[2]);
  A.permute(neworder);
  A.reshape(packedDims);
  A=Vdag*A; //newD dxDr
  if(newD!=dims[1]) dims[1]=newD;
  Indices newdims(dims[1],dims[0],dims[2]);
  A.reshape(newdims);
  A.permute(neworder);
  mps.setA(pos,A);
  // And the lambdas should be real
  for(int k=0;k<newD;k++){
    lambdas.push_back(real(S.getElement(Indices(k,k))));
  }
}

void siteGaugeL(mwArray& A,mwArray& leftterm){
  //Since A is dxDlxDr, I have to permute indices
  Indices neworder(2,1,3);
  A.permute(neworder);
  Indices dims=A.getDimensions(); // now Dl, d, Dr
  // And reshape to Dlx(d Dr)
  Indices packedDims(dims[0],dims[1]*dims[2]);
  wrapper::lq(reshape(A,packedDims),leftterm,A);
  int Dl=A.getDimension(0); // new Dl'
  dims[0]=Dl;
  A.reshape(dims); // Dl' d Dr
  A.permute(neworder);
  if(Dl==1){//fist site => do not include phase in the tensor
    complex_t aux=leftterm.getElement(0); // single element
    complex_t phase=(1./abs(aux))*aux;
    A.multiplyRight(mwArray(phase));
    leftterm.setElement((complex_t){abs(aux),0.},0);
  }
}

