#include <math.h>
#include <iomanip>

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "JoinedOperator.h"
#include "Properties.h"
#include "IsingHamiltonian.h"
#include <strstream>

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

#define MAXLEN 600

using namespace shrt;
using namespace std;


// Construct the MPO for the double commutator, directly.
void constructDoubleCommutatorIsing(MPO& mpo,int L,double g,double h);

void constructOperatorRange(MPS& opX,int range,const vector<complex_t>& coefs,double norm,complex_t phase,double cutoff);

void modulatedMPS(const MPS& opX,MPS& mps,int M);
int d=2;

/** 
    Take the result for  acertain range mdulated operator as produced
    by slowMPOhamilAnalysisHnOpt (program slowMPOHn, I believe) and
    compute the commutator for a large range of
    sies, to check the scaling. Since I do not use an exact
    modulation, I will interpolate the coefficients
*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  //  string directory="";
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

  int M1=props.getIntProperty("initM");
  int M2=props.getIntProperty("finM");
  int stepM=props.getIntProperty("stepM");

  int range=props.getIntProperty("range");
  // txt file containing the coefficients
  string inputfile=props.getProperty("inputfile");
  // the length from which to take the 
  //  int refM=props.getIntProperty("refM");

  string outfile=props.getProperty("outputfile");
  double gP=props.getDoubleProperty("gP");
  double hP=props.getDoubleProperty("hP");
  // minimum abs value of coefficient to be taken into account 
  double cutoff=props.getDoubleProperty("cutoff");

  //  string outfileIm=outfile+"_im";

  //  contractor.setConvTol(tol);
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  //  srandom(rSeed);

  // 1) Read the basic operator computed by the other program
  // In this case, it is  a list of terms of Hamiltonian powers, each
  // with a coefficient, so I have to construct the lists
  vector<vector<int> > indices;
  vector<complex_t> coeffs;
  ifstream input(inputfile.data());
  string aLine;
  double maxAbs=0.; //to keep track of the relative weight of terms
		    //and apply a cutoff
  complex_t phase=0.; // idem for phase
  while(input.good()){
    getline(input,aLine);
    //cout<<"Read line "<<aLine<<endl;
    int _n;
    vector<int> _code; // string of positions
    // if the proper range appears, take the line
    istrstream str(aLine.data());
    str>>_n;
    for(int k=0;k<_n;k++){
      // read one character
      char dig;
      str>>dig;
      _code.push_back(dig-'0');      
    }
    long double _re,_im;
    str>>_re>>_im;
    complex_t component={_re,_im};
    if(abs(component)>maxAbs){
      maxAbs=abs(component);
      phase=component/maxAbs;
    }
    coeffs.push_back(component);
    indices.push_back(_code);
    cout<<"Read line n="<<_n<<", indices="<<_code<<", value="<<component<<endl;
  }

  // Now I need to 

  exit(1);
  /**
  //  cout<<"Read from file vector of coeffs "<<coefs<<endl;
  // ********
  // Now write the op as a MPO
  MPS opX(range,1,d*d);
  constructOperatorRange(opX,range,coefs,norm,phase,cutoff);
  opX.gaugeCond('R',true);
  opX.gaugeCond('L',true);

  cout<<"Got the MPS for Opx: "<<opX<<endl;
  for(int m=M1;m<=M2;m+=stepM){
    MPO H2(m);
    constructDoubleCommutatorIsing(H2,m,gP,hP);

    //    cout<<"Got the double commutator for M="<<m<<endl;

    //For each size, I construct the MPS for the modulated operator
    MPS opM(m,1,d*d);
    modulatedMPS(opX,opM,m);
    opM.gaugeCond('R',1);
    //  cout<<"Constructed modulated MPS"<<opM<<endl;
    // and compute the lambda
    complex_t valH2=contractor.contract(opM,H2,opM);
    cout<<"valH2="<<valH2<<endl;
    // and the norm
    complex_t valN2=contractor.contract(opM,opM);
    cout<<"valN2="<<valN2<<endl;
    // and the trace
    MPS mpsId(m,1,d*d);
    {  
      MPO idMPO(m);
      mwArray id2=identityMatrix(2);
      for(int k=0;k<m;k++){
	idMPO.setOp(k,new Operator(reshape(id2,Indices(2,1,2,1))),true);
      }
      MPSfromMPO(idMPO,mpsId,d);
      mpsId.gaugeCond('R',true); // normalize properly
    }
    complex_t valT=contractor.contract(opM,mpsId);
    cout<<"valT="<<valT<<endl;
    complex_t lambda=valH2/(valN2-valT*conjugate(valT));
    ofstream* out;
    out=new ofstream(outfile.data(),ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfile<<
	" for output"<<endl;
      exit(1);
    }
    *out<<setprecision(10);
    *out<<m<<"\t"<<range<<"\t"<<real(lambda)<<endl;
    cout<<m<<"\t"<<range<<"\t"<<lambda<<endl;

    out->close();
    delete out;

  }

  */

}



double modulatedFunction(int M,int R,int pos){
    return cos(M_PIl*(pos+1.-(M-R+1)*.5)/(M-R));
}

void modulatedMPS(const MPS& opX,MPS& mps,int M){
  // I just need to expand the MPS I had 
  // Since I will have the sum of displaced MPS, I need a bond
  // dimension as big as the sum of bond dimensions of opX onall
  // links, with two extra values for D (0 and End). 
  int range=opX.getLength();
  int newD(0);
  for(int k=0;k<range;k++)
    newD+=opX.getA(k).getDr();
  mps=MPS(M,newD+1,d*d); // will just replace
  cout<<"The new dimensions of A matrices would be "<<newD<<"+1"<<endl;
  // the common part
  mwArray auxA(Indices(d*d,newD+1,newD+1));auxA.fillWithZero();
  // the block that I need to modulate
  mwArray auxFirst(Indices(d*d,newD+1,newD+1));auxFirst.fillWithZero();
  int offL=0;int offR=1; //for column 0 is for 0,0 el (Id)
  for(int p=0;p<range;p++){
    const mwArray& tmpA=opX.getA(p).getA();
    cout<<"Copying elements of A["<<p<<"]-"<<tmpA.getDimensions()<<endl;
    int dimL=tmpA.getDimension(1);
    int dimR=tmpA.getDimension(2);
    for(int i1=0;i1<d*d;i1++){
      for(int Dl=0;Dl<dimL;Dl++){
	for(int Dr=0;Dr<dimR;Dr++){
	  // cout<<"Setting element "<<Indices(i1,Dl+offL,Dr+offR)<<" from element "
	  //     <<Indices(i1,Dl+offL,Dr+offR)<<" of tensor in pos "<<p<<endl;
	  if(p>0){
	    auxA.setElement(tmpA.getElement(Indices(i1,Dl,Dr)),Indices(i1,Dl+offL,Dr+offR));
	  }
	  else{
	    auxFirst.setElement(tmpA.getElement(Indices(i1,Dl,Dr)),Indices(i1,Dl+offL,Dr+offR));
	  }
	}
      }
    }
    offL+=dimL;offR+=dimR;
    //    cout<<"After filling in tensor A["<<p<<"], offL="<<offL<<", offR="<<offR<<endl;
  }
  // and at the end,set the identity elements
  mwArray id=identityMatrix(d);id.reshape(Indices(d*d,1,1));
  for(int i1=0;i1<d*d;i1++){
    auxA.setElement(id.getElement(Indices(i1,0,0)),Indices(i1,0,0));
    auxA.setElement(id.getElement(Indices(i1,0,0)),Indices(i1,newD,newD));
  }  
  // the first will be the first row and the last the last column
  mwArray auxAL=auxA.subArray(Indices(-1,0,-1));auxAL.reshape(Indices(d*d,1,newD+1));
  mwArray auxFirstL=auxFirst.subArray(Indices(-1,0,-1));auxFirstL.reshape(Indices(d*d,1,newD+1));
  cout<<"Cut auxA(First) for L to dims "<<auxAL.getDimensions()<<endl;
  mwArray auxAR=auxA.subArray(Indices(-1,-1,newD));
  mwArray auxFirstR=auxFirst.subArray(Indices(-1,-1,newD));

  // and I need to introduce the modulation by multiplying a block
  for(int k=0;k<M;k++){
    mwArray newA;
    double cosFac=modulatedFunction(M,range-1,k);
    if(k==0)
      newA=auxAL+cosFac*auxFirstL;
    else
      if(k==M-1)
	newA=auxAR+cosFac*auxFirstR;
      else
	newA=auxA+cosFac*auxFirst;
    mps.replaceSite(k,newA,false);
  }
}

void sigmaProduct(const vector<int>& indices,mwArray& result){
  mwArray sig0=identityMatrix(d);sig0.reshape(Indices(d,1,d,1));
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);sigX.reshape(Indices(d,1,d,1));
  complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);sigY.reshape(Indices(d,1,d,1));
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);sigZ.reshape(Indices(d,1,d,1));
  Operator Op0(sig0),OpX(sigX),OpY(sigY),OpZ(sigZ);
  Operator* C[4]={&Op0,&OpX,&OpY,&OpZ};
  // I set a MPO with range sites and D=1
  int range=indices.size();
  MPO auxMPO(range);
  for(int ir=0;ir<range;ir++){
    auxMPO.setOp(ir,C[indices[ir]],false);
  }  
  // And expand it in result
  //  expandOper(auxMPO,result);
  //cout<<"MPO for "<<indices<<" expanded to "<<result<<endl;
  // Now I map it to MPS (just for convenience)
  MPS auxMPS(range,1,d*d);
  MPSfromMPO(auxMPO,auxMPS,true);
  expandVec(auxMPS,result);
}

// Return false in case it is at the end (all max - M-1)
// I increment WITHOUT the constraint that i1<=i2<=i3<=...iN
bool incrementRecursiveAll(vector<int>& indices,int M,int& ptr){
  // cout<<"incrementRecursive("<<indices<<") pos "<<ptr<<endl;
  if(indices[ptr]<M){
    indices[ptr]++;
    for(int k=ptr+1;k<indices.size();k++){
      indices[k]=0;
    }
    return true;
  }
  else{ // this one cannot be incremented => go left
    if(ptr>0)
      return incrementRecursiveAll(indices,M,--ptr);
    else
      return false;
  }
}

// Construct the range-body operator with the found coefficients and
// return it as a reference MPO
//void constructOperatorRange(mwArray& opX,int range,const vector<complex_t>& coefs,double norm,complex_t phase,double cutoff){
void constructOperatorRange(MPS& opX,int range,const vector<complex_t>& coefs,double norm,complex_t phase,double cutoff){
  mwArray operX;
  // add proper contributions one by one
  vector<int> indices(range,0); 
  bool done=false;
  int k=0;
  while(!done){
    if(abs(coefs[k])/norm>cutoff){
      // construct the full operator for this sequence, and add with
      // correct weight
      mwArray result;
      sigmaProduct(indices,result);
      if(operX.getDimension(0)==0) operX=((coefs[k]/phase)/norm)*result;
      else
	operX=operX+((coefs[k]/phase)/norm)*result;
    } // otherwise, ignore
    // increment index
    k++;
    int ptr=range-1;
    done=!incrementRecursiveAll(indices,3,ptr);
  }
  //cout<<"At the end, opX is "<<operX<<endl;
  // what if it is not Hermitian? Because coeffs could be complex!
  // reshape to MPS/MPO
  Indices dims;for(int p=0;p<range;p++)dims.push_back(d*d);
  //cout<<"dims="<<dims<<endl;
  opX=MPS(range,pow(d*d,floor(range/2)),d*d);
  opX.gaugeCond('R',1);  opX.gaugeCond('L',1);
  cout<<"Prepared MPS for oper: "<<opX<<endl;
  vecToMPS(operX,dims,opX,pow(d*d,floor(range/2)));
  cout<<"Computed MPS for oper: "<<opX<<endl;
}


// Copied from Liouvillian, should be somewhere else
void constructOperatorProduct(mwArray& result,const mwArray& opA,
			      const mwArray& opB){
#ifdef CHECKDIMS
  if(opA.getDimensions()!=opB.getDimensions()){
    cout<<"Error in LiouvillianXYEdges::constructOperatorProduct for"
	<<" A="<<opA<<" and B="<<opB<<endl;
    exit(1);
  }
#endif
  result=opA; 
  result.reshape(Indices(d*d,1));
  result.multiplyRight(reshape(opB,Indices(1,d*d)));
  result.reshape(Indices(d,d,d,d)); // order here: iji'j'
  result.permute(Indices(1,3,2,4));  // order in the MPO: ii',jj'  
  result.reshape(Indices(d*d,d*d));
}


void constructDoubleCommutatorIsing(MPO& mpo,int L,double g,double h){
  cout<<"constructDoubleCommutatorIsing(g="<<g<<",h="<<h<<")"<<endl;
  if(mpo.getLength()!=L)
    mpo.initLength(L);
  // Basic operators
  int nrOps=5;
  int dd=d*d;
  mwArray Z=mwArray(Indices(nrOps,dd,dd));
   // Basic pieces
  mwArray sig0=identityMatrix(d);//identity
  complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),datax);//sigmax
  complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  mwArray sigZ(Indices(d,d),dataz);//sigmaz

  // *** First of all, we need the identity on both
  mwArray term0=identityMatrix(dd);
  // *** Now the ones appearing in single-body terms
  // (1) sigX (x) Id
  mwArray term1;
  constructOperatorProduct(term1,sigX,sig0);
  // (2)  Id (x) sigX -> just a permutation of the previous one
  mwArray term2;
  constructOperatorProduct(term2,sig0,permute(sigX,Indices(2,1)));
  // (3) sigZ (x) Id
  mwArray term3;
  constructOperatorProduct(term3,sigZ,sig0);
  // (4)  Id (x) sigZ^T -> just a permutation of the previous one
  mwArray term4;
  constructOperatorProduct(term4,sig0,permute(sigZ,Indices(2,1)));
  // Fill in the operators in the Z array, to be contracted with C in
  // the construction of the MPO. The order is as listed above.
  for(int d1=0;d1<dd;d1++)
    for(int d2=0;d2<dd;d2++){
      Z.setElement(term0.getElement(Indices(d1,d2)),Indices(0,d1,d2));
      Z.setElement(term1.getElement(Indices(d1,d2)),Indices(1,d1,d2));
      Z.setElement(term2.getElement(Indices(d1,d2)),Indices(2,d1,d2));
      Z.setElement(term3.getElement(Indices(d1,d2)),Indices(3,d1,d2));
      Z.setElement(term4.getElement(Indices(d1,d2)),Indices(4,d1,d2));
    }
  Z.reshape(Indices(nrOps,dd*dd));
  //cout<<"Finished Z:"<<Z.getDimensions()<<endl;
  // now the coefficients
  int D=4; // bond dimension of the MPO
  mwArray C1(Indices(1,D,nrOps)); // the first site
  mwArray C(Indices(D,D,nrOps)); // the middle site
  mwArray CN(Indices(D,1,nrOps)); // the last site

  // Identity terms when nothing has happened yet
  C.setElement(ONE_c,Indices(0,0,0));
  C1.setElement(ONE_c,Indices(0,0,0));
  // Identity terms after the end
  C.setElement(ONE_c,Indices(D-1,D-1,0));
  CN.setElement(ONE_c,Indices(D-1,0,0));
  // Single body terms (one for each operator, with proper weights)
  C.setElement(h*ONE_c,Indices(0,D-1,3)); // h sigZ x I
  C.setElement(-h*ONE_c,Indices(0,D-1,4)); // -h I x sigZ^T
  C.setElement(g*ONE_c,Indices(0,D-1,1)); // g sigX x I
  C.setElement(-g*ONE_c,Indices(0,D-1,2)); // -g I x sigX^T
  C1.setElement(h*ONE_c,Indices(0,D-1,3)); // h sigZ x I
  C1.setElement(-h*ONE_c,Indices(0,D-1,4)); // -h I x sigZ^T
  C1.setElement(g*ONE_c,Indices(0,D-1,1)); // g sigX x I
  C1.setElement(-g*ONE_c,Indices(0,D-1,2)); // -g I x sigX^T
  CN.setElement(h*ONE_c,Indices(0,0,3)); // h sigZ x I
  CN.setElement(-h*ONE_c,Indices(0,0,4)); // -h I x sigZ^T
  CN.setElement(g*ONE_c,Indices(0,0,1)); // g sigX x I
  CN.setElement(-g*ONE_c,Indices(0,0,2)); // -g I x sigX^T
  //  Two body terms from H and H^T on the other side
  C.setElement(ONE_c,Indices(0,1,3)); //  sigZ x I
  C.setElement(ONE_c,Indices(1,D-1,3)); 
  C.setElement(-1.*ONE_c,Indices(0,2,4));  //  I x sigZ 
  C.setElement(ONE_c,Indices(2,D-1,4));
  C1.setElement(ONE_c,Indices(0,1,3)); //  sigZ x I
  C1.setElement(-1.*ONE_c,Indices(0,2,4));  //  I x sigZ 
  CN.setElement(ONE_c,Indices(1,0,3)); 
  CN.setElement(ONE_c,Indices(2,0,4));

  // Now set the operators in place
  // Now I reshape, contract with Z and set the indices in proper order
  C.reshape(Indices(D*D,nrOps));
  C1.reshape(Indices(1*D,nrOps));
  CN.reshape(Indices(D*1,nrOps));
  C.multiplyRight(Z);C.reshape(Indices(D,D,dd,dd));C.permute(Indices(3,1,4,2));

  C1.multiplyRight(Z);C1.reshape(Indices(1*D*dd,dd));
  // Special situation: first and last sites are contracted with Id as traced out
  mwArray id2=1/sqrt(2.)*identityMatrix(d);id2.reshape(Indices(dd,1));
  C1.multiplyRight(id2);
  C1.reshape(Indices(1,D,dd,1));C1.permute(Indices(3,1,4,2));
  CN.multiplyRight(Z);CN.reshape(Indices(D*1*dd,dd));
  CN.multiplyRight(id2);
  CN.reshape(Indices(D,1,dd,1));CN.permute(Indices(3,1,4,2)); // dd x D x 1 x 1

  // Now I build the MPO with required ops.
  //1) on the edges, I contract C1 and CN with their adjoint
  C1.reshape(Indices(dd,D));
  CN.reshape(Indices(dd,D)); // dd D
  mwArray tmp(C1);tmp.Hconjugate(); 
  C1.multiplyLeft(tmp); // Dl,Dr
  tmp=CN;tmp.Hconjugate(); 
  CN.multiplyLeft(tmp); //Dlu,Dld
  // 2) These, I contract with the middle ones 
  tmp=C;tmp.permute(Indices(2,1,3,4));tmp.reshape(Indices(D,dd*dd*D));
  mwArray edgeL=C1*tmp;edgeL.reshape(Indices(D*dd,1,dd,D)); // Ddd 1 dd D
  tmp=C;tmp.reshape(Indices(dd*D*dd,D));
  CN.permute(Indices(2,1)); // put first the index for 2nd MPO
  mwArray edgeR=tmp*CN;edgeR.reshape(Indices(dd,D,dd,D));
  edgeR.permute(Indices(1,4,2,3));edgeR.reshape(Indices(dd*D,D,dd,1)); // ddD D dd 1
  // Now for the adjoint I also need modifications
  mwArray Cd(C);Cd.permute(Indices(3,2,1,4),true); // adjoint of normal C
  mwArray edgeLd(Cd);edgeLd.reshape(Indices(dd,1,D*dd,D));
  mwArray edgeRd(Cd);edgeRd.reshape(Indices(dd,D,dd*D,1));
  Operator CdOp(Cd),COp(C);
  Operator edgeLOp(edgeL),edgeROp(edgeR);
  Operator edgeLdOp(edgeLd),edgeRdOp(edgeRd);
  const Operator* arrL[2]={&edgeLdOp,&edgeLOp};
  const Operator* arrR[2]={&edgeRdOp,&edgeROp};
  const Operator* arrC[2]={&CdOp,&COp};
  //cout<<"About to place operators in MPO"<<endl;
  // And now I set the operators one by one
  mpo.setOp(0,new JoinedOperator(2,arrL),true);
  mpo.setOp(L-1,new JoinedOperator(2,arrR),true);
  if(L>2){
    mpo.setOp(1,new JoinedOperator(2,arrC),true);
    for(int k=2;k<L-1;k++)
      mpo.setOp(k,&mpo.getOp(1),false);
  }
}
