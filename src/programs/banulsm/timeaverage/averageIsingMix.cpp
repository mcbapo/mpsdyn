#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "VariableOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

#include "IsingHamiltonian.h"

using namespace shrt;

/** averageIsing constructs (by hand) a transverse operator to try to 
    approximate the time average after long time, using
    Ising Hamiltonian \ref <IsingHamiltonian>.

    Receives arguments:
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <delta> (double) width of the time step used
    \param <M> (int) number of steps
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the energy
    \param <app> (int) whether the output file is to be kept (app==1)

    The first version has a problem. Although the MPO is right, it
    does not have a unique max. eigenvector. On the contrary, it has
    M+1 degenerate ones, all with lambda=1, corresponding to the
    transfer matrices of all the possible evolutions involved, and
    therefore the power method will not work. Could I try to estimate
    directly the trace?
    Or could I find two orthogonal eigenvectors?
    Here I try to initialize a symmetric state by buildInitialSymMPS
    And also with the symmetric proyector.

    In averageIsingMix, I forget about the separate subspaces and mix
    everything, to check what's the eigenvector of that.

*/

void swap(const MPS& mps,int pos,MPS& result);
void swapProj(const MPS& mps,int pos,MPS& result);
void prodProj(const MPS& mps,int pos,MPS& result);
void swapProj2(const MPS& mps,int pos,MPS& result,char dir);

void buildInitialSymMPS(MPS& mps,mwArray& dataState);
void buildSymMPO(MPO& mpo,int Dop,int len);
//void findMaxEigenvector(const MPO& ops,int D,
void findMaxEigenvector(const MPO& ops,const MPO& proj,int D,
			double* lambda,char dir,
			MPS& lastMPS,double convtol=1E-5);
complex_t computeAllProjections(int len,int pos,int Dop,MPS& rEV,MPS& lEV,
				MPO& oper,MPO& idMPO,ofstream* out);
void projOnSym(int len,int Dop,MPS& mps,char dir);

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  double J=atof(argv[++cntr]);
  double g=atof(argv[++cntr]);
  double h=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  if(D%2!=0){
    cout<<"Here I need even D (for the moment), increasing it in 1"<<endl;
    D++;
  }

  cout<<"Initialized arguments "
      <<", J="<<J
      <<", g="<<g
      <<", h="<<h
      <<", delta="<<delta
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% J="<<J<<", g="<<g<<", h="<<h
	<<" D="<<D<<endl;
    *out<<"% M\t delta\t t\t D\t avrg sigx"<<endl;
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }

   Contractor& contractor=Contractor::theContractor();
   cout<<"Initialized Contractor"<<endl;
   int d=2;

   // First: create the Hamiltonian
   int Dop=2; // bond dimension of Ising exponential mpo
   IsingHamiltonian hamI(3,d,J,g,h);
   cout<<"Created the Hamiltonian"<<endl;
   mwArray basicU=hamI.getUOperator(delta,0.,false,2);
   // This is the du+Dl+dd+Dr basic piece
   // I want to increase it with an identity
   // and combine with the adjoints for the ancillary side

   // TEST!!!!
   mwArray identOp=identityMatrix(d*d); // double identity
   mwArray basicUad(basicU);
   basicUad.conjugate(); // the adjoint for the ancillas

   mwArray doubleU(basicU);
   doubleU.reshape(Indices(-1,1));
   basicUad.reshape(Indices(1,-1));
   doubleU.multiplyRight(basicUad);
   doubleU.reshape(Indices(d,Dop,d,Dop,d,Dop,d,Dop));
   doubleU.permute(Indices(1,5,2,6,3,7,4,8));
   doubleU.reshape(Indices(d*d,Dop*Dop,d*d,Dop*Dop));

   // The coefficients of the "mpo"
   mwArray C(Indices(2,2,2));
   C.setElement(ONE_c,Indices(0,0,0));
   C.setElement(ONE_c,Indices(1,1,1));
   C.setElement(ONE_c,Indices(1,0,0));

   // The "operator basis"
   mwArray Z(Indices(2,d*d,Dop*Dop,d*d,Dop*Dop));
   for(int i=0;i<d*d;i++)
     for(int j=0;j<d*d;j++){
      for(int l=0;l<Dop*Dop;l++)
	 for(int r=0;r<Dop*Dop;r++){
	   if(i==j&&l==r)
	     Z.setElement(ONE_c,Indices(0,i,l,j,r));
 	   Z.setElement(doubleU.getElement(Indices(i,l,j,r)),
			Indices(1,i,l,j,r));
	 }
     }

   // Now construct the Operator itselt
   C.reshape(Indices(-1,2));
   Z.reshape(Indices(2,-1));
   mwArray data=C*Z;
   data.reshape(Indices(2,2,d*d,Dop*Dop,d*d,Dop*Dop));
   // this is as: duxDlxddxDr, now left is dd
   data.permute(Indices(2,3,4,1,5,6));
   data.reshape(Indices(2*d*d,Dop*Dop,2*d*d,Dop*Dop));
   cout<<"Created operator data "<<data.getDimensions()<<endl;
   data.permute(Indices(2,3,4,1)); // as transverse coming from right
   Operator* oper=new Operator(data);

   // now the edges (contracted with a middle op!)
   mwArray stateX(Indices(2,1));
   stateX.setElement(1./sqrt(2)*ONE_c,Indices(0,0));
   stateX.setElement(1./sqrt(2)*ONE_c,Indices(1,0));
   mwArray adjSt(Hconjugate(stateX));
   mwArray dataX=stateX*adjSt;
   mwArray edgeSt(Indices(2,d,d));
   for(int i=0;i<d;i++)
     for(int j=0;j<d;j++){
       //edgeSt.setElement(dataX.getElement(Indices(i,j)),Indices(0,i,j));
       edgeSt.setElement(dataX.getElement(Indices(i,j)),Indices(1,i,j));
     }
   edgeSt.reshape(Indices(1,2*d*d));
   cout<<"Created operator edge X "<<edgeSt.getDimensions()<<endl;
   edgeSt.multiplyRight(reshape(permute(data,Indices(2,1,3,4)),
				Indices(2*d*d,-1)));
   edgeSt.reshape(Indices(1,Dop*Dop,Dop*Dop,2*d*d));
   edgeSt.permute(Indices(2,1,3,4));
   cout<<"Created operator edge X "<<edgeSt.getDimensions()<<endl;
   Operator* operX=new Operator(edgeSt);

   // identity op
   mwArray lastId=identityMatrix(d);
   mwArray idData(Indices(2,d,d));
   mwArray dataAux(data); // to multiply with std oper, unless there is only 1
   if(M==1)
     dataAux=edgeSt;
   for(int i=0;i<d;i++)
     for(int j=0;j<d;j++){
       idData.setElement(lastId.getElement(Indices(i,j)),Indices(0,i,j));
       idData.setElement(lastId.getElement(Indices(i,j)),Indices(1,i,j));
     }
   idData.reshape(Indices(1,2*d*d,1,1));
   cout<<"Created operator edge Id "<<idData.getDimensions()<<endl;
   int Dl=dataAux.getDimension(0);
   int Dr=dataAux.getDimension(2);
   int dd=dataAux.getDimension(1);
   int du=dataAux.getDimension(3);
   idData.reshape(Indices(du,1));
   idData.multiplyLeft(reshape(dataAux,Indices(-1,du)));
   idData.reshape(Indices(Dl,dd,Dr,1));
   Operator* operId=new Operator(idData);

   // at the end op sigmaX (contracted with a middle op!
   complex_t sigmaX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
   mwArray sigX(Indices(d,d),sigmaX);
   mwArray opData(Indices(2,d,d));
   for(int i=0;i<d;i++)
     for(int j=0;j<d;j++){
       opData.setElement(sigX.getElement(Indices(i,j)),Indices(0,i,j));
       opData.setElement(sigX.getElement(Indices(i,j)),Indices(1,i,j));
     }
   opData.reshape(Indices(2*d*d,1));
   opData.multiplyLeft(reshape(dataAux,Indices(-1,du)));
   opData.reshape(Indices(Dl,dd,Dr,1));
   cout<<"Created operator edge opX "<<opData.getDimensions()<<endl;
   Operator* operOp=new Operator(opData);
   
   // An attempt to find a right/left eigenvector
   int len=M;
   MPO trMPO(len);
   MPO trMPOx(len);
   for(int k=1;k<len-1;k++){
     trMPO.setOp(k,oper);
     trMPOx.setOp(k,oper);
   }
   trMPO.setOp(0,operX);
   trMPO.setOp(len-1,operId);
   trMPOx.setOp(0,operX);
   trMPOx.setOp(len-1,operOp);
  // the operator with a sigx at the end

   if(M<=6){
     // to check, save the operator MPO as text
     trMPO.exportMPOtext("/home/banulsm/SVN/peps/projects/mpsdyn/src/bin/transv1.txt");
     trMPOx.exportMPOtext("/home/banulsm/SVN/peps/projects/mpsdyn/src/bin/transvX.txt");
   }

   int maxiter=1;
   MPS anEV(len,D,Dop*Dop);
   anEV.setRandomState();
   MPO proj(len);
   buildSymMPO(proj,Dop,len);
   // Make it "symmetric", somehow
   MPS aux(anEV);
   //   contractor.optimize(proj,aux,anEV);
   //projOnSym(len,Dop,anEV,'R');
   //projOnSym(len,Dop,aux,'L');

   std::vector<MPS> rightEV(maxiter,anEV);
   std::vector<MPS> leftEV(maxiter,aux);

   for(int iter=0;iter<maxiter;iter++){
     MPS& _rightEV(rightEV[iter]);
     MPS& _leftEV(leftEV[iter]);
     if(len>1&&0){
       buildInitialSymMPS(_rightEV,dataX);
       _rightEV.gaugeCond('R',1);
       buildInitialSymMPS(_leftEV,dataX);
       _leftEV.gaugeCond('R',1);
     }
     double lambda=0;
     cout<<"Starting the optimization"<<endl;
     contractor.findRightEigenvector(trMPO,D,&lambda,_rightEV);
     contractor.findLeftEigenvector(trMPO,D,&lambda,_leftEV);
     MPS::gaugeCond(_rightEV,_leftEV,'R');
   }

  if(M<=6){
     // to check, save the MPSs too as text
     rightEV[0].exportMPStext("/home/banulsm/SVN/peps/projects/mpsdyn/src/bin/rightEV.txt");
     leftEV[0].exportMPStext("/home/banulsm/SVN/peps/projects/mpsdyn/src/bin/leftEV.txt");
   }

   std::vector<complex_t> denom(maxiter);
   std::vector<complex_t> num(maxiter);
   for(int k=0;k<maxiter;k++){
     // should be just lambda, with the gauge above
     denom[k]=contractor.contract(rightEV[k],trMPO,leftEV[k]);
   }
   //complex_t denom=lambda*ONE_c;

   for(int k=0;k<maxiter;k++){
     num[k]=contractor.contract(rightEV[k],trMPOx,leftEV[k]);
   }

   double time=M*delta;

   cout<<"For total time "<<time<<", values= ";
   complex_t avrg=ZERO_c;
   for(int k=0;k<maxiter;k++){
     cout<<num[k]/denom[k]<<", ";
     avrg=avrg+num[k]/denom[k];
   }
   cout<<endl;


   *out<<"%"<<M<<"\t"<<delta<<"\t"<<time+delta<<"\t"<<D<<"\t"<<avrg/maxiter;
   for(int k=0;k<maxiter;k++)
     *out<<"\t"<<num[k]/denom[k];
   *out<<endl;

   // Now try the projections of the eigenvectors!
   MPS auxR(rightEV[0]);
   MPS auxL(leftEV[0]);
   MPO tmpMPO(len);
   MPO tmpMPOx(len);
   mwArray identD=identityMatrix(Dop*Dop);
   identD.reshape(Indices(Dop*Dop,1,Dop*Dop,1));
   Operator* idOp=new Operator(identD);
   for(int k=0;k<len;k++){
     tmpMPO.setOp(k,&trMPO.getOp(k));
     tmpMPOx.setOp(k,&trMPOx.getOp(k));
   }
   //   for(int pos=len-2;pos>=0;pos--){
   for(int pos=len-2;pos>0;pos--){
     cout<<"Computing \"projection\" for position "<<pos<<endl;
     *out<<pos<<"\t"<<pos*delta<<"\t";
     if(OPTIONA){ 
       // swap dos ultimos
       MPS aux1(auxR);
       swap(aux1,pos,auxR);
       MPS aux2(auxL);
       swap(aux2,pos,auxL);
     }
     // proyeccion
     MPS projR(auxR);
     MPS projL(auxL);
     if(OPTIONA){
       prodProj(auxR,pos,projR);
       prodProj(auxL,pos,projL);
       MPS::gaugeCond(projR,projL,'R');
     }     
     else{
       swapProj2(auxR,pos,projR,'R');
       swapProj2(auxL,pos,projL,'L');
     }

     tmpMPO.setOp(pos,&trMPO.getOp(len-1));
     tmpMPOx.setOp(pos,&trMPOx.getOp(len-1));
     tmpMPO.setOp(pos+1,idOp);
     tmpMPOx.setOp(pos+1,idOp);
     complex_t num=contractor.contract(projR,tmpMPOx,projL);
     complex_t den=contractor.contract(projR,tmpMPO,projL);
     complex_t val=num/den;
     *out<<val<<endl;
   }
   delete idOp;
   delete oper;
   delete operX;
   delete operId;
   delete operOp;
   out->close();
}

void swap(const MPS& mps,int pos,MPS& result){
  cout<<"Swaping "<<pos<<" and next"<<endl;
  int len=mps.getLength();
  if(pos<0||pos>=len-1){
    cout<<"Error: can only swap 0-"<<len-2<<" with next one"<<endl;
    exit(2);
  }
  int d=result.getA(pos).getd();
  static bool first=true;
  static mwArray swapL;
  static mwArray swapR;
  if(first){
    mwArray swap(Indices(d,d,d,d));
    for(int i=0;i<d;i++)
      for(int j=0;j<d;j++)
	swap.setElement(ONE_c,Indices(i,j,j,i));
    swap.permute(Indices(1,3,2,4));
    swap.reshape(Indices(d*d,d*d));
    mwArray S; // the diagonal values are all one
    wrapper::svd(swap,swapL,S,swapR);
    swapL.reshape(Indices(d,d,d*d));
    swapL.permute(Indices(1,3,2));
    swapL.reshape(Indices(d*d*d,d)); // up, right, down
    swapR.reshape(Indices(d*d*d,d)); // left, up, down
    
    first=false;
  }
  result=mps;
  mwArray left=result.getA(pos).getA();
  mwArray right=result.getA(pos+1).getA();
  int Dl=left.getDimension(1);
  int Dr=left.getDimension(2);
  left.reshape(Indices(d,-1));
  left.multiplyLeft(swapL);
  left.reshape(Indices(d,d*d,Dl,Dr));
  left.permute(Indices(1,3,2,4));
  left.reshape(Indices(d,Dl,d*d*Dr));
  int Drr=right.getDimension(2);
  right.reshape(Indices(d,-1));
  right.multiplyLeft(swapR);
  right.reshape(Indices(d*d,d,Dr,Drr));
  right.permute(Indices(2,1,3,4));
  right.reshape(Indices(d,d*d*Dr,Drr));
  result.replaceSite(pos,left,false);
  result.replaceSite(pos+1,right,false);
  result.gaugeCond(pos,'R',true);
}

void prodProj(const MPS& mps,int pos,MPS& result){
  cout<<"Projecting after "<<pos<<endl;
  mwArray right=result.getA(pos).getA();
  MPS result2(result);
  right.resize(Indices(right.getDimension(0),right.getDimension(1),1));
  result.replaceSite(pos,right,false); // do not check dimensions
  if(pos+1<=mps.getLength()-1){
    mwArray aux=result.getA(pos+1).getA();
    aux.resize(Indices(aux.getDimension(0),1,aux.getDimension(2)));
    result.replaceSite(pos+1,aux,false); // do not check dimensions
  }
  Contractor& contractor=Contractor::theContractor();
  vector<const MPS*> kets(1,&result2);
  vector<complex_t> beta(1,ONE_c);
  contractor.optimizeSum(kets,beta,result);
}

void swapProj(const MPS& mps,int pos,MPS& result){
  swap(mps,pos,result);
  MPS aux(result);
  prodProj(aux,pos+1,result);
}


void swapProj2(const MPS& mps,int pos,MPS& result,char dir){
  int len=mps.getLength();
  int d=result.getA(pos).getd();
  static bool first=true;
  static mwArray swapL;
  static mwArray swapR;
  static mwArray swapM;
  static mwArray swapMinv;
  static mwArray idM=identityMatrix(d);
  if(first){
    cout<<"Constructing matrices in swawProj2 for the first time"<<endl; 
    idM.reshape(Indices(d,1,d,1));
    mwArray swap(Indices(d,d,d,d));
    for(int i=0;i<d;i++)
      for(int j=0;j<d;j++)
	swap.setElement(ONE_c,Indices(i,j,j,i));
    swap.permute(Indices(1,3,2,4));
    swap.reshape(Indices(d*d,d*d));
    mwArray S; // the diagonal values are all one
    wrapper::svd(swap,swapL,S,swapR);
    swapL.reshape(Indices(d,d,d*d));
    swapL.permute(Indices(1,3,2));
    swapL.reshape(Indices(d*d*d,d)); // up, right, down
    swapR.reshape(Indices(d*d*d,d)); // left, up, down
    swapM=swapR*reshape(swapL,Indices(d,d*d*d));
    swapM.reshape(Indices(d*d,d,d*d,d));
    swapM.permute(Indices(2,1,4,3));
    swapMinv=reshape(swapR,Indices(d*d,d,d));
    swapMinv.permute(Indices(2,1,3));
    swapMinv.reshape(Indices(d,d*d*d));
    swapMinv.multiplyLeft(swapL);
    swapMinv.reshape(Indices(d,d*d,d*d,d));
    swapMinv.permute(Indices(1,3,4,2));
    swapR.reshape(Indices(d*d,d,d,1));
    swapR.permute(Indices(2,1,3,4));
    swapL.reshape(Indices(d,1,d*d,d));
    swapL.permute(Indices(1,2,4,3));
    cout<<"idOp:"<<idM.getDimensions()<<endl;
    cout<<"swapL:"<<swapL.getDimensions()<<endl;
    cout<<"swapR:"<<swapR.getDimensions()<<endl;
    cout<<"swapM:"<<swapM.getDimensions()<<endl;
    cout<<"swapMinv:"<<swapMinv.getDimensions()<<endl;
    first=false;
  }
  MPO swap(len);
  MPO swapB(len);
  Operator idOp(idM),sL(swapL),sR(swapR),sM(swapM),sMi(swapMinv);
  for(int k=0;k<pos;k++){
    swap.setOp(k,&idOp);
    swapB.setOp(k,&idOp);
  }
  swap.setOp(pos,&sL);
  swapB.setOp(pos,&sL);
  swap.setOp(len-1,&sR);
  swapB.setOp(len-1,&sR);
  for(int k=pos+1;k<len-1;k++){
    swap.setOp(k,&sM);
    swapB.setOp(k,&sMi);
  }

  MPS result2(mps);
  mwArray right=result.getA(pos).getA();
  right.resize(Indices(right.getDimension(0),right.getDimension(1),1));
  result.replaceSite(pos,right,false); // do not check dimensions
  if(pos+1<=mps.getLength()-1){
    mwArray aux=result.getA(pos+1).getA();
    aux.resize(Indices(aux.getDimension(0),1,aux.getDimension(2)));
    result.replaceSite(pos+1,aux,false); // do not check dimensions
  }
  //cout<<"New MPS result "<<result<<endl;
  Contractor& contractor=Contractor::theContractor();
  switch(dir){
  case 'R':
    cout<<"Optimizing swap on mps"<<endl;
    contractor.optimize(swap,mps,result);
    cout<<"Optimizing swapB on result"<<endl;
    contractor.optimize(swapB,result,result2);
    break;
  default:
    contractor.optimizeL(swapB,mps,result);
    contractor.optimizeL(swap,result,result2);
  }
  result=result2;
}


void buildInitialSymMPS(MPS& mps,mwArray& state){
  // the structure  
  mwArray Cmid(Indices(2,2,2));
  mwArray Cleft(Indices(1,2,2));
  mwArray Cright(Indices(2,1,2));
  // FOR 0001....
   // Cmid.setElement(ONE_c,Indices(0,0,0));
   // Cmid.setElement(ONE_c,Indices(1,1,1));
   // Cmid.setElement(ONE_c,Indices(0,1,1));
   // Cleft.setElement(ONE_c,Indices(0,1,1));
   // Cleft.setElement(ONE_c,Indices(0,0,0));
   // Cright.setElement(ONE_c,Indices(1,0,1));
   // Cright.setElement(ONE_c,Indices(0,0,0));
   // Cright.setElement(ONE_c,Indices(0,0,1));

  // FOR 1110000....
  Cmid.setElement(ONE_c,Indices(0,0,0));
  Cmid.setElement(ONE_c,Indices(1,1,1));
  Cmid.setElement(ONE_c,Indices(1,0,0));
  /*  Cleft.setElement(ONE_c,Indices(0,1,1));
  Cleft.setElement(ONE_c,Indices(0,0,0));
  Cright.setElement(ONE_c,Indices(0,0,0));
  Cright.setElement(ONE_c,Indices(1,0,1));*/

  int len=mps.getLength();
  for(int k=1;k<len-1;k++){
    mwArray C(Cmid);
    //if(k==0) C=Cleft;
    //else if(k==len-1) C=Cright;
    int lD=C.getDimension(0);
    int rD=C.getDimension(1);
    int d=mps.getA(k).getd();
    int Dl=mps.getA(k).getDl()/lD;
    int Dr=mps.getA(k).getDr()/rD;
    mwArray A(Indices(2,d,Dl,Dr));
    A.fillRandom();// but the 0 component is a delta
    for(int p=0;p<d;p++)
      for(int i=0;i<Dl;i++)
	for(int j=0;j<Dr;j++){
	  complex_t val0=(p==0&&(i==j))?ONE_c:ZERO_c;
	  //if((k==0||k==len-1)&&p==0)
	  //val0=(complex_t){random()*1./RAND_MAX,random()*1./RAND_MAX};
	  A.setElement(val0,Indices(0,p,i,j));
	  A.setElement(ZERO_c,Indices(1,0,i,j));
	}
    C.reshape(Indices(-1,2));
    //cout<<"C="<<C<<endl;
    A.reshape(Indices(2,-1));
    //cout<<"A="<<A<<endl;
    C.multiplyRight(A);
    C.reshape(Indices(lD,rD,d,Dl,Dr));
    C.permute(Indices(3,1,4,2,5));
    C.reshape(Indices(d,lD*Dl,rD*Dr));
    mps.setA(k,C);
  } 
  // Now the edges, which are special
  //left one
  int d=mps.getA(0).getd();
  int dp=state.getDimension(0); // physical dimension
  mwArray A(Indices(d,2,dp,dp));
  A.fillWithZero();
  // the 0 is state
  for(int l1=0;l1<dp;l1++)
    for(int l2=0;l2<dp;l2++){
      A.setElement(state.getElement(Indices(l1,l2)),Indices(0,0,l1,l2));
      // the 1 is random
      for(int k=1;k<d;k++){
	complex_t val=(complex_t){random()*1./RAND_MAX,random()*1./RAND_MAX};
	A.setElement(val,Indices(k,1,l1,l2));
      }
    }
  A.reshape(Indices(d,1,2*dp*dp));
  int realD=mps.getA(0).getDr();
  if(realD>2*dp*dp) A.resize(Indices(d,1,realD));
  mps.setA(0,A);
  // right edge
  A.reshape(Indices(d,2,dp,dp));
  A.fillWithZero();
  // the 0 is identity
  for(int l1=0;l1<dp;l1++){
    A.setElement(ONE_c,Indices(0,0,l1,l1));
    A.setElement(ONE_c,Indices(0,1,l1,l1));
    for(int l2=0;l2<dp;l2++){
      // the 1 is random
      for(int k=1;k<d;k++){
	complex_t val=(complex_t){random()*1./RAND_MAX,random()*1./RAND_MAX};
	A.setElement(val,Indices(k,1,l1,l2));
      }
    }
  }
  A.reshape(Indices(d,2*dp*dp,1));
  realD=mps.getA(len-1).getDl();
  if(realD>2*dp*dp) A.resize(Indices(d,realD,1));
  mps.setA(len-1,A);

}


void buildSymMPO(MPO& mpo,int Dop,int len){
  mpo.initLength(len);
  mwArray Zop(Indices(2,1+Dop*Dop,1+Dop*Dop)); // operators
  Zop.fillWithZero();
  Zop.setElement(ONE_c,Indices(0,0,0));
  for(int k=0;k<Dop*Dop;k++)
    Zop.setElement(ONE_c,Indices(1,k+1,k+1));
  Zop.reshape(Indices(2,-1));
  // coefficients
  if(len>2){
    mwArray Cop(Indices(2,2,2));
    Cop.setElement(ONE_c,Indices(0,0,0));
    Cop.setElement(ONE_c,Indices(1,0,0));
    Cop.setElement(ONE_c,Indices(1,1,1));
    mwArray op(Cop);
    op.reshape(Indices(2*2,2));
    op.multiplyRight(Zop);
    op.reshape(Indices(2,2,1+Dop*Dop,1+Dop*Dop));
    op.permute(Indices(3,1,4,2));
    Operator* theop=new Operator(op);
    mpo.setOp(1,theop,true);
    for(int k=2;k<len-1;k++){
      mpo.setOp(k,theop,false);    
    }
  }
  // edge operators
  mwArray CopL(Indices(1,2,2));
  CopL.setElement(ONE_c,Indices(0,0,0));
  CopL.setElement(ONE_c,Indices(0,1,1));
  mwArray op(CopL);
  op.reshape(Indices(1*2,2));
  op.multiplyRight(Zop);
  op.reshape(Indices(1,2,1+Dop*Dop,1+Dop*Dop));
  op.permute(Indices(3,1,4,2));
  mpo.setOp(0,new Operator(op),true);
  // and right
  mwArray CopR(Indices(2,1,2));
  CopR.setElement(ONE_c,Indices(0,0,0));
  CopR.setElement(ONE_c,Indices(1,0,1));
  CopR.setElement(ONE_c,Indices(1,0,0));
  mwArray opR(CopR);
  opR.reshape(Indices(2*1,2));
  opR.multiplyRight(Zop);
  opR.reshape(Indices(2,1,1+Dop*Dop,1+Dop*Dop));
  opR.permute(Indices(3,1,4,2));
  mpo.setOp(len-1,new Operator(opR),true);
}

#define MAXITER 2000 // max nr of trials for eigenvector
#define MAXPERT 1 // max nr of perturbations and retrials for eigenvector

void findMaxEigenvector(const MPO& ops,const MPO& proj,int D,
//void findMaxEigenvector(const MPO& ops,int D,
			double* lambda,char dir,
			MPS& lastMPS,double convtol){
  Contractor& contractor=Contractor::theContractor();
  if(D!=lastMPS.getBond())
    lastMPS.increaseBondDimension(D);
  lastMPS.gaugeCond('R',1);lastMPS.gaugeCond('L',1);
  // Repeat the following procedure:
  // (1) Apply Operator prj -> find new (unnormalized) MPS
  // (2) Check if normalization factor converges
  // (3) If not, normalize obtained state and repeat from (1)
  bool done=0;int count=0;
  *lambda=0.; // initialize
  int nrkicks=0,perturbcnt=0; // how many times I have retried and perturbed (todo: set a limit?)
  //cout<<"Contractor find MaxEigenvector("<<dir<<") with initial MPS "
  //<<lastMPS<<endl;
  while(!done){
    count++; // count number of attempts => set a limit
    cout<<"Contractor::findMaxEigenvector"<<dir<<"["<<count<<"] lamb="
   	<<*lambda<<"; "<<endl;
    // 6-10-2008 Try with approximation!!!
    MPS aux(lastMPS); // on which ops is acting (constant after optimize)
    if(count<2){
      switch(dir){
      case 'R':{
	lastMPS.approximate(ops,aux,D);
	break;
      }
      case 'L':{
	lastMPS.approximate(ops,aux,D,'U'); // op acting upwards
	break;
      }
      }
      cout<<"First round: approximation for initial guess"<<endl;
    }
    //    int Dop=(int)sqrt(ops.getOp(0).getd()-1);
    //projOnSym(ops.getLength(),Dop,lastMPS,dir);
    switch(dir){
    case 'R':
      contractor.optimize(proj,aux,lastMPS);
      lastMPS.gaugeCond('r',1);
      aux=lastMPS;
      contractor.optimize(ops,aux,lastMPS);
      break;
    case 'L':
      contractor.optimizeL(proj,aux,lastMPS);
      lastMPS.gaugeCond('r',1);
      aux=lastMPS;
      contractor.optimizeL(ops,aux,lastMPS);
      break;
    }
    
    double newLamb=lastMPS.getNormFact();
    // At the end, check convergence
    if(abs(newLamb-*lambda)<abs((*lambda)*convtol)){
      cout<<"Contractor::find "<<dir<<" MaxEigenvector converged after "<<count
	  <<" attempts with lambda="<<newLamb<<" and tol="<<convtol<<endl;
      done=1;
    }
    else{
      //cout<<"last var="<<abs(newLamb-*lambda)/abs(*lambda)<<"; ";
      if(count>MAXITER){
	cout<<"No convergence of findMaxEigenvector after "<<count<<" trials"<<endl;
	lastMPS.perturb(.01); // Perturbation on 1 out of 100 sites
	count=0;
	perturbcnt++;
	if(perturbcnt>=MAXPERT){
	  cout<<"Contract could not find an eigenvector after "<<perturbcnt
	      <<" trials (each of them = perturbation-"<<MAXITER<<" applications"
	      <<" of the operator=> FAILED"<<endl;
	  exit(1); // TODO: Seguir de todas formas??
	}
      }    
    }
    *lambda=newLamb;
    lastMPS.gaugeCond('R',1);lastMPS.gaugeCond('L',1); 
    // Check again?
     // cout<<"Now <vec "<<dir<<"|O|old vec>="<<contract(aux,ops,lastMPS)<<", and <vec|vec>="
//  	<<contract(lastMPS,lastMPS)<<", <oldV|oldV>="<<contract(aux,aux)<<", <vec|oldV>="
//  	<<contract(aux,lastMPS)<<endl;
  }

}

complex_t computeAllProjections(int len,int pos,int Dop,MPS& rEV,MPS& lEV,
				MPO& mpo,MPO& idMPO,ofstream* out){
  mwArray Zop(Indices(1+Dop*Dop,1,1+Dop*Dop,1)); // operators
  Zop.fillWithZero();
  Zop.setElement(ONE_c,Indices(0,0,0,0));
  static int cnt=0; // how many not present
  static Operator* op1=new Operator(Zop);
  Zop.setElement(ZERO_c,Indices(0,0,0,0));
  for(int k=0;k<Dop*Dop;k++)
    Zop.setElement(ONE_c,Indices(k+1,0,k+1,0));
  static Operator* opU=new Operator(Zop);
  static MPO proj(len);
  for(int k=0;k<len;k++)
    proj.setOp(k,op1,0);
  for(int k=0;k<pos;k++){
    proj.setOp(k,opU,0);
  }

  // Now compute projections and value
  Contractor& contractor=Contractor::theContractor();
  complex_t value=contractor.contract(rEV,proj,lEV);
  cout<<"Before contracting: <lEV|proj("<<pos<<")|rEV>="
      <<value<<endl;
  *out<<abs(value)<<"\t";
  if(abs(value)<1E-10){
    cnt++;
    cout<<"Not present: "<<cnt<<endl;
    return ZERO_c;
  }
  MPS rEVk(rEV);
  contractor.optimize(proj,rEV,rEVk);
  double overR=rEVk.getNormFact();
  MPS lEVk(lEV);
  contractor.optimizeL(proj,lEV,lEVk);
  double overL=lEVk.getNormFact();
  cout<<"Projections for pos "<<pos<<" have factors: "<<overR
      <<"(R) and "<<overL<<"(L)"<<endl;
  rEVk.gaugeCond('R',1);
  MPS::gaugeCond(rEVk,lEVk,'R');
  complex_t valD=contractor.contract(rEVk,idMPO,lEVk);
  // Now compute the expectation value
  complex_t valN=contractor.contract(rEVk,mpo,lEVk);
  cout<<"Denominator for pos "<<pos<<"="<<valD<<" and numerator="
      <<valN<<endl;
  //complex_t valN=contractor.contract(rEVk,mpo,lEV);
  // assume denominator is one (only need that lambda is one anbd they
  // really are eigenvectors)
  //return overL*overR*valN/valD;
  return valN/valD;
}

void projOnSym(int len,int Dop,MPS& mps,char dir){
  mwArray Zop(Indices(1+Dop*Dop,1,1+Dop*Dop,1)); // operators
  Zop.fillWithZero();
  Zop.setElement(ONE_c,Indices(0,0,0,0));
  static Operator* op1=new Operator(Zop);
  Zop.setElement(ZERO_c,Indices(0,0,0,0));
  for(int k=0;k<Dop*Dop;k++)
    Zop.setElement(ONE_c,Indices(k+1,0,k+1,0));
  static Operator* opU=new Operator(Zop);
  static MPO proj(len);
  for(int k=0;k<len;k++)
    proj.setOp(k,op1,0);
  std::vector<complex_t> values(len+1,ZERO_c);
  Contractor& contractor=Contractor::theContractor();
  int nUs=0; // projector with NUs
  while(nUs<=len){
    values[nUs]=contractor.contract(mps,proj,mps);
    cout<<"Size of the projection at "<<nUs<<"="<<values[nUs]
	<<endl;
    // set the site nUs to the U one except for the last
    if(nUs<len) proj.setOp(nUs,opU,0);    
    nUs++;
  }

  // now build the true MPO for "projection"
  static MPO mpo(len);
  static bool first(true);
  if(first){//initialize the MPO
    cout<<"Initializing the MPO"<<endl;
    mwArray Zop_(Indices(2,1+Dop*Dop,1+Dop*Dop)); // operators
    Zop_.fillWithZero();
    Zop_.setElement(ONE_c,Indices(0,0,0));
    for(int k=0;k<Dop*Dop;k++)
      Zop_.setElement(ONE_c,Indices(1,k+1,k+1));
    Zop_.reshape(Indices(2,-1));
    mwArray Cop(Indices(2,2,2));
    Cop.setElement(ONE_c,Indices(0,0,0));
    Cop.setElement(ONE_c,Indices(1,0,0));
    Cop.setElement(ONE_c,Indices(1,1,1));
    mwArray op(Cop);
    op.reshape(Indices(2*2,2));
    op.multiplyRight(Zop_);
    op.reshape(Indices(2,2,1+Dop*Dop,1+Dop*Dop));
    op.permute(Indices(3,1,4,2));
    for(int k=1;k<len-1;k++)
      mpo.setOp(k,new VariableOperator(op),true);
    // Set the edges too
    mwArray CopL(Indices(1,2,2));
    CopL.setElement(ONE_c,Indices(0,0,0));
    CopL.setElement(ONE_c,Indices(0,1,1));
    mwArray opL(CopL);
    opL.reshape(Indices(1*2,2));
    opL.multiplyRight(Zop_);
    opL.reshape(Indices(1,2,1+Dop*Dop,1+Dop*Dop));
    opL.permute(Indices(3,1,4,2));
    mpo.setOp(0,new VariableOperator(opL),true);
    mwArray CopR(Indices(2,1,2));
    CopR.setElement(ONE_c,Indices(0,0,0));
    CopR.setElement(ONE_c,Indices(1,0,1));
    CopR.setElement(ONE_c,Indices(1,0,0));
    mwArray opR(CopR);
    opR.reshape(Indices(2*1,2));
    opR.multiplyRight(Zop_);
    opR.reshape(Indices(2,1,1+Dop*Dop,1+Dop*Dop));
    opR.permute(Indices(3,1,4,2));
    mpo.setOp(len-1,new VariableOperator(opR),true);
    first=false;
  }
  cout<<"Initialization of MPO behind"<<endl;
  // Now reset the required components.
 // HACK MUUUUYYY FEO!!!!
 VariableOperator* ref=(VariableOperator*)dynamic_cast<const VariableOperator*>(&mpo.getOp(0));
 ref->setElement(ONE_c/values[0],Indices(0,0,0,0));
 ref=(VariableOperator*)dynamic_cast<const VariableOperator*>(&mpo.getOp(len-1));
 ref->setElement(ONE_c/values[len-1],Indices(0,1,0,0));
 for(int k=0;k<Dop*Dop;k++)
   ref->setElement(ONE_c/values[len],Indices(k,1,k,0));
 for(int k=1;k<len-1;k++){
   ref=(VariableOperator*)dynamic_cast<const VariableOperator*>(&mpo.getOp(k));
   ref->setElement(ONE_c/values[k],Indices(0,1,0,0));
 }
 // And finally, apply the operator
 MPS aux(mps);
 switch(dir){
 case 'L':
   contractor.optimizeL(mpo,aux,mps);
   break;
 default:
   contractor.optimize(mpo,aux,mps);
 }
 // And normalize again (?)
 mps.gaugeCond('R',1);
}

 
