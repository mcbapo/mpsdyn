#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

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
    Here I try to initialize a sym state and, at the end, I try to 
    project on each of the M+1) subspaces to get corresponding vectors

*/


void buildInitialSymMPS(MPS& mps);

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
   mwArray Z(Indices(2,d*d,Dop*Dop+1,d*d,Dop*Dop+1));
   for(int i=0;i<d*d;i++)
     for(int j=0;j<d*d;j++){
       Z.setElement(identOp.getElement(Indices(i,j)),
		    Indices(0,i,0,j,0));
       for(int l=0;l<Dop*Dop;l++)
	 for(int r=0;r<Dop*Dop;r++){
	   Z.setElement(doubleU.getElement(Indices(i,l,j,r)),
			Indices(1,i,l+1,j,r+1));
	 }
     }

   // Now construct the Operator itselt
   C.reshape(Indices(-1,2));
   Z.reshape(Indices(2,-1));
   mwArray data=C*Z;
   data.reshape(Indices(2,2,d*d,Dop*Dop+1,d*d,Dop*Dop+1));
   data.permute(Indices(1,3,4,2,5,6));
   data.reshape(Indices(2*d*d,Dop*Dop+1,2*d*d,Dop*Dop+1));
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
       edgeSt.setElement(dataX.getElement(Indices(i,j)),Indices(0,i,j));
       edgeSt.setElement(dataX.getElement(Indices(i,j)),Indices(1,i,j));
     }
   edgeSt.reshape(Indices(1,2*d*d));
   cout<<"Created operator edge X "<<edgeSt.getDimensions()<<endl;
   edgeSt.multiplyRight(reshape(permute(data,Indices(2,1,3,4)),
				Indices(2*d*d,-1)));
   edgeSt.reshape(Indices(1,Dop*Dop+1,Dop*Dop+1,2*d*d));
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
       //idData.setElement(lastId.getElement(Indices(i,j)),Indices(0,i,j));
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
       //opData.setElement(sigX.getElement(Indices(i,j)),Indices(0,i,j));
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
   for(int k=1;k<len-1;k++)
     trMPO.setOp(k,oper);
   trMPO.setOp(0,operX);
   trMPO.setOp(len-1,operId);

   if(M<=6){
     // to check, save the operator MPO as text
     trMPO.exportMPOtext("/home/banulsm/SVN/peps/projects/mpsdyn/src/bin/transv1.txt");
   }

   MPS anEV(len,D,Dop*Dop+1);
   // Make it "symmetric", somehow
   if(len>1)
     buildInitialSymMPS(anEV);
   MPS _rightEV(anEV);
   MPS _leftEV(anEV);
   double lambda=0;
   cout<<"Starting the optimization"<<endl;
   contractor.findRightEigenvector(trMPO,D,&lambda,_rightEV);
   contractor.findLeftEigenvector(trMPO,D,&lambda,_leftEV);
   MPS::gaugeCond(_rightEV,_leftEV,'R');

  if(M<=6){
     // to check, save the MPSs too as text
     rightEV[0].exportMPStext("/home/banulsm/SVN/peps/projects/mpsdyn/src/bin/rightEV.txt");
     leftEV[0].exportMPStext("/home/banulsm/SVN/peps/projects/mpsdyn/src/bin/leftEV.txt");
   }

   // // Now I would like a vector orthogonal to the first one.
   // MPS rightEV2(len,D,Dop*Dop+1);
   // rightEV2.setRandomState();
   // double lambda2=0;
   // complex_t overlap=contractor.contract(rightEV2,rightEV);
   // cout<<"Initialized second rightEV candidate (overlap "<<overlap
   //     <<")"<<endl;
   // contractor.findRightEigenvector(trMPO,D,&lambda2,rightEV2);
   // cout<<"After Optimizing, overlap "
   //     <<contractor.contract(rightEV2,rightEV)<<endl;
   // // for(int p=0;p<3;p++){
   // // {  
   // //   MPS aux(rightEV2);
   // //   std::vector<const MPS*> rhs(2);
   // //   rhs[0]=&aux;rhs[1]=&rightEV;
   // //   std::vector<complex_t> beta(2);
   // //   beta[0]=ONE_c;
   // //   beta[1]=-contractor.contract(rightEV2,rightEV);
   // //   contractor.optimizeSum(rhs,beta,rightEV2);
   // // }
   // // cout<<"After orthogonalizing for "<<p<<" time, overlap "
   // //     <<contractor.contract(rightEV2,rightEV)<<endl;
   // // }
   // MPS leftEV2(rightEV2);
   // contractor.findLeftEigenvector(trMPO,D,&lambda,leftEV2);
   //  // for(int p=0;p<3;p++){
   // // {  
   // //   MPS aux(leftEV2);
   // //   std::vector<const MPS*> rhs(2);
   // //   rhs[0]=&aux;rhs[1]=&leftEV;
   // //   std::vector<complex_t> beta(2);
   // //   beta[0]=ONE_c;
   // //   beta[1]=-contractor.contract(leftEV2,leftEV);
   // //   contractor.optimizeSum(rhs,beta,leftEV2);
   // // }
   // // cout<<"After orthogonalizing for "<<p<<" time, overlap "
   // //     <<contractor.contract(leftEV2,leftEV)<<endl;
   // // }
   // MPS::gaugeCond(rightEV2,leftEV2,'R');
   // cout<<"Second leftEV candidate (after gauge overlap "<<
   //   contractor.contract(leftEV2,leftEV)<<")"<<endl;
   
   std::vector<MPS> rightEV(M+1,_rightEV);
   std::vector<MPS> leftEV(M+1,_leftEV);

   std::vector<complex_t> denom(maxiter);
   std::vector<complex_t> num(maxiter);
   for(int k=0;k<maxiter;k++){
     // should be just lambda, with the gauge above
     denom[k]=contractor.contract(rightEV[k],trMPO,leftEV[k]);
   }
   //complex_t denom=lambda*ONE_c;
   
   // the operator with a sigx at the end
   trMPO.setOp(len-1,operOp);

   for(int k=0;k<maxiter;k++){
     num[k]=contractor.contract(rightEV[k],trMPO,leftEV[k]);
   }

   double time=M*delta;

   cout<<"For total time "<<time<<", values= ";
   complex_t avrg=ZERO_c;
   for(int k=0;k<maxiter;k++){
     cout<<num[k]/denom[k]<<", ";
     avrg=avrg+num[k]/denom[k];
   }
   cout<<endl;


   *out<<M<<"\t"<<delta<<"\t"<<time+delta<<"\t"<<D<<"\t"<<avrg/maxiter;
   for(int k=0;k<maxiter;k++)
     *out<<"\t"<<num[k]/denom[k];
   *out<<endl;

   delete oper;
   delete operX;
   delete operId;
   delete operOp;
   out->close();
}

void buildInitialSymMPS(MPS& mps){
  // the structure
  mwArray Cmid(Indices(2,2,2));
  Cmid.setElement(ONE_c,Indices(0,0,0));
  Cmid.setElement(ONE_c,Indices(1,1,1));
  //Cmid.setElement(ONE_c,Indices(1,0,0));
  Cmid.setElement(ONE_c,Indices(0,1,0));
  mwArray Cleft(Indices(1,2,2));
  Cleft.setElement(ONE_c,Indices(0,1,1));
  Cleft.setElement(ONE_c,Indices(0,0,0));
  Cleft.setElement(ONE_c,Indices(0,1,0));
  mwArray Cright(Indices(2,1,2));
  Cright.setElement(ONE_c,Indices(1,0,1));
  Cright.setElement(ONE_c,Indices(0,0,0));
  //Cright.setElement(ONE_c,Indices(1,0,0));
  int len=mps.getLength();
  for(int k=0;k<len;k++){
    mwArray C(Cmid);
    if(k==0) C=Cleft;
    else if(k==len-1) C=Cright;
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
	  if(k==0||k==len-1)
	    val0=(complex_t){random()*1./RAND_MAX,random()*1./RAND_MAX};
	  A.setElement(val0,Indices(0,p,i,j));
	  A.setElement(ZERO_c,Indices(1,0,i,j));
	}
    C.reshape(Indices(-1,2));
    A.reshape(Indices(2,-1));
    C.multiplyRight(A);
    C.reshape(Indices(lD,rD,d,Dl,Dr));
    C.permute(Indices(3,1,4,2,5));
    C.reshape(Indices(d,lD*Dl,rD*Dr));
    mps.setA(k,C);

  } 




}
