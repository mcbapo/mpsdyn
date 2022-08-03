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

    This has a problem. Although the MPO is right, it does not have a
    unique max. eigenvector. On the contrary, it has M+1 degenerate
    ones, all with lambda=1, corresponding to the transfer matrices of
    all the possible evolutions involved, and therefore the power
    method will not work.

*/



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

   MPS rightEV(len,D,Dop*Dop+1);
   rightEV.setRandomState();
   cout<<"Starting rightEV to random state"<<endl;
   rightEV.gaugeCond('R',1);
   MPS leftEV(rightEV);
   double lambda=0;
   cout<<"Starting the optimization"<<endl;
   contractor.findRightEigenvector(trMPO,D,&lambda,rightEV);
   contractor.findLeftEigenvector(trMPO,D,&lambda,leftEV);

   MPS::gaugeCond(rightEV,leftEV,'R');

   // should be just lambda, with the gauge above
   complex_t denom=contractor.contract(rightEV,trMPO,leftEV);
   //complex_t denom=lambda*ONE_c;

   // the operator with a sigx at the end
   trMPO.setOp(len-1,operOp);

   
   if(M<=6){
     // to check, save the operator MPO as text
     trMPO.exportMPOtext("/home/banulsm/SVN/peps/projects/mpsdyn/src/bin/transvSig.txt");
   }


   complex_t num=contractor.contract(rightEV,trMPO,leftEV);
   double time=M*delta;

   cout<<"For total time "<<time<<", value "<<(num/denom)<<endl;
   cout<<"num="<<num<<", denom="<<denom<<endl;
   *out<<M<<"\t"<<delta<<"\t"<<time+delta<<"\t"<<D<<"\t"<<num/denom<<endl;

   delete oper;
   delete operX;
   delete operId;
   delete operOp;
   out->close();
}
