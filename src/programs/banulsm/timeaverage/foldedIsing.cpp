/**
   \file foldedIsing.cpp
   Program that implements a folded contraction for the Ising Hamiltonian.
   
   \author Mari-Carmen Banuls
   \date 21/10/2011
*/

#include <math.h>
#include <iomanip>
#include <typeinfo>

#include "MPS.h"
#include "Contractor.h"
#include "VariableOperator.h"
#include "DoubleOperator.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;

/** foldedIsing constructs (by hand) a transverse operator to try to 
    compare the performance of the folded (std) technique with the 
    average MPO in averageIsing (since the old implementation is 
    not cooperative)
    only for Ising Hamiltonian \ref <IsingHamiltonian>.

    Receives arguments:
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <delta> (double) width of the time step used
    \param <M> (int) number of steps
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the energy
    \param <app> (int) whether the output file is to be kept (app==1)
*/

complex_t computeIntermediateValue(int pos,int Dop,MPS& rEV,MPS& lEV,
				   MPO& oper,int d);

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
   // This is the du x Dl x dd x Dr basic piece
   // I want to combine with the adjoints for the ancillary side

   // TEST!!!!
   mwArray basicUad(basicU);
   basicUad.conjugate(); // the adjoint for the ancillas

#ifdef TESTINGDOUB
   // Now I try a true double operator
   basicU.permute(Indices(2,3,4,1)); // as transverse coming from right
   basicUad.permute(Indices(2,3,4,1)); // as transverse coming from right, after folding
   DoubleOperator* oper=new DoubleOperator(basicU,basicUad);
   //FoldedOperator* oper=new FoldedOperator(basicU,basicUad);

   // And double operators, too, for states and observables
   // 1) the original state
   int Dst=1; // bond dimension of original state
   mwArray stateX(Indices(2,Dst,Dst));
   // particular case: product |X+>
   stateX.setElement(1./sqrt(2)*ONE_c,Indices(0,0,0));
   stateX.setElement(1./sqrt(2)*ONE_c,Indices(1,0,0));
   stateX.reshape(Indices(d,Dst,Dst,1));
   stateX.permute(Indices(2,4,3,1));
   mwArray adjSt(conjugate(stateX)); // The adjoint
   DoubleOperator* operX=new DoubleOperator(stateX,adjSt);
   //FoldedOperator* operX=new FoldedOperator(permute(stateX,Indices(1,4,3,2)),adjSt);
   // This one contains no time step yet!!

   //2) identity operator on the other edge: for now, this edge is
   //just a block, but one should include the structure in the
   //contractions, to make them more efficient (another option: using
   //DoubleOperator and an extra site with dimension 1

   // This constructs a blocked double operator
   mwArray doubleU(basicU);
   doubleU.reshape(Indices(-1,1));
   basicUad.reshape(Indices(1,-1));
   doubleU.multiplyRight(basicUad);
   doubleU.reshape(Indices(Dop,d,Dop,d,Dop,d,Dop,d));
   doubleU.permute(Indices(1,5,2,6,3,7,4,8));
   doubleU.reshape(Indices(Dop*Dop,d*d,Dop*Dop,d*d));
   // identity op
   mwArray lastId=identityMatrix(d);
   mwArray dataAux(doubleU); // to multiply with std oper, unless there is only 1 step, and the state comes here

   int Dl=dataAux.getDimension(0);
   int Dr=dataAux.getDimension(2);
   int dd=dataAux.getDimension(1);
   int du=dataAux.getDimension(3);
   lastId.reshape(Indices(du,1));
   lastId.multiplyLeft(reshape(dataAux,Indices(-1,du)));
   lastId.reshape(Indices(Dl,dd,Dr,1));
   Operator* operId=new Operator(lastId);

   // at the end op sigmaX (contracted with a middle op!
   complex_t sigmaX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
   mwArray sigX(Indices(d,d),sigmaX);
   sigX.reshape(Indices(d*d,1));
   sigX.multiplyLeft(reshape(dataAux,Indices(-1,du)));
   sigX.reshape(Indices(Dl,dd,Dr,1));
   Operator* operOp=new Operator(sigX);

#else
   // This constructs a blocked double operator
   mwArray doubleU(basicU);
   doubleU.reshape(Indices(-1,1));
   basicUad.reshape(Indices(1,-1));
   doubleU.multiplyRight(basicUad);
   doubleU.reshape(Indices(d,Dop,d,Dop,d,Dop,d,Dop));
   doubleU.permute(Indices(1,5,2,6,3,7,4,8));
   doubleU.reshape(Indices(d*d,Dop*Dop,d*d,Dop*Dop));

   // In this case, this is all the MPO, besides the edges
   doubleU.permute(Indices(2,3,4,1)); // as transverse coming from right
   Operator* oper=new Operator(doubleU);

   // now the edges (contracted with a middle op!)
   mwArray stateX(Indices(2,1));
   stateX.setElement(1./sqrt(2)*ONE_c,Indices(0,0));
   stateX.setElement(1./sqrt(2)*ONE_c,Indices(1,0));
   mwArray adjSt(Hconjugate(stateX));
   mwArray dataX=stateX*adjSt;
   dataX.reshape(Indices(1,d*d));
   dataX.multiplyRight(reshape(permute(doubleU,Indices(2,1,3,4)),
				Indices(d*d,-1)));
   dataX.reshape(Indices(1,Dop*Dop,Dop*Dop,d*d));
   dataX.permute(Indices(2,1,3,4));
   Operator* operX=new Operator(dataX);

   // identity op
   mwArray lastId=identityMatrix(d);
   mwArray dataAux(doubleU); // to multiply with std oper, unless there is only 1
   if(M==1)
     dataAux=dataX;
   int Dl=dataAux.getDimension(0);
   int Dr=dataAux.getDimension(2);
   int dd=dataAux.getDimension(1);
   int du=dataAux.getDimension(3);
   lastId.reshape(Indices(du,1));
   lastId.multiplyLeft(reshape(dataAux,Indices(-1,du)));
   lastId.reshape(Indices(Dl,dd,Dr,1));
   Operator* operId=new Operator(lastId);

   // at the end op sigmaX (contracted with a middle op!
   complex_t sigmaX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
   mwArray sigX(Indices(d,d),sigmaX);
   sigX.reshape(Indices(d*d,1));
   sigX.multiplyLeft(reshape(dataAux,Indices(-1,du)));
   sigX.reshape(Indices(Dl,dd,Dr,1));
   Operator* operOp=new Operator(sigX);
#endif   


   // An attempt to find a right/left eigenvector
#ifdef TESTINGDOUB
   int len=M+1; //ops + starting state
#else
   int len=M;
#endif
   MPO trMPO(len);
   MPO trMPOx(len);
   for(int k=1;k<len-1;k++){
     trMPO.setOp(k,oper);
     trMPOx.setOp(k,oper);
   }
   //cout<<"Set operators pos 1-"<<len-2<<" in trMPO(x) to "<<*oper<<endl;
   trMPO.setOp(0,operX);
   //cout<<"Set operator pos "<<0<<" in trMPO to "<<*operX<<endl;
   trMPO.setOp(len-1,operId);
   //cout<<"Set operator pos "<<len-1<<" in trMPO to "<<*operId<<endl;
   trMPOx.setOp(0,operX);
   //cout<<"Set operator pos "<<0<<" in trMPOx to "<<*operX<<endl;
   trMPOx.setOp(len-1,operOp);
   //cout<<"Set operator pos "<<len-1<<" in trMPOx to "<<*operOp<<endl;
  // the operator with a sigx at the end

   int maxiter=1;
   MPS anEV(len,D,trMPO.getOriginalDimensions());
   // cout<<"constructed default MPS anEV "<<endl;
   anEV.setRandomState();
   //cout<<"anEV set to random state "<<anEV<<endl;

   std::vector<MPS> rightEV(maxiter,anEV);
   std::vector<MPS> leftEV(maxiter,anEV);

   for(int iter=0;iter<maxiter;iter++){
     MPS& _rightEV(rightEV[iter]);
     MPS& _leftEV(leftEV[iter]);

     double lambda=0;
     //cout<<"Starting the optimization for right EV"<<endl;
     contractor.findRightEigenvector(trMPO,D,&lambda,_rightEV);
     //cout<<"Starting the optimization for left EV"<<endl;
     contractor.findLeftEigenvector(trMPO,D,&lambda,_leftEV);
     MPS::gaugeCond(_rightEV,_leftEV,'R');
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


#ifdef TESTINGDOUB
   *out<<M<<"\t"<<delta<<"\t"<<time<<"\t"<<D<<"\t"<<avrg/maxiter;
#else
   *out<<M<<"\t"<<delta<<"\t"<<time+delta<<"\t"<<D<<"\t"<<avrg/maxiter;
#endif
   for(int k=0;k<maxiter;k++)
     *out<<"\t"<<num[k]/denom[k];
   *out<<endl;

   for(int pos=0;pos<len-1;pos++){
#ifdef TESTINGDOUB
   // Compute also all the intermediate ones (from step 0 to M-1)
     *out<<pos<<"\t"<<delta<<"\t"<<delta*(pos)<<"\t"<<D<<"\t";
#else
   // Compute also all the intermediate ones (from step 1 to M-1)
     *out<<pos+1<<"\t"<<delta<<"\t"<<delta*(pos+1)<<"\t"<<D<<"\t";
#endif
     complex_t result=computeIntermediateValue(pos,Dop,rightEV[0],
					       leftEV[0],
					       trMPO,d);
     *out<<result/denom[0]<<"\t"<<0.<<endl;
   }

   delete oper;
   delete operX;
   delete operId;
   delete operOp;
   out->close();
}

complex_t computeIntermediateValue(int pos,int Dop,MPS& rEV,MPS& lEV,
				   MPO& oper,int d){
  if(pos<0||pos>=oper.getLength()-1){
    cout<<"Invalid position "<<pos<<", only from 0 to len-2 supported"<<endl;
    exit(2);
  }
  //cout<<"Function computeIntermediateValue("<<pos<<")"<<endl;
  const Operator& op=oper.getOp(pos);
  //  cout<<"Recovered old Operator, "<<op<<endl;
  mwArray oldData=op.getFullData();
  Indices dims=oldData.getDimensions();
  // Now multiply the operator
  complex_t sigmaX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),sigmaX);
  oldData.reshape(Indices(-1,d,d));
  oldData.permute(Indices(1,3,2));
  oldData.reshape(Indices(-1,d));
  oldData.multiplyRight(sigX);
  oldData.reshape(Indices(-1,d,d));
  oldData.permute(Indices(1,3,2));
  oldData.reshape(dims); // restore the original dimensions
  // set the modified operator in place
  oper.setOp(pos,new Operator(oldData),true);
  // compute the value
  Contractor& contractor=Contractor::theContractor();
  complex_t result=contractor.contract(rEV,oper,lEV);
  // restore the old one
  oper.setOp(pos,&op);
  return result;
}

// complex_t computeIntermediateValue(int pos,int Dop,MPS& rEV,MPS& lEV,
// 				   MPO& oper,int d){
//   if(pos<0||pos>=oper.getLength()-1){
//     cout<<"Invalid position "<<pos<<", only from 0 to len-2 supported"<<endl;
//     exit(2);
//   }
//   mwArray oldOp=oper.getOpData(pos);
//   mwArray restore(oldOp); // to restore it at the end!
//   Indices dims=oldOp.getDimensions();
//   // Now multiply the operator
//   complex_t sigmaX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
//   mwArray sigX(Indices(d,d),sigmaX);
//   oldOp.reshape(Indices(-1,d,d));
//   oldOp.permute(Indices(1,3,2));
//   oldOp.reshape(Indices(-1,d));
//   oldOp.multiplyRight(sigX);
//   oldOp.reshape(Indices(-1,d,d));
//   oldOp.permute(Indices(1,3,2));
//   oldOp.reshape(dims); // restore the original dimensions
//   // set the modified operator in place
//   oper.setOp(pos,new Operator(oldOp),true);
//   // compute the value
//   Contractor& contractor=Contractor::theContractor();
//   complex_t result=contractor.contract(rEV,oper,lEV);
//   // restore the old one
//   oper.setOp(pos,new Operator(restore),true);
//   return result;
// }
