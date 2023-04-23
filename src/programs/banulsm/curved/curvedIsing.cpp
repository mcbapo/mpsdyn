/**
   \file curvedIsing.cpp
   Evolve ground state of a deformed Ising Hamiltonian 
   after a local perturbation, and explore causal cone.
   
   \author Mari-Carmen Banuls
   \date 18/04/2012

*/

#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "HyperbolicIsingHamiltonian.h"

using namespace shrt;
using namespace std;

void computeCoefficients(int L,int L0,double lambda,vector<double>& ZZ,
			 vector<double>& X,vector<double>& Z);
void applyPerturbation(const MPS& current,MPS& result,int center);
void computeExpectation(MPS& state,const Operator& oper,ofstream& out);
void writeLocalEnergies(ostream *out,ofstream* out3,MPS& state,
			  vector<double>& coeffZZ,
			  vector<double>& coeffX,vector<double>& coeffZ);


/** curvedIsing runs the findGroundState routine with the MPO for the
    \ref <HyperbolicIsingHamiltonian>, then applies a local
    perturbation and evolves in real time for M steps using width delta.
    In the output file, we write the magnetization site by site 
    as a function of time.
    In the second file, entropies are written across every cut. 
    In the third one, we write two-body energies.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <L0> (int) length of the "flat" part in the center of the
                      chain (the deformation will act on the \f$L-L0\f$ 
                      external sites)
    \param <lambda> (double) parameter \f$\lambda\f$ of the deformation
    \param <D> (int) maximum bond dimension
    \param <M> (int) number of time steps
    \param <stepBlock> (int) frequency (number of time steps) of writing
    \param <delta> (double) step width
    \param <outfname> (char*) name of the output file for the magnetization
    \param <app> (int) whether the output file is to be kept (app==1)
    \param <outfname2> (char*) name of the output file for the entropy
    \param <app2> (int) whether the output file is to be kept (app2==1)
    \param <outfname3> (char*) name of the output file for the energy (2body)
    \param <app3> (int) whether the output file is to be kept (app3==1)
*/


int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  int L0=atoi(argv[++cntr]);
  double lambda=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  int stpBlk=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);
  const char* outfname2=argv[++cntr];
  int app2=atoi(argv[++cntr]);
  const char* outfname3=argv[++cntr];
  int app3=atoi(argv[++cntr]);

  if(L%2!=0){
    cout<<"WARNING: In curvedIsing, the length of the chain is not even"
	<<endl;
  }
  if(L%2!=L0%2){
    cout<<"Warning: the flat section in the center does not have the same "
	<<"parity as the total chain: the result will look asymmetric "<<endl;
  }
  // in the even case, the center is double, so this is the first one
  int center=L%2==0?L/2-1:(L-1)/2;

  cout<<"Initialized arguments: L="<<L
      <<", L0="<<L0
      <<", lambda="<<lambda
      <<", outfile="<<outfname
      <<", app="<<app
      <<", M="<<M
      <<", delta="<<delta
      <<", D="<<D<<endl;

  ofstream *out,*out2,*out3;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% L\t L0\t lambda\t D\t Energy"<<endl;
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }

  // Compute the coefficients for the curved problem
  vector<double> coeffZZ,coeffX,coeffZ;
  computeCoefficients(L,L0,lambda,coeffZZ,coeffX,coeffZ);
  // Construct the corresponding Hamiltonian
  HyperbolicIsingHamiltonian hamI(L,coeffZZ,coeffX,coeffZ);
  // HyperbolicIsingHamiltonian hamI(L,-1.,0.,0.,lambda);
  cout<<"Created the Hamiltonian"<<endl;

  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d=2;
  MPS gs(L,D,d);
  gs.setRandomState(); // the intial state, random
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  cout<<"Initialized random state, norm "<<contractor.contract(gs,gs)<<endl;
   
  // First: find the GS
  double E0=0.;
  const MPO& hamil=hamI.getHMPO();
  cout<<"Constructed the hamil MPO"<<endl;
  //  hamil.exportMPOtext("curvedHmpo.txt");
  cout<<"Initial value, with initial state"<<endl;
  cout<<contractor.contract(gs,hamil,gs)<<endl;

  contractor.findGroundState(hamil,D,&E0,gs);
  cout<<"Ground state found with eigenvalue "<<E0<<endl;
  *out<<"% "<<L<<"\t"<<L0<<"\t"<<lambda<<"\t"
      <<D<<"\t"<<E0<<endl;

  // ////// DEBUG
  // MPS gsFlat;
  // gsFlat.importMPS("gsflat.dat");    
  // cout<<"Overlap of GS: "<<contractor.contract(gs,gsFlat)<<endl;
  // cout<<"Norms: gs "<<contractor.contract(gs,gs)<<endl;
  // cout<<"Norms: gsF "<<contractor.contract(gsFlat,gsFlat)<<endl;
  // //// DEBUG

  // Now apply perturbation (and normalize)
  // Create an MPO for it
  MPS state(gs);
  complex_t sigmaX_[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,1,d,1),sigmaX_);
  Operator sigmaX(sigX);
  applyPerturbation(gs,state,center+1);

  // ////// DEBUG
  // MPS perFlat;
  // perFlat.importMPS("perflat.dat");    
  // cout<<"Overlap of pert GS: "<<contractor.contract(state,perFlat)<<endl;
  // cout<<"Norms: pert gs "<<contractor.contract(state,state)<<endl;
  // cout<<"Norms: pert gsF "<<contractor.contract(perFlat,perFlat)<<endl;
  // //// DEBUG


  *out<<"% Evolution in time, with "<<M<<" steps of width "<<delta<<endl;
  *out<<"% M\t t\t <sigx(pos)>"<<endl;
  out->close();
  
  if(!app2){
    out2=new ofstream(outfname2);
    *out2<<"% L\t L0\t lambda\t D\t Energy"<<endl;
    *out2<<"% "<<L<<"\t"<<L0<<"\t"<<lambda<<"\t"
	 <<D<<"\t"<<E0<<endl;
    *out2<<"% Evolution in time, with "<<M<<" steps of width "<<delta<<endl;
    *out2<<"% M\t t\t E \t <S(pos)>"<<endl;
  }
  else
    out2=new ofstream(outfname2,ios::app);
  if(!out2->is_open()){
    cout<<"Error: impossible to open file "<<outfname2<<
      " for output"<<endl;
    exit(1);
  }

  if(!app3){
    out3=new ofstream(outfname3);
    *out3<<"% L\t L0\t lambda\t D\t Energy"<<endl;
    *out3<<"% "<<L<<"\t"<<L0<<"\t"<<lambda<<"\t"
	 <<D<<"\t"<<E0<<endl;
    *out3<<"% Evolution in time, with "<<M<<" steps of width "<<delta<<endl;
    *out3<<"% M\t t\t E \t <E(pos,pos+1)>"<<endl;
  }
  else
    out3=new ofstream(outfname3,ios::app);
  if(!out3->is_open()){
    cout<<"Error: impossible to open file "<<outfname2<<
      " for output"<<endl;
    exit(1);
  }

  // Now run evolution with the exponential
  MPO evol(L);
  double t=0.;// H is constant
  int oT=2; // Trotter order
  bool imag=false; // no imaginary time evolution
  hamI.getUMPO(evol,delta,oT,imag);
  //  evol.exportMPOtext("operCurved.mat");
  int cnt=0;
  while(cnt<M){
    cout<<"Time evolution step nr "<<cnt+1<<", time="<<delta*(cnt+1)<<endl;
    // if(cnt==2){
    //   MPS evoFlat;
    //   evoFlat.importMPS("evolflat.dat");      
    //   cout<<"Overlap (modulus) of evol GS("<<cnt<<"): "
    // 	  <<abs(contractor.contract(state,evoFlat))<<endl;
    //   cout<<"Norms: evol gs "<<contractor.contract(state,state)<<endl;
    //   cout<<"Norms: evol gsF "<<contractor.contract(evoFlat,evoFlat)<<endl;
    // }
    MPS aux(state); // temporary copy
    contractor.optimize(evol,aux,state,D);
    if(cnt%stpBlk==0||cnt==M||cnt==0){
      double time=delta*(cnt+1); // after optimizing=> next step
      state.gaugeCond('R',1); // normalize
      double energy=real(contractor.contract(state,hamil,state));
      out2->open(outfname2,ios::app);
      *out2<<cnt<<"\t"<<time<<"\t"<<energy<<"\t";
      // Now I would also like to compute the entropy between the right 
      // part and the left one across every cut in the chain!
      for(int pl=1;pl<=L-1;pl++){
	double entr=contractor.getEntropy(state,pl);
	*out2<<entr<<"\t";
      }
      *out2<<endl;
      out2->close();
      out->open(outfname,ios::app);
      *out<<cnt<<"\t"<<time<<"\t";
      // done in writeLocalEnergies
      //computeExpectation(state,sigmaX,*out); *out<<endl;
      out3->open(outfname3,ios::app);
      *out3<<cnt<<"\t"<<time<<"\t"<<energy<<"\t";
      writeLocalEnergies(out,out3,state,coeffZZ,coeffX,coeffZ);
      *out<<endl;
      out->close();
      *out3<<endl;
      out3->close();
    }
    cnt++;
  }
  if(out->is_open())
    out->close();
  delete out;
}

void computeCoefficients(int L,int L0,double lambda,vector<double>& ZZ,
			 vector<double>& X,vector<double>& Z){
  ZZ.clear();X.clear();Z.clear();
  cout<<"Computing coefficients"<<endl;
  // limits of the flat region
  // first site of the flat part
  int flatXo=(L-L0)%2==0?(L-L0)/2:(L-L0-1)/2;
  // last site of the flat part
  int flatXl=flatXo+L0-1;

  for(int k=0;k<L;k++){
    double distk=0.;
    if(k<flatXo)
      distk=flatXo-k;
    else if(k>flatXl)
      distk=k-flatXl;
    else
      distk=0.;
    //    double value=exp(-lambda*distk);
    double value=lambda==0?1:(distk==0?1:1-(1-exp(lambda*distk))/(1-exp(lambda*(L-flatXl))));
    if(k<L-1)
      ZZ.push_back(value);
    X.push_back(value);
  }
}

#include "IsingHamiltonian.h";

void applyPerturbation(const MPS& current,MPS& result,int center){
  cout<<"Inside applyPerturbation"<<endl;
  int L=current.getLength();
  int d=current.getA(0).getd();
  // the MPO
  MPO perturbation(L);
  mwArray ident=identityMatrix(d);ident.reshape(Indices(d,1,d,1));
  Operator id(ident);
  for(int k=0;k<L;k++)
    if(k!=center&&k!=center+1)
      perturbation.setOp(k,&id,false);
  // The operator in the middle I do by hand, using a trick
  IsingHamiltonian auxH(2,2,0.,.5,0.);
  const MPO& hamil=auxH.getHMPO(0.);
  perturbation.setOp(center,&(hamil.getOp(0)),false);
  perturbation.setOp(center+1,&(hamil.getOp(1)),false);
  
  Contractor& contractor=Contractor::theContractor();
  // // The operators on sites center and center+1 by hand
  //  mwArray delta0(Indices(2,1),sigmaX_);  
  //  mwArray delta1(Indices(2,1),&sigmaX_[2]);
  //  sigX.reshape(Indices(1,d*d));ident.reshape(Indices(1,d*d));
  //  mwArray pertL=.5*delta0*sigX+delta1*ident;
  //  mwArray pertR=.5*delta1*sigX+delta0*ident;
  //  pertL.reshape(Indices(2,d,1,d));pertL.permute(Indices(2,3,4,1));
  //  pertR.reshape(Indices(2,d,d,1));pertR.permute(Indices(2,1,3,4));
  //  perturbation.setOp(center,new Operator(pertL),true);
  //  perturbation.setOp(center+1,new Operator(pertL),true);
  contractor.optimize(perturbation,current,result);
  // Normalize
  result.gaugeCond('R',true);   
  cout<<"Leaving applyPerturbation"<<endl;

}

void computeExpectation(MPS& state,const Operator& oper,ofstream& out){
  int len=state.getLength();
  int d=2;
  static MPO ops(len);
  static bool init(false);
  static Operator ident(Indices(d,1,d,1));
  if(!init){
    mwArray id_=identityMatrix(d);
    id_.reshape(Indices(d,1,d,1));
    ident.setData(&id_);
    // A trick to store copies of the mwArrays in the opers
    ident.permuteOp(Indices(1,2,3,4));
    init=true;
    for(int k=0;k<len;k++) ops.setOp(k,&ident,false);
  }
  Contractor& contractor=Contractor::theContractor();
  double norm=real(contractor.contract(state,state));
  for(int pos=0;pos<len;pos++){
    ops.setOp(pos,&oper,false);
    // Compute and write result
    complex_t valk=contractor.contract(state,ops,state);
    out<<real(valk)/norm<<"\t";
    ops.setOp(pos,&ident,false);
  }
}

void writeLocalEnergies(ostream *out,ofstream* out3,MPS& state,
			  vector<double>& coeffZZ,
			  vector<double>& coeffX,vector<double>& coeffZ){
  cout<<"Computing and writing local observables"<<endl;
  int len=state.getLength();
  int d=2;
  mwArray sigX,sigZ;
  static bool init(false);
  static MPO ops(len);
  if(!init){
    complex_t sigmaX_[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    sigX=mwArray(Indices(d,1,d,1),sigmaX_);
    complex_t sigmaZ_[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
    sigZ=mwArray(Indices(d,1,d,1),sigmaZ_);
  }
  static Operator ident(Indices(d,1,d,1));
  static Operator sigmaX(sigX), sigmaZ(sigZ); 
  if(!init){
    mwArray id_=identityMatrix(d);
    id_.reshape(Indices(d,1,d,1));
    ident.setData(&id_);
    // A trick to store copies of the mwArrays in the opers
    ident.permuteOp(Indices(1,2,3,4));
    init=true;
    for(int k=0;k<len;k++) ops.setOp(k,&ident,false);
  }
  Contractor& contractor=Contractor::theContractor();
  double Etot=0.;
  complex_t leftX=ZERO_c;
  double norm=real(contractor.contract(state,state));
  for(int k=0;k<len;k++){
    ops.setOp(k,&sigmaX,false);
    // Compute <sigx(pos)> result
    complex_t valXk=contractor.contract(state,ops,state);
    // Compute <sigz(pos)>
    ops.setOp(k,&sigmaZ,false);
    complex_t valZk=contractor.contract(state,ops,state);
    *out<<real(valXk)/norm<<"\t"<<real(valZk)/norm<<"\t";
    // Compute <sigz(pos) sigz(pos+1)>
    if(k<len-1){
      ops.setOp(k+1,&sigmaZ,false);
      complex_t valZZk=contractor.contract(state,ops,state);
      ops.setOp(k+1,&ident,false);
      // Write to file out3
      double energyK=real(coeffZZ[k]*valZZk+coeffX[k]*valXk)/norm;
      *out3<<energyK<<"\t";
      Etot+=energyK;
    }
    ops.setOp(k,&ident,false);
    leftX=valXk;
  }
  *out3<<Etot;
  cout<<"Total energy "<<Etot<<endl;
}
