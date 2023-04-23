/**
   \file nnpert.cpp
   Read a nearest-neighbour Hamiltonian, find its ground state, apply
   a local perturbation and evolve in time.
   
   \author Mari-Carmen Banuls
   \date 10/04/2012

*/

#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "MPO.h"
#include "Contractor.h"

#include "NNHamiltonian.h"

using namespace shrt;
using namespace std;

/** nnpert reads a Hamiltonian \ref <NNHamiltonian>
    from a text file, finds its ground state and uses as initial state
    for evolution after applying a local perturbation on the central
    site. If the total length is even, the perturbation is applied on
    the site to the right of the center of the chain.

    Receives arguments:
    \param <infname> (char*) name of the input file with NNHamiltonian
    \param <infnameP> (char*) name of the input file with perturbation
    \param <D> (int) maximum bond dimension to be used
    \param <stepBl> (int) periodicity for saving results of time evolution
    \param <timefname> (char*) name of the input file containing the
    time configuration (tot nr of steps, list of delta and nrs)
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)
*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infname=argv[++cntr];
  const char* infnameP=argv[++cntr];
  int D=atoi(argv[++cntr]);
  int stepBlk=atoi(argv[++cntr]);
  const char* timefname=argv[++cntr];
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<"% NNHamiltonian from "<<infname<<endl;
  *out<<"% D="<<D<<", stepBlk="<<stepBlk<<endl;
  *out<<"% Perturbation from "<<infnameP<<endl;
  *out<<"% Evolution configuration from "<<timefname<<endl;
  *out<<setprecision(10);

  *out<<"% D\t E(GS) \t "<<endl;

  int M;

  ifstream inlog(timefname);
  double t1=0;
  inlog>>M; // true final number of steps
  double* deltas=new double[M];
  // Initialize the step widths
  double tmp;int cntD;double normFact;
  int k=0; // actual step
  while(k<M-1){ // read all the file
    inlog>>cntD;
    inlog>>tmp;
    for(int l=0;l<cntD;l++){
      deltas[k++]=tmp;
      //cout<<"Step "<<k<<" is "<<tmp<<endl;
    }
  }
  inlog.close();
  cout<<"Read time configuration: M="<<M<<", delta[0]="<<deltas[0]<<endl;

  NNHamiltonian myH(infname);
  const MPO& hamil=myH.getHMPO();
  int L=myH.getLength();
  vector<int> dims(L);
  myH.getDimensions(dims);
  MPS gsP(L,D,dims);
  gsP.setRandomState(); // the initial state, random -> not to the GS!!
  gsP.gaugeCond('R',1);
  gsP.gaugeCond('L',1);
  double lambda=0.;
  Contractor& contractor=Contractor::theContractor();
  cout<<"Before starting, expectation value is="
      <<contractor.contract(gsP,hamil,gsP)<<", and norm="
      <<contractor.contract(gsP,gsP)<<endl;

  contractor.findGroundState(hamil,D,&lambda,gsP);

  cout<<"Ground state found with eigenvalue "<<lambda<<endl;
  complex_t normGS=contractor.contract(gsP,gsP);
  complex_t energyGS=contractor.contract(gsP,hamil,gsP);
  cout<<"After optimizing, expectation value is="
      <<energyGS<<", and norm="
      <<normGS<<endl;
  *out<<D<<"\t"<<energyGS/normGS<<endl;

  // Now read and apply the perturbation
  // Site on which perturbation acts: it length is odd, in the center
  int center=L%2==0?L/2:(L-1)/2;
  ifstream inputP(infnameP);
  if(!inputP.is_open()){
    cout<<"Error: cannot open file "<<infnameP<<" to read NNHamiltonian"<<endl;
    exit(1);
  }
  mwArray pert;
  pert.loadtext(inputP);
  inputP.close();
  // Apply the perturbation now and normalize the state
  cout<<"Applying perturbation "<<pert.getDimensions()<<" to site "
      <<center<<endl;
  MPS evolS(gsP);
  evolS.applyLocalOperator(center,pert,true);
  // Normalize
  //evolS.gaugeCond('R',true);
  cout<<"After applying perturbation, energy changed from "
      <<energyGS/normGS<<" to "<<contractor.contract(evolS,hamil,evolS)<<endl;

  *out<<"% After perturbation, energy "
      <<contractor.contract(evolS,hamil,evolS)<<endl;
  *out<<"% Overlap with GS "<<contractor.contract(evolS,gsP)<<endl;


  // Now run steps of time evolution

  // Get the exponential MPOs
  MPO expHe(L),expHo(L),expHe_2(L);
  double deltaR=deltas[0];
  myH.getExponentialMPOeven(expHe,(complex_t){0,-deltaR});
  myH.getExponentialMPOeven(expHe_2,(complex_t){0,-deltaR/2});
  myH.getExponentialMPOodd(expHo,(complex_t){0,-deltaR});

  int cnt=0;  double time=0.;
  *out<<"% %%%%%% Result of time evolution %%%% "<<endl;
  *out<<"% Nr step \t delta \t time \t Energy \t overlap GS \t Entropy (cut by cut)"
      <<endl;
  out->close();

  while(cnt<M){
    cout<<"Time evolution step nr "<<cnt<<", delta="<<deltaR<<endl;
    MPS aux(evolS);
    //double errE,errO;
    contractor.optimize(expHe_2,aux,evolS,D);
    //cout<<"After step "<<cnt<<" rel errE="<<errE<<endl;
    int ks=0;
    while((ks<stepBlk-1)&&(cnt<M-1)&&(deltas[cnt+1]==deltaR)){
      aux=evolS;
      //contractor.optimize(expHo,aux,evolS,D,&errO);
      contractor.optimize(expHo,aux,evolS,D);
      aux=evolS;
      //      contractor.optimize(expHe,aux,evolS,D,&errE);
      contractor.optimize(expHe,aux,evolS,D);
      ks++;cnt++;time+=deltaR;
      //cout<<"After step "<<cnt-1<<" rel errE="<<errE
      //  <<", rel errO="<<errO<<endl;
    }
    aux=evolS;
    //contractor.optimize(expHo,aux,evolS,D,&errO);
    contractor.optimize(expHo,aux,evolS,D);
    aux=evolS;
    //contractor.optimize(expHo,aux,evolS,D,&errO);
    contractor.optimize(expHe_2,aux,evolS,D);
    //cout<<"After step "<<cnt<<" rel errO="<<errO
    //	<<" rel errE="<<errE<<endl;
    cnt++;
    time+=deltaR;
    // Normalize, just in case
    evolS.gaugeCond('R',1);
    complex_t energyM=contractor.contract(evolS,hamil,evolS);
    complex_t normM=contractor.contract(evolS,evolS);
    complex_t overlM=contractor.contract(evolS,gsP);
    out->open(outfname,ios::app);
    *out<<cnt<<"\t"<<deltaR<<"\t"<<time
	<<"\t"<<real(energyM/normM)<<"\t"<<abs(overlM)<<"\t";
    // Now I would also like to compute the entropy between the right 
    // part and the left one across every cut in the chain!
    for(int pl=1;pl<=L-1;pl++){
      double entr=contractor.getEntropy(evolS,pl);
      *out<<entr<<"\t";
    }
    *out<<endl;
    out->close();
    if(deltaR!=deltas[cnt]){
      cout<<"Recomputing exponentials, due to changing step"<<endl;
      deltaR=deltas[cnt];
      myH.getExponentialMPOeven(expHe,(complex_t){0,-deltaR});
      myH.getExponentialMPOeven(expHe_2,(complex_t){0,-deltaR/2});
      myH.getExponentialMPOodd(expHo,(complex_t){0,-deltaR});
    }
  }
  if(out->is_open())
    out->close();
  delete out;
}

/**
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  int d=N+1;
  MPS gsP(L,D,d);

  BoseHubbardHamiltonian hBH0(L,N,t,U,mu,V,Vint);
  const MPO& hamil=hBH0.getHMPO();


  // Construct a MPO for the particle number on a given site
  mwArray identM=identityMatrix(d);
  identM.reshape(shrt::Indices(d,1,d,1));
  Operator ident(identM);
  double* values=new double[d];
  for(int k=0;k<d;k++){
    values[k]=k;
  }
  mwArray nrop=realDiag(d,values); // individual operator number
  nrop.reshape(shrt::Indices(d,1,d,1));
  Operator nrk(nrop);
  // Basic MPO: all Identity
  MPO mpoN(L);
  for(int k=0;k<L;k++)
    mpoN.setOp(k,&ident,0);
  // Now iterate for occupation number
  *out<<"% t="<<t<<endl;
  *out<<"% U="<<U<<endl;
  *out<<"% mu="<<mu<<endl;
  *out<<"% V="<<V<<endl;
  *out<<"% Vint="<<Vint<<endl;
  *out<<"%k\t <n(k)> \t (Delta nk)^2"<<endl;
  // and at the same time compute total number of particles
  double Ntot=0.;
  for(int k=0;k<L;k++){
    mpoN.setOp(k,&nrk,0);
    complex_t nrk_=contractor.contract(gsP,mpoN,gsP);
    complex_t nrk2=contractor.contract2(gsP,mpoN,gsP);
    nrk_=nrk_/normGS;
    nrk2=nrk2/normGS;
    *out<<k  //<<"\t"<<t<<"\t"<<U<<"\t"<<mu<<"\t"<<V<<"\t"<<Vint
	<<"\t"<<nrk_<<"\t"<<nrk2-nrk_*nrk_<<endl;
    // restore the identity
    mpoN.setOp(k,&ident,0);
    Ntot=Ntot+real(nrk_);
  }

  cout<<"Total number of particles N="<<Ntot<<endl;
  *out<<"%\n Ntot="<<Ntot<<endl;

  // Now apply creation operator on the central site
  mwArray aDagger(Indices(d,d));
  for(int n=0;n<N;n++){
    aDagger.setElement(ONE_c*sqrt(n+1),Indices(n+1,n));
  }
  aDagger.reshape(Indices(d,1,d,1));
  int center=L%2==0?L/2:(L-1)/2;
  cout<<"Cut located in "<<center<<endl;
  *out<<"%Before perturbation, S("<<center<<")="<<contractor.getEntropy(gsP,center)<<endl;
  Operator* aDagop=new Operator(aDagger);
  mpoN.setOp(center,aDagop,false);
  MPS init(gsP);
  {
    MPS tmp(init);
    contractor.optimize(mpoN,tmp,init,D);
    mpoN.setOp(center,&ident,false);
  }
  delete aDagop;

  complex_t normInit=contractor.contract(init,init);
  complex_t energyInit=contractor.contract(init,hamil,init);

  *out<<"% Perturbation applied, now S("<<center<<")="<<contractor.getEntropy(init,center)<<endl;
  *out<<"%k\t <n(k)> \t (Delta nk)^2"<<endl;
  Ntot=0.;
  for(int k=0;k<L;k++){
    mpoN.setOp(k,&nrk,0);
    complex_t nrk_=contractor.contract(init,mpoN,init);
    complex_t nrk2=contractor.contract2(init,mpoN,init);
    nrk_=nrk_/normInit;
    nrk2=nrk2/normInit;
    *out<<k  //<<"\t"<<t<<"\t"<<U<<"\t"<<mu<<"\t"<<V<<"\t"<<Vint
	<<"\t"<<nrk_<<"\t"<<nrk2-nrk_*nrk_<<endl;
    // restore the identity
    mpoN.setOp(k,&ident,0);
    Ntot=Ntot+real(nrk_);
  }

  cout<<"Total number of particles after creating one, N="<<Ntot<<endl;
  cout<<"New energy "<<energyInit/normInit<<endl;


  // Now start the evolution
  // Get the exponential MPOs
  MPO expHe(L),expHo(L),expHe_2(L);
  hBH0.getExponentialMPOeven(expHe,(complex_t){0,-delta});
  hBH0.getExponentialMPOeven(expHe_2,(complex_t){0,-delta*.5});
  hBH0.getExponentialMPOodd(expHo,(complex_t){0,-delta});

  //  expHe.exportMPOtext("checkExpHe.txt");
  //expHe_2.exportMPOtext("checkExpHe_2.txt");
  //expHo.exportMPOtext("checkExpHo.txt");
  //  exit(1);

  int cnt=0;
  *out<<"% Results of the time evolution "<<endl;
  *out<<"% cnt(delta="<<delta<<")\t ";
  for(int k=0;k<L;k++) *out<<"n("<<k<<")\t";
  *out<<"Ntot \t E \t S(L/2)"<<endl;

  MPS evolS(init);
  while(cnt<M){
    cout<<"Time evolution step "<<cnt<<endl;
    MPS aux(evolS);
    double errE,errO;
    contractor.optimize(expHe_2,aux,evolS,D,&errE);
    cout<<"After step "<<cnt<<" rel errE="<<errE<<endl;
    int ks=0;
    while((ks<stepBlk-1)&&(cnt<M-1)){
      aux=evolS;
      contractor.optimize(expHo,aux,evolS,D,&errO);
      //contractor.optimize(expHo,aux,evolS,D);
      aux=evolS;
      contractor.optimize(expHe,aux,evolS,D,&errE);
      //contractor.optimize(expHe,aux,evolS,D);
      ks++;cnt++;
      cout<<"After step "<<cnt-1<<" rel errE="<<errE
        <<", rel errO="<<errO<<endl;
    }
    aux=evolS;
    contractor.optimize(expHo,aux,evolS,D,&errO);
    //contractor.optimize(expHo,aux,evolS,D);
    aux=evolS;
    contractor.optimize(expHe_2,aux,evolS,D,&errE);
    //contractor.optimize(expHe_2,aux,evolS,D);
    cout<<"After step "<<cnt<<" rel errO="<<errO
    	<<" rel errE="<<errE<<endl;
    cnt++;
    double time=delta*cnt;
    // // MPS aux(evolS);
    // // contractor.optimize(expHe_2,aux,evolS,D);
    // // aux=evolS;
    // // contractor.optimize(expHo,aux,evolS,D);
    // // aux=evolS;
    // // contractor.optimize(expHe_2,aux,evolS,D);
    // After step of evolution, compute distrib of part.
    // Normalize, just in case
    evolS.gaugeCond('R',1);
    *out<<cnt<<"\t";
    double Ntot=0.;
    for(int k=0;k<L;k++){
      mpoN.setOp(k,&nrk,0);
      complex_t nrk=contractor.contract(evolS,mpoN,evolS);
      *out<<real(nrk)<<"\t";
      // restore the identity
      mpoN.setOp(k,&ident,0);
      Ntot=Ntot+real(nrk);
    }
    *out<<Ntot<<"\t"
	<<real(contractor.contract(evolS,hamil,evolS))<<"\t";
    // Now I would also like to compute the entropy between the right 
    // part and the left one
    double entr=contractor.getEntropy(evolS,center);
    *out<<entr<<endl;
  }


  out->close();


}

*/
