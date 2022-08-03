
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "FoldedOperator.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "SchwingerHamiltonian.h"

using namespace shrt;

/** schwingerThermUnfold tries to find the thermal state of the
    Schwinger Hamiltonian \ref <SchwingerHamiltonian> using
    purification and imaginary time evolution. In this case, for
    efficiency, the thermal state is written as an MPS with
    real and ancillary sites splitted, so that 0,2,4... are the
    real ones and odd numbered sites are for the conjugate dof.

    Receives arguments:
    \param <N> (int) length of the chain
    \param <mg> (double) parameter \f$m/g\f$ of the Hamiltonian
    \param <x> (double) parameter \f$x\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <beta> (double) inverse temperature
    \param <delta> (double) width of time step (actually, this is the 
                   coeff. in the exponential, and the time step turns 
		   out to be 2*delta)
    \param <stepBlk> (int) save results every stepBlk steps (means
            every delta*stepBlk*2 (4) in beta)
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)

*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int N=atoi(argv[++cntr]);
  double mg=atof(argv[++cntr]);
  double x=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  double beta=atof(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int stepBlk=atoi(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  double mu=2*mg*sqrt(x);

  double alpha=.5;
  double L0P=0.;
  //double L0M=-1.;

  cout<<"Initialized arguments: L="<<N
      <<", mu="<<mu
      <<", x="<<x
      <<", Tmin="<<1/beta
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% N\t x\t mg\t D\t delta\t nrsteps\t beta\t T\t Energy"
	<<endl;
    *out<<setprecision(10);
  }
  else
    out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);

  int M=beta/(2*delta);

  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  cout<<"Initialized Contractor"<<endl;
  int d0=2;int Dint=d0;
  MPS thS(2*N,Dint,d0);
  //cout<<"Initialized MPS "<<thS<<endl;
  // Tensor for the even sites (left part of the max. entgld)
  mwArray Aeven(shrt::Indices(d0,Dint,Dint));
  // Tensor for the odd sites (right part of the max. entgld)
  mwArray Aodd(shrt::Indices(d0,Dint,Dint));
  //  double val=1./sqrt(sqrt(d0));
  double val=1.;
  for(int k=0;k<d0;k++){
    Aeven.setElement(val,0.,shrt::Indices(k,0,k));
    Aodd.setElement(val,0.,shrt::Indices(k,k,0));
  }
  // Now set the sites in the MPS
  for(int k=0;k<N;k++){
    if(k==0){
      mwArray firstSite(Aeven);
      firstSite.resize(shrt::Indices(d0,1,Dint));
      thS.setA(2*k,firstSite);
    }
    else
      thS.setA(2*k,Aeven);
    if(k==N-1){
      mwArray lastSite(Aodd);
      lastSite.resize(shrt::Indices(d0,Dint,1));
      thS.setA(2*k+1,lastSite);
    }
    else
      thS.setA(2*k+1,Aodd);
  }
  cout<<"Initialized identity state. Check norm="
      <<contractor.contract(thS,thS)<<endl;
  // I keep it for later (identity contraction)
  MPS idState(thS);
 
  double offset=0.;
  SchwingerHamiltonian hSchP(N,mu,x,L0P,alpha,offset);
  //SchwingerHamiltonian hSchM(N,mu,x,L0M,alpha);
  const MPO& hamilP=hSchP.getHMPO();
  cout<<"Initialized basic Hamiltonian"<<endl;

  complex_t deltaC=(complex_t){-delta,0.}; // delta  as complex
  complex_t deltaC_2=(complex_t){-delta*.5,0.}; // delta/2

  // Get the exponential MPOs
  MPO expHxe_2(N),expHxo(N),expHL(N),expHL_2(N),expHxe(N),expHxo_2(N);
  hSchP.getExponentialMPOx(expHxe_2,deltaC_2,true);  
  hSchP.getExponentialMPOx(expHxe,deltaC,true);  
  cout<<"Constructed exponential MPO for Hxe_2"<<endl;
  hSchP.getExponentialMPOx(expHxo,deltaC,false);
  hSchP.getExponentialMPOx(expHxo_2,deltaC_2,false);
  cout<<"Constructed exponential MPO for Hxo"<<endl;
  hSchP.getExponentialMPOL(expHL,deltaC);  
  cout<<"Constructed exponential MPO for HL"<<endl;
  hSchP.getExponentialMPOL(expHL_2,deltaC_2);
  cout<<"Constructed all exponential MPOs"<<endl;

  cout<<"expHxe is "<<expHxe<<endl;
  cout<<"expHxo is "<<expHxo<<endl;
  //exit(1);


  // Now we want to extend these ones to act on doubled MPS, with half
  // the sites auxiliary. On those the MPO has to be identity in
  // vertical tensor identity in horizontal direction, but the
  // horizontal one will have the dimension of the MPO in between the
  // physical sites (sometimes Dop=d*d, sometimes 1, sometimes DopL=3)

  // An identity operator for intermediate sites
  int Dop=d0*d0; // the dimension in between 
  int DopL=3;
  int DopH=5;
  mwArray identPhys=identityMatrix(d0);
  mwArray identD=identityMatrix(DopH);
  identPhys.reshape(shrt::Indices(d0*d0,1)); 
  identD.reshape(shrt::Indices(1,DopH*DopH)); 
  mwArray _idOp=identPhys*identD;
  _idOp.reshape(shrt::Indices(d0,d0,DopH,DopH));
  _idOp.permute(shrt::Indices(1,3,2,4));
  Operator idOpH(_idOp);
  // Cut down the bond dim for the expHx operators!
  _idOp.resize(shrt::Indices(d0,Dop,d0,Dop));
  Operator idOp(_idOp);
  // Cut down again for the intermediate in expHL (DopL=3)
  _idOp.resize(shrt::Indices(d0,DopL,d0,DopL));
  Operator idOpL(_idOp);
  // the last one is different, because Dop should be one now
  identPhys.reshape(shrt::Indices(d0,1,d0,1));
  Operator idOpI(identPhys);
  cout<<"Constructed the auxiliary identities"<<endl;

  // The MPO for H:
  MPO unfoldedH(2*N);
  for(int k=0;k<N;k++){
    //cout<<"Constructing unfoldedH at "<<k<<endl;
    unfoldedH.setOp(2*k,&hamilP.getOp(k));
    if(k<N-1)
      unfoldedH.setOp(2*k+1,&idOpH);
    else
      unfoldedH.setOp(2*k+1,&idOpI);
  }  
  cout<<"Checking our long MPO for H "
      <<contractor.contract(thS,unfoldedH,thS)/contractor.contract(thS,thS)
      <<endl;

  // Now combine the operators, to write
  // the expH X conj(expH) on the doubled system
  // joined pairs of operators
  MPO unexpHe_2T(2*N),unexpHoT(2*N),unexpHLT(2*N),unexpHL_2T(2*N);

  // operators on physical systems
  MPO unexpHe_2(2*N),unexpHo(2*N),unexpHL(2*N),unexpHL_2(2*N);
  // And the operators on the ancillas
  MPO unexpHe_2A(2*N),unexpHoA(2*N),unexpHLA(2*N),unexpHL_2A(2*N);
  // another trick: operators Hxo_2 and Hxe
  MPO unexpHe(2*N),unexpHo_2(2*N),unexpHeA(2*N),unexpHo_2A(2*N);

  for(int k=0;k<N;k++){
    //cout<<"Constructing expHxe_2 at 2*"<<k<<endl;
    unexpHe_2.setOp(2*k,&expHxe_2.getOp(k));
    unexpHo.setOp(2*k,&expHxo.getOp(k));
    unexpHL.setOp(2*k,&expHL.getOp(k));
    unexpHL_2.setOp(2*k,&expHL_2.getOp(k));
    // the extra ones
    unexpHe.setOp(2*k,&expHxe.getOp(k));
    unexpHo_2.setOp(2*k,&expHxo_2.getOp(k));
    // for the ancillas
    unexpHe_2A.setConjugatedOp(2*k+1,&expHxe_2.getOp(k));
    unexpHoA.setConjugatedOp(2*k+1,&expHxo.getOp(k));
    unexpHLA.setConjugatedOp(2*k+1,&expHL.getOp(k));
    unexpHL_2A.setConjugatedOp(2*k+1,&expHL_2.getOp(k));
    // the extra ones
    unexpHeA.setConjugatedOp(2*k+1,&expHxe.getOp(k));
    unexpHo_2A.setConjugatedOp(2*k+1,&expHxo_2.getOp(k));
  
    if(k!=N-1){// the last one is out the MPO
      unexpHL.setOp(2*k+1,&idOpL);
      unexpHL_2.setOp(2*k+1,&idOpL);
    }
    else{
      unexpHL.setOp(2*k+1,&idOpI);
      unexpHL_2.setOp(2*k+1,&idOpI);
    }
    if(k!=0){// the first one is out the ancillary MPO
      unexpHLA.setOp(2*k,&idOpL);
      unexpHL_2A.setOp(2*k,&idOpL);
    }
    else{
      unexpHLA.setOp(2*k,&idOpI);
      unexpHL_2A.setOp(2*k,&idOpI);
    }

    if(k%2==0){
      if(k!=N-1){
	unexpHe_2.setOp(2*k+1,&idOp);
	unexpHe.setOp(2*k+1,&idOp);// extras
      }
      else{
	unexpHe_2.setOp(2*k+1,&idOpI);
	unexpHe.setOp(2*k+1,&idOpI);// extras
      }
      unexpHo.setOp(2*k+1,&idOpI);
      unexpHo_2.setOp(2*k+1,&idOpI);// extras

      unexpHe_2A.setOp(2*k,&idOpI);
      unexpHeA.setOp(2*k,&idOpI);// extras
      if(k!=0){
	unexpHoA.setOp(2*k,&idOp);
	unexpHo_2A.setOp(2*k,&idOp);// extras
      }
      else{
	unexpHoA.setOp(2*k,&idOpI);
	unexpHo_2A.setOp(2*k,&idOpI);// extras
      }
   }
    else{
      unexpHe_2.setOp(2*k+1,&idOpI);
      unexpHe.setOp(2*k+1,&idOpI);// extras
      if(k!=N-1){
	unexpHo.setOp(2*k+1,&idOp);
	unexpHo_2.setOp(2*k+1,&idOp);// extras
      }
      else{
	unexpHo.setOp(2*k+1,&idOpI);
	unexpHo_2.setOp(2*k+1,&idOpI);// extras
      }

      unexpHe_2A.setOp(2*k,&idOp);
      unexpHoA.setOp(2*k,&idOpI);
      // extras
      unexpHeA.setOp(2*k,&idOp);// extras
      unexpHo_2A.setOp(2*k,&idOpI);// extras
    }
  }
  cout<<"Constructed unfolded exponential operators "<<endl;


  // Now join them
  const MPO* aux[2];
  cout<<"Joining expHxe_2"<<endl;
  aux[0]=&unexpHe_2A;aux[1]=&unexpHe_2; 
  MPO::join(2,aux,unexpHe_2T);
  cout<<"Joining expHxo"<<endl;
  aux[0]=&unexpHoA;aux[1]=&unexpHo; 
  MPO::join(2,aux,unexpHoT);
  /**
  cout<<"Joining expHxe_2"<<endl;
  aux[0]=&unexpHo_2A;aux[1]=&unexpHe_2; 
  MPO::join(2,aux,unexpHe_2T);
  cout<<"Joining expHxo"<<endl;
  aux[0]=&unexpHeA;aux[1]=&unexpHo; 
  MPO::join(2,aux,unexpHoT);
  */
  cout<<"Joining expHL"<<endl;
  aux[0]=&unexpHLA;aux[1]=&unexpHL; 
  MPO::join(2,aux,unexpHLT);
  cout<<"Joining expHL_2"<<endl;
  aux[0]=&unexpHL_2A;aux[1]=&unexpHL_2; 
  MPO::join(2,aux,unexpHL_2T);
  cout<<"Constructed joined exponential operators "<<endl;

  // Ready to compute things

  int cnt=0;
  double betai=0;
  while(cnt<M){
    MPS aux(thS);
    contractor.optimize(unexpHL_2T,aux,thS,D);
    //cout<<"Optimized unexpHL_2 step "<<cnt<<endl;
    int k=0;
    while((k<stepBlk-1)&&(cnt<M-1)){
      aux=thS;
      contractor.optimize(unexpHe_2T,aux,thS,D);
      //cout<<"Optimized unexpHe_2 step "<<cnt<<endl;
      aux=thS;
      contractor.optimize(unexpHoT,aux,thS,D);
      //cout<<"Optimized unexpHo step "<<cnt<<endl;
      aux=thS;
      contractor.optimize(unexpHe_2T,aux,thS,D);
      //cout<<"Optimized second unexpHe_2 step "<<cnt<<endl;
      aux=thS;
      contractor.optimize(unexpHLT,aux,thS,D);
      //cout<<"Optimized unexpHL step "<<cnt<<endl;
      betai+=2*delta;
      k++;cnt++;
    }
    aux=thS;
    contractor.optimize(unexpHe_2T,aux,thS,D);
    //cout<<"Optimized unexpHe_2 step "<<cnt<<endl;
    aux=thS;
    contractor.optimize(unexpHoT,aux,thS,D);
    //Cout<<"Optimized unexpHo step "<<cnt<<endl;
    aux=thS;
    contractor.optimize(unexpHe_2T,aux,thS,D);
    //cout<<"Optimized unexpHe_2 step "<<cnt<<endl;
    aux=thS;
    contractor.optimize(unexpHL_2T,aux,thS,D);
    //cout<<"Optimized unexpHL_2 step "<<cnt<<endl;
    betai+=2*delta;
    cnt++;
    // Now I need to normalize
    thS.gaugeCond('R',1);
    cout<<"Result normalized"<<endl;

    // Now save some results
    double norm1=real(contractor.contract(thS,idState)); //trace of rho
    double norm2=real(contractor.contract(thS,thS)); //trace of rho^2
    cout<<"Computed norm "<<norm1<<" and norm2="<<norm2<<endl;
    complex_t energy=contractor.contract(thS,unfoldedH,idState);
    //cout<<"Computed energy "<<energy<<endl;
    complex_t energy2=contractor.contract(thS,unfoldedH,thS);
    //cout<<"Computed energy2 "<<energy2<<endl;
    *out<<N<<"\t"<<x<<"\t"<<mg<<"\t"<<D<<"\t"
	<<delta<<"\t"<<cnt<<"\t"<<betai<<"\t"<<1/betai<<"\t"
	<<real(energy)/norm1<<endl;
    *out<<N<<"\t"<<x<<"\t"<<mg<<"\t"<<D<<"\t"
	<<delta<<"\t"<<2*cnt<<"\t"<<2*betai<<"\t"<<.5/betai<<"\t"
	<<real(energy2)/norm2<<endl;
    cout<<"At the end of step "<<cnt<<", 2*beta="<<2*betai<<", T="<<.5/betai
	<<", energy="<<real(energy2)/norm2
	<<", energy="<<energy/norm1<<", norm="<<norm1<<endl;

  }
    MPO testPos(N); // test positivity?
    for(int pos=0;pos<N;pos++){
      // Construct the operator for site pos
      cout<<"Constructing op for site "<<pos<<endl;
      mwArray site=thS.getA(2*pos).getA();
      cout<<"Phys site mwArray="<<site.getDimensions()<<endl;
      int d=site.getDimension(0);
      int Dl=site.getDimension(1);
      int Dr=site.getDimension(2);
      site.reshape(Indices(d*Dl,Dr));
      mwArray ancilla=thS.getA(2*pos+1).getA();
      cout<<"ancilla mwArray="<<ancilla.getDimensions()<<endl;
      int Dl2=thS.getA(2*pos+1).getDl();
      int Dr2=thS.getA(2*pos+1).getDr();
      ancilla.permute(Indices(2,1,3));ancilla.reshape(Indices(Dl2,d*Dr2));
      site=site*ancilla;
      site.reshape(Indices(d,Dl,d,Dr2));
      testPos.setOp(pos,new Operator(site),true);
    }
    // Now find the min eigenvalue of this guy
    double minValRho=0.;
    MPS valRho(N,D,2);valRho.setRandomState();
    contractor.findGroundState(testPos,D,&minValRho,valRho);
    cout<<"Found min eigenval of rho="<<minValRho<<endl;
}
