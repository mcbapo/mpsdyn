
#include <math.h>
#include <iomanip>

#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "HeisenbergHamiltonian.h"

using namespace shrt;

/** testHeisenberg runs the findGroundState routine with the MPO for the
    Heisenberg Hamiltonian \ref <HeisenberggHamiltonian>, to check convergence
    and then runs imaginary time evolution to check also the
    implementation of the unitary evolution.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <Jx> (double) parameter \f$Jx\f$ of the Hamiltonian
    \param <Jy> (double) parameter \f$Jy\f$ of the Hamiltonian
    \param <Jz> (double) parameter \f$Jz\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <M> (int) maximum number of imaginary time steps 
           (if 0, no imaginary time phase is run)
    \param <delta> (double) step width for imaginary time evolution
    \param <outfname> (char*) name of the output file for the energy
    \param <app> (int) whether the output file is to be kept (app==1)
*/


int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  int L=atoi(argv[++cntr]);
  double Jx_=atof(argv[++cntr]);
  double Jy_=atof(argv[++cntr]);
  double Jz_=atof(argv[++cntr]);
  double h_=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  const char* outfname=argv[++cntr];
  int app=atoi(argv[++cntr]);

  cout<<"Initialized arguments: L="<<L
      <<", Jx="<<Jx_
      <<", Jy="<<Jy_
      <<", Jz="<<Jz_
      <<", h="<<h_
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname);
    *out<<"% N\t Jx\t Jy\t Jz\t h\t D\t Energy\t Energy/N"<<endl;
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
   MPS gs(L,D,d);
   gs.setRandomState(); // the intial state, random
   gs.gaugeCond('R',1);
   gs.gaugeCond('L',1);
   cout<<"Initialized random state, norm "<<contractor.contract(gs,gs)<<endl;
   cout<<setprecision(10);
   
   // First: put H and find the GS
   //int Dop=5; // bond dimension of Heisenbergg H
   double E0=0.;
   // Set random coefficients
   vector<double> Jx(L-1,Jx_);
   vector<double> Jy(L-1,Jy_);
   vector<double> Jz(L-1,Jz_);
   vector<double> h(L,h_);
   for(int k=0;k<L;k++){
     if(k<L-1){
       Jx[k]=(2.*rand()/RAND_MAX-1.)*Jx_;
       Jy[k]=(2.*rand()/RAND_MAX-1.)*Jy_;
       Jz[k]=(2.*rand()/RAND_MAX-1.)*Jz_;
       // Jx[k]=Jx_;
       // Jy[k]=Jy_;
       // Jz[k]=Jz_;
     }
     h[k]=(2.*rand()/RAND_MAX-1.)*h_;
     // h[k]=h_;
   }

   HeisenbergHamiltonian hamH(L,Jx,Jy,Jz,h,d);
   cout<<"Created the Hamiltonian"<<endl;
   const MPO& hamil=hamH.getHMPO();
   cout<<"Constructed the hamil MPO"<<endl;
   hamil.exportForMatlab("heisenbergHmpo.m");
   cout<<"Initial value, with initial state"<<endl;
   cout<<contractor.contract(gs,hamil,gs)<<endl;

   contractor.findGroundState(hamil,D,&E0,gs);
   cout<<"Ground state found with eigenvalue "<<E0
       <<" (Energy per particle="<<E0/L<<")"<<endl;
   // mwArray Amid=gs.getA(L/2).getA();
   // ofstream Amid_("Amid.dat");
   // Amid.savetext(Amid_);
   // Amid_.close();
   // //   putForMatlab(cout,Amid,"Amid");
   // cout<<endl;

   *out<<setprecision(10);
   *out<<"Jx=[";
   for(int k=0;k<L-1;k++) *out<<Jx[k]<<" ";
   *out<<"];"<<endl;
   *out<<"Jy=[";
   for(int k=0;k<L-1;k++) *out<<Jy[k]<<" ";
   *out<<"];"<<endl;
   *out<<"Jz=[";
   for(int k=0;k<L-1;k++) *out<<Jz[k]<<" ";
   *out<<"];"<<endl;
   *out<<"h=[";
   for(int k=0;k<L;k++) *out<<h[k]<<" ";
   *out<<"];"<<endl;

   *out<<L<<"\t"
       <<D<<"\t"<<E0<<"\t"<<E0/L<<"\t";
   *out<<endl;

   if(M>0){
     contractor.setConvTol(1E-8);
     // Now try evolution with the exponential, instead
     MPO expHe(L);
     MPO expHo(L);
     MPO expHe2(L);
     int oT=2; // Trotter order
     hamH.getExponentialMPOeven(expHe,-delta*ONE_c);
     expHe.exportForMatlab("heisenbergExpHe.m");
     hamH.getExponentialMPOeven(expHe2,-.5*delta*ONE_c);
     expHe2.exportForMatlab("heisenbergExpHe2.m");
     hamH.getExponentialMPOodd(expHo,-delta*ONE_c);
     expHo.exportForMatlab("heisenbergExpHo.m");
     *out<<endl;
     *out<<"% Imaginary time evolution!! "<<endl;
     *out<<"% t\t Energy\t Energy/N"<<endl;

     MPS imaGS(L,D,d);
     //     imaGS.setRandomState(); imaGS.gaugeCond('R',1);
     imaGS.exportForMatlab("initialV.m");
     double oldE=real(contractor.contract(imaGS,hamil,imaGS));
     double time=0.;
     int cnt=0;
     bool conver=0;
     while(cnt<M&&!conver){
       time+=delta;
       cout<<"Time evolution step nr "<<cnt+1<<", time="<<time
	   <<", E0="<<oldE<<endl;
       {MPS aux(imaGS); // temporary copy
	 //	 contractor.optimize(expHe2,aux,imaGS,D);
	 contractor.optimize(expHe2,aux,imaGS,D);
       }
       // cout<<"After expHe2, v'*H*v="<<contractor.contract(imaGS,hamil,imaGS)
       // 	   <<", v*v="<<contractor.contract(imaGS,imaGS)<<endl;
       {MPS aux(imaGS); 
	 contractor.optimize(expHo,aux,imaGS,D);
       }
       // cout<<"After expHo, v'*H*v="<<contractor.contract(imaGS,hamil,imaGS)
       // 	   <<", v*v="<<contractor.contract(imaGS,imaGS)<<endl;
       {MPS aux(imaGS); 
         contractor.optimize(expHe2,aux,imaGS,D);
       }
       // cout<<"After expHe2, v'*H*v="<<contractor.contract(imaGS,hamil,imaGS)
       // 	   <<", v*v="<<contractor.contract(imaGS,imaGS)<<endl;
       imaGS.gaugeCond('R',true);
       //cout<<"After normalizing, v'*v="<<contractor.contract(imaGS,imaGS)<<endl;
       E0=real(contractor.contract(imaGS,hamil,imaGS));
       
       
       *out<<time<<"\t"<<E0<<"\t"<<E0/L<<"\t";
       *out<<endl;
       if(cnt>0&&abs((E0-oldE)/oldE)<1E-8) conver=true;
       cnt++;
       oldE=E0;
     }

   //  evol.exportForMatlab("isingUmpo.m");
   // hamil.exportForMatlab("isingHmpo.m");
   }

   out->close();
}
