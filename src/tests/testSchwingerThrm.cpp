// - slightly chnaging to change distance automatically
//   from version of testSchwingerITTrtr_20131109.cpp
// - cahnging convergence without result of variational method 
//   Old version is putted to testSchwingerITTrtr_20131111.cpp
// - Thermal state version
#include <math.h>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string.h>

#include "MPS.h"
#include "Contractor.h"
#include "SchwingerHamiltonian.h"
#include "DoubleOperator.h"

#include <vector>
using namespace std;
using namespace shrt;

#define LOCALDIR "./"

/** SchwingerPD runs the findGroundState routine with the MPO for the
    Schwinger Hamiltonian \ref <SchwingerHamiltonian>, and saves
    values of order parameters, to construct a phase diagram. 

    To be used as test in the new implementation, compile make test9

    Receives arguments:
    \param <L> (int) length of the chain
    \param <mg> (double) parameter \f$m/g\f$ of the Hamiltonian
    \param <x> (double) parameter \f$x\f$ of the Hamiltonian
    \param <L0> (double) parameter \f$l_0\f$ of the Hamiltonian
    \param <alpha> (double) parameter \f$\alpha\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file
    \param <app> (int) whether the output file is to be kept (app==1)

    For it to work in Mac Os X, a thread has to be created which
    receives the main function as argument.
*/

int main(int argc,const char* argv[]){

  double delta[6];  int cntr=0;
  clock_t start, end;
  clock_t start2, end2;

  int L=atoi(argv[++cntr]);
  double mg=atof(argv[++cntr]);
  double x=atof(argv[++cntr]);
  double L0=atof(argv[++cntr]);
  double alpha=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  delta[0]=atof(argv[++cntr]);
  delta[1]=atof(argv[++cntr]);
  delta[2]=atof(argv[++cntr]);
  delta[3]=atof(argv[++cntr]);
  delta[4]=atof(argv[++cntr]);
  delta[5]=atof(argv[++cntr]);
  int init_MPS=atoi(argv[++cntr]);  // =0 :random start, !=0 :read start
  const char* outfname=argv[++cntr];
  const char* rdflnm_init=argv[++cntr];
  const char* bkdrctnm_init=argv[++cntr];
  int app=atoi(argv[++cntr]);

  double mu=2*mg*sqrt(x);
  int d=2; // physical dimension
  ofstream *out; ofstream *out_tm;
  ifstream *in_tm;
  char rdflnm[50]; char bkupflnm[50]; 
  char rdtmflnm[50]; char bkuptmflnm[50];

  out=new ofstream(outfname);

  cout<<"L : "<<L
      <<", mg : "<<mg
      <<", x : "<<x
      <<", L0 : "<<L0
      <<", alpha : "<<alpha
      <<", D : "<<D
      <<", M : "<<M
      <<", outfname : "<<outfname<<endl;

  cout<<"delta[0]="<<delta[0]<<endl;
  cout<<"delta[1]="<<delta[1]<<endl;
  cout<<"delta[2]="<<delta[2]<<endl;
  cout<<"delta[3]="<<delta[3]<<endl;
  cout<<"delta[4]="<<delta[4]<<endl;
  cout<<"delta[5]="<<delta[5]<<endl;

  *out<<"#L : "<<L
      <<", mg : "<<mg
      <<", x : "<<x
      <<", L0 : "<<L0
      <<", alpha : "<<alpha
      <<", D : "<<D
      <<", M : "<<M
      <<", outfname : "<<outfname<<endl;

  double E0=0.; double Norm=0.;
  double E0_init=0.; double Cond_init=0.;
  complex_t CondSch;
  double time_chngstp[5];

  Contractor& contractor = Contractor::theContractor();

  MPS gs(L, D, d);
  gs.setRandomState();
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);

  SchwingerHamiltonian hamSch(L, mu, x, L0, alpha);
  const MPO& hamil=hamSch.getHMPO();

  MPO CondMPO(L);
  hamSch.constructCondensateMPO(CondMPO);

  //=====================
  // GROUND STATE SEARCH
  //=====================
  if(init_MPS==0){ 
    start = clock();
    contractor.findGroundState(hamil, D, &E0_init, gs);
    Cond_init=real(contractor.contract(gs, CondMPO, gs));
    cout<<"Ground state :"<<scientific<<setprecision(12)<<E0_init<<endl;
    cout<<"Condensate :"<<scientific<<setprecision(12)<<Cond_init<<endl;
    end = clock();
    cout<<"TIMER OF VAR. : "
	<<(end-start)/(float)CLOCKS_PER_SEC<<endl;
  }

  gs.clear();

  //======================
  // TIME REVOLUTION PART
  //======================
  cout<<"=== TIME EVOLUTION ==="<<endl;
  if(M>0){
    start = clock();

    MPO evolxo(L), evolxe_half(L), evolL_half(L), evolL(L), Hdoub(L), Conddoub(L);  // For double opr
    MPO _evolxo(L), _evolxe_half(L), _evolL_half(L), _evolL(L);  // For single opr

    MPS imaGS(L,D,d*d);

    MPS Identity(L,D,d*d);
    Identity.setProductState(p_maxent); 
    Identity.gaugeCond('R',1);
    Identity.gaugeCond('L',1);

    //    double time_ckpnt=0.001;
    double time_ckpnt=0.0001;  // For N>=200
    double time=0.;    double time_ck=0.;
    //    double time_bk=0.0005;// N=50
    double time_bk=0.3;// N=20
    double distance;    double E0_old=1.;
    double Norm;
    
    int istp=0;
    int cnt=0;
    int init_HMPO=0;

    if(init_MPS==0){
      imaGS.setProductState(p_maxent); 
      imaGS.gaugeCond('R',1);
      imaGS.gaugeCond('L',1);
      cout<<"Random starting"<<endl;
    }else{
      strcpy(rdflnm, rdflnm_init);
      ostringstream ost;
      ost<<init_MPS;
      const char* flnm_tmp=ost.str().c_str();
      strcat(rdflnm, flnm_tmp);
      cout<<"Reading file : "<<rdflnm<<endl;
      imaGS.importMPS(rdflnm);  // For read start

      strcpy(rdtmflnm, rdflnm);
      strcat(rdtmflnm, ".sub");
      in_tm=new ifstream(rdtmflnm);
      *in_tm>>time;
      in_tm->close();
    }
      
    MPS aux;

    //    while(cnt<M && istp<6){
    while(cnt<M && istp<1){
      aux=imaGS;

      if(init_HMPO==0){
	start2 = clock();
	hamSch.getExponentialMPOx(_evolxe_half, -delta[istp]*.5*.5*ONE_c,true);
	hamSch.getExponentialMPOx(_evolxo, -delta[istp]*.5*ONE_c,false);
	hamSch.getExponentialMPOL(_evolL_half, -delta[istp]*.5*.5*ONE_c);
	hamSch.getExponentialMPOL(_evolL, -delta[istp]*.5*ONE_c);
	
	mwArray identPhys=identityMatrix(d);
	identPhys.reshape(Indices(d,1,d,1));
			  
	for(int k=0; k<L; k++){
	  mwArray aux=_evolxe_half.getOp(k).getFullData();
	  evolxe_half.setOp(k, new DoubleOperator(_evolxe_half.getOp(k).getFullData(),
						  permute(aux, Indices(3,2,1,4))), true);
	  aux=_evolxo.getOp(k).getFullData();
	  evolxo.setOp(k, new DoubleOperator(_evolxo.getOp(k).getFullData(),
					     permute(aux, Indices(3,2,1,4))), true);
	  aux=_evolL_half.getOp(k).getFullData();
	  evolL_half.setOp(k, new DoubleOperator(_evolL_half.getOp(k).getFullData(),
						 permute(aux, Indices(3,2,1,4))), true);      
	  aux=_evolL.getOp(k).getFullData();
	  evolL.setOp(k, new DoubleOperator(_evolL.getOp(k).getFullData(),
					    permute(aux, Indices(3,2,1,4))), true);      

	  Hdoub.setOp(k, new DoubleOperator(hamil.getOp(k).getFullData(),
					    identPhys), true);

	  Conddoub.setOp(k, new DoubleOperator(CondMPO.getOp(k).getFullData(),
					       identPhys), true);
	}

	end2 = clock();
	cout<<"TIMER OF Construct Double MPOs : "
	    <<(end2-start2)/(float)CLOCKS_PER_SEC<<endl;

	// For single opr
	_evolxo.clear();
	_evolxe_half.clear();
	_evolL_half.clear();
	_evolL.clear();
	
	//	hamil.clear();
	CondMPO.clear();
	
	init_HMPO++;
      }

      if(cnt==0){
	Norm=real(contractor.contract(imaGS, Identity));
	E0=real(contractor.contract(imaGS, Hdoub, Identity))/Norm;
	CondSch=contractor.contract(imaGS, Conddoub, Identity)/Norm;
	
	cout<<"Red MPS configuration OR Random start MPS"<<endl;
	cout<<"#"<<fixed<<time
	    <<"  "<<scientific<<setprecision(12)<<E0
	    <<"  "<<distance
	    <<"  "<<real(CondSch)
	    <<endl; 
      }

      contractor.optimize(evolL_half, aux, imaGS, D);
      //      while(time_ck<time_ckpnt){
      while(time_ck<time_ckpnt-delta[istp]){
	contractor.optimize(evolxe_half, imaGS, aux, D);
	contractor.optimize(evolxo, aux, imaGS, D);
	contractor.optimize(evolxe_half, imaGS, aux, D);
	contractor.optimize(evolL, aux, imaGS, D);
	time_ck+=delta[istp];  time+=delta[istp];
	cnt++;
      }

      contractor.optimize(evolxe_half, imaGS, aux, D);
      contractor.optimize(evolxo, aux, imaGS, D);
      contractor.optimize(evolxe_half, imaGS, aux, D);
      contractor.optimize(evolL_half, aux, imaGS, D);
      time_ck+=delta[istp];  time+=delta[istp];
      cnt++;

      imaGS.gaugeCond('R', true);
	
      Norm=real(contractor.contract(imaGS, Identity));
      E0=real(contractor.contract(imaGS, Hdoub, Identity))/Norm;
      CondSch=contractor.contract(imaGS, Conddoub, Identity)/Norm;
      
      if( cnt!=0 && E0_old==1.){ cout<<"Something wrong with E0"<<endl; }
      distance = abs(E0-E0_old)/abs(E0_old);
	
      if(cnt>0 && distance<1E-8){
	time_chngstp[istp]=time;
	istp++; 
	init_HMPO=0; 
      }

      cout<<fixed<<time
	  <<"  "<<scientific<<setprecision(12)<<E0
	  <<"  "<<distance
	  <<"  "<<real(CondSch)
	  <<endl; 
      *out<<fixed<<time
	  <<"  "<<scientific<<setprecision(12)<<E0
	  <<"  "<<distance
	  <<"  "<<real(CondSch)
	  <<endl;

      //      if(fmod(time, time_bk)<1e-13 && time<=1.2){
      if(time==0.3 || time==0.6 || time==0.9){
	//	cout<<"Bck up work at "<<time<<endl;
	strcpy(bkupflnm, bkdrctnm_init);
	ostringstream ost_b;
	ost_b<<time;
	const char* bflnm_tmp=ost_b.str().c_str();
	strcat(bkupflnm, bflnm_tmp);
	cout<<"Buckup file (intermid. time): "<<bkupflnm<<endl;
	imaGS.exportMPS(bkupflnm);
      }

      E0_old=E0;
      time_ck=0;
	
    }

    if(M==cnt){ cout<<"Maximum step achieved!! It should be checked!!"<<endl; }

    strcpy(bkupflnm, rdflnm_init);
    ostringstream ost_b;
    ost_b<<init_MPS+1;
    const char* bflnm_tmp=ost_b.str().c_str();
    strcat(bkupflnm, bflnm_tmp);
    cout<<"Buckup file : "<<bkupflnm<<endl;
    imaGS.exportMPS(bkupflnm);

    strcpy(bkuptmflnm, bkupflnm);
    strcat(bkuptmflnm, ".sub");
    out_tm=new ofstream(bkuptmflnm);
    *out_tm<<time<<endl;
    out_tm->close();

    end = clock();
    cout<<"TIMER OF TIME EVOL.: "
        <<(end-start)/(float)CLOCKS_PER_SEC<<endl;

  }

  //  cout<<"==== Time STP chang ====="<<endl;
  //  cout<<"const0="<<fixed<<time_chngstp[0]<<endl;
  //  cout<<"const1="<<fixed<<time_chngstp[1]<<endl;
  //  cout<<"const2="<<fixed<<time_chngstp[2]<<endl;
  //  cout<<"const3="<<fixed<<time_chngstp[3]<<endl;
  //  cout<<"const4="<<fixed<<time_chngstp[4]<<endl;
  //  cout<<"const5="<<fixed<<time_chngstp[5]<<endl;


  out->close();

  return 0;

}

