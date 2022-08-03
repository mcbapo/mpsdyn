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

#include "misc.h"
#include "MPS.h"
#include "Contractor.h"
#include "SchwingerHamiltonian.h"
#include "DoubleOperator.h"

#include <vector>
using namespace std;
using namespace shrt;

#define LOCALDIR "./"

int main(int argc,const char* argv[]){

  int cntr=0; 
  clock_t start, end;
  clock_t start2, end2;

  int L=atoi(argv[++cntr]);
  double mg=atof(argv[++cntr]);
  double x=atof(argv[++cntr]);
  double L0=atof(argv[++cntr]);
  double alpha=atof(argv[++cntr]);
  int D=atoi(argv[++cntr]);
  int M=atoi(argv[++cntr]);
  double delta=atof(argv[++cntr]);
  int init_MPS=atoi(argv[++cntr]);  // =0 :random start, !=0 :read start
  const char* outfname=argv[++cntr];
  const char* rdflnm_init=argv[++cntr];
  int Intrvl_bckup=atoi(argv[++cntr]);   
  double distanceBeta=atof(argv[++cntr]);  // For a case of large delta
  int app=atoi(argv[++cntr]);

  double mu=2*mg*sqrt(x);
  int d=2; // physical dimension
  ofstream *out; ofstream *out_tm;
  ifstream *in_tm;
  char rdflnm[150]; char bkupflnm[150]; 
  char rdtmflnm[150]; char bkuptmflnm[150];
  string lastTmpBackup,lastTmpBackupT;

  if(!app){
    out=new ofstream(outfname);
    *out<<"#L : "<<L
	<<", mg : "<<mg
	<<", x : "<<x
	<<", L0 : "<<L0
	<<", alpha : "<<alpha
	<<", delta : "<<delta
	<<", D : "<<D
	<<", M : "<<M
	<<", outfname : "<<outfname<<endl;
    out->close();
  }

  cout<<"L : "<<L
      <<", mg : "<<mg
      <<", x : "<<x
      <<", L0 : "<<L0
      <<", alpha : "<<alpha
      <<", D : "<<D
      <<", M : "<<M
      <<", outfname : "<<outfname<<endl;

  cout<<"delta="<<delta<<endl;
  cout<<"Intrvl_bckup="<<Intrvl_bckup<<endl;
  cout<<"distanceBeta="<<distanceBeta<<endl;


  double E0=0.; double Norm=0.;
  double E0_init=0.; double Cond_init=0.;
  complex_t CondSch;
  double time_chngstp[5];

  Contractor& contractor = Contractor::theContractor();

  // MPS gs(L, D, d);
  // gs.setRandomState();
  // gs.gaugeCond('R',1);
  // gs.gaugeCond('L',1);

  SchwingerHamiltonian hamSch(L, mu, x, L0, alpha);
  const MPO& hamil=hamSch.getHMPO();

  MPO CondMPO(L);
  hamSch.constructCondensateMPO(CondMPO);

  //======================
  // TIME REVOLUTION PART
  //======================
  cout<<"=== TIME EVOLUTION ==="<<endl;
  if(M>0){
    start = clock();

    MPO evolxo(L), evolxe(L), evolxe_half(L), evolL_half_S(L),evolL_half_A(L), Hdoub(L), Conddoub(L);  // For double opr
    MPO _evolxo(L), _evolxe(L), _evolxe_half(L), _evolL_half(L);  // For single opr

    MPS imaGS(L,D,d*d);

    MPS Identity(L,D,d*d);
    Identity.setProductState(p_maxent); 
    Identity.gaugeCond('R',1);
    Identity.gaugeCond('L',1);

    //double time_ckpnt=min(10*delta,.5*distanceBeta/(2*sqrt(x)));  // H.S. 2015.02.05
    double time_ckpnt=min(10*delta,distanceBeta/(2*sqrt(x)));  // H.S. 2015.02.05
    cout<<"time_ckpnt = "<<time_ckpnt<<endl;
    if(time_ckpnt==10*delta){
      cout<<"  (10*delta)"<<endl;
    }else{
      cout<<"  (distanceBeta/(2*sqrt(x)))"<<endl;
    }

    double Intrvl_comp=time_ckpnt/(2*delta);
    int Intrvl_comp_int=(int)Intrvl_comp;
    if( Intrvl_comp_int<Intrvl_comp ){ Intrvl_comp=Intrvl_comp+1.; }

    if( Intrvl_comp > Intrvl_bckup ){
      cout<<"  You're ordering backup files more than possible. Check the valuse of time_ckpnt, Intrvl_bckup etc."<<endl;
      exit(3);
    }

    double time=0.;    double time_ck=0.;
    //    double time_bk=0.0005;// N=50
    double time_bk=0.3;// N=20
    double distance;    double E0_old=1.;
    double Norm;
    
    int cnt=0;
    int cnt_bckup=1;
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
      cout<<"Imported MPS "<<imaGS<<"from file "<<rdflnm<<" for time "<<time<<endl;
    }
      
    MPS aux;

    while(cnt<M){
      aux=imaGS;

      if(init_HMPO==0){
	start2 = clock();
	hamSch.getExponentialMPOx(_evolxe_half, -delta*.5*.5*ONE_c,true);
	hamSch.getExponentialMPOx(_evolxe, -delta*.5*ONE_c,true);
	hamSch.getExponentialMPOx(_evolxo, -delta*.5*ONE_c,false);
	hamSch.getExponentialMPOL(_evolL_half, -delta*.5*.5*ONE_c);

	extendMPO(_evolL_half,evolL_half_S,d);
	extendTransposeMPO(_evolL_half,evolL_half_A,d);
	
	mwArray identPhys=identityMatrix(d);
	identPhys.reshape(Indices(d,1,d,1));
			  
	for(int k=0; k<L; k++){
	  mwArray aux=_evolxe_half.getOp(k).getFullData();
	  evolxe_half.setOp(k, new DoubleOperator(_evolxe_half.getOp(k).getFullData(),
						  permute(aux, Indices(3,2,1,4))), true);
	  aux=_evolxo.getOp(k).getFullData();
	  evolxo.setOp(k, new DoubleOperator(_evolxo.getOp(k).getFullData(),
					     permute(aux, Indices(3,2,1,4))), true);
	  //aux=_evolL_half.getOp(k).getFullData();
	  //evolL_half.setOp(k, new DoubleOperator(_evolL_half.getOp(k).getFullData(),
	  //					 permute(aux, Indices(3,2,1,4))), true);      

	  aux=_evolxe.getOp(k).getFullData();
	  evolxe.setOp(k, new DoubleOperator(_evolxe.getOp(k).getFullData(),
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
	_evolxe.clear();
	CondMPO.clear();
	
	init_HMPO++;
      }

      if(cnt==0){
	Norm=real(contractor.contract(imaGS, imaGS));
	E0=real(contractor.contract(imaGS, Hdoub, imaGS))/Norm;
	CondSch=contractor.contract(imaGS, Conddoub, imaGS)/Norm;

	distance = abs(E0-E0_old)/abs(E0_old);
	
	cout<<"Read MPS configuration OR Random start MPS"<<endl;
	cout<<"#"<<fixed<<time
	    <<"  "<<scientific<<setprecision(12)<<E0
	    <<"  "<<distance
	    <<"  "<<real(CondSch)
	    <<endl; 
      }

      contractor.optimize(evolxe_half, aux, imaGS, D);

      while(time_ck<time_ckpnt-2*delta){
	//	contractor.optimize(evolL_half, imaGS, aux, D);
	//	contractor.optimize(evolxo, aux, imaGS, D);
	//	contractor.optimize(evolL_half, imaGS, aux, D);
	//	contractor.optimize(evolxe, aux, imaGS, D);
	contractor.optimize(evolL_half_S, imaGS, aux, D);
	contractor.optimize(evolL_half_A, aux, imaGS, D);
	contractor.optimize(evolxo, imaGS,aux, D);
	contractor.optimize(evolL_half_S, aux,imaGS, D);
	contractor.optimize(evolL_half_A, imaGS, aux, D);
	contractor.optimize(evolxe, aux, imaGS, D);
	time_ck+=2*delta;  time+=2*delta;
	cnt++;
      }

      // contractor.optimize(evolL_half, imaGS, aux, D);
      // contractor.optimize(evolxo, aux, imaGS, D);
      // contractor.optimize(evolL_half, imaGS, aux, D);
      // contractor.optimize(evolxe_half, aux, imaGS, D);
      contractor.optimize(evolL_half_S, imaGS, aux, D);
      contractor.optimize(evolL_half_A, aux, imaGS, D);
      contractor.optimize(evolxo, imaGS,aux, D);
      contractor.optimize(evolL_half_S, aux,imaGS, D);
      contractor.optimize(evolL_half_A, imaGS, aux, D);
      contractor.optimize(evolxe_half, aux, imaGS, D);
      time_ck+=2*delta;  time+=2*delta;
      cnt++;
      // end2 = clock();
      // cout<<"TIMER OF "<<cnt-cnt0 <<" steps : "
      // 	  <<(end2-start2)/(float)CLOCKS_PER_SEC<<endl;

      imaGS.gaugeCond('R', true);

      Norm=real(contractor.contract(imaGS, imaGS));
      E0=real(contractor.contract(imaGS, Hdoub, imaGS))/Norm;
      CondSch=contractor.contract(imaGS, Conddoub, imaGS)/Norm;
      	
      distance = abs(E0-E0_old)/abs(E0_old);

      cout<<fixed<<time
	  <<"  "<<scientific<<setprecision(12)<<E0
	  <<"  "<<distance
	  <<"  "<<real(CondSch)
	  <<endl; 
      out=new ofstream(outfname,ios::app);
      *out<<fixed<<time
	  <<"  "<<scientific<<setprecision(12)<<E0
	  <<"  "<<distance
	  <<"  "<<real(CondSch)
	  <<endl;
      out->close();
      // Saving MPS Backup files 
      /*
      //  time dependent backup generations
      if(time==0.3 || time==0.6 || time==0.9){
        cout<<"Bck up work at "<<time<<endl;
	strcpy(bkupflnm, bkdrctnm_init);
	ostringstream ost_b;
	ost_b<<time;
	const char* bflnm_tmp=ost_b.str().c_str();
	strcat(bkupflnm, bflnm_tmp);
	cout<<"Buckup file (intermid. time): "<<bkupflnm<<endl;
	imaGS.exportMPS(bkupflnm);
      }
      */

      // iteration dependent backup generations
      // cnt_bckup = 1, 2, ..., M/Intrvl_bckup-1
      if( cnt>=Intrvl_bckup*cnt_bckup && cnt_bckup*Intrvl_bckup<M ){
        cout<<"  Generating Backup file at "<<time<<endl;
	strcpy(bkupflnm, rdflnm_init);

	ostringstream ost_b;
	ost_b<<init_MPS;
	const char* bflnm_tmp=ost_b.str().c_str();
	strcat(bkupflnm, bflnm_tmp);
       
	ostringstream ost_bckupID;
	ost_bckupID<<cnt_bckup;
	const char* bflnm_tmp2=ost_bckupID.str().c_str();
	strcat(bkupflnm, "-");
	strcat(bkupflnm, bflnm_tmp2);

	cout<<"  Backup file (intermid. time): "<<bkupflnm<<endl;
	imaGS.exportMPS(bkupflnm);

	// To Print out beta
	strcpy(bkuptmflnm, bkupflnm);
	strcat(bkuptmflnm, ".sub");
	out_tm=new ofstream(bkuptmflnm);
	*out_tm<<time<<endl;
	out_tm->close();

	// Remove the last ones, if they exist
	if(lastTmpBackup.length()>0){
	  if(file_exists(lastTmpBackup.data()))
	    remove(lastTmpBackup.data());
	  if(file_exists(lastTmpBackupT.data()))
	    remove(lastTmpBackupT.data());
	  cout<<"Removed old tmp files "<<lastTmpBackup<<" and "<<lastTmpBackupT<<endl;
	}
	// And store these names for next iteration
	lastTmpBackup=string(bkupflnm);
	lastTmpBackupT=string(bkuptmflnm);
	cnt_bckup++;
      }

      E0_old=E0;
      time_ck=0;

      if(real(CondSch)<0){ exit(2); }
	
    }

    //    if(M==cnt){ cout<<"Maximum step achieved!! It should be checked!!"<<endl; }

    strcpy(bkupflnm, rdflnm_init);
    ostringstream ost_b;
    ost_b<<init_MPS+1;
    const char* bflnm_tmp=ost_b.str().c_str();
    strcat(bkupflnm, bflnm_tmp);
    cout<<"Backup file : "<<bkupflnm<<endl;
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

