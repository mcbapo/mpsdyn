
#include <math.h>
#include <iomanip>

#include "Properties.h"
#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "SpinMPO.h"
#include "misc.h"


#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "ThirringHamiltonian.h"

using namespace shrt;

/** thirringP runs the findGroundState routine with the MPO for the
    Thirring Hamiltonian \ref <ThirringHamiltonian> with
    possibly a penalty term for a targetted total \f$S_z=S_{\mathrm{tot}}\f$.
    At the end, the total \f$S_z\f$ is computed.
    This time I group sites in pairs.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <ma> (double) adimensional mass parameter \f$ma\f$ of the Hamiltonian
    \param <g2> (double) adimensional coupling parameter \f$g^2\f$ of the Hamiltonian
    \param <gt2> (double) adimensional coupling parameter \f$\tilde{g}^2\f$
    of the Hamiltonian, to allow for a different coupling for even links
    \param <lambda> (double) penalty term \f$\lambda\f$ (\f$>0\f$)
    \param <mu> (double) chemical potential \f$\mu\f$ 
    \param <Starget> (int) targetted total \f$S_z\f$
    \param <D> (int) maximum bond dimension
    \param <outfname> (char*) name of the output file for the energy
    \param <app> (int) whether the output file is to be kept (app==1)
*/

int findTmpFile(int D,const string& mpsdir,const string& baseName,string& mpsfile);
const string mpsFileName(const string& mpsdir,const string& baseName,int D);

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  int L = props.getIntProperty("L");  
  if(L%2!=0){
    cout<<"ERROR: Nr of site L should be even"<<endl;
    exit(1);
  }
  double ma=props.getDoubleProperty("ma");  
  double gt2_=props.getDoubleProperty("gtilde2");
  double g2_=0; //props.getDoubleProperty("gtilde2"); // the non-alternating one
  double lambda_=props.getDoubleProperty("lambda");
  if(lambda_<0) lambda_=0.;
  double mu_=props.getDoubleProperty("mu");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  int Stot_=props.getIntProperty("Starget");
  int D=props.getIntProperty("D");
  int kmax=props.getIntProperty("maxlevel");
  if(kmax<0) kmax=0;
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir");
  string mpsfile=props.getProperty("mpsfile");
  string outfnameS=props.getProperty("outputSz"); // for entropies and polarizations
  bool saveS=(outfnameS.length()>0);
  bool app=(props.getIntProperty("append")!=0);

  cout<<"Initialized arguments: L="<<L
      <<", ma="<<ma
      <<", g2="<<g2_
      <<", gt2="<<gt2_
      <<", lambda="<<lambda_
      <<", mu="<<mu_
      <<", Starget="<<Stot_
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  cout<<"Initialized Contractor"<<endl;
  contractor.setConvTol(tol);
  int d=2;
  MPS gs(L/2,D,d*d);
  // Check if initial state was given
  string tmpMPSfile;
  if(findTmpFile(D,mpsdir,mpsfile,tmpMPSfile)){// found some tmp or smaller D file
    cout<<"Recovered initial guess from file "<<tmpMPSfile<<endl;
    if(gs.getLength()!=L/2){
      cout<<"Error! File "<<tmpMPSfile<<" does not contain valid MPS "<<endl;
      exit(1);
    }
    gs.importMPS(tmpMPSfile.data());
  }
  else{
    srandom(time(NULL));
    gs.setRandomState(); // the intial state, random
    cout<<"Initialized random state, norm "<<contractor.contract(gs,gs)<<endl;
  }
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  cout<<setprecision(10);
  
  // First: construct H
  double E0=0.;
  
  MPO SzmpoS(L);
  SpinMPO::getSzMPO(L,d,SzmpoS);
  MPO Szmpo(L/2);
 
  ThirringHamiltonian hamH(L,ma,g2_,lambda_,Stot_,d,gt2_,mu_);
  //cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamilS=hamH.getHMPO();
  MPO hamil(L/2);
  for(int k=0;k<L;k+=2){
    Operator* opK=new Operator(hamilS.getOp(k).getFullData());
    opK->contractRight(&hamilS.getOp(k+1));
    hamil.setOp(k/2,opK,true);
    Operator* opSz=new Operator(SzmpoS.getOp(k).getFullData());
    opSz->contractRight(&SzmpoS.getOp(k+1));
    Szmpo.setOp(k/2,opSz,true);
  }

  // if(L<12)
  //   hamil.exportForMatlab("thirringHmpo.m");
  //cout<<"Initial E value, with initial random state"<<endl;
  cout<<"E0="<<contractor.contract(gs,hamil,gs)<<endl;
  //cout<<"Initial Sz value, with initial random state"<<endl;
  cout<<"Sz="<<contractor.contract(gs,Szmpo,gs)<<endl;

  //  contractor.findGroundState(hamil,2,&E0,gs); // start sweeping with D=2

  // for(int pp=0;pp<3;pp++){
  //   for(int pos=0;pos<(L-1)/2+1;pos++){
  //     contractor.sweepPart(pos,pos+1,hamil,D,gs); // one right
  //     contractor.sweepPart(L-1-pos,L-1-pos-1,hamil,D,gs); // one left
  //   }
  //   double lambda=real(contractor.contract(gs,hamil,gs));
  //   cout<<"After symmetric sweep nr "<<pp+1<<"energy value "<<lambda<<endl;
  // }

  contractor.findGroundState(hamil,D,&E0,gs);
  // Save it!
  tmpMPSfile=mpsFileName(mpsdir,mpsfile,D);
  gs.exportMPS(tmpMPSfile.data());
  cout<<"Ground state saved to "<<tmpMPSfile<<endl;
  
  complex_t Sztot=contractor.contract(gs,Szmpo,gs);
  cout<<"Ground state found for S_tot="<<Stot_<<" with energy "<<E0
      <<" (Energy per particle="<<E0/L<<", plus offset from David->"
      <<E0/L+gt2_*.25+g2_*(L-1.)/(4.*L)<<")"
      <<" Total Sz="<<real(Sztot)<<endl;
  complex_t Sz2=contractor.contract2(Szmpo,gs);
  complex_t H2=contractor.contract2(hamil,gs);
  cout<<"DeltaH="<<real(H2)-E0*E0<<endl;
  cout<<"Delta S_z="<<Sz2-Sztot*Sztot<<endl;

  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    *out<<"% N\t ma\t gt2\t Starget\t lambda\t D\t Energy\t Delta H \t S_z\t Delta Sz\t E/L(David)"<<endl;
  }
  else
    out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(15);
  *out<<L<<"\t"<<ma<<"\t"<<gt2_<<"\t"<<Stot_<<"\t"<<lambda_<<"\t"
      <<D<<"\t 0 \t"<<E0<<"\t"<<real(H2)-E0*E0<<"\t"<<real(Sztot)<<"\t"<<real(Sz2-Sztot*Sztot)
      <<"\t"<<E0/L+gt2_*.25+g2_*(L-1.)/(4.*L);
  *out<<endl;

  ofstream* outS;
  if(saveS){
    if(!app){
      outS=new ofstream(outfnameS.data());
      *outS<<"% N\t ma\t g2\t D\t level\t Sz\t pos\t S_z\t S(0-pos:rest)"<<endl;
    }
    else
      outS=new ofstream(outfnameS.data(),ios::app);
    if(!outS->is_open()){
      cout<<"Error: impossible to open file "<<outfnameS<<
  	" for output"<<endl;
      exit(1);
    }
    *outS<<setprecision(15);
    // Compute entropies and Sz of all sites /cuts

    MPS aux(gs);aux.gaugeCond('R',1);
    complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
    mwArray sigz(Indices(d,d),dataZ); // sigmaZ operator on one site
    mwArray sig0=identityMatrix(d); // id
    // sigz on the left
    mwArray sigzL(sigz);sigzL.reshape(Indices(d*d,1));
    sigzL.multiplyRight(reshape(sig0,Indices(1,d*d)));
    sigzL.reshape(Indices(d,d,d,d));
    mwArray sigzR(sigzL);
    sigzL.permute(Indices(1,3,2,4));
    sigzL.reshape(Indices(d*d,d*d));
    sigzR.permute(Indices(3,1,4,2));
    sigzR.reshape(Indices(d*d,d*d));
    for(int pos=0;pos<L/2;pos++){
      aux.applyLocalOperator(pos,sigzL);
      complex_t valZl=.5*contractor.contract(aux,gs);
      // undo the change
      aux.applyLocalOperator(pos,sigzL);
      *outS<<L<<"\t"<<ma<<"\t"<<gt2_<<"\t"<<D<<"\t"
  	   <<0<<"\t"<<real(Sztot)<<"\t"
  	   <<2*pos<<"\t"<<real(valZl)<<"\t"<<0<<endl;      
      aux.applyLocalOperator(pos,sigzR);
      complex_t valZr=.5*contractor.contract(aux,gs);
      // undo the change
      aux.applyLocalOperator(pos,sigzR);
      // And the entropy
      double Spos=contractor.getEntropy(gs,pos+1);
      *outS<<L<<"\t"<<ma<<"\t"<<gt2_<<"\t"<<D<<"\t"
  	   <<0<<"\t"<<real(Sztot)<<"\t"
  	   <<2*pos+1<<"\t"<<real(valZr)<<"\t"<<Spos<<endl;      
    }
  }
  
  // vector<MPS*> levels;
  // levels.push_back(&gs); 
  // int k=1;
  // while(k<=kmax){
  //   cout<<"Computing now next excited state"<<endl;
  //   double E1=0.;
  //   MPS* exc=new MPS(L,D,d);
  //   exc->setRandomState(); // the initial state, random
  //   contractor.findNextExcitedState(hamil,D,levels,&E1,*exc);
  //   Sztot=contractor.contract(*exc,Szmpo,*exc);
  //   levels.push_back(exc); 
  //   Sz2=contractor.contract2(Szmpo,*exc);
  //   H2=contractor.contract2(hamil,*exc);
    
  //   cout<<k<<"-th excited state state found with eigenvalue "<<E1
  // 	<<" (Gap="<<E1-E0<<")"
  // 	<<" Total Sz="<<real(Sztot)<<endl;
  //   *out<<L<<"\t"<<"\t"<<ma<<"\t"<<gt2_<<"\t"<<Stot_<<"\t"<<lambda_<<"\t"
  // 	<<D<<"\t"<<k<<"\t"<<E1<<"\t"<<real(H2)-E1*E1<<"\t"<<real(Sztot)<<"\t"<<real(Sz2)-real(Sztot)*real(Sztot);
  //   *out<<endl;
  //   if(saveS){
  //     MPS aux(gs);aux.gaugeCond('R',1);
  //     complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
  //     mwArray sigz(Indices(d,d),dataZ); // sigmaZ operator on one site
  //     for(int pos=0;pos<L;pos++){
  // 	aux.applyLocalOperator(pos,sigz,0);
  // 	complex_t valZ=.5*contractor.contract(aux,gs);
  // 	// undo the change
  // 	aux.applyLocalOperator(pos,sigz,0);
  // 	// And the entropy
  // 	double Spos=contractor.getEntropy(gs,pos);
  // 	*outS<<L<<"\t"<<ma<<"\t"<<gt2_<<"\t"<<D<<"\t"
  // 	     <<k<<"\t"<<real(Sztot)<<"\t"
  // 	     <<pos<<"\t"<<real(valZ)<<"\t"<<Spos<<endl;      
  //     }
  //   }
  //   k++;
  // }
  
  out->close();
  delete out;
  if(saveS){
    outS->close();
    delete outS;
  }
}


const string mpsFileName(const string& mpsdir,const string& baseName,int D){
  stringstream s;
  s<<mpsdir<<"/"<<baseName;
  if(D>0) s<<"_"<<D;
  return s.str();
}

int findTmpFile(int D,const string& mpsdir,const string& baseName,string& mpsfile){
  mpsfile=mpsFileName(mpsdir,baseName,D)+"_tmp";
  bool found=0;
  if(file_exists(mpsfile)){
    cout<<"findTmpFile found mpsfile="<<mpsfile<<endl;
    found=true;
    return 1;
  }
  else{// search for smaller Ds
    int cnt=D;
    bool found=0;
    while(cnt>0&&!found){
      mpsfile=mpsFileName(mpsdir,baseName,cnt);
      if(file_exists(mpsfile)){
	cout<<"findTmpFile found mpsfile="<<mpsfile<<endl;
	found=true;
	return cnt;
      }
      else{
	mpsfile="";
	cnt--;
      }
    }
  }
  if(!found) return 0;
  cout<<"ERROR: Should not be here"<<endl;
  exit(1);
}
