
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

/** thirringP runs the findGroundState routine with the MPO for the Thirring
    Hamiltonian \ref <ThirringHamiltonian> with possibly a penalty term for a
    targetted total \f$S_z=S_{\mathrm{tot}}\f$.  At the end, the total
    \f$S_z\f$ is computed.  In the version thirringPenaltyEVs, expectation
    values of SxSx, SySy, SzSz and (-1)^n Sz are stored too, to be able to
    compute expectation values of Hamiltonians with other parameters.

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


int findTmpFile(int D,const string& mpsdir,const string& baseName,string& mpsfile,int level=0);
const string mpsFileName(const string& mpsdir,const string& baseName,int D,int level=0);
const string baseFileName(const Properties& props);
void check(const MPS& MPS,double gt2,double g2,double ma,double mu,const MPO& hamil,double lambda,double Stot); // debugging only
complex_t valZ(const MPS& gs,int pos);

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
  double g2_=props.getDoubleProperty("g2"); // the non-alternating one (Delta)
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
  //  string mpsfile=props.getProperty("mpsfile");
  //  string mpstostretch=props.getProperty("stretchmps");
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
  MPS gs(L,D,d);
    // Check if initial state was given
  string mpsfile=baseFileName(props); // base of file name containing these parameters

  string tmpMPSfile;
  if(findTmpFile(D,mpsdir,mpsfile,tmpMPSfile)){// found some tmp or smaller D file
    cout<<"Recovered initial guess from file "<<tmpMPSfile<<endl;
    if(gs.getLength()!=L){
      cout<<"Error! File "<<tmpMPSfile<<" does not contain valid MPS "<<endl;
      exit(1);
    }
    gs.importMPS(tmpMPSfile.data());
  }
  else{
    //    gs.setProductState(p_one);
    srandom(time(NULL));
    gs.setRandomState(); // the intial state, random
    // but make it random TI
    // mwArray randA=gs.getA(L/2).getA();
    // for(int k=0;k<L;k++){
    //   mwArray aux(randA);
    //   aux.resize(gs.getA(k).getDimensions());
    //   gs.replaceSite(k,aux);
    // }
    cout<<"Initialized random state, norm "<<contractor.contract(gs,gs)<<endl;
    //    gs.exportMPStext("randomMPS.txt");
  }
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
  cout<<setprecision(10);
  
  // First: construct H
  double E0=0.;
  
  MPO Szmpo(L);
  SpinMPO::getSzMPO(L,d,Szmpo);

 
  ThirringHamiltonian hamH(L,ma,g2_,lambda_,Stot_,d,gt2_,mu_);
  //cout<<"Created the Hamiltonian"<<endl;
  const MPO& hamil=hamH.getHMPO();
  if(L<12)
    hamil.exportForMatlab("thirringHmpo.m");
  //  cout<<"Initial E value, with initial random state"<<endl;
  cout<<contractor.contract(gs,hamil,gs)<<endl;
  //cout<<"Initial Sz value, with initial random state"<<endl;
  cout<<contractor.contract(gs,Szmpo,gs)<<endl;

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
  complex_t Sztot=contractor.contract(gs,Szmpo,gs);
  cout<<"Ground state saved to "<<tmpMPSfile<<endl;
  cout<<"Energy contraction is:"<<contractor.contract(gs,hamil,gs)<<" and Sz contraction "<<Sztot<<endl;
  

  cout<<"Ground state found for S_tot="<<Stot_<<" with energy "<<E0
      <<" (Energy per particle="<<E0/L<<", plus offset from David->"
      <<E0/L+gt2_*.25+g2_*(L-1.)/(4.*L)<<")"
      <<" Total Sz="<<real(Sztot)<<endl;

  complex_t Sz2=contractor.contract2(Szmpo,gs);
  complex_t H2=contractor.contract2(hamil,gs);
  cout<<"DeltaH="<<real(H2)-E0*E0<<endl;
  cout<<"Delta S_z="<<Sz2-Sztot*Sztot<<endl;

  cout<<"Computing pieces of H "<<endl;
  MPO Hterm(L);
  SpinMPO::getSxSxMPO(L,d,Hterm);
  complex_t SxSxtot=contractor.contract(gs,Hterm,gs);
  cout<<"Obtained <SxSx>="<<SxSxtot<<endl;
  SpinMPO::getSySyMPO(L,d,Hterm);
  complex_t SySytot=contractor.contract(gs,Hterm,gs);
  cout<<"Obtained <SySy>="<<SySytot<<endl;
  SpinMPO::getSzSzMPO(L,d,Hterm);
  complex_t SzSztot=contractor.contract(gs,Hterm,gs);
  cout<<"Obtained <SzSz>="<<SzSztot<<endl;
  SpinMPO::getSzMPO(L,d,Hterm,true);
  complex_t SzStaggered=contractor.contract(gs,Hterm,gs);
  cout<<"Obtained <(-1)^n Sz>="<<SzStaggered<<endl;

  //check(gs,gt2_,g2_,ma,mu_,hamil,lambda_,Stot_); // just a debugging piece to check energy terms

  // I will also save <Sz> in first and last sites
  

  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    *out<<"% N\t ma\t gt2\t g2\t Starget\t lambda\t D\t level\t Energy\t  H \t S_z\t Delta Sz\t E/L (David)\t<SxSx>\t<SySy>\t<SzSz>\t(-1)^n Sz>\t<Sz(0)>\t<Sz(L-1)>"<<endl;
  }
  else
    out=new ofstream(outfname.data(),ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(15);
  *out<<L<<"\t"<<ma<<"\t"<<gt2_<<"\t"<<g2_<<"\t"<<Stot_<<"\t"<<lambda_<<"\t"
      <<D<<"\t"<<0<<"\t"<<E0<<"\t"<<real(H2)-E0*E0<<"\t"<<real(Sztot)<<"\t"<<real(Sz2-Sztot*Sztot)
      <<"\t"<<E0/L+gt2_*.25+g2_*(L-1.)/(4.*L)
      <<"\t"<<real(SxSxtot)<<"\t"<<real(SySytot)<<"\t"<<real(SzSztot)<<"\t"<<real(SzStaggered)
      <<"\t"<<real(valZ(gs,0))<<"\t"<<real(valZ(gs,L-1));
  *out<<endl;

  // cout<<"Checking!!!!"<<endl;
  // for(int k=0;k<L;k++)
  //   cout<<k<<"\t"<<real(valZ(gs,k))<<endl;
  // exit(1);
  
  ofstream* outS;
  if(saveS){
    if(!app){
      outS=new ofstream(outfnameS.data());
      *outS<<"% N\t ma\t gt2\t g2\t D\t level\t Sz\t pos\t S_z\t S(0-pos:rest)"<<endl;
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
    for(int pos=0;pos<L;pos++){
      aux.applyLocalOperator(pos,sigz);
      complex_t valZ=.5*contractor.contract(aux,gs);
      // undo the change
      aux.applyLocalOperator(pos,sigz);
      // And the entropy
      double Spos=contractor.getEntropy(gs,pos+1);
      *outS<<L<<"\t"<<ma<<"\t"<<gt2_<<"\t"<<g2_<<"\t"<<D<<"\t"
	   <<0<<"\t"<<real(Sztot)<<"\t"
	   <<pos<<"\t"<<real(valZ)<<"\t"<<Spos<<endl;      
    }
    
  }
  
  vector<MPS*> levels;
  levels.push_back(&gs); 
  int k=1;
  while(k<=kmax){
    cout<<"Computing now next excited state"<<endl;
    double E1=0.;
    MPS* exc=new MPS(L,D,d);
    if(findTmpFile(D,mpsdir,mpsfile,tmpMPSfile,k)){// found some tmp or smaller D file
      cout<<"Recovered initial guess for excitation from file "<<tmpMPSfile<<endl;
      if(exc->getLength()!=L){
	cout<<"Error! File "<<tmpMPSfile<<" does not contain valid MPS "<<endl;
	exit(1);
      }
      exc->importMPS(tmpMPSfile.data());
    }
    else{
      exc->setRandomState(); // the initial state, random
    }
    contractor.findNextExcitedState(hamil,D,levels,&E1,*exc);

    tmpMPSfile=mpsFileName(mpsdir,mpsfile,D,k);
    exc->exportMPS(tmpMPSfile.data());
    cout<<"Level "<<k<<" state saved to "<<tmpMPSfile<<endl;

    Sztot=contractor.contract(*exc,Szmpo,*exc);
    levels.push_back(exc); 
    Sz2=contractor.contract2(Szmpo,*exc);
    H2=contractor.contract2(hamil,*exc);
    complex_t overl=contractor.contract(*exc,gs);
    
    cout<<k<<"-th excited state state found with eigenvalue "<<E1
	<<" (Gap="<<E1-E0<<", <gs|lev("<<k<<")>="<<overl<<")"
	<<" Total Sz="<<real(Sztot)<<endl;
    *out<<L<<"\t"<<"\t"<<ma<<"\t"<<gt2_<<"\t"<<g2_<<"\t"<<Stot_<<"\t"<<lambda_<<"\t"
	<<D<<"\t"<<k<<"\t"<<E1<<"\t"<<real(H2)-E1*E1<<"\t"<<real(Sztot)
	<<"\t"<<real(Sz2)-real(Sztot)*real(Sztot)
	<<"\t"<<E1/L+gt2_*.25+g2_*(L-1.)/(4.*L);
    *out<<endl;
    if(saveS){
      MPS aux(*exc);aux.gaugeCond('R',1);
      complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
      mwArray sigz(Indices(d,d),dataZ); // sigmaZ operator on one site
      for(int pos=0;pos<L;pos++){
	aux.applyLocalOperator(pos,sigz,0);
	complex_t valZ=.5*contractor.contract(aux,*exc);
	// undo the change
	aux.applyLocalOperator(pos,sigz,0);
	// And the entropy
	double Spos=contractor.getEntropy(*exc,pos+1);
	*outS<<L<<"\t"<<ma<<"\t"<<gt2_<<"\t"<<D<<"\t"
	     <<k<<"\t"<<real(Sztot)<<"\t"
	     <<pos<<"\t"<<real(valZ)<<"\t"<<Spos<<endl;      
      }
    }
    k++;
  }
  
  out->close();
  delete out;
  if(saveS){
    outS->close();
    delete outS;
  }
}

const string baseFileName(const Properties& props){
  int L=props.getIntProperty("L");
  double ma=props.getDoubleProperty("ma");  
  double gt2_=props.getDoubleProperty("gtilde2");
  double g2_=props.getDoubleProperty("g2"); // the non-alternating one
  double lambda_=props.getDoubleProperty("lambda");
  if(lambda_<0) lambda_=0.;
  double mu_=props.getDoubleProperty("mu");
  double tol=props.getDoubleProperty("tol");
  if(tol<0) tol=1E-6;
  int Stot_=props.getIntProperty("Starget");
  stringstream s;
  s<<"MPSthirring_L"<<L<<"_ma"<<ma<<"_gt2"<<gt2_<<"_g2"<<g2_<<"_mu"<<mu_<<"_Stot"<<Stot_;
  return s.str();
}

const string mpsFileName(const string& mpsdir,const string& baseName,int D,int level){
  stringstream s;
  s<<mpsdir<<"/"<<baseName;
  if(level>0) s<<"_"<<level;
  if(D>0) s<<"_"<<D;
  return s.str();
}

int findTmpFile(int D,const string& mpsdir,const string& baseName,string& mpsfile,int level){
  mpsfile=mpsFileName(mpsdir,baseName,D,level)+"_tmp";
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
      mpsfile=mpsFileName(mpsdir,baseName,cnt,level);
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

complex_t valZ(const MPS& gs,int pos){
  int d=2;
  complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  Contractor& contractor=Contractor::theContractor();
  MPS aux(gs);
  aux.applyLocalOperator(pos,sigZ,0);
  return contractor.contract(aux,gs);  
}


void check(const MPS& gs,double gt2,double g2,double ma,double mu,const MPO& hamil,double lambda,double Stot){
  int d=2;
  mwArray sig0=identityMatrix(d);
  complex_t dataX[]={ZERO_c,.5*ONE_c,.5*ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  complex_t dataY[]={ZERO_c,.5*I_c,-.5*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  complex_t dataZ[]={.5*ONE_c,ZERO_c,ZERO_c,-.5*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);


  vector<complex_t> valsXX;
  vector<complex_t> valsYY;
  vector<complex_t> valsZZ;
  vector<complex_t> valsZ;

  Contractor& contractor=Contractor::theContractor();
  int L=gs.getLength();
  MPS aux(gs);
  for(int i=0;i<L-1;i++){

    aux.applyLocalOperator(i,sigX,0);
    aux.applyLocalOperator(i+1,sigX,0);
    valsXX.push_back(contractor.contract(aux,gs));
    
    aux.applyLocalOperator(i,4*sigY*sigX,0);
    aux.applyLocalOperator(i+1,4*sigY*sigX,0);
    valsYY.push_back(contractor.contract(aux,gs));
    
    aux.applyLocalOperator(i,4*sigZ*sigY,0);
    aux.applyLocalOperator(i+1,4*sigZ*sigY,0);
    valsZZ.push_back(contractor.contract(aux,gs));
    
    aux.applyLocalOperator(i+1,4*sigZ,0);
    valsZ.push_back(contractor.contract(aux,gs));

    aux.applyLocalOperator(i,4*sigZ,0);
    
  }
  aux.applyLocalOperator(L-1,sigZ,0);
  valsZ.push_back(contractor.contract(aux,gs));

  
  complex_t totXX=ZERO_c;
  complex_t totYY=ZERO_c;
  complex_t totZZ=ZERO_c;
  complex_t totZZst=ZERO_c;
  complex_t totZ=ZERO_c;
  complex_t totZst=ZERO_c;
  for(int i=0;i<L;i++){
    int signPos=(i%2==0)?1:-1;
    if(i<L-1){
      totXX+=valsXX[i];
      totYY+=valsYY[i];
      totZZ+=valsZZ[i];
      totZZst+=signPos*valsZZ[i];
    }
    totZ+=valsZ[i];
    totZst+=signPos*valsZ[i];
  }

  cout<<"After checking: \n\ttotXX="<<totXX<<"\n\ttotYY="<<totYY<<"\n\ttotZZ="<<totZZ<<"\n\ttotZ="<<totZ<<"\n\ttotZst="<<totZst<<endl;

  cout<<"Energy should be "<<-1.*(totXX+totYY)+g2*totZZ+gt2*(totZZ+totZZst)+ma*totZst+(gt2+mu+g2)*totZ-g2*.5*(valsZ[0]+valsZ[L-1])<<endl;

  complex_t Ec=contractor.contract(gs,hamil,gs);
  complex_t normV=contractor.contract(gs,gs);
  cout<<"And contraction is E="<<Ec<<" (norm "<<normV<<")"<<endl;
  cout<<"Substracting penalty:"<<Ec-lambda*(totZ-Stot*ONE_c)*(totZ-Stot*ONE_c)<<endl;

  // Experimenting with other hamil
  ThirringHamiltonian hamH1(L,ma*0,g2*0,lambda*0,Stot,d,gt2*0,mu*0); //only hopping
  const MPO& hamil1=hamH1.getHMPO();
  complex_t valHop=contractor.contract(gs,hamil1,gs);
  cout<<"Expectation value of only hopping term: "<<valHop<<endl;
  cout<<"\t in my calculation "<<-1.*(totXX+totYY)<<endl;
  
  ThirringHamiltonian hamH2(L,ma,g2*0,lambda*0,Stot,d,gt2*0,mu*0); //only hopping
  const MPO& hamil2=hamH2.getHMPO();
  complex_t valHopMass=contractor.contract(gs,hamil2,gs);
  cout<<"Expectation value of only hopping and mass terms: "<<valHopMass<<", from which mass term is "<<valHopMass-valHop
      <<" and <(-1)^nSz>="<<(valHopMass-valHop)/ma<<endl;
  cout<<"\t in my calculation "<<ma*totZst<<endl;

  ThirringHamiltonian hamH3(L,ma*0,g2,lambda*0,Stot,d,gt2*0,mu*0); //only hopping
  const MPO& hamil3=hamH3.getHMPO();
  complex_t valHopG2=contractor.contract(gs,hamil3,gs);
  cout<<"Expectation value of only hopping and Delta terms: "<<valHopG2<<", from which contrib of G2 is "<<valHopG2-valHop
      <<" (without G2 "<<(valHopG2-valHop)/g2<<")"<<endl;
  cout<<"\t in my calculation "<<g2*totZZ-g2*.5*(valsZ[0]+valsZ[L-1])+g2*totZ<<endl;

  ThirringHamiltonian hamH4(L,ma*0,g2*0,lambda*0,Stot,d,gt2,mu*0); //only hopping
  const MPO& hamil4=hamH4.getHMPO();
  complex_t valHopGt2=contractor.contract(gs,hamil4,gs);
  cout<<"Expectation value of only hopping and gtilde terms: "<<valHopGt2<<", from which contrib of Gt2 is "<<valHopGt2-valHop
      <<" (without Gt2 "<<(valHopG2-valHop)/gt2<<")"<<endl;
  cout<<"\t in my calculation "<<gt2*(totZZ+totZZst)+gt2*totZ<<endl;

  cout<<"For these pieces, recompose totatl energy: "<<valHop+valHopMass+valHopG2+valHopGt2-3*valHop<<endl;
  
  
  //cout<<valsZ<<endl;

  // Now as theirs
  vector<complex_t> valsProjZ;

  vector<complex_t> valsProjZZ;
  mwArray projZ=sigZ+.5*sig0;
  for(int i=0;i<L;i++){    
    MPS aux(gs);
    aux.applyLocalOperator(i,projZ,0);
    valsProjZ.push_back(contractor.contract(aux,gs));

    if(i<L-1){
    aux.applyLocalOperator(i+1,projZ,0);
    valsProjZZ.push_back(contractor.contract(aux,gs));
    }
  }
  
  complex_t totProjZZ=ZERO_c;
  complex_t totProjZst=ZERO_c;
  for(int i=0;i<L;i++){
    int signPos=(i%2==0)?1:-1;
    if(i<L-1){
      totProjZZ+=valsProjZZ[i];
    }
    totProjZst+=signPos*valsProjZ[i];
  }
  
  cout<<"And with David's form: "<<-1.*(totXX+totYY)+g2*totProjZZ+g2*totProjZst<<endl;
  
}
