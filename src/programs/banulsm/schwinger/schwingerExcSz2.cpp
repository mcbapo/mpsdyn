
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"
//#include "HermitianTensorMultiplier.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "SchwingerHamiltonianSz.h"
#include "SpinMPO.h"
#include <cmath>

using namespace std;
using namespace shrt;

#ifdef MACOSX
#define LOCALDIR "/Users/banulsm/Desktop/SVN/projects/mpsdyn/src/bin/"
#else
#define LOCALDIR "/home/banulsm/ALaPorra/bin/"
#endif
#define MAXLEN 120


/** Auxiliary function which constructs a blocked version of a MPO */
void blockMPO(const MPO& orig,MPO& final,int nrSites,int Llow);

/** Construct a block version given a list of Operator pointers. Returns a new Operator. */
Operator* blockOps(Operator** list,int nrsites);

/** SchwingerExcSz tries to find the scalar mass gap by first finding
    the MPS approximation to the ground state, then 
    finding another one for the first excited state, but in every
    step, the Hamiltonian is projected onto the Sz=0 subspace.

    Receives as argument a Porperties file, such as
    config/excitProp.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.

    To replace some of the values in the file by a command line
    option, the additional arguments have to be of the form
    -propName=propValue
    (except for the vector properties)
*/

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  const char* infile=argv[++cntr];
  string directory="";
  // if(argc>2){
  //   const char* indir=argv[++cntr];
  //   directory=string(indir)+"/";
  // }
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  int L = props.getIntProperty("L");
  double mg=props.getDoubleProperty("mg");
  double alpha=props.getDoubleProperty("alpha");
  double x=props.getDoubleProperty("x");
  int D0 = props.getIntProperty("D0");
  int nLev = props.getIntProperty("nLev");
  directory=props.getProperty("outputdir");
  const string outfile=directory+"/"+props.getProperty("outfile");
  int app=props.getIntProperty("append");
  double tol=props.getDoubleProperty("tol");
  int unif=props.getIntProperty("unif");
  int Dall=props.getIntProperty("Dcommon");
  double zpen=props.getDoubleProperty("penalty");
  double offset=props.getDoubleProperty("offset");
  if(offset==-1) offset=0;

  vector<int> Ds=props.getVectorIntProperty("Ds");
  while(Ds.size()<nLev)
    Ds.push_back(Ds.back());
  //cout<<"Read vector property "<<Ds<<endl;
  if(unif==1){ // ignore the different bond dimensions
    for(int k=0;k<Ds.size();k++)
      Ds[k]=Dall;
  }
  // Parameters I do not use here: offset and chemical potential
  double nu=0;
  // I might want to switch off the gauge terms
  double gweight=props.getDoubleProperty("gweight");
  if(gweight==-1) gweight=1; // default: do not alter

  double mu=2*mg*sqrt(x);
  double L0=0.; // l0 parameter, which is irrelevant, as only appears
		// as alpha+l0

  //  double t0=.01; // todo: differently -> reference time

  cout<<"Initialized arguments: L="<<L
      <<", mu="<<mu
      <<", x="<<x
      <<", alpha="<<alpha
      <<", outfile="<<outfile
      <<", app="<<app
      <<", D0="<<D0<<endl;

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
    out=new ofstream(outfile.data());
  }
  else{
    out=new ofstream(outfile.data(),ios::app);
    *out<<"%%%%%%%%%%%%%%%%%%%"<<endl;
  }
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output (append="<<app<<")"<<endl;
    exit(1);
  }

  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  int d=2;
  MPS gs(L,D0,d);
  gs.setRandomState(); // the intial state, random
  gs.gaugeCond('R',1);
  gs.gaugeCond('L',1);
   
  // First: put H and find the GS
  int Dop=5; // bond dimension of Schwinger H
  double lambda=0.;
  SchwingerHamiltonianSz hSch(L,mu,x,L0,alpha,zpen,offset,nu,gweight);
  //  SchwingerHamiltonian hSch(L,mu,x,L0,alpha,offset,nu,gweight);
  const MPO& hamil=hSch.getHMPO();
  // To check!!
  //  hamil.exportForMatlab("schwHamil.m");
  //  exit(1);
  contractor.findGroundState(hamil,D0,&lambda,gs);
  cout<<"Ground state found with eigenvalue "<<lambda<<endl;
  *out<<"% L="<<L<<endl;
  *out<<"% mg="<<mg<<endl;
  *out<<"% x="<<x<<endl;
  *out<<"% alpha="<<alpha<<endl;
  *out<<"% D0="<<D0<<endl;
  *out<<setprecision(10);
  *out<<"% E0="<<lambda<<endl;
  
  { // try to compute the expectation value of the condensate, too
    MPO cond(L);
    hSch.constructCondensateMPO(cond);
    complex_t fC=contractor.contract(gs,cond,gs);
    *out<<"% <condensate>="<<real(fC)<<endl;
  }

  // Observables I will compute
  // Compute the expectation value of Sz, as it commutes with H
  MPO Szmpo(L);
  SpinMPO::getSzMPO(L,d,Szmpo);
  MPO Sz2mpo(L);
  SpinMPO::getSz2MPO(L,d,Sz2mpo);
  // The momentum operator
  MPO Pmpo(L);
  hSch.constructMomentumMPO(Pmpo);
  const MPO* oprsP[2]={&Pmpo,&Pmpo};
  MPO P2mpo(L);
  MPO::join(2,oprsP,P2mpo);
  // The momentum squared operator
  MPO PSqmpo(L);
  hSch.constructMomentumSquaredMPO(PSqmpo);

  MPO ZZempo(L),ZZompo(L);
  hSch.constructZZMPO(ZZempo,true);
  hSch.constructZZMPO(ZZompo,false);

  MPO SPmpo(L);
  hSch.constructShiftPairMPO(SPmpo);

  // the shift by one and rotate sigmax
  MPO SRmpo(L);
  hSch.constructShiftRotateMPO(SRmpo);
  SRmpo.exportForMatlab("SRmpo.m");
  const MPO* oprsSR[2]={&SRmpo,&SRmpo};
  MPO SR2mpo(L);
  MPO::join(2,oprsSR,SR2mpo);  

  // the shift by one and rotate sigmax, but to left
  MPO SRmpoL(L);
  hSch.constructShiftRotateMPO(SRmpoL,false);
  SRmpo.exportForMatlab("SRmpoL.m");

  // the shift by one and rotate sigmax, to the right but cyclic
  MPO TRmpo(L);
  hSch.constructCyclicTranslationRotateMPO(TRmpo,true);
  MPO TRmpoL(L); // idem to left
  hSch.constructCyclicTranslationRotateMPO(TRmpoL,false);
  TRmpo.exportForMatlab("TRmpo.m");



  // // A special (product) state: the vacuum for x=0, to check the coefficients
  // MPS sceVac(L,1,d);sceVac.setProductState(p_zero);
  // MPS sceVacC(L,1,d);sceVacC.setProductState(p_one); // the complementary one
  // for(int k=0;k<L;k++){
  //   if(k%2!=0)
  //     sceVac.setA(k,sceVacC.getA(k).getA());
  //   else
  //     sceVacC.setA(k,sceVac.getA(k).getA());
  // }


  vector<MPS*> levels;
  levels.push_back(&gs); 
  // For each excited state, repeat the following procedure
  for(int l=0;l<=nLev;l++){
    if(l>0){ // compute one new excitation)
      int D1=Ds[l-1];
      cout<<"Initializing computation of excitation nr "<<l<<" with D="<<D1<<endl;
      // Now find the first excited state
      MPS* exc=new MPS(L,D1,d);
      exc->setRandomState(); // the initial state, random
      double lambda1=0.;
      contractor.findNextExcitedState(hamil,D1,levels,&lambda1,*exc);
      levels.push_back(exc);
    }
    // Compute the expectation values
    const MPS& exc=*levels[l];
    //MPS reflected(exc);
    //hSch.applyReflectionRotation(exc,reflected);
    complex_t Hval=contractor.contract(exc,hamil,exc);
    complex_t H2val=contractor.contract2(hamil,exc);
    complex_t normN=contractor.contract(exc,exc);
    complex_t Szval=contractor.contract(exc,Szmpo,exc);
    complex_t Sz2val=contractor.contract(exc,Sz2mpo,exc);
    complex_t Pval=contractor.contract(exc,Pmpo,exc);
    complex_t P2val=contractor.contract(exc,P2mpo,exc);
    complex_t PSqval=contractor.contract(exc,PSqmpo,exc);
    complex_t ZZeval=contractor.contract(exc,ZZempo,exc);
    complex_t ZZoval=contractor.contract(exc,ZZompo,exc);
    complex_t SPval=contractor.contract(exc,SPmpo,exc);
    complex_t SRval=contractor.contract(exc,SRmpo,exc);
    complex_t SRvalL=contractor.contract(exc,SRmpoL,exc);
    complex_t SR2val=contractor.contract(exc,SR2mpo,exc);
    complex_t TRval=contractor.contract(exc,TRmpo,exc);
    complex_t TRvalL=contractor.contract(exc,TRmpoL,exc);
    //complex_t Parity=contractor.contract(exc,reflected);
    //    complex_t SCEcoeff=contractor.contract(exc,sceVac);
    // complex_t SCEcoeff2=contractor.contract(exc,sceVacC);

    *out<<l<<"\t"<<exc.getBond()<<"\t";
    *out<<real(normN)<<"\t";
    *out<<real(H2val)<<"\t"; //<<imag(H2)<<"\t";
    *out<<real(Hval-zpen*Sz2val)<<"\t"; //<<imag(Hval)<<"\t";
    //*out<<real(H2/normN-Hval*Hval/(normN*normN))<<"\t";
    *out<<real(Szval)<<"\t";
    //*out<<real(Pval)<<"\t";
    *out<<Pval<<"\t";
    *out<<real(P2val)<<"\t";
    *out<<real(PSqval)<<"\t";
    *out<<real(ZZeval)<<"\t";
    *out<<real(ZZoval)<<"\t";
    *out<<2*real(SPval)<<"\t";
    *out<<real(SRval)<<"\t"<<imag(SRval)<<"\t";
    *out<<real(SRvalL)<<"\t"<<imag(SRvalL)<<"\t";
    *out<<real(SR2val)<<"\t"<<imag(SR2val)<<"\t";
    *out<<real(TRval)<<"\t"<<imag(TRval)<<"\t";
    *out<<real(TRvalL)<<"\t"<<imag(TRvalL)<<"\t";
    //*out<<real(Parity)<<"\t"<<imag(Parity)<<"\t";
    //*out<<SCEcoeff<<"\t";
    //*out<<SCEcoeff2<<"\t";
    for(int s=0;s<=nLev;s++)
      if(s<l)
	*out<<contractor.contract(exc,*levels[s])<<"\t";
      else if(s==l)
	*out<<1<<"\t";
      else
	*out<<0<<"\t";
    *out<<endl;
  }

  out->close();
  delete out;

   // for(int z=1;z<L;z++){
   //   delete Skmpos[z];
   // }
   // delete []Skmpos;
}


void blockMPO(const MPO& orig,MPO& final,int midSites,int Llow){
  int L=orig.getLength();
  final.initLength(L-midSites+1);
   for(int k=0;k<L-midSites+1;k++){
     if(k<Llow){
       //final.setOp(k,&orig.getOp(k),false);
       final.setOp(k,new Operator(orig.getOp(k).getFullData()),true);
     }
     else if(k>Llow){
       // final.setOp(k,&orig.getOp(k+midSites-1),false);
       final.setOp(k,new Operator(orig.getOp(k+midSites-1).getFullData()),true);
     }
     else{ // the middle site!
       mwArray effSite=orig.getOpData(Llow); // du dl dd dr
       int du=effSite.getDimension(0);
       int dd=effSite.getDimension(2);
       int Dl=effSite.getDimension(1);int Dr=effSite.getDimension(3);
       for(int p=1;p<midSites;p++){
	 effSite.reshape(Indices(du*Dl*dd,Dr));
	 mwArray auxOp=orig.getOpData(Llow+p); // du Dr dd Drp
	 Indices newdims=auxOp.getDimensions();
	 // cout<<"Adding site "<<p+Llow<<" to effective operator, dims now are"
	 //     <<effSite.getDimensions()<<" and the new piece has "<<newdims<<endl;
	 auxOp.permute(Indices(2,1,3,4));
	 auxOp.reshape(Indices(newdims[1],newdims[0]*newdims[2]*newdims[3]));
	 effSite.multiplyRight(auxOp);
	 effSite.reshape(Indices(du,Dl,dd,newdims[0],newdims[2],newdims[3]));
	 effSite.permute(Indices(1,4,2,3,5,6));
	 du=du*newdims[0];dd=dd*newdims[2];
	 effSite.reshape(Indices(du,Dl,dd,Dr));
       }
       final.setOp(k,new Operator(effSite),true);
     }
   }
}

 Operator* blockOps(Operator** list,int nrsites){
   mwArray aux=list[0]->getFullData();
   int dimU=aux.getDimension(0);
   int Dl=aux.getDimension(1);
   int dimD=aux.getDimension(2);
   int Dr=aux.getDimension(3);
   for(int s=1;s<nrsites;s++){
     aux.reshape(Indices(dimU*Dl*dimD,Dr));
     mwArray aux2=list[s]->getFullData();
     Indices dim2=aux2.getDimensions();
     aux2.permute(Indices(2,1,3,4));
     aux2.reshape(Indices(dim2[1],dim2[0]*dim2[2]*dim2[3]));
     aux.multiplyRight(aux2);
     aux.reshape(Indices(dimU,Dl,dimD,dim2[0],dim2[2],dim2[3]));
     aux.permute(Indices(1,4,2,3,5,6));
     dimU*=dim2[0];
     dimD*=dim2[2];
     Dr=dim2[3];
     aux.reshape(Indices(dimU,Dl,dimD,Dr));
   }
   return new Operator(aux);
}
