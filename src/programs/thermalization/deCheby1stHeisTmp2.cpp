
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "HeisenbergHamiltonian.h"

using namespace shrt;

void sumMPO(MPO& result,const MPO& mpo1,const MPO& mpo2);
void computeTwoBodyCorrCenter(MPS& result,const MPS& rho,const mwArray& op1,
                              const mwArray& op2,int l);


double getJacksonKPMcoeff(int k,int M);
/** Construct a file name to save (tmp) results */
const string mpsFileName(const string& mpsdir,const string& baseName,const string& label,int cnt);

/** Save the tmp files */
void tmpSave(MPS& mps,double norm_MPS,MPS& chebyU_Nm1,double norm_Nm1,MPS& chebyU_N,
             double norm_N,const string& baseDir,const string& baseName,int cnt,vector<string>& tmpFiles);
//void getSchwingerInitState(MPS &mps,int N,InitState state);
bool constructInitialRho(int L,int d,int initSt,MPS& rho0);


/** Find tmp files and restore the status */
int tmpRestore(MPS& mps,double& norm_MPS,MPS& chebyU_Nm1,double& norm_Nm1,MPS& chebyU_N,
               double& norm_N,const string& baseDir,const string& baseName,int M,vector<string>& tmpFiles);
/** minVarChebyDelta constructs a Chebyshev poly approx to the delta
    function and uses it to filter out large energies, and thus
    decrease the variance.

    In this version, the expansion approximates directly the product
    of the Chebysev polynomial times the initial state, since the recursion is
    the same.

 */

int d=2;

int main(int argc,const char* argv[]){
  // Read input arguments
  int cntr=0;
  srandom(time(NULL));
  const char* infile=argv[++cntr];
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }

  int L=props.getIntProperty("L");
  double Jx_=props.getDoubleProperty("Jx");
  double Jy_=props.getDoubleProperty("Jy");
  double Jz_=props.getDoubleProperty("Jz");
  double h_=props.getDoubleProperty("h");
  double delta=props.getDoubleProperty("delta");
  if(delta<0) delta=0.02; 
  int D = props.getIntProperty("D");  
  string outfname=props.getProperty("output");
  string mpsdir=props.getProperty("mpsdir");
  string mpsfile=props.getProperty("mpsfile");
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  int M=props.getIntProperty("M"); // LARGEST ORDER
  double Emin=props.getDoubleProperty("Emin");
  double Emax=props.getDoubleProperty("Emax");
  string efname=props.getProperty("energiesfile");
  int initSt=props.getIntProperty("initSt"); // initial state to use:
  //                               +1=X+; -1=X-; (2,3 for Y, Z)
  //                               +4,+5,+6 are the staggered(+-) for X, Y and Z and -4,-5,-6 the -+
  //                               +7 is a (complex) random state and +8 a real random (in which cases, the random seed may be given)
  // These three parameters are obsolete
  bool initStag=(abs(initSt)>=4&&abs(initSt)<=6); //(props.getIntProperty("staggered")>0);
  bool initRand=(initSt>=7); //(props.getIntProperty("random")>0);
  bool initY=(initSt==+2); //(props.getIntProperty("initY")>0);
  bool realInit=(abs(initSt)!=2&&abs(initSt)!=5&&initSt!=7);
  int randSeed=props.getIntProperty("seed");


  int rate=props.getIntProperty("rate");
  if(rate<0) rate=1;
  int saveFreq=props.getIntProperty("saveFreq");
  if(saveFreq<0) saveFreq=M;
  //  saveFreq=lcm(rate,saveFreq);
  int lastSavedCnt=-1;

  if(saveFreq>0&&mpsdir.size()==0){
    cout<<"ERROR: Cannot save tmp results because no dir is given as argument!"<<endl;
    exit(1);
  }
 


  cout<<"Initialized arguments: L="<<L
      <<", Jx="<<Jx_
      <<", Jy="<<Jy_
      <<", Jz="<<Jz_
      <<", h="<<h_
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;

  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfname<<
	" for output"<<endl;
      exit(1);
    }
    *out<<"% cnt\t D\t real(tr)\t im(tr)\t re(E)\t im(E)\t re(Sx)\t re(Sy)\t re(Sz)\tcomm";
    *out<<endl;
    out->close();
    delete out;
  }
  cout<<setprecision(10);

  
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  
  // The original Hamiltonian (I'll also need the rescaled one for the Chebyshev expansion)
  HeisenbergHamiltonian hamH0(L,4*Jx_,4*Jy_,4*Jz_,2*h_);
  const MPO& hamil0=hamH0.getHMPO();
    // The original Hamiltonian to get the energy
  MPS hmps0(L,1,d*d);
  MPSfromMPO(hamil0,hmps0);


  double scale,offset;
  if(Emin==-1&&Emax==-1){// not provided
    cout<<"First need to estimate the energy band to rescale "
	<<"and shift the Hamiltonian, and accordingly the commutator"<<endl;
    // First rescale the Hamiltonian:
    //hamH.registerTimeDependence();
    const MPO& hamil=hamH0.getHMPO();
    // Get ground and highest state energy, for rescaling
    Emin=0;
    MPS gs(L,D,d);
    gs.setRandomState();gs.gaugeCond('R');gs.gaugeCond('L',1);
      contractor.findGroundState(hamil,D,&Emin,gs);
      cout<<"Found GS of H with E="<<Emin<<endl;
      // Same for largest eigenvalue:
      Emax=0;
      MPS exc(L,D,d);
      exc.setRandomState();exc.gaugeCond('R');exc.gaugeCond('L',1);
      {
	MPO hamilMinus(L);
	hamilMinus.setOp(0,new Operator(-1.*hamil.getOp(0).getFullData()),true);
	for(int k=1;k<L;k++){
	  hamilMinus.setOp(k,&hamil.getOp(k),false);
	}
	contractor.findGroundState(hamilMinus,D,&Emax,exc);
	Emax=-Emax;
      }
      cout<<"Found max Exc of H with E="<<Emax<<endl;
      if(efname.size()!=0){ // save for other cases
	ofstream* outE=new ofstream(efname.data());
	*outE<<setprecision(15);
	*outE<<Emin<<endl;
	*outE<<Emax<<endl;
	outE->close();delete outE;
      }
    }
  else{
    cout<<"Using for extreme energies the provided arguments Emin="<<Emin
	<<" Emax="<<Emax<<endl;
  }
  // Now the scale has to be such that E fits in [-1+delta,1+delta]
  scale=(1.-delta)/(Emax-Emin);
  offset=0.;
  
  // Now the rescaled H, such that the spectrum of [H,.] is in [-1+delta,1+delta]
  HeisenbergHamiltonian hamH(L,Jx_*scale,Jy_*scale,Jz_*scale,h_*scale);
  //hamH.registerTimeDependence();
  cout<<"Created the rescaled Hamiltonian: scale="<<scale<<", offset="<<offset<<endl;
  const MPO& hamil_=hamH.getHMPO();
  int Dh=3; // bond dimension of H, just in case

  //Commutator superoperator
  MPO commMPO(L);
  hamH.getCommutatorMPO(commMPO); //commutator constructed
  cout<<"Commutator constructed"<<endl;

    // Prepare the observables to compute
  // The identity, to compute the trace
  MPS idMPS(L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++)
    idMPS.setA(k,id2);
  // And I will also need the sigmaz/sigmay/sigmax operator on central site, which I could compute ony by one, but I assemble it here
  MPS SxCenter(idMPS);
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  sigX.reshape(Indices(d*d,1,1));
  SxCenter.setA(L/2,sigX);sigX.reshape(Indices(d,d));

  MPS SyCenter(idMPS);
  complex_t dataY[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  sigY.reshape(Indices(d*d,1,1));
  SyCenter.setA(L/2,sigY);sigY.reshape(Indices(d,d));

  MPS SzCenter(idMPS);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  sigZ.reshape(Indices(d*d,1,1));
  SzCenter.setA(L/2,sigZ);sigZ.reshape(Indices(d,d));





    // Initial state rho0
  MPS rho0(L,D,d*d);
  bool isReal=constructInitialRho(L,d,initSt,rho0);
  rho0.gaugeCond('R',0);rho0.gaugeCond('L',0);
    double normMPS=log2(rho0.getNormFact()); // hopefully 0
   rho0.setNormFact(1.);
  complex_t E0=contractor.contract(rho0,hmps0);
  complex_t tr0=contractor.contract(rho0,idMPS);
  cout<<"Initial state has trace="<<tr0<<", <H>="<<E0/tr0<<endl;
    //  <<", |[H,rho0]|^2="<<pow(2,2*normMPS)*contractor.contract2(commMPO,rho0)/(scale*scale)<<endl;    
 
  
  // Now start the Chebyshev expansion computing U-1(0) and U0(Id) on the state
  MPS* chebyU_Nm1=new MPS(idMPS); // poly(-1)=0
  // Set the first tensor to 0=> ignore norm factor in this case
  chebyU_Nm1->gaugeCond('R',0);chebyU_Nm1->gaugeCond('L',1);
  chebyU_Nm1->setA(0,ZERO_c*chebyU_Nm1->getA(0).getA());
  double norm_Nm1=0.;
  
  MPS* chebyU_N=new MPS(rho0); // normalized already
  double norm_N=normMPS;

  double coeffK=getJacksonKPMcoeff(0,M);
  normMPS+=log2(coeffK)-log2(M_PIl);



  MPS SSz1(idMPS);
  MPS SSz2(idMPS);


  computeTwoBodyCorrCenter(SSz1,rho0,sigZ,sigZ,1);
  //  MPO auxSSz1(L), mpoSSz1(L);
  //MPOfromMPS(SSz1,auxSSz1);
  // extendMPO(auxSSz1,mpoSSz1,d);


  computeTwoBodyCorrCenter(SSz2,rho0,sigZ,sigZ,2);
  //MPO auxSSz2(L), mpoSSz2(L);
  // MPOfromMPS(SSz2,auxSSz2);
  //extendMPO(auxSSz2,mpoSSz2,d);

  
 
  complex_t trR=contractor.contract(rho0,idMPS);
  complex_t Eev=contractor.contract(rho0,hmps0);
  complex_t SxCev=contractor.contract(rho0,SxCenter);
  complex_t SyCev=contractor.contract(rho0,SyCenter);
  complex_t SzCev=contractor.contract(rho0,SzCenter);
   complex_t valComm=contractor.contract2(commMPO,rho0);
   double Entropy=contractor.getEntropy(rho0);

  complex_t valSSz1=contractor.contract(rho0,SSz1);
  complex_t valSSz2=contractor.contract(rho0,SSz2);  


  out=new ofstream(outfname.data(),ios::app);
  *out<<0<<"\t"<<D<<"\t"<<scale<<"\t"<<normMPS
      <<"\t"<<real(trR)<<"\t"<<imag(trR)
      <<"\t"<<real(Eev)<<"\t"<<imag(Eev)
    //  <<"\t"<<real(Sxev)<<"\t"<<imag(Sxev)
    //  <<"\t"<<real(Szev)<<"\t"<<imag(Szev)
      <<"\t"<<real(SxCev)/real(trR)
    //  <<"\t"<<real(SyCev)
      <<"\t"<<real(SzCev)/real(trR)
      <<"\t"<<real(valComm)/(scale*scale)
      <<"\t"<<Entropy
      <<"\t"<<real(valSSz1)
      <<"\t"<<real(valSSz2)
      <<endl;
  out->close();
  // If intermediate steps are saved, I need mps with norm_MPS, chebyT_N with norm_N and chebyT_Nm1 with norm_Nm1                                                                                           
   vector<string> tmpFiles;
  int cnt=0;
  if(mpsdir.size()>0)
    cnt=tmpRestore(rho0,normMPS,*chebyU_Nm1,norm_Nm1,*chebyU_N,norm_N,mpsdir,mpsfile,M,tmpFiles);
  if(cnt>0){
    cout<<"After reading tmp files for step "<<cnt<<", normMPS="<<normMPS<<", norm_N="<<norm_N<<", norm_Nm1="<<norm_Nm1<<endl;
    cnt++; // do not want to repeat the iteration for cnt                                                                                                                                                   
  }
  
  // iteration
  for(int k=cnt;k<=M;k++){
    MPS* chebyU_Np1;double norm_Np1(0.);
    if(k==0){
      chebyU_Np1=chebyU_N; 
      norm_Np1=norm_N;
    
    }
    else{
      chebyU_Np1=new MPS(*chebyU_N);
      if(k==1){
        // for Type I, T_1 is simply x*T0=> only multiply by H                                                                                                                                              
        contractor.optimize(commMPO,*chebyU_N,*chebyU_Np1,D);
        chebyU_Np1->gaugeCond('R',0);chebyU_Np1->gaugeCond('L',0);
        norm_Np1=log2(chebyU_Np1->getNormFact())+norm_N;chebyU_Np1->setNormFact(1.);
      }

    else{
      // chebyU_Np1=new MPS(*chebyU_N);
      // TODO: Optimize the sum of ops acting on vecs!
      MPS aux(*chebyU_N);
      contractor.optimize(commMPO,*chebyU_N,aux,D);
      aux.gaugeCond('R',0);aux.gaugeCond('L',0);
      double auxN=log2(aux.getNormFact());aux.setNormFact(1.);
      vector<const MPS*> kets;
      kets.push_back(&aux);kets.push_back(chebyU_Nm1);
      vector<complex_t> coefs;coefs.push_back(2.*ONE_c);
      coefs.push_back(-pow(2,norm_Nm1-norm_N-auxN)*ONE_c);
      contractor.optimizeSum(kets,coefs,*chebyU_Np1,D);
      chebyU_Np1->gaugeCond('R',0);chebyU_Np1->gaugeCond('L',0);
      norm_Np1=log2(chebyU_Np1->getNormFact())+norm_N+auxN;chebyU_Np1->setNormFact(1.);
      // Shift the pointers and delete the discarded one 
    }
      delete chebyU_Nm1;
      chebyU_Nm1=chebyU_N;norm_Nm1=norm_N;
      chebyU_N=chebyU_Np1;norm_N=norm_Np1;      
    

    }


    // Now compute the estimation of the state and the expectation values 
    if(k>0&&k%2==0){ // o.w. already done, and odd terms have 0 coeff (only for the delta!)
      MPS aux(rho0);
      vector<const MPS*> kets;
      kets.push_back(&aux);kets.push_back(chebyU_Np1);
      vector<complex_t> coefs;coefs.push_back(ONE_c);
      coeffK=getJacksonKPMcoeff(k,M);
      if((k/2)%2!=0) coeffK=-coeffK;
      coefs.push_back((2./M_PIl)*coeffK*pow(2,norm_Np1-normMPS)*ONE_c);
      contractor.optimizeSum(kets,coefs,rho0,D);
      rho0.gaugeCond('R',0);rho0.gaugeCond('L',0);
      normMPS+=log2(rho0.getNormFact());rho0.setNormFact(1.);
    

      trR=contractor.contract(rho0,idMPS);
      Eev=contractor.contract(rho0,hmps0);
      SxCev=contractor.contract(rho0,SxCenter);
      SyCev=contractor.contract(rho0,SyCenter);
      SzCev=contractor.contract(rho0,SzCenter);
      valComm=contractor.contract2(commMPO,rho0);
      Entropy=contractor.getEntropy(rho0);

   valSSz1=contractor.contract(rho0,SSz1);
   valSSz2=contractor.contract(rho0,SSz2);  


      out=new ofstream(outfname.data(),ios::app);
      *out<<k<<"\t"<<D<<"\t"<<scale<<"\t"<<normMPS
	  <<"\t"<<real(trR)<<"\t"<<imag(trR)
	  <<"\t"<<real(Eev)<<"\t"<<imag(Eev)
	//  <<"\t"<<real(Sxev)<<"\t"<<imag(Sxev)
	//  <<"\t"<<real(Szev)<<"\t"<<imag(Szev)
	  <<"\t"<<real(SxCev)/real(trR)
	//  <<"\t"<<real(SyCev)
	  <<"\t"<<real(SzCev)/real(trR)
	  <<"\t"<<real(valComm)/(scale*scale)
	  <<"\t"<<Entropy
	  <<"\t"<<real(valSSz1)
		       <<"\t"<<real(valSSz2)
	  <<endl;
      out->close();

      cout<<"k="<<k<<", <E>="<<Eev/trR //*pow(2,normMPS)*scale
	  <<", Sx="<<SxCev/trR //*pow(2,normMPS)*scale
	  <<", comm="<<valComm/(trR*trR)/(scale*scale) // *pow(2,2*normMPS)
	  <<endl;
    } // If the nr of the step is multiple of the desired saving frequency, save tmp files                                                                                                                    
  if((saveFreq>0&&k%saveFreq==0&&k>0)||(M-k<=1)){
    cout<<"About to save tmp files for step "<<k<<", normMPS="<<normMPS<<", norm_N="<<norm_N<<", norm_Nm1="<<norm_Nm1<<endl;
    tmpSave(rho0,normMPS,*chebyU_Nm1,norm_Nm1,*chebyU_N,norm_N,mpsdir,mpsfile,k,tmpFiles);

      
    
  }
  }
  // at the very end, delete the two left
  delete chebyU_Nm1,chebyU_N;
  
  }

void sumMPO(MPO& result,const MPO& mpo1,const MPO& mpo2){
  int L=mpo1.getLength();
  if(mpo2.getLength()!=L){
    cout<<"Error: incompatible MPO lengths for sum"<<endl;
    exit(1);
  }
  if(L==1){
    cout<<"Error: cannot do sum for MPO of length 1!!"<<endl;
    exit(1);
  }
  result.initLength(L);
 
  cout<<"Summing MPOs "<<mpo1<<" and "<<mpo2<<endl;
  for(int k=0;k<L;k++){
    // take both Operators involved
    const mwArray& op1=mpo1.getOpData(k);
    const mwArray& op2=mpo2.getOpData(k);
    Indices dims1=op1.getDimensions();
    Indices dims2=op2.getDimensions();
    if((dims1[0]!=dims2[0])||(dims1[2]!=dims2[2])){
      cout<<"Error: incompatible physical dimensions at position "<<k<<endl;
      exit(1);
    }
    int Dl=dims1[1]+dims2[1];
    if(k==0) Dl=1;
    int Dr=dims1[3]+dims2[3];
    if(k==L-1) Dr=1;
    cout<<"sumMPO("<<k<<") op1:"<<dims1<<", op2:"<<dims2<<"; new Dl="<<Dl<<", Dr="<<Dr<<endl;
    // Now fill in element by element
    mwArray blockOp(Indices(dims1[0],Dl,dims1[2],Dr));
    for(int l=0;l<Dl;l++){
      for(int r=0;r<Dr;r++){
	bool copy1=0;bool copy2=0;
	if((l<dims1[1])&&(r<dims1[3])) copy1=1;
	if((l>=dims1[1])&&(r>=dims1[3])) copy2=1;
	// special: edges
	if((k==0)&&(r<dims1[3])) copy1=1;
	if((k==0)&&(r>=dims1[3])) copy2=1;
	if((k==L-1)&&(l<dims1[1])) copy1=1;
	if((k==L-1)&&(l>=dims1[1])) copy2=1;
	// (operator element Dl, Dr): if Dl, Dr<dims1[2,3]
	//cout<<"For element "<<l<<","<<r<<" copy1="<<copy1
	//  <<" copy2="<<copy2<<endl;
	if(copy1+copy2){ // sth to copy
	  for(int i1=0;i1<dims1[0];i1++){
	    for(int i2=0;i2<dims1[2];i2++){
	      if(copy1) blockOp.setElement(op1.getElement(Indices(i1,l,i2,r)),Indices(i1,l,i2,r));
	      if(copy2){
		int indl=k==0?0:l-dims1[1];
		int indr=k==L-1?0:r-dims1[3];
		blockOp.setElement(op2.getElement(Indices(i1,indl,i2,indr)),Indices(i1,l,i2,r));
	      }
	    }
	  }
	}
      }
    }
    result.setOp(k,new Operator(blockOp),true);
  }

}

double getJacksonKPMcoeff(int k,int M){
  return ((M+1.-k)*cos(M_PIl*k/(M+1))+sin(M_PIl*k/(M+1))*cos(M_PIl/(M+1))/sin(M_PIl/(M+1)))/(M+1);
}


const string mpsFileName(const string& mpsdir,const string& baseName,const string& label,int cnt){
  stringstream s;
  s<<mpsdir<<"/"<<baseName<<"_"<<label;
  if(cnt>0) s<<"_"<<cnt;
  return s.str();
}

void tmpSave(MPS& mps,double norm_MPS,MPS& chebyU_Nm1,double norm_Nm1,MPS& chebyU_N,
             double norm_N,const string& baseDir,const string& baseName,int cnt,vector<string>& tmpFiles){
  // first create new tmpFile names                                                                                                                                                                         
  vector<string> existing(tmpFiles); // copy to remove at the end                                                                                                                                           
  tmpFiles.clear();
  // now save ony by one                                                                                                                                                                                    
  string newTmp=mpsFileName(baseDir,baseName,"mps",cnt);
  double exNorm=mps.getNormFact();
  mps.setNormFact(norm_MPS+log2(exNorm));
  mps.exportMPS(newTmp.data()); // saved containing log of norm!                                                                                                                                            
  mps.setNormFact(exNorm); // restore                                                                                                                                                                       
  tmpFiles.push_back(newTmp);

  newTmp=mpsFileName(baseDir,baseName,"Nm1",cnt);
  exNorm=chebyU_Nm1.getNormFact();
  chebyU_Nm1.setNormFact(norm_Nm1+log2(exNorm));
  chebyU_Nm1.exportMPS(newTmp.data()); // saved containing log of norm!                                                                                                                                     
  chebyU_Nm1.setNormFact(exNorm); // restore                                                                                                                                                                
  tmpFiles.push_back(newTmp);

  newTmp=mpsFileName(baseDir,baseName,"N",cnt);
  exNorm=chebyU_N.getNormFact();
  chebyU_N.setNormFact(norm_N+log2(exNorm));
  chebyU_N.exportMPS(newTmp.data()); // saved containing log of norm!                                                                                                                                       
  chebyU_N.setNormFact(exNorm); // restore                                                                                                                                                                  
  tmpFiles.push_back(newTmp);

  cout<<"Saved tmp files for step "<<cnt<<" to ";
  for(int p=0;p<tmpFiles.size();p++)
    cout<<tmpFiles[p]<<", ";
  cout<<endl;
  // now remove the old one                                                                                                                                                                                 
  if(existing.size()>0){
    for(int k=0;k<existing.size();k++){
      string tmpF=existing[k];
      if(tmpF.length()>0&&file_exists(tmpF)){
        remove(tmpF.data());
      }
    }
  }
}



int tmpRestore(MPS& mps,double& norm_MPS,MPS& chebyU_Nm1,double& norm_Nm1,MPS& chebyU_N,
	       double& norm_N,const string& mpsdir,const string& baseName,int M,vector<string>& tmpFiles){
  tmpFiles.clear();
  // first find if there are tmp files                                                                                                                                                                      
  int cnt=M;
  bool found=0;
  while(cnt>0&&!found){
    string mpsfile=mpsFileName(mpsdir,baseName,"mps",cnt);
    if(file_exists(mpsfile)){
      cout<<"findTmpFile found mpsfile="<<mpsfile<<endl;
      // Check also for the others                                                                                                                                                                          
      string Nm1file=mpsFileName(mpsdir,baseName,"Nm1",cnt);
      string Nfile=mpsFileName(mpsdir,baseName,"N",cnt);
      if(file_exists(Nm1file)&&file_exists(Nfile)){
        // Read the files and return                                                                                                                                                                        
        mps.importMPS(mpsfile.data());
        chebyU_Nm1.importMPS(Nm1file.data());
        chebyU_N.importMPS(Nfile.data());
        norm_MPS=mps.getNormFact();mps.setNormFact(1.);
        norm_Nm1=chebyU_Nm1.getNormFact();chebyU_Nm1.setNormFact(1.);
        norm_N=chebyU_N.getNormFact();chebyU_N.setNormFact(1.);
        tmpFiles.push_back(mpsfile);
        tmpFiles.push_back(Nm1file);
        tmpFiles.push_back(Nfile);
        found=true;
        return cnt;
      }
      else{
        cout<<"ERROR: Found a file for MPS but not the corresponding ones for Nm1 and N=> no tmp file used"<<endl;
        return 0;
      }
    }
    else{
      cnt--;
    }
  }
  if(!found) return 0;
  cout<<"ERROR: Should not be here"<<endl;
  exit(1);
}

bool constructInitialRho(int L,int d,int initSt,MPS& rho0){
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  complex_t dataY[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sig0=identityMatrix(d);
  mwArray sigZ(Indices(d,d),dataZ);
  mwArray sigX(Indices(d,d),dataX);
  mwArray sigY(Indices(d,d),dataY);
  bool isReal=false;  

  // Product states
  if(initSt<9){
    mwArray state;
    switch(initSt){
    case 1: // X+
      {
	state=.5*(sig0+sigX);
	isReal=true;
	break;
      }
    case 2: // Y+
      {
	state=.5*(sig0+sigY);
	break;
      }



  						\

    case 3: // Z+
      {
	state=.5*(sig0+sigZ);
	isReal=true;
	break;
      }
    default:
      cout<<"Error: unknown type onf intSt="<<initSt<<endl;
      exit(1);
    }
    state.reshape(Indices(d*d,1,1));

    rho0=MPS(L,1,d*d);
    for(int k=0;k<L;k++)
      rho0.setA(k,state);
  }
  if(initSt==9){
    // Special case: |Y+>+|Y->
    MPS auxP(L,1,d);auxP.setProductState(p_yplus);
    MPS auxM(L,1,d);auxM.setProductState(p_yminus);
    Contractor& contr=Contractor::theContractor();
    MPS tmp(L,2,d);tmp.setRandomState();
    vector<const MPS*>vecs;vector<complex_t> coeffs(2,ONE_c);
    vecs.push_back(&auxP);vecs.push_back(&auxM);
    contr.optimizeSum(vecs,coeffs,tmp,2);
    tmp.gaugeCond('L',1);
    tmp.gaugeCond('R',1);
    // And now transform to double MPS
    for(int k=0;k<L;k++){
      mwArray tens=tmp.getA(k).getA(); //d x Dl xDr
      int Dl=tens.getDimension(1);int Dr=tens.getDimension(2);
      tens.reshape(Indices(d*Dl*Dr,1));
      tens.multiplyRight(Hconjugate(tens));
      tens.reshape(Indices(d,Dl,Dr,d,Dl,Dr));
      tens.permute(Indices(1,4,2,5,3,6));
      tens.reshape(Indices(d*d,Dl*Dl,Dr*Dr));
      rho0.replaceSite(k,tens,0);
    }
    isReal=true;
 
  
  return isReal;
  }


  if(initSt==10){

     mwArray spinsiteup(Indices(2,1,1));           spinsiteup.fillWithZero();
  mwArray spinsitedown(Indices(2,1,1));         spinsitedown.fillWithZero();

  spinsiteup.setElement(ONE_c,Indices(0,0,0));
  spinsitedown.setElement(ONE_c,Indices(1,0,0));

  //  if(state == is_updown)
    
  /*    cout << "Constructing state |u>|0>|d>|0>..." << endl;

      rho0.clear();
      rho0=MPS(L,1,2);

      for(int i=0; i<L; i++)
	{
	  if((i%2)!=0)
	    rho0.setA(i,spinsitedown);
	  else
	    rho0.setA(i,spinsiteup);
	
	    }} */

    // else if(state == is_bgfield_superpos)
  
  {
      //Little error check                                                                                                                                                                                    
      if((L%2)!=0)
	{
	  cout << "State can only be constructed for an even number of sites" << endl;
	  cout << "Program will be aborted..." << endl;
	  exit(666);
	}
      cout << "Constructing superposition state" << endl;

      rho0.clear();
      rho0=MPS(L,2,d);

      mwArray First,Odd,Even,Last;
      First=mwArray(Indices(2,1,2));      First.fillWithZero();
      Odd=mwArray(Indices(2,2,2));        Odd.fillWithZero();
      Even=mwArray(Indices(2,2,2));       Even.fillWithZero();
      Last=mwArray(Indices(2,2,1));       Last.fillWithZero();

      //Set elements                                                                                                                                                                                          
      First.setElement(ONE_c,Indices(0,0,0));
      First.setElement(ONE_c,Indices(1,0,1));

      Odd.setElement(ONE_c,Indices(0,0,0));
      Odd.setElement(ONE_c,Indices(0,0,1));

      Even.setElement(ONE_c,Indices(1,0,0));
      Even.setElement(ONE_c,Indices(1,1,0));

      Last.setElement(ONE_c,Indices(0,1,0));
      Last.setElement(ONE_c,Indices(1,0,0));

      rho0.setA(0,First);
      rho0.setA(L-1,Last);

      for(int i=2; i<L; i++)
	{
	  if((i%2)==0)
	    rho0.setA(i-1,Even);
	  else
	    rho0.setA(i-1,Odd);
	}
      for(int k=0;k<L;k++){
 mwArray tens=rho0.getA(k).getA(); //d x Dl xDr
      int Dl=tens.getDimension(1);int Dr=tens.getDimension(2);
      tens.reshape(Indices(d*Dl*Dr,1));
      tens.multiplyRight(Hconjugate(tens));
      tens.reshape(Indices(d,Dl,Dr,d,Dl,Dr));
      tens.permute(Indices(1,4,2,5,3,6));
      tens.reshape(Indices(d*d,Dl*Dl,Dr*Dr));
      rho0.replaceSite(k,tens,0);
  
      }}
 	isReal=true;
	return isReal;


 }


  
}





 void computeTwoBodyCorrCenter(MPS& result,const MPS& rho,const mwArray& op1,
			       const mwArray& op2,int l){

   int d0=2;
   int L=rho.getLength();
   MPS aux(L,1,d0*d0);
   mwArray basic=identityMatrix(d0);
   basic.reshape(Indices(d0*d0,1,1));
   for(int k=0;k<L;k++){
     aux.setA(k,basic);
   }

   mwArray op1_=reshape(op1,Indices(d0*d0,1,1));
   op1_.conjugate();
   mwArray op2_=reshape(op2,Indices(d0*d0,1,1));
   op2_.conjugate(); // gets conjugated in the scalar product!                                                                                                                                              \
 \                                                                                                                                                                                                          
                                                                                                                                                                                                           \


   // Now, site by site, replace the identity at k by the operator op1,                                                                                                                                     \
 \                                                                                                                                                                                                          
                                                                                                                                                                                                           \

   // the identity in k+l by op2,                                                                                                                                                                           \
 \                                                                                                                                                                                                          
                                                                                                                                                                                                           \


   aux.setA(L/2,op1_); // operator on site k                                                                                                                                                                \
                                                                                                                                                                                                            
 aux.setA(L/2+l,op2_); // operator on site k+l                                                                                                                                                            \
                                                                                                                                                                                                            
 result=aux;


 }


