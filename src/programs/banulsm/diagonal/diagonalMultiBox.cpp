
#include <math.h>
#include <iomanip>

#include "misc.h"
#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"
#include "Properties.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"

using namespace shrt;


void initialState(MPS& mps,int intSt,int D);

double getJacksonKPMcoeff(int k,int M);

/** diagonalMultiBox constructs a Chebyshev poly approx to the box
    function for the different intervals in the energy band, and uses
    this to estmate expectation values in the diagonal ensemble
    corresponding to a certain initial state. 

    \TODO This first version is not saving anything for later use, but
    it should be done.

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
  double J_=props.getDoubleProperty("J");
  double h_=props.getDoubleProperty("h");
  double g_=props.getDoubleProperty("g");
  double delta=props.getDoubleProperty("delta");
  if(delta<0) delta=0.2; 
  int D = props.getIntProperty("D");  
  string outfname=props.getProperty("output");
  //  string mpsdir=props.getProperty("mpsdir");
  // string mpsfile=props.getProperty("mpsfile");
  bool app=1;
  app=!(props.getIntProperty("append")==0); // replace only if 0 explicitly given (ow -1)
  int M=props.getIntProperty("M"); // LARGEST ORDER of Chebyshev poly
  double Emin=props.getDoubleProperty("Emin");
  double Emax=props.getDoubleProperty("Emax");
  string efname=props.getProperty("energiesfile");
  int initSt=props.getIntProperty("initSt"); // initial state to use:
  //                               +1=X+; -1=X-; (2,3 for Y, Z)
  //                               +4,+5,+6 are the staggered(+-) for X, Y and Z and -4,-5,-6 the -+
  //                               +7 is a (complex) random state and +8 a real random (in which cases, the random seed may be given)
  //  int randSeed=props.getIntProperty("seed");
  int nrBins=props.getIntProperty("Nbins");
  if(nrBins<0) nrBins=40; // default
  
  // int rate=props.getIntProperty("rate");
  // if(rate<0) rate=1;
  // int saveFreq=props.getIntProperty("saveFreq");
  // if(saveFreq<0) saveFreq=M;
  // saveFreq=lcm(rate,saveFreq);
  // int lastSavedCnt=-1;

  cout<<"Initialized arguments: L="<<L
      <<", J="<<J_
      <<", g="<<g_
      <<", h="<<h_
      <<", outfile="<<outfname
      <<", app="<<app
      <<", D="<<D<<endl;
  
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  
  // The original Hamiltonian (I'll also need the rescaled one for the Chebyshev expansion)
  IsingHamiltonian hamH0(L,d,J_,g_,h_);
  const MPO& hamil0=hamH0.getHMPO();
    // The original Hamiltonian to get the energy

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
  scale=(1.-delta)/max(abs(Emax),abs(Emin));
  offset=0.;
  
  // Now the rescaled H, such that the spectrum of [H,.] is in [-1+delta,1+delta]
  IsingHamiltonian hamH(L,d,J_*scale,g_*scale,h_*scale,offset);
  //hamH.registerTimeDependence();
  cout<<"Created the rescaled Hamiltonian: scale="<<scale<<", offset="<<offset<<endl;
  const MPO& hamil_=hamH.getHMPO();
  int Dh=3; // bond dimension of H, just in case

    // Initial state

  MPS mps(L,D,d);
  initialState(mps,initSt,D);
  cout<<"Initialized state ("<<initSt<<") from scratch, norm "<<contractor.contract(mps,mps)<<endl;
  mps.gaugeCond('R',1);
  mps.gaugeCond('L',1);
  
  double normMPS=0.;

  cout<<"Initial state has <H>="<<contractor.contract(mps,hamil0,mps)
      <<", <H^2>="<<contractor.contract2(hamil0,mps)<<endl;

  // Since there will be fewer bins than polynomials, it makes more
  // sense to compute all the boxes and keep them, although then for
  // different binning we'll need to repeat the full calculations.
  
  vector<MPS*> boxN;
  //  vector<double> normN; // norms of the resulting vectors, but saved as their log2

  // Now start the Chebyshev expansion computing T-1(0) and T0(Id) on the state
  MPS* chebyT_Nm1=new MPS(mps); // poly(-1)=0
  // Set the first tensor to 0=> ignore norm factor in this case
  chebyT_Nm1->gaugeCond('R',0);chebyT_Nm1->gaugeCond('L',1);
  chebyT_Nm1->setA(0,ZERO_c*chebyT_Nm1->getA(0).getA());
  double norm_Nm1=0.;
  
  MPS* chebyT_N=new MPS(mps); // normalized already by gaugeCond
  double norm_N=0.;//-log2(M_PIl); // But I need to put in front the coeff of T_0:1/pi=> I do it in the vector cN

  // iteration
  for(int k=0;k<=M;k++){
    MPS* chebyT_Np1;double norm_Np1(0.);
    if(k==0){
      chebyT_Np1=chebyT_N; 
      norm_Np1=norm_N;
    }
    else{
      chebyT_Np1=new MPS(*chebyT_N);
      // TODO: Optimize the sum of ops acting on vecs!
      MPS aux(*chebyT_N);
      contractor.optimize(hamil_,*chebyT_N,aux,D);
      aux.gaugeCond('R',0);aux.gaugeCond('L',0);
      double auxN=log2(aux.getNormFact());aux.setNormFact(1.);
      vector<const MPS*> kets;
      kets.push_back(&aux);kets.push_back(chebyT_Nm1);
      vector<complex_t> coefs;
      if(k==1) coefs.push_back(ONE_c); // special trick to get the T1
      else coefs.push_back(2.*ONE_c);
      coefs.push_back(-pow(2,norm_Nm1-norm_N-auxN)*ONE_c);
      contractor.optimizeSum(kets,coefs,*chebyT_Np1,D);
      chebyT_Np1->gaugeCond('R',0);chebyT_Np1->gaugeCond('L',0);
      norm_Np1=log2(chebyT_Np1->getNormFact())+norm_N+auxN;chebyT_Np1->setNormFact(1.);
      // Shift the pointers but do not delete, because I am keeping all of them
      delete chebyT_Nm1;
      chebyT_Nm1=chebyT_N;norm_Nm1=norm_N;
      chebyT_N=chebyT_Np1;norm_N=norm_Np1;      
    }
    //cout<<"Norm of poly nr "<<k<<" :"<<norm_Np1<<" ("<<pow(2,norm_Np1)<<")"<<endl;
    // Use the new poly to update all the boxes (for k=0, initialize them)
    double Nm=k==0?M_PIl:.5*M_PIl;
    double coeffKPMm=getJacksonKPMcoeff(k,M);
    double wBin=(Emax-Emin)/nrBins; // width of the bin
    for(int ib=0;ib<nrBins;ib++){
      double x1=ib==0?max(Emin-.5*wBin,-1./scale):Emin+ib*wBin;
      double x2=ib==nrBins-1?min(1./scale,Emax+.5*wBin):Emin+(ib+1)*wBin;
      x1*=scale;x2*=scale; // in energy scale
      //      double x1=.5*acos(scale*(Emin+ib*(Emax-Emin)/nrBins));
      // double x2=.5*acos(scale*(Emin+(ib+1)*(Emax-Emin)/nrBins));
      // if(ib==0){// for the first bin, put the lower limit further
      // 	x1=.5*acos(max(scale*(Emin-.5*(Emax-Emin)/nrBins),-1.));
      // }
      // if(ib==nrBins-1){
      // 	x2=.5*acos(min(scale*(Emax+.5*(Emax-Emin)/nrBins),1.));
      // }
      double cm;
      if(k==0){
	cm=acos(x1)-acos(x2);
	// cm=2.*(x1-x2)*coeffKPMm/Nm;
      }
      else{
	cm=2.*cos(.5*k*(acos(x1)+acos(x2)))*sin(.5*k*(acos(x1)-acos(x2)))/k;
	//cos(k*(x1+x2))*sin(k*(x1-x2))*2.*coeffKPMm/(Nm*k);
      }
      cm*=coeffKPMm/Nm;
      // cout<<"Step "<<k<<" of Cheby, treating box nr "<<ib<<" ["<<x1<<","<<x2<<"]"
      // 	  <<", cm="<<cm<<endl;    
      MPS* box;
      if(k==0){
	box=new MPS(*chebyT_Np1);
	// NOTICE: using a trick: keep log2 of norm instead of norm in
	// the factor, so I need to be careful when summing and
	// combining MPS
	box->setNormFact(norm_Np1+log2(chebyT_Np1->getNormFact())+log2(cm));
	// Factor in front: coef(m,k)*gk/Nk
	//box->setNormFact(pow(2,norm_Np1)*cm/Nm);
	// insert it in the vector
	boxN.push_back(box);
	//cout<<"For step "<<k<<" of box "<<ib<<" norm: "<<pow(2,boxN[ib]->getNormFact())<<endl;
      }
      else{
	box=boxN[ib];
	double logNb=box->getNormFact();box->setNormFact(1.);
	MPS aux(*box);
	vector<const MPS*> kets;
	kets.push_back(&aux);kets.push_back(chebyT_Np1);
	vector<complex_t> coefs;coefs.push_back(ONE_c);
	coefs.push_back(cm*pow(2,norm_Np1-logNb)*ONE_c);
	contractor.optimizeSum(kets,coefs,*box,D);
	box->gaugeCond('R',0);box->gaugeCond('L',0);
	logNb+=log2(box->getNormFact());box->setNormFact(logNb);
	//cout<<"For step "<<k<<" of box "<<ib<<" norm: "<<pow(2,logNb)<<", norm_Np1="<<norm_Np1<<endl;
	//if(k==1&&ib==nrBins-1) exit(1);
      }
    }
  }

  // Now that I have all ofthe boxes, compute the observables

  // Prepare the observables to compute
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  sigX.reshape(Indices(d,1,d,1)); Operator opX(sigX);
  complex_t dataY[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  sigY.reshape(Indices(d,1,d,1)); Operator opY(sigY);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  sigZ.reshape(Indices(d,1,d,1)); Operator opZ(sigZ);
  Operator idOp(reshape(identityMatrix(d),Indices(d,1,d,1)));

  MPO mpo(L); // I'll just set the op where I need it
  for(int p=0;p<L;p++) mpo.setOp(p,&idOp,false);

  double valX(0),valY(0),valZ(0); // estimates of local ev at the diag
  double totE(0),totProb(0);
  
  if(!file_exists(outfname.data())) app=0; // also if it does not exist (new)
  ofstream* out;
  if(!app){
    out=new ofstream(outfname.data());
    *out<<"% nBin\t D\t Emin\t Emax\t nrBins\t M\t lowE\t highE\t weightBin\t re(Sx_bin)\t re(Sy_bin)\t re(Sz_bin)";
    *out<<endl;
  }
  else{
    out=new ofstream(outfname.data(),ios::app);
  }
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  *out<<setprecision(10);
  for(int m=0;m<nrBins;m++){
    // weight of the box: the norm squared
    MPS* box=boxN[m];
    double logNb=box->getNormFact();box->setNormFact(1.);
    stringstream str;str<<"box"<<m+1<<".m";
    box->exportForMatlab(str.str().data());
    double probNb=pow(2,2*logNb)*real(contractor.contract(*box,*box));
    mpo.setOp(L/2,&opX,false);
    complex_t xNb=pow(2,2*logNb)*contractor.contract(*box,mpo,*box);
    mpo.setOp(L/2,&opY,false);
    complex_t yNb=pow(2,2*logNb)*contractor.contract(*box,mpo,*box);
    mpo.setOp(L/2,&opZ,false);
    complex_t zNb=pow(2,2*logNb)*contractor.contract(*box,mpo,*box);
    complex_t ENb=pow(2,2*logNb)*contractor.contract(*box,hamil0,*box);
    //mpo.setOp(L/2,&idOp,false);
    cout<<"For box "<<m<<", <Sx>="<<xNb<<", weight="<<pow(2,2*logNb)<<", "<<contractor.contract(*box,mpo,*box)<<endl;
    // Save values to file
    double loE=m==0?Emin-.5*(Emax-Emin)/nrBins:Emin+m*(Emax-Emin)/nrBins;
    double hiE=m==nrBins-1?Emax+.5*(Emax-Emin)/nrBins:Emin+(m+1)*(Emax-Emin)/nrBins;
    *out<<m<<"\t"<<D<<"\t"<<Emin<<"\t"<<Emax<<"\t"<<nrBins<<"\t"<<M<<"\t"<<loE<<"\t"<<hiE<<"\t";
    *out<<probNb<<"\t"<<real(xNb)<<"\t"<<real(yNb)<<"\t"<<real(zNb)<<"\t"<<real(ENb)<<endl;
    valX+=real(xNb);
    valY+=real(yNb);
    valZ+=real(zNb);
    totE+=real(ENb);
    totProb+=probNb;
  }
  cout<<"Final results for M="<<M<<", nrBins="<<nrBins<<": "<<endl;
  cout<<"\t<X>="<<valX<<endl;
			cout<<"\t<Y>="<<valY<<endl;
					      cout<<"\t<Z>="<<valZ<<endl;
cout<<"\t<E>="<<totE<<endl;
		      cout<<"\ttotProb="<<totProb<<endl;
  out->close();
  delete out;

  // And now I should clean up memory
  
}

double getJacksonKPMcoeff(int k,int M){
  return ((M+1.-k)*cos(M_PIl*k/(M+1))+sin(M_PIl*k/(M+1))*cos(M_PIl/(M+1))/sin(M_PIl/(M+1)))/(M+1);
}



const string mpsFileName(const string& mpsdir,const string& baseName,int cnt){
  stringstream s;
  s<<mpsdir<<"/"<<baseName;
  if(cnt>0) s<<"_"<<cnt;
  return s.str();
}

void tmpSave(MPS& rho,const string& baseDir,const string& baseName,int cnt,string& tmpFile){
  // first create a new tmpFile name
  string newTmp=mpsFileName(baseDir,baseName,cnt);
  rho.exportMPS(newTmp.data());
  // now remove the old one
  if(tmpFile.length()>0&&file_exists(tmpFile)){
    remove(tmpFile.data());
  }
  // and replace the string
  tmpFile=newTmp;
}

int findTmpFile(int M,const string& mpsdir,const string& baseName,string& mpsfile){
  int cnt=M;
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
  if(!found) return 0;
  cout<<"ERROR: Should not be here"<<endl;
  exit(1);
}

void initialState(MPS& mps,int initSt,int D){
  int L=mps.getLength();
  if(initSt==+1||abs(initSt)==4) mps.setProductState(p_xplus);
  else{
    if(initSt==-1) mps.setProductState(p_xminus);
    else{
      if(initSt==+2||abs(initSt)==5) mps.setProductState(p_yplus);
      else{
	if(initSt==-2) mps.setProductState(p_yminus);
	else{
	  if(initSt==+3||abs(initSt)==6) mps.setProductState(p_zero);
	  else{
	    if(initSt==-3) mps.setProductState(p_one);
	    else{ // any other number, I start random
	      cout<<"Initial state a non-translational random MPS with D="<<mps.getBond()<<endl;
	      mps.setRandomState(); // This creates it with complex values
	      // If real, discard imaginary parts
	      if(initSt==8){
		for(int pos=0;pos<L;pos++){
		  mwArray aux=mps.getA(pos).getA();
		  mps.replaceSite(pos,.5*(aux+conjugate(aux)));
		}
	      }
	      mps.gaugeCond('R',1);
	    }
	  }
	}
      }
    }
  }
  mps.increaseBondDimension(D);
  // Now for the staggered, change one every two by applying a single site op (sigma_z to x,y sigma_y to z)  
  if(abs(initSt)>=4&&abs(initSt)<=6){
    int firstPos=initSt<0?0:1; //+- changes first site 1
    complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
    mwArray oper;
    if(abs(initSt)<6) oper=mwArray(Indices(d,d),dataz);
    else oper=mwArray(Indices(d,d),datax);
    for(int k=firstPos;k<L;k+=2){
	mps.applyLocalOperator(k,oper);
    }
  }
}
