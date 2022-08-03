
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

void sumMPO(MPO& result,const MPO& mpo1,const MPO& mpo2);

void probeDoS(vector<double>& values,const string mpsprobelist,int L,int d,
	      const MPS& dosMPS,const MPS& wX,const MPO& hamil);


/** chebyDoS constructs an approximation to the DoS operator of a certain Ising Hamiltonian
    using the Chebyshev expansion. 
    Yet missing: probing the values with low variance states.
    TODO: Same for Heisenberg with disorder.
 */

int d=2;

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

  int L=props.getIntProperty("L");
  double J_=props.getDoubleProperty("J");
  double h_=props.getDoubleProperty("h");
  double g_=props.getDoubleProperty("g");
  double delta=props.getDoubleProperty("delta");
  if(delta<0) delta=0.02; 
  int D=props.getIntProperty("D");
  // int dpur=props.getIntProperty("dpur");
  // if(dpur<0) dpur=d;
  string outfname=props.getProperty("outputfile");
  string mpsfname=props.getProperty("mpsfile"); // file for the resulting DoS
  string mpsfnameW=props.getProperty("mpsfileW"); // file for the weight (can be reused!)
  bool app=(props.getIntProperty("append")!=0);
  int M=props.getIntProperty("M"); // LARGEST ORDER in Chebyshev exp.
  int Mw=props.getIntProperty("Mw"); // Largest order in the expansion of the weight
  if(Mw<0) Mw=M;
  double Emin=props.getDoubleProperty("Emin");
  double Emax=props.getDoubleProperty("Emax");

  string mpsprobelist=props.getProperty("mpsprobelist"); // file with a list of MPS to use to probe the DoS
  int freqProbe=props.getIntProperty("freqprobe"); // how often we compute the values with the MPS
  if(freqProbe<0) freqProbe=M/10; // default
  
  cout<<"Initialized arguments: L="<<L
      <<", J="<<J_
      <<", g="<<g_
      <<", h="<<h_
      <<", outfile="<<outfname
      <<", D="<<D<<endl;


  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;

  // The original Hamiltonian to get the energy
  IsingHamiltonian hamH0(L,d,J_,g_,h_);
  const MPO& hamil0=hamH0.getHMPO();
  double scale=1.;
  double offset=0.;
  if(Emin==-1&&Emax==-1){
    cout<<"First need to estimate the energy band to rescale "
	<<"and shift the Hamiltonian"<<endl;
    // First rescale the Hamiltonian:
    //hamH.registerTimeDependence();
    cout<<"Created the non-rescaled Hamiltonian"<<endl;
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
  }
    // Now the scale has to be such that E+alpha fits in [-1+delta,1+delta]
  //scale=2.*(1.-delta)/(Emax-Emin);
  //offset=-Emin*scale-1.+delta;
  scale=(1.-delta)/max(abs(Emax),abs(Emin));
  offset=0.;
  
  // Now the rescaled H, such that the spec is in [-1+delta,1+delta]
  IsingHamiltonian hamH(L,d,J_*scale,g_*scale,h_*scale,offset);
  //hamH.registerTimeDependence();
  cout<<"Created the rescaled Hamiltonian: scale="<<scale<<", offset="<<offset<<endl;
  const MPO& hamil_=hamH.getHMPO();
  int Dh=3; // bond dimension of H, just in case
  
  // Prepare the observables to compute

  MPS idMPS(L,1,d*d);
  for(int k=0;k<L;k++)
    idMPS.setA(k,reshape(identityMatrix(d),Indices(d*d,1,1)));

  MPS hamilMPS(L,Dh,d*d);
  MPSfromMPO(hamil_,hamilMPS);

  MPS hamil0MPS(L,Dh,d*d);
  MPSfromMPO(hamil0,hamil0MPS);

  MPS SxCenter(idMPS);
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  sigX.reshape(Indices(d*d,1,1));
  SxCenter.setA(L/2,sigX);//sigX.reshape(Indices(d,d));

  MPS SyCenter(idMPS);
  complex_t dataY[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  sigY.reshape(Indices(d*d,1,1));
  SyCenter.setA(L/2,sigY);//sigY.reshape(Indices(d,d));

  MPS SzCenter(idMPS);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  sigZ.reshape(Indices(d*d,1,1));
  SzCenter.setA(L/2,sigZ);//sigZ.reshape(Indices(d,d));
  
  MPO hamilId(L); // needed for the recursion
  extendMPO(hamil_,hamilId,d);
  //  hamilId.exportForMatlab("hId.m");

  MPS wx(idMPS);double normW(0.); // will be the log
  // first the weighting function
  if(mpsfnameW.size()>0&&file_exists(mpsfnameW)){ // read from file and use
    wx.importMPS(mpsfnameW.data());
    cout<<"Read weighting function from file "<<mpsfnameW<<endl;
    normW=wx.getNormFact();wx.setNormFact(1.); // trick: saving the log!
    wx.gaugeCond('R',0);wx.gaugeCond('L',0);
    normW+=log2(wx.getNormFact()); // hopefully unchanged
    wx.setNormFact(1.);
  }
  else{
    cout<<"First, approximate the weight as a polynomial to order "<<Mw<<endl;
    // first term id the 0-th power is just 1
    wx.gaugeCond('R',0);wx.gaugeCond('L',0);
    normW+=log2(wx.getNormFact()); // hopefully unchanged
    wx.setNormFact(1.);
    // I just apply powers of 2, so I prepare a MPO for this first
    MPO hamil2Id(L); // the double one to apply
    {
      const MPO* ptrs[]={&hamil_,&hamil_};
      MPO ham2(L); // auxiliary
      MPO::join(2,ptrs,ham2);// MPO for H^2
      extendMPO(ham2,hamil2Id,d);
    }
    MPS Hm(wx); // the first power (0) is Id
    double normM=normW;
    double cM=1.; // coefficient
    for(int m=2;m<=Mw;m+=2){ // only even powers
      // apply H^2
      MPS aux(Hm);
      contractor.optimize(hamil2Id,aux,Hm,D);
      Hm.gaugeCond('R',0);Hm.gaugeCond('L',0);
      normM+=log2(Hm.getNormFact());
      Hm.setNormFact(1.);
      // decide the right coefficient
      if(m==2) cM=-.5;
      else{ // use recurrence TODO: in log form, so I pull out the minus sign
	cM*=((m-3)/m);
      }
      // Now optimize the linear combination
      aux=wx;
      vector<const MPS*> vecs;vector<complex_t> coeffs;
      vecs.push_back(&aux);coeffs.push_back(ONE_c);
      vecs.push_back(&Hm);coeffs.push_back(pow(2,normM-normW)*cM*ONE_c);
      contractor.optimizeSum(vecs,coeffs,wx,D);
      wx.gaugeCond('R',0);wx.gaugeCond('L',0);
      normW+=log2(wx.getNormFact());
      wx.setNormFact(1.);
    }


    // at the end, save it
    if(mpsfnameW.size()>0){
      wx.setNormFact(normW); // trick to store the factor together!
      wx.exportMPS(mpsfnameW.data());
      wx.setNormFact(1.);
    }
  }



  //// TODO: What to write during the expansion????  
  // During the expansion, I only record, with certain frequency, the
  // result of probing the energy distribution with certain MPS, found in the given files.
  // This is done with method probeDoS
  
  // For output
  ofstream* out;
  if(mpsprobelist.size()>0){
    if(!file_exists(outfname)){
      out=new ofstream(outfname.data());
      *out<<"% Chebyshev expansion of the DoS, M="<<M<<endl;
      *out<<"% k\t J\t g\t h\t D\t E(k)\t <H^2(k)>\t <k|dos|k>, for k=1...."<<endl;
      if(!out->is_open()){
	cout<<"Error: impossible to open file "<<outfname<<
	  " for output"<<endl;
	exit(1);
      }
      out->close();delete out;
    }
  }
  
  // Now start the Chebyshev expansion computing U0(Id) and U1(X)
  MPS* chebyU_Nm1=new MPS(idMPS); 
  MPS* chebyU_N=new MPS(hamilMPS); // times 2, but I put this in the norm
  chebyU_Nm1->gaugeCond('R',0);chebyU_Nm1->gaugeCond('L',0);
  chebyU_N->gaugeCond('R',0);chebyU_N->gaugeCond('L',0);
  double norm_Nm1=log2(chebyU_Nm1->getNormFact());chebyU_Nm1->setNormFact(1.);
  double norm_N=1.+log2(chebyU_N->getNormFact());chebyU_N->setNormFact(1.);
  
  MPS rhoDoS(idMPS); // times Cs[0]=tr(Id)=2^L, but I put it in the norm
  rhoDoS.gaugeCond('R',0);rhoDoS.gaugeCond('L',0);
  double normDoS=L+log2(rhoDoS.getNormFact());
  rhoDoS.setNormFact(1.);
  
  // place for trace, energy, trace^2 and ev of X Y Z at the center
  // TODO!!!!! Sth else?????
  double trRho(0.);
  double trRho2(0.);
  double trRhoH(0.);
  double valX(0.);
  double valY(0.);
  double valZ(0.);
      
  // iteration
  for(int k=0;k<=M;k++){
    MPS* chebyU_Np1;double norm_Np1(0.);
    if(k==0){
      chebyU_Np1=chebyU_Nm1; 
      norm_Np1=norm_Nm1;
    }
    else if(k==1){
      chebyU_Np1=chebyU_N;
      norm_Np1=norm_N;
    }
    else{
      chebyU_Np1=new MPS(*chebyU_N);
      // TODO: Optimize the sum of ops acting on vecs!
      MPS aux(*chebyU_N);
      contractor.optimize(hamilId,*chebyU_N,aux,D);
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
      delete chebyU_Nm1;
      chebyU_Nm1=chebyU_N;norm_Nm1=norm_N;
      chebyU_N=chebyU_Np1;norm_N=norm_Np1;      
    }
    // Now compute the estimation of rhoDoS and maybe some expectation values?
    if(k>0){ // o.w. already done
      MPS aux(rhoDoS);
      vector<const MPS*> kets;
      kets.push_back(&aux);kets.push_back(chebyU_Np1);
      vector<complex_t> coefs;coefs.push_back(ONE_c);
      // The Ck coeff is the trace of the corresponding poly
      complex_t Ck=contractor.contract(*chebyU_Np1,idMPS); // times the norm factor!
      coefs.push_back(Ck*pow(2,2*norm_Np1-normDoS)); // included here (factor 2)
      contractor.optimizeSum(kets,coefs,rhoDoS,D);
      rhoDoS.gaugeCond('R',0);rhoDoS.gaugeCond('L',0);
      normDoS+=log2(rhoDoS.getNormFact());rhoDoS.setNormFact(1.);
      cout<<"Computed approximation for k="<<k<<endl;
    }

    // If appropriate, probe the distribution
    if(mpsprobelist.size()>0&&(k%freqProbe==0||k==M)){
      vector<double> values;
      probeDoS(values,mpsprobelist,L,d,rhoDoS,wx,hamil0);

      // and save values to file
      out=new ofstream(outfname.data(),ios::app);
      *out<<setprecision(15);
      *out<<k<<"\t"<<J_<<"\t"<<g_<<"\t"<<h_<<"\t"<<D<<"\t";
      for(int p=0;p<values.size();p++)
	*out<<values[p]<<"\t";
      *out<<endl;
      out->close();
      delete out;
    }

  }

  // at the very end, delete the two left
  delete chebyU_Nm1; delete chebyU_N;

  // And export the resulting rho, for comparison
  rhoDoS.setNormFact(normDoS); // trick!
  rhoDoS.exportMPS(mpsfname.data());
  
}


void probeDoS(vector<double>& values,const string mpsprobelist,int L,int d,
	      const MPS& dosMPS,const MPS& wX,const MPO& hamil){
  static vector<string> list;
  static vector<double> energies;
  static vector<double> vars;
  static MPO wXmpo(L);
  static int nr=0;
  static bool init=false;
  if(!init){ // read the names of the files
    ifstream input(mpsprobelist.data());
    if(!input.is_open()){
      cout<<"Error trying to read list of files from "<<mpsprobelist<<endl;
      exit(1);
    }
    string s;
    while(input.is_open()){
      getline(input,s);
      if(input.good()){
	if(s.length()>0&&s[0]!='%'){ // ignore comments and empty lines
	  list.push_back(s);
	  nr++;
	  cout<<"Read file name: "<<s<<endl;
	}
      }
      else{ // no more lines
	input.close();
	cout<<"Closed the list file"<<endl;
      }
    }
    MPOfromMPS(wX,wXmpo);
  }
  MPO dosMPO(L),auxMPO(L);
  MPOfromMPS(dosMPS,dosMPO);
  const MPO* ptrs[2]={&wXmpo,&dosMPO};
  MPO::join(2,ptrs,auxMPO);
  
  Contractor& contractor=Contractor::theContractor();
  for(int k=0;k<nr;k++){
    // Read and evaluate state nr
    MPS mps(L,1,d);
    if(file_exists(list[k].data())){
      mps.importMPS(list[k].data());
    }
    else{
      // Try tmp file with the same name?
      string tmpFilek=list[k]+".tmp";
      if(!file_exists(tmpFilek.data())){
	cout<<"Error! Cannot find file "<<list[k]<<" neither its tmp "<<tmpFilek<<endl;
	exit(1);	
      }
      mps.importMPS(tmpFilek.data());
    }
    complex_t val=contractor.contract(mps,auxMPO,mps);
    complex_t norm=contractor.contract(mps,mps);
    if(!init){
      complex_t Ek=contractor.contract(mps,hamil,mps);
      energies.push_back(real(Ek/norm));
      complex_t E2k=contractor.contract2(hamil,mps);
      vars.push_back(real(E2k/norm));
      cout<<"Probe state "<<k<<" has energy "<<Ek<<" and <H^2>="<<E2k<<", thus variance="
	  <<E2k-Ek*Ek<<endl;
    }
    values.push_back(energies[k]);
    values.push_back(vars[k]);
    values.push_back(real(val/norm));
  }
  init=true;
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
