#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "SpinMPO.h"
#include "Properties.h"
#include "misc.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#define OPTIONA 0

#include "IsingHamiltonian.h"

using namespace shrt;

/** Construct the initial state. Return true if the state is real! */
bool constructInitialRho(int L,int d,int initSt,MPS& rho0);

/** Prepare a name for the mps correspondng to step cnt */
const string mpsFileName(const string& mpsdir,const string& baseName,int cnt);

/** Save a tmp mps file, replacing the previous one */
void tmpSave(MPS& rho,const string& baseDir,const string& baseName,
	     int cnt,string& tmpFile);

/** Check for tmp files */
int findTmpFile(int M,const string& mpsdir,const string& baseName,string& mpsfile);

/** Build a sum MPO making blocks from two different MPOs */
void sumMPO(MPO& result,const MPO& mpo1,const MPO& mpo2);

/**

   diagonalEnsembleProj tries to find the infinite time limit of the time
   averaged RDM, by applying an operator that kills the off-diagonal
   terms in the energy basis.
   In this case, I iteratively substract the component in the direction of
   \f[[H,\cdot]^2 |\Psi\rangle\f]

   In this version, the ansatz for the diagonal ensemble is:
   |rho>=alpha |Id> + beta |H> +|rho'>
   where alpha and beta are fixed by the trace(1) and the energy value
   at the beginning, and the manipulation only affects rho'
   In order to do this, each step has to be done 
   with the constraint that the result, rho', is orthogonal to both
   |Id> and |H> vectors (in the exact case it would be, but the truncation 
is not expected to respect this property).

   Parameters are provided as a Properties file, where, in particular,
   the following should be included
   \param <L> (int) number of sites
   \param <J> (double) parameter \f$J\f$ of the Hamiltonian
   \param <g> (double) parameter \f$g\f$ of the Hamiltonian
   \param <h> (double) parameter \f$h\f$ of the Hamiltonian
   \param <initSt> (int) which initial rho to consider. Options are:
   1 - |X+>
   2 - |Y+>
   3 - |Z+>
   9 - |Y+>+|Y-> (real part of Y+)
   \param <M> (int) max number of steps
   \param <rate> (int) frequency (steps) to write results
   \param <D> (int) maximum bond dimension
   \param <Dpsi> (int) [optional] bond dimension used for the pure
                     state (default D)
   \param <outfname> (char*) name of the output file for the results
                     (will be replaced if it exists)
   \param <mpsdir> (char*) directory where to store the tmp (or final)
                     state, so that I can continue the evolution
   \param <mpsfile> (char*) base of the filename where to save or from
                     which to read (only if newInstance is false) the
		     tmp MPS. The full name will be changed to include
		     the time step it corresponds to. The basis of
		     thename should be such that we do not mix
		     evolutions with different parameters (e.g. cahnging step width).
   \param <savingFreq> (int) frequency (steps) to save the tmp MPS file
   \param <newInstance> (bool) whether to start a new evolution from
                         the inhomogeneous polarization distribution
			 (if not, the MPS in mpsfile is used)
 */

int d=2;

int main(int argc,const char* argv[]){
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
  int initSt=props.getIntProperty("initSt");
  int M=props.getIntProperty("M"); // max number of steps
  int rate=props.getIntProperty("rate");
  int D=props.getIntProperty("D");
  string outfname=props.getProperty("outputfile");
  string mpsdir=props.getProperty("mpsdir");
  string mpsfile=props.getProperty("mpsfile");
  int savingFreq=props.getIntProperty("savingFreq");
  int _aux_=props.getIntProperty("newInstance");
  bool newInstance=_aux_>0;

  // For output
  ofstream* out;
  Contractor& contractor=Contractor::theContractor();
  cout<<"Initialized Contractor"<<endl;
  IsingHamiltonian hamH(L,d,J_,g_,h_);
  //hamH.registerTimeDependence();
  cout<<"Created the Hamiltonian"<<endl;

  // Prepare the observables to compute
  // The identity, to compute the trace
  MPS Id(L,1,d*d);
  mwArray id2=identityMatrix(d);id2.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++)
    Id.setA(k,id2);
  // The Hamiltonian to get the energy
  MPS hmps(L,1,d*d);
  const MPO& hamil=hamH.getHMPO();
  MPSfromMPO(hamil,hmps);

  //Commutator superoperator
  MPO commMPO(L);
  hamH.getCommutatorMPO(commMPO,0.); //commutator constructed
  cout<<"Commutator constructed"<<endl;

  // Sum of sigmaz/sigmax operators 
  MPS SxSum(L,2,d*d),SzSum(L,2,d*d);
  {
    MPO Saux(L);
    SpinMPO::getSxMPO(L,d,Saux);
    MPSfromMPO(Saux,SxSum);
    SpinMPO::getSzMPO(L,d,Saux);
    MPSfromMPO(Saux,SzSum);
  }
  // And I will also need the sigmaz/sigmay/sigmax operator on central site, which I could compute ony by one, but I assemble it here
  MPS SxCenter(Id);
  complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
  mwArray sigX(Indices(d,d),dataX);
  sigX.reshape(Indices(d*d,1,1));
  SxCenter.setA(L/2,sigX);sigX.reshape(Indices(d,d));

  MPS SyCenter(Id);
  complex_t dataY[]={ZERO_c,I_c,-1.*I_c,ZERO_c};
  mwArray sigY(Indices(d,d),dataY);
  sigY.reshape(Indices(d*d,1,1));
  SyCenter.setA(L/2,sigY);sigY.reshape(Indices(d,d));

  MPS SzCenter(Id);
  complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-1.*ONE_c};
  mwArray sigZ(Indices(d,d),dataZ);
  sigZ.reshape(Indices(d*d,1,1));
  SzCenter.setA(L/2,sigZ);sigZ.reshape(Indices(d,d));

  // The operator to be applied is the square commutator, so I can construct it as a joined MPO

  MPO comm2(L);
  {
    const MPO* mpos[2]={&commMPO,&commMPO};
    MPO::join(2,mpos,comm2);
  }

  // Initial state rho0
  MPS rho0(L,1,d*d);
  // Now if we start from scratch, initialize the MPS. Otherwise, read
  // the tmp file
  
  // In this version, even if I read from file, I need the initial
  // value of the energy in order to compute the corresponding
  // parameter of the |H> vector in the final state. So I always
  // generate the initial state and later, if possible, replace it by
  // file data.
  bool isReal=constructInitialRho(L,d,initSt,rho0);
  cout<<"Created state from scratch"<<endl;
  complex_t E0=contractor.contract(rho0,hmps);
  complex_t trH2=contractor.contract(hmps,hmps);
  // cout<<"Right now, trace of rho:"<<contractor.contract(rho0,Id)<<", energy "<<E0
  //     <<" and comm "<<contractor.contract2(commMPO,rho0)<<endl;
  // cout<<", and trace(H^2) "<<trH2<<", while d^L="<<pow(d,L)<<", and (J^2*(L-1)+g^2*L+h^2*L) "
  //     <<(J_*J_*(L-1)+g_*g_*L+h_*h_*L)<<" so the product: "<<(J_*J_*(L-1)+g_*g_*L+h_*h_*L)*pow(d,L)<<endl;
  // cout<<" and <rho|rho>="<<contractor.contract(rho0,rho0)<<endl;
    
  string tmpMPSfile;
  int cnt=0;
  double normR; // the log of the norm
  if(!newInstance){
    cnt=findTmpFile(M,mpsdir,mpsfile,tmpMPSfile);
    if(cnt>0){
      rho0.importMPS(tmpMPSfile.data());
      cout<<"Recovered state for cnt="<<cnt<<" at file "<<tmpMPSfile<<endl;
      // I use a trick and store the log of the norm!
      normR=rho0.getNormFact();
      rho0.setNormFact(1.);
      //normR=log2(real(contractor.contract(rho0,rho0)));
    }
  }
  if(newInstance||cnt==0){ // init from scratch
    //    constructInitialRho(L,d,initSt,rho0);
    //cout<<"Created state from scratch"<<endl;
    // I need to prepare the orthogonal part of rho, i.e. substract |Id> and |H>
    vector<const MPS*> vecs;vector<complex_t> coeffs;
    MPS auxId(Id);MPS auxH(hmps);
    for(int pos=0;pos<L;pos++){
      auxId.setA(pos,(1./d)*auxId.getA(pos).getA());
      auxH.setA(pos,(1./d)*auxH.getA(pos).getA());
    }
    MPS aux0(rho0); // the full initial state
    vecs.push_back(&aux0);
    vecs.push_back(&auxId);vecs.push_back(&auxH);
    coeffs.push_back(ONE_c); 
    coeffs.push_back(-1.*ONE_c); // coeff of the Id, with the factor 1/d^L included
    coeffs.push_back(-1.*(E0/(J_*J_*(L-1)+g_*g_*L+h_*h_*L))*ONE_c); // coeff of the H, which is E/tr(H^2) with the factor 1/d^L included in the vector
    contractor.optimizeSum(vecs,coeffs,rho0,max(D,5)); // I want this to be exact        
    //    normR=real(contractor.contract(rho0,rho0));
    normR=.5*log2(real(contractor.contract(rho0,rho0)));

    // cout<<"And after orthogonalization trace of rho:"<<contractor.contract(rho0,Id)
    // 	<<", energy "<<contractor.contract(rho0,hmps)
    // 	<<" and comm "<<contractor.contract2(commMPO,rho0)
    // 	<<" and <rho|rho>="<<normR<<endl;
    out=new ofstream(outfname.data());
    // Write header
    *out<<"% The first line (cnt=-1) contains values of trace, E, and exp values of the component \n"
	<<"% of sigma which is fixed, (1/dim)*Id +E/(tr H^2)*H . The line for cnt=0 contains already \n"
	<<"% values for rho_tilde, the same as the remaining evolution. Notice that the (Euclidean) \n"
	<<"% norm of rho_tilde, written in the 5th column, is not conserved, but the trace of rho\n"
	<<"% and of rho_tilde is."<<endl;
    *out<<"%"<<endl;
    *out<<"% cnt\t D\t log(tr(rho^2))\t tr(rho)\t tr(rho H)\t tr(rho Sx(L/2))\t tr(rho(Sy(L/2))\t tr(rho Sz(L/2))\t <rho|Hc^T.Hc|rho>\t "<<endl;
    *out<<setprecision(12);
    // And now write the contributions of the substracted part to observables: it is the value at t=0, and then it needs to be added, as
    // this is always constant. So I will label this time with -1!!
    complex_t trR0=contractor.contract(aux0,Id);
    complex_t normR0=contractor.contract(aux0,aux0);
    complex_t Eev0=contractor.contract(aux0,hmps);
    // complex_t Sxev=contractor.contract(rho0,SxSum);
    // complex_t Szev=contractor.contract(rho0,SzSum);
    complex_t SxCev0=contractor.contract(aux0,SxCenter);
    complex_t SyCev0=contractor.contract(aux0,SyCenter);
    complex_t SzCev0=contractor.contract(aux0,SzCenter);
    complex_t valComm0=contractor.contract2(commMPO,aux0);

    complex_t trR=contractor.contract(rho0,Id);
    complex_t Eev=contractor.contract(rho0,hmps);
    // complex_t Sxev=contractor.contract(rho0,SxSum);
    // complex_t Szev=contractor.contract(rho0,SzSum);
    complex_t SxCev=contractor.contract(rho0,SxCenter);
    complex_t SyCev=contractor.contract(rho0,SyCenter);
    complex_t SzCev=contractor.contract(rho0,SzCenter);
    complex_t valComm=contractor.contract2(commMPO,rho0);
    
    out=new ofstream(outfname.data(),ios::app);
    *out<<setprecision(12);
    *out<<-1<<"\t"<<D<<"\t"<<log2(real(normR0))
	<<"\t"<<real(trR0-trR)<<"\t"<<imag(trR0-trR)
	<<"\t"<<real(Eev0-Eev)<<"\t"<<imag(Eev0-Eev)
      //  <<"\t"<<real(Sxev)<<"\t"<<imag(Sxev)
      //  <<"\t"<<real(Szev)<<"\t"<<imag(Szev)
	<<"\t"<<real(SxCev0-SxCev)
	<<"\t"<<real(SyCev0-SyCev)
	<<"\t"<<real(SzCev0-SzCev)
	<<"\t"<<real(valComm-valComm0)
	<<endl;
    out->close();
    delete out;
  }
  rho0.increaseBondDimension(D);
  rho0.gaugeCond('R',1);
  rho0.gaugeCond('L',1);


  // And now iterate the projections
  cout<<"=== STARTING TIME PROJECTIONS FROM STEP "<<cnt<<" ==="<<endl;
  // These guys should be zero!
  complex_t trR=contractor.contract(rho0,Id);
  complex_t Eev=contractor.contract(rho0,hmps);
  //  complex_t normR=contractor.contract(rho0,rho0);
  // complex_t Sxev=contractor.contract(rho0,SxSum);
  // complex_t Szev=contractor.contract(rho0,SzSum);
  complex_t SxCev=contractor.contract(rho0,SxCenter);
  complex_t SyCev=contractor.contract(rho0,SyCenter);
  complex_t SzCev=contractor.contract(rho0,SzCenter);
  complex_t valComm=contractor.contract2(commMPO,rho0);

  out=new ofstream(outfname.data(),ios::app);
  *out<<cnt<<"\t"<<D<<"\t"<<normR
      <<"\t"<<real(trR)<<"\t"<<imag(trR)
      <<"\t"<<real(Eev)<<"\t"<<imag(Eev)
    //  <<"\t"<<real(Sxev)<<"\t"<<imag(Sxev)
    //  <<"\t"<<real(Szev)<<"\t"<<imag(Szev)
      <<"\t"<<real(SxCev)
      <<"\t"<<real(SyCev)
      <<"\t"<<real(SzCev)
      <<"\t"<<real(valComm)
      <<endl;
  out->close();

  // For the orthogonalization
  vector<const MPS*> vecsOrth;vecsOrth.push_back(&Id);vecsOrth.push_back(&hmps);
  vector<double> penOrth(2,50); // we need to adjust this penalty!!
  
  int stepsLeft=M-cnt;
  clock_t start,finish;
  while(stepsLeft>0){
    MPS aux(rho0);
    // cout<<"About to call orthogonalize with state rho0:"<<rho0<<" and op comm2:"<<comm2<<
    //   " and orthog wrt Id:"<<Id<<" and H:"<<hmps<<endl;
    // Optimize the comm2 on the rhoTilde, with orthog to Id and H
    contractor.orthogonalize(aux,D,comm2,rho0,vecsOrth,penOrth);
    cout<<"Orthogonalized Op on state"<<endl;
    
    // Now substract the component:
    // Normalize and compute the corresponding proj of rho0
    aux.gaugeCond('R',1);    
    complex_t alphaP2=contractor.contract(rho0,aux);
    // Substract and orthogonalize
    MPS aux0(rho0); // for the sum
    vector<const MPS*> vecs;vecs.push_back(&aux0);vecs.push_back(&aux);
    vector<complex_t> coeffs(2,ONE_c);coeffs[1]=-alphaP2;    
    contractor.orthogonalizeSum(rho0,D,vecs,coeffs,vecsOrth,penOrth);
    //    cout<<"Orthogonalized wrt previous result"<<endl;
    // Second option: instead of substracting, directly orthogonalize and set the right norm factor!
    //vector<const MPS*> vecs;vecs.push_back(&Id);vecs.push_back(&hmps);vecs.push_back(&aux);
    //vector<double> pens(3,50);
    //contractor.orthogonalize(rho0,D,vecs,pens);

    
    cnt++;stepsLeft--;

    // Q: what should be the norm of the remaining part????
    rho0.gaugeCond('R',1);
    normR+=log2(1-real(conjugate(alphaP2)*alphaP2));
    //    normR*=.5*(1.+real(contractor.contract(aux1,aux2)));


    // Now save results
    trR=contractor.contract(rho0,Id);
    Eev=contractor.contract(rho0,hmps);
    //    Sxev=contractor.contract(rho0,SxSum);
    //    Szev=contractor.contract(rho0,SzSum);
    SxCev=contractor.contract(rho0,SxCenter);
    SyCev=contractor.contract(rho0,SyCenter);
    SzCev=contractor.contract(rho0,SzCenter);
    valComm=contractor.contract2(commMPO,rho0);
    out=new ofstream(outfname.data(),ios::app);
    *out<<setprecision(12);
    *out<<cnt<<"\t"<<D<<"\t"<<normR
	<<"\t"<<real(trR)<<"\t"<<imag(trR)
	<<"\t"<<real(Eev)<<"\t"<<imag(Eev)
      //	<<"\t"<<real(Sxev)<<"\t"<<imag(Sxev)
      //	<<"\t"<<real(Szev)<<"\t"<<imag(Szev)
    <<"\t"<<real(SxCev)
    <<"\t"<<real(SyCev)
    <<"\t"<<real(SzCev)
    <<"\t"<<real(valComm)
	<<endl;
    out->close();
    delete out;
    // And if necessary, store tmpfile

    if(cnt>0&&(cnt%savingFreq==0||cnt==M)){
      cout<<"Saving backup at cnt="<<cnt<<endl;
      rho0.setNormFact(normR); // to have it stored: TRICK!!!! It is the log of the norm squared!
      tmpSave(rho0,mpsdir,mpsfile,cnt,tmpMPSfile);
      rho0.setNormFact(1.); // to continue
 }
  }

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
  }
  
  return isReal;
}
