
#include <math.h>
#include <iomanip>

#include "MPS.h"
#include "Contractor.h"
#include "DoubleOperator.h"
#include "Properties.h"

#include "HeisenbergHamiltonian.h"
#include "MPOMultiplierHermitian.h"
#include "SpinMPO.h"
#include <cmath>
#include <unistd.h>

#include "misc.h"

//#include "quicksort.h"

using namespace std;
using namespace shrt;

#define MAXLEN 120

int d=2;
double tol=1E-8;

/** Auxiliary object, level-energy, to sort */
class Level{
public:
  int lev;
  double energy;
  Level(int l, double E):lev(l),energy(E){};
  Level(const Level& l2):lev(l2.lev),energy(l2.energy){};
  bool operator<=(const Level& l2){return energy<=l2.energy;}
  bool operator<(const Level& l2){return energy<l2.energy;}
  bool operator==(const Level& l2){return energy==l2.energy;}
  friend bool operator<=(const Level& l1,const Level& l2){return l1.energy<=l2.energy;}
  friend bool operator<(const Level& l1,const Level& l2){return l1.energy<l2.energy;}
  friend bool operator==(const Level& l1,const Level& l2){return l1.energy==l2.energy;}
};



/** Read a GS MPS from a file (needed input)
    and compute the spectrum of the effective local Hamiltonian for a certain number of sites over the whole chain.
    In the XY case, I am also computing the total Sz of the resulting "excitations".
    This version tries to compute the spectrum of Heff as a MPO itself.
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

  const string mpsfile=props.getProperty("mpsfile"); // full path expected
  // Hamiltonian parameters
  double Jx_=props.getDoubleProperty("Jx");
  double Jy_=props.getDoubleProperty("Jy");
  double Jz_=props.getDoubleProperty("Jz");
  double h_=props.getDoubleProperty("h");  
  int block = props.getIntProperty("block");
  if(block<=0) block=2;
  const string outfile=props.getProperty("outfile");  
  bool app=(props.getIntProperty("append")!=0); // default: append
  int nrEV=props.getIntProperty("numberEV"); // max nr of eigenvalues to record
  const string outfileS=props.getProperty("outfileS");  
  bool computeS=!outfileS.empty();
  const string outfileSz=props.getProperty("outfileSz");  
  bool computeSz=!outfileSz.empty();
  const string outfileH2=props.getProperty("outfileVariance");  
  bool computeH2=!outfileH2.empty(); // whether to compute the variance of the MPSs for the excitations  

  int pos1 = props.getIntProperty("initpos");
  int pos2 = props.getIntProperty("finalpos");
  int Dmps = props.getIntProperty("D");
  
  // Read the MPS
  if(!file_exists(mpsfile)){
    cout<<"ERROR: File for MPS "<<mpsfile<<" not found!"<<endl;
    exit(1);
  }

  MPS state(1,1,1);
  state.importMPS(mpsfile.data());
  cout<<"Read GS file "<<endl;
  int L=state.getLength();

  state.gaugeCond('L',1);



  
  if(pos1<=0) pos1=1;
  if(pos2<=0||pos2>=L-block-1) pos2=L-block-2;

  
  // Prepare output file
  ofstream* out;
  if(!app||!file_exists(outfile))
    out=new ofstream(outfile.data());
  else{
    out=new ofstream(outfile.data(),ios::app);
    // Write the header
  }
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfile<<
      " for output (append="<<app<<")"<<endl;
    exit(1);
  }
  out->close();delete out;

  if(computeS){
    if(!app||!file_exists(outfileS))
      out=new ofstream(outfileS.data());
    else
      out=new ofstream(outfileS.data(),ios::app);
    if(!out->is_open()){
      cout<<"Error: impossible to open file "<<outfileS<<
	" for output (append="<<app<<")"<<endl;
      exit(1);
    }
    out->close();delete out;
  }

  vector<double> Jx(L-1,Jx_);
  vector<double> Jy(L-1,Jy_);
  vector<double> Jz(L-1,Jz_);
  vector<double> h(L,h_);

  HeisenbergHamiltonian hamH(L,Jx,Jy,Jz,h,d);
  const MPO& mpoH=hamH.getHMPO();
  MPO Szmpo(L);
  SpinMPO::getSzMPO(L,d,Szmpo);
 
  Contractor& contractor=Contractor::theContractor();
  contractor.setConvTol(tol);
  contractor.setEigenSolver(primme);

  // Now it may be that my state was a purification: check
  bool purif=0; MPO doubleH(L);MPO doubleSz(L);
  int dphys=state.getA(0).getd();
  if(dphys==d*d){
    cout<<"The state read was a purification: extending the MPO of H"<<endl;
    purif=true;
    extendMPO(mpoH,doubleH,d);
    extendMPO(Szmpo,doubleSz,d);
  }

  
  complex_t initE;
  if(!purif)
    initE=contractor.contract(state,mpoH,state);
  else
    initE=contractor.contract(state,doubleH,state);
  cout<<"The energy of the read state is "<<initE<<endl;
  out=new ofstream(outfile.data(),ios::app);
  *out<<setprecision(15);
  *out<<"% E="<<real(initE)<<endl;
  *out<<endl;
  out->close();delete out;

  // max dim of the SVD, if I compute entropies
  int maxSVD=state.getA(L/2).getDr()*pow(d,(int)block/2);


  // First, impose canonical cond to right up to pos1 (to L was in all)
  int lastGaugeR=0;

  
  for(int k=pos1;k<=pos2;k++){
    // impose canonical rto R up to pos1
    for(int kg=lastGaugeR;kg<k;kg++)
      state.gaugeCond(kg,'R',true);
    lastGaugeR=k-1;
    
    MPOMultiplierHermitian Heff(block+2); // maybe could use non-Hermitian, as it is in fact Hermitian (factor 2 cheaper)
    MPOMultiplier Szeff(block+2); 
    if(!purif){
      contractor.getEffectiveOperatorMPOMultiplier(state,mpoH,block,k,Heff);
      contractor.getEffectiveOperatorMPOMultiplier(state,Szmpo,block,k,Szeff);
    }
    else{
      contractor.getEffectiveOperatorMPOMultiplier(state,doubleH,block,k,Heff);
      contractor.getEffectiveOperatorMPOMultiplier(state,doubleSz,block,k,Szeff);
    }
    
    if(purif){
      cout<<"******   ERROR: Behaviour for purification not yet clear with the MPO version of efectiveH !!!"<<endl;exit(1);
      // my Heff is a product of boundaries witht eh true H and
      // identities on ancillas. The identity is tensorized with the
      // rest and I want to remove it, as it only gives degeneracy
      int firstOpPos=1;
      int lastOpPos=block;
      for(int p=firstOpPos;p<=lastOpPos;p++){
	const FoldedOperator& op=(FoldedOperator&)Heff.getOp(p);
	//cout<<"Operator at p is "<<op<<endl;
	mwArray aux=op.getDataLeft();
	Heff.setOp(p,new Operator(aux),true);
	//cout<<"Substituted op "<<p<<" by simple"<<endl;
	const FoldedOperator& opZ=(FoldedOperator&)Szeff.getOp(p);
	//cout<<"Operator at p is "<<op<<endl;
	aux=opZ.getDataLeft();
	Szeff.setOp(p,new Operator(aux),true);
      }
    }

    // Copy the current GS as it is the right lowest eigenvector
    MPS gs(block+2,Dmps,d);
    for(int p=1;p<=block;p++){
      gs.replaceSite(p,state.getA(k+p-1).getA(),false);
    }
    int Dl=gs.getA(1).getDl();
    int Dr=gs.getA(block).getDr();
    gs.replaceSite(0,reshape(identityMatrix(Dl),Indices(Dl,1,Dl)),false);
    gs.replaceSite(block+1,reshape(identityMatrix(Dr),Indices(Dr,Dr,1)),false);
    
    vector<MPS*> levels; levels.push_back(&gs);
    vector<double> energies; energies.push_back(real(initE));
    vector<double> variances; variances.push_back(real(contractor.contract2(mpoH,state)-initE*conjugate(initE)));
    //vector<Level> energies; energies.push_back(Level(0,real(initE)));

    // Check the energy
    complex_t computedE0=contractor.contract(gs,Heff,gs);
    cout<<"Energy of the GS computed with Heff("<<k<<") "<<real(computedE0)<<endl;
    
    int lval=0; // last one computed
    while(lval<nrEV-1){
      double lambda=L;
      // compute one more
      MPS* exc=new MPS(gs);
      exc->setRandomState();
      exc->gaugeCond('R',1);
      exc->gaugeCond('L',1);

      //cout<<"In newly created MPS energy is "<<contractor.contract(*exc,Heff,*exc)<<endl;
      // orthogonalize it to gs
      //contractor.orthogonalize(*exc,Dmps,levels,vector<double>(levels.size(),1.));
      contractor.findNextExcitedState(Heff,Dmps,levels,&lambda,*exc,-L);
      levels.push_back(exc);
      energies.push_back(lambda);
      lval++;
      //cout<<"Found level "<<lval<<" with energy "<<lambda<<endl;
      if(computeH2){
	// need to reconstruct the full tensor, but will have to absorb tensors at 0 and block+1
	//exc->gaugeCond(0,'R',true); // this only changes the basis, but does not guarantee that they are back to identity
	//exc->gaugeCond(block+1,'L',true);
	MPS aux(state);
	for(int p=1;p<=block;p++){
	  //cout<<"Setting site "<<p<<" of small mps, in the tensor of site "<<k+p-1<<" of big one (for block="<<block<<", starting "<<k<<")"<<endl;
	  mwArray tens=exc->getA(p).getA(); // diff dims
	  if(p==1){
	    mwArray A0=exc->getA(0).getA(); // Dlx1xDl
	    Indices dims=tens.getDimensions();// d,Dl,Dp
	    tens.permute(Indices(2,1,3)); // Dl,d*Dp
	    tens.reshape(Indices(dims[1],dims[0]*dims[2]));
	    tens.multiplyLeft(reshape(A0,Indices(Dl,Dl)));
	    tens.reshape(Indices(Dl,dims[0],dims[2]));
	    tens.permute(Indices(2,1,3)); // d,Dl,Dp
	  }
	  if(p==block){
	    mwArray Ab=exc->getA(block+1).getA(); // DrxDrx1
	    Indices dims=tens.getDimensions();// d,Dp,Dr
	    tens.reshape(Indices(dims[0]*dims[1],dims[2]));
	    tens.multiplyRight(permute(reshape(Ab,Indices(Dr,Dr)),Indices(2,1)));
	    tens.reshape(Indices(dims[0],dims[1],Dr));
	  }
	  aux.replaceSite(k+p-1,tens);	  
	}
	// now check energy
	complex_t enerL=contractor.contract(aux,mpoH,aux);
	//cout<<"Checking: energy found "<<lambda<<", computed small: "<<contractor.contract(*exc,Heff,*exc)<<" full: "<<enerL<<", norm "<<contractor.contract(aux,aux)<<" and of exc "<<contractor.contract(*exc,*exc)<<endl;
	//cout<<"Overlap with gs: \nsmall :"<<contractor.contract(*exc,gs)<<"\nfull: "<<contractor.contract(aux,state)<<endl;
	complex_t valH2=contractor.contract2(mpoH,aux);
	variances.push_back(real(valH2-enerL*conjugate(enerL)));
	
      }
      
      if(computeS){
	vector<complex_t> sv_;
	contractor.getSchmidtValues(*exc,sv_);	
	out=new ofstream(outfileS.data(),ios::app);
	*out<<setprecision(15);
	*out<<k<<"\t"<<lval<<"\t";
	for(int ip=0;ip<sv_.size();ip++){
	  *out<<real(sv_[ip])<<"\t";
	}
	for(int ip=sv_.size();ip<maxSVD;ip++) // fill in
	  *out<<0.<<"\t";
	*out<<endl;
	out->close();delete out;
      }
      if(computeSz){
	//cout<<"Will now try to apply Szeff:"<<Szeff<<endl;
	complex_t valSz=contractor.contract(*exc,Szeff,*exc);	
	out=new ofstream(outfileSz.data(),ios::app);
	*out<<setprecision(15);
	*out<<k<<"\t"<<lval<<"\t"<<real(valSz)<<"\t"<<imag(valSz)<<endl;
	out->close();delete out;
      }
    }

    out=new ofstream(outfile.data(),ios::app);
    *out<<setprecision(15);
    *out<<k<<"\t";
    for(int p=0;p<min(nrEV,(int)energies.size());p++){
	*out<<energies[p]<<"\t";
    }
    if(energies.size()<nrEV) // fill in with 0
      for(int p=energies.size();p<nrEV;p++) *out<<0<<"\t";
    *out<<endl;
    out->close();delete out;

    if(computeH2){
      out=new ofstream(outfileH2.data(),ios::app);
      *out<<setprecision(15);
      *out<<k<<"\t";
      for(int p=0;p<min(nrEV,(int)energies.size());p++){
	*out<<variances[p]<<"\t";
      }
      if(variances.size()<nrEV) // fill in with 0
	for(int p=variances.size();p<nrEV;p++) *out<<0<<"\t";
      *out<<endl;
      out->close();delete out;
    }
  }
  
}
