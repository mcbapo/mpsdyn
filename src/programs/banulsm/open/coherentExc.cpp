
//#include "mwArray.h"
//#include "MPO.h"
#include "CoherentDissipationLiouvillian.h"
//#include "Operator.h"
#include "Contractor.h"
#include "Properties.h"
//#include "JoinedOperator.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <stdio.h>

#define MAXLEN 200  // Max nr of characters in a filename


/** 
    Program coherentExc: Compute the steady state of a certain
    Liouvillian
    as in \class <CoherentDissipationLiouvillian>,
    with g, gamma given as input parameters.
    In this case, I will compute (to a lower precision) several
    (excited) states

    To compile: make coherentEx

    Receives as argument a Properties file, such as
    config/excitProp.conf, from where the required parameters are
    read.  Additional parameters can be given to override the values
    in the Properties file.
    To replace some of the values in the file by a command line
    option, the additional arguments have to be of the form
    -propName=propValue

    The properties that will be required:
    + L (int) length of the chain
    + D (int) maximum bond dimension
    + g (double) parameter $g$ of the Hamiltonian
    + gamma (double) parameter $gamma$ of the dissipation
    + scale (double) scale factor to multiply Lmpo
    + tol (double) tolerance asked from Contractor
    + mpsfile (string) file where the resulting MPS will be stored.
    + oldmpsfile [OPTIONAL] (string) name of another (binary) file, containing 
                    an earlier MPS to be used as starting point
    + outdir (string) directory where to write the (automatically
                    named) file(s) containing observables.
    + nLev (int) number of states I will accumulate

*/

using namespace std;
using namespace shrt;

/** Given a state, compute all the expectation values (on each
    position) of a given single body operator. The results are
    written to the ostream (could be a file, or cout). */
void computeSingleBodyEV(const MPS& rho,const mwArray& op1,ostream& os);

/** Same as before, but compute the expectation value of the product
    of op1 acting on one site and op2 acting on the site plus l.  If
    maxLen!=0 is provided, it writes up to maxLen elements in the
    line, padding with 0 when needed.  
*/
void computeTwoBodyEV(const MPS& rho,const mwArray& op1,
		      const mwArray& op2,int l,ostream& os,
		      int maxLen=0);

/** Prepare header for a file */
void writeHeader(ofstream& out,const CoherentDissipationLiouvillian& L,
		 double lambda0,
		 int D,double scal,double costT,const string oldMPSfile);

/** Global matrices */
int d=2;
complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
mwArray sigX(Indices(d,d),dataX);
complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
mwArray sigY(Indices(d,d),dataY);
complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
mwArray sigZ(Indices(d,d),dataZ);

//double noise=1E-60;

int main(int argc,const char* argv[]){
  int cntr=0;
  const char* infile=argv[++cntr];
  string directory="";
  // Recover the necessary parameters from a Properties file
  Properties props(infile);
  if(argc>2){
    cntr++;
    cout<<"Some properties may be now replaced by command line arguments"
	<<endl;
    props.loadProperties(argc-cntr,&argv[cntr]);
  }


  int L=props.getIntProperty("L");
  int D=props.getIntProperty("D");
  double g=props.getDoubleProperty("g");
  double gamma=props.getDoubleProperty("gamma");
  double scal=props.getDoubleProperty("scale");
  double tol=props.getDoubleProperty("tol");
  const string newMPSfile=props.getProperty("mpsfile");
  const string oldMPSfile=props.getProperty("oldmpsfile");
  bool usingOldMPS=(!oldMPSfile.empty()&&oldMPSfile.length()>1);
  const string outdir=props.getProperty("outdir");
  int nLev=props.getIntProperty("nLev");

  cout<<"oldMPSfile is "<<oldMPSfile<<endl;
  cout<<"its size is "<<oldMPSfile.size()<<endl;
  cout<<"its length is "<<oldMPSfile.length()<<endl;


  // I construct here the name of the text output files, where I will
  // save the results
  char fileNameC[MAXLEN];
  // For the single body observables (x,y,z, site after site, each in
  // one line)
  sprintf(fileNameC,"coherentExcSingle_%d_%1.2f_%D.txt",L,gamma,D);
  const string out1bodyfile=outdir+"/"+fileNameC;
  // For the two body observables (xx,yy,zz, site after site, first
  // for distance l=1, until N/2, each (xx, yy, zz) for each l
  // in one line, whichstarts with an identifier and the value of l
  sprintf(fileNameC,"coherentExcDouble_%d_%1.2f_%D.txt",L,gamma,D);
  const string out2bodyfile=outdir+"/"+fileNameC;


  // Open here the output file: for normal output Repeat the same for
  // as many output files as required, but not necessarily here at the
  // beginning!
  ofstream* out1=new ofstream(out1bodyfile.data(),ios::app);
  ofstream* out2=new ofstream(out2bodyfile.data(),ios::app);
  if(!out1->is_open()){
    cout<<"Error: impossible to open file "<<out1bodyfile<<
      " for output"<<endl;
    exit(1);
  }
  if(!out2->is_open()){
    cout<<"Error: impossible to open file "<<out2bodyfile<<
      " for output"<<endl;
    exit(1);
  }
  cout<<"Writing output to files: \n\t"<<out1bodyfile
      <<"\n\t"<<out2bodyfile<<endl;

  //srandom(time(NULL));
  // Construct the Liouvillian
  CoherentDissipationLiouvillian model(L,scal*g,scal*gamma);
  
  // Liouvillian
  const MPO& Lop=model.getLMPO();  
  // Adjoint Liouvillian
  const MPO& Ldagger=model.getLAdjointMPO();  

  // // To debug
  // if(L<=5)
  //   Lop.exportForMatlab("myLmpo.m");


  // And now construct a Joined MPO for the L+ L
  const MPO* ptrs[2]={&Ldagger,&Lop};
  MPO LdaggerL(L);
  MPO::join(2,ptrs,LdaggerL);

  // And directly try for the GS!
  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  if(D<=12)
    contractor.setEigenSolver(fulleig);
  contractor.setConvTol(tol);

  MPS myGS(L,D,d*d);
  // Initialize the MPS with something already existing, maybe?
  if(usingOldMPS){
    myGS.importMPS(oldMPSfile.data());
    //    myGS.increaseBondDimensionWithNoise(D,noise);
    if(D>=myGS.getBond())
      myGS.increaseBondDimension(D);
    else
      myGS=MPS(myGS,D);
    cout<<"Using initial MPS from file "<<oldMPSfile<<endl;
  }
  else{
    myGS.setRandomState();
    cout<<"Starting from random MPS"<<endl;
  }

  myGS.gaugeCond('R',1);
  myGS.gaugeCond('L',1);

  double lambda=real(contractor.contract(myGS,LdaggerL,myGS));
  cout<<"Starting from an energy value "<<lambda<<endl;
  double offset=0.; // Contractor will set it alone
  //if(abs(lambda)<1E-5)
  //offset=10*lambda;

  // // Truquillo 1)
  // contractor.sweepPart(L-1,(L-1)/2+1,LdaggerL,D,myGS,offset);
  // lambda=real(contractor.contract(myGS,LdaggerL,myGS));
  // cout<<"After half sweep energy value "<<lambda<<endl;

  // // Truquillo 2)
  // {MPS aux(myGS);
  //   contractor.sweepPart(0,(L-1)/2+1,LdaggerL,D,aux,offset);
  //   // And change! To force reflection symmetry at the beginning
  //   for(int pos=(L-1)/2+1;pos<L-1;pos++){
  //     // get the tensor, flip it and set it to reflected place
  //     myGS.setRotatedA(pos,aux.getA(L-1-pos).getA(),Indices(1,3,2));
  //   }
  // }

  // Truquillo 3)
  for(int posL=(L+1)/2;posL>=0;posL--){
    contractor.sweepPart(posL,L-1-posL,LdaggerL,D,myGS,offset);
    contractor.sweepPart(L-1-posL,posL,LdaggerL,D,myGS,offset);
    cout<<"After sweep "<<posL<<"-"<<L-1-posL<<" energy value "<<contractor.contract(myGS,LdaggerL,myGS)<<endl;
  }

  // Truquillo 4) (more costly, but rigorous -?- )
  for(int pp=0;pp<3;pp++){
    for(int pos=0;pos<(L-1)/2+1;pos++){
      contractor.sweepPart(pos,pos+1,LdaggerL,D,myGS,offset); // one right
      contractor.sweepPart(L-1-pos,L-1-pos-1,LdaggerL,D,myGS,offset); // one left
    }
    lambda=real(contractor.contract(myGS,LdaggerL,myGS));
    cout<<"After symmetric sweep nr "<<pp+1<<"energy value "<<lambda<<endl;
  }

  clock_t start=clock();
  contractor.findGroundState(LdaggerL,D,&lambda,myGS,offset);
  clock_t end=clock();

  cout<<"Finished GS search with lambda="<<lambda<<endl;

  // Now I have the GS, I will write it to a file
  myGS.exportMPS(newMPSfile.data());


  // I write out some results to output file, and close it (instead, I
  // could start computing correaltions afterwards) Actually, I could
  // open out just now, because it does not need to be open all during
  // the calculation of Contractor
  double cost=(double)(end-start)/CLOCKS_PER_SEC;
  
  // Now compute observables and write results to files
  writeHeader(*out1,model,lambda,D,scal,cost,oldMPSfile);
  *out1<<setprecision(10);
  writeHeader(*out2,model,lambda,D,scal,cost,oldMPSfile);
  *out2<<setprecision(10);
  *out2<<"% Every line starts with an integer: 0 for the lowest eigenvalue "
       <<"(would be one for the next one). Then the expectation value of L+L in "
       <<"this state. It follows another integer: 1 for XX, "
       <<"2 for YY, 3 for ZZ correlators; the next value is the distance l "
       <<"(so correlators  <X(i)X(i+3)> will be in a line starting with 0 1"
       <<" 3 and then a series of double values from i=1 to the end"<<endl;

  vector<MPS*> levels;
  levels.push_back(&myGS);
  double offsetE=-1;
  bool stopping=false;
  for(int lev=0;lev<nLev&&!stopping;lev++){
    if(lev>0){
      MPS* exc=new MPS(L,D,d*d);
      exc->setRandomState();
      double lambda1=0.;
      contractor.findNextExcitedState(LdaggerL,D,levels,&lambda1,*exc,scal*scal*offsetE);
      if(lambda1+scal*scal*offsetE>0){
	cout<<"Stopping iteration as I got an energy above zero!"<<endl;
	stopping=true;
      }
      levels.push_back(exc);
    }
    const MPS& exc=*levels[lev];

    complex_t LLval=contractor.contract(exc,LdaggerL,exc);

    // Compute now expectation values of sigmaX,Y,Z
    *out1<<lev<<"\t"<<LLval<<"\t"<<1<<"\t";
    computeSingleBodyEV(exc,sigX,*out1);
    *out1<<endl;
    *out1<<lev<<"\t"<<LLval<<"\t"<<2<<"\t";
    computeSingleBodyEV(exc,sigY,*out1);
    *out1<<endl;
    *out1<<lev<<"\t"<<LLval<<"\t"<<3<<"\t";
    computeSingleBodyEV(exc,sigZ,*out1);
    *out1<<endl;

    // Compute now expectation values of sigmaXX,YY,ZZ at distance 1-3
    int maxLen=L-1; // the largest number of terms
    for(int l=1;l<=L/2;l++){
      *out2<<lev<<"\t"<<LLval<<"\t"<<1<<"\t"<<l<<"\t";
      computeTwoBodyEV(exc,sigX,sigX,l,*out2,maxLen);
      *out2<<endl;
      *out2<<lev<<"\t"<<LLval<<"\t"<<2<<"\t"<<l<<"\t";
      computeTwoBodyEV(exc,sigY,sigY,l,*out2,maxLen);
      *out2<<endl;
      *out2<<lev<<"\t"<<LLval<<"\t"<<3<<"\t"<<l<<"\t";
      computeTwoBodyEV(exc,sigZ,sigZ,l,*out2,maxLen);
      *out2<<endl;
    }
    
  }



  out1->close();
  delete out1;
  out2->close();
  delete(out2);

}

void computeSingleBodyEV(const MPS& rho,const mwArray& op1,ostream& os){
  int L=rho.getLength();
  MPS aux(L,1,d*d);
  mwArray basic=identityMatrix(d); // basic piece of the MPO
  basic.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++){
    aux.setA(k,basic);
  }

  Contractor& contractor=Contractor::theContractor();
  double traceRho=real(contractor.contract(rho,aux));

  // WARNING: No check of dimensions done!!
  mwArray op1_=reshape(op1,Indices(d*d,1,1));
  op1_.conjugate(); // It gets conjugated in the scalar product! 

  // Now, site by site, replace the identity by the operator,
  // contract, write the result to os, and replace the operator by the
  // identity again
  for(int k=0;k<L;k++){
    aux.setA(k,op1_); // operator on site k
    complex_t valK=contractor.contract(rho,aux);
    os<<real(valK)/traceRho<<"\t";
    aux.setA(k,basic); // operator removed
  }
}

void computeTwoBodyEV(const MPS& rho,const mwArray& op1,
		      const mwArray& op2,int l,ostream& os,
		      int maxLen){
  int L=rho.getLength();
  MPS aux(L,1,d*d);
  mwArray basic=identityMatrix(d); // basic piece of the MPO
  basic.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++){
    aux.setA(k,basic);
  }

  Contractor& contractor=Contractor::theContractor();
  double traceRho=real(contractor.contract(rho,aux));

  // WARNING: No check of dimensions done!!
  mwArray op1_=reshape(op1,Indices(d*d,1,1));
  op1_.conjugate(); // It gets conjugated in the scalar product! 
  mwArray op2_=reshape(op2,Indices(d*d,1,1));
  op2_.conjugate(); // It gets conjugated in the scalar product! 

  // Now, site by site, replace the identity at k by the operator op1,
  // the identity in k+l by op2,
  // contract, write the result to os, and replace the operators by the
  // identity again
  int cnt=0; // nr of terms written
  for(int k=0;k+l<L;k++){
    aux.setA(k,op1_); // operator on site k
    aux.setA(k+l,op2_); // operator on site k+l
    complex_t valK=contractor.contract(rho,aux);
    os<<real(valK)/traceRho<<"\t";
    cnt++;
    aux.setA(k,basic); // operator1 removed
    aux.setA(k+l,basic); // operator2 removed
  }
  if(maxLen>0)
    while(cnt<maxLen){
      os<<0.<<"\t";
      cnt++;
    }

}

void writeHeader(ofstream& out,const CoherentDissipationLiouvillian& L,
		 double lambda0,
		 int D,double scal,double costT,const string oldMPSfile){
  out<<"% L="<<L.getLength()<<endl;
  out<<"% g="<<L.getG()<<endl;
  out<<"% gamma="<<L.getGamma()<<endl;
  out<<"% D="<<D<<endl;
  out<<"% scale="<<scal<<endl;
  Contractor& contr=Contractor::theContractor();
  out<<"% Contractor tol="<<contr.getConvTol()<<endl;
  out<<"% Contractor solver="<<contr.getSolverName()<<endl;
  out<<"% Time (s) needed = "<<costT<<endl;
  if(!oldMPSfile.empty())
    out<<"% Using initial state from "<<oldMPSfile<<endl;
  else
    out<<"% Using random initial state "<<endl;
  out<<"% Lowest eigenvalue found to be "<<lambda0<<endl;
}

  // MPO H(L),LI(L),IL(L),Ha(L);
  // Nt=L;
  // MPO Ham(Nt);
  // MPO H6(Nt);
  // Operator* op1;
  //   for(int j=0;j<Nt;j++){
  //     cout<<"Doing double operator...at "<<j<<endl;
  // 	//op1=new Operator(L6.getOp(j).doubleop(0));    Do=L6.getOp(0).getDr();//bond dimension of MPO;
  // 	op1=new Operator(L5.getOp(j).doubleop(0));  Do=L5.getOp(0).getDr();
  // 	d=4;
  //     Ham.setOp(j,op1);
  //     double LM=1;
  //     LM=pow(LM,(double)1/L);
  //     if(1){// directly add an LM*Identity to the MPO L^\dag L;if this is taken then the findgroundstate should use Ha or H5_2
  // 	//    *out<<Ham.getOpData(j)<<endl;
  // 	if(j==0){
  //         mwArray* temp0=new mwArray(Indices(d,1,d,Do*Do+1));
  // 	  temp0->fillWithZero();
  // 	  for(int ju=0;ju<d;ju++){ 
  // 	    for(int jd=0;jd<d;jd++){ 
  // 	      for(int jl=0;jl<1;jl++){
  // 		for(int jr=0;jr<Do*Do;jr++){
  // 		  temp0->setElement(Ham.getOpData(j).getElement(Indices(ju,jl,jd,jr)),Indices(ju,jl,jd,jr));  
  // 		}
  //             }
  // 	      double ele;
  // 	      if(ju==jd) ele=LM;
  // 	      else ele=0.0;
  // 	      temp0->setElement(ele,0.0,Indices(ju,0,jd,Do*Do));
	      
  // 	    }
  //         }
  //         Operator* op0=new Operator(*temp0);  
  // 	  //H6.setOp(j,op0);  
  // 	  Ha.setOp(j,op0);

  // 	}     
  // 	else if(j==Nt-1)
  // 	  {
  // 	    mwArray* temp0=new mwArray(Indices(d,Do*Do+1,d,1));
  // 	    temp0->fillWithZero();
  // 	    for(int ju=0;ju<d;ju++){ 
  // 	      for(int jd=0;jd<d;jd++){ 
  // 		for(int jr=0;jr<1;jr++){
  // 		  for(int jl=0;jl<Do*Do;jl++){
  // 		    temp0->setElement(Ham.getOpData(j).getElement(Indices(ju,jl,jd,jr)),Indices(ju,jl,jd,jr));
  // 		  } 
  // 		}
  // 		double ele;
  // 		if(ju==jd) ele=LM;
  // 		else ele=0.0;
  // 		temp0->setElement(ele,0.0,Indices(ju,Do*Do,jd,0));
  // 	      }
  // 	    }
  // 	    Operator* op0=new Operator(*temp0);  
  // 	      Ha.setOp(j,op0);  
  // 	    //H6.setOp(j,op0);  
  // 	  }      
  // 	else {
  //         mwArray* temp0=new mwArray(Indices(d,Do*Do+1,d,Do*Do+1));
  // 	  temp0->fillWithZero();
  //         for(int ju=0;ju<d;ju++){ 
  // 	    for(int jd=0;jd<d;jd++){           
  // 	      for(int jr=0;jr<Do*Do;jr++){
  // 		for(int jl=0;jl<Do*Do;jl++){
  // 		  temp0->setElement(Ham.getOpData(j).getElement(Indices(ju,jl,jd,jr)),Indices(ju,jl,jd,jr));
  //               }
  // 	      } 
  // 	      double ele;
  // 	      if(ju==jd) ele=LM;
  // 	      else ele=0.0;
  // 	      temp0->setElement(ele,0.0,Indices(ju,Do*Do,jd,Do*Do));
  // 	    }
  // 	  }
  // 	  Operator* op0=new Operator(*temp0);  
  // 	    Ha.setOp(j,op0);
  // 	    //H6.setOp(j,op0);
  // 	} 
  //     }//end adding Identity
  //   }//end loop of each site
    
  //   //check the spectrum of MPO
  //   if(0){//Ham(for both split or not) / Ha(not split) or H5_2(split) 
  //     mwArray* res=new mwArray(Ham.getOpData(0));
  //     cout<<"the MPO L^dag L to check is "<<Ham<<endl;
  //     for(int i=1;i<Ham.getLength();i++){
  // 	mwArray temp=Ham.getOpData(i);
  // 	int utemp=temp.getDimension(0);
  // 	int ltemp=temp.getDimension(1);
  // 	int dtemp=temp.getDimension(2);
  // 	int rtemp=temp.getDimension(3);
  // 	int ures=(*res).getDimension(0);
  // 	int lres=(*res).getDimension(1);
  // 	int dres=(*res).getDimension(2);
  // 	int rres=(*res).getDimension(3);
  // 	(*res).reshape(Indices(ures*lres*dres,rres));
  // 	temp.permute(Indices(2,1,3,4));
  // 	temp.reshape(Indices(ltemp,utemp*dtemp*rtemp));
  // 	(*res)=(*res)*temp;
  // 	(*res).reshape(Indices(ures,lres,dres,utemp,dtemp,rtemp));
  // 	(*res).permute(Indices(1,4,2,3,5,6));
  // 	(*res).reshape(Indices(ures*utemp,lres,dres*dtemp,rtemp));
  //     }
  //     (*res).reshape(Indices((*res).getDimension(0),(*res).getDimension(2)));
  //      char filename[]="/afs/ipp/home/j/jcui/mpsdyn/src/bin/mpotocheck";
  //      ofstream outfile(filename);
  //      (*res).save(outfile);
  //     cout<<setprecision(16);
  //     cout<<"Check"<<(*res)<<endl;
  //     vector<complex_t> Dval;
  //     mwArray U;
  //     delete res;
  //   }
  //   // end checking the spectrum of MPO
    
  //   //begin calculation of finding ground state and/or saving ground state MPS
  //   MPS gsP;
  //   MPS gs(Nt,Dim,d);
  //   //      int D=10;
  //   double lambda=0.0;
  //   double gvalue=0.0;
  //   MyContractor& contractor=MyContractor::theContractor();
  //   Contractor& oldcontractor=Contractor::theContractor();
  //   double dd=d*d;
  //   if(1){      //begin calculation of finding ground state and/or saving ground state MPS
  //     int no=0;
  //     int z14=0;//number of Znaupd -14 error
  //     int incE=0;//number of energy increasing
  // 	//set random ground state mps or import mps
  // 	if (call==1){
	  
  // 	  //MPS gs;	
  // 	  char address[150];
  // 	  strcpy(address,"/afs/ipp/home/j/jcui/mpsdyn/src/bin/IO/IOin/IO_GSmps_gamma=");
  // 	  char gam2is[5];
  // 	  sprintf(gam2is,"%1.2f",gamma);
  // 	  strcat(address,gam2is);
  // 	  //strcat(address,"_r2=");
  // 	  //char fifth[5];
  // 	  //sprintf(fifth,"%1.2f",gamma2);
  // 	  //strcat(address,fifth);    
  // 	  strcat(address,"_N=");
  // 	  char third[3];
  // 	  sprintf(third,"%d",L);
  // 	  strcat(address,third);
  // 	  strcat(address,"_D=");
  // 	  char da[3];
  // 	  sprintf(da,"%d",Dim);
  // 	  strcat(address,da);
  // 	  char ends[]=".txt";
  // 	  strcat(address,ends);
  // 	  gs.importMPStext(address);
  // 	}
  // 	else{
  // 	  //MPS gs(Nt,Dd,d);
  // 	  gs.setRandomState();
  // 	  cout<<"starting findGSeig..."<<endl;
  // 	   round=contractor.findGSIOeig(Ha,Dim,&lambda,gs,energy_,correlationx,correlationy,correlationz);//eig,test converge 
  // 	}
	
  // 	*out<<setprecision(16);
  // 	*out<<"Initial energy is "<<contractor.contract(gs,Ha,gs)/contractor.contract(gs,gs)<<endl;	
  // 	  // round=contractor.findGroundStateneweigs(Ha,Dd,&lambda,gs,z14,incE);//new eigs(Heff+eigs)
  // 	  //round=contractor.findGSIOeigs(Ha,Dd,&lambda,gs,z14,incE,energy_,correlationx,correlationy,correlationz);//eigs,test converge(eigs,if z-14 use cg or cgeigs or eig and save Heff matrix) 
  // 	  //round=contractor.findGSIOeigs0(Ha,Dd,&lambda,gs,z14,incE,correlationzz,correlationx,correlationy,correlationz);//eigs,test converge(eigs,if z-14 increase reserved states in eigs) 
  // 	  //round=contractor.findGSIOeigs3(Ha,Dd,&lambda,gs,z14,incE,correlationzz,correlationx,correlationy,correlationz);//eigs,test converge(eigs,if z-14 redo findgs with Heff) 
  // 	  //	  round=contractor.findGSIOeigs2(Ha,Dd,&lambda,gs,z14,incE,correlationzz,correlationx,correlationy,correlationz);//eigs,test converge,if z14 happens skip that site, and/or save the effective Hamiltonian
  // 	  //contractor.findgstestiteration(Ha,Dd,&lambda,gs);
  // 	  //contractor.findGSaddLM(H6,Dd,&lambda,gs,z14,incE,energy_,correlationx,correlationy,correlationz,correlationxx,correlationyy,correlationzz);//eigs with LM added to Heff, test converge(eigs,if z-14 use eig)
  // 	  //round=contractor.findGSIOcg(H6,Dd,&lambda,gs,energy_,correlationx,correlationy,correlationz);//cg,test converge 
  // 	  //round=contractor.findGSIOcg2(Ha,Dd,&lambda,gs,energy_,correlationx,correlationy,correlationz);//cgeigs,test converge 
  // 	gs.increaseBondDimension(Dd);
  //      char filename[]="/afs/ipp/home/j/jcui/mpsdyn/src/bin/mpotocheck";	
  //      char Energy[]="/afs/ipp/home/j/jcui/mpsdyn/src/bin/IO/IO_energy";

  //      cout<<"starting findGScg..."<<endl;
  //      round=contractor.findGSIOcg2(Ha,Dd,&lambda,gs,energy_,correlationx,correlationy,correlationz,correlationxx,correlationyy,correlationzz);//cgeigs,test converge, save the correlation profile in each iteration or every 10 iterations
  //      //round=contractor.findGroundState(Ha,Dd,&lambda,gs); //eig
	

  //      {gsP=gs;D=Dd;
  // 	  cout<<"ground energy is "<<lambda<<", mps bound dimension D is  "<<Dd<<" and Iteration times are round="<<round<<". Znaupd -14 happened "<<z14<<" times."<<endl;
  // 	  *out<<"ground energy is "<<lambda<<", mps bound dimension D is  "<<Dd<<" and Iteration times are round="<<round<<". Znaupd -14 happened "<<z14<<" times."<<endl;
  // 	}
      //finish finding the ground state with proper D
       //      delete  op1;
      // cout<<endl<<"Ground state found with eigenvalue "<<lambda<<endl<<endl;
      // cout<<"system size is "<<L<<endl;  
      // cout<<"D is "<<D<<", and g is "<<g<<". Bond of ground state MPS is "<<gsP.getBond()<<endl;
      // cout<<contractor.contract(gsP,gsP)<<endl;
      // *out<<"system size L is "<<L<<endl;  
      // *out<<"D is "<<D<<", gamma="<<gamma1<<", and Iteration times are round= "<<round<<endl;
      // //      *out<<"Znaupd -14 happens "<<z14<<"times."<<endl;
      // *out<<"Ground state found with eigenvalue "<<lambda<<endl;
  //    }    //end calculation of finding ground state and/or saving ground state MPS
    
    
    // if(1){  
    //     // begin correlations 
    //   //   if(0){
    // 	//   mwArray sig0=identityMatrix(2);
    // 	complex_t dataI[]={ONE_c,ZERO_c,ZERO_c,ONE_c};
    // 	complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    // 	complex_t datay[]={ZERO_c,I_c,-1*I_c,ZERO_c};
    // 	complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-1*ONE_c};
	
    // 	mwArray sig0(Indices(4,1,1),dataI);
    // 	mwArray sig1(Indices(4,1,1),datax);
    // 	mwArray sig2(Indices(4,1,1),datay); 
    // 	mwArray sig3(Indices(4,1,1),dataz);
	
    // 	MPS ID(L,1,4);
    // 	for(int k=0;k<L;k++){
    // 	  ID.setA(k,sig0);
    // 	}
    // 	MPS x(L,1,4),x1(L,1,4),x2(L,1,4),y(L,1,4),y1(L,1,4), y2(L,1,4), z(L,1,4),z1(L,1,4),z2(L,1,4), xx(L,1,4), yy(L,1,4), zz(L,1,4);
    // 	complex_t norm=contractor.contract(ID,gsP);
    // 	cout<<"normalization factor of state is "<<norm<<endl;
    // 	*out<<"normalization factor of state is "<<norm<<endl;
    // 	//*  single site correlation
    // 	for(int i=0;i<L;i++){
	  
    // 	  x=ID;
    // 	  x.setA(i,sig1);
    // 	  y=ID;
    // 	  y.setA(i,sig2);
    // 	  z=ID;
    // 	  z.setA(i,sig3);
	  
    // 	  *outcox<<contractor.contract(x,gsP)/norm<<endl;
    // 	  *outcoy<<contractor.contract(y,gsP)/norm<<endl;
    // 	  *outcoz<<contractor.contract(z,gsP)/norm<<endl;  
    // 	}
	
    // 	//*two-site correlation
    //     LB1=L/2;
    // 	for(int i=LB1;i<L;i++){
	 	  
    // 	  xx=ID;
    // 	  x1=ID;
    // 	  x2=ID;
    // 	  xx.setA(LB1-1,sig1);
    // 	  xx.setA(i,sig1);
    // 	  x1.setA(LB1-1,sig1);
    // 	  x2.setA(i,sig1);
	  
    // 	  yy=ID;
    // 	  y1=ID;
    // 	  y2=ID;
    // 	  yy.setA(LB1-1,sig2);
    // 	  yy.setA(i,sig2);
    // 	  y1.setA(LB1-1,sig2);
    // 	  y2.setA(i,sig2);
	  
    // 	  zz=ID;
    // 	  z1=ID;
    // 	  z2=ID;
    // 	  zz.setA(LB1-1,sig3);
    // 	  zz.setA(i,sig3);
    // 	  z1.setA(LB1-1,sig3);
    // 	  z2.setA(i,sig3);
	  
    // 	   *outxx<<contractor.contract(xx,gsP)/norm-contractor.contract(x1,gsP)*contractor.contract(x2,gsP)/norm/norm<<endl;
    // 	   *outyy<<contractor.contract(yy,gsP)/norm-contractor.contract(y1,gsP)*contractor.contract(y2,gsP)/norm/norm<<endl;
    // 	   *outzz<<contractor.contract(zz,gsP)/norm-contractor.contract(z1,gsP)*contractor.contract(z2,gsP)/norm/norm<<endl;
    // 	}
    // 	//end correlations
    // }

    // //save the ground state MPS    
    //   char address[150];
    //   strcpy(address,"/afs/ipp/home/j/jcui/mpsdyn/src/bin/IO/IOout/IO_GSmps_gamma=");
    //   char gam2is[5];
    //   sprintf(gam2is,"%1.2f",gamma);
    //   strcat(address,gam2is);
    //   //strcat(address,"_r2=");
    //   //char fifth[5];
    //   //sprintf(fifth,"%1.2f",gamma2);
    //   //strcat(address,fifth);
    //   strcat(address,"_N=");
    //   char third[3];
    //   sprintf(third,"%d",L);
    //   strcat(address,third);
    //   strcat(address,"_D=");
    //   char da[3];
    //   sprintf(da,"%d",D);
    //   strcat(address,da);
    //   char ends[]=".txt";
    //   strcat(address,ends);
    //   gs.exportMPStext(address);
    
    
    /*
    //convert final state vector to density matrix
    // int Dop=d*d;
    int b[L];
    int c[L];
    int e[2];
    mwArray rho(Indices(pow(2,L),pow(2,L)));
    rho.fillWithZero();
    for (int j=0;j<pow(2*d,L);j++){
    quaternary(j,L,b);
    binary(b[0],d,e);
    int k1=e[0];
    int k2=e[1];
    mwArray res=gsP.getA(0).getA();
    // res=gsP.getA(0).getA();
    //	cout<<res<<endl;
    //	cout<<b[0]<<endl;
    mwArray m(Indices(res.getDimension(1),res.getDimension(2)));
    for (int t1=0;t1<res.getDimension(1);t1++){
    for (int t2=0;t2<res.getDimension(2);t2++){ 
    m.setElement(res.getElement(Indices(b[0],t1,t2)),Indices(t1,t2));
    
    }
    }
    for (int p=1;p<L;p++){
    binary(b[p],d,e);
    int j1=e[0];
    int j2=e[1];
    k1=k1*2+j1;
    k2=k2*2+j2;
    mwArray res=gsP.getA(p).getA();
    //res=gsP.getA(p).getA();
    mwArray te(Indices(res.getDimension(1),res.getDimension(2)));
    for (int t1=0;t1<res.getDimension(1);t1++){
    for (int t2=0;t2<res.getDimension(2);t2++){ 
    te.setElement(res.getElement(Indices(b[p],t1,t2)),Indices(t1,t2));
    }
    }
    m=m*te;
    }
    rho.setElement(m.getElement(Indices(0,0)),Indices(k1,k2));
    //    rho.setElement(ONE_c,Indices(k1,k2));
    }
    cout<<"The final density matrix is "<<rho<<endl;
     
    */

