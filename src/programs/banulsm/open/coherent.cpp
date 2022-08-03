
//#include "mwArray.h"
#include "MPO.h"
#include "CoherentDissipationLiouvillian.h"
#include "Operator.h"
#include "Contractor.h"
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

//#define LOCALDIR "/afs/ipp/home/j/jcui/mpsdyn/src/bin/"
//#include <fstream>

/** 
    Program coherent: Compute the steady state of a certain
    Liouvillian
    as in \class <CoherentDissipationLiouvillian>,
    with g, gamma given as input parameters.

    To compile: make coherent

    Receives arguments:
    \param<outfname> (char*) name of the output file
    \param<L> (int) length of the chain
    \param<D> (int) maximum bond dimension
    \param<g> (double) parameter $g$ of the Hamiltonian
    \param<gamma> (double) parameter $gamma$ of the dissipation
    \param<scal> (double) scale factor to multiply Lmpo
    \param<tol> (double) tolerance asked from Contractor
    \param<dir> (char*) directory where the output is to be written
    \param<oldfile> [OPTIONAL] (char*) name of another (binary) file, containing 
                    an earlier MPS to be used as starting point

*/

using namespace std;
using namespace shrt;

/** Given a state, compute all the expectation values (on each
    position) of a given single body operator. The results are
    written to the ostream (could be a file, or cout). */
void computeSingleBodyEV(const MPS& rho,const mwArray& op1,ostream& os);

/** Same as before, but compute the expectation value of the product
    of op1 acting on one site and op2 acting on the site plus l. */
void computeTwoBodyEV(const MPS& rho,const mwArray& op1,
		      const mwArray& op2,int l,ostream& os);

/** Global matrices */
int d=2;
complex_t dataX[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
mwArray sigX(Indices(d,d),dataX);
complex_t dataY[]={ZERO_c,I_c,-I_c,ZERO_c};
mwArray sigY(Indices(d,d),dataY);
complex_t dataZ[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
mwArray sigZ(Indices(d,d),dataZ);


int main(int argc,const char* argv[]){
  int count=0;
  int L=atoi(argv[++count]);
  int D=atoi(argv[++count]);
  double g=atof(argv[++count]); 
  double gamma=atof(argv[++count]); 
  double scal=atof(argv[++count]); 
  double tol=atof(argv[++count]); 
  const char* basedir=argv[++count];  
  bool usingOldMPS=false;
  const char* oldMPSfile="";
  if(argc>8){
    usingOldMPS=true;
    oldMPSfile=argv[++count];  
  }


  // I construct here the name of the MPS output file, where I will save the tensors
  char newMPSfile[MAXLEN];
  sprintf(newMPSfile,"%s/MPS_%d_%1.2f_%D",basedir,L,gamma,D);

  // Open here the out put file: for normal output
  // Repeat the same for as many output files as required, but not necessarily here at the beginning! 
  char outfname[200];
  // Constructing the filenames dynamically!
  sprintf(outfname,"%s/coherent_%d_%1.2f_%D.txt",basedir,L,gamma,D);
  ofstream* out=new ofstream(outfname,ios::app);
  if(!out->is_open()){
    cout<<"Error: impossible to open file "<<outfname<<
      " for output"<<endl;
    exit(1);
  }
  cout<<"Writing output to file "<<outfname<<endl;

  

  //srandom(time(NULL));
  CoherentDissipationLiouvillian model(L,scal*g,scal*gamma);
  
  // Liouvillian
  const MPO& Lop=model.getLMPO();  
  // Adjoint Liouvillian
  const MPO& Ldagger=model.getLAdjointMPO();  

  // To debug
  if(L<=5)
    Lop.exportForMatlab("myLmpo.m");


  // And now construct a Joined MPO for the L+ L
  const MPO* ptrs[2]={&Ldagger,&Lop};
  MPO LdaggerL(L);
  MPO::join(2,ptrs,LdaggerL);

  // And directly try for the GS!
  Contractor& contractor=Contractor::theContractor();
  contractor.setEigenSolver(primme);
  if(D<=10)
    contractor.setEigenSolver(fulleig);
  contractor.setConvTol(tol);

  MPS myGS(L,D,d*d);
  // Initialize the MPS with something already existing, maybe?
  if(usingOldMPS){
    myGS.importMPS(oldMPSfile);
    myGS.increaseBondDimension(D);
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

  clock_t start=clock();
  contractor.findGroundState(LdaggerL,D,&lambda,myGS,offset);
  clock_t end=clock();

  cout<<"Finished GS search with lambda="<<lambda<<endl;

  // Now I have the GS, I will write it to a file
  myGS.exportMPS(newMPSfile);


  // I write out some results to output file, and close it (instead, I
  // could start computing correaltions afterwards) Actually, I could
  // open out just now, because it does not need to be open all during
  // the calculation of Contractor
  double cost=(double)(end-start)/CLOCKS_PER_SEC;
  
  *out<<"% Summary of the L+ L calculation for "<<endl;
  *out<<"% L="<<L<<endl;
  *out<<"% gamma="<<gamma<<endl;
  *out<<"% D="<<D<<endl;
  if(usingOldMPS)
    *out<<"% Using initial state from "<<oldMPSfile<<endl;
  else
    *out<<"% Using random initial state "<<endl;
  *out<<"% And scale factor "<<scal<<endl; 
  *out<<"% Lowest eigenvalue found to be "<<lambda<<endl;
  *out<<"% Using convergence criterion "<<tol<<endl;
  *out<<"% time cost is "<<cost<<" seconds"<<endl;

  // Compute now expectation values of sigmaX,Y,Z
  *out<<"% X \t";
  computeSingleBodyEV(myGS,sigX,*out);
  *out<<endl;
  *out<<"% Y \t";
  computeSingleBodyEV(myGS,sigY,*out);
  *out<<endl;
  *out<<"% Z \t";
  computeSingleBodyEV(myGS,sigZ,*out);
  *out<<endl;

  // Compute now expectation values of sigmaXX,YY,ZZ at distance 1-3
  for(int l=1;l<=3;l++){
    *out<<"% X(k)X(k+"<<l<<") \t";
    computeTwoBodyEV(myGS,sigX,sigX,l,*out);
  *out<<endl;
  *out<<"% Y(k)Y(k+"<<l<<") \t";
  computeTwoBodyEV(myGS,sigY,sigY,l,*out);
  *out<<endl;
  *out<<"% Z(k)Z(k+"<<l<<") \t";
  computeTwoBodyEV(myGS,sigZ,sigZ,l,*out);
  *out<<endl;
  }


  out->close();
  delete out;
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
  complex_t traceRho=contractor.contract(rho,aux);

  // WARNING: No check of dimensions done!!
  mwArray op1_=reshape(op1,Indices(d*d,1,1));
  op1_.conjugate(); // It gets conjugated in the scalar product! 

  // Now, site by site, replace the identity by the operator,
  // contract, write the result to os, and replace the operator by the
  // identity again
  for(int k=0;k<L;k++){
    aux.setA(k,op1_); // operator on site k
    complex_t valK=contractor.contract(rho,aux);
    os<<valK/traceRho<<"\t";
    aux.setA(k,basic); // operator removed
  }
}

void computeTwoBodyEV(const MPS& rho,const mwArray& op1,
		      const mwArray& op2,int l,ostream& os){
  int L=rho.getLength();
  MPS aux(L,1,d*d);
  mwArray basic=identityMatrix(d); // basic piece of the MPO
  basic.reshape(Indices(d*d,1,1));
  for(int k=0;k<L;k++){
    aux.setA(k,basic);
  }

  Contractor& contractor=Contractor::theContractor();
  complex_t traceRho=contractor.contract(rho,aux);

  // WARNING: No check of dimensions done!!
  mwArray op1_=reshape(op1,Indices(d*d,1,1));
  op1_.conjugate(); // It gets conjugated in the scalar product! 
  mwArray op2_=reshape(op2,Indices(d*d,1,1));
  op2_.conjugate(); // It gets conjugated in the scalar product! 

  // Now, site by site, replace the identity at k by the operator op1,
  // the identity in k+l by op2,
  // contract, write the result to os, and replace the operators by the
  // identity again
  for(int k=0;k+l<L;k++){
    aux.setA(k,op1_); // operator on site k
    aux.setA(k+l,op2_); // operator on site k+l
    complex_t valK=contractor.contract(rho,aux);
    os<<valK/traceRho<<"\t";
    aux.setA(k,basic); // operator1 removed
    aux.setA(k+l,basic); // operator2 removed
  }

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

