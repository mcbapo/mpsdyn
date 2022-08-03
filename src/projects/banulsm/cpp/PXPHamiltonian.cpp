/**
  \file PXPHamiltonian.cpp
  Hamiltonian for a stochastic model

  \author Mari-Carmen Banuls
  \date 15/4/2019

  Add getUMPO function for obc
  Add coupling constant to do rescaling
  \modified Yilun Yang
  \date 14/5/2019
*/

#include "PXPHamiltonian.h"
#include "Contractor.h"


using namespace std;
using namespace shrt;

PXPHamiltonian::PXPHamiltonian(int L_,double J_,double lambda_,bool pbc_,double offset_):d(2),J(J_),lambda(lambda_),L(L_),hamil(L),offset(offset_),pbc(pbc_),proj(L){
    initZ();
    initHMPO();
    initProjMPO();
    //cout << localSigMPS->getA(k).getA();
}

PXPHamiltonian::~PXPHamiltonian(){};

void PXPHamiltonian::initZ(){
    int nrOps=(lambda!=0.)?4:3;
    // basic spin operators appearing
    mwArray sig0=identityMatrix(d);
    complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};
    mwArray sig1(Indices(d,d),datax);
    complex_t dataP[]={ONE_c,ZERO_c,ZERO_c,ZERO_c};
    mwArray operP(Indices(d,d),dataP);
    mwArray operQ=sig0-operP;
    Z=mwArray(Indices(nrOps,d,d));Z.fillWithZero();
    for(int i1=0;i1<d;i1++)
        for(int i2=0;i2<d;i2++){
            Z.setElement(sig0.getElement(Indices(i1,i2)),Indices(0,i1,i2));
            Z.setElement(operP.getElement(Indices(i1,i2)),Indices(1,i1,i2));
            Z.setElement(sig1.getElement(Indices(i1,i2)),Indices(2,i1,i2));
            if(nrOps==4)
                Z.setElement(operQ.getElement(Indices(i1,i2)),Indices(3,i1,i2));
        }
    Z.reshape(Indices(nrOps,d*d));
}

void PXPHamiltonian::initHMPO(){
    if(pbc){
        initHMPOpbc();return;
    }
    int D=(lambda==0.)?4:5;
    int nrOps=Z.getDimension(0); // 4 if lambda!=0
    for(int k=0;k<L;k++){
        // only first, second and last are distinct
        if(k>1&&k<L-1){
            hamil.setOp(k,&hamil.getOp(k-1));
        }
        else{
            //cout<<"Creating new Operator for site "<<k<<" of PXPH MPO"<<endl;
            int Dl=D; int Dr=D;
            if(k==0) Dl=1;
            if(k==L-1) Dr=1;
            mwArray C(Indices(Dl,Dr,nrOps));C.fillWithZero();
            // Set elements      
            if((k<L-1)&&(k>0)){ // normal ones
                C.setElement(ONE_c,Indices(0,0,0)); // nothing yet
                C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // After everything
                // First of the trio	
                C.setElement(ONE_c,Indices(0,1,1));
                // Second of the trio 
                C.setElement(ONE_c,Indices(1,2,2));
                // Last one
                C.setElement(-J*ONE_c,Indices(2,Dr-1,1));
                if(lambda!=0.){ // penalty term
                    C.setElement(lambda*ONE_c,Indices(0,3,3));
                    C.setElement(ONE_c,Indices(3,Dr-1,3));
                }
            }
            // Now the edges
            if(k==0){
                C.setElement(ONE_c,Indices(0,0,0)); // nothing yet
                // First of the trio	
                C.setElement(ONE_c,Indices(0,1,1));
                // Second of the trio (special)
                C.setElement(ONE_c,Indices(0,2,2));
                if(lambda!=0.){ // penalty term (first)
                    C.setElement(lambda*ONE_c,Indices(0,3,3));
                }
            }
            if(k==L-1){
                C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // After everything
                // Second of the trio (special)
                C.setElement(-J*ONE_c,Indices(1,Dr-1,2));
                // Last one
                C.setElement(-J*ONE_c,Indices(2,Dr-1,1));
                if(lambda!=0.){ // penalty term (second)
                    C.setElement(ONE_c,Indices(3,0,3));
                }
            }
            if(offset!=0)
                C.setElement((offset/L)*ONE_c,Indices(0,Dr-1,0)); // the offset term
            // Reshape, multiply operators and set in MPO
            C.reshape(Indices(Dl*Dr,nrOps));
            C.multiplyRight(Z);
            C.reshape(Indices(Dl,Dr,d,d));
            C.permute(Indices(3,1,4,2));
            hamil.setOp(k,new Operator(C),true);
        }
    }
}

void PXPHamiltonian::initHMPOpbc(){
    int D=(lambda==0)?6:8;
    int nrOps=Z.getDimension(0);
    for(int k=0;k<L;k++){
        //cout<<"Creating new Operator for site "<<k<<" of PXPH MPO"<<endl;
        int Dl=D; int Dr=D;
        if(k==0) Dl=1;
        if(k==L-1) Dr=1;
        mwArray C(Indices(Dl,Dr,nrOps));C.fillWithZero();
        // Set elements
        if(k==0){
            C.setElement(ONE_c,Indices(0,0,0)); // nothing yet
            C.setElement(ONE_c,Indices(0,1,1)); // first
            C.setElement(ONE_c,Indices(0,4,2)); // second
            C.setElement(-J*ONE_c,Indices(0,3,1)); // last
            if(lambda!=0.){ // penalty term
                C.setElement(lambda*ONE_c,Indices(0,5,3)); // first (regular)
                C.setElement(ONE_c,Indices(0,6,3)); // second (over boundary)
            }
        }
        else{
            if(k==L-1){
                C.setElement(-J*ONE_c,Indices(2,0,1)); // last
                C.setElement(ONE_c,Indices(4,0,1)); // first
                C.setElement(ONE_c,Indices(3,0,2)); // second
                C.setElement(ONE_c,Indices(Dl-1,0,0)); // nothing
                if(lambda!=0.){ // penalty term
                    C.setElement(lambda*ONE_c,Indices(6,0,3)); // first (over boundary)
                    C.setElement(ONE_c,Indices(5,0,3)); // second (regular)
                }
            }
            else{ // the rest, with two exceptions for sites 1 and L-2
                C.setElement(ONE_c,Indices(0,0,0)); // nothing yet
                if(k>1) C.setElement(ONE_c,Indices(Dl-1,Dr-1,0)); // After everything
                if(k<L-2) C.setElement(ONE_c,Indices(0,1,1)); // the first one
                C.setElement(ONE_c,Indices(1,2,2)); // the (normal) second one
                if(k>1) C.setElement(-J*ONE_c,Indices(2,Dr-1,1)); // the (normal) third one
                if(k>1) C.setElement(ONE_c,Indices(4,4,0)); // the one that starts at the last 
                else C.setElement(-J*ONE_c,Indices(4,4,1)); 
                if(k<L-2) C.setElement(ONE_c,Indices(3,3,0)); // the one that finishes at the first 
                else C.setElement(ONE_c,Indices(3,3,1)); 
                if(lambda!=0.){ // penalty term
                    C.setElement(lambda*ONE_c,Indices(0,5,3));
                    C.setElement(ONE_c,Indices(5,Dr-1,3));
                    C.setElement(ONE_c,Indices(6,6,0)); // the one that starts at the last 
                }
            }
        }
        if(offset!=0)
            C.setElement((offset/L)*ONE_c,Indices(0,Dr-1,0)); // the offset term
        // Reshape, multiply operators and set in MPO
        C.reshape(Indices(Dl*Dr,nrOps));
        C.multiplyRight(Z);
        C.reshape(Indices(Dl,Dr,d,d));
        C.permute(Indices(3,1,4,2));
        hamil.setOp(k,new Operator(C),true);
    }
}

void PXPHamiltonian::initProjMPO(){
    if(pbc){
        initProjMPOpbc();return;
    }
    int D=2;
    mwArray C(Indices(d,D,d,D));C.fillWithZero();
    C.setElement(ONE_c,Indices(0,0,0,0));
    C.setElement(ONE_c,Indices(1,0,1,1));
    C.setElement(ONE_c,Indices(0,1,0,0));
    proj.setOp(1,new Operator(C),true);
    for(int k=2;k<L-1;k++){
        proj.setOp(k,&proj.getOp(k-1));
    }
    // Edges
    mwArray Cl(Indices(d,1,d,D));Cl.fillWithZero();
    Cl.setElement(ONE_c,Indices(0,0,0,0));
    Cl.setElement(ONE_c,Indices(1,0,1,1));
    proj.setOp(0,new Operator(Cl),true);
    mwArray Cr(Indices(d,D,d,1));Cr.fillWithZero();
    Cr.setElement(ONE_c,Indices(0,1,0,0));
    Cr.setElement(ONE_c,Indices(0,0,0,0));
    Cr.setElement(ONE_c,Indices(1,0,1,0));
    proj.setOp(L-1,new Operator(Cr),true);
}

void PXPHamiltonian::initProjMPOpbc(){  
    int D=2;
    mwArray C(Indices(d,D,d,D));C.fillWithZero();
    mwArray auxInd=identityMatrix(D);auxInd.reshape(Indices(1,D*D));
    C.setElement(ONE_c,Indices(0,0,0,0));
    C.setElement(ONE_c,Indices(1,0,1,1));
    C.setElement(ONE_c,Indices(0,1,0,0));
    C.reshape(Indices(d*D*d*D,1));
    C.multiplyRight(auxInd);
    C.reshape(Indices(d,D,d,D,D,D));
    C.permute(Indices(1,2,5,3,4,6));
    C.reshape(Indices(d,D*D,d,D*D));
    proj.setOp(1,new Operator(C),true);
    for(int k=2;k<L-1;k++){
        proj.setOp(k,&proj.getOp(k-1));
    }
    // And the edges
    mwArray Cl(Indices(d,1,d,D,D));Cl.fillWithZero();
    Cl.setElement(ONE_c,Indices(0,0,0,0,0));
    Cl.setElement(ONE_c,Indices(1,0,1,1,1));
    Cl.reshape(Indices(d,1,d,D*D));
    proj.setOp(0,new Operator(Cl),true);
    mwArray Cr(Indices(d,D,D,d,1));Cr.fillWithZero();
    Cr.setElement(ONE_c,Indices(0,0,0,0,0));
    Cr.setElement(ONE_c,Indices(1,0,0,1,0));
    Cr.setElement(ONE_c,Indices(0,1,0,0,0));
    Cr.setElement(ONE_c,Indices(0,0,1,0,0));
    Cr.setElement(ONE_c,Indices(0,1,1,0,0));
    Cr.reshape(Indices(d,D*D,d,1));
    proj.setOp(L-1,new Operator(Cr),true);
}









void PXPHamiltonian::getThreeBodyTermExponential(mwArray &Ol, mwArray &Om, mwArray &Or,
        complex_t delta,
        int pos) const {
    mwArray H123;
    computeThreeBodyTerm(H123, pos);
    getThreeBodyTermExponential(Ol, Om, Or, H123, delta);
}

void PXPHamiltonian::getThreeBodyTermExponential(mwArray &Ol, mwArray &Om, mwArray &Or,
        const mwArray &H123,
        complex_t delta) const {
    // Now take the matrix exponential
    mwArray expH;
    wrapper::expm(H123, expH, delta);
    // putForMatlab(cout,expH,"expH12");
    // cout<<"Exponential of the h12="<<expH<<endl;

    // Obtained with indices (i'j'k')(ijk) => permute to (i'i)(j'jk'k) for first svd
    expH.reshape(Indices(d, d, d, d, d, d));
    expH.permute(Indices(1, 4, 2, 5, 3, 6));
    expH.reshape(Indices(d * d, d * d * d * d));
    // And compute the SVD
    mwArray S; // place for the singular values
    mwArray Omr;
    int nr = 0;
    double tol = 0.;
    // cout<<"Before decomposition, expH="<<expH<<endl;
    wrapper::svd(expH, tol, nr, Ol, S, Omr);
    // cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<", S="<<S<<endl;
    // redistribute the singular values
    S = sqrt(S);
    Ol.multiplyRight(S);
    Omr.multiplyLeft(S);
    Ol.reshape(Indices(d, 1, d, -1));
    Ol.permute(Indices(1, 2, 4, 3));
    // Result is in shape (i', 1, Dl, i)

    // Do second svd
    int Dm = Ol.getDimension(2);
    Omr.reshape(Indices(-1, d * d));
    wrapper::svd(Omr, tol, nr, Om, S, Or);
    S = sqrt(S);
    Om.multiplyRight(S);
    Or.multiplyLeft(S);
    Om.reshape(Indices(Dm, d, d, -1));
    Or.reshape(Indices(-1, d, d, 1));
    Om.permute(Indices(2, 1, 4, 3));
    Or.permute(Indices(2, 1, 4, 3));
    // Result is in shape (i', Dr, 1, i)
    // Dl, Dr is the bond dimension of Om
}

void PXPHamiltonian::computeThreeBodyTerm(mwArray &result, int pos) const {
    if (pos >= L - 3 || L <= 1) {
        cout << "ERROR: PXPHamiltonian::computeThreeBodyTerm for site "
            << pos << " when L=" << L << endl;
        exit(212);
    }
    computeThreeBodyTerm(result);
}

void PXPHamiltonian::computeThreeBodyTerm(mwArray &result) const {
    static mwArray PXP;
    static bool firstCall(true);
    if (firstCall) {        
        complex_t dataPXP[] = {
            ZERO_c, ZERO_c, ONE_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c,
            ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c,
            ONE_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c,
            ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c,
            ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c,
            ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c,
            ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c,
            ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c, ZERO_c
        };                   
        PXP = J * mwArray(Indices(d * d * d, d * d * d), dataPXP);
        firstCall = false;
    }
    result = PXP;
}


void PXPHamiltonian::getTwoBodyTermExponential(mwArray &Ol, mwArray &Or,
        complex_t delta,
        int pos) const {
    mwArray H12;
    computeTwoBodyTerm(H12, pos);
    getTwoBodyTermExponential(Ol, Or, H12, delta);
}

void PXPHamiltonian::getTwoBodyTermExponential(mwArray &Ol, mwArray &Or,
        const mwArray &H12,
        complex_t delta) const {
    // Now take the matrix exponential
    mwArray expH;
    wrapper::expm(H12, expH, delta);
    // putForMatlab(cout,expH,"expH12");
    // cout<<"Exponential of the h12="<<expH<<endl;

    // Obtained with indices (i'j')(ij) => permute to (i'i)(j'j) for svd
    expH.reshape(Indices(d, d, d, d));
    expH.permute(Indices(1, 3, 2, 4));
    expH.reshape(Indices(d * d, d * d));
    // And compute the SVD
    mwArray S; // place for the singular values
    int nr = 0;
    double tol = 0.;
    // cout<<"Before decomposition, expH="<<expH<<endl;
    wrapper::svd(expH, tol, nr, Ol, S, Or);
    // cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<", S="<<S<<endl;
    // redistribute the singular values
    S = sqrt(S);
    Ol.multiplyRight(S);
    Or.multiplyLeft(S);
    // reshape for operators (i'x 1 x i x alpha ))
    Ol.reshape(Indices(d, 1, d, -1));
    Ol.permute(Indices(1, 2, 4, 3));
    // Result is in shape (i', 1, D, i)
    Or.reshape(Indices(-1, d, d, 1));
    Or.permute(Indices(2, 1, 4, 3));
    // Result is in shape (i', D, 1, i)
    // cout<<"After decomposition, Ol="<<Ol<<", Or="<<Or<<endl;
    // putForMatlab(cout,Ol,"Ol");
    // putForMatlab(cout,Or,"Or");
}

void PXPHamiltonian::computeTwoBodyTerm(mwArray &result, int pos) const {
    if (pos != L-2 && pos != 0) {
        cout << "ERROR: PXPHamiltonian::computeTwoBodyTerm for site "
            << pos << " when L=" << L << endl;
        exit(212);
    }	
    computeTwoBodyTerm(pos == 0 ? 0 : 1, result);
}


void PXPHamiltonian::computeTwoBodyTerm(int lr, mwArray &result) const {
    static mwArray PX;
    static mwArray XP;
    static bool firstCall(true);
    if (firstCall) {
        complex_t dataXP[]={
            ZERO_c, ZERO_c, ONE_c, ZERO_c,
            ZERO_c, ZERO_c, ZERO_c, ZERO_c,
            ONE_c, ZERO_c, ZERO_c, ZERO_c,
            ZERO_c, ZERO_c, ZERO_c, ZERO_c
        };
        complex_t dataPX[]={
            ZERO_c, ONE_c, ZERO_c, ZERO_c,
            ONE_c, ZERO_c, ZERO_c, ZERO_c,
            ZERO_c, ZERO_c, ZERO_c, ZERO_c,
            ZERO_c, ZERO_c, ZERO_c, ZERO_c
        };
        XP = J * mwArray(Indices(d * d, d * d), dataXP);
        PX = J * mwArray(Indices(d * d, d * d), dataPX);
        firstCall = false;
    }
    if (lr == 0){
        result = XP;
    }
    else{
        result = PX;
    }
}

void PXPHamiltonian::getUMPO(MPO &expH, complex_t delta) const {
    mwArray PXPl, PXPm, PXPr, PXPlHalf, PXPmHalf, PXPrHalf;
    mwArray XPlHalf, XPrHalf;
    
    //cout << "Create local exp ops" << endl;
    getTwoBodyTermExponential(XPlHalf, XPrHalf, -delta/2, 0);
    getThreeBodyTermExponential(PXPl, PXPm, PXPr, -delta, 1);
    getThreeBodyTermExponential(PXPlHalf, PXPmHalf, PXPrHalf, -delta/2, 1);

    /*
    cout << "XPlHalf" << XPlHalf.getDimensions() << endl;
    cout << "XPrHalf" << XPrHalf.getDimensions() << endl;
    out << "PXPlHalf" << PXPlHalf.getDimensions() << endl;
    out << "PXPmHalf" << PXPmHalf.getDimensions() << endl;
    out << "PXPrHalf" << PXPrHalf.getDimensions() << endl;
    */
    
    expH.initLength(L);

    int PXPDl = PXPm.getDimension(1), PXPDr = PXPm.getDimension(2);
    int XPD = XPlHalf.getDimension(2);
    
    //cout << "Create UMPO" << endl;
    // Set first Op
    XPlHalf.reshape(Indices(-1, d));
    PXPl.reshape(Indices(d, -1));
    mwArray Op = XPlHalf * PXPl;
    Op.reshape(Indices(-1, d));
    XPlHalf.reshape(Indices(d, -1));
    Op = Op * XPlHalf;
    Op.reshape(Indices(d, -1, d, 1));
    Op.permute(Indices(1, 4, 3, 2));
    expH.setOp(0, new Operator(Op), true);
   
    // Set second Op
    XPrHalf.reshape(Indices(-1, d));
    PXPlHalf.reshape(Indices(d, -1));
    Op = XPrHalf * PXPlHalf;
    Op.reshape(Indices(-1, d));
    PXPm.reshape(Indices(d, -1));
    mwArray Op2 = Op * PXPm;
    Op2.reshape(Indices(-1, d));
    Op.reshape(Indices(d, -1));
    Op = Op2 * Op;
    Op.reshape(Indices(d, XPD, PXPDl, PXPDl, PXPDr, XPD, PXPDl, d));
    Op.permute(Indices(1, 2, 4, 6, 8, 3, 5, 7));
    Op.reshape(Indices(d, XPD*PXPDl*XPD, d, -1));
    expH.setOp(1, new Operator(Op), true);


    // Set Ops for odd k
    PXPrHalf.reshape(Indices(-1, d));
    PXPlHalf.reshape(Indices(d, -1));
    Op = PXPrHalf * PXPlHalf;
    Op.reshape(Indices(-1, d));
    PXPm.reshape(Indices(d, -1));
    Op2 = Op * PXPm;
    Op2.reshape(Indices(-1, d));
    Op.reshape(Indices(d, -1));
    Op = Op2 * Op;
    Op.reshape(Indices(d, PXPDr, PXPDl, PXPDl, PXPDr, PXPDr, PXPDl, d));
    Op.permute(Indices(1, 2, 4, 6, 8, 3, 5, 7));
    Op.reshape(Indices(d, PXPDr*PXPDl*PXPDr, d, -1));
    for (int k = 3; k < L - 2; k += 2){
        if (k==3){
            expH.setOp(k, new Operator(Op), true);
        }
        else{
            expH.setOp(k, &expH.getOp(k-2));
        }
    }
    
    // Set Ops for even k
    PXPmHalf.reshape(Indices(-1, d));
    PXPr.reshape(Indices(d, -1));
    Op = PXPmHalf * PXPr;
    Op.reshape(Indices(-1, d));
    PXPl.reshape(Indices(d, -1));
    Op = Op * PXPl;
    Op.reshape(Indices(-1, d));
    PXPmHalf.reshape(Indices(d, -1));
    Op = Op * PXPmHalf;
    Op.reshape(Indices(d, PXPDl, PXPDr, PXPDr, PXPDl, PXPDl, PXPDr, d));
    Op.permute(Indices(1, 2, 4, 6, 8, 3, 5, 7));
    Op.reshape(Indices(d, PXPDl*PXPDr*PXPDl, d, -1));
    for (int k = 2; k < L - 2; k += 2){
        if (k==2){
            expH.setOp(k, new Operator(Op), true);
        }
        else{
            expH.setOp(k, &expH.getOp(k-2));
        }
    }
    
    // Set last two Ops
    if (L % 2 == 0){
        mwArray PXl, PXr;
        getTwoBodyTermExponential(PXl, PXr, -delta, L-2);
        int PXD = PXl.getDimension(2);

        PXPmHalf.reshape(Indices(-1, d));
        PXPr.reshape(Indices(d, -1));
        Op = PXPmHalf * PXPr;
        Op.reshape(Indices(-1, d));
        PXl.reshape(Indices(d, -1));
        Op = Op * PXPl;
        Op.reshape(Indices(-1, d));
        PXPmHalf.reshape(Indices(d, -1));
        Op = Op * PXPmHalf;
        Op.reshape(Indices(d, PXPDl, PXPDr, PXPDr, PXD, PXPDl, PXPDr, d));
        Op.permute(Indices(1, 2, 4, 6, 8, 3, 5, 7));
        Op.reshape(Indices(d, PXPDl*PXPDr*PXPDl, d, -1));
        expH.setOp(L-2, new Operator(Op), true);

        PXPrHalf.reshape(Indices(-1, d));
        PXr.reshape(Indices(d, -1));
        mwArray Op = PXPrHalf * PXr;
        Op.reshape(Indices(-1, d));
        PXPrHalf.reshape(Indices(d, -1));
        Op = Op * PXPrHalf;
        Op.reshape(Indices(d, -1, d, 1));
        expH.setOp(L-1, new Operator(Op), true);
    }
    else{
        mwArray PXlHalf, PXrHalf;
        getTwoBodyTermExponential(PXlHalf, PXrHalf, -delta/2, L-2);
        int PXD = PXlHalf.getDimension(2);

        PXPrHalf.reshape(Indices(-1, d));
        PXlHalf.reshape(Indices(d, -1));
        Op = PXPrHalf * PXlHalf;
        Op.reshape(Indices(-1, d));
        PXPm.reshape(Indices(d, -1));
        Op2 = Op * PXPm;
        Op2.reshape(Indices(-1, d));
        Op.reshape(Indices(d, -1));
        Op = Op2 * Op;
        Op.reshape(Indices(d, PXPDr, PXD, PXPDl, PXPDr, PXPDr, PXD, d));
        Op.permute(Indices(1, 2, 4, 6, 8, 3, 5, 7));
        Op.reshape(Indices(d, PXPDr*PXPDl*PXPDr, d, -1));
        expH.setOp(L-2, new Operator(Op), true);


        PXrHalf.reshape(Indices(-1, d));
        PXPr.reshape(Indices(d, -1));
        mwArray Op = PXrHalf * PXPr;
        Op.reshape(Indices(-1, d));
        PXrHalf.reshape(Indices(d, -1));
        Op = Op * PXrHalf;
        Op.reshape(Indices(d, -1, d, 1));
        expH.setOp(L-1, new Operator(Op), true);
    }
    MPS expHMPS(L, 1, 4);
    MPSfromMPO(expH, expHMPS);
    MPS auxMPS(expHMPS);
    Contractor& contractor=Contractor::theContractor();
    contractor.optimizeMPS(expHMPS, auxMPS, 4);
    MPOfromMPS(auxMPS, expH);
    
}


