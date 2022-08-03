#include "HeisenbergAncillaryNoiseHamiltonian.h"
#include "DoubleOperator.h"

using namespace shrt;

HeisenbergAncillaryNoiseHamiltonian::HeisenbergAncillaryNoiseHamiltonian(int L_,double J_,double B_):
    HeisenbergHamiltonian(L_,2),J1x(1.),J1y(1.),J1z(1.),J2x(J_),J2y(J_),J2z(J_),B(B_),offset(0.),TI(true){
        d=2;
        // vectors of parameters here are all just equal to a single number
        // Jxi.push_back(J_);
        // Jyi.push_back(J_);
        // Jzi.push_back(J_);
        // hi.push_back(B_);
        initHMPO();
        // cout<<"Initialized H MPO for HeisenbergAncillaryNoise, with mpo "<<hamil<<endl;
    }

HeisenbergAncillaryNoiseHamiltonian::HeisenbergAncillaryNoiseHamiltonian(int L_,double J1x_,double J1y_,double J1z_,double J2x_,double J2y_,double J2z_,double B_, double offset_):
    HeisenbergHamiltonian(L_,2),J1x(J1x_),J1y(J1y_),J1z(J1z_),J2x(J2x_),J2y(J2y_),J2z(J2z_),B(B_),offset(offset_),TI(true){
        d=2;
        initHMPO();
    }

HeisenbergAncillaryNoiseHamiltonian::HeisenbergAncillaryNoiseHamiltonian(int L_,double J_,
        const vector<double>& B_):
    HeisenbergHamiltonian(L_,2),J1x(1.),J1y(1.),J1z(1.),J2x(J_),J2y(J_),J2z(J_),B(0),offset(0.),TI(false){
        d=2;
        int len=B_.size();
        if(len!=L){
            cout<<"Error: trying to construct non-TI HeisenbergAncillaryNoiseHamiltonian for L="<<L<<" with a vector of Bs="<<B_<<endl;
            exit(1);
        }
        for(int k=0;k<L;k++)
            Bs.push_back(B_[k]);
        initHMPOnonTI();
        cout<<"Initialized H MPO for HeisenbergAncillaryNoise, with mpo "<<hamil<<endl;
    }

HeisenbergAncillaryNoiseHamiltonian::~HeisenbergAncillaryNoiseHamiltonian(){
    if(!TI) Bs.clear();
};

void HeisenbergAncillaryNoiseHamiltonian::getExponentialMPOeven(MPO& expHe,complex_t delta,bool sigma)  const{
    cout<<"HeisenbergAncillaryNoiseHamiltonian::getExponentialMPOeven(delta="
        <<delta<<","<<(sigma?"sigma)":"tau)")<<endl;
    getExponentialMPO(expHe,delta,true,sigma); 
}

void HeisenbergAncillaryNoiseHamiltonian::getExponentialMPOodd(MPO& expHo,complex_t delta,bool sigma) const {
    cout<<"HeisenbergAncillaryNoiseHamiltonian::getExponentialMPOodd(delta="
        <<delta<<","<<(sigma?"sigma)":"tau)")<<endl;
    getExponentialMPO(expHo,delta,false,sigma); 
}

void HeisenbergAncillaryNoiseHamiltonian::getExponentialMPO(MPO& expH,complex_t delta,bool even,bool sigma) const {
    /** The even and odd terms for sigma and tau are the same, and
      identical to the Heisenberg ones, only stretched. */
    expH.initLength(2*L);
    double Jx_=sigma?J1x:J2x;
    double Jy_=sigma?J1y:J2y;
    double Jz_=sigma?J1z:J2z;
    double h=0.;
    mwArray idPhys=identityMatrix(d);
    idPhys.reshape(Indices(d,1,d,1)); // for the idOp
    Operator* idOpPhys=new Operator(idPhys);
    bool someterm=containsOperators(even,sigma);
    if(someterm){
        mwArray H12,Ol,Or; // Compute the exponential of a two-body Heisenberg term
        computeTwoBodyTerm(H12,Jx_,Jy_,Jz_,h,h);
        getTwoBodyTermExponential(Ol,Or,H12,delta); // Return the already SVDd exponential
        // Additionally I need the identity in the physical dof, and the
        // identity in both
        int D=Ol.getDimension(3); // max dim of the virtual bond
        mwArray idVir=identityMatrix(D);
        Operator* operOl=new Operator(Ol);
        Operator* operOr=new Operator(Or);
        idPhys.reshape(Indices(d*d,1));
        idVir.reshape(Indices(1,D*D));
        idPhys.multiplyRight(idVir);
        idPhys.reshape(Indices(d,d,D,D));idPhys.permute(Indices(1,3,2,4));
        Operator* idOpBoth=new Operator(idPhys);
        setExponentialOperators(expH,operOl,operOr,idOpPhys,idOpBoth,even,sigma);
    }
    else{
        cout<<"WARNING: constructing an exponentialMPO ("<<(even?"even":"odd")<<","<<(sigma?"sigma":"tau")
            <<") which is identically 1: would be better not to apply it at all"<<endl;
        setExponentialOperators(expH,NULL,NULL,idOpPhys,NULL,even,sigma);
    }
}

void HeisenbergAncillaryNoiseHamiltonian::getExponentialMPOmixed(MPO& expH,complex_t delta)  const {
    // The term is an ising interaction, sigmaz sigmaz. Therefore, I
    // know the exponential exactly, and can construct the Ol and Or by
    // hand (as in Ising).
    int D=2; // bond dimension (max)
    expH.initLength(2*L);
    mwArray C(Indices(1,D,4));
    if(TI){
        C.setElement(sqrt(cosh(B*.25*delta)),Indices(0,0,0)); // coeff with Id
        C.setElement(2*sqrt(sinh(B*.25*delta)),Indices(0,1,3)); // coeff with Z (in Z, there is Sz, not sigma_z)
        //cout<<"Coefficients in the exponent for site "<<pos<<" eps="<<eps<<endl;
        C.reshape(Indices(1*D,4));
        //cout<<"C tensor "<<reshape(resize(C,Indices(Dl*Dr,2)),Indices(Dl,Dr,2))<<endl;
        C.multiplyRight(Z); // this is now 1*D*d*d
        C.reshape(Indices(1,D,d,d));
        Operator* Ol=new Operator(C,Indices(3,1,4,2));
        Operator* Or=new Operator(C,Indices(3,2,4,1));
        // Now set them in the chain
        expH.initLength(2*L);
        expH.setOp(0,Ol,true);
        expH.setOp(1,Or,true);
        for(int k=1;k<L;k++){
            expH.setOp(2*k,Ol,false);
            expH.setOp(2*k+1,Or,false);
        }
    }
    else{ // if the Bs are not TI
        for(int k=0;k<L;k++){
            double Bx=Bs[k];
            C.setElement(sqrt(cosh(Bx*.25*delta)),Indices(0,0,0)); // coeff with Id
            C.setElement(2*sqrt(sinh(Bx*.25*delta)),Indices(0,1,3)); // coeff with Z (in Z, there is Sz, not sigma_z)
            //cout<<"Coefficients in the exponent for site "<<pos<<" eps="<<eps<<endl;
            mwArray Caux=reshape(C,Indices(1*D,4));
            //cout<<"C tensor "<<reshape(resize(C,Indices(Dl*Dr,2)),Indices(Dl,Dr,2))<<endl;
            Caux.multiplyRight(Z); // this is now 1*D*d*d
            Caux.reshape(Indices(1,D,d,d));
            Operator* Ol=new Operator(Caux,Indices(3,1,4,2));
            Operator* Or=new Operator(Caux,Indices(3,2,4,1));
            // Now set them in the chain
            expH.setOp(2*k,Ol,true);
            expH.setOp(2*k+1,Or,true);
        }
    }
}

void HeisenbergAncillaryNoiseHamiltonian::getDoubleExponentialMPOeven(MPO& expHe,complex_t delta,bool sigma)  const {
    // cout<<"HeisenbergAncillaryNoiseHamiltonian::getDoubleExponentialMPOeven(delta="
    //     <<delta<<","<<(sigma?"sigma)":"tau)")<<endl;
    getDoubleExponentialMPO(expHe,delta,true,sigma); 
}

void HeisenbergAncillaryNoiseHamiltonian::getDoubleExponentialMPOodd(MPO& expHo,complex_t delta,bool sigma)  const {
    // cout<<"HeisenbergAncillaryNoiseHamiltonian::getDoubleExponentialMPOodd(delta="
    //     <<delta<<","<<(sigma?"sigma)":"tau)")<<endl;
    getDoubleExponentialMPO(expHo,delta,false,sigma); 
}

void HeisenbergAncillaryNoiseHamiltonian::getDoubleExponentialMPO(MPO& expH,complex_t delta,bool even,bool sigma) const {
    /** The even and odd terms for sigma and tau are the same, and
      identical to the Heisenberg ones, only stretched. */
    cout<<"HeisenbergAncillaryNoiseHamiltonian::getDoubleExponentialMPO(delta="
        <<delta<<(even?",even,":",odd,")
        <<(sigma?"sigma)":"tau)")<<endl;
    expH.initLength(2*L);
    double Jx_=sigma?J1x:J2x;
    double Jy_=sigma?J1y:J2y;
    double Jz_=sigma?J1z:J2z;
    double h=0.;
    // Additionally I need the identity in the physical dof, and the identity in
    // both
    mwArray idPhys=identityMatrix(d);
    idPhys.reshape(Indices(d,1,d,1)); // for the idOp

    bool someterm=containsOperators(even,sigma);
    if(someterm){
        mwArray H12,Ol,Or; // Compute the exponential of a two-body Heisenberg term
        computeTwoBodyTerm(H12,Jx_,Jy_,Jz_,h,h);
        getTwoBodyTermExponential(Ol,Or,H12,delta); // Return the already
        // SVDd exponential
        mwArray OlA,OrA; // I need the complex ones acting on the ancillary
        // indices
        getTwoBodyTermExponential(OlA,OrA,H12,conjugate(delta)); 
        DoubleOperator* operOl=new DoubleOperator(Ol,permute(OlA,Indices(3,2,1,4)));
        DoubleOperator* operOr=new DoubleOperator(Or,permute(OrA,Indices(3,2,1,4)));
        DoubleOperator* idOpPhys=new DoubleOperator(idPhys,idPhys);
        int D=Ol.getDimension(3); // max dim of the virtual bond
        mwArray idVir=identityMatrix(D);
        idPhys.reshape(Indices(d*d,1));
        idVir.reshape(Indices(1,D*D));
        idPhys.multiplyRight(idVir);
        idPhys.reshape(Indices(d,d,D,D));idPhys.permute(Indices(1,3,2,4));
        DoubleOperator* idOpBoth=new DoubleOperator(idPhys,idPhys);
        setExponentialOperators(expH,operOl,operOr,idOpPhys,idOpBoth,even,sigma);
    }
    else{
        DoubleOperator* idOpPhys=new DoubleOperator(idPhys,idPhys);
        setExponentialOperators(expH,NULL,NULL,idOpPhys,NULL,even,sigma);
    }
}

void HeisenbergAncillaryNoiseHamiltonian::getDoubleExponentialMPOmixed(MPO& expH,complex_t delta) const{
    // The term is an ising interaction, sigmaz sigmaz. Therefoer, I
    // know the exponential exactly, and can construct the Ol and Or by
    // hand (as in Ising).
    int D=2; // bond dimension (max)
    mwArray C(Indices(1,D,4));
    mwArray CA(Indices(1,D,4)); // for the ancilla
    expH.initLength(2*L);
    if(TI){
        C.setElement(sqrt(cosh(B*.25*delta)),Indices(0,0,0)); // coeff with Id
        C.setElement(2*sqrt(sinh(B*.25*delta)),Indices(0,1,3)); // coeff with Sz (NOT sigmaz!!)
        CA.setElement(sqrt(cosh(B*.25*conjugate(delta))),Indices(0,0,0)); // coeff with Id
        CA.setElement(2*sqrt(sinh(B*.25*conjugate(delta))),Indices(0,1,3)); // coeff with Sz (NOT sigmaz!!)
        //cout<<"Coefficients in the exponent for site "<<pos<<" eps="<<eps<<endl;
        C.reshape(Indices(1*D,4));CA.reshape(Indices(1*D,4));
        //cout<<"C tensor "<<reshape(resize(C,Indices(Dl*Dr,2)),Indices(Dl,Dr,2))<<endl;
        C.multiplyRight(Z);CA.multiplyRight(Z); // this is now 1*D*d*d
        C.reshape(Indices(1,D,d,d));CA.reshape(Indices(1,D,d,d));
        // For the ancilla: transpose physical indices
        DoubleOperator* Ol=new DoubleOperator(permute(C,Indices(3,1,4,2)),
                permute(CA,Indices(4,1,3,2)));
        DoubleOperator* Or=new DoubleOperator(permute(C,Indices(3,2,4,1)),
                permute(CA,Indices(4,2,3,1)));
        // Now set them in the chain
        expH.setOp(0,Ol,true);
        expH.setOp(1,Or,true);
        for(int k=1;k<L;k++){
            expH.setOp(2*k,Ol,false);
            expH.setOp(2*k+1,Or,false);
        }
    }
    else{
        C.reshape(Indices(1*D,4));CA.reshape(Indices(1*D,4));
        for(int k=0;k<L;k++){
            double Bk=Bs[k];
            C.setElement(sqrt(cosh(Bk*.25*delta)),Indices(0,0)); // coeff with Id
            C.setElement(2*sqrt(sinh(Bk*.25*delta)),Indices(1,3)); // coeff with Sz (NOT sigmaz!!)
            CA.setElement(sqrt(cosh(Bk*.25*conjugate(delta))),Indices(0,0)); // coeff with Id
            CA.setElement(2*sqrt(sinh(Bk*.25*conjugate(delta))),Indices(1,3)); // coeff with Sz (NOT sigmaz!!)
            //cout<<"Coefficients in the exponent for site "<<pos<<" eps="<<eps<<endl;
            mwArray Cl=C*Z;Cl.reshape(Indices(1,D,d,d));
            mwArray Cr=CA*Z;Cr.reshape(Indices(1,D,d,d));
            // For the ancilla: transpose physical indices
            DoubleOperator* Ol=new DoubleOperator(permute(Cl,Indices(3,1,4,2)),
                    permute(Cr,Indices(4,1,3,2)));
            DoubleOperator* Or=new DoubleOperator(permute(Cl,Indices(3,2,4,1)),
                    permute(Cr,Indices(4,2,3,1)));
            expH.setOp(2*k,Ol,false);
            expH.setOp(2*k+1,Or,false);
        }
    }
}



void HeisenbergAncillaryNoiseHamiltonian::initHMPO(){
    //  cout<<"HeisenbergAncillaryNoiseHamiltonian::initHMPO()"<<endl;
    int D=8;
    // coefficients that appear, different from 1
    complex_t valJ1x=J1x*ONE_c;
    complex_t valJ1y=J1y*ONE_c;
    complex_t valJ1z=J1z*ONE_c;
    complex_t valJ2x=J2x*ONE_c;
    complex_t valJ2y=J2y*ONE_c;
    complex_t valJ2z=J2z*ONE_c;
    complex_t valB=B*ONE_c;
    // even and odd sites not on the edges
    mwArray Ce(Indices(D,D,4));
    mwArray Co(Indices(D,D,4));
    mwArray Ce0(Indices(1,D,4)); // the very first
    mwArray CeL; // the last even
    mwArray CoL(Indices(D,1,4)); // the very last
    
    Ce0.setElement(ONE_c,Indices(0,0,0)); // identity only
    Ce0.setElement(ONE_c,Indices(0,1,1)); // first part of XX
    Ce0.setElement(ONE_c,Indices(0,2,2)); // first part of YY
    Ce0.setElement(ONE_c,Indices(0,3,3)); // first part of ZZ (or B term)

    // Fill in the even tensors. 
    Ce.setElement(ONE_c,Indices(0,0,0)); // nothing yet
    Ce.setElement(ONE_c,Indices(0,1,1)); // sigx first part
    Ce.setElement(valJ1x,Indices(1,D-1,1)); // sigx second part
    Ce.setElement(ONE_c,Indices(0,2,2)); // sigy first part
    Ce.setElement(valJ1y,Indices(2,D-1,2)); // sigy second part
    Ce.setElement(ONE_c,Indices(0,3,3)); // sigz first part
    Ce.setElement(valJ1z,Indices(3,D-1,3)); // sigz second part
    Ce.setElement(ONE_c,Indices(4,4,0)); // transparent to a XX pair in taus
    Ce.setElement(ONE_c,Indices(5,5,0)); // transparent to a YY pair in taus
    Ce.setElement(ONE_c,Indices(6,6,0)); // transparent to a ZZ pair in taus
    Ce.setElement(ONE_c,Indices(D-1,D-1,0)); // identities until the end
    // The last even one has two different elements because it cannot
    // start a XX or YY pair, so I copy it and set them to zero by hand
    CeL=Ce;
    CeL.setElement(ZERO_c,Indices(0,0,0)); // only identities
    CeL.setElement(ZERO_c,Indices(0,1,1)); // sigx first part
    CeL.setElement(ZERO_c,Indices(0,2,2)); // sigy first part

    // Same for odd tensors
    Co.setElement(ONE_c,Indices(0,0,0)); // nothing yet
    Co.setElement(ONE_c,Indices(0,4,1)); // taux first part
    Co.setElement(valJ2x,Indices(4,D-1,1)); // taux second part
    Co.setElement(ONE_c,Indices(0,5,2)); // tauy first part
    Co.setElement(valJ2y,Indices(5,D-1,2)); // tauy second part
    Co.setElement(ONE_c,Indices(0,6,3)); // tauz first part
    Co.setElement(valJ2z,Indices(6,D-1,3)); // tauz second part
    Co.setElement(valB,Indices(3,D-1,3)); // tauz second part of B term!!
    Co.setElement(ONE_c,Indices(1,1,0)); // transparent to a XX pair in sigmas
    Co.setElement(ONE_c,Indices(2,2,0)); // transparent to a YY pair in sigmas
    Co.setElement(ONE_c,Indices(3,3,0)); // transparent to a ZZ pair in sigmas
    Co.setElement(ONE_c,Indices(D-1,D-1,0));  // identities until the end


    CoL.setElement(valJ2x,Indices(4,0,1)); // second part of XX
    CoL.setElement(valJ2y,Indices(5,0,2)); // second part of YY
    CoL.setElement(valJ2z,Indices(6,0,3)); // second part of ZZ
    CoL.setElement(valB,Indices(3,0,3)); // second part of B term
    CoL.setElement(ONE_c,Indices(D-1,0,0)); //final identity

    // Set offset
    if (offset!=0){
        Ce.setElement(offset/L/2.*ONE_c,Indices(0,D-1,0));     
        CeL.setElement(offset/L/2.*ONE_c,Indices(0,D-1,0));
        Co.setElement(offset/L/2.*ONE_c,Indices(0,D-1,0));
        CoL.setElement(offset/L/2.*ONE_c,Indices(0,0,0));
        Ce0.setElement(offset/L/2.*ONE_c,Indices(0,D-1,0));
    }

    // Reshape, multiply operators and set in MPO
    Ce.reshape(Indices(D*D,4));Ce.multiplyRight(Z);
    Ce.reshape(Indices(D,D,d,d));Ce.permute(Indices(3,1,4,2));
    Co.reshape(Indices(D*D,4));Co.multiplyRight(Z);
    Co.reshape(Indices(D,D,d,d));Co.permute(Indices(3,1,4,2));
    CeL.reshape(Indices(D*D,4));CeL.multiplyRight(Z);
    CeL.reshape(Indices(D,D,d,d));CeL.permute(Indices(3,1,4,2));
    Ce0.reshape(Indices(D,4));Ce0.multiplyRight(Z);
    Ce0.reshape(Indices(1,D,d,d));Ce0.permute(Indices(3,1,4,2));
    CoL.reshape(Indices(D*1,4));CoL.multiplyRight(Z);
    CoL.reshape(Indices(D,1,d,d));CoL.permute(Indices(3,1,4,2));

    hamil.initLength(2*L);
    // the edges
    hamil.setOp(0,new Operator(Ce0),true);
    hamil.setOp(2*(L-1),new Operator(CeL),true);
    hamil.setOp(2*(L-1)+1,new Operator(CoL),true);
    // the bulk (I will assume there are at least 4 sites (L>=2)
    Operator* Opo=new Operator(Co);
    hamil.setOp(1,Opo,true);
    if(L>2){
        Operator* Ope=new Operator(Ce);
        hamil.setOp(2,Ope,true);
        // now loop over the rest of the chain and set pointers
        for(int k=1;k<L;k++){
            if(k<L-1) // the very last is special
                hamil.setOp(2*k+1,Opo,false); // the odd ones
            if(k>1) // the first one is the original
                hamil.setOp(2*k,Ope,false); // the even ones
        }
    }
}


void HeisenbergAncillaryNoiseHamiltonian::initHMPOnonTI(){
    cout<<"HeisenbergAncillaryNoiseHamiltonian::initHMPOnonTI()"<<endl;
    int D=8;
    // coefficients that appear, different from 1
    complex_t valJ1x=sqrt(J1x)*ONE_c;double sgnJ1x=J1x<0?-1.:1.;
    complex_t valJ1y=sqrt(J1y)*ONE_c;double sgnJ1y=J1y<0?-1.:1.;
    complex_t valJ2x=sqrt(J2x)*ONE_c;double sgnJ2x=J2x<0?-1.:1.;
    complex_t valJ2y=sqrt(J2y)*ONE_c;double sgnJ2y=J2y<0?-1.:1.;
    complex_t valJ2z=sqrt(J2z)*ONE_c;double sgnJ2z=J2z<0?-1.:1.;
    // even and odd sites not on the edges
    mwArray Ce(Indices(D,D,4));
    mwArray Co(Indices(D,D,4));
    mwArray Ce0(Indices(1,D,4)); // the very first
    mwArray CeL; // the last even
    mwArray CoL(Indices(D,1,4)); // the very last
    // Fill in the even tensors. 
    Ce.setElement(ONE_c,Indices(0,0,0)); // nothing yet
    Ce.setElement(sgnJ1x*valJ1x,Indices(0,1,1)); // sigx first part
    Ce.setElement(valJ1x,Indices(1,D-1,1)); // sigx second part
    Ce.setElement(sgnJ1y*valJ1y,Indices(0,2,2)); // sigy first part
    Ce.setElement(valJ1y,Indices(2,D-1,2)); // sigy second part
    Ce.setElement(ONE_c,Indices(0,3,3)); // sigz first part
    Ce.setElement(J1z*ONE_c,Indices(3,D-1,3)); // sigz second part
    Ce.setElement(ONE_c,Indices(4,4,0)); // transparent to a XX pair in taus
    Ce.setElement(ONE_c,Indices(5,5,0)); // transparent to a YY pair in taus
    Ce.setElement(ONE_c,Indices(6,6,0)); // transparent to a ZZ pair in taus
    Ce.setElement(ONE_c,Indices(D-1,D-1,0)); // identities until the end
    // The last even one has two different elements because it cannot
    // start a XX or YY pair, so I copy it and set them to zero by hand
    CeL=Ce;
    CeL.setElement(ZERO_c,Indices(0,0,0)); // only identities
    CeL.setElement(ZERO_c,Indices(0,1,1)); // sigx first part
    CeL.setElement(ZERO_c,Indices(0,2,2)); // sigy first part

    // Same for odd tensors
    Co.setElement(ONE_c,Indices(0,0,0)); // nothing yet
    Co.setElement(sgnJ2x*valJ2x,Indices(0,4,1)); // taux first part
    Co.setElement(valJ2x,Indices(4,D-1,1)); // taux second part
    Co.setElement(sgnJ2y*valJ2y,Indices(0,5,2)); // tauy first part
    Co.setElement(valJ2y,Indices(5,D-1,2)); // tauy second part
    Co.setElement(sgnJ2z*valJ2z,Indices(0,6,3)); // tauz first part
    Co.setElement(valJ2z,Indices(6,D-1,3)); // tauz second part
    Co.setElement(ONE_c,Indices(1,1,0)); // transparent to a XX pair in sigmas
    Co.setElement(ONE_c,Indices(2,2,0)); // transparent to a YY pair in sigmas
    Co.setElement(ONE_c,Indices(3,3,0)); // transparent to a ZZ pair in sigmas
    Co.setElement(ONE_c,Indices(D-1,D-1,0));  // identities until the end
    // cout<<"Sizes of the different terms here: Ce"<<Ce.getDimensions()
    //     <<", Co "<<Co.getDimensions()
    //     <<", CeL "<<CeL.getDimensions()
    //     <<", Ce0 "<<Ce0.getDimensions()
    //     <<", CoL "<<CoL.getDimensions()
    //     <<", Z "<<Z.getDimensions()
    //     <<endl;
    // And the edges
    Ce0.setElement(ONE_c,Indices(0,0,0)); // identity only
    Ce0.setElement(sgnJ1x*valJ1x,Indices(0,1,1)); // first part of XX
    Ce0.setElement(sgnJ1y*valJ1y,Indices(0,2,2)); // first part of YY
    Ce0.setElement(ONE_c,Indices(0,3,3)); // first part of ZZ (or B term)

    CoL.setElement(valJ2x,Indices(4,0,1)); // second part of XX
    CoL.setElement(valJ2y,Indices(5,0,2)); // second part of YY
    CoL.setElement(valJ2z,Indices(6,0,3)); // second part of ZZ
    CoL.setElement(ONE_c,Indices(D-1,0,0)); //final identity
    CoL.setElement(Bs[L-1]*ONE_c,Indices(3,0,3)); // second part of B term for the last pair

    // Reshape, multiply operators and set in MPO
    Ce.reshape(Indices(D*D,4));Ce.multiplyRight(Z);
    Ce.reshape(Indices(D,D,d,d));Ce.permute(Indices(3,1,4,2));
    CeL.reshape(Indices(D*D,4));CeL.multiplyRight(Z);
    CeL.reshape(Indices(D,D,d,d));CeL.permute(Indices(3,1,4,2));
    Ce0.reshape(Indices(D,4));Ce0.multiplyRight(Z);
    Ce0.reshape(Indices(1,D,d,d));Ce0.permute(Indices(3,1,4,2));
    CoL.reshape(Indices(D*1,4));CoL.multiplyRight(Z);
    CoL.reshape(Indices(D,1,d,d));CoL.permute(Indices(3,1,4,2));

    hamil.initLength(2*L);
    // the edges
    hamil.setOp(0,new Operator(Ce0),true);
    hamil.setOp(2*(L-1),new Operator(CeL),true);
    hamil.setOp(2*(L-1)+1,new Operator(CoL),true);
    // the bulk (I will assume there are at least 4 sites (L>=2)
    // Set the even ones
    if(L>2){
        Operator* Ope=new Operator(Ce);
        hamil.setOp(2,Ope,true);
        cout<<"Set operator "<<2<<" to Ce"<<endl;
        // now loop over the rest of the chain and set pointers
        for(int k=2;k<L;k++){
            hamil.setOp(2*k,Ope,false); // the even ones
            cout<<"Set operator "<<2*k<<" to Ce"<<endl;
        }
    }

    // And now the odd, but these change
    for(int k=0;k<L-1;k++){
        // For the odd ones, the value of B changes for each site! 
        Co.setElement(Bs[k]*ONE_c,Indices(3,D-1,3)); // tauz second part of B term!!
        mwArray Co_=reshape(Co,Indices(D*D,4));Co_.multiplyRight(Z);
        Co_.reshape(Indices(D,D,d,d));Co_.permute(Indices(3,1,4,2));
        Operator* Opo=new Operator(Co_);
        hamil.setOp(2*k+1,Opo,true); // the odd ones
        cout<<"Set operator "<<2*k+1<<" to Co"<<endl;
    }
}


void HeisenbergAncillaryNoiseHamiltonian::setExponentialOperators(MPO& expH,
        Operator* operOl,
        Operator* operOr,
        Operator* idOpPhys,
        Operator* idOpBoth,
        bool even,
        bool sigma) const {
    // Now, depending on even or odd,and sigma or tau, place the operators and 
    // the required identities in the MPO
    int kfirst=even?0:1; // where it starts (in terms of the "double" chain)
    int klast=(L-1)-(L-kfirst)%2; // last occupied by the loop
    //  cout<<"kfirst is "<<kfirst<<", klast is "<<klast<<endl;
    // Now, if we are in a sigma term, everything happens on even sites
    int posOl1=sigma?2*kfirst:2*kfirst+1; // where the first Ol goes
    int posOrLast=sigma?2*klast:2*klast+1; // where the last Or goes
    //cout<<"posOl1 is "<<posOl1<<", posOrLast is "<<posOrLast<<endl;
    bool someterm=containsOperators(even,sigma);
    if(someterm){
        // Let's place the Ol operators
        //    DoubleOperator* operOl=new DoubleOperator(Ol,OlA);
        expH.setOp(posOl1,operOl,true); // first one owns the pointer
        //cout<<"first Ol->"<<posOl1<<endl;
        for(int pos=posOl1+4;pos<posOrLast;pos+=4){
            expH.setOp(pos,operOl,false);  
            //cout<<" ref to Ol->"<<pos<<endl;
        }
        // Let's place the Or operators
        //    DoubleOperator* operOr=new DoubleOperator(Or,OrA);
        expH.setOp(posOl1+2,operOr,true); // first one owns the pointer
        //cout<<"first Or->"<<posOl1+2<<endl;
        for(int pos=posOl1+6;pos<=posOrLast;pos+=4){
            expH.setOp(pos,operOr,false);  
            //cout<<" ref to Or->"<<pos<<endl;
        }
        // Now fill in with the proper identities
        // First simple idOp goes on 0 except for even sigma, when it goes on 3
        int firstPosId=(kfirst==1||!sigma)?0:3;
        expH.setOp(firstPosId,idOpPhys,true);
        //cout<<"first Id(d)->"<<firstPosId<<endl;
        int lastPosLeft=-1;
        for(int pos=1;pos<posOl1;pos++){ // simple identities at the beginning, if there are
            expH.setOp(pos,idOpPhys,false);
            //cout<<" ref to Id(d) (left)->"<<pos<<endl;
            lastPosLeft=pos;
        }
        // where does the second Id go=
        int secondIdPos=lastPosLeft>=0?lastPosLeft+4:firstPosId+4;
        for(int pos=secondIdPos;pos<2*L;pos+=4){ // periodically repeated Ids in between terms
            expH.setOp(pos,idOpPhys,false);
            //cout<<" ref to Id(d) (middle) ->"<<pos<<endl;
        }
        for(int pos=max(firstPosId+1,posOrLast+1);pos<2*L;pos++){ // and after the last term until the end of the chain (one wwill be probably set twice)
            expH.setOp(pos,idOpPhys,false);
            //cout<<" ref to Id(d) (right)->"<<pos<<endl;
        }
        // And now the composite Ids which go between each Ol-Or pair
        // first one after posOl1
        expH.setOp(posOl1+1,idOpBoth,true);
        //cout<<"first Id(d,D)->"<<posOl1+1<<endl;
        for(int pos=posOl1+1+4;pos<posOrLast;pos+=4){
            expH.setOp(pos,idOpBoth,false);
            //cout<<" ref to Id(d,D)->"<<pos<<endl;
        }
    }
    else{
        cout<<"WARNING: constructing an exponentialMPO ("<<(even?"even":"odd")<<","<<(sigma?"sigma":"tau")
            <<") which is identically 1: would be better not to apply it at all"<<endl;
        // the operator is just the identity
        //DoubleOperator* idOp=new DoubleOperator(idPhys,idPhys);
        expH.setOp(0,idOpPhys,true); // first one owns the pointer
        //cout<<"first Id(d)->"<<0<<endl;
        for(int k=1;k<2*L;k++){
            expH.setOp(k,idOpPhys,false); // all the rest hold a copy
            //cout<<" ref to Id(d)->"<<k<<endl;
        }
    }
}

bool HeisenbergAncillaryNoiseHamiltonian::containsOperators(bool even,bool sigma) const{
    int kfirst=even?0:1; // where it starts (in terms of the "double" chain)
    int klast=(L-1)-(L-kfirst)%2; // last occupied by the loop

    // Now, if we are in a sigma term, everything happens on even sites
    int posOl1=sigma?2*kfirst:2*kfirst+1; // where the first Ol goes
    int posOrLast=sigma?2*klast:2*klast+1; // where the last Or goes
    bool someterm=posOrLast<2*L; // does at least one term fit in the
    // chain?
    if(!someterm){
        cout<<"Detecting identity exponential for "<<(even?"even ":"odd ")
            <<(sigma?"sigma ":"tau ")<<"term, posOl1("<<posOl1
            <<","<<posOrLast<<")"<<endl;
    }
    return someterm;
}
