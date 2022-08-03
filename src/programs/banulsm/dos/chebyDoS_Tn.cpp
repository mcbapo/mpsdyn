/*  
 *  Calculate T_n(\tilde{H}) and save them.
 *  Output the truncation error of each term, Tr(T_n(H)) and Tr(T_n(H)O) for O in opList.
 *  */
#include <iomanip>
#include <math.h>
#include <sstream>
#include <string>

#include "Contractor.h"
#include "MPS.h"
#include "Properties.h"
#include "SpinMPO.h"
#include "misc.h"
#include "mwArray.h"

#define _USE_MATH_DEFINES // For PI
#define TMPSTORE 0		  // Using temporary storage for intermediate results
#define SSTR(x)                                                                \
    static_cast<std::ostringstream &>((std::ostringstream() << std::dec << x)) \
.str() 
// Convert variables to strings

// Included models
#include "HeisenbergHamiltonian.h"
#include "IsingHamiltonian.h"
#include "PXPHamiltonian.h"
#include "HeisenbergAncillaryNoiseHamiltonian.h"
int d = 2; // Spin 1/2

using namespace shrt;

int main(int argc, const char *argv[])
{
    // Read input arguments
    int cntr = 0;
    const char *infile = argv[++cntr];
    // Recover the necessary parameters from a Properties file
    Properties props(infile);
    if (argc > 2)
    {
        cntr++;
        cout << "Some properties may be now replaced by command line arguments"
            << endl;
        props.loadProperties(argc - cntr, &argv[cntr]);
    }

    // Model parameters
    int L = props.getIntProperty("L");		   // Length of system
    string model = props.getProperty("model"); 
    // Name of model. Support Ising, Heisenberg and PXP
    double J_ = props.getDoubleProperty("J");  
    double B_ = props.getDoubleProperty("B");
    double Jx_ = props.getDoubleProperty("Jx");
    double Jy_ = props.getDoubleProperty("Jy");
    double g_ = props.getDoubleProperty("g");
    double h_ = props.getDoubleProperty("h");

    // Chebyshev parameters
    int M = props.getIntProperty("M");				 // LARGEST ORDER in Chebyshev exp.
    double delta = props.getDoubleProperty("delta"); // Delta in Chebyshev exp.
    if (delta < 0 || delta > 1)
        delta = 0.05;    
    int saveTn = props.getIntProperty("saveTn");

    // MPS parameters
    double tol = props.getDoubleProperty("tol");
    int D = props.getIntProperty("D");
    int numD = props.getIntProperty("numD");

    // Output file names
    string outfname = props.getProperty("outputfile");
    string scalingfname = props.getProperty("scalingfile"); // Storing Emin and Emax


    cout << "Initialized arguments: L=" << L << ", D=" << D << ", M=" << M 
        << ", model=" << model << ", outfile=" << outfname <<endl;
    Contractor &contractor = Contractor::theContractor();
    contractor.setConvTol(tol);
    cout << "Initialized Contractor" << endl;
    cout << "First need to estimate the energy band to rescale "
        << "and shift the Hamiltonian" << endl;


    if (model == "HAN"){
        L *= 2;
    }

    // The original Hamiltonian to get the energy range
    MPO hamil0(L);
    if (model == "Ising"){
        IsingHamiltonian hamH0(L, d, J_, g_, h_);
        const MPO &tmphamil0 = hamH0.getHMPO();
        for (int k = 0; k<L; k++){
            hamil0.setOp(k, new Operator(tmphamil0.getOp(k).getFullData()), true);
        }
        cout << "J=" << J_ <<", g=" << g_ << ", h=" << h_ << endl;
    }
    else if (model == "Heisenberg")
    { 
        HeisenbergHamiltonian hamH0(L, Jx_, Jy_, 1.0, h_, d);
        const MPO &tmphamil0 = hamH0.getHMPO();
        for (int k = 0; k<L; k++){
            hamil0.setOp(k, new Operator(tmphamil0.getOp(k).getFullData()), true);
        }
        cout << "Jx=" << Jx_ <<", Jy=" << Jy_ << ", Jz=1.0, h=" << h_ << endl;

    }
    else if (model == "PXP"){
        PXPHamiltonian hamH0(L, 1., 0., 0, 0.);
        const MPO &tmphamil0 = hamH0.getHMPO();
        for (int k = 0; k<L; k++){
            hamil0.setOp(k, new Operator(tmphamil0.getOp(k).getFullData()), true);
        }
    }     
    else if (model == "HAN"){
        HeisenbergAncillaryNoiseHamiltonian hamH0(L/2, 1., 1., 0., J_, J_, 0., B_, 0.);
        const MPO &tmphamil0 = hamH0.getHMPO();
        for (int k = 0; k<L; k++){
            hamil0.setOp(k, new Operator(tmphamil0.getOp(k).getFullData()), true);
        }
        cout << "J=" << J_ <<", B=" << B_ << endl;

    }     
    cout << "Created the non-rescaled Hamiltonian" << endl;


    // First rescale the Hamiltonian:
    double scale = 1.;
    double offset = 0.;

    // Get ground and highest state energy, for rescaling
    double Emin = 0;
    MPS gs(L, D, d);
    gs.setRandomState();
    gs.gaugeCond('R');
    gs.gaugeCond('L', 1);
    contractor.findGroundState(hamil0, D, &Emin, gs);
    cout << "Found GS of H with E=" << Emin << endl;

    // Same for largest eigenvalue:
    double Emax = 0;
    MPS exc(L, D, d);
    exc.setRandomState();
    exc.gaugeCond('R');
    exc.gaugeCond('L', 1);
    {
        MPO hamilMinus(L);
        hamilMinus.setOp(
                0, new Operator(-1. * hamil0.getOp(0).getFullData()), true);
        for (int k = 1; k < L; k++)
        {
            hamilMinus.setOp(k, &hamil0.getOp(k), false);
        }
        contractor.findGroundState(hamilMinus, D, &Emax, exc);
        Emax = -Emax;
    }
    cout << "Found max Exc of H with E=" << Emax << endl;

    ofstream *fout;
    fout = new ofstream(scalingfname.data());
    *fout << Emin << "\t" << Emax << "\t" << endl;
    fout->close();
    delete fout;

    // Now rescale H such that the spec is in [-1+delta,1+delta]
    scale = 2. * (1. - delta) / (Emax - Emin);
    offset = -Emin * scale - 1. + delta;

    MPO hamil_(L);
    MPO proj(L);
    if (model == "Ising"){
        IsingHamiltonian hamH(L, d, J_*scale, g_*scale, h_*scale, offset);
        const MPO &tmphamil_ = hamH.getHMPO();
        for (int k = 0; k<L; k++){
            hamil_.setOp(k, new Operator(tmphamil_.getOp(k).getFullData()), true);
        }
    }
    else if (model == "Heisenberg")
    { 
        HeisenbergHamiltonian hamH(L, Jx_*scale, Jy_*scale, 
                1.0*scale, h_*scale, d, offset);
        const MPO &tmphamil_ = hamH.getHMPO();
        for (int k = 0; k<L; k++){
            hamil_.setOp(k, new Operator(tmphamil_.getOp(k).getFullData()), true);
        }
    }
    else if (model == "PXP"){
        PXPHamiltonian hamH(L, scale, 0., 0, offset);
        const MPO &tmphamil_ = hamH.getHMPO();
        for (int k = 0; k<L; k++){
            hamil_.setOp(k, new Operator(tmphamil_.getOp(k).getFullData()), true);
        }
        const MPO &tmpproj = hamH.getProjectorConstr();
        for (int k = 0; k<L; k++){
            proj.setOp(k, new Operator(tmpproj.getOp(k).getFullData()), true);
        }
    }
    else if (model == "HAN"){
        HeisenbergAncillaryNoiseHamiltonian hamH(L/2, scale, scale, 0., J_*scale, J_*scale, 0., B_*scale, offset);
        const MPO &tmphamil_ = hamH.getHMPO();
        for (int k = 0; k<L; k++){
            hamil_.setOp(k, new Operator(tmphamil_.getOp(k).getFullData()), true);
        }
    }
    cout << "Created the rescaled Hamiltonian: scale=" << scale
        << ", offset=" << offset << endl;



    // Prepare MPO and MPS needed for Chebyshev terms
    int Dh = 0;
    if (model == "Ising") {
        Dh = 3;
    } else if (model == "Heisenberg") {
        Dh = 5;
    } else if (model == "PXP") {
        Dh = 4;
    } else if (model == "HAN") {
        Dh = 8;
    }

    MPS idMPS(L, 1, d * d);
    for (int k = 0; k < L; k++)
        idMPS.setA(k, reshape(identityMatrix(d), Indices(d * d, 1, 1)));

    MPS hamilMPS(L, Dh, d * d);
    MPSfromMPO(hamil_, hamilMPS);
    MPS hamil0MPS(L, Dh, d * d);
    MPSfromMPO(hamil0, hamil0MPS);
    MPO hamilId(L); // needed for the recursion
    extendMPO(hamil_, hamilId, d);
    MPS projMPS(L, 2, d*d);
    MPO projId(L);
    if (model == "PXP"){
        MPSfromMPO(proj, projMPS);
        extendMPO(proj, projId, d);
    }


    // Prepare observables
    cout << "Prepare observables!" << endl;
    vector<const MPS *> opList;
    vector<string> opNameList;
    int numOp = 0;

    complex_t SigId[] = {ONE_c, ZERO_c, ZERO_c, ONE_c};
    complex_t SigX[] = {ZERO_c, ONE_c, ONE_c, ZERO_c};
    complex_t SigY[] = {ZERO_c, I_c, -1. * I_c, ZERO_c};
    complex_t SigZ[] = {ONE_c, ZERO_c, ZERO_c, -1. * ONE_c};
    vector<mwArray> PauliMatrices;
    mwArray *SigIdM = new mwArray(Indices(d * d, 1, 1), SigId);
    PauliMatrices.push_back(*SigIdM);
    mwArray *SigXM = new mwArray(Indices(d * d, 1, 1), SigX);
    PauliMatrices.push_back(*SigXM);
    mwArray *SigYM = new mwArray(Indices(d * d, 1, 1), SigY);
    PauliMatrices.push_back(*SigYM);
    mwArray *SigZM = new mwArray(Indices(d * d, 1, 1), SigZ);
    PauliMatrices.push_back(*SigZM);

    // DoS
    opList.push_back(&idMPS);
    opNameList.push_back("_DoS");
    numOp ++;

    // Hamiltonian
    opList.push_back(&hamil0MPS);
    opNameList.push_back("_H");
    numOp ++;

    // Sx,Sy,Sz in the middle
    MPS *Sx = new MPS(L, 1, d*d);
    MPS *Sy = new MPS(L, 1, d*d);
    MPS *Sz = new MPS(L, 1, d*d);
    MPS *Szp = new MPS(L, 1, d*d);
    for (int k = 0; k < L; k++){
        if (k == L / 2){
            Sx->setA(k, reshape((PauliMatrices[1]), Indices(d * d, 1, 1)));
            Sy->setA(k, reshape((PauliMatrices[2]), Indices(d * d, 1, 1)));
            Sz->setA(k, reshape((PauliMatrices[3]), Indices(d * d, 1, 1)));
        }
        else{
            Sx->setA(k, reshape((PauliMatrices[0]), Indices(d * d, 1, 1)));
            Sy->setA(k, reshape((PauliMatrices[0]), Indices(d * d, 1, 1)));
            Sz->setA(k, reshape((PauliMatrices[0]), Indices(d * d, 1, 1)));
        }
        if (k == L / 2+1){
            Szp->setA(k, reshape((PauliMatrices[3]), Indices(d * d, 1, 1)));
        }
        else{
            Szp->setA(k, reshape((PauliMatrices[0]), Indices(d * d, 1, 1)));
        }
    }
    opList.push_back(Sx);
    opList.push_back(Sy);
    opList.push_back(Sz);
    opList.push_back(Szp);

    opNameList.push_back("_SxMid");
    opNameList.push_back("_SyMid");
    opNameList.push_back("_SzMid");
    opNameList.push_back("_SzMidp");

    numOp += 4;

    // Sum Sz, Sum Sz Even, Sum Sz Odd
    SpinMPO spinMPO;
    MPO Sz2(L);
    spinMPO.getSz2MPO(L, 2, Sz2);
    MPO Sz2Even(L);
    spinMPO.getSz2EvenMPO(L, 2, Sz2Even);
    MPO Sz2Odd(L);
    spinMPO.getSz2OddMPO(L, 2, Sz2Odd);

    MPS *Sz2MPS = new MPS(L, 1, d*d);
    MPS *Sz2EvenMPS = new MPS(L, 1, d*d);
    MPS *Sz2OddMPS = new MPS(L, 1, d*d);
    MPSfromMPO(Sz2, *Sz2MPS);
    MPSfromMPO(Sz2Even, *Sz2EvenMPS);
    MPSfromMPO(Sz2Odd, *Sz2OddMPS);

    opList.push_back(Sz2MPS);
    opList.push_back(Sz2EvenMPS);
    opList.push_back(Sz2OddMPS);
    opNameList.push_back("_Sz2");
    opNameList.push_back("_Sz2Even");
    opNameList.push_back("_Sz2Odd");

    numOp += 3;



    // Clear the previous output files of ops
    for (int num = 0; num < numOp; num++){
        string tmpfname = outfname + opNameList[num];
        fout = new ofstream(tmpfname.data());
        fout->close();
        delete fout;
    }



    // Prepare first two Chebyshev polynomials
    MPS *chebyT_Nm1 = new MPS(idMPS);
    MPS *chebyT_N = new MPS(idMPS);
    chebyT_Nm1->gaugeCond('R', 0);
    chebyT_Nm1->gaugeCond('L', 0);
    chebyT_N->gaugeCond('R', 0);
    chebyT_N->gaugeCond('L', 0);
    double norm_Nm1 = log2(chebyT_Nm1->getNormFact());
    chebyT_Nm1->setNormFact(1.);
    double norm_N = log2(chebyT_N->getNormFact());
    chebyT_N->setNormFact(1.);


    // Truncate bond dimension of Tn to Ds to get truncation error
    fout = new ofstream((outfname+"_TnErr").data());	
    *fout << setprecision(10) << fixed;
    int Ds[numD];
    for (int k = 0; k < numD; k++){
        Ds[k] = D / (numD+1) * (k+1);
    }


    cout << "Start Chebyshev polynomials calculation" << endl;
    for (int k = 0; k < M; k++)
    {
        MPS *chebyT_Np1;
        double norm_Np1(0.);
        if (k == 0)
        {
            chebyT_Np1 = chebyT_Nm1;
            norm_Np1 = norm_Nm1;
        }
        else if (k == 1)
        {  
            MPS aux(*chebyT_Nm1);
            if (model == "PXP"){
                contractor.optimize(projId, *chebyT_N, aux, D);
            }
            contractor.optimize(hamilId, aux, *chebyT_N, D);
            if (model == "PXP"){
                contractor.optimize(projId, *chebyT_N, aux, D);
                delete chebyT_N;
                chebyT_N = new MPS(aux);
            }
            chebyT_N->gaugeCond('R', 0);
            chebyT_N->gaugeCond('L', 0);
            norm_N += log2(chebyT_N->getNormFact());
            chebyT_N->setNormFact(1.);     
            chebyT_Np1 = chebyT_N;
            norm_Np1 = norm_N;
        }
        else
        {
            chebyT_Np1 = new MPS(*chebyT_N);
            // TODO: Optimize the sum of ops acting on vecs!
            MPS aux1(*chebyT_N);
            contractor.optimize(hamilId, *chebyT_N, aux1, D);
            aux1.gaugeCond('R', 0);
            aux1.gaugeCond('L', 0);
            double auxN = log2(aux1.getNormFact());
            aux1.setNormFact(1.);
            MPS aux(aux1);

            if (model == "PXP"){
                cout << "Before apply proj" << endl;
                contractor.optimize(projId, aux1, aux, D);
                aux.gaugeCond('R', 0);
                aux.gaugeCond('L', 0);
                auxN += log2(aux.getNormFact());
                aux.setNormFact(1.);
                cout << "After apply proj" << endl;
            }
            vector<const MPS *> kets;
            kets.push_back(&aux);
            kets.push_back(chebyT_Nm1);
            vector<complex_t> coefs;
            coefs.push_back(2. * ONE_c);
            coefs.push_back(-pow(2, norm_Nm1 - norm_N - auxN) * ONE_c);
            contractor.optimizeSum(kets, coefs, *chebyT_Np1, D);
            chebyT_Np1->gaugeCond('R', 0);
            chebyT_Np1->gaugeCond('L', 0);
            norm_Np1 = log2(chebyT_Np1->getNormFact()) + norm_N + auxN;
            chebyT_Np1->setNormFact(1.);

            // Shift the pointers and delete the discarded one
            delete chebyT_Nm1;
            chebyT_Nm1 = chebyT_N;
            norm_Nm1 = norm_N;
            chebyT_N = chebyT_Np1;
            norm_N = norm_Np1;
        }

        // Store truncating error
        *fout << k << "\t";
        for (int num = 0; num < numD; num++){
            MPS tmpD(L, Ds[num], d*d);
            double err =1.;
            contractor.optimizeMPS(*chebyT_Np1, tmpD, Ds[num], &err);
            //cout << D << " " << Ds[num] << " " << tmpD.getBond() << endl;
            *fout << err << "\t";
        }
        *fout << endl;

        if (saveTn){
            chebyT_Np1->exportMPS((outfname+"_T"+SSTR(k)).data());
        }

        cout << "T_" << k << " has been calculated!" << endl;


        for (int num = 0; num < numOp; num++){
            complex_t Ck = contractor.contract(*(opList[num]), *chebyT_Np1);
            ofstream *foutOp;
            string tmpfname = outfname + opNameList[num];
            foutOp = new ofstream((tmpfname).data(), ios::app);
            *foutOp << setprecision(10) << fixed;
            *foutOp << k << '\t' << norm_Np1 << '\t' << Ck.re << endl;
            foutOp->close();
            delete foutOp;
        }

    }
    delete chebyT_Nm1;
    delete chebyT_N;
    fout->close();
    delete fout;
}
