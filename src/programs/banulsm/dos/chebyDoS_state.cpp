/* Calculate the local density of states. */
#include <iomanip>
#include <math.h>
#include <sstream>
#include <string>

#include "Contractor.h"
#include "MPS.h"
#include "Properties.h"
#include "misc.h"
#include "mwArray.h"

#define _USE_MATH_DEFINES // For PI
#define TMPSTORE 0		  // Using temporary storage for intermediate results
#define SSTR(x)                                                                \
    static_cast<std::ostringstream &>((std::ostringstream() << std::dec << x)) \
.str() // Convert variables to strings

// Included models
#include "HeisenbergHamiltonian.h"
#include "IsingHamiltonian.h"
#include "PXPHamiltonian.h"
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

    int L = props.getIntProperty("L");		   // Length of system
    string model = props.getProperty("model"); // Name of model
    // Ising: H = J_Σs_i^z*s_{i+1}^z + g_Σs_i^x + h_Σs_i^z
    // Heisenberg with disorder: H = J_Σs_i*s_{i+1} + h_Σs_i^z
    double J_ = props.getDoubleProperty("J");
    double Jy_ = props.getDoubleProperty("Jx");
    double Jx_ = props.getDoubleProperty("Jy");
    double g_ = props.getDoubleProperty("g");
    double h_ = props.getDoubleProperty("h");
    // Chebyshev parameters
    int M = props.getIntProperty("M");				 // LARGEST ORDER in Chebyshev exp.
    double delta = props.getDoubleProperty("delta"); // Delta in Chebyshev exp.
    if (delta < 0 || delta > 1)
        delta = 0.02;
    double Emin = props.getDoubleProperty("Emin");
    double Emax = props.getDoubleProperty("Emax");
    // MPS parameters
    double tol = props.getDoubleProperty("tol");
    int D = props.getIntProperty("D");
    // Output file names
    string outfname = props.getProperty("outputfile");
    string scalingfname =
        props.getProperty("scalingfile"); // Storing Emin and Emax

    cout << "Initialized arguments: L=" << L << ", J=" << J_ << ", g=" << g_
        << ", h=" << h_ << ", outfile=" << outfname << ", D=" << D << endl;

    Contractor &contractor = Contractor::theContractor();
    contractor.setConvTol(tol);
    cout << "Initialized Contractor" << endl;

    // The original Hamiltonian to get the energy
    MPO hamil0(L);
    if (model == "Ising"){
        IsingHamiltonian hamH0(L, d, J_, g_, h_);
        const MPO &tmphamil0 = hamH0.getHMPO();
        for (int k = 0; k<L; k++){
            hamil0.setOp(k, new Operator(tmphamil0.getOp(k).getFullData()), true);
        }
    }
    else if (model == "Heisenberg")
    { 
        HeisenbergHamiltonian hamH0(L, J_, J_, J_, h_, d);
        const MPO &tmphamil0 = hamH0.getHMPO();
        for (int k = 0; k<L; k++){
            hamil0.setOp(k, new Operator(tmphamil0.getOp(k).getFullData()), true);
        }
    }
    else if (model == "PXP"){
        PXPHamiltonian hamH0(L, 1., 0., 0, 0.);
        const MPO &tmphamil0 = hamH0.getHMPO();
        for (int k = 0; k<L; k++){
            hamil0.setOp(k, new Operator(tmphamil0.getOp(k).getFullData()), true);
        }
    }     

    double scale = 1.;
    double offset = 0.;
    if (Emin == -1 && Emax == -1)
    {
        cout << "First need to estimate the energy band to rescale "
            << "and shift the Hamiltonian" << endl;

        // First rescale the Hamiltonian:
        cout << "Created the non-rescaled Hamiltonian" << endl;

        // Get ground and highest state energy, for rescaling
        Emin = 0;
        MPS gs(L, D, d);
        gs.setRandomState();
        gs.gaugeCond('R');
        gs.gaugeCond('L', 1);
        contractor.findGroundState(hamil0, D, &Emin, gs);
        cout << "Found GS of H with E=" << Emin << endl;

        // Same for largest eigenvalue:
        Emax = 0;
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
    }

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



    cout << "Created the rescaled Hamiltonian: scale=" << scale
        << ", offset=" << offset << endl;

    // bond dimension of H, just in case
    int Dh = 0;
    if (model == "Ising") {
        Dh = 3;
    } else if (model == "Heisenberg") {
        Dh = 5;
    } else if (model == "PXP") {
        Dh = 4;
    }

    MPS idMPS(L, 1, d * d);
    for (int k = 0; k < L; k++)
        idMPS.setA(k, reshape(identityMatrix(d), Indices(d * d, 1, 1)));

    // Prepare the vectors to compute local density of states
    double sqrt2 = sqrt(2);
    complex_t xUp[] = {ONE_c / sqrt2, ONE_c / sqrt2};
    mwArray xUpVec(Indices(d, 1, 1), xUp);
    complex_t yUp[] = {ONE_c / sqrt2, I_c / sqrt2};
    mwArray yUpVec(Indices(d, 1, 1), yUp);
    complex_t zUp[] = {ONE_c, ZERO_c};
    mwArray zUpVec(Indices(d, 1, 1), zUp);
    complex_t zDn[] = {ZERO_c, ONE_c};
    mwArray zDnVec(Indices(d, 1, 1), zDn);

    vector<const MPS *> vecs;
    vector<string> names;

    MPS *xPlus = new MPS(L, 1, d);
    for (int k = 0; k < L; k++)
        xPlus->setA(k, xUpVec);
    vecs.push_back(xPlus);
    names.push_back("_Xp");

    MPS *yPlus = new MPS(L, 1, d);
    for (int k = 0; k < L; k++)
        yPlus->setA(k, yUpVec);
    vecs.push_back(yPlus);
    names.push_back("_Yp");

    MPS *zPlus = new MPS(L, 1, d);
    for (int k = 0; k < L; k++)
        zPlus->setA(k, zUpVec);
    vecs.push_back(zPlus);
    names.push_back("_Zp");

    int numRound = 10;
    for (int u = 0; u < numRound; u++){
        complex_t xToY[] = {ONE_c/sqrt2, ONE_c/sqrt2 * cos(M_PI/2/(numRound+1)*(u+1)) + 
            I_c/sqrt2 * sin(M_PI/2/(numRound+1)*(u+1)) };
        mwArray xToYVec(Indices(d,1,1), xToY);
        complex_t zToY[] = {ONE_c*cos(M_PI/4/(numRound+1)*(u+1)), 
            I_c*sin(M_PI/4/(numRound+1)*(u+1))};
        mwArray zToYVec(Indices(d,1,1), zToY);
        MPS *xY = new MPS(L, 1, d);
        for (int k = 0 ; k < L; k++)
            xY->setA(k, xToYVec);
        vecs.push_back(xY);
        names.push_back("_xtoY"+SSTR(u));
        MPS *zY = new MPS(L, 1, d);
        for (int k = 0 ; k < L; k++)
            zY->setA(k, zToYVec);
        vecs.push_back(zY);
        names.push_back("_ztoY"+SSTR(u));
    }

    MPS *zNeel = new MPS(L, 1, d);
    for (int k = 0; k < L; k++)
        if (k % 2 == 0)
        {
            zNeel->setA(k, zDnVec);
        }
        else
        {
            zNeel->setA(k, zUpVec);
        }
    vecs.push_back(zNeel);
    names.push_back("_Neel");



    for (int nvec = 0; nvec < vecs.size(); nvec++){
        MPS *chebyT_Nm1 = new MPS(*(vecs[nvec]));
        MPS *chebyT_N = new MPS(*chebyT_Nm1);
        MPS auxMPS(*chebyT_Nm1);
        if (model == "PXP"){
            contractor.optimize(proj, *chebyT_Nm1, auxMPS, D);
        }
        contractor.optimize(hamil_, auxMPS, *chebyT_N, D);
        if (model == "PXP"){
            contractor.optimize(proj, *chebyT_N, auxMPS, D);
            delete chebyT_N;
            chebyT_N = new MPS(auxMPS);
        }

        chebyT_Nm1->gaugeCond('R', 0);
        chebyT_Nm1->gaugeCond('L', 0);
        chebyT_N->gaugeCond('R', 0);
        chebyT_N->gaugeCond('L', 0);
        double norm_Nm1 = log2(chebyT_Nm1->getNormFact());
        chebyT_Nm1->setNormFact(1.);
        double norm_N = log2(chebyT_N->getNormFact());
        chebyT_N->setNormFact(1.);

        cout << setprecision(10) << fixed;
        fout = new ofstream((outfname + names[nvec]).data());
        *fout << setprecision(10) << fixed;

        // iteration
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
                chebyT_Np1 = chebyT_N;
                norm_Np1 = norm_N;
            }
            else
            {
                chebyT_Np1 = new MPS(*chebyT_N);
                MPS aux1(*chebyT_N);
                contractor.optimize(hamil_, *chebyT_N, aux1, D);
                aux1.gaugeCond('R', 0);
                aux1.gaugeCond('L', 0);
                double auxN = log2(aux1.getNormFact());
                aux1.setNormFact(1.);
                MPS aux(aux1);
                if (model == "PXP"){
                    contractor.optimize(proj, aux1, aux, D);
                    aux.gaugeCond('R', 0);
                    aux.gaugeCond('L', 0);
                    auxN += log2(aux.getNormFact());
                    aux.setNormFact(1.);
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

            *fout << k << '\t' << norm_Np1 << '\t';
            complex_t Ck = contractor.contract(*(vecs[nvec]), *chebyT_Np1);
            *fout << Ck.re << endl;
            cout << "Computed approximation for k=" << k << endl;
        }

        delete chebyT_Nm1;
        delete chebyT_N;
        fout->close();
        delete fout;
    }
}
