#define CATCH_CONFIG_MAIN

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision


#ifdef USING_PRIMME
extern "C" {
#include "primme.h"
}
#endif

#include "Indices.h"
#include "mwArray.h"
#include "catch.hpp"

using namespace std;
using namespace shrt;

bool areEqual(mwArray A, mwArray B, double eps) {
    if ((A.getRank() != B.getRank())
        || (A.getNrComponents() != B.getNrComponents())
        || (A.getDimensions() != B.getDimensions())
            ) {
        return false;
    }
    mwArray C = B - A;
    for (int i = 0; i < 2 * C.getNrComponents(); ++i) {
        if (abs(C.getComponents()[i]) > eps) {
            cout << i <<  "," << C.getComponents()[i] << endl;
            return false;
        };
    }

    return true;
}
/**
 * TODO: HConjugate, resize
 */



TEST_CASE("Creating an empty mwArray") {
    mwArray arr = mwArray();
    REQUIRE(arr.getRank() == 0);
    REQUIRE(arr.getNrComponents() == 0);
}

TEST_CASE("Basic mwArray Manipulations") {
    Indices dims(2, 4, 3);
    mwArray A2(dims), A3, A4, A5;
    double testDataDouble[dims[0] * dims[1] * dims[2] * 2];
    complex_t testDataComplex[dims[0] * dims[1] * dims[2]];
    for (int i = 0; i < dims[2]; i++) {
        for (int j = 0; j < dims[1]; j++) {
            for (int k = 0; k < dims[0]; k++) {
                testDataDouble[2 * (dims[0] * dims[1] * i + dims[0] * j + k)] = dims[0] * dims[1] * i + dims[0] * j + k;
                testDataDouble[2 * (dims[0] * dims[1] * i + dims[0] * j + k) + 1] =
                        dims[0] * dims[1] * i + dims[0] * j + k;
                testDataComplex[dims[0] * dims[1] * i + dims[0] * j + k] = {
                        (double) (dims[0] * dims[1] * i + dims[0] * j + k),
                        (double) (dims[0] * dims[1] * i + dims[0] * j + k)};

            }
        }
    }

    SECTION("Dimensions test") {
        REQUIRE(A2.getDimensions() == Indices(2, 4, 3));
    }


    SECTION("Filling mwArray with random values between 0 and 1") {
        A2.fillRandom();
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 3; k++) {
                    Indices ind(i, j, k);
                    REQUIRE(A2.getElement(ind).im <= 1);
                    REQUIRE(A2.getElement(ind).re <= 1);
                }
            }
        }
    }
    SECTION("== Operator") {
        A3 = mwArray(dims, testDataComplex);
        A4 = mwArray(dims, testDataComplex);
        REQUIRE(A3 == A4);
    }
    SECTION("Creation methods") {

        A3 = mwArray(dims, testDataDouble, false);
        for (int i = 0; i < dims[2]; i++) {
            for (int k = 0; k < dims[1]; k++) {
                for (int j = 0; j < dims[0]; j++) {
                    Indices ind(j, k, i);
                    REQUIRE(A3.getElement(ind).im == dims[0] * dims[1] * i + dims[0] * k + j);
                    REQUIRE(A3.getElement(ind).re == dims[0] * dims[1] * i + dims[0] * k + j);
                }
            }
        }
        A4 = mwArray(dims, testDataComplex);


        REQUIRE(A3 == A4);
        A5 = mwArray(A3);
        REQUIRE(A5 == A3);
    }

    SECTION("IO functions for mwArray") {
        char filename[] = "data.dat";

        /**
         * Write
         */

        ofstream outfile(filename, ios::binary);
        if (!outfile.is_open()) {
            cout << "Error: unable to open file to write" << filename << endl;
        } else {
            A3 = mwArray(dims, testDataDouble, false);
            A3.save(outfile);
            outfile.close();

        }

        /**
         * Read
         */

        ifstream infile(filename, ios::binary);
        if (!infile.is_open()) {
            cout << "Error: unable to open file " << filename << endl;
        } else {
            A4 = mwArray(infile);
            infile.close();
        }

        REQUIRE(A3 == A4);
    }
    SECTION("Subarray and permutations") {

        Indices subIdx(-1, 2, -1);
        A3 = mwArray(dims, testDataComplex);
        mwArray subA3 = A3.subArray(subIdx);
        complex_t subData[] = {
                {4.,  4.},
                {5.,  5.},
                {12., 12.},
                {13., 13.},
                {20., 20.},
                {21., 21.},

        };

        REQUIRE(subA3 == mwArray(Indices(2, 3), subData));
        A4 = mwArray(A3);
        Indices permut(2, 3, 1);
        A4.permute(permut);
        for (int k = 0; k < dims[1]; k++) {
            for (int j = 0; j < dims[0]; j++) {
                for (int i = 0; i < dims[2]; i++) {
                    Indices ind(k, i, j);
                    REQUIRE(A4.getElement(ind).im == dims[0] * dims[1] * i + dims[0] * k + j);
                    REQUIRE(A4.getElement(ind).re == dims[0] * dims[1] * i + dims[0] * k + j);
                }
            }
        }
    }
}

TEST_CASE("Basic matrix operations and per element amipulations") {
    Indices dim(2, 2);
    complex_t dataA1[] = {
            {1., 1.},
            {3., 3.},
            {2., 2.},
            {4., 4.}

    }, dataA1_1[] = {
            {1., 0.},
            {3., 3.},
            {2., 2.},
            {4., 4.}

    }, dataA1_2[] = {
            {0., 1.},
            {3., 3.},
            {2., 2.},
            {4., 4.}

    }, dataA2[] = {
            {2., 2.},
            {6., 6.},
            {4., 4.},
            {8., 8.}

    }, dataA3[] = {
            {0, 28},
            {0, 60},
            {0, 40},
            {0, 88}

    }, dataA5[] = {
            {1,  0},
            {4,  0},
            {9,  0},
            {16, 0}

    }, dataA6[] = {
            {1, 0},
            {2, 0},
            {3, 0},
            {4, 0}

    }, dataA4[] = {
            ZERO_c, ZERO_c,
            ZERO_c, ZERO_c
    };
    mwArray A1(dim, dataA1);
    mwArray A1_1(dim, dataA1_1);
    mwArray A1_2(dim, dataA1_2);
    mwArray A2(dim, dataA2);
    mwArray A3(dim, dataA3);
    mwArray A4(dim, dataA4);
    mwArray A5(dim, dataA5);
    mwArray A6(dim, dataA6);


    SECTION("+ Operator") {
        REQUIRE(A1 + A1 == A2);
    }
    SECTION("- Operator") {
        REQUIRE(A2 - A1 == A1);
    }
    SECTION("* Operator") {
        REQUIRE((complex_t) {2, 0} * A1 == A2);
    }
    SECTION("Elementwise square root") {
        REQUIRE(sqrt(A5) == A6);
    }
    SECTION("Multiple Operators") {
        REQUIRE(A2 - (complex_t) {2, 0} * A1 == A4);
        REQUIRE(A1 * A2 == A3);
    }


    /**
     * Per element selection
     */


    Indices pos(0, 0);
    SECTION("Set element with real and imaginary part") {
        A1.setElement(1, 0, pos);
        REQUIRE(A1 == A1_1);
    }

    SECTION("Set element with complex number") {
        A1.setElement(I_c, pos);
        REQUIRE(A1 == A1_2);
    }

    SECTION("Get element") {
        REQUIRE(A1_1.getElement(pos) == ONE_c);
    }
}

TEST_CASE("Creating diagonal matrices") {
    double dataA1[] = {1., 2., 3.};
    complex_t dataA2[] = {
            {1, 0},
            {0, 0},
            {0, 0},
            {0, 0},
            {2, 0},
            {0, 0},
            {0, 0},
            {0, 0},
            {3, 0}
    };
    mwArray A1 = realDiag(3, dataA1);
    mwArray A2(Indices(3, 3), dataA2);

    SECTION("Real") {
        REQUIRE(A1 == A2);
    }
    double dataA3[] = {1., 1., 2., 2., 3., 3.};
    complex_t dataA4[] = {
            {1, 1},
            {0, 0},
            {0, 0},
            {0, 0},
            {2, 2},
            {0, 0},
            {0, 0},
            {0, 0},
            {3, 3}
    };
    mwArray A3 = diag(3, dataA3);
    mwArray A4(Indices(3, 3), dataA4);

    SECTION("Complex(double pairs)") {
        REQUIRE(A3 == A4);
    }

    vector<complex_t> dataA5(3, ONE_c);
    complex_t dataA6[] = {
            {1, 0},
            {0, 0},
            {0, 0},
            {0, 0},
            {1, 0},
            {0, 0},
            {0, 0},
            {0, 0},
            {1, 0}
    };

    mwArray A5 = diag(dataA5);
    mwArray A6(Indices(3, 3), dataA6);
    SECTION("Complex(vector)") {
        REQUIRE(A5 == A6);
    }

    SECTION("Checks if matrix is diagonal") {
        REQUIRE(A5.isDiagonal());
        REQUIRE(A3.isDiagonal());
        A3.setElement(1, 0, Indices(0, 1));
        REQUIRE(!A3.isDiagonal());
    }

    double re, im;
    A3.trace(re, im);
    SECTION("Trace") {
        REQUIRE(A3.trace() == (complex_t) {6, 6});
        REQUIRE(re == 6);
        REQUIRE(im == 6);
    }


    complex_t dataA7[] = {
            {1. / 2.,  -1. / 2.},
            {0,        0},
            {0,        0},
            {0,        0},
            {2. / 8.,  -2. / 8.},
            {0,        0},
            {0,        0},
            {0,        0},
            {3. / 18., -3. / 18.}
    };
    complex_t dataA8[] = {
            {1,       0},
            {0,       0},
            {0,       0},
            {0,       0},
            {1. / 2., 0},
            {0,       0},
            {0,       0},
            {0,       0},
            {1. / 3., 0}
    };
    int aux;
    mwArray A7(Indices(3, 3), dataA7);
    mwArray A8(Indices(3, 3), dataA8);
    SECTION("Inverse of diagonal matrix(real)") {
        REQUIRE(invertDiag(A2, aux) == A8);
    }
    SECTION("Inverse of diagonal matrix(complex)") {
        REQUIRE(invertDiag(A4, aux) == A7);

    }

}

TEST_CASE("SVD") {
    mwArray A(Indices(3, 4));
    A.fillRandom();
    mwArray S, V, U;
    wrapper::svd(A, U, S, V);

    /**
     * Checks if results are the same within uncertainty epsilon
     */
    bool flag = true;
    double epsilon = 1e-15;
    mwArray diff = A - U * S * V;
    for (int i = 0; i < diff.getNrComponents(); ++i) {
        flag &= diff.getComponents()[2 * i] < epsilon;
        flag &= diff.getComponents()[2 * i + 1] < epsilon;
    }
    REQUIRE(flag);
}

TEST_CASE("Some") {
    mwArray I1 = identityMatrix(4);
    mwArray I2 = diag(vector<complex_t>(4, ONE_c));
    REQUIRE(I1 == I2);


    complex_t dataVector[] = {
            {1, 0},
            {2, 0},
            {3, 0},
            {4, 0}
    };
    mwArray vectorR(Indices(1, 4), dataVector);
    mwArray vectorC(Indices(4, 1), dataVector);
    mwArray vector0((complex_t) {30, 0});
    REQUIRE(!I1.isVector());
    REQUIRE(vectorR.isVector());
    REQUIRE(vectorC.isVector());
    REQUIRE(!vector0.isVector());


    complex_t dataM[] = {
            {1, 0},
            {2, 0},
            {3, 0},
            {4, 0},
            {1, 0},
            {2, 0},
            {3, 0},
            {4, 0},
            {1, 0},
            {2, 0},
            {3, 0},
            {4, 0},
            {1, 0},
            {2, 0},
            {3, 0},
            {4, 0}
    };
    mwArray M(Indices(4, 4), dataM);


    complex_t res1[] = {
            {10, 0},
            {20, 0},
            {30, 0},
            {40, 0}
    };

    complex_t res2[] = {
            {30, 0},
            {30, 0},
            {30, 0},
            {30, 0}
    };
    complex_t res3[] = {
            {1,  0},
            {2,  0},
            {3,  0},
            {4,  0},
            {2,  0},
            {4,  0},
            {6,  0},
            {8,  0},
            {3,  0},
            {6,  0},
            {9,  0},
            {12, 0},
            {4,  0},
            {8,  0},
            {12, 0},
            {16, 0}
    };


    SECTION("Right") {
        REQUIRE(M * vectorC == mwArray(Indices(4, 1), res1));
    }
    SECTION("Left") {
        REQUIRE(vectorR * M == mwArray(Indices(1, 4), res2));
    }
    SECTION("Vector (inner) Vector") {
        REQUIRE(vectorR * vectorC == vector0);
    }
    SECTION("Vector (outer) Vector") {
        REQUIRE(vectorC * vectorR == mwArray(Indices(4, 4), res3));
        REQUIRE(outerproduct(vectorR, vectorC) == mwArray(Indices(4, 4), res3));
    }
    complex_t dataA1[] = {
            {1, 1},
            {2, 2},
            {3, 3},
            {4, 4},
            {5, 5},
            {6, 6}
    }, dataA2[] = {
            {1, -1},
            {4, -4},
            {2, -2},
            {5, -5},
            {3, -3},
            {6, -6}
    };
    mwArray A1(Indices(3, 2), dataA1);
    mwArray A2(Indices(2, 3), dataA2);
    SECTION("Hermitian Conjugate") {
        REQUIRE(Hconjugate(A1) == A2);
        REQUIRE(Hconjugate(mwArray((complex_t) {1, 1})) == mwArray((complex_t) {1, -1}));
    }
    complex_t dataA3[] = {
            {0, 0},
            {1, 0},
            {1, 0},
            {0, 0}
    };
    complex_t dataEigVRes[] = {
            {-sqrt(2) / 2, 0},
            {sqrt(2) / 2,  0},
            {sqrt(2) / 2,  0},
            {sqrt(2) / 2,  0}
    };
    mwArray A3(Indices(2, 2), dataA3);
    mwArray eigVRes(Indices(2, 2), dataEigVRes);
    mwArray eigVectors1, eigVectors2;
    vector<complex_t> eigValues1, eigValues2;
    bool isEqEig = true, isEqEigs = true;
    double eps = 1e-15;
    wrapper::eig(A3, eigValues1, eigVectors1, 1);
    mwArray diff = eigVRes - eigVectors1;
    for (int i = 0; i < diff.getNrComponents(); ++i) {
        isEqEig &= abs(diff.getComponents()[2 * i]) < eps;
        isEqEig &= abs(diff.getComponents()[2 * i + 1]) < eps;
    }
    SECTION("Eigenvalues and Eigenvectors with 'eig' function") {
        REQUIRE(eigValues1[0] == (complex_t) {-1, 0});
        REQUIRE(eigValues1[1] == (complex_t) {1, 0});
        REQUIRE(isEqEig);
    }

    mwArray herm2(Indices(4, 4));
    herm2.fillRandom();
    herm2 = herm2 + Hconjugate(herm2);

    complex_t dataX[] = {
            {1., 1.},
            {2., 2.},
            {3., 3.},
            {4., 4.}
    };

    /**
     * TODO :  uncertainties
     */

    mwArray X(Indices(4, 1), dataX); // actual sol
    mwArray B = herm2 * X; // rhs
    mwArray sol(Indices(4, 1));
    wrapper::lsd(herm2, B, sol);


    SECTION("Solution with lsd") {
        REQUIRE(areEqual(sol, X, 1e-12));
    }
    wrapper::lss(herm2, B, sol);
    SECTION("Solution with lss") {
        REQUIRE(areEqual(sol, X, 1e-12));

    }
    wrapper::lslu(herm2, B, sol);
    SECTION("Solution with lslu") {
        REQUIRE(areEqual(sol, X, 1e-12));
    }


    mwArray testMatrix(Indices(4, 3));
    testMatrix.fillRandom();
    mwArray L, Uu, P;
    wrapper::lu(testMatrix, L, Uu, P);
    SECTION("LU Decomposition") {
        REQUIRE(areEqual(P * L * Uu, testMatrix, 1e-14));
    }
    mwArray Q1, R1;
    wrapper::qr(testMatrix, Q1, R1);
    SECTION("QR Decomposition") {
        REQUIRE(areEqual(Q1 * R1, testMatrix, 1e-14));
    }
    mwArray L1;
    wrapper::lq(testMatrix, L1, Q1);
    SECTION("LQ Decomposition") {
        REQUIRE(areEqual(L1 * Q1, testMatrix, 1e-14));
    }


    complex_t dataA4[] = {
            {1,  0},
            {2,  0},
            {3,  0},
            {4,  0},
            {5,  0},
            {6,  0},
            {7,  0},
            {8,  0},

            {3,  0},
            {6,  0},
            {9,  0},
            {12, 0},
            {2,  0},
            {4,  0},
            {6,  0},
            {8,  0},

            {3,  0},
            {6,  0},
            {9,  0},
            {12, 0},
            {4,  0},
            {8,  0},
            {12, 0},
            {16, 0}
    };
    mwArray VVT = outerproduct(vectorR, vectorC);
    mwArray A4(Indices(2, 4, 3), dataA4);

    SECTION("Shallow copy test") {

        A4.permute(Indices(2, 3, 1));
        mwArray Ashallow;
        Ashallow.setPointer(A4.getDimensions(), (double *) A4.getComponents());
        mwArray A80 = permute(A4, Indices(3, 1, 2));
        Ashallow.permute(Indices(3, 1, 2));
        REQUIRE(areEqual(A80, reshape(A4, A80.getDimensions()), 1e-15));

    }


    SECTION("Test new trasposition!") {
        mwArray A81(Indices(2, 5));
        A81.fillRandom();
        mwArray A82 = Hconjugate(A81);
        mwArray A83 = reshape(A81, Indices(2, 1, 5));
        A83.permute(Indices(3, 2, 1), true);
        A83.reshape(Indices(5, 2));
        REQUIRE(areEqual(A83, A82, 1e-15));
    }


    SECTION("Test expm") {
        complex_t dataLn[] = {
                {0, 0},
                {1, 0},
                {1, 0},
                {0, 0}
        }, dataExp[] = {
                {exp(1) + exp(-1), 0},
                {exp(1) - exp(-1), 0},
                {exp(1) - exp(-1), 0},
                {exp(1) + exp(-1), 0},
        };
        mwArray Ln(Indices(2, 2), dataLn);
        mwArray Exp(Indices(2, 2), dataExp);
        Exp = 0.5 * Exp;
        mwArray expH;

        wrapper::expm(Ln, expH);
        REQUIRE(areEqual(Exp, expH, 1e-15));
    }
    SECTION("Test inverse") {
        mwArray A90(Indices(5, 5));
        A90.fillRandom();
        int nrVals;
        mwArray A91 = inverse(A90, nrVals);
        REQUIRE(areEqual(A90 * A91, identityMatrix(A90.getDimension(0)), 1e-12));
        REQUIRE(areEqual(A91 * A90, identityMatrix(A90.getDimension(0)), 1e-12));
    }
    int dimC = 5;
    int d1(1), d2(4), d3(5), d4(5), d5(1), d6(4);
    Indices dims1(d1, d2, d3), order1(1, 4, 2, 3);
    int dim1L = d1 * d2 * d3;
    Indices dims2(d4, d5, d6), order2(2, 3, 1, 4);
    int dim2R = d4 * d5 * d6;

    SECTION("Test single leg contraction") {
        mwArray A92(Indices(dims1, dimC));
        A92.fillRandom();
        mwArray A93(Indices(dimC, dims2));
        A93.fillRandom();
        mwArray A94 = reshape(A92, Indices(dim1L, dimC)) * reshape(A93, Indices(dimC, dim2R));
        A94.reshape(Indices(dims1, dims2));
        // now reshuffle and contract
        A92.permute(Indices(order1));
        A93.permute(Indices(order2));
        mwArray A95 = contractLeg(A92, 2, A93, 3);
    }

    SECTION("Test multiple leg contraction") {
        Indices dimsC(2, 3, 2);
        int cTot = dimsC[0] * dimsC[1] * dimsC[2];
        mwArray A96(Indices(dim1L, cTot));
        A96.fillRandom();
        mwArray A97(Indices(cTot, dim2R));
        A97.fillRandom();
        mwArray A98 = A96 * A97;
        A98.reshape(Indices(dims1, dims2));
        A96.reshape(Indices(dims1, dimsC)); // ordered d1 d2 d3 dc1 dc2 dc3
        A97.reshape(Indices(dimsC, dims2)); // ordered dc1 dc2 dc3 d4 d5 d6
        order1 = Indices(1, 4, 5, 2, 6, 3);
        A96.permute(order1);
        order2 = Indices(4, 1, 5, 2, 6, 3);
        A97.permute(order2);
        mwArray A99 = contractLegs(A96, Indices(2, 3, 5), A97, Indices(2, 4, 6));
        REQUIRE(areEqual(A98, A99, 1e-12));

    }
}



/**
 * TODO:S
*/
/*   wrapper::eig(herm2,eigValues1,eigVectors1,1);
      wrapper::eigs(herm2,2,"SM",eigValues2,eigVectors2,1);

      SECTION("Eigenvalues and Eigenvectors with 'eigs' function"){
          bool areEigs1 = false;
          bool areEigs2 = false;
          double eps = 1e-15;

          for (int j = 0; j < 4; ++j) {
              if( abs((eigValues1[j] - eigValues2[0]).re) < eps
                      && abs((eigValues1[j] - eigValues2[0]).im) < eps) areEigs1 = true;
              if( abs((eigValues1[j] - eigValues2[1]).re) < eps
                      && abs((eigValues1[j] - eigValues2[1]).im) < eps) areEigs2 = true;
          }

          REQUIRE((areEigs1 && areEigs2));
          /** TODO: Eigenvectors ???
      }*/
/**
TODO: PRIMME
#ifdef USING_PRIMME
  int err=wrapper::eigs_primme(herm2,2,primme_smallest,eigValues2,eigVectors2,1);
  if(err!=0){
    cout<<"ERROR in the diagonalization!"<<endl;
    exit(err);
  }

SECTION("Eigenvalues and Eigenvectors with 'eigs' function(PRIMME)"){
    bool areEigs1 = false;
    bool areEigs2 = false;
    double eps = 1e-15;
    for (int j = 0; j < 4; ++j) {
        if( abs((eigValues1[j] - eigValues2[0]).re) < eps
                && abs((eigValues1[j] - eigValues2[0]).im) < eps) areEigs1 = true;
        if( abs((eigValues1[j] - eigValues2[1]).re) < eps
                && abs((eigValues1[j] - eigValues2[1]).im) < eps) areEigs2 = true;
    }

    REQUIRE((areEigs1 && areEigs2));
    //** TODO: Eigenvectors ???
}
#endif
*/
/**TODO:???
cout<<"A="<<VVT<<endl;
cout<<"X="<<X<<endl;
B=VVT*X; // rhs
wrapper::lsd(VVT,B,sol);
cout<<"After lsd, solution: "<<sol<<endl;
*/
// Test permutation identity
/*cout<<"Permute A4 "<<A4<<" as (2,1,3)"<<endl;
cout<<permute(A4,Indices(2,1,3))<<endl;*/
