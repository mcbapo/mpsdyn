#include "SpinMPO.h"

using namespace std;
using namespace shrt;

void SpinMPO::getSxMPO(int L, int spindim, MPO &Sx, bool staggered)
{
  if (spindim != 2)
  {
    cout << "ERROR: Only spin 1/2 supported by SpinMPO" << endl;
    exit(1);
  }
  cout << "SpinMPO::getSxMPO()" << endl;
  int d = spindim;
  mwArray sig0 = identityMatrix(d);
  complex_t datax[] = {ZERO_c, ONE_c, ONE_c, ZERO_c};
  mwArray sig1(Indices(d, d), datax);
  getSingleBodyMPO(L, sig0, .5 * sig1, Sx, staggered);
}

void SpinMPO::getSyMPO(int L, int spindim, MPO &Sy, bool staggered)
{
  if (spindim != 2)
  {
    cout << "ERROR: Only spin 1/2 supported by SpinMPO" << endl;
    exit(1);
  }
  // cout<<"SpinMPO::getSyMPO()"<<endl;
  int d = spindim;
  mwArray sig0 = identityMatrix(d);
  complex_t datay[] = {ZERO_c, I_c, -1. * I_c, ZERO_c};
  mwArray sig1(Indices(d, d), datay);
  getSingleBodyMPO(L, sig0, .5 * sig1, Sy, staggered);
}

void SpinMPO::getSzMPO(int L, int spindim, MPO &Sz, bool staggered)
{
  if (spindim != 2)
  {
    cout << "ERROR: Only spin 1/2 supported by SpinMPO" << endl;
    exit(1);
  }
  // cout<<"SpinMPO::getSzMPO()"<<endl;
  int d = spindim;
  mwArray sig0 = identityMatrix(d);
  complex_t dataz[] = {ONE_c, ZERO_c, ZERO_c, -1. * ONE_c};
  mwArray sig1(Indices(d, d), dataz);
  getSingleBodyMPO(L, sig0, .5 * sig1, Sz, staggered);
}

void SpinMPO::getSz2MPO(int L, int spindim, MPO &Sz2)
{
  if (spindim != 2)
  {
    cout << "ERROR: Only spin 1/2 supported by SpinMPO" << endl;
    exit(1);
  }
  // cout<<"SpinMPO::getSz2MPO()"<<endl;
  Sz2.initLength(L);
  int d = spindim;
  // The operators
  mwArray sig0 = identityMatrix(d);
  complex_t dataz[] = {ONE_c, ZERO_c, ZERO_c, -1. * ONE_c};
  mwArray sig1(Indices(d, d), dataz);
  mwArray Z(Indices(2, d, d));
  for (int i = 0; i < d; i++)
    for (int j = 0; j < d; j++)
    {
      Z.setElement(sig0.getElement(Indices(i, j)), Indices(0, i, j));
      Z.setElement(sig1.getElement(Indices(i, j)), Indices(1, i, j));
    }
  Z.reshape(Indices(2, d * d));
  // The coefficients (edges are different)
  int D = 3; // bond dimension
  mwArray C(Indices(D, D, 2));
  mwArray Cl(Indices(1, D, 2));              // left edge
  mwArray Cr(Indices(D, 1, 2));              // right edge
  C.setElement(ONE_c, Indices(0, D - 1, 0)); // identity terms
  C.setElement(ONE_c, Indices(0, 0, 0));     // before any operator
  Cl.setElement(ONE_c, Indices(0, D - 1, 0));
  Cl.setElement(ONE_c, Indices(0, 0, 0));
  Cr.setElement(ONE_c, Indices(0, 0, 0));
  C.setElement(ONE_c, Indices(1, 1, 0));               // intermediate
  C.setElement(sqrt(2) * ONE_c, Indices(0, 1, 1));     // left part of ZZ
  Cl.setElement(sqrt(2) * ONE_c, Indices(0, 1, 1));    // left part of ZZ
  C.setElement(sqrt(2) * ONE_c, Indices(1, D - 1, 1)); // right part of ZZ
  Cr.setElement(sqrt(2) * ONE_c, Indices(1, 0, 1));
  C.setElement(ONE_c, Indices(D - 1, D - 1, 0)); // after every operator
  Cr.setElement(ONE_c, Indices(D - 1, 0, 0));    // after every operator
  C.reshape(Indices(D * D, 2));
  Cl.reshape(Indices(D, 2));
  Cr.reshape(Indices(D, 2));
  C.multiplyRight(Z);
  Cl.multiplyRight(Z);
  Cr.multiplyRight(Z);
  C.reshape(Indices(D, D, d, d));
  C.permute(Indices(3, 1, 4, 2));
  Cl.reshape(Indices(1, D, d, d));
  Cl.permute(Indices(3, 1, 4, 2));
  Cr.reshape(Indices(D, 1, d, d));
  Cr.permute(Indices(3, 1, 4, 2));
  // Assign the operators to the MPO
  Sz2.setOp(0, new Operator(Cl), true);
  Sz2.setOp(L - 1, new Operator(Cr), true);
  if (L > 2)
    Sz2.setOp(1, new Operator(C), true);
  for (int l = 2; l < L - 1; l++)
    Sz2.setOp(l, &Sz2.getOp(1), false);
}

void SpinMPO::getSz2OddMPO(int L, int spindim, MPO &Sz2Odd)
{
  if (spindim != 2)
  {
    cout << "ERROR: Only spin 1/2 supported by SpinMPO" << endl;
    exit(1);
  }
  // cout<<"SpinMPO::getSz2MPO()"<<endl;
  Sz2Odd.initLength(L);
  int d = spindim;
  // The operators
  mwArray sig0 = identityMatrix(d);
  complex_t dataz[] = {ONE_c, ZERO_c, ZERO_c, -1. * ONE_c};
  mwArray sig3(Indices(d, d), dataz);
  mwArray Z(Indices(2, d, d));
  for (int i = 0; i < d; i++)
    for (int j = 0; j < d; j++)
    {
      Z.setElement(sig0.getElement(Indices(i, j)), Indices(0, i, j));
      Z.setElement(sig3.getElement(Indices(i, j)), Indices(1, i, j));
    }
  Z.reshape(Indices(2, d * d));
  // The coefficients (edges are different)
  int D = 3; // bond dimension
  mwArray C(Indices(D, D, 2));
  C.setElement(ONE_c, Indices(0, D - 1, 0));           // identity terms
  C.setElement(ONE_c, Indices(0, 0, 0));               // before any operator
  C.setElement(ONE_c, Indices(1, 1, 0));               // intermediate
  C.setElement(sqrt(2) * ONE_c, Indices(0, 1, 1));     // left part of ZZ=
  C.setElement(sqrt(2) * ONE_c, Indices(1, D - 1, 1)); // right part of ZZ
  C.setElement(ONE_c, Indices(D - 1, D - 1, 0));       // after every operator

  mwArray CId(Indices(D, D, 2));
  CId.setElement(ONE_c, Indices(0, 0, 0));
  CId.setElement(ONE_c, Indices(1, 1, 0));
  CId.setElement(ONE_c, Indices(2, 2, 0));

  mwArray Cl(Indices(1, D, 2)); // left edge
  Cl.setElement(ONE_c, Indices(0, 0, 0));

  mwArray Cr(Indices(D, 1, 2)); // right edge
  if (L % 2 == 0)
  {
    Cr.setElement(ONE_c, Indices(0, 0, 0));
    Cr.setElement(sqrt(2) * ONE_c, Indices(1, 0, 1));
    Cr.setElement(ONE_c, Indices(D - 1, 0, 0)); // after every operator
  }
  else
  {
    Cr.setElement(ONE_c, Indices(D - 1, 0, 0));
  }

  C.reshape(Indices(D * D, 2));
  CId.reshape(Indices(D * D, 2));
  Cl.reshape(Indices(D, 2));
  Cr.reshape(Indices(D, 2));

  C.multiplyRight(Z);
  CId.multiplyRight(Z);
  Cl.multiplyRight(Z);
  Cr.multiplyRight(Z);

  C.reshape(Indices(D, D, d, d));
  C.permute(Indices(3, 1, 4, 2));
  CId.reshape(Indices(D, D, d, d));
  CId.permute(Indices(3, 1, 4, 2));
  Cl.reshape(Indices(1, D, d, d));
  Cl.permute(Indices(3, 1, 4, 2));
  Cr.reshape(Indices(D, 1, d, d));
  Cr.permute(Indices(3, 1, 4, 2));

  // Assign the operators to the MPO
  Sz2Odd.setOp(0, new Operator(Cl), true);
  Sz2Odd.setOp(L - 1, new Operator(Cr), true);
  if (L > 2)
    Sz2Odd.setOp(1, new Operator(C), true);
  if (L > 3)
    Sz2Odd.setOp(2, new Operator(CId), true);

  for (int l = 3; l < L - 1; l++)
  {
    if (l % 2 == 1)
    {
      Sz2Odd.setOp(l, &Sz2Odd.getOp(1), false);
    }
    else
    {
      Sz2Odd.setOp(l, &Sz2Odd.getOp(2), false);
    }
  }
}

void SpinMPO::getSz2EvenMPO(int L, int spindim, MPO &Sz2Even)
{
  if (spindim != 2)
  {
    cout << "ERROR: Only spin 1/2 supported by SpinMPO" << endl;
    exit(1);
  }
  // cout<<"SpinMPO::getSz2MPO()"<<endl;
  Sz2Even.initLength(L);
  int d = spindim;
  // The operators
  mwArray sig0 = identityMatrix(d);
  complex_t dataz[] = {ONE_c, ZERO_c, ZERO_c, -1. * ONE_c};
  mwArray sig3(Indices(d, d), dataz);
  mwArray Z(Indices(2, d, d));
  for (int i = 0; i < d; i++)
    for (int j = 0; j < d; j++)
    {
      Z.setElement(sig0.getElement(Indices(i, j)), Indices(0, i, j));
      Z.setElement(sig3.getElement(Indices(i, j)), Indices(1, i, j));
    }
  Z.reshape(Indices(2, d * d));
  // The coefficients (edges are different)
  int D = 3; // bond dimension
  mwArray C(Indices(D, D, 2));
  C.setElement(ONE_c, Indices(0, D - 1, 0));           // identity terms
  C.setElement(ONE_c, Indices(0, 0, 0));               // before any operator
  C.setElement(ONE_c, Indices(1, 1, 0));               // intermediate
  C.setElement(sqrt(2) * ONE_c, Indices(0, 1, 1));     // left part of ZZ=
  C.setElement(sqrt(2) * ONE_c, Indices(1, D - 1, 1)); // right part of ZZ
  C.setElement(ONE_c, Indices(D - 1, D - 1, 0));       // after every operator

  mwArray CId(Indices(D, D, 2));
  CId.setElement(ONE_c, Indices(0, 0, 0));
  CId.setElement(ONE_c, Indices(1, 1, 0));
  CId.setElement(ONE_c, Indices(2, 2, 0));

  mwArray Cl(Indices(1, D, 2)); // left edge
  Cl.setElement(ONE_c, Indices(0, D - 1, 0));
  Cl.setElement(ONE_c, Indices(0, 0, 0));
  Cl.setElement(sqrt(2) * ONE_c, Indices(0, 1, 1)); // left part of ZZ

  mwArray Cr(Indices(D, 1, 2)); // right edge
  if (L % 2 == 1)
  {
    Cr.setElement(ONE_c, Indices(0, 0, 0));
    Cr.setElement(sqrt(2) * ONE_c, Indices(1, 0, 1));
    Cr.setElement(ONE_c, Indices(D - 1, 0, 0)); // after every operator
  }
  else
  {
    Cr.setElement(ONE_c, Indices(D - 1, 0, 0));
  }

  C.reshape(Indices(D * D, 2));
  CId.reshape(Indices(D * D, 2));
  Cl.reshape(Indices(D, 2));
  Cr.reshape(Indices(D, 2));

  C.multiplyRight(Z);
  CId.multiplyRight(Z);
  Cl.multiplyRight(Z);
  Cr.multiplyRight(Z);

  C.reshape(Indices(D, D, d, d));
  C.permute(Indices(3, 1, 4, 2));
  CId.reshape(Indices(D, D, d, d));
  CId.permute(Indices(3, 1, 4, 2));
  Cl.reshape(Indices(1, D, d, d));
  Cl.permute(Indices(3, 1, 4, 2));
  Cr.reshape(Indices(D, 1, d, d));
  Cr.permute(Indices(3, 1, 4, 2));

  // Assign the operators to the MPO
  Sz2Even.setOp(0, new Operator(Cl), true);
  Sz2Even.setOp(L - 1, new Operator(Cr), true);
  if (L > 2)
    Sz2Even.setOp(1, new Operator(CId), true);
  if (L > 3)
    Sz2Even.setOp(2, new Operator(C), true);

  for (int l = 3; l < L - 1; l++)
  {
    if (l % 2 == 1)
    {
      Sz2Even.setOp(l, &Sz2Even.getOp(1), false);
    }
    else
    {
      Sz2Even.setOp(l, &Sz2Even.getOp(2), false);
    }
  }
}
void SpinMPO::getCyclicTranslationMPO(int L, int dim, MPO &T)
{
  T.initLength(L);
  int d = dim;
  int D = d * d;
  // The basic operator on any middle site
  mwArray basicOp = identityMatrix(d * d * d);
  basicOp.reshape(Indices(d, d, d, d, d, d));
  basicOp.permute(Indices(5, 1, 2, 3, 4, 6));
  basicOp.reshape(Indices(d, d * d, d, d * d));
  // The operator on the left
  mwArray leftOp = identityMatrix(d * d);
  leftOp.reshape(Indices(d, 1, d, d * d));
  mwArray rightOp = identityMatrix(d * d);
  rightOp.reshape(Indices(d, d, d, d));
  rightOp.permute(Indices(4, 1, 2, 3));
  rightOp.reshape(Indices(d, d * d, d, 1));
  // Now set the operators in the chain
  T.setOp(0, new Operator(leftOp), true);
  T.setOp(L - 1, new Operator(rightOp), true);
  if (L > 2)
    T.setOp(1, new Operator(basicOp), true);
  for (int k = 2; k < L - 1; k++)
    T.setOp(k, &T.getOp(1), false);
}

void SpinMPO::getMomentumDistributionMPO(int L, int dim, int k, MPO &Nkmpo,
                                         bool prepared)
{
  if (dim != 2)
  {
    cout << "Wrong dimension for SpinMPO::getMomentumDistributionMPO, right now,"
         << " only dim=2 supported" << endl;
    exit(212);
  }
  // for the time being, prepare every time
  prepared = false; // TODO: change!
  if (!prepared)
  {
    Nkmpo.initLength(L);
    // the basic operators are identity, sigma+, sigma- and sigma+sigma-
    mwArray idOp = identityMatrix(dim);
    mwArray sigP(Indices(dim, dim));
    sigP.setElement(ONE_c, Indices(0, 1));
    mwArray sigM(sigP);
    sigM.Hconjugate();
    mwArray Nop(idOp);
    Nop.setElement(ZERO_c, Indices(1, 1));
    int D = 4; // bond dimension of the MPO
    // The operators
    mwArray Z(Indices(4, dim, dim));
    for (int i = 0; i < dim; i++)
      for (int j = 0; j < dim; j++)
      {
        Z.setElement(idOp.getElement(Indices(i, j)), Indices(0, i, j));
        Z.setElement(Nop.getElement(Indices(i, j)), Indices(3, i, j));
        Z.setElement(sigP.getElement(Indices(i, j)), Indices(1, i, j));
        Z.setElement(sigM.getElement(Indices(i, j)), Indices(2, i, j));
      }
    Z.reshape(Indices(4, dim * dim));

    // The coefficients
    complex_t kfac = exp(I_c * k * 2 * M_PIl / L);
    mwArray C(Indices(D, D, 4));
    mwArray Cl(Indices(1, D, 4)); // left edge
    mwArray Cr(Indices(D, 1, 4)); // right edge

    C.setElement(ONE_c, Indices(0, 0, 0)); // nothing yet
    Cl.setElement(ONE_c, Indices(0, 0, 0));
    C.setElement((1 / L) * ONE_c, Indices(0, D - 1, 3));    // local term
    Cl.setElement((1 / L) * ONE_c, Indices(0, D - 1, 3));   // local term
    Cr.setElement((1 / L) * ONE_c, Indices(0, 0, 3));       // local term
    C.setElement((1 / sqrt(L)) * ONE_c, Indices(0, 1, 1));  // sigP first part
    Cl.setElement((1 / sqrt(L)) * ONE_c, Indices(0, 1, 1)); // sigP first part

    C.setElement((1 / sqrt(L)) * kfac, Indices(1, D - 1, 2)); // sigM second part
    Cr.setElement((1 / sqrt(L)) * kfac, Indices(1, 0, 2));    // sigM second part

    C.setElement((1 / sqrt(L)) * ONE_c, Indices(0, 2, 2));               // sigM first part
    Cl.setElement((1 / sqrt(L)) * ONE_c, Indices(0, 2, 2));              // sigM first part
    C.setElement((1 / sqrt(L)) * conjugate(kfac), Indices(2, D - 1, 1)); // sigP second part
    Cr.setElement((1 / sqrt(L)) * conjugate(kfac), Indices(2, 0, 1));    // sigP second part
    // the intermediate terms carry the factors, too
    // These are the only elements that should change if k changes
    C.setElement(kfac, Indices(1, 1, 0));            // between sigP and sigM
    C.setElement(conjugate(kfac), Indices(2, 2, 0)); // between sigM and sigP
    C.reshape(Indices(D * D, 4));
    Cl.reshape(Indices(D, 4));
    Cr.reshape(Indices(D, 4));
    C.multiplyRight(Z);
    Cl.multiplyRight(Z);
    Cr.multiplyRight(Z);
    C.reshape(Indices(D, D, dim, dim));
    C.permute(Indices(3, 1, 4, 2));
    Cl.reshape(Indices(1, D, dim, dim));
    Cl.permute(Indices(3, 1, 4, 2));
    Cr.reshape(Indices(D, 1, dim, dim));
    Cr.permute(Indices(3, 1, 4, 2));
    // Assign the operators to the MPO
    Nkmpo.setOp(0, new Operator(Cl), true);
    Nkmpo.setOp(L - 1, new Operator(Cr), true);
    if (L > 2)
      Nkmpo.setOp(1, new Operator(C), true);
    for (int l = 2; l < L - 1; l++)
      Nkmpo.setOp(l, &Nkmpo.getOp(1), false);
  }
}

void SpinMPO::getSkMPO(int N, int dim, int z, MPO &Skmpo, bool prepared)
{
  if (dim != 2)
  {
    cout << "Wrong dimension for SpinMPO::getMomentumDistributionMPO, right now,"
         << " only dim=2 supported" << endl;
    exit(212);
  }
  // for the time being, prepare every time
  prepared = false; // TODO: change!
  if (!prepared)
  {
    Skmpo.initLength(N);
    // the basic operators are identity and sigma-
    mwArray idOp = identityMatrix(dim);
    mwArray sigM(Indices(dim, dim));
    sigM.fillWithZero();
    sigM.setElement(ONE_c, Indices(1, 0));

    mwArray Z = mwArray(Indices(2, dim, dim));
    Z.fillWithZero();
    for (int i1 = 0; i1 < dim; i1++)
      for (int i2 = 0; i2 < dim; i2++)
      {
        Z.setElement(idOp.getElement(Indices(i1, i2)), Indices(0, i1, i2));
        Z.setElement(sigM.getElement(Indices(i1, i2)), Indices(1, i1, i2));
      }
    Z.reshape(Indices(2, dim * dim));
    // Basic coefficients
    int D = 2;
    mwArray C(Indices(D, D, 2));
    mwArray Cl(Indices(1, D, 2)), Cr(Indices(D, 1, 2));
    C.setElement(ONE_c, Indices(0, 0, 0));
    Cl.setElement(ONE_c, Indices(0, 0, 0));
    C.setElement(ONE_c, Indices(D - 1, D - 1, 0));
    Cr.setElement(ONE_c, Indices(D - 1, 0, 0));
    C.setElement(ONE_c, Indices(0, D - 1, 1));
    Cl.setElement(ONE_c, Indices(0, D - 1, 1));
    Cr.setElement(ONE_c, Indices(0, 0, 1));
    C.reshape(Indices(D * D, 2));
    Cl.reshape(Indices(D, 2));
    Cr.reshape(Indices(D, 2));
    mwArray term = C * Z;
    term.reshape(Indices(D, D, dim, dim));
    term.permute(Indices(3, 1, 4, 2));
    mwArray termL = Cl * Z;
    termL.reshape(Indices(1, D, dim, dim));
    termL.permute(Indices(3, 1, 4, 2));
    mwArray termR = Cr * Z;
    termR.reshape(Indices(D, 1, dim, dim));
    termR.permute(Indices(3, 1, 4, 2));

    complex_t kval = (-2 * M_PIl * z / N) * I_c;
    for (int l = 0; l < N; l++)
    {
      complex_t coeff = exp(kval * l) / sqrt(N);
      if (l == 0)
      {
        termL.setElement(coeff, Indices(1, 0, 0, 1));
        Skmpo.setOp(l, new Operator(termL), true);
      }
      else if (l == N - 1)
      {
        termR.setElement(coeff, Indices(1, 0, 0, 0));
        Skmpo.setOp(l, new Operator(termR), true);
      }
      else
      {
        term.setElement(coeff, Indices(1, 0, 0, 1));
        Skmpo.setOp(l, new Operator(term), true);
      }
    }
  }
}
void SpinMPO::getOBCSkMPO(int N, int dim, int z, MPO &Skmpo, bool prepared)
{
  if (dim != 2)
  {
    cout << "Wrong dimension for SpinMPO::getMomentumDistributionMPO, right now,"
         << " only dim=2 supported" << endl;
    exit(212);
  }
  // for the time being, prepare every time
  prepared = false; // TODO: change!
  if (!prepared)
  {
    Skmpo.initLength(N);
    // the basic operators are identity and sigma-
    mwArray idOp = identityMatrix(dim);
    mwArray sigM(Indices(dim, dim));
    sigM.fillWithZero();
    sigM.setElement(ONE_c, Indices(1, 0));

    mwArray Z = mwArray(Indices(2, dim, dim));
    Z.fillWithZero();
    for (int i1 = 0; i1 < dim; i1++)
      for (int i2 = 0; i2 < dim; i2++)
      {
        Z.setElement(idOp.getElement(Indices(i1, i2)), Indices(0, i1, i2));
        Z.setElement(sigM.getElement(Indices(i1, i2)), Indices(1, i1, i2));
      }
    Z.reshape(Indices(2, dim * dim));
    // Basic coefficients
    int D = 2;
    mwArray C(Indices(D, D, 2));
    mwArray Cl(Indices(1, D, 2)), Cr(Indices(D, 1, 2));
    C.setElement(ONE_c, Indices(0, 0, 0));
    Cl.setElement(ONE_c, Indices(0, 0, 0));
    C.setElement(ONE_c, Indices(D - 1, D - 1, 0));
    Cr.setElement(ONE_c, Indices(D - 1, 0, 0));
    C.setElement(ONE_c, Indices(0, D - 1, 1));
    Cl.setElement(ONE_c, Indices(0, D - 1, 1));
    Cr.setElement(ONE_c, Indices(0, 0, 1));
    C.reshape(Indices(D * D, 2));
    Cl.reshape(Indices(D, 2));
    Cr.reshape(Indices(D, 2));
    mwArray term = C * Z;
    term.reshape(Indices(D, D, dim, dim));
    term.permute(Indices(3, 1, 4, 2));
    mwArray termL = Cl * Z;
    termL.reshape(Indices(1, D, dim, dim));
    termL.permute(Indices(3, 1, 4, 2));
    mwArray termR = Cr * Z;
    termR.reshape(Indices(D, 1, dim, dim));
    termR.permute(Indices(3, 1, 4, 2));

    double kval = (M_PIl * z / (N + 1));
    for (int l = 0; l < N; l++)
    {
      complex_t coeff = ONE_c * sin(kval * (l + 1)) * sqrt(2) / sqrt(N + 1);
      if (l == 0)
      {
        termL.setElement(coeff, Indices(1, 0, 0, 1));
        Skmpo.setOp(l, new Operator(termL), true);
      }
      else if (l == N - 1)
      {
        termR.setElement(coeff, Indices(1, 0, 0, 0));
        Skmpo.setOp(l, new Operator(termR), true);
      }
      else
      {
        term.setElement(coeff, Indices(1, 0, 0, 1));
        Skmpo.setOp(l, new Operator(term), true);
      }
    }
  }
}

void SpinMPO::getSingleBodyMPO(int L, const mwArray &idOp,
                               const mwArray &spinOp, MPO &Smpo, bool staggered)
{
  Smpo.initLength(L);
  int d = idOp.getDimension(0);
  //
  // mwArray sig2(Indices(d,d),datay);
  mwArray Z = mwArray(Indices(2, d, d));
  Z.fillWithZero();
  for (int i1 = 0; i1 < d; i1++)
    for (int i2 = 0; i2 < d; i2++)
    {
      Z.setElement(idOp.getElement(Indices(i1, i2)), Indices(0, i1, i2));
      Z.setElement(spinOp.getElement(Indices(i1, i2)), Indices(1, i1, i2));
    }
  Z.reshape(Indices(2, d * d));

  int D = 2;
  mwArray C(Indices(D, D, 2));
  mwArray Cl(Indices(1, D, 2)), Cr(Indices(D, 1, 2));
  Cl.setElement(ONE_c, Indices(0, 0, 0));
  Cr.setElement(ONE_c, Indices(D - 1, 0, 0));
  C.setElement(ONE_c, Indices(0, 0, 0));
  C.setElement(ONE_c, Indices(D - 1, D - 1, 0));
  C.setElement(ONE_c, Indices(0, D - 1, 1));
  Cl.setElement(ONE_c, Indices(0, D - 1, 1));
  int signR = 1;
  if (staggered && (L - 1) % 2 != 0)
    signR = -1;
  Cr.setElement(signR * ONE_c, Indices(0, 0, 1));
  // if staggered, I set a C for odd sites
  mwArray Codd(C);
  if (staggered)
    Codd.setElement(-1 * ONE_c, Indices(0, D - 1, 1)); // only this element changes sign
  C.reshape(Indices(D * D, 2));
  Codd.reshape(Indices(D * D, 2));
  Cl.reshape(Indices(D, 2));
  Cr.reshape(Indices(D, 2));
  mwArray term = C * Z;
  term.reshape(Indices(D, D, d, d));
  term.permute(Indices(3, 1, 4, 2));
  mwArray termOdd;
  if (!staggered)
    termOdd = term;
  else
  {
    termOdd = Codd * Z;
    termOdd.reshape(Indices(D, D, d, d));
    termOdd.permute(Indices(3, 1, 4, 2));
  }
  mwArray termL = Cl * Z;
  termL.reshape(Indices(1, D, d, d));
  termL.permute(Indices(3, 1, 4, 2));
  mwArray termR = Cr * Z;
  termR.reshape(Indices(D, 1, d, d));
  termR.permute(Indices(3, 1, 4, 2));

  Smpo.setOp(0, new Operator(termL), true);
  Smpo.setOp(L - 1, new Operator(termR), true);
  if (!staggered)
  {
    if (L > 2)
    {
      Smpo.setOp(1, new Operator(term), true);
      for (int k = 2; k < L - 1; k++)
      {
        Smpo.setOp(k, &Smpo.getOp(1), false);
      }
    }
  }
  else
  {
    if (L > 2)
    {
      Smpo.setOp(1, new Operator(termOdd), true);
      if (L > 3)
      {
        Smpo.setOp(2, new Operator(term), true);
        for (int k = 3; k < L - 1; k += 2)
        {
          Smpo.setOp(k, &Smpo.getOp(1), false);
          if (k + 1 < L - 1)
            Smpo.setOp(k + 1, &Smpo.getOp(2), false);
        }
      }
    }
  }
}

void SpinMPO::getProjectorTotalSz(int L, int spindim, MPO &PSz)
{
  if (spindim != 2)
  {
    cout << "ERROR: Only spin 1/2 supported by SpinMPO" << endl;
    exit(1);
  }
  cout << "SpinMPO::getProjectorSz0()" << endl;
  int d = spindim;
  PSz.initLength(L);
  // In this case, I construct the tensors by hand
  for (int k = 0; k < L / 2; k++)
  {
    mwArray termK(Indices(d, k + 1, d, k + 2));
    for (int j = 0; j < k + 1; j++)
    {
      termK.setElement(ONE_c, Indices(0, j, 0, j));
      termK.setElement(ONE_c, Indices(1, j, 1, j + 1));
    }
    // When I have to connect both sides, an auxiliary intermediate
    // matrix ensures that indices match (to total sum 0, although it
    // can be modified to accommodate a different value of Sz)
    // And I have to take proper care of the even and odd lengths
    if ((L % 2 != 0 && k == (L - 1) / 2) ||
        (L % 2 == 0 && k == L / 2 - 1))
    {
      if (L % 2 == 0)
      {
        PSz.setOp(L - k - 1, new Operator(permute(termK, Indices(1, 4, 3, 2))), true);
      }
      mwArray transf(Indices(L / 2 + 1, L / 2 + 1));
      for (int j = 0; j < L / 2 + 1; j++)
        transf.setElement(ONE_c, Indices(j, L / 2 - j));
      termK.reshape(Indices(d * (k + 1) * d, k + 2));
      termK.multiplyRight(transf);
      termK.reshape(Indices(d, (k + 1), d, k + 2));
      PSz.setOp(k, new Operator(termK), true);
    }
    else
    {
      PSz.setOp(k, new Operator(termK), true);
      termK.permute(Indices(1, 4, 3, 2));
      PSz.setOp(L - k - 1, new Operator(termK), true);
    }
  }
}

void SpinMPO::getSx2MPO(int L, int spindim, MPO &Sz2)
{
  if (spindim != 2)
  {
    cout << "ERROR: Only spin 1/2 supported by SpinMPO" << endl;
    exit(1);
  }
  // cout<<"SpinMPO::getSz2MPO()"<<endl;
  Sz2.initLength(L);
  int d = spindim;
  // The operators
  mwArray sig0 = identityMatrix(d);
  complex_t datax[] = {ZERO_c, ONE_c, ONE_c, ZERO_c};
  mwArray sig1(Indices(d, d), datax);
  mwArray Z(Indices(2, d, d));
  for (int i = 0; i < d; i++)
    for (int j = 0; j < d; j++)
    {
      Z.setElement(sig0.getElement(Indices(i, j)), Indices(0, i, j));
      Z.setElement(sig1.getElement(Indices(i, j)), Indices(1, i, j));
    }
  Z.reshape(Indices(2, d * d));
  // The coefficients (edges are different)
  int D = 3; // bond dimension
  mwArray C(Indices(D, D, 2));
  mwArray Cl(Indices(1, D, 2));              // left edge
  mwArray Cr(Indices(D, 1, 2));              // right edge
  C.setElement(ONE_c, Indices(0, D - 1, 0)); // identity terms
  C.setElement(ONE_c, Indices(0, 0, 0));     // before any operator
  Cl.setElement(ONE_c, Indices(0, D - 1, 0));
  Cl.setElement(ONE_c, Indices(0, 0, 0));
  Cr.setElement(ONE_c, Indices(0, 0, 0));
  C.setElement(ONE_c, Indices(1, 1, 0));               // intermediate
  C.setElement(sqrt(2) * ONE_c, Indices(0, 1, 1));     // left part of ZZ
  Cl.setElement(sqrt(2) * ONE_c, Indices(0, 1, 1));    // left part of ZZ
  C.setElement(sqrt(2) * ONE_c, Indices(1, D - 1, 1)); // right part of ZZ
  Cr.setElement(sqrt(2) * ONE_c, Indices(1, 0, 1));
  C.setElement(ONE_c, Indices(D - 1, D - 1, 0)); // after every operator
  Cr.setElement(ONE_c, Indices(D - 1, 0, 0));    // after every operator
  C.reshape(Indices(D * D, 2));
  Cl.reshape(Indices(D, 2));
  Cr.reshape(Indices(D, 2));
  C.multiplyRight(Z);
  Cl.multiplyRight(Z);
  Cr.multiplyRight(Z);
  C.reshape(Indices(D, D, d, d));
  C.permute(Indices(3, 1, 4, 2));
  Cl.reshape(Indices(1, D, d, d));
  Cl.permute(Indices(3, 1, 4, 2));
  Cr.reshape(Indices(D, 1, d, d));
  Cr.permute(Indices(3, 1, 4, 2));
  // Assign the operators to the MPO
  Sz2.setOp(0, new Operator(Cl), true);
  Sz2.setOp(L - 1, new Operator(Cr), true);
  if (L > 2)
    Sz2.setOp(1, new Operator(C), true);
  for (int l = 2; l < L - 1; l++)
    Sz2.setOp(l, &Sz2.getOp(1), false);
}

void SpinMPO::getSy2MPO(int L, int spindim, MPO &Sz2)
{
  if (spindim != 2)
  {
    cout << "ERROR: Only spin 1/2 supported by SpinMPO" << endl;
    exit(1);
  }
  // cout<<"SpinMPO::getSz2MPO()"<<endl;
  Sz2.initLength(L);
  int d = spindim;
  // The operators
  mwArray sig0 = identityMatrix(d);
  complex_t dataY[] = {ZERO_c, I_c, -I_c, ZERO_c};
  mwArray sig1(Indices(d, d), dataY);
  mwArray Z(Indices(2, d, d));
  for (int i = 0; i < d; i++)
    for (int j = 0; j < d; j++)
    {
      Z.setElement(sig0.getElement(Indices(i, j)), Indices(0, i, j));
      Z.setElement(sig1.getElement(Indices(i, j)), Indices(1, i, j));
    }
  Z.reshape(Indices(2, d * d));
  // The coefficients (edges are different)
  int D = 3; // bond dimension
  mwArray C(Indices(D, D, 2));
  mwArray Cl(Indices(1, D, 2));              // left edge
  mwArray Cr(Indices(D, 1, 2));              // right edge
  C.setElement(ONE_c, Indices(0, D - 1, 0)); // identity terms
  C.setElement(ONE_c, Indices(0, 0, 0));     // before any operator
  Cl.setElement(ONE_c, Indices(0, D - 1, 0));
  Cl.setElement(ONE_c, Indices(0, 0, 0));
  Cr.setElement(ONE_c, Indices(0, 0, 0));
  C.setElement(ONE_c, Indices(1, 1, 0));               // intermediate
  C.setElement(sqrt(2) * ONE_c, Indices(0, 1, 1));     // left part of ZZ
  Cl.setElement(sqrt(2) * ONE_c, Indices(0, 1, 1));    // left part of ZZ
  C.setElement(sqrt(2) * ONE_c, Indices(1, D - 1, 1)); // right part of ZZ
  Cr.setElement(sqrt(2) * ONE_c, Indices(1, 0, 1));
  C.setElement(ONE_c, Indices(D - 1, D - 1, 0)); // after every operator
  Cr.setElement(ONE_c, Indices(D - 1, 0, 0));    // after every operator
  C.reshape(Indices(D * D, 2));
  Cl.reshape(Indices(D, 2));
  Cr.reshape(Indices(D, 2));
  C.multiplyRight(Z);
  Cl.multiplyRight(Z);
  Cr.multiplyRight(Z);
  C.reshape(Indices(D, D, d, d));
  C.permute(Indices(3, 1, 4, 2));
  Cl.reshape(Indices(1, D, d, d));
  Cl.permute(Indices(3, 1, 4, 2));
  Cr.reshape(Indices(D, 1, d, d));
  Cr.permute(Indices(3, 1, 4, 2));
  // Assign the operators to the MPO
  Sz2.setOp(0, new Operator(Cl), true);
  Sz2.setOp(L - 1, new Operator(Cr), true);
  if (L > 2)
    Sz2.setOp(1, new Operator(C), true);
  for (int l = 2; l < L - 1; l++)
    Sz2.setOp(l, &Sz2.getOp(1), false);
}

void SpinMPO::getS2MPO(int L, int spindim, MPO &S2, double alpha)
{
  if (spindim != 2)
  {
    cout << "ERROR: Only spin 1/2 supported by SpinMPO" << endl;
    exit(1);
  }
  // cout<<"SpinMPO::getS2MPO()"<<endl;
  S2.initLength(L);
  int d = spindim;
  // The operators
  mwArray sig0 = identityMatrix(d);
  complex_t datax[] = {ZERO_c, .5 * ONE_c, .5 * ONE_c, ZERO_c};
  mwArray sig1(Indices(d, d), datax);
  complex_t datay[] = {ZERO_c, .5 * I_c, -.5 * I_c, ZERO_c};
  mwArray sig2(Indices(d, d), datay);
  complex_t dataz[] = {.5 * ONE_c, ZERO_c, ZERO_c, -.5 * ONE_c};
  mwArray sig3(Indices(d, d), dataz);
  mwArray Z(Indices(4, d, d));
  for (int i = 0; i < d; i++)
    for (int j = 0; j < d; j++)
    {
      Z.setElement(sig0.getElement(Indices(i, j)), Indices(0, i, j));
      Z.setElement(sig1.getElement(Indices(i, j)), Indices(1, i, j));
      Z.setElement(sig2.getElement(Indices(i, j)), Indices(2, i, j));
      Z.setElement(sig3.getElement(Indices(i, j)), Indices(3, i, j));
    }
  Z.reshape(Indices(4, d * d));
  // The coefficients (edges are different)
  int D = 5; // bond dimension
  mwArray C(Indices(D, D, 4));
  mwArray Cl(Indices(1, D, 4));                                  // left edge
  mwArray Cr(Indices(D, 1, 4));                                  // right edge
  C.setElement((.75 - alpha / L) * ONE_c, Indices(0, D - 1, 0)); // identity terms
  C.setElement(ONE_c, Indices(0, 0, 0));                         // before any operator
  Cl.setElement((.75 - alpha / L) * ONE_c, Indices(0, D - 1, 0));
  Cl.setElement(ONE_c, Indices(0, 0, 0));
  Cr.setElement((.75 - alpha / L) * ONE_c, Indices(0, 0, 0));
  C.setElement(ONE_c, Indices(1, 1, 0));                 // intermediate
  C.setElement(ONE_c, Indices(2, 2, 0));                 // intermediate
  C.setElement(ONE_c, Indices(3, 3, 0));                 // intermediate
  C.setElement((sqrt(2)) * ONE_c, Indices(0, 3, 3));     // left part of ZZ
  Cl.setElement((sqrt(2)) * ONE_c, Indices(0, 3, 3));    // left part of ZZ
  C.setElement((sqrt(2)) * ONE_c, Indices(3, D - 1, 3)); // right part of ZZ
  Cr.setElement((sqrt(2)) * ONE_c, Indices(3, 0, 3));
  C.setElement((sqrt(2)) * ONE_c, Indices(0, 2, 2));     // left part of YY
  Cl.setElement((sqrt(2)) * ONE_c, Indices(0, 2, 2));    // left part of YY
  C.setElement((sqrt(2)) * ONE_c, Indices(2, D - 1, 2)); // right part of YY
  Cr.setElement((sqrt(2)) * ONE_c, Indices(2, 0, 2));
  C.setElement((sqrt(2)) * ONE_c, Indices(0, 1, 1));     // left part of XX
  Cl.setElement((sqrt(2)) * ONE_c, Indices(0, 1, 1));    // left part of XX
  C.setElement((sqrt(2)) * ONE_c, Indices(1, D - 1, 1)); // right part of XX
  Cr.setElement((sqrt(2)) * ONE_c, Indices(1, 0, 1));
  C.setElement(ONE_c, Indices(D - 1, D - 1, 0)); // after every operator
  Cr.setElement(ONE_c, Indices(D - 1, 0, 0));    // after every operator
  C.reshape(Indices(D * D, 4));
  Cl.reshape(Indices(D, 4));
  Cr.reshape(Indices(D, 4));
  C.multiplyRight(Z);
  Cl.multiplyRight(Z);
  Cr.multiplyRight(Z);
  C.reshape(Indices(D, D, d, d));
  C.permute(Indices(3, 1, 4, 2));
  Cl.reshape(Indices(1, D, d, d));
  Cl.permute(Indices(3, 1, 4, 2));
  Cr.reshape(Indices(D, 1, d, d));
  Cr.permute(Indices(3, 1, 4, 2));
  // Assign the operators to the MPO
  S2.setOp(0, new Operator(Cl), true);
  S2.setOp(L - 1, new Operator(Cr), true);
  if (L > 2)
    S2.setOp(1, new Operator(C), true);
  for (int l = 2; l < L - 1; l++)
    S2.setOp(l, &S2.getOp(1), false);
}

void SpinMPO::getProjectorTotalSPair(int L, int spindim, MPO &PS1, bool even)
{
  // First, construct the basic operators
  int D = 4;
  int d = spindim;
  mwArray basL(Indices(d, d, D));
  PS1.initLength(L);
  // The operators
  mwArray sig0 = sqrt(3) * .5 * identityMatrix(d);
  complex_t datax[] = {ZERO_c, .5 * ONE_c, .5 * ONE_c, ZERO_c};
  mwArray sig1(Indices(d, d), datax);
  complex_t datay[] = {ZERO_c, .5 * I_c, -.5 * I_c, ZERO_c};
  mwArray sig2(Indices(d, d), datay);
  complex_t dataz[] = {.5 * ONE_c, ZERO_c, ZERO_c, -.5 * ONE_c};
  mwArray sig3(Indices(d, d), dataz);
  for (int d1 = 0; d1 < d; d1++)
    for (int d2 = 0; d2 < d; d2++)
    {
      basL.setElement(sig0.getElement(Indices(d1, d2)), Indices(d1, d1, 0));
      basL.setElement(sig1.getElement(Indices(d1, d2)), Indices(d1, d1, 1));
      basL.setElement(sig2.getElement(Indices(d1, d2)), Indices(d1, d1, 2));
      basL.setElement(sig3.getElement(Indices(d1, d2)), Indices(d1, d1, 3));
    }
  basL.reshape(Indices(d, 1, d, D));
  Operator *OL = new Operator(basL);
  mwArray basR(basL);
  basR.permute(Indices(1, 4, 3, 2));
  // These are already the first and last operators
  Operator *OR = new Operator(basR);
  // I might need also an identity operator, if even=false (for site
  // 0) and/or even and last is even (odd nr of sites) or not even and
  // last is odd (even nr of sites)
  int firstL = even ? 0 : 1;
  mwArray idOp = identityMatrix(d);
  idOp.reshape(Indices(d, 1, d, 1));
  bool usedId = false;
  if (firstL == 1)
  {
    PS1.setOp(0, new Operator(idOp), true);
    usedId = true;
  }
  int k = firstL;
  bool usedL = false;
  bool usedR = false;
  for (k = firstL; k < L - 1; k++)
  {
    if (usedL == false)
    {
      PS1.setOp(k, OL, true);
      usedL = true;
    }
    else
    {
      PS1.setOp(k, OL, false);
    }
    k++;
    if (usedR == false)
    {
      PS1.setOp(k, OR, true);
      usedR = true;
    }
    else
    {
      PS1.setOp(k, OR, false);
    }
  }
  if (k == L - 1)
  { // after the last OR was placed in L-2,
    // the counter is incremented
    PS1.setOp(k, new Operator(idOp), true);
    // could put a pointer to the first one, if there, but small overhead
  }
}

void SpinMPO::getModulatedSingleBodyMPO(int N, const mwArray &idOp,
                                        const mwArray &spinOp, double q, MPO &Smpo, bool obc)
{
  Smpo.initLength(N);
  int d = idOp.getDimension(0);
  //
  // mwArray sig2(Indices(d,d),datay);
  mwArray Z = mwArray(Indices(2, d, d));
  Z.fillWithZero();
  for (int i1 = 0; i1 < d; i1++)
    for (int i2 = 0; i2 < d; i2++)
    {
      Z.setElement(idOp.getElement(Indices(i1, i2)), Indices(0, i1, i2));
      Z.setElement(spinOp.getElement(Indices(i1, i2)), Indices(1, i1, i2));
    }
  Z.reshape(Indices(2, d * d));

  // Basic coefficients
  int D = 2;
  for (int l = 0; l < N; l++)
  {
    complex_t factorQ = exp(-l * q * M_PIl * I_c);
    if (obc)
      factorQ = sin(l * q * M_PIl) * ONE_c;
    //    cout<<"Phase factor for site "<<l<<" with q="<<q<<"=>"<<factorQ<<endl;
    int Dl = D;
    int Dr = D;
    if (l == 0)
      Dl = 1;
    if (l == N - 1)
      Dr = 1;
    mwArray C(Indices(Dl, Dr, 2));
    if (l < N - 1)
      C.setElement(ONE_c, Indices(0, 0, 0));
    if (l > 0)
      C.setElement(ONE_c, Indices(Dl - 1, Dr - 1, 0));
    C.setElement(factorQ, Indices(0, Dr - 1, 1));
    C.reshape(Indices(Dl * Dr, 2));
    mwArray term = C * Z;
    term.reshape(Indices(Dl, Dr, d, d));
    term.permute(Indices(3, 1, 4, 2));
    Smpo.setOp(l, new Operator(term), true);
  }
}

void SpinMPO::getSxSxMPO(int L, int spindim, MPO &SxSx)
{
  if (spindim != 2)
  {
    cout << "ERROR: Only spin 1/2 supported by SpinMPO" << endl;
    exit(1);
  }
  cout << "SpinMPO::getSxSxMPO()" << endl;
  int d = spindim;
  mwArray sig0 = identityMatrix(d);
  complex_t datax[] = {ZERO_c, ONE_c, ONE_c, ZERO_c};
  mwArray sig1(Indices(d, d), datax);
  getNearestNeighbourMPO(L, sig0, .5 * sig1, .5 * sig1, SxSx);
}

void SpinMPO::getSySyMPO(int L, int spindim, MPO &SySy)
{
  if (spindim != 2)
  {
    cout << "ERROR: Only spin 1/2 supported by SpinMPO" << endl;
    exit(1);
  }
  cout << "SpinMPO::getSySyMPO()" << endl;
  int d = spindim;
  mwArray sig0 = identityMatrix(d);
  complex_t datay[] = {ZERO_c, I_c, -1. * I_c, ZERO_c};
  mwArray sig1(Indices(d, d), datay);
  getNearestNeighbourMPO(L, sig0, .5 * sig1, .5 * sig1, SySy);
}

void SpinMPO::getSzSzMPO(int L, int spindim, MPO &SzSz)
{
  if (spindim != 2)
  {
    cout << "ERROR: Only spin 1/2 supported by SpinMPO" << endl;
    exit(1);
  }
  cout << "SpinMPO::getSzSzMPO()" << endl;
  int d = spindim;
  mwArray sig0 = identityMatrix(d);
  complex_t dataz[] = {ONE_c, ZERO_c, ZERO_c, -1. * ONE_c};
  mwArray sig1(Indices(d, d), dataz);
  getNearestNeighbourMPO(L, sig0, .5 * sig1, .5 * sig1, SzSz);
}

void SpinMPO::getNearestNeighbourMPO(int L, const mwArray &idOp,
                                     const mwArray &opA, const mwArray &opB,
                                     MPO &Smpo)
{
  // cout<<"SpinMPO::getSz2MPO()"<<endl;
  Smpo.initLength(L);
  int d = idOp.getDimension(0);
  mwArray Z(Indices(3, d, d));
  for (int i = 0; i < d; i++)
    for (int j = 0; j < d; j++)
    {
      Z.setElement(idOp.getElement(Indices(i, j)), Indices(0, i, j));
      Z.setElement(opA.getElement(Indices(i, j)), Indices(1, i, j));
      Z.setElement(opB.getElement(Indices(i, j)), Indices(2, i, j));
    }
  Z.reshape(Indices(3, d * d));
  // The coefficients (edges are different)
  int D = 3; // bond dimension
  mwArray C(Indices(D, D, 3));
  mwArray Cl(Indices(1, D, 3));          // left edge
  mwArray Cr(Indices(D, 1, 3));          // right edge
  C.setElement(ONE_c, Indices(0, 0, 0)); // before any operator
  Cl.setElement(ONE_c, Indices(0, 0, 0));
  C.setElement(ONE_c, Indices(0, 1, 1)); // left part of AB
  Cl.setElement(ONE_c, Indices(0, 1, 1));
  C.setElement(ONE_c, Indices(1, D - 1, 2)); // right part of AB
  Cr.setElement(ONE_c, Indices(1, 0, 2));
  C.setElement(ONE_c, Indices(D - 1, D - 1, 0)); // after every operator
  Cr.setElement(ONE_c, Indices(D - 1, 0, 0));
  C.reshape(Indices(D * D, 3));
  Cl.reshape(Indices(D, 3));
  Cr.reshape(Indices(D, 3));
  C.multiplyRight(Z);
  Cl.multiplyRight(Z);
  Cr.multiplyRight(Z);
  C.reshape(Indices(D, D, d, d));
  C.permute(Indices(3, 1, 4, 2));
  Cl.reshape(Indices(1, D, d, d));
  Cl.permute(Indices(3, 1, 4, 2));
  Cr.reshape(Indices(D, 1, d, d));
  Cr.permute(Indices(3, 1, 4, 2));
  // Assign the operators to the MPO
  Smpo.setOp(0, new Operator(Cl), true);
  Smpo.setOp(L - 1, new Operator(Cr), true);
  if (L > 2)
    Smpo.setOp(1, new Operator(C), true);
  for (int l = 2; l < L - 1; l++)
    Smpo.setOp(l, &Smpo.getOp(1), false);
}
