#define USE_FC_LEN_T
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "l.h"

#ifndef FCONE
# define FCONE
#endif

double xLx(double* Q, double* x, int noa, int b)
{
  double erg=0.0;

  for (int i=0; i<noa; i++)
    {
      for (int j=0; j<noa; j++)
	{
	  if (ABS(i-j) < b)
	    {
	      erg += x[i]*x[j]*Q[(int)(min(i,j)*b+ABS(i-j))]*Q[(int)(min(i,j)*b+ABS(i-j))];
	    }
	}
    }

  return erg;
}
double xMx(double* Q, double* x, int noa, int b)
{
  double erg=0.0;

  for (int i=0; i<noa; i++)
    {
      for (int j=0; j<noa; j++)
	{
	  if (ABS(i-j) < b)
	    {
	      erg += x[i]*x[j]*Q[(int)(min(i,j)*b+ABS(i-j))];
	    }
	}
    }

  return erg;
}
double berechneVarianzsumx_i(double* L,int noa,int bandw)
{
  int n = bandw+1;
  int lda = noa;
  int info = 0;
  int nn = n;

  double* L2=new double[noa*n];
  for (int i=0; i< (noa*n); i++)
    L2[i]=L[i];

  F77_CALL(dpotri)("L", &nn, L2, &lda, &info FCONE);

  double* x = new double[noa];
  for (int i=0; i<noa; i++)
    x[i]=1.0;
  double S=xMx(L2, x, noa, bandw+1);

  delete[] x;
  delete[] L2;

  if (S<0){S=1e-5;}
  return (1.0/S);
}


double det(double* A,int n)
{
  double erg=0.0;

  // Fortran-Methode: DSYEV

  if (n==1){erg=A[0];}
  if (n==2){erg=A[0]*A[3]-(A[1]*A[2]);}
  return erg;
}

void cholesky77(int* n, int* bw, double* matrix,  int* lda, int* info) {
  F77_CALL(dpbtrf)("L", n, bw, matrix, lda, info FCONE);
}

void loese77T(double* A, double* z, int* n, int* bw, int* lda, int* incx) {
  F77_CALL(dtbsv)("L","T","N", n, bw, A, lda, z, incx FCONE FCONE FCONE);
}
void loese77N(double* A, double* z, int* n, int* bw, int* lda, int* incx) {
  F77_CALL(dtbsv)("L","N","N", n, bw, A, lda, z, incx FCONE FCONE FCONE);
}



double* cholesky(int n, double* matrix,  int& bw)
{
  int info = 0;
  int nn = n;
  int bw0 = bw;
  int lda = bw+1;
  int sz = n * lda;

  cholesky77(&nn, &bw0, matrix, &lda, &info);

  if (info != 0) {
    /* Matrix not positive definite (common with RW2 at extreme hyperparameter
       values during tuning). Save original, add increasing diagonal jitter,
       and retry. Diagonal elements are at matrix[j*lda] in band storage. */
    double* orig = new double[sz];
    for (int i = 0; i < sz; i++) orig[i] = matrix[i];

    double jitter = 1e-6;
    for (int attempt = 0; attempt < 15 && info != 0; attempt++) {
      for (int i = 0; i < sz; i++) matrix[i] = orig[i];
      for (int j = 0; j < n; j++) matrix[j * lda] += jitter;
      info = 0;
      cholesky77(&nn, &bw0, matrix, &lda, &info);
      jitter *= 10.0;
    }
    delete[] orig;

    if (info != 0)
      error("Cholesky decomposition failed: matrix is not positive definite");
  }

  return matrix;
}


void loese(double* A, double* z, int &n, int& bw)
{
  int nn = n;
  int lda = bw+1;
  int bw0 = bw;
  int incx = 1;
  loese77T(A, z, &nn, &bw0, &lda, &incx);
  return;
}

void loese2(double* A, double* z, int &n, int& bw)
{
  int nn = n;
  int bw0 = bw;
  int lda = bw+1;
  int incx = 1;
  loese77N(A, z, &nn, &bw0, &lda, &incx);
  return;
}
