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
  int* lda=new int[1];
  lda[0] = noa;
  int* info=new int[1];
  info[0]=0;

  double* L2=new double[noa*n];
  for (int i=0; i< (noa*n); i++)
    {
      L2[i]=L[i];
    }
  int* nn=new int[1];
  nn[0]=n;
  F77_CALL(dpotri)("L", nn, L2, lda, info FCONE);

  //dpotri_(uplo, n, L2, lda, info);

  double* x = new double[noa];
  for (int i=0; i<noa; i++)
    {
      x[i]=1.0;
    }
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
int* info=new int[1];
info[0]=1;
int* nn=new int[1];
nn[0]=n;
int* bw0=new int[1];
bw0[0]=bw;
int* lda=new int[1];
lda[0]=bw+1;
cholesky77(nn, bw0, matrix, lda, info);
return matrix;
}


void loese(double* A, double* z, int &n, int& bw)
{
  int* nn=new int[1];
nn[0]=n;
int* lda=new int[1];
lda[0]=bw+1;
int* bw0=new int[1];
bw0[0]=bw;
int* incx=new int[1];
incx[0]=1;

loese77T(A, z, nn, bw0, lda, incx);
return;
}

void loese2(double* A, double* z, int &n, int& bw)
{
  int* nn=new int[1];
nn[0]=n;
int* bw0=new int[1];
bw0[0]=bw;
int* lda=new int[1];
lda[0]=bw+1;
int* incx=new int[1];
incx[0]=1;
loese77N(A, z, nn, bw0, lda, incx);
return;
}
