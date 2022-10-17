#include <stdlib.h>
#include <math.h>


//Methoden fuer Bandmatrizen


// macht eine Matrix zu einer Bandmatrix
double* mxsortieren(int n, int** matrix, int* permut, int& bandwidth);


//Macht Cholesky-Zerlegung A=LL^T
double* cholesky(int n, double* matrix,  int& bw);


//loest A^T * u = z
void loese(double* A, double* z, int &n, int& bw);


//loest A * u = z
void loese2(double* A, double* z, int &n, int& bw);

// berechne x^t M x fuer Bandmatrizen
double xMx(double* Q, double* x, int noa, int b);
double xLx(double* Q, double* x, int noa, int b);

// Determinante einer symmetrischen Matrix
double det(double* A,int n);


double berechneVarianzsumx_i(double* L,int noa,int bandw);
