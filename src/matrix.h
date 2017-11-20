//multipliziert A mit B; A ist (a x n), B ist (n x b)  
void multipliziere(double* A, double* B, int a, int n, int b, double* ergebnis);

// BERECHNET A-1, k ist matrixlaenge
void invers(double* A, int k);

//schreibt eine Matrix auf Console
void mxschreibe(double* A, int  a, int b);


//testet, ob Matrix symmetrisch und Zeilensumme null (1=fehler)
int mxcheck(int n, int** matrix);
