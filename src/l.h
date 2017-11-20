//#include "mystring.h"
#include <string>

void zentriere(double &my, double* theta, int noa);

void start(double* theta,double*  phi,double*  psi,int number_of_agegroups,int number_of_periods, int number_of_cohorts,int** y,int** n, int vielfaches_der_breite);
//Erzeugt den Kohortenindex (ab 1!!)
//  erste Werte alter, periode werden mit 0 ?bergeben
int coh(int altersgruppe, int periode,int number_of_agegroups, int vielfaches_der_breite);

//Berechnet den ksi-Wert aus psi am Anfang

void ksi_berechnen(double** ksi, double* psi,int vielfaches_der_breite, int noa, int nop);

// Berechnet den neuen z-Wert aus ksi

void z_aus_ksi_berechnen(double** z, double my, double** ksi, double* theta, double* phi, double* psi,int vielfaches_der_breite, int noa, int nop);

//sortiert matrix theta mit laenge noa*noe, jeweils in den Spalten
void sortieren(double* theta, int noa, int noe);

//Minimum
double min(double i, double j);

// Setzt Matrix yn gleich null
void nullmatrix(int** yn, int number_of_agegroups, int number_of_periods);

// TUNING-PROGRAMM alt
void tune(double** accepted, double** konstante, int &need, int number_of_agegroups, int number_of_periods);

//Maximum
int max(int a, int b);

//schreibt die Elemente q1 bis q5 aus dem Vector my (oder spalte index der Matrix my mit spaltenanzahl laenge) in die Datei namens mystring. (davor schreibt sie vorher, zum schlu? Zeilenumbruch
double ABS(double x);
