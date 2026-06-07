#include <stdlib.h>
#include <math.h>

//Erzeugt gleichverteilte Zufallszahlen auf [0,1]
double nulleins();

//Erzeugt normalverteilte Zufallszahlen
double normal(double mu, double sigma);


/* gamma verteilte Zufallszahlen; Bernardo Smith Notation */
double RNDGAM(double a, double b);

//int rbinom(double p,int n);
