#include <stdlib.h>
#include <math.h>

//Erzeugt gleichverteilte Zufallszahlen auf [0,1]

double nulleins();




//Erzeugt diskrete gleichverteilte Zufallszahlen auf 1,...,k

int einsk(int k);


//Erzeugt diskrete gleichverteilte Zufallszahlen auf 1,...,k
//gibt double zurueck
double d_einsk(int k);


//Erzeugt diskrete gleichverteilte Zufallszahlen auf 0,...,k

int nullk(int k);


//Erzeugt gleichverteilte Zufallszahlen auf [0,L]

double nullell(double l);


//Erzeugt gleichverteilte Zufallszahlen auf [a,b]

double gleichab(double a, double b);


//Erzeugt gleichverteilte Zufallszahlen auf [-1,1]

double reins();



//Erzeugt gleichverteilte Zufallszahlen auf [-0.5,0.5]

double halbeins();


//Erzeugt exonentialverteilte Zufallszahlen

double expo(double lambda);



//Erzeugt poissonverteilte Zufallszahlen

int poisson(double lambda);



//Erzeugt normalverteilte Zufallszahlen


double zahl(double mu, double sigma);
double normal(double mu, double sigma);


/* gamma verteilte Zufallszahlen; Bernardo Smith Notation */
double RNDGAM(double a, double b);

int rbinom(double p,int n);
  