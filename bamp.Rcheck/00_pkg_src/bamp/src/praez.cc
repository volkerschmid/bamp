
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "zufall.h"

// *********** Gibbs-Steps zur Erzeugung der Pr?zisionen (RW1/2/Space) *************
double hyper(int rw, double* theta, double k_a, double k_b, int n)
{

//double k_a = 1.0;
//double k_b = 0.00005;

double aa;
double bb;
double kappa_neu;

if (rw==1)
{
double summe = 0.0;



aa = k_a + 0.5 * double(n-1);

for (int i=1; i < n; i++)
{
    summe = summe + (theta[i] - theta[i-1]) * (theta[i] - theta[i-1]);
}

bb = k_b + 0.5 * summe;

kappa_neu = RNDGAM(aa,bb);
}
else
{

double dopp_diff;
double summe = 0.0;

aa = k_a + 0.5 * double(n-2);

for (int i=1; i < (n - 1); i++)
{
	dopp_diff = theta[i-1] - 2*theta[i] + theta[i+1];
        summe = summe + dopp_diff * dopp_diff;
}

bb = k_b + 0.5 * summe;

kappa_neu = RNDGAM(aa,bb);
}

return kappa_neu;
}



double hyper_a(double hyper2, int rw, double* theta, double k_a, double k_b, int n)
{

//double k_a = 1.0;
//double k_b = 0.00005;

double aa;
double bb;
double kappa_neu;

if (rw==1)
{
double summe = 0.0;



aa = k_a + 0.5 * double(n);

for (int i=1; i < n; i++)
{
    summe = summe + (theta[i] - theta[i-1]) * (theta[i] - theta[i-1]);
}

bb = k_b + (0.5 * summe / hyper2);

kappa_neu = RNDGAM(aa,bb);
}
else
{

double dopp_diff;
double summe = 0.0;

aa = k_a + 0.5 * double(n);

for (int i=1; i < (n - 1); i++)
{
	dopp_diff = theta[i-1] - 2*theta[i] + theta[i+1];
        summe = summe + dopp_diff * dopp_diff;
}

bb = k_b + 0.5 * summe;

kappa_neu = RNDGAM(aa,bb);
}

return kappa_neu;
}

double hyper2(double* z, double d_g, double d_h, int n)
{
//double d_g = 1.0;
//double d_h = 0.005;

double gg;
double hh;

double summe = 0.0;

double delta_neu;

gg = d_g + 0.5 * double(n);

for (int i = 0; i<n; i++)
{
	 summe = summe + (z[i] * z[i]);
}


hh = d_h + 0.5 * summe;

delta_neu = RNDGAM(gg,hh);

return delta_neu;
}

double delta_berechnen(double** z, double d_g, double d_h, int n_o_a, int n_o_p)
{
//double d_g = 1.0;
//double d_h = 0.005;

double gg;
double hh;

double summe = 0.0;

double delta_neu;

gg = d_g + 0.5 * double(n_o_a) * double(n_o_p);

for (int i = 0; i<n_o_a; i++)
{
	for (int j=0; j<n_o_p; j++)
	{
		 summe = summe + (z[i][j] * z[i][j]);
	}
}


hh = d_h + 0.5 * summe;
 // cout << endl<<hh<<":" << gg<< " ";

delta_neu = RNDGAM(gg,hh);

return delta_neu;
}

double delta_berechnen_S(double*** z, double d_g, double d_h, int n_o_a, int n_o_p, int n_o_r)
{
//double d_g = 1.0;
//double d_h = 0.005;

double gg;
double hh;

double summe = 0.0;

double delta_neu;

gg = d_g + 0.5 * double(n_o_a) * double(n_o_p);

for (int i = 0; i<n_o_a; i++)
{
	for (int j=0; j<n_o_p; j++)
	{
		 for (int k=0; k< n_o_r; k++)
		   {
		     summe = summe + (z[i][j][k] * z[i][j][k]);
		   }
	}
}


hh = d_h + 0.5 * summe;

delta_neu = RNDGAM(gg,hh);

return delta_neu;
}


double tau_berechnen(double* alpha, double k_a, double k_b, int** k_alpha, int number_of_regions)
{


double aa;
double bb;


double summe = 0.0;

double tau_neu;

aa = k_a + 0.5 * (double(number_of_regions)-1);

for (int i=0; i < (number_of_regions-1); i++)
{
	for (int j=i+1; j<number_of_regions; j++)
	{
		if (k_alpha[i][j]==-1)
			{summe = summe + (alpha[i] - alpha[j]) * (alpha[i] - alpha[j]);}
	}
}


bb = k_b + 0.5 * summe;

tau_neu = RNDGAM(aa,bb);

return tau_neu;
}

