#include "zufall.h"
#include "l.h"
#include "block.h"

void update_my_1(double &my, double** ksi, double* theta, double* phi, double* psi, int vielfaches_der_breite, int number_of_agegroups, int number_of_periods, double delta)
{
  my = 0.0;
 
 for (int i=0; i<number_of_agegroups; i++)
  {
  for (int t=0; t<number_of_periods; t++)
	  {
	  my+= ksi[i][t]-theta[i]-phi[t]-psi[coh(i,t,number_of_agegroups,vielfaches_der_breite)-1];
	}
  }
 
 
 my = my/number_of_periods;
 my = my/number_of_agegroups;
 
 my = my + normal(0,1)*sqrt(1/(delta*number_of_agegroups*number_of_periods));
 return;
}


void update_my_mh(double &my, double** ksi, double* theta, double* phi, double* psi, int vielfaches_der_breite, int number_of_agegroups, int number_of_periods, int lung_summe, int** n, int** y, int& ja_my)
{

//Update my
double q = 0.0;
double m = double(lung_summe);

for (int i=0; i< number_of_agegroups; i++)
{
  for (int t=0; t<number_of_periods; t++)
    {
     
 q += double(n[i][t])*taylor2(my+theta[i]+phi[t]+psi[coh(i,t,number_of_agegroups,vielfaches_der_breite)-1]);
 m += double(n[i][t])*my*taylor2(my+theta[i]+phi[t]+psi[coh(i,t,number_of_agegroups,vielfaches_der_breite)-1])-double(n[i][t])*taylor1(my+theta[i]+phi[t]+psi[coh(i,t,number_of_agegroups,vielfaches_der_breite)-1]);

    }
}

double prop = normal(m/q,1.0/q);

//  cout << m <<":"<< q<<" = "<<m/q<<endl;

double loglik_alt=0.0;
double loglik_neu=0.0;

for (int i=0; i< number_of_agegroups; i++)
{
  for (int t=0; t< number_of_periods; t++)
    {
      loglik_alt += my*double(y[i][t])-double(n[i][t])*log(1+exp(my+theta[i]+phi[t]+psi[coh(i,t,number_of_agegroups,vielfaches_der_breite)-1]));
      loglik_neu += prop*double(y[i][t])-double(n[i][t])*log(1+exp(prop+theta[i]+phi[t]+psi[coh(i,t,number_of_agegroups,vielfaches_der_breite)-1]));
    }
}

double q_alt = 0.0;
double m_alt = double(lung_summe);

for (int i=0; i< number_of_agegroups; i++)
{
  for (int t=0; t<number_of_periods; t++)
    {
     
 q_alt += double(n[i][t])*taylor2(prop+theta[i]+phi[t]+psi[coh(i,t,number_of_agegroups,vielfaches_der_breite)-1]);
 m_alt += double(n[i][t])*prop*taylor2(prop)-double(n[i][t])*taylor1(prop+theta[i]+phi[t]+psi[coh(i,t,number_of_agegroups,vielfaches_der_breite)-1]);

    }
}



double logdichte_alt=log(q_alt)/2.0-(my-m_alt)*(my-m_alt)*q_alt/2.0;
double logdichte_neu=log(q)/2.0-(prop-m)*(prop-m)*q/2.0;
// cout <<loglik_neu << logdichte_alt << loglik_alt<< logdichte_neu;


double alph=loglik_neu + logdichte_alt - loglik_alt - logdichte_neu;

//  cout << my << " -> "<<prop<< " : "<<alph << endl;
alph =exp(alph);


 if (alph>nulleins())
 {
  my=prop;
  ja_my++;
 }

return;

}
