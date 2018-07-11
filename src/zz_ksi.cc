#include <stdlib.h>
#include <math.h>

#include "l.h"
#include "zufall.h"

void ZZ_aus_fc_von_ksi0(double my, double* theta, double* phi, double* psi, double** ksi, double delta, int number_of_agegroups, int number_of_periods, int vielfaches_der_breite, int** y, int** n, int** yes, int** no, int** schalter)
{
double atemp= 1.0/((exp(1.0/delta))-1.0);
double btemp;
double proposal;
double alpha;
double offset;

  // ** Metropolis **
for (int j=0; j < number_of_agegroups; j++)
      {
	for (int k=0; k < number_of_periods; k++)
	  {

offset=my + theta[j] + phi[k] + psi[coh(j,k,number_of_agegroups,vielfaches_der_breite)-1];

//  mytemp=exp(offset+(0.5/delta));

//  sigma2temp=exp((2.0*offset)+(1.0/delta))*(exp(1.0/delta)-1.0);
	 
// atemp=mytemp*mytemp/sigma2temp;   
	    
// btemp= mytemp/sigma2temp;

// atemp=exp(4*offset+2/delta-2*offset*exp(1/delta)-(exp(1/delta)/delta));
// btemp=exp(3*offset+1.5/delta-2*offset*exp(1/delta)-(exp(1/delta)/delta));

   
btemp=atemp/exp(offset+(1.0/(2.0*delta)));
         
if(schalter[j][k]==0)
  {
proposal = log ( RNDGAM(  ((double)y[j][k]) + atemp ,((double)n[j][k]) + btemp ));
 
alpha = ( (proposal-ksi[j][k]) * ((delta*offset)-atemp) ) + ( delta * (ksi[j][k]*ksi[j][k] - proposal*proposal) /2.0) + ((btemp+n[j][k]) * (exp(proposal)-exp(ksi[j][k]))) + n[j][k] * (log(1+exp(ksi[j][k])) - log(1+exp(proposal)));
  }
else
{ 
proposal = log ( RNDGAM( atemp , btemp ));

alpha = ( (proposal-ksi[j][k]) * (y[j][k]+(delta*offset)-atemp) ) + ( delta * (ksi[j][k]*ksi[j][k] - proposal*proposal) /2.0) + btemp * (exp(proposal)-exp(ksi[j][k])) + n[j][k] * (log(1+exp(ksi[j][k])) - log(1+exp(proposal)));
}

alpha=exp(alpha);
 


//  cout << alpha << " ";
	    if (nulleins() <= alpha)
	      {
		ksi[j][k]=proposal;
		yes[j][k]+=1 ;
	      }
	    else { no[j][k]+=1;}
	    
	    
      	  } 
// 	cout << endl;
      }

//  cout << endl;

return;
}

		
