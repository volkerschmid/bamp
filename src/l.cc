//#include <mystring.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string>
#include "zufall.h"

double logit(double x)
{

  if (x==0)
    {
      return -20.0;
    }
  if (x==1)
    {
      return 20.0;
    }
  return (log(x/(1-x)));
}


void zentriere(double &my, double* theta, int noa)
{
  double sum=0.0;
  for (int i=0; i< noa; i++)
    {
      sum += theta[i];
    }
  for (int i=0; i< noa; i++)
    {
      theta[i] = (theta[i]-(sum/(double)noa));
    }
  my += sum/(double)noa;

  return;
}

//Erzeugt den Kohortenindex (ab 1)

int coh(int altersgruppe, int periode,int number_of_agegroups, int vielfaches_der_breite)
//  erste Werte alter, periode werden mit 0 ?bergeben
{

int cohorte;

cohorte = (number_of_agegroups - (altersgruppe+1))   * vielfaches_der_breite + (periode +1);

return cohorte;
}


void start(double* theta,double*  phi,double*  psi,int number_of_agegroups,int number_of_periods, int number_of_cohorts,int** y,int** n, int vielfaches_der_breite)
{

  int* nwerte=new int[number_of_cohorts];

  for (int i=0; i< number_of_periods; i++)
    {
      phi[i]=0.0;
    }
  for (int i=0; i<number_of_cohorts; i++)
    {
      psi[i]=0.0;
      nwerte[i]=0;
    }

  for (int i=0; i< number_of_agegroups; i++)
    {
      theta[i]=0.0;
      for (int t=0; t< number_of_periods; t++)
	{
	  theta[i] += double(y[i][t])/(double(n[i][t])*double(number_of_periods));
	  phi[t] += double(y[i][t])/double((n[i][t])*double(number_of_agegroups));
	  psi[coh(i,t,number_of_agegroups,vielfaches_der_breite)-1] += double(y[i][t])/double(n[i][t]);
	  nwerte[coh(i,t,number_of_agegroups,vielfaches_der_breite)-1] +=1;
	}
      theta[i]=logit(theta[i]);

    }
   for (int i=0; i<number_of_cohorts; i++)
    {
      psi[i]= logit(psi[i]/((double)nwerte[i]));
    }

  for (int i=0; i< number_of_periods; i++)
    {
      phi[i]=logit(phi[i]);
    }

   delete[] nwerte;

   return;
}

//Berechnet den neuen ksi-Wert

void ksi_berechnen(double** ksi, double* psi,int vielfaches_der_breite, int noa, int nop)
{

for (int i=0; i<noa; i++)
{
       for (int j=0; j<nop; j++)
       {
             ksi[i][j]=psi[coh(i,j,noa, vielfaches_der_breite)-1];
       }
}
return;
}

// Berechnet den neuen z-Wert

void z_aus_ksi_berechnen(double** z, double my, double** ksi, double* theta, double* phi, double* psi,int vielfaches_der_breite, int noa, int nop)
{
for (int i=0; i<noa; i++)
{
       for (int j=0; j<nop; j++)
       {
       	z[i][j]=ksi[i][j] - my -theta[i] - phi[j] - psi[coh(i,j, noa, vielfaches_der_breite)-1];
       }
}
return;
}

// Berechnet pr

// void  pr_aus_ksi_berechnen(matrix &pr,matrix &ksi)
// {
//
// for (int i = 0; i<ksi.rows() ; i++)
//   {
//       for (int j = 0; j < ksi.cols();j++)
// 	{
//                         pr.put(i,j, exp(ksi.get(i,j)) / (1+exp(ksi.get(i,j))) );
// 	}
//   }
// return;
// }



void sort(double* vector, int low, int high, int index, int laenge)
{
int i,j;
double x,y;

  i = low; j = high; x = vector[index+(((low+high)/2)*laenge)];
   do {
      while ( vector[index+(i*laenge)] < x ) { i++; }
      while ( vector[index+(j*laenge)] > x ) { if (j) j--; }
      if ( i <= j )
      {
         y = vector[index+(i*laenge)];
         vector[index+(i*laenge)] = vector[index+(j*laenge)];
         vector[index+(j*laenge)]= y ;
         i++; if(j) j--;
      }
   } while ( i <= j );

   if ( low < j )
      sort(vector, low, j, index, laenge);
   if ( i < high )
      sort(vector, i, high, index, laenge);
}


void sortieren(double* theta, int noa, int noe)
{
	for (int i=0; i< noa; i++)
	{
		sort(theta, 0, noe-1, i, noa);
	}
return;
}


double min(double i, double j)
{
 if (i<j)
 {
 return i;
 }
 else
 {
 return j;
 }
}
double betrag(double x)
{
if (x<0)
{
return -x;
}
else
{
return x;
}
}


void center(double* theta,double* phi,double* psi,int number_of_agegroups ,int number_of_periods ,int number_of_cohorts,int vdb,double x){


double sum = 0.0, adjust;

for (int l = 0; l < number_of_agegroups;l++){
  sum = sum + theta[l];
}

adjust = sum/double(number_of_agegroups);

for (int g = 0; g < number_of_agegroups; g++){
  theta[g] -=  adjust;
}

for (int g = 0; g < number_of_cohorts; g++){
  psi[g] += adjust;
}

sum = 0.0;

for (int l = 0; l < number_of_periods;l++){
  sum = sum + phi[l];
}

adjust = sum/double(number_of_periods);

for (int g = 0; g < number_of_periods; g++){
  phi[g] -= adjust;
}

for (int g = 0; g < number_of_cohorts; g++){
  psi[g] += adjust;
}
return;

}

void nullmatrix(int** yn, int number_of_agegroups, int number_of_periods)
  {
    for (int i=0;i<number_of_agegroups;i++)
      {
	for (int j=0;j< number_of_periods; j++)
	  {
	    yn[i][j]=0;
	  }
      }
return;
  }

// TUNING-PROGRAMM

void tune(double** accepted, double** konstante, int &need, int number_of_agegroups, int number_of_periods){

  //accepted = accepted * double(100);

  for (int i=0; i <number_of_agegroups; i++){
    for (int j=0; j <number_of_periods; j++){
      if (accepted[i][j] < .10){
        konstante[i][j]=konstante[i][j]/2;
        need +=1;
      }
      else if (accepted[i][j] < .20){
        konstante[i][j]= konstante[i][j]/1.4;
        need +=1;
      }
      else if (accepted[i][j] < .25){
        konstante[i][j]= konstante[i][j]/1.2;
        need +=1;
      }
            else if (accepted[i][j] < .30){
        konstante[i][j]= konstante[i][j]/1.1;
        need +=1;
      }
            else if (accepted[i][j] < .35){
        konstante[i][j]= konstante[i][j]/1.05;
      }

      if (accepted[i][j] > .90){
        konstante[i][j]= konstante[i][j]*2;
        need +=1;
      }
      else if (accepted[i][j] > .70) {
        konstante[i][j]= konstante[i][j]*1.5;
        need +=1;
      }
      else if (accepted[i][j] > .65){
        konstante[i][j]= konstante[i][j]*1.4;
        need +=1;
      }
      else if (accepted[i][j] > .60){
        konstante[i][j]= konstante[i][j]*1.3;
        need +=1;
      }
      else if (accepted[i][j] > .55){
        konstante[i][j]= konstante[i][j]*1.2;
        need +=1;
      }
      else if (accepted[i][j] > .50){
        konstante[i][j]= konstante[i][j]*1.1;
        need +=1;
      }
      else if (accepted[i][j] > .45){
        konstante[i][j]= konstante[i][j]*1.05;
      }
    }
  }

  return;
}

void rette(double* matrix, double* speicher, int laenge)
{

for (int i=0; i< laenge; i++)
  {
    speicher[i]=matrix[i];
  }
return;
}


void adjust(double* theta,double* phi,double* psi,int number_of_agegroups ,int number_of_periods ,int number_of_cohorts,int vdb,double x)
{

//double mu, summe = 0.0, faktor, alpha;
double noa = (double)(number_of_agegroups);
double nop = (double)(number_of_periods);
double noc = (double)(number_of_cohorts);



// for (int g = 0; g < number_of_cohorts; g++){
//   summe += psi[g];
// }
// mu = summe / noc;

// faktor = vdb * ((double)(noa+1)/2.0 - (double)noa) - (double)(nop+1)/2.0
// 				  + (double)(noc+1)/2.0;
//  cout << faktor << " ";
// alpha = (x - mu) / faktor;

//  cout << alpha << endl;

// // Adjustierung so dass mean(psi) = starwert


//   for (int i=1; i<= number_of_agegroups; i++){
//     theta[i-1]= theta[i-1] + alpha  * vdb * (double(i) - (noa + 1)/2);
//   }

//   for (int i=1; i<= number_of_periods; i++){
//     phi[i-1]= phi[i-1] - alpha * (double(i) - (nop + 1)/2);
//   }

//   for (int i=1; i<= number_of_cohorts; i++){
//     psi[i-1]= psi[i-1] + alpha * (double(i) - (noc + 1)/2) + alpha * faktor;

//   }

// //   // Test

// summe = 0.0;
// for (int g = 0; g < number_of_agegroups; g++){
//   summe += theta[g];
// }
//  cout << summe / noa << " ";

// summe = 0.0;
// for (int g = 0; g < number_of_periods; g++){
//   summe += phi[g];
// }
//  cout << summe / nop << " ";

// summe = 0.0;
// for (int g = 0; g < number_of_cohorts; g++){
//   summe += psi[g];
// }
//  cout << summe / noc << endl;



  // Leveladjustierung

// double sum = 0.0, adjust;

// for (int l = 0; l < number_of_agegroups;l++){
//   sum = sum + theta.get(l,0);
// }

// adjust = sum/double(number_of_agegroups);

// for (int g = 0; g < number_of_agegroups; g++){
//   theta.put(g,0, theta.get(g,0) - adjust);
// }

// for (int g = 0; g < number_of_cohorts; g++){
//   psi.put(g,0, psi.get(g,0) + adjust);
// }

// sum = 0.0;

// for (int l = 0; l < number_of_periods;l++){
//   sum = sum + phi.get(l,0);
// }

// adjust = sum/double(number_of_periods);

// for (int g = 0; g < number_of_periods; g++){
//   phi.put(g,0, phi.get(g,0) - adjust);
// }

// for (int g = 0; g < number_of_cohorts; g++){
//   psi.put(g,0, psi.get(g,0) + adjust);
// }

// double nny;

// nny=theta[1]-theta[0];

// for (int t1=0; t1<number_of_agegroups;t1++)
// 	{
// 	theta[t1]=theta[t1]-(t1*nny);
// 	}
// for (int t2=0; t2<number_of_periods; t2++)
// 	{
// 	phi[t2]=phi[t2]+(t2*nny);
// 	}
// for (int t3=0; t3<number_of_cohorts; t3++)
// 	{
// 	psi[t3]=psi[t3]-(t3*nny)+(nny*number_of_agegroups);
// 	}

// **Trendadjustierung**


double c;
c=2*(theta[0]+phi[0]+psi[0])/(noa-nop+noc-1.0);

for (int g=0; g < number_of_agegroups; g++){

  theta[g]=theta[g] + c*((g+1)-(number_of_agegroups+1)/2.0);
 }


for (int g=0; g < number_of_periods; g++){
  phi[g]=phi[g] - c*((g+1)-(number_of_periods+1)/2.0);
 }


for (int g=0; g < number_of_cohorts; g++){
  psi[g]=psi[g] + c*((g+1)-(number_of_cohorts+1)/2.0);
 }


// c=2*psi[0]/(noc-1.0);

// for (int g=0; g < number_of_agegroups; g++){

//   theta[g]=theta[g] + c*((g+1)-(number_of_agegroups+1)/2.0);
//  }


// for (int g=0; g < number_of_periods; g++){
//   phi[g]=phi[g] - c*((g+1)-(number_of_periods+1)/2.0);
//  }


// for (int g=0; g < number_of_cohorts; g++){
//   psi[g]=psi[g] + c*((g+1)-(number_of_cohorts+1)/2.0);
//  }

// c=2*phi[0]/(1.0-nop);

// for (int g=0; g < number_of_agegroups; g++){

//   theta[g]=theta[g] + c*((g+1)-(number_of_agegroups+1)/2.0);
//  }


// for (int g=0; g < number_of_periods; g++){
//   phi[g]=phi[g] - c*((g+1)-(number_of_periods+1)/2.0);
//  }


// for (int g=0; g < number_of_cohorts; g++){
//   psi[g]=psi[g] + c*((g+1)-(number_osinkt dadurch die Zahl der Senatoren auf unter neun, so w?hlen die Kandidaten der Liste einen Nachr?cker.f_cohorts+1)/2.0);
//  }

// c=2*theta[0]/(noa-1.0);

// for (int g=0; g < number_of_agegroups; g++){

//   theta[g]=theta[g] + c*((g+1)-(number_of_agegroups+1)/2.0);
//  }


// for (int g=0; g < number_of_periods; g++){
//   phi[g]=phi[g] - c*((g+1)-(number_of_periods+1)/2.0);
//  }


// for (int g=0; g < number_of_cohorts; g++){
//   psi[g]=psi[g] + c*((g+1)-(number_of_cohorts+1)/2.0);
//  }


// if (theta[0]>10e3)
//   {
//     for (int i=0; i< number_of_agegroups; i++)
//       {
// 	theta[i]=theta[i]/theta[number_of_agegroups-1];
//       }
//   }
// if (phi[0]>10e3)
//   {
//     for (int i=0; i< number_of_periods; i++)
//       {
// 	phi[i]=phi[i]/phi[number_of_periods-1];
//       }
//   }
// if (psi[0]>10e3)
//   {
//     for (int i=0; i< number_of_cohorts; i++)
//       {
// 	psi[i]=psi[i]/psi[number_of_cohorts-1];
//       }
//   }
return;
}


int max(int a, int b)
{
  if (a>b){return a;}
  else {return b;}
}

double ABS(double i)
{
double erg;
if (i>0.0)
  {
erg=i;
  }
else
  {
erg = -i;
  }
return erg;
}
