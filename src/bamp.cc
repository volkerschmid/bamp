#include <R.h>
#include <Rmath.h>
#include <time.h>
//#include "mystring.h"
//#include "l.h"
#include "update_my.h"
#include "zz_ksi.h"
#include "block.h"
#include "praez.h"
#include "zufall.h"
//#include "prog.h"
//#include "matrix.h"


extern "C"{

// cases: Faelle (y) in order=1(!)
// population: N
// blocks: A-P-C, 0 ohne, 1/2 RW 1/2, 3/4 RW 1/2 + unstrukturiert
// dims: nr agegroups, nr periods
// vdb vielfaches der breite (double!)
// mcnumbers: interations, burnin, thinning, tuning
// settings: prognose nein/ja, deviance nein/ja, overdispersion ja/nein
// ab: age_a, age_b, period_a, period_b, cohort_a, cohort_b, age_a2, age_b2, period_a2, period_b2, cohort_a2, cohort_b2, z_a, z_b
void bamp(int* cases, int* population, int* blocks, int* dims, double* vdb, int* mcnumbers, int* settings, double* ab,
          double* ttt, double* pph, double* pps, double* ttt2, double* pph2, double* pps2, double* ksisam,
          double* KK_delta, double* KK_th, double* KK_th2, double* KK_ph, double* KK_ph2, double* KK_ps, double* KK_ps2, double* my_off, 
          double* devsample, int* verbose0)
{
  int verbose=verbose0[0];
  if(verbose>=2){Rprintf("cc verbose level: %d\n\n", verbose);}
double age_a=ab[0];
double age_b=ab[1];
double period_a=ab[2];
double period_b=ab[3];
double cohort_a=ab[4];
double cohort_b=ab[5];
double age_a2=ab[6];
double age_b2=ab[7];
double period_a2=ab[8];
double period_b2=ab[9];
double cohort_a2=ab[10];
double cohort_b2=ab[11];
double z_a=ab[12];
double z_b=ab[13];
int age_block=blocks[0];
int period_block=blocks[1];
int cohort_block=blocks[2];
//int max_block=max(age_block,(max(period_block,cohort_block)));
int period_plus=0;
int cohort_plus=0;
int number_of_agegroups=dims[0];
int number_of_periods=dims[1];
double vielfaches_der_breite=vdb[0];
int number_of_predictions=dims[3];
int prognosis=settings[0];
int devi=settings[1];
int mode=settings[2];
int number_of_iterations=mcnumbers[0];
int burn_in=mcnumbers[1];
int abstand=mcnumbers[2];
int tun_konst=mcnumbers[3];
int number_of_cohorts=0;
int number_of_periods2=0;
int number_of_cohorts2=0;
int case_predict =0;
double temp1;
int max_number_of_ap_combinations;
double* period_data;
double* cohort_data;
int tuningcount=0;

//counters for sampling
int gen_c=0;
int ttt_c=0;
int pph_c=0;
int pps_c=0;
int ksisam_c=0;
int ttt2_c=0;
int pph2_c=0;
int pps2_c=0;

clock();
number_of_cohorts = (((int)vielfaches_der_breite)*(number_of_agegroups - 1))+number_of_periods;

if (prognosis == 0)
   {
     number_of_periods2=number_of_periods;
     number_of_cohorts2=number_of_cohorts;
     number_of_predictions = 0;
   }
 if (prognosis == 2)
   {
     number_of_periods2=number_of_periods+number_of_predictions;
     number_of_cohorts2=number_of_cohorts+number_of_predictions;
   }
 if (prognosis == 3)
 {
  prognosis = 1;
  case_predict = 1;
 }
 if (prognosis == 1)
   {
     number_of_periods2=number_of_periods;
     number_of_cohorts2=number_of_cohorts;
   }

 temp1 = ((number_of_periods2)/vielfaches_der_breite)+0.999;
 max_number_of_ap_combinations = (int)temp1;
 if (number_of_agegroups < max_number_of_ap_combinations)
   {
     max_number_of_ap_combinations = number_of_agegroups;
   }

 if (prognosis != 2 && devi==2)
   {
     devi=1;
   }

if (period_block==8 || period_block==9)
	{
	period_block = period_block-7;
	period_data = new double[number_of_periods+number_of_predictions];
	period_plus = 1;
	}
if (cohort_block==8 || cohort_block==9)
	{
	cohort_block = cohort_block-7;
	cohort_data = new double[number_of_cohorts+number_of_predictions];
	cohort_plus = 1;
	}

int** yes = new int*[number_of_agegroups];
int** no =new int*[number_of_agegroups];
double** akzeptanz = new double*[number_of_agegroups];
int** schalter = new int*[number_of_agegroups];

for (int i=0; i< number_of_agegroups; i++)
{
  yes[i] = new int[number_of_periods];
  no[i] =new int[number_of_periods];
  akzeptanz[i] = new double[number_of_periods];
  schalter[i] =new int[number_of_periods];

  for (int j=0; j< number_of_periods; j++)
    {
      yes[i][j]=0;
      no[i][j]=0;
      akzeptanz[i][j]=0.0;
      schalter[i][j]=0;
    }


}
// int need = 0;

// fuer blockupdate

double* ageQ = new double[number_of_agegroups*(age_block+1)];
double* perQ = new double[number_of_periods*(period_block+1)];
double* cohQ = new double[number_of_cohorts*(cohort_block+1)];



// ***** READING DATA *****
int** y = new int*[number_of_agegroups];
int** n =new int*[number_of_agegroups];

for (int i=0; i< number_of_agegroups; i++)
{
   y[i] = new int[number_of_periods2];
   if (case_predict==0){
   n[i] =new int[number_of_periods2];
   }
   else{ n[i] =new int[number_of_periods + number_of_predictions];}
}


int c=0;
for(int l=0; l<number_of_periods2;l++)
  {
    for(int m=0;m<number_of_agegroups;m++)
	  {
	    n[m][l] = population[c];
	    y[m][l] = cases[c];
	    c++;
	  }
  if(case_predict==1)
	 {
	   for (int l=number_of_periods2; l<(number_of_periods+number_of_predictions);l++)
	     {
	     for(int m=0;m<number_of_agegroups;m++)
		     {
		      n[m][l] = population[c];
		     }
	     }
	 }
}

// ***** CALCULATING STARTING-VALUE PSI *****

double bev_summe=0.0, lung_summe=0.0, startwert;


for (int i=0; i< number_of_agegroups; i++)
{
  for (int j=0; j< number_of_periods; j++)
    {
      bev_summe+=n[i][j];
      lung_summe+=y[i][j];
    }
}

startwert = lung_summe/bev_summe;
startwert = log(startwert/(1-startwert));


//  ***** SETTING VARIABLES *****

double* psi =new double[number_of_cohorts+number_of_predictions];
double* theta=new double[number_of_agegroups];
double* phi=new double[number_of_periods+number_of_predictions];
double* psi2 =new double[number_of_cohorts+number_of_predictions];
double* theta2=new double[number_of_agegroups];
double* phi2=new double[number_of_periods+number_of_predictions];
double* psitemp =new double[number_of_cohorts+number_of_predictions];
double* thetatemp=new double[number_of_agegroups];
double* phitemp=new double[number_of_periods+number_of_predictions];

for (int i=0; i<number_of_agegroups; i++)
  {theta[i]=0.0;theta2[i]=0.0;}
for (int i=0; i<number_of_periods2; i++)
  {
  phi[i]=0.0;phi2[i]=0.0;
  }
if (cohort_block>0)
{
for (int i=0; i<number_of_cohorts2; i++)
  {
  psi[i]=normal(0.0,1.0);psi2[i]=normal(0.0,1.0);
  }
}
if (cohort_block==0)
{
  for (int i=0; i<number_of_cohorts2; i++)
  {
    psi[i]=0.0;psi2[i]=0.0;
  }
}

double my=0.0;

double** z=new double*[number_of_agegroups];
for (int j=0;j<number_of_agegroups;j++)
{
z[j]=new double[number_of_periods];
}
double** ksi=new double*[number_of_agegroups];
for (int j=0;j<number_of_agegroups;j++)
{
ksi[j]=new double[number_of_periods+number_of_predictions];
}

//double z_prog;
ksi_berechnen(ksi, psi, vielfaches_der_breite, number_of_agegroups,number_of_periods);
double kappa = 10*age_a/age_b;
double lambda = 10*period_a/period_b;
double ny = 10*cohort_a/cohort_b;
double kappa2 = 10*age_a2/age_b2;
double lambda2 = 10*period_a2/period_b2;
double ny2 = 10*cohort_a2/cohort_b2;
double delta = 1.0;

//double mydach;
//double devsumme;
//double devianztemp;
int index=0;
//int tip = 0;

my=startwert;

// zaehle akzeptanz
int ja_age = 0;
int ja_period=0;
int ja_cohort=0;
int ja_my=0;

if (period_block==4||period_block==3) {delta=100;}
if (cohort_block==4||cohort_block==3) {delta=100;}
double itpercent=0.0;
double lastitpercent=-10.0;
int tempcounter=0;
int ididwritesomething=0;

for (int iteration=1; iteration<= number_of_iterations; iteration++)
{
  if (verbose>=3){Rprintf("Iteration: %d\n",iteration);ididwritesomething=0;}
  // progress bar
  if(verbose>0){
  //if (ididwritesomething==0){
    //for (tempcounter=0; tempcounter<=27; tempcounter++)
    //{Rprintf(" ");}
    //ididwritesomething=0;
  //}
  itpercent=floor(100.0*iteration/number_of_iterations);
  if (itpercent>lastitpercent)
  {
    for (tempcounter=0; tempcounter<=ididwritesomething; tempcounter++){Rprintf("\b");}
    for (tempcounter=0; tempcounter <= (itpercent/2.5); tempcounter++){Rprintf("=");}
    if (itpercent<100){Rprintf(">");}else{Rprintf(" ");}
    for (;tempcounter<=40; tempcounter++){Rprintf(" ");}
    if (itpercent<100){Rprintf(" ");}
    if (itpercent<10){Rprintf(" ");}
    Rprintf("%d %% ",(int)itpercent);
    lastitpercent=itpercent;
    ididwritesomething=47;
  }
}

    // ** Gibbs-Steps **

  if (mode==0)
    {

    
    if (age_block==1 || age_block==2)
	    {
	      blocken(1, age_block, kappa, delta, number_of_agegroups, number_of_periods, my, theta,  phi,  psi,  ageQ,  thetatemp,  vielfaches_der_breite, n, y, ja_age );
	    }
      if (age_block==3||age_block==4)
	    {
	      blocken2(1, age_block, kappa, kappa2, number_of_agegroups, number_of_periods, my, theta, theta2,  phi,  psi,  ageQ,  thetatemp,  vielfaches_der_breite, n, y, ja_age );
	    }

      if (period_block==1 || period_block==2)
	    {
	      blocken(2, period_block, lambda, delta, number_of_periods, number_of_agegroups, my, phi,  theta,  psi,  perQ,  phitemp,  vielfaches_der_breite, n, y, ja_period );
	    }
      if (period_block==3 || period_block==4)
	    {
	      blocken2(2, period_block, lambda, lambda2, number_of_periods, number_of_agegroups, my, phi, phi2,  theta,  psi,  perQ,  phitemp,  vielfaches_der_breite, n, y, ja_period );
	    }

      if (cohort_block==1 || cohort_block==2)
	    {
	      blocken(-number_of_agegroups, cohort_block, ny, delta, number_of_cohorts, number_of_periods, my, psi,  phi,  theta,  cohQ,  psitemp, vielfaches_der_breite, n, y, ja_cohort );
	    }
      if (cohort_block==3 || cohort_block==3)
	    {
	      blocken2(-number_of_agegroups, cohort_block, ny, ny2, number_of_cohorts, number_of_periods, my, psi, psi2,  phi,  theta, cohQ,  psitemp, vielfaches_der_breite, n, y, ja_cohort );
	    }

      update_my_mh(my, ksi, theta, phi,  psi, vielfaches_der_breite, number_of_agegroups, number_of_periods, lung_summe,  n,  y, ja_my);
      
      for (int i=0; i<number_of_agegroups; i++)
	    {
 	      for (int j=0; j<number_of_periods; j++)
 	        {
 	          ksi[i][j]=my+theta[i]+phi[j]+psi[coh(i,j, number_of_agegroups, vielfaches_der_breite)-1];
 	        }
 	    }

    } //mode==0

  if (mode==1)
    {
      if (age_block==1 || age_block==2)
    	{
	      blockupdate(1, age_block, kappa, delta, number_of_agegroups, number_of_periods,  ksi, my, theta,  phi,  psi,  ageQ,  thetatemp,  vielfaches_der_breite);
	    }
      if (age_block==3 || age_block==4)
	    {
	      blockupdate_a2(1, age_block-2, kappa, kappa2, delta, number_of_agegroups, number_of_periods,  ksi, my, theta,  theta2,  phi,  psi, thetatemp,  vielfaches_der_breite);
	    }
      if (period_plus == 0)
    	{
	      if (period_block==1 || period_block==2)
	        {
	          blockupdate(2, period_block, lambda, delta, number_of_periods, number_of_agegroups,  ksi, my, phi,  theta,  psi,  perQ,  phitemp,  vielfaches_der_breite);
    	    }
	     if (period_block==3||period_block==4)
	        {
	          blockupdate_a2(2, period_block-2, lambda, lambda2, delta, number_of_periods, number_of_agegroups,  ksi, my, phi,  phi2,  theta,  psi, phitemp,  vielfaches_der_breite);
	        }
	    }
      else
	    {
	      if (period_block==1 || period_block==2)
	        {
	          blockupdateplus(2, period_block, lambda, delta, number_of_periods, number_of_agegroups,  ksi, my, phi,  theta,  psi,  perQ,  phitemp,  vielfaches_der_breite, period_data);
	        }
	    }

      if (cohort_plus == 0)
	    {
	      if (cohort_block==1 || cohort_block==2)
	        {
	          blockupdate(-number_of_agegroups, cohort_block, ny, delta, number_of_cohorts, number_of_periods,  ksi, my, psi,  phi,  theta,  cohQ,  psitemp, vielfaches_der_breite);
	        }
	      if (cohort_block==3||cohort_block==4)
	        {
	          blockupdate_a2(-number_of_agegroups, cohort_block-2, ny,ny2,  delta, number_of_cohorts, number_of_periods,  ksi, my, psi, psi2, phi,  theta,  psitemp, vielfaches_der_breite);
	        }
	    }
      else
	    {
	      if (cohort_block==1 || cohort_block==2)
	        {
	          blockupdateplus(-number_of_agegroups, cohort_block, ny, delta, number_of_cohorts, number_of_periods,  ksi, my, psi,  phi,  theta,  cohQ,  psitemp, vielfaches_der_breite, cohort_data);
	        }
	     }

      ZZ_aus_fc_von_ksi0(my, theta,  phi,  psi,  ksi, delta,  number_of_agegroups,  number_of_periods, vielfaches_der_breite,  y,  n,yes,no,schalter);
      z_aus_ksi_berechnen(z, my, ksi, theta, phi, psi, vielfaches_der_breite,number_of_agegroups,number_of_periods) ;
      update_my_1(my, ksi,  theta, phi,  psi,  vielfaches_der_breite,  number_of_agegroups, number_of_periods, delta);
    } //mode==1


  // Updaten der Hyperparameter
      if (period_plus==1)
	{
	  for (int i=0; i<number_of_periods;i++)
	    {
	      phi[i]=phi[i]/period_data[i];
	    }
	}
      if (cohort_plus==1)
	{
	  for (int i=0; i<number_of_cohorts;i++)
	    {
	      psi[i]=psi[i]/cohort_data[i];
	    }
	}

  if (age_block==3 || age_block==4)
    {
      for (int i=0; i< number_of_agegroups; i++)
	{
	  theta[i]=theta[i]-theta2[i];
	}
    }
  if (period_block==3 || period_block==4)
    {
      for (int i=0; i< number_of_periods; i++)
	{
	  phi[i]=phi[i]-phi2[i];
	}
    }
  if (cohort_block==3 || cohort_block==4)
    {
      for (int i=0; i< number_of_cohorts; i++)
	{
	  psi[i]=psi[i]-psi2[i];
	}
    }

  if (age_block!=0)
    {
      kappa = hyper(age_block, theta, age_a, age_b,number_of_agegroups);
    }

  if (period_block!=0)
    {
      lambda = hyper(period_block, phi, period_a, period_b, number_of_periods);
    }

  if (cohort_block!=0)
    {
     ny = hyper(cohort_block, psi, cohort_a, cohort_b,number_of_cohorts);
    }


if (age_block==3 || age_block==4)
  {
    for (int i=0; i< number_of_agegroups; i++)
      {
	theta[i]=theta[i]+theta2[i];
	kappa2 = hyper2(theta2, age_a2, age_b2, number_of_agegroups);
      }
  }
if (period_block==3 || period_block==4)
  {
    for (int i=0; i< number_of_periods; i++)
      {
	phi[i]=phi[i]+phi2[i];
	lambda2 = hyper2(phi2, period_a2, period_b2, number_of_periods);
      }
  }
if (cohort_block==3 || cohort_block==4)
  {
    for (int i=0; i< number_of_cohorts; i++)
      {
	psi[i]=psi[i]+psi2[i];
	ny2 =  hyper2(psi2, cohort_a2, cohort_b2, number_of_cohorts);
      }
  }
      if (period_plus==1)
	{
	  for (int i=0; i<number_of_periods;i++)
	    {
	      phi[i]=phi[i]*period_data[i];
	    }
	}
      if (cohort_plus==1)
	{
	  for (int i=0; i<number_of_cohorts;i++)
	    {
	      psi[i]=psi[i]*cohort_data[i];
	    }
	}

if (mode==1)
  {
      delta = delta_berechnen(z, z_a, z_b, number_of_agegroups, number_of_periods);
  }

      //
      // TUNING 
      //
      // CHECK ACCEPTANCE RATE


 if (iteration == tun_konst && mode==1){
	   int back=0;

	   for (int l=0;l<number_of_agegroups;l++)
	     {
	       for (int m=0;m<number_of_periods;m++)
		      {
		      akzeptanz[l][m]=(double)(yes[l][m])/(double)(yes[l][m]+no[l][m]);
		      if (akzeptanz[l][m]==0)
		      {
		         back=1;
		         schalter[l][m]=1-schalter[l][m];
		      }

		    }
	     }
	   if (back==1)
	     {
	       iteration=0;
	       if(verbose>=1){
	         for (tempcounter=0; tempcounter<=ididwritesomething; tempcounter++){Rprintf("\b");}
	         
	          Rprintf("Tuning is necessary. Restarting iterations.");
	          ididwritesomething=42;
	          lastitpercent=0.0;
	          tuningcount++;
	       }
	       for (int l=0;l<number_of_agegroups;l++)
	       {
	         for (int m=0;m<number_of_periods;m++)
	         {
	           yes[l][m]=0;
	           no[l][m]=0;
	           akzeptanz[l][m]=0.0;
	           z[l][m]=0.0;
	         }
	       }
	       delta=1000;
	     }
 }


 if (iteration == tun_konst && mode==0){
	   int back=0;

	   if ((double)ja_my/(double)tun_konst<.3){back=1; my=0;}
	   if(age_block>0)if ((double)ja_age/(double)tun_konst<.3){back=1;}
	   if(period_block>0)if ((double)ja_period/(double)tun_konst<.3){back=1;}
	   if(cohort_block>0)if ((double)ja_cohort/(double)tun_konst<.3){back=1;}

	   if (back==1)
	     {
	       iteration=0;
	        if(verbose>=1){
	          for (tempcounter=0; tempcounter<=ididwritesomething; tempcounter++){Rprintf("\b");}}
	          
            Rprintf("Acceptance rate low. Restarting iterations.\n\n");
	          lastitpercent=0.0;
	          ididwritesomething=0;
	          kappa = 10*age_a/age_b;
	          lambda = 10*period_a/period_b;
	          ny = 10*cohort_a/cohort_b;
	          kappa2 = 10*age_a2/age_b2;
	          lambda2 = 10*period_a2/period_b2;
	          ny2 = 10*cohort_a2/cohort_b2;
	          for (int i=0; i<number_of_agegroups; i++){theta[i]=0.0;theta2[i]=0.0;}
	          for (int i=0; i<number_of_periods2; i++){phi[i]=0.0;phi2[i]=0.0;}}
	          if (cohort_block>0){for (int i=0; i<number_of_cohorts2; i++){psi[i]=normal(0.0,1.0);psi2[i]=normal(0.0,1.0);}}
	          if (cohort_block==0){for (int i=0; i<number_of_cohorts2; i++){psi[i]=0.0;psi2[i]=0.0;}}
	          my=startwert+normal(0,1.0);

	   if (verbose>=2)
	   {
	     Rprintf("%d ",ja_my);
	     Rprintf("%d ",ja_age);
	     Rprintf("%d ",ja_period);
	     Rprintf("%d\n\n",ja_cohort);
	   }

	  ja_age=0;
	  ja_period=0;
	  ja_cohort=0;
	  ja_my=0;

 }


      // save samples

 if(fmod((double)iteration,2.0*(double)abstand)==0.0 && iteration>tun_konst)
	{

	  // zentriere(my, psi, number_of_cohorts);

	  ja_age=0;
	  ja_period=0;
	  ja_cohort=0;
	  ja_my=0;

	 }

 if(iteration>burn_in && fmod((double)iteration,(double)abstand)==0.0)
	{
	  if (period_plus==1)
	    {
	      for (int i=0; i< number_of_periods; i++){phi[i]=phi[i]/period_data[i];}
	    }

	  if (cohort_plus==1)
	    {
	      for (int i=0; i< number_of_cohorts; i++){psi[i]=psi[i]/cohort_data[i];}
	    }




	      int ksi_cc=0;
	      double temp1=double(ksisam_c);
	      double temp2=temp1/(1.0+temp1);
	      temp1=1.0/(temp1+1.0);
	      for (int h=0; h< number_of_agegroups;h++){
	        for (int j=0; j< number_of_periods;j++){
	          ksisam[ksi_cc]=temp2*ksisam[ksi_cc]+temp1*ksi[h][j];
	          ksi_cc++;
	          //Rprintf("ksisam_c:%d",(int)ksisam_c);
	         }
	        }
	      ksisam_c++;
	      
	  for (int h=0; h< number_of_agegroups;h++)
	    {
	      ttt[ttt_c]=theta[h];
  	    ttt_c++;
	    }

	  	if (age_block==3 || age_block==4)
	    {
	      for (int h=0; h< number_of_agegroups;h++)
		    {
		      ttt2[ttt2_c]=theta2[h];
	        ttt2_c++;
    		}
	    }
	  if (prognosis==0)
	    {
      for (int h=0; h< number_of_periods;h++)
    		{
	  	  if (period_plus==1)
		      {
		        pph[pph_c]=(phi[h]*period_data[h]);
	  	      pph_c++;
		    }
		  else
		    {
		      pph[pph_c]=phi[h];
		      pph_c++;
		    }
		}
    if (period_block==3 || period_block==4)
		{
		  for (int h=0; h< number_of_periods;h++)
		    {
		      if (period_plus==1)
	    		{
			      pph2[pph2_c]=(phi2[h]*period_data[h]);
		        pph2_c++;
			    }
		      else
			    {
			      pph2[pph2_c]=phi2[h];
		        pph2_c++;
			    }
		    }
		}

  for (int h=0; h< number_of_cohorts;h++)
		{
		  if (cohort_plus==1)
		    {
		      pps[pps_c]=(psi[h]*cohort_data[h]);
		    }
		  else
		    {
		      pps[pps_c]=psi[h];
		    }
		  pps_c++;
		}
  if (cohort_block==3 || cohort_block==4)
		{
		  for (int h=0; h< number_of_cohorts;h++)
		    {
		      if (cohort_plus==1)
    			{
		    	  pps2[pps2_c]=(psi2[h]*cohort_data[h]);
			    }
		      else
    			{
		    	  pps2[pps2_c]=psi2[h];
			    }
		      pps2_c++;
		    }
		}

  KK_th[gen_c] = kappa;
  if (age_block==3 || age_block==4) KK_th2[gen_c] = kappa2;
  KK_ph[gen_c] = lambda;
	if (period_block==3 || period_block==4) KK_ph2[gen_c] = lambda2;
	KK_ps[gen_c] = ny;
	if (cohort_block==3 || cohort_block==4) KK_ps2[gen_c] = ny2;
	if (mode==1) KK_delta[gen_c] = delta;
  my_off[gen_c] = my;

	   if (period_plus==1)
	     {
	       for (int i=0; i< number_of_periods; i++){phi[i]=phi[i]*period_data[i];}
	     }

	   if (cohort_plus==1)
	     {
	       for (int i=0; i< number_of_cohorts; i++){psi[i]=psi[i]*cohort_data[i];}
	     }

	   index += 1;
	}


	  double temp=0.0;
	  devsample[gen_c]=0.0;
	  for (int l6=0;l6<number_of_agegroups;l6++)
	  {
	    for(int m=0;m<number_of_periods;m++)
	    {						       
  	   temp = exp(ksi[l6][m])/(1+exp(ksi[l6][m])); // pr
	     temp =temp*n[l6][m]; //mydach
	      if(y[l6][m]==0.0)
	      {   
	        devsample[gen_c]+=2*((n[l6][m]-y[l6][m])*log((n[l6][m]-y[l6][m])/(n[l6][m]-temp)));
	      }
	      else
	      {
	        devsample[gen_c]+=2*(y[l6][m]*log(y[l6][m]/temp)+(n[l6][m]-y[l6][m])*log((n[l6][m]-y[l6][m])/(n[l6][m]-temp)));
	      }
	    }	
	  }	
	  gen_c++;
	  
	}

 }//Ende Iterationsschleife
  if (verbose>=1){Rprintf("\n\n");}

  return;
}
}
