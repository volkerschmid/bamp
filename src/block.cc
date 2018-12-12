#include <R.h>
#include <iostream>
#include "matrix.h"
#include "l.h"
#include "mxs.h"
#include "zufall.h"


//Berechnet Matrix Q (f?r age, period)
void berechneQ(double* temp, int age_block, double kappa,int noa,int nop,double delta)
{

if (age_block==1)
  {
    int index=0;
    temp[index]=kappa+delta*nop;
    index++;
    temp[index]=-kappa;
    index++;
    for (int i=1; i<noa-1; i++)
      {
	temp[index]=2*kappa+delta*nop;
	index++;
	temp[index]=-kappa;
	index++;
      }
    temp[index]=kappa+delta*nop;
  }
if (age_block==2)
  {
    int index=0;
    temp[index]=kappa+delta*nop;
    index++;
    temp[index]=-2*kappa;
    index++;
    temp[index]=kappa;
    index++;
    temp[index]=5*kappa+delta*nop;
    index++;
    temp[index]=-4*kappa;
    index++;
    temp[index]=kappa;
    index++;
    for (int i=2; i<noa-2; i++)
      {
	temp[index]=6*kappa+delta*nop;
	index++;
	temp[index]=-4*kappa;
	index++;
	temp[index]=kappa;
	index++;
      }
    temp[index]=5*kappa+delta*nop;
    index++;
    temp[index]=-2*kappa;
    index++;
    //temp[index]=0;
    index++;
    temp[index]=kappa+delta*nop;
  }
return;
}

void berechneQcohort(int* n, double* temp, int age_block, double kappa,int noa,double delta)
{

if (age_block==1)
  {
    int index=0;
    temp[index]=kappa+delta*n[0];
    index++;
    temp[index]=-kappa;
    index++;
    for (int i=1; i<noa-1; i++)
      {
	temp[index]=2*kappa+delta*n[i];
	index++;
	temp[index]=-kappa;
	index++;
      }
    temp[index]=kappa+delta*n[noa-1];
  }
if (age_block==2)
  {
    int index=0;
    temp[index]=kappa+delta*n[0];
    index++;
    temp[index]=-2*kappa;
    index++;
    temp[index]=kappa;
    index++;
    temp[index]=5*kappa+delta*n[1];
    index++;
    temp[index]=-4*kappa;
    index++;
    temp[index]=kappa;
    index++;
    for (int i=2; i<noa-2; i++)
      {
	temp[index]=6*kappa+delta*n[i];
	index++;
	temp[index]=-4*kappa;
	index++;
	temp[index]=kappa;
	index++;
      }
    temp[index]=5*kappa+delta*n[noa-2];
    index++;
    temp[index]=-2*kappa;
    index++;
   // temp[index]=0;
    index++;
    temp[index]=kappa+delta*n[noa-1];
    index++;
   // temp[index]=0;
    index++;
    //temp[index]=0;
  }
return;
}
void berechneB(int swit, double* temp, double** ksi, double my,  double* phi, double* psi,int noa, int nop, int vdb, double delta)
{

switch(swit)
  {

  case 1:

for (int i=0; i< noa; i++)
  {
    temp[i]=0.0;
    for (int j=0; j< nop; j++)
      {
	temp[i]=temp[i]+ksi[i][j]-my-phi[j]-psi[coh(i,j,noa,vdb)-1];
      }
    temp[i]=temp[i]*delta;
  }
break;

  case 2:

for (int i=0; i< noa; i++)
  {
    temp[i]=0.0;
    for (int j=0; j< nop; j++)
      {
	temp[i]=temp[i]+ksi[j][i]-my-phi[j]-psi[coh(j,i,nop,vdb)-1];
      }
    temp[i]=temp[i]*delta;
  }
break;

}
return;
}
void berechneQplus(double* temp, int age_block, double kappa,int noa,int nop,double delta, double* data)
{

if (age_block==1)
  {
    int index=0;
    temp[index]=kappa+delta*nop*data[0]*data[0];
    index++;
    temp[index]=-kappa;
    index++;
    for (int i=1; i<noa-1; i++)
      {
	temp[index]=2*kappa+delta*nop*data[i]*data[i];
	index++;
	temp[index]=-kappa;
	index++;
      }
    temp[index]=kappa+delta*nop*data[noa-1]*data[noa-1];
  }
if (age_block==2)
  {
    int index=0;
    temp[index]=kappa+delta*nop*data[0]*data[0];
    index++;
    temp[index]=-2*kappa;
    index++;
    temp[index]=kappa;
    index++;
    temp[index]=5*kappa+delta*nop*data[1]*data[1];
    index++;
    temp[index]=-4*kappa;
    index++;
    temp[index]=kappa;
    index++;
    for (int i=2; i<noa-2; i++)
      {
	temp[index]=6*kappa+delta*nop*data[i]*data[i];
	index++;
	temp[index]=-4*kappa;
	index++;
	temp[index]=kappa;
	index++;
      }
    temp[index]=5*kappa+delta*nop*data[noa-2]*data[noa-2];
    index++;
    temp[index]=-2*kappa;
    index++;
    //temp[index]=0;
    index++;
    temp[index]=kappa+delta*nop*data[noa-1]*data[noa-1];
  }
return;
}

void berechneQcohortplus(int* n, double* temp, int age_block, double kappa,int noa,double delta, double* data)
{

if (age_block==1)
  {
    int index=0;
    temp[index]=kappa+delta*n[0]*data[0]*data[0];
    index++;
    temp[index]=-kappa;
    index++;
    for (int i=1; i<noa-1; i++)
      {
	temp[index]=2*kappa+delta*n[i]*data[i]*data[i];
	index++;
	temp[index]=-kappa;
	index++;
      }
    temp[index]=kappa+delta*n[noa-1]*data[noa-1]*data[noa-1];
  }
if (age_block==2)
  {
    int index=0;
    temp[index]=kappa+delta*n[0]*data[0]*data[0];
    index++;
    temp[index]=-2*kappa;
    index++;
    temp[index]=kappa;
    index++;
    temp[index]=5*kappa+delta*n[1]*data[1]*data[1];
    index++;
    temp[index]=-4*kappa;
    index++;
    temp[index]=kappa;
    index++;
    for (int i=2; i<noa-2; i++)
      {
	temp[index]=6*kappa+delta*n[i]*data[i]*data[i];
	index++;
	temp[index]=-4*kappa;
	index++;
	temp[index]=kappa;
	index++;
      }
    temp[index]=5*kappa+delta*n[noa-2]*data[noa-2]*data[noa-2];
    index++;
    temp[index]=-2*kappa;
    index++;
   // temp[index]=0;
    index++;
    temp[index]=kappa+delta*n[noa-1]*data[noa-1]*data[noa-1];
    index++;
   // temp[index]=0;
    index++;
    //temp[index]=0;
  }
return;
}
void berechneBplus(int swit, double* temp, double** ksi, double my, double* phi, double* psi,int noa, int nop, int vdb, double delta, double* data)
{

switch(swit)
  {

  case 1:

for (int i=0; i< noa; i++)
  {
    temp[i]=0.0;
    for (int j=0; j< nop; j++)
      {
	temp[i]=temp[i]+ksi[i][j]-my-phi[j]-psi[coh(i,j,noa,vdb)-1];
      }
    temp[i]=temp[i]*delta*data[i];
  }
break;

  case 2:

for (int i=0; i< noa; i++)
  {
    temp[i]=0.0;
    for (int j=0; j< nop; j++)
      {
	temp[i]=temp[i]+ksi[j][i]-my-phi[j]-psi[coh(j,i,nop,vdb)-1];
      }
    temp[i]=temp[i]*delta*data[i];
  }
break;

}
return;
}

void berechneB_S(int swit, double* temp, double*** ksi, double* phi, double* psi, double* alpha, int noa, int nop, int vdb, double delta, int nor)
{

switch(swit)
  {

  case 1:

for (int i=0; i< noa; i++)
  {
    temp[i]=0.0;
    for (int j=0; j< nop; j++)
      {
	for (int k=0; k< nor; k++)
	  {
	    temp[i]=temp[i]+ksi[i][j][k]-phi[j]-psi[coh(i,j,noa,vdb)-1]-alpha[k];
	  }
      }
    temp[i]=temp[i]*delta;
  }
break;

  case 2:

for (int i=0; i< noa; i++)
  {
    temp[i]=0.0;
    for (int j=0; j< nop; j++)
      {
	for (int k=0; k< nor; k++)
	  {
	    temp[i]=temp[i]+ksi[j][i][k]-phi[j]-psi[coh(j,i,nop,vdb)-1]-alpha[k];
	  }
      }
    temp[i]=temp[i]*delta;
  }
break;

  case 4:

    for (int i=0; i<noa; i++)
      {
	temp[i]=0.0;
	for (int j=0; j< nop; j++)
	  {
	    for (int k=0; k< nor; k++)
	      {
		temp[i]=temp[i]+ksi[k][j][i]-phi[j]-psi[coh(k,j,nop,vdb)-1]-alpha[k];
	      }
	  }

	temp[i]=temp[i]*delta;
      }
break;

}
return;
}
void berechneBcohort(int* n, double* temp, double** ksi, double my, double* phi, double* theta,int noa, int nop, int vdb, double delta, int noc)
{

for (int i=0; i< noc; i++)
  {
    temp[i]=0.0;
    n[i]=0;
  }

for (int i=0; i< noa; i++)
  {
    for (int j=0; j< nop; j++)
      {
	temp[coh(i,j,noa,vdb)-1]=temp[coh(i,j,noa,vdb)-1]+ksi[i][j]-my-phi[j]-theta[i];
	n[coh(i,j,noa,vdb)-1]=n[coh(i,j,noa,vdb)-1]+1;
      }

  }

for (int i=0; i< noc; i++)
  {
temp[i]=temp[i]*delta;
  }
return;
}
void berechneBcohortplus(int* n, double* temp, double** ksi, double my, double* phi, double* theta,int noa, int nop, int vdb, double delta, int noc, double* data)
{

for (int i=0; i< noc; i++)
  {
    temp[i]=0.0;
    n[i]=0;
  }

for (int i=0; i< noa; i++)
  {
    for (int j=0; j< nop; j++)
      {
	temp[coh(i,j,noa,vdb)-1]=temp[coh(i,j,noa,vdb)-1]+ksi[i][j]-my-phi[j]-theta[i];
	n[coh(i,j,noa,vdb)-1]=n[coh(i,j,noa,vdb)-1]+1;
      }

  }

for (int i=0; i< noc; i++)
  {
temp[i]=temp[i]*delta*data[i];
  }
return;
}

//Berechnet Vektor b f?r Kohorte in SAPC
void berechneBcohort_S(int* n, double* temp, double*** ksi, double* phi, double* theta, double* alpha, int noa, int nop, int vdb, double delta, int noc ,int nor)
{

for (int i=0; i< noc; i++)
  {
    temp[i]=0.0;
    n[i]=0;
  }

for (int i=0; i< noa; i++)
  {
    for (int j=0; j< nop; j++)
      {
	for (int k=0; k<nor; k++)
	  {
	    temp[coh(i,j,noa,vdb)-1]+=ksi[i][j][k]-phi[j]-theta[i]-alpha[k];
	    n[coh(i,j,noa,vdb)-1]=n[coh(i,j,noa,vdb)-1]+1;
	  }
      }

  }

for (int i=0; i< noc; i++)
  {
temp[i]=temp[i]*delta;
n[i]=n[i]*nor;
  }
return;
}

//Erzeugt Normalverteilten Zufallsvektor der L?nge noa
void gausssample(double* temp, int noa)
{

for (int i=0; i< noa; i++)
  {
    temp[i]=normal(0,1);
  }

return;
}

void berechneQspace(double tau, double* Qspace, double deltamalJI, int bandw, int nor)
{
int index=0;
for (int i=0;i<nor; i++)
  {
    Qspace[index]=tau*Qspace[index]+deltamalJI;
    index++;
    for (int j=1; j<bandw;j++)
	{
	  Qspace[index]=tau*Qspace[index];
	  index++;
	}
  }
return;
}

void blockupdate(int swit, int age_block, double kappa, double delta, int noa, int nop, double** ksi, double &my, double* theta, double* phi, double* psi, double* ageQ, double* thetatemp, int vdb)
{
int bandw=age_block;
int* nwerte;
if (swit<0)
  {
    nwerte=new int[noa];
    berechneBcohort(nwerte, thetatemp,ksi, my, phi, psi, -swit, nop,vdb,delta, noa);
    berechneQcohort(nwerte, ageQ, age_block, kappa, noa, delta);
    delete[] nwerte;
  }
else
  {
    berechneB(swit, thetatemp,ksi, my, phi, psi, noa, nop,vdb,delta);
    berechneQ(ageQ, age_block, kappa, noa, nop, delta);
  }


ageQ = cholesky(noa, ageQ, bandw);
double* L=new double[(age_block+1)*noa];
for (int i=0; i< (age_block+1)*noa; i++)
  {
    L[i]=ageQ[i];
  }

//Lv=b
loese2(L, thetatemp, noa, bandw);
//L^T my=v
loese(L,thetatemp, noa, bandw);

//z
gausssample(theta,noa);

//L^T y=z
loese(L,theta,noa,bandw);

//x=my + y
double sum=0;
for (int i=0; i< noa; i++)
  {
    theta[i]=theta[i]+thetatemp[i];
    sum+=theta[i];
  }

//Zentrieren
sum=sum/(double)noa;

for (int i=0; i< noa; i++)
  {
    theta[i]=theta[i]-sum;
  }
//my=my+sum;


delete[] L;


return;

}

void blockupdateplus(int swit, int age_block, double kappa, double delta, int noa, int nop, double** ksi, double &my, double* theta, double* phi, double* psi, double* ageQ, double* thetatemp, int vdb, double* data)
{
int bandw=age_block;
int* nwerte;
for (int i=0; i< noa; i++)
{
	theta[i]=theta[i]/data[i];
}

if (swit<0)
  {
    nwerte=new int[noa];
    berechneBcohortplus(nwerte, thetatemp,ksi, my, phi, psi, -swit, nop,vdb,delta, noa, data);
    berechneQcohortplus(nwerte, ageQ, age_block, kappa, noa, delta, data);
    delete[] nwerte;
  }
else
  {
    berechneBplus(swit, thetatemp,ksi, my, phi, psi, noa, nop,vdb,delta, data);
    berechneQplus(ageQ, age_block, kappa, noa, nop, delta, data);
  }


ageQ = cholesky(noa, ageQ, bandw);
double* L=new double[(age_block+1)*noa];
for (int i=0; i< (age_block+1)*noa; i++)
  {
    L[i]=ageQ[i];
  }

//Lv=b
loese2(L, thetatemp, noa, bandw);
//L^T my=v
loese(L,thetatemp, noa, bandw);

//z
gausssample(theta,noa);

//L^T y=z
loese(L,theta,noa,bandw);

//x=my + y
for (int i=0; i< noa; i++)
  {
    theta[i]=theta[i]+thetatemp[i];
  }
delete[] L;

for (int i=0; i< noa; i++)
{
	theta[i]=theta[i]*data[i];
}

double sum=0;
for (int i=0; i< noa; i++)
  {
    sum+=theta[i];
  }

//Zentrieren
sum=sum/(double)noa;
for (int i=0; i< noa; i++)
  {
    theta[i]=theta[i]-sum;
  }
// my=my+sum;


return;


}




void blockupdate_S(int swit, int age_block, double kappa, double delta, int noa, int nop, double*** ksi, double* theta, double* phi, double* psi, double* alpha, double* ageQ, double* thetatemp, int vdb, int nor)
{

int bandw=age_block;
int* nwerte;
if (swit<0)
  {
    nwerte=new int[noa];
    berechneBcohort_S(nwerte, thetatemp,ksi, phi, psi,alpha, -swit, nop,vdb,delta, noa, nor);
    berechneQcohort(nwerte, ageQ, age_block, kappa, noa, delta);
    delete[] nwerte;
  }
else
  {
    berechneB_S(swit, thetatemp,ksi, phi, psi, alpha, noa, nop,vdb,delta, nor);
    if (swit==4)
      {
	berechneQspace(kappa, ageQ, delta*nop*nor, age_block, noa);
      }
    else
      {
	berechneQ(ageQ, age_block, kappa, noa, nop*nor, delta);
      }
  }

ageQ = cholesky(noa, ageQ, bandw);
double* L=new double[(age_block+1)*noa];
for (int i=0; i< (age_block+1)*noa; i++)
  {
    L[i]=ageQ[i];
  }

loese2(L, thetatemp, noa, bandw);

loese(L,thetatemp, noa, bandw);

gausssample(theta,noa);
loese(L,theta,noa,bandw);
for (int i=0; i< noa; i++)
  {
    theta[i]=theta[i]+thetatemp[i];
  }
delete[] L;

return;
}

void bedinge(int age_block, int noa, double* theta, double* L, double* thetatemp, int k, double* A, double* b)
{

int bandw=age_block;
double* V;
V=new double[k*noa];
// double* Atemp;
// Atemp=new double[k*noa];
// for (int i=0; i< noa*k; i++)
//   {
//     Atemp[i]=A[i];
//   }


for (int i=0; i<k; i++)
  {
    // erzeuge A^T_i
    double* Atransi=new double[noa];
   for (int j=0; j<noa; j++)
     {
       Atransi[j]=A[(i*noa)+j];
// 	cout <<  Atransi[j]<< " ";
      }
    // Lu=A^T_i
    loese2(L, Atransi, noa, bandw);
    // L^T*V_i=u
    loese(L, Atransi, noa, bandw);

    //erzeuge V
    //int jj=0;
    // for (int
    for (int j=0; j<noa; j++)
      {

	V[j*k+i]=Atransi[j];
// 	cout <<  Atransi[j]<< " ";
      }
    delete[] Atransi;
//      cout << endl;
  }
//   cout << "--" << endl;
// mxschreibe(V,noa,k);
double* y=new double[k];
//Ax
multipliziere(A,theta,k,noa,1,y);
//Ax-b
for (int i=0; i<k; i++)
  {
    y[i]=y[i]-b[i];
  }

double* AV=new double[k*k];

//A*V
multipliziere(A,V,k,noa,k,AV);

//(AV)^-1
invers(AV,k);

double* W=new double[noa*k];
//W=V(AV)^-1
multipliziere(V,AV,noa,k,k,W);

//W(Ax-b)
multipliziere(W,y,noa,k,1,thetatemp);

//berechne x-W(Ax-b)
for (int i=0; i< noa; i++)
  {
    theta[i]=theta[i]-thetatemp[i];
  }

delete[] V;
delete[] y;
delete[] W;
delete[] AV;

return;

}


#include <iostream>
#include "matrix.h"
#include "l.h"
#include "mxs.h"
#include "zufall.h"

#include <fstream>
using namespace std;

void gausssample(double* temp, int noa);


//Berechnet Matrix Q (f?r age, period)
void berechneQ(double* temp, int age_block, double kappa,int noa,int nop,double delta);



void berechneQcohort(int* n, double* temp, int age_block, double kappa,int noa,double delta);




double taylor1(double x)
{
  return ( (exp(x)) /  ( 1+exp(x) ) );
}
double taylor2(double x)
{
  double temp=taylor1(x);
  return (temp-temp*temp);
}
double taylor0(double x)
{
  return (log ( 1+exp(x) ) );
}
void MausQtheta(double* Q, int b, int** daten_n, double* theta, double* phi, double* psi, int noa, int nop, int vdb, double my)
{

  for (int i=0; i< noa; i++)
  {
    for (int t=0; t < nop; t++)
    {
      Q[i*b] = Q[i*b] + (daten_n[i][t]*taylor2(theta[i] + phi[t] + psi[coh(i,t,noa,vdb)-1] + my));
    }
  }

  return;
}

void berechneQ2(double* temp, int age_block, double kappa1,int noa,int nop,double delta, double kappa2);

void machQ2(int swit, double* Q, double* Qnull, int age_block, int** daten_n, double* theta, double* phi, double* psi, int noa, int nop, int vdb, double my, double kappa2, double kappa)
{
  int b=(2*age_block)+1;
  berechneQ2(Q, age_block, kappa, noa, nop, 0.0, kappa2);
  //mxschreibe(Q, 2*noa,b);

  if (swit==1)
  {
    for (int i=0; i< noa; i++)
    {
      for (int t=0; t < nop; t++)
      {

        Q[((2*i)+1)*b] = Q[((2*i)+1)*b] + (daten_n[i][t]*taylor2(theta[i] + phi[t] + psi[coh(i,t,noa,vdb)-1] + my));
      }
    }
  }

  if (swit==2)
  {
    for (int i=0; i< noa; i++)
    {

      for (int t=0; t < nop; t++)
      {
        Q[((2*i)+1)*b] = Q[((2*i)+1)*b] + (daten_n[t][i]*taylor2(theta[i] + phi[t] + psi[coh(t,i,nop,vdb)-1] + my));
      }
    }
  }
  if (swit<0)
  {
    for (int i=0; i< (-swit); i++)
    {
      for (int t=0; t < nop; t++)
      {
        Q[(2*(coh(i,t,(-swit),vdb)-1)+1)*b] = Q[(2*(coh(i,t,(-swit),vdb)-1)+1)*b] + (daten_n[i][t]*taylor2(psi[i] + phi[t] + theta[coh(i,t,(-swit),vdb)-1] + my ));
      }
    }

  }

  return;
}
void MausQphi(double* Q, int b, int** daten_n, double* theta, double* phi, double* psi, int noa, int nop, int vdb, double my)
{

  for (int i=0; i< noa; i++)
  {
    for (int t=0; t < nop; t++)
    {
      Q[t*b] = Q[t*b] + ( daten_n[i][t]*taylor2(theta[i] + phi[t] + psi[coh(i,t,noa,vdb)-1] + my));
    }
  }

  return;
}
void MausQpsi(double* Q, int b, int** daten_n, double* theta, double* phi, double* psi, int noa, int nop, int vdb, double my)
{

  for (int i=0; i< noa; i++)
  {
    for (int t=0; t < nop; t++)
    {
      Q[(coh(i,t,noa,vdb)-1)*b] = Q[(coh(i,t,noa,vdb)-1)*b] + (daten_n[i][t]*taylor2(theta[i] + phi[t] + psi[coh(i,t,noa,vdb)-1] + my ));
    }
  }

  return;
}


void berechneBtaylorcohort(int* n, double* temp, double my, double* psi, double* phi, double* theta,int noa, int nop, int vdb, int noc, int** daten_n, int** daten_y)
{
  double offset;
  for (int i=0; i< noc; i++)
  {
    temp[i]=0.0;
    n[i]=0;
  }

  for (int i=0; i< noa; i++)
  {
    for (int j=0; j< nop; j++)
    {
      offset = my+phi[j]+theta[i]+psi[coh(i,j,noa,vdb)-1];
      temp[coh(i,j,noa,vdb)-1] += ((double)daten_y[i][j] - ((double)daten_n[i][j]*taylor1(offset)) + ((double)daten_n[i][j]*psi[coh(i,j,noa,vdb)-1]*taylor2(offset)));
      n[coh(i,j,noa,vdb)-1]++;
    }

  }

  return;
}


void berechneBtaylor(int swit, double* temp, double my, double* theta, double* phi, double* psi,int noa, int nop, int vdb, int** daten_n, int** daten_y)
{

  switch(swit)
  {
    double offset;
  case 1:

    for (int i=0; i< noa; i++)
    {
      temp[i]=0.0;
      for (int j=0; j< nop; j++)
      {
        offset=my+theta[i]+phi[j]+psi[coh(i,j,noa,vdb)-1];
        temp[i]+=(double)daten_y[i][j] - ((double)daten_n[i][j]*taylor1(offset)) + ((double)daten_n[i][j]*theta[i]*taylor2(offset));

      }
    }
    break;

  case 2:

    for (int i=0; i< noa; i++)
    {
      temp[i]=0.0;
      for (int j=0; j< nop; j++)
      {
        offset=my+theta[i]+phi[j]+psi[coh(j,i,nop,vdb)-1];
        temp[i]+=(double)daten_y[j][i] - ((double)daten_n[j][i]*taylor1(offset)) + ((double)daten_n[j][i]*theta[i]*taylor2(offset));
      }
    }
    break;

  }
  return;
}

double detQ(double*L, int n, int b)
{
  double det=0.0;

  for (int i=0; i<n; i++)
  {
    det += log(L[i*b]);
  }

  return det;
}
double loglikelihood(int swit, double my, double* theta, double* phi, double* psi, int** daten_y, int** daten_n,int age_block,int noa,int nop, int vdb, double kappa){

  double erg=0.0;

  if (swit==1)
  {
    for (int i=0; i<noa; i++)
    {
      for (int j=0; j<nop; j++)
      {
        erg += theta[i]*daten_y[i][j];
        erg -= daten_n[i][j] * log( 1+ exp(my+ theta[i]+phi[j]+psi[coh(i,j,noa,vdb)-1]));
      }
    }
  }
  if (swit==2)
  {
    for (int i=0; i< nop; i++)
    {
      for (int j=0; j<noa; j++)
      {
        erg += theta[j]*daten_y[i][j];
        erg -= daten_n[i][j]*log( 1+exp(my+  phi[i]+theta[j]+psi[coh(i,j,nop,vdb)-1] ) );
      }
    }
  }
  if (swit<0)
  {
    for (int i=0; i< (-swit); i++)
    {
      for (int j=0; j<nop; j++)
      {
        erg += theta[coh(i,j,-swit,vdb)-1]*daten_y[i][j];
        erg -= daten_n[i][j]*log( 1+exp (my+  psi[i]*phi[j]*theta[coh(i,j,-swit,vdb)-1]));
      }
    }
  }

  if (age_block==1)
  {
    for (int i=1; i< noa; i++)
    {
      erg -= (theta[i]-theta[i-1])*(theta[i]-theta[i-1])*kappa/2.0;
    }
  }

  if (age_block==2)
  {
    for (int i=2; i< noa; i++)
    {
      erg -= (theta[i]-2*theta[i-1]+theta[i-2])*(theta[i]-2*theta[i-1]+theta[i-2])*kappa/2.0;
    }
  }
  return erg;
}
double loglikelihood2o(int swit, double my, double* thetastern, double* theta2, double* phi, double* psi, int** daten_y, int** daten_n,int age_block,int noa,int nop, int vdb, double kappa, double kappa2){

  double erg=0.0;

  if (swit==1)
  {
    for (int i=0; i<noa; i++)
    {
      for (int j=0; j<nop; j++)
      {
        erg += thetastern[i]*daten_y[i][j];
        erg -= daten_n[i][j] * log( 1+ exp(my+ thetastern[i]+phi[j]+psi[coh(i,j,noa,vdb)-1]));
      }
    }
  }
  if (swit==2)
  {
    for (int i=0; i< nop; i++)
    {
      for (int j=0; j<noa; j++)
      {
        erg += thetastern[j]*daten_y[i][j];
        erg -= daten_n[i][j]*log( 1+exp(my+  phi[i]+thetastern[j]+psi[coh(i,j,nop,vdb)-1] ) );
      }
    }
  }
  if (swit<0)
  {
    for (int i=0; i< (-swit); i++)
    {
      for (int j=0; j<nop; j++)
      {
        erg += thetastern[coh(i,j,-swit,vdb)-1]*daten_y[i][j];
        erg -= daten_n[i][j]*log( 1+exp (my+  psi[i]*phi[j]*thetastern[coh(i,j,-swit,vdb)-1]));
      }
    }
  }

  for (int i=0; i< noa; i++)
  {
    erg -= theta2[i]*theta2[i]*kappa2/2.0;
  }

  if (age_block==1)
  {
    for (int i=1; i< noa; i++)
    {
      erg -= (thetastern[i]-theta2[i]-thetastern[i-1]+theta2[i-1])*(thetastern[i]-theta2[i]-thetastern[i-1]+theta2[i-1])*kappa/2.0;
    }
  }

  if (age_block==2)
  {
    for (int i=2; i< noa; i++)
    {
      erg -= (thetastern[i]-theta2[i]-2.0*thetastern[i-1]+2.0*theta2[i-1]+thetastern[i-2]-theta2[i-2])*(thetastern[i]-theta2[i]-2.0*thetastern[i-1]+2.0*theta2[i-1]+thetastern[i-2]-theta2[i-2])*kappa/2.0;
    }
  }
  return erg;
}
double loglikelihood2(int swit, double my, double* theta, double* phi, double* psi, int** daten_y, int** daten_n,int age_block,int noa,int nop, int vdb, double kappa, double kappa2){

  double erg=0.0;

  if (swit==1)
  {
    for (int i=0; i<noa; i++)
    {
      for (int j=0; j<nop; j++)
      {
        erg += theta[(2*i)+1]*daten_y[i][j];
        erg -= daten_n[i][j] * log( 1+ exp(my+ theta[(2*i)+1]+phi[j]+psi[coh(i,j,noa,vdb)-1]));
      }
    }
  }
  if (swit==2)
  {
    for (int i=0; i< nop; i++)
    {
      for (int j=0; j<noa; j++)
      {
        erg += theta[(2*j)+1]*daten_y[i][j];
        erg -= daten_n[i][j]*log( 1+exp(my+  phi[i]+theta[(2*j)+1]+psi[coh(i,j,nop,vdb)-1] ) );
      }
    }
  }
  if (swit<0)
  {
    for (int i=0; i< (-swit); i++)
    {
      for (int j=0; j<nop; j++)
      {
        erg += theta[((coh(i,j,-swit,vdb)-1)*2)+1]*daten_y[i][j];
        erg -= daten_n[i][j]*log( 1+exp (my+  psi[i]*phi[j]*theta[(2*(coh(i,j,-swit,vdb)-1))+1]));
      }
    }
  }

  for (int i=0; i< noa; i++)
  {
    erg -= ((theta[(2*i)+1]-theta[2*i])*(theta[(2*i)+1]-theta[2*i])*kappa2/2.0);
  }


  if (age_block==1)
  {
    for (int i=1; i< noa; i++)
    {
      erg -= (theta[2*i]-theta[2*(i-1)])*(theta[2*i]-theta[2*(i-1)])*kappa/2.0;
    }
  }

  if (age_block==2)
  {
    for (int i=2; i< noa; i++)
    {
      erg -= (theta[2*i]-2*theta[2*(i-1)]+theta[2*(i-2)])*(theta[2*i]-2*theta[2*(i-1)]+theta[2*(i-2)])*kappa/2.0;
    }
  }
  return erg;
}





double xstrichy(double* x,double* y,int n)
{
  double erg=0.0;
  for (int i=0; i<n; i++)
  {
    erg += x[i]*y[i];
  }
  return erg;
}


void bedinge_lik(int age_block, int noa, double* theta, double* L, double* mu, int k, double* A, double* b, double& lik)
{

  /*
  theta: unbedingtes sample
  thetatemp: ?
  A
  */

  lik=0.0;

  int bandw=age_block;
  double* V;
  V=new double[k*noa];
  double* temp1 = new double[k];
  double* temp2 = new double[k];
  double* temp3 = new double[1];
  double* thetatemp = new double[noa];


  for (int i=0; i<k; i++)
  {
    // erzeuge A^T_i
    double* Atransi=new double[noa];
    for (int j=0; j<noa; j++)
    {
      Atransi[j]=A[(i*noa)+j];
      // 	cout <<  Atransi[j]<< " ";
    }
    // Lu=A^T_i
    loese2(L, Atransi, noa, bandw);
    //mxschreibe(Atransi,1,noa);
    //mxschreibe(L,bandw,noa);    // L^T*V_i=u
    loese(L, Atransi, noa, bandw);
    //mxschreibe(Atransi,1,noa);
    //erzeuge V
    //int jj=0;
    // for (int
    for (int j=0; j<noa; j++)
    {

      V[j*k+i]=Atransi[j];
      //  	cout <<  Atransi[j]<< " ";
    }
    delete[] Atransi;
    //      cout << endl;
  }
  //   cout << "--" << endl;
  //  mxschreibe(V,noa,k);
  double* y=new double[k];
  //Ax
  multipliziere(A,theta,k,noa,1,y);
  //Ax-b
  for (int i=0; i<k; i++)
  {
    y[i]=y[i]-b[i];
  }

  double* AV=new double[k*k];

  //A*V
  multipliziere(A,V,k,noa,k,AV);
  //  cout << "1 "<<lik<< " " <<AV[0] << endl;
  lik -= det(AV,k)/2.0;
  //  cout << "2 "<< lik << endl;
  //(AV)^-1
  invers(AV,k);

  //V^T*?
  multipliziere(thetatemp,V,1,noa,k,temp1);

  //(V^T*?)^T*(AV)^-1*(V^T?)
  multipliziere(temp1,AV,1,k,k,temp2);
  multipliziere(temp2,temp1,1,k,1,temp3);

  lik -= temp3[0]/2.0;
  //  cout << lik << endl; cout << endl;
  double* W=new double[noa*k];
  //W=V(AV)^-1
  multipliziere(V,AV,noa,k,k,W);

  //W(Ax-b)
  multipliziere(W,y,noa,k,1,thetatemp);

  //berechne x-W(Ax-b)
  for (int i=0; i< noa; i++)
  {
    theta[i]=theta[i]-thetatemp[i];
  }

  delete[] V;
  delete[] y;
  delete[] W;
  delete[] AV;
  delete[] temp1;
  delete[] temp2;
  delete[] temp3;
  delete[] thetatemp;

  return;

}
double bedinge_lik2(int age_block, int noa, double* theta, double* Q, double* mu, int k, double* A, double* b)
{

  /*
  theta: unbedingtes sample
  thetatemp: ?
  A
  */

  double lik=0.0;

  int bandw=age_block;
  double* AV=new double[k*k];
  double* V=new double[noa*k];
  double* M=new double[noa*noa];
  double* temp1 = new double[k];
  double* temp2 = new double[k];
  double* temp3 = new double[1];
  double* thetatemp = new double[noa];
  double* W=new double[noa*k];
  double* y=new double[k];

  for (int i=0; i<noa; i++)
  {
    for (int j=0; j<noa; j++)
    {
      if ((ABS(i-j))<bandw)
      {
        M[(i*noa)+j]=Q[((int)min(i,j)*bandw)+(int)ABS(i-j)];
        // 	     cout << i<<":"<<j<<"="<<(i*noa)+j <<" / " << ((int)min(i,j)*bandw)+(int)ABS(i-j)<<endl;
      }
      else
      {
        M[(i*noa)+j]=0.0;
      }
      //  	 cout << M[(i*noa)+j]<<"  | ";
    }
    //      cout << endl;
  }
  //  cout << endl;


  invers(M,noa);



  multipliziere(M,A,noa,noa,k,V);

  //Ax
  multipliziere(A,theta,k,noa,1,y);
  //Ax-b
  for (int i=0; i<k; i++)
  {
    y[i]=y[i]-b[i];
  }



  //A*V
  multipliziere(A,V,k,noa,k,AV);
  //cout << "1 "<<lik<< " " <<AV[0] << endl;
  lik -= det(AV,k)/2.0;
  //cout << "2 "<< lik << endl;
  //(AV)^-1
  invers(AV,k);

  //V^T*?
  multipliziere(thetatemp,V,1,noa,k,temp1);

  //(V^T*?)^T*(AV)^-1*(V^T?)
  multipliziere(temp1,AV,1,k,k,temp2);
  multipliziere(temp2,temp1,1,k,1,temp3);

  lik -= temp3[0]/2.0;
  //cout << lik << endl; cout << endl;

  //W=V(AV)^-1
  multipliziere(V,AV,noa,k,k,W);

  //W(Ax-b)
  multipliziere(W,y,noa,k,1,thetatemp);

  //berechne x-W(Ax-b)
  for (int i=0; i< noa; i++)
  {
    theta[i]=theta[i]-thetatemp[i];
  }


  delete[] y;
  delete[] W;
  //  delete[] V;
  delete[] AV;
  delete[] M;
  delete[] temp1;
  delete[] temp2;
  delete[] temp3;
  delete[] thetatemp;

  return lik;

}

double lik_bedingt(int age_block, int noa, double* theta, double* Q, double* mu, int k, double* A, double* b)
{

  /*
  theta: unbedingtes sample
  thetatemp: ?
  A
  */

  double lik=0.0;

  int bandw=age_block;
  double* V;
  V=new double[k*noa];
  double* temp1 = new double[k];
  double* temp2 = new double[k];
  double* temp3 = new double[1];
  double* thetatemp = new double[noa];


  double* M;
  M=new double[noa*noa];

  for (int i=0; i<noa; i++)
  {
    for (int j=0; j<noa; j++)
    {
      if ((ABS(i-j))<bandw)
      {
        M[(i*noa)+j]=Q[((int)min(i,j)*bandw)+(int)ABS(i-j)];
      }
      else
      {
        M[(i*noa)+j]=0.0;
      }
      // 	 cout << M[(i*noa)+j]<<"  | ";
    }
    //      cout << endl;
  }
  //  cout << endl;

  invers(M,noa);



  multipliziere(M,A,noa,noa,k,V);

  //   cout << "--" << endl;
  // mxschreibe(V,noa,k);
  double* y=new double[k];
  //Ax
  multipliziere(A,theta,k,noa,1,y);
  //Ax-b
  for (int i=0; i<k; i++)
  {
    y[i]=y[i]-b[i];
  }

  double* AV=new double[k*k];

  //A*V
  multipliziere(A,V,k,noa,k,AV);

  lik -= det(AV,k)/2.0;

  //(AV)^-1
  invers(AV,k);

  //V^T*?
  multipliziere(thetatemp,V,1,noa,k,temp1);

  //(V^T*?)^T*(AV)^-1*(V^T?)
  multipliziere(temp1,AV,1,k,k,temp2);
  multipliziere(temp2,temp1,1,k,1,temp3);

  lik -= temp3[0]/2.0;



  delete[] V;
  delete[] y;
  delete[] M;
  delete[] AV;
  delete[] temp1;
  delete[] temp2;
  delete[] temp3;
  delete[] thetatemp;

  return lik;

}




void blocken(int swit, int age_block, double kappa, double delta, int noa, int nop, double &my, double* theta, double* phi, double* psi, double* ageQ, double* thetatemp, int vdb, int** daten_n, int** daten_y, int &ja)
{


  int bandw=age_block;
  int* nwerte;
  if (swit<0)
  {
    nwerte=new int[noa];
    berechneBtaylorcohort(nwerte, thetatemp, my, theta, phi, psi, -swit, nop,vdb, noa, daten_n, daten_y);
    berechneQcohort(nwerte, ageQ, age_block, kappa, noa, 0.0);
    delete[] nwerte;
    MausQpsi(ageQ, age_block+1, daten_n, psi, phi, theta, -swit, nop, vdb, my);
  }
  else
  {
    berechneBtaylor(swit, thetatemp, my, theta, phi, psi, noa, nop,vdb, daten_n, daten_y);
    berechneQ(ageQ, age_block, kappa, noa, nop, 0.0);
  }

  if (swit==1){MausQtheta(ageQ, age_block+1, daten_n, theta, phi, psi, noa, nop, vdb, my);}
  if (swit==2){MausQphi(ageQ, age_block+1, daten_n, phi, theta, psi, nop, noa, vdb, my);}

  //  cout << "B:"<< endl;
  //  mxschreibe(thetatemp,1,noa);

  //thetatemp == b
  //ageQ == M

  double* thetaproposal=new double[noa];

  ageQ = cholesky(noa, ageQ, bandw);

  //ageQ == L

  double* L=new double[(age_block+1)*noa];
  for (int i=0; i< (age_block+1)*noa; i++)
  {
    L[i]=ageQ[i];
  }

  //Lv=b
  loese2(L, thetatemp, noa, bandw);

  //L^T my=v
  loese(L,thetatemp, noa, bandw);

  //  cout << "?"<<endl;
  //  mxschreibe(thetatemp,1,noa);

  double* mu= new double[noa];
  for (int i=0; i<noa; i++)
  {
    mu[i]=thetatemp[i];
  }

  // mu==M^-1m

  //z
  gausssample(thetaproposal,noa);
  //L^T y=z
  loese(L,thetaproposal,noa,bandw);


  //x=my + y

  double sum=0.0;
  for (int i=0; i< noa; i++)
  {
    thetaproposal[i]=thetaproposal[i]+thetatemp[i];
    sum += thetaproposal[i];
  }

  // und zentrieren

  for (int i=0; i< noa; i++)
  {
    thetaproposal[i] = thetaproposal[i] - sum/double(noa);
  }

  double logalpha=0.0;

  logalpha += loglikelihood(swit, my,  thetaproposal, phi, psi, daten_y, daten_n, age_block, noa, nop, vdb, kappa);
  logalpha -= loglikelihood(swit, my, theta, phi, psi, daten_y, daten_n, age_block, noa, nop, vdb, kappa);
  logalpha += detQ(ageQ, noa, age_block+1);


  for (int i=0; i<noa; i++)
  {
    thetatemp[i]=thetaproposal[i]-mu[i];
  }
  logalpha -= 0.5*xLx(ageQ,thetatemp,noa,age_block+1);


  //Berechne m: Erwartungswert von \sum x_i
  //Berechne s: Varianz von \sum x_i

  double S_new=berechneVarianzsumx_i(L,noa,bandw);
  double m_new=0.0;

  for (int i=0; i< noa; i++)
  {
    m_new += mu[i];
  }

  logalpha -= 0.5*( log(S_new) + m_new*m_new/S_new );

  // Berechne M und m f?r altes theta

  if (swit<0)
  {
    nwerte=new int[noa];
    berechneBtaylorcohort(nwerte, thetatemp, my, thetaproposal, phi, psi, -swit, nop,vdb, noa, daten_n, daten_y);
    berechneQcohort(nwerte, ageQ, age_block, kappa, noa, 0.0);
    delete[] nwerte;
    MausQpsi(ageQ, age_block+1, daten_n, psi, phi, thetaproposal, -swit, nop, vdb, my);
  }
  else
  {
    berechneBtaylor(swit, thetatemp, my, thetaproposal, phi, psi, noa, nop,vdb, daten_n, daten_y);
    berechneQ(ageQ, age_block, kappa, noa, nop, 0.0);
  }

  if (swit==1){MausQtheta(ageQ, age_block+1, daten_n, thetaproposal, phi, psi, noa, nop, vdb, my);}
  if (swit==2){MausQphi(ageQ, age_block+1, daten_n, phi, thetaproposal, psi, nop, noa, vdb, my);}


  //Lv=b
  loese2(ageQ, thetatemp, noa, bandw);
    //L^T my=v
  loese(ageQ,thetatemp, noa, bandw);

  for (int i=0; i< noa; i++)
  {
    mu[i]=thetatemp[i];
    thetatemp[i]=theta[i]-mu[i];
  }


    // mu=M^-1m, thetatemp = theta - M^-1m, ageQ = L


  //  logalpha -= detQ(ageQ, noa, age_block+1);


  //Berechne m: Erwartungswert von \sum x_i
  //Berechne s: Varianz von \sum x_i

  double S_old=berechneVarianzsumx_i(ageQ,noa,bandw);
  double m_old=0.0;
  for (int i=0; i< noa; i++)
  {
    m_old += mu[i];
  }

  logalpha += 0.5*xLx(ageQ,thetatemp,noa,age_block+1);
  logalpha += 0.5*(log(S_old) + m_old*m_old/S_old);
  logalpha -= detQ(ageQ, noa, age_block+1);

  logalpha = exp(logalpha);

  if (logalpha > nulleins())
  {
    for (int i=0; i<noa; i++)
    {
      theta[i]=thetaproposal[i];
    }
    ja ++;
    //  cout << "OK: ";
  }

  //  cout << "Vorschlag: "<< lvdo-lvdp << "  Likelihood: "<<llhp-llho<<endl;


  //  else
  //    {
  //      blocken(swit, age_block, kappa, delta, noa, nop, my,  theta,  phi,  psi,  ageQ,  thetatemp,  vdb, daten_n, daten_y, ja);
  //    }

  //  mxschreibe(thetaproposal,1,noa);
  //  mxschreibe(theta,1,noa);
  //  if (swit<0)
  //    {
  //      cout << lvdo << " " <<lvdp << endl;
  //      ofstream dev2( (string("/home/schmidt/temp/proposal.dat")).c_str(),ios::app);
  //      for (int i=0; i< noa; i++)
  //        {
  // 	 dev2 << thetaproposal[i]<< " ";
  //        }
  //      dev2 << endl;
  //    }



  delete[] L;


  delete[] thetaproposal;
  delete[] mu;


  return;

}



void blocken2(int swit, int age_block, double kappa, double kappa2, int noa, int nop, double &my, double* theta, double* theta2, double* phi, double* psi, double* ageQ, double* thetatemp, int vdb, int** daten_n, int** daten_y, int &ja)
{

  double* thetatemp2 = new double[2*noa];
  double* ageQ2 = new double[2*noa*((2*age_block)+1)];

  age_block=age_block-2;
  int bandw=(2*age_block)+1;
  int bw=bandw-1;
  int* nwerte;
  int noa_2=2*noa;

  if (swit<0)
  {
    nwerte=new int[noa];
    berechneBtaylorcohort(nwerte, thetatemp, my, theta, phi, psi, -swit, nop,vdb, noa, daten_n, daten_y);
    delete[] nwerte;
  }
  else
  {
    berechneBtaylor(swit, thetatemp, my, theta, phi, psi, noa, nop,vdb, daten_n, daten_y);

  }


  machQ2(swit, ageQ2, ageQ, age_block, daten_n, theta, phi, psi, noa, nop, vdb, my, kappa2, kappa);


  //mxschreibe(ageQ2,2*noa,bandw);

  double* b= new double[2*noa];

  for (int i=0; i<noa; i++)
  {
    thetatemp2[2*i+1]=thetatemp[i];
  }
  for (int i=0; i<noa; i++)
  {
    thetatemp2[2*i]=0.0;
  }

  for (int i=0; i<(2*noa); i++)
  {
    b[i]=thetatemp2[i];
  }
  // cout << "B:"<< endl;
  //   mxschreibe(thetatemp2,1,2*noa);
  double* Q=new double[((2*age_block)+1)*2*noa];
  for (int i=0; i< ((2*age_block)+1)*2*noa; i++)
  {
    Q[i]=ageQ2[i];
  }
  // cout<< "Q:";
  // mxschreibe(Q,2*noa,bandw);
  double* thetaproposal=new double[2*noa];
  ageQ2 = cholesky(2*noa, ageQ2, bw);



  // mxschreibe(ageQ2,2*noa,bandw);
  double* L=new double[bandw*2*noa];
  for (int i=0; i< bandw*2*noa; i++)
  {
    L[i]=ageQ2[i];
  }
  // mxschreibe(L,2*noa,bandw);

  //Lv=b
  loese2(L, thetatemp2, noa_2, bw);
  // mxschreibe(thetatemp2,1,2*noa);
  //L^T my=v
  loese(L,thetatemp2, noa_2, bw);
  // mxschreibe(thetatemp2,1,2*noa);
  //  cout << "?"<<endl;
  //  mxschreibe(thetatemp,1,noa);

  double* mu= new double[noa*2];
  for (int i=0; i<2*noa; i++)
  {
    mu[i]=thetatemp2[i];
  }




  //z
  gausssample(thetaproposal,2*noa);
  //mxschreibe(thetaproposal,1,2*noa);
  double* z= new double[2*noa];
  for (int i=0; i<2*noa; i++)
  {
    z[i]=thetaproposal[i];
  }

  //L^T y=z
  loese(L,thetaproposal,noa_2,bw);
  // mxschreibe(thetaproposal,1,2*noa);
  //x=my + y




  //Berechne m: Erwartungswert von \sum x_i
  //Berechne s: Varianz von \sum x_i

  //  double S_new=berechneVarianzsumx_i(L,noa,bandw);

  //  double m_new=0.0;

  for (int i=0; i< 2*noa; i++)
  {
    //     m_new += thetatemp[i];
    thetaproposal[i]=thetaproposal[i]+thetatemp2[i];

  }
  //mxschreibe(thetaproposal,1,noa);


  // // nur \sum theta1 = 0

  int k=1;
  double* b_A=new double[k*noa_2];
  double* b_b=new double[k];


  b[0]=0.0;
  for (int i=0; i< noa; i++)
  {
    b_A[2*i]=
      b_A[(2*i)+1]=1.0;
  }
  // nur \sum theta1 = 0

  //  int k=2;
  //  double* b_A=new double[k*noa_2];
  //  double* b_b=new double[k];


  //  b[0]=0.0; b[1]=0.0;
  //  for (int i=0; i< noa; i++)
  //    {
  //      b_A[2*i]=1.0;
  //      b_A[(2*i)+1]=0.0;
  //      b_A[2*(i+noa)]=-1.0;
  //      b_A[(2*(i+noa))+1]=1.0;
  //    }







  double logalpha=0.0;


  //  mxschreibe(thetaproposal, 1, noa_2);

  //  logalpha += bedinge_lik2(bandw, noa_2, thetaproposal, Q, b, k, b_A, b_b);
  //  mxschreibe(thetaproposal, 1, noa_2);
  //  cout << endl;
  logalpha += loglikelihood2(swit, my,  thetaproposal, phi, psi, daten_y, daten_n, age_block, noa, nop, vdb, kappa, kappa2);

  logalpha -= loglikelihood2o(swit, my, theta, theta2, phi, psi, daten_y, daten_n, age_block, noa, nop, vdb, kappa, kappa2);

  //  cout << logalpha << endl;
  //Proposaldichte
  logalpha += detQ(ageQ2, 2*noa, bandw);
  //  cout <<":" <<logalpha << endl;
  for (int i=0; i<noa_2; i++)
  {

    thetatemp2[i]=thetaproposal[i]-mu[i];
    // thetaproposal[i]=thetaproposal[i]-sum;
  }
  logalpha -= xLx(ageQ2,thetatemp2,2*noa,bandw)/2.0;
  // cout << logalpha << endl;


  double* thetaproposal1=new double[noa];

  for (int i=0; i<noa; i++)
  {
    thetaproposal1[i]=thetaproposal[(2*i)+1];
  }

  if (swit<0)
  {
    nwerte=new int[noa];
    berechneBtaylorcohort(nwerte, thetatemp, my, thetaproposal1, phi, psi, -swit, nop,vdb, noa, daten_n, daten_y);

    delete[] nwerte;
  }
  else
  {
    berechneBtaylor(swit, thetatemp, my, thetaproposal1, phi, psi, noa, nop,vdb, daten_n, daten_y);

  }
  machQ2(swit, ageQ2, ageQ, age_block, daten_n, thetaproposal1, phi, psi, noa, nop, vdb, my, kappa2, kappa);



  for (int i=0; i<noa; i++)
  {
    thetatemp2[(2*i)+1]=thetatemp[i];
  }
  for (int i=0; i<noa; i++)
  {
    thetatemp2[(2*i)]=0.0;
  }
  for (int i=0; i<2*noa; i++)
  {
    b[i]=thetatemp2[i];
  }
  for (int i=0; i< bandw*2*noa; i++)
  {
    L[i]=ageQ2[i];
  }
  L = cholesky(2*noa, L, bw);

  //Lv=b
  loese2(L, thetatemp2, noa_2, bw);
  //L^T my=v
  loese(L,thetatemp2, noa_2, bw);

  // cout << logalpha << endl;

  logalpha -= detQ(L, 2*noa, bandw);

  for (int i=0; i< noa; i++)
  {
    thetatemp2[2*i] = theta[i]-thetatemp2[2*i];
    thetatemp2[(2*i)+1] = theta[i]+theta2[i]-thetatemp2[2*i];
  }

  logalpha += xLx(L, thetatemp2, 2*noa, bandw);

  // cout <<"?"<< logalpha << endl;

  for (int i=0; i< noa; i++)
  {
    thetatemp2[2*i]=theta[i];
    thetatemp2[(2*i)+1]=theta[i]+theta2[i];
  }

  //  logalpha -= lik_bedingt(bandw, noa_2, thetatemp2, ageQ2, b, k, b_A, b_b);


  // cout << logalpha << endl;
  logalpha = exp(logalpha);

  if (logalpha > nulleins())
  {
    for (int i=0; i<noa; i++)
    {
      theta[i]=thetaproposal[(2*i)+1];
      theta2[i]=theta[i]-thetaproposal[2*i];
    }
    ja ++;
    //  cout << "OK: ";
  }

  //  cout << "Vorschlag: "<< lvdo-lvdp << "  Likelihood: "<<llhp-llho<<endl;


  //  else
  //    {
  //      blocken(swit, age_block, kappa, delta, noa, nop, my,  theta,  phi,  psi,  ageQ,  thetatemp,  vdb, daten_n, daten_y, ja);
  //    }

  //  mxschreibe(theta,1,noa);
  //  mxschreibe(theta2,1,noa);
  //  cout << endl;

  double summe=0;
  double summe2=0;
  for (int i=0; i<noa; i++)
  {
    summe+=theta[i];
    summe2+=theta2[i];
  }
  for (int i=0; i<noa; i++)
  {
    theta[i]=theta[i]-(summe/double(noa));
    theta2[i]=theta2[i]-(summe2/double(noa));
  }



  delete[] L;
  delete[] b;
  delete[] Q;
  delete[] thetaproposal;
  delete[] thetaproposal1
    ; delete[] mu;
    delete[] z;
    delete[] thetatemp2;
    delete[] ageQ2;
    delete[] b_A;
    delete[] b_b;

    return;

}

#include <iostream>
using namespace std;
#include "matrix.h"
#include "l.h"
#include "mxs.h"
#include "zufall.h"

void berechneQ(double* temp, int age_block, double kappa,int noa,int nop,double delta);
void berechneQcohort(int* n, double* temp, int age_block, double kappa,int noa,double delta);

//Berechnet Matrix Q (f?r age, period)
void berechneQ2(double* temp, int age_block, double kappa1,int noa,int nop,double delta, double kappa2)
{

  if (age_block==1)
  {
    int index=0;
    temp[index]=kappa1+kappa2;
    index++;
    temp[index]=-kappa2;
    index++;
    temp[index]=-kappa1;
    index++;
    temp[index]=kappa2+delta*nop;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    for (int i=1; i<noa-1; i++)
    {
      temp[index]=2*kappa1+kappa2;
      index++;
      temp[index]=-kappa2;
      index++;
      temp[index]=-kappa1;
      index++;
      temp[index]=kappa2+delta*nop;
      index++;
      temp[index]=0.0;
      index++;
      temp[index]=0.0;
      index++;
    }
    temp[index]=kappa1+kappa2;
    index++;
    temp[index]=-kappa2;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=kappa2+delta*nop;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
  }

  // // TO:DO
  if (age_block==2)
  {
    //1.Zeile
    int index=0;
    temp[index]=kappa1+kappa2;
    index++;
    temp[index]=-kappa2;
    index++;
    temp[index]=-2*kappa1;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=kappa1;
    index++;
    //2.Zeile
    temp[index]=kappa2+delta*nop;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    //3.Zeile
    temp[index]=5*kappa1+kappa2;
    index++;
    temp[index]=-kappa2;
    index++;
    temp[index]=-4*kappa1;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=kappa1;
    index++;
    //4.Zeile
    temp[index]=kappa2+delta*nop;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;

    for (int i=2; i<noa-2; i++)
    {
      temp[index]=6*kappa1+kappa2;
      index++;
      temp[index]=-kappa2;
      index++;
      temp[index]=-4*kappa1;
      index++;
      temp[index]=0.0;
      index++;
      temp[index]=kappa1;
      index++;

      temp[index]=kappa2+delta*nop;
      index++;
      temp[index]=0.0;
      index++;
      temp[index]=0.0;
      index++;
      temp[index]=0.0;
      index++;
      temp[index]=0.0;
      index++;


    }
    //-4.Zeile

    temp[index]=5*kappa1+kappa2;
    index++;
    temp[index]=-kappa2;
    index++;
    temp[index]=-2*kappa1;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    //-3.Zeile
    temp[index]=kappa2+delta*nop;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    //-2.Zeile

    temp[index]=kappa1+kappa2;
    index++;
    temp[index]=-kappa2;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    //-1.Zeile
    temp[index]=kappa2+delta*nop;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;

  }
  return;
}

void berechneQcohort2(int* n, double* temp, int age_block, double kappa1,int noa,double delta, double kappa2)
{

  if (age_block==1)
  {
    int index=0;
    temp[index]=kappa1+kappa2;
    index++;
    temp[index]=-kappa2;
    index++;
    temp[index]=-kappa1;
    index++;
    temp[index]=kappa2+delta*n[0];
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    for (int i=1; i<noa-1; i++)
    {
      temp[index]=2*kappa1+kappa2;
      index++;
      temp[index]=-kappa2;
      index++;
      temp[index]=-kappa1;
      index++;
      temp[index]=kappa2+delta*n[i];
      index++;
      temp[index]=0.0;
      index++;
      temp[index]=0.0;
      index++;
    }
    temp[index]=kappa1+kappa2;
    index++;
    temp[index]=-kappa2;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=kappa2+delta*n[noa-1];
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
  }
  if (age_block==2)
  {
    int index=0;
    temp[index]=kappa1+kappa2;
    index++;
    temp[index]=-kappa2;
    index++;
    temp[index]=-2*kappa1;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=kappa1;
    index++;
    //2.Zeile
    temp[index]=kappa2+delta*n[0];
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    //3.Zeile
    temp[index]=5*kappa1+kappa2;
    index++;
    temp[index]=-kappa2;
    index++;
    temp[index]=-4*kappa1;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=kappa1;
    index++;
    //4.Zeile
    temp[index]=kappa2+delta*n[1];
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;

    for (int i=2; i<noa-2; i++)
    {
      temp[index]=6*kappa1+kappa2;
      index++;
      temp[index]=-kappa2;
      index++;
      temp[index]=-4*kappa1;
      index++;
      temp[index]=0.0;
      index++;
      temp[index]=kappa1;
      index++;

      temp[index]=kappa2+delta*n[i];
      index++;
      temp[index]=0.0;
      index++;
      temp[index]=0.0;
      index++;
      temp[index]=0.0;
      index++;
      temp[index]=0.0;
      index++;


    }
    //-4.Zeile

    temp[index]=5*kappa1+kappa2;
    index++;
    temp[index]=-kappa2;
    index++;
    temp[index]=-4*kappa1;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=kappa1;
    index++;
    //-3.Zeile
    temp[index]=kappa2+delta*n[noa-2];
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    //-2.Zeile

    temp[index]=kappa1+kappa2;
    index++;
    temp[index]=-kappa2;
    index++;
    temp[index]=-2*kappa1;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=kappa1;
    index++;
    //-1.Zeile
    temp[index]=kappa2+delta*n[noa-1];
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
  }
  return;
}
//Berechnet Matrix Q (f?r age, period)
void berechneQ3(double* temp, int age_block, double kappa1,int noa,int nop,double delta, double kappa2)
{

  if (age_block==1)
  {
    int index=0;
    temp[index]=kappa1+(delta*nop);
    index++;
    temp[index]=-(delta*nop);
    index++;
    temp[index]=-kappa1;
    index++;
    temp[index]=kappa2+(delta*nop);
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    for (int i=1; i<noa-1; i++)
    {
      temp[index]=2*kappa1+(delta*nop);
      index++;
      temp[index]=-(delta*nop);
      index++;
      temp[index]=-kappa1;
      index++;
      temp[index]=kappa2+(delta*nop);
      index++;
      temp[index]=0.0;
      index++;
      temp[index]=0.0;
      index++;
    }
    temp[index]=kappa1+(delta*nop);
    index++;
    temp[index]=-(delta*nop);
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=kappa2+(delta*nop);
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
  }

  // // TO:DO
  // if (age_block==2)
  //   {
  //     int index=0;
  //     temp[index]=kappa+delta*nop;
  //     index++;
  //     temp[index]=-2*kappa;
  //     index++;
  //     temp[index]=kappa;
  //     index++;
  //     temp[index]=5*kappa+delta*nop;
  //     index++;
  //     temp[index]=-4*kappa;
  //     index++;
  //     temp[index]=kappa;
  //     index++;
  //     for (int i=2; i<noa-2; i++)
  //       {
  // 	temp[index]=6*kappa+delta*nop;
  // 	index++;
  // 	temp[index]=-4*kappa;
  // 	index++;
  // 	temp[index]=kappa;
  // 	index++;
  //       }
  //     temp[index]=5*kappa+delta*nop;
  //     index++;
  //     temp[index]=-2*kappa;
  //     index++;
  //     //temp[index]=0;
  //     index++;
  //     temp[index]=kappa+delta*nop;
  //   }
  return;
}

void berechneQcohort3(int* n, double* temp, int age_block, double kappa1,int noa,double delta, double kappa2)
{

  if (age_block==1)
  {
    int index=0;
    temp[index]=kappa1+(delta*n[0]);
    index++;
    temp[index]=(delta*n[0]);
    index++;
    temp[index]=-kappa1;
    index++;
    temp[index]=kappa2+(delta*n[0]);
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
    for (int i=1; i<noa-1; i++)
    {
      temp[index]=2*kappa1+(delta*n[i]);
      index++;
      temp[index]=(delta*n[i]);
      index++;
      temp[index]=-kappa1;
      index++;
      temp[index]=kappa2+(delta*n[i]);
      index++;
      temp[index]=0.0;
      index++;
      temp[index]=0.0;
      index++;
    }
    temp[index]=kappa1+(delta*n[noa-1]);
    index++;
    temp[index]=(delta*n[noa-1]);
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=kappa2+(delta*n[noa-1]);
    index++;
    temp[index]=0.0;
    index++;
    temp[index]=0.0;
    index++;
  }
  // if (age_block==2)
  //   {
  //     int index=0;
  //     temp[index]=kappa+delta*n[0];
  //     index++;
  //     temp[index]=-2*kappa;
  //     index++;
  //     temp[index]=kappa;
  //     index++;
  //     temp[index]=5*kappa+delta*n[1];
  //     index++;
  //     temp[index]=-4*kappa;
  //     index++;
  //     temp[index]=kappa;
  //     index++;
  //     for (int i=2; i<noa-2; i++)
  //       {
  // 	temp[index]=6*kappa+delta*n[i];
  // 	index++;
  // 	temp[index]=-4*kappa;
  // 	index++;
  // 	temp[index]=kappa;
  // 	index++;
  //       }
  //     temp[index]=5*kappa+delta*n[noa-2];
  //     index++;
  //     temp[index]=-2*kappa;
  //     index++;
  //    // temp[index]=0;
  //     index++;
  //     temp[index]=kappa+delta*n[noa-1];
  //     index++;
  //    // temp[index]=0;
  //     index++;
  //     //temp[index]=0;
  //   }
  return;
}


void berechneB(int swit, double* temp, double** ksi, double my,  double* phi, double* psi,int noa, int nop, int vdb, double delta);

void berechneBcohort(int* n, double* temp, double** ksi, double my, double* phi, double* theta,int noa, int nop, int vdb, double delta, int noc);
//Erzeugt Normalverteilten Zufallsvektor der L?nge noa
void gausssample(double* temp, int noa);

void bedinge(int age_block, int noa, double* theta, double* L, double* thetatemp, int k, double* A, double* b);

void blockupdate_a2(int swit, int age_block, double kappa, double kappa2, double delta, int noa, int nop, double** ksi, double &my, double* theta, double* theta2, double* phi, double* psi, double* thetatemp, int vdb)
{
  int bandw=age_block*2;
  int* nwerte;

  int n=2*noa;

  double* L = new double[(bandw+1)*n];
  double* Lback = new double[(bandw+1)*n];
  if (swit<0)
  {
    nwerte=new int[noa];
    berechneBcohort(nwerte, thetatemp,ksi, my, phi, psi, -swit, nop,vdb,delta, noa);
    berechneQcohort2(nwerte, L, age_block, kappa, noa, delta, kappa2);
    delete[] nwerte;
  }
  else
  {
    berechneB(swit, thetatemp,ksi, my, phi, psi, noa, nop,vdb,delta);
    berechneQ2(L, age_block, kappa, noa, nop, delta, kappa2);
  }

  double* thetaschlangetemp= new double[n];
  double* thetaschlange= new double[n];

  //  for (int i=0; i<noa; i++)
  //    {
  //      thetaschlangetemp[(i*2)]=thetatemp[i];
  //      thetaschlangetemp[(i*2)+1]=0.0;
  //    }
  // mxschreibe(L, n,3);

  //oder so?
  for (int i=0; i<noa; i++)
  {
    thetaschlangetemp[(i*2)+1]=thetatemp[i];
    thetaschlangetemp[(i*2)]=0.0;
  }


  L = cholesky(n, L, bandw);

  for (int i=0; i< (bandw+1)*n; i++)
  {
    Lback[i]=L[i];
  }

  //Lv=b
  loese2(L, thetaschlangetemp, n, bandw);
  //L^T my=v
  loese(L,thetaschlangetemp, n, bandw);

  //z
  gausssample(thetaschlange,n);

  //L^T y=z
  loese(L,thetaschlange,n,bandw);

  //x=my + y
  double sum=0;
  for (int i=0; i< n; i++)
  {
    thetaschlange[i]=thetaschlange[i]+thetaschlangetemp[i];
    sum+=thetaschlange[i];
  }

  //Zentrieren ???
  // sum=sum/(double)n;

  // for (int i=0; i< noa; i++)
  //   {
  //     thetaschlange[i]=thetaschlange[i]-sum;
  //   }
  //my=my+sum;

  double* A=new double[2*n];
  double* b=new double[2];
  b[0]=0.0;
  b[1]=0.0;
  //  b[2]=0.0;
  for (int i=0; i< noa; i++)
  {
    A[2*i]=0.0;
    A[2*i+1]=1.0;
    A[2*(noa+i)]=1.0;
    A[2*(noa+i)+1]=0.0;
    //A[(4*noa)+(2*i)]=0.0;
    //  A[(4*noa)+(2*i)+1]=0.0;
  }
  //A[(4*noa)]=1.0;
  //A[(4*noa)+1]=-1.0;
  //   A[(6*noa)-2]=-1.0;
  //   A[(6*noa)-1]=1.0;
  //  mxschreibe(A, n,2);
  //  cout << thetaschlange[n-1]-thetaschlange[n-2] << endl;
  bedinge(bandw, n, thetaschlange, Lback, thetaschlangetemp, 2, A, b);

  //cout << thetaschlange[n-1]-thetaschlange[n-2] << endl<< endl;

  for (int i=0; i< noa; i++)
  {
    theta2[i] = thetaschlange[(2*i)+1] - thetaschlange[2*i];
    theta[i] = thetaschlange[(2*i)+1];
  }


  //   cout << endl;
  delete[] L;
  delete[] Lback;
  delete[] A;
  delete[] b;
  delete[] thetaschlange;
  delete[] thetaschlangetemp;

  return;

}

