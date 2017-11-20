
//multipliziert A mit B; A ist (a x n), B ist (n x b)  
void multipliziere(double* A, double* B, int a, int n, int b, double* ergebnis){

  //zu machen: wenn a=1 oder b=1 oder n=1 !

// if (a==1)
//   {
//     if (n==1)
//       {
// 	if (b==1)
// 	  {
// 	    ergebnis[0]=A[0]*B[0];
// 	  }
// 	else
// 	  {
// 	    for (int i=0; i< b; i++)
// 	      {
// 		ergebnis[i]=A[0]*B[i];
// 	      }
// 	  }
//       }
//     else
//       {
// 	if (b==1)
// 	  {
// 	    ergebnis[0]=0.0;
// 	    for (int i=0;i<n;i++)
// 	      {
// 		ergebnis[0]+=A[i]*B[i];
// 	      }
// 	  }
// 	else
// 	  {
// 	    for (int i=0; i<b; i++)
// 	      {
// 		ergebnis[i]=0.0;
// 		for (int j=0; j<n; j++)
// 		  {
// 		    ergebnis[i]+=A[j]*B[i*n+j];
// 		  }
// 	      }
// 	  }
//       }
//   }
// else
//   {
//     if (n==1)
//       {
// 	if (b==1)
// 	  {
// 	    for (int i=0; i<a; i++)
// 	      {
// 		ergebnis[i]=A[i]*B[0];
// 	      }
// 	  }
// 	else
// 	  {
// 	    for (int i=0; i<a; i++)
// 	      {
// 		for (int j=0; j<b; j++)
// 		  {
// 		    ergebnis[i*b+j]=A[i]*B[j];
// 		  }
// 	      }
// 	  }
//       }
//     else
//       {
// 	if (b==1)
// 	  {
// 	    for (int i=0; i<a; i++)
// 	      {
// 		ergebnis[i]=0.0;
// 		for (int j=0; j<n; j++)
// 		  {
// 		    ergebnis[i]+=A[i*n+j]*B[j];
// 		  }
// 	      }
// 	  }
// 	else

// 	  {
		


  //double* ergebnis=new double[a*b];
 for (int i=0; i<a; i++)
  {
    for (int j=0; j<b; j++)
      {
	ergebnis[i*b+j]=0.0;
	for (int k=0; k<n; k++)
	  {
	    ergebnis[i*b+j]=ergebnis[i*b+j]+A[i*n+k]*B[k*b+j];
	  }
      }
  }
// 	  }
//       }
//   }

return;

}



// BERECHNET A-1, k ist matrixlaenge
void invers(double* A,int k)
{

double* ergebnis=new double[k*k];
if (k==1)
  {
    ergebnis[0] = 1.0/A[0];
   
  }
if (k==2)
  {
    double det = A[0]*A[3]-A[1]*A[2];
    ergebnis[0] = A[3]/det;
    ergebnis[1] = 0.0-A[1]/det;
    ergebnis[2] = 0.0-A[2]/det;
    ergebnis[3] = A[0]/det;
  }
for (int i=0; i< k*k; i++)
 {
   A[i]=ergebnis[i];
 }

delete[] ergebnis;
return;

}


int mxcheck(int n, int** matrix)
{
int zs=0;
for (int i=0; i<n; i++)
  {
    zs=0;
    for (int j=0; j<n; j++)
      {
	zs+=matrix[i][j];
	if (matrix[i][j]!=matrix[j][i])
	  {
	    return 1;
	  }
      }
    if (zs != 0)
      {
        return 1;
      }
  }
return 0;
}
