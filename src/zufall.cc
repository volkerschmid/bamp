#include <R.h>
#include <Rmath.h>
#include <time.h>
#include "l.h"



// extern "C" {btpec_(int& N,double* PP,double& ISEED,int& JX);}
extern "C" {long ignbin(long n,float pp);}

//Erzeugt gleichverteilte Zufallszahlen auf [0,1[ mit linearem
//Kongruenzgenerator mit Multiplikator 25214903917 und Modulus 2**48


double nulleins(){
  GetRNGstate();
  double temp=runif(0.0, 1.0);
  PutRNGstate();
  return temp;
  }

//Erzeugt diskrete gleichverteilte Zufallszahlen auf 1,...,k

int einsk(int k){
  return int(floor(nulleins()*k)) + 1;
}

//Erzeugt gleichverteilte Zufallszahlen auf [-1,1]

double reins()
{
  return nulleins()*2.0 - 1.0;
}


//Erzeugt exonentialverteilte Zufallszahlen

double expo(double lambda)
{
  double u,v;
  do{
    u = nulleins();
  }while (u==0.0);   // Sicherstellen dass u>0

  v = (-1 / lambda) * log(u);
  return v;
}



//Erzeugt poissonverteilte Zufallszahlen

int poisson(double lambda)
{
	double sum = 0.0;
	int zaehler = 0;
	while(sum <= 1)
	{
		sum = sum + expo(lambda);
		zaehler++;
	}
	return (zaehler - 1);
}



//Erzeugt standardnormalverteilte Zufallszahlen
double zahl(double mu, double sigma)
{
  GetRNGstate();
  double temp=rnorm(0,1);
  PutRNGstate();

  return temp;
}
double zahl_old(double mu, double sigma)
{
	static int gespeichert=0;
	static double alteZahl;

	if(gespeichert==0)
	{
		gespeichert=1;
		double X, Y, r2;

		do
		{
			X=reins();
			Y=reins();
			r2=X*X+Y*Y;
		}while(r2>1.0);

		double norm=sqrt((-2.0*log(r2))/r2);
		alteZahl=Y*norm*sigma+mu;
		return X*norm*sigma+mu;
	}
	else
	{
		gespeichert=0;
		return alteZahl;
	}
}

//Erzeugt gleichverteilte Zufallszahlen auf [0,L]

double nullell(double l)
{
	return (l*nulleins());
}



//Erzeugt gleichverteilte Zufallszahlen auf [a,b]

double gleichab(double a, double b)
{
	return (nullell(b-a)+a);
}

//Erzeugt gleichverteilte Zufallszahlen auf [-0.5,0.5]

double halbeins()
{
  return nulleins()-0.5;
}






/* gamma verteilte Zufallszahlen; Bernardo Smith Notation */

// -------------------------------------
double RNDGAM(double a, double b){
// -------------------------------------
int accept=0;
double c,d,u,v,w,x,y,z;
if(a>1){                        /*   Algorithmus S.410 Devroye */
  c=a-1;
  d=3*a-3/4;
  while(accept==0){
    u=nulleins();
    v=nulleins();
    w=u*(1-u);
    y=sqrt(d/w)*(u-0.5);
    x=c+y;
    if(x >= 0){
      z=64*w*w*w*v*v;
      if( z<= (1-(2*y*y/x))) accept=1;
      if (accept==0)
        if (log(z)<=2*((c*log(x/c))-y)) accept=1;
      }
  }
}
else{                               /* Fall: a<=1; Stuart's theorem */
x = (double(pow(nulleins(),(1/a)))*RNDGAM(a+1,1));
}
return x/b;


}

int rbinom(double p,int n)
{
  return (int)ignbin((long)n,(float)p);
}


long ignbin(long n,float pp)

{
static float psave = -1.0;
static long nsave = -1;
static long ignbin,i,ix,ix1,k,m,mp,T1;
static float al,alv,amaxp,c,f,f1,f2,ffm,fm,g,p,p1,p2,p3,p4,q,qn,r,u,v,w,w2,x,x1,
    x2,xl,xll,xlr,xm,xnp,xnpq,xr,ynorm,z,z2;

    if(pp != psave) goto S10;
    if(n != nsave) goto S20;
    if(xnp < 30.0) goto S150;
    goto S30;
S10:
/*
*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
*/
    psave = pp;
    p = min(psave,1.0-psave);
    q = 1.0-p;
S20:
    xnp = n*p;
    nsave = n;
    if(xnp < 30.0) goto S140;
    ffm = xnp+p;
    m = ffm;
    fm = m;
    xnpq = xnp*q;
    p1 = (long) (2.195*sqrt(xnpq)-4.6*q)+0.5;
    xm = fm+0.5;
    xl = xm-p1;
    xr = xm+p1;
    c = 0.134+20.5/(15.3+fm);
    al = (ffm-xl)/(ffm-xl*p);
    xll = al*(1.0+0.5*al);
    al = (xr-ffm)/(xr*q);
    xlr = al*(1.0+0.5*al);
    p2 = p1*(1.0+c+c);
    p3 = p2+c/xll;
    p4 = p3+c/xlr;
S30:
/*
*****GENERATE VARIATE
*/
    u = nulleins()*p4;
    v = nulleins();
/*
     TRIANGULAR REGION
*/
    if(u > p1) goto S40;
    ix = xm-p1*v+u;
    goto S170;
S40:
/*
     PARALLELOGRAM REGION
*/
    if(u > p2) goto S50;
    x = xl+(u-p1)/c;
    v = v*c+1.0-ABS(xm-x)/p1;
    if(v > 1.0 || v <= 0.0) goto S30;
    ix = x;
    goto S70;
S50:
/*
     LEFT TAIL
*/
    if(u > p3) goto S60;
    ix = xl+log(v)/xll;
    if(ix < 0) goto S30;
    v *= ((u-p2)*xll);
    goto S70;
S60:
/*
     RIGHT TAIL
*/
    ix = xr-log(v)/xlr;
    if(ix > n) goto S30;
    v *= ((u-p3)*xlr);
S70:
/*
*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
*/
    k = ABS(ix-m);
    if(k > 20 && k < xnpq/2-1) goto S130;
/*
     EXPLICIT EVALUATION
*/
    f = 1.0;
    r = p/q;
    g = (n+1)*r;
    T1 = m-ix;
    if(T1 < 0) goto S80;
    else if(T1 == 0) goto S120;
    else  goto S100;
S80:
    mp = m+1;
    for(i=mp; i<=ix; i++) f *= (g/i-r);
    goto S120;
S100:
    ix1 = ix+1;
    for(i=ix1; i<=m; i++) f /= (g/i-r);
S120:
    if(v <= f) goto S170;
    goto S30;
S130:
/*
     SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))
*/
    amaxp = k/xnpq*((k*(k/3.0+0.625)+0.1666666666666)/xnpq+0.5);
    ynorm = -(k*k/(2.0*xnpq));
    alv = log(v);
    if(alv < ynorm-amaxp) goto S170;
    if(alv > ynorm+amaxp) goto S30;
/*
     STIRLING'S FORMULA TO MACHINE ACCURACY FOR
     THE FINAL ACCEPTANCE/REJECTION TEST
*/
    x1 = ix+1.0;
    f1 = fm+1.0;
    z = n+1.0-fm;
    w = n-ix+1.0;
    z2 = z*z;
    x2 = x1*x1;
    f2 = f1*f1;
    w2 = w*w;
    if(alv <= xm*log(f1/x1)+(n-m+0.5)*log(z/w)+(ix-m)*log(w*p/(x1*q))+(13860.0-
      (462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0+(13860.0-(462.0-
      (132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0+(13860.0-(462.0-(132.0-
      (99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0+(13860.0-(462.0-(132.0-(99.0
      -140.0/w2)/w2)/w2)/w2)/w/166320.0) goto S170;
    goto S30;
S140:
/*
     INVERSE CDF LOGIC FOR MEAN LESS THAN 30
*/
    qn = pow(q,(double)n);
    r = p/q;
    g = r*(n+1);
S150:
    ix = 0;
    f = qn;
    u = nulleins();
S160:
    if(u < f) goto S170;
    if(ix > 110) goto S150;
    u -= f;
    ix += 1;
    f *= (g/ix-r);
    goto S160;
S170:
    if(psave > 0.5) ix = n-ix;
    ignbin = ix;
    return ignbin;
}


double normal(double m, double s)
{
  return (zahl(0.0,1.0)*sqrt(s)+m);
}
