#include <R.h>
#include <Rmath.h>
#include "l.h"


//Erzeugt gleichverteilte Zufallszahlen auf [0,1[ mit linearem
//Kongruenzgenerator mit Multiplikator 25214903917 und Modulus 2**48


double nulleins(){
  GetRNGstate();
  double temp=runif(0.0, 1.0);
  PutRNGstate();
  return temp;
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
        if (std::log(z)<=2*((c*std::log(x/c))-y)) accept=1;
      }
  }
}
else{                               /* Fall: a<=1; Stuart's theorem */
x = (double(std::pow(nulleins(),(1/a)))*RNDGAM(a+1,1));
}
return x/b;
}

double normal(double m, double s)
{
  GetRNGstate();
  double temp=rnorm(0,1)*sqrt(s)+m;
  PutRNGstate();
  return(temp);
  }
