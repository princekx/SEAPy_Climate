#include <Stdio.h>
#include <Math.h>

#define MAXIT  20
#define TOLIT  1.0e-10
#define START  1.0e-5

double power_norm(double , double , int );


void regfal(double *kb, double *beta, double bmax)

{ 

    int nit=0;

    double k1, k2;
    double fa, fb;
    double kn=0.0, fn;
    double aa;

    aa = ((*kb - 1.0) / *kb);

    k1 = 1.0+START;
    k2 = *kb;

    if(k1 < 2.0 * *beta) *beta = bmax * k1;
    fa = power_norm(k1, *beta, 0) - aa;
    fb = power_norm(k2, *beta, 0) - aa;

    while(nit++ < MAXIT){


       kn = (k1 * fb - k2 * fa) / (fb - fa);



       if(kn < 2.0 * *beta) *beta = bmax * kn;

       fn = power_norm(kn, *beta, 0) - aa;


       if(fabs(fn) < TOLIT){

          *kb = kn;
          break;

       }

       if(fa * fn > 0.0){ k1 = kn; fa = fn;}
       else{ k2 = kn; fb = fn;}


    }

    if(nit >= MAXIT){

      printf("***WARNING***, maximum number of iterations exceeded trying to correct bandwidth support\n\n");

      *kb = kn;


    }

    *kb = kn;

    return;

}
