#include <Stdio.h>
#include <Math.h>


/* function to bracket a maximum for the line optimization,
   taken from Numerical Recipies.                             */

#define  GOLD               1.618034
#define  TINY               1.0e-20
#define  MAX(a, b)          (double)((a) > (b) ? (a) : (b))
#define  SIGN(a, b)         (double)((b) > 0.0 ? fabs(a) : -fabs(a))
#define  SHFT(a, b, c, d)   (a)=(b); (b)=(c); (c)=(d);

double func(double * , double * , int , int , int );

void bracket(double *a, double *b, double *mvec, double *fa, double *fb, double *suma, double *sumb, double *h, double lam, int dimv, int type)

{

    int i;

    double c, u, r, q, fc, fu, dum;
    double sumc, sumu, mv[dimv], dvec[dimv];

/*    if(*fb < *fa){
       SHFT(dum, *a, *b, dum)
       SHFT(dum, *fa, *fb, dum)
       SHFT(dum, *suma, *sumb, dum)
    } */

    c = (*b) + GOLD*(*b - *a);

    for(i=0; i < dimv; i++) mv[i] = mvec[i] + c * h[i];

    fc = func(mv, dvec, dimv, type, 1);

    sumc = 0.;

    for(i=0; i < dimv; i++) sumc += dvec[i] * h[i];

    while(*fb < fc){

         r = (*b - *a) * (*fb - fc);
         q = (*b - c) * (*fb - *fa);
         u = (*b) - ((*b - c) * q - (*b - *a) * r)/
             (2.0 * SIGN(MAX(fabs(q-r), TINY), q-r));

         if((*b - u) * (u - c) > 0.0) {

           for(i=0; i < dimv; i++) mv[i] = mvec[i] + u * h[i];

           fu = func(mv, dvec, dimv, type, 1);

           sumu = 0.;

           for(i=0; i < dimv; i++) sumu += dvec[i] * h[i];

           if(fu > fc){

              *a = *b;
              *b = c;
              *fa = *fb;
              *fb = fc;
              *suma = *sumb;
              *sumb = sumc;  
              return;

            }

            else if(fu < *fb){

               *b = u;
               *fb = fu;
               *sumb = sumu;
               return;

            }

            u = c + GOLD * (c - *b);

            for(i=0; i < dimv; i++) mv[i] = mvec[i] + u * h[i];

            fu = func(mv, dvec, dimv, type, 1);

            sumu = 0.;

            for(i=0; i < dimv; i++) sumu += dvec[i] * h[i];

         }

         else if((c - u) * (u - lam) > 0.0){


            for(i=0; i < dimv; i++) mv[i] = mvec[i] + u * h[i];

            fu = func(mv, dvec, dimv, type, 1);

            sumu = 0.;

            for(i=0; i < dimv; i++) sumu += dvec[i] * h[i];

            if(fu > fc){

               *b = c;
               c = u;
               u = c + GOLD * (c - *b);
               *fb = fc;
               fc = fu;
               
               for(i=0; i < dimv; i++) mv[i] = mvec[i] + u * h[i];

               fu = func(mv, dvec, dimv, type, 1);

               *sumb = sumc;
               sumc = sumu;

               sumu = 0.;

               for(i=0; i < dimv; i++) sumu += dvec[i] * h[i];

             }

           }

           else if((u - lam) * (lam - c) >= 0.0) {

               u = lam;

               for(i=0; i < dimv; i++) mv[i] = mvec[i] + u * h[i];

               fu = func(mv, dvec, dimv, type, 1);

               sumu = 0.;

               for(i=0; i < dimv; i++) sumu += dvec[i] * h[i];

           }

           else{

              u = c + GOLD * (c - *b);

              for(i=0; i < dimv; i++) mv[i] = mvec[i] + u * h[i];

              fu = func(mv, dvec, dimv, type, 1);

              sumu = 0.;

              for(i=0; i < dimv; i++) sumu += dvec[i] * h[i];

           }

           SHFT(*a, *b, c, u)
           SHFT(*fa, *fb, fc, fu);
           SHFT(*suma, *sumb, sumc, sumu);

    }

}
