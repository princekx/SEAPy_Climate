#include <Stdio.h>
#include <stdarg.h>
#include "mem_er.h"

#define  SMALL  1.0e-20

/* function to update the hessian type matrix h used in the conjugate
   gradient optimization procedure when a constraint is added or deleted 
   from the basis. */

void update_h(double *hh, double *n, int dimv, int dir, ...)

{

       double nth[dimv], hn[dimv];
       double np;
       double denom;
       double *pp;

       int i, j, ii, jj;

       va_list argptr;
       va_start(argptr, dir);

       if(dir == 1){

          for(i=0; i < dimv; i++){

              hn[i] = nth[i] = 0.;

              ii = i * dimv;

              for(j=0; j < dimv; j++){

                 jj = j*dimv + i;
                 np = *(n + j);

                 hn[i] += *(hh+jj) *np;
                 nth[i] += *(hh + ii + j) * np;

              }

          }

          denom = 0.;

          for(i=0; i < dimv; i++) denom += *(n+i) * hn[i]; 

          if(denom < SMALL) denom = SMALL;

          for(i=0; i < dimv; i++){

              for(j=0; j < dimv; j++) *(hh + j*dimv + i) -= hn[i] * nth[j] / denom;


          }

       }

       else{

          pp = va_arg(argptr, double * );
      
          for(i=0; i < dimv; i++){

              hn[i] = nth[i] = 0.;

              ii = i * dimv;

              for(j=0; j < dimv; j++){

                 jj = j*dimv + i;
                 np = *(n + j);
                 hn[i] += *(pp+jj) *np;
                 nth[i] += *(pp + ii + j) * np;

              }

          }

          denom = 0.;

          for(i=0; i < dimv; i++) denom += *(n+i) * hn[i];

          if(denom < SMALL) denom = SMALL;

          for(i=0; i < dimv; i++){

              for(j=0; j < dimv; j++) *(hh + j*dimv + i) += hn[i] * nth[j] / denom;

          }

       }

/*printf("H\n");
for(i=0; i<dimv; i++){
for(j=0; j < dimv; j++) printf("%f  ", *(hh + j*dimv + i));
printf("\n");
}printf("\n\n"); */


       va_end(argptr);


       return;

}
