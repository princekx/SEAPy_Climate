#include <Stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mem_er.h"

#define  TOLLD   0.00001

/* function to update the nn1 matrix used in the conjugate gradient
   optimization procedure when a constraint is added to the basis.  */

int update_nn1(double *nn1, double *n ,int q, int dimv, int dimn, int dir, ...)

{
   int i, j, k, jj=0, qd=q*dimv, qt, qc;
   int ld=0;

   double *nn=NULL, *r=NULL, *tt=NULL, *p=NULL;
   double b4;
   double pp, pn=0., div;

   va_list argptr;
   va_start(argptr, dir);


   if(dir == 1){

      nn = (double * )calloc(dimn, sizeof(double));
      mem_er((nn == NULL) ? 0 : 1, dimn * sizeof(double));

      r = (double * )calloc(q, sizeof(double));
      mem_er((r == NULL) ? 0 : 1, q * sizeof(double));

      for(j=0; j < q; j++){

         nn[j] = 0.;

         jj = j*dimv;

         for(i=0; i < dimv; i++) nn[j] += n[jj + i] * n[qd+i];

      }

      for(i=0; i < q; i++){

         r[i] = 0.;

         for(j=0; j < q; j++) r[i] += *(nn1 + j*dimn +i) * nn[j];

      }

      for(i=0; i < dimv; i++){

         pp = 0.0;

         for(j=0; j < q; j++)pp -= n[j*dimv + i] * r[j];

         pp += n[qd+i];

         pn += pp*pp;

      }

      if(pn < TOLLD) ld = 1;

      else{

        b4 = nn1[q*dimn + q] = 1.0 / pn;

        for(i=0; i < q; i++) {

          nn1[q*dimn + i] = nn1[i*dimn + q] = - r[i] * b4;

          for(j=0; j < q; j++) nn1[j*dimn + i] += r[i]*r[j]*b4;

        }

      }

      free(nn);
      free(r);
     
   }

   else {

      p = va_arg(argptr, double * );

      qt = q - 1;
      qc = qt * dimn;
      div = 1.0 / nn1[qc + qt];

      for(i=0; i < qt; i++){

          for(j=0; j < qt; j++) 

             nn1[j*dimn + i] -= nn1[qc + i]*nn1[qc + j] * div;

      }

      nn = (double * )calloc(qt*dimv, sizeof(double));
      mem_er((nn == NULL) ? 0 : 1, qt*dimv * sizeof(double));

/* compute new p */

      for(i=0; i < dimv; i++){

          for(j=0; j < qt; j++){

             tt = nn + i*qt + j;

             *tt = 0.;

             for(k=0; k < qt; k++) *tt += nn1[k*dimn+j] * n[k*dimv+i];

          }

      }

      for(i=0; i < dimv; i++){

          for(j=0; j < dimv; j++){

              tt = p + i*dimv + j;
              *tt = (i == j) ? 1. : 0.;

              for(k=0; k < qt; k++) *tt -= n[k*dimv+j] * nn[i*qt+k];

          }

      }

     free(nn);

/* printf("P\n");
for(i=0; i < dimv; i++){
for(j=0; j < dimv; j++)printf("%f  ", p[j*dimv+i]);
printf("\n");
}printf("\n"); */


   }

/*printf("N-1\n");
for(i=0; i < q+1; i++){
for(j=0; j < q+1; j++)printf("%f  ", nn1[j*dimn+i]);
printf("\n");
}printf("\n"); */

   va_end(argptr);

   return ld;

}
