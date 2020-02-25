#include <Stdio.h>

/* function to perform operations to correct initial feasable vector
   which violates the constraints for the conjugate gradient optimization. */

void init_vec(double *mv, double *nn1, double *n, double *alpha, int dimv, int dimnc, int q)

{

       int i, j;

       double nl[q];
       double sum;

       for(i=0; i < q; i++){

           nl[i] = 0.;

           for(j=0; j < q; j++) nl[i] += nn1[j*dimnc+i]*alpha[j];

       }

       for(i=0; i < dimv; i++){

          sum = 0.;

          for(j=0; j < q; j++) sum += n[j*dimv+i]*nl[j];

          mv[i] -= sum;

       }


       return;

}
