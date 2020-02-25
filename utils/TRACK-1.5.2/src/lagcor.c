#include <Stdio.h>
#include <Math.h>
#include "complex.h"

#define TOLMISS  1.0e-4


/* taken from Kay, S.M, 1989: Modern Spectral Estimation: Theory
   and Application, and adapted for missing data to produce the 
   unbiased autocorrelation or crosscorrelation.                 */

int missing(float, float, int );

void  lagcor(int nn, int nlag, complex *x, complex *y, complex *r, float fmiss, int imiss, int icmp)
{

   int i, j;
   int nk, nw;

   complex sum, conj, cm;

   if(imiss){

      for(i=0; i<nlag; i++){

          nk = nn - i;
          comp(0.0, 0.0, &sum);
          nw = 0;
          for(j=0; j < nk; j++){


             if(!missing((x + j)->real, fmiss, icmp) && 
                !missing((x + j + i)->real, fmiss, icmp)){

                conjg(*(x + j), &conj);
                cmul(conj, *(y + j + i), &cm); 
                cadd(sum, cm, &sum);
                ++nw;
             }

          }

          cmx(1.0 / (float) nw, sum, (r + i));
  
      }

   }

   else {

      for(i=0; i<nlag; i++){

          nk = nn - i;
          comp(0.0, 0.0, &sum);
          nw = 0;
          for(j=0; j < nk; j++){

             conjg(*(x + j), &conj);
             cmul(conj, *(y + j + i), &cm); 
             cadd(sum, cm, &sum);
             ++nw;

          }

          cmx(1.0 / (float) nw, sum, (r + i));
   
      }

   }

   return;

}
