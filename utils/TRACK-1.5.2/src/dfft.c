/* Discrete fast Fourier transform based on the Fortran program of Kay, 
   Modern Spectral Estimation: Theory and Application.                  */

#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "complex.h"
#include "m_values.h"
#include "mem_er.h"

int powi(int , int );

complex *zsave1=NULL, *zsave2=NULL;    
complex *yt=NULL;


void dfft(complex *xt, int nexp, int inv, int itrig)

{

   static int n=0;
   
   int ic;

   int i, j, l, k;
   int nn, ld, nd; 
   int np, nq;
   int nnex=0;

   double arg, ss, cc;
   double aa;
   
   complex w, sav;
   complex *zsave=NULL;

   complex zz;


   if(inv != 1){

      if(inv != -1){

         printf("****ERROR****, transform direction must be either foreward (1)\r\n"
                "               or backward (-1).                              \n\n");
         exit(1);

      }

   }

   aa = -(FPI2 * (float)inv);

   if(itrig > 0) n = powi(2, nexp);


   if(!yt){

     yt = (complex *)calloc(n, sizeof(complex));
     mem_er((yt == NULL) ? 0 : 1, n * sizeof(complex));

   }

   if(!zsave1 && inv == 1){
      nnex = nexp*powi(2, nexp-1);
      zsave1 = (complex *)calloc(nnex, sizeof(complex));
      mem_er((zsave1 == NULL) ? 0 : 1, nnex * sizeof(complex));

   }

   if(!zsave2 && inv == -1){
      nnex = nexp*powi(2, nexp-1);
      zsave2 = (complex *)calloc(nnex, sizeof(complex));
      mem_er((zsave2 == NULL) ? 0 : 1, nnex * sizeof(complex));

   }


   if(itrig < 0){free(zsave1); free(zsave2); free(yt); return;}

   if(inv == 1) zsave = zsave1;
   else zsave = zsave2;

   *yt = *xt;

   for(i=1; i < n; i++){

      nn = i;
      j = 0;
      l = n;

      for(k=0; k < nexp; k++){

         l /= 2;
         if(nn < l) continue;
         j += powi(2, k);
         nn -= l;
      }

      *(yt + i) = *(xt + j);
 

   }


   if(itrig > 0){     /* compute all trigonometric functions  */

      ld = 1;
      nd = n;
      ic = 0; 

      for(i=0; i< nexp; i++){

          ld *= 2;
          nd /= 2;

          for(j=0; j < nd; j++){

              for(k=0; k < ld/2; k++){

                  arg = aa * ((float)k / (float)ld);
                  sincos(arg, &ss, &cc);
                  comp(cc, ss, (zsave + ic));

                  ++ic;

              }

          }
        
      }

   }


   if(!zsave){
      printf("****ERROR****, no trigonometric evaluations avaliable\n\n");
      return;
   }

   ld = 1;
   nd = n;
   ic = 0;

   for(i=0; i< nexp; i++){

       ld *= 2;
       nd /= 2;

       for(j=0; j < nd; j++){

          for(k=0; k < ld/2; k++){

              np = k + ld * j;
              nq = np + ld / 2;
              w = *(zsave + ic);
              ++ic;
              cmul(w, *(yt+nq), &zz);
              cadd(*(yt+np), zz, &sav);
              csub(*(yt+np), zz, (yt+nq));
              *(yt + np) = sav;
 
          }

       }

    }

    if(inv == 1){

       for(i=0; i<n; i++) *(xt + i) = *(yt + i);

    }

    else {

       for(i=0; i<n; i++) {
          (xt + i)->real = (yt + i)->real / n;
          (xt + i)->imag = (yt + i)->imag / n;

       }


    }


    return;

}
