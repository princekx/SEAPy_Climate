#include <Stdio.h>
#include <Math.h>
#include "complex.h"
#include "m_values.h"


/* Blackman-Tukey estimate of the power spectra based on  
   Kay, Modern Spectral Estimation: Theory and Application. */

int powi(int , int );
void dfft(complex * , int , int , int );

void blackman_tukey(complex *timser, complex *rcor, int m, double *w, int nexp, double *pwr, int nf, int itrig)

{

   int i;

   int m1=m+1, mm1=m-1;
   int n2=nf/2;

   double arg;

   complex zz, znew;

   if(w){

      cmx(*(w + m), *rcor, (timser + m));
 
      for(i=0; i <m; i++){
          zz = *(rcor + i + 1);
          cmx(*(w + m1 + i), zz, (timser + m1 + i));
          conjg(*(rcor + i + 1), &znew);
          cmx(*(w + m1 + i), znew, (timser + mm1 - i));

      }

   }

   else {

      *(timser + m) = *rcor;
      for(i=0; i <m; i++){
          zz = *(rcor + i + 1);
          *(timser + m1 + i) = zz;
          conjg(zz, &znew);
          *(timser + mm1 - i) = znew;

      }

   }

/* zero pad time series */

   for(i=2*m+1; i < nf; i++) comp(0.0, 0.0, (timser+i));

   dfft(timser, nexp, 1, itrig);

   for(i=0; i < nf; i++){

      arg = 2.0 * FPI * (float) m * (float)i / (float) nf;
      comp(cos(arg), sin(arg), &zz);

      cmul(*(timser + i), zz, (timser + i)); 

   }

   

   for(i=0; i<n2; i++){

       *(pwr + i + n2) = (timser + i)->real / n2;
       *(pwr + i) = (timser + i + n2)->real / n2;


   }

   return;

}
