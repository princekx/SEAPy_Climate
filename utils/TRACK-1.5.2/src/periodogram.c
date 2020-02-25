/* Periodogram estimate of the power spectra based on the Fortran program
   of Kay, Modern Spectral Estimation: Theory and Application.            */

#include <Stdio.h>
#include <Math.h>
#include "complex.h"

int powi(int , int );
void dfft(complex * , int , int , int );


void periodogram(complex *xt, int n, double *ww, double wsum, int nexp, double *per, double *phase, int m, int itrig)
{

   int i, mm;

   float tnorm;

   double zabs;
/*   double wsum=0.0; */

/*   m = powi(2, nexp); */

   tnorm = (float)(n * (m / 2));

/* apply window to the data */

   if(ww){

      for(i=0; i<n; i++) {

         cmx(*(ww+i), *(xt+i), (xt+i));

/*         wsum += *(ww+i) * *(ww + i); */

      }

      tnorm = wsum  * (m / 2);


   }

/* zero padding */

   for(i=n; i<m; i++) comp(0.0, 0.0, (xt+i));


   dfft(xt, nexp, 1, itrig);


/* compute PSD */

   mm = m / 2;

   for(i=0; i<mm; i++){

       zabs = cabs2(*(xt+i));

       *(phase + i + mm) = atan2((xt+i)->imag, (xt+i)->real);

       *(per+i + mm) = zabs / tnorm;
       zabs = cabs2(*(xt + i + mm));

       *(phase + i) = atan2((xt+i+mm)->imag, (xt+i+mm)->real);

       *(per+i) = zabs / tnorm; 

   }


   return;

}
