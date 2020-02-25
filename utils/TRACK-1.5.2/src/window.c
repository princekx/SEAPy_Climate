#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "m_values.h"
#include "mem_er.h"


/* function to provide windows for PSD calculation, taken mostely from the
   Spectra package, Numerical Recipies and Kay (Modern Spectral Estimation) */


double *window(int m, int win_id, double *sum)

{

   int i, ii;
   int mm = 2 * m;

   double *w=NULL;

   *sum = 0.0;

   w = (double *)calloc(mm+1, sizeof(double));
   mem_er((w == NULL)? 0 : 1, (mm+1) * sizeof(double));

   if(win_id == 0){      /* Bartlett window */

      for(i=0; i <= mm; i++) {

         ii = i - m;
         *(w+i) = 1. - (((ii > 0) ? (float)ii : -(float)ii) / (float)m );

      }

   }

   else if(win_id == 1){   /* Hanning Window  */

      for(i=0; i <= mm; i++) 
         *(w+i) = 0.5 * (1.0 + cos(FPI * (float)(i - m) / (float) m));

   }

   else if(win_id == 2){   /* Hamming Window */

      for(i=0; i <= mm; i++) 
         *(w+i) = 0.54 + 0.46 * cos(FPI * (float)(i - m) / (float) m);

   }

   else 

       printf("****ERROR****, no windowing available for the chosen identifier\n\n");

   for(i=0; i <=mm; i++) *sum += *(w+i);
 
   return w;

}
