#include <Stdio.h>
#include "m_values.h"
#include "tessel.h"

#define   CNT   1      /* set to 1 if triangle centroids are required */

/* function to decompose an octant of the sphere using a projection
   onto a trianglular facet of an octahedron.                       */

int powi(int , int );


void tesselate_facet(int ll, struct cpt *cpos)

{

   int nn,nn1,n1=ll;
   int i, j, jj;

   float x, y;

   struct cpt *cp;

   nn1 = powi(2, n1);
   nn = 2 * nn1;

   cp = cpos;

   for(i=1; i <= nn1; i++) {

      jj = nn - 2 * i + 1;

      for(j=0; j < jj; j++){

          x = (float) (i + j);

          if(CNT) y = (float) (i - 1) + (float) (2 - (1 + powi(-1, j)) / 2) / 3.0;
          else y = 0.5 + (float) (i - 1);

          cp->yct = FPI * y / (nn * FP_PI);

          cp->xct = (FPI * (SQ_R3 * x - y)) / (4. * (SQ_R3 * nn1 - y));
          cp->xct /= FP_PI;

          ++cp;

      }

   }

   return;

}
