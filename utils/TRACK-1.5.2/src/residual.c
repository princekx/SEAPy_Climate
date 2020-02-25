#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "statistic.h"
#include "p_vecs.h"
#include "sqt.h"



/* function to compute residuals for M-smoothing */


void residual(double (*cvalln)(), struct dpt *dt , int dtn, double tsm, double *plt , float *wght, int ht, int im, struct cvecs *pcom, int nni, LEAF *lf, VEC *gv, int nlf, int sqt_c)

{

   int din=1;

   float smm;
   double stp;

   stp = tsm;

   printf("Do you want to apply a different amount of smoothing for \r\n"
          "checking for robustness, 'y' or 'n'?                     \n\n");

   scanf("\n");
   if(getchar() == 'y'){
      printf("What is the new smoothing parameter (>0)?\n\n");
      scanf("%f", &smm);
      tsm = smm;
   }

   if(tsm < 0.){
      printf("****ERROR****, invalid smoothing parameter for residual estimation \r\n"
             "               using original value.                               \n\n");
      tsm = fabs(stp);
   }


/*   printf("Do you want data point inclusion or exclusion in residual calculation, 'i' or 'e'\n\n");
   scanf("\n");
   din = (getchar() == 'i') ? 1 : 0; */

   switch(im){
      case 0: case 2: case 3:
          if(pcom) (*cvalln)(dt, dtn, tsm, plt, wght, 0, ht, 0, din, pcom, nni);
          else if(lf) (*cvalln)(dt, dtn, tsm, NULL, wght, 0, ht, 0, din, lf, gv, nlf, 1);
          else (*cvalln)(dt, dtn, tsm, NULL, wght, 0, ht, 0, din);
          if(im == 0 || im == 2) break;
      case -1:
          if(pcom) (*cvalln)(dt, dtn, tsm, plt, wght, -1, ht, 0, din, pcom, nni);
          else if(lf) (*cvalln)(dt, dtn, tsm, NULL, wght, -1, ht, 0, din, lf, gv, nlf, 1, sqt_c);
          else (*cvalln)(dt, dtn, tsm, NULL, wght, -1, ht, 0, din);
          break;
      default:
          printf("***ERROR***, illegal switch for type of statistic, see sp_kernal_estimate\n\n");
          exit(1);
   }

   return;

}
