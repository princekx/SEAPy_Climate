#include <Stdio.h>
#include <Math.h>

#define  TOLIND   1.0e-6
#define  R5       2.236067977
#define  R575     0.335410197
#define  BIW      0.937500000

#define  LVAL     1.0e+10

/* function to compute 1D kernel for teleconnection statistics */


double tele_kernel(double tst, double tval, double h, int ktyp)

{

   static double otval=LVAL;
   static double ker;
   double val, vabs;

   if(fabs(tval - otval)  < TOLIND) return ker;

   val = (tst - tval) / h;
   vabs = fabs(val);

   switch(ktyp) {
      case 0:        /* triangular  */
         if(vabs < 1.0) ker = (1.0 - vabs);
         else ker = 0.0;
         break;
      case 1:        /* rectangular */
         if(vabs < 1.0) ker = 0.5;
         else ker = 0.0;
         break;
      case 2:        /* Epanechnikov */
         if(vabs < R5) ker = R575 * (1.0 - 0.2 * val * val);
         else ker = 0.0;
         break;
      case 3:        /* Bi-weights  */
         if(vabs < 1.0) {
            val = (1.0 - val * val);
            ker = BIW * val * val;
         }
         else ker = 0.0;
         break;
   }

   otval = tval;

   return ker;

}
