#include <Stdio.h>

/* function to compute the intersection of two ordered sets of integers 
   for use in picking out tracks for statistical analysis               */

int tr_overlap(int *nn, int f1, int f2, int ts, int tf)

{

   int in=0;

   if(f1 <= ts && f2 >= tf) { in = 1; *nn = tf - ts; }

   else if(ts < f1 && tf > f2) { in = 1; *nn = f2 - f1 + 2;}

   else if(f1 <= ts && ts <= f2) {

      in = 1;
      *nn = f2 - ts + 1;

    }

    else if(f1 <= tf && tf <= f2) {

       in = 1;
       *nn = tf - f1 + 1;

    }

    else {in = 0; *nn = 0; }

   return in;

} 
