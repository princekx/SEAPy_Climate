#include <Stdio.h>

/* fumction to raise an integer i to the power of the integer l */

int powi(int i, int l)

{

   int r=1, k=0;

   while(k++ < l) r=r*i;

   return r;

}
