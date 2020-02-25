#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>

#define  TOLMISS      1.0e-6

int missing(float val, float mval, int icmp)
{

    int ims=0;

    switch(icmp){
          case 0:
            if(fabs(val - mval) < TOLMISS) ims = 1;
             break;
          case 1:
             if(val < mval) ims = 1;
             break;
          case 2:
             if(val > mval) ims = 1;
             break;
          default:
             printf("****ERROR****, no missing value set for checking.\n\n");
             exit(1);
    }

    return ims;

}
