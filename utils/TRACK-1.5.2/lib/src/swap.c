#include <Stdio.h>

/* function to swap constraints in basis and operators */

void swap(double *n, double *nn1, double *b, int *acon, int iq, int q, int dimv, int dimn)

{
    int i, ii, jj;

    double temp, *t1, *t2;
 

    ii = iq-1;
    jj = q-1;

    temp = acon[ii];
    acon[ii] = acon[jj];
    acon[jj] = temp;

    temp = b[ii];
    b[ii] = b[jj];
    b[jj] = temp;

    t1 = n + (iq-1) * dimv;
    t2 = n + (q-1) * dimv;

    for(i=0; i < dimv; i++) {

         temp = *t1;
         *t1 = *t2;
         *t2 = temp;
         ++t1;
         ++t2;

    }

    t1 = nn1 + (iq-1)*dimn;
    t2 = nn1 + (q-1)*dimn;

    for(i=0; i < q; i++){     /* swap columns */

         temp = *t1;
         *t1 = *t2;
         *t2 = temp;
         ++t1;
         ++t2;

    }

    for(i=0; i < q; i++){     /* swap rows */

        ii = i*dimn-1;
        t1 = nn1 + ii + iq;
        t2 = nn1 + ii + q;

        temp = *t1;
        *t1 = *t2;
        *t2 = temp;

    }

    return;

}
