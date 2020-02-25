#include <Stdio.h>
#include <Math.h>

#define TOLDIAG  1.0e-12


int cholesky(double **aa, int ncof)

{


    int i, j, k;

    double ad;

    for(i=0; i < ncof; i++){


        ad = *(aa[i] + i) = sqrt(*(aa[i] + i));

        if(ad < TOLDIAG) {

           printf("****ERROR****, normal equations matrix probably not \r\n"
                  "               positive to within set tolerance.    \r\n"
                  "               No least squares solution availavble.\n\n"); 
           return 1;

        }

        for(j=i+1; j < ncof; j++) *(aa[j] + i) /= ad;

        for(j=i+1; j < ncof; j++) {

            for(k=j; k < ncof; k++) *(aa[k] + j) -= *(aa[k] + i) * *(aa[j] + i);


        }


    }



    return 0;

}
