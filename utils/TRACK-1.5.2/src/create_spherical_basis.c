#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "mem_er.h"
#include "m_values.h"
#include "grid.h"

#define TOLG  1.0e-4

extern float period;

void pnm(double * , double * , double , int , int );

int create_spherical_basis(GRID *gcur, int ntr, int ider)

{


    int i, j;
    int ntr1 = ntr + 1;
    int n2=0;
    long int nd1, nd2;
    int nx;


    double *tmp=NULL, *tmp2=NULL;
    double sn, cn, arg, argn;

/* memory for spherical harmonic coefficients */

/* determine number of Lengendre functions for truncation */

    for(i=0; i <= ntr; i++)n2 += (ntr1 - i);

    gcur->nleng = n2;

    nd1 = (long int)(gcur->iy) * (long int)(n2);

    gcur->aleng = (double *)calloc(nd1, sizeof(double));
    mem_er((gcur->aleng == NULL)? 0 : 1, nd1 * sizeof(double));

    if(ider){
      gcur->dleng = (double *)calloc(nd1, sizeof(double));
      mem_er((gcur->dleng == NULL)? 0 : 1, nd1 * sizeof(double));
    }

    nx = (fabs(*(gcur->xgrid) + period - *(gcur->xgrid + gcur->ix - 1)) < TOLG) ? gcur->ix - 1 : gcur->ix;


    nd2 = (long int) nx * (long int) ntr1;

    gcur->sico = (complex *)calloc(nd2, sizeof(complex));
    mem_er((gcur->sico == NULL) ? 0 : 1, nd2 * sizeof(complex));



    for(i=0; i < nx; i++){

       arg = *(gcur->xgrid + i) * FP_PI;
       if(arg < 0.) arg = 0.;
       else if (arg > FPI2) arg = FPI2;

       for(j=0; j < ntr1; j++){

           argn = arg * (float) j;
           sincos(argn, &sn, &cn);
           comp(cn, sn, (gcur->sico+i * ntr1 + j));

       }

    }



    for(i=0; i < gcur->iy; i++){

       tmp = gcur->aleng + i * n2;
       tmp2 = gcur->dleng + i * n2;

       arg = FP_PI2 - *(gcur->ygrid + i) * FP_PI;

       if(arg < 0.) arg = 0.;
       else if(arg > FPI) arg = FPI;

       arg = cos(arg);

       pnm(tmp, tmp2, arg, ntr, ider);


    }


    return nx;    

}
