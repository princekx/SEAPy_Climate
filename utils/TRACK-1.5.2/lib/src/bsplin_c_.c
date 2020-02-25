#include <Stdio.h>

/* function to call the fortran subroutine which evaluates the b-splines */


#ifdef NOUNDERSCORE

void bsplin();

void bsplin_c(double *t, int *n, int *k, double *x, int *l, double *h)

{

    double hh[*k];

    bsplin(t, n, k, x, l, h, hh);

    return;

}

#else


void bsplin_();

#ifdef   G77
void bsplin_c__(double *t, int *n, int *k, double *x, int *l, double *h)
#else
void bsplin_c_(double *t, int *n, int *k, double *x, int *l, double *h)
#endif
{

    double hh[*k];

    bsplin_(t, n, k, x, l, h, hh);

    return;

}

#endif
