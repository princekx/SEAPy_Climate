#include <Stdio.h>
#include <Math.h>
#include <stdlib.h>
#include "mem_er.h"

#define SCTOL   1.0e-10

/* funtion to compute the Householder transformation of a matrix */

void householder(double **aa, int n, double *dd, double *ee)
{

    int l,k,j,i;
    double scale, hh, h, g, f;
 
    double **a, *d, *e;

    a = (double **)calloc(n, sizeof(double *));
    mem_er((a == NULL) ? 0 : 1, n * sizeof(double *));

    a -= 1;

    for(i=1; i<=n; i++) a[i] = &aa[i-1][0] - 1;
    d = dd - 1;
    e = ee - 1;

    for(i=n; i>=2; i--){

        l = i - 1;
        h = scale = 0.0;
        if(l > 1) {
           for(k=1; k <= l; k++) scale += fabs(a[i][k]);

           if(scale <= SCTOL) e[i] = a[i][l];

           else {

              for(k=1; k<=l; k++){


                  a[i][k] /= scale;
                  h += a[i][k] * a[i][k];

              }

              f = a[i][l];
              g = (f > 0) ? -sqrt(h) : sqrt(h);
              e[i] = scale * g;
              h -= f*g;
              a[i][l] = f - g;
              f = 0.0;
              for(j=1; j<=l; j++){
                  a[j][i] = a[i][j] / h;
                  g = 0.0;
                  for(k=1; k<=j; k++) g += a[j][k] * a[i][k];
                  for(k=j+1; k<=l; k++) g += a[k][j] * a[i][k];
                  e[j] = g / h;
                  f += e[j] * a[i][j];
              }

              hh = f / (h + h);
              for(j=1; j<=l; j++){

                  f =a[i][j];
                  e[j] = g = e[j] - hh * f;
                  for(k=1; k<=j; k++) a[j][k] -= (f * e[k] + g * a[i][k]);

              }

           }

        }
        else  e[i] = a[i][l];

        d[i] = h;
    }

    d[1] = 0.0;
    e[1] = 0.0;

    for(i=1; i <= n; i++){

        l = i - 1;
        if(d[i]) {

          for(j=1; j<=l; j++){

              g = 0.0;
              for(k=1; k<=l; k++) g+= a[i][k] * a[k][j];
              for(k=1; k<=l; k++) a[k][j] -= g * a[k][i];

          }

        }

        d[i] = a[i][i];
        a[i][i] = 1.0;
        for(j=1; j<=l; j++) a[j][i] = a[i][j] = 0.0;

    }

    a += 1;

    free(a);

    return;

}
