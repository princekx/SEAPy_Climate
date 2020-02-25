#include <Stdio.h>
#include <Math.h>
#include <stdlib.h>
#include "mem_er.h"

#define  MAXITER  30
#define  LTOL  1.0e-8

#define SIGN(a,b)  (((b) < 0 ) ? -fabs(a) : fabs(a))

/* function to determine the eigen-values and vectors of a reduced matrix
   from householder.c                                                      */

void ql(double *da, double *ee, int n, double **zz)

{
   int m, l, iter, i, k;
   double s, r, p, g, f, dd, c, b;
   double **z, *d, *e;

   z = (double **)calloc(n, sizeof(double *));
   mem_er((z == NULL) ? 0 : 1, n * sizeof(double *));
   z -= 1;

   for(i=1; i<=n; i++) z[i] = &zz[i-1][0] - 1;
   d = da - 1;
   e = ee - 1;


   for(i=2; i<=n; i++) e[i-1] = e[i];
   e[n] = 0.;
   for(l=1; l<=n; l++){
       iter = 0;
       do {
           for(m=l; m <= n-1; m++){
               dd = fabs(d[m]) + fabs(d[m+1]);
               if((fabs(e[m]) + dd) - dd < LTOL) break;
           }

           if(m != l){
              if(iter++ == MAXITER) {
                 printf("Number of iterations exceeded in QL routine in %s at line %d\n", __FILE__, __LINE__);
                 exit(1);
              }

              g = 0.5 * (d[l+1] - d[l]) / e[l];
              r = sqrt((g*g) + 1.);
              g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
              s = c = 1.0;
              p = 0.0;
              for(i=m-1; i>=l; i--){
                  f = s * e[i];
                  b = c * e[i];
                  if(fabs(f) >= fabs(g)){
                     c = g / f;
                     r = sqrt((c*c) + 1.0);
                     e[i+1] = f * r;
                     c *= (s = 1. / r);
                  }

                  else {
                     s = f / g;
                     r = sqrt((s * s) + 1.0);
                     e[i+1] = g * r;
                     s *= (c = 1. / r);
                  }

                  g = d[i+1] - p;
                  r = (d[i] - g) * s + 2. * c * b;
                  p = s * r;
                  d[i+1] = g + p;
                  g = c * r - b;

                  for(k=1; k<=n; k++){

                      f = z[k][i+1];
                      z[k][i+1] = s * z[k][i] + c * f;
                      z[k][i] = c * z[k][i] - s * f;

                  }

              }

              d[l] -= p;
              e[l] = g;
              e[m] = 0.0;

           }

       } while(m != l);


   }

   z += 1;

   free(z);

   return;

}

/* sort eigen-values/vectors in decending order */

void eigsort(double *dd, double **vv, int n)

{

    int i, j, k;
    double p;
    double **v, *d;

    v = (double **)calloc(n, sizeof(double *));
    mem_er((v == NULL) ? 0 : 1, n * sizeof(double *));
    v -= 1;
    for(i=1; i<=n; i++) v[i] = &vv[i-1][0] - 1;
    d = dd - 1;

    for(i=1; i<n; i++){

       p = d[k=i];
       for(j=i+1; j <=n; j++){
           if(d[j] >= p) p = d[k=j];
       }

       if(k != i){
          d[k] = d[i];
          d[i] = p;
          for(j=1; j<=n; j++){
              p = v[j][i];
              v[j][i] = v[j][k];
              v[j][k] = p;
          }

        }

    }

    v += 1;

    free(v);

    return;

}
