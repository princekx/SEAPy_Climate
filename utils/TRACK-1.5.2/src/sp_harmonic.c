#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "grid.h"

#define TOLG  1.0e-4

extern float period;


float sp_harmonic(double *cof, float *newf, GRID *gr, int ntr, int nx, float *ap, float *filt)

{

     int i, j, n, m;
     int ntr1 = ntr+1;
     int np;
     int iper=0;

     float rms=0.;
     float diff;

     double *cc=NULL;
     double ff, frst=0.0;
     double *leng=NULL, *cleng=NULL;
     double prfl;
     complex *lng=NULL;

     float *ftmp=NULL;

     ftmp = newf;


     if(fabs(*(gr->xgrid) + period - *(gr->xgrid + nx - 1)) < TOLG){--nx; iper=1;}


     if(!cof){

        printf("****ERROR****, no spectral coefficients for calculation in %s\n\n", __FILE__);

        exit(1);


     }

     for(i=0; i<gr->iy; i++){

        leng = gr->aleng + i * gr->nleng;

        for(j=0; j<nx; j++){

            cc = cof;

            ff = 0.0;

            lng = gr->sico + j * ntr1;

            np = 0;


            for(m=0; m<ntr1;m++){


                cleng = leng + np;

                for(n=m; n <ntr1; n++){

                   prfl = *cleng * *(filt + n);

                   ff += *(cc++) * prfl * lng->real; 

                   if(m > 0)ff += *(cc++) * prfl * lng->imag;
                   ++cleng;

                }

                ++lng;

                np += (ntr1 - m);

            }

            if(!newf){

               diff = *(ap+i * gr->ix + j) - ff;
               rms += diff * diff;

            }

            else *(ftmp++) = ff; 


            if(j == 0) frst = ff;

        }

        if(iper) *(ftmp++) = frst;

     }


     return sqrt(rms / (float)(gr->iy * nx));

}
