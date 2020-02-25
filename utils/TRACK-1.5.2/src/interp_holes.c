#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#include "mem_er.h"
#include "grid.h"

#define  HW           4

/* simple function to fill holes by local interpolation */

extern GRID *gr;

int missing(float, float, int);


void interp_holes(float *abuf, float mval, int ix, int iy, int ifr, int icmp)

{

    int i, j, k, l;
    int i1, i2, j1, j2;
    int ii=0, jn=0;
    int nb0=0, nb1=0, nb2=0;
    int iy0=0, iy1=0, iy2=0;
    int np=0;
    int gtyp=0;
    int nn=0;

    static int dim=0;
    static int *mask=NULL;

    float *aa=NULL, *at=NULL;
    float f0=0.0, f1=0.0;
    float sum=0.0;

    static float *alat=NULL;


    if(gr->igtyp == 0 || gr->igtyp == 1) gtyp = 1;

    if(!ifr){

       dim = ix * iy;
       alat = (float *)calloc(dim, sizeof(float));
       mem_er((alat == NULL) ? 0 : 1, dim * sizeof(float));

       mask = (int *)calloc(dim, sizeof(int));
       mem_er((mask == NULL) ? 0 : 1, dim * sizeof(int));

    }

    else if(ifr < 0){ free(alat); free(mask); return;}

    for(i=0; i < dim; i++) {*(alat + i) = mval; *(mask + i) = 0;}

    for(i=0; i < iy; i++){

       for(j=0; j < ix; j++){

           nb0 = nb1 = 0;

           f0 = f1 = 0.0;

           ii = i * ix + j;

           aa = abuf + ii;
           if(missing(*aa, mval,icmp)){

              *(mask + ii) = 1;

/* find extent in latitudinal range for this grid point */
              iy0 = i - 1;

              while(iy0 >= 0){

                   at = abuf + iy0 * ix + j;
                   if(!missing(*at, mval, icmp)) {f0 = *at; nb0 = 1; break;}
                   --iy0;
              }

              if(iy0 < 0)iy0 = 0;
              
              iy1 = i + 1;

              while(iy1 < iy){

                   at = abuf + iy1 * ix + j;
                   if(!missing(*at, mval, icmp)) {f1 = *at; nb1 = 1; break;}
                   ++iy1;

              }

              if(iy1 == iy) iy1 = iy - 1;


              if(nb0 && nb1){

                 np = iy1 - iy0;

                 for(k=1; k < np; k++){

                     at = alat + (iy0 + k) * ix + j;
                     *at = ((float)(np - k) * f0 + (float)k * f1) / (float) np;

                 }

              }


/* SH polar interpolation */

              else if(!nb0 && nb1 && gtyp){

                  jn = j + ix / 2;
                  if(jn >= ix) jn -= ix;

                  iy2 = 1;

                  while(iy2 < iy){

                      at = abuf + iy2 * ix + jn;
                      if(!missing(*at, mval, icmp)) {f0 = *at; nb2 = 1; break;}
                      ++iy2;

                  }

                  if(iy2 == iy) iy2 = iy - 1;

                  np = iy1 + iy2;

                  if(nb1 && nb2){

                     for(k=1; k <= iy1; k++){

                        at = alat + (iy1 - k) * ix + j;
                        *at = ((float)k * f0 + (float)(np - k) * f1) / (float) np;
                     }

                     for(k=1; k < iy2; k++){

                        at = alat + k * ix + jn;
                        *at = ((float)(iy1 + k) * f0 + (float)(np - iy1 - k) * f1) / (float) np;
                     }

                  } 


              } 

           } 

       }

    }



    for(i=0; i < iy; i++){

       for(j=0; j < ix; j++){

           ii = i * ix + j;
           aa = abuf + ii;
           at = alat + ii;
           if(*(mask + ii)) *aa = *at;

        }


    }


/* perform local averaging */


    for(i=0; i < iy; i++){

       for(j=0; j < ix; j++){ 

           ii = i * ix + j;
           aa = abuf + ii;
           at = alat + ii;

           if(*(mask + ii)){

              i1 = i - HW;
              if(i1 < 0) i1 = 0;
              i2 = i + HW;
              if(i2 >= iy) i2 = iy - 1;

              j1 = j - HW;
              if(j1 < 0) j1 = 0;
              j2 = j + HW;
              if(j2 >= ix) j2 = ix - 1;

              sum = 0.0;
              nn = 0;

              for(k=i1; k <= i2; k++){

                  for(l=j1; l <= j2; l++){

                      sum += *(abuf + k * ix + l);
                      ++nn;

                  }

              }

              *at = sum / (float) nn;

           }

       }

    }   

/* Replace missing values */

    for(i=0; i < iy; i++){

       for(j=0; j < ix; j++){

           ii = i * ix + j;
           aa = abuf + ii;
           at = alat + ii;
           if(*(mask + ii)) *aa = *at;
           
       }

    } 

    return;

}


