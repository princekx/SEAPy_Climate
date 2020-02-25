#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#include "mem_er.h"

/* simple function to fill holes */

int missing(float, float, int );

void fill_holes(float *abuf, float mval, float mdval, int ix, int iy, int ifr, int icmp)

{

    int i, j;
    int ii, i0;
    int nav=0;
    int ixx=ix+2, iyy=iy+2;

    static int dim;
    static int isd=0;

    float *aa=NULL, *at=NULL, *att=NULL;
    float sum=0.0;

    static float *atmp=NULL;

    if(!ifr){

       dim = ixx * iyy;
       atmp = (float *)calloc(dim, sizeof(float));
       mem_er((atmp == NULL) ? 0 : 1, dim * sizeof(float));

       printf("Do you want a foreward raster scan, '0' or both foreward and backward, '1' \n\n");
       scanf("%d", &isd);

    }

    else if (ifr < 0){free(atmp); return;}

    for(i=0; i<dim; i++) *(atmp + i) = mdval;

    for(i=0; i < iy; i++)
       memcpy(atmp + (i+1)*ixx + 1, abuf+i*ix, ix*sizeof(float));

    for(i=0; i < iy; i++){

       for(j=0; j < ix; j++){

           ii = i * ix + j;
           i0 = (i+1) * ixx + j + 1;
           aa = abuf + ii;
           at = atmp + i0;


           if(missing(*aa, mval, icmp)){
              nav = 0;
              sum = 0.0;
              att = at - 1;
              if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
              att = at + 1;
              if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
              att = at - ixx - 1;
              if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
              att = at - ixx;
              if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
              att = at - ixx + 1;
              if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
              att = at + ixx - 1;
              if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
              att = at + ixx;
              if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
              att = at + ixx + 1;
              if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
              if(nav > 0)*at = *aa = sum / nav;
           } 


       }

    }

    if(isd){

       for(i=iy-1; i >= 0; i--){

          for(j=ix-1; j >= 0; j--){

              ii = i * ix + j;
              i0 = (i+1) * ixx + j + 1;
              aa = abuf + ii;
              at = atmp + i0;


              if(missing(*aa, mval, icmp)){

                 nav = 0;
                 sum = 0.0;
                 att = at - 1;
                 if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
                 att = at + 1;
                 if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
                 att = at - ixx - 1;
                 if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
                 att = at - ixx;
                 if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
                 att = at - ixx + 1;
                 if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
                 att = at + ixx - 1;
                 if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
                 att = at + ixx;
                 if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
                 att = at + ixx + 1;
                 if(!missing(*att, mval, icmp)) {sum += *att; ++nav;}
                 if(nav > 0)*at = *aa = sum / nav;

              } 


          }

       }


    }

    return;

}
