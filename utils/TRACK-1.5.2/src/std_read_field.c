#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grid.h"
#include "pp.h"
#include "netcdf_info.h"

#define  INV      -1.0
#define  MAXCHAR   100

extern int tl, gof;
extern int form;
extern int nw_ln;
extern float offs, ilat;
extern float *abuf;
extern int shp;
extern int neqt, ieq;
extern int std_x, std_y;

extern long int frtime;

extern GRID *gr;

extern PP *pph;

/* function to read in each field of data and apply an offset if required */

int pp_read(FILE * , int );

int std_read_field(float *ap, float *ata, float scl, FILE *fdatin, int sign, int as, int of, int hemi, int tscl)

{

    int i, j, k;
    int pp, plus=0, ppt, pj=0;
    int ndim=0, dim;
    int iret=0;
    int nx, dd;
    int insp=0;
    int ii=0;

    char texin[MAXCHAR];
    char *ctim;

    float *at=NULL;
    float *gw=NULL;
    float yy;

    frtime = 0;

    switch(form){
       case 0:

         dim = std_x * std_y;
         fgets(texin, MAXCHAR, fdatin);

         if(feof(fdatin)){
            printf("****WARNING****, EOF has been detected.\n\n");
            iret = 1;
            break;
         }

         if(REPORT)printf("%s\n", texin);

         if((ctim = strstr(texin, "TIME"))){
            sscanf(ctim, "%*s %ld", &frtime);
         }

         if((ndim = fread(abuf, sizeof(float), dim, fdatin)) != dim){
            printf("***WARNING***, problem reading in field data, \r\n"
                   "               insufficient data, %d != %d\n\n", ndim, dim);
            iret = (!ndim) ? 1 : 0;
            break;
         }
         if(nw_ln == 'y') fscanf(fdatin, "%*c");

         break;
       case 1:

         dim = std_x * std_y;
         fgets(texin, MAXCHAR, fdatin);

         if(feof(fdatin)){
            printf("****WARNING****, EOF has been detected.\n\n");
            iret = 1;
            break;
         }

         if(REPORT)printf("%s\n", texin);

         if((ctim = strstr(texin, "TIME"))){
            sscanf(ctim, "%*s %ld", &frtime);
         }

         for(i=0; i< dim; i++) fscanf(fdatin, "%e", abuf+i);

         fgets(texin, MAXCHAR, fdatin);

         break;

       case 3:

         iret = pp_read(fdatin, 1);
         iret = (iret) ? 0 : 1;
         break;

       case 4:
         iret = netcdf_read_field((void *)fdatin);
         break;
    }

    if(iret) return iret;

    if(gr->gcen == 1){
    
       dd = std_x - gr->nneg;

       gw = gr->wrk + dd;

       for(j=0; j < std_y; j++){

          pj = j * std_x;
          at = abuf + pj;

/*          for(i=gr->nneg; i < std_x; i++) *(gr->wrk + i - gr->nneg) = *(abuf + pj + i);
          for(i=0; i< gr->nneg; i++) *(gr->wrk + gr->nneg + i) = *(abuf + pj + i);*/

          memcpy(gr->wrk, at + gr->nneg, dd * sizeof(float));
          memcpy(gw, at, gr->nneg * sizeof(float));
          memcpy(at, gr->wrk, std_x * sizeof(float));

       }
    }

    if(gr->h_inv){
       for(i=0; i < std_y / 2; i++){
           memcpy(gr->wrk, abuf + i * std_x, std_x * sizeof(float));
           memcpy(abuf + i * std_x, abuf + (std_y - i - 1) * std_x, std_x * sizeof(float));
           memcpy(abuf + (std_y - i - 1) * std_x, gr->wrk, std_x * sizeof(float));
       }

    }

    nx = (gr->igc) ? gr->ix - 1: gr->ix;

    insp = (shp == 'n') ? std_x : 0;
    ii = (shp == 'n') ? 1 : 0;

    for(j=0; j < gr->iy; j++){

       pp = j * gr->ix;
       ppt = j * std_x + insp;
       if(neqt == 'n' && ii >= ieq) ppt += std_x;
       yy = *(gr->ygrid + j);

       for(i=0; i< nx; i++) {

          if(tl == 'n') plus = pp+i;

          else if(gof >= 0){

             k = i + gof;
             if(k >= gr->ix) k -= nx;
             plus = pp + k;

          }

          else{

             k = i + gof;
             if(k < 0) k += nx;
             plus = pp + k;

          }

          at = ap+plus;
          *at = *(abuf + ppt + i);    

          if(tscl == 'y') *at *= scl;

          if(as == 'y') *at = *at - *(ata+plus); 

          if(of == 'y') *at -= offs;

          *at *= sign;

          if(hemi == 'n' && yy > ilat) *at *= INV;
          else if(hemi == 's' && yy < ilat) *at *= INV;

       }


       if(gr->igc){
         if(gof > 0) *(ap + pp) = *(ap + pp + gr->ix - 1);
         else if (gof <= 0) *(ap + pp + gr->ix - 1) = *(ap + pp);

       }

       ++ii;

    }

    return iret;

}
