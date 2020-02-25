#include <Stdio.h>
#include <stdlib.h>

#ifndef  REGULAR

void surfit()

{

   printf("***error***, surface fitting impossible unless correct libraries\r\n"
          "             are linked, see compilation options\n");
   exit(1);

}

#else

#include <stdarg.h>
#include "bisp.h"
#include "reg_dat.h"
#include "sphery_dat.h"
#include "mem_er.h"
#include "grid.h"

#define  IOPT  0

void smoopy_setup(int );
int smoopy_c(double * , double  , int );
void sphery_setup(int );
int sphery_c(double * , double  , int );

struct sp_dat *cokn=NULL;

extern struct rspline *rs;
extern struct sspline *ss;

extern struct savedat *sd;
extern struct savedat_sphy *ssd;

extern float *ap;
extern double *apd;
extern GRID *gr;

void surfit(double *s, int foro, int smty, ...)

{

    va_list smptr;

    int ier;
    int nmax;
    int jj;
    int i, j;

    double fpp;

    va_start(smptr, smty);

    cokn = va_arg(smptr, struct sp_dat *);

    if(smty){

          ss = va_arg(smptr, struct sspline *);
          ssd = va_arg(smptr, struct savedat_sphy *);

    }

    else {

          rs = va_arg(smptr, struct rspline *);
          sd = va_arg(smptr, struct savedat *);

    }

/* assign memory for spline knots and coefficients at first use */

    if(foro == 0){

       cokn->sm_type = smty;

/* choose the smoothing factor at first function call */

       printf("What smoothing factor do you wish to use with ");
       (smty) ? printf("%s\n", "SPHERY") : printf("%s\n", "SMOOPY");
       
       printf("    s = 0  interpolating B-spline surface \r\n"
              "    s > 0  least squares B-spline surface   \n");

       scanf("%lf", s);

/* perform setup assignments for spline fitting */

       if(cokn->sm_type) {

          sphery_setup(0); nmax = ss->nmax; cokn->ncof = ss->ncof;

       }
       else {

          smoopy_setup(0); nmax = rs->nmax; cokn->ncof = rs->ncof;

       }


       cokn->tx = (double *)calloc(nmax, sizeof(double));
       mem_er((cokn->tx == NULL) ? 0 : 1, nmax * sizeof(double));

       cokn->ty = (double *)calloc(nmax, sizeof(double));
       mem_er((cokn->ty == NULL) ? 0 : 1, nmax * sizeof(double));

       cokn->c = (double *)calloc(cokn->ncof, sizeof(double));
       mem_er((cokn->c == NULL) ? 0 : 1, cokn->ncof * sizeof(double));


    }

    else if(foro == -1){  /* free memory */

       free(cokn->tx);
       free(cokn->ty);
       free(cokn->c);

       if(cokn->sm_type) sphery_setup(1);
       else smoopy_setup(1);

       va_end(smptr); 

       return;

    }

    if(cokn->sm_type){

/* put data into correct order, may want to change this later */

       if(ap){
          for(j=0; j < ss->mu; j++){
              jj = ss->mu - j - 1;
              for(i=0; i < ss->mv; i++){
                  *(ss->r + jj * ss->mv + i) = *(ap + j * gr->ix + i);
              }
          } 
       }
       else if (apd){
          for(j=0; j < ss->mu; j++){
              jj = ss->mu - j - 1;
              for(i=0; i < ss->mv; i++){
                  *(ss->r + jj * ss->mv + i) = *(apd + j * gr->ix + i);
              }
          }
       }
       else{
          printf("***ERROR***, no data available for surface fitting in file %s.\n\n", __FILE__);
          exit(1);
       }

       ier = sphery_c(&fpp, *s, IOPT);
       cokn->ncof = ss->ncof;

    }

    else {

       if(ap){
          for(j=0; j < rs->mx; j++){
             for(i=0; i < rs->my; i++) {
                *(rs->z + j * rs->my + i) = *(ap + (rs->sy +i-1) * gr->ix + rs->sx + j - 1);
             }
          }
       }
       else if(apd){
          for(j=0; j < rs->mx; j++){
             for(i=0; i < rs->my; i++) {
                *(rs->z + j * rs->my + i) = *(apd + (rs->sy +i-1) * gr->ix + rs->sx + j - 1);
             }
          }
       }
       else{
          printf("***ERROR***, no data available for surface fitting in file %s.\n\n", __FILE__);
          exit(1);
       }

       ier = smoopy_c(&fpp, *s, IOPT);
       cokn->ncof = rs->ncof;


    }

    printf("%s ier = %d residual = %e\n", __FILE__, ier, fpp);

    va_end(smptr); 

    return;  

}

#endif
