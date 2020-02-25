#include <Stdio.h>
#include <Math.h>
#include <stdlib.h>
#include "sphery_dat.h"
#include "mem_er.h"
#include "m_values.h"
#include "grid.h"


/* function to perform the setup for interpolation on a spherical surface.
   converted from the sphery.f code of DIERCKX. Hard wired to bi-cubics.   
   Data is assumed to be in lat-long (degrees) and is converted to radians
   before use.                                                             */


void  error_interp(int , char * , int , int );

struct savedat_sphy *ssd;
struct sspline *ss;

extern float xmn, ymn, xmx, ymx;
extern int x1u, x2u, y1u, y2u;

extern int delb, bs;

extern GRID *gr;
extern float period;

void sphery_setup(int mem)

{

    int i;
    int nmaxu, nmaxv;

    double *uu;


    if(mem){    /* free memory if assigned */

       free(ssd->nrdatu);
       free(ssd->nrdatv);
       free(ss->u);
       free(ss->v);
       free(ss->r);
       return;

    }


    printf("////////////////////////////////////////////////////////////////\r\n"
           "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");

    printf("***INFORMATION***, performing setup for surface fitting on a   \r\n"
           "                   sphere using bi-cubic splines. Continuity   \r\n"
           "                   properties at the poles are specified by    \r\n"
           "                   the user by setting options (See sphery.f). \n\n");

    printf("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\n");

    printf("Do you you want to specify continuity properties at the poles 'y' or 'n'\n\n");

    printf("***INFORMATION***, for best results with the local optimization\n"
           "                   scheme choose C1 continuity at the poles  \n\n");

    ss->iopt[2] = ss->iopt[1] = 0;

    scanf("\n");
    if(getchar() == 'y'){

       printf("What continuity is required at the poles, for C1 input '1'\r\n"
              "otherwise no continuity is assumed\n\n");

       scanf("%d", &ss->iopt[1]);
       ss->iopt[2] = ss->iopt[1];
       if(ss->iopt[1] > 1 || ss->iopt[1] < 0) {

          printf("***ERROR***, incorrect option for polar continuity,  \n"
                 "             performing surface fit on a sphere with \n"
                 "             no continuity at the poles.             \n\n");

          ss->iopt[2] = ss->iopt[1] = 0;

       }

    }

    ss->iopt[0] = 0;    /* always start as a completely new calculation
                          i.e. discard knots from previous calculations. */
    ss->mu = gr->iy;
    ss->mv = (fabs(*(gr->xgrid + gr->ix - 1) - period) < TOLGRID) ? gr->ix - 1 : gr->ix;
    ss->muv = ss->mv * gr->iy;

/*  information at the poles currentely set to none (see sphery.f) */

    ss->idr[0] = -1;  /* if == 0, 1 must set r0 */
    ss->idr[1] = 0;
    ss->idr[2] = -1;  /* if == 0, 1 must set r1 */
    ss->idr[3] = 0;

    ss->r0 = ss->r1 = 0.;

    nmaxu = ss->mu + 6 + ss->iopt[1] + ss->iopt[2];
    nmaxv = ss->mv + 7;

    ss->nuest = nmaxu;
    ss->nvest = nmaxv;

    ss->mvv = ss->mv + ss->nvest;
    ss->nmax = (nmaxu > nmaxv) ? nmaxu : nmaxv;
    ss->ncof = (ss->nuest - 4) * (ss->nvest - 4);

/* perform some checks on input so far */

    if(ss->nuest < 8 || ss->nvest < 11) error_interp(10, __FILE__, __LINE__, 1);
    if(ss->iopt[0] < 0 || ss->iopt[0] > 1) error_interp(10, __FILE__, __LINE__, 1);
    if(ss->iopt[1] < 0 || ss->iopt[1] > 1) error_interp(10, __FILE__, __LINE__, 1);
    if(ss->iopt[2] < 0 || ss->iopt[2] > 1) error_interp(10, __FILE__, __LINE__, 1);
    if(ss->idr[0] < -1 || ss->idr[0] > 1) error_interp(10, __FILE__, __LINE__, 1);
    if(ss->idr[1] < 0 || ss->idr[1] > 1) error_interp(10, __FILE__, __LINE__, 1);
    if(ss->idr[2] < -1 || ss->idr[2] > 1) error_interp(10, __FILE__, __LINE__, 1);
    if(ss->idr[3] < 0 || ss->idr[3] > 1) error_interp(10, __FILE__, __LINE__, 1);

    ss->mumin = 4;

    if(ss->idr[0] >= 0 ) ss->mumin -= 1;
    if(ss->iopt[1] == 1 && ss->idr[1] == 1) ss->mumin -= 1;
    if(ss->idr[2] >= 0) ss->mumin -= 1;
    if(ss->iopt[2] == 1 && ss->idr[3] == 1) ss->mumin -= 1;
    if(!(ss->mumin)) ss->mumin = 1;
    if(ss->mu < ss->mumin || ss->mv < 4) error_interp(10, __FILE__, __LINE__, 1);


/* assign storage to pointers */

    ssd->nrdatu = (int * )calloc(ss->nmax, sizeof(int));
    mem_er((ssd->nrdatu == NULL) ? 0 : 1, ss->nmax * sizeof(int));

    ssd->nrdatv = (int * )calloc(ss->nmax, sizeof(int));
    mem_er((ssd->nrdatv == NULL) ? 0 : 1, ss->nmax * sizeof(int));

/* assign memory for grid data */


    ss->u = (double * )calloc(ss->mu, sizeof(double));
    mem_er((ss->u == NULL) ? 0 : 1, ss->mu * sizeof(double));

    ss->v = (double * )calloc(ss->mv, sizeof(double));
    mem_er((ss->v == NULL) ? 0 : 1, ss->mv * sizeof(double));

    ss->r = (double *)calloc(ss->muv, sizeof(double));
    mem_er((ss->r == NULL) ? 0 : 1, ss->muv * sizeof(double));

   ss->vb = 0.;
   ss->ve = FPI2;
   ss->ub = 0.;
   ss->ue = FPI;

/* initialize grid arrays */

   for(i=0; i < ss->mu; i++){

      uu = ss->u + ss->mu - i - 1;

      *uu = FP_PI2 - (double) *(gr->ygrid + i) * FP_PI;
      if(*uu > FPI) *uu = FPI;
      else if (*uu < 0.) *uu = 0.;

   }

   for(i=0; i < ss->mv; i++)

      *(ss->v + i) = (double) *(gr->xgrid + i) * FP_PI;
   

/* perform some more checks on the data */

   if(*(ss->u) <= 0. || *(ss->u + ss->mu - 1) >= FPI)
 
       error_interp(10, __FILE__, __LINE__, 1);


   if(ss->mu > 1){

      for(i=1; i < ss->mu; i++) {

          if(*(ss->u + i - 1) >= *(ss->u + i)) 

             error_interp(10, __FILE__, __LINE__, 1);

      }


   }


   if(*(ss->v) != 0. || *(ss->v + ss->mv - 1) >= 2.*FPI)

      error_interp(10, __FILE__, __LINE__, 1);


   for(i=1; i < ss->mv; i++){

       if(*(ss->v + i - 1) >= *(ss->v + i))

         error_interp(10, __FILE__, __LINE__, 1);

   }





    return;

}
