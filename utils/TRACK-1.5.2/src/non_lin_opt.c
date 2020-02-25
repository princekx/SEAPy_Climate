#include <Stdio.h>
#include <stdlib.h>

#ifndef  OPT

void non_lin_opt()

{

    printf("***error***, non linear optimization impossible unless correct\r\n"
           "              libraries are linked, see compilation options\n");
    exit(1);

}

#else

#include <Math.h>
#include "mem_er.h"
#include "constraint.h"
#include "file_handle.h"
#include "reg_dat.h"
#include "sphery_dat.h"
#include "st_obj.h"
#include "st_fo.h"
#include "bisp.h"
#include "m_values.h"
#include "grid.h"

#define  XDIM     2      /* dimension of object function */
#define  IOPT     0      /* flag for equality constraint operator storage */
#define  TOLSEP   0.0001 /* tolerance on feature point seperation */
#define  TOLEQ    0.0001 /* need to counteract linear dependance of constraints */
#define  TOLPER   1.0e-6 /* tolerance on points on periodic boundary. */

float measure(struct feature_pts * , struct feature_pts * );
double gdfp_optimize(double * , double , int , int , int , int * );

struct cnst cst;

double *ih=NULL, *inn=NULL;
double *norm=NULL;

int ocon='d';
int ncon;

extern struct rspline *rs;
extern struct sspline *ss;
extern struct sp_dat *cokn;

extern int x1u, y1u;
extern int pb;
extern int tf;

extern GRID *gr;

void non_lin_opt(struct frame_objs *ff, int forl, int icons)

{

    int ndim, i, j, k, ii, xdim = XDIM;
    int iopt=IOPT;
    int ier;

    double *np, xy;
    double mvec[XDIM];

    float tstr;

    struct object *ob;
    struct feature_pts *ptt;

    FILE *fcon=NULL;

/* assign memory for constraints */

    if(forl == 0){

       if(tf == 8) ncon = 'y';
       else {

          printf("do you want constrained optimization, 'y' or 'n'\n");
          scanf("\n");
          ncon = getchar();

       }

       if(ncon == 'y'){

          if(tf == 8) ocon = 'o';

          else {

             printf("***IMPORTANT***, see the file 'data/constraints.dat' for restrictions\n");

             printf("do you want to use the default constraints \r\n"
                    "or the object extent constraints, 'd' or 'o'.   \n");

             scanf("\n");
             ocon = getchar();

          }

          if(cokn->sm_type && ocon == 'o'){

             printf("***WARNING***, can only use default constraints with the\n"
                    "               spherical surface fit, changing option to\n"
                    "               default, check constraints are correct.  \n"
                    "               Constraints must be at least for         \n"
                    "               co-latitude and in radians.          \n\n");

             ocon = 'd';

          }


          if(icons == 0){

             printf("***INFORMATION***, using constarint file %s\n\n", CONSDAT_SMOOPY);

             fcon = open_file(CONSDAT_SMOOPY, "r");

          }

          else if(icons == 1){
             printf("***INFORMATION***, using constarint file %s\n\n", CONSDAT_SPHERY);

             fcon = open_file(CONSDAT_SPHERY, "r");

          }

          fscanf(fcon, "%d %d", &cst.numc, &cst.numeq);

          if(ocon == 'o' && (cst.numeq + 2*xdim > cst.numc)){

             printf("***error***, in %s incompatable number of constraints\n", __FILE__);
             exit(1);

          }

          if(cst.numc){

             ndim = xdim * cst.numc;
             cst.nn = (double * )calloc(ndim, sizeof(double));
             mem_er((cst.nn == NULL) ? 0 : 1, ndim * sizeof(double));
 
             cst.bb = (double * )calloc(cst.numc, sizeof(double));
             mem_er((cst.bb == NULL) ? 0 : 1, cst.numc * sizeof(double));

          
             norm = (double * )calloc(cst.numc, sizeof(double));
             mem_er((norm == NULL) ? 0 : 1, cst.numc * sizeof(double));


             for(i=0; i < xdim; i++)
  
                 for(j=0; j < cst.numc; j++)

                     fscanf(fcon, "%lf", cst.nn+j*xdim+i);

             i = 0;

             while(fscanf(fcon, "%lf", cst.bb+(i++))){};


          }


          if(icons == 0){
            close_file(fcon, CONSDAT_SMOOPY);
          }
          else if(icons == 1){
            close_file(fcon, CONSDAT_SPHERY);
          }


/* normalize constraints */

          for(i=0; i < cst.numc; i++){

              norm[i] = 0.;
              ii = i * xdim;

              for(j=0; j < xdim; j++){np = cst.nn + ii + j; norm[i] += *np * (*np);}

              norm[i] = 1. / sqrt(norm[i]);

              for(j=0; j < xdim; j++)*(cst.nn + ii + j) *= norm[i];

              *(cst.bb + i) *= norm[i];

          }
 

          if(cst.numeq){

             ih = (double * )calloc(xdim*xdim, sizeof(double));
             mem_er((ih == NULL) ? 0 : 1, xdim*xdim * sizeof(double));

             inn = (double * )calloc(cst.numc*cst.numc, sizeof(double));
             mem_er((inn == NULL) ? 0 : 1, cst.numc*cst.numc * sizeof(double));

          }

       }

    }

    else if(forl == -1){

          if(ncon == 'y'){

             if(cst.numc) {free(cst.nn); free(cst.bb); free(norm);}

             if(cst.numeq) {free(ih); free(inn);}

          }

          return;

    }



    if(ncon == 'y' && ocon == 'd'){

        j = cst.numeq;
        np = cst.bb + j;

/* IMPORTANT that rs->xb, rs->xe, rs->yb and rs->ye match the correct 
   constraints read from constraints.dat file.                       */

        if(cokn->sm_type){

          *(np++) = *(cokn->ty + 3) * norm[j++];
          *(np++) = *(cokn->tx + 3) * norm[j++];
          *(np++) = - *(cokn->ty + cokn->ny - 4) * norm[j++];
          *(np++) = - *(cokn->tx + cokn->nx - 4) * norm[j++];

        }

        else {

          *(np++) = *(cokn->tx + 3) * norm[j++];
          *(np++) = *(cokn->ty + 3) * norm[j++];
          *(np++) = - *(cokn->tx + cokn->nx - 4) * norm[j++];
          *(np++) = - *(cokn->ty + cokn->ny - 4) * norm[j++];

        }

        j = cst.numeq;       

        if(fabs(*(cst.bb+j) + *(cst.bb+j+2)) < TOLEQ   || 
           fabs(*(cst.bb+j+1) + *(cst.bb+j+3)) < TOLEQ){

           printf("***error***, attempt to use linearly dependant positivity constraints, \n"
                 "             if this is what you really want, use an equality constraint\n");
           exit(1);

        }

    }


    for(i=0; i < ff->obj_num; i++){

        ob = (ff->objs) + i;

        if(ocon == 'o'){

           j = cst.numeq;
           np = cst.bb + j;

/* check consistency between constraints and extent of B-spline approx. 
   and assign constraints accordingly.                                    */

           if(tf == 8){

              *(np++) = rs->xb * norm[j++];
              *(np++) = rs->yb * norm[j++];
              *(np++) = - rs->xe * norm[j++];
              *np = - rs->ye * norm[j];

           }

           else {

              xy = obj_xreal(ob->ext->x1 + x1u - 2);
              if(xy < rs->xb) xy = rs->xb;
              *(np++) = xy * norm[j++];

              xy = *(gr->ygrid + ob->ext->y1 + y1u - 2);
              if(xy < rs->yb) xy = rs->yb;
              *(np++) = xy * norm[j++];

              xy = obj_xreal(ob->ext->x2 + x1u - 2);
              if(xy > rs->xe) xy = rs->xe;
              *(np++) = - xy * norm[j++];

              xy = *(gr->ygrid + ob->ext->y2 + y1u - 2);
              if(xy > rs->ye) xy = rs->ye;
              *np = - xy * norm[j];

           }



/* this next check will be corrected in the optimization routine if it violates 
   any  equality constraints. Check for linear dependance in positivity 
   constraints.                                                                 */

           j = cst.numeq;
           if(fabs(*(cst.bb+j) + *(cst.bb+j+2)) < TOLEQ){

              if(fabs(*(cst.bb+j)-rs->xb) < TOLEQ) *(cst.bb+j+2) += TOLEQ;
              else *(cst.bb+j) -= TOLEQ;

           }

           if(fabs(*(cst.bb+j+1) + *(cst.bb+j+3)) < TOLEQ){

              if(fabs(*(cst.bb+j+1)-rs->yb) < TOLEQ) *(cst.bb+j+3) += TOLEQ;
              else *(cst.bb+j+1) -= TOLEQ;

           }

        }

        for(j=0; j < ob->fet->feature_num; j++){

            ptt = (ob->fet->fpt) + j;

            tstr = ptt->str;

            if(cokn->sm_type){

               mvec[0] = FP_PI2 - (double)(ptt->y).xy * FP_PI;
               if(mvec[0] < 0.) mvec[0] = 0.;
               else if(mvec[0] > FPI) mvec[0] = FPI;

               mvec[1] = (double)(ptt->x).xy * FP_PI;

               ptt->str = (float)gdfp_optimize(mvec, (double)ptt->str, xdim, 0, iopt, &ier);

               if(ocon == 'd' && pb == 'y'){

                 if(mvec[1] - *(gr->xgrid) * FP_PI < TOLPER){

                    mvec[1] = *(gr->xgrid + gr->ix - 1) * FP_PI ;
                    ptt->str = (float)gdfp_optimize(mvec, (double)ptt->str, xdim, 0, iopt, &ier);

                 }

                 else if(*(gr->xgrid + gr->ix -1) * FP_PI - mvec[1] < TOLPER){

                    mvec[1] = *(gr->xgrid) * FP_PI;
                    ptt->str = (float)gdfp_optimize(mvec, (double)ptt->str, xdim, 0, iopt, &ier);

                 }


               }

               if(!ier){

                  (ptt->x).xy = (float)(mvec[1] / FP_PI);
                  (ptt->y).xy = (float)((FP_PI2 - mvec[0]) / FP_PI);
                  if((ptt->y).xy > 90. ) (ptt->y).xy = 90.;
                  else if((ptt->y).xy < -90.) (ptt->y).xy = -90.;

               }

               else{

                  printf("****WARNING****, a problem was encountered during the    \r\n"
                         "                 during optimization for feature points. \r\n"
                         "                 Defaulting to grid point feature points.\n\n");

                  ptt->str = tstr;

               }


            }

            else {

               mvec[0] = (double)(ptt->x).xy;
               mvec[1] = (double)(ptt->y).xy;


               ptt->str = (float)gdfp_optimize(mvec, (double)ptt->str, xdim, 0, iopt, &ier);

               if(ocon == 'd' && pb == 'y'){

                 if(mvec[0] - *(gr->xgrid) < TOLPER){

                    mvec[0] = *(gr->xgrid + gr->ix - 1);
                    ptt->str = (float)gdfp_optimize(mvec, (double)ptt->str, xdim, 0, iopt, &ier);

                 }

                 else if(*(gr->xgrid + gr->ix -1) - mvec[0] < TOLPER){

                    mvec[0] = *(gr->xgrid);
                    ptt->str = (float)gdfp_optimize(mvec, (double)ptt->str, xdim, 0, iopt, &ier);

                 }


               }


               if(!ier){

                  (ptt->x).xy = (float)mvec[0];
                  (ptt->y).xy = (float)mvec[1];

               }


               else{

                  printf("****WARNING****, a problem was encountered during the    \r\n"
                         "                 during optimization for feature points. \r\n"
                         "                 Defaulting to grid point feature points.\n\n");

                  ptt->str = tstr;

               }


            }


            for(k=0; k < j; k++)

               if(measure(ptt, (ob->fet->fpt) + k) < TOLSEP) ptt->str = DUFF_PT;


        }

    }

    return;

}

#endif
