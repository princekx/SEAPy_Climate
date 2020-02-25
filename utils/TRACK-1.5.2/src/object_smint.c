#include <Stdio.h>
#include <stdlib.h>

#ifndef  REGULAR

void object_smint()

{

   printf("***error***, surface fitting and optimization impossible unless correct \r\n" 
          "             libraries are linked, see compilation options\n");
   exit(1);

}

#else

#include "st_obj.h"
#include "st_fo.h"
#include "bisp.h"
#include "reg_dat.h"
#include "mem_er.h"
#include "grid.h"



void error_interp(int , char * , int );
int smoopy_c(double * , double  , int );

#ifdef  NOUNDERSCORE

void bisp(double * , double * , double * , int * , double * , int * , 
           double * , int * , double * , double * , int * );

#else

void bisp_(double * , double * , double * , int * , double * , int * , 
           double * , int * , double * , double * , int * );

#endif

void object_local_maxs(struct frame_objs * , struct frame_objs * , int ,int );
void non_lin_opt(struct frame_objs * , int , int );

/* function to interpolate/smooth over a frame object */


extern struct savedat *sd;
extern struct rspline *rs;
extern struct sp_dat *cokn;

extern float *ap;

extern GRID *gr;
extern int x1u, x2u, y1u, y2u;
extern int delb, bs, pb;

extern float period;

void object_smint(struct frame_objs *fo, double *s, struct savedat *sdo, struct rspline *rso, struct sp_dat *cokno, int fstlst)

{

   int i, j, k;
   int nmaxx, nmaxy;
   int ix, iy;
   int err=0;
   int ier=0;
   int pxmx=0, pxmn=0, pymx=0, pymn=0;
   int xd1=x1u-delb, xd2=x2u+delb, yd1=y1u-delb, yd2=y2u+delb;
   int dd;
   int idif=0;
   int fn_tmp=0;
   int iter_sm=0;

   static int os='a';
   static int flin=1;
   static int fetty=0;
   static int dc;
   static int b_exc=0;

   double fpp;
   double da[4];

   double xx, yy, zz;
   double stemp;

   struct object *ob=NULL;
   struct point *ptt=NULL;

   struct frame_objs fp, ft;


/* if first entry assign spline degree data to structures */

   sd = sdo;
   rs = rso;
   cokn = cokno;


   if(fstlst == 0){

      rs->kx = KX;
      rs->ky = KY;
      rs->kx1 = KX + 1;
      rs->ky1 = KY + 1;
      rs->kx2 = KX + 2;
      rs->ky2 = KY + 2;
      rs->kmax = (rs->kx > rs->ky) ? rs->kx : rs->ky;
      rs->kmax2 = rs->kmax + 2;
      rs->k2max2 = 2 * rs->kmax + 2;

      printf("What is the value of the smoothing parameter ( >= 0)?\n");
      scanf("%lf", s);

      printf("Do you want to use all points within the object region for \r\n"
             "surface fitting or just the object points, 'a' or 'o'.     \n\n");

      scanf("\n");
      os = getchar();

      printf("Do you want grid maxima or off grid maxima, '0' or '1'?\n");
      scanf("%d", &fetty);

      printf("When finding grid local maxima do you want to exclude  \r\n"
             "boundary maxima, e.g.. maxima adjacent to threshold    \r\n"
             "boundary ormaxima adjacent to a land boundary for ocean\r\n"
             "problems?                                              \r\n"
             "Input '1' for exlusion, '0' otherwize.                 \n\n"); 

      scanf("%d", &b_exc);


      flin = 0;

   }

   else if (fstlst < 0){

      free(sd->nrdatx);
      free(sd->nrdaty);
      free(rs->x); 
      free(rs->y); 
      free(cokn->tx); 
      free(cokn->ty); 
      free(cokn->c);
      free(rs->z); 

      if(fetty) non_lin_opt(NULL, -1, 0);

      return;

   }

   else flin = 1;

   fo->tot_f_f_num = 0;

   if(fstlst == 0){

      printf("Do you want to perform a data check for each object, 'y' or 'n'\n");

      scanf("\n");
      dc = getchar();

   }

/* loop over all objects */

   for(j=0; j< fo->obj_num; j++){

      ob = fo->objs + j;

      if(!(ob->ext)){

        pxmx = pxmn = ob->pt->x;
        pymx = pymn = ob->pt->y;

        for(i=0; i < ob->point_num; i++){

            ptt = (ob->pt) + i;

            if(ptt->x > pxmx) pxmx = ptt->x;
            else if(ptt->x < pxmn) pxmn = ptt->x;

            if(ptt->y > pymx) pymx = ptt->y;
            else if(ptt->y < pymn) pymn = ptt->y;

        }

        ob->ext = (struct extent * )malloc_initl(sizeof(struct extent));
        mem_er((ob->ext == NULL) ? 0 : 1, sizeof(struct extent));

        ob->ext->x1 = pxmn;
        ob->ext->x2 = pxmx;
        ob->ext->y1 = pymn;
        ob->ext->y2 = pymx; 
 
         
      }

      rs->sx = x1u + ob->ext->x1 - 2;
      rs->sy = y1u + ob->ext->y1 - 2;
      rs->ex = x1u + ob->ext->x2;
      rs->ey = y1u + ob->ext->y2;

      rs->mx = rs->ex - rs->sx + 1;
      rs->my = rs->ey - rs->sy + 1;


/* Check region compatability and re-assign. */

      while(rs->mx < rs->kx1){
            --(rs->sx);
            ++(rs->ex);
            rs->mx += 2;
      }

      while(rs->my < rs->ky1){
            --(rs->sy);
            ++(rs->ey);
            rs->my += 2;
      }

      if(bs == 'y'){

         if(pb == 'y') dd = xd1;
         else dd = (xd1 > 1) ? xd1 : 1;

         if(rs->sx < dd && pb != 'y') { rs->ex += dd - rs->sx; rs->sx = dd;}

         if(pb == 'y') dd = xd2;
         else dd = (xd2 < gr->ix) ? xd2 : gr->ix;

         if(rs->ex > dd && pb != 'y') { rs->sx -= rs->ex - dd; rs->ex = dd;}

      }

      else {

         dd = x1u;
         if(rs->sx < dd){ rs->ex += dd - rs->sx; rs->sx = dd;}
         dd = x2u;
         if(rs->ex > dd){ rs->sx -= rs->ex - dd; rs->ex = dd;}

      }

      dd = (yd1 > 1) ? yd1 : 1;

      if(rs->sy < dd) { rs->ey += dd - rs->sy; rs->sy = dd;}

      dd = (yd2 < gr->iy) ? yd2 : gr->iy;

      if(rs->ey > dd) { rs->sy -= rs->ey - dd; rs->ey = dd;}


      rs->mxy = rs->mx * rs->my;
      rs->ncof = rs->mxy;

      cokn->ncof = rs->ncof;

      nmaxx = rs->nxest = rs->mx + rs->kx1;
      nmaxy = rs->nyest = rs->my + rs->ky1;
      rs->nmax = (nmaxx > nmaxy) ? nmaxx : nmaxy;
      ix = rs->sx - 1;
      if(ix < 0)
         rs->xb = (double) (* (gr->xgrid + gr->ix - 1 + ix) - period);
      else rs->xb = (double) * (gr->xgrid + ix);

      ix = rs->ex - 1;

      if(ix >= gr->ix)
        rs->xe = (double) (*(gr->xgrid + ix - gr->ix + 1) + period);
      else rs->xe = (double) * (gr->xgrid + ix);


      rs->yb = (double) * (gr->ygrid + rs->sy - 1);
      rs->ye = (double) * (gr->ygrid + rs->ey - 1);

/* memory allocation/re-allocation */

      if(!(sd->nrdatx)){

         sd->nrdatx = (int *)calloc(rs->nmax, sizeof(int));
         mem_er((sd->nrdatx == NULL) ? 0 : 1, rs->nmax * sizeof(int));

      }

      else {

         sd->nrdatx = (int *)realloc_n(sd->nrdatx, rs->nmax*sizeof(int));
         mem_er((sd->nrdatx == NULL) ? 0 : 1, rs->nmax*sizeof(int));

      }

      if(!(sd->nrdaty)){

         sd->nrdaty = (int *)calloc(rs->nmax, sizeof(int));
         mem_er((sd->nrdaty == NULL) ? 0 : 1, rs->nmax * sizeof(int));

      }

      else {

         sd->nrdaty = (int *)realloc_n(sd->nrdaty, rs->nmax*sizeof(int));
         mem_er((sd->nrdaty == NULL) ? 0 : 1, rs->nmax*sizeof(int));


      }

      if(!(rs->x)){

         rs->x = (double *)calloc(rs->mx, sizeof(double));
         mem_er((rs->x == NULL) ? 0 : 1, rs->mx * sizeof(double));

      }

      else {

         rs->x = (double *)realloc_n(rs->x, rs->mx*sizeof(double));
         mem_er((rs->x == NULL) ? 0 : 1, rs->mx*sizeof(double));

      }

      if(!(rs->y)){

         rs->y = (double *)calloc(rs->my, sizeof(double));
         mem_er((rs->y == NULL) ? 0 : 1, rs->my * sizeof(double));

      }

      else {

         rs->y = (double *)realloc_n(rs->y, rs->my*sizeof(double));
         mem_er((rs->y == NULL) ? 0 : 1, rs->my*sizeof(double));

      }

      if(!(cokn->tx)){

         cokn->tx = (double *)calloc(rs->nmax, sizeof(double));
         mem_er((cokn->tx == NULL) ? 0 : 1, rs->nmax * sizeof(double));

      }

      else {

         cokn->tx = (double *)realloc_n(cokn->tx, rs->nmax*sizeof(double));
         mem_er((cokn->tx == NULL) ? 0 : 1, rs->nmax*sizeof(double));

      }


      if(!(cokn->ty)){

         cokn->ty = (double *)calloc(rs->nmax, sizeof(double));
         mem_er((cokn->ty == NULL) ? 0 : 1, rs->nmax * sizeof(double));

      }

      else {

         cokn->ty = (double *)realloc_n(cokn->ty, rs->nmax*sizeof(double));
         mem_er((cokn->ty == NULL) ? 0 : 1, rs->nmax*sizeof(double));

      }

      if(!(cokn->c)){

         cokn->c = (double *)calloc(cokn->ncof, sizeof(double));
         mem_er((cokn->c == NULL) ? 0 : 1, cokn->ncof * sizeof(double));

      }

      else {

         cokn->c = (double *)realloc_n(cokn->c, cokn->ncof*sizeof(double));
         mem_er((cokn->c == NULL) ? 0 : 1, cokn->ncof*sizeof(double));

      }


      if(!(rs->z)){

         rs->z = (double *)calloc(rs->mxy, sizeof(double));
         mem_er((rs->z == NULL) ? 0 : 1, rs->mxy * sizeof(double));

      }

      else { 

         rs->z = (double *)realloc_n(rs->z, rs->mxy*sizeof(double));
         mem_er((rs->z == NULL) ? 0 : 1, rs->mxy*sizeof(double));
      }

      if(dc == 'y'){

         if(  (rs->kx <= 0 || rs->ky <= 0)                                               ||
              (rs->mx < rs->kx1 || rs->nxest < 2 * rs->kx1)                              ||
              (rs->my < rs->ky1 || rs->nyest < 2 * rs->ky1)                              ||
              (rs->xb > *(gr->xgrid + rs->sx - 1) || rs->xe < *(gr->xgrid + rs->ex - 1)) ||
              (rs->yb > *(gr->ygrid + rs->sy - 1) || rs->ye < *(gr->ygrid + rs->ey-1))      ){

              printf("***error***, file = %s, line = %d\n", __FILE__, __LINE__);
              printf("kx = %d <= 0 or ky = %d <= 0\n", rs->kx, rs->ky);
              printf("mx = %d < kx1 = %d or my = %d < ky1 = %d\n", rs->mx, rs->kx1, rs->my, rs->ky1);
              printf("defined region incompatable with data\n\n");
              error_interp(10, __FILE__, __LINE__);



         }


      }

/* assign grid arrays */

      for(i=0; i < rs->mx; i++){

         ix = rs->sx + i - 1;

         if(ix < 0)
           *(rs->x + i) = (double) (* (gr->xgrid + gr->ix - 1 + ix) - period);
         else if(ix >= gr->ix)
           *(rs->x + i) = (double) (* (gr->xgrid + ix - gr->ix + 1) + period);
         else
           *(rs->x + i) = (double) * (gr->xgrid + ix);

      }

      for(i=0; i < rs->my; i++) *(rs->y + i) = (double) * (gr->ygrid + rs->sy + i - 1);

      if(dc == 'y'){

         err = 0;

         for(i=1; i < rs->mx; i++)
             if(*(rs->x + i - 1) >= *(rs->x + i)) { err = 1; break;}

         if(err == 1 ) {

            printf("***error***, file = %s, line = %d\n", __FILE__, __LINE__);
            printf("***error*** in x-data %f %f\n\n", *(rs->x + i - 1), *(rs->x + i));
            error_interp(10, __FILE__, __LINE__);

        }

        err = 0;

        for(i=1; i < rs->my; i++)

           if(*(rs->y + i - 1) >= *(rs->y + i)) {err = 1; break;}


        if(err == 1 ){

          printf("***error***, file = %s, line = %d\n", __FILE__, __LINE__);
          printf("***error*** in y-data %f %f\n\n", *(rs->y + i - 1), *(rs->y + i));
          error_interp(10, __FILE__, __LINE__); 

        }  




      } 


/* assign field points */

      if(os == 'o'){

/* initialize field array */

         for(i=0; i < rs->mxy; i++) *(rs->z + i) = 0.0;

         for(i=0; i < ob->point_num; i++){

             ptt = (ob->pt) + i;

             ix = x1u + ptt->x - 2;


             if(ix < 0) 
               *(rs->z + (x1u + ptt->x - rs->sx - 1) * rs->my + (y1u + ptt->y - rs->sy - 1)) = *(ap + (y1u + ptt->y - 2) * gr->ix + gr->ix - 1 + ix);
             else if(ix >= gr->ix)
               *(rs->z + (x1u + ptt->x - rs->sx - 1) * rs->my + (y1u + ptt->y - rs->sy - 1)) = *(ap + (y1u + ptt->y - 2) * gr->ix + ix - gr->ix + 1);
             else 
               *(rs->z + (x1u + ptt->x - rs->sx - 1) * rs->my + (y1u + ptt->y - rs->sy - 1)) = *(ap + (y1u + ptt->y - 2) * gr->ix + x1u + ptt->x -2);

         }


      }

      else {

         for(i=0; i < rs->my; i++){

             iy = rs->sy + i - 1;

             for(k=0; k < rs->mx; k++){

                ix = rs->sx + k - 1;

                if(ix < 0)              
                   *(rs->z + k * rs->my + i) = *(ap + iy * gr->ix + gr->ix - 1 + ix);
                else if(ix >= gr->ix)
                   *(rs->z + k * rs->my + i) = *(ap + iy * gr->ix + ix - gr->ix + 1);
                else
                   *(rs->z + k * rs->my + i) = *(ap + iy * gr->ix + ix);

             }

         }

      }


/* surface fitting, if surface fitting error occurs, try more smoothing 
   upto NITER times.                                                     */

      ier = 100;
      iter_sm = 0;

      stemp = *s;

      while(ier > 0 && iter_sm < NITER){

         ier = 0;

         ier = smoopy_c(&fpp, stemp, 0);

         stemp *= SM_FAC;
         ++iter_sm;

         if(iter_sm >= NITER) break;

         if(ier > 0)
            printf("****WARNING****, a problem was encountered with surface fitting\r\n"
                   "                 for the smoothing %f, trying again with %f\n\n", stemp/SM_FAC, stemp);

      }
         


      cokn->ncof = rs->ncof;

      if(dc == 'y') printf("Object %d, ier = %d, residual = %e\n", j+1, ier, fpp);

/* re-assign object point values */


      if(ier == 0 || ier == -2){

         for(i=0; i< ob->point_num; i++){

             ptt = ob->pt + i;

             ix = x1u + ptt->x - 2;

             if(ix < 0)
                xx = (double) (* (gr->xgrid + gr->ix - 1 + ix) - period);
             else if(ix >= gr->ix)
                xx = (double) (* (gr->xgrid + ix - gr->ix + 1) + period);
             else
                xx = (double) * (gr->xgrid + ix);

             yy = (double) * (gr->ygrid + y1u + ptt->y - 2);

#ifdef  NOUNDERSCORE

             bisp(&zz, da, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof, &xx, &yy, &idif);

#else

             bisp_(&zz, da, cokn->tx, &cokn->nx, cokn->ty, &cokn->ny, cokn->c, &cokn->ncof, &xx, &yy, &idif);

#endif

             ptt->val = zz;
           

         }


      }


/* find grid point maxima for each object */

      fp.obj_num = 0;
      fp.objs = NULL;
      fp.tot_f_f_num = 0;

      ft.obj_num = 1;
      ft.objs = ob;
      ft.tot_f_f_num = 0;

      object_local_maxs(&ft, &fp, 8, b_exc);

      fn_tmp = ft.tot_f_f_num;

      if(fetty){

         non_lin_opt(&ft, flin, 0);

         if(!flin) flin = 1;

         fn_tmp = ft.tot_f_f_num;

      }

      fo->tot_f_f_num += fn_tmp;


   }

   return;

}

#endif
