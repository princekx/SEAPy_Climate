#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "st_fo.h"
#include "boundary.h"
#include "st_obj.h"
#include "mem_er.h"
#include "grid.h"
#include "m_values.h"

/* function to perform an 8-connected search for a boundary from a
   thresholded data set.                                           */


void dist_trans(int * , int , int , int );
void dir_search(int * , int * , int * , int * , int , int );
int powi(int , int );
void dfft(complex * , int , int , int );
float obj_xreal(int );
int inside(int , struct boundary_pt * , float , float );

extern int x1u, y1u;
extern int cc;

extern GRID *gr;
extern float period;
extern int tf;
extern int bs;
extern int pb;
extern int dfil;

extern float *ap;

int boundary_find(struct object *ob, struct boundary_cntl *bcntl, int *llb, int *rrb, int nfour)
{

    int nf=0, nexp=0;
    int ix1=0;
    int nmf=0;
    int ii;
    int nm;
    int ix=0, iy=0;

    int oset=0;
    int os = 1;

    int j, k;
    int pxmx, pxmn, pymx, pymn, ydim, xdim;
    int dim;
    int *aa=NULL, ab;
    int ixf=0, iyf=0;
    int ixx, iyy;
    int ifnd=0;
    int mxb=0;
    int dir=0;
    int nfet=0;
    int nly, nry;
    int oobptn=0;
    int pty;

    float *ob_field=NULL;
    float val=0.0;
    float fx;
    float xx, yy;

    struct point *ptt=NULL;
    struct boundary_pt *btt=NULL;
    struct feature_pts *fpt=NULL;

    complex *bbb=NULL;

    if(cc == 'e') mxb = 1;
    else if(cc == 'v') mxb = 2;
    else {
       printf("****ERROR****, connectivity tag %c not known,              \r\n"
              "               defaulting to 4-connected distance transform\r\n"
              "               for boundary search.\n\n", cc);
       cc = 'e';
       mxb = 1;

    }


    if(!(ob->point_num)) return nfet;
    if(ob->point_num == 1) return nfet;

    pxmx = pxmn = ob->pt->x;
    pymx = pymn = ob->pt->y;

    for(j=0; j < ob->point_num; j++){

        ptt = (ob->pt) + j;

        if(ptt->x > pxmx) pxmx = ptt->x;
        else if(ptt->x < pxmn) pxmn = ptt->x;

        if(ptt->y > pymx) pymx = ptt->y;
        else if(ptt->y < pymn) pymn = ptt->y;

    }

    xdim = pxmx - pxmn + 3;
    ydim = pymx - pymn + 3;

    dim = xdim*ydim;


    aa = (int *)calloc(dim, sizeof(int));
    mem_er((aa == NULL) ? 0 : 1, dim * sizeof(int));

    ob_field = (float *)calloc(dim, sizeof(float));
    mem_er((ob_field == NULL) ? 0 : 1, dim * sizeof(float));

    for(j=0; j < dim; j++) {*(aa + j) = 0; *(ob_field + j) = 0.0;}

    for(j=0; j < ob->point_num; j++){

       ptt = (ob->pt) + j;

       *(aa+(ptt->y - pymn + os)*xdim + (ptt->x - pxmn + os)) = 1;
       *(ob_field+(ptt->y - pymn + os)*xdim + (ptt->x - pxmn + os)) = ptt->val;

    }


    dist_trans(aa, xdim, ydim, 0);


    for(j=0; j<dim; j++) {

        if(*(aa+j) > mxb) *(aa+j) = 0;

    }



/* find the first boundary point */

    ob->bound = (struct boundary_pt *)malloc_initl(sizeof(struct boundary_pt));
    mem_er((ob->bound == NULL) ? 0 : 1, sizeof(struct boundary_pt));

    ifnd = 0;

    for(j=0; j < ydim; j++){

        for(k=0; k < xdim; k++){

            if(*(aa + j * xdim + k) == mxb){ 
              ixf = k; 
              iyf = j; 
              val = *(ob_field + j*xdim + k);
              ifnd = 1;
              break;
            }

        }

        if(ifnd) break;

    }

    (ob->bound->x).ixy = ixf;
    (ob->bound->y).ixy = iyf;
    ob->bound->val = val;
    ixx = ixf;
    iyy = iyf;
    ob->bound_num += 1;

    ixf = -1;
    iyf = -1;

    btt = ob->bound;

    dir = 1;

    while(!((ixx == ixf) && (iyy == iyf))){

          ixf = (btt->x).ixy;
          iyf = (btt->y).ixy;

          dir_search(aa, &ixf, &iyf, &dir, mxb, xdim);

          ++(ob->bound_num);

          ob->bound = (struct boundary_pt *)realloc_n(ob->bound, ob->bound_num * sizeof(struct boundary_pt));
          mem_er((ob->bound == NULL) ? 0 : 1, ob->bound_num * sizeof(struct boundary_pt));

          btt = ob->bound + ob->bound_num - 1;
          (btt->x).ixy = ixf;
          (btt->y).ixy = iyf;
          btt->val = *(ob_field + iyf*xdim + ixf);

          if(ixf == ixx && iyf == iyy) break;

    }

    for(j=0; j < ob->bound_num; j++){

        btt = ob->bound + j;

        (btt->x).ixy += pxmn - os;
        (btt->y).ixy += pymn - os;


    }

/* return after boundary determination if no Fourier analysis required */

    if(nfour) return 0;

/* produce fourier shape descriptors */

    nm = bcntl->nmode;
    oset=1+bcntl->obex;

    if(bcntl->fd){

/* convert boundary to float format */

      for(j=0; j < ob->bound_num; j++){

          btt = ob->bound + j;

          (btt->x).xy = obj_xreal((btt->x).ixy + x1u - 2);

          (btt->y).xy = *(gr->ygrid + (btt->y).ixy + y1u - 2);                 

      }

          

/* find number of frequencies */

      nexp = (int)(log((float)(ob->bound_num))/log(2.0));
      if(powi(2, nexp) != ob->bound_num){
         ++nexp;
      }

      if(bcntl->nexp > 0){

         if(nexp > bcntl->nexp) {

            printf("****ERROR****, boundary size to big for number of frequencies, choose a larger fixed number of frequencies.\n\n");
            exit(1);

         }

         nexp = bcntl->nexp;

      }

      else nexp += bcntl->nadd;

      nf = powi(2, nexp);

      if(nm < nf/2){

         ob->bound = (struct boundary_pt *)realloc_n(ob->bound, (nf+1)*sizeof(struct boundary_pt));
         mem_er((ob->bound == NULL) ? 0 : 1, (nf+1)*sizeof(struct boundary_pt));

         bbb = (complex *)calloc(nf, sizeof(complex));
         mem_er((bbb == NULL) ? 0 : 1, nf * sizeof(complex));

         comp((ob->bound->x).xy, (ob->bound->y).xy, bbb);

         for(j=1; j<nf; j++){

             fx = (float) (j * (ob->bound_num - 1)) / (float) nf;
             ix1 = (int) fx;

             btt = ob->bound + ix1;
             xx = (btt->x).xy + (fx - ix1) * (((btt+1)->x).xy - (btt->x).xy);
             yy = (btt->y).xy + (fx - ix1) * (((btt+1)->y).xy - (btt->y).xy);

             comp(xx, yy, bbb+j);
                            

         }

         dfft(bbb, nexp, 1, 1);


/* low pass filter fourier series */

         nmf = nf - 2 * nm - 1;

         for(j=0; j <nmf; j++){

            ii = nm + j + 1;
            cmx(0.0, *(bbb+ii), bbb+ii);

         }

         dfft(bbb, nexp, -1, 1);

         for(j=0; j<nf; j++){
            btt=ob->bound + j; 
            (btt->x).xy = (bbb+j)->real; 
            (btt->y).xy = (bbb+j)->imag; 
            btt->val = 0.0;

         } 

         btt=ob->bound + nf;
         (btt->x).xy = bbb->real;
         (btt->y).xy = bbb->imag;
         btt->val = 0.0;
 
         dfft(NULL, 0, 1, -1);

         ob->bound_num = nf + 1;

/* remove object points outside of new boundary */


         for(j=0; j < ob->point_num; j++){

            ptt = (ob->pt) + j;
            pty = ptt->y + y1u - 2;
            xx = obj_xreal(ptt->x + x1u - 2);
            if(fabs(xx - GERROR) < 1.0) continue;
            if(pty < 0 || pty >= gr->iy) continue;
            yy = *(gr->ygrid + pty);

            if(!inside(ob->bound_num, ob->bound, xx, yy)){

               for(k=j+1; k < ob->point_num; k++){

                   *ptt = *(ptt+1);
                   ++ptt;
               }

               --(ob->point_num);
               --j;
               --dfil;

           }

         } 


/* add points to object inside new boundary */

         if(bcntl->ofill == 'y'){


           xdim = pxmx - pxmn + 3 + 2*bcntl->obex;
           ydim = pymx - pymn + 3 + 2*bcntl->obex;


           dim = xdim*ydim;

           aa = (int *)realloc_n(aa, dim * sizeof(int));
           mem_er((aa == NULL) ? 0 : 1, dim * sizeof(int));

           for(j=0; j < dim; j++) *(aa + j) = 0;

           for(j=0; j < ob->point_num; j++){

               ptt = (ob->pt) + j;

               *(aa+(ptt->y - pymn + oset)*xdim + (ptt->x - pxmn + oset)) = 1;

           } 

 
           nly = nry = 0;

           oobptn = ob->point_num;    

           for(j=0; j < ydim; j++){

               iy =  j + pymn - oset;

               for(k=0; k < xdim; k++){


                   ix =  k + pxmn - oset;

                   ab = *(aa + j*xdim + k);

                   if(!ab){
                      pty = iy + y1u - 2;
                      xx = obj_xreal(ix + x1u - 2);
                      if(fabs(xx - GERROR) < 1.0) continue;
                      if(pty < 0 || pty >= gr->iy) continue;
                      yy = *(gr->ygrid + pty);
                      if(inside(ob->bound_num, ob->bound, xx, yy)){

                         ++(ob->point_num);
                         ++dfil;

                         *(aa + j*xdim + k) = 1;


                      }

                   }

               }

           }

           ob->mem_pt_size = ob->point_num;

           if(pb == 'y' && llb && rrb){

              for(j=0; j<ydim; j++){

                 iy = j + pymn - oset;

                 for(k=0; k <xdim; k++){

                     ix = k + pxmn - oset;

                     ab = *(aa + j*xdim + k);
  
                     if(ab && ix == 1){
                       if(ob->b_or_i == '\0') {ob->b_or_i = 'l';} 
                       else if(ob->b_or_i == 'r') {ob->b_or_i = 'b';}
                       *(llb + iy - 1) = ob->lab;
                       ++nly;
                     }

                     else if(ab && ix == gr->ix){
                       if(ob->b_or_i == '\0') {ob->b_or_i = 'r';} 
                       else if(ob->b_or_i == 'l') {ob->b_or_i = 'b';}
                       *(rrb + iy - 1) = ob->lab;
                       ++nry;
                     }
                  
                 }

              }

              if(ob->b_or_i == 'l') {ob->mem_pt_size += nly;}
              else if (ob->b_or_i == 'r'){ob->mem_pt_size += nry;}
              else if (ob->b_or_i == 'b'){ob->mem_pt_size = (nly > nry) ? nly : nry;}

           }            



           if(ob->mem_pt_size > oobptn){

              ob->pt = (struct point * )realloc_n(ob->pt, ob->mem_pt_size*sizeof(struct point));
              mem_er((ob->pt == NULL) ? 0 : 1, ob->mem_pt_size*sizeof(struct point));

           }


           ptt = ob->pt;

           for(j=0; j<ydim; j++){

              iy = j + pymn - oset;

              for(k=0; k <xdim; k++){

                  ix = k + pxmn - oset;

                  ab = *(aa + j*xdim + k);

                  if(ab){

                     ptt->x = ix;
                     ptt->y = iy;

                     if(ix < 0) ix += gr->ix;
                     else if(ix >= gr->ix) ix -= gr->ix;

                     ptt->val = *(ap +  (iy + y1u - 2) * gr->ix + ix + x1u - 2);
                     ++ptt; 

                  }


              }

           }

         }

         if(!ob->point_num){

            ob->bound_num = 0;
            free(ob->bound);

         }

/* remove feature points outside of new boundary */


         if(ob->fet){


            for(j=0; j < ob->fet->feature_num; j++){

                fpt = ((ob->fet->fpt) + j);
                if(tf == 3){

                   xx = obj_xreal((fpt->x).ixy + x1u - 2);
                   yy = *(gr->ygrid + (fpt->y).ixy + y1u - 2);

                }

                else {

                   xx = (fpt->x).xy;
                   yy = (fpt->y).xy;

                }

                if(!inside(ob->bound_num, ob->bound, xx, yy)){

                   for(k=j+1; k < ob->fet->feature_num; k++){

                       *fpt = *(fpt + 1);
                       ++fpt;

                   }

                   --(ob->fet->feature_num);
                   --j;


                }

                else ++nfet;


            }


         }


         free(bbb);

      }

    }

    else nfet = ob->fet->feature_num;

    free(aa);
    free(ob_field);

   
    return nfet;

}


void dir_search(int *aa, int *ix, int *iy, int *dir, int mxb, int xdim)
{

      switch(*dir){
          case 1:
            if(*(aa + *iy*xdim + *ix - 1) == mxb) {*ix -= 1; *dir = 7;}
            else if(*(aa + (*iy + 1)*xdim + *ix - 1) == mxb) {*ix -= 1; *iy += 1; *dir = 8;}
            else if(*(aa + (*iy + 1)*xdim + *ix) == mxb) {*iy += 1; *dir = 1;}
            else if(*(aa + (*iy + 1)*xdim + *ix + 1) == mxb) {*ix += 1; *iy += 1; *dir = 2;}
            else if(*(aa + *iy * xdim + *ix + 1) == mxb) {*ix += 1; *dir = 3;}
            else if(*(aa + (*iy - 1) * xdim + *ix + 1) == mxb) {*ix += 1; *iy -= 1; *dir = 4;}
            else if(*(aa + (*iy - 1) * xdim + *ix) == mxb) {*iy -= 1; *dir = 5;}
            else if(*(aa + (*iy - 1) * xdim + *ix - 1) == mxb) {*ix -= 1; *iy -= 1; *dir = 6;}
            else {*dir = 5; dir_search(aa, ix, iy, dir, mxb, xdim);}
            break;
          case 2:
            if(*(aa + (*iy + 1)*xdim + *ix - 1) == mxb) {*ix -= 1; *iy += 1; *dir = 8;}
            else if(*(aa + (*iy + 1)*xdim + *ix) == mxb) {*iy += 1; *dir = 1;}
            else if(*(aa + (*iy + 1)*xdim + *ix + 1) == mxb) {*ix += 1; *iy += 1; *dir = 2;}
            else if(*(aa + *iy * xdim + *ix + 1) == mxb) {*ix += 1; *dir = 3;}
            else if(*(aa + (*iy - 1) * xdim + *ix + 1) == mxb) {*ix += 1; *iy -= 1; *dir = 4;}
            else if(*(aa + (*iy - 1) * xdim + *ix) == mxb) {*iy -= 1; *dir = 5;}
            else if(*(aa + (*iy - 1) * xdim + *ix - 1) == mxb) {*ix -= 1; *iy -= 1; *dir = 6;}
            else if(*(aa + *iy*xdim + *ix - 1) == mxb) {*ix -= 1; *dir = 7;}
            else {*dir = 6; dir_search(aa, ix, iy, dir, mxb, xdim);};
            break;
          case 3:
            if(*(aa + (*iy + 1)*xdim + *ix) == mxb) {*iy += 1; *dir = 1;}
            else if(*(aa + (*iy + 1)*xdim + *ix + 1) == mxb) {*ix += 1; *iy += 1; *dir = 2;}
            else if(*(aa + *iy * xdim + *ix + 1) == mxb) {*ix += 1; *dir = 3;}
            else if(*(aa + (*iy - 1) * xdim + *ix + 1) == mxb) {*ix += 1; *iy -= 1; *dir = 4;}
            else if(*(aa + (*iy - 1) * xdim + *ix) == mxb) {*iy -= 1; *dir = 5;}
            else if(*(aa + (*iy - 1) * xdim + *ix - 1) == mxb) {*ix -= 1; *iy -= 1; *dir = 6;}
            else if(*(aa + *iy*xdim + *ix - 1) == mxb) {*ix -= 1; *dir = 7;}
            else if(*(aa + (*iy + 1)*xdim + *ix - 1) == mxb) {*ix -= 1; *iy += 1; *dir = 8;}
            else {*dir = 7; dir_search(aa, ix, iy, dir, mxb, xdim);};
            break;
          case 4:
            if(*(aa + (*iy + 1)*xdim + *ix + 1) == mxb) {*ix += 1; *iy += 1; *dir = 2;}
            else if(*(aa + *iy * xdim + *ix + 1) == mxb) {*ix += 1; *dir = 3;}
            else if(*(aa + (*iy - 1) * xdim + *ix + 1) == mxb) {*ix += 1; *iy -= 1; *dir = 4;}
            else if(*(aa + (*iy - 1) * xdim + *ix) == mxb) {*iy -= 1; *dir = 5;}
            else if(*(aa + (*iy - 1) * xdim + *ix - 1) == mxb) {*ix -= 1; *iy -= 1; *dir = 6;}
            else if(*(aa + *iy*xdim + *ix - 1) == mxb) {*ix -= 1; *dir = 7;}
            else if(*(aa + (*iy + 1)*xdim + *ix - 1) == mxb) {*ix -= 1; *iy += 1; *dir = 8;}
            else if(*(aa + (*iy + 1)*xdim + *ix) == mxb) {*iy += 1; *dir = 1;}
            else {*dir = 8; dir_search(aa, ix, iy, dir, mxb, xdim);};
            break;
          case 5:
            if(*(aa + *iy * xdim + *ix + 1) == mxb) {*ix += 1; *dir = 3;}
            else if(*(aa + (*iy - 1) * xdim + *ix + 1) == mxb) {*ix += 1; *iy -= 1; *dir = 4;}
            else if(*(aa + (*iy - 1) * xdim + *ix) == mxb) {*iy -= 1; *dir = 5;}
            else if(*(aa + (*iy - 1) * xdim + *ix - 1) == mxb) {*ix -= 1; *iy -= 1; *dir = 6;}
            else if(*(aa + *iy*xdim + *ix - 1) == mxb) {*ix -= 1; *dir = 7;}
            else if(*(aa + (*iy + 1)*xdim + *ix - 1) == mxb) {*ix -= 1; *iy += 1; *dir = 8;}
            else if(*(aa + (*iy + 1)*xdim + *ix) == mxb) {*iy += 1; *dir = 1;}
            else if(*(aa + (*iy + 1)*xdim + *ix + 1) == mxb) {*ix += 1; *iy += 1; *dir = 2;}
            else {*dir = 1; dir_search(aa, ix, iy, dir, mxb, xdim);};
            break;
          case 6:
            if(*(aa + (*iy - 1) * xdim + *ix + 1) == mxb) {*ix += 1; *iy -= 1; *dir = 4;}
            else if(*(aa + (*iy - 1) * xdim + *ix) == mxb) {*iy -= 1; *dir = 5;}
            else if(*(aa + (*iy - 1) * xdim + *ix - 1) == mxb) {*ix -= 1; *iy -= 1; *dir = 6;}
            else if(*(aa + *iy*xdim + *ix - 1) == mxb) {*ix -= 1; *dir = 7;}
            else if(*(aa + (*iy + 1)*xdim + *ix - 1) == mxb) {*ix -= 1; *iy += 1; *dir = 8;}
            else if(*(aa + (*iy + 1)*xdim + *ix) == mxb) {*iy += 1; *dir = 1;}
            else if(*(aa + (*iy + 1)*xdim + *ix + 1) == mxb) {*ix += 1; *iy += 1; *dir = 2;}
            else if(*(aa + *iy * xdim + *ix + 1) == mxb) {*ix += 1; *dir = 3;}
            else {*dir = 2; dir_search(aa, ix, iy, dir, mxb, xdim);};
            break;
          case 7:
            if(*(aa + (*iy - 1) * xdim + *ix) == mxb) {*iy -= 1; *dir = 5;}
            else if(*(aa + (*iy - 1) * xdim + *ix - 1) == mxb) {*ix -= 1; *iy -= 1; *dir = 6;}
            else if(*(aa + *iy*xdim + *ix - 1) == mxb) {*ix -= 1; *dir = 7;}
            else if(*(aa + (*iy + 1)*xdim + *ix - 1) == mxb) {*ix -= 1; *iy += 1; *dir = 8;}
            else if(*(aa + (*iy + 1)*xdim + *ix) == mxb) {*iy += 1; *dir = 1;}
            else if(*(aa + (*iy + 1)*xdim + *ix + 1) == mxb) {*ix += 1; *iy += 1; *dir = 2;}
            else if(*(aa + *iy * xdim + *ix + 1) == mxb) {*ix += 1; *dir = 3;}
            else if(*(aa + (*iy - 1) * xdim + *ix + 1) == mxb) {*ix += 1; *iy -= 1; *dir = 4;}
            else {*dir = 3; dir_search(aa, ix, iy, dir, mxb, xdim);};
            break;
          case 8:
            if(*(aa + (*iy - 1) * xdim + *ix - 1) == mxb) {*ix -= 1; *iy -= 1; *dir = 6;}
            else if(*(aa + *iy*xdim + *ix - 1) == mxb) {*ix -= 1; *dir = 7;}
            else if(*(aa + (*iy + 1)*xdim + *ix - 1) == mxb) {*ix -= 1; *iy += 1; *dir = 8;}
            else if(*(aa + (*iy + 1)*xdim + *ix) == mxb) {*iy += 1; *dir = 1;}
            else if(*(aa + (*iy + 1)*xdim + *ix + 1) == mxb) {*ix += 1; *iy += 1; *dir = 2;}
            else if(*(aa + *iy * xdim + *ix + 1) == mxb) {*ix += 1; *dir = 3;}
            else if(*(aa + (*iy - 1) * xdim + *ix + 1) == mxb) {*ix += 1; *iy -= 1; *dir = 4;}
            else if(*(aa + (*iy - 1) * xdim + *ix) == mxb) {*iy -= 1; *dir = 5;}
            else {*dir = 4; dir_search(aa, ix, iy, dir, mxb, xdim);};
            break;
      }

      return;

}

/* function to check if point is inside a boundary, 2D only. */

int inside(int nm, struct boundary_pt *bdpt, float xx, float yy)
{


    int i;
    int inb=0;

    float x1, x2, y1, y2;
    float cp, dp;
    float atot=0.0;

    struct boundary_pt *pt1, *pt2;


    for(i=0; i < nm-1; i++){

       pt1 = bdpt + i;
       pt2 = pt1 + 1;

       x1 = (pt1->x).xy - xx;
       x2 = (pt2->x).xy - xx;
       y1 = (pt1->y).xy - yy;
       y2 = (pt2->y).xy - yy;

       cp = x1 * y2 - x2 * y1;
       dp = x1 * x2 + y1 * y2;

       atot += atan2(cp, dp);

    }
 
    inb = (fabs(atot) > FPI) ? 1 : 0;

    return inb;

}
