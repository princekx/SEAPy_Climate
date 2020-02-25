#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_im.h"
#include "st_obj.h"
#include "st_fo.h"
#include "grid.h"

#define  LESS(A, B, C) (int)((A <= B) ? (C + 1) : C)


void border_obj_fet_filter(struct object * , double ** , struct frame_objs * );

/* void object_split_init(struct object * , struct frame_objs * ); */

void hierarc_segment(struct image ** , struct image * , int , struct frame_objs * , int , int , int );
int powi(int , int );

extern int x1u, y1u;
extern int bs;

extern GRID *gr;

/* function to calculate the local max's of each object in each frame */

void object_local_maxs(struct frame_objs *ff, struct frame_objs *fp, int tf, int b_exc) 

{

     int i, j, k, ln, dcount=0;
     int pxmx, pxmn, pymx, pymn, ydim, xdim;
     int xd;
     int lmax, level, p2;

     float *aa=NULL, *ab=NULL;
     float sx, sy;

     double **fuz=NULL;

     struct image *ia=NULL, **hierarchy=NULL;
     struct object *ob=NULL, *obf=NULL;
     struct point *ptt=NULL;
     struct feature_pts *fpts=NULL;

     for(i=0; i < ff->obj_num; i++){

         ob = (ff->objs)+i;

         ob->fet = (struct features *)malloc_initl(sizeof(struct features));
         mem_er((ob->fet == NULL) ? 0 : 1, sizeof(struct features));

         ob->fet->fpt = (struct feature_pts *)malloc_initl(sizeof(struct feature_pts));
         mem_er((ob->fet->fpt == NULL) ? 0 : 1, sizeof(struct feature_pts));

     }

     ff->tot_f_f_num = 0;

     for(i=0; i < ff->obj_num; i++){

         ob = (ff->objs) + i;
         ln = 0;

         if(!ob->ext){

            pxmx = pxmn = ob->pt->x;
            pymx = pymn = ob->pt->y;

            for(j=0; j < ob->point_num; j++){

                ptt = (ob->pt) + j;

                if(ptt->x > pxmx) pxmx = ptt->x;
                else if(ptt->x < pxmn) pxmn = ptt->x;

                if(ptt->y > pymx) pymx = ptt->y;
                else if(ptt->y < pymn) pymn = ptt->y;

             }

/* assign extent of objects if required */

             if(tf >= 7) {

               ob->ext = (struct extent * )malloc_initl(sizeof(struct extent));
               mem_er((ob->ext == NULL) ? 0 : 1, sizeof(struct extent));

               ob->ext->x1 = pxmn;
               ob->ext->x2 = pxmx;
               ob->ext->y1 = pymn;
               ob->ext->y2 = pymx;
  
             }

          }

          else {

               pxmn = ob->ext->x1;
               pxmx = ob->ext->x2;
               pymn = ob->ext->y1;
               pymx = ob->ext->y2;

          }


          xdim = pxmx - pxmn + 3;
          ydim = pymx - pymn + 3;

          aa = (float *)calloc(xdim*ydim, sizeof(float));
          mem_er((aa == NULL) ? 0 : 1, xdim*ydim * sizeof(float));

          for(j=0; j < ob->point_num; j++){

              ptt = (ob->pt) + j;

              *(aa+(ptt->y - pymn+1)*xdim + (ptt->x - pxmn +1)) = ptt->val;

          }
 
         for(j=1; j < ydim-1; j++){

             for(k=1; k < xdim-1; k++){

                ab = (aa+j*xdim+k);
                dcount = 0;

                if(*ab > 0.0){

                   if(b_exc){

                      if(*(ab-1) > 0.0)dcount = LESS(*(ab-1), *ab, dcount); 
                      if(*(ab+1) > 0.0)dcount = LESS(*(ab+1), *ab, dcount);
                      if(*(ab - xdim) > 0.0)dcount = LESS(*(ab - xdim), *ab, dcount);
                      if(*(ab + xdim) > 0.0)dcount = LESS(*(ab + xdim), *ab, dcount);
                      if(*(ab-xdim-1) > 0.0)dcount = LESS(*(ab-xdim-1), *ab, dcount);
                      if(*(ab+xdim+1) > 0.0)dcount = LESS(*(ab+xdim+1), *ab, dcount);
                      if(*(ab+xdim-1) > 0.0)dcount = LESS(*(ab+xdim-1), *ab, dcount);
                      if(*(ab-xdim+1)> 0.0)dcount = LESS(*(ab-xdim+1), *ab, dcount);

                   }

                   else {

                      dcount = LESS(*(ab-1), *ab, dcount); 
                      dcount = LESS(*(ab+1), *ab, dcount);
                      dcount = LESS(*(ab - xdim), *ab, dcount);
                      dcount = LESS(*(ab + xdim), *ab, dcount);
                      dcount = LESS(*(ab-xdim-1), *ab, dcount);
                      dcount = LESS(*(ab+xdim+1), *ab, dcount);
                      dcount = LESS(*(ab+xdim-1), *ab, dcount);
                      dcount = LESS(*(ab-xdim+1), *ab, dcount);

                   }

                   if(dcount == 8){

                      ++(ob->fet->feature_num);
                      ++(ff->tot_f_f_num);
 
                      if(ob->fet->feature_num > 1){
                          ob->fet->fpt = (struct feature_pts *)realloc_n(ob->fet->fpt, (ob->fet->feature_num)*sizeof(struct feature_pts));
                          mem_er((ob->fet->fpt == NULL) ? 0 : 1, (ob->fet->feature_num)*sizeof(struct feature_pts));

                      }

                      (((ob->fet->fpt)+ln)->x).ixy = (((ob->fet->fpt)+ln)->ox).ixy = k+pxmn-1;
                      (((ob->fet->fpt)+ln)->y).ixy = (((ob->fet->fpt)+ln)->oy).ixy = j+pymn-1;
                      ((ob->fet->fpt)+ln)->str = ((ob->fet->fpt)+ln)->ostr = *ab;
                      ((ob->fet->fpt)+ln)->track_id = 0;
		      ((ob->fet->fpt)+ln)->add_fld = NULL;
                      ++ln;

                   }


                } 


             }


          } 

          if(tf >= 4){

/* analyse object feature points for groupings and calculate the means */         

             if(ob->fet->feature_num > 1){

                xd = xdim;
                level = 1;
                lmax = (xdim >= ydim) ? xdim : ydim;
                p2 = powi(2, level);
 
                while(p2 <= lmax){
                      ++level;
                      p2 = powi(2, level);
                }

                if(p2 != lmax) ++level;

                xdim = ydim = p2;

/* assign memory to hierarchy if the relevent options are applicable */

                hierarchy = (struct image ** )calloc(level, sizeof(struct image *));
                mem_er((hierarchy == NULL) ? 0 : 1, level * sizeof(struct image *));

                ia = (struct image *)calloc(xdim*ydim, sizeof(struct image));
                mem_er((ia == NULL) ? 0 : 1, xdim*ydim * sizeof(struct image));         

                for(j=0; j < ob->fet->feature_num; j++){

                    fpts = (ob->fet->fpt) + j;

                    (ia+((fpts->y).ixy - pymn)*xdim + ((fpts->x).ixy - pxmn))->pval = 1;

                }

                fp->obj_num=0;
                fp->objs = NULL;
                fp->tot_f_f_num = 0;

                hierarc_segment(hierarchy, ia, level, fp, xdim, ydim, 'v');

                free(ia);

                ff->tot_f_f_num += (fp->obj_num - ob->fet->feature_num);
                ob->fet->feature_num = fp->obj_num;

                ob->fet->fpt = (struct feature_pts * )realloc_n(ob->fet->fpt, fp->obj_num*sizeof(struct feature_pts));
                mem_er((ob->fet->fpt == NULL) ? 0 : 1, fp->obj_num*sizeof(struct feature_pts));

                for(j=0; j < fp->obj_num; j++){

                   obf = (fp->objs)+j;

                   sx = 0.0;
                   sy = 0.0;

                   for(k=0; k < obf->point_num; k++){

                       sx += obj_xreal((((obf->pt)+k)->x) + pxmn +x1u - 3);
                       sy += *(gr->ygrid + (((obf->pt)+k)->y) + pymn +y1u - 3);

                   }

                   (((ob->fet->fpt)+j)->x).xy = (((ob->fet->fpt)+j)->ox).xy = sx/obf->point_num;
                   (((ob->fet->fpt)+j)->y).xy = (((ob->fet->fpt)+j)->oy).xy = sy/obf->point_num;
                   ((ob->fet->fpt)+j)->str = ((ob->fet->fpt)+j)->ostr = *(aa+ ((obf->pt->y))*xd + (obf->pt->x));

                   free(obf->pt);

                }
                   

                fp->obj_num = 0;
                fp->tot_f_f_num =0;
                free(fp->objs);
                fp->objs = NULL;
                free(hierarchy);

             }

             else if (ob->fet->feature_num == 1){

                ((ob->fet->fpt)->x).xy = ((ob->fet->fpt)->ox).xy = obj_xreal(((ob->fet->fpt)->x).ixy +x1u - 2);
                ((ob->fet->fpt)->y).xy = ((ob->fet->fpt)->oy).xy = *(gr->ygrid + ((ob->fet->fpt)->y).ixy +y1u - 2);
                (ob->fet->fpt)->ostr = (ob->fet->fpt)->str;

             }

         }

         free(aa);

         if(bs == 'y' && ( tf == 3 || tf == 4)) border_obj_fet_filter(ob, fuz, ff);
/*         if(tf == 5 || tf == 6) object_split_init(ob, ff); */
 

     }
  
     return;

}
