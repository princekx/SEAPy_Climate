#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_obj.h"

/* function to filter out objects below a user defined size */

int dfil;

extern float obr;
extern int bs, pb;

void object_filter(struct frame_objs *ff , int filt_pt_num, int xdim, int ydim)

{

     int i, j, k, l, bod;
     int ptx, pty;
     int sob=0, msiz=0;

     double div, totm, intm;

     struct object *ob=NULL, *ob1=NULL;
     struct point *ptt=NULL;

     dfil = 0;
     bod = 0;
     div = 0.0;

     for (i=0; i < ff->obj_num; i++){

         l=0;

         reset:

         div = 0.0;

         ob = (ff->objs) + i;

         if(pb == 'y')

           sob = (ob->b_or_i == 'l' || 
                  ob->b_or_i == 'r' || ob->b_or_i == 'b') ? 1 : 0;

         msiz = (sob) ? ob->mem_pt_size : ob->point_num;

         if(bs == 'y' && obr > 0.0){

            totm = intm = 0;
            bod = 0;

            if(ob->b_or_i != 'i' && ob->point_num != 0){

               for(j=0; j < ob->point_num; j++){

                   ptt = (ob->pt) + j;
                   ptx = (ptt->x - 1)*(xdim - ptt->x);
                   pty = (ptt->y - 1)*(ydim - ptt->y);

                   if(ptx >= 0 && pty >= 0) intm += (float)ptt->val;

                   totm += (float)ptt->val;

                }

                div = intm * 100.0/totm;

                if(div < obr) bod = 1;

            }

         }


         if(ob->point_num <= filt_pt_num || bod == 1){

            ++l;
            k=0;

            while(i+k < (ff->obj_num)-l){

                 ob = (ff->objs) + (i+k);
                 ob1 = ob + 1;

                 if(pb == 'y')

                   sob = (ob1->b_or_i == 'l' ||
                          ob1->b_or_i == 'r' || ob1->b_or_i == 'b' ) ? 1 : 0;

                 msiz = (sob) ? ob1->mem_pt_size : ob1->point_num;

                 ob->pt = (struct point *)realloc_n(ob->pt, msiz*sizeof(struct point));
                 mem_er((ob->pt == NULL) ? 0 : 1, msiz*sizeof(struct point));


                 for(j=0; j < ob1->point_num; j++)

                     *((ob->pt)+j) = *((ob1->pt)+j);

                 ob->point_num = ob1->point_num;
                 ob->mem_pt_size = ob1->mem_pt_size;

                 ob->lab = ob1->lab;

                 ob->b_or_i = ob1->b_or_i;

                 ++k;
            }

            if(i+l < ff->obj_num) goto reset;

          }

          else dfil += msiz;

          if(l > 0){

            ff->obj_num -= l;

            for(j=(ff->obj_num); j < (ff->obj_num)+l; j++){

                free(((ff->objs)+j)->pt);

            }

            ff->objs = (struct object *)realloc_n(ff->objs, (ff->obj_num)*sizeof(struct object));
            mem_er((ff->objs == NULL) ? 0 : 1, (ff->obj_num)*sizeof(struct object));

          }

      }

      return;

}
