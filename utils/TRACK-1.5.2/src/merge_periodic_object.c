#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_obj.h"
#include "grid.h"

/* function to merge pieces of a periodic object */

extern int x1u, y1u;

extern GRID *gr;

void merge_periodic_object(struct frame_objs *fo, int *lb, int *rb, int lab1, int lab2)

{

   int i;
   int lb1, lb2, spt, npt, tpt;
   int ixx, inx;
   int *ll=NULL;
   int ix1;

   struct object *ob1=NULL, *ob2=NULL, *ob=NULL, *cob=NULL;
   struct point *pt1=NULL, *pt2=NULL;

   ix1 = gr->ix - 1;

   ob1 = (fo->objs) + lab1 - 1;
   ob2 = (fo->objs) + lab2 - 1;

   if(ob1->point_num > ob2->point_num){

     ob = ob1;
     cob = ob2;
     lb1 = lab2;
     lb2 = lab1;
     ll = rb;
     inx = gr->ix;

     spt = ob1->point_num;
     npt = ob2->point_num;
     tpt = ob2->mem_pt_size;

   }

   else {
 
     ob = ob2;
     cob = ob1;
     lb1 = lab1;
     lb2 = lab2;
     ll = lb;
     inx = 1;

     spt = ob2->point_num;
     npt = ob1->point_num;
     tpt = ob1->mem_pt_size;

   }

   ob->point_num += npt;
   ob->mem_pt_size += tpt;
 
   ob->pt = (struct point * )realloc_n(ob->pt, (ob->mem_pt_size)*sizeof(struct point));
    mem_er((ob->pt == NULL) ? 0 : 1, (ob->mem_pt_size)*sizeof(struct point));

   pt1 = (ob->pt) + spt;
   pt2 = cob->pt;

   for(i=0; i < npt; i++){

       ixx = x1u + pt2->x - 1;

       if(ixx == inx) {*(ll + y1u + pt2->y - 2) = lb2; --(ob->point_num);}

       else {

         if(cob->b_or_i == 'r') pt2->x -= ix1;

         else if(cob->b_or_i == 'l') pt2->x += ix1;

         *(pt1++) = *pt2;

       }

       ++pt2;       

   }

   ob = (fo->objs) + lb1 - 1;

   ob->mem_pt_size = ob->point_num = 0;
   ob->b_or_i = 'n';

   return;

}
