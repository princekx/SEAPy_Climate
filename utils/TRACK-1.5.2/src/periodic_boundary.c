#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_obj.h"
#include "grid.h"

extern int x1u, x2u, y1u, y2u;

extern GRID *gr;

/* function to determine periodic objects in X-direction only using 
   boundary matching.                                                */

void merge_periodic_object(struct frame_objs * , int * , int * , int , int );

void periodic_boundary(struct frame_objs *fo, int *llb, int *rrb)

{


   int *lb=NULL, *rb=NULL;
   int l1, r1;
   int i, j, ii, ip, ipx;
   int ilb;

   struct object *ob=NULL;
   struct point *ptt=NULL;

/* assign initialized memory for boundary matching */

   lb = (int * )calloc(gr->iy, sizeof(int));
   mem_er((lb == NULL) ? 0 : 1, gr->iy * sizeof(int));

   rb = (int * )calloc(gr->iy, sizeof(int));
   mem_er((rb == NULL) ? 0 : 1, gr->iy * sizeof(int));

   for(i=0; i < fo->obj_num; i++) {

       ob = (fo->objs) + i;

       ii = i + 1;

       if(ob->b_or_i == 'x' || ob->b_or_i == 'b'){

          ob->mem_pt_size = ob->point_num;

          ilb = 0;

          for(j=0; j < ob->point_num; j++){

              ptt = (ob->pt) + j;

              ipx = x1u + ptt->x - 1;

              ip = y1u + ptt->y - 2;

              if(ipx == 1){

                *(llb + ip) = *(lb + ip) = ii;
                if(!ilb) {ob->b_or_i = 'l'; ilb = 1;}
                else if(ilb && ob->b_or_i == 'r') ob->b_or_i = 'b';

              }

              else if(ipx == gr->ix){

                *(rrb + ip) = *(rb + ip) = ii;
                if(!ilb) {ob->b_or_i = 'r'; ilb = 1;}
                else if(ilb && ob->b_or_i == 'l') ob->b_or_i = 'b';

              }

          }

       }

   }

   for(i=y1u-1; i<y2u; i++){

       l1 = *(lb+i);
       r1 = *(rb+i);

       if(l1*r1 && l1 != r1) merge_periodic_object(fo, lb, rb, l1, r1);

   }

/*   for(i=0; i < gr->iy; i++)printf("%d %d\n", *(lb+i), *(rb+i)); */

   free(rb);

   free(lb);

   return;

}
