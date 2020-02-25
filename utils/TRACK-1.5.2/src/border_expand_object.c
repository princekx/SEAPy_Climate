#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_obj.h"
#include "grid.h"

void border_search(struct frame_objs * ,int * ,int , int , float );

/* function to expand the object regions for those objects on a boundary */ 
 
extern int x1u, y1u, x2u, y2u;

extern GRID *gr;

void border_expand_object(struct frame_objs *fo, int xdim, int ydim, float thresh)

{

     int i, j, bi;

     int *ib, *ibb;

     struct object *ob;
     struct point *ptt;

/* assign memory for boundary searching block */

     ib = (int * )calloc(gr->ix * gr->iy, sizeof(int));
     mem_er((ib == NULL) ? 0 : 1, (gr->ix * gr->iy * sizeof(int)));   

     for(i=0; i < fo->obj_num; i++){

         ob = (fo->objs) + i;

         bi = 0;

/* is object on the boundary                                              
   x == X-boundary; y == Y-boundary; b == both boundarys */

         for(j=0; j < ob->point_num; j++){

             ptt = (ob->pt) + j;

             if(ptt->x == 1){

               if(bi == 0){ bi = 1; ob->b_or_i = 'x';}

               else if(bi == 1 && ob->b_or_i == 'y') ob->b_or_i = 'b';

             }

             else if(ptt->x == xdim){

                if(bi == 0) {bi = 1; ob->b_or_i = 'x'; }

                else if(bi == 1 && ob->b_or_i == 'y') ob->b_or_i = 'b';

             }

             if(ptt->y == 1){ 

                if(bi == 0) {bi = 1; ob->b_or_i = 'y'; }

                else if(bi == 1 && ob->b_or_i == 'x') ob->b_or_i = 'b';

             }

             else if(ptt->y == ydim){

                if(bi == 0) {bi = 1; ob->b_or_i = 'y'; }

                else if(bi == 1 && ob->b_or_i == 'x') ob->b_or_i = 'b';

             }

         }

         if(bi == 1){

            for(j=0; j < ob->point_num; j++){

               ptt = (ob->pt) + j;

               ibb = ib + (y1u + ptt->y -2) * gr->ix + (x1u + ptt->x -2);

               *ibb = i+1;

            }

         }

         else ob->b_or_i = 'i';

     }

/*     for(i=y1u-1; i < y2u; i++){
        for(j=x1u-1; j < x2u; j++){

           printf("%3d ", *(ib + i * gr->ix +j));

        }
        printf("\n");

     }
     printf("\n"); */

     for(j=y1u - 1; j < y2u; j++){

        if(*(ib + j * gr->ix + x1u-1) > 0)

           border_search(fo, ib, j, x1u-1, thresh); 

        if(*(ib + j * gr->ix + x2u-1) > 0)

           border_search(fo, ib, j, x2u-1, thresh); 

     }              

     for(j=x1u-1; j < x2u; j++){

        if(*(ib + (y1u-1) * gr->ix + j) > 0)

           border_search(fo, ib, y1u-1, j, thresh);

        if(*(ib + (y2u-1) * gr->ix + j) > 0)

           border_search(fo, ib, y2u-1, j, thresh); 

     }

/*     for(i=y1u-1; i < y2u; i++){
        for(j=x1u-1; j < x2u; j++){

           printf("%3d ", *(ib + i * gr->ix +j));

        }
        printf("\n");

     } */

     free(ib);

     return;

}
