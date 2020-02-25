#include <Stdio.h>
#include <stdlib.h>
#include "st_obj.h"
#include "grid.h"

void add_point_to_object(struct frame_objs * , int , int , int , float );
void merge_objects(struct frame_objs * , int * , int , int );
void border_search(struct frame_objs * , int * , int , int , float );

extern float *ap;

extern GRID *gr;

/* function for combining points with objects and objects with objects
   for objects wich intersect the boundarys or if periodic boundarys exist */

void search_adjacency(struct frame_objs *fo, int *ib, int lab1, int xcor, int ycor, float thresh)

{

     float *apt;

     int lab2, plus;

     plus = ycor * gr->ix + xcor;
     lab2 = *(ib+plus);

     if(lab2 == 0){

       apt = ap+plus;

       if(*apt >= thresh){

          *(ib+plus) = lab1;

          add_point_to_object(fo, ycor, xcor, lab1, *apt);

          border_search(fo, ib, ycor, xcor, thresh);

       } 

       else *(ib+plus) = -1;

     }

     else if(lab2 > 0 && lab2 != lab1) merge_objects(fo, ib, lab1, lab2);

     return;

}
