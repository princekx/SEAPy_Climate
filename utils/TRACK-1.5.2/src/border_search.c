#include <Stdio.h>
#include <stdlib.h>
#include "st_obj.h"
#include "grid.h"

void search_adjacency(struct frame_objs * ,int * ,int ,int ,int ,float ); 

/* function to find connected components outside the region of interest
   for objects on the boundary, using an omnidirectional region growing
   technique.
                                                           */
extern int x1u, x2u, y1u, y2u;
extern int cc;
extern int perh, dbp, dbm;
/*extern int perv; */

extern GRID *gr;

void border_search(struct frame_objs *fo, int *ib, int yc, int xc, float thresh)

{ 

     int xcor, ycor;
     int xp1, xp2, yp1, yp2, lab1;
     int prx, pry;

     if(x1u < dbp && xc > x2u+dbm) 
        prx = (xc - x1u + dbp - gr->ix + 1)*(x2u + dbm - xc + gr->ix - 1);

     else if(x2u > gr->ix - dbp && xc < x1u-dbm)
        prx = (xc - x1u + dbp + gr->ix - 1)*(x2u + dbm - xc - gr->ix + 1);

     else prx = (xc - x1u + dbp)*(x2u + dbm - xc);

     pry = (yc - y1u + dbp)*(y2u + dbm - yc);    

     if(prx > 0 && pry > 0){

        lab1 = *(ib + yc * gr->ix + xc);
        
        if(perh == 1){

          if(xc <= 0) xcor = gr->ix - 2;

          else  xcor = xc - 1;

          search_adjacency(fo, ib, lab1, xcor, yc, thresh);

        }

        else if(xc > 0 && xc <= gr->ix - 1){

             xcor = xc - 1;
             search_adjacency(fo, ib, lab1, xcor, yc, thresh);

        }

        else xcor = xc;

        xp1 = xcor;

        if(perh == 1){

          if(xc >= gr->ix - 1) xcor = 1;

          else  xcor = xc + 1;

          search_adjacency(fo, ib, lab1, xcor, yc, thresh);

        }

        else if(xc >= 0 && xc < gr->ix - 1){

             xcor = xc+1;
             search_adjacency(fo, ib, lab1, xcor, yc, thresh);

        }

        else xcor = xc;

        xp2 = xcor; 

        if(yc > 0 && yc <= gr->iy - 1){

           ycor = yc - 1;

           search_adjacency(fo, ib, lab1, xc, ycor, thresh);

        }
 
        else ycor = yc;

        yp1 = ycor;


        if(yc >= 0 && yc < gr->iy - 1){

          ycor = yc+1;

          search_adjacency(fo, ib, lab1, xc, ycor, thresh);

        }

        else ycor = yc;

        yp2 = ycor;

        if(cc == 'v'){

          search_adjacency(fo, ib, lab1, xp1, yp1, thresh);
          search_adjacency(fo, ib, lab1, xp1, yp2, thresh);
          search_adjacency(fo, ib, lab1, xp2, yp1, thresh);
          search_adjacency(fo, ib, lab1, xp2, yp2, thresh);


        }

     }

     return;

}
