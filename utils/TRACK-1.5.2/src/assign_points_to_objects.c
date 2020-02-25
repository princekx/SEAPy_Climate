#include <Stdio.h>
#include "st_obj.h"
#include "grid.h"

int powi(int , int );

/* function to assign point structures to object structures */

extern float *ap;
extern int x1u, y1u;
extern GRID *gr;

void assign_points_to_objects(struct object *ob, int cpt, int ll, int xp, int yp)

{

    int i, j, num=powi(2, ll);
    int ixx, iyy;

    struct point *pnt, *ptt;

    pnt=(ob->pt)+cpt;

    for(j=0; j < num; j++){

       iyy = num*yp+j;

       for(i=0; i < num; i++){

           ptt = pnt + j*num +i;

           ixx = num*xp+i;
           ptt->x = ixx+1;
           ptt->y = iyy+1;
           ptt->val = *(ap + (y1u+iyy-1) * gr->ix + (x1u+ixx -1));

/*           printf("%d %d %e\n", ptt->x, ptt->y, ptt->val); */


        }

     }

     return;

}
