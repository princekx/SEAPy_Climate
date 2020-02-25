#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_obj.h"
#include "grid.h"

extern int y1u, x1u;

extern GRID *gr;

/* function to merge objects on the boundarys */

void merge_objects(struct frame_objs *fo, int *ib, int lab1, int lab2)

{

    int i, spt, npt, lb1, lb2;

    struct object *ob1, *ob2, *ob, *cob;
    struct point *pt1, *pt2;

    ob1 = (fo->objs) + lab1 - 1;
    ob2 = (fo->objs) + lab2 - 1;

    if(lab1 > lab2){

       ob = ob2;
       cob = ob1;
       lb1 = lab1;
       lb2 = lab2;

       spt = ob2->point_num;
       npt = ob1->point_num;

       ob->point_num += npt;


    }

    else {

        ob = ob1;
        cob = ob2;
        lb1 = lab2;
        lb2 = lab1;

        spt = ob1->point_num;
        npt = ob2->point_num;

        ob->point_num += npt;

    }

    ob->pt = (struct point * )realloc_n(ob->pt, (ob->point_num)*sizeof(struct point));
    mem_er((ob->pt == NULL) ? 0 : 1, (ob->point_num)*sizeof(struct point));

    for(i=0; i < npt; i++){

        pt1 = (ob->pt) + spt + i;
        pt2 = (cob->pt) + i;

        *pt1 = *pt2;
        *(ib + (y1u + (pt2->y) - 2) * gr->ix + (x1u + (pt2->x) -2)) = lb2;

    }

    if(ob->b_or_i != cob->b_or_i) ob->b_or_i = 'b';

/* re-allocate storage to zero size for one of the pre-merged objects */

    ob = (fo->objs) + lb1 - 1;

    ob->point_num = 0;
    ob->b_or_i ='n';

    return;

}
