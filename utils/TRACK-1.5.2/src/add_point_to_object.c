#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"
#include "st_obj.h"

/* function to add points that are outside the boundarys but which belong 
   to the object.                                                         */

extern int x1u, y1u;

void add_point_to_object(struct frame_objs *fo, int yc, int xc, int lab, float val)

{

    struct object *ob=NULL;
    struct point *ptt=NULL;

    ob = (fo->objs) + lab -1;

    ++(ob->point_num);

    ob->pt = (struct point * )realloc_n(ob->pt, (ob->point_num)*sizeof(struct point));
    mem_er((ob->pt == NULL) ? 0 : 1, (ob->point_num)*sizeof(struct point));

    ptt = (ob->pt) + (ob->point_num) - 1;

    ptt->x = xc+2-x1u;
    ptt->y = yc+2-y1u;
    ptt->val = val;

    return;

}
