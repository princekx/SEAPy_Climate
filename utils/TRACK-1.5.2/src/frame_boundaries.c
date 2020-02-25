#include <Stdio.h>
#include "st_obj.h"
#include "st_fo.h"
#include "boundary.h"
#include "grid.h"

int boundary_find(struct object * , struct boundary_cntl * , int * , int * , int );

/* function to call boundary_find for all objects in a frame */

extern int pb;

extern GRID *gr;

void frame_boundaries(struct frame_objs *ff, struct boundary_cntl *bcntl, int *llb, int *rrb)
{

    int i;


    ff->b_state = (bcntl->fd) ? 1 : 0;

    ff->tot_f_f_num = 0;

    if(pb == 'y'){

      for(i=0; i < gr->iy; i++) *(llb + i) = *(rrb + i) = 0;

    }

    for(i=0; i< ff->obj_num; i++){

        ff->tot_f_f_num += boundary_find(ff->objs + i, bcntl, llb, rrb, 0);

    }


    return;

}
