#include <Stdio.h>
#include <stdlib.h>
#include "st_im.h"
#include "st_obj.h"

void assign_object(struct image * , struct frame_objs * , int , int , int );

/* function to partition labeled map into objects */

void form_objects(struct image *iq, struct frame_objs *ff, int ll, int xp, int yp)

{

     if(iq->pval == -1){

        form_objects(iq->nw, ff, ll-1, 2*xp, 2*yp);
        form_objects(iq->ne, ff, ll-1, 2*xp+1, 2*yp);
        form_objects(iq->sw, ff, ll-1, 2*xp, 2*yp+1);
        form_objects(iq->se, ff, ll-1, 2*xp+1, 2*yp+1);

     }

     else if(iq->pval == 1) 

          assign_object(iq, ff, ll, xp, yp);

     return;

}
