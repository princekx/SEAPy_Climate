#include <Stdio.h>
#include "st_im.h"

/* function to propogate equivelences through the tree for the
   hierarchical data structure.                                        */

struct eq_class *find(struct eq_class *);

void update(struct image *iq)

{

     if(iq->pval == -1){

        update(iq->nw);
        update(iq->ne);
        update(iq->sw);
        update(iq->se);

     }

     else if(iq->pval == 1) iq->label=find(iq->label);

     return;

}
