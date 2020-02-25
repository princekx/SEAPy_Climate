#include <Stdio.h>
#include <stdlib.h>
#include "st_im.h"

/* recursive function for assigning labels to data nodes in the 
   hierarchical tree.                                            */

struct image *edge_neighbour(struct image * , int);
struct image *vertex_neighbour(struct image * , int);
void label_adj(struct image * , struct image * , int , int );

int assign_label(struct image *ih, int nclass, int cc)

{

    struct image *iq;
    int q1, q2, val;

    val=ih->pval;

    if(val == -1) {

       nclass = assign_label(ih->nw, nclass, cc);
       nclass = assign_label(ih->ne, nclass, cc);
       nclass = assign_label(ih->sw, nclass, cc);
       nclass = assign_label(ih->se, nclass, cc);
    }

    else if(val == 1) {

/*       ih->label->name=0; */

       iq = edge_neighbour(ih, 'W');

       if(iq != NULL && iq->pval != 0)

          label_adj(ih, iq, q1=2, q2=4);

       iq = edge_neighbour(ih, 'N');

       if(iq != NULL && iq->pval != 0)

          label_adj(ih, iq, q1=3, q2=4);

       if(cc == 'v'){

          iq = vertex_neighbour(ih, 'V');

          if(iq != NULL && iq->pval != 0)

              label_adj(ih, iq, q1=0, q2=4);

          iq = vertex_neighbour(ih, 'U');

          if(iq != NULL && iq->pval != 0)

              label_adj(ih, iq, q1=3, q2=0);

       }

       if(ih->label->name == 0) {

          ++nclass;

          ih->label->name=nclass;

        }

     }

     return nclass;

}
