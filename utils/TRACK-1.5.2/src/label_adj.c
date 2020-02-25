#include <Stdio.h>
#include "st_im.h"

/* function to find all descendents of a node adjacent to the current node  */

struct eq_class *unify(struct eq_class * , struct eq_class * );
struct eq_class *find(struct eq_class * );

void label_adj(struct image *ip, struct image *ir, int q1, int q2)

{

      if(ir->pval == -1) {

         switch(q1) {
         case 2:
           label_adj(ip, ir->ne, q1, q2);
           break;
         case 3:
           label_adj(ip, ir->sw, q1, q2);
           break;
         case 0:
           break;
         }

         switch(q2){
         case 4:
           label_adj(ip, ir->se, q1, q2);
           break;
         case 0:
           break;
         }
  
      }

      else if(ir->pval == 1) 

          ip->label = unify(find(ip->label), find(ir->label));

      return;

}
