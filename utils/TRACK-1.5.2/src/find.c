#include <Stdio.h>
#include <stdlib.h>
#include "st_im.h"

/* function to manipulate equivelence classes for hierarchical data tree */

struct eq_class *find(struct eq_class *lab)

{

     struct eq_class *r, *temp;

     if(lab->farther == NULL ) return lab;
 
     else {

          r = lab;
          
          while(r->farther != NULL) r=r->farther;

          while(lab->farther != NULL) {

                temp = lab->farther;
                lab->farther = r;
                lab = temp;

           }

           return r;

      }

}
