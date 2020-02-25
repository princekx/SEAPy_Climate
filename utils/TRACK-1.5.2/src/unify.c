#include <Stdio.h>
#include <stdlib.h>
#include "st_im.h"

/* funtion to merge labels for the hierarchical data tree */

struct eq_class *unify(struct eq_class *lab1, struct eq_class *lab2)

{

     if(lab1 != NULL && lab1 != lab2)

        lab1->farther = lab2;

     return lab2;

}
