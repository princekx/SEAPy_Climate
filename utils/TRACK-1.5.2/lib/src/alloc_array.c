#include <stdlib.h>
#include "mem_er.h"

void *alloc_array(int size)

{

   void *ptr=NULL;

   ptr = (void *)malloc(size);

   mem_er((ptr == NULL) ? 0 : 1, size);

   memset(p, NULL, size);

   return ptr;

}
