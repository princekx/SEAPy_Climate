#include <Stdio.h>
#include <stdlib.h>
#include "mem_er.h"

void *realloc_n(void *ptr, int size)

{

      if(size < 0){

         printf("****ERROR****, in %s, memory block size  must be >=0\n\n", __FILE__);
         exit(1);

      }

      if(!ptr && !size) return ptr;

      if(!size) {free(ptr); ptr = NULL;}

      else
        ptr = (void *) realloc(ptr, size);

      return ptr;

}
