#include <stdlib.h>
#include <string.h>

/* function to call malloc and initialize the memory allocated */

void *malloc_initl(int size)

{

    void *p;

    p = (void *)malloc(size);

    memset(p, 0, size);

    return p;

}
