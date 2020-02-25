#include <Stdio.h>
#include <stdlib.h>
#include "file_handle.h"

FILE *open_file(char *filenm, char *mode)

{

   FILE *ptr;

   if(!(ptr = fopen(filenm, mode))){

     printf("***error*** opening file %s for %s \n", filenm, mode);
     perror(filenm);
     exit(1);

   }

   printf("\nFile %s \r\n opened for '%s'\n\n", filenm, mode);

   return ptr;

}
