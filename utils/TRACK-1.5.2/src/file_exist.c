#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "file_handle.h"


char *file_exist(char *filenm, char *mode, char *apps)

{

   FILE *ptr;

   while((ptr = fopen(filenm, mode))){

     printf("***WARNING*** file %s exists.\n", filenm);
     filenm = strcat(filenm, apps);
     printf("              new file name is %s \n\n", filenm);
     close_file(ptr, filenm);


   }

   return filenm;

}
