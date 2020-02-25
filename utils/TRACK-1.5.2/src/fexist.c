#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(HP_UX)
#define _INCLUDE_POSIX_SOURCE
#endif
#include <dirent.h>
#include "file_handle.h"


int fexist(char *filenm, char *mode)

{

   FILE *ptr;
   DIR *dptr;

   if((ptr = fopen(filenm, mode))){

     close_file(ptr, filenm);

     if((dptr = opendir(filenm))) {

        printf("***ERROR***, filename is a directory\n\n");

        closedir(dptr);
        return 0;

     }

     return 1;


   }

   return 0;

}
