#include <Stdio.h>
#include <stdlib.h>

void close_file(FILE *ptr, char *filenm)

{

    if(fclose(ptr) != 0){

       printf("***error*** closing file %s \n", filenm);
       perror(filenm);
       exit(1);

     }

     ptr = NULL;

     return;

}
