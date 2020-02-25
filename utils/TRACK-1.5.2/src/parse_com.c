#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "file_handle.h"
#include "file_path.h"

#define  NWFCHR     100

#define  FPATHI     Add(USER,PATHI,)

/* function to allow command line changes to input file and file extensions */

int parse_com(int argc, char **argv, char *filnam, char *fext, int maxchr, int mext)
{

    int narg = 1;
    int iext=0;

    char usage[] = "Usage: track -i <input file> -f <file extension> -u <URL file> -w <URL>";

    char wfile[NWFCHR];

    FILE *fw=NULL;


    while(narg < argc){

         if(!strncmp(argv[narg], "-i", 2)){
            strcpy(filnam, FPATHI);
            if(strlen(argv[narg]) == 2) {
               ++narg;
               if(!strncmp(argv[narg], "-", 1)){
                 printf("****ERROR****, incorrect command line input.\n\n");
                 printf("\n %s\n\n", usage);
                 exit(1);
               }

               if(strlen(argv[narg]) > maxchr) {
                  printf("****ERROR****, command line input to long for character array\n\n");
                  exit(1);
               }
               strcat(filnam, argv[narg]);

             }
             else {
               if(strlen(argv[narg]) > maxchr) {
                  printf("****ERROR****, command line input to long for character array\n\n");
                  exit(1);
               }
               strcat(filnam, argv[narg]+2);
             }
             ++narg;
         }

         else if(!strncmp(argv[narg], "-f", 2)){
            if(strlen(argv[narg]) == 2) {
               ++narg;
               if(!strncmp(argv[narg], "-", 1)){
                 printf("****ERROR****, incorrect command line input.\n\n");
                 printf("\n %s\n\n", usage);
                 exit(1);
               }

               if(strlen(argv[narg]) > mext) {
                  printf("****ERROR****, command line input to long for character array\n\n");
                  exit(1);
               }               
               sscanf(argv[narg], "%s", fext);

             }
             else {
               if(strlen(argv[narg]) > mext) {
                  printf("****ERROR****, command line input to long for character array\n\n");
                  exit(1);
               }
               sscanf(argv[narg]+2, "%s", fext);
             }

             ++narg;
             iext = 1;
         }

         else if (!strncmp(argv[narg], "-u", 2)){
             if(strlen(argv[narg]) == 2) {
               ++narg;
               if(!strncmp(argv[narg], "-", 1)){
                 printf("****ERROR****, incorrect command line input.\n\n");
                 printf("\n %s\n\n", usage);
                 exit(1);
               }

               if(strlen(argv[narg]) > maxchr) {
                  printf("****ERROR****, command line input to long for character array\n\n");
                  exit(1);
               }               
               sscanf(argv[narg], "%s", wfile);
               fw = open_file(wfile, "r");
                   fscanf(fw, "%s", filnam);
               close_file(fw, wfile);

             }
             else {
               if(strlen(argv[narg]) > maxchr) {
                  printf("****ERROR****, command line input to long for character array\n\n");
                  exit(1);
               }
               sscanf(argv[narg]+3, "%s", filnam);
             }

             ++narg;
         }


         else if(!strncmp(argv[narg], "-w", 2)){
            if(strlen(argv[narg]) == 2) {
               ++narg;
               if(!strncmp(argv[narg], "-", 1)){
                 printf("****ERROR****, incorrect command line input.\n\n");
                 printf("\n %s\n\n", usage);
                 exit(1);
               }

               if(strlen(argv[narg]) > maxchr) {
                  printf("****ERROR****, command line input to long for character array\n\n");
                  exit(1);
               }               
               sscanf(argv[narg], "%s", filnam);

             }

             else {
               if(strlen(argv[narg]) > maxchr) {
                  printf("****ERROR****, command line input to long for character array\n\n");
                  exit(1);
               }
               sscanf(argv[narg]+2, "%s", filnam);
             }

             ++narg;
         }

         else {
           printf("****ERROR****, incorrect command line input.\n\n");
           printf("\n %s\n\n", usage);
           exit(1);
         }

    }

    return iext;

}
