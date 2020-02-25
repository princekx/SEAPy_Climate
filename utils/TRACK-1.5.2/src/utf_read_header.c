#include <Stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grid.h"

#define  MAXCHAR  100


void utf3_read_header(float ,GRID * , int * , FILE * ,int * , int * );
void utf4_read_header(float ,GRID * , int * , FILE * ,int * , int * );

int utfv = 3;
int eqsw = 0;
int whemi = 0;

void utf_read_header(float ival, GRID *gr, int *frnum, FILE *fdatin, int *tl, int *gof)

{


    char flin[MAXCHAR];

    fgets(flin, MAXCHAR, fdatin);

    if(!strstr(flin, "UTF")){
       printf("***ERROR***, the input data file is not a UTF file, aborting program\n");
       exit(1);
    }

    if(strstr(flin, "1.4")) utfv = 4;

    switch(utfv){
        case 3:
          utf3_read_header(ival, gr, frnum, fdatin, tl, gof);
          break;
        case 4:
          utf4_read_header(ival, gr, frnum, fdatin, tl, gof);
          break;
        default:
          printf("***ERROR***, icorrect specifier for utf header\n");
    }

    return;

}
