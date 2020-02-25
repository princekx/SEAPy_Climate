#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zones.h"
#include "files_in.h"
#include "file_handle.h"
#include "mem_er.h"
#include "m_values.h"

#define  MAXCHR  200

int wc(char * );

extern int tom;

ADPT *read_adptp(int parn)

{

    int i;
    
    FILE  *fad=NULL;

    char cbuf[MAXCHR];
    char par_file[MAXCHR];
    char num[4];

    ADPT *add=NULL;

    sprintf(num, "%d", parn);
    strncpy(par_file, DATAD, MAXCHR);
    strcat(par_file, num);
    
    fad = open_file(par_file, "r");

    printf("****INFORMATION****, reading adaptive phimax values from file:\r\n"
           "                     %s\n\n", par_file);

    add = (ADPT *)malloc_initl(sizeof(ADPT));
    mem_er((add == NULL) ? 0 : 1, sizeof(ADPT));

    fgets(cbuf, MAXCHR, fad);

    if(wc(cbuf) != 4){

       printf("****ERROR****, insuficient input for adaptive tracking \r\n"
              "               check input file.                       \n\n");

       exit(1);

    }

    sscanf(cbuf, "%f %f %f %f", &(add->ct[0]), &(add->ct[1]), &(add->ct[2]), &(add->ct[3]));

    if(add->ct[0] > add->ct[1] || add->ct[1] > add->ct[2] || add->ct[2] > add->ct[3]){

       printf("****ERROR****, distance cutt-offs for adaptive phimax must be incremental\n\n");
       exit(1);

    }

    if(tom == 'g'){

       add->ct[0] *= FP_PI;
       add->ct[1] *= FP_PI;
       add->ct[2] *= FP_PI;
       add->ct[3] *= FP_PI;

    }

    fgets(cbuf, MAXCHR, fad);

    if(wc(cbuf) != 4){

       printf("****ERROR****, insuficient input for adaptive tracking \r\n"
              "               check input file.                       \n\n");

       exit(1);

    }

    sscanf(cbuf, "%f %f %f %f", &(add->phii[0]), &(add->phii[1]), &(add->phii[2]), &(add->phii[3]));

    close_file(fad, par_file);

    add->maxad = 0.0;
    
    for(i=0; i < 4; i++){
        if(add->phii[i] > add->maxad) add->maxad = add->phii[i];
    }


    add->psl[0] = (add->phii[1] - add->phii[0]) / (add->ct[1] - add->ct[0]);
    add->incp[0] = (add->ct[1] * add->phii[0] - add->ct[0] * add->phii[1]) / (add->ct[1] - add->ct[0]);

    add->psl[1] = (add->phii[2] - add->phii[1]) / (add->ct[2] - add->ct[1]);
    add->incp[1] = (add->ct[2] * add->phii[1] - add->ct[1] * add->phii[2]) / (add->ct[2] - add->ct[1]); 
    
    add->psl[2] = (add->phii[3] - add->phii[2]) / (add->ct[3] - add->ct[2]);
    add->incp[2] = (add->ct[3] * add->phii[2] - add->ct[2] * add->phii[3]) / (add->ct[3] - add->ct[2]); 

    return add;    

}


int wc(char *st)

{

    int nw=0, state=0;
    int i=0, c;

    while((c=st[i++]) != '\0'){

       if(c == ' ' || c == '\n' || c == '\t')

          state = 0;

       else if (!state){

          state = 1;
          ++nw;
       }


    }

    return nw;

}      
