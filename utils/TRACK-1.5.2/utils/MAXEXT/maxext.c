#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "splice.h"
#include "m_values.h"
#include "mem_er.h"

#define  LARGE   1.0e+20
#define  SMALL  -1.0e+20

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int *
);
void meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );
void *realloc_n(void *ptr, int size);
int maxmin(struct tot_tr * , int , int , int );

int noheader=0;

extern int nfld, nff;
extern int *nfwpos;

int main(int argc, char *argv[])
{

    int i=0, j=0;

    int gpr=0, ipr=0;
    int trnum=0;
    int ift=0, ifft=0;
    int ifrid=0;
    int inoin=0;
    int ilif=0;

    float alng=0.0, alat=0.0;
    float str=0.0;
    float *addtmp=NULL;

    struct tot_tr *all_tr=NULL, *atr=NULL;
    struct fet_pt_tr *fp1=NULL, *fp2=NULL; 

    FILE *fin=NULL;
    FILE *fout=NULL;

    char filin[MAXCHR];
    char filout[MAXCHR];

    if(argc != 4){
       printf("Usage: maxext [infile] [additional field] [lifecycle stage: 0 (genesis), 1 (lysis), 2 (max), 3 (min)]\n\n");
       inoin=1; 
    }

    if(inoin){
       printf("What is the input track file?\n\n");
       scanf("%s", filin);
       strcpy(filout, filin);
       strcat(filout, ".new"); 
    }
    else{
      strcpy(filout, argv[1]);
      strcat(filout, ".new");
      sscanf(argv[1], "%s", filin);
      sscanf(argv[2], "%d", &ift); 
      sscanf(argv[3], "%d", &ilif);
    }

    fin = fopen(filin, "r");
    if(!fin){

       printf("****ERROR****, no such file exists, %s\n\n", argv[1]);
       exit(1);


    }

    all_tr = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    fclose(fin);

    if(nff){
       printf("****INFORMATION****, there are %d additional fields associated with these tracks \r\n", nff);

       if(inoin){
          printf("Which additional field is required to extract lifecycle stage, input '0' for default?\n\n");
          scanf("%d", &ift);
       }
       if(ift < 0 || ift > nff){
          printf("****ERROR****, not a valid field Id.\n\n");
          exit(1);
       }
       ifft = 0;
       if(ift){
          for(i=0; i < ift; i++){
              if(*(nfwpos + i)) ifft += 3;
              else ifft += 1;
          }
          --ifft;
       } 

    }

    if(inoin){
       printf("What stage of lifecycle to test:- \r\n"
              "Input '0' for genesis             \r\n"
              "Input '1' for lysis               \r\n"
              "Input '2' for maximum intensity   \r\n"
              "Input '3' for minimum intensity   \n\n");
       scanf("%d", &ilif);
    }
    if(ilif < 0 || ilif > 3){
       printf("****ERROR****, incorrect option value input.\n\n");
       exit(1);
    }

    for(i=0; i < trnum; i++){
        atr = all_tr + i;
        fp1 = atr->trpt;

        if(atr->num == 0) continue;

        if(!ilif){
           if(ift)
             str = (atr->trpt)->add_fld[ifft];
           else
             str = (atr->trpt)->zf;

           if(str > ADD_CHECK) continue;

           ifrid = 0;
        }
        else if (ilif == 1){
           if(ift)
             str = (atr->trpt + atr->num - 1)->add_fld[ifft];
           else
             str = (atr->trpt + atr->num - 1)->zf;

           if(str > ADD_CHECK) continue;

           ifrid =  atr->num - 1;
        }
        else if(ilif > 1){
           ifrid = maxmin(atr, ilif, ift, ifft);
        }

        if(ifrid < 0) {
           if(nff){
              for(j=0; j < atr->num; j++) {
                 fp2 = atr->trpt + j;
                 free(fp2->add_fld);
              } 
           }
           free(atr->trpt);
           atr->num = 0;
        }
        else{

           if(ilif){
              addtmp = fp1->add_fld;
              fp2 = atr->trpt + ifrid;
              memcpy(fp1, fp2, sizeof(struct fet_pt_tr));
              fp1->add_fld = addtmp;

              if(nff){
                for(j=0; j < nfld; j++) fp1->add_fld[j] = fp2->add_fld[j];
              }
           }

           for(j=1; j < atr->num; j++) {
               fp2 = atr->trpt + j;
               free(fp2->add_fld);
           }

           atr->num = 1;

           atr->trpt = (struct fet_pt_tr *)realloc_n(atr->trpt, sizeof(struct fet_pt_tr));
           mem_er((atr->trpt == NULL) ? 0 : 1, sizeof(struct fet_pt_tr));  

        }


    }

    fout=fopen(filout, "w");
    if(!fin){
       printf("****ERROR****, can't open file %s\n", filout);
       exit(1);
    }

    meantrd(fout, all_tr , trnum, 's', gpr, ipr, alat, alng);

    fclose(fin);

    return 0;
}

int maxmin(struct tot_tr *atr, int ilif, int ift, int ifft)
{
    int i=0;
    int ifrd=-1;

    float str=0.0;
    float mxmn=0.0;

    struct fet_pt_tr *fp=NULL; 

    if(ilif == 2) mxmn = SMALL;
    else if(ilif == 3) mxmn = LARGE;
    else {
       printf("****ERROR****, incorect lifetime identifier.\n\n");
       exit(1);
    }

    for(i=0; i < atr->num; i++){
       fp = atr->trpt + i; 
       if(ift)
         str = fp->add_fld[ifft];
       else
         str = fp->zf;

       if(str > ADD_CHECK) continue;

       if(ilif == 2){
         if(str > mxmn) {mxmn = str; ifrd = i;}
       }
       else if(ilif == 3){
         if(str < mxmn) {mxmn = str; ifrd = i;}
       }
       else{
         printf("****ERROR****, option unknown.\n\n");
         exit(1);
       }

    }

    return ifrd;
}
