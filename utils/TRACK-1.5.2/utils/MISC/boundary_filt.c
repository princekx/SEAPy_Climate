#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "splice.h"
#include "mem_er.h"
#include "file_handle.h"

/* combine sets of tr_trs files */

struct tot_tr *read_tracks(FILE * , int * , int * , int * , int , float * , float * , float ** , int *);

int meantrd(FILE * , struct tot_tr * , int , int , int , int , float , float );

int noheader=0;

extern float sum_per;
extern int iper_num;
extern int aniso;
extern int nff, nfld;
extern int *nfwpos;


int main(int argc, char *argv[])

{
    int i, j;
    int np=0, ip=0, num=0, imiss=0;
    int gpr, ipr;
    int tr_count=0, trnum=0;
    int nf=0;
    int iadd=0, iff=0;

    float alat=0.0, alng=0.0;

    float x1in, x2in, y1in, y2in;
    float x1ou, x2ou, y1ou, y2ou;
    
    float xpf=0.0, ypf=0.0;

    char filin[MAXCHR];
    char filout[MAXCHR];
    char usage[MAXCHR]="bfilt [infile] [NFILT]";

    FILE *fin=NULL;
    FILE *fout=NULL;

    struct tot_tr *alltr=NULL, *atr=NULL, *att=NULL;
    struct fet_pt_tr *fp=NULL, *fptt=NULL;

    if(argc != 3){
      printf("%s\n", usage);
      exit(1);
    }

    fin = fopen("region.dat", "r");
    if(!fin){
       printf("****ERROR****, can't open file %s for read.\n\n", "region.dat");
       exit(1);
    }

    fscanf(fin, "%f %f %f %f", &x1in, &x2in, &y1in, &y2in);
    fscanf(fin, "%f %f %f %f", &x1ou, &x2ou, &y1ou, &y2ou);

    fclose(fin);

    sscanf(argv[1], "%s", filin);
    sscanf(argv[2], "%d", &nf);

    strncpy(filout, filin, MAXCHR);
    strcat(filout, ".new");

    fin = fopen(filin, "r");
    if(!fin){
       printf("****ERROR****, can't open file %s for read.\n\n", filin);
       exit(1);
    } 


    alltr = read_tracks(fin, &tr_count, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);

    fclose(fin);
    
    if(nff){
       printf("Which additional field do you want to use, input '0' for default.\n\n");
       scanf("%d", &iadd);
       
       if(iadd < 0 || iadd > nff){
          printf("****ERROR****, addtional field identifier not valid.\n\n");
	  exit(1);
       }
       
 
       
       if(iadd){
       
          if(! *(nfwpos + iadd - 1)){
             printf("****ERROR****, additional field does not have positional information.\n\n");
	     exit(1);
          }
	         
          iff = 0;
          for(i=0; i < iadd; i++){
              if(*(nfwpos + i)) iff += 3;
              else iff += 1;
          }
          --iff;
       }
    }

    trnum = tr_count;

    for(i=0; i < tr_count; i++){

        atr = alltr + i;
        imiss = 0;
        for(j=0; j < atr->num; j++){
            fp = atr->trpt + j;
	    
	    if(iadd){
	       xpf = fp->add_fld[iff-2];
	       ypf = fp->add_fld[iff-1];
	    }
	    else {
	       xpf = fp->xf;
	       ypf = fp->yf;
	    }
	   

            if((xpf < x1in || xpf > x2in) ||
               (ypf < y1in || ypf > y2in)    ) {fp->fr_id = -1; imiss = 1;}
        }

        if(imiss){

           ip = 0;
           num = atr->num;
           fptt = atr->trpt;

           while(ip < num){

              fp = fptt;
              np = 0;
              while(ip < num && fp->fr_id < 0){ 
                   ++np, ++ip;
                   ++fp;
              }
              
              fptt = fp;
              atr->num -= np;

              if(ip == num) break;

              np = 0;
              fp = fptt;

              while(ip < num && fp->fr_id > 0) {
                    ++ip; ++np;
                    ++fp;
              }

              atr->num -= np;

              if(np > 0){

                 ++trnum;

                 alltr = (struct tot_tr *)realloc_n(alltr, trnum * sizeof(struct tot_tr));
                 mem_er((alltr == NULL) ? 0 : 1, trnum * sizeof(struct tot_tr));

                 att = alltr + trnum - 1;
                 att->num = np;
                 att->trpt = NULL;
                 att->awt = 0;
                 att->time = 0;

                 atr = alltr + i;

                 att->trid = atr->trid;

                 att->trpt = (struct fet_pt_tr *)calloc(np, sizeof(struct fet_pt_tr));
                 mem_er((att->trpt == NULL) ? 0 : 1, np*sizeof(struct fet_pt_tr));

                 memcpy(att->trpt, fptt, np * sizeof(struct fet_pt_tr));

              }

              fptt = fp;

              if(ip == num) break;

           }

           if(! atr->num){free(atr->trpt); atr->trpt = NULL;}

        }

    }

    for(i=0; i < trnum; i++){
        atr = alltr + i;
        if(atr->num < nf){free(atr->trpt); atr->num = 0;}
    }

    fout = fopen(filout, "w");
    if(!fout){
       printf("****ERROR****, can't open file %s for write.\n\n", filout);
       exit(1);
    }

    meantrd(fout, alltr, trnum, 's', gpr, ipr, alat, alng);

    fclose(fout);

    return 0;
}
