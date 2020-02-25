#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "file_cat_out.h"
#include "mem_er.h"

/* function to reformat tracks into IBTrACS format */

struct tot_tr *read_tracks(FILE *, int *, int *, int *, int , float *, float * , float ** , int * );
void write_track_netcdf(struct tot_tr * , int , char * , int , char ** , int );


int noheader=0;

extern int nfld, nff;
extern int *nfwpos;

int main(int argc, char *argv[])
{
   int i=0;
   int trnum=0;
   int gpr=0, ipr=0;
   int itrtyp='s';
   int nmeta=0;
   
   float alat=0.0, alng=0.0;

   char com[]="Usage: tr2nc trfilein itrtype(s or v) metafile";
   
   char *meta=NULL;
   char **metad=NULL;
   char line[MAXCHR];
   
   FILE *fin=NULL;
   FILE *fmeta=NULL;
   
   struct tot_tr *tracks=NULL;

   if(argc < 3){

     printf("****ERROR***, %s \n", com);
     exit(1);
   }  

   sscanf(argv[2], "%d", &itrtyp);
   if(argc == 4) {
      meta = (char *)calloc(MAXCHR, sizeof(char));
      mem_er((meta == NULL) ? 0 : 1, MAXCHR*sizeof(char));
      sscanf(argv[3], "%s", meta);
      
/* read attribute data */
 
      fmeta = fopen(meta, "r");
      if(!fmeta){
          printf("****ERROR****, can't open file %s\n", meta);
          exit(1);
      }
      fgets(line, MAXCHR, fmeta);
      sscanf(line, "%d", &nmeta);
      metad = (char **)calloc(nmeta, sizeof(char *));
      mem_er((metad == NULL) ? 0 : 1, nmeta*sizeof(char *));
      for(i=0; i < nmeta; i++){
           metad[i] = calloc(MAXCHR, sizeof(char));
	   mem_er((metad[i] == NULL) ? 0 : 1, MAXCHR*sizeof(char));
	   fgets(line, MAXCHR, fmeta);
	   strncpy(metad[i], line, MAXCHR);
      }
       
      fclose(fmeta);    
        
   }

   fin = fopen(argv[1], "r");
   if(!fin){
      printf("***ERROR***, unable to open file %s for 'r'\n\n", argv[1]);
      exit(1);
   }

   tracks = read_tracks(fin, &trnum, &gpr, &ipr, 's', &alat, &alng, NULL, NULL);
   
   fclose(fin);   
   
   write_track_netcdf(tracks, trnum, argv[1], itrtyp, metad, nmeta);

   return 0;
}
