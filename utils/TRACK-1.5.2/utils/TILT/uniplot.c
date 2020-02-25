#include <stdio.h>
#include <strings.h>
#include "uniras.h"

#define MXLINE 2000

int main(int argc, char **argv)
{

   int i, j, k;
   int iseg;
   int ntr=0, nf=0;
   int trnpt=0;
   int dum;
   int legn[100];

   long int *tim=NULL;

   float xlim, ylim, xsize, ysize;
   float xmin=-10.0, xmax=10.0, ymin=900.0, ymax=200.0;
   float srat=0.8, xof=0.0, yof=0.0;

   float *lev=NULL, *levt=NULL;
   float *tilt=NULL;
   float color[3], hatch[5];

   FILE *fin=NULL;

   char trin[100];
   char device[50]="device=mx11;e";
   char line[MXLINE];
   char *ih=NULL;

   static char leg[100][20];

   color[1] = 50.0;
   color[2] = 50.0;

   printf("What is the track input file?\n\n");
   scanf("%s", trin);

   printf("Specify device for UNIRAS use.\n\n");
   scanf("%s", device);

   iseg = 1;
   finit_c(argc, argv);
   groute_c(device);
   gopen_c();
   gsegwk_c(0);
   glimit_c(xmin, xmax, ymin, ymax, 0.0, 0.0);
   grpsiz_c(&xlim, &ylim);
   xsize = 0.8 * ((xlim < ylim) ? xlim : ylim);
   ysize = srat * xsize;
   xof = 0.4 * (xlim - xsize);
   yof = 0.7 * (ylim - ysize);
   gvport_c(xof, yof, xsize, ysize);
   gwbox_c(xsize, ysize, 0.0);
   gscale_c();
   
   fin = fopen(trin, "r");
   if(!fin){
      printf("****ERROR****, can't open file %s for read.\n\n", trin);
      exit(1);
   }

   fgets(line, MXLINE, fin);
   sscanf(line, "%*s %d %*s %d", &ntr, &nf);
   
   lev = (float *)calloc(nf, sizeof(float));
   if(!lev){
      printf("****ERROR****, cannot allocate memory for level info.\n\n");
      exit(1);
   }

   levt = (float *)calloc(nf, sizeof(float));
   if(!levt){
      printf("****ERROR****, cannot allocate memory for level info.\n\n");
      exit(1);
   }

   tilt = (float *)calloc(nf, sizeof(float));
   if(!tilt){
      printf("****ERROR****, cannot allocate memory for level info.\n\n");
      exit(1);
   }

   fgets(line, MXLINE, fin);
   ih = strtok(line, " ");
   sscanf(ih, "%*s");
   ih = strtok(NULL, " ");
   for(i=0; i < nf; i++) {
       sscanf(ih, "%f", lev+i);
       ih = strtok(NULL, " ");
   }

   for(i=0; i < ntr; i++){
 
       fgets(line, MXLINE, fin);
       sscanf(line, "%*s %*d %*s %d", &trnpt);

       tim = (long int *)calloc(trnpt, sizeof(long int));
       if(!tim){
          printf("****ERROR****, cannot allocate memory for level info.\n\n");
          exit(1);
       }

       legn[0] = 0;

       gsegcr_c(i+1);
       raxis_c(1, 900.0, 3., 1);
       ruxlop_c(2);
       raxis_c(2, 0.0, 3., 1);

       for(j=0; j < trnpt; j++){
           color[0] = j * (int)(350 / trnpt);
           rcolor_c(0, 14+j, color, hatch);
           memcpy(levt, lev, nf*sizeof(float));
           fgets(line, MXLINE, fin);
           ih = strtok(line, " ");
           sscanf(ih, "%ld", tim+j);
           sprintf(&leg[j+1][0], "%ld", *(tim+j));
           legn[j +1] = 10;
           ih = strtok(NULL, " ");
           for(k=0; k < nf; k++){
               sscanf(ih, "%e", tilt+k);
               if(*(tilt+k) > 1.0e+10) {*(tilt+k) = 999.999; *(levt+k) = 999.999;}
               ih = strtok(NULL, " ");
           }

           blcol_c(14+j);
           blwid_c(1.0);
           bdimx_c(1, tilt, nf);
           blivec_c(levt, 7, 0, 0);

/*gempty_c();
printf("input dummy\n");
scanf("%*d"); */

       }

       bpleg_c(9., 500., legn, 20, leg, trnpt+1, 3.);

       gempty_c();
       gsegcl_c(i+1);
       gclear_c();


       free(tim); tim=NULL;

   }

   fclose(fin);

   gclose_c();

   free(lev);

   free(tilt);

   return 0;
}
