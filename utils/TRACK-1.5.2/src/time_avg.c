#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Math.h>
#include <sys/types.h>
#include <stdint.h>
#include "grid.h"
#include "files_out.h"
#include "file_handle.h"
#include "mem_er.h"
#include "utf.h"
#include "pp.h"
#include "m_values.h"
#include "periodogram.h"
#include "region.h"
#include "netcdf_info.h"



#define  LARGE_FREQ   1.0e+6
#define  TOLVAR       1.0e-10
#define  TOLWT        1.0e-6
#define  DEFVAL       0.0
#define  TOLSCL       1.0e-20


/* function to compute time average */

extern GRID *gr;
extern int frnum;
extern int form;
extern int eqsw;
extern int utfv;
extern int i_utf4x, i_utf4y;
extern int nw_ln;
extern int std_x, std_y;
extern int tom;

extern float *abuf;
extern float period;

extern PP *pph;

extern char *fext;

extern int iext;

float ssc;

extern char *chrfld;
int nchr=0;

char *ihead[NUMLINES];
int nlines=0;

int rf(FILE * ,FILE * , int , int );
void writef(FILE * , float * , int );
float *field_varience(FILE * , float * , int , int , int , int , int , int , off_t , float * , float , float * , int );

void var_filt_freq(float * , float * , FILE * , int , int , int , int , int , off_t , float * , float , int );

int sub_tim_avg(FILE * , float * ,int , int , int , int , off_t , float * , float , int );
int missing(float, float, int);

int time_avg(FILE *fdat, int fr1, int fri, int frl, int nvar)

{

   int i, nf=0;
   int nn=0;
   int dim=0;
   int numfr=0;
   int numwt=0;
   int imiss='n';
   int icmp=0;

   int wf='a';

   float *ata=NULL, *avar=NULL;
   float *atmp=NULL;
   float mval=MVAL, wt;
   float *wght=NULL;
   float *anum=NULL, *anumn=NULL;

   double *acmp=NULL, *acmpn=NULL;

   off_t chrnum=0, place1=0, place2=0;

   FILE *ftav=NULL;
   FILE *fwght=NULL;

   char fintav[MAXCHR];
   char foutav[MAXCHR];
   char finw[MAXCHR];
   char fsumw[MAXCHR];

   strncpy(fsumw, TAVGE, MAXCHR);
   if(iext) strcpy(strstr(fsumw, EXTENSION), fext);
   strcat(fsumw, "_wsum");

/* assign space for header strings */

   nlines = 0;

   for(i=0; i < NUMLINES; i++){
       ihead[i] = (char *)calloc(MAXCHR, sizeof(char));
       mem_er((ihead[i] == NULL) ? 0 : 1, MAXCHR * sizeof(char));
   }


   if(form < 2 || form == 4) dim = std_x * std_y;
   else if (form == 2){
      if(utfv == 3) dim=gr->ix * ((eqsw) ? gr->iy + 1 : gr->iy);
      else if(utfv == 4) dim = i_utf4x * i_utf4y;

   }

   else if(form == 3) dim = pph->lbrow * pph->lbnpt;

   printf("Do you want all frames in data file to contribute to the   \r\n"
          "time average or time average removal, or just those frames \r\n"
          " selected or all frames at the specified sampling rate from\r\n"
          " the starting frame, 'a', 's', or 'i'.                     \n\n");
 

   scanf("\n");
   wf = getchar();

   if(wf == 'a' ){fr1 = 1; fri = 1; frl = frnum;}

   else if(wf == 's'){
      if(fri < 1) fri = 1;
      if(frl > frnum) frl = frnum;
   }

   else if(wf == 'i'){
      if(fri < 1) fri = 1;
      frl = frnum;
   }

   else {

      printf("***WARNING***, incorrect decriptor for time average calculation\r\n"
             "               defaulting to 'a' (all).                        \n\n");

      wf = 'a';
      fr1 = 1; fri = 1; frl = frnum;

   }

   printf("***INFORMATION***, frame counters for computing time average or time average removal are:\r\n"
          "                   Start = %d, Interval = %d, Last = %d \n\n", fr1, fri, frl);

   ata = (float *)calloc(dim, sizeof(float));
   mem_er((ata == NULL) ? 0 : 1, dim * sizeof(float));

   if(!nvar){

      if(form == 4) {

        if(((NETCDF_INFO *)fdat)->imiss){
           printf("****INFORMATION****, missing data value is %e\n\n", ((NETCDF_INFO *)fdat)->missing);
           mval = ((NETCDF_INFO *)fdat)->missing;
           imiss = 'y';

        }
 
      }

      else {

         printf("Are there missing data values, delineated by some extreme value. 'y' or 'n'\n\n");

         scanf("\n");
         imiss = getchar();

      }

      if(imiss == 'y'){

         printf("What is the missing data value?\n\n");
         scanf("%f", &mval);

         printf("How do you want to compare with missing value?\r\n"
                "'0' -- equality.                              \r\n"
                "'1' -- less than.                             \r\n"
                "'2' -- greater than.                          \n\n");
         scanf("%d", &icmp);

         anum = (float *)calloc(dim, sizeof(float));
         mem_er((anum == NULL) ? 0 : 1, dim * sizeof(float));

      }

      else mval = MVAL;

      printf("Do you want to load and subtract an existing time average field\r\n"
             "and output the filtered fields, 'y' or 'n'.                    \n\n");

      scanf("\n");

      if(getchar() == 'y'){

         printf("****WARNING****, any time avarage file you load must be on the \r\n"
                "                 same grid as the data you want to filter. This\r\n"
                "                 is the users responsibility.                  \n\n");

         printf("What time average file do you want to load?\n\n");

         scanf("%s", fintav);

         if(form != 4) ftav = open_file(fintav, "r");
         else {
            ftav = (FILE *)nc_clone((NETCDF_INFO *)fdat, fintav, NC_OPEN_MODE);
            ((NETCDF_INFO *)ftav)->iframe = 0;
            ((NETCDF_INFO *)ftav)->nframe = 0;
         }


         atmp = abuf;
         abuf = ata;

         rf(ftav, NULL, dim, 1);

         abuf = atmp;

         if(form != 4) place1 = ftello(fdat);

         nn = sub_tim_avg(fdat, ata, fr1, fri, frl, dim, place1, anum, mval, icmp);
      
         if(form != 4) close_file(ftav, fintav);
         else netcdf_close((NETCDF_INFO *)ftav);

         free(ata);
         free(anum);

         return nn;

      }

   }

   acmp = (double *)calloc(dim, sizeof(double));
   mem_er((acmp == NULL) ? 0 : 1, dim * sizeof(double));

/* get weights for weighted statistics */

   printf("Do you want to computed weighted statistics, e.g. means,\r\n"
          "and variance 'y' or 'n'.                                \n\n");

   scanf("\n");
   if(getchar() == 'y'){

      numfr = ((frl - fr1) / fri) + 1;

      printf("What is the file containing the weights for each time step?\n\n");
      scanf("%s", finw);

      wght = (float *)calloc(numfr, sizeof(float));
      mem_er((wght == NULL) ? 0 : 1, numfr * sizeof(float));

      fwght = open_file(finw, "r");

      numwt=0;
      while(fscanf(fwght, "%f", wght+numwt) > 0) {
         ++numwt;
         if(numwt == numfr) break;
      }

      if(numwt != numfr){

         printf("****WARNING****, the file                            \r\n"
                "               %s\r\n"
                "               may contain insufficient values for the\r\n"
                "               number of time steps.               \n\n", finw);   

      }

      close_file(fwght, finw);

      if(!anum){
         anum = (float *)calloc(dim, sizeof(float));
         mem_er((anum == NULL) ? 0 : 1, dim * sizeof(float));

         mval = MVAL;

      }

      acmpn = (double *)calloc(dim, sizeof(double));
      mem_er((acmpn == NULL) ? 0 : 1, dim * sizeof(double));

      anumn = (float *)calloc(dim, sizeof(float));
      mem_er((anumn == NULL) ? 0 : 1, dim * sizeof(float));


   }

   printf("***INFORMATON***, Computing Time Average of Data....\r\n\n"
          "                                                    \r\n\n"
          " Please Wait, this may take some time ..............  \n\n");


/* initialize time average */

   for(i=0; i<dim; i++) *(acmp + i) = 0.0;

   if(anum){

      for(i=0; i<dim; i++) *(anum + i) = 0.0;

   }


/* open file for writing the time average */

   strncpy(foutav, TAVGE, MAXCHR);
   if(iext) strcpy(strstr(foutav, EXTENSION), fext);

   if(form != 4) ftav = open_file(foutav, "w");
   else {
     ftav = (FILE *)nc_define((NETCDF_INFO *)fdat, foutav);
     ((NETCDF_INFO *)ftav)->nframe = 0;
     ((NETCDF_INFO *)ftav)->iframe = 0;
   }

/* read in field data */

   while(nf <= frl){

        if(form != 4){

           if(!nf) {

              place1 = ftello(fdat);

              nlines = 0;

              rf(fdat, NULL, dim, 1);

              for(i=0; i<nlines; i++) fprintf(ftav, "%s", ihead[i]);

              if(fr1 == 1) ++nf;

              place2 = ftello(fdat);
              chrnum = place2 - place1;

              if(fr1 > 1){

                 fseeko(ftav, (off_t)0, FSTART);

                 fseeko(fdat, (fr1-2)*chrnum, ORIGIN);
                 nlines = 0;

                 rf(fdat, NULL, dim, 1);


                 for(i=0; i<nlines; i++) fprintf(ftav, "%s", ihead[i]);
   
                 nf = fr1;        

              }

              nf += fri;

           }

           else{

              fseeko(fdat, (fri-1)*chrnum, ORIGIN);

              if(rf(fdat, NULL, dim, 0)) break;

              nf += fri;

           }

        }

        else {

           if(!nf) nf = fr1;

           ((NETCDF_INFO *)fdat)->iframe = nf - 1;
           
           if(rf(fdat, NULL, dim, 1)) break; 
      
           nf += fri;

        }

        if(wght){

           if(nn+1 > numwt){
              printf("****ERROR****, insufficient weights for time series.\n\n");
              exit(1);

           }

           wt = *(wght + nn);


           for(i=0; i<dim; i++){


              if(!missing(*(abuf+i), mval, icmp)){

                *(acmpn + i) += *(abuf + i) * wt;
                *(acmp + i) += *(abuf + i);
                *(anum + i) += 1.0;
                *(anumn + i) += wt;

              }                    

           }

        }

        else if(anum){

           for(i=0; i<dim; i++){

              if(!missing(*(abuf+i), mval, icmp)){

                *(acmp + i) += *(abuf+i);
                *(anum + i) += 1;

              }  
                

           }

        }

        else{
           
           for(i=0; i<dim; i++) *(acmp + i) += *(abuf+i);

        }

        ++nn; 


        if(nf > frl) break;

   }

   printf("\n");

   printf("****INFORMATION****, number of frames for averaging is: %d\n\n", nn);


/* average the field */

   if(wght){

      for(i=0; i<dim; i++) {
         if(*(anumn+i) > TOLWT)
           *(ata + i) = *(acmpn + i) / *(anumn+i);
         else *(ata + i) = DEFVAL;
      }

      writef(ftav, ata, dim);

      for(i=0; i<dim; i++) {
          if(*(anum+i) > TOLWT)
            *(ata + i) = *(acmp + i) / *(anum+i);
          else *(ata + i) = DEFVAL;

      }



/* write sum of weights to file for use in combining statistics */

      if(form != 4) {
         fwght = open_file(fsumw, "w");
         for(i=0; i<nlines; i++) fprintf(fwght, "%s", ihead[i]);
      }
      else {
         fwght = (FILE *)nc_define((NETCDF_INFO *)fdat, fsumw);
         ((NETCDF_INFO *)fwght)->nframe = 0;
         ((NETCDF_INFO *)fwght)->iframe = 0;
      }

      writef(fwght, anumn, dim);

      if(form != 4) close_file(fwght, fsumw);
      else netcdf_close((NETCDF_INFO *)fwght);

      free(acmpn);
      free(anumn);

   }

   else if(anum){

      for(i=0; i<dim; i++) {
          if(*(anum+i) > TOLWT)
            *(ata + i) = *(acmp + i) / *(anum+i);
          else *(ata + i) = DEFVAL;
      }
      writef(ftav, ata, dim);

   }

   else{
      if(nn > 0){
         for(i=0; i<dim; i++) *(ata + i) = *(acmp + i) / (float) nn;
      }
      else{
         for(i=0; i<dim; i++) *(ata + i) = DEFVAL;
      }
      writef(ftav, ata, dim);

   }


   if(form != 4) close_file(ftav, foutav);
   else netcdf_close((NETCDF_INFO *)ftav);


   if(!nvar){

      printf("Do you want to output the field with time average removed?, 'y' or 'n'\n\n");
      scanf("\n");
      if(getchar() == 'y'){
         sub_tim_avg(fdat, ata, fr1, fri, frl, dim, place1, anum, mval, icmp); 

      }

      avar = field_varience(fdat, ata, fr1, fri, frl, dim, nn, wf, place1, anum, mval, wght, icmp);

      printf("Do you want to compute filtered varience fields, 'y' or 'n'?\n\n");

      scanf("\n");
      if(getchar() == 'y')var_filt_freq(ata, avar, fdat, fr1, fri, frl, dim, nn, place1, anum, mval, icmp);

   }

   if(chrfld) free(chrfld);

   free(ata);
   free(acmp);
   free(avar);
   free(anum);
   free(wght);

   for(i=0; i<NUMLINES; i++) free(ihead[i]);

   return nn;

}

void read_chr_field(FILE * , int , int , char * , int );
int pp_read(FILE * , int );
 

int rf(FILE *fdat, FILE *ftav, int dim, int hf)

{

   int i, k, ad=0;
   int br=0;
   int dum;
   int ic[CHMX];
   int isum;
   int ier=0;

   static int ixx=0, iyy=0;
   static int iwmis=0;

   float scl;
   float tmax, tmin;
   float arng = ARNG;
   float ffmx, ffmn, frng, fmax;

   char text[MAXCHR], txtmp[MAXCHR];
   char *ih=NULL;
   char *icp=NULL;
   char ii[CHMX];

   if(hf) nlines = 0;


   if(form == 0 || form == 1){

      fgets(text, MAXCHR, fdat);
      if(PTSTR > 0) printf("%s", text);

      if(hf){

/*        fprintf(ftav, "%s", text); */
        sprintf(ihead[nlines], "%s", text);
        ++nlines;

      }
      if(form == 0) {
         if(fread(abuf, sizeof(float), dim, fdat) != dim){
            printf("***WARNING***, problem reading field data in %s at %d, possible EOF.\n", __FILE__, __LINE__);
            return 1; 

         }

         if(nw_ln == 'y') fscanf(fdat, "%*c");


      }
      else {

         for(i=0; i<dim; i++)fscanf(fdat, "%e", (abuf+i));

         fgets(text, MAXCHR, fdat);

      }


   }

   else if(form == 2){

          for(i=0; i<3; i++){

             if(!fgets(text, MAXCHR, fdat)) {br = 1; break;}
             if(PTSTR > 0) printf("%s\n", text);
             if(strstr(text, EOD)) {br = 1; break;}
             if(hf){


                if(i < 2) {

/*                  fprintf(ftav, "%s", text); */

                  sprintf(ihead[nlines], "%s", text);

                  ++nlines;

                }
                else {

                   *(ihead[nlines]) = '\0';

                   ih = strtok(text, " ");

                   while(ih){

                       sscanf(ih, "%d", &dum);

                       ++ad;
                       switch(ad){
                          case 1:
                             ixx = dum;
                             break;
                          case 2:
                             iyy = dum;
                             break;
                          case 10:
                             nchr = dum;
                             if(nchr > CHMX){
                                printf("***ERROR***, to many charcters per number at decoding.\n");
                                exit(1);
                             }

                             else if(nchr <= 0){
                                printf("***ERROR***, number of characters for decoding UTF not set, what is the number of chracters per number?\n\n");
                                scanf("%d", &nchr);
                             }

                             break;
                       }

                       if(ad == 10) {
/*                         fprintf(ftav, "%d  ", nchr); */
                         sprintf(txtmp,"%d  ", nchr);
                         strcat(ihead[nlines], txtmp);
                       }
                       else {
/*                         fprintf(ftav, "%d  ", dum); */
                         sprintf(txtmp,"%d  ", dum);
                         strcat(ihead[nlines], txtmp);
                       }




                       if(ad == 15) break;

                       ih = strtok(NULL, " ");



                   }
                   ad = 0;
/*                   fprintf(ftav, "\n"); */
                   strcat(ihead[nlines], "\n");

                   ++nlines;
                }
 
             }


          }



          if(br) return 1;


          if(!chrfld){

             chrfld = (char *)calloc((ixx*iyy+1)*nchr, sizeof(char));
             mem_er((chrfld == NULL) ? 0 : 1, (ixx*iyy+1)*nchr * sizeof(char));

          }


          arng = pow((float)CNUM, (float)nchr) - 1.0;

          fscanf(fdat, "%e %e %e\n", &tmin, &tmax, &ssc);  

          ffmx = fabs(tmax);
          ffmn = fabs(tmin);

          fmax = (ffmx > ffmn) ? ffmx : ffmn;

          frng = tmax - tmin;

          if(frng <= fmax*FTOL)
             scl = 0.;
          else
             scl = frng/arng;        

          read_chr_field(fdat, ixx, iyy, chrfld, nchr);

          icp = chrfld;


          for(i=0; i<dim; i++){

            strncpy(ii, icp, nchr);
            icp += nchr;
            isum = 0;
            for(k=0; k < nchr; k++) {

               if(ii[k] == CMISS){
                  isum = -1; 
                  if(!iwmis){
                     printf("****WARNING****, missing data values have been encountered in this data.\n\n");
                     iwmis = 1;
                  }    
                  break;
               }

               isum *= CNUM;

               ic[k] = (int)(uintptr_t)strchr(lkup, ii[k])-(int)(uintptr_t)lkup; 
               if(ic[k] > CNUM || ic[k] < 0){

                 printf("***error***, illegal character in data field, character %c\n", ic[k]);
                 exit(1);

               }

               isum += ic[k];

            }

            if(isum < 0) *(abuf + i) = VMISS;

            else *(abuf + i) = (float)isum * scl + tmin;

         }

   }

   else if(form == 3) ier = (pp_read(fdat, -1)) ? 0 : 1;

   else if(form == 4) ier = netcdf_read_field((void *)fdat);

   if(!PTSTR) putchar('#');

   return ier;



}


void writef(FILE *ftav, float *ata, int dim)

{

   int i, j, ad=0;
   int ic1;
   int ifv;
   int bend;

   long int ss;

   float tmax, tmin;
   float scl, fscl;
   float fv;
   float arng=ARNG;

   char *cbuf=NULL;


/* write to file */

   if(form == 0){

      if(fwrite(ata, sizeof(float), dim, ftav) != dim){
         printf("***ERROR***, writing data, aborting\n");
         exit(1);
      }
      if(nw_ln == 'y') fprintf(ftav, "\n");

   }

   else if(form == 1){

      ad = 0;

      for(i=0; i<dim; i++){

          fprintf(ftav, "%12.5e ", *(ata + i));
          if((++ad) == 10){fprintf(ftav, "\n"); ad = 0;}

      }

      fprintf(ftav, "\n");

   }

   else if(form == 2){

/* compute max and min values of field */

       arng = pow((float)CNUM, (float)nchr) - 1.0;


       bend = ((RECL % nchr) ? 1 : 0) + (RECL / nchr);

       cbuf = (char *)calloc(bend*nchr+1, sizeof(char));
       mem_er((cbuf == NULL) ? 0 : 1, (bend*nchr+1) * sizeof(char)); 

       *(cbuf + bend*nchr) = '\0';     

       tmax = tmin = *ata;

       for(i=0; i< dim; i++){

           fv = *(ata + i);
           if(fv > tmax) tmax = fv;
           if(fv < tmin) tmin = fv;

       }

       if(fabs(tmax - tmin) < TOLSCL) {tmin = 0.5 * tmin; tmax = 1.5 *tmax;}

       fprintf(ftav, "%14.7e %17.7e %14.7e\n", tmin, tmax, ssc);

       if(fabs(tmax - tmin) < TOLSCL) {
          printf("****WARNING****, FMAX - FMIN too small.\n\n");
          scl = 0.0;
       }

       else scl = arng/ (tmax - tmin);

       ad = 0;

       for(i=0; i<dim; i++){


          fscl = (*(ata + i) - tmin) * scl;
          ifv = (int) (fscl >= 0.) ? fscl + 0.5 : fscl - 0.5;
          if(ifv > arng) ifv = arng;
          if(ifv < 0) ifv = 0;

          for(j=nchr-1; j>=0; j--){

             ic1 = (ifv % CNUM);
             cbuf[ad*nchr+j] = lkup[ic1];
             ifv /= CNUM;

          }


          if((++ad) == bend){fputs(cbuf, ftav); fprintf(ftav, "\n"); ad = 0;}

      }

      if(ad){
         *(cbuf + ad*nchr) = '\0';
         fputs(cbuf, ftav);
         fprintf(ftav, "\n");
      }

      free(cbuf)

   }

   else if(form == 3){

       fwrite(&ss,4,1,ftav);
       fwrite(pph,sizeof(PP),1,ftav);
       fwrite(&ss,4,1,ftav);
       fwrite(&ss,4,1,ftav);
       fwrite(ata, dim * sizeof(float),1,ftav);
       fwrite(&ss,4,1,ftav);



   }

   else if(form == 4) write_netcdf_field(ata, ftav);

   return;

}


float *field_varience(FILE *fdat, float *mean, int fr1, int fri, int frl, int dim, int nn, int wf, off_t pl1, float *anum, float mval, float *wght, int icmp)

{

   int i, nf=0;
   int nt=0;

   char varf[MAXCHR];

   float diff;
   float wt;

   float *ata=NULL;
   float *anumn=NULL;
   double *acmp=NULL;
   double *acmpn=NULL;


   off_t chrnum=0, place1, place2;

   FILE *ftav=NULL;


   strncpy(varf, TAVGE, MAXCHR);
   if(iext) strcpy(strstr(varf, EXTENSION), fext);
   strcat(varf, "_var");


   printf("***INFORMATON***, Computing Varience of Data........\r\n\n"
          "                                                    \r\n\n"
          " Please Wait, this may take some time ..............  \n\n");

   ata = (float *)calloc(dim, sizeof(float));
   mem_er((ata == NULL) ? 0 : 1, dim * sizeof(float));

   acmp = (double *)calloc(dim, sizeof(double));
   mem_er((acmp == NULL) ? 0 : 1, dim * sizeof(double));

   if(wght){

     acmpn = (double *)calloc(dim, sizeof(double));
     mem_er((acmpn == NULL) ? 0 : 1, dim * sizeof(double));

     anumn = (float *)calloc(dim, sizeof(float));
     mem_er((anumn == NULL) ? 0 : 1, dim * sizeof(float));


   }


/* initialize varience */

   for(i=0; i<dim; i++) *(acmp + i) = 0.0;


/* open file for writing the varience */

   if(form != 4) ftav = open_file(varf, "w");
   else {
      ftav = (FILE *)nc_define((NETCDF_INFO *)fdat, varf);
      ((NETCDF_INFO *)ftav)->nframe = 0;
      ((NETCDF_INFO *)ftav)->iframe = 0;
   }

/* position file pointer */

   if(form != 4) fseeko(fdat, pl1, FSTART);


   while(nf <= frl){


        if(form != 4){
           if(!nf) {

              place1 = ftello(fdat); 
              nlines = 0;

              rf(fdat, ftav, dim, 1);

         
              for(i=0; i<nlines; i++) fprintf(ftav, "%s", ihead[i]);

              if(fr1 == 1) ++nf;;
 
              place2 = ftello(fdat);
              chrnum = place2 - place1;

              if(fr1 > 1){

                 fseeko(ftav, (off_t)0, FSTART);
                 fseeko(fdat, (fr1-2)*chrnum, ORIGIN);
                 nlines = 0;

                 rf(fdat, ftav, dim, 1);

                 for(i=0; i<nlines; i++) fprintf(ftav, "%s", ihead[i]);

                 nf = fr1;           

              }

              nf += fri;

           }

           else {

              fseeko(fdat, (fri-1)*chrnum, ORIGIN);

              if(rf(fdat, ftav, dim, 0)) break;

              nf += fri;

           }

        }

        else {

           if(!nf) nf = fr1;

           ((NETCDF_INFO *)fdat)->iframe = nf - 1;
           
           if(rf(fdat, NULL, dim, 1)) break; 
      
           nf += fri;

        }

        if(wght){

           wt = *(wght + nt);

           for(i=0; i<dim; i++) {

              if(!missing(*(abuf+i), mval, icmp)){

                 diff = *(abuf+i) - *(mean + i);
                 *(acmp + i) += diff * diff;
                 *(acmpn + i) += diff * diff * wt;
                 *(anumn + i) += wt;

              }

           }

        }

        else if(anum){

           for(i=0; i<dim; i++) {

              if(!missing(*(abuf+i), mval, icmp)){

                 diff = *(abuf+i) - *(mean + i);
                 *(acmp + i) += diff * diff;

              }

           }

        }

        else {
           
           for(i=0; i<dim; i++) {

              diff = *(abuf+i) - *(mean + i);
              *(acmp + i) += diff * diff;

           }

        }

        ++nt;

        if(nf > frl) break;

   }

   printf("\n");

/* average the field */

   if(wght){

      for(i=0; i<dim; i++) {
         if(*(anumn + i) > TOLWT)
           *(ata + i) = sqrt(*(acmpn + i) / *(anumn + i));
         else *(ata + i) = DEFVAL;
      }

      writef(ftav, ata, dim);
      for(i=0; i<dim; i++) {
          if(*(anum + i) > TOLWT)
             *(ata + i) = sqrt(*(acmp + i) / *(anum + i));
          else *(ata + i) = DEFVAL;
      }

      free(acmpn);
      free(anumn);

   }

   else if(anum){

      for(i=0; i<dim; i++) {
          if(*(anum + i) > TOLWT)
             *(ata + i) = sqrt(*(acmp + i) / *(anum + i));
          else *(ata + i) = DEFVAL;
      }
      writef(ftav, ata, dim);

   }

   else {
   
      if(nn > 0){
         for(i=0; i<dim; i++) *(ata + i) = sqrt(*(acmp + i) / (float) nn);
      }
      else{
         for(i=0; i<dim; i++) *(ata + i) = DEFVAL;
      }
      writef(ftav, ata, dim);

   }


   if(form != 4) close_file(ftav, varf);
   else netcdf_close((NETCDF_INFO *)ftav);

   free(acmp);

   return ata;

}


/* function to compute the filtered varience of a time series of fields 
   in the frequency domain. Much of the code for this comes from:

   Kay, S. M., 1989: Modern Spectral Estimation: Theory and 
             Application, prentice-hall.

   and the SPECTRA package. 

   This has entailed a re-write of the code in C.                         */

double *window(int , int  , double * );
void periodogram(complex * , int , double * , double , int , double * , double * , int , int );
void dfft(complex * , int , int , int );
int powi(int , int );
void lagcor(int , int , complex * , complex * , complex * , float , int , int );
void blackman_tukey(complex * , complex * , int , double * , int , double * , int , int );
void rednoise(double , double , double * , double * , double , int );
void corcoff1(float * , int , float * , float , int , int );
float confcor(int , int );
void region_grp(REGION * , int , int * );
void write_region_data(FILE * , REGION * , float * , int );


void var_filt_freq(float *ata, float *avar, FILE *fdat, int fr1, int fri, int frl, int dim, int nn, off_t pl1, float *anum, float mval, int icmp)

{

    int i, j, k, l;
    int nband=0;
    int nfreq=0, nfreq2=0;
    int npoint=0;
    int tsamp=0;
    int nblock=0;
    int nf=0;
    int itim=0;
    int pstype=0;
    int m=0;
    int nexp=0;
    int iwin=0;
    int itrig=0;
    int chksum=0;
    int dimtim=0;
    int ntot = 0;
    int iout=0;
    int ospec='n';
    int ixtyp = 0;
    int gnum=0;
    int nseg=0, nlast, nst=0, nover;
    int ilat, ispec=0, irep=0;
    int ired='n';
    int ntype='w', ctype=0;
    int iuseg = 0;
    int orig_var=0;
    int ireg_av=0, ireg_num=0, in_reg=0;
    int av_typ='b';
    int isc_var=1;
    int nt=0;
    int n_filt='n';
    int gbad=0;

    int *index=NULL;
    int *ireg=NULL, *ig=NULL;

    off_t place1, place2, chrnum=0;

    char varf[MAXCHR];
    char gridspec[MAXCHR];
    char phase[MAXCHR];
    char spec_av[MAXCHR];

    char scans;

    float **bltim=NULL, *btim=NULL;
    float *speriod=NULL;
    float *freq=NULL;
    float *bll=NULL, *att=NULL, *avr=NULL;
    float *bands=NULL;
    float *about=NULL;
    float rph=0.0;
    float ptlng, ptlat;
    float corcof=0.0;
    float siglev;
    float wtt=0.0;

    double **aband=NULL, **atemp=NULL;

    double *w=NULL;
    double *pspec=NULL, *tspec=NULL;
    double *phs=NULL;
    double varsum;
    double ratio;
    double wsum=0.0;
    double diff=0.0;


    double scale=1.0, scale2=1.0;
    double fval;
    double var=1.0;
    double *redn=NULL;
    double *redc=NULL;

    WS *wa=NULL;
    REGION *reg_av=NULL, *rtmp=NULL;

    FILE *filvar=NULL;
    FILE *gspec=NULL;
    FILE *fspec_av=NULL;

    complex *timser=NULL, *bt_timser=NULL;
    complex *rcor=NULL;

    strncpy(varf, TAVGE, MAXCHR);
    strncpy(gridspec, TAVGE, MAXCHR);
    strncpy(phase, TAVGE, MAXCHR);
    strncpy(spec_av, TAVGE, MAXCHR);

    if(iext) {
      strcpy(strstr(varf, EXTENSION), fext);
      strcpy(strstr(gridspec, EXTENSION), fext);
      strcpy(strstr(phase, EXTENSION), fext);
      strcpy(strstr(spec_av, EXTENSION), fext);
    }

    strcat(varf, "_varfil");
    strcat(gridspec, "_gspec"); 
    strcat(phase, "_phase");
    strcat(spec_av, "_spec_av");

     
    printf("****WARNING****, computing the band pass filtered varience\r\n"
           "                 may require large amounts of memory.     \n\n");

    printf("****INFORMATION****, the dimension of the data is = %d\n\n", dim);


    printf("How many PSD evaluations per block, large number requires large memory,\r\n"
           "small number requires more I/O.                                        \n\n");

    scanf("%d", &npoint);

    if(npoint > dim) npoint = dim;

    nblock = (dim % npoint) ? (dim / npoint) + 1 : dim / npoint; 

/* assign memory for time series */

    bltim = (float **)calloc(npoint, sizeof(float *));
    mem_er((bltim == NULL) ? 0 : 1, npoint * sizeof(float *));

    for(i=0; i < npoint; i++){

       *(bltim + i) = (float *)calloc(nn, sizeof(float));
       mem_er((*(bltim + i) == NULL) ? 0 : 1, nn * sizeof(float));

    }

    printf("How many bands do you require for filtering?\n\n");

    while(nband <= 0){
       scanf("%d", &nband);
       if(nband <= 0)
         printf("****ERROR****, number of bands must be positive, try again.\n\n");

    }

    bands = (float *)calloc(nband+1, sizeof(float));
    mem_er((bands == NULL) ? 0 : 1, (nband+1) * sizeof(float));

    printf("Do you want any region averaged spectra? Input '1' for yes and '0' for no.\n\n");
    scanf("%d", &ireg_av);

    if(ireg_av){

       ireg = (int *)calloc(dim, sizeof(int));
       mem_er((ireg == NULL) ? 0 : 1, dim * sizeof(int));

       printf("How many regions do you want?\n\n");
       scanf("%d", &ireg_num);

       reg_av = (REGION *)calloc(ireg_num, sizeof(REGION));
       mem_er((reg_av == NULL) ? 0 : 1, ireg_num * sizeof(REGION));

       if(tom == 'g'){

          printf("Do you want Box averaging or Circular (geodesic) averaging, 'b' for box and 'c' for circular.\n\n");
          scanf("\n");
          av_typ = getchar();

       }

       else av_typ = 'b';

       if(!(av_typ == 'b' || av_typ == 'c')){

          printf("****WARNING****, the box averaging tag type %c is not valid, \r\n"
                 "                 defaulting to Box averaging.                \n\n", av_typ);
          av_typ = 'b';

       }


       for(i=0; i < ireg_num; i++){

           rtmp = (reg_av + i);
           rtmp->av_typ = av_typ;
           rtmp->reg_id = i+1;
           rtmp->nband = nband;
           rtmp->band_av_spec = (double *)calloc(nband, sizeof(double));
           mem_er((rtmp->band_av_spec == NULL) ? 0 : 1, nband * sizeof(double));
           rtmp->band_spec_av = (double *)calloc(nband, sizeof(double));
           mem_er((rtmp->band_spec_av == NULL) ? 0 : 1, nband * sizeof(double));
           rtmp->reg_tim = (float *)calloc(nn, sizeof(float));
           mem_er((rtmp->reg_tim == NULL) ? 0 : 1, nn * sizeof(float));
           if(anum){
              rtmp->nav = (int *)calloc(nn, sizeof(int));
              mem_er((rtmp->nav == NULL) ? 0 : 1, nn * sizeof(int));

           }

       }



       printf("****INFORMATION*****, spectra averaging type is %s\n\n", (av_typ == 'b') ? "BOX" : "CIRCULAR");

       printf("What regions do you require; for box, input two (x,y) pairs \r\n"
              "for the lower left and upper right coordinates; for circular\r\n"
              "input the cicle postion (x,y) and angular radius.           \n\n");

       if(av_typ == 'b'){

          for(i=0; i < ireg_num; i++){

              rtmp = (reg_av + i);
              printf("Input region %d coordinates, x1, y1, x2, y2\n", i+1);
              scanf("%f %f %f %f", &(rtmp->x1), &(rtmp->y1), &(rtmp->x2), &(rtmp->y2));
              if((rtmp->x1 > rtmp->x2) || (rtmp->y1 > rtmp->y2)){

                 printf("****ERROR****, box co-ordinates not entered correctely.\r\n"
                        "               Must have X1 < X2 and Y1 < Y2.          \n\n");
                 exit(1);

              }

          }

       }
       else{

          for(i=0; i < ireg_num; i++){

              rtmp = (reg_av + i);
              printf("Input region %d coordinates, x1, y1, rad\n", i+1);
              scanf("%f %f %f", &(rtmp->x1), &(rtmp->y1), &(rtmp->rad));

          }

       }

       printf("\n\n");

/* determine grid points for region */

       for(i=0; i < ireg_num; i++) {
          (reg_av + i)->igrp=NULL; 
          (reg_av + i)->n_grp=0; 
          region_grp(reg_av + i, av_typ, ireg);
          printf("Number of points found in Region %d = %d\n", i+1, (reg_av + i)->n_grp);
       }

       printf("\nFields are normalized to unit varience before calculating spectra,\r\n"
              "do you want to scale individual spectra back to actual varience,  \r\n"
              "before output, input '1' for yes and '0' for no.                  \n\n");

       scanf("%d", &isc_var);

    }


    printf("What is the number of data samples per time unit (e.g. day) for period?\n\n");
    scanf("%d", &tsamp);

    if(tsamp < 1){

       printf("****ERROR****, samples per unit time must be positive.\n\n");
       exit(1);

    }

    printf("Input the period boundaries for filtering, e.g., 0 10 100,\r\n"
           " for 2 bands between 0-10 time units and 10-100 time units\r\n"
           " ordered as high pass to low pass in frequency.           \r\n"
           "Ouput is in the order low pass to high pass.              \n\n");

    i = 0;

    while(i < (nband + 1)) {
         printf("Boundary %d =", i+1);
         scanf("%f", (bands+i));
         printf("\n");
         ++i;
    }

    printf("\n");

/* assign memory for band filtered fields. */

    
    aband = (double **)calloc(nband, sizeof(double *));
    mem_er((aband == NULL) ? 0 : 1, nband * sizeof(double *));

    atemp = (double **)calloc(nband, sizeof(double *));
    mem_er((atemp == NULL) ? 0 : 1, nband * sizeof(double *));

    about = (float *)calloc(dim, sizeof(float));
    mem_er((about == NULL) ? 0 : 1, dim * sizeof(float));

    for(i=0; i < nband; i++){

       *(aband + i) = (double *)calloc(dim, sizeof(double));
       mem_er((*(aband + i) == NULL) ? 0 : 1, dim * sizeof(double));


    }

    printf("What type of PSD estimate do you require?                     \r\n"
           "Input '1' for standard Periodogram,                           \r\n"
           "Input '2' for Welch averaged Periodogram (including weighted),\r\n"
           "Input '3' for Blackman-Tukey,                                 \r\n"
           "Maybe more later, e.g. Maximum Entropy.                       \n\n");

    scanf("%d", &pstype);

    if(pstype > 3 || pstype < 1){

       printf("****WARNING*****, option %d is not available, using option %d\n", pstype, 1);
       pstype = 1;

    }

    printf("How many frequencies are required for the PSD, usualy some power of 2?\n\n");

    if(pstype == 1)
       printf("****INFORMATION****, for the periodogram method number of frequencies \r\n"
              "                     must be >= number of data points %d, it is reset \r\n" 
              "                     otherwize.\n\n", nn);
    else if(pstype == 2)
       printf("****INFORMATION****, for the averaged periodogram method number of frequencies \r\n"
              "                     must be >= window width, it is reset otherwize\n\n");


    while(nfreq <= 0){
         scanf("%d", &nfreq);
         if(nfreq <= 0)
           printf("****ERROR****, number of frequencies must be positive, try again.\n\n");

    }

    if(pstype == 1 || pstype == 2 || pstype == 3){

       if(pstype == 1 && nfreq < nn) nfreq = nn;

       else if(pstype == 2){ 

          printf("Do you want to specify the sections for the averaged\r\n"
                 "periodogram indvidually, together with weights,     \r\n"
                 "'y' or 'n'.                                         \n\n");

          scanf("\n");
          if(getchar() == 'y'){

             iuseg = 1;
             wtt = 0.0;

             printf("How many sections do you want?\n\n");
             scanf("%d", &nseg);

             wa = (WS *)calloc(nseg, sizeof(WS));
             mem_er((wa == NULL) ? 0 : 1, nseg * sizeof(WS));

             printf("Input section data including weights for averaging\r\n"
                    "(set weights to 1 for no weighting.)              \r\n"
                    "Input as: START END WEIGHT                        \n\n");

             for(i=0; i < nseg; i++){

                 printf("Data for section %d :", i+1);
                 scanf("%d %d %f", &((wa+i)->istart), &((wa+i)->iend), &((wa+i)->wt));
                 printf("\n");

                 if((wa+i)->istart < 1 || (wa+i)->iend < 1){

                    printf("****ERROR****, chosen section start or end are out of range, resetting to range.\n\n");
                    if((wa+i)->istart < 1) (wa+i)->istart = 1;
                    if((wa+i)->iend < 1) (wa+i)->iend = nn;

                 }

                 if((wa+i)->istart > nn || (wa+i)->iend > nn){

                    printf("****ERROR****, chosen section start or end are out of range, resetting to range.\n\n");
                    if((wa+i)->istart > nn) (wa+i)->istart = 1;
                    if((wa+i)->iend > nn) (wa+i)->iend = nn;


                 }

                 printf("START = %d, END = %d, WEIGHT = %f\n\n", (wa+i)->istart, (wa+i)->iend, (wa+i)->wt);

                 (wa+i)->m = (wa+i)->iend - (wa+i)->istart + 1;

                 if((wa+i)->m > m) m = (wa+i)->m;

                 wtt += (wa+i)->wt;


             }

             printf("****INFORMATION****, sum of weights for weighted \r\n"
                    "                     averaged periodogram = %f\n\n", wtt);
             
             if(nfreq < m) nfreq = m;

          }

          else {


             printf("What is the window width required for the averaged periodogram.\n\n");
             scanf("%d", &m);
             if(m > nn || m < 4) {

                printf("****WARNING****, window width %d not possible, using %d\n\n", m, nn);

                m = nn;

             }

             printf("How much overlap is required for the sections < window width.\n\n");
             scanf("%d", &nover);

             if(nover >= m - 1 || nover < 0) {
               printf("****WARNING****, overlap %d not allowed, using %d\n\n", nover, 0);
               nover = 0;
             }

             if(nfreq < m) nfreq = m;

/* compute number of overlapping sections */

             nseg = 0;
             nlast = m;
             nst = m - nover;

             while(nlast < nn){

                ++nseg;
                nlast = nst + (m - 1);
                nst = nlast - nover;

             }

             if(nlast == nn) ++nseg;

             nst = m - nover - 1;


          }

          printf("****INFORMATION****, the number of segments used for averaged periodogram is %d\n\n", nseg);

       }
       else if(pstype == 3) {


          printf("****INFORMATION****, for the Blackman-Tukey method number of lags \r\n"
                 "                     should be less than number of data points,   \n\n"
                 "                     ~n/5\n\n");


          printf("What is the maximum lag required for Blackman-Tukey\n\n");
          scanf("%d", &m);

          if(m > nn){

            printf("***WARNING***, too many lags for Blackman-Tukey,\r\n"
                   "               resetting to n/5.                \n\n");
            m = nn / 5;


          }

          if(nfreq < 2*m+2) nfreq = 2 * m + 2;

/* assign memory for lag correlation */

          rcor = (complex *)calloc(m+1, sizeof(complex));
          mem_er((rcor == NULL) ? 0 : 1, (m+1) * sizeof(complex));

          bt_timser = (complex *)calloc(nn, sizeof(complex));
          mem_er((bt_timser == NULL) ? 0 : 1, nn * sizeof(complex));


       }

       nexp = (int)(log((float)nfreq)/log(2.0)); 

       if(powi(2, nexp) != nfreq){
          ++nexp;
          nfreq = powi(2, nexp);

          printf("****WARNING****, number of frequencies must be a power of 2,\r\n" 
                 "                 resetting to %d\n\n", nfreq);

       }

       timser = (complex *)calloc(nfreq, sizeof(complex));
       mem_er((timser == NULL) ? 0 : 1, nfreq * sizeof(complex));


       dimtim = nfreq;


    }

    nfreq2 = nfreq / 2;

    printf("****INFORMATION****, caution should be excersised depending on      \r\n"
           "                     what is required from the time series analysis \r\n"
           "                     if using a data window.                        \r\n"
           "                     If band pass filtered varience is required,    \r\n"
           "                     as well as smoothed spectra then it is best to \r\n"
           "                     use the Bartlett window as this ensures        \r\n"
           "                     positivity. The other windows do not           \r\n"
           "                     ensure a strictley positive spectrum.          \r\n"
           "                     These other windows should only be used to     \r\n"
           "                     explore the nature of the periodicities of the \r\n"
           "                     power spectra.                                 \n\n");


    printf("Do you want tapered windowing of data to reduce sidelobe contanimation, 'y' or 'n'.\n\n");

    scanf("\n");
    if(getchar() == 'y'){

      printf("What window function do you want to use:   \r\n"
             "Input '0' for Bartlett.                    \r\n"
             "Input '1' for Hanning.                     \r\n"
             "Input '2' for Hamming.                     \n\n");

      scanf("%d", &iwin);

      wsum = 0.0;

      if(pstype == 1) w = window(nn/2, iwin, &wsum);

      else if(pstype == 2) {

         if(iuseg){

            for(i=0; i < nseg; i++) 
               (wa+i)->w = window((wa+i)->m / 2, iwin, &((wa+i)->wsum));

         }
         else 
            w = window(m/2, iwin, &wsum);

      }

      else if(pstype == 3) w = window(m, iwin, &wsum);


    }


    printf("****INFORMATION****, The spectral decomposition should sum       \r\n"
           "                     to the total varience, this is certainly    \r\n"
           "                     the case for the periodogram with no data   \r\n"
           "                     windowing. For the other estimators this    \r\n"
           "                     may not be the cse.                         \r\n"
           "                     However, if a weighted averaged periodogram \r\n"
           "                     is required, it is probably not sensible to \r\n"
           "                     do this unless unit weights are used.       \r\n\n"
           " Do you want to normalize to the total varience?, 'y' or 'n'. \n\n");

    scanf("\n");
    orig_var = (getchar() == 'y') ? 1 : 0;



    speriod = (float *)calloc(nfreq2, sizeof(float));
    mem_er((speriod == NULL) ? 0 : 1, nfreq2 * sizeof(float));

    index = (int *)calloc(nfreq2, sizeof(int));
    mem_er((index == NULL) ? 0 : 1, nfreq2 * sizeof(int));

    freq = (float *)calloc(nfreq2, sizeof(float));
    mem_er((freq == NULL) ? 0 : 1, nfreq2 * sizeof(float));

    pspec = (double *)calloc(nfreq, sizeof(double));
    mem_er((pspec == NULL) ? 0 : 1, nfreq * sizeof(double));

    if(pstype == 2){

       tspec = (double *)calloc(nfreq, sizeof(double));
       mem_er((tspec == NULL) ? 0 : 1, nfreq * sizeof(double));

    }

    phs = (double *)calloc(nfreq, sizeof(double));
    mem_er((phs == NULL) ? 0 : 1, nfreq * sizeof(double));

    printf("Do you want to compare and filter with a rednoise         \r\n" 
           "or whitenoise spectrum, 'y' or 'n'.                       \r\n"
           "****WARNING****, this is not yet reliable for all options.\n\n");

    scanf("\n");
    if((ired=getchar()) == 'y'){

       printf("Do you want band filtering with respect to the background spectra significance level, 'y' or 'n'\n\n");
       scanf("\n");
       n_filt = getchar();

       redn = (double *)calloc(nfreq2, sizeof(double));
       mem_er((redn == NULL) ? 0 : 1, nfreq2 * sizeof(double));

       redc = (double *)calloc(nfreq2, sizeof(double));
       mem_er((redc == NULL) ? 0 : 1, nfreq2 * sizeof(double));

       printf("Do you want rednoise or whitenoise, 'r' or 'w'\n\n");
       scanf("\n");
       ntype=getchar();

       if(ntype != 'r' && ntype != 'w'){

         printf("****ERROR****, noise type unknown, defaults to whitenoise\n\n");
         ntype = 'w';

       }

       if(ntype == 'r'){

          printf("Do you want a chosen lag1 correlation coefficient,\r\n"
                 "one computed from the time series, '0' for chosen, \r\n" 
                 "or '1' for computed.\n\n");

          scanf("%d", &ctype);

          if(ctype != 0 && ctype != 1){

             printf("****ERROR****, carrelation type unknown, defaulting to computed.\n\n");
             ctype = 1;


          }

          else if(!ctype){

              printf("What lag1 correaltion coefficient do you want?\n\n");
              scanf("%f", &corcof);

          }

       }

       else corcof = 0.0;

       printf("What significance level do you require, i.e 0.05 is \r\n"
              "the 95%% convidence level, significant at 5%%.        \n\n");

       scanf("%f", &siglev);

       if(siglev < 0 || siglev > 0.1){

          printf("****WARNING****, the significance level is outside of the prescibed\r\n"
                 "                 range of (0, 0.1), using default value of 0.05    \n\n");
          siglev = 0.05;


       }

    }

    if(ireg_av){

       for(i=0; i<ireg_num; i++){

           rtmp = (reg_av + i);
           rtmp->nfreq2 = nfreq2;
           rtmp->av_spec = (double *)calloc(nfreq2, sizeof(double));
           mem_er((rtmp->av_spec == NULL) ? 0 : 1, nfreq2 * sizeof(double));
           rtmp->spec_av = (double *)calloc(nfreq2, sizeof(double));
           mem_er((rtmp->spec_av == NULL) ? 0 : 1, nfreq2 * sizeof(double));

           if(ired == 'y'){

             rtmp->mos_redn = (double *)calloc(nfreq2, sizeof(double));
             mem_er((rtmp->mos_redn == NULL) ? 0 : 1, nfreq2 * sizeof(double));

             rtmp->mos_redc = (double *)calloc(nfreq2, sizeof(double));
             mem_er((rtmp->mos_redc == NULL) ? 0 : 1, nfreq2 * sizeof(double));

             rtmp->som_redn = redn;
             rtmp->som_redc = redc;

           }


       }

    }

/* calculate period */

    for(i=0; i < nfreq2; i++)

        *(freq + i) = i / (float)nfreq;


    *speriod = LARGE_FREQ;

    for(i=1; i < nfreq2; i++)
        *(speriod + i) = 1.0 / ((float)tsamp * *(freq + i));

/* set band indexs for periods */

    for(i=0; i < nfreq2; i++){

        *(index+i) = -1;

        for(j=1; j <= nband; j++){

           if(*(speriod + i) >= *(bands + j - 1) && *(speriod + i) < *(bands + j)) *(index + i) = j - 1;


        }

    }

    printf("Do you want the field scaled, useful for such fields as vorticity, 'y' or 'n'\n\n");
    scanf("\n");
    scans = getchar();

    if(scans == 'y'){

       printf("Input a scaling factor, must be positive.\n\n");
       scanf("%lf", &scale);

       if(scale < 0.){

          printf("****WARNING****, scaling is not positive, no scaling avalibale\n\n");

       }

       else {

          scale2 = scale * scale;

          for(i=0; i < dim; i++) {

              *(ata + i) *= scale;
              *(avar + i) *= scale;

          }

       }


    }

    if(ireg_av){

       printf("do you want to ouput the spectra in terms of frequency or period, \r\n"
              "Input '0' for frequency,                                          \r\n"
              "      '1' for period.                                             \r\n"
              "      '2' for log(period)                                         \n\n");

       scanf("%d", &ixtyp);
       if(ixtyp < 0 || ixtyp >2){
         printf("*****WARNING****, incorrect type specified for spectra output.\r\n"
                 "                  defaulting to frequency.                    \n\n");
         ixtyp = 0;

       }


       printf("Do you want to output the spectra for the grid points of the chosen regions, 'y' or 'n'\n");
       scanf("\n");
       if((ospec = getchar()) == 'y') gspec = open_file(gridspec, "w");
    

    }
      
    printf("***INFORMATON***, Computing Filtered Varience of Data........\r\n\n"
           "                                                             \r\n\n"
           " Please Wait, this may take some time ..............           \n\n");



/* put data into time series for chosen blocks and compute PSD */

   chksum = 0;
   itrig = 1;

   gnum = 0;


   for(i=0; i < nblock; i++){

/* position file pointer */

      chksum = ((i + 1) * npoint < dim) ? npoint : dim - i * npoint;

      if(form != 4) fseeko(fdat, pl1, FSTART);

      bll = abuf + i * npoint;
      att = ata + i * npoint;
      avr = avar + i * npoint;
      if(ireg_av) ig = ireg + i *npoint;

      for(k=0; k<nband; k++) atemp[k] = *(aband + k) + i * npoint;

      itim = 0;
      nf = 0;

/* construct time series */


      while(nf <= frl){

           if(form != 4){
              if(!nf) {

                 place1 = ftello(fdat); 

                 rf(fdat, NULL, dim, 0);         

                 if(fr1 == 1) ++nf;;
 
                 place2 = ftello(fdat);
                 chrnum = place2 - place1;

                 if(fr1 > 1){

                    fseeko(fdat, (fr1-2)*chrnum, ORIGIN);
                    rf(fdat, NULL, dim, 0);

                    nf = fr1;           

                 }

                 nf += fri;


              }

              else {

                 fseeko(fdat, (fri-1)*chrnum, ORIGIN);

                 if(rf(fdat, NULL, dim, 0)) break;
           
                 nf += fri;

              }

           }

           else {

              if(!nf) nf = fr1;

              ((NETCDF_INFO *)fdat)->iframe = nf - 1;
           
              if(rf(fdat, NULL, dim, 1)) break; 
      
              nf += fri;

           }

           if(anum){


              for(j=0; j < chksum; j++){

                  btim = *(bltim + j) + itim;

/* center and normalize time series */

                  fval = *(bll + j);


                  if(missing(fval, mval, icmp)) {

                     *btim = (pstype == 3) ? mval : 0.0;

                  }


                  else {

                     if(scans == 'y') fval *= scale;

                     *btim = (*(avr + j) > TOLVAR) ? (fval - *(att + j)) / *(avr + j) : 0.0;

                  }
 
              }


           }

           else {

              for(j=0; j < chksum; j++){



                  btim = *(bltim + j) + itim;

/* center and normalize time series */

                  fval = *(bll + j);


                  if(scans == 'y') fval *= scale;

                  *btim = (*(avr + j) > TOLVAR) ? (fval - *(att + j)) / *(avr + j) : 0.0;

              }

           }

           if(ireg_av && i == 0){
              for(j=0; j < ireg_num; j++){
                 rtmp = reg_av + j;
                 for(k=0; k < rtmp->n_grp; k++){
                     fval = *(abuf + *(rtmp->igrp + k));
                     if(anum){
                        if(!missing(fval, mval, icmp)){
                          if(scans == 'y') fval *= scale;
                          *(rtmp->reg_tim + itim) += fval;
                          *(rtmp->nav + itim) += 1;

                        }
  
                     }
                     else {
                        if(scans == 'y') fval *= scale;
                        *(rtmp->reg_tim + itim) += fval;
                     }

                 }

              }

           }

           ++itim;

           if(nf > frl) break;

      }

      printf("\n\n");


/* calculate PSD for time series in current block */

      for(j=0; j < chksum; j++){

         if(COUNT) printf("Current Point = %6d\r", i * npoint + j);

         btim = *(bltim + j);

/* put each time series into the complex 1-D array */

         if(pstype == 1) {

             for(k=0; k<dimtim; k++)(timser+k)->real = (timser+k)->imag = 0.0;
           
             for(k=0; k < nn; k++) comp((*(btim + k)), 0.0, (timser+k));

             periodogram(timser, nn, w, wsum, nexp, pspec, phs, nfreq, itrig);

         }

         else if(pstype == 2) {


             for(k=0; k< nfreq; k++) *(pspec + k) = 0.0;

             if(iuseg){

                for(k=0; k < nseg; k++){

                   for(l=0; l < (wa+k)->m; l++) comp((*(btim + l + (wa+k)->istart - 1)), 0.0, (timser+l));

                   periodogram(timser, (wa+k)->m, (wa+k)->w, (wa+k)->wsum, nexp, tspec, phs, nfreq, itrig); 

                   for(l=0; l < nfreq; l++) *(pspec + l) += *(tspec + l) * (wa+k)->wt;


                }

                if(wtt > TOLWT){

                   for(k=0; k< nfreq; k++) *(pspec + k) /= wtt;
                }

                else {

                   for(k=0; k< nfreq; k++) *(pspec + k) = 0.0;;
                }


             }

             else { 
  
                for(k=0; k < nseg; k++){

                   for(l=0; l < m; l++) comp((*(btim + l)), 0.0, (timser+l));
                   btim += nst; 

                   periodogram(timser, m, w, wsum, nexp, tspec, phs, nfreq, itrig); 

                   for(l=0; l < nfreq; l++) *(pspec + l) += *(tspec + l);              

                }

                for(k=0; k< nfreq; k++) *(pspec + k) /= (float) nseg;


             }

         }

         else if(pstype == 3) {

           
             for(k=0; k < nn; k++) comp((*(btim + k)), 0.0, (bt_timser+k));

             lagcor(nn, m+1, bt_timser, bt_timser, rcor, mval, (anum) ? 1 : 0, icmp);

             blackman_tukey(timser, rcor, m, w, nexp, pspec, nfreq, itrig);
 

         }

         itrig = 0;

         varsum = 0.0;


/* normalize back to unit variance */

         if(orig_var){

            for(k=nfreq2; k< nfreq; k++){ 
               varsum += *(pspec+k);
            }

            ratio = (varsum > TOLVAR) ? 1.0 / varsum : 0.0;

            for(k=nfreq2; k< nfreq; k++) *(pspec+k) *= ratio;

         }

/* compute noise spectra if required */

         if(ired == 'y'){


            if(ntype == 'r' && ctype == 1) {

               corcoff1(btim, nn, &corcof, mval, (anum) ? 1 : 0, icmp);


            }

            rednoise(1.0, corcof, redn, redc, (double)siglev, nfreq2);

         }


/* sum spectra and output individual spectra to file */


/* determine if current point is within a required region */

         if(ireg_av){

            if(*(ig + j)){

               for(k=0; k < ireg_num; k++) {

                   rtmp = (reg_av + k);

                   in_reg = 0;

                   for(l=0; l < rtmp->n_grp; l++){

                       if(gnum == *(rtmp->igrp + l) ) {in_reg = k+1; break;}

                   }

                   if(in_reg){

                      rtmp = reg_av + in_reg - 1;

                      ++(rtmp->ntmp);

                      var = 1.0;

                      if(isc_var) var = *(avr + j) * *(avr + j);

                      if(ired == 'y') {
                         for(l=0; l < nfreq2; l++){
                             *(rtmp->av_spec + l) += var * *(pspec + l + nfreq2);
                             *(rtmp->mos_redn + l) += var * *(redn + l);
                             *(rtmp->mos_redc + l) += var * *(redc + l);
                             if(*(index + l) >= 0 && *(pspec+l+nfreq2) >= *(redc + l))
                               *(rtmp->band_av_spec + nband - *(index + l) - 1) += var * *(pspec+l+nfreq2);
                         }

                      }
                      else{

                         for(l=0; l < nfreq2; l++) {
                            *(rtmp->av_spec + l) += var * *(pspec + l + nfreq2);
                            if(*(index + l) >= 0)
                               *(rtmp->band_av_spec + nband - *(index + l) - 1) += var * *(pspec+l+nfreq2);

                         }

                      }


                      if(ospec == 'y'){

                         ilat = gnum / gr->ix;
                         ptlat = *(gr->ygrid + ilat);
                         ptlng = *(gr->xgrid + gnum - ilat * gr->ix);


                         fprintf(gspec, "REGION (Id., Type, Extent) %d %c %f %f %f %f\n", in_reg, av_typ, rtmp->x1, rtmp->y1, rtmp->x2, rtmp->y2);
                         fprintf(gspec, "Position (Long., Lat.) %f %f\n", ptlng, ptlat);
                         fprintf(gspec, "%f %f\n", corcof, confcor(nn, 2));

                         if(!ixtyp){ 

                            rph = 0.0;
                            for(l=0; l < nfreq2; l++) {

                               rph =  *(phs + l + nfreq2);
                               if(rph < 0.0) rph += FPI2;
                               rph /= FPI2;

                               if(ired == 'y')
                                 fprintf(gspec, "%f %f %f %f\n", *(freq + l), var * *(pspec + l + nfreq2), var * *(redn + l), var * *(redc + l));
                               else
                                 fprintf(gspec, "%f %f\n", *(freq + l), var * *(pspec + l + nfreq2));

                            }

                         }

                         else if(ixtyp == 1){

                            for(l=0; l < nfreq2; l++) {

                               rph = *(phs + l + nfreq2);
                               if(rph < 0.0) rph += FPI2;
                               rph = rph / FPI2;

                               if(ired == 'y')
                                  fprintf(gspec, "%f %f %f %f\n", *(speriod + l), var * *(pspec + l + nfreq2), var * *(redn + l), var * *(redc + l));
                               else
                                  fprintf(gspec, "%f %f\n", *(speriod + l), var * *(pspec + l + nfreq2));

                            }

                         }

                         else if(ixtyp == 2){

                            for(l=0; l < nfreq2; l++) {

                               rph = *(phs + l + nfreq2);
                               if(rph < 0.0) rph += FPI2;
                               rph = rph / FPI2;

                               if(ired == 'y')
                                  fprintf(gspec, "%f %f %f %f\n", log(*(speriod + l)), var * *(pspec + l + nfreq2), var * *(redn + l), var * *(redc + l));
                               else
                                fprintf(gspec, "%f %f\n", log(*(speriod + l)), var * *(pspec + l + nfreq2));

                            }

                         }

                      }

                   }

               }

            }

         }


         var = *(avr + j) * *(avr + j);

         ispec = 0;

         for(k=nfreq2; k< nfreq; k++) {
            if(*(pspec+k) < 0.0) ispec = 1;
            *(pspec+k) *= var;
         }



         if(ired == 'y'){

            for(k=0; k< nfreq2; k++){ *(redc + k) *= var; *(redn + k) *= var;}

         }


         if(ispec == 1){

            ++gbad;

            if(!irep){
               printf("****WARNING****, a possible non-positive spectrum has been detected.\r\n"
                   "                 power spectra are not reliable.                    \n\n");

               irep = 1;
            }

         }


/* filter varience into bands */

         for(k=0; k < nband; k++) *atemp[k] = 0.0;

         if(ired == 'y' && n_filt == 'y'){

            for(k=0; k< nfreq2; k++){

               if(*(index + k) >= 0 && *(pspec+k+nfreq2) >= *(redc + k))

                  *atemp[nband - *(index + k) - 1] += *(pspec+k+nfreq2);

            }


         }

         else{

            for(k=0; k< nfreq2; k++){

               if(*(index + k) >= 0)

                  *atemp[nband - *(index + k) - 1] += *(pspec+k+nfreq2);

            }

         }

         for(k=0; k < nband; k++) ++atemp[k];

         ++gnum;

      }

      ntot += chksum;

      printf("Completed %d points so far.\n", ntot);

    }

    printf("\n");

    printf("****INFORMATION****, number of grid points with a non-positive \r\n"
           "                     spectrum is %d out of a total of %d\n\n", gbad, dim);


/* write to file the average of the spectra and spectra of the average */


    if(ireg_av){

      if(ospec == 'y') close_file(gspec, gridspec);

      fspec_av = open_file(spec_av, "w");


/* compute spectra of region averaged time series */

/* average MOS data for each region, compute SOM and generate SOM varience if required */

      fprintf(fspec_av, "!\n!\n!");
      fprintf(fspec_av, "%50s\n!\n", "Region Averaged Power Spectra");


      for(i=0; i < ireg_num; i++){

          rtmp = reg_av + i;

          if(rtmp->ntmp != rtmp->n_grp){

             printf("****WARNING****, number of spectra to average over \r\n"
                    "                 not equal to chosen number of grid\r\n"
                    "                 points for region %d, n1 = %d, n2 = %d\n\n", i+1, rtmp->ntmp, rtmp->n_grp);

          }

          for(k=0; k < nband; k++) *(rtmp->band_av_spec + k) /= rtmp->n_grp;

          for(k=0; k<nfreq2; k++) *(rtmp->av_spec + k) /= rtmp->n_grp;

          if(ired == 'y'){

             for(k=0; k<nfreq2; k++){
              *(rtmp->mos_redn + k) /= rtmp->n_grp;
              *(rtmp->mos_redc + k) /= rtmp->n_grp;
             }

          }


/* compute Spectra of Mean */

          if(anum){

             nt = 0;

             for(j=0; j < nn; j++) {
 
                if(*(rtmp->nav + j)){
                   ++nt;
                   *(rtmp->reg_tim + j) /= *(rtmp->nav + j);
                   rtmp->av_mn += *(rtmp->reg_tim + j);

                }


             }


             rtmp->av_mn = (nt) ? rtmp->av_mn / nt : 0.0;
    

             for(j=0; j < nn; j++) {

                if(*(rtmp->nav + j)){

                  diff = *(rtmp->reg_tim + j) - rtmp->av_mn;
                  rtmp->av_var += diff * diff;

                }

             }

             rtmp->av_var = (nt) ? rtmp->av_var / nt : 0.0;
             rtmp->av_var = sqrt(rtmp->av_var);


             if(rtmp->av_var < TOLVAR){

                printf("****WARNING****, zero varience, no spectra SOM for Region %d\n\n", i+1);

             }

             else {

                for(j=0; j < nn; j++){
                   if(*(rtmp->nav + j)){
                     *(rtmp->reg_tim + j) = (*(rtmp->reg_tim + j) - rtmp->av_mn) / rtmp->av_var;
                   }
                   else *(rtmp->reg_tim + j) = 0.0;

                }

             }

          }

          else {
             for(j=0; j < nn; j++) {
                *(rtmp->reg_tim + j) /= rtmp->n_grp;
                rtmp->av_mn += *(rtmp->reg_tim + j);
             }

             rtmp->av_mn /= nn;

             for(j=0; j < nn; j++){
                diff = *(rtmp->reg_tim + j) - rtmp->av_mn;
                rtmp->av_var += diff * diff;
             } 

             rtmp->av_var /= nn;
             rtmp->av_var = sqrt(rtmp->av_var);

             if(rtmp->av_var < TOLVAR){

                printf("****WARNING****, zero varience, no spectra SOM for Region %d\n\n", i+1);

             }

             else {
                for(j=0; j < nn; j++)
                   *(rtmp->reg_tim + j) = (*(rtmp->reg_tim + j) - rtmp->av_mn) / rtmp->av_var;
             }

          }


/* still need to finish with missing values */


          if(pstype == 1) {

             for(k=0; k<dimtim; k++)(timser+k)->real = (timser+k)->imag = 0.0;
           
             for(k=0; k < nn; k++) comp((*(rtmp->reg_tim + k)), 0.0, (timser+k));

             periodogram(timser, nn, w, wsum, nexp, pspec, phs, nfreq, itrig);

          }

          else if(pstype == 2) {


             for(k=0; k< nfreq; k++) *(pspec + k) = 0.0;

             if(iuseg){

                for(k=0; k < nseg; k++){

                   for(l=0; l < (wa+k)->m; l++) comp((*(rtmp->reg_tim + l + (wa+k)->istart - 1)), 0.0, (timser+l));

                   periodogram(timser, (wa+k)->m, (wa+k)->w, (wa+k)->wsum, nexp, tspec, phs, nfreq, itrig); 

                   for(l=0; l < nfreq; l++) *(pspec + l) += *(tspec + l) * (wa+k)->wt;


                }

                if(wtt > TOLWT){

                   for(k=0; k< nfreq; k++) *(pspec + k) /= wtt;
                }

                else {

                   for(k=0; k< nfreq; k++) *(pspec + k) = 0.0;;
                }


             }

             else { 
  
                for(k=0; k < nseg; k++){

                   for(l=0; l < m; l++) comp((*(rtmp->reg_tim + l)), 0.0, (timser+l));
                   rtmp->reg_tim += nst; 

                   periodogram(timser, m, w, wsum, nexp, tspec, phs, nfreq, itrig); 

                   for(l=0; l < nfreq; l++) *(pspec + l) += *(tspec + l);              

                }

                for(k=0; k< nfreq; k++) *(pspec + k) /= (float) nseg;

             }

          }

          else if(pstype == 3) {

           
             for(k=0; k < nn; k++) comp((*(rtmp->reg_tim + k)), 0.0, (bt_timser+k));

             lagcor(nn, m+1, bt_timser, bt_timser, rcor, mval, (anum) ? 1 : 0, icmp);

             blackman_tukey(timser, rcor, m, w, nexp, pspec, nfreq, itrig);
 

          }

          varsum = 0.0;

          ispec = 0;

/* normalize back to unit variance */


          if(orig_var){

            for(k=nfreq2; k< nfreq; k++){ 
               if(*(pspec+k) < 0.0) ispec = 1;
               varsum += *(pspec+k);
            }
            ratio = (varsum > TOLVAR) ? 1.0 / varsum : 0.0;

            for(k=nfreq2; k< nfreq; k++) *(pspec+k) *= ratio;

          }

          if(ispec == 1){

            printf("****WARNING****, a possible non-positive spectrum has been detected.\r\n"
                   "                 for the Spectra of the Mean, power spectra are not \r\n"
                   "                 reliable.                                          \n\n");

          }

          var = 1.0;

          if(isc_var) var = rtmp->av_var * rtmp->av_var;


/* compute noise spectra if required */

          if(ired == 'y'){


            if(ntype == 'r' && ctype == 1) {

               corcoff1(rtmp->reg_tim, nn, &corcof, mval, (anum) ? 1 : 0, icmp);


            }

            rednoise(1.0, corcof, redn, redc, (double)siglev, nfreq2);

            for(k=0; k< nfreq2; k++){ *(redc + k) *= var; *(redn + k) *= var;}

          }



          if(ired == 'y' && n_filt == 'y') {
             for(k=0; k < nfreq2; k++){
                 *(pspec + k + nfreq2) *= var;
                 *(rtmp->spec_av + k) += *(pspec + k + nfreq2);
                 if(*(index + k) >= 0 && *(pspec+k+nfreq2) >= *(redc + k))
                   *(rtmp->band_spec_av + nband - *(index + k) - 1) += *(pspec+k+nfreq2);
             }

          }
          else{

             for(k=0; k < nfreq2; k++) {
                *(pspec + k + nfreq2) *= var;
                *(rtmp->spec_av + k) += *(pspec + k + nfreq2);
                if(*(index + k) >= 0)
                   *(rtmp->band_spec_av + nband - *(index + k) - 1) += *(pspec+k+nfreq2);

             }

          }


/* write region averaged spectra to file */


          if(!ixtyp) write_region_data(fspec_av , rtmp, freq, ixtyp); 
          else if(ixtyp == 1 || ixtyp == 2) write_region_data(fspec_av , rtmp, speriod, ixtyp);


      }

      close_file(fspec_av, spec_av);
    }


    if(pstype == 1 || pstype == 2) dfft(NULL, 0, 1, -1);

/* write to file the various band filtered fields */

    printf("What kind of output do you require?             \r\n"
           "Input '0' for varience (power).                 \r\n"
           "Input '1' for normalized varience.              \r\n"
           "Input '2' for RMS sqrt(varience).               \r\n"
           "Input '3' for Power in DB (10*log_10(varience)).\n\n");

    scanf("%d", &iout);

    if(form != 4) filvar = open_file(varf, "w");
    else {
       filvar = (FILE *)nc_define((NETCDF_INFO *)fdat, varf);
       ((NETCDF_INFO *)filvar)->iframe = 0;
       ((NETCDF_INFO *)filvar)->nframe = 0;
    }

    for(i=0; i < nband; i++){

       switch(iout){
         case 0:
          if(scans == 'y'){
            for(j=0; j < dim; j++) *(about + j) = *(*(aband + i) + j) / scale2;
          }
          else {
            for(j=0; j < dim; j++) *(about + j) = *(*(aband + i) + j);
          }
          break;
         case 1:
          for(j=0; j < dim; j++) *(about + j) = *(*(aband + i) + j) / (*(avar + j) * *(avar + j)) ;
          break;
         case 2:  
          if(scans == 'y'){
             for(j=0; j < dim; j++) *(about + j) = sqrt(*(*(aband + i) + j)) / scale;
          }
          else{
             for(j=0; j < dim; j++) *(about + j) = sqrt(*(*(aband + i) + j));
          }
          break;
         case 3:
          if(scans == 'y'){
             for(j=0; j < dim; j++) *(about + j) = 10.0 * log10(*(*(aband + i) + j) / scale2);
          }
          else{
             for(j=0; j < dim; j++) *(about + j) = 10.0 * log10(*(*(aband + i) + j));
          }
          break;
         default:
          printf("****WARNING****, incorrect identifier, defaulting to varience output.\n\n");
          if(scans == 'y'){
            for(j=0; j < dim; j++) *(about + j) = *(*(aband + i) + j) / scale2;
          }
          else {
            for(j=0; j < dim; j++) *(about + j) = *(*(aband + i) + j);
          }

       }

       if(form != 4){
          for(j=0; j< nlines; j++)fprintf(filvar, "%s", ihead[j]);
       }
       else
          ((NETCDF_INFO *)filvar)->nframe = i;

       writef(filvar, about, dim);


    }

    if(form != 4) close_file(filvar, varf);
    else netcdf_close((NETCDF_INFO *)filvar);

    if(ireg_av){
       free(ireg);
       for(i=0; i < ireg_num; i++) {
          rtmp = reg_av + i;
          free(rtmp->igrp);
          free(rtmp->band_av_spec);
          free(rtmp->band_spec_av);
          free(rtmp->av_spec);
          free(rtmp->spec_av);
          free(rtmp->reg_tim);
          free(rtmp->mos_redn);
          free(rtmp->mos_redc);
          free(rtmp->nav);
       }
    }
    free(reg_av);
    free(speriod);
    free(freq);
    for(i=0; i < nband; i++) free(*(aband + i));
    free(aband);
    free(atemp);
    for(i=0; i < npoint; i++) free(*(bltim + i));
    free(bltim);
    free(timser);
    free(bt_timser)
    free(pspec);
    free(w);
    free(bands);
    free(index);
    free(about);
    free(tspec);
    free(rcor);
    free(redn);
    free(redc);
    if(wa){
      for(i=0; i< nseg; i++) free((wa+i)->w);
    }
    free(wa);
    free(phs);


    return;

}




int sub_tim_avg(FILE *fdat, float *mean, int fr1, int fri, int frl, int dim, off_t pl1, float *anum, float mval, int icmp)

{

   int i, nf=0;
   int nn=0;

   char tsub[MAXCHR];

   off_t chrnum=0, place1, place2;

   FILE *fsub=NULL;


   strncpy(tsub, TAVGE, MAXCHR);
   if(iext) strcpy(strstr(tsub, EXTENSION), fext);
   strcat(tsub, "_sub");


   printf("***INFORMATON***, Subtracting time mean from all fields........\r\n\n"
          "                                                               \r\n\n"
          " Please Wait, this may take some time ..............             \n\n");




/* open file for writing the varience */

   if(form != 4) {
      fsub = open_file(tsub, "w");

/* position file pointer */

      fseeko(fdat, pl1, FSTART);

   }
   else fsub = (FILE *)nc_define((NETCDF_INFO *)fdat, tsub);


   while(nf <= frl){

        if(form != 4) {

           if(!nf) {

              place1 = ftello(fdat); 
              nlines = 0;

              rf(fdat, NULL, dim, 1);
       
              for(i=0; i<nlines; i++) fprintf(fsub, "%s", ihead[i]);

              if(fr1 == 1) ++nf;
 
              place2 = ftello(fdat);
              chrnum = place2 - place1;

              if(fr1 > 1){

                 fseeko(fsub, (off_t)0, FSTART);
                 fseeko(fdat, (fr1-2)*chrnum, ORIGIN);
                 nlines = 0;

                 rf(fdat, NULL, dim, 1);

                 for(i=0; i<nlines; i++) fprintf(fsub, "%s", ihead[i]);

                 nf = fr1;           

              }

              nf += fri;

           }
  
           else {

              fseeko(fdat, (fri-1)*chrnum, ORIGIN);

              if(rf(fdat, NULL, dim, 1)) break;

              for(i=0; i<nlines; i++) fprintf(fsub, "%s", ihead[i]);

              nf += fri;

           }

        }

        else {

           if(!nf) nf = fr1;

           ((NETCDF_INFO *)fdat)->iframe = nf - 1;
           
           if(rf(fdat, NULL, dim, 1)) break; 

           ((NETCDF_INFO *)fsub)->iframe = nf - 1;
           ((NETCDF_INFO *)fsub)->nframe = nn;
      
           nf += fri;

        }

        if(anum){

           for(i=0; i < dim; i++){

               if(!missing(*(abuf+i), mval, icmp)) *(abuf+i) -= *(mean + i);

           }


        }

        else {
           
           for(i=0; i<dim; i++) {

              *(abuf+i) -= *(mean + i);

           }

        }

        writef(fsub, abuf, dim);

        ++nn;

        if(nf > frl) break;

   }

   printf("\n");

   if(form != 4) close_file(fsub, tsub);
   else netcdf_close((NETCDF_INFO *)fsub);

   return nn;

}


void region_grp(REGION *reg, int type, int *ireg)
{


    int i,j;
    int *ig=NULL;
    int iper=0;
    int ix=0, iy=0;
    float gv;
    float yy, xx;

    float x1, y1, z1;
    float x2, y2, z2;
    float dprd=0.0;
    double alat, alng, con;
    double s1, c1, s2, c2;

    if(form < 2 || form == 4) {ix = std_x; iy = std_y;}
    else if (form == 2){
      if(utfv == 3) {ix = gr->ix; iy = ((eqsw) ? gr->iy + 1 : gr->iy);}
      else if(utfv == 4) {ix = i_utf4x; iy = i_utf4y;}

    }

    else if(form == 3) {ix = pph->lbnpt; iy = pph->lbrow;}


    if(fabs(*(gr->xgrid + ix - 1) - *(gr->xgrid) - period) < TOLGRID) iper = gr->ix - 1;
    else iper = ix;

    if(type == 'b'){

       for(i=0; i<iy; i++){
           gv = *(gr->ygrid + i);
           yy = (reg->y2 - gv) * (gv - reg->y1);

           if(yy >= 0.0){

              for(j=0; j < iper; j++){
                  gv = *(gr->xgrid + j);
                  xx = (reg->x2 - gv) * (gv - reg->x1);
                  if(xx >= 0.0){
                     if(!reg->igrp){

                        reg->igrp = (int *)malloc_initl(sizeof(int));
                        mem_er((reg->igrp == NULL) ? 0 : 1, sizeof(int));
                        if(form == 2) *(reg->igrp) = (iy - i - 1) * gr->ix + j;
                        else *(reg->igrp) = i * gr->ix + j;
                        *(ireg + *(reg->igrp)) = 1;
                        reg->n_grp = 1;

                     }
                     else {
                        ++(reg->n_grp);
                        reg->igrp = (int *)realloc_n(reg->igrp, reg->n_grp * sizeof(int));
                        mem_er((reg->igrp == NULL) ? 0 : 1, reg->n_grp * sizeof(float));
                        ig = reg->igrp +  reg->n_grp - 1;
                        if(form == 2) *ig = (iy - i - 1) * gr->ix + j;
                        else *ig = i * gr->ix + j;
                        *(ireg + *ig) = 1;

                     }


                  }

              }

           }

       }

    }

    else if(type == 'c'){

       alat = FP_PI2 - reg->y1 * FP_PI;
       alng = reg->x1 * FP_PI;

       con = 1.0 / cos(reg->rad * FP_PI);

       sincos(alng, &s1, &c1);
       sincos(alat, &s2, &c2);

       x1 = s2 * c1;
       y1 = s2 * s1;
       z1 = c2;


       for(i=0; i<iy; i++){

           alat = FP_PI2 - *(gr->ygrid + i) * FP_PI;
           sincos(alat, &s2, &c2);

           for(j=0; j < iper; j++){

               alng = *(gr->xgrid + j) * FP_PI;
               sincos(alng, &s1, &c1);

               x2 = s2 * c1;
               y2 = s2 * s1;
               z2 = c2;

               dprd = con * (x1 * x2 + y1 * y2 + z1 * z2) - 1.0;

               if(dprd >= 0.0){
                  if(!reg->igrp){

                     reg->igrp = (int *)malloc_initl(sizeof(int));
                     mem_er((reg->igrp == NULL) ? 0 : 1, sizeof(int));
                     if(form == 2) *(reg->igrp) = (iy - i - 1) * gr->ix + j;
                     else *(reg->igrp) = i * gr->ix + j;
                     *(ireg + *(reg->igrp)) = 1;
                     reg->n_grp = 1;

                  }
                  else {
                     ++(reg->n_grp);
                     reg->igrp = (int *)realloc_n(reg->igrp, reg->n_grp * sizeof(int));
                     mem_er((reg->igrp == NULL) ? 0 : 1, reg->n_grp * sizeof(float));
                     ig = reg->igrp +  reg->n_grp - 1;
                     if(form == 2) *ig = (iy - i - 1) * gr->ix + j;
                     else *ig = i * gr->ix + j;
                     *(ireg + *ig) = 1;


                  }

               }

           }

       }

    }

    else {
       printf("****ERROR****, averaging specifier not a valid type, exiting\n\n");
       exit(1);

    }

    return;

}


void write_region_data(FILE *fspec , REGION *rtmp, float *ff, int ixtyp)

{

    int i;

    char cfq[MAXCHR];

    float frq=0.0;

    if(!ixtyp) strncpy(cfq, "Freq.", MAXCHR);
    else if (ixtyp == 1) strncpy(cfq, "Period", MAXCHR);
    else if (ixtyp == 2) strncpy(cfq, "Log(Period)", MAXCHR);


    if(rtmp->av_typ == 'b') 
      fprintf(fspec, "!Region %3d; X1 = %f, Y1 = %f, X2 = %f, Y2 = %f\n!\n", rtmp->reg_id, rtmp->x1, rtmp->y1, rtmp->x2, rtmp->y2);
    else if(rtmp->av_typ == 'c')
        fprintf(fspec, "!Region %3d; Xc = %f, Yc = %f, Rad = %f\n!\n", rtmp->reg_id, rtmp->x1, rtmp->y1, rtmp->rad);

    fprintf(fspec, "!=========================================================================\n!\n");

    fprintf(fspec, "!Mean of Spectra (MOS)\r\n"
                   "!=====================\n!\n");

    fprintf(fspec, "!Band Filter :-\r\n"
                   "!--------------\n!\n");

    fprintf(fspec, "!BAND:- ");
    for(i=0; i < rtmp->nband; i++) fprintf(fspec, "%9s%03d ", "Band", i+1);
    fprintf(fspec, "\n");
    fprintf(fspec, "!--------------------------------------------------------------------------\n");
    fprintf(fspec, "           ");
    for(i=0; i < rtmp->nband; i++) fprintf(fspec, "%12e ", *(rtmp->band_av_spec + i));
    fprintf(fspec, "\n!\n");


    if(rtmp->mos_redn){

      fprintf(fspec, "!Spectra:-%12s    %6s      %12s  %12s\n", cfq, "MOS", "Red Noise", "Confidence");
      fprintf(fspec, "!--------------------------------------------------------------------------\n");
      for(i=0; i<rtmp->nfreq2; i++){
          if(ixtyp == 2) frq = log(*(ff + i));
          else frq = *(ff + i);
          fprintf(fspec, "           %12e  %12e  %12e  %12e\n", frq, *(rtmp->av_spec + i), *(rtmp->mos_redn + i), *(rtmp->mos_redc + i));
      }
      fprintf(fspec, "!\n!\n");

    }

    else{

      fprintf(fspec, "!Spectra:-%12s      %6s\n", cfq, "MOS");
      fprintf(fspec, "!--------------------------------------------------------------------------\n");
      for(i=0; i<rtmp->nfreq2; i++){
          if(ixtyp == 2) frq = log(*(ff + i));
          else frq = *(ff + i);
          fprintf(fspec, "             %12e  %12e\n", frq, *(rtmp->av_spec + i));
      }

    }
    fprintf(fspec, "!\n!\n");


    fprintf(fspec, "!=========================================================================\n!\n");

    fprintf(fspec, "!Spectra of Mean (SOM)\r\n"
                   "!=====================\n!\n");

    fprintf(fspec, "!Band Filter :-\r\n"
                   "!--------------\n!\n");

    fprintf(fspec, "!BAND:- ");
    for(i=0; i < rtmp->nband; i++) fprintf(fspec, "%9s%03d ", "Band", i+1);
    fprintf(fspec, "\n");
    fprintf(fspec, "!--------------------------------------------------------------------------\n");
    fprintf(fspec, "           ");
    for(i=0; i < rtmp->nband; i++) fprintf(fspec, "%12e ", *(rtmp->band_spec_av + i));
    fprintf(fspec, "\n!\n");

    if(rtmp->mos_redn){

      fprintf(fspec, "!Spectra:-%12s   %6s      %12s  %12s\n", cfq, "SOM", "Red Noise", "Confidence");
      fprintf(fspec, "!--------------------------------------------------------------------------\n");
      for(i=0; i<rtmp->nfreq2; i++){
          if(ixtyp == 2) frq = log(*(ff + i));
          else frq = *(ff + i);
          fprintf(fspec, "           %12e  %12e  %12e  %12e\n", frq, *(rtmp->spec_av + i), *(rtmp->som_redn + i), *(rtmp->som_redc + i));
      }
      fprintf(fspec, "!\n!\n");

    }

    else{

      fprintf(fspec, "!Spectra:-%12s      %6s\n", cfq, "SOM");
      fprintf(fspec, "!--------------------------------------------------------------------------\n");
      for(i=0; i<rtmp->nfreq2; i++){
          if(ixtyp == 2) frq = log(*(ff + i));
          else frq = *(ff + i);
          fprintf(fspec, "             %12e  %12e\n", frq, *(rtmp->spec_av + i));
      }
      fprintf(fspec, "!\n!\n");

    }

    return;

}

/* function to combine weighted time averaged statistics */ 

void time_tele_avg(FILE *fdat)
{

   int i, j, k;
   int dim=0;
   int ifrt=0;

   char infil[5][MAXCHR], *ifl=NULL;
   char outfl[4][MAXCHR];

   FILE *fstat[5]={NULL, NULL, NULL, NULL, NULL};

   float *avg[4]={NULL, NULL, NULL, NULL};
   double *avgt[4]={NULL, NULL, NULL, NULL};

   double specsum=0.0;

   double a1, a2, a3, a4, a5, diff;

   float *specwt=NULL;

   printf("****WARNING****, data should come from a previous run of the time averaging and filtering\r\n"
          "                 part of the code, unpredictable results may occur otherwise.            \n\n");


/* assign space for header strings */

   nlines = 0;

   for(i=0; i < NUMLINES; i++){
       ihead[i] = (char *)calloc(MAXCHR, sizeof(char));
       mem_er((ihead[i] == NULL) ? 0 : 1, MAXCHR * sizeof(char));
   }

   if(form < 2 || form == 4) dim = std_x * std_y;
   else if (form == 2){
      if(utfv == 3) dim=gr->ix * ((eqsw) ? gr->iy + 1 : gr->iy);
      else if(utfv == 4) dim = i_utf4x * i_utf4y;

   }

   else if(form == 3) dim = pph->lbrow * pph->lbnpt;

   if(form == 2){
      printf("How mant frames are there?\n\n");
      scanf("%d", &frnum);
   }

/* assign memory to arrays */

   for(i=0; i < 4; i++){
      avg[i] = (float *)calloc(dim, sizeof(float));
      mem_er((avg[i] == NULL) ? 0 : 1, dim * sizeof(float));
      avgt[i] = (double *)calloc(dim, sizeof(double));
      mem_er((avgt[i] == NULL) ? 0 : 1, dim * sizeof(double));
   }

   specwt = (float *)calloc(frnum, sizeof(float));
   mem_er((specwt == NULL) ? 0 : 1, frnum * sizeof(float));

/* open files */

   ifl = &infil[0][0];

   printf("What is the file with time averages?\n\n");
   scanf("%s", ifl);

   if(form != 4) fstat[0] = open_file(ifl, "r");
   else {
     fstat[0] = (FILE *)nc_clone((NETCDF_INFO *)fdat, ifl, NC_OPEN_MODE);
     ((NETCDF_INFO *)fstat[0])->iframe = 0;
     ((NETCDF_INFO *)fstat[0])->nframe = 0;
   }

   ifl = &infil[1][0];

   printf("What is the file with the STD's (not variances)?\n\n");
   scanf("%s", ifl);

   if(form != 4) fstat[1] = open_file(ifl, "r");
   else {
      fstat[1] = (FILE *)nc_clone((NETCDF_INFO *)fdat, ifl, NC_OPEN_MODE);
      ((NETCDF_INFO *)(fstat[1]))->iframe = 0;
      ((NETCDF_INFO *)(fstat[1]))->nframe = 0;
   }

   ifl = &infil[2][0];

   printf("What is the file with filtered variances?\n\n");
   scanf("%s", ifl);

   if(form != 4) fstat[2] = open_file(ifl, "r");
   else {
      fstat[2] = (FILE *)nc_clone((NETCDF_INFO *)fdat, ifl, NC_OPEN_MODE);
      ((NETCDF_INFO *)(fstat[2]))->iframe = 0;
      ((NETCDF_INFO *)(fstat[2]))->nframe = 0;
   }

   ifl = &infil[3][0];

   printf("What is the file with gridded sum of weights?\n\n");
   scanf("%s", ifl);

   if(form != 4) fstat[3] = open_file(ifl, "r");
   else {
      fstat[3] = (FILE *)nc_clone((NETCDF_INFO *)fdat, ifl, NC_OPEN_MODE);
      ((NETCDF_INFO *)(fstat[3]))->iframe = 0;
      ((NETCDF_INFO *)(fstat[3]))->nframe = 0;
   }

   ifl = &infil[4][0];

   printf("What is the file with period sum of weights?\n\n");
   scanf("%s", ifl);

   fstat[4] = open_file(ifl, "r");

   for(i=0; i < frnum; i++) fscanf(fstat[4], "%f", specwt + i);

   close_file(fstat[4], ifl);

   for(i=0; i < 4; i++){
       strcpy(&outfl[i][0], &infil[i][0]);
       strcat(&outfl[i][0], "_comb");
   }

   nlines = 0;
   ifrt = 1;

   for(i=0; i < frnum; i++){

       for(j=0; j < 4; j++){
          abuf = avg[j];
          rf(fstat[j], NULL, dim, ifrt);

          ifrt = 0;
       }

       a5 = specsum + *(specwt + i);

       if(i == 0){
          for(j=0; j < 4; j++){
              for(k=0; k < dim; k++) *(avgt[j] + k) = (double)(*(avg[j] + k));
          }

       }

       else {
          for(j=0; j < dim; j++){
              a1 = *(avgt[3] + j) * *(avgt[0] + j) + *(avg[3] + j) * *(avg[0] + j);
              a2 = *(avgt[3] + j) * *(avgt[1] + j) * *(avgt[1] + j) + *(avg[3] + j) * *(avg[1] + j) * *(avg[1] + j);
              diff = *(avgt[0] + j) - *(avg[0] + j);
              a3 = *(avgt[3] + j) * *(avg[3] + j) * diff * diff;
              a4 = *(avgt[3] + j) + *(avg[3] + j);

              *(avgt[0] + j) = a1 / a4;
              *(avgt[1] + j) = sqrt((a2 / a4) + (a3 / (a4 * a4)));
              *(avgt[3] + j) = a4;


              *(avgt[2] + j) = (specsum * *(avgt[2] + j) + *(specwt + i) * *(avg[2] + j)) / a5;

          }
       }

       specsum = a5;

       if(form == 4){
          for(j=0; j < 4; j++) ((NETCDF_INFO *)(fstat[j]))->iframe += 1;
       }

   }

   for(i=0; i < 4; i++){
      if(fstat + i){
         if(form != 4) close_file(fstat[i], &infil[i][0]);
         else netcdf_close((NETCDF_INFO *)(fstat[i]));
      }
   }

/* write averages to file */

   for(i=0; i < 4; i++){

       ifl = &outfl[i][0];
       if(form != 4) fstat[0] = open_file(ifl, "w");
       else {
         fstat[0] = (FILE *)nc_define((NETCDF_INFO *)fdat, ifl);
         ((NETCDF_INFO *)(fstat[0]))->iframe = 0;
         ((NETCDF_INFO *)(fstat[0]))->nframe = 0;
       }

       for(j=0; j < dim; j++) *(avg[i] + j) = (float)(*(avgt[i] + j));

       writef(fstat[0], avg[i], dim);

       if(form != 4) close_file(fstat[0], ifl);
       else netcdf_close((NETCDF_INFO *)(fstat[0]));

   }

   for(i=0; i < 4; i++) {free(avg[i]); free(avgt[i]);}
   free(specwt);

   if(chrfld) free(chrfld);
   for(i=0; i<NUMLINES; i++) free(ihead[i]);


   return;

}
