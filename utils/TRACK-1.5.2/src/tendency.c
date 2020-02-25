#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Math.h>
#include <sys/types.h>
#include "grid.h"
#include "files_out.h"
#include "file_handle.h"
#include "mem_er.h"
#include "utf.h"
#include "pp.h"


/* function to compute the field tendencies, i.e. X(t+h) - X(t) */


extern GRID *gr;
extern int frnum;
extern int form;
extern int eqsw;
extern int utfv;
extern int i_utf4x, i_utf4y;
extern int std_x, std_y;

extern float *abuf;

extern PP *pph;

extern char *fext;

extern int iext;

extern char *chrfld;

extern char *ihead[NUMLINES];
extern int nlines;

int rf(FILE * ,FILE * , int , int );
void writef(FILE * , float * , int );

void tendency(FILE *fdat, int fr1, int frl)
{

   int i;
   int idt;
   int nf=0, nn=0;
   int dim=0;

   float *buf2=NULL, *atemp=NULL;
   float ftemp;

   off_t chrnum=0, place1, place2;

   FILE *ftend=NULL;

   char tend[MAXCHR];

/* assign space for header strings */

   nlines = 0;

   for(i=0; i < NUMLINES; i++){
       ihead[i] = (char *)calloc(MAXCHR, sizeof(char));
       mem_er((ihead[i] == NULL) ? 0 : 1, MAXCHR * sizeof(char));
   }


   if(form < 2) dim = std_x * std_y;
   else if (form == 2){
      if(utfv == 3) dim=gr->ix * ((eqsw) ? gr->iy + 1 : gr->iy);
      else if(utfv == 4) dim = i_utf4x * i_utf4y;

   }

   else if(form == 3) dim = pph->lbrow * pph->lbnpt;

   buf2 = (float *)calloc(dim, sizeof(float));
   mem_er((buf2 == NULL) ? 0 : 1, dim * sizeof(float)); 

   strncpy(tend, TENDENCY, MAXCHR);
   if(iext) strcpy(strstr(tend, EXTENSION), fext);

   ftend = open_file(tend, "w");

   fprintf(ftend, "%10c\n", ' ');

re_try_idt:


   printf("What is the number of time steps to difference over?\n\n");
   scanf("%d", &idt);

   if(idt <= 0){
      printf("****ERROR****, must be a positive number of frames to difference over.\r\n"
             "               input number of frames again.                          \n\n");
      goto re_try_idt;

   }

   printf("***INFORMATON***, Computing Tendency of Data....... \r\n\n"
          "                                                    \r\n\n"
          " Please Wait, this may take some time ..............  \n\n");

   while(nf < frl){


      if(!nf){

         atemp = abuf;
         abuf = buf2;

         place1 = ftello(fdat);

         if(rf(fdat, ftend, dim, 1)) break;


         if(fr1 == 1){++nf; ++nn;}

         place2 = ftello(fdat);


         chrnum = place2 - place1;

         if(fr1 > 1){

              fseeko(fdat, (fr1-2)*chrnum, ORIGIN);

              rf(fdat, ftend, dim, 0);

              nf = fr1;
              ++nn;

         }

         abuf = atemp;
         atemp = NULL;

      }


      if(nf > frl) break;


      fseeko(fdat, (idt - 1)*chrnum, ORIGIN);
      nlines = 0;

      if(rf(fdat, ftend, dim, 1))break;

      for(i=0; i<nlines; i++) fprintf(ftend, "%s", ihead[i]);

      nf += idt;

      ++nn;
      for(i=0; i<dim; i++) {

          ftemp = *(abuf + i);
          *(abuf + i) -= *(buf2 + i); 

          *(buf2 + i) = ftemp;

      }
 
      writef(ftend, abuf, dim);



   }

   fprintf(ftend, "%s\n", EOD);

   fseeko(ftend, (off_t)0, FSTART);

   fprintf(ftend, "%4d %4d", nn, idt);


   close_file(ftend, tend);

   free(buf2);
   if(chrfld) free(chrfld);

   printf("\n");

   for(i=0; i<NUMLINES; i++) free(ihead[i]);

   return;

}
