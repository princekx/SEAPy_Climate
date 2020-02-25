#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "pp.h"
#include "mem_er.h"
#include "grid.h"

/* function to read PP header */

extern int init;
extern float period;
extern float *abuf;
extern void *databuf;

extern FILE *finit;
extern int nhp, shp;
extern int std_x, std_y;

PP *pph=NULL;

void grid_trans(GRID * , float * , int * , int * );

int pp_header(GRID *gr, FILE *fpp, int *tl, int *gof, int ftm, int pr)
{

  int i;
  int iin='n';

  int ix=0, iy=0;

  long int ss;

  static int ftyp=0;
  static int imsg=0;
  static int fr=0;

  float *xgtmp=NULL;

  if(ftm > 0){

     pph = (PP *)malloc_initl(sizeof(PP));
     mem_er((pph == NULL) ? 0 : 1, sizeof(PP));


  }

  else if(ftm < 0) { free(pph); free(databuf); imsg = 0; return 0;};

  fread(&ss,4,1,fpp);

  fread(pph,sizeof(PP),1,fpp);

  fread(&ss,4,1,fpp);

  if(ferror(fpp)){ 

     printf("****ERROR****, an error has been encountered reading PP format header.\n\n");
     return 0;

  }

  else if(feof(fpp)){

     printf("****WARNING****, an EOF has been encountered.\n\n");

     return 0;

  }

  if(pr > 0){

      printf("Field date/time: %d %d %d %d\n",
              pph->lbyr,pph->lbmon,pph->lbdat,pph->lbhr);

      printf("Grid code type: %d\n", pph->lbcode);
      printf("Field code type: %d\n", pph->lbfc);
      printf("Level type code: %d\n", pph->lbvc);
      printf("Level value: %f\n", pph->blev);

  }

  if(!fr){
     if(pph->lbcode == 1 || pph->lbcode == 2){}
     else{

        printf("****ERROR****, data may not be on a regular lat-long grid\r\n"
            "               do you want to continue, 'y' or 'n'.      \n\n");

        if(init) {
          fscanf(finit,"\n");
          if(getc(finit) == 'n') exit(1);
        }
        else {
          scanf("\n");
          iin = getchar();
          fprintf(finit, "%c\n", iin);
          if(iin == 'n') exit(1);
        }

     }

     fr = 1;

  } 

  if(!ftyp) ftyp = pph->lbfc;
  else if(ftyp != pph->lbfc && !imsg){

     printf("****WARNING****, encountered field of different type\n");
     imsg = 1;


  }


  if(pph->lbpack){

     printf("****ERROR****, PP field is in packed format, no un-packing available yet.\n");
     exit(1);

  }

  if(ftm){

     nhp = shp = 0;

     ix = gr->ix = pph->lbnpt;

     iy = gr->iy = pph->lbrow;

     std_x = pph->lbnpt;
     std_y = pph->lbrow;

     printf("the current grid dimensions are %d * %d \n\n", gr->ix, gr->iy);


/* assign memory for grid */

     gr->xgrid = (float * )calloc(gr->ix + 1, sizeof(float));
     mem_er((gr->xgrid == NULL) ? 0 :1, (gr->ix + 1) * sizeof(float));

     xgtmp = (float * )calloc(gr->ix, sizeof(float));
     mem_er((xgtmp == NULL) ? 0 :1, gr->ix * sizeof(float));


     gr->ygrid = (float * )calloc(gr->iy, sizeof(float));
     mem_er((gr->ygrid == NULL) ? 0 : 1, gr->iy * sizeof(float));

/* assign memory for reading field data */

     abuf = (float *)calloc(pph->lbnpt * pph->lbrow, sizeof(float));
     mem_er((abuf == NULL) ? 0 : 1, pph->lbnpt * pph->lbrow * sizeof(float));

     databuf = (void *)malloc_initl(pph->lbnpt * pph->lbrow * sizeof(float));
     mem_er((databuf == NULL) ? 0 : 1, pph->lbnpt * pph->lbrow * sizeof(float));

     *(xgtmp) = pph->bzx + pph->bdx;


     for(i=1; i < gr->ix; i++){

         *(xgtmp + i) =  *(xgtmp + i - 1) + pph->bdx;

     }



     if(pph->bdy > 0.){

        *(gr->ygrid) = pph->bzy + pph->bdy;

        for(i=1; i < gr->iy; i++){

            *(gr->ygrid + i) =  *(gr->ygrid + i - 1) + pph->bdy;

        }

      }

      else {

        *(gr->ygrid + gr->iy - 1) = pph->bzy + pph->bdy;

        for(i=gr->iy - 2; i >= 0; i--){

            *(gr->ygrid + i) =  *(gr->ygrid + i + 1) + pph->bdy;

        }


      }

      grid_trans(gr, xgtmp, tl, gof);


      free(xgtmp);


  }

/*  else{
    if(ix != pph->lbnpt || iy != pph->lbrow){
       printf("****ERROR****, fields have different grid dimensions.       \r\n"
              "               All fields must have the same X-Y dimensions.\n\n");
       exit(1);
    }

  } */



  return 1;
 

}
