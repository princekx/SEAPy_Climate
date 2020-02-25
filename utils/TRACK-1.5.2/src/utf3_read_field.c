#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "grid.h"
#include "utf.h"

#define  INV      -1.0

extern int tl, gof;
extern float offs, ilat;

extern GRID *gr;
extern int eqsw;
extern int whemi;

/* function to read in each field of data and apply an offset if required */

int utf3_read_field(float *ap, float *ata, float scl, FILE *fdatin, int sign, int as, int of, int hemi, int tscl)

{

    int i, j, k;
    int pp=0, plus;
    int nft=NFT, ad;
    int ic1, ic2;
    int iyh;

    float fmx, fmn, fsc;
    float scale;

    char i1, i2;
    char text[MAXCHR];

    float *at;
    float yy;

    if(gr->gcen){

      printf("****WARNING****, no field translation available for UTF's\n\n");
      gr->gcen = 0;

    }

    iyh = gr->iy / 2;

    if(whemi){

       if(whemi == 1) iyh = gr->iy;
       else iyh = 0;

    }

    if(!fgets(text, MAXCHR, fdatin)) return 1;
    if(strstr(text, EOD)) {
      printf("***WARNING***, encountered end of data statement\n\n");
      return 1;
    }

    if(REPORT) printf("\n%s\n", text);

    fgets(text, MAXCHR, fdatin);

    if(REPORT) printf("\n%s\n", text);

    for(i=0; i < MHD; i++) fscanf(fdatin, "%*d");

    fscanf(fdatin, "%e %e %e\n", &fmn, &fmx, &fsc);

    scale = (fmx - fmn)/ARNG;

    ad = 0;
    
    for(j=0; j < gr->iy; j++){

       pp = (gr->iy - j - 1) * gr->ix;
       yy = *(gr->ygrid + j) * (-1.0);

       if(j == iyh && eqsw) {

         for(i=0; i < gr->ix; i++){
            fscanf(fdatin, "%*c %*c");
            if((++ad) == nft){fscanf(fdatin, "%*c"); ad = 0;}
         }

       }

       for(i=0; i< gr->ix; i++) {

          if(tl == 'n') plus = pp + i;

          else if(gof >= 0){

             k = i + gof;
             if(k >= gr->ix) k -= (gr->ix - 1);
             plus = pp + k;

          }

          else{

             k = i + gof;
             if(k < 0) k += (gr->ix - 1);
             plus = pp + k;

          }

          at = ap+plus;

          fscanf(fdatin, "%c %c", &i1, &i2);


          ic1 = (int)(uintptr_t)strchr(lkup, i1)-(int)(uintptr_t)lkup;
          ic2 = (int)(uintptr_t)strchr(lkup, i2)-(int)(uintptr_t)lkup;
          
          if(ic1 > CNUM || ic2 > CNUM){

            printf("***error***, illegal character in data field\n");
            exit(1);

          }

          *at = (float)(ic1*CNUM + ic2) * scale + fmn;                

          if(tscl == 'y') *at *= scl;

          if(as == 'y') *at = *at - *(ata+plus); 

          if(of == 'y') *at -= offs;

          *at *= sign;

          if(hemi == 'n' && yy > ilat) *at *= INV;
          else if(hemi == 's' && yy < ilat) *at *= INV;

          if((++ad) == nft){fscanf(fdatin, "%*c"); ad = 0;}

       }

       if(gof > 0) *(ap + pp) = *(ap + pp + gr->ix - 1);
       else if (gof < 0) *(ap + pp + gr->ix - 1) = *(ap + pp);
    }

    if(iyh == gr->iy && eqsw) {

      for(i=0; i < gr->ix; i++){
         fscanf(fdatin, "%*c %*c");

         if((++ad) == nft){fscanf(fdatin, "%*c"); ad = 0;}
      }

    }

    if(ad)fscanf(fdatin, "%*c");

    return 0;

}
