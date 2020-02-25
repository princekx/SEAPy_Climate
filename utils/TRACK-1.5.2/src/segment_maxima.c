#include <Stdio.h>
#include <stdlib.h>
#include <mem_er.h>


/* function to detect maxima in a digital image */

#define  LESS(A, B, C) (int)((A <= B) ? (C + 1) : C)

int xdl, ydl;

void collapse(int * , int * , int , int );

int *segment_maxima(int *oba, int xd, int yd)
{

    int i, j;
    int count=0, dim=xd*yd;
    int *ac=NULL, *aa=NULL, *ab=NULL;
    int *at=NULL, *att=NULL;

    xdl = xd;
    ydl = yd;

    at = (int *)calloc(dim, sizeof(int));
    mem_er((at == NULL) ? 0 : 1, dim * sizeof(int));

    for(i=0; i < dim; i++) *(at+i) = *(oba+i);

    for(i=yd-2; i>0; i--){

       ac = oba + (i+1) * xd + 1;
       aa = oba + i * xd + 1;
       ab = oba + (i-1) * xd + 1;

       att = at + i * xd + 1;

       for(j=1; j < xd-1; j++){

           count = 0;

           if(*att){

             count = LESS(*(aa-1), *aa, count); 
             count = LESS(*(aa+1), *aa, count);
             count = LESS(*ab, *aa, count);
             count = LESS(*(ab - 1), *aa, count);
             count = LESS(*(ab + 1), *aa, count);
             count = LESS(*ac, *aa, count);
             count = LESS(*(ac - 1), *aa, count);
             count = LESS(*(ac + 1), *aa, count); 

             if(count != 8) collapse(aa, att, j, i);            

           }

           ++aa; ++ab; ++ac, ++att;




       }


    }


    free(oba);

    return at;

}


void collapse(int *aa, int *at, int j, int i)

{

    int *ab, *ac;
    int *at1, *at2;
    int att = *at;

    *at  = 0;

    if(i == 1 || i == ydl-1) return;
    else if(j == 1 || j == xdl-1) return;

    if(att){

      if(*aa == *(aa - 1)) collapse(aa - 1, at - 1, j-1, i);
      if(*aa == *(aa + 1)) collapse(aa + 1, at + 1, j+1, i);
      ab = aa - xdl;
      at1 = at - xdl;
      ac = aa + xdl;
      at2 = at + xdl;
      if(*aa == *ac) collapse(ac, at2, j, i-1);
      if(*aa == *(ac - 1)) collapse(ac - 1, at2 - 1, j-1, i+1);
      if(*aa == *(ac + 1)) collapse(ac + 1, at2 + 1, j+1, i+1);
      if(*aa == *ab) collapse(ab, at1, j, i+1);
      if(*aa == *(ab - 1)) collapse(ab - 1, at1 - 1, j-1, i-1);
      if(*aa == *(ab + 1)) collapse(ab + 1, at1 + 1, j+1, i-1);

    }

    return;

}
