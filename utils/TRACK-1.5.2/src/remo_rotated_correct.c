#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "grid.h"
#include "st_fo.h"
#include "st_obj.h"
#include "m_values.h"
#include "mem_er.h"


/* Routine to correct for rotated coordinate systems. e.g. when using data from a LAM */

extern float period;
extern GRID *gr;
extern int x1u, x2u, y1u, y2u;

extern int nf, nfld;
extern int *nfwpos;

void remo_rotated_correct(struct frame_objs *fo, float cset, int *tf, int tl, int frn)
{


    int i, j, k;
    int nf0=0, nfld0=0;
    
    float nplat=0.0, nplng=0.0;
    
    double xx=0.0, yy=0.0;
    double arg1, arg2;
    double nlats, nlatc, nlngs, nlngc;
    double s1, c1, s2, c2;
   
    struct frame_objs *ft=NULL;
    struct object *ob=NULL;
    struct feature_pts *fp=NULL;
    
/* set added field settings for one added field */

    nfld0 = nfld;
    nf0 = nf;
    
    ++nf;
    nfld += 3;
    nfwpos = (int *)realloc_n(nfwpos, nf*sizeof(int));
    mem_er((nfwpos == NULL) ? 0 : 1, nf * sizeof(int));
    *(nfwpos + nf - 1) = 1;
    
    printf("****INFORMATION****, assuming the poles have been rotated, correct back to geographical coordinates.\n\n");
    printf("What are the longitude, latitude location of the rotated pole, npphi, nptheta?\n\n");
    scanf("%f %f", &nplng, &nplat);
    
    printf("Has grid been translated in the data pre-processing, 'y' or 'n'?\n\n");
    scanf("\n");
    if(getchar() == 'y'){
       tl = 'y';
       printf("How much is the translation?\n\n");
       scanf("%f", &cset);
    }
    
    nplat *= FP_PI;
    nplng *= FP_PI;
    sincos(nplat, &nlats, &nlatc);
    sincos(nplng, &nlngs, &nlngc);

    for(i=0; i < frn; i++){

        ft = fo + i;

        for(j=0; j < ft->obj_num; j++){

            ob = (ft->objs) + j;

            if(!(ob->fet)) continue;
	    
	    for(k=0; k < ob->fet->feature_num; k++){

               fp = (ob->fet->fpt) + k;

               if(!nf0) fp->add_fld = NULL;

               fp->add_fld = (float *)realloc_n(fp->add_fld, nfld*sizeof(float));
               mem_er((fp->add_fld == NULL) ? 0 : 1, nfld*sizeof(float));
	       
	       if(*tf == 3){

                 xx = *(gr->xgrid + ((fp->x).ixy)+x1u-2);
                 yy = *(gr->ygrid + ((fp->y).ixy)+y1u-2);
		 
		 (fp->x).xy = xx;
		 (fp->y).xy = yy;

               }

               else {

                 xx = (fp->x).xy;
                 yy = (fp->y).xy;

               }	       

               if(gr->gcen == 1) {
	          printf("****ERROR****, grid transformation not yet supported for rotated grids.\n\n");
		  exit(1);
	       }
	       else if(gr->gcen == 2) xx -= gr->f_off;
	       
	       if(tl == 'y') xx += cset;
	       
	       if(xx > 180.0) xx -= 360.0;
	       else if(xx < -180.0) xx += 360.0;

	       xx *= FP_PI;
               sincos(xx, &s1, &c1);
	       yy *= FP_PI;
	       sincos(yy, &s2, &c2);
	       
               arg1 = nlats * s2 + nlatc * c2 * c1;	       
	       yy = asin(arg1) / FP_PI;	       
	       
	       arg1 = nlngs * (-nlats * c1 * c2 + nlatc * s2) - nlngc * s1 * c2;
	       arg2 = nlngc * (-nlats * c1 * c2 + nlatc * s2) + nlngs * s1 * c2;
	       xx = atan2(arg1, arg2) / FP_PI;

               if(xx < 0.0) xx += period;
	       else if(xx > period) xx -= period;
	       
               *(fp->add_fld + nfld0) = (fp->x).xy;
	       *(fp->add_fld + nfld0 + 1) = (fp->y).xy;
	       *(fp->add_fld + nfld0 + 2) = fp->str;
               
	       (fp->x).xy = xx;
	       (fp->y).xy = yy;

	    }
	    
	}

    }

    if(*tf == 3) *tf = 4;

    return;
   
}
