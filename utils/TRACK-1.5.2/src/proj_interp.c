#include <Stdio.h>
#include <stdlib.h>

#ifndef  REGULAR

void proj_interp()

{

   printf("***error***, surface fitting impossible unless correct libraries\r\n"
          "             are linked, see compilation options\n");
   exit(1);

}

#else

#include <Math.h>
#include <stdarg.h>
#include "grid.h"
#include "reg_dat.h"
#include "sphery_dat.h"
#include "bisp.h"
#include "proj.h"
#include "m_values.h"


/* function to interpolate field onto new grid */

void surfit(double * , int , int , ... );

#ifdef  NOUNDERSCORE

void bisp(double * , double * , double * , int * , double * , int * , 
           double * , int * , double * , double * , int * );

#else

void bisp_(double * , double * , double * , int * , double * , int * , 
           double * , int * , double * , double * , int * );

#endif

extern int x1u, x2u, y1u, y2u;
extern float xmn, ymn, xmx, ymx;

extern GRID *gr, *gt;

extern float *ap;


void proj_interp(float *ap1, double *s, int fb, int smty, ... )

{

    va_list sptr;

    int i, j, idif=0;
    int iret=0;
    int ift=0;
    
    float xx, yy, yt;

    float *apt=NULL;
    float *az=NULL;

    double zz, dd[4];
    double xa, ya;

    GRID *gtemp=NULL;

    struct sp_dat *ct=NULL;
    struct rspline *rt=NULL;
    struct savedat *st=NULL;
    struct sspline *ss=NULL;
    struct savedat_sphy *ssd=NULL;

    va_start(sptr, smty);

    x1u = gr->ox1u;
    x2u = gr->ox2u;
    y1u = gr->oy1u;
    y2u = gr->oy2u;

    xmn = *(gt->xgrid + x1u-1);
    ymn = *(gt->ygrid + y1u-1);
    xmx = *(gt->xgrid + x2u-1);
    ymx = *(gt->ygrid + y2u-1);

    gtemp = gr;
    gr = gt;

    apt = ap;
    ap = ap1;

    ct = va_arg(sptr, struct sp_dat * );

    if(smty) {
       ss = va_arg(sptr, struct sspline *);
       ssd = va_arg(sptr, struct savedat_sphy *);
       surfit(s, fb, smty, ct, ss, ssd);
    }
    else {

       rt = va_arg(sptr, struct rspline *);
       st = va_arg(sptr, struct savedat *);
       surfit(s, fb, smty, ct, rt, st);

    }

 
    gr = gtemp;
    ap = apt;
    apt = NULL;

    if(fb >= 0){

       for(i=0; i < gr->iy; i++){

          yt = yy = *(gr->ygrid + i);
          for(j=0; j < gr->ix; j++){

              xx = *(gr->xgrid + j);

              az = ap + i * gr->ix + j;
	      
	      switch(gr->prgr){
                 case 1:
	            iret = azimuthal(&yy, &xx, gr->alat, gr->alng, 0, gr->prty);
		    break;
		 case 2:
		    iret = conic(&yy, &xx, gr->alat, gr->alng, gr->sp1, gr->sp2, 0, gr->prty, &ift);
		    break;
                 case 3:
                    iret = rotated_cyln(&yy, &xx, gr->alat, gr->alng, 0, &ift);
                    break;
                 default:
                    printf("****ERROR****, unknown projection type.\n\n");
                    exit(1); 
	      }


/*              if(azimuthal(&yy, &xx, gr->alat, gr->alng, 0, gr->prty) && 
                 (xmx - xx) * (xx - xmn) >= 0.0 && (ymx - yy) * (yy - ymn) >= 0.0) { */
		 
	      if(iret && (xmx - xx) * (xx - xmn) >= 0.0 && (ymx - yy) * (yy - ymn) >= 0.0) {

                 xa = xx;
                 ya = yy;

                 if(smty) {

                    xa *= FP_PI;
                    ya = FP_PI2 - ya * FP_PI;

#ifdef  NOUNDERSCORE

                    bisp(&zz, dd, ct->tx, &ct->nx, ct->ty, &ct->ny, ct->c, &ct->ncof, &ya, &xa, &idif);

#else

                    bisp_(&zz, dd, ct->tx, &ct->nx, ct->ty, &ct->ny, ct->c, &ct->ncof, &ya, &xa, &idif);

#endif


                 }

                 else

#ifdef  NOUNDERSCORE 

                   bisp(&zz, dd,ct->tx, &ct->nx, ct->ty, &ct->ny, ct->c, &ct->ncof, &xa, &ya, &idif);

#else

                   bisp_(&zz, dd,ct->tx, &ct->nx, ct->ty, &ct->ny, ct->c, &ct->ncof, &xa, &ya, &idif);

#endif


                 *az = (float) zz;


              }

              else *az = 0.0; 

              yy = yt;

           }

       }
                 
    }


 
    x1u = 1;
    x2u = gr->ix;
    y1u = 1;
    y2u = gr->iy;

    xmn = *(gr->xgrid + x1u-1);
    ymn = *(gr->ygrid + y1u-1);
    xmx = *(gr->xgrid + x2u-1);
    ymx = *(gr->ygrid + y2u-1);

    va_end(sptr);

    return;

}

#endif
