#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "grid.h"
#include "m_values.h"
#include "mem_er.h"


/* Discrete 2D fast cosine transform based on Jain, Fundamentals of Digital
   Image Processing.                                                     */
   
void dfct(complex * , int , int , int , int , int );
   
extern complex *yct;
extern complex *wct;

extern complex *zsave1, *zsave2;    
extern complex *yt;

extern GRID *gr;

extern int *njx, *njy;

extern double *dtrgx, *dtrgy;

extern double *dtrig;

extern int *nj;

void dfct2d(complex *ac, int nx, int ny, int nexpx, int nexpy, int inv, int itrig, int ifftyp)
{
    int i=0, j=0;
    int itt=0;

    static complex *zs1x=NULL, *zs2x=NULL;
    static complex *zs1y=NULL, *zs2y=NULL;
    static complex *ytx=NULL, *yty=NULL;
    static complex *yctx=NULL, *ycty=NULL;
    static complex *wct1x=NULL, *wct2x=NULL;
    static complex *wct1y=NULL, *wct2y=NULL;
    
    static complex *yy=NULL;
    
    complex *ctp=NULL;
    
    if(itrig > 0) itt = 1;
    
    if(itrig < 0){
       zsave1 = zs1x;
       zsave2 = zs2x;
       yt = ytx;
       wct = wct1x;
       yct = yctx;
       dfct(NULL, 0, 0, 1, -1, 0); 
       wct = wct2x;
       dfct(NULL, 0, 0, 1, -1, 0);
       
       zsave1 = zs1y;
       zsave2 = zs2y;
       yt = yty;
       wct = wct1y;
       yct = ycty;
       dfct(NULL, 0, 0, 1, -1, 0); 
       wct = wct2y;
       dfct(NULL, 0, 0, 1, -1, 0);
       
       free(yy);
       
       return;        
    }
    
    if(inv != 1){

      if(inv != -1){

         printf("****ERROR****, transform direction must be either foreward (1)\r\n"
                "               or backward (-1).                              \n\n");
         exit(1);

      }

    }
    
    if(!yy){
       yy = (complex *)calloc(ny, sizeof(complex));
       mem_er((yy == NULL) ? 0 : 1, ny * sizeof(complex));
    }
    
    if(inv == 1) wct = wct1x;
    else wct = wct2x;
    yct = yctx;    
    
    if(!ifftyp) {
       zsave1 = zs1x; 
       zsave2 = zs2x;      
       yt = ytx;
    }
    else {
       dtrig = dtrgx;
       nj = njx;
    }
    
    for(i=0; i < ny; i++){
       
        dfct(ac + i * nx, nx, nexpx, inv, itrig, ifftyp);
             
        if(itrig > 0){
	   zs1x = zsave1;
	   zs2x = zsave2;
	   ytx = yt;
	   yctx = yct;
	   if(inv == 1) wct1x = wct; 
	   else wct2x = wct;
	   itrig = 0; 
	}	   
    }    
   
    if(inv == 1) wct = wct1y;
    else wct = wct2y;
    yct = ycty;
 
    if(!ifftyp) {      
       zsave1 = zs1y;
       zsave2 = zs2y;
       yt = yty;
    }
    
    else {
       dtrig = dtrgy;
       nj = njy;
    }

    if(itt) itrig = 1;
       
    for(i=0; i < nx; i++){
       
        for(j=0; j < ny; j++) {ctp = ac + j * nx + i; comp(ctp->real, 0.0, (yy + j));}
	
	dfct(yy, ny, nexpy, inv, itrig, ifftyp);
	   
	if(itrig > 0){
	   zs1y = zsave1;
	   zs2y = zsave2;
	   yty = yt;
	   ycty = yct;
	   if(inv == 1) wct1y = wct;
	   else wct2y = wct; 
	   itrig = 0;  
	}
	   
	for(j=0; j < ny; j++) {ctp = yy + j; comp(ctp->real, 0.0, (ac + j * nx + i));}

    }
    
    return;
}
