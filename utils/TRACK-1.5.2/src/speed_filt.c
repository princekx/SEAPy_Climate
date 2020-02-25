#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "grid.h"
#include "splice.h"
#include "st_fo.h"
#include "geo_values.h"

#define  SPDLRG   1.0e+6
#define  SPDSML  -1.0e+6

/* function to filter tracks according to their propogation speed */

float measure(struct feature_pts * , struct feature_pts * );

extern int x1u, y1u, x2u, y2u;
extern GRID *gr;

extern int nf;

void speed_filt(struct tot_tr *all_tr, int trackn, int *totf)

{

    int i=0, j=0;
    int isc=0;
    int spty=0;
    int inreg=1;
    int inrg=0;
    int tn=0;
    int ttint=0;
    int ireg=0;
    int igm=0;
    
    float spd=0.0;
    float tstep=0.0;
    float sc=1.0;
    float spth1=0.0, spth2=0.0;
    float lon1=0.0, lon2=0.0, lat1=0.0, lat2=0.0;
    float spext=0.0;
    
    struct tot_tr *altr;
    struct fet_pt_tr *at, *at1;
    struct feature_pts f1, f2;
    
    printf("What is the time step in the required units?\n\n");
    scanf("%f", &tstep);
    
    printf("do you want to scale the phase speeds, 'y' or 'n'\r\n"
           "e.g. speed relative to the Earths surface, scale with radius of the Earth\n");

    scanf("\n");
    if(getchar() == 'y'){

      printf("What type of value do you want      \r\n"
             "Input '0' for user value            \r\n"
             "Input '1' for Earths radius in Km   \r\n");
      scanf("%d", &isc);
      
      if(isc < 0 || isc > 1){
         printf("****ERROR****, incorrect input for scaling choice.\n\n");
	 exit(1);
      }

      if(isc == 0 ) {

         printf("input scaling\n");
         scanf("%f", &sc);

      }

      else if (isc == 1) sc = EARTH_RADIUS;

    }
    
    printf("Do you want, the systems must attain speed in range, input '0'\r\n"
           "             system minimum speed must be in range,  input '1'\r\n"
           "             system maximum speed must be in range,  input '2'\n\n");

    scanf("%d", &spty);

    printf("What are the lower and upper speed thresholds?\n\n");
    scanf("%f %f", &spth1, &spth2);
    
    printf("Restrict to a Lon-Lat or sub region of domain, 'y' or 'n'.\n\n");
    scanf("\n");
    if(getchar() == 'y'){
       ireg = 1;
       printf("Input lon1, lon2, lat1, lat2 or x1, x2, y1, y2\n\n");
       scanf("%f %f %f %f", &lon1, &lon2, &lat1, &lat2);
       if(lon1 > lon2) igm = 1;
    }
    
    *totf = 0;
    
    for(i=0; i < trackn; i++){
    
       altr = all_tr + i;
       at = altr->trpt;
       at1 = at + 1;
       
       inrg = 0;
       
       if(spty == 1) spext = SPDLRG;
       else if(spty == 2) spext = SPDSML;
       
       for(j=1; j<altr->num; j++){
	  
	  if(ireg){
	     inreg = 0;
	     if(igm){
	        if((lat2 - at->yf) * (at->yf - lat1) >= 0.0 &&
		  ((lon2 - at->xf) * (at->xf - 0.0) >= 0.0 ||
		   (360.0 - at->xf) * (at->xf - lon1) >= 0.0)) inreg = 1;
                if((lat2 - at1->yf) * (at1->yf - lat1) >= 0.0 &&
		  ((lon2 - at1->xf) * (at1->xf - 0.0) >= 0.0 ||
		   (360.0 - at1->xf) * (at1->xf - lon1) >= 0.0)) inreg = 1;		 	     
	     }
	     else {
	        if((lon2 - at->xf) * (at->xf - lon1) >= 0.0 && 
	            (lat2 - at->yf) * (at->yf - lat1) >= 0.0) inreg = 1;
	        if((lon2 - at1->xf) * (at1->xf - lon1) >= 0.0 && 
		    (lat2 - at1->yf) * (at1->yf - lat1) >= 0.0) inreg = 1;
             }
	  }
	  
	   
          if(inreg) {
	  
             (f1.x).xy = at->xf;
             (f1.y).xy = at->yf;

             (f2.x).xy = at1->xf;
             (f2.y).xy = at1->yf;
	     ttint = (at->nfm + 1);
	  

             spd = sc * measure(&f1, &f2)/(ttint * tstep); 
	  
	     switch(spty){
	         case 0:
	             if((spth2 - spd) * (spd - spth1) >= 0.0) inrg=1;
	             break;
                 case 1:
	             if(spd < spext) spext = spd;
	             break;
	         case 2:
	             if(spd > spext) spext = spd;	      
	             break;
	  
	     }
	  
	  }
	  
	  	  
	  ++at;
          ++at1;
	  
        }
	
	if(spty > 0){
	  if((spth2 - spext) * (spext - spth1) >= 0.0) inrg = 1;
	}
	
	if(!inrg){
	
	   if(nf){
             for(j=0; j < altr->num; j++) free((altr->trpt + j)->add_fld);
	   }
           altr->num = 0;
           free(altr->trpt);
           altr->trpt = NULL;	   
	
	}
	
	else ++tn;
	
	*totf += altr->num;
     
    }

    printf("***INFORMATION***, current number of tracks is %d\n\n", tn);

    return;
}

