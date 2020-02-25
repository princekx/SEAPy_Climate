#include <Stdio.h>
#include <stdlib.h>
#include "st_obj.h"
#include "st_fo.h"
#include "grid.h"
#include "bisp.h"
#include "m_values.h"

/* function to add an additional field value onto the feature points */

#ifdef  REGULAR

#ifdef  NOUNDERSCORE

void bisp(double * , double * , double * , int * , double * , int * , 
          double * , int * , double * , double * , int * );
	  
#else

void bisp_(double * , double * , double * , int * , double * , int * , 
           double * , int * , double * , double * , int * );

#endif

#endif

extern int tf;

extern int x1u, x2u, y1u, y2u;
extern GRID *gr;

void adjust_fpt(struct frame_objs *ff, struct sp_dat *ct, float *apl, int aob)
{

   int i, j;
   int idif=0;

   double zz, dd[4];
   double xa, ya;

   struct object *ob=NULL; 
   struct point *ptt=NULL;
   struct feature_pts *fpts=NULL; 
    
   for(i=0; i < ff->obj_num; i++){

       ob = ff->objs + i;
       
       if(aob == 'y'){
       
          for(j=0; j < ob->point_num; j++){ 
	  
	      ptt=(ob->pt)+j;
	      	          	      
	      ptt->val += *(apl + (ptt->y + y1u - 2) * gr->ix + ptt->x + x1u - 2);  
	  
	  }
       
       }

       for(j=0; j < ob->fet->feature_num; j++){
       
          fpts = ob->fet->fpt + j;
	  
          if(fpts->str > CHECK_PT) {
	  	  
             if(tf == 3){
	     
	        fpts->str += *(apl + ((fpts->y).ixy +y1u-2) * gr->ix + ((fpts->x).ixy)+x1u-2);
	  
	     }
             else {
  
		
                if(ct->sm_type){
	           ya = (fpts->x).xy;
                   xa = (fpts->y).xy;
		   ya *= FP_PI;
		   xa = FP_PI2 - xa * FP_PI;
		   if(xa < 0.) xa = 0.0;
		   else if(xa > FPI) xa = FPI;		
		}
		
		else {
	           xa = (fpts->x).xy;
                   ya = (fpts->y).xy;		
		}

#ifdef REGULAR		  

#ifdef  NOUNDERSCORE

                bisp(&zz, dd, ct->tx, &ct->nx, ct->ty, &ct->ny, ct->c, &ct->ncof, &xa, &ya, &idif);

#else

                bisp_(&zz, dd, ct->tx, &ct->nx, ct->ty, &ct->ny, ct->c, &ct->ncof, &xa, &ya, &idif);

#endif

                fpts->str += zz;

#else

                printf("***ERROR***, surface fitting impossible unless correct libraries\r\n"
                       "             are linked, see compilation options.               \n\n");
                exit(1);

#endif
		
             }
	    
          } 

       }
	
   }
      
   return;
   
}
