#include <Stdio.h>
#include <Math.h>
#include <stdlib.h>
#include "m_values.h"
#include "vec.h"
#include "proj.h"

extern float period;
extern PROJ *ppm;

/* Rotated Cylindrical projection */

int rotated_cyln(float *lat, float *lng, float latp, float lngp, int itp, int *ift)


{

   int icon=0;

   static int cntset=0;
   
   float x=0., y=0.;
   float per2=period*0.5;
   
   double xt=0.0, yt=0.0;
   
   double sthp=0.0, cthp=0.0;
   double sphp=0.0, cphp=0.0;
   double sx=0.0, cx=0.0;
   double sy=0.0, cy=0.0;
   
   static VEC rotf[3]={{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
   static VEC rotb[3]={{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
   
   
   VEC ptt={0.0, 0.0, 0.0}, nptt={0.0, 0.0, 0.0};
   
/* compute constants on first use */

   if(!(*ift) && !cntset) {
      printf("****WARNING****, constants for projection are not set, they will be recalculated.\n\n");
      *ift = 1;
   }
   
   if(*ift){
      lngp *= FP_PI;
      latp *= FP_PI; 
      sincos(latp, &sthp, &cthp); 
      sincos(lngp, &sphp, &cphp); 
      
/* forward rotation matrix */
      
      rotf[0].x = cphp * cthp;
      rotf[0].y = -sphp * cthp;
      rotf[0].z = sthp;
      rotf[1].x = sphp;
      rotf[1].y = cphp;
      rotf[2].x = -cphp * sthp;
      rotf[2].y = sphp * sthp;
      rotf[2].z = cthp; 
      
/* backward rotation matrix, transpose of forward */
      
      rotb[0].x = rotf[0].x;
      rotb[0].y = rotf[1].x;
      rotb[0].z = rotf[2].x;
      rotb[1].x = rotf[0].y;
      rotb[1].y = rotf[1].y;
      rotb[1].z = rotf[2].y;
      rotb[2].x = rotf[0].z;
      rotb[2].z = rotf[2].z;   
      
      *ift = 0;
      cntset = 1;
      
   }

   if(itp){                                  /* forward */

      xt = *lng * FP_PI;
      yt = *lat * FP_PI;

      sincos(xt, &sx, &cx);
      sincos(yt, &sy, &cy);

      ptt.x = cx * cy;
      ptt.y = sx * cy;
      ptt.z = sy;
      
      nptt.x = dotp(&rotf[0], &ptt);
      nptt.y = dotp(&rotf[1], &ptt);
      nptt.z = dotp(&rotf[2], &ptt);    
  
      x = atan2(nptt.y, nptt.x) / FP_PI;
      if(x < 0.0) x = 360.0 + x;
      y = (FP_PI2 - acos(nptt.z)) / FP_PI;
      if(y > 90.0) y = 90.0;
      else if(y < -90.0) y = -90.0;
      
      if(x > per2) x -= period;
      
      if(ppm->prj2 != NULL) {
        (ppm->prj2)(&x, 'x', 1);
        (ppm->prj2)(&y, 'y', 1);
      }

   }

   else {
   
      x = *lng;
      y = *lat;

      if(ppm->prj2 != NULL) {
        (ppm->prj2)(&x, 'x', 0);
        (ppm->prj2)(&y, 'y', 0);  
      }

      x *= FP_PI;
      y *= FP_PI;

      sincos(x, &sx, &cx);
      sincos(y, &sy, &cy);

      ptt.x = cx * cy;
      ptt.y = sx * cy;
      ptt.z = sy;        

      nptt.x = dotp(&rotb[0], &ptt);
      nptt.y = dotp(&rotb[1], &ptt);
      nptt.z = dotp(&rotb[2], &ptt);
      
      x = atan2(nptt.y, nptt.x) / FP_PI;
      if(x < 0.0) x = 360.0 + x;
      y = (FP_PI2 - acos(nptt.z)) / FP_PI;
      if(y > 90.0) y = 90.0;
      else if(y < -90.0) y = -90.0;      

   }

   *lng = x;
   *lat = y;
   
   icon = 1;     

   return icon;

}
