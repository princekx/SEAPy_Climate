#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "m_values.h"

#define  TOLSP  0.001
#define  TOLRHO 1.0e-6
#define  TOLPOL 1.0e-5

double invtan(double , double );

extern float period;

/* subroutine to perform forward and backward conic projections. */


int conic(float *lat, float *lng, float lat1, float lng1, float sp1, float sp2, int itp, int prty, int *ift)

{
   int icon=0;
   
   int onesp=0;
   
   static int cntset=0;

   double x=0., y=0.;

   double xt=0.0, yt=0.0;

   double dl=0.0;
   double rho=0.0;
   double tht=0.;
   
   double nearpole=0;
   
   static double rho0=0.0, cc=0.0, nn=0.0;

   double s1, s2, s3, c1, c2, c3, slng, clng;
   
/* compute constants on first use */

   if(!(*ift) && !cntset) {
      printf("****WARNING****, constants for projection are not set, they will be recalculated.\n\n");
      *ift = 1;
   }

   if(*ift){
   
      if(fabs(sp1 - sp2) < TOLSP) onesp = 1;
      
      sp1 *= FP_PI;  
      sp2 *= FP_PI;
      
      sincos(sp1, &s1, &c1);
      sincos(sp2, &s2, &c2);
      
      
      lng1 *= FP_PI;
      lat1 *= FP_PI; 
      sincos(lat1, &s3, &c3);
      
      switch(prty){
          case 1:                     /* Albers Equal Area */
	    if(onesp) nn = s1;
            else nn = 0.5 * (s1 + s2);
	    cc = c1 * c1 + 2.0 * nn * s1;
            rho0 = sqrt(cc - 2.0 * nn * s3) / nn;
	    break;
	  case 2:                     /* Lambert Conformal */
	    if(onesp) {nn = s1;}
	    else nn = log(c1 / c2) / log(tan(FP_PI4 + 0.5 * sp2) / tan(FP_PI4 + 0.5 * sp1));
	    cc = c1 * pow(tan(FP_PI4 + 0.5 * sp1), nn) / nn;
	    nearpole = cos(FP_PI2 - fabs(lat1));
	    if(1.0 - nearpole < TOLPOL) rho0 = 0.0;
	    else rho0 = cc / pow(tan(FP_PI4 + 0.5 * lat1), nn);
	    break;
          case 3:                     /* Equidistant Conic */
	    if(onesp) nn = s1; 
	    else nn = (c1 - c2) / (sp2 - sp1);
	    cc= sp1 + (c1 / nn);
	    rho0 = cc - lat1;            
            break;
          default:
            printf("***ERROR***, projection identifier not known\n");
            exit(1);
      }   
      
      *ift = 0;
      cntset = 1;
      
      lng1 /= FP_PI;
      lat1 /= FP_PI;  

   }

   xt = *lng;
   yt = *lat;

   if(itp){

      dl = xt - lng1;
      if(dl > 180.0) dl -= 360.0;
      else if (dl < -180.0) dl += 360;
      xt *= FP_PI;
      yt *= FP_PI; 
      lng1 *= FP_PI;
      dl = nn * dl * FP_PI;
      sincos(yt, &s1, &c1);

      sincos(dl, &slng, &clng);
      
      switch(prty){
          case 1:
             rho = sqrt(cc - 2.0 * nn * s1) / nn;
             break;
          case 2:
             rho = cc / pow(tan(FP_PI4 + 0.5 * yt), nn);
             break;
          case 3:
             rho = cc - yt;
             break;
          default:
             printf("***ERROR***, projection identifier not known\n");
             exit(1);
      }

      x = rho * slng;
      y = rho0 - rho * clng;
      icon = 1;
   }

   else {

      dl = rho0 - yt;
      rho = sqrt(xt *xt + dl * dl);
      if(nn < 0.0){dl *= -1.0; xt *= -1.0; yt *= -1.0;} 
      tht = atan2(xt,dl);
      
      switch(prty){
          case 1:
             y = asin((cc - rho * rho * nn * nn) / (2.0 * nn));
             break;
          case 2:
             rho *= (nn != 0.0) ? ((nn > 0.0) ? 1.0 : -1.0) : 0.0; 
	     if(fabs(rho) < TOLRHO){
	        y = (nn > 0.0) ? FP_PI2 : -FP_PI2;
	     }
	     else { 
                y = 2.0 * invtan(pow(cc / rho, 1.0 / nn), 1.0) - FP_PI2;
             } 
             break;
          case 3:
             rho *= (nn != 0.0) ? ((nn > 0.0) ? 1.0 : -1.0) : 0.0;
	     y = (cc - rho);
             break;
          default:
             printf("***ERROR***, projection identifier not known\n");
             exit(1);
      }

      x = lng1 + (tht / (nn * FP_PI));
      
      y /= FP_PI;

      if(fabs(y) < 1.0e-5) y = 0.;
      if(x > 360.) x -= period;
      if (x < 0.0) x += period;
      icon = 1;

   }
   

   *lng = x;
   *lat = y;
   
   return icon;

}
