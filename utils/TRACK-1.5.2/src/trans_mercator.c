#include <Stdio.h>
#include <Math.h>
#include "m_values.h"

/* Transverse Mercator projection */

void trans_mercator(float *gridx, float *gridy, float k0, float lng0, float lat0, int itp)


{

   float lng, lat;
   double *sthet, *cthet, *sphi, *cphi;
   double tanthet;
   double b, d;

   lng = *gridx;
   lat = *gridy;

   if(itp){                                  /* forward */

       lng *= FP_PI;                         /* convert to radians */
       lng -= lng0;
       lat *= FP_PI;

       sincos((double)lat, sthet, cthet);
       sincos((double)lng, sphi, cphi);

       tanthet = *sthet / *cthet;

       b = *cthet * (*sphi);

       *gridx = 0.5 * k0 * (float) log((1.0 + b) / (1.0 - b));
       *gridy = k0 * ((float) atan(tanthet / *cphi) - lat0);


   }

   else {

       b = lng / k0;
       d = (lat / k0) + lat0;

       *gridy = (float) asin(sin(d) / cosh(b)) / FP_PI;
       *gridx = lng0 + (float) atan(sinh(b) / cos(d)) / FP_PI;

   }

   return;

}
