#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include "m_values.h"

double invtan(double , double );

extern float period;

/* subroutine to perform forward and backward azimuthal projections. */


int azimuthal(float *lat, float *lng, float lat1, float lng1, int itp, int prty)

{
   int icon=0;

   double x=0., y=0.;

   double xt, yt;

   double c=0., dl, cs, ss;
   double rho;
   double den;
   double kk=0.;
   double ac=0.;

   double s1, s2, c1, c2, slng, clng;

   xt = *lng;
   yt = *lat;

   lng1 *= FP_PI;
   lat1 *= FP_PI; 
   sincos(lat1, &s2, &c2);

   if(itp){

      xt *= FP_PI;
      yt *= FP_PI; 
      dl = xt - lng1;
      sincos(yt, &s1, &c1);

      sincos(dl, &slng, &clng);

      c = s1 * s2 + c1 * c2 * clng;

      if(c >= 0.) {

         x = c1 * slng;
         y = c2 * s1 - s2 * c1 * clng;

         switch(prty){
             case 1:
                kk = 1.0;
                break;
             case 2:
                kk = 2.0 / (1.0 + s1 * s2 + c1 * c2 * clng);
                break;
             case 3:
                kk = 1.0 / c;
                break;
             case 4:
                kk = sqrt(2.0 / (1.0 + s1 * s2 + c1 * c2 * clng));
                break;
             case 5:
                ac = acos(c);
                kk = ac / sin(ac);
                break;
             default:
                printf("***ERROR***, projection identifier not known\n");
                exit(1);
         }

         x *= kk;
         y *= kk;
         icon = 1;

      }

      else{ x = xt / FP_PI; y = yt / FP_PI;}


   }

   else {


      rho = sqrt(xt * xt + yt * yt);

      switch(prty){
          case 1:
             if(rho <= 1.0) {c = asin(rho); icon = 1;}
             break;
          case 2:
             if(rho <= 2.0) {c = 2.0 * atan(rho / 2.0); icon = 1;}
             break;
          case 3:
             c = atan(rho);
             icon = 1;
             break;
          case 4:
             if(rho <= SQ_R2) {c = 2.0 * asin(rho / 2.0); icon =1;}
             break;
          case 5:
             if(rho <= FP_PI2) {c = rho; icon = 1;}
             break;
          default:
             printf("***ERROR***, projection identifier not known\n");
             exit(1);
      }

      sincos(c, &ss, &cs);

      if(rho > 0.){

         if(icon){

            y = asin(cs * s2 + (yt * ss * c2 /rho));

            den = rho * c2 *cs - yt * s2 * ss;

            xt *= ss;

            x = invtan(xt, den);

            x += lng1;

         }

      }

      else { x = lng1; y = lat1; icon = 1;}

      x /= FP_PI;
      y /= FP_PI; 

      if(fabs(y) < 1.0e-5) y = 0.;
      if(x > 360.) x -= period;

   }

   *lng = x;
   *lat = y;


   return icon;

}



/* function to return the inverse tangent in the range 0 -- 2pi */

double  invtan(double yy, double xx)

{

    double th=0.;

    if(yy >= 0.)

       th = atan2(yy , xx);

    else if(yy < 0.)

       th = 2. * FPI + atan2(yy , xx);


    return th;

}
