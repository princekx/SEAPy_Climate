#include <Stdio.h>
#include <Math.h>
#include "m_values.h"
#include "complex.h"

#define  TOLCOR  1.0e-6

/* rednoise with confidence cuttoff computed for 2 degree's of freedom,
   this needs to be ammended for other spectral estimators other than
   periodogram.                                                       */

#ifdef NOUNDERSCORE

void cdfchi(int * , double * , double * , double * , double * , int * , double * );

void rednoise(double var, double alph, double *redn, double *redc, double siglev, int nf)
{
    int i;
    int iw=2;
    int status;

    double a2, a22;
    double pp=1.0-siglev;
    double dof=2.0;
    double bound;

    a2 = alph * alph;
    a22 = 2.0 * alph;

    var /= nf;

    for(i=0; i < nf; i++){

       *(redn + i)= var * (1.0 - a2) / (1.0 + a2 - a22 * cos(FPI * (float) i / (float) nf));

       cdfchi(&iw, &pp, &siglev, redc+i, &dof, &status, &bound);

       if(status < 0){

          printf("****WARNING****, possible problem with Chi Square calculation    \r\n"
                 "                 of confidence cuttoff, error code returned =%d. \r\n"
                 "                 See lib/chisqr.src/cdfchi.f for further details.\n\n", status);

       }

       *(redc+i) = *(redn + i) * *(redc+i) / dof;

    }

    return;

}



#else

void cdfchi_(int * , double * , double * , double * , double * , int * , double * );

void rednoise(double var, double alph, double *redn, double *redc, double siglev, int nf)
{
    int i;
    int iw=2;
    int status;

    double a2, a22;
    double pp=1.0-siglev;
    double dof=2.0;
    double bound;

    a2 = alph * alph;
    a22 = 2.0 * alph;

    var /= nf;

    for(i=0; i < nf; i++){

       *(redn + i)= var * (1.0 - a2) / (1.0 + a2 - a22 * cos(FPI * (float) i / (float) nf));

       cdfchi_(&iw, &pp, &siglev, redc+i, &dof, &status, &bound);


       if(status < 0){

          printf("****WARNING****, possible problem with Chi Square calculation    \r\n"
                 "                 of confidence cuttoff, error code returned =%d. \r\n"
                 "                 See lib/chisqr.src/cdfchi.f for further details.\n\n", status);

       }

       *(redc+i) = *(redn + i) * *(redc+i) / dof;


    }

    return;

}



#endif

int missing(float, float, int );

void corcoff1(float *timser, int n, float *lag, float fmiss, int imiss, int icmp)

{

    int i, n1;
    int nn=0;

    float difx, dify;

    double xy1, x1, y1;

    xy1 = x1 = y1 = 0.0;

    n1 = n -1 ;

    if(imiss){

        nn = 0;

        for(i=0; i < n1; i++){


           difx = *(timser + i);
           dify = *(timser + i + 1);

           if(!missing(difx, fmiss, icmp) && !missing(dify, fmiss, icmp))  {

                 xy1 += difx * dify;;
                 x1 += difx * difx;
                 y1 += dify * dify;
                 ++nn;

           }

        }

  


    }

    else {

        for(i=0; i < n1; i++){

            difx = *(timser + i);
            dify = *(timser + i + 1);
            x1 += difx * difx;
            y1 += dify * dify;
            xy1 += difx * dify;

        }

        nn = n;

    }

    if(y1 <= TOLCOR || x1 <= TOLCOR) *lag = 0.0;
    else *lag = xy1 / sqrt(x1 * y1);



    return;

}

/* confidence levels for correlation coefficient for the null hypothesis
   cor = 0                                                               */


float conf90[27]={0.9877, 0.900, 0.805, 0.729, 0.669, 0.621, 0.582, 0.549,
                  0.521,  0.497, 0.476, 0.457, 0.441, 0.426, 0.412, 0.400,
                  0.389,  0.378, 0.369, 0.360, 0.323, 0.296, 0.257, 0.231,
                  0.211,  0.183, 0.164};

float conf95[27]={0.9969, 0.950, 0.878, 0.811, 0.754, 0.707, 0.666, 0.632,
                  0.602,  0.576, 0.553, 0.532, 0.514, 0.497, 0.482, 0.468, 
                  0.456,  0.444, 0.433, 0.423, 0.381, 0.349, 0.304, 0.273, 
                  0.250,  0.217, 0.195};

float conf99[27]={0.9999, 0.990, 0.959, 0.917, 0.875, 0.834, 0.798, 0.765,
                  0.735,  0.708, 0.684, 0.661, 0.641, 0.623, 0.606, 0.590,
                  0.575,  0.561, 0.549, 0.537, 0.487, 0.449, 0.393, 0.354,
                  0.325,  0.283, 0.254};

float confcor(int nn, int icut)

{

   int dof = nn - 2;
   int ddof; 

   if(dof < 1){

      printf("****ERROR****, insufficient degrees of freedom for confidence\r\n"    
             "               cuttoff for correlation coefficient, possible \r\n"
             "               cause, too few points!\n\n");

   }

   if(dof <= 20) ddof = dof - 1;
   else if(dof <= 25) ddof = 20;
   else if(dof <= 30) ddof = 21;
   else if(dof <= 100) ddof = 22 + ((dof - 30) / 10);
   else ddof = 26;


   switch(icut){
     case 0:
       return conf90[ddof];
       break;
     case 1:
       return conf95[ddof];
       break;
     case 2:
       return conf99[ddof];
       break;
     default:
       printf("****ERRORR****, confidence level not know for Id=%d, no confidence cuttoff available\n\n", icut);
       return -1;
   }


}
