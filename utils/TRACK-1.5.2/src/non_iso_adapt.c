#include <Stdio.h>
#include <Math.h>
#include <stdlib.h>
#include "mem_er.h"
#include "statistic.h"
#include "p_vecs.h"
#include "m_values.h"


/* function to compute adaptive parameters for non-isotropic kernel (power). */


void estimate_vectors(double (*density)(), double , struct cvecs * , double * , float * ,double  ,int );
double const_density(struct dpt * , struct dpt * , double );
double linear_density(struct dpt * , struct dpt * , double );
double quad_density(struct dpt * , struct dpt * , double );

double non_iso_adapt(struct dpt *dt, int dtn, double smplt, double *plt, float *wght, struct cvecs *pcom, int ind)

{

    int i;
    int id=1;
    static int ks;

    double dnorm=0.0;
    double tsm, cn, cd=0.;
    double llh=0.0;
    double nn=0;

    struct cvecs *pp;
    struct dpt *dd;

    static double (*density)(struct dpt * , struct dpt * , double )=NULL;

/* put input data into first vector position */


    for(i=0; i< dtn; i++){

       dd = dt + i;
       pp = pcom + i;
       pp->p1[0] = dd->xdt;
       pp->p1[1] = dd->ydt;
       pp->p1[2] = dd->zdt;


    }

    if(!wght){
      cn = 1.0 / (2.0 * FPI * (float)dtn);
    }
    else { 
      nn = 0.0;
      for(i=0; i < dtn; i++) nn += *(wght + i);
      if(nn > TOLWT) cn = 1.0 / (2.0 * FPI * nn);
      else cn = 0.0;

    }

    if(smplt < SMMIN){tsm = -1.0; dnorm = 1.0; id = 0;}
    else {
       tsm = (1.0 + smplt) / smplt;
       cn *= tsm;
       cd = tsm - 1.0;
    }

    if(ind){

       printf("which smoothing isotropic kernal do you require for pilot estimate,\r\n"
              "        input '1', for constant (spherical cap)                    \r\n"
              "        input '2', for linear,                                     \r\n"  
              "        input '3', for quadratic.                                  \n\n");

       scanf("%d", &ks);

    }


    switch(ks){
       case 1:
         density = const_density;
         if(id) dnorm = cn / cd;
         break;
       case 2:
         density = linear_density;
         if(id) dnorm = 2.0 * cn / (cd * cd);
         break;
       case 3:
         density = quad_density;
         if(id) dnorm = 3.0 * cn / (cd * cd * (tsm + 2.0));
         break;
       default:
         printf("***error***, no kernel available for this index %s\n\n", __FILE__);
         exit(1);
    }


    estimate_vectors(density, dnorm, pcom, plt, wght, tsm, dtn);

    if(plt){

      for(i=0; i < dtn; i++) llh += log(*(plt+i));

    } 

    return llh;

}
