#include <Stdio.h>
#include <Math.h>
#include "p_vecs.h"
#include "statistic.h"

#define DENSS  0.0

void householder(double **, int , double * , double * );
void ql(double * , double * , int , double ** );
void eigsort(double *, double **, int );


/* function to estimate principle component eigen-vectors and eigen-values
   of a scatter matrix.                                                   */


void shape_vecs(struct dpt *pta, int npt, struct cvecs *pp, int direc)

{

    int j,k,l;

    double mm[3][3], yy[3], sc[3][3], *a[3];
    double d[3], e[3];
    double *cp=NULL;
    double sum=0.;

    struct dpt *pp1=NULL;

    for(j=0; j<3; j++) a[j] = &sc[j][0];

    for(j=0; j < 3; j++) {

       for(k=0; k < 3; k++)sc[j][k] = 0.0;
       d[j] = e[j] = 0.0;

    }

    cp = pp->p1;
    *(cp+3) = 0.0; 

/* Compute projection operator */

    if(direc){

       for(j=0; j < 3; j++){

          for(k=0; k < 3; k++){

             mm[j][k] = (j == k) ? 1.0 : 0.0;
 
          }

       }

    }

    else {

       for(j=0; j < 3; j++){

          for(k=0; k < 3; k++){

             mm[j][k] = ((j == k) ? 1.0 : 0.0 ) - *(cp + j) * *(cp + k);
 
          }

       }

    }


/* Estimate Scatter matrix */

   for(j=0; j < npt; j++){

       pp1 = pta + j;

       if(pp1->ic) continue;

       for(k=0; k < 3; k++){

            yy[k] = 0.0;

            yy[k] += mm[k][0] * pp1->xdt;
            yy[k] += mm[k][1] * pp1->ydt;
            yy[k] += mm[k][2] * pp1->zdt;

        }


/* Compute contribution to scatter matrix */

        for(k=0; k < 3; k++){

            for(l=0; l < 3; l++){

                 sc[k][l] += yy[k] * yy[l];

            }


        }


    }


/* Compute eigenvalues and eigenvectors */


    householder(a, 3, d, e);
    ql(d, e, 3, a);
    eigsort(d, a, 3);


/* sort eigenvalues and eigenvectors for a pair of zero degenerate eigenvalues */

    if(!direc){

       if(fabs(d[1]) < EIGTOL && fabs(d[2]) < EIGTOL){

         d[1] = d[2] = 0.;

         cp = pp->p1;

         sum = 0.;

         for(j=0; j<3; j++) sum += *(cp + j) * sc[j][1]; 


         if(fabs(sqrt(sum) - 1.0) < EIGTOL) {

            for(j=0; j<3; j++){

                sum = sc[j][1];
                sc[j][1] = sc[j][2];
                sc[j][2] = sum;

            }

         }

      }

    }

        
    cp = pp->p3;

    for(j=0; j<3; j++) *(cp+j) = sc[j][0];
    *(cp+3) = d[0];


    cp = pp->p2;
    for(j=0; j<3; j++) *(cp+j) = sc[j][1];
    *(cp+3) = d[1]; 


    cp = pp->p1;
    for(j=0; j<3; j++) *(cp+j) = sc[j][2];
    *(cp+3) = d[2];

    return;

}
