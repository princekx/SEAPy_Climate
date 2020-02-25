#include <Stdio.h>
#include <stdlib.h>
#include <reg_dat.h>
#include "bisp.h"
#include "grid.h"
#include "mem_er.h"


/* calling program for smoopy.f to produce either the interpolating
   or least squares smoothing B-spline.                                       */

/* This function and the dependent functions need re-writing to reduce the size
   of the parameter lists using correspondance between structures and common 
   blocks.                                                                       */ 

void error_interp(int , char * , int , int );

#ifdef  NOUNDERSCORE

void smoopy(double *, int *, double *, int *, double *, int *, double *, double *,
            double * , double *, int * , int * , double *, int * , double *, int *,
            double *, int *, double *, int * ,int *, int *, int *, int *, double *,
            double *, int *, int *, int *, double *, double *, double *, double *,
            double *, int *, int *);

#else

void smoopy_(double *, int *, double *, int *, double *, int *, double *, double *,
             double * , double *, int * , int * , double *, int * , double *, int *,
             double *, int *, double *, int * ,int *, int *, int *, int *, double *,
             double *, int *, int *, int *, double *, double *, double *, double *,
             double *, int *, int *);

#endif

extern struct savedat *sd;
extern struct rspline *rs;
extern struct sp_dat *cokn;


extern GRID *gr;

int *nrx=NULL, *nry=NULL;
double *spx=NULL, *spy=NULL;

int smoopy_c(double *fp, double s, int iopt)

{

    int ier;

    double fpintx[rs->nmax], fpinty[rs->nmax];

    *fp = 0.;

    if(s < 0. ) error_interp(10, __FILE__, __LINE__, 1);

    nrx = (int *)calloc(rs->mx, sizeof(int));
    mem_er((nrx == NULL) ? 0 : 1, rs->mx * sizeof(int));

    nry = (int *)calloc(rs->my, sizeof(int));
    mem_er((nry == NULL) ? 0 : 1, rs->my * sizeof(int));

    spx = (double *)calloc(rs->kx1*rs->mx, sizeof(double));
    mem_er((spx == NULL) ? 0 : 1, rs->kx1*rs->mx * sizeof(double));

    spy = (double *)calloc(rs->ky1*rs->my, sizeof(double));
    mem_er((spy == NULL) ? 0 : 1, rs->ky1*rs->my * sizeof(double));


#ifdef  NOUNDERSCORE

    smoopy(rs->x, &rs->mx, rs->y, &rs->my, rs->z, &rs->mxy, &rs->xb, &rs->xe, 
           &rs->yb, &rs->ye, &rs->kx, &rs->ky, &s, &cokn->nx, cokn->tx, &cokn->ny,
           cokn->ty, &rs->nmax, cokn->c, &rs->ncof, &rs->nxest, &rs->nyest, 
           sd->nrdatx, sd->nrdaty, fpintx, fpinty, &sd->nplusx, &sd->nplusy,
           &sd->lastdi, &sd->fp0, &sd->fpold, &sd->reducx, &sd->reducy, fp, 
           &iopt, &ier);

#else

    smoopy_(rs->x, &rs->mx, rs->y, &rs->my, rs->z, &rs->mxy, &rs->xb, &rs->xe, 
            &rs->yb, &rs->ye, &rs->kx, &rs->ky, &s, &cokn->nx, cokn->tx, &cokn->ny,
            cokn->ty, &rs->nmax, cokn->c, &rs->ncof, &rs->nxest, &rs->nyest, 
            sd->nrdatx, sd->nrdaty, fpintx, fpinty, &sd->nplusx, &sd->nplusy,
            &sd->lastdi, &sd->fp0, &sd->fpold, &sd->reducx, &sd->reducy, fp, 
            &iopt, &ier);

#endif


    free(nrx);
    free(nry);
    free(spx);
    free(spy);

    if(ier > 0) error_interp(ier, __FILE__, __LINE__, 0);

    return ier;
}
