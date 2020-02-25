#include <Stdio.h>
#include "sphery_dat.h"


/* C calling function for fortran subroutine sqgrid */

#ifdef NOUNDERSCORE

void sqgrid();


extern struct cmm vv;
extern struct sspline *ss;

void sqgrid_c(double *r, int *nu, int *nv, double *c, double *fp, double *fpintu, double *fpintv)

{

    sqgrid(&ss->mu,&ss->mv,r,&ss->muv,nu,nv,c,&ss->ncof,fp,fpintu,fpintv,
            &ss->nmax,vv.spu, vv.spv, vv.bu, vv.bv, vv.nru, vv.nrv);

    return;

}

#else

void sqgrid_();


extern struct cmm vv;
extern struct sspline *ss;

#ifdef  G77
void sqgrid_c__(double *r, int *nu, int *nv, double *c, double *fp, double *fpintu, double *fpintv)
#else
void sqgrid_c_(double *r, int *nu, int *nv, double *c, double *fp, double *fpintu, double *fpintv)
#endif
{

    sqgrid_(&ss->mu,&ss->mv,r,&ss->muv,nu,nv,c,&ss->ncof,fp,fpintu,fpintv,
            &ss->nmax,vv.spu, vv.spv, vv.bu, vv.bv, vv.nru, vv.nrv);

    return;

}

#endif
