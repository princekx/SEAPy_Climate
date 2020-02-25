#include <Stdio.h>
#include "sphery_dat.h"


/* C function to call grisph.f with variable  array sizes */

#ifdef NOUNDERSCORE

void grisph();

extern struct sspline *ss;
extern struct cmm vv;

void grisph_c(double *r, double *dr, int *iop0, int *iop1, 
               double *tu, int *nu, double *tv, int *nv, double *p,
               int *iflag, int *ib, double *c, double *sq)


{

     int mnumv, mvv, mvnu;

     double d[*nu], dd[*nv];
     double right[(mnumv = (*nu > (mvv = ss->mv + *nv)) ? *nu : mvv)];
     double q[(mvnu=(*nu) * mvv)], au[*nu][4], av1[*nv][6], av2[*nv][4], c0[*nv];
     double a0[ss->mv][2], b0[ss->mv][2], c1[*nv], a1[ss->mv][2], b1[ss->mv][2];


/* call fortran program */

     grisph(ss->u, &ss->mu, ss->v, &ss->mv, r, &ss->muv, dr, iop0, iop1, 
             tu, nu, tv, nv, p, iflag, ib, c, &ss->ncof, sq, right, &mnumv,
             q, &mvnu, au, av1, av2, a0, a1, b0, b1, c0, c1, d, dd,
             &ss->nmax, vv.spu, vv.spv, vv.bu, vv.bv, vv.cosi, vv.nru, vv.nrv);



     return;


}

#else

void grisph_();

extern struct sspline *ss;
extern struct cmm vv;

#ifdef   G77
void grisph_c__(double *r, double *dr, int *iop0, int *iop1, 
               double *tu, int *nu, double *tv, int *nv, double *p,
               int *iflag, int *ib, double *c, double *sq)
#else
void grisph_c_(double *r, double *dr, int *iop0, int *iop1, 
               double *tu, int *nu, double *tv, int *nv, double *p,
               int *iflag, int *ib, double *c, double *sq)
#endif
{

     int mnumv, mvv, mvnu;

     double d[*nu], dd[*nv];
     double right[(mnumv = (*nu > (mvv = ss->mv + *nv)) ? *nu : mvv)];
     double q[(mvnu=(*nu) * mvv)], au[*nu][4], av1[*nv][6], av2[*nv][4], c0[*nv];
     double a0[ss->mv][2], b0[ss->mv][2], c1[*nv], a1[ss->mv][2], b1[ss->mv][2];


/* call fortran program */

     grisph_(ss->u, &ss->mu, ss->v, &ss->mv, r, &ss->muv, dr, iop0, iop1, 
             tu, nu, tv, nv, p, iflag, ib, c, &ss->ncof, sq, right, &mnumv,
             q, &mvnu, au, av1, av2, a0, a1, b0, b1, c0, c1, d, dd,
             &ss->nmax, vv.spu, vv.spv, vv.bu, vv.bv, vv.cosi, vv.nru, vv.nrv);



     return;


}


#endif
