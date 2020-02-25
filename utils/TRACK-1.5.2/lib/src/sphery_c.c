#include <Stdio.h>
#include "sphery_dat.h"
#include "bisp.h"
#include "grid.h"


/* calling program for sphery.f to produce either the interpolating
   or least squares smoothing B-spline.                                     */  




void error_interp(int , char * , int , int );

#ifdef  NOUNDERSCORE

void sphery();

#else

void sphery_();

#endif

struct cmm vv;

extern struct savedat_sphy *ssd;
extern struct sspline *ss;
extern struct sp_dat *cokn;


extern GRID *gr;

int sphery_c(double *fp, double s, int iopt)

{

    int ier=0;

    int nru[ss->mu], nrv[ss->mv];

    double fpintu[ss->nmax], fpintv[ss->nmax];
    double spu[4][ss->mu], spv[4][ss->mv];
    double bu[5][ss->nmax], bv[5][ss->nmax], cosi[ss->nmax][2];

    vv.spu = &spu[0][0];
    vv.spv = &spv[0][0];
    vv.bu = &bu[0][0];
    vv.bv = &bv[0][0];
    vv.cosi = &cosi[0][0];
    vv.nru = nru;
    vv.nrv = nrv;

    *fp = 0.;

    ss->iopt[0] = iopt;

    if(s < 0. ) error_interp(10, __FILE__, __LINE__, 1);

/* nx, tx == u(lat.),  ny, ty == v(long.)  */


#ifdef  NOUNDERSCORE

    sphery(ss->u, &ss->mu, ss->v, &ss->mv, ss->r, &ss->muv, &ss->vb, &ss->ve, 
           &ss->mumin, &ss->r0, &ss->r1, &s, &cokn->nx, cokn->tx, &cokn->ny, 
           cokn->ty, &ss->nmax, cokn->c, &ss->ncof, &ss->nuest, &ss->nvest,
           ssd->nrdatu, ssd->nrdatv, fpintu, fpintv, &ssd->nplusu,&ssd->nplusv,
           &ssd->lastdi, &ssd->lastu0, &ssd->lastu1, &ssd->fp0, &ssd->fpold,
           &ssd->reducu, &ssd->reducv, ssd->step, ssd->dr, fp, ss->iopt, 
           ss->idr, &ier);

#else

    sphery_(ss->u, &ss->mu, ss->v, &ss->mv, ss->r, &ss->muv, &ss->vb, &ss->ve, 
            &ss->mumin, &ss->r0, &ss->r1, &s, &cokn->nx, cokn->tx, &cokn->ny, 
            cokn->ty, &ss->nmax, cokn->c, &ss->ncof, &ss->nuest, &ss->nvest,
            ssd->nrdatu, ssd->nrdatv, fpintu, fpintv, &ssd->nplusu,&ssd->nplusv,
            &ssd->lastdi, &ssd->lastu0, &ssd->lastu1, &ssd->fp0, &ssd->fpold,
            &ssd->reducu, &ssd->reducv, ssd->step, ssd->dr, fp, ss->iopt, 
            ss->idr, &ier);

#endif
        

    if(ier > 0) error_interp(ier, __FILE__, __LINE__, 0);



    return ier;
}
