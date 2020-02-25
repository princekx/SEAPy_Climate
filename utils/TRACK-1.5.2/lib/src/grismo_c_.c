#include <Stdio.h>
#include <reg_dat.h>

/* C calling routine for grismo.f, allows variable array lengths to be used */

#ifdef NOUNDERSCORE

void grismo();
void sqgrid();

extern struct rspline *rs;

extern int *nrx, *nry;
extern double *spx, *spy;

void grismo_c(double *z, double *tx, int *nx,
              double *ty, int *ny, double *p, double *c, double *fp,
              double *fpintx, double *fpinty, int *iflag)

{

   int mnxy, mnxmy, mynx;

   double h[rs->kmax2], right[(mnxmy=(*nx > rs->my) ? *nx : rs->my)];
   double q[(mynx=*nx * rs->my)];
   double ax[rs->kx1][*nx], ay[rs->ky1][*ny], d[(mnxy=(*nx > *ny) ? *nx : *ny)];
/*   double spx[rs->kx1][rs->mx], spy[rs->ky1][rs->my]; */
   double bx[rs->kx2][mnxy], by[rs->ky2][mnxy];

/*   int nrx[rs->mx], nry[rs->my]; */


   grismo(rs->x, &rs->mx, rs->y, &rs->my, z, &rs->mxy, tx, nx, ty, ny, p, c, 
           &rs->ncof, iflag, h, right, &mnxmy, q, &mynx, ax, ay, d, &mnxy, spx,
           spy, bx, by, nrx, nry);

/*  calculate for each knot interval tx(j+kx)<=x<=tx(j+kx+1) (ty(j+ky)<=
    y<=ty(j+ky+1)) the sum of squared residuals fpintx(j),j=1,2,...,
    nx-2*kx-1 (fpinty(j),j=1,2,...,ny-2*ky-1) for the data points having
    their absciss (ordinate) belonging to that interval.
    fp gives the total sum of squared residuals.                            */

    sqgrid(&rs->mx,&rs->my,z,&rs->mxy,nx,ny,c,&rs->ncof,fp,fpintx,fpinty,
            &mnxmy,spx,spy,bx,by,nrx,nry);

   return;
}


#else

void grismo_();
void sqgrid_();

extern struct rspline *rs;

extern int *nrx, *nry;
extern double *spx, *spy;

#ifdef   G77
void grismo_c__(double *z, double *tx, int *nx,
              double *ty, int *ny, double *p, double *c, double *fp,
              double *fpintx, double *fpinty, int *iflag)
#else
void grismo_c_(double *z, double *tx, int *nx,
              double *ty, int *ny, double *p, double *c, double *fp,
              double *fpintx, double *fpinty, int *iflag)
#endif
{

   int mnxy, mnxmy, mynx;

   double h[rs->kmax2], right[(mnxmy=(*nx > rs->my) ? *nx : rs->my)];
   double q[(mynx=*nx * rs->my)];
   double ax[rs->kx1][*nx], ay[rs->ky1][*ny], d[(mnxy=(*nx > *ny) ? *nx : *ny)];
/*   double spx[rs->kx1][rs->mx], spy[rs->ky1][rs->my]; */
   double bx[rs->kx2][mnxy], by[rs->ky2][mnxy];

/*   int nrx[rs->mx], nry[rs->my]; */


   grismo_(rs->x, &rs->mx, rs->y, &rs->my, z, &rs->mxy, tx, nx, ty, ny, p, c, 
           &rs->ncof, iflag, h, right, &mnxmy, q, &mynx, ax, ay, d, &mnxy, spx,
           spy, bx, by, nrx, nry);

/*  calculate for each knot interval tx(j+kx)<=x<=tx(j+kx+1) (ty(j+ky)<=
    y<=ty(j+ky+1)) the sum of squared residuals fpintx(j),j=1,2,...,
    nx-2*kx-1 (fpinty(j),j=1,2,...,ny-2*ky-1) for the data points having
    their absciss (ordinate) belonging to that interval.
    fp gives the total sum of squared residuals.                            */

    sqgrid_(&rs->mx,&rs->my,z,&rs->mxy,nx,ny,c,&rs->ncof,fp,fpintx,fpinty,
            &mnxmy,spx,spy,bx,by,nrx,nry);

   return;
}

#endif
