/* header file for data saved from previous call of reg_surf */

#define  NITER   3
#define  SM_FAC  2.

struct savedat{           /* data found from a previous call */
       int *nrdatx;
       int *nrdaty;
       int nplusx;
       int nplusy;
       int lastdi;
       double fp0;
       double fpold;
       double reducx;
       double reducy;
};

struct rspline {           /* knot data */
       int nxest;
       int nyest;
       int nmax;
       int ncof;

       int kx;           /* spline degree data */
       int ky;
       int kx1;
       int ky1;
       int kx2;
       int ky2;
       int kmax;
       int kmax2;
       int k2max2;

       double *x;         /* region data */
       double *y; 
       double *z;       
       double xb;
       double xe;
       double yb;
       double ye;
       int sx;
       int sy;
       int ex;
       int ey;
       int mx;
       int my;
       int mxy;

};


#define     KX     3      /* degree of the spline approx. in the x-direction */
#define     KY     3      /* degree of the spline approx. in the y-direction */

