/* header file for data saved from previous call of sph_surf */

struct savedat_sphy{           /* data found from a previous call */
     int *nrdatu;
     int *nrdatv;
     int nplusu;
     int nplusv;
     int lastdi;
     int lastu0;
     int lastu1;
     double fp0;
     double fpold;
     double reducu;
     double reducv;
     double step[2];
     double dr[6];
};

struct sspline {
     int nuest;
     int nvest;
     int nmax;
     int ncof;
     int mvv;

     double *u;
     double *v;
     double *r;

     double vb;
     double ve;
     double ue;
     double ub;

     double r0;
     double r1;

     int mu;
     int mv;
     int muv;
     int mumin;

     int iopt[3];
     int idr[4];

};


/* sphery uses only bicubic splines (hard wired)*/


struct cmm {
     double *spu;
     double *spv;
     double *bu;
     double *bv;
     double *cosi;
     int *nru;
     int *nrv;
};
