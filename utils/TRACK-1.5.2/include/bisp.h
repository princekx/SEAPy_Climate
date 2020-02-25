
/* include file to define structure for B-spline data */

struct sp_dat {
    int sm_type;
    int nx;
    int ny;
    int ncof;
    double *tx;
    double *ty;
    double *c;
};
