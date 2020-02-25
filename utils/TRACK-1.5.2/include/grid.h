#include "complex.h"

#define  MVAL     1.0e+20      /* generic missing value  */

#define REPORT    1
#define TOLGRID   0.0005
#define GERROR    9999.9999
#define MAXFRM    10000
#define MXPRCH    50
#define TOLPOLE   1.0e-4

/* structure for reading in grid data from the data file header */

typedef struct rgrid {
    float *xgrid;
    float *ygrid;
    float *coslat;    /* cosine of latitude for area integration */
    char prgr_nm[MXPRCH]; /* grid projection group name          */
    char prty_nm[MXPRCH]; /* grid projection type name           */
    int ix;
    int iy;
    int igtyp;        /* initial grid type, global/periodic,  '0' */
                      /*       global/non-periodic            '1' */
                      /*       non-global                     '2' */
                      /*       non-uniform in X (longitude)   '3' */
                      /*       not a latitude/longitude grid  '4' */ 
    int ixfrm;        /* uniform '0' in X or not '1'              */
    int iyfrm;        /* uniform '0' in Y or not '1'              */
    int igc;          /* flag for correction to periodic grid     */
    int h_inv;        /* flag data ordered differently            */
    int gcen;         /* flag for grid longitude center,          */
                      /*         [   0, 360]                 '0'  */
                      /*         [-180, 180]                 '1'  */
                      /*       some other grid               '2'  */
    int nneg;         /* offset for correcting grid.              */
    int ox1u;         /* old grid region extents                  */
    int ox2u;
    int oy1u;
    int oy2u;
    int prgr;         /* projection group type                    */
    int prty;         /* individual projection type               */
    int iaz;          /* data on a non-Plate Caree grid at input  */
    float f_off;      /* simple translation offset                */
    float *wrk;       /* pointer to work array                    */
    float alat;       /* origin of projection if any              */
    float alng;
    float sp1;        /* standard parallels of projection, if any */
    float sp2;
    int nleng;        /* number of legendre functions for truncation */
    double *aleng;    /* pointer to array of values of associated legendre 
                         functions at latitude nodes.                      */
    double *dleng;    /* pointer to array of values of derivatives of 
                         associated legendre functions at latitude nodes.  */
    complex *sico;    /* array of complex numbers for fourier components 
                         of longitude nodes.                               */
}GRID;

/* structure for country map data */

typedef struct country {
    float *cmxg;
    float *cmyg;
    int *cmi;
    int dcm;
}CNTRY;
