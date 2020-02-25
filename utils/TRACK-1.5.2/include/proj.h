
/* include file for projection function prototypes */


#define  TOLPROJ    1.0e-4


void proj_report(int , int );

/* prototypes for cylindrical projections */

void eq_rect(float * , int , int );
void eq_area(float * , int , int );
void urmaev(float * , int , int );
void mercator(float * , int , int );
void stereog(float * , int , int );
void miller(float * , int , int );

/* typedef for cylindrical projections */

typedef void (*PRFP)(float * , int , int );

/* prototype for azimuthal projection */

int azimuthal(float * , float * , float , float , int , int );

int conic(float * , float * , float , float , float , float , int , int , int * );

int rotated_cyln(float * , float * , float , float , int , int * );

typedef struct proj_ba{
     PRFP  prj1;
     PRFP  prj2;
}PROJ;

typedef struct grid_param{
     int prgr;
     int prty;
     float alng;
     float alat;
     float sp1;
     float sp2;
}GRDP;


PRFP proj_assign(int );
void trans_proj(PRFP , int );
PROJ map_project(int * , int );
void hemi_grid(int );
void conic_grid(int );
void rotated_cyln_grid(int );
PROJ proj_group(int ,int );

int proj_point(PROJ , float * , float * , GRDP * , GRDP * );
void vect_proj(float, float, float * , float * , PROJ , GRDP * , GRDP *);
