#define  DUFF_PT  -1.0e+12
#define  CHECK_PT -1.0e+10

/* union definition of feature points */

union u_pt{
       int ixy;
       float xy;
};


/* structure definition for feature points */

struct feature_pts {
       union u_pt x;
       union u_pt y;
       union u_pt ox;
       union u_pt oy;
       double gwk[3];     /* work space for geodesic measures */
       float str;
       float ostr;
       float area;
       float r_sh;
       float ornt[2];     /* [0] long. component, [1] latitude cpmponent */
       float *add_fld;
       int track_id;
       int gwky;          /* switch for assigning values to work space */
       int iob;           /* iob and ifet are id's for matching re-read */
       int ifet;          /* object and feature point data.             */
};

/* structure definition for pointer to feature points and total number of 
   feature points.                                                         */

struct features {
       struct feature_pts *fpt;
       int feature_num;
};
