/* data structure for combined track data */

#define  ADD_UNDEF   1.0e+25
#define  ADD_CHECK   1.0e+20

#ifndef MAXCHR
#define MAXCHR 500
#endif

struct fet_pt_tr {
       int fr_id;             /* Frame id.                                   */
       int obj_id;            /* Object Id. for track merging.               */
       long int time;         /* Frame time.                                 */
       int nfm;               /* Number of missing frames from current frame */
       int npt;               /* Number of points for averaging tracks in EPS*/
       int ix;                /* Nearest X-grid position                     */
       int iy;                /* Nearest Y-grid position                     */
       int imp;               /* Flag for merged point                       */
       float xf;              /* Longitude position                          */
       float yf;              /* Latitude position                           */
       float zf;              /* Strength                                    */
       double pp[3];          /* Position in Cartesian space                 */
       float gwthr;           /* Growth/decay rate, (1/f)df/dt               */
       float tend;            /* Tendency, df/dt                             */
       float vec[2];          /* Vector for velocities                       */
       float sh_an;           /* Anisotropy value                            */
       float or_vec[2];       /* Orientation vector                          */
       float area;            /* Area associated with point                  */
       float wght;            /* Additional weighting for statistics         */
       float *add_fld;        /* pointer to additional field values          */
       float *atmp;           /* temporary pointer to additional fields      */
};

struct tot_tr {
       int num;                /* number of points on track.                          */
       int trid;               /* track ID.                                           */
       int tr_sp;
       int eofs;
       int awt;                /* flag for weighted tracks.                         */
       int imt;                /* flag for merged tracks, '0' for no merging.       */
                               /*                         '1' for averaged merge.   */
                               /*                         '2' for no averaged merge */
       int ntt;                /* number of tracks combined for principle track.    */
       long int time;          /* time of start of track                            */
       struct fet_pt_tr *trpt; /* pointer to track points                           */
/*       complex *fr_cof;    */    /* Fourier modes for track.                      */
};
