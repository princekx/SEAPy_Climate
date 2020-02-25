/* generic data structure for spatial statistics */

#define   STNM    17              /* number of statistical measures */
#define   NFLD    19              /* number of statistical fields   */
#define   MXCHR   300             /* max. string array length       */

#define   TOLBOUND 0.001          /* tolerence for kernel-boundary overlap */

#define   STATMISS 1.0e+10        /* missing value for NETCDF output */ 


struct tr_flux{
     float xcomp;
     float ycomp;
};

struct tr_m_var{
     float mean;
     float var;
};


struct  pt_stat {
     int ptyp;                /* projection variable, 0 or 1          */
     int apos;                /* array position of estimation point   */
     float num;               /* computational variable               */
     float xs;                /* X-sample position                    */
     float ys;                /* Y-sample position                    */
     struct tr_m_var stat1;   /* mean strength                        */
     struct tr_m_var stat2;   /* mean speed                           */
     float stat3;             /* phenomena density                    */
     float stat4;             /* genesis                              */
     float stat5;             /* lysis                                */
     float stat6;             /* track density                        */
     struct tr_flux stat7;    /* track flux                           */
     float stat8;             /* lifetime                             */
     float stat9;             /* growth/decay rate, (1/f)df/dt        */
     float stat10;            /* anisotropy                           */
     struct tr_flux stat11;   /* orientation vector                   */
     float stat12;            /* tendency, df/dt                      */ 
     float stat13;            /* mean area                            */
     float spare1;            /* space for additional or derived data */
     float spare2;            /* space for additional or derived data */
};  

struct tot_stat {
     int ptnum;               /* number of estimation points             */
     int scden;               /* density type tag:                       */
                              /*      '0' pdf's except for track density */
                              /*      '1' traditional number densities   */
                              /*      '2' all pdf's including track      */
                              /*      '3' all number densities scaled    */
                              /*          from pdf's.                    */
     int kern[STNM];
     int ibound;
     int ibb[4];
     float add_den_sc;
     float datnm[STNM];
     float xa1;
     float xa2;
     float ya1;
     float ya2;
     float sm[STNM];
     struct  pt_stat *ptst;
};


/* structure for data input to density estimation */

struct dpt {
     int ic;
     float lng;
     float lat;
     double xdt;
     double ydt;
     double zdt;
     float sdt;
     float dsdt;
     float vec[2];
     float dvec[2];
};


#define   DBAND     20.0000000   /* smoothing parameter at which evaluation 
                                    of exponential of fisher density changes. */

#define   EXPON     40.0         /* limiting exponent for fisher density */

#define  SMMIN    1.0e-6         /* minimum value of smoothing parameter */

#define   DTOL    1.0e-4

#define  TOLWT    1.0e-10         /* tolerance for sum of additional weights */
