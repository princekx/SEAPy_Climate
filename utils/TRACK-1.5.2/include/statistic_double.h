/* generic data structure for spatial statistics */

struct tr_flux_d{
     double xcomp;
     double ycomp;
};

struct tr_m_var_d{
     double mean;
     double var;
};


struct  pt_stat_d {
     int ptyp;                /* projection variable, 0 or 1          */
     int apos;                /* array position of estimation point   */
     double num;               /* computational variable               */
     double xs;                /* X-sample position                    */
     double ys;                /* Y-sample position                    */
     struct tr_m_var_d stat1;   /* mean strength                        */
     struct tr_m_var_d stat2;   /* mean speed                           */
     double stat3;             /* phenomena density                    */
     double stat4;             /* genesis                              */
     double stat5;             /* lysis                                */
     double stat6;             /* track density                        */
     struct tr_flux_d stat7;    /* track flux                           */
     double stat8;             /* lifetime                             */
     double stat9;             /* growth/decay rate, (1/f)df/dt        */
     double stat10;            /* anisotropy                           */
     struct tr_flux_d stat11;   /* orientation vector                   */
     double stat12;            /* tendency, df/dt                      */ 
     double stat13;            /* mean area                            */
     double spare1;            /* space for additional or derived data */
     double spare2;            /* space for additional or derived data */
};  

struct tot_stat_d {
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
     double add_den_sc;
     double datnm[STNM];
     double xa1;
     double xa2;
     double ya1;
     double ya2;
     double sm[STNM];
     struct  pt_stat_d *ptst;
};

