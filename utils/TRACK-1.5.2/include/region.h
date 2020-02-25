/* structure definition for region averaged spectra */


typedef struct region {
   int n_grp;             /* number of grid points in region            */
   int ntmp;              /* grid point number check counter            */
   int *igrp;             /* array position of in region grid points    */
   int *nav;              /* number to average by for missing data      */
   int av_typ;            /* region averaging type                      */
   int nband;             /* number of bands for filtered power spectra */
   int reg_id;            /* region ID.                                 */
   int nfreq2;            /* number of frequencies                      */
   float x1;              /* X coordinate for region                    */
   float y1;              /* Y coordinate for region                    */
   float x2;              /* X coordinate for region                    */
   float y2;              /* Y coordinate for region                    */
   float rad;             /* radius for region                          */
   double av_var;         /* varience for SOM                           */
   double av_mn;          /* mean for SOM                               */
   double *band_av_spec;  /* band pass power for MOS                    */
   double *band_spec_av;  /* band pass power for SOM                    */
   double *av_spec;       /* MOS                                        */
   double *mos_redn;      /* MOS-red noise                              */
   double *mos_redc;      /* MOS-red noise confidence                   */
   double *som_redn;      /* SOM-red noise                              */
   double *som_redc;      /* SOM-red noise confidence                   */
   double *spec_av;       /* SOM, also used as tempory storage          */
   float *reg_tim;       /* region averaged time series                 */
} REGION;
