#include <Stdio.h>
#include "statistic.h"
#include "grid.h"
#include "files_out.h"

#ifndef NETCDF

void netcdf_write_stats(struct tot_stat *stats, GRID *gstat, char *statsf, char *trfil, int off)
{

   printf("****ERROR****, writing statistics in netcdf format not possible\r\n"
          "               without compiling and linking with netcdf.      \n\n");
   return;

}

#else

#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include "mem_er.h"


/* function to write statistics in netcdf format */

void handle_error(int , char * , int );

void netcdf_write_stats(struct tot_stat *stats, GRID *gstat, char *statsf, char *trfil, int snum, int off)
{

    int i;
    int dim;
    int ncid;
    int ierr=NC_NOERR;
    int xid, yid, xvid, yvid;
    int dimid[2];
    int varid[NFLD];

    float *field=NULL;
    float statmiss=STATMISS;

    struct pt_stat *pst=NULL, *pss=NULL;

    char info[MAXCHR];

    static char title1[]="Tracking Statistics: densities (pdf's) and mean attributes";
    static char title2[]="Tracking Statistics: number densities and mean attributes";
    static char *slname[NFLD]={"Mean Intensity", "Std of Intensity", "Mean Speed",
                               "Std of Speed", "Feature Density", "Genesis Density",
                               "Lysis Density", "Track Density", "X-component of Mean Velocity",
                               "Y-component of Mean Velocity", "Mean Lifetime", 
                               "Mean Growth/Decay Rate", "Mean Anisotropy", 
                               "X-component of Mean Orientation Vector",
                               "Y-component of Mean Orientation Vector",
                               "Mean Tendency", "Mean Area", "Spare1", "Spare2"};
    static char *ssname[NFLD]={"mstr", "stdstr", "msp", "stdsp", "fden", "gden", "lden",
                               "tden", "xvel", "yvel", "mlif", "mgdr", "miso", "xor", "yor",
                               "mten", "marea", "spare1", "spare2"};  



    dim = gstat->ix * gstat->iy;

/* assign statistics array */

    field = (float *)calloc(dim, sizeof(float));
    mem_er((field == NULL) ? 0 : 1, dim * sizeof(float));

    ierr = nc_create(statsf, NC_CLOBBER, &ncid);

    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    if(stats->scden)
      ierr = nc_put_att_text(ncid, NC_GLOBAL, "title", strlen(title2), title2);
    else 
      ierr = nc_put_att_text(ncid, NC_GLOBAL, "title", strlen(title1), title1);

    sprintf(info, "Track File: %s; Projection group: %s; Projection type: %s; Projection centre: longitude %f, latitude %f; Density Scaling: %f", trfil, gstat->prgr_nm, gstat->prty_nm, gstat->alng, gstat->alat, stats->add_den_sc);

    ierr = nc_put_att_text(ncid, NC_GLOBAL, "history", strlen(info), info);

/* define coordinates */


    if(!gstat->prty){
/* Grid is Latitude, Longitude */
       ierr = nc_def_dim(ncid, "long", gstat->ix, &xid);
       ierr = nc_def_var (ncid, "long", NC_FLOAT, 1, &xid, &xvid);
       ierr = nc_put_att_text(ncid, xvid, "long_name", strlen("longitude"), "longitude");
       ierr = nc_def_dim(ncid, "lat", gstat->iy, &yid);
       ierr = nc_def_var (ncid, "lat", NC_FLOAT, 1, &yid, &yvid);
       ierr = nc_put_att_text(ncid, yvid, "long_name", strlen("latitude"), "latitude");
    }
    else {
/* Grid is X-Y on some projection */
       ierr = nc_def_dim(ncid, "X-pos", gstat->ix, &xid);
       ierr = nc_def_var (ncid, "X-pos", NC_FLOAT, 1, &xid, &xvid);
       ierr = nc_put_att_text(ncid, xvid, "long_name", strlen("x-pos"), "x-pos");
       ierr = nc_def_dim(ncid, "Y-pos", gstat->iy, &yid);
       ierr = nc_def_var (ncid, "Y-pos", NC_FLOAT, 1, &yid, &yvid);
       ierr = nc_put_att_text(ncid, yvid, "long_name", strlen("y-pos"), "y-pos");
    }

    dimid[0] = yid;
    dimid[1] = xid;

/* define variables */

    for(i=0; i < NFLD; i++){
      ierr = nc_def_var (ncid, ssname[i], NC_FLOAT, 2, dimid, &varid[i]);
      ierr = nc_put_att_text(ncid, varid[i], "long_name", strlen(slname[i]), slname[i]);
      ierr = nc_put_att_float(ncid, varid[i], "missing_value", NC_FLOAT, 1, &statmiss);
    }

    ierr = nc_enddef(ncid);

/* Write data for variables */

    ierr = nc_put_var_float(ncid, xvid, gstat->xgrid);
    ierr = nc_put_var_float(ncid, yvid, gstat->ygrid);

    for(i=0; i < dim; i++) *(field + i) = STATMISS;

    pst = stats->ptst + off;

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = (pss->stat1).mean; }
    ierr = nc_put_var_float(ncid, varid[0], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = (pss->stat1).var; }
    ierr = nc_put_var_float(ncid, varid[1], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = (pss->stat2).mean; }
    ierr = nc_put_var_float(ncid, varid[2], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = (pss->stat2).var; }
    ierr = nc_put_var_float(ncid, varid[3], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = pss->stat3; }
    ierr = nc_put_var_float(ncid, varid[4], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = pss->stat4; }
    ierr = nc_put_var_float(ncid, varid[5], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = pss->stat5; }
    ierr = nc_put_var_float(ncid, varid[6], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = pss->stat6; }
    ierr = nc_put_var_float(ncid, varid[7], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = (pss->stat7).xcomp; }
    ierr = nc_put_var_float(ncid, varid[8], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = (pss->stat7).ycomp; }
    ierr = nc_put_var_float(ncid, varid[9], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = pss->stat8; }
    ierr = nc_put_var_float(ncid, varid[10], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = pss->stat9; }
    ierr = nc_put_var_float(ncid, varid[11], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = pss->stat10; }
    ierr = nc_put_var_float(ncid, varid[12], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = (pss->stat11).xcomp; }
    ierr = nc_put_var_float(ncid, varid[13], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = (pss->stat11).ycomp; }
    ierr = nc_put_var_float(ncid, varid[14], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = pss->stat12; }
    ierr = nc_put_var_float(ncid, varid[15], field);

    for(i=0; i < snum; i++) {pss = pst + i; *(field + pss->apos) = pss->stat13; }
    ierr = nc_put_var_float(ncid, varid[16], field);

    ierr = nc_close(ncid);

    free(field);

    return;

}

#endif
