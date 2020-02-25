#include <Stdio.h>

#ifndef NETCDF

float *write_fields_netcdf(char *ffl, double **avg, float *slng, float *slat, float *levdat, int nlng, int nlat, int nlev, int iext, int imxmn, long int ctime)
{

   printf("****ERROR****, writing statistics in netcdf format not possible\r\n"
          "               without compiling and linking with netcdf.      \n\n");
   return NULL;

}

#else

#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include "mem_er.h"
#include "splice.h"

/* function to write composite fields to netcdf */

void handle_error(int , char * , int );
float *leveldata(int );

float *write_fields_netcdf(char *ffl, double **avg, float *slng, float *slat, float *levdat, int nlng, int nlat, int nlev, int iext, int imxmn, int iwind, int nf, long int ctime)
{
    int i;
    int ncid;
    int ierr=NC_NOERR;
    int lngid, latid, levid;
    int lngvid, latvid, levvid;
    int varid;
    int dimid[3];

    size_t len1[3]={0, 0, 0};
    size_t len2[3]={1, 0, 0};

    static char title[MAXCHR];
    static char varnm[MAXCHR];
    char info[MAXCHR];
    char varnml[MAXCHR];

    float missval=ADD_UNDEF;

    len2[1] = nlat;
    len2[2] = nlng;

    ierr = nc_create(ffl, NC_CLOBBER, &ncid);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    if(nf == 0){
       gets(title);
       printf("Provide a title for this output\n\n");
       gets(title);
    }

    ierr = nc_put_att_text(ncid, NC_GLOBAL, "title", strlen(title), title);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    if(!iext){
       if(!imxmn) sprintf(info, "Composite at minimum intensity.");
       else sprintf(info, "Composite at maximum intensity.");
    }
    else if(iext == 1){
       if(!imxmn) sprintf(info, "Composite at minimum tendency.");
       else sprintf(info, "Composite at maximum tendency.");
    }
    else if(iext == 2) sprintf(info, "Composite at maximum vertical gradient.");
    else if(iext == 3) sprintf(info, "Composite at time %ld, real time or timestep.", ctime);
    else if(iext == 4) sprintf(info, "Composite at time of genesis.");

    else {
       printf("****ERROR****, stage Id. %d not valid.\n\n", iext);
       exit(1);
    }

    ierr = nc_put_att_text(ncid, NC_GLOBAL, "history", strlen(info), info);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    ierr = nc_def_dim(ncid, "long", nlng, &lngid);
    ierr = nc_def_var (ncid, "long", NC_FLOAT, 1, &lngid, &lngvid);
    ierr = nc_put_att_text(ncid, lngvid, "long_name", strlen("longitude"), "longitude");
    ierr = nc_put_att_text(ncid, lngvid, "units", 12, "degrees_east");
    ierr = nc_def_dim(ncid, "lat", nlat, &latid);
    ierr = nc_def_var (ncid, "lat", NC_FLOAT, 1, &latid, &latvid);
    ierr = nc_put_att_text(ncid, latvid, "long_name", strlen("latitude"), "latitude");
    ierr = nc_put_att_text(ncid, latvid, "units", 13, "degrees_north");
    ierr = nc_def_dim(ncid, "level", nlev, &levid);
    ierr = nc_def_var (ncid, "level", NC_FLOAT, 1, &levid, &levvid);
    ierr = nc_put_att_text(ncid, levvid, "long_name", strlen("Pressure"), "Pressure");
    ierr = nc_put_att_text(ncid, levvid, "units", 8, "millibar");

    if(!levdat) levdat = leveldata(nlev);

    dimid[0] = levid;
    dimid[1] = latid;
    dimid[2] = lngid;

    if(!iwind){
       printf("What is the short name of the composite variable?\n\n");
       scanf("%s", varnm);
       printf("What is the long name of the composite variable?\n\n");
       fflush(stdin);
       gets(varnml);
    }
    else {
       if(!nf) {sprintf(varnm, "wind"); sprintf(varnml, "Wind Speed");}
       else if(nf == 1) {sprintf(varnm, "wtan"); sprintf(varnml, "Tangential or Zonal Wind Speed");}
       else if(nf == 2) {sprintf(varnm, "wrad"); sprintf(varnml, "Radial or Meridional Wind Speed");}
       else if(nf == 3) {sprintf(varnm, "wrad"); sprintf(varnml, "Tangential Flux");}
       else if(nf == 4) {sprintf(varnm, "wrad"); sprintf(varnml, "Radial Flux");}
       else if(nf == 5) {sprintf(varnm, "wrad"); sprintf(varnml, "Advection");}
    }

    

    ierr = nc_def_var (ncid, varnm, NC_DOUBLE, 3, dimid, &varid);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);
    ierr = nc_put_att_text(ncid, varid, "long_name", strlen(varnml), varnml);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);
    ierr = nc_put_att_float(ncid, varid, "missing_value", NC_FLOAT, 1, &missval);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    ierr = nc_enddef(ncid);

/* Write data for variables */

    ierr = nc_put_var_float(ncid, lngvid, slng);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);
    ierr = nc_put_var_float(ncid, latvid, slat);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);
    ierr = nc_put_var_float(ncid, levvid, levdat);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    for(i=0; i < nlev; i++) {

        len1[0] = i;
        ierr = nc_put_vara_double(ncid, varid, len1, len2, *(avg + i));
        if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);
    }

    ierr = nc_close(ncid);

    return levdat;
}


void handle_error(int status, char *file, int line) {

   printf("%s in file %s at line %d\n", nc_strerror(status), file, line);
   exit(1);
   return;

}


#endif
