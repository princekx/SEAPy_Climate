#include <Stdio.h>

#ifndef NETCDF

void write_intflux_netcdf(char *ffl, double **linteg, float *slat, float *levdat, int nlat, int nlev, int imxmn, int iext, int ityp, int ctime)
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
#include "splice.h"

/* function to write composite fields to netcdf */

void handle_error(int , char * , int );

void write_intflux_netcdf(char *ffl, double **linteg, float *slat, float *levdat, int nlat, int nlev, int imxmn, int iext, int ityp, long int ctime)
{
    int i;
    int ncid, ierr=NC_NOERR;
    int latid, levid;
    int latvid, levvid;
    int dimid[2];
    int varid;

    float missval=ADD_UNDEF;

    char title[]="Azimuthaly integrated fluxes for chosen variable";
    char info[MAXCHR];

    size_t start[2]={0, 0}, count[2]={0, 0};

    count[0] = 1;
    count[1] = nlat;

    ierr = nc_create(ffl, NC_CLOBBER, &ncid);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

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

    if(ityp == 0)
       strcat(info, "combined fields integrated azimuthally");
    else if(ityp == 1)
       strcat(info, " tangential component integrated azimuthally");
    else if(ityp == 2)
       strcat(info, " radial component integrated azimuthally");
    else if(ityp == 3)
       strcat(info, " advection integrated azimuthally");
    else {
       printf("****ERROR****, component Id. not known.\n\n");
       exit(1);
    }

    ierr = nc_put_att_text(ncid, NC_GLOBAL, "history", strlen(info), info);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);


    ierr = nc_def_dim(ncid, "lat", nlat, &latid);
    ierr = nc_def_var (ncid, "lat", NC_FLOAT, 1, &latid, &latvid);
    ierr = nc_put_att_text(ncid, latvid, "long_name", strlen("latitude"), "latitude");
    ierr = nc_put_att_text(ncid, latvid, "units", 13, "degrees_north");
    ierr = nc_def_dim(ncid, "level", nlev, &levid);
    ierr = nc_def_var (ncid, "level", NC_FLOAT, 1, &levid, &levvid);
    ierr = nc_put_att_text(ncid, levvid, "long_name", strlen("Pressure"), "Pressure");
    ierr = nc_put_att_text(ncid, levvid, "units", 8, "millibar");

    dimid[0] = levid;
    dimid[1] = latid;

    ierr = nc_def_var (ncid, "fluxint", NC_DOUBLE, 2, dimid, &varid);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);
    if(ityp == 0){
       ierr = nc_put_att_text(ncid, varid, "long_name", strlen("Integrated Combined Fields"), "Integrated Combined Fields");
       if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);
    }
    else if(ityp == 1) {
       ierr = nc_put_att_text(ncid, varid, "long_name", strlen("Integrated Flux"), "Integrated Tangential Flux");
       if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);
    }
    else if(ityp == 2){
       ierr = nc_put_att_text(ncid, varid, "long_name", strlen("Integrated Flux"), "Integrated Radial Flux");
       if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);
    }
    ierr = nc_put_att_float(ncid, varid, "missing_value", NC_FLOAT, 1, &missval);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    ierr = nc_enddef(ncid);

    ierr = nc_put_var_float(ncid, latvid, slat);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);
    ierr = nc_put_var_float(ncid, levvid, levdat);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    start[0] = 0;
    for(i=0; i < nlev; i++){
       start[0] = i;
       ierr = nc_put_vara_double(ncid, varid, start, count, linteg[i]);
       if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);
    }

    ierr = nc_close(ncid);

    return;
}

void write_profile_netcdf(char *ffl, double *rintg, float *levdat, int nlev, int iext, int imxmn, long int ctime)
{
    int ncid, ierr=NC_NOERR;

    int levid;
    int levvid;
    int dimid[1];
    int varid;

    float missval=ADD_UNDEF;

    char title[]="Vertical profile of fluxes or combined fields, integrated azimuthaly and radialy, for chosen variable";
    char info[MAXCHR];

    size_t start[1]={0}, count[1]={0};

    count[0] = nlev;

    ierr = nc_create(ffl, NC_CLOBBER, &ncid);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

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

    else {
       printf("****ERROR****, stage Id. %d not valid.\n\n", iext);
       exit(1);
    }

    ierr = nc_def_dim(ncid, "level", nlev, &levid);
    ierr = nc_def_var (ncid, "level", NC_FLOAT, 1, &levid, &levvid);
    ierr = nc_put_att_text(ncid, levvid, "long_name", strlen("Pressure"), "Pressure");
    ierr = nc_put_att_text(ncid, levvid, "units", 8, "millibar");

    dimid[0] = levid;

    ierr = nc_def_var (ncid, "flx_prof", NC_DOUBLE, 1, dimid, &varid);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    ierr = nc_put_att_text(ncid, varid, "long_name", strlen("Flux Vertical Profile"), "Flux Vertical Profile");
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    ierr = nc_put_att_float(ncid, varid, "missing_value", NC_FLOAT, 1, &missval);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    ierr = nc_enddef(ncid);

    ierr = nc_put_var_float(ncid, levvid, levdat);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    ierr = nc_put_vara_double(ncid, varid, start, count, rintg);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    ierr = nc_close(ncid);

    return;

}

#endif
