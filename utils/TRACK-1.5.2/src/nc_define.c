#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grid.h"
#include "mem_er.h"
#include "netcdf_info.h"


extern int std_x, std_y;
extern int frnum;

NETCDF_INFO *nc_define(NETCDF_INFO *fdatin, char *nfile)
{

    NETCDF_INFO *fdat=NULL;

#ifdef NETCDF

    int ierr = NC_NOERR;
    int vard[5];
    int ndim=0;
    int timed=0;

    int adddim, timdim, levdim, latdim, lngdim;

    size_t adim[1];

    adddim = timdim = levdim = latdim = lngdim = -1;

    fdat = nc_clone(fdatin, nfile, NC_CREATE_MODE);

/* define dimensions and attributes */

    if(fdat->addid >= 0){
       if((ierr = nc_def_dim(fdat->ncid, fdat->dimn[fdat->addind], 1, &adddim)) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);
    }

    if(fdat->timid >= 0){
       timed = (fdat->addid >= 0) ? frnum : NC_UNLIMITED;
       if((ierr = nc_def_dim(fdat->ncid, fdat->dimn[fdat->timind], timed, &timdim)) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);
    }

    if(fdat->levid >= 0){
       if((ierr = nc_def_dim(fdat->ncid, fdat->dimn[fdat->levind], 1, &levdim)) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);
    }

    if((ierr = nc_def_dim(fdat->ncid, fdat->dimn[fdat->latind], std_y, &latdim)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    if((ierr = nc_def_dim(fdat->ncid, fdat->dimn[fdat->lngind], std_x, &lngdim)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);


/* define variables */

    if(fdat->addid >= 0){

       vard[0] = adddim;

       if((ierr = nc_def_var(fdat->ncid, fdat->dimn[fdat->addind], NC_FLOAT, 1, vard, &fdat->addid)) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);
    }


    if(fdat->timid >= 0){

       vard[0] = timdim;

       if((ierr = nc_def_var(fdat->ncid, fdat->dimn[fdat->timind], NC_DOUBLE, 1, vard, &fdat->timid)) != NC_NOERR)
           handle_error(ierr, __FILE__, __LINE__);
    }

    if(fdat->levid >= 0){

       vard[0] = levdim;

       if((ierr = nc_def_var(fdat->ncid, fdat->dimn[fdat->levind], NC_FLOAT, 1, vard, &fdat->levid)) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);

    }

    vard[0] = latdim;

    if((ierr = nc_def_var(fdat->ncid, fdat->dimn[fdat->latind], NC_FLOAT, 1, vard, &fdat->latid)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    vard[0] = lngdim;

    if((ierr = nc_def_var(fdat->ncid, fdat->dimn[fdat->lngind], NC_FLOAT, 1, vard, &fdat->longid)) != NC_NOERR) 
       handle_error(ierr, __FILE__, __LINE__);

    nc_put_att_text(fdat->ncid, fdat->timid, "units", NC_MAX_NAME, "");
    nc_copy_att(fdatin->ncid, fdatin->timid, "units", fdat->ncid, fdat->timid);

    if(fdat->addid >= 0) vard[fdat->addind] = adddim;
    vard[fdat->timind] = timdim;
    if(fdat->levid >= 0) vard[fdat->levind] = levdim;
    vard[fdat->latind] = latdim;
    vard[fdat->lngind] = lngdim;

    if(fdat->addid >= 0){
      ndim = (fdat->levid >= 0) ? 5 : 4;
    }
    else ndim = (fdat->levid >= 0) ? 4 : 3; 


    if((ierr = nc_def_var(fdat->ncid, fdat->name, NC_FLOAT, ndim, vard, &fdat->ifield)) != NC_NOERR) 
       handle_error(ierr, __FILE__, __LINE__);

/* assign attributes */

    if((ierr = nc_put_att_text(fdat->ncid, fdat->ifield, "long_name", NC_MAX_NAME, fdat->lname)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    if((ierr = nc_put_att_text(fdat->ncid, fdat->ifield, "units", NC_MAX_NAME, fdat->units)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    if(fdat->imiss){
       if((ierr = nc_put_att_float(fdat->ncid, fdat->ifield, "missing_value", NC_FLOAT, 1, &(fdat->missing))) != NC_NOERR)
           handle_error(ierr, __FILE__, __LINE__);
    }


/* end define mode */

    if((ierr = nc_enddef(fdat->ncid)) != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

/* write dimension data */

    if(fdat->addid >= 0){

       adim[0] = 0;
       if((ierr = nc_put_var1_float(fdat->ncid, fdat->addid, adim, &(fdat->addval))) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);

    }

    if(fdat->levid >= 0){

       adim[0] = 0;
       if((ierr = nc_put_var1_float(fdat->ncid, fdat->levid, adim, &(fdat->levval))) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);

    }

    if((ierr = nc_put_var_float(fdat->ncid, fdat->latid, fdat->latgr)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);
    

    if((ierr = nc_put_var_float(fdat->ncid, fdat->longid, fdat->lnggr)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    fdat->dattyp = NC_FLOAT;


#else

    return NULL;

#endif

    return fdat;

}
