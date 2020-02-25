#include <Stdio.h>
#include "grid.h"
#include "netcdf_info.h"

/* function to add a field to a netcdf file, written in a curse form */


void write_netcdf_field(float *arr, void *fdat)
{

#ifndef NETCDF

    printf("****WARNING****, netcdf support not available, field not written.\n\n");
    return;

}

#else

    size_t tind[1];
    size_t len1[NC_MY_MAX_DIM]={0, 0, 0, 0, 0};
    int ierr=0;

    NETCDF_INFO *fd=NULL;

    fd = (NETCDF_INFO *)fdat;

    tind[0] = fd->nframe;

    if((ierr = nc_put_var1_double(fd->ncid, fd->timid, tind, fd->timval + fd->iframe)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    fd->timval[fd->nframe] = fd->timval[fd->iframe];

    len1[fd->timind] = fd->nframe;

    if((ierr = nc_put_vara_float(fd->ncid, fd->ifield, len1, fd->len2, arr)) != NC_NOERR)
       handle_error(ierr,  __FILE__, __LINE__); 

    return;

}

#endif
