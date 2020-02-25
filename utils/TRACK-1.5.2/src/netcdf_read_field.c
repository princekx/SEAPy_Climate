#include <Stdio.h>
#include "grid.h"
#include "netcdf_info.h"


/* function to read a netcdf field */

extern void *databuf;
extern float *abuf;
extern int std_x, std_y;

int netcdf_read_field(void *fdatin){

#ifndef NETCDF

    return 1;
}

#else
   int i;
   int ierr=NC_NOERR;
   int dim=std_x*std_y;

   NETCDF_INFO *fdat=NULL;

   fdat = (NETCDF_INFO *) fdatin;

   if(fdat->timind >= 0) {
     fdat->len1[fdat->timind] = fdat->iframe;
     printf("Frame %d, time %lf of Variable %d, %s\n", fdat->iframe + 1, fdat->timval[fdat->iframe], fdat->ifield, fdat->lname);
   }
   else {
     fdat->len1[0] = fdat->iframe;
     printf("Frame %d, Variable %d, %s\n", fdat->iframe + 1, fdat->ifield, fdat->lname);
   }


   ierr = nc_read_data(fdat->ncid, fdat->ifield, dim, fdat->dattyp, abuf, databuf, fdat->len1, fdat->len2);

   if(fdat->ioff || fdat->iscl){
      if(fdat->ioff && fdat->iscl){
         for(i=0; i < dim; i++) *(abuf + i) = *(abuf + i) * fdat->scale_fac + fdat->fld_offset;
      }
      else if(fdat->ioff) {
         for(i=0; i < dim; i++) *(abuf + i) += fdat->fld_offset;
      }
      else {
         for(i=0; i < dim; i++) *(abuf + i) *= fdat->scale_fac;
      }
   }

   if(ierr != NC_NOERR) {
      printf("%s\n", nc_strerror(ierr)); 
      printf("****WARNING****, possible End of File encountered.\n\n");
      return 1;

   }

   return 0;

}

#endif
