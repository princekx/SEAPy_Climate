#include <Stdio.h>
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#include "grid.h"
#include "mem_er.h"
#include "netcdf_info.h"

extern int std_x, std_y;
extern int frnum;

NETCDF_INFO *nc_clone(NETCDF_INFO *fdatin, char *nfile, int mode)
{

    NETCDF_INFO *fdat=NULL;

#ifdef NETCDF

    int i, j;
    int ierr=NC_NOERR;
    int nvars=0, nvdim=0;
    int ndims=0, natts=0;
    int dims[NC_MAX_VAR_DIMS];
    int nfield=0;
    int isvar=0;
    int levid=0, addid=0;
    int dum=0;
    size_t nlev=0, atlen=0, nadd=0;
    size_t len1, len2;

    char **dim_name=NULL;
    char **var_name=NULL;

    float *llev=NULL, *addd=NULL;
    float alev=0.0;

    void *gtmp=NULL;

    nc_type type;

    fdat = (NETCDF_INFO *)malloc_initl(sizeof(NETCDF_INFO));
    mem_er((fdat == NULL) ? 0 : 1, sizeof(NETCDF_INFO));

    if(fdatin && nfile){
       printf("****INFORMATION****, NETCDF file to open is:- \r\n"
              "                     %s\n\n", nfile);
    }

    fdat->invar = fdatin->invar;

    if(fdatin || mode == NC_SAME) {

       memcpy(fdat, fdatin, sizeof(NETCDF_INFO));

       fdat->lnggr = (float *)calloc(std_x, sizeof(float));
       mem_er((fdat->lnggr == NULL) ? 0 : 1, std_x * sizeof(float));
       memcpy(fdat->lnggr, fdatin->lnggr, std_x * sizeof(float));

       fdat->latgr = (float *)calloc(std_y, sizeof(float));
       mem_er((fdat->latgr == NULL) ? 0 : 1, std_y * sizeof(float));
       memcpy(fdat->latgr, fdatin->latgr, std_y * sizeof(float));

       fdat->timval = (double *)calloc(frnum, sizeof(double));
       mem_er((fdat->timval == NULL) ? 0 : 1, frnum * sizeof(double));
       memcpy(fdat->timval, fdatin->timval, frnum * sizeof(double));

    }

    else {
       printf("****WARNING****, no netcdf info to clone, opening file \r\n"
              "%s anyway.\n\n", nfile);

    }


    if(mode == NC_CREATE_MODE){

       if((ierr = nc_create(nfile, NC_CLOBBER|NC_SHARE, &(fdat->ncid))) != NC_NOERR) 
          handle_error(ierr, __FILE__, __LINE__);

       printf("****INFORMATION****, NETCDF file %s created.\n\n", nfile);

       fdat->iscl = 0;
       fdat->ioff = 0;
       fdat->dattyp = NC_FLOAT;

    }

    else if(mode == NC_OPEN_MODE || mode == NC_SAME) {

       if(mode == NC_OPEN_MODE){
          if((ierr = nc_open(nfile, NC_NOWRITE, &(fdat->ncid))) != NC_NOERR) 
             handle_error(ierr, __FILE__, __LINE__);

          printf("****INFORMATION****, NETCDF file %s opened for 'r'.\n\n", nfile);

       }

       else if(mode == NC_SAME){
         fdat->ncid = fdatin->ncid;
       }

       nc_inq_nvars(fdat->ncid, &nvars);
       nc_inq_ndims(fdat->ncid, &ndims);

       dim_name = (char **)calloc(ndims, sizeof(char *));
       mem_er((dim_name == NULL) ? 0 : 1, ndims * sizeof(char *));

       for(i=0; i < ndims; i++){
           dim_name[i] = (char *)calloc(NC_MAX_NAME, sizeof(char));
           mem_er((dim_name[i] == NULL) ? 0 : 1, NC_MAX_NAME * sizeof(char));
           nc_inq_dimname(fdat->ncid, i, dim_name[i]);
       }

       var_name = (char **)calloc(nvars, sizeof(char *));
       mem_er((var_name == NULL) ? 0 : 1, nvars * sizeof(char *));

       printf("Available fields are:-   \n\n");

       for(i=0; i < nvars; i++){
          var_name[i] = (char *)calloc(NC_MAX_NAME, sizeof(char));
          mem_er((var_name[i] == NULL) ? 0 : 1, NC_MAX_NAME * sizeof(char));
          nc_inq_varname(fdat->ncid, i, var_name[i]);
          isvar = 0;
          for(j=0; j < ndims; j++) {
              if(!strcmp(dim_name[j], var_name[i])){isvar = 1; break;}
              else if(!strcmp(var_name[i], dim_name[j])){isvar = 1; break;}
          }
          if(!isvar){ 
             printf("Field Id. %d is %s\n", i, var_name[i]);
             fdat->ifield = i;
             ++nfield;
          }

       }

       if(nfield > 1) {
         printf("\nWhich field is required?\n\n");
         if(!(fdat->invar)) scanf("%d", &(fdat->ifield));
         else {
            scanf("%s", fdat->name);
            nc_inq_varid(fdat->ncid, fdat->name, &(fdat->ifield));
         }
       }

       ierr = nc_inq_var(fdat->ncid, fdat->ifield, fdat->name, &type, &nvdim, dims, &natts);
       ierr = nc_inq_attlen(fdat->ncid, fdat->ifield, "long_name", &atlen);
       if(atlen < NC_MAX_NAME && atlen){
          nc_get_att_text(fdat->ncid, fdat->ifield, "long_name", fdat->lname);
          fdat->lname[atlen] = '\0';
       }

       printf("****INFORMATION****, using field Id. %d, %s\n\n", fdat->ifield, fdat->name);

/* check for levels */
       if(fdat->levid >= 0){
          nc_inq_varname(fdat->ncid, fdat->levid, var_name[0]);
          nc_inq_dimid (fdat->ncid, var_name[0], &levid);

          nc_inq_dimlen(fdat->ncid, levid, &nlev);

          if(nlev == 1) fdat->len1[fdat->levind] = 0;
          else {
             llev = (float *)calloc(nlev, sizeof(float));
             mem_er((llev == NULL) ? 0 : 1, nlev * sizeof(float));
             gtmp = (void * )calloc(nlev, sizeof(double));
             mem_er((gtmp == NULL) ? 0 :1, nlev * sizeof(double));
             nc_inq_vartype(fdat->ncid, fdat->levid, &type);
             len1 = 0;
             len2 = nlev;
             nc_read_data(fdat->ncid, fdat->levid, nlev, type, llev, gtmp, &len1, &len2);
             printf("Available levels are:\n\n");
             for(i=0; i < nlev; i++) printf("Level %d is %12.4e\n", i, *(llev + i));
             printf("Which level is required?\n\n");
             if(!(fdat->invar)){
                scanf("%d", &dum);
             }
             else {
                scanf("%e", &alev);
                for(i=0; i < nlev; i++){
                   if(fabs(alev - *(llev + i)) < 1.0e-5) {dum = i; break;}
                }
                
             }
             fdat->len1[fdat->levind] = (size_t)dum;

             printf("****INFORMATION*****, level %e chosen.\n\n", *(llev + dum));
             free(gtmp);
             free(llev);
          }
       }

       if(fdat->addid >= 0){
          nc_inq_varname(fdat->ncid, fdat->addid, var_name[0]);
          nc_inq_dimid (fdat->ncid, var_name[0], &addid);

          nc_inq_dimlen(fdat->ncid, addid, &nadd);

          if(nadd == 1) fdat->len1[fdat->addind] = 0;
          else {
             addd = (float *)calloc(nadd, sizeof(float));
             mem_er((addd == NULL) ? 0 : 1, nadd * sizeof(float));
             gtmp = (void * )calloc(nadd, sizeof(double));
             mem_er((gtmp == NULL) ? 0 :1, nadd * sizeof(double));
             nc_inq_vartype(fdat->ncid, fdat->addid, &type);
             len1 = 0;
             len2 = nadd;
             nc_read_data(fdat->ncid, fdat->addid, nadd, type, addd, gtmp, &len1, &len2);
             printf("Available values of %s are:\n\n", var_name[0]);
             for(i=0; i < nadd; i++) printf("%s %d is %12.4e\n", var_name[0], i, *(addd + i));
             printf("Which value of %s is required?\n\n", var_name[0]);
             if(!(fdat->invar)){ 
                scanf("%d", &dum);
             }
             else {
                scanf("%e", &alev);
                for(i=0; i < nadd; i++){
                   if(fabs(alev - *(addd + i)) < 1.0e-5) {dum = i; break;}
                }
             }
             fdat->len1[fdat->addind] = (size_t)dum;
             printf("****INFORMATION*****, value %e chosen.\n\n", *(addd + dum));
             free(gtmp);
             free(addd);
          }
       }

/* attributes maybe different for new file so check for attributes */

       if(nc_read_att_value(fdat->ncid, fdat->ifield, "missing_value", &(fdat->missing)) == NC_NOERR)
         fdat->imiss = 1;
       else if(nc_read_att_value(fdat->ncid, fdat->ifield, "_FillValue", &(fdat->missing)) == NC_NOERR)
         fdat->imiss = 1;
       else fdat->imiss = 0;

       if(fdat->imiss){
          printf("****INFORMATION****, data has missing values,     \r\n"
                 "                     missing value is %e.            \n\n", fdat->missing);
       }

/* Check for data scaling factor */

       if(nc_read_att_value(fdat->ncid, fdat->ifield, "scale_factor", &(fdat->scale_fac)) == NC_NOERR)
          fdat->iscl = 1;
       else fdat->iscl = 0;

       if(fdat->iscl){
          printf("****INFORMATION****, data has a scaling factor ,     \r\n"
                 "                     scaling factor is %e.           \n\n", fdat->scale_fac);
       }

/* Check for data offset factor */

       if(nc_read_att_value(fdat->ncid, fdat->ifield, "add_offset", &(fdat->fld_offset)) == NC_NOERR)
          fdat->ioff = 1;
       else fdat->ioff = 0;

       if(fdat->ioff){
          printf("****INFORMATION****, data has an offset factor ,     \r\n"
                 "                     scaling factor is %e.           \n\n", fdat->fld_offset);
       }

    }

    else {
       printf("****ERROR****, don't know how to open file with mode %d\n\n", mode);
       exit(1);

    }

    for(i=0; i < ndims; i++) free(dim_name[i]);
    free(dim_name);

    for(i=0; i < ndims; i++) free(var_name[i]);
    free(var_name);


#else

    printf("****WARNING****, no netcdf support, can't open file.\n\n");

#endif

    return fdat;

}
