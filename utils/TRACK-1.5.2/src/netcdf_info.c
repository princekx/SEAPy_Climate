#include <Stdio.h>
#include "grid.h"

#include "netcdf_info.h"

#ifndef NETCDF
NETCDF_INFO *netcdf_info(GRID *gr, int *frnum, char *filnam, int *tl, int *gof)
{
    printf("****INFORMATION****, need to compile with NETCDF support.\n\n");
    return NULL;
}


void netcdf_close(NETCDF_INFO *fdat){}

void handle_error(int status, char *file, int line) {}


#else
#include <stdlib.h>
#include <Math.h>
#include <string.h>
#include "mem_er.h"


#define  ASSUME_COARDS  0

void grid_trans(GRID * , float * , int * , int * );

/* function to get information on the contents of a netcdf file, equivelent to
   reading a header. 
                                                                             */

extern int init;
extern float period;
extern float *abuf;
extern void *databuf;
extern int nhp, shp;
extern int std_x, std_y;

extern FILE *finit;

NETCDF_INFO *netcdf_info(GRID *gr, int *frnum, char *filnam, int *tl, int *gof)
{



    int i=0, j=0;
    int *did=NULL;
    int ipsum='y';
    int mxdim=0;
    int icoards='y';
    int maxdim=0;
    int nfield=0;
    int nd=0;

    int ierr=0;
    int ndims, nvars, ngatts, unlimdimid;
    int ilev=-1;
    int isvar=0;
    int irecord=0;
    int iadd=-1;

    nc_type *var_type=NULL;                   /* variable type */
    int *var_ndims=NULL;                      /* number of dims */
    int **var_dims=NULL;     /* variable shape */
    int *var_natts=NULL;
    char **var_name=NULL;
    char **dim_name=NULL;
    char ***var_dim_name=NULL;

    char lngname[NC_MAX_NAME], units[NC_MAX_NAME];

    void *gtmp=NULL;
    float *xgtmp=NULL;
    float *llev=NULL;
    float *addd=NULL;

    float scale_fac=0.0, goffs=0.0;
    float alev=0.0;

    size_t **var_dim_len=NULL;
    size_t len1[NC_MAX_VAR_DIMS], len2[NC_MAX_VAR_DIMS];
    size_t atlen=0;
    size_t nlev=0, nadd=0;

    NETCDF_INFO *ffield=NULL;

    printf("Print file summary information, 'y' or 'n'\n\n");
    if(init){
       fscanf(finit, "\n");
       ipsum = getc(finit);
    }
    else {
       scanf("\n");
       ipsum = getchar();
       fprintf(finit, "%c\n", ipsum);
    }

/* assign memory for field information */

    ffield = (NETCDF_INFO *)malloc_initl(sizeof(NETCDF_INFO));
    mem_er((ffield == NULL) ? 0 : 1, sizeof(NETCDF_INFO));

    ffield->lngind = ffield->latind = ffield->levind = ffield->timind = ffield->addind = -1;
    ffield->longid = ffield->latid = ffield->levid = ffield->timid = ffield->addid = -1;

/* what type of input for identifying fields and dimension values */

    printf("Use netcdf id's or variable names and dimension values, input '0' or '1'\n\n");
    if(init) fscanf(finit, "%d", &(ffield->invar));
    else {
       scanf("%d", &(ffield->invar));
       fprintf(finit, "%d\n", ffield->invar);
    }
    if(ffield->invar < 0 || ffield->invar > 1) {
      printf("****ERROR****, incorrect input, exiting.\n\n");
      exit(1);
    }

/* open file for read */

    ierr = nc_open(filnam, NC_NOWRITE, &(ffield->ncid));
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__); 

    ierr = nc_inq(ffield->ncid, &ndims, &nvars, &ngatts, &unlimdimid);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

/* read global attributes */

    if(ipsum == 'y' || ipsum == 'Y'){

       atlen = 0;
       lngname[0] = '\0';
       nc_inq_attlen(ffield->ncid, NC_GLOBAL, "title", &atlen);
       if(atlen < NC_MAX_NAME && atlen) {
         nc_get_att_text(ffield->ncid, NC_GLOBAL, "title", lngname);
         lngname[atlen] = '\0';
         printf("Title:- %s\n\n", lngname);
       }

       atlen = 0;
       lngname[0] = '\0';
       nc_inq_attlen(ffield->ncid, NC_GLOBAL, "history", &atlen);
       if(atlen < NC_MAX_NAME && atlen) {
         nc_get_att_text(ffield->ncid, NC_GLOBAL, "history", lngname);
         lngname[atlen] = '\0';
         printf("History:- %s\n\n", lngname);
       }

       atlen = 0;
       lngname[0] = '\0';
       nc_inq_attlen(ffield->ncid, NC_GLOBAL, "Conventions", &atlen);
       if(atlen < NC_MAX_NAME && atlen) {
         nc_get_att_text(ffield->ncid, NC_GLOBAL, "Conventions", lngname);
         lngname[atlen] = '\0';
         printf("Conventions:- %s\n\n", lngname);
       }

    }

/* Is this a COARDS file */

    if(!ASSUME_COARDS){

       printf("Is the data organized according to the COARDS convention, 'y' or 'n'\n\n");
       if(init){
          fscanf(finit, "\n");
          icoards = getc(finit);
       }
       else{
          scanf("\n");
          icoards = getchar();
          fprintf(finit, "%c\n", icoards);
       }
    }

/* assign memory for each dimension name */

    dim_name = (char **)calloc(ndims, sizeof(char *));
    mem_er((dim_name == NULL) ? 0 : 1, ndims * sizeof(char *));

    for(i=0; i < ndims; i++){
        dim_name[i] = (char *)calloc(NC_MAX_NAME, sizeof(char));
        mem_er((dim_name[i] == NULL) ? 0 : 1, NC_MAX_NAME * sizeof(char));
        nc_inq_dimname(ffield->ncid, i, dim_name[i]);
    }


/* assign memory for each variables attributes */

    var_type = (nc_type *)calloc(nvars, sizeof(nc_type));
    mem_er((var_type == NULL) ? 0 : 1, nvars * sizeof(nc_type));

    var_ndims = (int *)calloc(nvars, sizeof(int));
    mem_er((var_ndims == NULL) ? 0 : 1, nvars * sizeof(int));

    var_natts = (int *)calloc(nvars, sizeof(int));
    mem_er((var_natts == NULL) ? 0 : 1, nvars * sizeof(int));

    var_dims = (int **)calloc(nvars, sizeof(int *));
    mem_er((var_dims == NULL) ? 0 : 1, nvars * sizeof(int *));

    var_dim_len = (size_t **)calloc(nvars, sizeof(size_t *));
    mem_er((var_dim_len == NULL) ? 0 : 1, nvars * sizeof(size_t *));

    var_dim_name = (char ***)calloc(nvars, sizeof(char *));
    mem_er((var_dim_name == NULL) ? 0 : 1, nvars * sizeof(char *));

    var_name = (char **)calloc(nvars, sizeof(char *));
    mem_er((var_name == NULL) ? 0 : 1, nvars * sizeof(char *));


    for(i=0; i < nvars; i++){
        var_dims[i] = NULL;
        var_dims[i] = (int *)calloc(NC_MAX_VAR_DIMS, sizeof(int));
        mem_er((var_dims[i] == NULL) ? 0 : 1, NC_MAX_VAR_DIMS * sizeof(int));

        var_name[i] = NULL;
        var_name[i] = (char *)calloc(NC_MAX_NAME, sizeof(char));
        mem_er((var_name[i] == NULL) ? 0 : 1, NC_MAX_NAME * sizeof(char));

    }


    for(i=0; i < nvars; i++){

        ierr = nc_inq_var(ffield->ncid, i, var_name[i], var_type + i, var_ndims + i, var_dims[i], var_natts + i);

        if(var_ndims[i] > maxdim) maxdim = var_ndims[i];


        var_dim_len[i] = NULL;
        var_dim_len[i] = (size_t *)calloc(var_ndims[i], sizeof(size_t));
        mem_er((var_dim_len[i] == NULL) ? 0 : 1, var_ndims[i] * sizeof(size_t));

        var_dim_name[i] = NULL;
        var_dim_name[i] = (char **)calloc(var_ndims[i], sizeof(char *));
        mem_er((var_dim_name[i] == NULL) ? 0 : 1, var_ndims[i] * sizeof(char *));

        for(j=0; j < var_ndims[i]; j++){
            *(var_dim_name[i] + j) = NULL;
            *(var_dim_name[i] + j) = (char *)calloc(NC_MAX_NAME, sizeof(char));
            mem_er((*(var_dim_name[i] + j) == NULL) ? 0 : 1, NC_MAX_NAME * sizeof(char));
            ierr = nc_inq_dim(ffield->ncid, *(var_dims[i] + j), *(var_dim_name[i] + j), var_dim_len[i] + j);
       }

/* write variable summary */

       if(ipsum == 'y' || ipsum == 'Y'){

          var_summary(ffield->ncid, i, var_name[i], var_type[i], var_ndims[i], var_dims[i], var_natts[i], unlimdimid, var_dim_name[i], var_dim_len[i]);

       }

    }

    printf("Available fields are:-   \n\n");

    nfield = 0;


    for(i=0; i < nvars; i++){

       isvar = 0;

       for(j=0; j < ndims; j++) {

           if(!strcmp(dim_name[j], var_name[i])){isvar = 1; break;}
           else if(!strcmp(var_name[i], dim_name[j])){isvar = 1; break;}

       }

       if(!isvar){

          atlen = 0;
          lngname[0] = '\0';
          nc_inq_attlen(ffield->ncid, i, "long_name", &atlen);
          if(atlen < NC_MAX_NAME && atlen){
             nc_get_att_text(ffield->ncid, i, "long_name", lngname);
             lngname[atlen] = '\0';
          }

          atlen = 0;
          units[0] = '\0';
          nc_inq_attlen(ffield->ncid, i, "units", &atlen);
          if(atlen < NC_MAX_NAME && atlen) {
             nc_get_att_text(ffield->ncid, i, "units", units);
             units[atlen] = '\0';
          }

          printf("Field Id. %d is %s %s %s\n", i, var_name[i], lngname, units);

          ffield->ifield = i;

          ++nfield;

       }

    }

    printf("\n\n");

    if(nfield > 1) {
    
       printf("\nWhich field is required?\n\n");

       if(!(ffield->invar)){
          if(init) fscanf(finit, "%d", &(ffield->ifield));
          else{ 
            scanf("%d", &(ffield->ifield));
            fprintf(finit, "%d\n", ffield->ifield);
          }

       }
       else {
          if(init) fscanf(finit, "%s", ffield->name);
          else{
            scanf("%s", ffield->name);
            fprintf(finit, "%s\n", ffield->name);
          }
          nc_inq_varid(ffield->ncid, ffield->name, &(ffield->ifield));
       }
    }

    strncpy(ffield->name, var_name[ffield->ifield], NC_MAX_NAME);

    printf("****INFORMATION*****, field %s chosen.\n\n", ffield->name);

    atlen = 0;
    ffield->lname[0] = '\0';
    nc_inq_attlen(ffield->ncid, ffield->ifield, "long_name", &atlen);
    if(atlen < NC_MAX_NAME && atlen) {
       nc_get_att_text(ffield->ncid, ffield->ifield, "long_name", ffield->lname);
       ffield->lname[atlen] = '\0';
    }

    for(i=0; i < var_ndims[ffield->ifield]; i++)
       strncpy(ffield->dimn[i],  *(var_dim_name[ffield->ifield] + i), NC_MAX_NAME);


    ffield->dattyp = var_type[ffield->ifield];

    atlen = 0;
    ffield->units[0] = '\0';
    nc_inq_attlen(ffield->ncid, ffield->ifield, "units", &atlen);
    if(atlen < NC_MAX_NAME && atlen) { 
       nc_get_att_text(ffield->ncid, ffield->ifield, "units", ffield->units);
       ffield->units[atlen] = '\0';
    }

    if(var_ndims[ffield->ifield] < 2 || var_ndims[ffield->ifield] > 5){
       printf("****ERROR****, the chosen field does not have the correct dimensions, NDIMS = %d.\n\n", var_ndims[ffield->ifield]);
       exit(1);
    }

/* check for NCO record dimension */

/*
    for(i=0; i < var_ndims[ffield->ifield]; i++){
        if(strstr(*(var_dim_name[ffield->ifield] + i), "record") || 
           strstr(*(var_dim_name[ffield->ifield] + i), "Record") ||
           strstr(*(var_dim_name[ffield->ifield] + i), "RECORD")   ) {
               irecord = 1;

               nc_inq_dimid(ffield->ncid, *(var_dim_name[ffield->ifield] + i), &(ffield->addind));
               break;

        }
        
    }
*/

/* NCO record dimension not currently supported */
/*
    if(irecord) {
       printf("****ERROR****, NCO record dimension not currently supported.\n\n");
       exit(1);

    } 
*/

    if(icoards == 'y'){

      if(var_ndims[ffield->ifield] == 2){
         printf("****WARNING****, chosen field has two dimensions only, single field.\n\n");

         if(irecord){
            printf("One of the dimensions is a record dimension, not allowed for single 2D field.\n\n");
            exit(1);
         }

         if(nc_inq_varid(ffield->ncid, *(var_dim_name[ffield->ifield]), &(ffield->latid)) != NC_NOERR) ffield->latid = -1;
         if(nc_inq_varid(ffield->ncid, *(var_dim_name[ffield->ifield]+1), &(ffield->longid)) != NC_NOERR) ffield->longid = -1;

         ffield->lngind = 1;
         ffield->latind = 0;

      }

      else if(var_ndims[ffield->ifield] == 3){
         printf("****WARNING****, chosen field is not 4-dimensional, assuming 3rd dimension is time.\n\n");
         if(nc_inq_varid(ffield->ncid, *var_dim_name[ffield->ifield], &(ffield->timid)) != NC_NOERR) ffield->timid = -1;
         if(nc_inq_varid(ffield->ncid, *(var_dim_name[ffield->ifield]+1), &(ffield->latid)) != NC_NOERR) ffield->latid = -1;
         if(nc_inq_varid(ffield->ncid, *(var_dim_name[ffield->ifield]+2), &(ffield->longid)) != NC_NOERR) ffield->longid = -1;

/* support for NCO record not available yet */

/*         if(irecord && ffield->timid == -1){
            printf("One of the dimensions is a record dimension, do you want to use this \r\n"
                   "as the time dimension, 'y' or 'n'.                                   \n\n");
            scanf("\n");
            if(getchar() == 'y'){}
            else exit(1);
         }

*/


         ffield->lngind = 2;
         ffield->latind = 1;
         ffield->timind = 0;


      }

      else if(var_ndims[ffield->ifield] == 4){
         if(nc_inq_varid(ffield->ncid, *var_dim_name[ffield->ifield], &(ffield->timid)) != NC_NOERR) ffield->timid = -1;
         if(nc_inq_varid(ffield->ncid, *(var_dim_name[ffield->ifield]+1), &(ffield->levid)) != NC_NOERR) ffield->levid = -1;
         if(nc_inq_varid(ffield->ncid, *(var_dim_name[ffield->ifield]+2), &(ffield->latid)) != NC_NOERR) ffield->latid = -1;
         if(nc_inq_varid(ffield->ncid, *(var_dim_name[ffield->ifield]+3), &(ffield->longid)) != NC_NOERR) ffield->longid = -1;
         ffield->lngind = 3;
         ffield->latind = 2;
         ffield->levind = 1;
         ffield->timind = 0;
      }

      else if(var_ndims[ffield->ifield] == 5){
         if(nc_inq_varid(ffield->ncid, *var_dim_name[ffield->ifield], &(ffield->addid)) != NC_NOERR) ffield->addid = -1;
         if(nc_inq_varid(ffield->ncid, *(var_dim_name[ffield->ifield]+1), &(ffield->timid)) != NC_NOERR) ffield->timid = -1;
         if(nc_inq_varid(ffield->ncid, *(var_dim_name[ffield->ifield]+2), &(ffield->levid)) != NC_NOERR) ffield->levid = -1;
         if(nc_inq_varid(ffield->ncid, *(var_dim_name[ffield->ifield]+3), &(ffield->latid)) != NC_NOERR) ffield->latid = -1;
         if(nc_inq_varid(ffield->ncid, *(var_dim_name[ffield->ifield]+4), &(ffield->longid)) != NC_NOERR) ffield->longid = -1;

         ffield->lngind = 4;
         ffield->latind = 3;
         ffield->levind = 2;
         ffield->timind = 1;
         ffield->addind = 0;

         printf("****INFORMATION****, 5th dimension of field is %s\n", *(var_dim_name[ffield->ifield] + ffield->addind));
      }
      else {
         printf("****ERROR****, to many dimensions for current support.\n\n");
         exit(1);

      }


    }


/* Check how many dimensions have been allocated */

    if(ffield->timid >= 0) ++nd;
    if(ffield->levid >= 0) ++nd;
    if(ffield->latid >= 0) ++nd;
    if(ffield->longid >= 0) ++nd;
    if(ffield->addid >= 0) ++nd;

/* Identify longitude, latitude, level and time variable ID's if still unknown */

    if(var_ndims[ffield->ifield] != nd){

      for(i=0; i < var_ndims[ffield->ifield]; i++){

        isvar = 0;
        did = &(ffield->addid);
        strncpy(lngname, *(var_dim_name[ffield->ifield] + i), NC_MAX_NAME);
	
	if(ffield->timid < 0 && (strstr(lngname, "TIM") || strstr(lngname, "Tim")  || 
                         strstr(lngname, "tim") || lngname[0] == 't'        ||
                         lngname[0] == 'T')){
            isvar = 0;
            did = &(ffield->timid);
            if(nc_inq_varid(ffield->ncid, lngname, &(ffield->timid)) != NC_NOERR) ffield->timid = -1;
            else {ffield->timind = i; isvar = 1;}
        }

	else if(ffield->levid < 0 && (strstr(lngname, "LEV") || strstr(lngname, "Lev") || strstr(lngname, "lev") ||
	                              strstr(lngname, "PRESSURE") || strstr(lngname, "Pressure") || strstr(lngname, "pressure"))){
            isvar = 0;
            did = &(ffield->levid);
            if(nc_inq_varid(ffield->ncid, lngname, &(ffield->levid)) != NC_NOERR) ffield->levid = -1;
            else {ffield->levind = i; isvar = 1;}
        }

        else if(ffield->latid < 0 && (strstr(lngname, "LAT") || strstr(lngname, "Lat") || 
                              strstr(lngname, "lat") || strstr(lngname, "y") || strstr(lngname, "Y"))){

            if(strstr(lngname, "y") || strstr(lngname, "Y")) gr->igtyp = 4;
            isvar = 0;
            did = &(ffield->latid);
            if(nc_inq_varid(ffield->ncid, lngname, &(ffield->latid)) != NC_NOERR) ffield->latid = -1;
            else {ffield->latind = i; isvar = 1;} 
        }

        else if(ffield->longid < 0 && (strstr(lngname, "LON") || strstr(lngname, "Lon") || 
                               strstr(lngname, "lon") || strstr(lngname, "x") || strstr(lngname, "X"))){
            if(strstr(lngname, "x") || strstr(lngname, "X")) gr->igtyp = 4;
            isvar = 0;
            did = &(ffield->longid);
            if(nc_inq_varid(ffield->ncid, lngname, &(ffield->longid)) != NC_NOERR) ffield->longid = -1;
            else {ffield->lngind = i; isvar = 1;}
        }

        if(!isvar){
           if(ffield->addid < 0){
             nc_inq_varid(ffield->ncid, lngname, did);	     
             if(nc_inq_varid(ffield->ncid, lngname, &(ffield->addid)) != NC_NOERR) ffield->addid = -1;
             else{
               ffield->addind = i; 
               isvar = 1;
               printf("****INFORMATION****, field has an additional dimension of %s\n", lngname);
             }
           }
           else{
             printf("****ERROR****, additional dimension already assigned, cannot continue.\n\n");
             exit(1);
           }

        }

      }

    }

    if(ffield->longid < 0 || ffield->latid < 0){
      printf("****ERROR****, longitude and latitude dimension ID's unresolved\r\n");
      exit(1);
    }

/* If time and level variable Id. still unresolved */

    if(ffield->timid < 0 || ffield->levid < 0){
       if(var_ndims[ffield->ifield] == 5) {
          printf("****WARNING*****, assuming variable dimension 1 is additional dimension, dimension 2 is time and dimension 3 is level.\n\n");
          ffield->addind = 0;
          ffield->timind = 1;
          ffield->levind = 2;
       }
       else if(var_ndims[ffield->ifield] == 4){
          if(ffield->addind >= 0){
             printf("****WARNING*****, assuming variable dimension 1 is additional and dimension 2 is time.\n\n");
             ffield->addind = 0;
             ffield->timind = 1;
          }
          else {
             printf("****WARNING*****, assuming variable dimension 1 is time and dimension 2 is level.\n\n");
             ffield->timind = 0;
             ffield->levind = 1;
          }
       }

       else if (var_ndims[ffield->ifield] == 3) {
          printf("****WARNING*****, assuming variable dimension 1 is time.\n\n");
          ffield->timind = 0;
       }

    }

    if(var_ndims[ffield->longid] != 1 || var_ndims[ffield->latid] != 1 || 
       (ffield->timid >= 0) ? ((var_ndims[ffield->timid] != 1) ? 1 : 0) : 0 ||
       (ffield->levid >= 0) ? ((var_ndims[ffield->levid] != 1) ? 1 : 0) : 0){
       printf("****WARNING****, possible problem with grid, latitude/longitude not a single dimension.\n\n"); 

    }

/* set number of time frames */

    if(ffield->timind >= 0)
      *frnum = *(var_dim_len[ffield->ifield] + ffield->timind);

    gr->ix = std_x = *(var_dim_len[ffield->ifield] + ffield->lngind);
    gr->iy = std_y = *(var_dim_len[ffield->ifield] + ffield->latind);

    if(ffield->levid >= 0) nlev = *(var_dim_len[ffield->ifield] + ffield->levind);
    if(ffield->addid >= 0) nadd = *(var_dim_len[ffield->ifield] + ffield->addind);

    mxdim = (std_x > std_y) ? std_x : std_y;
    if(mxdim < nlev) mxdim = nlev;
    if(mxdim < *frnum) mxdim = *frnum;
    if(mxdim < nadd) mxdim = nadd;

/* assign memory for grid */

    gr->xgrid = (float * )calloc(gr->ix + 1, sizeof(float));
    mem_er((gr->xgrid == NULL) ? 0 :1, (gr->ix + 1) * sizeof(float));

    gr->ygrid = (float * )calloc(gr->iy, sizeof(float));
    mem_er((gr->ygrid == NULL) ? 0 : 1, gr->iy * sizeof(float));

/* assign memory for input buffers */

    gtmp = (void * )calloc(mxdim, sizeof(double));
    mem_er((gtmp == NULL) ? 0 :1, mxdim * sizeof(double));

    xgtmp = (void * )calloc(std_x, sizeof(float));
    mem_er((xgtmp == NULL) ? 0 :1, std_x * sizeof(float));

/* assign memory for reading field data */

    abuf = (float *)calloc(std_x * std_y, sizeof(float));
    mem_er((abuf == NULL) ? 0 : 1, std_x * std_y * sizeof(float));

/* read times */

    ffield->timval = (double *)calloc(*frnum, sizeof(double));
    mem_er((ffield->timval == NULL) ? 0 : 1, *frnum * sizeof(double));

    len1[0] = 0;
    len2[0] = *frnum;

    if(ffield->timid >= 0){
       nc_read_data_dbl(ffield->ncid, ffield->timid, *frnum, var_type[ffield->timid], ffield->timval, gtmp, len1, len2);
    }
    else {
       for(i=0; i < *frnum; i++) *(ffield->timval + i) = i + 1; 
    }

/* read longitudes */


    len1[0] = 0;
    len2[0] = std_x;

    nc_read_data(ffield->ncid, ffield->longid, std_x, var_type[ffield->longid], xgtmp, gtmp, len1, len2);
    ffield->lnggr = (float *) calloc(std_x, sizeof(float));
    mem_er((ffield->lnggr == NULL) ? 0 : 1, std_x * sizeof(float));
    if(nc_read_att_value(ffield->ncid, ffield->longid, "scale_factor", &scale_fac) == NC_NOERR){
       for(i=0; i < std_x; i++) *(xgtmp + i) *= scale_fac;
    }
    if(nc_read_att_value(ffield->ncid, ffield->longid, "add_offset", &goffs) == NC_NOERR){
       for(i=0; i < std_x; i++) *(xgtmp + i) += goffs;
    }
    memcpy(ffield->lnggr, xgtmp, std_x * sizeof(float));

/* read latitudes */

    len1[0] = 0;
    len2[0] = std_y;

    nc_read_data(ffield->ncid, ffield->latid, std_y, var_type[ffield->latid], gr->ygrid, gtmp, len1, len2);
    ffield->latgr = (float *) calloc(std_y, sizeof(float));
    mem_er((ffield->latgr == NULL) ? 0 : 1, std_y * sizeof(float));
    if(nc_read_att_value(ffield->ncid, ffield->latid, "scale_factor", &scale_fac) == NC_NOERR){
       for(i=0; i < std_y; i++) *(gr->ygrid + i) *= scale_fac;
    }
    if(nc_read_att_value(ffield->ncid, ffield->latid, "add_offset", &goffs) == NC_NOERR){
       for(i=0; i < std_y; i++) *(gr->ygrid + i) += goffs;
    }
    memcpy(ffield->latgr, gr->ygrid, std_y * sizeof(float));

/* perform grid transformations */


    if(gr->igtyp < 4) grid_trans(gr, xgtmp, tl, gof);
    else {
       printf("****INFORMATION****, grid is not a latitude/longitude grid.\n\n");
       for(i= 0; i< gr->ix; i++) *((gr->xgrid)+i) = *(xgtmp + i);
       
       printf("the current grid dimensions are %d * %d \n\n", gr->ix, gr->iy);
    }


/* read levels */


    if(ffield->levid >= 0){

       llev = (float *)calloc(nlev, sizeof(float));
       mem_er((llev == NULL) ? 0 : 1, nlev * sizeof(float));

       len1[0] = 0;
       len2[0] = nlev;

       nc_read_data(ffield->ncid, ffield->levid, nlev, var_type[ffield->levid], llev, gtmp, len1, len2);


       atlen = 0;
       units[0] = '\0';
       nc_inq_attlen(ffield->ncid, ffield->levid, "units", &atlen);
       if(atlen < NC_MAX_NAME && atlen) {
         nc_get_att_text(ffield->ncid, ffield->levid, "units", units);
         units[atlen] = '\0';
       }

       if(nlev > 1){

          printf("The levels are:-     \n\n");
          for(i=0; i < nlev; i++){
              printf("Level %3d is %12.4e %s\n", i, *(llev + i), units);
          }

          printf("\nWhich level is required?\n\n");

          if(!(ffield->invar)){
             if(init) fscanf(finit, "%d", &ilev);
             else{ 
                scanf("%d", &ilev);
                fprintf(finit, "%d\n", ilev);
             }
          }
          else {
             if(init) fscanf(finit, "%e", &alev);
             else {
                scanf("%e", &alev);
                fprintf(finit, "%e\n", alev);
             }
             for(i=0; i < nlev; i++){
                if(fabs(alev - *(llev + i)) < 1.0e-5) {ilev = i; break;}
             } 
          }

          if(ilev < 0 || ilev >= nlev){
             printf("****ERROR****, no such level Id.\n\n");
             exit(1);
          }

       }

       else ilev = 0;

       ffield->levval = *(llev + ilev);

       printf("****INFORMATION*****, level %e chosen.\n\n", ffield->levval);

       free(llev);

    }

    if(ffield->addid >= 0){

       addd = (float *)calloc(nadd, sizeof(float));
       mem_er((addd == NULL) ? 0 : 1, nadd * sizeof(float));

       len1[0] = 0;
       len2[0] = nadd;

       nc_read_data(ffield->ncid, ffield->addid, nadd, var_type[ffield->addid], addd, gtmp, len1, len2);

       atlen = 0;
       units[0] = '\0';
       nc_inq_attlen(ffield->ncid, ffield->addid, "units", &atlen);
       if(atlen < NC_MAX_NAME && atlen) {
         nc_get_att_text(ffield->ncid, ffield->addid, "units", units);
         units[atlen] = '\0';
       }

       if(nadd > 1){

	  nc_inq_varname(ffield->ncid, ffield->addid, lngname);
          printf("The values of dimension %s are:-     \n\n", lngname);
          for(i=0; i < nadd; i++){
              printf("%s %3d is %12.4e %s\n", lngname, i, *(addd + i), units);
          }

          printf("\nWhich %s value is required?\n\n", lngname);

          if(!(ffield->invar)){
             if(init) fscanf(finit, "%d", &iadd);
             else{ 
                scanf("%d", &iadd);
                fprintf(finit, "%d\n", iadd);
             }
          }
          else {
             if(init) fscanf(finit, "%e", &alev);
             else {
                scanf("%e", &alev);
                fprintf(finit, "%e\n", alev);
             }
             for(i=0; i < nadd; i++){
                if(fabs(alev - *(addd + i)) < 1.0e-5) {iadd = i; break;}
             }
          }

          if(iadd < 0 || iadd >= nadd){
             printf("****ERROR****, no such level Id.\n\n");
             exit(1);
          }

       }
       else iadd = 0;

       ffield->addval = *(addd + iadd);

       printf("****INFORMATION*****, added dimension value %e chosen.\n\n", ffield->addval);

       free(addd);

    }

/* update array counters */

    ffield->len1[ffield->timind] = 0;
    ffield->len1[ffield->latind] = 0;
    ffield->len1[ffield->lngind] = 0;

    ffield->len2[ffield->timind] = 1;
    ffield->len2[ffield->latind] = std_y;
    ffield->len2[ffield->lngind] = std_x;

    if(ffield->levid >= 0){
       ffield->len1[ffield->levind] = ilev;
       ffield->len2[ffield->levind] = 1;
    }

    if(ffield->addid >= 0){
       ffield->len1[ffield->addind] = iadd;
       ffield->len2[ffield->addind] = 1;
    }

/* Check for missing data value */

    if(nc_read_att_value(ffield->ncid, ffield->ifield, "missing_value", &(ffield->missing)) == NC_NOERR)
       ffield->imiss = 1;
    else if(nc_read_att_value(ffield->ncid, ffield->ifield, "_FillValue", &(ffield->missing)) == NC_NOERR)
       ffield->imiss = 1;
    else ffield->imiss = 0;

    if(ffield->imiss){
       printf("****INFORMATION****, data has missing values,     \r\n"
              "                     missing value is %e.            \n\n", ffield->missing);
    }

/* Check for data scaling factor */

    if(nc_read_att_value(ffield->ncid, ffield->ifield, "scale_factor", &(ffield->scale_fac)) == NC_NOERR)
       ffield->iscl = 1;
    else ffield->iscl = 0;
    if(ffield->iscl){
       printf("****INFORMATION****, data has a scaling factor ,     \r\n"
              "                     scaling factor is %e.           \n\n", ffield->scale_fac);
    }

/* Check for data offset factor */

    if(nc_read_att_value(ffield->ncid, ffield->ifield, "add_offset", &(ffield->fld_offset)) == NC_NOERR)
       ffield->ioff = 1;
    else ffield->ioff = 0;

    if(ffield->ioff){
       printf("****INFORMATION****, data has an offset factor ,     \r\n"
              "                     scaling factor is %e.           \n\n", ffield->fld_offset);
    }

/* assign memory for reading field data */

    abuf = (float *)calloc(std_x * std_y, sizeof(float));
    mem_er((abuf == NULL) ? 0 : 1, std_x * std_y * sizeof(float));

    databuf = (void *)malloc_initl(std_x * std_y * sizeof(double));
    mem_er((databuf == NULL) ? 0 : 1, std_x * std_y * sizeof(double));    

    free(gtmp);
    free(xgtmp);

    free(var_type);
    free(var_natts);

    for(i=0; i <ndims; i++) free(dim_name[i]);
    free(dim_name);

    for(i=0; i < nvars; i++){

        free(var_dims[i]);
        free(var_name[i]);
        free(var_dim_len[i]);
        for(j=0; j < var_ndims[i]; j++) free(*(var_dim_name[i] + j));
        free(var_dim_name[i]);

    }

    free(var_ndims);
    free(var_dims);
    free(var_name);
    free(var_dim_len);
    free(var_dim_name);

    return ffield;

}


void handle_error(int status, char *file, int line) {

   printf("%s in file %s at line %d\n", nc_strerror(status), file, line);
   exit(1);
   return;

}


void var_summary(int ncid, int ivar, char *name, nc_type type, int ndims, int *dims, int natts, int unlim, char **dimnam, size_t *dimlen)
{
    int i;
    int vid;

    size_t atlen = 0;
    nc_type attype;
    char dlen[NC_MAX_NAME], units[NC_MAX_NAME], lngnam[NC_MAX_NAME];
    char typstr[NC_MAX_NAME];

    float val=0.0;

    printf("============================================================\n");
    printf("------------------%s------------------\n", name);
    printf("Variable ID = %d\n", ivar);
    printf("Variable name = %s\n", name);
    vartyp(type, typstr);
    printf("Variable type  is %s\n", typstr);

    printf("------------Dimensions-------------\n");
    printf("Number of dimensions = %d\n", ndims);
    for(i=0; i < ndims; i++) {
        if(*(dims + i) == unlim) sprintf(dlen, "UNLIMITED, Current value = %d", (int)(*(dimlen + i)));
        else sprintf(dlen, "%d", (int)(*(dimlen + i)));
        vid = -1;
        nc_inq_varid(ncid, dimnam[i], &vid);
        printf("  Dimension %d is dimension ID %d, variable ID %d, name %s, length %s \n", i, *(dims+i), vid, dimnam[i], dlen);
    }
    printf("Number of variable attributes is %d\n", natts);
    printf("------------Attributes-------------\n");
    printf("Available attributes:- \n\n");
    for(i=0; i <natts; i++){
       nc_inq_attname(ncid, ivar, i, lngnam);
       nc_inq_atttype(ncid, ivar, lngnam, &attype);
       vartyp(attype, typstr);     
       printf("Att. %2d is %s of type %s\r\n", i, lngnam, typstr);
    }
    atlen = 0;
    lngnam[0]='\0';
    nc_inq_attlen(ncid, ivar, "long_name", &atlen);
    if(atlen < NC_MAX_NAME && atlen) {
      nc_get_att_text(ncid, ivar, "long_name", lngnam);
      lngnam[atlen] = '\0';
      printf("Long Name:- %s\n", lngnam);
    }
    atlen = 0;
    units[0]='\0';
    nc_inq_attlen(ncid, ivar, "units", &atlen);
    if(atlen < NC_MAX_NAME && atlen) {
      nc_get_att_text(ncid, ivar, "units", units);
      units[atlen] = '\0';
      printf("Units:- %s\n", units);
    }

    if(nc_read_att_value(ncid, ivar, "missing_value", &val) == NC_NOERR)
       printf("%s:- %e\n", "missing_value", val);


    printf("\n\n\n\n");

    return;

}


void vartyp(nc_type attype, char *type)
{

    switch(attype){
        case NC_BYTE:
           sprintf(type, "%s", "BYTE");
           break;
        case NC_CHAR:
           sprintf(type, "%s", "CHAR");
           break;
        case NC_SHORT:
           sprintf(type, "%s", "SHORT INTEGER");
           break;
         case NC_INT:
           sprintf(type, "%s", "INTEGER");
           break;
         case NC_FLOAT:
           sprintf(type, "%s", "FLOAT");
           break;
         case NC_DOUBLE:
           sprintf(type, "%s", "DOUBLE");
           break;
         default:
           printf("Variable type %d is unknown\n", attype); 
           break;
    }

}

int nc_read_data(int ncid, int var_id, int dim, int vartype, float *arr, void *buff, size_t *len1, size_t *len2)
{

       int i;
       int ierr=NC_NOERR;

       switch(vartype){
             case NC_SHORT:
               ierr = nc_get_vara_short(ncid, var_id, len1, len2, (short int *)buff);
               for(i=0; i < dim; i++) *(arr + i) = (float)(*((short int *)buff + i));
               break;
             case NC_INT:
               ierr = nc_get_vara_int(ncid, var_id, len1, len2, (int *)buff);
               for(i=0; i < dim; i++) *(arr + i) = (float)(*((int *)buff + i));
               break;
             case NC_FLOAT:
               ierr = nc_get_vara_float(ncid, var_id, len1, len2, (float *)buff);
               for(i=0; i < dim; i++) *(arr + i) = *((float *)buff + i);
               break;
             case NC_DOUBLE:
               ierr = nc_get_vara_double(ncid, var_id, len1, len2, (double *)buff);
               for(i=0; i < dim; i++) *(arr + i) = (float)(*((double *)buff + i));
               break;
             default:
               printf("****ERROR****, data type %d not allowed for data Id. %d\n\n", vartype, var_id);
               exit(1);
       }



    return ierr;

}

int nc_read_data_dbl(int ncid, int var_id, int dim, int vartype, double *arr, void *buff, size_t *len1, size_t *len2)
{

       int i;
       int ierr=NC_NOERR;

       switch(vartype){
             case NC_SHORT:
               ierr = nc_get_vara_short(ncid, var_id, len1, len2, (short int *)buff);
               for(i=0; i < dim; i++) *(arr + i) = (double)(*((short int *)buff + i));
               break;
             case NC_INT:
               ierr = nc_get_vara_int(ncid, var_id, len1, len2, (int *)buff);
               for(i=0; i < dim; i++) *(arr + i) = (double)(*((int *)buff + i));
               break;
             case NC_FLOAT:
               ierr = nc_get_vara_float(ncid, var_id, len1, len2, (float *)buff);
               for(i=0; i < dim; i++) *(arr + i) = (double)(*((float *)buff + i));
               break;
             case NC_DOUBLE:
               ierr = nc_get_vara_double(ncid, var_id, len1, len2, (double *)buff);
               for(i=0; i < dim; i++) *(arr + i) = *((double *)buff + i);
               break;
             default:
               printf("****ERROR****, data type %d not allowed for data Id. %d\n\n", vartype, var_id);
               exit(1);
       }



    return ierr;

}


int nc_read_att_value(int ncid, int var_id, char *att, float *val)
{

    int ierr=NC_NOERR;

    nc_type attype;
    size_t atlen=0;

    short int sp;
    int ip;
    float fp;
    double dp;

    atlen = 0;
    nc_inq_att(ncid, var_id, att, &attype, &atlen);
    if(atlen){
       
       switch(attype){
           case NC_SHORT:
              ierr = nc_get_att_short(ncid, var_id, att, &sp);
              *val = (float) sp;
              break;
            case NC_INT:
              ierr = nc_get_att_int(ncid, var_id, att, &ip);
              *val = (float) ip;
              break;
            case NC_FLOAT:
              ierr = nc_get_att_float(ncid, var_id, att, &fp);
              *val = fp;
              break;
            case NC_DOUBLE:
              ierr = nc_get_att_double(ncid, var_id, att, &dp);
              *val = (float) dp;
              break;
            default:
              printf("Variable type %d is unknown\n", attype); 
              ierr = 1;
              break;
       }

     }

     else ierr = 1;

     if(atlen > 1){
        printf("****WARNING****, NETCDF:- this attribute has more than one value.\n\n");

     }


     return ierr;

}

void netcdf_close(NETCDF_INFO *fdat)
{

    int ierr=0;

    ierr = nc_close(fdat->ncid);
    if(ierr != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);
    free(fdat->lnggr);
    free(fdat->latgr);
    free(fdat->timval);
    free(fdat);

    return;

}

#endif

