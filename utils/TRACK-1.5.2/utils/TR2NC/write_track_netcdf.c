#include <Stdio.h>


#ifndef NETCDF
#include "splice.h"
void write_track_netcdf(struct tot_tr *trr , int trnum)
{
    printf("****INFORMATION****, need to compile with NETCDF support.\n\n");
    return;
}

#else

#include <netcdf.h>
#include <stdlib.h>
#include <string.h>
#include "file_cat_out.h"
#include "splice.h"
#include "mem_er.h"

#define  RANK_TR  1
#define  RANK_VAR 1

void handle_error(int , char * , int );
void write_attribute(char ** , char * , int , int , int , long int * , int * );
long int new_time(long int , int );
void get_var_name(char ** , char * , int , char * );
long int julian(long int );

extern int nfld, nff;
extern int *nfwpos;

void write_track_netcdf(struct tot_tr *trr , int trnum, char *filnm, int itrtyp, char **metad, int nmeta)
{
    int i, j, k;

    int ncid=0;
    int itotrec=0;
    int ierr=NC_NOERR;
    int drecid=0, dtrid=0;
    int varid[1]={0};
    int irmax=0;
    int trdef[3]={0, 0, 0};
    int trid[9]={0, 0, 0, 0};
    int itim_typ=0;
    int recpos=0;
    int ifldc=0;
    int ichnm=0;
    int tstep=0;
    int maxtim=0, ttim=0;

    size_t strt[1]={0}, stct[1]={0};

    int *add_def=NULL;

    int *itrstart=NULL, *itrnum=NULL;
    int *itrid=NULL;
    int *dtimi=NULL;
    int *itrindx=NULL;
    
    long int tstart=0;
    long int *ttimes=NULL;
    long int ntim=0;
    long int jstart=0;

    char nfilout[MAXCHR], num[5];
    char **add_var_nm=NULL;
    char varnm[MAXCHR];
    char vartag[MAXCHR];

    float *fdump[5]={NULL, NULL, NULL, NULL, NULL};
    float *dlat=NULL, *dlng=NULL;

    float missval=ADD_UNDEF;

    double *dtimf=NULL;
    double *dtimc=NULL;
    double ftp=0.0;

    struct tot_tr *atr=NULL;
    struct fet_pt_tr *fpt=NULL;

    strncpy(nfilout, filnm, MAXCHR);
    strcat(nfilout, ".nc");

    if(itrtyp != 's' && itrtyp !='v'){
       printf("****ERROR****, unknown option '%c' for writing file.\n\n", itrtyp);
       exit(1);
    }

    printf("****INFORMATION****, writing netcdf file %s\n", nfilout);

/*    printf("Do you want to change default field variable names, 'y' or 'n'.\n\n");
    scanf("\n");
    if(getchar() == 'y') ichnm=1; */

    if(trr->time) itim_typ = 1;

    if(nff){
       add_def = (int *)calloc(nfld, sizeof(int));
       mem_er((add_def == NULL) ? 0 : 1, nfld*sizeof(int));

       add_var_nm = (char **)calloc(nfld, sizeof(char *));
       mem_er((add_var_nm == NULL) ? 0 : 1, nfld*sizeof(char *));

       for(i=0; i < nfld; i++) {
          add_var_nm[i] = (char *)calloc(MAXCHR, sizeof(char));
          mem_er((add_var_nm[i] == NULL) ? 0 : 1, MAXCHR*sizeof(char));
       }
    }

/* determine max record length */

    itrstart = (int *)calloc(trnum, sizeof(int));
    mem_er((itrstart == NULL) ? 0 : 1, trnum*sizeof(int));

    itrnum = (int *)calloc(trnum, sizeof(int));
    mem_er((itrnum == NULL) ? 0 : 1, trnum*sizeof(int));

    itrid = (int *)calloc(trnum, sizeof(int));
    mem_er((itrid == NULL) ? 0 : 1, trnum*sizeof(int));

    
    for(i=0; i < trnum; i++) {
        atr = trr + i;
        *(itrid + i) = i;
        if(atr->num > irmax) irmax = atr->num;
        *(itrstart + i) = itotrec;
        *(itrnum + i) = atr->num;
        itotrec += atr->num;
	if(!itim_typ){
	   ttim = (atr->trpt + atr->num - 1)->fr_id;
	   if(ttim > maxtim) maxtim = ttim;
	}
    }

    itrindx = (int *)calloc(irmax, sizeof(int));
    mem_er((itrindx == NULL) ? 0 : 1, itotrec*sizeof(int));

    if(itrtyp == 's'){
      fdump[0] = (float *)calloc(irmax, sizeof(float));
      mem_er((fdump[0] == NULL) ? 0 : 1, itotrec*sizeof(float));
    }
    else if(itrtyp == 'v'){
      for(i=0; i < 5; i++){
         fdump[i] = (float *)calloc(irmax, sizeof(float));
         mem_er((fdump[i] == NULL) ? 0 : 1, itotrec*sizeof(float));
      }
    }

    dlat = (float *)calloc(irmax, sizeof(float));
    mem_er((dlat == NULL) ? 0 : 1, itotrec*sizeof(float));

    dlng = (float *)calloc(irmax, sizeof(float));
    mem_er((dlng == NULL) ? 0 : 1, itotrec*sizeof(float));
    
    if(metad){
       dtimc = (double *)calloc(irmax, sizeof(double));
       mem_er((dtimc == NULL) ? 0 : 1, itotrec*sizeof(double)); 
       
       if(!itim_typ){    
          ttimes = (long int *)calloc(maxtim, sizeof(long int));
          mem_er((ttimes == NULL) ? 0 : 1, maxtim*sizeof(double));     
       }
           
    }    
    else{

       if(itim_typ) {
          dtimf = (double *)calloc(irmax, sizeof(double));
          mem_er((dtimf == NULL) ? 0 : 1, itotrec*sizeof(double));
       }
       else {
          dtimi = (int *)calloc(irmax, sizeof(int));
          mem_er((dtimi == NULL) ? 0 : 1, itotrec*sizeof(int));
       }
    }
    
/* open and create file */

    if((ierr = nc_create(nfilout,  NC_CLOBBER, &ncid)) != NC_NOERR) 
       handle_error(ierr, __FILE__, __LINE__);
       
/* write global attributes */
    if(metad) write_attribute(metad, "GLOBAL-ATTRIBUTES", ncid, nmeta, NC_GLOBAL, NULL, NULL);

/* define dimensions */

    if((ierr = nc_def_dim(ncid, "tracks", trnum, &dtrid)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    if((ierr = nc_def_dim(ncid, "record", NC_UNLIMITED, &drecid)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

/* define variables */

    varid[0] = dtrid;
    if((ierr = nc_def_var(ncid, "TRACK_ID", NC_INT, RANK_TR, varid, &trdef[0])) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    if((ierr = nc_def_var(ncid, "FIRST_PT", NC_INT, RANK_TR, varid, &trdef[1])) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    if((ierr = nc_def_var(ncid, "NUM_PTS", NC_INT, RANK_TR, varid, &trdef[2])) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

/* attributes */

    if((ierr = nc_put_att_int(ncid, trdef[0], "add_fld_num", NC_INT, 1, &nff)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    if((ierr = nc_put_att_int(ncid, trdef[0], "tot_add_fld_num", NC_INT, 1, &nfld)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    if((ierr = nc_put_att_int(ncid, trdef[0], "loc_flags", NC_INT, nff, nfwpos)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);
       
    if(metad) write_attribute(metad, "TRACK_ID", ncid, nmeta, trdef[0], NULL, NULL);   
    
    if(metad) write_attribute(metad, "NUM_PTS", ncid, nmeta, trdef[2], NULL, NULL);  
    

/* main track info */

    varid[0] = drecid;

    if((ierr = nc_def_var(ncid, "index", NC_INT, RANK_VAR, varid, &trid[0])) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);
    
    if(metad){
    
       if((ierr = nc_def_var(ncid, "time", NC_DOUBLE, RANK_VAR, varid, &trid[1])) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);
       
       write_attribute(metad, "TIME", ncid, nmeta, trid[1], &tstart, &tstep);
       
       jstart = julian(tstart);
       if(!itim_typ){
          if(!tstart || !tstep){
	     printf("****ERROR****, start time and time step not specified in time attribute input.\n\n");
	     exit(1);
	  }
	  else {
	     *ttimes = tstart;

	     for(i=1; i < maxtim; i++) *(ttimes + i) = new_time(*(ttimes + i -1), tstep);
	  }
                
       }
       
    }
    else{
       if(itim_typ){
          if((ierr = nc_def_var(ncid, "time", NC_DOUBLE, RANK_VAR, varid, &trid[1])) != NC_NOERR)
             handle_error(ierr, __FILE__, __LINE__);
       }
       else {
          if((ierr = nc_def_var(ncid, "time", NC_INT, RANK_VAR, varid, &trid[1])) != NC_NOERR)
             handle_error(ierr, __FILE__, __LINE__);
       }    
    }

    if((ierr = nc_def_var(ncid, "longitude", NC_FLOAT, RANK_VAR, varid, &trid[2])) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);
    
    if(metad) write_attribute(metad, "LONGITUDE", ncid, nmeta, trid[2], NULL, NULL);

    if((ierr = nc_def_var(ncid, "latitude", NC_FLOAT, RANK_VAR, varid, &trid[3])) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    if(metad) write_attribute(metad, "LATITUDE", ncid, nmeta, trid[3], NULL, NULL);

    if(itrtyp == 's'){        /* scaler file */

       if(ichnm){
          printf("What is the tracking variable name?\n\n");
          scanf("%s", varnm);
       }
       else {strcpy(varnm, "intensity");}
       
       if(metad) get_var_name(metad, "INTENSITY", nmeta, varnm);		 

       if((ierr = nc_def_var(ncid, varnm, NC_FLOAT, RANK_VAR, varid, &trid[4])) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);

       if(metad) write_attribute(metad, "INTENSITY", ncid, nmeta, trid[4], NULL, NULL);	  
       
    }
    else {   /* vector file */

       if((ierr = nc_def_var(ncid, "speed", NC_FLOAT, RANK_VAR, varid, &trid[4])) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);
	  
       if(metad) write_attribute(metad, "SPEED", ncid, nmeta, trid[4], NULL, NULL);

       if((ierr = nc_def_var(ncid, "gwthr", NC_FLOAT, RANK_VAR, varid, &trid[5])) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);
	  
       if(metad) write_attribute(metad, "GROWTH-RATE", ncid, nmeta, trid[5], NULL, NULL);

       if((ierr = nc_def_var(ncid, "tend", NC_FLOAT, RANK_VAR, varid, &trid[6])) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);
	  
       if(metad) write_attribute(metad, "TENDENCY", ncid, nmeta, trid[6], NULL, NULL);

       if((ierr = nc_def_var(ncid, "velX", NC_FLOAT, RANK_VAR, varid, &trid[7])) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);
	 
       if(metad) write_attribute(metad, "X-VEL", ncid, nmeta, trid[7], NULL, NULL);

       if((ierr = nc_def_var(ncid, "velY", NC_FLOAT, RANK_VAR, varid, &trid[8])) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);
	  
       if(metad) write_attribute(metad, "Y-VEL", ncid, nmeta, trid[8], NULL, NULL);

    }

/* additional fields */

    if(nff){
       for(i=0; i < nff; i++){

/* generate variable names */

           sprintf(num, "_%d", i+1);

           if(*(nfwpos + i)){
              strcpy(add_var_nm[ifldc], "longitude");
              strcat(add_var_nm[ifldc], num);
              if((ierr = nc_def_var(ncid, add_var_nm[ifldc], NC_FLOAT, RANK_VAR, varid, &add_def[ifldc])) != NC_NOERR)
                 handle_error(ierr, __FILE__, __LINE__);
              if((ierr = nc_put_att_float(ncid, add_def[ifldc], "missing_value", NC_FLOAT, 1, &missval)) != NC_NOERR)
                 handle_error(ierr, __FILE__, __LINE__);
	      if(metad) write_attribute(metad, "LONGITUDE", ncid, nmeta, add_def[ifldc], NULL, NULL);	      
              ++ifldc;
              strcpy(add_var_nm[ifldc], "latitude");
              strcat(add_var_nm[ifldc], num);
              if((ierr = nc_def_var(ncid, add_var_nm[ifldc], NC_FLOAT, RANK_VAR, varid, &add_def[ifldc])) != NC_NOERR)
                 handle_error(ierr, __FILE__, __LINE__);
              if((ierr = nc_put_att_float(ncid, add_def[ifldc], "missing_value", NC_FLOAT, 1, &missval)) != NC_NOERR)
                 handle_error(ierr, __FILE__, __LINE__);
              if(metad) write_attribute(metad, "LATITUDE", ncid, nmeta, add_def[ifldc], NULL, NULL);
              ++ifldc;

              if(ichnm){
                 printf("What is the additional variable name?\n\n");
                 scanf("%s", varnm);
                 strcpy(add_var_nm[ifldc], varnm);
                 strcat(add_var_nm[ifldc], num);
              }
              else {
                 strcpy(add_var_nm[ifldc], "addfld");
                 strcat(add_var_nm[ifldc], num);
              }
    
	      if(metad) {
	         strcpy(vartag, "VAR");
                 strcat(vartag, num);
	         get_var_name(metad, vartag, nmeta, add_var_nm[ifldc]);		 
	      }

              if((ierr = nc_def_var(ncid, add_var_nm[ifldc], NC_FLOAT, RANK_VAR, varid, &add_def[ifldc])) != NC_NOERR)
                 handle_error(ierr, __FILE__, __LINE__);
              if((ierr = nc_put_att_float(ncid, add_def[ifldc], "missing_value", NC_FLOAT, 1, &missval)) != NC_NOERR)
                 handle_error(ierr, __FILE__, __LINE__);
              if(metad) write_attribute(metad, vartag, ncid, nmeta, add_def[ifldc], NULL, NULL);
              ++ifldc;


           }
           else{
              strcpy(add_var_nm[ifldc], "addfld");
              strcat(add_var_nm[ifldc], num);
	      
	      if(metad) {
	         strcpy(vartag, "VAR");
                 strcat(vartag, num);
	         get_var_name(metad, vartag, nmeta, add_var_nm[ifldc]);		 
	      }	      
	      
              if((ierr = nc_def_var(ncid, add_var_nm[ifldc], NC_FLOAT, RANK_VAR, varid, &add_def[ifldc])) != NC_NOERR)
                 handle_error(ierr, __FILE__, __LINE__);
              if((ierr = nc_put_att_float(ncid, add_def[ifldc], "missing_value", NC_FLOAT, 1, &missval)) != NC_NOERR)
                 handle_error(ierr, __FILE__, __LINE__);		 
              if(metad) write_attribute(metad, vartag, ncid, nmeta, add_def[ifldc], NULL, NULL);
              ++ifldc;
           }


       }
    }


/* end define mode */

    if((ierr = nc_enddef(ncid)) != NC_NOERR) handle_error(ierr, __FILE__, __LINE__);

    strt[0] = 0;
    stct[0] = trnum;
    if((ierr = nc_put_vara_int(ncid, trdef[0], strt, stct, itrid)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    if((ierr = nc_put_vara_int(ncid, trdef[1], strt, stct, itrstart)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    if((ierr = nc_put_vara_int(ncid, trdef[2], strt, stct, itrnum)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    for(i=0; i < trnum; i++){
       atr = trr + i;

       for(j=0; j < atr->num; j++){
           fpt = atr->trpt + j;
           *(itrindx + j) = j;
	   if(metad){
	      if(itim_typ){ 
/*	         ftp = (double)(fpt->time / 100);
		 *(dtimc + j) = ftp + (((double)(fpt->time)  - ftp * 100.0) / 24.0); */
		 *(dtimc + j) = (double)(julian(fpt->time) - jstart);	 		 
	      }
	      else {
	         ntim = *(ttimes + fpt->fr_id - 1);
/*		 ftp = (double)(ntim / 100);
		 *(dtimc + j) = ftp + (((double)(ntim)  - ftp * 100.0) / 24.0); */
		 *(dtimc + j) = (double)(julian(ntim) - jstart);
	      }
	   }
	   else{
              if(itim_typ) *(dtimf + j) = fpt->time;
              else *(dtimi + j) = fpt->fr_id;
	   }
           *(dlng + j) = fpt->xf;
           *(dlat + j) = fpt->yf;
           *(fdump[0] + j) = fpt->zf;
           if(itrtyp == 'v') {
              *(fdump[1] + j) = fpt->gwthr;
              *(fdump[2] + j) = fpt->tend;
              *(fdump[3] + j) = fpt->vec[0];
              *(fdump[4] + j) = fpt->vec[1];
           }
       }

       strt[0] = recpos;
       stct[0] = atr->num;

       if((ierr = nc_put_vara_int(ncid, trid[0], strt, stct, itrindx)) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);

       if(metad){
             if((ierr = nc_put_vara_double(ncid, trid[1], strt, stct, dtimc)) != NC_NOERR)
                handle_error(ierr, __FILE__, __LINE__);       
       }
       else{
          if(itim_typ){
             if((ierr = nc_put_vara_double(ncid, trid[1], strt, stct, dtimf)) != NC_NOERR)
                handle_error(ierr, __FILE__, __LINE__);
          }
          else {
             if((ierr = nc_put_vara_int(ncid, trid[1], strt, stct, dtimi)) != NC_NOERR)
                handle_error(ierr, __FILE__, __LINE__);
          }
       }

       if((ierr = nc_put_vara_float(ncid, trid[2], strt, stct, dlng)) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);

       if((ierr = nc_put_vara_float(ncid, trid[3], strt, stct, dlat)) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);

       if((ierr = nc_put_vara_float(ncid, trid[4], strt, stct, fdump[0])) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);

       if(itrtyp == 'v'){
          if((ierr = nc_put_vara_float(ncid, trid[5], strt, stct, fdump[1])) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);
          if((ierr = nc_put_vara_float(ncid, trid[5], strt, stct, fdump[2])) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);
          if((ierr = nc_put_vara_float(ncid, trid[5], strt, stct, fdump[3])) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);
          if((ierr = nc_put_vara_float(ncid, trid[5], strt, stct, fdump[4])) != NC_NOERR)
          handle_error(ierr, __FILE__, __LINE__);
       }

/* additional fields */

       if(nff){

          ifldc = 0;

          for(j=0; j < nff; j++){
              if(*(nfwpos + j)){
                 for(k=0; k < atr->num; k++){
                    fpt = atr->trpt + k;
                    *(dlng + k) = fpt->add_fld[ifldc];
                    *(dlat + k) = fpt->add_fld[ifldc + 1];
                    *(fdump[0] + k) = fpt->add_fld[ifldc + 2];    
                 }

                 if((ierr = nc_put_vara_float(ncid, add_def[ifldc], strt, stct, dlng)) != NC_NOERR)
                     handle_error(ierr, __FILE__, __LINE__);
                 if((ierr = nc_put_vara_float(ncid, add_def[ifldc+1], strt, stct, dlat)) != NC_NOERR)
                     handle_error(ierr, __FILE__, __LINE__);
                 if((ierr = nc_put_vara_float(ncid, add_def[ifldc+2], strt, stct, fdump[0])) != NC_NOERR)
                     handle_error(ierr, __FILE__, __LINE__);

                 ifldc += 3;

              }
              else{
                 for(k=0; k < atr->num; k++){
                    fpt = atr->trpt + k;
                    *(fdump[0] + k) = fpt->add_fld[ifldc];
                 }

                 if((ierr = nc_put_vara_float(ncid, add_def[ifldc], strt, stct, fdump[0])) != NC_NOERR)
                     handle_error(ierr, __FILE__, __LINE__);

                 ifldc += 1;
              }

          }

       }

       recpos += atr->num;

    }  
    

    if((ierr = nc_close(ncid)) != NC_NOERR)
       handle_error(ierr, __FILE__, __LINE__);

    if(itrtyp == 's') {
       free(fdump[0]);
    }
    else {
      for(i=0; i < 5; i++) free(fdump[i]);
    }

    free(itrstart);
    free(itrnum);
    free(itrid);
    free(dlat);
    free(dlng);
    if(itim_typ) {free(dtimf);}
    else {free(dtimi);}
    if(ttimes)free(ttimes)
    if(dtimc) free(dtimc);
    if(nff) {
      free(add_def);
      for(i=0; i < nfld; i++) free(add_var_nm[i]);
      free(add_var_nm);
    }

    return;
}

void handle_error(int status, char *file, int line) {

   printf("%s in file %s at line %d\n", nc_strerror(status), file, line);
   exit(1);
   return;

}

void write_attribute(char **att, char *catt, int ncid, int natt, int varid, long int *tstart, int *tstep)
{
    int i=0;
    int iarray=0;
    int ncatt=0;
    int ierr=NC_NOERR;
    int len=0;
    
    char *name=NULL;
    char *value=NULL;
    char *sfnd=NULL;
    
    char strcp[MAXCHR];
    
    double scale=0.0;

    for(i=0; i < natt; i++){
        if((sfnd=strstr(att[i], catt)) != NULL){
	   iarray = i + 2;
	   sscanf(att[i + 1], "%d", &ncatt);
	   break;
	}
    }
    
    if(!sfnd){
       printf("*****WARNING****, no attributes for variable %s.\n", catt); 
       return;   
    }
    
    for(i=0; i < ncatt; i++){
        strncpy(strcp, att[iarray + i], MAXCHR);
        name = strtok(att[iarray + i], ":");
	value = strtok(NULL, ":");
	len = strlen(value) - 1;
	if(strstr(name, "start")){
	   sscanf(value, "%ld", tstart);
	}
	if(strstr(name, "step")){
	   sscanf(value, "%d", tstep);
	}
	
	if(strstr(name, "scale_factor")){
	   sscanf(value, "%lf", &scale);
           if((ierr = nc_put_att_double(ncid, varid, name, NC_DOUBLE, 1, &scale)) != NC_NOERR)
              handle_error(ierr, __FILE__, __LINE__);	   
	}
	else{	
           if((ierr = nc_put_att_text(ncid, varid, name, len, value)) != NC_NOERR)
              handle_error(ierr, __FILE__, __LINE__);
	} 
	
        strncpy(att[iarray + i], strcp, MAXCHR);
    }

    return;
}

void get_var_name(char **att, char *catt, int natt, char *vname)
{
    int i=0;
    int iarray=0;
    int ncatt=0;
    int len=0;
    
    char *name=NULL;
    char *value=NULL;
    char *sfnd=NULL;
    
    char strcp[MAXCHR];

    for(i=0; i < natt; i++){
        if((sfnd=strstr(att[i], catt)) != NULL){
	   iarray = i + 2;
	   sscanf(att[i + 1], "%d", &ncatt);
	   break;
	}
    }
    
    if(!sfnd){
       printf("*****WARNING****, no attributes for variable %s.\n", catt); 
       return;   
    }
    
    for(i=0; i < ncatt; i++){
        strncpy(strcp, att[iarray + i], MAXCHR);
        name = strtok(att[iarray + i], ":");
	value = strtok(NULL, ":");
	len = strlen(value) - 1;
	if(strstr(name, "standard_name")){
	   sscanf(value, "%s", vname);
	}
	
	strncpy(att[iarray + i], strcp, MAXCHR);

    }

    return;
}


#endif




