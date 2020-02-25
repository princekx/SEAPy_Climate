#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>

#define MAXCHR 100

#  define RANK_var 4
#  define RANK_lat 1
#  define RANK_lon 1
#  define RANK_level 1
#  define RANK_time 1

long int new_time(long int , int );

void
check_err(const int stat, const int line, const char *file) {
    if (stat != NC_NOERR) {
	   (void) fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
        exit(1);
    }
}


void main(int argc, char *argv[])

{

    int i, j;
    int dim;
    int ia;
    int ierr=0;
    int ig1, ig2;
    int ntime;
    int stat;
    int ict=0;

    int  ncid;

    int var_id;
    int lat_id;
    int lon_id;
    int level_id;
    int time_id;
    int ts=0;

   /* dimension ids */
    int lat_dim;
    int lon_dim;
    int level_dim;
    int time_dim;

    int var_dims[RANK_var];
    int lat_dims[RANK_lat];
    int lon_dims[RANK_lon];
    int level_dims[RANK_level];
    int time_dims[RANK_time];

   /* dimension lengths */
    int level_len = 1;
    int time_len, lon_len, lat_len;

    static float level[] = {1};

    FILE *fin=NULL;
    FILE *nc1=NULL;

    char str[MAXCHR];
    char outfil[MAXCHR];


    float *lon=NULL, *lat=NULL;
    float *field=NULL;

    size_t time_start[RANK_time];
    size_t time_count[RANK_time];

    size_t len1[NC_MAX_DIMS]={0, 0, 0, 0, 0};
    size_t len2[NC_MAX_DIMS]={1, 1, 0, 0, 0};

    long int *time=NULL, st=0;


    if(argc < 3 || argc > 4){

       printf("***Usage***, bin2drs [input file] [output file] [Add Time (1), optional]\n\n");

       exit(1);

    }

    strcpy(outfil, argv[2]);
    strcat(outfil, ".nc");

    if(argc == 4) sscanf(argv[3], "%d", &ict);

    if((nc1=fopen(outfil, "r"))){

       printf("****ERROR****, netCDF file already exist for output filename:- \r\n"
              "               %s \r\n" 
              "               rename or remove files:-                      \r\n"
              "               %s \r\n"
              "               and try again.\n\n", argv[2], outfil);

       if(nc1) fclose(nc1);
       exit(1);

    }

    fin = fopen(argv[1], "r");

    if(!fin){

       printf("****ERROR***, unable to open file:- \r\n"
              "              %s\n\n", argv[1]);
       exit(1);

    }

    fgets(str, 100, fin);
    sscanf(str, "%d %d %d", &lon_len, &lat_len, &ntime);
    dim = lon_len * lat_len;
    time_len = ntime;

    len2[2] = lat_len;
    len2[3] = lon_len;


    lon = (float *)calloc(lon_len, sizeof(float));
    if(!lon){

       printf("****ERROR****, unable to assign memory.\n\n");
       exit(1);

    }

    lat = (float *)calloc(lat_len, sizeof(float));
    if(!lat){

       printf("****ERROR****, unable to assign memory.\n\n");
       exit(1);

    }

    field = (float *)calloc(dim, sizeof(float));
    if(!field){

       printf("****ERROR****, unable to assign memory.\n\n");
       exit(1);

    }

    time = (long int *)calloc(ntime, sizeof(long int));
    if(!time){

       printf("****ERROR****, unable to assign memory.\n\n");
       exit(1);

    }

    if(ict){
       printf("What is the starting time and time step?\n");
       scanf("%ld %d", &st, &ts);
       for(i=0; i < ntime; i++) {
           time[i] = st;
printf("%ld\n", time[i]);
           st = new_time(st, ts);
       }
    }
    else
       for(i=0; i < ntime; i++) time[i] = i+1;


    ia = 0;
    for(i=0; i < lon_len; i++) {
        fscanf(fin, "%f", lon+i);
        ++ia;
        if(ia == 10){fgets(str, 100, fin); ia = 0;}
    }
    if(ia) fgets(str, 100, fin);
    ia = 0;
    for(i=0; i < lat_len; i++) {
        fscanf(fin, "%f", lat+i);
        ++ia;
        if(ia == 10){fgets(str, 100, fin); ia = 0;}
    }
    if(ia) fgets(str, 100, fin);


    /* enter define mode */
/*    stat = nc_create(outfil, NC_CLOBBER | NC_64BIT_OFFSET, &ncid); */
    stat = nc_create(outfil, NC_CLOBBER, &ncid);
    check_err(stat,__LINE__,__FILE__);

    /* define dimensions */
    stat = nc_def_dim(ncid, "lat", lat_len, &lat_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "lon", lon_len, &lon_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "level", level_len, &level_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "time", time_len, &time_dim);
    check_err(stat,__LINE__,__FILE__);

    var_dims[0] = time_dim;
    var_dims[1] = level_dim;
    var_dims[2] = lat_dim;
    var_dims[3] = lon_dim;

    lon_dims[0] = lon_dim;
    stat = nc_def_var(ncid, "lon", NC_FLOAT, RANK_lon, lon_dims, &lon_id);
    check_err(stat,__LINE__,__FILE__);

    lat_dims[0] = lat_dim;
    stat = nc_def_var(ncid, "lat", NC_FLOAT, RANK_lat, lat_dims, &lat_id);
    check_err(stat,__LINE__,__FILE__);

    level_dims[0] = level_dim;
    stat = nc_def_var(ncid, "level", NC_FLOAT, RANK_level, level_dims, &level_id);
    check_err(stat,__LINE__,__FILE__);

    time_dims[0] = time_dim;
    stat = nc_def_var(ncid, "time", NC_INT, RANK_time, time_dims, &time_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "var", NC_FLOAT, RANK_var, var_dims, &var_id);
    check_err(stat,__LINE__,__FILE__);


   /* assign attributes */
    stat = nc_put_att_text(ncid, var_id, "long_name", 8, "variable");
    check_err(stat,__LINE__,__FILE__);
     stat = nc_put_att_text(ncid, lat_id, "units", 7, "degrees");
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, lon_id, "units", 7, "degrees");
    check_err(stat,__LINE__,__FILE__);

   /* leave define mode */
    stat = nc_enddef (ncid);

    check_err(stat,__LINE__,__FILE__);


    stat = nc_put_var_float(ncid, level_id, level);
    check_err(stat,__LINE__,__FILE__);


    time_start[0] = 0;
    time_count[0] = time_len;
    stat = nc_put_vara_long(ncid, time_id, time_start, time_count, time);
    check_err(stat,__LINE__,__FILE__);   

    stat = nc_put_var_float(ncid, lat_id, lat);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_put_var_float(ncid, lon_id, lon);
    check_err(stat,__LINE__,__FILE__);

    time_len = 1;			/* number of records of time data */
    time_start[0] = 1;
    time_count[0] = ntime;

    for(i=0; i<ntime; i++){


        fgets(str, 100, fin);
        printf("%s\n", str);

        fread(field, dim*sizeof(float), 1, fin);

        len1[0] = i;
        stat = nc_put_vara_float(ncid, var_id, len1, len2, field);
        check_err(stat,__LINE__,__FILE__);

        if(!fgets(str, 100, fin))break;

    }


    stat = nc_close(ncid);
    check_err(stat,__LINE__,__FILE__);


    fclose(fin);

    free(lon);
    free(lat);
    free(field);

    return;

}
