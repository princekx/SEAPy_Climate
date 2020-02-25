/* Structure definition for netcdf info */

#define NC_CREATE_MODE  0
#define NC_OPEN_MODE    1
#define NC_SAME         2
#define  NC_MY_MAX_DIM  5

#ifdef NETCDF

#include <netcdf.h>

typedef struct netinfo {
    char name[NC_MAX_NAME];     /* variable name string.           */
    char lname[NC_MAX_NAME];    /* variable long name.             */
    char units[NC_MAX_NAME];    /* units of field                  */
    char dimn[NC_MY_MAX_DIM][NC_MAX_NAME];  /* dimension names     */
    int ncid;                   /* netcdf file descriptor          */
    int lngind;                 /* longitude array index           */
    int latind;                 /* latitude array index            */
    int levind;                 /* level array index               */
    int timind;                 /* time array index                */
    int addind;                 /* additional dim. array index     */
    int timid;                  /* Id. for time variable.          */
    int levid;                  /* Id. for level variable.         */
    int longid;                 /* Id. for longitude variable.     */
    int latid;                  /* Id. for latitude variable.      */
    int addid;                  /* Id. for additional dimension var.*/
    size_t len1[NC_MY_MAX_DIM]; /* start positions of field array  */
    size_t len2[NC_MY_MAX_DIM]; /* length positions of field array */
    int iframe;                 /* current frame Id.               */
    int nframe;                 /* actual frame Id.                */
    int ifield;                 /* required field Id.              */
    int dattyp;                 /* data type                       */
    int imiss;                  /* missing value flag              */
    int iscl;                   /* data scaling factor flag        */
    int ioff;                   /* data offset flag                */
    int invar;                  /* input type, id's or names/values*/
    float missing;              /* missing value                   */
    float *latgr;               /* pointer to latitude grid        */
    float *lnggr;               /* pointer to longitude grid       */
    float levval;               /* leval value                     */
    float addval;               /* additional dimension value      */
    float scale_fac;            /* scaling factor                  */
    float fld_offset;           /* field offset                    */
    double *timval;             /* pointer to time values          */
} NETCDF_INFO;


void var_summary(int , int , char * , nc_type , int , int * , int , int , char ** , size_t * );
void vartyp(nc_type , char * );
int nc_read_data(int , int , int , int , float * , void * , size_t * , size_t * );
int nc_read_data_dbl(int , int , int , int , double * , void * , size_t * , size_t * );
int nc_read_att_value(int , int , char * , float * );

#else

/* dummy structure if no NETCDF */

typedef struct netinfo {
    int ncid;    
    int lngind;
    int latind;
    int levind; 
    int timind; 
    int addind;
    int timid; 
    int levid;  
    int longid; 
    int latid;  
    int addid;    
    int len1[NC_MY_MAX_DIM];    
    int len2[NC_MY_MAX_DIM];  
    int iframe;     
    int nframe;
    int imiss;      
    float missing;  
    float levval;  
} NETCDF_INFO;


#endif

NETCDF_INFO *netcdf_info(GRID * , int * , char * , int * , int * );
int netcdf_read_field(void * );
void handle_error(int , char * , int );
void write_netcdf_field(float * , void * );
NETCDF_INFO *nc_define(NETCDF_INFO * , char * );
void netcdf_close(NETCDF_INFO * );
NETCDF_INFO *nc_clone(NETCDF_INFO * , char * , int );

