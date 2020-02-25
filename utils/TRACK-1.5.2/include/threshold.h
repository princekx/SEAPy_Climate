/* function prototypes for threshold */

void arrayd(struct image *, float *, float *, FILE * , float ,int , int , int );
void border_expand_object(struct frame_objs * , int , int , float );
/* void centroid(struct frame_objs * ); */
void extnpts_to_domain(struct frame_objs * , int *, int * );
void filtd(struct frame_objs * , int , int , int , float );
void hierarc_segment(struct image ** , struct image * , int , struct frame_objs * , int , int , int );
void objectd(struct frame_objs *, FILE * ,int , int );
void object_filter(struct frame_objs * , int , int , int );
void object_local_maxs(struct frame_objs * , struct frame_objs * , int , int );
void object_realf(struct frame_objs * , float , int , int , int , int *);
void periodic_boundary(struct frame_objs * , int * , int * );
int powi(int , int );
void proj_interp(float * , double * , int  , int , ...);
int query_interp(int );
int read_field(float * ,float * , float , FILE * , int , int , int , int , int );
struct frame_objs *read_obd(FILE * , float * , float * , int * , int * );
void surfit(double * , int , int , ... );
void non_lin_opt(struct frame_objs * , int , int );
int time_avg(FILE * , int , int , int , int );
void object_dist_trans(struct frame_objs * , struct frame_objs * , int * , int * ,int , int );
void object_smint(struct frame_objs * , double * , struct savedat * , struct rspline * , struct sp_dat * , int );
void anisotropy(struct frame_objs * , int );
void mfilt_obj(struct frame_objs * );
void fpt_intensity(struct frame_objs * , float , float );
void frame_boundaries(struct frame_objs * , struct boundary_cntl * , int * , int *);
void shape_setup(struct boundary_cntl * , int);
void feature_to_object(struct frame_objs * , struct frame_objs * , int );
void putback_area(struct frame_objs * , struct frame_objs * );
float sys_area(struct boundary_pt * , int , int , int , float , int * );
void correct_fpt(struct frame_objs * , int );
void adjust_fpt(struct frame_objs * , struct sp_dat * , float * , int );
long int new_time(long int , int );
struct frame_objs *mfilt_obj_dist(struct frame_objs * , float , int );
void write_header(GRID * , FILE * );
