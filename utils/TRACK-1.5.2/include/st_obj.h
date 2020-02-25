/* structure definition for object points */

#define  OBJWRT   1        /* determines whether the object point data is written to file, */
                           /* and whether it is read from file, can save on memory.        */
                           /* '0' all object data written, memory free'd after plotting    */
                           /* '1' prompted for whether to write/read object point data     */

#define  FPTWRT   0        /* determines whether the feature point data is free'd after    */
                           /* plotting, can save on memory while finding feature points    */
                           /* in subsequent frames. Feature points read in for tracking    */
                           /* following completion of feature finding for all frames.      */

struct extent{
       int x1;
       int x2;
       int y1;
       int y2;
};

struct point{
       int x;
       int y;
       float val;
};



/* structure definition for objects */

struct object{
      struct point *pt;
      struct features *fet;
      struct extent *ext;
      struct boundary_pt *bound;
      int point_num;
      int bound_num;
      int mem_pt_size;        /* only used in periodic boundary situations */
      int lab;
      char b_or_i;
};

/* structure definition for object groups per frame */

struct frame_objs {
      struct object *objs;
      int obj_num;
      int tot_f_f_num;
      int b_state;            /* boundary state, float or integer         */
      int frame_id;
      int nmiss;              /* number of missing frames to next frame   */
};

/* structure definition for frame dates and times
   this is associated with old missing frame checking. */

struct times{
      int year;
      int month;
      int day;
      int hour;
};


float obj_xreal(int );  /* function prototype for real grid point evaluation */
