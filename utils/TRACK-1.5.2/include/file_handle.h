/* header for file handing function prototypes */

#define   APPS    "_"        /* file name append character */

#define  ORIGIN  SEEK_CUR    /* random file access from current file position */
#define  FSTART  SEEK_SET    /* random file access from start of file         */

#define  NUMLINES     5

/* Next parameter controls what is printed as an input file is read */

#define  PTSTR        1   /* >0 to print header summaries, */
                          /*  0 to print '#',              */
                          /* <0 to print nothing.          */


FILE *open_file(char * , char * );
void close_file(FILE * , char * );
char *file_exist(char *, char * , char * );
int fexist(char *, char * );
