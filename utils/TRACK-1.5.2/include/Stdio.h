#if defined(CYGWIN) || defined(MAC)
#include "/usr/include/stdio.h"
#else
#include <stdio.h>
#endif

#if defined(SUNOS5) || defined(AIX)
#if !defined(_LARGEFILE_SOURCE)
#define ftello ftell
#define fseeko fseek
#endif
#if !defined(_LARGE_FILES)
#define ftello ftell
#define fseeko fseek
#endif
#endif

#if defined(SUNOS4)

#include <stdarg.h>

extern  int   printf(char *, ...);
extern  int   fprintf(FILE *, char *, ...);
extern  int   scanf(char *, ...);
extern  int   fscanf(FILE *, char *, ...);
extern  int   sscanf(char *, char *, ...);
extern  void  perror(char *);
extern  int   fclose(FILE *);
extern  int   fseek(FILE *, long, int);
extern  char  _filbuf(FILE *);
extern  int   fflush(FILE *);
extern  int   fread(void *, int, int, FILE *);
extern  int   fwrite(void *, int, int, FILE *);
extern  char  memset(char *, int, int);
extern  void  setbuf(FILE * , char *);

#endif


