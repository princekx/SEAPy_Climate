/* #include <stdlib.h> */

#ifdef MDB
  extern int  malloc_verify();
  extern int malloc_debug(int );
#define MEM_TYPE_DB    1       /* switch for type of memory debugging */

#if MEM_TYPE_DB              /* if true execute */

#define mem_verify()  {if(!(malloc_verify())){(void)fprintf(stderr,"***error***, bad memory block found: file \"%s\", line %d\n", __FILE__, __LINE__);exit(1);}}
#define mem_debug(A) {}

#else

#define mem_debug(A)  (int)malloc_debug(A)
#define mem_verify()  {}

#endif                       /* MEM_TYPE_DB */
#else
#define mem_debug(A) {}
#define mem_verify()  {}
#endif                       /* MDB */

#ifndef NDEBUG

/* Old memory error traps, now rationalised to a single trap */

/*# define mem_er(ex)	{mem_verify();if (!(ex)){(void)fprintf(stderr,"***error***, memory allocation failed: file \"%s\", line %d\n", __FILE__, __LINE__);exit(1);}}
# define mem_er_realloc(ex, size)	{mem_verify();if (!(ex) && size != 0){(void)fprintf(stderr,"***error***, memory allocation failed: file \"%s\", line %d\n", __FILE__, __LINE__);exit(1);}}  */

#define mem_er(ex, size)	{mem_verify();if (!(ex) && (size) != 0){(void)fprintf(stderr,"***error***, memory allocation of %ld bytes failed: file \"%s\", line %d\n", (long int)(size), __FILE__, __LINE__);exit(1);}}

#endif

#undef  free

#define free(p)  {if((p)){free(p); p=NULL;}}

/* prototype for malloc_init */

void *malloc_initl(int );
void *realloc_n(void * , int );

