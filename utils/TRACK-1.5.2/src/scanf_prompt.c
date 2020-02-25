#include <stdio.h>
#include <stdarg.h>

#define  PROMPT  "track->"

/* wrapper for scanf to provide a prompt */


int scanf_prompt(char *format, ...)
{

   va_list argl;
   int iret=0;

   va_start(argl, format);

   printf("%s", PROMPT);
   iret= scanf(format, argl);   

   va_end(argl);
   return iret;

}
