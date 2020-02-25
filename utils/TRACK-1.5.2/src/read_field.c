#include <Stdio.h>

int std_read_field(float * , float * , float , FILE * , int , int , int , int, int );
int utf3_read_field(float * , float * , float , FILE * , int , int , int , int, int );
int utf4_read_field(float * , float * , float , FILE * , int , int , int , int, int );

extern int form;
extern int utfv;

/* function to choose and read different formated fields */

int read_field(float *ap, float *ata, float scale, FILE *fdatin, int sign, int as, int of, int hemi, int tscl)

{

     int ir=0;

     if(form == 0 || form == 1 || form == 3 || form == 4) 

       ir = std_read_field(ap, ata, scale, fdatin, sign, as, of, hemi, tscl);

     else if(form == 2){

       if(utfv == 3) 
         ir = utf3_read_field(ap, ata, scale, fdatin, sign, as, of, hemi, tscl);
       else
         ir = utf4_read_field(ap, ata, scale, fdatin, sign, as, of, hemi, tscl);

     }

     return ir;

}

