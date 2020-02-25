#include <Stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Math.h>
#include <sys/types.h>
#include "mem_er.h"
#include "grid.h"
#include "file_handle.h"
#include "files_out.h"
#include "pp.h"
#include "utf.h"
#include "netcdf_info.h"

#define  NCHRB  30

/* function to perform a converion to the binary formatted field type */

void write_header(GRID * , FILE * );
int read_field(float * ,float * , float , FILE * , int , int , int , int , int );

extern GRID *gr;
extern char *chrfld;

extern char *fext;

extern int iext;
extern int form;

void convert(FILE *fdat, int fr1, int fri, int frl)

{

    int dim=gr->ix * gr->iy;
    int fram=0;

    int nf;

    off_t place1, place2, chrnum=0;
 
    float *ap=NULL;

    FILE *filspec=NULL;

    char conv[MAXCHR];

    ap=(float *)calloc(dim, sizeof(float));
    mem_er((ap == NULL) ? 0: 1, dim * sizeof(float));

    strncpy(conv, CONVERT, MAXCHR);
    if(iext) strcpy(strstr(conv, EXTENSION), fext);


    filspec = open_file(conv, "w");

/* write header */

    write_header(gr, filspec);

    
    nf = 0;


/* read and process fields write out in simple binary format for now */

    if(form == 4){
       if(((NETCDF_INFO *)fdat)->levid >=0) ((NETCDF_INFO *)fdat)->len1[((NETCDF_INFO *)fdat)->levind] = 0;
       if(((NETCDF_INFO *)fdat)->addid >=0) ((NETCDF_INFO *)fdat)->len1[((NETCDF_INFO *)fdat)->addind] = 0;
    }

    while(nf <= frl){

         if(form != 4){

            if(!nf) {

               place1 = ftello(fdat); 

 
               if(read_field(ap, NULL, 1.0, fdat, 1.0, 'n', 'n', '0', 'n'))break;        

               if(fr1 == 1) ++nf;
 
               place2 = ftello(fdat);
               chrnum = place2 - place1;

               if(fr1 > 1){

                  fseeko(fdat, (fr1-2)*chrnum, ORIGIN);
                  if(read_field(ap, NULL, 1.0, fdat, 1.0, 'n', 'n', '0', 'n'))break;

                  nf = fr1;  

         

               }

               nf += fri;

/* need to reset header strings for new grid */

            }

            else {

               fseeko(fdat, (fri-1)*chrnum, ORIGIN);
               if(read_field(ap, NULL, 1.0, fdat, 1.0, 'n', 'n', '0', 'n'))break;
          
               nf += fri;

            }


         }

         else {

            if(!nf) nf = fr1;

            ((NETCDF_INFO *)fdat)->iframe = nf - 1;
 
            if(read_field(ap, NULL, 1.0, fdat, 1.0, 'n', 'n', '0', 'n'))break;        
            nf += fri;

         }


/* write out new field */

         ++fram;


         fprintf(filspec, "FRAME %6d\n", fram);

         fwrite(ap, dim * sizeof(float), 1, filspec);
         fprintf(filspec, "\n");


         if(nf > frl) break;


    }


    fseeko(filspec, (off_t)0, FSTART);
    fprintf(filspec, "%6d %6d %6d\n", gr->ix, gr->iy, fram);

    close_file(filspec, conv);

    if(chrfld) free(chrfld);

    free(ap);

    return;

}
